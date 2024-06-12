/****************************************

  simlib_coadd.c  
  Created Oct 2, 2007 by R.Kessler

  Convert simlib with many exposures per measurement into a
  "COADD" library. Note that his program has the option to
  sum or average the exposures.  The default 'SUM' option 
  is recommended.

  Unlike SDSS, most surveys take multiple exposures in the same 
  passband and same night, and then combine the exposures into a 
  single co-added measurement for the passband/night.  The SIMLIB 
  files therefore can have many repeat measurements for each 
  effective measurement, resulting hundreds of measures on a
  lightcurve.  To avoid memory problems with the simulation & 
  fitter, this program 'coadds' these SIMLIB measures into a 
  single effective measurement.  

  The first use of this progam is for the SNLS/first-season SIMLIB
  [since it's not fair to ask Renault to do more work], using the
  AVG option. However, this program could also be useful to coadd
  simlibs for DES, LSST, ESSENCE.

  For each sequence within a night/passband, and defining N=1 if 
  images are summed and N = Nexposure if images are averaged, the 
  following SIMLIB quantities are modified:

   MJD       -> average MJD
   CCD gain  -> [sum G_i] * N/Nexpose
   CCD noise -> sqrt (sum of squares) / N
   SKYSIG    -> sqrt (sum of squares) / N
   PSF1      -> average PSF
   ZPT       -> 2.5 * log10 [ \sum 10^(0.4*ZPT_i) ] # coadd
   ZPTSIG    -> average ZPTSIG


  Usage (default is to sum images):
    simlib_coadd <simlib_file>

  Output is a new file called <simlib_file>.COADD

  Options:
    simlib_coadd <simlib_file> MWEBV      (compute MWEBV from Schlagel maps)
    simlib_coadd <simlib_file> --TDIF <TDIF>    (default is .4 days)
    simlib_coadd <simlib_file> --MINOBS <MINOBS>  (default is 2 epochs)

    simlib_coadd <simlib_file> --MJD_MIN <MJD_MIN>
    simlib_coadd <simlib_file> --MJD_MAX <MJD_MAX>
    simlib_coadd <simlib_file> --LIBID_MAX <LIBID_MAX>
    simlib_coadd <simlib_file> --LIBID_MIN <LIBID_MIN>

    simlib_coadd <simlib_file>  MJD_DMP
    simlib_coadd <simlib_file>  SUM   (add fluxes; default)
    simlib_coadd <simlib_file>  AVG   (average fluxes)
    simlib_coadd <simlib_file>  SNLS  (add Vega zeropts)

    simlib_coadd <simlib_file> SORT_BAND (sort by band before coadd) 
         # e.g., g,r,i,g,r,i -> gg,rr,ii so that each band is coadded.

  History
  ---------


 Jun 20 2017: 
    + define double MJD[MXMJD] and stop using float MJD

 Apr 3 2018:
   + long int IDEXPT -> char STRINGID[20] to allow ID*NEXPOSE
     e.g., 123438*2

 Jun 29 2018
   + fix logic appending COMMENT lines to original comments
   + final call to insert_NLIBID to add "NLIBID: <NLIBID>" key
     in global header.

 Jan 07 2021
   + sum NEXPOSE for IDEXPT*NEXPOSE column (for LSST DDF)
   + enable reading gzipped SIMLIB file (no need to speficy .gz)
   + replace local MADABORT with errmsg function in sntools.c  
   + print summary info including min/max NOBS and min/max MJD
   + misc code cleanup.

 Feb 25 2021:
   + refactor to require DOCANA keys, and append keys inside this
     YAML block
   + abort if there is no DOCANA block
   + abort on legacy COMMENT: keys

 Mar 01 2021: adapt to work for NEA column replacing 3 PSF columns.

 Jul 15 2021: 
   + MXCHAR_LINE-> 200 (was 120)
   + fix string-cat bug for ptrhead (caught by segfault in debug mode)

 Dec 6 2021: fix bug setting NEA or PSF label for each LIBID.
             Only impacts human-readability; no effect on sim.

 Mar 2022 refactor to enable co-adding bands that are no sequential
   + use SIMLIB_CONTENTS_DEF typedef 
   + all obs indices starts at 0 instead of 1
   + add SORT_BAND option to sort by filter so that non-sequential bands
     within a night are co-added.

 Jun 24 2022
   + fix bug computing min/max MJD for in update_summary_info();
     no impact on SIMLIB contents.

 Jun 12 2024: 
   + MXLINE_HEADER -> 200 (was 30)
   + define SIMLIB_INPUT_TEMP as global instead of local to avoid
     too much local memory inside function (was causing crash)

***************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h" 
#include "simlib_tools.c"

#define MXMJD     100000    // max MJDs per LIBID
#define MXLIBID   6000
#define MWEBV_MAX    2.0    // reject fields with such large MWEBV
#define MXLINE_HEADER 200    // max lines for simlib header
#define MXCHAR_LINE   200     

// global variables.

FILE *fp_simlib_input;
FILE *fp_simlib_output;

int  NLINE_HEADER ;        // number of header lines
int  NLINE_HEADER_ADD ;    // number of comment-header lines to add
bool PSF_NEA_UNIT;

char HEADER[MXLINE_HEADER][MXCHAR_LINE];
char HEADER_ADD[MXLINE_HEADER][MXCHAR_LINE];

// define user-input options

struct INPUTS {
  double MJD_MIN ;     // min MJD to process
  double MJD_MAX ;     // max MJD to process
  double MAXTDIF_COMBINE ;  // max time-diff to combine
  int    LIBID_MAX ;   // max LIBID to process
  int    LIBID_MIN ;   // min LIBID to process
  int    MINOBS_ACCEPT ;  // min # epochs to accept

  int   OPT_MWEBV ;     // compute MWEBV from Schlagel maps
  int   OPT_MJD_DMP ;   // dump flag
  int   OPT_SUM ;       // SUM exposures (default)
  int   OPT_AVG ;       // average of exposures (for strange SNLS lib?)
  int   OPT_SNLS;       // SNLS options

  int   OPT_SORT_BAND;  // Mar 2022
} INPUTS ;

// -----------------------
int  NLIBID_FOUND;
int  NLIBID_COADD;

char SIMLIB_FILTERS[100];  // from FILTERS: <filtlist> header

int  LTRACE;
bool REQUIRE_DOCANA_COADD ;

typedef struct {
  char  FILE[MXPATHLEN];

  // define contents for 1 LIBID entry in SIMLIB/cadence file
  int   LIBID ;
  int   NOBS ;
  int   NOBS_ACCEPT ;
  char  FIELDNAME[20];
  float INFO_HEAD[NPAR_HEAD];  // RA,DECL, MWEBV, PIXSIZE, Z, PEAKMJD

  char STRING_IDEXPT[MXMJD][20]; // Apr 2018
  int  NEXPOSE_IDEXPT[MXMJD];    // Jan 2021
  int  IDEXPT[MXMJD];

  char BAND[MXMJD][2];
  int  IFILT[MXMJD];

  // info is: CCDGAIN, CCDNOISE, SKYSIG, PSF[0-2], ZPTAVG, ZPTSIG, MAG
  double INFO_OBS[MXMJD][NPAR_OBS];

} SIMLIB_CONTENTS_DEF ;


SIMLIB_CONTENTS_DEF SIMLIB_INPUT ;
SIMLIB_CONTENTS_DEF SIMLIB_INPUT_TEMP ;
SIMLIB_CONTENTS_DEF SIMLIB_OUTPUT ;

// Jan 2021: store info for summary
struct {
  double MJD_MAX, MJD_MIN;
  int    NOBS_MIN, NOBS_MAX;
} SUMMARY_INFO;

// declare functions

void  print_simlib_coadd_help(void);
void  parse_args(int argc, char **argv);
void  SIMLIB_open_read();
void  SIMLIB_read(int *RDSTAT);
void  SIMLIB_sort_band(void);
void  SIMLIB_coadd(void);
void  insert_NLIBID(void);

void  init_summary_info(void);
void  update_summary_info(int obs);
void  print_summary_info(void);

void dmp_trace_main(char *string);

// define SNLS-specific functions
double ZPTOFF_SNLS(char *cfilt);

// define Galactic extiction
void MWgaldust(double RA, double DECL, double *XMW, double *MWEBV ); 

void copy_SIMLIB_CONTENTS(SIMLIB_CONTENTS_DEF *CONTENTS_INP,
			  SIMLIB_CONTENTS_DEF *CONTENTS_OUT );

void copy_SIMLIB_CONTENTS_OBS(SIMLIB_CONTENTS_DEF *CONTENTS_INP,
			      SIMLIB_CONTENTS_DEF *CONTENTS_OUT, 
			      int o_inp,  int o_out ) ;

// ****************************************
int main(int argc, char **argv) {

  int LIBID, obs, ifilt, NOBS ;
  int RDSTAT ;
  float XN, XNMOD=100.;
  char cfilt[2];

  // --------------- BEGIN --------

  if (argc < 2) { print_simlib_coadd_help();  exit(0); }

  parse_args(argc, argv );

  REQUIRE_DOCANA_COADD = true ; // Feb 25 2021

  if ( LTRACE > 0 ) dmp_trace_main("01");

  // open and read header info
  SIMLIB_open_read(); 

  if ( LTRACE > 0 ) dmp_trace_main("02");

  init_summary_info();

  RDSTAT  = 2 ;

  while ( RDSTAT != EOF ) {

  // set pointer used in simlib_tools
    FPLIB = fp_simlib_input ;

    SIMLIB_read(&RDSTAT);  // read next LIBID

    LIBID = SIMLIB_INPUT.LIBID ;

    if ( INPUTS.OPT_SORT_BAND && LIBID>=0 )  { SIMLIB_sort_band(); ; }

    XN = (float)NLIBID_FOUND;
    if ( fmodf(XN,XNMOD) == 0.0 && LIBID >= 0 ) 
      { printf("  Process LIBID %4d \n", LIBID); }

    sprintf(BANNER,"Loop 03: LIBID=%d", LIBID);
    if ( LTRACE > 0 ) dmp_trace_main(BANNER);

    if ( LIBID < 0   ) continue ;

    // process valid LIBIDs
    SIMLIB_coadd();

    // apply MINOBS requirement to output libid

    NOBS = SIMLIB_OUTPUT.NOBS ;
    if ( NOBS < INPUTS.MINOBS_ACCEPT ) {
      printf("\t Skipping LIBID %d : only %d compact exposures. \n", 
	     LIBID, NOBS );
      continue ;
    }

    sprintf(BANNER,"Loop 04: LIBID=%d", LIBID);
    if ( LTRACE > 0 ) dmp_trace_main(BANNER);

      // write to output [compact] SIMLIB file

    FPLIB = fp_simlib_output ;
    simlib_add_header(0 
		      ,SIMLIB_OUTPUT.LIBID
		      ,SIMLIB_OUTPUT.NOBS
		      ,SIMLIB_OUTPUT.FIELDNAME
		      ,SIMLIB_OUTPUT.INFO_HEAD
		      );

    sprintf(BANNER,"Loop 05: LIBID=%d", LIBID);
    if ( LTRACE > 0 ) dmp_trace_main(BANNER);
    
    for ( obs=0; obs < SIMLIB_OUTPUT.NOBS; obs++ ) {
      update_summary_info(obs);
      simlib_add_mjd(
		     1            // 1=>search info;  2=> template info
		     ,SIMLIB_OUTPUT.INFO_OBS[obs]
		     ,SIMLIB_OUTPUT.STRING_IDEXPT[obs]
		     ,SIMLIB_OUTPUT.BAND[obs]
		     );
    }
    
    sprintf(BANNER,"Loop 06: LIBID=%d", LIBID);
    if ( LTRACE > 0 ) dmp_trace_main(BANNER);

    // leave end-of-LIBID marker
    simlib_add_header(-1
		      ,SIMLIB_OUTPUT.LIBID
		      ,SIMLIB_OUTPUT.NOBS
		      ,SIMLIB_OUTPUT.FIELDNAME
		      ,SIMLIB_OUTPUT.INFO_HEAD
		      );

    NLIBID_COADD++;

    sprintf(BANNER,"Loop 07: LIBID=%d", LIBID);
    if ( LTRACE > 0 ) dmp_trace_main(BANNER);

  } // end of while loop

  fclose(fp_simlib_input);

  if ( LTRACE > 0 ) dmp_trace_main("after fclose");

  FPLIB = fp_simlib_output ;
  simlib_close_write();

  insert_NLIBID();

  print_summary_info();

}  // end of main


// ********************************************
void  print_simlib_coadd_help(void) {

  // Created Sep 1 2022

  static char *help[] = {

    "",
    "\t ***** simlib_coadd.exe help menu *****",
    "",
    "# Usage: ",
    "    simlib_coadd <simlib_file> <optional args> ",
    "#     (there is no config file)",
    "",
    "# Optional arguments: ",
    "",
    "MWEBV               # compute MWEBV from Schlagel maps",
    "--TDIF <TDIF>       # coadd within TDIF, days (default=0.4)", 
    "--MINOBS <MINOBS>  # reject LIBIDs with few than MINOBS (default=2)",
    "",
    "--MJD_MIN <MJD_MIN>      # select MJDs above MJD_MIN",
    "--MJD_MAX <MJD_MAX>      # select MJDs below MJD_MAX",
    "--LIBID_MAX <LIBID_MAX>    ",
    "--LIBID_MIN <LIBID_MIN>    ",
    "",
    "MJD_DMP     # print MJDs to stdout",
    "SUM         # compute ZP to add fluxes; default ",
    "AVG         # compute ZP to average fluxes instead of add",
    "SNLS        # original hack to add Vega zeropts for Astier06",
    ""
    "SORT_BAND   # sort by band before coadd",
    "#   (e.g., g,r,i,g,r,i -> gg,rr,ii so that each band is coadded)",
    0
  };

  int i;
  // ---------- BEGIN ---------
  for (i = 0; help[i] != 0; i++)
    { printf ("%s\n", help[i]); }
  
  return;

} // end print_simlib_coadd_help

// ********************************************
void  parse_args(int argc, char **argv) {

  int i, i1, N;
  char *ptrhead, copt[20];
  char fnam[] = "parse_args" ;

  // --------------- BEGIN --------------

  LTRACE      = 0 ;
  INPUTS.OPT_MJD_DMP = 0 ;
  INPUTS.OPT_AVG     = 0 ;
  INPUTS.OPT_SUM     = 1 ;
  sprintf(copt,"SUMMED") ;
  INPUTS.OPT_SNLS    = 0 ;
  INPUTS.OPT_MWEBV   = 0 ;
  INPUTS.OPT_SORT_BAND = 0 ;

  // combine consecutive exposures in same filter 
  // within this time-diff (days)
  INPUTS.MAXTDIF_COMBINE  = 0.4  ; 

  INPUTS.MJD_MIN       = 20000. ;
  INPUTS.MJD_MAX       = 80000. ;
  INPUTS.LIBID_MAX     = 1000000 ;
  INPUTS.LIBID_MIN     = 0 ;
  INPUTS.MINOBS_ACCEPT = 3 ;

  sprintf(SIMLIB_INPUT.FILE, "UNKNOWN" );

  if ( argc >= 2 ) {
    sscanf( argv[1], "%s", SIMLIB_INPUT.FILE );
    sprintf(SIMLIB_OUTPUT.FILE, "%s.COADD", SIMLIB_INPUT.FILE );
  }
  else {
    sprintf(c1err, "Must give SIMLIB file year as argument." ); 
    c2err[0] = 0;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  } 


  for ( i = 1; i < argc; i++ ) {
    i1 = i + 1 ;

    if ( strcmp(argv[i],"--MJD_MIN") == 0 ) 
      { sscanf ( argv[i1], "%le", &INPUTS.MJD_MIN ); }

    if ( strcmp(argv[i],"--MJD_MAX") == 0 ) 
      { sscanf ( argv[i1], "%le", &INPUTS.MJD_MAX ); }

    if ( strcmp(argv[i],"--LIBID_MAX") == 0 ) 
      { sscanf ( argv[i1], "%d", &INPUTS.LIBID_MAX ); }

    if ( strcmp(argv[i],"--LIBID_MIN") == 0 ) 
      { sscanf ( argv[i1], "%d", &INPUTS.LIBID_MIN ); }

    if ( strcmp(argv[i], "--TDIF" ) == 0 ) 
      { sscanf ( argv[i1], "%le", &INPUTS.MAXTDIF_COMBINE ); }

    if ( strcmp(argv[i], "--MINOBS" ) == 0 ) 
      { sscanf ( argv[i1], "%d", &INPUTS.MINOBS_ACCEPT ); }

    if ( strcmp(argv[i], "MWEBV" ) == 0 ) 
      { INPUTS.OPT_MWEBV = 1; }

    if ( strcmp(argv[i], "SORT_BAND" ) == 0 ) 
      { INPUTS.OPT_SORT_BAND = 1; }

    if ( strcmp(argv[i],"MJD_DMP") == 0 ) 
      { INPUTS.OPT_MJD_DMP = 1;  }

    if ( strcmp(argv[i],"TRACE") == 0 ) 
      { LTRACE = 1;  }

    if ( strcmp(argv[i],"AVG") == 0 ) 
      {  INPUTS.OPT_AVG = 1;   sprintf(copt,"AVERAGED"); }

    if ( strcmp(argv[i],"SUM") == 0 ) 
      { INPUTS.OPT_SUM = 1; sprintf(copt,"SUMMED"); }

    if ( strcmp(argv[i],"SNLS") == 0 ) 
      { INPUTS.OPT_SNLS = 1;  }

  }

  N=0;

  ptrhead = HEADER_ADD[N]; N++ ;
  sprintf(ptrhead,"This coadd SIMLIB was created with command");

  ptrhead = HEADER_ADD[N]; N++;
  sprintf(ptrhead,"   '"  );
  for ( i = 0; i < argc; i++ )  
    { strcat(ptrhead," "); strcat(ptrhead,argv[i]); }

  strcat(ptrhead," ' ");

  // now parse what's happening

  //  ptrhead = HEADER_ADD[N]; N++ ;
  //  sprintf(ptrhead," Select  MJD   between %9.2f and %9.2f " );

  ptrhead = HEADER_ADD[N]; N++ ;
  sprintf(ptrhead,"and the following co-add cuts/options:" );

  ptrhead = HEADER_ADD[N]; N++ ;
  sprintf(ptrhead,"   + Select LIBID between %d and %d", 
    INPUTS.LIBID_MIN, INPUTS.LIBID_MAX );

 ptrhead = HEADER_ADD[N]; N++ ;
 sprintf(ptrhead,"   + Reject LIBID with < %d exposures", 
    INPUTS.MINOBS_ACCEPT);

 ptrhead = HEADER_ADD[N]; N++ ;
 sprintf(ptrhead,"   + Combine consecutive exposures within %4.3f days", 
    INPUTS.MAXTDIF_COMBINE);

 if ( INPUTS.OPT_SNLS == 1 ) {
   ptrhead = HEADER_ADD[N]; N++ ;
   sprintf(ptrhead," USE SNLS OPTION => ADD VEGA OFFSETS TO ZPTs " );
   INPUTS.OPT_AVG = 1;
   INPUTS.OPT_SUM = 0;
 }

  if ( INPUTS.OPT_MWEBV > 0 ) {
    ptrhead = HEADER_ADD[N]; N++ ;
    sprintf(ptrhead,"   + Get MWEBV from Schlagel dust maps "
	    "(reject MWEBV > %3.1f). \n", MWEBV_MAX );
  }

  ptrhead = HEADER_ADD[N]; N++ ;
  sprintf(ptrhead,"   + Multiple exposures are '%s' ", copt);

  if ( N >= MXLINE_HEADER ) {
    sprintf(c1err,"%d HEADER_ADD lines exceeds array bound of %d.", 
	   N, MXLINE_HEADER );
    sprintf(c2err,"See MXLINE_HEADER in *.h files");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  NLINE_HEADER_ADD = N;

  // dump compact-header info to stdout (screen)
  for ( i=0; i < N ; i++ ) 
    { printf("HEADER_ADD: %s\n", HEADER_ADD[i] ); }

  return ;

} // end of function parse_args



// *************************************
void SIMLIB_open_read(void) {

  // Copied/modified from snlc_sim.c
  // Open input SIMLIB file and read header.
  // Jan 7 2021: use snana_open to allow reading gzipped SIMLIB.

  char cline[MXPATHLEN], c_get[MXPATHLEN], c_tmp[MXPATHLEN];
  char clast[MXPATHLEN], key[MXPATHLEN];
  char fullName[MXPATHLEN] ;

  bool FOUND_COMMENT = false, FOUND_DOCANA = false, FOUND_FILTERS=false ;
  int  READHEAD, i, gzipFlag, iwd, NWD, add ;
  int  langC = 0;
  char fnam[] = "SIMLIB_open_read" ;
  // ---------------- BEGIN --------------

  printf("\n SIMLIB_open_read(): \n ");

  ENVreplace(SIMLIB_INPUT.FILE, fnam, 1);

  int OPTMASK = 1;  // 1=verbose

  fp_simlib_input = snana_openTextFile (OPTMASK, "", SIMLIB_INPUT.FILE,
					fullName,  &gzipFlag ); 
  if ( !fp_simlib_input ) {
    sprintf(c1err,"cannot open input simlib file: ");
    sprintf(c2err," '%s' ", SIMLIB_INPUT.FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  rewind(fp_simlib_input);

  printf("\t Opened %s \n", SIMLIB_INPUT.FILE );
  fflush(stdout);

  // read header keywords. Stop when we reach "BEGIN"

  NLINE_HEADER = 0;
  READHEAD   = 1;
  PSF_NEA_UNIT = false;

  while( READHEAD > 0 ) {

    fgets(cline, 80, fp_simlib_input) ;
    printf(" Found header line: %s", cline );     fflush(stdout);

    sprintf( HEADER[NLINE_HEADER], "%s", cline );
    NLINE_HEADER++ ;

    if ( strlen(cline) < 2 ) { continue; }

    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,cline, fnam);
    iwd=0; get_PARSE_WORD(langC, iwd, key); 

    if ( strcmp(key,KEYNAME_DOCANA_REQUIRED) == 0 )  
      { FOUND_DOCANA = true; }

    // when we find DOCUMENTATION_END key, insert a few things above END key
    if ( strcmp(key,KEYNAME2_DOCANA_REQUIRED) == 0 )  {
      NLINE_HEADER-- ;
      sprintf(HEADER[NLINE_HEADER], "    SIMLIB_COADD:\n" );
      NLINE_HEADER++; 
      for (add = 0; add < NLINE_HEADER_ADD; add++ ) {
	sprintf(HEADER[NLINE_HEADER], "    - %s\n", HEADER_ADD[add] );
	NLINE_HEADER++; 
      }
      sprintf(HEADER[NLINE_HEADER], "%s", KEYNAME2_DOCANA_REQUIRED) ;
      NLINE_HEADER++; 

      sprintf(HEADER[NLINE_HEADER], "\n\n") ;
      NLINE_HEADER++; 
    }
    
    if ( strcmp(key,"COMMENT:") == 0 )  { 
      FOUND_COMMENT = true ; 
      if ( REQUIRE_DOCANA_COADD ) {
	sprintf(c1err,"COMMENT: keys no longer allowed." ) ;
	sprintf(c2err,"Add %s YAML block before SURVEY key.",
		KEYNAME_DOCANA_REQUIRED );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }      
    }


    // look for header end with "BEGIN" or "#"

    // remove pre-existing NLIBID key because a new NLIBID key
    // is added at the end.
    if ( strcmp(key,"NLIBID:") == 0 )
      { NLINE_HEADER--; continue ; }

    if ( strcmp(key,"PSF_UNIT:") == 0 ) {
      iwd++ ; get_PARSE_WORD(langC, iwd, c_tmp ); 
      if ( strcmp(c_tmp,"NEA_PIXEL") == 0 ) { PSF_NEA_UNIT = true; }
    }

    if ( strcmp(key,"BEGIN")   == 0 ) { READHEAD = 0 ; } // Feb 2021

    if ( strcmp(key,"LIBID:")  == 0 ) { READHEAD = 0; }


    // look for FILTERS keyword

    while ( !FOUND_FILTERS && iwd < NWD-1 ) {
      iwd++ ; get_PARSE_WORD(langC, iwd, key);
      if (strcmp(key,"FILTERS:") == 0 ) {
	iwd++ ; get_PARSE_WORD(langC, iwd, SIMLIB_FILTERS );
	printf(" Found SIMLIB_FILTERS: %s \n", SIMLIB_FILTERS );
	fflush(stdout);
	FOUND_FILTERS = true;
      }
    }

  } // end of while 


  // -----------------
  if ( REQUIRE_DOCANA_COADD && !FOUND_DOCANA ) {
    sprintf(c1err,"Missing required %s block.", 
	    KEYNAME_DOCANA_REQUIRED) ;
    sprintf(c2err,"Add %s YAML block before SURVEY key.",
	    KEYNAME_DOCANA_REQUIRED) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);

  }

  if ( !FOUND_FILTERS ) {
    sprintf(c1err,"Missing required FILTERS: key.");
    sprintf(c2err,"Check global header at top of file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);

  }

  // open output file manually instead of using
  // simlib_open_write() since we just have 80-char lines
  // instead of the detailed info needed by the open function
  // Dump header info to output file.

  if ( (fp_simlib_output = fopen(SIMLIB_OUTPUT.FILE, "wt")) == NULL ) {   
    sprintf(c1err,"Cannot open output simlib file:");
    sprintf(c2err," '%s'", SIMLIB_OUTPUT.FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }


  printf("\t Opened %s \n", SIMLIB_OUTPUT.FILE );

  for ( i=0; i < NLINE_HEADER; i++ ) {
    fprintf(fp_simlib_output, "%s", HEADER[i] );
  }

  // add a few  blank spaces
  fprintf(fp_simlib_output, "\n\n\n" );
  fflush(fp_simlib_output);

  return ;

}  // end of SIMLIB_open_read


// ************************************
void SIMLIB_read(int *RDSTAT) {

  /* Copied/modified from snlc_sim.c    
    Read next LIBID in simlib. Fill SIMLIB structure.
    Note that fp_simlib_input has already been opened

    Dec 6 2021: set INFO_HEAD[IPAR_NEA_UNIT] 
  */

  char  c_get[80], c_tmp[80],cfilt[4], STRING_IDEXPT[20] ;
  int   IDEXPT, NEXPOSE;

  int i, o, NOBS_READ, NOBS_ACCEPT, LIBID ;
  int OPTLINE, NOBS_EXPECT, ENDLIB, OKLIBID, iwd, NWD    ;

  double MJD, CCDGAIN, CCDNOISE, SKYSIG, NEA, PSF[3], ZPT[2] ;
  double MAG, RA, DEC, XMW[20], MWEBV    ;

  char  fnam[] = "SIMLIB_read"  ;

  // ---------------- BEGIN ------------------

  NOBS_READ    =  0  ;
  NOBS_ACCEPT  =  0  ;
  ENDLIB       =  0  ;
  *RDSTAT      =  2  ;

  SIMLIB_INPUT.LIBID = -9 ;

  while( *RDSTAT != EOF  ) {

    *RDSTAT = fscanf(fp_simlib_input, "%s", c_get) ;

    if ( strcmp(c_get,"END_OF_SIMLIB:") == 0 ) {    }

    if ( strcmp(c_get,"LIBID:") == 0 ) {
      NLIBID_FOUND++ ;
      readint ( fp_simlib_input, 1, &LIBID );
      SIMLIB_INPUT.LIBID = LIBID;

      // init info for this LIBID

      SIMLIB_INPUT.NOBS = 0 ;
      for ( i=0; i < NPAR_HEAD ; i++ )
	SIMLIB_INPUT.INFO_HEAD[i] = -9.0 ;
    }

    if ( strcmp(c_get,"NOBS:") == 0 ) {
      readint ( fp_simlib_input, 1, &NOBS_EXPECT );
      SIMLIB_INPUT.NOBS = NOBS_EXPECT ;
    }

    if ( strcmp(c_get,"RA:") == 0 ) 
      { readfloat(fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_RA] ); }

    if ( strcmp(c_get,"DECL:") == 0 || strcmp(c_get,"DEC")==0 ) 
      { readfloat(fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_DEC] ); }

    if ( strcmp(c_get,"MWEBV:") == 0  ) 
      { readfloat(fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_MWEBV] );}

    if ( strcmp(c_get,"PIXSIZE:") == 0 ) 
      { readfloat(fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_PIXSIZE]);}

    if ( strcmp(c_get,"REDSHIFT:") == 0 ) 
      { readfloat(fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_Z] );}

    if ( strcmp(c_get,"PEAKMJD:") == 0 ) 
      { readfloat(fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_PKMJD] );}

    if ( strcmp(c_get,"FIELD:") == 0 ) 
      { fscanf(fp_simlib_input,"%s", SIMLIB_INPUT.FIELDNAME ); }

    // Search for MJD-epoch info

    if ( strcmp(c_get,"S:") == 0 ) 
      { OPTLINE = 1 ; }
    else if ( strcmp(c_get,"T:") == 0 ) 
      { OPTLINE = 2 ; }
    else
      { OPTLINE = 0 ; }

    if ( OPTLINE > 0 ) {

      NOBS_READ++ ;

      // read line into temporary variables.
      readdouble(fp_simlib_input, 1, &MJD );
      readchar(fp_simlib_input, STRING_IDEXPT );
      parse_SIMLIB_IDplusNEXPOSE(STRING_IDEXPT, &IDEXPT, &NEXPOSE);

      fscanf(fp_simlib_input, "%s", cfilt );
      readdouble(fp_simlib_input, 1, &CCDGAIN  );
      readdouble(fp_simlib_input, 1, &CCDNOISE );
      readdouble(fp_simlib_input, 1, &SKYSIG );

      if ( PSF_NEA_UNIT ) 
	{ readdouble(fp_simlib_input, 1, &NEA ); }
      else 
	{ readdouble(fp_simlib_input, 3, PSF ); }

      readdouble(fp_simlib_input, 2, ZPT );

      if ( OPTLINE == 1 ) {
	readdouble ( fp_simlib_input, 1, &MAG );
	if ( MAG < 5.0 ) { MAG = 99.0 ; }
      }
      else {
	MAG = 99.0 ; 
      } 

      NOBS_ACCEPT++ ;
      SIMLIB_INPUT.NOBS_ACCEPT = NOBS_ACCEPT ;
      o = NOBS_ACCEPT-1;  // o index starts at zero (Mar 7 2022)

      if ( NOBS_ACCEPT >= MXMJD ) {
	sprintf(c1err,"NOBS=%d exceeds array bound at LIBID=%d. ", 
	       NOBS_ACCEPT, LIBID);
	sprintf(c2err,"MXMJD = %d", MXMJD);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      sprintf(SIMLIB_INPUT.STRING_IDEXPT[o], "%s", STRING_IDEXPT) ;
      SIMLIB_INPUT.IDEXPT[o]         = IDEXPT ;
      SIMLIB_INPUT.NEXPOSE_IDEXPT[o] = NEXPOSE ;

      sprintf(SIMLIB_INPUT.BAND[o], "%s", cfilt);

      SIMLIB_INPUT.INFO_OBS[o][IPAR_MJD]       = MJD ;
      SIMLIB_INPUT.INFO_OBS[o][IPAR_CCDGAIN]   = CCDGAIN ;
      SIMLIB_INPUT.INFO_OBS[o][IPAR_CCDNOISE]  = CCDNOISE ;
      SIMLIB_INPUT.INFO_OBS[o][IPAR_SKYSIG]    = SKYSIG   ;

      if ( PSF_NEA_UNIT ) {
	SIMLIB_INPUT.INFO_OBS[o][IPAR_PSF0]    = NEA ;
	SIMLIB_INPUT.INFO_OBS[o][IPAR_PSF0+1]  = -999.0 ;
	SIMLIB_INPUT.INFO_OBS[o][IPAR_PSF0+2]  = -999.0 ;

      }
      else {
	SIMLIB_INPUT.INFO_OBS[o][IPAR_PSF0]    = PSF[0] ;
	SIMLIB_INPUT.INFO_OBS[o][IPAR_PSF0+1]  = PSF[1] ;
	SIMLIB_INPUT.INFO_OBS[o][IPAR_PSF0+2]  = PSF[2] ;
      }

      SIMLIB_INPUT.INFO_OBS[o][IPAR_ZPT0]    = ZPT[0] ;
      SIMLIB_INPUT.INFO_OBS[o][IPAR_ZPT0+1]  = ZPT[1] ;
      SIMLIB_INPUT.INFO_OBS[o][IPAR_MAG]     = MAG ;

      if ( MJD < INPUTS.MJD_MIN ||  MJD > INPUTS.MJD_MAX ) {
	NOBS_ACCEPT-- ;
	SIMLIB_INPUT.NOBS_ACCEPT = NOBS_ACCEPT ;
      }

    } // end of OPTLINE if-block


    if ( strcmp(c_get,"END_LIBID:") == 0         ) ENDLIB = 1;
    if ( strcmp(c_get,"#") == 0 && NOBS_READ > 3 ) ENDLIB = 1;

    if ( ENDLIB == 1 ) {

      if ( NOBS_READ != NOBS_EXPECT && INPUTS.OPT_MJD_DMP == 0 ) {
	sprintf(c1err,"NOBS_READ = %d, but expected NOBS=%d (LIBID=%d)",
	       NOBS_READ, NOBS_EXPECT, LIBID);
	c2err[0] = 0 ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      OKLIBID = -9 ;
      if ( LIBID >= INPUTS.LIBID_MIN && LIBID <= INPUTS.LIBID_MAX ) 
	{ OKLIBID = 1; }

      //      printf("  LIBID=%d  OKLIBID = %d \n", LIBID, OKLIBID );
      if ( OKLIBID < 0 ) { SIMLIB_INPUT.LIBID = -9 ; }


      // if there are too few accepted exposures, set flag to skip LIBID

      if ( NOBS_ACCEPT < INPUTS.MINOBS_ACCEPT ) {
	printf("\t Skipping LIBID %d : only %d exposures. \n", 
	       LIBID, NOBS_ACCEPT);
	fflush(stdout);
	SIMLIB_INPUT.LIBID = -9;
      }

      // check option to compute MWEBV from Schlagel maps
      if ( INPUTS.OPT_MWEBV > 0 ) {
	RA   = SIMLIB_INPUT.INFO_HEAD[IPAR_RA] ;
	DEC  = SIMLIB_INPUT.INFO_HEAD[IPAR_DEC] ;
	MWgaldust(RA, DEC, &XMW[1], &MWEBV );  // returns XMW & MWEBV
	SIMLIB_INPUT.INFO_HEAD[IPAR_MWEBV] = (float)MWEBV ;
      }

      // skip really larger MW extinctions
      MWEBV = SIMLIB_INPUT.INFO_HEAD[IPAR_MWEBV];
      if ( MWEBV > MWEBV_MAX ) {
	printf("\t Skipping LIBID %d : MWEBV=%6.1f \n", LIBID, MWEBV );
	fflush(stdout);
	SIMLIB_INPUT.LIBID = -9;
      }
      return ;
    }

  }  // end of while fscanf loop

  SIMLIB_INPUT.INFO_HEAD[IPAR_NEA_UNIT] = (float)PSF_NEA_UNIT ;

  return ;
  
} // end of function SIMLIB_read


// *******************************************
void SIMLIB_sort_band(void) {

  // Created Mar 7 2022
  // Sort SIMLIB by band so that co-add includes
  // non-sequential bands.

  // xxx mark delete  SIMLIB_CONTENTS_DEF SIMLIB_INPUT_TEMP;
  int  NOBS_orig = SIMLIB_INPUT.NOBS;
  int  ifilt, NFILT, o, NOBS_copy=0 ;
  int  LDMP = 0;
  char band[2];
  char fnam[] = "SIMLIB_sort_band" ;

  // ------------ BEGIN --------------

  NFILT = strlen(SIMLIB_FILTERS);

  if ( LDMP ) {
    printf("xxx %s Sort obs by band NFILT=%d for %s  LIBID=%d\n", 
	   fnam, NFILT, SIMLIB_FILTERS, SIMLIB_INPUT.LIBID);
  }

  copy_SIMLIB_CONTENTS(&SIMLIB_INPUT, &SIMLIB_INPUT_TEMP);

  for(ifilt=0; ifilt < NFILT; ifilt++ ) {
    sprintf(band, "%c", SIMLIB_FILTERS[ifilt] );

    if ( LDMP ) 
      { printf("\t xxx collect %s band \n", band); fflush(stdout); }

    for(o=0 ; o < NOBS_orig; o++ ) {
      if (strcmp(band,SIMLIB_INPUT_TEMP.BAND[o]) == 0 ) {
	copy_SIMLIB_CONTENTS_OBS(&SIMLIB_INPUT_TEMP, 
				 &SIMLIB_INPUT, o, NOBS_copy);
	NOBS_copy++ ;
      }

    } // end obs loop
  } // end ifilt

  if ( LDMP ) {
    printf("\t xxx NOBS[orig,copy] = %d, %d \n", NOBS_orig, NOBS_copy);
  }
  return;

} // end SIMLIB_sort_band

// *********************
void SIMLIB_coadd(void) {

  //  transfer SIMLIB_input -> SIMLIB_OUTPUT structure,
  //  where the output is in compact form.

  // May 20, 2009: special fix for taking MJD average without roundoff error
  // Jun 20, 2017: MJD -> double instead of float
  // Jan 07, 2021: sum NEXPOSE and write proper IDEXPT string

  int  i, j, obs, ipar, NOBS_IN, NMEASURE, OVPFILT, IDEXPT, NEXPOSE ;
  int  OBSMIN[MXMJD], OBSMAX[MXMJD]    ;
  char *cfilt, *cfilt_last ;
  double MJD, MJD_LAST, MJD_DIF, XIN, XSUM, XN, XNOPT, ARG, ZPTOFF;
  double *PTR_INFO_INPUT, *PTR_INFO_OUTPUT  ;

  char fnam[] = "SIMLIB_coadd";

  // ------------- BEGIN ----------


  // transfer LIBID & header info without any changes; 

  SIMLIB_OUTPUT.LIBID = SIMLIB_INPUT.LIBID ;

  for ( i=0; i<NPAR_HEAD; i++ ) {
    SIMLIB_OUTPUT.INFO_HEAD[i] = SIMLIB_INPUT.INFO_HEAD[i] ;
  }
  sprintf(SIMLIB_OUTPUT.FIELDNAME, "%s", SIMLIB_INPUT.FIELDNAME ); 

  // ----------------------
  // first loop through and identify MJD-ranges to combine.

  NOBS_IN  = SIMLIB_INPUT.NOBS_ACCEPT ;
  MJD_LAST = -9.0 ;
  NMEASURE = 0 ;

  obs = 0;
  cfilt      = SIMLIB_INPUT.BAND[obs] ;
  cfilt_last = SIMLIB_INPUT.BAND[obs] ;

  for ( obs=0; obs < NOBS_IN; obs++ ) {
    MJD     = SIMLIB_INPUT.INFO_OBS[obs][IPAR_MJD] ;
    MJD_DIF = fabs(MJD - MJD_LAST);
    cfilt   = SIMLIB_INPUT.BAND[obs] ;
    OVPFILT = strcmp(cfilt,cfilt_last);  

    if ( MJD_DIF < INPUTS.MAXTDIF_COMBINE && OVPFILT == 0 ) {
      //  different MJD, same filter -> update OBSMAX
      OBSMAX[NMEASURE-1] = obs ;  
    } 
    else {
      // next band
      OBSMIN[NMEASURE] = obs ;
      OBSMAX[NMEASURE] = obs ;
      NMEASURE++ ;   
    }

    cfilt_last = cfilt ;
    MJD_LAST   = MJD ;

  }  // end of obs loop


  // -----------------
  // Now loop through and combine exposures and take appropriate
  // averages for SIMLIB_OUTPUT structure

  SIMLIB_OUTPUT.NOBS  = NMEASURE ;
  for ( i = 0; i < NMEASURE; i++ ) {

    // for filter, copy element from 1st exposure to OUTPUT measurement,
    obs = OBSMIN[i] ;

    cfilt = SIMLIB_INPUT.BAND[obs] ;
    sprintf(SIMLIB_OUTPUT.BAND[i], "%s", cfilt);

    NEXPOSE = 0;
    IDEXPT  = SIMLIB_INPUT.IDEXPT[obs] ;

    // setup pointer to output INFO array
    PTR_INFO_OUTPUT = &SIMLIB_OUTPUT.INFO_OBS[i][0] ;

    // init output INFO array 
    for ( ipar=0; ipar < NPAR_OBS; ipar++ ) 
      { PTR_INFO_OUTPUT[ipar] = 0.0 ; }

    // take sums or sums-of-squares over INPUT 'obs' observations

    XN = 0.0 ;

    for ( obs = OBSMIN[i]; obs <= OBSMAX[i]; obs++ ) {

      NEXPOSE += SIMLIB_INPUT.NEXPOSE_IDEXPT[obs];

      XN += 1.0 ;  // number of exposures for this measurement.

      PTR_INFO_INPUT   = &SIMLIB_INPUT.INFO_OBS[obs][0] ;

      XIN = PTR_INFO_INPUT[IPAR_MJD];
      SIMLIB_OUTPUT.INFO_OBS[i][IPAR_MJD] += XIN ;

      XIN = PTR_INFO_INPUT[IPAR_CCDGAIN] ;
      PTR_INFO_OUTPUT[IPAR_CCDGAIN] += XIN ;

      XIN = PTR_INFO_INPUT[IPAR_CCDNOISE] ;
      PTR_INFO_OUTPUT[IPAR_CCDNOISE] += (XIN * XIN) ;

      XIN = PTR_INFO_INPUT[IPAR_SKYSIG] ;
      PTR_INFO_OUTPUT[IPAR_SKYSIG] += (XIN * XIN) ;

      // - - - -
      for( j=0; j < 3; j++ ) {
	XIN = PTR_INFO_INPUT[IPAR_PSF0+j] ;
	PTR_INFO_OUTPUT[IPAR_PSF0+j] += XIN ;
      }

      // - - - 

      XIN = PTR_INFO_INPUT[IPAR_ZPT0] ;
      ARG = 0.4 * (XIN - 25.0) ;
      PTR_INFO_OUTPUT[IPAR_ZPT0] += powf(TEN,ARG); 

      XIN = PTR_INFO_INPUT[IPAR_ZPT0+1] ;
      PTR_INFO_OUTPUT[IPAR_ZPT0+1] += XIN ;

      XIN = PTR_INFO_INPUT[IPAR_MAG] ;
      PTR_INFO_OUTPUT[IPAR_MAG] += XIN ;

    } // end of obs loop

    // now divide by NMEASURE, take sqrt, etc ... to get
    // appropriate result for single measure
    if ( XN == 0.0 ) {
      sprintf(c1err,"Nexposure=0 for Meaure=%d OBS=%d-%d, filt=%s ",
	     i, OBSMIN[i], OBSMAX[i], cfilt );
      sprintf(c2err," MJD = %f", SIMLIB_INPUT.INFO_OBS[OBSMIN[i]][0] ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    sprintf(SIMLIB_OUTPUT.STRING_IDEXPT[i], "%d*%d", IDEXPT, NEXPOSE);

    if ( INPUTS.OPT_SUM > 0 ) 
      { XNOPT = 1.0 ; }
    else
      { XNOPT = XN ; }

    SIMLIB_OUTPUT.INFO_OBS[i][IPAR_MJD] /= XN ;

    PTR_INFO_OUTPUT[IPAR_CCDGAIN] *= (XNOPT/XN) ;

    XSUM = PTR_INFO_OUTPUT[IPAR_CCDNOISE] ;
    PTR_INFO_OUTPUT[IPAR_CCDNOISE] = sqrtf(XSUM) / XNOPT;

    XSUM = PTR_INFO_OUTPUT[IPAR_SKYSIG] ;
    PTR_INFO_OUTPUT[IPAR_SKYSIG] = sqrtf(XSUM) / XNOPT;

    for(j=0; j < 3; j++ ) 
      { PTR_INFO_OUTPUT[IPAR_PSF0+j] /= XN ; }

    XSUM = PTR_INFO_OUTPUT[IPAR_ZPT0] ;
    PTR_INFO_OUTPUT[IPAR_ZPT0] = 2.5 * log10f(XSUM/XNOPT) + 25.0 ;

    PTR_INFO_OUTPUT[IPAR_ZPT0+1] /= XN ;
    
    if ( INPUTS.OPT_SNLS == 1 ) {
      ZPTOFF = ZPTOFF_SNLS(cfilt);           // hard-coded hack
      PTR_INFO_OUTPUT[IPAR_ZPT0+0] += ZPTOFF;
      PTR_INFO_OUTPUT[IPAR_SKYSIG] *= 1.5 ;  // special hack
    }


    PTR_INFO_OUTPUT[IPAR_MAG] /= XN ;

    if ( INPUTS.OPT_MJD_DMP == 1 ) {
      MJD = SIMLIB_OUTPUT.INFO_OBS[i][IPAR_MJD];
      printf(" %f \n", MJD );      fflush(stdout);
    }

  }  // end of i-loop over NMEASURE

  return ;

} // end of SIMLIB_coadd


// ************************************************************
void copy_SIMLIB_CONTENTS_OBS(SIMLIB_CONTENTS_DEF *CONTENTS_INP,
			      SIMLIB_CONTENTS_DEF *CONTENTS_OUT, 
			      int o_inp,  int o_out ) {

  // Created Mar 2022

  int ipar;
  // --------- BEGIN -----------

  sprintf(CONTENTS_OUT->BAND[o_out], "%s", 
	  CONTENTS_INP->BAND[o_inp]);
  sprintf(CONTENTS_OUT->STRING_IDEXPT[o_out], "%s", 
	  CONTENTS_INP->STRING_IDEXPT[o_inp]);
  
  CONTENTS_OUT->IFILT[o_out] = CONTENTS_INP->IFILT[o_inp];
  CONTENTS_OUT->NEXPOSE_IDEXPT[o_out] = CONTENTS_INP->NEXPOSE_IDEXPT[o_inp];
  
  for(ipar=0; ipar < NPAR_OBS; ipar++ ) {  
    CONTENTS_OUT->INFO_OBS[o_out][ipar] = 
      CONTENTS_INP->INFO_OBS[o_inp][ipar]; 
  }
  
  return;

} // end copy_SIMLIB_CONTENTS_OBS

// ******************************************
void copy_SIMLIB_CONTENTS(SIMLIB_CONTENTS_DEF *CONTENTS_INP,
			  SIMLIB_CONTENTS_DEF *CONTENTS_OUT ) {

  // Created Mar 7 2022
  // Copy contents of CONTENTS_INP into CONTENTS_OUT.

  int NOBS = CONTENTS_INP->NOBS;
  int ipar, obs;
  char fnam[] = "copy_SIMLIB_CONTENTS";

  // ------------ BEGIN -------------

  sprintf(CONTENTS_OUT->FILE,      "%s", CONTENTS_INP->FILE);
  sprintf(CONTENTS_OUT->FIELDNAME, "%s", CONTENTS_INP->FIELDNAME);

  CONTENTS_OUT->NOBS        = CONTENTS_INP->NOBS;
  CONTENTS_OUT->NOBS_ACCEPT = CONTENTS_INP->NOBS_ACCEPT ;
  CONTENTS_OUT->LIBID       = CONTENTS_INP->LIBID;

  for(ipar=0; ipar < NPAR_HEAD; ipar++ ) 
    { CONTENTS_OUT->INFO_HEAD[ipar] = CONTENTS_INP->INFO_HEAD[ipar]; }

  for(obs=0; obs < NOBS; obs++ ) {
    copy_SIMLIB_CONTENTS_OBS(CONTENTS_INP, CONTENTS_OUT, obs,obs);
  }

  return ;

} // end copy_SIMLIB_CONTENTS

// ******************************************
void init_summary_info(void) {

  SIMLIB_INPUT.LIBID = 0 ;
  NLIBID_FOUND = NLIBID_COADD = 0 ;

  SUMMARY_INFO.MJD_MIN = +1.0E9;
  SUMMARY_INFO.MJD_MAX = 0.0;

  SUMMARY_INFO.NOBS_MIN = 9999999;
  SUMMARY_INFO.NOBS_MAX = 0;

  return ;
  
} // end init_var

// ***********************************
void update_summary_info(int obs) {

  // Jun 23 2022;
  // pass obs argument to check all observations for min/max MJD
  // rather than assuming that the first/last obs is min/max MJD.
  // This fix only impacts the printed MJD range; does NOT impact
  // the SIMLIB contents.

  int NOBS = SIMLIB_OUTPUT.NOBS;
  double MJD ;

  if ( obs == 0 ) {
    if ( NOBS < SUMMARY_INFO.NOBS_MIN ) { SUMMARY_INFO.NOBS_MIN=NOBS; }
    if ( NOBS > SUMMARY_INFO.NOBS_MAX ) { SUMMARY_INFO.NOBS_MAX=NOBS; }
  }

  MJD = SIMLIB_OUTPUT.INFO_OBS[obs][IPAR_MJD];
  if ( MJD < SUMMARY_INFO.MJD_MIN ) { SUMMARY_INFO.MJD_MIN = MJD; }
  if ( MJD > SUMMARY_INFO.MJD_MAX ) { SUMMARY_INFO.MJD_MAX = MJD; }
    
  return;
  
} // end UPDATE_SUMMARY_INFO

// ***********************************
void print_summary_info(void) {

  printf("  Coadd NOBS Range: %d to %d \n",
	 SUMMARY_INFO.NOBS_MIN, SUMMARY_INFO.NOBS_MAX );
  printf("  Coadd MJD Range: %.1f to %.1f \n",
	 SUMMARY_INFO.MJD_MIN, SUMMARY_INFO.MJD_MAX );

  printf("\n Done coadding %d LIBIDs (%d read) \n", 
	 NLIBID_COADD, NLIBID_FOUND );
    
  fflush(stdout);
  
  return ;
  
} // end


// ***********************************
double ZPTOFF_SNLS(char *cfilt) {

  // for Astier 2006
  double zptoff;

  zptoff = 0.0 ;
  if ( strcmp(cfilt,"g") == 0 ) 
    zptoff = +0.115 ;
  else if ( strcmp(cfilt,"r") == 0 ) 
    zptoff = -0.13 ;
  else if ( strcmp(cfilt,"i") == 0 ) 
    zptoff = -0.355 ;
  else if ( strcmp(cfilt,"z") == 0 ) 
    zptoff = -0.521 ;

  return zptoff ;

} // end of ZPTOFF_SNLS





// ******************************
void dmp_trace_main(char *string) {
  printf(" >>>>> TRACE_MAIN-%s  <<<<<< \n", string );
} // end of dmp_trace_main


// *******************************
void  insert_NLIBID(void) {

  // insert NLIBID key just SURVEY line.
  // Uses system command with 'sed' utility to 
  // do the string insert.
  //
  // Jan 10 2019: \ -> \\ in sedCmd to avoid warning on Midway2

  char INSERT_LINE[100], sedCmd[400];
  char fnam[] = "insert_NLIBID" ;
  
  // ------------- BEGIN -------------

  sprintf(INSERT_LINE,"NLIBID: %d", NLIBID_COADD);

  printf("\n Insert '%s' after SURVEY line in header.\n", INSERT_LINE);
  fflush(stdout);

  sprintf(sedCmd,"sed -i '/SURVEY/a\\%s' %s", 
	  INSERT_LINE, SIMLIB_OUTPUT.FILE );

  system(sedCmd);


  return ;

} // end insert_NLIBID
