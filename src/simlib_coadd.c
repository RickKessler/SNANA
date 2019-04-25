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
   ZPTAVG    -> 2.5 * log10 [ \sum 10^(0.4*ZPT_i) ] / N
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

    simlib_coadd <simlib_file>  ASTIER06  (idem SNLS, but
                                           use only Astier06 SN)


  History
  ---------

  Nov 15, 2007: add "ASTIER06" option to use exactly those
                SN that were published. See new function
                LOAD_LIBID_ASTIER06(). Note that 2 of the 71
                published SN are missing from the SIMLIB.
 
  Feb 27, 2008: WARNING, use only for SNLS.
                This code combines measurements in the same night
                without checking filters ... works for SNLS because
                their library is ordered by filter, but does NOT 
                work for ESSENCE since their library has R & I
                interlaced.  For now, the multiple ESSENCE exposures
                [within same passband] is so rare that it's not
                worth fixing this code yet ... but should do so in
                the future.
 

  May 14, 2009: fixed above problem so that multiple filters withing
                the same night are processed OK. However, if your
                library alternates  u, g, u, g ... then nothing
                gets combined. Must have uu gg etc ... to combine
                measurements.
                
  May 20, 2009: add TRACE option to trace seg-fault. 
                Fix dumb bug in SIMLIB_compact; set obs=1 for initial cfilt.
                Change default from AVG to SUM

                MXMJD -> 4000 (was 2000)

                SIMLIB_coadd:
                 special fix for taking MJD average without roundoff error

  May 28, 2009: 
     - add option --TDIF <MAXTDIF_COMBINE> 
     - add option --MINOBS <MINOBS_ACCEPT>
     - add option MWEBV to compute MWEBV (if zero) from Schlagel maps

     - add screen-dump comment to header in output simlib; 
       includes exact command to created compact simlib

 Jun 8, 2009: define FIELDNAME as part of header

 Jan 07, 2011: increase crazy limit on SKYSIG from 400 to 1000;
               this change is in simlib_tools.c (not here)
               Also updated commentary above since it confused David.

 Aug 14, 2012: increase MXMJD and  MXLIBID  to 6000 (from 4000 and 5000).
               Needed by Cinabro for LSST SIMLIBS.

 Aug 31, 2012: MXMJD -> 100,000 (for 10 year LSST simlibs)

 Aug 12, 2014: INPUTS.LIBID_MAX -> 100,000 (was 5,000)

 Jun 20 2017: 
    + define double MJD[MXMJD] and stop using float MJD

 Apr 3 2018:
   + long int IDEXPT -> char STRINGID[20] to allow ID*NEXPOSE
     e.g., 123438*2

 Jun 29 2018
   + fix logic appending COMMENT lines to original comments
   + final call to insert_NLIBID to add "NLIBID: <NLIBID>" key
     in global header.

***************************************/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "simlib_tools.c"
#include "sntools.h" 

#define MXMJD     100000    // max MJDs per LIBID
#define MXLIBID   6000
#define MWEBV_MAX  2.0     // reject fields with such large MWEBV
#define MXLINE_HEADER 30   // max lines for simlib header
#define MXCHAR_LINE 120 

// define index params for INFO_HEAD arra
#define NPAR_HEAD  6
#define IPAR_RA       0
#define IPAR_DECL     1
#define IPAR_MWEBV    2 
#define IPAR_PIXSIZE  3
#define IPAR_Z        4
#define IPAR_PKMJD    5 

// define index parameters for INFO_MJD array
#define NPAR_MJD     10
#define IPAR_CCDGAIN  1
#define IPAR_CCDNOISE 2
#define IPAR_SKYSIG   3
#define IPAR_PSF0     4
#define IPAR_ZPT0     7
#define IPAR_MAG      9



// global variables.

char msgerr[100];
//char BANNER[100];

FILE *fp_simlib_input;
FILE *fp_simlib_output;

int  NLINE_HEADER ;        // number of header lines
int  NLINE_HEADER_ADD ;    // number of comment-header lines to add

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

} INPUTS ;

// -----------------------
// xxx mark delete 6.29.2018 int  LIBID_USE[MXLIBID];
int  NLIBID_FOUND;
int  NLIBID_COADD;

char SIMLIB_FILTERS[20];  // from FILTERS: <filtlist> header

int LTRACE;

struct SIMLIB_INPUT {
  char  FILE[120];
  int   LIBID ;
  int   NLIBID ;
  int   NOBS ;
  int   NMJD_ACCEPT ;
  char  FIELDNAME[20];
  float INFO_HEAD[NPAR_HEAD];  // RA,DECL, MWEBV, PIXSIZE, Z, PEAKMJD

  double   MJD[MXMJD];   // June 2018
  char STRINGID[MXMJD][20]; // Apr 2018
  char FILTNAME[MXMJD][2];
  int  IFILT[MXMJD];

  // info is: MJD, CCDGAIN, CCDNOISE, SKYSIG, PSF[0-2], ZPTAVG, ZPTSIG, MAG
  float INFO_MJD[MXMJD][NPAR_MJD];

} SIMLIB_INPUT ;



struct SIMLIB_OUTPUT {
  char  FILE[120];
  int   LIBID ;
  int   NOBS ;
  char  FIELDNAME[20];
  float INFO_HEAD[NPAR_HEAD];  // RA,DECL, MWEBV, PIXSIZE, Z, PEAKMJD

  double   MJD[MXMJD];  //  June 2017
  // xxx mark delete   long int IDEXPT[MXMJD];
  char STRINGID[MXMJD][20]; // April 2018
  char FILTNAME[MXMJD][2];
 
  // info is: MJD, CCDGAIN, CCDNOISE, SKYSIG, PSF[0-2], ZPTAVG, ZPTSIG, MAG
  float INFO_MJD[MXMJD][NPAR_MJD];

} SIMLIB_OUTPUT ;



// declare functions

void  parse_args(int argc, char **argv);
void  SIMLIB_open();
void  SIMLIB_read(int *RDSTAT);
void  SIMLIB_coadd(void);
void  insert_NLIBID(void);

void dmp_trace_main(char *string);

// define SNLS-specific functions
double ZPTOFF_SNLS(char *cfilt);

// define Galactic extiction

void MWgaldust(double RA, double DECL, double *XMW, double *MWEBV ); 

// ****************************************
int main(int argc, char **argv) {

  int LIBID, obs, ifilt, NOBS ;
  int RDSTAT ;
  float XN, XNMOD=100.;
  char cfilt[2];

  // --------------- BEGIN --------

  parse_args(argc, argv );

  if ( LTRACE > 0 ) dmp_trace_main("01");

  // open and read header info
  SIMLIB_open(); 

  if ( LTRACE > 0 ) dmp_trace_main("02");

  SIMLIB_INPUT.LIBID = 0 ;
  NLIBID_FOUND = NLIBID_COADD = 0 ;
  RDSTAT  = 2 ;


  while ( RDSTAT != EOF ) {

  // set pointer used in simlib_tools
    FPLIB = fp_simlib_input ;

    SIMLIB_read(&RDSTAT);  // read next LIBID

    LIBID = SIMLIB_INPUT.LIBID ;

    XN = (float)NLIBID_FOUND;
    if ( fmodf(XN,XNMOD) == 0.0 && LIBID >= 0 ) 
      printf("  Process LIBID %4d \n", LIBID);

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

    for ( obs=1; obs <= SIMLIB_OUTPUT.NOBS; obs++ ) 
      simlib_add_mjd(
		     1            // 1=>search info;  2=> template info
		     ,SIMLIB_OUTPUT.MJD[obs]
		     ,SIMLIB_OUTPUT.STRINGID[obs]
		     ,SIMLIB_OUTPUT.FILTNAME[obs]
		     ,SIMLIB_OUTPUT.INFO_MJD[obs]
		     );

    
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
  simlib_close();

  insert_NLIBID();

  printf("\n Done compacting %d LIBIDs (%d read) \n", 
	 NLIBID_COADD, NLIBID_FOUND );

}  // end of main


// ********************************************
void  parse_args(int argc, char **argv) {

  int i, i1, N;
  char *ptrhead, BUFFER[20], copt[20];

  // --------------- BEGIN --------------

  LTRACE      = 0 ;
  INPUTS.OPT_MJD_DMP = 0 ;
  INPUTS.OPT_AVG     = 0 ;
  INPUTS.OPT_SUM     = 1 ;
  sprintf(copt,"SUMMED") ;
  INPUTS.OPT_SNLS    = 0 ;
  INPUTS.OPT_MWEBV   = 0 ;

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
    printf("\n ERROR: Must give SIMLIB file year as argument. \n" );
    MADABORT();
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

    /* xxx mark delete Jun 29 2018 xxxxxxx
    if ( strcmp(argv[i],"ASTIER06") == 0 ) {
      INPUTS.OPT_SNLS = 1;  
      INPUTS.LIBID_MAX   = -1 ;  INPUTS.LIBID_MIN   = -1;
      LOAD_LIBID_ASTIER06(); 
    }
    xxxxxxxx */
  }

  N=0;
  sprintf(BUFFER,"COMMENT:");

  // first comment is exact user-command
  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"# - - - - - - - - - - - - - - - - - - - - - - - - \n");

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s This coadd SIMLIB was created with command \n", BUFFER );

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s   '", BUFFER  );
  for ( i = 0; i < argc; i++ )  sprintf(ptrhead,"%s %s", ptrhead, argv[i] ); 
  sprintf(ptrhead,"%s ' \n", ptrhead ); 

  // now parse what's happening

  /* xxxx mark delete Jun 29 2018 xxxxxxxxx
  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s Process input  SIMLIB file : %s \n", 
	  BUFFER, SIMLIB_INPUT.FILE );
  // ptrhead,"%s

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s Create  ouptut SIMLIB file : %s \n", 
	  BUFFER, SIMLIB_OUTPUT.FILE );
  xxxxxxxx */

  //  N++;  ptrhead = HEADER_ADD[N];
  //  sprintf(ptrhead,"%s Select  MJD   between %9.2f and %9.2f \n", 

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s and the following co-add cuts/options: \n", BUFFER );

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s   + Select LIBID between %d and %d \n", 
	  BUFFER, INPUTS.LIBID_MIN, INPUTS.LIBID_MAX );

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s   + Reject LIBID with < %d exposures \n", 
	  BUFFER, INPUTS.MINOBS_ACCEPT);

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s   + Combine consecutive exposures within %4.3f days\n", 
	  BUFFER, INPUTS.MAXTDIF_COMBINE);

  if ( INPUTS.OPT_SNLS == 1 ) {
    N++;  ptrhead = HEADER_ADD[N];
    sprintf(ptrhead,"%s USE SNLS OPTION => ADD VEGA OFFSETS TO ZPTs \n",
	    BUFFER );
    INPUTS.OPT_AVG = 1;
    INPUTS.OPT_SUM = 0;
  }

  if ( INPUTS.OPT_MWEBV > 0 ) {
    N++;  ptrhead = HEADER_ADD[N];
    sprintf(ptrhead,"%s   + Get MWEBV from Schlagel dust maps (reject MWEBV > %3.1f). \n", BUFFER, MWEBV_MAX );
  }

  N++;  ptrhead = HEADER_ADD[N];
  sprintf(ptrhead,"%s   + Multiple exposures are '%s' \n", BUFFER, copt);

  if ( N >= MXLINE_HEADER ) {
    printf("\n ERROR: %d HEADER_ADD lines exceeds array bound of %d.\n", 
	   N, MXLINE_HEADER );
    MADABORT();
  }

  NLINE_HEADER_ADD = N;

  // dump compact-header info to stdout (screen)
  for ( i=1; i <=N ; i++ ) 
    { printf("%s", HEADER_ADD[i] ); }

} // end of function parse_args



// *************************************
void SIMLIB_open(void) {

  // Copied/modified from snlc_sim.c
  // Open input SIMLIB file and read header.

  char fnam[20] = "SIMLIB_open" ;
  char cline[80], c_get[80], c_tmp[80], clast[80] ;
  char *ptrtok;

  int READHEAD, i;

  // ---------------- BEGIN --------------

  printf("\n SIMLIB_open(): \n ");

  if ( (fp_simlib_input = fopen(SIMLIB_INPUT.FILE, "rt")) == NULL ) {   
    printf("\n ERROR: cannot open input simlib file: \n");
    printf("\t '%s' \n", SIMLIB_INPUT.FILE );
    MADABORT();
  }


  printf("\t Opened %s \n", SIMLIB_INPUT.FILE );

  // read header keywords. Stop when we reach "BEGIN"

  NLINE_HEADER = 0;
  READHEAD   = 1;

  while( READHEAD > 0 ) {

    fgets(cline, 80, fp_simlib_input) ;
    printf(" Found header line: %s", cline );

    NLINE_HEADER++ ;
    sprintf( HEADER[NLINE_HEADER], "%s", cline );

    // look for header end with "BEGIN" or "#"
    ptrtok = strtok(cline," ") ; // split string

    // remove pre-existing NLIBID key because a new NLIBID key
    // is added at the end.
    if ( strcmp( ptrtok,"NLIBID:") == 0 )
      { NLINE_HEADER--; continue ; }


    if ( strcmp( ptrtok,"BEGIN")   == 0 ) {
      READHEAD = 0;
      char SAVELINE[100];
      sprintf(SAVELINE, "%s", HEADER[NLINE_HEADER] );
      NLINE_HEADER-- ;
      for ( i=1; i <= NLINE_HEADER_ADD; i++ ) {
	NLINE_HEADER++ ;
	sprintf( HEADER[NLINE_HEADER], "%s", HEADER_ADD[i] );
      }
      NLINE_HEADER++ ;
      sprintf( HEADER[NLINE_HEADER], "%s", SAVELINE );
    }

    if ( strcmp( ptrtok,"LIBID:")  == 0 ) READHEAD = 0;

    /* xxxxxxxx mark delete Jun 29 2018 xxxxxxxx
    // look for place to insert extra header stuff
    if ( strcmp( ptrtok,"COMMENT:")  == 0 ) {
      for ( i=1; i <= NLINE_HEADER_ADD; i++ ) {
	NLINE_HEADER++ ;
	sprintf( HEADER[NLINE_HEADER], "%s", HEADER_ADD[i] );
      }
    }
    xxxxxxxxxxx end delete xxxxxxxxxx */


    // look for FILTERS keyword

    sprintf(clast,"XXX");
    while ( ptrtok != NULL  ) {

      if ( strcmp(clast,"FILTERS:") == 0 ) {
	sprintf(SIMLIB_FILTERS, "%s", ptrtok );
	printf(" Found SIMLIB_FILTERS: %s \n",SIMLIB_FILTERS );
      }
      sprintf(clast,"%s", ptrtok);
      ptrtok = strtok(NULL, " ");	
    }
    

  } // end of while 



  // open output file manually instead of using
  // simlib_open() since we just have 80-char lines
  // instead of the detailed info needed by the open function
  // Dump header info to output file.

  if ( (fp_simlib_output = fopen(SIMLIB_OUTPUT.FILE, "wt")) == NULL ) {   
    printf("\n ERROR: cannot open output simlib file: \n");
    printf("\t '%s' \n", SIMLIB_OUTPUT.FILE );
    MADABORT();
  }


  printf("\t Opened %s \n", SIMLIB_OUTPUT.FILE );

  for ( i=1; i <= NLINE_HEADER; i++ ) {
    fprintf(fp_simlib_output, "%s", HEADER[i] );
  }

  // add a few  blank spaces
    fprintf(fp_simlib_output, "\n\n\n" );


}  // end of SIMLIB_open


// ************************************
void SIMLIB_read(int *RDSTAT) {

  /* Copied/modified from snlc_sim.c    
    Read next LIBID in simlib. Fill SIMLIB structure.

  */

  char  c_get[80], c_tmp[80],cfilt[4], STRINGID[20] ;

  int i, NMJD, NMJD_READ, NMJD_ACCEPT, LIBID ;
  int OPTLINE, NOBS, ENDLIB, OKLIBID    ;

  double MJD, CCDGAIN, CCDNOISE, SKYSIG, PSF[3], ZPT[2] ;
  double MAG, RA, DECL, XMW[20], MWEBV    ;

  char  fnam[20] = "SIMLIB_read"  ;

  // ---------------- BEGIN ------------------



  NMJD_READ    =  0  ;
  NMJD_ACCEPT  =  0  ;
  NOBS         =  0  ;
  ENDLIB       =  0  ;
  *RDSTAT      =  2  ;

  SIMLIB_INPUT.LIBID = -9 ;

  while( *RDSTAT != EOF  ) {

    *RDSTAT = fscanf(fp_simlib_input, "%s", c_get) ;

    if ( strcmp(c_get,"END_OF_SIMLIB:") == 0 ) {
    }

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
      readint ( fp_simlib_input, 1, &NOBS );
      SIMLIB_INPUT.NOBS = NOBS ;
    }

    if ( strcmp(c_get,"RA:") == 0 ) 
      readfloat ( fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_RA] );

    if ( strcmp(c_get,"DECL:") == 0 ) 
      readfloat ( fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_DECL] );


    if ( strcmp(c_get,"MWEBV:") == 0  ) 
      readfloat ( fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_MWEBV] );

    if ( strcmp(c_get,"PIXSIZE:") == 0 ) 
      readfloat ( fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_PIXSIZE] );

    if ( strcmp(c_get,"REDSHIFT:") == 0 ) 
      readfloat ( fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_Z] );

    if ( strcmp(c_get,"PEAKMJD:") == 0 ) 
      readfloat ( fp_simlib_input, 1, &SIMLIB_INPUT.INFO_HEAD[IPAR_PKMJD] );

    if ( strcmp(c_get,"FIELD:") == 0 ) 
        fscanf(fp_simlib_input,"%s", SIMLIB_INPUT.FIELDNAME );  


    // Search for MJD-epoch info

    if ( strcmp(c_get,"S:") == 0 ) 
      OPTLINE = 1 ;
    else if ( strcmp(c_get,"T:") == 0 ) 
      OPTLINE = 2 ;
    else
      OPTLINE = 0 ;

    if (OPTLINE > 0 ) {

      NMJD_READ++ ;

      // read line into temporary variables.
      readdouble ( fp_simlib_input, 1, &MJD );
      readchar(fp_simlib_input, STRINGID );
      fscanf(fp_simlib_input, "%s", cfilt );
      readdouble ( fp_simlib_input, 1, &CCDGAIN  );
      readdouble ( fp_simlib_input, 1, &CCDNOISE );
      readdouble ( fp_simlib_input, 1, &SKYSIG );
      readdouble ( fp_simlib_input, 3, PSF );
      readdouble ( fp_simlib_input, 2, ZPT );

      if ( OPTLINE == 1 ) {
	readdouble ( fp_simlib_input, 1, &MAG );
	if ( MAG < 5.0 ) MAG = 99.0 ;
      }
      else
	{ MAG = 99.0 ; } 


      NMJD_ACCEPT++ ;
      SIMLIB_INPUT.NMJD_ACCEPT = NMJD_ACCEPT ;
      NMJD  = NMJD_ACCEPT ;

      if ( NMJD >= MXMJD ) {
	printf("  ERROR: NMJD=%d exceeds array bound at LIBID=%d. \n", 
	       NMJD, LIBID);
	MADABORT();
      }

      SIMLIB_INPUT.MJD[NMJD]      = MJD ;  // June 2017, use double
      sprintf(SIMLIB_INPUT.STRINGID[NMJD], "%s", STRINGID) ;
      sprintf(SIMLIB_INPUT.FILTNAME[NMJD], "%s", cfilt);

      SIMLIB_INPUT.INFO_MJD[NMJD][1]    = CCDGAIN ;
      SIMLIB_INPUT.INFO_MJD[NMJD][2]    = CCDNOISE ;
      SIMLIB_INPUT.INFO_MJD[NMJD][3]    = SKYSIG   ;
      SIMLIB_INPUT.INFO_MJD[NMJD][4]    = PSF[0] ;
      SIMLIB_INPUT.INFO_MJD[NMJD][5]    = PSF[1] ;
      SIMLIB_INPUT.INFO_MJD[NMJD][6]    = PSF[2] ;
      SIMLIB_INPUT.INFO_MJD[NMJD][7]    = ZPT[0] ;
      SIMLIB_INPUT.INFO_MJD[NMJD][8]    = ZPT[1] ;
      SIMLIB_INPUT.INFO_MJD[NMJD][9]    = MAG ;

      if ( MJD < INPUTS.MJD_MIN ||  MJD > INPUTS.MJD_MAX ) {
	NMJD_ACCEPT-- ;
	SIMLIB_INPUT.NMJD_ACCEPT = NMJD_ACCEPT ;
      }


    } // end of OPTLINE if-block


    if ( strcmp(c_get,"END_LIBID:") == 0         ) ENDLIB = 1;
    if ( strcmp(c_get,"#") == 0 && NMJD_READ > 3 ) ENDLIB = 1;

    if ( ENDLIB == 1 ) {

      if ( NMJD_READ != NOBS && INPUTS.OPT_MJD_DMP == 0 ) {
	printf("  ERROR: NMJD_READ = %d, but expected NOBS=%d (LIBID=%d) \n",
	       NMJD, NOBS, LIBID);
	//	MADABORT();
      }

      OKLIBID = -9 ;
      if ( LIBID >= INPUTS.LIBID_MIN && LIBID <= INPUTS.LIBID_MAX ) 
	{ OKLIBID = 1; }

      //      printf("  LIBID=%d  OKLIBID = %d \n", LIBID, OKLIBID );
      if ( OKLIBID < 0 ) { SIMLIB_INPUT.LIBID = -9 ; }


      // if there are too few accepted exposures, set flag to skip LIBID

      if ( NMJD_ACCEPT < INPUTS.MINOBS_ACCEPT ) {
	printf("\t Skipping LIBID %d : only %d exposures. \n", 
	       LIBID, NMJD_ACCEPT);
	SIMLIB_INPUT.LIBID = -9;
      }

      // check option to compute MWEBV from Schlagel maps
      if ( INPUTS.OPT_MWEBV > 0 ) {
	RA   = SIMLIB_INPUT.INFO_HEAD[IPAR_RA] ;
	DECL = SIMLIB_INPUT.INFO_HEAD[IPAR_DECL] ;
	MWgaldust(RA, DECL, &XMW[1], &MWEBV );  // returns XMW & MWEBV
	SIMLIB_INPUT.INFO_HEAD[IPAR_MWEBV] = (float)MWEBV ;
      }

      // skip really larger MW extinctions
      MWEBV = SIMLIB_INPUT.INFO_HEAD[IPAR_MWEBV];
      if ( MWEBV > MWEBV_MAX ) {
	printf("\t Skipping LIBID %d : MWEBV=%6.1f \n", LIBID, MWEBV );
	SIMLIB_INPUT.LIBID = -9;
      }
      return ;
    }

  }  // end of while fscanf loop



} // end of function SIMLIB_read



// *********************
void SIMLIB_coadd(void) {

  //  transfer SIMLIB_input -> SIMLIB_OUTPUT structure,
  //  where the output is in compact form.

  // May 20, 2009: special fix for taking MJD average without roundoff error
  // Jun 20, 2017: MJD -> double instead of float

  int  i, obs, ipar, NOBS_IN, NMEASURE, OVPFILT ;
  int  OBSMIN[MXMJD], OBSMAX[MXMJD]    ;

  char *cfilt, *cfilt_last ;

  double MJD, MJD_LAST, MJD_DIF, XIN, XSUM, XN, XNOPT, ARG, ZPTOFF ;
 
  float *PTR_INFO_OFFSET, *PTR_INFO_INPUT, *PTR_INFO_OUTPUT  ;

  // ------------- BEGIN ----------


  // transfer LIBID & header info without any changes; 

  SIMLIB_OUTPUT.LIBID = SIMLIB_INPUT.LIBID ;

  for ( i=0; i<NPAR_HEAD; i++ ) {
    SIMLIB_OUTPUT.INFO_HEAD[i] = SIMLIB_INPUT.INFO_HEAD[i] ;
  }
  sprintf(SIMLIB_OUTPUT.FIELDNAME, "%s", SIMLIB_INPUT.FIELDNAME ); 

  // ----------------------
  // first loop through and identify MJD-ranges to combine.

  NOBS_IN  = SIMLIB_INPUT.NMJD_ACCEPT ;
  MJD_LAST = -9.0 ;
  NMEASURE = 0 ;

  obs = 1;
  cfilt      = SIMLIB_INPUT.FILTNAME[obs] ;
  cfilt_last = SIMLIB_INPUT.FILTNAME[obs] ;


  for ( obs=1; obs <= NOBS_IN; obs++ ) {
    MJD     = SIMLIB_INPUT.MJD[obs] ;
    MJD_DIF = fabs(MJD - MJD_LAST);
    cfilt   = SIMLIB_INPUT.FILTNAME[obs] ;
    OVPFILT = strcmp(cfilt,cfilt_last);  // 5/14/2009

    if ( MJD_DIF < INPUTS.MAXTDIF_COMBINE && OVPFILT == 0 ) {
      OBSMAX[NMEASURE] = obs ;  //  different MJD, same filter
    } 
    else {
      NMEASURE++ ;                  
      OBSMIN[NMEASURE] = obs ;
      OBSMAX[NMEASURE] = obs ;
    }

    cfilt_last = cfilt ;
    MJD_LAST   = MJD ;

  }  // end of obs loop


  // -----------------
  // Now loop through and combine exposures and take appropriate
  // averages for SIMLIB_OUTPUT structure

  for ( i = 1; i <= NMEASURE; i++ ) {

    SIMLIB_OUTPUT.NOBS  = NMEASURE ;

    // for IDEXPT and filter, copy element from 1st exposure
    // to OUTPUT measurement,

    obs = OBSMIN[i] ;
    sprintf(SIMLIB_OUTPUT.STRINGID[i], "%s", SIMLIB_INPUT.STRINGID[obs]);

    cfilt = SIMLIB_INPUT.FILTNAME[obs] ;
    sprintf(SIMLIB_OUTPUT.FILTNAME[i], "%s", cfilt);

    // setup pointer to output INFO array
    PTR_INFO_OUTPUT = &SIMLIB_OUTPUT.INFO_MJD[i][0] ;

    // init output INFO array  
    SIMLIB_OUTPUT.MJD[i] = 0.0 ;
    for ( ipar=0; ipar < NPAR_MJD; ipar++ ) 
      *(PTR_INFO_OUTPUT + ipar) = 0.0 ;

    // take sums or sums-of-squares over INPUT 'obs' observations

    XN = 0.0 ;
    PTR_INFO_OFFSET  = &SIMLIB_INPUT.INFO_MJD[1][0] ;

    for ( obs = OBSMIN[i]; obs <= OBSMAX[i]; obs++ ) {

      XN += 1.0 ;  // number of exposures for this measurement.

      PTR_INFO_INPUT   = &SIMLIB_INPUT.INFO_MJD[obs][0] ;

      // special case for MJD; subtract offset to avoid float-roundoff
      // when taking average. Add offset back after taking average
      // of residuals.


      XIN = SIMLIB_INPUT.MJD[obs] - SIMLIB_INPUT.MJD[1];
      SIMLIB_OUTPUT.MJD[i] += XIN ;

      XIN = *(PTR_INFO_INPUT + IPAR_CCDGAIN) ;
      *(PTR_INFO_OUTPUT + IPAR_CCDGAIN) += XIN ;

      XIN = *(PTR_INFO_INPUT + IPAR_CCDNOISE) ;
      *(PTR_INFO_OUTPUT + IPAR_CCDNOISE) += (XIN * XIN) ;

      XIN = *(PTR_INFO_INPUT + IPAR_SKYSIG) ;
      *(PTR_INFO_OUTPUT + IPAR_SKYSIG) += (XIN * XIN) ;

      XIN = *(PTR_INFO_INPUT + IPAR_PSF0 + 0) ;
      *(PTR_INFO_OUTPUT + IPAR_PSF0+0) += XIN ;
      XIN = *(PTR_INFO_INPUT + IPAR_PSF0 + 1) ;
      *(PTR_INFO_OUTPUT + IPAR_PSF0+1) += XIN ;
      XIN = *(PTR_INFO_INPUT + IPAR_PSF0 + 2) ;
      *(PTR_INFO_OUTPUT + IPAR_PSF0+2) += XIN ;

      XIN = *(PTR_INFO_INPUT + IPAR_ZPT0 + 0) ;
      ARG = 0.4 * (XIN - 25.0) ;
      *(PTR_INFO_OUTPUT + IPAR_ZPT0+0) += powf(TEN,ARG); 

      XIN = *(PTR_INFO_INPUT + IPAR_ZPT0 + 1) ;
      *(PTR_INFO_OUTPUT + IPAR_ZPT0+1) += XIN ;

      XIN = *(PTR_INFO_INPUT + IPAR_MAG) ;
      *(PTR_INFO_OUTPUT + IPAR_MAG) += XIN ;

    } // end of obs loop

    // now divide by NMEASURE, take sqrt, etc ... to get
    // appropriate result for single measure
    if ( XN == 0.0 ) {
      printf("\n ERROR: Nexposure=0 for Meaure=%d OBS=%d-%d, filt=%s \n",
	     i, OBSMIN[i], OBSMAX[i], cfilt );
      printf(" MJD = %f\n", SIMLIB_INPUT.INFO_MJD[OBSMIN[i]][0] ) ;
      MADABORT();
    }

    if ( INPUTS.OPT_SUM > 0 ) 
      XNOPT = 1 ;
    else
      XNOPT = XN ;

    SIMLIB_OUTPUT.MJD[i] /= XN ;
    SIMLIB_OUTPUT.MJD[i] += SIMLIB_INPUT.MJD[1] ;

    *(PTR_INFO_OUTPUT + IPAR_CCDGAIN) *= (XNOPT/XN) ;

    XSUM = *(PTR_INFO_OUTPUT + IPAR_CCDNOISE) ;
    *(PTR_INFO_OUTPUT + IPAR_CCDNOISE) = sqrtf(XSUM) / XNOPT;

    XSUM = *(PTR_INFO_OUTPUT + IPAR_SKYSIG) ;
    *(PTR_INFO_OUTPUT + IPAR_SKYSIG) = sqrtf(XSUM) / XNOPT;

    *(PTR_INFO_OUTPUT + IPAR_PSF0+0) /= XN ;

    *(PTR_INFO_OUTPUT + IPAR_PSF0+1) /= XN ;

    *(PTR_INFO_OUTPUT + IPAR_PSF0+2) /= XN ;

    XSUM = *(PTR_INFO_OUTPUT + IPAR_ZPT0+0) ;
    *(PTR_INFO_OUTPUT + IPAR_ZPT0+0) = 2.5 * log10f(XSUM/XNOPT) + 25.0 ;

    *(PTR_INFO_OUTPUT + IPAR_ZPT0+1) /= XN ;
    
    if ( INPUTS.OPT_SNLS == 1 ) {
      ZPTOFF = ZPTOFF_SNLS(cfilt);
      *(PTR_INFO_OUTPUT + IPAR_ZPT0+0) += ZPTOFF;

      *(PTR_INFO_OUTPUT + IPAR_SKYSIG) *=1.5 ;
    }


    *(PTR_INFO_OUTPUT + IPAR_MAG) /= XN ;

    if ( INPUTS.OPT_MJD_DMP == 1 ) {
      MJD = SIMLIB_OUTPUT.MJD[i];
      printf(" %f \n", MJD );
    }


  }  // end of i-loop over NMEASURE

  return ;

} // end of SIMLIB_coadd


// *************************
double ZPTOFF_SNLS(char *cfilt) {

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

  printf("\n Insert '%s' after SURVEY line\n", INSERT_LINE);

  sprintf(sedCmd,"sed -i '/SURVEY/a\\ %s' %s", 
	  INSERT_LINE, SIMLIB_OUTPUT.FILE );

  system(sedCmd);
  //  printf("\n xxx sedCmd = \n %s \n", sedCmd);


  return ;

} // end insert_NLIBID
