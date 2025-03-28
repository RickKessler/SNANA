/************************************************
  Created Dec 2004 by R.Kessler

  Collect and compute all calibration and Kcor info and store
  in a single fits file.  

  SALTSHAKER SYNERGY WARNING: 
    The SALT3 training code (SALTshaker; https://github.com/djones1040/SALTShaker)
    reads the output from this code, so make sure to coordinate any changes in the
    output.

  K corrections  are computed on a grid of redshift, rest-frame epoch,
  and AVWARP parameter, where the latter is needed to warp the
  SN spectrum to match observed colors.

  The following are input via user input file:
    - files names for filter responses (lambda transmission)
    - file name for SN Ia template vs. day 
    - file name for primary SED
    - redshfit range for K correction grid

  Lambda binning notes.
  The SN templates with 10 A bins defines the nominal output binning. The filter 
  responses are stored in the original binning, which can be non-uniform, and
  therefore interpolation is used to multiply filter transmission with the same 
  binning as the SN SED.  The primary ref (AB, Vega, BD17..) can have can have 
  different (non-uniform) binning because it is immediately interpolated to have 
  the same binning as the SN. 

  Usage:
    
    kcor.exe  <inFile>   (pass name of user input file)

    kcor.exe  <inFile>  FILTER_LAMSHIFT  g 0.5  r -1.2  i 3.3
         [each number is a shift in A]

    kcor.exe  <inFile>  BD17_SED   <bd17File>
    kcor.exe  <inFile>  VEGA_SED   <vegaFile>
    kcor.exe  <inFile>  PRIMARY_SED  BD17 <bd17File>

    kcor.exe  <inFile> FILTPATH_REPLACE XXXXX Bessell90  
         (replace XXXXX" with "Bessell90")

    kcor.exe  <inFile> MAGSYSTEM_IGNORE BD17
         (ignores BD17 magsystem: allows multiple systems in one file.)

    kcor.exe <inFile> DUPLICATE_LAMSHIFT_GLOBAL 10
       (prep duplicate filter set with lambda shift)

  - - - - - - - - - - - - - - 
  INPUT FILE KEYS:

     SN_SED:  <file with SN flux vs. lam and phase>

     MAGSYSTEM: AB
     MAGSYSTEM: BD17
     MAGSYSTEM: BD17->AB  # read BD17 mag, but transform internally to AB

     # associate each filter set with survey(s) from SURVEY.DEF
     SURYEY: SDSS
     SURVEYS: FOUNDATION,PS1MD

     # optional ZPOFF file to override default ZPOFF.DAT (Jan 2021)
     ZPOFF_FILE: ZPOFF_UPDATED.DAT
       (if no slash in file name, check FILTPATH)

     FILTSYSTEM: COUNT    # most moder systems are count
     FILTSYSTEM: ENERGY   # older Bessell system may be energy

     FILTPATH:  <path>   # location of filter trans files

     FILTER: <stringName>  <fileName>  <mag>
       [note that last char of stringName is used in sim and snlc_fit]

     LAMBDA_RANGE: <lammin> <lammax>
     OUTFILE:  <outFile>              
      [this outFile is input arg for KCOR_FILE in snlc_sim and snlc_fit.]

   To define K-corrections:

       KCOR:  <restFilter>  <obsFilter>   K_ro

       REDSHIFT_RANGE: <zmin> <zmax>  # z-range to store kcor table
       REDSHIFT_BINSIZE: <binSize>    # bin size of kcor table

   AV-warp params for k-corrs:
       RV:          3.1 
       AV_RANGE:   -6.0  6.0  
       AV_BINSIZE:  0.5    # increase for quicker tests
       AV_OPTION:   2      # 2->proper integration over filter



  HISTORY
  ~~~~~~~~~~


 Jan 15 2021: new ZPOFF_FILE input key (for each FILTPATH) to override
              default ZPOFF.DAT

 Nov 3 2021:
   For spectrograph, extend stored wavelength range of SEDs to that
   of spectraograph. See new function set_store_lambda_range().
   Goal is to enable simulated spectra to go beyond wave range of filters.

 Mar 2 2023: abort on invalide filter char; see call to INTFILTER()

 July 2023: write CWD into header (to help find kcor-input file)

 Mar 28 2025; allow command-line override for SPECTROGRAPH

****************************************************/

#include "fitsio.h"
#include "sntools.h"    // defines some general tools
#include "sntools_calib.h"

#include "kcor.h"       // kcor-specific definitions 
#include "MWgaldust.h"
//#include "sntools_cosmology.h"
//#include "genmag_SEDtools.h"
#include "sntools_spectrograph.h"


// =================================================
// =================== MAIN ========================
// =================================================

int main(int argc, char **argv) {

  int i ;

  set_EXIT_ERRCODE(EXIT_ERRCODE_KCOR);

  if (argc < 2) { print_kcor_help();  exit(0); }

  // get name of user input file

  if ( argc >= 2 ) {
    sprintf(INPUTS.inFile_input, "%s", argv[1] );

    NARGV_LIST = argc ;
    for ( i = 0; i < NARGV_LIST ; i++ ) {
      ARGV_LIST[i] = (char*) malloc( MXPATHLEN*sizeof(char) );      
      sprintf( ARGV_LIST[i], "%s", argv[i] );
      USE_ARGV_LIST[i] = 0 ;
    }
    USE_ARGV_LIST[0] = 1 ; // program name
    USE_ARGV_LIST[1] = 1 ; // input file
  }
  else {
    sprintf(INPUTS.inFile_input, "kcor.input" );
  }

  // ------------------------------------------  
  // extract SNDATA_ROOT from env
  sprintf( PATH_SNDATA_ROOT,   "%s",    getenv("SNDATA_ROOT") );
  sprintf( PATH_SNDATA_FILTER, "%s/%s", PATH_SNDATA_ROOT, subdir_filter );
  printf(" SNANA_DIR   = %s \n", getenv("SNANA_DIR") );
  printf(" SNDATA_ROOT = %s \n", PATH_SNDATA_ROOT );
  
  set_FILTERSTRING(FILTERSTRING); // Mar 2023

  // -----------------------------------------

  // read survey name <-> integer map (Nov 2020)
  read_SURVEYDEF();

  // read user input file for directions
  if ( rd_input() != SUCCESS ) { madend(stdout,1) ; }

  // read SN spectra and filter responses
  if ( kcor_ini() != SUCCESS ) { madend(stdout,1) ;  }

  // allocate memory for multi-dimensional arrays.
  if ( malloc_ini() != SUCCESS ) { madend(stdout,1) ; }

  // determine SN color mag vs. day for each filter 
  if ( snmag() != SUCCESS ) { madend(stdout,1) ; }

  //  do K-cor grid vs. redshifts, and days 
  if ( kcor_grid() != SUCCESS ) { madend(stdout,1) ;  }

  // write output 
  if ( kcor_out() != SUCCESS ) { madend(stdout,1) ; }

  // end it all 
  
  happyend();

  return(0);

} // end of main


// ***********************************
void  print_kcor_help(void) {


  static char *help[] = {

    "",
    "\t ***** kcor.exe help menu *****",
    "",
    "# Original kcor program was designed for mlcs2k2 that uses explicit",
    "# k-corrections, but this code has evolved into storing filters ",
    "# and calibration info (without K-corrections) for SALT2/3 models.",
    "# Thus the code name 'kcor' is mis-leading since this is really a",
    "# calibration-storage program.",
    "# The program output is a fits-formatted file (see OUTDIR key below)",
    "# that is the input argument for the KCOR_FILE key in both the ",
    "# snlc_sim.exe and snlc_fit programs.",
    "",
    "# Usage:",
    "   kcor.exe <InputFile>",
    "",
    "# InputFile keys: ",
    "",
    "SN_SED: <file>   # SN flux vs. lam and phase (for k-cors)",
    "",
    "# primary SEDs",
    "BD17_SED:    $SNDATA_ROOT/standards/bd_17d4708_stisnic_003.dat ",
    "VEGA_SED:    $SNDATA_ROOT/standards/alpha_lyr_stis_005.dat  ",
    "AB_SED:      $SNDATA_ROOT/standards/flatnu.dat ",
    "",
    "# Example filter system for DES:",
    "MAGSYSTEM: AB        # AB, BD17, VEGA",
    "FILTSYSTEM: COUNT    # COUNT or ENERGY",
    "FILTPATH:   $SNDATA_ROOT/filters/DES/DES-SN3YR_DECam",
    "SURVEY: DES             # must be defined in $SNDATA_ROOT/SURVEY.DEF",   
    "FILTER:  DES-g   DECam_g.dat   0.0   # stringName  file  ZPoff",
    "FILTER:  DES-r   DECam_r.dat   0.0  ",
    "FILTER:  DES-i   DECam_i.dat   0.0  ",
    "FILTER:  DES-z   DECam_z.dat   0.0  ",    
    "#    (last char of stringName is used in snlc_sim and snlc_fit)",
    "",
    "# other magsystem options:",
    "MAGSYSTEM: BD17      # define BD17 system",
    "MAGSYSTEM: BD17->AB  # read BD17 mag, but transform internally to AB",
    "",
    "# mutiple surveys can match a filter set; e.g., ",
    "SURVEYS: FOUNDATION,PS1MD  ",
    "",
    "# optional ZPOFF file to override default ZPOFF.DAT, ",
    "# or to adjust published mag-system offsets",
    "ZPOFF_FILE: ZPOFF_UPDATED.DAT  ",
    "#   (if no slash in file name, check FILTPATH)",
    "",
    "LAMBDA_RANGE: <lammin> <lammax>   # wave-range for filters and primary",
    "OUTFILE:  <outFile>               # output fits file",
    "#   (this outFile is input arg for KCOR_FILE in snlc_sim and snlc_fit)",
    "",
    "# To define K-corrections: ",
    "KCOR:  <restFilter>  <obsFilter>   K_ro  # ",
    "",
    "REDSHIFT_RANGE: <zmin> <zmax>  # z-range to store kcor table ",
    "REDSHIFT_BINSIZE: <binSize>    # z-bin size of kcor table",
    "",
    "# AV-warp params for k-corrs:"
    "RV:         3.1       " ,
    "AV_RANGE:   -6.0  6.0 " , 
    "AV_BINSIZE:  0.5    # increase for faster kcor generation ",
    "AV_OPTION:   2      # 2->proper integration over filter",
    0
  };

  int i;

  // ----------- BEGIN ------------                                            
  for (i = 0; help[i] != 0; i++)
    { printf ("%s\n", help[i]); }
  
  return;

} // end print_kcor_help

// ******************************
int rd_input(void) {

  // Open and read user input file
  // Apr 25, 2009: initialize NFILT_SYSTEM = 0
  // May 25 2017: read optional LAMSHIFT_GLOBAL
  // Apr 10 2020: check for optional MAGSYST_INPUT->MAGSYS
  // 
  char 
    c_get[80]
    ,c_flt[2][40]
    ,kcorname[40]    // local name for filter 
    ,magcheck_name[40]
    ,MAGSYSTEM_TMP[60]
    ,FILTPATH_LIST[MXFILTDEF+1][MXCHAR_FILENAME]
    ,*ptr_tmp
    ,TXT_MAGREF[60]
    ;   

  INPUT_FILTER_DEF  INPUT_FILTER ;
  MAGSYSTEM_DEF     MAGSYSTEM ;
  FILTSYSTEM_DEF    FILTSYSTEM ;

  FILE *fp_input ;  // name of user input file 

  int iprim, ifilt, itmp, i, gzipFlag, INDX, INDX_INPUT ;
  int NFILT_SYSTEM, LANDOLT_OPT, FILTER_IGNORE   ;

  float xlim4[2];

  double magcheck_value[MXFILTDEF+1]    ;

  char fnam[] = "rd_input"  ;

  /* ----------------- BEGIN ----------------------- */

  // init some stuff to be safe 

   
  NFILTDEF  = NFILTPATH = NFILTDEF_SYN = 0;
  NFILT_SYSTEM  = 0 ;

  INPUTS.NPRIMARY = 0;
  INPUTS.NBIN_REDSHIFT = 1;
  INPUTS.REDSHIFT_MIN  = 0.0 ;
  INPUTS.REDSHIFT_MAX  = 0.0 ;
  INPUTS.REDSHIFT_BINSIZE  = 0.0 ;

  INPUTS.NBIN_AV       = 1;
  INPUTS.AV_MIN        = 0.0 ;
  INPUTS.AV_MAX        = 0.0 ;
  INPUTS.AV_BINSIZE    = 0.0 ;

  INPUTS.AV_OPTION     = 2;

  INPUTS.LEAKAGE_CUT = 0.0 ; // May 20 2021

  INPUTS.IRD_ZPOFF = 1; // read ZPOFF.DAT file by default

  NKCOR         = 0;
  SNSED.NBIN_LAMBDA   = 0;

  NZPOFF = 0;

  INPUTS.DUMP_SNMAG = 0;
  INPUTS.SN_MAGOFF  = 0.0;
  sprintf(INPUTS.SN_TYPE, "NULLTYPE" );

  SNSED.TREST_MIN = -30. ;
  SNSED.TREST_MAX = 100. ;

  // set default to range for ugrizY
  SNSED.LAMBDA_MIN =   2000. ;
  SNSED.LAMBDA_MAX =  12000. ;

  for ( iprim=0; iprim < MXPRIMARY; iprim++ ) {
    PRIMARYSED[iprim].NBIN_LAMBDA = 0;
    PRIMARYSED[iprim].NBIN_LAMBDA_RAW = 0;
    PRIMARYSED[iprim].USE  = 0;
    sprintf(PRIMARYSED[iprim].MAGSYSTEM_NAME,"NULL");
    INPUTS.inFile_PRIMARY[iprim][0] = 0 ;   
  }

  INPUTS.N_OUTFILE = 0 ;
  sprintf(INPUTS.OUTFILE[0],     "NULL" );
  sprintf(INPUTS.OUTFILE[1],     "NULL" );
  FILTSYSTEM.NAME[0]      = 0 ;
  FILTSYSTEM.INDX       = 0;

  INPUTS.PLOTFLAG_SN      = 1 ; // default => plot SNSED and SN mags
  INPUTS.PLOTFLAG_PRIMARY = 1 ; // default => plot primary spec

  sprintf(INPUTS.inFile_snsed,        "NULL" );
  sprintf(INPUTS.inFile_snsed_fudge,  "NULL" );
  INPUTS.snsed_powerlaw = -99.0;

  sprintf(INPUTS.inFile_spectrograph,    "NULL" );
  sprintf(INPUTS.stringOpt_spectrograph, "NULL" );

  sprintf(INPUTS.FILTPATH_replace1,   "NULL" );
  sprintf(INPUTS.FILTPATH_replace2,   "NULL" );
  sprintf(INPUTS.MAGSYSTEM_REPLACE1,     "NULL" );
  sprintf(INPUTS.MAGSYSTEM_REPLACE2,     "NULL" );
  sprintf(INPUTS.MAGSYSTEM_IGNORE,       "NULL" );


  SPECTROGRAPH_USEFLAG = 0 ;
  INPUTS_SPECTRO.SYN_FILTERLIST_BAND[0]=0;

  FILTER_IGNORE         = 0;

  MAGSYSTEM.NAME[0]       = 0 ;
  MAGSYSTEM.NAME_INPUT[0] = 0 ;
  MAGSYSTEM.INDX        = 0 ;
  MAGSYSTEM.INDX_INPUT  = 0 ;
  MAGSYSTEM.DO_TRANSFORM = 0 ;

  INPUTS.OPT_MWCOLORLAW = OPT_MWCOLORLAW_ODON94 ; 
  for(i=0; i<10; i++ ) { INPUTS.PARLIST_MWCOLORLAW[i] = 0.0 ; }
  
  INPUTS.FASTDEBUG    = 0 ;
  INPUTS.SKIPKCOR     = 0 ;
  INPUTS.FLUXERR_FLAG = 0 ;

  VEGAOFF_FLAG_ASTIER06 = 0;

  INPUTS.TREF_EXPLODE = -19.0 ;

  INPUTS.NLAMBIN_FT = 0;

  for ( ifilt=0; ifilt < MXFILTDEF; ifilt++ ) {
    FILTER[ifilt].MASKFRAME   = 0;
    FILTER[ifilt].NBIN_LAMBDA = 0;
    FILTER[ifilt].FILTSYSTEM_INDX = -9 ;
    FILTER[ifilt].IFLAG_DUPLICATE = 0   ;
    FILTER[ifilt].OOB.LAMRANGE[0] = -9.0 ;
    FILTER[ifilt].OOB.LAMRANGE[1] = -9.0 ;
    FILTER[ifilt].OOB.TRANS_MAX_RATIO = 0.0 ;

    ZPOFF_STORE.FILTVALUES[ifilt] = 0.0 ;
    INPUTS.FILTER_LAMSHIFT[ifilt] = 0.0 ; // NOT for duplicates
  }
  INPUTS.NFILTER_LAMSHIFT = 0 ;
  INPUTS.LAMSHIFT_GLOBAL  = 0.0 ; // Angstroms, for duplicates

  // check for FILTPATH overrides BEFORE parsing input file.
  kcor_input_override(1);

  if ( (fp_input = fopen(INPUTS.inFile_input, "rt"))==NULL ) { 
      sprintf ( c1err, "Cannot open input file :" );
      sprintf ( c2err," '%s' ", INPUTS.inFile_input);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      return ERROR;
  }

  printf(" --------------------------------------------------------\n");
  printf("  Read user input file: %s \n", INPUTS.inFile_input );


    /* ---------------------------------------------------    
             start reading input from text file.
       --------------------------------------------------- */

  while( (fscanf(fp_input, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"SN_SED:")==0 )  {
      readchar ( fp_input, INPUTS.inFile_snsed );
    }  
    if ( strcmp(c_get,"SN_SED_POWERLAW:")==0 )  {
      readdouble ( fp_input, 1, &INPUTS.snsed_powerlaw );
    }  

    if ( strcmp(c_get,"SN_SED_FUDGE:")==0 )  {
      readchar ( fp_input, INPUTS.inFile_snsed_fudge );
    }  
    if ( strcmp(c_get,"PLOTFLAG_SN:")==0 )  {
      readint ( fp_input, 1, &INPUTS.PLOTFLAG_SN );
    }  
    if ( strcmp(c_get,"PLOTFLAG_PRIMARY:")==0 )  {
      readint ( fp_input, 1, &INPUTS.PLOTFLAG_PRIMARY );
    }  
    if ( strcmp(c_get,"READ_ZPOFF:")==0 )  {
      readint ( fp_input, 1, &INPUTS.IRD_ZPOFF );
      if ( INPUTS.IRD_ZPOFF == 0 ) { NZPOFF = 0; }
    } 

    if ( strcmp(c_get,"NLAMBIN_FT:")==0 )  {
      readint ( fp_input, 1, &INPUTS.NLAMBIN_FT );
    }  

    if ( strcmp(c_get,"SN_TYPE:")==0 )  {
      readchar ( fp_input, INPUTS.SN_TYPE );
    }  
    if ( strcmp(c_get,"SN_MAGOFF:")==0 )  {
      readdouble ( fp_input, 1, &INPUTS.SN_MAGOFF );
    }  
   
    if ( strcmp(c_get,"LEAKAGE_CUT:")==0 )  {
      readdouble ( fp_input, 1, &INPUTS.LEAKAGE_CUT );
    }  

    if ( strcmp(c_get,"DUPLICATE_LAMSHIFT_GLOBAL:")==0 )  {
      readdouble ( fp_input, 1, &INPUTS.LAMSHIFT_GLOBAL );
    }  

    // ------ primary SED -----
    if ( strcmp(c_get,"VEGA_SED:")==0 )  {
      INPUTS.NPRIMARY++ ;
      iprim = INPUTS.NPRIMARY ;
      sprintf(INPUTS.name_PRIMARY[iprim], "VEGA"  ) ;
      readchar ( fp_input, INPUTS.inFile_PRIMARY[iprim] ) ;
      ENVreplace(INPUTS.inFile_PRIMARY[iprim],fnam,1) ;  
    }  
    if ( strcmp(c_get,"BD17_SED:")==0 )  {
      INPUTS.NPRIMARY++ ;
      iprim = INPUTS.NPRIMARY ;
      sprintf(INPUTS.name_PRIMARY[iprim], "BD17"  ) ;
      readchar ( fp_input, INPUTS.inFile_PRIMARY[iprim] ) ;
      ENVreplace(INPUTS.inFile_PRIMARY[iprim],fnam,1) ;
    }  
    if ( strcmp(c_get,"AB_SED:")==0 )  {
      INPUTS.NPRIMARY++ ;
      iprim = INPUTS.NPRIMARY ;
      sprintf(INPUTS.name_PRIMARY[iprim], "AB"  ) ;
      readchar ( fp_input, INPUTS.inFile_PRIMARY[iprim] ) ;
      ENVreplace(INPUTS.inFile_PRIMARY[iprim],fnam,1) ;
    }  

    if ( strcmp(c_get,"PRIMARY_SED:")==0 )  {
      INPUTS.NPRIMARY++ ;
      iprim = INPUTS.NPRIMARY ;
      readchar ( fp_input, INPUTS.name_PRIMARY[iprim]   ) ;
      readchar ( fp_input, INPUTS.inFile_PRIMARY[iprim] ) ;
      ENVreplace(INPUTS.inFile_PRIMARY[iprim],fnam,1) ;

      if ( INPUTS.NPRIMARY >= MXPRIMARY ) {
	// abort ??
      }
    }

    // -----

    if ( strcmp(c_get,"TREF_EXPLODE:") == 0 )  {
      readdouble ( fp_input, 1, &INPUTS.TREF_EXPLODE );
    }

    if ( strcmp(c_get,"FILTSYSTEM:")==0 )  {
      readchar ( fp_input,  FILTSYSTEM.NAME );

      NFILT_SYSTEM = 0 ;
      LANDOLT_OPT  = 0 ;
      sprintf(magcheck_name, "NULL" );
      for ( i=0; i<=MXFILTDEF; i++ ) 
	magcheck_value[i] = -99.0 ;


      if ( strcmp(FILTSYSTEM.NAME,"COUNT") == 0 ) 
	FILTSYSTEM.INDX = FILTSYSTEM_COUNT;
      else if ( strcmp(FILTSYSTEM.NAME,"ENERGY") == 0 ) 
	FILTSYSTEM.INDX = FILTSYSTEM_ENERGY;
      else {
	sprintf(c1err,"Do not recognize '%s' FILTSYSTEM", FILTSYSTEM.NAME );
	errmsg(SEV_FATAL, 0, fnam, c1err, "Check kcor-input file" ); 
      }
    }


    if ( strcmp(c_get,"MAGSYSTEM:")==0 )  {
      readchar ( fp_input,  MAGSYSTEM_TMP );
      parse_MAGSYSTEM(MAGSYSTEM_TMP, &MAGSYSTEM) ;

      // check for command-line override
      if ( strcmp(MAGSYSTEM.NAME,INPUTS.MAGSYSTEM_REPLACE1) == 0 ) 
	{ sprintf(MAGSYSTEM.NAME,"%s", INPUTS.MAGSYSTEM_REPLACE2); }

      if ( strcmp(MAGSYSTEM.NAME,INPUTS.MAGSYSTEM_IGNORE) == 0 )
	{ FILTER_IGNORE = 1 ; continue ; }
      else
	{ FILTER_IGNORE = 0 ; }

      INDX_INPUT   = index_primary(MAGSYSTEM.NAME_INPUT); 
      INDX         = index_primary(MAGSYSTEM.NAME);
      MAGSYSTEM.INDX_INPUT = INDX_INPUT ;
      MAGSYSTEM.INDX       = INDX ;
      PRIMARYSED[INDX_INPUT].USE = 1; 
      PRIMARYSED[INDX].USE       = 1; 
      MAGSYSTEM.OFFSET_INPUT     = 2.5 * log10 ( FNU_AB );
      MAGSYSTEM.OFFSET           = 2.5 * log10 ( FNU_AB );

      PRIMARYSED[INDX].IS_AB       = (strcmp(MAGSYSTEM.NAME,"AB")==0) ;
      PRIMARYSED[INDX_INPUT].IS_AB = (strcmp(MAGSYSTEM.NAME_INPUT,"AB")==0) ;

      printf("\n\t Found MAGSYSTEM '%s' with offset = %8.3f (INDX=%d->%d)\n",
	     MAGSYSTEM_TMP, MAGSYSTEM.OFFSET,
	     MAGSYSTEM.INDX_INPUT, MAGSYSTEM.INDX );
      sprintf(MAGSYSTEM.SURVEY_NAMES, "NONE" );
      sprintf(MAGSYSTEM.ZPOFF_FILE,   "NONE" );
    }  

    if ( strcmp(c_get,"SURVEY:") == 0 || strcmp(c_get,"SURVEYS:") == 0 ) {
      readchar(fp_input, MAGSYSTEM.SURVEY_NAMES);  // Nov 2020
      check_valid_survey_names(MAGSYSTEM.SURVEY_NAMES);
    }

    if ( strcmp(c_get,"ZPOFF_FILE:") == 0 ) {
      readchar(fp_input, MAGSYSTEM.ZPOFF_FILE);  // Jan 2021
    }

    if ( strcmp(c_get,"FILTPATH:")==0 && (FILTER_IGNORE == 0) )  {

      NFILTPATH++ ;
      readchar ( fp_input, INPUTS.FILTPATH );
      sprintf(INPUTS.FILTPATH_ORIG, "%s", INPUTS.FILTPATH);
      ENVreplace(INPUTS.FILTPATH,fnam,1);

      if ( strcmp(INPUTS.FILTPATH,INPUTS.FILTPATH_replace1) == 0 ) 
	{  sprintf(INPUTS.FILTPATH,"%s", INPUTS.FILTPATH_replace2 ); }

      rd_ZPOFF(INPUTS.FILTPATH, MAGSYSTEM.ZPOFF_FILE);  

      // store filter path to check for duplicates
      sprintf(FILTPATH_LIST[NFILTPATH], "%s", INPUTS.FILTPATH);

    }

    // - - - - - - - - - - - - - - -
   
    // July 2016: check for spectrograph
    if ( strcmp(c_get,"SPECTROGRAPH:") == 0 ) {      
      readchar   ( fp_input, INPUTS.inFile_spectrograph );

      // extract option in (). E.g., if inFile = myspec.data(blabla),
      // then inFile->myspec.dat and stringOpt->blabla.
      extractStringOpt(INPUTS.inFile_spectrograph,      // (I,O)
		       INPUTS.stringOpt_spectrograph);  // (O)
      ENVreplace(INPUTS.inFile_spectrograph,fnam,1);

      SPECTROGRAPH_USEFLAG = 1; 
    }

    // -------------------- 
    if ( strcmp(c_get,"FILTER:")==0 &&  (FILTER_IGNORE == 0) )  {

      readchar   ( fp_input, INPUT_FILTER.FILTNAME );
      readchar   ( fp_input, INPUT_FILTER.FILENAME );
      readchar   ( fp_input, TXT_MAGREF);
      parse_MAGREF(INPUT_FILTER.FILTNAME, TXT_MAGREF, 
		   &INPUT_FILTER.MAGREF ); // <== returned

      INPUT_FILTER.IFLAG_SYN = 0;     
      NFILTDEF++ ;      // increment total no. filters
      NFILT_SYSTEM++ ;  // increment NFILT for current system
      storeFilterInfo( &INPUT_FILTER, &MAGSYSTEM, &FILTSYSTEM );   
    }  // end of FITLER: read


    // --------------------
    // check for synthetic filter (after SPECTROGRAPH key)
    if ( strcmp(c_get,"SYN_FILTER:")==0  ) {

      readchar   ( fp_input, INPUT_FILTER.FILTNAME );
      readdouble ( fp_input, 2, INPUT_FILTER.LAMRANGE_SYN );
      readchar   ( fp_input, TXT_MAGREF);
      parse_MAGREF(INPUT_FILTER.FILTNAME, TXT_MAGREF, 
		   &INPUT_FILTER.MAGREF ); // <== returned

      if ( SPECTROGRAPH_USEFLAG == 0 ) {
	sprintf(c1err,"Cannot define SYN_FILTER %s (%.0f - %.0f A)",
		INPUT_FILTER.FILTNAME, 
		INPUT_FILTER.LAMRANGE_SYN[0],
		INPUT_FILTER.LAMRANGE_SYN[1]);
	sprintf(c2err,"because SPECTROGRAPH key is not defined.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      INPUT_FILTER.IFLAG_SYN = 1;
      sprintf(INPUT_FILTER.FILENAME,"NULL");
      
      NFILTDEF_SYN++ ;
      NFILTDEF++ ;      // increment total no. filters
      NFILT_SYSTEM++ ;  // increment NFILT for current system
      
      storeFilterInfo( &INPUT_FILTER, &MAGSYSTEM, &FILTSYSTEM );
    
    }  // end of IFU: read
       
    
    // - - - - - - - - - - - - - - - 
    
    // look for option to compare SN mags at t=0
    if ( strcmp(c_get,"SNMAG0_CHECK:") == 0 )  {
      readdouble ( fp_input, NFILT_SYSTEM, magcheck_value );
      readchar   ( fp_input, magcheck_name );
    }

    // store items for filters from most recent system

    i=0;
    for ( ifilt=NFILTDEF-NFILT_SYSTEM+1; ifilt<=NFILTDEF; ifilt++ ) {
      FILTER[ifilt].SNMAG0_CHECK_VALUE   = magcheck_value[i] ;
      sprintf(FILTER[ifilt].SNMAG0_CHECK_NAME, "%s", magcheck_name );
      i++ ;
    }

    if ( strcmp(c_get,"REDSHIFT_RANGE:")==0 )  {
      readfloat ( fp_input, 2, xlim4 );
      INPUTS.REDSHIFT_MIN = xlim4[0];
      INPUTS.REDSHIFT_MAX = xlim4[1];
    }

    if ( strcmp(c_get,"TREST_RANGE:")==0 )  {
      readfloat ( fp_input, 2, xlim4 );
      SNSED.TREST_MIN = xlim4[0];
      SNSED.TREST_MAX = xlim4[1];
    }

    if ( strcmp(c_get,"REDSHIFT_BINSIZE:")==0 )  {
      readdouble ( fp_input, 1, &INPUTS.REDSHIFT_BINSIZE );   
    }  

    if ( strcmp(c_get,"AV_RANGE:")==0 )  {
      readfloat ( fp_input, 2, xlim4 );
      INPUTS.AV_MIN = xlim4[0] ;
      INPUTS.AV_MAX = xlim4[1] ;
    }

    if ( strcmp(c_get,"AV_BINSIZE:")==0 )  
      { readdouble ( fp_input, 1, &INPUTS.AV_BINSIZE );   }  

    if ( strcmp(c_get,"OPT_MWCOLORLAW:")==0 )  
      { readint(fp_input, 1, &INPUTS.OPT_MWCOLORLAW );    }

    if ( strcmp(c_get,"RV:")==0 )  
      { readdouble ( fp_input, 1, &INPUTS.RV_MWCOLORLAW );  }  
    if ( strcmp(c_get,"RV_MWCOLORLAW:")==0 )  
      { readdouble ( fp_input, 1, &INPUTS.RV_MWCOLORLAW );  }  

    if ( strcmp(c_get,"AV_OPTION:")==0 )  
      { readint ( fp_input, 1, &INPUTS.AV_OPTION );  }  


    if ( strcmp(c_get,"LAMBDA_RANGE:")==0 )  {
      readfloat ( fp_input, 2, xlim4 );
      SNSED.LAMBDA_MIN = xlim4[0];
      SNSED.LAMBDA_MAX = xlim4[1]; 
    }  


    // check for adding out-of-band transmission
    // OOB: [bandList] MINLAM MAXLAM Trans/TranMax
    if ( strcmp(c_get,"FILTER_OOB:") == 0 ) {
      char bandList[40];  double LAMRANGE[2], RATIO ;
      readchar(fp_input, bandList);  
      readdouble(fp_input, 2, LAMRANGE);
      readdouble(fp_input, 1, &RATIO);
      parse_OOB(bandList, LAMRANGE, RATIO);      
    } 

    if ( strcmp(c_get,"KCOR:")==0  &&  INPUTS.SKIPKCOR == 0 )  {
       readchar ( fp_input, c_flt[0] );
       readchar ( fp_input, c_flt[1] );
       readchar ( fp_input, kcorname );

       NKCOR++;

       if ( NKCOR > MXKCOR ) {
	  sprintf(c1err, "%d K-cor matrices exceeds array bound.", 
		  NKCOR );
	  errmsg(SEV_FATAL, 0, fnam, c1err, " "); 
	  return ERROR;                  
       }

       sprintf ( KCORLIST[NKCOR][0], "%s", c_flt[0] );
       sprintf ( KCORLIST[NKCOR][1], "%s", c_flt[1] );
       sprintf ( KCORSYM[NKCOR],     "%s", kcorname );

       printf("\t  Do KCOR for rest '%s' and redshifted '%s' (%s) \n", 
	      KCORLIST[NKCOR][0], KCORLIST[NKCOR][1], KCORSYM[NKCOR] );

    }


    // -----------

    if ( strcmp(c_get,"OUTFILE:")==0 ) {
      i = INPUTS.N_OUTFILE ;
      readchar ( fp_input, INPUTS.OUTFILE[i] ); 
      INPUTS.N_OUTFILE++ ;
    }
    // -----------

    if ( strcmp(c_get,"DUMP_SNMAG:")==0 ) 
      { readint ( fp_input, 1, &INPUTS.DUMP_SNMAG ); }


  }  // end of fscanf while


  // ==================================================
  // check for command-line overrides
  // ==================================================

  kcor_input_override(2);
  check_argv();

  // -----------------------------------------------------
  // check option to add duplicate filters with LAMSHIFT
  ADDFILTERS_LAMSHIFT_GLOBAL();

  // Get number of Z and AV bins

  get_NZBIN();
  get_NAVBIN();


  // check option to read spectrograph file
  if ( SPECTROGRAPH_USEFLAG ) {
    char spectro_fileName[MXPATHLEN], SPECTRO_FILENAME[MXPATHLEN];
    char *inFile = INPUTS.inFile_spectrograph ;
    FILE *fp ;

    int OPENMASK = OPENMASK_VERBOSE;
    fp = snana_openTextFile(OPENMASK, INPUTS.FILTPATH, inFile,
			    SPECTRO_FILENAME, &gzipFlag );

    if ( !fp ) {
      sprintf(c1err, "Could not open SPECTROGRAPH file");
      sprintf(c2err, "%s", spectro_fileName );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    fclose(fp);
    init_spectrograph(SPECTRO_FILENAME,INPUTS.stringOpt_spectrograph);
    
  }

  // ==================================================
  //   parrot the input after command-line overrides
  // ==================================================

  printf("\t NPRIMARY:          %d \n", INPUTS.NPRIMARY);
  printf("\t SN Ia templates:   %s \n", INPUTS.inFile_snsed );
  printf("\t SN SED Fudge file: %s \n", INPUTS.inFile_snsed_fudge );
  printf("\t TREF_EXPLODE = %6.1f days (rest-frame)\n", 
	 INPUTS.TREF_EXPLODE);

  printf("\t REDSHIFT range from %5.3f to %5.3f \n", 
	 INPUTS.REDSHIFT_MIN, INPUTS.REDSHIFT_MAX ) ;
  printf("\t REDSHIFT binsize : %5.3f (NZBIN=%d) \n", 
	 INPUTS.REDSHIFT_BINSIZE, INPUTS.NBIN_REDSHIFT );
  printf("\t AV range from %4.2f to %4.2f \n", 
	 INPUTS.AV_MIN, INPUTS.AV_MAX ) ;
  printf("\t AV binsize : %4.2f  (NAVBIN=%d) \n", 
	 INPUTS.AV_BINSIZE,  INPUTS.NBIN_AV );

  printf("\t RV = A(V)/E(B-V) : %4.2f  \n", INPUTS.RV_MWCOLORLAW );
    
  char *ptrStr = INPUTS.STR_MWCOLORLAW ;
  text_MWoption("COLORLAW", INPUTS.OPT_MWCOLORLAW, ptrStr, fnam );
  printf("\t Galactic extinction law: %s \n", ptrStr );


  printf("\t AV option : %d  => ", INPUTS.AV_OPTION );

  if ( INPUTS.AV_OPTION == 1 ) 
    printf(" multiply Kcor x A(<lam>)  \n");
  else if ( INPUTS.AV_OPTION == 2 ) 
    printf(" convolve Kcor * A(lam)  \n");
  else {
    printf("\n");
    sprintf(c1err, "Illegal AV_OPTION = %d", INPUTS.AV_OPTION );
    errmsg(SEV_FATAL, 0, fnam, c1err, "        "); 
  }

  printf("\t MW Extinction slope measured at E(B-V)=%5.2f, RVMW=%.3f \n",
	 MWEBV_LIST[1], INPUTS.RV_MWCOLORLAW );

  printf("\t User SN LAMBDA range from %6.0f to %6.0f Angstroms \n",
	 SNSED.LAMBDA_MIN, SNSED.LAMBDA_MAX );

  printf("\t User SN Trest range from %6.0f to %6.0f days \n",
	 SNSED.TREST_MIN, SNSED.TREST_MAX );

  if ( INPUTS.LEAKAGE_CUT > 1.0E-12 ) {
    printf("\t Reject out-of-band (OOB) leakage with "
	   "Trans/TransMax < %9.2le\n", INPUTS.LEAKAGE_CUT );
  }

  // ==================================================
  // ==================================================
  // compute a few extra things from the input


  GET_SNSED_FUDGE((double)0.0);

  if ( INPUTS.SKIPKCOR == 1 ) { NKCOR = 0; }

  fclose ( fp_input );

  printf("\n  Done reading user input file.\n\n" );

  fflush(stdout);

  return SUCCESS;
}  // end of rd_input


// *********************************************************
void  check_valid_survey_names(char *SURVEYS) {

  // Created Nov 23 2020
  // Check each survey in comma-sep list of *SURVEYS;
  // abort if any survey name is invalid.

  int  n_survey, i, ID ;
  char **survey_list, *survey ;
  char fnam[] = "check_valid_survey_names" ;

  // ------------- BEGIN ---------------

  parse_commaSepList("SURVEYS", SURVEYS, 10, MXCHAR_FILENAME,
		     &n_survey, &survey_list );

  for(i=0; i < n_survey; i++ ) {
    survey = survey_list[i];
    ID = get_IDSURVEY(survey);
    if ( ID < 0 ) {
      sprintf(c1err, "Invalid survey = '%s'  "
	      "(check SURVEY keys in kcor-input).", survey);
      sprintf(c2err, "SURVEY name must exist in $SNDATA_ROOT/SURVEY.DEF");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  return ;

} // end check_valid_survey_names

// *********************************************************
void parse_OOB(char *bandList, double *LAMRANGE, double RATIO) {

  // Parse OOB info and store it in structure.

  int NFTMP, i, ifilt ;
  char band[2];
  //  char fnam[] = "parse_OOB" ;

  // -------------

  NFTMP = strlen(bandList);
  for(i=0; i < NFTMP; i++ ) {
    sprintf(band,"%c", bandList[i]);
    ifilt = INTFILTER_kcor(band);
    //	printf(" xxx load OOB(%s): LAMRANGE=%.1f to %.1f, RATIO=%le\n",
    //     band, LAMRANGE[0], LAMRANGE[1], RATIO);
    FILTER[ifilt].OOB.LAMRANGE[0] = LAMRANGE[0] ;
    FILTER[ifilt].OOB.LAMRANGE[1] = LAMRANGE[1] ;
    FILTER[ifilt].OOB.TRANS_MAX_RATIO = RATIO ;
  }
  
  return ;

} // end parse_OOB

// ***********************
void  parse_MAGREF(char *FILTNAME, char *TXT_MAGREF, double *MAGREF ) {

  // Created Jun 1 2017
  // Parse input TXT_MAGREF and return MAGREF.
  // FILTNAME is used only for error message.
  //
  // Examples:
  // TXT_MAGREF = 9.824           ---> MAGREF = 9.824
  // TXT_MAGREF = 9.824+.011      ---> MAGREF = 9.835
  // TXT_MAGREF = 9.824-.011      ---> MAGREF = 9.813
  // TXT_MAGREF = 9.824+.011+.02  ---> MAGREF = 9.855
  //
  // Motivation: allows easier adjusting ref mags in kcor-input file.
  //
  // Mar 2 2023: abort if FILTNAME is not valid

  int LENTXT  = strlen(TXT_MAGREF);
  int LENFILT = strlen(FILTNAME);
  int i, ISPLUS, ISMINUS, ISIGN, LASTCHAR, IFILTDEF ;
  double MAGREF_LOCAL, TMPNUM, XSIGN ;
  char   c1[2], cnum[40] ;
  char fnam[] = "parse_MAGREF" ;

  // ------------- BEGIN ---------------

  MAGREF_LOCAL = 0.0 ;
  XSIGN = +1.0 ;  cnum[0] = 0;

  //  Mar 2023: check for valid filter char at end of FILTNAME
  IFILTDEF = INTFILTER(FILTNAME);
  if ( IFILTDEF <= 0 ) {
    sprintf(c1err, "last filter char of '%s' is invalid", FILTNAME);
    sprintf(c2err, "Check valid chars with "
	    "'grep FILTERSTRING_DEFAULT sntools.h'");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  for(i=0; i < LENTXT; i++ ) {
    sprintf(c1,"%c", TXT_MAGREF[i] );
    ISPLUS  = ( i>0 && c1[0] == '+' ) ;
    ISMINUS = ( i>0 && c1[0] == '-' ) ;
    LASTCHAR = ( i == LENTXT-1 );

    // if NOT + or - then increment number string cnum
    if ( ISPLUS==0 || ISMINUS==0 ) { strcat(cnum,c1); }

    // when reacing + or - or last character,
    // increment MAGREF with last cnum value using last XSIGN
    if ( ISPLUS || ISMINUS || LASTCHAR ) { 
      sscanf(cnum, "%le", &TMPNUM);       
      MAGREF_LOCAL += (XSIGN*TMPNUM);
      
      // set XSIGN for next cnum
      ISIGN = ISPLUS - ISMINUS;  XSIGN = (double)ISIGN;
      cnum[0] = 0;  // reset cnum string
    }
      
  }

  *MAGREF = MAGREF_LOCAL ;

  return ;

} // end parse_MAGREF


// **************************************
void parse_MAGSYSTEM(char *MAGSYSTEM_ARG, MAGSYSTEM_DEF *MAGSYSTEM) {

  // Created Aril 2020 by R.Kessler
  // 
  // Examples:
  // Input MAGSYSTEM_ARG = VEGA    
  // Output
  //    MAGSYSTEM->NAME_INPUT = VEGA
  //    MAGSYSTEM->NAME       = VEGA
  //    MAGSYSTEM->DO_TRANSFORM = 0 (false)
  //
  // Input MAGSYSTEM_ARG = VEGA->AB 
  // Output
  //    MAGSYSTEM->NAME_INPUT = VEGA
  //    MAGSYSTEM->NAME       = AB
  //    MAGSYSTEM->DO_TRANSFORM = 1 (true)

  int i, len, jdash;
  char TMP_ARG[60];
  //  char fnam[] = "parse_MAGSYSTEM" ;

  // ----------- BEGIN -----------

  MAGSYSTEM->NAME_INPUT[0] = MAGSYSTEM->NAME[0] = 0 ;
  MAGSYSTEM->DO_TRANSFORM  = 0 ;

  if ( strstr(MAGSYSTEM_ARG,"->") == NULL ) {
    sprintf(MAGSYSTEM->NAME_INPUT, "%s", MAGSYSTEM_ARG);
    sprintf(MAGSYSTEM->NAME,       "%s", MAGSYSTEM_ARG);
  }
  else {
    MAGSYSTEM->DO_TRANSFORM  = 1 ;
    len = strlen(MAGSYSTEM_ARG);
    jdash = -9 ;
    for(i=0; i < len; i++ ) 
      { if ( MAGSYSTEM_ARG[i] == '-' ) {jdash=i;} }
    
    sprintf(TMP_ARG,"%s", MAGSYSTEM_ARG);

    snprintf(MAGSYSTEM->NAME_INPUT, jdash+1, "%s", TMP_ARG);
    // xxx merde !!  strncpy(MAGSYSTEM->NAME_INPUT,  TMP_ARG, jdash);
    sprintf(MAGSYSTEM->NAME, "%s", &MAGSYSTEM_ARG[jdash+2] );

    /*
    printf(" xxx %s: jdash=%d  NAME = '%s' -> '%s' \n", 
	   fnam, jdash, MAGSYSTEM->NAME_INPUT, MAGSYSTEM->NAME );
    */
    
  }

  return ;

} // end parse_MAGSYSTEM

// ******************************
void  storeFilterInfo(INPUT_FILTER_DEF *INPUT_FILTER,
		      MAGSYSTEM_DEF    *MAGSYSTEM,
		      FILTSYSTEM_DEF   *FILTSYSTEM) {

  // Created Dec 28 2015 by R.Kessler
  // Code moved from rd_input to here as part of refactor to
  // to prepare for reading IFU info.
  // Store filter info in global arrays.
  //

  int  INDX        = MAGSYSTEM->INDX;
  int  INDX_INPUT  = MAGSYSTEM->INDX_INPUT ;
  char *NAME       = MAGSYSTEM->NAME;
  char *NAME_INPUT = MAGSYSTEM->NAME_INPUT ;
  char *SURVEY     = MAGSYSTEM->SURVEY_NAMES ; // Nov 2020
  char *ZPOFF_FILE = MAGSYSTEM->ZPOFF_FILE;    // Jan 2021
  double OFFSET       = MAGSYSTEM->OFFSET;
  double OFFSET_INPUT = MAGSYSTEM->OFFSET_INPUT;

  int NF, lenf, IFLAG_SYN  ;
  char FILENAME[MXPATHLEN], band[4];
  char fnam[] = "storeFilterInfo" ;

  // strip inputs into local variables
  char *filtName = INPUT_FILTER->FILTNAME ;
  char *fileName = INPUT_FILTER->FILENAME ;  
  double magRef  = INPUT_FILTER->MAGREF ;
  
  // ------------- BEGIN -----------------

  // check if this is a synthetic filter from the spectrograph
  IFLAG_SYN = INPUT_FILTER->IFLAG_SYN ;

  // load file name
  sprintf(FILENAME, "%s/%s", INPUTS.FILTPATH, fileName);
  
  NF = NFILTDEF ;
  
  if ( NF >= MXFILTDEF ) {			      
    sprintf(c1err, "%d filters exceeds array bound of %d.", 
	    NF, MXFILTDEF );
    errmsg(SEV_FATAL, 0, fnam, c1err, " "); 
  }

  if ( fabs(magRef) > 22.0 ) {
    sprintf(c1err, "Crazy magRef = %f for filter = %s ", magRef, filtName );
    errmsg(SEV_FATAL, 0, fnam, c1err, " "); 
  }

  // get one-band char representation.
  lenf = strlen(filtName);
  sprintf(band, "%c", filtName[lenf-1] ); // last char is band
  sprintf ( FILTER[NF].band, "%s" , band );
  sprintf ( FILTER[NF].name, "%s" , filtName );
  sprintf ( FILTER[NF].file, "%s" , FILENAME ) ;  


  if ( INDX <= 0 || INDX > 10 ) {
    sprintf(c1err, "MAGSYSTEM = '%s' is not defined for FILTER='%s'",
	    NAME, filtName );
    sprintf(c2err,"Check input file: %s ", INPUTS.inFile_input);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  if ( FILTSYSTEM->INDX <= 0 || FILTSYSTEM->INDX > 10 ) {
    sprintf(c1err, "FILTSYSTEM = '%s' is not defined for FILTER='%s'",
	    FILTSYSTEM->NAME, filtName );
    sprintf(c2err,"Check input file: %s ", INPUTS.inFile_input);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // assign current MAGSYSTEM parameters
  
  FILTER[NF].MAGSYSTEM_OFFSET       = OFFSET ;
  FILTER[NF].MAGSYSTEM_INDX         = INDX ;
  FILTER[NF].MAGSYSTEM_INDX_INPUT   = INDX_INPUT;

  /*
  printf(" xxx %s: NF=%d load INDX = %d -> %d \n",
	 fnam, NF, INDX_INPUT, INDX );
  */

  sprintf(FILTER[NF].MAGSYSTEM_NAME,"%s", NAME ) ;
  FILTER[NF].MAGFILTER_REF    = magRef ;
  FILTER[NF].FILTSYSTEM_INDX  = FILTSYSTEM->INDX ;

  sprintf(FILTER[NF].FILTSYSTEM_NAME,"%s", FILTSYSTEM->NAME ) ;
  sprintf(FILTER[NF].PATH,     "%s", INPUTS.FILTPATH ) ;
  sprintf(FILTER[NFILTPATH].PATH_ORIG,"%s", INPUTS.FILTPATH_ORIG ) ;  
  FILTER[NF].IPATH = NFILTPATH ;  
  
  FILTER[NF].MAGFILTER_ZPOFF  = get_ZPOFF(filtName,NFILTPATH) ;

  FILTER[NF].SNMAG0_CHECK_VALUE  = -99.0 ; // fill with NULL
  
  // 4/25/2009 - load PRIMARYSED struct as well.
  PRIMARYSED[INDX].MAGSYSTEM_OFFSET = OFFSET ;
  sprintf(PRIMARYSED[INDX].MAGSYSTEM_NAME,"%s", NAME ) ;  
  sprintf(PRIMARYSED[INDX].MAGSYSTEM_SEDFILE,"%s", 
	  INPUTS.inFile_PRIMARY[INDX] ) ;

  PRIMARYSED[INDX_INPUT].MAGSYSTEM_OFFSET = OFFSET_INPUT ;
  sprintf(PRIMARYSED[INDX_INPUT].MAGSYSTEM_NAME,"%s", NAME_INPUT ) ;
  sprintf(PRIMARYSED[INDX_INPUT].MAGSYSTEM_SEDFILE,"%s", 
	  INPUTS.inFile_PRIMARY[INDX_INPUT] ) ;

  if ( REQUIRE_SURVEY_KCOR && IGNOREFILE(SURVEY) ) {
    sprintf(c1err,"Must associate SURVEY with filter %s", filtName);
    sprintf(c2err,"Add  'SURVEY: <surveyList>' in kcor-input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(FILTER[NF].SURVEY_NAMES, "%s", SURVEY);     // Nov 2020
  sprintf(FILTER[NF].ZPOFF_FILE,   "%s", ZPOFF_FILE); // Jan 2021

  if ( IFLAG_SYN ) {
    double L0 = INPUT_FILTER->LAMRANGE_SYN[0] ;
    double L1 = INPUT_FILTER->LAMRANGE_SYN[1] ;
    FILTER[NF].IFLAG_SYN  = IFLAG_SYN ;
    FILTER[NF].LAMBDA_MIN = L0 ;
    FILTER[NF].LAMBDA_MAX = L1 ;
    printf("    (%2d) %s  %s-flux => SYN_FILTER(%8.0f-%8.0f) \n" ,
	   NF, filtName, FILTSYSTEM->NAME, L0, L1 ) ;

    //update list of filter-bands defined by spectrograph
    char *ptrBand = INPUTS_SPECTRO.SYN_FILTERLIST_BAND ;
    sprintf(ptrBand,"%s%s", ptrBand, band);
  } 
  else {
    FILTER[NF].IFLAG_SYN = 0 ;
    printf("    (%2d) %s  %s-flux => %s \n" ,
	   NF, filtName, FILTSYSTEM->NAME, FILENAME );
  }

  return ;
  
} // end storeFilterInfo

// ******************
void get_NZBIN(void) {

  double tmpdif, xbin;
  char fnam[] = "get_NZBIN" ;

  // ----------- BEIGN ----------

  if ( INPUTS.REDSHIFT_BINSIZE > 0.0 ) { 
    tmpdif = INPUTS.REDSHIFT_MAX - INPUTS.REDSHIFT_MIN + 0.000001 ;
    xbin   = tmpdif / INPUTS.REDSHIFT_BINSIZE + 1.0;
    INPUTS.NBIN_REDSHIFT = (int)xbin;
  }
  else 
    { xbin = 1.0 ; INPUTS.NBIN_REDSHIFT = 1; }
  

  tmpdif = (double)INPUTS.NBIN_REDSHIFT - xbin ;


  if ( fabs(tmpdif) > 0.001 ) {
    sprintf(c1err, "redshift binsize (%5.3f) and range (%5.3f - %5.3f)"
	    ,INPUTS.REDSHIFT_BINSIZE
	    ,INPUTS.REDSHIFT_MIN
	    ,INPUTS.REDSHIFT_MAX );
    sprintf(c2err, "are not compatible. Check input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( INPUTS.NBIN_REDSHIFT >= MXREDSHIFT ) {
    sprintf(c1err, "%d redshift bins exceeds array bound.", 
	    INPUTS.NBIN_REDSHIFT );
    errmsg(SEV_FATAL, 0,fnam, c1err, " "); 
  }

} // end of get_NZBIN

// ******************
void get_NAVBIN(void) {

  double xbin, xn;
  char fnam[] = "get_NAVBIN" ;

  // ------------- BEGIN -----------

  if ( INPUTS.AV_BINSIZE  > 0.0 ) {
    xbin    = (INPUTS.AV_MAX - INPUTS.AV_MIN) / INPUTS.AV_BINSIZE + 1.0001 ;
    INPUTS.NBIN_AV = (int)xbin ;
  } else
    { xbin = 1.0 ; INPUTS.NBIN_AV = 1; }


  xn = (double)INPUTS.NBIN_AV ;
  if ( xn < xbin-0.001 ) {
    sprintf(c1err, "av binsize (%5.3f) and range (%5.3f - %5.3f)",
	    INPUTS.AV_BINSIZE, INPUTS.AV_MIN, INPUTS.AV_MAX );
    sprintf(c2err, "are not compatible. NBIN=%f  xbin=%f. Check input file.", 
	    xn, xbin );

    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( INPUTS.NBIN_AV >= MXAV ) {
    sprintf(c1err, "%d AV bins exceeds array bound.", 
	    INPUTS.NBIN_AV );
    errmsg(SEV_FATAL, 0, fnam, c1err, " "); 
  }

} // end of get_NAVBIN



// *********************************
void kcor_input_override(int OPT) {

  // Check for command-line overrides.
  // argv[0] is the name of the program.
  // argv[1] is the name of the input file.
  // start parsing at argv[2].
  //
  // OPT=1 --> parse before reading text file;
  // OPT=2 --> parse after reading text file
  //
  //
  // Feb 2019: read FILTER_OOB
  // Jan 15 2021: check ZPOFF_FILE
  // Mar 28 2025: check SPECTROGRAPH

  int i, ilast, iuse ;
  char tmpName[60], tmpFile[MXPATHLEN] ;
  char fnam[] = "kcor_input_override" ;

  // ------------ BEGIN -----------

  i = 2; ilast = 2 ;

  while ( i < NARGV_LIST ) {

    // start parsing command-line args that need to be known
    // before reading text file.

    // for FILTPATH, 1st arg is what is in file; 2nd arg is override
    if ( strcmp( ARGV_LIST[i], "FILTPATH_REPLACE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.FILTPATH_replace1 );
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.FILTPATH_replace2 );
    }

    if ( strcmp( ARGV_LIST[i], "MAGSYSTEM_IGNORE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.MAGSYSTEM_IGNORE );
    }

    if ( OPT == 1 ) { goto SKIPPY ; }

    // after reading  text file, parse all command-line args.
    printf("  PROCESS COMMAND LINE ARG: %s \n", ARGV_LIST[i] );

    // -----------------
    if ( strcmp( ARGV_LIST[i],"OUTFILE") == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.OUTFILE[0] ); 
    }
    // ----------


    if ( strcmp( ARGV_LIST[i], "DUMP_SNMAG" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.DUMP_SNMAG ); 
    }

    if ( strcmp( ARGV_LIST[i], "SKIPKCOR" ) == 0 ) 
      { INPUTS.SKIPKCOR = 1; USE_ARGV_LIST[i] = 1; }

    if ( strcmp( ARGV_LIST[i], "FASTDEBUG" ) == 0 ) 
      { INPUTS.FASTDEBUG = 1; USE_ARGV_LIST[i] = 1; }


    if ( strcmp( ARGV_LIST[i], "FLUXERR" ) == 0 ) 
      { INPUTS.FLUXERR_FLAG = 1; USE_ARGV_LIST[i] = 1; }
    

    if ( strcmp( ARGV_LIST[i], "SN_SED" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.inFile_snsed ); 
    }

    if ( strcmp( ARGV_LIST[i], "SN_SED_POWERLAW" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.snsed_powerlaw ); 
    }

    // ----------
    // check for primary override (Apr 2015)
    if ( strcmp( ARGV_LIST[i], "BD17_SED" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", tmpFile);
      primary_override("BD17", tmpFile);
    }
    if ( strcmp( ARGV_LIST[i], "VEGA_SED" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", tmpFile);
      primary_override("VEGA", tmpFile);
    }
    if ( strcmp( ARGV_LIST[i], "PRIMARY_SED" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", tmpName);
      i++ ; sscanf(ARGV_LIST[i] , "%s", tmpFile);
      primary_override(tmpName, tmpFile);
    }

    // ----------

    if ( strcmp( ARGV_LIST[i], "PLOTFLAG_SN" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PLOTFLAG_SN ); 
    }
    if ( strcmp( ARGV_LIST[i], "PLOTFLAG_PRIMARY" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.PLOTFLAG_PRIMARY ); 
    }

    if ( strcmp( ARGV_LIST[i], "READ_ZPOFF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.IRD_ZPOFF ); 
      if ( INPUTS.IRD_ZPOFF == 0 ) { NZPOFF = 0; }
    }


    if ( strcmp( ARGV_LIST[i], "NLAMBIN_FT" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.NLAMBIN_FT ); 
    }

    if ( strcmp( ARGV_LIST[i], "SN_TYPE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.SN_TYPE ); 
    }
    if ( strcmp( ARGV_LIST[i], "SN_MAGOFF" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.SN_MAGOFF ); 
    }

    if ( strcmp( ARGV_LIST[i], "SN_SED_FUDGE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.inFile_snsed_fudge ); 
    }

    if ( strcmp( ARGV_LIST[i], "TREF_EXPLODE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.TREF_EXPLODE ); 
    }

    if ( strcmp( ARGV_LIST[i], "REDSHIFT_RANGE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.REDSHIFT_MIN );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.REDSHIFT_MAX );
    }
    if ( strcmp( ARGV_LIST[i], "REDSHIFT_BINSIZE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.REDSHIFT_BINSIZE );
    }


    if ( strcmp( ARGV_LIST[i], "TREST_RANGE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &SNSED.TREST_MIN );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &SNSED.TREST_MAX );
    }


    if ( strcmp( ARGV_LIST[i], "SPECTROGRAPH" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%s", INPUTS.inFile_spectrograph );
    }

    
    // filter lam-shifts entered via command-line only
    // (input file only defines nominal filter-trans)
    if ( strcmp( ARGV_LIST[i], "FILTER_LAMSHIFT" ) == 0 ) {
      parse_FILTER_LAMSHIFT(&i);
    }
    if ( strcmp( ARGV_LIST[i], "DUPLICATE_LAMSHIFT_GLOBAL" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.LAMSHIFT_GLOBAL );
    }

    if ( strcmp( ARGV_LIST[i], "OPT_MWCOLORLAW" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%d", &INPUTS.OPT_MWCOLORLAW );
    }

    if ( strcmp( ARGV_LIST[i], "RV" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RV_MWCOLORLAW );
    }
    if ( strcmp( ARGV_LIST[i], "RV_MWCOLORLAW" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.RV_MWCOLORLAW );
    }

    if ( strcmp( ARGV_LIST[i], "AV_RANGE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.AV_MIN );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.AV_MAX );
    }
    if ( strcmp( ARGV_LIST[i], "AV_BINSIZE" ) == 0 ) {
      i++ ; sscanf(ARGV_LIST[i] , "%le", &INPUTS.AV_BINSIZE );
    }

    if ( strcmp( ARGV_LIST[i], "FILTER_OOB" ) == 0 ) { // out-of-band trans
      char BANDLIST[40];   double LAMRANGE[2], RATIO;
      i++ ; sscanf(ARGV_LIST[i] , "%s",  BANDLIST );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LAMRANGE[0] );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &LAMRANGE[1] );
      i++ ; sscanf(ARGV_LIST[i] , "%le", &RATIO );
      parse_OOB(BANDLIST, LAMRANGE, RATIO);
    }

  SKIPPY:

    if ( i > ilast ) {
      for ( iuse = ilast; iuse <= i; iuse++ ) 
	USE_ARGV_LIST[iuse] = 1;
    }

    i++ ;  ilast=i;

  }


} // end of kcor_input_override


// ************************************
void  primary_override(char *primName, char *fileName) {

  // Created Apr 20, 2015
  // 
  // Overwrite array PRIMARY array with command-line override.

  int  iprim, IPRIM  ;
  char oldFile[MXPATHLEN], *ptrFile ;

  IPRIM = -9 ;

  // first check of primName already exists:
  for(iprim=1; iprim <= INPUTS.NPRIMARY ; iprim++ ) {
    if ( strcmp(INPUTS.name_PRIMARY[iprim],primName) == 0 ) {
      ptrFile = INPUTS.inFile_PRIMARY[iprim];
      sprintf(oldFile, "%s", ptrFile);
      sprintf(ptrFile, "%s", fileName);
      printf("\n PRIMARY-OVERRIDE:\n    Replace %s file "
	     "'%s'\n\t\t with '%s'  (iprim=%d) \n\n",
	     primName, oldFile, ptrFile, iprim);
      IPRIM = iprim ;
    }
  }

  // if we get here then add to primary list
  if ( IPRIM < 0 ) {
    INPUTS.NPRIMARY++ ;
    IPRIM = INPUTS.NPRIMARY;  
    sprintf(INPUTS.name_PRIMARY[IPRIM],   "%s", primName);
    sprintf(INPUTS.inFile_PRIMARY[IPRIM], "%s", fileName);
    printf("\n PRIMARY-OVERRIDE: \n\t Add %s file '%s' \n\n",
	   primName, fileName);
  }


  // load another array.
  sprintf(PRIMARYSED[IPRIM].MAGSYSTEM_SEDFILE,"%s", fileName);

} // end of primary_override



// *********************************
void  parse_FILTER_LAMSHIFT(int *indx_ARGV) {

  // Nov 2014
  // parse INPUTS.FILTER_LAMSHIFT_STRING to get lambda shift
  // in each filter. This is designed for systematic tests
  //
  // *indx_ARGV is the index of the FILTER_LAMSHIFT key.
  // Note that *indx_ARGV is incremented here so that the
  // main parsing function can continue.
  //
  // Jan 28 2021: fix index bugs

  int i, ifilt, NBAND ;
  int LDMP = 0 ;
  double lamshift ;
  char *cfilt ;
  char fnam[] = "parse_FILTER_LAMSHIFT" ;

  // ------------ BEGIN --------------

  i = *indx_ARGV ;

  NBAND = 0 ;

  while ( i > 0  && i < NARGV_LIST-1) {

    if ( LDMP ) {
      printf(" xxx ------------------------------ \n");
      printf(" 1. xxx %s: ARG[%d of %d] = %s  \n", 
	     fnam, i, NARGV_LIST, ARGV_LIST[i] );    fflush(stdout);
    }

    i++ ; 
    cfilt = ARGV_LIST[i];   
    ifilt = INTFILTER_kcor(cfilt);

    if ( LDMP ) {
      printf(" 2. xxx %s: ARG[%d of %d] = %s  \n", 
	     fnam, i, NARGV_LIST, ARGV_LIST[i] );    fflush(stdout);
    }

    if ( ifilt > 0  ) {
      NBAND++ ;
      i++ ; sscanf( ARGV_LIST[i] , "%le", &lamshift ); 
      INPUTS.FILTER_LAMSHIFT[ifilt] = lamshift ;
      INPUTS.NFILTER_LAMSHIFT++ ;
    }
    else   {
      i-- ; 
      // bail when we get a string that is clearly not a filter name
      goto DONE ; 
    }

    printf("\t LAMSHIFT = %6.2f A  for filter = '%s' (%2d) \n",
	   lamshift, cfilt, ifilt );

    fflush(stdout);
  }


 DONE:
  *indx_ARGV = i ;

} // end of parse_FILTER_LAMSHIFT

// ****************************************
int INTFILTER_kcor(char *filterName) {

  // Created Nov 16 2014
  // return filter index for FILTER struct
  // Note that returned IFILT is a sparse filter index,
  // and is NOT the absolute index as in sntools function INTFILTER.
  //
  // First compare with full filter name (e.g., SDSS-g).
  // Then try just the single-char band name (e.g., 'g').
  // If the band appears more than once, abort.

  int ifilt, IFILT, NMATCH ;
  char fnam[] = "INTFILTER_kcor" ;
  // ------------- BEGIN -------------
  

  NMATCH = 0;  IFILT = -9 ;
  for(ifilt=1; ifilt <= NFILTDEF; ifilt++ ) {
    if ( strcmp(filterName,FILTER[ifilt].name) == 0  )
      { IFILT = ifilt;  NMATCH++ ; }
  }

  if ( NMATCH == 1 ) { goto DONE_CHECK ; }

  // ------------------------------------
  // if we get here, check single-char band.

  NMATCH = 0;
  for(ifilt=1; ifilt <= NFILTDEF; ifilt++ ) {
    if ( strcmp(filterName,FILTER[ifilt].band) == 0  )
      { IFILT = ifilt;  NMATCH++ ; }
  }

  // -----------

 DONE_CHECK:
  if ( NMATCH > 1 ) { 
    sprintf(c1err,"NMATCH=%d for filterName='%s'", 
	    NMATCH, filterName);
    sprintf(c2err,"Check filter names in input file: %s", 
	    INPUTS.inFile_input);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return IFILT ;

} // end of INTFILTER_kcor

// ********************************
void rd_ZPOFF(char *sdir, char *zpoff_file_override) {

  // Created Mar 25,  2010 by R.Kessler
  // Read ZPOFF.DAT file from *sdir and store AB offsets.
  // If zpoff_file_override is set, read this instead of ZPOFF.DAT
  // If ZPOFF.DAT does not exist, just return since this file is optional.
  //
  // Jan 2021; pass option zpoff_file_override

  double zpoff;
  FILE *fp;
  int gzipFlag, REQUIRE_ZPOFF_FILE ;
  char  cfilt[40], zpoff_File[MXPATHLEN], ZPOFF_FILE[MXPATHLEN] ;
  char  fnam[] = "rd_ZPOFF"     ;

  // ------------- BEGIN ------------

  if ( IGNOREFILE(zpoff_file_override) ) {
    REQUIRE_ZPOFF_FILE = 0;
    sprintf(zpoff_File,"%s/%s",  sdir, ZPOFF_FILE_DEFAULT);
  }
  else {
    // check user-define ZPOFF file (Jan 2021)
    REQUIRE_ZPOFF_FILE = 1 ; 
    // if there is a slash in the file name, assume it has full path.
    // Otherwise, glue on sdir to name.
    if ( strchr(zpoff_file_override,'/') == NULL ) 
      { sprintf(zpoff_File,"%s/%s",  sdir, zpoff_file_override); }
    else
      { sprintf(zpoff_File,"%s",  zpoff_file_override); }
  }


  fp = snana_openTextFile(0, PATH_SNDATA_FILTER, zpoff_File, 
			  ZPOFF_FILE, &gzipFlag );

  // if file does not exist ...
  if ( fp == NULL ) {
    if ( REQUIRE_ZPOFF_FILE ) {
      print_preAbort_banner(fnam);
      printf("\t Tried to open ZPOFF_FILE = \n %s\n", ZPOFF_FILE);
      sprintf(c1err,"User-requested ZPOFF file does not exist.");
      sprintf(c2err,"Check ZPOFF_FILE key");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    else {
      printf("\t Could not find default %s => no ZP offsets. \n", zpoff_File);
      return ;
    }
  }

  // file exists, so read it. 
  // Format is 'filterchar value', repeat for each filter

  printf("\t Found ZP offsets in %s  \n", ZPOFF_FILE );
  printf("\t BEWARE: ZP offsets are applied in analysis (not here)!\n");
  fflush(stdout);

  if ( INPUTS.IRD_ZPOFF == 0 ) { 
    printf("\t READ_ZPOFF=0 => ignore ZP offsets. \n" );
    return ; 
  }


  while( (fscanf(fp, "%s", cfilt )) != EOF) {
    readdouble(fp, 1, &zpoff );

    NZPOFF++;
    sprintf(ZPOFF_STORE.FILTLIST,"%s%s", 
	    ZPOFF_STORE.FILTLIST, cfilt );
    ZPOFF_STORE.FILTVALUES[NZPOFF] = zpoff ;
    ZPOFF_STORE.IFILTPATH[NZPOFF]  = NFILTPATH ;
  }

  fclose(fp);

} // end  of rd_ZPOFF


// ********************************
double get_ZPOFF(char *filtname, int IPATH) {

  // Inputs:
  // *filtname = full name of filter
  // IPATH     = index to path to distinguish same filter char in 
  //             different subdirs
  //
  // Dec 2012: major fix, pass IPATH to require IPATH match
  //           as well as filter-char match. For example, avoids setting
  //           SDSS-U offset to Bessell-U
  //

  double ZPOFF;
  int i, lnam, MATCH_PATH, MATCH_CHAR;
  char filtname1[2], ctmp1[2];

  // ------------ BEGIN ------------

  ZPOFF = 0.0;
  if( INPUTS.IRD_ZPOFF == 0 ) { return ZPOFF ; }

  // strip off last character of filtname to look for ZPOFF match
  lnam = strlen(filtname);
  sprintf(filtname1, "%c", *(filtname+lnam-1) );

  for ( i=1; i <= NZPOFF; i++ ) {
    sprintf(ctmp1, "%c", ZPOFF_STORE.FILTLIST[i-1] );

    MATCH_CHAR = ( strcmp(ctmp1,filtname1) == 0  );
    MATCH_PATH = ( IPATH == ZPOFF_STORE.IFILTPATH[i] ) ;

    if ( MATCH_CHAR && MATCH_PATH ) 
      { ZPOFF = ZPOFF_STORE.FILTVALUES[i]; }

  }

  return ZPOFF;

} // end of get_ZPOFF


// ******************************************
void  ADDFILTERS_LAMSHIFT_GLOBAL(void) {

  // May 2017
  // if LAMSHIFT_GLOBAL != 0 then DUPLICATE each filters
  // with this LAMSHIFT, and set flag for analysis to
  // pick either nominal or lam-shifted filters.
  // This allows running LAMSHIFT systematics (in snlc_fit)
  // with a single kcor file, rather than defining a separate
  // kcor file for each LAMSHIFTed passband.
  //
  // The duplicate filter name has an asterisk (*) in front;
  // e.g. SDSS-r duplicate is *SDSS-r
  //
  // May 22 2020: load MAGSYSTEM_INDX_INPUT

  int NF = NFILTDEF ;
  int ifilt, ifilt2, INDX ;
  double LAMSHIFT = INPUTS.LAMSHIFT_GLOBAL ;
  char fnam[] = "ADDFILTERS_LAMSHIFT_GLOBAL" ;

  // ------------- BEGIN ------------------

  if ( fabs(LAMSHIFT) < 0.0001 ) { return ; }

  if ( NKCOR > 0 ) {
    sprintf(c1err,"NKCOR=%d, but k-corrections are not", NKCOR);
    sprintf(c2err,"compatible with DUPLICATE_LAMSHIFT_GLOBAL=%.1f",
	    LAMSHIFT);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  for(ifilt=1; ifilt <= NF; ifilt++ ) {
    NFILTDEF++ ;
    ifilt2 = NF + ifilt; // index of duplicate
    FILTER[ifilt2].IFLAG_DUPLICATE = 1; 
    INPUTS.FILTER_LAMSHIFT[ifilt2] = LAMSHIFT ;
    INPUTS.NFILTER_LAMSHIFT++ ;

    sprintf ( FILTER[ifilt2].band, "%s" ,  FILTER[ifilt].band );
    sprintf ( FILTER[ifilt2].name, "*%s" , FILTER[ifilt].name );
    sprintf ( FILTER[ifilt2].file, "%s" ,  FILTER[ifilt].file );    

    FILTER[ifilt2].MAGSYSTEM_OFFSET = FILTER[ifilt].MAGSYSTEM_OFFSET ;

    INDX = FILTER[ifilt].MAGSYSTEM_INDX;
    FILTER[ifilt2].MAGSYSTEM_INDX        = INDX;
    FILTER[ifilt2].MAGSYSTEM_INDX_INPUT  = INDX; // May 22 2020

    sprintf(FILTER[ifilt2].MAGSYSTEM_NAME, "%s",FILTER[ifilt].MAGSYSTEM_NAME);
    sprintf(FILTER[ifilt2].FILTSYSTEM_NAME,"%s",FILTER[ifilt].FILTSYSTEM_NAME);
    sprintf(FILTER[ifilt2].PATH,           "%s",FILTER[ifilt].PATH );

    FILTER[ifilt2].MAGFILTER_REF    = FILTER[ifilt].MAGFILTER_REF ;
    FILTER[ifilt2].FILTSYSTEM_INDX  = FILTER[ifilt].FILTSYSTEM_INDX ;

    FILTER[ifilt2].IPATH              = FILTER[ifilt].IPATH ;
    FILTER[ifilt2].MAGFILTER_ZPOFF    = FILTER[ifilt].MAGFILTER_ZPOFF ;
    FILTER[ifilt2].SNMAG0_CHECK_VALUE = FILTER[ifilt].SNMAG0_CHECK_VALUE ;

    INDX = FILTER[ifilt].MAGSYSTEM_INDX; //   = MAGSYSTEM->INDX ;

    FILTER[ifilt2].IFLAG_SYN  = FILTER[ifilt].IFLAG_SYN ;
    FILTER[ifilt2].LAMBDA_MIN = FILTER[ifilt].LAMBDA_MIN ;
    FILTER[ifilt2].LAMBDA_MAX = FILTER[ifilt].LAMBDA_MAX ;
  }

} // end ADDFILTERS_LAMSHIFT_GLOBAL

// ********************************
int kcor_ini(void) {

/***

  - read SN SED, primary SED, and filter responses.
  - rebin primary to SN lambda bins.
 

  Nov 12, 2010: determine which filters are rest/obs frame

  Jul 11, 2011: add extra K_XX if BX filter is defined.

  Dec 8, 2011: remove RDFLAG loop over iprim

  Feb 11, 2013: check and set output formats.

 ***/
 
   int i, istat, iprim, ikcor, MSKTMP, ifilt, ifilt_rest, ifilt_obs ;
   int OVP_OBS, OVP_REST, lenf, LREST, LX     ;

   char cfilt[2] ;
   char fnam[] = "kcor_ini";

   /*  ----------------- BEGIN ------------------------- */

   // set kcor format for each file (Feb 2013)
   set_kcorFile_format();   

   printf("  ***** READ FILTER TRANSMISSIONS ***** \n" );

   FILTER_LAMBDA_MAX = -9999999.0;
   FILTER_LAMBDA_MIN = +9999999.0;

   for ( i=1;  i <= NFILTDEF ; i++ ) {
         istat = rd_filter(i);
   }
   fflush(stdout);

   // set min/max wavelength range to store SEDs
   set_store_lambda_range();

   if ( rd_snsed() != SUCCESS ) { return ERROR; }


   // be careful to loop over each 'iprim' in each stage so that
   // it works regardless of the order that PRIMARY refs are defined
   // in the kcor-input file.

   // read & rebin all USED primary refs.   
   for( iprim = 1; iprim <= INPUTS.NPRIMARY; iprim++ ) {
     if ( PRIMARYSED[iprim].USE == 0 ) { continue ; }
     istat = rd_primary(iprim, subdir_standards );
     if ( istat != SUCCESS ) { return ERROR; }
   }
 
   // compute primary mags and zp per filter
   for( iprim = 1; iprim <= INPUTS.NPRIMARY; iprim++ )
     { primarymag_zp(iprim) ; }

   // check option to transform primary; e.g., VEGA->AB
   for( iprim = 1; iprim <= INPUTS.NPRIMARY; iprim++ )
     { primarymag_zp2(iprim) ; }


   // print summary for each primary and filter
   for( iprim = 1; iprim <= INPUTS.NPRIMARY; iprim++ )
     { primarymag_summary(iprim) ; }
   

   // loop over Kcors to see which filters are rest/obs
   //
   for ( ikcor = 1; ikcor <= NKCOR; ikcor++ ) {
     index_filter ( ikcor, &ifilt_rest, &ifilt_obs );
     MSKTMP = FILTER[ifilt_rest].MASKFRAME ;
     FILTER[ifilt_rest].MASKFRAME = ( MSKTMP | MSKFRAME_REST );
     MSKTMP = FILTER[ifilt_obs].MASKFRAME ;
     FILTER[ifilt_obs].MASKFRAME = ( MSKTMP | MSKFRAME_OBS );
   }

   // now figure out rest-frame [non-obsframe] filters and make up
   // bogus K-corretions so that we get obs-mag for the
   // rest-frame filters (needed in snana at z=0).
   // Note that these extra KCOR tables are NOT written out,
   // so only the MAGOBS tables are kept.

   NKCOR_EXTRA = 0;

   for ( ifilt=1; ifilt <= NFILTDEF; ifilt ++ ) {
     MSKTMP = FILTER[ifilt].MASKFRAME ;

     OVP_OBS  = (MSKTMP & MSKFRAME_OBS ) ;
     OVP_REST = (MSKTMP & MSKFRAME_REST) ;

     lenf = strlen(FILTER[ifilt].name);
     sprintf(cfilt,"%c", FILTER[ifilt].name[lenf-1] );

     LREST = ( OVP_OBS == 0 && OVP_REST > 0 ) ;
     LX    = ( strcmp(cfilt,"X") == 0 && NKCOR > 0 ) ;

     if ( LREST || LX ) {

       NKCOR_EXTRA++;
       ikcor = NKCOR + NKCOR_EXTRA;  // absolute KCOR index
       if ( ikcor > MXKCOR ) {
	 sprintf(c1err,"ikcor = %d exceeds bound.", ikcor);
	 sprintf(c2err,"NKCOR=%d  NKCOR_EXTRA=%d", NKCOR, NKCOR_EXTRA);
	 errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
       }
       sprintf(KCORLIST[ikcor][0],"%s", FILTER[ifilt].name );
       sprintf(KCORLIST[ikcor][1],"%s", FILTER[ifilt].name );
       sprintf(KCORSYM[ikcor],"K_%s%s", cfilt, cfilt );     

     }
   }

   printf("\n Finished  %s  \n\n", fnam );

   fflush(stdout);
   return SUCCESS;
}

// ***************************************************
void  set_store_lambda_range(void) {

  // Created Nov 3 2021 by R.Kessler
  // Set wave range to store SEDs.
  // Default is min/max wavelength of bluest/reddest filters.
  // If spectrograph goes bluer/redder than filters, store
  // extended wvae range.
  //
  // Output is global STORE_LAMBDA_MIN and STORE_LAMBDA_MAX
  //
  // TEMP: code still uses FILTER_LAMBDA_MIN[MAX] until new
  //      STORE_LAMBDA_MIN[MAX] are verified.
  //

  int  LEGACY = 0;  // set True to restore using FILTER_LAMBDA_MIN[MAX]
  char fnam[] = "set_store_lambda_range" ;

  // ------------ BEGIN ------------

  printf("\n");
  printf(" Global wave range of all filters: %d to %d A \n",
	 (int)FILTER_LAMBDA_MIN, (int)FILTER_LAMBDA_MAX );

  STORE_LAMBDA_MIN = FILTER_LAMBDA_MIN;
  STORE_LAMBDA_MAX = FILTER_LAMBDA_MAX;

  if ( SPECTROGRAPH_USEFLAG && !LEGACY ) {
    printf(" Global wave range of spectrograph: %d to %d A\n",
	   (int)INPUTS_SPECTRO.LAM_MIN, (int)INPUTS_SPECTRO.LAM_MAX);

    if ( INPUTS_SPECTRO.LAM_MIN < STORE_LAMBDA_MIN ) 
      { STORE_LAMBDA_MIN = INPUTS_SPECTRO.LAM_MIN; }

    if ( INPUTS_SPECTRO.LAM_MAX > STORE_LAMBDA_MAX ) 
      { STORE_LAMBDA_MAX = INPUTS_SPECTRO.LAM_MAX; }
  }
  
  printf(" Final wavelength storage range: %d to %d A \n",
	 (int)STORE_LAMBDA_MIN, (int)STORE_LAMBDA_MAX );
  
  fflush(stdout);

  return;

} // end set_store_lambda_range

// ***************************************************
void set_kcorFile_format(void) {

  // Created Feb 11, 2013
  // set INPUTS.FORMAT_FLAG[ifile] for each OUTFILE.
  // If CERNLIB flag is not set and hbook file is requested, ABORT.

  int  i, FLAG ;
  char *SUFFIX, *suffix, *ptrFile;
  char CFORMAT[4][20] ;
  char fnam[] = "set_kcorFile_format" ;

  // ------------- BEGIN ---------


  sprintf( CFORMAT[FORMAT_FITS],  "FITS");
  sprintf( CFORMAT[FORMAT_HBOOK], "HBOOK");

  for ( i=0; i < INPUTS.N_OUTFILE; i++ ) {

    ptrFile = INPUTS.OUTFILE[i] ;
    SUFFIX  = strstr(ptrFile, ".HIS");
    suffix  = strstr(ptrFile, ".his");
    if ( suffix == NULL && SUFFIX == NULL ) 
      { FLAG = FORMAT_FITS ; }
    else
      { FLAG = FORMAT_HBOOK ; }

    INPUTS.FORMAT_FLAG[i] = FLAG ;

    printf("\t OUTFILE: %s  ->  %s format\n", ptrFile, CFORMAT[FLAG] );


    // if CERNLIB pre-processor flag is NOT set and hbook is requested,
    // then abort with error message.
#ifndef CERNLIB
    if ( FLAG == FORMAT_HBOOK ) {
      sprintf(c1err,"HBOOK format no longer supported.");
      sprintf(c2err,"Use .fits (.FITS) extension for OUTFILE argument.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
    }
#endif

  }

  printf("\n");

} // end of set_kcorFile_format


// ***************************************************
int  malloc_ini(void) {

  // Created Mar 2011 
  // allocate multi-dimensional arrays to avoid 
  // large static declarations leading to 1.6 GB program size.
  // Now the program size is well under 1 GB.
  //
  // Oct 25, 2011: replace hard-wired i4 and i8 with
  //               sizeof(float*) and sizeof(double*)
  //               to work on both 32-bit and 64-bit
  //
  //
  // Nov 11, 2011: 
  //   - set NBIN_FILT = NFILTDEF instead of MXFILTDEF
  //   - sizeof(float*) -> sizeof(float)
  //     and same with double.
  //
  // Feb 13, 2013: fix subtle bug that screws 64-bit 
  //               size(float) -> sizeof(float *)
  //               to get size of pointer. See i4 and i8 calc.
  //
  // ----------------------------------------------

  int i4, i4p, i8, i8p ;
  int i_ebv, i_filt, i_av, i_z, i_ep, i_kcor;
  int NBIN_EP, NBIN_Z, NBIN_AV, NBIN_FILT, NBIN_KCOR ;
  int MEM_R4MAG, MEM_MW, MEM_KCOR ;

  // ------------- BEGIN -------------

  if ( NKCOR <= 0 ) { return SUCCESS ; }

  i8  = sizeof(double  ) ;
  i8p = sizeof(double *) ;
  i4  = sizeof(float   ) ; 
  i4p = sizeof(float  *) ;  // i4p != i4  on 64 bit

  NBIN_Z    = INPUTS.NBIN_REDSHIFT ;
  NBIN_AV   = INPUTS.NBIN_AV ;
  NBIN_EP   = SNSED.NEPOCH ;
  NBIN_FILT = NFILTDEF ;
  NBIN_KCOR = NKCOR + NKCOR_EXTRA ;

  printf(" Allocate large memory chunks: \n");
  printf("\t NBIN(KCOR)   = %d + %d \n", NKCOR, NKCOR_EXTRA );
  printf("\t NBIN(AV)     = %d \n", NBIN_AV );
  printf("\t NBIN(Z)      = %d \n", NBIN_Z );
  printf("\t NBIN(EPOCH)  = %d \n", NBIN_EP );
  printf("\t NBIN(MWEBV)  = %d \n", MXMWEBV+1 );
  printf("\t NBIN(FILT)   = %d \n", NBIN_FILT );
  printf("\t sizeof(float,float*   double,double*) = %d,%d  %d,%d\n", 
	 i4, i4p,  i8,i8p );
  fflush(stdout);

  MEM_R4MAG = 0 ;
  MEM_MW    = 0 ;
  MEM_KCOR  = 0 ;

  // start wit R4MAG_OBS ...

  SNSED.R4MAG_OBS = (float *****)malloc(i4p*(MXMWEBV+1) );  // 0
  for ( i_ebv=0; i_ebv <= MXMWEBV; i_ebv++ ) {    // 1
    SNSED.R4MAG_OBS[i_ebv] = (float ****)malloc(i4p*(NBIN_FILT+1));

    for ( i_filt=0; i_filt <= NBIN_FILT; i_filt++ ) {  // 2
      SNSED.R4MAG_OBS[i_ebv][i_filt] = (float ***)malloc(i4p*(NBIN_AV+1));

      for ( i_av=0; i_av <= NBIN_AV; i_av++ ) {      // 3
	SNSED.R4MAG_OBS[i_ebv][i_filt][i_av] = 
	  (float **)malloc(i4p*(NBIN_Z+1));

	for ( i_z=0; i_z <= NBIN_Z; i_z++ ) {          // 4
	  SNSED.R4MAG_OBS[i_ebv][i_filt][i_av][i_z] = 
	    (float *)malloc(i4*(NBIN_EP+1));

	  for ( i_ep=0; i_ep <= NBIN_EP; i_ep++ ) {
	    SNSED.R4MAG_OBS[i_ebv][i_filt][i_av][i_z][i_ep] = NULLVAL ;
	    MEM_R4MAG += i4 ;
	  } // i_ep

	}  // i_z
      } // i_av
    } // i_filt
  } // i_ebv

  // --------------------------------------------
  // Now allocate 4-dim array for MW_dXT_dEBV

  SNSED.MW_dXT_dEBV = (double ****)malloc(i8p*(NBIN_FILT+1));   // 1

  for ( i_filt=0; i_filt <= NBIN_FILT; i_filt++ ) {  // 2
    SNSED.MW_dXT_dEBV[i_filt] = (double ***)malloc(i8p*(NBIN_AV+1));

    for ( i_av=0; i_av <= NBIN_AV; i_av++ ) {      // 3
      SNSED.MW_dXT_dEBV[i_filt][i_av] = (double **)malloc(i8p*(NBIN_Z+1));

      for ( i_z=0; i_z <= NBIN_Z; i_z++ ) {          // 4
	SNSED.MW_dXT_dEBV[i_filt][i_av][i_z] = (double*)malloc(i8*(NBIN_EP+1));

	for ( i_ep=0; i_ep <= NBIN_EP; i_ep++ ) {
	  SNSED.MW_dXT_dEBV[i_filt][i_av][i_z][i_ep] = NULLVAL ;
	  MEM_MW += i8 ;
	} // i_ep
	
      }  // i_z
    } // i_av
  } // i_filt


  // ----------------
  // R4KCOR_GRID[MXKCOR][MXAV][MXREDSHIFT][MXEP] ;

  R4KCOR_GRID.VALUE    = (float****)malloc(i4p*(NBIN_KCOR+1));   // 1  
  R4KCOR_GRID.REDSHIFT = (float****)malloc(i4p*(NBIN_KCOR+1));   // 1  
  R4KCOR_GRID.EPOCH    = (float****)malloc(i4p*(NBIN_KCOR+1));   // 1  

  

  for ( i_kcor=0; i_kcor <= NBIN_KCOR; i_kcor++ ) {  // 2
    R4KCOR_GRID.VALUE[i_kcor]    = 
      (float***)malloc(i4p*(NBIN_KCOR+1));   // 1  
    R4KCOR_GRID.REDSHIFT[i_kcor] = 
      (float***)malloc(i4p*(NBIN_KCOR+1));   // 1  
    R4KCOR_GRID.EPOCH[i_kcor]    = 
      (float***)malloc(i4p*(NBIN_KCOR+1));   // 1  

    for ( i_av=0; i_av <= NBIN_AV; i_av++ ) {      // 3
      R4KCOR_GRID.VALUE[i_kcor][i_av]    = 
	(float**)malloc(i4p*(NBIN_Z+1));
      R4KCOR_GRID.REDSHIFT[i_kcor][i_av] = 
	(float**)malloc(i4p*(NBIN_Z+1));
      R4KCOR_GRID.EPOCH[i_kcor][i_av]    = 
	(float**)malloc(i4p*(NBIN_Z+1));

      for ( i_z=0; i_z <= NBIN_Z; i_z++ ) {          // 4

	R4KCOR_GRID.VALUE[i_kcor][i_av][i_z] = 
	  (float*)malloc(i4*(NBIN_EP+1));
	R4KCOR_GRID.REDSHIFT[i_kcor][i_av][i_z] = 
	  (float*)malloc(i4*(NBIN_EP+1));
	R4KCOR_GRID.EPOCH[i_kcor][i_av][i_z] = 
	  (float*)malloc(i4*(NBIN_EP+1));

	for ( i_ep=0; i_ep <= NBIN_EP; i_ep++ ) {
	  R4KCOR_GRID.VALUE[i_kcor][i_av][i_z][i_ep]    = NULLVAL ;
	  R4KCOR_GRID.REDSHIFT[i_kcor][i_av][i_z][i_ep] = NULLVAL ;
	  R4KCOR_GRID.EPOCH[i_kcor][i_av][i_z][i_ep]    = NULLVAL ;
	  MEM_KCOR += 3*i4p ;
	} // i_ep
	
      }  // i_z
    } // i_av
  } // i_kcor

  

  // --------
  // print summary of memory allocations

  printf("\t Allocated %6.2f MB to SNSED.R4MAG_OBS \n",
	 1.0E-6 * (float)MEM_R4MAG );
  printf("\t Allocated %6.2f MB to SNSED.MW_dXT_dEBV \n",
	 1.0E-6 * (float)MEM_MW );
  printf("\t Allocated %6.2f MB to R4KCOR_GRID \n",
	 1.0E-6 * (float)MEM_KCOR );


  //  debugexit("malloc"); // xxxxxx

  fflush(stdout);
  return SUCCESS ;

} // end of malloc_ini

// *********************************
int kcor_out(void) {

  int i, FLAG ;
  char *ptrFile ;

  // check hbook option

  for ( i=0; i < INPUTS.N_OUTFILE; i++ ) {
   
    ptrFile = INPUTS.OUTFILE[i] ;
    FLAG    = INPUTS.FORMAT_FLAG[i] ;

    if ( FLAG == FORMAT_FITS  ) { wr_fits(ptrFile);  }  // fits option
  }

  // check option to text-dump SN mags.
  wrsnmag_text();

  return SUCCESS;
}  // end of kcor_out






// ******************************************
int rd_filter ( int ifilt ) {

/****
   - Open filter transmission file corresponding to "ifilt"
   - Read filter response from FILE.
   - File must have two columns: wavelen (A) and transmission
   - binning can be non-uniform
   - Load FILTER structure.

   - compute LAMAVG (Feb 15, 2006)

  May 2014: fix dlam calc that was removed back on Sep 12 2013.
            used only for internal SUMTOT --> no harm in previously
            generated kcor files.

  Nov 16 2014: add  INPUTS.FILTER_LAMSHIFT[ifilt] to each lambda value
                 (intended for systematic tests). LAMSHIFT values are 
                 entered via command-line  key "FILTER_LAMSHIFT ..."

  July 2016: check IFU option to strip off spectrograph info.

  Apr 14 2020
    Fix calc of fitler-averages to account for non-uniform binning.
    In particular, compute dlam for each wavelength bin.
 
  May 7 2020: 
    fix refactor bug from Apr 14 2020; LAMSHIFT now applied in 
    separate ilam loop.

 ****/

   double LAMSHIFT =  INPUTS.FILTER_LAMSHIFT[ifilt] ; 
   FILE *fp;

   double lambda, trans, tsum, wtsum ;
   double lammin, lammax, lamavg, dlam, dlam_last     ;

   int NBIN, IFLAG_SYN, ilam, gzipFlag ;
     
   char  *ptr_name, *ptr_file, txtFilter[60] ;
   char  FILTFILE_FULLNAME[MXCHAR_FILENAME] ;
   char  fnam[] = "rd_filter"    ;

  /* -------------------- BEGIN ------------------ */

   ptr_file   = FILTER[ifilt].file ;
   ptr_name   = FILTER[ifilt].name ;
   IFLAG_SYN  = FILTER[ifilt].IFLAG_SYN ;

   if ( IFLAG_SYN ) {
     double LMIN = FILTER[ifilt].LAMBDA_MIN ;
     double LMAX = FILTER[ifilt].LAMBDA_MAX ;    
     get_FILTERtrans_spectrograph(&LMIN, &LMAX, MXLAM_FILT, &NBIN, 
				  &FILTER[ifilt].LAMBDA[1], 
     				  &FILTER[ifilt].TRANS[1] ) ;

     FILTER[ifilt].LAMBDA_MIN = LMIN ;  // note that LMIN, LMAX are
     FILTER[ifilt].LAMBDA_MAX = LMAX ;  // are slightly extended
     sprintf(txtFilter, "%s", INPUTS_SPECTRO.INSTRUMENT_NAME);
   }
   else {
     // read 2-column ascii file
     int OPENMASK = OPENMASK_VERBOSE + OPENMASK_IGNORE_DOCANA ;
     fp = snana_openTextFile(OPENMASK,PATH_SNDATA_FILTER, ptr_file, 
			     FILTFILE_FULLNAME, &gzipFlag );
   
     // if we get here, abort because file cannot be found.
     if ( fp == NULL ) {
       sprintf ( c1err, "Cannot find file : %s", ptr_file );
       sprintf ( c2err," Check local dir and %s", PATH_SNDATA_FILTER );
       errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
     }     
     fclose(fp);
     rd2columnFile(FILTFILE_FULLNAME, MXLAM_FILT, &NBIN 
		   ,&FILTER[ifilt].LAMBDA[1]
		   ,&FILTER[ifilt].TRANS[1], 0 );

     sprintf(txtFilter,"Filter");
   }


   // --------------------------
   
   FILTER[ifilt].NBIN_LAMBDA = NBIN ;

   if ( NBIN < 1 ) {
     sprintf(c1err,"Found only %d bin for ifilt=%d (%s)", 
	     NBIN, ifilt, ptr_name );
     sprintf(c2err,"File='%s'", ptr_file);
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
   }
   
   if ( FILTER[ifilt].OOB.TRANS_MAX_RATIO > 0.0 ) {
     addOOBTrans_filter(ifilt);  // add optional out-of-band transmission
     NBIN = FILTER[ifilt].NBIN_LAMBDA ;
   }

   // May 20 2021: check option to remove OOB
   if ( INPUTS.LEAKAGE_CUT > 1.0E-12 ) {
     cutOOBTrans_filter(ifilt);
     NBIN = FILTER[ifilt].NBIN_LAMBDA ;
   }

   tsum = 0.0;
   wtsum = 0.0;
   dlam_last = 0.0 ;
   dlam      = 0.0 ;

   if ( LAMSHIFT != 0.0 ) {
     for ( ilam = 1; ilam <= NBIN; ilam++ ) 
       { FILTER[ifilt].LAMBDA[ilam] += LAMSHIFT ; }
   }

   // - - - -
   for ( ilam = 1; ilam <= NBIN; ilam++ ) {

     lambda = FILTER[ifilt].LAMBDA[ilam];
     trans  = FILTER[ifilt].TRANS[ilam] ;

     if ( lambda < SNSED.LAMBDA_MIN || lambda > SNSED.LAMBDA_MAX ) {
       sprintf(c1err,"lambda(%s) = %6.0f is outside SED-range of", 
	       ptr_name, lambda);
       sprintf(c2err,"%6.0f - %6.0f : check LAMBDA_RANGE: keyword",
	       SNSED.LAMBDA_MIN, SNSED.LAMBDA_MAX );
       errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
     }

      // divide-by-lambda for energy-based filter system
     if ( FILTER[ifilt].FILTSYSTEM_INDX == FILTSYSTEM_ENERGY ) {
	trans = trans * 1000. / lambda ;
	FILTER[ifilt].TRANS[ilam]  = trans ;
     }

     if ( ilam == 1 ) 
       { dlam = FILTER[ifilt].LAMBDA[ilam+1] - FILTER[ifilt].LAMBDA[ilam]; }
     else if ( ilam == NBIN ) 
       { dlam = FILTER[ifilt].LAMBDA[NBIN] - FILTER[ifilt].LAMBDA[NBIN-1]; }
     else {
       dlam = FILTER[ifilt].LAMBDA[ilam+1] - FILTER[ifilt].LAMBDA[ilam-1]; 
       dlam /= 2.0 ;
     }

      wtsum += dlam * trans * lambda ;
      tsum  += dlam * trans ;

   }  // end of fscanf loop   


   FILTER[ifilt].LAMAVG = wtsum / tsum ;
   FILTER[ifilt].SUMTOT = tsum ; // check later that filter covers SED

   lammin  = FILTER[ifilt].LAMBDA[1] ;
   lammax  = FILTER[ifilt].LAMBDA[NBIN] ;
   lamavg  = FILTER[ifilt].LAMAVG ;   

   FILTER[ifilt].LAMBDA_MIN     = lammin;
   FILTER[ifilt].LAMBDA_MAX     = lammax;

   if ( lammin < FILTER_LAMBDA_MIN ) { FILTER_LAMBDA_MIN = lammin ; }
   if ( lammax > FILTER_LAMBDA_MAX ) { FILTER_LAMBDA_MAX = lammax ; }
   

   printf("\t Filter-%2.2d %10s : "
	  "%3d values from %5.0f to %5.0f A, <lam>=%5.0f \n", 
	  ifilt, ptr_name, NBIN, lammin, lammax, lamavg );

   if ( LAMSHIFT  != 0.0 ) {
     printf("\t\t\t LAMSHIFT = %.1f \n", LAMSHIFT );
   }

   fflush(stdout);  
   
   return SUCCESS;

}  // end of rd_filter function


 // *************************************
void addOOBTrans_filter(int ifilt) {

  // Feb 2019
  // Add out-of-band (OOB) transmission to filter 'ifilt'.
  // It's a bit tricky because the wave-range is extended
  // with more wave bins.

  double OOB_RATIO  = FILTER[ifilt].OOB.TRANS_MAX_RATIO ;
  double OOB_MINLAM = FILTER[ifilt].OOB.LAMRANGE[0] ;
  double OOB_MAXLAM = FILTER[ifilt].OOB.LAMRANGE[1] ;
  char  *NAME       = FILTER[ifilt].name; 
  int    NBIN_LAM_ORIG = FILTER[ifilt].NBIN_LAMBDA ;
  double MINLAM_ORIG   = FILTER[ifilt].LAMBDA[1] ;
  double MAXLAM_ORIG   = FILTER[ifilt].LAMBDA[NBIN_LAM_ORIG] ;
  double LAM_ORIG[MXLAM_FILT], TRANS_ORIG[MXLAM_FILT];

  double TransMax = 0.0 ;
  double trans, lam, OOB_TRANS, xlam;
  double LAMBIN_BLUE, LAMBIN_RED, MINLAM_NEW, MAXLAM_NEW ;
  int    NBIN_LAM_NEW, ilam, ilam_orig, NBINADD_BLUE, NBINADD_RED ;
  char fnam[] = "addOOBTrans_filter" ;
  int LDMP = 0 ; 

  // ----------- BEGIN ------------

  if ( OOB_RATIO == 0.0 ) { return; }

  // first find max trans
  for(ilam=1; ilam <= NBIN_LAM_ORIG; ilam++ ) {
    lam   = FILTER[ifilt].LAMBDA[ilam];
    trans = FILTER[ifilt].TRANS[ilam];
    if ( trans > TransMax ) { TransMax = trans; }
    LAM_ORIG[ilam]   = lam;
    TRANS_ORIG[ilam] = trans ;
  }

  // get lam binsize on blue edge
  LAMBIN_BLUE = FILTER[ifilt].LAMBDA[2] - FILTER[ifilt].LAMBDA[1] ;
  LAMBIN_RED  = FILTER[ifilt].LAMBDA[NBIN_LAM_ORIG] - 
    FILTER[ifilt].LAMBDA[NBIN_LAM_ORIG-1] ;

  
  NBINADD_BLUE = (int)((MINLAM_ORIG - OOB_MINLAM)/LAMBIN_BLUE) ;  
  if ( NBINADD_BLUE < 0 ) { NBINADD_BLUE = 0 ; }
  MINLAM_NEW = MINLAM_ORIG - LAMBIN_BLUE*(double)NBINADD_BLUE ;

  NBINADD_RED = (int)((OOB_MAXLAM-MAXLAM_ORIG)/LAMBIN_RED) ; 
  if ( NBINADD_RED < 0 ) { NBINADD_RED = 0 ; }
  MAXLAM_NEW = MAXLAM_ORIG + LAMBIN_RED*(double)NBINADD_RED ;

  NBIN_LAM_NEW = NBIN_LAM_ORIG + NBINADD_BLUE + NBINADD_RED ;

  if ( LDMP ) {
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx TransMax(%s) = %3f  LAMBIN[BLUE,RED]=%.1f,%.1f\n",
	   NAME, TransMax, LAMBIN_BLUE, LAMBIN_RED );
    printf(" xxx MINLAM = %.1f -> %.1f \n", MINLAM_ORIG, MINLAM_NEW);   
    printf(" xxx MAXLAM = %.1f -> %.1f \n", MAXLAM_ORIG, MAXLAM_NEW);
    printf(" xxx NBINLAM = %d -> %d \n", NBIN_LAM_ORIG, NBIN_LAM_NEW);
    fflush(stdout);
  }

  if ( NBIN_LAM_NEW > MXLAM_FILT ) {
    sprintf(c1err,"NBIN_LAM_NEW=%d exceeds MXLAM_FILT=%d",
	    NBIN_LAM_NEW, MXLAM_FILT);
    sprintf(c2err,"NBIN_LAM[ORIG, addBlue, addRed] = %d, %d, %d ",
	    NBIN_LAM_ORIG, NBINADD_BLUE, NBINADD_RED);
    errmsg(SEV_FATAL, 0, fnam,  c1err, c2err); 
  }


  for(ilam=1; ilam <= NBIN_LAM_NEW; ilam++ ) {

    ilam_orig = ilam - NBINADD_BLUE ;
    trans     = OOB_TRANS; 

    if ( ilam <= NBINADD_BLUE ) { 
      xlam  = (double)ilam; 
      lam   = MINLAM_NEW + LAMBIN_BLUE*(xlam-1.0);  
    }
    else if ( ilam <= (NBIN_LAM_ORIG+NBINADD_BLUE) ) {
      if ( ilam_orig < 1 || ilam_orig > NBIN_LAM_ORIG ) {
	sprintf(c1err,"ilam_orig=%d outside expected range %d to %d",
		ilam_orig, 1, NBIN_LAM_ORIG);
	sprintf(c2err,"Check index logic.");
	errmsg(SEV_FATAL, 0, fnam,  c1err, c2err); 
      }
      lam    = LAM_ORIG[ilam_orig]; 
      trans += TRANS_ORIG[ilam_orig] ; 
    }
    else {
      xlam = (double)(ilam - (NBIN_LAM_ORIG+NBINADD_BLUE) );
      lam = MAXLAM_ORIG + LAMBIN_RED*(xlam);  
    }

    FILTER[ifilt].LAMBDA[ilam] = lam ;
    FILTER[ifilt].TRANS[ilam]  = trans ;
    
  } // end ilam loop

  FILTER[ifilt].NBIN_LAMBDA = NBIN_LAM_NEW ;

  // write out filter trans file for crosschecks
  char filterFile[MXPATHLEN];
  FILE *fp ;
  sprintf(filterFile, "%s+OOB.dat", NAME);  
  fp = fopen(filterFile, "wt") ;
  printf("\t Write OOB trans-file: %s \n", filterFile);
  for(ilam=1; ilam <= NBIN_LAM_NEW; ilam++ ) {
    fprintf(fp,"%9.2f  %.4le\n", 
	    FILTER[ifilt].LAMBDA[ilam], FILTER[ifilt].TRANS[ilam] ) ;
  }
  fclose(fp);

  
  return;

} // end addOOBTrans_filter


// *************************************
void cutOOBTrans_filter(int ifilt) {

  // May 2021
  // Remove OOB usage input LEAKAGE_CUT << 1.
  // LEAKAGE_CUT is applied to trans/Transmax.
  //
  // It's a bit tricky because the wave-range is reduced
  // and uniformity must be preserved.

  double LEAKAGE_CUT   = INPUTS.LEAKAGE_CUT ;
  char  *NAME          = FILTER[ifilt].name; 
  int    NBIN_LAM_ORIG = FILTER[ifilt].NBIN_LAMBDA ;
  double MINLAM_ORIG   = FILTER[ifilt].LAMBDA[1] ;
  double MAXLAM_ORIG   = FILTER[ifilt].LAMBDA[NBIN_LAM_ORIG] ;
  double LAM_ORIG[MXLAM_FILT], TRANS_ORIG[MXLAM_FILT];

  double TransMax = 0.0 ;
  double trans, lam, OOB_TRANS, xlam;
  double MINLAM_NEW, MAXLAM_NEW ;
  int    NBIN_LAM_NEW, ilam, ilam_new, ilam_orig, NBINCUT_BLUE, NBINCUT_RED ;
  bool   LCUT, FOUND_TRANSMAX=false;
  char fnam[] = "cutOOBTrans_filter" ;
  int LDMP = 0 ; 

  // ----------- BEGIN ------------

  if ( LEAKAGE_CUT == 0.0 ) { return; }

  // first find max trans
  for(ilam=1; ilam <= NBIN_LAM_ORIG; ilam++ ) {
    lam   = FILTER[ifilt].LAMBDA[ilam];
    trans = FILTER[ifilt].TRANS[ilam];
    if ( trans > TransMax ) { TransMax = trans; }
    LAM_ORIG[ilam]   = lam;
    TRANS_ORIG[ilam] = trans ;
  }

  ilam_new=0;

  NBINCUT_BLUE = NBINCUT_RED = 0 ;
  MINLAM_NEW   = 9.0E9; MAXLAM_NEW = 0.0 ;

  for(ilam_orig=1; ilam_orig <= NBIN_LAM_ORIG; ilam_orig++ ) {

    trans = TRANS_ORIG[ilam_orig] / TransMax ;
    lam   = LAM_ORIG[ilam_orig] ;

    // xxx    if ( trans > 0.999 ) { FOUND_TRANSMAX = true ; }
    LCUT = ( trans < LEAKAGE_CUT );
    if ( LCUT ) { continue; }

    ilam_new++ ;
    FILTER[ifilt].LAMBDA[ilam_new] = lam ;
    FILTER[ifilt].TRANS[ilam_new]  = trans ;

    if ( lam < MINLAM_NEW ) { MINLAM_NEW = lam; }
    if ( lam > MAXLAM_NEW ) { MAXLAM_NEW = lam; }

  } // end ilam loop

  NBIN_LAM_NEW = ilam_new;
  FILTER[ifilt].NBIN_LAMBDA = NBIN_LAM_NEW ;

  printf("  Remove %s OOB leakage: WAVE-RANGE = [%.1f,%.1f] -> [%.1f,%.1f]\n",
	 NAME, MINLAM_ORIG, MAXLAM_ORIG, MINLAM_NEW, MAXLAM_NEW);

  /* xxxx
  // write out filter trans file for crosschecks
  char filterFile[MXPATHLEN];
  FILE *fp ;
  sprintf(filterFile, "TMP_%s-OOB.dat", NAME);  
  fp = fopen(filterFile, "wt") ;
  printf("\t Write trans-file: %s \n", filterFile);
  for(ilam=1; ilam <= NBIN_LAM_NEW; ilam++ ) {
    fprintf(fp,"%9.2f  %.4le\n", 
	    FILTER[ifilt].LAMBDA[ilam], FILTER[ifilt].TRANS[ilam] ) ;
  }
  fclose(fp);
  xxxxxx */
  
  return;

} // end cutOOBTrans_filter


// ****************************************************
int rd_snsed ( void ) {

  /***
   Read SN template spectra from file.
   Input format is

    day  lambda(A)  flux(erg/s/cm^2/A)

  Jan 16 2017: allow missing SN-SED file if NKCOR=0

  July 28 2018: call ENVreplace(INPUTS.inFile_snsed)

  ****/

   FILE *fp;

   double POWERLAW = INPUTS.snsed_powerlaw; 
   double
     day, dayoff, daymin, daymax
     ,DAYrange[2], LAMrange[2], DAYREAD[MXEP], LAMREAD[MXLAM_SN]
     ,DAYSTEP, LAMSTEP, LAMMIN, LAMMAX, wfilt, wflux, trans
     ,flux, lam, fnu, fcount
     ,FMIN = 1.0E-19
     ;

   int NDAY, NLAM,ilam0, iep0, jflux0, ilam1, iep1, ifilt, ilam, iep  ;
   int FOUND_SNSEDFILE, nflux_nan ;
   char sedFile[MXPATHLEN], sedcomment[40], SNPATH[MXPATHLEN] ;
   char fnam[] = "rd_snsed" ;

   //   --------------------- BEGIN --------------------------

   printf("\n\n  ***** READ Supernova SED TEMPLATES ***** \n" );

   // ===================================
   // first check for sedFile in local area:

   FOUND_SNSEDFILE = 0 ;

   ENVreplace(INPUTS.inFile_snsed,fnam,1);
   sprintf(sedFile,"%s", INPUTS.inFile_snsed );
   if ( (fp = fopen(sedFile, "rt")) != NULL ) 
     { fclose(fp);  FOUND_SNSEDFILE=1 ; }
   else {
     // next check official area
     sprintf(SNPATH, "%s/%s", PATH_SNDATA_ROOT, subdir_snsed );
     sprintf(sedFile,"%s/%s", SNPATH, INPUTS.inFile_snsed );
     if ( (fp = fopen(sedFile, "rt")) != NULL ) 
       { fclose(fp);  FOUND_SNSEDFILE=1 ; }
   }
   
   // - - - - - - SN-SED file is NOT defined  - - - - - - 

   // // Jan 2017 allow missing SN-SED file if no k-corrections are defined
   if ( NKCOR == 0 && FOUND_SNSEDFILE == 0 ) {
     hardWire_snsed_bins();
     return(SUCCESS) ; 
   } 

   if ( FOUND_SNSEDFILE == 0 ) {
     sprintf(c1err,"Cannot find SED file '%s' .", INPUTS.inFile_snsed );
     sprintf(c2err,"Check %s for available SEDs", SNPATH );
     errmsg(SEV_FATAL, 0, fnam,  c1err, c2err); 
   }

   //  RD_SEDFLUX:

   sprintf(sedcomment,"K-cor" );

   // set limits to read in sed file
   DAYrange[0]   = SNSED.TREST_MIN ;
   DAYrange[1]   = SNSED.TREST_MAX ;
   LAMrange[0]   = SNSED.LAMBDA_MIN ;
   LAMrange[1]   = SNSED.LAMBDA_MAX ;

   // read sed file.
   rd_sedFlux(sedFile, sedcomment, DAYrange, LAMrange,
	      MXEP, MXLAM_SN, INPUTS.FLUXERR_FLAG,
	      &NDAY, DAYREAD, &DAYSTEP,
	      &NLAM, LAMREAD, &LAMSTEP,
	      DUMMYFLUX, DUMMYERR, &nflux_nan );

   // check array sizes
   check_snsed_bins();

   // extract arrays and load kcor structures

   dayoff = 0.0;
   
   for ( jflux0 = 0; jflux0 < NDAY*NLAM ; jflux0++ ) {

     iep0  = jflux0 / NLAM;              
     iep1  = iep0+1;           
     ilam0 = jflux0 - NLAM * iep0 ;    
     ilam1 = ilam0 + 1;

     flux    = *(DUMMYFLUX+jflux0);
     day     = *(DAYREAD+iep0);
     lam     = *(LAMREAD+ilam0);

     // check debug option to force power-law SN-vs-lam flux
     if ( POWERLAW > -90.0 ) 
       { flux = 1.0E-10 * pow( (lam/5000.) , POWERLAW ); }
     

     if ( flux == 0.0 ) { flux = FMIN; }  // avoid zero flux

     // check for optional flux-fudge
     flux  *= GET_SNSED_FUDGE(lam);

     flux_converter(lam, flux, &fnu, &fcount ); // returns fnu & fcount

     // allow for t=0 to be explosion day instead of max-luminosity
     if ( jflux0 == 0 && day >= 0.0 ) 
       { dayoff = INPUTS.TREF_EXPLODE ; }

     SNSED.EPOCH[iep1]             = day + dayoff ;
     SNSED.LAMBDA[iep1][ilam1]     = lam ;     // A
     SNSED.FLUX_WAVE[iep1][ilam1]  = flux ;    // erg/s/cm^2/A         
     SNSED.FLUX_NU[iep1][ilam1]    = fnu  ;    // erg/s/cm^2/Hz

   } // end of jflux loop


   // load more structure inf
   SNSED.LAMBDA_MIN  = SNSED.LAMBDA[NDAY][1] ;
   SNSED.LAMBDA_MAX  = SNSED.LAMBDA[NDAY][NLAM] ;
   SNSED.NBIN_LAMBDA = NLAM ;
   SNSED.LAMBDA_BINSIZE = LAMSTEP;

   SNSED.NEPOCH      = NDAY ;
   SNSED.TREST_MIN   = SNSED.EPOCH[1] ;      // Jan 2017
   SNSED.TREST_MAX   = SNSED.EPOCH[NDAY] ;   // Jan 2017

   daymin  = SNSED.EPOCH[1] ;
   daymax  = SNSED.EPOCH[NDAY] ;

   printf("  Finished reading %s \n", INPUTS.inFile_snsed );

   printf("  Found %d epochs from day %4.1f to %4.1f . \n", 
           NDAY, daymin, daymax );

   printf("  Each epoch has %d lambda bins with binsize = %4.2f A . \n",
	   SNSED.NBIN_LAMBDA, SNSED.LAMBDA_BINSIZE );

   printf("  Spectra stored between %6.0f and %6.0f A . \n", 
	  SNSED.LAMBDA_MIN, SNSED.LAMBDA_MAX );

   fflush(stdout);


   // Jan 2010: load FILTER[ifilt].SSUM_SN arrays

   iep = 1;

   for ( ifilt=1;  ifilt <= NFILTDEF ; ifilt++ ) {

     FILTER[ifilt].SSUM_SN = 0.0 ;
     LAMMIN = FILTER[ifilt].LAMBDA_MIN;
     LAMMAX = FILTER[ifilt].LAMBDA_MAX;

     for ( ilam=1; ilam <= SNSED.NBIN_LAMBDA; ilam++ ) {

	lam = SNSED.LAMBDA[iep][ilam] ;
	if ( lam < LAMMIN ) { continue ; }
	if ( lam > LAMMAX ) { continue ; }

	// fetch info for this lambda bin
	magflux_info( 0, ifilt, ilam, iep, 
		      &lam, &flux, &wflux, &wfilt ) ;

	trans  = filter_trans( lam, ifilt, 0);
	 
	 if ( trans > 0.0 ) 
	   { FILTER[ifilt].SSUM_SN  += (trans * wfilt); }
	 
      }   // end of ilam loop

   }  // end of 'ifilt' loop


   return SUCCESS;

} // end of rd_snsed


// *******************************************
void check_snsed_bins(void) {

  // Created Jan 2017: abort if NBIN(LAM,DAY) exceeds bound.

  int NLAM    = SNSED.NBIN_LAMBDA ;
  int NDAY    = SNSED.NEPOCH ;
  char fnam[] = "check_snsed_bins" ;

   // check array size
   if ( NLAM >= MXLAM_SN ) {
     sprintf ( c1err, "NLAM = %d exceeds array bound of MXLAM_SN=%d.", NLAM, MXLAM_SN );
     sprintf ( c2err, "See LAMBDA_RANGE: in input file.");
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
   }
   if ( NDAY >= MXEP ) {
     sprintf ( c1err, "NEPOCH = %d exceeds array bound.", NDAY );
     sprintf ( c2err, "Check SED file.");
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
   }

  return ;
}   // end check_snsed_bins

// *******************************************
void  hardWire_snsed_bins(void) {

  // Created Jan 16 2017 by R.Kessler
  //
  // without SN SED, hard-wire reasonable binning since filterrs 
  // and primary lambda bins are defined by the SN SED bins.
  // This allows NOT defining an SN spectral time series via
  // the SN_SED key.
  // Jan 2023: HARDWIRE_LAMBIN -> 10 (was 20)
  
  int  ilam, iday, NLAM,  NDAY ;
  double lamBin, dayBin, tmpRange, xi  ;
  char fnam[] = "hardWire_snsed_bins" ;

#define HARDWIRE_LAMBIN 10.0 // Angstroms
#define HARDWIRE_DAYBIN  1.0 // days

  // ------------- BEGIN --------------

  printf("\n %s: no SN SED --> hard wire bin sizes: \n", fnam );
  printf("\t binsize(lam) = %.1f \n", HARDWIRE_LAMBIN );
  printf("\t binsize(day) = %.1f \n", HARDWIRE_DAYBIN );
  fflush(stdout);
	 
  // start with lambda bins
  tmpRange       = (SNSED.LAMBDA_MAX-SNSED.LAMBDA_MIN) ;
  lamBin         = HARDWIRE_LAMBIN ;
  NLAM           = (int)(tmpRange/lamBin) + 1 ;
  SNSED.LAMBDA_BINSIZE = lamBin ;  
  SNSED.NBIN_LAMBDA    = NLAM ;

  // now day bins
  tmpRange  = (SNSED.TREST_MAX - SNSED.TREST_MIN) ;
  dayBin    = HARDWIRE_DAYBIN ;
  NDAY       = (int)(tmpRange/dayBin) + 1 ;
  SNSED.NEPOCH      = NDAY ;

  check_snsed_bins();

  // load each lambda and day bin value
  for ( iday=1; iday <= NDAY ; iday++ ) {      

    xi = (double)(iday-1);
    SNSED.EPOCH[iday] = SNSED.TREST_MIN + dayBin*xi ;

    for ( ilam=1; ilam <= NLAM; ilam++ ) {
      xi = (double)(ilam-1);
      SNSED.LAMBDA[iday][ilam] = SNSED.LAMBDA_MIN + lamBin*xi ;      
    }
  }
  
  return ;

} // end hardWire_snsed_bins

// *******************************************
int index_primary( char *name) {

  // Dec 8, 2011
  // Search INPUTS.name_primary[iprim] and return iprim.

  char fnam[] = "index_primary";
  int  iprim ;
  char *cptr;

  // -------------BEGIN ------------
  
  for ( iprim = 1; iprim <= INPUTS.NPRIMARY ; iprim++ ) {
    cptr = INPUTS.name_PRIMARY[iprim];
    if ( strcmp(name,cptr) == 0 ) { return(iprim) ;  }
  }


  // if system is AB and not yet defined, then define it here.
  // This is the only mag-system that is internally defined;
  // Users must specify the others using the PRIMARY_SED key
  // (for Vega, BD17, etc ...)
  if ( strcmp(name,"AB") == 0 ) {
    INPUTS.NPRIMARY++ ;
    iprim = INPUTS.NPRIMARY ;
    sprintf(INPUTS.name_PRIMARY[iprim], "%s",   name );

    // if user does not specify AB_SED in the input file,
    // then hard wire flatnu.dat
    char *inFile = INPUTS.inFile_PRIMARY[iprim] ;
    if ( strlen(inFile) == 0 ) { sprintf(inFile, "flatnu.dat" ); }
    return(iprim) ;
  }

  // if we get here then abort on error
  sprintf(c1err,"Could not find primary for '%s'", name);
  sprintf(c2err,"Check PRIMARY_SED key(s) in kcor-input file.");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 

  return(-9);
  
} // end of index_primary

// *************************************************
int rd_primary ( int INDX, char *subdir ) {

  /****************
    Read primary SED from file with format

     lambda(A)  FLux(erg/s/cm^2/A)

    Note that lambda binning can be non-uniform.
    The raw  SED is read into a "_RAW" arrays,
    which is then re-binned to have same bins as SN 
    for global structure.

   Nov 10 2021: replace FILTER_LAMBDA_MAX -> STORE_LAMBDA_MAX
                to store primary SED out to max of filter or spectrograph.

   Nov 8 2022: skip comment lines in SED-primary file


  ****************/

   FILE  *fp;
   double lambda, flam, fnu, fcount, LAMMIN_READ, LAMMAX_READ  ;
   int  ilam, NBIN, ibin, gzipFlag;
   char SNPATH[MXCHAR_FILENAME], fullName[MXCHAR_FILENAME];
   char *refName, *sedFile, line[200];
   char fnam[] = "rd_primary" ;

   /* --------------------- BEGIN -------------------------- */

   refName = PRIMARYSED[INDX].MAGSYSTEM_NAME;
   sedFile = PRIMARYSED[INDX].MAGSYSTEM_SEDFILE ;

   // idiot checks
   if ( strlen(refName) == 0 ) {
     sprintf(c1err,"Undefined MAGSYSTEM_NAME for INDX=%d", INDX);
     sprintf(c2err,"Check PRIMARYSED[INDX].MAGSYSTEM_NAME");
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
   }
   if ( strlen(sedFile) == 0 ) {
     sprintf(c1err,"Undefined MAGSYSTEM_SEDFILE for INDX=%d", INDX);
     sprintf(c2err,"Check PRIMARYSED[INDX].MAGSYSTEM_SEDFILE");
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
   }

   printf("\n\n  ***** READ PRIMARY %s SED (INDX=%d) ***** \n", 
	  refName, INDX );
 
   sprintf(SNPATH, "%s/%s", PATH_SNDATA_ROOT, subdir );

   int OPENMASK = OPENMASK_IGNORE_DOCANA ;
   fp = snana_openTextFile (OPENMASK,SNPATH, sedFile, fullName, &gzipFlag );

   if ( fp == NULL ) {
     sprintf(c1err,"%s", "Could not open file");
     sprintf(c2err,"%s", sedFile);
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
   }
   else {
     printf("\t Opened : %s \n", fullName );
   }

   // init a few things.

   PRIMARYSED[INDX].NBIN_LAMBDA_RAW = 0;
   ilam     = 0;
   LAMMIN_READ = 999999;
   LAMMAX_READ = 0;
   int NRAW=0;

   while( fgets(line, 200, fp)  != NULL ) {

     if ( commentchar(line) ) { continue; }
     sscanf(line, "%le %le", &lambda, &flam );


     if ( lambda < STORE_LAMBDA_MAX + 20. ) {

       ilam++;  NRAW=ilam;
  
       // keep track of min and max lambda read from file

       if ( lambda < LAMMIN_READ ) 	LAMMIN_READ = lambda;
       if ( lambda > LAMMAX_READ ) 	LAMMAX_READ = lambda;

       if ( ilam < MXLAM_PRIMARY ) {
	 PRIMARYSED[INDX].LAMBDA_RAW[ilam]     = lambda;
	 PRIMARYSED[INDX].FLUX_WAVE_RAW[ilam]  = flam;    // erg/s/cm^2/A
       }
     }
      
   }   /* end while */

   if ( NRAW >= MXLAM_PRIMARY ) { 
     sprintf(c1err,"NBIN(%s)=%d exceeds bound of MXPRIMARY=%d",
	     refName, NRAW, MXLAM_PRIMARY );
     sprintf(c2err,"check %s", sedFile);
     errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
   }

   PRIMARYSED[INDX].NBIN_LAMBDA_RAW = NRAW ;
   PRIMARYSED[INDX].LAMBDA_MIN_RAW  = LAMMIN_READ ;
   PRIMARYSED[INDX].LAMBDA_MAX_RAW  = LAMMAX_READ ;

   printf("   Read %d lambda bins from %d to %d A  \n"
	  , PRIMARYSED[INDX].NBIN_LAMBDA_RAW 
	  , (int)LAMMIN_READ
	  , (int)LAMMAX_READ
	  );

   // --------------------------
   // re-bin primary SED to uniform lambda grid

   rebin_primary ( 
		   PRIMARYSED[INDX].NBIN_LAMBDA_RAW    // (I) # raw lambda bins
		 ,&PRIMARYSED[INDX].LAMBDA_RAW[1]     // (I) raw lambda array 
		 ,&PRIMARYSED[INDX].FLUX_WAVE_RAW[1]   // (I) raw flux array
		 ,&NBIN                            // (O) # lambda bins
		 ,&PRIMARYSED[INDX].LAMBDA[1]       // (O) lam array
		 ,&PRIMARYSED[INDX].FLUX_WAVE[1]    // (O) re-binned flux
		       );

   PRIMARYSED[INDX].NBIN_LAMBDA = NBIN;
   PRIMARYSED[INDX].LAMBDA_MIN  = PRIMARYSED[INDX].LAMBDA[1] ;
   PRIMARYSED[INDX].LAMBDA_MAX  = PRIMARYSED[INDX].LAMBDA[NBIN] ;
   PRIMARYSED[INDX].LAMBDA_BINSIZE = 
     PRIMARYSED[INDX].LAMBDA[2] - PRIMARYSED[INDX].LAMBDA[1] ;


   ibin = (int)PRIMARYSED[INDX].LAMBDA_BINSIZE;
   printf("   Re-bin %s to uniform (%d A) lambda bins from %d to %d A\n", 
	  refName, ibin ,
	  (int)PRIMARYSED[INDX].LAMBDA_MIN,
	  (int)PRIMARYSED[INDX].LAMBDA_MAX
	  );

   //  Convert energy flux to count flux in each lambda bin.
   for ( ilam=1; ilam<=PRIMARYSED[INDX].NBIN_LAMBDA; ilam++ ) {
     lambda     = PRIMARYSED[INDX].LAMBDA[ilam];
     flam       = PRIMARYSED[INDX].FLUX_WAVE[ilam];
     flux_converter(lambda, flam, &fnu, &fcount ); // returns fnu & fcount
     PRIMARYSED[INDX].FLUX_NU[ilam]    = fnu;

     if ( ilam==1 ) {
       printf("\t Flam(%s,LAM=%.1f) = %le \n", refName, lambda, flam);
       fflush(stdout);
     }

   }


   fflush(stdout);

   fclose ( fp );
   return SUCCESS;

}  // end of rd_primary


// ***********************************
double GET_SNSED_FUDGE(double lam) {

  // Jan 18, 2009 R.Kessler
  // Return SN flux-fudge for this lambda.
  // If lambda = 0, read fudge-file and store.
  //
  // Jun 4,2010: if lambda is outside defined fudge-range,
  //              return 1.0 instead of aborting.

  int N, ilam_near, ilam,  i ;
  double lamtmp, fudge, maxfudge, lamatmax, a_lam[4], a_fud[4] ;
  FILE *fp;
  char fnam[] = "GET_SNSED_FUDGE";

  // ------------ BEGIN -----------

  if ( lam == 0.0 ) {

    SNSED_FUDGE.NBIN_LAMBDA = 0;
    if ( (fp = fopen(INPUTS.inFile_snsed_fudge, "rt")) == NULL ) 
      { fudge = 1.0 ; return fudge ; }

    printf("\n\t Read SNSED flux-Fudge from : %s \n", 
	   INPUTS.inFile_snsed_fudge );

    maxfudge = 0.0; lamatmax = 0.0;
    while( (fscanf(fp, "%le %le", &lamtmp, &fudge )) != EOF) {
      SNSED_FUDGE.NBIN_LAMBDA++ ;
      N = SNSED_FUDGE.NBIN_LAMBDA ;
      SNSED_FUDGE.LAMBDA[N] = lamtmp ;
      SNSED_FUDGE.SCALE[N]  = fudge ;

      if ( N == 1 )  SNSED_FUDGE.LAMBDA_MIN = lamtmp ;
      SNSED_FUDGE.LAMBDA_MAX = lamtmp ;

      if ( fudge > maxfudge ) {
	maxfudge = fudge ; lamatmax = lamtmp ;
      }

    } // end of while

    printf("\t MAX-FUDGE=%6.2f at LAM=%6.0f A \n", maxfudge, lamatmax );

    fclose(fp);

    return fudge ;

  } // end of lambda=0


  // -------------

  fudge = 1.0 ;
  if ( SNSED_FUDGE.NBIN_LAMBDA == 0 ) return fudge ;

  // June 2010: do nothing if we are outside the fudge range
  if ( lam < SNSED_FUDGE.LAMBDA_MIN ) return fudge ;
  if ( lam > SNSED_FUDGE.LAMBDA_MAX ) return fudge ;


  // find nearest lambda bin
  ilam_near = -9 ;
  N = SNSED_FUDGE.NBIN_LAMBDA ;
  for ( ilam=1; ilam <= N; ilam++ ) {
    lamtmp = SNSED_FUDGE.LAMBDA[ilam];
    if ( lamtmp < lam ) { ilam_near = ilam ; }
  }

  if  ( ilam_near < 0 ) {
    sprintf ( c1err, "ilam_near=%d for lam=%.2f", ilam_near, lam);
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING ); 
  }

  // printf(" xxxx lam=%f,  lamtmp=%f at ilam=%d \n", lam, lamtmp, ilam_near );

  // define lambda and fudge for three nearest bins
  for ( i=1; i <= 3; i++ ) {
     ilam = ilam_near + i - 1 ; 
     a_lam[i] = SNSED_FUDGE.LAMBDA[ilam] ;   
     a_fud[i] = SNSED_FUDGE.SCALE[ilam] ;   
   }

   // interpolate to get flux-fudge  at wavelen = "lam" 

   fudge = interp_1DFUN ( OPT_INTERP_FILTER, lam, NLAMBIN_INTERP, 
			  &a_lam[1], &a_fud[1], "fudge" );

   return fudge ;

}   // end of GET_SNSED_FUDGE


// ********************************************************
void flux_converter ( 
		      double  lambda       // (I) lambda (A)
		     ,double  flux_wave    // (I) dE/dlam 
		     ,double *flux_nu      // (O) dE/dnu
		     ,double *flux_count   // (O) dN/dlam
		        ) {

  /****
   convert input flux_wave (erg/s/cm^2/A) to
        - flux_nu     (erg/s/cm^2/Hz)
        - flux_count  (/s/cm^2/A)

   Note that "hc" has units of "erg * A", 
   and LIGHT_A is the speed of light in A/s.
  
  *****/

  double ephot, sqlam;

  // --------------------- BEGIN -------------------

  // ephot = h nu = hc/lambda

  ephot       = hc/lambda;
  sqlam       = lambda * lambda;

  *flux_count = flux_wave / ephot;
  *flux_nu    = flux_wave * sqlam / LIGHT_A;
  
}  // end of flux_converter




// ****************************************
void rebin_primary ( int  nblam_in,  double *lam_in,  double *flux_in, 
		     int *nblam_out, double *lam_out, double *flux_out ) {

  /***

    Jan 2010:
    re-bin primary SED to use the SNSED lambda binsize and the
    min/max lambda range for the filters.  This is needed
    because the primary (input) SED does NOT have uniform lambda bins;
    the output of this function has uniform lambda bins.
    Use linear interpolation for lam_in values closest
    to SNSED.LAMBDA grid point.

    Nov 3 2021: replace FILTER_LAMBDA_MIN[MAX] with STORE_LAMBDA_MIN[MAX]
             (to allow for spectrograph with broader wave range)

   May 31 2024: replace hard wired DLAM=10 with DLAM = SNSED.LAMBDA_BINSIZE

  ****/

  int NBLAM, ilam, idump=0;
  double  DLAM, LAM, LAM0,  LAM1, F0, F1, slope, FLUX_OUT  ;
  char fnam[] = "rebin_primary" ;

  /* ---------------------------- BEGIN -------------------- */

  DLAM = SNSED.LAMBDA_BINSIZE;
  
  LAM  = LAM0 = LAM1 = F0 = F1 = 0.0 ; 
  NBLAM = 0 ;

  while ( LAM < STORE_LAMBDA_MAX ) {
    LAM += DLAM ;
    if ( LAM < STORE_LAMBDA_MIN ) continue ;

    NBLAM++ ;

    // find input lambda(s) that bracket this LAM
    ilam = 0;
    while ( *(lam_in + ilam) < LAM+DLAM ) {
      LAM0 = *(lam_in + ilam) ;
      LAM1 = *(lam_in + ilam + 1 ) ;
      F0   = *(flux_in + ilam) ;
      F1   = *(flux_in + ilam + 1 ) ;

      if ( LAM >= LAM0 && LAM <= LAM1 ) { goto INTERPFLUX; }
      ilam++ ;
    }

  INTERPFLUX:

    slope = (F1 - F0)/(LAM1 - LAM0);
    FLUX_OUT = F0 + slope * (LAM  - LAM0);

    *(lam_out+NBLAM-1)  = LAM ;
    *(flux_out+NBLAM-1) = FLUX_OUT ;

    if ( idump ) {
      printf(" xxx LAM=%7.2f LAM0=%7.2f LAM1=%7.2f \n", LAM, LAM0, LAM1 ); 
      printf(" xxx FOUT=%le F0=%le F1=%le \n", FLUX_OUT, F0, F1 ); 
      debugexit("rebin");
    }

  } // LAM

  *nblam_out = NBLAM;


}  // end of rebin_primary


// ****************************************************
int kcor_grid (void) {

  // Shell to call kcor_eval  at each grid point.
  // Feb 15, 2006: add loop over AV index
  // Mar  7, 2006: add magobs8 output arg to kcor_eval(); 
  //                fill SNSED.SDSSMAG_OBS
  //
  // Jun 8, 2009; check kcor_eval output for nan
  //
  // Jun 9, 2009: 
  //  * change all float -> double
  //  * if kcor is outside valid range (KCORMIN - KCORMAX), 
  //    then set it to really crazy NULLVAL so that sim & fitter 
  //    know to ignore it
  // 
  // Nov 12, 2010: loop over NKCOR+KCOR_EXTRA to get synthetic
  //               'magobs' for the rest-frame filters that are
  //               needed by snana.
  // -------------------------------------------------

   char ctmp[20]
     ,  fnam[] = "kcor_grid"
     ;

   int  ikcor, OPT, ifilt_rest, ifilt_obs ;
   int  i_epoch, i_z, i_av, i_ebv, FLAG_MAGOBS, NZBIN ;
   
   double 
     z, epoch, av, dum, kcor, kcormin, kcormax
     ,err, ovp, magobs[MXMWEBV+2], magtmp, dxt, debv
     ; 

   /* -------------------- BEGIN ------------------ */

   OPT = 0;

   printf("\n  ***** START LOOPING for KCOR GRID ***** \n" );


   for ( ikcor=1; ikcor <= NKCOR + NKCOR_EXTRA ; ikcor++ ) {
       
     // check which filter indices

      ifilt_rest  = -1;
      ifilt_obs   = -1;

      index_filter ( ikcor, &ifilt_rest, &ifilt_obs );

      if ( ikcor <= NKCOR ) {
	NZBIN  = INPUTS.NBIN_REDSHIFT;
	ctmp[0]=0;
      }
      else {
	NZBIN = 1;
	sprintf(ctmp, "%s", "EXTRA");
      }

      printf("  Compute %s %s for '%s' (rest) => '%s' (obs) \n",
	     ctmp, KCORSYM[ikcor], KCORLIST[ikcor][0], KCORLIST[ikcor][1] );
      kcormin = 999999. ;
      kcormax = -99999. ;

      printf("\t AV = ");


      // now loop over redshift and epoch 

      for ( i_av=1;  i_av<=INPUTS.NBIN_AV;   i_av++ ) {

	  dum    = (double)(i_av-1) ;
	  av     = INPUTS.AV_MIN + dum * INPUTS.AV_BINSIZE;

          printf("%4.2f ", av);
	  fflush(stdout);

	for ( i_z=1;   i_z <= NZBIN ; i_z++ ) {

	  dum   = (double)(i_z-1) ;
	  z     = INPUTS.REDSHIFT_MIN + dum * INPUTS.REDSHIFT_BINSIZE;

	  for ( i_epoch=1; i_epoch<= SNSED.NEPOCH; i_epoch++ ) {

	    epoch = SNSED.EPOCH[i_epoch];  

	    R4KCOR_GRID.REDSHIFT[ikcor][i_av][i_z][i_epoch]  = (float)z ;
	    R4KCOR_GRID.EPOCH[ikcor][i_av][i_z][i_epoch]     = (float)epoch ;

	    // check if these obs mags have already been computed
	    if ( SNSED.R4MAG_OBS[0][ifilt_obs][i_av][i_z][i_epoch] == NULLVAL )
	      { FLAG_MAGOBS = 1 ; }
	    else
	      { FLAG_MAGOBS = 0; }

            kcor_eval( OPT
		       ,av, z, epoch
		       ,ifilt_rest, ifilt_obs 
		       ,FLAG_MAGOBS
		       ,&kcor, &err, &ovp, magobs        // return values
		       );

	    if ( kcor > kcormax ) { kcormax = kcor ; }
	    if ( kcor < kcormin ) { kcormin = kcor ; }

	    // if kcor is outside valid range, then set it to really
	    // crazy NULLVAL so that sim & fitter know to ignore it
	    if ( kcor > KCORMAX_VALID ) { kcor = NULLVAL ; }
	    if ( kcor < KCORMIN_VALID ) { kcor = NULLVAL ; }

	    // 6/08/2009: check for nan 
	    if ( isnan(kcor) ) {
	      sprintf(c1err,"kcor=%f  for z=%6.3f T=%6.3f  av=%6.3f",
		      kcor, z, epoch, av);
	      sprintf(c2err,"ifilt_[rest,obs]=%d,%d (%s,%s) FLAG_MAGOBS=%d"
		      ,ifilt_rest, ifilt_obs
		      ,FILTER[ifilt_rest].name
		      ,FILTER[ifilt_obs].name
		      ,FLAG_MAGOBS);
	      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
	    }

	    R4KCOR_GRID.VALUE[ikcor][i_av][i_z][i_epoch] = (float)kcor ;

	    // Feb 2007: store observer mags with array of MW E(B-V)
	    if ( FLAG_MAGOBS > 0 ) {
	      for ( i_ebv = 0; i_ebv <= MXMWEBV; i_ebv++ ) {
		magtmp = magobs[i_ebv];
		if ( isnan(magtmp) ) {
		  sprintf(c1err,"magobs=%f for i_ebv=%d z=%6.3f T=%6.2f",
			  magtmp, i_ebv, z, epoch );
		  sprintf(c2err,"ifilt_[rest,obs]=%d,%d", 
			  ifilt_rest, ifilt_obs);
		  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  

		}

		SNSED.R4MAG_OBS[i_ebv][ifilt_obs][i_av][i_z][i_epoch] = 
		  (float)magtmp;
	      }
	      // store d(mag)/d(xtmw) based on first two bins
	      dxt   = *(magobs + 1) - *(magobs + 0) ;
	      debv = MWEBV_LIST[1] -  MWEBV_LIST[0]  ;
	      SNSED.MW_dXT_dEBV[ifilt_obs][i_av][i_z][i_epoch] = 
		(dxt/debv);
	    }
	    

	  } // end of i_epoch loop 
	}  // end of i_z loop 
      }   // end if i_av loop

      printf(" \n\t %s min/max = %6.3f/%6.3f \n", 
	     KCORSYM[ikcor], kcormin, kcormax);

   }     // end of ikcor loop 
     

   return SUCCESS;


} // end of kcor_grid



// *************************************************
void kcor_eval(int opt                // (I) K cor option ("E" or "N")
               ,double av             // (I) host AV = RV * E(B-V)
	       ,double redshift       // (I) redshift
	       ,double epoch          // (I) SN epoch, days
	       ,int ifilt_rest        // (I) rest filter index
	       ,int ifilt_obs         // (I) observer filter index
	       ,int FLAG_MAGOBS       // (I) non-zer => compute *mag_obs
	       ,double *kcor_value   // (O) K correction value
	       ,double *kcor_error   // (O) error on above
	       ,double *overlap      // (O) rest-observer flux overlap
               ,double *mag_obs       // (O) observer-flux in ifilt_obs
   ) {

/***

  Jan 2005
  Here we compute the K correction for 'redshift' and 'epoch'.
  Integrate over SN template spectrum.

  For now, only kcor_value is returned; 
  kcor8_error is just a place-holder for future work.

   Feb 15, 2006: add "av" as argument and include dust extinction
                 as part of KCOR correction

  Mar 7, 2006: add *mag_obs output arg
  
  Feb 14, 2007: *mag_obs is now returned as an array of dimension
                 MXMWEBV with obserber-mag at each MW E(B-V)

                 Add FLAG_MAGOBS to avoid redundant calc of the same
                 observer-frame mags that get repeated for different
                 K-corrections

  May 15, 2009: when flux < 0, mag=9 instead of 30 to avoid
                abort in simulation & fitter. Found negative
                fluxes for badly interpolated non1a-003 SED.


  Jun 9, 2009: all floats -> double

 ***/

  int   
    ilam_sn
    ,NBIN
    ,iepoch
    ,iebv
    ;

   double 
     LAM
     , LAMMIN_FILT
     , LAMMAX_FILT
     , lam              // lambda (A)
     , flux_sn_rest   // SN flux, rest
     , flux_sn_obs    // SN flux in redshifted frame
     , trans_rest     // filter transmission, rest frame filter 
     , trans_obs      // filter trans in redshifted filter 
     , trans_min
     , conv_sn_rest   // rest frame convolution
     , conv_sn_obs    // redshifted convolution 
     , conv_sn_ovp
     , filtsum_rest
     , filtsum_obs
     , zp_rest
     , zp_obs
     , zp_dif
     , fluxsum_ratio
     , filtsum_ratio
     , oneplusz        // 1+z
     , tmp, arg
     , flux
     , flux_obs[MXMWEBV+1]
     , ftmp, fcount
     , wflux, wfilt
     , mwav, mwxt
     , kcortmp
     , LAMZ
     , RV
     , zero = 0.0
     , ten = 10.0
     ;

   char fnam[] = "kcor_eval";

   /* ------------------ BEGIN ------------------- */

   // init output 

   *kcor_value = NULLVAL ;  
   *kcor_error = NULLVAL ;
   *overlap    = NULLVAL ;

   for ( iebv=0; iebv <= MXMWEBV; iebv++ ) {
     mag_obs[iebv]        = NULLVAL ;
     flux_obs[iebv]       = 0.0 ;
   }

   if ( INPUTS.FASTDEBUG ) { return ; }

   RV = INPUTS.RV_MWCOLORLAW ;
   oneplusz   = ( 1.0 + redshift ) ;

   // get integer epoch index from "epoch" in days.
   iepoch  =  index_epoch ( epoch ) ;  

   // retrieve integrals of Filter-response

   filtsum_rest = FILTER[ifilt_rest].SSUM_SN ;
   filtsum_obs  = FILTER[ifilt_obs].SSUM_SN ;

   // retrieve zero points
   zp_rest = FILTER[ifilt_rest].MAGFILTER_ZP +
             FILTER[ifilt_rest].MAGSYSTEM_OFFSET ; 

   zp_obs  = FILTER[ifilt_obs].MAGFILTER_ZP +
             FILTER[ifilt_obs].MAGSYSTEM_OFFSET ; 

   NBIN    =  SNSED.NBIN_LAMBDA; // number of lambda bins

   /************************************************
         compute SNSED * FILTER(rest) * LAMBDA
   ************************************************/ 

   LAMMIN_FILT = FILTER[ifilt_rest].LAMBDA_MIN ;
   LAMMAX_FILT = FILTER[ifilt_rest].LAMBDA_MAX ;
   conv_sn_rest   = 0.0 ;
   conv_sn_ovp    = 0.0 ;

   for ( ilam_sn=1; ilam_sn <= NBIN; ilam_sn++ ) {

     LAM    = SNSED.LAMBDA[1][ilam_sn]; // get lambda from epoch=1

     if ( LAM >= LAMMIN_FILT && LAM <= LAMMAX_FILT ) {

       trans_rest  = filter_trans( LAM, ifilt_rest, 0 );
     
	if ( trans_rest > 0.0 ) {
	  flux_sn_rest  = snflux ( epoch, LAM, zero, av );   // flux at z=0 
	  conv_sn_rest += flux_sn_rest * trans_rest * LAM ;

	  // June 6, 2008 compute overlap function
	  LAMZ       = LAM * oneplusz ;
          trans_obs  = filter_trans( LAMZ, ifilt_obs, 0 );

	  if  ( trans_obs < trans_rest)  
	    trans_min = trans_obs;
	  else
	    trans_min = trans_rest;

	  conv_sn_ovp += flux_sn_rest * trans_min * LAM;

	}  // end positive trans if-block

     }    // LAMMIN,MAX if-check
   }   // end of ilam loop



   /*********************************************************

       compute conv_obs = 
          SED * FILTER(redshifted) * EXTINCTION 

   **********************************************************/ 

   LAMMIN_FILT  =  FILTER[ifilt_obs].LAMBDA_MIN ;
   LAMMAX_FILT  =  FILTER[ifilt_obs].LAMBDA_MAX ;
   conv_sn_obs = 0.0 ;


   for ( ilam_sn=1; ilam_sn <= NBIN; ilam_sn++ ) {

     LAM        = SNSED.LAMBDA[1][ilam_sn];

     if ( LAM >= LAMMIN_FILT && LAM <= LAMMAX_FILT ) {

       trans_obs   = filter_trans( LAM, ifilt_obs, 0 ); // filter trans

       if ( trans_obs > 0.0 ) {

	 // get redshifted flux needed for K-cor
	 flux_sn_obs  = snflux ( epoch, LAM, redshift, av ); 
	 conv_sn_obs += flux_sn_obs * trans_obs * LAM ;

	 if ( flux_sn_obs == NULLVAL ) { return ; }

	 // get info needed for observed mag
	 
	 if ( FLAG_MAGOBS == 0 ) { continue ; }

	 // get obsserved "flux"
	 flux_converter( LAM, flux_sn_obs, &flux, &fcount ); 

	 // get integration weight "wflux" for this filter
	 magflux_info( 1, ifilt_obs, ilam_sn, iepoch, 
		       &lam, &ftmp, &wflux, &wfilt ) ;

	 for ( iebv=0; iebv <= MXMWEBV; iebv++ ) {
	   mwav = INPUTS.RV_MWCOLORLAW * MWEBV_LIST[iebv] ;
	   tmp  = 0.4 * GALextinct ( RV, mwav, lam,
				     INPUTS.OPT_MWCOLORLAW, INPUTS.PARLIST_MWCOLORLAW, fnam );
	   mwxt = 1./pow(ten,tmp) ;
	   flux_obs[iebv]  += mwxt * wflux * flux * trans_obs  ; 
	   flux_obs[iebv]  += 0.1E-8;
	 }

       }

     }  // end of LAM if-block

   } // end of ilam_sn loop 


/**************************************************

       Final computation

 **************************************************/


   if ( conv_sn_rest <= zero || conv_sn_obs <= zero ) {
     *kcor_value = zero ;
   } 
   else {
      fluxsum_ratio = (conv_sn_rest / conv_sn_obs) ;
      filtsum_ratio = (filtsum_rest / filtsum_obs ) ; 
      tmp           = oneplusz * fluxsum_ratio / filtsum_ratio ;
      zp_dif        = zp_obs - zp_rest ;

      if ( tmp <= zero ) 
	{ kcortmp = zero  ; }
      else
      	{ kcortmp = +2.5 * log10( tmp ) + zp_dif ; }


      *kcor_value = kcortmp ;

      if ( isnan(kcortmp) ) {
	printf("\n\n ================================================== \n");
	printf(" ==> kcor isnan  at T=%6.2f ifilt_rest=%d(%s)\n", 
	       epoch, ifilt_rest, FILTER[ifilt_rest].name );
	printf(" xxxx kcortmp = %le  tmp = %le \n", kcortmp, tmp );
	printf(" xxxx conv_sn_rest = %le   conv_sn_obs=%le \n",
	       conv_sn_rest, conv_sn_obs );
	printf(" xxxx filtsum_rest = %le   filtsum_obs=%le \n",
	       filtsum_rest, filtsum_obs );
	printf(" xxxx zp_dif = %7.4f(obs) - %7.4f(rest) = %7.4f \n",
	       zp_obs, zp_rest, zp_dif );
      }


   }


   // convert observer flux to observer magnitude

   for ( iebv=0 ; iebv <= MXMWEBV; iebv++ ) {
     if ( flux_obs[iebv] <= zero ||  filtsum_obs <= zero) {
       mag_obs[iebv] = 9.0 ;
     }
     else {
       tmp = flux_obs[iebv] / filtsum_obs ;
       mag_obs[iebv] = -2.5 * log10 ( tmp ) + zp_obs ;
     }
   }


   // compute overlap 
   arg = conv_sn_rest * conv_sn_obs ;
   tmp = sqrt ( fabs(arg) );
   if ( tmp > 0.0 ) { *overlap = conv_sn_ovp / tmp ; }

   return;

}  // end of kcor_eval function




double filter_trans( double lam, int ifilt, int idump ) {

/**************************
  Returns interpolated filter transmission at
  wavelength = "lam".  "ifilt" specifies which filter.

  Nominal interpolation includes parabolic fit to 
  three points ... be careful at boundary.

  Sep 12 2013: fix to work with non-uniform filter binning.

 **************************/

  int NBIN; 
  double trans, *ptrLam, *ptrTrans ;

   /* --------------------- BEGIN ------------------- */


   // return 0 if lam is outside defined filter range 
   trans = 0.0 ;
   if ( lam < FILTER[ifilt].LAMBDA_MIN ) { return trans ; }
   if ( lam > FILTER[ifilt].LAMBDA_MAX ) { return trans ; }

   // get ilam_filt = lambda index for filter that has 
   //  closest lambda to SN "LAM" 

   NBIN      = FILTER[ifilt].NBIN_LAMBDA ; 
   ptrLam    = &FILTER[ifilt].LAMBDA[1] ;   
   ptrTrans  = &FILTER[ifilt].TRANS[1] ;   

   trans = interp_1DFUN ( OPT_INTERP_FILTER, lam, NBIN, 
			  ptrLam, ptrTrans, "filter-trans" );
      
   return trans;

} // end of filter_trans



// ***********************************************************************
double snflux( double epoch, double lambda, double redshift, double av ) {

  /*******

  Returns SN flux dE/dlam for this "epoch" and at this 
  "lambda/(1+redshift)".
  
  Use interpolation between SN bins of lambda.

  Feb 7, 2007: apply MW extinction of redshift is non-zero;
               valid for mags, but makes K-cors invalid.

  Jun 8, 2009: don't let f2=0 (see FMIN) to avoid /0 for non1a that
               have flux=0 in some regions.
 
               Change input args from float to double.

 ********/

   double     
     flux         // output value 
     , LMIN        // min lambda stored for SN spectrum (A) 
     , LMAX        // max ... 
     , LAMZ        // lambda / 1+z  = rest-frame lambda
     , binsize     // bin size for SN spectra 
     , xlam
     , a_lam[4]    // lambda array for interpolation 
     , a_flux[4]   // flux array for interp 
     , a_fnorm[4]  // normalized flux to avoid number problems 
     , dum
     , tmp  
     , warp
     , RV
     , renorm
     ;

   int iepoch, ilambda, ilam1, ilam, i;
   int idebug = 0;

   char fnam[] = "snflux";

   /* --------------------- BEGIN -------------------- */

   flux    = 0.0 ;             // init output 
   LMIN    = SNSED.LAMBDA_MIN;
   LMAX    = SNSED.LAMBDA_MAX;
   binsize = SNSED.LAMBDA_BINSIZE ;

   if ( lambda < LMIN || lambda > LMAX ) {
      return flux;
   }

   //   if ( av == 2. && epoch == -1888. && lambda == 3000.  )   idebug = 1;   

   // find "ilambda" bin for this "lambda"  
   LAMZ    = lambda / ( 1.0 + redshift ) ;
   xlam    = (LAMZ - LMIN) / binsize + 0.5 ;
   ilambda = (int)xlam + 1;

   iepoch  = index_epoch ( epoch );


  // find "ilam1" = 1st of three lambda bins to use for interpolation 


   if ( ilambda <= 1 )
     { ilam1 = 1 ; }      // watch lower boundary 
   else if ( ilambda >= SNSED.NBIN_LAMBDA ) 
     { ilam1 = SNSED.NBIN_LAMBDA - 2 ; }
   else
     { ilam1 = ilambda - 1 ; }


   // generate flux vs. lambda triplet to use for interpolation 

   for ( i=1;  i <= 3;  i++ ) {
      ilam = ilam1 + i - 1;
      a_lam[i]  = SNSED.LAMBDA[iepoch][ilam];
      a_flux[i] = SNSED.FLUX_WAVE[iepoch][ilam];
   }

   /* re-normalize to middle flux to avoid potential
      problem with very small numbers 
      6/08/2009: don't let renorm=0 to avoid divide-by-zero
   */

   if ( a_flux[2] == 0.0 ) 
     { renorm = 1.0 ; }
   else 
     { renorm = a_flux[2] ; }

   a_fnorm[1] = a_flux[1]/renorm ;
   a_fnorm[2] = a_flux[2]/renorm ;
   a_fnorm[3] = a_flux[3]/renorm ;

/* idiot check: make sure that a_lam[2] is within 
               "binsize" of requested lambda */

   dum = fabs(a_lam[2] - LAMZ);  // abs -> fabs (Feb 20 2014)
   if ( dum > binsize ) {
     return (double)NULLVAL ;
   }


   // do the interpolation 

   flux = interp_1DFUN ( OPT_INTERP_SNFLUX, LAMZ, NLAMBIN_INTERP, 
			  &a_lam[1], &a_fnorm[1], "SNFLUX" );
   flux *= renorm ;

   if ( idebug == 1  || isnan(flux) ) {
      printf("\n ------ LAMBDA = %6.0f, Z=%4.2f  Trest=%6.2f ----------- \n",
	     lambda, redshift, epoch );

      printf("  LAMBDA/(1+z) = %6.0f  (ilambda=%d) \n", LAMZ, ilambda );
      printf("  LMIN-LMAX = %5.0f - %5.0f  (xlam=%f binsize=%f) \n", 
	     LMIN, LMAX, xlam, binsize);
      printf("  a_lam  = %6.0f %6.0f %6.0f \n",
	     a_lam[1], a_lam[2], a_lam[3] );

      printf("  a_flux = %e  %e  %e \n",
	     a_flux[1], a_flux[2], a_flux[3] );

      printf("  interp SN flux = %e \n", flux );

   }

   // Apply warping parameter using AV and CCM89 law.
   // Note opt=94 => ODonell update for optical/NIR

   if ( INPUTS.AV_OPTION == 2 ) {
     RV          = INPUTS.RV_MWCOLORLAW ;
     tmp         = 0.4 * GALextinct ( RV, av, LAMZ,
				      INPUTS.OPT_MWCOLORLAW, INPUTS.PARLIST_MWCOLORLAW, fnam);
     warp        = 1.0/pow(TEN,tmp) ;
     flux       *= warp ;
   }

   if ( idebug == 1  )  { printf("  Final SN flux = %e \n", flux ); }

   return flux;

}  // end of snflux



double primaryflux( int iprim, double lambda ) {

/********************************

  Returns primary flux for this primary 'iprim' and wavelength LAMBDA.
  Use interpolation between  bins of lambda.
  Note that flux unit is  erg/s/cm**2/A .  

  Feb 20, 2014: abs -> fabs in a few places  (caught with c++)

  May 31 2024: replace hard-wired LBIN=10 with LBIN=SNSED.LAMBDA_BINSIZE

********************************/

   double     
     flux      // output value 
     , LMIN     // min lambda stored for SN spectrum (A) 
     , LMAX     // max ... 
     , LBIN     // bin size
     , xlam
     , a_lam[4]    // lambda array for interpolation 
     , a_flux[4]   // flux array for interp 
     , a_fnorm[4]  // normalized flux to avoid number problems 
     , dum
     , renorm
     ;

   int NBIN, ilambda, ilam1, ilam, i;

   int idebug = 0;

   char fnam[] = "primaryflux";
   
   /* --------------------- BEGIN -------------------- */

   flux    = 0.0 ;             // init output 
   LMIN    = PRIMARYSED[iprim].LAMBDA_MIN;
   LMAX    = PRIMARYSED[iprim].LAMBDA_MAX;
   LBIN    = SNSED.LAMBDA_BINSIZE;
   NBIN    = PRIMARYSED[iprim].NBIN_LAMBDA ;
   
   if ( lambda < LMIN || lambda > LMAX ) {
      return flux;
   }

  /* find "ilambda" bin for this "lambda"  */

   xlam    = ((lambda - LMIN)/LBIN) + 0.5 ;
   ilambda = (int)xlam + 1 ;

   // find "ilam1" = 1st of three lambda bins to use for interpolation 

   if ( ilambda == 1 )
     ilam1 = 1 ;      // watch lower boundary 
   else if ( ilambda == NBIN ) 
      ilam1 = NBIN - 2 ;
   else
      ilam1 = ilambda - 1 ;


   // generate flux vs. lambda triplet to use for interpolation

   for ( i=1;  i <= 3;  i++ ) {
      ilam = ilam1 + i - 1;
      a_lam[i]  = PRIMARYSED[iprim].LAMBDA[ilam];
      a_flux[i] = PRIMARYSED[iprim].FLUX_WAVE[ilam];
   }

   /* re-normalize to middle flux to avoid potential
      problem with very small numbers */

   if ( a_flux[2] == 0.0 ) 
     renorm = 1.0;
   else 
     renorm = a_flux[2] ;

   a_fnorm[1] = a_flux[1]/renorm ;
   a_fnorm[2] = a_flux[2]/renorm ;
   a_fnorm[3] = a_flux[3]/renorm ;

/* idiot check: make sure that a_lam[2] is within 
               "LBIN" of requested lambda */

   dum = fabs(a_lam[2] - lambda);
   if ( dum > LBIN ) {
      sprintf(c1err,"a_lam = %7.2f %7.2f %7.2f, but lam=%7.2f .",
               a_lam[1], a_lam[2], a_lam[3], lambda );

      sprintf(c2err,"LMIN=%7.2f, LBIN=%5.2f,  ilam=%d .",  
	      LMIN, LBIN, ilambda );

      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
   }

   // do the interpolation 

   //   flux_old = interp8 ( OPT_INTERP_VEGAFLUX, a_lam, a_fnorm, lambda ) ;
   flux = interp_1DFUN ( OPT_INTERP_VEGAFLUX, lambda, NLAMBIN_INTERP, 
			  &a_lam[1], &a_fnorm[1], "VEGAFLUX" );

   flux *= renorm ;

   if ( idebug == 1 ) {
      printf("\n -------- LAMBDA = %6.0f ------------- \n",
	     lambda );

      printf("  LAMBDA = %6.0f  (ilambda=%d) \n", lambda, ilambda );

      printf("  a_lam  = %6.0f %6.0f %6.0f \n",
	     a_lam[1], a_lam[2], a_lam[3] );

      printf("  a_flux = %e  %e  %e \n",
	     a_flux[1], a_flux[2], a_flux[3] );


      printf("  interp SN flux = %e \n", flux );

   }

   return flux;

} // end of primaryflux


// *******************************************
int index_epoch ( double epoch ) {

  //  Convert "epoch" to index in SNSED structure.
  //  For speed, jump every 10 bins at first,
  //  then check every bin.

   int iepoch, iepoch_best, iepoch_quick, I1, I2;
   double diff, diffmin, tmp, binsize;

   int JUMP = 10;

   /* ------------------ BEGIN ------------------- */


   binsize = SNSED.EPOCH[2] - SNSED.EPOCH[1];
   
   tmp = (epoch - SNSED.EPOCH[1]) / binsize;

   iepoch_quick = (int)(tmp+0.5) + 1;

   return iepoch_quick;
   

   // below is a more general (but slower) routine for arbitrary EPOCH binning

   diffmin     = 1.0E8;
   iepoch_best = 999;

   for ( iepoch=1; iepoch<= SNSED.NEPOCH; iepoch+=JUMP ) {
      tmp  = SNSED.EPOCH[iepoch];       
      diff = fabs(tmp - epoch);

      if ( diff < diffmin ) {
	 iepoch_best = iepoch;
	 diffmin = diff;
      }
   }


   I1 = iepoch_best - JUMP;
   I2 = iepoch_best + JUMP;

   for ( iepoch=I1; iepoch<= I2; iepoch++ ) {
      tmp  = SNSED.EPOCH[iepoch];       
      diff = fabs(tmp - epoch);

      if ( diff < diffmin ) {
	 iepoch_best = iepoch;
	 diffmin = diff;
      }
   }

 
   
   return iepoch_best;

}


void index_filter ( int ikcor, int *ifilt_rest, int *ifilt_obs ) {

/****************************
  for input "ikcor", returns "ifilt_rest" and "ifilt_obs"
  indices corresponding to rest-frame and observer frame
  filter for the given "ikcor".

******************************/

   int ifilt;

/* ---------------------- BEGIN ---------------- */

   *ifilt_rest = -1;
   *ifilt_obs  = -1;

   for ( ifilt=1; ifilt <= NFILTDEF; ifilt++ ) {
      if ( strcmp( KCORLIST[ikcor][0],FILTER[ifilt].name )==0 )
	 *ifilt_rest = ifilt ;

      if ( strcmp( KCORLIST[ikcor][1],FILTER[ifilt].name )==0 )
	 *ifilt_obs = ifilt ;
   }

  // now do error checking

   if ( *ifilt_rest < 0 ) {
      sprintf(c1err, "Could not find rest frame filter '%s' ", 
              KCORLIST[ikcor][0] );
      sprintf(c2err, "ikcor = %d", ikcor );
      errmsg(SEV_FATAL, 0, "index_filter", c1err, c2err ); 
   }

   if ( *ifilt_obs < 0 ) {
      sprintf(c1err, "Could not find redshifted filter '%s' ", 
              KCORLIST[ikcor][1] );
      sprintf(c2err, "ikcor = %d", ikcor );
      errmsg(SEV_FATAL, 0, "index_filter", c1err, c2err ); 
   }

} // end of index_filter



// *************************************************
int snmag(void) {

  /***
    fill 
        - SNSED.FLUXSUM[iepoch][ifilt]  (erg/s/cm^2)
        - SNSED.MAG_RST[iepoch][ifilt]

  Mar 11, 2009:  fix bug so that IEPOCH_SNPEAK  is evaluated
                 correctly when Trest has non-integer values.

  Jun 9, 2009: float --> double

  Jan 21, 2010: if filtsum_check (integral over SED-lambda range)
                does not match integral over entire filter range
                (to within 20%),
                then set mag=666 as a flag that it's not reliable.

  Feb 01, 2010: include DM15 in printout 
                 (instead of obsolete CURRENT-CHECK DIFF)

                Compute peak and day15 for each filter instead
                of just using iepoch of last filter.

  ****/


  int ifilt, iepoch, ilam, iepoch_peak, iepoch_15day ;

  double epoch, epoch_peak, epoch_15day, near_peak, near_15day ;
  double arg, lam, trans, flux,  fluxsum_sn, filtsum ;
  double filtsum_check, filter_check, zp, mag, wflux, wfilt ;
  double LAMMIN, LAMMAX, LAMAVG, dm15, mag15 ;
  double zero = 0.0 ;
  double TPEAK  = 0.0;
  double T15DAY = 15.0;
  char fnam[] = "snmag";

  /* ------------------- BEGIN --------------------- */

   if ( SNSED.NBIN_LAMBDA <= 0 ) { return SUCCESS ; }

   printf("\n ***** Compute SN magnitude in each filter/epoch ***** \n" );
   printf("\n");
  

   for ( ifilt=1;  ifilt<=NFILTDEF; ifilt++ ) {

     // get global properties of this filter

     zp      = FILTER[ifilt].MAGFILTER_ZP + 
               FILTER[ifilt].MAGSYSTEM_OFFSET ;

     filtsum = FILTER[ifilt].SSUM_SN ; 
     LAMMIN  = FILTER[ifilt].LAMBDA_MIN ;
     LAMMAX  = FILTER[ifilt].LAMBDA_MAX ;
     LAMAVG  = FILTER[ifilt].LAMAVG;   

     near_peak   = 9999. ; // Trest nearest peak
     near_15day  = 9999. ; //

     for ( iepoch=1; iepoch<=SNSED.NEPOCH; iepoch++ ) {

       fluxsum_sn     = 0.0;
       filtsum_check  = 0.0 ;

       for ( ilam=1; ilam  <= SNSED.NBIN_LAMBDA; ilam++ ) {

	 lam = SNSED.LAMBDA[iepoch][ilam] ;
	 if ( lam < LAMMIN ) { continue ; }
	 if ( lam > LAMMAX ) { continue ; }

	 // get info for this labmda bin

	 magflux_info( 0, ifilt, ilam, iepoch, 
		       &lam, &flux, &wflux, &wfilt ) ;

	 trans  = filter_trans( lam, ifilt, 0 );

	 /* xxxxxxxxxxxxxxxxx
	 if ( iepoch==10 ) {
	   printf(" xxx ifilt=%d  lam=%.1f  flux=%10.3le trans=%.4f \n",
		  ifilt, lam, flux, trans); fflush(stdout);
	 }
	 xxxxxxxxx */

	 if ( trans > 1.0E-20 ) {
	   fluxsum_sn     += (trans * wflux * flux)  ;
	   filtsum_check  += (trans * SNSED.LAMBDA_BINSIZE) ;
	 }

       }  // end of ilam loop

	// compute magnitudes. Zero points implicitly applied since
	// we normalize to prmary star spectrum)

       if ( iepoch == -999 ) {
	 printf(" xxxx ifilt=%d   filtsum=%f  fluxsum_sn=%f \n",
		ifilt, filtsum,fluxsum_sn );
       }

       if ( fluxsum_sn > zero && filtsum > zero ) {
	 arg     = (fluxsum_sn / filtsum) ;
	 mag     = -2.5 * log10(arg)  + zp ;
       }
       else
	 { mag = 0.0; }

       filter_check = filtsum_check/FILTER[ifilt].SUMTOT - 1.0;

       // if filtsum_check (integral over SED-lambda range)
       // does not match integral over entire filter range (to within 20%),
       // set MAG_UNDEFINED to flag problem
       if ( fabs(filter_check) > 0.05  ) { 
	 printf(" WARNING(%s): filter_check(%s) = %.3f  ep=%d\n",
		fnam, FILTER[ifilt].name, filter_check, iepoch );
	 mag = MAG_UNDEFINED ; 
       }

       SNSED.MAG_RST[iepoch][ifilt]  = mag;
       
       epoch = SNSED.EPOCH[iepoch];
       if ( fabs(epoch-TPEAK) < fabs(near_peak)  ) {
	 IEPOCH_SNPEAK[ifilt] = iepoch ;
	 near_peak = epoch-TPEAK ;
	 epoch_peak = epoch ;
       }

       if ( fabs(epoch-T15DAY) < fabs(near_15day)  ) {
	 IEPOCH_SN15DAY[ifilt] = iepoch ;
	 near_15day  = epoch - T15DAY ;
	 epoch_15day = epoch ;
       }


     }  // end of "iepoch" loop

     
     /*
     printf(" xxx FILT=%s  epoch(peak)=%f  epoch(15day)=%f \n",
	    FILTER[ifilt].name, epoch_peak, epoch_15day );
     */


   }  // end of 'ifilt' loop


   // #############################
   // print results to screen
   // Warning: T=0 and T=15 really mean closest epochs;
   // they could by -0.5 and 14.5 days, so don't over-interpret.
   // These mags are used only for screen dump.


   printf(" \n");
   printf("                                         Synth MAG   Synth    \n");
   printf("          filter              (system)   at T=0      DM15    \n" );
   printf("# ------------------------------------------------------------ \n");


   for ( ifilt=1;  ifilt<=NFILTDEF; ifilt++ ) {

     iepoch_peak   = IEPOCH_SNPEAK[ifilt];
     iepoch_15day  = IEPOCH_SN15DAY[ifilt];

     mag   = SNSED.MAG_RST[iepoch_peak][ifilt] ;
     mag15 = SNSED.MAG_RST[iepoch_15day][ifilt] ;
     dm15  = mag15 - mag;

     printf("SNMAG:   %-20.20s (%6s) %8.4f  %8.4f"
	    ,FILTER[ifilt].name
	    ,FILTER[ifilt].MAGSYSTEM_NAME
	    ,mag, dm15
	    );
     
     printf("\n");

   }


   printf("# ---------------------------------------------------------- \n");

   fflush(stdout);
   return SUCCESS;

}  // end of snmag function




// *****************************************************
void primarymag_zp(int iprim ) {

  /***
   fill

   -  FILTER[ifilt].SSUM_PRIM
   -  FILTER[ifilt].MAGFILTER_ZP
   -  PRIMARYSED[iprim].FLUXSUM[ifilt] (erg/s/cm^2)
   -  PRIMARYSED[iprim].ZP[ifilt] 

  Nov 10, 2006: add MAGFILTER_REF to zero points.

  May 29, 2008: print %20s instead of %10s for filter name
                (to allow those LONG HST names )

  Apr 25, 2009: evaluate only the filters in this MAGSYSTEM_INDX = iprim

  ****/

  int ifilt, ilam, INDX, INDX_INPUT ;
  int idump = 0;

  double  arg, lam, trans, flux, wflux, wfilt, mag, ftmp, LAMMIN, LAMMAX;
  double  fluxsum[MXFILTDEF+1], filtsum[MXFILTDEF+1] ;
  bool    USE_FILT ;
  char fnam[] = "primarymag_zp" ;

  /* ------------------- BEGIN --------------------- */

   printf("\n ***** Compute %s flux and ZP in each filter ***** \n", 
	  PRIMARYSED[iprim].MAGSYSTEM_NAME ) ;
   fflush(stdout);

   for ( ifilt=1;  ifilt <= NFILTDEF ; ifilt++ ) {

     INDX       = FILTER[ifilt].MAGSYSTEM_INDX ;
     INDX_INPUT = FILTER[ifilt].MAGSYSTEM_INDX_INPUT;  // 4/2020
     USE_FILT   = (iprim==INDX || iprim == INDX_INPUT);

     /*  xxxx
     printf(" xxx %s: iprim=%d ifilt=%2d INDX=%d->%d USE=%d \n",
	    fnam, iprim, ifilt, INDX_INPUT, INDX, USE_FILT);
     xxxx */

     if ( !USE_FILT ) { continue; }

     fluxsum[ifilt]   = 0.0;
     filtsum[ifilt]   = 0.0;

     LAMMIN = FILTER[ifilt].LAMBDA_MIN;
     LAMMAX = FILTER[ifilt].LAMBDA_MAX;

     ftmp = 0.0;
     for ( ilam=1; ilam <= PRIMARYSED[iprim].NBIN_LAMBDA; ilam++ ) {
       
	lam = PRIMARYSED[iprim].LAMBDA[ilam] ;

	if ( lam < LAMMIN ) { continue ; }
	if ( lam > LAMMAX ) { continue ; }

	// fetch info for this lambda bin
	magflux_info( iprim, ifilt, ilam, 0, 
		      &lam, &flux, &wflux, &wfilt ) ;

	trans  = filter_trans( lam, ifilt, idump );
	 
	 if ( trans > 0.0 ) {
	   fluxsum[ifilt]  += (trans * wflux * flux) ;
	   filtsum[ifilt]  += (trans * wfilt) ;
	 }
      }   // end of ilam loop

      // May 29, 2008: abort if flux is zero
      if ( fluxsum[ifilt] <= 0.0 ) {
	sprintf(c1err,"%s flux = %f for filter = '%s' "
		, PRIMARYSED[iprim].MAGSYSTEM_NAME
		, fluxsum[ifilt], FILTER[ifilt].name ) ; 
	errmsg(SEV_FATAL, 0, fnam, c1err, " ") ;  
      }

      arg   = fluxsum[ifilt] / filtsum[ifilt] ;
      if ( arg <= 0.0 ) {
	sprintf(c1err,"Bad arg=%le for ifilt=%d", arg, ifilt);
	errmsg(SEV_FATAL, 0, fnam, c1err, " ") ;  
      }
      mag   = -2.5*log10(arg) + FILTER[ifilt].MAGSYSTEM_OFFSET; 

      FILTER[ifilt].SSUM_PRIM            = filtsum[ifilt] ;
      PRIMARYSED[iprim].FLUXSUM[ifilt]   = fluxsum[ifilt] ;
      PRIMARYSED[iprim].SYNMAG[ifilt]    = +mag ;
      PRIMARYSED[iprim].ZP[ifilt]        = -mag ;

      FILTER[ifilt].MAGFILTER_ZP = FILTER[ifilt].MAGFILTER_REF ;

      // adjust filter zeropoint except for AB system because
      // since AB mags are zero
      if ( !PRIMARYSED[iprim].IS_AB ) 
	{ FILTER[ifilt].MAGFILTER_ZP += PRIMARYSED[iprim].ZP[ifilt] ; }

   }  // end of 'ifilt' loop

   return;

}  // end of primarymag_zp

// *****************************************************
void primarymag_zp2(int iprim ) {

  // Created Apri 2020 by R.Kessler
  // If magsys input is different than final magsys,
  // apply zp transformation. This applies to kcor-inputs
  // such as 
  //   MAGSYSTEM:  VEGA->AB
  //   MAGSYSTEM:  BD17->AB
  //

  int  ifilt, INDX, INDX_INPUT ;
  double ZP_INPUT, ZP ;
  //  char fnam[] = "primarymag_zp2" ;

  // ---------- BEGIN ------------

  for ( ifilt=1;  ifilt <= NFILTDEF ; ifilt++ ) {
    
    INDX       = FILTER[ifilt].MAGSYSTEM_INDX ;
    INDX_INPUT = FILTER[ifilt].MAGSYSTEM_INDX_INPUT; 

    if ( INDX != iprim      ) { continue; }
    if ( INDX == INDX_INPUT ) { continue; }

    // here do the transform from INDX_INPUT to INDX
    
    ZP       = PRIMARYSED[INDX].ZP[ifilt] ; 
    ZP_INPUT = PRIMARYSED[INDX_INPUT].ZP[ifilt] ; 
    //    FILTER[ifilt].MAGFILTER_ZP  += (ZP_INPUT - ZP);
    FILTER[ifilt].MAGFILTER_REF += (ZP_INPUT - ZP);
    
    /*xxxx
    printf(" xxx %s: iprim=%d ifilt=%2d(%s)  ZP=%6.3f->%6.3f \n",
	   fnam, iprim, ifilt, FILTER[ifilt].name,  ZP_INPUT, ZP);
    */
  }

  return;

} // end primarymag_zp2

// *****************************************************
void primarymag_summary(int iprim) {

  // print table summarizing primary mag and ZP

  char *NAME = PRIMARYSED[iprim].MAGSYSTEM_NAME ;
  int  ifilt, INDX, INDX_INPUT;
  bool USE_FILT;
  //  char fnam[] = "primarymag_summary" ;

  // ----------- BEGIN -----------

  // #############################
  // print results to screen
  
  printf("                                            syn                    \n");
  printf("                                flux        %s       system   final  \n",
	 NAME );
  printf("         filter      (system) (Nph/s/cm^2)  mag       ZPoff     ZP    \n" );
  printf("  ------------------------------------------------------------------ \n");
  
  for ( ifilt=1;  ifilt <= NFILTDEF ; ifilt++ ) {
    
    INDX       = FILTER[ifilt].MAGSYSTEM_INDX ;
    INDX_INPUT = FILTER[ifilt].MAGSYSTEM_INDX_INPUT;  // 4/2020
    USE_FILT = (iprim==INDX);
    if ( !USE_FILT ) { continue; }
    
    printf("%20s (%6s)   %9.3le %8.4f %8.4f %8.4f  \n"  // .xyz
	   ,FILTER[ifilt].name
	   ,FILTER[ifilt].MAGSYSTEM_NAME
	   ,PRIMARYSED[iprim].FLUXSUM[ifilt]
	   ,PRIMARYSED[iprim].SYNMAG[ifilt]
	   ,FILTER[ifilt].MAGFILTER_REF
	   ,FILTER[ifilt].MAGFILTER_ZP
	   );
    
  }

  printf("  ----------------------------------------------------------------- \n");


  fflush(stdout);

  return ;

} // end primarymag_summary


// *****************************************************
void magflux_info( 
		  int iprim         // index of primary
		  ,int ifilt        // filter index
		  ,int ilam         // (I) lambda index
		  ,int iepoch        // (I) epoch index (0=> primary; else SN)
		  ,double *lambda    // (O) lambda (A)
		  ,double *flux      // (O) flux
		  ,double *wflux     // (O) mag wgt for numerator (see below)
		  ,double *wfilt     // (O) mag wgt for denominator
		  ) {
	 
  /***

  The weights wflux and wfilt  for MAGSYSTEM are defined
  such that


         int [ Flux(lam,nu) * Trans(lam,nu) * wflux ]
  mag = -------------------------------------------
         int [  Trans(lam,nu) * wfilt ]


  Dec 8, 2011: remove primary-index options; everything is NU system.

  ***/

  //  char fnam[] = "magflux_info" ;
  int INDX, IDBUG ;
  double  dlam, lam, lam1, lam2, nu, dnu    ;

  // -------------- BEGIN function --------------

  IDBUG = 0;
  if ( ifilt == -222 && ilam > 1500 && ilam < 1510 ) IDBUG = 1;


  // init output args.
  *wflux   = 0.0 ;
  *wfilt   = 0.0 ;
  *flux    = 0.0;

  if ( iepoch == 0 ) {
    *lambda  = PRIMARYSED[iprim].LAMBDA[ilam];
    dlam     = PRIMARYSED[iprim].LAMBDA_BINSIZE ;
  }
  else {
    *lambda  = SNSED.LAMBDA[iepoch][ilam];
    dlam     = SNSED.LAMBDA_BINSIZE ;
  }

  lam      = *lambda ;  // get local "lam" variable for convenience

  // get generic binning info 

  lam1     = lam - 0.5 * dlam ;
  lam2     = lam + 0.5 * dlam ;
  nu       = LIGHT_A / lam ;
  dnu      = LIGHT_A * (1.0/lam1 - 1.0/lam2);

  INDX = FILTER[ifilt].MAGSYSTEM_INDX ;

  *wflux = dnu / ( PLANCK * nu ) ;   // divide by E/photon to get counts
  *wfilt = dnu / ( PLANCK * nu ) ;

  if ( iepoch == 0 ) 
    { *flux = PRIMARYSED[iprim].FLUX_NU[ilam]; }    // erg/s/cm&2/Hz
  else
    { *flux = SNSED.FLUX_NU[iepoch][ilam]; }     // erg/s/cm^2/Hz


}  // end of magsum_weights


// **********************************************
void wr_fits(char *ptrFile) {

  // Created Dec, 2012 by R.Kessler: 
  // write output to FITS file.
  //
  // July 2016: write optional spectrograph table

  int istat ;
  long  NAXIS = 1, NAXES = 0    ;
  fitsfile *fp ;

  char  clobberFile[MXPATHLEN];
  char fnam[] = "wr_fits"  ;

  // ------------- BEGIN --------------

  printf("\n %s: WRITE CALIB/KCOR TO '%s' \n",  fnam, ptrFile);
  fflush(stdout);
  
  sprintf(clobberFile, "!%s", ptrFile);

  istat = 0 ;
  fits_create_file(&fp, clobberFile, &istat) ;
  sprintf(c1err,"Open %s ", ptrFile );
  wr_fits_errorCheck(c1err,istat);

  // ---------------------------

  // create mandatory primary image (length=0)
  istat = 0 ;
  fits_create_img(fp, FLOAT_IMG, NAXIS, &NAXES, &istat) ;
  sprintf(c1err,"Create zero-len primary image" ) ;
  wr_fits_errorCheck(c1err, istat) ;
  
  // load forms
  sprintf(STRFITS.F4,    "1E");
  sprintf(STRFITS.D8,    "1D");
  sprintf(STRFITS.U4,    "1U");
  sprintf(STRFITS.C20,   "20A");
  STRFITS.blank[0]=0;

  // ---------
  wr_fits_HEAD(fp);
  wr_fits_ZPT(fp);          // zpt info vs. filter
  wr_fits_SNSED(fp);        // SN SED: flux vs. Trest
  wr_fits_KCOR(fp);         // KCOR vs. AV,Z,Trest
  wr_fits_MAGS(fp);         // mags and MWXT vs. AV,Z,Trest
  wr_fits_FilterTrans(fp);  // Filter transmissions
  wr_fits_PRIMARY(fp);      // primary SEDs
  wr_fits_SPECTROGRAPH(fp); // optional spectrograph

  fits_close_file(fp, &istat);

  return ;

} // end of wr_fits


// =====================================
void wr_fits_HEAD(fitsfile *fp) {

  // Nov 15 2020: write SURVEY=%s in comment field for each filter.

  int istat, iprim, ifilt, ikcor, N, jbinsize, IVER ;
  char  KEYNAME[40], KEYVAL[MXPATHLEN], MSG[200] ;
  char  fnam[] = "wr_fits_HEAD" ;

  // ------------ BEGIN -----------

  istat = 0;

  printf("\t %s: write header info.\n", fnam);
  fflush(stdout);

  // start with internal version in case we need
  // VERSION-dependent parsing code later.
  IVER = VERSION_KCOR ;
  fits_update_key(fp, TINT, "VERSION",
		  &IVER, "Internal KCOR version", &istat );

  sprintf(c1err,"Write VERSION key in header" ) ;
  wr_fits_errorCheck(c1err, istat) ;

  // Jul 2020: write name of input kcor file
  istat = 0 ;
  sprintf(KEYNAME,"INPUT_FILE");
  sprintf(KEYVAL,"%s", INPUTS.inFile_input);
  fits_update_key(fp, TSTRING, KEYNAME, KEYVAL,
		  "Name of kcor-input file", &istat ); 

  // July 2023: write cwd as well
  istat = 0 ;
  sprintf(KEYNAME,"CWD");
  getcwd(KEYVAL,MXPATHLEN);
  //  sprintf(KEYVAL,"%s", INPUTS.inFile_input);
  fits_update_key(fp, TSTRING, KEYNAME, KEYVAL,
		  "current work dir", &istat ); 

  // -----------------------------
  // write names of primary refs

  N = INPUTS.NPRIMARY ;
  fits_update_key(fp, TINT, "NPRIM",
		  &N, "Number of primary refs", &istat );
  sprintf(c1err,"Write NPRIM key in header" ) ;
  wr_fits_errorCheck(c1err, istat) ;

  for ( iprim = 1; iprim <= N ; iprim++ ) {

    sprintf(KEYNAME,"PRIM%3.3d", iprim);

    jbinsize = (int)SNSED.LAMBDA_BINSIZE ; 
    sprintf(KEYVAL, "%s SED (erg/s/cm^2! per %d A)", 
	    PRIMARYSED[iprim].MAGSYSTEM_NAME, jbinsize );

    istat = 0 ;
    fits_update_key(fp, TSTRING, KEYNAME, KEYVAL,
		    "Name of primary ref", &istat );    
   }


  // write Filter names into global header keys

  fits_update_key(fp, TINT, "NFILTERS",
		  &NFILTDEF, "Number of filters", &istat );
  sprintf(c1err,"Write NFILTERS key in header" ) ;
  wr_fits_errorCheck(c1err, istat) ;


  // July 2020 write filter paths so other codes can find DOCANA notes
  for(ifilt=1; ifilt <= NFILTPATH; ifilt++ ) {
    sprintf(KEYNAME,"FILTPATH%d", ifilt);
    istat = 0 ;
    fits_update_key(fp, TSTRING, KEYNAME, FILTER[ifilt].PATH_ORIG,
		    "Filter PATH", &istat );    
  }

  for ( ifilt = 1; ifilt <= NFILTDEF ; ifilt++ ) {
    sprintf(KEYNAME,"FILT%3.3d", ifilt);
  
    istat = 0 ;
    // xxxx sprintf(MSG,"Filter name; SURVEY=%s", FILTER[ifilt].SURVEY_NAMES);
    sprintf(MSG,"name; SURVEY=%s", FILTER[ifilt].SURVEY_NAMES);
    fits_update_key(fp, TSTRING, KEYNAME, FILTER[ifilt].name,
		    MSG, &istat );    
   }

  if ( NKCOR > 0 ) {
    // RV used for MW extinction and warping
    float RV = (float)INPUTS.RV_MWCOLORLAW;
    fits_update_key(fp, TFLOAT, "RV",
		  &RV, "RV used for CCM warping", &istat );
    
    // Sep 22 2013: write color law option
    fits_update_key(fp, TINT, "OPT_MWCOLORLAW",
		    &INPUTS.OPT_MWCOLORLAW, INPUTS.STR_MWCOLORLAW, &istat );
  }

  // write KCOR names into global header keys

  fits_update_key(fp, TINT, "NKCOR",
		  &NKCOR, "Number of KCOR tables", &istat );

  for ( ikcor = 1; ikcor <= NKCOR ; ikcor++ ) {
    sprintf(KEYNAME,"KCOR%3.3d", ikcor);
  
    sprintf(KEYVAL,"Kcor %s for rest %s to obs %s",
	    KCORSYM[ikcor], KCORLIST[ikcor][0], KCORLIST[ikcor][1] );

    istat = 0 ;
    fits_update_key(fp, TSTRING, KEYNAME, KEYVAL,
		    "KCOR definition", &istat );    
   }


  // =====================================
  // write binning info
  float FBIN, FMIN, FMAX ;
  int  NBIN ;

  // --- LAMBDA ------
  NBIN = SNSED.NBIN_LAMBDA  ;
  FBIN = SNSED.LAMBDA_BINSIZE ;
  FMIN = SNSED.LAMBDA_MIN;
  FMAX = SNSED.LAMBDA_MAX;

  fits_update_key(fp, TINT, "NBL", &NBIN,
		  "Number of lambda bins", &istat );
  fits_update_key(fp, TFLOAT, "LBIN", &FBIN,
		  "lambda bin size (A)", &istat );
  fits_update_key(fp, TFLOAT, "LMIN", &FMIN,
		  "min lambda (A)", &istat );
  fits_update_key(fp, TFLOAT, "LMAX", &FMAX,
		  "max lambda (A)", &istat );


  // --- Trest (Epoch)------

  NBIN = SNSED.NEPOCH ;
  FMIN = SNSED.TREST_MIN;  
  FMAX = SNSED.TREST_MAX; 
  FBIN = (FMAX-FMIN)/(float)(NBIN-1) ;

  fits_update_key(fp, TINT, "NBT", &NBIN,
		  "Number of Trest bins", &istat );
  fits_update_key(fp, TFLOAT, "TBIN", &FBIN,
		  "Trest bin size (days)", &istat );
  fits_update_key(fp, TFLOAT, "TMIN", &FMIN,
		  "min Trest (days)", &istat );
  fits_update_key(fp, TFLOAT, "TMAX", &FMAX,
		  "max Trest (days)", &istat );


  // --- redshift ------

  NBIN = INPUTS.NBIN_REDSHIFT ;
  FMIN = INPUTS.REDSHIFT_MIN  ;
  FMAX = INPUTS.REDSHIFT_MAX  ;
  if ( NBIN > 1 ) 
    {  FBIN = (FMAX-FMIN)/(float)(NBIN-1) ; }
  else
    { FBIN = 0.0 ; }

  fits_update_key(fp, TINT, "NBZ", &NBIN,
		  "Number of redshift bins", &istat );
  fits_update_key(fp, TFLOAT, "ZBIN", &FBIN,
		  "redshift bin size", &istat );
  fits_update_key(fp, TFLOAT, "ZMIN", &FMIN,
		  "min redshift", &istat );
  fits_update_key(fp, TFLOAT, "ZMAX", &FMAX,
		  "max redshift", &istat );


  // --- AVwarp ------

  NBIN = INPUTS.NBIN_AV ;
  FMIN = INPUTS.AV_MIN  ;
  FMAX = INPUTS.AV_MAX  ;
  FBIN = INPUTS.AV_BINSIZE ;

  fits_update_key(fp, TINT, "NBAV", &NBIN,
		  "Number of AV-warp bins", &istat );
  fits_update_key(fp, TFLOAT, "AVBIN", &FBIN,
		  "AV-warp bin size", &istat );
  fits_update_key(fp, TFLOAT, "AVMIN", &FMIN,
		  "min AV-warp", &istat );
  fits_update_key(fp, TFLOAT, "AVMAX", &FMAX,
		  "max AV-warp", &istat );
  
  // Nov 16 2014: write non-zero lambda shifts for filters.
  // These are just comments and are not read by other programs.

  N = INPUTS.NFILTER_LAMSHIFT ;
  fits_update_key(fp, TINT, "NFILTER_LAMSHIFT", &N,
		  "Nfilt with non-zero lamba-shift", &istat );

  float lamshift ;
  for ( ifilt = 1; ifilt <= NFILTDEF ; ifilt++ ) {
    lamshift = (float)INPUTS.FILTER_LAMSHIFT[ifilt];
    if ( lamshift == 0.0 ) { continue ; }

    sprintf(KEYNAME,"LAMSHIFT(%s)", FILTER[ifilt].name );
  
    istat = 0 ;
    fits_update_key(fp, TFLOAT, KEYNAME, &lamshift,
		    "wavelength shift (Angstroms)", &istat );
   }


  // write optional spectrograph instrument
  if ( SPECTROGRAPH_USEFLAG ) {
    istat = 0 ;

    sprintf(KEYNAME,"SPECTROGRAPH_INSTRUMENT" );  
    sprintf(KEYVAL,"%s", INPUTS_SPECTRO.INSTRUMENT_NAME);
    fits_update_key(fp, TSTRING, KEYNAME, KEYVAL,
		    "Spectrograph name", &istat );    

    sprintf(KEYNAME,"SPECTROGRAPH_NBLAM" );  
    fits_update_key(fp, TINT, KEYNAME, &INPUTS_SPECTRO.NBIN_LAM,
		    "Spectrograph NBLAM", &istat );    

    sprintf(KEYNAME,"SPECTROGRAPH_FILTERLIST" );  
    sprintf(KEYVAL,"%s", INPUTS_SPECTRO.SYN_FILTERLIST_BAND );
    fits_update_key(fp, TSTRING, KEYNAME, KEYVAL,
		    "", &istat );    
  }

  return ;

} // end of wr_fits_HEAD


// ========================================
void wr_fits_SNSED(fitsfile *fp) {
  
  char 
    fnam[] = "wr_fits_SNSED"
    ,LABEL_FLUX[]  = "SN Flux (erg/s/cm^2/A)"
    ,TBLname[]     = "SN SED";

  long NROW = 0;
  int i, ncol, istat;

  int NBLAM, NTREST, ilam, iday, icol ;
  int firstrow, firstelem, nrow ;
  float lam, epoch, flux, FBIN, FMIN, FMAX ;

  // --------------- BEGIN ------------

  printf("\t %s: write SN SED \n", fnam );
  fflush(stdout);

  ncol = 1;  istat = 0;
  /*
  STRFITS.tName[0] = LABEL_LAM   ;
  STRFITS.tName[1] = LABEL_TREST ;
  */
  STRFITS.tName[0] = LABEL_FLUX  ;

  for(i=0; i<ncol; i++ ) {
    STRFITS.tForm[i] = STRFITS.F4 ;
    STRFITS.tUnit[i] = STRFITS.blank ;
  }
    
  fits_create_tbl(fp, BINARY_TBL, NROW, ncol
		  ,&STRFITS.tName[0]
		  ,&STRFITS.tForm[0]
		  ,&STRFITS.tUnit[0]
		  ,TBLname, &istat );
    
  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  // write binning in header to double-check the binning
  // in the main header.


  // write lambda binning
  NBLAM  = SNSED.NBIN_LAMBDA;
  iday   = 1;
  FMIN   = SNSED.LAMBDA[iday][1];
  FMAX   = SNSED.LAMBDA[iday][NBLAM];
  FBIN   = SNSED.LAMBDA_BINSIZE ;

  fits_update_key(fp, TINT, "NBL",
		  &NBLAM, "Nubmer of Lambda bins", &istat );
  sprintf(c1err,"write NBLAM header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  fits_update_key(fp, TFLOAT, "LMIN",
		  &FMIN, "min Lambda", &istat );
  sprintf(c1err,"write LMIN header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  fits_update_key(fp, TFLOAT, "LMAX",
		  &FMAX, "max Lambda", &istat );
  sprintf(c1err,"write LMAX header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  fits_update_key(fp, TFLOAT, "LBIN",
		  &FBIN, "Lambda bin size", &istat );
  sprintf(c1err,"write LBIN header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  // write TREST binning

  NTREST = SNSED.NEPOCH ;

  FMIN   = SNSED.EPOCH[1] ;
  FMAX   = SNSED.EPOCH[NTREST] ;
  if ( NTREST > 1 ) 
    { FBIN = (FMAX-FMIN)/(float)(NTREST-1); }
  else
    { FBIN = 0.0 ; }


  fits_update_key(fp, TINT, "NBT",
		  &NTREST, "Nubmer of Trest bins", &istat );
  sprintf(c1err,"write NTREST header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  fits_update_key(fp, TFLOAT, "TMIN",
		  &FMIN, "min Trest", &istat );
  sprintf(c1err,"write TMIN header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  fits_update_key(fp, TFLOAT, "TMAX",
		  &FMAX, "max Trest", &istat );
  sprintf(c1err,"write TMAX header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  fits_update_key(fp, TFLOAT, "TBIN",
		  &FBIN, "Trest bin size", &istat );
  sprintf(c1err,"write TBIN header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  //------ fill SED table

  firstelem = nrow = 1;
  firstrow  = istat = 0;

  for ( iday=1; iday <= NTREST ; iday++ ) {      
    for ( ilam=1; ilam <= NBLAM; ilam++ ) {
      lam   = SNSED.LAMBDA[iday][ilam];
      epoch = SNSED.EPOCH[iday];
      flux  = (float)SNSED.FLUX_WAVE[iday][ilam] ;
      flux *= SNSED.LAMBDA_BINSIZE;     
      firstrow++ ;
      icol = 1 ;
      fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		     &flux, &istat);        
    }
   }   

  return ;

} // end of wr_fits_SNSED


// ========================================
void wr_fits_ZPT(fitsfile *fp) {

  // write zeropoint information.
  // ZPoff(prim) is the zeropoint offset computed internally
  // (mag_native - mag_synth). It is stored in the fits file only
  // to be printed by snana; snana.exe does not use this information 
  // since it is used to already to compute K-corrections.
  //
  // ZPoff(SNphot) are optional filter-dependent offsets defined in the
  // ZPOFF.DAT file in the /filters/xxx subdir. Example is the SDSS AB 
  // offsets. These offsets are NOT used here in the K-cor program;
  // instead they are stored so that snana can use them.
  //

  int istat, ncol, icol, ifilt ;
  int firstrow, firstelem, nrow ;

  long NROW      = 0 ;
 
  char 
    fnam[] = "wr_fits_ZPT"
    ,LABEL_FILT[]         = "Filter Name"
    ,LABEL_PRIMARY_NAME[] = "Primary Name"
    ,LABEL_PRIMARY_MAG[]  = "Primary Mag"
    ,LABEL_ZPOFF_REF[]    = "ZPoff(Primary)"  // used here, but not in snana
    ,LABEL_ZPOFF_SN[]     = "ZPoff(SNpot)"    // ignore here, used in snana
    ,TBLname[] = "ZPoff" 
    ;

  //     iprim =  FILTER[ifilt].MAGSYSTEM_INDX ;

  // ------------ BEGIN ---------

  printf("\t %s: write ZPT info \n", fnam );
  fflush(stdout);

  ncol = 5 ; istat = 0 ;

  STRFITS.tName[0] = LABEL_FILT ;
  STRFITS.tName[1] = LABEL_PRIMARY_NAME ;
  STRFITS.tName[2] = LABEL_PRIMARY_MAG ;
  STRFITS.tName[3] = LABEL_ZPOFF_REF ;
  STRFITS.tName[4] = LABEL_ZPOFF_SN ;



  for(icol=0; icol < ncol; icol++ ) {
    STRFITS.tForm[icol] = STRFITS.F4 ;
    STRFITS.tUnit[icol] = STRFITS.blank ;
  }
  STRFITS.tForm[0] = STRFITS.C20 ;
  STRFITS.tForm[1] = STRFITS.C20 ;


  fits_create_tbl(fp, BINARY_TBL, NROW, ncol
		  ,&STRFITS.tName[0]
		  ,&STRFITS.tForm[0]
		  ,&STRFITS.tUnit[0]
		  ,TBLname, &istat );

  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  // load table.

  firstelem = nrow = 1;
  firstrow = 0 ;

  float mag_prim, zptoff_ref, zptoff_filt ;
  char  *name_filt, *name_prim ;

  for ( ifilt = 1; ifilt <= NFILTDEF ; ifilt++) {

    firstrow++ ;  icol=0;
    name_filt    = FILTER[ifilt].name ;
    name_prim    = FILTER[ifilt].MAGSYSTEM_NAME ; 
    mag_prim     = FILTER[ifilt].MAGFILTER_REF ;     // primary mag
    zptoff_ref   = FILTER[ifilt].MAGFILTER_ZP ;      // applied filter zp
    zptoff_filt  = FILTER[ifilt].MAGFILTER_ZPOFF;    // stored AB  offset
   
    icol++;     
    fits_write_col(fp, TSTRING, icol, firstrow, firstelem, nrow,
		   &name_filt, &istat);  
    sprintf(c1err,"write filter '%s' into %s", name_filt, TBLname );
    wr_fits_errorCheck(c1err, istat) ;

    icol++;     
    fits_write_col(fp, TSTRING, icol, firstrow, firstelem, nrow,
		   &name_prim, &istat);  
    sprintf(c1err,"write primary '%s' into %s", name_prim, TBLname );
    wr_fits_errorCheck(c1err, istat) ;

    icol++ ;     
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &mag_prim, &istat);
    sprintf(c1err,"write MAG_PRIM into %s", TBLname );
    wr_fits_errorCheck(c1err, istat) ;  

    icol++ ;     
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &zptoff_ref, &istat);
    sprintf(c1err,"write ZPTOFF_REF into %s", TBLname );
    wr_fits_errorCheck(c1err, istat) ;    

    icol++ ;
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &zptoff_filt, &istat);  
    sprintf(c1err,"write ZPTOFF_FILT into %s", TBLname );
    wr_fits_errorCheck(c1err, istat) ;  

  }

  return ;

} // end of wr_fits_ZPT



// ========================================
void wr_fits_PRIMARY(fitsfile *fp) {

  // Mar 2019: for unused primary, fix to write NAME instead of blank
  //   so that AstroPy doesn't choke. Also write zero flux.


  int istat, ncol, iprim, ilam, NBLAM, NPRIM ;
  int icol,  firstrow, firstelem, nrow ;

  long NROW      = 0 ;
  char 
    fnam[]       = "wr_fits_PRIMARY" 
    ,LABEL_LAM[] = "wavelength (A)"
    ,TBLname[]   = "PrimarySED"  
    ;

  char *NAME;
  float lam, flux ;
  double lam8, flux8 ;

  // ----------- BEGIN ------------

  printf("\t %s: write Primary SED(s) \n", fnam);
  fflush(stdout);

  // -----------------------------------------
  // create table, first column is LAM, followed by 
  // SED-flux columns

  NPRIM = INPUTS.NPRIMARY;

  ncol = NPRIM+1 ;  istat = 0;
  for(iprim = 0; iprim <= NPRIM; iprim++ ) {
    NAME = INPUTS.name_PRIMARY[iprim];

    if ( iprim == 0 ) 
      {  STRFITS.tName[iprim] = LABEL_LAM ; }
    else 
      {  STRFITS.tName[iprim] = NAME ; }

    STRFITS.tForm[iprim] = STRFITS.F4 ;
    STRFITS.tUnit[iprim] = STRFITS.blank ;
  }

  fits_create_tbl(fp, BINARY_TBL, NROW, ncol
		  ,&STRFITS.tName[0]
		  ,&STRFITS.tForm[0]
		  ,&STRFITS.tUnit[0]
		  ,TBLname, &istat );

  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  // -----------------------------------------
  // fill table
  firstelem = firstrow = nrow = 1 ;
  NBLAM = SNSED.NBIN_LAMBDA;  

  for ( ilam=1; ilam <= NBLAM; ilam++ ) {
    lam8 = SNSED.LAMBDA[1][ilam];
    lam  = (float)lam8 ;

    firstrow = ilam ;

    // start by filling wavelength
    icol = 1 ;
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &lam, &istat);  

    sprintf(c1err,"writing primary ilam=%d (lam=%f)", ilam, lam );
    wr_fits_errorCheck(c1err, istat) ;

    // now fill SED flux for each primary ref.
    for(iprim=1 ; iprim <= NPRIM; iprim++ ) {
      icol++ ;           

      if ( PRIMARYSED[iprim].USE ) 
	{ flux8  = primaryflux ( iprim, lam8 );  }
      else
	{ flux8 = 0.0 ; }

      flux   = (float)flux8 ;  // per A.

      fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		     &flux, &istat) ;

      sprintf(c1err,"writing primary flux for iprim=%d (lam=%f)", 
	      iprim, lam );
      wr_fits_errorCheck(c1err, istat) ;
      
    } // iprim

  }  // ilam

  return ;

} // end of wr_fits_PRIMARY



// ========================================
void wr_fits_FilterTrans(fitsfile *fp) {

  // Write filter trans on same wavelength grid as SNSED.
  // July 2016: write IFLAG_SYN for synthetic filter from spectrograph

  int istat, ncol, ifilt, ilam, NBLAM ;
  int icol,  firstrow, firstelem, nrow ;

  long NROW      = 0 ;
  char 
    fnam[]       = "wr_fits_FilterTrans" 
    ,LABEL_LAM[] = "wavelength (A)"
    ,TBLname[]   = "FilterTrans"  
    ;

  float lam4, trans4 ;
  double lam8, trans8 ;

  // ----------- BEGIN ------------


  printf("\t %s: write filter transmissions. \n", fnam);
  fflush(stdout);

  // -----------------------------------------
  // create filter table. First column is LAMBDA, followed by 
  // transmission columns (1 column per filter)


  ncol = NFILTDEF+1 ;  istat = 0;
  for(ifilt = 0; ifilt <= NFILTDEF; ifilt++ ) {
    if ( ifilt == 0 ) 
      {  STRFITS.tName[ifilt] = LABEL_LAM ; }
    else 
      {  STRFITS.tName[ifilt] = FILTER[ifilt].name ; }

    STRFITS.tForm[ifilt] = STRFITS.F4 ;
    STRFITS.tUnit[ifilt] = STRFITS.blank ;
  }

  fits_create_tbl(fp, BINARY_TBL, NROW, ncol
		  ,&STRFITS.tName[0]
		  ,&STRFITS.tForm[0]
		  ,&STRFITS.tUnit[0]
		  ,TBLname, &istat );

  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;


  // -----------------------------------------
  // fill table
  firstelem = firstrow = nrow = 1 ;
  NBLAM = SNSED.NBIN_LAMBDA;
  int NONZERO = 0 ;

  for ( ilam=1; ilam <= NBLAM; ilam++ ) {
    lam8 = SNSED.LAMBDA[1][ilam];
    lam4 = (float)lam8 ;

    if ( lam8 > 0 ) { NONZERO++ ; }
    firstrow = ilam ;

    // start by filling wavelength
    icol = 1 ;
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &lam4, &istat);  

    sprintf(c1err,"writing filter ilam=%d (lam=%f)", ilam, lam4 );
    wr_fits_errorCheck(c1err, istat) ;

    // now fill transmission for each filter
    for(ifilt=1 ; ifilt <= NFILTDEF; ifilt++ ) {
      icol++ ;           
      trans8 = filter_trans( lam8, ifilt, 0 );
      trans4 = (float)trans8 ;
            
      fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		     &trans4, &istat) ;
      
      sprintf(c1err,"writing Trans for  ifilt=%d (lam=%f)", ifilt, lam4 );
      wr_fits_errorCheck(c1err, istat) ;

    } // ifilt

  }  // ilam


  if ( NONZERO == 0 ) {
    sprintf(c1err,"All lam=0 for filter trans");
    sprintf(c2err,"Something is really messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err,  c2err) ;  
  }

  return ;

} // end of wr_fits_FilterTrans



// ========================================
void wr_fits_KCOR(fitsfile *fp) {

  // Dec 2012
  // create and fill kcor table containing
  // AV Z Trest  K_X1Y2  K_X2Y2  K_X3Y3 ...

  int istat, ncol, ikcor, iday, iz, iav ;
  int icol,  firstrow, firstelem, nrow ;
  int NCOFF = 3;

  long NROW      = 0 ;
  char 
    fnam[]       = "wr_fits_KCOR" 
    ,LABEL_T[]   = "Trest"
    ,LABEL_Z[]   = "Redshift"
    ,LABEL_AV[]  = "AVwarp"
    ,TBLname[]   = "KCOR"  
    ;

  // ----------- BEGIN ------------


  printf("\t %s: write KCOR tables \n", fnam);
  fflush(stdout);

  // -----------------------------------------
  // create  table, first three columns are T,Z,Av;
  // then the K_ZY columns follow.


  ncol = NKCOR + NCOFF ;  istat = 0;


  STRFITS.tName[1] = LABEL_T  ; 
  STRFITS.tName[2] = LABEL_Z  ; 
  STRFITS.tName[3] = LABEL_AV ; 

  for(ikcor = 1; ikcor <= ncol; ikcor++ ) {
    STRFITS.tForm[ikcor] = STRFITS.F4 ;
    STRFITS.tUnit[ikcor] = STRFITS.blank ;
  }

  for(ikcor = 1; ikcor <= NKCOR; ikcor++ ) 
    { STRFITS.tName[NCOFF+ikcor] = KCORSYM[ikcor] ; }

  fits_create_tbl(fp, BINARY_TBL, NROW, ncol
		  ,&STRFITS.tName[1]
		  ,&STRFITS.tForm[1]
		  ,&STRFITS.tUnit[1]
		  ,TBLname, &istat );


  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  if ( NKCOR <= 0 ) { return ; }
  // ----------
  // fill table

  float av, z, epoch, dum, kcor ;

  firstelem = nrow = 1 ;
  firstrow = 0;

  
  for ( iav = 1; iav <= INPUTS.NBIN_AV; iav++ ) {

    dum    = (double)(iav-1) ;
    av     = INPUTS.AV_MIN + dum * INPUTS.AV_BINSIZE;

    for ( iz=1; iz<=INPUTS.NBIN_REDSHIFT; iz++ ) {
      for ( iday=1; iday<=SNSED.NEPOCH; iday++ ) {

	z        = R4KCOR_GRID.REDSHIFT[1][iav][iz][iday] ;
	epoch    = R4KCOR_GRID.EPOCH[1][iav][iz][iday] ;

	firstrow++ ;

	icol = 1 ;
	fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		       &epoch, &istat);  
	sprintf(c1err,"writing EPOCH=%f (iav,iz,iday=%d,%d,%d)", 
		epoch, iav, iz, iday);
	wr_fits_errorCheck(c1err, istat) ;

	icol = 2 ;
	fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		       &z, &istat);  
	sprintf(c1err,"writing REDSHIFT=%f (iav,iz,iday=%d,%d,%d)", 
		z, iav, iz, iday);
	wr_fits_errorCheck(c1err, istat) ;

	icol = 3 ;
	fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		       &av, &istat); 
	sprintf(c1err,"writing AV=%f (iav,iz,iday=%d,%d,%d)", 
		av, iav, iz, iday);
	wr_fits_errorCheck(c1err, istat) ;

	// and now the kcor values
	for(ikcor = 1; ikcor <= NKCOR; ikcor++ ) {
	  icol++ ;
	  kcor   = R4KCOR_GRID.VALUE[ikcor][iav][iz][iday] ;
	  fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
			 &kcor, &istat); 

	  sprintf(c1err,"writing KCOR=%f (iav,iz,iday,ikcor=%d,%d,%d,%d)", 
		  kcor, iav, iz, iday,ikcor);
	  wr_fits_errorCheck(c1err, istat) ;

	} // ikcor
	
      } // iday
    }   // iz
  } // iav

  return ;
} // end of wr_fits_KCOR



// ========================================
void wr_fits_MAGS(fitsfile *fp) {

  // Dec 2012
  // create and fill mag table containing
  // AV Z Trest  mag(filt1) dmag/dMWX(filt1) ... etc for each filter
  // The AV,Z,Trest binning is the same as for the KCOR table.

  int istat, ncol, ifilt, iday, iz, iav ;
  int icol,  firstrow, firstelem, nrow, lenf ;
  int NCOFF = 3 ;

  long NROW      = 0 ;
  char 
    fnam[]       = "wr_fits_MAGS" 
    ,LABEL_T[]   = "Trest"
    ,LABEL_Z[]   = "Redshift"
    ,LABEL_AV[]  = "AVwarp"
    ,TBLname[]   = "MAG+MWXTCOR" 
    ,LABELS[2][MXFILTINDX+1][80]
    ,MAGSYS[60], CFILT[2]
    ;

  // ----------- BEGIN ------------


  printf("\t %s: write MAG+MWXT tables \n", fnam);
  fflush(stdout);

  // -----------------------------------------
  // create  table, first three columns are T,Z,Av;
  // then the K_ZY columns follow.


  ncol = 2*NFILTDEF + NCOFF ;  istat = 0;

  STRFITS.tName[1] = LABEL_T  ; 
  STRFITS.tName[2] = LABEL_Z  ; 
  STRFITS.tName[3] = LABEL_AV ; 

  for(ifilt = 1; ifilt <= ncol; ifilt++ ) { 
    STRFITS.tForm[ifilt] = STRFITS.F4 ;
    STRFITS.tUnit[ifilt] = STRFITS.blank ;
  }

  for(ifilt = 1; ifilt <= NFILTDEF; ifilt++ ) {

    lenf = strlen(FILTER[ifilt].name);
    sprintf(CFILT, "%c", FILTER[ifilt].name[lenf-1]   ) ;
    sprintf(MAGSYS,"%s", FILTER[ifilt].MAGSYSTEM_NAME ) ;

    // define MAGOBS
    icol = NCOFF + ifilt ;
    sprintf(LABELS[0][ifilt], "MAGOBS_%s(%s)", CFILT, MAGSYS);
    STRFITS.tName[icol] = LABELS[0][ifilt];

    // now the MWXT corrections
    icol = NCOFF + ifilt + NFILTDEF;
    sprintf(LABELS[1][ifilt], "MWXT_SLOPE_%s(%s)", CFILT, MAGSYS);
    STRFITS.tName[icol] = LABELS[1][ifilt];
  }

  fits_create_tbl(fp, BINARY_TBL, NROW, ncol
		  ,&STRFITS.tName[1]
		  ,&STRFITS.tForm[1]
		  ,&STRFITS.tUnit[1]
		  ,TBLname, &istat );

  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  if ( NKCOR <= 0 ) { return ; }

  // ----------
  // fill table

  float av, z, epoch, dum, MAGOBS, MWXTSLP ;

  firstelem = nrow = 1 ;
  firstrow = 0;

  for ( iav = 1; iav <= INPUTS.NBIN_AV; iav++ ) {

    dum    = (double)(iav-1) ;
    av     = INPUTS.AV_MIN + dum * INPUTS.AV_BINSIZE;

    for ( iz=1; iz<=INPUTS.NBIN_REDSHIFT; iz++ ) {
      for ( iday=1; iday<=SNSED.NEPOCH; iday++ ) {

	z        = R4KCOR_GRID.REDSHIFT[1][iav][iz][iday] ; 
	epoch    = R4KCOR_GRID.EPOCH[1][iav][iz][iday];

	firstrow++ ;

	icol = 1 ;
	fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		       &epoch, &istat);  
	sprintf(c1err,"writing EPOCH=%f (iav,iz,iday=%d,%d,%d)", 
		epoch, iav, iz, iday);
	wr_fits_errorCheck(c1err, istat) ;

	icol = 2 ;
	fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		       &z, &istat);  
	sprintf(c1err,"writing REDSHIFT=%f (iav,iz,iday=%d,%d,%d)", 
		z, iav, iz, iday);
	wr_fits_errorCheck(c1err, istat) ;

	icol = 3 ;
	fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		       &av, &istat); 
	sprintf(c1err,"writing AV=%f (iav,iz,iday=%d,%d,%d)", 
		av, iav, iz, iday);
	wr_fits_errorCheck(c1err, istat) ;

	// the MAGOBS values
	for(ifilt = 1; ifilt <= NFILTDEF; ifilt++ ) {
	  icol++ ;
	  MAGOBS      = SNSED.R4MAG_OBS[0][ifilt][iav][iz][iday];
	  fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
			 &MAGOBS, &istat); 
	  sprintf(c1err,"writing MAGOBS=%f (iav,iz,iday,ifilt=%d,%d,%d,%d)", 
		  MAGOBS, iav, iz, iday,ifilt);
	  wr_fits_errorCheck(c1err, istat) ;


	  // xxxxxxxxxx
	  if ( iday==-20 && iz==4 && ifilt==2) {
	    printf(" xxx %s: ifilt=%2d  av=%6.3f  MAGOBS=%7.3f \n", 
		   fnam, ifilt, av, MAGOBS);
	    fflush(stdout);
	  }
	  //xxxxxxxxx


	} // ifilt

	// the MWXTSLP values
	for(ifilt = 1; ifilt <= NFILTDEF; ifilt++ ) {
	  icol++ ;
	  MWXTSLP  = SNSED.MW_dXT_dEBV[ifilt][iav][iz][iday];

	  fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
			 &MWXTSLP, &istat); 
	  sprintf(c1err,"writing MWXTSLP=%f (iav,iz,iday,ifilt=%d,%d,%d,%d)", 
		  MWXTSLP, iav, iz, iday,ifilt);
	  wr_fits_errorCheck(c1err, istat) ;
	} // ifilt
	
      } // iday
    }   // iz
  } // iav

  return ;

} // end of wr_fits_MAGS


// ==================================================
void wr_fits_SPECTROGRAPH(fitsfile *fp) {

  // Created July 2016
  // Columns are
  // LAMMIN  LAMMAX  ZP1 SQSKYSIG1  ZP2 SQSKYSIG2  etc ...
  //

  int NBL  = INPUTS_SPECTRO.NBIN_LAM ;
  int NBT  = INPUTS_SPECTRO.NBIN_TEXPOSE ;
  int NCOL = NCOL_noSNR + NBT*2 ;

  char TBLname[40] ;

  long NROW = 0;
  int i, istat, l,t ;

  int firstrow, firstelem, nrow, icol, ifilt ;
  float LAM_f ;
  char  keyName[40], *name;
  char  LABEL_ZP[MXTEXPOSE_SPECTROGRAPH][40];
  char  LABEL_SQ[MXTEXPOSE_SPECTROGRAPH][40];
  char  LABEL_LAMMIN[] = "LAMBIN_MIN" ;
  char  LABEL_LAMMAX[] = "LAMBIN_MAX" ;
  char  LABEL_LAMRES[] = "LAMRES" ;
  
  char fnam[]     = "wr_fits_SPECTROGRAPH" ;

  // ----------- BEGIN ------------

  if ( NBL == 0 ) { return ; }

  printf("\t %s: write SPECTROGRAPH tables \n", fnam);
  fflush(stdout);

  sprintf(TBLname, "%s", FITSTABLE_NAME_SPECTROGRAPH );

  istat = icol = 0 ;

  // set name of each column
  icol++ ; 
  STRFITS.tName[icol] = LABEL_LAMMIN ;
  STRFITS.tForm[icol] = STRFITS.D8 ; 

  icol++ ;
  STRFITS.tName[icol] = LABEL_LAMMAX ;
  STRFITS.tForm[icol] = STRFITS.D8 ; 

  icol++ ;
  STRFITS.tName[icol] = LABEL_LAMRES ;
  STRFITS.tForm[icol] = STRFITS.D8 ; 

  if ( icol != NCOL_noSNR ) {
    sprintf(c1err,"icol=%d, but expect ", icol);
    sprintf(c2err,"%d columns before SNR values", NCOL_noSNR);
    errmsg(SEV_FATAL, 0, fnam, c1err,  c2err) ;  
  }

  for(t=0; t < NBT; t++ ) {
    icol++ ;

    sprintf(LABEL_ZP[t], "ZP%2.2d", t  );
    STRFITS.tName[icol] = LABEL_ZP[t];
    STRFITS.tForm[icol] = STRFITS.F4 ; 

    icol++ ; 
    sprintf(LABEL_SQ[t], "SQSKYSIG%2.2d", t );
    STRFITS.tName[icol] = LABEL_SQ[t] ;
    STRFITS.tForm[icol] = STRFITS.F4 ; 
  }

  // load cast for each column
  for(i=1; i <= NCOL ; i++ ) { STRFITS.tUnit[i] = STRFITS.blank ;  }   

  // create table
  fits_create_tbl(fp, BINARY_TBL, NROW, NCOL
		  ,&STRFITS.tName[1]
		  ,&STRFITS.tForm[1]
		  ,&STRFITS.tUnit[1]
		  ,TBLname, &istat );
    
  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  // write a few table-header keys
  fits_update_key(fp, TINT, "NBL",
		  &NBL, "Nubmer of Lambda bins", &istat );
  sprintf(c1err,"write NBL header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  fits_update_key(fp, TINT, "NBT",
		  &NBT, "Nubmer of Texpose bins", &istat );
  sprintf(c1err,"write NBT header key for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  for(t=0 ; t < NBT; t++ ) {
    float t_f = INPUTS_SPECTRO.TEXPOSE_LIST[t]; 

    sprintf(keyName,"TEXPOSE%2.2d", t);
    fits_update_key(fp, TFLOAT, keyName,
		  &t_f, "Texpose (sec)", &istat );
    sprintf(c1err,"write TEXPOSE%2.2d for %s", t, TBLname );
    wr_fits_errorCheck(c1err, istat) ;
  }

  // fill table
  firstelem = nrow = 1;
  firstrow  = istat = 0;
  double LAM0_d, LAM1_d, LAMRES_d ;
  float  ZP_f, SQSIGSKY_f;
  for(l=0; l < NBL; l++ ) {
    
    firstrow++ ;

    icol = 1 ;
    LAM0_d  = INPUTS_SPECTRO.LAMMIN_LIST[l] ;
    fits_write_col(fp, TDOUBLE, icol, firstrow, firstelem, nrow,
		   &LAM0_d, &istat);   
    sprintf(c1err,"writing LAMMIN_LIST=%f", LAM0_d);
    wr_fits_errorCheck(c1err, istat) ;

    icol = 2 ;
    LAM1_d  = INPUTS_SPECTRO.LAMMAX_LIST[l] ;
    fits_write_col(fp, TDOUBLE, icol, firstrow, firstelem, nrow,
		   &LAM1_d, &istat) ;            
    sprintf(c1err,"writing LAMMAX_LIST=%f", LAM1_d);
    wr_fits_errorCheck(c1err, istat) ;

    icol = 3 ;  // write LAMRES (Oct 14 2016)
    LAMRES_d  = INPUTS_SPECTRO.LAMSIGMA_LIST[l] ;
    fits_write_col(fp, TDOUBLE, icol, firstrow, firstelem, nrow,
		   &LAMRES_d, &istat) ;            
    sprintf(c1err,"writing LAMSIGMA_LIST=%f", LAMRES_d);
    wr_fits_errorCheck(c1err, istat) ;

    for(t=0; t < NBT; t++ ) {

      ZP_f       = (float)INPUTS_SPECTRO.ZP[l][t];
      SQSIGSKY_f = (float)INPUTS_SPECTRO.SQSIGSKY[l][t];

      icol++ ;
      fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		     &ZP_f, &istat);  
      sprintf(c1err,"writing ZP%2.2d = %f", t, ZP_f);
      wr_fits_errorCheck(c1err, istat) ;

      icol++ ;
      fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		     &SQSIGSKY_f, &istat);   
      sprintf(c1err,"writing SQSIGKSY%2.2d = %f", t, SQSIGSKY_f);
      wr_fits_errorCheck(c1err, istat) ;

    } // end t loop over T_EXPOSE bins
  }   // end l lopo over lambda 


  // - - - - - - - - - - - - - - - - - - - - - - - -
  // write another table of LAMRANGE vs. filter 
  // - - - - - - - - - - - - - - - - - - - - - - - -

  sprintf(TBLname, "SYN_FILTER_SPECTROGRAPH" );

  istat = icol = 0 ;

  // set name and cast of each column
  icol++ ; 
  sprintf(STRFITS.tName[icol],"NAME"); 
  STRFITS.tForm[icol] = STRFITS.C20 ;

  icol++ ; 
  sprintf(STRFITS.tName[icol],"LAMMIN"); 
  STRFITS.tForm[icol] = STRFITS.F4 ;

  icol++ ; 
  sprintf(STRFITS.tName[icol],"LAMMAX"); 
  STRFITS.tForm[icol] = STRFITS.F4 ;

  NCOL = icol;

  for(i=1; i <= NCOL ; i++ ) {  STRFITS.tUnit[i] = STRFITS.blank ;  }
    
  // create table
  fits_create_tbl(fp, BINARY_TBL, NROW, NCOL
		  ,&STRFITS.tName[1]
		  ,&STRFITS.tForm[1]
		  ,&STRFITS.tUnit[1]
		  ,TBLname, &istat );
    
  sprintf(c1err,"fits_create_tbl for %s", TBLname );
  wr_fits_errorCheck(c1err, istat) ;

  // make ordered list of SYN_FILTER info since SYN_FILTERS
  // could be mixed with regular filters

  firstelem = nrow = 1;
  firstrow  = istat = icol = 0;

  for(ifilt=1; ifilt <= NFILTDEF; ifilt++ ) {
    if ( FILTER[ifilt].IFLAG_SYN == 0 ) { continue ; }

    name = FILTER[ifilt].name ;
    firstrow++ ;

    icol=1;
    fits_write_col(fp, TSTRING, icol, firstrow, firstelem, nrow,
		   &name, &istat); 
    sprintf(c1err,"writing SYN FILTER NAME '%s' ", name );
    wr_fits_errorCheck(c1err, istat) ;

    icol=2;  LAM_f = FILTER[ifilt].LAMBDA_MIN ;
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &LAM_f, &istat); 
    sprintf(c1err,"writing SYN_FILTER LAMMIN=%.2f for  %s ", LAM_f, name );
    wr_fits_errorCheck(c1err, istat) ;

    icol=3;  LAM_f = FILTER[ifilt].LAMBDA_MAX ;
    fits_write_col(fp, TFLOAT, icol, firstrow, firstelem, nrow,
		   &LAM_f, &istat); 
    sprintf(c1err,"writing SYN_FILTER LAMMAX=%.2f for  %s ", LAM_f, name );
    wr_fits_errorCheck(c1err, istat) ;
  }

  return ;

} // end wr_fits_SPECTROGRAPH


// ==================================================
void wr_fits_errorCheck(char *comment, int status) {
  // Print out cfitsio error messages and exit program
  char fnam[] = "wr_fits_errorCheck" ;
  if (status) {
    fits_report_error(stderr, status); /* print error report */
    errmsg(SEV_FATAL, 0, fnam, comment, "Check cfitsio routines." ); 
  }
  return;
}  // wr_fits_errorCheck


// ********************************
void dmp_filtersummary(void) {

  // Dump one-line filter summary with ZP and SNMAG at peak

  int ifilt ;
  double snmag0, zp ;
  char *cfilt ;

  // --------------- BEGIN ---------

  printf("\n FILTER SUMMARY \n" );
  printf("            zero     peak   \n" );
  printf("            point    SNmag  \n" );

   for ( ifilt=1;  ifilt<=NFILTDEF; ifilt++ ) {

     cfilt  = FILTER[ifilt].name ;
     snmag0 = SNSED.MAG_RST[IEPOCH_SNPEAK[ifilt]][ifilt] ;
     zp     = FILTER[ifilt].MAGFILTER_ZP;

     printf(" FILTER-%s  %6.3f  %8.3f  \n", cfilt, zp, snmag0 );
   }

} // end of dmp_filtersummary



// **********************************************
void wrsnmag_text(void) {

  /*****
    Created Apr 30, 2009 by R.Kessler
    
    Text dump SN mags in file with hard-wired
    name  SNMAG_[FILTPATH].TXT. Note that there
    is a separate file for each filter path such
    as SDSS, SNLS, etc ...

    File columns are Trest followed by the mag for
    each filter.

   DUMP_SNMAG = 1 => generic table with columns.
                     First row is table header.

   DUMP_SNMAG = 2 => table formatted for non1a preparation
                     (called by sednon1a_install.cmd)

  *****/

  int 
    ifilt,  ifiltmin, ifiltmax, iep, ipath, NPATH
    ,IFILTDEF_PATH[20][2]
    ;

  double snmag,epoch;

  char 
    dumpFile[MXCHAR_FILENAME]
    ,path[40], last_path[40]
    ,fnam[] = "wrsnmag_text" 
    ;

  FILE *fp;

  // ------- BEGIN --------


  if ( INPUTS.DUMP_SNMAG == 0 ) return ;

  printf("\n -------- wrsnmag_text: DUMP SNMAGS  ------------- \n" );

  // find out which filters belong to which path
  sprintf(last_path,"NULL"); NPATH = 0;

  for ( ifilt=1; ifilt <= NFILTDEF; ifilt++ ) {
    sprintf(path, "%s", FILTER[ifilt].PATH ) ;
    if ( strcmp(path,last_path) != 0 ) {
      NPATH++ ;
      IFILTDEF_PATH[NPATH][0] = ifilt;
    }
      IFILTDEF_PATH[NPATH][1] = ifilt;
    sprintf(last_path,"%s", path);
  } // end of ifilt loop


  for ( ipath=1; ipath <= NPATH; ipath++ ) {
    ifiltmin = IFILTDEF_PATH[ipath][0];
    ifiltmax = IFILTDEF_PATH[ipath][1];

    sprintf(dumpFile, "SNMAG_%s.TXT", FILTER[ifiltmin].PATH );

    if ( (fp = fopen(dumpFile, "wt"))==NULL ) { 
      sprintf ( c1err, "Cannot open SNMAG dump-file :" );
      sprintf ( c2err," '%s' ", dumpFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    printf("    Dump SNMAGS to %s (ifiltdef= %d - %d) \n", 
	   dumpFile, ifiltmin, ifiltmax );

    // print table header 
    if ( INPUTS.DUMP_SNMAG == 2 ) {  // for non1a prep

      fprintf(fp, "NEPOCH: %d \n",  SNSED.NEPOCH );
      fprintf(fp, "SNTYPE: %s \n",  INPUTS.SN_TYPE );
     
      for ( ifilt=ifiltmin; ifilt <= ifiltmax ; ifilt++ ) 
	fprintf(fp, "FILTER:  %s  $SNDATA_ROOT/filters/%s \n", 
		FILTER[ifilt].name, FILTER[ifilt].file );
    }
    else {
      fprintf(fp, " EPOCH  " );
      for ( ifilt=ifiltmin; ifilt <= ifiltmax ; ifilt++ ) 
	fprintf(fp, "%8s ", FILTER[ifilt].name );

      fprintf(fp,"\n" );     
    }

    // loop over epochs for dump
    for ( iep=1; iep <= SNSED.NEPOCH ; iep++ ) {
      epoch = SNSED.EPOCH[iep];

      if ( INPUTS.DUMP_SNMAG == 2 )  { fprintf(fp,"EPOCH: "); }

      fprintf(fp,"%9.4f   ", epoch );

      // loop over filters
      for ( ifilt=ifiltmin; ifilt <= ifiltmax ; ifilt++ ) {
	snmag = SNSED.MAG_RST[iep][ifilt] + INPUTS.SN_MAGOFF ;
	fprintf(fp,"%8.3f ", snmag );
      }      
	fprintf(fp, "\n" ); 

    }  // end of iep loop
     
    if ( INPUTS.DUMP_SNMAG == 2 ) { fprintf(fp, "END: \n" ); }


    fclose(fp);

  }
  
} // end of wrsnmag_text

