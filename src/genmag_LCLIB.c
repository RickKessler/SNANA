/************************************************
  Created July 2017 
  Pre-computed Light Curve library, mainly intended
  for stellar-related variables that are not associated
  with redshift.  The LC library contains pre-computed
  magnitudes on a time grid, and is interpolated to   
  a given MJD.

  These functions are used when user select LCLIB model
  with Sim-input  

    GENMODEL:  LCLIB  <lcLibFile>  <stringTemplateEpochs>

 WARNING:
   + set CUTWIN_MJD in sim-input file
   + remove GENRANGE_PEAKMJD and GENRANGE_TREST from sim-input file
      (these are set internally based on CUTWIN_MJD)

  May 21 2018: fix some l,b,GLAT,GLON stuff.

  June 9 2018: 
    + new function set_randomStart_LCLIB is called if sim passes
      IFLAG_RANSTART>0. Ensures full LCLIB sampling for sim-jobs
      split among multiple cores.

  July 24 2018:
    + in function magSearch_LCLIB(), fix bug setting mag for
      non-recurring events when DAY > DAYMAX
      
    + load FIRSTMAG and LASTMAG. For NON-RECUR, check that
      first and last mags are the same.

 Aug 26 2018: 
   + for periodic events, shift library LC with random phase.
     See new function ranPhase_PERIODIC_LCLIB().

 Dec 27 2018: if LCLIB_DEBUG.ZERO_TEMPLATE_FLUX>0, then zero template flux

 Feb 03 2021: if OPTMASK & 8, switch RA,DEC coords to those of LCLIB.
 Jun 25 2021: skip COMMENT lines and DOCANA block reading LCLIB
 Jun 29 2021: call coord_translate_LCLIB before anglematch cut
 Nov 15 2021: fix read_GLOBAL_HEADER_LCLIB bugs from refactor 
                back on Jun 23 2021
*************************************************/

#include "sntools.h"           // community tools
#include "genmag_LCLIB.h" 
#include "MWgaldust.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// ==========================================================
void init_genmag_LCLIB(char *lcLibFile, char *STRING_TEMPLATE_EPOCHS, 
		       char *SURVEY, char *FILTERS, 
		       double TOBS_RANGE_MAX, int NGENTOT, 
		       int IFLAG_RANSTART, int OPTMASK ) {

  // Inputs:
  //   lcLibFile : library file
  //   STRING_TEMPLATE_EPOCHS : e.g., 2+3+4
  //   SURVEY    : name of survey to check what is in LCLIBFILE
  //   FILTERS   : list of bands to check what is in LCLIBFILE
  //   TOBS_RANGE_MAX : max Tobs
  //   NGENTOT  = number of sim events to generate;
  //              used to compute NREPEAT per EVENT
  //
  //  IFLAG_RANSTART: 
  //     logical to start reading at random location in library ...
  //     if NEVENT is specified in LCLIB header. Useful for short
  //     batch jobs to make sure entire library is sampled.
  //
  //  OPTMASK: user bit-mask options passed from sim-input 
  //           GENMODEL_MSKOPT: <OPTMASK>
  //     += 1   --> ignore ANGLEMATCH cut 
  //     += 8   --> use LCLIB coordinates (Feb 2021)
  //     += 512 --> DEGUB/REFACTOR flag
  //
  //  If LCLIBFILE's SURVEY and FILTERLIST does not match input,
  //  abort on error.
  //
  //
  // HISTORY
  // Mar 26 2019: pass OPTMASK arg.
  // Feb    2021: option to use LCLIB coordinates

  char fnam[] = "init_genmag_LCLIB" ;

  // -------------- BEGIN --------------

  print_banner(fnam);

  parse_TEMPLATE_EPOCHS_LCLIB(STRING_TEMPLATE_EPOCHS);

  LCLIB_INFO.TOBS_RANGE_MAX = TOBS_RANGE_MAX;
  LCLIB_INFO.NGENTOT = NGENTOT ;
  LCLIB_INFO.OPTMASK = OPTMASK ;

  LCLIB_INFO.DO_ANGLEMATCH = true;
  if ( OPTMASK & OPTMASK_LCLIB_IGNORE_ANGLEMATCH ) 
    { LCLIB_INFO.DO_ANGLEMATCH = false; }
  if ( OPTMASK & OPTMASK_LCLIB_useRADEC ) 
    { LCLIB_INFO.DO_ANGLEMATCH = false; }

  // -------------------------------
  open_LCLIB(lcLibFile);
  read_GLOBAL_HEADER_LCLIB();

  // -----------------------------------------------------
  // check that LCLIB survey and filters are consistent
  // with SIMLIB/cadence file:
  check_LCLIB(SURVEY,FILTERS);
  fflush(stdout);

  // - - - - - - 

  LCLIB_EVENT.LAST_EXTERNAL_ID  = -999. ;
  LCLIB_EVENT.NROW_MAX    = 0 ;
  LCLIB_EVENT.MEMLC_MAX   = 0 ;
  LCLIB_EVENT.NEVENT_READ = 0 ;
  LCLIB_EVENT.NREPEAT     = LCLIB_INFO.NREPEAT ;
  LCLIB_EVENT.NforceTemplateRows = 0 ;
  LCLIB_EVENT.NROWADD_NONRECUR = 0 ;
  LCLIB_EVENT.NCHAR_ROW = 0;

  printf("\n");

  
  if ( IFLAG_RANSTART ) { set_randomStart_LCLIB(); }

  // -------------------------------
  // check for DEBUG options
  if ( LCLIB_DEBUG.DAYSCALE != 1.0 ) {
    printf("\t DEBUG OPTION: DAYSCALE = %.3f \n",
	   LCLIB_DEBUG.DAYSCALE ) ;
  }
  if ( LCLIB_DEBUG.TOBS_OFFSET_RANGE[1] != 0.0 ) {
    printf("\t DEBUG OPTION: Select random TOBS_OFFSET from "
	   "%.2f - %.2f  days.\n",
	   LCLIB_DEBUG.TOBS_OFFSET_RANGE[0], 
	   LCLIB_DEBUG.TOBS_OFFSET_RANGE[1] );
  }

  if ( (OPTMASK & OPTMASK_LCLIB_IGNORE_ANGLEMATCH)>0 ) {
    printf("\t Ignore ANGLEMATCH \n" );
  }

  if ( (OPTMASK & OPTMASK_LCLIB_useRADEC)>0 ) {
    char line[] = 
      "!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=" ;
    printf("# %s\n", line);
    printf("   WARNING: GENMODEL_MSKOPT = %d includes %d-bit --> \n",
	   OPTMASK, OPTMASK_LCLIB_useRADEC ); 
    printf("\t risky option to overwrite RA,DEC from LCLIB. \n" );
    printf("\t RA,DEC distribution belongs with SIMLIB, not LCLIB. \n" );
    printf("# %s\n", line);
    fflush(stdout);
  }

  //  debugexit(fnam); // xxxx REMOVE

  return ;

} // end  init_genmag_LCLIB


// ===============================
void open_LCLIB(char *lcLibFile) {

  // Sep 30 2020:
  // remove check for ENV PRIVATE_MODELPATH_NAME, and instead
  // check user input PATH_USER_INPUT.

  FILE *fp;
  int  OPTMASK_OPEN = 1;  // 1=verbose
  char MODELPATH_LIST[2*MXPATHLEN];
  char *LCLIB_FILE = LCLIB_INFO.FILENAME ;
  char fnam[] = "open_LCLIB" ;

  // ------------ BEGIN ---------------

  sprintf( MODELPATH_LIST, "%s %s/models/LCLIB", 
	   PATH_USER_INPUT, PATH_SNDATA_ROOT);

  printf("\n");
  fp = snana_openTextFile(OPTMASK_OPEN, MODELPATH_LIST, lcLibFile, 
			  LCLIB_FILE, &LCLIB_INFO.GZIPFLAG ); // <=== returned
  
  if ( fp == NULL ) {
    sprintf(c1err,"Could not open LCLIB file:");
    sprintf(c2err,"%s", LCLIB_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  LCLIB_INFO.FP = fp;  // store file pointer in global

  return;

} // end open_LCLIB

// ====================================================
void read_GLOBAL_HEADER_LCLIB(void) {

  // Jun 2 2018: call parse_PARNAMES_LCLIB
  // Jun 23 2021: refactor to skip comment lines
  // Nov 15 2021: fix refac bugs reading a few doubles (was missing &)

  int NRD_ABORT = 1000 ; // abort after this many words and no FIRST_EVT
  int NRD       = 0 ;
  int FIRST_EVT = 0 ;
  int ipar, NPAR, NFILT, iwd, NWD, NLINE=0; 
  int MSKOPT = MSKOPT_PARSE_WORDS_STRING + MSKOPT_PARSE_WORDS_IGNORECOMMA;
  bool IS_DOCANA = false, IS_COMMENT=false;
  FILE *fp = LCLIB_INFO.FP;
  char wd0[80], wd1[80], wd2[80], LINE[200], tmpString[60], comment[60] ;
  char fnam[] = "read_GLOBAL_HEADER_LCLIB" ;

  // -------------- BEGIN ---------------

  LCLIB_INFO.NFILTERS          = 0 ;

  LCLIB_INFO.IFLAG_RECUR_CLASS       = 0 ;
  LCLIB_INFO.STRING_RECUR_CLASS[0]   = 0 ;

  LCLIB_INFO.SURVEY[0]         = 0 ;
  LCLIB_INFO.FILTERS[0]        = 0 ;
  LCLIB_INFO.NEVENT            = 0 ;
  LCLIB_INFO.NREPEAT           = 1 ;

  LCLIB_INFO.DEBUGFLAG_RANMAG   = 0 ;
  LCLIB_INFO.GENRANGE_RANMAG[0] = 0.0 ;
  LCLIB_INFO.GENRANGE_RANMAG[1] = 0.0 ;
  LCLIB_INFO.GENRANGE_DIFMAG[0] = 0.0 ;
  LCLIB_INFO.GENRANGE_DIFMAG[1] = 0.0 ;

  LCLIB_INFO.IPAR_REDSHIFT = -9;
  LCLIB_INFO.IPAR_MWEBV    = -9;
  LCLIB_INFO.REDSHIFT_RANGE[0] = 0.0 ; 
  LCLIB_INFO.REDSHIFT_RANGE[1] = 0.0 ;

  // June 2 2018: hard wire photo-z error for now;
  // maybe later read this in from LCLIB header
  LCLIB_INFO.ZPHOTZ1ERR = 0.05 ; // error on Zphot/(1+z)

  // rewind to identify and skip documentation block
  snana_rewind(fp, LCLIB_INFO.FILENAME, LCLIB_INFO.GZIPFLAG);

  // - - - - - - - - - - 

  while ( !FIRST_EVT ) { 

    if ( NRD > NRD_ABORT ) {
      sprintf(c1err,"Could not find first event after %d words",NRD);
      sprintf(c2err,"Check global header in LCLIB file." );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
    }

    fgets(LINE, 200, fp ) ;  NLINE++ ;

    if ( commentchar(LINE) ) { continue; }

    NWD = store_PARSE_WORDS(MSKOPT,LINE);
    if ( NWD == 0 ) { continue; }

    NRD += NWD;

    // check first word in line
    iwd=0;  get_PARSE_WORD(0, iwd,  wd0);

    // skip parsing DOCANA and obsolete comment lines
    IS_COMMENT = ( strcmp(wd0,"COMMENT:") == 0 );
    if ( strcmp(wd0,KEYNAME_DOCANA_REQUIRED)  == 0 ) { IS_DOCANA=true;  }
    if ( strcmp(wd0,KEYNAME2_DOCANA_REQUIRED) == 0 ) { IS_DOCANA=false; }
    if ( IS_DOCANA || IS_COMMENT ) { continue; }

    // - - - - - 
    for ( iwd=0; iwd < NWD; iwd++ ) {

      wd0[0] = wd1[0] = 0;
      get_PARSE_WORD(0, iwd,  wd0);
      if ( iwd < NWD-1 ) { get_PARSE_WORD(0, iwd+1, wd1); }

      if ( strcmp(wd0,"NEVENT:") == 0 ) { 
	sscanf(wd1, "%d", &LCLIB_INFO.NEVENT ); iwd++ ;
      }

      if ( strcmp(wd0,"SURVEY:") == 0 ) {
	sscanf(wd1, "%s", LCLIB_INFO.SURVEY ); iwd++ ;
      }

      if ( strcmp(wd0,"FILTERS:") == 0 ) { 
	sscanf(wd1, "%s", LCLIB_INFO.FILTERS ); iwd++ ;
	LCLIB_INFO.NFILTERS = strlen(LCLIB_INFO.FILTERS);
      }

      if ( strcmp(wd0,"RECUR_CLASS:") == 0 )  { 
	sscanf(wd1, "%s", LCLIB_INFO.STRING_RECUR_CLASS ); iwd++ ;
      } 

      if ( strcmp(wd0,"RECUR_TYPE:") == 0 ) { // allow obsolete key
	sscanf(wd1, "%s", LCLIB_INFO.STRING_RECUR_CLASS );  iwd++ ;
      } 
     
      if ( strcmp(wd0,"MODEL:") == 0 )  { 
	sscanf(wd1, "%s", LCLIB_INFO.NAME_MODEL ); iwd++ ;
	if ( strcmp(LCLIB_INFO.NAME_MODEL,MODEL_RANMAG_LCLIB) == 0 ) 
	  { LCLIB_INFO.DEBUGFLAG_RANMAG = 1;  }
	continue ;
      }

      if ( strcmp(wd0,"MODEL_PARNAMES:")   == 0  ||
	   strcmp(wd0,"MODEL_PARAMETERS:") == 0 )  { 
	sscanf(wd1, "%s", tmpString ); iwd++ ;
	parse_PARNAMES_LCLIB(tmpString);
      }

      // read redshift range to pass back to snlc_sim for initializing
      // HOSTLIB.
      if ( strcmp(wd0,"REDSHIFT_RANGE:") == 0 ) { 
	get_PARSE_WORD(0, iwd+2, wd2);
	sscanf(wd1, "%le", &LCLIB_INFO.REDSHIFT_RANGE[0] ); iwd++ ;
	sscanf(wd2, "%le", &LCLIB_INFO.REDSHIFT_RANGE[1] ); iwd++ ;
      } 

      if ( LCLIB_INFO.DEBUGFLAG_RANMAG ) {
	get_PARSE_WORD(0, iwd+2, wd2);
	
	// check DEBUG mag ranges
	if ( strcmp(wd0,"GENRANGE_RANMAG:") == 0 )  {  
	  sscanf(wd1, "%le", &LCLIB_INFO.GENRANGE_RANMAG[0] ); iwd++ ;
	  sscanf(wd2, "%le", &LCLIB_INFO.GENRANGE_RANMAG[1] ); iwd++ ;
	}      
	else if ( strcmp(wd0,"GENRANGE_DIFMAG:") == 0 )  {
	  sscanf(wd1, "%le", &LCLIB_INFO.GENRANGE_DIFMAG[0] ); iwd++ ;
	  sscanf(wd2, "%le", &LCLIB_INFO.GENRANGE_DIFMAG[1] ); iwd++ ;
	}
      }     

      if ( strcmp(wd0,"START_EVENT:") == 0 ) { 
	snana_rewind(fp, LCLIB_INFO.FILENAME, LCLIB_INFO.GZIPFLAG);
	FIRST_EVT = 1 ;  iwd++ ;
      }         

    } // end iwd loop over words in LINE
  } // end while over FIRST_EVT


  // tack on extra model params: EVENT_ID and template mags
  NPAR  = LCLIB_INFO.NPAR_MODEL ;
  NFILT = LCLIB_INFO.NFILTERS ;

  sprintf(LCLIB_INFO.PARNAME_MODEL[NPAR], "EVENT_ID");    NPAR++ ;

  LCLIB_INFO.NPAR_MODEL_STORE = NPAR ;


  // translate RECUR_CLASS string to integer IFLAG_RECUR_CLASS ...
  LCLIB_INFO.IFLAG_RECUR_CLASS = 
    set_IFLAG_RECUR_LCLIB(LCLIB_INFO.STRING_RECUR_CLASS);

  // print model info
  printf("\t MODEL: %s  (%s) \n", 
	 LCLIB_INFO.NAME_MODEL, LCLIB_INFO.STRING_RECUR_CLASS ) ;

  for(ipar=0; ipar < LCLIB_INFO.NPAR_MODEL_STORE ; ipar++ ) {
    comment[0] = 0 ;
    if ( ipar == LCLIB_INFO.IPAR_REDSHIFT ) 
      { sprintf(comment, "<== flagged as redshift"); }

    printf("\t   ModelParam-%2.2d(%s)   %s \n",
	   ipar, LCLIB_INFO.PARNAME_MODEL[ipar], comment );
  }

  fflush(stdout);
  //  debugexit(fnam);

  return ;

} // end read_GLOBAL_HEADER_LCLIB


// ===========================================
void parse_PARNAMES_LCLIB(char *parNameString) {

  // Created Jun 1 2018
  // parse comma-separated list of PARNAMES, and store.
  // Also check for REDSHIFT.

  int ipar, NPAR ;
  char *ptrSplit[MXPAR_LCLIB], *parName ;
  char fnam[] = "parse_PARNAMES_LCLIB" ;

  // -------------- BEGIN --------------

  for(ipar=0; ipar<MXPAR_LCLIB; ipar++ ) {
    ptrSplit[ipar] = LCLIB_INFO.PARNAME_MODEL[ipar]; 
    LCLIB_EVENT.PARVAL_MODEL[ipar] = -999.0 ; 
  }
  splitString(parNameString, ",", MXPAR_LCLIB, &NPAR, ptrSplit);
  LCLIB_INFO.NPAR_MODEL = NPAR ;

  // check for redshift 

  for(ipar=0; ipar < NPAR ; ipar++ ) {
    parName = LCLIB_INFO.PARNAME_MODEL[ipar] ;
    if ( strchr(parName,'#') != NULL ) {
      sprintf(c1err,"Invalid PARNAME='%s' ", parName);
      sprintf(c2err,"hash (#) not allowed");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
    }

    if (strcmp_ignoreCase(parName, PARNAME_REDSHIFT_LCLIB) == 0) 
      { LCLIB_INFO.IPAR_REDSHIFT = ipar;}

    if (strcmp_ignoreCase(parName, PARNAME_MWEBV_LCLIB) == 0) { 
      LCLIB_INFO.IPAR_MWEBV    = ipar;
      sprintf(parName,"%s_LCLIB", PARNAME_MWEBV_LCLIB); // avoid name conflict 
    }
  }
  
  return ;
  
} // end parse_PARNAMES_LCLIB

// ====================================================
int  set_IFLAG_RECUR_LCLIB(char *string_recur_class) {

  int IFLAG = -9 ;
  char fnam[] = "set_IFLAG_RECUR_LCLIB" ;

  // ----------- BEGIN -----------

  if ( strcmp(string_recur_class,"RECUR-PERIODIC") == 0 ) 
    {  IFLAG = IFLAG_RECUR_PERIODIC ;  }
  else if ( strcmp(string_recur_class,"RECUR_PERIODIC") == 0 ) 
    {  IFLAG = IFLAG_RECUR_PERIODIC ;  }
  else if ( strcmp(string_recur_class,"PERIODIC") == 0 )
    {  IFLAG = IFLAG_RECUR_PERIODIC ;  }

  else if ( strcmp(string_recur_class,"RECUR-NONPERIODIC") == 0 ) 
    {  IFLAG = IFLAG_RECUR_NONPERIODIC ;  }
  else if ( strcmp(string_recur_class,"RECUR_NONPERIODIC") == 0 ) 
    {  IFLAG = IFLAG_RECUR_NONPERIODIC ;  }
  else if ( strcmp(string_recur_class,"NONPERIODIC") == 0 )
    {  IFLAG = IFLAG_RECUR_NONPERIODIC ;  }

  else if ( strcmp(string_recur_class,"NON-RECUR") == 0 )
    {  IFLAG = IFLAG_RECUR_NONRECUR ;  }
  else if ( strcmp(string_recur_class,"NON_RECUR") == 0 )
    {  IFLAG = IFLAG_RECUR_NONRECUR ;  }
  else if ( strcmp(string_recur_class,"NONRECUR") == 0 )
    {  IFLAG = IFLAG_RECUR_NONRECUR ;  }

  else {
    sprintf(c1err,"Invalid RECUR_TYPE: '%s'", string_recur_class);
    sprintf(c2err,"Check LCLIB header.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  return(IFLAG) ;

} // end set_IFLAG_RECUR_LCLIB

// ====================================================
void parse_TEMPLATE_EPOCHS_LCLIB(char *STRING_TEMPLATE_EPOCHS) {

  // Break up input STRING_TEMPLATE_EPOCHS into list of floats.
  // Example:
  //   STRING_TEMPLATE_EPOCHS = '1.2+1.3+4+6.3' -->
  //   Template epochs are 1.2  1.3  4.0  6.3 days.
  //   These epochs are relative to the template range
  //   for each event. So if the template range is 10-20 days
  //   for a particular event,
  //   the above TEMPLATE_EPOCHS would be shifted by 10 days.

  int  NEP, ep;
  int  MXEP = MXOBS_TEMPLATE_LCLIB ;
  char tmpString[80];
  char *ptrSplit[MXOBS_TEMPLATE_LCLIB];
  char string_EPLIST[MXOBS_TEMPLATE_LCLIB][12];
  //  char fnam[] = "parse_TEMPLATE_EPOCHS_LCLIB" ;

  // ------------------- BEGIN --------------------


  for(ep=0; ep < MXEP; ep++ ) 
    { ptrSplit[ep] = string_EPLIST[ep] ; }
      
  // copy input string into tmpString since tmpString
  // will be destroyed by splitString.

  sprintf(tmpString, "%s", STRING_TEMPLATE_EPOCHS);
  splitString(tmpString, "+", MXEP, &NEP, ptrSplit);

  LCLIB_INFO.NEP_TEMPLATE = NEP ;

  for(ep=0; ep < NEP; ep++ ) {
    sscanf(string_EPLIST[ep], "%le", &LCLIB_INFO.EPLIST_TEMPLATE[ep] );
    printf("\t Template epoch %2d : %.3f days \n",
	   ep, LCLIB_INFO.EPLIST_TEMPLATE[ep] );
  }

  LCLIB_INFO.EPRANGE_TEMPLATE =
    LCLIB_INFO.EPLIST_TEMPLATE[NEP-1] - 
    LCLIB_INFO.EPLIST_TEMPLATE[0] ;

  printf("\t Template epoch range: %.2f days \n",
	 LCLIB_INFO.EPRANGE_TEMPLATE );

  fflush(stdout);

  return ;


} // end parse_TEMPLATE_EPOCHS_LCLIB

// ==================================================
void check_LCLIB(char *SURVEY, char *FILTERS) {
  
  // Abort if user-defined SURVEY or FILTERS is not compatiable
  // with LCLIB values.  Inputs are passed from sim-input.
  //
  // Also store map:  ifiltmap_LCLIB(ifilt_obs) =ifilt 

  int  ifilt_obs;
  char cfilt[2];
  char fnam[] = "check_LCLIB" ;

  // ------------- BEGIN -----------

  // check for crazy-large TOBS_RANGE_MAX
  double TOBS_RANGE_MAX = LCLIB_INFO.TOBS_RANGE_MAX;
  if ( TOBS_RANGE_MAX > 36500.0 ) {
    sprintf(c1err,"TOBS_RANGE_MAX = %.1f day > 100 years !", TOBS_RANGE_MAX );
    sprintf(c2err,"Check GENRANGE_MJD in sim-input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // make sure that each band in FILTERS (from GENFILTERS sim-input key)
  // is included in the LCLIB.
  // Be careful with ifilt: ifilt_GEN is the sparse index for
  // the user-defined GENFILTERS, while ifilt_LCLIB is the sparse
  // index for the LCLIB.

  int NFILT_GEN = strlen(FILTERS);
  int ifilt_GEN, ifilt_LCLIB, NERR=0;
  char *c ;
  for(ifilt_GEN=0; ifilt_GEN < NFILT_GEN; ifilt_GEN++ ) {
    sprintf(cfilt, "%c", FILTERS[ifilt_GEN] );
    c = strstr(LCLIB_INFO.FILTERS,cfilt) ;
    if ( c == NULL ) {
      NERR++ ;
      printf(" **** LCLIB ERROR *** user band %s is not in LCLIB \n",cfilt);
    }

    // store sparse ifilt in map
    ifilt_obs = INTFILTER(cfilt);  // absolute filter index
    ifilt_LCLIB = (int)( c - LCLIB_INFO.FILTERS );
    ifiltmap_LCLIB[ifilt_obs] = ifilt_LCLIB;

    if ( ifilt_LCLIB < 0 || ifilt_LCLIB >= LCLIB_INFO.NFILTERS ) {
      sprintf(c1err,"Crazy ifilt_LCLIB=%d is outside expected range: %d-%d",
	      ifilt_LCLIB, 0, LCLIB_INFO.NFILTERS-1);
      sprintf(c2err,"cfilt='%s'  GENFILTERS='%s'  LCLIB FILTERS='%s' ",
	      cfilt, FILTERS, LCLIB_INFO.FILTERS );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
    // printf(" xxx ifiltmap_LCLIB[%d] = %d \n", ifilt_obs, ifilt_LCLIB);
  }

  if ( NERR > 0 ) {
    sprintf(c1err,"%d user bands are not in LCLIB", NERR);
    sprintf(c2err,"GENFILTERS='%s' but LCLIB FILTERS='%s' ",
	    FILTERS, LCLIB_INFO.FILTERS );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  return ;

} // end check_LCLIB


// ============================================ 
void set_randomStart_LCLIB(void) {

  // Created Jun 9 2018
  // If NEVENT(header) > NEVENT_MIN, then pick randon event 
  // and move to that event in the library.

  int NEVENT_MIN = 20 ; // at least this many to set random start
  int NEVENT = LCLIB_INFO.NEVENT ;
  
  int    IEVT_START, ievt ;
  double XEVT_START ;
  char   LINE[200];
  char KEY_SEARCH[] = "END_EVENT:" ;
  //  char fnam[] = "set_randomStart_LCLIB" ;

  // -------------- BEGIN ---------------

  if ( NEVENT < NEVENT_MIN ) { return ; }

  XEVT_START = unix_getRan_Flat1(0) * (double)(NEVENT-NEVENT_MIN) ;
  IEVT_START = (int)XEVT_START ;

  printf("\t Skip to random LCLIB event %d of %d ... ",
	 IEVT_START, NEVENT ); fflush(stdout);

  // start reading until reading IEVT_START'th event
  ievt=0;

  while ( ievt < IEVT_START ) {
    fgets(LINE, 40, LCLIB_INFO.FP ) ;
    if ( commentchar(LINE) ) { continue; }
    if ( strstr(LINE,KEY_SEARCH) != NULL )  { ievt++ ; }
  }
  printf("arrived.\n");  fflush(stdout);

  return ;

} // end set_randomStart_LCLIB

// ============================================ 
void set_TobsRange_LCLIB(double *TobsRange) {

  // Created Oct 18 2017
  // main sim routine calls this function to pass global TobsRange,
  // the min and max Tobs over all filters. These global values
  // are used to determine TOBS_OFFSET to map the LCLIB day range
  // into the simulated Tobs range.
  //
  // Input:
  //   TobsRange[0,1] = Tobs_min , Tobs_max
  //
  
  double Tobs_range = TobsRange[1] - TobsRange[0] ;
  char fnam[] = "set_TobsRange_LCLIB" ;

  // ---------------- BEGIN ---------------

  LCLIB_EVENT.TOBS_MIN   = TobsRange[0] ;
  LCLIB_EVENT.TOBS_MAX   = TobsRange[1] ;
  LCLIB_EVENT.TOBS_RANGE = Tobs_range;

  if ( Tobs_range > LCLIB_INFO.TOBS_RANGE_MAX ) {
    print_preAbort_banner(fnam);
    printf("\t TOBS_MIN,MAX = %.1f , %.1f \n",
	   LCLIB_EVENT.TOBS_MIN, LCLIB_EVENT.TOBS_MAX );
    printf("\t TOBS_RANGE = %.1f \n", LCLIB_EVENT.TOBS_RANGE);
    sprintf(c1err,"Tobs_range exceeds expected max");
    sprintf(c2err,"Check SIMLIB and GENRANGE_MJD");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
  }

  /*
  printf(" xxx %s: store TOBS_MIN/MAX=%.2f/%.2f \n",
	 fnam, LCLIB_EVENT.TOBS_MIN, LCLIB_EVENT.TOBS_MAX);
  */

  return ;

} // end set_TobsRange_LCLIB

// ===================================================
void genmag_LCLIB ( int EXTERNAL_ID     // (I) external ID 
		    ,int ifilt_obs       // (I) absolute filter index
		    ,double *RA          // (I,O) RA, deg
		    ,double *DEC         // (I,O) DEC, deg
		    ,double *mwebv       // (I,0) MW E(B-V)
		    ,double lamFilt     // (I) mean filt wave for MWEBV
                    ,int     Nobs       // (I) Number of observations
		    ,double *TobsList   // (I) obs-frame MJD-MJDMIN
		    ,double *magList_S  // (O) obs-frame mag per epoch
		    ,double *mag_T      // (O) template mag
		    ,double *TobsPeak   // (O) day at peak brightness
		    ) 
{

  // Return light curve magList_S at eatch epoch, TobsList.
  // Also return template mag_T.
  // TobsList = 0 corresponds to randomly selected PEAKMJD.
  //
  // July 26 2018: set LAST_EXTERNAL_ID before first return statement.
  // July 30 2018: new output arg, TobsPeak
  // Aug  23 2018: remove obsolete output args ztrue, zphot, zphoterr
  // Aug  29 2018: finally implmenent MWEBV 
  // Feb  03 2021: if OPTMASK & 8, return RA & DEC from LCLIB

  int  IFLAG_NONRECUR = (LCLIB_INFO.IFLAG_RECUR_CLASS==IFLAG_RECUR_NONRECUR);
  LCLIB_EVENT.MWEBV   = *mwebv ;

  double AV_MW, XT_MW;
  int obs, ifilt, NEXT_SIMEVENT, NEXT_LCLIBEVENT ;
  double Tobs, Tobs_shifted, mag_S ;
  //  char fnam[] = "genmag_LCLIB" ;

  // ------------ BEGIN ------------
  
  // init output arguments.
  for(obs=0; obs < Nobs; obs++ ) { magList_S[obs] = MAG_ZEROFLUX ; }
  *mag_T = MAG_ZEROFLUX ;
  *TobsPeak = TobsList[0] - 99.0 ;

  /*
  printf("\n xxx =========== IFILT_OBS=%d  NREAD=%d  ================\n",
	 ifilt_obs, LCLIB_EVENT.NEVENT_READ );
  */

  // if EXTERNAL_ID has changed then read next event;
  // otherwise use same event for different band or different epoch range.
  NEXT_SIMEVENT   = ( LCLIB_EVENT.LAST_EXTERNAL_ID  != EXTERNAL_ID );  
  NEXT_LCLIBEVENT = ( LCLIB_EVENT.NREPEAT == LCLIB_INFO.NREPEAT ||
		      LCLIB_EVENT.NEVENT_READ == 0 ) ;

  if ( NEXT_SIMEVENT && NEXT_LCLIBEVENT  ) {     
    readNext_LCLIB(RA,DEC);  // read next event
    set_NREPEAT_LCLIB(); 
  }

  LCLIB_EVENT.LAST_EXTERNAL_ID  = EXTERNAL_ID ; // July 26 2018
  
  if ( LCLIB_INFO.NREPEAT == 0 ) { return ; }


  if ( NEXT_SIMEVENT ) {  
    LCLIB_EVENT.NREPEAT++ ;

    if ( IFLAG_NONRECUR && LCLIB_EVENT.NROWADD_NONRECUR > 0 ) 
      { addTemplateReset_NONRECUR(); }

    if ( IFLAG_NONRECUR==0 && LCLIB_EVENT.NforceTemplateRows > 0 ) 
      { forceTemplateReset_LCLIB(); }

    set_TOBS_OFFSET_LCLIB();

    addTemplateRows_LCLIB(); // check to add template rows

    dumpEvent_LCLIB();       // user dump
  }


  // get sparse filter index for LCLIB
  ifilt = ifiltmap_LCLIB[ifilt_obs] ;

  // compute galactic extinction
  *mwebv  = LCLIB_EVENT.MWEBV ; // return func arg
  AV_MW   = RV_MWDUST * LCLIB_EVENT.MWEBV ;
  XT_MW   = GALextinct ( RV_MWDUST, AV_MW, lamFilt, 94 );

  // get template mag
  store_magTemplate_LCLIB(EXTERNAL_ID,ifilt,XT_MW);  
  *mag_T  =  LCLIB_EVENT.magTemplate[ifilt] ; // return arg.

  // loop over each epoch
  for(obs=0; obs < Nobs; obs++ ) {
    Tobs         = TobsList[obs] ;
    Tobs_shifted = Tobs + LCLIB_EVENT.TOBS_OFFSET ; // map into LCLIB DAYRANGE
    mag_S        = magSearch_LCLIB(ifilt,Tobs_shifted);
    mag_S       += XT_MW;    // Galactic extinction
    magList_S[obs] = mag_S ;
        
    if ( ifilt_obs == -1 ) {
      char cfilt[2];
      sprintf(cfilt,"%c", LCLIB_INFO.FILTERS[ifilt] );
      printf(" xxx -------------- ID(EXTERN,LCLIB) = %d,%lld ----------- \n",
	     EXTERNAL_ID, LCLIB_EVENT.ID );
      printf(" xxx ifilt=%d(%s)  Tobs[%2d]=%.3f -> %.3f  (TOBS_OFF=%.2f) \n",
	     ifilt, cfilt, obs, Tobs, Tobs_shifted,  
	     LCLIB_EVENT.TOBS_OFFSET);
      printf(" xxx mag[S,T] = %.4f, %.4f  (magDif = %.4f)\n", 
	     mag_S, *mag_T, mag_S - *mag_T );
      fflush(stdout);
      } 
  }


  *TobsPeak = LCLIB_EVENT.PEAKDAY_S - LCLIB_EVENT.TOBS_OFFSET ;

  return ;

}  // end genmag_LCLIB


// ============================================
void readNext_LCLIB(double *RA, double *DEC) {

  // read next LCLIB event.
  // If sky-dependence is part of model then find event 
  // close to input RA,DEC. Otherwise ignore RA,DEC and
  // read next event from the library
  // 
  // Jul 13 2018: MXWD -> += 2 in case NPAR=0
  // Aug 26 2018: call ranPhase_PERIODIC_LCLIB().
  // Feb 03 2021: check OPTMASK&8 option to pass RA,DEC back from HOSTLIB

  double RA_LOCAL  = *RA  ;
  double DEC_LOCAL = *DEC ;

  int  OPTMASK      = LCLIB_INFO.OPTMASK;
  bool switch_RADEC = (OPTMASK & OPTMASK_LCLIB_useRADEC) > 0 ;
  bool REFAC        = (OPTMASK & OPTMASK_LCLIB_DEBUG   ) > 0 ;

  FILE *fp   = LCLIB_INFO.FP ;
  int NPAR   = LCLIB_INFO.NPAR_MODEL ;
  int NFILT  = LCLIB_INFO.NFILTERS;
  int IFLAG_PERIODIC = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_PERIODIC);
  int MXWD   = NFILT + NPAR + 2 ;
  bool FIRST_EVENT = LCLIB_EVENT.NEVENT_READ==0;
  int START_EVENT, END_EVENT, ISROW_T, ISROW_S ;
  int ipar, ifilt, NROW_FOUND, NROW_EXPECT, KEEP, REJECT ;
  int NWD, iwd, NLINE_READ, NLINE_SKIP, Nfread, NCHAR_ROW, istat ;
  long int NCHAR_SKIP;
  bool  VALID_RADEC ;
  char  WDLIST[MXFILTINDX+MXPAR_LCLIB][100], *WD0, *WD1; 
  char *ptrWDLIST[MXFILTINDX+MXPAR_LCLIB];
  double GalLat, GalLong ;
  char key[60], LINE[200], tmpLINE[200];

  char fnam[] = "readNext_LCLIB" ;

  // ------------- BEGIN --------------

  LDUMP_EVENT_LCLIB = 0 ;

  if ( LCLIB_INFO.DEBUGFLAG_RANMAG ) { return; } // nothing to read for debug

  VALID_RADEC = (RA_LOCAL < 900.0 && DEC_LOCAL < 900.0 );
  if ( VALID_RADEC && !switch_RADEC ) 
    { slaEqgal ( RA_LOCAL, DEC_LOCAL, 
		 &GalLong,  &GalLat ); } // return GalLat/Long
  else
    { GalLat = GalLong = 999.0 ; }

  for(iwd=0; iwd < MXWD; iwd++ )  { ptrWDLIST[iwd] = WDLIST[iwd] ; }

 NEXT_EVENT:

  START_EVENT = END_EVENT = 0 ; 
  REJECT = NLINE_READ = NLINE_SKIP = Nfread = 0;
  LCLIB_EVENT.NROW = LCLIB_EVENT.NROW_S = LCLIB_EVENT.NROW_T = 0 ;
  LCLIB_EVENT.RA     = LCLIB_EVENT.DEC     = 999. ;
  LCLIB_EVENT.GLAT   = LCLIB_EVENT.GLON  = 999. ;
  LCLIB_EVENT.DAYRANGE_S[0] = LCLIB_EVENT.DAYRANGE_T[0] = +9.E9 ;
  LCLIB_EVENT.DAYRANGE_S[1] = LCLIB_EVENT.DAYRANGE_T[1] = -9.E9 ;
  LCLIB_EVENT.DAYCOVER_S = LCLIB_EVENT.DAYCOVER_T = 0.0 ;
  LCLIB_EVENT.DAYCOVER_ALL = 0.0 ;
  LCLIB_EVENT.FIRSTROW_T = LCLIB_EVENT.LASTROW_T = 0 ;
  LCLIB_EVENT.FIRSTROW_S = LCLIB_EVENT.LASTROW_S = 99999 ;
  LCLIB_EVENT.PEAKMAG_S  =  MAG_ZEROFLUX ;
  LCLIB_EVENT.PEAKDAY_S  = -99999.0 ;
  LCLIB_EVENT.ANGLEMATCH = LCLIB_EVENT.ANGLEMATCH_b = 99999. ;

  for(ipar=0; ipar < LCLIB_INFO.NPAR_MODEL_STORE; ipar++ ) 
    { LCLIB_EVENT.PARVAL_MODEL[ipar] = -999.0 ; }

  for(ifilt=0; ifilt<NFILT; ifilt++ ) {
    LCLIB_EVENT.FIRSTMAG[ifilt] = 99.0 ;
    LCLIB_EVENT.LASTMAG[ifilt]  = 99.0 ;
  }

  LCLIB_EVENT.REDSHIFT = LCLIB_EVENT.ZPHOT = LCLIB_EVENT.ZPHOTERR = 0.0 ;

  // init local var
  NROW_FOUND = NROW_EXPECT = 0 ;

  //  - - - - - - - - 
  while ( END_EVENT == 0 ) {
    
    LINE[0] = 0 ;
    if ( fgets(LINE, 200, fp ) == NULL ) 
      { snana_rewind(fp, LCLIB_INFO.FILENAME, LCLIB_INFO.GZIPFLAG);  }
    
    if ( commentchar(LINE) ) { continue; } 
    NLINE_READ++;

    // don't parse LINE if we know this event has been rejected
    // and we have not read all of the expected ROWs.
    if ( REJECT  &&  (NLINE_SKIP < NROW_EXPECT-5) )  { 
      NCHAR_ROW = LCLIB_EVENT.NCHAR_ROW;

      if ( REFAC && NCHAR_ROW > 0 ) {
	NCHAR_SKIP  = (NROW_EXPECT-5) * NCHAR_ROW ;
	NLINE_SKIP += (NROW_EXPECT-5) ;
	istat = fseek(fp, NCHAR_SKIP, SEEK_CUR);
      }
      else {
	NLINE_SKIP++ ;
      }
      continue ;
    }

    // split LINE into individual words
    sprintf(tmpLINE, "%s", LINE);
    splitString2(tmpLINE, " ", MXWD, &NWD, ptrWDLIST);
    if ( NWD < 2 ) { continue ; }

    if ( strcmp(WDLIST[0],"START_EVENT:") == 0 ) {
      sscanf(WDLIST[1],"%lld", &LCLIB_EVENT.ID); 
      START_EVENT = 1;
      NLINE_SKIP = REJECT = Nfread = NLINE_READ = 0;
    }       

    if ( REJECT           ) { NLINE_SKIP++ ; continue ; }
    if ( START_EVENT == 0 ) { continue ; }

    iwd=0;
    sprintf(key,"%s", WDLIST[0] );
    ISROW_T = ( strcmp(key,"T:") == 0 );
    ISROW_S = ( strcmp(key,"S:") == 0 );

    while(iwd < NWD-1 ) {

      WD0 = WDLIST[iwd+0];
      WD1 = WDLIST[iwd+1];

      if ( iwd==0 && (ISROW_T || ISROW_S) ) {

	if ( NROW_FOUND==0 ) {
	  if ( LCLIB_EVENT.NEVENT_READ>0 ) { malloc_LCLIB_EVENT(-1); }
	  malloc_LCLIB_EVENT(+1);
	}

	NCHAR_ROW = strlen(LINE);
	LCLIB_EVENT.NCHAR_ROW = NCHAR_ROW;

	if ( NROW_FOUND < NROW_EXPECT ) 
	  { read_ROW_LCLIB(NROW_FOUND,key,ptrWDLIST); }
	NROW_FOUND++ ;	
	iwd += (NFILT+1) ;   // DAY + filters
      }     
      else if ( strcmp(WD0,"NOBS:") == 0  || strcmp(WD0,"NROW:") == 0 ) { 
	sscanf(WD1, "%d", &NROW_EXPECT); iwd++; 
	LCLIB_EVENT.NROW = NROW_EXPECT ;
      }
      
      else if ( strcmp(WD0,"RA:") == 0 ) 
	{
		sscanf(WD1,"%le", &LCLIB_EVENT.RA); iwd++;
        	coord_translate_LCLIB(RA,DEC);
	} 
      else if ( strcmp(WD0,"DEC:") == 0 ) 
	{
		sscanf(WD1,"%le", &LCLIB_EVENT.DEC); iwd++;
  		coord_translate_LCLIB(RA,DEC);
	} 

      else if ( strcmp(WD0,"l:") == 0  || strcmp(WD0,"GLON:")==0 ) 
	{ sscanf(WD1,"%le", &LCLIB_EVENT.GLON ); iwd++; } 
      else if ( strcmp(WD0,"b:") == 0  || strcmp(WD0,"GLAT:")==0  ) 
	{ sscanf(WD1,"%le", &LCLIB_EVENT.GLAT); iwd++ ; } 
  
      else if ( strcmp(WD0,"PARVAL:") == 0 )  { 
	read_PARVAL_LCLIB(LINE);  iwd++ ; 
	START_EVENT = keep_PARVAL_LCLIB() ;
	REJECT = (START_EVENT==0);
      }
    
      else if ( strcmp(WD0,"ANGLEMATCH:") == 0 ) { 
	sprintf(c1err,"ANGLEMATCH not yet implemented.");
	sprintf(c2err,"Try ANGLEMATCH_b");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
      }
      
      else if ( strcmp(WD0,"ANGLEMATCH_b:") == 0 )  { 
	sscanf(WD1,"%le", &LCLIB_EVENT.ANGLEMATCH_b);  iwd++ ;
	LCLIB_INFO.NREPEAT = 1 ; // turn off this feature
	START_EVENT = keep_ANGLEMATCH_LCLIB(GalLat,GalLong);
	REJECT = (START_EVENT==0);
      }
      
      else if ( strcmp(WD0,"END_EVENT:") == 0 ) {
	END_EVENT = 1;  iwd++ ;
      }
      else {
	iwd++ ;
      }
    } // end while over words

  } // end fgets while over lines


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  // if inputs are RA, DEC then translate RA,DEC -> b,l
  coord_translate_LCLIB(RA,DEC);

  // check PARVAL cuts
  KEEP = keep_PARVAL_LCLIB(); 
  if ( KEEP == 0 ) { goto NEXT_EVENT ; }

  // check ANGLEMATCH cut
  KEEP = keep_ANGLEMATCH_LCLIB(GalLat,GalLong);
  if ( KEEP == 0 ) { goto NEXT_EVENT ; }

  // -----------------------------------
  LCLIB_EVENT.DAYCOVER_S = 
    LCLIB_EVENT.DAYRANGE_S[1] - LCLIB_EVENT.DAYRANGE_S[0] ;
  LCLIB_EVENT.DAYCOVER_T = 
    LCLIB_EVENT.DAYRANGE_T[1] - LCLIB_EVENT.DAYRANGE_T[0] ;

  LCLIB_EVENT.DAYCOVER_ALL = 
    LCLIB_EVENT.DAYCOVER_S + LCLIB_EVENT.DAYCOVER_T ;


  LCLIB_EVENT.NEVENT_READ++ ;


  //make sure we find expected number of obs
  if ( NROW_FOUND != NROW_EXPECT ) {
    sprintf(c1err,"NROW_FOUND=%d (T+S rows), but  NROW_EXPECT=%d", 
	    NROW_FOUND, NROW_EXPECT );
    sprintf(c2err,"Check  Event ID=%lld in LCLIB", LCLIB_EVENT.ID);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
  }

  if ( IFLAG_PERIODIC ) { ranPhase_PERIODIC_LCLIB();  }

  return ;

} // end void readNext_LCLIB


// =========================================
void read_ROW_LCLIB(int ROW, char *KEY, char **ptrWDLIST) {

  // Read one ROW from LCLIB.
  // This row has already been identified with
  // with input KEY equal to ether  T: or S:
  //
  // Inputs:
  //   ROW = index to load LCLIB_EVENT.DAY[MAGLIST]
  //   KEY = T: or S:
  //   **ptrWDLIST = words for line (includes KEY)

  int NFILT  = LCLIB_INFO.NFILTERS;
  int ifilt, I2MAG ;
  //  int IFLAG_NONRECUR=(LCLIB_INFO.IFLAG_RECUR_CLASS==IFLAG_RECUR_NONRECUR);
  double DAY, MAG, MAGLIST[MXFILTINDX];
  double *ptrDAYRANGE ;
  char fnam[] = "read_ROW_LCLIB" ;

  // -------------- BEGIN ---------------

  // parse input list of string-words to get DAY and MAGLIST.
  sscanf( ptrWDLIST[1], "%le", &DAY);
  DAY *= LCLIB_DEBUG.DAYSCALE ; // default=1. user sets !=1 for debug

  for(ifilt=0; ifilt < NFILT; ifilt++ )  { 
    sscanf( ptrWDLIST[2+ifilt], "%le", &MAGLIST[ifilt] );  

    if ( ROW==0 ) { LCLIB_EVENT.FIRSTMAG[ifilt] = MAGLIST[ifilt] ; }
    LCLIB_EVENT.LASTMAG[ifilt] = MAGLIST[ifilt] ;
  }


  if ( KEY[0] == 'T' ) { 

    if (LCLIB_EVENT.NROW_T==0 ) { LCLIB_EVENT.FIRSTROW_T=ROW; }
    LCLIB_EVENT.LASTROW_T = ROW;

    LCLIB_EVENT.NROW_T++ ; 
    ptrDAYRANGE = LCLIB_EVENT.DAYRANGE_T ;

    // abort if any T key is after an S key
    if ( ROW > LCLIB_EVENT.FIRSTROW_S ) {
      sprintf(c1err,"Cannot mix T: and S: rows in EVENT ID=%lld "
	      "(FIRSTROW_S=%d)",
	      LCLIB_EVENT.ID, LCLIB_EVENT.FIRSTROW_S );
      sprintf(c2err,"All T: rows must come before S: rows");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
  }
  else  { 
    if (LCLIB_EVENT.NROW_S == 0 ) {  LCLIB_EVENT.FIRSTROW_S = ROW;  }
    LCLIB_EVENT.LASTROW_S = ROW;
    LCLIB_EVENT.NROW_S++ ; 
    ptrDAYRANGE = LCLIB_EVENT.DAYRANGE_S ;
  }


  LCLIB_EVENT.STRING_FLAG[ROW] = KEY[0];
  LCLIB_EVENT.DAY[ROW] = DAY ;

  if ( ROW > 0  &&  DAY < LCLIB_EVENT.DAY[ROW-1] ) {
    print_preAbort_banner(fnam);
    printf("\t ROW = %d \n", ROW);
    printf("\t Previous DAY = %f\n", LCLIB_EVENT.DAY[ROW-1] );
    sprintf(c1err,"Invalid DAY=%f because previous day is larger", DAY);    
    sprintf(c2err,"Check DAY sequence for EVENT ID=%lld", LCLIB_EVENT.ID );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( DAY < ptrDAYRANGE[0] ) { ptrDAYRANGE[0] = DAY; }
  if ( DAY > ptrDAYRANGE[1] ) { ptrDAYRANGE[1] = DAY; }

  for(ifilt=0; ifilt < NFILT; ifilt++   ) { 
    MAG   = MAGLIST[ifilt];

    if ( MAG > 31.0 ) { MAG=31.0; } // avoid I2 overflow (7.15.2018)

    // find epoch of peak Mag (any band) to use as reference
    // for non-recurring transients
    if ( MAG < LCLIB_EVENT.PEAKMAG_S ) 
      { LCLIB_EVENT.PEAKMAG_S = MAG; LCLIB_EVENT.PEAKDAY_S=DAY; }

    // sanity check
    if ( MAG < MAGBRIGHT_LCLIB || MAG>MAGFAINT_LCLIB ) {	    
      sprintf(c1err,"Invalid MAG(%c)=%.3f is outside %.2f to %.2f",
	      LCLIB_INFO.FILTERS[ifilt],
	      MAG, MAGBRIGHT_LCLIB, MAGFAINT_LCLIB);
      sprintf(c2err,"Check  Event ID=%lld, DAY=%f", 
	      LCLIB_EVENT.ID, DAY);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
    }

    // store MAG*1000. as short int (to save memory)
    I2MAG = (int)(I2FLOAT_LCLIB * MAG + 0.5) ; 
    LCLIB_EVENT.I2MAG[ifilt][ROW] = I2MAG ;

  } // end ifilt loop

 
  return;

} // end read_ROW_LCLIB


// =========================================
void read_PARVAL_LCLIB(char *LINE) {

  // "PARVAL:" key has already been read, so now 
  // read list of parameter values after "PARVAL:" key.
  // List can be either commma-separated of space-separated,
  // so it's a little tricky.
  //
  // Input LINE is entire PARVAL line including PARVAL key.

  int NPAR   = LCLIB_INFO.NPAR_MODEL ;
  int MXSPLIT = MXPAR_LCLIB ;
  int ipar, NSPLIT ;
  double PARVAL;
  char  tmpLine[200], sepKey[4] = " " ;
  char *ptrPARVAL[MXPAR_LCLIB], stringPARVAL[MXPAR_LCLIB][20];
  char fnam[] = "read_PARVAL_LCLIB" ;

  // ------------ BEGIN -------------

  // create tmpLine = LINE so that tmpLine can be destroyed
  sprintf(tmpLine, "%s", &LINE[8]); // avoid PARVAL key
  if ( strstr(tmpLine,",") != NULL )  { sepKey[0] = ',' ; }

  for(ipar=0; ipar < MXSPLIT; ipar++ ) 
    { ptrPARVAL[ipar] = stringPARVAL[ipar] ; }

  splitString2(tmpLine, sepKey, MXSPLIT, &NSPLIT, ptrPARVAL);
  if ( NSPLIT != NPAR ) {
    print_preAbort_banner(fnam);
    //    printf("\t PARVAL line = '%s' \n", tmpLine);
    printf("\t sepKey = '%s' \n", sepKey );
    for(ipar=0; ipar<NSPLIT; ipar++ ) {
      printf("\t PARVAL[%2d] = '%s'\n", ipar, stringPARVAL[ipar]); 
    }

    sprintf(c1err,"Found %d PARVAL values for ID=%lld",
	    NSPLIT, LCLIB_EVENT.ID );
    sprintf(c2err,"but expected NPAR=%d values from header.", NPAR);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
  }

  for(ipar=0; ipar < NPAR; ipar++ ) {
    sscanf(stringPARVAL[ipar], "%le", &PARVAL );
    LCLIB_EVENT.PARVAL_MODEL[ipar] = PARVAL ;
  }
  LCLIB_EVENT.PARVAL_MODEL[NPAR+0] = (double)LCLIB_EVENT.ID ;  

  // check for optional redshift
  set_REDSHIFT_LCLIB();

  return ;

} // end read_PARVAL_LCLIB

// =========================================
void coord_translate_LCLIB(double *RA, double *DEC) {

  // Created Jun 29 2021
  // wrapper to convert RA,DEC -> b,l

  bool switch_RADEC = (LCLIB_INFO.OPTMASK & OPTMASK_LCLIB_useRADEC) > 0 ;
  char fnam[] = "coord_translate_LCLIB";

  if ( LCLIB_EVENT.RA < 900 && LCLIB_EVENT.DEC < 900.0 ) {
    // convert J2000 into Galactic coords
    slaEqgal ( LCLIB_EVENT.RA, LCLIB_EVENT.DEC, 
	       &LCLIB_EVENT.GLON,  &LCLIB_EVENT.GLAT );  // returned

    if ( switch_RADEC ) {
      // update Galactic extinction for LCLIB coords
      *RA = LCLIB_EVENT.RA;   *DEC=LCLIB_EVENT.DEC;
      double MWEBV = gen_MWEBV(*RA,*DEC);
      LCLIB_EVENT.MWEBV = MWEBV;
    }
  }
  else if ( LCLIB_EVENT.GLAT < 900 && LCLIB_EVENT.GLON < 900 ) {
    // Convert Galactic coords into J2000
    // crap, don't have translator function in SNANA ?!?!?!?
    LCLIB_EVENT.RA = LCLIB_EVENT.DEC = 66.666 ;
  }

  return;

} // end coord_translate_LCLIB

// =========================================
int keep_PARVAL_LCLIB(void) {

  // Oct 31 2017
  // check optional cuts on PARVAL
  // Returns 1 to keep, 0 to reject.

  int KEEP = 1;
  int NCUTWIN = LCLIB_CUTS.NCUTWIN ;
  int NPAR    = LCLIB_INFO.NPAR_MODEL ;

  int ipar, icut ;
  char *PARNAME, *CUTNAME ;
  double dval, DMIN, DMAX;

  // ------------- BEGIN ---------------
  if ( NCUTWIN == 0 ) { return(KEEP); }

  for(ipar=0; ipar < NPAR; ipar++ ) {
    PARNAME = LCLIB_INFO.PARNAME_MODEL[ipar];
    for(icut=0; icut<NCUTWIN; icut++ ) {
      CUTNAME = LCLIB_CUTS.PARNAME[icut];
      if ( strcmp(PARNAME,CUTNAME) == 0 ) {
	dval = LCLIB_EVENT.PARVAL_MODEL[ipar];
	DMIN = LCLIB_CUTS.CUTWIN[icut][0] ;
	DMAX = LCLIB_CUTS.CUTWIN[icut][1] ;
	if ( dval < DMIN ) { KEEP = 0 ; }
	if ( dval > DMAX ) { KEEP = 0 ; }

	/*
	if ( KEEP == 0 ) {
	  printf(" xxx REJECT EVENT_ID=%8d for %s = %f  (range=%.1f to %.1f)\n",
		 LCLIB_EVENT.ID, PARNAME, dval, DMIN, DMAX);
	} 
	*/

      }
    }
  }


  return(KEEP) ;

} // end keep_PARVAL_LCLIB

// ============================================
int keep_ANGLEMATCH_LCLIB(double b, double l) {

  // Dec 29 2017
  // Apply optional ANGLEMATCH_b  cut
  // ANGLEMATCH (radius) cut not defined.
  // Input b,l are the Galactic coords of the simulated event
  // Galactic coords of LCLIB event are stored in 
  // LCLIB_EVENT.GLAT[GLON].
  
  int KEEP=1 ;
  double b_SIM    = fabs(b);
  double b_LCLIB  = fabs(LCLIB_EVENT.GLAT);
  char fnam[] = "keep_ANGLEMATCH_LCLIB" ;

  // ------------- BEGIN ------------

  // check option to skip ANGLEMATCH cut
  if ( !LCLIB_INFO.DO_ANGLEMATCH ) { return(KEEP); }

  // bail if not cut is defined.
  if ( LCLIB_EVENT.ANGLEMATCH_b > 500.0 ) { return(KEEP); }

  // apply cut
  if ( fabs(b_SIM - b_LCLIB ) > LCLIB_EVENT.ANGLEMATCH_b ) 
    { KEEP = 0 ; }

  /*
  printf(" xxx b(SIM,LCLIB) = %7.1f, %7.1f   KEEP=%d \n",
	 b, LCLIB_EVENT.GLAT, KEEP);
  */

  return(KEEP);

} // end keep_ANGLEMATCH_LCLIB


// ====================================
void ranPhase_PERIODIC_LCLIB(void) {

  // Aug 26 2018
  // Called for periodic events, apply random phase shift to light curve.
  // This prevents phase-locking Search & template epochs.
  // 

#define I2MAG_UNFILLED 6

  int    NFILT    = LCLIB_INFO.NFILTERS;
  int    NROW_S   = LCLIB_EVENT.NROW_S ;
  double PERIOD   = LCLIB_EVENT.DAYCOVER_S ;
  double RANPHASE = PERIOD * getRan_Flat1(1) ;
  int    ID       = LCLIB_EVENT.ID;

  double DT ;
  int  NERR, ifilt, row, row_new, row_old, IROW_START = -9 ;
  int   *ROW_ORIG ;
  short int **I2MAGTMP, I2MAG, I2MAG_NEXT ;
  int  MEMI2 = sizeof(short int);
  int  DMPROW, DMPROW_LAST, LDMP = 0 ;
  char star[4], CLINE[100], CLINE_LAST[100] ;
  //  char fnam[] = "ranPhase_PERIODIC_LCLIB" ;

  // ----------------- BEGIN ----------------

  // allocate local memory to hold original event mags
  I2MAGTMP  = (short int**) malloc ( NFILT*sizeof(short int*) );
  for(ifilt=0; ifilt < NFILT; ifilt++ ) 
    { I2MAGTMP[ifilt]  = (short int*) malloc(NROW_S*MEMI2) ;  }
  ROW_ORIG = (int*) malloc( NROW_S * sizeof(int) ) ;


  // store original I2MAGs, and find epoch of random phase
  for(row=0; row < NROW_S; row++ ) {    
    for(ifilt=0; ifilt < NFILT; ifilt++ ) { 
      I2MAGTMP[ifilt][row] = LCLIB_EVENT.I2MAG[ifilt][row] ; 
      LCLIB_EVENT.I2MAG[ifilt][row] = I2MAG_UNFILLED ;
    }   
    DT = LCLIB_EVENT.DAY[row] - LCLIB_EVENT.DAY[0] ;
    if ( DT > RANPHASE && IROW_START < 0 ) { IROW_START = row;  }
  }
  if ( IROW_START >= NROW_S ) { IROW_START = NROW_S-1; }
 

  // re-load I2MAGs with random phase offset.
  // It's a little tricky because last original epoch
  // is a duplicate of first epoch, so exclude this.
  int NROW_copy = NROW_S - 1 ;
  for(row_new=0; row_new < NROW_copy; row_new++ ) {
    
    row_old = row_new + IROW_START;
    if ( row_old >= NROW_copy ) { row_old -= NROW_copy; }

    ROW_ORIG[row_new] = row_old;
    for(ifilt=0; ifilt < NFILT; ifilt++ ) { 
      LCLIB_EVENT.I2MAG[ifilt][row_new] = I2MAGTMP[ifilt][row_old] ;
    }
  }

  // make sure last epoch is same as first epoch
  for(ifilt=0; ifilt < NFILT; ifilt++ ) 
    { LCLIB_EVENT.I2MAG[ifilt][NROW_S-1] = LCLIB_EVENT.I2MAG[ifilt][0] ; }

  
  // - - - - - -
  if ( LDMP ) {
    NERR = DMPROW_LAST = DMPROW = 0 ;
    printf(" xxx ----------------------- \n");
    printf(" xxx "
	   "ID=%d  NROW_S=%4d  PERIOD=%8.2f   RanPhase=%8.2f  IROW_START=%d\n",
	   ID, NROW_S, PERIOD, RANPHASE, IROW_START );
    for(row_new=0; row_new < NROW_S; row_new++ ) {

      for(ifilt=0; ifilt < NFILT; ifilt++ ) {
	I2MAG = LCLIB_EVENT.I2MAG[ifilt][row_new] ;
	if ( I2MAG == I2MAG_UNFILLED ) { NERR++ ; }
      }

      star[0] = DMPROW = 0;
      sprintf(CLINE," xxx DAY=%.2f I2MAG=", LCLIB_EVENT.DAY[row_new] );
      for(ifilt=0; ifilt < NFILT; ifilt++ ) {
	I2MAG      = LCLIB_EVENT.I2MAG[ifilt][row_new] ;
	I2MAG_NEXT = LCLIB_EVENT.I2MAG[ifilt][row_new+1] ;
	sprintf(CLINE,"%s %5d", CLINE, I2MAG ) ;
	if ( abs(I2MAG-I2MAG_NEXT) > 200 && row_new<NROW_S-1) 
	  { sprintf(star,"*"); DMPROW=1; }
      }
      sprintf(CLINE, "%s %s (ROW_ORIG=%d)\n", 
	      CLINE, star, ROW_ORIG[row_new] ); fflush(stdout);

      if ( DMPROW_LAST ) {
	printf("%s\n", CLINE_LAST );
	printf("%s\n", CLINE );
      }
      
      DMPROW_LAST=DMPROW; sprintf(CLINE_LAST,"%s", CLINE);
    }

    
    printf(" xxx NERR(unfilled mag): %d \n", NERR);
    
  } // end LDMP

  
  // free I2MAGTMP
  for(ifilt=0; ifilt < NFILT; ifilt++ ) { free(I2MAGTMP[ifilt]); }
  free(I2MAGTMP);
  free(ROW_ORIG);

  return;

} // end ranPhase_PERIODIC_LCLIB

// ====================================
void  set_REDSHIFT_LCLIB(void) {

  // Created Jun 2 2018
  // if REDSHIFT is a parameter, this set global redshift
  // and photo-z info.

  int IPAR = LCLIB_INFO.IPAR_REDSHIFT ;
  double ZTRUE;

  // ------------ BEGIN -----------

  if ( IPAR < 0 ) { return ; }

  ZTRUE = LCLIB_EVENT.PARVAL_MODEL[IPAR] ;
  LCLIB_EVENT.REDSHIFT   = ZTRUE ;

 
  return ;
  
} // end set_REDSHIFT_LCLIB

// =========================================
void set_NREPEAT_LCLIB(void) {

  // set LCLIB_INFO.NREPEAT.
  // For recurring events, do nothing (i.e., leave NREPEAT to its
  // initialized value based on NEVENT)/
  // For non-recurring events, use event duration and survey duration
  // to set NREPEAT. For example, if event duration is 100 years
  // and survey duration is 10 years, then <NREPEAT>=10 for a
  // Poisson generator to pick NREPEAT.

  int  IFLAG_NONRECUR =(LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR);
  double XREPEAT ;
  int  LDMP ;
  char fnam[] = "set_NREPEAT_LCLIB"; 

  // --------------- BEGIN -------------

  if ( IFLAG_NONRECUR ) {

    LCLIB_EVENT.NROWADD_NONRECUR = 0 ;

    // total TGENRAGE depends on max possible TOBS_RANGE_MAX
    // because we haven't yet read each individual TOBS_RANGE.
    LCLIB_EVENT.TGENRANGE_NONRECUR = 
      LCLIB_EVENT.DAYCOVER_ALL + LCLIB_INFO.TOBS_RANGE_MAX ;


    // compute Poisson mean
    XREPEAT = LCLIB_EVENT.TGENRANGE_NONRECUR/LCLIB_EVENT.TOBS_RANGE;

    // compute random Poisson integer for number of LCLIB-EVENT repeats
    if ( LCLIB_EVENT.NEVENT_READ==1 ) { getRan_Poisson(0.0); }
    LCLIB_INFO.NREPEAT = getRan_Poisson(XREPEAT);

    if ( LCLIB_DEBUG.FORCE_NREPEAT > 0 ) 
      { LCLIB_INFO.NREPEAT = LCLIB_DEBUG.FORCE_NREPEAT; }

    LDMP = 0 ;
    if ( LDMP ) {
      printf(" xxx --------------------------------- \n");
      printf(" xxx %s DUMP  \n", fnam);
      printf(" xxx RANGE[survey,LCLIB,ALL] = %.1f, %.1f, %.1f  \n",
	     LCLIB_INFO.TOBS_RANGE_MAX, 
	     LCLIB_EVENT.DAYCOVER_ALL, 
	     LCLIB_EVENT.TGENRANGE_NONRECUR);     
      printf(" xxx XREPEAT=%.2f -> NREPEAT=%d\n", 
	     XREPEAT, LCLIB_INFO.NREPEAT );

      //      printf(" 77777  %d\n", LCLIB_INFO.NREPEAT );
      fflush(stdout);
    }
  }

  // always reset repeat-counter
  LCLIB_EVENT.NREPEAT = 0 ;

  return;

} // end set_NREPEAT_LCLIB


// =========================================
void addTemplateRows_LCLIB(void) {

  // Add template rows based on location of event within survey window.

  int NFILT  = LCLIB_INFO.NFILTERS;
  int PERIODIC = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_PERIODIC );
  int NONRECUR = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR );
  int ifilt ;
  char fnam[] = "addTemplateRows_LCLIB" ;

  // ------------- BEGIN ----------------

  /*
  printf(" xxx -------------------------- \n");
  printf(" 1. xxx %s: Nforce=%d \n",
	 fnam, LCLIB_EVENT.NforceTemplateRows ); fflush(stdout);
  */


  LCLIB_EVENT.NforceTemplateRows = 0;


  if ( PERIODIC ) {

    // For PERIODIC events, check to make templates.
    // First, make sure that first and last mag are equal  
    for(ifilt=0; ifilt < NFILT; ifilt++ ) 
      { checkMag_PERIODIC_LCLIB(ifilt);  } // end ifilt

    addTemplateRows_PERIODIC(); // add templates before search
  }
  else if ( NONRECUR ||  (LCLIB_EVENT.NROW_T == 0) ) {

    if ( NONRECUR ) { addTemplateRows_NONRECUR(); }

    // no template rows for non-periodoc events 
    // --> convert search rows into template rows.
    forceTemplateRows_LCLIB();

  }

  // for non-peridic events, must have explicit template rows
  if ( LCLIB_EVENT.NROW_T == 0 && PERIODIC==0 ) {
    sprintf(c1err,"No template (T:) rows found in LCLIB EVENT ID=%lld",
	    LCLIB_EVENT.ID);
    sprintf(c2err,"Must have at least 1 template row");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
  }


  return ;

} // end addTemplateRows_LCLIB


// ===============================================
void forceTemplateReset_LCLIB(void) {

  // reset forced template rows back into Search rows
  int FIRSTROW_T  = LCLIB_EVENT.FIRSTROW_T ;
  int LASTROW_T   = LCLIB_EVENT.LASTROW_T ;
  int row;
  //  char fnam[] = "forceTemplateReset_LCLIB" ;

  // ------------------- BEGIN -----------------

  for(row=FIRSTROW_T; row <= LASTROW_T; row++ ) 
    { LCLIB_EVENT.STRING_FLAG[row] = 'S' ; }

  LCLIB_EVENT.NROW_T      =  0 ;
  LCLIB_EVENT.FIRSTROW_T  =  0 ;
  LCLIB_EVENT.LASTROW_T   =  0 ;
  LCLIB_EVENT.FIRSTROW_S  =  0 ;

  return ;

} // end forceTemplateReset_LCLIB


// ===============================================
void forceTemplateRows_LCLIB(void) {

  // If Search time-window covers Tobs range, and there
  // are no template rows, convert some Search rows into
  // template rows. This is designed for non-recurring transients
  // that start before the survey and thus the true template mag
  // is not observed. Can also be used for recurring events
  // for which templates are not defined.
  // 

  char fnam[] = "forceTemplateRows_LCLIB" ;

  // -------------- BEGIN --------------

  if ( LCLIB_EVENT.NROW_T > 0 ) { return ; }

  // --------------------------------------
  // --------------------------------------

  int FIRSTROW_S     = LCLIB_EVENT.FIRSTROW_S ;
  int NEW_FIRSTROW_S = -1 ;
  int NEW_LASTROW_S  = LCLIB_EVENT.LASTROW_S ; // doesn't change
  int NEW_FIRSTROW_T = -1 ;
  int NEW_LASTROW_T  = -1 ;
  int row, NREAD, NBACK_ITER=0;
  int ERRFLAG=0, LDMP=0 ;
  double DAYBACK_T = DAYBACK_TEMPLATE_LCLIB ; 
  double dayback_T ;
  double DAYMIN_T, DAY, DAYRANGE_ACTUAL ;

  // start by finding last Search DAY that occurs
  // before LCLIB_EVENT.DAY_RANDOM.

  //  for(row=FIRSTROW_S; row < NROW; row++ ) {
  for(row=FIRSTROW_S; row < NEW_LASTROW_S; row++ ) {
    if ( LCLIB_EVENT.DAY[row] < LCLIB_EVENT.DAY_RANDOM ) 
      { NEW_FIRSTROW_S = row; }
  }


  // now go back in time to cover template time range that occurs
  // just before Search range. So if library spans 1000 years,
  // the template and search range should be close in time,
  // not separated by hundreds of years.
  
  /*
  printf(" 7. xxx %s: NROW=%d  NROW_S=%d \n",
	   fnam, LCLIB_EVENT.NROW, LCLIB_EVENT.NROW_S);
  */

  dayback_T = 0.0 ;    NEW_LASTROW_T = NEW_FIRSTROW_S ;
  NREAD = LCLIB_EVENT.NEVENT_READ ;

  while ( dayback_T < DAYBACK_T ) {
    NEW_LASTROW_T-- ;  NBACK_ITER++ ;

    if ( NEW_LASTROW_T < 1 ) {
      print_preAbort_banner(fnam);
      printf("\t NREAD=%d \n", NREAD );
      printf("\t NEW_FIRSTROW_S=%d --> NEW First DAY_S=%.3f \n",
	     NEW_FIRSTROW_S, LCLIB_EVENT.DAY[NEW_FIRSTROW_S] );
      printf("\t DAY_RANDOM = %.2f \n", LCLIB_EVENT.DAY_RANDOM);
      printf("\t TOBS_MIN/MAX = %.2f/%.2f (RANGE=%.2f)\n",
	     LCLIB_EVENT.TOBS_MIN, LCLIB_EVENT.TOBS_MAX, 
	     LCLIB_EVENT.TOBS_RANGE);
      printf("\t TOBS_OFFSET = %.2f \n", LCLIB_EVENT.TOBS_OFFSET);
      printf("\t DAYCOVER_ALL = %.2f = %.2f(S) + %.2f(T) \n",
	     LCLIB_EVENT.DAYCOVER_ALL, 
	     LCLIB_EVENT.DAYCOVER_S, LCLIB_EVENT.DAYCOVER_T) ;
      printf("\t DAYRANGE_S = %.2f, %.2f \n",
	     LCLIB_EVENT.DAYRANGE_S[0], LCLIB_EVENT.DAYRANGE_S[1] );

      sprintf(c1err,"Could not find NEW_LASTROW_T after %d iterations",
	      NBACK_ITER);
      sprintf(c2err,"Check EVENT_ID=%lld", LCLIB_EVENT.ID);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
    }
    dayback_T = 
      LCLIB_EVENT.DAY[NEW_FIRSTROW_S] - LCLIB_EVENT.DAY[NEW_LASTROW_T] ;
  }


  DAYMIN_T = LCLIB_EVENT.DAY[NEW_LASTROW_T] - LCLIB_INFO.EPRANGE_TEMPLATE;
  for(row=FIRSTROW_S; row < NEW_LASTROW_T; row++ ) {
    if ( LCLIB_EVENT.DAY[row] < DAYMIN_T ) 
      { NEW_FIRSTROW_T = row; } 
  }

  if ( NEW_FIRSTROW_S < 0 ) { ERRFLAG=1; }
  if ( NEW_FIRSTROW_T < 0 ) { ERRFLAG=2; }
  if ( NEW_LASTROW_T  < 0 ) { ERRFLAG=3; }
  if ( NEW_LASTROW_S < NEW_FIRSTROW_S ) { ERRFLAG=4; }

  if ( LDMP || ERRFLAG ) {
    printf("\n xxx --------- %s DUMP for ID=%lld --------------- \n",
	   fnam, LCLIB_EVENT.ID );

    printf(" xxx TOBS_OFFSET = %.3f\n", LCLIB_EVENT.TOBS_OFFSET );

    printf(" xxx DAYMIN_T = %.3f - %.3f = %.3f \n",
	   LCLIB_EVENT.DAY[NEW_LASTROW_T], 
	   LCLIB_INFO.EPRANGE_TEMPLATE, DAYMIN_T );

    DAY = LCLIB_EVENT.DAY[NEW_FIRSTROW_T];
    printf(" xxx FIRSTROW_T = %d -> %d  (DAY=%.3f)\n",
	   LCLIB_EVENT.FIRSTROW_T, NEW_FIRSTROW_T, DAY  );

    DAY = LCLIB_EVENT.DAY[NEW_LASTROW_T];
    printf(" xxx LASTROW_T  = %d -> %d  (DAY=%.3f\n",
	   LCLIB_EVENT.LASTROW_T, NEW_LASTROW_T, DAY ) ;

    DAY = LCLIB_EVENT.DAY[NEW_FIRSTROW_S];
    printf(" xxx FIRSTROW_S = %d -> %d  (DAY=%.3f) \n",
	   LCLIB_EVENT.FIRSTROW_S, NEW_FIRSTROW_S, DAY );

    DAY = LCLIB_EVENT.DAY[NEW_LASTROW_S];
    printf(" xxx LASTROW_S  = %d -> %d  (DAY=%.3f) \n",
	   LCLIB_EVENT.LASTROW_S, NEW_LASTROW_S, DAY );

    DAYRANGE_ACTUAL = 
      LCLIB_EVENT.DAY[NEW_LASTROW_T] - 
      LCLIB_EVENT.DAY[NEW_FIRSTROW_T] ;
    printf(" xxx Template time-range (Request,Actual): %.1f , %.1f \n",
	   LCLIB_INFO.EPRANGE_TEMPLATE, DAYRANGE_ACTUAL);

  }

  if (ERRFLAG ) {
    sprintf(c1err,"Could not force template rows for EVENT_ID=%lld",
	    LCLIB_EVENT.ID);
    sprintf(c2err,"ERRFLAG=%d; see above", ERRFLAG);    
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
  }

  //  debugexit(fnam);

  LCLIB_EVENT.NROW_S     = NEW_LASTROW_S - NEW_FIRSTROW_S + 1;
  LCLIB_EVENT.NROW_T     = NEW_LASTROW_T - NEW_FIRSTROW_T + 1;
  LCLIB_EVENT.NforceTemplateRows = LCLIB_EVENT.NROW_T;

  LCLIB_EVENT.FIRSTROW_T = NEW_FIRSTROW_T ;
  LCLIB_EVENT.LASTROW_T  = NEW_LASTROW_T ;
  LCLIB_EVENT.FIRSTROW_S = NEW_FIRSTROW_S ;
  LCLIB_EVENT.LASTROW_S  = NEW_LASTROW_S ;
  for(row = NEW_FIRSTROW_T; row <= NEW_LASTROW_T; row++ ) 
    { LCLIB_EVENT.STRING_FLAG[row] = 'T' ; }

  LCLIB_EVENT.DAYRANGE_T[0] = LCLIB_EVENT.DAY[NEW_FIRSTROW_T];
  LCLIB_EVENT.DAYRANGE_T[1] = LCLIB_EVENT.DAY[NEW_LASTROW_T];
  LCLIB_EVENT.DAYRANGE_S[0] = LCLIB_EVENT.DAY[NEW_FIRSTROW_S];
  LCLIB_EVENT.DAYRANGE_S[1] = LCLIB_EVENT.DAY[NEW_LASTROW_S];

  LCLIB_EVENT.DAYCOVER_T = 
    LCLIB_EVENT.DAYRANGE_T[1] - LCLIB_EVENT.DAYRANGE_T[0];
  LCLIB_EVENT.DAYCOVER_S = 
    LCLIB_EVENT.DAYRANGE_S[1] - LCLIB_EVENT.DAYRANGE_S[0];
  LCLIB_EVENT.DAYCOVER_ALL = 
    LCLIB_EVENT.DAYCOVER_S + LCLIB_EVENT.DAYCOVER_T ;


  return ;

} // end forceTemplateRows_LCLIB


// ===============================================
void checkMag_PERIODIC_LCLIB(int ifilt) {

  // ABORT if first and last mag are not the same.

  int  NROW       = LCLIB_EVENT.NROW ;
  int  FIRSTROW_S = LCLIB_EVENT.FIRSTROW_S ;
  int  IMAG_FIRST, IMAG_LAST ;
  float MAG_FIRST, MAG_LAST ;
  char fnam[] = "checkMag_PERIODIC_LCLIB";

  IMAG_FIRST = LCLIB_EVENT.I2MAG[ifilt][FIRSTROW_S] ;
  IMAG_LAST  = LCLIB_EVENT.I2MAG[ifilt][NROW-1] ;
  MAG_FIRST  = (float)IMAG_FIRST/I2FLOAT_LCLIB ;
  MAG_LAST   = (float)IMAG_LAST /I2FLOAT_LCLIB ;
  if ( fabsf(MAG_FIRST-MAG_LAST) > 0.0011 ) {
    sprintf(c1err,"Invalid PERIODIC Event (ID=%lld, filt=%c)", 
	    LCLIB_EVENT.ID, LCLIB_INFO.FILTERS[ifilt] );
    sprintf(c2err,"First/last MAG are different: %.4f/%.4f",
	    MAG_FIRST, MAG_LAST);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );      
  }

  return ;

} // end checkMag_PERIODIC_LCLIB

// ===============================================
void dumpEvent_LCLIB(void) {

  int NFILT  = LCLIB_INFO.NFILTERS;
  int ifilt, row ;
  char fnam[] = "dumpEvent_LCLIB" ;

  // ------------ BEGIN ---------------

  LDUMP_EVENT_LCLIB = ( LCLIB_EVENT.ID == -1 )  ;
  if  ( LDUMP_EVENT_LCLIB ) {
    printf("\n");
    printf(" xxx -------- %s DUMP for EVENT ID=%lld --------------- \n", 
	   fnam, LCLIB_EVENT.ID );
    printf(" xxx RA=%f    DEC=%f   FILTERS=%s \n", 
	   LCLIB_EVENT.RA, LCLIB_EVENT.DEC, LCLIB_INFO.FILTERS);
    printf(" xxx NROW_T=%5d   FIRSTROW_T=%5d   LASTROW_T=%5d \n",
	   LCLIB_EVENT.NROW_T, LCLIB_EVENT.FIRSTROW_T, LCLIB_EVENT.LASTROW_T );
    printf(" xxx NROW_S=%5d   FIRSTROW_S=%5d   LASTROW_S=%5d \n",
	   LCLIB_EVENT.NROW_S, LCLIB_EVENT.FIRSTROW_S, LCLIB_EVENT.LASTROW_S );

    printf(" xxx TOBS_OFFSET = %.3f   DAY_RANDOM=%.3f\n", 
	   LCLIB_EVENT.TOBS_OFFSET, LCLIB_EVENT.DAY_RANDOM );

    for(row=LCLIB_EVENT.FIRSTROW_T; row < LCLIB_EVENT.NROW; row++ ) {
      printf(" xxx DAY(%c-%3.3d)=%8.3f:  I2MAG = ", 
	     LCLIB_EVENT.STRING_FLAG[row], row, LCLIB_EVENT.DAY[row] );
      for(ifilt=0; ifilt < NFILT; ifilt++ ) 
	{ printf("%d ", LCLIB_EVENT.I2MAG[ifilt][row] ); }
      printf("\n"); fflush(stdout);
    } 
    printf(" xxx ----------- END DUMP in %s ---------------- \n\n", fnam );
    fflush(stdout);
  }

  return ;
} // end dumpEvent_LCLIB

// =========================================
void  malloc_LCLIB_EVENT(int OPT) {

  // OPT>0 --> allocate memory
  // OPT<0 --> free mem

  int  NROW  = LCLIB_EVENT.NROW +1000 ; 
  int  MEMD  = NROW * sizeof(double);
  int  MEMF  = NROW * sizeof(float);
  int  MEMI2 = NROW * sizeof(short int) ;
  int  MEMC  = NROW * sizeof(char);
  int  MEMTOT = 0 ;

  int  NFILT  = LCLIB_INFO.NFILTERS;
  int  NEP_T  = 0 ;
  int  ifilt; 
  char fnam[] = "malloc_LCLIB_EVENT" ;

  int IFLAG_PERIODIC = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_PERIODIC) ;
  int IFLAG_NONRECUR = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR );

  // ------------- BEGIN ---------------

  /*
  printf(" xxx %s malloc with OPT=%3d  (EVENT_ID=%d)\n",
	 fnam, OPT, LCLIB_EVENT.ID ); fflush(stdout);
  */
  
  if ( OPT > 0 ) {

    if ( NROW == 0 || NROW >= MXROW_LCLIB ) {
      sprintf(c1err,"Invalid NROW=%d", NROW);
      sprintf(c2err,"Check LCLIB and MXROW_LCLIB=%d", MXROW_LCLIB);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }


    if ( IFLAG_PERIODIC ) {  NEP_T  = LCLIB_INFO.NEP_TEMPLATE ; }
    if ( IFLAG_NONRECUR ) {  NEP_T  = 500 ; } // BEWARE fragile malloc

    if ( NEP_T > 0 ) {
      // increase allocation to include template rows      
      NROW   += NEP_T ;
      MEMD   += NEP_T * sizeof(double) ;
      MEMF   += NEP_T * sizeof(float) ;
      MEMI2  += NEP_T * sizeof(short int);
      MEMC   += NEP_T * sizeof(char) ;
    }

    MEMTOT += MEMD; LCLIB_EVENT.DAY          = (double*) malloc ( MEMD );
    MEMTOT += MEMC; LCLIB_EVENT.STRING_FLAG  = (char  *) malloc ( MEMC );

    LCLIB_EVENT.I2MAG  = (short int**) malloc ( NFILT*sizeof(short int*) );
    for(ifilt=0; ifilt < NFILT; ifilt++ ) {
      MEMTOT += MEMI2; 
      LCLIB_EVENT.I2MAG[ifilt]  = (short int*) malloc(MEMI2) ;
    }

    if ( MEMTOT > LCLIB_EVENT.MEMLC_MAX ) 
      { LCLIB_EVENT.MEMLC_MAX=MEMTOT; LCLIB_EVENT.NROW_MAX=NROW; }

  }
  else {
    free(LCLIB_EVENT.DAY) ;
    free(LCLIB_EVENT.STRING_FLAG) ;
    for(ifilt=0; ifilt < NFILT; ifilt++ ) 
      { free(LCLIB_EVENT.I2MAG[ifilt]);  }
    free(LCLIB_EVENT.I2MAG);
  }
  
  return ;

} // end malloc_LCLIB_EVENT

// ======================================
void  addTemplateRows_PERIODIC(void) {

  // For PERIODIC LCLIB, add template rows if none are given
  // in the LCLIB. Beware that template rows are added after
  // the search rows.
  //

  int    NFILT     = LCLIB_INFO.NFILTERS;
  int    NEP_T     = LCLIB_INFO.NEP_TEMPLATE ;
  double EPMIN_T   = LCLIB_INFO.EPLIST_TEMPLATE[0] ;
  double EPMAX_T   = LCLIB_INFO.EPLIST_TEMPLATE[NEP_T-1] ;
  double EPRANGE   = EPMAX_T - EPMIN_T ;
  double mag;
  double DAY_S, DAYFIRST_S, DAYLAST_S, NDAY_S ; 
  double DAY_T, DAYFIRST_T, DAYLAST_T, NDAY_T ; 
  
  int  FIRSTROW_S = LCLIB_EVENT.FIRSTROW_S ;
  int  LASTROW_S  = LCLIB_EVENT.LASTROW_S ;
  int  FIRSTROW_T ;
  int  NROW_S     = LCLIB_EVENT.NROW_S ;
  int  NROW_T     = LCLIB_EVENT.NROW_T ;
  int  irow, IROW, NCYCLE_S, ifilt ;
  int  LDMP = 0 ;
  char fnam[]     = "addTemplateRows_PERIODIC" ;

  // -------------- BEGIN -------------

  if ( NROW_T > 0 ) { return ; }

  DAYFIRST_S = LCLIB_EVENT.DAY[FIRSTROW_S] ;
  DAYLAST_S  = LCLIB_EVENT.DAY[LASTROW_S] ;
  NDAY_S     = DAYLAST_S - DAYFIRST_S ;

  DAYFIRST_T = DAYFIRST_S - EPRANGE - 10.0 ;
  DAYLAST_T  = DAYFIRST_T + EPRANGE ;
  NDAY_T     = DAYLAST_T - DAYFIRST_T ;
  NROW_T     = NEP_T ;

  if ( LDMP ) {
    printf("\n xxx -------- %s DUMP EVENT %lld -------------- \n", 
	   fnam, LCLIB_EVENT.ID  );
    printf(" xxx Construct %d PERIODIC template rows from %.2f to %.2f days\n",
	   NROW_T, DAYFIRST_T, DAYLAST_T );
  }

  FIRSTROW_T             = LASTROW_S + 1;
  
  for(irow=0; irow < NROW_T; irow++ ) {
    IROW  = FIRSTROW_T + irow ;
    DAY_T = DAYFIRST_T + (LCLIB_INFO.EPLIST_TEMPLATE[irow] - EPMIN_T);

    // find correspond search-day, which is really the phase
    NCYCLE_S = (int)((DAY_T-DAYFIRST_T)/NDAY_S) ;
    DAY_S = DAYFIRST_S + (DAY_T-DAYFIRST_T) - ((float)NCYCLE_S)*NDAY_S ;

    if ( DAY_S < DAYFIRST_S || DAY_S > DAYLAST_S ) {
      print_preAbort_banner(fnam);
      printf("\t DAYRANGE(S) = %8.2f to %8.2f \n", DAYFIRST_S, DAYLAST_S);
      printf("\t DAYRANGE(T) = %8.2f to %8.2f \n", DAYFIRST_T, DAYLAST_T);
      printf("\t irow=%d  IROW=%d \n", irow, IROW);
      sprintf(c1err,"Invalid Search Phase = %.2f for EVENT ID=%lld", 
	      DAY_S, LCLIB_EVENT.ID );
      sprintf(c2err,"DAY_T=%.2f   NCYCLE_S=%d", 
	      DAY_T, NCYCLE_S );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    LCLIB_EVENT.DAY[IROW]         =  DAY_T ;
    LCLIB_EVENT.STRING_FLAG[IROW] =  'T' ;

    // get template mag from search phase (DAY_S)    
    for(ifilt=0; ifilt < NFILT; ifilt++ ) {
      mag  = magInterp_LCLIB(DAY_S, NROW_S, 
			     &LCLIB_EVENT.DAY[FIRSTROW_S], 
			     &LCLIB_EVENT.I2MAG[ifilt][FIRSTROW_S] );
      LCLIB_EVENT.I2MAG[ifilt][IROW] = (int)( mag*I2FLOAT_LCLIB + 0.5 );
    }

    if ( LDMP ) {
      printf("\t xxx Add Template DAY[%2d] = %6.2f --> SearchPhase=%.2f\n", 
	     IROW, DAY_T, DAY_S );
    }

  } // end row loop


  LCLIB_EVENT.FIRSTROW_T = FIRSTROW_T ;;
  LCLIB_EVENT.LASTROW_T  = FIRSTROW_T + NROW_T -1 ;
  LCLIB_EVENT.NROW_T  = NROW_T ;
  LCLIB_EVENT.NROW   += NROW_T ;

  return ;

} // end addTemplateRows_PERIODIC


// ========================================
void addTemplateRows_NONRECUR(void) {

  // Createe Dec 2017
  // If there is only 1 template row for NONRECUR event,
  // add extra S rows to cover the survey-time range.
  // Note that the added rows are S , not T, but these S rows
  // have the template (T) mag values. This allows
  // forceTemplateRows to work the same wasy as for
  // PERIODIC events.  And the existing T row -> S row.
  //
  // July 24 2018: check that first and last mags are the same

  int  NROW_ORIG    = LCLIB_EVENT.NROW ;
  int  NFILT        = LCLIB_INFO.NFILTERS ;
  int  NROW_T       = LCLIB_EVENT.NROW_T ;

  int    ifilt, irow, irow2, NROW_S, NFERR;
  int    LDMP, NROWADD, I2MAG_STORE[MXFILTINDX] ;
  double TobsRange, MAGDIF, m0, m1 ;
  double DAYMIN_TEMPLATE, DAYADD_TEMP, DAYBIN, DAYRANGE[2] ;
  char fnam[] = "addTemplateRows_NONRECUR" ;

  // ------------- BEGIN ----------

  // if more than 1 T: row in library, bail.
  if ( NROW_T != 1 ) { return ; }
  

  TobsRange = LCLIB_INFO.TOBS_RANGE_MAX + 3.0*DAYBACK_TEMPLATE_LCLIB;

  LCLIB_EVENT.NROWADD_NONRECUR = 
    (int)(TobsRange/DAYBACK_TEMPLATE_LCLIB) + 2 ;  
  NROWADD = LCLIB_EVENT.NROWADD_NONRECUR ;

  DAYMIN_TEMPLATE = 
    LCLIB_EVENT.DAYRANGE_T[0] - TobsRange ;

  DAYBIN = TobsRange/(double)(NROWADD) ;

  // store template mag in each filter, and check that
  // first and last mags are the same.

  NFERR=0;
  for(ifilt=0; ifilt < NFILT; ifilt++   ) {
    I2MAG_STORE[ifilt] = LCLIB_EVENT.I2MAG[ifilt][0] ; 

    m0 = LCLIB_EVENT.FIRSTMAG[ifilt] ;
    m1 = LCLIB_EVENT.LASTMAG[ifilt] ;
    MAGDIF = fabs(m0-m1) ;
    if ( MAGDIF > 0.02 ) {
      printf(" WARNING: EVENT_ID=%lld first/last mag(%c) = %.3f/%.3f "
	     "(DIF=%.3f)\n",
	     LCLIB_EVENT.ID, LCLIB_INFO.FILTERS[ifilt], m0, m1, MAGDIF );
      NFERR++ ;
    } 
  }

  if ( NFERR > 0 ) {
    sprintf(c1err,"Invalid NONRECUR EVENT (ID=%lld)", LCLIB_EVENT.ID );
    sprintf(c2err,"First and Last mags do not match");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // move original rows up to span NROWADD -  NROWADD+NROW_ORIG
  // Run loop backward to avoid clobbering
  for(irow=NROW_ORIG-1; irow >= 0 ; irow-- ) {
    irow2 = irow + NROWADD ;
    LCLIB_EVENT.DAY[irow2] = LCLIB_EVENT.DAY[irow]  ;
    for(ifilt=0; ifilt < NFILT; ifilt++   ) { 
      LCLIB_EVENT.I2MAG[ifilt][irow2] = LCLIB_EVENT.I2MAG[ifilt][irow] ;
    }
    LCLIB_EVENT.STRING_FLAG[irow2] = 'S' ;  // T -> S
  }

  // now add the extra rows to cover survey window
  for(irow=0; irow < NROWADD; irow++ ) {
    DAYADD_TEMP = DAYMIN_TEMPLATE + DAYBIN*(double)irow;
    LCLIB_EVENT.DAY[irow] = DAYADD_TEMP ;
    LCLIB_EVENT.STRING_FLAG[irow] = 'S' ; 
    for(ifilt=0; ifilt < NFILT; ifilt++   ) { 
      LCLIB_EVENT.I2MAG[ifilt][irow] = I2MAG_STORE[ifilt] ;
    }
  }

  // update globals
  NROW_S                    =  NROW_ORIG + NROWADD ;
  DAYRANGE[0]               =  LCLIB_EVENT.DAY[0] ;
  DAYRANGE[1]               =  LCLIB_EVENT.DAY[NROW_S-1] ;
  LCLIB_EVENT.NROW_S        =  NROW_S ;
  LCLIB_EVENT.NROW          =  NROW_S ;
  LCLIB_EVENT.FIRSTROW_S    =  0 ;
  LCLIB_EVENT.LASTROW_S     =  NROW_S - 1 ;
  LCLIB_EVENT.DAYRANGE_S[0] =  DAYRANGE[0];
  LCLIB_EVENT.DAYRANGE_S[1] =  DAYRANGE[1];
  LCLIB_EVENT.DAYCOVER_S    =  DAYRANGE[1] - DAYRANGE[0];
  LCLIB_EVENT.DAYCOVER_ALL  =  DAYRANGE[1] - DAYRANGE[0];

  LCLIB_EVENT.NROW_T        =  0 ;
  LCLIB_EVENT.FIRSTROW_T    = -9 ;
  LCLIB_EVENT.LASTROW_T     = -9 ;
  LCLIB_EVENT.DAYCOVER_T    =  0.0 ;

  
  /*
  printf(" xxx NROW_ORIG=%d  NROWADD=%d  DAYBIN=%.2f   DAYMIN=%.2f \n", 
	 NROW_ORIG, NROWADD, DAYBIN, DAYMIN_TEMPLATE );
  */

  // check for monotonic day
  for (irow=1; irow < NROW_S; irow++ ) {
    if ( LCLIB_EVENT.DAY[irow] <= LCLIB_EVENT.DAY[irow-1] ) {

      print_preAbort_banner(fnam);
      for (irow2=0; irow2 < NROW_S; irow2++ ) 
	{ printf("\t DAY[%2d] = %.2f \n", irow2, LCLIB_EVENT.DAY[irow2] ); }

      sprintf(c1err,"DAY not monotonic at row=%d  ID=%lld", 
	      irow, LCLIB_EVENT.ID );
      sprintf(c2err,"DAY[%d]=%.2f  DAY[%d]=%.2f  "
	      "(NROW_ORIG=%d, NROWADD=%d)\n",
	      irow-1, LCLIB_EVENT.DAY[irow-1] ,
	      irow, LCLIB_EVENT.DAY[irow],
	      NROW_ORIG, NROWADD);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
  }

  LDMP = 0 ;
  if ( LDMP ) {
    printf(" xxx ------------------------------ \n");
    printf(" xxx %s DUMP for EVENT ID=%lld \n", fnam, LCLIB_EVENT.ID);
    printf("\t xxx NROW = %d -> %d  (NROWADD=%d) \n", 
	   NROW_ORIG, NROW_S, NROWADD );
    for(irow=0; irow < NROW_S; irow++ ) {
      printf("\t xxx New S-DAY[%2d] = %6.2f \n", irow, LCLIB_EVENT.DAY[irow]);
    }
    //    debugexit(fnam);
  }

  return  ;

} // end addTemplateRows_NONRECUR


// ===========================================
void addTemplateReset_NONRECUR(void) {

  // Restore LCLIB event back to how it was read from file.
  // Make sure to use NROW_S = LASTROW_S+1 instead of .NROW_S

  int NFILT      = LCLIB_INFO.NFILTERS ;
  int NROWADD    = LCLIB_EVENT.NROWADD_NONRECUR  ;
  int NROW_S     = LCLIB_EVENT.LASTROW_S + 1 ; // current NROW_S = NROW
  int NROW_ORIG  = NROW_S - NROWADD ;  // includes 1 template row
  int ifilt, irow, irow2, LDMP ;
  double DAYRANGE[2];
  char fnam[] = "addTemplateReset_NONRECUR" ;

  // ------------- BEGIN -----------

  for(irow=0; irow < NROW_ORIG; irow++ ) {   
    irow2 = irow + NROWADD ;
    LCLIB_EVENT.DAY[irow] = LCLIB_EVENT.DAY[irow2];
    for(ifilt=0; ifilt < NFILT; ifilt++ ) { 
      LCLIB_EVENT.I2MAG[ifilt][irow] = LCLIB_EVENT.I2MAG[ifilt][irow2] ;
    }

    if ( irow == 0 ) 
      { LCLIB_EVENT.STRING_FLAG[irow] = 'S' ;  }
    else
      { LCLIB_EVENT.STRING_FLAG[irow] = 'T' ;  }
  }

  DAYRANGE[0] = LCLIB_EVENT.DAY[1];
  DAYRANGE[1] = LCLIB_EVENT.DAY[NROW_ORIG-1];

  LDMP = 0;
  if ( LDMP ) {
    printf(" xxx %s: NROW_ORIG=%d  NROWADD=%d  NROW=%d  NROW_S=%d \n",
	   fnam, NROW_ORIG, NROWADD, 
	   LCLIB_EVENT.NROW, LCLIB_EVENT.NROW_S);
  }

  if ( NROW_ORIG < 0 ) {
    print_preAbort_banner(fnam);
    printf("\t NROW_ORIG=%d  NROWADD=%d  NROW=%d  NROW_S=%d \n",
	   NROW_ORIG, NROWADD, LCLIB_EVENT.NROW, LCLIB_EVENT.NROW_S);
    printf("\t ROW_S[FIRST,LAST] = %d, %d \n", 
	   LCLIB_EVENT.FIRSTROW_S, LCLIB_EVENT.LASTROW_S );
    printf("\t DAYRANGE = %.1f, %.1f \n", DAYRANGE[0], DAYRANGE[1]);
    sprintf(c1err,"NROW_ORIG < 0 --> Something is really messed up.");
    c2err[0] = 0 ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  
  LCLIB_EVENT.NROW_S     = NROW_ORIG-1; // subtract template row
  LCLIB_EVENT.FIRSTROW_S = 1;
  LCLIB_EVENT.LASTROW_S  = NROW_ORIG-1;
  LCLIB_EVENT.DAYRANGE_S[0] = DAYRANGE[0] ;
  LCLIB_EVENT.DAYRANGE_S[1] = DAYRANGE[1] ;
  LCLIB_EVENT.DAYCOVER_S    = DAYRANGE[1] - DAYRANGE[0] ; 
  LCLIB_EVENT.DAYCOVER_ALL  = DAYRANGE[1] - DAYRANGE[0] ;

  LCLIB_EVENT.NROW_T        = 1 ;
  LCLIB_EVENT.FIRSTROW_T    = 0 ;
  LCLIB_EVENT.LASTROW_T     = 0 ;
  LCLIB_EVENT.DAYRANGE_T[0] = LCLIB_EVENT.DAY[0];
  LCLIB_EVENT.DAYRANGE_T[1] = LCLIB_EVENT.DAY[1];
  LCLIB_EVENT.DAYCOVER_T    = 0.0 ;

  LCLIB_EVENT.NROW = LCLIB_EVENT.NROW_S + LCLIB_EVENT.NROW_T  ;

  return ;

} // end addTemplateReset_NONRECUR

// ===========================================
void set_TOBS_OFFSET_LCLIB(void) {
  
  // Created Oct 18 2017
  // Set TOBS_OFFSET offset to shift LCLIB map onto
  // simulated Tobs range.
  //
  // Note that Tobs refers to simulated range, while DAY
  // refers to LCLIB range.
  //

  int    NROW_T    = LCLIB_EVENT.NROW_T ;
  int    PERIODIC  = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_PERIODIC) ;
  int    NONRECUR  = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR) ;
  double TGENRANGE_NONRECUR = LCLIB_EVENT.TGENRANGE_NONRECUR ;

  double  Tobs_min  = LCLIB_EVENT.TOBS_MIN ;
  double  Tobs_max  = LCLIB_EVENT.TOBS_MAX ;
  double  TobsRange = LCLIB_EVENT.TOBS_RANGE ;

  double TOFF0 = LCLIB_DEBUG.TOBS_OFFSET_RANGE[0] ;
  double TOFF1 = LCLIB_DEBUG.TOBS_OFFSET_RANGE[1] ;

  int    FIRSTROW_S, row ;
  int    LDMP=0 ; 
  double TOBS_OFFSET = 0.0 ;
  double FlatRan ;
  double DAYRANGE_SPARE, DAY_RANDOM, DAY_FIRST, DAY_BUFFER;
  double DAY, DAY_FIRST_SHIFT ;
  char fnam[] = "set_TOBS_OFFSET_LCLIB";

  // -------------- BEGIN ------------

  FlatRan = getRan_Flat1(1);
  DAY_BUFFER = 0.0 ;

  // no offset for DEBUG
  if ( LCLIB_INFO.DEBUGFLAG_RANMAG ) { 
    // do nothing
  }
  else if ( NONRECUR ) {
    // non-recurring --> random time

    double Tstart_max = LCLIB_EVENT.DAYRANGE_S[1] ;    
    double Tstart_min = Tstart_max - TGENRANGE_NONRECUR;
    DAY_RANDOM = Tstart_min + (FlatRan*TGENRANGE_NONRECUR) ;

    /*
    printf(" xxx - - - - - - - - - - - \n");
    printf(" xxx %s: Tstart[min,max] = %.2f, %.2f  (DAY_RAN=%.2f) \n",
	   fnam, Tstart_min, Tstart_max, DAY_RANDOM );
    printf(" xxx %s: TGENRANGE_NONRECUR = %.2f \n", 
	   fnam, TGENRANGE_NONRECUR ); 
    */


    if ( TOFF0 != TOFF1 ) {
      DAY_RANDOM = TOFF0 + (TOFF1-TOFF0) * FlatRan ;
    }

    // set DAY_RANDOM for forceTemplates
    LCLIB_EVENT.DAY_RANDOM = DAY_RANDOM ;
    TOBS_OFFSET = DAY_RANDOM - Tobs_min ;
  }
  else {
    // must be recurring  here

    // make sure that LCLIB is big enough to cover Tobs
    if ( TobsRange > LCLIB_EVENT.DAYCOVER_S &&  PERIODIC == 0 ) {
      sprintf(c1err,"LCLIB DAYRANGE=%.1f days (ID=%lld) cannot cover",
	      LCLIB_EVENT.DAYCOVER_S, LCLIB_EVENT.ID ) ;
      sprintf(c2err,"TobsRange = %.1f days  (%.2f to %.2f)", 
	      TobsRange, Tobs_min, Tobs_max );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    // compute extra LCLIB days to cover TobsRange, and random offset
    
    if ( PERIODIC ) 
      { DAYRANGE_SPARE = LCLIB_EVENT.DAYCOVER_S ; }     
    else 
      { DAYRANGE_SPARE = LCLIB_EVENT.DAYCOVER_S  - TobsRange ; }

    // - - - - - - - 
    // if there are no template rows, leave enough BUFFER 
    // room to convert Search rows into template rows.

    int NBUF=0,  ibuf ; 
    double DAY_BUFFER_LIST[10];
    if ( NROW_T == 0 && PERIODIC==0 ) { 
      // user template  range
      DAY_BUFFER_LIST[NBUF] = LCLIB_INFO.EPRANGE_TEMPLATE; NBUF++ ;
      // time between template and Search
      DAY_BUFFER_LIST[NBUF] = DAYBACK_TEMPLATE_LCLIB;  NBUF++ ;
      // padding
      DAY_BUFFER_LIST[NBUF] = 10.0;  NBUF++ ;

      for(ibuf=0; ibuf < NBUF; ibuf++ ) 
	{ DAY_BUFFER += DAY_BUFFER_LIST[ibuf]; }
      DAYRANGE_SPARE -= DAY_BUFFER ;
    }


    if ( DAYRANGE_SPARE < 0.0 || DAYRANGE_SPARE > LCLIB_EVENT.DAYCOVER_S ) {
      sprintf(c1err,"Invalid DAYRANGE_SPARE=%.2f for EVENT_ID=%lld", 
	      DAYRANGE_SPARE, LCLIB_EVENT.ID );
      sprintf(c2err,"DAYCOVER_S=%.2f  DAY_BUFFER=%.2f", 
	      LCLIB_EVENT.DAYCOVER_S, DAY_BUFFER );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    // get start DAY in LCLIB
    FIRSTROW_S   = LCLIB_EVENT.FIRSTROW_S;
    DAY_FIRST    = LCLIB_EVENT.DAY[FIRSTROW_S];
    DAY_FIRST_SHIFT = DAY_FIRST + DAY_BUFFER ;

    // make sure that DAY_FIRST_SHIFT is past a real epoch which is
    // more than DAY_BUFFER past DAY_FIRST
    
    for(ibuf=0; ibuf < NBUF; ibuf++  ){
      row=0;   DAY=LCLIB_EVENT.DAY[row];
      while(DAY < DAY_FIRST_SHIFT) {
       	if ( row >= LCLIB_EVENT.NROW ) {
	  sprintf(c1err,"row=%d >= NROW=%d", row, LCLIB_EVENT.NROW);
	  sprintf(c2err,"Cannot compute DAY_FIRST_SHIFT for ID=%lld",
		  LCLIB_EVENT.ID );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}
	DAY = LCLIB_EVENT.DAY[row];  
	row++ ; 	
      }
      if ( DAY_FIRST_SHIFT<DAY ) 
	{ DAY_FIRST_SHIFT = DAY + DAY_BUFFER_LIST[ibuf] ; }
    }

    // - - - - - - - - - - - - - - - - 

    if ( DAYRANGE_SPARE > DAY_FIRST_SHIFT ) 
      { DAYRANGE_SPARE -= DAY_FIRST_SHIFT ; }
    else
      { DAYRANGE_SPARE = 0.0 ; }

    DAY_RANDOM  = DAY_FIRST_SHIFT + (DAYRANGE_SPARE * FlatRan) ;
    LCLIB_EVENT.DAY_RANDOM = DAY_RANDOM ;

    TOBS_OFFSET = DAY_RANDOM - Tobs_min ;
    
    if ( LDMP ) {
      printf(" xxx ----------------------------------------- \n");
      printf(" xxx TOBS_OFFSET DUMP for EVENT_ID=%lld  \n", 
	     LCLIB_EVENT.ID);
      printf(" xxx FlatRan = %f \n", FlatRan );
      printf(" xxx DAY_FIRST = %.2f->%.2f    DAY_RANDOM=%.2f of %.2f  "
	     "DAYCOVER=%.1f\n",
	     DAY_FIRST, DAY_FIRST_SHIFT, 
	     DAY_RANDOM, DAYRANGE_SPARE, 
	     LCLIB_EVENT.DAYCOVER_S);
      printf(" xxx DAY_BUFFER = %.3f \n", DAY_BUFFER);
      printf(" xxx Tobs(min,max) = %.2f,%.2f  TobsRange=%.2f \n",
	     Tobs_min, Tobs_max, TobsRange);
      printf(" xxx TOBS_OFFSET = %.2f \n", TOBS_OFFSET );
    }

  }

  LCLIB_EVENT.TOBS_OFFSET = TOBS_OFFSET ;

  return ;

} // end TOBS_OFFSET_LCLIB


// =================================================
void  get_TobsMinMax(int NOBS, double *TobsList, 
		     double *Tobs_min, double *Tobs_max) {
  // For input TobsList,
  // return Tobs_min and Tobs_max 

  double Tobs, Tmin= 99999999.0, Tmax = -99999.0 ;
  int obs;
  // ------------ BEGIN --------------

  *Tobs_min  = *Tobs_max = -99999 ;

  for(obs=0; obs < NOBS; obs++ ) {
    Tobs         = TobsList[obs] ;
    if ( Tobs < Tmin ) { Tmin = Tobs; }
    if ( Tobs > Tmax ) { Tmax = Tobs; }
  }

  *Tobs_min = Tmin ;
  *Tobs_max = Tmax ;

  return ;
} // end get_TobsMinMax

// ======================================
double magSearch_LCLIB(int ifilt, double Tobs) {

  // determine search-mag for input filter 'ifilt' and epoch Tobs.
  // Tobs has already been shifted to lie within LCLIB.DAYRANGE.
  // Be careful to treat RECURRING and NON-RECURRING.
  //
  // July 24 2018: fix bug setting NON-RECURR mag for DAY_LCLIB >= DAYMAX_S
  //

  int IFLAG_PERIODIC = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_PERIODIC) ;
  int IFLAG_NONRECUR = (LCLIB_INFO.IFLAG_RECUR_CLASS == IFLAG_RECUR_NONRECUR ) ;

  int FIRSTROW_S   = LCLIB_EVENT.FIRSTROW_S ;
  int NROW_S       = LCLIB_EVENT.NROW_S ;
  double mag_T     = LCLIB_EVENT.magTemplate[ifilt];
  double DAYMIN_S  = LCLIB_EVENT.DAYRANGE_S[0];
  double DAYMAX_S  = LCLIB_EVENT.DAYRANGE_S[1];
  double NDAY_S    = DAYMAX_S - DAYMIN_S ;
  double mag, DAY_LCLIB ;
  int    NCYCLE_S;
  char fnam[] = "magSearch_LCLIB" ;

  // ------------- BEGIN ---------------

  if ( LCLIB_INFO.DEBUGFLAG_RANMAG ) {
    double DIF0 = LCLIB_INFO.GENRANGE_DIFMAG[0];
    double DIF1 = LCLIB_INFO.GENRANGE_DIFMAG[1];
    mag    = mag_T + DIF0 + (DIF1-DIF0) * getRan_Flat1(1);
    return(mag) ;
  }

  DAY_LCLIB = Tobs ;
  
  if ( IFLAG_NONRECUR == 0  ) {
    // RECURRING

    if ( IFLAG_PERIODIC ) {
      NCYCLE_S   = (int)(( Tobs-DAYMIN_S)/NDAY_S ) ;
      DAY_LCLIB -= ( NDAY_S * (float)NCYCLE_S );
      
      /*
      printf(" xxx NCYCLE=%3d : DAY_LCLIB=%7.3f -> %7.3f \n", 
      NCYCLE_S, Tobs, DAY_LCLIB);   */
    }

    // make sure that DAY_LCLIB is covered by LCLIB; otherwise abort.
    if ( DAY_LCLIB < DAYMIN_S || DAY_LCLIB > DAYMAX_S ) {
      sprintf(c1err,"Cannot get RECURRING mag(%c) at Tobs=%.3f "
	      "for EVENT_ID=%lld",
	      LCLIB_INFO.FILTERS[ifilt], DAY_LCLIB, LCLIB_EVENT.ID);
      sprintf(c2err,"LCLIB.DAYRANGE = %.2f to %.2f is too narrow", 
	      DAYMIN_S, DAYMAX_S);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
  }
  else {
    // NON-RECURRING (SN-like):
    // Before LCLIB epoch range, set mag_S = mag_T
    //    (beware that mag_T is not necessarily quiescent mag)
    // After  LCLIB epoch range. set mag_S = last [quiescent] mag
    // (i.e., the quiescent value)

    if ( DAY_LCLIB <= DAYMIN_S ) { return(mag_T) ; }
    if ( DAY_LCLIB >= DAYMAX_S ) { DAY_LCLIB = DAYMAX_S - 1.0E-5 ; }
  }

  // ----------------------------------------------
  //linear interpolate MAG-vs-DAY grid.
  mag  = magInterp_LCLIB(DAY_LCLIB, NROW_S, 
			 &LCLIB_EVENT.DAY[FIRSTROW_S], 
			 &LCLIB_EVENT.I2MAG[ifilt][FIRSTROW_S] );
  
  return(mag) ;

} // end magSearch_LCLIB

// ======================================
double magTemplate_LCLIB(int EXTERNAL_ID, int ifilt) {

  // Return template mag for this sparse LCLIB-filter index
  // EXTERNAL_ID is for debug.

  int NEP_TEMPLATE   = LCLIB_INFO.NEP_TEMPLATE ;
  int IFLAG_PERIODIC = (LCLIB_INFO.IFLAG_RECUR_CLASS== IFLAG_RECUR_PERIODIC);

  int NROW_T        = LCLIB_EVENT.NROW_T ;
  int FIRSTROW_T    = LCLIB_EVENT.FIRSTROW_T ;
  int LASTROW_T     = LCLIB_EVENT.LASTROW_T ;
  int Nforce_T      = LCLIB_EVENT.NforceTemplateRows ;

  double Tep, fluxSum=0.0, flux, fluxAvg, arg, mag ;
  int    ep, I2MAG ;
  int    LDMP = (EXTERNAL_ID < -1 );
  char fnam[] =  "magTemplate_LCLIB";

  // ----------------- BEGIN ------------------

  if ( LCLIB_INFO.DEBUGFLAG_RANMAG ) {
    double MAG0 = LCLIB_INFO.GENRANGE_RANMAG[0];
    double MAG1 = LCLIB_INFO.GENRANGE_RANMAG[1];
    mag = MAG0 + ( MAG1 - MAG0 ) * getRan_Flat1(1) ;
    //    printf(" xxx mag_T = %f (%f - %f)\n", mag, MAG0, MAG1 );
    return(mag);
  }

  if ( LDMP ) { 
    printf(" xxx ------------------------------------- \n");
    printf(" xxx %s: EXTERN_ID=%d  DAY_RANDOM = %.2f \n", 
	   fnam, EXTERNAL_ID, LCLIB_EVENT.DAY_RANDOM); 
  }

  // check for non-recurring transient that has just
  // one template epoch
  if ( NROW_T == 1 ) {
    I2MAG = LCLIB_EVENT.I2MAG[ifilt][FIRSTROW_T] ;
    mag = ((double)I2MAG) / I2FLOAT_LCLIB ;
  }
  else {
    // get flux-average among template epochs, and then
    // convert flux-average back to mag.

    double TMIN = (double)LCLIB_EVENT.DAY[FIRSTROW_T] ;
    double TMAX = (double)LCLIB_EVENT.DAY[LASTROW_T] ;

    fluxSum = 0.0 ;
    for(ep=0; ep < NEP_TEMPLATE; ep++ ) {
      Tep = LCLIB_INFO.EPLIST_TEMPLATE[ep];

      if ( IFLAG_PERIODIC || Nforce_T > 0) 
	{ Tep -= (LCLIB_INFO.EPLIST_TEMPLATE[0] - TMIN); }

      if(Tep < TMIN-1.0E-9 || Tep > TMAX+1.0E-9 ) {
	sprintf(c1err,"Invalid user Tep(template)=%.2f for ID=%lld",
		Tep, LCLIB_EVENT.ID );
	sprintf(c2err,"LCLIB-Template epoch range: %.2f to %.2f",
		TMIN, TMAX ) ;
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      mag = magInterp_LCLIB(Tep, NROW_T, 
			    &LCLIB_EVENT.DAY[FIRSTROW_T], 
			    &LCLIB_EVENT.I2MAG[ifilt][FIRSTROW_T] );

      arg    = 0.4*(mag-ZEROPOINT_FLUXCAL_DEFAULT);
      flux   = pow(10.0,-arg);
      fluxSum += flux ;

      if ( LDMP ) {
	printf(" xxx magT: ID=%lld  ifilt=%d  mag_T=%.3f  "
	       "Tep=%.3f (DAY0=%.3f) \n", 
	       LCLIB_EVENT.ID, ifilt, mag, 
	       Tep, LCLIB_EVENT.DAY[FIRSTROW_T] );
      }

    }
  }

  fluxAvg = fluxSum / (double)NEP_TEMPLATE ;
  mag     = ZEROPOINT_FLUXCAL_DEFAULT - 2.5*log10(fluxAvg) ;

  //  printf(" xxx %s: ifilt=%d magTemplate = %f\n", fnam,ifilt, mag);

  return(mag) ;

} // end magTemplate_LCLIB


// ==================================
void store_magTemplate_LCLIB(int EXTERNAL_ID, int ifilt, double XT_MW) {

  // shell to call magTemplate_LCLIB and store it.
  // EXTERNAL_ID is for debug only.

  double mag_T ;

  if ( LCLIB_DEBUG.ZERO_TEMPLATE_FLUX ) {
    mag_T = MAG_ZEROFLUX ; // Dec 27 2018
  }
  else {
    mag_T  = magTemplate_LCLIB(EXTERNAL_ID,ifilt); 
    mag_T += XT_MW ;
  }

  LCLIB_EVENT.magTemplate[ifilt] = mag_T ; // global storage

} // end store_magTemplate_LCLIB 

// =========================================================
double magInterp_LCLIB(double T, int NROW, double *DAYLIST, 
		       short int *I2MAG) {

  // interpolate mag at epoch T.
  
  double mag, DAY, DAYFRAC, DAYSTEP, m0, m1 ;
  int    row, ROW ;
  char fnam[] = "magInterp_LCLIB" ;

  // --------------- BEGIN ---------------

  // find DAYLIST bin containing T
  row=0;    DAY = DAYLIST[row] ;
  while ( DAY < (T+1.0E-5) && row < NROW-1 ) 
    { row++ ;   DAY = DAYLIST[row];   }

  ROW = row-1;
  if ( ROW < 0 || ROW >= NROW ) {
    sprintf(c1err,"Interp problem: ROW=%d  not within 0-%d",
	    ROW, NROW-1 );
    sprintf(c2err,"T=%f  ID=%lld  NROW=%d", T, LCLIB_EVENT.ID, NROW );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  DAYSTEP  = DAYLIST[ROW+1] - DAYLIST[ROW] ;
  DAYFRAC  = (T-DAYLIST[ROW]) / DAYSTEP ;
  m0       = ((double)I2MAG[ROW+0])/ I2FLOAT_LCLIB ;
  m1       = ((double)I2MAG[ROW+1])/ I2FLOAT_LCLIB ;

  mag      = m0 + DAYFRAC*(m1-m0);

  /*
  printf(" xxx %s: T=%.2f  ROW=%d  DAYFRAC=%.2f  mag=%.3f \n",
	 fnam, T, ROW, DAYFRAC, mag );
  */
  
  return(mag) ;

} // end magInter_LCLIB
