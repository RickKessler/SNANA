/************************************
  Jan 17 2014 R.Kessler
 
  Pull out simulation trigger codes (PIPELINE and SPEC)
  into separate sntools_trigger.o to allow for other
  [non-sim] programs to use these trigger codes.
  
  Usage:
    Fill INPUTS_SEARCHEFF struct
    init_SEARCHEFF(survey);

    for each SN,
    * fill SEARCHEFF_DATA struct
    * gen_SEARCHEFF()


                HISTORY

 March 2018: 
   + lots of refactoring and new PHOTPROB map
   + read PHOTFLAG_TRIGGER from searcheff-pipeline file.

 Mar 28 2018:
    fix index bug for SEARCHEFF_DATA.SBMAG[ifiltobs];
    index is filter-index, not obs-index. Affect new PHOTPROB feature.
 
  Jun 19 2018: 
    + read & apply  optional "CUTWIN_SNRMAX: xx yy" in zHOST effic file.

  Mar 2019
   + inclue sntools_host.h to access HOSTLIB variables (for zHOST)
   + major refactor & update for SEARCHEFF_zHOST; multi-D map
     as an arbitrary function of HOSTLIB params.

  Feb 2 1 2020:
   add PATH_USER_INPUT to path list for snana_openTextFile,
   and if not file then call abort_openTextFile.

  Feb 25 2020: store README LINES from 0 to N-1 instead of 1 to N.

  Aug 26 2020: 
    + pass OPTMASK to snana_openTextFile to check DOCANA; If the new
      sim-input "REQUIRE_DOCANA: 1" is set, then code aborts if first
      key is NOT a DOCUMENTATION key.

  Feb 05 2021:
    + use new MATCH_SEARCHEFF_FIELD(field_map) function to handle
      overlaps.

************************************/

#include "sntools.h"
#include "sntools_trigger.h"
#include "sntools_host.h" 

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>


// ************************************
void init_SEARCHEFF(char *SURVEY, int APPLYMASK_SEARCHEFF ) {

  /***********************************
    Created Jan 7, 2008 by R.Kessler

    Initialization for search-efficiency.
    Read SEARCHEFF_SNR_[filter].DAT files, or
    Read SEARCHEFF_MAG_[filter].DAT files,
    and store in SEARCH_EFF structure.
    If search-eff file does not exist, then leave
    NBIN=0 as a flag to set EFF = 1 for any SNR value.

  Inputs:
    SURVEN:   
       name of survey
    APPLYMASK_SEARCHEFF: 
       argument of sim-input APPLY_SEARCHEFF_OPT (only for error checking)


        HISTORY
   
  Mar 06 2018: 
    + new input argument APPLYMASK_SEARCHEFF

  *************/

  //  char  fnam[] = "init_SEARCHEFF"  ;
  int  NMAP=0 ;
  
  // ------------- BEGIN ----------

  sprintf(BANNER, "Initialize SEARCH EFFICIENCY for '%s' \n", SURVEY );
  print_banner( BANNER );

  INPUTS_SEARCHEFF.NMAP_DETECT   = 0 ;
  INPUTS_SEARCHEFF.NMAP_PHOTPROB = 0 ;
  INPUTS_SEARCHEFF.NMAP_SPEC     = 0 ;
  INPUTS_SEARCHEFF.NMAP_zHOST    = 0 ;
  INPUTS_SEARCHEFF.NREDUCED_CORR_PHOTPROB = 0 ;
  INPUTS_SEARCHEFF.NPHOTPROB_DUMP = 0 ;

  INPUTS_SEARCHEFF.CUTWIN_SNRMAX_zHOST[0] = -99.9 ;
  INPUTS_SEARCHEFF.CUTWIN_SNRMAX_zHOST[1] = +1.0E8 ;

  SEARCHEFF_FLAG        = 0 ;
  SEARCHEFF_LOGIC.NMJD  = 0;

  NONZERO_SEARCHEFF_SPEC  = 0 ;
  NONZERO_SEARCHEFF_zHOST = 0 ;

  sprintf(PATH_SEARCHEFF, "%s %s/models/searcheff", 
	  PATH_USER_INPUT, PATH_SNDATA_ROOT );

  if ( INPUTS_SEARCHEFF.FUNEFF_DEBUG ) {
    printf("\t Use FUNEFF_DEBUG = %d \n", INPUTS_SEARCHEFF.FUNEFF_DEBUG );
    //    return ;
  }

  // read single file to get pipeline efficiency vs. SNR or vs. MAG.
  // Returns number of maps that have an efficiency curve defined.
  NMAP = init_SEARCHEFF_PIPELINE(SURVEY) ;

  if ( NMAP > 0 ) { init_SEARCHEFF_LOGIC(SURVEY); }  // read detection logic


  // init prob of get spec-confirmation
  init_SEARCHEFF_SPEC(SURVEY); 

  // init prob of getting zHOST for unconfirmed SN
  init_SEARCHEFF_zHOST(SURVEY); // May 2014 

  // Mar 2018: check that user SEARCHEFF mask is possible
  check_APPLYMASK_SEARCHEFF(SURVEY,APPLYMASK_SEARCHEFF);


}  // end of init_SEARCHEFF


// *************************************
void  check_APPLYMASK_SEARCHEFF(char *SURVEY, int APPLYMASK_SEARCHEFF_USER) {

  // Created March 2018
  // Compute SEARCHEFF-mask that is possible, 
  // and check that it is compatible with user-input
  // APPLYMASK_SEARCHEFF that is the argument of APPLY_SEARCHEFF_OPT.
  // Abort if APPLYMASK_SEARCHEFF_USER is impossible to obtain.
  //
  // July 6 2018: fix bug related to "ZERO" option for spec efficiency.
  // Sep 4 2019: write COMMENT_README_SEARCHEFF

  int  APPLYMASK_ALLOWED, OVP ;
  int  IFLAG_SPEC_EFFZERO = INPUTS_SEARCHEFF.IFLAG_SPEC_EFFZERO;
  int  NMAP_SPEC          = INPUTS_SEARCHEFF.NMAP_SPEC ;
  int  LSPEC              = IFLAG_SPEC_EFFZERO || NMAP_SPEC>0 ;

  char fnam[] = "check_APPLYMASK_SEARCHEFF" ;

  // ---------- BEGIN --------------

  // pipeline detection always allowed.
  APPLYMASK_ALLOWED = APPLYMASK_SEARCHEFF_PIPELINE ;

  if ( NONZERO_SEARCHEFF_SPEC ) 
    { APPLYMASK_ALLOWED += APPLYMASK_SEARCHEFF_SPEC ; }

  if ( NONZERO_SEARCHEFF_zHOST && LSPEC ) 
    { APPLYMASK_ALLOWED += APPLYMASK_SEARCHEFF_zHOST; }

  OVP = ( APPLYMASK_ALLOWED & APPLYMASK_SEARCHEFF_USER) ;


  if ( OVP != APPLYMASK_SEARCHEFF_USER) {
    print_preAbort_banner(fnam);
    printf("\t APPLY_SEARCHEFF_OPT += %d --> detection pipeline.\n",
	   APPLYMASK_SEARCHEFF_PIPELINE);
    printf("\t APPLY_SEARCHEFF_OPT += %d --> SPEC confirmed.\n",
	   APPLYMASK_SEARCHEFF_SPEC );
    printf("\t APPLY_SEARCHEFF_OPT += %d --> zHOST.\n",
	   APPLYMASK_SEARCHEFF_zHOST );
    printf("\t NONZERO_SEARCHEFF_SPEC=%d   LSPEC=%d \n", 
	   NONZERO_SEARCHEFF_SPEC, LSPEC);
    printf("\t NONZERO_SEARCHEFF_zHOST=%d \n", NONZERO_SEARCHEFF_zHOST);

    sprintf(c1err,"Invalid user-input 'APPLY_SEARCHEFF_OPT: %d' (OVP=%d) ",
	    APPLYMASK_SEARCHEFF_USER, OVP) ;
    sprintf(c2err,"Simulation can produce SEARCHEFF mask %d \n",
	    APPLYMASK_ALLOWED);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  // - - - - - - 
  // Sep 4 2019: prepare README comment for trigger.

  sprintf(COMMENT_README_SEARCHEFF, "APPLY_SEARCHEFF_OPT=%d --> ",
	  APPLYMASK_SEARCHEFF_USER);

  int REQ_PIPE, REQ_SPEC, REQ_zHOST, NUSE=0;
  char ctmp[100], clist[100], cplus[] = "+" ;
  REQ_PIPE  = (APPLYMASK_SEARCHEFF_USER & APPLYMASK_SEARCHEFF_PIPELINE );
  REQ_SPEC  = (APPLYMASK_SEARCHEFF_USER & APPLYMASK_SEARCHEFF_SPEC );
  REQ_zHOST = (APPLYMASK_SEARCHEFF_USER & APPLYMASK_SEARCHEFF_zHOST );

  ctmp[0] = clist[0] = 0; 
  if ( REQ_PIPE ) 
    { strcat(clist,"PIPELINE"); NUSE++; }

  if ( REQ_SPEC ) {
    if ( NUSE >0 ) { strcat(clist,cplus); }
    strcat(clist,"SPEC");  NUSE++ ;
  }
  if ( REQ_zHOST ) { 
    if ( NUSE > 0 ) { strcat(clist,cplus); }
    strcat(clist,"zHOST");   NUSE++ ;
  }

  if ( NUSE > 0 ) 
    { sprintf(ctmp,"Require EFF(%s)", clist); }
  else
    { sprintf(ctmp,"No trigger requirements"); }

  strcat(COMMENT_README_SEARCHEFF,ctmp);

  
  return ;

} // end  check_APPLYMASK_SEARCHEFF


// *******************************************
int init_SEARCHEFF_PIPELINE(char *survey) {
 
  // Read and init two kinds of maps:
  //   1. DETECTION efficiency map vs. SNR or MAG
  //   2. PHOTPROB map vs. SNR, host properties, etc ...
  //
  // Abort if user-requested file does not exist.
  // If default survey-dependent file does not exist
  // then just return 0.
  //
  // Nov 23, 2014: read optional PHOTGLAG_DETECT key and set
  //                SEARCHEFF_PHOTFLAG_DETECT
  //
  // Mar 7 2018: refactor to include MAPNAME & FIELD keys,
  //             and new SEARCHEFF_DETECT structure.
  //
  // Aug 26 2020: pass OPTMASK to snana_open to check for DOCANA
  //

  int   OPTMASK = INPUTS_SEARCHEFF.OPTMASK_OPENFILE ;
  FILE *fp;
  int   IREQUIRE, gzipFlag, imap, NMAP=0 ;
  int   FOUNDMAP_DETECT=0, FOUNDMAP_PHOTPROB=0 ;
  char  file_local[MXPATHLEN], c_get[60], *ptrFile_user, *ptrFile_final ;   
  char  fnam[] = "init_SEARCHEFF_PIPELINE"  ;
    

  // ---------------- BEGIN ----------------

  SEARCHEFF_FLAG  = 0 ;
  MAPVERSION_SEARCHEFF_DETECT   = 0 ;
  MAPVERSION_SEARCHEFF_PHOTPROB = 0 ;

  for ( imap=0; imap < MXMAP_SEARCHEFF_DETECT; imap++ ) { 
    SEARCHEFF_DETECT[imap].NBIN  = 0 ; 
    SEARCHEFF_DETECT[imap].MAPNAME[0] = 0 ; 
  }

  ptrFile_user  = INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE ;
  ptrFile_final = INPUTS_SEARCHEFF.PIPELINE_EFF_FILE ;

  if ( strcmp(ptrFile_user,"DEFAULT") == 0 )   { 
    sprintf(file_local, "SEARCHEFF_PIPELINE_%s.DAT",  survey); 
    IREQUIRE = 0 ;
  }
  else if ( IGNOREFILE(ptrFile_user) ) {
    sprintf(file_local, "%s",  ptrFile_user); 
    IREQUIRE = 0 ;
  }
  else {
    sprintf(file_local, "%s",  ptrFile_user); 
    IREQUIRE = 1 ;
  }

  // use utility to check local dir and path.
  fp = snana_openTextFile(OPTMASK, PATH_SEARCHEFF, file_local, 
			  ptrFile_final, &gzipFlag ); // returned

  if ( fp == NULL ) { 

    if ( IREQUIRE ) {
      abort_openTextFile("SEARCHEFF_PIPELINE_EFF_FILE", 
			 PATH_SEARCHEFF, file_local, fnam);
    }
    else { 
      printf("\n  Optional SEARCHEFF_PIPELINE_FILE not found -> skip. \n");
      fflush(stdout);
      return 0 ; 
    }
  }


  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[0], "\t %s\n", 
	  ptrFile_final ); 

  SEARCHEFF_FLAG = 0 ;

  while( (fscanf(fp, "%s", c_get )) != EOF) {

    if ( strcmp(c_get,"PHOTFLAG_DETECT:") == 0 && 
	 INPUTS_SEARCHEFF.PHOTFLAG_DETECT == 0 ) 
      { readint(fp, 1, &INPUTS_SEARCHEFF.PHOTFLAG_DETECT );  }

    if ( strcmp(c_get,"PHOTFLAG_TRIGGER:") == 0 && 
	 INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER == 0 ) 
      { readint(fp, 1, &INPUTS_SEARCHEFF.PHOTFLAG_TRIGGER );  }

    // check reading detection efficiency map
    if ( !FOUNDMAP_PHOTPROB ) 
      { FOUNDMAP_DETECT = readMap_SEARCHEFF_DETECT(fp,c_get); }    

    
    // check PHOTPROB map
    if ( !FOUNDMAP_DETECT )
      { FOUNDMAP_PHOTPROB = readMap_SEARCHEFF_PHOTPROB(fp, c_get); }    
    

  } // end of fscan

  fclose(fp); // done reading

  // check info for each map and set README comments
  NMAP = INPUTS_SEARCHEFF.NMAP_DETECT;
  for ( imap=0; imap < NMAP; imap++ )  { check_SEARCHEFF_DETECT(imap); }

  NMAP = INPUTS_SEARCHEFF.NMAP_PHOTPROB;
  for ( imap=0; imap < NMAP; imap++ )  { check_SEARCHEFF_PHOTPROB(imap); }
  
  //  debugexit(fnam); // xxx REMOVE

  return(INPUTS_SEARCHEFF.NMAP_DETECT);

} // end of  init_SEARCHEFF_PIPELINE

// ************************************************
int readMap_SEARCHEFF_DETECT  (FILE *fp,  char *key) {

  // Created Mar 7 2018
  // Code moved from init_SEARCHEFF_PIPELINE as part of refactor,
  // to allow reading PHOTPROB map in future upgrades.
  //
  //

  double VAL, EFF;
  int  IREAD, imap, NBIN ;
  char cfilt[MXFILTINDX] ;
  //  char fnam[] = "readMap_SEARCHEFF_DETECT";

  // ------------ BEGIN -------------

  if ( strcmp(key,"ENDMAP:")==0 ) { 
    MAPVERSION_SEARCHEFF_DETECT = 0 ;
    imap = INPUTS_SEARCHEFF.NMAP_DETECT-1 ; 
    return(0) ;  
  }

  if ( strcmp(key,"MAPNAME_DETECT:")==0 ) {
    imap = malloc_NEXTMAP_SEARCHEFF_DETECT();
    readchar(fp, SEARCHEFF_DETECT[imap].MAPNAME );
    MAPVERSION_SEARCHEFF_DETECT = 2 ;
    return(1) ;
  }

  if ( strcmp(key,"FILTER:")==0 || strcmp(key,"BAND:")==0 ) {
    if ( MAPVERSION_SEARCHEFF_DETECT == 2 ) { 
      // new map style has MAPNAME_DETECT and ENDMAP keys
      imap = INPUTS_SEARCHEFF.NMAP_DETECT-1 ; 
    } 
    else {
      // old-style map starts with FILTER key
      imap = malloc_NEXTMAP_SEARCHEFF_DETECT(); 
      MAPVERSION_SEARCHEFF_DETECT = 1 ;
    }
    readchar(fp,cfilt);
    sprintf(SEARCHEFF_DETECT[imap].FILTERLIST, "%s", cfilt);
    return(1);
  }
  

  IREAD = 0;
  if ( strcmp(key,"SNR:") == 0 ) 
    { SEARCHEFF_FLAG = FLAG_EFFSNR_DETECT ; IREAD = 1; }    
  else if ( strcmp(key,"ABS(SNR):") == 0 ) 
    { SEARCHEFF_FLAG = FLAG_EFFABSSNR_DETECT ; IREAD = 1; }    
  else if ( strcmp(key,"MAG:") == 0 ) 
    { SEARCHEFF_FLAG = FLAG_EFFMAG_DETECT ; IREAD = 1 ; }
  

  if ( IREAD ) {
    imap  = INPUTS_SEARCHEFF.NMAP_DETECT-1 ;
    readdouble(fp, 1, &VAL );  // SNR or MAG
    readdouble(fp, 1, &EFF );  // efficiency
    NBIN = SEARCHEFF_DETECT[imap].NBIN ;
    SEARCHEFF_DETECT[imap].VAL[NBIN]  = VAL;
    SEARCHEFF_DETECT[imap].EFF[NBIN]  = EFF;
    SEARCHEFF_DETECT[imap].NBIN++ ; 
    return(1);
  }
  
  return(0);

} // end  readMap_SEARCHEFF_DETECT 


// ************************************************
int malloc_NEXTMAP_SEARCHEFF_DETECT(void) {
  // Created March 2018
  // malloc next map, and return imap index.
  int imap, MEMD ;
  //  char fnam[] = "malloc_NEXTMAP_SEARCHEFF_DETECT" ;
  // ------------- BEGIN -------------
  INPUTS_SEARCHEFF.NMAP_DETECT++ ; 
  imap  = INPUTS_SEARCHEFF.NMAP_DETECT-1;
  MEMD  = MXROW_SEARCHEFF_DETECT * sizeof(double);
  SEARCHEFF_DETECT[imap].VAL = (double*)malloc(MEMD);
  SEARCHEFF_DETECT[imap].EFF = (double*)malloc(MEMD);
  return(imap);
} // end malloc_NEXTMAP_SEARCHEFF_DETECT


// ************************************************
int  readMap_SEARCHEFF_PHOTPROB(FILE *fp,  char *key) {

  // Created March 2018
  // Read optional map used to determine PHOTPROB.
  // This map is expected to be based on fakes 
  // processed with a real imaging pipeline.
  //
  // Mar 14 2019: 
  //   + warning, do NOT use the new read_GRIDMAP() utility.
  //   + remove malloc_SEARCHEFF_TMPMAP2D and malloc TMPMAP2D locally.

  int LEGACY = ( MAPVERSION_SEARCHEFF_DETECT == 1 );
  int MXROW  = MXROW_SEARCHEFF_PHOTPROB ;
  int LDMP = 0 ;
  int imap, ivar, NVAR_TOT, NVAR_READ, MEMVAR, MEMROW;
  int NVAR_MAP, NFUN, NROW, IDMAP,  IVARABS, ibin ;
  double tmpVal, TMPVAL[MXVAR_SEARCHEFF_PHOTPROB];
  char VARNAME[20], *MAPNAME;
  char fnam[] = "readMap_SEARCHEFF_PHOTPROB" ;

  // --------------- BEGIN --------------

  //  printf(" xxx %s check key = %s \n", fnam, key);

  if ( strcmp(key,"MAPNAME_PHOTPROB:")==0 ) {

    if ( LEGACY ) {
      sprintf(c1err,"Cannot mix legacy DETECT map with PHOTPROB map.");
      sprintf(c2err,"DETECT map must have MAPNAME_DETECT and ENDMAP keys.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
    INPUTS_SEARCHEFF.NMAP_PHOTPROB++ ;
    MAPVERSION_SEARCHEFF_PHOTPROB = 1 ;
    imap  = INPUTS_SEARCHEFF.NMAP_PHOTPROB-1 ;
    readchar(fp, SEARCHEFF_PHOTPROB[imap].NAME);

    if ( LDMP ) {
      printf(" xxx --------------------------------------- \n");
      printf(" xxx %s: key='%s' NMAP=%d \n",  fnam, key, imap+1);
    }

    SEARCHEFF_PHOTPROB[imap].NVAR_TOT       = 0 ;
    SEARCHEFF_PHOTPROB[imap].NVAR_MAP       = 0 ;
    SEARCHEFF_PHOTPROB[imap].NFUN_CDF       = 0 ;
    SEARCHEFF_PHOTPROB[imap].NROW           = 0 ;
    SEARCHEFF_PHOTPROB[imap].REQUIRE_DETECTION = 0 ;
    SEARCHEFF_PHOTPROB[imap].CUTVAL         = -999.0 ;
    SEARCHEFF_PHOTPROB[imap].REDUCED_CORR      = 0 ;
    sprintf( SEARCHEFF_PHOTPROB[imap].FIELDLIST, "ALL") ;
    sprintf( SEARCHEFF_PHOTPROB[imap].FILTERLIST,"ALL") ;

    if ( imap == 0 ) {
      sprintf(VARDEF_SEARCHEFF_PHOTPROB[IVARABS_PHOTPROB_SNR],     "SNR"   );
      sprintf(VARDEF_SEARCHEFF_PHOTPROB[IVARABS_PHOTPROB_LOGSNR],  "LOGSNR");
      sprintf(VARDEF_SEARCHEFF_PHOTPROB[IVARABS_PHOTPROB_SBMAG],   "SBMAG" );
      sprintf(VARDEF_SEARCHEFF_PHOTPROB[IVARABS_PHOTPROB_GALMAG],  "GALMAG" );
      sprintf(VARDEF_SEARCHEFF_PHOTPROB[IVARABS_PHOTPROB_REDSHIFT],"REDSHIFT");
    }
    return(1) ;
  }

  // - - - - - - - - -

  if ( MAPVERSION_SEARCHEFF_PHOTPROB == 0 ) { return(0); }

  // - - - - - - - -  -

  imap  = INPUTS_SEARCHEFF.NMAP_PHOTPROB-1 ;

  if ( strcmp(key,"FIELD:")==0 || strcmp(key,"FIELDLIST:")==0 ) 
    { readchar(fp, SEARCHEFF_PHOTPROB[imap].FIELDLIST );  }

  else if ( strcmp(key,"FILTER:")==0 || strcmp(key,"BAND:")==0 ) 
    { readchar(fp, SEARCHEFF_PHOTPROB[imap].FILTERLIST );  }

  else if ( strcmp(key,"NVAR:")==0 )   { 
    readint(fp, 1, &NVAR_TOT );
    SEARCHEFF_PHOTPROB[imap].NVAR_TOT  = NVAR_TOT ;  

    MEMVAR = NVAR_TOT  * sizeof(double*);
    MEMROW = MXROW     * sizeof(double );
    SEARCHEFF_TMPMAP2D = (double**) malloc(MEMVAR); // global array
    for(ivar=0; ivar < NVAR_TOT; ivar++ )
      { SEARCHEFF_TMPMAP2D[ivar] = (double*)malloc(MEMROW);  }
  }

  else if ( strcmp(key,"REQUIRE_DETECTION:")==0 )   { 
    readint(fp, 1, &SEARCHEFF_PHOTPROB[imap].REQUIRE_DETECTION );
  }
  else if ( strcmp(key,"CUTVAL:")==0 )   { 
    readdouble(fp, 1, &SEARCHEFF_PHOTPROB[imap].CUTVAL );
  }
  else if ( strcmp(key,"REDUCED_CORR:")==0 )   { 
    readdouble(fp, 1, &SEARCHEFF_PHOTPROB[imap].REDUCED_CORR );
    INPUTS_SEARCHEFF.NREDUCED_CORR_PHOTPROB++ ;
  }
  else if ( strcmp(key,"NDUMP:")==0 )   { 
    readint(fp, 1, &INPUTS_SEARCHEFF.NPHOTPROB_DUMP );
  }
  else if ( strcmp(key,"VARNAMES:")==0 ) {
    NVAR_READ = 0 ;  VARNAME[0] = 0 ;
    while ( strcmp(VARNAME,"PHOTPROB_WGTLIST") != 0 ) {
      readchar(fp, VARNAME);
      IVARABS = IVARABS_SEARCHEFF_PHOTPROB(VARNAME);
      if ( IVARABS < 0 ) {
	print_preAbort_banner(fnam);
	printf("   VALID PHOTPROB MAP VARIABLES: \n");
	for(ivar=0; ivar < MXDEF_VARNAMES_PHOTPROB; ivar++ ) {
	  printf("\t %s \n", VARDEF_SEARCHEFF_PHOTPROB[ivar]);
	}

	MAPNAME       = SEARCHEFF_PHOTPROB[imap].NAME ;
	sprintf(c1err,"Invalid VARNAME=%s", VARNAME);
	sprintf(c2err,"Check PHOTPROB map = '%s' ", MAPNAME);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }
      sprintf(SEARCHEFF_PHOTPROB[imap].VARNAMES[NVAR_READ],"%s",VARNAME);
      SEARCHEFF_PHOTPROB[imap].IVARABS[NVAR_READ] = IVARABS ;
      SEARCHEFF_PHOTPROB[imap].VALMIN[NVAR_READ] = +1.0E9 ;
      SEARCHEFF_PHOTPROB[imap].VALMAX[NVAR_READ] = -1.0E9 ;
      if ( LDMP ) {
	printf(" xxx %s: read VARNAME(%d) = '%s'\n", fnam,NVAR_READ,VARNAME);
      }
      NVAR_READ++ ;
      if ( NVAR_READ >= MXVAR_SEARCHEFF_PHOTPROB ) {
	sprintf(c1err,"NVAR_READ=%d exceeds bound.", NVAR_READ);
	sprintf(c2err,"Make sure PHOTPROB_WGTLIST is in VARNAMES list.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }	
    } // end while
    NVAR_TOT = SEARCHEFF_PHOTPROB[imap].NVAR_TOT;
    NVAR_MAP = NVAR_READ - 1;             // NVAR of multi-dim map
    NFUN     = NVAR_TOT  - NVAR_MAP; // number of photprob wgts
    SEARCHEFF_PHOTPROB[imap].NVAR_MAP      = NVAR_MAP ;
    SEARCHEFF_PHOTPROB[imap].NFUN_CDF      = NFUN;

    // compute PHOTPROB_CDFBINS; there are 1+NFUN bins including zero.
    double binSize = 1.0 / (double)NFUN ;
    SEARCHEFF_PHOTPROB[imap].BINSIZE = binSize;
    for(ibin=0; ibin <= NFUN; ibin++ ) {
      SEARCHEFF_PHOTPROB[imap].PHOTPROB_CDFBINS[ibin] = 
	binSize*((double)ibin) ;
    }

    if ( LDMP ) {
      printf(" xxx %s: imap=%d load NVAR(TOT,MAP,FUN) = %d, %d, %d \n",
	     fnam, imap, NVAR_TOT, NVAR_MAP, NFUN);
    }

  } // end VARNAMES

  else if ( strcmp(key,"ROW:")==0 ) {
    NVAR_TOT      = SEARCHEFF_PHOTPROB[imap].NVAR_TOT ; 
    NVAR_MAP      = SEARCHEFF_PHOTPROB[imap].NVAR_MAP ; 
    NFUN          = SEARCHEFF_PHOTPROB[imap].NFUN_CDF ; 
    NROW          = SEARCHEFF_PHOTPROB[imap].NROW;
    SEARCHEFF_PHOTPROB[imap].NROW++ ;
    readdouble(fp, NVAR_TOT, TMPVAL);       
    if ( NROW >= MXROW_SEARCHEFF_PHOTPROB ) { return(1); }

    LOAD_PHOTPROB_CDF( NFUN, &TMPVAL[NVAR_MAP] ); // returns CDF

    for(ivar=0; ivar < NVAR_TOT; ivar++ ) { 
      tmpVal = TMPVAL[ivar]; 
      SEARCHEFF_TMPMAP2D[ivar][NROW] = tmpVal ; 

      // keep track of min & max value for each map variable
      if ( ivar < NVAR_MAP ) {
	if ( tmpVal < SEARCHEFF_PHOTPROB[imap].VALMIN[ivar] ) 
	  { SEARCHEFF_PHOTPROB[imap].VALMIN[ivar] = tmpVal; }
	if ( tmpVal > SEARCHEFF_PHOTPROB[imap].VALMAX[ivar] ) 
	  { SEARCHEFF_PHOTPROB[imap].VALMAX[ivar] = tmpVal ; }
      }

    }  // end ivar loop
  }  // end ROW:


  else if ( strcmp(key,"ENDMAP:")==0 ) { 
    IDMAP          = IDGRIDMAP_PHOTPROB_OFFSET + imap ;
    NROW           = SEARCHEFF_PHOTPROB[imap].NROW ;
    MAPNAME        = SEARCHEFF_PHOTPROB[imap].NAME ;
    NVAR_TOT       = SEARCHEFF_PHOTPROB[imap].NVAR_TOT ;
    NVAR_MAP       = SEARCHEFF_PHOTPROB[imap].NVAR_MAP ;
    NFUN           = SEARCHEFF_PHOTPROB[imap].NFUN_CDF ;
    if ( NROW >= MXROW_SEARCHEFF_PHOTPROB ) {
      sprintf(c1err,"NROW=%d exceeds bound of MXROW_SEARCHEFF_PHOTPROB=%d",
	      NROW, MXROW_SEARCHEFF_PHOTPROB ) ;
      sprintf(c2err,"Check MAPNAME_PHOTROB: %s", MAPNAME);	      
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    if ( LDMP ) {
      printf(" xxx ENDMAP(%s): imap=%d  IDMAP=%d  "
	     "NVAR(TOT,MAP,FUN)=%d,%d,%d  NROW=%d \n",
	     MAPNAME, imap, IDMAP, 
	     NVAR_TOT, NVAR_MAP, NFUN, NROW) ; 
      fflush(stdout);
    }

    init_interp_GRIDMAP(IDMAP, MAPNAME, NROW, NVAR_MAP, NFUN, 0, 
			SEARCHEFF_TMPMAP2D, 
			&SEARCHEFF_TMPMAP2D[NVAR_MAP],
			&SEARCHEFF_PHOTPROB[imap].GRIDMAP  );  // <== returned


    for(ivar=0; ivar<NVAR_TOT; ivar++ ) { free(SEARCHEFF_TMPMAP2D[ivar]); }
    free(SEARCHEFF_TMPMAP2D);

    MAPVERSION_SEARCHEFF_PHOTPROB = 0 ; 
    return(0) ; 
  }
  else {
    return(1) ; 
  }

  return(0);

} // end readMap_SEARCHEFF_PHOTPROB


// *******************************************
int IVARABS_SEARCHEFF_PHOTPROB(char *VARNAME) {

  // Return absolute IVAR index for VARNAME
  // Return -9 if not defined.

  int ivar, IVARABS = -9;
  char *varTmp ;
  //  char fnam[] = "IVARABS_SEARCHEFF_PHOTPROB" ;

  for(ivar=0; ivar < MXDEF_VARNAMES_PHOTPROB; ivar++ ) {
    varTmp = VARDEF_SEARCHEFF_PHOTPROB[ivar];
    if ( strcmp(VARNAME,varTmp) == 0 ) { IVARABS = ivar; }
  }

  if ( strcmp(VARNAME,"PHOTPROB_WGTLIST") == 0 ) { IVARABS=99; }

  return(IVARABS);

} // end  IVARABS_SEARCHEFF_PHOTPROB


// *********************************************
void check_SEARCHEFF_DETECT(int imap) {

  // Created Jul 6, 2011 
  // Check structure SEARCHEFF_PIPELINE[ifilt_obs]
  // Load text for README file containing SNR/MAG
  // at which EFF = 1 and 0.5.
  //
  // Mar 7 2018: 
  //  + rename check_SEARCHEFF_FILTER -> check_SEARCHEFF_DETECT.
  //  + pass imap arg instead of ifilt_obs
  //  + float -> double
  //
  // Mar 13 2018:
  //  + fix bug looping ibin=1 to NBIN; now loops 0 to NBIN-1

  int NBIN    = SEARCHEFF_DETECT[imap].NBIN ;
  char *cfilt = SEARCHEFF_DETECT[imap].FILTERLIST ;

  int ibin, i, IBIN_HALF, IBIN_ONE ;
  double  VAL, VAL_LAST, EFF, EFFDIF, EFFDIF_ONE, EFFDIF_HALF ;

  char 
    *ptr_effname 
    ,cline[MXPATHLEN]
    ,fnam[] = "check_SEARCHEFF_DETECT" 
    ;

  // ------------ BEGIN -------------

  if ( imap==0 ) 
    { printf("   Read Detection-efficiency curves: \n"); }

  if ( NBIN >= MXROW_SEARCHEFF_DETECT ) {
    sprintf(c1err,"%d DETECT-SNR bins exceeds MXBIN_DETECT=%d", 
	    NBIN, MXROW_SEARCHEFF_DETECT);
    sprintf(c2err,"Check filter = %s", cfilt );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  ptr_effname    = SEARCHEFF_PARNAME[SEARCHEFF_FLAG];

  SEARCHEFF_DETECT[imap].NLINE_README = 0 ; 


  VAL_LAST   = 0.0 ;
  EFFDIF_ONE = EFFDIF_HALF = 999.0;
  IBIN_HALF  = IBIN_ONE = 0;

  for ( ibin=0; ibin < NBIN; ibin++ ) {
    //  for ( ibin=1; ibin <= NBIN; ibin++ ) {

    VAL =  SEARCHEFF_DETECT[imap].VAL[ibin] ;
    EFF =  SEARCHEFF_DETECT[imap].EFF[ibin] ;

    //    printf(" xxx bin=%d  VAL=%f  EFF=%f \n", ibin, VAL, EFF);

    // keep track of bin in which efficiency is closest to 1 or .5
    EFFDIF = fabs(1.0-EFF) ;
    if ( EFFDIF < EFFDIF_ONE ) 
      { EFFDIF_ONE = EFFDIF ; IBIN_ONE = ibin ; }
      
    EFFDIF = fabs(0.5-EFF) ;
    if ( EFFDIF < EFFDIF_HALF ) 
      { EFFDIF_HALF = EFFDIF ; IBIN_HALF = ibin ; }    

  } // NBIN 

  // ------
  // print VAL (SNR or MAG) when eff is closest to  1.
  
  //      printf(" xxx IBIN[ONE,HALF] = %d  %d \n", IBIN_ONE, IBIN_HALF);
  VAL = SEARCHEFF_DETECT[imap].VAL[IBIN_ONE];
  EFF = SEARCHEFF_DETECT[imap].EFF[IBIN_ONE];
  sprintf(cline, "\t Epoch SEARCH_EFF(%s) = %5.2f at %s = %5.2f ", 
	  cfilt, EFF, ptr_effname, VAL );
  
  printf("%s", cline);
  i = SEARCHEFF_DETECT[imap].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[imap].README[i], "%s", cline);
  SEARCHEFF_DETECT[imap].NLINE_README++ ;

  // ------
  // print VAL (SNR or MAG) when eff is closest to  0.5

  VAL = SEARCHEFF_DETECT[imap].VAL[IBIN_HALF];
  EFF = SEARCHEFF_DETECT[imap].EFF[IBIN_HALF];
  sprintf(cline, "\t Epoch SEARCH_EFF(%s) = %5.2f at %s = %5.2f ", 
	  cfilt, EFF, ptr_effname, VAL );

  printf("%s", cline);
  i = SEARCHEFF_DETECT[imap].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[imap].README[i], "%s", cline);
  SEARCHEFF_DETECT[imap].NLINE_README++ ;

} // end of check_SEARCHEFF_DETECT


// *********************************************
void check_SEARCHEFF_PHOTPROB(int imap) {
  
  // set README comment(s).
  //  char fnam[] = "check_SEARCHEFF_PHOTPROB" ;

  // -------------- BEGIN ----------------

  if ( imap==0 ) {
    printf("\n");
    printf("        PHOTPROB \n");
    printf("         MAPNAME             FIELD(s)     BAND(s)   "
	   "NFUN   NROW \n");
    printf("# -------------------------------"
	   "------------------------------ \n");
  }

  char *NAME      = SEARCHEFF_PHOTPROB[imap].NAME;
  char *FIELDLIST = SEARCHEFF_PHOTPROB[imap].FIELDLIST;
  char *BANDLIST  = SEARCHEFF_PHOTPROB[imap].FILTERLIST;
  int  NFUN       = SEARCHEFF_PHOTPROB[imap].NFUN_CDF ;
  int  NROW       = SEARCHEFF_PHOTPROB[imap].NROW;

  printf("  %d) %16.16s   %10.10s   %10.10s       %d   %4d \n",
	 imap, NAME, FIELDLIST, BANDLIST, NFUN, NROW );
  fflush(stdout) ;

  SEARCHEFF_PHOTPROB[imap].NLINE_README = 0 ; 
  int i = SEARCHEFF_PHOTPROB[imap].NLINE_README ;
  char *cptr = SEARCHEFF_PHOTPROB[imap].README[i] ;
  sprintf(cptr, "  PHOTPROB MAP %s: BANDS=%s  FIELD=%s \n",
	  NAME, BANDLIST, FIELDLIST );
  SEARCHEFF_PHOTPROB[imap].NLINE_README++ ;

  return ;

} // end check_SEARCHEFF_PHOTPROB

// *******************************************
void  init_SEARCHEFF_LOGIC(char *survey) {

  // Nov 29 2014: check local dir first, then official area
  //               use snana_openTextFile utility 
  //
  // Dec 27 2015: check for default logic file name, or user-input
  //
  // Aug 26 2020: pass OPTMASK to snana_open to check DOCANA

  int OPTMASK = INPUTS_SEARCHEFF.OPTMASK_OPENFILE ;
  int NMJD, i, gzipFlag ;

  FILE *fp ;

  char 
     cline[MXPATHLEN]
    ,logic[60]
    ,c_get[60]
    ,surveykey[60]
    ,logicFile_Default[] = "SEARCHEFF_PIPELINE_LOGIC.DAT"
    ,logicFile[MXPATHLEN]
    ,*ptrFile_user
    ,*ptrFile_final
    ,fnam[] = "init_SEARCHEFF_LOGIC" 
    ;

  // -------------- BEGIN ----------------


  ptrFile_user  = INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE ;
  ptrFile_final = INPUTS_SEARCHEFF.PIPELINE_LOGIC_FILE ;

  if ( strcmp(ptrFile_user,"DEFAULT") == 0 ) 
    { sprintf(logicFile,"%s", logicFile_Default); }
  else
    { sprintf(logicFile,"%s", ptrFile_user);  }


  //  printf(" xxx logic file -> '%s' \n", logicFile); 

  fp = snana_openTextFile(OPTMASK, PATH_SEARCHEFF, logicFile,
			  ptrFile_final, &gzipFlag ); // returned

  if ( !fp ) {
    abort_openTextFile("SEARCHEFF_PIPELINE_LOGIC_FILE",
		       PATH_SEARCHEFF, logicFile, fnam );
  }


  sprintf(cline, "\n   Fetch SOFTWARE SEARCH-LOGIC from : "); 
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", cline);
  printf("%s", cline);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;


  sprintf(cline, "\t %s", ptrFile_final ); 
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", cline);
  printf("%s", cline);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;

  
  NMJD = 0;
  sprintf(surveykey, "%s:", survey);
  while( (fscanf(fp, "%s", c_get )) != EOF) {
    if ( strcmp(surveykey,c_get) == 0 ) {
      readint(fp, 1, &NMJD);
      readchar(fp, logic ) ;
    }
  }
  fclose(fp);


  if ( NMJD == 0 ) {
    sprintf(c1err,"Could NOT find search pipeline logic for %s", survey );
    sprintf(c2err,"Check %s", ptrFile_final ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }
  
  parse_search_eff_logic(survey,NMJD,logic);



} // end of  init_SEARCHEFF_LOGIC


// *******************************************
void parse_search_eff_logic(char *survey, int NMJD, char *logic) {

  // Created may 22, 2008 by R.Kessler
  // Parse logic string;
  // example: 'gr+ri+gi' => require two of three gri filters
  //
  // Jan 13 2017: replace fortran function ifiltindx_ with INTFILTER()

  int len, i, ifiltdef, NOR, NAND ;
  char ctmp[4];
  char fnam[] = "parse_search_eff_logic" ;

  // ------------- BEGIN ----------------

  // store function args in global structure
  SEARCHEFF_LOGIC.NMJD = NMJD ;
  sprintf(SEARCHEFF_LOGIC.INPUT_STRING, "%s", logic);

  NAND  = SEARCHEFF_LOGIC.NMASK = 0;
  NOR   = 1;

  // init logic array to zero.
  for ( i=0; i< MXMASK_SEARCHEFF_LOGIC; i++ ) {
    SEARCHEFF_LOGIC.IFILTDEF_MASK[i] = 0;
  }

  
  // use c1err for message, even though it's not an error
  sprintf(c1err, "\t Logic: %d MJDs require filters=%s ", NMJD, logic);
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", c1err);
  printf("%s", c1err);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;

  sprintf(c1err, "\t Trigger epoch contains all obs withing %.3f days",
	  INPUTS_SEARCHEFF.TIME_SINGLE_DETECT);
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", c1err);
  printf("%s", c1err);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;

  len = strlen(logic);

  for ( i=0; i<len; i++ ) {
    sprintf(ctmp, "%c", *(logic+i) ) ;

    if ( strcmp(ctmp,"+")==0 )
      { ifiltdef = 0;   NOR++;   NAND=0; }
    else {

      ifiltdef = INTFILTER(ctmp) ;
      NAND++ ;
      if ( ifiltdef <= 0 ) {
	sprintf(c1err,"Invalid filter='%s' in search-eff logic.", ctmp);
	sprintf(c2err,"Check %s: %s", survey, logic );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
      }
    }

    if ( NAND > 0 ) 
      { SEARCHEFF_LOGIC.IFILTDEF_MASK[NOR-1] |= (1 << ifiltdef); }

    /*
    printf("\t logic-char(%d) = %s => MASK[%d] =%d \n", 
	   i, ctmp, NOR, SEARCHEFF_LOGIC.IFILTDEF_MASK[NOR] );
    */

  } // end of char-loop (i)


  if ( NOR >= MXMASK_SEARCHEFF_LOGIC ) {
    print_preAbort_banner(fnam);
    printf("   PIPELINE_LOGIC_FILE: %s\n", 
	   INPUTS_SEARCHEFF.PIPELINE_LOGIC_FILE);
    printf("   LOGIC STRING: '%s' ", logic);
    sprintf(c1err,"NOR=%d exceeds bound of %d", NOR, MXMASK_SEARCHEFF_LOGIC);
    sprintf(c2err,"Check LOGIC file and string above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  printf("\n"); fflush(stdout);

  SEARCHEFF_LOGIC.NMASK = NOR;

} // end of parse_search_eff_logic

// *******************************************
void  init_SEARCHEFF_SPEC(char *survey) {

  // Created July 2011 by R.Kessler
  // If optional EFF(SPEC-confirm) table exists, read and store.
  // Table name must be SEARCHEFF_[survey]_SPEC.DAT;
  // this table can reside in current [user] directory
  // or in PATH_SEARCHEFF
  //
  // Table VARNAMES include
  // * REDSHIFT
  // * any of the filters (i.e., 'g' => peak g-band mag).
  // * SPECEFF   (required last varname)
  //
  //
  // If user-requested file (SEARCHEFF_SPEC_FILE: <file>)
  // does not exist then abort. If default SPEC-eff file
  // does not exist then just return.
  //
  // Aug 26, 2011 fix aweful bug
  //     MAPSIZE    = SEARCHEFF_SPEC.MAPSIZE[NMAP] ;    
  //  -> MAPSIZE    = SEARCHEFF_SPEC.MAPSIZE[imap] ; 
  //
  // Jun 23 2016: check for non-uniform binning.
  //
  // Nov  6 2017: set SEARCHEFF_SPEC.FLAG_PEAKMAG_ONLY
  // Mar 11 2018: refactor
  // Jun 02 2018: check SPEC_FILE = ZERO option
  // Mar 14 2019: refactor to use read_GRIDMAP().
  //

  char 
     effspec_file_local[MXPATHLEN]
    ,*ptrFile_user, *ptrFile_final, *cptr
    ,c_get[60], LINE[100]
    ,fnam[] = "init_SEARCHEFF_SPEC" 
    ;

  FILE *fp ;

  int   NMAP, imap, NROW, ivar, NVAR, FOUND_VARNAMES ;
  int   ID, NDIM, NFUN, N, IREQUIRE, gzipFlag, ISFIELD ;
  char  KEY_ROW[]  = "SPECEFF:";
  char  KEY_STOP[] = "" ;
  char  *VARNAME, *VARLIST ; 
  char  FIELDLIST[100] ;

  int    OPTMASK       = INPUTS_SEARCHEFF.OPTMASK_OPENFILE ;
  double SPECEFF_SCALE = INPUTS_SEARCHEFF.USER_SPECEFF_SCALE ;

  // ------------ BEGIN --------------

  N = 0 ;  cptr = SEARCHEFF_SPEC_INFO.README[N] ;
  sprintf(cptr, "\t %s", "No spec-eff option specified ==> 100% efficiency.");
  N++;  SEARCHEFF_SPEC_INFO.NLINE_README = N ;  

  INPUTS_SEARCHEFF.IFLAG_SPEC_EFFZERO=0;
  SEARCHEFF_SPEC_INFO.IVARTYPE_MASK = 0 ;
  SEARCHEFF_SPEC_INFO.FLAG_PEAKMAG_ONLY = 0 ;
  SEARCHEFF_SPEC_INFO.BOOLEAN_OR        = 1 ;
  SEARCHEFF_SPEC_INFO.BOOLEAN_AND       = 0 ; // default logic is AND

  for( imap=0; imap < MXMAP_SEARCHEFF_SPEC; imap++ ) {
    SEARCHEFF_SPEC[imap].IVAR           = -9 ;
    SEARCHEFF_SPEC[imap].IVAR_REDSHIFT  = -9 ;
    SEARCHEFF_SPEC[imap].IVAR_PEAKMJD   = -9 ;
    SEARCHEFF_SPEC[imap].IVAR_DTPEAK    = -9 ;
    SEARCHEFF_SPEC[imap].IVAR_SALT2mB   = -9 ;
    sprintf(SEARCHEFF_SPEC[imap].FIELDLIST, "ALL" ); // default is all fields

    for ( ivar=0 ; ivar < MXVAR_SEARCHEFF_SPEC; ivar++ ) {
      SEARCHEFF_SPEC[imap].IVARTYPE[ivar]              = -9 ;
      SEARCHEFF_SPEC[imap].NFILTLIST_PEAKMAG[ivar]     =  0 ;
      SEARCHEFF_SPEC[imap].IFILTOBS_HOSTMAG[ivar]      = -9 ;
      SEARCHEFF_SPEC[imap].IFILTOBS_SBMAG[ivar]        = -9 ;
      SEARCHEFF_SPEC[imap].IFILTOBS_PEAKCOLOR[ivar][0] = -9 ;
      SEARCHEFF_SPEC[imap].IFILTOBS_PEAKCOLOR[ivar][1] = -9 ;  
    }
  } // imap


  // define effspec filename form user input or from default name

  ptrFile_user  = INPUTS_SEARCHEFF.USER_SPEC_FILE ;
  ptrFile_final = INPUTS_SEARCHEFF.SPEC_FILE ;


  if ( IGNOREFILE(ptrFile_user) ) {
    sprintf(effspec_file_local, "%s",  ptrFile_user); 
    IREQUIRE = 0 ;
  }
  else if ( strcmp(ptrFile_user,"ZERO") == 0  ) {
    INPUTS_SEARCHEFF.IFLAG_SPEC_EFFZERO = 1 ;
    sprintf(cptr,"\t SEARCHEFF_SPEC_FILE = ZERO -> Force EFF_SPEC=0. \n");
    printf("\n %s\n", cptr);
    fflush(stdout);
    return ;
  }
  else { 
    sprintf(effspec_file_local, "%s",  ptrFile_user); 
    IREQUIRE = 1 ;
  }


  // use utility to check local dir and path.
  fp = snana_openTextFile(OPTMASK, PATH_SEARCHEFF, effspec_file_local, 
			  ptrFile_final, &gzipFlag ); // returned

  if ( fp == NULL ) { 
    if ( IREQUIRE ) {
      abort_openTextFile("SEARCHEFF_SPEC_FILE",
			 PATH_SEARCHEFF, effspec_file_local, fnam );
    }
    else  { 
      //      printf("\n  Optional SEARCHEFF_SPEC_FILE not specified -> skip. \n");
      printf("\n  Optional SEARCHEFF_SPEC_FILE not specified -> EFF_SPEC=1 \n");
      fflush(stdout);
      NONZERO_SEARCHEFF_SPEC++; // default EFFspec = 100%
      return ; 
    }
  }

  printf("\n");
  printf("   Reading spectroscopic efficiency from \n\t %s\n", 
	 ptrFile_final );

  NVAR = NROW = NMAP = 0;   sprintf(FIELDLIST,"ALL");

  sprintf(c2err,"%s", "Check spec-eff file above");

  while( (fscanf(fp, "%s", c_get )) != EOF) {

    if ( strcmp(c_get,"NVAR:")==0 ) { warn_NVAR_KEY(ptrFile_final); }

    // - - - - - -
    // check optional key to associate map with particular field[list]
    ISFIELD = (strcmp(c_get,"FIELD:")==0 || strcmp(c_get,"FIELDLIST:")==0 );
    if ( ISFIELD ) { readchar(fp, FIELDLIST); }

    FOUND_VARNAMES = ( strcmp(c_get,"VARNAMES:")==0 ) ;

    if ( FOUND_VARNAMES && NMAP < MXMAP_SEARCHEFF_SPEC ) {
      fgets(LINE, 100, fp ); // scoop up varnames
      NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
      sprintf(SEARCHEFF_SPEC[NMAP].FIELDLIST,"%s", FIELDLIST );

      for ( ivar=0; ivar < NVAR; ivar++ ) {
	VARNAME = SEARCHEFF_SPEC[NMAP].VARNAMES[ivar] ;
	VARLIST = SEARCHEFF_SPEC[NMAP].GRIDMAP.VARLIST ;
	get_PARSE_WORD(0,ivar,VARNAME);
	assign_SPECEFF(NMAP,ivar,VARNAME); 
	if ( ivar == 0 ) 
	  { sprintf(VARLIST,"%s", VARNAME ); }
	else
	  { sprintf(VARLIST,"%s %s", VARLIST, VARNAME ); }
      }
      
      ID = IDGRIDMAP_SPECEFF_OFFSET + NMAP;   NDIM = NVAR-1; NFUN=1;
      read_GRIDMAP(fp, KEY_ROW, KEY_ROW, KEY_STOP, ID, NDIM, NFUN, 0,
		   MXROW_SEARCHEFF_SPEC, fnam,
		   &SEARCHEFF_SPEC[NMAP].GRIDMAP ) ;
      printf("\t for FIELDLIST=%s\n", FIELDLIST);
      NONZERO_SEARCHEFF_SPEC++; 
    }

    if ( FOUND_VARNAMES ) { NMAP++; }

  }  // end of read-line loop

 
  fclose(fp); // done reading SPECEFF map(s)

  INPUTS_SEARCHEFF.NMAP_SPEC = NMAP;

  if ( NMAP > MXMAP_SEARCHEFF_SPEC ) {
    sprintf(c1err,"NMAP=%d  exceeds  MXMAP_SEARCHEFF_SPEC=%d",
	    NMAP, MXMAP_SEARCHEFF_SPEC );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  // -------------------------------------------------------
  // Nov 2017: check if this map depends ONLY on PEAKMAG;
  //   if so, can make sim run faster.
  int IVARTYPE_MASK_PEAKMAG = 0, OVP;
  IVARTYPE_MASK_PEAKMAG += (1<< IVARTYPE_SPECEFF_PEAKMAG ) ;
  IVARTYPE_MASK_PEAKMAG += (1<< IVARTYPE_SPECEFF_COLOR   ) ;
  OVP = ( IVARTYPE_MASK_PEAKMAG & SEARCHEFF_SPEC_INFO.IVARTYPE_MASK ) ;
  if ( SEARCHEFF_SPEC_INFO.IVARTYPE_MASK == OVP ) 
    { SEARCHEFF_SPEC_INFO.FLAG_PEAKMAG_ONLY=1 ; }

  /* xxx
  printf(" xxx FLAG_PEAKMAG_ONLY = %d (OVP=%d, MASKREF=%d, MASK=%d) \n", 
	 SEARCHEFF_SPEC.FLAG_PEAKMAG_ONLY, 
	 OVP, IVARTYPE_MASK_PEAKMAG, SEARCHEFF_SPEC.IVARTYPE_MASK );
  */

  fflush(stdout) ;

  N = 0 ;
  cptr = SEARCHEFF_SPEC_INFO.README[N] ;  N++ ;  
  sprintf(cptr, "\t %s", "GRID-MAP read from ");
  cptr = SEARCHEFF_SPEC_INFO.README[N] ;  N++ ;  
  sprintf(cptr, "\t %s", ptrFile_final );

  if ( fabs(SPECEFF_SCALE-1.0) > 0.0001 ) {
    cptr = SEARCHEFF_SPEC_INFO.README[N] ;    N++ ;  
    sprintf(cptr, "\t EffSpec from file is scaled by %le", SPECEFF_SCALE );
  }

  SEARCHEFF_SPEC_INFO.NLINE_README = N ;

  return ;

} // end of init_SEARCHEFF_SPEC


// ***************************************
void init_SEARCHEFF_zHOST(char *survey) {

  // Created May 8 2014 by R.Kessler
  //
  // Read optional efficiency map for obtaining a host-galaxy
  // redshift for UNCONFIRMED SNe. Note that non-uniform 
  // redshift binning is allowed.
  //
  // Map is of the form
  // FIELD:    <bla>             # optional to give separate map per field
  // FIELDLIST <bla1+bla2+bla3>  # optional OR of several fields.
  // HOSTEFF:  <z>  <eff>        # required
  // HOSTEFF:  <z>  <eff>        # idem
  // HOSTEFF:  <z>  <eff>
  // etc ...
  //
  // Apr 11 2016: if no HOSTEFF keys, then read as 2-column file
  //              (same format as HOSTLIB_ZPHOTEFF_FILE)
  //
  // Mar 11 2018: refactor with new struct SEARCHEFF_zHOST
  // Jul 19 2018: read optional CUTWIN_SNRMAX 
  // Jul 30 2018: set NONZERO_SEARCHEFF_zHOST when reading 2-column format

  FILE *fp ;
  char fnam[] = "init_SEARCHEFF_zHOST" ;

  // --------------- BEGIN ------------

  // check file format version: LEGACY or MULTI-D
  // If init_HOSTLIB already checked, then 
  // INPUTS_SEARCHEFF.IVERSION_zHOST is already set > 0.
  if ( INPUTS_SEARCHEFF.IVERSION_zHOST == 0 ) {
    fp = open_zHOST_FILE(1);
    if ( fp != NULL ) { read_VARNAMES_zHOST(fp); fclose(fp); }
  }

  // ------------------------------------
  
  // check VARNAMES key to see if legacy (z-only) file
  int LEGACY = ( INPUTS_SEARCHEFF.IVERSION_zHOST == IVERSION_zHOST_LEGACY );

  fp = open_zHOST_FILE(1);
  if ( fp == NULL ) { return ; }

  if ( LEGACY ) 
    { read_zHOST_FILE_LEGACY(fp); }
  else
    { read_zHOST_FILE(fp); }

  fclose(fp);

  return;
} // end of init_SEARCHEFF_zHOSTEFF


// ************************************
FILE *open_zHOST_FILE(int OPT) {

  // Mar 20 2019
  // Open zHOST file and return file pointer 
  //
  // OPT = +1 -> normal init called from trigger code
  // OPT = -1 -> called from init_HOSTLIB to get VARNAMES

  int LPRINT = ( OPT > 0 ) ; // stdout printing
  int IREQUIRE, gzipFlag ;
  char *ptrFile_user ;
  char *ptrFile_final ;
  char localFile[MXPATHLEN];
  FILE *fp ;
  char fnam[] = "open_zHOST_FILE" ;

  // --------------- BEGIN ------------------

  fp = NULL ;

  INPUTS_SEARCHEFF.IFLAG_zHOST_EFFZERO = 0;
  INPUTS_SEARCHEFF.NMAP_zHOST = 0 ;

  ptrFile_user  = INPUTS_SEARCHEFF.USER_zHOST_FILE ; // from sim-input
  ptrFile_final = INPUTS_SEARCHEFF.zHOST_FILE ;      // includes full path
  
  if( IGNOREFILE(ptrFile_user) ) {
    // NULL or NONE, etc ... -> EFF=1
    if ( LPRINT )   { 
      printf("\n  Optional SEARCHEFF_zHOST_FILE not specified "
	     "-> Eff=1.0 \n"); 
    }
    return(fp);
  }
  else if ( strcmp(ptrFile_user,"ZERO") == 0 ) {
    // no file to read, but set all EFF_zHOST=0  (Jun 2018)
    INPUTS_SEARCHEFF.IFLAG_zHOST_EFFZERO = 1;
    if ( LPRINT ) 
      { printf("\n  SEARCHEFF_zHOST_FILE=ZERO -> Force EFF_zHOST=0 \n"); }
    return(fp) ;
  }
  else { 
    sprintf(localFile, "%s",  ptrFile_user); 
    IREQUIRE = 1 ;
  }


  // use utility to check local dir and path.
  int OPTMASK_OPEN  = INPUTS_SEARCHEFF.OPTMASK_OPENFILE ;
  fp = snana_openTextFile(OPTMASK_OPEN, PATH_SEARCHEFF, localFile, 
			  ptrFile_final, &gzipFlag); // returned
  
  // examine if there is no zHOST file
  if ( fp == NULL ) { 
    if ( IREQUIRE ) {
      abort_openTextFile("SEARCHEFF_zHOST_FILE",
			 PATH_SEARCHEFF, localFile, fnam );
    }
  }

  // --------------------------------------------

  if ( LPRINT ) {
    printf("\n");
    printf("   Read zHOST-for-unconfirmed effic map from \n\t %s\n", 
	   ptrFile_final );
    fflush(stdout) ;
  }

  return(fp);

} // end open_zHOST_FILE

// ***************************************
void read_zHOST_FILE(FILE *fp) {

  // March 2019:
  // Read zHOST effic file with arbitrary dependence on HOSTLIB properties.
  // If using FIELDLIST, each map must start with a new VARNAMES key 
  // 
  // Dec 3 2019: fix bug by setting KEY_STOP = ""
  // Jul 13 2020: read optional PEAKMJD

  int  OPT_EXTRAP = 0 ;
  int  NTAB=0;
  int  IDMAP, imap, NMAP, ivar, NVAR, NDIM, NFUN ;
  int  FOUND_VARNAMES;
  char c_get[60], FIELDLIST[100] ;
  char *ptr_VARNAMES[MXVAR_SEARCHEFF_zHOST], *VARLIST;
  char VARNAME_HOSTLIB_TMP[MXVAR_SEARCHEFF_zHOST][40];
  int  IVAR_HOSTLIB_TMP[MXVAR_SEARCHEFF_zHOST];
  double PEAKMJD_RANGE[2];
  char KEY_ROW[]   = "HOSTEFF:" ;
  char KEY_STOP[]  = "" ;
  char fnam[] = "read_zHOST_FILE" ;

  // ------------ BEGIN ----------

  for(imap=0; imap < MXMAP_SEARCHEFF_zHOST; imap++ ) {
    sprintf(SEARCHEFF_zHOST[imap].FIELDLIST,"NONE" );
  }
  for(ivar=0; ivar < MXVAR_SEARCHEFF_zHOST; ivar++ ) {
    ptr_VARNAMES[ivar] = VARNAME_HOSTLIB_TMP[ivar]; 
  }
  NMAP = NVAR = 0;    
  sprintf(FIELDLIST,"ALL");
  PEAKMJD_RANGE[0] = 10000.0;
  PEAKMJD_RANGE[1] = 90000.0;

  // -----------------------------
  while( (fscanf(fp, "%s", c_get )) != EOF) {

    if ( c_get[0] == '\t' ) { NTAB++; }
    
    if ( strcmp(c_get,"NVAR:")==0 ) 
      { warn_NVAR_KEY(INPUTS_SEARCHEFF.zHOST_FILE); }

    FOUND_VARNAMES = 0 ;

    if ( strcmp(c_get,"OPT_EXTRAP:") == 0 ) 
      { readint(fp,1,&OPT_EXTRAP); }
    if ( strcmp(c_get,"FIELDLIST:" ) == 0 ) 
      { readchar(fp,FIELDLIST); }
    if ( strcmp(c_get,"PEAKMJD:") == 0 || strcmp(c_get,"PEAKMJD_RANGE:")==0 ) 
      { readdouble(fp,2,PEAKMJD_RANGE);}

    if ( strcmp(c_get,"VARNAMES:"  ) == 0 ) { FOUND_VARNAMES=1; }

    // parse VARNAMES in map
    if ( FOUND_VARNAMES && NMAP < MXMAP_SEARCHEFF_zHOST ) {
      VARLIST = SEARCHEFF_zHOST[NMAP].GRIDMAP.VARLIST ;
      NVAR = parse_VARNAMES_zHOST(fp, IVAR_HOSTLIB_TMP, ptr_VARNAMES,VARLIST);
      NDIM = NVAR-1; NFUN=1;  IDMAP=IDGRIDMAP_zHOST_OFFSET+NMAP ;

      // store a few things in global
      sprintf(SEARCHEFF_zHOST[NMAP].FIELDLIST,"%s", FIELDLIST);
      SEARCHEFF_zHOST[NMAP].PEAKMJD_RANGE[0] = PEAKMJD_RANGE[0];
      SEARCHEFF_zHOST[NMAP].PEAKMJD_RANGE[1] = PEAKMJD_RANGE[1];

      for(ivar=0; ivar<NVAR; ivar++ ) {
	SEARCHEFF_zHOST[NMAP].IVAR_HOSTLIB[ivar] = IVAR_HOSTLIB_TMP[ivar];
	sprintf(SEARCHEFF_zHOST[NMAP].VARNAMES_HOSTLIB[ivar],"%s",
		ptr_VARNAMES[ivar] ) ;
      }

      read_GRIDMAP(fp,KEY_ROW,KEY_ROW, KEY_STOP, IDMAP, NDIM, NFUN,OPT_EXTRAP,
		   MXROW_SEARCHEFF_zHOST, fnam,
		   &SEARCHEFF_zHOST[NMAP].GRIDMAP ) ;

      printf("\t FIELDLIST = %s \n", FIELDLIST); 
      if ( PEAKMJD_RANGE[0] > 10001.0 ) {
	printf("\t PEAKMJD_RANGE = %d to %d \n", 
	       (int)PEAKMJD_RANGE[0], (int)PEAKMJD_RANGE[1] );
      }
      fflush(stdout);

      NONZERO_SEARCHEFF_zHOST++ ;
    } // end VARNAMES check

    // always increment NMAP, even if past bound.
    if ( FOUND_VARNAMES ) { NMAP++ ; }
   
  } // end while

  INPUTS_SEARCHEFF.NMAP_zHOST = NMAP ;


  if ( NMAP > MXMAP_SEARCHEFF_zHOST ) {
    sprintf(c1err,"NMAP=%d exceeds MXMAP_SEARCHEFF_zHOST=%d",
	    NMAP, MXMAP_SEARCHEFF_zHOST );
    sprintf(c2err,"Too many zHOST maps; check SEARCHEFF_zHOST_FILE arg.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;       
  }

  tabs_ABORT(NTAB, INPUTS_SEARCHEFF.zHOST_FILE, fnam);

  return ;

} // end read_zHOST_FILE


// ***************************************************
int parse_VARNAMES_zHOST(FILE *fp, int *ivar_HOSTLIB, 
			 char **varName_HOSTLIB, char *varNameList) {

  // Mar 2019
  // parse variable names after VARNAMES key; return(NVAR)
  //
  // Inputs:
  //    fp   file pointer to read
  // Outputs:
  //   *ivar_HOSTLIB     : HOSTLIB ivar for each variable
  //  **varName_HOSTLIB  : pointer to each varName
  //   *varNameList      : space-separate list (for printing)
  //

  int NVAR, ivar, *ivar_H;
  char LINE[100], *varName_H;
  char fnam[] = "parse_VARNAMES_zHOST" ;
  int  NERR = 0 ;
  int  LDMP = 0 ;
  // ----------- BEGIN ------------

  fgets(LINE, 100, fp ); // scoop up varnames
  if ( LDMP ) {
    printf(" xxx --------------------------------- \n" );
    printf(" xxx VARNAMES LINE = '%s' \n", LINE);
  }

  varNameList[0] = 0 ;

  NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
  for ( ivar=0; ivar < NVAR; ivar++ ) {
    varName_H = varName_HOSTLIB[ivar];
    ivar_H    = &ivar_HOSTLIB[ivar];

    get_PARSE_WORD(0,ivar,varName_H);
    checkAlternateVarNames_HOSTLIB(varName_H); // 10.03.2020
    *ivar_H  = IVAR_HOSTLIB(varName_H,0);
    
    if (LDMP ) {
      printf(" xxx varName[%d] = '%s' (IVAR_HOSTLIB=%d)\n", 
	     ivar, varName_H , *ivar_H);
    }

    // sanity checks
    if ( ivar < NVAR-1  && *ivar_H < 0 ) {
      NERR++;
      printf("\t ERROR: unknown zHOST variable '%s' not in HOSTLIB\n",
	     varName_H);
    }

    /*
    if ( ivar == NVAR-1 && strcmp(varName_H,"EFF")!=0 ) {
      NERR++ ;
      printf("\t ERROR: last zHOST variable is '%s, but must be EFF\n", 
	     varName_H);   
    }
    */

    if ( ivar==0 ) 
      { sprintf(varNameList,"%s", varName_H); }
    else
      { sprintf(varNameList,"%s %s", varNameList, varName_H); }

  }
  
  if ( NERR > 0 ) {
    sprintf(c1err,"%d zHOST variables are invalid", NERR);
    sprintf(c2err,"See ERROR messages above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;       
  }

  return(NVAR) ;

} // end parse_VARNAMES_zHOST

// ***************************************
void read_zHOST_FILE_LEGACY(FILE *fp) {

  // March 2019
  // This legacy code is separated from init_SEARCHEFF_zHOST.
  // Here the efficiency depends only on redshift.

  char c_get[60], FIELDLIST[100], *ptrField ;
  int  NROW = 0, NMAP=0, imap ;
  int  MEMD = MXROW_SEARCHEFF_zHOST * sizeof(double);
  double VAL[2];
  char fnam[] = "read_zHOST_FILE_LEGACY" ;

  for(imap=0; imap < MXMAP_SEARCHEFF_zHOST; imap++ ) {
    SEARCHEFF_zHOST_LEGACY[imap].NROW = 0 ;
    sprintf(SEARCHEFF_zHOST_LEGACY[imap].FIELDLIST,"NONE" );
  }

  // allocate memory for first map
  SEARCHEFF_zHOST_LEGACY[0].REDSHIFT = (double*)malloc(MEMD);
  SEARCHEFF_zHOST_LEGACY[0].EFF      = (double*)malloc(MEMD);

  while( (fscanf(fp, "%s", c_get )) != EOF) {

    if( strcmp(c_get,"FIELD:") == 0 || strcmp(c_get,"FIELDLIST:") == 0 ) {
      readchar(fp,FIELDLIST);
      sprintf(SEARCHEFF_zHOST_LEGACY[NMAP].FIELDLIST,"%s", FIELDLIST);

      if ( NMAP > 0 ) {
	SEARCHEFF_zHOST_LEGACY[NMAP].REDSHIFT = (double*)malloc(MEMD);
	SEARCHEFF_zHOST_LEGACY[NMAP].EFF      = (double*)malloc(MEMD);
      }
      NMAP++ ;
    }

    if ( strcmp(c_get,"CUTWIN_SNRMAX:") == 0 ) 
      { readdouble(fp, 2, INPUTS_SEARCHEFF.CUTWIN_SNRMAX_zHOST ); }
    

    if(strcmp(c_get,"HOSTEFF:") == 0 ) {

      // logic is tricky because FIELD key is optional
      imap = NMAP ;
      ptrField = SEARCHEFF_zHOST_LEGACY[imap].FIELDLIST ; 

      // if FIELDLIST not given, then associate ALL fields with map
      if ( NMAP==0 && strcmp(ptrField,"NONE")==0 ) 
	{ sprintf(ptrField,"ALL"); NMAP++ ; }

      imap = NMAP-1 ; 
      if ( NMAP < MXMAP_SEARCHEFF_zHOST ) {
	readdouble(fp, 2, VAL);
	NROW = SEARCHEFF_zHOST_LEGACY[imap].NROW ;
	SEARCHEFF_zHOST_LEGACY[imap].REDSHIFT[NROW] = VAL[0] ;
	SEARCHEFF_zHOST_LEGACY[imap].EFF[NROW]      = VAL[1] ;
	SEARCHEFF_zHOST_LEGACY[imap].NROW++ ;

	if ( VAL[1] > 1.0E-9 ) { NONZERO_SEARCHEFF_zHOST++ ; }
      }
    }    
  }  // end while

  INPUTS_SEARCHEFF.NMAP_zHOST = NMAP;

  if ( NMAP >= MXMAP_SEARCHEFF_zHOST ) {
    sprintf(c1err,"NMAP=%d exceeds bound.", NMAP );
    sprintf(c2err,"Check bound MXMAP_SEARCHEFF_zHOST=%d", 
	    MXMAP_SEARCHEFF_zHOST );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  printf("\t Stored %d redshift bins for SEARCHEFF_zHOST map. \n",
	 NROW );  fflush(stdout);

  NROW = SEARCHEFF_zHOST_LEGACY[0].NROW ;
  if ( NROW > 0 ) { 
    return ; 
  }
  else {
    sprintf(c1err,"No redshift bins found for SEARCHEFF_zHOST map.");
    sprintf(c2err,"Check SEARCHEFF_zHOST_FILE key");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;    
  }

  return ;

} // end of read_zHOST_FILE_LEGACY


// ********************************
void read_VARNAMES_zHOST(FILE *fp) {

  // Mar 20 2019
  //
  // * Read VARNAMES keys from *fp
  // * set INPUTS_SEARCHEFF.IVERSION_zHOST
  // * fill SEARCHEFF_zHOST[0].VARNAMES_HOSTLIB[ivar]
  //

  int MXVAR = MXVAR_SEARCHEFF_zHOST;
  int ivar, NVAR, NKEY, UNIQUE[MXVAR_SEARCHEFF_zHOST];
  char *VARNAMES[MXVAR_SEARCHEFF_zHOST];
  char fnam[] = "read_VARNAMES_zHOST" ;

  // -------------- BEGIN ---------------

  SEARCHEFF_zHOST[0].NVAR = 0;
  if ( INPUTS_SEARCHEFF.IVERSION_zHOST > 0 ) { return ; }

  for(ivar=0; ivar < MXVAR; ivar++ )
    { VARNAMES[ivar] = SEARCHEFF_zHOST[0].VARNAMES_HOSTLIB[ivar]; }

  read_VARNAMES_KEYS(fp, MXVAR, 1, fnam,&NVAR, &NKEY, UNIQUE, VARNAMES) ;
  SEARCHEFF_zHOST[0].NVAR = NVAR;

  if ( NVAR == 0 )
    { INPUTS_SEARCHEFF.IVERSION_zHOST = IVERSION_zHOST_LEGACY; }
  else
    { INPUTS_SEARCHEFF.IVERSION_zHOST = IVERSION_zHOST_MULTID; }


  return ;

} // end  read_VARNAMES_zHOST


// ***************************************
int gen_SEARCHEFF ( int ID                  // (I) identifier 
		    ,double *EFF_SPEC       // (O) spec eff
		    ,double *EFF_zHOST      // (O) EFF(zHOST)
		    ,MJD_DETECT_DEF *MJD_DETECT 
		    ) {

  /******
     Created Jan 7,2008 by R.Kessler

     Returns bitmask as follows:
      bit 1(1) => passes software search based on EFF vs. SNR or MAG
      bit 2(2) => pass human search using hard-wired function
               in human_search_eff.c
      bit 3(4) => failed spec trigger but has spec zHost.

     LFIND1,LFIND2 refer to software,human efficiencies, respectively.

   Beware that return arg *EFF_SPEC is really just the 
   spectroscopic efficiency because there is currently
   no way to calculate the pipeline efficiency.


  May 8 2014: call gen_SEARCHEFF_zHOST() and set 4-bit of MASK

  Aug 22, 2014: pass ID (integer identifier)

  Jan 28 2015: check spec eff only if found by pipeline
               Clarify variable names of LFIND* logicals.

  Mar 2018: 
   + use APPLYMASK_SEARCHEFF_xxx parmaeters to set return MASK
   + add argument *EFF_zHOST

  *****/

  int  LFIND1_PIPELINE, LFIND2_SPEC, LFIND3_zHOST, MASK ;
  char fnam[]  = "gen_SEARCHEFF" ;

  // ----------------- BEGIN -------------


  // init function args
  *EFF_SPEC = *EFF_zHOST = 0.0;  
  MJD_DETECT->TRIGGER = 1.0E6 ;
  MJD_DETECT->FIRST   = 1.0E6 ;
  MJD_DETECT->LAST    = 1.0E6 ;
  MASK = 0 ;  // init return arg

  // 5/26/2009: minimum trigger is at least 2 measurements
  if ( SEARCHEFF_DATA.NOBS < INPUTS_SEARCHEFF.MINOBS ) { return MASK ; }

  LFIND1_PIPELINE = LFIND2_SPEC = LFIND3_zHOST = 0 ;

  // --------------------------
  // if nothing is defined for software efficieincy,
  // just set it to FOUND and skip to SPEC eff.

  if ( INPUTS_SEARCHEFF.NMAP_DETECT == 0 ) {
    LFIND1_PIPELINE = 1;      // set flag that search software finds SN
  }
  else {
    // check check software/pipeline finds this SN
    LFIND1_PIPELINE = gen_SEARCHEFF_PIPELINE(ID,MJD_DETECT);
  }


  // ------------------------------------
  // if detected by pipeline, check if spec-confirmed
  if ( LFIND1_PIPELINE )  { 
    LFIND2_SPEC = gen_SEARCHEFF_SPEC(ID, EFF_SPEC) ;  // return EFF_SPEC
  }


  // if not spec-confirmed, check for spec zHOST
  if ( LFIND1_PIPELINE && LFIND2_SPEC == 0 )  { 
    LFIND3_zHOST = gen_SEARCHEFF_zHOST(ID,EFF_zHOST) ; // return EFF_zHOST
  }

  // --- set return bit-MASK ---- 

  if ( LFIND1_PIPELINE  ) 
    { MASK += APPLYMASK_SEARCHEFF_PIPELINE ; } // software trigger 

  if ( LFIND2_SPEC      ) 
    { MASK += APPLYMASK_SEARCHEFF_SPEC  ; } // got spec-confirmation 

  if ( LFIND3_zHOST     ) 
    { MASK += APPLYMASK_SEARCHEFF_zHOST ; } // got host-gal redshift (May 2014)

  return(MASK) ;

} // end of gen_SEARCHEFF



// ***************************************
int gen_SEARCHEFF_PIPELINE(int ID, MJD_DETECT_DEF *MJD_DETECT) {

  // Created Jul 3, 2011 by R.Kessler
  //
  // returns 1 if this SN is found by pipeline; returns 0 otherwise.
  // (I) integer identifier for debug dump and abort
  // (O) MJD_DETECT->TRIGGER is the MJD where the pipeline trigger 
  //      is first satisfied.
  // (O) MJD_DETECT-FIRST[LAST] is MJD of first and last detection
  //
  // Global outputs per observation:
  //    SEARCHEFF_DATA.detectFlag[obs] 
  //         += 1 for detection
  //         += 2 on epoch satisfying trigger
  //
  //    SEARCHEFF_DATA.PHOTPROB[obs] 
  //
  // Jan 2014: use SEARCHEFF_DATA struct instead of GENLC
  //
  // Oct 7 2017: fix buggy logic for MJD_TRIGGER and possibly the
  //             logic for determining trigger logic
  //             Bug impacts triggers formed within 1 night.
  //
  // Mar 14 2018: 
  //   + check to compute PHOTPROB from map.
  //   + on epoch satisfying trigger, set 2-bit for detectFlag
  //   + Refactor and fix logic to ignore MAG_UNDEFINED epochs
  //   + Refactor and fix logic so that MJD_TRIGGER is correct
  //
  // Feb 17 2020
  //  + refactor to implement CUTVAL cut on PHOTPROB and use
  //    it as part of detection.
  //
  // Oct 18 2021: load MJD_DETECT-FIRST[LAST]
  //

  int NMJD_DETECT, NDETECT, imask, NOBS, MARK, DETECT_MARK, IMAP ;
  int IFILTOBS, obs, OVP, obsLast, istore, LFIND, FIRST=0;
  int IFILTOBS_MASK, IFILTDEF_MASK, NEXT_DETECT, DETECT_FLAG ;
  int FOUND_TRIGGER=0,  FOUND_DETECT_FIRST=0, LCUT_PHOTPROB ;
  int OBSMARKER_DETECT[MXOBS_TRIGGER];
  double  RAN, EFF, MJD, MJD_LAST, MJD_DIF, TDIF_NEXT, SNR,MAG;
  double  PHOTPROB, CUTVAL ;
  char CFILT[4];
  int LDMP  = (ID == -39 ); 
  char fnam[] = "gen_SEARCHEFF_PIPELINE";

  // ------------- BEGIN -------------

  if ( INPUTS_SEARCHEFF.FUNEFF_DEBUG ) {
    RAN   = SEARCHEFF_RANDOMS.FLAT_PIPELINE[0] ;
    return gen_SEARCHEFF_DEBUG("PIPELINE", RAN, &EFF) ;
  }


  NMJD_DETECT = 0;
  MJD_LAST    = SEARCHEFF_DATA.MJD[0] ;
  NMJD_DETECT = 0 ;
  TDIF_NEXT   = INPUTS_SEARCHEFF.TIME_SINGLE_DETECT ;
  NOBS        = SEARCHEFF_DATA.NOBS ;

  if (LDMP ) {
    printf("\n xxx \n");
    printf(" xxx ------ DUMP CID=%d  NOBS=%d  "
	   "TDIF_NEXT=%.3f -------------- \n", 	ID, NOBS, TDIF_NEXT );    
  }

  // fill OBSMARKER_DETECT to mark end of each detection period.
  obsLast = -9;
  for(obs = 0 ; obs < NOBS; obs++ ) {
    OBSMARKER_DETECT[obs] = 0 ;
    MAG  = SEARCHEFF_DATA.MAG[obs] ;
    MJD  = SEARCHEFF_DATA.MJD[obs] ;

    if ( MAG == MAG_UNDEFINED ) { continue ; }
    if ( !FIRST ) {  FIRST=1; obsLast=obs; MJD_LAST= MJD; continue;  }       
    MJD_DIF     = fabs(MJD-MJD_LAST);
    NEXT_DETECT = ( MJD_DIF > TDIF_NEXT );

    if ( NEXT_DETECT && obsLast >= 0 ) 
      { OBSMARKER_DETECT[obsLast] = 1 ; MJD_LAST=MJD; }

    obsLast = obs;
  }  // end obs loop
  OBSMARKER_DETECT[obsLast] = 1 ;  // always set marker for last obs 

  // ---------
  
  NDETECT = IFILTOBS_MASK = LFIND = DETECT_MARK = 0 ;
  OBS_PHOTPROB.NSTORE = 0 ;

  // loop over each epoch and determine if there is a detection,
  // and also if there is a PHOTPROB measurement.
  for(obs = 0 ; obs < SEARCHEFF_DATA.NOBS; obs++ ) {
    OBS_PHOTPROB.IMAP_LIST[obs]    = -9;
    SEARCHEFF_DATA.detectFlag[obs] = 0; 

    MAG      = SEARCHEFF_DATA.MAG[obs] ;
    if ( MAG == MAG_UNDEFINED ) { continue ; }

    IFILTOBS = SEARCHEFF_DATA.IFILTOBS[obs] ;
    RAN      = SEARCHEFF_RANDOMS.FLAT_PIPELINE[obs] ;
    EFF      = GETEFF_PIPELINE_DETECT(obs); // compute effic
    DETECT_FLAG =  ( RAN < EFF ) ;
    if ( DETECT_FLAG ) 
      { SEARCHEFF_DATA.detectFlag[obs] += DETECT_MASK_SNR; }

    setObs_for_PHOTPROB(DETECT_FLAG,obs); // set obs for PHOTPROB

    // all epochs have "good PHOTPROB" by default
    SEARCHEFF_DATA.detectFlag[obs] += DETECT_MASK_PHOTPROB ; 

    if ( LDMP ) {
      printf(" xxx detectFlag[%2d] = %d (IFILTOBS=%d) EFF=%.3f RAN=%.3f\n",
	     obs, SEARCHEFF_DATA.detectFlag[obs], IFILTOBS, EFF, RAN);
      fflush(stdout);
    }

  }  // end obs loop

  // - - - - - - - - - - - - - - -
  // Compute photprob AFTER finding all PHOTPROB epochs so that
  // PHOTPROB covariance can be included.
  if ( OBS_PHOTPROB.NSTORE > 0 ) {
    setRan_for_PHOTPROB();
    for(istore=0; istore < OBS_PHOTPROB.NSTORE ; istore++ ) {
      obs      = OBS_PHOTPROB.OBS_LIST[istore];
      IMAP     = OBS_PHOTPROB.IMAP_LIST[istore] ;
      PHOTPROB = get_PIPELINE_PHOTPROB(istore); 
      // xxx if ( (ID%3) == 0 ) { PHOTPROB=0.7; } else { PHOTPROB=0.3; }

      SEARCHEFF_DATA.PHOTPROB[obs] = PHOTPROB ;

      // if PHOTPROB fails cut, then turn off DETECTFLAG bit
      CUTVAL = SEARCHEFF_PHOTPROB[IMAP].CUTVAL ;
      LCUT_PHOTPROB = ( PHOTPROB > CUTVAL ) ;
      if ( !LCUT_PHOTPROB ) 
	{ SEARCHEFF_DATA.detectFlag[obs] -= DETECT_MASK_PHOTPROB ; }
    } 
    dumpLine_PIPELINE_PHOTPROB(); // optional dump
  } // end PHOTPROB if-block


  // - - - - - - - - - - - - - - -
  // loop again over observations to determine trigger.
  // Check both detection and PHOTPROB.

  for(obs = 0 ; obs < SEARCHEFF_DATA.NOBS; obs++ ) {

    SNR      = SEARCHEFF_DATA.SNR[obs] ;
    MAG      = SEARCHEFF_DATA.MAG[obs] ;

    if ( MAG == MAG_UNDEFINED ) { continue ; }

    IFILTOBS    = SEARCHEFF_DATA.IFILTOBS[obs] ;
    MJD         = SEARCHEFF_DATA.MJD[obs] ;
    DETECT_FLAG = SEARCHEFF_DATA.detectFlag[obs];
    MARK        = OBSMARKER_DETECT[obs] ;

    // set filter-detect mask if both detection and PHOTPROB are satisfied.
    if ( (DETECT_FLAG & 1) > 0  && (DETECT_FLAG & 4)>0 ) { 
      IFILTOBS_MASK |= ( 1 << IFILTOBS )  ; 

      if ( MJD_DETECT->FIRST > 0.99E6 ) { MJD_DETECT->FIRST = MJD; }
      MJD_DETECT->LAST = MJD ;
    }

   
    // Check for trigger if trigger not yet formed, and at least
    // one obs-group marker has passed.  Check IFILTOBS_MASK for detection(s).

    FOUND_TRIGGER      = (MJD_DETECT->TRIGGER < 0.99E6 ) ;
 
    if ( !FOUND_TRIGGER  ) {

      NDETECT = 0 ;  // reset number of detections 
      for ( imask=0; imask < SEARCHEFF_LOGIC.NMASK; imask++ ) {
	IFILTDEF_MASK = SEARCHEFF_LOGIC.IFILTDEF_MASK[imask];
	OVP = IFILTDEF_MASK & IFILTOBS_MASK ;
	if ( OVP == IFILTDEF_MASK ) { NDETECT++ ; }
      }

      if ( NDETECT>0 && !DETECT_MARK ) { NMJD_DETECT++;  DETECT_MARK=1; }

      if ( NMJD_DETECT >= SEARCHEFF_LOGIC.NMJD ) {
	LFIND = 1 ;   MJD_DETECT->TRIGGER = MJD; 
	SEARCHEFF_DATA.detectFlag[obs] += DETECT_MASK_MJD_TRIGGER ; ;
	
	if ( LDMP ) {
	  printf(" xxx \t NMJD_DETECT=%d at MJD_TRIGGER=%.4f \n",
		 NMJD_DETECT, MJD );	   
	  printf(" xxx \t detectFlag=%d  MAG(%s)=%.2f  SNR=%.2f\n",
		 SEARCHEFF_DATA.detectFlag[obs], CFILT, MAG,SNR );
	}
      }

    } // end  !FOUND_TRIGGER 

    // reset MASK after each detection period
    if ( MARK ) { IFILTOBS_MASK = DETECT_MARK = 0 ; } 

  } // end obs loop

    //  if ( LDMP ) {  debugexit(fnam); } // check DUMP ID = 94

  return LFIND ;

} // end of gen_SEARCHEFF_PIPELINE



// ***************************************
void dumpLine_PIPELINE_PHOTPROB(void) {

  // Created Apr 2018
  // One-line dump of PHOTPROB(s) as follows:
  //
  //  "777777  [SNID]  [PHOTPROB-1] [PHOTPROB-2] ... [PHOTPROB-NDUMP]"
  //
  // The string 777777 allows using grep to extract from stdout.
  // The first "NDUMP" PHOTPROB>0 values are dumped;
  // PHOTPROB=0 observations are skipped.
  //
  int  NDUMP  = INPUTS_SEARCHEFF.NPHOTPROB_DUMP ;
  int  i, obs, NFOUND=0;
  double PHOTPROB, PHOTPROB_LIST[MXOBS_PHOTPROB];
  //  char fnam[] = "dumpLine_PIPELINE_PHOTPROB";
  // --------------- BEGIN --------------

  if ( NDUMP == 0 ) { return ; }

  // loop over all PHOTPROB and stop when we get "NDUMP" of them.

  for(obs=0; obs < SEARCHEFF_DATA.NOBS; obs++ ) {
    PHOTPROB = SEARCHEFF_DATA.PHOTPROB[obs];
    if ( PHOTPROB > 0.0 && NFOUND < NDUMP ) {
      PHOTPROB_LIST[NFOUND] = PHOTPROB;
      NFOUND++ ;
    }
  }
      
  if ( NFOUND < NDUMP ) { return ; }

  char LINE[200], cval[20];
  sprintf(LINE," 777777 %d  ", SEARCHEFF_DATA.CID );
  for(i=0; i < NDUMP; i++ ) {
    sprintf(cval, "%6.3f ", PHOTPROB_LIST[i] ); 
    strcat(LINE,cval);
  }

  // print line to stdout
  printf("%s\n", LINE);

  return;

}  // dumpLine_PIPELINE_PHOTPROB


// ***************************************
double GETEFF_PIPELINE_DETECT(int obs) {

  // Return pipeline/detection search efficiency for this obs.
  // Note that obs specifies both MJD and filter.
  //
  // Jan 3 2018: never use saturated epoch for trigger (see NPE_SAT)
  // Mar 7 2018: refactor to use new SEARCHEFF_DETECT struct.
  // Sep   2018: fix logic when NMAP_FOUND=0. Effects SDSS u & z bands
  //             which don't have efficiency maps, whild gri do.
  // Feb 15 2022: add more info for isnan abort.

  double
    MAG, SNR, MJD, EFF, FIX_EFF
    ,EFF_atmax, EFF_atmin, VAL_atmax, VAL_atmin, VAL
    ,ZERO = 0.0
    ,ONE  = 1.0
    ;

  int NMAP = INPUTS_SEARCHEFF.NMAP_DETECT ;

  int 
    CID, ifilt_obs, NPE_SAT
    ,OPT_INTERP  = 1   // 1=linear;  2=quadratic
    ,NBIN_EFF, imap, IMAP, NMAP_FOUND=0
    ;

  char cfilt[4];
  char fnam[] ="GETEFF_PIPELINE_DETECT" ;

  // ---------- BEGIN ---------

  // check debugging option with fixed effic.
  FIX_EFF = INPUTS_SEARCHEFF.FIX_EFF_PIPELINE ;
  if ( FIX_EFF > 0.0 ) { return FIX_EFF ; }

  EFF       = 0.0 ;

  // find map corresponding to filter
  ifilt_obs = SEARCHEFF_DATA.IFILTOBS[obs] ;  IMAP=-9;
  sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );
  for(imap=0; imap < NMAP; imap++ ) {
    if ( strstr(SEARCHEFF_DETECT[imap].FILTERLIST,cfilt) != NULL ) 
      {  IMAP = imap;   NMAP_FOUND++; }
  }

  // if no maps are found for this filter, there are two possibilities:
  // 1) there are no maps at all --> return EFF=1
  // 2) there are maps for other bands -> return EFF=0
  if ( NMAP_FOUND == 0 ) { 
    if ( NMAP == 0 ) { return(ONE); }  else { return(ZERO); }
  }


  if ( NMAP_FOUND > 1 ) {
    sprintf(c1err,
	    "Found %d PIPELINE/DETECT maps for ifilt_obs=%d(%s)",
	    NMAP_FOUND, ifilt_obs, cfilt);
    sprintf(c2err,"Check EFF maps in %s", 
	    INPUTS_SEARCHEFF.PIPELINE_EFF_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  CID       = SEARCHEFF_DATA.CID ;
  NBIN_EFF  = SEARCHEFF_DETECT[IMAP].NBIN;
  SNR       = SEARCHEFF_DATA.SNR[obs] ;
  MAG       = SEARCHEFF_DATA.MAG[obs] ;
  MJD       = SEARCHEFF_DATA.MJD[obs] ;
  NPE_SAT   = SEARCHEFF_DATA.NPE_SAT[obs] ; // Npe above saturation

  if ( isnan(SNR) ) {
    print_preAbort_banner(fnam);

    double PEAKMJD = SEARCHEFF_DATA.PEAKMJD ;
    double z       = SEARCHEFF_DATA.REDSHIFT ;
    double Trest   = (MJD-PEAKMJD)/(1.0+z);
    double FLUX    = SEARCHEFF_DATA.FLUX[obs];
    double FLUXERR = SEARCHEFF_DATA.FLUXERR[obs];

    printf("\t MAG=%.3f  MJD=%.3f  PEAKMJD=%.3f  Trest=%.3f  z=%.3f\n", 
	   MAG, MJD, PEAKMJD, Trest, z );
    printf("\t NPE_SAT   = %d\n", NPE_SAT);
    printf("\t FLUX(ADU) = %f +_ %f \n",  FLUX, FLUXERR);
    printf("\t MWEBV = %.3f \n",   SEARCHEFF_DATA.MWEBV);
    printf("\t LIBID     = %d \n", SEARCHEFF_DATA.SIMLIB_ID );

    sprintf(c1err,"Invalid SNR=isNan for CID=%d, ifilt_obs=%d(%s)",
	    CID, ifilt_obs, cfilt );
    sprintf(c2err,"See info above in pre-abort dump.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }
  
       
  if ( NBIN_EFF == 0   ) { return(ONE) ; } // nothing defined 

  // skip negative SNR unless ABSSNR key is used in PIPELINE-eff file
  if ( SNR <= 0.0 &&  SEARCHEFF_FLAG != FLAG_EFFABSSNR_DETECT ) 
    { return(ZERO) ; } 

  // Jan 3 2018: never use saturated epoch for trigger 
  if ( NPE_SAT > 0 ) { return(ZERO); } 

  EFF_atmax = SEARCHEFF_DETECT[IMAP].EFF[NBIN_EFF-1] ;
  EFF_atmin = SEARCHEFF_DETECT[IMAP].EFF[0] ;
  VAL_atmax = SEARCHEFF_DETECT[IMAP].VAL[NBIN_EFF-1] ;
  VAL_atmin = SEARCHEFF_DETECT[IMAP].VAL[0] ;

  // - - - - - - - - - - - - - -
  // NOTES FOR FUTURE OPTION TO TRIGGER ON SINGLE-IMAGE DETECTIONS
  // Jun 2 2022 .xyz 
  // for option to simulated single-exposure detections, 
  // strip off NEXPOSE here and SNR -> SNR / sqrt(NEXPOSE). 
  // For NEXPOSE chances to get a detection, 
  //    EFF -> 1-(1-EFF_1)^NEXPOSE
  // E.g., NEXPOSE=10 and EFF_1 = 0.02 -> EFF = 1 - 0.98^10 = 0.183
  // E.g., NEXPOSE=10 and EFF_1 = 0.05 -> EFF = 1 - 0.95^10 = 0.401
  // E.g., NEXPOSE=10 and EFF_1 = 0.50 -> EFF = 1 - 0.50^10 = 0.999
  // Above assumes stat-independent detection prob .... but is this
  // correct, or are the probabilities correlated?
  // - - - - - - - - - - - - - -

  /*   

   bool LDETCT_SINGLE = TRUE;
   if LDETECT_SINGLE{
    SNR *= 1.0/(double) NEXPOSE; 

}


 */


  if ( SEARCHEFF_FLAG == FLAG_EFFSNR_DETECT ) {
    VAL = SNR ;
  }
  else if ( SEARCHEFF_FLAG == FLAG_EFFABSSNR_DETECT ) {
    VAL = fabs(SNR) ;
  }
  else {
    VAL = MAG ;
  }

  if ( VAL > VAL_atmax ) { return EFF_atmax ; }
  if ( VAL < VAL_atmin ) { return EFF_atmin ; }

  // interpolate
  EFF = interp_1DFUN (OPT_INTERP, VAL, NBIN_EFF, 
		      SEARCHEFF_DETECT[IMAP].VAL,
		      SEARCHEFF_DETECT[IMAP].EFF, fnam);

  // printf(" xxxx EFF(VAL=%f) = %f   (inear = %d) \n", VAL, EFF, inear );

  return EFF ;

} // end of  GETEFF_PIPELINE_DETECT


// *************************************
void setObs_for_PHOTPROB(int DETECT_FLAG, int obs) {

  // Created Apr 2018
  // Store observations for which PHOTPROB will be computed.
  //
  // Inputs:
  //    DETECT_FLAG = 1 if there is a detection for this epoch
  //    DETECT_FLAG = 0 if no detection
  //    obs         = observer index for SEARCHEFF_DATA array
  //   
  
  int   NMAP       = INPUTS_SEARCHEFF.NMAP_PHOTPROB ;
  int   IFILTOBS   = SEARCHEFF_DATA.IFILTOBS[obs] ;
  char  *FIELD     = SEARCHEFF_DATA.FIELDNAME ; 

  int  NSTORE = OBS_PHOTPROB.NSTORE;
  int  imap, IMAP, NMATCH;
  bool MATCH_FIELD, MATCH_FILT ;
  char FILT[2], *FIELD_TMP, *FILT_TMP;
  char fnam[]      = "setObs_for_PHOTPROB" ;

  // ------------ BEGIN ------------

  if ( NMAP == 0 ) { return ; }


  sprintf(FILT, "%c", FILTERSTRING[IFILTOBS] );

  // find map for this filter and field.
  NMATCH = 0 ;  IMAP=-9;
  for(imap=0; imap < NMAP; imap++ ) {
    MATCH_FIELD = MATCH_FILT = false ;
    FIELD_TMP  = SEARCHEFF_PHOTPROB[imap].FIELDLIST;
    FILT_TMP   = SEARCHEFF_PHOTPROB[imap].FILTERLIST ;

    MATCH_FIELD = MATCH_SEARCHEFF_FIELD(FIELD_TMP); // Feb 2021

    if ( strcmp(FILT_TMP,"ALL")==0 )       { MATCH_FILT=1; }
    if ( strstr(FILT_TMP,FILT ) != NULL  ) { MATCH_FILT=1; }

    if ( MATCH_FIELD && MATCH_FILT ) { IMAP = imap;  NMATCH++ ; }

  } // end imap loop

  if(NMATCH==0 )  { return; }

  if(NMATCH > 1 ) {
    sprintf(c1err,"%d matches to PHOTPROB map invalid.", NMATCH );
    sprintf(c2err,"FIELD='%s'  FILT='%s' ", FIELD, FILT);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  // check if PHOTPROB map requires a detection
  if ( DETECT_FLAG==0 && SEARCHEFF_PHOTPROB[IMAP].REQUIRE_DETECTION )
    { return; }

  
  // if we get here, store info
  if ( NSTORE < MXOBS_PHOTPROB ) {
    OBS_PHOTPROB.OBS_LIST[NSTORE]  = obs ;
    OBS_PHOTPROB.OBSINV_LIST[obs]  = NSTORE ;
    OBS_PHOTPROB.IMAP_LIST[NSTORE] = IMAP ;
  }
  OBS_PHOTPROB.NSTORE++ ;

  return ;

}  // end setObs_for_PHOTPROB" ;


// ***************************************************
void setRan_for_PHOTPROB(void) {

  // Created April 2018
  // Set uniform randoms, U[0,1], for each PHOTPROB observation
  // and store in OBS_PHOTPROB.RAN_LIST[istore=0,NSTORE-1] 
  //
  // REDUCED_CORR == 0 -->
  //  Use already-determined randoms from SEARCHEFF_RANDOMS.PHOTPROB[obs]
  //
  // REDUCED_CORR != 0 --> 
  //   Compute correlated randoms. Start with Gaussian randoms in
  //   SEARCHEFF_RANDOMS.PHOTPROB and compute correlated Gaussian
  //   randoms (GRAN) using Cholesky decomp. Convert each GRAN
  //   into U[0,1] = Gaussian integral from -9.9 to GRAN.
  //   (-9.9 is an approx for -infinity due to integral grid)
  //
  // Beware of two distinct sets of indices.
  // "istore = 0, NSTORE-1" is a sparse index over observations
  // for which PHOTPROB is computed (e.g., those with detections).
  // "obs" is the observation index in SEARHEFF_DATA struct.
  //

  int  NSTORE     = OBS_PHOTPROB.NSTORE ;
  int  NMAP       = INPUTS_SEARCHEFF.NMAP_PHOTPROB ;
  double GAURAN_MAX = GAUSS_INTEGRAL_STORAGE.XMAX - 0.1 ;
  double COVDIAG = 1.0;  // Cov matrix has sigma=1

  CHOLESKY_DECOMP_DEF DECOMP;
  int  MEMD, istore, imap, imap1, irow, irow1, NCOV, obs ;
  double CORR, CORR1, FLATRAN, GAURAN ;
  double GAURAN_LIST[MXOBS_PHOTPROB], GAURANCORR_LIST[MXOBS_PHOTPROB];
  char   fnam[] = "setRan_for_PHOTPROB" ;


  // ------------ BEGIN -----------

  if ( NMAP == 0 ) { return ; }

  if ( NSTORE >= MXOBS_PHOTPROB ) {
    sprintf(c1err,"NSTORE = %d exceed bound MXOBS_PHOTRPOB=%d", 
	    NSTORE, MXOBS_PHOTPROB );
    sprintf(c2err,"Check MXOBS_PHOTPROB in sntools_trigger.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  if ( INPUTS_SEARCHEFF.NREDUCED_CORR_PHOTPROB == 0 ) { 
    // no photprob correlations --> use already-stored flatran[0,1]
    for(istore=0; istore < NSTORE; istore++ ) {
      obs  = OBS_PHOTPROB.OBS_LIST[istore] ;
      OBS_PHOTPROB.RAN_LIST[istore] = SEARCHEFF_RANDOMS.FLAT_PHOTPROB[obs] ; 
    }
    return ; 
  }



  // if we get here, compute correlated randoms.
  // Start with correlated Gaussians, then transform to 
  // uniform [0,1] distribution.


  // construct 1D cov matrix for unit-width Gaussians
  MEMD  = NSTORE * NSTORE * sizeof(double);
  DECOMP.COVMAT1D = (double*) malloc (MEMD); 
  DECOMP.MATSIZE = NSTORE ;

  // xxx  COVTMP1D = (double*) malloc (MEMD);  // REMOVE LATER
  NCOV  = 0;
  for (irow=0; irow < NSTORE; irow++) {
    for (irow1 = 0; irow1 < NSTORE ; irow1++) {       
      imap  = OBS_PHOTPROB.IMAP_LIST[irow];
      imap1 = OBS_PHOTPROB.IMAP_LIST[irow1];
      CORR  = SEARCHEFF_PHOTPROB[imap].REDUCED_CORR;
      CORR1 = SEARCHEFF_PHOTPROB[imap1].REDUCED_CORR;
      if ( irow == irow1 ) 
	{ DECOMP.COVMAT1D[NCOV] = COVDIAG; }
      else
	{ DECOMP.COVMAT1D[NCOV] = sqrt(CORR*CORR1); }

      // xxx      COVTMP1D[NCOV] = DECOMP.COVMAT1D[NCOV] ; // REMOVE LATER
      NCOV++ ;
    }
  }
 
  init_Cholesky(+1, &DECOMP);

  for(irow=0; irow < NSTORE; irow++ ) {
    obs    = OBS_PHOTPROB.OBS_LIST[irow] ;
    GAURAN_LIST[irow] = SEARCHEFF_RANDOMS.GAUSS_PHOTPROB[obs];
  }
  getRan_GaussCorr(&DECOMP, GAURAN_LIST, // (I)
		   GAURANCORR_LIST);     // (O)

  double x0=0.0 ;
  for(irow=0; irow < NSTORE; irow++ ) {

    GAURAN = GAURANCORR_LIST[irow] ;

    // keep GAURAN within limits of pre-defined Gauss integral grid
    if ( GAURAN > +GAURAN_MAX ) { GAURAN = +GAURAN_MAX; }
    if ( GAURAN < -GAURAN_MAX ) { GAURAN = -GAURAN_MAX; }

    // now convert random correlated Gaussian into flatran[0,1]
    // using Gaussian CDF. Beware that GaussIntegral integrates
    // from zero to positive NSIGMA, not from -infinity.
    if ( GAURAN >= 0.0 ) 
      { FLATRAN = 0.5 + GaussIntegral(x0,GAURAN); } 
    else
      { FLATRAN = 0.5 - GaussIntegral(x0,fabs(GAURAN)); } 

    if ( FLATRAN < 0.0 || FLATRAN > 1.0 ) {
      sprintf(c1err,"Invalid FLATRAN=%f from GUARAN=%f", FLATRAN, GAURAN);
      sprintf(c2err,"NSTORE = %d", NSTORE);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
    GAURAN_LIST[irow]            = GAURAN; // store for debug only
    OBS_PHOTPROB.RAN_LIST[irow]  = FLATRAN ;
  } // end irow loop


  if ( NSTORE < -3 ) {
    printf(" 777777 %6.3f %6.3f %6.3f    %6.3f %6.3f %6.3f \n",
	   OBS_PHOTPROB.RAN_LIST[0], 
	   OBS_PHOTPROB.RAN_LIST[1],
	   OBS_PHOTPROB.RAN_LIST[2],
	   GAURAN_LIST[0], GAURAN_LIST[1], GAURAN_LIST[2] );
  }


  // free memory
  init_Cholesky(-1, &DECOMP);

  return;

} // end setRan_for_PHOTPROB


// ***************************************************
double get_PIPELINE_PHOTPROB(int istore) {

  // Created April 2018
  // Determine random PHOTPROB for this 'istore'.
  // Before calling this function, must call 
  // setObs_for_PHOTPROB & setRan_for_PHOTPROB

  int    obs    = OBS_PHOTPROB.OBS_LIST[istore];  // SEARCHEFF_DATA index
  int    IMAP   = OBS_PHOTPROB.IMAP_LIST[istore]; // select PHOTPROB map
  double RANCDF = OBS_PHOTPROB.RAN_LIST[istore];  // retreive random [0,1]

  int    NVAR_MAP  = SEARCHEFF_PHOTPROB[IMAP].NVAR_MAP;
  int    NFUN      = SEARCHEFF_PHOTPROB[IMAP].NFUN_CDF ;
  int    CID       = SEARCHEFF_DATA.CID ;
  int    IFILTOBS  = SEARCHEFF_DATA.IFILTOBS[obs] ;
  double SNR       = SEARCHEFF_DATA.SNR[obs] ;
  double SBMAG     = SEARCHEFF_DATA.SBMAG[IFILTOBS] ;  
  //  double GALMAG    = SEARCHEFF_DATA.HOSTMAG[IFILTOBS] ;  
  double MJD       = SEARCHEFF_DATA.MJD[obs] ;

  int    istat, ivar, LDMP ;
  double PHOTPROB_CDF[MXVAR_SEARCHEFF_PHOTPROB];
  double VARDATA[MXVAR_SEARCHEFF_PHOTPROB];
  char  *VARNAME, cFILT[4];
  double PHOTPROB = 0.0 ;

  char fnam[] = "get_PIPELINE_PHOTPROB" ;

  // ------------ BEGIN --------------

  
  for(ivar=0; ivar < NVAR_MAP; ivar++ ) 
    { VARDATA[ivar] = LOAD_PHOTPROB_VAR(obs,IMAP,ivar);   }
  
  PHOTPROB_CDF[0] = 0.0 ;
  istat = interp_GRIDMAP(&SEARCHEFF_PHOTPROB[IMAP].GRIDMAP, VARDATA, 
			 &PHOTPROB_CDF[1] );  // <== returned  

  // ---------------------------------
  // Finally, use linear interpolation to pick random PHOTPROB 
  // from PHOTPROB CDF.  Note that NFUN+1 includes the zero bin
  // because the zero bin is not stored in the GRIDMAP.
  
  double *ptrCDF  = PHOTPROB_CDF ;
  double *ptrVAL  = SEARCHEFF_PHOTPROB[IMAP].PHOTPROB_CDFBINS ;
  PHOTPROB = interp_1DFUN(1, RANCDF, NFUN+1, ptrCDF, ptrVAL, fnam);


  LDMP = 0 ; // ( fabs(SEARCHEFF_DATA.SNR[obs]-8.0) < 0.1 );
  if ( LDMP ) {
    sprintf(cFILT, "%c", FILTERSTRING[IFILTOBS] );
    printf("\n xxx ------------------------------------------ \n");
    printf(" xxx obs=%d  CID=%d   IFILTOBS=%d(%s)\n", 
	   obs, CID, IFILTOBS, cFILT );
    printf(" xxx MJD=%.3f  SNR=%.2f  SBMAG=%.2f\n", MJD, SNR, SBMAG);
    printf(" xxx IMAP = %d   for  IDMAP=%d \n",
	   IMAP, SEARCHEFF_PHOTPROB[IMAP].GRIDMAP.ID );
    for(ivar=0; ivar < NVAR_MAP; ivar++ ) {
      VARNAME = SEARCHEFF_PHOTPROB[IMAP].VARNAMES[ivar]; 
      printf(" xxx ivar=%d : %s = %f \n", ivar, VARNAME, VARDATA[ivar] );
    }

    printf(" xxx istat(interp) = %d  RANCDF=%.3f --> PHOTPROB=%.3f\n", 
	   istat, RANCDF, PHOTPROB );

    /*
    printf(" xxx BINSIZE=%.3f   ibin=%d   CDF_bin=%.3f,%.3f\n", 
    BINSIZE, ibin, CDF_bin[0], CDF_bin[1] ); */

    printf(" xxx PHOTPROB_BINS: ");
    for(ivar=0; ivar <= NFUN ; ivar++ )  
      { printf(" %.3f",  SEARCHEFF_PHOTPROB[IMAP].PHOTPROB_CDFBINS[ivar]); }
    printf("\n");

    printf(" xxx PHOTPROB_CDF : ");
    for(ivar=0; ivar <= NFUN ; ivar++ )  
      {  printf(" %.3f", PHOTPROB_CDF[ivar] ); }
    printf("\n");
    debugexit(fnam);
  }

  return(PHOTPROB);

}  // end get_PIPELINE_PHOTPROB


// ***************************************************
int gen_SEARCHEFF_SPEC(int ID, double *EFF_SPEC) {

  // Return 1 if this SN is spec-confirmed; return 0 otherwise
  // (O) *EFF_SPEC is the spec-efficiency
  //
  // Oct 1, 2012: check  FIELD
  // May 10 2018: EFF *= USER_SPECEFF_SCALE
  // Jun 02 2018: check opton to force EFF_SPEC=0
  // Jun 08 2018: check BOOLEAN OR/AND logic for EFF.

  int  NMAP        = INPUTS_SEARCHEFF.NMAP_SPEC ;
  int  BOOLEAN_OR  = SEARCHEFF_SPEC_INFO.BOOLEAN_OR  ;
  int  BOOLEAN_AND = SEARCHEFF_SPEC_INFO.BOOLEAN_AND ;

  int  imap, ivar, NVAR, LFIND, istat ;
  bool MATCH ;
  double PnoSpec_OR, Pspec_AND, EFF, RAN, VARDATA[MXVAR_SEARCHEFF_SPEC];
  char *fld_gen, *fld_map ;
  char fnam[] = "gen_SEARCHEFF_SPEC" ;

  // ----------- BEGIN --------

  LFIND = 1 ;
  EFF   = 1.0 ;

  if ( INPUTS_SEARCHEFF.FUNEFF_DEBUG ) {
    RAN   = SEARCHEFF_RANDOMS.FLAT_SPEC[0] ;
    return gen_SEARCHEFF_DEBUG("SPEC", RAN, EFF_SPEC) ;
  }


  // if there is no spec-eff file, return 100% eff and True
  if ( NMAP == 0 )  { 
    if ( INPUTS_SEARCHEFF.IFLAG_SPEC_EFFZERO ) 
      { *EFF_SPEC = 0.0 ; LFIND=0 ; }
    else
      { *EFF_SPEC = EFF ; LFIND=1 ; }

    return LFIND ;  
  }

  // do multi-dimensional interpolation if SPEC-EFF map
  // is specified. Take the logical-OR of each map,
  // which means that EFF = 1 - product(1-Eff_imap)

  PnoSpec_OR = 1.0 ;
  Pspec_AND  = 1.0 ;
  fld_gen = SEARCHEFF_DATA.FIELDNAME ;

  for ( imap=0; imap < NMAP; imap++ ) {

    // check if current field goes with this map
    fld_map = SEARCHEFF_SPEC[imap].FIELDLIST ;

    MATCH = MATCH_SEARCHEFF_FIELD(fld_map); // Feb 2021

    // determine list of variables
    NVAR = SEARCHEFF_SPEC[imap].GRIDMAP.NDIM ;
    for ( ivar=0; ivar < NVAR; ivar++ )
      { VARDATA[ivar] = LOAD_SPECEFF_VAR(imap,ivar); }

    istat = interp_GRIDMAP(&SEARCHEFF_SPEC[imap].GRIDMAP, VARDATA, &EFF );
    PnoSpec_OR *= (1.0 - EFF);
    Pspec_AND  *= EFF ;

  } // imap for SPECEFF

  //  printf(" xxx BOOLEAN(OR,AND) = %d %d \n",  BOOLEAN_OR, BOOLEAN_AND);

  if ( NMAP > 0  ) { 
    if ( BOOLEAN_OR  ) { EFF = 1.0 - PnoSpec_OR ; }
    if ( BOOLEAN_AND ) { EFF = Pspec_AND        ; }
  }

  if ( EFF < 0.0 && EFF > -0.001 ) { EFF = 0.0 ; } // allow round-off

  // sanity check
  if ( EFF < 0.0 || EFF > 1.0001 ) {
    sprintf(c1err,"Invalid Spec effic= %f for ID=%d", EFF, ID);
    sprintf(c2err,"zHel=%6.3f  PEAKMJD=%9.2f ", 
	    SEARCHEFF_DATA.REDSHIFT, SEARCHEFF_DATA.PEAKMJD );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  // apply user scale (May 10 2018)
  EFF *= INPUTS_SEARCHEFF.USER_SPECEFF_SCALE ;

  // load output argument
  *EFF_SPEC = EFF ;  

  // compare spec-eff to random number
  RAN = SEARCHEFF_RANDOMS.FLAT_SPEC[0] ;
  if ( RAN > EFF ) { LFIND = 0 ; }

  return LFIND ;

} // end of gen_SEARCHEFF_SPEC


// ************************************************
int gen_SEARCHEFF_zHOST(int ID, double *EFF_zHOST) {

  // Created May 8 2014:
  // Check efficiency of finding host galaxy redshift;
  // return 0, 1 if not found, found.
  //
  // Jan 22 2017: check for field-dependent HOSTEFF map.
  // Mar 11 2018: refactor
  // Jun 19 2018: check CUTWIN_SNRMAX

  int NMAP = INPUTS_SEARCHEFF.NMAP_zHOST ;
  double  *CUTWIN_SNRMAX = INPUTS_SEARCHEFF.CUTWIN_SNRMAX_zHOST ;
  int     LFIND ;
  double  SNRMAX, RAN, EFF ;
  //  char fnam[] = "gen_SEARCHEFF_zHOST" ;

  // ----------- BEGIN -----------

  // check debug option
  if ( INPUTS_SEARCHEFF.FUNEFF_DEBUG ) {
    RAN = SEARCHEFF_RANDOMS.FLAT_SPEC[50] ;
    return gen_SEARCHEFF_DEBUG("zHOST", RAN, EFF_zHOST) ;
  }

  LFIND = 1 ;  *EFF_zHOST = 0.0 ;

  // if there is no spec-eff file, return 100% eff and True
  if ( NMAP == 0 )  { 
    if ( INPUTS_SEARCHEFF.IFLAG_zHOST_EFFZERO ) 
      { *EFF_zHOST = 0.0;  LFIND=0; }
    else
      { *EFF_zHOST = 1.0 ; LFIND=1; }

    return LFIND ;  
  }

  // check optional SNR cut 
  SNRMAX = SEARCHEFF_DATA.SNRMAX ;
  if ( SNRMAX < CUTWIN_SNRMAX[0] ) { return(0) ; }
  if ( SNRMAX > CUTWIN_SNRMAX[1] ) { return(0) ; }

  if(  INPUTS_SEARCHEFF.IVERSION_zHOST == IVERSION_zHOST_LEGACY ) 
    { EFF = interp_SEARCHEFF_zHOST_LEGACY(); }
  else
    { EFF = interp_SEARCHEFF_zHOST(); }


  *EFF_zHOST = EFF ;  // load function arg.

  // borrow one of the already-allocated SPEC randoms.
  RAN = SEARCHEFF_RANDOMS.FLAT_SPEC[50] ;
  if ( RAN > EFF ) { LFIND = 0 ; }

  /*
  printf(" xxxx %s : FIELD=%s z = %.3f  EFF=%.3f   LFIND=%d \n",
	 fnam, field_data, z, EFF, LFIND); fflush(stdout); // xxxx
  */

  return LFIND ;

}   // end of gen_SEARCHEFF_zHOST


// *******************************************
double interp_SEARCHEFF_zHOST(void) {

  // Mar 12 2019
  // Interpolate multi-D map to get EFF(HOSTLIB properties)
  //
  // July 2020: check PEAKMJD too
  // Oct 15 2020: clarify error message on interp failure.
  // May 27 2021: refactor to take logical OR of multiple maps

  int NMAP = INPUTS_SEARCHEFF.NMAP_zHOST ;
  int IMAP=0, istat, imap, NVAR, ivar, ivar_HOSTLIB, IGAL, NMATCH=0;
  double VARDATA[MXVAR_SEARCHEFF_zHOST];
  double EFF = 0.0, Pnoz=1.0, EFF_TMP, PEAKMJD, *PEAKMJD_RANGE ;
  char *field_map, *field_data, *varName ;  
  bool MATCH_FIELD, MATCH_PEAKMJD ;

  int LDMP = 0;
  char fnam[] = "interp_SEARCHEFF_zHOST" ;

  // ---------------- BEGIN -------------

  // determine which map based on FIELD
  for(imap=0; imap < NMAP;  imap++ ) {
    MATCH_FIELD = MATCH_PEAKMJD = false ;

    field_map  = SEARCHEFF_zHOST[imap].FIELDLIST ;
    field_data = SEARCHEFF_DATA.FIELDNAME ;
    MATCH_FIELD = MATCH_SEARCHEFF_FIELD(field_map);

    PEAKMJD_RANGE = SEARCHEFF_zHOST[imap].PEAKMJD_RANGE ;
    PEAKMJD       = SEARCHEFF_DATA.PEAKMJD ;
    if ( PEAKMJD >= PEAKMJD_RANGE[0] && PEAKMJD <= PEAKMJD_RANGE[1] ) 
      { MATCH_PEAKMJD = true; }

    if ( MATCH_FIELD && MATCH_PEAKMJD ) {
      // load VARDATA from HOSTLIB
      NMATCH++ ;
      NVAR = SEARCHEFF_zHOST[imap].GRIDMAP.NDIM ;
      IGAL = SNHOSTGAL.IGAL ;
      for(ivar=0; ivar < NVAR; ivar++ ) {
	ivar_HOSTLIB  = SEARCHEFF_zHOST[imap].IVAR_HOSTLIB[ivar] ;
	VARDATA[ivar] = HOSTLIB.VALUE_ZSORTED[ivar_HOSTLIB][IGAL] ;
      }
      
      istat = interp_GRIDMAP(&SEARCHEFF_zHOST[imap].GRIDMAP, VARDATA, 
			     &EFF_TMP );        // <== returned  

      Pnoz *= ( 1.0 - EFF_TMP ); // prob if NOT getting zHOST
    }

  } // end imap loop

  // E.g., EFF = 0, 0.7; Pnoz=(1-0)*(1-0.7) = 0.3; EFF = 1-0.3 = 0.7
  EFF = 1.0 - Pnoz ;

  if ( NMATCH == 0 ) {
    sprintf(c1err, "Invalid NMATCH=%d for", NMATCH );
    sprintf(c2err, "field = '%s'  PEAKMJD=%.3f", field_data, PEAKMJD );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  if ( LDMP || istat != SUCCESS ) {
    long long GALID = SNHOSTGAL.GALID ;
    double VALMIN, VALMAX, VAL ;
    char tmpStr[20];
    printf(" \n");
    printf(" xxx --------------------------------------------------- \n");
    printf(" xxx DUMP %s : \n", fnam);
    printf(" xxx   IMAP=%d  FIELD=%s  NVAR=%d  IGAL=%d  GALID=%lld\n", 
	   IMAP, field_data, NVAR, IGAL, GALID ) ;
    for(ivar=0; ivar < NVAR; ivar++ ) {
      varName      = SEARCHEFF_zHOST[IMAP].VARNAMES_HOSTLIB[ivar];
      ivar_HOSTLIB = SEARCHEFF_zHOST[IMAP].IVAR_HOSTLIB[ivar] ;
      VALMIN       = SEARCHEFF_zHOST[IMAP].GRIDMAP.VALMIN[ivar] ;
      VALMAX       = SEARCHEFF_zHOST[IMAP].GRIDMAP.VALMAX[ivar] ;
      VAL          = VARDATA[ivar];
      sprintf(tmpStr,"satisfies");
      if ( VAL < VALMIN || VAL > VALMAX ) { sprintf(tmpStr,"FAILS"); }
      printf(" xxx   %s = %f  (%s map bound: %f to %f) \n", 
	     varName, VAL, tmpStr, VALMIN, VALMAX );
    }
    printf(" xxx   EFF=%f  istat(interp_GRIDMAP)=%d \n", EFF, istat);
    fflush(stdout);
  }


  if ( istat != SUCCESS ) {
    sprintf(c1err, "Could not get zHOST_EFF from map");
    sprintf(c2err, "See info dump above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  return(EFF);

} // end of interp_SEARCHEFF_zHOST


bool MATCH_SEARCHEFF_FIELD(char *field_map) {

  // Created Feb 5 2021
  // Return true of SEARCHEFF_DATA.FIELDLIST_OVP matches field_map.
  // 
  // Ideally all overlap fields are checked, but this causes
  // downstream abort if an event overlaps two FIELDS.
  // Here we only check first field loaded in FIELDLIST_OVP;
  // maybe somebody we'll have an algorithm for choosing
  // among multiple FIELD-dependent maps.
  //
  // Previous logic had checked last overlap field (since each
  // subsequent field had clobbered previous field), and here we
  // check first field ... so new logic (Feb 5 2021) can result
  // in using a different map for field overlaps.

  int  NFIELD_OVP = SEARCHEFF_DATA.NFIELD_OVP;
  int  i;
  char *field_data;
  // ---------- BEGIN ----------

  if ( strcmp(field_map,"ALL")      == 0    ) { return true ; }

  // xxx  for(i=0; i < NFIELD_OVP; i++ ) { // maybe someday ??

  for(i=0; i < 1; i++ ) { // only check first FIELD among overlaps
    field_data = SEARCHEFF_DATA.FIELDLIST_OVP[i];
    if ( strstr(field_map,field_data) != NULL ) { return true ; }
  }

  // if we get here, there is no match -> return false
  return false ;

} // end MATCH_SEARCHEFF_FIELD


// *******************************************
double interp_SEARCHEFF_zHOST_LEGACY(void) {

  // Mar 12 2019
  // Code moved from gen_SEARCHEFF_zHOST.
  // Interpolate map to get z-dependent zSpec efficiency.

  int NMAP = INPUTS_SEARCHEFF.NMAP_zHOST ;
  int OPT_INTERP = 1; // linear interp

  double z, EFF;
  int NMATCH=0, IMAP=-9, imap;   char *field_map, *field_data ;
  bool MATCH_FIELD ;
  char fnam[] = "interp_SEARCHEFF_zHOST_LEGACY" ;

  z          = SEARCHEFF_DATA.REDSHIFT ; 
  field_data = SEARCHEFF_DATA.FIELDNAME ; 

  for(imap=0; imap < NMAP;  imap++ ) {

    field_map   = SEARCHEFF_zHOST_LEGACY[imap].FIELDLIST ;
    MATCH_FIELD = MATCH_SEARCHEFF_FIELD(field_map);
    if ( MATCH_FIELD ) { IMAP = imap ; NMATCH++ ; }
  }

  if ( NMATCH != 1 ) {
    print_preAbort_banner(fnam);
    for(imap=0; imap < NMAP;  imap++ ) {
      field_map = SEARCHEFF_zHOST_LEGACY[imap].FIELDLIST ;
      printf("\t FIELD_MAP[%2d] = %s", imap, field_map );
    }    
    sprintf(c1err, "Invalid NMATCH=%d for", NMATCH );
    sprintf(c2err, "field = '%s'", field_data );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  EFF = interp_1DFUN(OPT_INTERP, z   
		     , SEARCHEFF_zHOST_LEGACY[IMAP].NROW
		     , SEARCHEFF_zHOST_LEGACY[IMAP].REDSHIFT
		     , SEARCHEFF_zHOST_LEGACY[IMAP].EFF,  fnam );

  return(EFF);

} // end interp_SEARCHEFF_zHOST_LEGACY

// *********************************************
void LOAD_PHOTPROB_CDF(int NVAR_CDF, double *VAL ) {

  // Created Mar 2018
  // for input WGT values (*VAL), replace *VAL with CDF.
  // Note hat *VAL is input and output.

  int ivar;
  double SUM=0.0, sum=0.0, TMPVAL[MXVAR_SEARCHEFF_PHOTPROB];
  //  char fnam[] = "LOAD_PHOTPROB_CDF" ;

  // ------------- BEGIN -------------

  for(ivar=0; ivar < NVAR_CDF; ivar++ ) { 
    TMPVAL[ivar] = VAL[ivar]; 
    SUM         += VAL[ivar]; 
    VAL[ivar] = 0.0 ;
  }

  for(ivar=0; ivar < NVAR_CDF; ivar++ ) { 
    sum  += TMPVAL[ivar]; 
    VAL[ivar] = sum/SUM ;
  }

  return;

} // end LOAD_PHOTPROB_CDF


// *************************************
double LOAD_PHOTPROB_VAR(int OBS, int IMAP, int IVAR) {

  // Created Mar 2018
  // Return value of variable corresponding to IVAR 
  // in the PHOTPROB map 'IMAP'
  // OBS is the observation index.
  //
  // Feb 14 2020: add REDSHIFT dependence

  double VALMIN  = SEARCHEFF_PHOTPROB[IMAP].VALMIN[IVAR] ;
  double VALMAX  = SEARCHEFF_PHOTPROB[IMAP].VALMAX[IVAR] ;
  double safety  = 1.0E-9*(VALMAX-VALMIN);

  int    IFILTOBS = SEARCHEFF_DATA.IFILTOBS[OBS] ;
  double SNR      = SEARCHEFF_DATA.SNR[OBS]; 
  double SBMAG    = SEARCHEFF_DATA.SBMAG[IFILTOBS]; 
  double GALMAG   = SEARCHEFF_DATA.HOSTMAG[IFILTOBS]; 
  double REDSHIFT = SEARCHEFF_DATA.REDSHIFT ; 
  int    IVARABS ;
  double VAL=0.0 ;
  char *VARNAME;
  char fnam[] = "LOAD_PHOTPROB_VAR" ;

  // -------------- BEGIN --------------

  if ( SNR < 0.1 ) { SNR=0.1; }

  VARNAME = SEARCHEFF_PHOTPROB[IMAP].VARNAMES[IVAR] ;
  IVARABS = IVARABS_SEARCHEFF_PHOTPROB(VARNAME) ;

  if ( IVARABS == IVARABS_PHOTPROB_SNR ) 
    { VAL = SNR ; }
  else if ( IVARABS == IVARABS_PHOTPROB_LOGSNR ) 
    { VAL = log10(SNR) ; }
  else if ( IVARABS == IVARABS_PHOTPROB_SBMAG ) 
    { VAL = SBMAG ; }
  else if ( IVARABS == IVARABS_PHOTPROB_GALMAG ) 
    { VAL = GALMAG ; }
  else if ( IVARABS == IVARABS_PHOTPROB_REDSHIFT )  // FEB 2020
    { VAL = REDSHIFT ; }
  else {
    sprintf(c1err,"Invalid PHOTPROB-map variable '%s'", VARNAME);
    sprintf(c2err,"OBS=%d  IMAP=%d  IVAR=%d  IVARABS=%d", 
	    OBS, IMAP, IVAR, IVARABS );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  // avoid going off the grid
  if ( VAL < VALMIN ) { VAL = VALMIN + safety ; }
  if ( VAL > VALMAX ) { VAL = VALMAX - safety ; }

  return(VAL);

} // end LOAD_PHOTPROB_VAR

// ************************************************
double LOAD_SPECEFF_VAR(int imap, int ivar) {

  // return value of variable corresponding to this SPECEFF-ivar index
  // June 23 2016: check HOSTMAG and SBMAG 

  int    ifilt_obs, ifilt,  IVARTYPE, NFILTLIST ;
  double mag, MAG, mag0, mag1, color;
  char *varName = SEARCHEFF_SPEC[imap].VARNAMES[ivar] ;
  char fnam[] = "LOAD_SPECEFF_VAR" ;

  // ------------ BEGIN ------------------

  IVARTYPE = SEARCHEFF_SPEC[imap].IVARTYPE[ivar]; 

  if ( IVARTYPE == IVARTYPE_SPECEFF_REDSHIFT ) {
    return  SEARCHEFF_DATA.REDSHIFT ; 
  }
  else if ( IVARTYPE == IVARTYPE_SPECEFF_PEAKMJD ) {
    return  SEARCHEFF_DATA.PEAKMJD ; 
  }
  else if ( IVARTYPE == IVARTYPE_SPECEFF_DTPEAK ) {
    return  SEARCHEFF_DATA.DTPEAK_MIN ; // smallest abs(T-Tpeak)
  }
  else if ( IVARTYPE == IVARTYPE_SPECEFF_SALT2mB ) {
    return  SEARCHEFF_DATA.SALT2mB ; 
  }
  else if ( IVARTYPE == IVARTYPE_SPECEFF_PEAKMAG ) {
    // take average of peakMags that are OR'ed in this map.
    int NTMP = 0 ;
    MAG = 0.0 ;
    NFILTLIST = SEARCHEFF_SPEC[imap].NFILTLIST_PEAKMAG[ivar] ;
    for(ifilt=0; ifilt < NFILTLIST; ifilt++ ) {
      ifilt_obs = SEARCHEFF_SPEC[imap].IFILTLIST_PEAKMAG[ivar][ifilt] ;
      mag       = SEARCHEFF_DATA.PEAKMAG[ifilt_obs] ;
      /*
      printf(" xxx ifilt=%d of %d,  ifilt_obs=%2d   mag=%.3f \n", 
	     ifilt, NFILTLIST, ifilt_obs, mag); fflush(stdout);
      */
      if ( mag < 100. && mag > 0.0 ) { MAG += mag; NTMP++; }
    }
    if ( NTMP > 0 ) 
      { MAG /= (float)NTMP ; }
    else
      { MAG = -9.0; }
    //    check_magUndefined(mag,varName,fnam);  // abort if mag is undefined
    MAG  += INPUTS_SEARCHEFF.MAGSHIFT_SPECEFF ;
    return(MAG) ;
  }   
  else if ( IVARTYPE == IVARTYPE_SPECEFF_COLOR ) {
    ifilt_obs = SEARCHEFF_SPEC[imap].IFILTOBS_PEAKCOLOR[ivar][0] ;
    mag0 = SEARCHEFF_DATA.PEAKMAG[ifilt_obs]  ;
    check_magUndefined(mag0,varName,fnam);  

    ifilt_obs = SEARCHEFF_SPEC[imap].IFILTOBS_PEAKCOLOR[ivar][1] ;
    mag1 = SEARCHEFF_DATA.PEAKMAG[ifilt_obs] ;
    check_magUndefined(mag1,varName,fnam);  

    color = mag0 - mag1 ;
    return  color ;
  }          

  else if ( IVARTYPE == IVARTYPE_SPECEFF_HOSTMAG ) {
    ifilt_obs = SEARCHEFF_SPEC[imap].IFILTOBS_HOSTMAG[ivar] ;
    mag       = SEARCHEFF_DATA.HOSTMAG[ifilt_obs] ;
    check_magUndefined(mag,varName,fnam);
    return mag ;
  }          
  else if ( IVARTYPE == IVARTYPE_SPECEFF_SBMAG ) {
    ifilt_obs = SEARCHEFF_SPEC[imap].IFILTOBS_SBMAG[ivar] ;
    mag       = SEARCHEFF_DATA.SBMAG[ifilt_obs] ;
    check_magUndefined(mag,varName,fnam);
    return mag ;
  }          


  // if we get here then abort.
  sprintf(c1err,"Could not find ivar=%d for imap=%d ", ivar, imap);
  sprintf(c2err,"Check %s", INPUTS_SEARCHEFF.SPEC_FILE );
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 

  return(-9.0); // should never get here

} // end of LOAD_SPECEFF_VAR



// ***************************************
void assign_SPECEFF(int imap, int ivar, char *VARNAME) {
  
  // Assign SPECEFF VARNAME to one of
  //  * SEARCHEFF_SPEC.IVAR
  //  * SEARCHEFF_SPEC.IVAR_REDSHIFT  
  //  * SEARCHEFF_SPEC.IVAR_PEAKMJD
  //  * SEARCHEFF_SPEC.IVAR_DTPEAK
  //  * SEARCHEFF_SPEC.IFILTOBS_PEAKMAG
  //  * SEARCHEFF_SPEC.IFILTOBS_PEAKCOLOR
  //  * SEARCHEFF_SPEC.IFILTOBS_HOSTMAG 
  //  * SEARCHEFF_SPEC.IFILTOBS_SBMAG
  //
  // If assignment can't be made then abort.
  //
  // Oct 1 2012: check  for PEAKMJD
  //
  // Jan 18 2014: replace IFILTOBS_fromMAP(c0); with INTFILTER function.
  //              Remove obsolete IFILTOBS_fromMAP function.
  //
  // Jun 23 2016: check for HOSTMAG_[band] and SBMAG_[band]
  //
  // Feb 20 2017: check for OR of peakMag bands:
  //                a+b+c+d --> peakMag from any abcd band.
  //
  // Nov 6 2017: set SEARCHEFF_SPEC.IVARTYPE_MASK
  // May 13 2021: allow HOSTMAG or HOST_MAG
  //
  int  ifilt_obs, ifilt2_obs, ifiltlist[MXFILTINDX], LENVAR, ic, i ;
  int  ISPEAKMAG, ISCOLOR, LFF1, LFF2, LMNS ;

  char c0[2], c1[2], c2[2];
  char minus[2] = "-" ;
  char fnam[] = "assign_SPECEFF";

  // ----------- BEGIN ---------

  // check for the easy ones first.
  if ( strcmp(VARNAME,"SPECEFF") == 0 )  {
    SEARCHEFF_SPEC[imap].IVAR = ivar ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF ;
    return ;
  } 
  else if ( strcmp(VARNAME,"REDSHIFT") == 0 )  {
    SEARCHEFF_SPEC[imap].IVAR_REDSHIFT = ivar ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_REDSHIFT ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_REDSHIFT );
    return ;
  }
  else if ( strcmp(VARNAME,"PEAKMJD") == 0 )  {
    SEARCHEFF_SPEC[imap].IVAR_PEAKMJD = ivar ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_PEAKMJD ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_PEAKMJD );
    return ;
  }
  else if ( strcmp(VARNAME,"DTPEAK") == 0 )  {
    SEARCHEFF_SPEC[imap].IVAR_DTPEAK = ivar ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_DTPEAK ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_DTPEAK );
    return ;
  }
  else if ( strcmp(VARNAME,"mB") == 0 )  {
    SEARCHEFF_SPEC[imap].IVAR_SALT2mB = ivar ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] =  IVARTYPE_SPECEFF_SALT2mB ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_SALT2mB );
    return ;
  }

  // check for mag or color.
  // Color must be of the form [filt1]-[filt2]
  // where [filt1] and [filt2] are single-char strings.
  // LENVAR =1 for simple mag; e.g., 'g', 'r'
  // LENVAR =3 for color, e.g., 'g-r'

  ISPEAKMAG  = ISCOLOR = 0;
  LENVAR = strlen(VARNAME);
  ifilt_obs = ifilt2_obs = -9;

  if ( strstr(VARNAME,"+") != NULL ) {
    // OR of multiple peakMags (Feb 2017)
    // Store filter-list for any char that is not a plus (+)
    for(ic=0; ic < LENVAR; ic++ ) {
      sprintf(c0, "%c", VARNAME[ic] );
      if ( c0[0] != '+' ) {
	ifiltlist[ISPEAKMAG] = INTFILTER(c0);
	ISPEAKMAG++ ;
      }
    }
  }
  else if ( LENVAR == 1 ) {
    sprintf(c0, "%c", VARNAME[0] );
    ifilt_obs  = INTFILTER(c0);
    if ( ifilt_obs > 0 ) { ifiltlist[0]=ifilt_obs; ISPEAKMAG  = 1; }
  }
  else if ( LENVAR == 3 ) {
    sprintf(c0, "%c", VARNAME[0] );
    sprintf(c1, "%c", VARNAME[1] );
    sprintf(c2, "%c", VARNAME[2] );

    ifilt_obs  = INTFILTER(c0);
    ifilt2_obs = INTFILTER(c2);

    LMNS  = ( strcmp(c1,minus) == 0  ) ;           // require minus sign
    LFF1  = ( ifilt_obs != ifilt2_obs ) ;          // different filters
    LFF2  = ( ifilt_obs > 0 && ifilt2_obs > 0 ) ;  // valid filters
    if ( LMNS && LFF1 && LFF2 ) { ISCOLOR = 1; }
  }

  if ( ISPEAKMAG ) { 
    SEARCHEFF_SPEC[imap].NFILTLIST_PEAKMAG[ivar]  = ISPEAKMAG ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_PEAKMAG ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_PEAKMAG );
    for(i=0; i < ISPEAKMAG; i++ ) {
      SEARCHEFF_SPEC[imap].IFILTLIST_PEAKMAG[ivar][i] = ifiltlist[i] ;
    }
    return ; 
  }

  if ( ISCOLOR ) { 
    SEARCHEFF_SPEC[imap].IFILTOBS_PEAKCOLOR[ivar][0] = ifilt_obs ;
    SEARCHEFF_SPEC[imap].IFILTOBS_PEAKCOLOR[ivar][1] = ifilt2_obs ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_COLOR ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_COLOR );
    return ; 
  }


  // June 2016: Check for variables of the form XXX_[band]
  ifilt_obs = IFILTOBS_SPECEFF_VAR(VARNAME,"PEAKMAG");
  if ( ifilt_obs >=0 ) {
    SEARCHEFF_SPEC[imap].NFILTLIST_PEAKMAG[ivar]    = 1 ;
    SEARCHEFF_SPEC[imap].IFILTLIST_PEAKMAG[ivar][0] = ifilt_obs ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_PEAKMAG ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_PEAKMAG );
    return ;
  }

  ifilt_obs = IFILTOBS_SPECEFF_VAR(VARNAME,"HOSTMAG HOST_MAG");
  if ( ifilt_obs >=0 ) {
    SEARCHEFF_SPEC[imap].IFILTOBS_HOSTMAG[ivar] = ifilt_obs ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_HOSTMAG ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_HOSTMAG );
    return ;
  }

  ifilt_obs = IFILTOBS_SPECEFF_VAR(VARNAME,"SBMAG");
  if ( ifilt_obs >=0 ) {
    SEARCHEFF_SPEC[imap].IFILTOBS_SBMAG[ivar] = ifilt_obs ;
    SEARCHEFF_SPEC[imap].IVARTYPE[ivar] = IVARTYPE_SPECEFF_SBMAG ;
    SEARCHEFF_SPEC_INFO.IVARTYPE_MASK |= ( 1 << IVARTYPE_SPECEFF_SBMAG );
    return ;
  }

  // ------------------------------------------------------
  // if we get here then abort.
  sprintf(c1err,"Unknown SPECEFF VARNAME: '%s'", VARNAME);
  sprintf(c2err,"Check manual for list of options.");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 

  return ;

} // end of assign_SPECEFF

// ============================================================
int IFILTOBS_SPECEFF_VAR(char *VARNAME, char *PREFIX) {

  // if SPECEFF *VARNAME is of the form *PREFIX_[band],
  // then return integer filter index for [band].
  // If PREFIX is not part of VARNAME, return -9.
  //
  // Example:  VARNAME = HOSTMAG_r and PREFIX = HOSTMAG
  //           --> returns IFILTOBS=3
  //
  // May 13 2021:
  //   Allow multiple prefix names separated by space; e.g., 
  //    PREFIX = "HOSTMAG HOST_MAG"
  //

#define MXPREFIX 4
  int  NPREFIX=0;
  int  i, LENV, IFILTOBS = -9;
  char cfilt[2], PREFIX_LIST[MXPREFIX][80], *ptr_PREFIX;
  char *ptrList[MXPREFIX];
  char sepKey[] = " ";
  //  char fnam[] = "IFILTOBS_SPECEFF_VAR" ;

  // -------------- BEGIN ---------------

  for(i=0; i < MXPREFIX; i++ ) { ptrList[i] = PREFIX_LIST[i]; }

  // check how many space-separated PREFIX names to check
  splitString2(PREFIX, sepKey, MXPREFIX,
	       &NPREFIX, ptrList ); // <== returned

  for ( i=0; i < NPREFIX; i++ ) {
    ptr_PREFIX = PREFIX_LIST[i];
    if ( strstr(VARNAME,ptr_PREFIX) != NULL ) {
      LENV      = strlen(VARNAME);
      sprintf(cfilt,"%c", VARNAME[LENV-1] );
      IFILTOBS  = INTFILTER(cfilt);
    }
  }

  return(IFILTOBS);

} // end IFILTOBS_SPECEFF_VAR

// ============================================================
int  gen_SEARCHEFF_DEBUG(char *WHAT, double RAN, double *EFF) {

  // Created Jan 30 2015
  // DEBUG/hack efficiency function for specialized tests
  // (never for a real analysis).
  //
  // Input *WHAT = "PIPELINE", "SPEC", or "zHOST".
  // Input RAN is a random number from 0 to 1
  // Output *EFF is the efficiency

  char fnam[] = "gen_SEARCHEFF_DEBUG" ;
  int    IDEBUG ;
  double RAN_LOCAL, EFF_LOCAL, zTau ;

  IDEBUG = INPUTS_SEARCHEFF.FUNEFF_DEBUG ;
  if ( IDEBUG == 0 ) { return 0 ; }
  RAN_LOCAL = -9.0 ;

  if ( IDEBUG == 1 )  { 
    RAN_LOCAL = 0.0 ; 
    *EFF      = 1.0 ;
    return 1 ;
  }  
  else if ( IDEBUG == 2 ) {
    { RAN_LOCAL = RAN ; }  // apply effic-fun below
  }
  else {
    sprintf(c1err,"Invalid FUNEFF_DEBUG = %d", 
	    INPUTS_SEARCHEFF.FUNEFF_DEBUG );
    sprintf(c2err,"Must be 1 or 2");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  zTau = 999.0 ;
  if ( strcmp(WHAT,"PIPELINE") == 0 ) 
    { zTau = 0.4 ; }
  else if ( strcmp(WHAT,"SPEC") == 0 ) 
    { zTau = 0.2 ; }
  else if ( strcmp(WHAT,"zHOST") == 0 ) 
    { zTau = 0.3 ; }
  else {
    sprintf(c1err,"Invalid  WHAT input = '%s' ", WHAT);
    sprintf(c2err,"Must be PIPELINE, SPEC or zHOST");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  EFF_LOCAL   = exp( -SEARCHEFF_DATA.REDSHIFT/zTau);
  *EFF        = EFF_LOCAL ;  // load output argument

  //  printf(" xxx EFF(%s) = %f  (RAN=%f) \n", WHAT, EFF_LOCAL, RAN);

  if ( EFF_LOCAL > RAN_LOCAL ) { return 1; } else { return 0; }

} // end of  gen_SEARCHEFF_DEBUG
