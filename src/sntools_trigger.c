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

 Aug 09 2024: add logic to read and apply optional REQUIRE key in SPECEFF map.
 Aug 19 2024: replace c_get[60] with c_get[100] in a few places
 Mar 07 2026: set APPLYMASK_SEARCHEFF_zSPEC 

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
    Begin Refactor Mar 8 2026 [use same SEARCHEFF_MAP struct for SPECID and zHOST maps]

  Inputs:
    SURVEY:       
       name of survey
    APPLYMASK_SEARCHEFF: 
       argument of sim-input APPLY_SEARCHEFF_OPT (only for error checking)


        HISTORY
   
  *************/

  int  NMAP=0 ;
  char  fnam[] = "init_SEARCHEFF"  ;  (void)fnam;
 

  // ------------- BEGIN ----------

  sprintf(BANNER, "%s: Initialize SEARCH EFFICIENCY for '%s' \n", fnam, SURVEY );
  print_banner( BANNER );


  INPUTS_SEARCHEFF.NMAP_DETECT   = 0 ;
  INPUTS_SEARCHEFF.NMAP_PHOTPROB = 0 ;
  INPUTS_SEARCHEFF.NREDUCED_CORR_PHOTPROB = 0 ;
  INPUTS_SEARCHEFF.NPHOTPROB_DUMP = 0 ;

  INPUTS_SEARCHEFF.NMAP_SPECID   = 0 ;
  INPUTS_SEARCHEFF.NMAP_zHOST    = 0 ;

  SEARCHEFF_FLAG        = 0 ;
  SEARCHEFF_LOGIC.NMJD  = 0;

  NONZERO_SEARCHEFF_SPECID  = 0 ;
  NONZERO_SEARCHEFF_zHOST   = 0 ;

  sprintf(PATH_SEARCHEFF, "%s %s/models/searcheff", 
	  PATH_USER_INPUT, PATH_SNDATA_ROOT );

  if ( INPUTS_SEARCHEFF.FUNEFF_DEBUG ) {
    printf("\t Use FUNEFF_DEBUG = %d \n", INPUTS_SEARCHEFF.FUNEFF_DEBUG );
  }

  // read single file to get pipeline efficiency vs. SNR or vs. MAG.
  // Returns number of maps that have an efficiency curve defined.
  NMAP = init_SEARCHEFF_PIPELINE(SURVEY) ;
  if ( NMAP > 0 ) { init_SEARCHEFF_LOGIC(SURVEY); }  // read detection logic


  // init prob of getting specID -confirmation
  init_SEARCHEFF_SPECID(SURVEY);


  // init prob of getting zHOST for unconfirmed SN
  init_SEARCHEFF_zHOST(SURVEY);
  
  // Mar 2018: check that user SEARCHEFF mask is possible
  check_APPLYMASK_SEARCHEFF(SURVEY, APPLYMASK_SEARCHEFF);



}  // end of init_SEARCHEFF


// *************************************
void  check_APPLYMASK_SEARCHEFF(char *SURVEY, int APPLYMASK_SEARCHEFF_USER) {

  // Created March 2018
  // Compute SEARCHEFF-mask that is possible, 
  // and check that it is compatible with user-input
  // APPLYMASK_SEARCHEFF that is the argument of APPLY_SEARCHEFF_OPT.
  // Abort if APPLYMASK_SEARCHEFF_USER is impossible to obtain.
  //
  // Mar 7 2026: accommodate logic for zSPEC (8-bit)
  // Mar 9 2026: refactored 

  int  APPLYMASK_ALLOWED=0 , OVP ;
  int  IFLAG_SPECID_EFFZERO = SEARCHEFF_INFO_SPECID.IFLAG_EFFZERO;
  int  IFLAG_SPECID_EFFONE  = SEARCHEFF_INFO_SPECID.IFLAG_EFFONE ;  
  int  NMAP_SPECID          = SEARCHEFF_INFO_SPECID.NMAP ;
  int  LSPECID   = IFLAG_SPECID_EFFZERO || IFLAG_SPECID_EFFONE || NMAP_SPECID>0 ;

  int  NONZERO_SPECID = SEARCHEFF_INFO_SPECID.NONZERO_SEARCHEFF;
  int  NONZERO_zHOST  = SEARCHEFF_INFO_zHOST.NONZERO_SEARCHEFF;

  char fnam[] = "check_APPLYMASK_SEARCHEFF" ;  (void)fnam;

  // ---------- BEGIN --------------

  // pipeline detection always allowed.
  APPLYMASK_ALLOWED += APPLYMASK_SEARCHEFF_PIPELINE ;

  if ( NONZERO_SPECID )  { 
    APPLYMASK_ALLOWED += APPLYMASK_SEARCHEFF_SPECID ; 
    APPLYMASK_ALLOWED |= APPLYMASK_SEARCHEFF_zSPEC ;  // Mar 2026
  }

  if ( NONZERO_zHOST && LSPECID ) {
    APPLYMASK_ALLOWED += APPLYMASK_SEARCHEFF_zHOST; 
    APPLYMASK_ALLOWED |= APPLYMASK_SEARCHEFF_zSPEC ;  // Mar 2026
  }

  OVP = ( APPLYMASK_ALLOWED & APPLYMASK_SEARCHEFF_USER) ;

  if ( OVP != APPLYMASK_SEARCHEFF_USER) {
    print_preAbort_banner(fnam);
    printf("\t APPLY_SEARCHEFF_OPT += %d --> detection pipeline.\n",
	   APPLYMASK_SEARCHEFF_PIPELINE);

    printf("\t APPLY_SEARCHEFF_OPT += %d --> SPECID confirmed.\n",
	   APPLYMASK_SEARCHEFF_SPECID );

    printf("\t APPLY_SEARCHEFF_OPT += %d --> zHOST(spec).\n",
	   APPLYMASK_SEARCHEFF_zHOST );

    printf("\t APPLY_SEARCHEFF_OPT += %d --> zSPEC(SN or HOST).\n",
	   APPLYMASK_SEARCHEFF_zSPEC );

    printf("\t NONZERO_SPECID=%d   LSPECID=%d \n",  NONZERO_SPECID, LSPECID);
    printf("\t NONZERO_zHOST=%d \n", NONZERO_zHOST);

    sprintf(c1err,"Invalid user-input 'APPLY_SEARCHEFF_OPT: %d' (OVP=%d) ",
	    APPLYMASK_SEARCHEFF_USER, OVP) ;
    sprintf(c2err,"Simulation can produce SEARCHEFF mask %d \n",
	    APPLYMASK_ALLOWED);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  // - - - - - - 
  // Sep 4 2019: prepare README comment for trigger.

  COMMENT_README_SEARCHEFF[0][0] = 0;
  COMMENT_README_SEARCHEFF[1][0] = 0;

  sprintf(COMMENT_README_SEARCHEFF[0], "APPLY_SEARCHEFF_OPT=%d --> ",
	  APPLYMASK_SEARCHEFF_USER);

  int REQ_PIPE, REQ_SPEC, REQ_zHOST, NUSE=0;
  char ctmp[100], clist[100], cplus[] = "+" ;
  REQ_PIPE  = (APPLYMASK_SEARCHEFF_USER & APPLYMASK_SEARCHEFF_PIPELINE );
  REQ_SPEC  = (APPLYMASK_SEARCHEFF_USER & APPLYMASK_SEARCHEFF_SPECID );
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

  strcat(COMMENT_README_SEARCHEFF[0],ctmp);

  if ( INPUTS_SEARCHEFF.APPLY_DETECT_SINGLE ) {
    sprintf(COMMENT_README_SEARCHEFF[1],
	    "SPECIAL OPTION: "
	    "EFF(PIPE) applied to single-exposures instead of co-add");
  }
  
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
  FILE *fp = NULL;
  int   REQUIRE_EFF_FILE=0, gzipFlag, imap, NMAP=0 ;
  int   FOUNDMAP_DETECT=0, FOUNDMAP_PHOTPROB=0 ;
  char  file_local[MXPATHLEN], c_get[100], *ptrFile_user, *ptrFile_final ;   
  char  fnam[] = "init_SEARCHEFF_PIPELINE"  ;  (void)fnam;
    

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
    REQUIRE_EFF_FILE = 0 ;
  }
  else if ( IGNOREFILE(ptrFile_user) ) {
    sprintf(file_local, "%s",  ptrFile_user); 
    REQUIRE_EFF_FILE = 0 ;
  }
  else {
    sprintf(file_local, "%s",  ptrFile_user); 
    REQUIRE_EFF_FILE = 1 ;
  }


  // use utility to check local dir and path.
  if ( REQUIRE_EFF_FILE ) {
    fp = snana_openTextFile(OPTMASK, PATH_SEARCHEFF, file_local, 
			    ptrFile_final, &gzipFlag ); // returned
  }

  
  if ( fp == NULL ) { 

    if ( REQUIRE_EFF_FILE ) {
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
  // Function returns 1 if map is found; 0 otherwise
  //
  // Nov 30 2022: read FIELD key
  //

  double VAL, EFF;
  int  IREAD, imap, NBIN ;
  char ctmp[MXFILTINDX] ;
  char fnam[] = "readMap_SEARCHEFF_DETECT";  (void)fnam;

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
    readchar(fp,ctmp);
    sprintf(SEARCHEFF_DETECT[imap].FILTERLIST, "%s", ctmp);
    return(1);
  }
  
  if ( strcmp(key,"FIELD:")==0 ) {           // Nov 30 2022
    imap = INPUTS_SEARCHEFF.NMAP_DETECT-1 ; 
    readchar(fp,ctmp);
    sprintf(SEARCHEFF_DETECT[imap].FIELDLIST, "%.60s", ctmp);
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
  char fnam[] = "malloc_NEXTMAP_SEARCHEFF_DETECT" ;  (void)fnam;
  // ------------- BEGIN -------------
  INPUTS_SEARCHEFF.NMAP_DETECT++ ; 
  imap  = INPUTS_SEARCHEFF.NMAP_DETECT-1;
  MEMD  = MXROW_SEARCHEFF_DETECT * sizeof(double);
  SEARCHEFF_DETECT[imap].VAL = (double*)malloc(MEMD);
  SEARCHEFF_DETECT[imap].EFF = (double*)malloc(MEMD);
  SEARCHEFF_DETECT[imap].FILTERLIST[0] = 0 ;
  sprintf(SEARCHEFF_DETECT[imap].FIELDLIST, "%s", ALL);

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
  char fnam[] = "readMap_SEARCHEFF_PHOTPROB" ; (void)fnam;

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
    sprintf( SEARCHEFF_PHOTPROB[imap].FIELDLIST,  "%s", ALL) ;
    sprintf( SEARCHEFF_PHOTPROB[imap].FILTERLIST, "%s", ALL) ;

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
  char fnam[] = "IVARABS_SEARCHEFF_PHOTPROB" ; (void)fnam;

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

  int     ibin, i, IBIN_HALF, IBIN_ONE ;
  double  VAL, EFF, EFFDIF, EFFDIF_ONE, EFFDIF_HALF ;

  char  *ptr_effname, *MAPNAME, cline[MXPATHLEN] ;
  char fnam[] = "check_SEARCHEFF_DETECT" ; (void)fnam;

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


  EFFDIF_ONE = EFFDIF_HALF = 999.0;
  IBIN_HALF  = IBIN_ONE = 0;

  for ( ibin=0; ibin < NBIN; ibin++ ) {

    VAL =  SEARCHEFF_DETECT[imap].VAL[ibin] ;
    EFF =  SEARCHEFF_DETECT[imap].EFF[ibin] ;

    // keep track of bin in which efficiency is closest to 1 or .5
    EFFDIF = fabs(1.0-EFF) ;
    if ( EFFDIF < EFFDIF_ONE ) 
      { EFFDIF_ONE = EFFDIF ; IBIN_ONE = ibin ; }
      
    EFFDIF = fabs(0.5-EFF) ;
    if ( EFFDIF < EFFDIF_HALF ) 
      { EFFDIF_HALF = EFFDIF ; IBIN_HALF = ibin ; 
	// printf(" xxx imap=%d  bin=%d  VAL=%.3f  EFF=%.4f  EFFDIF=%f \n", 
	//   imap, ibin, VAL, EFF, EFFDIF ); fflush(stdout);
      }    

  } // NBIN 

  
  // ------
  // print VAL (SNR or MAG) when eff is closest to  0.5

  VAL     = SEARCHEFF_DETECT[imap].VAL[IBIN_HALF];
  EFF     = SEARCHEFF_DETECT[imap].EFF[IBIN_HALF];
  MAPNAME = SEARCHEFF_DETECT[imap].MAPNAME;
  sprintf(cline, "\t EFF_DETECT(%s) = %5.3f at %s = %5.2f (%s)", 
	  cfilt, EFF, ptr_effname, VAL, MAPNAME );

  printf("%s\n", cline); fflush(stdout);
  i = SEARCHEFF_DETECT[imap].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[imap].README[i], "%s", cline);
  SEARCHEFF_DETECT[imap].NLINE_README++ ;

} // end of check_SEARCHEFF_DETECT


// *********************************************
void check_SEARCHEFF_PHOTPROB(int imap) {
  
  // set README comment(s).
  char NAME[100], FIELDLIST[100], BANDLIST[MXFILTINDX];
  char fnam[] = "check_SEARCHEFF_PHOTPROB" ; (void)fnam;

  // -------------- BEGIN ----------------

  if ( imap==0 ) {
    printf("\n");
    printf("        PHOTPROB \n");
    printf("         MAPNAME             FIELD(s)     BAND(s)   "
	   "NFUN   NROW \n");
    printf("# -------------------------------"
	   "------------------------------ \n");
  }

  sprintf(NAME,      "%s", SEARCHEFF_PHOTPROB[imap].NAME);
  sprintf(FIELDLIST, "%s",  SEARCHEFF_PHOTPROB[imap].FIELDLIST);
  sprintf(BANDLIST,  "%s",  SEARCHEFF_PHOTPROB[imap].FILTERLIST);
  int  NFUN       = SEARCHEFF_PHOTPROB[imap].NFUN_CDF ;
  int  NROW       = SEARCHEFF_PHOTPROB[imap].NROW;

  printf("  %d) %16.16s   %10.10s   %10.10s       %d   %4d \n",
	 imap, NAME, FIELDLIST, BANDLIST, NFUN, NROW );
  fflush(stdout) ;

  SEARCHEFF_PHOTPROB[imap].NLINE_README = 0 ; 
  int i = SEARCHEFF_PHOTPROB[imap].NLINE_README ;

  if ( i < 4 ) {
    char *cptr = SEARCHEFF_PHOTPROB[imap].README[i] ;
    sprintf(cptr, "  PHOTPROB MAP %.40s: BANDS=%.80s  FIELD=%.40s \n",
	    NAME, BANDLIST, FIELDLIST );
  }

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
     cline[MXPATHLEN+40]
    ,logic[100]
    ,c_get[100]
    ,surveykey[60]
    ,logicFile_Default[] = "SEARCHEFF_PIPELINE_LOGIC.DAT"
    ,logicFile[MXPATHLEN]
    ,*ptrFile_user
    ,*ptrFile_final
    ;

  char fnam[] = "init_SEARCHEFF_LOGIC";   (void)fnam;
   

  // -------------- BEGIN ----------------


  ptrFile_user  = INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE ;
  ptrFile_final = INPUTS_SEARCHEFF.PIPELINE_LOGIC_FILE ;

  if ( strcmp(ptrFile_user,"DEFAULT") == 0 ) 
    { sprintf(logicFile,"%s", logicFile_Default); }
  else
    { sprintf(logicFile,"%s", ptrFile_user);  }

  
  fp = snana_openTextFile(OPTMASK, PATH_SEARCHEFF, logicFile,
			  ptrFile_final, &gzipFlag ); // returned

  if ( !fp ) {
    abort_openTextFile("SEARCHEFF_PIPELINE_LOGIC_FILE",
		       PATH_SEARCHEFF, logicFile, fnam );
  }


  sprintf(cline, "\n   Fetch SOFTWARE SEARCH-LOGIC from : "); 
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", cline);
  printf("%s\n", cline);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;


  sprintf(cline, "\t %s", ptrFile_final ); 
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", cline);
  printf("%s\n", cline); 
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;

  
  fflush(stdout);

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
    sprintf(c2err,"Check %.*s", MXCHAR_MSGERR, ptrFile_final ) ;
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
  char fnam[] = "parse_search_eff_logic" ; (void)fnam;

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
  printf("%s\n", c1err);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;

  sprintf(c1err, "\t Trigger epoch contains all obs within %.3f days",
	  INPUTS_SEARCHEFF.TIME_SINGLE_DETECT);
  i = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ;
  sprintf(SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[i], "%s", c1err);
  printf("%s\n", c1err);
  SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README++ ;

  fflush(stdout);
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
void  init_SEARCHEFF_SPECID(char *survey) {

  // Created July 2011 by R.Kessler
  // If optional EFF(SPEC-confirm) table exists, read and store.
  // Table name must be SEARCHEFF_[survey]_SPEC.DAT;
  // this table can reside in current [user] directory
  // or in PATH_SEARCHEFF
  //
  // Refactor March 2026
  //
  //
  // If user-requested file (SEARCHEFF_SPEC_FILE: <file>)
  // does not exist then abort. If default SPEC-eff file
  // does not exist then just return with message that EFFSPEC=1.
  //

  char fnam[] = "init_SEARCHEFF_SPECID";  (void)fnam;

  // ------------ BEGIN --------------

  init_searcheff_map(MAPTYPE_SEARCHEFF_SPECID, &SEARCHEFF_INFO_SPECID) ;

  read_searcheff_map(INPUTS_SEARCHEFF.USER_SPEC_FILE, &SEARCHEFF_INFO_SPECID) ;


  // temp assignment for snlc_sim ... can remove this later after REFACTOR is released
  INPUTS_SEARCHEFF.NMAP_SPECID = SEARCHEFF_INFO_SPECID.NMAP; 


} // end of init_SEARCHEFF_SPECID


// **************************************
void init_searcheff_map(char *MAPTYPE, SEARCHEFF_INFO_DEF *SEARCHEFF_INFO) {

  // Created Mar 2026 [as part of refactor]
  // init generic INFO including multi-Dimensional MAP_LIST

  int ivar, imap;
  char fnam[] = "init_searcheff_map" ; (void)fnam;

  // ------------ BEGIN ------------

  printf(" -------------------------------------------------- \n");
  printf("  Init SEARCHEFF_%s Efficiency map(s): \n", MAPTYPE); fflush(stdout);

  sprintf(SEARCHEFF_INFO->MAPTYPE, "%s", MAPTYPE);

  SEARCHEFF_INFO->NMAP                   = 0 ;
  SEARCHEFF_INFO->OPT_EXTRAP             = 0 ;
  SEARCHEFF_INFO->OPT_FIELDMATCH_REQUIRE = 1 ;
  SEARCHEFF_INFO->IVARTYPE_MASK     = 0 ;
  SEARCHEFF_INFO->FLAG_PEAKMAG_ONLY = 0 ;
  SEARCHEFF_INFO->BOOLEAN_OR        = 1 ; // default logic is OR
  SEARCHEFF_INFO->BOOLEAN_AND       = 0 ; 

  SEARCHEFF_INFO->NLINE_README      = 0;

  SEARCHEFF_INFO->NONZERO_SEARCHEFF = 0;
  SEARCHEFF_INFO->IFLAG_EFFZERO          = 0 ;
  SEARCHEFF_INFO->IFLAG_EFFONE           = 0 ;

  for ( imap=0; imap < MXMAP_SEARCHEFF_MAP;  imap++ ) {

    SEARCHEFF_INFO->MAP_LIST[imap].NVAR_TOT       =  0 ; 
    SEARCHEFF_INFO->MAP_LIST[imap].NVAR_HOST      =  0 ; 
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR           = -9 ; 
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_REDSHIFT  = -9 ;
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_PEAKMJD   = -9 ;
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_DTPEAK    = -9 ;
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_DTSEASON_PEAK = -9 ;
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_LOGMASS   = -9 ;    // Aug 27 2024
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_SALT2mB   = -9 ;
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_SALT2x1   = -9 ;     // Aug 27 2024    
    SEARCHEFF_INFO->MAP_LIST[imap].IVAR_SALT2c    = -9 ;     // Aug 27 2024
    SEARCHEFF_INFO->MAP_LIST[imap].REQUIRE        = 0 ; 
    SEARCHEFF_INFO->MAP_LIST[imap].MAGSHIFT       = 0.0 ;

    sprintf(SEARCHEFF_INFO->MAP_LIST[imap].FIELDLIST, "%s", ALL ); // default is all fields

    
    for ( ivar=0 ; ivar < MXVAR_SEARCHEFF_MAP; ivar++ ) {
      SEARCHEFF_INFO->MAP_LIST[imap].IVARTYPE[ivar]              = -9 ;

      SEARCHEFF_INFO->MAP_LIST[imap].FLAG_MAG[ivar]   = 0 ;
      SEARCHEFF_INFO->MAP_LIST[imap].SHIFT[ivar]      = 0 ;

      SEARCHEFF_INFO->MAP_LIST[imap].NFILTLIST_PEAKMAG[ivar]     =  0 ;
      SEARCHEFF_INFO->MAP_LIST[imap].IFILTLIST_PEAKMAG[ivar][0]  = -9 ;
      SEARCHEFF_INFO->MAP_LIST[imap].IFILTLIST_PEAKMAG[ivar][1]  = -9 ;

      SEARCHEFF_INFO->MAP_LIST[imap].NFILTLIST_HOSTMAG[ivar]     =  0 ;
      SEARCHEFF_INFO->MAP_LIST[imap].IFILTLIST_HOSTMAG[ivar][0]  = -9 ;
      SEARCHEFF_INFO->MAP_LIST[imap].IFILTLIST_HOSTMAG[ivar][1]  = -9 ;

      SEARCHEFF_INFO->MAP_LIST[imap].NFILTLIST_SBMAG[ivar]     =  0 ;
      SEARCHEFF_INFO->MAP_LIST[imap].IFILTLIST_SBMAG[ivar][0]  = -9 ;
      SEARCHEFF_INFO->MAP_LIST[imap].IFILTLIST_SBMAG[ivar][1]  = -9 ;

      SEARCHEFF_INFO->MAP_LIST[imap].IVAR_HOSTLIB[ivar] = -9 ;
    }
    
  } // end imap


  return;

} // end init_searcheff_map


void  init_searcheff_shifts(int imap, SEARCHEFF_INFO_DEF *SEARCHEFF_INFO) {

  // Created Mar 20 2026
  // Define systematic shift per map variable (default shift =0).
  // See sim input keys
  //   SEARCHEFF_SPEC_SHIFT([varname]): [value]
  //   SEARCHEFF_zHOST_SHIFT([varname]): [value]
  //

  char *MAPTYPE = SEARCHEFF_INFO->MAPTYPE ;
  int    ivar, ivar2, NSHIFT ;
  double MAGSHIFT ;
  double *ptr_SHIFT_VALUES;
  char   *ptr_SHIFT_VARNAMES[MXVAR_SEARCHEFF_MAP] ;
  SEARCHEFF_MAP_DEF *MAP = &SEARCHEFF_INFO->MAP_LIST[imap];
  char fnam[] = "init_searcheff_shifts" ;  (void)fnam;

  // ------------ BEGIN -----------

  // select inputs based on SPEC or zHOST map type
  if ( strcmp(MAPTYPE,MAPTYPE_SEARCHEFF_SPECID)==0 ) {
    NSHIFT           = INPUTS_SEARCHEFF.NSHIFT_SPEC ;
    MAGSHIFT         = INPUTS_SEARCHEFF.MAGSHIFT_SPECEFF; 
    ptr_SHIFT_VALUES = INPUTS_SEARCHEFF.SHIFT_VALUES_SPEC ;
    for(ivar=0; ivar < NSHIFT; ivar++ ) 
      { ptr_SHIFT_VARNAMES[ivar] = INPUTS_SEARCHEFF.SHIFT_VARNAMES_SPEC[ivar]; }
  }
  else if ( strcmp(MAPTYPE,MAPTYPE_SEARCHEFF_zHOST)==0 ) {
    NSHIFT           = INPUTS_SEARCHEFF.NSHIFT_zHOST ;
    MAGSHIFT         = INPUTS_SEARCHEFF.MAGSHIFT_zHOSTEFF; 
    ptr_SHIFT_VALUES = INPUTS_SEARCHEFF.SHIFT_VALUES_zHOST ;
    for(ivar=0; ivar < NSHIFT; ivar++ ) 
      { ptr_SHIFT_VARNAMES[ivar] = INPUTS_SEARCHEFF.SHIFT_VARNAMES_zHOST[ivar]; }
  }
  else {
    NSHIFT   = 0;
    MAGSHIFT = 0.0 ;
  }

  // assign same MAGSHIFT to all mags (later applied only to MAG type map);
  MAP->MAGSHIFT = MAGSHIFT ; 

  // - - - - - - -
  /* xxx mark delete May 4 2026 ... SHIFT unclear ?? xxxxx
  if ( MAGSHIFT != 0.0 ) {
    printf("\t set SHIFT = %8.2f for %s PEAK & HOST mags \n",
	   SHIFT, MAPTYPE ); fflush(stdout);
  }    
  xxxxxxxxxx end mark xxxxxxx  */



  if ( NSHIFT == 0 ) { return ; }  

  char *VARNAME, *VARNAME2;
  double SHIFT;
  // assign more specific shift based on VARNAME
  for(ivar=0; ivar < MAP->NVAR_TOT-1; ivar++ ) {
    VARNAME = MAP->VARNAMES[ivar];
    for(ivar2=0; ivar2 < NSHIFT; ivar2++ ) {
      VARNAME2 = ptr_SHIFT_VARNAMES[ivar2] ;
      SHIFT    = ptr_SHIFT_VALUES[ivar2];
      if ( strcmp(VARNAME,VARNAME2) == 0 ) {
	MAP->SHIFT[ivar] = SHIFT;
	printf("\t set SHIFT = %8.2f for %s VARNAME = %s \n",
	       SHIFT, MAPTYPE, VARNAME); fflush(stdout);
      }
    }
  }

  return;
} // end init_searcheff_shifts

void read_searcheff_map(char *USER_MAP_FILE, SEARCHEFF_INFO_DEF *SEARCHEFF_INFO) {

  // Created Mar 2026
  // Inputs:
  //    MAPTYPE:   type of map; SPECID or zHOST
  //    USER_MAP_FILE : full path name OR base name OR NONE or ZERO or ONE
  //
  // Outputs:
  //    SEARCHEFF_INFO: typedef struct with info per map
  //

  char *MAPTYPE = SEARCHEFF_INFO->MAPTYPE;
  FILE *fp = NULL ;

  int   NMAP, NROW, ivar, NVAR, FOUND_VARNAMES ;
  int   ID, NDIM, NFUN, N, gzipFlag ;
  int   REQUIRE_EFF_FILE = 0, REQUIRE_MAP, OPT_EXTRAP=0;
  char  KEY_ROW[20], KEYNAME_MAP_FILE[40], KEYNAME_EFF[40];
  char  KEY_STOP[] = "" ;
  char  *VARNAME, *VARLIST, FIELDLIST[100] ;
  char  eff_file_local[MXPATHLEN], *ptrFile_user, *ptrFile_final, *cptr, *fg;
  char  c_get[100], LINE[MXPATHLEN];
  double PEAKMJD_RANGE[2];
  int   OPTMASK       = INPUTS_SEARCHEFF.OPTMASK_OPENFILE ;

  int IDGRIDMAP_OFFSET=-99, IVARTYPE_MASK = 0  ;
  char fnam[] = "read_searcheff_map" ;  (void)fnam;

  // ----------- BEGIN -------------

  
  if ( strcmp(MAPTYPE,MAPTYPE_SEARCHEFF_SPECID)==0 )  { 
    IDGRIDMAP_OFFSET = IDGRIDMAP_SPECEFF_OFFSET ; 
    sprintf(KEY_ROW,"SPECEFF:");
  }
  else if ( strcmp(MAPTYPE,MAPTYPE_SEARCHEFF_zHOST)==0 ) { 
    IDGRIDMAP_OFFSET = IDGRIDMAP_zHOST_OFFSET ; 
    sprintf(KEY_ROW,"HOSTEFF:");
  }
  else {
    sprintf(c1err,"Unknown MAPTYPE = '%s' (expect %s or  %s)", 
	    MAPTYPE, MAPTYPE_SEARCHEFF_SPECID, MAPTYPE_SEARCHEFF_zHOST);
    sprintf(c2err,"for USER_MAP_FILE = %s", USER_MAP_FILE); 
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;   
  }

  sprintf(KEYNAME_MAP_FILE,"SEARCHEFF_%s_FILE", MAPTYPE); // user-input key name
  sprintf(KEYNAME_EFF,     "EFF_%s",            MAPTYPE); // internal string for prints

  
  // define effspec filename form user input or from default name

  N = 0 ;  cptr = SEARCHEFF_INFO->README[N] ;
  N++;  SEARCHEFF_INFO->NLINE_README = N ;  

  ptrFile_user  = USER_MAP_FILE ;             // passed from user input
  ptrFile_final = SEARCHEFF_INFO->MAP_FILE;   // full path of eff file if defined


  if ( IGNOREFILE(ptrFile_user) ) {
    SEARCHEFF_INFO->IFLAG_EFFONE = 1;
    sprintf(cptr, "\t No %s ==> Force %s = 1",  KEYNAME_MAP_FILE, KEYNAME_EFF);
    printf("\n %s\n", cptr);    fflush(stdout);
    SEARCHEFF_INFO->NONZERO_SEARCHEFF++ ;  // pretend there is a map
  }

  else if ( strcmp(ptrFile_user,"ONE") == 0  ) {
    SEARCHEFF_INFO->IFLAG_EFFONE = 1;
    sprintf(cptr,"\t %s = ONE -> Force %s = 1. \n", KEYNAME_MAP_FILE, KEYNAME_EFF);
    printf("\n %s\n", cptr);    fflush(stdout);
    SEARCHEFF_INFO->NONZERO_SEARCHEFF++ ;  // pretend there is a map
  }

  else if ( strcmp(ptrFile_user,"ZERO") == 0  ) {
    SEARCHEFF_INFO->IFLAG_EFFZERO = 1;
    sprintf(cptr,"\t %s = ZERO -> Force %s = 0 (refac) \n", 
	    KEYNAME_MAP_FILE, KEYNAME_EFF);
    printf("\n %s\n", cptr);    fflush(stdout);
  }
  else { 
    sprintf(eff_file_local, "%s",  ptrFile_user); 
    REQUIRE_EFF_FILE = 1;
  }

  // - - - - -
  // use utility to check local dir and path.
  if ( REQUIRE_EFF_FILE ) {
    fp = snana_openTextFile(OPTMASK, PATH_SEARCHEFF, eff_file_local, 
			    ptrFile_final, &gzipFlag ); // returned

    if ( !fp ) {
      abort_openTextFile(KEYNAME_MAP_FILE,
			 PATH_SEARCHEFF, eff_file_local, fnam );
    }
  }
  else {
    return ;
  }

  printf("\n   Reading %s efficiency map from \n\t %s\n", MAPTYPE, ptrFile_final );

  NVAR = NROW = NMAP = 0;

  REQUIRE_MAP = 0 ;  sprintf(FIELDLIST,"%s", ALL);
  PEAKMJD_RANGE[0] = 10000;  PEAKMJD_RANGE[1] = 90000;  

  sprintf(c2err, "Check %s-eff file above", MAPTYPE );


  while( (fscanf(fp, "%s", c_get )) != EOF) {
    
    if ( commentchar(c_get) ) { fg = fgets(LINE, 100, fp ); continue ; }

    if ( strcmp(c_get,"NVAR:")==0 ) { warn_NVAR_KEY(ptrFile_final); }

    if ( strcmp(c_get,"OPT_EXTRAP:") == 0 ) 
      { readint(fp, 1, &SEARCHEFF_INFO->OPT_EXTRAP); }

    if ( strcmp(c_get,"OPT_FIELDMATCH_REQUIRE:") == 0 ) 
      { readint(fp, 1, &SEARCHEFF_INFO->OPT_FIELDMATCH_REQUIRE); }  

    // - - - - - -
    // check optional key to associate map with particular field[list]
    if ( strcmp(c_get,"FIELD:")==0 || strcmp(c_get,"FIELDLIST:")==0 ) 
      { readchar(fp, FIELDLIST); }

    if ( strcmp(c_get,"PEAKMJD:") == 0 || strcmp(c_get,"PEAKMJD_RANGE:")==0 ) 
      { readdouble(fp, 2, PEAKMJD_RANGE);}

    if ( strcmp(c_get,"REQUIRE:") == 0 ) 
      { readint(fp, 1, &REQUIRE_MAP); }
    
    FOUND_VARNAMES = ( strcmp(c_get,"VARNAMES:")==0 ) ;

    if ( FOUND_VARNAMES && NMAP < MXMAP_SEARCHEFF_MAP ) {
      fg = fgets(LINE, 100, fp ); // scoop up varnames

      NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING, LINE, fnam );

      SEARCHEFF_INFO->MAP_LIST[NMAP].NVAR_TOT = NVAR ;
      SEARCHEFF_INFO->MAP_LIST[NMAP].REQUIRE  = REQUIRE_MAP;      
      SEARCHEFF_INFO->MAP_LIST[NMAP].PEAKMJD_RANGE[0] = PEAKMJD_RANGE[0];
      SEARCHEFF_INFO->MAP_LIST[NMAP].PEAKMJD_RANGE[1] = PEAKMJD_RANGE[1];
      sprintf(SEARCHEFF_INFO->MAP_LIST[NMAP].FIELDLIST,"%s", FIELDLIST );

      for ( ivar=0; ivar < NVAR; ivar++ ) {
	VARNAME = SEARCHEFF_INFO->MAP_LIST[NMAP].VARNAMES[ivar] ;
	VARLIST = SEARCHEFF_INFO->MAP_LIST[NMAP].GRIDMAP.VARLIST ;
	get_PARSE_WORD(0, ivar, VARNAME, fnam);

	IVARTYPE_MASK |= assign_MAP_VARNAME(MAPTYPE, ivar, VARNAME, 
					    &SEARCHEFF_INFO->MAP_LIST[NMAP] ); 

	catVarList_with_sep(VARLIST, VARNAME, " ");
      }

      SEARCHEFF_INFO->IVARTYPE_MASK = IVARTYPE_MASK;

      OPT_EXTRAP = SEARCHEFF_INFO->OPT_EXTRAP;
      ID = IDGRIDMAP_OFFSET + NMAP;   NDIM = NVAR-1; NFUN=1;
      read_GRIDMAP(fp, KEY_ROW, KEY_ROW, KEY_STOP, ID, NDIM, NFUN, OPT_EXTRAP,
		   MXROW_SEARCHEFF_MAP, fnam,
		   &SEARCHEFF_INFO->MAP_LIST[NMAP].GRIDMAP ) ;
      
      printf("\t for FIELDLIST=%s   and  %.0f < PEAKMJD < %.0f \n", 
	     FIELDLIST, PEAKMJD_RANGE[0], PEAKMJD_RANGE[1]) ;

      if ( REQUIRE_MAP ) 
	{ printf("\t This %s-EFF map is required for event selection.\n", MAPTYPE); }
	
      init_searcheff_shifts(NMAP, SEARCHEFF_INFO);	

      SEARCHEFF_INFO->NONZERO_SEARCHEFF++; 
      NMAP++;  // increment map counter

      // reset map attributes
      REQUIRE_MAP = 0 ;  sprintf(FIELDLIST,"%s", ALL);
      PEAKMJD_RANGE[0] = 10000;  PEAKMJD_RANGE[1] = 90000;  

    } // end FOUND_VARNAMES

  }  // end of read-line loop

  // - - - - - - - - - - - - - - - - - - - - - 
  fclose(fp); // done reading EFF map(s)

  SEARCHEFF_INFO->NMAP       = NMAP ;

  printf("\n\t Finished preparing %d %s maps\n", 
	 NMAP, MAPTYPE);
  printf("\t OPT_EXTRAP             = %d \n", 
	 SEARCHEFF_INFO->OPT_EXTRAP );
  printf("\t OPT_FIELDMATCH_REQUIRE = %d \n", 
	 SEARCHEFF_INFO->OPT_FIELDMATCH_REQUIRE );
  fflush(stdout) ;

  if ( NMAP > MXMAP_SEARCHEFF_MAP ) {
    sprintf(c1err,"NMAP=%d  exceeds  MXMAP_SEARCHEFF=%d", NMAP, MXMAP_SEARCHEFF_MAP );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }



  // -------------------------------------------------------
  // Nov 2017: check if this map depends ONLY on PEAKMAG;
  //   if so, set flag to make sim run faster.
  int IVARTYPE_MASK_PEAKMAG = 0, OVP;
  IVARTYPE_MASK_PEAKMAG += (1<< IVARTYPE_EFFMAP_PEAKMAG ) ;
  IVARTYPE_MASK_PEAKMAG += (1<< IVARTYPE_EFFMAP_PEAKCOLOR   ) ;
  OVP = ( IVARTYPE_MASK_PEAKMAG & SEARCHEFF_INFO->IVARTYPE_MASK ) ;
  if ( SEARCHEFF_INFO->IVARTYPE_MASK == OVP ) 
    { SEARCHEFF_INFO->FLAG_PEAKMAG_ONLY=1 ; }

  fflush(stdout) ;

  N = 0 ;
  cptr = SEARCHEFF_INFO->README[N] ;  N++ ;  
  sprintf(cptr, "\t %s", "GRID-MAP read from ");
  cptr = SEARCHEFF_INFO->README[N] ;  N++ ;  
  sprintf(cptr, "\t %s", ptrFile_final );

  /* xxx mark delete Mar 9 2026 xxxxxxxx
  if ( fabs(SPECEFF_SCALE-1.0) > 0.0001 ) {
    cptr = SEARCHEFF_SPECID_INFO.README[N] ;    N++ ;  
    sprintf(cptr, "\t EffSpec from file is scaled by %le", SPECEFF_SCALE );
  }
  xxxx */

  SEARCHEFF_INFO->NLINE_README = N ;

  (void)fg; 

  return ;

} // end read_searcheff_map


// ***************************************
int assign_MAP_VARNAME(char *MAPTYPE, int ivar, char *VARNAME, 
			SEARCHEFF_MAP_DEF *MAP) {
  
  // Refactored Mar 2026 to handle both SPEC and zHOST maps.
  // 
  // Inputs:
  //   MAPTYPE   :  SPECID or zPHOT (for diangostic messages)
  //   ivar      :  variable index to fill in *MAP.IVARTYPE
  //  *VARNAME   :  name of map variable (table column) to examine
  //
  // Outputs:
  //   *MAP       :  search eff map struct to load.
  //
  // If assignment can't be made then abort.
  // Function returns mask of assigned variables.
  //
  // May 13 2021: allow HOSTMAG or HOST_MAG
  // Aug 27 2024: store IVAR_SALT2x1, IVAR_SALT2c, IVAR_LOGMASS
  // Mar 02 2026: store IVAR_SNRSUM_REST_V
  //

  bool USE_HOSTLIB = (HOSTLIB.NGAL_STORE > 0 ) ; 
  int  IVARTYPE_MASK = 0; 

  char fnam[] = "assign_MAP_VARNAME";  (void)fnam;

  // ----------- BEGIN ---------
  

  // - - - - -
  // check for the easy ones first.
  if ( strcmp(VARNAME,"SPECEFF") == 0 )  {
    MAP->IVAR           = ivar ;
    MAP->IVARTYPE[ivar] = IVARTYPE_EFFMAP ;
    return IVARTYPE_MASK ;
  } 

  else if ( strcmp(VARNAME,"HOSTEFF") == 0 )  {
    MAP->IVAR           = ivar ;
    MAP->IVARTYPE[ivar] = IVARTYPE_EFFMAP ; 
    return IVARTYPE_MASK;
  } 
  else if ( strcmp(VARNAME,"REDSHIFT") == 0 )  {
    MAP->IVAR_REDSHIFT  = ivar ;
    MAP->IVARTYPE[ivar] = IVARTYPE_EFFMAP_REDSHIFT ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_REDSHIFT );
    return IVARTYPE_MASK;
  }

  else if ( strcmp(VARNAME,"PEAKMJD") == 0 )  {
    MAP->IVAR_PEAKMJD = ivar ;
    MAP->IVARTYPE[ivar] = IVARTYPE_EFFMAP_PEAKMJD ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_PEAKMJD );
    return IVARTYPE_MASK;
  }
  else if ( strcmp(VARNAME,"DTPEAK") == 0 )  {
    MAP->IVAR_DTPEAK = ivar ;
    MAP->IVARTYPE[ivar] = IVARTYPE_EFFMAP_DTPEAK ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_DTPEAK );
    return IVARTYPE_MASK;
  }
  else if ( strcmp(VARNAME,"DTSEASON_PEAK") == 0 )  {
    MAP->IVAR_DTSEASON_PEAK = ivar ;
    MAP->IVARTYPE[ivar] = IVARTYPE_EFFMAP_DTSEASON_PEAK ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_DTSEASON_PEAK );
    return IVARTYPE_MASK ;
  }
  else if ( strcmp(VARNAME,"SALT2mB")==0 )  {
    MAP->IVAR_SALT2mB = ivar ;
    MAP->IVARTYPE[ivar] =  IVARTYPE_EFFMAP_SALT2mB ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_SALT2mB );
    return IVARTYPE_MASK ;
  }
  else if ( strcmp(VARNAME,"SALT2c")==0 )  {
    MAP->IVAR_SALT2c = ivar ;
    MAP->IVARTYPE[ivar] =  IVARTYPE_EFFMAP_SALT2c ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_SALT2c );
    return IVARTYPE_MASK ;
  }
  else if ( strcmp(VARNAME,"SALT2x1")==0 )  {
    MAP->IVAR_SALT2x1 = ivar ;
    MAP->IVARTYPE[ivar] =  IVARTYPE_EFFMAP_SALT2x1 ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_SALT2x1 );
    return IVARTYPE_MASK ;
  }
  else if ( strcmp(VARNAME,"LOGMASS") ==0 )  {
    MAP->IVAR_LOGMASS = ivar ;
    MAP->IVARTYPE[ivar] =  IVARTYPE_EFFMAP_LOGMASS ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_LOGMASS );
    return IVARTYPE_MASK ;
  }
  else if ( strcmp(VARNAME,"SNRSUM_REST_V") ==0 )  { // Mar 2026
    MAP->IVAR_SNRSUM_REST_V = ivar ;
    MAP->IVARTYPE[ivar] =  IVARTYPE_EFFMAP_SNRSUM_REST_V ;
    IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_SNRSUM_REST_V );
    return IVARTYPE_MASK ;
  }



  // - - - - - - - - - - - - - - -
  // check for filter-dependent mag or color.

  IVARTYPE_MASK = assign_MAP_VARNAME_FILTERS(MAPTYPE, ivar, VARNAME, MAP);
  if ( IVARTYPE_MASK > 0 ) { return IVARTYPE_MASK; }


  // check for other HOSTLIB variables that are not filter-dependent
  if ( USE_HOSTLIB ) {
    int NVAR_HOST, ivar_tmp;
    // check if VARNAME is in HOSTLIB 
    ivar_tmp = IVAR_HOSTLIB_STORE(VARNAME,0, fnam); 
    if ( ivar_tmp > 0 ) {
      //printf("\t Found %-20.20s in HOSTLIB for %s map \n", VARNAME, MAPTYPE);
      NVAR_HOST = MAP->NVAR_HOST;
      MAP->IVAR_HOST[NVAR_HOST]    = ivar ;     // SEARCHEFF map ivar  
      MAP->IVAR_HOSTLIB[ivar]      = ivar_tmp ; // hostlib ivar/column
      MAP->IVARTYPE[ivar]          = IVARTYPE_EFFMAP_HOSTLIB ; 
      IVARTYPE_MASK = ( 1 << IVARTYPE_EFFMAP_HOSTLIB );
      MAP->NVAR_HOST++ ;
      return IVARTYPE_MASK ;
    }
  }



  // ------------------------------------------------------
  // if we get here then abort.
  sprintf(c1err,"Unknown VARANME='%s' in SEARCHEFF_%s_FILE", VARNAME, MAPTYPE);
  sprintf(c2err,"Check manual for list of options.");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 

  return IVARTYPE_MASK ;


} // end of assign_MAP_VARNAME


int assign_MAP_VARNAME_FILTERS(char *MAPTYPE, int ivar, char *VARNAME, 
			       SEARCHEFF_MAP_DEF *MAP) {

  int IVARTYPE_MASK = 0 ;
  int IVARTYPE      = 0 ;

  bool HAS_PLUS  = strstr(VARNAME,PLUS)  != NULL ; // logical-OR among bands; eg l+p
  bool HAS_MINUS = strstr(VARNAME,MINUS) != NULL ; // color; e.g., g-r
  bool FOUND_MATCH_PREFIX = false ;
  bool IS_COLOR = HAS_MINUS ; 
  bool IS_MAG   = !HAS_MINUS ;

  int MXSUBSTR = MXSUBSTR_SEARCHEFF_MAP;

  char **ptrLIST, SUB[40], sep[2];
  float fmem;   (void)fmem;
  int  ipre, isub, NSUBSTR, NMATCH, ichar_band, ifiltlist[20], NFILT=0  ;
  int  i, ifilt_obs;

  // for PEAKMAG, allow color = PEAKMAG_[band0]-PEAKMAG_[band1], etc ...
#define NPREFIX_VARNAME 4
  char PREFIX_LIST[NPREFIX_VARNAME][20] = {  
    "PEAKMAG_",  // 0
    "HOSTMAG_",  // 1
    "_obs",      // 2
    "SBMAG_"     // 3
  } ;

  int IPREFIX_PEAKMAG = 0; 
  int IPREFIX_obs     = 2;

  int  IVARMAG_EFFMAP_LIST[NPREFIX_VARNAME] = { 
    IVARTYPE_EFFMAP_PEAKMAG, 
    IVARTYPE_EFFMAP_HOSTMAG, 
    IVARTYPE_EFFMAP_HOSTMAG, 
    IVARTYPE_EFFMAP_SBMAG 
  } ;

  int  IVARCOLOR_EFFMAP_LIST[NPREFIX_VARNAME] = { 
    IVARTYPE_EFFMAP_PEAKCOLOR,  
    IVARTYPE_EFFMAP_HOSTCOLOR, 
    IVARTYPE_EFFMAP_HOSTCOLOR, 
    IVARTYPE_EFFMAP_SBCOLOR 
  } ;

  char PREFIX[40], cband[2];
  int LDMP = 0 ;

  char fnam[] = "assign_MAP_VARNAME_FILTERS" ; (void)fnam;
  
  // -------------- BEGIN ---------------

  fmem = malloc_strlist(+1, MXSUBSTR, 40, &ptrLIST );

  for(i=0; i < 10; i++ ) { ifiltlist[i] = -9; }

  // check how many space-separated SUBSTR names to check
  // Split by both '+' and '-'
  nsplitString(VARNAME, "+-",  fnam, MXSUBSTR, &NSUBSTR, ptrLIST, sep );

  if ( LDMP ) {
    printf(" xxx ---------------------------------------------- \n");
    printf(" xxx %s DUMP for %s-MAP  ivar=%d  VARNAME=%s \n", 
	   fnam, MAPTYPE, ivar, VARNAME);
    printf(" xxx   HAS[+ -] = [ %d %d ]  IS_[MAG,COLOR] = [ %d %d ] \n", 
	   HAS_PLUS, HAS_MINUS, IS_MAG, IS_COLOR );
    fflush(stdout);
  }

  // - - - - 
  // check which prefix matches; return if none
  for ( ipre=0;  ipre < NPREFIX_VARNAME; ipre++ ) {
    NMATCH = 0 ;
    sprintf(PREFIX, "%s", PREFIX_LIST[ipre] );

    for ( isub=0; isub < NSUBSTR; isub++ ) {
      FOUND_MATCH_PREFIX = false;
      sprintf(SUB, "%s", ptrLIST[isub] );
      ichar_band = strlen(SUB) - 1; // default band location is last char of VARNAME

      if ( strstr(SUB,PREFIX) != NULL ) {  FOUND_MATCH_PREFIX = true ;  }

      // check for single band representation of PEAKMAG; e.g. r is same as PRAKMAG_r
      if ( ipre == IPREFIX_PEAKMAG && strlen(SUB)==1 ) 
	{  FOUND_MATCH_PREFIX = true;  PREFIX[0] = 0; }

      // check for HOSTLIB name [band]_obs ... treat same as HOSTMAG_[band]
      if ( ipre == IPREFIX_obs && FOUND_MATCH_PREFIX ) {
	if ( strlen(SUB) == 5 ) { ichar_band = 0; }  
	else                    { FOUND_MATCH_PREFIX = false; }  // fragile alert
      }

      if ( FOUND_MATCH_PREFIX ) { 
	sprintf(cband, "%c", SUB[ichar_band] );
	ifilt_obs         = INTFILTER(cband);
	ifiltlist[NMATCH] = ifilt_obs ;
	NMATCH++ ; 
	NFILT = NMATCH;
	if ( LDMP ) {
	  printf(" xxx   NMATCH=%d for SUBSTR='%s'  PREFIX='%s'  "
		 "ichar_band=%d  band='%s'  ifilt_obs=%d\n",
		 NMATCH, SUB, PREFIX, ichar_band, cband, ifilt_obs );
	  fflush(stdout);
	}
      }

    } // end isub loop


    if ( NMATCH > 0 && NMATCH != NSUBSTR ) {
      sprintf(c1err,"Found %d prefix matches among %d subtrings in VARNAME=%s",
	      NMATCH, NSUBSTR, VARNAME);
      sprintf(c2err,"Check %s in %s map", VARNAME, MAPTYPE);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    // if we have matches, set pointer to identified prefix
    if ( NMATCH == NSUBSTR ) { 
      if ( HAS_MINUS ) {
	IVARTYPE      = IVARCOLOR_EFFMAP_LIST[ipre] ;
      }
      else {
	IVARTYPE      = IVARMAG_EFFMAP_LIST[ipre] ;
      }	
    }

  } // end loop over possible prefixes

  // free memory
  fmem = malloc_strlist(+1, MXSUBSTR, 40, &ptrLIST );
  // - - - - - - - - - - -  -

  if ( IVARTYPE > 0 ) {
    IVARTYPE_MASK = ( 1 << IVARTYPE ) ;
    MAP->IVARTYPE[ivar] = IVARTYPE;

    if ( IS_MAG   ) { MAP->FLAG_MAG[ivar] =  FLAG_EFFMAP_MAG   ; }
    if ( IS_COLOR ) { MAP->FLAG_MAG[ivar] =  FLAG_EFFMAP_COLOR ; }

    if ( LDMP ) {
      printf(" xxx   IVARTYPE=%d  ivar=%d  IS_[MAG,COLOR] = [%d,%d]  FLAG_MAG=%d\n",
	     IVARTYPE, ivar, IS_MAG, IS_COLOR, MAP->FLAG_MAG[ivar] );
      printf(" xxx\t   (FLAG_EFFMAP_MAG=%d  FLAG_EFFMAP_COLOR=%d) \n",
	     FLAG_EFFMAP_MAG, FLAG_EFFMAP_COLOR);
      fflush(stdout);
    }
  }
  else {
    if ( LDMP ) { printf(" xxx   return with no matches. \n xxx \n"); fflush(stdout);  }
    return IVARTYPE_MASK ;
  }



  bool IS_PEAKMAG = ( IVARTYPE==IVARTYPE_EFFMAP_PEAKMAG || IVARTYPE==IVARTYPE_EFFMAP_PEAKCOLOR); 
  bool IS_HOSTMAG = ( IVARTYPE==IVARTYPE_EFFMAP_HOSTMAG || IVARTYPE==IVARTYPE_EFFMAP_HOSTCOLOR); 
  bool IS_SBMAG   = ( IVARTYPE==IVARTYPE_EFFMAP_SBMAG   || IVARTYPE==IVARTYPE_EFFMAP_SBCOLOR  ); 

  if ( LDMP ) { 
    printf(" xxx   IS[PEAKMAG, HOSTMAG, SBMAG] = [ %d %d %d ] \n",
	   IS_PEAKMAG, IS_HOSTMAG, IS_SBMAG ); fflush(stdout);
  }

  if ( IS_PEAKMAG ) {
    MAP->NFILTLIST_PEAKMAG[ivar]  = NFILT ;
    if ( LDMP ) { printf(" xxx   load NFILT=%d for ivar=%d  IFILT_OBS = %d %d\n", 
			 NFILT, ivar, ifiltlist[0], ifiltlist[1]); }
    for(i=0; i < NFILT; i++ ) { MAP->IFILTLIST_PEAKMAG[ivar][i] = ifiltlist[i] ; }
  }

  if ( IS_HOSTMAG ) {
    MAP->NFILTLIST_HOSTMAG[ivar]  = NFILT ;
    for(i=0; i < NFILT; i++ ) { MAP->IFILTLIST_HOSTMAG[ivar][i] = ifiltlist[i] ; }
  }

  if ( IS_SBMAG ) {
    MAP->NFILTLIST_SBMAG[ivar]  = NFILT ;
    for(i=0; i < NFILT; i++ ) { MAP->IFILTLIST_SBMAG[ivar][i] = ifiltlist[i] ; }
  }

  if ( LDMP ) { printf(" xxx \n"); fflush(stdout); }

  return IVARTYPE_MASK ; 

  
} // end assign_MAP_VARNAME_FILTERS





// ***************************************
void init_SEARCHEFF_zHOST(char *survey) {

  // Refactored March 2026
  //
  // Read optional efficiency map for obtaining a host-galaxy
  // redshift for UNCONFIRMED SNe.
  //

  char fnam[] = "init_SEARCHEFF_zHOST" ; (void)fnam;

  // --------------- BEGIN ------------
  
  init_searcheff_map(MAPTYPE_SEARCHEFF_zHOST, &SEARCHEFF_INFO_zHOST) ;

  read_searcheff_map(INPUTS_SEARCHEFF.USER_zHOST_FILE, &SEARCHEFF_INFO_zHOST) ;

  INPUTS_SEARCHEFF.NMAP_zHOST   = SEARCHEFF_INFO_zHOST.NMAP; // for snlc_sim

  return;
} // end of init_SEARCHEFF_zHOSTEFF


void read_searcheff_raw_varnames(char *SEARCHEFF_FILE, int *OPEN_STATUS, 
				 int *NVAR_RAW, char **VARNAMES_RAW) {

  // Created Mar 27 2026
  // For input EFFMAP_FILE, read VARNAMES and break up variables that have
  // + or -.
  //
  // Example: SPECEFF file has  VARNAMES:  r  g-i
  //   This functions returns *NVAR_RAW = 3  and VARNAMES_RAW = 'r', 'g', 'i'
  //   

  int MXVAR = MXVAR_SEARCHEFF_MAP ;
  int NSKIP = 1; // skip HOSTEFF of SPECEFF column
  int NKEY, NVAR, UNIQUE[MXVAR_SEARCHEFF_MAP], ivar, NVAR_RAW_LOCAL = 0 ;
  int MXSUBSTR=4, NSUBSTR, isub ;
  FILE *fp;
  float fmem; 
  char **VARNAMES_LOCAL, **SUBSTRING_LIST, SUB[40], sep[2], VARLIST[80] ;
  int LDMP = 0 ;
  char fnam[] = "read_searcheff_raw_varnames"; (void)fnam;

  // ------------ BEGIN -----------

  *NVAR_RAW    = 0 ;
  *OPEN_STATUS = 0;  // not SUCCESS, nor ERROR
  VARLIST[0]   = 0;
  if ( strlen(SEARCHEFF_FILE) < 6 ) { return; } // avoid options like ONE or ZERO or NONE

  int OPTMASK = 1;  // +=1(verbose)  +=4(don't check DOC)
  int gzipFlag;
  char PATH_LIST[] = "",   SEARCHEFF_FULLPATH[MXPATHLEN];
  fp  = snana_openTextFile (OPTMASK, PATH_LIST, SEARCHEFF_FILE, SEARCHEFF_FULLPATH, &gzipFlag);

  fmem += malloc_strlist(+1, MXVAR_SEARCHEFF_MAP, 40, &VARNAMES_LOCAL );
  fmem += malloc_strlist(+1, MXVAR_SEARCHEFF_MAP, 40, &SUBSTRING_LIST );

  read_VARNAMES_KEYS(fp, MXVAR, NSKIP, fnam, &NVAR, &NKEY, UNIQUE, VARNAMES_LOCAL);
  
  for(ivar=0; ivar < NVAR; ivar++ ) {

    // if VARNAME_LOCAL is g-r, split into 'g' and 'r' to report both RAW variables
    nsplitString(VARNAMES_LOCAL[ivar], "+-",  fnam, MXSUBSTR, &NSUBSTR, SUBSTRING_LIST, sep );
    for (isub=0; isub < NSUBSTR; isub++ ) {
      sprintf(SUB, "%s", SUBSTRING_LIST[isub] );
      sprintf(VARNAMES_RAW[NVAR_RAW_LOCAL], "%s", SUB );
      strcat(VARLIST," ");      strcat(VARLIST,SUB);
      if ( LDMP ) {
	printf(" xxx %s: VARNAMES_RAW[%d] = %s \n", fnam, NVAR_RAW_LOCAL, SUB); fflush(stdout);
      }
      NVAR_RAW_LOCAL++ ;
    } // end isub loop
  } // end ivar loop


  *NVAR_RAW = NVAR_RAW_LOCAL;

  printf("\t Found %d RAW variables: '%s' \n", NVAR_RAW_LOCAL, VARLIST);
  fflush(stdout);

  fmem = malloc_strlist(-1, MXVAR_SEARCHEFF_MAP, 40, &VARNAMES_LOCAL );
  fmem = malloc_strlist(-1, MXVAR_SEARCHEFF_MAP, 40, &SUBSTRING_LIST );
  fclose(fp);

  return;

} // end read_searcheff_raw_varnames


// ************************************
FILE *open_zHOST_FILE(int VBOSE) {

  // Mar 20 2019
  // Open zHOST file and return file pointer 
  //
  // VBOSE = +1 -> print to stdout
  // VBOSE =  0 -> do NOT print to stdout
  //
  // Apr 28 2026: just read; don't set or print flags 

  int  REQUIRE_zHOST_FILE=0, gzipFlag ;
  char *ptrFile_user ;
  char *ptrFile_final ;
  char localFile[MXPATHLEN];
  FILE *fp = NULL;
  char fnam[] = "open_zHOST_FILE" ; (void)fnam;

  // --------------- BEGIN ------------------

  fp = NULL ;

  /* xxxxxx mark delete 
  INPUTS_SEARCHEFF.IFLAG_zHOST_EFFZERO = 0;
  INPUTS_SEARCHEFF.IFLAG_zHOST_EFFONE  = 0;  
  INPUTS_SEARCHEFF.NMAP_zHOST = 0 ;
  xxxxxx end mark */

  ptrFile_user  = INPUTS_SEARCHEFF.USER_zHOST_FILE ; // from sim-input
  ptrFile_final = INPUTS_SEARCHEFF.zHOST_FILE ;      // includes full path
  
  if( IGNOREFILE(ptrFile_user) ) {
    // NULL or NONE, etc ... -> EFF=1
    return(fp);
  }
  else if ( strcmp(ptrFile_user,"ZERO") == 0 ) {

    /* xxxx mark delete 4.28.2026 xxxxxxx
    // no file to read, but set all EFF_zHOST=0  (Jun 2018)
    INPUTS_SEARCHEFF.IFLAG_zHOST_EFFZERO = 1;
    if ( LPRINT ) 
      { printf("\n  SEARCHEFF_zHOST_FILE=ZERO -> Force EFF_zHOST=0 \n"); }
    xxxxxxxxx end mark xxxx */

    return(fp) ;
  }
  else if ( strcmp(ptrFile_user,"ONE") == 0 ) {
    /* xxxx mark delete 4.28.2026 xxxxxxx
    // no file to read, but set all EFF_zHOST=1  (Nov 2024)
    INPUTS_SEARCHEFF.IFLAG_zHOST_EFFONE = 1;
    if ( LPRINT ) 
      { printf("\n  SEARCHEFF_zHOST_FILE=ONE -> Force EFF_zHOST=1 \n"); }
    xxxxxxxx end mark xxxxxx */
    return(fp) ;
  }
  else { 
    sprintf(localFile, "%s",  ptrFile_user); 
    REQUIRE_zHOST_FILE = 1 ;
  }


  // use utility to check local dir and path.
  int OPTMASK_OPEN  = INPUTS_SEARCHEFF.OPTMASK_OPENFILE ;
  if ( REQUIRE_zHOST_FILE ) {
    fp = snana_openTextFile(OPTMASK_OPEN, PATH_SEARCHEFF, localFile, 
			    ptrFile_final, &gzipFlag); // returned

    // examine if there is no zHOST file
    if ( !fp  ) {    
      abort_openTextFile("SEARCHEFF_zHOST_FILE",
			 PATH_SEARCHEFF, localFile, fnam );
    }    
  }
  

  // --------------------------------------------

  // always print message for opened file
  printf("   Open zHOST-for-unconfirmed effic map from \n\t %s\n", ptrFile_final );
  fflush(stdout) ;

  return(fp);

} // end open_zHOST_FILE




// ***************************************
int gen_SEARCHEFF ( int ID                  // (I) identifier 
		    ,double *EFF_SPECID     // (O) specID eff
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

   Beware that return arg *EFF_SPECID is really just the 
   spectroscopic efficiency because there is currently
   no way to calculate the pipeline efficiency.


  Mar 2018: 
   + use APPLYMASK_SEARCHEFF_xxx parmaeters to set return MASK
   + add argument *EFF_zHOST

  Mar 2026: set APPLYMASK_SEARCHEFF_zSPEC

  *****/

  int  LFIND1_PIPELINE, LFIND2_SPECID, LFIND3_zHOST, MASK ;
  char fnam[]  = "gen_SEARCHEFF" ; (void)fnam;

  // ----------------- BEGIN -------------


  // init function args
  *EFF_SPECID = *EFF_zHOST = 0.0;  
  MJD_DETECT->TRIGGER = 1.0E6 ;
  MJD_DETECT->FIRST   = 1.0E6 ;
  MJD_DETECT->LAST    = 1.0E6 ;
  MASK = 0 ;  // init return arg


  // 5/26/2009: minimum trigger is at least 2 measurements
  if ( SEARCHEFF_DATA.NOBS < INPUTS_SEARCHEFF.MINOBS ) { return MASK ; }

  LFIND1_PIPELINE = LFIND2_SPECID = LFIND3_zHOST = 0 ;

  // --------------------------
  // if nothing is defined for software efficieincy,
  // just set it to FOUND and skip to SPEC eff.

  if ( INPUTS_SEARCHEFF.NMAP_DETECT == 0 ) {
    LFIND1_PIPELINE = 1;      // set flag that search software finds SN
  }
  else {
    // check if software/pipeline detection finds this SN
    LFIND1_PIPELINE = gen_SEARCHEFF_PIPELINE(ID,MJD_DETECT);
  }


  // ------------------------------------
  // if detected by pipeline, check if spec-confirmed
  if ( LFIND1_PIPELINE )  { 
    LFIND2_SPECID = gen_SEARCHEFF_SPECID(ID, EFF_SPECID) ; 
  }


  // if not spec-confirmed, check for spec zHOST
  if ( LFIND1_PIPELINE  &&  LFIND2_SPECID == 0 )  { 
     LFIND3_zHOST = gen_SEARCHEFF_zHOST(ID,EFF_zHOST) ;   // return EFF_zHOST
  }

  // --- set return bit-MASK ---- 

  if ( LFIND1_PIPELINE  ) {
    MASK += APPLYMASK_SEARCHEFF_PIPELINE ;  // software trigger 
  }

  if ( LFIND2_SPECID    ) {
    MASK += APPLYMASK_SEARCHEFF_SPECID  ;  // got spec-confirmation 
    MASK |= APPLYMASK_SEARCHEFF_zSPEC   ;  // Mar 7 2026: got zSPEC
  }
  if ( LFIND3_zHOST     ) {
    MASK += APPLYMASK_SEARCHEFF_zHOST ;  // got host-gal redshift
    MASK |= APPLYMASK_SEARCHEFF_zSPEC   ;  // Mar 7 2026: got zSPEC
  }

  return(MASK) ;

} // end of gen_SEARCHEFF



// ******************************************************
int gen_SEARCHEFF_SPECID(int ID, double *EFF_SPECID) {

  int  IS_SPECID ;
  char fnam[] = "gen_SEARCHEFF_SPECID";  (void)fnam;

  // ---------- BEGIN -------------

  // return EFF_SPECID
  IS_SPECID = gen_searcheff_map(ID, &SEARCHEFF_INFO_SPECID, EFF_SPECID) ; 


  return IS_SPECID ;

} // end gen_SEARCHEFF_SPECID

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
  int FOUND_TRIGGER=0, LCUT_PHOTPROB ;
  int OBSMARKER_DETECT[MXOBS_TRIGGER];
  double  RAN, EFF, MJD, MJD_LAST, MJD_DIF, TDIF_NEXT, SNR,MAG;
  double  PHOTPROB, CUTVAL, SEP_NEAREST_SRC, PSFSIG, MINSEP_DETECT ;
  char CFILT[4];
  int LDMP  = (ID == -39 ); 
  char fnam[] = "gen_SEARCHEFF_PIPELINE"; (void)fnam;

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
  SEP_NEAREST_SRC  = SEARCHEFF_DATA.SEP_NEAREST_SRC;

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

    // Jul 2022 check resolving nearby source 
    if ( SEP_NEAREST_SRC < 5.0 ) {
      PSFSIG = SEARCHEFF_DATA.PSFSIG[obs]; // PSF sigma, arcsec
      // if less than 1 sigma separation, disable detection
      MINSEP_DETECT = INPUTS_SEARCHEFF.NPSFSIGMA_MINSEP_DETECT * PSFSIG ;
      if ( SEP_NEAREST_SRC < MINSEP_DETECT ) { DETECT_FLAG = 0; }
    }

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

    SNR      = SEARCHEFF_DATA.SNR_CALC[obs] ;
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
  char fnam[] = "dumpLine_PIPELINE_PHOTPROB";  (void)fnam;
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
  // Jun 15 2022: check opt for single-exposure detections instead of coadd
  // Nov 30 2022: check for FIELD dependence

  int NMAP                = INPUTS_SEARCHEFF.NMAP_DETECT ;
  int APPLY_DETECT_SINGLE = INPUTS_SEARCHEFF.APPLY_DETECT_SINGLE ;

  double MAG, SNR, MJD, EFF, FIX_EFF, XNEXPOSE;
  double EFF_atmax, EFF_atmin, VAL_atmax, VAL_atmin, VAL ;
  double ZERO = 0.0, ONE  = 1.0 ;

  int CID, ifilt_obs, NPE_SAT, NBIN_EFF, imap, IMAP, NMAP_FOUND=0;
  int OPT_INTERP  = 1;   // 1=linear;  2=quadratic
  bool MATCH_FILTER, MATCH_FIELD;

  char cfilt[4], *field_map, *filt_map;
  char fnam[] ="GETEFF_PIPELINE_DETECT" ;  (void)fnam;

  // ---------- BEGIN ---------

  // check debugging option with fixed effic.
  FIX_EFF = INPUTS_SEARCHEFF.FIX_EFF_PIPELINE ;
  if ( FIX_EFF > 0.0 ) { return FIX_EFF ; }

  EFF       = 0.0 ;

  // find map corresponding to filter and [optional] FIELD
  ifilt_obs = SEARCHEFF_DATA.IFILTOBS[obs] ;  IMAP=-9;
  sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] );

  for(imap=0; imap < NMAP; imap++ ) {
    field_map     = SEARCHEFF_DETECT[imap].FIELDLIST;
    filt_map      = SEARCHEFF_DETECT[imap].FILTERLIST ;
    MATCH_FILTER  = ( strstr(filt_map,cfilt) != NULL );
    if ( strlen(field_map) > 0 ) 
      { MATCH_FIELD   = MATCH_SEARCHEFF_FIELD(field_map); }
    else
      { MATCH_FIELD = true; }

    if ( MATCH_FILTER && MATCH_FIELD ) 	
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
    sprintf(c2err,"Check EFF maps in %.*s", 
	    MXCHAR_MSGERR, INPUTS_SEARCHEFF.PIPELINE_EFF_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  CID       = SEARCHEFF_DATA.CID ;
  NBIN_EFF  = SEARCHEFF_DETECT[IMAP].NBIN;
  SNR       = SEARCHEFF_DATA.SNR_CALC[obs] ;
  MAG       = SEARCHEFF_DATA.MAG[obs] ;
  MJD       = SEARCHEFF_DATA.MJD[obs] ;
  NPE_SAT   = SEARCHEFF_DATA.NPE_SAT[obs] ; // Npe above saturation
  XNEXPOSE  = (double)SEARCHEFF_DATA.NEXPOSE[obs] ; 

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
  // Check option to trigger on single-exposures instead of coadd.
  // Note assumption that each exposure has same SNR.
  if ( APPLY_DETECT_SINGLE ) { SNR  /= sqrt(XNEXPOSE); }

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

  // - - - - - 
  if ( APPLY_DETECT_SINGLE && XNEXPOSE > 1.0 ) {
    // EFF -> 1 - (1-EFF_1)^NEXPOSE
    // E.g., NEXPOSE=10 and EFF_1 = 0.02 -> EFF = 1 - 0.98^10 = 0.183
    // E.g., NEXPOSE=10 and EFF_1 = 0.05 -> EFF = 1 - 0.95^10 = 0.401
    // E.g., NEXPOSE=10 and EFF_1 = 0.50 -> EFF = 1 - 0.50^10 = 0.999
    double EFF_1 = EFF, INEFF_1 = 1.0-EFF_1 ;
    EFF = 1.0 - pow(INEFF_1,XNEXPOSE);
  }

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
  char fnam[]      = "setObs_for_PHOTPROB" ;  (void)fnam;

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

    if ( strcmp(FILT_TMP,ALL) ==0 )        { MATCH_FILT=1; }
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
  char   fnam[] = "setRan_for_PHOTPROB" ;  (void)fnam;


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
  double SNR       = SEARCHEFF_DATA.SNR_CALC[obs] ;
  double SBMAG     = SEARCHEFF_DATA.SBMAG[IFILTOBS] ;  
  //  double GALMAG    = SEARCHEFF_DATA.HOSTMAG[IFILTOBS] ;  
  double MJD       = SEARCHEFF_DATA.MJD[obs] ;

  int    istat, ivar, LDMP ;
  double PHOTPROB_CDF[MXVAR_SEARCHEFF_PHOTPROB];
  double VARDATA[MXVAR_SEARCHEFF_PHOTPROB];
  char  *VARNAME, cFILT[4];
  double PHOTPROB = 0.0 ;

  char fnam[] = "get_PIPELINE_PHOTPROB" ;  (void)fnam;

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


  LDMP = 0 ; // ( fabs(SEARCHEFF_DATA.SNR_CALC[obs]-8.0) < 0.1 );
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
int gen_searcheff_map(int ID, SEARCHEFF_INFO_DEF *SEARCHEFF_INFO, double *EFF) {

  //
  // Return 1 if this SN is selected by MAPs defined by SEARCHEFF_INFO; 
  // return 0 otherwise.
  // General utility for SPECID and zHOST efficiency maps.
  //
  //  [Mar 2026 Refactor from original/obsolete gen_SEARCHEFF_SPECID]
  //
  // Inputs:
  //   ID:        event ID (for diagnostic/error messages) 
  //   MAPTYPE:   SPECID or zHOST (for diagnostic messages)
  //   SEARCHEFF_INFO : struct containing map info
  //
  // Outputs
  //    *EFF efficiency from map
  //

  char *MAPTYPE       = SEARCHEFF_INFO->MAPTYPE ;
  int  NMAP           = SEARCHEFF_INFO->NMAP ;
  int  BOOLEAN_OR     = SEARCHEFF_INFO->BOOLEAN_OR  ;
  int  BOOLEAN_AND    = SEARCHEFF_INFO->BOOLEAN_AND ;
  int  IFLAG_EFFZERO  = SEARCHEFF_INFO->IFLAG_EFFZERO ;
  int  IFLAG_EFFONE   = SEARCHEFF_INFO->IFLAG_EFFONE ;
  int  OPT_FIELDMATCH_REQUIRE = SEARCHEFF_INFO->OPT_FIELDMATCH_REQUIRE;

  int    imap, ivar, NVAR, LFIND, istat, NMATCH = 0, INDX_RAN ;
  bool   MATCH_FIELD, MATCH_PEAKMJD, REQUIRE_PASS = true ;
  double PROB_all_FAIL, PROB_all_PASS;  
  double EFF_LOCAL, RAN, VARDATA[MXVAR_SEARCHEFF_MAP], PEAKMJD=-9.0 ;
  double *ptr_FLAT_RANDOMS;
  char   *FIELD_MAP, *VARNAME ;
  int    LDMP = 0 ; // (ID == 7) ;
  char   fnam[] = "gen_searcheff_map" ;  (void)fnam;

  // ----------- BEGIN --------

  LFIND      = 1 ;
  EFF_LOCAL  = 1.0 ;

  if ( INPUTS_SEARCHEFF.FUNEFF_DEBUG ) {
    RAN   = SEARCHEFF_RANDOMS.FLAT_SPEC[0] ; 
    return gen_SEARCHEFF_DEBUG(MAPTYPE, RAN, EFF) ;
  }

  if ( strcmp(MAPTYPE,MAPTYPE_SEARCHEFF_SPECID) == 0 ) {
    ptr_FLAT_RANDOMS = SEARCHEFF_RANDOMS.FLAT_SPEC;
    INDX_RAN = 0 ;
  }
  else if ( strcmp(MAPTYPE,MAPTYPE_SEARCHEFF_zHOST) == 0 ) {
    ptr_FLAT_RANDOMS = SEARCHEFF_RANDOMS.FLAT_zHOST;
    INDX_RAN = 50 ;
  }
  else {
    return 0;
  }

  // check option of there is no SEARCHEFF map file
  if ( NMAP == 0 )  { 
    if ( IFLAG_EFFZERO ) 
      { EFF_LOCAL = 0.0 ; LFIND=0 ; }
    else if ( IFLAG_EFFONE ) 
      { EFF_LOCAL = 1.0 ; LFIND=1 ; }
    else
      { EFF_LOCAL = 1.0 ; LFIND=1 ; }

    *EFF = EFF_LOCAL;
    return LFIND ;  
  }

  // do multi-dimensional interpolation if SPEC-EFF mapis specified. 
  // Take the logical-OR of each map,
  // which means that EFF = 1 - product(1-Eff_imap)

  if ( LDMP ) {
    printf(" xxx ================================================================= \n");
    printf(" xxx %s: %s MAP DUMP for CID = %d  \n", 
	   fnam, MAPTYPE, ID);
    printf(" xxx %s: BOOLEAN_[OR,AND] = %d, %d \n",
	   fnam, BOOLEAN_OR, BOOLEAN_AND );
    fflush(stdout);
  }
  
  // - - - - - - - - - - - - -
  PROB_all_FAIL = PROB_all_PASS = 1.0;  
  SEARCHEFF_MAP_DEF *MAP;

  for ( imap=0; imap < NMAP; imap++ ) {

    MAP = &SEARCHEFF_INFO->MAP_LIST[imap];

    // check if current FIELD & PEAKMJD goes with this map
    FIELD_MAP   = MAP->FIELDLIST ;
    MATCH_FIELD = MATCH_SEARCHEFF_FIELD(FIELD_MAP);

    PEAKMJD = SEARCHEFF_DATA.PEAKMJD ;
    MATCH_PEAKMJD = ( PEAKMJD >= MAP->PEAKMJD_RANGE[0] && 
		      PEAKMJD <= MAP->PEAKMJD_RANGE[1] );

    if ( LDMP ) {
      printf(" xxx - - - - - - - - - - - - - - - - - - \n");
      printf(" xxx %s: FIELD=%s  PEAKMJD=%.0f  MATCH[FIELD,PEAKMJD] = %d, %d \n",
	     fnam, SEARCHEFF_DATA.FIELDNAME, PEAKMJD, MATCH_FIELD, MATCH_PEAKMJD ); fflush(stdout);
    }
    if ( !(MATCH_PEAKMJD && MATCH_FIELD)  ) { continue ; }

    NMATCH++ ; // count how many maps are matched to field/peakmjd
    
    // determine list of variables
    NVAR = MAP->GRIDMAP.NDIM ;
    for ( ivar=0; ivar < NVAR; ivar++ )  {
      VARDATA[ivar] = LOAD_SEARCHEFF_VAR(MAPTYPE, MAP, ivar); 
      VARDATA[ivar] += MAP->SHIFT[ivar]; // optional systematic shift
      if ( LDMP ) {
	VARNAME = MAP->VARNAMES[ivar] ;	
	printf(" xxx %s MATCH: imap=%d \n", fnam, imap);
	printf(" xxx %s \t load ivar=%d:  %s = %f \n",
	       fnam, ivar, VARNAME, VARDATA[ivar]); fflush(stdout);
	printf(" xxx %s \t FIELD=%s  PEAKMJD_RANGE=%.0f-%.0f \n",
	       fnam, MAP->FIELDLIST, MAP->PEAKMJD_RANGE[0], MAP->PEAKMJD_RANGE[1]);
	fflush(stdout);
      }
    }

    istat = interp_GRIDMAP(&MAP->GRIDMAP, VARDATA, &EFF_LOCAL );  (void)istat;

    if ( LDMP ) {
      printf(" xxx %s \t EFF(interp)=%.3f \n", fnam, EFF_LOCAL );
      fflush(stdout);
    }
    
    PROB_all_FAIL *= (1.0 - EFF_LOCAL); // prob that all maps fail
    PROB_all_PASS *= EFF_LOCAL ;        // prob that all maps pass
    
    // Aug 2024: check for required map; e.g.. PEAKMJD or DTPEAK range
    if ( MAP->REQUIRE ) {
      RAN = ptr_FLAT_RANDOMS[imap+10] ;
      if ( RAN > EFF_LOCAL ) { REQUIRE_PASS = false ; } 
    }
    
  } // imap 

  // - - - - - -
  if ( NMATCH == 0 ) {
    if ( OPT_FIELDMATCH_REQUIRE > 0 ) {
      sprintf(c1err, "Invalid NMATCH=%d  %s maps for", NMATCH, MAPTYPE);
      sprintf(c2err, "FIELD = '%s'  PEAKMJD=%.3f", SEARCHEFF_DATA.FIELDNAME, PEAKMJD );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
    else {
      return 0.0 ;  
    }
  }

  // - - - - -

  if ( NMAP > 0  ) { 
    if ( BOOLEAN_OR  )
      { EFF_LOCAL = 1.0 - PROB_all_FAIL ; }
    else if ( BOOLEAN_AND )
      { EFF_LOCAL = PROB_all_PASS  ; }
    else {
      sprintf(c1err,"Undefined BOOL logic for SPECEFF maps.");
      sprintf(c2err,"Expecting OR or AND logic.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;       
    }
    
  }

  // allow round-off
  if ( EFF_LOCAL < 0.0 && EFF_LOCAL > -0.001 ) { EFF_LOCAL = 0.0 ; } 
  if ( EFF_LOCAL > 1.0 && EFF_LOCAL <  1.0003) { EFF_LOCAL = 0.99999; }
  // sanity check
  if ( EFF_LOCAL < 0.0 || EFF_LOCAL > 1.0003 ) {
    sprintf(c1err,"Invalid EFF(%s) = %.4f for ID=%d", MAPTYPE, EFF_LOCAL, ID);
    sprintf(c2err,"zHel=%6.3f  PEAKMJD=%9.2f  FIELD=%s", 
	    SEARCHEFF_DATA.REDSHIFT, 
	    SEARCHEFF_DATA.PEAKMJD, SEARCHEFF_DATA.FIELDNAME );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  // apply user scale (May 10 2018)
  // ??? EFF_LOCAL *= INPUTS_SEARCHEFF.USER_SPECEFF_SCALE ;

  // load output argument
  *EFF = EFF_LOCAL ;  

  // compare spec-eff to random number
  // xxx mark  RAN = SEARCHEFF_RANDOMS.FLAT_SPEC[0] ;
  RAN = ptr_FLAT_RANDOMS[INDX_RAN];
  if ( RAN > EFF_LOCAL ) { LFIND = 0 ; }

  if ( LDMP ) {
    printf(" xxx %s: PROB_all_FAIL=%.4f  PROB_all_PASS=%.4f \n",
	   fnam, PROB_all_FAIL, PROB_all_PASS);
    printf(" xxx %s: RAN=%.3f  EFF(FINAL) = %.3f  LFIND=%d \n",
	   fnam, RAN, EFF_LOCAL, LFIND); 
    fflush(stdout);
  }
  // Aug 2024: if required map is not satisfied, reject this event
  if ( !REQUIRE_PASS ) { LFIND = 0 ; }
  
  return LFIND ;

} // end of gen_searcheff_map

// ************************************************
double LOAD_SEARCHEFF_VAR(char *MAPTYPE, SEARCHEFF_MAP_DEF *MAP, int ivar) {

  // Refactored Mar 2026 to work for both SPECID and zHOST maps,
  // and to include optional HOSTLIB columns.
  // Function returns value of variable corresponding to this MAP & ivar index

  int    IVARTYPE ;
  double mag; 
  char *varName = MAP->VARNAMES[ivar] ;
  char fnam[] = "LOAD_SEARCHEFF_VAR" ;  (void)fnam;

  // ------------ BEGIN ------------------

      
  IVARTYPE = MAP->IVARTYPE[ivar]; 
 

  if ( IVARTYPE == IVARTYPE_EFFMAP_HOSTLIB ) {        // March 2026
    int IGAL         = SNHOSTGAL.IGAL ;
    int ivar_HOSTLIB = MAP->IVAR_HOSTLIB[ivar];
    double val       = HOSTLIB.VALUE_ZSORTED[ivar_HOSTLIB][IGAL] ;
    //    printf(" xxx %s: ivar=%d  ivar_HOSTLIB=%d(%s)  value=%f \n",
    //	   fnam, ivar, ivar_HOSTLIB, MAP->VARNAMES[ivar], val) ; fflush(stdout);
    return val;
  }

  else if ( IVARTYPE == IVARTYPE_EFFMAP_REDSHIFT ) {
    return  SEARCHEFF_DATA.REDSHIFT ; 
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_PEAKMJD ) {
    return  SEARCHEFF_DATA.PEAKMJD ; 
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_DTPEAK ) {
    return  SEARCHEFF_DATA.DTPEAK_MIN ; // smallest abs(T-Tpeak)
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_DTSEASON_PEAK ) {
    return  SEARCHEFF_DATA.DTSEASON_PEAK ; // smallest abs(T-Tpeak)
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_LOGMASS ) {
    return  SEARCHEFF_DATA.LOGMASS ; 
  }  
  else if ( IVARTYPE == IVARTYPE_EFFMAP_SNRSUM_REST_V ) {
    return  SEARCHEFF_DATA.SNRSUM_REST_V ; 
  }  
  else if ( IVARTYPE == IVARTYPE_EFFMAP_SALT2mB ) {
    return  SEARCHEFF_DATA.SALT2mB ; 
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_SALT2x1 ) {
    return  SEARCHEFF_DATA.SALT2x1 ; 
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_SALT2c ) {
    return  SEARCHEFF_DATA.SALT2c ; 
  }

  else if ( IVARTYPE == IVARTYPE_EFFMAP_PEAKMAG || IVARTYPE == IVARTYPE_EFFMAP_PEAKCOLOR) {
    mag = get_searcheff_mag(MAPTYPE, MAP->FLAG_MAG[ivar], MAP->MAGSHIFT,
			    MAP->VARNAMES[ivar],
			    MAP->NFILTLIST_PEAKMAG[ivar], MAP->IFILTLIST_PEAKMAG[ivar], 
			    SEARCHEFF_DATA.PEAKMAG );
    return mag;
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_HOSTMAG || IVARTYPE == IVARTYPE_EFFMAP_HOSTCOLOR) {
    mag = get_searcheff_mag(MAPTYPE, MAP->FLAG_MAG[ivar], MAP->MAGSHIFT,
			    MAP->VARNAMES[ivar],
			    MAP->NFILTLIST_HOSTMAG[ivar], MAP->IFILTLIST_HOSTMAG[ivar], 
			    SEARCHEFF_DATA.HOSTMAG );
    return mag;
  }
  else if ( IVARTYPE == IVARTYPE_EFFMAP_SBMAG || IVARTYPE == IVARTYPE_EFFMAP_SBCOLOR) {
    double MAGSHIFT = 0.0 ; // no shift for SBmag
    mag = get_searcheff_mag(MAPTYPE, MAP->FLAG_MAG[ivar], MAGSHIFT,
			    MAP->VARNAMES[ivar],
			    MAP->NFILTLIST_SBMAG[ivar], MAP->IFILTLIST_SBMAG[ivar], 
			    SEARCHEFF_DATA.SBMAG );
    return mag;
  }


  // if we get here then abort.
  sprintf(c1err,"Could not find ivar=%d for varName = %s ", ivar, varName);
  sprintf(c2err,"Check %s map", MAPTYPE );
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 


  return(-9.0); // should never get here

} // end of LOAD_SEARCHEFF_VAR


double get_searcheff_mag(char *MAPTYPE, int FLAG_MAG, double MAGSHIFT, char *VARNAME,
			 int NFILT, int *IFILTLIST, double *MAG_DATA) {

  // Created Mar 27 2026
  // utility to return mag among list of bands, or color from list of bands.
  // Inputs:
  //   MAPTYPE    : SPECID or zHOST
  //   FLAG_MAG   : FLAG_EFFMAP_MAG (for mag) or FLAG_EFFMAP_COLOR
  //   MAGSHIFT   : for systematics
  //   VARNAME    : var name in map (for diagnostic only)
  //   NFILT      : number of filters passed
  //   IFILTLIST  : list of NFILT absolute filter indices
  //   MAG_DATA   : list of NFILT mag values from sim data
  //
  // Function returns either mag or color based on FLAG_MAG value.

  int ifilt, ifilt_obs  ;
  int LDMP = 0 ;
  char fnam[] = "get_searcheff_mag" ;  (void)fnam;

  // ------------ BEGIN --------------

  if ( NFILT <= 0 || NFILT > 20 )  {
    sprintf(c1err,"Invalid NFILT=%d for VARNAME=%s in %s-MAP", NFILT, VARNAME, MAPTYPE);
    sprintf(c2err,"FLAG_MAG = %d ", FLAG_MAG);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  if ( LDMP ) {
    printf(" xxx --------------------------------------------------- \n");
    printf(" xxx %s(%s) DUMP for %s-MAP FLAG_MAG=%d NFILT=%d \n", 
	   fnam, VARNAME, MAPTYPE, FLAG_MAG, NFILT );
    fflush(stdout);
  }

  if ( FLAG_MAG == FLAG_EFFMAP_MAG ) {

    double mag=0.0,  mag_tmp=0.0;   int nmag_valid=0 ;
    for(ifilt=0; ifilt < NFILT; ifilt++ ) {
      ifilt_obs = IFILTLIST[ifilt] ;
      mag_tmp   = MAG_DATA[ifilt_obs] ;
      if ( mag_tmp < 40. && mag_tmp > 0.0 ) { mag += mag_tmp; nmag_valid++; }

      if ( LDMP ) {
	printf(" xxx %s(%s): ifilt_obs = %d mag_tmp = %.3f \n", 
	       fnam, VARNAME, ifilt_obs, mag_tmp);	fflush(stdout);
      }
    }
    if ( nmag_valid > 0 ) { 
      mag /= (double)nmag_valid ;
      mag += MAGSHIFT; // systematic shift; see user keys MAGSHIFT_[SPECEFF,zHOSTEFF]
    }
    else
      { mag = -9.0; }
    

    if ( LDMP )
      { printf(" xxx %s(%s): return  mag = %.3f \n",  fnam, VARNAME, mag); fflush(stdout); }

    return mag;
  }
  else if ( FLAG_MAG == FLAG_EFFMAP_COLOR ) {

    double mag0, mag1, color ;   int ifilt_obs0, ifilt_obs1;
    ifilt_obs0 = IFILTLIST[0] ;    mag0 = MAG_DATA[ifilt_obs0];
    ifilt_obs1 = IFILTLIST[1] ;    mag1 = MAG_DATA[ifilt_obs1];
    check_magUndefined(mag0, VARNAME, fnam); 
    check_magUndefined(mag1, VARNAME, fnam);  
	  
    color = mag0 - mag1 ;

    if ( LDMP ) {
      printf(" xxx %s(%s): ifilt_obs= %d,%d   color = %.3f - %.3f = %.3f \n",
	     fnam, VARNAME, ifilt_obs0, ifilt_obs1, mag0, mag1, color); fflush(stdout);
    }

    return color ;
  }

  else {
    sprintf(c1err,"Invalid FLAG_MAG = %d",  FLAG_MAG);
    sprintf(c2err,"Expecting FLAG_MAG = %d(mag) or %d(color)",
	    FLAG_EFFMAP_MAG, FLAG_EFFMAP_COLOR);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  // should never get here
  return MAG_UNDEFINED;

} // end get_searcheff_mag

// ************************************************
int gen_SEARCHEFF_zHOST(int ID, double *EFF_zHOST) {

  // Created May 8 2014:
  // Check efficiency of finding host galaxy redshift;
  // return 0, 1 if not found, found.
  //
  // Jan 22 2017: check for field-dependent HOSTEFF map.
  // Mar 11 2018: refactor
  // Oct 15 2025: 
  //    + bail if hostless event (can't get zSpec_host if host is not found)

  bool  USE_HOSTLIB   = (HOSTLIB.NGAL_STORE > 0 ) ; // Oct 2025
  int     LFIND ;
  char fnam[] = "gen_SEARCHEFF_zHOST" ;  (void)fnam;

  // ----------- BEGIN -----------


  *EFF_zHOST = 0.0 ;
  if ( INPUTS_SEARCHEFF.RESTORE_DES5YR ) {
    // restore D5YR bug by ignoring hostless events
  }
  else {
    // implement host-less fix Oct 15 2025 
    // (original motivation: for DESC all photo-z/re-analysis of DES-SN5YR)
    if ( USE_HOSTLIB && SNHOSTGAL.NNBR_DDLRCUT == 0 ) { return 0; }
  }


  LFIND = gen_searcheff_map(ID, &SEARCHEFF_INFO_zHOST, EFF_zHOST) ; 

  return LFIND ;

}   // end of gen_SEARCHEFF_zHOST



// ===============================================
bool MATCH_SEARCHEFF_FIELD(char *field_map) {

  // Created Feb 5 2021
  // Return true iff SEARCHEFF_DATA.FIELDLIST_OVP matches field_map.
  // 
  // Ideally all overlap fields are checked, but this causes
  // downstream abort if an event overlaps two FIELDS.
  // Here we only check first field loaded in FIELDLIST_OVP;
  // maybe somebody will have an algorithm for choosing
  // among multiple FIELD-dependent maps.
  //
  // Previous logic had checked last overlap field (since each
  // subsequent field had clobbered previous field), and here we
  // check first field ... so new logic (Feb 5 2021) can result
  // in using a different map for field overlaps.
  //

  int  i;
  bool OVP;
  char *field_data;
  char fnam[] = "MATCH_SEARCHEFF_FIELD";  (void)fnam;

  // ---------- BEGIN ----------

  if ( strcmp(field_map,ALL)  == 0    ) { return true ; }

  // xxx  for(i=0; i < NFIELD_OVP; i++ ) { // maybe someday ??

  for(i=0; i < 1; i++ ) { // only check first FIELD among overlaps
    field_data = SEARCHEFF_DATA.FIELDLIST_OVP[i];

    // xxx mark delete Mar 3 2026 if ( strstr(field_map,field_data) != NULL ) { return true ; }
    OVP = ( strstr(field_map,field_data) != NULL );

    if ( OVP ) { return true ; }
  }

  // if we get here, there is no match -> return false
  return false ;

} // end MATCH_SEARCHEFF_FIELD


// *********************************************
void LOAD_PHOTPROB_CDF(int NVAR_CDF, double *VAL ) {

  // Created Mar 2018
  // for input WGT values (*VAL), replace *VAL with CDF.
  // Note hat *VAL is input and output.

  int ivar;
  double SUM=0.0, sum=0.0, TMPVAL[MXVAR_SEARCHEFF_PHOTPROB];
  char fnam[] = "LOAD_PHOTPROB_CDF" ;  (void)fnam;

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
  double SNR      = SEARCHEFF_DATA.SNR_CALC[OBS]; 
  double SBMAG    = SEARCHEFF_DATA.SBMAG[IFILTOBS]; 
  double GALMAG   = SEARCHEFF_DATA.HOSTMAG[IFILTOBS]; 
  double REDSHIFT = SEARCHEFF_DATA.REDSHIFT ; 
  int    IVARABS ;
  double VAL=0.0 ;
  char *VARNAME;
  char fnam[] = "LOAD_PHOTPROB_VAR" ;  (void)fnam;

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




// ============================================================
int  gen_SEARCHEFF_DEBUG(char *WHAT, double RAN, double *EFF) {

  // Created Jan 30 2015
  // DEBUG/hack efficiency function for specialized tests
  // (never for a real analysis).
  //
  // Input *WHAT = "PIPELINE", "SPEC", or "zHOST".
  // Input RAN is a random number from 0 to 1
  // Output *EFF is the efficiency

  int    IDEBUG ;
  double RAN_LOCAL, EFF_LOCAL, zTau ;
  char fnam[] = "gen_SEARCHEFF_DEBUG" ;  (void)fnam;
  // --------- BEGIN ------------

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



// ======= END: =========
