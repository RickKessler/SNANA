/***********************************************************
  sntools_hpak_legacy.cpp

  Legacy SNLCPAK (light curve pack) and SPECPAK (spectra pack)
  functions.

  Nov 09 2026: REFACTOR ONLY -- moved verbatim from sntools_output.cpp
               in preparation for removing SNLCPAK/SPECPAK.
               No functional changes.

  This is NOT a stand-alone compilation unit: it is an include file
  that is #included from sntools_output.cpp at the location where
  these functions used to be defined, so that it sees the same
  macros, globals and package-specific functions (e.g., 
  SNLCPAK_FILL_TEXT, SPECPAK_FILL_MARZ) as before.
  Do not add it to the list of objects to compile.
***********************************************************/

// ==========================================
//    SNLCPAK functions
// ==========================================

void HPAK_INIT(void) {
  // Created Nov 9 2026
  // one-time global init 

  SNLCPAK_USE_ROOT = SNLCPAK_USE_TEXT = SNLCPAK_USE_NPZ = 0 ;
  SPECPAK_USE_ROOT = SPECPAK_USE_TEXT = SPECPAK_USE_NPZ = 0 ;
  SPECPAK_USE_MARZ  = 0;

  sprintf(SNLCPAK_OUTPUT.SURVEY,             "NULL" );
  sprintf(SNLCPAK_OUTPUT.VERSION_PHOTOMETRY, "NULL" );
  sprintf(SNLCPAK_OUTPUT.VERSION_SNANA,      "NULL" );
  sprintf(SNLCPAK_OUTPUT.SURVEY_FILTERS,     "NULL" );
  SNLCPAK_OUTPUT.TEXT_FORMAT[0] = 0;
  
} // end HPAK_INIT


void SNLCPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *VERSION_SNANA,
		  char *SURVEY_FILTERS, int SIMFLAG, 
		  int NFIT_PER_SN, char *TEXT_FORMAT) {

  // One-time global init called before opening output file.
  // NFIT = Number of fits per SN
  //
  // Sep 7 2014: add TEXT_FORMAT (used only if TEXT option is selected)
  // Dec 19 2019: pass SIMFLAG argument.

  char BANNER[100];
  char fnam[] = "SNLCPAK_INIT" ;

  // ---------------- BEGIN ------------

  sprintf(BANNER,"%s for SURVEY=%s  (%d fits per SN)",
	  fnam, SURVEY, NFIT_PER_SN);

  print_banner(BANNER) ;

  SNLCPAK_USE_ROOT  = 0 ; 
  SNLCPAK_USE_TEXT  = 0 ; // Sep 7 2014
  SNLCPAK_USE_HDF5  = 0 ; 
  NCALL_SNLCPAK_FILL = 0 ;

  set_FILTERSTRING(FILTERSTRING);

  // store info in global
  sprintf(SNLCPAK_OUTPUT.SURVEY,             "%s", SURVEY         );
  sprintf(SNLCPAK_OUTPUT.VERSION_PHOTOMETRY, "%s", VERSION_PHOT   );
  sprintf(SNLCPAK_OUTPUT.VERSION_SNANA,      "%s", VERSION_SNANA  );
  sprintf(SNLCPAK_OUTPUT.SURVEY_FILTERS,     "%s", SURVEY_FILTERS );
  sprintf(SNLCPAK_OUTPUT.TEXT_FORMAT,        "%s", TEXT_FORMAT    );

  SNLCPAK_OUTPUT.NFILTDEF_SURVEY = strlen(SNLCPAK_OUTPUT.SURVEY_FILTERS) ;
  SNLCPAK_OUTPUT.SIMFLAG = SIMFLAG ;

  // store max possible fits per SN (can be fewer)
  SNLCPAK_OUTPUT.NFIT_PER_SN = NFIT_PER_SN ;
  if ( NFIT_PER_SN <= 0 ) {
    sprintf(MSGERR1,"Invalid NFIT_PER_SN = %d", NFIT_PER_SN);
    sprintf(MSGERR2,"NFIT_PER_SN must be > 0") ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  // clear arrays to allow for SNLCPAK_XXX calls.
  SNLCPAK_CLEAR_PLOT();     

  // ---------------------------------------------
  // set flag that init has been done.
  SNLCPAK_OUTPUT.INITDONE = SNLCPAK_INITDONE ;


}  // end of SNLCPAK_INIT

void snlcpak_init__(char *SURVEY, char *VER_PHOT, char *VER_SNANA,
		    char *SURVEY_FILTERS, int *SIMFLAG, 
		    int *NFIT_PER_SN, char *TEXTFMT) {
  SNLCPAK_INIT(SURVEY, VER_PHOT, VER_SNANA, SURVEY_FILTERS, 
	       *SIMFLAG, *NFIT_PER_SN, TEXTFMT );
}


// ======================================================
void SNLCPAK_NFIT_PER_SN(int NFIT_PER_SN) {
  // override value passed to SNLCPAK_INIT
  //  char fnam[] = "SNLCPAK_NFIT_PER_SN" ;
  printf("\t SNLCPAK_NFIT_PER_SN:  NFIT_PER_SN -> %d\n", NFIT_PER_SN);
  SNLCPAK_OUTPUT.NFIT_PER_SN = NFIT_PER_SN ;
}

void snlcpak_nfit_per_sn__(int *NFIT_PER_SN) {
  SNLCPAK_NFIT_PER_SN(*NFIT_PER_SN) ;
}

// ======================================================
void SNLCPAK_SURVEY(void) {

  // one-time global init after output file is opened.
  // Store global info such as name of survey, filters, etc ...

  //  char fnam[] = "SNLCPAK_SURVEY" ;

#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT ) { SNLCPAK_SURVEY_ROOT(); }
#endif


} // end of SNLCPAK_SURVEY

void snlcpak_survey__(void) {
  SNLCPAK_SURVEY();
}


void SNLCPAK_CLEAR_PLOT(void) {

  // Called after filling each subdir;
  // clear SNLCPAK counters to prepare for next SUBDIR.
  
  int  iflag, ifilt, i ;
  //  char fnam[] = "SNLCPAK_CLEAR_PLOT" ;

  // init counters
  SNLCPAK_OUTPUT.MXTEXT   = MXTEXT_SNLCPAK ;
  SNLCPAK_OUTPUT.NTEXT    = 0 ;
  SNLCPAK_OUTPUT.NOBS_MAX = 0 ;

  for ( iflag=0; iflag < MXFLAG_SNLCPAK_EPOCH; iflag++ )  {   

    //    if ( NOBS == 0 ) { continue ; }

    SNLCPAK_OUTPUT.NOBS[iflag] = 0 ;   
    free(SNLCPAK_OUTPUT.MJD[iflag]);
    free(SNLCPAK_OUTPUT.TOBS[iflag]);
    free(SNLCPAK_OUTPUT.EPDATA[iflag]);
    free(SNLCPAK_OUTPUT.EPDATA_ERR[iflag]);
    free(SNLCPAK_OUTPUT.IFILTOBS[iflag]);
    free(SNLCPAK_OUTPUT.IFILT[iflag]);

    for ( ifilt=0; ifilt<MXFILTINDX; ifilt++ )  { 
      SNLCPAK_OUTPUT.NOBS_FILT[iflag][ifilt] = 0; 
    }
  }


  for ( iflag=0; iflag < MXFLAG_SNLCPAK_FILTER; iflag++ )  {   
    for ( ifilt=0; ifilt<MXFILTINDX; ifilt++ )  { 
      SNLCPAK_OUTPUT.FILTDATA[iflag][ifilt] = 0.0 ; 
      SNLCPAK_OUTPUT.FILTDATA_ERR[iflag][ifilt] = 0.0 ; 
    }
  }

  // set CCID to blank
  SNLCPAK_OUTPUT.CCID[0] = 0 ;

  // reset plot-filter info
  SNLCPAK_OUTPUT.NFILTDEF_PLOT = 0;
  SNLCPAK_OUTPUT.FILTLIST_PLOT[0] = 0 ;

  // rest dislay text
  for ( i=0; i < MXTEXT_SNLCPAK; i++ ) 
    { sprintf( SNLCPAK_OUTPUT.DISPLAYTEXT[i],"NULL" ) ; }

} // end of SNLCPAK_CLEAR_PLOT

void snlcpak_clear_plot__(void) {
  SNLCPAK_CLEAR_PLOT();
}

void SNLCPAK_CLEAR_SN(void) {
  // May 2014: !!! OBSOLETE !!!
  SNLCPAK_OUTPUT.NLCPAK      = 0 ; // incremental counter 
  SNLCPAK_OUTPUT.NLCPAK_TOT  = 1 ; // expected total number of fits/subdirs
  printf(" xxx SNLCPAK_CLEAR_SN: NLCPAK_TOT -> 1 \n" ); fflush(stdout);
}


void SNLCPAK_DISPLAYTEXT(char *CCID, char *DISPLAYTEXT) {

  // Store text strings to display on plot.
  // NOTHING to do with text output in TEXTFILE_OPEN.

  int N;
  char fnam[] = "SNLCPAK_DISPLAYTEXT" ;

  // ------------- BEGIN -------------

  SNLCPAK_CHECK(CCID,fnam);

  SNLCPAK_OUTPUT.NTEXT++ ;
  N = SNLCPAK_OUTPUT.NTEXT ;

  if ( N >= MXTEXT_SNLCPAK ) {
    sprintf(MSGERR1,"NTEXT exceeds bound of MXTEXT_SNLCPAK=%d",
	    MXTEXT_SNLCPAK);
    sprintf(MSGERR2,"CCID=%s, TEXT='%s'", CCID, DISPLAYTEXT);
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  sprintf( SNLCPAK_OUTPUT.DISPLAYTEXT[N],"%s", DISPLAYTEXT) ;

} // end of SNLCPAK_TEXT

void snlcpak_displaytext__(char *CCID, char *DISPLAYTEXT) {
  SNLCPAK_DISPLAYTEXT(CCID,DISPLAYTEXT);
}

void SNLCPAK_NFIT(int NFIT) {
  // store number of light curve fits (NFIT) for this SN;
  // note that NFIT <= NFIT_PER_SN
  //  char fnam[] = "SNLCPAK_NFIT" ;
  SNLCPAK_OUTPUT.NLCPAK_TOT = NFIT;
  SNLCPAK_OUTPUT.NLCPAK     = 0 ; // reset counter (5/01/2014)
} 
void snlcpak_nfit__(int *NFIT) {
  SNLCPAK_NFIT(*NFIT);
}

void SNLCPAK_DATA(char *CCID, int NOBS, double *MJD, double *TOBS, 
		  double *DATA, double *DATA_ERR, int *IFILTOBS, int FLAG) {

  // store NOBS data values (epochs) passed from main program.
  // FLAG indicates data type:  DATA, REJECT-FLAG, FIT-CHI2, FIT-FUNCTION.
  // The first 3  must have NOBS and TOBS corresponding to the data;
  // the last one (FIT-FUNCTION) has different NOBS and TOBS, usually
  // in very fine TOBS-bins to show a smooth best-fit model function.
  //
  // May 20 2016: new FLAG for simulated flux
  // Jan 28 2026: new FLAG for fit flux at each MJD

  int i, NFILT, IFILT, ifiltobs, flag, MEMD, MEMI ;
  char cfilt[2], comment[80];
  char fnam[] = "SNLCPAK_DATA";

  // -------- BEGIN ---------

  if ( NOPEN_TABLEFILE == 0 ) {
    sprintf(comment, "Cannot store light curve");
    TABLEFILE_notOpen_ABORT(fnam, comment );
  }

  // increment plot-counter for FLUX-DATA since this is the
  // only required FLAG
  if ( FLAG  == SNLCPAK_EPFLAG_FLUXDATA )  
    { SNLCPAK_OUTPUT.NLCPAK++ ; }

  sprintf(comment,"%s(FLAG=%d)", fnam, FLAG);
  SNLCPAK_CHECK(CCID,comment);

  if ( NOBS <= 0  ) {
    sprintf(MSGERR1,"Invalid NOBS = %d for CCID=%s and FLAG=%d . ", 
	    NOBS, CCID, FLAG);
    sprintf(MSGERR2,"NOBS must be > 0. " );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }


  // keep track of max NOBS among all FLAG values
  if ( NOBS >  SNLCPAK_OUTPUT.NOBS_MAX ) 
    { SNLCPAK_OUTPUT.NOBS_MAX = NOBS ; }


  if ( FLAG < MXFLAG_SNLCPAK_EPOCH ) {
    // allocate memory
    MEMD = (NOBS+10) * sizeof(double);
    MEMI = (NOBS+10) * sizeof(int) ;  
    SNLCPAK_OUTPUT.TOBS[FLAG]            = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.MJD[FLAG]             = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.EPDATA[FLAG]          = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.EPDATA_ERR[FLAG]      = (double*)malloc(MEMD) ;
    SNLCPAK_OUTPUT.IFILTOBS[FLAG]        = (int*)malloc(MEMI) ;
    SNLCPAK_OUTPUT.IFILT[FLAG]           = (int*)malloc(MEMI) ;

    SNLCPAK_OUTPUT.NOBS[FLAG] = NOBS ;
    for ( i=0; i < NOBS; i++ ) {
      SNLCPAK_OUTPUT.MJD[FLAG][i]         = MJD[i] ;
      SNLCPAK_OUTPUT.TOBS[FLAG][i]        = TOBS[i] ;
      SNLCPAK_OUTPUT.EPDATA[FLAG][i]      = DATA[i] ;
      SNLCPAK_OUTPUT.EPDATA_ERR[FLAG][i]  = DATA_ERR[i] ;
      SNLCPAK_OUTPUT.IFILTOBS[FLAG][i]    = IFILTOBS[i] ;
    }
  }
  else {
    NFILT = NOBS ;  // one value per filter
    for (i=0; i < NFILT; i++ ) {
      ifiltobs = IFILTOBS[i] ;
      sprintf(cfilt,"%c", FILTERSTRING[ifiltobs] );
      IFILT = strcspn( SNLCPAK_OUTPUT.SURVEY_FILTERS, cfilt) ;
      flag  = FLAG - 100; // storage index

      if ( flag < 0 || flag >= MXFLAG_SNLCPAK_FILTER ) {
	sprintf(MSGERR1,"Invalid flag=%d for FLAG=%d", flag, FLAG);
	sprintf(MSGERR2,"IFILT=%d  IFILTOBS=%d(%s)", IFILT, ifiltobs, cfilt);
	errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
      }

      SNLCPAK_OUTPUT.FILTDATA[flag][IFILT]      = DATA[i] ; 
      SNLCPAK_OUTPUT.FILTDATA_ERR[flag][IFILT]  = DATA_ERR[i] ;

    }    // i loop
  }


}  // end of SNLCPAK_DATA

void snlcpak_data__(char *CCID, int *NOBS, double *MJD, double *TOBS, 
		    double *DATA, double *DATA_ERR, int *IFILTOBS,int *FLAG) {
  SNLCPAK_DATA(CCID, *NOBS, MJD, TOBS, DATA, DATA_ERR, IFILTOBS, *FLAG);
}

// ========================================
void SNLCPAK_FILL(char *CCID) {
 
  char fnam[] = "SNLCPAK_FILL" ;

  // ----------- BEGIN ------------

  SNLCPAK_CHECK(CCID,fnam);

  // do some global stuff needed for any output format.
  SNLCPAK_FILL_PREP();
  NCALL_SNLCPAK_FILL++ ;  

#ifdef USE_ROOT
  if ( SNLCPAK_USE_ROOT  ) { SNLCPAK_FILL_ROOT(); }
#endif


#ifdef USE_TEXT
  if ( SNLCPAK_USE_TEXT  ) { SNLCPAK_FILL_TEXT(); }
#endif

 
  SNLCPAK_CLEAR_PLOT();

}  // end  SNLCPAK_FILL

void snlcpak_fill__(char *CCID) {
  SNLCPAK_FILL(CCID);
}

// ====================================
void SNLCPAK_FILL_PREP() {

  // compute quantities useful for for the _FILL_ routines.

  int  IFILT, FLAG, obs, NFILT, USE, NOBS, IFILTOBS, IFMIN, IFMAX ;
  char cfilt[2], *SURVEY_FILTERS ;
  //  char fnam[] = "SNLCPAK_FILL_PREP" ;

  // ------------- BEGIN --------------

  for ( FLAG=0; FLAG < MXFLAG_SNLCPAK_EPOCH; FLAG++ ) {
    for ( IFILT=0; IFILT < MXFILTINDX; IFILT++ )  { 
      SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] = 0 ; 
    }
  }


  // use required data array to check for which filters to plot.
  // Filter order is determined by SURVEY_FILTERS.

  SURVEY_FILTERS = SNLCPAK_OUTPUT.SURVEY_FILTERS ;

  for ( FLAG=0; FLAG < MXFLAG_SNLCPAK_EPOCH; FLAG++ ) {
    NOBS      = SNLCPAK_OUTPUT.NOBS[FLAG] ;
    for(obs=0; obs<NOBS; obs++ ) {   

      IFILTOBS = SNLCPAK_OUTPUT.IFILTOBS[FLAG][obs] ;
      sprintf(cfilt,"%c", FILTERSTRING[IFILTOBS] );
      IFILT = strcspn( SURVEY_FILTERS, cfilt);
      SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] = 1 ;  
      
      SNLCPAK_OUTPUT.IFILT[FLAG][obs] = IFILT ;  // sparse IFILT for this obs.
      SNLCPAK_OUTPUT.NOBS_FILT[FLAG][IFILT]++ ; //  NOBS for this filter

    }  // obs loop
  }  // FLAG loop


  // count how many defined filters for DATA, and set FILTLIST_PLOT
  NFILT  = 0 ;
  IFMIN  = IFMAX = -9 ;
  FLAG   = SNLCPAK_EPFLAG_FLUXDATA ;
  for ( IFILT=0; IFILT < MXFILTINDX; IFILT++ ) {
    USE = SNLCPAK_OUTPUT.USEFILT[FLAG][IFILT] ;
    sprintf(cfilt,"%c", SURVEY_FILTERS[IFILT] );
    if ( USE ) { 
      NFILT++ ; 
      strcat(SNLCPAK_OUTPUT.FILTLIST_PLOT, cfilt); 

      if ( IFMIN < 0 ) { IFMIN = IFILT ; }
      IFMAX = IFILT ;
    }
  }

  SNLCPAK_OUTPUT.NFILTDEF_PLOT = NFILT ;
  SNLCPAK_OUTPUT.IFILT_MIN     = IFMIN ;
  SNLCPAK_OUTPUT.IFILT_MAX     = IFMAX ;

} // end of SNLCPAK_FILL_PREP



// ===========================================
void SNLCPAK_CHECK(char *CCID, char *callFun) {

  // Either store CCID in global array,
  // or check the CCID matches what is in global array.
  // This CHECK prevents switching CIDs before calling SNLCPAK_FILL.

  char *ptrCCID ;

  char fnam0[] = "SNLCPAK_CHECK";
  char fnam[100];
  concat_callfun_plus_fnam(callFun, fnam0, fnam);
  // ------------ BEGIN -----------

  // make sure that SNLCPAK_INIT was called
  if ( SNLCPAK_OUTPUT.INITDONE != SNLCPAK_INITDONE ) {
    sprintf(MSGERR1,"%s called before SNLCPAK_INIT.", callFun);
    sprintf(MSGERR2,"Must call SNLCPAK_INIT first.");
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }


  ptrCCID = SNLCPAK_OUTPUT.CCID ;

  // store CCID or make sure that CCID is correct
  if ( strlen(ptrCCID) == 0 ) {
    // define CCID
    sprintf(ptrCCID,"%s", CCID);    
  }
  else if ( strcmp(CCID,ptrCCID) != 0 ) {
    // abort on error
    sprintf(MSGERR1,"Invalid CCID = '%s' passed to %s", CCID, callFun);
    sprintf(MSGERR2,"Expected CCID = '%s'", ptrCCID) ;
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

} // end of SNLCPAK_CHECK



// ==========================================
//    SPECPAK functions (April 2019)
//    Only for TEXT output to feed external plotter. Ignore ROOT
// ==========================================


void SPECPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *TEXT_FORMAT) {

  // April 2019
  // One-time global init called before opening output file.


  char BANNER[100];
  char fnam[] = "SPECPAK_INIT" ;

  // ---------------- BEGIN ------------

  sprintf(BANNER,"%s for SURVEY=%s ", fnam, SURVEY);

  print_banner(BANNER) ;

  SPECPAK_USE_ROOT  = 0 ; 
  SPECPAK_USE_TEXT  = 0 ; // Sep 7 2014
  SPECPAK_USE_MARZ  = 0 ;
  NCALL_SPECPAK_FILL = 0 ;

  // store info in global
  sprintf(SPECPAK_OUTPUT.SURVEY,             "%s", SURVEY         );
  sprintf(SPECPAK_OUTPUT.VERSION_PHOTOMETRY, "%s", VERSION_PHOT   );
  sprintf(SPECPAK_OUTPUT.TEXT_FORMAT,        "%s", TEXT_FORMAT    );

  // clear arrays
  SPECPAK_CLEAR_PLOT();     

  // ---------------------------------------------
  // set flag that init has been done.
  SPECPAK_OUTPUT.INITDONE = SPECPAK_INITDONE ;


}  // end of SPECPAK_INIT

void specpak_init__(char *SURVEY, char *VER_PHOT, char *TEXTFMT) {
  SPECPAK_INIT(SURVEY, VER_PHOT, TEXTFMT );
}


// =========================================
void SPECPAK_CLEAR_PLOT(void) {

  // May 13 2019:
  //  + fix bug: NCALL_SPECPAK_FILL>10 --> >0 so that
  //    allocated memory is freed. Fixes valgrind errors.

  //  char fnam[] = "SPECPAK_CLEAR_PLOT" ;

  //  printf(" xxx %s  NCALL=%d \n", fnam, NCALL_SPECPAK_FILL ) ;

  SPECPAK_OUTPUT.NSPEC = 0 ;
  SPECPAK_OUTPUT.NLAMBIN_TOT = 0 ;

  if ( NCALL_SPECPAK_FILL > 0 ) {
    free(SPECPAK_OUTPUT.ID    ) ;
    free(SPECPAK_OUTPUT.LAMMIN) ;
    free(SPECPAK_OUTPUT.LAMMAX) ;
    free(SPECPAK_OUTPUT.FLAM) ;
    free(SPECPAK_OUTPUT.FLAMERR) ;
  }


} // end SPECPAK_CLEAR_PLOT

void specpak_clear_plot__(void) {  SPECPAK_CLEAR_PLOT();  }


// =======================================================
void SPECPAK_DATA(char *CCID, int ID, double MJD, double Tobs, double Texpose,
		  int NLAMBIN, double *LAMMIN, double *LAMMAX, 
		  double *FLAM, double *FLAMERR)
{

  // pack/store multiple spectra for one event.
  // Note that MJD = Tobs = -9 for host galaxy.

  int MEMD, MEMI, i, ii, NSPEC, NLAMTOT, ILAM_START ;
  char fnam[] = "SPECPAK_DATA" ;

  // ------------ BEGIN ------------

  if ( NLAMBIN <= 0  ) {
    sprintf(MSGERR1,"Invalid NLAMBIN = %d for CCID=%s and ID=%d . ", 
	    NLAMBIN, CCID, ID);
    sprintf(MSGERR2,"NLAMBIN must be > 0. " );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);    
  }

  NSPEC      = SPECPAK_OUTPUT.NSPEC;
  ILAM_START = SPECPAK_OUTPUT.NLAMBIN_TOT;
  SPECPAK_OUTPUT.NLAMBIN_TOT += NLAMBIN ;
  NLAMTOT    = SPECPAK_OUTPUT.NLAMBIN_TOT ;

  sprintf(SPECPAK_OUTPUT.CCID, "%s", CCID );
  SPECPAK_OUTPUT.ID_LIST[NSPEC]       = ID ;
  SPECPAK_OUTPUT.NLAMBIN_LIST[NSPEC]  = NLAMBIN ;
  SPECPAK_OUTPUT.MJD_LIST[NSPEC]      = MJD;
  SPECPAK_OUTPUT.TOBS_LIST[NSPEC]     = Tobs ;
  SPECPAK_OUTPUT.TEXPOSE_LIST[NSPEC]  = Texpose ;

  MEMD    = NLAMTOT * sizeof(double);
  MEMI    = NLAMTOT * sizeof(int) ;  

  // malloc wave-dependent storage
  if ( NSPEC == 0 ) {
    SPECPAK_OUTPUT.ID      = (int   *)malloc(MEMI);
    SPECPAK_OUTPUT.LAMMIN  = (double*)malloc(MEMD); // array of lammin in each bin
    SPECPAK_OUTPUT.LAMMAX  = (double*)malloc(MEMD); // array of lammax in each bin
    SPECPAK_OUTPUT.FLAM    = (double*)malloc(MEMD); // array of flam in each bin
    SPECPAK_OUTPUT.FLAMERR = (double*)malloc(MEMD); 
  }
  else {
    // realloc for more wave bins
    SPECPAK_OUTPUT.ID      = (int   *)realloc(SPECPAK_OUTPUT.ID,     MEMI );
    SPECPAK_OUTPUT.LAMMIN  = (double*)realloc(SPECPAK_OUTPUT.LAMMIN, MEMD );
    SPECPAK_OUTPUT.LAMMAX  = (double*)realloc(SPECPAK_OUTPUT.LAMMAX, MEMD );
    SPECPAK_OUTPUT.FLAM    = (double*)realloc(SPECPAK_OUTPUT.FLAM,   MEMD );
    SPECPAK_OUTPUT.FLAMERR = (double*)realloc(SPECPAK_OUTPUT.FLAMERR,MEMD );
  }

  // load SPECPAK_OUTPUT
  for(i=0; i < NLAMBIN; i++ ) {
    ii = i + ILAM_START ;
    SPECPAK_OUTPUT.ID[ii]      = ID ;
    SPECPAK_OUTPUT.LAMMIN[ii]  = LAMMIN[i];
    SPECPAK_OUTPUT.LAMMAX[ii]  = LAMMAX[i];
    SPECPAK_OUTPUT.FLAM[ii]    = FLAM[i];
    SPECPAK_OUTPUT.FLAMERR[ii] = FLAMERR[i];
  }

  SPECPAK_OUTPUT.NSPEC++ ;
  if ( SPECPAK_OUTPUT.NSPEC >= MXSPEC_SPECPAK ) {
    sprintf(MSGERR1, "NSPEC exceeds MXSPEC_SPECPAK=%d bound.", 
	    MXSPEC_SPECPAK);
    sprintf(MSGERR2, "Reduce NSPEC or increase MXSPEC_SPECPAK"); 
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2);
  }

  return ;

} // end SPECPAK_DATA


void specpak_data__(char *CCID, int *ID, double *MJD,double *Tobs,
		    double *Texpose,int *NLAMBIN,double *LAMMIN,double *LAMMAX,
		    double *FLAM, double *FLAMERR) {
  SPECPAK_DATA(CCID, *ID, *MJD, *Tobs, *Texpose,
	       *NLAMBIN, LAMMIN, LAMMAX, FLAM, FLAMERR);
}


// ========================================
void SPECPAK_FILL(char *CCID) {
 
  //  char fnam[] = "SPECPAK_FILL" ;

  // ----------- BEGIN ------------

  NCALL_SPECPAK_FILL++ ; 

  // do nothing for ROOT ... only for TEXT.

#ifdef USE_TEXT
  if ( SPECPAK_USE_TEXT  ) { SPECPAK_FILL_TEXT(); }
#endif

#ifdef USE_MARZ
  if ( SPECPAK_USE_MARZ  ) { SPECPAK_FILL_MARZ(); }
#endif


  SPECPAK_CLEAR_PLOT(); 

}  // end  SPECPAK_FILL

void specpak_fill__(char *CCID) {  SPECPAK_FILL(CCID); }

// END
