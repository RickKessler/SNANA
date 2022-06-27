/******************************************
  Created Dec 10 2021 by R.Kessler

  Data-handling tools that don't depend on format.
  Format-specific tools are in 
     sntools_dataformat_fits.c
     sntools_dataformat_text.c

  Initial motivation is to move data-override utility out of
  fortran code and into here, and to make it more flexible so 
  that hard-coding variable match-list isn't needed.

******************************************/

#include  "sntools.h"
#include  "sntools_spectrograph.h"  
#include  "sntools_data.h"
#include  "sntools_output.h"
#include  "sntools_cosmology.h"

//#include  "sntools_dataformat_text.h"
#include  "sntools_host.h" 
//#include  "sntools_trigger.h" 


// =======================================================
bool IS_SIMKEY_SNDATA(char *key) {

  // Created Mar 8 2022
  // Return true if input *key is from simulation,
  // where key is a data key in either FITS or TEXT format.

  bool IS_KEYSIM = false;
  // ----------- BEGIN ------------

  if ( strncmp(key,"SIM",3)  ==   0 ) { IS_KEYSIM = true; }
  if ( strstr(key,"SIMSED") != NULL ) { IS_KEYSIM = true; }
  if ( strstr(key,"LCLIB" ) != NULL ) { IS_KEYSIM = true; }

  return IS_KEYSIM ;

} // end IS_SIMKEY_SNDATA


// ******************************************
void copy_int(int copyFlag, double *DVAL0, int *IVAL1) {
  if   ( copyFlag > 0)  { *IVAL1 = (int)(*DVAL0);  }  
  else                  { *DVAL0 = (double)(*IVAL1);  }
}

void copy_lli(int copyFlag, double *DVAL0, long long *IVAL1) {
  if   ( copyFlag > 0)  { *IVAL1 = (long long)(*DVAL0);  }  
  else                  { *DVAL0 = (double)(*IVAL1);  }
}

void copy_flt(int copyFlag, double *DVAL0, float *FVAL1) {
  if   ( copyFlag > 0)  { *FVAL1 = (float)(*DVAL0);  }  
  else                  { *DVAL0 = (double)(*FVAL1);  }
}

void copy_dbl(int copyFlag, double *DVAL0, double *DVAL1) {
  if   ( copyFlag > 0)  { *DVAL1 = (double)(*DVAL0);  }  
  else                  { *DVAL0 = (double)(*DVAL1);  }
}

void copy_str(int copyFlag, char *STR0, char *STR1) {
  if   ( copyFlag > 0)  { sprintf(STR1, "%s", STR0); }
  else                  { sprintf(STR0, "%s", STR1); }
}


// ===================================================
void copy_SNDATA_GLOBAL(int copyFlag, char *key, int NVAL, 
			char *stringVal, double *parVal ) {

  // Created Feb 2021
  //
  // Copy *key value to or from SNDATA struct.
  // Note that event has already been read, but is not stored
  // in useful arrays.
  //
  // GLOBAL refers to params that do not depend on event.
  //
  // Inputs:
  //  copyFlag : 
  //     +1 -> copy from string or parVal to SNDATA (prep  data write)
  //     -1 -> copy from SNDATA to string or parVal (after read data)
  //
  //   *key  : name of variable to copy to/from SNDATA struct
  //   NVAL  : number of values to copy
  // 
  // Output:
  //    *stringVal  : string value if *key points to string
  //    *parVal     : double value if *key points to double,float,int
  //
  // Apr 24 2021: add SIM_BIASCOR_MASK
  // Oct 08 2021: add SIM_MODEL_INDEX

  bool ISKEY_PRIVATE = ( strstr (key,"PRIVATE")   != NULL ) ;
  bool ISKEY_BYOSED  = ( strncmp(key,"BYOSED",6)  == 0 ) ;
  bool ISKEY_SNEMO   = ( strncmp(key,"SNEMO",5)   == 0 ) ;
  bool ISKEY_SIMSED  = ( strncmp(key,"SIMSED",6)  == 0 ) ;
  bool ISKEY_LCLIB   = ( strncmp(key,"LCLIB",5)   == 0 ) ;
  bool ISKEY_SIM     = ( strncmp(key,"SIM",3)     == 0 && !ISKEY_SIMSED) ;
  bool ISKEY_ZPHOT_Q = ( strstr (key,"ZPHOT_Q")  != NULL ) ;

  int ivar, NVAR, ipar, PCT ;
  char fnam[] = "copy_SNDATA_GLOBAL" ;

  // --------------- BEGIN --------------

  if(copyFlag<0) { sprintf(stringVal,"NOTSET");  parVal[0] = -999.0; }

  if ( strcmp(key,"SNANA_VERSION") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA.SNANA_VERSION );  }

  else if ( strcmp(key,"SURVEY") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA.SURVEY_NAME );  }

  else if ( strcmp(key,"SUBSURVEY_FLAG") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.SUBSURVEY_FLAG );  }

  else if ( strcmp(key,"FILTERS") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA_FILTER.LIST );  }

  else if ( strcmp(key,"DATATYPE") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA.DATATYPE );  }

  else if ( strcmp(key,"PySEDMODEL") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA.PySEDMODEL_NAME );  }

  else if ( strcmp(key,"NXPIX") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.NXPIX );  }
  else if ( strcmp(key,"NYPIX") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.NYPIX );  }

  else if ( ISKEY_PRIVATE  ) {
    if ( strcmp(key,"NVAR_PRIVATE") == 0 ) { 
      copy_int(copyFlag, parVal, &SNDATA.NVAR_PRIVATE );  
    }
    else {
      sscanf(&key[7], "%d", &ivar);  // PRIVATEnn
      copy_str(copyFlag, stringVal, SNDATA.PRIVATE_KEYWORD[ivar] ); 
    }
  }

  else if ( ISKEY_ZPHOT_Q ) {

    if ( strcmp(key,"NZPHOT_Q") == 0 ) {
      copy_int(copyFlag, parVal, &SNDATA.HOSTGAL_NZPHOT_Q );
    }
    else if ( strstr(key,"PERCENTILE_ZPHOT_Q") != NULL ) {
      int lentmp = strlen(key) - 2 ; // PERCENTILE_ZPHOT_Qnn
      sscanf(&key[lentmp], "%d", &ivar); 
      PCT = SNDATA.HOSTGAL_PERCENTILE_ZPHOT_Q[ivar] ;
      copy_int(copyFlag, parVal, &PCT ); 
    }
  }

  else if ( ISKEY_SIMSED  ) {

    if ( strcmp(key,"SIMSED_NPAR") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.NPAR_SIMSED );  }
    else {
      sscanf(&key[10], "%d", &ivar);  // 'SIMSED_PARnn'
      copy_str(copyFlag, stringVal, SNDATA.SIMSED_KEYWORD[ivar] ); 
    }
  }

  else if ( ISKEY_LCLIB  ) {
    if ( strcmp(key,"LCLIB_NPAR") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.NPAR_LCLIB );  }
    else {
      sscanf(&key[9], "%d", &ivar); // LCLIB_PARnn
      copy_str(copyFlag, stringVal, SNDATA.LCLIB_KEYWORD[ivar] ); 
    }
  }

  else if ( ISKEY_BYOSED  ) {
    if ( strcmp(key,"BYOSED_NPAR") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.NPAR_PySEDMODEL );  }
    else {
      sscanf(&key[10], "%d", &ivar); // BYOSED_PARnn
      copy_str(copyFlag, stringVal, SNDATA.PySEDMODEL_KEYWORD[ivar] ); 
    }
  }

  else if ( ISKEY_SNEMO  ) {
    if ( strcmp(key,"SNEMO_NPAR") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.NPAR_PySEDMODEL );  }
    else {
      sscanf(&key[9], "%d", &ivar); // SNEMO_PARnn
      copy_str(copyFlag, stringVal, SNDATA.PySEDMODEL_KEYWORD[ivar] ); 
    }
  }

  else if ( ISKEY_SIM ) {   

    if ( strcmp(key,"SIMLIB_FILE") == 0 ) 
      { copy_str(copyFlag, stringVal, SNDATA.SIMLIB_FILE ); }

    else if ( strcmp(key,"SIMLIB_MSKOPT") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIMLIB_MSKOPT ); }

    else if ( strcmp(key,"SIMOPT_MWCOLORLAW") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIMOPT_MWCOLORLAW ); }

    else if ( strcmp(key,"SIMOPT_MWEBV") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIMOPT_MWEBV ); }

    else if ( strcmp(key,"SIM_MWRV") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_MWRV ); }

    else if ( strcmp(key,"SIM_VARNAME_SNRMON") == 0 ) 
      { copy_str(copyFlag, stringVal, SNDATA.VARNAME_SNRMON ); }

    else if ( strcmp(key,"SIM_HOSTLIB_NPAR") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.NPAR_SIM_HOSTLIB ); }

    else if ( strncmp(key,"SIM_HOSTLIB_PAR",15) == 0 ) {
      for(ipar=0; ipar < SNDATA.NPAR_SIM_HOSTLIB; ipar++ ) {
	copy_str(copyFlag, stringVal, SNDATA.SIM_HOSTLIB_KEYWORD[ipar]);
      }
    }
    else if ( strcmp(key,"SIM_SL_FLAG") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_SL_FLAG ); }

    else if ( strcmp(key,"SIM_BIASCOR_MASK") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_BIASCOR_MASK ); }

    else if ( strcmp(key,"SIM_MODEL_INDEX") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_MODEL_INDEX );  }

    else {
      // error message
      sprintf(c1err,"Unknown SIM key = %s (copyFlag=%d)", key, copyFlag);
      sprintf(c2err,"stringVal='%s'  parVal=%f", stringVal, parVal[0] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

  }
  else {
    // error message
    sprintf(c1err,"Unknown key = %s  (copyFlag=%d)", key, copyFlag);
    sprintf(c2err,"stringVal='%s'  parVal=%f", stringVal, parVal[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return;

} // end copy_SNDATA_GLOBAL


// = = = = = = = = = = = = = = = = = = = = = = = =

void copy_SNDATA_HEAD(int copyFlag, char *key, int NVAL, 
		      char *stringVal, double *parVal ) {

  // Created Feb 2021:
  //
  // Copy *key value to or from SNDATA struct.
  // Note that event has already been read, but is not stored
  // in useful arrays.
  //
  // HEAD refers SN-dependent HEADER that updates for each event.
  //   (i.e., depends on SN; does NOT depend on MJD)

  // Inputs:
  //  copyFlag : 
  //     +1 -> copy from string or parVal to SNDATA (prep  data write)
  //     -1 -> copy from SNDATA to string or parVal (after read data)
  //
  //   *key  : name of variable to copy to/from SNDATA struct
  //   NVAL  : number of values to copy
  // 
  // Output:
  //    *stringVal  : string value if *key points to string
  //    *parVal     : double value if *key points to double,float,int
  //
  // History:
  // Sep 9 2021 A. Gagliano: Load igal dimension of SNDATA.SIM_HOSTLIB_PARVAL

  int NFILT = SNDATA_FILTER.NDEF ;
  char *PySEDMODEL_NAME = SNDATA.PySEDMODEL_NAME ; // BYOSED or SNEMO
  int  len_PySEDMODEL   = strlen(PySEDMODEL_NAME);
  int  ncmp_PySEDMODEL  = strncmp(key,PySEDMODEL_NAME,len_PySEDMODEL) ;
  int igal, NGAL, ifilt, ifilt_obs, NVAR, ivar, ipar, q, PCT;
  double DVAL;
  char PREFIX[40], KEY_TEST[60], cfilt[2] ;
  char fnam[] = "copy_SNDATA_HEAD" ;

  // ------------- BEGIN ------------

  if ( copyFlag < 0 ) 
    { sprintf(stringVal,"NOTSET");  parVal[0] = -999.0; }


  if ( strcmp(key,"SUBSURVEY") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA.SUBSURVEY_NAME ); }

  else if ( strcmp(key,"SNID") == 0 ) {
     copy_str(copyFlag, stringVal, SNDATA.CCID ); 
     SNDATA.NOBS_STORE = 0 ; // must call select_MJD_SNDATA to set this
  }

  else if ( strcmp(key,"IAUC") == 0 ) 
    { copy_str(copyFlag, stringVal, SNDATA.IAUC_NAME ); }

  else if ( strcmp(key,"FAKE") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.FAKE ); }

  else if ( strcmp(key,"MASK_FLUXCOR_SNANA") == 0 )
    { copy_int(copyFlag, parVal, &SNDATA.MASK_FLUXCOR ); } 

  else if ( strcmp(key,"RA") == 0 ) 
    { copy_dbl(copyFlag, parVal, &SNDATA.RA ); }
  else if ( strcmp(key,"DEC") == 0 ) 
    { copy_dbl(copyFlag, parVal, &SNDATA.DEC ); }

  else if ( strcmp(key,"PIXSIZE") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.PIXSIZE ); }
  else if ( strcmp(key,"NXPIX") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.NXPIX ); } 
  else if ( strcmp(key,"NYPIX") == 0 )
    { copy_int(copyFlag, parVal, &SNDATA.NYPIX ); } 

  else if ( strcmp(key,"CCDNUM") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.CCDNUM[0] ); } // should be obsolete

  else if ( strcmp(key,"SNTYPE") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.SNTYPE ); } 

  else if ( strcmp(key,"NOBS") == 0 ) { 
    copy_int(copyFlag, parVal, &SNDATA.NOBS ); 
    SNDATA.NOBS_STORE = -1 ; // must call select_MJD_SNDATA to set this
  } 

  else if ( strcmp(key,"MWEBV") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.MWEBV ); } 
  else if ( strcmp(key,"MWEBV_ERR") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.MWEBV_ERR ); } 

  else if ( strcmp(key,"REDSHIFT_HELIO") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.REDSHIFT_HELIO ); } 
  else if ( strcmp(key,"REDSHIFT_HELIO_ERR") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.REDSHIFT_HELIO_ERR ); } 

  else if ( strcmp(key,"REDSHIFT_FINAL") == 0 || 
	    strcmp(key,"REDSHIFT_CMB"  ) == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.REDSHIFT_FINAL ); } 
  else if ( strcmp(key,"REDSHIFT_FINAL_ERR") == 0 ||
	    strcmp(key,"REDSHIFT_CMB_ERR"  ) == 0 )  
    { copy_flt(copyFlag, parVal, &SNDATA.REDSHIFT_FINAL_ERR ); } 
  else if ( strcmp(key,"REDSHIFT_QUALITYFLAG") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.REDSHIFT_QUALITYFLAG ); } 

  else if ( strcmp(key,"VPEC") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.VPEC ); } 
  else if ( strcmp(key,"VPEC_ERR") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.VPEC_ERR ); } 

  else if ( strncmp(key,"HOSTGAL",7) == 0 ) {

    if ( strcmp(key,"HOSTGAL_NMATCH") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.HOSTGAL_NMATCH[0] ); } 
    else if ( strcmp(key,"HOSTGAL_NMATCH2") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.HOSTGAL_NMATCH[1] ); } 
    else if ( strcmp(key,"HOSTGAL_CONFUSION") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_CONFUSION ); } 

    if ( strstr(key,"SB_FLUXCAL") != NULL ) {
      for(ifilt=0; ifilt < NFILT; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST,"HOSTGAL_SB_FLUXCAL_%c", FILTERSTRING[ifilt_obs]); 
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_SB_FLUXCAL[ifilt]); } 
      }
    }

    NGAL = MXHOSTGAL ;
    for(igal=0; igal < NGAL; igal++ ) {
      sprintf(PREFIX,"HOSTGAL");
      if ( igal > 0 ) { sprintf(PREFIX,"HOSTGAL%d",igal+1); }

      sprintf(KEY_TEST,"%s_OBJID", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_lli(copyFlag, parVal, &SNDATA.HOSTGAL_OBJID[igal] ); } 

      sprintf(KEY_TEST,"%s_FLAG", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_int(copyFlag, parVal, &SNDATA.HOSTGAL_FLAG[igal] ); } 

      sprintf(KEY_TEST,"%s_PHOTOZ", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_PHOTOZ[igal] ); } 
      sprintf(KEY_TEST,"%s_PHOTOZ_ERR", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_PHOTOZ_ERR[igal] ); } 

      sprintf(KEY_TEST,"%s_SPECZ", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_SPECZ[igal] ); } 
     
      sprintf(KEY_TEST,"%s_SPECZ_ERR", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_SPECZ_ERR[igal] ); } 

      sprintf(KEY_TEST,"%s_RA", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_dbl(copyFlag, parVal, &SNDATA.HOSTGAL_RA[igal] ); } 
      sprintf(KEY_TEST,"%s_DEC", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_dbl(copyFlag, parVal, &SNDATA.HOSTGAL_DEC[igal] ); } 

      sprintf(KEY_TEST,"%s_SNSEP", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_SNSEP[igal] ); } 

      sprintf(KEY_TEST,"%s_DDLR", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_DDLR[igal] ); } 

      sprintf(KEY_TEST,"%s_LOGMASS", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_LOGMASS_OBS[igal] ); } 
      sprintf(KEY_TEST,"%s_LOGMASS_ERR", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_LOGMASS_ERR[igal] ); } 

      sprintf(KEY_TEST,"%s_LOGSFR", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_LOGSFR_OBS[igal] ); } 
      sprintf(KEY_TEST,"%s_LOGSFR_ERR", PREFIX); 
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_LOGSFR_ERR[igal] ); } 

      sprintf(KEY_TEST,"%s_LOGsSFR", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 )
        { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_LOGsSFR_OBS[igal] ); }
      sprintf(KEY_TEST,"%s_LOGsSFR_ERR", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 )
        { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_LOGsSFR_ERR[igal] ); }

      sprintf(KEY_TEST,"%s_COLOR", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 )
        { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_COLOR_OBS[igal] ); }
      sprintf(KEY_TEST,"%s_COLOR_ERR", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 )
        { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_COLOR_ERR[igal] ); }

      for(ifilt=0; ifilt < NFILT; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST,"%s_MAG_%c", PREFIX, FILTERSTRING[ifilt_obs]); 
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_MAG[igal][ifilt]); } 
      }

      if ( strstr(key,PREFIX_ZPHOT_Q) != NULL ) {
	for(q=0; q < SNDATA.HOSTGAL_NZPHOT_Q; q++ ) {
	  PCT = SNDATA.HOSTGAL_PERCENTILE_ZPHOT_Q[q] ;
	  sprintf(KEY_TEST,"%s_%s%3.3d", PREFIX, PREFIX_ZPHOT_Q, PCT);
	  if ( strcmp(key,KEY_TEST) == 0 ) 
	    { copy_flt(copyFlag, parVal, &SNDATA.HOSTGAL_ZPHOT_Q[igal][q]);  } 
	}
      }

    } // end igal

  }  // end "HOSTGAL"

  else if ( strcmp(key,"PEAKMJD") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.SEARCH_PEAKMJD) ; }  

  else if ( strcmp(key,"MJD_TRIGGER") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.MJD_TRIGGER) ; }  
  else if ( strcmp(key,"MJD_DETECT_FIRST") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.MJD_DETECT_FIRST) ; }  
  else if ( strcmp(key,"MJD_DETECT_LAST") == 0 ) 
    { copy_flt(copyFlag, parVal, &SNDATA.MJD_DETECT_LAST) ; }  

  else if ( strcmp(key,"SEARCH_TYPE") == 0 ) 
    { copy_int(copyFlag, parVal, &SNDATA.SEARCH_TYPE) ; }  

  else if ( strncmp(key,"PRIVATE",7) == 0 ) {
    NVAR = SNDATA.NVAR_PRIVATE ;
    for(ivar=1; ivar <= NVAR; ivar++ ) {
      if ( strcmp(key,SNDATA.PRIVATE_KEYWORD[ivar]) == 0 ) 
	{ copy_dbl(copyFlag, parVal, &SNDATA.PRIVATE_VALUE[ivar]) ; }  
    }    
  }

  else if ( strncmp(key,"SIM",3) == 0 ) {
    // ------ SIM ------- 

    if ( strcmp(key,"SIM_MODEL_NAME") == 0 ) 
      { copy_str(copyFlag, stringVal, SNDATA.SIM_MODEL_NAME ); }
    else if ( strcmp(key,"SIM_TYPE_NAME") == 0 ) 
      { copy_str(copyFlag, stringVal, SNDATA.SIM_TYPE_NAME ); }

    else if ( strcmp(key,"SIM_MODEL_INDEX") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_MODEL_INDEX) ; }  

    else if ( strcmp(key,"SIM_TYPE_INDEX") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_TYPE_INDEX) ; }  

    else if ( strcmp(key,"SIM_TEMPLATE_INDEX") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_TEMPLATE_INDEX) ; }  

    else if ( strcmp(key,"SIM_SUBSAMPLE_INDEX") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SUBSAMPLE_INDEX) ; }  

    else if ( strcmp(key,"SIM_LIBID") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_LIBID) ; }  

    else if ( strcmp(key,"SIM_NGEN_LIBID") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_NGEN_LIBID) ; }  

    else if ( strcmp(key,"SIM_NOBS_UNDEFINED") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_NOBS_UNDEFINED) ; }  

    else if ( strcmp(key,"SIM_SEARCHEFF_MASK") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_SEARCHEFF_MASK) ; }  

    else if ( strcmp(key,"SIM_REDSHIFT_HELIO") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_REDSHIFT_HELIO) ; }  
    else if ( strcmp(key,"SIM_REDSHIFT_CMB") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_REDSHIFT_CMB) ; }  
    else if ( strcmp(key,"SIM_REDSHIFT_HOST") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_REDSHIFT_HOST) ; }  
    else if ( strcmp(key,"SIM_REDSHIFT_FLAG") == 0 ) 
      { copy_int(copyFlag, parVal, &SNDATA.SIM_REDSHIFT_FLAG) ; }  
    else if ( strcmp(key,"SIM_VPEC") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_VPEC) ; }  

    else if ( strcmp(key,"SIM_HOSTLIB_GALID") == 0 ) 
      { copy_lli(copyFlag, parVal, &SNDATA.SIM_HOSTLIB_GALID) ; }  

    else if ( strncmp(key,"SIM_HOSTLIB",11) == 0 ) {
      for(igal=0; igal < MXHOSTGAL; igal++ ) {
      	for(ipar=0; ipar < SNDATA.NPAR_SIM_HOSTLIB; ipar++ ) {
		copy_flt(copyFlag, parVal, &SNDATA.SIM_HOSTLIB_PARVAL[ipar][igal]) ; 
        }
      }
    }
    else if ( strcmp(key,"SIM_DLMU") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_DLMU) ; }  

    else if ( strcmp(key,"SIM_LENSDMU") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_LENSDMU) ; }  

    else if ( strcmp(key,"SIM_RA") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_RA) ; }  

    else if ( strcmp(key,"SIM_DEC") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_DEC) ; }  

    else if ( strcmp(key,"SIM_MWEBV") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_MWEBV) ; }  

    else if ( strcmp(key,"SIM_PEAKMJD") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_PEAKMJD) ; }  

    else if ( strcmp(key,"SIM_MAGSMEAR_COH") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_MAGSMEAR_COH) ; }  

    else if ( strcmp(key,"SIM_AV") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_AV) ; }  
    else if ( strcmp(key,"SIM_RV") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_RV) ; }  

    else if ( strncmp(key,"SIM_SALT2",9) == 0 ) {
      if ( strcmp(key,"SIM_SALT2x0") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2x0) ; }  
      else if ( strcmp(key,"SIM_SALT2x1") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2x1) ; }  
      else if ( strcmp(key,"SIM_SALT2c") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2c) ; }  
      else if ( strcmp(key,"SIM_SALT2mB") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2mB) ; }  
      else if ( strcmp(key,"SIM_SALT2alpha") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2alpha) ; }  
      else if ( strcmp(key,"SIM_SALT2beta") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2beta) ; }  
      else if ( strcmp(key,"SIM_SALT2gammaDM") == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.SIM_SALT2gammaDM) ; }  
    }
    else if ( strcmp(key,"SIM_DELTA") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_DELTA) ; }  

    else if ( strcmp(key,"SIM_STRETCH") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_STRETCH) ; }  

    else if ( strcmp(key,"SIM_DM15") == 0 ) 
      { copy_flt(copyFlag, parVal, &SNDATA.SIM_DM15) ; }  

    else if ( strncmp(key,"SIMSED",6) == 0 ) {
      for(ipar=0; ipar < SNDATA.NPAR_SIMSED; ipar++ ) { 
	sprintf(KEY_TEST, "%s", SNDATA.SIMSED_KEYWORD[ipar]) ;
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.SIMSED_PARVAL[ipar]) ; } 
      }
    }
    else if ( strncmp(key,"SIM_PEAKMAG",11) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_PEAKMAG_%c", FILTERSTRING[ifilt_obs]);
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.SIM_PEAKMAG[ifilt]) ; }
      }
    }
    else if ( strncmp(key,"SIM_TEMPLATEMAG",15) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];    
	sprintf(KEY_TEST, "SIM_TEMPLATEMAG_%c", FILTERSTRING[ifilt_obs]);
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.SIM_TEMPLATEMAG[ifilt]) ; }
      }
    }
    else if ( strncmp(key,"SIM_EXPOSURE",12) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_EXPOSURE_%c", FILTERSTRING[ifilt_obs]);
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.SIM_EXPOSURE_TIME[ifilt]) ; }
      }
    }
    else if ( strncmp(key,"SIM_GALFRAC",11) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_GALFRAC_%c", FILTERSTRING[ifilt_obs]);
	if ( strcmp(key,KEY_TEST) == 0 ) 
	  { copy_flt(copyFlag, parVal, &SNDATA.SIM_GALFRAC[ifilt]) ; }
      }
    }

    else if ( strncmp(key,"SIM_STRONGLENS",14) == 0 ) {
      // continue with SIM_STROGNLENS ...
      sprintf(PREFIX, "SIM_STRONGLENS") ;

      sprintf(KEY_TEST,"%s_ID", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_int(copyFlag, parVal, &SNDATA.SIM_SL_IDLENS); }

      sprintf(KEY_TEST,"%s_z", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_dbl(copyFlag, parVal, &SNDATA.SIM_SL_zLENS); }

      sprintf(KEY_TEST,"%s_TDELAY", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_dbl(copyFlag, parVal, &SNDATA.SIM_SL_TDELAY); }

      sprintf(KEY_TEST,"%s_MAGSHIFT", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_dbl(copyFlag, parVal, &SNDATA.SIM_SL_MAGSHIFT); }

      sprintf(KEY_TEST,"%s_NIMG", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_int(copyFlag, parVal, &SNDATA.SIM_SL_NIMG); }

      sprintf(KEY_TEST,"%s_IMGNUM", PREFIX);
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_int(copyFlag, parVal, &SNDATA.SIM_SL_IMGNUM); }

    }
    else {
      sprintf(c1err,"Unknown SIM key = %s", key);
      sprintf(c2err,"stringVal='%s'  parVal=%f", stringVal, parVal[0] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }

  }

  // other SIM keys that don't start with SIM
  else if ( ncmp_PySEDMODEL==0 && len_PySEDMODEL > 2) {
    for(ipar=0; ipar < SNDATA.NPAR_PySEDMODEL; ipar++ ) { 
      sprintf(KEY_TEST, "%s", SNDATA.PySEDMODEL_KEYWORD[ipar]) ;
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.PySEDMODEL_PARVAL[ipar]) ;  } 
    }
  }

  else if ( strncmp(key,"LCLIB_PARAM",11) == 0 || 
	    strncmp(key,"LCLIB(",6) == 0 ) {    // read legacy PLASTICC data from 2018
    for(ipar=0; ipar < SNDATA.NPAR_LCLIB; ipar++ ) { 
      sprintf(KEY_TEST, "%s", SNDATA.LCLIB_KEYWORD[ipar]) ;
      if ( strcmp(key,KEY_TEST) == 0 ) 
	{ copy_flt(copyFlag, parVal, &SNDATA.LCLIB_PARVAL[ipar]) ; } 
    }
  }
  
  else {
    // error message
    sprintf(c1err,"Unknown key = %s (copyFlag=%d)", key, copyFlag);
    sprintf(c2err,"stringVal='%s'  parVal=%f", stringVal, parVal[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  return ;

}  // end copy_SNDATA_HEAD

// = = = = = = = = = = = = = = = = = = = = = = = = 
int select_MJD_SNDATA(double *CUTWIN_MJD) {

  // Created Feb 11 2021
  // Use input CUTWIN_MJD to select a subset of observations
  // for copy_SNDATA_OBS.  The sparse list OBS_STORE_LIST
  // enables fast copy by only looping over elements to keep.
  // Note that copy_SNDATA_OBS will abort if this function is not 
  // called for each event.
  //
  // Function returns SNDATA.NOBS_STORE
  //
  // Example: for a 10-year survey, a typical MJD window contains
  // 1 of the 10 seasons, and thus copy_SNDATA_OBS returns just 1 
  // season of data instead of all 10 seasons.

  int  o, NOBS_STORE  = 0;
  int  NOBS = SNDATA.NOBS ;
  bool LDMP = false ;
  double MJD;
  char fnam[] = "select_MJD_SNDATA" ;
  // ---------- BEGIN -------

  for(o=1; o <= NOBS; o++ ) {  // note fortran-like index
    MJD = SNDATA.MJD[o];
    if ( MJD < CUTWIN_MJD[0] ) { continue; }
    if ( MJD > CUTWIN_MJD[1] ) { continue; }
    SNDATA.OBS_STORE_LIST[NOBS_STORE] = o ;
    NOBS_STORE++ ;
  }

  if ( LDMP ) {
    printf(" xxx %s: CUTWIN_MJD = %.1f to %.1f  NOBS_STORE=%d\n", 
	   fnam, CUTWIN_MJD[0], CUTWIN_MJD[1], NOBS_STORE);
  }

  SNDATA.NOBS_STORE = NOBS_STORE ;
  return(NOBS_STORE) ;
} // end select_MJD_SNDATA

int select_mjd_sndata__(double *MJD_WINDOW) 
  {  return select_MJD_SNDATA(MJD_WINDOW); }

// = = = = = = = = = = = = = = = = = = = = = = = =
void host_property_list_sndata(char *HOST_PROPERTY_LIST) {
  // Created Mar 14 2022
  // Return list of stored host property strings based on
  // true values, e.g,
  //   'LOGMASS,LOGsSFR'
  // Original intent is to tell analysis codes which properties
  // to store in output tables.

  double NOVAR = HOSTLIB_PROPERTY_UNDEFINED + 1.0;
  char fnam[] = "host_property_list_sndata";
  char TMPLIST[MXPATHLEN];
  HOST_PROPERTY_LIST[0] = TMPLIST[0] = 0 ;

  if ( SNDATA.HOSTGAL_LOGMASS_TRUE[0] > NOVAR ) 
    { catVarList_with_comma(TMPLIST,HOSTGAL_PROPERTY_BASENAME_LOGMASS); }

  if ( SNDATA.HOSTGAL_LOGSFR_TRUE[0] > NOVAR ) 
    { catVarList_with_comma(TMPLIST,HOSTGAL_PROPERTY_BASENAME_LOGSFR); }

  if ( SNDATA.HOSTGAL_LOGsSFR_TRUE[0] > NOVAR ) 
    { catVarList_with_comma(TMPLIST,HOSTGAL_PROPERTY_BASENAME_LOGsSFR); }

  if ( SNDATA.HOSTGAL_COLOR_TRUE[0] > NOVAR ) 
    { catVarList_with_comma(TMPLIST,HOSTGAL_PROPERTY_BASENAME_COLOR); }

  sprintf(HOST_PROPERTY_LIST,"%s", TMPLIST);

  return;
} // end host_property_list_sndata

void host_property_list_sndata__(char *HOST_PROPERTY_LIST) 
{ host_property_list_sndata(HOST_PROPERTY_LIST); }

// = = = = = = = = = = = = = = = = = = = = = = = =
void copy_SNDATA_OBS(int copyFlag, char *key, int NVAL, 
		       char *stringVal, double *parVal ) {

  // Created Feb 2021:
  // For input *key, copy observations to/from SNDATA struct.
  // Note that select_MJD_SNDATA must be called first to 
  // select MJD subset; if not, this function aborts.
  //
  // Inputs:
  //  copyFlag : 
  //     +1 -> copy from string or parVal to SNDATA (prep  data write)
  //     -1 -> copy from SNDATA to string or parVal (after read data)
  //
  //   *key  : name of variable to copy to/from SNDATA struct
  //   NVAL  : number of values to copy
  // 

  int  NOBS       = SNDATA.NOBS ;
  int  NOBS_STORE = SNDATA.NOBS_STORE ;
  int  obs, OBS, NSPLIT ;
  char **str2d ;
  char fnam[] = "copy_SNDATA_OBS" ;

  // ------------- BEGIN ------------

  if(copyFlag<0) { sprintf(stringVal,"NOTSET");  parVal[0] = -999.0; }
  
  if ( NOBS_STORE < 0 ) {
    sprintf(c1err,"Must call select_MJD_SNDATA to set MJD window");
    sprintf(c2err,"for which obs to copy (CID=%s)", SNDATA.CCID );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( strcmp(key,"MJD") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {  // sparse index
      OBS = SNDATA.OBS_STORE_LIST[obs];     // index to full list
      copy_dbl(copyFlag, &parVal[obs], &SNDATA.MJD[OBS]) ; 
    }  
  } 
  else if ( strcmp(key,"FLT") == 0 || strcmp(key,"BAND") == 0 ) {

    if ( copyFlag > 0 ) {
      sprintf(c1err,"key = %s doesn't work with copyFlag=%d", key, copyFlag);
      sprintf(c2err,"Needs a code fix here");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // FILTCHAR_1D includes every observation; here we pick out subset
    // subset of FILTCHAR_1D that are on STORE_LIST

    splitString(SNDATA.FILTCHAR_1D, COMMA, MXEPOCH,    // inputs    
		&NSPLIT, &SNDATA.FILTCHAR[1] );            // outputs 

    stringVal[0] = 0 ;
    for(obs=0; obs < NOBS_STORE; obs++ ) { 
      OBS = SNDATA.OBS_STORE_LIST[obs]; // back to C index    
      catVarList_with_comma(stringVal, SNDATA.FILTCHAR[OBS] );
    }

  }
  else if ( strcmp(key,"FIELD") == 0 ) {
    if ( copyFlag > 0 ) {
      sprintf(c1err,"key = %s doesn't work with copyFlag=%d", key, copyFlag);
      sprintf(c2err,"Needs a code fix here");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    splitString(SNDATA.FIELDNAME_1D, COMMA, MXEPOCH,    // inputs    
		&NSPLIT, &SNDATA.FIELDNAME[1] );            // outputs 

    stringVal[0] = 0 ;
    for(obs=0; obs < NOBS_STORE; obs++ ) { 
      OBS = SNDATA.OBS_STORE_LIST[obs] ;    
      catVarList_with_comma(stringVal, SNDATA.FIELDNAME[OBS] );
    }

  }
  else if ( strcmp(key,"CCDNUM") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_int(copyFlag, &parVal[obs], &SNDATA.CCDNUM[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"IMGNUM") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_int(copyFlag, &parVal[obs], &SNDATA.IMGNUM[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"PHOTFLAG") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_int(copyFlag, &parVal[obs], &SNDATA.PHOTFLAG[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"PHOTPROB") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.PHOTPROB[OBS]) ; 
    }  
  }

  else if ( strcmp(key,"FLUXCAL") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.FLUXCAL[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"FLUXCALERR") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.FLUXCAL_ERRTOT[OBS]) ; 
    }  
  }

  else if ( strcmp(key,"PSF_SIG1") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.PSF_SIG1[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"PSF_SIG2") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.PSF_SIG2[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"PSF_RATIO") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.PSF_RATIO[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"PSF_NEA") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.PSF_NEA[OBS]) ; 
    }  
  }

  else if ( strcmp(key,"SKY_SIG") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.SKY_SIG[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"SKY_SIG_T") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.SKY_SIG_T[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"ZEROPT") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.ZEROPT[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"ZEROPT_ERR") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.ZEROPT_ERR[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"GAIN") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.GAIN[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"XPIX") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.XPIX[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"YPIX") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.YPIX[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"SIM_MAGOBS") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.SIMEPOCH_MAG[OBS]) ; 
    }  
  }
  else if ( strcmp(key,"SIM_FLUXCAL_HOSTERR") == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.SIMEPOCH_FLUXCAL_HOSTERR[OBS]) ; 
    }  
  }

  else if ( strcmp(key,SNDATA.VARNAME_SNRMON) == 0 ) {
    for(obs=0; obs < NOBS_STORE; obs++ ) {
      OBS = SNDATA.OBS_STORE_LIST[obs];  
      copy_flt(copyFlag, &parVal[obs], &SNDATA.SIMEPOCH_SNRMON[OBS]) ; 
    }  
  }

  else {
    // error message
    sprintf(c1err,"Unknown key = %s (copyFlag=%d)", key, copyFlag);
    sprintf(c2err,"stringVal='%s'  parVal=%f", stringVal, parVal[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if  ( NVAL != NOBS_STORE ) {
    sprintf(c1err,"Copied %d values (NOBS_STORE)", NOBS_STORE);
    sprintf(c2err,"but expected NVAL=%d", NVAL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  return ;

} // end copy_SNDATA_OBS


// ==========================================
void copy_GENSPEC(int copyFlag, char *key, int ispec, double *parVal ) {

  // Created Feb 18 2021
  // Return contents of GENSPEC structure; 
  // intended as fortran interface to read spectra from data.

  int  NBLAM, ilam ;
  char fnam[] = "copy_GENSPEC" ;

  // ------------ BEGIN ----------

  if ( copyFlag<0 ) {  parVal[0] = -999.0; }  // init outpuit

  if ( ispec >= 0 ) { NBLAM = GENSPEC.NBLAM_VALID[ispec]; }

  if ( strcmp(key,"NSPECTRA") == 0 ) 
    { copy_int(copyFlag, parVal, &GENSPEC.NMJD_PROC );  }

  // ispec-dependent info

  else if ( strcmp(key,"ID") == 0 ) 
    { copy_int(copyFlag, parVal, &GENSPEC.ID_LIST[ispec] );  }

  else if ( strcmp(key,"MJD") == 0 ) 
    { copy_dbl(copyFlag, parVal, &GENSPEC.MJD_LIST[ispec] );  }

  else if ( strcmp(key,"NBLAM") == 0 ) 
    { copy_int(copyFlag, parVal, &GENSPEC.NBLAM_VALID[ispec] );  }

  else if ( strcmp(key,"TEXPOSE") == 0 ) 
    { copy_dbl(copyFlag, parVal, &GENSPEC.TEXPOSE_LIST[ispec] );  }

  // lam-dependent arrays

  else if ( strcmp(key,"LAMMIN") == 0 ) {
    for(ilam=0; ilam  < NBLAM; ilam++ ) 
      { copy_dbl(copyFlag, &parVal[ilam], &GENSPEC.LAMMIN_LIST[ispec][ilam] ); }
  }
  else if ( strcmp(key,"LAMMAX") == 0 ) {
    for(ilam=0; ilam  < NBLAM; ilam++ ) 
      { copy_dbl(copyFlag, &parVal[ilam], &GENSPEC.LAMMAX_LIST[ispec][ilam] ); }
  }
  else if ( strcmp(key,"FLAM") == 0 ) {
    for(ilam=0; ilam  < NBLAM; ilam++ ) 
      { copy_dbl(copyFlag, &parVal[ilam], &GENSPEC.FLAM_LIST[ispec][ilam] ); }
  }
  else if ( strcmp(key,"FLAMERR") == 0 ) {
    for(ilam=0; ilam  < NBLAM; ilam++ ) 
      { copy_dbl(copyFlag, &parVal[ilam], &GENSPEC.FLAMERR_LIST[ispec][ilam]); }

  }
  else if ( strcmp(key,"SIM_FLAM") == 0 ) {
    for(ilam=0; ilam  < NBLAM; ilam++ ) 
      { copy_dbl(copyFlag, &parVal[ilam], &GENSPEC.GENFLAM_LIST[ispec][ilam]); }

  }

  return;

}  // end copy_GENSPEC

// ==============================================
void RD_OVERRIDE_INIT(char *OVERRIDE_FILE) {

  // read and store columns from comma-sep list of  override files 
  // to override values in data headers (not photometry).
  // Allows float/double/int, but not strings
  // (e.g., cannot override SNID, FIELD, .. )

  int NROW, ivar, ifile, NFILE = 0;
  int OPTMASK_SNTABLE = 4; // append next file
  char **file_list, *ptrFile;
  char TABLE_NAME[] = "OVERRIDE";
  char VARLIST[]    = "ALL";
  char fnam[]       = "RD_OVERRIDE_INIT";

  // ----------- BEGIN -----------

  RD_OVERRIDE.USE = false;
  if ( IGNOREFILE(OVERRIDE_FILE) ) { return; }

  print_banner(fnam);

  // split comma-sep OVERRIDE_FILE 
  parse_commaSepList(fnam, OVERRIDE_FILE, MXFILE_OVERRIDE, MXPATHLEN,
		     &NFILE, &file_list ); // <== returned
  
  SNTABLE_AUTOSTORE_RESET(); // May 2022
  for(ifile=0; ifile < NFILE; ifile++ ) {
    ptrFile = file_list[ifile] ;
    ENVreplace(ptrFile, fnam, 1);
    NROW = SNTABLE_AUTOSTORE_INIT(ptrFile, TABLE_NAME, VARLIST,
				  OPTMASK_SNTABLE );
    printf("   Stored %d rows of header-override data\n", NROW);

    if ( NROW == 0 ) {
      sprintf(c1err,"NROW=0 in HEADER_OVERRIDE_FILE");
      sprintf(c2err,"Check %s", ptrFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  } // end ifile

  // - - - - - - - - - - - - 
  RD_OVERRIDE.USE    = true ;
  RD_OVERRIDE.NFILE  = NFILE;
  for(ivar=0; ivar < MXVAR_OVERRIDE; ivar++ )
    { RD_OVERRIDE.N_PER_VAR[ivar] = 0 ; }

  // - - - - - - - 
  // set z logicals in case zHEL <-> zCMB needs to be recomputed
  RD_OVERRIDE.FOUND_zCMB = false;
  RD_OVERRIDE.FOUND_zHEL = false;
  if ( EXIST_VARNAME_AUTOSTORE("REDSHIFT_FINAL") ) 
    { RD_OVERRIDE.FOUND_zCMB = true; }
  if ( EXIST_VARNAME_AUTOSTORE("REDSHIFT_CMB") ) 
    { RD_OVERRIDE.FOUND_zCMB = true; }
  if ( EXIST_VARNAME_AUTOSTORE("REDSHIFT_HELIO") ) 
    { RD_OVERRIDE.FOUND_zHEL = true; }

  // if both zCMB and zHEL are on header-override list,
  // turn them off since there is no need to recompute.
  if ( RD_OVERRIDE.FOUND_zCMB && RD_OVERRIDE.FOUND_zHEL ) 
    { RD_OVERRIDE.FOUND_zCMB = RD_OVERRIDE.FOUND_zHEL = false; }

  return ;


} // end RD_OVERRIDE_INIT

int RD_OVERRIDE_FETCH(char *CCID, char *VARNAME, double *DVAL) {

  // Created Dec 2021
  // If SNID and VARNAME is on override list, return DVAL
  // and function returns 1.
  // Function returns 0 if there is no override.

  int  ISTAT, NRD, IVAR, NTMP ;
  char STRDUM[20];
  char fnam[] = "RD_OVERRIDE_FETCH";

  // ----------- BEGIN -----------
  *DVAL = 0.0;
  if ( !RD_OVERRIDE.USE ) { return 0; }

  IVAR = IVAR_VARNAME_AUTOSTORE(VARNAME);
  if ( IVAR < 0 ) { return 0; }

  // read from override table; ISTAT and DVALare returned
  SNTABLE_AUTOSTORE_READ(CCID, VARNAME, &ISTAT, DVAL, STRDUM);  

  // *ISTAT =  0  if CCID is found; 
  // *ISTAT = -1  if CCID is NOT found    
  // *ISTAT = -2  if VARNAME is NOT found

  if ( ISTAT == 0 ) {
    //    printf(" xxx %s: CID=%s  '%s'=%f  ISTAT=%d \n",
    //	   fnam, CCID, VARNAME, *DVAL, ISTAT); fflush(stdout);
    NRD = 1;
    if ( RD_OVERRIDE.N_PER_VAR[IVAR] == 0 ) 
      { printf("\t Found override for %s\n", VARNAME );  fflush(stdout); }

    RD_OVERRIDE.N_PER_VAR[IVAR]++ ;
  }
  else {
    NRD = 0 ;
  }

  return NRD ;

} // end RD_OVERRIDE_FETCH

// =====================================
void RD_OVERRIDE_POSTPROC(void) {

  char fnam[] = "RD_OVERRIDE_POSTPROC" ;
  
  // ------------ BEGIN --------------

  if ( !RD_OVERRIDE.USE ) { return; }

  // for text format, check for variables that are not in
  // the data files and thus header_override is an addition.
  if ( FORMAT_SNDATA_READ == FORMAT_SNDATA_TEXT ) 
    { rd_override_append(); }

  // check for redshift_helio update that forces zcmb to also change.
  rd_override_zcalc();

  return ;

} // end RD_OVERRIDE_POSTPROC


void rd_override_append(void) {

  // Called only for TEXT format that might be missing
  // some variables.

#define NVAR_OVERRIDE_CHECK 16

  char VARNAME_CHECK[NVAR_OVERRIDE_CHECK][40] = {
    "VPEC", "VPEC_ERR", 
    "REDSHIFT_HELIO",  "REDSHIFT_HELIO_ERR",
    "REDSHIFT_FINAL",  "REDSHIFT_FINAL_ERR",
    "REDSHIFT_CMB",    "REDSHIFT_CMB_ERR",
    "HOSTGAL_LOGMASS", "HOSTGAL_LOGMASS_ERR", 
    "HOSTGAL_LOGSFR",  "HOSTGAL_LOGSFR_ERR", 
    "HOSTGAL_LOGsSFR", "HOSTGAL_LOGsSFR_ERR", 
    "HOSTGAL_COLOR",   "HOSTGAL_COLOR_ERR"
  };
  float *ptr_SNDATA[NVAR_OVERRIDE_CHECK] = {
    &SNDATA.VPEC,            &SNDATA.VPEC_ERR, 
    &SNDATA.REDSHIFT_HELIO,  &SNDATA.REDSHIFT_HELIO_ERR,
    &SNDATA.REDSHIFT_FINAL,  &SNDATA.REDSHIFT_FINAL_ERR,
    &SNDATA.REDSHIFT_FINAL,  &SNDATA.REDSHIFT_FINAL_ERR,
    &SNDATA.HOSTGAL_LOGMASS_OBS[0], &SNDATA.HOSTGAL_LOGMASS_ERR[0],
    &SNDATA.HOSTGAL_LOGSFR_OBS[0],  &SNDATA.HOSTGAL_LOGSFR_ERR[0],
    &SNDATA.HOSTGAL_LOGsSFR_OBS[0], &SNDATA.HOSTGAL_LOGsSFR_ERR[0],
    &SNDATA.HOSTGAL_COLOR_OBS[0],   &SNDATA.HOSTGAL_COLOR_ERR[0]
  } ;

  int ivar;
  double DVAL ;
  char *varName ;
  char fnam[] = "rd_override_append" ;

  // --------- BEGIN --------

  for (ivar=0; ivar < NVAR_OVERRIDE_CHECK; ivar++ ) {
    varName = VARNAME_CHECK[ivar] ;
    if ( EXIST_VARNAME_AUTOSTORE(varName) ) { 
      RD_OVERRIDE_FETCH(SNDATA.CCID, varName, &DVAL) ;
      *ptr_SNDATA[ivar] = (float)DVAL;

      if ( strstr(varName,"HOSTGAL") != NULL ) {
	SNDATA.HOSTGAL_NMATCH[0] = 1;
      }

    }  // end exist varname in header overrid
  } // end ivar loop

  return;
} // end rd_override_append

void rd_override_zcalc(void) {

  // If either zCMB or zHEL is on override list; recompute the other.
  
  double RA, DEC, zCMB, zHEL ;
  bool FOUND_z= ( RD_OVERRIDE.FOUND_zCMB || RD_OVERRIDE.FOUND_zHEL );
  char eq[]   = "eq" ;
  char fnam[] = "rd_override_zcalc" ;
  // ---------- BEGIN -------------

  if ( !FOUND_z ) { return; }
  RA  = SNDATA.RA;  
  DEC = SNDATA.DEC ;

  if ( RD_OVERRIDE.FOUND_zCMB ) {
    zCMB = (double)SNDATA.REDSHIFT_FINAL;
    zHEL = zhelio_zcmb_translator(zCMB,RA,DEC,eq,-1); 
    SNDATA.REDSHIFT_HELIO = (float)zHEL ;
  }
  else if ( RD_OVERRIDE.FOUND_zHEL ) {
    zHEL = (double)SNDATA.REDSHIFT_HELIO ;
    zCMB = zhelio_zcmb_translator(zHEL,RA,DEC,eq,+1);
    SNDATA.REDSHIFT_FINAL = (float)zCMB ;

    /* xxx mark delete xxxxxxxx
    if ( strcmp(SNDATA.CCID,"1032") == 0 ) {
      printf(" xxx %s: ov zHEL=%.4f -> zCMB=%.4f (RA,DEC=%f,%f)\n",
	     fnam, zHEL, zCMB, RA, DEC); fflush(stdout);
    }
    xxxxxx end mark xxx */
  }

  return ;

} // end rd_override_zcalc


// - - - - - - - 
// mangled fortran functions

void copy_sndata_global__(int *copyFlag, char *key, int *NVAL, 
			  char *stringVal, double *parVal ) 
{ copy_SNDATA_GLOBAL(*copyFlag, key, *NVAL, stringVal, parVal); }

void copy_sndata_head__(int *copyFlag, char *key, int *NVAL, 
			char *stringVal, double *parVal ) 
{ copy_SNDATA_HEAD(*copyFlag, key, *NVAL, stringVal, parVal); }

void copy_sndata_obs__(int *copyFlag, char *key, int *NVAL, 
		       char *stringVal, double *parVal ) 
{ copy_SNDATA_OBS(*copyFlag, key, *NVAL, stringVal, parVal); }

void copy_genspec__(int *copyFlag, char *key, int *ispec, double *parVal ) 
{ copy_GENSPEC(*copyFlag, key, *ispec, parVal); }

void rd_override_init__(char *override_file)
{ RD_OVERRIDE_INIT(override_file); }
