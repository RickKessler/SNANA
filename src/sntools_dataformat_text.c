/************************************************
  Created Feb 2021

  Refactor write/read for text format to allow extracting
  text files from FITS format. Longer term, data will be
  written directly to FITS format because we can't create 
  millions of intermediate text files ... but sometimes 
  it's useful to extract a few events in text format for 
  debugging. 

  Another goal is to cleanup so really old code that has
  so much legacy crap that it's hard to follow.

  Analogous to sntools_dataformat_fits.c[h], here we use
  the SNDATA structure to pass the data.

*************************************************/

#include  "sntools.h"
#include  "sntools_dataformat_text.h"
#include  "sntools_host.h" 
#include  "sntools_trigger.h" 
#include  "sntools_spectrograph.h"


void WR_SNTEXTIO_DATAFILE(char *OUTFILE) {

  // Createed Feb 2021
  // Driver function to write single data file for single event.
  //

  FILE *fp ;
  char fnam[] = "WR_SNTEXTIO_DATAFILE" ;

  // ----------- BEGIN ------------

  fp = fopen(OUTFILE, "wt");
  if ( !fp ) {
    sprintf(c1err,"Could not open data text file");
    sprintf(c2err,"%s", OUTFILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  wr_dataformat_text_HEADER(fp);

  wr_dataformat_text_SNPHOT(fp);

  wr_dataformat_text_SNSPEC(fp);
    
  fclose(fp);

  return ;

} // end WR_SNTEXTIO_DATAFILE

void wr_sntextio_datafile__(char *OUTFILE)  
{ WR_SNTEXTIO_DATAFILE(OUTFILE); }

// =====================================================
void  wr_dataformat_text_HEADER(FILE *fp) {

  char comment[80];

  char COMMENT_FAKEFLAG[4][80] = {
    "data",
    "fakes overlaid on images",
    "simulated with snlc_sim.exe" ,
    "BLIND-TEST simulation"
  } ;

  char SURVEY_ARG[100];
  char fnam[] = "wr_dataformat_text_HEADER" ;

  // ------------ BEGIN -----------

  // write either "SURVEY: SURVEY" or "SURVEY: SURVEY(SUBSAMPLE)"
  int LENS = strlen(SNDATA.SUBSURVEY_NAME);
  int OVP  = strcmp(SNDATA.SURVEY_NAME,SNDATA.SUBSURVEY_NAME);
  if ( LENS > 0 && OVP!=0 ) 
    { sprintf(SURVEY_ARG, "%s(%s)", 
	      SNDATA.SURVEY_NAME, SNDATA.SUBSURVEY_NAME ); }
  else
    { sprintf(SURVEY_ARG, "%s", SNDATA.SURVEY_NAME); }
  fprintf(fp,"SURVEY:   %s\n", SURVEY_ARG);


  fprintf(fp,"SNID:     %s\n", SNDATA.CCID);
  fprintf(fp,"IAUC:     %s\n", SNDATA.IAUC_NAME);
  // fprintf(fp,"SNTYPE:   %d\n", SNDATA.SEARCH_TYPE);
  fprintf(fp,"SNTYPE:   %d\n", SNDATA.SNTYPE);
  fprintf(fp,"RA:       %.6f  # deg\n", SNDATA.RA);
  fprintf(fp,"DEC:      %.6f  # deg\n", SNDATA.DEC);
  fprintf(fp,"FILTERS:  %s\n", SNDATA_FILTER.LIST);
  fprintf(fp,"PIXSIZE:  %.4f  # arcsec \n", SNDATA.PIXSIZE);
  fprintf(fp,"CCDNUM:   %d   \n", SNDATA.CCDNUM[0] );

  int FAKE = SNDATA.FAKE ;
  if ( FAKE < 0 ) {
    sprintf(c1err,"Invalid FAKE = %d (check FAKE key in header)", FAKE);
    sprintf(c2err,"For valid values, grep FAKEFLAG sndata.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  fprintf(fp, "FAKE:     %d  # %s \n", FAKE, COMMENT_FAKEFLAG[FAKE] );

  fprintf(fp, "MWEBV:    %6.4f +- %.4f  # MW E(B-V) \n", 
	  SNDATA.MWEBV, SNDATA.MWEBV_ERR);
  if ( SNDATA.APPLYFLAG_MWEBV ) {
    fprintf(fp, "MWEBV_APPLYFLAG: %d    # FLUXCAL corrected for MWEBV \n",
	    SNDATA.APPLYFLAG_MWEBV ) ;
  }

  fprintf(fp, "PEAKMJD:  %9.3f     # estimate for LC fit codes\n", 
	  SNDATA.SEARCH_PEAKMJD ); 

  // redshift info
  fprintf(fp, "\n");

  fprintf(fp,"REDSHIFT_HELIO:   %.6f +- %.6f  # Helio \n"
	  ,SNDATA.REDSHIFT_HELIO, SNDATA.REDSHIFT_HELIO_ERR );

  fprintf(fp,"REDSHIFT_FINAL:   %.6f +- %.6f  # CMB (no VPEC corr) \n"
	  ,SNDATA.REDSHIFT_FINAL, SNDATA.REDSHIFT_FINAL_ERR );

  fprintf(fp,"VPEC:             %8.2f +- %.2f      # v_pec correction\n", 
	  SNDATA.VPEC, SNDATA.VPEC_ERR );
  
  // HOST galaxy info
  fprintf(fp, "\n");
  wr_dataformat_text_HOSTGAL(fp);

  // --------------------------------------------------
  // write optional AUXHEADER_LINES (June 2012)
  // at end of header, but before the SIM_XXX keys
  int j ;
  for(j=0; j < SNDATA.NLINES_AUXHEADER; j++ ) 
    {  fprintf(fp,"%s \n", SNDATA.AUXHEADER_LINES[j]); }


  // - - - - - - - - -  -
  // SIM_ info
  if ( SNDATA.WRFLAG_BLINDTEST      ) { return; } // skip for BLIND test
  if ( SNDATA.FAKE == FAKEFLAG_DATA ) { return; } // skip for real data

  wr_dataformat_text_SIMPAR(fp);

  return ;
} // end wr_dataformat_text_HEADER

// ========================================= 
void wr_dataformat_text_SIMPAR(FILE *fp) {

  // write sim-truth variables to data file header (text format)

  int NTMP, NPAR, ipar, ifilt, ifilt_obs, iscat ;
  char key[80], ctmp[100] ;  
  float parval ;
  char fnam[] = "wr_dataformat_text_SIMPAR" ;

  // ------------ BEGIN -----------
  // write SIM-truth params

  fprintf(fp, "\n" );
  fprintf(fp, "# SIM-truth values\n" );
  
  fprintf(fp, "SIM_MODEL_NAME:      %s  \n",
	  SNDATA.SIM_MODEL_NAME  ) ;

  fprintf(fp, "SIM_MODEL_INDEX:     %d  # grep MODEL_ sntools.h | grep def\n",
	  SNDATA.SIM_MODEL_INDEX  ) ;

  fprintf(fp, "SIM_TYPE_INDEX:      %d  # true GENTYPE \n",
	  SNDATA.SIM_TYPE_INDEX );

  fprintf(fp, "SIM_TYPE_NAME:       %s \n",
	  SNDATA.SIM_TYPE_NAME );

  fprintf(fp, "SIM_TEMPLATE_INDEX:  %d   # template index for SIMSED,NONIa\n", 
	  SNDATA.SIM_TEMPLATE_INDEX ) ;

  // Jun 2017: SUBSAMPLE_INDEX                                                
  if ( SNDATA.SUBSAMPLE_INDEX >= 0 ) {
    fprintf(fp,"SIM_SUBSAMPLE_INDEX: %d \n", SNDATA.SUBSAMPLE_INDEX);
  }
  //  fprintf(fp, "SIM_COMMENT:  %s  \n", SNDATA.SIM_COMMENT  ) ;

  fprintf(fp, "SIM_LIBID:           %d    # LIBID index in SIMLIB file\n",
	  SNDATA.SIM_LIBID ) ;

  fprintf(fp, "SIM_NGEN_LIBID:      %d  \n", 
	  SNDATA.SIM_NGEN_LIBID ) ;

  fprintf(fp, "SIMLIB_MSKOPT:       %d  \n", 
	  SNDATA.SIMLIB_MSKOPT ) ;

  fprintf(fp, "SIM_NOBS_UNDEFINED:  %d  \n", 
	  SNDATA.SIM_NOBS_UNDEFINED );

  fprintf(fp, "SIM_REDSHIFT_HELIO:  %.5f  \n",
	  SNDATA.SIM_REDSHIFT_HELIO );

  fprintf(fp, "SIM_REDSHIFT_CMB:    %.5f  \n",
	  SNDATA.SIM_REDSHIFT_CMB );

  fprintf(fp, "SIM_REDSHIFT_HOST:   %.5f  \n",
	  SNDATA.SIM_REDSHIFT_HOST );

  fprintf(fp,"SIM_REDSHIFT_FLAG:   %d  # %s\n", 
	  SNDATA.SIM_REDSHIFT_FLAG, SNDATA.SIM_REDSHIFT_COMMENT );

  fprintf(fp, "SIM_VPEC:            %.1f   # km/sec) \n", 
	  SNDATA.SIM_VPEC ) ;

  // - - - -  SIM_HOSTLIB_

  fprintf(fp,"SIM_HOSTLIB_GALID:   %lld \n", SNDATA.SIM_HOSTLIB_GALID);
  NPAR = SNDATA.NPAR_SIM_HOSTLIB;
  fprintf(fp, "SIM_HOSTLIB_NPAR:    %d \n", NPAR);
  for(ipar=0; ipar < NPAR; ipar++ ) {
    sprintf(key,"%s:", SNDATA.SIM_HOSTLIB_KEYWORD[ipar] );
    parval = SNDATA.SIM_HOSTLIB_PARVAL[ipar] ;
    fprintf(fp, "%-28.28s  %.3f \n", key, parval);
  }

  fprintf(fp, "SIM_DLMU:            %.4f   # mag   [ -5*log10(10pc/dL) ]\n", 
	  SNDATA.SIM_DLMU );

  fprintf(fp, "SIM_LENSDMU:         %.4f   # mag \n", 
	  SNDATA.SIM_LENSDMU );

  fprintf(fp, "SIM_RA:              %f     # deg  \n", 
	  SNDATA.SIM_RA );
  fprintf(fp, "SIM_DEC:             %f      # deg  \n", 
	  SNDATA.SIM_DEC );

  fprintf(fp, "SIM_MWRV:            %.3f    # MilkyWay RV \n", 
	  SNDATA.SIM_MWRV );

  fprintf(fp, "SIM_MWEBV:           %.4f   # MilkyWay E(B-V) \n", 
	  SNDATA.SIM_MWEBV );

  fprintf(fp, "SIMOPT_MWCOLORLAW:   %d   # MW color-law option \n", 
	  SNDATA.SIMOPT_MWCOLORLAW ) ;

  fprintf(fp, "SIMOPT_MWEBV:        %d   # MWEBV option \n",
	  SNDATA.SIMOPT_MWEBV ) ;

  if ( SNDATA.SIM_AV != NULLFLOAT ) {

    fprintf(fp, "SIM_AVTAU:         %.3f     # dN/dAV = exp(-AV/AVTAU) \n",
	    SNDATA.SIM_AVTAU ) ;

    fprintf(fp, "SIM_AV:            %.3f     # host extinction at 5510 A\n",
	    SNDATA.SIM_AV ) ;

    fprintf(fp, "SIM_RV:            %.3f     # CCM89 extinction param \n",
	    SNDATA.SIM_RV ) ;
  }

  fprintf(fp, "SIM_MAGSMEAR_COH:    %.3f     # coh model scatter \n", 
	  SNDATA.SIM_MAGSMEAR_COH ) ;

  fprintf(fp, "SIM_PEAKMJD:         %.3f    #  true PEAKMJD, days \n", 
	  SNDATA.SIM_PEAKMJD ) ;


  if ( SNDATA.SIM_DELTA != NULLFLOAT ) {
    fprintf(fp, "SIM_DELTA:         %.3f    # MLCS lumi-par \n", 
	    SNDATA.SIM_DELTA) ;
  }
  if ( SNDATA.SIM_DM15 != NULLFLOAT ) {
    fprintf(fp, "SIM_DM15:          %.3f    # DM15 lumi-par \n",
	    SNDATA.SIM_DM15 ) ;
  }

  if ( SNDATA.SIM_SALT2alpha != NULLFLOAT ) {
    fprintf(fp, "SIM_SALT2alpha:      %.3f \n", SNDATA.SIM_SALT2alpha ) ;
    fprintf(fp, "SIM_SALT2beta:       %.3f \n", SNDATA.SIM_SALT2beta ) ;
    fprintf(fp, "SIM_SALT2gammaDM:    %.3f \n", SNDATA.SIM_SALT2gammaDM ) ;
  }

  if ( SNDATA.SIM_SALT2x0 != NULLFLOAT ) {
    fprintf(fp, "SIM_SALT2x0:         %.4e \n", SNDATA.SIM_SALT2x0 ) ;
  }
  if ( SNDATA.SIM_SALT2x1 != NULLFLOAT ) {
    fprintf(fp, "SIM_SALT2x1:         %.4f \n", SNDATA.SIM_SALT2x1 ) ;
  }
  if ( SNDATA.SIM_SALT2c != NULLFLOAT ) {
    fprintf(fp, "SIM_SALT2c:          %.4f \n", SNDATA.SIM_SALT2c ) ;
  }
  if ( SNDATA.SIM_SALT2mB != NULLFLOAT ) {
    fprintf(fp, "SIM_SALT2mB:         %.4f \n", SNDATA.SIM_SALT2mB ) ;
  }

  if ( SNDATA.SIM_STRETCH > 0 ) {
    fprintf(fp, "SIM_STRETCH:         %.3f\n", SNDATA.SIM_STRETCH ) ;
  }

  // covmat-scatter info                                                      
  if ( SNDATA.SIMFLAG_COVMAT_SCATTER ) {
    for ( iscat=0; iscat < 3; iscat++ ) {
      fprintf(fp, "SIM_SCATTER[%s]:     %.4f \n"
	      , SNDATA.SIM_COVMAT_SCATTER_NAME[iscat]
	      , SNDATA.SIM_COVMAT_SCATTER[iscat] ) ;
    }
  }

  sprintf(ctmp,"bits 1,2=> found by software,spec" );
  fprintf(fp, "SIM_SEARCHEFF_MASK:  %d  # %s \n", 
	  SNDATA.SIM_SEARCHEFF_MASK, ctmp );

  sprintf(ctmp,"spectro-search efficiency");
  fprintf(fp, "SIM_SEARCHEFF_SPEC:  %6.4f  # %s \n", 
	  SNDATA.SIM_SEARCHEFF_SPEC, ctmp );

  if ( SNDATA.SIM_TRESTMIN != NULLFLOAT ) {
    fprintf(fp, "SIM_TRESTMIN:        %.2f  # days \n", SNDATA.SIM_TRESTMIN );
    fprintf(fp, "SIM_TRESTMAX:        %.2f  # days \n", SNDATA.SIM_TRESTMAX );
  }

  fprintf(fp, "SIM_RISETIME_SHIFT:  %.1f   # days \n", 
	  SNDATA.SIM_RISETIME_SHIFT ) ;
  fprintf(fp, "SIM_FALLTIME_SHIFT:  %.1f   # days \n", 
	  SNDATA.SIM_FALLTIME_SHIFT ) ;

  // SIMSED info                                                              
  if ( SNDATA.NPAR_SIMSED > 0 ) {
    fprintf(fp,"\n");
    fprintf(fp,"SIMSED_NPAR:   %d \n", SNDATA.NPAR_SIMSED );
    for ( ipar = 0; ipar < SNDATA.NPAR_SIMSED; ipar++ ) {
      fprintf(fp,"%s:    %f\n"
	      ,SNDATA.SIMSED_KEYWORD[ipar]
	      ,SNDATA.SIMSED_PARVAL[ipar] );
    }
  }

  // PySEDMODEL info for BYOSED, SNEMO                                        
  if ( SNDATA.NPAR_PySEDMODEL > 0 ) {
    fprintf(fp,"\n");
    fprintf(fp,"%s_NPAR:   %d \n",
	    SNDATA.SIM_MODEL_NAME, SNDATA.NPAR_PySEDMODEL );
    for ( ipar = 0; ipar < SNDATA.NPAR_PySEDMODEL; ipar++ ) {
      fprintf(fp,"%s:    %f \n"
	      ,SNDATA.PySEDMODEL_KEYWORD[ipar]
	      ,SNDATA.PySEDMODEL_PARVAL[ipar] );
    }
  }

  // LCLIB info (Sep 8 2017)                                                  
  if ( SNDATA.NPAR_LCLIB > 0 ) {
    fprintf(fp,"\n");
    fprintf(fp,"LCLIB_NPAR:   %d \n", SNDATA.NPAR_LCLIB );
    for ( ipar = 0; ipar < SNDATA.NPAR_LCLIB; ipar++ ) {      
      fprintf(fp,"%s:    %f \n"
	      ,SNDATA.LCLIB_KEYWORD[ipar]
	      ,SNDATA.LCLIB_PARVAL[ipar] );
    }
  }

  // strong lens info (July 20 2019)                                          
  if ( SNDATA.SIM_SL_FLAG ) {
    fprintf(fp,"\n");
    fprintf(fp,"SIM_STRONGLENS_ID:        %d   \n", SNDATA.SIM_SL_IDLENS );
    fprintf(fp,"SIM_STRONGLENS_z:         %.3f \n", SNDATA.SIM_SL_zLENS  );
    fprintf(fp,"SIM_STRONGLENS_DELAY:     %.3f  # days \n",
	    SNDATA.SIM_SL_TDELAY );
    fprintf(fp,"SIM_STRONGLENS_XIMG:      %.3f  # arcsec\n",
	    SNDATA.SIM_SL_XIMG );
    fprintf(fp,"SIM_STRONGLENS_YIMG:      %.3f  # arcsec\n",
	    SNDATA.SIM_SL_YIMG );
    fprintf(fp,"SIM_STRONGLENS_MAGSHIFT:  %.3f \n", SNDATA.SIM_SL_MAGSHIFT );
    fprintf(fp,"SIM_STRONGLENS_NIMG:      %d   \n", SNDATA.SIM_SL_NIMG    );
    fprintf(fp,"SIM_STRONGLENS_IMGNUM:    %d   \n", SNDATA.SIM_SL_IMGNUM  );
  }

  // - - - - - 
  // filter-dependent quantities
  
  fprintf(fp,"\n");
  fprintf(fp,"# %s filter-dependent quantities vs. \n", SNDATA_FILTER.LIST);

  // gal/SN flux-fraction                                                     
  fprintf(fp, "SIM_GALFRAC: "); NTMP = 0;
  for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
    ifilt_obs = SNDATA_FILTER.MAP[ifilt];
    fprintf(fp," %6.3f", SNDATA.SIM_GALFRAC[ifilt_obs] ) ;
    NTMP++ ;
    if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
  }
  fprintf(fp,"  # F_gal/F_SNpeak for PSF=1''\n");

  fprintf(fp, "SIM_PEAKMAG: "); NTMP=0;
  for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
    ifilt_obs = SNDATA_FILTER.MAP[ifilt];
    fprintf(fp," %6.2f",SNDATA.SIM_PEAKMAG[ifilt_obs] ) ;
    NTMP++ ;
    if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
  }
  fprintf(fp,"  # \n");

  if ( SNDATA.SIM_MODEL_INDEX == MODEL_LCLIB ) {
    fprintf(fp, "SIM_TEMPLATEMAG: "); NTMP=0;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.2f",SNDATA.SIM_TEMPLATEMAG[ifilt_obs] ) ;
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp,"  # \n");
  }

  fprintf(fp, "SIM_EXPOSURE: ");  NTMP=0;
  for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
    ifilt_obs = SNDATA_FILTER.MAP[ifilt];
    fprintf(fp," %6.1f",SNDATA.SIM_EXPOSURE_TIME[ifilt_obs] ) ;
    NTMP++ ;
    if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
  }
  fprintf(fp,"  # \n");

  return ;

} // end wr_dataformat_text_SIMPAR


// ===================================
void wr_dataformat_text_HOSTGAL(FILE *fp) {

  // Created Feb 6 2021
  // Copy from wr_HOSTGAL() in sntools.c so that wr_HOSTGAL can be removed.

  int ifilt, ifilt_obs, NTMP, igal, NGAL ;
  char PREFIX[20] = "HOSTGAL";
  char filtlist[MXFILTINDX], ctmp[100] ;
 
  char fnam[] = "wr_dataformat_text_HOSTGAL" ;

  // ------------- BEGIN -------------

  sprintf(filtlist,"%s", SNDATA_FILTER.LIST );

  NGAL = SNDATA.HOSTGAL_NMATCH[0];
  if ( NGAL > MXHOSTGAL ) { NGAL = MXHOSTGAL ; }

  fprintf(fp, "%s_NMATCH:      %d  \n",  
	  PREFIX, SNDATA.HOSTGAL_NMATCH[0] );
  fprintf(fp, "%s_NMATCH2:     %d  \n",  
	  PREFIX, SNDATA.HOSTGAL_NMATCH[1] );

  for(igal=0; igal < NGAL; igal++ ) {

    if ( igal > 0 ) { sprintf(PREFIX,"HOSTGAL%d", igal+1); }

    fprintf(fp, "%s_OBJID:       %lld  \n",  
	    PREFIX, SNDATA.HOSTGAL_OBJID[igal] );

    fprintf(fp, "%s_PHOTOZ:      %.4f  +- %.4f \n", PREFIX,
	    SNDATA.HOSTGAL_PHOTOZ[igal], 
	    SNDATA.HOSTGAL_PHOTOZ_ERR[igal]);

    fprintf(fp, "%s_SPECZ:       %.4f  +- %.4f \n", PREFIX,
	    SNDATA.HOSTGAL_SPECZ[igal], SNDATA.HOSTGAL_SPECZ_ERR[igal] ); 

    fprintf(fp, "%s_RA:          %.6f    # deg \n", 
	    PREFIX, SNDATA.HOSTGAL_RA[igal] );
    fprintf(fp, "%s_DEC:         %.6f    # deg \n", 
	    PREFIX, SNDATA.HOSTGAL_DEC[igal] );

    fprintf(fp, "%s_SNSEP:       %.3f       # arcsec \n", 
	    PREFIX, SNDATA.HOSTGAL_SNSEP[igal] );
    fprintf(fp, "%s_DDLR:        %.3f       # SNSEP/DLR  \n", 
	    PREFIX, SNDATA.HOSTGAL_DDLR[igal] );
    
    if ( igal==0 ) {
      fprintf(fp, "HOSTGAL_CONFUSION:  %6.3f  \n", 
	      SNDATA.HOSTGAL_CONFUSION );
    }

    if ( SNDATA.HOSTGAL_LOGMASS_OBS[igal] > 0.0 ) {
      fprintf(fp, "%s_LOGMASS:     %.3f +- %.3f   # log10(Mgal/Msolar)\n", 
	      PREFIX, 
	      SNDATA.HOSTGAL_LOGMASS_OBS[igal], 
	      SNDATA.HOSTGAL_LOGMASS_ERR[igal] );

      fprintf(fp, "%s_sSFR:        %.3e +- %.3e  # specific SFR\n",
	      PREFIX, 
	      SNDATA.HOSTGAL_sSFR[igal], 
	      SNDATA.HOSTGAL_sSFR_ERR[igal] );
    }

    // if MAGOBS has been read for any filter, then write MAGOBS
    // for all filters.

    if ( SNDATA.HOSTLIB_NFILT_MAGOBS > 0 ) {
      fprintf(fp, "%s_MAG:    ", PREFIX ); NTMP=0;    
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt] ;
	fprintf(fp,"%6.2f ", SNDATA.HOSTGAL_MAG[igal][ifilt] );
	NTMP++ ;
	if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
      }
      fprintf(fp,"# %s\n", filtlist) ;
    }

    if ( SNDATA.HOSTLIB_NFILT_MAGOBS > 0 ) {
      fprintf(fp, "%s_MAGERR: ", PREFIX ); NTMP=0;    
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs = SNDATA_FILTER.MAP[ifilt] ;
	fprintf(fp,"%6.2f ", SNDATA.HOSTGAL_MAGERR[igal][ifilt] );
	NTMP++ ;
	if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
      }
      fprintf(fp,"# %s\n", filtlist) ;
    }
    
    fprintf(fp,"\n");

  } // end igal loop

  // ---------- surface brightness -----------

  sprintf(ctmp,"%s/asec^2",filtlist );
  if ( (SNDATA.HOSTGAL_USEMASK & 4) > 0 ) {
    fprintf(fp,"HOSTGAL_SB_FLUXCAL:  " ); NTMP=0 ;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %.2f",SNDATA.HOSTGAL_SB_FLUXCAL[ifilt] ) ;
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp,"  # %s\n", ctmp );
  }

  // and the uncertainties ...
  if ( (SNDATA.HOSTGAL_USEMASK & 8) > 0 ) {
    fprintf(fp,"HOSTGAL_SB_FLUXCAL_ERR:    " ); NTMP=0 ;
    for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
      ifilt_obs = SNDATA_FILTER.MAP[ifilt];
      fprintf(fp," %6.2f",SNDATA.HOSTGAL_SB_FLUXCALERR[ifilt_obs] ) ;
      NTMP++ ;
      if ( NTMP == 10 ) { fprintf(fp,"\n    ");  NTMP=0; }
    }
    fprintf(fp," # %s \n", ctmp );
  }

  return;

} // end wr_dataformat_text_HOSTGAL


// =====================================================
void  wr_dataformat_text_SNPHOT(FILE *fp) {

  // Created Feb 2021
  // write photometry rows
  
  char OBSKEY[] = "OBS:" ;
  bool ISMODEL_FIXMAG    = ( SNDATA.SIM_MODEL_INDEX == MODEL_FIXMAG );
  bool WRFLAG_BLINDTEST  = SNDATA.WRFLAG_BLINDTEST ;
  bool WRFLAG_PHOTPROB   = SNDATA.WRFLAG_PHOTPROB ;
  bool WRFLAG_PHOTFLAG   = true;
  bool WRFLAG_SKYSIG_T   = SNDATA.WRFLAG_SKYSIG_T ; // template sky sig
  bool WRFLAG_SIM_MAGOBS = ( SNDATA.FAKE > 0 && !SNDATA.WRFLAG_BLINDTEST );
  bool WRFLAG_TRIGGER    = (SNDATA.MJD_TRIGGER < 0.99E6 && 
			    SNDATA.MJD_TRIGGER > 1000.0 );

  double MJD ;
  int  ep, NVAR, NVAR_EXPECT, NVAR_WRITE;
  char VARLIST[200], cvar[40], cval[40], LINE_EPOCH[200] ;
  char fnam[] = "wr_dataformat_text_SNPHOT" ;

  // ------------ BEGIN -----------

  VARLIST[0] = NVAR = 0;
  NVAR++ ;  strcat(VARLIST,"MJD ");  
  NVAR++ ;  strcat(VARLIST,"FLT ");
  NVAR++ ;  strcat(VARLIST,"FIELD ");
  NVAR++ ;  strcat(VARLIST,"FLUXCAL ");
  NVAR++ ;  strcat(VARLIST,"FLUXCALERR ");

  if ( WRFLAG_PHOTFLAG )  { NVAR++ ;  strcat(VARLIST,"PHOTFLAG "); }
  if ( WRFLAG_PHOTPROB )  { NVAR++ ;  strcat(VARLIST,"PHOTPROB "); }

  NVAR++ ;  strcat(VARLIST,"GAIN ");
  NVAR++ ;  strcat(VARLIST,"ZPT ");
  NVAR++ ;  strcat(VARLIST,"PSF ");
  NVAR++ ;  strcat(VARLIST,"SKY_SIG ");
  if ( WRFLAG_SKYSIG_T ) { NVAR++ ;  strcat(VARLIST,"SKY_SIG_T "); }

  if ( WRFLAG_SIM_MAGOBS )
    { NVAR++ ;  strcat(VARLIST,"SIM_MAGOBS "); }

  if ( SNDATA.MAGMONITOR_SNR  ) {
    sprintf(cvar, "SIM_SNRMAG%2.2d", SNDATA.MAGMONITOR_SNR );
    NVAR++ ;  strcat(VARLIST," " );   strcat(VARLIST,cvar);
  }

  /* xxx maybe later ???
  if ( APPEND_MAGREST )
    { strcat(VARLIST,"  FLT_REST  SIM_MAGREST AVWARP KCOR"); }
  xxxxxxx */

  NVAR_EXPECT = NVAR;

  fprintf(fp,"\n\n# ============================================ \n");
  fprintf(fp,"# TEXT LIGHT CURVE OUTPUT: \n#\n");
  fprintf(fp,"NOBS: %d \n", SNDATA.NOBS  );
  fprintf(fp,"NVAR: %d \n", NVAR );
  fprintf(fp,"VARLIST: %s\n", VARLIST);

  for ( ep = 1; ep <= SNDATA.NEPOCH; ep++ ) {

    if ( !SNDATA.OBSFLAG_WRITE[ep] )  { continue ; }

    sprintf(LINE_EPOCH,"%s ", OBSKEY);
    NVAR_WRITE = 0 ;

    sprintf(cval, "%10.4f ",  SNDATA.MJD[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%s ",  SNDATA.FILTCHAR[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%s ",  SNDATA.FIELDNAME[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%11.4le ",  SNDATA.FLUXCAL[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%10.3le ",  SNDATA.FLUXCAL_ERRTOT[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    if ( WRFLAG_PHOTFLAG ) {
      sprintf(cval, "%4d ",  SNDATA.PHOTFLAG[ep] ); 
      NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);
    }
    if ( WRFLAG_PHOTPROB ) {
      sprintf(cval, "%5.3f ",  SNDATA.PHOTPROB[ep] ); 
      NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);
    }

    sprintf(cval, "%6.3f ",  SNDATA.GAIN[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%6.3f ",  SNDATA.ZEROPT[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%5.2f ",  SNDATA.PSF_SIG1[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    sprintf(cval, "%.3le ",  SNDATA.SKY_SIG[ep] ); 
    NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);

    if ( WRFLAG_SKYSIG_T ) {
      sprintf(cval, "%.3le ",  SNDATA.SKY_SIG_T[ep] ); 
      NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);
    }

    if ( WRFLAG_SIM_MAGOBS ) {
      sprintf(cval, "%8.4f ",  SNDATA.SIMEPOCH_MAG[ep] ); 
      NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);
    }

    if ( SNDATA.MAGMONITOR_SNR  ) { 
      sprintf(cval, "%6.1f ",  SNDATA.SIMEPOCH_SNRMON[ep] ); 
      NVAR_WRITE++ ;    strcat(LINE_EPOCH,cval);
    }

    /* xxx ??? maybe someday 
    if ( APPEND_MAGREST ) {    }
    xxxx */

    // check for trigger marker
    MJD = SNDATA.MJD[ep] ;
    if ( WRFLAG_TRIGGER && MJD > SNDATA.MJD_TRIGGER+1.0E-5  ) {
      fprintf(fp,"TRIGGER:  trigger logic satisfied at MJD_TRIGGER=%.4f\n",
              SNDATA.MJD_TRIGGER);
      WRFLAG_TRIGGER = false ;
    }

    if ( NVAR_WRITE != NVAR_EXPECT ) {
      sprintf(c1err,"Prepared to write %d variables for MJD=%f (ep=%d of %d)", 
	      NVAR_WRITE, MJD, ep, SNDATA.NEPOCH );      
      sprintf(c2err,"but NVAR_EXPECT = %d", NVAR_EXPECT );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    // write epoch line
    fprintf(fp, "%s\n", LINE_EPOCH);
  }

  fprintf(fp, "END_PHOTOMETRY: \n\n");
  fflush(fp);

  return ;

} // end wr_dataformat_text_SNPHOT

// =====================================================
void  wr_dataformat_text_SNSPEC(FILE *fp) {

  bool WRFLAG_SIM = (SNDATA.FAKE == FAKEFLAG_LCSIM);
  // xxx mark delete Feb 24 2021 int  NMJD       = GENSPEC.NMJD_TOT ;
  int  NMJD       = GENSPEC.NMJD_PROC ;  // Feb 24 2021
  int  NBLAM_TOT  = GENSPEC.NBLAM_TOT ;
  
  int  NBLAM_VALID, NBLAM_WR, IDSPEC, IS_HOST, NVAR, NVAR_EXPECT ;
  int  imjd, ilam ;
  double L0, L1, LCEN, FLAM, FLAMERR, GENFLAM, GENMAG, WARP ;

  char VARLIST[200], tmpLine[200], cval[40] ;
  char fnam[] = "wr_dataformat_text_SNSPEC" ;

  // ------------ BEGIN -----------

  if ( NMJD == 0 ) { return; }

  VARLIST[0] = NVAR = 0 ;

  NVAR++ ; strcat(VARLIST,"LAMMIN ");
  NVAR++ ; strcat(VARLIST,"LAMMAX ");
  NVAR++ ; strcat(VARLIST,"FLAM ");
  NVAR++ ; strcat(VARLIST,"FLAMERR ");

  if ( WRFLAG_SIM ) {
    NVAR++ ; strcat(VARLIST,"SIM_GENFLAM ");
    NVAR++ ; strcat(VARLIST,"SIM_GENMAG ");
    if ( GENSPEC.USE_WARP )  { NVAR++;  strcat(VARLIST,"SIM_WARP "); }
  }

  // write header info                                                          
  fprintf(fp,"\n# ============================================= \n");
  fprintf(fp,"NSPECTRA:   %d \n\n",  NMJD);
  fprintf(fp,"NVAR_SPEC:  %d \n",    NVAR );
  fprintf(fp,"VARNAMES_SPEC: %s \n", VARLIST);

  for(imjd=0; imjd < NMJD; imjd++ ) {
    IDSPEC = imjd + 1 ;  // start at 1                                          
    NBLAM_VALID = GENSPEC.NBLAM_VALID[imjd] ;
    IS_HOST     = GENSPEC.IS_HOST[imjd];

    if ( NBLAM_VALID == 0 ) { return; } // suppress legacy bug (Aug 23 2017)

    fprintf(fp,"SPECTRUM_ID:       %d  \n", IDSPEC ) ;

    fprintf(fp,"SPECTRUM_MJD:      %9.3f            ", 
	    GENSPEC.MJD_LIST[imjd]);

    if ( IS_HOST )
      { fprintf(fp, "# HOST \n"); }
    else
      { fprintf(fp, "# Tobs = %8.3f \n", GENSPEC.TOBS_LIST[imjd] ); }

    fprintf(fp,"SPECTRUM_TEXPOSE:  %9.1f            "
            "# seconds\n",
            GENSPEC.TEXPOSE_LIST[imjd] );

    fprintf(fp,"SPECTRUM_SNR_COMPUTE:  %9.3f        "
            "# from user SNR request\n",
            GENSPEC.SNR_COMPUTE_LIST[imjd] );

    fprintf(fp,"SPECTRUM_LAMOBS_SNR:   %7.1f %7.1f  "
            "# range for SNR_COMPUTE\n",
            GENSPEC.LAMOBS_SNR_LIST[imjd][0],
            GENSPEC.LAMOBS_SNR_LIST[imjd][1] );

    fprintf(fp,"SPECTRUM_NLAM:       %4d  %4d         "
            "# Number of wave bins: VALID  TOTAL\n",
            NBLAM_VALID, NBLAM_TOT );

    /* ??
    if ( INPUTS.NHOST_TAKE_SPECTRUM > 0 && !IS_HOST ) {
      fprintf(fp,"SPECTRUM_HOSTFRAC:   %.2f               "
              "# host-frac flux contamination\n",
              INPUTS.TAKE_SPECTRUM_HOSTFRAC );
    }
    xxx*/

    NBLAM_WR = 0 ;
    NVAR_EXPECT = NVAR ;

    for(ilam=0; ilam < NBLAM_TOT; ilam++ ) {
      GENFLAM    = GENSPEC.GENFLAM_LIST[imjd][ilam];
      GENMAG     = GENSPEC.GENMAG_LIST[imjd][ilam];
      FLAM       = GENSPEC.FLAM_LIST[imjd][ilam];
      FLAMERR    = GENSPEC.FLAMERR_LIST[imjd][ilam];
      WARP       = GENSPEC.FLAMWARP_LIST[imjd][ilam];

      if ( FLAMERR <= 0.0 ) { continue ; } // skip unphysical values            
      L0      = INPUTS_SPECTRO.LAMMIN_LIST[ilam];
      L1      = INPUTS_SPECTRO.LAMMAX_LIST[ilam];
      LCEN    = INPUTS_SPECTRO.LAMAVG_LIST[ilam];

      NVAR = 0; sprintf(tmpLine,"SPEC: ");

      sprintf(cval, "%8.2f ", L0);
      NVAR++ ; strcat(tmpLine,cval);
      sprintf(cval, "%8.2f ", L1);  
      NVAR++ ; strcat(tmpLine,cval);

      sprintf(cval, "%10.3le ", FLAM);  
      NVAR++ ; strcat(tmpLine,cval);
      sprintf(cval, "%10.3le ", FLAMERR);  
      NVAR++ ; strcat(tmpLine,cval);

      if ( WRFLAG_SIM ) {
	sprintf(cval, "%10.3le ", GENFLAM);  
	NVAR++ ; strcat(tmpLine,cval);

	sprintf(cval, "%.2f ", GENMAG);  
	NVAR++ ; strcat(tmpLine,cval);

	if ( GENSPEC.USE_WARP ) {
	  sprintf(cval,"%6.3f ", WARP );
	  NVAR++ ; strcat(tmpLine,cval);
	}
	  
      } // end WRFLAG_SIM

      fprintf(fp,"%s \n", tmpLine);
      NBLAM_WR++ ;
 
    } // end ilam
  
    fprintf(fp,"SPECTRUM_END:  \n\n" );
  
    if ( NBLAM_WR != NBLAM_VALID ) {
      sprintf(c1err,"Wrote %d lamBins, but expected %d",
	      NBLAM_WR, NBLAM_VALID);
      sprintf(c2err,"CID=%s  SpecID=%d MJD=%.3f",
	      SNDATA.CCID, IDSPEC, GENSPEC.MJD_LIST[imjd] );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    fflush(fp); 
  }  // end imjd loop 

  return ;

} // end wr_dataformat_text_SNSPEC



/*************************************************************

   READ UTILITIES

*************************************************************/

void RD_SNTEXTIO_INIT(void) {
  // Feb 2021: one-time init
  SNTEXTIO_VERSION_INFO.NVERSION        = 0 ;
  SNTEXTIO_VERSION_INFO.NFILE           = 0 ;
  SNTEXTIO_VERSION_INFO.PHOT_VERSION[0] = 0 ;
  SNTEXTIO_VERSION_INFO.DATA_PATH[0]    = 0 ;
  check_head_sntextio(0);

  init_SNDATA_GLOBAL();
  init_SNDATA_EVENT();
  init_GENSPEC_GLOBAL();

  return ;
} // end RD_SNTEXTIO_INIT




// =================================================
int RD_SNTEXTIO_PREP(int MSKOPT, char *PATH, char *VERSION) {

  // Read LIST of files and global info from first data file.
  //
  // Inputs
  //   MSKOPT  
  //     += 8  -> DUMP
  //   PATH    : optional private data path to check
  //   VERSION : name of data version = folder name

  int  istat, NFILE = 0;
  int  LDMP   = (MSKOPT & 8)> 0;
  char fnam[] = "RD_SNTEXTIO_PREP" ;

  // ------------- BEGIN ------------

  sprintf(BANNER,"%s: Prepare reading text format for %s",
	  fnam, VERSION);
  print_banner(BANNER);

  if ( LDMP) {
    printf(" xxx %s: PATH = '%s' \n", fnam, PATH);
    printf(" xxx %s: VERS = '%s' \n", fnam, VERSION);
  }

  sprintf(SNTEXTIO_VERSION_INFO.DATA_PATH,   "%s", PATH );
  sprintf(SNTEXTIO_VERSION_INFO.PHOT_VERSION,"%s", VERSION);

  istat = 
    getInfo_PHOTOMETRY_VERSION(SNTEXTIO_VERSION_INFO.PHOT_VERSION,
			       SNTEXTIO_VERSION_INFO.DATA_PATH,
			       SNTEXTIO_VERSION_INFO.LIST_FILE,
			       SNTEXTIO_VERSION_INFO.README_FILE );

  if ( istat == ERROR ) {
    sprintf(c1err,"Cannot find VERSION_PHOTOMETRY = '%s'", VERSION);
    sprintf(c2err,"Check $SNDTA_ROOT/lcmerge, PATH_SNDATA_SIM, "
	    "PRIVATE_DATA_PATH");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // read & store list of data files 
  NFILE = rd_sntextio_list() ;
  if ( NFILE < 0 ) { return -1 ; } // not TEXT format

  // read/store global info from first file
  

  if ( SNTEXTIO_VERSION_INFO.NVERSION == 0 ) 
    { rd_sntextio_global(); }

  SNTEXTIO_VERSION_INFO.NVERSION++ ; 

  return NFILE;

} //end RD_SNTEXTIO_PREP


// - - - - - 
void rd_sntextio_init__(void) { RD_SNTEXTIO_INIT(); }

int rd_sntextio_prep__(int *MSKOPT, char *PATH, char *VERSION)
{ return RD_SNTEXTIO_PREP(*MSKOPT, PATH,VERSION); }


// ===========================================
int rd_sntextio_list(void) {

  // Created Feb 2021
  // Read auxilary LIST file and count how many text files are specified.
  // Return NFILE = -1 if FITS extension is detected.
  // Open first file and read header to store global info 
  //  e..g, SURVEY, FILTERS, ...

  char *LIST_FILE = SNTEXTIO_VERSION_INFO.LIST_FILE ;
  char *DATA_PATH = SNTEXTIO_VERSION_INFO.DATA_PATH ;
  int   MSKOPT     = MSKOPT_PARSE_TEXT_FILE ;

  int  langC  = LANGFLAG_PARSE_WORDS_C ;
  int  iwd, NFILE = 0 ;
  int  RETCODE_FITS = -1;
  int  NFILE_LAST = SNTEXTIO_VERSION_INFO.NFILE ;

  char firstFile[MXPATHLEN], FIRSTFILE[MXPATHLEN];
  char fnam[] = "rd_sntextio_list" ;

  // ------------- BEGIN --------------

  NFILE = store_PARSE_WORDS(MSKOPT, LIST_FILE);

  // read first file
  iwd=0;  get_PARSE_WORD(langC, iwd, firstFile);

  sprintf(FIRSTFILE, "%s/%s", DATA_PATH, firstFile);	 

  //  printf(" xxx %s: firstFile = \n\t '%s' \n\t '%s' \n", 
  //	 fnam, firstFile, FIRSTFILE );

  // return -1 if this has fits extension
  if ( strstr(firstFile,".FITS") != NULL ) { return(RETCODE_FITS); }
  if ( strstr(firstFile,".fits") != NULL ) { return(RETCODE_FITS); }

  // not fits, so store list of text file names to read.
  // Just store base file name in VERSION.LIST ...
  // path is appended later when the data file is read.
  bool DO_FREE = ( SNTEXTIO_VERSION_INFO.NVERSION > 0 ) ;
  if ( DO_FREE ) { rd_sntextio_malloc_list(-1, NFILE_LAST); }
  rd_sntextio_malloc_list(+1, NFILE);

  for(iwd=0; iwd < NFILE; iwd++ ) {
    get_PARSE_WORD(langC, iwd, 
		   SNTEXTIO_VERSION_INFO.DATA_FILE_LIST[iwd] );
  }

  SNTEXTIO_VERSION_INFO.NFILE = NFILE ; 

  return NFILE;

} // end  rd_sntextio_list


// =============================================
void rd_sntextio_malloc_list(int OPT, int NFILE) {
  
  int i;
  char fnam[] = "rd_sntextio_malloc_list" ;
  // ---------- BEGIN ----------


  if ( OPT < 0 ) {
    // free mem
    for(i=0; i < NFILE; i++ )
      { free(SNTEXTIO_VERSION_INFO.DATA_FILE_LIST[i]);   }
    free(SNTEXTIO_VERSION_INFO.DATA_FILE_LIST);
  }
  else {
    // malloc
    int MEMC1 = NFILE * sizeof(char*);
    int MEMC0 = 80    * sizeof(char);
    SNTEXTIO_VERSION_INFO.DATA_FILE_LIST = (char**) malloc(MEMC1);
    for(i=0; i < NFILE; i++ )
      { SNTEXTIO_VERSION_INFO.DATA_FILE_LIST[i] =  (char*)malloc(MEMC0); }    
  }

  return;
} // end rd_sntextio_malloc_list

// =======================================
void  rd_sntextio_global(void) {

  // Created Feb 2021
  // Open first text file and read info that is global;
  // skip SN-dependent info. Stop reading at first "OBS:" key.

  int   MSKOPT     = MSKOPT_PARSE_TEXT_FILE ;
  int  NVERSION    = SNTEXTIO_VERSION_INFO.NVERSION ;
  char *firstFile  = SNTEXTIO_VERSION_INFO.DATA_FILE_LIST[0];
  char *DATA_PATH  = SNTEXTIO_VERSION_INFO.DATA_PATH ;
  char FIRSTFILE[MXPATHLEN] ;
  int  langC  = LANGFLAG_PARSE_WORDS_C ;
  int  NWD, iwd, ITMP, LENKEY, NVAR, NPAR ;
  bool HAS_COLON, HAS_PARENTH, IS_TMP, IS_SIM, FOUND_FAKEKEY=false ;
  bool IS_PRIVATE, IS_SIMSED, IS_LCLIB, IS_BYOSED, IS_SNEMO ;
  char word0[100], word1[100], word2[100];    
  char fnam[] = "rd_sntextio_global" ;
  int  LDMP = 0 ;
  // ---------- BEGIN ----------

  sprintf(FIRSTFILE, "%s/%s", DATA_PATH, firstFile);	 
  NWD = store_PARSE_WORDS(MSKOPT,FIRSTFILE);
  
  if ( LDMP ) {
    printf(" xxx %s: store %d words from \n\t %s\n", fnam, NWD, FIRSTFILE);
  }

  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(langC, iwd, word0 );
    
    HAS_COLON    =  strstr(word0,COLON) != NULL  ;
    HAS_PARENTH  =  strstr(word0,"("  ) != NULL  ;
    IS_TMP       =  HAS_COLON && HAS_PARENTH ;
    IS_PRIVATE   =  strstr(word0,"PRIVATE")       != NULL && IS_TMP ;
    IS_SIMSED    =  strstr(word0,"SIMSED_")       != NULL && IS_TMP ;
    IS_LCLIB     =  strstr(word0,"LCLIB_PARAM")   != NULL && IS_TMP ;
    IS_BYOSED    =  strstr(word0,"BYOSED_PARAM")  != NULL && IS_TMP ;
    IS_SNEMO     =  strstr(word0,"SNEMO_PARAM")   != NULL && IS_TMP ;

    IS_SIM       =  strncmp(word0,"SIM",3)   == 0    && HAS_COLON ;

    if( strcmp(word0,"SNANA_VERSION:") == 0 ) {
      iwd++; get_PARSE_WORD(langC, iwd, SNDATA.SNANA_VERSION );
    }
    else if ( strcmp(word0,"SURVEY:") == 0 ) {
      iwd++; get_PARSE_WORD(langC, iwd, SNDATA.SURVEY_NAME );

      // check for SURVEY(SUBSURVEY); e.g., LOWZ_COMBINED(CFA3)
      extractStringOpt(SNDATA.SURVEY_NAME, SNDATA.SUBSURVEY_NAME); 
    }

    else if( strcmp(word0,"FILTERS:") == 0 ) {
      iwd++; get_PARSE_WORD(langC, iwd, word1);
      set_SNDATA_FILTER(word1);
    }
    else if ( strcmp(word0,"FAKE:") == 0 ) {
      iwd++; get_PARSE_WORD_INT(langC, iwd, &ITMP );
      SNDATA.FAKE = ITMP ;
      FOUND_FAKEKEY = true ;
      if ( ITMP == FAKEFLAG_DATA ) 
	{ sprintf(SNDATA.DATATYPE, "%s", DATATYPE_DATA); } 
      else if ( ITMP == FAKEFLAG_LCSIM ) 
	{ sprintf(SNDATA.DATATYPE, "%s", DATATYPE_SIM_SNANA); } 
      else if ( ITMP == FAKEFLAG_FAKES ) 
	{ sprintf(SNDATA.DATATYPE, "%s", DATATYPE_SIM_MAGOBS); } 
    }
    else if ( IS_PRIVATE ) {
      SNDATA.NVAR_PRIVATE++ ; // note fortran-like index
      NVAR = SNDATA.NVAR_PRIVATE ;
      copy_keyword_nocolon(word0,SNDATA.PRIVATE_KEYWORD[NVAR]);
    }
    else if ( IS_SIMSED ) {
      NPAR = SNDATA.NPAR_SIMSED;
      if ( strcmp(word0,"SIMSED_NPAR:") == 0 ) { continue; }  
      copy_keyword_nocolon(word0,SNDATA.SIMSED_KEYWORD[NPAR]);
      //      sprintf(SNDATA.SIMSED_KEYWORD[NPAR], "%s", word0);
      NPAR++ ;    SNDATA.NPAR_SIMSED = NPAR ;
    }
    else if ( IS_LCLIB ) {
      NPAR = SNDATA.NPAR_LCLIB ;
      if ( strcmp(word0,"LCLIB_NPAR:") == 0 ) { continue; }  
      copy_keyword_nocolon(word0,SNDATA.LCLIB_KEYWORD[NPAR]);
      //      sprintf(SNDATA.LCLIB_KEYWORD[NPAR], "%s", word0);
      NPAR++ ;    SNDATA.NPAR_LCLIB = NPAR ;
    }
    else if ( IS_BYOSED || IS_SNEMO ) {
      NPAR = SNDATA.NPAR_PySEDMODEL ;
      if ( strcmp(word0,"BYOSED_NPAR:") == 0 ) { continue; }  
      if ( strcmp(word0,"SNEMO_NPAR:" ) == 0 ) { continue; }  
      copy_keyword_nocolon(word0,SNDATA.PySEDMODEL_KEYWORD[NPAR]);
      NPAR++ ;    SNDATA.NPAR_PySEDMODEL = NPAR ;
      if ( IS_BYOSED ) { sprintf(SNDATA.PySEDMODEL_NAME, "BYOSED"); }
      if ( IS_SNEMO  ) { sprintf(SNDATA.PySEDMODEL_NAME, "SNEMO" ); }
    } 

    
    else if ( IS_SIM ) {
      if ( strcmp(word0,"SIMLIB_MSKOPT:") == 0 ) {
	iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIMLIB_MSKOPT) ;
      }
      else if ( strcmp(word0,"SIMOPT_MWCOLORLAW:") == 0 ) {
	iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIMOPT_MWCOLORLAW ) ;
      }
      else if ( strcmp(word0,"SIMOPT_MWEBV:") == 0 ) {
	iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIMOPT_MWEBV ) ;
      }

      else if ( strcmp(word0,"SIM_MWRV:") == 0 ) {
	iwd++; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_MWRV ) ;
      }
      else if ( strstr(word0,"SIM_STRONGLENS") != NULL ) {
	SNDATA.SIM_SL_FLAG = 1 ;
      }

    } // end IS_SIM

    // - - - - - - OBS info - - - - -
    else if ( strcmp(word0,"NVAR:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNTEXTIO_FILE_INFO.NVAROBS );
    }
    else if ( strcmp(word0,"VARLIST:") == 0 ) {
      rd_sntextio_varlist_obs(&iwd);
    }

    else if ( strcmp(word0,"OBS:") == 0 ) {

      if ( !FOUND_FAKEKEY) {
	printf("\n     WARNING: no FAKE key -> assume real data\n\n");
	SNDATA.FAKE = FAKEFLAG_DATA;
	sprintf(SNDATA.DATATYPE, "%s", DATATYPE_DATA); 
      }
      return ;  // done reading global info; bye bye
    }    

  } // end iwd loop
  
  return;

} // end rd_sntextio_global

// =============================
void copy_keyword_nocolon(char *key_in, char *key_out) {
  int lenkey = strlen(key_in);
  sprintf(key_out, "%s", key_in) ;
  if ( strstr(key_out,COLON) != NULL ) 
    { sprintf(&key_out[lenkey-1],"") ; }
  return;
} // end copy_keyword_nocolon

// ==============================================
void rd_sntextio_varlist_obs(int *iwd_file) {

  // Created Feb 2021
  // read OBS-varlist elements to prepare for reading OBS later.


  int  iwd    = *iwd_file;
  int  langC  = LANGFLAG_PARSE_WORDS_C ;
  int  NVAR, ivar ;
  char *varName ;
  char fnam[] = "rd_sntextio_varlist_obs" ;

  // ---------- BEGIN -------

  IVAROBS_SNTEXTIO.MJD = IVAROBS_SNTEXTIO.BAND = IVAROBS_SNTEXTIO.FIELD = -9 ;
  IVAROBS_SNTEXTIO.FLUXCAL = IVAROBS_SNTEXTIO.FLUXCALERR = -9 ;
  IVAROBS_SNTEXTIO.ZPFLUX = IVAROBS_SNTEXTIO.ZPERR = IVAROBS_SNTEXTIO.PSF = -9;
  IVAROBS_SNTEXTIO.SKYSIG = IVAROBS_SNTEXTIO.SKYSIG_T = -9;
  IVAROBS_SNTEXTIO.GAIN = -9;
  IVAROBS_SNTEXTIO.PHOTFLAG = IVAROBS_SNTEXTIO.PHOTPROB = -9 ;
  IVAROBS_SNTEXTIO.XPIX = IVAROBS_SNTEXTIO.YPIX = IVAROBS_SNTEXTIO.CCDNUM = -9;
  IVAROBS_SNTEXTIO.SIMEPOCH_MAG -9 ;

  NVAR = SNTEXTIO_FILE_INFO.NVAROBS ;
  if ( NVAR < 5 || NVAR >= MXVAROBS_TEXT ) {
    sprintf(c1err,"Invalid NVAR=%d (MXVAROBS_TEXT=%d)", 
	    NVAR, MXVAROBS_TEXT ) ;
    sprintf(c2err,"Check NVAR and VARLIST keys");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // printf(" xxx %s: read %d VARLIST names\n", fnam, SNTEXTIO_FILE_INFO.NVAROBS );

  for(ivar=0; ivar < NVAR; ivar++ ) {
    varName = SNTEXTIO_FILE_INFO.VARNAME_OBS_LIST[ivar];
    iwd++ ; get_PARSE_WORD(langC, iwd, varName);
    // printf(" xxx %s: varName[%2d] = %s \n", fnam, ivar, varName);

    if ( strcmp(varName,"MJD") == 0 ) 
      { IVAROBS_SNTEXTIO.MJD = ivar; }

    else if ( strcmp(varName,"BAND") == 0 ) 
      { IVAROBS_SNTEXTIO.BAND = ivar; }
    else if ( strcmp(varName,"FLT") == 0 ) 
      { IVAROBS_SNTEXTIO.BAND = ivar; }

    else if ( strcmp(varName,"FIELD") == 0 ) 
      { IVAROBS_SNTEXTIO.FIELD = ivar; }
    else if ( strcmp(varName,"FLUXCAL") == 0 )  
      { IVAROBS_SNTEXTIO.FLUXCAL = ivar; }
    else if ( strcmp(varName,"FLUXCALERR") == 0 ) 
      { IVAROBS_SNTEXTIO.FLUXCALERR = ivar; }

    // check obsolete keys left in very old data files
    else if ( strcmp(varName,"MAG")    == 0 )  { ; }  
    else if ( strcmp(varName,"MAGERR") == 0 )  { ; }
    else if ( strcmp(varName,"SNR")    == 0 )  { ; } // in SNLS3year_MEGACAM

    else if ( strcmp(varName,"PHOTFLAG") == 0 ) 
      { IVAROBS_SNTEXTIO.PHOTFLAG = ivar; }  
    else if ( strcmp(varName,"PHOTPROB") == 0 ) 
      { IVAROBS_SNTEXTIO.PHOTPROB = ivar; }  


    else if ( strcmp(varName,"PSF") == 0 ) 
      { IVAROBS_SNTEXTIO.PSF = ivar; }   

    else if ( strcmp(varName,"ZPFLUX") == 0 ||
	      strcmp(varName,"ZPT")    == 0 ||
	      strcmp(varName,"Zpt")    == 0 ) 
      { IVAROBS_SNTEXTIO.ZPFLUX = ivar; }  

    else if ( strcmp(varName,"ZPERR") == 0 ) 
      { IVAROBS_SNTEXTIO.ZPERR = ivar; }  
    else if ( strcmp(varName,"SKY_SIG") == 0 || strcmp(varName,"SKYSIG")==0 ) 
      { IVAROBS_SNTEXTIO.SKYSIG = ivar; }  
    else if (strcmp(varName,"SKY_SIG_T") == 0 ||strcmp(varName,"SKYSIG_T")==0) 
      { IVAROBS_SNTEXTIO.SKYSIG_T = ivar; }

    else if ( strcmp(varName,"GAIN") == 0 ) 
      { IVAROBS_SNTEXTIO.GAIN = ivar; }  

    else if ( strcmp(varName,"XPIX") == 0 ) 
      { IVAROBS_SNTEXTIO.XPIX = ivar; }  
    else if ( strcmp(varName,"YPIX") == 0 ) 
      { IVAROBS_SNTEXTIO.YPIX = ivar; }  
    else if ( strcmp(varName,"CCDNUM") == 0 ) 
      { IVAROBS_SNTEXTIO.CCDNUM = ivar; }  

    else if ( strcmp(varName,"SIM_MAGOBS") == 0 ) 
      { IVAROBS_SNTEXTIO.SIMEPOCH_MAG = ivar; }  

    else {
      sprintf(c1err,"Invalid varName = %s (ivar=%d of %d)", 
	      varName, ivar, NVAR );
      sprintf(c2err,"Check VARLIST args.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  } // end ivar


  // - - - - - - - - - - - - - - - 
  // check that required columns are defined
  int IVAR_MIN  = 0;
  int IVAR_MAX  = MXVAROBS_TEXT-1 ;
  int NVAL=1;
  checkval_I("IVAROBS_MJD", NVAL, &IVAROBS_SNTEXTIO.MJD, 
	     IVAR_MIN, IVAR_MAX);
  checkval_I("IVAROBS_BAND", NVAL, &IVAROBS_SNTEXTIO.BAND, 
	     IVAR_MIN, IVAR_MAX);
  checkval_I("IVAROBS_FLUXCAL", NVAL, &IVAROBS_SNTEXTIO.FLUXCAL, 
	     IVAR_MIN, IVAR_MAX);
  checkval_I("IVAROBS_FLUXCALERR", NVAL, &IVAROBS_SNTEXTIO.FLUXCALERR, 
	     IVAR_MIN, IVAR_MAX);

  *iwd_file = iwd;

  return;
} // end rd_sntextio_varList_obs


// ==============================================
void rd_sntextio_varlist_spec(int *iwd_file) {

  // Created Feb 2021
  // read SPEC-varlist elements to prepare for reading OBS later.


  int  iwd    = *iwd_file;
  int  langC  = LANGFLAG_PARSE_WORDS_C ;
  int  NVAR, ivar ;
  char *varName ;
  char fnam[] = "rd_sntextio_varlist_spec" ;

  // ---------- BEGIN -------

  IVARSPEC_SNTEXTIO.LAMMIN = IVARSPEC_SNTEXTIO.LAMMAX = -9;
  IVARSPEC_SNTEXTIO.LAMAVG = -9;
  IVARSPEC_SNTEXTIO.LAMMIN = IVARSPEC_SNTEXTIO.FLAMERR = -9;
  IVARSPEC_SNTEXTIO.SIM_GENFLAM = IVARSPEC_SNTEXTIO.SIM_GENMAG = -9;


  NVAR = SNTEXTIO_FILE_INFO.NVARSPEC ;
  if ( NVAR < 3 || NVAR >= MXVAROBS_TEXT ) {
    sprintf(c1err,"Invalid NVAR=%d (MXVARSPEC_TEXT=%d) for CID=%s", 
	    NVAR, MXVAROBS_TEXT, SNDATA.CCID ) ;
    sprintf(c2err,"Check NVAR_SPEC and VARNAMES_SPEC keys");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }



  for(ivar=0; ivar < NVAR; ivar++ ) {
    varName = SNTEXTIO_FILE_INFO.VARNAME_SPEC_LIST[ivar];
    iwd++ ; get_PARSE_WORD(langC, iwd, varName);
    // printf(" xxx %s: varName[%2d] = %s \n", fnam, ivar, varName);

    if ( strcmp(varName,"LAMAVG") == 0 ) 
      { IVARSPEC_SNTEXTIO.LAMAVG = ivar; }

    else if ( strcmp(varName,"LAMMIN") == 0 ) 
      { IVARSPEC_SNTEXTIO.LAMMIN = ivar; }

    else if ( strcmp(varName,"LAMMAX") == 0 ) 
      { IVARSPEC_SNTEXTIO.LAMMAX = ivar; }

    else if ( strcmp(varName,"FLAM") == 0 ) 
      { IVARSPEC_SNTEXTIO.FLAM = ivar; }

    else if ( strcmp(varName,"FLAMERR") == 0 ) 
      { IVARSPEC_SNTEXTIO.FLAMERR = ivar; }

    else if ( strcmp(varName,"SIM_GENFLAM") == 0 ) 
      { IVARSPEC_SNTEXTIO.SIM_GENFLAM = ivar; }

    else if ( strcmp(varName,"SIM_GENMAG") == 0 ) 
      { IVARSPEC_SNTEXTIO.SIM_GENMAG = ivar; }

    else if ( strcmp(varName,"DQ")       == 0 ) { ; } // do nothing
    else if ( strcmp(varName,"SPECFLAG") == 0 ) { ; } // do nothing

    else {
      sprintf(c1err,"Invalid varName = %s (ivar=%d of %d)", 
	      varName, ivar, NVAR );
      sprintf(c2err,"Check VARLIST_SPEC args for CID=%s", SNDATA.CCID);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  } // end ivar


  // - - - - - - - - - - - - - - - 
  // check that required columns are defined
  int IVAR_MIN  = 0;
  int IVAR_MAX  = MXVAROBS_TEXT-1 ;
  int NVAL=1;
  checkval_I("IVAROBS_FLAM", NVAL, &IVAROBS_SNTEXTIO.MJD, 
	     IVAR_MIN, IVAR_MAX);
  checkval_I("IVAROBS_FLAMERR", NVAL, &IVAROBS_SNTEXTIO.BAND, 
	     IVAR_MIN, IVAR_MAX);

  *iwd_file = iwd;

  return;
} // end rd_sntextio_varlist_spec

// ==============================================
void RD_SNTEXTIO_EVENT(int OPTMASK, int ifile_inp) {

  // Created Feb 2021
  // Read event for DATA_FILE_LIST[ifile]
  // OPTMASK controls reading of header, obs, spec,
  // so that header cuts can be applied before read OBS.
  // Beware that input ifile_inp runs from 1 to NFILE;
  // so define file = ifile_inp-1 to use as C index


  int  ifile      = ifile_inp - 1; // convert to C index starting at 0
  bool LRD_HEAD   = (OPTMASK & OPTMASK_TEXT_HEAD) > 0 ;
  bool LRD_OBS    = (OPTMASK & OPTMASK_TEXT_OBS ) > 0 ;
  bool LRD_SPEC   = (OPTMASK & OPTMASK_TEXT_SPEC) > 0 ;
  int  MSKOPT     = MSKOPT_PARSE_TEXT_FILE ;
  int  NFILE_TOT  = SNTEXTIO_VERSION_INFO.NFILE ;
  char *DATA_PATH = SNTEXTIO_VERSION_INFO.DATA_PATH ;
  char *fileName  = SNTEXTIO_VERSION_INFO.DATA_FILE_LIST[ifile];
  
  char FILENAME[MXPATHLEN]; 
  int  NWD, iwd; 
  bool LRD_NEXT = false;
  char fnam[] = "RD_SNTEXTIO_EVENT";

  // ------------ BEGIN ----------

  if ( ifile < 0 || ifile >= NFILE_TOT ) {
    sprintf(c1err,"Invalid ifile=%d", ifile);
    sprintf(c2err,"Valid ifile is 0 to %d", NFILE_TOT-1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if ( LRD_HEAD ) {
    sprintf(FILENAME, "%s/%s", DATA_PATH, fileName);
    NWD = store_PARSE_WORDS(MSKOPT,FILENAME);

    SNTEXTIO_FILE_INFO.NWD_TOT    = NWD ;
    SNTEXTIO_FILE_INFO.IPTR_READ  = 0 ;
    SNTEXTIO_FILE_INFO.NOBS_READ  = 0 ;
    SNTEXTIO_FILE_INFO.NSPEC_READ = 0 ;
    SNTEXTIO_FILE_INFO.NVAR_PRIVATE_READ = 0 ;
    init_SNDATA_EVENT();
    init_GENSPEC_EVENT(-1,-1);
    check_head_sntextio(1);

    iwd = 0;  LRD_NEXT = true ;
    while ( LRD_NEXT && iwd < NWD ) {
      LRD_NEXT = parse_SNTEXTIO_HEAD(&iwd);
      SNTEXTIO_FILE_INFO.IPTR_READ = iwd;
      iwd++ ;  
    }

    // run header sanity checks to catch common user mistakes
    // when making text formatted data files.
    check_head_sntextio(2);

  } // end LRD_HEAD


  // - - - - - - - 
  if ( LRD_OBS ) {

    NWD = SNTEXTIO_FILE_INFO.NWD_TOT ;   // restore from LRD_HEAD
    iwd = SNTEXTIO_FILE_INFO.IPTR_READ;  // restore from LRD_HEAD

    /*
    printf(" xxx %s: read obs for %s (iwd=%d of %d)\n", 
    fnam, SNDATA.CCID, iwd, NWD); */
    LRD_NEXT = true ;
    while ( LRD_NEXT && iwd < NWD ) {
      LRD_NEXT = parse_SNTEXTIO_OBS(&iwd);
      SNTEXTIO_FILE_INFO.IPTR_READ = iwd;
      iwd++ ;  
    }

    int NOBS_READ   =  SNTEXTIO_FILE_INFO.NOBS_READ ;
    int NOBS_EXPECT =  SNDATA.NOBS; 
    if ( NOBS_READ != NOBS_EXPECT ) {
      sprintf(c1err,"Read %d OBS rows for CID=%s", 
	      NOBS_READ, SNDATA.CCID );
      sprintf(c2err,"but expected %d rows from NOBS key.", NOBS_EXPECT);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);      
    }

  } // end LRD_OBS

  // - - - -
  if ( LRD_SPEC ) {

    NWD = SNTEXTIO_FILE_INFO.NWD_TOT ;   // restore from LRD_OBS
    iwd = SNTEXTIO_FILE_INFO.IPTR_READ;  // restore from LRD_OBS
    /*
    printf(" xxx %s: read obs for %s (iwd=%d of %d)\n", 
    fnam, SNDATA.CCID, iwd, NWD); */
    LRD_NEXT = true ;
    while ( LRD_NEXT && iwd < NWD ) {
      LRD_NEXT = parse_SNTEXTIO_SPEC(&iwd);
      SNTEXTIO_FILE_INFO.IPTR_READ = iwd;
      iwd++ ;  
    }
  }     // end LRD_SPEC

  //  debugexit(fnam); // xxx REMOVE

  return;

} // end RD_SNTEXTIO_EVENT

void rd_sntextio_event__(int *OPTMASK, int *ifile)
{ RD_SNTEXTIO_EVENT(*OPTMASK,*ifile); }


bool parse_SNTEXTIO_HEAD(int *iwd_file) {

  // Created Feb 15 2021
  //
  // Data file has already been read and each word is stored.
  // Here, examine word(*iwd_file) for standard key.
  // If standard key, read argument (next word after key)
  // and load SNDATA struct.
  // Also increment and return *iwd_file.
  //
  // Abort if SURVEY and FILTERS arg changes w.r.t. first data file.
  //
  // Function returns true to keep reading;
  // returns false when end of header is reached by finding NOBS key.

  int  NFILT     = SNDATA_FILTER.NDEF;
  int  langC     = LANGFLAG_PARSE_WORDS_C ;
  char *PySEDMODEL_NAME = SNDATA.PySEDMODEL_NAME ;
  int  len_PySEDMODEL   = strlen(PySEDMODEL_NAME);
  int  ncmp_PySEDMODEL ;

  int  iwd       = *iwd_file ;
  int  igal, ivar, NVAR, ipar, NPAR, ifilt, ifilt_obs ;
  double DVAL; 
  bool IS_PRIVATE ;
  char word0[100], PREFIX[40], KEY_TEST[80], ARG_TMP[80];
  char fnam[] = "parse_SNTEXTIO_HEAD" ;

  // ------------ BEGIN -----------g

  get_PARSE_WORD(langC, iwd, word0);

  // bail when reaching first obs
  if ( strcmp(word0,"OBS:") == 0 ) { return false; }

  IS_PRIVATE =  
    strncmp(word0,"PRIVATE",7) == 0  &&  strstr(word0,COLON) != NULL ;

  ncmp_PySEDMODEL  = strncmp(word0,PySEDMODEL_NAME,len_PySEDMODEL) ;

  // - - - - - - -
  // parse keys for data or sim
  if ( strcmp(word0,"SNID:") == 0 ) {
    iwd++ ; get_PARSE_WORD(langC, iwd, SNDATA.CCID);
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_SNID] = true ;
  }
  else if ( strcmp(word0,"SURVEY:") == 0 ) {
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_SURVEY] = true ;

    // check for SURVEY(SUBSURVEY); e.g., LOWZ_COMBINED(CFA3)
    iwd++ ; get_PARSE_WORD(langC, iwd, ARG_TMP);
    extractStringOpt(ARG_TMP, SNDATA.SUBSURVEY_NAME); 
    if ( strcmp(ARG_TMP,SNDATA.SURVEY_NAME) != 0 ) {
      sprintf(c1err,"Invalid 'SURVEY: %s' for CID=%s", ARG_TMP, SNDATA.CCID);
      sprintf(c2err,"Expected 'SURVEY: %s' from first data file",
	      SNDATA.SURVEY_NAME);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }
  else if ( strcmp(word0,"FILTERS:") == 0 ) {
    // already parsed in global; check here that FILTERS arg doesn't change
    iwd++; get_PARSE_WORD(langC, iwd, ARG_TMP);
    if ( strcmp(ARG_TMP,SNDATA_FILTER.LIST) != 0 ) {
      sprintf(c1err,"Invalid 'FILTERS: %s' for CID=%s", ARG_TMP, SNDATA.CCID);
      sprintf(c2err,"Expected 'FILTERS: %s' from first data file",
	      SNDATA_FILTER.LIST);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_FILTERS] = true ;
  }
  else if ( strcmp(word0,"IAUC:") == 0 ) {
    iwd++ ; get_PARSE_WORD(langC, iwd, SNDATA.IAUC_NAME);
  }
  else if ( strcmp(word0,"FAKE:") == 0 ) {
    iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.FAKE);
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_FAKE] = true ;
  }
  else if ( strcmp(word0,"MASK_FLUXCOR_SNANA:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.MASK_FLUXCOR );
  }

  else if ( strcmp(word0,"RA:") == 0 ) {
    iwd++; get_PARSE_WORD_DBL(langC, iwd, &SNDATA.RA ) ;
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_RA] = true ;
  }
  else if ( strcmp(word0,"DEC:") == 0 || strcmp(word0,"DECL:") == 0 ) {
    iwd++; get_PARSE_WORD_DBL(langC, iwd, &SNDATA.DEC );
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_DEC] = true ;
  }
  else if ( strcmp(word0,"PIXSIZE:") == 0 ) {
    iwd++; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.PIXSIZE );
  }
  else if ( strcmp(word0,"NXPIX:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.NXPIX );
  }
  else if ( strcmp(word0,"NYPIX:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.NYPIX );
  }
  else if ( strcmp(word0,"CCDNUM:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.CCDNUM[1] );
  }
  else if ( strcmp(word0,"TYPE:")==0 || strcmp(word0,"SNTYPE:")==0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SNTYPE );
  }

  else if ( strstr(word0,"MWEBV") != NULL ) {
    parse_plusminus_sntextio(word0, "MWEBV", &iwd, 
			     &SNDATA.MWEBV, &SNDATA.MWEBV_ERR );
    //    printf(" xxx %s: MWEBV = %f +- %f (CID=%s) \n",
    //	   fnam, SNDATA.MWEBV, SNDATA.MWEBV_ERR, SNDATA.CCID );
  }
  
  else if ( strstr(word0,"REDSHIFT_HELIO") != NULL ) {
    parse_plusminus_sntextio(word0, "REDSHIFT_HELIO", &iwd, 
			     &SNDATA.REDSHIFT_HELIO, 
			     &SNDATA.REDSHIFT_HELIO_ERR );
  }
  else if ( strstr(word0,"REDSHIFT_FINAL") != NULL ) {
    parse_plusminus_sntextio(word0, "REDSHIFT_FINAL", &iwd, 
			     &SNDATA.REDSHIFT_FINAL, 
			     &SNDATA.REDSHIFT_FINAL_ERR );
    SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[HEAD_REQUIRE_z] = true ;
  }
  
  else if ( strstr(word0,"VPEC") != NULL ) {
    parse_plusminus_sntextio(word0, "VPEC", &iwd, 
			     &SNDATA.VPEC, &SNDATA.VPEC_ERR );
    //	printf(" xxx %s: VPEC = %f +- %f \n",
    //     fnam, SNDATA.VPEC, SNDATA.VPEC_ERR );
  }

  else if ( strncmp(word0,"HOSTGAL",7) == 0 ) {
    
    if ( strcmp(word0,"HOSTGAL_NMATCH:") == 0 ) {
      iwd++; get_PARSE_WORD_INT(langC,iwd, &SNDATA.HOSTGAL_NMATCH[0] );
    }
    else if ( strcmp(word0,"HOSTGAL_NMATCH2:") == 0 ) {
      iwd++; get_PARSE_WORD_INT(langC,iwd, &SNDATA.HOSTGAL_NMATCH[1] );
    } 
    else if ( strcmp(word0,"HOSTGAL_CONFUSION:") == 0 ) {
      iwd++; get_PARSE_WORD_FLT(langC,iwd, &SNDATA.HOSTGAL_CONFUSION );
    }
    
    else if ( strcmp(word0,"HOSTGAL_SB_FLUXCAL:") == 0 ) {
      for(ifilt=0; ifilt < NFILT; ifilt++ ) {
	iwd++; get_PARSE_WORD_FLT(langC, iwd, 
				  &SNDATA.HOSTGAL_SB_FLUXCAL[ifilt]); 
      } 
    } // end HOSTGAL_SB_FLUXCAL
    
    for(igal=0; igal < MXHOSTGAL; igal++ ) {
      sprintf(PREFIX,"HOSTGAL");
      if ( igal > 0 ) { sprintf(PREFIX,"HOSTGAL%d",igal+1); }
      
      sprintf(KEY_TEST,"%s_OBJID:", PREFIX); 
      if ( strcmp(word0,KEY_TEST) == 0 ) {
	iwd++; get_PARSE_WORD_DBL(langC,iwd, &DVAL);
	SNDATA.HOSTGAL_OBJID[igal] = (long long) DVAL ;
      }

      sprintf(KEY_TEST,"%s_PHOTOZ", PREFIX); 
      if ( strstr(word0,KEY_TEST) != NULL ) {
	parse_plusminus_sntextio(word0, KEY_TEST, &iwd, 
				 &SNDATA.HOSTGAL_PHOTOZ[igal], 
				 &SNDATA.HOSTGAL_PHOTOZ_ERR[igal] );
      }
      sprintf(KEY_TEST,"%s_SPECZ", PREFIX); 
      if ( strstr(word0,KEY_TEST) != NULL ) {
	parse_plusminus_sntextio(word0, KEY_TEST, &iwd, 
				 &SNDATA.HOSTGAL_SPECZ[igal], 
			     &SNDATA.HOSTGAL_SPECZ_ERR[igal] );
      }

      sprintf(KEY_TEST,"%s_RA:", PREFIX); 
      if ( strcmp(word0,KEY_TEST) == 0 ) {
	iwd++; get_PARSE_WORD_DBL(langC,iwd, &SNDATA.HOSTGAL_RA[igal]); 
      }
      sprintf(KEY_TEST,"%s_DEC:", PREFIX); 
      if ( strcmp(word0,KEY_TEST) == 0 ) {
	iwd++; get_PARSE_WORD_DBL(langC,iwd, &SNDATA.HOSTGAL_DEC[igal]); 
      }

      sprintf(KEY_TEST,"%s_SNSEP:", PREFIX); 
      if ( strcmp(word0,KEY_TEST) == 0 ) {
	iwd++; get_PARSE_WORD_FLT(langC,iwd, &SNDATA.HOSTGAL_SNSEP[igal]); 
      }
      sprintf(KEY_TEST,"%s_DDLR:", PREFIX); 
      if ( strcmp(word0,KEY_TEST) == 0 ) {
	iwd++; get_PARSE_WORD_FLT(langC,iwd, &SNDATA.HOSTGAL_DDLR[igal]); 
      }

      sprintf(KEY_TEST,"%s_LOGMASS", PREFIX); 
      if ( strstr(word0,KEY_TEST) != NULL ) {
	parse_plusminus_sntextio(word0, KEY_TEST, &iwd, 
				 &SNDATA.HOSTGAL_LOGMASS_OBS[igal], 
				 &SNDATA.HOSTGAL_LOGMASS_ERR[igal] );
      }
      sprintf(KEY_TEST,"%s_sSFR", PREFIX); 
      if ( strstr(word0,KEY_TEST) != NULL ) {
	parse_plusminus_sntextio(word0, KEY_TEST, &iwd, 
				 &SNDATA.HOSTGAL_sSFR[igal], 
				 &SNDATA.HOSTGAL_sSFR_ERR[igal] );
      }

      sprintf(KEY_TEST,"%s_MAG:", PREFIX); 
      if ( strcmp(word0,KEY_TEST) == 0 ) {
	for(ifilt=0; ifilt < NFILT; ifilt++ ) {
	  iwd++ ; get_PARSE_WORD_FLT(langC, iwd, 
				     &SNDATA.HOSTGAL_MAG[igal][ifilt]); 
	}
      }  
      
    } // end igal      

  } // end HOSTGAL
      
  else if ( strcmp(word0,"PEAKMJD:") == 0 || 
	    strcmp(word0,"SEARCH_PEAKMJD:") == 0 ) {
    iwd++; get_PARSE_WORD_FLT(langC,iwd, &SNDATA.SEARCH_PEAKMJD);    
  }
  else if ( strcmp(word0,"SEARCH_TYPE:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC,iwd, &SNDATA.SEARCH_TYPE );
  }
  
  else if ( IS_PRIVATE ) {
    char key_with_colon[100];
    SNTEXTIO_FILE_INFO.NVAR_PRIVATE_READ++ ;
    NVAR = SNDATA.NVAR_PRIVATE ;
    for(ivar=1; ivar <= NVAR; ivar++ ) {
      sprintf(key_with_colon,"%s:", SNDATA.PRIVATE_KEYWORD[ivar]) ;
      if ( strcmp(word0,key_with_colon) == 0 ) {
	iwd++; get_PARSE_WORD_DBL(langC,iwd,&SNDATA.PRIVATE_VALUE[ivar]);
      }  
    }    
  }

  else if ( strcmp(word0,"NOBS:") == 0 ) {
    iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.NOBS );
    SNDATA.NEPOCH = SNDATA.NOBS; // goofy logic here
  }

  // ---------------------
  // !!! SIM !!! 

  if ( SNDATA.FAKE == FAKEFLAG_LCSIM ) { 

    if ( strcmp(word0,"SIM_MODEL_NAME:") == 0 ) {
      iwd++ ; get_PARSE_WORD(langC, iwd, SNDATA.SIM_MODEL_NAME );
    }
    else if ( strcmp(word0,"SIM_TYPE_NAME:") == 0 ) {
      iwd++ ; get_PARSE_WORD(langC, iwd, SNDATA.SIM_TYPE_NAME );
    }

    if ( strcmp(word0,"SIM_MODEL_INDEX:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_MODEL_INDEX );
    }
    else if ( strcmp(word0,"SIM_TYPE_INDEX:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_TYPE_INDEX );
    }
    else if ( strcmp(word0,"SIM_TEMPLATE_INDEX:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_TEMPLATE_INDEX );
    }
    else if ( strcmp(word0,"SIM_SUBSAMPLE_INDEX:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SUBSAMPLE_INDEX );
    }
    else if ( strcmp(word0,"SIM_LIBID:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_LIBID );
    }
    else if ( strcmp(word0,"SIM_NGEN_LIBID:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_NGEN_LIBID );
    }
    else if ( strcmp(word0,"SIM_NOBS_UNDEFINED:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_NOBS_UNDEFINED );
    }
    else if ( strcmp(word0,"SIM_SEARCHEFF_MASK:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_SEARCHEFF_MASK );
    }
    else if ( strcmp(word0,"SIM_REDSHIFT_HELIO:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_REDSHIFT_HELIO );
    }

    else if ( strcmp(word0,"SIM_REDSHIFT_CMD:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_REDSHIFT_CMB );
    }
    else if ( strcmp(word0,"SIM_REDSHIFT_HOST:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_REDSHIFT_HOST );
    }
    else if ( strcmp(word0,"SIM_REDSHIFT_FLAG:") == 0 ) {
      iwd++ ; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_REDSHIFT_FLAG );
    }
    else if ( strcmp(word0,"SIM_VPEC:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_VPEC );
    }
    else if ( strcmp(word0,"SIM_HOSTLIB_GALID:") == 0 ) {
      iwd++ ; get_PARSE_WORD_DBL(langC, iwd, &DVAL );
      SNDATA.SIM_HOSTLIB_GALID = (long long)DVAL ;
    }
    else if ( strncmp(word0,"SIM_HOSTLIB",11) == 0 ) {
      for(ipar=0; ipar < SNDATA.NPAR_SIM_HOSTLIB; ipar++ ) {
	sprintf(KEY_TEST,"%s:", SNDATA.SIM_HOSTLIB_KEYWORD[ipar]) ;
	if ( strcmp(word0,KEY_TEST) == 0 ) {
	  iwd++ ; get_PARSE_WORD_FLT(langC, iwd, 
				     &SNDATA.SIM_HOSTLIB_PARVAL[ipar]) ; 
	}
      }
    }

    else if ( strcmp(word0,"SIM_DLMU:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_DLMU );
    }
    else if ( strcmp(word0,"SIM_LENSDMU:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_LENSDMU );
    }
    else if ( strcmp(word0,"SIM_RA:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_RA );
    }
    else if ( strcmp(word0,"SIM_DEC:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_DEC );
    }
    else if ( strcmp(word0,"SIM_MWEBV:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_MWEBV );
    }
    else if ( strcmp(word0,"SIM_PEAKMJD:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_PEAKMJD );
    }
    else if ( strcmp(word0,"SIM_MAGSMEAR_COH:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_MAGSMEAR_COH );
    }
    else if ( strcmp(word0,"SIM_AV:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_AV );
    }
    else if ( strcmp(word0,"SIM_RV:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_RV );
    }
    else if ( strncmp(word0,"SIM_SALT2",9) == 0 ) {
      if ( strcmp(word0,"SIM_SALT2x0:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2x0 );
      }
      else if ( strcmp(word0,"SIM_SALT2x1:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2x1 );
      }
      else if ( strcmp(word0,"SIM_SALT2c:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2c );
      }
      else if ( strcmp(word0,"SIM_SALT2mB:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2mB );
      }
      else if ( strcmp(word0,"SIM_SALT2alpha:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2alpha );
      }
      else if ( strcmp(word0,"SIM_SALT2beta:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2beta );
      }
      else if ( strcmp(word0,"SIM_SALT2gammaDM:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_SALT2gammaDM );
      }

    } // end SIM_SALT2
    else if ( strcmp(word0,"SIM_DELTA:") == 0 ) {
	iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_DELTA );
    }
    else if ( strcmp(word0,"SIM_STRETCH:") == 0 ) {
      iwd++ ; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_STRETCH );
    }
    else if ( strncmp(word0,"SIMSED",6) == 0 ) {
      for(ipar=0; ipar < SNDATA.NPAR_SIMSED; ipar++ ) { 
	sprintf(KEY_TEST, "%s:", SNDATA.SIMSED_KEYWORD[ipar]) ;
	if ( strcmp(word0,KEY_TEST) == 0 ) { 
	  iwd++; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIMSED_PARVAL[ipar]);
	}
      }
    }
    
    else if ( strncmp(word0,"SIM_PEAKMAG",11) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_PEAKMAG_%c:", FILTERSTRING[ifilt_obs]);
	if ( strcmp(word0,KEY_TEST) == 0 )  {
	  iwd++; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_PEAKMAG[ifilt]);
	}
      }
    }

    else if ( strncmp(word0,"SIM_TEMPLATEMAG",15) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_TEMPLATEMAG_%c:", FILTERSTRING[ifilt_obs]);
	if ( strcmp(word0,KEY_TEST) == 0 )  {
	  iwd++; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_TEMPLATEMAG[ifilt]);
	}
      }
    }    

    else if ( strncmp(word0,"SIM_EXPOSURE",12) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_EXPOSURE_%c:", FILTERSTRING[ifilt_obs]);
	if ( strcmp(word0,KEY_TEST) == 0 )  {
	  iwd++; get_PARSE_WORD_FLT(langC,iwd,&SNDATA.SIM_EXPOSURE_TIME[ifilt]);
	}
      }
    }

    else if ( strncmp(word0,"SIM_GALFRAC",11) == 0 ) {
      for ( ifilt=0; ifilt < SNDATA_FILTER.NDEF; ifilt++ ) {
	ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
	sprintf(KEY_TEST, "SIM_GALFRAC_%c:", FILTERSTRING[ifilt_obs]);
	if ( strcmp(word0,KEY_TEST) == 0 )  {
	  iwd++; get_PARSE_WORD_FLT(langC, iwd, &SNDATA.SIM_GALFRAC[ifilt]);
	}
      }
    }

    else if ( strncmp(word0,"SIM_STRONGLENS",14) == 0 ) {
      sprintf(PREFIX, "SIM_STRONGLENS") ;

      sprintf(KEY_TEST,"%s_ID:", PREFIX);
      if ( strcmp(word0,KEY_TEST) == 0 ) 
	{ iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_SL_IDLENS); }

      sprintf(KEY_TEST,"%s_z:", PREFIX);
      if ( strcmp(word0,KEY_TEST) == 0 ) 
	{ iwd++; get_PARSE_WORD_DBL(langC, iwd, &SNDATA.SIM_SL_zLENS); }

      sprintf(KEY_TEST,"%s_TDELAY:", PREFIX);
      if ( strcmp(word0,KEY_TEST) == 0 ) 
	{ iwd++; get_PARSE_WORD_DBL(langC, iwd, &SNDATA.SIM_SL_TDELAY); }

      sprintf(KEY_TEST,"%s_MAGSHIFT:", PREFIX);
      if ( strcmp(word0,KEY_TEST) == 0 ) 
	{ iwd++; get_PARSE_WORD_DBL(langC, iwd, &SNDATA.SIM_SL_MAGSHIFT); }

      sprintf(KEY_TEST,"%s_NIMG", PREFIX);
      if ( strcmp(word0,KEY_TEST) == 0 ) 
	{ iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_SL_NIMG); }

      sprintf(KEY_TEST,"%s_IMGNUM", PREFIX);
      if ( strcmp(word0,KEY_TEST) == 0 ) 
	{ iwd++; get_PARSE_WORD_INT(langC, iwd, &SNDATA.SIM_SL_IMGNUM); }
      
    }  // end SIM_STRONGLENS

    else if ( ncmp_PySEDMODEL == 0 && len_PySEDMODEL > 2 ) {
      for(ipar=0; ipar < SNDATA.NPAR_PySEDMODEL; ipar++ ) { 
	sprintf(KEY_TEST, "%s:", SNDATA.PySEDMODEL_KEYWORD[ipar]) ;
	if ( strcmp(word0,KEY_TEST) == 0 ) {
	  iwd++; get_PARSE_WORD_FLT(langC, iwd, 
				    &SNDATA.PySEDMODEL_PARVAL[ipar]); 
	}
      }
    } // end PySEDMODEL

    else if ( strncmp(word0,"LCLIB_PARAM",11) == 0 ) {
      for(ipar=0; ipar < SNDATA.NPAR_LCLIB; ipar++ ) { 
	sprintf(KEY_TEST, "%s:", SNDATA.LCLIB_KEYWORD[ipar]) ;
	if ( strcmp(word0,KEY_TEST) == 0 ) {
	  iwd++; get_PARSE_WORD_FLT(langC, iwd, 
				    &SNDATA.LCLIB_PARVAL[ipar]); 
	}
      }
    }  // end LCLIB

  } // end SIM_LCSIM


  /*
  printf(" xxx %s check word0 = '%s'  cnt=%d\n", 
	 fnam, word0, strncmp(word0,"LCLIB_PARAM",11) );
  */
  

  *iwd_file = iwd ;

  return true ;

} // end parse_SNTEXTIO_HEAD

// ========================================================
void parse_plusminus_sntextio(char *word, char *key, int *iwd_file, 
			      float *PTR_VAL, float *PTR_ERR) {

  // Created Feb 15 2021
  //
  // Parse text file keys of the form
  //   [key]: <value> +- <error>
  //      or
  //   [key]:     <value> 
  //   [key_ERR]: <error>
  //
  // Inputs:
  //   *word   : current word in text file
  //   *key    : desired key name without colon; e.g., 'MWEBV'
  //   *iwd_file : current iwd index while reading text file
  //
  // Outputs:
  //   *iwd_file : final iwd index reading file
  //   *PTR_VAL  : float value associated with key
  //   *PTR_ERR  : float value for error
  //
  // Beware that input *iwd_file is modified here.

  int iwd = *iwd_file;
  char next_word[60], KEY[60], KEY_ERR[60];
  int  langC  =  LANGFLAG_PARSE_WORDS_C ;
  int  lenkey = strlen(key);
  char fnam[] = "parse_plusminus_sntextio" ;

  // ------------ BEGIN --------

  if ( strstr(word,key) == NULL ) { return; }

  // construct key names with colons
  sprintf(KEY,     "%s:",     key);
  sprintf(KEY_ERR, "%s_ERR:", key);

  if ( strcmp(word,KEY) == 0 ) { 

    // read value after KEY
    iwd++; get_PARSE_WORD_FLT(langC, iwd, PTR_VAL );

    // check +- format
    iwd++; get_PARSE_WORD(langC, iwd, next_word );
    if ( strcmp(next_word,"+-") == 0 || strcmp(next_word,"+_") == 0 ) 
      { iwd++; get_PARSE_WORD_FLT(langC, iwd, PTR_ERR ); }
    else 
      { iwd--; }
  }

  // check explicit error key of the form [KEY]_ERR: 
  if ( strcmp(word,KEY_ERR) == 0 ) { 
    iwd++; get_PARSE_WORD_FLT(langC, iwd, PTR_ERR );
  }

  // update iwd pointer in text file
  *iwd_file = iwd;

  return;

} // end parse_plusminus_sntextio


// ================================
void check_head_sntextio(int OPT) {

  // OPT = 0 -> one time global init
  // OPT = 1 -> init for each event
  // OPT = 2 -> do the check

  int i;
  bool HEAD_EXIST ;
  char *CCID = SNDATA.CCID ;
  char *KEY;
  char fnam[] = "check_head_sntextio" ;

  // ------------- BEGIN -------------

  if ( OPT == 0 ) {
    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_SURVEY];
    sprintf(KEY,"SURVEY") ;

    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_SNID];
    sprintf(KEY,"SNID") ;

    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_z];
    sprintf(KEY,"REDSHIFT_FINAL");

    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_FILTERS];
    sprintf(KEY,"FILTERS") ;

    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_FAKE];
    sprintf(KEY,"FAKE") ;

    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_RA];
    sprintf(KEY,"RA") ;

    KEY = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[HEAD_REQUIRE_DEC];
    sprintf(KEY,"DEC") ;
  }

  else if ( OPT == 1 ) {
    // for each event: init exist logicals to false
    for(i=0; i < NHEAD_REQUIRE; i++ )
      { SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[i] = false; }

  } 
  else {
    // if FAKE key not found, assume real data and give warning
    i = HEAD_REQUIRE_FAKE ;
    HEAD_EXIST = SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[i];
    if ( !HEAD_EXIST ) {
      SNDATA.FAKE = FAKEFLAG_DATA;
      sprintf(SNDATA.DATATYPE, "%s", DATATYPE_DATA); 
      SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[i] = true ;       
    }
    
    // abort if any sanity check fails
    int NERR = 0 ;  bool HEAD_EXIST;
    for(i=0; i < NHEAD_REQUIRE; i++ ) {
      HEAD_EXIST = SNTEXTIO_FILE_INFO.HEAD_EXIST_REQUIRE[i];
      KEY        = SNTEXTIO_FILE_INFO.HEAD_KEYNAME_REQUIRE[i];
      if ( !HEAD_EXIST ) {
	NERR++;
	printf("\n ERROR: missing required header key %s:", KEY);     
      }
    } // end i loop over NHEAD
    if ( NERR > 0 ) {
      sprintf(c1err,"Missing %d required header keys (see above)", NERR);
      sprintf(c2err,"Check SNID=%s", CCID );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }

    int NVAR_EXPECT = SNDATA.NVAR_PRIVATE ;
    int NVAR_FOUND  = SNTEXTIO_FILE_INFO.NVAR_PRIVATE_READ ;
    if ( NVAR_FOUND != NVAR_EXPECT ) {
      sprintf(c1err,"Found %d PRIVATE variables in CID=%s", 
	      NVAR_FOUND, CCID );
      sprintf(c2err,"but expected %d from first data file", NVAR_EXPECT);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
    }

  }

  return ;

} // end check_head_sntextio

// =====================================
bool parse_SNTEXTIO_OBS(int *iwd_file) {

  // Created Feb 15 2021
  //
  // Data file has already been read and each word is stored.
  // Here, examine word(*iwd_file) for standard key.
  // If standard key, read argument (next word after key)
  // and load SNDATA struct.
  // Also increment and return *iwd_file.
  //
  // Function returns true to keep reading;
  // returns false when end of header is reached by 
  // finding NOBS key.

  int  langC     = LANGFLAG_PARSE_WORDS_C ;
  int  iwd       = *iwd_file ;
  int  ep, ivar, NVAR = SNTEXTIO_FILE_INFO.NVAROBS ;
  bool DONE_OBS = false ;
  char word0[100], PREFIX[40], KEY_TEST[80], *varName, *str;
  char STRING[40];
  char fnam[] = "parse_SNTEXTIO_OBS";

  // ------------ BEGIN -----------
 
  get_PARSE_WORD(langC, iwd, word0) ;

  if ( strstr(word0,"END")       != NULL ) { DONE_OBS = true; }
  if ( strcmp(word0,"NSPECTRA:") == 0    ) { DONE_OBS = true; }
  if ( DONE_OBS ) { *iwd_file--; return false; }


  if ( strcmp(word0,"OBS:") == 0 ) {
    for(ivar=0; ivar < NVAR; ivar++ ) {
      iwd++ ;   get_PARSE_WORD(langC, iwd, STRING);
      if ( strcmp(STRING,"OBS:") == 0 ) {
	sprintf(c1err,"Found OBS key at ivar=%d of %d (last MJD=%.3f)",
		ivar, NVAR, SNDATA.MJD[ep]);
	sprintf(c2err,"Data file for SNID=%s is messed up.", SNDATA.CCID);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
      }
      sprintf(SNTEXTIO_FILE_INFO.STRING_LIST[ivar],"%s", STRING);
      
    }

    SNTEXTIO_FILE_INFO.NOBS_READ++ ;
    ep = SNTEXTIO_FILE_INFO.NOBS_READ ; // ep starts at 1 for SNDATA struct

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.MJD] ;
    sscanf(str, "%le", &SNDATA.MJD[ep] );
    SNDATA.OBSFLAG_WRITE[ep] = true ; 

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.BAND] ;
    sprintf(SNDATA.FILTCHAR[ep], "%s", str);
    catVarList_with_comma(SNDATA.FILTCHAR_1D,str);

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.FIELD] ;
    sprintf(SNDATA.FIELDNAME[ep], "%s", str);
    catVarList_with_comma(SNDATA.FIELDNAME_1D,str);

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.FLUXCAL] ;
    sscanf(str, "%f", &SNDATA.FLUXCAL[ep] );

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.FLUXCALERR] ;
    sscanf(str, "%f", &SNDATA.FLUXCAL_ERRTOT[ep] );
    
    if ( IVAROBS_SNTEXTIO.PHOTFLAG >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.PHOTFLAG] ;
      sscanf(str, "%d", &SNDATA.PHOTFLAG[ep] );
    }
    if ( IVAROBS_SNTEXTIO.PHOTPROB >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.PHOTPROB] ;
      sscanf(str, "%f", &SNDATA.PHOTPROB[ep] ) ;
    }

    if ( IVAROBS_SNTEXTIO.ZPFLUX >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.ZPFLUX] ;
      sscanf(str, "%f", &SNDATA.ZEROPT[ep] );
    }
    if ( IVAROBS_SNTEXTIO.ZPERR >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.ZPERR] ;
      sscanf(str, "%f", &SNDATA.ZEROPT_ERR[ep] );
    }

    if ( IVAROBS_SNTEXTIO.PSF >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.PSF] ;
      sscanf(str, "%f", &SNDATA.PSF_SIG1[ep]);
    }

    if ( IVAROBS_SNTEXTIO.SKYSIG >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.SKYSIG] ;
      sscanf(str, "%f", &SNDATA.SKY_SIG[ep] );
    }
    if ( IVAROBS_SNTEXTIO.SKYSIG_T >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.SKYSIG_T] ;
      sscanf(str, "%f", &SNDATA.SKY_SIG_T[ep] );
    }

    if ( IVAROBS_SNTEXTIO.GAIN >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.GAIN] ;
      sscanf(str, "%f", &SNDATA.GAIN[ep] );
    }

    if ( IVAROBS_SNTEXTIO.XPIX >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.XPIX] ;
      sscanf(str, "%f", &SNDATA.XPIX[ep] );
    }
    if ( IVAROBS_SNTEXTIO.YPIX >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.YPIX] ;
      sscanf(str, "%f", &SNDATA.YPIX[ep] );
    }

    if ( IVAROBS_SNTEXTIO.CCDNUM >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.CCDNUM] ;
      sscanf(str, "%d", &SNDATA.CCDNUM[ep] );
    }

    // - - -
    if ( IVAROBS_SNTEXTIO.SIMEPOCH_MAG >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVAROBS_SNTEXTIO.SIMEPOCH_MAG] ;
      sscanf(str, "%f", &SNDATA.SIMEPOCH_MAG[ep] );
    }

  }     // end OBS key

  // - - - - -

  *iwd_file = iwd;

  return true ;

} // end parse_SNTEXTIO_OBS


// =====================================
bool parse_SNTEXTIO_SPEC(int *iwd_file) {

  // Created Feb 17 2021
  // Look for SPECTRUM keys, and load GENSPEC struct.

  int  iwd     = *iwd_file ;
  int  langC   = LANGFLAG_PARSE_WORDS_C ;
  int  ID, ivar, NBLAM, ilam ; 
  int  ISPEC   = SNTEXTIO_FILE_INFO.NSPEC_READ-1 ;
  double MJD ;
  char word0[100], *str ;
  char fnam[]  = "parse_SNTEXTIO_SPEC" ;

  // ------------ BEGIN -------------

  get_PARSE_WORD(langC, iwd, word0) ;

  if ( strcmp(word0,"NSPECTRA:") == 0 ) {
      iwd++ ;  get_PARSE_WORD_INT(langC, iwd, &GENSPEC.NMJD_PROC );
  }

  else if ( strcmp(word0,"NVAR_SPEC:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &SNTEXTIO_FILE_INFO.NVARSPEC );
  }

  else if ( strcmp(word0,"VARNAMES_SPEC:") == 0 || 
	    strcmp(word0,"VARLIST_SPEC:" ) == 0  ) {
    rd_sntextio_varlist_spec(&iwd);
  }

  else if ( strcmp(word0,"SPECTRUM_ID:") == 0 ) {
    SNTEXTIO_FILE_INFO.NSPEC_READ++ ;
    SNTEXTIO_FILE_INFO.NLAM_READ = 0 ;
    ISPEC = SNTEXTIO_FILE_INFO.NSPEC_READ - 1 ;
    iwd++; get_PARSE_WORD_INT(langC, iwd, &GENSPEC.ID_LIST[ISPEC] );
  }

  else if ( strcmp(word0,"SPECTRUM_MJD:") == 0 ) {
    iwd++ ;  get_PARSE_WORD_DBL(langC, iwd, &MJD );
    GENSPEC.MJD_LIST[ISPEC] = MJD ;
    if ( MJD > 0.0 ) 
      { GENSPEC.IS_HOST[ISPEC] = 0 ; } // is SN spectrum
    else
      { GENSPEC.IS_HOST[ISPEC] = 1 ; }  // is HOST spectrum
  }

  else if ( strcmp(word0,"SPECTRUM_TEXPOSE:") == 0 ) {
    iwd++; get_PARSE_WORD_DBL(langC, iwd, &GENSPEC.TEXPOSE_LIST[ISPEC] );
  }

  else if ( strcmp(word0,"SPECTRUM_SNR_COMPUTE:") == 0 ) {
    iwd++; get_PARSE_WORD_FLT(langC, iwd, &GENSPEC.SNR_COMPUTE_LIST[ISPEC] );
  }

  else if ( strcmp(word0,"SPECTRUM_LAMOBS_SNR:") == 0 ) {
    iwd++; get_PARSE_WORD_DBL(langC, iwd, &GENSPEC.LAMOBS_SNR_LIST[ISPEC][0]);
    iwd++; get_PARSE_WORD_DBL(langC, iwd, &GENSPEC.LAMOBS_SNR_LIST[ISPEC][1]);
  }

  else if ( strcmp(word0,"SPECTRUM_NLAM:") == 0 ) {
    iwd++; get_PARSE_WORD_INT(langC, iwd, &NBLAM );
    init_GENSPEC_EVENT(ISPEC,NBLAM);    // malloc GENSPEC arrays vs. ilam
    GENSPEC.NBLAM_VALID[ISPEC] = NBLAM ;
  }

  else if ( strcmp(word0,"SPEC:") == 0 ) {
    for ( ivar=0; ivar < SNTEXTIO_FILE_INFO.NVARSPEC ; ivar++ ) {
      iwd++; get_PARSE_WORD(langC, iwd, SNTEXTIO_FILE_INFO.STRING_LIST[ivar] );
    }

    ilam = SNTEXTIO_FILE_INFO.NLAM_READ ; // ilam starts at 0
    SNTEXTIO_FILE_INFO.NLAM_READ++ ;

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVARSPEC_SNTEXTIO.LAMMIN] ;
    sscanf(str, "%le", &GENSPEC.LAMMIN_LIST[ISPEC][ilam]) ;

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVARSPEC_SNTEXTIO.LAMMAX] ;
    sscanf(str, "%le", &GENSPEC.LAMMAX_LIST[ISPEC][ilam]) ;

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVARSPEC_SNTEXTIO.FLAM] ;
    sscanf(str, "%le", &GENSPEC.FLAM_LIST[ISPEC][ilam]) ;

    str = SNTEXTIO_FILE_INFO.STRING_LIST[IVARSPEC_SNTEXTIO.FLAMERR] ;
    sscanf(str, "%le", &GENSPEC.FLAMERR_LIST[ISPEC][ilam]) ;

    if ( IVARSPEC_SNTEXTIO.SIM_GENFLAM >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVARSPEC_SNTEXTIO.SIM_GENFLAM] ;
      sscanf(str, "%le", &GENSPEC.GENFLAM_LIST[ISPEC][ilam]) ;
    }
    if ( IVARSPEC_SNTEXTIO.SIM_GENMAG >= 0 ) {
      str = SNTEXTIO_FILE_INFO.STRING_LIST[IVARSPEC_SNTEXTIO.SIM_GENMAG] ;
      sscanf(str, "%le", &GENSPEC.GENMAG_LIST[ISPEC][ilam]) ;
    }
    //
  }


  *iwd_file = iwd;
  return true ;

} // end parse_SNTEXTIO_SPEC


