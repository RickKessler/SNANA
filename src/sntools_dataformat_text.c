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

#include "sntools.h"
#include "sntools_dataformat_text.h"
#include "sntools_host.h" 
#include "sntools_trigger.h" 
#include "sntools_spectrograph.h"


void WR_DATAFILE_TEXT(void) {

  // Createed Feb 2021
  // Driver function to write single data file for single event.
  //
  // simFlag bits defined in 'grep WRITE_MASK sndata.h'

  char OUTFILE[MXPATHLEN];
  FILE *fp ;
  char fnam[] = "WR_DATAFILE_TEXT" ;

  // ----------- BEGIN ------------

  sprintf(OUTFILE, "%s", SNDATA.SNFILE_OUTPUT );

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

} // end WR_DATAFILE_TEXT


// =====================================================
void  wr_dataformat_text_HEADER(FILE *fp) {

  char comment[80];

  char COMMENT_FAKEFLAG[4][80] = {
    "data",
    "fakes overlaid on images",
    "simulated with snlc_sim.exe" ,
    "BLIND-TEST simulation"
  } ;

  char fnam[] = "wr_dataformat_text_HEADER" ;

  // ------------ BEGIN -----------

  fprintf(fp,"SURVEY:   %s\n", SNDATA.SURVEY_NAME);
  fprintf(fp,"SNID:     %s\n", SNDATA.CCID);
  fprintf(fp,"IAUC:     %s\n", SNDATA.IAUC_NAME);
  fprintf(fp,"SNTYPE:   %d\n", SNDATA.SEARCH_TYPE);
  fprintf(fp,"RA:       %.6f  # deg\n", SNDATA.RA);
  fprintf(fp,"DEC:      %.6f  # deg\n", SNDATA.DEC);
  fprintf(fp,"FILTERS:  %s\n", SNDATA_FILTER.LIST);
  fprintf(fp,"PIXSIZE:  %.4f  # arcsec \n", SNDATA.PIXSIZE);
  fprintf(fp,"CCDNUM:   %d   \n", SNDATA.CCDNUM[0] );

  int FAKE = SNDATA.FAKE ;
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

  fprintf(fp, "SIM_VPEC:            %.1f (km/sec) \n", 
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


  //.xyz
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
      fprintf(fp," %.2f",SNDATA.HOSTGAL_SB_FLUX[ifilt] ) ;
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
      fprintf(fp," %6.2f",SNDATA.HOSTGAL_SB_FLUXERR[ifilt_obs] ) ;
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
  bool WRFLAG_TRIGGER    = (SNDATA.MJD_TRIGGER < 0.99E6);

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
  int  NMJD       = GENSPEC.NMJD_TOT ;
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

  //.xyz

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

    fprintf(fp,"SPECTRUM_NLAM:     %d (of %d)         "
            "# Number of valid wave bins\n",
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




