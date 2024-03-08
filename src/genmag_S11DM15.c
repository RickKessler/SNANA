/***********************************************
  Created Aug 2013 by R.Kessler


  Compute DM15/stretch model used in Sako et al, [2008, 2011, 2013]
  for photometric classification. (S08,S11,S13).
  The SNIa mag-model is from Appendic C of S08 and is used in the 
  succeeding publications. The mag-error model is from Eq. 6 of S11, 
  and the 8 templates (not part of this code) are also from S11. 
  Hence we use 'S11' as a generic reference to this model that 
  includes mags, errors and templates.
 
  Note that host RV=3.1 in S08, and RV=2.2 in S11;
  you pick the value with sim-input key  
         GENRANGE_RV: <RVlo> <RVhi>

  Apr 12 2017: call malloc_S11DM15_SEDMODEL().

************************************************/

/*
#include <stdio.h> 
#include <math.h>     // log10, pow, ceil, floor
#include <stdlib.h>   // includes exit(),atof()
*/

#include "sntools.h"           // community tools
#include "MWgaldust.h"
#include "genmag_SEDtools.h"
#include "genmag_S11DM15.h" 


// ================= BEGIN FUNCTIONS ========================

int init_genmag_S11DM15(char *VERSION, int OPTMASK ) {

  // Aug 2013
  // initialize S11DM15 model; read INFO file and SED.
  // This function is largely copied from init_genamg_SALT2.
  //
  // Feb 22 2017: allow input *VERSOIN to include path

  int retval = 0 ;     // return 0 if no errors
  int NDAY, NLAM, NZBIN ;

  char  fnam[] = "init_genmag_S11DM15"  ;
  char  BANNER[120], version[60] ;

  // -------------- BEGIN -------------

  sprintf(BANNER, "%s : Initialize %s", fnam, version );
  print_banner(BANNER) ;

  sprintf(S11DM15_INFO_FILE, "S11DM15.INFO" );

  if ( NFILT_SEDMODEL == 0 ) {
    sprintf(c1err,"No filters defined ?!?!?!? " );
    sprintf(c2err,"Need to call init_filter_SEDMODEL");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // summarize filter info
  filtdump_SEDMODEL();

    
  // ==========================================
  // construct path to model area
  
  extract_MODELNAME(VERSION,
		    S11DM15_MODELPATH, version); // returned

  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    sprintf( S11DM15_MODELPATH, "%s/%s", 
	     getenv(PRIVATE_MODELPATH_NAME), version );
  }
  else if ( strlen(S11DM15_MODELPATH) > 0 ) {
    // do nothing for user-defined path/model
  }
  else {
    // default location under $SNDATA_ROOT
    sprintf( S11DM15_MODELPATH, "%s/models/S11DM15/%s", 
	     getenv("SNDATA_ROOT"), version );
  }
  
  // misc inits
  SEDMODEL.NSURFACE   =  1   ;
  SEDMODEL_MWEBV_LAST = -999. ;

  RVMW_S11DM15        =  3.1 ;  // beware hard-wire
  S11DM15_LAST.z      = -9.0 ; 
  S11DM15_LAST.dm15   = -9.0 ; 
  S11DM15_LAST.AV     = -9.0 ;
  S11DM15_LAST.RV     = -9.0 ;
  
  // read info file.
  read_S11DM15_INFO_FILE();

  // set extreme ranges to read any SED
  int     nflux_nan;
  double  Trange[2], Lrange[2] ;
  char    tmpFile[MXPATHLEN], sedComment[60];

  Trange[0] = -20. ;  // rest-frame days
  Trange[1] = 200. ;
  Lrange[0] = 100. ;  // rest-frame wavelength (A)
  Lrange[1] = 30000. ;

  sprintf(tmpFile, "%s/SED.DAT", S11DM15_MODELPATH );
  sprintf(sedComment,"S11DM15");

  // Apr 12 2017: fix old refactor and malloc S11DM15_SEDMODEL
  NZBIN            = NZBIN_SEDMODEL_DEFAULT ;
  init_redshift_SEDMODEL(NZBIN,
			 ZMIN_SEDMODEL_DEFAULT, ZMAX_SEDMODEL_DEFAULT);
  malloc_S11DM15_SEDMODEL();   
  
  rd_sedFlux(tmpFile, sedComment, Trange, Lrange
	     ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0
	     ,&S11DM15_SEDMODEL.NDAY
	     ,S11DM15_SEDMODEL.DAY, &S11DM15_SEDMODEL.DAYSTEP
	     ,&S11DM15_SEDMODEL.NLAM
	     ,S11DM15_SEDMODEL.LAM, &S11DM15_SEDMODEL.LAMSTEP
	     ,S11DM15_SEDMODEL.FLUX, S11DM15_SEDMODEL.FLUXERR
	     ,&nflux_nan );
  
  // store useful quantities
  NDAY = S11DM15_SEDMODEL.NDAY ;
  NLAM = S11DM15_SEDMODEL.NLAM ;
  S11DM15_SEDMODEL.LAMMIN = S11DM15_SEDMODEL.LAM[0] ;
  S11DM15_SEDMODEL.LAMMAX = S11DM15_SEDMODEL.LAM[NLAM-1] ;
  S11DM15_SEDMODEL.DAYMIN = S11DM15_SEDMODEL.DAY[0] ;
  S11DM15_SEDMODEL.DAYMAX = S11DM15_SEDMODEL.DAY[NDAY-1] ;

  // store SEDMODL struct values that are used for error checking 
  SEDMODEL.LAMMIN_ALL = S11DM15_SEDMODEL.LAMMIN ;
  SEDMODEL.LAMMAX_ALL = S11DM15_SEDMODEL.LAMMAX ;

  printf("\t LAM(MIN,MAX,STEP)= %.0f, %.0f, %.0f \n"
	 ,S11DM15_SEDMODEL.LAMMIN
	 ,S11DM15_SEDMODEL.LAMMAX
	 ,S11DM15_SEDMODEL.LAMSTEP );

  printf("\t DAY(MIN,MAX,STEP)= %.0f, %.0f, %.1f \n"
	 ,S11DM15_SEDMODEL.DAYMIN
	 ,S11DM15_SEDMODEL.DAYMAX
	 ,S11DM15_SEDMODEL.DAYSTEP );

  fflush(stdout);

  // find Trest(B band max) of template
  Trest_Bmax_S11DM15 = get_Tmax_S11DM15();
  
  // print stretch vs. dm15 and extinction vs. lambda
  dumpStuff_S11DM15();

  return retval ;

} // end of init_genmag_S11DM15


// =====================================
void malloc_S11DM15_SEDMODEL(void) {

  // Created Apr 12 2017 (modified from malloc_TEMP_SEDMODEL)
  
  int MEMD = sizeof(double);
  S11DM15_SEDMODEL.N_FINEBIN = 0 ; 
  S11DM15_SEDMODEL.DAY     = (double*) malloc ( MEMD * MXBIN_DAYSED_SEDMODEL );
  S11DM15_SEDMODEL.LAM     = (double*) malloc ( MEMD * MXBIN_LAMSED_SEDMODEL );
  S11DM15_SEDMODEL.FLUX    = (double*) malloc ( MEMD * MXBIN_SED_SEDMODEL );
  S11DM15_SEDMODEL.FLUXERR = (double*) malloc ( MEMD * MXBIN_SED_SEDMODEL );

  return ;
}

// =============================================
void read_S11DM15_INFO_FILE(void) {

  // read info file and fill S11DM15.INFO structre

  char 
    infoFile[MXPATHLEN]
    ,c_get[60]
    ,fnam[] = "read_S11DM15_INFO_FILE"
    ;

  int i ;

  FILE *fp ;

  // --------- BEGIN ----------

  sprintf(infoFile, "%s/%s", S11DM15_MODELPATH, S11DM15_INFO_FILE );

  if (( fp = fopen(infoFile, "rt")) == NULL ) {
    sprintf(c1err,"Could not open info file:");
    sprintf(c2err," %s", infoFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  printf("\n  Read S11DM15 model parameters from \n  %s\n", infoFile );

  S11DM15_INFO.DM15_REF = -9.0 ;
  for(i=0; i<=2; i++) { S11DM15_INFO.a_LAMPAR[i]        = -9.0 ; }
  for(i=0; i<=3; i++) { S11DM15_INFO.TAUPOLY_STRETCH[i] = -9.0 ; }
  S11DM15_INFO.b_PAR = -9.0 ;
  S11DM15_INFO.RESTLAMBDA_RANGE[0] = -9.0 ;
  S11DM15_INFO.RESTLAMBDA_RANGE[1] = -9.0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get, "MBPEAK:") == 0 ) 
      { readdouble(fp, 1, &S11DM15_INFO.MBPEAK );  }

    if ( strcmp(c_get, "DM15_REF:") == 0 ) 
      { readdouble(fp, 1, &S11DM15_INFO.DM15_REF );  }

    if ( strcmp(c_get, "a_LAMPAR:") == 0 ) 
      { readdouble(fp, 3, S11DM15_INFO.a_LAMPAR );  }
    if ( strcmp(c_get, "b_PAR:") == 0 ) 
      { readdouble(fp, 1, &S11DM15_INFO.b_PAR );  }

    if ( strcmp(c_get, "TAUPOLY_STRETCH:") == 0 ) 
      { readdouble(fp, 4, S11DM15_INFO.TAUPOLY_STRETCH );  }

    if ( strcmp(c_get, "RESTLAMBDA_RANGE:") == 0 ) 
      { readdouble(fp, 2, S11DM15_INFO.RESTLAMBDA_RANGE );  }

  } // end of fscanf

  fclose(fp);


  // print INFO to screen

  printf("  S11DM15.INFO \n");

  printf("\t RESTLAMBDA_RANGE:  %6.0f - %6.0f A\n"
	 ,S11DM15_INFO.RESTLAMBDA_RANGE[0]
	 ,S11DM15_INFO.RESTLAMBDA_RANGE[1] );

  printf("\t X = DM15 - %.2f \n", S11DM15_INFO.DM15_REF);

  printf("\t a(lam) = ") ;
  for(i=0; i<=2; i++) { printf("%f*LAM^%d  ", S11DM15_INFO.a_LAMPAR[i], i); }
  printf("\n");

  printf("\t b_PAR = %f \n", S11DM15_INFO.b_PAR );

  printf("\t TAUPOLY_STRETCH = ");
  for(i=0; i<=3; i++) { printf("%le  ", S11DM15_INFO.TAUPOLY_STRETCH[i] ); }
  printf("\n");

  fflush(stdout);

} // end of read_S11DM15_INFO_FILE


// ========================================
void  dumpStuff_S11DM15(void) {

  // print summary table of
  // - stretch vs. dm15
  // - color law vs. lambda
  // These tables are just for visual checking and
  // are not used anywhere else.

  double dm15, s, lam, AV, RV, arg, XT_MAG, XT_FRAC;

  // ------- BEGIN ---------

  printf("\n");

  for(dm15 = 0.6; dm15<=1.6; dm15+= 0.2 ) {
    s = 15.0/get_tau_S11DM15(dm15) ;
    printf("\t stretch(dm15=%3.1f) = %.3f \n", dm15, s);
  }


  AV = 1.0 ;
  RV = 3.1 ;
  printf("\n  Host Transmission with RV=%3.1f and AV=%3.1f : \n", RV, AV);

  for ( lam=1000.0; lam <= 10000.0 ; lam += 1000.0 ) {
    XT_MAG   = GALextinct ( RV, AV, lam, 89 ); // host XT
    arg      = -0.4*XT_MAG ;
    XT_FRAC  = pow(TEN,arg);    // flux-fraction thru host
    printf("\t Host-Trans(Lam=%7.1f) = %.5f \n", lam, XT_FRAC );
  }

  fflush(stdout);

} // end of dumpStuff_S11DM15


// ==========================================
void genmag_S11DM15( int ifilt_obs  // (I) absolute filter index
		     ,double dm15    // (I) Phillips dm15 parameter
		     ,double AV      // (I) host extiction in V band
		     ,double RV      // (I) RV parameter for CCM89 
		     ,double mwebv   // (I) Galactic E(B-V)
		     ,double z       // (I) redshift
		     ,double mu      // (I) distance modulus
		     ,int Nobs       // (I) Number of obs passed below
		     ,double *Tobs_list    // (I) List of Tobs
		     ,double *magobs_list  // (O) return list of mags
		     ,double *magerr_list  // (O) return list of mag errors
		     ) {


  int    epobs, ifilt ;
  double z1, Tobs, Trest, Tobs_template, s, Flux, ZP ;
  char   fnam[] = "genmag_S11DM15" ;

  // ----------------- BEGIN ----------------

  z1    = 1.0 + z ;
  ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
  ZP    = FILTER_SEDMODEL[ifilt].ZP ;

  // make sure lambda range is OK for this filter
  checkLamRange_SEDMODEL(ifilt,z,fnam) ;


  // store info for Galactic extinction (genmag_SEDtools.c)
  fill_TABLE_MWXT_SEDMODEL(RVMW_S11DM15,mwebv);

  // store product of flux-scale and host extinction vs. lambda
  fill_TABLE_FLUXSCALE_S11DM15(dm15,RV,AV,z);

  // compute rest-frame stretch as in Appendic C of S08
  s = 15.0/get_tau_S11DM15(dm15) ;

  for ( epobs=0; epobs < Nobs; epobs++ ) {

    magobs_list[epobs] = 99.0 ;  // init value
    magerr_list[epobs] =  9.0 ;  

    Tobs  = Tobs_list[epobs] ;

    // we want Tobs above to be relative to template Tmax,
    // so make internal shift here
    Tobs_template = Tobs  - (z1*Trest_Bmax_S11DM15) ;
    Trest         = Tobs_template/z1 ;

    // bail if outside model range
    if ( Trest/s < S11DM15_SEDMODEL.DAYMIN+.0001 ) { continue ; }
    if ( Trest/s > S11DM15_SEDMODEL.DAYMAX-.0001 ) { continue ; }

    // integrate over filter to get obs-frame Flux
    Flux = INTEG_SEDFLUX_S11DM15(ifilt_obs, z, Tobs_template/s );

    magobs_list[epobs] = 
      -2.5*log10(Flux)       // integrated flux over filte trans
      + ZP                   // zeroppint
      + mu                   // distance modulus
      + S11DM15_INFO.MBPEAK  // normalization
      - 0.2                  // fudge to avoid internal Bpeak calc/norm.
      ;
    
    magerr_list[epobs] = get_magerr_S11DM15(Trest);

  } // epobs


} // end of genmag_S11DM15

// ============================================
void  fill_TABLE_FLUXSCALE_S11DM15(double dm15, double RV, 
				   double AV, double z) {

  // Aug 24 2013
  // store FLUXSCALE_S11DM15 = FLUXSCALE * XTHOST vs. wavelength 
  // for each filter to speed integrations. Note that FLUXSCALE
  // is computed for rest-frame LAMREST, but it is stored vs.
  // obs-frame lambda index for filters. Note use of 1+z to
  // translate between frames.

  int  NLAMFILT, ilam, I8, I8p, ifilt, MEM ;

  double 
    LAMOBS, LAMREST
    ,XT_MAG
    ,XT_FRAC
    ,arg, z1, Fscale
    ;

  //  char fnam[] = "fill_TABLE_XTMW_SEDMODEL";
  
  // ------------- BEGIN ------------------

  I8  = sizeof(double) ;
  I8p = sizeof(double*) ;
  
  // allocate memory for each filter only once
  if ( S11DM15_LAST.z < 0.0 ) {

    printf("\t Malloc S11DM15 TABLES for FLUXSCALE & XTHOST \n");
    fflush(stdout);

    MEM = I8p * (NFILT_SEDMODEL+1) ;
    S11DM15_TABLE_FLUXSCALE   = (double**)malloc(MEM); 
    S11DM15_TABLE_XTHOST_FRAC = (double**)malloc(MEM); 

    for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
      NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
      MEM       = I8 * NLAMFILT ;
      S11DM15_TABLE_FLUXSCALE[ifilt]   = (double*)malloc(MEM); 
      S11DM15_TABLE_XTHOST_FRAC[ifilt] = (double*)malloc(MEM); 
    }
  }


  if ( z    != S11DM15_LAST.z    ) { goto UPDATE ; }
  if ( dm15 != S11DM15_LAST.dm15 ) { goto UPDATE ; }
  if(  AV   != S11DM15_LAST.AV   ) { goto UPDATE ; }
  if(  RV   != S11DM15_LAST.RV   ) { goto UPDATE ; }
  return ; 

 UPDATE:
  z1 = 1.0 + z;

  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
    NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;

    for ( ilam=0; ilam < NLAMFILT; ilam++ ) {
      LAMOBS     = FILTER_SEDMODEL[ifilt].lam[ilam] ;
      LAMREST    = LAMOBS/z1 ;

      XT_MAG   = GALextinct ( RV, AV, LAMREST, 89 ); // host XT
      arg      = -0.4*XT_MAG ;
      XT_FRAC  = pow(TEN,arg);    // flux-fraction thru MW

      Fscale   = get_FluxScale_S11DM15(LAMREST,dm15);

      S11DM15_TABLE_XTHOST_FRAC[ifilt][ilam]  = XT_FRAC ;
      S11DM15_TABLE_FLUXSCALE[ifilt][ilam]    = Fscale ;

    } // ilam

  } // ifilt
  
  
  S11DM15_LAST.z    = z    ;
  S11DM15_LAST.dm15 = dm15 ;
  S11DM15_LAST.AV   = AV   ;
  S11DM15_LAST.RV   = RV   ;

} // end of fill_TABLE_FLUXSCALE_S11DM15


// =================================================
double  INTEG_SEDFLUX_S11DM15(int ifilt_obs, double z, double Tobs ) {

  //  Extract rest-frame spectrum (Flux vs. wavelength) from
  //  SED vs. Trest, using the passed Trest value. Also apply
  //  the flux scale vs. wavelength. The returned *flux vs. alambda
  //  is ready to integrate to get the observed flux.

  double z1, Trest, FSCALE_REST, FSCALE_TOT ;
  double DAYSTEP, DAYMIN, DAYDIF ;
  double TRANS, LAMOBS, LAMSED, LAMDIF, LAMSTEP ;
  double FRAC_INTERP_DAY, FRAC_INTERP_LAMSED, MWXT_FRAC, XTHOST_FRAC ;
  double VAL0, VAL1, FTMP, FDIF, FSED[2], FLUX_INTEG ;

  int    ifilt, IDAY, NDAY, NLAMFILT, ilamobs, ilamsed, iday ;
  int    LABORT, NLAMSED, jflux ;

  char *cfilt ;
  char fnam[] = "INTEG_SEDFLUX_S11DM15";

  int LDMP ;

  // ------------- BEGIN ----------------

  LDMP = (ifilt_obs == -3 && fabs(Tobs) < 10.0);

  FLUX_INTEG = 0.0 ; // output 

  ifilt     = IFILTMAP_SEDMODEL[ifilt_obs] ;
  cfilt     = FILTER_SEDMODEL[ifilt].name ;
  z1        = 1.0 + z ;
  Trest     = Tobs/z1 ;
  NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
  NLAMSED   = S11DM15_SEDMODEL.NLAM ;
  LAMSTEP   = S11DM15_SEDMODEL.LAMSTEP ;

  // S11DM15_SEDMODEL
  DAYSTEP = S11DM15_SEDMODEL.DAYSTEP ;  
  DAYMIN  = S11DM15_SEDMODEL.DAY[0]  ;
  DAYDIF  = Trest - DAYMIN ;
  IDAY    = (int)(DAYDIF/DAYSTEP);  
  NDAY    = S11DM15_SEDMODEL.NDAY ;
  DAYDIF  = Trest - S11DM15_SEDMODEL.DAY[IDAY] ;
  FRAC_INTERP_DAY = DAYDIF/DAYSTEP ;  

  if ( LDMP ) {
    printf("\n\n xxx ########################################## \n");
    printf(" xxx      %s DUMP-%s  \n", fnam, cfilt );
    printf(" xxx Tobs = %f  z=%f \n",  Tobs, z );
    printf(" xxx IDAY = %d/%d   FRAC_INTERP_DAY=%f\n", 
	   IDAY, NDAY, FRAC_INTERP_DAY );
  }

  // -----------------------------------

  for ( ilamobs=0; ilamobs < NLAMFILT; ilamobs++ ) {

    TRANS  = FILTER_SEDMODEL[ifilt].transSN[ilamobs] ;
    if ( TRANS < 1.0E-12 ) { continue ; } 

    LAMOBS     = FILTER_SEDMODEL[ifilt].lam[ilamobs] ;
    LAMSED     = LAMOBS / z1 ;  // rest-frame lambda

    MWXT_FRAC   = SEDMODEL_TABLE_MWXT_FRAC[ifilt][ilamobs] ;  // Galactic XT
    XTHOST_FRAC = S11DM15_TABLE_XTHOST_FRAC[ifilt][ilamobs];  // host XT
    FSCALE_REST = S11DM15_TABLE_FLUXSCALE[ifilt][ilamobs] ;   // model scale
    FSCALE_TOT  = MWXT_FRAC * XTHOST_FRAC * FSCALE_REST ;     // product

    // get rest-frame lambda index for SED space
    LAMDIF  = LAMSED - S11DM15_SEDMODEL.LAMMIN ;
    ilamsed = (int)(LAMDIF/LAMSTEP);
    LAMDIF  = LAMSED - S11DM15_SEDMODEL.LAM[ilamsed] ;
    FRAC_INTERP_LAMSED = LAMDIF / LAMSTEP ; // 0-1

    LABORT = ( FRAC_INTERP_LAMSED < -1.0E-8 || 
	       FRAC_INTERP_LAMSED > 1.0000000001 ) ;
    if ( LABORT ) { 
      sprintf(c1err,"Invalid FRAC_INTERP_LAMSED=%le ", 
	      FRAC_INTERP_LAMSED );
      sprintf(c2err,"check Tobs(%s)=%6.2f at z=%5.3f ",
	      cfilt, Tobs, z );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    
    // interpolate SED Fluxes to LAMSED for each day-grid
    for ( iday=0; iday<2; iday++ ) {
      jflux  = NLAMSED*(IDAY+iday) + ilamsed ;
      VAL0   = S11DM15_SEDMODEL.FLUX[jflux];
      VAL1   = S11DM15_SEDMODEL.FLUX[jflux+1];
      FSED[iday] = VAL0 + (VAL1-VAL0)*FRAC_INTERP_LAMSED ;
    } // iday
      
    // interpolate DAY
    FDIF = FSED[1] - FSED[0] ;
    FTMP = FSED[0] + FDIF * FRAC_INTERP_DAY ; 

    // increment flux
    FLUX_INTEG += (FTMP * LAMSED * TRANS * FSCALE_TOT );

    if ( LDMP && ilamobs < 90 && ilamobs > 80 ) {
      printf(" xxx ----------------------------------------- \n");
      printf(" xxx ilamobs = %3d   LAMOBS=%.0f   LAMSED=%.0f  TRANS=%.4f \n", 
	     ilamobs, LAMOBS, LAMSED, TRANS );
      printf(" xxx SCALE(XTMW*XTHOST*DM15 = TOT) = %f x %f x %f = %f \n",
	     MWXT_FRAC, XTHOST_FRAC, FSCALE_REST, FSCALE_TOT );
      printf(" xxx FSED = %le, %le   FTMP=%le\n", FSED[0], FSED[1], FTMP );
    }

  } // ilamobs

  
  if ( LDMP ) { fflush(stdout); debugexit(fnam); }

  // apply final normalization:
  double FNORM, hc8;
  hc8   = (double)hc ;
  FNORM = LAMSTEP  / hc8 ;
  FLUX_INTEG *= FNORM ;

  return FLUX_INTEG ;

} // end of INTEG_SEDFLUX_S11DM15


// =======================================================
double get_Tmax_S11DM15(void) {

  // loop over Trest(template) and find max B flux.
  // Get B-filter transmission from utility filterTrans_BessB(lam)
  // in genmag_SEDtools.c.

  double Trest, Tmax, BFLUX, BFLUXMAX  ;
  //  char fnam[] = "get_Tmax_S11DM15" ;

  // ------------- BEGIN ---------

  Tmax     = -999.0 ;
  BFLUXMAX = -999.0 ;

  for ( Trest = -5.0; Trest <= 5.0; Trest+=1.0 ) {

    BFLUX = INTEG_BFLUX_S11DM15(Trest);
    
    if ( BFLUX > BFLUXMAX ) {
      Tmax = Trest ;
      BFLUXMAX = BFLUX ;
    }
    //    printf(" xxxx Trest = %.1f -> BFLUX = %le \n", Trest, BFLUX);
  }

  printf("\t Max B-flux at Trest = %f days \n", Tmax); 

  return Tmax ;

} // end of get_Tmax_S11DM15

// ===================================================
double INTEG_BFLUX_S11DM15(double Trest) {

  // Seperate function to integrate rest-frame B-flux
  // Cannot use generic INTEG_S11DM15_SEDFLUX because
  // the Bessell-B is not necessarily defined.

  int    jflux, NLAMSED, ilam, iday ;
  double FLUX, TRANS, FTMP, LAMSED, LAMDIF, DAYDIF ;
  double DAYSTEP, LAMSTEP ;

  // -------------- BEGIN ------------
  FLUX = 0.0 ;

  NLAMSED   = S11DM15_SEDMODEL.NLAM ;
  DAYSTEP   = S11DM15_SEDMODEL.DAYSTEP ;  
  LAMSTEP   = S11DM15_SEDMODEL.LAMSTEP ;

  DAYDIF    = Trest - S11DM15_SEDMODEL.DAYMIN ; 
  iday      = (int)( (DAYDIF+0.0001)/DAYSTEP);  

  for(LAMSED=3600.0 ; LAMSED <= 5600.0 ; LAMSED += 10.0 ) {

    LAMDIF  = LAMSED - S11DM15_SEDMODEL.LAMMIN ;
    ilam    = (int)( (LAMDIF+0.001)/LAMSTEP);

    TRANS  = filterTrans_BessB(LAMSED);
    jflux  = NLAMSED*iday + ilam ;
    FTMP   = S11DM15_SEDMODEL.FLUX[jflux];
    FLUX  += (FTMP * LAMSED * TRANS );
  }

  return FLUX ;

} // end of INTEG_BFLUX_S11DM15

// =======================================================
double get_FluxScale_S11DM15(double lambda, double dm15) {

  // Returns 10**[-0.4*(a*x + b*x^2)],
  // Eq C3 in Sako 2008

  int i ;
  double arg_a, arg_b, arg, alam, b, atmp, X , lamPow ;

  // ----------- begin ------------

  X = dm15 - S11DM15_INFO.DM15_REF ;

  alam = 0.0 ;
  for(i=0; i<3; i++ ) {
    lamPow = pow( lambda, (double)i );
    atmp   = S11DM15_INFO.a_LAMPAR[i];
    alam  += (atmp * lamPow) ;
  }


  b     = S11DM15_INFO.b_PAR ;
  arg_a = alam *  X ;
  arg_b = b    * (X*X) ;
  arg   = -0.4 * ( arg_a + arg_b );

  return  pow(TEN,arg) ;

} // end of FluxScale_S11DM15


// ====================================
double get_tau_S11DM15(double dm15) {

  // return Eq C4 of Sako 2008

  int i ;
  double pow_dm15, tau, coeff ;

  tau = 0.0 ;
  for(i=0; i<4; i++ ) {
    pow_dm15 = pow( dm15, (double)i ) ;
    coeff    = S11DM15_INFO.TAUPOLY_STRETCH[i] ;
    tau     += (coeff * pow_dm15);    
  }

  return tau ;

} // end of tau_S11DM15


// ============================================
double get_magerr_S11DM15(double Trest) {

  // Return hard-wired error From Eq. 6 in S11:
  // magerr = 
  //    0.08 + 0.04 × (|t|/20) |t| < 20 days,
  //    0.12 + 0.08 × ((|t| − 20)/60) |t| ≥ 20 days.
  //
  // ----------- BEGIN -------------

  double magerr, TT ;

  TT = fabs(Trest);

  if ( TT < 20.0 ) {
    magerr = 0.08 + (0.04 * TT) ;
  }
  else {
    magerr = 0.12 + 0.08 * (TT-20.0)/60.0 ;
  }

  return magerr ;

} // end of get_magerr_S11DM15
