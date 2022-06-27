/*******************************************
 Created Oct 14, 2009 by R.Kessler


 Generic strutures and tools to use for SED-based models
 such as SALT2 or SN-explosion models. Includes storage
 of primary reference, filters, and SN SED, along with
 zeropoint calculations.

 calling sequence from external function (C or fortran):

   init_primary_SEDMODEL(name, ...)
   init_filter_SEDMODEL(ifilt_obs, ...) ; // repeat for each filter
 
 Then from init_genmag_XXXX, 

   rd_sedFlux(...)
   init_flux_SEDMODEL(ifilt_obs,ISED) ; // for each filter, ISED
   
  WARNING: SEDMODEL.TEMPFLUX is over-written for each SED to save memory.
          Only the filter flux-integrals are stored for each SED.

 Finally, from genmag_XXX

   interp_flux_SEDMODEL(ISED, ifilt_obs, ... ) 


            HISTORY

  July 17 2016: new function  fill_TABLE_HOSTXT_SEDMODEL(...)
                to store host galaxy extinction.

  Jul 30 2016: add SPECTROGRAPH functions

  Jan 10 2017: fix init_primary_SEDMODEL() to return(SUCCESS)

  Apr 12 2018: write FLUX*FLUXSCALE to avoid floating overflow
               whe FLUX ~ E39 and FLUSCALE ~ E-41

  May 2018: new function UVLAM_EXTRAPFLUX_SEDMODEL

  Jul 20, 2018: 
    + new function get_flux_SEDMODEL() to call interp_flux_SEDMODEL
      and take care of extrapolating Trest outside model range.
   
  Aug 23 2019:
    in interp_primaryFlux_SEDMODEL(), add lam-edge protection.

  Jul 13 2020:
    + fix dumb logic but in fill_TABLE_HOSTXT_SEDMODEL() so that
      extinction table is no longer created for every iteration.
      Fits now go almost x10 faster ... same speed as in March 2020

********************************************/

#include "sntools.h"           // community tools
#include "sntools_spectrograph.h"
#include "genmag_SEDtools.h"   // SED tools
#include "MWgaldust.h"         // GALextinct is here

// ******************************
int reset_SEDMODEL(void) {

  // Nov 10, 2010: do one-time inits
  int ifilt;

  SEDMODEL.NSURFACE        =  0 ;
  SEDMODEL.FLUXSCALE       = -9.0 ; // require user to set this later
  NVAR_FIRSTBIN_SEDMODEL   =  0 ;
  NFILT_SEDMODEL           =  0 ;
  FILTLIST_SEDMODEL[0]     =  0 ;
  ISED_SEDMODEL            = -9 ; // Mar 6 2017

  SEDMODEL.DAYMIN_ALL = +9999999.0 ;
  SEDMODEL.DAYMAX_ALL = -9999999.0 ;

  for(ifilt=0; ifilt < MXFILT_SEDMODEL; ifilt++ ) {
    FILTER_SEDMODEL[ifilt].name[0]    = 0;  // Nov 2020
    FILTER_SEDMODEL[ifilt].survey[0]  = 0;  // Nov 2020
    IFILTMAP_SEDMODEL[ifilt]          = -9 ;
    FILTER_SEDMODEL[ifilt].ifilt_obs  = -9 ;
    FILTER_SEDMODEL[ifilt].magprimary = 0.0 ;
    FILTER_SEDMODEL[ifilt].lamshift   = 0.0 ;
  }

  // set default redshift range and NZBIN
  init_redshift_SEDMODEL(NZBIN_SEDMODEL_DEFAULT, 
			 ZMIN_SEDMODEL_DEFAULT, ZMAX_SEDMODEL_DEFAULT);

  return SUCCESS;

} // end of reset_SEDMODEL

// ******************************
int init_primary_SEDMODEL(char *refname  // (I) name of primary ref
		       ,int NLAM         // (I) number of lambda bins
		       ,double *LAMLIST  // (I) lambda-list for spectrum
		       ,double *FLUXLIST // (I) flux vs. lambda (erg/s/cm^2/A)
		       ) {

  int ilam;
  char fnam[] = "init_primary_SEDMODEL" ;

  // -------------- BEGIN ---------------

  if ( NLAM >= MXBIN_PRIMARY_SEDMODEL ) {
    sprintf(c1err,"NLAM(%s) = %d  exceeds array bound of %d", 
	    refname, NLAM, MXBIN_PRIMARY_SEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  if ( NLAM <= 10 ) {
    sprintf(c1err,"NLAM = %d is too small for primary = '%s' ", 
	    NLAM, refname );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  PRIMARY_SEDMODEL.NLAM = NLAM;

  for ( ilam=0; ilam < NLAM; ilam++ ) {
    PRIMARY_SEDMODEL.lam[ilam]  = LAMLIST[ilam];
    PRIMARY_SEDMODEL.flux[ilam] = FLUXLIST[ilam];
  }

  PRIMARY_SEDMODEL.lammin  = PRIMARY_SEDMODEL.lam[0] ;
  PRIMARY_SEDMODEL.lammax  = PRIMARY_SEDMODEL.lam[NLAM-1] ;
  PRIMARY_SEDMODEL.lamstep = PRIMARY_SEDMODEL.lam[1] - PRIMARY_SEDMODEL.lam[0];

  sprintf(PRIMARY_SEDMODEL.name, "%s", refname);

  return(SUCCESS) ;

} // end of init_primary_SEDMODEL


// ********************************
double interp_primaryFlux_SEDMODEL(double lam){

  double 
    lamstep, lammin
    ,flux
    ,a_lam[3]
    ,a_flux[3]
    ;

  int ilam, NLAM ;
  int NBIN_INTERP = 3;

  char fnam[] = "interp_primaryFlux_SEDMODEL";

  // -------- BEGIN --------

  flux = 0.0;
  lamstep = PRIMARY_SEDMODEL.lamstep;
  lammin  = PRIMARY_SEDMODEL.lammin;
  NLAM    = PRIMARY_SEDMODEL.NLAM;

  ilam = (int)((lam - lammin)/lamstep + 1.0E-8) ;
  if ( ilam >= NLAM-1 ) { ilam = NLAM-2; } // edge protect, Aug 23 2019

  a_lam[0] = PRIMARY_SEDMODEL.lam[ilam-1] ;
  a_lam[1] = PRIMARY_SEDMODEL.lam[ilam] ;
  a_lam[2] = PRIMARY_SEDMODEL.lam[ilam+1] ;

  a_flux[0] = PRIMARY_SEDMODEL.flux[ilam-1] ;
  a_flux[1] = PRIMARY_SEDMODEL.flux[ilam] ;
  a_flux[2] = PRIMARY_SEDMODEL.flux[ilam+1] ;

  // idiot check
  if ( lam < a_lam[0] || lam > a_lam[2] ) {
    print_preAbort_banner(fnam);
    printf("\t Expected lam-range: %6.1f - %6.1f A  ilam=%d/%d . \n",
	   a_lam[0], a_lam[2], ilam, PRIMARY_SEDMODEL.NLAM );
    sprintf(c1err,"lam=%6.1f is outside expected range", lam);
    sprintf(c2err,"Check LAMBDA_RANGE key in kcor-input file."); 	    
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // interpolate
  flux = interp_1DFUN (OPT_INTERP_LINEAR, lam, 
		       NBIN_INTERP, a_lam, a_flux, fnam ) ;

  return flux;

} // end of interp_primaryFlux_SEDMODEL



// =================================================
double interp_primaryMag_SEDMODEL(double lam) {

  // Aug 1 2016
  // Interpolate/extrapolate primary mag among broadband
  // filters to estimate mag at wavelength "lam".

  int    NFILT = NFILT_SEDMODEL ;
  double mag;
  char fnam[] = "interp_primaryMag_SEDMODEL" ;

  // ------------ BEGIN -------------

  if ( lam < PRIMARY_SEDMODEL.LAM_SORTED[1] ) 
    { mag = PRIMARY_SEDMODEL.MAG_SORTED[1] ; }
  else if ( lam > PRIMARY_SEDMODEL.LAM_SORTED[NFILT] ) 
    { mag = PRIMARY_SEDMODEL.MAG_SORTED[NFILT] ; }
  else {
    mag = interp_1DFUN (OPT_INTERP_LINEAR, lam, 
			NFILT_SEDMODEL, 
			&PRIMARY_SEDMODEL.LAM_SORTED[1],
			&PRIMARY_SEDMODEL.MAG_SORTED[1] , fnam);   
  }

  return(mag);

} // end interp_primaryMag_SEDMODEL

// ***********************************************
int init_filter_SEDMODEL(
			 int ifilt_obs        // (I) obs filter index
			 ,char   *filter_name // (I) filter name
			 ,char   *survey_name // (I) name of survey
			 ,double  magprimary  // (I) primary mag
			 ,int     NLAM        // (I) Number of lambda bins
			 ,double *LAMLIST     // (I) array of lambda
			 ,double *TRANSSNLIST // (I) array of SN filt-trans
			 ,double *TRANSREFLIST // (I) idem for ref
			 ,double  LAMSHIFT     // (I) shift filter curve
			 )  {

  // Utility to pass filter-response information to SEDMODEL.
  // Call this function (from external program) 
  // for each filter before calling init_genmag_SEDMODEL
  // Don't print anything here since the filter-summary
  // is printed from init_genmag_SEDMODEL().
  //
  // Nov 10 2020: allow modifying already existing filter.

  int ilam, ifilt ;
  char fnam[] = "init_filter_SEDMODEL" ;

  double
    lam, lamstep, transSN, transREF, transSN_MAX, transREF_MAX
    ,fluxREF, fluxREF_sum, fluxSN_sum, transSN_sum, transREF_sum
    ,lamtransSN_sum, lamtransREF_sum, hc8
    ;

  // to debug filter-trans
  int LDMP_FILT = 0 ; 
  FILE *fp_filt;
  char filtFile[80];
  char cfilt1[2];

  // ----- BEGIN -----

  /*  
  printf(" init_filter_SEDMODEL: ifilt_obs=%d  filtname = %s  NLAM=%d \n", 
	 ifilt_obs, filter_name, NLAM );  
  */

  if ( IFILTMAP_SEDMODEL[ifilt_obs] < 0 ) {
    NFILT_SEDMODEL++ ;
    ifilt = NFILT_SEDMODEL; // sparse filter index
  }
  else {
    // Nov 10 2020: update already defined filter
    ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
  }

  if ( ifilt >= MXFILT_SEDMODEL ) {
    filtdump_SEDMODEL();
    sprintf(c1err, "NFILT_SEDMODEL = %d exceeds bound.", NFILT_SEDMODEL);
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }


  if ( NLAM >= MXBIN_LAMFILT_SEDMODEL ) {    
    sprintf(c1err,"NLAM(%s) = %d  exceeds array bound of %d", 
	    filter_name, NLAM, MXBIN_LAMFILT_SEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  lamstep     = LAMLIST[1] - LAMLIST[0] ;
  fluxREF_sum = transREF_sum = lamtransREF_sum = 0.0 ;
  fluxSN_sum  = transSN_sum  = lamtransSN_sum  = 0.0 ;
  transSN_MAX = transREF_MAX = 0.0 ;

  for ( ilam=0; ilam < NLAM; ilam++ ) {
    lam       = LAMLIST[ilam]  + LAMSHIFT ;
    transSN   = TRANSSNLIST[ilam] ;
    transREF  = TRANSREFLIST[ilam] ;
    FILTER_SEDMODEL[ifilt].lam[ilam]   = lam ;
    FILTER_SEDMODEL[ifilt].transSN[ilam] = transSN ; 
    FILTER_SEDMODEL[ifilt].transREF[ilam] = transREF ; 

    if ( transSN  > transSN_MAX  ) { transSN_MAX  = transSN  ; }
    if ( transREF > transREF_MAX ) { transREF_MAX = transREF ; }

    // interpolate primary reference at this lambda.
    fluxREF = interp_primaryFlux_SEDMODEL(lam) ;
    fluxREF_sum  += transREF * fluxREF * lam ;

    // sums used to define mean filter wavelength
    // Could have used SN or REF trans, we choose to pick SN
    lamtransSN_sum += lam * transSN;
    transSN_sum    += transSN ;      

  } // ilam

  if ( transSN_sum < 0.  ) {
    sprintf(c1err,"transSN_sum = %f for ifilt_obs=%d (%s) \n",
	    transSN_sum, ifilt_obs, filter_name );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }


  // load FILTER structure
  IFILTMAP_SEDMODEL[ifilt_obs]      = ifilt ;
  FILTER_SEDMODEL[ifilt].ifilt_obs  = ifilt_obs ;
  FILTER_SEDMODEL[ifilt].magprimary = magprimary ;
  FILTER_SEDMODEL[ifilt].lamshift   = LAMSHIFT ;


  if ( transSN_sum > 0.0 ) 
    {   FILTER_SEDMODEL[ifilt].mean  = lamtransSN_sum/transSN_sum; }
  else
    {   FILTER_SEDMODEL[ifilt].mean  = 0.0 ; }

  FILTER_SEDMODEL[ifilt].lamstep   =   lamstep ;
  FILTER_SEDMODEL[ifilt].NLAM      =   NLAM ;
  FILTER_SEDMODEL[ifilt].lammin    =   FILTER_SEDMODEL[ifilt].lam[0];
  FILTER_SEDMODEL[ifilt].lammax    =   FILTER_SEDMODEL[ifilt].lam[NLAM-1];
  sprintf(FILTER_SEDMODEL[ifilt].name,   "%s", filter_name);
  sprintf(FILTER_SEDMODEL[ifilt].survey, "%s", survey_name);

  FILTER_SEDMODEL[ifilt].transSN_MAX   = transSN_MAX ;
  FILTER_SEDMODEL[ifilt].transREF_MAX  = transREF_MAX ;

  // strip off last char of filtername to make filter-string list
  int len = strlen(filter_name);
  sprintf(cfilt1, "%c", filter_name[len-1] ) ;

  strcat(FILTLIST_SEDMODEL,cfilt1);

  hc8 = (double)hc ;
  fluxREF_sum *= (lamstep/hc8) ;  // Jan 2010: divide by hc 

  if( fluxREF_sum != 0.0 ) 
    { FILTER_SEDMODEL[ifilt].ZP = 2.5*log10(fluxREF_sum) + magprimary ; }
  else
    { FILTER_SEDMODEL[ifilt].ZP = 0.0 ; }

  // dump filter response to text file (DEBUG only)
  if ( LDMP_FILT == 1 ) {
    sprintf(filtFile, "genmag_SEDMODEL_%s-trans.dat", 
	    FILTER_SEDMODEL[ifilt].name );
    fp_filt = fopen(filtFile, "wt");
    printf("\t dump filter-trans to %s )\n", filtFile );
    
    for ( ilam=0; ilam < FILTER_SEDMODEL[ifilt].NLAM ; ilam++ ) {
      lam       =  FILTER_SEDMODEL[ifilt].lam[ilam] ;
      transSN   =  FILTER_SEDMODEL[ifilt].transSN[ilam] ;
      transREF  =  FILTER_SEDMODEL[ifilt].transREF[ilam] ;
      fprintf(fp_filt,"%8.2f  %8.4f %8.4f\n", lam, transSN, transREF );
    }
    fclose(fp_filt);
    
  }  // end of LDMP_FILT 


  return 0;

}  // end of init_filter_SEDMODEL


// ***********************************************
void filtdump_SEDMODEL(void) {

  int ifilt, ifilt_obs;
  char *name;
  char dashLine[] = 
    "----------------------------------------------------------------------" ;
  char fnam[] = "filtdump_SEDMODEL";

  // --------- BEGIN ----

  name  = PRIMARY_SEDMODEL.name ;
  printf("   %s: Read primary ref '%s' with %d lambda bins \n",
	 fnam, name, PRIMARY_SEDMODEL.NLAM );

  printf("\n");
  printf(
	 "   Defined                                 Lambda           "
	 "%s     %s \n", 
	 name, name );
  printf(
	 "   Filter                  <LAMBDA>        Range(step)      "
	 "mag    ZP \n" );
  printf( "   %s\n", dashLine);

  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {

    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs ;

    printf("   id=%2.2d (%12s)   (%7.1f +%2.0f)  %5.0f-%5.0f(%3.0f)  "
	   "%6.3f  %6.3f  \n"
	   ,FILTER_SEDMODEL[ifilt].ifilt_obs
	   ,FILTER_SEDMODEL[ifilt].name 
	   ,FILTER_SEDMODEL[ifilt].mean - FILTER_SEDMODEL[ifilt].lamshift
	   ,FILTER_SEDMODEL[ifilt].lamshift
	   ,FILTER_SEDMODEL[ifilt].lammin
	   ,FILTER_SEDMODEL[ifilt].lammax
	   ,FILTER_SEDMODEL[ifilt].lamstep
	   ,FILTER_SEDMODEL[ifilt].magprimary
	   ,FILTER_SEDMODEL[ifilt].ZP
	   );

    fflush(stdout);
  }  // end of ifilt loop

  printf( "   %s\n", dashLine);


  fflush(stdout);

} // end of filtdump_SEDMODEL



// ***********************************
void init_MWXT_SEDMODEL(int OPT_COLORLAW, double RV) {

  // Created Sep 2013
  // ------------------------
  
  MWXT_SEDMODEL.OPT_COLORLAW = OPT_COLORLAW ;
  MWXT_SEDMODEL.RV = RV ;

} // end of init_MWXT_SEDMODEL

// **************************************
void init_redshift_SEDMODEL(int NZbin, double Zmin, double Zmax) {

  // Nov 2010: optional function to set redshift range of SED grid
  // Jul 2017: fill REDSHIFT_SEDMODEL.ZTABLE[iz]
  // Jan 11 2022: 
  //    Fix bug so that last bin is at zmax, not zmax-zbin
  //    For z > zmax-zbin, flux was based on SED(zmax-zbin) + MU
  //

  char fnam[] = "init_redshift_SEDMODEL";

  // ---------- BEGIN ---------

  if ( NZbin >= MXZBIN_SEDMODEL ) {
    sprintf(c1err,"NZbin=%d exceeds array bound of", NZbin);
    sprintf(c2err,"MXZBIN_SEDMODEL = %d.  Zmin/Zmax=%6.4f/%6.4f", 
	    MXZBIN_SEDMODEL, Zmin, Zmax );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  REDSHIFT_SEDMODEL.NZBIN    = NZbin ;
  REDSHIFT_SEDMODEL.ZMIN     = Zmin  ;
  REDSHIFT_SEDMODEL.ZMAX     = Zmax  ;
  REDSHIFT_SEDMODEL.LOGZMIN  = log10(Zmin)  ;
  REDSHIFT_SEDMODEL.LOGZMAX  = log10(Zmax)  ;

  // setup table. 
  double LOGZMIN, LOGZMAX, LOGZBIN, logz, z;
  int iz;
  LOGZMIN = log10(Zmin);
  LOGZMAX = log10(Zmax);
  // xxx mark delete Jan 11 2022 LOGZBIN = (LOGZMAX-LOGZMIN) / (double)NZbin ;
  LOGZBIN = (LOGZMAX-LOGZMIN) / (double)(NZbin-1) ;
  logz, z;
  for( iz = 0; iz <= NZbin ; iz++ ) { 
    if ( iz == 0 )  
      { z = 0.0; logz = -999. ; } 
    else {                                                                     
      logz = LOGZMIN + LOGZBIN * (double)(iz-1) ;
      z    = pow(10.0 , logz); 
    }
    REDSHIFT_SEDMODEL.ZTABLE[iz]    = z ; 
    REDSHIFT_SEDMODEL.LOGZTABLE[iz] = logz ; 
  } // iz                     

} // end of init_redshift_SEDMODEL

// ******************************************
void malloc_FLUXTABLE_SEDMODEL( int NFILT, int NZBIN, int NLAMPOW, 
				int NDAY, int NSED ) {

  // Nov 24, 2008: allocate flux-integral memory for NSED & NZBIN
  // Jan 30, 2010: switch from fancy 5-dim pointer to 1d pointer
  // Dec 15, 2021: fix isize=sizeof(float) instead of pointer size.

  int isize;
  char fnam[] = "malloc_FLUXTABLE_SEDMODEL" ;

  // -------------- BEGIN -----------------

  // check args against limits

  if ( NFILT <= 0 || NFILT > MXFILT_SEDMODEL ) {
    sprintf(c1err,"NFILT=%d is invalid", NFILT );
    sprintf(c2err,"Check filter-init");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( NZBIN <= 0 || NZBIN > MXZBIN_SEDMODEL ) {
    sprintf(c1err,"NZBIN=%d  is invalid.",  NZBIN );
    sprintf(c2err,"Check redshift bins");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( NLAMPOW < 0 || NLAMPOW > MXLAMPOW_SEDMODEL ) {
    sprintf(c1err,"NLAMPOW=%d  is invalid.",  NLAMPOW );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  }

  if ( NDAY <=0 || NDAY > MXBIN_DAYSED_SEDMODEL ) {
    sprintf(c1err,"NDAY=%d  is invalid (MXBIN_DAY=%d).", 
	    NDAY, MXBIN_DAYSED_SEDMODEL );
    sprintf(c2err,"Check epoch bins");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( NSED <=0 || NSED > MXSEDMODEL ) {
    sprintf(c1err,"NSED=%d  is invalid.",  NSED );
    sprintf(c2err,"Check number of SEDs");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  NBIN_SEDMODEL_FLUXTABLE[IDIM_SEDMODEL_FILTER]   = NFILT ;
  NBIN_SEDMODEL_FLUXTABLE[IDIM_SEDMODEL_REDSHIFT] = NZBIN ;
  NBIN_SEDMODEL_FLUXTABLE[IDIM_SEDMODEL_LAMPOW]   = NLAMPOW ;
  NBIN_SEDMODEL_FLUXTABLE[IDIM_SEDMODEL_DAY]      = NDAY ;
  NBIN_SEDMODEL_FLUXTABLE[IDIM_SEDMODEL_SED]      = NSED ;

  N1DBINOFF_SEDMODEL_FLUXTABLE[0] = 
    (NSED+1) * (NDAY+1) * (NLAMPOW+1) * (NZBIN+1) * (NFILT+1) ;
  N1DBINOFF_SEDMODEL_FLUXTABLE[1] = 
    (NSED+1) * (NDAY+1) * (NLAMPOW+1) * (NZBIN+1) ;
  N1DBINOFF_SEDMODEL_FLUXTABLE[2] = 
    (NSED+1) * (NDAY+1) * (NLAMPOW+1) ;
  N1DBINOFF_SEDMODEL_FLUXTABLE[3] = 
    (NSED+1) * (NDAY+1) ;
  N1DBINOFF_SEDMODEL_FLUXTABLE[4] = 
    (NSED+1) ;
  N1DBINOFF_SEDMODEL_FLUXTABLE[5] = 
    1;


  sprintf(VARNAME_SEDMODEL_FLUXTABLE[1],"NFILT");
  sprintf(VARNAME_SEDMODEL_FLUXTABLE[2],"NZBIN");
  sprintf(VARNAME_SEDMODEL_FLUXTABLE[3],"NLAMPOW");
  sprintf(VARNAME_SEDMODEL_FLUXTABLE[4],"NDAY");
  sprintf(VARNAME_SEDMODEL_FLUXTABLE[5],"NSED");

  //  isize = sizeof(PTR_SEDMODEL_FLUXTABLE);
  isize = sizeof(float);
  NBTOT_SEDMODEL_FLUXTABLE = N1DBINOFF_SEDMODEL_FLUXTABLE[0] ;
  ISIZE_SEDMODEL_FLUXTABLE = NBTOT_SEDMODEL_FLUXTABLE * isize ;

  PTR_SEDMODEL_FLUXTABLE =  (float*)malloc(ISIZE_SEDMODEL_FLUXTABLE);

  printf("  %s : allocate %6.2f Mb of memory for integral-flux tables. \n", 
	 fnam, 1.E-6*(double)ISIZE_SEDMODEL_FLUXTABLE );

  //  printf("\t\t Tables include lambda powers up to %d .\n",  NLAMPOW );
  printf("\t Table bins include %3d DAYs. \n",    NDAY);
  printf("\t Table bins include %3d FILTERs. \n", NFILT);
  printf("\t Table bins include %3d SEDs. \n",    NSED);
  printf("\t Table bins include %3d log10(Z): %.5f <= Z <= %.5f .\n"
	 ,REDSHIFT_SEDMODEL.NZBIN
	 ,REDSHIFT_SEDMODEL.ZMIN
	 ,REDSHIFT_SEDMODEL.ZMAX );

  // print few z-bins at low and high end
  int iz;  double z, logz ;
  for(iz=1; iz <= NZBIN; iz++ ) {
    if ( iz <= 3 || iz >= NZBIN-3 ) {
      z    = REDSHIFT_SEDMODEL.ZTABLE[iz] ;
      logz = REDSHIFT_SEDMODEL.LOGZTABLE[iz] ;
      printf("\t\t ZTABLE[%3d] = %8.5f   (LOGZ=%.4f)\n", 
	     iz, z,  logz ); 
    }
  }


  zero_flux_SEDMODEL();

  // - - - - - - - 
 

  printf("\n");
  fflush(stdout);

  return ;

} // end of malloc_FLUXTABLE_SEDMODEL


void malloc_SEDFLUX_SEDMODEL(SEDMODEL_FLUX_DEF *SEDMODEL_FLUX,
			     int NBIN_DAY_USER, int NBIN_LAM_USER, 
			     int NBIN_SED_USER ) {

  // Created Mar 7 2017
  // May 2018: 
  //  + if any input NBIN > 0, use it instead of default
  //  + pass structure to modify instead of operating on TEMP_SEDMODEL
  //  + rename malloc_TEMP_SEDMODEL -> malloc_SEDFLUX_SEDMODEL

  int MEMD      = sizeof(double);
  int NBIN_DAY  = MXBIN_DAYSED_SEDMODEL;
  int NBIN_LAM  = MXBIN_LAMSED_SEDMODEL;
  int NBIN_SED  = MXBIN_SED_SEDMODEL ;

  // ------------- BEGIN ----------------
  if ( NBIN_DAY_USER  > 0 ) { NBIN_DAY  = NBIN_DAY_USER; }
  if ( NBIN_LAM_USER  > 0 ) { NBIN_LAM  = NBIN_LAM_USER; }
  if ( NBIN_SED_USER  > 0 ) { NBIN_SED  = NBIN_SED_USER; }


  SEDMODEL_FLUX->N_FINEBIN = 0 ; // Nov 17 2016
  SEDMODEL_FLUX->DAY     = (double*) malloc ( MEMD * NBIN_DAY );
  SEDMODEL_FLUX->LAM     = (double*) malloc ( MEMD * NBIN_LAM );
  SEDMODEL_FLUX->FLUX    = (double*) malloc ( MEMD * NBIN_SED );
  SEDMODEL_FLUX->FLUXERR = (double*) malloc ( MEMD * NBIN_SED );

  int i;
  for(i=0; i < MXBIN_SED_SEDMODEL; i++ ) {
    SEDMODEL_FLUX->FLUX[i] = 0.0 ;
    SEDMODEL_FLUX->FLUXERR[i] = 0.0 ;
  }
  return ;

} // end malloc_SEDFLUX_SEDMODEL

// ***************************************
void zero_flux_SEDMODEL(void) {
  int  index ;
  long int N;

  // zero the entire table
  N = N1DBINOFF_SEDMODEL_FLUXTABLE[0] ;
  printf("\t\t Zero entire SEDMODEL flux table (%ld entries). \n", N );
  for ( index = 0 ; index < N; index++ )
    { PTR_SEDMODEL_FLUXTABLE[index] = 0.0 ; }


} // endof zero_flux_SEDMODEL



// ********************************
int NSED_SEDMODEL(void) {
  return SEDMODEL.NSURFACE ;
}

int NPAR_SEDMODEL(void) {
  return SEDMODEL.NPAR;
}

// *********************************
int IPAR_SEDMODEL(char *parName) {

  // for inpuar *parName, return ipar index.
  // If *parName is not defined, return -9.

  int IPAR, ipar  ;

  IPAR = -9 ;

  for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
    if ( strcmp(parName,SEDMODEL.PARNAMES[ipar]) == 0 ) { IPAR=ipar; }
  }

  return IPAR ;

  
} // end of IPAR_SEDMODEL 


// ********************************************************
void fetch_parVal_SEDMODEL(int ISED, int IPAR, double *PARVAL) {
  // Jan 2 2019
  // For input ISED and IPAR index, return *PARVAL 
  // ------------- BEGIN -------------
  *PARVAL = SEDMODEL.PARVAL[ISED][IPAR];
  return ; 
} // end fetch_parVal_SEDMODEL


int fetch_parinfo_sedmodel__(int *ipar,char *parname,int *NBIN,double *range)
{ return( fetch_parInfo_SEDMODEL(*ipar, parname, NBIN, range) ); }

void fetch_parval_sedmodel__(int *ISED, int *IPAR, double *PARVAL)
{ fetch_parVal_SEDMODEL(*ISED,*IPAR,PARVAL); }


// ****************************************************
int fetch_parInfo_SEDMODEL(int ipar, char *parname,int *NBIN,double *PARLIM){

  // Dec 20, 2009
  // for input ipar, return *parname, number of bins,
  // and min/max value
  // Function returns SUCCESS if ipar is defined; ERROR otherwise

  double parval, parval_min, parval_max, dif ;
  int ised, i2, NBIN_LOCAL, USED, I8;
  double *TMPLIST ;

  char fnam[] = "fetch_parInfo_SEDMODEL" ;

  // ----------------- BEGIN ----------------

  // make sure that ipar is  valid
  if ( ipar < 0 || ipar >= SEDMODEL.NPAR ) {
    sprintf(c1err,"invalid arg: ipar=%d  (SEDMODEL.NPAR=%d)", 
	    ipar, SEDMODEL.NPAR );
    sprintf(c2err,"Check SIMSED model parameters." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
   return ERROR ;
  }

  sprintf(parname, "%s", SEDMODEL.PARNAMES[ipar] );

  // get the min/max value for this ipar

  NBIN_LOCAL = 0 ;
  parval_min = +1.0e+19 ;
  parval_max = -1.0e+19 ;
  I8         = sizeof(double);
  TMPLIST    = (double*)malloc(I8*SEDMODEL.NSURFACE+I8);

  for ( ised = 1; ised <= SEDMODEL.NSURFACE ; ised++ ) {
    parval = SEDMODEL.PARVAL[ised][ipar];
    if ( parval > parval_max ) { parval_max = parval ; }
    if ( parval < parval_min ) { parval_min = parval ; }

    // keep track on how many unique parameter values (NBINS)    
    USED = 0;    
    for ( i2=1; i2 <= NBIN_LOCAL; i2++ ) {
      dif = parval - TMPLIST[i2] ;
      if ( fabs(dif) < 1.0E-9 ) { USED = 1 ; }
    } 
    if ( USED == 0 ) {
      NBIN_LOCAL++ ;
      TMPLIST[NBIN_LOCAL] = parval ;
    }


  } // ised

  free(TMPLIST);

  PARLIM[0] = parval_min ;
  PARLIM[1] = parval_max ;
  *NBIN    = NBIN_LOCAL ;

  return SUCCESS ;

} // end of fetch_parInfo_SEDMODEL



// *****************************************************
void init_flux_SEDMODEL(int ifilt_obs, int ised) {

  // Apr 11, 2009 RSK - 
  // for input filter and surface index, 
  // compute integrals for all redshifts & epochs.
  //
  // May 5, 2011
  // Major change: integrate over filter-LAMOBS instead
  // of over LAMSED in the rest-frame. The LAMSED integration
  // does not fully sample the filter-bins at high-redshift.
  // Note that 'ilampow' is now an observer-frame lambda-power.
  //
  // Since the Trest 'quadInterp' is the slowest part,
  // create a very fine-grained SED in lambda space ...
  // then use fast linear interp inside the integration loop.
  //
  // Jun 9, 2011: compute SEDMODEL.LAMMIN_ALL and SEDMODEL.MAXLAM_ALL
  //
  // Mar 22 2017: F->0 if any part of filter trans is not contained by model.
  //
  // Apr 30 2018: store MINDAY_ALL and MAXDAY_ALL
  // Nov 15 2020: protect ifilt for ifilt_obs==0

  int  ilampow, iep, ifilt, ilamfilt, iz ;
  int  NLAMFILT, NLAMSED, EPMIN, EPMAX, N, NZBIN, index ;

  double 
    lamsed, lamobs, lampow, trans, day, z, z1, logzdif
    ,LOGZMIN, LOGZMAX, LOGZBIN, FLUX, tmp, mag, mutmp, x0tmp
    ,LAMOBS_MIN, LAMOBS_MAX, LAMOBS_STEP, SEDMODELNORM
    ,LAMTPOW[MXLAMPOW_SEDMODEL+1]
    ,hc8, TDUM=-99.9
    ;


  char *cfilt ;
  char fnam[] = "init_flux_SEDMODEL" ;

  // -------- BEGIN --------

  //  check flag to print SED-summary
  if ( ifilt_obs == 0 && ised == 0 ) {

    printf("\n");

    printf("  SED-epoch  range: %6.2f <= Trest <= %6.2f (<STEP>=%4.2f day)\n",
	TEMP_SEDMODEL.DAYMIN, TEMP_SEDMODEL.DAYMAX, TEMP_SEDMODEL.DAYSTEP );
    
    printf("  SED-lambda range: %5.0f <= LAM <= %5.0f (<STEP>=%4.0f A) \n",
	TEMP_SEDMODEL.LAMMIN, TEMP_SEDMODEL.LAMMAX, TEMP_SEDMODEL.LAMSTEP );

    printf("  FILTER-lambda range: %5.0f <= LAM <= %5.0f  \n",
	   SEDMODEL.RESTLAMMIN_FILTERCEN, SEDMODEL.RESTLAMMAX_FILTERCEN );

    if ( TEMP_SEDMODEL.DAYMIN < SEDMODEL.DAYMIN_ALL ) 
      { SEDMODEL.DAYMIN_ALL = TEMP_SEDMODEL.DAYMIN ; }
    if ( TEMP_SEDMODEL.DAYMAX > SEDMODEL.DAYMAX_ALL ) 
      { SEDMODEL.DAYMAX_ALL = TEMP_SEDMODEL.DAYMAX ; }

    return ;
  }


  if ( ifilt_obs > 0 )
    { ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ; }
  else
    { ifilt = 0; }

  cfilt     = FILTER_SEDMODEL[ifilt].name ;
  
  if ( SEDMODEL.NSURFACE <= 2 ) {
    printf("  Initialize SED flux-integrals for  %s  and  ISED=%d \n",
	   cfilt, ised);
  }

  // compute misc. stuff
  N = TEMP_SEDMODEL.NLAM ;
  TEMP_SEDMODEL.LAMMIN  = TEMP_SEDMODEL.LAM[0] ;
  TEMP_SEDMODEL.LAMMAX  = TEMP_SEDMODEL.LAM[N-1] ;

  N = TEMP_SEDMODEL.NDAY ;
  TEMP_SEDMODEL.DAYMIN  = TEMP_SEDMODEL.DAY[0] ;
  TEMP_SEDMODEL.DAYMAX  = TEMP_SEDMODEL.DAY[N-1] ;
  
  if ( TEMP_SEDMODEL.DAYMIN < SEDMODEL.DAYMIN_ALL ) 
    { SEDMODEL.DAYMIN_ALL = TEMP_SEDMODEL.DAYMIN ; }
  if ( TEMP_SEDMODEL.DAYMAX > SEDMODEL.DAYMAX_ALL ) 
    { SEDMODEL.DAYMAX_ALL = TEMP_SEDMODEL.DAYMAX ; }


  // store binning info for this "ised" (TEMP_SEDMODEL -> SEDMODEL)
  // This is needlessly repeated for each filter, but it's OK.
  SEDMODEL.NDAY[ised]    = TEMP_SEDMODEL.NDAY ;
  SEDMODEL.DAYMIN[ised]  = TEMP_SEDMODEL.DAYMIN  ;
  SEDMODEL.DAYMAX[ised]  = TEMP_SEDMODEL.DAYMAX ;
  SEDMODEL.DAYSTEP[ised] = TEMP_SEDMODEL.DAYSTEP ;

  SEDMODEL.NLAM[ised]    = TEMP_SEDMODEL.NLAM ;
  SEDMODEL.LAMMIN[ised]  = TEMP_SEDMODEL.LAMMIN ;
  SEDMODEL.LAMMAX[ised]  = TEMP_SEDMODEL.LAMMAX ;
  SEDMODEL.LAMSTEP[ised] = TEMP_SEDMODEL.LAMSTEP ;

  // Aug 12 2017 : 
  // Check flag to use explicit DAY array that could
  // have non-uniform bins. Array is malloced elsewhere
  if ( (SEDMODEL.OPTMASK & OPTMASK_DAYLIST_SEDMODEL)>0 ) {
    for(iep=0; iep < TEMP_SEDMODEL.NDAY; iep++ ) 
      { SEDMODEL.DAY[ised][iep] = TEMP_SEDMODEL.DAY[iep] ;  } 
  }


  // keep track of min and max LAM valid for ALL SEDs ...
  // needed for the IFILTSTAT_SEDMODEL function.
  if ( ised == 1 ) {
    SEDMODEL.LAMMIN_ALL =  0.0  ;      // init
    SEDMODEL.LAMMAX_ALL =  9999999. ;  // init
  }
  if ( TEMP_SEDMODEL.LAMMIN > SEDMODEL.LAMMIN_ALL ) 
    {  SEDMODEL.LAMMIN_ALL = TEMP_SEDMODEL.LAMMIN ; }

  if ( TEMP_SEDMODEL.LAMMAX < SEDMODEL.LAMMAX_ALL ) 
    {  SEDMODEL.LAMMAX_ALL = TEMP_SEDMODEL.LAMMAX ; }

  
  // local variables.
  NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
  NLAMSED   = TEMP_SEDMODEL.N_FINEBIN ;

  // error checks.

  if ( SEDMODEL.FLUXSCALE <= 0.0 ) {
    sprintf(c1err,"SEDMODEL.FLUXSCALE = %f has not been set.", 
	    SEDMODEL.FLUXSCALE );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  if ( ifilt < 0 || ifilt >= MXFILT_SEDMODEL ) {
    sprintf(c1err,"invalid ifilt=%d for ifilt_obs=%d", ifilt, ifilt_obs);
    sprintf(c2err,"array bound is MXFILT_SEDMODEL=%d", MXFILT_SEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( ised < 0 || ised >= MXSEDMODEL ) {
    sprintf(c1err,"invalid ised=%d . Array bound is %d", ised, MXSEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }

  if ( TEMP_SEDMODEL.NDAY >= MXBIN_DAYSED_SEDMODEL ) {
    sprintf(c1err,"NDAY(ISED=%d) = %d exceeds array bound.", 
	    ised, TEMP_SEDMODEL.NDAY );
    sprintf(c2err,"Check MXBIN_DAYSED_SEDMODEL = %d", MXBIN_DAYSED_SEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  if ( TEMP_SEDMODEL.NLAM >= MXBIN_LAMSED_SEDMODEL ) {
    sprintf(c1err,"NLAM(ISED=%d) = %d exceeds array bound.", 
	    ised, TEMP_SEDMODEL.NLAM );
    sprintf(c2err,"Check MXBIN_LAMSED_SEDMODEL = %d", MXBIN_LAMSED_SEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

    

  // ==========================================

  // fill Z-table
  LOGZMIN = REDSHIFT_SEDMODEL.LOGZMIN;
  LOGZMAX = REDSHIFT_SEDMODEL.LOGZMAX;
  NZBIN   = REDSHIFT_SEDMODEL.NZBIN ;


  logzdif = (LOGZMAX - LOGZMIN ) ;

  if ( NZBIN > 1 ) 
    { LOGZBIN = logzdif/(double)(NZBIN-1) ; }
  else
    { LOGZBIN = logzdif ; }

  REDSHIFT_SEDMODEL.LOGZBIN = LOGZBIN ; // store in global array


  if ( NZBIN >= MXZBIN_SEDMODEL ) {
    sprintf(c1err,"NZBIN_SEDMODEL=%d exceeds array bound of %d", 
	    NZBIN, MXZBIN_SEDMODEL );
    errmsg(SEV_FATAL, 0, fnam, c1err, ""); 
  }


  // check option to skip flux-table calc (if binary file is read)
  if ( ifilt_obs == 0 ) { return ; }

  EPMIN = 0 ; 
  EPMAX = TEMP_SEDMODEL.NDAY - 1 ;

  LAMOBS_STEP = FILTER_SEDMODEL[ifilt].lamstep ;
  LAMOBS_MIN  = FILTER_SEDMODEL[ifilt].lammin ;
  LAMOBS_MAX  = FILTER_SEDMODEL[ifilt].lammax ;

  hc8 = (double)hc;
  SEDMODELNORM = SEDMODEL.FLUXSCALE * LAMOBS_STEP / hc8 ;

  // =========================================
  // start  main loop with SED lambda in rest-frame
  for ( iz=0; iz <= NZBIN; iz++ ) {

    z         = REDSHIFT_SEDMODEL.ZTABLE[iz];
    z1        = 1. + z ;

    // Mar 22 2017: 
    // bail if any part of filter trans it outside of model range
    if ( LAMOBS_MIN/z1 < SEDMODEL.LAMMIN[ised] ) { continue ; }
    if ( LAMOBS_MAX/z1 > SEDMODEL.LAMMAX[ised] ) { continue ; }

    // normalization for integrals
    for ( ilamfilt=0; ilamfilt < NLAMFILT; ilamfilt++ ) {

      lamobs = FILTER_SEDMODEL[ifilt].lam[ilamfilt];      
      trans  = FILTER_SEDMODEL[ifilt].transSN[ilamfilt]; 
      lamsed = lamobs/z1 ;

      // bail of lamsed is outside SED range
      if ( lamsed < SEDMODEL.LAMMIN[ised] ) { continue ; }
      if ( lamsed > SEDMODEL.LAMMAX[ised] ) { continue ; }
      
      if ( trans <= 0.0 ) { continue ; }
  
      lampow = 1.0 ; 
      // store powers of obs-lambda before looping over other indices
      for ( ilampow=0; ilampow <= NLAMPOW_SEDMODEL; ilampow++ ) {
	LAMTPOW[ilampow] = (SEDMODELNORM * lamsed * trans * lampow) ;
	lampow *= lamobs ;
      } 

      // loop over epochs, and compte flux
      for ( iep=EPMIN ; iep <= EPMAX ; iep++ ) {

	// get flux at this epoch and obs-wavelength
	FLUX = getFluxLam_SEDMODEL(ised, iep, TDUM, lamobs, z, fnam) ; 

	for ( ilampow=0; ilampow <= NLAMPOW_SEDMODEL; ilampow++ ) {
	  tmp   = LAMTPOW[ilampow] * FLUX ;
	  index = INDEX_SEDMODEL_FLUXTABLE(ifilt,iz,ilampow,iep,ised);
	  PTR_SEDMODEL_FLUXTABLE[index] += tmp ;
	} 
	  
      }  // end iep loop over epochs
    } // end ilam loop over SED-lambda bins in rest-frame
  } // end of iz loop


  if ( ised == -9 ) {
    iep   = 20; iz=0 ;
    day   = TEMP_SEDMODEL.DAY[iep]; 
    z     = REDSHIFT_SEDMODEL.ZTABLE[iz];
    mutmp = 0.0 * 23.158781 ;
    x0tmp = pow(10.0,12-0.4*mutmp);
    index = INDEX_SEDMODEL_FLUXTABLE(ifilt,iz,0,iep,ised);
    FLUX  = x0tmp * PTR_SEDMODEL_FLUXTABLE[index] ;
    
    mag   = FILTER_SEDMODEL[ifilt].ZP - 2.5*log10(FLUX) ;
    printf(" xxx peak SNmag(%s) = %f  (day=%f,z=%le) \n", 
	   cfilt, mag, day, z ) ;
  }

  return ;

} // end of init_flux_SEDMODEL


// *********************************************
double getFluxLam_SEDMODEL(int ISED, int IEP, double TOBS, double LAMOBS, 
			   double z, char *funCall ) {

  // Nov 2016
  // Return rest-frame SED Flam for inputs
  //  ised = SED index
  //  iep  = epoch index (if >= 0 )
  //  Tobs = T - Tpeak (if iep<0)
  //  lamobs  = observer wavelength (A)
  //  funCall = call function, and used only for abort message.
  // 
  // Note that the epoch can be input as an index (iep)
  // or a double Tobs (if iep<0). If Tobs is input, then nearest
  // iep index is used (no Tobs interpolation).
  //
  // This function is part of refactor needed so that NON1ASED
  // model can return a spectrum for the SPECTROGRAPH option.
  //
  // BEWARE: 
  //   This works for NON1ASED using TEMP_SEDMODEL struct.
  //   but won't work for SIMSED without switching to SEDMODEL struct.
  //
  // Jan 19 2017: fix to work if just one DAY (i.e., one spectrum)
  //

  double FLUX   = 0.0 ;
  double FLUXTMP[2], fluxTmp[2];
  double FRAC_INTERP_LAM, FRAC_INTERP_DAY ;
  double LAMDIF, LAMDIF2, LAMSED, TREST, z1 ;
  int    ilamsed, i, iep, jflux, ONEDAY ;
  int    NLAMSED = TEMP_SEDMODEL.N_FINEBIN ;
  char   fnam[] = "getFluxLam_SEDMODEL" ;

  double DAYMIN  = TEMP_SEDMODEL.DAYMIN ;
  double DAYMAX  = TEMP_SEDMODEL.DAYMAX ;
  double DAYSTEP = TEMP_SEDMODEL.DAYSTEP ; 
  int    NDAY    = TEMP_SEDMODEL.NDAY ; 

  //  int LDMP = ( fabs(TOBS)<.1 ) ;
  // ----------- BEGIN ------------

  z1 = 1.0 + z ;
  LAMSED = LAMOBS / z1 ;
  
  if ( LAMSED < SEDMODEL.LAMMIN[ISED] ) { return(FLUX); }
  if ( LAMSED > SEDMODEL.LAMMAX[ISED] ) { return(FLUX); }

  ONEDAY = ( NDAY == 1 ) ;

  if ( IEP >= 0 ) 
    { iep = IEP ; FRAC_INTERP_DAY=0.0; }
  else {
    // compute iep from TOBS
    TREST = TOBS/z1 ;

    if ( TREST < DAYMIN ) { return(FLUX); }
    if ( TREST > DAYMAX ) { return(FLUX); }

    get_DAYBIN_SEDMODEL(ISED, TREST, &iep, &FRAC_INTERP_DAY);

    // TREST in SED file can be truncated, so allow a little
    // slop for FRAC_INTERP_DAY
    if ( FRAC_INTERP_DAY < -0.001 || FRAC_INTERP_DAY > 1.001 ) {  
      print_preAbort_banner(fnam);
      printf("\t DAYMIN=%f   DAYMAX=%f  DAYSTEP=%f \n",
             DAYMIN, DAYMAX, DAYSTEP);
      printf("\t DAY[iep=%d] = %f   Trest=%f\n",
             iep, TEMP_SEDMODEL.DAY[iep], TREST );

      sprintf(c1err,"FRAC_DAY=%f for Trest=%.2f, z=%3f",
	      FRAC_INTERP_DAY, TREST, z );
      sprintf(c2err,"Called by %s (ISED=%d, iep=%d)", funCall, ISED, iep);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
    }    
  }

  // get ilamsed that bounds observer-lamobs
  LAMDIF    = LAMSED - SEDMODEL.LAMMIN[ISED] ;
  ilamsed   = (int)(LAMDIF/TEMP_SEDMODEL.FINEBIN_LAMSTEP) ;
  if ( ilamsed > NLAMSED - 1 ) { ilamsed = NLAMSED - 1 ; }
  LAMDIF2 = LAMSED - TEMP_SEDMODEL.FINEBIN_LAM[ilamsed] ;
  FRAC_INTERP_LAM  = LAMDIF2 / TEMP_SEDMODEL.FINEBIN_LAMSTEP ;

  if ( FRAC_INTERP_LAM < -1.0E-8 || FRAC_INTERP_LAM > 1.000000001 ) {
    print_preAbort_banner(fnam);
    printf("\t LAMDIF  = %.2f - %.2f = %.2f \n",
	   LAMSED, SEDMODEL.LAMMIN[ISED], LAMDIF );
    printf("\t ilamsed = (int) %.2f / %.3f  = %d \n", 
	   LAMDIF, TEMP_SEDMODEL.FINEBIN_LAMSTEP, ilamsed);
    printf("\t LAMDIF2 = %.2f - %.2f = %.2f \n",
	   LAMSED, TEMP_SEDMODEL.FINEBIN_LAM[ilamsed], LAMDIF2);
    printf("\t FRAC_INTERP = %.3f / %f = %f \n",
	   LAMDIF2, TEMP_SEDMODEL.FINEBIN_LAMSTEP, FRAC_INTERP_LAM);

    sprintf(c1err,"FRAC_LAM=%f for LAMSED=%6.2f(%d), z=%3f",
	    FRAC_INTERP_LAM, LAMSED, ilamsed, z );
    sprintf(c2err,"Called by %s (ISED=%d, iep=%d)", funCall, ISED, iep);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  }

  // do linear interp for fine-binned flux
  for ( i=0; i < 2; i++ ) {
    jflux      = NLAMSED*(iep+0) + (ilamsed+i) ;
    fluxTmp[0] = TEMP_SEDMODEL.FINEBIN_FLUX[jflux];
    if ( IEP >= 0 || ONEDAY ) 
      { FLUXTMP[i] = fluxTmp[0] ; } // no day-interpolation
    else {
      // use day-interpolation
      jflux      = NLAMSED*(iep+1) + (ilamsed+i) ;
      fluxTmp[1] = TEMP_SEDMODEL.FINEBIN_FLUX[jflux];
      FLUXTMP[i] = fluxTmp[0] + FRAC_INTERP_DAY*(fluxTmp[1]-fluxTmp[0]);
    }

  }

  // interpolate in lambda space
  FLUX  = FLUXTMP[0] + FRAC_INTERP_LAM*(FLUXTMP[1] - FLUXTMP[0]) ; 

  return(FLUX);

} // end getFluxLam_SEDMODEL


// *********************************************
void init_FINEBIN_SEDMODEL(int ised) {

  // Created May 5, 2011
  // fill TEMP_SEDMODEL.FINBIN_XXX variables.
  // This defines an SED with 2 A lambda bins determined
  // here with 'quadInterp'. Allows using fast linear interp
  // in the integration loops in init_flux_SEDMODEL(), 
  // saving significant CPU time for the init stage.
  //
  // ised = 0,1 => process this sed
  // ised =  -1 => free memory
  //
  // Mar 11 2016: LAMSTEP_FIND -> LAMSTE/5 instead of 2.0 to ensure
  //              integer number of divisions
  //
  // Nov 17 2016: require N_FINEBIN>0 to free memory
  //
  int 
     ilam_fine, ilam_orig, jflux_fine,jflux_orig
    ,N_FINEBIN,  N2D, N_REBIN, EPMIN, EPMAX, iep, i, I8
    ,MXLAM_ORIG
    ;

  double 
    LAMSTEP_RATIO
    , LAM_FINE
    , LAM_DIF
    , LAM_MIN
    , FRAC, xlam
    , FLUXTMP[4], LAMTMP[4], FLUX_FINE
    , LAMSTEP_FINE
    ;

  char fnam[] = "init_FINEBIN_SEDMODEL" ;

  int LDMP = 0 ;

  // ----------- BEGIN --------------

  if ( ised < 0 &&  TEMP_SEDMODEL.N_FINEBIN > 0) {
    free(TEMP_SEDMODEL.FINEBIN_LAM);
    free(TEMP_SEDMODEL.FINEBIN_FLUX);
    if ( LDMP ) { printf(" xxx %s : free memory \n", fnam ); }
  }

  if (ised < 0 ) { return ; }

  LAMSTEP_FINE = TEMP_SEDMODEL.LAMSTEP/5.0 ; // Mar 2016

  if ( LAMSTEP_FINE < TEMP_SEDMODEL.LAMSTEP )   
    { TEMP_SEDMODEL.FINEBIN_LAMSTEP = LAMSTEP_FINE ; }
  else
    { TEMP_SEDMODEL.FINEBIN_LAMSTEP = TEMP_SEDMODEL.LAMSTEP ; }


  LAMSTEP_RATIO = TEMP_SEDMODEL.LAMSTEP /TEMP_SEDMODEL.FINEBIN_LAMSTEP ;

  N_REBIN   = (int)LAMSTEP_RATIO ;  // NLAM -> *= N_REBIN
  N_FINEBIN = N_REBIN * TEMP_SEDMODEL.NLAM ;

  if ( LAMSTEP_RATIO != (double)N_REBIN ) {
    sprintf(c1err,"SED LAMSTEP=%3.2f is not a multiple of "
	    "FINEBIN_LAMSTEP=%2.1f",
	    TEMP_SEDMODEL.LAMSTEP, TEMP_SEDMODEL.FINEBIN_LAMSTEP );
    sprintf(c2err,"check LAMBDA binning of SED model.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  }

  /*
  printf("\t Init %d  %2.1f-A lambda bins for SED-%4.4d \n", 
	 N_FINEBIN, TEMP_SEDMODEL.FINEBIN_LAMSTEP, ised);
  */

  TEMP_SEDMODEL.N_FINEBIN = N_FINEBIN ;

  // get size of fine-binned SED
  N2D = N_FINEBIN * TEMP_SEDMODEL.NDAY ;
  I8  = sizeof(double);
  TEMP_SEDMODEL.FINEBIN_LAM  = (double*)malloc(I8*N_FINEBIN);
  TEMP_SEDMODEL.FINEBIN_FLUX = (double*)malloc(I8*N2D);
  MXLAM_ORIG = TEMP_SEDMODEL.NLAM - 3 ;

  LAM_MIN = TEMP_SEDMODEL.LAM[0] ;
  EPMIN = 0 ;
  EPMAX = TEMP_SEDMODEL.NDAY - 1 ;

  // ------- loop over fine-binned lambda -----------

  for ( ilam_fine=0; ilam_fine < N_FINEBIN ; ilam_fine++ ) {
    
    xlam = (double)ilam_fine ;
    LAM_FINE  = LAM_MIN + (xlam * TEMP_SEDMODEL.FINEBIN_LAMSTEP) ;
    TEMP_SEDMODEL.FINEBIN_LAM[ilam_fine] = LAM_FINE ;

    // get ilamsed that bounds observer-lamobs
    LAM_DIF    = LAM_FINE - LAM_MIN ;
    ilam_orig  = (int)(LAM_DIF/TEMP_SEDMODEL.LAMSTEP) ;
    if ( ilam_orig > MXLAM_ORIG ) { ilam_orig = MXLAM_ORIG ; }
    LAM_DIF    = LAM_FINE - TEMP_SEDMODEL.LAM[ilam_orig] ;
    FRAC       = LAM_DIF / TEMP_SEDMODEL.LAMSTEP ;
    if ( FRAC < 0.5 && ilam_orig > 0 ) { ilam_orig-- ; }
    
    
    // loop over epochs, and compte flux
    for ( iep=EPMIN ; iep <= EPMAX ; iep++ ) {
      // do quadratic interp for fine-binned flux
      for ( i=0; i <= 2; i++ ) {
	jflux_orig = TEMP_SEDMODEL.NLAM*iep + (ilam_orig+i) ;
	FLUXTMP[i] = TEMP_SEDMODEL.FLUX[jflux_orig];
	LAMTMP[i]  = TEMP_SEDMODEL.LAM[ilam_orig+i] ;
      }
      FLUX_FINE  = quadInterp( LAM_FINE, LAMTMP, FLUXTMP, fnam );
      
      // store in global TEMP array
      jflux_fine = N_FINEBIN*iep + ilam_fine ;
      TEMP_SEDMODEL.FINEBIN_FLUX[jflux_fine] = FLUX_FINE ;

      if ( iep == -20 && ilam_fine == 122 && ised == 1 ) {
	printf(" XXXX --------------------------------------------- \n");

	printf(" XXXX iep = %d    Trest = %6.1f \n",
	       iep, TEMP_SEDMODEL.DAY[iep] );

	printf(" XXXX FLUX_FINE = %le  at  LAM_FINE(%d)=%6.1f \n",
	       FLUX_FINE, ilam_fine, LAM_FINE );

	printf(" XXXX LAM_ORIG  = %6.1f  %6.1f  %6.1f \n",
	       LAMTMP[0], LAMTMP[1], LAMTMP[2] );

	printf(" XXXX FLUX_ORIG = %le  %le  %le \n",
	       FLUXTMP[0], FLUXTMP[1], FLUXTMP[2] );

	printf("\n");
      }

    } // ep

  } // ilam_fine

  if ( LDMP ) { 
    printf(" xxx %s : allocate memory with N_FINEBIN=%d\n", 
	   fnam, N_FINEBIN ) ; 
  }

} // end of init_FINEBIN_SEDMODEL


// ************************************
double get_flux_SEDMODEL( int ISED, int ilampow, int ifilt_obs,
			  double z, double Trest) {

  // July 20 2018
  // Shell to call interp_flux_SEDMODEL when MINDAY < Trest < MAXDAY.
  // If Trest < MINDAY, return 0.
  // If Trest > MAXDAY the extrapolate magnitude.
  //
  // Aug 25 2018: 
  //  for late-time extrap, compute least-squares fit for last
  //  5 days of model. Should be more robust for noisy models.
  // 

#define NDAYFIT_SLOPE 6

  double MINDAY   = SEDMODEL.DAYMIN[ISED] ;
  double MAXDAY   = SEDMODEL.DAYMAX[ISED] ;
  double MINSLOPE = INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE;

  double XNDAY  = (double)NDAYFIT_SLOPE ;
  double Tref1, Fref1, m1, m, slope;
  double Tlist[NDAYFIT_SLOPE];
  double Flist[NDAYFIT_SLOPE];
  double Mlist[NDAYFIT_SLOPE];
  double T, F, Tsum, Msum, TTsum, TMsum, flux=0.0 ;  
  char fnam[] = "get_flux_SEDMODEL" ;
  int LDMP = 0; 
  int ilist, N_FLUXNEG=0 ; 
  
  // ------------- BEGIN --------------

    
  if ( Trest < MINDAY ) {
    flux = 0.0 ;
  }
  else if ( Trest > MAXDAY-1.0 ) {
    // extraploate late-time

    Tsum = Msum = TTsum = TMsum = 0.0;
    for(ilist=0; ilist < NDAYFIT_SLOPE; ilist++ ) {
      T            = MAXDAY-0.01 - (double)(NDAYFIT_SLOPE-ilist-1) ;
      F            = interp_flux_SEDMODEL(ISED,ilampow,ifilt_obs,z,T);
      Tlist[ilist] = T ;
      Flist[ilist] = F ;
      Mlist[ilist] = 99.0 ;
      if ( F > 1.0E-80 ) {
	m            = 30.0 - 2.5*log10(F); // ZP=30 is arbitrary
	Mlist[ilist] = m ;
	Tsum += T; Msum += m; TTsum += (T*T); TMsum+=(T*m);
      }
      else {
	N_FLUXNEG++ ;
      }
    } // end ilist

    if ( N_FLUXNEG==0 ) {
      slope = ( TMsum - (Tsum*Msum)/XNDAY) / ( TTsum - Tsum*Tsum/XNDAY );
      slope = fabs(slope); // fabs --> force getting dimmer with time
      if ( slope < MINSLOPE ) { slope = MINSLOPE; }
      Fref1 = Flist[NDAYFIT_SLOPE-1];
      Tref1 = Tlist[NDAYFIT_SLOPE-1];
      m1    = Mlist[NDAYFIT_SLOPE-1];
      m     = m1 + (Trest-Tref1)*slope ;
      flux  = pow(TEN,-0.4*(m-30.0));
    }
    else {
      m = slope = flux = 0.0 ;
    }


    if ( LDMP ) {
      printf(" xxx ----------%s DUMP ---------------- \n", fnam);
      printf(" xxx ISED=%d  ifilt_obs=%d  Trest=%f  MAXDAY=%f \n",
	     ISED, ifilt_obs, Trest, MAXDAY);

      for(ilist=0; ilist < NDAYFIT_SLOPE; ilist++ ) {
	printf(" xxx ilist=%d: Trest=%.2f  Fref=%le  Mref=%.3f\n",
	       ilist, Tlist[ilist], Flist[ilist], Mlist[ilist] );
      }
      printf(" xxx slope=%le  m=%.4f  flux=%le  (N_NEG=%d)\n", 
	     slope, m, flux, N_FLUXNEG );
      fflush(stdout);
    }

    
  }
  else {
    // interpolate DAY grid
    flux = interp_flux_SEDMODEL(ISED,ilampow,ifilt_obs,z,Trest);
  }
  
  return (flux);
  
} // end get_flux_SEDMODEL


// ***************************************
double interp_flux_SEDMODEL(
			    int ISED       // SED surface id 1-MXSED
			    ,int ilampow    // lambda power of integral
			    ,int ifilt_obs  // observer  filter index
			    ,double z       // redshift
			    ,double Trest   // Tobs (days; 0=peak)
			    ) {

  // Created Oct 11, 2009 by R.Kessler
  // Return interpolated flux-integral 'S' 
  // from pre-tabulated tables
  //
  // Aug 17 2015: fix long-standing bug; replace TEMP_SEDMODEL.DAY with
  //              DAYLIST from this SED. Previously was using last read
  //              DAYLIST.
  //
  // Jan 19 2017: 
  //   + add NEPBIN_SPLINE and ONEDAY check for single-epoch spectrum
  //
  // July 19 2018: 
  //   + replace quadInterp with linear interp to avoid pathological
  //     parabolic interp that can result in negative flux
  //             
  int 
    ifilt, NDAY, EPMAX, index, LDMP, NZBIN, IZLO, EPLO, ep, iz, NZTMP
    ,NBIN_SPLINE = 3
    ,NZBIN_SPLINE, NEPBIN_SPLINE
    ;

  int  USE_DAYSTEP = 
    ( (SEDMODEL.OPTMASK & OPTMASK_DAYLIST_SEDMODEL)==0 );

  double 
    logz
    ,DAYMIN, DAYSTEP, DAYMAX
    ,DAYLIST_INTERP[4]
    ,LOGZMIN, LOGZMAX, LOGZBIN, FRAC, LOGZ_TMP
    ,S2DTMP[10][10]  // [iz][iday]
    ,SZTMP[10], S, LAMOBS_MIN, LAMOBS_MAX
    ,*ptr_EP, *ptr_LOGZ    
    ,z1 = 1.0 + z
    ;

  char fnam[] = "interp_flux_SEDMODEL" ;
  
  // ---------- BEGIN -------------

  S = 0.0 ;  // init flux to be returned.

  if ( ISED > SEDMODEL.NSURFACE ) { return(S); };

  ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;

  NDAY  = SEDMODEL.NDAY[ISED];
  EPMAX = SEDMODEL.NDAY[ISED] - 1; // because ep starts at zero

  NZBIN   = REDSHIFT_SEDMODEL.NZBIN ;
  LOGZMIN = REDSHIFT_SEDMODEL.LOGZMIN ;
  LOGZMAX = REDSHIFT_SEDMODEL.LOGZMAX ;
  LOGZBIN = REDSHIFT_SEDMODEL.LOGZBIN ;

  NZBIN_SPLINE = NEPBIN_SPLINE = NBIN_SPLINE ; // default

  if ( NZBIN == 1 )  { NZBIN_SPLINE  = 1; }
  if ( NDAY  == 1 )  { NEPBIN_SPLINE = 1; }


  if ( EPMAX < 0 || EPMAX > MXBIN_DAYSED_SEDMODEL-1 ) {
    sprintf(c1err,"invalid EPMAX = %d ", EPMAX );
    sprintf(c2err,"ISED=%d  NDAY=%d z=%.3f Trest=%.3f \n", 
	    ISED, NDAY, z, Trest );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  }

  // find nearest SED-epoch bins that sandwhich Trest
  DAYMIN  = SEDMODEL.DAYMIN[ISED] ;
  DAYMAX  = SEDMODEL.DAYMAX[ISED] ;
  DAYSTEP = SEDMODEL.DAYSTEP[ISED] ;


  if ( Trest < DAYMIN ) { return(S) ; }
  if ( Trest > DAYMAX ) { return(S) ; }


  // if filter trans is not contained in model wavelength range,
  // return undefined flux
  LAMOBS_MIN  = FILTER_SEDMODEL[ifilt].lammin ;
  LAMOBS_MAX  = FILTER_SEDMODEL[ifilt].lammax ;
  if ( LAMOBS_MIN/z1 < SEDMODEL.LAMMIN[ISED] ) { return(FLUX_UNDEFINED) ; }
  if ( LAMOBS_MAX/z1 > SEDMODEL.LAMMAX[ISED] ) { return(FLUX_UNDEFINED) ; }

  if ( ifilt_obs == -3 ) {
    printf(" xxx ifilt_obs=%d  Trest=%6.2f  "
	   "DAY[MIN,MAX,STEP]=%6.2f,%6.2f,%f \n",
	   ifilt_obs, Trest, DAYMIN, DAYMAX, DAYSTEP ); fflush(stdout);
  }


  get_DAYBIN_SEDMODEL(ISED, Trest, &EPLO, &FRAC); // return EPLO & FRAC
  if ( EPLO > 0 && FRAC < 0.5      ) { EPLO-- ; }
  if ( EPLO > NDAY - NEPBIN_SPLINE ) { EPLO = NDAY - NEPBIN_SPLINE ; }

  // construct local DAY-array used for interpation
  if ( USE_DAYSTEP ) {
    for(ep=0; ep < 3; ep++ ) 
      { DAYLIST_INTERP[ep] = DAYMIN + (double)(EPLO+ep) * DAYSTEP; }
  }
  else {
    for(ep=0; ep < 3; ep++ ) 
      { DAYLIST_INTERP[ep] = SEDMODEL.DAY[ISED][EPLO+ep] ; }
  }
  

  // now check redshift.
  if ( z < 1.0E-10 ) 
    {  logz = log10( z + 1.0E-10); }
  else
    { logz = log10(z); }

  if ( z < 0.0 || logz > LOGZMAX+0.00001 ) {
    char cfilt[2];  sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] ) ;
    print_preAbort_banner(fnam);
    printf("\t Filter=%s logz = %f   LOGZMAX=%f \n", cfilt, logz, LOGZMAX );
    printf("\t ISED=%d  Trest=%.2f   \n", ISED,Trest);
    sprintf(c1err,"z=%f is outside SEDMODEL_FLUXTABLE lookup table", z);
    sprintf(c2err,"ZTABLE range is 0 to %f  (logz-LOGZMAX=%le)", 
	    REDSHIFT_SEDMODEL.ZTABLE[NZBIN], logz-LOGZMAX );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // find nearest redshift bin in log10 space
  if ( z < 3.E-9 ) {
    IZLO = 0;
  }
  else if ( logz <= LOGZMIN ) {
    IZLO = 1 ;
  }
  else {
    IZLO     = (int)( (logz - LOGZMIN + 1.0E-9) / LOGZBIN ) + 1 ;
    if( NZBIN > 1 ) {
      LOGZ_TMP = REDSHIFT_SEDMODEL.LOGZTABLE[IZLO];    
      FRAC     = (logz - LOGZ_TMP)/LOGZBIN;    
      if ( IZLO > 1 && FRAC < 0.5       ) { IZLO-- ; }

      // avoid going past upper redshift range
      NZTMP = NZBIN - NZBIN_SPLINE + 1 ;
      if ( IZLO > NZTMP  ) { IZLO = NZTMP ; }
    }
  }



  ptr_EP   = DAYLIST_INTERP ;
  ptr_LOGZ = &REDSHIFT_SEDMODEL.LOGZTABLE[IZLO] ;

  // for now do quad-interp in each dimension ...
  // later should use a real spline-interp
 
  // create local LOGZ x TREST grid 
  for ( iz=0; iz < NZBIN_SPLINE; iz++ ) {
    for ( ep=0; ep < NEPBIN_SPLINE; ep++ ) {
      index = INDEX_SEDMODEL_FLUXTABLE(ifilt,IZLO+iz,ilampow,EPLO+ep,ISED);
      S2DTMP[iz][ep] = PTR_SEDMODEL_FLUXTABLE[index] ;
      
    } // ep

    // interpolate across epoch to get SZ
    if ( NDAY > 1 ) {
	// linear interp to avoid crazy parabolic fits
      SZTMP[iz] = interp_1DFUN (1, Trest,
				  NBIN_SPLINE, ptr_EP, S2DTMP[iz], fnam ) ;
    }
    else
      { SZTMP[iz] = S2DTMP[iz][0] ; }  // no interp of just 1 spectrum

  } // iz

  // Now interpolate across LOGZ

  if ( NZBIN == 1 ) 
    { S = SZTMP[0]; } // 1 z-bin => no interp needed.
  else 
    { S = quadInterp( logz, ptr_LOGZ, SZTMP, fnam ); }

  LDMP = ( EPLO == 20 &&  ifilt_obs == -2  ) ;
  if ( LDMP ) {
    printf(" SSSS ------------------------------------------------------ \n");

    printf(" SSSS z=%5.3f logz=%6.4f (NZBIN=%d,%d)  Trest=%6.2f \n", 
	   z, logz, NZBIN, NZBIN_SPLINE, Trest );
    printf(" SSSS IZLO=%d  ptr_LOGZ = %f, %f, %f \n",
	   IZLO,   *(ptr_LOGZ+0), *(ptr_LOGZ+1), *(ptr_LOGZ+2) );

    printf(" SSSS EPLO=%d  ptr_EP   = %f, %f, %f \n",
	   EPLO, *(ptr_EP+0), *(ptr_EP+1), *(ptr_EP+2) );

    printf(" SSSS S2DTMP[0] = %le, %le, %le \n",
	   S2DTMP[0][0], S2DTMP[0][1], S2DTMP[0][2] );
    printf(" SSSS SZTMP = %le, %le, %le \n",
	   SZTMP[0], SZTMP[1], SZTMP[2] );
    printf(" SSSS Z=%le \n", S );

    //    debugexit("interp_SEDMODEL"); // xxxxxxxx
  } // end of LDMP


  return(S) ;

} // end of interp_flux_SEDMODEL



// *******************************************
long int INDEX_SEDMODEL_FLUXTABLE(int ifilt, int iz, 
			     int ilampow, int iep, int  ised) {

  // Created Jan 30, 2010 by R.Kessler
  // Return 1D index corresponding to the five indices
  // of the FLXUTABLE. 

  long int INDEX;

  // ------- BEGIN --------

  INDEX = 0;

  INDEX += N1DBINOFF_SEDMODEL_FLUXTABLE[1] * ifilt ;
  INDEX += N1DBINOFF_SEDMODEL_FLUXTABLE[2] * iz ;
  INDEX += N1DBINOFF_SEDMODEL_FLUXTABLE[3] * ilampow ;
  INDEX += N1DBINOFF_SEDMODEL_FLUXTABLE[4] * iep ;
  INDEX += N1DBINOFF_SEDMODEL_FLUXTABLE[5] * ised ;

  return INDEX ;

} // end of INDEX_SEDMODEL_FLUXTABLE

// ************************************
double get_magerr_SEDMODEL( int ISED, int ifilt_obs,
			     double z, double Trest) {

  int    NLAM    = TEMP_SEDMODEL.NLAM ;
  int    NDAY    = TEMP_SEDMODEL.NDAY ; 
  double DAYMIN  = TEMP_SEDMODEL.DAYMIN ;
  double DAYMAX  = TEMP_SEDMODEL.DAYMAX ;
  int    ifilt, EP, ILAM, jflux ;
  double FRAC, LAMSED, LAMDIF, FLUX, FLUXERR, magerr = 0.10 ;
  //  char fnam[] = "get_magerr_SEDMODEL";

  // -------------- BEGIN --------------

  // get epoch index for Trest
  if ( Trest <= DAYMIN ) 
    { EP = 0 ; }
  else if ( Trest >= DAYMAX ) 
    { EP = NDAY-1 ; }
  else
    { get_DAYBIN_SEDMODEL(ISED, Trest, &EP, &FRAC); }

  // get lambda bin
  ifilt     = IFILTMAP_SEDMODEL[ifilt_obs] ;
  LAMSED    = FILTER_SEDMODEL[ifilt].mean / (1.0 + z) ;
  LAMDIF    = LAMSED - SEDMODEL.LAMMIN[ISED] ;
  ILAM      = (int)(LAMDIF/TEMP_SEDMODEL.LAMSTEP) ;
  if ( ILAM < 0        ) { ILAM = 0 ; }
  if ( ILAM > NLAM - 1 ) { ILAM = NLAM - 1 ; }

  jflux      = NLAM*EP + ILAM ;
  FLUX       = TEMP_SEDMODEL.FLUX[jflux];
  FLUXERR    = TEMP_SEDMODEL.FLUXERR[jflux];

  if ( FLUX > 1.0E-15 ) 
    { magerr  = 2.5*log10(1.0+FLUXERR/FLUX); }
  else
    { magerr = 5.0 ; }


  return(magerr) ;

} // end get_magerr_SEDMODEL


// *******************************************
int IFILTSTAT_SEDMODEL(int ifilt_obs, double z) { 

  // Created Aug 26, 2008 by R.Kessler

  // returns 1 if flux can be defined in  observer-filter;
  // returns 0 otherwise
  // Allows checking for valid filters before calling
  // genmag_XXXX and getting an abort.
  //
  // Jun 9, 2011: also check min and max to make sure that
  //              entire transmission is defined.
  //

  int ifilt;
  double LMEAN, LMIN, LMAX, z1;
  //  char fnam[] = "IFILTSTAT_SALT2";

  // ---------- BEGIN ----------

  // translate absolute filter index into sparse index

  z1 = 1. + z;
  ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
  LMEAN = FILTER_SEDMODEL[ifilt].mean / z1 ;
  LMIN  = FILTER_SEDMODEL[ifilt].lammin / z1 ;
  LMAX  = FILTER_SEDMODEL[ifilt].lammax / z1 ;

  //make sure filter-mean is in range
  if ( LMEAN < SEDMODEL.RESTLAMMIN_FILTERCEN ) { return 0 ; }
  if ( LMEAN > SEDMODEL.RESTLAMMAX_FILTERCEN ) { return 0 ; }

  // make sure that entire filter-lambda range is covered by SED.
  if ( LMIN  < SEDMODEL.LAMMIN_ALL  ) { return 0 ; }
  if ( LMAX  > SEDMODEL.LAMMAX_ALL  ) { return 0 ; }

  // if we get here, then this filter is OK
  return 1 ;


}  // end  of IFILTSTAT_SEDMODEL


void get_LAMTRANS_SEDMODEL(int ifilt, int ilam, double *LAM, double *TRANS) {

  // Created Mar 23 2021
  // For input ifilt and ilam, return LAM and TRANS.
  // Use different array for SPECTROGRAPH to hold more wave bins.

  if ( ifilt == JFILT_SPECTROGRAPH ) {
    *LAM   = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] ;
    *TRANS = 1.0 ;
  }
  else {
    *LAM   = FILTER_SEDMODEL[ifilt].lam[ilam];
    *TRANS = FILTER_SEDMODEL[ifilt].transSN[ilam];
  }

}// end get_LAMTRANS_SEDMODEL

// ==============================================
void get_LAMRANGE_SEDMODEL(int opt, double *lammin, double *lammax) {
  // this function is used by external programs to 
  // return LAMRANGE of filter-center.
  // opt=1 -> return range of central filter
  // opt=2 -> return wider range of SED
  
  char fnam[] = "get_LAMRANGE_SEDMODEL" ;

  if ( opt == 1 ) {
    *lammin = SEDMODEL.RESTLAMMIN_FILTERCEN ;
    *lammax = SEDMODEL.RESTLAMMAX_FILTERCEN ;
  }
  else if ( opt == 2 ) {
    *lammin = SEDMODEL.LAMMIN_ALL ;
    *lammax = SEDMODEL.LAMMAX_ALL ;
  }
  else {
    sprintf(c1err, "Invalid opt = %d", opt );
    sprintf(c2err, "Cannot set wavelength range.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return;

} // end of get_LAMRANGE_SEDMODEL


// *************************************
void  checkLamRange_SEDMODEL(int ifilt, double z, char *callFun) {

  // May 2013:
  // Abort if any part of the filter transmission lies outside
  // the SED lambda range.
  // [Aug 2013: moved from genamg_SALT2]
  //
  // Jan 2014: add new input arg *callFun to include name
  //          of calling function in the abort message.
  //

  int    LAMFAIL_LO, LAMFAIL_HI ;
  double LAMREST_MIN, LAMREST_MAX, z1 ;
  char  *cfilt ;

  char fnam[] = "checkLamRange_SEDMODEL" ;

  // ------------------ BEGIN -----------------
  
  z1 = 1.0 + z ;

  LAMREST_MIN  = FILTER_SEDMODEL[ifilt].lammin / z1 ;
  LAMREST_MAX  = FILTER_SEDMODEL[ifilt].lammax / z1 ;

  LAMFAIL_LO   = (LAMREST_MIN < SEDMODEL.LAMMIN_ALL ) ;
  LAMFAIL_HI   = (LAMREST_MAX > SEDMODEL.LAMMAX_ALL ) ;

  if ( LAMFAIL_LO == 0 && LAMFAIL_HI == 0 ) { return ; }

  // if we get here then abort
  cfilt  = FILTER_SEDMODEL[ifilt].name ;
  print_preAbort_banner(fnam);
  printf("  Calling function '%s' passed ifilt=%d('%s') and z=%.4f \n", 
	 callFun, ifilt, cfilt, z);

  printf("  Rest-frame Lambda-Range (LAMREST_MIN/MAX)  = %f to %f \n",
	 LAMREST_MIN, LAMREST_MAX) ;

  printf("  SED Lambda-Range: %f to %f \n",
	 SEDMODEL.LAMMIN_ALL, SEDMODEL.LAMMAX_ALL );

  printf("  LAMFAIL_LO/HI = %d, %d \n", LAMFAIL_LO, LAMFAIL_HI);

  sprintf(c1err,"%s lambda range is outside defined SED range.",
	  cfilt);
  sprintf(c2err,"Consider using &FITINP variable RESTLAMBDA_FITRANGE.");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 

  return ;

} // end of   checkLamRange_SEDMODEL


// ==============================================
void get_DAYBIN_SEDMODEL(int ISED, double DAY, int *IDAY, double *FRAC) {

  // Created Aug 2017
  // For input DAY, return IDAY bin and *FRAC within bin.
  // Use SEDMODEL struct as additional inputs.
  // If DAY-grid is uniform, use DAYSTEP to compute IDAY and FRAC.
  // If DAY-grid is not uniform, use explicit DAY array.
  //

  int    USE_DAYSTEP = 
    ( (SEDMODEL.OPTMASK & OPTMASK_DAYLIST_SEDMODEL)==0 );
  int    NDAY     = SEDMODEL.NDAY[ISED] ;
  double DAYMIN   = SEDMODEL.DAYMIN[ISED] ;
  //  double DAYMAX   = SEDMODEL.DAYMAX[ISED] ;
  double DAYSTEP  = SEDMODEL.DAYSTEP[ISED] ;

  int    IDAY_LOCAL;
  double FRAC_LOCAL, DAYREF0=0.0, DAYREF1=0.0 ;
  char fnam[] = "get_DAYBIN_SEDMODEL" ;

  // --------------- BEGIN ------------

  *IDAY = -9;  *FRAC = -9.0 ;

  if ( NDAY == 1 ) 
    { IDAY_LOCAL = 0 ;  FRAC_LOCAL = 0.0 ;  goto DONE ;   }


  if ( USE_DAYSTEP  ) {
    // use fixed DAYSTEP 
    IDAY_LOCAL = (int)( (DAY - DAYMIN) / DAYSTEP );
    DAYREF0    = DAYMIN + (double)IDAY_LOCAL * DAYSTEP;
    FRAC_LOCAL = (DAY - DAYREF0)/DAYSTEP ;
  }
  else {
    // find where DAY lands in SEDMODEL.DAY array
    IDAY_LOCAL = quickBinSearch(DAY, NDAY, SEDMODEL.DAY[ISED], fnam);
    DAYREF0    = SEDMODEL.DAY[ISED][IDAY_LOCAL];
    DAYREF1    = SEDMODEL.DAY[ISED][IDAY_LOCAL+1];
    DAYSTEP    = DAYREF1-DAYREF0 ;
    FRAC_LOCAL = (DAY - DAYREF0)/DAYSTEP ;
  }

  if ( FRAC_LOCAL < -0.0001 || FRAC_LOCAL > 1.0001 ) {
    sprintf(c1err,"Invalid  FRAC_LOCAL=%le  at IDAY_LOCAL=%d",
	    FRAC_LOCAL, IDAY_LOCAL ) ;
    sprintf(c2err,"DAY=%f  DAYREF = %f,%f", 
	    DAY, DAYREF0, DAYREF1 ) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
 DONE: 
  *IDAY = IDAY_LOCAL ;
  *FRAC = FRAC_LOCAL ;

  return ;

} // end get_DAYBIN_SEDMODEL


// ======== MANGLED FUNCTIONS for FORTRAN ==================

int reset_sedmodel__(void) {
  int istat;
  istat = reset_SEDMODEL();
  return istat;
}

int init_primary_sedmodel__(char *refname, int *NLAM, 
			    double *LAMLIST, double *FLUXLIST ) {
  int istat;
  istat = init_primary_SEDMODEL(refname, *NLAM, LAMLIST, FLUXLIST );
  return istat;
}


// =======================================================
int init_filter_sedmodel__(int *ifilt_obs, char *filter_name, 
			   char *survey_name, double *magprimary,
			   int *NLAM,  double *LAMLIST, 
			   double *TRANSSNLIST, 
			   double *TRANSREFLIST, 
			   double *LAMSHIFT) {
  int istat;
  istat = init_filter_SEDMODEL(*ifilt_obs, filter_name, survey_name,
			       *magprimary, *NLAM, 
			       LAMLIST, TRANSSNLIST,TRANSREFLIST, *LAMSHIFT );
  return istat;
}

void init_redshift_sedmodel__(int *NZbin, double *Zmin, double *Zmax) {
  init_redshift_SEDMODEL(*NZbin, *Zmin, *Zmax);
}

void init_mwxt_sedmodel__(int *OPT_COLORLAW, double *RV) {
  init_MWXT_SEDMODEL(*OPT_COLORLAW,*RV)  ;
}

void get_lamrange_sedmodel__(int *opt, double *lammin, double *lammax) {
  get_LAMRANGE_SEDMODEL(*opt,lammin,lammax);
}



// ==================================================
void pack_SEDBINARY(int OPT) {

  // Created Jan 13, 2010 by R.Kessler
  // OPT = +1 => transfer SEDMODEL to SEDBINARY array
  // OPT = -1 => transfer SEDBINARY array to SEDMODEL
  // Handle variables read by rd_sedFlux().
  //
  // April 12 2018
  //  + write IVERSION 
  //  + write PADWORD after header info
  //  + write FLUX_SCALE in header
  //  + FLUX *= FLUXSCALE to avoid float>E39
  //
  // April 30 2018
  //   + compress long list of zeros using ZEROLIST_SEDBINARY
  //
  // Jun 6 2022: fix aweful index bug restoring flux from SEDBINARY
  //           (N++ was after instead of before)

  int N, NZLEN, NZLEN_LAST, NFLUX, j, IVERSION ;
  double tmpFlux, FLUXSCALE_LOCAL, PADWORD ;
  char fnam[] = "pack_SEDBINARY" ;

  // -------------- BEGIN ----------

  if ( OPT > 0 ) {
    N=0;

    IVERSION = IVERSION_SEDBINARY; // Apr 12 2018
    PADWORD  = PADWORD_SEDBINARY ;
    N++; SEDBINARY[N] = (float)IVERSION           ; // 4/12/2018
    N++; SEDBINARY[N] = (float)TEMP_SEDMODEL.NDAY ;
    N++; SEDBINARY[N] = (float)TEMP_SEDMODEL.DAYSTEP ;
    N++; SEDBINARY[N] = (float)SEDMODEL.FLUXSCALE ; // 4/12/2018 
    N++; SEDBINARY[N] = (float)PADWORD ;            // 4/12/2018

    for ( j=0; j < TEMP_SEDMODEL.NDAY; j++ )
      { N++;  SEDBINARY[N] = TEMP_SEDMODEL.DAY[j]; }

    N++; SEDBINARY[N] = (float)TEMP_SEDMODEL.NLAM ;
    N++; SEDBINARY[N] = (float)TEMP_SEDMODEL.LAMSTEP ;
    for ( j=0; j < TEMP_SEDMODEL.NLAM; j++ )
      { N++;  SEDBINARY[N] = TEMP_SEDMODEL.LAM[j]; }

    NZLEN = NZLEN_LAST = 0; // number of consecutive zeros
    NFLUX = TEMP_SEDMODEL.NDAY * TEMP_SEDMODEL.NLAM ;
    for ( j=0; j<NFLUX; j++ ) { 
      tmpFlux = (TEMP_SEDMODEL.FLUX[j]*SEDMODEL.FLUXSCALE); 
      if ( tmpFlux < 0.0 ) { tmpFlux = 0.0 ; }
      if ( fabs(tmpFlux) < 1.0E-30 ) { tmpFlux=0.0; }
      N++ ; SEDBINARY[N] = (float)tmpFlux;
    } 


    /*
    // xxxxxxxxxxxxxxxx
    float XMEM = (float)N * 4/1.0E3;
    printf(" xxx %s: NDAY=%d  NLAM=%d  NFLUX=%d  N=%d (%.1f KB)\n",
	   fnam, TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.NLAM, NFLUX, N, XMEM );
    // xxxxxxxxxxxxxxxx
    */

    NSEDBINARY = N;
  }
  else {
    N=0;

    N++; IVERSION              = (int)SEDBINARY[N];
    N++; TEMP_SEDMODEL.NDAY    = (int)SEDBINARY[N];
    N++; TEMP_SEDMODEL.DAYSTEP = (double)SEDBINARY[N];
    N++; FLUXSCALE_LOCAL       = (double)SEDBINARY[N];
    N++; PADWORD               = (double)SEDBINARY[N]; // should be 77777

    if ( PADWORD != PADWORD_SEDBINARY ) {
      sprintf(c1err,"PADWORD=%.3f but expected %.3f",
	      PADWORD, PADWORD_SEDBINARY);
      sprintf(c2err,"IVERSION=%d. Try re-making BINARY file.", 
	      IVERSION);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    for ( j=0; j < TEMP_SEDMODEL.NDAY; j++ )
      { N++;  TEMP_SEDMODEL.DAY[j] = SEDBINARY[N]; }

    N++; TEMP_SEDMODEL.NLAM    = (int)SEDBINARY[N];
    N++; TEMP_SEDMODEL.LAMSTEP = SEDBINARY[N];
    for ( j=0; j < TEMP_SEDMODEL.NLAM; j++ )
      { N++;  TEMP_SEDMODEL.LAM[j] = SEDBINARY[N]; }

    NFLUX = TEMP_SEDMODEL.NDAY * TEMP_SEDMODEL.NLAM ;
    for ( j=0; j<NFLUX; j++ ) { 
      N++; // bug fix, Jun 2022
      tmpFlux = (double)SEDBINARY[N];
      // apply 1/flux-scale to return to originally stored flux
      // xxx mark delete N++ ; 
      TEMP_SEDMODEL.FLUX[j] = tmpFlux/FLUXSCALE_LOCAL ; 
    } 

    if ( N != NSEDBINARY ) {
      sprintf(c1err,"Unpacked %d words from SEDBINARY", N);
      sprintf(c2err,"But expected to unpack %d words", NSEDBINARY);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  return ;

} // end of pack_SEDBINARY


// =======================================================
int get_SEDMODEL_INDICES( int IPAR, double LUMIPAR, 
			 int *ILOSED, int *IHISED) {

  // for input paramater IPAR with value LUMIPAR,
  // return SED indices (I0SED and I1SED) that 
  // bracket this LUMIPAR value.
  //

  int ised;
  double     parval ,dif1, mindif1, dif0, mindif0  ;
  //  char fnam[] = "get_SEDMODEL_INDICES" ;

  // ------------ BEGIN ----------------

  *ILOSED = *IHISED = -9 ;
  mindif0 = mindif1 = 9999999.;

  for ( ised = 1; ised <= SEDMODEL.NSURFACE ; ised++ ) {
    parval = SEDMODEL.PARVAL[ised][IPAR];

    dif1 = parval - LUMIPAR;
    if ( dif1 >= 0.0  &&  dif1 < mindif1 ) {
      mindif1 = dif1 ;
      *IHISED = ised ;
    }
    dif0 = LUMIPAR - parval ;
    if ( dif0 >= 0.0  &&  dif0 < mindif0 ) {
      mindif0 = dif0 ;
      *ILOSED = ised ;
    }
  }


  if ( *ILOSED <= 0 || *IHISED <= 0 )
    { return ERROR; }
  else
    { return SUCCESS ; }


} // end of get_SIMSED_INDICES


// ============================================================
double gridval_SIMSED(int ipar, int ibin) {
  
  // Return PARVAL for ipar and bin 'ibin'.
  // ibin = 0, 1, 2 ...
  // If ibin exceeds NBIN for this ipar, then take fmod
  // so that the parameter bin wraps around.

  double BIN, XN, PMIN,  PMAX, PARVAL, xbin, x ;
  //  char fnam[] = "gridval_SIMSED" ;

  // ----------- BEGIN -------------

  PMIN   = SEDMODEL.PARVAL_MIN[ipar] ;
  PMAX   = SEDMODEL.PARVAL_MAX[ipar] ;
  BIN    = SEDMODEL.PARVAL_BIN[ipar] ; 
  XN     = (double)SEDMODEL.NBIN_PARVAL[ipar]  ;

  xbin   = (double)ibin ; 
  x      = fmod(xbin,XN) ;
  PARVAL = PMIN + x*BIN;

  if ( PARVAL <= PMIN ) { PARVAL = PMIN + 1.0E-6*(PMAX-PMIN) ; }
  if ( PARVAL >= PMAX ) { PARVAL = PMAX - 1.0E-6*(PMAX-PMIN) ; }
  
  /*  
  printf(" xxxx ibin = %2d  XN=%2.0f  x=%f  PARVAL=%f \n",
	 ibin, XN, x, PARVAL );
  */

  return PARVAL ;

} // end of gridval_SIMSED


// ============================================================
double nearest_gridval_SIMSED (int ipar, double lumipar ) {

  // Created  Jun 2010
  // return value of parval at grid node nearest to lumipar.
  //
  // July 28 2017: bugfix: check for invalid I0SED and I1SED

  int I0SED, I1SED, istat, NSED ;
  double  parval0, parval1, frac, lumigrid ;

  // -------------- BEGIN --------------

  // get SED indices that bound the input lumipar
  istat = get_SEDMODEL_INDICES( ipar, lumipar, &I0SED, &I1SED ); 

  // check for invalid ISED index
  if ( I0SED < 0 ) { I0SED = I1SED; }
  if ( I1SED < 0 ) { I1SED = I0SED; }

  // extract parval at the nearest nodes.
  parval0 = SEDMODEL.PARVAL[I0SED][ipar]; 
  parval1 = SEDMODEL.PARVAL[I1SED][ipar];

  // snap to closest parval
  frac    = (lumipar - parval0)/(parval1-parval0);
  if ( frac < 0.5 ) 
    { lumigrid = parval0 ; }
  else
    { lumigrid = parval1 ; }

  /*
  printf(" xxxx --------------------------------- \n" );
  printf(" xxxx DUMP nearest_gridval_SIMSED \n");
  printf(" xxxx lumipar=%6.3f  ISED[0,1]=%d,%d  parval[0,1] = %6.3f,%6.3f \n", 
	 lumipar, I0SED, I1SED, parval0, parval1 );
  printf(" xxxx frac=%f  lumigrid = %f \n", frac, lumigrid );
  */

  NSED = SEDMODEL.NSURFACE ; 

  double range = SEDMODEL.PARVAL[NSED][ipar] - SEDMODEL.PARVAL[1][ipar] ;
  if ( lumigrid <= SEDMODEL.PARVAL[1][ipar] )    { lumigrid += range*1.0E-8 ; }
  if ( lumigrid >= SEDMODEL.PARVAL[NSED][ipar] ) { lumigrid -= range*1.0E-8 ; }
  
  return(lumigrid) ;

} // end of nearest_gridval_SIMSED


// ============================================================
void check_sedflux_bins(int ised        // (I) sed index
			,char *VARNAME  // (I) name of variable to check
			,int NBIN       // (I) Number of bins
			,double VAL0    // (I) first value
			,double BINSIZE // (I) bin size
			)  {

  /*
   Created Apr 24, 2010 by R.Kessler
   Store binning from 1st SED, and check for consistency
   for each successive SED. Abort if binning changes.

   Note that ISED can start at 0 or 1, so explicitly count
   NSED locally to know which one is first.
  */

  char fnam[] = "check_sedflux_bins" ;

  int ivar, IVAR, NVAR;
  char *ptrvar;

  int    NBIN_REF;
  double VAL0_REF, BINSIZE_REF;
  char   *ptrvar_ref;

  // -------------- BEGIN -------------

  NVAR = NVAR_FIRSTBIN_SEDMODEL;

  // first check if this is a new or old VARNAME
  IVAR = 0;

  if ( NVAR > 100 ) {
    sprintf(c1err,"Invalid NVAR_FIRSTBIN_SEDMODEL=%d for", NVAR);
    sprintf(c2err,"VARNAME='%s  ised=%d  NBIN=%d", VARNAME, ised, NBIN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  for ( ivar = 1; ivar <= NVAR; ivar++ ) {
    ptrvar = FIRSTBIN_SEDMODEL[ivar].VARNAME;
    if ( strcmp(VARNAME,ptrvar) == 0 ) { IVAR = ivar; }
  }

  // if we have a new variable, increment NVAR
  if ( IVAR == 0 )  { 
    NVAR++ ; 
    IVAR = NVAR; 
    FIRSTBIN_SEDMODEL[IVAR].NSED = 1;
  }
  else
    FIRSTBIN_SEDMODEL[IVAR].NSED++ ;


  if ( FIRSTBIN_SEDMODEL[IVAR].NSED == 1 ) {
    /*
    printf(" xxx store IVAR=%d VARNAME=%s : NBIN=%d VAL0=%f BINSIZE=%f \n",
	   IVAR, VARNAME, NBIN, VAL0, BINSIZE );
    */

    sprintf(FIRSTBIN_SEDMODEL[ivar].VARNAME,"%s", VARNAME);
    FIRSTBIN_SEDMODEL[IVAR].NBIN    = NBIN;
    FIRSTBIN_SEDMODEL[IVAR].VAL0    = VAL0;
    FIRSTBIN_SEDMODEL[IVAR].BINSIZE = BINSIZE;
  }
  else {
    ptrvar_ref  =  FIRSTBIN_SEDMODEL[IVAR].VARNAME ;
    NBIN_REF    =  FIRSTBIN_SEDMODEL[IVAR].NBIN ;
    VAL0_REF    =  FIRSTBIN_SEDMODEL[IVAR].VAL0;
    BINSIZE_REF =  FIRSTBIN_SEDMODEL[IVAR].BINSIZE;

    if ( NBIN != NBIN_REF ) {
      sprintf(c1err,"NBIN(%s)=%d at ISED=%d .", ptrvar_ref, NBIN, ised);
      sprintf(c2err,"but expected NBIN=%d from 1st SED", NBIN_REF);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( VAL0 != VAL0_REF ) {
      sprintf(c1err,"VAL0(%s)=%f at ISED=%d .", ptrvar_ref, VAL0, ised);
      sprintf(c2err,"but expected VAL0=%f from 1st SED", VAL0_REF);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( BINSIZE != BINSIZE_REF ) {
      sprintf(c1err,"BINSIZE(%s)=%f at ISED=%d .", ptrvar_ref, BINSIZE, ised);
      sprintf(c2err,"but expected BINSIZE=%f from 1st SED", BINSIZE_REF);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  }

  NVAR_FIRSTBIN_SEDMODEL = NVAR;

} // end of check_sedflux_bins


void check_surveyDefined_SEDMODEL(void) {

  // Nov 24 2020
  // Abort if survey is not defined for any filter.

  int ifilt, NERR=0;
  char *survey, *name ;
  char fnam[] = "check_surveyDefined_SEDMODEL";

  // ---------- BEGIN ---------

  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
    name   = FILTER_SEDMODEL[ifilt].name  ;
    survey = FILTER_SEDMODEL[ifilt].survey ;
    if (  IGNOREFILE(survey) ) { 
      printf(" ERROR: missing SURVEY name for filter = %s\n", name);
      fflush(stdout);
      NERR++; 
    }
  }

  if ( NERR > 0 ) {
    sprintf(c1err,"Missing survey name for %d of %d filters (see above).",
	    NERR, NFILT_SEDMODEL);
    sprintf(c2err,"Add SURVEY key(s) to kcor/calib input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  return;

} // end check_surveyDefined_SEDMODEL



// ***********************************************
void fill_TABLE_MWXT_SEDMODEL(double RV, double mwebv) {

  // [Aug 23 2013: moved from genmag_SALT2.c]
  // Fill SEDMODEL_TABLE_MWXT_FRAC[ifilt][ilam] = flux-fraction
  // through our Galaxy. Fill lambda-grid on same lambda-grid
  // as the filters. The table is filled only if mwebv has changed
  // so that the table is fixed within each SN processed.
  // 
  // July 24 2016: if spectrograph option is set, load IFILT_SPECTROGRPAPH
  //

  int  NLAMFILT, NBSPEC, ilam, I8, I8p, ifilt, ifilt_min ;
  int  OPT_COLORLAW ;
  double LAMOBS, AV, XT_MAG, XT_FRAC, arg    ;

  //  char fnam[] = "fill_TABLE_MWXT_SEDMODEL";
  
  // ------------- BEGIN ------------------

  I8  = sizeof(double) ;
  I8p = sizeof(double*) ;

  // allocate memory for each filter only once
  if ( SEDMODEL_MWEBV_LAST < 0.0 ) {
    SEDMODEL_TABLE_MWXT_FRAC  = (double**)malloc(I8p*(NFILT_SEDMODEL+1)); 
    for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
      NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
      SEDMODEL_TABLE_MWXT_FRAC[ifilt]  = (double*)malloc(I8*NLAMFILT); 
    }


    // check optional ifilt=0 for spectrograph
    NBSPEC =  FILTER_SEDMODEL[JFILT_SPECTROGRAPH].NLAM ;
    if ( NBSPEC > 0 ) {
      SEDMODEL_TABLE_MWXT_FRAC[JFILT_SPECTROGRAPH] = 
	(double*)malloc(I8*NBSPEC);
    }
  }

  if ( mwebv == SEDMODEL_MWEBV_LAST ) { return ; }

  AV  = RV * mwebv ; 
  OPT_COLORLAW = MWXT_SEDMODEL.OPT_COLORLAW ;
  NBSPEC = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;
  if ( NBSPEC>0 ) { ifilt_min=0; } else { ifilt_min=1; }

  for(ifilt=ifilt_min; ifilt <= NFILT_SEDMODEL; ifilt++) {
   
    NLAMFILT = FILTER_SEDMODEL[ifilt].NLAM ; 

    for ( ilam=0; ilam < NLAMFILT; ilam++ ) {
      LAMOBS     = FILTER_SEDMODEL[ifilt].lam[ilam] ;
      XT_MAG     = GALextinct ( RV, AV, LAMOBS, OPT_COLORLAW ) ;
      arg        = -0.4*XT_MAG ;
      XT_FRAC    = pow(TEN,arg);    // flux-fraction thru MW
      SEDMODEL_TABLE_MWXT_FRAC[ifilt][ilam]  = XT_FRAC ;

    } // ilam
  } // ifilt
  

  // --------
  SEDMODEL_MWEBV_LAST = mwebv ;

  int LDMP = 0; // ( NBSPEC > 0 ) ;
  if ( LDMP ) {
    ifilt = 0 ;  ilam=102 ; 
    if ( NBSPEC == 0 ) { LAMOBS = FILTER_SEDMODEL[ifilt].lam[ilam] ; }
    if ( NBSPEC >  0 ) { LAMOBS = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] ; }
    XT_FRAC = SEDMODEL_TABLE_MWXT_FRAC[ifilt][ilam] ;
    printf(" xxx AV=%5.3f, ifilt=%d   MWXT_FRAC[lam=%6.1f] = %.3f \n", 
	   AV, ifilt, LAMOBS,  XT_FRAC ) ;
    fflush(stdout);
  }

  return ;

} // end of fill_TABLE_MWXT_SEDMODEL


// ***********************************************
void fill_TABLE_HOSTXT_SEDMODEL(double RV, double AV, double z) {

  // Created July 2016
  //
  // Analogous to fill_TABLE_MWXT_SEDMODEL (for Galactic extinct),
  // Fill SEDMODEL_TABLE_HOSTXT_FRAC[ifilt][ilam] = flux-fraction
  // through host Galaxy. Fill lambda-grid on same lambda-grid
  // as the filters. The table is filled only if AV or z has changed
  // so that the table is fixed within each SN processed.
  //
  // Beware that extinction is computed for rest-frame wavelength,
  // but stored as a function of observer-frame wavelength.
  //
  // Mar 2 2017: fix to work with SPECTROGRAPH 
  // Mar 18 2020: DJB added logic for changing RV
  // Jul 13 2020: 
  //   + return if !update_hostxt (bug from v10_76c, Mar 2020)
  //     fixes silly bug causing fit to be almost x10 slower.
  //

  int  NLAMFILT, ilam, I8, I8p, ifilt, ifilt_obs, ifilt_min ;
  int  OPT_COLORLAW, NBSPEC ;

  double  LAMOBS, LAMREST, XT_MAG, XT_FRAC, arg    ;

  char fnam[] = "fill_TABLE_HOSTXT_SEDMODEL";
  
  // ------------- BEGIN ------------------
  
  if ( AV < 0.0 || RV < 0.0 ) { return ; }

  I8  = sizeof(double) ;
  I8p = sizeof(double*) ;

  // allocate memory for each filter only once
  if ( SEDMODEL_HOSTXT_LAST.AV < 0.0 ) {
    SEDMODEL_TABLE_HOSTXT_FRAC  = (double**)malloc(I8p*(NFILT_SEDMODEL+1)); 
    for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
      NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
      SEDMODEL_TABLE_HOSTXT_FRAC[ifilt]  = (double*)malloc(I8*NLAMFILT); 
    }

    // check optional ifilt=0 for spectrograph
    NBSPEC =  FILTER_SEDMODEL[JFILT_SPECTROGRAPH].NLAM ;
    if ( NBSPEC > 0 ) {
      SEDMODEL_TABLE_HOSTXT_FRAC[JFILT_SPECTROGRAPH] = 
	(double*)malloc(I8*NBSPEC);
    }

  }

  bool update_hostxt = false;
  if ( AV != SEDMODEL_HOSTXT_LAST.AV ) { update_hostxt = true; }
  if ( RV != SEDMODEL_HOSTXT_LAST.RV ) { update_hostxt = true; }
  if ( z  != SEDMODEL_HOSTXT_LAST.z  ) { update_hostxt = true; }
  if ( !update_hostxt ) { return; } // put back, July 13 2020

  OPT_COLORLAW = MWXT_SEDMODEL.OPT_COLORLAW ;
  NBSPEC = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;
  if ( NBSPEC>0 ) { ifilt_min=0; } else { ifilt_min=1; }

  for(ifilt=ifilt_min; ifilt <= NFILT_SEDMODEL; ifilt++) {
    NLAMFILT  = FILTER_SEDMODEL[ifilt].NLAM ;
    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs ;

    for ( ilam=0; ilam < NLAMFILT; ilam++ ) {
      LAMOBS     = FILTER_SEDMODEL[ifilt].lam[ilam] ;
      LAMREST    = LAMOBS/(1.0 + z);
      XT_MAG     = GALextinct ( RV, AV, LAMREST, OPT_COLORLAW ) ;
      arg        = -0.4*XT_MAG ;
      XT_FRAC    = pow(TEN,arg);    // flux-fraction thru host
      SEDMODEL_TABLE_HOSTXT_FRAC[ifilt][ilam]  = XT_FRAC ;
    } // ilam

  } // ifilt

  SEDMODEL_HOSTXT_LAST.RV = RV ; // 7.14.2020
  SEDMODEL_HOSTXT_LAST.AV = AV ;
  SEDMODEL_HOSTXT_LAST.z  = z ;


  return ;

} // end fill_TABLE_HOSTXT_SEDMODEL

// ************************************
double filterTrans_BessB(double lam) {

  // Aug 25 2013: moved from SIMSED_fudge.c to here.
  //
  // Return approximate Bessell B transmission at this lambda.
  // Transmission is energy fraction, not count fraction.
  // Hard-wire the B-transmission vs. lambda to avoid more I/O.

#define NBIN_Btrans 21
  double Btrans_LAMBDA[NBIN_Btrans] = { 3600., 3700., 3800., 3900., 4000., 4100., 4200., 4300., 4400., 4500., 4600., 4700., 4800., 4900., 5000., 5100., 5200., 5300., 5400., 5500., 5600. } ;

  double Btrans_VALUE[NBIN_Btrans]  = { 0.0,   0.03,  0.134, 0.567, 0.920, 0.978, 1.00, 0.978, 0.935, 0.853, 0.740, 0.640, 0.536, 0.424, 0.325, 0.235, 0.150, 0.095, 0.043, 0.009, 0.00 } ;


  double trans, lam0, lam1, frac, tr0,tr1 ;
  int ilam;

  char fnam[] = "filterTrans_BessB" ;

  //  ----------- BEGIN -----------

  trans = 0.0 ;

  if ( lam < Btrans_LAMBDA[1] )             return trans ;
  if ( lam > Btrans_LAMBDA[NBIN_Btrans-1] ) return trans ;

  for ( ilam=0; ilam < NBIN_Btrans-1 ;  ilam++ ) {
    lam0 = Btrans_LAMBDA[ilam] ;
    lam1 = Btrans_LAMBDA[ilam+1] ;

    tr0  = Btrans_VALUE[ilam];
    tr1  = Btrans_VALUE[ilam+1];

    if ( lam >= lam0 && lam <= lam1 ) {
      frac = (lam - lam0)/(lam1-lam0);
      trans = tr0 + frac*(tr1-tr0);
      return trans;
    }

  } // ilam

  sprintf(c1err,"Could not evaluate Btrans at lam=%f", lam);
  errmsg(SEV_FATAL, 0, fnam, c1err, "" ); 
  return(-9.0); // should never get here.

} // end of filterTrans_BessB


// =============================================
void T0shiftPeak_SEDMODEL(SEDMODEL_FLUX_DEF *SEDFLUX, int vboseFlag) {

  // Created Aug 11 2015 
  // [moved from shift_NON1A_DAY in genmag_NON1ASED to work for NON1A & SIMSED]
  // Shift TEMP_SEDMODEL.DAY so that peak bolometric flux is at DAY=0.
  // This only affects the MAGT0_[band] parameters and has no effect 
  // on the light curves.  
  //
  // May 15 2018: 
  // + change fun name:  shiftPeakDay_SEDMODEL -> T0shiftPeak_SEDMODEL
  // + pass SEDFLUX struct as arg instead of modifying global struct.
  //
  // July 10 2018: set SEDFLUX->TSHIFT

  int iday, ilam, NDAY, NLAM, IDAY_PEAK, jflux ;
  double FLUXTMP, FLUXSUM, FLUXSUM_MAX, Tshift ;
  char fnam[] = "T0shiftPeak_SEDMODEL" ;

  // ------------- BEGIN ----------

  NDAY = SEDFLUX->NDAY ;
  NLAM = SEDFLUX->NLAM ;
  FLUXSUM_MAX = 0.0 ;
  IDAY_PEAK = -9 ;

  for(iday=0; iday < NDAY; iday++ ) {
    FLUXSUM = 0.0 ;
    for(ilam=0; ilam < NLAM; ilam++ ) {
      jflux      = NLAM*iday + ilam ;
      FLUXTMP    = SEDFLUX->FLUX[jflux];
      FLUXSUM += FLUXTMP ;
    } // end ilam           

    if ( FLUXSUM > FLUXSUM_MAX ) {
      FLUXSUM_MAX = FLUXSUM ;
      IDAY_PEAK   = iday ;
    }
  } // end iday 


  if ( IDAY_PEAK >= 0 ) 
    {  Tshift = -SEDFLUX->DAY[IDAY_PEAK]; }
  else
    { Tshift = 0.0 ; } // all fluxes are zero; temp fix ??

  SEDFLUX->TSHIFT = Tshift; // July 10 2018

  if( vboseFlag ) {
    printf("    %s: shift DAY-array by %.2f days so that "
	   "maxFlux is at DAY=0\n",
	   fnam, Tshift); fflush(stdout);
  }

  for(iday=0; iday < NDAY; iday++ )
    { SEDFLUX->DAY[iday] += Tshift ;  }
  
  return ;

} // end of T0shiftPeak_SEDMODEL

// =================================================
void T0shiftExplode_SEDMODEL(int OPTMASK, SEDMODEL_FLUX_DEF *SEDFLUX, 
			     int vboseFlag) {

  // Created May 15 2018 by R.Kessler
  // Shift DAY array so that T=0 corresponds to explosion time
  // instead of time of peak brightness.
  //
  // Inputs:
  //   OPTMASK   = option mask (see below)
  //   SEDFLUX   = structure to analyse & modify
  //   vboseflag = verbosity flag
  //
  // Outputs:
  //   SEDFLUX->DAY[IDAY]  are modified so that DAY=0 at explosion
  //
  //
  // Option to determine explosion time
  //   OPTMASK & 1  --> Explosion  T=0 when flux *= E-2
  //   OPTMASK & 2  --> Explosion  T=0 when flux *= E-3
  //   OPTMASK & 32 --> add 32 days from peak (idiot test)
  //
  // Jun 4 2018: fix sign error in Tshift.
  //
  int  jflux, NDAY, NLAM, IDAY_MAX=0, ILAM_MAX, iday, ilam ;
  double Tshift, FLUXTMP, FLUXMAX, DAY, FLUX_EXPLODE=0.0 ;
  char fnam[] = "T0shiftExplode_SEDMODEL" ;

  // ---------------- BEGIN ---------------

  if ( OPTMASK == 0 ) { return ; }

  NDAY = SEDFLUX->NDAY ;
  NLAM = SEDFLUX->NLAM ;

  // find max flux for reference
  FLUXMAX = -9.0 ;
  for(iday=0; iday < NDAY; iday++ ) {
    for(ilam=0; ilam < NLAM; ilam++ ) {
      jflux      = NLAM*iday + ilam ;
      FLUXTMP    = SEDFLUX->FLUX[jflux];
      if ( FLUXTMP > FLUXMAX ) {
	FLUXMAX = FLUXTMP;  IDAY_MAX=iday; ILAM_MAX=ilam;
      }
    } // end ilam 
  } // end iday                                                   

  if ( FLUXMAX < 0.0 ) {
    sprintf(c1err, "Could not find FLUXMAX for OPTMASK=%d", OPTMASK);
    sprintf(c2err, "Something really weird??");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( (OPTMASK & 1)>0 ) {
    FLUX_EXPLODE= 0.01 * FLUXMAX ;
  }
  else if ( (OPTMASK & 2)>0 ) {
    FLUX_EXPLODE= 0.001 * FLUXMAX ;
  }
  else if ( (OPTMASK & 32)>0 ) {
    T0shiftPeak_SEDMODEL(SEDFLUX,0);
    Tshift = 32.0 ;
  }
  else {
    sprintf(c1err, "Invalid OPTMASK=%d", OPTMASK);
    sprintf(c2err, "Check options in manual or in genmag_SEDtools.c");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  /*
  printf(" xxx OPTMASK=%d FLUXMAX=%le  FLUX_EXPLODE=%le \n", 
	 OPTMASK, FLUXMAX, FLUX_EXPLODE ); fflush(stdout);
  */

  if ( FLUX_EXPLODE > 1.0E-12 ) {
    Tshift = 9999. ;
    for(ilam=0; ilam < NLAM; ilam++ ) {

      // start at 2nd bin to allow checking previous day bin.
      for(iday=1; iday < IDAY_MAX; iday++ ) { 
	jflux      = NLAM*iday + ilam ;
	FLUXTMP    = SEDFLUX->FLUX[jflux];
	DAY        = SEDFLUX->DAY[iday-1];  // note previous day bin
	if ( FLUXTMP > FLUX_EXPLODE  &&  DAY < Tshift ) {
	  Tshift = DAY;
	}
      } // end ilam 
    } // end iday                                                       
  }  // end FLUX_EXPLODE if block


  if ( vboseFlag ) {
    printf("    Shift %.1f days -> T=0 at explosion (OPTMASK=%d)\n",
	   -Tshift, OPTMASK); fflush(stdout);
  }


  for(iday=0; iday < NDAY; iday++ )
    { SEDFLUX->DAY[iday] -= Tshift ;  }    

  SEDFLUX->TSHIFT = -Tshift; // July 10 2018

  return;

} // end T0shiftExplode_SEDMODEL


// ======================================================
void FLUX_SCALE_SEDMODEL(double SCALE, SEDMODEL_FLUX_DEF *SEDFLUX) {

  // Apr 22 2019
  // multiply every flux by SCALE

  int NDAY  = SEDFLUX->NDAY;
  int NLAM  = SEDFLUX->NLAM;
  int NBTOT = NLAM * NDAY;
  int jflux ;
  // -------------- BEGIN -------------

  for(jflux=0; jflux < NBTOT; jflux++ ) 
    {  SEDFLUX->FLUX[jflux] *= SCALE ;   }
  
  return ;

} // end FLUX_SCALE_SEDMODEL

// ===============================================
bool found_fluxerr_SEDMODEL(char *sedFile) {

  // Created Feb 6 2020
  // Returns true if 4th column (fluxerr) is found in the
  // input SED file.

  bool found = false;
  int  gzipFlag, NRDWORD, i ;
  FILE *fp;
  char line[200], *stringVal[MXWORDLINE_FLUX], space[] = " ";
  char fnam[] = "found_fluxerr_SEDMODEL" ;

  // -------------- BEGIN ------------

  fp = open_TEXTgz(sedFile, "rt", &gzipFlag);
  if ( !fp ) {
    sprintf(c1err,"Could not open sedFile:");
    sprintf(c2err,"%s", sedFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  for(i=0; i < MXWORDLINE_FLUX; i++ )
    { stringVal[i] = (char*) malloc(60*sizeof(char) ) ;  }

  while ( fgets (line, 200, fp ) != NULL  ) {
    if ( strlen(line) < 2  ) { continue ; }
    if ( commentchar(line) ) { continue ; }

    splitString(line, space, MXWORDLINE_FLUX,
                &NRDWORD, stringVal ) ;  // returned                         

    if ( NRDWORD >= 4 ) { found = true; }
    goto DONE;
  }

 DONE:
  // free mem
  for(i=0; i < MXWORDLINE_FLUX; i++ )    { free(stringVal[i]); }
  fclose(fp);
  return(found);

} // end found_fluxerr_SEDMODEL

// ======================================
void print_ranges_SEDMODEL(SEDMODEL_FLUX_DEF *SEDFLUX) {

  // Creatd June 2022 by R.Kessler  
  // Print range of Trest and LAMMIN for SEDFLUX.
  // Also check if average flux at DAYMIN and DAYMAX is below 2% of peak;
  // if edge fluxes are above 2%, give warning.     

  int    NDAY   = SEDFLUX->NDAY;
  int    NLAM   = SEDFLUX->NLAM;
  double DAYMIN = SEDFLUX->DAYMIN;
  double DAYMAX = SEDFLUX->DAYMAX;
  double LAMMIN = SEDFLUX->LAMMIN;
  double LAMMAX = SEDFLUX->LAMMAX;

  int jflux, iday, ilam, iday_sparse;
  int IDAY_PEAK=-9, IDAY_EDGE_MIN=0, IDAY_EDGE_MAX=NDAY-1;
  double day, lam, dif, difmin, FLUXTMP ;

  int    NDAY_SPARSE = 3;
  int    IDAY_LIST[3];
  double FLUXSUM_LIST[3]; // peak, minEdge, maxEdge
  char fnam[] = "print_ranges_SEDMODEL";

  // --------- BEGIN ----------

  printf("\t Trest: %6.2f to %6.2f     LAMBDA: %5.0f to %5.0f A \n"
         ,DAYMIN, DAYMAX, LAMMIN, LAMMAX);
  fflush(stdout) ;

  // assume peak is near DAY=0 and find IDAY_PEAK
  difmin = 99999.0 ;
  for(iday=0; iday < NDAY; iday++ ) {
    day = SEDFLUX->DAY[iday];
    if ( fabs(day) < difmin ) { difmin=fabs(day); IDAY_PEAK=iday; }
  }

  // get bolometric flux at edges and peak
  
  IDAY_LIST[0] = IDAY_PEAK;
  IDAY_LIST[1] = IDAY_EDGE_MIN;
  IDAY_LIST[2] = IDAY_EDGE_MAX;

  for(iday_sparse=0; iday_sparse < NDAY_SPARSE; iday_sparse++ ) {
    iday = IDAY_LIST[iday_sparse];
    FLUXSUM_LIST[iday_sparse] = 0.0 ;

    for(ilam=0; ilam < NLAM; ilam++ ) {
      jflux      = NLAM*iday + ilam ;
      FLUXTMP    = SEDFLUX->FLUX[jflux];
      FLUXSUM_LIST[iday_sparse] += FLUXTMP ;
    } // end ilam           

  } // end iday_sparse

  // - - - - - -
  double EDGE_FLUX_RATIO_WARNING = 0.04 ;
  double fluxratio_edge_min = FLUXSUM_LIST[1] / FLUXSUM_LIST[0];
  double fluxratio_edge_max = FLUXSUM_LIST[2] / FLUXSUM_LIST[0];

  if ( fluxratio_edge_min > EDGE_FLUX_RATIO_WARNING ) {
    printf("\t WARNING: Flux(Trestmin)/Flux(peak) = %.3f -> "
	   "beware of extrapolation\n",
           fluxratio_edge_min);
  }
  if ( fluxratio_edge_max > EDGE_FLUX_RATIO_WARNING ) {
    printf("\t WARNING: Flux(Trestmax)/Flux(peak) = %.3f -> "
	   "beware of extrapolation\n",
           fluxratio_edge_max);
  }

  
  fflush(stdout) ;

  return ;

} // end  print_ranges_SEDMODEL

// ======================================================
void UVLAM_EXTRAPFLUX_SEDMODEL(double UVLAM, SEDMODEL_FLUX_DEF *SEDFLUX) {

  // Created May 21 2018
  //
  // For inputs UV wavelenth "UVLAM" (Angstroms), extrapolate SEDFLUX 
  // down to UVLAM. This routine is called from sim-input key
  //  UVLAM_EXTRAPFLUX:  <UVLAM>
  //
  // Note that SEDFLUX structure is modified.
  // Extrapolted flux is linear starting at LAMMIN_ORIG,
  // and going down to zero ad LAMMIN_NEW.
  //
  // If flux=0 at LAMMIN_ORIG, then find smallest MINLAM with non-zero
  // flux and extrapolate flux from this MINLAM.
  //

  int NDAY           = SEDFLUX->NDAY;
  int NLAM_ORIG      = SEDFLUX->NLAM;
  double LAMMIN_ORIG = SEDFLUX->LAM[0];
  double LAMSTEP     = SEDFLUX->LAMSTEP;

  int    NLAM_NEW   = NLAM_ORIG ;
  double LAMMIN_NEW = LAMMIN_ORIG ;

  double FLUXTMP, *FLUXORIG, LAM_ADD, LAM_EXTRAP, LAM, LAMMIN_TMP ;
  double *LAMMIN_LIST ;  // minLam for each DAY

  int NBTOT_ORIG = NLAM_ORIG * NDAY;
  int NLAM_ADD, NBTOT_NEW, *ILAMMIN_LIST, ILAM_TMP;
  int iday, ilam_new, ilam_orig, jflux_orig, jflux_new ;
  char fnam[] = "UVLAM_EXTRAPFLUX_SEDMODEL" ;

  // -------------- BEGIN ---------------

  if ( UVLAM < 0.0         ) { return ; }
  if ( UVLAM > LAMMIN_ORIG ) { return ; }

  // figure out new LAMMIN on original LAMSTEP grid
  while ( LAMMIN_NEW > UVLAM ) { LAMMIN_NEW -= LAMSTEP ;  NLAM_NEW++ ; }

  NLAM_ADD  = NLAM_NEW - NLAM_ORIG; // number of LAM bins to add
  NBTOT_NEW = NLAM_NEW * NDAY ;
  LAM_ADD   = LAMMIN_ORIG - LAMMIN_NEW;

  printf("    Extrap UVLAM: %.1f -> %.1f A  (NLAM=%d->%d) \n",
	 LAMMIN_ORIG, LAMMIN_NEW,  NLAM_ORIG, NLAM_NEW);
  fflush(stdout);

  // make sure extended lambda range does not exceed bound.
  if ( NLAM_NEW >= MXBIN_LAMSED_SEDMODEL  ) {
    sprintf(c1err, "NLAM_NEW=%d exceeds bound", NLAM_NEW);
    sprintf(c2err, "Check bound MXBIN_LAMSED_SEDMODEL = %d",
	    MXBIN_LAMSED_SEDMODEL  );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // allocate storage for original flux
  FLUXORIG = (double*) malloc( sizeof(double) * NBTOT_ORIG );
  for(jflux_orig=0; jflux_orig < NBTOT_ORIG; jflux_orig++ ) {  
    FLUXORIG[jflux_orig] = SEDFLUX->FLUX[jflux_orig]; 
    SEDFLUX->FLUX[jflux_orig] = 0.0 ;
  }

  // allocate storage for minLam 
  LAMMIN_LIST  = (double*) malloc( NDAY * sizeof(double) );
  ILAMMIN_LIST = (int   *) malloc( NDAY * sizeof(int) );
  
  
  // to make room for UV wavelengths, push current fluxes 
  // to higher LAM indices. 

  for(iday=0; iday < NDAY; iday++ ) {

    LAMMIN_TMP = 1.0E12;  ILAM_TMP=-9;

    for(ilam_orig=0; ilam_orig<NLAM_ORIG ; ilam_orig++ ) {
      jflux_orig    = NLAM_ORIG*iday + ilam_orig ;
      FLUXTMP       = FLUXORIG[jflux_orig];

      ilam_new  = ilam_orig + NLAM_ADD;
      jflux_new = NLAM_NEW*iday + ilam_new ;
      SEDFLUX->FLUX[jflux_new] = FLUXTMP;

      // keep track of minLam with non-zerp flux
      LAM = SEDFLUX->LAM[ilam_orig];
      if ( FLUXTMP > 1.0E-12 && LAM < LAMMIN_TMP ) 
	{ LAMMIN_TMP = LAM; ILAM_TMP=ilam_orig; }

      if ( jflux_new >= NBTOT_NEW ) {
	sprintf(c1err, "jflux_new=%d exceeds bound NBTOT_NEW=%d",
		jflux_new, NBTOT_NEW);
	sprintf(c2err, "iday=%d ilam_orig=%d  ilam_new=%d",
		iday, ilam_orig, ilam_new);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
    } // end ilam

    if ( ILAM_TMP < 0 ) {
      LAMMIN_LIST[iday]  = LAMMIN_ORIG ;
      ILAMMIN_LIST[iday] = 0 ;
    } else {
      LAMMIN_LIST[iday]  = LAMMIN_TMP ;
      ILAMMIN_LIST[iday] = ILAM_TMP;
    }

    /*
    printf(" xxx iday=%2d: MINLAM=%.1f  (ilam=%d) \n",
    iday, LAMMIN_LIST[iday], ILAMMIN_LIST[iday] ); */

  } // end iday


  // update entire wavelength array
  SEDFLUX->NLAM   = NLAM_NEW ;
  for(ilam_new=0; ilam_new<NLAM_NEW ; ilam_new++ ) {
    SEDFLUX->LAM[ilam_new] = LAMMIN_NEW + LAMSTEP*(double)(ilam_new);
  }


  // now fill in UV range for added wavelength range,
  // plus original wavelengths that have zero flux.

  double FLUXEDGE ;
  int NLAM_EXTRAP ;

  for(iday=0; iday < NDAY; iday++ ) {

    LAMMIN_TMP  = LAMMIN_LIST[iday] ;
    NLAM_EXTRAP = NLAM_ADD + ILAMMIN_LIST[iday];
    LAM_EXTRAP  = LAMMIN_TMP - LAMMIN_NEW;

    /* xxxxxxxxx DEBUG ONLY xxxxx
    if ( UVLAM < 499.99999 ) { 
      NLAM_EXTRAP = NLAM_ADD; 
      LAM_EXTRAP  = LAM_ADD ;
      LAMMIN_TMP  = LAMMIN_ORIG ;
    } // xxx REMOVE
    xxxxxxxxxxxxxxxxxxxxxxxxxx */

    jflux_new = NLAM_NEW*iday + NLAM_EXTRAP ;
    FLUXEDGE  = SEDFLUX->FLUX[jflux_new];


    for(ilam_new=0; ilam_new<NLAM_EXTRAP ; ilam_new++ ) {
      jflux_new = NLAM_NEW*iday + ilam_new ;
      LAM       = SEDFLUX->LAM[ilam_new];
      SEDFLUX->FLUX[jflux_new] = FLUXEDGE * (LAM-LAMMIN_NEW)/LAM_EXTRAP ;
    }
  }


  free(FLUXORIG);  free(LAMMIN_LIST);  free(ILAMMIN_LIST);
  return ;

} // end UVLAM_EXTRAPFLUX_SEDMODEL


// ==================================================
void set_UVLAM_EXTRAPFLUX_SEDMODEL(float UVLAM_MIN) {
  // --------------- BEGIN ------------------
  INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX = (double)UVLAM_MIN;
  return ;
} // end set_UVLAM_EXTRAPFLUX_SEDMODEL

void set_uvlam_extrapflux_sedmodel__(float *UVLAM_MIN)
  {  set_UVLAM_EXTRAPFLUX_SEDMODEL(*UVLAM_MIN); }


// =========================================================
// ============== SPECTROGRAPH FUNCTIONS ===================
// =========================================================


void INIT_SPECTROGRAPH_SEDMODEL(char *MODEL_NAME, int NBLAM,
                                double *LAMMIN_LIST, double *LAMMAX_LIST) {

  // Aug 2016
  // called by simulation to initialize stuff for SPECTROGRAPH.
  //
  // Nov 10 2016: fix sorting bug after sortDouble call
  // Jul 12 2019: 
  //  + store FILTER_SEDMODEL[IFILT].lammin/lammax (for BYOSED)
  //

  int  IFILT  = JFILT_SPECTROGRAPH ;
  int  MEMD   = NBLAM * sizeof(double);
  int  ilam, DUMPFLAG_ZP ;
  double L0, L1, LAVG;
  char fnam[] = "INIT_SPECTROGRAPH_SEDMODEL" ;

  // ---------- BEGIN ----------

  sprintf(BANNER, "%s : prep %s spectra for %s (NBLAM=%d)",
          fnam, MODEL_NAME, INPUTS_SPECTRO.INSTRUMENT_NAME, NBLAM );
  print_banner(BANNER);


  SPECTROGRAPH_SEDMODEL.NBLAM_TOT   = NBLAM ;
  SPECTROGRAPH_SEDMODEL.LAMMIN_LIST = (double*) malloc(MEMD) ;
  SPECTROGRAPH_SEDMODEL.LAMMAX_LIST = (double*) malloc(MEMD) ;
  SPECTROGRAPH_SEDMODEL.LAMAVG_LIST = (double*) malloc(MEMD) ;
  SPECTROGRAPH_SEDMODEL.ZP_LIST     = (double*) malloc(MEMD) ;

  for(ilam=0; ilam < NBLAM; ilam++ ) {
    L0 = LAMMIN_LIST[ilam] ;
    L1 = LAMMAX_LIST[ilam] ;
    LAVG = (L0+L1)/2.0 ;
    SPECTROGRAPH_SEDMODEL.LAMMIN_LIST[ilam] = L0 ;
    SPECTROGRAPH_SEDMODEL.LAMMAX_LIST[ilam] = L1 ;
    SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ilam] = LAVG ;

    // get zero point needed to get source mag in each lambin
    DUMPFLAG_ZP = 0; // ( fabs(L0-4210.0) < 2.0 ) ;
    SPECTROGRAPH_SEDMODEL.ZP_LIST[ilam]
      = getZP_SPECTROGRAPH_SEDMODEL(L0,L1,DUMPFLAG_ZP);
  }


  // store global min/max wavelength  
  L0 = SPECTROGRAPH_SEDMODEL.LAMMIN_LIST[0];
  L1 = SPECTROGRAPH_SEDMODEL.LAMMAX_LIST[NBLAM-1];
  SPECTROGRAPH_SEDMODEL.LAMMIN = L0;
  SPECTROGRAPH_SEDMODEL.LAMMAX = L1;

  FILTER_SEDMODEL[IFILT].lammin =   L0;
  FILTER_SEDMODEL[IFILT].lammax =   L1;

  FILTER_SEDMODEL[IFILT].NLAM      = NBLAM ;
  IFILTMAP_SEDMODEL[IFILT]         = IFILT ;
  FILTER_SEDMODEL[IFILT].ifilt_obs = IFILT ;
  FILTER_SEDMODEL[IFILT].lamstep   = 10.0   ;  // anything not zero

  char   *cfilt      = FILTER_SEDMODEL[IFILT].name ;
  sprintf(cfilt, "%s", "SPECTROGRAPH");

  printf("\t %d wavelength bins, %.1f to %.1f A \n", NBLAM,
         SPECTROGRAPH_SEDMODEL.LAMMIN, SPECTROGRAPH_SEDMODEL.LAMMAX );
  fflush(stdout);

  // --------------------------------------------
  // sort primary mags vs. increasing wavelength, to interpolate
  // primary mag at arbitrary wavelength.


  int ifilt, isort, INDEX_SORT[MXFILTINDX] ;
  int ORDER = +1 ;
  int NFILT = NFILT_SEDMODEL ;
  double  LAMMIN_FILT, LAMMAX_FILT, LAM_ARRAY[MXFILTINDX] ;
  char   *name = PRIMARY_SEDMODEL.name ;

  for ( ifilt=1; ifilt <= NFILT; ifilt++ ) {
    LAM_ARRAY[ifilt]            = FILTER_SEDMODEL[ifilt].mean ;
    PRIMARY_SEDMODEL.MAG[ifilt] = FILTER_SEDMODEL[ifilt].magprimary; 
  }

  sortDouble(NFILT, &LAM_ARRAY[1], ORDER, &INDEX_SORT[1] ) ;

  for ( isort=1; isort <= NFILT; isort++ ) {
    ifilt = INDEX_SORT[isort] + 1 ;
    PRIMARY_SEDMODEL.MAG_SORTED[isort] =  PRIMARY_SEDMODEL.MAG[ifilt] ;
    PRIMARY_SEDMODEL.LAM_SORTED[isort] =  LAM_ARRAY[ifilt];
  }

  // print summary of primaryMag interpolation.
  // Get min,max lambda for broad band filters where primary mag
  // is well defined.
  LAMMIN_FILT = PRIMARY_SEDMODEL.LAM_SORTED[1]; 
  LAMMAX_FILT = PRIMARY_SEDMODEL.LAM_SORTED[NFILT];
  
  int NLINE=0, iline ;
  double LAMMIN_PRIMARY[6], LAMMAX_PRIMARY[5];
  char   str_mag[6][40];

  if ( SPECTROGRAPH_SEDMODEL.LAMMIN < LAMMIN_FILT ) {
    LAMMIN_PRIMARY[NLINE] = SPECTROGRAPH_SEDMODEL.LAMMIN ;
    LAMMAX_PRIMARY[NLINE] = LAMMIN_FILT ;
    sprintf(str_mag[NLINE], "%.3f", PRIMARY_SEDMODEL.MAG_SORTED[1]);
    NLINE++ ;
  }

  LAMMIN_PRIMARY[NLINE] = LAMMIN_FILT ;
  LAMMAX_PRIMARY[NLINE] = LAMMAX_FILT ;
  sprintf(str_mag[NLINE], "interpolated value" );
  NLINE++ ;

  if ( SPECTROGRAPH_SEDMODEL.LAMMIN < LAMMIN_FILT ) {
    LAMMIN_PRIMARY[NLINE] = LAMMAX_FILT ;
    LAMMAX_PRIMARY[NLINE] = SPECTROGRAPH_SEDMODEL.LAMMAX ;
    sprintf(str_mag[NLINE], "%.3f", PRIMARY_SEDMODEL.MAG_SORTED[NFILT]);
    NLINE++ ;
  }


  for(iline=0; iline < NLINE; iline++ ) {
    printf("\t magSpecBin(%s, %7.1f - %7.1f A) = %s  \n", 
	   name, LAMMIN_PRIMARY[iline], LAMMAX_PRIMARY[iline],
	   str_mag[iline] );
  }

  return ;

} // end INIT_SPECTROGRAPH_SEDMODEL


double getZP_SPECTROGRAPH_SEDMODEL(double LAMMIN, double LAMMAX, 
				   int DUMPFLAG) {

  // Return zeropoint for spectrograph bin
  // bounded by input LAMMIN & LAMMIN.
  //
  // Jun 18 2021: abort of primary lam range does not cover spectrograph range.

  double ZP, lam0, lamCen, lamStep, fluxSum, flux, magPrimary ;
  double hc8     = (double)hc ;
  double LAMSTEP = SEDMODEL.LAMSTEP[0] ;
  int    LDMP    = DUMPFLAG ;
  int    NLOOP   = 0;
  char   fnam[]  = "getZP_SPECTROGRAPH_SEDMODEL" ;

  // ------------ BEGIN -----------

  // for NON1SED, LAMSTEP is not read yet, so just use 10 A
  if ( LAMSTEP == 0.0  ) { LAMSTEP = 10.0 ; } 

  fluxSum = 0.0 ;  // for primary

  if ( LDMP ) {
    printf(" xxx ------- %s DUMP  --------- \n", fnam);
    printf(" xxx LAMRANGE: %.3f - %.3f   LAMSTEP=%.3f \n",
	   LAMMIN, LAMMAX, LAMSTEP);
    fflush(stdout);
  }
  
  // Jun 18 2021: check that primary spectrum is defined over wave range  
  int NLAM = PRIMARY_SEDMODEL.NLAM;
  double LAMMIN_PRIM = PRIMARY_SEDMODEL.lam[0];
  double LAMMAX_PRIM = PRIMARY_SEDMODEL.lam[NLAM-1];
  if ( LAMMIN < LAMMIN_PRIM || LAMMAX > LAMMAX_PRIM ) {
    print_preAbort_banner(fnam);
    printf("\t Spectrograph bin wave range: %.1f to %.1f \n",
           LAMMIN, LAMMAX);
    printf("\t Primary wave range: %.1f to %.1f \n",
           LAMMIN_PRIM, LAMMAX_PRIM);
    sprintf(c1err,"Primary SED wave range does not cover spectrograph.");
    sprintf(c2err,"Check KCOR-input, and allow for extended spectro bins.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  // - - - - - - -

  for(lam0=LAMMIN; lam0 < LAMMAX; lam0 += LAMSTEP ) {

    // get lambda bin size; beware near edge
    if ( lam0+LAMSTEP < LAMMAX )
      { lamStep = LAMSTEP ; } 
    else
      { lamStep = (LAMMAX-lam0) ; }

    lamCen   = lam0 + lamStep/2.0 ;
    flux     = interp_primaryFlux_SEDMODEL(lamCen) ; 
    fluxSum += (flux * lamCen * lamStep/hc8); 

    if ( LDMP ) {
      printf(" xxx lam0=%.2f  lamStep=%4.2f  lamCen=%.2f flux=%10.3le\n",
	     lam0, lamStep, lamCen, flux); fflush(stdout);
    }
    NLOOP++ ;

    if ( NLOOP > 100 ) {
      sprintf(c1err, "NLOOP=%d is too much", NLOOP);
      sprintf(c2err, "LAMMIN=%.0f  LAMMAX=%.0f  LAMSTEP=%.1f",
	      LAMMIN, LAMMAX, lamStep);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  } // end lam

  lamCen     = 0.5 * ( LAMMIN + LAMMAX );
  magPrimary = interp_primaryMag_SEDMODEL(lamCen);

  if( fluxSum != 0.0 ) 
    { ZP = 2.5*log10(fluxSum) + magPrimary ; }
  else
    { ZP = 0.0 ; }


  if( LDMP ) {
    printf(" xxx fluxSum=%le magPrim=%.3f -->  ZP=%f \n",
	   fluxSum, magPrimary, ZP); fflush(stdout);
  }

  return(ZP);

} // end getZP_SPECTROGRAPH_SEDMODEL

// *****************************************************
void getSpec_SEDMODEL(int ised,
		      double MWEBV, double RV_host, double AV_host,
		      double z, double MU, double Tobs, double MAGOFF,
		      double *GENFLUX_LIST, double *GENMAG_LIST ) {

  // Return integrated flux in each wavelength bin (not FLAM=dF/dlam),
  // and corresponding magnitude.
  // Conversion to Flam is done in GENSPEC_DRIVER.
  //
  // WARNING: Works for NON1ASED, but does NOT work for SIMSED
  // because TEMP_SEDMODEL array is not kept for each SED.
  //
  // Mar  6 2017: pass extinction arguments MWEBV,RV_host,AV_host
  // Dec 12 2018: replace x0 argument with MU; compute x0 below.
  // Mar 29 2019: SEDMODELNORM is separate for MAG and SPEC (hc factor)

  double x0     = pow(TEN,-0.4*MU);
  double hc8    = (double)hc;
  double LAMBIN = 10.0 ;  // Angstromgs, integration binsize
  int    IFILT  = JFILT_SPECTROGRAPH ;
  int    NBSPEC = FILTER_SEDMODEL[IFILT].NLAM ;

  double LAMTMP_OBS, LAMTMP_REST, LAM0, LAM1, LAMAVG, lamBin, lam ;
  double ZP, MAG, FLUXGEN_forSPEC, FLUXGEN_forMAG;
  double FTMP, MWXT_FRAC, z1, x0fac ;
  double SEDNORM_forSPEC, SEDNORM_forMAG;
  int    ispec, NBLAM ;
  int    NDMP_SKIP=0;
  int    LDMP  = ( fabs(Tobs) < -5.0 ) ; 

  char   fnam[] = "getSpec_SEDMODEL" ;

  // --------------- BEGIN ---------------

  if ( ised < 0 ) {
    sprintf(c1err, "Undefined  ised=%d", ised);
    sprintf(c2err, "z=%.3f  Tobs=%.2f", z, Tobs );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  SEDNORM_forMAG  = SEDMODEL.FLUXSCALE / hc8 ;
  SEDNORM_forSPEC = SEDMODEL.FLUXSCALE ;
  z1 = 1.0 + z ;

  if ( LDMP ) {
    printf("\n xxx %s DUMP for ised=%d z=%.3f Tobs=%.2f  NBSPEC=%d \n",
	   fnam, ised, z, Tobs, NBSPEC );
    printf(" xxx FLUXSCALE=%f  SEDNORM(MAG)=%le  x0=%le\n",
	   SEDMODEL.FLUXSCALE, SEDNORM_forMAG, x0 );
  }


  // loop over spectrograph wavelength bins
  for ( ispec=0; ispec< NBSPEC; ispec++ ) {

    GENFLUX_LIST[ispec] = 0.0 ;
    GENMAG_LIST[ispec]  = MAG_UNDEFINED ;

    // get wavelength range in this spectrograph bin
    LAM0      = INPUTS_SPECTRO.LAMMIN_LIST[ispec]; // obs-frame
    LAM1      = INPUTS_SPECTRO.LAMMAX_LIST[ispec];
    LAMAVG    = INPUTS_SPECTRO.LAMAVG_LIST[ispec];
    ZP        = SPECTROGRAPH_SEDMODEL.ZP_LIST[ispec] ;

    NBLAM = (int)((LAM1-LAM0)/LAMBIN) ;
    if ( NBLAM < 3 ) { NBLAM=3; }
    lamBin = (LAM1-LAM0)/(double)NBLAM ;

    // sub loop in finer observer-lambda bins
    FLUXGEN_forSPEC = FLUXGEN_forMAG = 0.0 ;
    for(lam=LAM0; lam < LAM1-0.1 ; lam+=lamBin) {
      LAMTMP_OBS   = lam + lamBin/2.0 ;
      LAMTMP_REST  = LAMTMP_OBS/z1 ;
      FTMP  = getFluxLam_SEDMODEL(ised, -9, Tobs, LAMTMP_OBS, z, fnam) ;
      FLUXGEN_forMAG  += (FTMP * lamBin * LAMTMP_REST); 
      FLUXGEN_forSPEC += (FTMP * lamBin );
    }

    MWXT_FRAC   = SEDMODEL_TABLE_MWXT_FRAC[IFILT][ispec] ;
    x0fac       = (x0 *MWXT_FRAC) ;
    FLUXGEN_forSPEC  *= (x0fac*SEDNORM_forSPEC) ;
    FLUXGEN_forMAG   *= (x0fac*SEDNORM_forMAG ) ; 

    MAG  = MAG_UNDEFINED ;
    if ( ZP > 0.0 && FLUXGEN_forMAG > 1.0E-50 ) 
      { MAG = ZP - 2.5*log10(FLUXGEN_forMAG) ; }

    // load function output 
    GENFLUX_LIST[ispec] = FLUXGEN_forSPEC ;
    GENMAG_LIST[ispec]  = MAG;

    // --------- check DUMP ---------

    NDMP_SKIP++ ;
    if ( LDMP && FLUXGEN_forSPEC > 0.0 && NDMP_SKIP==10 ) {
      printf(" xxx ispec=%3d  LAMAVG=%8.1f  ZP=%.2f "
	     " GEN[FLUX,MAG] = %.3le , %.3f \n",
	     ispec, LAMAVG, ZP, FLUXGEN_forSPEC, MAG );
      NDMP_SKIP=0 ;
    }

  } // end ispec

  // ------------------------------------------------
  // ------- apply host galaxy extinction -----------

  double MAGOFF_XT, FRAC, FRAC_XT, LAM ;
  int  DOXT = ( AV_host > 1.0E-9 ) ;
  // check before removing: int NBSPEC = SPECTROGRAPH_SEDMODEL.NBLAM_TOT ;

  if ( MAGOFF == 0.0  &&  DOXT==0 ) { return ; }

  FRAC = pow(TEN, -0.4*MAGOFF);
  MAGOFF_XT=0.0 ; FRAC_XT=1.0 ;

  for(ispec=0; ispec < NBSPEC; ispec++ ) {
  
    // get extinction from host in rest-frame (Jan 15, 2012)
    if ( DOXT ) {
      LAM       = SPECTROGRAPH_SEDMODEL.LAMAVG_LIST[ispec] ;
      MAGOFF_XT = GALextinct ( RV_host, AV_host, LAM, 94 ); 
      FRAC_XT   = pow(TEN, -0.4*MAGOFF_XT);
    }

    GENMAG_LIST[ispec]  += ( MAGOFF + MAGOFF_XT ) ;
    GENFLUX_LIST[ispec] *= ( FRAC * FRAC_XT );
  }


  // - - - - - - - - - - - - - - - -  
  if ( LDMP ) { debugexit(fnam); }

  return ;

}  // getSpec_SEDMODEL
