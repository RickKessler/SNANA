/******************************************

 Created Sep 2023 by R.Kessler

 Utilities to extrapolate LC mags to phases beyond nominal phase range of model.
 This is useful for simulations so that there are no sharp cutoffs in the
 light curve.
 These extrap models are designed to be fairly generic so that they can be 
 applied to multiple models.

     UNDER CONSTRUCTION: NOT YET IMPLEMENTED 

*****************************************/

#include "sntools.h"
#include "genmag_extrap.h"


// ***************************************************
double modelflux_extrap(
			double T          // (I) Time to extrapolate to
			,double Tref      // (I) nearest time with model flux
			,double fluxref   // (I) flux at Tref
			,double fluxslope // (I) dflux/dT  at Tref
			,int LDMP         // (I) print to screen if true
			) {

  /****
    Created Apr 12, 2009 by R.Kessler
       [moved from sntools.c[h] to here, Sep 2023]

    Extrapolates model flux from Tref to time T
    (rest-frame), where T is outside range of model.
    T and Tref are in days, defined so that T=0 at peak brightness.

    T < Tref => 
    Rises like  F0*(T-T0)^2 x (T-T1) where the quadratic term 
    is a simple fireball model, and the 2nd (linear) term allows 
    matching both the function and its derivative at T=Tref.

    T > Tref => on the tail; use linear mag-extrap, which
                is a power-law flux-extrap.  If fluxref < 0,
                then just hard-wire flux=1E-30 since model
                is unphysical and extrapolation makes no sense.
  
  ****/

  double flux, slope, arg, F0, FP, T1 ;
  double DTref0, SQDTref0, CUBEDTref0 ,DT0, DT1 ;
  double T0  = -20.0 ; // time of explosion, wher T=0 at peak

  char fnam[] = "modelflux_extrap" ;

  // -------- BEGIN ----------

  // always use same sign to ensure that flux-> 0 at T->infinity
  slope = fabs(fluxslope) ;

  if ( T > Tref ) {

    // linear mag-extrap at late times => power-law extrap for flux
    if ( fluxref > 0.0 ) {
      arg  = -(T-Tref)*slope / (LNTEN * fluxref) ;
      flux = fluxref * pow(TEN,arg) ;
    }
    else { 
      flux = 1.0E-30 ; 
    }
  }
  else {

    if ( T < T0 ) { flux = 0.0 ; return flux ; }

    // solve for F0 and T1 to describe pre-max rise
    FP         = fluxslope ;
    DTref0     = Tref - T0;
    SQDTref0   = DTref0*DTref0;
    CUBEDTref0 = SQDTref0*DTref0;

    F0  = (FP * DTref0 - 2.*fluxref)/CUBEDTref0 ;
    T1  = Tref - fluxref / (F0*SQDTref0) ;

    DT0 = T - T0 ;
    DT1 = T - T1 ;
    flux = F0 * DT0*DT0 * DT1 ;
  }

  if ( LDMP > 0 ) {
    printf("FLUX_EXTRAP: F(T=%5.1f)=%10.3le", T, flux );
    printf(" [Fref(Tref=%3.0f)=%10.3le & dF/dT=%10.3le] \n",
	   Tref, fluxref, fluxslope );
  }

  return flux ;

} // end of modelflux_extrap


// ==========================================================
void init_extrap_latetime_Ia(void) {

  // Created June 25 2018
  // Init optional mag-extrapolation for epochs later than
  // what is defined in SALT2 model. 
  // Note that default extrapolation is in SED-flux space,
  // extrapolation last few days, but this extrap can have
  // large errors.
  //
  // Input: model_extrap_latetime --> fileName with model info

  char *fileName = INPUT_EXTRAP_LATETIME_Ia.FILENAME ;
  char fnam[] = "init_extrap_latetime_Ia" ;

  int    ipar, ilam, NLAMBIN=0;
  int    NPAR_READ = NPAR_EXTRAP_LATETIME_Ia ;
  double DAYMIN, LAM, TAU1, TAU2, EXPRATIO, MAGSLOPE1, MAGSLOPE2, DAYPIVOT;
  double TMPVAL[10];

  FILE *fp;
  char c_get[60];

  // -------------- BEGIN ----------

  // xxxx  printf(" xxx %s: filename = %s \n", fnam, fileName); fflush(stdout);

  INPUT_EXTRAP_LATETIME_Ia.NLAMBIN = 0 ;

  if ( IGNOREFILE(fileName) ) { return ; }
  ENVreplace(fileName,fnam,1);

  fp = fopen(fileName,"rt");
  if ( !fp ) {
    sprintf(c1err,"Could not open MODEL_EXTRAP_LATETIME:");
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  printf("\n   Read EXTRAP_LATETIME(SNIa) parameters from :\n");
  printf("\t %s \n", fileName);
  fflush(stdout);

  INPUT_EXTRAP_LATETIME_Ia.NLAMBIN   = 1 ;
  INPUT_EXTRAP_LATETIME_Ia.DAYMIN    = 0.0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,"EXTRAP_DAYMIN:") == 0 ) 
      { readdouble(fp, 1, &INPUT_EXTRAP_LATETIME_Ia.DAYMIN ); }

    if ( strcmp(c_get,"EXTRAP_PARLIST:") == 0 ) { 
      readdouble(fp, NPAR_READ, TMPVAL );
      if ( NLAMBIN < MXLAMBIN_EXTRAP_LATETIME_Ia ) {
	for(ipar=0; ipar < NPAR_READ; ipar++ ) 
	  { INPUT_EXTRAP_LATETIME_Ia.PARLIST[ipar][NLAMBIN] = TMPVAL[ipar]; } 
      }
      NLAMBIN++ ;
      INPUT_EXTRAP_LATETIME_Ia.NLAMBIN = NLAMBIN ;
    }
  }

  fclose(fp);


  if ( NLAMBIN >= MXLAMBIN_EXTRAP_LATETIME_Ia ) {
    sprintf(c1err,"NLAMBIN=%d exceeds bound of %d", 
	    NLAMBIN, MXLAMBIN_EXTRAP_LATETIME_Ia);
    sprintf(c2err,"Check MXLAMBIN_EXTRAP_LATETIME_Ia = %d",
	    MXLAMBIN_EXTRAP_LATETIME_Ia);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // -----------------------------------------------------
  // prep/check stuff stuff


  NLAMBIN = INPUT_EXTRAP_LATETIME_Ia.NLAMBIN ;
  DAYMIN  = INPUT_EXTRAP_LATETIME_Ia.DAYMIN  ;

  if ( DAYMIN < 10.0 ) { 
    sprintf(c1err,"Invalid DAYMIN=%.2f (too small)", DAYMIN);
    sprintf(c2err,"Check EXTRAP_DAYMIN key for SNIa");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // compute a few quantities and print info for each lambin

  printf("\t DAYMIN_EXTRAP = %.1f \n", DAYMIN );
  printf("\n\t FLUX_EXTRAP(t) ~ [ exp(t/TAU1) + RATIO*exp(t/TAU2) ] \n");

  printf("                                TAU1     TAU2     DAY when\n");
  printf("   LAM    TAU1   TAU2   RATIO  mag/day  mag/day    F1=F2 \n");
  printf("   --------------------------------------------------------\n");

  for(ilam=0; ilam < NLAMBIN; ilam++ ) {  
    LAM  = INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_LAM_Ia][ilam] ;
    TAU1 = INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_TAU1_Ia][ilam] ;
    TAU2 = INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_TAU2_Ia][ilam] ;

    if ( TAU2 < TAU1 ) {
      sprintf(c1err,"Invalid TAU2(%.2f) < TAU1(%.2f)", TAU2, TAU1);
      sprintf(c2err,"Check EXTRAP_PARLIST with lam=%.1f", LAM);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    EXPRATIO = INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_EXPRATIO_Ia][ilam] ;

    MAGSLOPE1 = 1.086/TAU1;
    MAGSLOPE2 = 1.086/TAU2;
    if ( EXPRATIO > 1.0E-9 && TAU1 > 0.0 && TAU2 > 0.0 ) 
      { DAYPIVOT  = log(1.0/EXPRATIO) / (1.0/TAU1 - 1.0/TAU2); }
    else
      { DAYPIVOT  = 1.0E4; }

    INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_MAGSLOPE1_Ia][ilam] = MAGSLOPE1;
    INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_MAGSLOPE2_Ia][ilam] = MAGSLOPE2;
    INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_DAYPIVOT_Ia][ilam]  = DAYPIVOT;

    printf(" %7.1f %6.2f %6.2f %6.4f  %6.3f   %6.3f     %.0f \n",
	   LAM, TAU1, TAU2, EXPRATIO,    MAGSLOPE1, MAGSLOPE2, DAYPIVOT);
    fflush(stdout);
  }
  printf("   --------------------------------------------------------\n");

  return ;

} // end init_extrap_latetime_Ia

// ===============================================
double genmag_extrap_latetime_Ia(double mag_daymin, double day, double lam ) {

  // Created Jun 25 2018
  // for input mag_daymin, return extrapolated magnitude.
  //
  // Inputs:
  //   mag_daymin = mag (obs or rest) to extrapolate from 'day' 
  //   day        = rest-frame day (day=0 at peak brightness)
  //   lam        = rest-frame wavelength of filter
  //

  int    NLAMBIN = INPUT_EXTRAP_LATETIME_Ia.NLAMBIN ;  
  double DAYMIN  = INPUT_EXTRAP_LATETIME_Ia.DAYMIN ;  
  double mag_extrap = mag_daymin;

  double arg, F_DAYMIN, F_EXTRAP, VAL, PARLIST[MXPAR_EXTRAP_LATETIME_Ia];
  double *ptrLam, *ptrVal;
  int    ipar;
  int    NPAR = NPAR_EXTRAP_LATETIME_Ia ;
  int    OPT_INTERP = 1; // linear
  int    LDMP = 0, ABORT=0 ;
  char   fnam[] = "genmag_extrap_latetime_Ia" ;

  // ----------- BEGIN ---------

  if ( day < DAYMIN ) {
    sprintf(c1err,"Invalid day=%.2f is < DAYMIN=%.2f", day, DAYMIN);
    sprintf(c2err,"day must be > DAYMIN");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // compute flux at daymin
  arg      = 0.4*(mag_daymin - ZEROPOINT_FLUXCAL_DEFAULT);
  F_DAYMIN = pow(10.0,-arg);

  // interpolate each extrap parameter vs. wavelength
  for(ipar=1; ipar < NPAR; ipar++ ) { // skip LAM parameter
   
    ptrLam = INPUT_EXTRAP_LATETIME_Ia.PARLIST[IPAR_EXTRAP_LAM_Ia] ;
    ptrVal = INPUT_EXTRAP_LATETIME_Ia.PARLIST[ipar] ;
    if ( lam < ptrLam[0] ) {
      VAL = ptrVal[0];
    }
    else if ( lam > ptrLam[NLAMBIN-1] ) {
      VAL = ptrVal[NLAMBIN-1];
    }
    else {
      VAL = interp_1DFUN(OPT_INTERP, lam, 
			 NLAMBIN, ptrLam, ptrVal, fnam );
    }
    PARLIST[ipar] = VAL ;
  }

  // ----------
  
  double TAU1  = PARLIST[IPAR_EXTRAP_TAU1_Ia] ;
  double TAU2  = PARLIST[IPAR_EXTRAP_TAU2_Ia] ;
  double RATIO = PARLIST[IPAR_EXTRAP_EXPRATIO_Ia] ;
  double DAYDIF = 0.0 ;
  double FTMP, FNORM ;

  // get reference extrap flux at DAYDIF = DAY-DAYMIN =0
  FTMP  = FLUXFUN_EXTRAP_LATETIME_Ia(DAYDIF,TAU1,TAU2,RATIO);
  FNORM = F_DAYMIN / FTMP;

  DAYDIF   = day - DAYMIN ;
  F_EXTRAP = FNORM * FLUXFUN_EXTRAP_LATETIME_Ia(DAYDIF,TAU1,TAU2,RATIO);

  mag_extrap = ZEROPOINT_FLUXCAL_DEFAULT - 2.5*log10(F_EXTRAP);

  if ( mag_extrap > 40.0 ) { mag_extrap = MAG_ZEROFLUX; }

  if ( mag_extrap < 0.0 || mag_extrap > 99. ) { ABORT=1; }
  if ( LDMP || ABORT ) {
    printf(" xxx \n");
    printf(" xxx -------- DUMP   %s  ---------- \n", fnam);
    printf(" xxx INPUTS: mag_daymin=%.3f  day=%.3f  lam=%.1f \n",
	   mag_daymin, day, lam);
    printf(" xxx TAU1=%.3f  TAU2=%.3f  RATIO=%.5f \n", 
	   TAU1, TAU2, RATIO );
    printf(" xxx F_DAYMIN = %f   FLUXFUN_EXTRAP(0)=%f \n",
	   F_DAYMIN, FTMP);
    printf(" xxx DAYDIF=%.2f  F_EXTRAP=%f  --> mag_extrap=%.3f \n",
	   DAYDIF, F_EXTRAP, mag_extrap);

    if ( ABORT ) {
      sprintf(c1err,"Crazy mag_extrap = %le", mag_extrap);
      sprintf(c2err,"Check above DUMP");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
    }

  }

  return(mag_extrap);

} // end genmag_extrap_latetime_Ia

double FLUXFUN_EXTRAP_LATETIME_Ia(double t, double tau1, double tau2, 
				  double ratio) {
  double F1 = exp(-t/tau1);
  double F2 = ratio * exp(-t/tau2);
  double F  = F1 + F2;
  return(F);
} 
