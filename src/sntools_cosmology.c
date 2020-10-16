
/*********************************************************
**********************************************************

      Cosmology and SFR Utilities 
      [moved from sntools.c, Oct 2020]

**********************************************************
**********************************************************/

#include "sntools.h"
#include "sntools_cosmology.h"

// ***********************************
void init_HzFUN_INFO(double *cosPar, char *fileName, HzFUN_INFO_DEF *HzFUN) {

  // initialize input struct HzFUN with either
  //  * wCDM params in cosPar (analytic cosmology theory)
  //  * Hz vs. z map read from fileName (used for interpolation)

  int ipar;
  char fnam[] = "init_HzFUN_INFO";

  // ----------- BEGIN ------------

  print_banner(fnam);

  for(ipar=0; ipar < NCOSPAR_HzFUN; ipar++ ) 
    { HzFUN->COSPAR_LIST[ipar] = -9.0; }

  HzFUN->Nzbin_MAP = 0;

  // - - - - - - 
  HzFUN->USE_MAP = !IGNOREFILE(fileName) ;
  if ( HzFUN->USE_MAP ) {

    
  }
  else {
    for(ipar=0; ipar < NCOSPAR_HzFUN; ipar++ ) 
      { HzFUN->COSPAR_LIST[ipar] = cosPar[ipar]; }
    
    double H0 = HzFUN->COSPAR_LIST[ICOSPAR_HzFUN_H0] ;
    double OM = HzFUN->COSPAR_LIST[ICOSPAR_HzFUN_OM] ;
    double OL = HzFUN->COSPAR_LIST[ICOSPAR_HzFUN_OL] ;
    double w0 = HzFUN->COSPAR_LIST[ICOSPAR_HzFUN_w0] ;
    double wa = HzFUN->COSPAR_LIST[ICOSPAR_HzFUN_wa] ;
    double Ok = 1.0 - OM - OL ;
    printf("\t H0         = %.2f      # km/s/Mpc \n",  H0);
    printf("\t OM, OL, Ok = %7.5f, %7.5f, %7.5f \n", OM, OL, Ok );
    printf("\t w0, wa     = %6.3f, %6.3f \n", w0, wa);
    fflush(stdout) ;
  }

  return ;

} // end init_HzFUN_INFO


// ****************************
double SFR_integral( double H0, double OM, double OL, double W, double Z ) {

  /***
   Integrate SFR(t) from 0 to current time.
   For convenience, integrate  over 'a' instead of redshift
   [since z goes to infinity]
 
        /a
       |   SFR(a') 
     c |  ----------- da'
       |  a' * H(a')
       /0

  ***/

  int    ia, NABIN = 100 ;
  double AMIN, AMAX, ABIN, atmp, ztmp, xa, tmp  ;
  double sum, sfr, aH ;

  double SECONDS_PER_YEAR = 3600. * 24. * 365. ;

  // ---------- BEGIN ------------

  AMIN = 0.0 ;
  AMAX = 1. / (1. + Z) ;
  ABIN = (AMAX - AMIN) / (double)NABIN ;

  sum = 0.0 ;

  for ( ia=1; ia <= NABIN; ia++ ) {
    xa   = (double)ia ;
    atmp = AMIN + ABIN * ( xa - 0.5 ) ;
    ztmp = (1. / atmp) - 1.0 ;

    sfr = SFRfun_BG03(H0,ztmp);
    aH  = atmp * Hzfun(H0,OM,OL,W,ztmp) ;
    sum += sfr / aH ;
  }

  // convert H (km/s/Mpc) to H(/year)

  tmp = (1.0E6 * PC_km) / SECONDS_PER_YEAR ;
  sum *= (ABIN * tmp) ;
  return sum ;

}  // end of function SFR_integral

// ****************************
double SFRfun_BG03(double H0, double z) {

  /***
      Compute sfr(z) = (a+b*z)/(1+(z/c)^d)*h solar_mass/yr/Mpc^3
      where H0 unit is km/s/Mpc and  h = H0/(100 km/s/Mpc)
     
      Oct 16 2020: rename SFRfun -> SFRfun_BG03
  ***/

  // -- Baldry and Glazebrook IMF --
  // https://ui.adsabs.harvard.edu/abs/2003ApJ...593..258B

  double a = 0.0118 ;
  double b = 0.08 ;
  double c = 3.3 ;
  double d = 5.2 ;
  double tmp1, tmp2, zc, h, SFRLOC;

  // ------------ BEGIN --------------
  zc   = z/c ;
  tmp1 = a + b*z ;
  tmp2 = 1.0 + pow(zc,d) ;
  h    = H0 / 100. ; 
  SFRLOC = h * tmp1 / tmp2 ;
  return(SFRLOC) ;
}  // end SFRfun_BG03

// *******************************************
double SFRfun_MD14(double z, double *params) {

  // Created Dec 2016 by R.Kessler
  // use function from  Madau & Dickoson 2014,
  // that was also used in Strolger 2015 for CC rate.
  // This function intended for CC rate, so beware
  // using for other purposes (e.g., no H0 factor here).

  double A      = params[0];
  double B      = params[1];
  double C      = params[2];
  double D      = params[3];
  double z1     = 1.0 + z;
  double top    = A*pow(z1,C);
  double bottom = 1.0 + pow( (z1/B), D );
  double SFR    = top/bottom; // also see Eq 8+9 of Strolger 2015 
  return(SFR);

} // end SFRfun_MD14

// *******************************************
double dVdz_integral
( 
  double H0    // (I) km/s per MPc
  ,double OM    // (I) Omega_matter
  ,double OL    // (I) Omega_lamba
  ,double W     //  (I) w = rho/p
  ,double Zmax  //  (I) integrate up to this redshift
  ,int wgtopt   //  (I) weight integral by z^wgtopt
 ) {

  //
  // return integral of dV/dz = r(z)^2/H(z) dz
  // wgtopt = 0:  returns standard volume integral
  // wgtopt = 1:  returns  z-wgted integral

  double sum, tmp, dz, Ztmp, wz, xz ;
  int NZbin, iz;

  // ---- BEGIN ----------

  // compute exact integral

  NZbin = (int)( Zmax * 1000.0 ) ;
  if ( NZbin < 10 ) { NZbin = 10 ; }
  dz   = Zmax / (float)NZbin ;   // integration binsize
  sum = 0.0;

  for ( iz=1; iz <= NZbin; iz++ ) {
    xz   = (double)iz ;
    Ztmp = dz * (xz - 0.5) ;
    tmp  = dVdz ( H0, OM, OL, W, Ztmp );

    wz = pow(Ztmp,(double)wgtopt);
    sum += wz * tmp;

  }

  sum *= dz ;

  //  printf(" xxxx dVdz_integral = %e (approx=%e) \n", sum, sumtmp  );

  return sum ;

}  // end of dVdz_integral


double dvdz_integral__(double *H0, double *OM, double *OL, double *W,
		       double *Zmax, int *wgtopt) {
  return dVdz_integral(*H0,*OM,*OL,*W,*Zmax,*wgtopt);
}



// **********************************
double dVdz 
( 
  double H0    // (I) km/s per MPc
 ,double OM    // (I) Omega_matter
 ,double OL    // (I) Omega_lamba
 ,double W    //  (I) w = rho/p
 ,double Z    //  (I) redshift
 ) {

  // returns dV/dz = r(z)^2 / H(z)

  double r, H, tmp ;
  double zero = 0.0;

  r = Hzinv_integral ( H0, OM, OL, W, zero, Z );
  H = Hzfun ( H0, OM, OL, W, Z );

  tmp = LIGHT_km * r * r / H ;

  return tmp;

}  // end of dVdz


// ******************************************
double Hzinv_integral 
( 
 double H0      // (I) km/s per MPc
 ,double OM     // (I) Omega_matter
 ,double OL     // (I) Omega_lamba
 ,double W      //  (I) w = rho/p
 ,double Zmin   //  (I) min integ bin
 ,double Zmax   //  (I) integrate up to this redshift
 ) {

  // 
  // Jun 2016: bug fix, (float)NZbin -> (double)NZbin

  int iz, NZbin ;
  double dz, Hz, xz, Ztmp, sum, Hzinv, KAPPA, SQRT_KAPPA ; 

  // ------ return integral c*r(z) = int c*dz/H(z) -------------
  // Note that D_L = (1+z)*Hzinv_integral

  sum = 0.0;

  NZbin = (int)( (Zmax-Zmin) * 1000.0 ) ;
  if ( NZbin < 10 ) { NZbin = 10 ; }
  dz  = (Zmax-Zmin) / (double)NZbin ;      // integration binsize

  for ( iz=1; iz <= NZbin; iz++ ) {
    xz   = (double)iz ;
    Ztmp = Zmin + dz * (xz - 0.5) ;
    Hz   = Hzfun ( H0, OM, OL, W, Ztmp );
    sum += (1.0/Hz) ;
  }

  // remove H0 factor from inetgral before checking curvature.

  sum *= (dz * H0) ;

  // check for curvature
  KAPPA      = 1.0 - OM - OL ; 
  SQRT_KAPPA = sqrt(fabs(KAPPA));

  if ( KAPPA < -0.00001 ) 
    { Hzinv = sin( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else if ( KAPPA > 0.00001 ) 
    { Hzinv = sinh( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else
    { Hzinv = sum ; }


  // return Hzinv with c/H0 factor
  return (Hzinv * LIGHT_km / H0 ) ;

} // end of Hzinv_integral


// ******************************************
double Hainv_integral 
( 
 double H0     // (I) km/s per MPc
 ,double OM     // (I) Omega_matter
 ,double OL     // (I) Omega_lamba
 ,double W      //  (I) w = rho/p
 ,double amin   //  (I) min integ bin
 ,double amax   //  (I) integrate up to this redshift
 ) {

  // May 29, 2008: same as Hzinv_integral, but integrate over a
  // instead of over z.
  // dz/E(z) :  z=1/a-1   dz = -da/a^2

  int ia, Nabin ;
  double da, Hz, xa, atmp, Ztmp, sum, Hzinv, KAPPA, SQRT_KAPPA ; 

  // ------ return integral c*r(z) = int c*dz/H(z) -------------
  // Note that D_L = (1+z)*Hzinv_integral

  sum = 0.0;

  Nabin = (int)( (amax-amin) * 1000.0 ) ;
  if ( Nabin < 10 ) { Nabin = 10 ; }
  da   = (amax-amin) / (double)Nabin ;   // integration binsize

  for ( ia=1; ia <= Nabin; ia++ ) {
    xa   = (double)ia ;
    atmp = amin + da * (xa - 0.5) ;
    Ztmp = 1./atmp - 1.0 ;
    Hz   = Hzfun ( H0, OM, OL, W, Ztmp );
    sum += 1.0/( Hz * atmp * atmp) ;
  }

  // remove H0 factor from inetgral before checking curvature.

  sum *= (da * H0) ;

  // check for curvature
  KAPPA      = 1.0 - OM - OL ; 
  SQRT_KAPPA = sqrt(fabs(KAPPA));

  if ( KAPPA < -0.00001 ) 
    { Hzinv = sin( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else if ( KAPPA > 0.00001 ) 
    { Hzinv = sinh( SQRT_KAPPA * sum ) / SQRT_KAPPA ; }
  else
    { Hzinv = sum ; }


  // return Hzinv with c/H0 factor
  return (Hzinv * LIGHT_km / H0 ) ;

} // end of Hainv_integral


// ******************************************
double Hzfun ( double H0, double OM, double OL, double W, double Z ) {

  //  int iz, NZbin ;
  double sqHz, Hz, ZZ, Z2, Z3, ZL, WW, KAPPA ;

  // ------ returns H(z) -------------

  KAPPA = 1.0 - OM - OL ;  // curvature

  ZZ  = 1.0 + Z ;
  Z2  = ZZ * ZZ ;  // avoid pow fun
  Z3  = Z2 * ZZ ; 

  WW  = 3.0 * (1.0 + W) ;
  ZL  = pow(ZZ,WW) ;

  sqHz = OM*Z3  +  OL*ZL  + KAPPA*Z2 ;

  Hz = H0 * sqrt ( sqHz ) ;

  return Hz ;

} // end of Hzfun


// ******************************************
double dLmag ( double H0, double OM, double OL, double W, 
	       double zCMB, double zHEL ) {

  // returns luminosity distance in mags:
  //   dLmag = 5 * log10(DL/10pc)
  //
  // BEWARE that input H0 is 1/seconds (not km/s/Mpc!)
  // Jan 5 2015: pass zCMB and zHEL
  //
  double rz, dl, arg, mu ;
  double zero = 0.0 ;
  // ----------- BEGIN -----------
  rz     = Hzinv_integral ( H0, OM, OL, W, zero, zCMB );
  dl     = ( 1.0 + zHEL ) * rz ;
  arg    = (double)10.0 * PC_km / dl ;
  mu     = -5.0 * log10( arg );
  return mu ;
}  // end of dLmag

/* xxxx mark delete xxxxx
double dlmag_ (double *H0, double *OM, double *OL, double *W,
               double *zCMB, double *zHEL ) {
  double mu = dLmag(*H0,*OM,*OL,*W,*zCMB,*zHEL);
  return(mu);
}
xxxxxxxx */

double zcmb_dLmag_invert( double H0, double OM, double OL, double W, double MU ) {

  // Created Jan 4 2018
  // for input distance modulus (MU), solve for zCMB.
  // Beware that H0 unit is km/s/pc (not per Mpc)

  double zCMB, zCMB_start, dmu, DMU, mutmp, DL ;
  double DMU_CONVERGE = 1.0E-4 ;
  int    NITER=0;
  char fnam[] = "zcmb_dLmag_invert" ;

  // ---------- BEGIN ----------

  // use naive Hubble law to estimate zCMB_start
  DL    = pow( 10.0,(MU/5.0) ) * 1.0E-5 ; // Mpc
  zCMB_start = (70.0*DL)/LIGHT_km ;
  zCMB_start *= exp(-zCMB_start/6.0); // very ad-hoc estimate

  zCMB = zCMB_start ;
  DMU = 9999.0 ;
  while ( DMU > DMU_CONVERGE ) {
    mutmp  = dLmag(H0,OM,OL,W, zCMB, zCMB); // MU for trial zCMB
    dmu    = mutmp - MU ;             // error on mu
    DMU    = fabs(dmu);
    zCMB  *= exp(-dmu/2.0); 

    NITER++ ;
    if ( NITER > 500 ) {
      sprintf(c1err,"Could not solve for zCMB after NITER=%d", NITER);
      sprintf(c2err,"MU=%f  dmu=%f  ztmp=%f", MU, dmu, zCMB);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  } // end dz                                                                   

  int LDMP=0 ;
  if ( LDMP ) {
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx MU=%.4f -> DL = %.2f Mpc  zCMB(start) = %.5f \n",
	   MU, DL, zCMB_start );
    printf(" xxx zCMB -> %.5f  DMU=%f after %d iterations \n",
	   zCMB, DMU, NITER);
  }

  return(zCMB);

} // end zcmb_dLmag_invert



// ************************************************
double zhelio_zcmb_translator (double z_input, double RA, double DEC, 
			       char *coordSys, int OPT ) {

  /**********************

   Created Dec 2013 by R.Kessler
   General redshift-translator function to between helio and cmb frame.
   [replaces Z2CMB that translated only in 1 direction]

   OPT > 0  -> z_input = z_helio, return z_out = zcmb
   OPT < 0  -> z_input = z_cmb,   return z_out = zhelio
   RA,DEC   = sky coordinates
   coordSys = coordinate system; e.g., 'eq' or 'gal' or 'J2000'


   l = longitude = RA (deg)
   b = lattitue  = DEC  (deg)

   Use exact formuala,
   
    1 + z_cmb = ( 1 + z_helio ) / ( 1 - V0 . nhat/c )

   where V0 is the CMB velocity vector, 
   nhat is the unit vector between us and the SN,
   and c = 3E5 km/s.

   Note that the NED-calculator used in JRK07 is 
   an approximation, z_cmb = z_helio + V0.nhat/c,
   that is OK when z_helio << 1.

 ****************/

  double 
     ra_gal, dec_gal
    ,ss, ccc, c1, c2, c3, vdotn, z_out
    ;

  char fnam[] = "zhelio_zcmb_translator" ;

  // --------------- BEGIN ------------

  // on negative redshift, just return input redshift with
  // no calculation. Allows flags such as -9 to be unperturbed.
  if ( z_input < 1.0E-10 ) { return z_input ; }

  if ( strcmp(coordSys,"eq"   ) == 0 || 
       strcmp(coordSys,"J2000") == 0 ) {

    // input and output in degrees
    slaEqgal( RA, DEC, &ra_gal, &dec_gal ) ;
  }
  else if ( strcmp(coordSys,"gal") == 0 ) {
    ra_gal  = RA ;
    dec_gal = DEC ;
  }
  else {
    sprintf(c1err,"Invalid coordSys = '%s' ", coordSys );
    sprintf(c2err,"OPT=%d z_in=%f RA=%f DEC=%f", OPT, z_input, RA, DEC);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // get projection

  ss  = sin(RADIAN*dec_gal) * sin(RADIAN*CMBapex_b);
  c1  = cos(RADIAN*dec_gal) ;
  c2  = cos(RADIAN*CMBapex_b) ;
  c3  = cos(RADIAN*(ra_gal-CMBapex_l));
  ccc = c1 * c2 * c3 ;
  vdotn = CMBapex_v * ( ss + ccc ) / LIGHT_km ;

  z_out = -9.0 ;

  if ( OPT > 0 ) {
    z_out  = ( 1. + z_input ) / ( 1. - vdotn ) - 1. ; 
  }
  else if ( OPT < 0 )  {
    z_out  = ( 1. + z_input) * ( 1. - vdotn ) - 1.0 ;  
  }
  else if ( OPT == 0 ) {
    sprintf(c1err,"Invalid OPT=0" );
    sprintf(c2err,"z_input=%f  RA=%f  DEC=%f", z_input,RA,DEC);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
   
  return(z_out) ;


} // end of zhelio_zcmb_translator 


double zhelio_zcmb_translator__ (double *z_input, double *RA, double *DEC,
                                 char *coordSys, int *OPT ) {
  return zhelio_zcmb_translator(*z_input, *RA, *DEC, coordSys, *OPT) ;
}



