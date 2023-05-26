/**************************************** 
  Created May 2023 by R.Kessler

  Tools to simulate atmosphere effects such as DCR(coords) and 
  DCR(PSF-shape) ; 
  Motivated by Le et al, 2023: https://arxiv.org/abs/2304.01858

  These functions lean heavily on INPUTS and GENLC structures
  in snlc_sim.h, so these tools are not easily plugged into
  other codes.

 ****************************************/

#include "sntools.h"
#include "sntools_cosmology.h"  // clumsy to require this include
#include "sntools_stronglens.h" // clumsy to ...

#include "snlc_sim.h"
#include "sntools_sim_atmosphere.h"

// ***************************************
void GEN_ATMOSPHERE_DRIVER(void) {

  // Created May 2023 by R.Kessler
  // Driver routine to simulate atmoshperic effects

  int NEPOCH = GENLC.NEPOCH ;
  int ep;
  char fnam[] = "GEN_ATMOSPHERE_DRIVER" ;

  // ------------ BEGIN ----------

  if ( INPUTS.ATMOSPHERE_OPTMASK == 0 ) { return; }

  for ( ep = 1; ep <= GENLC.NEPOCH; ep++ )  {  
    gen_airmass(ep);
    gen_dcr_coordShift(ep);
    genSmear_coords(ep);
  }
  
  return;

} // end GEN_ATMOSPHERE_DRIVER

 
// ==========================================
void gen_airmass(int epoch) {

  // Created May 2023 by R.Kessler and Jason Lee (U.Penn)
  // Compute airmass for this epoch
  // 

  double MJD      = GENLC.MJD[epoch];
  double RA       = GENLC.RA ;  // true RA
  double DEC      = GENLC.DEC ; // true DEC
  int    IDSURVEY = GENLC.IDSURVEY;

  double geoLAT, geoLON;
  get_geoSURVEY(IDSURVEY, &geoLAT, &geoLON);

  double RAD = RADIAN ;
  double airmass  = 1.11 ;
  double sin_alt, ang_zenith_rad, ang_zenith_deg, h_hr, h_deg, COS_h ;
  double GLAT, GLON;
    
  // test example from ESO calculator:
  // https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
  bool DO_TEST = false ; 
  double test_geoLAT  = -29.257 ;    // La Silla
  double test_geoLON  = -70.738 ;
  double test_MJD    = 59583.2409 ; // from LSST minion simlib
  double test_RA     = 149.;
  double test_DEC    = 2.2 ;

  char fnam[] = "gen_airmass" ;

  // ------------ BEGIN ------------

  GENLC.AIRMASS[epoch] = -9.0 ;


  if ( DO_TEST ) {
    MJD       = test_MJD ;
    geoLAT    = test_geoLAT ;
    geoLON    = test_geoLON ;
    RA  = GENLC.RA  = test_RA;
    DEC = GENLC.DEC = test_DEC ;
    slaEqgal(RA, DEC, &GLON, &GLAT ); // return GLON, GLAT in degrees
    GENLC.SIN_GLON = sin(GLON*RAD);
    GENLC.COS_GLON = cos(GLON*RAD);
    GENLC.SIN_DEC = sin(DEC*RAD);
    GENLC.COS_DEC = cos(DEC*RAD);

    SURVEY_INFO.sin_geoLAT[IDSURVEY] = sin(geoLAT*RAD) ;
    SURVEY_INFO.cos_geoLAT[IDSURVEY] = cos(geoLAT*RAD) ;

    printf("\n xxx %s: prep comparison with ESO calculator: \n", fnam);
    printf("\t xxx geo(LAT,LON) = %f , %f \n", geoLAT, geoLON);
    printf("\t xxx RA, DEC = %f , %f \n", RA, DEC);
    printf("\t xxx MJD = %f \n", MJD);
    fflush(stdout);
  }

  // if geo coords are not available, return
  if ( geoLAT > 1000.0 || geoLON > 1000.0 ) { return; }    

  double JD2000   = 2451545.0;
  int    iMJD     = (int)MJD;       // MJD of previous midnight
  double JD       = MJD + 2400000.5 ;
  double JD0      = (double)iMJD + 2400000.5 ;
  double D_UT     = JD0 - JD2000 ;
  double T_U      = (JD - JD2000)/36525.0 ; //Number of Julian Centuries since J2000.0

  // compute h = hourAngle = Local Siderial Time (LST) - RA
  // xxx double GMST_deg = fmod(18.697375 + 24.065709824279*D_UT, 24.0) * 360.0/24.0;
  double GMST_deg = (fmod(24110.54841 + 8640184.812866 * T_U + 0.093104 * T_U * T_U, 86400.0)/86400.0 + 1.0027379 * fmod(MJD, 1.0) ) * 360.0 ; // Grabbed this from http://www.astro.sunysb.edu/metchev/AST443/times.html#:~:text=GAST%20(Greenwich%20Apparent%20Sidereal%20Time,%2B0.8%20to%20%2B1.2%20seconds - GMST equation + adding number of hours from 0h UT.
  double LST_deg  = (geoLON + GMST_deg);
  h_deg           = LST_deg - GENLC.RA ;
  h_hr            = h_deg * 24.0/360.0 ;
  COS_h           = cos(h_deg*RAD) ;

  // compute airmass ...
  double SIN_geoLAT = SURVEY_INFO.sin_geoLAT[IDSURVEY];
  double COS_geoLAT = SURVEY_INFO.cos_geoLAT[IDSURVEY];
  double SIN_LON    = GENLC.SIN_GLON ;
  double COS_LON    = GENLC.COS_GLON ;

  // xx  double SIN_DEC    = sin(DEC*RAD) ; // GENLC.SIN_DEC ;
  // xx  double COS_DEC    = cos(DEC*RAD) ; // GENLC.COS_DEC ; 

  // avoid re-computing trig functions for each obs
  double SIN_DEC    = GENLC.SIN_DEC ;
  double COS_DEC    = GENLC.COS_DEC ; 

  sin_alt = 
    (SIN_geoLAT * SIN_DEC) + 
    (COS_geoLAT * COS_DEC * COS_h);

  ang_zenith_rad = 0.25*TWOPI - asin(sin_alt) ; // zenight angle, radians
  ang_zenith_deg = ang_zenith_rad / RAD ;

  airmass = 1.0/cos(ang_zenith_rad) ;

  GENLC.AIRMASS[epoch]    = airmass;
  GENLC.ANG_ZENITH[epoch] = ang_zenith_deg ;
  if ( DO_TEST ) {

    printf("\n xxx quantities computed by sim function %s: \n", fnam);
    printf("\t xxx GLON, GLAT = %f, %f \n", GLON, GLAT);
    printf("\t xxx GMST, LST = %f , %f deg \n", GMST_deg, LST_deg);
    printf("\t xxx hour angle h = %f deg = %f hr\n", h_deg, h_hr);
    printf("\t xxx ang_zenith = %f deg / %f rad  (RADIAN=%f)\n", 
	   ang_zenith_deg, ang_zenith_rad, RAD );
    printf("\t xxx airmass = %f \n", airmass);

    printf("\n xxx ESO calculator results:\n");
    printf("\t xxx Galactic coords = 235°.94 ,  41°.21 \n");
    printf("\t xxx Hour Angle HA = 22:03:14\n");
    printf("\t xxx Target az =  46°.66  alt =  47°.93 \n");
    printf("\t xxx Zenith distance =  42°.07 \n");
    printf("\t xxx Airmass =  1.347 \n");

    fflush(stdout);
    debugexit(fnam);
  }

  return;
} // end gen_airmass

// ==================================
void genSmear_coords(int epoch) {

  // Created May 2023: 
  // determined measured RA,DEC for this epoch
  // Start with silly model for testing ...

  double SNR_REF          = 100.0;
  double ANGRES_REF_mASEC = 0.010; // 10 milli-asec res for SNR_REF
  double ANGRES_REF_DEG   = ANGRES_REF_mASEC/3600.0;
  double cosDEC  = GENLC.cosDEC;
  double trueSNR = GENLC.trueSNR[epoch] ;
  if ( trueSNR < 0.01 ) { trueSNR = 0.01; }

  double ANGRES, ran_RA, ran_DEC, WGT ;
  char fnam[] = "genSmear_coords" ;

  // ----------- BEGIN ---------

  ran_RA  = getRan_Gauss(1);
  ran_DEC = getRan_Gauss(1);

  ANGRES = ANGRES_REF_DEG * sqrt(SNR_REF/trueSNR);

  GENLC.RA_OBS[epoch]  = GENLC.RA  + (ANGRES * ran_RA)/cosDEC;
  GENLC.DEC_OBS[epoch] = GENLC.DEC + (ANGRES * ran_DEC) ;

  // ?? what about uncertainty ???

  // update wgted-avg among all epochs
  WGT = (ANGRES_REF_DEG*ANGRES_REF_DEG) / (ANGRES*ANGRES);
  GENLC.RA_SUM     += WGT * GENLC.RA_OBS[epoch] ;
  GENLC.DEC_SUM    += WGT * GENLC.DEC_OBS[epoch] ;
  GENLC.RA_WGTSUM  += WGT ;
  GENLC.DEC_WGTSUM += WGT ;

  GENLC.RA_AVG      = GENLC.RA_SUM  / GENLC.RA_WGTSUM ;
  GENLC.DEC_AVG     = GENLC.DEC_SUM / GENLC.DEC_WGTSUM ;

  //  printf(" xxx %s: DEC_AVG -> %f for CID=%d \n",
  //	 fnam, GENLC.DEC_AVG, GENLC.CID);

  return;

} // end genSmear_coords


// ========================================
void gen_dcr_coordShift(int ep) {

  // Created May 2023
  // Compute DCR astrometric shift for RA and DEC.

  double AIRMASS    = GENLC.AIRMASS[ep];
  double ANG_ZENITH = GENLC.ANG_ZENITH[ep];
  
  double RA       = GENLC.RA ;  // true RA 
  double DEC      = GENLC.DEC ; // true DEC   
  double LON      = GENLC.GLON;
  double LAT      = GENLC.GLAT ;
  int    IDSURVEY = GENLC.IDSURVEY;
  double geoLAT, geoLON;
  get_geoSURVEY(IDSURVEY, &geoLAT, &geoLON);

  // .xyz SN-SED top of atmosphere (not thru filter),
  //    PSF ... effective <lambda> for SED !!!!
  //  integral[ l*SED(l)*Filt(l) / integral[SED(l)*Filt(l)]

  double wave_sed_wgted;
  char fnam[] = "gen_dcr_coordShift" ;

  // -------------- BEGIN -----------

  wave_sed_wgted = gen_wave_sed_wgted(ep);

  GENLC.RA_dcr_shift[ep]  = 0.0 ; // true shift; degrees
  GENLC.DEC_dcr_shift[ep] = 0.0 ;

  return ;

} // end gen_dcr_coordShift

// ================================
double gen_wave_sed_wgted(int ep) {

  // Created May 2023
  // Return mean wavelegnth in filter corresponding to epoch 'ep';
  // use SED * FILTER_TRANS weight.
  double wave = 0.0 ;

  char fnam[] = "gen_wave_sed_wgted";

  // --------- BEGIN ---------

  return wave;

} // end gen_wave_sed_wgted


