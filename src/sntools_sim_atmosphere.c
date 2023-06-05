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
#include "genmag_SEDtools.h"

#include "snlc_sim.h"
#include "sntools_sim_atmosphere.h"
#include "sntools_spectrograph.h"
#include "sntools_trigger.h"

//#include "fitsio.h"
//#include "sntools_calib.h"

#define UNIT_TEST_AIRMASS      false 
#define UNIT_TEST_COMPUTE_DCR  false

// *****************************************
void INIT_ATMOSPHERE(void) {

  // Created Jun 2023
  // One-time init to prepare to simulate DCR effects on 
  // coordinates and PSF-fitted mags.

  int  ID = GENLC.IDSURVEY ;

  int ifilt, ifilt_obs ;
  double lamavg, n_calstar;
  char *cfilt ;
  char fnam[] = "INIT_ATMOSPHERE" ;

  // ------------ BEGIN ------------

  sprintf(BANNER,"%s to model DCR effects on RA, DEC, MAG \n", fnam);
  print_banner(BANNER);

  print_SURVEY(ID);

  ATMOS_INFO.PRESSURE    = SURVEY_INFO.pressure_atmos[ID] ;
  ATMOS_INFO.TEMPERATURE = SURVEY_INFO.temperature_atmos[ID] ;
  ATMOS_INFO.PWV         = SURVEY_INFO.pwv_atmos[ID] ;

  ATMOS_INFO.SNRMIN = 3.0 ;

  // for avg stellar wave per band, start with mean filter wave
  // (flat SED) 

  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) { 
    ATMOS_INFO.LAMAVG_CALSTAR[ifilt] = -9.0; 
    ATMOS_INFO.n_CALSTAR[ifilt] = -9.0; 
  }


  printf("   Mean wavelength & index of refraction per band:\n");
  // ifilt=0 is reserved for spectrograph, so passband ifilt starts at 1
  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++ ) {
    cfilt     = FILTER_SEDMODEL[ifilt].name ;
    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs;
    lamavg    = FILTER_SEDMODEL[ifilt].mean ;
    n_calstar = compute_index_refrac_atmos(lamavg, 0);

    ATMOS_INFO.LAMAVG_CALSTAR[ifilt_obs] = lamavg;
    ATMOS_INFO.n_CALSTAR[ifilt_obs] = n_calstar;

    printf("\t CalStar(%s) : <lam>=%7.1f A   <n-1>=%le\n",
	   cfilt, lamavg, n_calstar-1.0 );
    fflush(stdout);
  }

  //  debugexit(fnam);
  return;

} // end INIT_ATMOSPHERE

// ***************************************
void GEN_ATMOSPHERE_DRIVER(void) {

  // Created May 2023 by R.Kessler
  // Driver routine to simulate atmoshperic effects

  int NEPOCH = GENLC.NEPOCH ;
  int ep;
  char fnam[] = "GEN_ATMOSPHERE_DRIVER" ;

  // ------------ BEGIN ----------

  if ( INPUTS.ATMOSPHERE_OPTMASK == 0 ) { return; }

  if ( UNIT_TEST_COMPUTE_DCR ) { test_compute_dcr(); }

  // reset avg sums and other misc.
  INIT_EVENT_ATMOSPHERE();

  for ( ep = 1; ep <= GENLC.NEPOCH; ep++ )  {  
    if ( !GENLC.OBSFLAG_GEN[ep]  )  { continue ; }
    gen_airmass(ep);
    gen_dcr_coordShift(ep);
    genSmear_coords(ep);
  }

  // determine magShift after obs-weighted RA,DEC are determined.
  for ( ep = 1; ep <= GENLC.NEPOCH; ep++ )  {  
    if ( !GENLC.OBSFLAG_GEN[ep]  )  { continue ; }
    gen_dcr_magShift(ep);
  }
  
  return;

} // end GEN_ATMOSPHERE_DRIVER

// ====================================
void INIT_EVENT_ATMOSPHERE(void) {

  // Init stuff for each event
  int ifilt, ifilt_obs;
  char fnam[] = "INIT_EVENT_ATMOSPHERE";

  // ------------ BEGIN ----------

  reset_COORD_AVG(&ATMOS_INFO.COORD_RA);
  reset_COORD_AVG(&ATMOS_INFO.COORD_DEC);
  reset_COORD_AVG(&ATMOS_INFO.COORD_SIM_RA);
  reset_COORD_AVG(&ATMOS_INFO.COORD_SIM_DEC);

  return;

} // end reset_ATMOSPHERE_DRIVER

// ==============
void reset_COORD_AVG(COORD_AVG_DEF *COORD_AVG) {

  int ifilt, ifilt_obs ;
  char fnam[] = "reset_COORD_AVG" ;
  // ---------- BEGIN -------

  COORD_AVG->AVG = COORD_AVG->SUM = COORD_AVG->WGTSUM = 0.0 ;
  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++ ) {
    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs;
    
    COORD_AVG->AVG_BAND[ifilt_obs] = 0.0 ;
    COORD_AVG->SUM_BAND[ifilt_obs] = 0.0 ;
    COORD_AVG->WGTSUM_BAND[ifilt_obs] = 0.0 ;
  }

  return ;
} // end reset_COORD_AVG

void sum_COORD_AVG(COORD_AVG_DEF *COORD, 
		   double RA_OBS, double WGT, int IFILT_OBS) {

  // Increment sums for wgtd avg in COORD struct.
  char fnam[] = "sum_COORD_AVG" ;

  // --------- BEGIN -------

  COORD->SUM     += WGT * RA_OBS ;
  COORD->WGTSUM  += WGT ;

  COORD->SUM_BAND[IFILT_OBS]     += WGT * RA_OBS ;
  COORD->WGTSUM_BAND[IFILT_OBS]  += WGT;

  // re-compute AVG here after each obs so that we don't need a separate
  // code snippet at the end to compute final avg.
  COORD->AVG  = 
    COORD->SUM  / COORD->WGTSUM ;

  COORD->AVG_BAND[IFILT_OBS]  = 
    COORD->SUM_BAND[IFILT_OBS] / COORD->WGTSUM_BAND[IFILT_OBS] ;


  return ;

} // end sum_COORD_AVG

// ==========================================
void gen_airmass(int epoch) {

  // Created May 2023 by R.Kessler and Jason Lee (U.Penn)
  // Compute airmass for this epoch
  // 

  double MJD      = GENLC.MJD[epoch];
  double RA       = GENLC.RA ;  // true RA
  double DEC      = GENLC.DEC ; // true DEC
  int    IDSURVEY = GENLC.IDSURVEY;
  
  double geoLAT   = SURVEY_INFO.geoLAT[IDSURVEY] ;
  double geoLON   = SURVEY_INFO.geoLON[IDSURVEY] ;

  double RAD = RADIAN ;
  double airmass  = 0.01 ;
  double alt_rad, sin_alt, ang_zenith_rad, ang_zenith_deg;
  double h_hr, h_deg, cos_h ;
  double GLAT, GLON;
    
  // test example from ESO calculator:
  // https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
  bool DO_TEST = UNIT_TEST_AIRMASS ;
  double test_geoLAT  = -29.257 ;    // La Silla
  double test_geoLON  = -70.738 ;
  double test_MJD    = 59583.2409 ; // from LSST minion simlib
  double test_RA     = 149.;
  double test_DEC    = 2.2 ;

  char fnam[] = "gen_airmass" ;

  // ------------ BEGIN ------------

  GENLC.AIRMASS[epoch] = -9.0 ;
  GENLC.ALTITUDE[epoch] = -9.0 ;

  if ( DO_TEST ) {
    MJD       = test_MJD ;
    geoLAT    = test_geoLAT ;
    geoLON    = test_geoLON ;
    RA  = GENLC.RA  = test_RA;
    DEC = GENLC.DEC = test_DEC ;
    slaEqgal(RA, DEC, &GLON, &GLAT ); // return GLON, GLAT in degrees
    GENLC.sin_GLON = sin(GLON*RAD);
    GENLC.cos_GLON = cos(GLON*RAD);
    GENLC.sin_DEC  = sin(DEC*RAD);
    GENLC.cos_DEC  = cos(DEC*RAD);

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
  cos_h           = cos(h_deg*RAD) ;

  // compute airmass ...
  double sin_geoLAT = SURVEY_INFO.sin_geoLAT[IDSURVEY];
  double cos_geoLAT = SURVEY_INFO.cos_geoLAT[IDSURVEY];
  double sin_LON    = GENLC.sin_GLON ;
  double cos_LON    = GENLC.cos_GLON ;

  // avoid re-computing trig functions for each obs
  double sin_DEC    = GENLC.sin_DEC ;
  double cos_DEC    = GENLC.cos_DEC ; 

  sin_alt = 
    (sin_geoLAT * sin_DEC) + 
    (cos_geoLAT * cos_DEC * cos_h);

  alt_rad = asin(sin_alt);

  ang_zenith_rad = 0.25*TWOPI - alt_rad ; // zenight angle, radians
  ang_zenith_deg = ang_zenith_rad / RAD ;

  airmass = 1.0/cos(ang_zenith_rad) ;

  // store lots of useful info to avoid repeat compuations later
  GENLC.ALTITUDE[epoch]   = alt_rad/RAD; // store altitude in degrees
  GENLC.sin_ALT[epoch]    = sin_alt;
  GENLC.cos_ALT[epoch]    = cos(alt_rad);
  GENLC.AIRMASS[epoch]    = airmass;
  GENLC.ANG_ZENITH[epoch] = ang_zenith_deg ;
  GENLC.tan_ZENITH[epoch] = tan(ang_zenith_rad); // needed for DCR calc later

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

  double SNR_OBS = SEARCHEFF_DATA.SNR[epoch-1];

  int IFILT_OBS  = GENLC.IFILT_OBS[epoch];
  int detectFlag = SEARCHEFF_DATA.detectFlag[epoch-1] ; // regular C index


  bool VALID_DCR_SHIFT = ( GENLC.RA_dcr_shift[epoch] < COORD_SHIFT_NULL_DEG );

  double RA_OBS, DEC_OBS, RA_TRUE, DEC_TRUE ;
  double ANGRES, ran_RA, ran_DEC, WGT ;
  char fnam[] = "genSmear_coords" ;

  // ----------- BEGIN ---------

  ran_RA  = getRan_Gauss(1);
  ran_DEC = getRan_Gauss(1);

  ANGRES = ANGRES_REF_DEG * sqrt(SNR_REF/trueSNR);
  if ( !VALID_DCR_SHIFT ) { ANGRES = 0.0; } // nothing to smear

  // get true coords with DCR shift
  RA_TRUE  = GENLC.RA  + GENLC.RA_dcr_shift[epoch];
  DEC_TRUE = GENLC.DEC + GENLC.DEC_dcr_shift[epoch];
  
  // apply random smear to get observed coords
  RA_OBS  = RA_TRUE  + (ANGRES * ran_RA)/cosDEC;
  DEC_OBS = DEC_TRUE + (ANGRES * ran_DEC) ;

  // store observed and true coords
  GENLC.RA_OBS[epoch]  = RA_OBS;
  GENLC.DEC_OBS[epoch] = DEC_OBS;

  GENLC.RA_TRUE[epoch]  = RA_TRUE;
  GENLC.DEC_TRUE[epoch] = DEC_TRUE;
  
  // ?? what about uncertainty ???
  
  // update wgted-avg among all detctions
  bool USE_OBS   = SNR_OBS > ATMOS_INFO.SNRMIN ;
  if ( USE_OBS ) {

    if ( ANGRES > 0.0 ) 
      { WGT = (ANGRES_REF_DEG*ANGRES_REF_DEG) / (ANGRES*ANGRES); }
    else
      { WGT = 1.0E-20; }

    sum_COORD_AVG(&ATMOS_INFO.COORD_RA,  RA_OBS,  WGT, IFILT_OBS);
    sum_COORD_AVG(&ATMOS_INFO.COORD_DEC, DEC_OBS, WGT, IFILT_OBS);

    sum_COORD_AVG(&ATMOS_INFO.COORD_SIM_RA,  RA_TRUE,  WGT, IFILT_OBS);
    sum_COORD_AVG(&ATMOS_INFO.COORD_SIM_DEC, DEC_TRUE, WGT, IFILT_OBS);
  }


  return;

} // end genSmear_coords


// ========================================
void gen_dcr_coordShift(int ep) {

  // Created May 2023
  // Compute DCR astrometric shift for RA and DEC.

  double ALTITUDE   = GENLC.ALTITUDE[ep];
  double AIRMASS    = GENLC.AIRMASS[ep];
  double ANG_ZENITH = GENLC.ANG_ZENITH[ep];
  double tan_ZENITH = GENLC.tan_ZENITH[ep];
  int    IFILT_OBS  = GENLC.IFILT_OBS[ep];
  
  double RA         = GENLC.RA ;  // true RA 
  double DEC        = GENLC.DEC ; // true DEC   
  double LON        = GENLC.GLON; // SN coord
  double LAT        = GENLC.GLAT ;
  double sin_DEC    = GENLC.sin_DEC;
  double cos_DEC    = GENLC.cos_DEC ;

  int IDSURVEY      = GENLC.IDSURVEY ;
  double geoLAT     = SURVEY_INFO.geoLAT[IDSURVEY];  // telescope geo coord
  double geoLON     = SURVEY_INFO.geoLON[IDSURVEY] ;
  double geoALT     = SURVEY_INFO.geoALT[IDSURVEY] ;
  double sin_geoLAT = SURVEY_INFO.sin_geoLAT[IDSURVEY] ;
  double cos_geoLAT = SURVEY_INFO.cos_geoLAT[IDSURVEY] ;

  // define null RA,DEC shift if there is no SED model;
  // e.g., pre-explosion or at late times where model-mags are 
  // extrapolated and thus there is no SED
  double SHIFT_NULL = COORD_SHIFT_NULL_DEG ;
  double wave_sed_wgted;
  char fnam[] = "gen_dcr_coordShift" ;

  // -------------- BEGIN -----------

  GENLC.RA_dcr_shift[ep]  = SHIFT_NULL ;
  GENLC.DEC_dcr_shift[ep] = SHIFT_NULL ;

  // compute <wave> = integeral[lam*SED*Trans] / integ[SED*Trans]
  wave_sed_wgted = gen_wave_sed_wgted(ep);
  GENLC.LAMAVG_SED_WGTED[ep] = wave_sed_wgted ;

  if ( wave_sed_wgted < 0.01 ) { return; } // no model SED --> bail

  // set number of processed spectra (SEDs) to zero so that spectra/SED
  // are not written to data files. NMJD_PROC is reset next event.
  GENSPEC.NMJD_PROC = 0 ;

  // - - - - - - - - - - - - - - 
  // begin computatin of RA & DEC shifts
  double DCR, DCR_deg, sin_ALT, cos_ALT, q, cos_q, sin_q, cos_product;
  int DUMPFLAG = 0 ;

  // start with DCR angle shift in arcsec
  DCR = compute_DCR_angle(wave_sed_wgted, tan_ZENITH, IFILT_OBS, DUMPFLAG);
  DCR_deg = DCR / 3600.0 ;

  // compute cos and sin of paralatic angle
  sin_ALT = GENLC.sin_ALT[ep];
  cos_ALT = GENLC.cos_ALT[ep];

  cos_product = cos_DEC*cos_ALT ;
  if ( cos_product != 0.0 ) 
    { cos_q   = (sin_geoLAT - sin_DEC*sin_ALT) / cos_product; }
  else
    { cos_q   = 0.0 ; }

  q       = acos(cos_q);
  sin_q   = sin(q);

  // take projection for RA and DEC shifts in degrees
  GENLC.RA_dcr_shift[ep]  = DCR_deg * sin_q ;
  GENLC.DEC_dcr_shift[ep] = DCR_deg * cos_q ;

  return ;

} // end gen_dcr_coordShift

// ================================
double gen_wave_sed_wgted(int ep) {

  // Created May 2023
  // Return mean wavelegnth in filter corresponding to epoch 'ep';
  //    <wave> = integeral[lam*SED*Trans] / integ[SED*Trans]
  //

  double wave         = 0.0 ;
  int    NMJD_TOT     = GENSPEC.NMJD_TOT ;
  double MJD          = GENLC.MJD[ep];
  double TOBS         = GENLC.epoch_obs[ep];
  int    IFILT_OBS    = GENLC.IFILT_OBS[ep];
  int    IFILT        = IFILTMAP_SEDMODEL[IFILT_OBS] ;
  char   *FILTER_NAME = FILTER_SEDMODEL[IFILT].name ;
  int    imjd, IMJD = -9;
  
  int    ilam, NLAM_FILT;
  double MJD_SPEC, MJD_DIF_MIN, MJD_DIF ;
  char fnam[] = "gen_wave_sed_wgted";
  int LDMP    = 0 ;
 
  // --------- BEGIN ---------

  MJD_DIF_MIN = 1.0E9;
  for(imjd = 0; imjd < NMJD_TOT; imjd++ ) {
    MJD_SPEC = GENSPEC.MJD_LIST[imjd] ;
    MJD_DIF  = fabs(MJD_SPEC-MJD);
    if ( MJD_DIF < MJD_DIF_MIN ) {
      IMJD = imjd;  MJD_DIF_MIN = MJD_DIF;
    }
  } // end imjd

  if ( IMJD < 0 ) {
    sprintf(c1err,"Unable to find SED IMJD for MJD=%f" );
    sprintf(c2err,"ep=%d  MJD=%.4f  Tobs=%.2f", ep, MJD, TOBS);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  int    ifilt_obs_check;
  ifilt_obs_check = FILTER_SEDMODEL[IFILT].ifilt_obs;
  if ( ifilt_obs_check != IFILT_OBS ) {
    sprintf(c1err,"filter index mis-match");
    sprintf(c2err,"epoch IFILT_OBS=%d but FILTER_SEDMODEL(ifilt_obs)=%d",
	    IFILT_OBS, ifilt_obs_check);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // - - - - - - - -
  int NLAM_FILTER = FILTER_SEDMODEL[IFILT].NLAM;
  int NLAM_SED    = INPUTS_SPECTRO.NBIN_LAM;
  double lam, trans, sedFlux, ST ;
  double *ptr_SEDFLUX = GENSPEC.GENFLUX_LIST[IMJD];
  double *ptr_SEDLAM  = INPUTS_SPECTRO.LAMAVG_LIST;
  double sum0=0.0, sum1=0.0 ;
  for(ilam=0; ilam < NLAM_FILTER; ilam++ ) {
    lam   = FILTER_SEDMODEL[IFILT].lam[ilam];
    trans = FILTER_SEDMODEL[IFILT].transSN[ilam];
    sedFlux = interp_1DFUN(1, lam, NLAM_SED, ptr_SEDLAM, ptr_SEDFLUX, fnam);
    ST      = sedFlux * trans ;
    sum0 += ( ST ) ;
    sum1 += ( ST * lam );
  }

  if ( sum0 > 0.0 ) {  wave = sum1 / sum0; }

  if ( LDMP ) {
    printf(" xxx ---------------------------------- \n");
    printf(" xxx %s DUMP for CID=%d  NMJD_TOT=%d \n", 
	   fnam, GENLC.CID, NMJD_TOT );

    printf(" xxx MJD=%.3f  Tobs=%.3f  IFILTOBS=%d IFILT=%d"
	   "(ep=%d IMJD=%d) \n",
	   MJD, TOBS, IFILT_OBS, IFILT, ep, IMJD); fflush(stdout);
    printf(" xxx NLAM[FILTER,SED] = %d, %d \n",
	   NLAM_FILTER, NLAM_SED );
    printf(" xxx %s <wave> = %f \n",
	   FILTER_NAME, wave );
    fflush(stdout);

    if ( GENLC.CID > 2 ) { debugexit(fnam); }
  }


  return wave;

} // end gen_wave_sed_wgted


// =======================================
double compute_DCR_angle(double LAM, double tan_ZENITH, int IFILT_OBS, int DUMPFLAG) {

  // Created Jun 2023
  // Compute DCR from Eq 4 in Fillipenko 1982,
  //   https://articles.adsabs.harvard.edu/full/1982PASP...94..715F
  //
  // Inputs:
  //  LAM        : mean SED-weighted wavelength in passband
  //  tan_ZENITH :  tan(zenith angle)
  //  IFILT_OBS  : absolute passband index, used to find reference n-1
  //  DUMPFLAG   : optional dump flag

  double DCR = 0.0 ;
  double n_ref = ATMOS_INFO.n_CALSTAR[IFILT_OBS];
  double n_tele ;  // index of refrac for transient

  char fnam[] = "compute_DCR_angle" ;

  // ------------ BEGIN ----------

  n_tele = compute_index_refrac_atmos(LAM, DUMPFLAG);
  
  DCR = 206265.0 * ( n_tele - n_ref ) * tan_ZENITH ; // arcsec

  if ( DUMPFLAG ) {
    double z = atan(tan_ZENITH);
    double airmass = 1.0/cos(z);
    printf(" xxx %s: DCR = %f  (airmass=%f  tan_ZENITH = %.3f)\n", 
	   fnam, DCR, airmass, tan_ZENITH);
    fflush(stdout);
  }

  return DCR; // arcsec

} // end compute_DCR_angle


// ==========================
void test_compute_dcr(void) {

  // Created Jun 2023
  // Compute DCR on a grid of airmass and LAM to compute with 
  // Table 1 in Fillipenko 1982.

  int    IDSURVEY = 12; // LSST
  double LAMMIN_TEST=3000.0, LAMMAX_TEST=10000.0, LAMBIN_TEST=1000.0;
  double airmass, tanz, z, lam, dcr ;
  int    DUMPFLAG = 0 ;
  char fnam[] = "test_compute_dcr" ;

  // ---------- BEGIN ---------

  print_banner(fnam);
  // write table header
  printf("# Airmass  ");
  for(lam=LAMMIN_TEST ; lam < LAMMAX_TEST; lam+=LAMBIN_TEST ) 
    { printf(" %6.0f ", lam); }
  printf("\n# --------------------------------------"
	 "------------------------- \n");
  fflush(stdout);

  // - - - - - 
  for ( airmass = 1.0; airmass < 4.0; airmass += 1.0 ) {
    z = acos(1.0/airmass);
    tanz = tan(z);
    printf(" %6.3f    ", airmass);

    for(lam=LAMMIN_TEST ; lam < LAMMAX_TEST; lam+=LAMBIN_TEST ) {
      dcr = compute_DCR_angle(lam, tanz, 2, DUMPFLAG);
      printf(" %6.3f ", dcr); fflush(stdout);
    }
    printf("\n");
    fflush(stdout);

  } // end airmass loop

  debugexit(fnam);

  return;
} // end test_compute_dcr


// ========================================
double compute_index_refrac_atmos(double LAM, int DUMPFLAG) {
  
  // Created Jun 2023
  // Compute index of refraction for wavelength LAM in Angstroms,
  // and altitude ALT (meters) ... not to be confused with angle alt later.
  // Use Eqs 1,2,3 in Fillipenko 1982,                                         
  //   https://articles.adsabs.harvard.edu/full/1982PASP...94..715F 
  //
  double INVLAMSQ = 1.0E8/(LAM*LAM); // convert to inverse micron^2
  double n_0, n_1, n_tele;  // index of refrac at sea level, 2km height, +water

  // hard-code telescope altitude of 2km, but maybe later
  // need to add altitude argument to geo key in SURVEY.DEF
  double P_tele   = ATMOS_INFO.PRESSURE; // atmos pressure, mm Hg 
  double T_tele   = ATMOS_INFO.TEMPERATURE;   // temperature, Celsius
  double PWV_tele = ATMOS_INFO.PWV ;   // water vapor pressure, mm Hg

  double denom_T  = 1.0 + 0.003661*T_tele ;
  double tmp0, tmp1, tmp2;
  double ONE = 1.0;
  char fnam[] = "compute_index_refrac_atmos" ;

  // ---------- BEGIN ----------

  tmp0 = 64.328;
  tmp1 = 29498.1 / ( 146.0 - INVLAMSQ );
  tmp2 = 255.4 / (41.0 - INVLAMSQ);
  n_0  = ONE + (tmp0 + tmp1 + tmp2) * 1.0E-6 ;

  // correct for telescope altitude
  tmp0 = ( n_0 - ONE ) ;
  tmp1 = P_tele * (ONE + (1.049-0.0157*T_tele)*1.0E-6*P_tele);
  tmp2 = 720.883 * denom_T;
  n_1 = ONE + tmp0 * (tmp1/tmp2);

  // correct for water vapor
  tmp1 = (0.0624 - 0.000680*INVLAMSQ) * PWV_tele / denom_T ;
  n_tele = ONE + (n_1-ONE) - tmp1*1.0E-6 ;


  if ( DUMPFLAG && LAM != 5000.0 ) {
    printf(" xxx ----------- \n");
    printf(" xxx %s dump for LAM = %.1f A  (INVLAMSQ=%f)\n", 
	   fnam, LAM, INVLAMSQ);
    printf(" xxx n_0-1, n_1-1, n_tele-1 = %le %le %le (see level, 2km, +PWVcor) \n",
	   n_0-ONE, n_1-ONE, n_tele-ONE );
    fflush(stdout);
  }

  return n_tele ;

} // end compute_index_refrac_atmos

// ========================================
void gen_dcr_magShift(int ep) {

  // Compute mag shift for PSF-fitted flux where PSF centerl location
  // is offset from the band-average center.

  bool VALID_DCR_SHIFT = ( GENLC.RA_dcr_shift[ep] < COORD_SHIFT_NULL_DEG );
  char fnam[] = "gen_dcr_magShift" ;

  // ---------- BEGIN -------------

  GENLC.mag_dcr_shift[ep] = 0.0 ;
  if ( !VALID_DCR_SHIFT ) { return; }

  // - - - - 
  GENLC.mag_dcr_shift[ep] = 1.0E-4 ;

  return;

} // end gen_dcr_magShift

