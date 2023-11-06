/****************************************  
  Created May 2023 by R.Kessler 

  Tools to simulate atmosphere effects such as DCR(coords) and 
  DCR(PSF-shape) ; 
  Motivated by Le et al, 2023: https://arxiv.org/abs/2304.01858

  These functions lean heavily on INPUTS and GENLC structures
  in snlc_sim.h, so these tools are not easily plugged into
  other codes.

  To DO:
  - determine reference wavelength from calib stars;
    computation or input numbers in the sim-input file ?
  - fluctuations in pressure, temperature, pwv 
  - noise model for RA,DEC precision

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
  int  OPTMASK = INPUTS_ATMOSPHERE.OPTMASK;
  int ifilt, ifilt_obs, NLAM, ilam, MEMD ;
  double lamavg_flat, lamavg_calstar, n_calstar, lam, n_site;
  char *cfilt ;
  char fnam[] = "INIT_ATMOSPHERE" ;

  // ------------ BEGIN ------------

  sprintf(BANNER,"%s to model DCR effects on RA, DEC, MAG", fnam);
  print_banner(BANNER);

  INPUTS_ATMOSPHERE.DO_DCR_COORD    = (OPTMASK & ATMOSPHERE_OPTMASK_DCR_COORD) > 0;
  INPUTS_ATMOSPHERE.DO_DCR_PSFSHAPE = (OPTMASK & ATMOSPHERE_OPTMASK_DCR_PSFSHAPE) > 0;

  printf("\t DO_DCR_COORD    = %d \n", INPUTS_ATMOSPHERE.DO_DCR_COORD);
  printf("\t DO_DCR_PSFSHAPE = %d \n", INPUTS_ATMOSPHERE.DO_DCR_PSFSHAPE);
  fflush(stdout);

  print_SURVEY(ID);


  if ( SURVEY_INFO.geoALT[ID] < 0.001 ) {
    sprintf(c1err,"Survey ID=%d missing geoLat:geoLONG:ALT", ID);
    sprintf(c2err,"See SURVEY.DEF file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );          
  }

  ATMOS_INFO.PRESSURE_AVG    = SURVEY_INFO.pressure_atmos[ID] ;
  ATMOS_INFO.TEMPERATURE_AVG = SURVEY_INFO.temperature_atmos[ID] ;
  ATMOS_INFO.PWV_AVG         = SURVEY_INFO.pwv_atmos[ID] ;

  ATMOS_INFO.SNRMIN = 3.0 ;

  printf("\t Sigma(temperature/Pressure/PWV) = %.1f C / %.1f mmHg / %.1f mmHg\n\n",
	 INPUTS_ATMOSPHERE.SIGMA_SITE_TEMP,
	 INPUTS_ATMOSPHERE.SIGMA_SITE_BP,
	 INPUTS_ATMOSPHERE.SIGMA_SITE_PWV ); fflush(stdout);
  INPUTS_ATMOSPHERE.APPLY_SIGMA_SITE = 
    ( INPUTS_ATMOSPHERE.SIGMA_SITE_TEMP > 0.0 ||
      INPUTS_ATMOSPHERE.SIGMA_SITE_BP   > 0.0 ||
      INPUTS_ATMOSPHERE.SIGMA_SITE_PWV  > 0.0 );


  read_stellar_sed_atmos();

  // for avg stellar wave per band, start with mean filter wave
  // (flat SED) 

  for(ifilt=0; ifilt < MXFILTINDX; ifilt++ ) { 
    ATMOS_INFO.LAMAVG_CALSTAR[ifilt] = -9.0; 
    ATMOS_INFO.n_CALSTAR_AVG[ifilt]  = -9.0; 
  }


  printf("\t                              mean   \n");
  printf("\t         flatSED  calStar    calStar \n");
  printf("\t  band    <lam>    <lam>      <n-1>  \n");
  printf("\t# ------------------------------------------------- \n");

  // ifilt=0 is reserved for spectrograph, so passband ifilt starts at 1
  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++ ) {
    cfilt     = FILTER_SEDMODEL[ifilt].name ;
    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs;

    init_stellar_sed_atmos(ifilt_obs);

    lamavg_flat       = FILTER_SEDMODEL[ifilt].mean ;    
    lamavg_calstar    = ATMOS_INFO.LAMAVG_CALSTAR[ifilt_obs] ;
    n_calstar         = compute_index_refrac_atmos(lamavg_calstar, 0);

    ATMOS_INFO.n_CALSTAR_AVG[ifilt_obs]  = n_calstar ;

    printf("\t %s   %7.1f  %7.1f  %le\n",
	   cfilt, lamavg_flat, lamavg_calstar, n_calstar-1.0 );
    fflush(stdout);

    // if there are no site fluctuations, store n_site[ifilt_obs][ilam]    
    if ( !INPUTS_ATMOSPHERE.APPLY_SIGMA_SITE ) {
      NLAM = FILTER_SEDMODEL[ifilt].NLAM;
      MEMD = NLAM * sizeof(double);
      ATMOS_INFO.n_SITE_LIST[ifilt_obs] = (double*)malloc(MEMD);
      // xxx printf("\t xxx Store n(lam) for %s-band \n", cfilt); fflush(stdout);
      for(ilam=0; ilam < NLAM; ilam++ ) {
	lam    = FILTER_SEDMODEL[ifilt].lam[ilam];
	n_site = compute_index_refrac_atmos(lam, 0);    
	ATMOS_INFO.n_SITE_LIST[ifilt_obs][ilam] = n_site;
      }
    }

  } // end ifilt



  // - - - - - - -
  // check polynomial functions to characterize a few things,
  // and print some values to stdout for visual inspection.

  bool REQUIRE_RESPOLY = INPUTS_ATMOSPHERE.DO_DCR_COORD ;
  bool REQUIRE_MAGPOLY = INPUTS_ATMOSPHERE.DO_DCR_COORD ;

  if ( REQUIRE_RESPOLY ) { 

    double SNR, ANGRES, STATRES;
    GENPOLY_DEF *RESPOLY = &INPUTS_ATMOSPHERE.DCR_COORDRES_POLY;

    if ( RESPOLY->ORDER < 0 ) {
      sprintf(c1err,"Missing required coord res vs. PSF/SNR");
      sprintf(c2err,"Set sim-input key %s", KEYNAME_ATMOSPHERE_DCR_COORDRES_POLY );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );       
    }

    // summarize polyFun to describe astrometry resolution vs. 1/sqrt(SNR)
    printf("\n");
    print_GENPOLY(RESPOLY); 
    for (SNR=10.0; SNR <= 100.0; SNR+= 30.0 ) {
      double PSF_FWHM = 1.0 ; // arcsec
      double x = PSF_FWHM/SNR;
      STATRES = eval_GENPOLY(x, RESPOLY, fnam);
      printf("\t STATRES = %7.4f arcsec for SNR = %4.0f (PSF_FWHM=%.1f asec)\n", 
	     STATRES, SNR, PSF_FWHM); fflush(stdout);
    }
  }  // end REQUIRE_RESPOLY

  if ( REQUIRE_MAGPOLY ) {

    double fracPSF, mag_shift;
    GENPOLY_DEF *MAGPOLY = &INPUTS_ATMOSPHERE.DCR_MAGSHIFT_POLY;
    if ( MAGPOLY->ORDER < 0 ) {
      sprintf(c1err,"Missing required coord mag vs. PSF-shift-fraction.");
      sprintf(c2err,"Set sim-input key %s", KEYNAME_ATMOSPHERE_DCR_MAGSHIFT_POLY );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );       
    }

    // summarize polyFun to describe mag offset vs. PSF-fraction-shift
    printf("\n");
    print_GENPOLY(MAGPOLY);
    for(fracPSF = 0.0; fracPSF <= 0.2; fracPSF += 0.04 ) {
      mag_shift = eval_GENPOLY(fracPSF, MAGPOLY, fnam);
      printf("\t mag_shift = %7.4f mag for PSFshift/PSF = %.4f \n", 
	     mag_shift, fracPSF); fflush(stdout);
      
    }
  } // end REQUIRE_MAGPOLY



  printf("\n\t Finished %s \n", fnam); fflush(stdout);

  //  debugexit(fnam);
  return;

} // end INIT_ATMOSPHERE

// ********************************************
void  read_stellar_sed_atmos(void) {

  // Created JUn 2023
  // open stellar SED file and read/store contents

  char *ptrFile = INPUTS_ATMOSPHERE.SEDSTAR_FILE;
  int    MEMD = sizeof(double);
  char fnam[] = "read_stellar_sed_atmos" ;

  // ----------- BEGIN -----------

  printf("   Read average calStar SED from : %s\n", ptrFile);

  ENVreplace(ptrFile, fnam, 1);

  ATMOS_INFO.LAM_ARRAY_CALSTAR  = (double*) malloc(MXBIN_LAMSED_SEDMODEL * MEMD);
  ATMOS_INFO.FLUX_ARRAY_CALSTAR = (double*) malloc(MXBIN_LAMSED_SEDMODEL * MEMD);

  rd2columnFile(ptrFile, MXBIN_LAMSED_SEDMODEL, &ATMOS_INFO.NBINLAM_CALSTAR,
		ATMOS_INFO.LAM_ARRAY_CALSTAR, ATMOS_INFO.FLUX_ARRAY_CALSTAR );

  int NB = ATMOS_INFO.NBINLAM_CALSTAR ;
  printf("\t Found %d wave bins from %.0f to %.0f A \n",  NB,
	 ATMOS_INFO.LAM_ARRAY_CALSTAR[0],
	 ATMOS_INFO.LAM_ARRAY_CALSTAR[NB-1] );

  printf("\n");
  fflush(stdout);

  return;

} // end read_stellar_sed_atmos

// ********************************************
void init_stellar_sed_atmos(int ifilt_obs) {
  
  // store FLUX-vs-lam on filter-lam grid to speed up integrals.
  // Also store diagnostic info such as <lam>.

  int    ifilt        = IFILTMAP_SEDMODEL[ifilt_obs];
  int    NLAM_FILTER  = FILTER_SEDMODEL[ifilt].NLAM;
  int    NLAM_CALSTAR = ATMOS_INFO.NBINLAM_CALSTAR ;
  int    MEMD         = NLAM_FILTER * sizeof(double);
  int ilam;
  double lamavg, lam, trans, flux_star, sum0=0.0, sum1=0.0 ;
  char fnam[] = "init_stellar_sed_atmos" ;

  // ---------- BEGIN ----------

  lamavg = 0.0 ;

  ATMOS_INFO.FLUX_CALSTAR[ifilt_obs] = (double*) malloc(MEMD);

  for(ilam=0; ilam < NLAM_FILTER; ilam++ ) {
    lam   = FILTER_SEDMODEL[ifilt].lam[ilam];
    trans = FILTER_SEDMODEL[ifilt].transSN[ilam];
    flux_star = interp_1DFUN(1, lam, NLAM_CALSTAR,
			     ATMOS_INFO.LAM_ARRAY_CALSTAR,
			     ATMOS_INFO.FLUX_ARRAY_CALSTAR, fnam);

    ATMOS_INFO.FLUX_CALSTAR[ifilt_obs][ilam] = flux_star;

    sum0 += (flux_star * trans);
    sum1 += (flux_star * trans * lam);
  }

  lamavg = sum1/sum0;

  ATMOS_INFO.LAMAVG_CALSTAR[ifilt_obs] = lamavg;

  return ;

} // end init_stellar_sed_atmos

// ***************************************
void GEN_ATMOSPHERE_DRIVER(void) {

  // Created May 2023 by R.Kessler
  // Driver routine to simulate atmoshperic effects

  int NEPOCH = GENLC.NEPOCH ;
  int ep;
  char fnam[] = "GEN_ATMOSPHERE_DRIVER" ;

  // ------------ BEGIN ----------

  if ( INPUTS_ATMOSPHERE.OPTMASK == 0 ) { return; }

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
  double h_hr, h_deg, cos_h, sin_h ;
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
  sin_h           = sin(h_deg*RAD) ;

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
  GENLC.sin_h[epoch]      = sin_h; // needed for ALT shift to RA shift later

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
  // Determine measured RA,DEC for this epoch, using coord-smear 
  // resolution defined by sim-input key ATMOSPHERE_DCR_COORDRES_POLY and ATMOSPHERE_DCR_COORDRES_FLOOR

  double cosDEC  = GENLC.cosDEC;
  double SNR_TRUE = SEARCHEFF_DATA.SNR_CALC[epoch-1];
  if ( SNR_TRUE < 0.01 ) { SNR_TRUE = 0.01; }
  double SNR_OBS  = SEARCHEFF_DATA.SNR_OBS[epoch-1];
  double PSF_FWHM = SIMLIB_OBS_GEN.PSF_FWHM[epoch];

  int IFILT_OBS  = GENLC.IFILT_OBS[epoch];
  int detectFlag = SEARCHEFF_DATA.detectFlag[epoch-1] ; // regular C index


  bool VALID_DCR_SHIFT = ( GENLC.RA_dcr_shift[epoch] < COORD_SHIFT_NULL_DEG );

  GENPOLY_DEF  *DCR_COORDRES_POLY = &INPUTS_ATMOSPHERE.DCR_COORDRES_POLY;
  double DCR_COORDRES_FLOOR = INPUTS_ATMOSPHERE.DCR_COORDRES_FLOOR;
  double RA_OBS, DEC_OBS, RA_TRUE, DEC_TRUE ;
  double STATRES_TRUE_asec, STATRES_OBS_asec;
  double ANGRES_TRUE_asec, ANGRES_TRUE_deg, ran_RA, ran_DEC, WGT ;
  double ANGRES_OBS_asec,  ANGRES_OBS_deg;
  char fnam[] = "genSmear_coords" ;

  // ----------- BEGIN ---------

  ran_RA  = getRan_Gauss(1);
  ran_DEC = getRan_Gauss(1);


  double x_TRUE = PSF_FWHM/SNR_TRUE;
  double x_OBS  = PSF_FWHM/SNR_OBS;

  STATRES_TRUE_asec = eval_GENPOLY(x_TRUE, DCR_COORDRES_POLY, fnam);
  STATRES_OBS_asec = eval_GENPOLY(x_OBS, DCR_COORDRES_POLY, fnam);
  ANGRES_TRUE_asec = sqrt(pow(STATRES_TRUE_asec, 2) + pow(DCR_COORDRES_FLOOR, 2)); //sqrt(sigma_syst**2 + sigma_stat**2)
  ANGRES_OBS_asec  = sqrt(pow(STATRES_OBS_asec, 2) + pow(DCR_COORDRES_FLOOR, 2));
  if ( !VALID_DCR_SHIFT ) { ANGRES_TRUE_asec = 0.0; } // nothing to smear

  // convert ANGRES to degrees and divide by sqrt(2) for 
  // projected resolution on separate RA and DEC coords.
  double unit_convert = 1.0 / (3600.0*1.4142);
  ANGRES_TRUE_deg = ANGRES_TRUE_asec * unit_convert;
  ANGRES_OBS_deg  = ANGRES_OBS_asec  * unit_convert;

  // get true coords with DCR shift
  RA_TRUE  = GENLC.RA  + GENLC.RA_dcr_shift[epoch];
  DEC_TRUE = GENLC.DEC + GENLC.DEC_dcr_shift[epoch];
  
  // apply random smear to get observed coords
  RA_OBS  = RA_TRUE  + (ANGRES_TRUE_deg * ran_RA)/cosDEC;
  DEC_OBS = DEC_TRUE + (ANGRES_TRUE_deg * ran_DEC) ;

  // store observed and true coords
  GENLC.RA_OBS[epoch]  = RA_OBS;
  GENLC.DEC_OBS[epoch] = DEC_OBS;

  GENLC.RA_TRUE[epoch]  = RA_TRUE;
  GENLC.DEC_TRUE[epoch] = DEC_TRUE;
  
  // ?? what about uncertainty ???
  
  // update wgted-avg among all detctions
  bool USE_OBS   = SNR_OBS > ATMOS_INFO.SNRMIN ;
  if ( USE_OBS ) {

    if ( ANGRES_TRUE_deg > 0.0 ) 
      { WGT = 1.0E-6 / (ANGRES_OBS_deg*ANGRES_OBS_deg); }
    else
      { WGT = 1.0E-20; }

    ATMOS_INFO.COORDRES[epoch] = ANGRES_OBS_asec;

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

  int  OPTMASK      = INPUTS_ATMOSPHERE.OPTMASK;
  bool DO_DCR_COORD = (OPTMASK & ATMOSPHERE_OPTMASK_DCR_COORD) > 0 ;

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
  double sin_h      = GENLC.sin_h[ep];

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

  if ( !DO_DCR_COORD ) {
    GENLC.dcr_shift[ep]     = 0.0 ;
    GENLC.RA_dcr_shift[ep]  = 0.0 ;
    GENLC.DEC_dcr_shift[ep] = 0.0 ;
    return ;
  }

  // init dcr coord shifts to crazy value to indicate SED model not available.
  GENLC.dcr_shift[ep]     = SHIFT_NULL ;
  GENLC.RA_dcr_shift[ep]  = SHIFT_NULL ;
  GENLC.DEC_dcr_shift[ep] = SHIFT_NULL ;


  // - - - - - - - - - - - - - - 
  // begin computatin of RA & DEC shifts
  double DCR, DCR_deg, sin_ALT, cos_ALT, cos_q, sin_q, cos_product, sin_AZ;
  int DUMPFLAG = 0 ;

  // start with DCR angle shift in arcsec
  DCR = compute_DCR_angle(ep, DUMPFLAG);
  DCR_deg = DCR / 3600.0 ;

  if ( DCR >=  999.0 ) { return; }

  // compute cos and sin of paralatic angle
  sin_ALT = GENLC.sin_ALT[ep];
  cos_ALT = GENLC.cos_ALT[ep];

  if ( cos_ALT != 0.0 )
    { sin_AZ   = (-sin_h*cos_DEC) / cos_ALT; }
  else
    { sin_AZ   = 0.0 ; }

  if ( cos_DEC != 0.0 )
    { sin_q   = (-cos_geoLAT*sin_AZ) / cos_DEC; }
  else
    { sin_q   = 0.0 ; }

  cos_product = cos_DEC*cos_ALT ;
  if ( cos_product != 0.0 ) 
    { cos_q   = (sin_geoLAT - sin_DEC*sin_ALT) / cos_product; }
  else
    { cos_q   = 0.0 ; }

  // take projection for RA and DEC shifts in degrees, update 10/3/2023 J.L. dividing RA shift by cos_DEC
  GENLC.dcr_shift[ep]     = DCR_deg ;
  GENLC.RA_dcr_shift[ep]  = DCR_deg * sin_q / cos_DEC ;
  GENLC.DEC_dcr_shift[ep] = DCR_deg * cos_q ;

  // set number of processed spectra (SEDs) to zero so that spectra/SED
  // are not written to data files ... unless sim-input explicitly
  // requests SED output. NMJD_PROC is reset next event.
  if ( (INPUTS.WRITE_MASK & WRITE_MASK_SPECTRA) == 0 ) 
    { GENSPEC.NMJD_PROC = 0 ; }

  return ;

} // end gen_dcr_coordShift



// =======================================
double compute_DCR_angle(int ep, int DUMPFLAG) {

  // Created Jun 2023
  // Compute DCR from Eq 4 in Fillipenko 1982,
  //   https://articles.adsabs.harvard.edu/full/1982PASP...94..715F
  // and integrate over filter wavelengths as in Eq 1 of 
  //   Plazas & Bernstein 2012,
  //   https://iopscience.iop.org/article/10.1086/668294/meta
  //
  // Inputs:
  //  ep         : epoch index
  //  DUMPFLAG   : optional dump flag

  double DCR_SED_WGTED    = 999.0 ; // init output
  double LAMAVG_SED_WGTED =   0.0 ; // init output

  int    IFILT_OBS      = GENLC.IFILT_OBS[ep];
  int    IFILT          = IFILTMAP_SEDMODEL[IFILT_OBS] ;
  int    NMJD_TOT       = GENSPEC.NMJD_TOT ;
  double MJD            = GENLC.MJD[ep];
  double TOBS           = GENLC.epoch_obs[ep];
  double tan_ZENITH     = GENLC.tan_ZENITH[ep];
  double FAC_DCR        = 206265.0;

  int    IMJD ;
  double n_site;  // index of refrac at site

  char fnam[] = "compute_DCR_angle" ;

  // ------------ BEGIN ----------

  // determine which SED based on MJD matching
  IMJD = IMJD_GENSPEC(MJD);

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

  int NLAM_FILTER = FILTER_SEDMODEL[IFILT].NLAM;
  int NLAM_SED    = INPUTS_SPECTRO.NBIN_LAM;
  double *ptrFLUXSED    = GENSPEC.GENFLUX_LIST[IMJD];   // transient SED
  double *ptrLAMSED     = INPUTS_SPECTRO.LAMAVG_LIST;
  double  LAMSED_BINSIZE = INPUTS_SPECTRO.LAMBIN_LIST[0];

  double sum0_SED=0.0, sum1_SED=0.0 ;
  double sum0_CAL=0.0, sum1_CAL=0.0 ;
  double sum0_LAM=0.0, sum1_LAM=0.0 ;
  double lam, trans, sedFlux, calFlux, ST, DCR_lam; 
  int    ilam ;

  for(ilam=0; ilam < NLAM_FILTER; ilam++ ) {
    lam   = FILTER_SEDMODEL[IFILT].lam[ilam];
    trans = FILTER_SEDMODEL[IFILT].transSN[ilam];

    // interpolate SED flux 
    sedFlux = interp_1DFUN(1, lam, NLAM_SED, ptrLAMSED, ptrFLUXSED, fnam);

    // calstar flux already interpolated and stored on filter-lam grid
    calFlux = ATMOS_INFO.FLUX_CALSTAR[IFILT_OBS][ilam]; 

    if ( INPUTS_ATMOSPHERE.APPLY_SIGMA_SITE ) 
      { n_site = compute_index_refrac_atmos(lam, DUMPFLAG);  }
    else
      { n_site = ATMOS_INFO.n_SITE_LIST[IFILT_OBS][ilam] ; }

    DCR_lam = FAC_DCR * (n_site-1.0) * tan_ZENITH ; // arcsec

    ST      = sedFlux * trans ;
    sum0_SED += ( ST ) ;
    sum1_SED += ( ST * DCR_lam );

    sum0_LAM += ( ST );
    sum1_LAM += ( ST * lam);

    ST      = calFlux * trans ;
    sum0_CAL += ( ST ) ;
    sum1_CAL += ( ST * DCR_lam );

  } // end ilam loop
  
  if ( sum0_SED > 0.0 && sum0_CAL > 0.0 ) { 
    double DCR_SED    = sum1_SED/sum0_SED ;  // transient 
    double DCR_CAL    = sum1_CAL/sum0_CAL;   // avg calib star
    DCR_SED_WGTED     = DCR_SED - DCR_CAL ; 
    LAMAVG_SED_WGTED  = sum1_LAM / sum0_LAM ;

    /* xxx
    printf(" xxx %s: DCR[SED,CAL,net] = %f %f %f  \n",
	   fnam, DCR_SED, DCR_CAL, DCR ); fflush(stdout);
    debugexit(fnam);
    */
  }


  if ( DUMPFLAG ) {
    double z = atan(tan_ZENITH);
    double airmass = 1.0/cos(z);
    printf(" xxx %s: DCR = %f  (airmass=%f  tan_ZENITH = %.3f)\n", 
	   fnam, DCR_SED_WGTED, airmass, tan_ZENITH);
    fflush(stdout);
  }

  GENLC.LAMAVG_SED_WGTED[ep] = LAMAVG_SED_WGTED ;// diagnostic; not used for calc.
  return DCR_SED_WGTED ; // arcsec

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
      //      dcr = compute_DCR_angle(ep, DUMPFLAG); // what is ep ??
      dcr = -9.0 ;
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

  // Telescope params
  double TEMP_tele   = ATMOS_INFO.TEMPERATURE_AVG;   // temperature, Celsius
  double BP_tele     = ATMOS_INFO.PRESSURE_AVG; // atmos pressure, mm Hg 
  double PWV_tele    = ATMOS_INFO.PWV_AVG ;   // water vapor pressure, mm Hg

  double denom_T, tmp0, tmp1, tmp2, r ;
  double ONE = 1.0;
  char fnam[] = "compute_index_refrac_atmos" ;

  // ---------- BEGIN ----------


  // This site condition model is way too extreme because it doesn't
  // account for weather correlations.
  if ( INPUTS_ATMOSPHERE.APPLY_SIGMA_SITE ) {
    r = getRan_GaussClip(1,-3.0, 3.0);
    TEMP_tele += r * INPUTS_ATMOSPHERE.SIGMA_SITE_TEMP ;

    r = getRan_GaussClip(1,-3.0, 3.0);
    BP_tele += r * INPUTS_ATMOSPHERE.SIGMA_SITE_BP ;

    r = getRan_GaussClip(1,-3.0, 3.0);
    PWV_tele += r * INPUTS_ATMOSPHERE.SIGMA_SITE_PWV ;
  }

  denom_T  = 1.0 + 0.003661*TEMP_tele ;

  tmp0 = 64.328 ;
  tmp1 = 29498.1 / ( 146.0 - INVLAMSQ );
  tmp2 = 255.4 / (41.0 - INVLAMSQ);
  n_0  = ONE + (tmp0 + tmp1 + tmp2) * 1.0E-6 ;

  // correct for telescope altitude
  tmp0 = ( n_0 - ONE ) ;
  tmp1 = BP_tele * (ONE + (1.049-0.0157*TEMP_tele)*1.0E-6*BP_tele);
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

  // Compute mag shift for PSF-fitted flux where PSF center location
  // is offset from the band-average center.

  double dcr_shift_deg  = GENLC.dcr_shift[ep]; 
  double dcr_shift_asec = dcr_shift_deg*3600.0 ;
  bool VALID_DCR_SHIFT  = (  dcr_shift_deg < COORD_SHIFT_NULL_DEG );

  int    LDMP = ( GENLC.CID == -2 ) ;
  double fracPSF, PSF_FWHM, mag_shift ;
  char fnam[] = "gen_dcr_magShift" ;

  // ---------- BEGIN -------------

  GENLC.mag_dcr_shift[ep] = 0.0 ;
  if ( !VALID_DCR_SHIFT ) { return; }

  // - - - - - - - - - - - -
  PSF_FWHM = SIMLIB_OBS_GEN.PSF_FWHM[ep];
  fracPSF  = fabs(dcr_shift_asec) / PSF_FWHM ;

  if ( INPUTS_ATMOSPHERE.DO_DCR_COORD ) {
    // compute mag shift from polynominal fit to PSF-fitted flux using galsim.
    GENPOLY_DEF *MAGPOLY = &INPUTS_ATMOSPHERE.DCR_MAGSHIFT_POLY;
    mag_shift            = eval_GENPOLY(fracPSF, MAGPOLY, fnam);
    GENLC.mag_dcr_shift[ep] += mag_shift ; // store global
  }

  if ( LDMP ) {
    int    IFILT_OBS      = GENLC.IFILT_OBS[ep];
    double MJD            = GENLC.MJD[ep];
    printf(" xxx ---------------------- \n");
    printf(" xxx %s: CID=%d  MJD=%.4f  IFILTOBS=%d  PSF_FWHM=%.3f asec\n",
	   fnam, GENLC.CID, MJD, IFILT_OBS, PSF_FWHM );
    printf("\t xxx PSF_FWHM=%5.3f  dcr_shift=%7.4f  frac=%.4f  magShift=%.4f\n",
	   PSF_FWHM, dcr_shift_asec, fracPSF, mag_shift);
    fflush(stdout);
  }

  return;

} // end gen_dcr_magShift

