//
// Collection of test functions for simulations. These test
// functions call other functions with controlled inputs
// and print results. Do not put dump-utilities here.
// This file is included in snlc_sim.c, so there is no
// sntools_sim_testfuns.h and no object file.
//
// [moved out of snlc_sim.c on Nov 19 2022]
//
//

// ================================================

void unit_test_driver(char *UNIT_TEST_NAME);

void test_kcor_utils(void);
void test_kcor_table_lcmag(void);
void test_kcor_table_mwxt(void);
void test_kcor_table_avwarp(void);
void test_kcor_table_kcor(void);
void test_kcor_nearfilt_rest(void);
void test_GET_KCOR_DRIVER(void);

void test_fortran(void);
void test_igm(void);
void test_ran(void);
void test_PARSE_WORDS(void);
void test_zcmb_dLmag_invert(void);

void test_getRan_funVal(char *FUNVAL_NAME);
void load_test_GENGAUSS(GENGAUSS_ASYM_DEF *GENGAUSS );
void load_test_GENEXP(GEN_EXP_HALFGAUSS_DEF *GENEXP);

char TEST_REFAC[]  = "REFAC";
char TEST_LEGACY[] = "LEGACY";

#define STRING_FUNVAL_GENGAUSS     "FUNVAL_GENGAUSS"
#define STRING_FUNVAL_GENEXP       "FUNVAL_GENEXP"
#define STRING_FUNVAL_GENPDF       "FUNVAL_GENPDF"

// ************************


void unit_test_driver(char *UNIT_TEST_NAME) {

  char fnam[] = "unit_test_driver" ;

  // ------ BEGIN ----------

  if ( strlen(UNIT_TEST_NAME) == 0 ) { return; }

  if ( strcmp(UNIT_TEST_NAME,"IGM") == 0 ) 
    { test_igm(); }

  else if ( strstr(UNIT_TEST_NAME,"FUNVAL") != NULL ) 
    { test_getRan_funVal(UNIT_TEST_NAME); }

  else {
    sprintf(c1err,"Undefined UNIT_TEST: %s", UNIT_TEST_NAME);
    sprintf(c2err,"Check UNIT_TEST key in sim-input file");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

} // end unit_test_driver


// =====================================
void test_kcor_utils(void) {
  
  // Created Nov 2022

  char fnam[] = "test_kcor_utils" ;

  // --------------- BEGIN ------------

  // test_kcor_table_lcmag();
  // test_kcor_table_mwxt();
  // test_kcor_nearfilt_rest();

  test_kcor_table_avwarp();
  test_kcor_table_kcor();
  test_GET_KCOR_DRIVER();

  debugexit(fnam);
  return;

} // end test_kcor_utils

void test_kcor_table_lcmag(void) {
  
  // Created Nov 2022

  FILTERCAL_DEF *FILTERCAL_REST = &CALIB_INFO.FILTERCAL_REST ;
  int    NFILTDEF_REST          = FILTERCAL_REST->NFILTDEF ;  

  double Trest=-7.7, z=ZAT10PC, AVwarp=0.0, MAG ;
  int  ifilt_r, ifiltdef_r, istat;
  char *filtName_r, *text  ;
  char fnam[] = "test_kcor_table_lcmag" ;

  // --------------- BEGIN ------------

  printf("\n");
  printf(" xxx -------------------------------------------- \n");
  printf(" xxx %s: \n",fnam);
  printf("\n  Trest    AVwarp       Filter     LCMAG   istat\n");
  for(ifilt_r=0; ifilt_r < NFILTDEF_REST; ifilt_r++ ) {
    ifiltdef_r = FILTERCAL_REST->IFILTDEF[ifilt_r];
    filtName_r = FILTERCAL_REST->FILTER_NAME[ifilt_r];
    istat = 0;

    for ( AVwarp=-6.0; AVwarp <= 6.0; AVwarp+=3.0 ) {


      MAG = eval_kcor_table_LCMAG(ifiltdef_r, Trest, z, AVwarp) ;
      text = TEST_REFAC ;
      printf("  %6.1f  %6.3f  %12s  %7.4f  %d  (%s)\n",
             Trest, AVwarp, filtName_r, MAG, istat, text );
      fflush(stdout);
    }
  }

  return;

} // end test_kcor_table_lcmag


void test_kcor_table_mwxt(void) {
  
  // Created Nov 2022

  FILTERCAL_DEF *FILTERCAL_OBS  = &CALIB_INFO.FILTERCAL_OBS ;
  int    NFILTDEF_OBS           = FILTERCAL_OBS->NFILTDEF ;  

  double Trest=9.4, z=ZAT10PC, AVwarp=0.0, MWXT ;
  double MWEBV = 0.05;
  double RV    = 3.1;
  int    OPT_MWCOLORLAW = 94;

  int  ifilt_o, ifiltdef_o, istat;
  char *filtName_o, *text  ;
  char fnam[] = "test_kcor_table_mwxt" ;

  // --------------- BEGIN ------------

  printf("\n");
  printf(" xxx -------------------------------------------- \n");
  printf(" xxx %s: \n",fnam);
  printf("\n  Trest   AVwarp       Filter     MWXT   istat\n");
  for(ifilt_o=0; ifilt_o < NFILTDEF_OBS; ifilt_o++ ) {
    ifiltdef_o = FILTERCAL_OBS->IFILTDEF[ifilt_o];
    filtName_o = FILTERCAL_OBS->FILTER_NAME[ifilt_o] ;
    istat = 0;

    for ( AVwarp=-1.087; AVwarp < 3.0; AVwarp+=1.6 ) {


      MWXT = eval_kcor_table_MWXT(ifiltdef_o, Trest, z, AVwarp, 
				  MWEBV, RV, OPT_MWCOLORLAW) ;
      text = TEST_REFAC ;

      printf("  %6.1f  %6.3f  %12s  %7.4f  %d  (%s)\n",
             Trest, AVwarp, filtName_o, MWXT, istat, text ); 
      fflush(stdout);
    }
  }

  return;

} // end test_kcor_table_mwxt


// =================================================
void test_kcor_nearfilt_rest(void) {

  FILTERCAL_DEF *FILTERCAL_OBS  = &CALIB_INFO.FILTERCAL_OBS ;
  FILTERCAL_DEF *FILTERCAL_REST = &CALIB_INFO.FILTERCAL_REST ;
  int NFILTDEF_OBS = FILTERCAL_OBS->NFILTDEF ;
  int OPT_REFAC    = MASK_FRAME_OBS + 8;
  int OPT_LEGACY   = OPT_REFAC ;
  double z = 0.24;

  int ifiltdef, ifilt, ifiltdef_rest, ifilt_r, irank ;
  double lamdif_min; 
  char *filter_name_obs, *filter_name_rest, *TEXT ;
  char fnam[] = "test_kcor_nearfilt_rest" ;

  // ----------- BEGIN ----------

  printf("\n %s with z=%.3f\n", fnam, z );

  printf("                                 nearest                     \n");
  printf("  rank  ifiltdef_obs:name   ifiltdef_rest:name    lamdif_min \n");
  
  for(irank=1; irank<=3; irank++ ) {
    for(ifilt=0; ifilt < NFILTDEF_OBS; ifilt++ ) {
      ifiltdef = FILTERCAL_OBS->IFILTDEF[ifilt] ;
      ifiltdef_rest = nearest_ifiltdef_rest(OPT_REFAC, ifiltdef, irank, z,
					    fnam, &lamdif_min );

      ifilt_r = FILTERCAL_REST->IFILTDEF_INV[ifiltdef_rest];
      filter_name_obs  = FILTERCAL_OBS->FILTER_NAME[ifilt];
      filter_name_rest = FILTERCAL_REST->FILTER_NAME[ifilt_r];

      printf(" rank=%d   %2d:%-12s       %2d:%-12s   %7.1f \n",  irank, 
	     ifiltdef, filter_name_obs,  ifiltdef_rest,filter_name_rest, lamdif_min);
      fflush(stdout);
    } 
  } // end irank_want

  return;

} // end test_kcor_nearfilt_rest


void test_kcor_table_avwarp(void) {

  int i, NTEST=0;
  int ifiltdef_a_list[10], ifiltdef_b_list[10];
  double C_list[10];
  double Trest = 4.5 ;
  double z = ZAT10PC;

  char fnam[] = "test_kcor_table_avwarp" ;

  // ------------ BEGIN -----------

  printf("\n %s: \n", fnam);

  ifiltdef_a_list[NTEST] = INTFILTER("B");
  ifiltdef_b_list[NTEST] = INTFILTER("V");
  C_list[NTEST] = -0.34 ;
  NTEST++ ;

  ifiltdef_a_list[NTEST] = INTFILTER("B");
  ifiltdef_b_list[NTEST] = INTFILTER("V");
  C_list[NTEST] = +0.34 ;
  NTEST++ ;

  ifiltdef_a_list[NTEST] = INTFILTER("R");
  ifiltdef_b_list[NTEST] = INTFILTER("I");
  C_list[NTEST] = -0.20 ;
  NTEST++ ;

  ifiltdef_a_list[NTEST] = INTFILTER("U");
  ifiltdef_b_list[NTEST] = INTFILTER("B");
  C_list[NTEST] = -5.4 ;
  NTEST++ ;

  printf("  filtdef_a      ifiltdef_b      Color     AVwarp    istat\n");
  FILTERCAL_DEF  *FILTERCAL_REST = &CALIB_INFO.FILTERCAL_REST ;
  char *TEXT;
  double AVWARP = 0.0, C, mag_a, mag_b ; ;
  int ifiltdef_a, ifiltdef_b, ifilt_a, ifilt_b, istat ;
  char *name_a, *name_b ;

  for(i=0; i < NTEST; i++ ) {

    ifiltdef_a = ifiltdef_a_list[i];
    ifiltdef_b = ifiltdef_b_list[i];
    C          = C_list[i];

    mag_a      = 0.0;
    mag_b      = mag_a - C;  // C = mag_a - mag_b
    ifilt_a    = FILTERCAL_REST->IFILTDEF_INV[ifiltdef_a];
    ifilt_b    = FILTERCAL_REST->IFILTDEF_INV[ifiltdef_b];
    name_a     = FILTERCAL_REST->FILTER_NAME[ifilt_a];
    name_b     = FILTERCAL_REST->FILTER_NAME[ifilt_b];

    AVWARP = eval_kcor_table_AVWARP(ifiltdef_a, ifiltdef_b, 
				    mag_a, mag_b, Trest, &istat);

    printf("  %2d:%-10s  %2d:%-10s  %7.3f, %8.4f  %2d \n",
	   ifiltdef_a, name_a,  ifiltdef_b, name_b, C, 
	   AVWARP, istat );
    fflush(stdout);
  }

  return ;

} // end test_kcor_table_avwarp


void test_kcor_table_kcor(void) {

  int ifiltdef_obs_list[20], ifiltdef_rest_list[20];
  int ifiltdef_o, ifiltdef_r, ifilt_o, ifilt_r ;
  double T_list[20], z_list[20], AV_list[20];
  double T, z, AV, kcor;
  int  i, NTEST, len_o, len_r ;
  char *TEXT, *name_o, *name_r, kcor_sym[20] ;
  char fnam[] = "test_kcor_table_kcor" ;

  // ------------- BEGIN -----------

  printf("\n %s: \n", fnam);

  ifiltdef_rest_list[NTEST] = INTFILTER("B");
  ifiltdef_obs_list[NTEST]  = INTFILTER("g");
  T_list[NTEST]  = 5.44 ;
  z_list[NTEST]  = 0.12 ;
  AV_list[NTEST] = 1.2; 
  NTEST++ ;

  ifiltdef_rest_list[NTEST] = INTFILTER("V");
  ifiltdef_obs_list[NTEST]  = INTFILTER("r");
  T_list[NTEST]  = 5.44 ;
  z_list[NTEST]  = 0.29 ;
  AV_list[NTEST] = -0.4 ; 
  NTEST++ ;

  printf("KCOR:  Trest      z       AV    kcorSym   kcorVal  \n");

  for(i=0; i < NTEST; i++ ) {
    ifiltdef_r = ifiltdef_rest_list[i];
    ifiltdef_o = ifiltdef_obs_list[i];
    T = T_list[i]; z= z_list[i];  AV=AV_list[i];

    ifilt_o = CALIB_INFO.FILTERCAL_OBS.IFILTDEF_INV[ifiltdef_o];
    ifilt_r = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef_r];
    name_o  = CALIB_INFO.FILTERCAL_OBS.FILTER_NAME[ifilt_o];
    name_r  = CALIB_INFO.FILTERCAL_REST.FILTER_NAME[ifilt_r];
    len_o   = strlen(name_o);
    len_r   = strlen(name_r);

    kcor = eval_kcor_table_KCOR(ifiltdef_r, ifiltdef_o, T, z, AV);

    sprintf(kcor_sym, "K_%c%c", name_r[len_r-1], name_o[len_o-1] );
    printf("KCOR:  %5.2f    %5.3f   %6.2f   %s   %7.4f \n",
	   T, z, AV, kcor_sym, kcor );
    fflush(stdout);
  }



  return;
} // end test_kcor_table_kcor



void test_GET_KCOR_DRIVER(void) {

  int IFILTDEF_OBS = INTFILTER("r");
  int IFILTDEF_REST_LIST[6];
  double Trest = -3.3 ;
  double z     = 0.20 ;
  double z0    = ZAT10PC ;
  int OPT_NEAREST   = MASK_FRAME_OBS + 8;

  double lamdif_min, LAMDIF[6], LAMAVG[6], MAG_REST_LIST[6], AVWARP_LIST[6] ;
  double kcor, lamavg, AVwarp = 0.1 ;
  int inear, ifiltdef_rest, ifilt_r;
  char *TEXT ;
  char fnam[] = "test_GET_KCOR_DRIVER" ;

  // ------------- BEGIN ----------

  printf("\n xxx %s: \n", fnam);

  for(inear=0; inear <3; inear++ ) {
    int rank = inear + 1;
    ifiltdef_rest =
      nearest_ifiltdef_rest(OPT_NEAREST, IFILTDEF_OBS, rank, z,
			    fnam, &lamdif_min );    
    IFILTDEF_REST_LIST[inear] = ifiltdef_rest ;
    MAG_REST_LIST[inear] = eval_kcor_table_LCMAG(ifiltdef_rest, Trest, z0, AVwarp) ;

    ifilt_r = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef_rest];
    lamavg  = CALIB_INFO.FILTERCAL_REST.LAMMEAN[ifilt_r];
    LAMAVG[inear] = lamavg;
    LAMDIF[inear] = lamavg - LAMAVG[0];

    printf("  inear=%d  ifiltdef_rest=%d  MAG=%.2f LAMAVG=%.1f LAMDIF=%.1f \n",
	   inear, ifiltdef_rest, MAG_REST_LIST[inear], lamavg, LAMDIF[inear]);
    fflush(stdout);

  }

  
  kcor = GET_KCOR_DRIVER(IFILTDEF_OBS, IFILTDEF_REST_LIST, MAG_REST_LIST,
			 LAMDIF, Trest, z, AVWARP_LIST );

  printf("   kcor=%.4f  AVwarp=%.4f \n",	
	 kcor, AVWARP_LIST[1] );

  return ;

} // end test_GET_KCOR_DRIVER

// =================================
#ifdef USE_KCOR_FORTRAN
void test_fortran(void) {

  char kcorFile[80] = "  " ;
  int IERR, ifilt_rest, ifilt_obs, iz;
  double t8, z8, av8, kcor8 ;

  // ----------- BEGIN -----------

  print_banner ( " TEST FORTRAN INTERFACE " );

  init_snvar__(&IERR);
  sprintf(kcorFile, "%s", INPUTS.KCOR_FILE );
  rdkcor_(kcorFile, &IERR, strlen(kcorFile) );

  av8 = 0.0;
  z8  = 0.20;
  t8  = 0.0;

  ifilt_rest = 7 ;
  ifilt_obs  = 3 ;


  for ( iz=0; iz<=10; iz++ ) {
    z8 = 0.04 * (double)iz ;
    printf(" --------------------------------------- \n");
    kcor8  = get_kcor8__( &ifilt_rest, &ifilt_obs, &t8, &z8, &av8 );
    printf(" xxx z=%4.2f  Tr=%5.2f AV=%4.2f   K_(%d->%d) = %f \n",
	   z8, t8, av8, ifilt_rest, ifilt_obs, kcor8 );
  }

  exit(1);
  
}  // end of test_fortran
#endif

// ******************************
void test_zcmb_dLmag_invert(void) {

  char fnam[] = "test_zcmb_dLmag_invert" ;
  double MU, zCMB=0.0 ;
  // ----------- RETURN ------------
  for(MU=32.0; MU < 49.0; MU+=1.0 ) {
    zCMB = zcmb_dLmag_invert(MU, &INPUTS.HzFUN_INFO, &INPUTS.ANISOTROPY_INFO );
  }
  debugexit(fnam);
  return ;
} // end test_zcmb_dLmag_invert

// ******************************
void test_PARSE_WORDS(void) {

  char fnam[] = "test_PARSE_WORDS" ;
  char dumString[] = "I walk to the beach";
  char dumFile[]   = "SIMGEN_SALT2_SDSS_FITS.input" ;
  //char dumFile[]   = "DES_SVA2_Y2.HOSTLIB" ;
  char tmpWord[80];
  int NWD, iwd;

  store_PARSE_WORDS(-1,"", fnam);

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE,dumFile, fnam);

  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(0,iwd,tmpWord);
    printf(" xxx word[%2d] = '%s' \n", iwd, tmpWord);
  }

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,dumString, fnam);
  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(0,iwd,tmpWord);
    printf(" xxx word[%2d] = '%s' \n", iwd, tmpWord);
  }

  debugexit(fnam);

  return ;
} // end test_PARSE_WORDS


// *********************
void test_ran(void) {
  int i, NTMP=0 ;
  double x, y; 
  double x0 = 0.0 ;
  // ------------- BEGIN --------
  for ( i = 1; i<=100000; i++ ) {
    // now smear with sigma=1
    NTMP++ ;
    if ( NTMP==100 ) {  fill_RANLISTs();  NTMP = 0 ; }
    x = getRan_Gauss(1);
    if ( x >=  0.0 ) 
      { y = 0.5 + GaussIntegral(x0,x);     }
    else 
      { y = 0.5 - GaussIntegral(x0,fabs(x) );     }
    printf(" 7777777 %f %f \n", x, y );
  }
  exit(1);
}  // end of test_ran

// ***********************
void test_igm(void) {

  int ilam;
  double lrest, lobs, z, z_a, z_b, tau_a, tau_b;
  char PATH_IGM_PARAM[MXPATHLEN];
  //  char fnam[] = "test_igm";

  // --------------- BEGIN --------------

  print_banner ( " Initialize igm correction for Spectra " );

  sprintf(PATH_IGM_PARAM,"%s/models/inoue_igm", PATH_SNDATA_ROOT);
  sprintf(DLA_FILE,"%s/DLAcoeff.txt", PATH_IGM_PARAM);
  sprintf(LAF_FILE,"%s/LAFcoeff.txt", PATH_IGM_PARAM);

  read_Inoue_coeffs();

  z_a = 2.0;
  z_b = 5.0;

  for (ilam=500; ilam <=6000; ilam+=500 ) {
    lrest = (double)ilam ;
    lobs  = lrest*(1.+z);

    z     = z_a ;
    tau_a = tLSLAF(z,lobs) + tLCLAF(z,lobs) + tLSDLA(z,lobs) + tLCDLA(z,lobs); 

    z     = z_b ;
    tau_b = tLSLAF(z,lobs) + tLCLAF(z,lobs) + tLSDLA(z,lobs) + tLCDLA(z,lobs); 

    printf(" xxx lrest=%6.0f tau(z=%.1f)=%le   tau(z=%.1f)=%le\n", 
	   lrest, z_a, tau_a,  z_b, tau_b) ;
    fflush(stdout);
  }
  
  exit(1);

} // end test_igm 


// ===============================
void test_getRan_funVal(char *FUNVAL_NAME) {

  // Created Dec 2023
  // Test asymGauss and expHalf by generating randoms from
  // distribution and comparing random distribution to 
  // funVal-computed distribution

  double range[2], xval[10];
  char   parName[40], parName2[40] ;
  int    IDMAP ;
  bool IS_LOGPARAM;
  GENGAUSS_ASYM_DEF     GENGAUSS ;
  GEN_EXP_HALFGAUSS_DEF GENEXP   ;

  char fnam[] = "test_getRan_funVal" ;

  // ----------- BEGIN ---------

  print_banner(fnam);

  GENGAUSS.USE = GENEXP.USE = false;

  if ( strcmp(FUNVAL_NAME,STRING_FUNVAL_GENGAUSS) == 0 ) { 
    load_test_GENGAUSS(&GENGAUSS );
    range[0]  = GENGAUSS.RANGE[0] ;
    range[1]  = GENGAUSS.RANGE[1] ;
  }
  else if ( strcmp(FUNVAL_NAME,STRING_FUNVAL_GENEXP) == 0 ) {
    load_test_GENEXP(&GENEXP);
    range[0]  = GENEXP.RANGE[0] ;
    range[1]  = GENEXP.RANGE[1] ;
  }
  else if ( strcmp(FUNVAL_NAME,STRING_FUNVAL_GENPDF) == 0 ) {
    int  GENPDF_OPTMASK   = 1 ;
    char GENPDF_FILE[200] = 
      "$DES5YR/populations/FINAL_forDES5yr/DES5YR_S3W22N21_GENPDF.DAT"; 
    ENVreplace(GENPDF_FILE, fnam, 1);
    init_genPDF(GENPDF_OPTMASK, NULL, GENPDF_FILE, BLANK_STRING) ;

    // sprintf(parName,"SALT2x1"); 
    // sprintf(parName,"SALT2c"); 
    sprintf(parName,"EBV"); 
    SNHOSTGAL.IGAL = 20;   
    IDMAP          = IDMAP_GENPDF(parName, &IS_LOGPARAM) ;
    range[0]       = GENPDF[IDMAP].GRIDMAP.VALMIN[0] ;
    range[1]       = GENPDF[IDMAP].GRIDMAP.VALMAX[0] ; 

    // restrict range to make better use of binned results below 
    // SALT2x1 range is ok
    if ( strcmp(parName,"SALT2c") == 0 ) { range[0] = -0.3; range[1] = 0.3; }
    if ( strcmp(parName,"EBV")    == 0 ) { range[0] =  0.0; range[1] = 0.9; }

  }
  else {
    print_preAbort_banner(fnam);
    printf("  Valid FUNVAL names for UNIT_TEST argument: \n");
    printf("\t %s\n", STRING_FUNVAL_GENGAUSS);
    printf("\t %s\n", STRING_FUNVAL_GENEXP);
    printf("\t %s\n", STRING_FUNVAL_GENPDF);

    sprintf(c1err,"Undefined FUNVAL_NAME = '%s' ", FUNVAL_NAME);
    sprintf(c2err,"Check valid FUNVAL unit test names above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  printf("\n - - -  Done with UNIT_TEST init - - - - \n");

  printf("\t Range is %.3f to %.3f \n", 
	 range[0], range[1]); fflush(stdout);

  if ( NMAP_GENPDF > 0 ) {
    printf("\t Select GENPDF map for %s (IDMAP=%d) \n",
	   parName, IDMAP);

  }
  // ---------------------------------
  int    NRANGEN = 500000;
  int    NBIN_HIST = 600; //600;
  double range_tot = range[1] - range[0];
  double BINSIZE_HIST = range_tot / (double)NBIN_HIST ;
  double r, x, xmin, xmax, funVal, genVal, genVal_err, nsig;
  double funVal_edge[2], funVal_integral=0.0 ;

  int MEMD = NBIN_HIST * sizeof(double);
  int MEMI = NBIN_HIST * sizeof(int);
  double *funVal_PER_BIN = (double *) malloc(MEMD) ;
  double *x_PER_BIN      = (double *) malloc(MEMD) ;
  int    *NGEN_PER_BIN   = (int    *) malloc(MEMI) ;
  int i, bin, MX_PER_BIN=0, N2SIG=0, N3SIG=0, NBIN_NONZERO=0 ;

  printf("\t Generate %d random numbers from %s:\n",
	 NRANGEN, FUNVAL_NAME);

  fflush(stdout);

  for(bin=0; bin < NBIN_HIST; bin++ ) { 
    NGEN_PER_BIN[bin] = 0 ; 
    x_PER_BIN[bin] = range[0] + ((double)bin + 0.5) * BINSIZE_HIST; // bin center
  }

  int NWRAP, NWRAP_LAST = -9;

  for(i=0; i < NRANGEN; i++ ) {

    NWRAP = GENRAN_INFO.NWRAP[1] ;
    if ( NWRAP != NWRAP_LAST ) { printf("\t xxx increment NWRAP=%d\n", NWRAP); }
    NWRAP_LAST = NWRAP ;

    if ( i % 200 == 0 ) { fill_RANLISTs(); }

    if ( GENGAUSS.USE ) 
      { r   = getRan_GENGAUSS_ASYM(&GENGAUSS) ; }
    else if ( GENEXP.USE ) 
      { r   = getRan_GEN_EXP_HALFGAUSS(&GENEXP) ; }
    else if ( NMAP_GENPDF > 0 ) 
      { r = getRan_genPDF(parName, &GENGAUSS); }

    if ( r < range[0] || r > range[1] ) { continue; }

    bin = (int)( (r - range[0])/BINSIZE_HIST );
    NGEN_PER_BIN[bin]++;
    if ( i < -20 ) {
      printf(" xxx i=%2d r=%f bin=%d \n", i, r, bin);
    }
  }

  // - - - - - 
  // get funVal and integral 
  for(bin=0; bin < NBIN_HIST; bin++ ) {
    x      = x_PER_BIN[bin];
    xmin   = x - 0.5*BINSIZE_HIST ;
    xmax   = x + 0.5*BINSIZE_HIST ;

    if ( GENGAUSS.USE ) {
      funVal_edge[0] = funVal_GENGAUSS_ASYM(xmin, &GENGAUSS) ;
      funVal_edge[1] = funVal_GENGAUSS_ASYM(xmax, &GENGAUSS) ;
    }
    else if ( GENEXP.USE ) {
      funVal_edge[0] = funVal_GEN_EXP_HALFGAUSS(xmin, &GENEXP) ;
      funVal_edge[1] = funVal_GEN_EXP_HALFGAUSS(xmax, &GENEXP) ;
    }
    else if ( NMAP_GENPDF > 0 ) {
      xval[0] = xmin; 
      funVal_edge[0] = funVal_genPDF(parName, xmin, &GENGAUSS) ;
      xval[0] = xmax; 
      funVal_edge[1] = funVal_genPDF(parName, xmax, &GENGAUSS) ;
    }

    funVal         = 0.5 * ( funVal_edge[0] + funVal_edge[1] );
    funVal_integral    += funVal ;
    funVal_PER_BIN[bin] = funVal ;
  }

  // - - - - - 
  double XNORM = funVal_integral / (double)NRANGEN;
  printf("\t XNORM(funVal_integral / %d) = %le \n", NRANGEN, XNORM);

  printf("\n");
  printf("#    x      funVal      genVal  \n");
  for(bin=0; bin < NBIN_HIST; bin++ ) {
    x      = x_PER_BIN[bin];
    funVal = funVal_PER_BIN[bin];
    genVal = (double)NGEN_PER_BIN[bin] * XNORM ;

    genVal_err = nsig = 0.0 ;
    bool USE = ( funVal > 0.01 );
    if ( USE ) {      
      genVal_err = sqrt(funVal/XNORM) * XNORM; // use model to get error
      nsig       = fabs(genVal-funVal) / genVal_err ;
      NBIN_NONZERO++ ;
      if ( nsig > 2.0 ) { N2SIG++ ; }
      if ( nsig > 3.0 ) { N3SIG++ ; }
    }

    if ( USE ) {
      printf("   %6.3f   %8.5f  %8.5f +_ %8.5f  (%4.2f sigma) \n", 
	     x, funVal, genVal, genVal_err, nsig );
    }
    fflush(stdout);
    // N_PER_BIN[bin] = 0 ;
  }
  
  // - - - - -
  double frac2, frac3;
  printf("\n");
  frac2 = (double)N2SIG / (double)(NBIN_NONZERO);
  frac3 = (double)N3SIG / (double)(NBIN_NONZERO);
  printf(" Number of 2-sigma outliers: %d of %d  (frac = %.3f)\n",
	 N2SIG, NBIN_NONZERO, frac2 );
  printf(" Number of 3-sigma outliers: %d of %d  (frac = %.3f)\n",
	 N3SIG, NBIN_NONZERO, frac3 );

  fflush(stdout);

  // summarize variables for GENPDF
  if ( NMAP_GENPDF ) {
    int IVAR_HOSTLIB, ivar;
    int IGAL = SNHOSTGAL.IGAL ;
    long long GALID = get_GALID_HOSTLIB(IGAL);;
    char *tmpName;
    printf("\n    GENPDF varname: %s  from MAP=%s  (IDMAP=%d)\n", 
	   parName, GENPDF[IDMAP].MAPNAME, IDMAP);
    printf("\t GALID = %lld   (IGAL = %d)\n", GALID, IGAL);
    for(ivar=1; ivar < GENPDF[IDMAP].GRIDMAP.NDIM ; ivar++ ) {
      IVAR_HOSTLIB   = GENPDF[IDMAP].IVAR_HOSTLIB[ivar];
      x              = get_VALUE_HOSTLIB(IVAR_HOSTLIB, IGAL);
      tmpName        = GENPDF[IDMAP].VARNAMES[ivar];
      printf("\t %-12s = %f \n", tmpName, x);
      fflush(stdout);
    }
  }

  debugexit(fnam);

  return;
} // end test_getRan_funVal


// =========
void load_test_GENGAUSS(GENGAUSS_ASYM_DEF *GENGAUSS ) {

  // Created Dec 2023:
  // hard-wire test values into GENGAUSS structure.

  // --------- BEGIN ----------

  init_GENGAUSS_ASYM(GENGAUSS, (double)0.0 ); 
  GENGAUSS->USE      = true;
  sprintf(GENGAUSS->NAME,"UNIT_TEST");

  GENGAUSS->PEAK     =  0.0;
  GENGAUSS->SIGMA[0] =  0.2; 
  GENGAUSS->SIGMA[1] =  0.4 ;
  GENGAUSS->RANGE[0] = -1.0 ;
  GENGAUSS->RANGE[1] =  3.0 ;  // check 1.4 and 3

  GENGAUSS->PROB2     = 0.3 ;
  GENGAUSS->PEAK2     = 1.0 ; 
  GENGAUSS->SIGMA2[0] = 0.1 ;
  GENGAUSS->SIGMA2[1] = 0.3 ;

  double prob_expon_rewgt         = 0.5 ;
  GENGAUSS->PROB_EXPON_REWGT      = prob_expon_rewgt;
  GENGAUSS->SQRT_PROB_EXPON_REWGT = sqrt(prob_expon_rewgt) ;

  dump_GENGAUSS_ASYM(GENGAUSS);

  return ;

} // end load_test_GENGAUSS


void load_test_GENEXP(GEN_EXP_HALFGAUSS_DEF *GENEXP) {

  double tau = 0.3;
  double range[2] = { 0.0, 2.0 } ;

  // --------- BEGIN ----------

  init_GEN_EXP_HALFGAUSS(GENEXP, (double)0.0 );

  sprintf(GENEXP->NAME,"UNIT_TEST");  
  set_GEN_EXPON(tau, range, GENEXP);

  // define half-Gaussian core
  bool DO_GAUSS = true ;
  if ( DO_GAUSS ) {
    GENEXP->RATIO = 0.3 ;  
    GENEXP->PEAK  = 0.0 ;
    GENEXP->SIGMA = 0.05 ;
  }

  double prob_expon_rewgt       = 0.5 ;
  GENEXP->PROB_EXPON_REWGT      = prob_expon_rewgt ;
  GENEXP->SQRT_PROB_EXPON_REWGT = sqrt(prob_expon_rewgt) ;

  dump_GEN_EXP_HALFGAUSS(GENEXP);

  return ;

} // end load_test_GENEXP
