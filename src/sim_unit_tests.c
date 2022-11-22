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

void test_kcor_utils(void);
void test_kcor_table_lcmag(void);
void test_kcor_table_mwxt(void);
void test_kcor_table_avwarp(void);
void test_kcor_table_kcor(void);
void test_kcor_nearfilt_rest(void);

void test_fortran(void);
void test_igm(void);
void test_ran(void);
void test_PARSE_WORDS(void);
void test_zcmb_dLmag_invert(void);

char TEST_REFAC[]  = "REFAC";
char TEST_LEGACY[] = "LEGACY";

// ************************

void test_kcor_utils(void) {
  
  // Created Nov 2022

  char fnam[] = "test_kcor_utils" ;

  // --------------- BEGIN ------------

  // test_kcor_table_lcmag();
  // test_kcor_table_mwxt();
  // test_kcor_nearfilt_rest();

  test_kcor_table_avwarp();
  test_kcor_table_kcor();

  //.xyz
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

      if ( INPUTS.USE_KCOR_REFACTOR == 2 ) {
	MAG = eval_kcor_table_LCMAG(ifiltdef_r, Trest, z, AVwarp) ;
	text = TEST_REFAC ;
      }
      else {
	MAG   = get_maglc8__(&ifiltdef_r, &Trest, &z, &AVwarp );
	text  = TEST_LEGACY ;
      }
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

      if ( INPUTS.USE_KCOR_REFACTOR == 2 ) {
	MWXT = eval_kcor_table_MWXT(ifiltdef_o, Trest, z, AVwarp, 
			     MWEBV, RV, OPT_MWCOLORLAW) ;
	text = TEST_REFAC ;
      }
      else {
	MWXT  = get_mwxt8__(&ifiltdef_o, &Trest, &z, &AVwarp, 
			    &MWEBV, &RV, &OPT_MWCOLORLAW );
	text  = TEST_LEGACY ;
      }
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

  if (INPUTS.USE_KCOR_REFACTOR==2 ) 
    { TEXT= TEST_REFAC; }
  else
    { TEXT = TEST_LEGACY; }

  printf("\n %s with z=%.3f (%s)\n", fnam, z, TEXT );

  printf("                                 nearest                     \n");
  printf("  rank  ifiltdef_obs:name   ifiltdef_rest:name    lamdif_min \n");
  
  for(irank=1; irank<=3; irank++ ) {
    for(ifilt=0; ifilt < NFILTDEF_OBS; ifilt++ ) {
      ifiltdef = FILTERCAL_OBS->IFILTDEF[ifilt] ;

      if ( INPUTS.USE_KCOR_REFACTOR == 2 ) {
	ifiltdef_rest = nearest_ifiltdef_rest(OPT_REFAC, ifiltdef, irank, z,
					      fnam, &lamdif_min );
      }
      else {
	float zf = (float)z ;
	float  lamdif_min_f ;
	ifiltdef_rest = nearest_ifilt_rest__(&OPT_LEGACY, &ifiltdef, &irank, &zf, 
					     &lamdif_min_f);
	lamdif_min = (double)lamdif_min_f ;
      }

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

    if ( INPUTS.USE_KCOR_REFACTOR == 2 ) {
      //  AVWARP = fit_AVWARP(ifiltdef_a, ifiltdef_b, Trest, C);
      AVWARP = eval_kcor_table_AVWARP(ifiltdef_a, ifiltdef_b, mag_a, mag_b, Trest, &istat);
      TEXT = TEST_REFAC;
    }
    else {  
      AVWARP = get_avwarp8__(&Trest, &z, &mag_a, &mag_b, 
			     &ifiltdef_a, &ifiltdef_b, &istat);
      TEXT = TEST_LEGACY ;
    }

    printf("  %2d:%-10s  %2d:%-10s  %7.3f, %8.4f  %2d (%s)\n",
	   ifiltdef_a, name_a,  ifiltdef_b, name_b, C, 
	   AVWARP, istat, TEXT);
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

    if ( INPUTS.USE_KCOR_REFACTOR == 2 ) {
      kcor = eval_kcor_table_KCOR(ifiltdef_r, ifiltdef_o, T, z, AV);
      TEXT = TEST_REFAC ;
    }
    else {
      kcor = get_kcor8__(&ifiltdef_r, &ifiltdef_o, &T, &z, &AV);
      TEXT = TEST_LEGACY ;
    }

    sprintf(kcor_sym, "K_%c%c", name_r[len_r-1], name_o[len_o-1] );
    printf("KCOR:  %5.2f    %5.3f   %6.2f   %s   %7.4f  (%s) \n",
	   T, z, AV, kcor_sym, kcor, TEXT);
    fflush(stdout);
  }



  return;
} // end test_kcor_table_kcor

// =================================
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

// ******************************
void test_zcmb_dLmag_invert(void) {

  char fnam[] = "test_zcmb_dLmag_invert" ;
  double MU, zCMB=0.0 ;
  // ----------- RETURN ------------
  for(MU=32.0; MU < 49.0; MU+=1.0 ) {
    zCMB = zcmb_dLmag_invert(MU, &INPUTS.HzFUN_INFO);
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

  store_PARSE_WORDS(-1,"");

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_FILE,dumFile);

  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(0,iwd,tmpWord);
    printf(" xxx word[%2d] = '%s' \n", iwd, tmpWord);
  }

  NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,dumString);
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


