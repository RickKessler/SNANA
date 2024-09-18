/******************
 May 2017
 Implement weak lensing effects on distance modulus
 using map defined by sim-input 
   WEAKLENS_PROBMAP_FILE:  <file with: z deltaMU prob>

 HOSTLIB-WEAKLENS_DMU is treated elsewhere (not here).
 The map in WEAKLENS_PROBMAP_FILE overrides optional 
 HOSTLIB-WEAKLENS_DMU.

 Oct 2021: require DOCANA at top of weak lensing map file.

 ******************/

#include "sntools.h"
#include "sntools_weaklens.h"

// =============================================
void init_lensDMU(char *mapFileName, float dsigma_dz) {

  // Created May 1 2017
  // Read lensing map of Prob(z,dMu) from input *mapFileName.
  // Compute inverse map, DMU(z,Prob) for fast selection of
  // random DMU from random number.
  //
  // If dsigma_dz > 0 then just print message that model
  // is symmetric Gaussian.
  //
  // If file does not exist, return but don't abort.
  //
  // Oct 10 2021: require DOCANA
  // Apr 12 2024: apply bound checks on NZ and NDUM (number of bins)
  
  FILE *FPMAP;
  int    OPENMASK = OPENMASK_VERBOSE + OPENMASK_REQUIRE_DOCANA ;
  int    MEMD, NROW_TOT, NROW_TABLE, irow;
  int    jj, iz, imu, gzipFlag, NZ=0, NDMU=0 ;  
  bool   DOCANA_END = false;
  double Prob, dmu, ztmp, z=0.0, zLAST=-9.9 ;
  double SUM_WGT, SUM_dmu, dmu_avg, SUM_Prob ;

  char PATH_DEFAULT[2*MXPATHLEN], MAPFILENAME[MXPATHLEN];
  char tmpLine[MXPATHLEN], tmpLine_copy[MXPATHLEN];
  int  LEN, iwd, NWD ;
  char *ptrWORD[10]; int MXWD=10;
  char *firstWord, sep[] = " ";
  char fnam[] = "init_lensDMU" ;

  // ---------------- BEGIN -----------------

  LENSING_PROBMAP.USEFLAG = 0 ;

  if ( dsigma_dz > 1.0E-8 ) {
    print_banner(fnam);
    printf("\t Symmetric Gaussian model: sigma = %.3f * z \n", 
	   dsigma_dz );
    return ;
  }

  // check option to ignore lensing if mapFileName = 'NULL', 'NONE', 'BLANK'
  if ( IGNOREFILE(mapFileName) ) { return; }

  print_banner(fnam);

  sprintf(PATH_DEFAULT,"%s %s/models/lensing", 
	  PATH_USER_INPUT, PATH_SNDATA_ROOT);
  
  // check mapFileName in user dir, then under PATH_DEFAULT
  FPMAP = snana_openTextFile(OPENMASK, PATH_DEFAULT, mapFileName, 
			     MAPFILENAME, &gzipFlag );

  if ( FPMAP == NULL ) {
    abort_openTextFile("WEAKLENS_PROBMAP_FILE", 
		       PATH_DEFAULT, mapFileName, fnam);
  }


  // ---------------------------------
  LENSING_PROBMAP.USEFLAG = 1 ;
  // ---------------------------------

  for(iwd=0; iwd < MXWD; iwd++ ) 
    { ptrWORD[iwd] = (char*)malloc(200*sizeof(char)) ; }

  // we have a mapFile
  // get NROW for malloc

  MEMD = sizeof(double);
  NROW_TOT = nrow_read(MAPFILENAME,fnam); // includes DOCANA rows
  printf("\t Found %d total rows in file. \n", NROW_TOT );

  // malloc temp 1D arrays to read everything
  double *z_TMP1D    = (double*) malloc( MEMD * NROW_TOT );
  double *dmu_TMP1D  = (double*) malloc( MEMD * NROW_TOT );
  double *Prob_TMP1D = (double*) malloc( MEMD * NROW_TOT );

  NROW_TABLE = 0;
  for(irow=0; irow < NROW_TOT; irow++ ) {

    fgets(tmpLine, MXPATHLEN, FPMAP);
    if ( commentchar(tmpLine) ) { continue ; }

    if ( !DOCANA_END ) { 
      sprintf(tmpLine_copy,"%s", tmpLine);
      LEN = strlen(tmpLine_copy);
      firstWord = strtok(tmpLine_copy, sep);
      if ( LEN > 0 ) { firstWord[LEN-1] = 0; } // remove termination
      if ( strcmp(firstWord,KEYNAME2_DOCANA_REQUIRED)==0 ) 
	{ DOCANA_END = true; }
      continue ;
    }

    splitString(tmpLine, " ", fnam, MXWD, &NWD, ptrWORD);
    if ( NWD != 3 ) { continue; }

    /* xxxxxx
    if ( NROW_TABLE == 0 ) {
      printf(" xxx %s: first table line: '%s' \n", fnam, tmpLine);
    }
    xxx */

    // we have a line after DOCANA that is not a comment line; parse it
    sscanf(ptrWORD[0], "%le", &z_TMP1D[NROW_TABLE] );
    sscanf(ptrWORD[1], "%le", &dmu_TMP1D[NROW_TABLE] );
    sscanf(ptrWORD[2], "%le", &Prob_TMP1D[NROW_TABLE] );
    
    z = z_TMP1D[NROW_TABLE] ;
    if ( z > zLAST ) { NZ++ ; }
    zLAST = z;
    NROW_TABLE++ ;
  }
  fclose(FPMAP);

  // - - - - - 
  if ( NZ == 0 ) {
    sprintf(c1err,"Number of weak-lensing z bins is zero.");
    sprintf(c2err,"Something is really messed up.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }
  
  // convert TMP1D arrays to 2D map 
  NDMU = NROW_TABLE/NZ ;

  printf("\t Convert %d lens-table rows into %d(z) x %d(DMU) array\n",
	 NROW_TABLE, NZ, NDMU); fflush(stdout);

  if ( NZ >= MXBIN_LENSING_z ) {
    sprintf(c1err,"NZ(bins) = %d exceeds bound of MXBIN_LENSING_z = %d", NZ, MXBIN_LENSING_z);
    sprintf(c2err,"Either reduce NZ or increase MXBIN_LENSING_z");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );        
  }

  if ( NDMU >= MXBIN_LENSING_dmu ) {
    sprintf(c1err,"NDMU(bins) = %d exceeds bound of MXBIN_LENSING_dmu = %d",
	    NDMU, MXBIN_LENSING_dmu);
    sprintf(c2err,"Either reduce NDMU or increase MXBIN_LENSING_dmu");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );        
  }

  // allocate memory for LENSING structure
  LENSING_PROBMAP.NBIN_z   = NZ;
  LENSING_PROBMAP.NBIN_dmu = NDMU;
  LENSING_PROBMAP.z_LIST     = (double*) malloc(MEMD * NZ );
  LENSING_PROBMAP.dmu_LIST   = (double*) malloc(MEMD * NDMU );

  LENSING_PROBMAP.PROB       = (double**) malloc(sizeof(double*) * NZ );
  LENSING_PROBMAP.FUNPROB    = (double**) malloc(sizeof(double*) * NZ );
  LENSING_PROBMAP.FUNDMU     = (double**) malloc(sizeof(double*) * NZ );
  for(iz=0; iz<NZ; iz++ ) { 
    LENSING_PROBMAP.PROB[iz]    = (double*)malloc(MEMD*NDMU) ; 
    LENSING_PROBMAP.FUNPROB[iz] = (double*)malloc(MEMD*NDMU) ; 
    LENSING_PROBMAP.FUNDMU[iz]  = (double*)malloc(MEMD*NDMU) ; 
  }


  // transfer 1D contents to 2D PROBMAP
  irow=0;
  for(iz=0; iz < NZ; iz++ ) {
    SUM_WGT = SUM_Prob = 0.0 ;
    for(imu=0; imu < NDMU; imu++ ) {
      z     = z_TMP1D[irow] ;
      dmu   = dmu_TMP1D[irow] ;
      Prob  = Prob_TMP1D[irow] ;
      SUM_WGT  += (Prob*dmu) ;
      SUM_Prob += (Prob);

      LENSING_PROBMAP.z_LIST[iz]    = z;
      LENSING_PROBMAP.dmu_LIST[imu] = dmu;
      LENSING_PROBMAP.PROB[iz][imu] = Prob ;

      LENSING_PROBMAP.FUNPROB[iz][imu] = SUM_Prob ;
      LENSING_PROBMAP.FUNDMU[iz][imu]  = dmu ;

      irow++ ;
    } // end imu

    // normalize SUM_prob to have max=1
    SUM_Prob = LENSING_PROBMAP.FUNPROB[iz][NDMU-1];
    LENSING_PROBMAP.FUNPROB[iz][0] = 0.0 ; // set small prob to 0
    for(imu=0; imu < NDMU; imu++ ) 
      { LENSING_PROBMAP.FUNPROB[iz][imu] /= SUM_Prob; }

    // compute <dmu> to check how close to zero
    SUM_dmu = LENSING_PROBMAP.dmu_LIST[NDMU-1] - LENSING_PROBMAP.dmu_LIST[0] ;
    dmu_avg = SUM_WGT/SUM_dmu ;

    for(jj=1; jj<=4; jj++ ) {
      ztmp = (double)jj;
      if ( fabs(z-ztmp) < .02 ) 
	{ printf("\t\t <dmu>[z=%5.3f] = %12.4le \n", z, dmu_avg); }
    }

  } // end iz

  // free TMP1D arrays
  free(z_TMP1D);    free(dmu_TMP1D);    free(Prob_TMP1D);

  // store min/max redshift and dmu
  LENSING_PROBMAP.zMIN   =  LENSING_PROBMAP.z_LIST[0];
  LENSING_PROBMAP.zMAX   =  LENSING_PROBMAP.z_LIST[NZ-1];
  LENSING_PROBMAP.dmuMIN =  LENSING_PROBMAP.dmu_LIST[0] ;
  LENSING_PROBMAP.dmuMAX =  LENSING_PROBMAP.dmu_LIST[NDMU-1] ;

  // print summary 
  printf("\t Done initializing Prob(%.3f < DeltaMU < %.3f) \n",
	 LENSING_PROBMAP.dmu_LIST[0], LENSING_PROBMAP.dmu_LIST[NDMU-1] );
  printf("\t in %d redshfit bins from %.3f to %.3f \n",
	 NZ, LENSING_PROBMAP.z_LIST[0], LENSING_PROBMAP.z_LIST[NZ-1] );
  fflush(stdout);

  for(iwd=0; iwd < MXWD; iwd++ )  { free(ptrWORD[iwd]); }

  return ;

} // end init_lensDMU


// ============================================
double gen_lensDMU( double z, double ran1, int DUMP_FLAG ) {
  // Created Apr 2017 
  // Return DMU from random lensing magnification.
  // DMU = -2.5*log10(magnification)
  //
  // Inputs:
  //  z = redshift
  //  ran1 = random number from 0 to 1
  //
  // Feb 16 2022; pass DUMP_FLAG arg

  double lensDMU = 0.0 ;  // default  
  double zMIN = LENSING_PROBMAP.zMIN ;
  double zMAX = LENSING_PROBMAP.zMAX ;
  int    NBIN_dmu = LENSING_PROBMAP.NBIN_dmu ;
  int    NBIN_z   = LENSING_PROBMAP.NBIN_z ;

  int OPT_INTERP = OPT_INTERP_LINEAR ;
  double ztmp, z0, z1, zScale = 0.0 ;
  double DMU0=-9.0, DMU1=-9.0, zFac=-9.0 ;
  int    iz, iz0, iz1;
  char fnam[] = "gen_lensDMU" ;

  // -------------- BEGIN ------------- 

  if ( LENSING_PROBMAP.USEFLAG == 0 ) { return(lensDMU) ;}

  //  printf(" xxx %s for z=%.4f and ran1=%f \n", fnam, z, ran1);

  iz0 = iz1 = -9;

  if ( z <= zMIN ) {
    iz0 = iz1 = 0 ;
    zScale = z/zMIN ;
  }
  else if ( z >= zMAX ) {
    iz0 = iz1 = NBIN_z - 1 ;
    zScale = z/zMAX ;
  }
  else {
    zScale = -9.0 ;
    for(iz=0; iz < NBIN_z-1; iz++ ) {
      ztmp = LENSING_PROBMAP.z_LIST[iz] ;
      if ( ztmp < z ) { iz0=iz; iz1=iz0+1; }
    }
  }


  if ( iz0 < 0 || iz1 < 0 ) {
    sprintf(c1err,"Could not find z-bin for z=%f", z);
    sprintf(c2err,"Check lensing DMU-map");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  z0 = LENSING_PROBMAP.z_LIST[iz0] ;
  z1 = LENSING_PROBMAP.z_LIST[iz1] ;

  double *ptrPROB[2], *ptrDMU[2] ;
  ptrPROB[0] = LENSING_PROBMAP.FUNPROB[iz0] ;
  ptrPROB[1] = LENSING_PROBMAP.FUNPROB[iz1] ;
  ptrDMU[0]  = LENSING_PROBMAP.FUNDMU[iz0] ;
  ptrDMU[1]  = LENSING_PROBMAP.FUNDMU[iz1] ;

  // - - - - - - - - - - - - - - - - - - - 
  if ( zScale > 0.0  ) {
    // scale lensDMU by zScale
    lensDMU = interp_1DFUN(OPT_INTERP, ran1, NBIN_dmu,
			   ptrPROB[0], ptrDMU[0], fnam );
    lensDMU *= zScale ;
  }
  else {
    // interpolate lensDMU between two z bins
    DMU0 = interp_1DFUN(OPT_INTERP, ran1, NBIN_dmu,
			ptrPROB[0], ptrDMU[0], fnam );
    DMU1 = interp_1DFUN(OPT_INTERP, ran1, NBIN_dmu,
			ptrPROB[1], ptrDMU[1], fnam );

    zFac    = (z - z0)/(z1-z0);
    lensDMU = DMU0 + (DMU1-DMU0)*zFac ;
  }

  if ( DUMP_FLAG ) {
    printf(" xxx ------------------------------ \n");
    printf(" xxx %s: dump for z=%f and ran1=%f \n", 
	   fnam, z, ran1);
    printf(" xxx %s: NBIN_dmu=%d  zMIN,zMAX = %f %f \n", 
	   fnam, NBIN_dmu, zMIN, zMAX);
    printf(" xxx %s: zScale = %f \n", 
	   fnam, zScale);
    printf(" xxx %s: z=%f z0=%f z1=%f zFac=%f \n", 
	   fnam, z, z0, z1, zFac);
    printf(" xxx %s: DMU0=%f  DMU1=%f  lensDMU = %f \n", 
	   fnam, DMU0, DMU1, lensDMU );
    fflush(stdout);
  }

  return(lensDMU) ; 

} // end gen_magLens 

