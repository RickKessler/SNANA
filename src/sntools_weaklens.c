/******************
 May 2017
 Implement weak lensing effects on distance modulus.

 ******************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

  FILE *FPMAP;

  int MEMD, NROW, irow, jj, iz, imu, gzipFlag, NZ=0, NDMU=0 ;  
  double VALTMP[3], Prob, dmu, ztmp, z=0.0, zLAST=-9.9 ;
  double SUM_WGT, SUM_dmu, dmu_avg, SUM_Prob ;

  char PATH_DEFAULT[2*MXPATHLEN], MAPFILENAME[MXPATHLEN] ;
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
  FPMAP = snana_openTextFile(1,PATH_DEFAULT, mapFileName, 
			     MAPFILENAME, &gzipFlag );

  if ( FPMAP == NULL ) {
    abort_openTextFile("WEAKLENS_PROBMAP_FILE", 
		       PATH_DEFAULT, mapFileName, fnam);

    /* xxxxxxxxx mark delete Feb 1 2020 xxxxxxxxxxx
    sprintf(c1err, "Could not find lensing mapFile");
    sprintf(c2err, " '%s' ", mapFileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
    xxxxxxxxxxxxx */
  }


  // ---------------------------------
  LENSING_PROBMAP.USEFLAG = 1 ;
  // ---------------------------------

  // we have a mapFile
  // get NROW for malloc

  MEMD = sizeof(double);
  NROW = nrow_read(MAPFILENAME,fnam);
  printf("\t Found %d rows to read. \n", NROW);

  // malloc temp 1D arrays to read everything
  double *z_TMP1D    = (double*) malloc( MEMD * NROW );
  double *dmu_TMP1D  = (double*) malloc( MEMD * NROW );
  double *Prob_TMP1D = (double*) malloc( MEMD * NROW );

  
  for(irow=0; irow < NROW; irow++ ) {
    readdouble(FPMAP, 3, VALTMP);
    z_TMP1D[irow]    = VALTMP[0] ;
    dmu_TMP1D[irow]  = VALTMP[1];
    Prob_TMP1D[irow] = VALTMP[2];

    z = z_TMP1D[irow] ;
    if ( z > zLAST ) { NZ++ ; }
    zLAST = z;
  }
  fclose(FPMAP);


  // convert TMP1D arrays to 2D map 
  NDMU = NROW/NZ ;

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
  printf("\t Done initializing Prpb(%.3f < DeltaMU < %.3f) \n",
	 LENSING_PROBMAP.dmu_LIST[0], LENSING_PROBMAP.dmu_LIST[NDMU-1] );
  printf("\t in %d redshfit bins from %.3f to %.3f \n",
	 NZ, LENSING_PROBMAP.z_LIST[0], LENSING_PROBMAP.z_LIST[NZ-1] );
  fflush(stdout);

  //  debugexit(fnam); // xxx REMOVE

  return ;

} // end init_lensDMU


// ============================================
double gen_lensDMU( double z, double ran1 ) {
  // Created Apr 2017 
  // Return DMU from random lensing magnification.
  // DMU = -2.5*log10(magnification)
  //
  // Inputs:
  //  z = redshift
  //  ran1 = random number from 0 to 1

  double lensDMU = 0.0 ;  // default  
  double zMIN = LENSING_PROBMAP.zMIN ;
  double zMAX = LENSING_PROBMAP.zMAX ;
  int    NBIN_dmu = LENSING_PROBMAP.NBIN_dmu ;
  int    NBIN_z   = LENSING_PROBMAP.NBIN_z ;

  int OPT_INTERP = OPT_INTERP_LINEAR ;
  double ztmp, z0, z1, zScale = 0.0 ;
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
    double DMU0, DMU1, zFac ;
    DMU0 = interp_1DFUN(OPT_INTERP, ran1, NBIN_dmu,
			ptrPROB[0], ptrDMU[0], fnam );
    DMU1 = interp_1DFUN(OPT_INTERP, ran1, NBIN_dmu,
			ptrPROB[1], ptrDMU[1], fnam );

    zFac    = (z - z0)/(z1-z0);
    lensDMU = DMU0 + (DMU1-DMU0)*zFac ;
  }

  return(lensDMU) ; 

} // end gen_magLens 

