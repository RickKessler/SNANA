/***********************************************************
 Created Dec 19 2015 by R.Kessler

 Model wrong host fraction and wrong redshift.
 This code is independent of the HOSTLIB code in sntools_host.c,
 but some links might be introduced later.

 Sim-input WRONGHOST_FILE specifies to things to read and store:
  1) wrongHost prob vs. ZTRUE(SN), defined as a 3rd order polynomial
     function of redshift

  2) list of "ZTRUE ZMATCH" when there is  mis-matched host.
     Code below finds a ZTRUE close to ZSN, and returns ZMATCH-ZTRUE
     to get a corrupted redshift.

 Apr 12 2019: 
   + refactor to allow comma-separated coefficients after the
     PROB_WRONGHOST_POLY key, and to also allow legacy input
     of 3rd-order poly with space-sep coefficients.

***********************************************************/

/*
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>
*/

#include "sntools.h"
#include "sntools_wronghost.h"

// =====================================

void INIT_WRONGHOST(char *inFile, double ZMIN, double ZMAX) {

  // read inFile and store WRONGHOST map between ZMIN and ZMAX
  //
  // Read the following
  //   PROB_WRONGHOST_POLY: [a0] [a1] [a2] [a3]
  //      or
  //   PROB_WRONGHOST_POLY: [a0],[a1],...  (arbitrary poly-order)
  //
  //   and then read list of ZTRUE,ZMATCH pairs after the
  //   comment lines.
  //
  //  Apr 2019: refactor to read comma-sep poly coefficients.

  FILE *fp;
  int  NROW_TOT, NROW_READ, MEMD, NLIST, isort, iz, *INDEX_SORT ;
  char c_get[60];
  char fnam[] = "INIT_WRONGHOST" ;

  // ---------- BEGIN ------------

  WRONGHOST.NLIST = 0 ;
  if ( IGNOREFILE(inFile) ) { return ; }

  WRONGHOST.USER_zMIN = ZMIN;
  WRONGHOST.USER_zMAX = ZMAX;

  if ( ZMAX < 1.0E-8 || ZMAX <= ZMIN ) {
    sprintf(c1err,"Invalid ZMIN=%f, ZMAX=%f", ZMIN, ZMAX);
    sprintf(c2err,"Check xxx" );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  fp = fopen(inFile, "rt") ;
  if ( !fp ) {
    sprintf(c1err,"Could not open WRONGHOST_FILE");
    sprintf(c2err,"%s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  sprintf(BANNER,"%s: Init wrongHost prob and zmatch-ztrue " , fnam );
  print_banner(BANNER);
  printf("\t Read WRONGHOST info from \n\t %s\n", inFile);

  // search top of file for POLY function of ZTRUE that
  // defines wrongHost PROB vs. ZTRUE

  char tmpWord[100];
  char key_prob[] = "PROB_WRONGHOST_POLY:" ;
  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,key_prob) == 0 ) {
      readchar(fp,tmpWord);
      if ( strstr(tmpWord,",") == NULL )  { 

	// if no comma, read legacy hard-coded 3rd order poly
	// with space-separated coefficients
	sscanf(tmpWord, "%le", &WRONGHOST.PROB_POLY[0]);
	readdouble(fp, ORDER_PROB_WRONGHOST_POLY-1, &WRONGHOST.PROB_POLY[1]);
	sprintf(tmpWord,"%f,%f,%f,%f",
		WRONGHOST.PROB_POLY[0], WRONGHOST.PROB_POLY[1],
		WRONGHOST.PROB_POLY[2], WRONGHOST.PROB_POLY[3]   );
      }

      // if last char of tmpWord is a comma, abort since space-separated
      // values are not allowed.
      int lentmp = strlen(tmpWord);
      if ( tmpWord[lentmp-1] == ',' ) {
	sprintf(c1err, "Invalid input %s %s", key_prob, tmpWord);
	sprintf(c2err, "pad spaces not allowed with commas.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      parse_GENPOLY(tmpWord, "z", &WRONGHOST.ZPOLY_PROB, fnam);
      goto NEXTREAD ;
    }

  }


  // if we get here, abort on error
  sprintf(c1err,"Could not find requred key = '%s'", key_prob);
  sprintf(c2err,"in WRONGHOST_FILE = '%s'", inFile);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);


 NEXTREAD:

  // since the file is open, read every line to estimate
  // how much memory to allocate
  rewind(fp);
  NROW_TOT = 0;
  while( fgets(c_get, 60, fp)  != NULL ) { NROW_TOT++ ; }
  fclose(fp);

  printf("\t Allocate ZTRUE & ZMATCH arrays of len %d\n", NROW_TOT); 
  fflush(stdout);

  double *tmpZTRUE, *tmpZMATCH ;
  MEMD = NROW_TOT * sizeof(double);

  INDEX_SORT   = (int*) malloc ( NROW_TOT * sizeof(int) );
  tmpZTRUE     = (double*) malloc ( MEMD );
  tmpZMATCH    = (double*) malloc ( MEMD );
  WRONGHOST.ZTRUE_LIST   = (double*) malloc ( MEMD );
  WRONGHOST.ZMATCH_LIST  = (double*) malloc ( MEMD );

 
  rd2columnFile(inFile, NROW_TOT,         // inputs
		&NROW_READ, tmpZTRUE, tmpZMATCH); // outputs

  printf("\t Read %d ZTRUE-ZMATCH pairs \n", NROW_READ);
  printf("\t Require %.3f < ZTRUE & ZMATCH < %.3f \n", ZMIN, ZMAX);

  // sort tmpZTRUE and store WRONGHOST.ZTRUE[ZMATCH]
  // as a ZTRUE-sorted list

  double ZTRUE, ZMATCH, DZ, DZMAX;
  int    iz_atDZMAX ;
  int ORDER_SORT = +1 ;        // increasing order                                  
  sortDouble( NROW_READ, tmpZTRUE, ORDER_SORT, INDEX_SORT ) ;

  // load global WRONGHOST arrays with ZTRUE-sorted list
  NLIST = 0 ;
  for(iz = 0; iz < NROW_READ; iz++ ) {
    isort  = INDEX_SORT[iz];
    ZTRUE  = tmpZTRUE[isort];
    ZMATCH = tmpZMATCH[isort];

    if ( ZTRUE   < ZMIN ) { continue ; }
    if ( ZTRUE   > ZMAX ) { continue ; }
    if ( ZMATCH  < ZMIN ) { continue ; }
    if ( ZMATCH  > ZMAX ) { continue ; }

    WRONGHOST.ZTRUE_LIST[NLIST]  = ZTRUE ;
    WRONGHOST.ZMATCH_LIST[NLIST] = ZMATCH ;
    NLIST++ ;

    /*
    if ( i < 20 ) {
      printf(" xxx i=%3d  NLIST=%3d  isort=%4d ZTRUE = %.4f \n",
	     i, NLIST, isort, tmpZTRUE[isort] );
	     } */

  }

  WRONGHOST.NLIST = NLIST ; 

  printf("\t Stored %d ZTRUE-ZMATCH pairs, sorted by ZTRUE. \n",  NLIST);

  // check max gap between adjacent ZTRUE
  DZMAX = -9.0 ;
  iz_atDZMAX = 0 ;
  for(iz=0; iz < NLIST-1 ; iz++ ) {
    DZ = WRONGHOST.ZTRUE_LIST[iz+1] - WRONGHOST.ZTRUE_LIST[iz] ;
    if ( DZ > DZMAX ) { DZMAX = DZ;   iz_atDZMAX = iz ; }
  }

  iz = iz_atDZMAX ;
  printf("\t Max ZTRUE gap is %f (ZTRUE=%f to %f)\n",
	 DZMAX, WRONGHOST.ZTRUE_LIST[iz], WRONGHOST.ZTRUE_LIST[iz+1] );
  printf("\t WRONGHOST range of ZTRUE(min,max): %f  to %f \n",
	 WRONGHOST.ZTRUE_LIST[0], WRONGHOST.ZTRUE_LIST[NLIST-1] );

  fflush(stdout);

  // free local memory
  free(tmpZTRUE);
  free(tmpZMATCH);
  free(INDEX_SORT);

  return;

} // end init_WRONGHOST


double gen_zHOST(int CID, double zSN, int *hostMatch) {

  // Dec 2015
  // use input zSN to generate host redshift (CID is for error msg only).
  // For correct host match, returns zHOST = zSN.
  // For mis-matched host, return zHOST based on user-input
  // WRONGHOST_FILE. 
  // Note that init_WRONGHOST must be called first.
  // 
  // Inputs:
  //    CID = integer candidate ID
  //    zSN = SN redshift
  //
  // Output:
  //    function returns zHOST to generate from HOSTLIB
  //    hostMatch = 1 for correct match (zHOST=zSN)
  //    hostMatch = 0 for mis-match
  //
  // Mar  6 2016: fix bug loop over iz < NLIST intead of NZMAX
  //
  // Apr 12 2019: 
  //   + refactor to use parse_GENPOLY
  //   + minor fix to iz loop to find iz_start & iz_end
  //

  int    NLIST = WRONGHOST.NLIST ;
  double zHOST = zSN ;

  // always burn the randoms to stay synced
  double FlatRan_prob   = getRan_Flat1(1);
  double FlatRan_zmatch = getRan_Flat1(1);

  double P_wrongHost, XNRAN ;
  int i ;

  char fnam[] = "gen_zHOST";

  // ------------ BEGIN --------------

  *hostMatch = 1;   // default is correct match
  if ( NLIST == 0 ) { return zHOST ; }

  // compute wrongHost prob at this zSN using polynom coeff.

  P_wrongHost = eval_GENPOLY(zSN,&WRONGHOST.ZPOLY_PROB,fnam);      

  // check for correct host match
  if ( P_wrongHost < FlatRan_prob ) { return zHOST ; }

  // here we have the wrong host, so select zMatch to wrong host.

#define NZJUMP_WRONGHOST 50  // junp this many z bins for rough ZTRUE check

  int PAST_END = 0 ;
  int iz, iz_start, iz_end, N_NEAR;
  int iz_list[5*NZJUMP_WRONGHOST] ;
  double ZTRUE, ZMATCH, DZ ;

  iz_start = iz_end = -9 ; iz=0;
  //  for ( iz=0; iz < NLIST ; iz+= NZJUMP_WRONGHOST ) {
  while ( iz < NLIST ) {
    ZTRUE = WRONGHOST.ZTRUE_LIST[iz] ;

    // xxxxxxx
    if ( CID == -52304 && ZTRUE > 0.44 ) {
      printf(" xxx iz=%d of %d: ZTRUE=%.5f  zSN=%.5f \n",
	     iz, NLIST, ZTRUE, zSN); fflush(stdout);
    }
    // xxxxx

    if ( ZTRUE > zSN ) {
      iz_start = iz - 2*NZJUMP_WRONGHOST ;
      iz_end   = iz + 1*NZJUMP_WRONGHOST ;
      goto SEARCH2_ZTRUE ;
    }
    
    iz += NZJUMP_WRONGHOST ;
    if ( iz >= NLIST && !PAST_END )  { iz = NLIST-1;  PAST_END++; }

  } // end iz

 SEARCH2_ZTRUE:
  
  if ( iz_start == -9 || iz_end == -9 ) {
    print_preAbort_banner(fnam);
    printf("\t user gen zMIN/zMAX = %f/%f \n", 
	   WRONGHOST.USER_zMIN,WRONGHOST.USER_zMAX) ;
    printf("\t WrongHost map ZMIN/ZMAX = %f/%f \n",
	   WRONGHOST.ZTRUE_LIST[0], WRONGHOST.ZTRUE_LIST[NLIST-1] );
    sprintf(c1err,"Invalid iz[start,end] = %d, %d", iz_start, iz_end);
    sprintf(c2err,"CID=%d  zSN=%f", CID, zSN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if ( iz_start < 0     )  { iz_start = 0 ; }
  if ( iz_end  >= NLIST )  { iz_end = NLIST-1; }

  // count how many have ZTRUE close enough
  N_NEAR = 0 ;
  for( iz=iz_start; iz <= iz_end; iz++ ) {
    ZTRUE = WRONGHOST.ZTRUE_LIST[iz] ;
    DZ    = fabs(ZTRUE-zSN);
    if ( DZ < DZTRUEMAX_WRONGHOST ) {
      iz_list[N_NEAR] = iz ;
      N_NEAR++ ;      
    }
  }

  if ( N_NEAR < 1 ) {
    sprintf(c1err,"Found %d WRONGHOST entries with dz<%.3f",
	    N_NEAR, DZTRUEMAX_WRONGHOST );
    sprintf(c2err,"CID=%d  zSN=%f  iz[start,end]=%d,%d", 
	    CID, zSN, iz_start, iz_end );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // pick random element is iz_list
  XNRAN = FlatRan_zmatch * (double)N_NEAR ;
  i = (int)XNRAN ;
  iz = iz_list[i];
  if ( iz < 0 || iz >= NLIST ) {
    sprintf(c1err,"Invalid iz=%d(i=%d) for N_NEAR=%d and FlatRan=%f",
	    iz, i, N_NEAR, FlatRan_zmatch);
    sprintf(c2err,"CID=%d and zSN=%f", CID, zSN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  ZTRUE  = WRONGHOST.ZTRUE_LIST[iz] ;
  ZMATCH = WRONGHOST.ZMATCH_LIST[iz] ;

  // set outputs
  zHOST  = zSN + ( ZMATCH - ZTRUE ); 

  // make sure to respect HOSTLIB z range
  if ( zHOST >= WRONGHOST.USER_zMAX ) { zHOST = WRONGHOST.USER_zMAX - 0.001 ; }
  if ( zHOST <= WRONGHOST.USER_zMIN ) { zHOST = WRONGHOST.USER_zMIN + 0.001 ; }

  *hostMatch = 0;

  /*
  printf(" xxx CID=%d : zSN=%.3f --> zHOST = %.3f \n",
  CID, zSN, zHOST);  fflush(stdout); */

  return(zHOST) ;

} // end gen_zHOST

