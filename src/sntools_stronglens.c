/***************************************
  Created July 2019 by R.Kessler and J.Peirel

  Model multiple images strong lensed images for SN.

 ***************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_stronglens.h"

// ==========================================
void init_stronglens(char *MODEL_FILE) {

  // Initialize strong lens model.
  FILE *fp;
  char fnam[] = "init_stronglens";

  // --------------- BEGIN ---------------

  INPUTS_STRONGLENS.USE_FLAG = 0 ;
  INPUTS_STRONGLENS.NCALL    = 0 ;
  if ( IGNOREFILE(MODEL_FILE) ) { return; }

  sprintf(BANNER,"%s", fnam);
  print_banner(BANNER);


  // open input file for reading
  fp = fopen(MODEL_FILE,"rt");

  if ( fp == NULL ) {
    sprintf(c1err, "Could not open strong-lens model file");
    sprintf(c2err, " '%s' ", MODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  INPUTS_STRONGLENS.USE_FLAG = 1;
  printf("\t Read strong lens model from\n\t %s\n", MODEL_FILE);
  fflush(stdout);


  // estimate NLENS with upper bound; e.g., read number of lines in file,
  // or read NLENS key from library.  Then allocate memory.
  int NLENS_APPROX = 50 ;
  malloc_stronglens(NLENS_APPROX);


  // read here ...



  fclose(fp);

  return ;

} // end init_stronglens

// ===================================
void malloc_stronglens(int NLENS) {

  // malloc memory to hold strong lens library.
  // Access 2D angle arrays as
  //  INPUTS_STRONGLENS.angSep[libid][image]
  //  INPUTS_STRONGLENS.phi[libid][image]
  //

  int MEMFF = NLENS * sizeof(float*);
  int MEMF  = NLENS * sizeof(float);
  int MEMI  = NLENS * sizeof(int);
  int i;
  //  char fnam[] = "malloc_stronglens";

  // ------------ BEGIN --------------

  INPUTS_STRONGLENS.IDLENS = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.zLENS  = (float*) malloc(MEMF);
  INPUTS_STRONGLENS.Nimage = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.angSep = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.phi    = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.tdelay = (float**)malloc(MEMFF);
  
  int memf = MXIMG_STRONGLENS * sizeof(float);
  for(i=0; i < NLENS; i++ ) {
    INPUTS_STRONGLENS.angSep[i] = (float*)malloc(memf);
    INPUTS_STRONGLENS.phi[i]    = (float*)malloc(memf);
    INPUTS_STRONGLENS.tdelay[i] = (float*)malloc(memf);
  }

  return ;

} // end malloc_stronglens

// ==========================================
int get_stronglens(double zSN, double *hostpar, 
		   double *zLENS, int *blend_flag, int *Nimage, 
		   double *tdelay, double *mu, double *angSep, double *phi) {

  // Inputs:
  //   zSN       redshift of SN
  //   hostpar   placeholder for host properties (future upgrade?)
  //
  // Ouptuts:
  //  zLENS       redshift of lens galaxy
  //  blend_flag  1 if blended into single LC; 0 if each image is resolved
  //  Nimage      Number of images
  //  tdelay      list of Nimage time delays (days)
  //  mu          list of Nimage magnifications (not distance modulus)
  //  angSep      list of Nimage angular separations from unlensed SN (arcsec)
  //  phi         list of Nimage phi angle, degrees (phi=0 along RA)
  //
  // Function returns integer ID of lens galaxy. This ID is written
  // to data files, and is used to link multiple light curves to 
  // same lensing event.

  int    IDLENS, Nimage_local=0, img;
  double FlatRan, GauRan, zLENS_local;
  char fnam[] = "get_stronglens" ;

  // ---------------- BEGIN ---------------

  *Nimage = 0 ;
  if ( !INPUTS_STRONGLENS.USE_FLAG ) { return(0); }

  INPUTS_STRONGLENS.NCALL++ ;

  // illustrate how to generate randoms
  FlatRan = FlatRan1(2);   // flat between 0 and 1
  GauRan  = GaussRan(2);   // Gaussian, sigma=1

  zLENS_local = zSN * ( 0.1 + 0.8*FlatRan ) ;
 
  Nimage_local = (int)(FlatRan1(2) * 5.0);
  for(img=0; img < Nimage_local; img++ ) {
    tdelay[img] = 100.0 * (FlatRan1(2)-0.5) ;
    angSep[img] = fabs(GaussRan(2));
    phi[img]    = 360.0 * FlatRan1(2);
    mu[img]     = 1.0 + 0.2*GaussRan(2) ;
    if ( mu[img] < .1 ) { mu[img] = 0.1; }    
  }

  // load return arguments
  IDLENS  = INPUTS_STRONGLENS.NCALL ;
  *zLENS  = zLENS_local ;
  *Nimage = Nimage_local ;


  int  LDMP = (IDLENS < 0 ) ;
  if ( LDMP ) {
    printf(" xxx ------ %s DUMP -------- \n", fnam );
    printf(" xxx input zSN = %.3f \n", zSN);
    printf(" xxx output IDLENS=%d at zLENS=%.3f \n", IDLENS, *zLENS);
    for(img=0 ; img < Nimage_local; img++ ) {
      printf(" xxx output image-%d: mu=%.3f  angSep=%.2f  phi=%.1f\n",
	     img, mu[img], angSep[img], phi[img] );
    }
    printf(" xxx \n");    fflush(stdout);
  }

  return(IDLENS);

} // end get_stronglens


