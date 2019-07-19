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
  char c_get[60];
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
  int NROWS = nrow_read(MODEL_FILE,fnam);
  int NLENS_APPROX = 2*NROWS ;
  malloc_stronglens(NLENS_APPROX);

  int STOP=0.;
  int n_spaces = 0, i;
  char ** indices = NULL;
  
  fscanf(fp,"%[^\n]",c_get);
  char *p = strtok(c_get," ");
  while(p){
    indices = realloc (indices, sizeof (char*) * ++n_spaces);
    if (indices == NULL)
      exit (-1);
    indices[n_spaces-1] = p;
  
    p = strtok (NULL, " ");
  }
  indices = realloc (indices, sizeof (char*) * (n_spaces+1));
  indices[n_spaces] = 0;
  printf("%s\n",indices[0]);
  printf("%s\n",indices[1]);
  //STOP = 1;
  //char ** templine = NULL;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  float test[n_spaces] = {0.0};

  //char empty = ;
  //int n_spaces_temp = 0;
  //while((read = getline(&line,&len,fp)) != -1){
  //  if(strlen(line)>=n_spaces){
  //    for(i = 0; i<(n_spaces);++i){
  //	//p=strtok(line," ");
  //	printf("%s",strtok(line," "));
  //	//p = strtok (NULL, " ");
  //  }
  //}
    //n_spaces_temp=0;
    
    //while(p){
    //  templine = realloc (templine, sizeof (char*) * ++n_spaces_temp);
    //templine[n_spaces_temp-1] = p;
    //p = strtok(NULL," ");
    //}
    //templine = realloc (indices, sizeof (char*) * (n_spaces+1));
    //templine[n_spaces] = 0;
    //for(i = 0; i<(n_spaces);++i){
    //  printf("%s",templine[i]);
    //}
  //}
  
  while(!STOP){
    STOP=1;
      //if ( strcmp(c_get,"VARNAMES:") == 0 ) { ; }
    fscanf(fp,"%[^\n]",c_get);
    printf("%s\n",c_get);
    if (c_get == NULL){
      printf("%s\n","Reached end");
      STOP = 1.;
    } else{
      char *p = strtok(c_get," ");
      n_spaces = 0;
      for(i = 0; i< (n_spaces);++i){
	printf("%s\n",p);
	  //templine[i] = p;
	  //++n_spaces;

	p = strtok (NULL, " ");
      }
      //for (i = 0; i < (n_spaces+1); ++i)
      //  if (res[i] == 'IDLENS')
	    
      //    }
      //printf("%s\n",templine[0]);
    //printf("%s\n",templine[4]);
    //printf("%s\n",p);
    }
  }
  

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
  char fnam[] = "malloc_stronglens";

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
  if ( !INPUTS_STRONGLENS.USE_FLAG ) { return; }

  INPUTS_STRONGLENS.NCALL++ ;
 

  // illustrate how to generate randoms
  FlatRan = FlatRan1(2);   // flat between 0 and 1
  GauRan  = GaussRan(2);   // Gaussian, sigma=1

  zLENS_local = zSN * ( 0.1 + 0.8*FlatRan ) ;
 
  Nimage_local = (int)(FlatRan1(2) * 5.0);
  for(img=0; img < Nimage_local; img++ ) {
    angSep[img] = fabs(GaussRan(2));
    phi[img]    = 360.0 * FlatRan1(2);
    mu[img]     = 1.0 + 0.2*GaussRan(2) ;
    if ( mu[img] < .1 ) { mu[img] = 0.1; }    
  }

  // load return arguments
  IDLENS  = INPUTS_STRONGLENS.NCALL ;
  *zLENS  = zLENS_local ;
  *Nimage = Nimage_local ;


  int  LDMP = (IDLENS > 0 ) ;
  if ( LDMP ) {
    printf(" xxx ------ %s DUMP -------- \n", fnam );
    printf(" xxx input zSN = %.3f \n", zSN);
    printf(" xxx output IDLENS=%d at zLENS=%.3f \n", IDLENS, zLENS);
    for(img=0 ; img < Nimage_local; img++ ) {
      printf(" xxx output image-%d: mu=%.3f  angSep=%.2f  phi=%.1f\n",
	     img, mu[img], angSep[img], phi[img] );
    }
    printf(" xxx \n");    fflush(stdout);
  }

  return(IDLENS);

} // end get_stronglens


