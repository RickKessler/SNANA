/***************************************
  Created July 2019 by R.Kessler and J.Peirel

  Model multiple images strong lensed images for SN.

  Aug 7 2019 RK - pass DUMPFLAG argument to get_stronglens

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
  char cline[200],tmpWord[64];
  int MSKOPT_PARSE = MSKOPT_PARSE_WORDS_STRING + MSKOPT_PARSE_WORDS_IGNORECOMMA;
  int iwd,NWD,i,j,k,NVARS,NIMG,Nsplit ;
  char *cptr[MXIMG_STRONGLENS];
  char comma[] = ",";
  char MISSING_VAR[40];
  int  MEMC  = 64*sizeof(char);
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



  sprintf(INPUTS_STRONGLENS.VARNAME_LENSID,"LENSID");
  sprintf(INPUTS_STRONGLENS.VARNAME_ZSRC,  "ZSRC"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_ZLENS, "ZLENS" );
  sprintf(INPUTS_STRONGLENS.VARNAME_NIMG,  "NIMG"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_XIMG,  "XIMG"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_YIMG,  "YIMG"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_MAG,    "MAG"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_DELAY,"DELAY"  );
  

  INPUTS_STRONGLENS.ICOL_LENSID = -9;
  INPUTS_STRONGLENS.ICOL_ZLENS  = -9;
  INPUTS_STRONGLENS.ICOL_ZSRC   = -9;
  INPUTS_STRONGLENS.ICOL_NIMG   = -9;
  INPUTS_STRONGLENS.ICOL_XIMG   = -9;
  INPUTS_STRONGLENS.ICOL_YIMG   = -9;
  INPUTS_STRONGLENS.ICOL_MAG    = -9;
  INPUTS_STRONGLENS.ICOL_DELAY  = -9;

  // - - - - - - - - - - - - - - - -
  // read LENS library below.
  
  // allocate pointers to strip comma-separated values; e.g., Ximg, Yimg ...
  for(k=0;k<MXIMG_STRONGLENS;++k){
    cptr[k] = (char*)malloc(MEMC); 
  }
    

  
  NVARS=-1;
  while(NVARS==-1 && fgets(cline, 200, fp)  != NULL ){
    NWD = store_PARSE_WORDS(MSKOPT_PARSE,cline);
    if(NWD>1){
      get_PARSE_WORD(0,0,tmpWord);
      if(strcmp(tmpWord,"VARNAMES:")==0){
	       NVARS=NWD-1;
      }
    }
  }

  // make sure that we found VARNAMES
  if(NVARS == -1){
    sprintf(c1err,"Could not find required VARNAMES key"); 
    sprintf(c2err,"Check %s", MODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);
  }

  // read and store each VARNAME
  char VARLIST[NVARS][40];   
  for(k=0;k<NVARS+1;++k){
    get_PARSE_WORD(0,k,VARLIST[k]);
    //  printf(" xxx %s: found VARLIST[%d] = '%s' \n", fnam, k, VARLIST[k]);

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_LENSID) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_LENSID = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_ZLENS) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_ZLENS = k; }    

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_ZSRC) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_ZSRC = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_NIMG) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_NIMG = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_XIMG) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_XIMG = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_YIMG) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_YIMG = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_MAG) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_MAG = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_DELAY) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_DELAY = k; }

  }

  

  if( INPUTS_STRONGLENS.ICOL_LENSID < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_LENSID); }
  else if( INPUTS_STRONGLENS.ICOL_ZLENS < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_ZLENS); }
  else if( INPUTS_STRONGLENS.ICOL_ZSRC < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_ZSRC); }
  else if( INPUTS_STRONGLENS.ICOL_NIMG < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_NIMG); }
  else if( INPUTS_STRONGLENS.ICOL_XIMG < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_XIMG); }
  else if( INPUTS_STRONGLENS.ICOL_YIMG < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_YIMG); }
  else if( INPUTS_STRONGLENS.ICOL_MAG < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_MAG); } 
  else if( INPUTS_STRONGLENS.ICOL_DELAY < 0 ){ sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_DELAY); }
  else{ sprintf(MISSING_VAR,"NONE"); }
  
  if( strcmp(MISSING_VAR,"NONE") != 0 ){
    sprintf(c1err,"Could not find required key: %s",MISSING_VAR); 
    sprintf(c2err,"Check %s", MODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);
  }
 
  i=0; // represents the library entry

  while( fgets(cline, 200, fp)  != NULL ){

      NWD = store_PARSE_WORDS(MSKOPT_PARSE,cline);
      iwd = 0;
      get_PARSE_WORD(0,iwd,tmpWord); // read first word

      if ( strcmp(tmpWord,"LENS:") != 0 ) { continue ; } // RK

      if ( NWD-1 != NVARS )  { // RK
	sprintf(c1err,"Found %i strings after LENS: key",NWD-1);
    	sprintf(c2err,"but expected %d strings", NVARS);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
  
      NIMG = -9;

      while(iwd<NVARS+1) {
      	get_PARSE_WORD(0,iwd,tmpWord);

      	if ( iwd == INPUTS_STRONGLENS.ICOL_LENSID ) 
      	  { sscanf(tmpWord, "%d", &INPUTS_STRONGLENS.IDLENS[i]); }

      	else if ( iwd == INPUTS_STRONGLENS.ICOL_ZLENS ) 
      	  { sscanf(tmpWord, "%f", &INPUTS_STRONGLENS.ZLENS[i]); }

        else if ( iwd == INPUTS_STRONGLENS.ICOL_ZSRC ) 
          { sscanf(tmpWord, "%f", &INPUTS_STRONGLENS.ZSRC[i]); }

      	else if ( iwd == INPUTS_STRONGLENS.ICOL_NIMG ) {
      	  sscanf(tmpWord, "%d", &INPUTS_STRONGLENS.NIMG[i]) ; 
      	  NIMG = INPUTS_STRONGLENS.NIMG[i] ;
      	}

      	else if ( iwd == INPUTS_STRONGLENS.ICOL_XIMG ) 	{ 
          splitString(tmpWord,comma,MXIMG_STRONGLENS,&Nsplit,cptr);
      	  if ( Nsplit != NIMG )  { 
             sprintf(c1err,"Found %d images but expected %d in line %d",Nsplit,NIMG,i);
             sprintf(c2err,"of %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }

      	  for(j=0; j<NIMG; ++j) {
                  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.XIMG[i][j] );
      	  }
      	 
      	}

        else if ( iwd == INPUTS_STRONGLENS.ICOL_YIMG )  { 
          splitString(tmpWord,comma,MXIMG_STRONGLENS,&Nsplit,cptr);
          if ( NIMG < 0 )  { 
             sprintf(c1err,"NIMG must come before variables with multiple images (e.g. MAG, DELAY, etc.)");
             sprintf(c2err,"in %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }
          else if ( Nsplit != NIMG )  { 
             sprintf(c1err,"Found %d images but expected %d in line %d",Nsplit,NIMG,i);
             sprintf(c2err,"of %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }

          for(j=0; j<NIMG; ++j) {
                  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.YIMG[i][j] );
          }
         
        }

        else if ( iwd == INPUTS_STRONGLENS.ICOL_MAG )  { 
          splitString(tmpWord,comma,MXIMG_STRONGLENS,&Nsplit,cptr);
          if ( NIMG < 0 )  { 
             sprintf(c1err,"NIMG must come before variables with multiple images (e.g. MAG, DELAY, etc.)");
             sprintf(c2err,"in %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }
          else if ( Nsplit != NIMG )  { 
             sprintf(c1err,"Found %d images but expected %d in line %d",Nsplit,NIMG,i);
             sprintf(c2err,"of %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }

          for(j=0; j<NIMG; ++j) {
                  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.MAG[i][j] );
          }
         
        }
      	
        else if ( iwd == INPUTS_STRONGLENS.ICOL_DELAY )  { 
	  splitString(tmpWord,comma,MXIMG_STRONGLENS,&Nsplit,cptr);
          if ( NIMG < 0 )  { 
             sprintf(c1err,"NIMG must come before variables with multiple images (e.g. MAG, DELAY, etc.)");
             sprintf(c2err,"in %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }
          else if ( Nsplit != NIMG )  { 
             sprintf(c1err,"Found %d images but expected %d in line %d",Nsplit,NIMG,i);
             sprintf(c2err,"of %s", MODEL_FILE);
             errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
          }

          for(j=0; j<NIMG; ++j) {
                  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.DELAY[i][j] );
          }
         
        }
      	iwd++ ;
            } // end loop over columns

          i++; // increment library entry

      } // end while over input lines 
    
  
  fclose(fp);

  INPUTS_STRONGLENS.NLENS = i-1;
  //printf("%i\n",INPUTS_STRONGLENS.NLENS);
  return ;

} // end init_stronglens

// ===================================
void malloc_stronglens(int NLENS) {

  // malloc memory to hold strong lens library.
  // Access 2D angle arrays as
  //  INPUTS_STRONGLENS.Ximg[libid][image]
  //  INPUTS_STRONGLENS.Yimg[libid][image]
  //

  int MEMFF = NLENS * sizeof(float*);
  int MEMF  = NLENS * sizeof(float);
  int MEMI  = NLENS * sizeof(int);
    
  int i;  
  //  char fnam[] = "malloc_stronglens";

  // ------------ BEGIN --------------

  INPUTS_STRONGLENS.IDLENS = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.ZLENS  = (float*) malloc(MEMF);
  INPUTS_STRONGLENS.NIMG = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.XIMG   = (float**)malloc(MEMFF); //RA Offset
  INPUTS_STRONGLENS.YIMG   = (float**)malloc(MEMFF); //RA*cos(DEC) offset
  INPUTS_STRONGLENS.DELAY = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.MAG     = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.ZSRC   = (float*) malloc(MEMF);


  int memf = MXIMG_STRONGLENS * sizeof(float);
  for(i=0; i < NLENS; i++ ) {
    INPUTS_STRONGLENS.XIMG[i]    = (float*)malloc(memf);
    INPUTS_STRONGLENS.YIMG[i]    = (float*)malloc(memf);
    INPUTS_STRONGLENS.DELAY[i]  = (float*)malloc(memf);
    INPUTS_STRONGLENS.MAG[i]      = (float*)malloc(memf);
  }

  return ;

} // end malloc_stronglens

// ==========================================
void get_stronglens(double zSN, double *hostpar, int DUMPFLAG,
		    int *IDLENS, double *ZLENS, int *blend_flag, int *NIMG, 
		    double *DELAY, double *MAG, double *XIMG, double *YIMG) {

  // Inputs:
  //   zSN       redshift of SN
  //   hostpar   placeholder for host properties (future upgrade?)
  //   DUMPFLAG  dump flag (Aug 7 2019, RK)
  //
  // Ouptuts:
  //  IDLENS      integer identifier for galaxy lens
  //  ZLENS       redshift of lens galaxy
  //  blend_flag  1 if blended into single LC; 0 if each image is resolved
  //  NIMG      Number of images
  //  DELAY      list of NIMG time delays (days)
  //  MAG          list of NIMG magnifications 
  //  XIMG        list of NIMG X separations (arcsec)
  //  YIMG        list of NIMG Y separations (arcsec)
  //
  int    IDLENS_local, NIMG_local=0, img,i,j;
  double FlatRan, GauRan, zLENS_local;
  char fnam[] = "get_stronglens" ;

  // ---------------- BEGIN ---------------

  // always burn random
  FlatRan = FlatRan1(2);   // flat between 0 and 1
  // GauRan  = GaussRan(2);   // Gaussian, sigma=1
  
  *NIMG = 0 ;
  if ( !INPUTS_STRONGLENS.USE_FLAG ) { return ; }

  INPUTS_STRONGLENS.NCALL++ ;

  
  int numLens = 0;
  for(i=0;i<INPUTS_STRONGLENS.NLENS;++i){
    if(INPUTS_STRONGLENS.ZSRC[i]>=zSN-0.05 && INPUTS_STRONGLENS.ZSRC[i]<=zSN+0.05){
      ++numLens;
    }

  }
  if(numLens==0){
    //errmsg(SEV_FATAL, 0, fnam, "No Lenses in your library matching your source redshift."," ");
    printf("No lenses found in library matching your source redshift %f\n",zSN);
    return;
  }

  int possible_lenses[numLens];
  j=0;
  for(i=0;i<INPUTS_STRONGLENS.NLENS;++i){
    if(INPUTS_STRONGLENS.ZSRC[i]>=zSN-0.05 && INPUTS_STRONGLENS.ZSRC[i]<=zSN+0.05){
      possible_lenses[j]=i;
      ++j;
    }
  }

  int random_lens_index = possible_lenses[ (int)( FlatRan*(numLens-1) ) ];

  IDLENS_local  = INPUTS_STRONGLENS.IDLENS[random_lens_index];
  zLENS_local   = (double)INPUTS_STRONGLENS.ZLENS[random_lens_index];  
  NIMG_local    = INPUTS_STRONGLENS.NIMG[random_lens_index];

  
  // strip off image-dependent quantities, and recast to double
  for(img=0 ; img < NIMG_local; img++ ) {
    XIMG[img]   = (double)INPUTS_STRONGLENS.XIMG[random_lens_index][img]; 
    YIMG[img]   = (double)INPUTS_STRONGLENS.YIMG[random_lens_index][img]; 
    MAG[img]     = (double)INPUTS_STRONGLENS.MAG[random_lens_index][img];     
    DELAY[img] = (double)INPUTS_STRONGLENS.DELAY[random_lens_index][img]; 
  }

  // load return argument scalars  
  *IDLENS = IDLENS_local;
  *ZLENS  = zLENS_local ;
  *NIMG = NIMG_local;


  if ( DUMPFLAG ) {
    printf(" xxx \n");
    printf(" xxx ------ %s DUMP -------- \n", fnam );
    printf(" xxx input zSN = %.3f \n", zSN);
    printf(" xxx output IDLENS=%d at zLENS=%.3f \n", 
	   IDLENS_local, zLENS_local );

    for(img=0 ; img < NIMG_local; img++ ) {
      printf(" xxx output image-%d: mu=%.3f deltaT=%.2f  X,Y-offset=%f,%f\n",
	     img, MAG[img],DELAY[img], XIMG[img], YIMG[img] );
    }
    printf(" xxx \n");    fflush(stdout);
  }

  return ;

} // end get_stronglens


