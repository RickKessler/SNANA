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
  char cline[200];
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

  //int STOP=0.;
  //int n_spaces = 0, i;
  //char ** indices = NULL;
  
  //fscanf(fp,"%[^\n]",c_get);
  //char *p = strtok(c_get," ");
  //while(p){
  //  indices = realloc (indices, sizeof (char*) * ++n_spaces);
  //  if (indices == NULL)
  //    exit (-1);
  //  indices[n_spaces-1] = p;
  
  //  p = strtok (NULL, " ");
  //}
  //indices = realloc (indices, sizeof (char*) * (n_spaces+1));
  //indices[n_spaces] = 0;
  //printf("%s\n",indices[0]);
  //printf("%s\n",indices[1]);
  //STOP = 1;
  //char ** templine = NULL;
  //char * line = NULL;
  //size_t len = 0;
  //ssize_t read;
  //float test[n_spaces] = {0.0};

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
  
  //while(!STOP){
  //  STOP=1;
      //if ( strcmp(c_get,"VARNAMES:") == 0 ) { ; }
  //  fscanf(fp,"%[^\n]",c_get);
  //  printf("%s\n",c_get);
  //  if (c_get == NULL){
  //    printf("%s\n","Reached end");
  //    STOP = 1.;
  //  } else{
  //    char *p = strtok(c_get," ");
  //    n_spaces = 0;
  //    for(i = 0; i< (n_spaces);++i){
  //	printf("%s\n",p);
	  //templine[i] = p;
	  //++n_spaces;

  //	p = strtok (NULL, " ");
  //    }
      //for (i = 0; i < (n_spaces+1); ++i)
      //  if (res[i] == 'IDLENS')
	    
      //    }
      //printf("%s\n",templine[0]);
    //printf("%s\n",templine[4]);
    //printf("%s\n",p);
    //}
  //}
  int iwd,NWD,i,j,k,NVARS,STOP,NIMG,Nsplit,VARLINE;
  char tmpWord[64];
  char *cptr;
  char semicolon[] = ";";
  int    MEMC  = 64*sizeof(char);

  i=0;
  NVARS=-1;
  while(NVARS==-1 && fgets(cline, 200, fp)  != NULL ){
    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,cline);
    if(NWD>1){
      get_PARSE_WORD(0,0,tmpWord);
      if(strcmp(tmpWord,"VARNAMES:")==0){
	NVARS=NWD;
    }
  }
  }
  if(NVARS==-1){
    errmsg(SEV_FATAL, 0, fnam, "Error reading strong lens file, make sure you have a 'VARNAMES:' row.", c2err );
  }
  const char VARLIST[NVARS][20];
  
  for(k=0;k<NVARS;++k){
    get_PARSE_WORD(0,k,VARLIST[k]);
      }
  while( fgets(cline, 200, fp)  != NULL ){
    if(i>0){
      NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,cline);
      if(NWD==NVARS){
	STOP = 0;
      }else{
	STOP = 1;
      }
    }else{STOP = 0;}
    
    
    iwd=0;
    NIMG=-1;
    char *cptr[12];
    for(k=0;k<12;++k){
      cptr[k] = (char*)malloc(MEMC);
    }
    while(iwd<NVARS && !STOP){
      
      if(iwd==0){
	get_PARSE_WORD(0,iwd,tmpWord);
	if(strcmp(tmpWord,"LENS:")!=0){
	  STOP = 1;
	}
      }else{
	get_PARSE_WORD(0,iwd,tmpWord);
	if(strcmp(VARLIST[iwd],"LENSID")==0){
	  INPUTS_STRONGLENS.IDLENS[i] = atoi(tmpWord);
	}else if(strcmp(VARLIST[iwd],"NIMG")==0){
	  INPUTS_STRONGLENS.Nimage[i] = atoi(tmpWord);
	  NIMG = atoi(tmpWord);
	  
	}else if(strcmp(VARLIST[iwd],"ZLENS")==0){
	  INPUTS_STRONGLENS.zLENS[i] = atof(tmpWord);
	}else if(strcmp(VARLIST[iwd],"ZSRC")==0){
          INPUTS_STRONGLENS.zSRC[i] = atof(tmpWord);
        }else if(strcmp(VARLIST[iwd],"XIMG")==0){
	  if(NIMG==-1){
            errmsg(SEV_FATAL, 0, fnam, "Error reading strong lens file, make sure your NIMG and LENSID come first.", c2err );
          }
	  
	  
	  
	  
          splitString(tmpWord,semicolon,NIMG,&Nsplit,cptr);
	  for(j=0;j<NIMG;++j){
            INPUTS_STRONGLENS.Ximg[i][j] = atof(cptr[j]);
	  }
	}else if(strcmp(VARLIST[iwd],"YIMG")==0){
	  if(NIMG==-1){
            errmsg(SEV_FATAL, 0, fnam, "Error reading strong lens file, make sure your NIMG and LENSID come first.", c2err );
          }
	  
          splitString(tmpWord,";",NIMG,&Nsplit,cptr);
	  for(j=0;j<Nsplit;++j){
            INPUTS_STRONGLENS.Yimg[i][j] = atof(cptr[j]);
	  }
	}else if(strcmp(VARLIST[iwd],"MAG")==0){
	  if(NIMG==-1){
            errmsg(SEV_FATAL, 0, fnam, "Error reading strong lens file, make sure your NIMG and LENSID come first.", c2err );
          }
	  
          splitString(tmpWord,";",NIMG,&Nsplit,cptr);
	  for(j=0;j<Nsplit;++j){
            INPUTS_STRONGLENS.mu[i][j] = atof(cptr[j]);
	  } 
	}else if(strcmp(VARLIST[iwd],"DELAY")==0){
	    if(NIMG==-1){
	      errmsg(SEV_FATAL, 0, fnam, "Error reading strong lens file, make sure your NIMG and LENSID come first.", c2err );
	    }
	    splitString(tmpWord,";",NIMG,&Nsplit,cptr);
	    for(j=0;j<Nsplit;++j){
	      INPUTS_STRONGLENS.tdelay[i][j] = atof(cptr[j]);
	    }
	}
      }
	++iwd;
	
    }
    
    ++i;
  }
    
  
  fclose(fp);

  INPUTS_STRONGLENS.NLENS = i-2;
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
  INPUTS_STRONGLENS.NLENS = (int)malloc(sizeof(int));
  INPUTS_STRONGLENS.IDLENS = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.zLENS  = (float*) malloc(MEMF);
  INPUTS_STRONGLENS.Nimage = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.Ximg = (float**)malloc(MEMFF); //RA Offset
  INPUTS_STRONGLENS.Yimg    = (float**)malloc(MEMFF); //RA*cos(DEC) offset
  INPUTS_STRONGLENS.tdelay = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.mu = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.zSRC = (float*) malloc(MEMF);


  int memf = MXIMG_STRONGLENS * sizeof(float);
  for(i=0; i < NLENS; i++ ) {
    INPUTS_STRONGLENS.Ximg[i] = (float*)malloc(memf);
    INPUTS_STRONGLENS.Yimg[i]    = (float*)malloc(memf);
    INPUTS_STRONGLENS.tdelay[i] = (float*)malloc(memf);
    INPUTS_STRONGLENS.mu[i] = (float*)malloc(memf);
  }

  return ;

} // end malloc_stronglens

// ==========================================
int get_stronglens(double zSN, double *hostpar, 
		   double *zLENS, int *blend_flag, int *Nimage, 
		   double *tdelay, double *mu, double *Ximg, double *Yimg) {

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
  //  Ximg      list of Nimage RA separation from center of lens (arcsec)
  //  Yimg         list of Nimage RA*cost(DEC) separation from center of lens (arcsec)
  //
  // Function returns integer ID of lens galaxy. This ID is written
  // to data files, and is used to link multiple light curves to 
  // same lensing event.

  int    IDLENS, Nimage_local=0, img,i,j;
  double FlatRan, GauRan, zLENS_local;
  char fnam[] = "get_stronglens" ;

  // ---------------- BEGIN ---------------

  // always burn random
  FlatRan = FlatRan1(2);   // flat between 0 and 1
  // GauRan  = GaussRan(2);   // Gaussian, sigma=1
  
  *Nimage = 0 ;
  if ( !INPUTS_STRONGLENS.USE_FLAG ) { return(0); }

  INPUTS_STRONGLENS.NCALL++ ;

  
  int numLens = 0;
  for(i=0;i<INPUTS_STRONGLENS.NLENS;++i){
    if(INPUTS_STRONGLENS.zSRC[i]>=zSN-0.05 && INPUTS_STRONGLENS.zSRC[i]<=zSN+0.05){
      ++numLens;
    }

  }
  if(numLens==0){
    errmsg(SEV_FATAL, 0, fnam, "No Lenses in your library matching your source redshift."," ");
  }
  const int *possible_lenses[numLens];
  j=0;
  for(i=0;i<INPUTS_STRONGLENS.NLENS;++i){
    if(INPUTS_STRONGLENS.zSRC[i]>=zSN-0.05 && INPUTS_STRONGLENS.zSRC[i]<=zSN+0.05){
      possible_lenses[j]=i;
      ++j;
    }
  }
  int random_lens_index=(int)( FlatRan*(numLens-1) );
  //zLENS_local = zSN * ( 0.1 + 0.8*FlatRan ) ;
  zLENS = INPUTS_STRONGLENS.zLENS[random_lens_index];
  
  Nimage_local = INPUTS_STRONGLENS.Nimage[random_lens_index];
  Ximg = INPUTS_STRONGLENS.Ximg[random_lens_index];
  Yimg = INPUTS_STRONGLENS.Yimg[random_lens_index];
  mu = INPUTS_STRONGLENS.mu[random_lens_index];
  tdelay = INPUTS_STRONGLENS.tdelay[random_lens_index];
  // load return arguments
  
  IDLENS = INPUTS_STRONGLENS.IDLENS[random_lens_index];
  
  int  LDMP = 1 ;
  if ( LDMP ) {
    printf(" xxx ------ %s DUMP -------- \n", fnam );
    printf(" xxx input zSN = %.3f \n", zSN);
    printf(" xxx output IDLENS=%d at zLENS=%.3f \n", IDLENS, INPUTS_STRONGLENS.zLENS[random_lens_index]);

    for(img=0 ; img < Nimage_local; img++ ) {
      printf(" xxx output image-%d: mu=%.3f  X-offset=%.5f  Y-offset=%.5f\n",
	     img, INPUTS_STRONGLENS.mu[random_lens_index][img], INPUTS_STRONGLENS.Ximg[random_lens_index][img], INPUTS_STRONGLENS.Yimg[random_lens_index][img] );
    }
    printf(" xxx \n");    fflush(stdout);
  }

  return(IDLENS);

} // end get_stronglens


