/***************************************
  Created July 2019 by R.Kessler and J.Peirel

  Model multiple images strong lensed images for SN.

  Aug 7 2019 RK - pass DUMPFLAG argument to get_stronglens

  Jul 16 2020 JP - Added zSRC to DUMPFLAG

  Jun 30 2022 RK - 
    start adding code to read LOGMASS_LENS and LOGMASS_ERR_LENS
    to enable selecting appropriate LENS galaxy as observed host.

 ***************************************/


#include "sntools.h"
#include "sntools_stronglens.h"
#include "sntools_host.h" 

double prob_stronglens(double z) {

  // Created Jul 5 2022 by R.Kessler and J.Pierel
  // Use Eq 31 (Fig 5) approx from https://arxiv.org/pdf/1907.06830.pdf 
  //
  double zcube        = z * z * z;
  double z11          = pow(z,1.1) ;
  double top          = 5.0E-4 * zcube;
  double bottom_fac   = (1.0 + 0.41*z11);
  double bottom       = pow(bottom_fac,2.7) ;
  double prob         = top/bottom;
  return prob ;

} // end prob_stronglens

// ==========================================
void init_stronglens(char *MODEL_FILE) {

  // Initialize strong lens model by reading MODEL_FILE and
  // storing contents.

  FILE *fp;
  char fnam[] = "init_stronglens";
  int MXCHAR_LINE = 300;
  char cline[MXCHAR_LINE],tmpWord[100];
  int MSKOPT_PARSE = MSKOPT_PARSE_WORDS_STRING + MSKOPT_PARSE_WORDS_IGNORECOMMA;
  int iwd,NWD,i,j,k,NVARS,NIMG,Nsplit, gzipFlag ;
  char *cptr[MXIMG_STRONGLENS];
  char comma[] = ",";
  char MISSING_VAR[40], FULLNAME_MODEL_FILE[MXPATHLEN] ;
  int  MEMC  = 64*sizeof(char);
  // --------------- BEGIN ---------------

  INPUTS_STRONGLENS.USE_FLAG = 0 ;
  INPUTS_STRONGLENS.NCALL    = 0 ;
  if ( IGNOREFILE(MODEL_FILE) ) { return; }

  sprintf(BANNER,"%s", fnam);
  print_banner(BANNER);


  // open input file for reading
  fp = snana_openTextFile(1,PATH_USER_INPUT, MODEL_FILE,
			  FULLNAME_MODEL_FILE, &gzipFlag );

  if ( fp == NULL ) {
    abort_openTextFile("STRONGLENS_FILE",
		       PATH_USER_INPUT, MODEL_FILE, fnam );
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
  sprintf(INPUTS_STRONGLENS.VARNAME_LOGMASS_LENS,     "LOGMASS_LENS" );
  sprintf(INPUTS_STRONGLENS.VARNAME_LOGMASS_ERR_LENS, "LOGMASS_ERR_LENS" );
  sprintf(INPUTS_STRONGLENS.VARNAME_NIMG,     "NIMG"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_XIMG_SRC, "XIMG_SRC"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_YIMG_SRC, "YIMG_SRC"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_XGAL_SRC, "XGAL_SRC"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_YGAL_SRC, "YGAL_SRC"  );
  sprintf(INPUTS_STRONGLENS.VARNAME_MAGNIF,   "MAGNIF"    );
  sprintf(INPUTS_STRONGLENS.VARNAME_DELAY,    "DELAY"     );
  

  INPUTS_STRONGLENS.ICOL_LENSID = -9;
  INPUTS_STRONGLENS.ICOL_ZLENS  = -9;
  INPUTS_STRONGLENS.ICOL_LOGMASS_LENS      = -9;
  INPUTS_STRONGLENS.ICOL_LOGMASS_ERR_LENS  = -9;
  INPUTS_STRONGLENS.ICOL_ZSRC       = -9;
  INPUTS_STRONGLENS.ICOL_NIMG       = -9;
  INPUTS_STRONGLENS.ICOL_XIMG_SRC   = -9;
  INPUTS_STRONGLENS.ICOL_YIMG_SRC   = -9;
  INPUTS_STRONGLENS.ICOL_XGAL_SRC   = -9;
  INPUTS_STRONGLENS.ICOL_YGAL_SRC   = -9;
  INPUTS_STRONGLENS.ICOL_MAGNIF     = -9;
  INPUTS_STRONGLENS.ICOL_DELAY      = -9;

  // - - - - - - - - - - - - - - - -
  // read LENS library below.
  
  // allocate pointers to strip comma-separated values; e.g., Ximg, Yimg ...
  for(k=0;k<MXIMG_STRONGLENS;++k)
    {  cptr[k] = (char*)malloc(MEMC);   }
      
  // read until VARNAMES key is found
  NVARS = -1;
  while( NVARS==-1 && fgets(cline, MXCHAR_LINE, fp)  != NULL ){
    NWD = store_PARSE_WORDS(MSKOPT_PARSE,cline, fnam );
    if ( NWD > 1 ) {
      get_PARSE_WORD(0,0,tmpWord);
      if ( strcmp(tmpWord,"VARNAMES:") ==0  ) { NVARS=NWD-1; }
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
  for(k=0; k < NVARS+1; k++ ) {
    get_PARSE_WORD(0,k,VARLIST[k]);
    //  printf(" xxx %s: found VARLIST[%d] = '%s' \n", fnam, k, VARLIST[k]);

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_LENSID) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_LENSID = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_ZLENS) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_ZLENS = k; }    

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_LOGMASS_LENS) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_LOGMASS_LENS = k; }    
    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_LOGMASS_ERR_LENS) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_LOGMASS_ERR_LENS = k; }    

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_ZSRC) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_ZSRC = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_NIMG) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_NIMG = k; }

    // - - - -
    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_XIMG_SRC) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_XIMG_SRC = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_YIMG_SRC) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_YIMG_SRC = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_XGAL_SRC) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_XGAL_SRC = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_YGAL_SRC) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_YGAL_SRC = k; }

    // - - - 
    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_MAGNIF) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_MAGNIF = k; }

    if ( strcmp(VARLIST[k],INPUTS_STRONGLENS.VARNAME_DELAY) == 0 ) 
      { INPUTS_STRONGLENS.ICOL_DELAY = k; }

  }

  

  if( INPUTS_STRONGLENS.ICOL_LENSID < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_LENSID); }
  else if( INPUTS_STRONGLENS.ICOL_ZLENS < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_ZLENS); }

  // June 30 LOGMASS_LENS is optional 

  else if( INPUTS_STRONGLENS.ICOL_ZSRC < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_ZSRC); }
  else if( INPUTS_STRONGLENS.ICOL_NIMG < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_NIMG); }
  else if( INPUTS_STRONGLENS.ICOL_XIMG_SRC < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_XIMG_SRC); }
  else if( INPUTS_STRONGLENS.ICOL_YIMG_SRC < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_YIMG_SRC); }

  // XGAL_SRC and YGAL_SEC are optional

  else if( INPUTS_STRONGLENS.ICOL_MAGNIF < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_MAGNIF); } 
  else if( INPUTS_STRONGLENS.ICOL_DELAY < 0 )
    { sprintf(MISSING_VAR,INPUTS_STRONGLENS.VARNAME_DELAY); }
  else{ sprintf(MISSING_VAR,"NONE"); }
  
  if( strcmp(MISSING_VAR,"NONE") != 0 ){
    sprintf(c1err,"Could not find required key: %s",MISSING_VAR); 
    sprintf(c2err,"Check %s", MODEL_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err,c2err);
  }
 
  if ( INPUTS_STRONGLENS.ICOL_LOGMASS_LENS < 0 ) {
    printf("\t WARNING: LOGMASS_LENS not in SL library \n");
    fflush(stdout);
  }


  i=0; // represents the library entry

  while( fgets(cline, MXCHAR_LINE, fp)  != NULL ){

    NWD = store_PARSE_WORDS(MSKOPT_PARSE,cline, fnam );
    iwd = 0;
    get_PARSE_WORD(0,iwd,tmpWord); // read first word

    if ( strcmp(tmpWord,"LENS:") != 0 ) { continue ; } // RK

    if ( NWD-1 != NVARS )  { // RK
      sprintf(c1err,"Found %i strings after LENS: key",NWD-1);
      sprintf(c2err,"but expected %d strings", NVARS);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
  
    NIMG = -9 ;
    INPUTS_STRONGLENS.IDLENS[i]           = -9 ;
    INPUTS_STRONGLENS.ZLENS[i]            = -9.0 ;
    INPUTS_STRONGLENS.LOGMASS_LENS[i]     = -9.0 ;
    INPUTS_STRONGLENS.LOGMASS_ERR_LENS[i] = -9.0 ;
    INPUTS_STRONGLENS.ZSRC[i]             = -9.0 ;
    INPUTS_STRONGLENS.NIMG[i]             =  0;

    while ( iwd < NVARS+1 ) {
      get_PARSE_WORD(0,iwd,tmpWord);

      if ( iwd == INPUTS_STRONGLENS.ICOL_LENSID ) 
	{ sscanf(tmpWord, "%lld", &INPUTS_STRONGLENS.IDLENS[i]); }
		
      else if ( iwd == INPUTS_STRONGLENS.ICOL_ZLENS ) 
	{ sscanf(tmpWord, "%f", &INPUTS_STRONGLENS.ZLENS[i]); }

      else if ( iwd == INPUTS_STRONGLENS.ICOL_LOGMASS_LENS ) 
	{ sscanf(tmpWord, "%f", &INPUTS_STRONGLENS.LOGMASS_LENS[i]); }
      else if ( iwd == INPUTS_STRONGLENS.ICOL_LOGMASS_ERR_LENS ) 
	{ sscanf(tmpWord, "%f", &INPUTS_STRONGLENS.LOGMASS_ERR_LENS[i]); }

      else if ( iwd == INPUTS_STRONGLENS.ICOL_ZSRC ) 
	{ sscanf(tmpWord, "%f", &INPUTS_STRONGLENS.ZSRC[i]); }

      else if ( iwd == INPUTS_STRONGLENS.ICOL_NIMG ) {
	sscanf(tmpWord, "%d", &INPUTS_STRONGLENS.NIMG[i]) ; 
	NIMG = INPUTS_STRONGLENS.NIMG[i] ;
      }

      else if ( iwd == INPUTS_STRONGLENS.ICOL_XIMG_SRC ) 	{ 
	splitString(tmpWord, comma, fnam, MXIMG_STRONGLENS,&Nsplit,cptr);
	if ( Nsplit != NIMG )  { 
	  sprintf(c1err,"Found %d images but expected %d in line %d",
		  Nsplit,NIMG,i);
	  sprintf(c2err,"of %s", MODEL_FILE);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}
	
	for(j=0; j<NIMG; ++j) {
	  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.XIMG_SRC[i][j] );
	}
      	 
      }

      else if ( iwd == INPUTS_STRONGLENS.ICOL_YIMG_SRC )  { 
	splitString(tmpWord, comma, fnam, MXIMG_STRONGLENS,&Nsplit,cptr);
	if ( NIMG < 0 )  { 
	  sprintf(c1err,"NIMG must be defined before variables with "
		  "multiple images (e.g. MAGNIF, DELAY, etc.)");
	  sprintf(c2err,"in %s", MODEL_FILE);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}
	else if ( Nsplit != NIMG )  { 
	  sprintf(c1err,"Found %d images but expected %d in line %d",
		  Nsplit,NIMG,i);
	  sprintf(c2err,"of %s", MODEL_FILE);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}

	for(j=0; j<NIMG; ++j) {
	  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.YIMG_SRC[i][j] );
	}
         
      }

      else if ( iwd == INPUTS_STRONGLENS.ICOL_MAGNIF )  { 
	splitString(tmpWord, comma, fnam, MXIMG_STRONGLENS,&Nsplit,cptr);
	if ( NIMG < 0 )  { 
	  sprintf(c1err,"NIMG must be defined before variables with "
		  "multiple images (e.g. MAGNIF, DELAY, etc.)");
	  sprintf(c2err,"in %s", MODEL_FILE);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}
	else if ( Nsplit != NIMG )  { 
	  sprintf(c1err,"Found %d images but expected %d in line %d",
		  Nsplit,NIMG,i);
	  sprintf(c2err,"of %s", MODEL_FILE);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}

	for(j=0; j<NIMG; ++j) {
	  sscanf(cptr[j],"%f", &INPUTS_STRONGLENS.MAGNIF[i][j] );
	}
         
      }
      	
      else if ( iwd == INPUTS_STRONGLENS.ICOL_DELAY )  { 
	splitString(tmpWord, comma, fnam, MXIMG_STRONGLENS,&Nsplit,cptr);
	if ( NIMG < 0 )  { 
	  sprintf(c1err,"NIMG must be defined before variables with "
		  "multiple images (e.g. MAGNIF, DELAY, etc.)");
	  sprintf(c2err,"in %s", MODEL_FILE);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
	}
	else if ( Nsplit != NIMG )  { 
	  sprintf(c1err,"Found %d images but expected %d in line %d",
		  Nsplit,NIMG,i);
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
  char fnam[] = "malloc_stronglens";

  // ------------ BEGIN --------------

  INPUTS_STRONGLENS.IDLENS            = (long long int  *) malloc(MEMI);
  INPUTS_STRONGLENS.ZLENS             = (float*) malloc(MEMF);
  INPUTS_STRONGLENS.LOGMASS_LENS      = (float*) malloc(MEMF);
  INPUTS_STRONGLENS.LOGMASS_ERR_LENS  = (float*) malloc(MEMF);
  INPUTS_STRONGLENS.XGAL_SRC          = (float*) malloc(MEMF); 
  INPUTS_STRONGLENS.YGAL_SRC          = (float*) malloc(MEMF); 
  INPUTS_STRONGLENS.NIMG              = (int  *) malloc(MEMI);
  INPUTS_STRONGLENS.XIMG_SRC   = (float**)malloc(MEMFF); 
  INPUTS_STRONGLENS.YIMG_SRC   = (float**)malloc(MEMFF); 
  INPUTS_STRONGLENS.DELAY      = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.MAGNIF     = (float**)malloc(MEMFF);
  INPUTS_STRONGLENS.ZSRC       = (float*) malloc(MEMF);


  int memf = MXIMG_STRONGLENS * sizeof(float);
  for(i=0; i < NLENS; i++ ) {
    INPUTS_STRONGLENS.XIMG_SRC[i]    = (float*)malloc(memf);
    INPUTS_STRONGLENS.YIMG_SRC[i]    = (float*)malloc(memf);
    INPUTS_STRONGLENS.DELAY[i]       = (float*)malloc(memf);
    INPUTS_STRONGLENS.MAGNIF[i]      = (float*)malloc(memf);
  }

  return ;

} // end malloc_stronglens

// ==========================================
void get_stronglens(double zSN, double *hostpar, int DUMPFLAG,
		    EVENT_STRONGLENS_DEF *EVENT_SL ) {

  //
  // For input redshift zSN, return information about strong lens.
  //
  // Inputs:
  //   zSN       redshift of SN
  //   hostpar   placeholder for host properties (future upgrade?)
  //   DUMPFLAG  dump flag (Aug 7 2019, RK)
  //
  // Ouptuts in EVENT_SL struct:
  //  IDLENS        random integer identifier for galaxy lens
  //  ZLENS         redshift of lens galaxy
  //  LOGMASS_LENS  logmass of lens galaxy
  //  LOGMASS_ERR_LENS uncertainty on logmass
  //  NIMG      Number of images
  //  DELAY      list of NIMG time delays (days)
  //  MAGNIF     list of NIMG magnifications 
  //  XIMG_SRC   list of NIMG XSRC-LENS separations (arcsec)
  //  YIMG_SRC   list of NIMG YSRC-LENS separations (arcsec)
  //  XGAL_SRC   XGAL-LENS sep
  //  YGAL_SRC   
  //
  // July 1 2022: pass new output args LOGMASS[_ERR]_LENS
  //

  int    NLENS_LIB = INPUTS_STRONGLENS.NLENS;
  int    NIMG_local=0, img,i,j, numLens ;
  long long int IDLENS_local;
  double FlatRan, ZSRC, ZSRC_MINTOL, ZSRC_MAXTOL ;
  double zLENS_local, LOGMASS_local, LOGMASS_ERR_local, zSRC_local ;

  char fnam[] = "get_stronglens" ;

  // ---------------- BEGIN ---------------

  // always burn random
  FlatRan = getRan_Flat1(2);   // flat between 0 and 1

  
  EVENT_SL->NIMG = 0 ;
  if ( !INPUTS_STRONGLENS.USE_FLAG ) { return ; }

  INPUTS_STRONGLENS.NCALL++ ;

  ZSRC_MINTOL = 0.99 * zSN ;
  ZSRC_MAXTOL = 1.01 * zSN ;

  numLens = 0;
  for(i=0; i < NLENS_LIB; i++) {
    ZSRC = INPUTS_STRONGLENS.ZSRC[i] ;
    if ( ZSRC>=ZSRC_MINTOL && ZSRC<=ZSRC_MAXTOL ) { numLens++; }
  }

  if ( numLens==0 ) {
    // printf("WARNING: No lens in library matches zSRC= %f\n", zSN);
    return;
  }


  int possible_lenses[numLens];
  j = 0;
  for(i=0; i < NLENS_LIB; ++i) {
    ZSRC = INPUTS_STRONGLENS.ZSRC[i];
    if( ZSRC >= ZSRC_MINTOL && ZSRC <= ZSRC_MAXTOL) {
      possible_lenses[j] = i ;
      j++ ;
    }
  }

  int random_lens_index = possible_lenses[ (int)( FlatRan*(numLens-1) ) ];

  IDLENS_local  = INPUTS_STRONGLENS.IDLENS[random_lens_index];
  zLENS_local   = (double)INPUTS_STRONGLENS.ZLENS[random_lens_index];
  LOGMASS_local = (double)INPUTS_STRONGLENS.LOGMASS_LENS[random_lens_index];
  LOGMASS_ERR_local = (double)INPUTS_STRONGLENS.LOGMASS_ERR_LENS[random_lens_index];
  zSRC_local    = (double)INPUTS_STRONGLENS.ZSRC[random_lens_index];
  NIMG_local    = INPUTS_STRONGLENS.NIMG[random_lens_index];
  
  
  // strip off image-dependent quantities, and recast to double
  for(img=0 ; img < NIMG_local; img++ ) {

    EVENT_SL->XIMG_SRC_LIST[img]   = 
      (double)INPUTS_STRONGLENS.XIMG_SRC[random_lens_index][img]; 
    EVENT_SL->YIMG_SRC_LIST[img]   = 
      (double)INPUTS_STRONGLENS.YIMG_SRC[random_lens_index][img]; 

    EVENT_SL->MAGNIF_LIST[img]        = 
      (double)INPUTS_STRONGLENS.MAGNIF[random_lens_index][img]; 

    EVENT_SL->DELAY_LIST[img]      = 
      (double)INPUTS_STRONGLENS.DELAY[random_lens_index][img]; 
  }

  EVENT_SL->XGAL_SRC   = 
    (double)INPUTS_STRONGLENS.XGAL_SRC[random_lens_index]; 
  EVENT_SL->YGAL_SRC   = 
    (double)INPUTS_STRONGLENS.YGAL_SRC[random_lens_index]; 
  
  // load return argument scalars  
  EVENT_SL->IDLENS      = IDLENS_local ; 
  EVENT_SL->zLENS       = zLENS_local ;
  EVENT_SL->LOGMASS     = LOGMASS_local ;
  EVENT_SL->LOGMASS_ERR = LOGMASS_ERR_local ;
  EVENT_SL->NIMG        = NIMG_local;
  
  DUMPFLAG = 0 ;
  if ( DUMPFLAG ) {
    printf(" xxx \n");
    printf(" xxx -------- %s DUMP -------- \n", fnam );
    printf(" xxx input zSN = %.3f \n", zSN);
    printf(" xxx output IDLENS=%d at zLENS=%.3f and zSRC=%.3f \n", 
	   IDLENS_local, zLENS_local, zSRC_local );
    printf(" xxx output LOGMASS_LENS = %.3f +_ %.3f \n",
	   LOGMASS_local, LOGMASS_ERR_local );
    for(img=0 ; img < NIMG_local; img++ ) {
      printf(" xxx output image-%d: mu=%.3f deltaT=%.2f  X,Y-offset=%f,%f\n",
	     img, 
	     EVENT_SL->MAGNIF_LIST[img], 
	     EVENT_SL->DELAY_LIST[img], 
	     EVENT_SL->XIMG_SRC_LIST[img], EVENT_SL->YIMG_SRC_LIST[img] );
    }
    printf(" xxx \n");    fflush(stdout);
  }

  return ;

} // end get_stronglens


