/****************************************************

  Created Jun 12 2020 by R.Kessler

  Generic utility to draw random numbers from probability distribution
  functions (PDF) based on analytical function (asym Gauss), or 
  multi-dimensional map. The map can depend on redshift and any
  combination of HOSTLIB parameters.

  Initial use is for SALT2 parameters (c,x1,RV,AV).

  USE_SUBPROCESS is a pre-processor flag so that SALT2mu.c can
  include only init_genPDF (to read and store map), while not
  including the extra baggage from sntools_host. The simulation
  does NOT define USE_SUBPROCESS, and therefore all of the code
  below is used for the simulation.

  Oct 23 2020: use get_VAL_RANGE_genPDF() 
  Dec 23 2020: improve algorithmn in get_VAL_RANGE_genPDF
  Mar 26 2021: allow LOGparName in GENPDF maps; e.g., LOGEBV, LOGAV.
  May 26 2021: new function free_memory_genPDF() to free GRIDMAP memory.
  Aug 10 2021: check priority of GENPDF vs. GENGAUSS

  Mar 05 2024: minor refac; rename IDMAP var -> IMAP for sparse index;
               reserve IDMAP as absolute index. Refac is to avoid confusion.

  Jun 4 2025: enable reading GEN[PEAK,SIGMA,RANGE]_SALT2c so that asyn Gauss
              profile can go into genPDF map file.

 ****************************************************/

#ifndef USE_SUBPROCESS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"
#include "sntools_genPDF.h"
#include "sntools_host.h"

//#include "sntools_genGauss_asym.h"

#endif

// =====================================================
void init_genPDF(int OPTMASK, FILE *FP, char *fileName, char *ignoreList) {

  // June 2020
  // Parse map in fileName. 
  //
  // Inputs
  //   OPTMASK -> bit-mask of options
  //         1 : allow extrapolation outside range of map
  //               (e.g, LOGMASS in HOSTLIB extends beyone map)
  //         8 : use already opened FP
  //
  //   fileName -> name of file with map(s)
  //   ignoreList -> comma-separated list of map(s) to ignore 
  // 
  //
  // Inside fileName, each VARNAMES key corresponds to a separate map.
  // First variable after VARNAMES is the variable for which we wish 
  // to generate: VARGEN (e.g., SALT2x1, SALT2c, RV_HOST, etc ...)
  // The last variable should be PROB. Variables between VARGEN and 
  // PROB are are addition multi-dimensional dependences; can be 
  // any HOSTLIB variable.
  //
  // Examples of fileName contents:
  //  VARNAMES: SALT2c PROB  
  //     PROB is a 1D function of SALT2c (color)
  //
  //  VARNAMES: SALT2c LOGMASS PROB
  //     PROB is a 2D function of SALT2c and LOGMASS
  //
  //  VARNAMES: SALT2c ZTRUE LOGMASS PROB
  //     PROB is a 3D function of SALT2c, redshift and LOGMASS
  //     [note that REDSHIFT can be used instead of ZTRUE]
  //
  //  HISTORY:
  //   Jun 24 2021 Dillon - including alpha beta asym gauss
  //   Aug 10 2021 Brodie and Rick fixed bug setting NVAR properly
  //   Nov  4 2021 RK - check OPTMASK_GENPDF_KEYSOURCE_ARG
  //   Nov 10 2022 RK - skip parsing comment lines
  //   May 20 2024 RK - read MAG_OFFSET
  // -----------

  FILE *fp;
  int gzipFlag,  NDIM, NFUN, i, NWD=0 ;
  int NMAP=0, NVAR=0, NITEM=0, ivar, IDMAP=0;
  bool IGNORE_MAP ;
  int  MXWD_TMPLINE = 20;
  int  MXCHAR_TMPLINE = 800;
  int  MEMC_TMPLINE   = MXCHAR_TMPLINE * sizeof(char);

  char c_get[400], fileName_full[MXPATHLEN];
  char *LINE, *TMPLINE ;
  char *MAPNAME, *ptrVar;
  char KEY_ROW[]  = "PDF:", KEY_STOP[] = "", PATH[] = "" ;

  GENGAUSS_ASYM_DEF  gengauss_SALT2ALPHA;
  GENGAUSS_ASYM_DEF  gengauss_SALT2BETA;
  GENGAUSS_ASYM_DEF  gengauss_SALT2c;  // June 2025
  bool  USE_SALT2 = false;

  char **ptr_ITEMLIST;
  int KEYSOURCE = 1; // default source is from file

  char fnam[]     = "init_genPDF";
  
  // ------------- BEGIN -------------

  OPTMASK_GENPDF = OPTMASK ; // Dec 23 2020

  NMAP_GENPDF = NCALL_GENPDF = 0;
  MAG_OFFSET_GENPDF = 0.0 ;
  
  if ( IGNOREFILE(fileName) ) { return; }

  sprintf(BANNER,"%s with OPTMASK=%d", fnam, OPTMASK);
  print_banner(BANNER);

  // init optional asymGauss params for SALT2alpha and beta
  init_GENGAUSS_ASYM(&gengauss_SALT2ALPHA, 0.0 );
  init_GENGAUSS_ASYM(&gengauss_SALT2BETA,  0.0 );
  init_GENGAUSS_ASYM(&gengauss_SALT2c,     0.0 );

  ptr_ITEMLIST = (char**)malloc( MXWD_TMPLINE*sizeof(char*));
  for(i=0; i<MXWD_TMPLINE; i++) 
    { ptr_ITEMLIST[i] = (char*)malloc(80*sizeof(char)); }

  LINE    = (char*) malloc(MEMC_TMPLINE);     // allow for long comments
  TMPLINE = (char*) malloc(MEMC_TMPLINE+100); 

#ifndef USE_SUBPROCESS
  if ( HOSTLIB_WGTMAP.N_SNVAR > 0 ) {
    sprintf(c1err,"Found SN params in WGTMAP");
    sprintf(c2err,"-> not compatible with GENPDF_FILE");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
#endif

  // check OPTMASK options
  OPT_EXTRAP_GENPDF = 0 ;
  if ( (OPTMASK & OPTMASK_GENPDF_EXTRAP) > 0 ) { 
    printf("\t Enable extrapolation beyond map ranges.\n");
    OPT_EXTRAP_GENPDF = 1; 
  }

  if ( (OPTMASK & OPTMASK_GENPDF_SLOW) > 0 ) { 
    printf("\t skip speed-trick: select regardless of PROB. \n");
  }
  else {
    printf("\t use speed-trick: select from range bounded by PROB ~ 0. \n");
  }

  if ( (OPTMASK & OPTMASK_GENPDF_KEYSOURCE_ARG) > 0 ) {
    printf("\t GENPDF_FILE is from command line -> allow overrides.\n");
    KEYSOURCE = 2; // genpdf_file was from command line
  }
  
  if ( strlen(ignoreList) > 0 ) 
    { printf("\t Ignore PDF map(s): %s \n", ignoreList); }
  
  // - - - - - - - - - - - - - 
  // open file and read it

  if ( (OPTMASK & OPTMASK_GENPDF_EXTERNAL_FP)> 0 ) { 
    fp = FP; 
    printf("  Read already opened file,\n\t %s\n", fileName);
  }
  else  { 
    fp = snana_openTextFile(1, PATH, fileName, fileName_full, &gzipFlag); 
  }


  if ( !fp ) {
    sprintf(c1err,"Could not open GENPDF_FILE");
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  bool IS_VARNAMES, IS_SALT2, IS_SIM ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    // check for asymmetric gaussian for alpha,beta
    //if ( strstr(c_get,"SALT2") != NULL ) {  // SALT2 is in c_get
    fgets(LINE,200,fp);

    if ( commentchar(c_get) ) { continue; }

    sprintf(TMPLINE,"%s %s", c_get, LINE);
    splitString(TMPLINE, " ", fnam, MXWD_TMPLINE,  // inputs             
		&NITEM, ptr_ITEMLIST );     // outputs

    IS_VARNAMES = (strcmp(c_get,"VARNAMES:") == 0 );
    IS_SALT2    = (strstr(c_get,"SALT2")     != NULL );
    IS_SIM      = (strstr(c_get,"SIM_")      != NULL );

    if ( strcmp(c_get,"MAG_OFFSET:") == 0 ) {
      sscanf(ptr_ITEMLIST[1], "%le", &MAG_OFFSET_GENPDF);
      printf("    Load MAG_OFFSET = %.3f \n", MAG_OFFSET_GENPDF);
      fflush(stdout);
    }
    // try sim-input keys
    if ( IS_SALT2 ) {

      NWD = parse_input_GENGAUSS("SALT2ALPHA", ptr_ITEMLIST, KEYSOURCE, 
				 &gengauss_SALT2ALPHA);      
      if ( NWD > 0 ) { USE_SALT2 = true; }

      NWD = parse_input_GENGAUSS("SALT2BETA", ptr_ITEMLIST, KEYSOURCE, 
				 &gengauss_SALT2BETA);
      if ( NWD > 0 ) { USE_SALT2 = true; }

      NWD = parse_input_GENGAUSS("SALT2c", ptr_ITEMLIST, KEYSOURCE, 
				 &gengauss_SALT2c);      
      if ( NWD > 0 ) { USE_SALT2 = true; }

    }

    // try column name in FITRES file
    if ( IS_SIM ) {
      NWD = parse_input_GENGAUSS("SIM_beta", ptr_ITEMLIST, KEYSOURCE, 
				 &gengauss_SALT2BETA);
      NWD = parse_input_GENGAUSS("SIM_alpha", ptr_ITEMLIST, KEYSOURCE, 
				 &gengauss_SALT2ALPHA);      
    }

    if ( IS_VARNAMES ) {

      NVAR = NITEM - 1;  // avoid VARNAMES key
      GENPDF[NMAP].NVAR = NVAR;
      for(ivar=0; ivar < NVAR; ivar++ )	{ 
	ptrVar = ptr_ITEMLIST[ivar+1]; 
	assign_VARNAME_GENPDF(NMAP, ivar, ptrVar);
	GENPDF[NMAP].IVAR_HOSTLIB[ivar] = -9 ; // init for below
      }

      // check options to abort or ignore map(s)
      ptrVar = GENPDF[NMAP].VARNAMES[0] ;
      checkAbort_VARNAME_GENPDF(ptrVar);
      IGNORE_MAP = false;
      if ( strlen(ignoreList) > 0 ) 
	{ if ( strstr(ignoreList,ptrVar) != NULL ) { IGNORE_MAP=true;}  }
      if ( IGNORE_MAP ) { 
	printf("\t IGNORE %s \n", GENPDF[NMAP].MAPNAME );
	continue; 
      }

      IDMAP  = IDGRIDMAP_GENPDF + NMAP ;  // absolute ID
      MAPNAME = GENPDF[NMAP].MAPNAME;
      NFUN   = 1  ; // for now, assume only one PROB column
      NDIM   = NVAR - NFUN ;
      read_GRIDMAP(fp, MAPNAME, KEY_ROW, KEY_STOP, 
		   IDMAP, NDIM, NFUN, OPT_EXTRAP_GENPDF, 
		   MXROW_GENPDF, fnam, &GENPDF[NMAP].GRIDMAP );

      GENPDF[NMAP].N_CALL      = 0 ;
      GENPDF[NMAP].N_ITER_SUM  = 0 ;
      GENPDF[NMAP].N_ITER_MAX  = 0 ;
      GENPDF[NMAP].PROB_EXPON_REWGT = 1.0 ; // default is no rewgt 
      /*
      int NROW = GENPDF[NMAP].GRIDMAP.NROW;
      char *VARLIST = GENPDF[NMAP].GRIDMAP.VARLIST ;
      printf(" Found PROB(%s)  NDIM=%d, NROW=%d \n",  VARLIST, NDIM, NROW);
      */

      NMAP++ ;
    }
  }

  // - - - - - - -
  if ( USE_SALT2 ) {

    if ( gengauss_SALT2ALPHA.USE ) {
      init_genPDF_from_GenGauss(NMAP, &gengauss_SALT2ALPHA) ; 
      NMAP++ ; 
    }
    if ( gengauss_SALT2BETA.USE ) {
      init_genPDF_from_GenGauss(NMAP, &gengauss_SALT2BETA) ;
      NMAP++ ;
    }
    if ( gengauss_SALT2c.USE ) {
      init_genPDF_from_GenGauss(NMAP, &gengauss_SALT2c) ;
      NMAP++ ;
    }

    /*
    dump_GENGAUSS_ASYM(&gengauss_SALT2ALPHA);
    dump_GENGAUSS_ASYM(&gengauss_SALT2BETA);
    debugexit(fnam); // xxx remove me!
    */
  } // end USE_SALT2
  

  NMAP_GENPDF = NMAP ;

#ifndef USE_SUBPROCESS
  // - - - - - - - -
  // loop thru maps again and check that extra variables (after 1st column)
  // exist in HOSTLIB
  bool IS_LOGPARAM;
  int  ivar_hostlib, imap, imap_tmp, NVAR_HOSTLIB=0 ;
  int  ABORTFLAG = 0 ;
  char *VARNAME;

  printf("\n   Check for GENPDF variables in HOSTLIB:\n"); fflush(stdout);
  
  for(imap=0; imap < NMAP; imap++ ) {
    NVAR = GENPDF[imap].NVAR ;

    if ( GENPDF[imap].GRIDMAP.NDIM == 1 ) { continue; }

    NVAR_HOSTLIB = 0 ;
    for(ivar=1; ivar < NVAR-1; ivar++ )	{ 
      VARNAME = GENPDF[imap].VARNAMES[ivar] ;
      checkAlternateVarNames_HOSTLIB(VARNAME);
      ivar_hostlib = IVAR_HOSTLIB(VARNAME,ABORTFLAG);
      if ( ivar_hostlib < 0 ) {
	sprintf(c1err,"Could not find HOSTLIB variable '%s'", VARNAME);
	sprintf(c2err,"Check HOSTLIB and check GENPDF_FILE");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      printf("\t Found HOSTLIB IVAR=%2d for VARNAME='%s' (%s) \n",
	     ivar_hostlib, VARNAME, GENPDF[imap].MAPNAME ); fflush(stdout);
      GENPDF[imap].IVAR_HOSTLIB[ivar] = ivar_hostlib;
      NVAR_HOSTLIB++ ;
    }
    printf("    Found %d HOSTLIB variables in GENPDF map %s\n",
	   NVAR_HOSTLIB, GENPDF[imap].MAPNAME);    fflush(stdout);
  }

#endif

  // free memory
  for(i=0; i<MXWD_TMPLINE; i++)   { free(ptr_ITEMLIST[i]); }
  free(ptr_ITEMLIST);
  free(LINE); free(TMPLINE);

  //  debugexit(fnam);
  return;

} // end init_genPDF

// =======================================
void init_genPDF_from_GenGauss(int IMAP, GENGAUSS_ASYM_DEF *GENGAUSS) {

  // Created July 16 2021
  // Compute 1D GRIDMAP from GENGAUSS and load GENPDF[IMAP].GRIDMAP
  // as if it were read from GENPDF map file.
  // Initial use for SALT2 alpha/beta in SALT2mu SUBPROCESS
  // WARNING! Currently only works for 1D maps
  //
  // Mar 5 2024: reduce RANGE if more than 10 sigma away from peak;
  //   needed to avoid downstream abort from not finding prob > 0.001
  //   in course-grid search.

  // May 2024: avoid sigma = 0

  double PEAK    = GENGAUSS->PEAK ;
    
  double sig_tiny = PEAK/1.0E6;
  if ( GENGAUSS->SIGMA[0] == 0.0 )
    { GENGAUSS->SIGMA[0] = sig_tiny;  GENGAUSS->RANGE[0] = PEAK - 2.0*sig_tiny; }
  if ( GENGAUSS->SIGMA[1] == 0.0 )
    { GENGAUSS->SIGMA[1] = sig_tiny;  GENGAUSS->RANGE[1] = PEAK + 2.0*sig_tiny; }

  char   *NAME   = GENGAUSS->NAME ;
  double *RANGE  = GENGAUSS->RANGE;
  double siglo   = GENGAUSS->SIGMA[0];
  double sighi   = GENGAUSS->SIGMA[1];
  double sigavg  = 0.5*(siglo+sighi);

  int   IDMAP  = IDGRIDMAP_GENPDF + IMAP ;
  int   NBIN_PER_SIGMA = 20; // hard-wired guess
  double NSIGMA_MAX    = 8.0; // range must not be more than this away from peak
  
  int   NBIN;
  double XNBIN, xnsig_tmp, x_orig;
  char fnam[] = "init_genPDF_from_GenGauss" ; 

  // ---------- BEGIN ----------
 
  // if RANGE is to wide compared to sigma, reduce RANGE so that
  // downstream selection does not fail in finding a hard-to-find
  // delta-function 
  xnsig_tmp = (PEAK - RANGE[0]) / siglo;
  if ( xnsig_tmp > NSIGMA_MAX ) {
    x_orig = RANGE[0] ; // for comment
    RANGE[0] = (PEAK - NSIGMA_MAX*siglo) ;
    printf("  %s: [PEAK-min(%s)]/siglo=%.1f; change RANGE[0]=%.3f to %.3f\n", 
	   fnam, NAME, xnsig_tmp,   x_orig, RANGE[0] );
  }

  xnsig_tmp = (RANGE[1] - PEAK ) / sighi;
  if ( xnsig_tmp > NSIGMA_MAX ) {
    x_orig = RANGE[1] ; // for comment
    RANGE[1] = (PEAK + NSIGMA_MAX*sighi) ;
    printf("  %s: [max(%s)-PEAK]/sighi=%.1f; change RANGE[1]=%.3f to %.3f\n", 
	   fnam, NAME, xnsig_tmp, x_orig, RANGE[1] );
  }

  // - - - -


  XNBIN = (float)NBIN_PER_SIGMA * (RANGE[1] - RANGE[0]) / sigavg ;
  NBIN  = (int)(XNBIN+0.5);

  // call utility in sntools_genGauss_asym.c
  compute_genGauss_GRIDMAP(GENGAUSS, NAME, IDMAP, OPT_EXTRAP_GENPDF, 
			   NBIN, RANGE, fnam,
			   &GENPDF[IMAP].GRIDMAP ); // <== returned

  GENPDF[IMAP].PROB_EXPON_REWGT = 1.0 ; // Mar 14, 2024

  assign_VARNAME_GENPDF(IMAP, 0, NAME );

  bool debugflag = false ;
  if (debugflag) {
    int nbin_test = 4;
    double binsize = (RANGE[1] - RANGE[0]) / (float)nbin_test ;
    double x, funVal;
    for (x = RANGE[0]; x <= RANGE[1]; x+= binsize) {
      interp_GRIDMAP(&GENPDF[IMAP].GRIDMAP, &x, &funVal) ;
      printf("xxx %s: %s=%f ----> funVal = %f \n", fnam, NAME, x, funVal);
    }
    debugexit(fnam);
  }

  return;

} // end init_genPDF_from_GenGauss


// =======================================
void free_memory_genPDF(void) {

  // Created May 26 2021
  // Use malloc_GRIDMAP utility to free memory of all genPDF maps.
  // E.g., useful for using SALT2mu as iterative SUBPROCESS

  int  ivar, NVAR, NDIM, imap, NMAP = NMAP_GENPDF; 
  int  NFUN = 1;
  int  MAPSIZE = MXROW_GENPDF ;
  char fnam[] = "free_memory_genPDF" ;

  // --------- BEGIN --------

  for ( imap=0; imap < NMAP; imap++ )  {
    NDIM = GENPDF[imap].GRIDMAP.NDIM ;
    NVAR = GENPDF[imap].NVAR;

    malloc_GRIDMAP(-1, &GENPDF[imap].GRIDMAP, NFUN, NDIM, MAPSIZE);
    for(ivar=0; ivar<NVAR; ivar++ )  { free(GENPDF[imap].VARNAMES[ivar]); }
  }

  return;
} // end free_memory_genPDF

// =======================================
void assign_VARNAME_GENPDF(int imap, int ivar, char *varName) {

  char *MAPNAME = GENPDF[imap].MAPNAME ;
  char *VARLIST = GENPDF[imap].GRIDMAP.VARLIST;
  //  char fnam[] = "assign_VARNAME_GENPDF";
  // ---------- BEGIN ------------

  if ( ivar == 0 )  { 
    sprintf(MAPNAME, "%s_PDF", varName);
    sprintf(VARLIST,"%s", varName);
  } 
  else {
    strcat(VARLIST,",");    strcat(VARLIST,varName);
  }

  GENPDF[imap].VARNAMES[ivar] = (char*) malloc(40*sizeof(char));
  sprintf(GENPDF[imap].VARNAMES[ivar], "%s", varName);
  
  return ;

} // end assign_VARNAME_GENPDF


// ===================================
void checkAbort_VARNAME_GENPDF(char *varName) {

  // Created Aug 11 2021
  // Abort on easy/common mistake

  char fnam[] = "checkAbort_VARNAME_GENPDF" ;
  if ( strcmp(varName,"c")==0 || strcmp(varName,"x1") == 0 ) {
    sprintf(c1err,"Invalid GENPDF variable '%s'", varName);
    sprintf(c2err,"Try 'SALT2%s' instead", varName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return;

} // end checkAbort_VARNAME_GENPDF

#ifndef USE_SUBPROCESS

// ===================================
double funVal_genPDF(char *parName, double x, GENGAUSS_ASYM_DEF *GENGAUSS) {

  // Created Dec 19 2023
  // For input parName and value x, return genPDF function value.
  // For multi-dim genPDF map, additional HOSTLIB-dependent values
  // are internally included.
  //

  int    IMAP = -9 ;
  double prob  = -9.0 ;
  bool   IS_LOGPARAM = false ;
  char fnam[] = "funVal_genPDF" ;

  // -------------- BEGIN ------------

  if ( NMAP_GENPDF > 0 ) {
    int    istat, NDIM, ivar, IVAR_HOSTLIB;
    double xval[MXVAR_GENPDF], EXPON_REWGT ;
    IMAP        = IMAP_GENPDF(parName, &IS_LOGPARAM) ;
    NDIM        = GENPDF[IMAP].GRIDMAP.NDIM ;
    EXPON_REWGT = GENPDF[IMAP].PROB_EXPON_REWGT ;

    xval[0] = x;
    for(ivar=1; ivar < NDIM ; ivar++ ) {
      IVAR_HOSTLIB   = GENPDF[IMAP].IVAR_HOSTLIB[ivar];
      xval[ivar]     = get_VALUE_HOSTLIB(IVAR_HOSTLIB, SNHOSTGAL.IGAL);
    }

    istat = interp_GRIDMAP(&GENPDF[IMAP].GRIDMAP, xval, &prob);

    if ( EXPON_REWGT != 1.0 )  { prob = pow(prob,EXPON_REWGT); }
  }

  if  ( GENGAUSS->USE ) {
    prob = funVal_GENGAUSS_ASYM(x, GENGAUSS);
  }

  if ( prob < -1.0E-6 || prob > 1.000001 ) {
    sprintf(c1err,"Invalid prob = %12.9f for x=%f  parName=%s (IMAP=%d)", 
	    prob, x, parName, IMAP);
    sprintf(c2err,"Need either GENPDF map or GENGAUSS");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return prob ;

} // end funVal_genPDF

// =====================================================
double getRan_genPDF(char *parName, GENGAUSS_ASYM_DEF *GENGAUSS) {
	
  // Return random number from GENPDF map corresponding to input *parName.
  // If *parName does not match a GENPDF map, use input *GENGAUSS instead.
  // 
  // Dec 2023: implement PROB_EXPON_REWGT

  int    KEYSOURCE_GENGAUSS = GENGAUSS->KEYSOURCE ;
  int    IGAL               = SNHOSTGAL.IGAL;

  int    ILIST_RAN = 1;
  int    N_ITER=0, MAX_ITER  = MXITER_GENPDF ;
  int    N_EVAL = 0, IMAP, ivar, NDIM, istat, itmp, IVAR_HOSTLIB;
  double val_inputs[MXVAR_GENPDF], prob_ref, prob, r = 0.0 ;
  double VAL_RANGE[2], FUNMAX, EXPON_REWGT, prob_ratio ;
  int    LDMP = 0 ;
  bool   DO_GENGAUSS; 
  bool   IS_LOGPARAM = false ; // true -> param stored as LOGparam
  char   *MAPNAME, *VARNAME ;
  char fnam[] = "getRan_genPDF";
  
  // ------------- BEGIN -----------

  // check for map in GENPDF_FILE argument of sim-input file
  if ( NMAP_GENPDF > 0 ) {
    IMAP = IMAP_GENPDF(parName, &IS_LOGPARAM);
    if ( IMAP >= 0 ) {
      N_EVAL++ ; NCALL_GENPDF++ ;
      MAPNAME       = GENPDF[IMAP].MAPNAME ;
      VARNAME       = GENPDF[IMAP].VARNAMES[0] ;
      FUNMAX        = GENPDF[IMAP].GRIDMAP.FUNMAX[0] ;
      NDIM          = GENPDF[IMAP].GRIDMAP.NDIM ;
      EXPON_REWGT   = GENPDF[IMAP].PROB_EXPON_REWGT ; 
      if ( EXPON_REWGT != 1.0 ) { FUNMAX = pow(FUNMAX,EXPON_REWGT); }
      prob_ref=1.0; prob=0.0;

      // tack on optional dependence on HOSTLIB
      // Leave var_inputs[0] to be filled below inside while loop
      for(ivar=1; ivar < NDIM; ivar++ ) {
	IVAR_HOSTLIB     = GENPDF[IMAP].IVAR_HOSTLIB[ivar];
	val_inputs[ivar] = get_VALUE_HOSTLIB(IVAR_HOSTLIB,IGAL);
      }

      // get min/max VALUE range for random selection;
      // function here returns VAL_RANGE

      if ( IS_LOGPARAM ) {
	VAL_RANGE[0]  = GENPDF[IMAP].GRIDMAP.VALMIN[0] ;
	VAL_RANGE[1]  = GENPDF[IMAP].GRIDMAP.VALMAX[0] ;
      }
      else {
	get_VAL_RANGE_genPDF(IMAP, val_inputs, VAL_RANGE, 0 );
      }

      // - - - - - -
      LDMP = 0; // (NCALL_GENPDF < 5 );
      if ( LDMP ) {
	printf(" xxx \n");
	printf(" xxx %s -------- (%.3f < %s < %.3f )--------------- \n", 
	       fnam, VAL_RANGE[0], parName, VAL_RANGE[1] );	       
      }

      while ( prob < prob_ref ) {
	prob_ref      = getRan_Flat1(ILIST_RAN); 
	r             = getRan_Flat(ILIST_RAN, VAL_RANGE);
	val_inputs[0] = r ;  

	istat = interp_GRIDMAP(&GENPDF[IMAP].GRIDMAP, val_inputs, &prob);
	if ( EXPON_REWGT != 1.0 )  { prob = pow(prob,EXPON_REWGT); }

	if ( istat != SUCCESS ) {
	  print_preAbort_banner(fnam);
	  for(ivar=0; ivar < NDIM; ivar++ ) {
	    printf("   %s = %f  [mapRange: %.2f to %.2f]\n", 
		   GENPDF[IMAP].VARNAMES[ivar], val_inputs[ivar],
		   GENPDF[IMAP].GRIDMAP.VALMIN[ivar],
		   GENPDF[IMAP].GRIDMAP.VALMAX[ivar]	 );
	  }

	  sprintf(c1err,"interp_GRIDMAP returned istat=%d", istat);
	  sprintf(c2err,"Value probably outside GENPDF map range");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	} // end abort check

	prob /= FUNMAX; // normalize to max prob = 1.0

	if ( LDMP ) {
	  printf(" xxx %s: %s=%8.4f, m=%6.2f  prob(PDF,ref)=%f,%f \n",
		 fnam, parName, r, val_inputs[1],  prob, prob_ref); 
	  fflush(stdout);
	}

	TMPSTORE_PROB_REF_GENPDF[N_ITER] = (float)prob_ref;
	TMPSTORE_PROB_GENPDF[N_ITER]     = (float)prob;
	TMPSTORE_RAN_GENPDF[N_ITER]      = (float)r ;
	N_ITER++ ;
	if ( N_ITER >= MAX_ITER ) {
	  print_preAbort_banner(fnam);
	  printf("   %s(%s) \n", MAPNAME, GENPDF[IMAP].GRIDMAP.VARLIST );
	  
	  for(ivar=1; ivar < NDIM; ivar++ ) {
	    IVAR_HOSTLIB = GENPDF[IMAP].IVAR_HOSTLIB[ivar];
	    val_inputs[ivar] = get_VALUE_HOSTLIB(IVAR_HOSTLIB,IGAL);
	    printf("\t %s = %f \n", 
		   HOSTLIB.VARNAME_STORE[IVAR_HOSTLIB], val_inputs[ivar] );
	  }

	  for(itmp = 0; itmp < N_ITER; itmp++ ) {
	    r          = TMPSTORE_RAN_GENPDF[itmp] ;
	    prob       = TMPSTORE_PROB_GENPDF[itmp];
	    prob_ref   = TMPSTORE_PROB_REF_GENPDF[itmp];
	    prob_ratio = prob/prob_ref;
	    printf("   %4d: PROB=%8.5f, PROB_REF=%8.5f, "
		   "RATIO=%.4f, %s=%8.5f \n", 
		   itmp, prob, prob_ref, prob_ratio, parName, r);
	  }	  
	  get_VAL_RANGE_genPDF(IMAP, val_inputs, VAL_RANGE, 1 );
	  sprintf(c1err,"N_ITER=%d exceeds bound (prob_ref=%f)", 
		  N_ITER, prob_ref );
	  sprintf(c2err,"Check %s or increase MXITER_GENPDF", MAPNAME );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	} // end abort check

      } // end while loop over prob

    } // end IDMAP >= 0

    // track N_ITER stats
    GENPDF[IMAP].N_CALL++ ;
    GENPDF[IMAP].N_ITER_SUM += (double)N_ITER;
    if ( N_ITER > GENPDF[IMAP].N_ITER_MAX ) 
      { GENPDF[IMAP].N_ITER_MAX = N_ITER; }
    
  } // end genPDF 

  // - - - - - -
  // check explicit asymGauss function 
  DO_GENGAUSS = GENGAUSS->USE ;

  if  ( GENGAUSS->USE ) {
    N_EVAL++ ;
    r = getRan_GENGAUSS_ASYM(GENGAUSS) ;
  }

  // - - - - - -
  if ( N_EVAL != 1 ) {
    sprintf(c1err,"Found %d PDF options for parName='%s'", N_EVAL, parName);
    sprintf(c2err,"Requires 1 PDF option; check AsymGauss and GENPDF_FILE.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
 
  if ( LDMP ) 
    { printf(" xxx %s: return %s = %f \n", fnam, parName, r); }

  if ( IS_LOGPARAM ) { r = pow(TEN,r); }

  return(r);

} // end getRan_genPDF



// ========================================================
void get_VAL_RANGE_genPDF(int IMAP, double *val_inputs, 
			  double *VAL_RANGE, int dumpFlag ) {

  // Created Oct 23 2020
  // Return VAL_RANGE for random selection.
  // Default is the extreme min & max of map.
  // For better efficiency, however, this function checks the
  // val_inputs[0] range for which PROB > 0, and returns
  //
  //  VAL_RANGE[0] = min(val_inputs[0]) for which PROB > 0
  //  VAL_RANGE[1] = max(val_inputs[0]) for which PROB > 0
  //
  // Inputs:
  //   IMAP  : sparse integer ID of GRIDMAP
  //   val_inputs : includes HOSTLIB and PROB elements 1 to NDIM
  //                 val_inputs[0] is varied to see where PROB>0.
  //
  // The full VALUE range is divided int NBIN_CHECKPROB0 bins,
  // and min,max VALUE is stored for which PROB > 0.
  //
  // Dec 23 2020: 
  //   + add and subtract 1*BINSIZE to range.
  //   + check OPTMASK for SLOW option using full range.
  //

  int    NBIN_CHECKPROB0 = 10 ;
  double XNBIN = (double)NBIN_CHECKPROB0;
  
  double VAL_RANGE_PROB[2], VAL_TMP, VAL_BINSIZE, prob ;
  int itmp, istat;
  char fnam[] = "get_VAL_RANGE_genPDF" ;

  // ----------- BEGIN ------------

  VAL_RANGE[0]  = GENPDF[IMAP].GRIDMAP.VALMIN[0] ;
  VAL_RANGE[1]  = GENPDF[IMAP].GRIDMAP.VALMAX[0] ;

  if ( OPTMASK_GENPDF & OPTMASK_GENPDF_SLOW ) { return; }

  VAL_RANGE_PROB[0] = 9.0E12;  VAL_RANGE_PROB[1] = -9.0E12 ;
  VAL_BINSIZE = (VAL_RANGE[1] - VAL_RANGE[0]) / XNBIN; 
  for(itmp=0; itmp <= NBIN_CHECKPROB0; itmp++ ) {
    VAL_TMP       = VAL_RANGE[0] + VAL_BINSIZE*(double)itmp ;
    val_inputs[0] = VAL_TMP;
    istat = interp_GRIDMAP(&GENPDF[IMAP].GRIDMAP, val_inputs, &prob);

    if ( dumpFlag ) {
      printf("\t %s: prob(%.3f) = %le \n", fnam, VAL_TMP, prob);
      fflush(stdout);
    }

    if ( prob > PROBMAX_REJECT_GENPDF ) {
      // increment min range only once
      if (VAL_RANGE_PROB[0] > 8.0E12 )  { VAL_RANGE_PROB[0] = VAL_TMP; }

      // always increment max range
      VAL_RANGE_PROB[1] = VAL_TMP ; 
    }      
  } // end itmp

  // - - - - -
  if ( VAL_RANGE_PROB[0] > 8.0E12 || VAL_RANGE_PROB[1] < -8.0E12 ) {
    
    print_preAbort_banner(fnam);
    printf("\t IMAP(sparse)=%d  IDMAP=%d \n", IMAP, GENPDF[IMAP].GRIDMAP.ID);
    printf("\t NBIN_CHECKPROB0 = %d \n", NBIN_CHECKPROB0);
    printf("\t PROBMAX_REJECT_GENPDF = %le \n", PROBMAX_REJECT_GENPDF);
    sprintf(c1err,"Invalid VAL_RANGE_PROB (%.3le to %.3le) for %s=%.3f",
            VAL_RANGE_PROB[0], VAL_RANGE_PROB[1], val_inputs[0] );

    sprintf(c2err,"Consider GENPDF_OPTMASK to extrapolate.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // Dec 23 2020: add and subtract one binsize  
  VAL_RANGE[0] = VAL_RANGE_PROB[0] - VAL_BINSIZE ;
  VAL_RANGE[1] = VAL_RANGE_PROB[1] + VAL_BINSIZE ;

  // avoid going past boundaries            
  if ( VAL_RANGE[0] < GENPDF[IMAP].GRIDMAP.VALMIN[0] )
    { VAL_RANGE[0] = GENPDF[IMAP].GRIDMAP.VALMIN[0]; }

  if ( VAL_RANGE[1] > GENPDF[IMAP].GRIDMAP.VALMAX[0] )
    { VAL_RANGE[1] = GENPDF[IMAP].GRIDMAP.VALMAX[0]; }

  return;
} // end get_VAL_RANGE_genPDF

// ==========================================================
int IMAP_GENPDF(char *parName, bool *FOUND_LOGPARAM) {

  // return sparse IMAP for this input parName.
  // Match against first varName in each map.
  // Mar 26 2021: if LOGparName exists, return LOGPARAM=True.

  int IMAP_NONE =-9, imap;
  char *tmpName, LOGparName[60];

  *FOUND_LOGPARAM = false; 
  sprintf(LOGparName,"LOG%s", parName);  // Mar 26 2021
  for(imap=0; imap < NMAP_GENPDF; imap++ ) {
    tmpName = GENPDF[imap].VARNAMES[0];

    if ( strcmp(tmpName,parName)==0 ) 
      { return(imap); }

    if ( strcmp(tmpName,LOGparName)==0 ) 
      { *FOUND_LOGPARAM=true; return(imap); }
  }

  return(IMAP_NONE);

} // end IMAP_GENPDF

// =============================
void iter_summary_genPDF(void) {

  // Created Oct 23 2020
  // for each genPDF map, summarize mean number of iteration,
  // and max number.

  int imap, N_CALL, N_ITER_SUM, N_ITER_MAX ;
  double SUM, XNCALL, MEAN ;
  char *MAPNAME ;

  if ( NMAP_GENPDF == 0 ) { return; }

  printf("\n   N_ITER Summary for GENPDF_FILE:\n");

  for ( imap=0; imap < NMAP_GENPDF; imap++ ) {
    MAPNAME     = GENPDF[imap].MAPNAME ;
    N_ITER_MAX  = GENPDF[imap].N_ITER_MAX ;
    N_CALL      = GENPDF[imap].N_CALL ;
    N_ITER_SUM  = GENPDF[imap].N_ITER_SUM ;
    if ( N_CALL > 0 ) 
      { MEAN        = (double)N_ITER_SUM / (double)N_CALL ; }
    else
      { MEAN        = 0.0 ; }

    printf("\t%-16.16s : N_ITER[MEAN,MAX] = %.1f, %d \n",
	   MAPNAME, MEAN, N_ITER_MAX );

    fflush(stdout);
  }

} // end iter_summary_genPDF

#endif

