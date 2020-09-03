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
  // -----------

  FILE *fp;
  int gzipFlag, NMAP=0, NVAR=0, ivar, imap=-9, IDMAP=0, NDIM, NFUN ;
  bool IGNORE_MAP ;
  char c_get[60], fileName_full[MXPATHLEN], LINE[100], varName[60];
  char *MAPNAME, *ptrVar;
  char KEY_ROW[]  = "PDF:", KEY_STOP[] = "", PATH[] = "" ;
  char fnam[]     = "init_genPDF";

  // ------------- BEGIN -------------

  NMAP_GENPDF = NCALL_GENPDF = 0;
  if ( IGNOREFILE(fileName) ) { return; }

  print_banner(fnam);

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

  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"VARNAMES:") == 0 ) {
      fgets(LINE,100,fp); // scoop up variable names
      NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
      for(ivar=0; ivar < NVAR; ivar++ )	{ 
	get_PARSE_WORD(0,ivar,varName); 
	assign_VARNAME_GENPDF(NMAP,ivar,varName);
	GENPDF[NMAP].IVAR_HOSTLIB[ivar] = -9 ; // init for below
      }


      // check options to rejet or ignore map(s)
      ptrVar = GENPDF[NMAP].VARNAMES[0] ;
      IGNORE_MAP = false;
      if ( strlen(ignoreList) > 0 ) 
	{ if ( strstr(ignoreList,ptrVar) != NULL ) { IGNORE_MAP=true;}  }
      if ( IGNORE_MAP ) { 
	printf("\t IGNORE %s \n", GENPDF[NMAP].MAPNAME );
	continue; 
      }

      IDMAP  = IDGRIDMAP_GENPDF + NMAP ;
      MAPNAME = GENPDF[NMAP].MAPNAME;
      NFUN   = 1  ; // for now, assume only one PROB column
      NDIM   = NVAR - NFUN ;
      read_GRIDMAP(fp, MAPNAME, KEY_ROW, KEY_STOP, 
		   IDMAP, NDIM, NFUN, OPT_EXTRAP_GENPDF, 
		   MXROW_GENPDF, fnam, &GENPDF[NMAP].GRIDMAP );

      /*
      int NROW = GENPDF[NMAP].GRIDMAP.NROW;
      char *VARLIST = GENPDF[NMAP].GRIDMAP.VARLIST ;
      printf(" Found PROB(%s)  NDIM=%d, NROW=%d \n",  VARLIST, NDIM, NROW);
      */

      NMAP++ ;
    }
  }

  NMAP_GENPDF = NMAP ;

#ifndef USE_SUBPROCESS
  // - - - - - - - -
  // loop thru maps again and check that extra variables (after 1st column)
  // exist in HOSTLIB
  int ivar_hostlib ;
  int ABORTFLAG = 0 ;
  char *VARNAME;
  for(imap=0; imap < NMAP; imap++ ) {

    if ( GENPDF[imap].GRIDMAP.NDIM == 1 ) { continue; }

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
	     ivar_hostlib, VARNAME, GENPDF[imap].MAPNAME );
      GENPDF[imap].IVAR_HOSTLIB[ivar] = ivar_hostlib;
    }
  }
#endif

  //  debugexit(fnam);
  return;

} // end init_genPDF


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


#ifndef USE_SUBPROCESS

// =====================================================
double get_random_genPDF(char *parName, GENGAUSS_ASYM_DEF *GENGAUSS) {
			 
  int    ILIST_RAN = 1;
  int    N_ITER=0, MAX_ITER  = MXITER_GENPDF ;
  int    N_EVAL = 0, IDMAP, ivar, NDIM, istat, itmp ;
  int    IVAR_HOSTLIB, IGAL = SNHOSTGAL.IGAL;
  double val_inputs[MXVAR_GENPDF], prob_ref, prob, r = 0.0 ;
  double PAR_RANGE[2], FUNMAX, prob_ratio ;
  int    LDMP = 0 ;
  char   *MAPNAME ;
  char fnam[] = "get_random_genPDF";
  
  // ------------- BEGIN -----------

  // check for map in GENPDF_FILE argument of sim-input file
  if ( NMAP_GENPDF > 0 ) {
    IDMAP = IDMAP_GENPDF(parName);
    if ( IDMAP >= 0 ) {
      N_EVAL++ ; NCALL_GENPDF++ ;
      MAPNAME       = GENPDF[IDMAP].MAPNAME ;
      PAR_RANGE[0]  = GENPDF[IDMAP].GRIDMAP.VALMIN[0] ;
      PAR_RANGE[1]  = GENPDF[IDMAP].GRIDMAP.VALMAX[0] ;
      FUNMAX        = GENPDF[IDMAP].GRIDMAP.FUNMAX[0] ;
      NDIM          = GENPDF[IDMAP].GRIDMAP.NDIM ;
      prob_ref=1.0; prob=0.0;

      // LDMP = (NCALL_GENPDF > 10 && NCALL_GENPDF < 15 );
      if ( LDMP ) {
	printf(" xxx \n");
	printf(" xxx %s -------- (%.3f < %s < %.3f )--------------- \n", 
	       fnam, PAR_RANGE[0], parName, PAR_RANGE[1] );	       
      }

      while ( prob < prob_ref ) {
	prob_ref      = FlatRan1(ILIST_RAN); 
	r             = FlatRan( ILIST_RAN, PAR_RANGE);
	val_inputs[0] = r ;  

	// tack on optional dependence on HOSTLIB
	for(ivar=1; ivar < NDIM; ivar++ ) {
	  IVAR_HOSTLIB = GENPDF[IDMAP].IVAR_HOSTLIB[ivar];
	  val_inputs[ivar] = get_VALUE_HOSTLIB(IVAR_HOSTLIB,IGAL);
	}

	istat = interp_GRIDMAP(&GENPDF[IDMAP].GRIDMAP, val_inputs, &prob);
	if ( istat != SUCCESS ) {
	  print_preAbort_banner(fnam);
	  for(ivar=0; ivar < NDIM; ivar++ ) {
	    printf("   %s = %f  [mapRange: %.2f to %.2f]\n", 
		   GENPDF[IDMAP].VARNAMES[ivar], val_inputs[ivar],
		   GENPDF[IDMAP].GRIDMAP.VALMIN[ivar],
		   GENPDF[IDMAP].GRIDMAP.VALMAX[ivar]	 );
	  }

	  sprintf(c1err,"interp_GRIDMAP returned istat=%d", istat);
	  sprintf(c2err,"Value probably outside GENPDF map range");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	}

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
	  printf("   %s(%s) \n", MAPNAME, GENPDF[IDMAP].GRIDMAP.VARLIST );
	  
	  for(ivar=1; ivar < NDIM; ivar++ ) {
	    IVAR_HOSTLIB = GENPDF[IDMAP].IVAR_HOSTLIB[ivar];
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
	  sprintf(c1err,"N_ITER=%d exceeds bound (prob_ref=%f)", 
		  N_ITER, prob_ref );
	  sprintf(c2err,"Check %s or increase MAX_ITER", MAPNAME );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	}

      } // end while loop over prob

    } // end IDMAP >= 0
  } // end genPDF 

  // - - - - - -
  // check explicit asymGauss function
  if  ( GENGAUSS->USE ) {
    N_EVAL++ ;
    r = exec_GENGAUSS_ASYM(GENGAUSS) ;
  }

  // - - - - - -
  if ( N_EVAL != 1 ) {
    sprintf(c1err,"Found %d PDF options for parName='%s'", N_EVAL, parName);
    sprintf(c2err,"Requires 1 PDF option; check AsymGauss and GENPDF_FILE.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
 
  if ( LDMP ) 
    { printf(" xxx %s: return %s = %f \n", fnam, parName, r); }

  return(r);

} // end get_random_genPDF

// ====================================
int IDMAP_GENPDF(char *parName) {

  // return IDMAP for this input parName.
  // Match against first varName in each map.

  int ID=-9, imap;
  char *tmpName;
  for(imap=0; imap < NMAP_GENPDF; imap++ ) {
    tmpName = GENPDF[imap].VARNAMES[0];
    if ( strcmp(tmpName,parName)==0 ) { return(imap); }
  }

  return(ID);

} // end IDMAP_GENPDF

#endif

