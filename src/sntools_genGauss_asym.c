// =============================
//   sntools_genGauss_asym.h
//
//  Jun 12 2020: remove obsolete/unused skewnormal code
// =============================

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "sntools.h"
#include "sntools_genGauss_asym.h"


// ******************************
void init_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss, double VAL ) {

  // Created Apr 13 2016 by R.Kessler
  // Init all genGasuss parameter values to VAL .

  sprintf(genGauss->NAME,"NULL");

  genGauss->USE       = false ;
  genGauss->PEAK      = VAL ;

  genGauss->RANGE[0]  = VAL ;
  genGauss->RANGE[1]  = VAL ;

  genGauss->SIGMA[0]  = VAL ;
  genGauss->SIGMA[1]  = VAL ;

  // SKEW params = 0 regardless of VAL
  genGauss->SKEW[0]       = 0.0 ; 
  genGauss->SKEW[1]       = 0.0 ;  

  genGauss->NGRID     = 0 ;
  genGauss->FUNINDEX  = 0;

  // Mar 29 2017: 2nd peak
  genGauss->PROB2     = 0.0 ;
  genGauss->PEAK2     = 0.0 ;
  genGauss->SIGMA2[0] = 0.0 ;
  genGauss->SIGMA2[1] = 0.0 ;

  genGauss->RMS       = 0.0 ;

} // end init_GENGAUSS_ASYM

// ************************************
void copy_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss1,
			GENGAUSS_ASYM_DEF *genGauss2) {

  // Aug 30 2016
  // copy contents of genGauss1 into genGauss2
  int i ;

  genGauss2->USE = genGauss1->USE ;

  sprintf(genGauss2->NAME,"%s", genGauss1->NAME );

  genGauss2->PEAK      =  genGauss1->PEAK ;
  genGauss2->NGRID     =  genGauss1->NGRID ;

  for(i=0; i < 2; i++ ) {
    genGauss2->RANGE[i]  = genGauss1->RANGE[i] ;
    genGauss2->SIGMA[i]  = genGauss1->SIGMA[i] ;
    genGauss2->SKEW[i]   = genGauss1->SKEW[i] ;
  }

  genGauss2->FUNINDEX = genGauss1->FUNINDEX ;

  // 2nd peak
  genGauss2->PROB2      =  genGauss1->PROB2 ;
  genGauss2->PEAK2      =  genGauss1->PEAK2 ;
  genGauss2->SIGMA2[0]  =  genGauss1->SIGMA2[0] ;
  genGauss2->SIGMA2[1]  =  genGauss1->SIGMA2[1] ;

  genGauss2->RMS        =  genGauss1->RMS ;

  return ;

} // end copy_GENGAUSS_ASYM

// ======================================
void prepIndex_GENGAUSS(char *varName, GENGAUSS_ASYM_DEF *genGauss ) {

  // Created Sep 2 2016
  // Store NAME and increment index.
  // Called by readFile routine and command-line read function.
  //
  // Jun 11 2020: moved from snlc_sim.h, and set USE=true.

  char *ptrName = genGauss->NAME;
  //  char fnam[] = "prepIndex_GENGAUSS" ;
  // ---------- BEGIN ---------

  // if genGauss name is not set, then set name and FUNINDEX
  if ( strcmp(ptrName,varName) != 0 ) {
    genGauss->USE      = true ;
    genGauss->FUNINDEX = NFUN_GENGAUSS_ASYM ;
    NFUN_GENGAUSS_ASYM++;  
    sprintf(genGauss->NAME, "%s", varName);  
  }

  // copy each GENGAUSS_ASYM struct into master list in case
  // some kind of global operation or init is needed.

  int FUNINDEX = genGauss->FUNINDEX ;
  copy_GENGAUSS_ASYM( genGauss, &GENGAUSS_ASYM_LIST[FUNINDEX] );

} // end prepIndex_GENGAUSS


// **********************************
double exec_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss) {

  // Created May 2012: 
  // return random Guassian number with above params
  // If siglo & sighi > 100*(hi-lo), then just assume a flat
  // distribution between lo & hi.
  // Same as BIGAUSRAN, but here a struct is passed with all the args.
  //
  // Mar 16 2014: use new SKEW parameter
  // Apr 14 2016: sigmax -> 10*range [was 100*range]
  // Apr 20 2016: check NGRID option to snap to grid
  // Aug 30 2016: add call to skewNormal().
  // Oct 02 2016: compute rangeDif after DO_GRID if-block (bug fix)
  // Mar 29 2017: check for 2nd peak
  // Feb 28 2018: check for correlated randoms (see 'redCor')
  // Jun 12 2020: 
  //   + remove skewnormal (never worked)
  //   + return -9 if !USE
  // 
  
  double peak, lo, hi, siglo, sighi, skewlo, skewhi, xlo, xhi ; 
  double gridsize, grid0;
  int NTRY, DO_SKEWSIGMA, DO_GRID;
  int NGRID, FUNINDEX, j ;
  int USE_PEAK1=1, MXTRY = 1000, LDMP=0 ;
  double ranval=-9.0, rangeDif, RANGE[2], sigmax, ran1, ran2, PROB2 ;
  char *NAME  = genGauss->NAME;
  char fnam[] = "exec_GENGAUSS_ASYM" ;

  // ---------- BEGIN -------------

  //  LDMP = (strstr(NAME,"SALT2") != NULL ); // xxx REMOVE

  // always burn random to stay synced.
  ran1 = FlatRan1(1) ;

  if ( !genGauss->USE ) {  return(ranval);  }

  /*
  if ( !genGauss->USE ) {  
    printf(" xxx %s: undefined '%s' PEAK=%f  RANGE = %f to %f\n", 
    	   fnam, NAME, genGauss->PEAK,
	   genGauss->RANGE[0], genGauss->RANGE[1] ); 
    return(ranval); 
  }
  */

  // check optional 2nd peak (Mar 2017)
  PROB2 = genGauss->PROB2 ;
  if ( PROB2 > 0.0000001 ) {
    ran2 = FlatRan1(1) ;
    if ( ran2 < PROB2 ) {
      peak        = genGauss->PEAK2 ;
      siglo       = genGauss->SIGMA2[0] ;
      sighi       = genGauss->SIGMA2[1] ;
      skewlo = skewhi = 0.0 ;
      USE_PEAK1   = 0;
    }
  }

  if ( USE_PEAK1 ) {
    peak        = genGauss->PEAK ;
    siglo       = genGauss->SIGMA[0] ;
    sighi       = genGauss->SIGMA[1] ;
    skewlo      = genGauss->SKEW[0] ;
    skewhi      = genGauss->SKEW[1] ;
  }

  lo          = genGauss->RANGE[0] ;
  hi          = genGauss->RANGE[1] ;
  NGRID       = genGauss->NGRID ;
  FUNINDEX    = genGauss->FUNINDEX ;

  RANGE[0]=lo; RANGE[1]=hi;

  // abort on crazy NGRID value
  if ( NGRID < 0 || NGRID > 100 ) {
    print_preAbort_banner(fnam);
    printf("\t peak = %f \n", peak);
    printf("\t range(lo,hi) = %f, %f \n", lo, hi );
    printf("\t sigma(lo,hi) = %f, %f \n", siglo, sighi );
    sprintf(c1err,"Crazy NGRID=%d", NGRID );
    errmsg(SEV_FATAL, 0, fnam, c1err, "" );
  }


  gridsize = grid0 = -9.0 ;
  DO_GRID = ( NGRID >=2 && hi>lo && siglo>0.0 && sighi>0.0 );
  if ( DO_GRID ) {
    // compute gridsize here to allow command-line override of RANGE
    gridsize = (hi-lo)/(double)(NGRID-1) ;
    grid0    = lo ; // store first grid location
    lo   -= gridsize/2.0 ;
    hi   += gridsize/2.0 ;
  } 

  rangeDif = hi-lo;
  NTRY  = 0;
  sigmax = 10.*rangeDif ;
  DO_SKEWSIGMA  = ( fabs(skewlo)  > 1.0E-9 || fabs(skewhi) > 1.0E-9 ) ;

  
  if ( LDMP ) {
    printf("\t xxx ----------------- %s ------------ \n", NAME);
    printf("\t xxx peak = %f \n", peak);
    printf("\t xxx range(lo,hi) = %f, %f \n", lo, hi );
    printf("\t xxx sigma(lo,hi) = %f, %f \n", siglo, sighi );
    printf("\t xxx DO_GRID=%d  sigmax=%f \n", DO_GRID, sigmax);
    printf("\t xxx rangeDif=%f   ran1=%f \n", rangeDif, ran1);
    fflush(stdout);
  }

  // always burn random to stay synced.

  if ( lo == hi ) {
    ranval = lo ;    // delta function
  }
  else if ( siglo > sigmax && sighi > sigmax ) {
    ranval = lo + rangeDif*ran1 ;  // flat distribution
  }
  else {
  GENVAL:
    NTRY++ ;
    if ( NTRY > MXTRY ) {
      print_preAbort_banner(fnam);
      dump_GENGAUSS_ASYM(genGauss);
      printf(" DO_SKEW[SIGMA] = %d, %d \n",  DO_SKEWSIGMA);

      sprintf(c1err,"Could not find %s RANDOM after %d tries ", 
	      NAME, NTRY );
      sprintf(c2err,"Something is crazy.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
    
    if ( DO_SKEWSIGMA  )  {
      xlo = lo - peak ;
      xhi = hi - peak ;
      ranval = peak + skewGaussRan(xlo,xhi,siglo,sighi,skewlo,skewhi); 
    }
    else { 
      ranval = peak + biGaussRan(siglo,sighi) ; 
    }
    
    if ( ranval < lo ) { goto GENVAL;  } 
    if ( ranval > hi ) { goto GENVAL;  } 
  }


  // April 2016: check option to snap to grid
  if ( DO_GRID ) {
    double ranval_orig = ranval ;
    int    ibin = (int)((ranval_orig-lo)/gridsize) ;
    double xbin = (double)ibin ;
    ranval      = grid0 + (xbin*gridsize) ;

    /*
    printf(" xxx ranval=%.3f->0.3f  ibin=%d  grid0=%.3f\n",
    ranval_orig, ranval, ibin, grid0); fflush(stdout); */
  }

  if ( LDMP ) {
    printf("\t xxx ranval = %f \n", ranval); 
    fflush(stdout);
  }

  return(ranval) ;

} // end of exec_GENGAUSS_ASYM


// ****************************************
void dump_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss) {

  // Created Sep 2 2016 
  // screen-Dump contents of genGauss

  double *ptrVal ;
  char fnam[] = "dump_GENGAUSS_ASYM" ;

  // ---------- BEGIN -----------

  printf("\n");
  printf("# ------------------------------------------- \n");
  printf(" START %s  for '%s' (INDEX=%d) \n",
	 fnam, genGauss->NAME, genGauss->FUNINDEX);

  printf("\t PEAK = %.3f  \n", genGauss->PEAK);

  ptrVal = genGauss->RANGE ;
  printf("\t RANGE = %.3f to %.3f \n", 	ptrVal[0], ptrVal[1] );

  ptrVal = genGauss->SIGMA; 
  printf("\t SIGMA(-/+) = %.3f / %.3f \n", ptrVal[0], ptrVal[1] );
	  
  printf(" END %s \n", fnam );
  printf("# ------------------------------------------- \n");

  fflush(stdout);

  return ;

} // end dump_GENGAUSS_ASYM


