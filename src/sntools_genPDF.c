/****************************************************

  Created Jun 12 2020 by R.Kessler

  Generic utility to draw random numbers from probability distribution
  functions (PDF) based on analytical function (asym Gauss), or 
  multi-dimensional map. The map can depend on redshift and any
  combination of HOSTLIB parameters.

  Initial use is for SALT2 parameters (c,x1,RV,AV).

 ****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"
#include "sntools_genPDF.h"

// =====================================================
void init_genPDF(char *fileName) {

  char fnam[] = "init_genPDF";

  // ------------- BEGIN -------------


  return;

} // end init_genPDF

// =====================================================
double get_random_genPDF(char *parName, double zCMB, 
			 GENGAUSS_ASYM_DEF *GENGAUSS){

  double r = 0.0 ;
  char fnam[] = "get_random_genPDF";

  // ------------- BEGIN -----------

  return(r);

} // end get_random_genPDF


