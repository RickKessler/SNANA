#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

// SNANA stuff
#include "sntools.h"
#include "sntools_devel.h"

void  init_genSmear_devel(int opt) {

  // Invoke with
  // snlc_sim.exe <myInput>  GENMAG_SMEAR_SCALE(DEVEL) <opt>
  //
  char fnam[] = "init_genSmear_devel" ;

  OPT_DEVEL = opt; // store global
  printf("\n %s: OPT_DEVEL = %d \n", fnam, opt);
  fflush(stdout);

  return;

} // end init_genSmear_devel

// =========================================
double get_genSmear_devel(double *parList) {
  // June 1 2021: return scale for G10 intrinsic scatter model:
  // Input parList[0:1] = Trest, x1, c, logmass
  //
  // Invoke with
  //   snlc_sim.exe <myInput>  GENMAG_SMEAR_SCALE(DEVEL) <opt>

  double Trest = parList[0]; // (MJD-PEAKMJD)/(1+z)
  double x1    = parList[1]; // SALT2 x1
  double c     = parList[2]; // SALT2 c
  double m     = parList[3]; // logMass

  double SCALE = 1.0 ;
  char fnam[] = "get_genSmear_devel" ;

  // ------------ BEGIN -------------

  /* 
  printf(" xxx %s: Trest=%.2f  x1=%.3f  c=%.3f  m=%.2f \n",
	 fnam, Trest, x1, c, m);  fflush(stdout);
  */

  return SCALE;

} // end get_genSmear_devel


