#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

// SNANA stuff
#include "sntools.h"
#include "sntools_devel.h"

void  init_genSmear_devel(void) {

  // Invoke with
  // snlc_sim.exe <myInput>  GENMAG_SMEAR_SCALE DEVEL
  //

  char fnam[] = "init_genSmear_devel" ;

  printf("\n %s: hello \n", fnam);
  fflush(stdout);

  return;

} // end init_genSmear_devel

// ================ 
double get_genSmear_devel(double *parList) {
  // June 1 2021: return scale for G10 intrinsic scatter model:
  // Input parList[0:1] = Trest, x1, c, logmass
  //
  // Invoke with
  //   snlc_sim.exe <myInput>  GENMAG_SMEAR_SCALE DEVEL

  double Trest = parList[0];
  double x1    = parList[1];
  double c     = parList[2];
  double m     = parList[3]; // logMass

  double SCALE = 1.0 ;
  char fnam[] = "get_genSmear_devel" ;

  // ------------ BEGIN -------------

  printf(" xxx %s: Trest=%.2f  x1=%.3f  c=%.3f  m=%.2f \n",
	 fnam, Trest, x1, c, m);  fflush(stdout);

  return SCALE;

} // end get_genSmear_devel


