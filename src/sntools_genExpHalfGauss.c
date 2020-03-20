// =============================
//   sntools_genExpHalfGauss.h
// =============================

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"
#include "sntools_genExpHalfGauss.h"


// ******************************
void init_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS, double VAL ) {
  // Created March 20 2020 by R.Kessler and D.Brout
  // Init all genGasuss parameter values to VAL .
  
  char fnam[] = "init_GEN_EXP_HALFGAUSS" ;
  
  printf("xxx Hello from %s\n",fnam);

  gen_EXP_HALFGAUSS->USE = false;

  gen_EXP_HALFGAUSS->PEAK = VAL;
  gen_EXP_HALFGAUSS->SIGMA = VAL;
  gen_EXP_HALFGAUSS->EXP_TAU = VAL;
  gen_EXP_HALFGAUSS->RANGE[0]  = VAL ;
  gen_EXP_HALFGAUSS->RANGE[1]  = VAL ;
  gen_EXP_HALFGAUSS->RATIO = 0;
  

}
