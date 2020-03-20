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

  gen_EXP_HALFGAUSS->NAME[0] = 0;

  gen_EXP_HALFGAUSS->PEAK = VAL;
  gen_EXP_HALFGAUSS->SIGMA = VAL;
  gen_EXP_HALFGAUSS->EXP_TAU = VAL;
  gen_EXP_HALFGAUSS->RANGE[0]  = VAL ;
  gen_EXP_HALFGAUSS->RANGE[1]  = VAL ;
  gen_EXP_HALFGAUSS->RATIO = 0;
  
  

}

void setUseFlag_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS, char *name ) {
  // Created by D.Brout on March 20 2020
  // This checks if you have set proper profile args 
  // and sets USE flag accordingly.
  // also sets NAME if not null

  char fnam[] = "setUseFlag_GEN_EXP_HALFGAUSS";

  double valmax = gen_EXP_HALFGAUSS->RANGE[1] ;
  double sigma = gen_EXP_HALFGAUSS->SIGMA ;
  double exptau = gen_EXP_HALFGAUSS->EXP_TAU ;

  bool DO_PROFILE = (sigma > 1.0E-9 || exptau > 1.0E-9 );
  bool DO_RANGE = (valmax > 1.0E-9) ;

  if (DO_PROFILE && !DO_RANGE) { 
    print_preAbort_banner(fnam);
    printf("SIGMA = %f",sigma); printf("EXP_TAU = %f",exptau); printf("RANGE[1] = %f",valmax);

    sprintf ( c1err, "INVALID PROFILE FOR %s .", name );
    sprintf ( c2err, "CHECK SIM INPUT FILE. RANGE MIGHT NOT BE SET.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if (DO_RANGE && DO_PROFILE){
    gen_EXP_HALFGAUSS->USE = true ;
  }

  


  if (strlen(name)>0){
    sprintf(gen_EXP_HALFGAUSS->NAME,"%s",name);
  } 

}
