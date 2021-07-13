// =============================
//   sntools_genExpHalfGauss.h
// =============================

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "sntools.h"
//#include "sntools_genExpHalfGauss.h"


// ******************************
void init_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS, double VAL ) {
  // Created March 20 2020 by R.Kessler and D.Brout
  // Init all genGasuss parameter values to VAL .
  
  //  char fnam[] = "init_GEN_EXP_HALFGAUSS" ;
  
  gen_EXP_HALFGAUSS->USE       = false;
  gen_EXP_HALFGAUSS->NAME[0]   = 0;
  gen_EXP_HALFGAUSS->PEAK      = VAL;
  gen_EXP_HALFGAUSS->SIGMA     = VAL;
  gen_EXP_HALFGAUSS->EXP_TAU   = VAL;
  gen_EXP_HALFGAUSS->RANGE[0]  = VAL ;
  gen_EXP_HALFGAUSS->RANGE[1]  = VAL ;
  gen_EXP_HALFGAUSS->RATIO     = 0;
  gen_EXP_HALFGAUSS->INDEX     = -999;

} //  end init_GEN_EXP_HALFGAUSS

void set_GEN_EXPON(double tau, double *range,
		   GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS) {

  // Jul 8 2020
  // Utility to set EXPONENT-only for EXP_HALFGAUSS ... ignore Gaussian.
  gen_EXP_HALFGAUSS->USE      = true ;
  gen_EXP_HALFGAUSS->EXP_TAU  = tau ;
  gen_EXP_HALFGAUSS->RANGE[0] = range[0];
  gen_EXP_HALFGAUSS->RANGE[1] = range[1];

  // ignore Gaussian part
  gen_EXP_HALFGAUSS->PEAK     = 0.0 ;
  gen_EXP_HALFGAUSS->SIGMA    = 0.0 ;
  gen_EXP_HALFGAUSS->RATIO    = 0.0 ;
}


void copy_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *inp_EXP_HALFGAUSS, GEN_EXP_HALFGAUSS_DEF *out_EXP_HALFGAUSS){
  out_EXP_HALFGAUSS->USE = inp_EXP_HALFGAUSS->USE;

  sprintf(out_EXP_HALFGAUSS->NAME,"%s",inp_EXP_HALFGAUSS->NAME);

  out_EXP_HALFGAUSS->PEAK     = inp_EXP_HALFGAUSS->PEAK;
  out_EXP_HALFGAUSS->SIGMA    = inp_EXP_HALFGAUSS->SIGMA;
  out_EXP_HALFGAUSS->EXP_TAU  = inp_EXP_HALFGAUSS->EXP_TAU;
  out_EXP_HALFGAUSS->RANGE[0] = inp_EXP_HALFGAUSS->RANGE[0];
  out_EXP_HALFGAUSS->RANGE[1] = inp_EXP_HALFGAUSS->RANGE[1];
  out_EXP_HALFGAUSS->RATIO    = inp_EXP_HALFGAUSS->RATIO;
  out_EXP_HALFGAUSS->INDEX    = inp_EXP_HALFGAUSS->INDEX;
}


void setUseFlag_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS, char *name ) {
  // Created by D.Brout on March 20 2020
  // This checks if you have set proper profile args 
  // and sets USE flag accordingly.
  // also sets NAME if not null

  double valmax = gen_EXP_HALFGAUSS->RANGE[1] ;
  double sigma  = gen_EXP_HALFGAUSS->SIGMA ;
  double exptau = gen_EXP_HALFGAUSS->EXP_TAU ;

  bool DO_PROFILE = (sigma > 1.0E-9 || exptau > 1.0E-9 );
  bool DO_RANGE   = (valmax > 1.0E-9) ;

  char fnam[] = "setUseFlag_GEN_EXP_HALFGAUSS";

  // ------------ BEGIN --------------

  /* xxxxxxxxxxx mark delete Jul 12 2021 xxxxxxxx
  if (DO_PROFILE && !DO_RANGE) { 
    print_preAbort_banner(fnam);
    printf("SIGMA = %f \n",sigma); 
    printf("EXP_TAU = %f \n",exptau); 
    printf("RANGE[1] = %f \n",valmax);

    sprintf ( c1err, "INVALID PROFILE FOR %s .", name );
    sprintf ( c2err, "CHECK SIM INPUT FILE. RANGE MIGHT NOT BE SET.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  xxxxxx end mark xxxxxx*/


  if (DO_RANGE && DO_PROFILE){
    gen_EXP_HALFGAUSS->USE = true ;
  }
  
  if ( strlen(name) > 0 ){
    sprintf(gen_EXP_HALFGAUSS->NAME,"%s",name);
  } 

} // end setUseFlag_GEN_EXP_HALFGAUSS

// ************************************************************
double funVal_GEN_EXP_HALFGAUSS(double x, GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS){

  //WARNING!! Half Gaussian not yet coded!
  double tau = gen_EXP_HALFGAUSS->EXP_TAU ;
  double funVal = NULLDOUBLE ;
  double sigma = gen_EXP_HALFGAUSS->SIGMA ; 
  char fnam[] = "funVal_GEN_EXP_HALFGAUSS";

  if (sigma > 0.) {
    sprintf(c1err,"Half Gaussian not yet implemented!");
    sprintf(c2err,"TAU = %f, SIGMA = %f", tau, sigma);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
  }

  funVal = exp(-x/tau) ; 


  return funVal;

} // end funVal_GEN_EXP_HALFGAUSS

// ************************************************************
double getRan_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS){
  // March 20 2020
  // generate a random value following profile consisting of 
  // EXP & Half Gaussian
  // TO DO: z dependence of parameters.

  char fnam[] = "getRan_GEN_EXP_HALFGAUSS";
  
  double sig    = gen_EXP_HALFGAUSS->SIGMA ;//half gaussian sigma
  double peak   = gen_EXP_HALFGAUSS->PEAK ; //half gaussian peak
  double tau    = gen_EXP_HALFGAUSS->EXP_TAU;//exponential
  double ratio  = gen_EXP_HALFGAUSS->RATIO;//ratio GAUSS/EXP at zero
  double *range = gen_EXP_HALFGAUSS->RANGE;//generate random within this range
  char *name    = gen_EXP_HALFGAUSS->NAME;

  bool LDMP = false;
  
  bool DOFUN_EXPON = false;
  bool DOFUN_GAUSS = false;

  // always burn randoms to stay synced.                                                                                        
  double ran_EXPON = getRan_Flat1(1) ;                                 
  double ran_GAUSS = getRan_Gauss(1);
  double ran_WGT   = getRan_Flat1(1) ;  
  double epsilon = 1.0E-14;
  double ranval = -9.0 ; //output random value

  // saveforlater if(INPUTS.WV07_GENAV_FLAG) { AV=GENAV_WV07();return(AV); }

  // ----- BEGIN ---------------------

  
  // sanity check
  if ( (range[1] > range[0]) && (tau<epsilon && sig<epsilon) ) {
    sprintf(c1err,"GENRANGE(%s) = %.3f to %.3f, but EXP_TAU=0 and SIG=0", 
	    name, range[0], range[1]);
    sprintf(c2err,"EXP_TAU and/or SIG must not be zero");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if (LDMP) { 
    printf("xxx ----------------------------------------------\n");
    printf("xxx %s DUMP\n",fnam); 
    printf("xxx %s tau=%f sig=%f peak=%f ratio=%f\n",fnam,tau,sig,peak,ratio);
    printf("xxx %s ran_EXPON=%f ran_GAUSS=%f ran_WGT=%f\n",fnam,ran_EXPON,ran_GAUSS,ran_WGT);
  }
  
  if ( range[0] == range[1] )       // delta-function
    { ranval = range[0]; goto DONE ; }
  if ( tau > 50. && ratio == 0.0 )    // flat distribution     
    { ranval = range[0] + (range[1] - range[0]) * ran_EXPON ; goto DONE ; }


  // check for pure exponential or pure Gaussian
  DOFUN_EXPON = ( tau > epsilon &&  sig < epsilon );
  DOFUN_GAUSS = ( sig > epsilon &&  tau < epsilon );

  // check for mixed. 
  if ( tau > epsilon && sig > epsilon && ratio > epsilon ) {
    double WGT_EXPON, WGT_GAUSS, WGT_SUM ;
    WGT_EXPON = 1.0 / tau ;
    WGT_GAUSS = 0.5*ratio / sqrt(TWOPI * sig*sig);
    WGT_SUM   = WGT_EXPON + WGT_GAUSS ;
    if ( ran_WGT < WGT_EXPON/WGT_SUM )
      { DOFUN_EXPON = true; DOFUN_GAUSS = false ;}
    else
      { DOFUN_GAUSS = true; DOFUN_EXPON = false ;}
  }

  if (DOFUN_EXPON && DOFUN_GAUSS) { 
    sprintf(c1err,"Can't have DOFUN_EXPON=T and DOFUN_GAUSS=T\n");     
    sprintf(c2err,"Code logic problem for %s\n",name);                                                              
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);                                                                             
  }

  // pure exponential                                                                                           
  if ( DOFUN_EXPON ) {
    double expmin,expmax,expdif;
    expmin = expf(-range[0]/tau) ;  // note that expmin > expmax !!!                       
    expmax = expf(-range[1]/tau) ;
    expdif = expmin - expmax ;
    ranval = -tau * log( expmin - expdif*ran_EXPON ) ;
    if (LDMP) {
      printf("xxx %s %s=%f (pure exponential)\n",fnam,name,ranval);      
    }
    goto DONE ;
  }


  // pure Guassian (AV > 0 only)                                                                                                                                                                                                              
  if ( DOFUN_GAUSS ) {
    int MAXTRY_ABORT = 10000  ;  // abort after this many tries                                                                                                                                                                               
    int itry = 0 ;
    ranval = -9999.0 ;
    while ( ranval < range[0] || ranval > range[1] ) {
      if ( itry > MAXTRY_ABORT ) {
        sprintf(c1err,"Can't find Gauss-%s between %.2f and %.2f",
                name, range[0], range[1]);
        sprintf(c2err,"after %d tries (sigma=%.2f)\n", itry, sig );
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }
      if ( itry  > 1 ) { ran_GAUSS = getRan_Gauss(1); }
      if (peak > 0.0001 )
        { ranval = sig * ran_GAUSS + peak; }
      else
        { ranval = sig * fabs(ran_GAUSS); }
      itry++ ;
    }
    goto DONE ;
  }


  // if we get here then abort on confusion.                                                                                                                                                                                                  
  sprintf(c1err,"Could not determine %s from Expon. or Gaussian ??", name);
  sprintf(c2err,"tau=%f  sig=%f  ratio=%f", tau, sig, ratio);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  
 DONE:
  return(ranval) ;

} // end of getRan_GEN_EXP_HALFGAUSS()


int parse_input_EXP_HALFGAUSS(char *VARNAME, char **WORDS, int keySource,
                              GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS ) {

  // Created Jul 2021
  // Utility to read gen_EXP_HALFGAUSS structure.

  int N = 0 ;
  char KEYNAME[60];
  char fnam[] = "parse_input_EXP_HALFGAUSS" ;

  // -------------- BEGIN --------------

  // first make sure that VARNAME is contained in WORDS[0] 
  if ( strstr(WORDS[0],VARNAME) == NULL ) { return(N); }

 
  // check Gauss core params
  sprintf(KEYNAME, "GENGAUPEAK_%s",  VARNAME );
  if ( keyMatchSim(1, KEYNAME, WORDS[0], keySource) )  {
    N++; sscanf(WORDS[N], "%le", &gen_EXP_HALFGAUSS->PEAK );
  }

  sprintf(KEYNAME, "GENSIG_%s GENSIGMA_%s",  VARNAME, VARNAME );
  if ( keyMatchSim(1, KEYNAME, WORDS[0], keySource) )  {
    N++; sscanf(WORDS[N], "%le", &gen_EXP_HALFGAUSS->SIGMA );
  }

  // check Exponential param
  sprintf(KEYNAME, "GENTAU_%s",  VARNAME );
  if ( keyMatchSim(1, KEYNAME, WORDS[0], keySource) )  {
    N++; sscanf(WORDS[N], "%le", &gen_EXP_HALFGAUSS->EXP_TAU );
  }


  // Gauss(0)/EXP(0)
  sprintf(KEYNAME, "GENRATIO_%s GENRATIO_%s0",  VARNAME, VARNAME );
  if ( strcmp(VARNAME,"EBV_HOST") ==0 ) 
    { strcat(KEYNAME," EBV0_HOST"); } // back-compatible
  if ( keyMatchSim(1, KEYNAME, WORDS[0], keySource) )  {
    N++; sscanf(WORDS[N], "%le", &gen_EXP_HALFGAUSS->RATIO );
  }

  sprintf(KEYNAME, "GENRANGE_%s",  VARNAME );
  if ( keyMatchSim(1, KEYNAME, WORDS[0], keySource) )  {
    N++; sscanf(WORDS[N], "%le", &gen_EXP_HALFGAUSS->RANGE[0] );
    N++; sscanf(WORDS[N], "%le", &gen_EXP_HALFGAUSS->RANGE[1] );
  }

  // set NAME and USE flag
  if ( N > 0 ) 
    {  setUseFlag_GEN_EXP_HALFGAUSS(gen_EXP_HALFGAUSS,VARNAME); }

  return N ;

} // end parse_input_EXP_HALFGAUSS
