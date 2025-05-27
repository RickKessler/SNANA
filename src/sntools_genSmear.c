/********************************
 Created April 2012 by R.Kessler

 models of intrinsic brightness variation as a function of wavelength.
 The init functions depend on the model, but there is just one 
 function to evaluate the smear. Flags to set in the sim-input 

 GENMAG_SMEAR_MODELNAME: G10      # SALT2 model disp and SIGCOH
 GENMAG_SMEAR_MODELNAME: Guy10    # idem
 GENMAG_SMEAR_MODELNAME: G10Fig8  # Fig 8 only (SIGCOH=0, FUDGE=1)

 GENMAG_SMEAR_MODELNAME: C11      # or Chotard11 (farUV uncorrelated)
 GENMAG_SMEAR_MODELNAME: C11_0    # or Chotard11 (farUV uncorrelated)
 GENMAG_SMEAR_MODELNAME: C11_1    # farUV 100% correlated with U
 GENMAG_SMEAR_MODELNAME: C11_2    # farUV 100% anti-correlated with U

 GENMAG_SMEAR_MODELNAME: CCM89    # modify color law with CCM89  
 GENMAG_SMEAR_MODELNAME: PRIVATE  
 GENMAG_SMEAR_MODELNAME: COH           # SIGCOH=.13  (Feb 2014)
 GENMAG_SMEAR_MODELNAME: COH(0.12)     # SIGCOH=.12  (Nov 2019)
 GENMAG_SMEAR_MODELNAME: COH(0.9+0.2)  # 2 sigma-scat terms (Nov 2019)
 GENMAG_SMEAR_MODELNAME: VCR      # Velocit-Color Relation (MFK14)
 GENMAG_SMEAR_MODELNAME: BIMODAL_UV  # sort'of like Milne 2015
 GENMAG_SMEAR_USRFUN:  <list of 8 params> 

 GENMAG_SMEAR_MODELNAME: OIR    # Optical+IR using CSP+CfA (under development)
 GENMAG_SMEAR_MODELNAME: [path]/COVSED.[version] # cov(wave x phase), Oct 2019
 
 External program must set
   NSMEARPAR_OVERRIDE = 0 

 then call
   init_genSmear_FLAGS() ;

 Then call one of
   init_genSmear_Guy10
   init_genSmear_Guy10Fig8  | Fig 8 only (no lam-fudge and SIGCOH=0)
   init_genSmear_Chotard11
   init_genSmear_private
   init_genSmear_COH
   init_genSmear_OIR
   init_genSmear_USRFUN
   etc ...


 To see if anything is set
   istat = istat_genSmear() = 0 (nothing set)
                            = 1 (a model is set)

 Then to get magSmear array for input *Lam array,
   get_genSmear(...)


  HISTORY
 ~~~~~~~~~~~

 Jan 24 2013: Fixed init_genSmear_G10 bug preventing local 
              smear file from loading (JLM)

 May 01, 2014: add SALT2.JLA smearing model; very similar to G10
               but slightly larger sigma_int.

 July 30 2015: add BIMODAL_UV model

 April 6 2016: fix awful bug in get_genSmear_COH(), and .12 -> 0.13

 Jun 14 2016: set global MAGSMEAR_COH so that it ends up in output table
              as SIM_MAGSMEAR_COH

 July 29 2016: new utility exec_genSmear_override().

 Aug 30 2019: RK,DJ - start adding infrastructure for OIR model.
 Sep 03 2019: 
   + include genmag_SALT2[SEDtools].h
   + after reading salt2 dispersion file, clip wavelengths outside
     SALT2 wave range. Avoids abort at very low redshifts.

 Oct 18 2019: begin COVSED model [based on SNF]

 Oct 21 2019: 
   + significant refactor for speed improvements using repeat_genSmear().
     BEWARE to check repeat_genSmear for Trest-dependent smear models.

 Nov 30 2019: 
   + refactor and upgrade COH model to pass 1 or 2 sigma values.
      Default is still sigma=0.13 mag.
 
 Feb 17 2020
   + refactor correlated Gauss randoms to use init_Cholesky and
     GaussRanCorr utilities.

**********************************/

#include "fitsio.h"
#include "sntools.h"
#include "genmag_SEDtools.h"
// xxx #include "genmag_extrap.h"
#include "genmag_SALT2.h"
#include "sntools_genSmear.h"
#include "sntools_devel.h"
#include "MWgaldust.h"
#include "sntools_dataformat_fits.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>


// ===========================================================
//
//    Begin functions
//
// ===========================================================

int  istat_genSmear(void) {
  return(GENSMEAR.NUSE);
}

void  init_genSmear_FLAGS(int MSKOPT, char *SCALE_STRING) {

  // Mar 22 2020: 
  //  replace double SCALE with string that may have polynom funct

  char fnam[] = "init_genSmear_FLAGS" ;

  GENSMEAR.NUSE           = 0 ; // number of GENSMEAR init calls
  GENSMEAR.NCALL          = 0 ; // number of get_genSmear calls.

  GENSMEAR_USRFUN.USE     = 0 ;
  GENSMEAR_SALT2.USE      = 0 ;
  GENSMEAR_C11.USE        = 0 ;
  GENSMEAR_CCM89.USE      = 0 ;
  GENSMEAR_VCR.USE        = 0 ;
  GENSMEAR_BIMODAL_UV.USE = 0 ;
  GENSMEAR_OIR.USE        = 0 ;
  GENSMEAR_COVSED.USE     = 0 ;
  GENSMEAR_PHASECOR.USE   = 0 ;
  GENSMEAR_SCALE.USE      = 0 ;

  GENSMEAR.NSET_RANGauss  = 0 ;
  GENSMEAR.NSET_RANFlat   = 0 ;
  GENSMEAR.MSKOPT         = MSKOPT ; // Oct 2019

  // hard-wire wavelengths to monitor COVARIANCE between 
  // 2 arbitrary wavelengths
  GENSMEAR.SUMSMEAR_CHECK[0]   =  0.0 ;
  GENSMEAR.SUMSMEAR_CHECK[1]   =  0.0 ;
  GENSMEAR.SQSUMSMEAR_CHECK[0] =  0.0 ;
  GENSMEAR.SQSUMSMEAR_CHECK[1] =  0.0 ;
  GENSMEAR.SUMCROSS            =  0.0 ;
  GENSMEAR.NCHECK              =  0 ;

  GENSMEAR.CID_LAST    = -9;
  GENSMEAR.TREST_LAST  = 99999. ;
  GENSMEAR.NLAM_LAST   = -9 ;
  GENSMEAR.LAMMIN_LAST = -9.0 ;
  GENSMEAR.LAMMAX_LAST = -9.0 ;
    
  // abort on unitialized NSMEARPAR_OVERRIDE
  if ( NSMEARPAR_OVERRIDE < 0 || NSMEARPAR_OVERRIDE >= MXSMEARPAR_OVERRIDE ) {
    sprintf(c1err,"Unitialized NSMEARPAR_OVERRIDE=%d", NSMEARPAR_OVERRIDE);
    sprintf(c2err,"External calling routine must set NSMEARPAR_OVERRIDE=0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // Jan 2020: allocate magsmear-vs-wave array
  int MEMD = MXLAM_GENSMEAR_SALT2*sizeof(double);
  GENSMEAR.MAGSMEAR_LIST = (double*) malloc(MEMD);


  init_genSmear_SCALE(SCALE_STRING);

  return ;

}  // end of init_genSmear_FLAGS


// *************************************
void init_genSmear_SCALE(char *SCALE_STRING) {
  
  // Created Mar 22 2020 by R.Kessler
  // SCALE_STRING = '[VARNAME] [POLYSTRING]'

  char VARNAME[40], cPOLY[40], *cptr[2];
  char space[] = " " ;
  int  NTMP ;
  char fnam[] = "init_genSmear_SCALE" ;
    
  // ---------- BEGIN ----------

  GENSMEAR_SCALE.GLOBAL = 1.0 ; // default
  if ( strlen(SCALE_STRING) == 0 ) { return ; }

  GENSMEAR_SCALE.USE = 1;

  // split string to get variable name and polynom string
  cptr[0] = VARNAME;
  cptr[1] = cPOLY;

  splitString(SCALE_STRING, space, fnam, 2, &NTMP, cptr);

  // xxx  printf(" xxx VARNAME='%s'  cPOLY = '%s' \n", VARNAME, cPOLY);

  parse_GENPOLY(cPOLY, VARNAME, &GENSMEAR_SCALE.POLY, fnam);
    sprintf(GENSMEAR_SCALE.VARNAME, "%s", VARNAME);

  if ( GENSMEAR_SCALE.POLY.ORDER == 0 ) {
    GENSMEAR_SCALE.GLOBAL = GENSMEAR_SCALE.POLY.COEFF_RANGE[0][0];
    printf("\t Global GENMAG_SMEAR_SCALE: %.3f \n", GENSMEAR_SCALE.GLOBAL);
  }
  else {
    printf("\t GENMAG_SMEAR_SCALE: poly(%s) = %s \n", VARNAME, cPOLY);
  }

  // xxxxxxxxx DEVEL OPTION, May 31 2021
  if ( strcmp(VARNAME,"DEVEL") == 0 ) { 
    int opt_devel = (int)GENSMEAR_SCALE.GLOBAL ;
    init_genSmear_devel(opt_devel); 
    return;  
  }
  // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  return;

} // end of init_genSmear_SCALE


// *******************************************************
void init_genSmear_randoms(int NRANGauss, int NRANFlat) {

  // Created Oct 21 2019
  // Store number of Gauss and Flat randoms, and
  // malloc memory to hold Gauss and flat randoms for each event.

  int MEMD;
  char fnam[] = "init_genSmear_randoms";

  // ------------ BEGIN --------------

  GENSMEAR.NGEN_RANGauss = NRANGauss ;
  GENSMEAR.NGEN_RANFlat  = NRANFlat  ;

  MEMD = (NRANGauss + MXFILTINDX + 1) * sizeof(double);
  GENSMEAR.RANGauss_LIST = (double*) malloc(MEMD);

  MEMD = (NRANFlat + MXFILTINDX + 1) * sizeof(double);
  GENSMEAR.RANFlat_LIST = (double*) malloc(MEMD);

  printf("   Number of Gaussian  randoms per SN: %d \n", 
         NRANGauss);
  printf("   Number of Flat(0-1) randoms per SN: %d \n", 
         NRANFlat);  

  fflush(stdout);

  return;
} // end init_genSmear_randoms


void load_genSmear_randoms(int CID, double gmin, double gmax, double RANFIX) {

  // Created Oct 21 2019
  // load Gauss and Flat randoms into global GENSMEAR struct.
  // Note that at least MXFILTINDX randoms are stored to 
  // preserve random synch ... may change this later.
  //
  // Inputs:
  //    CID to track new event.
  //    gmin = min Gauran to clip Gaussian
  //    gmax = max Gauran to clip Gaussian
  //    RANFIX is a user option to fix random Gaussians for debug

  int  NRANGauss = GENSMEAR.NGEN_RANGauss ;
  int  NRANFlat  = GENSMEAR.NGEN_RANFlat ;
  int  ILIST_RAN = 2 ; // list to use for genSmear randoms  
  int  LDMP = 0 ;

  int iran;
  char fnam[] = "load_genSmear_randoms" ;

  // ------------- BEGIN ------------

  GENSMEAR.CID = CID;

  // generate Guassian randoms for intrinsic scatter [genSmear] model
  if ( NRANGauss < MXFILTINDX-1 ) { NRANGauss = MXFILTINDX-1; } 
  for(iran=0; iran < NRANGauss; iran++ ) {                 
    GENSMEAR.RANGauss_LIST[iran] = getRan_GaussClip(ILIST_RAN,gmin,gmax); 
  }

 
  // repeat for 0-1 [flat] randoms
  if ( NRANFlat < MXFILTINDX-1 ) { NRANFlat = MXFILTINDX-1; }
  for ( iran=0; iran < NRANFlat; iran++ ) 
    {  GENSMEAR.RANFlat_LIST[iran]  = getRan_Flat1(ILIST_RAN);  } 
  
  if ( LDMP ) {
    double GFIRST = GENSMEAR.RANGauss_LIST[0];
    double GLAST  = GENSMEAR.RANGauss_LIST[NRANGauss-1];
    double FFIRST = GENSMEAR.RANFlat_LIST[0];
    double FLAST  = GENSMEAR.RANFlat_LIST[NRANFlat-1];
    printf(" xxx load_ran: CID=%d  First/Last = %f/%f(G), %f/%f(F) \n",
	   CID, GFIRST,GLAST, FFIRST,FLAST);  fflush(stdout);
  }


  if ( RANFIX > -99.0 ) { 
    for(iran=0; iran < NRANGauss; iran++ ) {                 
      GENSMEAR.RANGauss_LIST[iran] = RANFIX ;
    }
  }

  // Feb 12 2020: check randoms for phaseCor model
  if ( GENSMEAR_PHASECOR.USE  ) {  
    int NBIN = GENSMEAR_PHASECOR.NBIN ;
    for(iran=0; iran < NBIN; iran++ ) {                 
      GENSMEAR_PHASECOR.RANGauss_LIST[iran] = 
	getRan_GaussClip(ILIST_RAN,gmin,gmax); 
    }
  }

  return;

} // end load_genSmear_randoms


// ********************************
void get_genSmear(double *parList, int NLam, double *Lam, 
		  double *magSmear) {

  // April 2012
  // Figure out which smearing model was initialized,
  // and return magSmear at each *Lam.
  //
  // Inputs:
  //     parList[0:3] : Trest, x1, c, logMass
  //     NLam   : number of wave bins
  //    *Lam    : array of rest-frame wavelenths
  //
  // Output
  //   magSmear : magSmear at each *Lam bin.
  //
  //
  // Oct 9 2018: check option to scale the magSmear values
  // Oct 21 2019: add CID argument
  // Nov 30 2019: MAGSMEAR_COH -> MAGSMEAR_COH[2]
  // Feb 17 2020: add c & x1 input args
  // May 31 2021: refactor to pass parList that includes logMass

  double Trest   = parList[0];
  double x1      = parList[1];
  double c       = parList[2];
  double logMass = parList[3];

  int ilam, repeat ;
  char fnam[] = "get_genSmear" ;

  // -------------- BEGIN -----------

  GENSMEAR.NCALL++ ;

  repeat = repeat_genSmear(Trest,NLam,Lam);
  if ( repeat ) {  goto SET_LAST; }

  // abort if more than one model has been initialized.
  if ( GENSMEAR.NUSE > 1 ) {
    sprintf(c1err,"%d GENSMEAR models initilialized", GENSMEAR.NUSE);
    sprintf(c2err,"but only one model allowed.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  GENSMEAR.MAGSMEAR_COH[0] = 0.0 ;
  GENSMEAR.MAGSMEAR_COH[1] = 0.0 ;

  for(ilam=0; ilam < NLam; ilam++ ) { magSmear[ilam] = 0.0 ; }

  if ( GENSMEAR_USRFUN.USE ) {
    get_genSmear_USRFUN(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_SALT2.USE ) {
    get_genSmear_SALT2(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_C11.USE ) {
    get_genSmear_Chotard11(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_VCR.USE ) {
    get_genSmear_VCR(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_CCM89.USE ) {
    get_genSmear_CCM89(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_COH.USE ) {
    get_genSmear_COH(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_PRIVATE.USE ) {
    get_genSmear_private(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_BIMODAL_UV.USE ) {
    get_genSmear_biModalUV(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_OIR.USE ) {
    get_genSmear_OIR(Trest, NLam, Lam, magSmear) ;
  }
  else if ( GENSMEAR_COVSED.USE ) {
    get_genSmear_COVSED(Trest, NLam, Lam, magSmear) ;
  }
  else {
    sprintf(c1err,"Unknown smear model.");
    sprintf(c2err,"Check GENMAG_SMEAR_MODELNAME key.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // Mar 2020: check option to scale the smearing vs. c & x1
  if ( GENSMEAR_SCALE.USE ) {    
    double SCALE = get_genSmear_SCALE(parList);
    for(ilam=0; ilam < NLam; ilam++ ) { magSmear[ilam] *= SCALE ; }
  }


 SET_LAST:
  GENSMEAR.CID_LAST    = GENSMEAR.CID ;
  GENSMEAR.TREST_LAST  = Trest;
  GENSMEAR.NLAM_LAST   = NLam ;
  GENSMEAR.LAMMIN_LAST = Lam[0]; 
  GENSMEAR.LAMMAX_LAST = Lam[NLam-1];

  return ;

} // end of get_genSmear


// ********************************************************
int repeat_genSmear(double Trest, int NLam, double *Lam) {
  
  // Created Oct 21 2019
  // Return 0 --> get_genSmear returns new set of magSmear values
  // Return 1 --> get_genSmear returns same magSmear values
  //
  // If sim-input key
  //   GENMAG_SMEAR_MSKOPT:  32
  // then always return NEW (no repeats).
  //

  int REPEAT  = 1;
  int NEW     = 0;
  int LDMP    = 0 ;
  double LAMMIN = Lam[0];
  double LAMMAX = Lam[NLam-1];
  char fnam[] = "repeat_genSmear";

  // ------------ BEGIN ------------

  // for 1 lam bin, it's central filter wavelength,
  // so always re-calculate.
  if ( NLam == 1 ) { return(NEW); }

  if ( (GENSMEAR.MSKOPT & 32 )>0 ) { return(NEW); }

  if ( LDMP ) {
    printf(" xxx ---------------------------------- \n");
    printf(" xxx CID(%d,%d)  Trest(%.2f,%.2f) \n",
	   GENSMEAR.CID, GENSMEAR.CID_LAST,  Trest, GENSMEAR.TREST_LAST ) ;
    printf(" xxx LAMMIN(%.1f,%.1f)  LAMMAX(%.1f,%.1f) \n"
	   ,LAMMIN, GENSMEAR.LAMMIN_LAST
	   ,LAMMAX, GENSMEAR.LAMMAX_LAST );
    fflush(stdout);
  }
  
  if ( GENSMEAR.CID != GENSMEAR.CID_LAST    ) { return(NEW); }
  if ( LAMMIN       != GENSMEAR.LAMMIN_LAST ) { return(NEW); }
  if ( LAMMAX       != GENSMEAR.LAMMAX_LAST ) { return(NEW); }

  if ( LDMP ) { printf(" xxx REPEAT \n"); fflush(stdout); }

  return(REPEAT);

} // end repeat_genSmear


// ***********************************
void init_genSmear_COVLAM_debug(double *lam, double COVMAT[2][2]) {

  // Inputs
  //   lam[0], lam[1] = two wavelengths to monitor COV
  //   COVMAT = expected 2x2 COV matrix 
  //
  // Simply print expected RMS and reduced covar,
  // to be compared during generation with empirical 
  // values.
  // 
  double SIG0 = sqrt( COVMAT[0][0] );
  double SIG1 = sqrt( COVMAT[1][1] );
  double RHO  = COVMAT[0][1] / (SIG0*SIG1);
  char fnam[] = "init_genSmear_COVLAM_debug" ;

  // ----------- BEGIN ---------

  printf("\n");
  printf(" %s: Expected RMS(lam=%.0f) = %.5f \n", fnam, lam[0], SIG0 );
  printf(" %s: Expected RMS(lam=%.0f) = %.5f \n", fnam, lam[1], SIG1 );
  printf(" %s: Expected Reduced Covar = %.5f \n", fnam, RHO );
  printf("\n");
  fflush(stdout);

  return ;

} // end init_genSmear_COVLAM_debug


void update_genSmear_COVLAM_debug(double *magSmear) {

  // Inputs:
  //  magSmear[0], magSmear[1] = magSmear values for each wavelength.

  double mtmp, XN, SUM, SQSUM, MEAN[2], RMS[2], COV_01, RHO ;
  int i, N;
  char fnam[] = "update_genSmear_COVLAM_debug" ;

  // ----------- BEGIN -------

  for(i=0; i < 2; i++ ) {
    mtmp = magSmear[i];
    GENSMEAR.SUMSMEAR_CHECK[i]   +=  mtmp;
    GENSMEAR.SQSUMSMEAR_CHECK[i] += (mtmp*mtmp);    
  }
  GENSMEAR.SUMCROSS += (magSmear[0]*magSmear[1]);
  GENSMEAR.NCHECK++ ;

  if ( (GENSMEAR.NCHECK % 100) == 0 ) {
    N     = GENSMEAR.NCHECK;
    XN    = (double)N;

    for(i=0; i < 2; i++ ) {
      SUM      = GENSMEAR.SUMSMEAR_CHECK[i];
      SQSUM    = GENSMEAR.SQSUMSMEAR_CHECK[i];
      MEAN[i]  = SUM/XN;
      RMS[i]   = STD_from_SUMS(N, SUM, SQSUM);
    }

    COV_01 = 
      GENSMEAR.SUMCROSS + 
      GENSMEAR.SUMSMEAR_CHECK[0]*MEAN[1] +
      GENSMEAR.SUMSMEAR_CHECK[1]*MEAN[0] +
      MEAN[0]*MEAN[1] 
      ;
    COV_01 /= XN;
    RHO = COV_01/(RMS[0]*RMS[1]);

    printf(" COV_debug: N=%4d  <smear>=%6.3f/%6.3f  "
	   "RMS=%.4f/%.4f  RHO=%5.3f\n",
	   N, MEAN[0], MEAN[1], RMS[0], RMS[1], RHO );
    fflush(stdout);
  }

  return ;

} // end  update_genSmear_COVLAM_debug


// *********************************************
void init_genSmear_USRFUN(int NPAR, double *parList, double *LAMRANGE) {

  /*
    Created Mar 6, 2012 by R.Kessler
    init function for intrinsic mag-Smearing.
    *parList =
    0. Coherent
    1. Amp_5500
    2. lambda interval (L1)
    3. lambda phase (L0: L0=0 => random phase)
    4. tau_lam 
    5. tau_day
    6. lam_corr    (wavelength correlation, A)
    7. epoch_corr  (epoch correlation, days)

    LAMRANGE[2] is the wavelength range for which the smearing
    model will be used.

    The random mag-nodes are splined together with a sin function
    so that the mag-function derivative is 0 at each node.

    Returns number of Gaussian randoms to generate for each SN.

    Apr 14 2014: fix bug setting LAMDIF (bug found by Gautham)
  */

  int    NRANGEN, i   ;
  double LAMDIF ;

  char fnam[] = "init_genSmear_USRFUN" ;

  // -------------- BEGIN -----------

  if ( NPAR != NPAR_GENSMEAR_USRFUN_REQUIRED ) {
    sprintf(c1err,"passed NPAR=%d, but require NPAR=%d",
	    NPAR, NPAR_GENSMEAR_USRFUN_REQUIRED ) ;
    sprintf(c2err,"%s", "Check parameters");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


  // load global struct
  GENSMEAR_USRFUN.USE   = 1 ;  GENSMEAR.NUSE++ ;

  i = -1 ;
  
  i++ ; GENSMEAR_USRFUN.SIGCOH       = parList[i] ;
  i++ ; GENSMEAR_USRFUN.A5500        = parList[i] ;
  i++ ; GENSMEAR_USRFUN.LAMINTERVAL  = parList[i] ;
  i++ ; GENSMEAR_USRFUN.LAMPHASE     = parList[i] ;
  i++ ; GENSMEAR_USRFUN.TAU_LAM      = parList[i] ;
  i++ ; GENSMEAR_USRFUN.TAU_EPOCH    = parList[i] ;
  i++ ; GENSMEAR_USRFUN.CORR_LAM     = parList[i] ;
  i++ ; GENSMEAR_USRFUN.CORR_EPOCH   = parList[i] ;


  GENSMEAR_USRFUN.LAMREF = 5500.0 ;  // hard wired

  GENSMEAR_USRFUN.LAMRANGE[0] = LAMRANGE[0];
  GENSMEAR_USRFUN.LAMRANGE[1] = LAMRANGE[1];

  // first get number of randoms to generate.
  // A random for each lambda interval, plus a spare.

  //  LAMDIF  = LAMRANGE[1] - LAMRANGE[0] ; 
  LAMDIF  = LAMRANGE[1]  ; // fixed Apr 14 2014
  NRANGEN = (int)(LAMDIF / GENSMEAR_USRFUN.LAMINTERVAL) + 3 ;

  /*
  printf(" xxx LAMDIF = %f - %f = %f \n",
	 LAMRANGE[1], LAMRANGE[0], LAMDIF);
  printf(" xxx LAMINTERVAL = %f \n", GENSMEAR_USRFUN.LAMINTERVAL) ;
  printf(" xxx NRANGEN=%d \n", NRANGEN); fflush(stdout);
  */


  // add random number for (optional) random phase
  if ( GENSMEAR_USRFUN.LAMPHASE <= 0.00001 ) { NRANGEN++ ; }

  // add random for coherent part
  NRANGEN++ ;

  // store NRANGEN in global struct
  init_genSmear_randoms(NRANGEN,0);

  sprintf(BANNER,"%s", fnam);
  print_banner(BANNER);

  printf("\t Coherent = %7.3f \n", 
	 GENSMEAR_USRFUN.SIGCOH );
  printf("\t Amp5500 = %6.3f * EXP(-LAM/%6.0f) * EXP(-Trest/%5.1f) \n"
	 ,GENSMEAR_USRFUN.A5500
	 ,GENSMEAR_USRFUN.TAU_LAM
	 ,GENSMEAR_USRFUN.TAU_EPOCH );
  printf("\t LAMBDA(phase,interval) = %7.0f, %7.0f \n"
	 ,GENSMEAR_USRFUN.LAMPHASE
	 ,GENSMEAR_USRFUN.LAMINTERVAL );

  printf("\t Requires %d Gaussian random numbers per SN \n", NRANGEN ) ;
  printf("\t Function defined from %7.0f < LAMBDA < %7.0f \n",
	 LAMRANGE[0], LAMRANGE[1] );

  if ( NRANGEN >= MXRAN_GENSMEAR ) {
    sprintf(c1err,"NRANGEN=%d exceeds bound of %d",
	    NRANGEN, MXRAN_GENSMEAR ) ;
    sprintf(c2err,"%s", "Check FUNPAR values and  LAMRANGE.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }


} // end of init_genSmear_USRFUN



void get_NRAN_genSmear(int *NRANGauss, int *NRANFlat) {
  // interface to external calling program to tell how 
  // many randoms to generate.
  *NRANGauss = GENSMEAR.NGEN_RANGauss ;
  *NRANFlat  = GENSMEAR.NGEN_RANFlat ;
} 




void SETSNPAR_genSmear(double shape, double color, double redshift) {

  // Jan 2014: 
  // SN parameters passed from simulation or external program
  // Just store these parameters in case they are needed by
  // one of the get_genSmear_XXX functions.
  
  // ------------- BEGIN -----------

  GENSMEAR.SHAPE    = shape ;
  GENSMEAR.COLOR    = color ;
  GENSMEAR.REDSHIFT = redshift ;

  // store powers of redshift 
  int i;  double zpow ;
  zpow = 1.0 ;
  for(i=0; i < 8; i++ ) {
    GENSMEAR.ZPOW[i] = zpow ;
    zpow *= redshift;
  }
 

  // -------------------------------------
  // do once-per-SN tasks
  if ( GENSMEAR_CCM89.USE )  { GENSMEAR_CCM89.RV = get_CCM89_RV(); }

}  // end of SETSNPAR_genSmear


// ************************************
double get_genSmear_SCALE(double *parList) {

  // Created Mar 22 2020 by R.Kessler
  // Return magSmear scale based on inpout parList:
  //   parList[0:3] = Trest, x1, c, logMass

  double x1 = parList[1];
  double c  = parList[2];
  double m  = parList[3];

  double SCALE  = 1.0 ;
  double VAL;
  char *varName = GENSMEAR_SCALE.VARNAME ;
  char fnam[] = "get_genSmear_SCALE" ;

  // -------------- BEGIN ------------

  if ( !GENSMEAR_SCALE.USE ) { return(SCALE); }

  // xxxxxxxx DEVEL option, May 31 2021 xxxxxx
  if ( strcmp(varName,"DEVEL") == 0 ) {
    SCALE = get_genSmear_devel(parList);
    return(SCALE);
  }
  // xxxxxxxxxxxxxxxxx

  if ( strcmp(varName,"c") == 0 ) 
    { VAL = c; }
  else if ( strcmp(varName,"x1") == 0 ) 
    { VAL = x1 ; }
  else {
    VAL = -9999.0 ;
    sprintf(c1err,"Invalid VARNAME='%s' for poly fun", varName);
    sprintf(c2err,"Check sim-input GENMAG_SMEAR_SCALE");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  SCALE = eval_GENPOLY(VAL, &GENSMEAR_SCALE.POLY, fnam);

  return(SCALE);

} // end get_genSmear_SCALE


// ************************************************
void get_genSmear_USRFUN(double Trest, int NLam, double *Lam, 
			 double *magSmear) {

  // Created Mar 6, 2012 by R.Kessler
  // Return magSmear array in magnitudes (of length NLam).
  // Inputs are the rest-frame epochs in days, with Trest=0 at peak
  // and wavelength array (*Lam) in Angstroms.

  int Nran, ilam, NNODE, INODE, LAST_INODE  ;
  double 
    SIGCOH
    ,A5500
    ,LAMRANGE[2]
    ,LAMPHASE
    ,LAMINTERVAL, LAMDIF, MAGCOH
    ,LAM, lam, rr
    ,LAM_NODE[2], MAG_NODE[2]
    ,TAU_LAM, TAU_EPOCH
    ,EXP_LAM, EXP_EPOCH, AMPL
    ;

  char fnam[] = "get_genSmear_USRFUN" ;

  // ------------- BEGIN -------------------

  // init all magSmear values to zero.
  for ( ilam=0; ilam < NLam; ilam++ ) { magSmear[ilam] = 0.0 ; }

  if ( GENSMEAR_USRFUN.USE == 0 ) { return ; }

  SIGCOH      = GENSMEAR_USRFUN.SIGCOH ;
  A5500       = GENSMEAR_USRFUN.A5500 ;
  LAMINTERVAL = GENSMEAR_USRFUN.LAMINTERVAL ;
  TAU_LAM     = GENSMEAR_USRFUN.TAU_LAM ;
  TAU_EPOCH   = GENSMEAR_USRFUN.TAU_EPOCH ;
  LAMRANGE[0] = GENSMEAR_USRFUN.LAMRANGE[0];
  LAMRANGE[1] = GENSMEAR_USRFUN.LAMRANGE[1];

  Nran  = -1 ; // random counter

  Nran++ ; rr = GENSMEAR.RANGauss_LIST[Nran] ;
  MAGCOH = rr * SIGCOH ;


  // determine wavelength phase.
  if ( GENSMEAR_USRFUN.LAMPHASE > 0 ) 
    { LAMPHASE = GENSMEAR_USRFUN.LAMPHASE ; }
  else {
    Nran++ ; rr = GENSMEAR.RANGauss_LIST[Nran] ;
    LAMPHASE = LAMINTERVAL * rr ;
  }


  // Now keep subtracting or adding LAMINTERVCAL from LAMPHASE
  // until we are below LAMRANGE[0] and still above 0.

  while ( LAMPHASE > LAMRANGE[0]  &&   LAMPHASE > LAMINTERVAL ) {
    LAMPHASE -= LAMINTERVAL ;
  }

  // loop over lambda intervals (nodes) and set random magSmear 

  NNODE = 0 ;


  for ( LAM=LAMPHASE; LAM < LAMRANGE[1]+LAMPHASE; LAM += LAMINTERVAL ) {

    if ( LAM < LAMRANGE[0] - LAMPHASE ) { continue ; }

    NNODE++ ;
    Nran++ ; rr = GENSMEAR.RANGauss_LIST[Nran] ;

    LAMDIF    = LAM - GENSMEAR_USRFUN.LAMREF ;
    EXP_LAM   = exp(-LAMDIF/TAU_LAM) ;
    EXP_EPOCH = exp(-Trest/TAU_EPOCH) ;
    AMPL  = MAGCOH + rr * A5500 * EXP_LAM * EXP_EPOCH ;

    GENSMEAR_USRFUN.AMPL_NODE[NNODE-1] = AMPL ;      
    GENSMEAR_USRFUN.LAM_NODE[NNODE-1]  = LAM ;

    /*    
    printf(" xxxxx Ampl[%2d: LAM=%6.0f] = %f \n", 
	   NNODE, LAM, AMPL );
    */
    
    if ( Nran > GENSMEAR.NGEN_RANGauss ) {
      sprintf(c1err,"Nran=%d, but expected list has only %d",
	      Nran, GENSMEAR.NGEN_RANGauss ) ;
      sprintf(c2err,"Something is really fishy.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
  }

  GENSMEAR_USRFUN.NNODE = NNODE;

  LAST_INODE = -99;


  // loop over user-passed lambda array
  for ( ilam=0; ilam < NLam; ilam++ ) {
    lam = Lam[ilam];
    
    // find nodes that  bound this 'lam'
    INODE = INODE_LAMBDA(lam, NNODE, GENSMEAR_USRFUN.LAM_NODE );

    if ( INODE < 0 ) {
      sprintf(c1err,"Could not find INODE for lam=%7.1f", lam);
      sprintf(c2err,"NNODE=%d  LAMPHASE=%7.1f LAMINTERVAL=%7.1f",
	      NNODE, LAMPHASE, LAMINTERVAL);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }


    if ( LAST_INODE != INODE ) {
      LAM_NODE[0] = GENSMEAR_USRFUN.LAM_NODE[INODE] ;
      LAM_NODE[1] = GENSMEAR_USRFUN.LAM_NODE[INODE+1] ;
      MAG_NODE[0] = GENSMEAR_USRFUN.AMPL_NODE[INODE];
      MAG_NODE[1] = GENSMEAR_USRFUN.AMPL_NODE[INODE+1];
    }
    
    // note that the sin(arg) function slows the
    // generation by about 50% !!!

    magSmear[ilam] = interp_SINFUN(lam, LAM_NODE, MAG_NODE, fnam );
    
    /*
    printf(" xxxx %4.4d : lam=%7.1f  magSmear(old,new) = %f/%f \n", 
	   ilam, Lam[ilam], magSmear[ilam], mtest ); // DDDDDDDDDDD
    */

    LAST_INODE = INODE ;

  } // ilam


}  // end of get_genSmear_USRFUN


// ********************************
void init_genSmear_CCM89(double *LAMRANGE) {

  // LAMRANGE[2] is the input wavelength range to
  // store the CCM map. The MAP speeds up generation by
  // about 50%, but this model still slows down the
  // overall generation speed by about a factor of
  // 2 compared to no smear model.

  int    NLAM, I8 ;
  double sig, REF, LAM, LAMBIN;
  char fnam[] = "init_genSmear_CCM89" ;

  // --------- BEGIN ------------

  GENSMEAR_CCM89.USE  = 1 ;  GENSMEAR.NUSE++ ;

  GENSMEAR_CCM89.SIGCOH = 0.09 ;

  REF = 2.5 ;
  sig = 0.6 ;
  GENSMEAR_CCM89.RV_REF   = REF ;
  GENSMEAR_CCM89.RV_sigma = sig ;
  GENSMEAR_CCM89.RV_range[0] = REF - 2.0*sig ; // sigma clip
  GENSMEAR_CCM89.RV_range[1] = REF + 2.0*sig ;

  printf("   Init CCM89-ColorLaw smear: <RV>=%3.1f  sigma_RV=%3.1f \n",	 
	 REF, sig );
  printf("\t SIGCOH = %5.2f \n", GENSMEAR_CCM89.SIGCOH);

  // allocate map to store magsmear vs. lambda at start of each SN.
  I8     = sizeof(double);
  NLAM   = 0 ;
  LAMBIN = 20.0 ; // hard-wire 20 A bins
  NLAM   = (int)((LAMRANGE[1] - LAMRANGE[0])/LAMBIN) + 4 ; // approx
  GENSMEAR_CCM89.LAM_MAP      = (double*)malloc(NLAM*I8);
  GENSMEAR_CCM89.XTDIF_MAP    = (double*)malloc(NLAM*I8);

  NLAM = 0; // reset NLAM
  for ( LAM = LAMRANGE[0]; LAM < LAMRANGE[1]+1.0; LAM+=LAMBIN ) {
    NLAM++ ;
    GENSMEAR_CCM89.LAM_MAP[NLAM-1] = LAM ;
  }
  GENSMEAR_CCM89.NLAM_MAP = NLAM;

  printf("    CCM89 map allocated for %d LAMBDA bins (%6.0f - %6.0f)\n",
	 NLAM, LAMRANGE[0], LAMRANGE[1]);

  // allocate Gaussian random numbers for SIGCOH(1) and random RV(2).
  init_genSmear_randoms(2,0);

} // end of init_genSmear_CCM89


// ***********************************
double get_CCM89_RV(void) {
  
  // called from SETRAN to return random RV for this SN.
  // Also fill MAGSMEAR_MAP for CPU speed.

  double RVmin, RVmax, RV0, RVsig, RV, RVsmear, ran ;
  char fnam[] = "get_CCM89_RV" ;

  // ---------- BEGIN ----------

  ran    = GENSMEAR.RANGauss_LIST[1] ;
  RV0    = GENSMEAR_CCM89.RV_REF ;
  RVsig  = GENSMEAR_CCM89.RV_sigma ;
  RVmin  = GENSMEAR_CCM89.RV_range[0];
  RVmax  = GENSMEAR_CCM89.RV_range[1];

  RVsmear = ran*RVsig ;
  RV   = RV0 + RVsmear ;

  // to apply clipping, keep dividing the fluctuation by 2
  // until it's withing range ... 
  // avoids burning extra random numbers.

  while ( RV > RVmax || RV < RVmin ) {
    RVsmear /= 2.0 ;
    RV   = RV0 + RVsmear ;
  }

  // load/store magSmear corrections at all lambda

  int ilam;
  int OPT = 89 ;
  double PARDUM=0.0 ;
  double LAM, Color, XT, XT0, XTDIF ;
  Color = GENSMEAR.COLOR ;  // SALT2 'c' parameter

  for ( ilam=0; ilam < GENSMEAR_CCM89.NLAM_MAP; ilam++ ) {
    LAM  = GENSMEAR_CCM89.LAM_MAP[ilam] ;
    XT   = GALextinct(RV,  Color, LAM, OPT, &PARDUM, fnam ) ;
    XT0  = GALextinct(RV0, Color, LAM, OPT, &PARDUM, fnam ) ;
    XTDIF = XT - XT0 ;
    GENSMEAR_CCM89.XTDIF_MAP[ilam]  = XTDIF ;
  }

  return RV ;

} // end of get_CCM89_RV
 
// ******************************
void get_genSmear_CCM89(double Trest, int NLam, double *Lam, 
			double *magSmear) {

  double  lam, XTDIF, SMEAR0 ;
  double *ptrLam, *ptrXTDIF ;

  int ilam, NLAM;

  char fnam[] = "get_genSmear_CCM89";

  // ---------- BEGIN ----------

  NLAM     = GENSMEAR_CCM89.NLAM_MAP ;
  ptrLam   = GENSMEAR_CCM89.LAM_MAP ;
  ptrXTDIF = GENSMEAR_CCM89.XTDIF_MAP ;

  // get coherent piece
  SMEAR0 = GENSMEAR_CCM89.SIGCOH * GENSMEAR.RANGauss_LIST[0] ;

  for ( ilam=0; ilam < NLam; ilam++ ) {
    lam  = Lam[ilam];
    // interpolate to get magShifts.
    XTDIF  = interp_1DFUN(1, lam, NLAM, ptrLam, ptrXTDIF, fnam );
    magSmear[ilam] = SMEAR0 + XTDIF ;
  }

} // end of get_genSmear_CCM89


// ***********************************************
void  init_genSmear_SALT2(char *versionSALT2, char *smearModel, 
			  double SIGCOH, double *GENRANGE_REDSHIFT) {

  // Created May 1 2014 by R.Kessler
  // General init for SALT2 intrinsic scatter.
  // [replaces init_genSmear_G10* functions]
  //
  // Inputs:
  //   versionSALT2 = SALT2 version, including PATH
  //   smearModel = name of smear model (G10, G10Fig8, G10FUDGE ..)
  //   SIGCOH     = coherent luminosity scatter; if negative use default.
  //   GENRANGE_REDSHIFT = zmin, zmax (to check for errors)
  //
  // Aug 26 2015: 
  //  Sort'of fix problem when MAXLAM(NODE) < MAXLAM(SED) and it fails
  //  later when lambda is between the two MAXLAM values.
  //  The fix here is to define an extra NODE at MAXLAM(SED), even though
  //  the lambda spacing to the last node is not correct. It's a kluge
  //  to prevent aborting.
  //
  // May 29 2018: input versionSALT2 now includes full path
  //
  // May 30 2018: check SIGCOH=-8 --> command line override of SIGMA_COH
  //
  // Apr 11 2019:
  //   Tried to work with BYOSED, but two problems:
  //    1. SIGCOH cannot be read from SALT2.INFO file
  //    2. dispFile unknown.
  //  Need to pass a full SALT2 model path to work with BYOSED.
  //
  // + remove call to getFileName_SALT2colorDisp(dispFile), and
  //   construct dispFile using versionSALT2 ... needed to work
  //   for BYOSED & SALT2 models.
  //
  // Aug 28 2019:
  //  + add GENRANGE_REDSHIFT argument, and abort on potential gen-errors.
  //  + add logic for last NODE bin
  //
  // Nov 4 2020: allow for salt2 or salt3 in file name for 
  //             salt2[3]_color_dispersions
  //
  // Dec 7 2020: for SIGCOH>0, call read_genSmear_SALT2sigcoh
  //             (fix needed for G10 scatter model with BYOSED)
  //
  // Dec 28 2020: check for SALT3
  // Mar 29 2023: flag crazy SIGMA values and abort.

  double SED_LAMMIN =  SALT2_TABLE.LAMMIN;
  double SED_LAMMAX =  SALT2_TABLE.LAMMAX;  
  double COLOR_DISP_MAX = 5.0;
  double zmin = GENRANGE_REDSHIFT[0];
  int   NSIGMA_CRAZY = 0;
  double SIGMA_CRAZY = 2.0 ;
  char dispFile[MXPATHLEN] ;  
  char fnam[] = "init_genSmear_SALT2" ;

  // ----------------- BEGIN --------------

  GENSMEAR_SALT2.USE  = 1 ;    GENSMEAR.NUSE++ ;
  GENSMEAR_SALT2.SIGCOH  = 0.0 ;

  // note that SIGCOH=-8 is a flag that user passed SIGMA_COH
  // via command-line override, so ignore this key in SALT2.INFO file.
  if ( SIGCOH != -8.0 ) { GENSMEAR_SALT2.SIGCOH_LAM.NBIN = 0 ; }

  if ( strstr(smearModel,"G10Fig8") != NULL ) {
    GENSMEAR_SALT2.SIGCOH = 0.0 ;
    GENSMEAR_SALT2.FUDGE_LAMSWITCH      = 0.0 ;
    GENSMEAR_SALT2.FUDGE_dSmear_dLAM[0] = 0.0 ; // lam < 2800
    GENSMEAR_SALT2.FUDGE_dSmear_dLAM[1] = 0.0 ; // lam > 2800
  }
  else {
    GENSMEAR_SALT2.FUDGE_LAMSWITCH      =  2157.3 ;
    GENSMEAR_SALT2.FUDGE_dSmear_dLAM[0] =  0.0 ;     // lam < LAMSWITCH
    GENSMEAR_SALT2.FUDGE_dSmear_dLAM[1] = +1.08E-4 ; // lam > LAMSWITCH
  }

  // Dec 28 2020 : check if SALT2 or SALT3
  // setFlags sets global logicals ISMODEL_SALT2 & ISMODEL_SALT3
  setFlags_ISMODEL_SALT2(versionSALT2); //

  if ( ISMODEL_SALT2 ) {
    sprintf(SALT2_INFO_FILE,"SALT2.INFO"); 
    sprintf(dispFile, "%s/salt2_color_dispersion.dat", versionSALT2);
  }
  else if ( ISMODEL_SALT3 ) {
    sprintf(SALT2_INFO_FILE,"SALT3.INFO"); 
    sprintf(dispFile, "%s/salt3_color_dispersion.dat", versionSALT2);
  }
  else {
    sprintf(c1err,"Unknown ISMODEL for version = ");
    sprintf(c2err,"%s", versionSALT2);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // --------------------------------------
  // read SIGMA_INT from SALT2.INFO file
  
  if ( SIGCOH < -8.1 )  {  
    read_genSmear_SALT2INFO(versionSALT2, &GENSMEAR_SALT2.SIGCOH_LAM,
			    &COLOR_DISP_MAX); 
  }
  else if ( SIGCOH >= 0.0 ) {
    read_genSmear_SALT2INFO(versionSALT2, &GENSMEAR_SALT2.SIGCOH_LAM,
			    &COLOR_DISP_MAX); 
  }
  else {
    // do nothing because user input includes
    //   SIGMA_COH(lamList):   sigcoh-list
    // to override whatever is in SALT2.INFO file
  }

  // --------------------------------------------
  // read dispersion vs. lambda : 
  // Use already-determine filename that allows for salt2 or salt3

  read_genSmear_SALT2disp(dispFile, COLOR_DISP_MAX) ;

  // ----

  int    NRANGEN, ilam, NLAM, NNODE, LAST=0 ;
  double LAM, LAM2=0.0, MINLAM, MAXLAM, MAXLAM_LOOP, DLAM, SIG ;
  double F0, F1, FUDGE, L0;

  // apply sigma-scale fudge so that lambda-interpolated model
  // gives correct broadband  magsmear at all wavelens. 
  // Note that Fudge is 1.0 at L0, and increases linearly 
  // for smaller and larger lambda.

  NLAM = GENSMEAR_SALT2.NLAM ;
  F0 = GENSMEAR_SALT2.FUDGE_dSmear_dLAM[0];
  F1 = GENSMEAR_SALT2.FUDGE_dSmear_dLAM[1];
  L0 = GENSMEAR_SALT2.FUDGE_LAMSWITCH ;
  for(ilam=0; ilam < NLAM; ilam++ ) {
    LAM = GENSMEAR_SALT2.LAM[ilam] ; 
    if ( LAM < L0 ) 
      { FUDGE = 1.0 + (LAM-L0)*F0 ;  }
    else 
      { FUDGE = 1.0 + (LAM-L0)*F1 ;  }

    GENSMEAR_SALT2.SIGMA[ilam] *= FUDGE ;
  }


  // print info
  int NBIN = GENSMEAR_SALT2.SIGCOH_LAM.NBIN ;
  char tmpLine[100];

  if ( NBIN == 1 ) {
    sprintf(tmpLine, "SIGCOH = %5.3f",  GENSMEAR_SALT2.SIGCOH_LAM.YVAL[0] );
  }
  else {
    sprintf(tmpLine, "SIGCOH = " );
    for(ilam=0; ilam < NBIN; ilam++ ) 
      { sprintf(tmpLine, "%s %5.3f", 
		tmpLine, GENSMEAR_SALT2.SIGCOH_LAM.YVAL[ilam] ); }
    sprintf(tmpLine, "%s for LAM =", tmpLine );
    for(ilam=0; ilam < NBIN; ilam++ ) 
      { sprintf(tmpLine, "%s %.1f", 
		tmpLine, GENSMEAR_SALT2.SIGCOH_LAM.XVAL[ilam] );} 
  }
  printf("    %s\n", tmpLine);

  if ( L0 > 0.0 ) {
    printf("    SMEAR *= [ 1 + %9.6f*(LAM-%4.0f)]   LAM < %4.0f \n",
	   F0, L0, L0 );
    printf("    SMEAR *= [ 1 + %9.6f*(LAM-%4.0f)]   LAM > %4.0f \n",
	   F1, L0, L0 );
    fflush(stdout);
  }


  // count randoms
  NRANGEN = 0 ;

  NRANGEN++ ; // for SIGCOH
  
  // Count new random every 800 A and store SIGMAs
  GENSMEAR_SALT2.LAMSEP_NODE = 800.0 ; // lambda node separation (A)
  DLAM   = GENSMEAR_SALT2.LAMSEP_NODE ; // temp variable
  NNODE  = 0 ;
  MINLAM = GENSMEAR_SALT2.MINLAM ;
  MAXLAM = GENSMEAR_SALT2.MAXLAM ;
  MAXLAM_LOOP = MAXLAM + DLAM ; 

  // Apr 2022
  // make sure that lam-range in salt2_color_dispersion.dat covers
  // lam-range of SED template. Beware UV extrap cause failure on blue side.
  //  if ( MINLAM > SED_LAMMIN || MAXLAM < SED_LAMMAX ) {
  if (  MAXLAM < SED_LAMMAX ) {
    print_preAbort_banner(fnam);
    char cdisp_file[] = "salt2_color_dispersion.dat";
    printf("  %s MINLAM / MAXLAM = %.2f / %.2f\n", 
	   cdisp_file, MINLAM, MAXLAM);
    printf("  SALT2 SED MINLAM / MAXLAM = %.2f / %.2f\n", 
	   SED_LAMMIN, SED_LAMMAX);

    sprintf(c1err,"wave-range for %s is too small;", cdisp_file); 
    sprintf(c2err,"Check MINLAM/MAXLAM ranges above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );    
  }

  for ( LAM = MINLAM ; LAM <= MAXLAM; LAM+=DLAM ) {

    LAM2 = LAM ;
    // protect bounds
    if ( LAM < MINLAM ) { LAM2 = MINLAM+0.001 ; }

    SIG = interp_1DFUN( 1, LAM2
			,GENSMEAR_SALT2.NLAM
			,GENSMEAR_SALT2.LAM 
			,GENSMEAR_SALT2.SIGMA
			,fnam );    

    if ( SIG > COLOR_DISP_MAX_DEFAULT ) {
      SIG = COLOR_DISP_MAX_DEFAULT - 1.0E-8 ;
    }


    printf("\t Set LAM-node %2d at %7.1f A : SIGMA=%6.3f \n", 
	   NNODE, LAM2, SIG ); fflush(stdout);

    GENSMEAR_SALT2.LAM_NODE[NNODE] = LAM2 ; // lambda at each node
    GENSMEAR_SALT2.SIG_NODE[NNODE] = SIG  ; // sigma at each node.

    if ( NNODE > 80 ) {
      sprintf(c1err,"NNODE=%d is way too big. What the ?@!?", NNODE);
      sprintf(c2err,"GENSMEAR_SALT2.LAM[min,max] = %7.1f  %7.1f",
	      GENSMEAR_SALT2.MINLAM, GENSMEAR_SALT2.MAXLAM );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    NNODE++ ;    NRANGEN++ ;

    if ( MAXLAM-LAM < DLAM  && LAM < SED_LAMMAX && !LAST )
      { LAST=1; DLAM = MAXLAM-LAM-.001; } // Sep 3 2019

  } // end LAM loop
 
  GENSMEAR_SALT2.NNODE = NNODE ;

  if ( NSIGMA_CRAZY > 0 ) {
    sprintf(c1err,"Found %d crazy SIGMA values", NSIGMA_CRAZY);
    sprintf(c2err,"Either trim rest-frame wave range or fix model dispersion");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  // Aug 28 2019: make sure last node is covered by wavelength range
  double LAMCHECK = LAM2 * (1.0+zmin);
  if ( LAMCHECK < MAXLAM && LAMCHECK < SED_LAMMAX ) { // Sep 19 2019
    double zmin_suggest = MAXLAM/LAM2 - 1.0 ;
    sprintf(c1err,"genSmear model can't handle zmin=%.3f "
	    "(see GENRANGE_REDSHIFT)", zmin);
    sprintf(c2err,"Suggest increasing zmin to be > %.3f\n", zmin_suggest);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  // -------------------

  // store number of randoms in global.
  init_genSmear_randoms(NRANGEN,0);

  return ;

} // end of init_genSmear_SALT2


// ********************************************
void read_genSmear_SALT2INFO(char *versionSALT2, GRIDMAP1D_DEF *SIGCOH_LAM,     
			     double  *COLOR_DISP_MAX) 
{ 

  // Minor refactor Oct 21 2022 to read and return COLOR_DISP_MAX.
  // Read SALT2[3].INFO file for information needed for genSmear.
  // Input *versionSALT2 is the name of the SALT2/SALT3 version,
  // including full path
  //
  // Fix INFO_FILE to work with SALT2 or SALT3

  FILE *fp ;
  char INFO_FILE[MXPATHLEN], *keyName, c_get[60], keyArg[80] ;
  char KEY1[]  = "SIGMA_INT" ;
  char KEY2[]  = "SIGMA_COH" ;
  char fnam[] = "read_genSmear_SALT2sigcoh" ;

  // --------- BEGIN -------

  // xxx  sprintf(INFO_FILE, "%s/SALT2.INFO", versionSALT2);
  sprintf(INFO_FILE, "%s/%s", versionSALT2, SALT2_INFO_FILE ); // 11.04.2020


  if ( (fp = fopen(INFO_FILE, "rt")) == NULL ) { 
      sprintf(c1err,"Cannot open infoFile " );
      sprintf(c2err,"%s", INFO_FILE);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;    
  }

  SIGCOH_LAM->NBIN = 0 ;  // SIGCOH vs. lam array

  while ( (fscanf(fp, "%s", c_get)) != EOF ) {

    if ( c_get[0] == '#' ) { continue; }
    if ( strstr(c_get,KEY1) != NULL || 	strstr(c_get,KEY2) != NULL ) {
      keyName = c_get;
      readchar(fp, keyArg); 
      parse_SIGCOH_SALT2(keyName,keyArg, SIGCOH_LAM);
    }

    if ( strcmp(c_get,"COLOR_DISP_MAX:") == 0 ) 
      { readdouble(fp, 1, COLOR_DISP_MAX); }
  }

  fclose(fp);
  
  // abort if SIGCOH not found.
  int NBIN = SIGCOH_LAM->NBIN ;
  if ( NBIN == 0  ){
    sprintf(c1err,"Could not find '%s' or %s in SALT2 info-file:", 
	    KEY1, KEY2);
    sprintf(c2err,"%s", INFO_FILE) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ; 
  }

  //  debugexit(fnam);

  return;

} // end of read_genSmear_SALT2sigcoh


// ************************************
void  parse_SIGCOH_SALT2(char *KEYNAME, char *KEYARG, 
			 GRIDMAP1D_DEF *SIGCOH_LAM) {

  // Created May 30 2018
  //
  // Parse inputs KEYNAME and KEYARG and fill output SIGCOH_LAM.
  // Examples:
  //   
  //   KEYNAME                       KEYARG
  //   SIGMA_INT:                   0.10
  //   SIGMA_INT(1000,9000,25000):  0.12,0.10,0.05
  //
  // Note that SIGMA_INT or SIGMA_COH are allowed keys.
  //
  // Output struct is
  //   SIGCOH_LAM->NBIN = number of bins in map
  //   SIGCOH_LAM->XVAL = wavelength array
  //   SIGCOH_LAM->YVAL = sigma_coh array
  //

#define MXLAMBIN_SIGCOH 10

  double SIGCOH;
  char key[40], stringLam[80], stringArg[80] ;
  char fnam[] = "parse_SIGCOH_SALT2" ;

  // --------------- BEGIN ----------------
  
  //  printf(" xxx %s for '%s'  '%s'  \n", fnam, KEYNAME, KEYARG);

  SIGCOH_LAM->XVAL = (double*) malloc ( MXLAMBIN_SIGCOH * sizeof(double) );
  SIGCOH_LAM->YVAL = (double*) malloc ( MXLAMBIN_SIGCOH * sizeof(double) );

  // check for lambda values in ()
  sprintf(key, "%s", KEYNAME);
  extractStringOpt(key,stringLam);
  if ( strlen(stringLam) == 0 ) { // KEYARG scalar SIGCOH
    sscanf( KEYARG, "%le", &SIGCOH);
    SIGCOH_LAM->NBIN    = 1;
    SIGCOH_LAM->XVAL[0] = 0.0 ; // no particular wavelength
    SIGCOH_LAM->YVAL[0] = SIGCOH ;
    return ;
  }


  char comma[] = "," ;
  char *ptrSplit[MXLAMBIN_SIGCOH], stringList[MXLAMBIN_SIGCOH][20];
  int  ilam, NLAM, NSIGCOH ;

  // parse comma-separated wavelengths in stringLam.  
  for(ilam=0; ilam < MXLAMBIN_SIGCOH; ilam++ ) 
    { ptrSplit[ilam] = stringList[ilam]; }

  splitString(stringLam, comma, fnam, MXLAMBIN_SIGCOH,
	      &NLAM, ptrSplit );          // <=== returned
  for(ilam=0; ilam < NLAM; ilam++ ) 
    { sscanf(ptrSplit[ilam], "%le", &SIGCOH_LAM->XVAL[ilam] ); }


  // parse comma-separated SIGCOH values in KEYARG
  sprintf(stringArg, "%s", KEYARG);
  splitString(stringArg, comma, fnam, MXLAMBIN_SIGCOH,
	      &NSIGCOH, ptrSplit );          // <=== returned
  for(ilam=0; ilam < NSIGCOH; ilam++ ) 
    { sscanf(ptrSplit[ilam], "%le", &SIGCOH_LAM->YVAL[ilam] ); }

  if ( NLAM != NSIGCOH ) {
    sprintf(c1err,"N(LAM)=%d but N(SIGCOH)=%d", NLAM, NSIGCOH );
    sprintf(c2err,"Check '%s %s' ", KEYNAME, KEYARG);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }
  else
    { SIGCOH_LAM->NBIN = NLAM; }


  int LDMP = 0 ;
  if ( LDMP ) {
    for(ilam=0; ilam < NSIGCOH; ilam++ ) {
      printf(" xxx ilam=%d : LAM = %8.1f   SIGCOH=%6.4f \n",
	     ilam, SIGCOH_LAM->XVAL[ilam], SIGCOH_LAM->YVAL[ilam] );	     
    }
  }


  return ;

} // end  parse_SIGCOH_SALT2


// ************************************
void read_genSmear_SALT2disp(char *smearFile, double COLOR_DISP_MAX ) {

  // Created May 1 2014
  // read dispersion vs. wavelength from file.
  // Abort if input *smearFile  == ""
  //
  // Dec 29 2017: use open_TEXTgz() to allow for gzipped file.
  // Aug 28 2019: fill GENSMEAR_SALT2.MINLAM[MAXLAM]
  // Sep 03 2019: GENSMEAR_SALT2.LAM[SIGMA] arrays start at 0, not 1
  // Oct 21 2022: apply COLOR_DISP_MAX passed as new arg

  char fnam[] = "read_genSmear_SALT2disp" ;
  FILE *fp ;
  int  ilam, NLAM,  GZIPFLAG ;

  // -------------- BEGIN ----------------

  if ( strlen(smearFile) == 0 ) {
    sprintf(c1err,"input smearFile = '' is not defined ???");
    sprintf(c2err,"Something is really wrong.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  fp = open_TEXTgz(smearFile, "rt", 0, &GZIPFLAG, fnam ) ;

  if ( fp == NULL ) {
      sprintf(c1err,"Cannot open smearFile " );
      sprintf(c2err,"%s", smearFile);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;    
  }
  fclose(fp);

  sprintf(GENSMEAR_SALT2.FILE,"%s", smearFile);

  printf("  Read SALT2-smear file: \n  '%s' \n", GENSMEAR_SALT2.FILE );
  
  rd2columnFile( GENSMEAR_SALT2.FILE,  MXLAM_GENSMEAR_SALT2 // input
		 ,&NLAM                       // returned
		 ,GENSMEAR_SALT2.LAM        // returned
		 ,GENSMEAR_SALT2.SIGMA, 0  );  // returned

  GENSMEAR_SALT2.NLAM   = NLAM;  
  GENSMEAR_SALT2.MINLAM = GENSMEAR_SALT2.LAM[0];
  GENSMEAR_SALT2.MAXLAM = GENSMEAR_SALT2.LAM[NLAM-1];

  // Oct 2022: put a cap on color smear, to avoid crazy mags in far UV
  for(ilam=0; ilam < NLAM; ilam++ ) {
    if ( GENSMEAR_SALT2.SIGMA[ilam] > COLOR_DISP_MAX ) 
      { GENSMEAR_SALT2.SIGMA[ilam] = COLOR_DISP_MAX ; }
  }


  if ( NLAM >= MXLAM_GENSMEAR_SALT2 ) {
    sprintf(c1err,"G10 NLAM=%d exceeds array bound of %d",
	    NLAM, MXLAM_GENSMEAR_SALT2 );
    sprintf(c2err,"Either increase MXLAM_GENSMEAR_SALT2 or reduce NLAM");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  printf("    Read %d SALT2-Smear values from lambda=%6.0f to %6.0f \n",
	 NLAM, GENSMEAR_SALT2.MINLAM, GENSMEAR_SALT2.MAXLAM );

  fflush(stdout);

  return ;

} // end of void read_genSmear_SALT2disp



// *******************************************
void get_genSmear_SALT2(double Trest, int NLam, double *Lam, 
			double *magSmear) {

  // Returns magSmear array corresponding to input *Lam array.
  // Jul 27 2016: protect against undefined model (see LAMMIN, LAMMAX)
  // May 30 2018: interpolate SIGCOH vs. LAMBDA (see NBCOH)

  int NBCOH         = GENSMEAR_SALT2.SIGCOH_LAM.NBIN ;
  double *ptrLAMCOH = GENSMEAR_SALT2.SIGCOH_LAM.XVAL ;
  double *ptrSIGCOH = GENSMEAR_SALT2.SIGCOH_LAM.YVAL ;
  int OPT_INTERP=1;

  int    ilam, INODE, N, LDMP, NLAM ;
  double lam, rCOH, r0, r1, SMEAR0=0.0, SMEAR, MINLAM, MAXLAM ;
  double LAM_NODE[2], MAG_NODE[2];
  char   fnam[] = "get_genSmear_SALT2";
  char   comment[40];

  // ---------- BEGIN ----------

  LDMP = ( fabs(Trest) > 2.0E9 ) ;
  //  LDMP = ( fabs(Trest) < 5.0 ); /// xxx REMOVE

  if ( LDMP ) {
    printf(" xxx %s : Trest=%4.1f  NLAM=%d  lam=%7.1f  SIGCOH=%.3f\n",
	   fnam, Trest, NLam, Lam[0], ptrSIGCOH[0] );
  }

  NLAM   = GENSMEAR_SALT2.NLAM ;  
  MINLAM = GENSMEAR_SALT2.MINLAM ;
  MAXLAM = GENSMEAR_SALT2.MAXLAM ;
  

  N      = GENSMEAR.NGEN_RANGauss;
  rCOH   = GENSMEAR.RANGauss_LIST[N-1] ;  // use last random 

  if ( NBCOH == 1 ) {
    SMEAR0 = rCOH * ptrSIGCOH[0];
    GENSMEAR.MAGSMEAR_COH[0] = SMEAR0; // load global for SNTABLE (Jun 14 2016)
  }


  for ( ilam=0; ilam < NLam; ilam++ ) {
    lam = Lam[ilam];    
    if ( NBCOH >  1111111 ) {
      SMEAR0 = rCOH * interp_1DFUN(OPT_INTERP, lam, NBCOH,
				   ptrLAMCOH, ptrSIGCOH, fnam) ; 
    }

    magSmear[ilam] = SMEAR0 ;

    if ( lam <= (MINLAM+0.001) ) { continue ; }
    if ( lam >= (MAXLAM-0.001) ) { continue ; }

    INODE = INODE_LAMBDA(lam, GENSMEAR_SALT2.NNODE, GENSMEAR_SALT2.LAM_NODE);
    if ( INODE < 0 || INODE >= GENSMEAR_SALT2.NNODE ) {      
      print_preAbort_banner(fnam); 
      printf("  MINLAM / MAXLAM = %.2f / %.2f \n", MINLAM, MAXLAM);
      printf("  ilam = %d of %d   Trest=%.2f\n",  ilam, NLam, Trest);
      sprintf(c1err,"Could not find INODE for lam=%7.1f", lam);
      sprintf(c2err,"NNODE=%d  INODE=%d", GENSMEAR_SALT2.NNODE, INODE) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    LAM_NODE[0] = GENSMEAR_SALT2.LAM_NODE[INODE];
    LAM_NODE[1] = GENSMEAR_SALT2.LAM_NODE[INODE+1];

    r0          = GENSMEAR.RANGauss_LIST[INODE] ;
    r1          = GENSMEAR.RANGauss_LIST[INODE+1] ;

    MAG_NODE[0] = r0 * GENSMEAR_SALT2.SIG_NODE[INODE] ;
    MAG_NODE[1] = r1 * GENSMEAR_SALT2.SIG_NODE[INODE+1] ;

    sprintf(comment,"%s[INODE=%d]", fnam, INODE);
    SMEAR = interp_SINFUN(lam, LAM_NODE, MAG_NODE, comment );

    magSmear[ilam] = SMEAR0 + SMEAR ;

    // DDDDDDDDDDDDDDDDDDD
    if ( fabs(lam+3000.0)<0.0 && fabs(Trest) < 2.0 ) {
      printf(" xxx ----------------------------------------- \n" );
      printf(" xxx lam=%7.1f  SMEAR=%f  (r0,r1=%5.3f,%5.3f) \n", 
	     lam, SMEAR, r0,r1 );
      printf(" xxx \t NODE0: LAM=%7.1f SIG=%6.3f\n", 
	     LAM_NODE[0], MAG_NODE[0] );
      printf(" xxx \t NODE1: LAM=%7.1f SIG=%6.3f\n", 
	     LAM_NODE[1], MAG_NODE[1] );
    }
    // DDDDDDDDDDDDDD


  } // end of ilam loop

  return ;

} // end of  get_genSmear_SALT2


// ************************************
void  init_genSmear_Chotard11(int OPT_farUV) {

  // Created Apr 06, 2012
  // Based on UBVRI [anti]-correlations.
  // Copy/modify Cholesky stuff from INIT_COVMAT_SCATTER
  //
  // Hard-wire reduced matrix from Chotard's thesis, page 185
  // Use upper triangle corresponding to EWSi 4131 only
  // since Julien claims this best matches the SALT2 stretch corr.
  //
  // Chotard's 2-digit rounding leads to a gsl matrix error of
  // 'not positive definite'. Instead we use John M's matrix
  // that is slightly tweaked, along with adding 1.0E-9 to 
  // the diagonal covariances (see COV_DIAG_FUDGE).
  //
  // Chotard gives a passband covariance model, but here we translate 
  // into a wavelength model by picking the 5 correlated randoms at
  // the central UBVRI wavelengths, and then splining the points
  // together with a SIN function (interp_SINFUN). 
  // Lastly, COV_SCALE multiplies all covariances so that the 
  // wavelength-dependent model gives the more-or-less correct 
  // passband covariances. For example, the input/recovered diagonal 
  // terms from a simulation are 
  // RMS-UBVRI(input/recoverd) = 
  //    0.06/0.0586  0.04/0.0344   0.05/0.0492   0.04/0.0408  0.08/0.0783
  //
  // May 23, 2012: far-UV update.
  //  Include far UV filter (v) at 2500 A using G10 sigma(2500)=0.59
  //               with 3 options determined by input OPT_farUV
  //               0 -> no correlation
  //               1 -> 100% correlation with U
  //               2 -> 100% anti-correlation with U
  //
  //
  // Feb 17 2020: refactor to use CHOLESLY_DECOMP_DEF struct.
  //
  // -------------------

  char FILTERS_C11[NBAND_C11+1] = "vUBVRI" ;
  int  IFILT_v = 0 ;
  int  IFILT_U = 1 ;

  // define upper triangle of Table 14.2: roundoff -> not invertible
  /*    
  double COVAR_reduced_EWSi[NBAND_C11][NBAND_C11] = {
    +1.00, -0.12, -0.77, -0.91, -0.22,
    -0.12, +1.00, +0.57, -0.24, -0.89,
    -0.77, +0.57, +1.00, +0.53, -0.40,
    -0.91, -0.24, +0.53, +1.00, +0.49,
    -0.22, -0.89, -0.40, +0.49, +1.00
  } ;
  double COVAR_diag_EWSi[NBAND_C11] = { 0.06, 0.04, 0.05, 0.04, 0.08  };
  double COV_DIAG_FUDGE = 1.0E-5 ; 
  */

  
  // Fudged by Marriner to be invertible
  // Default for 'v' band is no correlation with other bands.
  double COVAR_reduced_EWSi[NBAND_C11][NBAND_C11] = 
    {
      {+1.000, 0.000,     0.000,     0.000,     0.000,     0.000},
      {0.000, +1.000000, -0.118516, -0.768635, -0.908202, -0.219447},
      {0.000, -0.118516, +1.000000, +0.570333, -0.238470, -0.888611},
      {0.000, -0.768635, +0.570333, +1.000000, +0.530320, -0.399538},
      {0.000, -0.908202, -0.238470, +0.530320, +1.000000, +0.490134},
      {0.000, -0.219447, -0.888611, -0.399538, +0.490134, +1.000000}
    } ;

  double COVAR_diag_EWSi[NBAND_C11] = 
    { 0.5900, 0.06001, 0.040034, 0.050014, 0.040017, 0.080007  };

  double COV_DIAG_FUDGE = 1.0E-9 ; // needed to be invertible
  double COV_SCALE = 1.3 ; // to get correct filter-COV from lambda model


  /*
  // define lower triangle of Table 14.2
  double COVAR_reduced_EWSiCa[NBAND_C11][NBAND_C11] = {
    +1.00, +0.63, -0.73, -0.97, -0.52,
    +0.63, +1.00, +0.07, -0.80, -0.99,
    -0.73, +0.07, +1.00, +0.53, -0.20,
    -0.97, -0.80, +0.53, +1.00, +0.72,
    -0.52, -0.99, -0.20, +0.72, +1.00
  } ;
  double COVAR_diag_EWSiCa[NBAND_C11] = { 0.03, 0.02, 0.03, 0.03, 0.06  };
  */

  double COVAR1[NBAND_C11* NBAND_C11] ;
  double COVAR2[NBAND_C11][NBAND_C11] ;

  double CC, COVred, tmp, covscale_v ;
  int i,j, N ;
  int MEMD = NBAND_C11 * NBAND_C11 * sizeof(double);
  //  char fnam[] = "init_genSmear_Chotard11" ;

  // ------------ BEGIN -----------

  GENSMEAR_C11.USE = 1;    GENSMEAR.NUSE++ ;
  GENSMEAR_C11.OPT_farUV = OPT_farUV;
  GENSMEAR_C11.DECOMP.COVMAT1D = (double*) malloc(MEMD);
  GENSMEAR_C11.DECOMP.MATSIZE  = NBAND_C11 ;

  printf("\t Initialize Chotard11/SNF model of %s correlations\n", 
	 FILTERS_C11 );
  printf("\t OPT_farUV = %d \n", OPT_farUV);

  // set scaling for off-diagonal elements for farUV 'v'
  if ( OPT_farUV == 1 ) 
    { covscale_v = 0.999 ; }
  else if ( OPT_farUV == 2 ) 
    { covscale_v = -0.999 ; } 
  else
    { covscale_v = 0.0 ; }


  // translate reduced covariance into covariances
  N = 0 ;
  printf("\n\t Reduced UBVRI covariances: \n" );
  for (i =0; i < NBAND_C11; i++){
    printf("\t  "); fflush(stdout);
    for (j = 0; j < NBAND_C11 ; j++){      

      COVred       = COVAR_reduced_EWSi[i][j] ;

      // check option for farUV correlation
      if ( i == IFILT_v && j != IFILT_v ) 
	{ COVred = covscale_v * COVAR_reduced_EWSi[IFILT_U][j]; }
      if ( i != IFILT_v && j == IFILT_v ) 
	{ COVred = covscale_v * COVAR_reduced_EWSi[i][IFILT_U]; }

      CC   = COVAR_diag_EWSi[i] * COVAR_diag_EWSi[j];
      if ( i == j ) { CC += COV_DIAG_FUDGE ; }

      COVAR2[i][j] = COVred * CC ; // xxx REMOVE
      GENSMEAR_C11.DECOMP.COVMAT1D[N] = COVred * CC * COV_SCALE;

      N++ ;  
      COVAR1[N-1] = COVAR2[i][j] * COV_SCALE ; // xxx REMOVE

      printf(" %7.4f ", COVred );
    }
    printf("\n"); fflush(stdout);
  }

  init_Cholesky(+1, &GENSMEAR_C11.DECOMP); 

  // print Cholesky matrix
  printf("\n\t Cholesky Decomp: \n" );
  for (i =0; i < NBAND_C11; i++){
    printf("\t  d%c(Ran) = ", FILTERS_C11[i] ); fflush(stdout);
    for (j = 0; j < NBAND_C11 ; j++){      
      //      tmp = GENSMEAR_C11.Cholesky[j][i] ; // xxx REMOVE
      tmp = GENSMEAR_C11.DECOMP.CHOLESKY2D[j][i] ;
      printf("+ %7.4f*R%d ", tmp, j );
    }    
    printf("\n"); fflush(stdout);
  }

  //  debugexit(fnam); // DDDDDDD

  init_genSmear_randoms(NBAND_C11,0);

  return;

} // end of   init_genSmear_Chotard11


// ****************************************
void get_genSmear_Chotard11(double Trest, int NLam, double *Lam, 
			    double *magSmear) {

  
  int    ilam, IFILT ;
  double lam, tmp, SCATTER_VALUES[NBAND_C11] ;

  double LAMCEN[NBAND_C11] = 
    { 2500.0, 3560.0, 4390.0, 5490.0, 6545.0, 8045.0 } ;

  char fnam[] = "get_genSmear_Chotard11" ;

  // ---------- BEGIN -------

  // Feb 17 2020: use new utility for correlated randoms
  getRan_GaussCorr(&GENSMEAR_C11.DECOMP, GENSMEAR.RANGauss_LIST, // (I)
		   SCATTER_VALUES );            // (O)

  // -------------
  for ( ilam=0; ilam < NLam; ilam++ ) {

    lam = Lam[ilam];  // exact lambda
    
    // interpolate SCATTER values vs. wavelength
    if( lam < LAMCEN[0] ) {
      tmp   = SCATTER_VALUES[0] ;        // extend blueward of UV
    }
    else if ( lam > LAMCEN[NBAND_C11-1] ) {      
      tmp   = SCATTER_VALUES[NBAND_C11-1];  // extend redward of I band
    }
    else {
      IFILT = INODE_LAMBDA(lam, NBAND_C11, LAMCEN);
      if ( IFILT < 0 ) {
	sprintf(c1err,"Could not find UBVRI band for lam=%7.1f", lam);
	sprintf(c2err,"ilam = %d", ilam);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      // interpolate with sin function
      tmp = interp_SINFUN( lam, &LAMCEN[IFILT], &SCATTER_VALUES[IFILT], fnam );
    }

    magSmear[ilam] = tmp ;
  }


  return ;

} // end of get_genSmear_Chotard11



// ***************************************
void init_genSmear_VCR(char *VCR_version, int index_SNmodel) {

  // Jan 2014: Velocity-Color-Relation (VCR) model
  //
  // Inputs:
  //   *VCR_version: specific version for the intrinsic-scatter [VCR]
  //                  model; e.g., 'VCR.MFK14'
  // index_SNmodel: lightcurve model index from snana, 
  //                indicating SALT2, MLCS2k2, snoopy, etc ...
  //                  (grep MODEL_  sntools.h)
  //

  char *path ;
  //  char fnam[] = "init_genSmear_VCR" ;

  // ---------- BEGIN -------------

  GENSMEAR_VCR.USE = 1;   GENSMEAR.NUSE++ ;

  printf("\t Init Velocity-Color-Relation model: %s\n", VCR_version );
  fflush(stdout);

  // ------------------------------------
  // get path and filenames to parse
  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
  path = GENSMEAR_VCR.MODELPATH ;
  sprintf(path,"%s/models/VCR/%s", PATH_SNDATA_ROOT, VCR_version ) ;  
  sprintf(GENSMEAR_VCR.INFO_FILE, "%s/VCR.INFO", path);
  sprintf(GENSMEAR_VCR.VSI_FILE,  "%s/v_Si.dat", path);
  printf("\t Read model info from :\n\t\t %s\n", path); 
  fflush(stdout);

  // ------- read v_Si distribution ------------
  read_VCR_VSI();

  // -------------------------------------
  // read VCR.INFO file 
  read_VCR_INFO();

  // sort VCR bands by wavelength
  sort_VCR_BANDS();

  // prepare covariance matrix
  prep_VCR_COVAR();

  // SALT2-specific preparation
  if ( index_SNmodel == MODEL_SALT2 ) { prep_VCR_forSALT2(); }
  
  // ----------------------------
  // store number of randoms to generate.

  int NRANGauss = GENSMEAR_VCR.NCOLOR + 1; // each color + sigma_Coh
  int NRANFlat  = 1 ;  // pick v_Si
  init_genSmear_randoms(NRANGauss,NRANFlat);

} // end of init_genSmear_VCR



// **********************************
void read_VCR_VSI(void) {

  // read distribution of silicon velocities, with units of 10^3 km/sec/
  // Load distribution into GENSMEAR_VCR structure.

  int MAXLEN, MEMD, i, NBVSI ;
  double VSI, PROB, SUMPROB, SUMVSI ;

  MAXLEN               = 100;
  MEMD                 = sizeof(double);
  GENSMEAR_VCR.VSI     = (double*) malloc( MAXLEN * MEMD );
  GENSMEAR_VCR.PROB    = (double*) malloc( MAXLEN * MEMD );
  GENSMEAR_VCR.SUMPROB = (double*) malloc( MAXLEN * MEMD );

  rd2columnFile(GENSMEAR_VCR.VSI_FILE, MAXLEN, 
		&GENSMEAR_VCR.NBIN_VSI,  // return number of bins read
		GENSMEAR_VCR.VSI,        // return v_Si axis
		GENSMEAR_VCR.PROB, 0 );     // return prob axis
  
  NBVSI = GENSMEAR_VCR.NBIN_VSI;

  // get running sum needed for random selection
  SUMPROB = SUMVSI = 0.0 ;
  for(i=0; i < NBVSI; i++ ) {
    VSI      = GENSMEAR_VCR.VSI[i] ;
    PROB     = GENSMEAR_VCR.PROB[i] ;
    GENSMEAR_VCR.SUMPROB[i] = SUMPROB ;  // previous SUMPROB
    SUMPROB += PROB ;
    SUMVSI  += (VSI*PROB);
  }

  // renormalize SUMPROBs so that max prob = 1
  SUMPROB = GENSMEAR_VCR.SUMPROB[NBVSI-1] ;
  for(i=0; i < NBVSI; i++ ) { GENSMEAR_VCR.SUMPROB[i] /= SUMPROB ; }
  GENSMEAR_VCR.VSI_MEAN = SUMVSI/SUMPROB ;

  printf("\t v_Si range: %.2f - %.2f  x10^3 km/sec \n",
	 GENSMEAR_VCR.VSI[0], GENSMEAR_VCR.VSI[NBVSI-1]);

  //  printf("\t <v_Si> : %.3f  x10^3 km/sec \n",  GENSMEAR_VCR.VSI_MEAN );

  fflush(stdout) ;

  /*
  for(i=0; i < NBVSI; i++ ) { 
    printf(" XXX bin %3d: vsi=%.3f  PROB=%f  SUMPROB=%f \n", i,
	   GENSMEAR_VCR.VSI[i],
	   GENSMEAR_VCR.PROB[i],
	   GENSMEAR_VCR.SUMPROB[i]     );  fflush(stdout);  }    
  */ 

} // end of read_VCR_VSI


// **********************************
void read_VCR_INFO(void) {

  // read model params from VCR.INFO file
  //
  FILE *fp ;
  int  ic, ic2, NC, irowmat, iband, IFILTDEF, j ;
  char c_get[60], band[2];
  double LAMCEN ;
  char *infoFile = GENSMEAR_VCR.INFO_FILE ;
  char fnam[]    = "read_VCR_INFO";

  // ------------- BEGIN ------------

  if ( (fp = fopen(infoFile,"rt") ) == NULL ) {
    sprintf(c1err,"Cannot open info file:");
    sprintf(c2err,"%s", infoFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  GENSMEAR_VCR.NCOLOR = NC = 0 ;
  GENSMEAR_VCR.SIGMACOH_MB = NULLDOUBLE ;

  for(ic=0; ic < MXCOLOR_VCR; ic++ ) {
    GENSMEAR_VCR.COLOR_STRING[ic][0] = 0 ;
    GENSMEAR_VCR.COLOR_SLOPE[ic]  = NULLDOUBLE ;
    GENSMEAR_VCR.COLOR_SIGMA[ic]  = NULLDOUBLE ;

    for(ic2=0; ic2 < MXCOLOR_VCR; ic2++ )
      { GENSMEAR_VCR.COLOR_CORMAT[ic][ic2] = NULLDOUBLE ; }
  }

  for(iband=0; iband < 2*MXCOLOR_VCR; iband++  )
    { GENSMEAR_VCR.LAMCEN_BAND[iband] = -999.0 ; }


  GENSMEAR_VCR.COLOR_SIGMA_SCALE = 1.0 ;
  GENSMEAR_VCR.SPECSHIFT_SCALE = 1.0 ;


  GENSMEAR_VCR.USE_ZVAR = 0 ;
  for(j=0; j < 4; j++ ) 
    { GENSMEAR_VCR.ZVARIATION_POLY_VSI[j] = 0.0 ; }

  irowmat = 0; // number of matrix rows read.

  while( (fscanf(fp, "%s", c_get)) != EOF ) {
    
    if ( strcmp(c_get,"NCOLOR:") == 0  ) { 
      readint(fp, 1, &NC); GENSMEAR_VCR.NCOLOR = NC;  
      if ( NC > MXCOLOR_VCR ) {
	sprintf(c1err,"NCOLOR = %d  exceeds MXCOLOR_VCR=%d",
		NC, MXCOLOR_VCR ) ;
	sprintf(c2err,"Check NCOLOR key and #define MXCOLOR_VCR");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
      }
    }

    if ( strcmp(c_get,"COLOR_DEFINE:") == 0  ) {
      for ( ic=0; ic < NC; ic++ ) 
	{ readchar(fp, GENSMEAR_VCR.COLOR_STRING[ic] ); }
    }

    if ( strcmp(c_get,"COLOR_SLOPE:") == 0  ) 
      { readdouble(fp, NC, GENSMEAR_VCR.COLOR_SLOPE ); }

    if ( strcmp(c_get,"COLOR_SIGMA:") == 0  ) 
      { readdouble(fp, NC, GENSMEAR_VCR.COLOR_SIGMA ); }

    if ( strcmp(c_get,"COLOR_SIGMA_SCALE:") == 0  ) 
      { readdouble(fp, 1, &GENSMEAR_VCR.COLOR_SIGMA_SCALE ); }

    if ( strcmp(c_get,"COLOR_CORMAT:") == 0  ) {
      readdouble(fp, NC, GENSMEAR_VCR.COLOR_CORMAT[irowmat] ); 
      irowmat++ ;
    }

    if ( strcmp(c_get,"SIGMACOH_MB:") == 0  ) 
      { readdouble(fp, 1, &GENSMEAR_VCR.SIGMACOH_MB); }


    if ( strcmp(c_get,"SPECSHIFT_SCALE:") == 0  ) 
      { readdouble(fp, 1, &GENSMEAR_VCR.SPECSHIFT_SCALE); }

    if ( strcmp(c_get,"LAMCEN_BAND:") == 0  ) {
      readchar(fp, band);
      readdouble(fp, 1, &LAMCEN );
      IFILTDEF = INTFILTER(band);
      GENSMEAR_VCR.LAMCEN[IFILTDEF] = LAMCEN ;
    }


    if ( strcmp(c_get,"ZVARIATION_POLY_VSI:") == 0  )   { 
      readdouble(fp, 4, GENSMEAR_VCR.ZVARIATION_POLY_VSI); 
      GENSMEAR_VCR.USE_ZVAR = 1 ;
    }

  } // end while


  fflush(stdout);

  // ------------------------------------------
  // check for sim-input override parameters
  int NVAL, ipar;
  
  if ( NSMEARPAR_OVERRIDE > 0 ) { printf("\n"); }
  for(ipar=0; ipar < NSMEARPAR_OVERRIDE; ipar++ ) {

    NVAL = exec_genSmear_override(ipar,"COLOR_SIGMA_SCALE", 
				  &GENSMEAR_VCR.COLOR_SIGMA_SCALE );

    NVAL = exec_genSmear_override(ipar,"ZVARIATION_POLY_VSI", 
				  GENSMEAR_VCR.ZVARIATION_POLY_VSI );
    if ( NVAL > 0 ) { GENSMEAR_VCR.USE_ZVAR = 1; }


  } // end ipar loop
  

  // -------------------------------------------------------
  // sanity checks; because sane people do insane things.

  if ( irowmat != NC ) {
    sprintf(c1err,"Read %d rows of corr. matrix; expected %d rows.",
	    irowmat, NC);
    sprintf(c2err,"Check COLOR_CORMAT keys in VCR.INFO file.") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  if ( GENSMEAR_VCR.NCOLOR == 0 ) 
    {  missingKey_ABORT("NCOLOR:", infoFile, fnam);  }

  if ( GENSMEAR_VCR.SIGMACOH_MB == NULLDOUBLE ) 
    {  missingKey_ABORT("SIGMACOH_MB:", infoFile, fnam);  }

  // -----
  if ( GENSMEAR_VCR.USE_ZVAR ) {
    double tmp ;
    printf("\n   Redshift dependence:\n");
    printf("\t v_Si -> v_Si ");
    for(j=0; j<4; j++ ) {
      tmp = GENSMEAR_VCR.ZVARIATION_POLY_VSI[j] ;
      if( tmp != 0.0 ) { printf(" + (%.2f * z^%d)", tmp, j); }
    }
    printf("\n");
    fflush(stdout);
  }

  // -------------------------------------------------------------
  // find mean color to use as offset so that 
  // mean color shift is always zero. Print summary.
  // Note that sigma is the scatter about the average linear trend,
  // while RMS is the scatter about 0; thus RMS > sigma.

  double OFFSET, RMS ;
  printf("\n\t ColorDef  dColor/dv_Si   sigma   rmsTot    OFFSET \n" );
  fflush(stdout);


  for(ic=0; ic < NC; ic++ ) {

    // apply user option to multiply all of the COLOR_SMEAR terms.
    GENSMEAR_VCR.COLOR_SIGMA[ic] *= GENSMEAR_VCR.COLOR_SIGMA_SCALE ;

    parse_VCR_colorString(ic) ;
    get_VCR_colorOffset(ic, &OFFSET, &RMS); // return OFFSET and RMS

    GENSMEAR_VCR.COLOR_OFFSET[ic] = -OFFSET ;

    printf("\t   %s-%s       %6.3f       %5.3f   %5.3f   %6.3f \n"
	   , GENSMEAR_VCR.COLOR_BAND[ic][0]
	   , GENSMEAR_VCR.COLOR_BAND[ic][1]
	   , GENSMEAR_VCR.COLOR_SLOPE[ic] 
	   , GENSMEAR_VCR.COLOR_SIGMA[ic] 
	   , RMS
	   , GENSMEAR_VCR.COLOR_OFFSET[ic] );
    fflush(stdout);
  }

  return ;

} // end of read_VCR_INFO


// **********************************************
void get_VCR_colorOffset(int ic, double *OFFSET, double *RMS ) {

  // for color index ic, compute OFFSET so that average
  // color is 0. Also return RMS.

  int i, NBVSI ;
  double VSI, PROB, color, slope, sig ;
  double sum_PxColor, sum_PROB, PxColor, colorDif, SQSUM;
  //  char fnam[] = "get_VCR_colorOffset" ;

  // ----------- BEGIN --------------

  slope = GENSMEAR_VCR.COLOR_SLOPE[ic] ;  // dColor/dv_Si
  NBVSI = GENSMEAR_VCR.NBIN_VSI ;
  sum_PxColor = sum_PROB = 0.0 ;
  for( i=0; i < NBVSI; i++ ) { 
    VSI           = GENSMEAR_VCR.VSI[i]  ; // note that VSI = |v_Si|
    PROB          = GENSMEAR_VCR.PROB[i] ;
    color         = slope*VSI  ;  // 
    PxColor       = PROB * color ;
    sum_PxColor  += PxColor ;
    sum_PROB     += PROB ; 
  }
  *OFFSET = sum_PxColor / sum_PROB ; // average color shift
  
  // now estimate RMS by brute force, and include COLOR_SIGMA.
  // RMS is just for screen print and is not used for anything.
  SQSUM = 0.0 ;
  for( i=0; i < NBVSI; i++ ) { 
    VSI         = GENSMEAR_VCR.VSI[i] ;
    PROB        = GENSMEAR_VCR.PROB[i] ;
    color       = slope*VSI ;
    colorDif    = ( color - *OFFSET ) ;
    SQSUM      += PROB * (colorDif * colorDif ) ;
  }
  sig = GENSMEAR_VCR.COLOR_SIGMA[ic];
  *RMS = sqrt(SQSUM/sum_PROB + sig*sig) ;

} // end of  get_VCR_colorOffset


// ****************************************
void sort_VCR_BANDS(void) {

  // use already filled GENSMEAR_VCR.LAMCEN[IFILTDEF] array
  // to define the list of bands, then sort them by lambda 
  // and fill ORDERED arrays.

  int ORDER =  1; // --> increasing order
  int NBAND, IFILTDEF, IFILTDEF_B, iband, isort, ic, NC ;
  int IFILTDEF_LIST[MXFILTINDX];
  int INDSORT[MXFILTINDX];
  double LAMCEN, LAMCEN_LIST[MXFILTINDX];
  char fnam[] = "sort_VCR_BANDS" ;

  // ------------ BEGIN -------------

  // loop over every possible filter index,
  // and pick out the ones with LAMCEN > 0
  NBAND = 0 ;
  for( IFILTDEF=0; IFILTDEF < MXFILTINDX; IFILTDEF++ ) {
    LAMCEN = GENSMEAR_VCR.LAMCEN[IFILTDEF] ;
    if ( LAMCEN > 0.0 ) {
      LAMCEN_LIST[NBAND]   = LAMCEN ;
      IFILTDEF_LIST[NBAND] = IFILTDEF ;
      NBAND++ ;
    }
  }
  GENSMEAR_VCR.NFILTDEF = NBAND ;
  
  // sort by increasing wavelength
  sortDouble(NBAND, LAMCEN_LIST, ORDER, INDSORT);
  NC = GENSMEAR_VCR.NCOLOR ;
  for(ic=0; ic < NC; ic++ ) { GENSMEAR_VCR.COLOR_IFILT[ic] = -9 ; }
      
  IFILTDEF_B = INTFILTER("B");

  // store in global ORDERED arrays.
  for(iband=0; iband < NBAND; iband++ ) {
    isort    = INDSORT[iband];
    LAMCEN   = LAMCEN_LIST[isort];
    IFILTDEF = IFILTDEF_LIST[isort];

    GENSMEAR_VCR.ORDERED_IFILTDEF[iband] = IFILTDEF ;
    GENSMEAR_VCR.ORDERED_LAMCEN[iband]   = LAMCEN ;

    if ( IFILTDEF == IFILTDEF_B ) { GENSMEAR_VCR.IFILT_B = iband; }

    for(ic=0; ic < NC; ic++ ) {
      if ( GENSMEAR_VCR.COLOR_IFILTDEF[ic][1] == IFILTDEF ) 
	{  GENSMEAR_VCR.COLOR_IFILT[ic] = iband ;  }
    }

    printf("\t\t VCR_LAMCEN(%c) = %7.0f \n", FILTERSTRING[IFILTDEF], LAMCEN);
    fflush(stdout);  
  }


  // make sure that each COLOR_IFILT is valid
  for(ic=0; ic < NC; ic++ ) {
    iband = GENSMEAR_VCR.COLOR_IFILT[ic] ;
    if ( iband < 0 || iband >= NBAND ) {
      sprintf(c1err,"Invalid COLOR_IFILT[%d] = %d", ic, iband);
      sprintf(c2err,"Check VCR model bands.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     	  
    }
  }


} // end of  sort_VCR_BANDS


// ****************************************
void  prep_VCR_COVAR(void) {

  int ic1, ic2, NC, N;
  double S1, S2, rho, COV ;
  //  char fnam[] = "prep_VCR_COVAR" ;

  // ------------- BEGIN ----------

  NC = GENSMEAR_VCR.NCOLOR ;
  N  = 0 ;

  for(ic1=0; ic1<NC; ic1++ ) {
    for(ic2=0; ic2<NC; ic2++ ) {
      rho = GENSMEAR_VCR.COLOR_CORMAT[ic1][ic2] ;
      S1  = GENSMEAR_VCR.COLOR_SIGMA[ic1] ;
      S2  = GENSMEAR_VCR.COLOR_SIGMA[ic2] ;
      COV = (S1 * S2 * rho) ;
      GENSMEAR_VCR.COLOR_COVAR2[ic1][ic2] = COV ;

      N++ ;
      GENSMEAR_VCR.COLOR_COVAR1[N-1] = COV ;

    }  // ic2
  }    // ic1


  gsl_matrix_view chk;  

  chk  = gsl_matrix_view_array ( GENSMEAR_VCR.COLOR_COVAR1, NC, NC); 
  gsl_linalg_cholesky_decomp ( &chk.matrix )  ;  

  for (ic1 =0; ic1 < NC ; ic1++){
    for (ic2 = 0; ic2 < NC ; ic2++) { 
      if ( ic2 >= ic1 ) 
	{ GENSMEAR_VCR.Cholesky[ic1][ic2] = 
	    gsl_matrix_get(&chk.matrix,ic1,ic2); }
      else
	{ GENSMEAR_VCR.Cholesky[ic1][ic2] = 0.0 ; }
    }
  }

} // end of   prep_VCR_COVAR

// ****************************************
void prep_VCR_forSALT2(void) {
  // Mar 2014  
  //  char fnam[] = "prep_VCR_forSALT2" ;
  // --------------- BEGIN -------------
  // ??
} // end of prep_VCR_forSALT2



// *****************************************
void  parse_VCR_colorString(int ic) {

  // Mar 2014
  // store filter index for each color,
  // and store unique list of all requried filters.

  char *cstr, c2[2] ;
  int len ;
  char fnam[] = "parse_VCR_colorString" ;

  // ------------ BEGIN -------------------

  // idiot checks; because sometimes smart people do dumb things.
  cstr = GENSMEAR_VCR.COLOR_STRING[ic] ;
  len = strlen(cstr);
  if ( len != 3 ) {
    sprintf(c1err,"Invalid VCR color definition: '%s' (len=%d !=3)", 
	    cstr, len );
    sprintf(c2err,"Color def must be of the form 'X-Y'");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  // make sure we have a dash as 2nd element.
  sprintf(c2,"%c", cstr[1] );
  if (strcmp(c2,"-") != 0 ) {
    sprintf(c1err,"Invalid VCR color definition: '%s'  (2nd char ne '-'", 
	    cstr);
    sprintf(c2err,"Color def must be of the form 'X-Y'");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  // strip bands from X-Y
  int  iband, IFILTDEF, USE;
  char bands[2][2], *BAND ;

  sprintf(bands[0],"%c", cstr[0] );
  sprintf(bands[1],"%c", cstr[2] );

  for(iband=0; iband <=1; iband++ ) {
    BAND = bands[iband] ;
    sprintf(GENSMEAR_VCR.COLOR_BAND[ic][iband],"%s", BAND );
    IFILTDEF = INTFILTER(BAND) ;
    GENSMEAR_VCR.COLOR_IFILTDEF[ic][iband] = IFILTDEF ;

    // abort if LAMCEN (i.e., node) is not defined for this band
    USE = ( GENSMEAR_VCR.LAMCEN[IFILTDEF] > 0.0 ) ;
    if ( USE == 0 ) { 
      sprintf(c1err,"LAMCEN_BAND undefined for band = '%s'", BAND);
      sprintf(c2err,"Check LAMCEN_BAND keys in VCR.INFO file");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
    } // USE

  } // end iband


} // end of   parse_VCR_colorString


// ****************************************
void get_genSmear_VCR(double Trest, int NLam, double *Lam, 
		      double *magSmear) {

  // 
  // Created Jan 2014 by R.Kessler
  // Model based on MFK14 paper (Mandel, Foley, Kirshner)
  // Implementation here is analogous to the Chotard 2011 model,
  // except that here there is another component related
  // to the Si velocity.
  //

  int ilam, NC, i, j, IFILTDEF, iband, NBAND ;

  double 
    VSI, VSI_zShift, SUMPROB, tmp, lam, ranB, sigB, magSmearB
    ,colorSmear[MXFILTINDX]
    ,LAMCEN[MXFILTINDX] 
    ,MAGSMEAR[MXFILTINDX]
    ;

  char fnam[] = "get_genSmear_VCR" ;
  
  // -------- BEGIN ------------

  for(iband=0; iband < GENSMEAR_VCR.NFILTDEF; iband++ ) 
    { MAGSMEAR[iband] = 9999. ; }
  
  if ( GENSMEAR_VCR.USE == 1 ) {  
    // select random v_Si
    int OPT_INTERP=1;
    SUMPROB = GENSMEAR.RANFlat_LIST[0] ; 
    VSI = interp_1DFUN(OPT_INTERP, SUMPROB, GENSMEAR_VCR.NBIN_VSI,
		       GENSMEAR_VCR.SUMPROB, GENSMEAR_VCR.VSI,  fnam) ;    
  }
  else if ( GENSMEAR_VCR.USE == 2 ) {
    // use externally set VSI (e.g., from SIMSED_fudge)
    VSI = GENSMEAR_VCR.GENRAN_VSI ;  
  }
  else {
    // abort
    VSI = -9.0 ;
    sprintf(c1err,"Invalid GENSMEAR_VCR.USE = %d", GENSMEAR_VCR.USE) ;
    sprintf(c2err,"Only 1 or 2 is valid.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }


  // check for redshift dependence 
  VSI_zShift = 0.0 ;
  if ( GENSMEAR_VCR.USE_ZVAR ) {
    for(j=0; j<4; j++ ) {
      tmp = GENSMEAR_VCR.ZVARIATION_POLY_VSI[j] ;
      VSI_zShift += (GENSMEAR.ZPOW[j] * tmp);
    }
  }
  VSI += VSI_zShift ;

  // load VSI to global
  GENSMEAR_VCR.GENRAN_VSI = VSI ;

  // -------------------------
  NC = GENSMEAR_VCR.NCOLOR ;

  // get coherent/luminosity scatter in B band;
  // note that randoms 0 - NC-1 are reserved for the Cholesky matrix below.
  ranB = GENSMEAR.RANGauss_LIST[NC] ;
  sigB = GENSMEAR_VCR.SIGMACOH_MB ;
  magSmearB = ranB * sigB ;
  MAGSMEAR[GENSMEAR_VCR.IFILT_B] = magSmearB ;

  //Matrix Multiply to get color-scatter values = ch^T normalvector

  for (i = 0 ; i < NC ; i++) {
    colorSmear[i] = 0.0 ;  
    for (j = 0 ; j < NC ; j++){
      //transpose cholesky matrix
      tmp = GENSMEAR_VCR.Cholesky[j][i] ;      
      colorSmear[i] += tmp * GENSMEAR.RANGauss_LIST[j] ;
    }
  }


  // get color shifts
  double slope, offset, colorShift, magShift ;

  for (i = 0 ; i < NC ; i++) {
    slope  = GENSMEAR_VCR.COLOR_SLOPE[i]; 
    offset = GENSMEAR_VCR.COLOR_OFFSET[i]; 

    colorShift = offset + (slope * VSI) + colorSmear[i] ;
    GENSMEAR_VCR.GENRAN_COLORSHIFT[i] = colorShift ; // un-scaled

    // get absolute and sparse filter index of band to smear
    IFILTDEF = GENSMEAR_VCR.COLOR_IFILTDEF[i][1] ; 
    iband    = GENSMEAR_VCR.COLOR_IFILT[i] ;

    colorShift *= GENSMEAR_VCR.SPECSHIFT_SCALE ;
    magShift = magSmearB - colorShift ;
    MAGSMEAR[iband] = magShift ;
    
  } // end i loop over colors
  

  NBAND = GENSMEAR_VCR.NFILTDEF ; // number of VCR model bands
  for (i = 0 ; i < NBAND ; i++) 
    { LAMCEN[i] = GENSMEAR_VCR.ORDERED_LAMCEN[i] ;  }


  // idiot check that each band has a valid MAGSMEAR
  for(iband=0; iband < NBAND; iband++ ) { 
    if ( MAGSMEAR[iband] > 9990.0 ) {
      IFILTDEF = GENSMEAR_VCR.ORDERED_IFILTDEF[iband] ;
      sprintf(c1err,"Undefined MAGSMEAR for band = '%c' (%d) ", 
	      FILTERSTRING[IFILTDEF], iband);
      sprintf(c2err,"Something is messed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );       
    }

  } // end of iband


  // -------------
  for ( ilam=0; ilam < NLam; ilam++ ) {

    lam = Lam[ilam];  // exact lambda
    
    // interpolate smear values vs. wavelength
    if( lam < LAMCEN[0] ) {
      tmp   = MAGSMEAR[0] ;        // extend blueward of UV
    }
    else if ( lam > LAMCEN[NBAND-1] ) {      
      tmp   = MAGSMEAR[NBAND-1] ;  // extend redward 
    }
    else {
      iband = INODE_LAMBDA(lam, NBAND, LAMCEN);
      if ( iband < 0 || iband >= NBAND ) {
	sprintf(c1err,"Could not find band for lam=%7.1f", lam);
	sprintf(c2err,"ilam = %d", ilam);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      // interpolate with sin function
      tmp = interp_SINFUN(lam, &LAMCEN[iband], &MAGSMEAR[iband], fnam);
    }
    

    magSmear[ilam] = tmp ;
  }


} // end of get_genSmear_VCR


// ***************************************
void init_genSmear_COH(char *stringArg) {

  // stringArg == "COH" 
  //     --> hard code coherent sigma = 0.13 mag
  //
  // stringArg = COH(0.11) 
  //     --> sigma = 0.11 mag
  //
  // stringArg = COH(0.55+.11) 
  //     --> sigma = 0.55 and sigma=0.11 (2 variations),
  //         and stored separate MAGSMEAR_COH values.
  //
  // The two scatter values can be output to SIMGEN_DUMP using
  // MAGSMEAR_COH and MAGSMEAR_COH2.
  //
  // Nov 30 2019: major refactor/update, and pass stringArg input
  //

  int  Nsigma = 0, i, ISBAD ;
  int  MEMC   = sizeof(char) * 20 ;
  double SIGMA;
  char stringLocal[60], stringSigma[60], *ptrSigma[2] ;
  char plus[] = "+" ;
  char fnam[] = "init_genSmear_COH" ;

  // -------------- BEGIN ----------------
  GENSMEAR_COH.USE = 1;    GENSMEAR.NUSE++ ;

  GENSMEAR_COH.MAGSIGMA[0] = -9.0;
  GENSMEAR_COH.MAGSIGMA[1] = -9.0;

  sprintf(stringLocal, "%s", stringArg);
  extractStringOpt(stringLocal, stringSigma);

  // only + is allowed between floats, so abort on invalid separators
  ISBAD=0;
  if ( strchr(stringSigma,':') != NULL )  { ISBAD = 1; }
  if ( strchr(stringSigma,'!') != NULL )  { ISBAD = 1; }
  if ( strchr(stringSigma,',') != NULL )  { ISBAD = 1; }
  if ( strchr(stringSigma,'%') != NULL )  { ISBAD = 1; }

  if ( ISBAD || strcmp(stringLocal,"COH") != 0 ) {
    sprintf(c1err,"Invalid GENMAG_SMEAR_MODELNAME: %s", stringArg);
    sprintf(c2err,"Must be COH or COH([sig]) or COH([sig0]+[sig1])");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }



  if ( strlen(stringSigma) == 0 ) {
    GENSMEAR_COH.MAGSIGMA[0] = 0.13 ; // hard-wired
    Nsigma = 1;
  }
  else {
    // examine sigma(s) in parentheses
    ptrSigma[0] = (char*) malloc(MEMC);
    ptrSigma[1] = (char*) malloc(MEMC);

    splitString(stringSigma, plus, fnam, 2, &Nsigma, ptrSigma);
    for(i=0; i < Nsigma; i++ ) 
      { sscanf(ptrSigma[i], "%le", &GENSMEAR_COH.MAGSIGMA[i])  ; }

    free(ptrSigma[0]); free(ptrSigma[1]);
  }

  GENSMEAR_COH.NSIGMA = Nsigma;
  for(i=0; i < Nsigma; i++ ) {
    SIGMA = GENSMEAR_COH.MAGSIGMA[i];
    printf("\t Coherent MAGSMEAR SIGMA = %.3f \n", SIGMA );
  }

  fflush(stdout);

  // set number of Gaussian randoms neede per SN
  init_genSmear_randoms(Nsigma,0);

  //  debugexit(fnam); // xxxx REMOVE

} // end init_genSmear_COH

void get_genSmear_COH(double Trest, int NLam, double *Lam, 
		      double *magSmear) {

  // Nov 20 2019: refactor and allow up to two scatter terms.
  int    Nsigma =   GENSMEAR_COH.NSIGMA;
  int    ilam, isig ; 
  double SIG, MAGSMEAR, SUM_MAGSMEAR = 0.0 ;
  //  char fnam[] = "get_genSmear_COH";

  // ------------ BEGIN -----------

  for(isig=0; isig < Nsigma; isig++ ) {
    SIG      = GENSMEAR_COH.MAGSIGMA[isig] ;
    MAGSMEAR = SIG*GENSMEAR.RANGauss_LIST[isig] ;
    GENSMEAR.MAGSMEAR_COH[isig] = MAGSMEAR;
    SUM_MAGSMEAR += MAGSMEAR ; // note linear sum, not in quadrature
  }

  /*
  printf(" xxx %s: SUM_MAGSMEAR = %f  (NLAM=%d) \n",
	 fnam, SUM_MAGSMEAR, NLam ); fflush(stdout);
  */

  for ( ilam=0; ilam < NLam; ilam++ )  
    { magSmear[ilam] = SUM_MAGSMEAR ; }

} // end of get_genSmear_COH


// ***************************************
void init_genSmear_biModalUV(void) {

  // Motivated by claim of NUV sub classes in Milne et al, 2015
  // July 21 2016: MAGU_SPLIT -> 0.55 (was 0.65) from D. Cinabro
  // July 29 2016: check for override params
  //
  // Also see "Search for NUV Sub-classes" in
  //   http://adsabs.harvard.edu/abs/2016arXiv161102984C
  //    (arXiv:1611.02984)
  //

  double LAMU_CEN = 3600.0 ;
  double LAMU_MIN = LAMU_CEN - 800.0 ;
  double LAMU_MAX = LAMU_CEN + 800.0 ;
  //  char fnam[] = "init_genSmear_biModalUV" ;
  // ------------ BEGIN ------------

  GENSMEAR_BIMODAL_UV.USE = 1;  GENSMEAR.NUSE++ ;

  GENSMEAR_BIMODAL_UV.LAMU_MIN   = LAMU_MIN ; // U-shift falls to
  GENSMEAR_BIMODAL_UV.LAMU_MAX   = LAMU_MAX ; // zero at MIN,MAX  
  GENSMEAR_BIMODAL_UV.MAGU_SPLIT = 0.55 ; // mag-split at LAMU_CEN
  GENSMEAR_BIMODAL_UV.MAGU_SIGMA = 0.04 ; // relative sigma on MAGU_SPLIT
  GENSMEAR_BIMODAL_UV.SIGMACOH   = 0.08 ; // coherent part

  int NGEN_RANGauss = 2; // sigmaCOH and sigmaUV
  int NGEN_RANFlat  = 1; // pick + or - mag shift

  init_genSmear_randoms(NGEN_RANGauss,NGEN_RANFlat);

  // July 29 2016: check for GENMAG_SMEARPAR_OVERRIDE
  int ipar ;
  if ( NSMEARPAR_OVERRIDE > 0 ) { printf("\n"); }
  for(ipar=0; ipar < NSMEARPAR_OVERRIDE; ipar++ ) {
    exec_genSmear_override(ipar,"BIMODAL_LAMU_MIN",    // smear->0 here
			   &GENSMEAR_BIMODAL_UV.LAMU_MIN);
    exec_genSmear_override(ipar,"BIMODAL_LAMU_MAX",    // smear->0 here
			   &GENSMEAR_BIMODAL_UV.LAMU_MAX);
    exec_genSmear_override(ipar,"BIMODAL_MAGU_SPLIT",   // UV peak sep, mag
			   &GENSMEAR_BIMODAL_UV.MAGU_SPLIT);
    exec_genSmear_override(ipar,"BIMODAL_MAGU_SIGMA",  // width of split peaks 
			   &GENSMEAR_BIMODAL_UV.MAGU_SIGMA);
    exec_genSmear_override(ipar,"BIMODAL_SIGMACOH",    // coherent smear
			   &GENSMEAR_BIMODAL_UV.SIGMACOH);
  } // end ipar


  // compute central wavelength from min,max
  GENSMEAR_BIMODAL_UV.LAMU_CEN = 
    0.5*(GENSMEAR_BIMODAL_UV.LAMU_MIN + GENSMEAR_BIMODAL_UV.LAMU_MAX);

  printf("\n BIMODAL_UV model of intrinsic scatter: \n" );

  printf("\t Coherent sigmaCOH: %.3f mag \n", 
	 GENSMEAR_BIMODAL_UV.SIGMACOH );

  printf("\t U_MAG split: %.3f with relative-sigma = %.3f \n",
	 GENSMEAR_BIMODAL_UV.MAGU_SPLIT,
	 GENSMEAR_BIMODAL_UV.MAGU_SIGMA );

  printf("\t U_MAG shift is max at <lamU> = %.0f A\n", LAMU_CEN);

  printf("\t U_MAG shift --> 0 at %0.f and %0.f A \n",
	 GENSMEAR_BIMODAL_UV.LAMU_MIN,
	 GENSMEAR_BIMODAL_UV.LAMU_MAX );

  printf("\n");
  fflush(stdout);

  return ;

} // end of  init_genSmear_biModalUV

void get_genSmear_biModalUV(double Trest, int NLam, double *Lam, 
			    double *magSmear) {

  // motivated by Milne et al, 2015.
  // This model is intended to approximately reproduce the 
  // Milne results, with a bi-model  UV scatter so that the 
  // analysis/fitting program can verify
  // sensitivity to such an effect.

  char fnam[] = "get_genSmear_biModalUV" ;
  int ilam;
  double SIGMA_COH, SHIFT_UV, SIGMA_UV;
  double LAM_NODE[2], MAG_NODE[2];
  double lam, LAMU_CEN, LAMU_MIN, LAMU_MAX, MAGU_SPLIT ;
  double ranFlat, ranGau0, ranGau1,  splitSign ; 

  int LDMP = 0; // ( fabs(Trest) < 0.1 );

  // -------------- BEGIN ----------------

  LAMU_CEN = GENSMEAR_BIMODAL_UV.LAMU_CEN ;
  LAMU_MIN = GENSMEAR_BIMODAL_UV.LAMU_MIN ;
  LAMU_MAX = GENSMEAR_BIMODAL_UV.LAMU_MAX ;
  MAGU_SPLIT = GENSMEAR_BIMODAL_UV.MAGU_SPLIT ;

  ranGau0 = GENSMEAR.RANGauss_LIST[0] ; 
  ranGau1 = GENSMEAR.RANGauss_LIST[1] ; 
  ranFlat = GENSMEAR.RANFlat_LIST[0] ; 

  SIGMA_COH = GENSMEAR_BIMODAL_UV.SIGMACOH ;
  SIGMA_UV  = GENSMEAR_BIMODAL_UV.MAGU_SIGMA ;

  if ( ranFlat < 0.5 ) { splitSign = -0.5; } else { splitSign = 0.5 ; }

  if ( LDMP ) {
    printf(" xxx --------------------------------------------------- \n" );
    printf(" xxx SIGMA(COH,UV) = %.3f,%.3f   splitSign=%.2f \n", 
	   SIGMA_COH, SIGMA_UV, splitSign);
    printf(" xxx ranGau0=%.3f  ranGau1=%.3f  ranFlat=%.3f \n",
	   ranGau0, ranGau1, ranFlat);
  }

  // return zero magSmear by default.
  for ( ilam=0; ilam < NLam; ilam++ )  { 

    lam = Lam[ilam];  // exact lambda

    SHIFT_UV = 0.0 ;
    if ( lam < LAMU_CEN && lam > LAMU_MIN ) {
      LAM_NODE[0] = LAMU_MIN ;  MAG_NODE[0] = 0.0 ;
      LAM_NODE[1] = LAMU_CEN ;  MAG_NODE[1] = MAGU_SPLIT * splitSign ;
      SHIFT_UV    = interp_SINFUN(lam, LAM_NODE, MAG_NODE, fnam );
    }
    else if ( lam > LAMU_CEN && lam < LAMU_MAX ) {
      LAM_NODE[0] = LAMU_CEN ;  MAG_NODE[0] = MAGU_SPLIT * splitSign ;;
      LAM_NODE[1] = LAMU_MAX ;  MAG_NODE[1] = 0.0 ;
      SHIFT_UV    = interp_SINFUN(lam, LAM_NODE, MAG_NODE, fnam );
    }

    magSmear[ilam] = 
      ranGau0*SIGMA_COH  +  // coherent term
      SHIFT_UV * ( 1.0 + ranGau1*SIGMA_UV) ; 

    // xxxxxxxxxxxxxxxxxxxxxx
    if ( LDMP && fabs(lam-3600.0) < 600.0  ) {
      printf(" xxx lam=%.0f    SHIFT_UV=%.3f  magSmear=%.3f \n",
	     lam, SHIFT_UV, magSmear[ilam] );
      fflush(stdout);
    }
    // xxxxxxxxxxxxxxxxxxxxxx

  } // end of ilam loop

} // end of get_genSmear_biModalUV


// ***************************************
void init_genSmear_private(void) {

  GENSMEAR_PRIVATE.USE = 1;  GENSMEAR.NUSE++ ;

  // set number of Gaussian randoms neede per SN
  init_genSmear_randoms(1,1);
}

void get_genSmear_private(double Trest, int NLam, double *Lam, 
			   double *magSmear) {

  int ilam; 

  // return zero magSmear by default.
  for ( ilam=0; ilam < NLam; ilam++ )  { magSmear[ilam] = 0.0 ; }
} // end of get_genSmear_private


// ***************************************
void init_genSmear_OIR(char *VERSION) {

  // Created Aug 30 2019 by R.Kessler and D.Jones
  // Optical+IR smear model based on CfA and CSP.
  // Very similar to Chotard11 method
  //
  // version: e.g., $PATH/OIR.v1
  //
  // Apr 6 2020 RK - pass VERSION argument 
  //

  char *PATH = GENSMEAR_OIR.MODELPATH ;
  char version[100];
  char FILTERS_OIR[NBAND_OIR+1] = "BgriYJH" ;
  
  double COV_DIAG_FUDGE = 1.0E-9 ; // needed to be invertible
  double COV_SCALE = 1.3 ; // to get correct filter-COV from lambda model
  double COVAR1[NBAND_OIR* NBAND_OIR] ;
  double COVAR2[NBAND_OIR][NBAND_OIR] ;

  double CC, COVred, tmp ;
  int i,j, N ;
  int LDMPCOV = 0; 
  gsl_matrix_view chk;  

  //  char fnam[] = "init_genSmear_OIR";

  // --------------- BEGIN ---------------

  GENSMEAR_OIR.USE = 1;    GENSMEAR.NUSE++ ;

  // ------------------------------------
  // get path and filenames to parse.
  // Note that VERSION may include full path.
  // Returned version is just the OIR.ABC part.

  extract_MODELNAME(VERSION,          // input may include full path
		    PATH, version);   // returned
  
  // if no path is given, use default path under SNDATA_ROOT
  if ( strlen(PATH) == 0 ) {
    sprintf(PATH,"%s/models/OIR/%s", PATH_SNDATA_ROOT, version ) ;  
  }

  // read in covmat from file - using utilities from VCR model
  printf("\t Init OIR smear model: %s\n", version );
  fflush(stdout);


  sprintf(GENSMEAR_OIR.INFO_FILE, "%s/OIR.INFO", PATH );
  printf("\t Read model info from :\n\t\t %s\n", GENSMEAR_OIR.INFO_FILE ); 
  fflush(stdout);

  // -------------------------------------
  // read OIR.INFO file 
  read_OIR_INFO();

  // sort OIR bands by wavelength
  sort_OIR_BANDS();

  // prepare covariance matrix
  prep_OIR_COVAR();

  printf("\t Initialize OIR model of %s correlations (A19) \n", 
	 FILTERS_OIR );

  // translate reduced covariance into covariances
  N = 0 ;

  if ( LDMPCOV ) 
    {  printf("\n\t xxx COVred(BgriYJH): \n" ); fflush(stdout); }

  for (i =0; i < NBAND_OIR; i++){

    if ( LDMPCOV )  { printf("\t xxx "); fflush(stdout); }
    for (j = 0; j < NBAND_OIR ; j++){      

      COVred       = GENSMEAR_OIR.COLOR_CORMAT[i][j] ;
      if ( LDMPCOV ) { printf(" %7.4f\n",COVred); }
      
      CC           = GENSMEAR_OIR.COLOR_SIGMA[i] * GENSMEAR_OIR.COLOR_SIGMA[j];
      
      if ( i == j ) { CC += COV_DIAG_FUDGE ; }
      COVAR2[i][j] = COVred * CC ;

      // fill 1D array for gsl argument below.
      N++ ;  COVAR1[N-1] = COVAR2[i][j] * COV_SCALE ;

    }
    if ( LDMPCOV ) { printf("\n"); fflush(stdout); fflush(stdout); }
  }

  chk  = gsl_matrix_view_array ( COVAR1, NBAND_OIR, NBAND_OIR); 
  gsl_linalg_cholesky_decomp ( &chk.matrix )  ;  

  for (i =0; i < NBAND_OIR ; i++){
    for (j = 0; j < NBAND_OIR ; j++) { 
      if ( j >= i ) 
	{ GENSMEAR_OIR.Cholesky[i][j] = gsl_matrix_get(&chk.matrix,i,j); }
      else
	{ GENSMEAR_OIR.Cholesky[i][j] = 0.0 ; }
    }
  }


  // print Cholesky matrix
  printf("\n\t Cholesky Decomp: \n" );
  for (i =0; i < NBAND_OIR; i++){
    printf("\t  d%c(Ran) = ", FILTERS_OIR[i] ); fflush(stdout);
    for (j = 0; j < NBAND_OIR ; j++){
      tmp = GENSMEAR_OIR.Cholesky[j][i] ;
      printf("+ %7.4f*R%d ", tmp, j );
    }    
    printf("\n"); fflush(stdout);
  }

  //  debugexit(fnam); // DDDDDDD

  // ----------------------------
  // store number of randoms to generate.


  int NRANGauss = NBAND_OIR + 1 ; //NCOLOR  + COH
  int NRANFlat  = 0;
  init_genSmear_randoms(NRANGauss,NRANFlat);

  return;

} // end init_genSmear_OIR



// =====================================================
void get_genSmear_OIR(double Trest, int NLam, double *Lam, 
		      double *magSmear) {

  // Created Aug 30 2019 by R.Kessler and D.Jones
  //

  int NBAND           = GENSMEAR_OIR.NBAND_OIRDEF ;
  double RANGauss_COH = GENSMEAR.RANGauss_LIST[NBAND];

  int    ilam, i, j, IFILT ;
  double lam, tmp, SCATTER_VALUES[NBAND_OIR], LAMCEN[NBAND_OIR] ;
  double SCATTER_COH;
  char fnam[] = "get_genSmear_OIR";

  // ---------------- BEGIN -----------------

  // RK - start with coherent scatter
  SCATTER_COH = GENSMEAR_OIR.SIGMACOH * RANGauss_COH ;

  //uBgriYJH filter wavelengths from 
  //   https://csp.obs.carnegiescience.edu/data/filters
  //from Swope site2 and RetroCam

  for (i = 0 ; i < NBAND ; i++) 
    { LAMCEN[i] = GENSMEAR_OIR.ORDERED_LAMCEN[i] ;  }


  //Matrix Multiply
  //scatter_values = ch^T normalvector
  for (i = 0 ; i < NBAND_OIR ; i++) {
    SCATTER_VALUES[i] = 0.0 ;  
    for (j = 0 ; j < NBAND_OIR ; j++){
      //transpose cholesky matrix
      tmp = GENSMEAR_OIR.Cholesky[j][i] ;      
      SCATTER_VALUES[i] += tmp * GENSMEAR.RANGauss_LIST[j] ;
    }
  }


  // -------------
  for ( ilam=0; ilam < NLam; ilam++ ) {

    lam = Lam[ilam];  // exact lambda
    
    // interpolate SCATTER values vs. wavelength
    if( lam < LAMCEN[0] ) {
      tmp   = SCATTER_VALUES[0] ;        // extend blueward of UV
    }
    else if ( lam > LAMCEN[NBAND_OIR-1] ) {      
      tmp   = SCATTER_VALUES[NBAND_OIR-1];  // extend redward of I band
    }
    else {
      IFILT = INODE_LAMBDA(lam, NBAND_OIR, LAMCEN);
      if ( IFILT < 0 ) {
	sprintf(c1err,"Could not find UBVRI band for lam=%7.1f", lam);
	sprintf(c2err,"ilam = %d", ilam);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      // interpolate with sin function
      //printf("xxxxxxxxxx %d\n",&SCATTER_VALUES[IFILT]);
      tmp = interp_SINFUN( lam, &LAMCEN[IFILT], &SCATTER_VALUES[IFILT], fnam );

    }

    magSmear[ilam] = tmp + SCATTER_COH ;

  }

  return;

} // end get_genSmear_OIR

// **********************************
void read_OIR_INFO(void) {

  // read model params from OIR.INFO file
  //
  FILE *fp ;
  int  ic, ic2, NC, irowmat, iband, IFILTDEF ;
  char c_get[60], band[2];
  double LAMCEN ;
  char *infoFile = GENSMEAR_OIR.INFO_FILE ;
  char fnam[]    = "read_OIR_INFO";

  // ------------- BEGIN ------------

  NC = 7;
  
  if ( (fp = fopen(infoFile,"rt") ) == NULL ) {
    sprintf(c1err,"Cannot open info file:");
    sprintf(c2err,"%s", infoFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }

  for(ic=0; ic < NBAND_OIR; ic++ ) {
    GENSMEAR_OIR.COLOR_SIGMA[ic]  = NULLDOUBLE ;

    for(ic2=0; ic2 < NBAND_OIR; ic2++ )
      { GENSMEAR_OIR.COLOR_CORMAT[ic][ic2] = NULLDOUBLE ; }
  }

  for(iband=0; iband < NBAND_OIR; iband++  )
    { GENSMEAR_OIR.LAMCEN_BAND[iband] = -999.0 ; }

  irowmat = 0; // number of matrix rows read.

  while( (fscanf(fp, "%s", c_get)) != EOF ) {
    if ( c_get[0] == '#' ) { continue; }

    if ( strcmp(c_get,"COLOR_SIGMA_SCALE:") == 0  ) 
      { readdouble(fp, 1, &GENSMEAR_OIR.COLOR_SIGMA_SCALE ); }

    if ( strcmp(c_get,"COLOR_CORMAT:") == 0  ) {
      readdouble(fp, NC, GENSMEAR_OIR.COLOR_CORMAT[irowmat] ); 
      irowmat++ ;
    }

    if ( strcmp(c_get,"LAMCEN_BAND:") == 0  ) {
      readchar(fp, band);
      readdouble(fp, 1, &LAMCEN );
      IFILTDEF = INTFILTER(band);
      GENSMEAR_OIR.LAMCEN[IFILTDEF] = LAMCEN ;
    }

    if ( strcmp(c_get,"COLOR_SIGMA:") == 0  ) 
      { readdouble(fp, NC, GENSMEAR_OIR.COLOR_SIGMA ); }

    if ( strcmp(c_get,"COLOR_SIGMA_SCALE:") == 0  ) 
      { readdouble(fp, 0, &GENSMEAR_OIR.COLOR_SIGMA_SCALE ); }
    
    if ( strcmp(c_get,"SIGMACOH_MB:") == 0  ||
	 strcmp(c_get,"SIGMACOH:"   ) == 0 ) 
      { readdouble(fp, 1, &GENSMEAR_OIR.SIGMACOH); }

    
  } // end while


  fflush(stdout);

  
  // -------------------------------------------------------
  // sanity checks; because sane people do insane things.

  if ( irowmat != NC ) {
    sprintf(c1err,"Read %d rows of corr. matrix; expected %d rows.",
	    irowmat, NC);
    sprintf(c2err,"Check COLOR_CORMAT keys in OIR.INFO file.") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );     
  }
  
  return ;

} // end of read_OIR_INFO

// ****************************************
void sort_OIR_BANDS(void) {

  // use already filled GENSMEAR_OIR.LAMCEN[IFILTDEF] array
  // to define the list of bands, then sort them by lambda 
  // and fill ORDERED arrays.

  int ORDER =  1; // --> increasing order
  int NBAND, IFILTDEF, IFILTDEF_B, iband, isort, NC ;
  int IFILTDEF_LIST[MXFILTINDX];
  int INDSORT[MXFILTINDX];
  double LAMCEN, LAMCEN_LIST[MXFILTINDX];
  //  char fnam[] = "sort_OIR_BANDS" ;

  // ------------ BEGIN -------------

  // loop over every possible filter index,
  // and pick out the ones with LAMCEN > 0
  NBAND = 0 ;
  for( IFILTDEF=0; IFILTDEF < MXFILTINDX; IFILTDEF++ ) {
    LAMCEN = GENSMEAR_OIR.LAMCEN[IFILTDEF] ;
    if ( LAMCEN > 0.0 ) {
      LAMCEN_LIST[NBAND]   = LAMCEN ;
      IFILTDEF_LIST[NBAND] = IFILTDEF ;
      NBAND++ ;
    }
  }
  GENSMEAR_OIR.NBAND_OIRDEF = NBAND ;
  
  // sort by increasing wavelength
  sortDouble(NBAND, LAMCEN_LIST, ORDER, INDSORT);
  NC = GENSMEAR_OIR.NCOLOR ;
      
  IFILTDEF_B = INTFILTER("B");

  // store in global ORDERED arrays.
  for(iband=0; iband < NBAND; iband++ ) {
    isort    = INDSORT[iband];
    LAMCEN   = LAMCEN_LIST[isort];
    IFILTDEF = IFILTDEF_LIST[isort];

    GENSMEAR_OIR.ORDERED_IFILTDEF[iband] = IFILTDEF ;
    GENSMEAR_OIR.ORDERED_LAMCEN[iband]   = LAMCEN ;

    if ( IFILTDEF == IFILTDEF_B ) { GENSMEAR_OIR.IFILT_B = iband; }

    printf("\t\t OIR_LAMCEN(%c) = %7.0f \n", FILTERSTRING[IFILTDEF], LAMCEN);
    fflush(stdout);  
  }


} // end of  sort_OIR_BANDS

// ****************************************
void  prep_OIR_COVAR(void) {

  int ic1, ic2, NC, N;
  double S1, S2, rho, COV ;
  //  char fnam[] = "prep_OIR_COVAR" ;

  // ------------- BEGIN ----------

  NC = NBAND_OIR ;
  N  = 0 ;
  for(ic1=0; ic1<NC; ic1++ ) {
    for(ic2=0; ic2<NC; ic2++ ) {
      rho = GENSMEAR_OIR.COLOR_CORMAT[ic1][ic2] ;
      S1  = GENSMEAR_OIR.COLOR_SIGMA[ic1] ;
      S2  = GENSMEAR_OIR.COLOR_SIGMA[ic2] ;
      COV = (S1 * S2 * rho) ;
      // xxx      printf("hiiiiiiiiiii8 %f\n",S1);
      N++ ;

    }  // ic2
  }    // ic1

} // end of   prep_OIR_COVAR



// ***************************************
void init_genSmear_COVSED(char *version, int OPTMASK) {

  // Created Oct 18 2019 by R.Kessler.
  // Covariance model in wavelength x phase bins.
  // Initial model is from Sugar/Pierre-Francoise Leget.
  //
  // Feb 17 2020: use CHOLESKY_DECOMP_DEF struct 


  double COV_DIAG_FUDGE = 1.0E-9 ; // needed to be invertible
  char *ptrPATH = GENSMEAR_COVSED.MODEL_PATH;
  char *ptrVERS = GENSMEAR_COVSED.VERSION;
  char *ptrFILE = GENSMEAR_COVSED.COVSED_FILE ;
  char fnam[] = "init_genSmear_COVSED" ;

  // --------------- BEGIN ---------------

  GENSMEAR_COVSED.USE  = 1;  GENSMEAR.NUSE++;
  GENSMEAR_COVSED.OPTMASK = OPTMASK;

  printf("\t %s \n", fnam);

  // separate input version into path and version
  ENVreplace(version, fnam, 1);
  extract_MODELNAME(version,
		    ptrPATH, ptrVERS);  // <== returned

  if ( strlen(ptrPATH) == 0 ) {
    // default path
    sprintf(ptrPATH,"%s/models/COVSED/%s", 
	    getenv("SNDATA_ROOT"), ptrVERS);
  }

  printf("   model path: %s\n",     ptrPATH);
  printf("   model version: %s \n", ptrVERS);

  parse_COVSED_INFO_FILE();

  // read contents of FITS file (cov matrix, wave bins ...)
  readFits_genSmear_COVSED(ptrFILE);

  // - - - - - - - - - - - 
  // prepare Cholesky decomp for correlated random in each WAVE bin

  int NBIN = GENSMEAR_COVSED.NBIN_WAVExEPOCH ;
  int i, jj;

  // fudge diagonal COV to ensure positive-definite
  for (i =0; i < NBIN ; i++){
    jj = NBIN*i + i;
    GENSMEAR_COVSED.DECOMP.COVMAT1D[jj] += COV_DIAG_FUDGE;
  }



  // - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Check option to monitor COV at two wavelengths
  // Make this check before any gsl calls modify COVMAT1D

  double *LAMPAIR   = GENSMEAR_COVSED.LAMPAIR_DEBUG;
  double  COVMAT[2][2] ; 
  double *WAVE_LIST = GENSMEAR_COVSED.WAVE;
  int    NBIN_WAVE  = GENSMEAR_COVSED.NBIN_WAVE;
  int i0, i1 ;
  if ( LAMPAIR[0] > 0.0 ) {
    i0 = quickBinSearch(LAMPAIR[0], NBIN_WAVE, WAVE_LIST, fnam);
    i1 = quickBinSearch(LAMPAIR[1], NBIN_WAVE, WAVE_LIST, fnam);
    LAMPAIR[0] = WAVE_LIST[i0] ; // snap to nearest bin
    LAMPAIR[1] = WAVE_LIST[i1] ; // snap to nearest bin
    COVMAT[0][0] = GENSMEAR_COVSED.DECOMP.COVMAT1D[NBIN*i0+i0] ;
    COVMAT[1][1] = GENSMEAR_COVSED.DECOMP.COVMAT1D[NBIN*i1+i1] ;
    COVMAT[0][1] = GENSMEAR_COVSED.DECOMP.COVMAT1D[NBIN*i0+i1] ;
    COVMAT[1][0] = GENSMEAR_COVSED.DECOMP.COVMAT1D[NBIN*i1+i0] ;
    init_genSmear_COVLAM_debug(LAMPAIR,COVMAT);
  }

  
  printf("\t Prepare Cholesky decomp. \n"); fflush(stdout);
  init_Cholesky(+1, &GENSMEAR_COVSED.DECOMP);

  // alocate array to use in gen_magSmear_COV
  GENSMEAR_COVSED.SCATTER_VALUES = (double*) malloc(NBIN*sizeof(double) );

  init_genSmear_randoms(NBIN,0);

  // free original COVMAT from memory
  free(GENSMEAR_COVSED.DECOMP.COVMAT1D);


  return;

} // end init_genSmear_COVSED


// =====================================================
void parse_COVSED_INFO_FILE(void) {

  char *ptrPATH    = GENSMEAR_COVSED.MODEL_PATH;
  char info_file[] = "COVSED.INFO";
  char *INFO_FILE  = (char*) malloc ( MXPATHLEN * sizeof(char) );
  char c_get[80], covsed_file[100];
  FILE *fp;
  char fnam[] = "parse_COVSED_INFO_FILE";

  // ----------- BEGIN ---------
  printf("   Read %s\n", info_file);
  sprintf(INFO_FILE, "%s/%s", ptrPATH, info_file);


  GENSMEAR_COVSED.COVSED_FILE[0]   = 0;
  GENSMEAR_COVSED.REBIN_LAM        = 1;
  GENSMEAR_COVSED.LAMPAIR_DEBUG[0] = -999.0;
  GENSMEAR_COVSED.LAMPAIR_DEBUG[1] = -999.0;

  if ( (fp = fopen(INFO_FILE, "rt")) == NULL ) { 
      sprintf(c1err,"Cannot open required info File " );
      sprintf(c2err,"%s", INFO_FILE);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;    
  }
  

  while ( (fscanf(fp, "%s", c_get)) != EOF ) {

    if ( strcmp(c_get,"COVSED_FILE:") == 0 )  { 
      readchar(fp, covsed_file);
      sprintf(GENSMEAR_COVSED.COVSED_FILE, "%s/%s", ptrPATH, covsed_file);
      printf("\t COVSED_FILE: %s \n", covsed_file);
    }
            
    if ( strcmp(c_get,"REBIN_LAM:") == 0 ) { 
      readint(fp, 1, &GENSMEAR_COVSED.REBIN_LAM); 
      printf("\t REBIN_LAM: %d \n", GENSMEAR_COVSED.REBIN_LAM);
    }

    if ( strcmp(c_get,"LAMPAIR_DEBUG:") == 0 ) { 
      double *ptrPair = GENSMEAR_COVSED.LAMPAIR_DEBUG ;
      readdouble(fp, 2, ptrPair);
      printf("\t LAMPAIR_DEBUG: %.1f and %.1f \n", 
	     ptrPair[0], ptrPair[1] );
    }

  }
  
  fflush(stdout);
  fclose(fp);
  return ;

} // end parse_COVSED_INFO_FILE


// =====================================================
void readFits_genSmear_COVSED(char *fileName) {

  // Created Oct 2019
  // for genSmear-COV model, read COV matrix from FITS fileName,
  // and store in structure.

  fitsfile *fp;
  long NCOLUMN; 
  int  icol, hdutype, lenstr, nmove=1, istat=0, NBIN ;
  int  NBIN_WAVE, NBIN_EPOCH, NBIN_COVMAT ;
  int  ICOL_WAVE, ICOL_EPOCH, ICOL_COVMAT;
  char keyname[40], comment[200], TTYPE[60], TFORM[60], tform[60] ;
  char fnam[] = "readFits_genSmear_COVSED" ;

  // -------------- BEGIN -------------

  printf("   %s: \n", fnam);

  // open file
  fits_open_file(&fp, fileName, READONLY, &istat );
  sprintf(c1err,"Open %s", fileName);
  snfitsio_errorCheck(c1err, istat);

  // move to table
  fits_movrel_hdu( fp, nmove, &hdutype, &istat );
  sprintf(c1err,"move to table" ) ;
  snfitsio_errorCheck(c1err, istat);

  sprintf(keyname, "TFIELDS" );
  fits_read_key(fp, TLONG, keyname,  &NCOLUMN, comment, &istat );
  sprintf(c1err, "read %s key", keyname);
  snfitsio_errorCheck(c1err, istat); 
  //  printf("\t Found %d columns \n", NCOLUMN); fflush(stdout);

  NBIN_WAVE = NBIN_EPOCH = NBIN_COVMAT = 0;
  ICOL_WAVE = ICOL_EPOCH = ICOL_COVMAT = -9;

  for ( icol = 1; icol <= NCOLUMN; icol++ ) {

    istat = 0 ;
    sprintf(keyname,"TTYPE%d", icol );
    fits_read_key(fp, TSTRING, keyname, TTYPE, comment, &istat );
    sprintf(c1err, "read %s key", keyname);
    snfitsio_errorCheck(c1err, istat);

    istat = 0 ;
    sprintf(keyname,"TFORM%d", icol );
    fits_read_key(fp, TSTRING, keyname, TFORM, comment, &istat );
    sprintf(c1err, "read %s key", keyname);
    snfitsio_errorCheck(c1err, istat);

    sprintf(tform, "%s", TFORM);
    lenstr = strlen(tform); tform[lenstr-1] = 0 ;
    sscanf(tform, "%d", &NBIN );

    if ( strstr(TTYPE,"WAVELEN") != NULL ) 
      { NBIN_WAVE = NBIN ;  ICOL_WAVE = icol;  }
    else if ( strstr(TTYPE,"EPOCH") != NULL ) 
      { NBIN_EPOCH = NBIN ; ICOL_EPOCH = icol;   }
    else if ( strstr(TTYPE,"MATRIX") != NULL )
      { NBIN_COVMAT = NBIN; ICOL_COVMAT = icol; }

    printf("\t Column %d:   TTYPE=%s   TFORM=%s  NBIN=%d\n",
	   icol, TTYPE, TFORM, NBIN ); fflush(stdout);
  }


  // - - - - - - - - - - - - - 
  // sanity checks
  sprintf(c1err,"Problem reading FITS file");
  if ( NBIN_WAVE <= 0 ) {
    sprintf(c2err,"NBIN_WAVE = %d", NBIN_WAVE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  if ( NBIN_COVMAT <= 0 ) {
    sprintf(c2err,"NBIN_COVMAT = %d", NBIN_COVMAT);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( NBIN_EPOCH <= 0 ) {
    NBIN_EPOCH = 1 ;  // maybe change later
  }

  int NBIN_TOT = NBIN_WAVE * NBIN_EPOCH;
  if ( NBIN_TOT*NBIN_TOT != NBIN_COVMAT ) {
    sprintf(c1err,"NBIN_TOT = %d(WAVE) x %d(EPOCH) = %d",
	    NBIN_WAVE, NBIN_EPOCH, NBIN_TOT);
    sprintf(c2err,"but expected NBIN_TOT = sqrt(NBIN_COVMAT) = %d",
	    (int)sqrt((double)NBIN_COVMAT) );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // malloc memory to store contents
  int MEMD1 = sizeof(double);

  GENSMEAR_COVSED.NBIN_WAVE   = NBIN_WAVE;
  GENSMEAR_COVSED.NBIN_EPOCH  = NBIN_EPOCH ;
  GENSMEAR_COVSED.NBIN_COVMAT = NBIN_COVMAT ;
  GENSMEAR_COVSED.NBIN_WAVExEPOCH = NBIN_WAVE * NBIN_EPOCH;

  GENSMEAR_COVSED.WAVE     = (double*) malloc ( MEMD1 * NBIN_WAVE    );
  GENSMEAR_COVSED.EPOCH    = (double*) malloc ( MEMD1 * NBIN_EPOCH   );
  GENSMEAR_COVSED.DECOMP.COVMAT1D = (double*) malloc ( MEMD1 * NBIN_COVMAT  );
  GENSMEAR_COVSED.DECOMP.MATSIZE  = NBIN_COVMAT;


  long FIRSTROW=1, FIRSTELEM=1, NREAD;
  int  anynul;

  // read wave column
  NREAD = NBIN_WAVE;  istat=0;
  fits_read_col_dbl(fp, ICOL_WAVE, FIRSTROW, FIRSTELEM, NREAD, NULL_1D,
		    GENSMEAR_COVSED.WAVE, &anynul, &istat);
  snfitsio_errorCheck("read WAVE bins", istat);
  
  // read covmat into 1D array
  NREAD = NBIN_COVMAT;  istat=0;
  fits_read_col_dbl(fp, ICOL_COVMAT, FIRSTROW, FIRSTELEM, NREAD, NULL_1D,
		    GENSMEAR_COVSED.DECOMP.COVMAT1D, &anynul, &istat);
  snfitsio_errorCheck("read COVMAT1D", istat);


  // - - - - - - - - - - - - - 
  // print summary

  printf("\t Read %d Wave bins: %.1f  to %.1f A \n",
	 NBIN_WAVE, GENSMEAR_COVSED.WAVE[0], GENSMEAR_COVSED.WAVE[NBIN_WAVE-1]);

  int iwave, jcov;
  double WAVE, COV, WAVE_LAST = -9.0 ;
  printf("\t Read %d COVMAT values: \n", NBIN_COVMAT); fflush(stdout);
  for(iwave=0; iwave < NBIN_WAVE; iwave++ ) {
    jcov = iwave*NBIN_WAVE + iwave;
    WAVE = GENSMEAR_COVSED.WAVE[iwave];
    COV  = GENSMEAR_COVSED.DECOMP.COVMAT1D[jcov];
    if ( (WAVE - WAVE_LAST) > 500.0 ) {
      printf("\t\t WAVE[%4d] = %6.0f : sqrt[ DIAG ] = %6.3f \n", 
	     iwave, WAVE, sqrt(COV) );      
      WAVE_LAST = WAVE ;
    }
  }

  fflush(stdout);
  fits_close_file(fp, &istat);

  return ;

} // end readFits_genSmear_COVSED


// =====================================================
void get_genSmear_COVSED(double Trest, int NWAVE, double *WAVE, 
			 double *magSmear) {

  // Created Oct 18 2019 by R.Kessler.
  //

  int    NBIN_WAVE      = GENSMEAR_COVSED.NBIN_WAVE;
  double WAVE_MIN       = GENSMEAR_COVSED.WAVE[0];
  double WAVE_MAX       = GENSMEAR_COVSED.WAVE[NBIN_WAVE-1];
  double *LAMPAIR_DEBUG = GENSMEAR_COVSED.LAMPAIR_DEBUG ;

  int    ISPEAK     = (fabs(Trest)<0.001);
  int    DEBUG      = ISPEAK && (LAMPAIR_DEBUG[0]>0.0) ;
  int    NDEBUG_STORE = 0 ;
  int    iwave;
  double tmp, tmpWave, magSmear_debug[2];
  char fnam[] = "get_genSmear_COVSED";

  // ---------------- BEGIN -----------------


  // Feb 17 2020: new utility to fetch correlated randoms
  getRan_GaussCorr(&GENSMEAR_COVSED.DECOMP, GENSMEAR.RANGauss_LIST, // (I)
		   GENSMEAR_COVSED.SCATTER_VALUES );        // (O)

  // -------------
  for ( iwave=0; iwave < NWAVE; iwave++ ) {

    tmpWave = WAVE[iwave];  // exact lambda
    
    // interpolate SCATTER values to current tmpWave

    if( tmpWave <= WAVE_MIN ) {
      tmp   = GENSMEAR_COVSED.SCATTER_VALUES[0] ; // extend blueward 
    }
    else if ( tmpWave >= WAVE_MAX ) {      
      tmp   = GENSMEAR_COVSED.SCATTER_VALUES[NBIN_WAVE-1];  // extend redward
    }
    else {
      // interpolate
      tmp = interp_1DFUN(1,tmpWave, NBIN_WAVE, GENSMEAR_COVSED.WAVE, 
			 GENSMEAR_COVSED.SCATTER_VALUES, fnam);
    }

    // -------

    if ( DEBUG ) {
      double wdif0 = fabs(tmpWave-LAMPAIR_DEBUG[0]);
      double wdif1 = fabs(tmpWave-LAMPAIR_DEBUG[1]);
      if ( wdif0 < 5.0 && NDEBUG_STORE==0 )
	{ magSmear_debug[0] = tmp;  NDEBUG_STORE++ ; }
      if ( wdif1 < 5.0 && NDEBUG_STORE==1 )
	{ magSmear_debug[1] = tmp;  NDEBUG_STORE++ ; }
    }
    
    magSmear[iwave] = tmp ;
  }


  /*
  printf(" xxx %s: DEBUG=%d Trest=%.3f  LAMPAIR=%.1f\n", 
	 fnam, DEBUG, Trest, LAMPAIR_DEBUG[0] ); fflush(stdout);
  */
 
  if ( NDEBUG_STORE == 2 ) 
    { update_genSmear_COVLAM_debug(magSmear_debug); }

  return;

} // end get_genSmear_COVSED


// *********************************************************
int INODE_LAMBDA(double LAM, int NNODE, double *LAM_NODES) {

  // April 12, 2012
  // Return INODE for which 
  // LAM_NODES[INODE] <= LAM <= LAM_NODES[INODE+1]
  // If not found then return -9

  int INODE, i ;
  //  char fnam[] = "INODE_LAMBDA" ;

  // -------------- BEGIN -----------

  INODE = -9;
  for (i=0; i < NNODE-1 ; i++ ) {
    if ( LAM >= LAM_NODES[i] && LAM <= LAM_NODES[i+1] ) 
      {  INODE = i ; }      
  }

  return INODE ;

} // end of INODE_LAMBDA


// ===============================================
void extraFilters_4genSmear(char *modelName,      // (I) name of scatter model
			    int *NFILT_extra,    // (O) Number of extra filt
			    int *IFILTOBS_extra) // (O) list of indices
{

  // Jan 23, 2014
  // For input scatter modelName, return list of 
  // extra filters needed to implement the model.
  // Called by snlc_sim to include these extra filter in the init.

  int NFILT  ;
  char s3[8];
  // ---------------- BEGIN --------------

  NFILT = 0 ; // default is no extra filters.


  snprintf(s3, 4, "%s ", modelName);
  if ( strcmp(s3,"VCR") == 0 ) {
    NFILT++ ; IFILTOBS_extra[NFILT-1] = INTFILTER("B");
    NFILT++ ; IFILTOBS_extra[NFILT-1] = INTFILTER("V");
    NFILT++ ; IFILTOBS_extra[NFILT-1] = INTFILTER("R");
    NFILT++ ; IFILTOBS_extra[NFILT-1] = INTFILTER("I");
  }



  *NFILT_extra = 0 ; // NFILT ; // load output arg
  
  //  printf(" xxxx s3 = '%s'  NFILT = %d \n",s3,NFILT); fflush(stdout);

} // end of 


// =====================
int  nval_genSmear_override(char *inputKey, char *parName) {

  // Mar 22 2014
  // Utility to parse override key and return number of values
  // associated with this parameter.
  //
  // Example for sim-input file :
  // GENMAG_SMEARPAR_OVERRIDE:  MYPAR      443
  // GENMAG_SMEARPAR_OVERRIDE:  MYLIST[4]  3 4 8 9
  //
  // Results in calling this function twice,
  //   nval = nval_genSmear_override("MYPAR",  *parName) ; 
  //   nval = nval_genSmear_override("MYLIST", *parName) ; 
  //
  // The 1st call returns nval=1 and parName = "MYPAR"
  // The 2nd call returns nval=4 and parName = "MYLIST"
  //
  
  int nval, br0, br1, LENKEY, i ;
  char c1[4];
  char cbr0[] = "[" ;
  char cbr1[] = "]" ;
  //  char fnam[] = "nval_genSmear_override" ;

  // ---------------- BEGIN --------------------

  nval = 1; // defalt return arg.
  parName[0] = 0 ;

  // check for optional brackets containing nval.
  
  LENKEY = strlen(inputKey);
  br0 = strcspn(inputKey, cbr0 );
  br1 = strcspn(inputKey, cbr1 );

  if ( br1 > br0 ) { 
    sscanf(&inputKey[br0+1], "%c", c1);
    sscanf(c1, "%d", &nval);
  }

  // strip off [] to get parName
  for(i=0; i < br0; i++)  { sprintf(&parName[i], "%c", inputKey[i]); }

  /*
  printf(" xxx -------------------------------- \n") ;
  printf(" xxx pass key = '%s' \n", inputKey); fflush(stdout);
  printf(" xxx br0,1 = %d,%d   LEN=%d \n", br0, br1, LENKEY );
  printf(" xxx nval = %d   parName = '%s' \n", nval, parName);
  fflush(stdout);
  debugexit(fnam); // xxxxx
  */


  return nval ;

} // end of nval_genSmear_override


// =============================
void  store_genSmear_override(char *parName, int NVAL, double *tmpList) {

  // Mar 22 2014
  // Call nval_genSmear_override to get parName and NVAL.
  // Then read list of params externally. Then call this function
  // to store parameter values.

  int N, i ;
  char cval[20];
  char fnam[] = "store_genSmear_override" ;

  // ------------- BEGIN ----------

  N = NSMEARPAR_OVERRIDE ;
  if ( N >= MXSMEARPAR_OVERRIDE || N < 0 ) {
    sprintf(c1err,"Invalid NSMEAR_OVERRIDE=%d (MXSMEAR_OVERRIDE=%d)", 
	    N, MXSMEARPAR_OVERRIDE );
    sprintf(c2err,"with parName= '%s' and NVAL=%d", parName, NVAL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(GENMAG_SMEARPAR_OVERRIDE[N].NAME,"%s", parName);
  GENMAG_SMEARPAR_OVERRIDE[N].NVAL = NVAL ;

  printf("\t %s : store %2d values for '%s' \n", fnam, NVAL, parName);
  fflush(stdout);
  for(i=0 ; i < NVAL; i++ ) 
    { GENMAG_SMEARPAR_OVERRIDE[N].VALUE[i] = tmpList[i] ; }
  NSMEARPAR_OVERRIDE++ ;

  return;

} // end of store_genSmear_override

int  exec_genSmear_override(int IPAR, char *PARNAME, double *VAL) {

  // execute override and load value(s) into pointer *VAL.
  // Function returns NVAL override values.

  char   *tmpName, parString[60];
  double *newVal, oldVal[MXSMEARPAR_OVERRIDE] ;
  int    NVAL, ival ;
  //  char fnam[] = "exec_genSmear_override" ;

  // ------------- BEGIN --------------

  tmpName = GENMAG_SMEARPAR_OVERRIDE[IPAR].NAME ;
  newVal  = GENMAG_SMEARPAR_OVERRIDE[IPAR].VALUE ;
  NVAL    = GENMAG_SMEARPAR_OVERRIDE[IPAR].NVAL ;

  if ( strcmp(tmpName,PARNAME) != 0  ) { return 0 ; }

  for(ival=0; ival < NVAL; ival++ ) {
    oldVal[ival] = VAL[ival] ;
    VAL[ival] = newVal[ival] ;

    if ( NVAL == 1 ) 
      { sprintf(parString,"%s", PARNAME) ; }
    else
      { sprintf(parString,"%s[%d]", PARNAME, ival) ; }

    printf("\t Override %s:  %.3f -> %.3f \n", 
	   parString, oldVal[ival], newVal[ival] );  
    fflush(stdout);
  }
  
  return NVAL ;

} // end exec_genSmear_override



// ==================================================
void  init_genSmear_phaseCor(double magSmear, double expTau) {

  // Created Feb 11 2020 by R. Kessler
  // Initialize phase-correlation model where 
  //  + magSmear = Gaussiam smear (few percent)
  //  + reduced correlation between any two epochs is 
  //    exp(-phaseDif/exptau)
  //
  // This smear-model should only be noticable for very
  // high-SNR epochs where SNR ~ few x 10.
  // Motivation is low FITPROB excess (compared to sims)
  // hard code GRID binning; maybe later add this to sim-input

  int NBIN_PHASECOR    =   60 ;
  double TMIN_PHASECOR = -18.0 ;
  double TBIN_PHASECOR =   1.0 ;

  int NBIN, N, MEMD, i, j ;
  double phase;
  char fnam[] = "init_genSmear_phaseCor" ;

  // ------------- BEGIN --------------

  if ( magSmear == 0.0 ) { return; } 

  printf("\n %s: magSmear=%.3f, rho=exp(-Tdif/%.1f) \n",
	 fnam, magSmear, expTau);

  init_genSmear_randoms(2, 0);
    
  fflush(stdout);

  
  GENSMEAR_PHASECOR.USE = 1;
  GENSMEAR_PHASECOR.INPUT_MAGSMEAR = magSmear ;  
  GENSMEAR_PHASECOR.INPUT_EXPTAU   = expTau;
  GENSMEAR_PHASECOR.CID_LAST       = -9 ;

  GENSMEAR_PHASECOR.NCHECK       = 0 ;
  GENSMEAR_PHASECOR.NSUM         = 0  ;
  GENSMEAR_PHASECOR.sumCHECK     = 0.0 ;
  GENSMEAR_PHASECOR.SUMCHECK[0]  = 0.0 ;
  GENSMEAR_PHASECOR.SUMCHECK[1]  = 0.0 ;
  GENSMEAR_PHASECOR.SUMCHECK[2]  = 0.0 ;
  GENSMEAR_PHASECOR.SUMCHECK[3]  = 0.0 ;
  
  NBIN = NBIN_PHASECOR ;
  MEMD = NBIN * sizeof(double) ;
  GENSMEAR_PHASECOR.GRID_PHASE    = (double*)malloc(MEMD);
  GENSMEAR_PHASECOR.GRID_MAGSMEAR = (double*)malloc(MEMD);
  GENSMEAR_PHASECOR.RANGauss_LIST = (double*)malloc(MEMD);

  for(i=0; i < NBIN; i++ ) {
    phase = TMIN_PHASECOR + TBIN_PHASECOR*(double)i ;
    GENSMEAR_PHASECOR.GRID_PHASE[i] = phase ;
  }

  GENSMEAR_PHASECOR.NBIN      = NBIN ;
  GENSMEAR_PHASECOR.RANGE[0]  = GENSMEAR_PHASECOR.GRID_PHASE[0] ;
  GENSMEAR_PHASECOR.RANGE[1]  = GENSMEAR_PHASECOR.GRID_PHASE[NBIN-1] ;
  GENSMEAR_PHASECOR.BINSIZE   = TBIN_PHASECOR;

  MEMD = NBIN * NBIN * sizeof(double);
  GENSMEAR_PHASECOR.DECOMP.COVMAT1D = (double *) malloc(MEMD);
  GENSMEAR_PHASECOR.DECOMP.MATSIZE  = NBIN ;

  // construct COV matrix vs. phase
  double Ti, Tj, TDIF, COV, RHO ;
  N=0;
  for(i=0; i < NBIN; i++ ) {
    for(j=0; j < NBIN; j++ ) {
      Ti   = GENSMEAR_PHASECOR.GRID_PHASE[i];
      Tj   = GENSMEAR_PHASECOR.GRID_PHASE[j];
      TDIF = fabs(Ti-Tj);
      RHO  = exp(-TDIF/expTau);
      COV  = magSmear * magSmear * RHO ;
      GENSMEAR_PHASECOR.DECOMP.COVMAT1D[N] = COV;
      N++ ;
    }
  }

  // printf("\t Prepare PhaseCor Cholesky decomp. \n"); fflush(stdout);
  init_Cholesky(+1, &GENSMEAR_PHASECOR.DECOMP);

  //  debugexit(fnam);

  return ;

} // end init_genSmear_phaseCor

// ----------------------------
void  get_genSmear_phaseCor(int CID, double phase, double *magSmear ) {

  int CID_LAST = GENSMEAR_PHASECOR.CID_LAST ;
  int NBIN     = GENSMEAR_PHASECOR.NBIN ;
  double magSmear_local = 0.0 ;
  char fnam[] = "get_genSmear_phaseCor" ;
  int  LDMP = 0 ;

  // ---------- BEGIN ---------

  *magSmear = 0.0 ;
  if ( !GENSMEAR_PHASECOR.USE ) { return; }

  // for each new event, prepare new list of correlated randoms
  if ( CID != CID_LAST ) {

    if ( LDMP ) 
      { printf(" xxx ----------- CID = %d -------------- \n", CID); }

    getRan_GaussCorr(&GENSMEAR_PHASECOR.DECOMP, 
		     GENSMEAR_PHASECOR.RANGauss_LIST,
		     GENSMEAR_PHASECOR.GRID_MAGSMEAR);

    check_genSmear_phaseCor();    
	
  } // end new CID

  GENSMEAR_PHASECOR.CID_LAST = CID;

  // interpolate magSmear at phase

  if ( phase <= GENSMEAR_PHASECOR.RANGE[0] ) {
    magSmear_local = GENSMEAR_PHASECOR.GRID_MAGSMEAR[0] ;
  }
  else if ( phase >= GENSMEAR_PHASECOR.RANGE[1] ) {
    magSmear_local = GENSMEAR_PHASECOR.GRID_MAGSMEAR[NBIN-1] ;
  }
  else {
    magSmear_local = interp_1DFUN(1, phase, NBIN, 
				  GENSMEAR_PHASECOR.GRID_PHASE, 
				  GENSMEAR_PHASECOR.GRID_MAGSMEAR, fnam);
  }
  
  *magSmear = magSmear_local ;
  return ;

} // end get_genSmear_phaseCor


// - - - 
void  check_genSmear_phaseCor(void) {

  // Feb 2020
  // Utility to verify the reduced covariance form the 
  // ADDPHASECOR model.

  int    IBIN_CHECK = -10 ;
  int    NBIN       = GENSMEAR_PHASECOR.NBIN ;
  double EXPTAU     = GENSMEAR_PHASECOR.INPUT_EXPTAU;
  double DT_CHECK[10], T0, m0, m1, RHO_PREDICT, RHO_CHECK;
  char fnam[] = "check_genSmear_phaseCor" ;

  // ------------- BEGIN ------------

  if ( IBIN_CHECK < 0 ) { return; }

  GENSMEAR_PHASECOR.NCHECK = 3; // dT=0, 1, max
  GENSMEAR_PHASECOR.NSUM++ ;

  DT_CHECK[0] = 0.0;
  DT_CHECK[1] = GENSMEAR_PHASECOR.BINSIZE ;
  DT_CHECK[2] = GENSMEAR_PHASECOR.RANGE[1] - GENSMEAR_PHASECOR.RANGE[0];

  // check neighboring bins
  m0 = GENSMEAR_PHASECOR.GRID_MAGSMEAR[IBIN_CHECK];
  GENSMEAR_PHASECOR.SUMCHECK[0] += (m0*m0);
  GENSMEAR_PHASECOR.sumCHECK    += m0 ;

  // check neighboring bins
  m0 = GENSMEAR_PHASECOR.GRID_MAGSMEAR[IBIN_CHECK];
  m1 = GENSMEAR_PHASECOR.GRID_MAGSMEAR[IBIN_CHECK+1];
  GENSMEAR_PHASECOR.SUMCHECK[1] += (m0*m1);
  
  // check first/last bin with min correlation
  m0 = GENSMEAR_PHASECOR.GRID_MAGSMEAR[0];
  m1 = GENSMEAR_PHASECOR.GRID_MAGSMEAR[NBIN-1] ;
  GENSMEAR_PHASECOR.SUMCHECK[2] += (m0*m1);

  int icheck, NSUM;
  double XN, SUMCHECK, sumCHECK, ARG, DT, SIG, SQSIG, RMS ;
  NSUM = GENSMEAR_PHASECOR.NSUM; 
  

  if ( NSUM > 0 && (NSUM % 100) == 0 ) {

    printf(" xxx --------------------------------- \n");
    printf(" xxx %s at N = %d: \n", fnam, NSUM );
    
    XN      = (double)NSUM ;
    SIG     = GENSMEAR_PHASECOR.INPUT_MAGSMEAR ;
    SQSIG   = SIG*SIG;
    T0      = GENSMEAR_PHASECOR.GRID_PHASE[IBIN_CHECK];

    sumCHECK    = GENSMEAR_PHASECOR.sumCHECK ;    // sum(m0)
    SUMCHECK    = GENSMEAR_PHASECOR.SUMCHECK[0] ; // sum(m0^2)
    RMS = STD_from_SUMS(NSUM, sumCHECK, SUMCHECK);

    for(icheck = 0; icheck < GENSMEAR_PHASECOR.NCHECK; icheck++ ) {
      DT  = DT_CHECK[icheck];       ARG = DT/EXPTAU;
      SUMCHECK    = GENSMEAR_PHASECOR.SUMCHECK[icheck] ;
      RHO_PREDICT = exp(-ARG);
      RHO_CHECK   = SUMCHECK / (XN*SQSIG);
      printf(" xxx RHO(dT=%4.1f) = %6.3f (expect RHO=%6.3f) ",
	     DT, RHO_CHECK, RHO_PREDICT);

      if ( icheck == 0 ) { printf("  RMS=%6.4f", RMS ); }
      printf("\n");

      fflush(stdout);
    }
    
    //    printf(" xxx NEVT=%5d: phaseRHO(dT=%.0f)=%.3f, 
  }

} // end check_genSmear_phaseCor
