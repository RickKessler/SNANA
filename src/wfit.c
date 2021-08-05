#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "fitsio.h"
#include "longnam.h"

#include "wfit.h"
#include "sntools.h"
#include "sntools_cosmology.h"
#include "sntools_output.h"

/***************************************************************************** 
   WFIT: Read in Hubble diagram and produce constraints in 
   Omega_m/w plane assuming flatness.

   - Read in redshifts and distances from fitres file.

   - Calculate chi2 for observed data and errors over a grid of 
      values in omega_m, w0 and H0.

   - Convert to likelihood.

   - Marginalize over H0.  

   - Create a probability distribution and optionally fold 
   in an Omega prior (default: 0.30 +/- 0.04 as in Tegmark '04)
   or BAO constraints (Eisenstein et al, '06)

   - Normalize and write out the distributions.
   
   - Marginalize over omega_m to get the uncertainty in w.

   Gajus Miknaitis, Fermilab

   ------------------
    CHANGES/HISTORY
   -----------------
   2006-07-18, GM: 

   * Added BAO constraints.  To avoid doing unnecessary integration,
    this involved a rewrite of the luminosity distance calculation as well, 
    "lumdist" has been replaced by "codist".  See comments in codist.c.

   * Also rewrote numerical integrator ("simpint") to allow for the passing 
     of additional parameters via a structure.  


   June 2, 2007 RSK: combine wfit.c & codist.c & simpint.c into one file.

   Nov 27, 2007 GM: fix fitter to run on large data sets.
                    Do analytic marginalization for H0.
     
  May 8, 2008: fix computation of chi2, and write chi2 & dof to .cospar file.

  May 21, 2008 RK - read cid instead of dummystr and write CID in
                    extra column in .resid file.

  May 29, 2008 RK - add -cmb option for CMB prior.
                    See "Rcmb" and zcmb variables.

             Replace codist() with Hainv_integral (sntools.c)
             for cmb integral, but NOT for other integrals.
             Hainv_integral is much faster for large z,
             but codist is faster for smaller z ?
             For 265 SN, wfit takes 4 minutes using codist for the
             CMB integral, and 1 minute using Hainv_integral;
             the w-results agree to within 0.001.

             Finally, note that wfit is linked to sntools.c
             and "sntool.h" is included.  There are only
             a few sntools function being used, so it would be
             easy to remove the sntools dependence if needed.

  June 10, 2008 RK:  write grid-chi2 to text file: char gridchi2file

  June 13, 2008 RK:  define mu_siginv so that mu_sig retains its meaning
                     Also define mu_sqsiginv to speed up loops.

                     New function "read_fitres" read old-style (LEGACY)
                     format or new self-documented format.

 Sep 11, 2008: write chi2tot with 14.6f (instead of 10.2f)
               to get better precision for contour-calculation.

 Sep 25, 2008: compute SALT2-style sigmu_int and add to cospar output.

 Nov 9, 2008: use Komatsu RWMAP values (1.710 +- 0.019)
              instead of values from Kowalski 2008.

 Jan 14, 2009 (RK): write SN-only contour to separate chi2grid file.
                    See chi2gridfile_SN

 Jan 16, 2009: print both mean and mode for w.
               Init w0_prob[i] = 0.0 in marginalization loop;
               lucky that this fix did not change results.

 Jan 23, 2009: 
    - replace fprintf with printf

    - make new function get_chi2wOM so that refined chi2grid
      can be used to get w,OM at min chi2. "-minchi2" option
      reports values at minchi2 instead of marginalized values.

 Jan 30, 2009: usemarg=0 by default => report values at chi2 min

 Apr 23, 2009 RSK - accept MU or DLMAG, MUERR or DLMAGERR (for v8_02)
                  - fix dumb bug in get_minwOM.  minimized values
                    were off more as initial grid is courser.

 May 1, 2009: 
   -  change default omSteps from 51 to 81.
   - Add -marg option to marginalize.
   - default result is UNBLINDED !  Use -blind to hide result.


 May 20, 2009: read distance covariances with option -covar <filename>
               apply to chi2 in get_chi2wOM().

 Jun 23, 2009: fix dumb bug writing to residFile: add mu_sig and sqsnrms
               in quadrature.

               Fix read_fitres to work if CID is specified in fitres file.

 Aug 20, 2009: when reading fitres file, allow either Z[ERR] or ZPHOT[ERR]
               to be used as the redshift.

 Apr 29, 2011: allow IDSURVEY or IDTEL key

 Aug 4, 2011 - fix declarations that gave compile warnings:
               void * -> Cosparam * here and in wfit.h.


 Sep 12, 2011: 
      MXSN -> 100,000 (from 10,000)
      New selection options: 
         -GENTYPE <gentype>       
         -SNTYPE <sntype>

 Mar 1, 2012:
      Add -Rcmb and -sigma_Rcmb flags to allow user-selected priors. 
      To use this option, -cmb flag must also be set. [JLM]

 Jun 13, 2012:
      Add -refit option. If this flag is included, wfit will fit cosmo twice,
      first time as is, second time using intrinsic scatter found 
      in first fit. [JLM]

 Jun 27, 2012:
      Fixed dumb bug - only refit cosmo data if intrinsic scatter > snrms. [JLM]
      Made getname subroutine to remove duplications in code.
      Since intrinsic scatter includes any input snrns, set refit snrms to int. scat.

 Jul 16, 2012:
      Add -errscale option to allow user to scale prior uncertainties by sample size. 
      Usage: -errscale <comparison SN sample size>
      Code will determine N = input SN sample size / comparison sample size,
      then adjust prior bao and/or cmb uncertainty as unc/sqrt(N). [JLM]

 Aug 20, 2012 - JLM
      Use fitcount variable to prevent errscale bug 
      (if -errscale AND -refit, errors were being rescaled twice).

 May 24, 2013 - JLM
      Fixed bug that hung code for empty fitres file.

 OCT 08, 2013 - JLM
      Adjusted code to allow proper fitting of binned SN data. 
      If NBIN column exists in fitres file ( i.e. SALT2train_biascorr.pl output),
      prior uncertainties are weighted by the sum of NBIN values rather than the
      number of "SNE". NB: intrinsic scatter output has not been adjusted to take
      binning into account. 
      Also, altered the chi-sq minimization "refined grid" stepsize to very as a
      function of input grid stepsize (was hardcoded to be 0.001 before, is now
      10x smaller than input grid stepsize). If default stepsizes are used,
      refined grid stepsize is unchanged.
      Using -marg flag now causes output w0 uncertainty to be the marg value.


 Aug 27 2015: read zHD and zHDERR with higher priority over z & zERR

 Apr 15 2016: if there is a MUDIF & MUDIFERR columm, then interpret
              as MU-MU(REF) based on wref and omref. See -wref and -omref
              options to specify MU(REF).  

 Jun 29 2016: in read_fitres(), float -> double for reading  inputs.
          
 Sep 20 2016: 
   + H0->70 (was 65)
   + for fitswrite option, write ascii in FITRES format using ROW keys

 June 2017:
   + optional argument "-label <label>" to label cospar outpout.
   + read MUERR_INCLUDE flags from fitres file to avoud
     double-counting MUERR contributions.

 Apr 3 2019: 
   + fix parsing to work withot NVAR key.
   + remove LEGACY option to read input file without keys.
     Now it only reads FITRES-formatted options.

 May 07, 2019: use MUREF column if it's there (for MUDIF option)

 Jun 13, 2019:
   + input syntax for cospar outfile is
     -cospar <fileName>  or
     -outfile_cospar <fileName>
   + same for resid file with key -resid or -outfile_resid
   + -csv command writes cospar in csv format.
   + -outfile_resid NONE --> no resid file

 Mar 17 2020: 
    + H0=70 replaced with H0 = H0_SALT2 (param in sntools.h)
    + omm_prior = 0.3 -> OMEGA_MATTER_DEFAULT (paramn in sntools.h)

 Mar 18 2020:
    + DJB added cmb_sim and bao_sim flags to set cosmology of priors to 
      same as default cosmology in sntools

 Aug 15 2020
   + new input  -cospar_yaml <outFile> to write results in yaml format.
     Motivation is for easier parsing with batch scripts.

 Sep 21 2020: few fixes for YAML output.
 Oct 17 2020: pass HzFUN typedef to refactored Hainv_integral function.
 
 Mar 27 2021: new -om_sim option to change OM for -cmb_sim.

*****************************************************************************/

int compare_double_reverse (const void *, const void *);
void writemodel(char*, float, float, float);
void printerror( int status);
void read_fitres(char *inFile);
void parse_VARLIST(FILE *fp);
void read_mucovar(char *inFile);
void invert_mucovar(double sqmurms_add);
void   get_chi2wOM(double w0, double wa, double OM, double sqmurms_add,
		   double *mu_off, double *chi2sn, double *chi2tot );// AM w0 --> w0
void getname(char *basename, char *tempname, int nrun);

double get_minwOM( double *w0_atchimin, double *wa_atchimin, double *OM_atchimin ); // AM w0 --> w0
void set_priors(void);

void set_HzFUN_for_wfit(double H0, double OM, double OE, double w0, double wa,
                        HzFUN_INFO_DEF *HzFUN ) ;



typedef struct {
  double  w0_out, w0sig, w0sig_upper, w0sig_lower ;
  double  wa_out, wasig, wasig_upper, wasig_lower ;
  double  omm_out, omm_sig,  chi2_final ;
  double  sigmu_int, w0rand, warand,  ommrand ;
  char    label_cospar[40];
  int     Ndof ;
} RESULTS_DEF ;  // Aug 15 2020

double *mu, *mu_sig, *mu_ref, *mu_sqsig, *z, *z_sig;
int    *snnbin; // to allow for binning of SNe. default is 1. [JLM]

double mu_offset ;
double omm_stepsize, w0_stepsize, wa_stepsize, h_stepsize;// w0_stepsoze --> w0_stepsize
int    cidindex(char *cid);

void write_output_cospar(FILE *fp, RESULTS_DEF *RESULTS, 
			 int usemarg, int format_cospar );
void write_output_resid(void);
void write_output_contour(void);

struct CUTSELECT {
  double ZMIN ;
  double ZMAX ;
  int    GENTYPE_SIM  ;
  int    SNTYPE       ;
  int    MXSNFIT ;
} CUTSELECT ;


// define CERN LIB function for matrix inversion
// extern void dfinv_( int *NCOVSN, double *covtmp, int *NMAX, int *IR ) ;

// =========== variables ================

// xxxdouble H0       = 70.0 ;  // approx value only (changed Sep 2016)
double H0       = H0_SALT2;      //
double H0SIG    = 1000.0 ;        // error used in prior
double SIG_MUOFF ;
double SQSIG_MUOFF ;

#define MXSN 100000 // max number of SN to read & fit
char CIDLIST[MXSN][20];

int NCIDLIST, NCUT ;
int NSNE_NBIN ; // for errscale 
int *tid, *GENTYPE_SIM, *SNTYPE ;
int snana = 0 ;
int snanasim = 0;


int fitnumber = 1;
int fitcount = 0;
double snrms   = 0.0 ;
double sqsnrms = 0.0 ;

int usebao  = 0;
int usecmb  = 0;
int usewa = 0 ; //AM
int csv_out = 0; // optional csv format for output cospar and resids
int mudif_flag = 0 ; // Apri 2016
int MUERR_INCLUDE_zERR;    // True if zERR already included in MUERR
int MUERR_INCLUDE_LENS ;   // True if lensing sigma already included

  /* WMAP + LSS, from SDSS, Tegmark et al, astro-ph/0310723  */

double OMEGA_MATTER_SIM = OMEGA_MATTER_DEFAULT ;
double omm_prior        = OMEGA_MATTER_DEFAULT ;      // OM prior
double omm_prior_sig    = 0.04;  /* 1 sigma uncertainty on prior: 10% */

// Apr 2016: reference cosmology to fit  MUDIF = mu- mu(ref), 
// only if MUDIF column is present in the input file
double wref  = -1.0 ;
double waref =  0.0 ;
double omref =  0.3 ;

char label_cospar[40] = "none" ; // string label for cospar file.

  /* Cosparam cosmological parameter structure */
Cosparam cpar;

double c_light = 299792.458 ;  /* speed of light, km/s */

  /* BAO parameters */
double abest = 0.469;
double sigma_a = 0.017;
double z_bao = 0.35;
double a1 ;

  // CMB params (May 29, 2008 RSK)
double zcmb = 1089. ;
double acmb ;
  //  double sigma_Rcmb = .021;  // used in Kowalski ??
  // double Rcmb_best = 1.715 ;  // idem

double sigma_Rcmb = 0.019;   // from Komatsu
double Rcmb_best  = 1.710 ;  // 

double TWOTHIRD = 2./3. ;
double TWO = 2. ;
double NEGTHIRD = -1./3. ;
double ONE = 1.0;
double ZERO = 0.0;


// May 20, 2009: mu-covariance structures


#define MXCOVSN   400                   // max # of SN with MU-covariance 
#define MXCOVPAIR MXCOVSN*(MXCOVSN-1)/2 // max off-diag covariance pairs 

int NCOVPAIR = 0 ;
int NCOVSN   = 0 ; 
int INDEX_COVSN_MAP[MXCOVSN];      // CIDLIST index vs. NCOVSN index
int INDEX_COVSN_INVMAP[MXSN];      // NCOVSN index vs. CIDLIST

struct MUCOV_INPUT {
  char   CID[2][40]; // name of each SN
  int    INDEX[2];   // index in CIDLIST
  double COV ;   // covariance
} MUCOV_INPUT[MXCOVPAIR];


struct MUCOV_STORE {
  int    INDEX[2];   // index in CIDLIST
  double COV ;     // covariance
  double COVINV;   // inverse
} MUCOV_STORE[MXCOVSN][MXCOVSN];

double *snprob, ***snprob3d, *extprob, ***extprob3d, *snchi, ***snchi3d, *extchi, ***extchi3d;

// =============================
int main(int argc,char *argv[]){

  /* Command line help */
  
  static char *help[] = {		/* instructions */
    "",
    "WFIT - Usage: wfit [hubble diagram] (options)",
    "",
    " Default for input file is 3 columns: z, mu, sigma_mu, but see below",
    "",
    "  Options:",
    "   -wa\t\tfit for w0-wa from Chevallier & Polarski, 2001 [Int.J.Mod.Phys.D10,213(2001)]",
    "   -ompri\tcentral value of omega_m prior [default: Planck2018]", 
    "   -dompri\twidth of omega_m prior [default: 0.04]",
    "   -bao\t\tuse BAO constraints from Eisenstein et al, 2006",
    "   -bao_sim\tuse BAO constraints with simulated cosmology and E06 formalism",
    "   -cmb\t\tuse CMB constraints from 5-year WMAP (2008)",
    "   -cmb_sim\tuse CMB constraints with simulated cosmology and WMAP formalism",
    "   -om_sim\tOmega_M for cmb_sim (default from sntools.h)",
    "   -minchi2\tget w and OM from minchi2 instead of marginalizing",
    "   -marg\tget w and OM from marginalizing",
    "   -Rcmb\tCMB comstraints: R = Rcmb +/- sigma_Rcmb [= 1.710 +/- 0.019]",
    "   -sigma_Rcmb\tUncertainty on Rcmb",
    "   -abest\tBAO constraints: A = abest +/- sigma_a [= 0.469 +/- 0.017]",
    "   -sigma_a\tSee Eisensten et al., '06, eqn. 4 for details",
    "   -z_bao\tBAO average redshift, z_bao=0.35",
    "   \t\tcurrently user must choose either Omega_m or BAO, not both",
    "   -czerr\terror in c*redshift (km/s), used for all SNe",
    "   -zerr\terror in redshift, used for all SNe",
    "   -snrms\tadd this to reported errors (mags) [default: 0]",
    "   -dz\t\twfit expects input file to be: z, sigma_z, mu, sigma_mu",
    "   -zmin\tFit only data with z>zmin",
    "   -zmax\tFit only data with z>zmax",    
    "   -GENTYPE\tSelect GENTYPE (for sim only)",
    "   -snana\tinput is snana .fitres format",
    "   -snanasim\tinput is snana .fitres format for simulated data",
    "   \t\tThis overrides redshift error given by czerr",
    "   -blind\tIf set, results are obscured with sin(random) ",
    "   \t\tlarge random integer is added to results to hide their values.",
    "   -fitswrite\tWrite 2D likelihoods to output fits file.",
    "   -gridwrite\tWrite 2D likelihoods as output file.",
    "   -mucovar\tUse distance covariances from this input file.",
    "   -refit\tfit once for sigint then refit with snrms=sigint.", 
    "   -errscale\trescale prior errors as 1/sqrt(N).. must include comp SN sample size",
    "",
    "   -wref \tfit for w-wref   (reads MUDIF column from SALT2mu)",    
    "   -omref\tfit for om-omref (reads MUDIF column from SALT2mu).",
    "",
    " Grid spacing:",
    " wCDM Fit:",
    "   -hmin/-hmax/-hsteps\t\tH0 grid [40,100,121]",
    "   -wmin/-wmax/-wsteps\t\tw grid  [0,-2,201]",
    "   -ommin/-ommax/-omsteps\tOM grid [0,1.0,81]",
    "",
    " w0-wa Fit:",
    "   -hmin/-hmax/-hsteps\t\tH0 grid [40,100,121]",
    "   -w0min/-w0max/-w0steps\tw0 grid [-3.0,1.0,201]",
    "   -wamin/-wamax/-wasteps\twa grid [-4.0,4.0,301]",
    "   -ommin/-ommax/-omsteps\tOM grid [0,1.0,81]",
    "",


    " Output:",
    "   -cospar\tname of output LEGACY file with fit cosmological params",
    "   -cospar_yaml\tname of output YAML file with fit cosmo params",
    "   -resid\tname of output file with mu-residuals",
    "   -chi2grid\tname of output file containing chi2-grid",
    "   -label\tstring-label for cospar file.",
    "",
    0
  };

  int i, j,kk, iarg ;
  FILE *fpcospar, *fpresid;
  char infile[2000];
  char cosparfile[1000] = "";
  char cosparfilevar[1000] = "";
  int  format_cospar = 1;  // 1=legacy, 2=YAML (RK Aug 15, 2020)
  char residfile[1000]= "";
  char residfilevar[1000]="";
  char chi2gridfile[1000]="";
  char chi2gridfilevar[1000]="";
  char chi2gridfile_SN[1000]="";
  char mucovarfile[1000]="";
  char tempfilename1[1000];
  char tempfilename2[1000];

  
  double *w0_prob,*wa_prob,  *w0_sort, *wa_sort,  *hvec, *omm_prob;
  double *chitmp;
  //double *snprob, ***snprob3d, *extprob, ***extprob3d, *snchi, ***snchi3d, *extchi, ***extchi3d;
  double f_memory; // allocated memory in Mb 
  double w0,wa, omm, ld_cos, mu_cos, coeff ;
  int  omm_steps, w0_steps,wa_steps,  h_steps;//Ayan Mitra
  double h_min, h_max, omm_max, omm_min, w0_max, w0_min, wa_min, wa_max;
  double w0_sum, w0_sort_1sigma, wa_sum, wa_sort_1sigma;
  double snchi_min=1.e20, extchi_min=1.e20;
  double w0_probsum, w0_probmax, w0_mean, wa_probsum, wa_probmax, wa_mean;
  //double w_probsum, w_probmax, w_mean, w_probsum, w_probmax, w_mean;
  double w0sig, w0sig_upper=0.0, w0sig_lower=0.0;
  double wasig, wasig_upper=0.0, wasig_lower=0.0;
  //double w0sig, w0sig_upper=0.0, w0sig_lower=0.0;
  double w0_out, wa_out, omm_out ;
  double Pt1mp,Pt2mp, P1max, P2max ;
  double omm_probsum, omm_mean, omm_sig;
  double delta_w0,delta_wa,  mu_dif ;
  double w0_atchimin,wa_atchimin,  OM_atchimin, chi2atmin=0.0;
  int iw0_mean=0, iwa_mean = 0;
  int imin=0, kmin=0, jmin=0;
  double snprobtot=0, snprobmax, extprobtot=0, extprobmax;

  double zerr=0.0, mu_sig_z;
  float  dif, mindif ;
  double sigmu_int=0.0, sigmu_int1=0.0, sigmu_tmp, chi2tmp ;
  double chi_approx, chi2_final;
  double chi2tot, chi2sn ;
  double w0rand, warand,ommrand;
  int Ndof;

  RESULTS_DEF RESULTS ;
  double muoff_tmp, sqmusig_tmp, musig_tmp, snchi_tmp, extchi_tmp, chidif; 
  
  /* flags */
  int dz=0;
  int blind   = 0;
  int fitsflag= 0;
  int gridflag= 0;
  int writechi= 0;
  int usemarg = 0;
  int usemucovar = 0;


  FILE *FILEPTR_TOT, *FILEPTR_SN ;

  /* Cfitsio things, for writing out probabilities*/

  fitsfile *fptr;
  /*   char snfitsname[] = "!snprob.fits"; */
  /*   char ommfitsname[] = "!extprob.fits"; */
  char snfitsname[1000], ommfitsname[1000], snfitsstr[1000], ommfitsstr[1000];
  long naxis=2;
  long fpixel[2], nelements, naxes[2];
  int status;

  /* for option -errscale */
  int Nerrscale=0;
  int NSNEerrscale;
  float SNsampleratio;

  // ----------------- BEGIN MAIN ----------------

  /** Give help if no arguments **/

  set_EXIT_ERRCODE(EXIT_ERRCODE_wfit);

  if (argc < 2) {
    for (i = 0; help[i] != 0; i++)
      printf ("%s\n", help[i]);
    exit(0);
  }

  // Aug 15 2020: print full command 
  printf(" Full command: \n   ");
  for(i=0; i < argc; i++ ) { printf(" %s", argv[i] );  }
  printf("\n\n"); fflush(stdout);

  /** Initialize parameters **/
  /* Cosmological params */

  //  test_codist();  only to test r(z) integral


  SIG_MUOFF   = 5.0 * log10(1. + H0SIG/H0);
  SQSIG_MUOFF = SIG_MUOFF * SIG_MUOFF ;

  double rz ;

  /* 2df prior, as used by JT's wcont */

  /*   We specify the *number* of steps for omega_m, w0, H0. */

  a1   = 1./(1. + z_bao);
  acmb = 1./(1.+zcmb);

  /* Omega_m & w0 grid parameters */
 
 /* number of steps across omega_m range */


  w0_steps = 201; w0_min = -2.0 ;  w0_max = 0. ;  
  wa_steps = 1;   wa_min = 0. ;  wa_max = 0. ;
  omm_steps = 81; omm_min = 0.0; omm_max = 1.0;
  h_steps = 121; h_min = 40; h_max = 100;    /* number of steps across H0 range */




  
  CUTSELECT.ZMIN = 0.0 ;
  CUTSELECT.ZMAX = 9.9 ;
  CUTSELECT.GENTYPE_SIM = -9 ;
  CUTSELECT.SNTYPE      = -9 ;
  CUTSELECT.MXSNFIT     = MXSN+1 ;

  /* Range of Hubble parameters to marginalize over */

  /** Parse the args **/
  strcpy(infile,argv[1]);

  for (iarg=2; iarg<argc; iarg++) {
    if (argv[iarg][0]=='-') {
      if (strcasecmp(argv[iarg]+1,"ompri")==0) { 
	omm_prior = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"dompri")==0) { 
	omm_prior_sig = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"bao")==0) { 
	usebao=1;
      } else if (strcasecmp(argv[iarg]+1,"bao_sim")==0) {
        usebao=2;
      } else if (strcasecmp(argv[iarg]+1,"csv")==0) { 
	csv_out=1;
      } else if (strcasecmp(argv[iarg]+1, "Rcmb")==0) {
	Rcmb_best = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1, "sigma_Rcmb")==0) {
	sigma_Rcmb = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"abest")==0) { 
	abest = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"sigma_a")==0) { 
	sigma_a = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"z_bao")==0) { 
	z_bao = atof(argv[++iarg]); a1 = 1./(1. + z_bao);
      } else if (strcasecmp(argv[iarg]+1,"cmb")==0) { 
	usecmb=1;
      } else if (strcasecmp(argv[iarg]+1,"wa")==0) {
        usewa=1; 
	wa_steps = 301;w0_steps = 201;   wa_min = -4. ;  wa_max = 4. ; w0_min = -3.; w0_max = 1; //-3,1 Later might add inputs to alter. 

      } else if (strcasecmp(argv[iarg]+1,"cmb_sim")==0) {
        usecmb=2;
	omm_prior     = OMEGA_MATTER_SIM ;
	omm_prior_sig = 0.5;  // turn off omm prior
      } else if (strcasecmp(argv[iarg]+1,"om_sim")==0) {
	OMEGA_MATTER_SIM = atof(argv[++iarg]); 
	omm_prior     = OMEGA_MATTER_SIM ;
      } else if (strcasecmp(argv[iarg]+1,"minchi2")==0) { 
	usemarg=0;
      } else if (strcasecmp(argv[iarg]+1,"marg")==0) { 
	usemarg=1;

      } else if (strcasecmp(argv[iarg]+1,"wref")==0) { 
	wref = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"omref")==0) { 
	omref = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"czerr")==0) { 
	zerr = atof(argv[++iarg])/c_light; /* take out c */
      } else if (strcasecmp(argv[iarg]+1,"zerr")==0) { 
	zerr = atof(argv[++iarg]); 
      } else if (strcasecmp(argv[iarg]+1,"snrms")==0) { 
	snrms = atof(argv[++iarg]);
	sqsnrms = snrms * snrms ;
      } else if (strcasecmp(argv[iarg]+1,"refit")==0){
        fitnumber = fitnumber - 1;
      } else if (strcasecmp(argv[iarg]+1,"errscale")==0){
	Nerrscale = atoi(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"dz")==0) { 
	dz=1;
      } else if (strcasecmp(argv[iarg]+1,"zmin")==0) { 
	CUTSELECT.ZMIN = atof(argv[++iarg]); 

      } else if (strcasecmp(argv[iarg]+1,"zmax")==0) { 
	CUTSELECT.ZMAX = atof(argv[++iarg]); 	

      } else if (strcasecmp(argv[iarg]+1,"SNTYPE")==0) { 
	CUTSELECT.SNTYPE = atoi(argv[++iarg]); 

      } else if (strcasecmp(argv[iarg]+1,"GENTYPE")==0) { 
	CUTSELECT.GENTYPE_SIM = atoi(argv[++iarg]); 

      } else if (strcasecmp(argv[iarg]+1,"MXSNFIT")==0) { 
	CUTSELECT.MXSNFIT = atoi(argv[++iarg]); 

      } else if (strcasecmp(argv[iarg]+1,"blind")==0) { 
	blind=1;
      } else if (strcasecmp(argv[iarg]+1,"snana")==0) { 
	snana=1;
      } else if (strcasecmp(argv[iarg]+1,"snanasim")==0) { 
	snana=1;
	snanasim=1;
	
      } else if (strcasecmp(argv[iarg]+1,"fitswrite")==0) {   
	/* Output likelihoods as fits files & text files */
	fitsflag=1;  gridflag = 1;

      } else if (strcasecmp(argv[iarg]+1,"gridwrite")==0) {   
	/* Output likelihoods as fits files */
	gridflag=1; 

      } else if (strcasecmp(argv[iarg]+1,"mucovar")==0) {   
  	strcpy(mucovarfile,argv[++iarg]);
	usemucovar =1 ;

      /* Change H0 grid parameters? */
      } else if (strcasecmp(argv[iarg]+1,"hmin")==0) {   
	h_min=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"hmax")==0) {   
	h_max=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"hsteps")==0) {   
	h_steps = (int)atof(argv[++iarg]);
      }

      /* Change Omega_m grid parameters? */
      else if (strcasecmp(argv[iarg]+1,"ommin")==0) {   
	omm_min=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"ommax")==0) {   
	omm_max=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"omsteps")==0) {   
	omm_steps = (int)atof(argv[++iarg]);
      }

      /* Change w0 grid parameters? */
      else if (strcasecmp(argv[iarg]+1,"w0min")==0) {   
	w0_min=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"w0max")==0) {   
	w0_max=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"w0steps")==0) {   
	w0_steps = (int)atof(argv[++iarg]);
      }


      /* Change wa grid parameters? */
      else if (strcasecmp(argv[iarg]+1,"wamin")==0) {
        wa_min=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"wamax")==0) {
        wa_max=atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"wasteps")==0) {
        wa_steps = (int)atof(argv[++iarg]);
	}//*/

      
      
      /* Write out chi2 instead probabilities */
      else if (strcasecmp(argv[iarg]+1,"writechi")==0) {   
	writechi=1;
      }

      /* Output filenames */
      else if (strcasecmp(argv[iarg]+1,"outfile_cospar")==0) 
 	{ strcpy(cosparfilevar, argv[++iarg]); }      
      else if (strcasecmp(argv[iarg]+1,"cospar")==0) 
 	{ strcpy(cosparfilevar, argv[++iarg]); }      

      else if (strcasecmp(argv[iarg]+1,"cospar_yaml")==0) 
 	{ strcpy(cosparfilevar, argv[++iarg]); format_cospar=2; }      

      else if (strcasecmp(argv[iarg]+1,"label")==0) 
 	{ sprintf(label_cospar,"%s", argv[++iarg]); }      

      else if (strcasecmp(argv[iarg]+1,"outfile_resid")==0)
	{ strcpy(residfilevar,argv[++iarg]); }
      else if (strcasecmp(argv[iarg]+1,"resid")==0) 
	{ strcpy(residfilevar,argv[++iarg]);  }

      else if (strcasecmp(argv[iarg]+1,"chi2grid")==0)  
	{ strcpy(chi2gridfilevar,argv[++iarg]); }

      else if (strcasecmp(argv[iarg]+1,"outfile_chi2grid")==0)  
	{ strcpy(chi2gridfilevar,argv[++iarg]); }
      
      else {
	printf("Bad arg: %s\n", argv[iarg]);
	exit(EXIT_ERRCODE_wfit);
      }
    }
  }


  /*** Set up vectors and matrices ***/
  /* Data vectors */
  mu          = (double *)calloc(MXSN,sizeof(double));
  mu_sig      = (double *)calloc(MXSN,sizeof(double));
  mu_ref      = (double *)calloc(MXSN,sizeof(double));
  mu_sqsig    = (double *)calloc(MXSN,sizeof(double));

  
  z      = (double *)calloc(MXSN,sizeof(double));
  z_sig  = (double *)calloc(MXSN,sizeof(double));
  chitmp = (double *)calloc(MXSN,sizeof(double));
  tid    = (int    *)calloc(MXSN,sizeof(int));
  GENTYPE_SIM = (int*)calloc(MXSN,sizeof(int));
  SNTYPE      = (int*)calloc(MXSN,sizeof(int));
  snnbin      = (int*)calloc(MXSN,sizeof(int));

  w0_prob   = (double *)calloc(w0_steps,sizeof(double));
  wa_prob   = (double *)calloc(wa_steps,sizeof(double));
  w0_sort   = (double *)calloc(w0_steps,sizeof(double));
  wa_sort   = (double *)calloc(wa_steps,sizeof(double));
  hvec     = (double *)calloc(h_steps,sizeof(double));
  omm_prob = (double *)calloc(omm_steps,sizeof(double));

  /* Probability arrays */
  extprob = (double *) calloc(omm_steps*w0_steps,sizeof(double));
  snprob  = (double *) calloc(omm_steps*w0_steps,sizeof(double));
  extchi  = (double *) calloc(omm_steps*w0_steps,sizeof(double));
  snchi   = (double *) calloc(omm_steps*w0_steps,sizeof( double));

  f_memory = 0;
 
  f_memory += malloc_double3D(+1, w0_steps, wa_steps, omm_steps,&snchi3d );
  f_memory += malloc_double3D(+1, w0_steps, wa_steps, omm_steps,&extchi3d );
  f_memory += malloc_double3D(+1, w0_steps, wa_steps, omm_steps,&extprob3d );
  f_memory += malloc_double3D(+1, w0_steps, wa_steps, omm_steps,&snprob3d );


  for(i=0; i < w0_steps; i++ ) {
    for(kk=0; kk < wa_steps; kk ++ ) {
      for (j=0; j < omm_steps; j++ ) {
	snchi3d[i][kk][j] = 0.0	;
	extchi3d[i][kk][j] = 0.0 ;
	extprob3d[i][kk][j] = 0.0 ;
	snprob3d[i][kk][j] = 0.0 ;  
      }
    }
  }
  
  printf(" =============================================== \n");
  printf(" SNANA_DIR = %s \n", getenv("SNANA_DIR") );

  while (fitnumber <= 1) {

    /************************/
    /*** Read in the data ***/
    /************************/

    read_fitres(infile);  // local function

    /*******************************************/
    /*** Adjust prior uncerts (if specified) ***/
    /*******************************************/

    if (Nerrscale > 0 && fitcount == 0){

      // IF input data file has an NBIN column (as is produced by SALT2train_biascorr.pl, 
      // which averages distances in redshift bins before applying bias corrections) 
      // the sum of NBIN column values will be used to weight the prior. Otherwise, 
      // the number of SNe in the data sample (NCIDLIST) will be used for weighting. 
      printf("\n\nAdjusting priors -- \n");
      printf("\t\t  NCIDLIST: %d, NSNE_NBIN: %d (based on NBIN column values) \n", 
	     NCIDLIST, NSNE_NBIN);     
      if (NSNE_NBIN == 0 ) {
	NSNEerrscale = NCIDLIST; 
      } else { 
	NSNEerrscale = NSNE_NBIN;
      }
      printf("\n\n RESCALING BAO, CMB PRIOR UNC'S BY SAMPLE SIZE \n");
      SNsampleratio = (1.0*NSNEerrscale)/(1.0*Nerrscale);
      printf("   Nerrscale = %i\n", Nerrscale);
      printf("   NSNe      = %i\n", NSNEerrscale);
      if (usebao){
        printf("   Rescale bao sigma_a by sqrt(%0.3f)\n\n", SNsampleratio);
	sigma_a = sigma_a/sqrt(SNsampleratio);
      } else {
	printf("   Rescale flat prior uncertainty by sqrt(%0.3f)\n\n", SNsampleratio);
	omm_prior_sig = omm_prior_sig/sqrt(SNsampleratio);
      }
      if (usecmb){
        printf("   Rescale cmb sigma_Rcmb by sqrt(%0.3f)\n\n", SNsampleratio);
	sigma_Rcmb = sigma_Rcmb/sqrt(SNsampleratio);
      }
    }

    /********************************************/
    /*** Add redshift error to distance error ***/
    /********************************************/

    coeff = 5./log(10.);
    for (i=0; i < NCIDLIST; i++){

      // Add any additional z error, adjust for binning if any 
      if (zerr > 1.0E-12 && MUERR_INCLUDE_zERR ){
	z_sig[i] = sqrt(z_sig[i]*z_sig[i] + zerr*zerr/snnbin[i]);
	
	if (i==0) 
	  printf(" Adding redshift uncertainty of %7.4f.\n",zerr);
      }

      // Skip if zERR is already included in MUERR (Jun 2017)
      if ( MUERR_INCLUDE_zERR == 0 ) {	
	mu_sig_z    = coeff * ((1+z[i])/(z[i]*(1+z[i]/2))) * z_sig[i];
	mu_sig[i]   = sqrt(mu_sig[i]*mu_sig[i] + mu_sig_z*mu_sig_z);
      }
      mu_sqsig[i] = mu_sig[i] * mu_sig[i] ;
    }
    

    // April 2016
    // if fitting MUDIF, then add reference cosmology to each MUDIF
    if ( mudif_flag ) {
      Cosparam cparref ;     
      cparref.omm = omref ;
      cparref.ome = 1.0 - omref ;
      cparref.w0  = wref ;
      cparref.wa  = waref; //0.0 ;
    
      double rz, ld_cos, mudif, muref ;
      for (i=0; i < NCIDLIST; i++){
	rz       = codist(z[i], &cparref);
	ld_cos   = (1+z[i]) *  rz * c_light / H0;
	muref    = mu_ref[i]; // May 7 2019
	mudif    = mu[i] ;
	mu[i]   += muref ;
	/*
	printf(" xxx %d: z=%.3f  MU=%.5f + %.5f = %.5f\n",
	       i, z[i], mudif, muref, mu[i] ); fflush(stdout);
	*/
      }
    }

    // Set BAO and CMB priors
    set_priors();
  
  /********************************************************/
    /******* read MU-covariances (if specified) *************/
    /********************************************************/

    for ( i=0; i<MXSN; i++ )  { INDEX_COVSN_INVMAP[i] = -9 ; }
  
    if ( usemucovar > 0 ) {
      read_mucovar(mucovarfile);

      printf(" Invert MU-covariance matrix with CERN's dfinv. \n");
      invert_mucovar(sqsnrms);
    }

    fflush(stdout) ;

    /****************************************/
    /*** Calculate chi-squared  on a grid ***/
    /****************************************/
    //wa_steps = 51;wa_max = 1; wa_min = -1; // XXX check
    snprobtot=0.0;
    snprobmax=0.0;
    extprobtot=0.0;
    extprobmax=0.0;
    w0_stepsize = wa_stepsize = omm_stepsize = h_stepsize = 0. ;

    if (w0_steps > 1) {w0_stepsize = (w0_max-w0_min)/(w0_steps-1);}
    if (wa_steps > 1) {wa_stepsize = (wa_max-wa_min)/(wa_steps-1);}
    if (omm_steps > 1) {omm_stepsize = (omm_max-omm_min)/(omm_steps-1);}
    if (h_steps > 1) {h_stepsize = (h_max-h_min)/(h_steps-1);}
    
 
    printf("\n\n****************************************\n");
    printf("****************************************\n");
    printf("   Applying  w(z) = w0 + wa(1-a) \n");
    printf("   w ~ w0. wa implemented.       \n");
    printf("****************************************\n");
    printf("****************************************\n\n");
    printf("   --------   Grid parameters   -------- \n");
    printf("  w0_min: %6.2f   w0_max: %6.2f  %5i steps of size %8.5f\n",
	    w0_min,w0_max,w0_steps,w0_stepsize);

    if(usewa){
    printf("  wa_min: %6.2f   wa_max: %6.2f  %5i steps of size %8.5f\n",
            wa_min,wa_max,wa_steps,wa_stepsize);
    }
    printf(" omm_min: %6.2f  omm_max: %6.2f  %5i steps of size %8.5f\n",
	    omm_min,omm_max,omm_steps,omm_stepsize);
    printf("   h_min: %6.2f   h_max: %6.2f  %5i steps of size %8.5f\n",
	      h_min,  h_max,  h_steps,  h_stepsize);
    printf("   ------------------------------------\n");

    /* Report priors being used */

    if ( H0SIG < 100. ) 
      printf("Fit with H0 prior: sigH0/H0=%4.1f/%4.1f  => sig(MUOFF)=%5.3f \n", 
	     H0SIG, H0, SIG_MUOFF );

    if (usebao) {
      if (usebao == 1) {
	printf("Fit with BAO Eisenstein 2006 prior: "
	       "A(BAO) =%5.3f +- %5.3f \n", abest, sigma_a);
      }
      if (usebao == 2) {
        printf("Fit with BAO prior in sim cosmology: "
	       "OM=%5.3f, w=%5.3f, A(BAO) =%5.3f +- %5.3f \n" ,
               OMEGA_MATTER_SIM, w0_DEFAULT, abest, sigma_a);
      } 
    } else {
      printf("Fit data with Omega_m prior: %5.3f +/- %5.3f\n",
	      omm_prior,omm_prior_sig);
    }

    if (usecmb) {
      if ( usecmb == 1 ) {
	printf("Fit data with CMB (WMAP) prior:  R=%5.3f +- %5.3f  \n" ,
	        Rcmb_best, sigma_Rcmb);
      }
      if ( usecmb == 2 ) {
        printf("Fit data with CMB (WMAP) prior in sim cosmology: "
	       "OM=%5.3f, w=%5.3f, R=%5.3f +- %5.3f  \n" ,
               OMEGA_MATTER_SIM, w0_DEFAULT, Rcmb_best, sigma_Rcmb);
      } 
      //debugexit('hello');
    }

    if(usewa){
      if ( usemarg != 0 ) 
	printf("Will MARGINALIZE for final (w0, wa) & OM \n");
      else
	printf("Will MINIMIZE for final (w0, wa) & OM \n");
    }
    else {
      if ( usemarg != 0 )
	printf("Will MARGINALIZE for final w0 & OM \n");
      else
	printf("Will MINIMIZE for final w0 & OM \n");
    }

    printf("Add %4.2f mag-error (snrms) in quadrature to MU-error \n", snrms);
    printf("---------------------------------------\n");
    fflush(stdout);

    /* Get approximate expected minimum chi2 (= NCIDLIST - 3 dof),
       used to keep numbers small in the chi2 loop.  */

    //mark delete line below
    //Ndof = NCIDLIST - 3 + usebao + usecmb ;
    
    Ndof = NCIDLIST - 3 ;
    if ( usebao ) 
      Ndof++ ;
    if ( usecmb ) 
      Ndof++ ;

    

    chi_approx = (double)(Ndof);

    for( i=0; i < w0_steps; i++){
      cpar.w0 = w0_min + i*w0_stepsize;
      for( int kk=0; kk < wa_steps; kk++){    // indent
	cpar.wa = (wa_min + kk*wa_stepsize);// we will simply incremenet wa without 'a' at the moment *(1-(a_min + k*z_stepsize))	
	for(j=0; j < omm_steps; j++){
	  cpar.omm = omm_min + j*omm_stepsize; // AM look here
	  cpar.ome = 1 - cpar.omm;

	  get_chi2wOM ( cpar.w0,cpar.wa, cpar.omm, sqsnrms, 
			&muoff_tmp, &snchi_tmp, &extchi_tmp );   // return args. inputs  AM --> insert cpar.wa ?

	/*printf("xxxchi %d   %0.4f   %0.4f   %0.4g   %0.4g   %0.4g   \n", 
	  i, cpar.w0, cpar.omm, snchi_tmp, extchi_tmp, diff); */

	  //snchi[i*omm_steps+j]  = snchi_tmp ;
	  //extchi[i*omm_steps+j] = extchi_tmp ;
	  snchi3d[i][kk][j]  = snchi_tmp ; // AM
	  extchi3d[i][kk][j] = extchi_tmp ;
	  
        /* Keep track of minimum chi2 */
	  if(snchi_tmp < snchi_min) 
	    snchi_min = snchi_tmp ;

	  if(extchi_tmp < extchi_min) 
	    { extchi_min = extchi_tmp ;  imin=i; jmin=j;kmin=kk; }
	} // j loop
      }  // end of k-loop
    }  // end of i-loop
  
    // get w,OM at min chi2 by using more refined grid
    // Pass approx w,OM,  then return w,OM at true min
    //  printf("  Get minimized w,OM from refined chi2grid \n");
    w0_atchimin  = w0_min   + imin*w0_stepsize ;
    wa_atchimin = wa_min  + kmin*wa_stepsize ;
    OM_atchimin = omm_min + jmin*omm_stepsize ;

    fflush(stdout);

    if ( usemarg == 0 ) 
      chi2atmin   = get_minwOM( &w0_atchimin,&wa_atchimin,  &OM_atchimin );

    fflush(stdout);

    /*******************************************/
    /*** Normalize probability distributions ***/
    /*******************************************/
    
    /* First, convert chi2 to likelihoods */
    // for(i=0; i<w0_steps; i++){
    // for(j=0; j<omm_steps; j++){

       ///* Probability distribution from SNe alone */
    // chidif = snchi[i*omm_steps+j] - snchi_min ;
    //  snprob[i*omm_steps+j] = exp(-0.5*chidif);
    //  snprobtot += snprob[i*omm_steps+j];
      
	/*Probability distribution of SNe + external prior */
    //  chidif = extchi[i*omm_steps+j] - extchi_min ;
    //  extprob[i*omm_steps+j] = exp(-0.5*chidif);
    //	printf("\n XXX1 extprob %f, chidif = %f, extprob3d = %f\n",extchi[i*omm_steps+j], chidif, exp(-0.5*chidif));
    //	extprobtot += extprob[i*omm_steps+j];
    // }
    // }

    /* Now normalize so that these sum to 1.  Note that if the
       grid does not contain all of the probability, 
       then these normalizations are incorrect */

    // for(i=0; i<w0_steps; i++){
    //  for(j=0; j<omm_steps; j++){
    //  snprob[i*omm_steps+j] /= snprobtot;
    //  extprob[i*omm_steps+j] /= extprobtot;
	/*printf("xxxprob %d   %0.4f   %0.4f   %0.4g   %0.4g   \n", 
	  i, w0_min + i*w0_stepsize, omm_min + j*omm_stepsize, 
	  snprob[i*omm_steps+j], extprob[i*omm_steps+j]); */
    //    }
    // }
       //----------------- New block below


        /* First, convert chi2 to likelihoods */
    for(i=0; i<w0_steps; i++){
      for(kk = 0; kk<wa_steps; kk++){
	for(j=0; j<omm_steps; j++){

	  /* Probability distribution from SNe alone */
	  chidif = snchi3d[i][kk][j] - snchi_min ;
	  snprob3d[i][kk][j] = exp(-0.5*chidif);
	  //printf("\n XXX snchi3d %f, chidif = %f, snprob3d = %f\n",snchi3d[i][kk][j], chidif, exp(-0.5*chidif));
	  snprobtot += snprob3d[i][kk][j];
      
	  /* Probability distribution of SNe + external prior */
	  chidif = extchi3d[i][kk][j] - extchi_min ;
	  extprob3d[i][kk][j] = exp(-0.5*chidif);
	  //	  printf("\n XXX extprob %f, chidif = %f, extprob3d = %f\n",extchi3d[i][kk][j], chidif, exp(-0.5*chidif));
	  extprobtot += extprob3d[i][kk][j];
	}
      }
    }

    /* Now normalize so that these sum to 1.  Note that if the
       grid does not contain all of the probability, 
       then these normalizations are incorrect */

    for(i=0; i<w0_steps; i++){
      for(kk = 0; kk<wa_steps; kk++){
	for(j=0; j<omm_steps; j++){
	  snprob3d[i][kk][j] /= snprobtot;
	  extprob3d[i][kk][j] /= extprobtot;
	/*printf("xxxprob %d   %0.4f   %0.4f   %0.4g   %0.4g   \n", 
	  i, w0_min + i*w0_stepsize, omm_min + j*omm_stepsize, 
	  snprob[i*omm_steps+j], extprob[i*omm_steps+j]); */
	}
      }
    }



    /***************************************/
    /*** Marginalize over omega_m, get w ***/
    /***************************************/    
    printf("---------------------------------------\n");
    printf("Marginalizing over omega_m...\n");
    w0_probsum = 0.0 ;wa_probsum = 0.0 ;
    w0_mean    = 0.0 ;wa_mean    = 0.0 ;
    w0_probmax = 0.0 ;wa_probmax = 0.0 ;
    P1max      = 0.0 ;P2max = 0.0;
    

    // Get Marginalise w0
    for(i=0; i < w0_steps; i++){
      w0 = w0_min + i*w0_stepsize;
      w0_prob[i] = 0. ; // initialize (RK, Jan 16, 2009)                                                                                    
      for(kk=0; kk < wa_steps; kk++){
	for(j=0; j < omm_steps; j++){
	  Pt1mp       = extprob3d[i][kk][j] ;
	  w0_prob[i] += Pt1mp ;
	  if ( Pt1mp > P1max ) {
	    P1max = Pt1mp ;
	    w0_probmax = w0 ; // for test only                                                                                        
	  } // if Loop
	} // j
      } // kk
	// printf(" %f  %f  \n", w0, 40.*w0_prob[i]) ;                                                                                       
      w0_mean    += w0_prob[i]*w0;
      w0_probsum += w0_prob[i];
      //printf("\n XXX w0_probsm = %f\n",w0_probsum);
      
    } // i
    /** Mean of w0 distribution: this is our estimate of w0 **/
    w0_mean /= w0_probsum;
	
    // Get Marginalise wa 	
	
    if(usewa){
      for(kk=0; kk < wa_steps; kk++){
	wa = wa_min + kk*wa_stepsize;
	wa_prob[kk] = 0. ; 
	for(i=0; i < w0_steps; i++){
	  for(j=0; j < omm_steps; j++){
	    Pt2mp       = extprob3d[i][kk][j] ;
	    wa_prob[kk] += Pt2mp ;
	    if ( Pt2mp > P2max ) {
	      P2max = Pt2mp ;
	      wa_probmax = wa ; 
	    } // if Loop
	  } // j
	} // i
	  //    printf(" %f  %f  \n", wa, 40.*wa_prob[i]) ;                                                                                
	wa_mean    += wa_prob[kk]*wa;
	wa_probsum += wa_prob[kk];
      } // kk
    /** Mean of wa distribution: this is our estimate of wa **/
      wa_mean /= wa_probsum;   
    }

    //---------------------
    //printf("XXX w0_mean = %f\n", w0_mean);
    if ( !blind ){
      if (usemarg){
       printf("  CHECK:  w0(marg) = %f \n",  w0_mean);
	if (usewa)
	  {printf("  CHECK:  wa(marg) = %f \n",  wa_mean );}
      }
      else
	{
	  printf("  CHECK:  ChiSq(min) = %f \n", chi2atmin);
	  printf("  CHECK:  w0(min)    = %f, %f\n",  w0_atchimin);
	  if (usewa)
	    {printf("  CHECK:  wa(min)    = %f, %f \n",  wa_atchimin );}  
	}
    }
    printf("\n" );


    /** Get std dev for the weighted mean.  This is a reasonable                                                                             
       measure of uncertainty in w if distribution is ~Gaussian **/
    w0sig=0;
    for(i=0; i<w0_steps; i++){
      w0 = w0_min + i*w0_stepsize;
      w0_prob[i] /= w0_probsum;   /** normalize the distribution **/
      w0sig += w0_prob[i]*pow((w0-w0_mean),2);
      //printf("XXX  w0sig=%f, w0_prob=%f, w0=%f, w0_mean=%f\n ",w0sig, w0_prob[i],w0, w0_mean);    
      w0_sort[i] = w0_prob[i];        /** make a copy to use later **/
      //      printf("\n XXX w0_sort = %f\n",w0_sort[i]);
    }
    w0sig = sqrt(w0sig/w0_probsum);
    printf("marg w0 err estimate = %f\n", w0sig);
    
    if (usewa){
      
      wasig=0;
      for(kk=0; kk<wa_steps; kk++){
	wa = wa_min + kk*wa_stepsize;
	wa_prob[kk] /= wa_probsum;   // normalize the distribution 
	wasig += wa_prob[kk]*pow((wa-wa_mean),2);
	wa_sort[kk] = wa_prob[kk];        // make a copy to use later 
      }
      wasig = sqrt(wasig/wa_probsum);
      printf("marg wa err estimate = %f\n", wasig); 
    } 
       
    
    /** Find location of mean in the prob vector **/
    delta_w0=1e3;

    for (i=0; i<w0_steps; i++){
      w0 = w0_min + i*w0_stepsize;
      if (fabs(w0 - w0_mean) < delta_w0){
        delta_w0 = fabs(w0 - w0_mean);
        iw0_mean = i;
      }
    }

   
    if (usewa){
      delta_wa=1e3;
      for (kk=0; kk<wa_steps; kk++){
	wa = wa_min + kk*wa_stepsize;
	if (fabs(wa - wa_mean) < delta_wa){
	  delta_wa = fabs(wa - wa_mean);
	  iwa_mean = kk;
	}
      }
    }
    
    /** Sort probability **/
    qsort(w0_sort,w0_steps,sizeof(double),compare_double_reverse);
    if (usewa) {qsort(wa_sort,wa_steps,sizeof(double),compare_double_reverse);}
    
    /** Add up probability until you get to 68.3%, use that value to 
       determine the upper/lower limits on w **/
    w0_sum=0.;
    w0_sort_1sigma=-1;
    for (i=0; i<w0_steps; i++){
      w0_sum += w0_sort[i];
      //     printf("\n XXX w0_sum = %f, %i \n",w0_sum,i);
      if (w0_sum > 0.683) {
	//printf("\n XXX w0_sort[i] = %f\n",w0_sort[i]);
        w0_sort_1sigma = w0_sort[i];  /** w0_sort value corresponding to 68.3%  **/
        break;
      }
    }

    if (usewa){
      wa_sum=0.;
      wa_sort_1sigma=-1;
      for (kk=0; kk<wa_steps; kk++){
	wa_sum += wa_sort[kk];
	if (wa_sum > 0.683) {
	  //printf("\n XXX wa_sort[i] = %f\n",wa_sort[kk]);
	  wa_sort_1sigma = wa_sort[kk];  /** w0_sort value corresponding to 68.3%  **/
	  break;
	}
      }
    }
    
    if (w0_sort_1sigma < 0 ){
      printf("ERROR: w0 grid does not enclose 68.3 percent of probability!\n");
      exit(EXIT_ERRCODE_wfit);
    }

    if (usewa){
      if (wa_sort_1sigma < 0 ){
	printf("ERROR: wa grid does not enclose 68.3 percent of probability!\n");
	exit(EXIT_ERRCODE_wfit);
      }
    }

   
    /** Prob array is in order of ascending w.  Count from i=iw0_mean
       towards i=0 to get lower 1-sigma bound on w **/
    for (i=iw0_mean; i>=0; i--){
      if (w0_prob[i] < w0_sort_1sigma){
        w0sig_lower = w0_mean - (w0_min + i*w0_stepsize);
        break;
      }
    }
    if(usewa){
      for (kk=iwa_mean; kk>=0; kk--){
	//printf( "XXX %f\n",wa_prob[kk]);
	if (wa_prob[kk] < wa_sort_1sigma){
	  //printf( "XXX %f",wa_prob[kk]);
	  wasig_lower = wa_mean - (wa_min + kk*wa_stepsize);
	  break;
	}
      }
    }
    
    /** Error checking **/
    if(i==0){
      printf("WARNING: lower 1-sigma limit outside range explored\n");
      w0sig_lower = 100;
    }

    if (usewa){
      if(kk==0){
	printf("WARNING: lower 1-sigma limit outside range explored\n");
	wasig_lower = 100;
      }
    }
    if (w0sig_lower <= w0_stepsize){                                              
      printf("WARNING: 1. w0 grid is too coarse to resolve lower 1-sigma limit\n");
      w0sig_lower = w0_stepsize;
    }
    if (usewa){
      if (wasig_lower <= wa_stepsize){                                                                    
      printf("WARNING: 1. wa grid is too coarse to resolve lower 1-sigma limit\n");
      wasig_lower = wa_stepsize;
      }
    }

    /** Count from i=iw0_mean towards i=w0_steps to get                                                                                                      
        upper 1-sigma bound on w **/
    for (i=iw0_mean; i<w0_steps; i++){
      if (w0_prob[i] < w0_sort_1sigma){
        w0sig_upper = (w0_min + i*w0_stepsize) - w0_mean;
        break;
      }
    }
    if(usewa){
    for (kk=iwa_mean; kk<wa_steps; kk++){
      if (wa_prob[kk] < wa_sort_1sigma){
        wasig_upper = (wa_min + kk*wa_stepsize) - wa_mean;
        break;
      }
    }
   }
    /** Error checking **/
    if(i==(w0_steps-1)){
      printf("WARNING: upper 1-sigma limit outside range explored\n");
      w0sig_lower = 100;
    } 
  

    if (w0sig_upper <= w0_stepsize){
      printf("WARNING: 2. w0 grid is too coarse to resolve upper 1-sigma limit\n %f, %f\n", w0sig_upper, w0_stepsize);
      w0sig_upper = w0_stepsize; 
    }

    if(usewa){
      if(kk==(wa_steps-1)){
	printf("WARNING: upper 1-sigma limit outside range explored\n");
	wasig_lower = 100;
	
      }
      if (wasig_upper <= wa_stepsize){
	printf("WARNING: 2. wa grid is too coarse to resolve upper 1-sigma limit\n %f, %f\n", wasig_upper, wa_stepsize);
	wasig_upper = wa_stepsize;
      }
    }    
    printf("\n---------------------------------------\n");
    printf("probw0-err estimates: lower = %f, upper = %f\n", w0sig_lower, w0sig_upper);
    if(usewa){
      //printf("probw0-err estimates: lower = %f, upper = %f\n", w0sig_lower, w0sig_upper);
      printf("probwa-err estimates: lower = %f, upper = %f\n", wasig_lower, wasig_upper);
    }

    /***************************************/
    /*** Marginalize over w, get omega_m ***/
    /***************************************/
    omm_probsum=0;
    omm_mean=0;
    /*  // XXX Mark Delete Ayan Mitra 21/07/2021
    for(j=0; j<omm_steps; j++){
      omm = omm_min + j*omm_stepsize;

      for(i=0; i<w0_steps; i++){
	printf("XXX om probsum = %f\n", extprob[i*omm_steps+j]);
        omm_prob[j] += extprob[i*omm_steps+j];
      }

      omm_mean += omm_prob[j]*omm;
      omm_probsum += omm_prob[j];
      }*/  // XXX End Mark DELETE  

    
      for(j=0; j<omm_steps; j++){
	omm = omm_min + j*omm_stepsize;
	for(kk=0; kk<wa_steps; kk++){
	  for(i=0; i<w0_steps; i++){
	    //printf("XXX om probsum = %f\n", extprob[i*omm_steps+j]);
	    omm_prob[j] +=extprob3d[i][kk][j];// extprob[i*omm_steps+j];
	      }
	}
      omm_mean += omm_prob[j]*omm;
      omm_probsum += omm_prob[j];
      //printf("2. XXX omm mean = %f, omm_prob = %f, omm = %f\n", omm_mean, omm_prob[j], omm);

      }
    


    /** Mean of omega_m distribution: this is our estimate of omega_m **/

      omm_mean /= omm_probsum;

    /** Get std dev for the weighted mean. This should be a sufficient estimate
        of the uncertainty in omega_m. **/
    omm_sig=0;
    for(i=0; i<omm_steps; i++){ 	
      omm = omm_min + i*omm_stepsize;    
      omm_prob[i] /= omm_probsum;	/** normalize the distribution **/     
      omm_sig += omm_prob[i]*pow((omm-omm_mean),2);
    }

    omm_sig = sqrt(omm_sig/omm_probsum);
    printf("marg ommerr estimate = %f\n", omm_sig);
    printf("\n---------------------------------------\n");

    /********************************************************/
    /*** Write out distance residuals from best cosmology ***/
    /********************************************************/

    /* Kludge: Because the H0 marginalization happens right away in the 
       above, we can't do the full marginalization here.  Instead, 
       we can analytically calculate the value of the 'offset' between
       mu_theory and the data for the fit values of the parameters.  
       Note that this 'offset' term is not directly interprettable as
       -5*log10(H0) because the distances reported by MLCS are based on
       assumed values for H0 and the peak mag of an SN. */
 
    if ( usemarg ) 
      { w0_out = w0_mean ;  wa_out = wa_mean;   omm_out = omm_mean ; }
    else
      { w0_out = w0_atchimin ;  wa_out = wa_atchimin ; omm_out = OM_atchimin ; }

    cpar.omm = omm_out;
    cpar.w0  = w0_out;
    cpar.ome = 1. - cpar.omm;
    cpar.wa  = wa_out;//cpar.wa; //0.0 ;

    // get mu-offset from final parameters.
    get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, sqsnrms,  // inputs
	  	&mu_offset, &snchi_tmp, &chi2_final );   // return args

    /* Print residuals to file */
    if ( !IGNOREFILE(residfilevar) ) {
      if ( strlen(residfilevar) == 0 ) {
	strcpy(residfilevar, infile);
	strcat(residfilevar, ".resid");
      }

      getname(residfilevar, tempfilename1, fitnumber);
      strcpy(residfile, tempfilename1);

      printf("Write MU-residuals to %s with MU_offset=%6.3f \n",
	     residfile, mu_offset);
      fpresid=fopen(residfile, "w");
      
      if (fpresid == NULL){
	printf("ERROR: couldn't open %s\n",residfile);
	exit(EXIT_ERRCODE_wfit);
      }
      
      fprintf(fpresid,"z, mu_dif,  mu_sig,  tel_id,  CID \n"); // csv format
      for (i=0; i<NCIDLIST; i++){
	rz     = codist(z[i], &cpar) ;
	ld_cos = (1+z[i]) *  rz * c_light / H0;
	mu_cos =  5.*log10(ld_cos) + 25. ;
	mu_dif =  mu_cos - mu[i] - mu_offset ;
	sqmusig_tmp  =  mu_sqsig[i] + sqsnrms ;
	musig_tmp    = sqrt(sqmusig_tmp);
	
	fprintf(fpresid,"%7.4f %7.4f %7.4f %4i  %6s \n",
		z[i], mu_dif, musig_tmp, tid[i], CIDLIST[i]);
      }
      fclose(fpresid);
    } // end resid-write 

    // --------------------------------------------------
    // Sep 25, 2008: determine sigma_mu^int the same way as SALT2

    mindif = 999999.;
    for ( i = 0; i< 20; i++ ) {
      sigmu_tmp = (double)i * .01 ;
      sqmusig_tmp = sigmu_tmp * sigmu_tmp ;
      invert_mucovar(sqmusig_tmp);
      get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, sqmusig_tmp,  // inputs
	  	  &muoff_tmp, &snchi_tmp, &chi2tmp );   // return args

      dif = chi2tmp/(double)Ndof - 1.0 ;
      if ( fabsf(dif) < mindif ) {
        sigmu_int1 = sigmu_tmp ; mindif = dif;
      }
    }

    // search again in .001 sigmu_int bins

    mindif = 9999999. ;
    for ( i = -10; i < 10; i++ ) {
      sigmu_tmp = sigmu_int1 + (double)i * .0002 ;
      sqmusig_tmp = sigmu_tmp * sigmu_tmp ;
      invert_mucovar(sqmusig_tmp);

      get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, sqmusig_tmp,  // inputs
		  &muoff_tmp, &snchi_tmp, &chi2tmp );   // return args

      dif = chi2tmp/(double)Ndof - 1.0 ;

      if ( fabsf(dif) < mindif ) {
        sigmu_int = sigmu_tmp ; mindif = dif;
      }
    }


    invert_mucovar(sigmu_int);
    get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, sigmu_int*sigmu_int,  // inputs
	  	&muoff_tmp, &snchi_tmp, &chi2tmp );   // return args

    //  printf(" xxxx chi2(sigmu_int=%7.4f) = %8.2f  (Ndof=%d) \n", 
    //	 sigmu_int, chi2tmp, Ndof );


    /*************************************************************/
    /*** Output estimated cosmological parameters, obscured by ***/
    /*** sine of a large integer if blindness is enabled.      ***/
    /*************************************************************/
    if(usewa){

      if (blind){
      //    srand(time(NULL));
	srand(48901);
	w0rand   = floor(1e6*rand()/RAND_MAX);
	warand   = floor(1e6*rand()/RAND_MAX);
	ommrand = floor(1e6*rand()/RAND_MAX);
      } else {
	w0rand   = 0.0 ;
	warand   = 0.0 ;
	ommrand  = 0.0 ;
      }

      if ( strlen(cosparfilevar) == 0){
	strcpy(cosparfilevar,infile);
	strcat(cosparfilevar,".cospar");
      } 
    
      getname(cosparfilevar, tempfilename1, fitnumber);
      strcpy(cosparfile, tempfilename1);
      printf("Write cosmo params to %s \n",cosparfile);
      fpcospar=fopen(cosparfile, "w");
      if (fpcospar == NULL){
	printf("ERROR: couldn't open %s\n",cosparfile);
	exit(EXIT_ERRCODE_wfit);
      }

    // modified format to sneak in chi2 & Ndof (RK May 8, 2008)
    // modified format to include w0_marg sigma (JLM Sep 26, 2013)
    // modified format to include label (RK, Jun 2017)

      w0_out   += sin(w0rand);
      wa_out   += sin(warand);
      omm_out += sin(ommrand);
      

    // Jun 2019: check for csv format
      char sep[] = " ";
      if ( csv_out ) { sprintf(sep,","); }
      RESULTS.w0_out       = w0_out ;
      RESULTS.w0sig        = w0sig  ;
      RESULTS.w0sig_upper  = w0sig_upper  ;
      RESULTS.w0sig_lower  = w0sig_lower  ;
      RESULTS.wa_out       = wa_out ;
      RESULTS.wasig        = wasig  ;
      RESULTS.wasig_upper  = wasig_upper  ;
      RESULTS.wasig_lower  = wasig_lower  ;
      RESULTS.omm_out     = omm_out;
      RESULTS.omm_sig     = omm_sig;
      RESULTS.chi2_final  = chi2_final;
      RESULTS.Ndof        = Ndof ;
      RESULTS.sigmu_int   = sigmu_int ;
      RESULTS.w0rand       = w0rand;
      RESULTS.warand       = warand;
      RESULTS.ommrand     = ommrand ;
      sprintf(RESULTS.label_cospar, "%s", label_cospar);
      write_output_cospar(fpcospar, &RESULTS, usemarg, format_cospar);
      fclose(fpcospar);
    }

    else{
      if (blind){
      //    srand(time(NULL));
	srand(48901);
	w0rand   = floor(1e6*rand()/RAND_MAX);
	ommrand = floor(1e6*rand()/RAND_MAX);
      } else {
	w0rand   = 0.0 ;
	ommrand = 0.0 ;
      }

      if ( strlen(cosparfilevar) == 0){
	strcpy(cosparfilevar,infile);
	strcat(cosparfilevar,".cospar");
      } 
    
      getname(cosparfilevar, tempfilename1, fitnumber);
      strcpy(cosparfile, tempfilename1);
      printf("Write cosmo params to %s \n",cosparfile);
      fpcospar=fopen(cosparfile, "w");
      if (fpcospar == NULL){
	printf("ERROR: couldn't open %s\n",cosparfile);
	exit(EXIT_ERRCODE_wfit);
      }

    // modified format to sneak in chi2 & Ndof (RK May 8, 2008)
    // modified format to include w0_marg sigma (JLM Sep 26, 2013)
    // modified format to include label (RK, Jun 2017)

      w0_out   += sin(w0rand);
      omm_out += sin(ommrand);


    // Jun 2019: check for csv format
      char sep[] = " ";
      if ( csv_out ) { sprintf(sep,","); }
      RESULTS.w0_out       = w0_out ;
      RESULTS.w0sig        = w0sig  ;
      RESULTS.w0sig_upper  = w0sig_upper  ;
      RESULTS.w0sig_lower  = w0sig_lower  ;
      RESULTS.omm_out     = omm_out;
      RESULTS.omm_sig     = omm_sig;
      RESULTS.chi2_final  = chi2_final;
      RESULTS.Ndof        = Ndof ;
      RESULTS.sigmu_int   = sigmu_int ;
      RESULTS.w0rand       = w0rand;
      RESULTS.ommrand     = ommrand ;
      sprintf(RESULTS.label_cospar, "%s", label_cospar);
      write_output_cospar(fpcospar, &RESULTS, usemarg, format_cospar);
      fclose(fpcospar);
    }
    /*********************************************/
    /*** Write out the prob dist using cfitsio ***/
    /*********************************************/
	  
    if (fitsflag){
    
      naxes[0] = omm_steps;
      naxes[1] = w0_steps;
      status = 0;
      nelements = naxes[0] * naxes[1];    /* number of pixels to write */
      fpixel[0]=1;			/* first pixel to write */
      fpixel[1]=1;
      
      if (writechi) {
        /* subtract off minchi from each, use prob arrays for storage */
        for(i=0; i<w0_steps; i++){
	  for(kk=0; kk<wa_steps;kk++){
	    for(j=0; j<omm_steps; j++){
	      //	      snprob[i*omm_steps+j]  = snchi[i*omm_steps+j] - snchi_min;
	      //	      extprob[i*omm_steps+j] = extchi[i*omm_steps+j] - extchi_min;
              snprob3d[j][kk][i]  = snchi3d[j][kk][i] - snchi_min;
              extprob3d[j][kk][i] = extchi3d[j][kk][i] - extchi_min;
	    }
	  }
        }
      
        strcpy(snfitsname,infile);
        strcat(snfitsname,"_snchi.fits");

	getname(snfitsname, tempfilename1, fitnumber);

        strcpy(ommfitsname,infile);
        strcat(ommfitsname,"_extchi.fits");

	getname(ommfitsname, tempfilename2, fitnumber);

        printf("writing out chi2 distributions:\n");
        printf("  SNe only:      %s\n",tempfilename1);
        printf("  SNe + Omega_m: %s\n",tempfilename2);

      } else {
        /* just write out probabilities */
        strcpy(snfitsname,infile);
        strcat(snfitsname,"_snprob.fits");

	getname(snfitsname, tempfilename1, fitnumber);

        strcpy(ommfitsname,infile);
        strcat(ommfitsname,"_extprob.fits");
    
	getname(ommfitsname, tempfilename2, fitnumber);

        printf("writing out prob distributions:\n");
        printf("  SNe only:      %s\n",tempfilename1);
        printf("  SNe + Omega_m: %s\n",tempfilename2);
      }

      /* prepend with "!" so cfitsio clobbers old files */
      sprintf(snfitsstr,"!%s",tempfilename1);
      sprintf(ommfitsstr,"!%s",tempfilename2);
    
      /* Write out SN-only distribution */
      if (fits_create_file(&fptr, snfitsstr, &status))  printerror( status );
      if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) 
        printerror( status );          
      if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, snprob, &status))
        printerror( status );
      if ( fits_close_file(fptr, &status))                /* close the file */
        printerror( status );           
    
      /* Write out SN+Omega_m distribution */
      if (fits_create_file(&fptr, ommfitsstr, &status))  printerror( status );
      if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status))
        printerror( status );          
      if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, extprob, &status))
        printerror( status );
      if ( fits_close_file(fptr, &status))                /* close the file */
        printerror( status );           

      // 6/10/2008 RK write chi2 to text file (like Andy's)

      if (chi2gridfilevar[0] == 0){
        strcpy(chi2gridfilevar, infile);
      }
      strcpy(chi2gridfile, chi2gridfilevar); 
      strcat(chi2gridfile,".chi2grid");
      sprintf(chi2gridfile_SN,"%s_SN", chi2gridfile);
      

      getname(chi2gridfile, tempfilename1, fitnumber);
      getname(chi2gridfile_SN, tempfilename2, fitnumber);
      
      printf("\n");
      printf("\t CHI2GRID(total)   dump to: '%s' \n", tempfilename1 );
      printf("\t CHI2GRID(SN-only) dump to: '%s' \n", tempfilename2 );
      printf("\t CHI2GRID format : i_OM  i_w  OM  w  chi2 \n" );

      int irow=0; char CVARDEF[100];
      FILEPTR_TOT = fopen(tempfilename1, "wt");
      FILEPTR_SN  = fopen(tempfilename2, "wt");

      sprintf(CVARDEF,"NVAR: 6 \nVARNAMES: ROW iM iw OM w chi2\n");
      fprintf(FILEPTR_TOT,"%s", CVARDEF);
      fprintf(FILEPTR_SN, "%s", CVARDEF);

        for(i=0; i<w0_steps; i++){
	  for(j=0; j<omm_steps; j++){
	    w0   = w0_min   + i*w0_stepsize;
	    omm = omm_min + j*omm_stepsize;
	    chi2tot = extchi[i*omm_steps+j] ;
	    chi2sn  = snchi[i*omm_steps+j] ;

	    irow++ ;
	    fprintf( FILEPTR_TOT, "ROW: %5d  %5d %5d  %8.4f  %8.4f  %14.6f \n",
		     irow, j, i, omm, w0, chi2tot  );

	    fprintf( FILEPTR_SN, "ROW: %5d  %5d %5d  %8.4f  %8.4f  %14.6f \n",
		     irow, j, i, omm, w0, chi2sn  );

	  }
        }
      fclose(FILEPTR_TOT);
      fclose(FILEPTR_SN);
    }
    if ( fitnumber == 0){
      if (sigmu_int > snrms ) {
	printf("REFITTING DATA WITH UPDATED INTRINSIC SCATTER\n\n");
	printf("\t input snrms was %0.4f \n", snrms);
	printf("\t refitting with snrms %0.4f \n", sigmu_int);
	snrms = sigmu_int;
	sqsnrms = snrms*snrms;
	fitcount = fitcount+1;
      } else {
	printf("SKIPPING REFIT - NO ADDITIONAL INTRINSIC SCATTER\n\n");
	printf("\t initial snrms was %0.4f \n", snrms);
	printf("\t final int scatter was %0.4f \n\n", sigmu_int);
	printf("\t moving inital files to final files.\n");

	getname(residfilevar, tempfilename1, 0);
	getname(residfilevar, tempfilename2, 1);
	printf("\t\t moving %s to %s \n", tempfilename1, tempfilename2);
	rename(tempfilename1, tempfilename2);
	
	getname(cosparfilevar, tempfilename1, 0);
	getname(cosparfilevar, tempfilename2, 1);
	printf("\t\t moving %s to %s \n", tempfilename1, tempfilename2);
	rename(tempfilename1, tempfilename2);

	getname(snfitsname, tempfilename1, 0);	
	getname(snfitsname, tempfilename2, 1);	
	printf("\t\t moving %s to %s \n", tempfilename1, tempfilename2);
	rename(tempfilename1, tempfilename2);

	getname(ommfitsname, tempfilename1, 0);
	getname(ommfitsname, tempfilename2, 1);
	printf("\t\t moving %s to %s \n", tempfilename1, tempfilename2);
	rename(tempfilename1, tempfilename2);

	getname(chi2gridfile, tempfilename1, 0);
	getname(chi2gridfile, tempfilename2, 1);
	printf("\t\t moving %s to %s \n", tempfilename1, tempfilename2);
	rename(tempfilename1, tempfilename2);

	getname(chi2gridfile_SN, tempfilename1, 0);
	getname(chi2gridfile_SN, tempfilename2, 1);
	printf("\t\t moving %s to %s \n", tempfilename1, tempfilename2);
	rename(tempfilename1, tempfilename2);

	fitnumber = fitnumber+1;
      }
    }
    fitnumber = fitnumber+1;
  } /* end refit while loop */

  /**************************************/
  /*** Free up the memory and go home ***/
  /**************************************/

  free(mu);
  free(mu_sig);
  free(mu_ref);
  free(z);
  free(z_sig);
  free(snnbin);
  free(chitmp);
  free(tid);
  
  free(w0_prob);
  free(w0_sort);
  free(wa_prob);
  free(wa_sort);
  free(hvec); 
  free(omm_prob);

  free(extprob3d);
  free(snprob3d);
  free(extchi3d);
  free(snchi3d);


  printf("DONE. \n");
  return(0);
}  


// ==================================
void read_fitres(char *inFile) {

  /*************
    Created Jun 13, 2008 by R.Kessler

    Read snana fitres file. Read either old-style (LEGACY)
    format, or the newer self-documented format that uses
    "NVAR:", "VARNAMES:" and "SN:" keywords.

    With the newer format, we do NOT need the -snana or -snanasim
    flags. Only need the -snanasim flag for the LEGACY format.

   Aug 20, 2009: allow Z[ERR] or ZPHOT[ERR]

   June 30 2016: float -> double 

   Jun 21 2017: init IWD_NBIN=-9 to fix bug
   Jun 25 2017: 
     + optional -label <label> arg to label the .cospar output
     + check MUERR_INCLUDE_xxx flags
     + skip events with MUERR>100 (unfilled bins from SALT2mu M0DIF output)

   Apr 1 2019: fix to work without NVAR key; beware that code is fragile.
   Apr 3 2019: remove LEGACY check; only reads FITRES-formatted files.

   Sep 15 2020: use snana_openTextFile to handle gzipped inputs

  *************/

  char ctmp[80] ,ctmp2[80], VARLIST[200] ;
  int NVAR, IWD ;
  int IWD_CID, IWD_Z, IWD_ZERR, IWD_MU, IWD_MUERR, IWD_MUREF=-9, IWD_TID ;
  int IWD_MUDIF, IWD_MUDIFERR ;
  int IWD_GTYPE, IWD_SNTYPE, IWD_NBIN, IWD_NFIT ;

  char CID[12], inFile_opened[200];
  double Z, ZERR, MU, MUERR, MUREF ;
  int TID, GTYPE, STYPE, LCUT, i, NRDTOT, NBIN, NFIT, ISROWKEY, gzipFlag ;

  FILE *fp;
  char fnam[] = "read_fitres";

  // -------- BEGIN --------

  printf(" Open fitres file: %s \n", inFile);

  // xxx mark delete fp = fopen(inFile,"rt");
  int OPENMASK = OPENMASK_VERBOSE ;
  fp = snana_openTextFile(OPENMASK, "", inFile, 
			  inFile_opened, &gzipFlag );

  NCIDLIST = NCUT = NRDTOT = NSNE_NBIN = NBIN = 0;

  /* xxxx mark delete 
  printf(" Reading fitres file: NCIDLIST %d, NSNE_NBIN %d \n", 
	 NCIDLIST, NSNE_NBIN);
  */

  MUERR_INCLUDE_zERR=0;
  MUERR_INCLUDE_LENS=0;

  NVAR = 0;
  IWD_CID   =  0 ;  // CID is always the first word
  IWD_Z = IWD_ZERR = IWD_MU = IWD_MUERR = IWD_TID = -9;
  IWD_MUDIF = IWD_MUDIFERR = -9 ;
  IWD_GTYPE = IWD_SNTYPE = IWD_NBIN = IWD_NFIT = -9 ;

  IWD = 0;


  // read header info and stop once VARNAMES key is reached.
  while ( IWD < 200 ) {

    IWD++ ;

    readchar ( fp, ctmp) ; 

    if ( strcmp(ctmp,"MUERR_INCLUDE:") == 0 )  {
      readchar(fp,ctmp2);
      if ( strcmp(ctmp2,"zERR") == 0 ) 
	{ MUERR_INCLUDE_zERR = 1; }
      else if ( strcmp(ctmp2,"ZERR") == 0 ) 
	{ MUERR_INCLUDE_zERR = 1; }
      else if ( strcmp(ctmp2,"zerr") == 0 ) 
	{ MUERR_INCLUDE_zERR = 1; }
      else if (strstr(ctmp2,"SIGMA_LENS") != NULL ) 
	{ MUERR_INCLUDE_LENS = 1; }
    } 

    if ( strcmp(ctmp,"VARNAMES:") == 0 ) { goto READ_VARLIST;  }

  } //end of while (IWD < 200) loop

  sprintf(c1err,"Could not find required VARNAMES key");
  sprintf(c2err,"Check %s", inFile);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);


 READ_VARLIST:

  fgets( VARLIST, 100, fp );
  NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,VARLIST);

  for(IWD = 0 ; IWD < NVAR; IWD++ ) {

    get_PARSE_WORD(0, IWD, ctmp);

    if ( strcmp(ctmp,"CID")      == 0 ) { IWD_CID = IWD ; }
    if ( strcmp(ctmp,"cid")      == 0 ) { IWD_CID = IWD ; }
    if ( strcmp(ctmp,"ROW")      == 0 ) { IWD_CID = IWD ; }
    
    // note that zHD takes priority over z (Aug 2015)
    if ( strcmp(ctmp,"zHD")        == 0 ) IWD_Z     = IWD;
    if ( strcmp(ctmp,"zHDERR")     == 0 ) IWD_ZERR  = IWD;

    if ( IWD_Z < 0 ) {
      if ( strcmp(ctmp,"Z")        == 0 ) IWD_Z     = IWD;
      if ( strcmp(ctmp,"z")        == 0 ) IWD_Z     = IWD;
    }
    if ( IWD_ZERR < 0 ) {
      if ( strcmp(ctmp,"ZERR")     == 0 ) IWD_ZERR  = IWD;	
      if ( strcmp(ctmp,"zERR")     == 0 ) IWD_ZERR  = IWD;
    }

    // if there is a photoz, then use it.
    if ( strcmp(ctmp,"ZPHOT")    == 0 ) IWD_Z     = IWD;
    if ( strcmp(ctmp,"ZPHOTERR") == 0 ) IWD_ZERR  = IWD;

    // allow DLMAG or MU
    if ( strcmp(ctmp,"MU")       == 0 ) IWD_MU    = IWD;
    if ( strcmp(ctmp,"MUERR")    == 0 ) IWD_MUERR = IWD;
    if ( strcmp(ctmp,"DLMAG")    == 0 ) IWD_MU    = IWD;
    if ( strcmp(ctmp,"DLMAGERR") == 0 ) IWD_MUERR = IWD;
    

    if ( strcmp(ctmp,"MUDIF")    == 0 ) {
      IWD_MUDIF    = IWD;
      printf("\n !!! Found MUDIF --> fit for w-wref and OM-OMref \n\n");
      ZERR       = 0.0 ;
      mudif_flag = 1 ;
      MUERR_INCLUDE_zERR = 1;
    }
    if ( strcmp(ctmp,"MUDIFERR") == 0 ) IWD_MUDIFERR = IWD;
    if ( strcmp(ctmp,"MUREF")    == 0 ) IWD_MUREF    = IWD;
    
    if ( strcmp(ctmp,"IDTEL")    == 0 ) IWD_TID   = IWD;
    if ( strcmp(ctmp,"IDSURVEY") == 0 ) IWD_TID   = IWD;

    if ( strcmp(ctmp,"GENTYPE") == 0 ) IWD_GTYPE   = IWD;
    if ( strcmp(ctmp,"SNTYPE")  == 0 ) IWD_SNTYPE  = IWD;

    if ( strcmp(ctmp,"NBIN")  == 0 )   IWD_NBIN  = IWD;

    if ( strcmp(ctmp,"NFIT")  == 0 )   IWD_NFIT  = IWD;
    
  } // end of IWD-varname read loop


  int ERR = EXIT_ERRCODE_wfit;
  if ( IWD_Z < 0 ) 
    { printf(" ABORT ON ERROR: IWD_Z = %d \n", IWD_Z ); exit(ERR); }

  if ( IWD_ZERR < 0 && mudif_flag == 0 ) 
    { printf(" ABORT ON ERROR: IWD_ZERR = %d \n", IWD_ZERR ); exit(ERR); }

  if ( IWD_MU < 0  && IWD_MUDIF<0 ) 
    { printf(" ABORT ON ERROR: IWD_MU = %d \n", IWD_MU ); exit(ERR); }

  if ( IWD_MUERR < 0 && IWD_MUDIFERR < 0 ) 
    { printf(" ABORT ON ERROR: IWD_MUERR = %d \n", IWD_MUERR ); exit(ERR); }

  if ( IWD_TID < 0 && mudif_flag==0 ) 
    { printf(" ABORT ON ERROR: IWD_TID = %d \n", IWD_TID ); exit(ERR); }

  //  printf(" xxx IWD(z,MUDIF,MUDIFERR) = %d   %d %d  \n",
  //	 IWD_Z, IWD_MUDIF, IWD_MUDIFERR );

  IWD = 0;
  
  

  while( (fscanf(fp, "%s", ctmp)) != EOF) {

    ISROWKEY = ( strcmp(ctmp,"SN:")==0  || strcmp(ctmp,"ROW:")==0 );

    if ( ISROWKEY )  { 
      IWD = 0; 
      MU = MUERR = MUREF = Z = ZERR = 0.0 ;
      continue ; 
    }

    if ( IWD == IWD_CID      ) sscanf ( ctmp, "%s",  CID    );
    if ( IWD == IWD_Z        ) sscanf ( ctmp, "%le", &Z     );
    if ( IWD == IWD_ZERR     ) sscanf ( ctmp, "%le", &ZERR  );
    if ( IWD == IWD_MU       ) sscanf ( ctmp, "%le", &MU    );
    if ( IWD == IWD_MUERR    ) sscanf ( ctmp, "%le", &MUERR );

    if ( IWD == IWD_MUDIF    ) sscanf ( ctmp, "%le", &MU    );
    if ( IWD == IWD_MUDIFERR ) sscanf ( ctmp, "%le", &MUERR );
    if ( IWD == IWD_MUREF    ) sscanf ( ctmp, "%le", &MUREF );

    if ( IWD == IWD_TID      ) sscanf ( ctmp, "%i", &TID   );

    if ( IWD == IWD_GTYPE)   sscanf ( ctmp, "%i", &GTYPE );
    if ( IWD == IWD_SNTYPE)  sscanf ( ctmp, "%i", &STYPE );
    if ( IWD == IWD_NBIN)    sscanf ( ctmp, "%i", &NBIN );
    if ( IWD == IWD_NFIT)    sscanf ( ctmp, "%i", &NFIT );

    IWD++ ; 

    // after reading all NVAR variables, increment "NCIDLIST"
    // and store variables needed for cosmology fit

    if ( IWD == NVAR ) {

      // check cuts
      LCUT = 1;
      if ( Z < CUTSELECT.ZMIN ) LCUT = 0 ;
      if ( Z > CUTSELECT.ZMAX ) LCUT = 0 ;

      if ( CUTSELECT.GENTYPE_SIM > 0 && 
	   GTYPE != CUTSELECT.GENTYPE_SIM ) { LCUT = 0; }

      if ( CUTSELECT.SNTYPE > 0 && 
	   STYPE != CUTSELECT.SNTYPE ) { LCUT = 0; }

      
      if ( NCIDLIST >= CUTSELECT.MXSNFIT ) { LCUT = 0 ; }

      if ( MUERR > 100.0 ) { LCUT = 0; } // Jun 26 2017

      if ( mudif_flag && NFIT <= 1 ) { LCUT=0; } // Sep 21 2020


      // sanity checks
      if ( LCUT ) {
	if ( (Z < 0.0 && fabs(Z+9.0)>.001 )  || Z > 5.0 ) 
	  { printf(" Found INSANE Z = %f  ==> ABORT \n", Z ); exit(ERR); }

	if ( ZERR < 0.0 || ZERR > 10.0 ) 
	  { printf(" Found INSANE ZERR = %f  ==> ABORT \n",ZERR ); exit(ERR);}

	if ( mudif_flag == 0 ) {
	  if ( MU < 0.0 || MU > 100.0 ) 
	    { printf(" Found INSANE MU = %f  ==> ABORT \n", MU ); exit(ERR); }
	}
	else {
	  if ( fabs(MU) > 2. && fabs(MU/MUERR) > 3.0 ) 
	    { printf(" Found INSANE MUDIF = %f += %f  ==> ABORT \n",
		     MU, MUERR ); exit(ERR); }
	}
      } // end LCUT 

      // xxx move check cuts above (9.21.2020) xxxx

      NRDTOT++ ;

      if ( LCUT  ) {
	NCIDLIST++ ;  i = NCIDLIST-1;
	sprintf(CIDLIST[i],"%s", CID);
	z[i]       = Z;
	z_sig[i]   = ZERR;
	mu[i]      = MU ;
	mu_sig[i]  = MUERR;
	mu_ref[i]  = MUREF ; // 5.2019
	tid[i]     = TID ;
	if ( NBIN ) {
	  NSNE_NBIN = NSNE_NBIN + NBIN; 
	  snnbin[i] = NBIN;
	} else {
	  snnbin[i] = 1;
	}

      }
      else {
	NCUT++;
      }

      /***
      printf(" xxx NCIDLIST=%3d  CID=%5s Z=%6.4f +- %6.4f  MU = %7.3f +- %6.3f   TID = %d \n"
	     ,NCIDLIST, CID
	     ,z[i], z_sig[i]
	     , mu[i], mu_sig[i]
	     , tid[i] 
	     );
      *****/

    }

  }

  fclose(fp);

  printf(" Read %d SNe from %s \n", NRDTOT, inFile);
  printf(" Select %d SNe for fitting.\n", NCIDLIST );
  printf(" MUERR_INCLUDE(zERR,LENS) = %d, %d \n",
	 MUERR_INCLUDE_zERR, MUERR_INCLUDE_LENS );
  printf(" Done reading file -- NCIDLIST: %d, NSNE_NBIN: %d \n", 
	 NCIDLIST, NSNE_NBIN);

  fflush(stdout);

  return ;

} // end of read_fitres



// ==================================
void read_mucovar(char *inFile) {

  /*************
    Created May 20, 2009 by R.Kessler

    Read & load off-diagonal distance-modulus covariances.

  *************/

  char ctmp[80], SN[2][12], locFile[1000] ;
  char fnam[] = "read_mucovar" ;   
  double cov;
  int N, N0, N1, i0, i1, j;
  FILE *fp;

  // -------- BEGIN --------


  // first try to open file in current user area

  sprintf(locFile, "%s", inFile );
  fp = fopen(locFile,"rt");
  if ( fp != NULL ) goto READIT ;

  // now try reading from official (public) area

  sprintf(locFile, "%s/models/mucovar/%s", 
	  getenv("SNDATA_ROOT"), inFile );
  fp = fopen(locFile,"rt");
  if ( fp == NULL ) {
    printf(" Could not open mucovar file: %s\n", inFile );
    printf(" ***** ABORT ****** \n");
    exit(EXIT_ERRCODE_wfit);
  }


 READIT:

  // init stuff
  for ( N0=0; N0 < MXCOVSN; N0++ ) {
    for ( N1=0; N1 < MXCOVSN; N1++ ) {
      MUCOV_STORE[N0][N1].COV      = -9.0 ;
      MUCOV_STORE[N0][N1].INDEX[0] = -9   ;
      MUCOV_STORE[N0][N1].INDEX[1] = -9   ;
    }
  }


  printf("\n Open mucovar file: \n  %s \n", locFile);

  while( (fscanf(fp, "%s", ctmp)) != EOF) {

    if ( strcmp(ctmp,"COV:") == 0  ) {

      readchar(fp, SN[0] );
      readchar(fp, SN[1] );
      readdouble(fp, 1, &cov );

      i0 = cidindex(SN[0]); // global sparse index for fitres input
      i1 = cidindex(SN[1]);

      // store off-diag elements only
      if ( i0 >= 0 && i1 >= 0 && i0 != i1 ) {

	NCOVPAIR++;  N = NCOVPAIR;

	if ( N >= MXCOVPAIR ) {
	  printf("\n %s FATAL ERROR: NCOVPAIR=%d exceeds bound.\n", 
		 fnam, NCOVPAIR);
	  printf(" ***** ABORT ***** \n");
	  exit(EXIT_ERRCODE_wfit);
	}

	sprintf ( MUCOV_INPUT[N].CID[0], "%s", SN[0] ) ;
	sprintf ( MUCOV_INPUT[N].CID[1], "%s", SN[1] ) ;
	MUCOV_INPUT[N].COV      = cov ;
	MUCOV_INPUT[N].INDEX[0] = i0 ;
	MUCOV_INPUT[N].INDEX[1] = i1 ;

	if ( INDEX_COVSN_INVMAP[i0] < 0 ) {
	  NCOVSN++ ; 
	  INDEX_COVSN_MAP[NCOVSN]  = i0 ; 
	  INDEX_COVSN_INVMAP[i0]   = NCOVSN ; 
	}

	if ( INDEX_COVSN_INVMAP[i1] < 0 ) {
	  NCOVSN++ ; 
	  INDEX_COVSN_MAP[NCOVSN]  = i1 ; 
	  INDEX_COVSN_INVMAP[i1]   = NCOVSN ; 
	}

	// get sparse indices for MUCOV_STORE array
	N0 = N1 = 0;
	for ( j=1; j <= NCOVSN; j++ ) {
	  if ( INDEX_COVSN_MAP[j] == i0 ) N0 = j;
	  if ( INDEX_COVSN_MAP[j] == i1 ) N1 = j;
	}


	if ( N0 <= 0 || N1 <= 0 ) {
	  printf("\n %s FATAL ERROR: N0=%d  N1=%d \n", fnam, N0, N1);
	  printf(" ***** ABORT \n");
	  exit(EXIT_ERRCODE_wfit);
	}

	MUCOV_STORE[N0][N1].COV      = cov;
	MUCOV_STORE[N0][N1].INDEX[0] = i0 ;
	MUCOV_STORE[N0][N1].INDEX[1] = i1 ;

	MUCOV_STORE[N1][N0].COV      = cov;
	MUCOV_STORE[N1][N0].INDEX[0] = i0 ;
	MUCOV_STORE[N1][N0].INDEX[1] = i1 ;

      }

      /*
      printf(" xxx Read MUCOVAR for SN pair: %s and %s (indices %d %d) \n"
	     ,MUCOV_INPUT[N].CID[0],   MUCOV_INPUT[N].CID[1]
	     ,MUCOV_INPUT[N].INDEX[0], MUCOV_INPUT[N].INDEX[1]
	     );
      */

    }
  } // end of read loop


  printf(" Store %d off-diagonal MU-covariances for %d SNe. \n", 
	 NCOVPAIR, NCOVSN );


  // check for missing cov elements

  for ( N0=1; N0 <= NCOVSN; N0++ ) {
    for ( N1=1; N1 <= NCOVSN; N1++ ) {

      if ( N0 == N1 ) continue ;

      i0  = INDEX_COVSN_MAP[N0];
      i1  = INDEX_COVSN_MAP[N1];
      cov = MUCOV_STORE[N0][N1].COV ;

      if ( cov < -8. ) {

	if ( N1 > N0 ) 
	  printf("\t WARNING: cov(%s,%s) is undefined. Set to zero. \n",
		 CIDLIST[i0], CIDLIST[i1] );

	MUCOV_STORE[N0][N1].COV      = 0.0 ;
	MUCOV_STORE[N0][N1].INDEX[0] = i0  ;
	MUCOV_STORE[N0][N1].INDEX[1] = i1  ;
      }

    }  // N1
  } // N0

  fclose(fp);



} // end of read_mucovar()


//===================================
void set_priors(void) {

  //  char fnam[]="set_priors";

  double rz, tmp1, tmp2;
  double OM = OMEGA_MATTER_SIM ;
  double OE = 1 - OM ;
  double w0 = w0_DEFAULT ;
  double wa= wa_DEFAULT ;
  Cosparam cparloc;

  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa; //0.0 ;

  // ===== BEGIN ==========
  
  if ( usebao == 2 ){ 
    //recompute abest

    rz = codist(z_bao, &cparloc);
    tmp1 = pow( EofZ(z_bao, &cparloc), NEGTHIRD) ;
    tmp2 = pow( (1./z_bao) * rz, TWOTHIRD );
    abest = sqrt(OM) * tmp1 * tmp2 ;

  }

  if ( usecmb == 2) {
    //recompute Rbest
    HzFUN_INFO_DEF HzFUN;
    set_HzFUN_for_wfit(ONE, OM, OE, w0, wa, &HzFUN);
    // xxx mark delete rz=Hainv_integral(ONE,OM,OE,w0,acmb, ONE ) / LIGHT_km;
    rz = Hainv_integral ( acmb, ONE, &HzFUN ) / LIGHT_km;
    Rcmb_best = sqrt(OM) * rz ;
  }

  return;

} //end of set_priors()

// ==================================
void set_HzFUN_for_wfit(double H0, double OM, double OE, double w0, double wa,
			HzFUN_INFO_DEF *HzFUN ) {

  double COSPAR_LIST[10];  int VBOSE=0;
  // ---------- BEGIN ----------

  COSPAR_LIST[ICOSPAR_HzFUN_H0] = 1.0;
  COSPAR_LIST[ICOSPAR_HzFUN_OM] = OM ;
  COSPAR_LIST[ICOSPAR_HzFUN_OL] = OE ;
  COSPAR_LIST[ICOSPAR_HzFUN_w0] = w0 ;
  COSPAR_LIST[ICOSPAR_HzFUN_wa] = wa; //0.0 ;
  init_HzFUN_INFO(VBOSE, COSPAR_LIST, "", HzFUN);
    
} // end set_HzFUN_for_wfit

// ==================================
void invert_mucovar(double sqmurms_add) {

  // May 20, 2009
  // add diagonal elements to covariance matrix, then invert.
  //
  // Feb 2013: replace CERNLIB's dfact,dfinv with invertMatrix based on gsl.
  //
  int N0, N1, i0, NMAX;
  double covtmp[MXCOVSN][MXCOVSN];

  // =================================

  if ( NCOVSN <= 0 ) return ;

  // set diagonal terms to errors computed before this function was called.
  for ( N0=1; N0 <= NCOVSN; N0++ ) {
    i0  = INDEX_COVSN_MAP[N0];
    MUCOV_STORE[N0][N0].INDEX[0] = i0 ;
    MUCOV_STORE[N0][N0].INDEX[1] = i0 ;
    MUCOV_STORE[N0][N0].COV      = mu_sqsig[i0] + sqmurms_add ;
  }


  // load full COV matrix into temp array that starts at zero.

  for ( N0=1; N0 <= NCOVSN; N0++ ) {
    for ( N1=1; N1 <= NCOVSN; N1++ ) {
      covtmp[N0-1][N1-1] = MUCOV_STORE[N0][N1].COV ;
    }
  }

  NMAX = MXCOVSN ;

  invertMatrix( NMAX, NCOVSN, &covtmp[0][0] ) ;

  // store inverted matrix in global struct

  for ( N0=1; N0 <= NCOVSN; N0++ ) {
    for ( N1=1; N1 <= NCOVSN; N1++ ) {
      MUCOV_STORE[N0][N1].COVINV = covtmp[N0-1][N1-1];
    }
  }


} // end of invert_mucovar



// ===========================
void get_chi2wOM ( 
	      double w0             // (I) 
	      ,double wa           // (I)
	      ,double OM           // (I) 
	      ,double sqmurms_add  // (I) anomalous mu-error squared
	      ,double *mu_off      // (O) distance offset
	      ,double *chi2sn      // (O) SN-only chi2
	      ,double *chi2tot     // (O) SN+prior chi2
	      ) {

  // Jul 10, 2021: Return chi2 at this w0, wa, OM and z value
  // Jan 23, 2009: return chi2 at this w,OM value
  // Apr 23, 2009: fix bug: use cparloc instead of global cpar
  // May 22, 2009: add args *mu_off and sqmurms_add
  //               include nonflat H0 prior term using SQSIG_MUOFF
  //

  double OE ;

  double 
    a, rz, sqmusig, sqmusiginv, Bsum, Csum
    ,chi_hat, dchi_hat, ld_cos
    ,mu_cos[MXSN]
    ,tmp1, tmp2
    ,Rcmb
    ,nsig
    ,dmu
    ,sqdmu
    ,covinv
    ;

  Cosparam cparloc;
  int k, k0, k1, N0, N1, icov;

  // --------- BEGIN --------

  OE = 1.0 - OM ;
  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa ;


  Bsum = Csum = chi_hat = 0.0 ;

  /* Loop over all data and calculate chi2 */
  for (k=0; k < NCIDLIST; k++){

    sqmusig     = mu_sqsig[k] + sqmurms_add ;
    sqmusiginv  = 1./sqmusig ;
    
    rz        = codist(z[k], &cparloc); // cparloc -- input. w, om etc. 
    ld_cos    = (1+z[k]) *  rz * c_light / H0;
    mu_cos[k] =  5.*log10(ld_cos) + 25. ;


    dmu  = mu_cos[k] - mu[k] ;

    // add chi2 only for SNe that have no off-diag terms

    icov = INDEX_COVSN_INVMAP[k] ; 

    if ( icov < 0 )  {
      Bsum    += sqmusiginv * dmu ;       // Eq. A.11 of Goliath 2001
      Csum    += sqmusiginv ;             // Eq. A.12 of Goliath 2001
      chi_hat += sqmusiginv * dmu*dmu ;
    }

    /*
    printf("  xxxx dmu(%4s) = %6.2f - %6.2f = %6.2f   sigmu=%f  \n", 
	   CIDLIST[k], mu_cos[k], mu[k], dmu, mu_sig[k] );
    */
  }


  // check for SNe with off-diagonal terms
  // Note that below includes both diag & off-diag for 
  // the SN-subset with covariances.

  if ( NCOVPAIR > 0 ) {
    for ( N0=1; N0 <= NCOVSN; N0++ ) {
      for ( N1=1; N1 <= NCOVSN; N1++ ) {

	k0     = INDEX_COVSN_MAP[N0]  ;
	k1     = INDEX_COVSN_MAP[N1]  ;
	covinv = MUCOV_STORE[N0][N1].COVINV ;

	sqdmu    = (mu_cos[k0] - mu[k0]) * (mu_cos[k1] - mu[k1]) ;
	dchi_hat = sqdmu * covinv ;
	chi_hat += dchi_hat ;

	// increment B & C terms
	dmu    = mu_cos[k0] - mu[k0] ;
	Bsum  += covinv * dmu ;
	Csum  += covinv ;      

	if ( k0 >= 9999990 ) 
	  printf("  xxxx covini(%4s,%4s) = %10.4f   dchi=%12.4f  chi=%f \n", 
		 CIDLIST[k0], CIDLIST[k1], covinv, dchi_hat, chi_hat );

      } // N1
    } // N0
  }


  *mu_off  = Bsum/Csum ;  // load function output before adding H0-prior corr

  /* Analytic marginalization over H0.  
     See appendices in Goliath et al, A&A, 2001 ,
     and add contribution from PRIOR(H0) = Gaussian instead of flat.
  */
  
  Csum    += 1./SQSIG_MUOFF ;  // H0 prior term
  *chi2sn  =  chi_hat - Bsum*Bsum/Csum ;

  // Ignore constant term:
  //  logCsum  = log(Csum/(2.*M_PI));
  //  *chi2sn += logCsum;

  *chi2tot = *chi2sn ;     // load intermediate function output

  /* Compute likelihood with external prior: Omega_m or BAO */
  if (usebao) {
    /* Use BAO constraints from Eisenstein et al, 2006 */

    rz = codist(z_bao, &cparloc);
    tmp1 = pow( EofZ(z_bao, &cparloc), NEGTHIRD) ;
    tmp2 = pow( (1./z_bao) * rz, TWOTHIRD );
    a = sqrt(OM) * tmp1 * tmp2 ;
    nsig = (a-abest)/sigma_a ;
    *chi2tot =  *chi2sn + pow( nsig, TWO );

  } else {
    /* Use Omega_m constraints */
    nsig = (OM-omm_prior)/omm_prior_sig ;
    *chi2tot =  *chi2sn + pow( nsig, TWO );
  }

  // May 29, 2008 RSK - add CMB prior if requested
  if (usecmb) {

    HzFUN_INFO_DEF HzFUN;
    set_HzFUN_for_wfit(ONE, OM, OE, w0, wa, &HzFUN);
    // xxx mark rz = Hainv_integral ( ONE, OM, OE, w0,wa,  acmb, ONE ) / LIGHT_km;
    rz = Hainv_integral ( acmb, ONE, &HzFUN ) / LIGHT_km;
    Rcmb = sqrt(OM) * rz ;
    nsig = (Rcmb - Rcmb_best)/sigma_Rcmb ;
    *chi2tot += pow( nsig, TWO );
  }


}  // end of get_chi2wOM


// ==============================================
// Look here AM, 3rd Arg
double get_minwOM( double *w0_atchimin, double *wa_atchimin, double *OM_atchimin ) {

  // Jan 23, 2009: 
  // search min chi2 on refined grid with smaller step size.
  // Return minimized w and OM
  // Note that input w,OM are close to min based on coarse grid
  //
  // Apr 23, 2009: fix refined step (w & OM) and compute
  // nbw and nbm = number of bins.
  // Also fix bugs: nbm instead of nbw 
  //  (did not affect last version since nbm=nbw)
  // and fix bug in w,om loops
  //
  // Oct 09, 2013: Fixed refined step (w & OM) so always 10x smaller
  // than starting grid step. If using default grid, this change
  // preserves the prior hard-coded refined grid steps (=0.001).
  // Also, moved hardcoded nb_factor (used in nbw, nbm calculations)
  // to variable declaration region. 

  double 
     w0cen_tmp ,omcen_tmp, wacen_tmp
    ,w0_tmp, wa_tmp, om_tmp
    ,snchi_tmp, extchi_tmp
    ,snchi_min, extchi_min
    ,muoff_tmp
    ;

  double w0step_tmp  = w0_stepsize/10.0; // pre-Oct 2013: 0.001 ;
  double wastep_tmp  = wa_stepsize/10.0;
  double omstep_tmp = omm_stepsize/10.0; // pre-Oct 2013: 0.001 ;
  int    nbw0, nbwa, nbm,i, j, kk,imin=0,kmin =0,  jmin=0    ;  // Added AM

  double nb_factor = 4.0; //JLM SEP24

  // ---------- BEGIN ------------

  snchi_min = 1.e20;  extchi_min = 1.e20 ;

  nbw0 = (int)(nb_factor*w0_stepsize/w0step_tmp) ;
  nbwa= (int)(nb_factor*wa_stepsize/wastep_tmp) ;
  //nbwa= (int)(nb_factor*w0_stepsize/w0step_tmp) ; // Addedd AM : simpint error
  nbm = (int)(nb_factor*omm_stepsize/omstep_tmp) ;

  printf("   Minimize with refined w0step=%6.4f (%d bins, w0=%7.4f) \n", 
	 w0step_tmp, 2*nbw0, *w0_atchimin);
  if(usewa){
  printf("   Minimize with refined wastep=%6.4f (%d bins, wa=%7.4f) \n",
         wastep_tmp, 2*nbwa, *wa_atchimin);
  }
  printf("   Minimize with refined OMstep=%6.4f (%d bins, OMcen=%6.4f) \n", 
	 omstep_tmp, 2*nbm, *OM_atchimin );
  printf("---------------------------------------\n");

  w0cen_tmp  = *w0_atchimin ;
  wacen_tmp    = *wa_atchimin ;
  omcen_tmp = *OM_atchimin ;

  for ( i = -nbw0; i <= nbw0; i++ ) {
    w0_tmp = w0cen_tmp + (double)i*w0step_tmp ;
    // break up w into w0,wa AM

    for ( kk = -nbwa; kk <= nbwa; kk++ ) {
      wa_tmp = wacen_tmp + (double)kk*wastep_tmp ;

      for ( j = -nbm; j < nbm; j++ ) {
	om_tmp = omcen_tmp + (double)j*omstep_tmp ;

	get_chi2wOM(w0_tmp, wa_tmp, om_tmp, sqsnrms, 
		  &muoff_tmp, &snchi_tmp, &extchi_tmp );

	if ( extchi_tmp < extchi_min ) 
	  { extchi_min = extchi_tmp ;  imin=i; kmin = kk; jmin=j; }
      }
    } // end j
  } // end i
  //for (int ii=0; ii < NCIDLIST; ii++){
    //printf("XXX z=%6.4f\n",z[ii]);}

  // change input values with final w,OM

  *w0_atchimin  = w0cen_tmp + (double)imin * w0step_tmp ;
  *wa_atchimin = wacen_tmp    + (double)kmin * wastep_tmp ;
  *OM_atchimin = omcen_tmp + (double)jmin * omstep_tmp ;
  // printf("XXX om_chi min=%6.4f\n",*OM_atchimin);
  return extchi_min;

} // end of get_minwOM

// ==================================
void printerror(int status) {
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/
  
	

  if (status)    {
      fits_report_error(stderr, status); /* print error report */
      exit(EXIT_ERRCODE_wfit);
  }
  return;
}


int compare_double_reverse (const void *a, const void *b)
{ 				/* qsort in decreasing order */
  int rv = 0;

  if (  *((double *) a)  <  *((double *) b)   )
    {
      rv = 1;
    }
  else if (  *((double *) a)  >  *((double *) b)   )
    {
      rv = -1;
    }
  return rv;
}




/****************************************************************************
 CODIST: This routine computes the dimensionless comoving distance, which 
   turns out to be the generically useful quantity for various cosmological
   distance measures.  The dimensionless comoving distance is defined here as:
   
   d(z) = Integral{0,z}( dz' / E(z') )

   where
   
   E(z) = sqrt( Omega_m*(1+z)^3 + Omega_k*(1+z)^2 + 
                Omega_e*(1+z)^3*(1+w0+wa) * exp(-3*wa*(z/(1+z))) )
	    
   Hogg, astro-ph/9905116, relates E(z) to various cosmlogical distance 
   measures; I've added a parametrized dark energy equation of state parameter
   following Linder, 2003, PRL, v90:
   
     w(a) = w0 + wa*(1-a)
          = w0 + wa*z/(1+z)

   Two uses of the dimensionless comoving distance of interest to me are:

   * Luminosity distances:
     
     Dl(z) = (1+z)*c/H0*sqrt(Omega_k) * sinn{ sqrt(Omega_k)*d(z) }
           = (1+z)*c/H0 * d(z) for flat cosmologies

     See Hogg, eqns 16 and 21 for more.

   * Baryonic Acoustic Oscillations:

     See Eisenstein et al, 2006, ApJ, section 4.5, where they introduce the 
     function A:

     A = sqrt(Omega_m) * E(z_bao)^(-1/3) * [(1/z_bao) * d(z_bao)]^(2/3)

     where z_bao=0.35 and A is constrained by their measurements to be
     A = 0.469 +/- 0.017.

  For convenience, E(z), 1/E(z) and d(z) are provided here.

  Integration over redshift done using Simpson's Rule via 'simpint', which is 
  derived from Numerical Recipces 'qsimp', but allows for passing of additional
  parameters.

  Note that the integral is divergent for some values of Omega_m,Omega_e,
  e.g. "Big Bounce" cosmologies.  No checking is done here to avoid these 
  regions, so watch your step.

  Gajus Miknaitis, Fermilab, Aug '06
******************************************************************************/


double codist(double z, Cosparam *cptr)
/* Returns dimensionless comoving distance. */
{
  double zero = 0.0 ;
  return simpint(one_over_EofZ, zero, z, cptr);
}

double one_over_EofZ(double z, Cosparam *cptr){
  /* This is actually the function that we pass to the integrator. */
  return 1./EofZ(z, cptr);
}

double EofZ(double z, Cosparam *cptr){
  double E;
  double omm, omk, ome, w0, wa;

  omm = cptr->omm;
  ome = cptr->ome;
  w0  = cptr->w0;
  wa  = cptr->wa;

  omk = 1. - omm - ome;

  E =  sqrt( omm*pow((1+z),3.) + 
	     omk*pow((1+z),2.) + 
	     ome*pow((1+z),3.*(1+w0+wa)) * exp(-3*wa*(z/(1+z))) );



  return E;
}



/***************************************************************************
 * SIMPINT is a minor modification of the Numerical Recipes 'qsimp'        *
 * routine, which integrates a function numerically by Simpson's rule.     *
 * 'simpint' now allows for the passing of additional arguments to the     *
 * function via a pointer to a structure.  Also, the use of a static       *
 * variable in the workhose routine 'trapzd' has been eliminated.          *
 *                                                                         *
 * Integration by trapezoidal rule is also provided via 'trapint'. For     *
 * smooth functions, 'simpint' should usually be more accurate, whereas    *
 * 'trapint' may provide better results when integrating noisy data.       *
 * See Press, Teukolsky, Vetterling and  Flannery, "Numerical Recipes in   *
 * C", Chapter 4, section 2 for more detail.                               *
 *                                                                         *
 * Gajus Miknaitis, Fermilab, Aug '06                                      *
 ***************************************************************************/

double simpint(double (*func)(double, Cosparam *), double x1, 
	       double x2, Cosparam *vptr)
/* based on NR 'qsimp' */
{
  int j;
  double s, st;
  double ost=0.0, os=0.0;

  for(j=1; j<=JMAX; j++){
    st = trapezoid(func,x1,x2,j,ost,vptr);
    s = (4.0*st-ost)/3.0;
    if (j > 5)  		/* Avoid spurious early convergence */
      if (fabs(s-os) < EPS*fabs(os) || (s==0.0 && os==0.0))  return s;
    os = s;
    ost = st;
  }

  /* If we got to here, the integral did not converge in JMAX iterations */
  printf("ERROR: simpint did not converge\n");
  return 0.0;
}

double trapint(double (*func)(double, Cosparam *), double x1, 
	       double x2, Cosparam *vptr)
/* based on NR 'qtrap' */
{
  int j;
  double s,olds=0.0;

  for(j=1; j<=JMAX; j++){
    s = trapezoid(func,x1,x2,j,olds,vptr);
    if (j > 5)  		/* Avoid spurious early convergence */
      if (fabs(s-olds) < EPS*fabs(olds) || (s==0.0 && olds==0.0))  return s;
    olds = s;
  }

  /* If we got to here, the integral did not converge in JMAX iterations */
  printf("ERROR: trapint did not converge\n");
  return 0.0;
}


double trapezoid(double (*func)(double, Cosparam *), double x1, double x2, 
		 int n, double s, Cosparam *vptr)
/* Based on NR 'trapzd'.  Here, s must be fed back in from previous 
   iterations, rather than using a static variable as in the NR version. */
{
  double x,tnm,sum,del;
  int it,j;

  if (n == 1) {
    return 0.5*(x2-x1)*( (*func)(x1,vptr) + (*func)(x2,vptr) );
  } else {
    for(it=1,j=1; j<n-1; j++) it <<= 1;
    tnm = it;
    del = (x2-x1)/tnm;
    x = x1+0.5*del;
    for (sum=0.0,j=1; j<=it; j++,x+=del) sum += (*func)(x,vptr);
    s = 0.5*(s+(x2-x1)*sum/tnm);
    return s;
  }
}


void test_codist() {

  // test variables
  int iz;
  double Ztmp, atmp, amax, Zmin, rz, ra;
  double rcodist;
  Cosparam cpar;

  // test integration routines for accuracy and speed

  cpar.omm  =  0.3;
  cpar.ome  =  0.7;
  cpar.w0   = -1.0;
  cpar.wa   =  0.0;

  HzFUN_INFO_DEF HzFUN;
  set_HzFUN_for_wfit(H0, cpar.omm, cpar.ome, cpar.w0, cpar.wa, &HzFUN);

  printf(" test_codist with O(M,E) = %4.2f %4.2f,  w=%4.2f \n, wa=%4.2f \n",
	 cpar.omm, cpar.ome, cpar.w0, cpar.wa );

  for ( iz = 100;  iz <= 1100; iz+=1 ) {

    Ztmp = (double)iz;
    Ztmp = 1089. ;
    Zmin = 0.0;
    atmp = 1./(1. + Ztmp);  amax = 1.0;

    /*
    rz = Hzinv_integral ( H0, cpar.omm, cpar.ome, cpar.w0
			  ,Zmin, Ztmp );
    */
    rz = 1.;

    // xxx ra = Hainv_integral(H0, cpar.omm, cpar.ome, cpar.w0, atmp, amax );
    ra = Hainv_integral ( atmp, amax, &HzFUN );

    rz *= H0/LIGHT_km ;
    ra *= H0/LIGHT_km ;

    rcodist = 2. * ra;
    rcodist = codist(Ztmp, &cpar);

    if ( iz == 100 || iz == 500 || iz == 1100 ) {
      printf("    Z=%7.2f  ra/rz = %f   ra/codist = %f \n", 
	     Ztmp, ra/rz, ra/rcodist );
    }


  }

  printf("  Exit gracefully after test_codist \n");
  exit(0); // xxxxx


} //

// =================================
void getname(char *basename, char *tempname, int nrun) {


   //get clean string name
   strcpy(tempname, basename);

   if ( nrun == 0 ) {
     // works for residfile
     // works for cosparfile
     strcat(tempname, ".init");
     
   } 


} // end of getname

// =================================
int cidindex(char *cid) {

  int i;
  for ( i=0; i < NCIDLIST; i++ ) {
    if ( strcmp(CIDLIST[i],cid) == 0 ) return i;
  }

  return -1;

} // end of cidindex


// ********************************
void write_output_cospar(FILE *fp, RESULTS_DEF *RESULTS, 
			 int usemarg, int format_cospar ) {

  // Created Aug 15 2020
  // format_cospar = 1 : legacy csv format
  // format_cospar = 2 : YAML format

  char sep[] = " " ;
  // ----------- BEGIN -------------

 

  if(usewa){

    if ( format_cospar == 1 ) {
    // legacy format
      if ( !usemarg ) {
	fprintf(fp,"# w0%s w0sig_marg%s wa%s wasig_marg%s OM%s OM_sig%s chi2%s Ndof%s "
		"sigint%s wran%s OMran%s label \n",
		sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep );
	fprintf(fp,"%8.4f%s %7.4f%s %8.4f%s %7.4f%s %7.4f%s %7.4f%s %8.1f%s "
	      "%5d%s %6.3f%s %.2f%s %.2f%s %s\n"
	      , RESULTS->w0_out, sep
	      , RESULTS->w0sig,  sep
	      , RESULTS->wa_out, sep
	      , RESULTS->wasig,  sep
	      , RESULTS->omm_out, sep
	      , RESULTS->omm_sig, sep 
	      , RESULTS->chi2_final, sep
	      , RESULTS->Ndof, sep
	      , RESULTS->sigmu_int, sep
	      , RESULTS->w0rand,   sep
	      , RESULTS->ommrand, sep
	      , RESULTS->label_cospar );
      } else {
	fprintf(fp,"# w0%s w0sig_up%s wa%s wasig_up%s w0sig_low%s OM%s OM_sig%s chi2%s "
	      "Ndof%s sigint%s wran%s OMran%s label\n",
	      sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep, sep );
      
	fprintf(fp,"%8.4f%s %7.4f%s %7.4f%s %8.4f%s %7.4f%s %7.4f%s %7.4f%s %7.4f%s %8.1f%s "
	      "%5d%s %6.3f%s  %.2f%s %.2f%s  %s\n"
	      , RESULTS->w0_out, sep
	      , RESULTS->w0sig_upper, sep
	      , RESULTS->w0sig_lower, sep
              , RESULTS->wa_out, sep
              , RESULTS->wasig_upper, sep
              , RESULTS->wasig_lower, sep
	      , RESULTS->omm_out, sep
	      , RESULTS->omm_sig, sep
	      , RESULTS->chi2_final, sep
	      , RESULTS->Ndof, sep
	      , RESULTS->sigmu_int, sep
	      , RESULTS->w0rand, sep
	      , RESULTS->ommrand, sep
	      , RESULTS->label_cospar);
      }
    }
    else {
    // YAML format
      fprintf(fp, "w0:        %.4f \n", RESULTS->w0_out );
      fprintf(fp, "w0_sig:    %.4f \n", RESULTS->w0sig  );
      fprintf(fp, "wa:        %.4f \n", RESULTS->wa_out );
      fprintf(fp, "wa_sig:    %.4f \n", RESULTS->wasig  );
      fprintf(fp, "omm:      %.4f \n", RESULTS->omm_out );
      fprintf(fp, "omm_sig:  %.4f \n", RESULTS->omm_sig );    
      fprintf(fp, "chi2:     %.1f \n", RESULTS->chi2_final ); 
      fprintf(fp, "sigint:   %.4f \n", RESULTS->sigmu_int );    
      fprintf(fp, "w0rand:    %.4f \n", RESULTS->w0rand );    
      fprintf(fp, "ommrand:  %.4f \n", RESULTS->ommrand );  
      fprintf(fp, "label:    %s \n", RESULTS->label_cospar );    
    }

    
  }
  else {
    if ( format_cospar == 1 ) {
    // legacy format
      if ( !usemarg ) {
	fprintf(fp,"# w%s w0sig_marg%s OM%s OM_sig%s chi2%s Ndof%s "
	      "sigint%s wran%s OMran%s label \n",
	      sep, sep, sep, sep, sep, sep, sep, sep, sep );
	fprintf(fp,"%8.4f%s %7.4f%s %7.4f%s %7.4f%s %8.1f%s "
	      "%5d%s %6.3f%s %.2f%s %.2f%s %s\n"
	      , RESULTS->w0_out, sep
	      , RESULTS->w0sig,  sep
	      , RESULTS->omm_out, sep
	      , RESULTS->omm_sig, sep 
	      , RESULTS->chi2_final, sep
	      , RESULTS->Ndof, sep
	      , RESULTS->sigmu_int, sep
	      , RESULTS->w0rand,   sep
	      , RESULTS->ommrand, sep
	      , RESULTS->label_cospar );
      } else {
	fprintf(fp,"# w%s w0sig_up%s w0sig_low%s OM%s OM_sig%s chi2%s "
	      "Ndof%s sigint%s wran%s OMran%s label\n",
	      sep, sep, sep, sep, sep, sep, sep, sep, sep, sep );
      
	fprintf(fp,"%8.4f%s %7.4f%s %7.4f%s %7.4f%s %7.4f%s %8.1f%s "
	      "%5d%s %6.3f%s  %.2f%s %.2f%s  %s\n"
	      , RESULTS->w0_out, sep
	      , RESULTS->w0sig_upper, sep
	      , RESULTS->w0sig_lower, sep
	      , RESULTS->omm_out, sep
	      , RESULTS->omm_sig, sep
	      , RESULTS->chi2_final, sep
	      , RESULTS->Ndof, sep
	      , RESULTS->sigmu_int, sep
	      , RESULTS->w0rand, sep
	      , RESULTS->ommrand, sep
	      , RESULTS->label_cospar);
      }
    }
    else {
    // YAML format
    fprintf(fp, "w0:        %.4f \n", RESULTS->w0_out );
    fprintf(fp, "w0_sig:    %.4f \n", RESULTS->w0sig  );
    fprintf(fp, "omm:      %.4f \n", RESULTS->omm_out );
    fprintf(fp, "omm_sig:  %.4f \n", RESULTS->omm_sig );    
    fprintf(fp, "chi2:     %.1f \n", RESULTS->chi2_final ); 
    fprintf(fp, "sigint:   %.4f \n", RESULTS->sigmu_int );    
    fprintf(fp, "w0rand:    %.4f \n", RESULTS->w0rand );    
    fprintf(fp, "ommrand:  %.4f \n", RESULTS->ommrand );  
    fprintf(fp, "label:    %s \n", RESULTS->label_cospar );    
    }


  }


  
  return ;

} // end write_output_cospar

// ********************************
void write_output_resid(void) {

  // ----------- BEGIN -------------

  return ;
} // end write_output_resid

// ********************************
void write_output_contour(void) {

  // ----------- BEGIN -------------
  return ;

} // end write_output_contour
