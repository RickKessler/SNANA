/***************************************************************************** 
   wfit [refactored release, Oct 4 2021] 

   Read Hubble diagram in SNANA's "FITRES" format with VARNAMES key,
   and optionally read cov_syst matrix from create_covariance;
   then fit for constraints on Omega_m & w  or  Omega_m and [w0,wa].
   Fit assumes flatness. Each parameter constraint is marginalized
   (default) or chi2-minimized value; uncertainty is +_1 sigma estimate. 

   Optional priors:
    + Gaussian Omega_M
    + CMB R-shift param (E.g. Eq 69 of Komatsu 2009)
    + BAO A-param as in Eisenstein 2006

   To apply CMB prior using sim cosmology and similar sensitivity 
   to Planck 2016,
      wfit.exe <input_HD>  -cmb_sim -sigma_Rcmb 0.007

   Few warnings: 
     - no MCMC sampling
     - no output contours (but could be added)
     - no Planck constraints (see above for similar CMB constraint)
     - no gaurantees

                --------- Code History ---------
   This code was originally developed within SNLS (Conley et al, 2011, 
   arxiv:1104.1443), then passed to SDSS where it was installed into
   SNANA. The wfit program has been superseded by better MCMC sampling 
   codes elsewhere (e.g., CosmoMC, CosmoSIS), but wfit has remained useful 
   for very fast (< 1 minute) cosmology checks on sims and systematics. 
   Many wfit updates were made during the 2010s, but the code was poorly 
   maintained since it was used for w-shifts and crosschecks, but not for
   final cosmology results. The utility of wfit is illustrated by its 
   use in these published method-papers:
     Kessler  2013 (1209.2482)   w-syst from intrinsic scatter models
     Mosher   2014 (1401.4065)   w-syst from SALT2 training on sim samples
     Kessler  2017 (1610.04677)  test w-bias on sims using BBC 
     Brout    2019 (1811.02377)  constraint validation on 200 sim samples
     Lasker   2019 (1811.02380)  w-syst with Chromatic corr
     Popovic  2019 (1910.05228)  w-syst from host association in SDSS
     Popovic  2021 (2102.01776)  w-syst from SN-host correlations
     Vincenzi 2021 (2111.10382)  w-syst for SNCC contamination in DES
     Sanchez  2021 (2111.06858)  cosmology result for LSST-DC2-SNIa analysis
     Brout    2022 (2202.04077)  Pantheon+ systematics
     Dai      2022 (in prep)     SALT3 training syst
     Mitra    2022 (in prep)     SNIa cosmology with photo-z/PLASTICC
     
  In Sep/Oct 2021, R.Kessler and A.Mitra made a few major updates:
    + completely refactor/overhaul code for easier & proper maintainance.
    + include option to float wa in w0waCDM model
    + read and implement cov_syst matrix produced by create_covarance.py

  Regression test of refactored vs. old wfit: 
     There is no change from reading MUDIF output from BBC. If reading 
     unbinned FITRES file with zHDERR (or zERR) column, the incorrect 
     mu_sig_z = [dmu/dz]*zERR term is no longer added and thus the new
     w results changes by ~E-3.

  If using wfit in a publication, please indicate usage with something 
  like "fast cosmology grid-search program in SNANA," and consider 
  optional footnote with "wfit" code name.



     HISTORY (RK=R.Kessler, AM=A.Mitra, ...)
   --------------------------------------------
   
 Nov 24 2021 RK :
   + BOA prior from Alam 2020 is the new default for -bao option.
   + disable computation of sigint unless refit option is used.
     (for large NSN, calc was slow because of matrix inversion each step)

 Dec 01 2021 RK :
   in get_chi2wOM(), skip off-diag computation if chi2(diag) is 
   > 30 sigma above naive chi2=Ndof.
   With 1800^2 cov matrix (unbinned data), speed increase is 
   about x2.6 for wCDM fit.

 Dec 03 2021 RK : new input option -lcdm to fit for OM with w=-1

 Jan 26 2022 RK 
  + fix speed_trick to set nsig = nsig(orig) * 30 instead of 30.
  + write "COSPAR: <results>" to stdout so that regression test
     can scoop up w, werr om, omerr on one line.

 Feb 18 2022 RK - print status updates inside chi2 min loop

 Feb 21 2022 RK 
    + fix subtle bug with speed_trick; apply analytic H0 marg to 
       handle Pantheon distances that are missing M0.
    + inside chi2-min loop, add screen updates with timing info.
    + write NWARNINGS to yaml output

 Mar 09 2022 RK
    + fix chi2_bao_prior to use BAO_PRIOR structure rather than DEFAULTs;
      fixes bug in which -bao_sim was same as -bao.

 Apr 26 2022 RK
    + speed_flag_chi2  +=1 -> offDiag trick
                       +=2 -> interpolation trick
      set default speed_flag_chi2=3 [use both speed tricks]

 May 19 2022 RK - check COV*COVINV if "-debug_flag 1000"

*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "fitsio.h"
#include "longnam.h"

#include "sntools.h"
#include "sntools_cosmology.h"
#include "sntools_output.h"

// ======== global params ==========

#define MXSN 100000 // max number of SN to read & fit
#define MEMC_FILENAME 1000*sizeof(char)

// bit-mask options for speed_flag_chio2
#define SPEED_MASK_OFFDIAG 1  // skip off-diag calc if chi2(diag)>threshold
#define SPEED_MASK_INTERP  2  // interplate r(z) and mu_cos(z)
#define SPEED_FLAG_CHI2_DEFAULT  SPEED_MASK_OFFDIAG + SPEED_MASK_INTERP

// ======== global structures ==========

struct INPUTS {

  // misc flags
  int fitsflag ;
  int blind;  // blind cosmology results by adding sin(big number) 
  int debug_flag ;

  int   speed_flag_chi2; // default = 1; set to 0 to disable
  bool  USE_SPEED_OFFDIAG; // internal: skip off-diag calc if chi2(diag)>threshold
  bool  USE_SPEED_INTERP;  // internal: intero r(z) and mu(z)

  int fitnumber;   // default=1; legacy for iterative fit after sigint calc

  char infile[1000];          // input hubble diagram in fitres format
  char outFile_cospar[1000] ; // output name of cospar file
  char outFile_resid[1000] ;
  char outFile_chi2grid[1000];
  char mucov_file[1000] ;  // input cov matrix; e.g., from create_covariance
  char label_cospar[40]  ;   // string label for cospar file.
  int  ndump_mucov ; // dump this many column/rows

  // grid ranges and steps
  int    h_steps, omm_steps, w0_steps, wa_steps ;
  double h_min,   omm_min,   w0_min,   wa_min   ;    
  double h_max,   omm_max,   w0_max,   wa_max   ;

  int dofit_wcdm ;    // default is true
  int dofit_lcdm ;    // option with -lcdm
  int dofit_w0wa ;    // option with -wa

  int use_bao, use_cmb ; // flags for bao and cmb priors
  //  int use_wa;            // flag to fit for wa (added Sep 2021)
  int use_marg;          // flag to marginalize (instead of only minimize)
  int use_mucov;         // flag to read cov matrix
  double zmin, zmax;               // redshift selection

  int csv_out ; // optional csv format for output cospar and resids

  // Gauss OM prior
  double omm_prior, omm_prior_sig; 

  double snrms, sqsnrms;
  double OMEGA_MATTER_SIM ;

  // computed from user inputs
  double h_stepsize, omm_stepsize, w0_stepsize, wa_stepsize;
  int    format_cospar;

} INPUTS ;


// define workspace
struct  {

  // define variables interpolate rz for large samples
  int     n_exec_interp; // number of interpolate calls for r(z) and mu_cos(z)
  int     n_logz_interp; // number of logz bins for interpolation
  double *logz_list_interp, *z_list_interp, logz_bin_interp; 
  double *rz_list_interp, *mucos_list_interp;
  double  rz_dif_max ;

  // - - - -
  double *w0_prob, *wa_prob, *omm_prob;
  double *w0_sort, *wa_sort;

  double   *snprob,     *extprob,      *snchi,      *extchi ;
  double ***snprob3d, ***extprob3d,  ***snchi3d,  ***extchi3d;

  int    Ndof, Ndof_prior, Ndof_data;
  double sig_chi2min_naive; // sqrt(2*Ndof)
  double nsig_chi2min_skip; // skip off-diag of chi2diag > this value (Jan 26 2022)

  double snchi_min, extchi_min ;
  double snprobtot, snprobmax, extprobtot, extprobmax;
  double w0_atchimin, wa_atchimin,  omm_atchimin, chi2atmin;

  double w0_probsum, wa_probsum,  omm_probsum ;
  double w0_mean,    wa_mean,     omm_mean ; 
  double w0_probmax, wa_probmax,  omm_probmax;

  double w0_sig_marg,  w0_sig_upper,  w0_sig_lower;
  double wa_sig_marg,  wa_sig_upper,  wa_sig_lower;
  double omm_sig_marg, omm_sig_upper, omm_sig_lower;
  double cov_w0wa, cov_w0omm;
  double rho_w0wa, rho_w0omm; 
  
  double *MUCOV; // 1D array for cov matrix
  int    NCOV_NONZERO ;

  double w0_ran,   wa_ran,   omm_ran;
  double w0_final, wa_final, omm_final, chi2_final ;
  double sigmu_int, muoff_final, FoM_final ;

  int NWARN;

} WORKSPACE ;


// define structure to hold hubble diagram (Oct 1 2021)
struct {
  int    NSN, NSN_ORIG ; // NSN after cuts, before cuts
  char   **cid;
  bool   *pass_cut;
  double *mu, *mu_sig, *mu_ref, *mu_sqsig, *z, *logz, *z_sig ;
  int    *nfit_perbin;  
  double zmin, zmax;
} HD;


typedef struct  {
  // Cosmological parameter structure used by codist.c 
  double omm;  // Omega_matter 
  double ome;  // Omega_energy
  double w0;   // equation of state of Dark Energy: 
  double wa;   //    w(z) = w0 + wa(1-a)  
} Cosparam ;



struct {
  double R, sigR;    // CMB R shift parameter and error
  double z, a ;      // redshift and 1/(1+z)
  char comment[200];
} CMB_PRIOR;



// define SDSS4-eBOSS prior/structure
#define NZBIN_BAO_SDSS4 7

double rd_DEFAULT = 150.0; // from ???
// define default values from Table 3 of Alam 2020
double  z_sdss4_DEFAULT[NZBIN_BAO_SDSS4] =
  { 0.38,   0.51,  0.70,  0.85,  1.48,  2.33, 2.33 } ;
double DMrd_sdss4_DEFAULT[NZBIN_BAO_SDSS4] = 
  { 10.27, 13.38, 17.65, 19.50, 30.21, 37.60, 37.30 } ;
double sigDMrd_sdss4_DEFAULT[NZBIN_BAO_SDSS4] = 
  {  0.15,  0.18,  0.30,  1.00,  0.79,  1.90,  1.70 } ;
double DHrd_sdss4_DEFAULT[NZBIN_BAO_SDSS4] = 
  { 24.89, 22.43, 19.78, 19.60, 13.23, 8.93, 9.08 };   
double sigDHrd_sdss4_DEFAULT[NZBIN_BAO_SDSS4] = 
  {  0.58,  0.48,  0.46,  2.10,  0.47, 0.28, 0.34 } ;
double COV_BAO_SDSS4[NZBIN_BAO_SDSS4*2][NZBIN_BAO_SDSS4*2]; 
double COVINV_BAO_SDSS4[NZBIN_BAO_SDSS4*2][NZBIN_BAO_SDSS4*2];

struct {

  bool use_sdss;  // from Eisenstein 2005  astro-ph/0501171 
  bool use_sdss4; // from SDSS-IV-eBOSS    arxiv:2007.08991

  // define params for Eisenstein 2005  astro-ph/0501171 
  double z_sdss, a_sdss, siga_sdss; 

  // define params for Alam 2010 (SDSS-IV-eBOSS)  arxiv:2007.08991 
  double  rd_sdss4; 
  double  z_sdss4[NZBIN_BAO_SDSS4] ;
  double  DMrd_sdss4[NZBIN_BAO_SDSS4] ;
  double  DHrd_sdss4[NZBIN_BAO_SDSS4] ;
  double  sigDMrd_sdss4[NZBIN_BAO_SDSS4] ;
  double  sigDHrd_sdss4[NZBIN_BAO_SDSS4] ;

  // comment for stdout
  char comment[200];
	      
} BAO_PRIOR;

// =========== global variables ================


double H0SIG     = 1000.0 ;       // error used in prior
double c_light   = 299792.458 ;   // speed of light, km/s 
double c_sound   = 155776.730694; 
double TWOTHIRD  = 2./3. ;
double TWO       = 2. ;
double NEGTHIRD  = -1./3. ;
double ONE       = 1.0;
double ZERO      = 0.0;

double SIG_MUOFF ;  // computed from H0 and H0SIG
double SQSIG_MUOFF ;

double NSIG_CHI2MIN_SAFETY = 30.0; // skip off-diag if chi2_diag - chi2min > 30

// define cosmo varnames used for consistent output names
char varname_w[4];          // either w for wCDM or w0 for w0wa model
char varname_wa[4]  ;
char varname_omm[4] ; 

time_t t_start, t_end_init, t_end_fit ;


int MUERR_INCLUDE_zERR;    // True if zERR already included in MUERR
int MUERR_INCLUDE_LENS ;   // True if lensing sigma already included
#define MXCOVSN   400                   // max # of SN with MU-covariance 
#define MXCOVPAIR MXCOVSN*(MXCOVSN-1)/2 // max off-diag covariance pairs 
int NCOVPAIR = 0 ;
int NCOVSN   = 0 ; 

int INDEX_COVSN_MAP[MXCOVSN];      // CIDLIST index vs. NCOVSN index
int INDEX_COVSN_INVMAP[MXSN];      // NCOVSN index vs. CIDLIST

// =========== function prototypes ================

void print_help(void);
void init_stuff(void);
void parse_args(int argc, char **argv) ;
int  compare_double_reverse (const void *, const void *);
void writemodel(char*, float, float, float);
void printerror( int status);
void read_fitres(char *inFile);
void malloc_HDarrays(int opt, int NSN);
void malloc_workspace(int opt);
void parse_VARLIST(FILE *fp);
void read_mucov_sys(char *inFile);
void dump_MUCOV(char *comment);
void invert_mucovar(double sqmurms_add);
double check_invertMatrix(int N, double *COV, double *COVINV );
void set_stepsizes(void);
void set_Ndof(void);
void init_rz_interp(void);
void exec_rz_interp(int k, Cosparam *cospar, double *rz, double *dmu);
void check_refit(void);

void wfit_minimize(void);
void prep_speed_offdiag(double extchi_tmp);
void wfit_normalize(void);
void wfit_marginalize(void);
void wfit_uncertainty(void);
void wfit_final(void);
void wfit_FoM(void);
void wfit_Covariance(void);


void get_chi2wOM(double w0, double wa, double OM, double sqmurms_add,
		   double *mu_off, double *chi2sn, double *chi2tot );
void getname(char *basename, char *tempname, int nrun);

double get_DMU_chi2wOM(double z, double rz, double mu); 
double get_mu_cos(double z, double rz) ;

double get_minwOM( double *w0_atchimin, double *wa_atchimin, 
		   double *OM_atchimin ); 

void   set_priors(void);
void   init_bao_prior(int OPT) ;
void   init_bao_cov(void);
double rd_bao_prior(Cosparam *cpar) ;
double DM_bao_prior(double z, Cosparam *cpar);
double DH_bao_prior(double z, Cosparam *cpar);

void   init_cmb_prior(int OPT) ;
double chi2_bao_prior(Cosparam *cpar);
double chi2_cmb_prior(Cosparam *cpar);

void set_HzFUN_for_wfit(double H0, double OM, double OE, double w0, double wa,
                        HzFUN_INFO_DEF *HzFUN ) ;

int  cidindex(char *cid);

void WRITE_OUTPUT_DRIVER(void);
void write_output_resid(int fitnum);
void write_output_contour(void);
void write_output_cospar(void);
void write_output_fits(void);

void CPU_summary(void);

// cosmology functions
double EofZ(double z, Cosparam *cptr);
double one_over_EofZ(double z, Cosparam *cptr);
double codist(double z, Cosparam *cptr); // maybe change name to Ezinv_integral
void   test_codist(void);

// Nov 1 2021 R.Kessler - new set of functions to integrate over "a"
double EofA(double a, Cosparam *cptr);
double one_over_EofA(double a, Cosparam *cptr);
double Eainv_integral(double amin, double amax, Cosparam *cptr) ;

// from simpint.h

#define EPS  1.0e-6
#define JMAX 20

double simpint(double (*func)(double, Cosparam *), 
	       double x1, double x2, Cosparam *vptr);

double trapint(double (*func)(double, Cosparam *), 
	       double x1, double x2, Cosparam *vptr);

double trapezoid(double (*func)(double, Cosparam *), 
		 double x1, double x2, int n, double s, Cosparam *vptr);


// =============================================
// ================ MAIN =======================
// =============================================

int main(int argc,char *argv[]){

  // ----------------- BEGIN MAIN ----------------

  t_start = time(NULL);

  set_EXIT_ERRCODE(EXIT_ERRCODE_wfit);

  // Give help if no arguments
  if (argc < 2) { print_help();  exit(0);  }

  // init variables
  init_stuff(); 

  //  test_codist();  only to test r(z) integral

  // Parse command line args 
  parse_args(argc,argv);

  if ( INPUTS.debug_flag == -10 ) { malloc_HDarrays(+1,MXSN); }  // legacy

  malloc_workspace(+1);

  printf("# =============================================== \n");
  printf(" SNANA_DIR = %s \n", getenv("SNANA_DIR") );

  while (INPUTS.fitnumber <= 1) {

    /************************/
    /*** Read in the data ***/
    /************************/
    
    read_fitres(INPUTS.infile); 

    // for large samples, setup logz grid to interpolate rz(z)
    init_rz_interp();

    // Set BAO and CMB priors
    set_priors();
  
    // read optional mu-cov matrix (e..g, Cov_syst)
    read_mucov_sys(INPUTS.mucov_file);
   
    // compute grid step size per floated variable
    set_stepsizes();

    // compute number of degrees of freedom
    set_Ndof(); 


    t_end_init = time(NULL);
 
    // minimize chi2 on a grid
    wfit_minimize(); 
    
    // Normalize probability distributions 
    wfit_normalize();

    // Marginalize 
    wfit_marginalize();

    // get uncertainties
    wfit_uncertainty();

    // Compute covriance with fitted parameters
    wfit_Covariance();
    
    // determine "final" quantities, including sigma_mu^int
    wfit_final();

    // estimate FoM
    wfit_FoM();

    t_end_fit = time(NULL);

    // call driver routine for output(s)
    WRITE_OUTPUT_DRIVER();

    // --------------------------------------------------
    // check option to repeat fit with updated snrms = sigmu_int
    check_refit();


    INPUTS.fitnumber++ ;

  } /* end refit while loop */


  // ----------------------------
  // free memory
  malloc_HDarrays(-1,0);
  malloc_workspace(-1);
  CPU_summary();

  printf("DONE. \n"); fflush(stdout);

  return(0);
}   // end main

// ================================
void init_stuff(void) {
  
  char fnam[] = "init_stuff" ;

  // ------------ BEGIN -----------

  INPUTS.blind = INPUTS.fitsflag = INPUTS.debug_flag = 0;

  INPUTS.speed_flag_chi2 = SPEED_FLAG_CHI2_DEFAULT ;

  INPUTS.OMEGA_MATTER_SIM = OMEGA_MATTER_DEFAULT ;

  INPUTS.mucov_file[0]       = 0 ;
  INPUTS.ndump_mucov         = 0 ;
  INPUTS.outFile_cospar[0]   = 0 ;
  INPUTS.outFile_resid[0]    = 0 ;
  INPUTS.outFile_chi2grid[0] = 0 ;
  INPUTS.use_mucov           = 0 ;
  sprintf(INPUTS.label_cospar,"none");
  INPUTS.format_cospar = 1; // csv like format
  INPUTS.fitnumber = 1;

  // Gauss OM params
  INPUTS.omm_prior        = OMEGA_MATTER_DEFAULT ;  // mean
  INPUTS.omm_prior_sig    = 0.04;  // sigma
  
  init_bao_prior(-9); 
  init_cmb_prior(-9);  

  INPUTS.dofit_wcdm = 1 ;
  INPUTS.dofit_lcdm = 0 ;
  INPUTS.dofit_w0wa = 0 ;

  INPUTS.use_bao  = INPUTS.use_cmb ;
  INPUTS.use_marg = 1;


  INPUTS.w0_steps  = 201; INPUTS.w0_min  = -2.0 ;  INPUTS.w0_max  = 0. ;  
  INPUTS.wa_steps  = 101; INPUTS.wa_min  = -4.0 ;  INPUTS.wa_max  = 4.0 ;
  INPUTS.omm_steps =  81; INPUTS.omm_min =  0.0 ;  INPUTS.omm_max = 1.0;
  INPUTS.h_steps   = 121; INPUTS.h_min   = 40.0 ;  INPUTS.h_max   = 100;   

  INPUTS.zmin = 0.0;  INPUTS.zmax = 9.0;

  INPUTS.snrms  = INPUTS.sqsnrms = 0.0 ;


  // - - - - - - - - - - -  - -
  // init misc variables
  SIG_MUOFF   = 5.0 * log10(1. + H0SIG/H0_SALT2);
  SQSIG_MUOFF = SIG_MUOFF * SIG_MUOFF ;

  sprintf(varname_w,   "w"  );  // default wCDM model
  sprintf(varname_wa,  "wa" );
  sprintf(varname_omm, "OM" );

  // - - - - -
  // WORKSPACE
  WORKSPACE.snchi_min  = 1.0e20 ;
  WORKSPACE.extchi_min = 1.0e20 ;

  WORKSPACE.snprobtot  = 0.0 ;
  WORKSPACE.snprobmax  = 0.0 ; 
  WORKSPACE.extprobtot = 0.0 ;
  WORKSPACE.extprobmax = 0.0 ;
  
  WORKSPACE.w0_atchimin =  0.0 ;
  WORKSPACE.wa_atchimin =  0.0 ;
  WORKSPACE.omm_atchimin =  0.0 ;
  WORKSPACE.chi2atmin   =  0.0 ;
  WORKSPACE.NWARN = 0 ;


  return ;

} // end init_stuff


// ==================================
void print_help(void) {

  // Created Oct 1 2021
  // [moved from main]

  static char *help[] = {		/* instructions */
    "",
    "WFIT - Usage: wfit [hubble diagram] (options)",
    "",
    " Default for input file is 3 columns: z, mu, sigma_mu, but see below",
    "",
    "  Options:",
    "   -wa\t\tfit for w0 & wa [w=w0+wa*(1-a)] from ",
    "     \t\t Chevallier & Polarski, 2001 [Int.J.Mod.Phys.D10,213(2001)]",
    "   -lcdm\tfit for OM only (fix w=-1)" ,
    "   -ompri\tcentral value of omega_m prior [default: Planck2018]", 
    "   -dompri\twidth of omega_m prior [default: 0.04]",
    "   -bao\t\tuse BAO prior from Eisenstein et al, 2006",
    "   -bao_sim\tuse BAO prior with simulated cosmology and E06 formalism",
    "   -cmb\t\tuse CMB prior from 5-year WMAP (2008)",
    "   -cmb_sim\tuse CMB prior with simulated cosmology and WMAP formalism",
    "   -om_sim\tOmega_M for cmb_sim (default from sntools.h)",
    "   -minchi2\tget w and OM from minchi2 instead of marginalizing",
    "   -marg\tget w and OM from marginalizing (default)",
    "   -Rcmb\tCMB comstraints: R = Rcmb +/- sigma_Rcmb [= 1.710 +/- 0.019]",
    "   -sigma_Rcmb\tUncertainty on Rcmb",
    "   -snrms\tadd this to reported errors (mags) [default: 0]",
    "   -zmin\tFit only data with z>zmin",
    "   -zmax\tFit only data with z<zmax",    
    "   -blind\tIf set, results are obscured with sin(large random) ",
    "   -fitswrite\tWrite 2D likelihoods to output fits file.",
    "   -mucov_file\tfile with COV_syst e.g., from create_covariance",
    "   -ndump_mucov\t dump this many rows/columns of MUCOV and MUCOVINV",
    "   -mucovar\t\t [Legacy key for previous]",
    "   -refit\tfit once for sigint then refit with snrms=sigint.", 
    "   -speed_flag_chi2   +=1->offdiag trick, +=2->interp trick"
    "\n",
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
    "   -outfile_cospar\tname of output file with fit cosmo params",
    "   -cospar\t\t  [same as previous]",
    "   -cospar_yaml\tname of output YAML file with fit cosmo params",
    "   -outfile_resid\tname of output file with mu-residuals",
    "   -resid\t\t  [same as previous]" ,
    "   -outfile_chi2grid\tname of output file containing chi2-grid",
    "   -label\tstring-label for cospar file.",
    "",
    0
  };

  int i;
  // ----------- BEGIN ------------
  
  for (i = 0; help[i] != 0; i++)
    { printf ("%s\n", help[i]); }

  return;

} // end print_help


// ==================================
void parse_args(int argc, char **argv) {

  // Created Oct 1 2021 [code moved from main]

  int iarg;
  char fnam[] = "parse_args" ;

  // ------------ BEGIN ------------

  // Aug 15 2020: print full command 
  printf(" Full command: \n   ");
  for(iarg=0; iarg < argc; iarg++ ) { printf(" %s", argv[iarg] );  }
  printf("\n\n"); fflush(stdout);

  strcpy(INPUTS.infile,argv[1]); // positional arg is HD

  for (iarg=2; iarg<argc; iarg++) {
    if (argv[iarg][0]=='-') {
      if (strcasecmp(argv[iarg]+1,"ompri")==0) { 
	INPUTS.omm_prior = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"dompri")==0) { 
	INPUTS.omm_prior_sig = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"bao")==0) { 
	INPUTS.use_bao=1;
      } else if (strcasecmp(argv[iarg]+1,"bao_sim")==0) {
        INPUTS.use_bao=2;
      } else if (strcasecmp(argv[iarg]+1,"csv")==0) { 
	INPUTS.csv_out = 1;

      } else if (strcasecmp(argv[iarg]+1, "Rcmb")==0) {
	CMB_PRIOR.R = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1, "sigma_Rcmb")==0) {
	CMB_PRIOR.sigR = atof(argv[++iarg]);

      } else if (strcasecmp(argv[iarg]+1,"cmb")==0) { 
	INPUTS.use_cmb = 1;
	INPUTS.omm_prior_sig = 0.9;  // turn off omm prior (9/2021)
      } else if (strcasecmp(argv[iarg]+1,"wa")==0) {
        INPUTS.dofit_w0wa = 1 ;
	INPUTS.dofit_wcdm = 0 ;
	sprintf(varname_w,"w0");
      }
      else if (strcasecmp(argv[iarg]+1,"lcdm")==0) {
	INPUTS.dofit_lcdm = 1 ;
	INPUTS.dofit_wcdm = 0 ;
      }
      else if (strcasecmp(argv[iarg]+1,"cmb_sim")==0) {
        INPUTS.use_cmb = 2 ;
	INPUTS.omm_prior     = INPUTS.OMEGA_MATTER_SIM ;
	INPUTS.omm_prior_sig = 0.9;  // turn off Gauss omm prior
      } else if (strcasecmp(argv[iarg]+1,"om_sim")==0) {
	INPUTS.OMEGA_MATTER_SIM = atof(argv[++iarg]); 
	INPUTS.omm_prior        = INPUTS.OMEGA_MATTER_SIM ;

      } else if (strcasecmp(argv[iarg]+1,"minchi2")==0) { 
	INPUTS.use_marg = 0 ;
      } else if (strcasecmp(argv[iarg]+1,"marg")==0) { 
	INPUTS.use_marg=1;

      } else if (strcasecmp(argv[iarg]+1,"snrms")==0) { 
	INPUTS.snrms = atof(argv[++iarg]);
	INPUTS.sqsnrms = INPUTS.snrms * INPUTS.snrms ;

      } else if (strcasecmp(argv[iarg]+1,"refit")==0){
        INPUTS.fitnumber -= 1;  

      } else if (strcasecmp(argv[iarg]+1,"zmin")==0) { 
	INPUTS.zmin = atof(argv[++iarg]); 
      } else if (strcasecmp(argv[iarg]+1,"zmax")==0) { 
	INPUTS.zmax = atof(argv[++iarg]); 	

      } else if (strcasecmp(argv[iarg]+1,"blind")==0) { 
	INPUTS.blind=1;
	
      } else if (strcasecmp(argv[iarg]+1,"fitswrite")==0) {   
	/* Output likelihoods as fits files & text files */
	INPUTS.fitsflag=1; 

      } else if (strcasecmp(argv[iarg]+1,"mucov_file")==0) {   
  	strcpy(INPUTS.mucov_file,argv[++iarg]);
	INPUTS.use_mucov =1 ;
      } else if (strcasecmp(argv[iarg]+1,"mucovar")==0) {    // legacy key
  	strcpy(INPUTS.mucov_file,argv[++iarg]);
	INPUTS.use_mucov =1 ;
       
      } else if (strcasecmp(argv[iarg]+1,"ndump_mucov")==0 ) { 
	INPUTS.ndump_mucov = atoi(argv[++iarg]);      

      /* Change H0 grid parameters? */
      } else if (strcasecmp(argv[iarg]+1,"hmin")==0) {   
	INPUTS.h_min = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"hmax")==0) {   
	INPUTS.h_max = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"hsteps")==0) {   
	INPUTS.h_steps = (int)atof(argv[++iarg]);
      }

      /* Change Omega_m grid parameters? */
      else if (strcasecmp(argv[iarg]+1,"ommin")==0) {   
	INPUTS.omm_min = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"ommax")==0) {   
	INPUTS.omm_max = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"omsteps")==0) {   
	INPUTS.omm_steps = (int)atof(argv[++iarg]);
      }

      /* accept w or w0 */
      else if (strcasecmp(argv[iarg]+1,"w0min")==0) {   
	INPUTS.w0_min = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"w0max")==0) {   
	INPUTS.w0_max = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"w0steps")==0) {   
	INPUTS.w0_steps = (int)atof(argv[++iarg]);
      }

      else if (strcasecmp(argv[iarg]+1,"wmin")==0) {   
	INPUTS.w0_min = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"wmax")==0) {   
	INPUTS.w0_max = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"wsteps")==0) {   
	INPUTS.w0_steps = (int)atof(argv[++iarg]);
      }

      /* Change wa grid parameters? */
      else if (strcasecmp(argv[iarg]+1,"wamin")==0) {
        INPUTS.wa_min = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"wamax")==0) {
        INPUTS.wa_max = atof(argv[++iarg]);
      } else if (strcasecmp(argv[iarg]+1,"wasteps")==0) {
        INPUTS.wa_steps = (int)atof(argv[++iarg]);
      }
           
      /* Output filenames */
      else if (strcasecmp(argv[iarg]+1,"outfile_cospar")==0) 
 	{ strcpy(INPUTS.outFile_cospar, argv[++iarg]); }      
      else if (strcasecmp(argv[iarg]+1,"cospar")==0) 
 	{ strcpy(INPUTS.outFile_cospar, argv[++iarg]); }      

      else if (strcasecmp(argv[iarg]+1,"cospar_yaml")==0) 
 	{ strcpy(INPUTS.outFile_cospar, argv[++iarg]); INPUTS.format_cospar=2; }  

      else if (strcasecmp(argv[iarg]+1,"label")==0) 
 	{ sprintf(INPUTS.label_cospar,"%s", argv[++iarg]); }      

      else if (strcasecmp(argv[iarg]+1,"outfile_resid")==0)
	{ strcpy(INPUTS.outFile_resid,argv[++iarg]); }
      else if (strcasecmp(argv[iarg]+1,"resid")==0) 
	{ strcpy(INPUTS.outFile_resid,argv[++iarg]);  }

      else if (strcasecmp(argv[iarg]+1,"outfile_chi2grid")==0)  
	{ strcpy(INPUTS.outFile_chi2grid,argv[++iarg]); }
      else if (strcasecmp(argv[iarg]+1,"chi2grid")==0)  
	{ strcpy(INPUTS.outFile_chi2grid,argv[++iarg]); }

      else if (strcasecmp(argv[iarg]+1,"debug_flag")==0)  
	{ INPUTS.debug_flag = atoi(argv[++iarg]);  } 

      else if (strcasecmp(argv[iarg]+1,"speed_flag_chi2")==0)
	{ INPUTS.speed_flag_chi2 = atoi(argv[++iarg]); }      

      else {
	printf("Bad arg: %s\n", argv[iarg]);
	exit(EXIT_ERRCODE_wfit);
      }
    }
  }


  // - - - - - - 
  char var_w[20], str_marg_min[20];
	   
  INPUTS.USE_SPEED_OFFDIAG = (INPUTS.speed_flag_chi2 & SPEED_MASK_OFFDIAG)>0;

  printf(" ****************************************\n");
  if ( INPUTS.dofit_w0wa )  { 
    printf("   Fit w0waCDM  model:  w(z) = w0 + wa(1-a) \n"); 
    sprintf(var_w,"%s,%s", varname_w, varname_wa);
  }
  else if ( INPUTS.dofit_lcdm ) {
    printf("  Fit LCDM model for OM  with w=-1 \n" );
    sprintf(var_w,"%s", varname_w);
    INPUTS.w0_steps = 1;   INPUTS.w0_min = INPUTS.w0_max = -1.0 ;
    INPUTS.wa_steps = 1;   INPUTS.wa_min = INPUTS.wa_max =  0.0 ;
  }
  else {
    printf("   Fit wCDM model:  w(z) = w0  \n"); 
    sprintf(var_w,"%s", varname_w);
    INPUTS.wa_steps = 1;   INPUTS.wa_min = INPUTS.wa_max = 0.0 ;
  }

  if ( INPUTS.use_marg ) 
    { sprintf(str_marg_min,"MARGINALIZE"); }
  else
    { sprintf(str_marg_min,"MINIMIZE"); }

  printf("   %s for final %s,%s \n",
	 str_marg_min, var_w, varname_omm);

  if ( INPUTS.snrms > 0.0  ) {
    printf("Add %4.2f mag-error (snrms) in quadrature to MU-error \n", 
	   INPUTS.snrms);
  }
  printf("---------------------------------------\n\n");
  fflush(stdout);

  return;

} // end parse_args


// ==================================
void  malloc_workspace(int opt) {

  // Created Oct 1 2021 [code moved from main]
  // opt > 0 -> malloc
  // opt < 0 -> free

  int memd = sizeof(double);
  int i, kk, j;
  float f_mem;
  char fnam[] = "malloc_workspace" ;
  // ------------ BEGIN ------------

  if ( opt > 0 ) {
    WORKSPACE.w0_prob   = (double *)calloc(INPUTS.w0_steps,  memd);
    WORKSPACE.wa_prob   = (double *)calloc(INPUTS.wa_steps,  memd);
    WORKSPACE.omm_prob  = (double *)calloc(INPUTS.omm_steps, memd);

    WORKSPACE.w0_sort   = (double *)calloc(INPUTS.w0_steps, memd);
    WORKSPACE.wa_sort   = (double *)calloc(INPUTS.wa_steps, memd);

    /* Probability arrays */
    
    f_mem = 0; 
    f_mem += malloc_double3D(+1, INPUTS.w0_steps, INPUTS.wa_steps, 
			      INPUTS.omm_steps, &WORKSPACE.snchi3d );
    f_mem += malloc_double3D(+1, INPUTS.w0_steps, INPUTS.wa_steps, 
			      INPUTS.omm_steps, &WORKSPACE.extchi3d );
    f_mem += malloc_double3D(+1, INPUTS.w0_steps, INPUTS.wa_steps, 
			      INPUTS.omm_steps, &WORKSPACE.extprob3d );
    f_mem += malloc_double3D(+1, INPUTS.w0_steps, INPUTS.wa_steps, 
			      INPUTS.omm_steps, &WORKSPACE.snprob3d );

    for(i=0; i < INPUTS.w0_steps; i++ ) {
      for(kk=0; kk < INPUTS.wa_steps; kk ++ ) {
	for (j=0; j < INPUTS.omm_steps; j++ ) {
	  WORKSPACE.snchi3d[i][kk][j]   = 0.0 ;
	  WORKSPACE.extchi3d[i][kk][j]  = 0.0 ;
	  WORKSPACE.extprob3d[i][kk][j] = 0.0 ;
	  WORKSPACE.snprob3d[i][kk][j]  = 0.0 ;  
	}
      }
    }

  }
  else {
    free(WORKSPACE.w0_prob);
    free(WORKSPACE.wa_prob);
    free(WORKSPACE.omm_prob);
    free(WORKSPACE.w0_sort);
    free(WORKSPACE.wa_sort);
    
    free(WORKSPACE.extprob);
    free(WORKSPACE.snprob);
    free(WORKSPACE.extchi);
    free(WORKSPACE.snchi);

  }

  return;

} // end malloc_workspace

// ==================================
void  malloc_HDarrays(int opt, int NSN) {

  // Created Oct 01 2021 by R.Kessler
  // malloc arrays to store Hubble diagram data
  // opt > 0 -> malloc
  // opt < 0 -> free

  int i;
  // --------- BEGIN --------

  if ( opt > 0 ) {
    HD.cid = (char**) malloc( NSN * sizeof(char*) );
    for(i=0; i < NSN; i++ ) { HD.cid[i] = (char*)malloc( 20*sizeof(char) ); }

    HD.pass_cut    = (bool *)calloc(NSN,sizeof(bool));
    HD.mu          = (double *)calloc(NSN,sizeof(double));
    HD.mu_sig      = (double *)calloc(NSN,sizeof(double));
    HD.mu_ref      = (double *)calloc(NSN,sizeof(double));
    HD.mu_sqsig    = (double *)calloc(NSN,sizeof(double));   
    HD.z           = (double *)calloc(NSN,sizeof(double));
    HD.logz        = (double *)calloc(NSN,sizeof(double));
    HD.z_sig       = (double *)calloc(NSN,sizeof(double));
    HD.nfit_perbin = (int*) calloc(NSN,sizeof(int));
  }
  else {
    free(HD.pass_cut);
    free(HD.mu); 
    free(HD.mu_sig); 
    free(HD.mu_ref); 
    free(HD.mu_sqsig); 
    free(HD.z); 
    free(HD.logz); 
    free(HD.z_sig); 
    free(HD.nfit_perbin); 
    for(i=0; i < NSN; i++ ) { free(HD.cid[i]); } 
    free(HD.cid);
  }

} // end malloc_HDarrays


// ==================================
void read_fitres(char *inFile) {

  // Created Oct 1 2021 by R.Kessler
  // Refactored routine to read Hubble diagram from fitres-formatted file
  // using SNANA read utilities.

  int IVAR_ROW=-8, IVAR_MU=-8, IVAR_MUERR=-8, IVAR_MUREF;
  int IVAR_zHD=-8, IVAR_zHDERR=-8, IVAR_NFIT=-8 ;
  int IFILETYPE, NVAR_ORIG, LEN, NROW, irow ;  
  int VBOSE = 1;
  char TBNAME[] = "HD" ;  // table name is Hubble diagram
  char fnam[] = "read_fitres" ;

  // --------------- BEGIN --------------


  TABLEFILE_INIT();
  NROW = SNTABLE_NEVT(inFile,TBNAME);
  IFILETYPE = TABLEFILE_OPEN(inFile,"read");
  NVAR_ORIG = SNTABLE_READPREP(IFILETYPE,TBNAME);

  malloc_HDarrays(+1,NROW); 

  IVAR_ROW   = SNTABLE_READPREP_VARDEF("CID:C*20 ROW:C*20",
				       HD.cid, NROW, VBOSE);

  IVAR_MU    = SNTABLE_READPREP_VARDEF("MU:D DLMAG:D MUDIF:D",    
				       HD.mu,     NROW, VBOSE) ;
  IVAR_MUERR = SNTABLE_READPREP_VARDEF("MUERR:D DLMAGERR:D MUDIFERR:D", 
				       HD.mu_sig, NROW, VBOSE) ;
  
  IVAR_MUREF = SNTABLE_READPREP_VARDEF("MUREF:D", 
				       HD.mu_ref, NROW, VBOSE);
  IVAR_NFIT  = SNTABLE_READPREP_VARDEF("NFIT:I", 
				       HD.nfit_perbin, NROW, VBOSE );
  // - - - -
  IVAR_zHD    = SNTABLE_READPREP_VARDEF("zHD:D zCMB:D z:D Z:D",    
					HD.z,     NROW, VBOSE) ;
  IVAR_zHDERR = SNTABLE_READPREP_VARDEF("zHDERR:D zCMBERR:D zERR:D ZERR:D", 
					HD.z_sig, NROW, VBOSE) ;

  // check for required elements
  if ( IVAR_MU < 0 ) {
    sprintf(c1err,"Could not find required distance column:");
    sprintf(c2err,"MU or DLMAG or MUDIF", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }
  if ( IVAR_MUERR < 0 ) {
    sprintf(c1err,"Could not find required distance-uncertainty column:");
    sprintf(c2err,"MUERR or DLMAGERR or MUDIFERR", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }
  if ( IVAR_zHD < 0 ) {
    sprintf(c1err,"Could not find required redshift column:");
    sprintf(c2err,"zHD or z or Z or zCMB", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  // - - - - - - - - - 
  // read table ; note that HD.NSN_ORIG = NROW
  HD.NSN_ORIG = SNTABLE_READ_EXEC();

  // for MUDIF output from BBC, MU is actually MUDIF,
  // so set MU += MUREF
  if ( IVAR_MUREF > 0 ) {
    printf("\t Input is MUDIF File: MU -> MU_REF + MUDIF \n");
    for(irow=0; irow < NROW; irow++ ) { HD.mu[irow] += HD.mu_ref[irow]; }
  }

  // - - - - - -

  int    PASSCUTS,  NFIT, NROW2=0 ;  
  double ztmp;
  HD.zmin = 99999.0;  HD.zmax = -9999.90 ;

  for(irow=0; irow < NROW; irow++ ) {

    NFIT = 999;
    if( IVAR_NFIT > 0 ) { NFIT = HD.nfit_perbin[irow] ; }

    ztmp = HD.z[irow] ;
    PASSCUTS = 
      (ztmp >= INPUTS.zmin) &&
      (ztmp <= INPUTS.zmax) &&
      (NFIT > 1) ;
    
    HD.pass_cut[irow] = false;
    if ( PASSCUTS ) {
      HD.pass_cut[irow]  = true;
      sprintf(HD.cid[NROW2], "%s", HD.cid[irow] );
      HD.mu[NROW2]       = HD.mu[irow];
      HD.mu_sig[NROW2]   = HD.mu_sig[irow];
      HD.z[NROW2]        = ztmp;
      HD.logz[NROW2]     = log10(ztmp); // Apr 22 2022
      HD.z_sig[NROW2]    = HD.z_sig[irow];
      HD.mu_sqsig[NROW2] = HD.mu_sig[irow]*HD.mu_sig[irow];
      if ( ztmp < HD.zmin ) { HD.zmin = ztmp; }
      if ( ztmp > HD.zmax ) { HD.zmax = ztmp; }

      NROW2++ ;

    } // end PASSCUTS
  } // end irow

  HD.NSN  = NROW2;
  printf("   Keep %d of %d z-bins with z cut\n", 
	 HD.NSN, NROW);  

  fflush(stdout);
  //  debugexit(fnam);

  return;

} // end read_fitres



// ==================================
void read_mucov_sys(char *inFile){

  // Created October 2021 by A.Mitra and R.Kessler
  // Read option Cov_syst matrix, and invert it.

#define MXSPLIT_mucov 20

  int MSKOPT_PARSE = MSKOPT_PARSE_WORDS_STRING+MSKOPT_PARSE_WORDS_IGNORECOMMA;
  char ctmp[200], SN[2][12], locFile[1000] ;
  int NSPLIT, NROW_read=0, NDIM_ORIG = 0, NDIM_STORE = HD.NSN, NMAT_read=0,  NMAT_store = 0;
  float f_MEM;
  double cov;
  int N, N0, N1, i0, i1, j, iwd, NWD, i, k0, k1, kk,gzipFlag ;
  char  **ptrSplit;
  FILE *fp;
  char fnam[] = "read_mucov_sys" ;

  // ---------- BEGIN ----------------

  // init legacy INVMAP ... should be able to remove this
  // after refactor.
  for ( i=0; i<HD.NSN; i++ )  { INDEX_COVSN_INVMAP[i] = -9 ; }  

  WORKSPACE.NCOV_NONZERO = 0;

  // bail if there is no cov matrix to read
  if ( !INPUTS.use_mucov  ) { return; }

  printf("\n# ======================================= \n");
  printf("  Process MUCOV systematic file  \n"); 

  // Open File using the utility
  int OPENMASK = OPENMASK_VERBOSE + OPENMASK_IGNORE_DOCANA ;
  fp = snana_openTextFile(OPENMASK, "", inFile,
			  locFile, &gzipFlag );

  if ( !fp ) {
    sprintf(c1err,"Cannot open mucov_file") ;
    sprintf(c2err,"%s", inFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // allocate strings to read each line ... in case there are comments
  ptrSplit = (char **)malloc(MXSPLIT_mucov*sizeof(char*));
  for(j=0;j<MXSPLIT_mucov;j++){
    ptrSplit[j]=(char *)malloc(200*sizeof(char));
  }
  
  i0 = i1 = 0 ;
  k0 = k1 = 0 ;
  while ( fgets(ctmp, 100, fp) != NULL ) {
    // ignore comment lines 
    if ( commentchar(ctmp) ) { continue; }
    
    // break line into words
    splitString(ctmp, " ",MXSPLIT_mucov,  // (I)
		&NSPLIT, ptrSplit);       // (O)

    if ( NROW_read == 0 ) {
      sscanf(ptrSplit[0],"%d",&NDIM_ORIG);
      printf("\t Found COV dimension %d\n", NDIM_ORIG);
      printf("\t Store COV dimension %d\n", NDIM_STORE);
      
      if ( NDIM_ORIG != HD.NSN_ORIG ) {
	sprintf(c1err,"NDIM(COV)=%d does not match NDIM(HD)=%d ??",
		NDIM_ORIG, HD.NSN_ORIG);
	sprintf(c2err,"Above NDIM are before cuts.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      int MEMD = NDIM_STORE * NDIM_STORE * sizeof(double);
      WORKSPACE.MUCOV = (double*) malloc(MEMD);
    }
    else {
      NMAT_read++ ;
      sscanf( ptrSplit[0],"%le",&cov);      
      if(HD.pass_cut[i0] && HD.pass_cut[i1] ) {
	kk = k1*NDIM_STORE + k0;
	WORKSPACE.MUCOV[kk] = cov;
	if ( cov != 0.0 ) { WORKSPACE.NCOV_NONZERO++ ; }
	NMAT_store += 1;
	k0++;
	if ( k0 == NDIM_STORE ) { k1++; k0=0; }
      }
      i0++;
      if ( i0 == NDIM_ORIG ) { i1++; i0=0; }
    }
    
    NROW_read += 1 ;

  } // end of read loop


  printf("\t Read %d non-zero COV elements.\n", WORKSPACE.NCOV_NONZERO );
  fflush(stdout);

  // if all COV elements are zero, this is a request for stat-only fit,
  // so disable cov
  if ( WORKSPACE.NCOV_NONZERO == 0 ) { INPUTS.use_mucov = 0; return; }

  // - - - - - - - - - - - - - - - - -
  // sanify check
  if ( NMAT_read != NDIM_ORIG*NDIM_ORIG )  {
    sprintf(c1err,"Read %d cov elements, but expected %d**2=%d",
	    NMAT_read, NDIM_ORIG, NDIM_ORIG*NDIM_ORIG);
    sprintf(c2err,"Check %s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if ( NMAT_store != NDIM_STORE*NDIM_STORE )  {
    sprintf(c1err,"Store %d cov elements, but expected %d**2=%d",
            NMAT_store, NDIM_STORE, NDIM_STORE*NDIM_STORE);
    sprintf(c2err,"Check %s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // Add diagonal errors from the Hubble diagram
  double COV_STAT ;
  for ( i=0; i<HD.NSN; i++ )  {
    kk = i*NDIM_STORE + i;
    COV_STAT = HD.mu_sqsig[i] ;
    WORKSPACE.MUCOV[kk] += COV_STAT ;
  }

  int LDMP_MUCOV = (INPUTS.ndump_mucov>0);
  if ( LDMP_MUCOV ) { dump_MUCOV("MUCOV"); }
  
  // - - - - -
  // invert MUCOV matrix ... this will ovewrite WORKSPACE.MUCOV
  // with its inverse.
  invert_mucovar(INPUTS.sqsnrms);

  if ( LDMP_MUCOV ) { dump_MUCOV("MUCOV^-1"); }

  return ; 
}
// end of read_mucov_sys


void dump_MUCOV(char *comment ) {

  char fnam[]="dump_MUCOV";
  int i0, i1, NROW, kk ;
  int MAX_ROW = INPUTS.ndump_mucov ;
  if(HD.NSN < MAX_ROW){NROW = HD.NSN;}
  else{NROW= MAX_ROW;}
  
  // dump
  printf("\n DUMP %s \n", comment);
  
  for (i0=0; i0 < NROW; i0++)  {
      for(i1=0; i1 < NROW; i1++) {
	kk = i0*NROW + i1;
	printf("%9.5f ", WORKSPACE.MUCOV[kk] );
      }
      printf("\n");
  }
  printf("\n"); 
}// end of dump_MUCOV




//===================================
void set_priors(void) {

  bool noprior = true;

  double OM    = INPUTS.OMEGA_MATTER_SIM ;
  double OE    = 1 - OM ;
  double w0    = w0_DEFAULT ;
  double wa    = wa_DEFAULT ;

  Cosparam cparloc;
  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa ; 

  char fnam[]="set_priors";

  // =========== BEGIN ============
  
  if ( INPUTS.use_bao )   {
    init_bao_prior(1); }

  if ( INPUTS.use_cmb )   { init_cmb_prior(1); }

  // - - - - - - - - - - - -
  // Report priors to stdout

  printf("\n# ============================================ \n");
  printf(" PRIOR(s): \n");
  // - - - -
  if ( H0SIG < 100. ) {
    printf("\t Fit with H0 prior: sigH0/H0=%4.1f/%4.1f  "
	   "=> sig(MUOFF)=%5.3f \n",  H0SIG, H0_SALT2, SIG_MUOFF );
    noprior = false;
  }
 
  
  // - - - - - - - -
  if ( INPUTS.use_bao ) {
    noprior = false;
    char *comment = BAO_PRIOR.comment ; 
    printf("   %s\n", comment) ;
  } 
  else if ( INPUTS.omm_prior_sig < 1. ) {
    noprior = false;
    printf("   Gaussian Omega_m prior: %5.3f +/- %5.3f\n",
	   INPUTS.omm_prior, INPUTS.omm_prior_sig);
  }

  // - - - -
  if ( INPUTS.use_cmb ) {
    noprior = false;
    char *comment = CMB_PRIOR.comment ; 
    printf("   %s\n", comment) ; 
  }

  if ( noprior ) { printf("\t None.\n"); }

  fflush(stdout);

  return;

} //end of set_priors()


// =========================================
void init_cmb_prior(int OPT) {

  // Created Oct 2021 by R.Kessler
  // OPT= -9 --> set all params to -9 
  //                (before reading user input)
  // OPT= +1 --> set params to measured or sim values
  //                (after reading user input)

  double rz, a, z;
  double OM = INPUTS.OMEGA_MATTER_SIM ;
  double OE = 1 - OM ;
  double w0 = w0_DEFAULT ;
  double wa = wa_DEFAULT ;

  Cosparam cparloc;
  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa ; 

  char *comment = CMB_PRIOR.comment;
  char fnam[] = "init_cmb_prior" ;

  // ---------- BEGIN ------------

  if ( OPT == -9 ) {
    CMB_PRIOR.z      = 1089.0 ;
    CMB_PRIOR.a      = 1.0/(1.0 + CMB_PRIOR.z) ;
    CMB_PRIOR.R      = -9.0 ;
    CMB_PRIOR.sigR   = .019 ;
    comment[0]       =  0;
    return ;
  }

  
  // - - - - - - -
  if ( INPUTS.use_cmb == 1 ) {
    // Set WMAP default from Komatsu 2009, Eq 69 and Table 10,
    // unless user has already specified value on command line input
    if(CMB_PRIOR.R    < 0.0) { CMB_PRIOR.R     = 1.710; }
    if(CMB_PRIOR.sigR < 0.0) { CMB_PRIOR.sigR  = 0.019; }
    sprintf(comment, "CMB WMAP-prior:  R=%5.3f +- %5.3f " ,
	   CMB_PRIOR.R, CMB_PRIOR.sigR);
  }
  else if ( INPUTS.use_cmb == 2 ) {
    //recompute R from sim
    HzFUN_INFO_DEF HzFUN ;
    a = CMB_PRIOR.a ;  // 1/(1+z)
    set_HzFUN_for_wfit(ONE, OM, OE, w0, wa, &HzFUN) ;
    rz = Hainv_integral ( a, ONE, &HzFUN ) / LIGHT_km;
    CMB_PRIOR.R = sqrt(OM) * rz ;
    sprintf(comment, "CMB sim-prior:  R=%5.3f +- %5.3f " ,
	    CMB_PRIOR.R, CMB_PRIOR.sigR);
  }

  return;

} // end init_cmb_prior

// =========================================
void init_bao_prior(int OPT) {

  // Created Oct 2021 by R.Kessler
  // OPT= -9 --> set all params to -9
  //               (before reading user input)
  // OPT= +1 --> set params to measured or sim values.
  //               (after reading user input)
  //
  // Nov 24 2021: BAO from SDSSIV is new default.
  //              legacy prior used if debug_flag = -7
  //
  // Mar 09 2022: remove LEGACY flag

  int iz;
  double rz, tmp1, tmp2, z;
  double OM = INPUTS.OMEGA_MATTER_SIM ;
  double OE = 1 - OM ;
  double w0 = w0_DEFAULT ;
  double wa = wa_DEFAULT ;

  Cosparam cparloc;
  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa ; 

  char *comment = BAO_PRIOR.comment;
  char fnam[] = "init_bao_prior" ;

  // -----------------BEGIN ------------------

  if ( OPT == -9 ) {
    BAO_PRIOR.use_sdss  = false;
    BAO_PRIOR.use_sdss4 = false;

    BAO_PRIOR.a_sdss    = -9.0 ;
    BAO_PRIOR.siga_sdss = -9.0 ;
    BAO_PRIOR.z_sdss    = -9.0 ;

    BAO_PRIOR.rd_sdss4    = -9.0 ;
    for(iz=0; iz < NZBIN_BAO_SDSS4; iz++ ) {
      BAO_PRIOR.z_sdss4[iz]  = -9.0 ;
      BAO_PRIOR.DMrd_sdss4[iz] = -9.0 ;
      BAO_PRIOR.DHrd_sdss4[iz] = -9.0 ;
    }

    BAO_PRIOR.comment[0] = 0;
    
    return;
  }
  
  // - - - - - -
  BAO_PRIOR.use_sdss4   = true ;

  // Default values from SDSS4
  BAO_PRIOR.rd_sdss4 = -9.0 ;
  for(iz=0; iz < NZBIN_BAO_SDSS4; iz++ ) {
    BAO_PRIOR.z_sdss4[iz]       = z_sdss4_DEFAULT[iz];
    BAO_PRIOR.DMrd_sdss4[iz]    = DMrd_sdss4_DEFAULT[iz];
    BAO_PRIOR.DHrd_sdss4[iz]    = DHrd_sdss4_DEFAULT[iz];
    BAO_PRIOR.sigDMrd_sdss4[iz] = sigDMrd_sdss4_DEFAULT[iz];
    BAO_PRIOR.sigDHrd_sdss4[iz] = sigDHrd_sdss4_DEFAULT[iz];
  }
  // keep loading all 7 z ranges ...
  sprintf(comment,"BAO prior from SDSS-IV (Alam 2020)" );
  
  // - - - - - -
  // check option to compute BAO params from sim cosmology
  if ( INPUTS.use_bao == 2 ) { 
    double rd,z, DM, DH;
    HzFUN_INFO_DEF HzFUN;
    rd = rd_bao_prior(&cparloc);
    for (iz=0; iz < NZBIN_BAO_SDSS4; iz++ ) {
      z = BAO_PRIOR.z_sdss4[iz];
      DM = DM_bao_prior(z, &cparloc);
      DH = DH_bao_prior(z, &cparloc);
      BAO_PRIOR.DMrd_sdss4[iz] = DM/rd ;
      BAO_PRIOR.DHrd_sdss4[iz] = DH/rd ;
    }
    sprintf(comment,"BAO prior from SDSS-IV using sim cosmology" );
  }
  init_bao_cov();
  return ;

} // end init_bao_prior

void init_bao_cov(void){
  char fnam[] = "init_bao_cov";


  return;

}
double rd_bao_prior( Cosparam *cpar) {
  // rd ~ 150 Mpc                  Pg. 9, Aubourg et al. [1411.1074]
  // rd = int_0^inf [c_s /H(z)];                       Eq. 13, Alam 2020.
  // c_s= 1/ (sqrt( 3 * (1 + 3*Om_b / 4*Om_gamma) ) )  Davis et al, Page 4.
  // Om_b ~ 0.02/h^2                                   Davis T. Note
  // Om_g 2.469 * 10e-5 * T_CMB / 2.725                Davis et al, after eq. 15 
  double H0      = H0_Planck;
  double c_sound = 0.9 * c_light / sqrt(3.0); // Davis internal note
  double h       = H0 / 100.0 ;
  double Om_b = 0.02 / (h * h), Om_gamma = 2.46*0.00001;
  //c_sound = 1/ (sqrt( 3 * (1 + 3*Om_b / 4*Om_gamma) ) ); 
  double z_d    = 1060.0 ; 
  double amin   = 1.0E-6,  amax = 1.0/(1+z_d);
  double Hinv_integ, Einv_integ, rd = 1.0;
  HzFUN_INFO_DEF HzFUN_INFO;
  bool DO_INTEGRAL = true ; 
  int  LDMP = 0;
  char fnam[] = "rd_bao_prior" ;
  // ---------- BEGIN ----------  
  rd = rd_DEFAULT;
  if ( DO_INTEGRAL  ) {
    Einv_integ = Eainv_integral(amin, amax, cpar);
    Hinv_integ = Einv_integ / H0 ;
    rd = c_sound * Hinv_integ  ; 
  }
  return rd;
}

double DM_bao_prior(double z, Cosparam *cpar){
  // Ayan Mitra Nov, 2021
  double DM = 1.0,  DA, H0 = H0_Planck;
  double amin   = 1.,  amax = 1.0/(1+z);
  double rd = codist(z, cpar);
  double Hinv_integ, Einv_integ;
  HzFUN_INFO_DEF HzFUN_INFO;
  bool DO_INTEGRAL = true ;
  int  LDMP = 0;
  char fnam[] = "DM_bao_prior" ;
  
  if ( DO_INTEGRAL  ) {
    Einv_integ = Eainv_integral(amax, amin, cpar);
    Hinv_integ = Einv_integ / H0 ;
    DM = c_light * Hinv_integ  ;
  }  
  return DM;
}


double DH_bao_prior(double z, Cosparam *cpar){
  double H0 = H0_Planck, H, DH;
  HzFUN_INFO_DEF HzFUN;
  set_HzFUN_for_wfit(H0, cpar->omm, cpar->ome, cpar->w0, cpar->wa, &HzFUN);
  H  = Hzfun(z, &HzFUN);
  DH = (c_light / H);
  return DH;
}



// =========================================
void set_stepsizes(void) {

  // Created Oct 1 2021
  // Compute grid size for each cosmoPar dimension

  char fnam[] = "set_stepsizes";

  // ------------ BEGIN ------------

  INPUTS.w0_stepsize = INPUTS.wa_stepsize = 0.0;
  INPUTS.omm_stepsize = INPUTS.h_stepsize = 0. ;

  if (INPUTS.w0_steps  > 1) 
    { INPUTS.w0_stepsize  = (INPUTS.w0_max-INPUTS.w0_min)/(INPUTS.w0_steps-1);}
  if (INPUTS.wa_steps  > 1) 
    { INPUTS.wa_stepsize  = (INPUTS.wa_max-INPUTS.wa_min)/(INPUTS.wa_steps-1);}
  if (INPUTS.omm_steps > 1) 
    { INPUTS.omm_stepsize = (INPUTS.omm_max-INPUTS.omm_min)/(INPUTS.omm_steps-1);}
  if (INPUTS.h_steps   > 1) 
    { INPUTS.h_stepsize   = (INPUTS.h_max-INPUTS.h_min)/(INPUTS.h_steps-1);}
  
   
  printf("   --------   Grid parameters   -------- \n");
  printf("  %s_min: %6.2f   %s_max: %6.2f  %5i steps of size %8.5f\n",
	 varname_w, varname_w, 
	 INPUTS.w0_min, INPUTS.w0_max, INPUTS.w0_steps, INPUTS.w0_stepsize);
  
  if( INPUTS.dofit_w0wa ){
    printf("  %s_min: %6.2f   %s_max: %6.2f  %5i steps of size %8.5f\n",
	   varname_wa, varname_wa, 
	   INPUTS.wa_min, INPUTS.wa_max, INPUTS.wa_steps, INPUTS.wa_stepsize);
  }
  printf("  %s_min: %6.2f   %s_max: %6.2f  %5i steps of size %8.5f\n",
	 varname_omm, varname_omm, 
	 INPUTS.omm_min, INPUTS.omm_max, INPUTS.omm_steps, INPUTS.omm_stepsize);
  printf("  h_min:  %6.2f   h_max:  %6.2f  %5i steps of size %8.5f\n",
	 INPUTS.h_min,  INPUTS.h_max,  INPUTS.h_steps,  INPUTS.h_stepsize);
  
  printf("   ------------------------------------\n");

  fflush(stdout);

  return;
} // end set_stepsizes


// ==================
void set_Ndof(void) {
  // Compute number of degrees of freedom for wfit
  int Ndof;
  int Ndof_data = HD.NSN - 3 ; // h, w, om = 3 fitted parameters
  int Ndof_prior = 0;
  if ( INPUTS.dofit_w0wa  ) { Ndof_data-- ; } 
  if ( INPUTS.use_bao     ) { Ndof_prior += 2*NZBIN_BAO_SDSS4 ; }
  if ( INPUTS.use_cmb     ) { Ndof_prior += 1 ; }

  Ndof = Ndof_data + Ndof_prior;
  
  WORKSPACE.Ndof       = Ndof ;
  WORKSPACE.Ndof_data  = Ndof_data;
  WORKSPACE.Ndof_prior = Ndof_prior;
  WORKSPACE.sig_chi2min_naive = sqrt(2.0*(double)Ndof);

  printf("\n# ====================================== \n");
  printf("   Ndof(data,prior,final) = %d, %d, %d \n",
	 Ndof_data, Ndof_prior, Ndof);
  fflush(stdout);

  return ;

} // end set_Ndof


// ==================================
void init_rz_interp(void) {

  // Created Apr 22 2022
  // For large SN sample, initialize lists to interpolate rz(z) in chi2 fun.
  double zmin = HD.zmin;
  double zmax = HD.zmax;
  int    NSN  = HD.NSN;
  int MEMD, n_logz=0, nz_check=0;
  double logz_bin, logz_min, logz_max, logz, z ;
  char fnam[] = "init_rz_interp";

  // ----------- BEGIN ------------

  WORKSPACE.n_exec_interp = 0;

  if ( NSN > 1000 ) 
    { INPUTS.USE_SPEED_INTERP  = (INPUTS.speed_flag_chi2 & SPEED_MASK_INTERP )>0; }
  else
    { INPUTS.USE_SPEED_INTERP = false; }

  if ( !INPUTS.USE_SPEED_INTERP ) { return; }

  logz_min = log10(zmin);
  logz_max = log10(zmax);

  n_logz   = (int)(200.0 * (zmax - zmin));

  if ( n_logz > NSN/2 ) { INPUTS.USE_SPEED_INTERP=false; return; }

  logz_bin = ( logz_max - logz_min ) / (double)(n_logz-1) ;
  MEMD     = n_logz * sizeof(double);

  WORKSPACE.n_logz_interp     = n_logz;
  WORKSPACE.logz_bin_interp   = logz_bin ;
  WORKSPACE.logz_list_interp  = (double*)malloc(MEMD);
  WORKSPACE.z_list_interp     = (double*)malloc(MEMD);
  WORKSPACE.rz_list_interp    = (double*)malloc(MEMD);
  WORKSPACE.mucos_list_interp = (double*)malloc(MEMD);
  
  printf("\n# ========================================================= \n");
  printf(" load %d logz bins (%.5f <= z <= %.5f) to interpolate rz(z)\n", 
	 n_logz, zmin, zmax);
  fflush(stdout);
  
  for(logz=logz_min; logz <= logz_max+1.0E-6; logz += logz_bin ) {
    WORKSPACE.logz_list_interp[nz_check] = logz;
    WORKSPACE.z_list_interp[nz_check]    = pow(10.0,logz);
    nz_check++ ;
  }

  if ( nz_check != n_logz ) {
    sprintf(c1err,"Expected to load %d logz bins", n_logz);
    sprintf(c2err,"but loaded %d logz bins", nz_check);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return;

} // end init_rz_interp

void exec_rz_interp(int k, Cosparam *cparloc, double *rz, double *mucos) {

  // Created Apr 22 2022
  // return interpolated rz and dmu for SN index k
  // cparloc is used only as a diagnostic for first few chi2 loops.

  int    n_logz   = WORKSPACE.n_logz_interp;
  double logz_min = WORKSPACE.logz_list_interp[0];
  double logz_bin = WORKSPACE.logz_bin_interp;

  double logz, rz0, rz1, frac, mucos0, mucos1, rz_dif, rz_exact ;
  double rz_loc, mucos_loc;
  int  iz;
  char fnam[] = "exec_rz_interp";

  // ----------------- BEGIN --------------

  logz = HD.logz[k];
  iz   = (int)((logz - logz_min)/logz_bin );
  if ( iz < 0 || iz >= n_logz ) {
    sprintf(c1err,"Invalid iz=%d (range is 0 to %d)", iz, n_logz);
    sprintf(c2err,"z=%f  logz=%f  for CID=%s", 
	    HD.z[k], HD.logz[k], HD.cid[k] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  if ( iz < n_logz-1 ) {
    frac = (logz - WORKSPACE.logz_list_interp[iz])/logz_bin;

    rz0    = WORKSPACE.rz_list_interp[iz];
    rz1    = WORKSPACE.rz_list_interp[iz+1];
    rz_loc = rz0 + frac*(rz1-rz0); 

    mucos0    = WORKSPACE.mucos_list_interp[iz];
    mucos1    = WORKSPACE.mucos_list_interp[iz+1];
    mucos_loc = mucos0 + frac*(mucos1-mucos0); 
  }
  else {
    rz_loc    = WORKSPACE.rz_list_interp[iz]; // last z-bin
    mucos_loc = WORKSPACE.mucos_list_interp[iz]; // last z-bin
  }
  
  // load output function args
  *rz    = rz_loc;
  *mucos = mucos_loc;

  if ( k == HD.NSN-1 ) { WORKSPACE.n_exec_interp++ ; }

  // print diagnostic for first few events.
  int LDMP = (WORKSPACE.n_exec_interp < 5) ;
  if ( LDMP ) {
    if ( k==0 ) { WORKSPACE.rz_dif_max = 0.0 ; }    
    rz_exact   = codist(HD.z[k], cparloc);
    rz_dif     = fabs(rz_loc / rz_exact - 1.0);
    if ( rz_dif > WORKSPACE.rz_dif_max ) 
      { WORKSPACE.rz_dif_max = rz_dif; }  

    if ( k == HD.NSN-1 ) {
      printf("\t    Diagnostic: max|rz(interp)/rz(exact) - 1| "
	     "= %8.2le \n",  fnam, WORKSPACE.rz_dif_max); fflush(stdout);
    }
  }

  return;
} // end exec_rz_interp

// ==================================
void check_refit(void) {
  
  // Created Oct 6 2021
  // [move old JLM code from main to here]

  double sigmu_int = WORKSPACE.sigmu_int ;

  // ----------- BEGIN --------------

  if ( INPUTS.fitnumber == 0 ) {
    if ( sigmu_int > INPUTS.snrms ) {
      printf("REFITTING DATA WITH UPDATED INTRINSIC SCATTER\n\n");
      printf("\t input snrms = %.4f --> %.4f\n", INPUTS.snrms, sigmu_int);
      INPUTS.snrms   = sigmu_int;
      INPUTS.sqsnrms = sigmu_int * sigmu_int;
    } else {
      printf("SKIPPING REFIT - NO ADDITIONAL INTRINSIC SCATTER\n\n");
      printf("\t initial snrms   = %.4f \n", INPUTS.snrms);
      printf("\t final sigma_int = %0.4f \n\n", sigmu_int);
      INPUTS.fitnumber++ ;
    }
  }

  return ;

} // end check_refit

// ==================================
void wfit_minimize(void) {

  // Created Oct 2 2021
  // Driver function to minimize chi2 on grid,
  // Outputs loaded to WORKSPACE struct.

  int Ndof                 = WORKSPACE.Ndof;
  double sig_chi2min_naive = WORKSPACE.sig_chi2min_naive ;
  bool   USE_SPEED_OFFDIAG = INPUTS.USE_SPEED_OFFDIAG ;
  
  int    use_mucov         = INPUTS.use_mucov;
  Cosparam cpar;
  Cosparam cpar_fixed;
  double snchi_tmp, extchi_tmp, muoff_tmp;
  bool UPDATE_STDOUT;
  int  i, kk, j;
  int  imin = -9, kmin = -9, jmin = -9;
  char fnam[] = "wfit_minimize" ;

  // ---------- BEGIN --------------

  int NBTOT = INPUTS.w0_steps * INPUTS.wa_steps * INPUTS.omm_steps;
  int NB=0;

  printf("\n# ======================================= \n");
  printf(" Get prob at %d grid points, and approx mimimized values: \n", 
	 NBTOT );
  printf("\t USE_SPEED_OFFDIAG = %d \n", INPUTS.USE_SPEED_OFFDIAG);
  printf("\t USE_SPEED_INTERP  = %d \n", INPUTS.USE_SPEED_INTERP);
  fflush(stdout);
    
  // prep speed trick
  if ( use_mucov && USE_SPEED_OFFDIAG ) {
    printf("   Prepare speed trick to reduce off-diagonal calculation\n");

    INPUTS.USE_SPEED_OFFDIAG = false; // disable speed flag for approx min chi2 
    cpar_fixed.w0 = -1.0; cpar_fixed.wa=0.0; cpar_fixed.omm=0.3; 
    get_chi2wOM(cpar_fixed.w0,cpar_fixed.wa, cpar_fixed.omm, INPUTS.sqsnrms, 
		&muoff_tmp, &snchi_tmp, &extchi_tmp ); 
    printf("    Very approx chi2min(SNonly,SN+prior: w0=-1,wa=0,omm=0.3) "
	   "= %.1f %.1f \n", snchi_tmp, extchi_tmp);
    fflush(stdout);

    INPUTS.USE_SPEED_OFFDIAG = true ; // restore speed flag
    prep_speed_offdiag(extchi_tmp);
  }

  // - - - - - - - - 
  time_t t0 = time(NULL);  // monitor time to build prob grid

  for( i=0; i < INPUTS.w0_steps; i++){
    cpar.w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    for( kk=0; kk < INPUTS.wa_steps; kk++){    
      cpar.wa = (INPUTS.wa_min + kk*INPUTS.wa_stepsize);
      for(j=0; j < INPUTS.omm_steps; j++){
	cpar.omm = INPUTS.omm_min + j*INPUTS.omm_stepsize; 
	cpar.ome = 1 - cpar.omm;
	
	get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, INPUTS.sqsnrms, 
		      &muoff_tmp, &snchi_tmp, &extchi_tmp ); 

	WORKSPACE.snchi3d[i][kk][j]  = snchi_tmp ; 
	WORKSPACE.extchi3d[i][kk][j] = extchi_tmp ;
	  
	// Keep track of minimum chi2 
	if(snchi_tmp < WORKSPACE.snchi_min) 
	  { WORKSPACE.snchi_min = snchi_tmp ; }
	
	if(extchi_tmp < WORKSPACE.extchi_min)  { 
	  WORKSPACE.extchi_min = extchi_tmp ;  
	  imin=i; jmin=j; kmin=kk; 
	}

	// stdout update with timing information
	NB++;
	if ( NB < 1000 ) 
	  { UPDATE_STDOUT = ( NB % 100 == 0 ); }
	else if ( NB < 10000 ) 
	  { UPDATE_STDOUT = ( NB % 1000 == 0 ); }
	else
	  { UPDATE_STDOUT = ( NB % 10000 == 0 ); }

	if ( UPDATE_STDOUT ) {
	  time_t t_update = time(NULL);
	  double dt = t_update - t0 ;
	  printf("\t finished chi2 bin %8d of %8d  (%.0f sec)\n",
		 NB, NBTOT, dt); fflush(stdout) ;
	}


      } // j loop
    }  // end of k-loop
  }  // end of i-loop


  // get w,OM at min chi2 by using more refined grid
  // Pass approx w,OM,  then return w,OM at true min
  //  printf("  Get minimized w,OM from refined chi2grid \n");
  WORKSPACE.w0_atchimin = INPUTS.w0_min  + imin * INPUTS.w0_stepsize ;
  WORKSPACE.wa_atchimin = INPUTS.wa_min  + kmin * INPUTS.wa_stepsize ;
  WORKSPACE.omm_atchimin= INPUTS.omm_min + jmin * INPUTS.omm_stepsize ;
  
  if ( !INPUTS.use_marg  ) {
    WORKSPACE.chi2atmin  = 
      get_minwOM( &WORKSPACE.w0_atchimin, 
		  &WORKSPACE.wa_atchimin,  
		  &WORKSPACE.omm_atchimin ); 
    printf("   Final chi2min(SN+prior) = %f \n", WORKSPACE.chi2atmin);
    fflush(stdout);
  }


  return;

} // end wfit_minimize

// =============================
void prep_speed_offdiag(double chi2min_approx) {

  // Prepare speed trick for computing chi2 by ignoring 
  // off-diagonal terms when diagonal sum is greater than threshold.
  // Input is appriximate chi2min.

  bool USE_SPEED_OFFDIAG   = INPUTS.USE_SPEED_OFFDIAG ;
  int Ndof                 = WORKSPACE.Ndof;
  double sig_chi2min_naive = WORKSPACE.sig_chi2min_naive ;

  double nsig, Xdof=(double)Ndof;
  char fnam[] = "prep_speed_offdiag";

  // -------------- BEGIN -----------

  if ( !USE_SPEED_OFFDIAG ) { return; }

  // compute nsig_chi2 above naive chi2min=Ndof and compare
  // with nsig_chi2_skip used for speed trick
  nsig = (chi2min_approx - Xdof) / sig_chi2min_naive;

  printf("\t Naive nsig(chi2min) = %.1f  # (chi2min-Ndof)/sqrt(2*Ndof)\n",  
	 nsig);
  
  if ( INPUTS.use_mucov ) {
    double nsig_chi2min_skip;
    double NSIG_MULTIPLIER = 10.0 ;
    if ( nsig < 1.0 )
      { nsig_chi2min_skip = NSIG_MULTIPLIER; } 
    else
      { nsig_chi2min_skip = nsig * NSIG_MULTIPLIER ; } // Jan 26 2026

    WORKSPACE.nsig_chi2min_skip = nsig_chi2min_skip;

    printf("\t Skip off-diag chi2-calc if nsig(diag) > %.1f\n", 
	  nsig_chi2min_skip );

  }
  fflush(stdout);

  return;

} // end prep_speed_offdiag

// ==================================
void wfit_normalize(void) {

  // Created Oct 2 2021
  // Convert chi2map to normalized probability map.

  int i, kk, j;
  double chidif; 
  char fnam[] = "wfit_normalize" ;

  // ----------- BEGIN ------------

  /* First, convert chi2 to likelihoods */
  for(i=0; i< INPUTS.w0_steps; i++){
    for(kk = 0; kk < INPUTS.wa_steps; kk++){
      for(j=0; j < INPUTS.omm_steps; j++){
	
	// Probability distribution from SNe alone 
	chidif = WORKSPACE.snchi3d[i][kk][j] - WORKSPACE.snchi_min ;
	WORKSPACE.snprob3d[i][kk][j] = exp(-0.5*chidif);
	WORKSPACE.snprobtot += WORKSPACE.snprob3d[i][kk][j];
      
	// Probability distribution of SNe + external prior
	chidif = WORKSPACE.extchi3d[i][kk][j] - WORKSPACE.extchi_min ;
	WORKSPACE.extprob3d[i][kk][j] = exp(-0.5*chidif);
	WORKSPACE.extprobtot += WORKSPACE.extprob3d[i][kk][j];
      }
    }
  }

  /* Now normalize so that these sum to 1.  Note that if the
     grid does not contain all of the probability, 
     then these normalizations are incorrect */
  
  for(i=0; i  < INPUTS.w0_steps; i++){
    for(kk = 0; kk < INPUTS.wa_steps; kk++){
      for(j=0; j < INPUTS.omm_steps; j++){
	WORKSPACE.snprob3d[i][kk][j] /= WORKSPACE.snprobtot;
	WORKSPACE.extprob3d[i][kk][j] /= WORKSPACE.extprobtot;
      }
    }
  }
  
  return;
} // end wfit_normalize


// ==================================
void wfit_marginalize(void) {

  // Created Oct 2 2021
  // use prob map to marginalize; store means in WORKSPACE struct.

  double P1max, P2max, Pt1mp, Pt2mp, w0, wa, omm, chi2 ;
  int    i, kk, j;
  char fnam[] = "wfit_marginalize" ;

  // ----------- BEGIN ------------

  printf("\n# =================================================== \n");
  printf(" Get Marginalized values (speed_flag_chi2=%d): \n", 
	 INPUTS.speed_flag_chi2) ;
  fflush(stdout);

  WORKSPACE.w0_probsum = WORKSPACE.wa_probsum = WORKSPACE.omm_probsum = 0.0 ;
  WORKSPACE.w0_mean    = WORKSPACE.wa_mean    = WORKSPACE.omm_mean    = 0.0 ;
  WORKSPACE.w0_probmax = WORKSPACE.wa_probmax = WORKSPACE.omm_probmax = 0.0 ;
  P1max = P2max = 0.0;
    
  // Get Marginalise w0
  printf("\t Get Marginalized %s\n", varname_w ); fflush(stdout);
  for(i=0; i < INPUTS.w0_steps; i++){
    w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    WORKSPACE.w0_prob[i] = 0. ;                 
    for(kk=0; kk < INPUTS.wa_steps; kk++){
      for(j=0; j < INPUTS.omm_steps; j++){
	Pt1mp       = WORKSPACE.extprob3d[i][kk][j] ;
	WORKSPACE.w0_prob[i] += Pt1mp ;
	if ( Pt1mp > P1max ) {
	  P1max = Pt1mp ;
	  WORKSPACE.w0_probmax = w0 ;
	} // if Loop

      } // j
    } // kk

    WORKSPACE.w0_mean    += WORKSPACE.w0_prob[i]*w0;
    WORKSPACE.w0_probsum += WORKSPACE.w0_prob[i] ; 
  }  // i

  // - - - -
  // Mean of w0 distribution: this is our estimate of w0 
  WORKSPACE.w0_mean /= WORKSPACE.w0_probsum;
	
  // - - - - - - - -
  // Get Marginalise wa 	
  if( INPUTS.dofit_w0wa ) {
    printf("\t Get Marginalized %s\n", varname_wa ); fflush(stdout);
    for(kk=0; kk < INPUTS.wa_steps; kk++){
      wa = INPUTS.wa_min + kk*INPUTS.wa_stepsize;
      WORKSPACE.wa_prob[kk] = 0. ; 
      for(i=0; i < INPUTS.w0_steps; i++){
	for(j=0; j < INPUTS.omm_steps; j++){
	  Pt2mp       = WORKSPACE.extprob3d[i][kk][j] ;
	  WORKSPACE.wa_prob[kk] += Pt2mp ;
	  if ( Pt2mp > P2max ) {
	    P2max = Pt2mp ;
	    WORKSPACE.wa_probmax = wa ; 
	  } // if Loop
	} // j
      } // i

      WORKSPACE.wa_mean    += WORKSPACE.wa_prob[kk]*wa;
      WORKSPACE.wa_probsum += WORKSPACE.wa_prob[kk];
    } // kk

    // - - -
    // Mean of wa distribution: this is our estimate of wa
    WORKSPACE.wa_mean /= WORKSPACE.wa_probsum;   
  } // end use_wa if block


  // - - - - - - - - 
  // Marginalize over w,wa to get omega_m 
  printf("\t Get Marginalized %s\n", varname_omm ); fflush(stdout);
  for(j=0; j < INPUTS.omm_steps; j++){
    omm = INPUTS.omm_min + j*INPUTS.omm_stepsize;
    for(kk=0; kk < INPUTS.wa_steps; kk++){
      for(i=0; i < INPUTS.w0_steps; i++){
	WORKSPACE.omm_prob[j] += WORKSPACE.extprob3d[i][kk][j]; 
      }
    }
    WORKSPACE.omm_mean    += WORKSPACE.omm_prob[j]*omm;
    WORKSPACE.omm_probsum += WORKSPACE.omm_prob[j];
  }
  
  // Mean of omega_m distribution: this is our estimate of omega_m
  WORKSPACE.omm_mean /= WORKSPACE.omm_probsum;
  
  //---------------------
  // print update 
  if ( INPUTS.blind ) { return ; }

  char str[20];
  if ( INPUTS.use_marg ) {
    sprintf(str,"marg");
    w0  = WORKSPACE.w0_mean;  wa = WORKSPACE.wa_mean;
    omm = WORKSPACE.omm_mean;
  }
  else {
    sprintf(str,"min") ;
    w0  = WORKSPACE.w0_atchimin ; wa = WORKSPACE.wa_atchimin ;
    omm = WORKSPACE.omm_atchimin ;
    chi2 = WORKSPACE.chi2atmin;
  }

  printf("  CHECK:  %s(%s) = %f \n",  varname_w, str, w0);
  if ( INPUTS.dofit_w0wa )
    { printf("  CHECK:  %s(%s) = %f \n",  varname_wa, str, wa);  }

  printf("  CHECK:  %s(%s) = %f \n",  varname_omm, str, omm);

  if(!INPUTS.use_marg) { printf("  CHECK:  ChiSq(min) = %f \n", chi2); }

  printf("\n" );

  fflush(stdout);

  return;

} // end wfit_marginalize


// =========================================
void wfit_uncertainty(void) {

  // Created Oct 2 2021
  // Estimate uncertainties

  double w0, wa, omm, delta, sqdelta ;
  int i, kk, j ;
  char fnam[] = "wfit_uncertainty" ;

  // ------------- BEGIN ------------

  printf("\n# ============================================ \n");
  printf(" Estimate uncertainties: \n");
  fflush(stdout);

  WORKSPACE.w0_sig_marg = 0.0 ;
  WORKSPACE.w0_sig_upper = WORKSPACE.w0_sig_lower = 0.0;

  WORKSPACE.wa_sig_marg = 0.0;
  WORKSPACE.wa_sig_upper = WORKSPACE.wa_sig_lower = 0.0;

  WORKSPACE.omm_sig_marg = 0.0;
  WORKSPACE.omm_sig_upper = WORKSPACE.omm_sig_lower = 0.0;

  // Get std dev for the weighted mean.  This is a reasonable   
  //  measure of uncertainty in w if distribution is ~Gaussian
  sqdelta = 0.0 ;
  for(i=0; i < INPUTS.w0_steps; i++){
    w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    WORKSPACE.w0_prob[i] /= WORKSPACE.w0_probsum;  
    delta       = w0 - WORKSPACE.w0_mean ;
    sqdelta    += WORKSPACE.w0_prob[i] * (delta*delta);
    WORKSPACE.w0_sort[i] = WORKSPACE.w0_prob[i];  
  }
  WORKSPACE.w0_sig_marg = sqrt(sqdelta/WORKSPACE.w0_probsum);
  printf("\t marg %s_sig estimate = %.4f\n", 
	 varname_w, WORKSPACE.w0_sig_marg);
    

  if ( INPUTS.dofit_w0wa ) {
    sqdelta = 0.0 ;
    for(kk=0; kk < INPUTS.wa_steps; kk++){
      wa = INPUTS.wa_min + kk*INPUTS.wa_stepsize;
      WORKSPACE.wa_prob[kk] /=WORKSPACE.wa_probsum; 
      delta    = wa - WORKSPACE.wa_mean ;
      sqdelta += WORKSPACE.wa_prob[kk]*(delta*delta);
      WORKSPACE.wa_sort[kk] = WORKSPACE.wa_prob[kk];        // make a copy 
    }
    WORKSPACE.wa_sig_marg = sqrt(sqdelta/WORKSPACE.wa_probsum);
    printf("\t marg %s_sig estimate = %.4f\n", 
	   varname_wa, WORKSPACE.wa_sig_marg);
  } // end dofit_w0wa 
     
  // omm ...  
  sqdelta = 0.0;
  for(i=0; i < INPUTS.omm_steps; i++){ 	
    omm = INPUTS.omm_min + i*INPUTS.omm_stepsize;    
    WORKSPACE.omm_prob[i] /= WORKSPACE.omm_probsum;  
    delta    = omm - WORKSPACE.omm_mean ;
    sqdelta += WORKSPACE.omm_prob[i] * (delta*delta);
  }
  WORKSPACE.omm_sig_marg = sqrt(sqdelta/WORKSPACE.omm_probsum);
  printf("\t marg %s_sig estimate = %.4f\n", 
	 varname_omm, WORKSPACE.omm_sig_marg);

  // - - - - - - - - -
  double delta_w0=1.0E3, delta_wa=1.0E3 ;
  double w0_sum, w0_sort_1sigma, wa_sum, wa_sort_1sigma;
  int iw0_mean, iwa_mean, memd=sizeof(double) ;
  // Find location of mean in the prob vector


  for (i=0; i < INPUTS.w0_steps; i++){
    w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    if (fabs(w0 - WORKSPACE.w0_mean) < delta_w0){
      delta_w0 = fabs(w0 - WORKSPACE.w0_mean);
      iw0_mean = i;
    }
  }

  if ( INPUTS.dofit_w0wa ) {
    for (kk=0; kk < INPUTS.wa_steps; kk++){
      wa = INPUTS.wa_min + kk*INPUTS.wa_stepsize;
      if (fabs(wa - WORKSPACE.wa_mean) < delta_wa){
	delta_wa = fabs(wa - WORKSPACE.wa_mean);
	iwa_mean = kk;
      }
    }
  }

    // Sort probability
  qsort(WORKSPACE.w0_sort, INPUTS.w0_steps, memd,
	compare_double_reverse);
  if (INPUTS.dofit_w0wa) {
    qsort(WORKSPACE.wa_sort, INPUTS.wa_steps, memd,
	  compare_double_reverse);
  }

  // Add up probability until you get to 68.3%, use that value to 
  // determine the upper/lower limits on w 
  w0_sum=0.;
  w0_sort_1sigma=-1;
  for (i=0; i < INPUTS.w0_steps; i++){
    w0_sum += WORKSPACE.w0_sort[i];
    if (w0_sum > 0.683) {
      w0_sort_1sigma = WORKSPACE.w0_sort[i];  // w0_sort value for 68.3%
      break;
    }
  }

  if (INPUTS.dofit_w0wa){
    wa_sum=0.;
    wa_sort_1sigma=-1;
    for (kk=0; kk < INPUTS.wa_steps; kk++){
      wa_sum += WORKSPACE.wa_sort[kk];
      if (wa_sum > 0.683) {
	wa_sort_1sigma = WORKSPACE.wa_sort[kk]; // w0_sort value for 68.3
	break;
      }
    }
  }

  // - - - - - - - -
  // ERROR checking ... 
 
  if (w0_sort_1sigma < 0 ){
    printf("ERROR: w0 grid doesn't enclose 68.3% of probability!\n");
    exit(EXIT_ERRCODE_wfit);
  }

  if ( INPUTS.dofit_w0wa ){
    if (wa_sort_1sigma < 0 ){
      printf("ERROR: wa grid doesn't enclose 68.3% of probability!\n");
      exit(EXIT_ERRCODE_wfit);
    }
  }


  // Prob array is in order of ascending w.  Count from i=iw0_mean
  // towards i=0 to get lower 1-sigma bound on w 
  for (i=iw0_mean; i>=0; i--){
    if (WORKSPACE.w0_prob[i] < w0_sort_1sigma){
      WORKSPACE.w0_sig_lower = WORKSPACE.w0_mean - 
	(INPUTS.w0_min + i*INPUTS.w0_stepsize);
      break;
    }
  }

  if ( INPUTS.dofit_w0wa ) {
    for (kk=iwa_mean; kk>=0; kk--){
      if (WORKSPACE.wa_prob[kk] < wa_sort_1sigma){
	WORKSPACE.wa_sig_lower = WORKSPACE.wa_mean - 
	  (INPUTS.wa_min + kk*INPUTS.wa_stepsize);
	break;
      }
    }
  }

  if ( INPUTS.dofit_lcdm ) { return; }

  // Error checking 
  if(i==0){
    printf("WARNING: lower 1-sigma limit outside range explored\n");
    WORKSPACE.w0_sig_lower = 100;
    WORKSPACE.NWARN++ ;
  }
  
  if ( INPUTS.dofit_w0wa ) {
    if(kk==0){
      printf("WARNING: lower 1-sigma limit outside range explored\n");
      WORKSPACE.wa_sig_lower = 100;
      WORKSPACE.NWARN++ ;
    }
  }

  if ( WORKSPACE.w0_sig_lower <= INPUTS.w0_stepsize ) {   
    printf("WARNING: 1. w0 grid is too coarse to resolve "
	   "lower 1-sigma limit\n");
    WORKSPACE.w0_sig_lower = INPUTS.w0_stepsize;
    WORKSPACE.NWARN++ ;
  }
  if (INPUTS.dofit_w0wa){
    if (WORKSPACE.wa_sig_lower <= INPUTS.wa_stepsize) {  
      printf("WARNING: 1. wa grid is too coarse to resolve "
	     "lower 1-sigma limit\n");
      WORKSPACE.wa_sig_lower = INPUTS.wa_stepsize;
      WORKSPACE.NWARN++ ;
    }
  }

    
  // Count from i=iw0_mean towards i=w0_steps to get   
  //  upper 1-sigma bound on w 
  for (i=iw0_mean; i < INPUTS.w0_steps; i++){
    if (WORKSPACE.w0_prob[i] < w0_sort_1sigma){
      WORKSPACE.w0_sig_upper = 
	(INPUTS.w0_min + i*INPUTS.w0_stepsize) - WORKSPACE.w0_mean;
      break;
    }
  }

  if ( INPUTS.dofit_w0wa ) {
    for (kk=iwa_mean; kk < INPUTS.wa_steps; kk++){
      if (WORKSPACE.wa_prob[kk] < wa_sort_1sigma){
        WORKSPACE.wa_sig_upper = 
	  (INPUTS.wa_min + kk*INPUTS.wa_stepsize) - WORKSPACE.wa_mean;
        break;
      }
    }
  }

    
  // Error checking 
  if(i==(INPUTS.w0_steps-1)){
    printf("WARNING: upper 1-sigma limit outside range explored\n");
    WORKSPACE.w0_sig_lower = 100;
    WORKSPACE.NWARN++ ;
  } 
  

  if ( WORKSPACE.w0_sig_upper <= INPUTS.w0_stepsize){
    printf("WARNING: 2. w0 grid is too coarse to resolve "
	   "upper 1-sigma limit\n %f, %f\n", 
	   WORKSPACE.w0_sig_upper, INPUTS.w0_stepsize);
    WORKSPACE.w0_sig_upper = INPUTS.w0_stepsize; 
    WORKSPACE.NWARN++ ;
  }

  if ( INPUTS.dofit_w0wa ) {
    if(kk==(INPUTS.wa_steps-1)){
      printf("WARNING: upper 1-sigma limit outside range explored\n");
      WORKSPACE.wa_sig_lower = 100;	
      WORKSPACE.NWARN++ ;
    }
    if (WORKSPACE.wa_sig_upper <= INPUTS.wa_stepsize){
      printf("WARNING: 2. wa grid is too coarse to resolve "
	     "upper 1-sigma limit\n %f, %f\n", 
	     WORKSPACE.wa_sig_upper, INPUTS.wa_stepsize);
      WORKSPACE.wa_sig_upper = INPUTS.wa_stepsize;
      WORKSPACE.NWARN++ ;
    }
  }    

    
  //printf("\n---------------------------------------\n");
  printf("  Prob %s-err estimates: lower = %f, upper = %f\n", 
	 varname_w, WORKSPACE.w0_sig_lower, WORKSPACE.w0_sig_upper);
  if ( INPUTS.dofit_w0wa ) {
    printf("Prob %s-err estimates: lower = %f, upper = %f\n", 
	   varname_wa, WORKSPACE.wa_sig_lower, WORKSPACE.wa_sig_upper);
  }
  
  fflush(stdout);

  return ;

} // end wfit_uncertainty


// ==========================================x 
void wfit_Covariance(void){

  // Created Oct 18 2021 by A.Mitra and R.Kessler
  // Compute the fitted covariance Cov(w0,Om).
  // and if wa option is set, also compute cov(w0,wa)/
  // Store covariances in WORKSPACE.cov_* and store
  // reduced covar in WORKSPACE.rho*
  //
  int i,kk, j ;
  Cosparam cpar;  
  double cov_w0wa=0., cov_w0omm =0., rho_w0wa=0., rho_w0omm=0.;
  double probsum = 0., prob=0;
  double diff_w0, diff_wa, diff_omm, sig_product ;
  bool dofit_w0wa = INPUTS.dofit_w0wa;
  
  // ----------- BEGIN --------------
  
  for( i=0; i < INPUTS.w0_steps; i++) {     
    cpar.w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    for( kk=0; kk < INPUTS.wa_steps; kk++) { 
      cpar.wa = (INPUTS.wa_min + kk*INPUTS.wa_stepsize); 
      for(j=0; j < INPUTS.omm_steps; j++) {
	cpar.omm = INPUTS.omm_min + j*INPUTS.omm_stepsize;    

	prob     = WORKSPACE.extprob3d[i][kk][j];
	probsum += prob;
	diff_w0 = cpar.w0 - WORKSPACE.w0_mean;
	diff_omm = cpar.omm - WORKSPACE.omm_mean;
	cov_w0omm +=  diff_w0 * diff_omm * prob;

	if ( dofit_w0wa ) {
	  diff_wa   = cpar.wa - WORKSPACE.wa_mean;
	  cov_w0wa +=  diff_w0 * diff_wa * prob;
	}

      } // end j
    } // end kk 
  } // end i

  // - - - - - -
  // Compute cov and reduced covariances
  sig_product = (WORKSPACE.w0_sig_marg * WORKSPACE.omm_sig_marg);
  cov_w0omm /= probsum;
  rho_w0omm = cov_w0omm / sig_product;

  if( dofit_w0wa ) { 
    sig_product = (WORKSPACE.w0_sig_marg * WORKSPACE.wa_sig_marg) ;
    cov_w0wa  /= probsum ;
    rho_w0wa = cov_w0wa / sig_product ;
  }

  WORKSPACE.cov_w0omm  = cov_w0omm ;
  WORKSPACE.cov_w0wa   = cov_w0wa ;
  WORKSPACE.rho_w0omm  = rho_w0omm ;
  WORKSPACE.rho_w0wa   = rho_w0wa ;

  printf("# ========================================== \n");
  printf(" Reduced Covariance, w0omm = %.4f \n", WORKSPACE.rho_w0omm);
  if ( dofit_w0wa )
    { printf(" Reduced Covariance, w0wa = %.4f \n", WORKSPACE.rho_w0wa); }

  fflush(stdout);

  return;  
} // end wfit_covariance




// ==========================================x
void wfit_final(void) {

  // Created Oct 2 2021
  // Compute "final" cosmology quantities, 
  // and also sigma_int to get chi2/dof=1.
  // Store in WORKSPACE.

  int    Ndof = WORKSPACE.Ndof ;
  Cosparam cpar;
  double sigmu_int=0.0, sigmu_int1=0.0, sigmu_tmp, sqmusig_tmp ;
  double dif, mindif, muoff_tmp, snchi_tmp, chi2_tmp ;
  double w0_final, wa_final, omm_final ;
  double w0_ran=0.0,  wa_ran=0.0, omm_ran=0.0 ;
  double muoff_final, chi2_final, sigint_binsize ;
  int i;
  char fnam[] = "wfit_final" ;

  // ----------- BEGIN -----------

  printf("# =================================== \n");
  printf(" Extract final cosmology quantities \n");
  fflush(stdout);

  // load final cosmology parameters based on marg or minimize option
  if ( INPUTS.use_marg ) {
    w0_final  = WORKSPACE.w0_mean ;  
    wa_final  = WORKSPACE.wa_mean;   
    omm_final = WORKSPACE.omm_mean ; }
  else {
    w0_final  = WORKSPACE.w0_atchimin ; 
    wa_final  = WORKSPACE.wa_atchimin ; 
    omm_final = WORKSPACE.omm_atchimin ; 
  }

  // store final results before adding random offsets
  cpar.w0  = w0_final;
  cpar.wa  = wa_final;
  cpar.omm = omm_final;
  cpar.ome = 1. - cpar.omm;
  
  // add unknown offset if blind option 
  if ( INPUTS.blind ) {
    srand(48901);
    w0_ran     = floor(1e6*rand()/RAND_MAX);
    wa_ran     = floor(1e6*rand()/RAND_MAX);
    omm_ran    = floor(1e6*rand()/RAND_MAX);
    w0_final   += sin(w0_ran);
    wa_final   += sin(wa_ran);
    omm_final  += sin(omm_ran);
  }

  // store final results and random numbers
  WORKSPACE.w0_final  = w0_final ;
  WORKSPACE.wa_final  = wa_final ; 
  WORKSPACE.omm_final = omm_final ; 

  WORKSPACE.w0_ran  = w0_ran ;
  WORKSPACE.wa_ran  = wa_ran ;
  WORKSPACE.omm_ran = omm_ran ;
    
  // get final chi2 and mu-offset from final parameters.
  get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, INPUTS.sqsnrms,  // inputs
		&muoff_final, &snchi_tmp, &chi2_final );   // return args

  WORKSPACE.chi2_final  = chi2_final;
  WORKSPACE.muoff_final = muoff_final ;

  // - - - -  sigmu_int - - - - 
  WORKSPACE.sigmu_int = 0.0;
  if ( INPUTS.fitnumber == 1 ) { return; } // Nov 24 2021

  sigint_binsize = 0.01;
  printf("\t search for sigint in bins of %.4f \n", sigint_binsize);
  fflush(stdout);
  mindif = 999999.;
  for ( i = 0; i< 20; i++ ) {
    sigmu_tmp   = (double)i * sigint_binsize ;
    sqmusig_tmp = sigmu_tmp * sigmu_tmp ;
    invert_mucovar(sqmusig_tmp);
    get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, sqmusig_tmp,  // inputs
	  	  &muoff_tmp, &snchi_tmp, &chi2_tmp );   // return args
    
    dif = chi2_tmp/(double)Ndof - 1.0 ;
    if ( fabsf(dif) < mindif ) {
      sigmu_int1 = sigmu_tmp ; mindif = dif;
    }
  } // end i loop

  // search again in .001 sigmu_int bins
  sigint_binsize = 0.0002;
  printf("\t search for sigint in bins of %.4f \n", sigint_binsize);
  fflush(stdout);
  mindif = 9999999. ;
  for ( i = -10; i < 10; i++ ) {
    sigmu_tmp   = sigmu_int1 + (double)i * sigint_binsize ;
    sqmusig_tmp = sigmu_tmp * sigmu_tmp ;
    invert_mucovar(sqmusig_tmp);
    
    get_chi2wOM ( cpar.w0, cpar.wa, cpar.omm, sqmusig_tmp,  // inputs
		  &muoff_tmp, &snchi_tmp, &chi2_tmp );   // return args
    
    dif = chi2_tmp/(double)Ndof - 1.0 ;
    if ( fabsf(dif) < mindif ) {
      sigmu_int = sigmu_tmp ; mindif = dif;
    }
  }
 
  WORKSPACE.sigmu_int   = sigmu_int ;
  

  return;
} // end wfit_final

// ==================================
void wfit_FoM(void) {
  // estimate FoM for w0wa model ... need to finish ...
  int Ndof = WORKSPACE.Ndof;
  double extchi_min = WORKSPACE.extchi_min;
  int i, kk, j;
  double extchi, extchi_dif, chi_approx, snchi_tmp, extchi_tmp, muoff_tmp;;
  double sig_product, rho;
  Cosparam cpar;
  char fnam[] = "wfit_FoM" ;
  // --------------BEGIN --------------
  /*
    1. Run a loop over w0, wa, Om
    2. Call getchi2Om function and compute the Chi-sq for each bin
    3. Compute delta_Chi = [Chi-sq] - [Chi-sq_min]
    4. Check if delta_Chi is less than 4.6
    5. Store it in a variable (eg. Chi_Tot) if Yes
    6. Return the Sum of Chi_Tot
   */
  
  WORKSPACE.FoM_final = 0.0 ;

  if ( !INPUTS.dofit_w0wa ) {return ; }

  sig_product = (WORKSPACE.w0_sig_marg * WORKSPACE.wa_sig_marg);


  rho = WORKSPACE.rho_w0wa;
  if(fabs(rho)>=1.){
    sprintf(c1err,"Invalid Rho = %f \n",rho);
    sprintf(c2err,"Check Covariance Calculation ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
}

  sig_product *= sqrt(1.0- rho*rho);

  
  if (sig_product > 0. ) 
    { WORKSPACE.FoM_final = 1.0/sig_product; }
  else 
    { WORKSPACE.FoM_final = -9.0; }

  printf("# ====================================== \n");
  printf(" FOM = %.2f\n", WORKSPACE.FoM_final);
  fflush(stdout);

  return ;
} // emd wfit_FoM

// ==================================
void set_HzFUN_for_wfit(double H0, double OM, double OE, double w0, double wa,
			HzFUN_INFO_DEF *HzFUN ) {

  double COSPAR_LIST[10];  int VBOSE=0;
  // ---------- BEGIN ----------

  COSPAR_LIST[ICOSPAR_HzFUN_H0] = H0 ;
  COSPAR_LIST[ICOSPAR_HzFUN_OM] = OM ;
  COSPAR_LIST[ICOSPAR_HzFUN_OL] = OE ;
  COSPAR_LIST[ICOSPAR_HzFUN_w0] = w0 ;
  COSPAR_LIST[ICOSPAR_HzFUN_wa] = wa ; 
  init_HzFUN_INFO(VBOSE, COSPAR_LIST, "", HzFUN);
    
} // end set_HzFUN_for_wfit

// ==================================
void invert_mucovar(double sqmurms_add) {

  // May 20, 2009 + 04 OCt, 2021
  // add diagonal elements to covariance matrix, then invert.
  //
  int  i, NSN = HD.NSN;
  time_t t0, t1;
  bool check_inverse = (INPUTS.debug_flag == 1000);
  double *MUCOV_ORIG ;

  // ---------------- BEGIN --------------

  printf("\t Invert %d x %d mucov matrix with COV_DIAG += %f \n", 
	 NSN, NSN, sqmurms_add);

  if ( sqmurms_add != 0.0 ) 
    { printf("\t\t xxx WARNING: fix bug and add sqmurms ... \n"); }

  fflush(stdout);

  if ( INPUTS.use_mucov) {
    t0 = time(NULL);

    if ( check_inverse ) {
      int MEMD = NSN * NSN * sizeof(double);
      MUCOV_ORIG = (double*) malloc(MEMD);
      for(i=0; i < NSN*NSN; i++ ) { MUCOV_ORIG[i] = WORKSPACE.MUCOV[i]; }
    }

    invertMatrix( NSN, NSN, WORKSPACE.MUCOV ) ;
    t1 = time(NULL);
    double t_invert = (t1-t0);
    printf("\t Time to invert mucov matrix: %.1f seconds.\n", t_invert);

    if ( check_inverse ) 
      { check_invertMatrix(NSN,MUCOV_ORIG,WORKSPACE.MUCOV); }

  }
  fflush(stdout);

  return ;

} // end of invert_mucovar


// =========================================
double check_invertMatrix(int N, double *COV, double *COVINV ) {

  int i, j, k;
  double val, valinv;
  double prod, absprod, prodmax_offdiag=0.0, prodmax_diag=0.0 ;
  double *ptr_prodmax;
  char fnam[] = "check_invertMatrix" ;

  // ---------- BEGIN -----------

  for(i=0; i < N; i++ ) {
    for(j=0; j < N; j++ ) {
      prod = 0.0 ;
      for(k=0; k < N; k++ ) {
	val    = COV[i*N+k] ; 
	valinv = COVINV[k*N+j] ; 
	prod  += val * valinv;
	//c[i][j]+=a[i][k]*b[k][j];
      } // end k

      if ( i==j ) 
	{ absprod = fabs(prod-1.0); ptr_prodmax = &prodmax_diag; }
      else
	{ absprod = fabs(prod);     ptr_prodmax = &prodmax_offdiag; }
      
      if ( absprod > *ptr_prodmax ) { *ptr_prodmax = absprod; }
    }
  }

  printf("\t max[diag-1  C*Cinv] = %le\n", prodmax_diag);
  printf("\t max[offdiag C*Cinv] = %le\n", prodmax_offdiag);
  fflush(stdout);

  return;
} // end check_invertMatrix

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
  // Oct     2021: refactor using new cov matrix formalism
  // Dec 01  2021: enable speed trick using nsig_chi2_skip (R.Kessler)
  // 
  // Apr 22 2022:
  //   + use dmu_list to avoid redundant log10 calculations in get_DMU_chi2wOM
  //   + implement rz-interpolation option

  bool USE_SPEED_OFFDIAG = INPUTS.USE_SPEED_OFFDIAG ;
  bool USE_SPEED_INTERP  = INPUTS.USE_SPEED_INTERP ;
  int  use_mucov = INPUTS.use_mucov ;
  int  NSN       = HD.NSN;
  int  Ndof      = WORKSPACE.Ndof ;
  double sig_chi2min_naive = WORKSPACE.sig_chi2min_naive;
  double nsig_chi2min_skip = WORKSPACE.nsig_chi2min_skip;  
  double chi_hat_naive     = (double)Ndof;

  double OE, rz, sqmusig, sqmusiginv, Bsum, Csum ;
  double nsig_chi2, chi_hat, chi_tmp ;
  double dmu, dmu0, dmu1, mu_cos  ;
    
  double  chi2_prior = 0.0 ;
  double *rz_list  = (double*) malloc(NSN * sizeof(double) );
  double *dmu_list = (double*) malloc(NSN * sizeof(double) );
  Cosparam cparloc;
  int k, k0, k1, N0, N1, k1min, n_count=0 ;
  
  // rz-interp variables
  int n_logz, iz;
  double z ;

  char fnam[] = "get_chi2wOM";

  // --------- BEGIN --------

  OE = 1.0 - OM ;
  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa ;

  Bsum = Csum = chi_hat = 0.0 ;

  // Apr 2022: check option to interpolate rz(z) [speed trick]
  if ( USE_SPEED_INTERP ) {
    n_logz   = WORKSPACE.n_logz_interp;
    for(iz=0; iz < n_logz; iz++ ) {
      z   = WORKSPACE.z_list_interp[iz];
      rz  = codist(z, &cparloc); 
      mu_cos = get_mu_cos(z,rz);  // theory mu
      WORKSPACE.rz_list_interp[iz]    = rz;
      WORKSPACE.mucos_list_interp[iz] = mu_cos ;
    }
  }

  // Compute diag part first and precompute rz in each z bin to 
  // avoid redundant calculations when using covariance matrix.
  for(k=0; k < NSN; k++ )  { 

    if ( USE_SPEED_INTERP )  { 
      exec_rz_interp(k, &cparloc, &rz, &mu_cos); 
    }
    else { 
      rz     = codist(HD.z[k], &cparloc);
      mu_cos = get_mu_cos(HD.z[k], rz) ;
    }

    rz_list[k]  = rz ;
    dmu_list[k] = HD.mu[k] - mu_cos;

    n_count++ ;
    if ( use_mucov ) {
      sqmusiginv = WORKSPACE.MUCOV[k*(NSN+1)]; 
    }
    else {
      sqmusig     = HD.mu_sqsig[k] + sqmurms_add ;
      sqmusiginv  = 1.0 / sqmusig ;
    }

    dmu         = dmu_list[k] ;
    Bsum       += sqmusiginv * dmu ;       // Eq. A.11 of Goliath 2001
    Csum       += sqmusiginv ;             // Eq. A.12 of Goliath 2001
    chi_hat    += sqmusiginv * dmu*dmu ;

  } // end k

  
  // - - - - - -
  // check for adding off-diagonal terms from cov matrix.
  // If chi_hat(diag) is already > 10 sigma above naive chi2 -> 
  // skip off-diag computation to save time.
  bool do_offdiag = false;
  if ( use_mucov ) {
    if ( USE_SPEED_OFFDIAG ) {
      chi_tmp     = chi_hat - Bsum*Bsum/Csum ;
      nsig_chi2  = (chi_tmp - chi_hat_naive ) / sig_chi2min_naive ;
      do_offdiag = nsig_chi2 < nsig_chi2min_skip ;

      /* xxxx
      if ( fabs(w0+1.0) < 0.06  &&  fabs(OM-0.3)<0.02 ) {
	printf(" xxx %s: ------------------------------ \n", fnam);
	printf(" xxx %s: w0=%f  wa=%f  OM=%f \n",
	       fnam, w0, wa, OM);
	printf(" xxx %s: nsig=%.1f chi(hat,naive)=%.1f,%.1f  "
	       "sig_naive=%.1f \n",
	       fnam, nsig_chi2, chi_tmp, chi_hat_naive, sig_chi2min_naive);
	debugexit(fnam); // xxxx remove
	} xxxxx*/

    }
    else {
      do_offdiag = true ; 
    }
  }

  // add off-diag elements if using cov matrix
  if ( do_offdiag ) {
    for ( k0=0; k0 < NSN-1; k0++) {
      k1min = k0 + 1;
      for ( k1=k1min; k1 < NSN; k1++)  {
	k = k0*NSN + k1;
	sqmusiginv = WORKSPACE.MUCOV[k]; // Inverse of the matrix 

	dmu0     = dmu_list[k0];
	dmu1     = dmu_list[k1];

	chi_tmp  = (sqmusiginv * dmu0 * dmu1 );

	Bsum    += sqmusiginv*(dmu0+dmu1); // Eq. A.11 of Goliath 2001  
	Csum    += (2.0*sqmusiginv);       // Eq. A.12 of Goliath 2001
	chi_hat += (2.0*chi_tmp);

	/* xxxx
	Bsum    += (sqmusiginv * dmu0) ;   // Eq. A.11 of Goliath 2001  
	Csum    += sqmusiginv ;          // Eq. A.12 of Goliath 2001
	chi_hat += chi_tmp ;

	Bsum    += (sqmusiginv * dmu1) ;
	Csum    += sqmusiginv ;
	chi_hat += chi_tmp ; 
	xxxx */

	n_count++ ;

      } // end k1
    } // end k0
  }  // end do_offdiag

  //  - - - - - -
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
  if ( INPUTS.use_bao ) {
    chi2_prior += chi2_bao_prior(&cparloc);
  } 
  else {
    // Gaussian Omega_m prior
    double nsig;
    nsig = (OM-INPUTS.omm_prior)/INPUTS.omm_prior_sig ;
    chi2_prior += nsig*nsig;
  }

  // CMB prior
  if ( INPUTS.use_cmb ) {
    chi2_prior += chi2_cmb_prior(&cparloc);
  }

  *chi2tot += chi2_prior;

  free(rz_list);
  free(dmu_list);

  return ;

}  // end of get_chi2wOM


// ======================================
double chi2_cmb_prior(Cosparam *cpar) {

  // Created Oct 2021
  // Return chi2 for cmb prior

  double R_prior    = CMB_PRIOR.R;
  double sigR_prior = CMB_PRIOR.sigR;
  double a          = CMB_PRIOR.a ;

  double rz, nsig, R_calc, chi2=0.0 ;
  HzFUN_INFO_DEF HzFUN;

  // ------------- BEGIN --------------

  set_HzFUN_for_wfit(ONE, cpar->omm, cpar->ome, cpar->w0, cpar->wa, &HzFUN);
  rz         = Hainv_integral ( a, ONE, &HzFUN ) / LIGHT_km;
  R_calc     = sqrt(cpar->omm) * rz ;
  nsig       = (R_calc - R_prior) / sigR_prior ;
  chi2       = nsig*nsig;

  return chi2;

} // end chi2_cmb_prior

// ======================================
double chi2_bao_prior(Cosparam *cpar) {

  // Created Oct 2021
  // Return chi2 for bao prior.
  // Mar 09 2022 RK - fix bug to allow -bao_sim

  double OM        = cpar->omm ;
  double rz, tmp1, tmp2, a, z, nsig, E, chi2 = 0.0 ;
  double DM_rd_measured,  DH_rd_measured,  DMvar, DHvar;
  double rd_model,dif_M, dif_H, DM_rd_var, DH_rd_var,DH_rd_model,DM_rd_model;
  int    iz;
  bool   DEBUG_TINYERR = (INPUTS.debug_flag == 309 );
  char   fnam[] = "chi2_bao_prior" ;

  // ------------ BEGIN -------------

  if ( BAO_PRIOR.use_sdss4 ) {

    // Alam 2020, SDSS-IV;  Ayan Mitra Nov, 2021    
    double *z_sdss4 = BAO_PRIOR.z_sdss4;
    double *DMrd    = BAO_PRIOR.DMrd_sdss4 ;
    double *DHrd    = BAO_PRIOR.DHrd_sdss4 ;
    double *sigDMrd = BAO_PRIOR.sigDMrd_sdss4;
    double *sigDHrd = BAO_PRIOR.sigDHrd_sdss4;

    rd_model      = rd_bao_prior( cpar);
    for(iz=0; iz < NZBIN_BAO_SDSS4; iz++ ) {
      
      // Measured
      z               = z_sdss4[iz];
      DH_rd_measured  = DHrd[iz];
      DM_rd_measured  = DMrd[iz]; 
      DM_rd_var       = (sigDMrd[iz] * sigDMrd[iz]); 
      DH_rd_var       = (sigDHrd[iz] * sigDHrd[iz]);

      if ( DEBUG_TINYERR ) {
	DM_rd_var *= 0.1 ;
	DH_rd_var *= 0.1 ;
      }

      // Model
      DH_rd_model   = DH_bao_prior(z, cpar)/rd_model ;
      DM_rd_model   = DM_bao_prior(z, cpar)/rd_model ;

      // measured - model
      dif_H         = (DH_rd_measured - DH_rd_model);
      dif_M   	    = (DM_rd_measured - DM_rd_model);
     
      chi2 += (dif_H * dif_H) / DH_rd_var;
      chi2 += (dif_M * dif_M) / DM_rd_var; 
    }
  }
  else if ( BAO_PRIOR.use_sdss ) {
    // Eisen 2006
    double z_sdss    = BAO_PRIOR.z_sdss;
    double a_sdss    = BAO_PRIOR.a_sdss;
    double siga_sdss = BAO_PRIOR.siga_sdss;

    rz   = codist(z_sdss, cpar);
    E    = EofZ(z_sdss, cpar) ;
    tmp1 = pow( E, NEGTHIRD) ;
    tmp2 = pow( (rz/z_sdss),   TWOTHIRD );
    a    = sqrt(OM) * tmp1 * tmp2 ;
    nsig = (a - a_sdss) / siga_sdss ;
    chi2 = nsig * nsig ;
  }
  return chi2;

} // end chi2_bao_prior

// ===============================================
double get_DMU_chi2wOM(double z, double rz, double mu)  {

  // Created oct, 2021. Mitra, Kessler.
  // For the chi sq. function to evaluate the Hubble residual.
  // Inputs :
  //   z  : redshift
  //   rz : codist
  //   mu : HD.mu[k]
  // Output :
  //   Function returns mu_obs-mu_theory 
  
  //  double z      = HD.z[k];
  //  double mu_obs = HD.mu[k] ;
  double H0     = H0_SALT2 ;
  double ld_cos, mu_cos, dmu;
  
  ld_cos      = (1.0 + z) *  rz * c_light / H0;
  mu_cos      = 5.0 * log10(ld_cos) + 25. ;
  dmu         = mu_cos - mu ;
  return dmu ;

}  // end get_DMU_chi2wOM

// ===============================================
double get_mu_cos(double z, double rz)  {

  // Created oct, 2021. Mitra, Kessler.
  // For the chi sq. function to evaluate the Hubble residual.
  // Inputs :
  //   z  : redshift
  //   rz : codist
  // Output :
  //   Function returns mu_cos
  
  //  double z      = HD.z[k];
  //  double mu_obs = HD.mu[k] ;
  double H0     = H0_SALT2 ;
  double ld_cos, mu_cos ;
  
  ld_cos      = (1.0 + z) *  rz * c_light / H0;
  mu_cos      = 5.0 * log10(ld_cos) + 25. ;
  return mu_cos ;

}  // end get_mu_cos


// ==============================================
double get_minwOM( double *w0_atchimin, double *wa_atchimin, 
		   double *omm_atchimin ) {

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
    ,snchi_tmp, extchi_tmp, snchi_min, extchi_min, muoff_tmp
    ;


  double w0step_tmp  = INPUTS.w0_stepsize / 10.0; 
  double wastep_tmp  = INPUTS.wa_stepsize / 10.0;
  double omstep_tmp  = INPUTS.omm_stepsize/ 10.0; 
  int    nbw0, nbwa, nbm,i, j, kk, imin=0, kmin =0, jmin=0 ;
  double nb_factor = 4.0;

  int  LDMP = 0;
  char fnam[] = "get_minwOM";

  // ---------- BEGIN ------------

  snchi_min = 1.e20;  extchi_min = 1.e20 ;

  nbw0 = (int)(nb_factor*INPUTS.w0_stepsize/w0step_tmp) ;
  nbwa = (int)(nb_factor*INPUTS.wa_stepsize/wastep_tmp) ;
  nbm  = (int)(nb_factor*INPUTS.omm_stepsize/omstep_tmp) ;

  printf("   Minimize with refined %sstep=%6.4f (%d bins, %s=%7.4f) \n", 
	 varname_w, w0step_tmp, 2*nbw0, varname_w, *w0_atchimin);
  if(INPUTS.dofit_w0wa){
    printf("   Minimize with refined wastep=%6.4f (%d bins, wa=%7.4f) \n",
	   wastep_tmp, 2*nbwa, *wa_atchimin);
  }
  printf("   Minimize with refined OMstep=%6.4f (%d bins, OMcen=%6.4f) \n", 
	 omstep_tmp, 2*nbm, *omm_atchimin );
  printf("---------------------------------------\n");

  w0cen_tmp  = *w0_atchimin ;
  wacen_tmp  = *wa_atchimin ;
  omcen_tmp  = *omm_atchimin ;

  for ( i = -nbw0; i <= nbw0; i++ ) {
    w0_tmp = w0cen_tmp + (double)i*w0step_tmp ;
    // break up w into w0,wa AM

    for ( kk = -nbwa; kk <= nbwa; kk++ ) {
      wa_tmp = wacen_tmp + (double)kk*wastep_tmp ;

      for ( j = -nbm; j < nbm; j++ ) {
	om_tmp = omcen_tmp + (double)j*omstep_tmp ;

	if ( LDMP ) {
	  fflush(stdout);
	}

	get_chi2wOM(w0_tmp, wa_tmp, om_tmp, INPUTS.sqsnrms, 
		  &muoff_tmp, &snchi_tmp, &extchi_tmp );

	if ( extchi_tmp < extchi_min ) 
	  { extchi_min = extchi_tmp ;  imin=i; kmin = kk; jmin=j; }
      }
    } // end j
  } // end i

  // - - - - - -
  // change input values with final w,OM
  *w0_atchimin  = w0cen_tmp + (double)imin * w0step_tmp ;
  *wa_atchimin  = wacen_tmp + (double)kmin * wastep_tmp ;
  *omm_atchimin = omcen_tmp + (double)jmin * omstep_tmp ;
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


double codist(double z, Cosparam *cptr) {
  // Returns dimensionless comoving distance.
  double zero = 0.0 ;
  return simpint(one_over_EofZ, zero, z, cptr);
}

double Eainv_integral(double amin, double amax, Cosparam *cptr) {
  // Created Oct 2021 by R.Kessler
  // return integral of 1/E..
  // This is the analog of codist, but integrating over a instead of z.
  // Note that the one_over_EofA function includes 1/a^2 jacobian factor.
  return simpint(one_over_EofA, amin, amax, cptr);
}

double one_over_EofZ(double z, Cosparam *cptr){
  // This is actually the function that we pass to the integrator.
  return 1./EofZ(z, cptr);
}
double one_over_EofA(double a, Cosparam *cptr) {
  double Einv     = 1./EofA(a, cptr);
  double Jacobian = 1.0/(a*a);  // dz = da/a^2, Jac=1/a^2
  return Einv * Jacobian;
}

double EofZ(double z, Cosparam *cptr){

  // Oct 31 2021 R.Kessler
  //  Implement omr (radiation) erm such that setting omr=0.0 reduce
  //  to previous code.
  //  omr includes photons+neutrinos for early-universe integrals.
  //  This is inaccurate at late times when neutrinos become non-relativistic,
  //  but the late time omr contribution is so small that the approximation
  //  is insignificant.

  double E;
  double arg, argr, omm, omk, ome, w0, wa;
  double z1       = 1.0 + z;
  double z1_pow2  = z1*z1;
  double z1_pow3  = z1_pow2 * z1;
  double z1_pow4;
  double omr      = 0.9E-4 ; // photon+neutrino, Oct 31 2021, RK

  omm = cptr->omm;
  ome = cptr->ome;
  omk = 1.0 - omm - ome;
  w0  = cptr->w0;
  wa  = cptr->wa;

  if ( omr > 0.0 ) {
    // omm + ome + omr = 1 => preserve flat universe 
    omm *= (1.0 - omr);
    ome *= (1.0 - omr);
    z1_pow4  = z1_pow2 * z1_pow2;
    argr     = omr * z1_pow4;
  }

  arg = 
    omm * z1_pow3 + 
    omk * z1_pow2 + 
    ome * pow(z1, 3.*(1+w0+wa))*exp(-3.0*wa*(z/z1)) ;

  if ( omr > 0.0 )  { arg += (omr * z1_pow4); }

  E = sqrt(arg);

  return E;
}

double EofA(double a, Cosparam *cptr) {
  double z = 1.0/a - 1.0 ;
  double E = EofZ(z, cptr);
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
  double H0 = H0_SALT2;
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
  exit(0); 


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
  for ( i=0; i < HD.NSN; i++ ) {
    if ( strcmp(HD.cid[i],cid) == 0 ) return i;
  }

  return -1;

} // end of cidindex


// ************************************
void WRITE_OUTPUT_DRIVER(void) {

  char fnam[] = "WRITE_OUTPUT_DRIVER" ;

  // ------------ BEGIN ------------

  printf("# ==================================== \n");
  printf(" %s\n", fnam);

  // primary output is the cosmo params
  write_output_cospar();

  // Print residuals (not used for long time; beware)
  write_output_resid(INPUTS.fitnumber);
 
  // write chi2 & prob maps to fits files;
  // not used for VERY long time ... beware.
  if ( INPUTS.fitsflag ) { write_output_fits(); }

  return;
} // end WRITE_OUTPUT_DRIVER

// ********************************
void write_output_cospar(void) {

  // Created Aug 15 2020
  // Refactor Aug 2021 to avoid too much redundant code.
  //
  // format_cospar = 1 : legacy csv format
  // format_cospar = 2 : YAML format

  int  use_marg  = INPUTS.use_marg;
  int  dofit_w0wa= INPUTS.dofit_w0wa ;
  int  format    = INPUTS.format_cospar;
  int  blind     = INPUTS.blind ;
  char *outFile  = INPUTS.outFile_cospar;
  double dt_fit  = (double)(t_end_fit  - t_start)/ 60.0 ;

#define MXVAR_WRITE 20
  int    ivar, NVAR = 0;
  FILE *fp ;
  char   VALUES_LIST[MXVAR_WRITE][20];
  char   VARNAMES_LIST[MXVAR_WRITE][20], LINE_STRING[200] ;
  char   ckey[40], cval[40];
  char sep[] = " " ;

  // ----------- BEGIN -------------

  if ( strlen(outFile) == 0){
    sprintf(outFile, "%s.cospar", INPUTS.infile);
  } 
    
  printf("  Write cosmo params to %s \n", outFile);
  fp = fopen(outFile, "w");
  if (fp == NULL){
    printf("ERROR: couldn't open %s\n", outFile);
    exit(EXIT_ERRCODE_wfit);
  }

  // - - - -
  // define variables to write out based on use_marg and dofit_w0wa flags.

  sprintf(VARNAMES_LIST[NVAR],"%s", varname_w); 
  sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.w0_final ) ;
  NVAR++ ;

  if ( use_marg ) {
    sprintf(VARNAMES_LIST[NVAR],"%ssig_marg", varname_w); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.w0_sig_marg ) ;
    NVAR++ ;
  }
  else {
    sprintf(VARNAMES_LIST[NVAR],"%ssig_lo", varname_w); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.w0_sig_lower ) ;
    NVAR++ ;
    sprintf(VARNAMES_LIST[NVAR],"%ssig_up", varname_w); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.w0_sig_upper ) ;
    NVAR++ ;
  }

  // - - - 
  if ( dofit_w0wa ) {
    sprintf(VARNAMES_LIST[NVAR],"%s", varname_wa); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.wa_final ) ;
    NVAR++ ;
    if ( use_marg ) {
      sprintf(VARNAMES_LIST[NVAR],"%ssig_marg", varname_wa); 
      sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.wa_sig_marg ) ;
      NVAR++ ;
    }
    else {
      sprintf(VARNAMES_LIST[NVAR],"%ssig_lo", varname_wa); 
      sprintf(VALUES_LIST[NVAR], "%7.4f",  WORKSPACE.wa_sig_lower ) ;
      NVAR++ ;
      sprintf(VARNAMES_LIST[NVAR],"%ssig_up", varname_wa); 
      sprintf(VALUES_LIST[NVAR], "%7.4f",  WORKSPACE.wa_sig_upper ) ;
      NVAR++ ;
    }
  } // end dofit_w0wa


  // - - - 

  sprintf(VARNAMES_LIST[NVAR],"%s", varname_omm); 
  sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.omm_final ) ;
  NVAR++ ;
  if ( use_marg ) {
    sprintf(VARNAMES_LIST[NVAR],"%ssig_marg", varname_omm); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.omm_sig_marg ) ;
    NVAR++ ;
  }
  else {
    // WARNING: we don't have omm_sig_upper/lower, so just
    // write omm_sig_marg
    sprintf(VARNAMES_LIST[NVAR],"%ssig_marg", varname_omm); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.omm_sig_marg ) ;
    NVAR++ ;
  }

  if ( dofit_w0wa ) {
    sprintf(VARNAMES_LIST[NVAR],"FoM" );
    sprintf(VALUES_LIST[NVAR], "%5.1f",  WORKSPACE.FoM_final ) ;
    NVAR++ ;

    sprintf(VARNAMES_LIST[NVAR],"Rho" );
    sprintf(VALUES_LIST[NVAR], "%6.3f",  WORKSPACE.rho_w0wa ) ;
    NVAR++;
  }




  sprintf(VARNAMES_LIST[NVAR],"chi2" );
  sprintf(VALUES_LIST[NVAR], "%.1f",  WORKSPACE.chi2_final ) ;
  NVAR++ ;

  sprintf(VARNAMES_LIST[NVAR],"Ndof" );
  sprintf(VALUES_LIST[NVAR], "%5d",  WORKSPACE.Ndof ) ;
  NVAR++ ;

  sprintf(VARNAMES_LIST[NVAR],"sigint" );
  sprintf(VALUES_LIST[NVAR], "%6.3f",  WORKSPACE.sigmu_int ) ;
  NVAR++ ;

  sprintf(VARNAMES_LIST[NVAR],"%sran", varname_w );
  sprintf(VALUES_LIST[NVAR], "%d",  (int)WORKSPACE.w0_ran ) ;
  NVAR++ ;

  if ( dofit_w0wa ) {
    sprintf(VARNAMES_LIST[NVAR],"%sran", varname_wa );
    sprintf(VALUES_LIST[NVAR], "%d",  (int)WORKSPACE.wa_ran ) ;
    NVAR++ ;
  }

  sprintf(VARNAMES_LIST[NVAR],"%sran", varname_omm );
  sprintf(VALUES_LIST[NVAR], "%d",  (int)WORKSPACE.omm_ran ) ;
  NVAR++ ;



  sprintf(VARNAMES_LIST[NVAR],"label" );
  sprintf(VALUES_LIST[NVAR], "%s",  INPUTS.label_cospar ) ;
  NVAR++ ;
  
  // - - - - - - -
  // List of variable names and values are stored;
  // now write them out based on format option;
  // csv-like or YAML


  // legacy format in csv-like format with default sep=" "
  // Always write this to stdout; write to file if format==1.

  LINE_STRING[0] = 0;
  sprintf(LINE_STRING,"# ");
  for(ivar=0; ivar < NVAR; ivar++ ) {
    strcat(LINE_STRING, VARNAMES_LIST[ivar] );
    strcat(LINE_STRING, sep );
  }
  if ( format == 1 ) { fprintf(fp,"%s\n", LINE_STRING); }  // write VARNAMES
  printf("COSPAR:  %s\n", &LINE_STRING[1]);  // write VARNAMES
    
  // now print values
  LINE_STRING[0] = 0;
  for(ivar=0; ivar < NVAR; ivar++ ) {
    strcat(LINE_STRING, VALUES_LIST[ivar] );
    strcat(LINE_STRING, sep );
  }
  if ( format == 1 ) {fprintf(fp,"%s\n", LINE_STRING); } // write values
  printf("COSPAR:  %s\n", LINE_STRING); // write values
 

  // - - --  - check YAML format - - - - 
  if ( format == 2 ) {
    // YAML format
    for(ivar=0; ivar < NVAR; ivar++ ) {
      sprintf(ckey, "%s:",  VARNAMES_LIST[ivar]);
      sprintf(cval, "%s",   VALUES_LIST[ivar]);
      fprintf(fp, "%-14s  %s \n", ckey, cval);
    }
    fprintf(fp,"NWARNINGS:      %d \n", WORKSPACE.NWARN);
    fprintf(fp,"ABORT_IF_ZERO:  %d   # same as Ndof \n", WORKSPACE.Ndof);
    fprintf(fp,"CPU_MINUTES:    %.2f \n", dt_fit);
    fprintf(fp,"BLIND:          %d \n", blind);
  }

  fclose(fp);

  return ;

} // end write_output_cospar


// ********************************
void write_output_resid(int fitnum) {

  Cosparam cpar;
  FILE *fpresid;
  char *tempfilename1 = (char*) malloc(MEMC_FILENAME);
  char *residfile     = (char*) malloc(MEMC_FILENAME);
  char fnam[] = "write_output_resid" ;

  // ----------- BEGIN -------------

  if ( IGNOREFILE(INPUTS.outFile_resid) ) { return; }

  getname(INPUTS.outFile_resid, tempfilename1, fitnum);
  strcpy(residfile, tempfilename1);

  printf("Write MU-residuals to %s  \n",
	 residfile );
  fpresid=fopen(residfile, "w");
      
  if (fpresid == NULL) {
    printf("ERROR: couldn't open %s\n",residfile);
    exit(EXIT_ERRCODE_wfit);
  }
    
  int i;
  double rz, ld_cos, mu_cos, mu_dif, sqmusig_tmp, musig_tmp;
  double H0 = H0_SALT2;
  fprintf(fpresid,"z, mu_dif,  mu_sig,  tel_id,  CID \n"); // csv format
  for (i=0; i<HD.NSN; i++){
    rz     = codist( HD.z[i], &cpar) ;
    ld_cos = ( 1.0 + HD.z[i]) *  rz * c_light / H0;
    mu_cos =  5.*log10(ld_cos) + 25. ;
    mu_dif =  mu_cos - HD.mu[i];
    sqmusig_tmp  =  HD.mu_sqsig[i] + INPUTS.sqsnrms ;
    musig_tmp    = sqrt(sqmusig_tmp);
	
    fprintf(fpresid,"%7.4f %7.4f %7.4f  %6s \n",
	    HD.z[i], mu_dif, musig_tmp, HD.cid[i]);
  }
  fclose(fpresid);

  return ;

} // end write_output_resid

// ********************************
void write_output_contour(void) {
  // ----------- BEGIN -------------
  return ;
} // end write_output_contour

// =======================================
void write_output_fits(void) {

  // Created Oct 2 2021
  // Write to fits files.
  // Ancient code moved out of main.

  FILE *FILEPTR_TOT, *FILEPTR_SN ;
  fitsfile *fptr;
  char *snfitsname  = (char*) malloc(MEMC_FILENAME);
  char *ommfitsname = (char*) malloc(MEMC_FILENAME);
  char *snfitsstr   = (char*) malloc(MEMC_FILENAME);
  char *ommfitsstr  = (char*) malloc(MEMC_FILENAME);
  char *cosparfile  = (char*) malloc(MEMC_FILENAME);
  char *residfile   = (char*) malloc(MEMC_FILENAME);
  char *chi2gridfile= (char*) malloc(MEMC_FILENAME);
  char *chi2gridfile_SN = (char*) malloc(MEMC_FILENAME);
  char *tempfilename1   = (char*) malloc(MEMC_FILENAME);
  char *tempfilename2   = (char*) malloc(MEMC_FILENAME);

  long naxis=2;
  long fpixel[2], nelements, naxes[2];
  int status, i, kk, j ;

  // ----------- BEGIN -----------

  naxes[0] = INPUTS.omm_steps;
  naxes[1] = INPUTS.w0_steps;
  status = 0;
  nelements = naxes[0] * naxes[1];    /* number of pixels to write */
  fpixel[0]=1;			/* first pixel to write */
  fpixel[1]=1;
      
  int writechi=0; // ??
  if (writechi) {
    /* subtract off minchi from each, use prob arrays for storage */
    for(i=0; i < INPUTS.w0_steps; i++){
      for(kk=0; kk < INPUTS.wa_steps;kk++){
	for(j=0; j < INPUTS.omm_steps; j++){
	  WORKSPACE.snprob3d[j][kk][i] = 
	    WORKSPACE.snchi3d[j][kk][i] - WORKSPACE.snchi_min;
	  WORKSPACE.extprob3d[j][kk][i] = 
	    WORKSPACE.extchi3d[j][kk][i] - WORKSPACE.extchi_min;
	}
      }
    }
      
    strcpy(snfitsname,INPUTS.infile);
    strcat(snfitsname,"_snchi.fits");

    getname(snfitsname, tempfilename1, INPUTS.fitnumber);

    strcpy(ommfitsname,INPUTS.infile);
    strcat(ommfitsname,"_extchi.fits");

    getname(ommfitsname, tempfilename2, INPUTS.fitnumber);

    printf("writing out chi2 distributions:\n");
    printf("  SNe only:      %s\n",tempfilename1);
    printf("  SNe + Omega_m: %s\n",tempfilename2);

  } else {
    /* just write out probabilities */
    strcpy(snfitsname,INPUTS.infile);
    strcat(snfitsname,"_snprob.fits");

    getname(snfitsname, tempfilename1, INPUTS.fitnumber);
    
    strcpy(ommfitsname,INPUTS.infile);
    strcat(ommfitsname,"_extprob.fits");
    
    getname(ommfitsname, tempfilename2, INPUTS.fitnumber);
    
    printf("writing out prob distributions:\n");
    printf("  SNe only:      %s\n",tempfilename1);
    printf("  SNe + Omega_m: %s\n",tempfilename2);
  }

  // - - - -
  /* prepend with "!" so cfitsio clobbers old files */
  sprintf(snfitsstr,"!%s",tempfilename1);
  sprintf(ommfitsstr,"!%s",tempfilename2);
    
  /* Write out SN-only distribution */
  if (fits_create_file(&fptr, snfitsstr, &status))  
    { printerror( status ); }

  if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status)) 
    { printerror( status ); }

  if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, 
		     WORKSPACE.snprob, &status))
    { printerror( status ); }

  if ( fits_close_file(fptr, &status))                /* close the file */
    { printerror( status ); }           
    
  /* Write out SN+Omega_m distribution */
  if (fits_create_file(&fptr, ommfitsstr, &status))  
    { printerror( status ); }

  if (fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status))
    { printerror( status ); }

  if (fits_write_pix(fptr, TDOUBLE, fpixel, nelements, 
		     WORKSPACE.extprob, &status))
    { printerror( status ); }

  if ( fits_close_file(fptr, &status) )                /* close the file */
    { printerror( status ); }           

  // 6/10/2008 RK write chi2 to text file (like Andy's)
  
  if ( strlen(INPUTS.outFile_chi2grid) == 0 ) {
    sprintf(INPUTS.outFile_chi2grid,"%s", INPUTS.infile);
  }
  strcpy(chi2gridfile, INPUTS.outFile_chi2grid); 
  strcat(chi2gridfile,".chi2grid");
  sprintf(chi2gridfile_SN,"%s_SN", chi2gridfile);
      
  getname(chi2gridfile,    tempfilename1, INPUTS.fitnumber);
  getname(chi2gridfile_SN, tempfilename2, INPUTS.fitnumber);
      
  printf("\n");
  printf("\t CHI2GRID(total)   dump to: '%s' \n", tempfilename1 );
  printf("\t CHI2GRID(SN-only) dump to: '%s' \n", tempfilename2 );
  printf("\t CHI2GRID format : i_OM  i_w  OM  w  chi2 \n" );

  int irow=0; char CVARDEF[100];
  double w0, omm, chi2tot, chi2sn ;
  FILEPTR_TOT = fopen(tempfilename1, "wt");
  FILEPTR_SN  = fopen(tempfilename2, "wt");

  sprintf(CVARDEF,"NVAR: 6 \nVARNAMES: ROW iM iw OM w chi2\n");
  fprintf(FILEPTR_TOT,"%s", CVARDEF);
  fprintf(FILEPTR_SN, "%s", CVARDEF);

  for(i=0; i < INPUTS.w0_steps; i++){
    for(j=0; j < INPUTS.omm_steps; j++){
      w0      = INPUTS.w0_min  + i*INPUTS.w0_stepsize;
      omm     = INPUTS.omm_min + j*INPUTS.omm_stepsize;
      chi2tot = WORKSPACE.extchi[i*INPUTS.omm_steps+j] ;
      chi2sn  = WORKSPACE.snchi[i*INPUTS.omm_steps+j] ;

      irow++ ;
      fprintf( FILEPTR_TOT, "ROW: %5d  %5d %5d  %8.4f  %8.4f  %14.6f \n",
	       irow, j, i, omm, w0, chi2tot  );

      fprintf( FILEPTR_SN, "ROW: %5d  %5d %5d  %8.4f  %8.4f  %14.6f \n",
	       irow, j, i, omm, w0, chi2sn  );
      
    }
  }
  fclose(FILEPTR_TOT);
  fclose(FILEPTR_SN);

  return;
} // end write_output_fits


// ==========================
void CPU_summary(void) {

  double dt_init = (double)(t_end_init - t_start)   / 60.0 ;
  double dt_fit  = (double)(t_end_fit  - t_end_init)/ 60.0 ;

  // ----------- BEGIN -----------

  printf("# =================================== \n");
  printf(" CPU Summary \n");
  printf("\t init/fit: %.2f / %.2f minutes \n",	 dt_init, dt_fit);
  fflush(stdout);

  return;
} // end CPU_summary
