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
    + BAO prior from DESI_2024 or SDSS4_2020

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
     Dai      2022 (2212.06879)  SALT3 training syst
     Mitra    2022 (2210.07560)  SNIa cosmology with photo-z/PLASTICC
     Armstron 2023 (2307.13862)  contraint validation
     Kessler  2023 (2306.05819)  Binning-is-sinning redemption     
     Qu       2023 (2307.13696)  Biases from host mis-match
     Vincenzi 2024 (2401.02945)  DES-SN5YR analysis & syst


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
   + disable computation of sigint unless refit option is used.
     (for large NSN, calc was slow because of matrix inversion each step)

 Dec 01 2021 RK :
   in get_chi2_fit(), skip off-diag computation if chi2(diag) is 
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

 July 19 2022 RK
   + new flags -blind_seed and -blind_auto
   + -blind_auto reads ISDATA_REAL flag from HD file and sets
     blinding on for real data, or off for sim data.

 July 26 2022 RK - implement outfile_chi2grid ... not well tested.

 Sep 12 2022 RK - abort if ommin<0 and -cmb
                  in EofZ, if arg < 0.01; arg=0.01 (to allow om<0)

 Sep 27 2022 RK - write rho_wom
 Dec 6 2022 RK - remove obsolete fitswrite option

 Jan 31 2023 RK 
    - enable -outfile_resid feature 
    - in write_output_chi2grid, add blind offsets

 Feb 6 2023 - fix sntools_output.c to select 1st varname on list
              instead of last. E.g., for 'zHD zCMB', use zHD instead of
              zCMB.

 Mar 13 2023: add inputs w0_sim and wa_sim to alter cmb prior in the same
              way as om_sim input.

 Mar 14 2023: begin prep for HDIBC method:
                 Hubble Diagram Interpolated with Bias-cor Cosmology

 Mar 29 2023:
   +  "-debug_flag 91"  to dump cospar comparison between wfit and sim
   +  -muerr_ideal <value> to replace mu with mu_true + Gauss(0,muerr_ideal).

 Aug 03 2023: 
   + use print_cputime(..) utility to grep stdout for standard CPUTIME string.

 Oct 17 2023: 
   + new option ranseed_Rcmb to fluctuate Rcmb value 
     (intended for error validation on many sim-data samples)

 Dec 3 2023: replace fixed muerr_ideal with polyFun(z); 
             e.g., muerr_ideal = .11,0.01,0.07

 Apr 23 2024: Update to more general BAO prior for SDSS4 or DESI based on these args:
      -bao DESI_2024
      -bao SDSS4_2020
      -bao_sim DESI_2024   (use sim cosmology for means + DESI cov)
      -bao_sim SDSS4_2020  (use sim cosmology for means + SDSS4 cov)
   that read z-dependent means and cov from $SNDATA_ROOT/models/BAO

 Apr 25 2024:
   + include optional Gauss prior on rd (to mimic Camilari et al)
   + define separate chi2_prior variables for OM, CMB, BOA, rd  
     and print each chi2-prior contribution separately at end of fit;
     see get_chi2_priors(...)
  
  Apr 26 2024:
    + refactor wfit_uncertainty to simplify code.
    + replace default std error with average of lower+upper at 68% confidence;
      this new default error is larger and FoM will therefore be smaller.
      Revert back to std errors with argument "-sig_type std"
    + fix lower/upper 68% confidence calc using interpolation to 
      avoid original nearest bin approximation and round-off error on error. 
    + add -mucovtot_inv_file option to read already inverted cov matrix.
      This goes with create_covariancy.py update to write covtot_inv_[nnn].txt

  Oct 23 2024:
    write mu_obs in resid file, where mu_obs is from HDIBC method.

  Nov 25 2024:
    +  use long long int for size of covmat malloc to handle
       matrices with more than few GB size.
    +  when reading cov, print update every 2 million elements

 Dec 7 2024
   + mu_obs -= B/C in residual output file 

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

// bit-mask options for speed_flag_chio2
#define SPEED_MASK_OFFDIAG 1  // skip off-diag calc if chi2(diag)>threshold
#define SPEED_MASK_INTERP  2  // interplate r(z) and mu_cos(z)
#define SPEED_FLAG_CHI2_DEFAULT  SPEED_MASK_OFFDIAG + SPEED_MASK_INTERP

#define PROBSUM_1SIGMA  0.683

// Define variable names to read in hubble diagram file.
// VARLIST_DEFAULT_XXX means that any variable is valid for XXX.
#define  VARLIST_DEFAULT_CID     "CID:C*20 ROW:C*20"
#define  VARLIST_DEFAULT_MUERR   "MUERR:D DLMAGERR:D MUDIFERR:D"
#define  VARLIST_DEFAULT_MU      "MU:D DLMAG:D MUDIF:D"
#define  VARLIST_DEFAULT_zHD     "zHD:D zCMB:D z:D Z:D"
#define  VARLIST_DEFAULT_zHDERR  "zHDERR:D zCMBERR:D zERR:D ZERR:D"
#define  VARLIST_DEFAULT_MUREF   "MUREF:D"
#define  VARLIST_DEFAULT_NFIT    "NFIT:I"

#define  KEYNAME_ISDATA_REAL     "ISDATA_REAL:"

// define INFO_YML file name and sim-keys to read cosmo params for biasCor sim
#define INFO_YML_FILENAME     "INFO.YML"  // produced by create_covariance
#define NSIMKEY_COSPAR         5
#define SIMKEY_MUSHIFT        "MUSHIFT:"
#define SIMKEY_OMEGA_LAMBDA   "OMEGA_LAMBDA:"
#define SIMKEY_OMEGA_MATTER   "OMEGA_MATTER:"
#define SIMKEY_w0_LAMBDA      "w0_LAMBDA:"
#define SIMKEY_wa_LAMBDA      "wa_LAMBDA:"


// Define default grid for each fitted parameter.
// They can be modified via command line args; for help type "wfit.exe"
#define DEFAULT_omm_steps  81   // number of steps
#define DEFAULT_omm_min    0.0
#define DEFAULT_omm_max    1.0

#define DEFAULT_w0_steps   201
#define DEFAULT_w0_min    -2.0
#define DEFAULT_w0_max    +0.0

#define DEFAULT_wa_steps   101
#define DEFAULT_wa_min    -4.0
#define DEFAULT_wa_max    +4.0


#define OPT_RD_CALC     0  // 0=Planck; 1=compute with constant c_s, 2=c_s(z)

// ======== global structures ==========

struct INPUTS {

  // misc flags
  int fitsflag ;
  bool blind;      // blind cosmology results by adding sin(big number) 
  bool blind_auto; // automatically blind data and unblind sim
  int  blind_seed; // used to pick large random number for sin arg
  int  debug_flag ;

  GENPOLY_DEF zpoly_muerr_ideal; // define muerr = polynom(z)
  char string_muerr_ideal[100];

  int   speed_flag_chi2; // default = 1; set to 0 to disable
  bool  USE_SPEED_OFFDIAG; // internal: skip off-diag calc if chi2(diag)>threshold
  bool  USE_SPEED_INTERP;  // internal: intero r(z) and mu(z)

  int fitnumber;   // default=1; legacy for iterative fit after sigint calc

  char **HD_infile_list ;
  bool USE_HDIBC;
  int  NHD; // number of HDs; 1 or 2
  char outFile_cospar[MXCHAR_FILENAME] ; // output name of cospar file
  char outFile_resid[MXCHAR_FILENAME] ;
  char outFile_chi2grid[MXCHAR_FILENAME];
  char outFile_mucovtot_inv[MXCHAR_FILENAME];
  
  char **mucov_file ;  // input cov matrix(es); e.g., produced by create_cov
  int    NMUCOV ; // 1 or 2 cov matrices; 2 for HDIBC method
  
  char label_cospar[40]  ;   // string label for cospar file.
  int  ndump_mucov ; // dump this many column/rows
  char varname_muerr[40] ; // name of muerr variable, default MUERR

  // grid ranges and steps
  int    omm_steps, w0_steps, wa_steps ;
  double omm_min,   w0_min,   wa_min   ;    
  double omm_max,   w0_max,   wa_max   ;

  int dofit_wcdm ;    // default is true
  int dofit_lcdm ;    // option with -lcdm
  int dofit_w0wa ;    // option with -wa

  int use_bao, use_cmb ; // flags for bao and cmb priors
  char bao_sample[20];   // e.g., DESI_2024 or SDSS4_2020
  double rd_prior, rd_prior_sig ;
  
  int use_marg;          // flag to marginalize (instead of only minimize)
  int use_mucov;         // flag to read cov matrix
  double zmin, zmax;               // redshift selection

  int csv_out ; // optional csv format for output cospar and resids

  // Gauss OM prior
  double omm_prior, omm_prior_sig; 

  double snrms, sqsnrms, muerr_force ;
  double OMEGA_MATTER_SIM ; // for priors
  double w0_SIM, wa_SIM ;   // idem

  // computed from user inputs
  double omm_stepsize, w0_stepsize, wa_stepsize;
  int    format_cospar;

  char sig_type[40];
  bool use_sig_std, use_sig_68;
  
} INPUTS ;


typedef struct COVMAT {
  char    fileName[MXCHAR_FILENAME]; // file name that was read
  double *ARRAY1D ;   // 1D representation of matrix
  int     N_NONZERO ; // number of non-zero elements
  int     NDIM ;      // dimension size  
} COVMAT_DEF ;

// define workspace
struct  {

  // define variables interpolate rz for large samples
  int     n_exec_interp; // number of interpolate calls for r(z) and mu_cos(z)
  int     n_logz_interp; // number of logz bins for interpolation
  double *logz_list_interp, *z_list_interp, logz_bin_interp; 
  double *rz_list_interp, *mucos_list_interp;
  double  rz_dif_max ;

  // - - - -
  double *omm_val,  *w0_val,  *wa_val;
  double *omm_prob, *w0_prob, *wa_prob;
  double *omm_cdf,  *w0_cdf,  *wa_cdf;
  double *omm_sort, *w0_sort, *wa_sort;

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

  double omm_sig_std, omm_sig_upper, omm_sig_lower, omm_sig_final;
  double w0_sig_std,  w0_sig_upper,  w0_sig_lower,  w0_sig_final ;
  double wa_sig_std,  wa_sig_upper,  wa_sig_lower,  wa_sig_final ;
  double cov_w0wa, cov_w0omm;
  double rho_w0wa, rho_w0omm; 

  COVMAT_DEF MUCOV[2]; // up to two cov matrices
  COVMAT_DEF MUCOV_FINAL ;
  
  double w0_ran,   wa_ran,   omm_ran;
  double w0_final, wa_final, omm_final, chi2_final ;
  double sigmu_int, muoff_final, FoM_final ;

  int NWARN;

} WORKSPACE ;


typedef struct  {
  // Cosmological parameter structure used by codist.c 
  double omm;  // Omega_matter 
  double ome;  // Omega_energy
  double w0;   // equation of state of Dark Energy: 
  double wa;   //    w(z) = w0 + wa(1-a)  
  double mushift;  // for HDIBC method only
} Cosparam ;


// define structure to hold hubble diagram (Oct 1 2021)
typedef struct {
  int    NSN, NSN_ORIG ; // NSN after cuts, before cuts
  char   **cid ;
  bool   *pass_cut;
  double *mu, *mu_sig, *mu_ref, *mu_sim, *mu_sqsig, *z, *logz, *z_sig ;
  int    *nfit_perbin;  
  double zmin, zmax;
  bool   ISDATA_REAL;
  double muresid_avg;
  
  // for HDIBC method (Mar 2023)
  Cosparam cospar_biasCor;
  double *mu_cospar_biascor; // store distance(z_data) using biasCor cosmology
  double *f_interp;
  double *z_orig ;
  
} HD_DEF ;

HD_DEF HD_LIST[2];
HD_DEF HD_FINAL ;  // differs for HDIBC

Cosparam COSPAR_SIM ; // cosmo params from simulated data 
Cosparam COSPAR_LCDM; // w0,wa = -1,0

double *temp0_list ;
double *temp1_list ;
double *temp2_list ;

struct {
  double R, sigR;    // CMB R shift parameter and error
  double z, a ;      // redshift and 1/(1+z)
  char comment[200];

  int ranseed_R;     // random seed used to fluctuate R (1->internally compute)
} CMB_PRIOR;


#define MXDATA_BAO_PRIOR 40
#define IVAR_BAO_DM_over_rd   1
#define IVAR_BAO_DH_over_rd   2
#define IVAR_BAO_DV_over_rd   3

struct {

  int    NDATA;
  double REDSHIFT[MXDATA_BAO_PRIOR];
  double MEAN[MXDATA_BAO_PRIOR];
  double *COV1D, *COVINV1D ;            // store cov matrix in 1D flattened array
  char   VARNAME[MXDATA_BAO_PRIOR][20]; // e.g., DM_over_rd
  int    IVAR[MXDATA_BAO_PRIOR];        // e.g., = IVAR_BAO_DH_over_rd
  char   comment[200];
  char   comment_rd_prior[200];
} BAO_PRIOR;

// =========== global variables ================


double H0SIG     = 1000.0 ;       // error used in prior
double c_light   = 299792.458 ;   // speed of light, km/s 
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

// =========== function prototypes ================

void print_wfit_help(void);
void init_stuff(void);
void parse_args(int argc, char **argv) ;
int  compare_double_reverse (const void *, const void *);
void writemodel(char*, float, float, float);
void printerror( int status);
void read_HD(int index_HD, char *inFile, HD_DEF *HD);
bool read_ISDATA_REAL(char *inFile);
void compute_mu_biascor(char *inFile, HD_DEF *HD) ;
void read_cospar_biascor(char *info_yml_file, Cosparam *cospar);
void malloc_HDarrays(int opt, int NSN, HD_DEF *HD);
void malloc_COVMAT(int opt, COVMAT_DEF *COV);
void malloc_workspace(int opt);
void parse_VARLIST(FILE *fp);
void read_mucov(char *inFile, int imat, COVMAT_DEF *COV);
int  applyCut_HD(bool *PASS_CUT_LIST, HD_DEF *HD);
int  applyCut_COVMAT(bool *PASS_CUT_LIST, COVMAT_DEF *MUCOV);
void dump_MUCOV(COVMAT_DEF *COV, char *comment);
void sync_HD_LIST(HD_DEF *HD_REF, HD_DEF *HD, COVMAT_DEF *MUCOV);
void sync_HD_redshifts(HD_DEF *HD0, HD_DEF *HD1) ;
void compute_MUCOV_FINAL();
void invert_mucovar(COVMAT_DEF *COV, double sqmurms_add);
void check_invertMatrix(int N, double *COV, double *COVINV );
void set_stepsizes(void);
void set_Ndof(void);
void init_rz_interp(HD_DEF *HD);
void exec_rz_interp(int k, Cosparam *cospar, double *rz, double *dmu);
void check_refit(void);

void wfit_minimize(void);
void prep_speed_offdiag(double extchi_tmp);
void wfit_normalize(void);
void wfit_marginalize(void);
void wfit_uncertainty(void);
void wfit_uncertainty_fitpar(char *varname);
void wfit_final(void);
void wfit_FoM(void);
void wfit_Covariance(void);


void get_chi2_fit(double w0, double wa, double OM, double sqmurms_add,
		  double *z_list, double *mu_list, double *f_interp_list,
		  double *muresid_avg, double *chi2sn, double *chi2tot );
void getname(char *basename, char *tempname, int nrun);

double get_DMU_chi2wOM(double z, double rz, double mu); 
double get_mu_cos(double z, double rz) ;

double get_minimized_cospar( double *w0_atchimin, double *wa_atchimin, 
			       double *OM_atchimin ); 

void   set_priors(void);
void   init_bao_prior(int OPT) ;
void   init_bao_cov(void);
double rd_bao_prior(int OPT, Cosparam *cpar) ;
double c_sound(int OPT, double z, double H0);
double DM_bao_prior(double z, Cosparam *cpar);
double DH_bao_prior(double z, Cosparam *cpar);
double DV_bao_prior(double z, Cosparam *cpar);
void   dump_bao_prior(void);

void   init_cmb_prior(int OPT) ;
double get_Rcmb_ranshift(void) ;
double chi2_bao_prior(Cosparam *cpar, double *rd);
double chi2_cmb_prior(Cosparam *cpar);
void   get_chi2_priors(Cosparam *cpar, double *chi2_om, double *chi2_cmb,
		       double *chi2_bao, double *chi2_rd);

int    ISEED_UNIQUE(void);

void set_HzFUN_for_wfit(double H0, double OM, double OE, double w0, double wa,
                        HzFUN_INFO_DEF *HzFUN ) ;

int  cidindex(char *cid);

void WRITE_OUTPUT_DRIVER(void);
void write_output_resid(void);
void write_output_cospar(void);
void write_output_chi2grid(void);

void CPU_summary(void);

// cosmology functions
double EofZ(double z, Cosparam *cptr);
double one_over_EofZ(double z, Cosparam *cptr);
double codist(double z, Cosparam *cptr); // maybe change name to Ezinv_integral
void   test_codist(void);
void   test_c_sound(void);

// Nov 1 2021 R.Kessler - new set of functions to integrate over "a"
double EofA(double a, Cosparam *cptr);
double one_over_EofA(double a, Cosparam *cptr);
double cs_over_EofA(double a, Cosparam *cptr);
double Eainv_integral(double amin, double amax, Cosparam *cptr) ;
double cs_over_E_integral(double amin, double amax, Cosparam *cptr);
  
// from simpint.h
#define EPS_CONVERGE_POSomm  1.0e-6
#define EPS_CONVERGE_NEGomm  1.0e-3
#define JMAX 20

double simpint(double (*func)(double, Cosparam *), 
	       double x1, double x2, Cosparam *vptr);

double trapint(double (*func)(double, Cosparam *), 
	       double x1, double x2, Cosparam *vptr);

double trapezoid(double (*func)(double, Cosparam *), 
		 double x1, double x2, int n, double s, Cosparam *vptr);

void test_cospar(void);


// =============================================
// ================ MAIN =======================
// =============================================

int main(int argc,char *argv[]){

  int f ;
  // ----------------- BEGIN MAIN ----------------

  print_full_command(stdout, argc, argv);

  t_start = time(NULL);

  set_EXIT_ERRCODE(EXIT_ERRCODE_wfit);

  // Give help if no arguments
  if (argc < 2) { print_wfit_help();  exit(0);  }

  // init variables
  init_stuff(); 

  //  test_codist();  only to test r(z) integral
  // test_c_sound();
  
  // Parse command line args 
  parse_args(argc,argv);

  if ( INPUTS.debug_flag == 91 ) { test_cospar(); }

  malloc_workspace(+1);

  printf("# =============================================== \n");
  printf(" SNANA_DIR = %s \n", getenv("SNANA_DIR") );

  while (INPUTS.fitnumber <= 1) {

    /************************/
    /*** Read in the data ***/
    /************************/
    
    read_HD(0, INPUTS.HD_infile_list[0], &HD_LIST[0]); 
    if ( INPUTS.USE_HDIBC ) 
      { read_HD(1, INPUTS.HD_infile_list[1], &HD_LIST[1]); }

    // for large samples, setup logz grid to interpolate rz(z)
    init_rz_interp(&HD_LIST[0]);
    
    // Set BAO and CMB priors
    set_priors();
  
    // read optional mu-covSys or mucovtot_inv matrix
    for(f=0; f < INPUTS.NMUCOV; f++ ) 
      { read_mucov(INPUTS.mucov_file[f], f, &WORKSPACE.MUCOV[f] ); }

    // sync HD events and MUCOV if there are two HDs
    if ( INPUTS.USE_HDIBC ) { 
      sync_HD_LIST(&HD_LIST[0],
		   &HD_LIST[1], &WORKSPACE.MUCOV[1] );  // <== modified

      sync_HD_LIST(&HD_LIST[1],
		   &HD_LIST[0], &WORKSPACE.MUCOV[0] );  // <== modified

      sync_HD_redshifts(&HD_LIST[0], &HD_LIST[1]); 
    }

    if ( INPUTS.use_mucov ) {
      compute_MUCOV_FINAL();
      invert_mucovar(&WORKSPACE.MUCOV_FINAL, INPUTS.sqsnrms);
    }
    
    // compute grid step size per floated variable
    set_stepsizes();

    // compute number of degrees of freedom
    set_Ndof(); 

    printf("\n# ======================================= \n");
    print_cputime(t_start, STRING_CPUTIME_INIT, UNIT_TIME_SECOND, 0);
    t_end_init = time(NULL);
 
    // minimize chi2 on a grid
    wfit_minimize(); 
    
    // Normalize probability distributions 
    wfit_normalize();

    wfit_marginalize();  // marginalize

    // get uncertainties
    wfit_uncertainty();
      
    // determine "final" quantities, including sigma_mu^int
    wfit_final();

    // Compute covriance with fitted parameters
    wfit_Covariance();
    
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
  malloc_workspace(-1);
  CPU_summary();

  printf("DONE. \n"); fflush(stdout);

  return(0);
}   // end main

// ================================
void init_stuff(void) {
  
  char fnam[] = "init_stuff" ;

  // ------------ BEGIN -----------

  INPUTS.USE_HDIBC = false;
  INPUTS.NHD       = 0 ;

  INPUTS.blind = INPUTS.blind_auto = INPUTS.fitsflag = INPUTS.debug_flag = 0;
  INPUTS.blind_seed = 48901 ;

  init_GENPOLY(&INPUTS.zpoly_muerr_ideal);
  INPUTS.string_muerr_ideal[0] = 0 ;

  INPUTS.speed_flag_chi2 = SPEED_FLAG_CHI2_DEFAULT ;

  INPUTS.OMEGA_MATTER_SIM = OMEGA_MATTER_DEFAULT ;
  INPUTS.w0_SIM           = w0_DEFAULT ;
  INPUTS.wa_SIM           = wa_DEFAULT ;

  INPUTS.NMUCOV              = 0 ;
  INPUTS.ndump_mucov         = 0 ;
  INPUTS.outFile_cospar[0]   = 0 ;
  INPUTS.outFile_resid[0]    = 0 ;
  INPUTS.outFile_chi2grid[0] = 0 ;
  INPUTS.outFile_mucovtot_inv[0] = 0 ;
  INPUTS.use_mucov           = 0 ;
  sprintf(INPUTS.label_cospar,"none");
  INPUTS.format_cospar = 1; // csv like format
  INPUTS.fitnumber = 1;
  INPUTS.varname_muerr[0] = 0 ;

  INPUTS.use_sig_std  = false;   // legacy option
  INPUTS.use_sig_68   = true;  
   
  // Gauss OM params
  INPUTS.omm_prior        = OMEGA_MATTER_DEFAULT ;  // mean
  INPUTS.omm_prior_sig    = 99.0;  // sigma

  INPUTS.rd_prior         = 147.0 ;
  INPUTS.rd_prior_sig     = 1.0E8 ;

  
  init_bao_prior(-9); 
  init_cmb_prior(-9);  

  INPUTS.dofit_wcdm = 1 ;
  INPUTS.dofit_lcdm = 0 ;
  INPUTS.dofit_w0wa = 0 ;

  INPUTS.use_bao  = INPUTS.use_cmb = 0 ;
  INPUTS.bao_sample[0] = 0 ;
  INPUTS.use_marg = 1;

  INPUTS.omm_steps =  DEFAULT_omm_steps;
  INPUTS.omm_min   =  DEFAULT_omm_min ;
  INPUTS.omm_max   =  DEFAULT_omm_max;
  INPUTS.w0_steps  =  DEFAULT_w0_steps;
  INPUTS.w0_min    =  DEFAULT_w0_min ;
  INPUTS.w0_max    =  DEFAULT_w0_max ;  
  INPUTS.wa_steps  =  DEFAULT_wa_steps ;
  INPUTS.wa_min    =  DEFAULT_wa_min ;
  INPUTS.wa_max    =  DEFAULT_wa_max ;

  INPUTS.zmin = 0.0;  INPUTS.zmax = 9.0;

  INPUTS.snrms  = INPUTS.sqsnrms = 0.0 ;
  INPUTS.muerr_force = -9.0 ;

  // - - - - - - - - - - -  - -
  // init misc variables
  SIG_MUOFF   = 5.0 * log10(1. + H0SIG/H0_SALT2);
  SQSIG_MUOFF = SIG_MUOFF * SIG_MUOFF ;

  sprintf(varname_w,   "w"  );  // default wCDM model
  sprintf(varname_wa,  "wa" );
  sprintf(varname_omm, "OM" );

  // - - - - -
  // WORKSPACE

  WORKSPACE.MUCOV[0].NDIM = 0;
  WORKSPACE.MUCOV[1].NDIM = 0;

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
void print_wfit_help(void) {

  // Created Oct 1 2021
  // [moved from main]

  static char *help_main[] = {		/* instructions */
    "",
    "WFIT - Usage: wfit [hubble diagram] (options)",
    "WFIT - Usage: wfit [HD0,HD1] (options)   for HDIBC method",
    "",
    " Default for input file is 3 columns: z, mu, sigma_mu, but see below",
    "",
    "  Options:",
    "   -wa\t\tfit for w0 & wa [w=w0+wa*(1-a)] from ",
    "     \t\t Chevallier & Polarski, 2001 [Int.J.Mod.Phys.D10,213(2001)]",
    "   -lcdm\tfit for OM only (fix w=-1)" ,
    "   -ompri\tcentral value of omega_m prior [default: Planck2018]", 
    "   -dompri\tGauss sigma of omega_m prior [default: 99]",
    //
    "   -bao DESI_2024\tBAO prior (mean+cov) from DESI 2024",
    "   -bao_sim DESI_2024\tBAO cov from DESI 2024, mean from sim cosmology params",
    "   -bao SDSS4_2020\tBAO prior (mean+cov) from SDSS4 2020",
    "   -bao_sim SDSS4_2020\tBAO cov from SDSS4 2020, mean from sim cosmology params",
    "   -rd_prior \t rd prior (with BAO prior only)  [default=147]",
    "   -rd_prior_sig \t Gauss sigma for rd prior (with BAO prior only) [default=999.]",
    //
    "   -cmb\t\tuse CMB prior from 5-year WMAP (2008)",
    "   -cmb_sim\tuse CMB prior with simulated cosmology and WMAP formalism",
    "   -om_sim\tOmega_M for cmb_sim (default from sntools.h)",
    "   -w0_sim\tw0 for cmb_sim (default from sntools.h)",
    "   -wa_sim\twa for cmb_sim (default from sntools.h)",
    "   -minchi2\tget w and OM from minchi2 instead of marginalizing",
    "   -marg\tget w and OM from marginalizing (default)",
    "   -sig_type 68 \tUncertainty from 68 percent confidence (default)",
    "   -sig_type std\tUncertainty from STD of marginalized pdf",
    "   -Rcmb\tCMB comstraints: R = Rcmb +/- sigma_Rcmb [= 1.710 +/- 0.019]",
    "   -sigma_Rcmb\tUncertainty on Rcmb",
    "   -ranseed_Rcmb random seed to fluctuate Rcmb (1-> internally compute seed)",
    "   -snrms\tadd this to reported errors (mags) [default: 0]",
    "   -muerr_force\tforce this mu_sig on all events",
    "   -zmin\tFit only data with z>zmin",
    "   -zmax\tFit only data with z<zmax",    
    "   -blind\tIf set, results are obscured with sin(large random) ",
    "   -blind_auto\tBlind data, unblind sim; requires ISDATA_REAL in HD file",
    "   -blind_seed\tSeed to pick large random numbers for sin arg ",
    "   -mucovsys_file\tfile with COV_syst e.g., from create_covariance",
    "   -mucov_file   \tlegacy key for mucovsys_file",
    "   -mucovtot_inv_file\tfile with inverse of COVTOT",
    "   -ndump_mucov\t dump this many rows/columns of MUCOV and MUCOVINV",
    "   -varname_muerr\t column name with distance errors (default=MUERR)",
    "   -refit\tfit once for sigint then refit with snrms=sigint.", 
    "   -speed_flag_chi2   +=1->offdiag trick, +=2->interp trick",
    "   -debug_flag 91\t compare calc mu(wfit) vs. mu(sim)",
    "   -muerr_ideal  replace all mu with mu_true + Gauss(0,muerr);",
    "                 e.g.,  muerr_ideal 0.1,0.01,0.05 -> "
    "muerr = .1 + .01*z + .05*z^2",
    "                 e.g.,  muerr_ideal muerr -> use original muerr from data",
    "",
    0
  };

  char help_grid[1000];
  sprintf(help_grid,
	  " Grid spacing:\n"
	  " wCDM Fit: \n"
	  "   -ommin/-ommax/-omsteps   OM default grid [%4.1f / %4.1f / %3d] \n"
	  "   -wmin /-wmax /-wsteps    w  default grid [%4.1f / %4.1f / %3d] \n"
	  "\n"
	  " w0-wa Fit:\n"
	  "   -ommin/-ommax/-omsteps   OM default grid [%4.1f / %4.1f / %3d] \n"
	  "   -w0min/-w0max/-w0steps   w0 default grid [%4.1f / %4.1f / %3d] \n"
	  "   -wamin/-wamax/-wasteps   wa default grid [%4.1f / %4.1f / %3d] \n" 
	  "\n",
	  DEFAULT_omm_min, DEFAULT_omm_max, DEFAULT_omm_steps,
	  DEFAULT_w0_min,  DEFAULT_w0_max,  DEFAULT_w0_steps,
	  //
	  DEFAULT_omm_min, DEFAULT_omm_max, DEFAULT_omm_steps,
	  DEFAULT_w0_min,  DEFAULT_w0_max,  DEFAULT_w0_steps,
	  DEFAULT_wa_min,  DEFAULT_wa_max,  DEFAULT_wa_steps
	  );

  static char *help_output[] = {
    " Output:",
    "   -outfile_cospar\tname of output file with fit cosmo params",
    "   -cospar\t\t  [same as previous]",
    "   -cospar_yaml\tname of output YAML file with fit cosmo params",
    "   -outfile_resid\tname of output file with mu-residuals",
    "   -resid\t\t  [same as previous]" ,
    "   -outfile_chi2grid\tname of output file containing chi2-grid",
    "   -label\tstring-label for cospar file.",
    "   -outfile_mucovtot_inv\t Write mucovtot_inv to file",
    "",
    0
  };

  int i; 
  
  // ----------- BEGIN ------------
  
  for (i = 0; help_main[i] != 0; i++)
    { printf ("%s\n", help_main[i]); }


  printf("%s", help_grid);

  for (i = 0; help_output[i] != 0; i++)
    { printf ("%s\n", help_output[i]); }    

  return;

} // end print_wfit_help


// ==================================
void parse_args(int argc, char **argv) {

  // Created Oct 1 2021 [code moved from main]
  // 
  int N_HDfile ;
  int iarg;
  char fnam[] = "parse_args" ;

  // ------------ BEGIN ------------

  // parse HD as 1st positional arg
  parse_commaSepList("HD_infile", argv[1], 2, 2*MXCHAR_FILENAME,
		     &INPUTS.NHD, &INPUTS.HD_infile_list );

  if ( INPUTS.NHD == 2 ) { INPUTS.USE_HDIBC = true; }

  for (iarg=2; iarg<argc; iarg++) {
    if (argv[iarg][0]=='-') {

      if (strcasecmp(argv[iarg]+1,"ompri")==0) { 
	INPUTS.omm_prior = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"dompri")==0) { 
	INPUTS.omm_prior_sig = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"bao")==0) { 
	sprintf(INPUTS.bao_sample,"%s", argv[++iarg]);
	INPUTS.use_bao=1;
      }
      else if (strcasecmp(argv[iarg]+1,"bao_sim")==0) {
	sprintf(INPUTS.bao_sample,"%s", argv[++iarg]);
        INPUTS.use_bao=2;
      }
      else if (strcasecmp(argv[iarg]+1,"rd_prior")==0) {
	INPUTS.rd_prior     = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"rd_prior_sig")==0) {
	INPUTS.rd_prior_sig = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"csv")==0) { 
	INPUTS.csv_out = 1;
      }
      else if (strcasecmp(argv[iarg]+1, "Rcmb")==0) {
	CMB_PRIOR.R = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1, "sigma_Rcmb")==0) {
	CMB_PRIOR.sigR = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1, "ranseed_Rcmb")==0) {
	CMB_PRIOR.ranseed_R = atoi(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"cmb")==0) { 
	INPUTS.use_cmb = 1;
	INPUTS.omm_prior_sig = 0.9;  // turn off omm prior (9/2021)
      }
      else if (strcasecmp(argv[iarg]+1,"wa")==0) {
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

      }
      else if (strcasecmp(argv[iarg]+1,"om_sim")==0) {
	INPUTS.OMEGA_MATTER_SIM = atof(argv[++iarg]); 
	INPUTS.omm_prior        = INPUTS.OMEGA_MATTER_SIM ;
      }
      else if (strcasecmp(argv[iarg]+1,"w0_sim")==0) {
	INPUTS.w0_SIM = atof(argv[++iarg]); 
      }
      else if (strcasecmp(argv[iarg]+1,"wa_sim")==0) {
	INPUTS.wa_SIM = atof(argv[++iarg]); 
      }
      else if (strcasecmp(argv[iarg]+1,"minchi2")==0) { 
	INPUTS.use_marg = 0 ;
      }
      else if (strcasecmp(argv[iarg]+1,"marg")==0) { 
	INPUTS.use_marg=1;
      }
      else if (strcasecmp(argv[iarg]+1,"sig_type")==0) {
	sprintf(INPUTS.sig_type,"%s", argv[++iarg]);
	if ( strcmp(INPUTS.sig_type,"std")== 0 )
	  { INPUTS.use_sig_std = true; INPUTS.use_sig_68= false; }
	else if ( strcmp(INPUTS.sig_type,"68")== 0 )
	  { INPUTS.use_sig_std = false; INPUTS.use_sig_68= true; }
	else {
	  sprintf(c1err,"Invalid sig_type = '%s'", INPUTS.sig_type);
	  sprintf(c2err,"Valid sig_types are std and 68");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
	}

      }
      else if (strcasecmp(argv[iarg]+1,"snrms")==0) { 
	INPUTS.snrms = atof(argv[++iarg]);
	INPUTS.sqsnrms = INPUTS.snrms * INPUTS.snrms ;
      }
      else if (strcasecmp(argv[iarg]+1,"muerr_force")==0) { 
	INPUTS.muerr_force = atof(argv[++iarg]); // March 2023
      }
      else if (strcasecmp(argv[iarg]+1,"refit")==0){
        INPUTS.fitnumber -= 1;  
      }
      else if (strcasecmp(argv[iarg]+1,"zmin")==0) { 
	INPUTS.zmin = atof(argv[++iarg]); 
      }
      else if (strcasecmp(argv[iarg]+1,"zmax")==0) { 
	INPUTS.zmax = atof(argv[++iarg]); 	

      }
      else if (strcasecmp(argv[iarg]+1,"blind")==0) { 
	INPUTS.blind = true ;

      } else if (strcasecmp(argv[iarg]+1,"blind_auto")==0) { 
	INPUTS.blind_auto = true ;

      } else if (strcasecmp(argv[iarg]+1,"blind_seed")==0) { 
	INPUTS.blind_seed = atoi(argv[++iarg]); 	
	
      } else if (strcasecmp(argv[iarg]+1,"fitswrite")==0) {   
	/* Output likelihoods as fits files & text files */
	INPUTS.fitsflag = 1; 

      }
      else if ( strcasecmp(argv[iarg]+1,"mucovsys_file")==0 ||
		strcasecmp(argv[iarg]+1,"mucov_file"   )==0    ) {

	parse_commaSepList("mucov_file", argv[++iarg], 2, 2*MXCHAR_FILENAME,
			   &INPUTS.NMUCOV, &INPUTS.mucov_file );
	INPUTS.use_mucov =1 ;  // flag that mucovsys has been read

      }
      else if (strcasecmp(argv[iarg]+1,"mucovtot_inv_file")==0) {
	parse_commaSepList("mucov_file", argv[++iarg], 2, 2*MXCHAR_FILENAME,
			   &INPUTS.NMUCOV, &INPUTS.mucov_file );
	INPUTS.use_mucov = 2 ;  // flag that mucovtot_inv has been read
	
      }
      else if (strcasecmp(argv[iarg]+1,"varname_muerr")==0) {
        strcpy(INPUTS.varname_muerr,argv[++iarg]);
       
      }
      else if (strcasecmp(argv[iarg]+1,"ndump_mucov")==0 ) { 
	INPUTS.ndump_mucov = atoi(argv[++iarg]);      

      /* Change Omega_m grid parameters? */
      }
      else if (strcasecmp(argv[iarg]+1,"ommin")==0) {   
	INPUTS.omm_min = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"ommax")==0) {   
	INPUTS.omm_max = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"omsteps")==0) {   
	INPUTS.omm_steps = (int)atof(argv[++iarg]);
      }

      /* accept w or w0 */
      else if (strcasecmp(argv[iarg]+1,"w0min")==0) {   
	INPUTS.w0_min = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"w0max")==0) {   
	INPUTS.w0_max = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"w0steps")==0) {   
	INPUTS.w0_steps = (int)atof(argv[++iarg]);
      }

      else if (strcasecmp(argv[iarg]+1,"wmin")==0) {   
	INPUTS.w0_min = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"wmax")==0) {   
	INPUTS.w0_max = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"wsteps")==0) {   
	INPUTS.w0_steps = (int)atof(argv[++iarg]);
      }

      /* Change wa grid parameters? */
      else if (strcasecmp(argv[iarg]+1,"wamin")==0) {
        INPUTS.wa_min = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"wamax")==0) {
        INPUTS.wa_max = atof(argv[++iarg]);
      }
      else if (strcasecmp(argv[iarg]+1,"wasteps")==0) {
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
      // xxx else if (strcasecmp(argv[iarg]+1,"chi2grid")==0)  
      // xxx { strcpy(INPUTS.outFile_chi2grid,argv[++iarg]); }

      else if (strcasecmp(argv[iarg]+1,"outfile_mucovtot_inv")==0)  
	{ strcpy(INPUTS.outFile_mucovtot_inv,argv[++iarg]); }
      
      else if (strcasecmp(argv[iarg]+1,"debug_flag")==0)  
	{ INPUTS.debug_flag = atoi(argv[++iarg]);  } 

      else if (strcasecmp(argv[iarg]+1,"muerr_ideal")==0)  { 
	char *string_muerr_ideal = INPUTS.string_muerr_ideal;
	strcpy(string_muerr_ideal,argv[++iarg]); 
	if ( strcmp(string_muerr_ideal,"muerr") == 0 ) {   
	  INPUTS.zpoly_muerr_ideal.ORDER = 100; // flag to use data muerr
	}
	else {
	  parse_GENPOLY(string_muerr_ideal, "zpoly_muerr_ideal",
			&INPUTS.zpoly_muerr_ideal, fnam ) ;
	}  
      }

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

  COSPAR_SIM.omm = INPUTS.OMEGA_MATTER_SIM ;
  COSPAR_SIM.ome = 1.0 - INPUTS.OMEGA_MATTER_SIM ;
  COSPAR_SIM.w0  = INPUTS.w0_SIM;
  COSPAR_SIM.wa  = INPUTS.wa_SIM;
  COSPAR_SIM.mushift  = 0.0 ;

  COSPAR_LCDM.omm =  INPUTS.OMEGA_MATTER_SIM ;
  COSPAR_LCDM.ome =  1.0 - INPUTS.OMEGA_MATTER_SIM ;
  COSPAR_LCDM.w0  = -1.0 ;
  COSPAR_LCDM.wa  =  0.0 ;
  COSPAR_LCDM.mushift  = 0.0 ;

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
    WORKSPACE.omm_val  = (double *)calloc(INPUTS.omm_steps, memd);
    WORKSPACE.w0_val   = (double *)calloc(INPUTS.w0_steps,  memd);
    WORKSPACE.wa_val   = (double *)calloc(INPUTS.wa_steps,  memd);

    WORKSPACE.omm_prob  = (double *)calloc(INPUTS.omm_steps, memd);
    WORKSPACE.w0_prob   = (double *)calloc(INPUTS.w0_steps,  memd);
    WORKSPACE.wa_prob   = (double *)calloc(INPUTS.wa_steps,  memd);
    
    WORKSPACE.omm_sort  = (double *)calloc(INPUTS.omm_steps, memd);
    WORKSPACE.w0_sort   = (double *)calloc(INPUTS.w0_steps,  memd);
    WORKSPACE.wa_sort   = (double *)calloc(INPUTS.wa_steps,  memd);

    WORKSPACE.omm_cdf  = (double *)calloc(INPUTS.omm_steps, memd);
    WORKSPACE.w0_cdf   = (double *)calloc(INPUTS.w0_steps,  memd);
    WORKSPACE.wa_cdf   = (double *)calloc(INPUTS.wa_steps,  memd);    

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
    free(WORKSPACE.omm_val);
    free(WORKSPACE.w0_val);
    free(WORKSPACE.wa_val);

    free(WORKSPACE.omm_prob);
    free(WORKSPACE.w0_prob);
    free(WORKSPACE.wa_prob);

    free(WORKSPACE.omm_cdf);
    free(WORKSPACE.w0_cdf);
    free(WORKSPACE.wa_cdf);

    free(WORKSPACE.omm_sort);
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
void  malloc_HDarrays(int opt, int NSN, HD_DEF *HD) {

  // Created Oct 01 2021 by R.Kessler
  // malloc arrays to store Hubble diagram data
  // opt > 0 -> malloc
  // opt < 0 -> free

  int i;
  char fnam[] = "malloc_HDarrays" ;
  
  // --------- BEGIN --------

  if ( opt > 0 ) {
    HD->cid = (char**) malloc( NSN * sizeof(char*) );
    for(i=0; i < NSN; i++ ) { HD->cid[i] = (char*)malloc( 20*sizeof(char) ); }

    HD->pass_cut    = (bool   *)calloc(NSN,sizeof(bool));
    HD->mu          = (double *)calloc(NSN,sizeof(double));
    HD->mu_sig      = (double *)calloc(NSN,sizeof(double));
    HD->mu_ref      = (double *)calloc(NSN,sizeof(double));
    HD->mu_sim      = (double *)calloc(NSN,sizeof(double));    
    HD->mu_sqsig    = (double *)calloc(NSN,sizeof(double));   
    HD->z           = (double *)calloc(NSN,sizeof(double));
    HD->logz        = (double *)calloc(NSN,sizeof(double));
    HD->z_sig       = (double *)calloc(NSN,sizeof(double));
    HD->f_interp    = (double *)calloc(NSN,sizeof(double));
    HD->z_orig      = (double *)calloc(NSN,sizeof(double));        
    
    HD->nfit_perbin = (int    *)calloc(NSN,sizeof(int));
    if ( INPUTS.USE_HDIBC ) 
      { HD->mu_cospar_biascor = (double*)calloc(NSN,sizeof(double)); }
  }
  else {
    free(HD->pass_cut);
    free(HD->mu); 
    free(HD->mu_sig); 
    free(HD->mu_ref);
    free(HD->mu_sim);     
    free(HD->mu_sqsig); 
    free(HD->z); 
    free(HD->logz); 
    free(HD->z_sig);
    free(HD->f_interp);
    free(HD->z_orig);    
    free(HD->nfit_perbin); 
    for(i=0; i < NSN; i++ ) { free(HD->cid[i]); } 
    free(HD->cid);

    if ( INPUTS.USE_HDIBC )  { free(HD->mu_cospar_biascor); }
  }

} // end malloc_HDarrays

// ==========================
void malloc_COVMAT(int opt, COVMAT_DEF *COVMAT) {

  int NDIM = COVMAT->NDIM ;
  long long int MEMD = NDIM * NDIM * sizeof(double);
  // ------------ BEGIN ------------

  if ( opt > 0 )  {
    COVMAT->ARRAY1D = (double*) malloc(MEMD);
  }
  else {
    free(COVMAT->ARRAY1D);
  }
  
  return ;
} // end malloc_COVMAT

// ==================================
void read_HD(int index_HD, char *inFile, HD_DEF *HD) {

  // Created Oct 1 2021 by R.Kessler
  // Refactored routine to read Hubble diagram from fitres-formatted file
  // using SNANA read utilities.
  // Mar 2023: refactor to pass HD struct to enable reading 1 or 2 HDs

  int IVAR_ROW=-8, IVAR_MU=-8, IVAR_MUERR=-8, IVAR_MUREF;
  int IVAR_zHD=-8, IVAR_zHDERR=-8, IVAR_NFIT=-8 ;
  int IFILETYPE, NVAR_ORIG, LEN, NROW, irow ;    
  int VBOSE = 1;
  double rz, mu_cos, ztmp, sigtmp ;
  bool ISDATA_REAL = false ;
  char TBNAME[] = "HD" ;  // table name is Hubble diagram
  char fnam[] = "read_HD" ;

  // --------------- BEGIN --------------

  if ( INPUTS.blind_auto ) { ISDATA_REAL = read_ISDATA_REAL(inFile); }
  HD->ISDATA_REAL = ISDATA_REAL ;

  TABLEFILE_INIT();
  NROW      = SNTABLE_NEVT(inFile,TBNAME);
  IFILETYPE = TABLEFILE_OPEN(inFile,"read");
  NVAR_ORIG = SNTABLE_READPREP(IFILETYPE,TBNAME);

  malloc_HDarrays(+1, NROW, HD); 

  IVAR_ROW   = SNTABLE_READPREP_VARDEF(VARLIST_DEFAULT_CID,
				       HD->cid, NROW, VBOSE);

  IVAR_MU    = SNTABLE_READPREP_VARDEF(VARLIST_DEFAULT_MU, 
				       HD->mu,  NROW, VBOSE) ;

  char STRING_MUERR[100] ;
  if ( strlen(INPUTS.varname_muerr) > 0 )  // command line override
    { sprintf(STRING_MUERR,"%s:D",INPUTS.varname_muerr); }
  else
    { sprintf(STRING_MUERR,"%s", VARLIST_DEFAULT_MUERR) ;  }

  IVAR_MUERR = SNTABLE_READPREP_VARDEF(STRING_MUERR,
                                       HD->mu_sig, NROW, VBOSE) ;
  
  IVAR_MUREF = SNTABLE_READPREP_VARDEF(VARLIST_DEFAULT_MUREF,
				       HD->mu_ref, NROW, VBOSE);
  IVAR_NFIT  = SNTABLE_READPREP_VARDEF(VARLIST_DEFAULT_NFIT,
				       HD->nfit_perbin, NROW, VBOSE );
  // - - - -
  IVAR_zHD    = SNTABLE_READPREP_VARDEF(VARLIST_DEFAULT_zHD, 
					HD->z,     NROW, VBOSE) ;
  IVAR_zHDERR = SNTABLE_READPREP_VARDEF(VARLIST_DEFAULT_zHDERR,
					HD->z_sig, NROW, VBOSE) ;

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
  HD->NSN_ORIG = SNTABLE_READ_EXEC();


  // for MUDIF output from BBC, MU is actually MUDIF,
  // so set MU += MUREF
  if ( IVAR_MUREF > 0 ) {
    printf("\t Input is MUDIF File: MU -> MU_REF + MUDIF \n");
    for(irow=0; irow < NROW; irow++ ) { HD->mu[irow] += HD->mu_ref[irow]; }
  }

  // - - - - - -
  // apply cuts and store mu_sim per event
  
  int    PASSCUTS,  NFIT, NROW2=0 ;  
  HD->zmin = 99999.0;  HD->zmax = -9999.90 ;

  for(irow=0; irow < NROW; irow++ ) {

    NFIT = 999;
    if( IVAR_NFIT > 0 ) { NFIT = HD->nfit_perbin[irow] ; }

    ztmp = HD->z[irow] ;
    PASSCUTS = 
      (ztmp >= INPUTS.zmin) &&
      (ztmp <= INPUTS.zmax) &&
      (NFIT > 1) ;

    HD->pass_cut[irow]  = PASSCUTS ;

    rz                = codist(ztmp, &COSPAR_SIM);
    HD->mu_sim[irow]  = get_mu_cos(ztmp, rz);
      
  } // end irow


  // - - - - - - - - - - - -
  // compute mu_biascor before applutCut_HD re-shuffles lists.
  if ( INPUTS.USE_HDIBC ) {  compute_mu_biascor(inFile, HD);  } 

  
  // apply cuts by removing rejected events from HD lists
  // and from cov matrix.
  HD->NSN = applyCut_HD(HD->pass_cut, HD);

  printf("   Keep %d of %d z-bins with z cut\n", 
	 HD->NSN, NROW);  

  // - - - - - 
  // Mar 29 2023: test with ideal mu and mu_err 
  double muerr_ideal, gran ;
  GENPOLY_DEF *zpoly_muerr_ideal = &INPUTS.zpoly_muerr_ideal;
  int         ORDER              = zpoly_muerr_ideal->ORDER ;
  int ISEED ;

  if ( ORDER >= 0.0 ) {

    // compute SEED using first MU-uncertainty so that
    // each data set has a unique randome seed.
    ISEED = ISEED_UNIQUE();
    init_random_seed(ISEED,1);
 
    printf("\n");
    printf("   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("\n\t TEST: Replace MU -> MU_SIM + N(0,muerr) \n");

    if ( ORDER < 10 )
	{ print_GENPOLY(zpoly_muerr_ideal); }
    else
      { printf("\t\t muerr = original muerr\n"); }
	  
    printf("\t ISEED = %d \n", ISEED );
    printf("\n   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("\n");

    for(irow=0; irow < HD->NSN; irow++ ) {
      ztmp               = HD->z[irow] ;
      sigtmp             = HD->mu_sig[irow];
      gran               = unix_getRan_Gauss(0);
      rz                 = codist(ztmp, &COSPAR_SIM);
      mu_cos             = get_mu_cos(ztmp, rz);

      if ( ORDER < 10 ) 
	{ muerr_ideal  = eval_GENPOLY(ztmp, zpoly_muerr_ideal, fnam); }
      else
	{ muerr_ideal = sigtmp; } // use original muerr

      HD->mu[irow]       = mu_cos + (muerr_ideal * gran) ;
      HD->mu_sig[irow]   = muerr_ideal;
      HD->mu_sqsig[irow] = muerr_ideal * muerr_ideal ;
    }
  } // end muerr_ideal

  fflush(stdout);
  //  debugexit(fnam);


  // - - - -
  if ( index_HD == 0 ) {
    int MEMD = (NROW+10) * sizeof(double) ;
    printf("\t malloc temp[0,1,2]_list (NROW=%d) for get_chi2_fit calls.\n",
	   NROW);
    temp0_list  = (double*)malloc(MEMD);
    temp1_list  = (double*)malloc(MEMD);
    temp2_list  = (double*)malloc(MEMD);
  }
  
  printf("\n"); fflush(stdout);

  return;

} // end read_HD


// ============================================================
void compute_mu_biascor(char *inFile, HD_DEF *HD) {

  // for HDIBC method, read cospar_biascor from INFO.YML in same
  // directory as HD; then compute mu_biascor for each event;
  // these distances are used later to interpolate.

  int NSN = HD->NSN_ORIG ; // orignal number before cuts
  char *e, path[MXCHAR_FILENAME], info_yml_file[MXCHAR_FILENAME];
  int  irow, jslash = -9 ;
  double ztmp, rz, mu_cos ;
  char fnam[] = "compute_mu_biascor" ;

  // ------------- BEGIN ------------

  if ( strchr(inFile,'/') != NULL ) {
      e      = strrchr(inFile, '/');
      jslash = (int)(e - inFile);
      sprintf(path,"%s", inFile); path[jslash] = 0;
      sprintf(info_yml_file,"%s/%s", path, INFO_YML_FILENAME );
    }
    else {
      sprintf(info_yml_file,"%s", INFO_YML_FILENAME);
    }

    read_cospar_biascor(info_yml_file, &HD->cospar_biasCor);

    // for each data event, store mu(z,cospar_biascor)
    for(irow=0; irow < NSN; irow++ ) {
      ztmp    = HD->z[irow] ;
      rz      = codist(ztmp, &HD->cospar_biasCor);
      mu_cos  = get_mu_cos(ztmp, rz);
      HD->mu_cospar_biascor[irow] = mu_cos;
    }
  
  return;
} // end compute_mu_biascor

// ============================================================
void read_cospar_biascor(char *info_yml_file, Cosparam *cospar) {

  // Created Mar 2023
  // Open and read info_yml_file, and scoop up cosmology parameters
  // that were used to create biasCor sim; return cospar.

  FILE *fp ;
  char c_get[MXCHAR_FILENAME];
  char SIMKEY_LIST[NSIMKEY_COSPAR][40] = {
    SIMKEY_MUSHIFT, SIMKEY_OMEGA_MATTER, SIMKEY_OMEGA_LAMBDA,
    SIMKEY_w0_LAMBDA, SIMKEY_wa_LAMBDA } ;
  char SIMKEY_FOUND_STRING[100] = "" ;

  char fnam[] = "read_cospar_biascor" ;

  // --------------- BEGIN --------------

  fp = fopen(info_yml_file,"rt");
  if ( !fp ) {
    sprintf(c1err,"Cannot open INFO file") ;
    sprintf(c2err,"%s", info_yml_file );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  printf("\n   Read cospar_biascor from \n\t %s \n", info_yml_file);
  fflush(stdout);

  cospar->omm = cospar->ome = cospar->w0 = cospar->wa = -9.0 ;
  cospar->mushift = 0.0 ;

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    if ( strcmp(c_get,SIMKEY_MUSHIFT) == 0 ) { 
      fscanf(fp, "%le", &cospar->mushift ); 
      catVarList_with_comma(SIMKEY_FOUND_STRING,SIMKEY_MUSHIFT);
    }  

    else if ( strcmp(c_get,SIMKEY_OMEGA_MATTER) == 0 ) {
      fscanf(fp, "%le", &cospar->omm ); 
      catVarList_with_comma(SIMKEY_FOUND_STRING,SIMKEY_OMEGA_MATTER);
    }  

    else if ( strcmp(c_get,SIMKEY_OMEGA_LAMBDA) == 0 ) {
      fscanf(fp, "%le", &cospar->ome ); 
      catVarList_with_comma(SIMKEY_FOUND_STRING,SIMKEY_OMEGA_LAMBDA);
    }  

    else if ( strcmp(c_get,SIMKEY_w0_LAMBDA) == 0 ) {
      fscanf(fp, "%le", &cospar->w0 ); 
      catVarList_with_comma(SIMKEY_FOUND_STRING,SIMKEY_w0_LAMBDA);
    }  

    else if ( strcmp(c_get,SIMKEY_wa_LAMBDA) == 0 ) {
      fscanf(fp, "%le", &cospar->wa ); 
      catVarList_with_comma(SIMKEY_FOUND_STRING,SIMKEY_wa_LAMBDA);
    }  
    
  } // end while

  fclose(fp);

  printf("   Found cosmology parameters used to generate sim biasCor: \n");
  printf("      OM,OL,w0,wa,mushift = %.3f %.3f %.3f %.3f %.3f \n",
	 cospar->omm, cospar->ome, cospar->w0,  cospar->wa, cospar->mushift );
  
  // check for missing cospar
  int j, NERR=0;
  char *key ;
  for(j=0; j < NSIMKEY_COSPAR; j++ ) {
    key = SIMKEY_LIST[j];
    if ( strstr(SIMKEY_FOUND_STRING,key) == NULL ) {
      printf(" ERROR: missing required cospar key '%s' \n", key);
      NERR++ ;
    }
  }

  if ( NERR > 0 ) {
    sprintf(c1err,"%d missing cospar keys; see above.", NERR) ;
    sprintf(c2err,"check %s ", INFO_YML_FILENAME );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  fflush(stdout);

  return;

} // end read_cospar_biascor

// ====================================
bool read_ISDATA_REAL(char *inFile) {
  // Created July 19 2022 
  // read comment-header lines to search for
  //   # ISDATA_REAL: <val>
  // If ISDATA_REAL is not found, abort 

  FILE *fp;
  int gzipFlag;
  bool FOUND_KEY_ISDATA = false, ISDATA_REAL ;
  char locFile[MXPATHLEN];
  char fnam[] = "read_ISDATA_REAL" ;

  // ------------- BEGIN ------------

  int OPENMASK = OPENMASK_VERBOSE + OPENMASK_IGNORE_DOCANA ;
  fp = snana_openTextFile(OPENMASK, "", inFile,
			  locFile, &gzipFlag );

  if ( !fp ) {
    sprintf(c1err,"Cannot open HD file") ;
    sprintf(c2err,"%s", inFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  
  char c_get[200];
  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"VARNAMES:") == 0 ) { break; }

    if ( strcmp(c_get,KEYNAME_ISDATA_REAL) == 0 )  { 
      fscanf(fp, "%d", &ISDATA_REAL);  
      FOUND_KEY_ISDATA = true;
      break; 
    }
  } // end while

  // - - - - 
  if ( FOUND_KEY_ISDATA ) {
    // adjust blind flag based on real or sim data
    char ctype[40];
    if ( ISDATA_REAL ) 
      { INPUTS.blind = true; sprintf(ctype,"blind REAL"); }
    else
      { INPUTS.blind = false; sprintf(ctype,"unblind SIM"); }
    
    printf("\t Found ISDATA_REAL = %d --> %s data\n", 
	   ISDATA_REAL, ctype);
  }
  else {
    sprintf(c1err,"Could not find required ISDATA_REAL in HD file.") ;
    sprintf(c2err,"Either add flag in HD or re-run create_cov stage." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);   
  }

  fclose(fp);

  return ISDATA_REAL;

} // end read_ISDATA_REAL

// =================================================================
void read_mucov(char *inFile, int imat, COVMAT_DEF *MUCOV ){

  // Created October 2021 by A.Mitra and R.Kessler
  // Read option Cov_syst matrix, and invert it.
 
#define MXSPLIT_mucov 20

  int MSKOPT_PARSE = MSKOPT_PARSE_WORDS_STRING+MSKOPT_PARSE_WORDS_IGNORECOMMA;
  int NSN_STORE    = HD_LIST[imat].NSN; // number passing cuts
  int NSN_ORIG     = HD_LIST[imat].NSN_ORIG; // total number read from HD file
  int NDIM_STORE   = NSN_STORE ;
  int NMAT_READ_UPDATE = 5000000;  // 5 million
  
  char ctmp[200], SN[2][12], locFile[1000] ;
  int NSPLIT, NROW_read=0, NDIM_ORIG = 0, NMAT_ORIG=0 ;
  int NMAT_read=0,  NMAT_store = 0;
  float f_MEM;
  double cov, XM, XMTOT, dt_read;
  int N, N0, N1, i0, i1, j, iwd, NWD, i, k0, k1, kk,gzipFlag ;
  char  **ptrSplit, covtype[60];
  FILE *fp;
  bool UPDSTD;
  time_t t_start_read, t_read;
  char fnam[] = "read_mucov" ;

  // ---------- BEGIN ----------------

  MUCOV->N_NONZERO = 0;

  printf("\n# ======================================= \n");

  if ( INPUTS.use_mucov == 1 ) 
    { sprintf(covtype, "MUCOVSYS");  }
  else if ( INPUTS.use_mucov == 2 )
    { sprintf(covtype, "MUCOVTOT^{-1}");  }    

  printf("  Process %s file \n", covtype);   fflush(stdout);
  sprintf(MUCOV->fileName, "%s", inFile);
  
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
  t_start_read = time(NULL);
  
  while ( fgets(ctmp, 100, fp) != NULL ) {
    // ignore comment lines 
    if ( commentchar(ctmp) ) { continue; }
    
    // break line into words
    splitString(ctmp, " ", fnam, MXSPLIT_mucov,  // (I)
		&NSPLIT, ptrSplit);       // (O)

    if ( NROW_read == 0 ) {
      sscanf(ptrSplit[0],"%d",&NDIM_ORIG);
      printf("\t Found COV dimension %d\n", NDIM_ORIG);
      printf("\t Store COV dimension %d\n", NDIM_STORE);
      
      if ( NDIM_ORIG != HD_LIST[imat].NSN_ORIG ) {
	sprintf(c1err,"NDIM(COV)=%d does not match NDIM(HD)=%d ??",
		NDIM_ORIG, HD_LIST[imat].NSN_ORIG);
	sprintf(c2err,"Above NDIM are before cuts.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
      }

      MUCOV->NDIM    = NDIM_ORIG; // read entire COV without cuts
      malloc_COVMAT(+1,MUCOV);

      NMAT_ORIG = NDIM_ORIG * NDIM_ORIG ;
      XMTOT = (double)(NMAT_ORIG) * 1.0E-6 ;
    }
    else {
      // store entire cov matrix (no cuts yet)
      sscanf( ptrSplit[0],"%le",&cov);
      MUCOV->ARRAY1D[NMAT_read] = cov;
      NMAT_read++ ;
      UPDSTD = ( (NMAT_read % NMAT_READ_UPDATE) == 0 ||
		 NMAT_read == NMAT_ORIG);
      if ( UPDSTD ) {
	t_read  = time(NULL);
	dt_read = t_read - t_start_read;
	XM  = (double)NMAT_read * 1.0E-6 ;
	printf("\t Finished reading %.3f of %.3f million "
	       "%s elements (%.0f seconds)\n",
	       XM, XMTOT, covtype, dt_read );
	fflush(stdout);
      }
    }
    
    NROW_read += 1 ;

  } // end of read loop


  // - - - - - - - - - - -  -
  // Trim cov matrix using pass_cuts; MUCOV->ARRAY1D gets overwritten
  int NDIM_CHECK = applyCut_COVMAT(HD_LIST[imat].pass_cut, MUCOV);

  NMAT_store = NDIM_STORE*NDIM_STORE;

  printf("\t Read %d non-zero %s elements in %.0f seconds.\n",
	 MUCOV->N_NONZERO, covtype, dt_read );
  fflush(stdout);

  // if all COV elements are zero, this is a request for stat-only fit,
  // so disable cov
  if ( MUCOV->N_NONZERO == 0 ) { 
    printf("\t -> disable off-diag COV computations.\n");
    fflush(stdout);
    INPUTS.use_mucov = 0; 
    INPUTS.USE_SPEED_OFFDIAG = false; // disable speed flag for approx min chi2 
    return; 
  }

  // - - - - - - - - - - - - - - - - -
  // sanity check
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


  if ( INPUTS.use_mucov == 1 ) {
    // Add diagonal errors (from Hubble diagram) to covsys
    double COV_STAT ;
    for ( i=0; i<NDIM_STORE; i++ )  {
      kk = i*NDIM_STORE + i;
      COV_STAT = HD_LIST[imat].mu_sqsig[i] ;
      MUCOV->ARRAY1D[kk] += COV_STAT ;
    }
  } 

  
  return ; 
}
// end of read_mucov



// ===================================
int applyCut_HD(bool *PASS_CUT_LIST, HD_DEF *HD) {

  // Use PASS_CUT_LIST to update HD arrays for elements passing cuts.
  //
  // Oct 31 2024:
  //  + adjust mu_ref  ... bug fixes
  //  + adjust new mu_sim
  
  int NSN_ORIG = HD->NSN_ORIG;
  int irow, NSN_STORE=0;
  double ztmp, mu_sig ;
  char *cid;
  char fnam[] = "applyCut_HD";

  // ------------ BEGIN ------------

  for(irow = 0; irow < NSN_ORIG; irow++ ) {

    if ( !PASS_CUT_LIST[irow] ) { continue; }

    ztmp        = HD->z[irow];
    mu_sig      = HD->mu_sig[irow];
    cid         = HD->cid[irow] ;
    
    sprintf(HD->cid[NSN_STORE], "%s", cid );
    HD->mu[NSN_STORE]       = HD->mu[irow];      
  
    if ( INPUTS.muerr_force > 0.0 ) { mu_sig = INPUTS.muerr_force; }    
    HD->mu_sig[NSN_STORE]   = mu_sig;
    HD->mu_sqsig[NSN_STORE] = mu_sig*mu_sig ;

    HD->z[NSN_STORE]        = ztmp;
    HD->logz[NSN_STORE]     = log10(ztmp); 
    HD->z_sig[NSN_STORE]    = HD->z_sig[irow];

    HD->mu_ref[NSN_STORE]   = HD->mu_ref[irow];
    HD->mu_sim[NSN_STORE]   = HD->mu_sim[irow];


    if ( INPUTS.USE_HDIBC )
      {  HD->mu_cospar_biascor[NSN_STORE] = HD->mu_cospar_biascor[irow]; }
    
    if ( ztmp < HD->zmin ) { HD->zmin = ztmp; }
    if ( ztmp > HD->zmax ) { HD->zmax = ztmp; }

    NSN_STORE++ ;
  }


  HD->NSN = NSN_STORE ;

  return(NSN_STORE);

} // end applyCut_HD

// =========================================================
int applyCut_COVMAT(bool *PASS_CUT_LIST, COVMAT_DEF *MUCOV) {

  // Created Mar 2023
  // Overwrite input MUCOV->ARRAY1D with matrix containing only
  // the elements passing 1D-array of cuts in PASS_CUT_LIST.
  // Function Return final size of matrix = number of True elements in 
  // PASS_CUT_LIST.

  int  i0, i1, k0, k1, kk, j, NDIM_ORIG = MUCOV->NDIM;
  int  NDIM_STORE = 0 ;
  double cov;
  char fnam[] = "applyCut_COVMAT" ;

  // ----------- BEGIN --------------

  i0 = i1 = k0 = k1 = 0 ;
  MUCOV->N_NONZERO = 0 ;

  // first find out how many pass cuts
  for(j=0; j < NDIM_ORIG; j++ ) {
    if ( PASS_CUT_LIST[j] ) { NDIM_STORE++ ; }
  }

  for(j=0; j < NDIM_ORIG*NDIM_ORIG; j++ ) {
    cov = MUCOV->ARRAY1D[j];
    if( PASS_CUT_LIST[i0] && PASS_CUT_LIST[i1] ) {
	kk = k1*NDIM_STORE + k0;
	MUCOV->ARRAY1D[kk] = cov;
	if ( cov != 0.0 ) { MUCOV->N_NONZERO++ ; }
	k0++;
	if ( k0 == NDIM_STORE ) { k1++; k0=0; }
      }
      i0++;
      if ( i0 == NDIM_ORIG ) { i1++; i0=0; }
  }

  MUCOV->NDIM = NDIM_STORE;
  return(NDIM_STORE) ;

} // end applyCut_COVMAT


// ==================================
void dump_MUCOV(COVMAT_DEF *MUCOV, char *comment ) {

  char fnam[]="dump_MUCOV";
  int i0, i1, NROW, kk ;
  int MAX_ROW = INPUTS.ndump_mucov ;
  int NSN     = MUCOV->NDIM ;
  if (NSN < MAX_ROW ) { NROW = NSN; }   else { NROW =  MAX_ROW; }
  
  // dump
  printf("\n DUMP %s \n", comment);
  
  for (i0=0; i0 < NROW; i0++)  {
      for(i1=0; i1 < NROW; i1++) {
	kk = i0*NROW + i1;
	printf("%9.5f ", MUCOV->ARRAY1D[kk] );
      }
      printf("\n");
  }
  printf("\n");
  fflush(stdout);
  return;
}// end of dump_MUCOV

// ==================================
void sync_HD_LIST(HD_DEF *HD_REF, 
		  HD_DEF *HD,  COVMAT_DEF *MUCOV) {

  // Created Mar 2023
  // Called for HDIBC method where we have two HDs and possibly two COVMATs.
  // Re-define HD and MUCOV to include only those elements in HD_REF.
  //

  int NSN_REF = HD_REF->NSN ;
  int MEMB    = sizeof(bool) * (NSN_REF + 100);
  bool  *SYNC_LIST = (bool*) malloc(MEMB);
  int ilist0 = 0, ilist1=1, NREMOVE=0 ;
  int iref, i,  idum;
  char *cid;
  char fnam[] = "sync_HD_LIST" ;
  
  // ------------ BEGIN -----------

  match_cid_hash("", -1,0); // reset hash table

  for(iref=0; iref < NSN_REF; iref++ ) {
    cid = HD_REF->cid[iref] ;
    idum = match_cid_hash(cid,ilist0,iref);  // load hash table
  }
  
  // loop over  HD and set SYNC_LIST
  for(i=0; i < HD->NSN; i++ ) {
    cid  = HD->cid[i] ;
    iref = match_cid_hash(cid,ilist1,i);
    SYNC_LIST[i] = false; 
    if ( iref >= 0 )  { SYNC_LIST[i] = true; }
    else              { NREMOVE++ ; }
  }

  printf("\t %s: remove %d SN from HD and MUCOV \n", fnam, NREMOVE);
  fflush(stdout);

  // update HD with CID selection ...
  HD->NSN_ORIG = HD->NSN;
  applyCut_HD(SYNC_LIST, HD);

  // update MUCOV with CID selection ...
  if ( MUCOV->NDIM > 0 ) 
    { applyCut_COVMAT(SYNC_LIST, MUCOV); }

  return ;
} // end sync_HD_LIST

// ===================================================
void sync_HD_redshifts(HD_DEF *HD0, HD_DEF *HD1) {

  // Created Dec 18 2023
  // For binned HD, the redshifts for each HD can be very
  // slightly different because of how BBC computes <z> in each bin.
  // This function forces HD1 redshifts to match those in HD0,
  // and make a slight MU-correction to distance.
  // Report average diff and abort on crazy-large average.
  
  int NSN = HD0->NSN;
  int MEMD = NSN * sizeof(double);
  
  int i, NDIF = 0 ;
  double z0, z1, zdif, zdif_sum=0, zdif_max=0.0, zdif_avg ;
  double rz1_orig, rz1_new, mu1_dif, mu1_orig;
  double *zdif_list   = (double*)malloc(MEMD);
  double *mu1dif_list = (double*)malloc(MEMD);

  bool CHEAT_COSPAR_SIM = false;  // default is false; use only for debugging
  char fnam[] = "sync_HD_redshifts" ;

  // ------------- BEGIN ----------

  if ( NSN > 100 ) { return; } // bail on unbinned HD
  
  print_banner(fnam);
  
  for(i=0; i < NSN; i++ ) {
    z0 = HD0->z[i];
    z1 = HD1->z[i];
    HD0->z_orig[i]   = z0;
    HD1->z_orig[i]   = z1;         
    zdif = z1 - z0 ;

    if ( zdif == 0.0 ) { continue; }
    NDIF++ ;
    
    zdif_sum += zdif;
    if ( fabs(zdif) > fabs(zdif_max) ) { zdif_max = zdif; }

    // compute Mu1 shift corresponding to z1 shift
    if ( CHEAT_COSPAR_SIM ) {
      // cheat mode for testing only
      rz1_orig     = codist(z1, &COSPAR_SIM);
      rz1_new      = codist(z0, &COSPAR_SIM);
    }
    else {
      // nominal default is to use cospar from biasCor.
      rz1_orig     = codist(z1, &HD1->cospar_biasCor);
      rz1_new      = codist(z0, &HD1->cospar_biasCor);
    }
    
    mu1_orig     = HD1->mu[i];
    
    mu1_dif = (get_mu_cos(z0, rz1_new) - get_mu_cos(z1, rz1_orig) ) ;
    
    HD1->z[i]   = z0;  // force same redshift
    HD1->mu[i] += mu1_dif;
    
    zdif_list[i]   = zdif;
    mu1dif_list[i] = mu1_dif;
  }

  // - - - - - - - - - 
  // print z & mu change for each redshift
  if ( NDIF > 0 ) {
    for(i=0; i < NSN; i++ ) {
      z0 = HD0->z[i];
      printf("\t z=%.5f: dz1 = %9.6f -> mu1 += %8.5f\n",
	     z0, zdif_list[i], mu1dif_list[i] );
      fflush(stdout);
    }
  } // end NDIF
  
  zdif_avg = zdif_sum / (double)NSN;
  printf("\t <zdif>= %le  zdif_max=%le\n", fnam, zdif_avg, zdif_max);
  fflush(stdout);

  free(zdif_list); free(mu1dif_list);
  
  //  debugexit(fnam);
  
  return;
} // sync_HD_redshifts

// ===================================
void compute_MUCOV_FINAL(void) {

  // Ceated Mar 2023
  // If there is just one COV, then MUCOV_FINAL = MUCOV[0].
  // If there are two COVs, MUCOV_FINAL is the average for each element.
 
  int  j, icov ;
  int  NSN    = HD_LIST[0].NSN;
  int  NMUCOV = INPUTS.NMUCOV; // 1 or 2
  if ( NMUCOV == 0 ) { return; }
  double fac = 1.0 / (double)NMUCOV ;
  double covsum ;

  // ---------- BEGIN ------------

  WORKSPACE.MUCOV_FINAL.NDIM = NSN ;
  malloc_COVMAT(+1, &WORKSPACE.MUCOV_FINAL);
  
  for(j=0; j < NSN*NSN; j++ ) {
    covsum = 0.0 ;
    for(icov = 0; icov < INPUTS.NMUCOV; icov++ )
      { covsum += WORKSPACE.MUCOV[icov].ARRAY1D[j]; }

    WORKSPACE.MUCOV_FINAL.ARRAY1D[j] = fac * covsum ;
  }

    return ;
} // end compute_MUCOV_FINAL

//===================================
void set_priors(void) {

  char *comment ;
  bool noprior = true;
  char fnam[]="set_priors";

  // =========== BEGIN ============
  
  if ( INPUTS.use_bao )  { init_bao_prior(1); }

  if ( INPUTS.use_cmb )  { init_cmb_prior(1); }

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
    comment = BAO_PRIOR.comment ; 
    printf("   %s\n", comment) ;

    comment = BAO_PRIOR.comment_rd_prior ; 
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
    comment = CMB_PRIOR.comment ; 
    printf("   %s\n", comment) ; 
  }

  if ( noprior ) 
    { printf("\t None.\n"); }
  else {
    printf("\t OM,w0,wa for priors: %.3f  %.3f  %.3f \n",
	   INPUTS.OMEGA_MATTER_SIM, INPUTS.w0_SIM, INPUTS.wa_SIM );
  }

  // print a few computed mu resids at simCospar values
  double z, rz, mu_sim, mu_lcdm;
  printf("\n   Ref sim distance resids:\n");
  for (z=0.5; z <=3.0; z+=0.5) {
    rz      =  codist( z, &COSPAR_SIM) ;
    mu_sim  =  get_mu_cos(z,rz) ;
    rz      =  codist( z, &COSPAR_LCDM) ;
    mu_lcdm =  get_mu_cos(z,rz) ;
    printf("\t z=%.2f  mu(sim) - mu(LCDM) = %.3f - %3f = %.3f \n",
	   z, mu_sim, mu_lcdm, mu_sim-mu_lcdm );
  }
  
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
  double w0 = INPUTS.w0_SIM ;
  double wa = INPUTS.wa_SIM ;
  
  char *comment = CMB_PRIOR.comment;
  char fnam[] = "init_cmb_prior" ;

  // ---------- BEGIN ------------

  if ( OPT == -9 ) {
    CMB_PRIOR.z      = 1089.0 ;
    CMB_PRIOR.a      = 1.0/(1.0 + CMB_PRIOR.z) ;
    CMB_PRIOR.R      = -9.0 ;
    CMB_PRIOR.sigR   = .019 ;
    CMB_PRIOR.ranseed_R = 0 ;
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
    double Rcmb_ranshift = get_Rcmb_ranshift();
    char str_ranshift[40]; str_ranshift[0] = 0 ;

    a = CMB_PRIOR.a ;  // 1/(1+z)
    set_HzFUN_for_wfit(ONE, OM, OE, w0, wa, &HzFUN) ;
    rz = Hainv_integral ( a, ONE, &HzFUN ) / LIGHT_km;
    CMB_PRIOR.R = sqrt(OM) * rz ;

    if ( Rcmb_ranshift != 0.0 ) {
      CMB_PRIOR.R += Rcmb_ranshift;
      sprintf(str_ranshift,"(Rcmb_ranshift=%.4f)", Rcmb_ranshift);
    }
    sprintf(comment, "CMB sim-prior:  R=%5.3f +- %5.3f  %s" ,
	    CMB_PRIOR.R, CMB_PRIOR.sigR, str_ranshift );
  }

  
  if ( INPUTS.omm_min < 0.0 ) {
    sprintf(c1err,"Negative omm (ommin=%.3f) not allowed with CMB prior",
	    INPUTS.omm_min);
    sprintf(c2err,"Remove CMB prior or set omm_min >=0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);

  }

  return;

} // end init_cmb_prior


// ==============================
double get_Rcmb_ranshift(void) {

  // Created Oct 17 2023
  // Return random Rcmb shift (if user-selects this option)
  // The final Rcmb is determined externally from this function.
  
  int NSTREAM = 1, istream=0 ;
  int iseed;
  double gran, ranshift = 0.0 ;
  char fnam[] = "get_Rcmb_ranshift";

  // ---------- BEGIN -------------

  if (CMB_PRIOR.ranseed_R <= 0) { return ranshift; }

  if ( CMB_PRIOR.ranseed_R == 1 ) {
    // internally compute unique seed based on HD values
    iseed = ISEED_UNIQUE();
  }
  else {
    // seed is explicitly passed by user
    iseed = CMB_PRIOR.ranseed_R; 
  }

  init_random_seed(iseed, NSTREAM);
  gran      = unix_getRan_Gauss(istream);
  ranshift  = gran * CMB_PRIOR.sigR;

  printf(" %s: random Rcmb shift = %.4f (ISEED=%d, Gaussran=%.4f)\n",
	 fnam, iseed, ranshift, iseed, gran); fflush(stdout);

  return ranshift;

} // end get_Rcmb_ranshift


// ===============================
int ISEED_UNIQUE(void) {

  // Created Oct 17 2023.
  // return ISEED that is unique for each data set.
  // Used for user-input options:
  //   -ranseed_Rcmb 1
  //   -muerr_ideal <muerr>
  //
  HD_DEF *HD = &HD_LIST[0];
  int ISEED  = 0 ;
  int k, NUSE=0, NUSE_MAX = 5;
  double product = 1.0 ;
  char fnam[] = "ISEED_UNIQUE" ;

  // --------- BEGIN --------

  k = 0 ;
  while ( NUSE < NUSE_MAX &&  k < HD->NSN-1 ) {

    // make sure that z-bin has a real MU value 
    if ( HD->mu_sig[k] < 8.0 ) {
      product *= ( HD->mu[k]/40.0 ) ;  // ~ unity 
      NUSE++ ;
    }
    k++ ;
  }

  // product is of order unity, so multiply by a very big number
  ISEED = (int)(product * 8224532.0);
  
  return ISEED;

} // end ISEED_UNIQUE

// =========================================
void init_bao_prior(int OPT) {

  // Created Oct 2021 by R.Kessler
  // Re-written April 2024 to read BAO params from data file
  //
  // OPT= -9 --> set all params to -9
  //               (before reading user input)
  // OPT= +1 --> set params to measured or sim values.
  //               (after reading user input)
  //

  int i, IVAR ;
  char *comment = BAO_PRIOR.comment;
  char bao_file[MXPATHLEN], VARNAME[40];
  double z, MEAN, COV;
  FILE *fp;
  int LDMP = 0 ;

  double frac_H0err   = (double)H0err_Planck/(double)H0_Planck ;
  double frac_rderr   = (double)rderr_Planck/(double)rd_Planck ;
  double sqfrac_H0err = frac_H0err * frac_H0err ;
  double sqfrac_rderr = frac_rderr * frac_rderr ;
  double mm;
  
  char fnam[] = "init_bao_prior" ;

  // -----------------BEGIN ------------------


  if ( OPT == -9 ) {
    for(i=0; i < MXDATA_BAO_PRIOR; i++ ) {
      BAO_PRIOR.REDSHIFT[i]   = -9.0 ;
      BAO_PRIOR.MEAN[i]       = -9.0 ;
      BAO_PRIOR.VARNAME[i][0] =  0  ;
      BAO_PRIOR.IVAR[i]       = -9;
    }
    BAO_PRIOR.comment[0] = 0;    
    return;
  } 
  
  // - - - - - - - - - - - - -
  // read redshift, means, quantity 
  sprintf(PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
  sprintf(bao_file, "%s/models/BAO/%s_mean.dat", PATH_SNDATA_ROOT, INPUTS.bao_sample);
  fp = fopen(bao_file,"rt");
  if ( !fp ) {
    sprintf(c1err,"Could not open %s", bao_file);
    sprintf(c2err,"Check -bao argument");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  char **ptrSplit, cline[200];
  int j, NWD, ndata=0;
  ptrSplit = (char **)malloc(MXDATA_BAO_PRIOR*sizeof(char*));
  for(i=0; i<MXDATA_BAO_PRIOR; i++) { ptrSplit[i]=(char *)malloc(40*sizeof(char));  }

  while ( fgets(cline, 200, fp) != NULL ) {
    // ignore comment lines 
    if ( commentchar(cline) ) { continue; }
    
    // break line into words
    splitString(cline, " ", fnam, MXDATA_BAO_PRIOR,  // (I)
		&NWD, ptrSplit);       // (O)

    sscanf(ptrSplit[0], "%le", &z );
    sscanf(ptrSplit[1], "%le", &MEAN );
    sscanf(ptrSplit[2], "%s",  VARNAME );

    if (strstr(VARNAME,"DM_over") != NULL ) 
      { IVAR = IVAR_BAO_DM_over_rd; }
    else if (strstr(VARNAME,"DH_over") != NULL ) 
      { IVAR = IVAR_BAO_DH_over_rd; }
    else if (strstr(VARNAME,"DV_over") != NULL ) 
      { IVAR = IVAR_BAO_DV_over_rd; }
    else {
      sprintf(c1err,"Could not decode BAO variable '%s'", VARNAME);
      sprintf(c2err,"from %s", bao_file);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
    }

    // store bao info for this data point
    BAO_PRIOR.REDSHIFT[ndata] = z;
    BAO_PRIOR.MEAN[ndata]     = MEAN;
    BAO_PRIOR.IVAR[ndata]     = IVAR;
    sprintf(BAO_PRIOR.VARNAME[ndata],"%s", VARNAME);   
    ndata++ ;   
  }
  
  BAO_PRIOR.NDATA = ndata;
  fclose(fp);

  
  // next read cov matrix
  int MEMD  = ndata * ndata * sizeof(double) ;
  BAO_PRIOR.COV1D    = (double*) malloc(MEMD);
  BAO_PRIOR.COVINV1D = (double*) malloc(MEMD);
  
  sprintf(bao_file, "%s/models/BAO/%s_cov.dat", PATH_SNDATA_ROOT, INPUTS.bao_sample);
  fp = fopen(bao_file,"rt");
  if ( !fp ) {
    sprintf(c1err,"Could not open %s", bao_file);
    sprintf(c2err,"Check -bao argument");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  int nrow = 0, ncov=0 ;
  double cov_planck;
  while ( fgets(cline, 200, fp) != NULL ) {
    // ignore comment lines 
    if ( commentchar(cline) ) { continue; }
    
    // break line into words
    splitString(cline, " ", fnam, MXDATA_BAO_PRIOR,  // (I)
		&NWD, ptrSplit);       // (O) 

    if ( NWD != ndata ) {
      sprintf(c1err,"Expected %d cov elements in row %d, but found %d",
	      ndata, nrow, NWD);
      sprintf(c2err,"Check %s", bao_file);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
    }

    
    for(i=0; i < NWD; i++ ) {
      // H0 & rd error is an external (Planck) error; here we tack it onto
      // BAO measurement error instead of floating these parameters in the fit.
      if ( OPT_RD_CALC == 0 ) {
	mm          = BAO_PRIOR.MEAN[i] * BAO_PRIOR.MEAN[nrow] ;
	cov_planck  = (sqfrac_H0err + sqfrac_rderr) * mm ;
      }
      else {
	cov_planck = 0.0 ;
      }
      
      sscanf(ptrSplit[i], "%le", &BAO_PRIOR.COV1D[ncov] );
      BAO_PRIOR.COV1D[ncov]   += cov_planck;
      BAO_PRIOR.COVINV1D[ncov] = BAO_PRIOR.COV1D[ncov]  ; // invert below
      ncov++ ;
    }
    nrow++ ;
  
  } // end while

  invertMatrix( ndata, ndata, BAO_PRIOR.COVINV1D );
  sprintf(comment,"BAO prior from %s data (n=%d)", INPUTS.bao_sample, ndata );
  
  // - - - - - -
  if ( INPUTS.use_bao == 2 ) { 
    double rd ; 
    HzFUN_INFO_DEF HzFUN;
    rd = rd_bao_prior(OPT_RD_CALC, &COSPAR_SIM);
    
    //printf(" xxx %s rd  = %f \n", fnam, rd);
    for (i=0; i < ndata; i++ ) {
      z    = BAO_PRIOR.REDSHIFT[i];
      IVAR = BAO_PRIOR.IVAR[i];

      if ( IVAR == IVAR_BAO_DM_over_rd ) 
	{ MEAN = DM_bao_prior(z, &COSPAR_SIM) / rd ; }
      else if ( IVAR == IVAR_BAO_DH_over_rd )
	{ MEAN = DH_bao_prior(z, &COSPAR_SIM) / rd ; }
      else if ( IVAR == IVAR_BAO_DV_over_rd )
      	{ MEAN = DV_bao_prior(z, &COSPAR_SIM) / rd ; }
      
      BAO_PRIOR.MEAN[i]    = MEAN ;
    }
    sprintf(comment,"BAO prior from %s using sim cosmology (n=%d)",
	    INPUTS.bao_sample, ndata );
  }

  sprintf(BAO_PRIOR.comment_rd_prior,"rd Gauss prior: mean=%.2f  sig=%.2f\n",
	  INPUTS.rd_prior, INPUTS.rd_prior_sig);

  if ( LDMP ) { dump_bao_prior(); }
  //  debugexit(fnam);
  
  return ;
  
} // end init_bao_prior

void dump_bao_prior(void) {

  double z, mean, err, cov;
  char *varname;
  int i, kk, NDATA = BAO_PRIOR.NDATA;
  char fnam[] = "dump_bao_prior";
  // --------- BEGIN ---------

  print_banner(fnam);
  printf("\t xxx %s\n", BAO_PRIOR.comment);
  for(i=0; i < NDATA; i++ ) {
    z        = BAO_PRIOR.REDSHIFT[i];
    mean     = BAO_PRIOR.MEAN[i];
    varname  = BAO_PRIOR.VARNAME[i];

    kk = i*NDATA + i; // diagonal element of COV
    cov      = BAO_PRIOR.COV1D[kk] ;
    err      = sqrt(cov);
    
    printf("\t xxx %s:  z=%5.3f   mean(%s) = %7.4f +- %7.4f \n",
	   fnam, z, varname, mean, err);
    fflush(stdout);
  }

  return;
} // end dump_bao_prior


double rd_bao_prior(int OPT, Cosparam *cpar) {

  // OPT = 0 -> Planck value
  // OPT = 1 -> compute integral with constant c_s
  // OPT = 2 -> compute integral with c_s(z)
  
  // Compute drag epochs rd.
  // rd ~ 150 Mpc                  Pg. 9, Aubourg et al. [1411.1074]
  // rd = int_0^inf [c_s /H(z)];                       Eq. 13, Alam 2020.
  // c_s= c_light/sqrt[3*(1 + 3*Om_b / 4*Om_gamma)]    Davis et al, Page 4.
  // Om_b ~ 0.02/h^2                                   Davis T. Note
  // Om_g ~ 2.469 * 10e-5 * T_CMB / 2.725              Davis et al, after eq. 15

  double z_d    = 1060.0 ;  
  double H0     = H0_Planck;

  int OPT_c_sound  = 1;  // 1 = simple approx, 2= Om_[b,g] dependence
  double amin   = 1.0E-6,  amax = 1.0/(1.0+z_d);
  double Hinv_integ, Einv_integ, c_s, rd = 1.0;
  HzFUN_INFO_DEF HzFUN_INFO;
  int  LDMP = 0;
  char fnam[] = "rd_bao_prior" ;
  
  // ---------- BEGIN ----------  

  if ( OPT == 0 ) {
    rd = rd_Planck ;
  }
  
  else if ( OPT == 1 ) {
    // approx c_s = contant
    c_s = c_sound(OPT_c_sound, z_d, H0_Planck);
    Einv_integ = Eainv_integral(amin, amax, cpar);
    Hinv_integ = Einv_integ / H0 ;
    rd = c_s * Hinv_integ  ; 
  }
  else {
    // full integral
    rd = cs_over_E_integral(amin, amax, cpar) / H0;
  }
  
  return rd;
}


// ======================================
double c_sound(int OPT, double z, double H0) {

  // OPT=1 -> simplest approx
  // OPT=2 -> include dependence om Om_b(z)/Om_g(z)
  
  double c_s;
  double root3_inv = 0.57735 ;  // 1/sqrt(3)
  double Om_b0     = 0.0225;    // Omega_baryon, today (1/h^2 cancels with Om_g)
  double Om_g0     = 2.47E-5;   // Omega_gamma, today    
  double z1        = 1.0 + z;
  double a         = 1.0/z1 ;
  char fnam[] = "c_sound" ;

  // ---------- BEGIN ----------

  c_s  = c_light * root3_inv ; // common factors for both methods
  
  if (  OPT == 1 ) {
    c_s *= 0.9 ; // Davis internal note
  }
  else if ( OPT == 2 ) {
    c_s /= sqrt( 1.0 + 0.75*a*(Om_b0/Om_g0) )  ;
  }
  else {
    sprintf(c1err,"Invalid OPT = %d", OPT);
    sprintf(c2err,"Valid OPT = 1 or 2");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  return c_s;
  
} // end c_cound

double DM_bao_prior(double z, Cosparam *cpar){
  // Ayan Mitra Nov, 2021
  double DM = 1.0, H0 = H0_Planck;
  double amin   = 1.,  amax = 1.0/(1+z);
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


double DV_bao_prior(double z, Cosparam *cpar){
  double DM, DH, DV;
  DM = DM_bao_prior(z,cpar);
  DH = DH_bao_prior(z,cpar);
  DV = cbrt( z * DM*DM * DH ) ;
  return DV;
}



// =========================================
void set_stepsizes(void) {

  // Created Oct 1 2021
  // Compute grid size for each cosmoPar dimension

  char fnam[] = "set_stepsizes";

  // ------------ BEGIN ------------

  INPUTS.w0_stepsize = INPUTS.wa_stepsize =  INPUTS.omm_stepsize = 0 ;

  if (INPUTS.w0_steps  > 1) 
    { INPUTS.w0_stepsize  = (INPUTS.w0_max-INPUTS.w0_min)/(INPUTS.w0_steps-1);}
  if (INPUTS.wa_steps  > 1) 
    { INPUTS.wa_stepsize  = (INPUTS.wa_max-INPUTS.wa_min)/(INPUTS.wa_steps-1);}
  if (INPUTS.omm_steps > 1) 
    { INPUTS.omm_stepsize = (INPUTS.omm_max-INPUTS.omm_min)/(INPUTS.omm_steps-1);}
  

  printf("\n");
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

  printf("   ------------------------------------\n");

  fflush(stdout);

  return;
} // end set_stepsizes


// ==================
void set_Ndof(void) {
  // Compute number of degrees of freedom for wfit
  int Ndof;
  int Ndof_data = HD_LIST[0].NSN - 3 ; // h, w, om = 3 fitted parameters
  int Ndof_prior = 0;
  if ( INPUTS.dofit_w0wa  ) { Ndof_data-- ; } 
  if ( INPUTS.use_bao     ) { Ndof_prior += BAO_PRIOR.NDATA ; }
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
void init_rz_interp(HD_DEF *HD) {

  // Created Apr 22 2022
  // For large SN sample, initialize lists to interpolate rz(z) in chi2 fun.
  double zmin = HD->zmin;
  double zmax = HD->zmax;
  int    NSN  = HD->NSN;
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
  HD_DEF *HD0 = &HD_LIST[0];
  HD_DEF *HD1 = &HD_LIST[1];
  char fnam[] = "exec_rz_interp";

  // ----------------- BEGIN --------------

  logz = HD0->logz[k]; 
  if ( INPUTS.USE_HDIBC ) { logz = 0.5*(HD0->logz[k] + HD1->logz[k]); }

  iz   = (int)((logz - logz_min)/logz_bin );
  if ( iz < 0 || iz >= n_logz ) {
    sprintf(c1err,"Invalid iz=%d (range is 0 to %d)", iz, n_logz);
    sprintf(c2err,"z=%f  logz=%f  for CID=%s", 
	    HD0->z[k], HD0->logz[k], HD0->cid[k] );
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
    rz_loc    = WORKSPACE.rz_list_interp[iz];    // last z-bin
    mucos_loc = WORKSPACE.mucos_list_interp[iz]; // last z-bin
  }
  
  // load output function args
  *rz    = rz_loc;
  *mucos = mucos_loc;

  if ( k == HD0->NSN-1 ) { WORKSPACE.n_exec_interp++ ; }

  // print diagnostic for first few events.
  int LDMP = (WORKSPACE.n_exec_interp < 5) ;
  if ( LDMP ) {
    if ( k==0 ) { WORKSPACE.rz_dif_max = 0.0 ; }    
    rz_exact   = codist(HD0->z[k], cparloc);
    rz_dif     = fabs(rz_loc / rz_exact - 1.0);
    if ( rz_dif > WORKSPACE.rz_dif_max ) 
      { WORKSPACE.rz_dif_max = rz_dif; }  

    if ( k == HD0->NSN-1 ) {
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
  double snchi_tmp, extchi_tmp, mures_tmp ;

  bool UPDATE_STDOUT;
  int  i, kk, j;
  int  imin = -9, kmin = -9, jmin = -9;
  char fnam[] = "wfit_minimize" ;

  // ---------- BEGIN --------------

  int NBTOT = INPUTS.w0_steps * INPUTS.wa_steps * INPUTS.omm_steps;
  int NB=0;

  printf("\n# ======================================= \n");
  printf(" %s: Get prob at %d grid points, and approx mimimized values \n", 
	 fnam, NBTOT );
  printf("\t USE_SPEED_OFFDIAG = %d \n", INPUTS.USE_SPEED_OFFDIAG);
  printf("\t USE_SPEED_INTERP  = %d \n", INPUTS.USE_SPEED_INTERP);
  fflush(stdout);
    
  // prep speed trick
  if ( use_mucov && USE_SPEED_OFFDIAG ) {
    printf("   Prepare speed trick to reduce off-diagonal calculation\n");

    INPUTS.USE_SPEED_OFFDIAG = false; // disable speed flag for approx min chi2 
    cpar_fixed.w0 = -1.0; cpar_fixed.wa=0.0; cpar_fixed.omm=0.3; 
    cpar_fixed.mushift = 0.0 ;
    get_chi2_fit(cpar_fixed.w0,cpar_fixed.wa, cpar_fixed.omm, INPUTS.sqsnrms, 
		 temp0_list, temp1_list, temp2_list,
		 &mures_tmp, &snchi_tmp, &extchi_tmp ); 
    printf("    Very approx chi2min(SNonly,SN+prior: w0=-1,wa=0,omm=0.3) "
	   "= %.1f %.1f \n", snchi_tmp, extchi_tmp);
    fflush(stdout);

    INPUTS.USE_SPEED_OFFDIAG = true ; // restore speed flag
    prep_speed_offdiag(extchi_tmp);
  }

  // - - - - - - - - 
  time_t t0 = time(NULL);  // monitor time to build prob grid

  cpar.mushift = 0.0;

  for( i=0; i < INPUTS.w0_steps; i++){
    cpar.w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    for( kk=0; kk < INPUTS.wa_steps; kk++){    
      cpar.wa = (INPUTS.wa_min + kk*INPUTS.wa_stepsize);
      for(j=0; j < INPUTS.omm_steps; j++){
	cpar.omm = INPUTS.omm_min + j*INPUTS.omm_stepsize; 
	cpar.ome = 1 - cpar.omm;
	
	get_chi2_fit ( cpar.w0, cpar.wa, cpar.omm, INPUTS.sqsnrms, 
		       temp0_list, temp1_list, temp2_list,
		       &mures_tmp, &snchi_tmp, &extchi_tmp ); 

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

	if ( UPDATE_STDOUT || NB==NBTOT ) {
	  char comment[60];
	  sprintf(comment, "chi2 bin %8d of %8d", NB, NBTOT); 
	  print_elapsed_time(t0, comment, UNIT_TIME_SECOND);
	}

      } // j loop
    }  // end of k-loop
  }  // end of i-loop


  // get w,OM at min chi2 by using more refined grid
  // Pass approx w,OM,  then return w,OM at true min

  WORKSPACE.w0_atchimin = INPUTS.w0_min  + imin * INPUTS.w0_stepsize ;
  WORKSPACE.wa_atchimin = INPUTS.wa_min  + kmin * INPUTS.wa_stepsize ;
  WORKSPACE.omm_atchimin= INPUTS.omm_min + jmin * INPUTS.omm_stepsize ;
  
  if ( !INPUTS.use_marg  ) {
    //printf("  Repeat minimization with refined parameter grid \n");
    fflush(stdout);
    WORKSPACE.chi2atmin  = 
      get_minimized_cospar( &WORKSPACE.w0_atchimin, 
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
    for(kk=0; kk < INPUTS.wa_steps; kk++) {
      for(j=0; j < INPUTS.omm_steps; j++) {
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
      for(i=0; i < INPUTS.w0_steps; i++) {
	for(j=0; j < INPUTS.omm_steps; j++) {
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
    WORKSPACE.omm_prob[j] = 0.0 ;
    for(kk=0; kk < INPUTS.wa_steps; kk++) {
      for(i=0; i < INPUTS.w0_steps; i++) {
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
  
  return;

} // end wfit_marginalize


// =========================================
void wfit_uncertainty(void) {

  // Created Apr 26 2024 by R.Kessler
  // Refactor to process each variable exactly the same way
  // without repeating so much code.
  
  printf("\n# ============================================ \n");
  printf(" Estimate uncertainties: \n");
  fflush(stdout);

  wfit_uncertainty_fitpar(varname_omm);
  wfit_uncertainty_fitpar(varname_w);

  if ( INPUTS.dofit_w0wa ) 
    { wfit_uncertainty_fitpar(varname_wa); }

  return ;
}  // end wfit_uncertainty


void wfit_uncertainty_fitpar(char *varname) {

  // Created Apr 2024
  // Compute and store two kinds of fitted uncertainty for parameter 'varname';
  // 1. standard devistion (std) of marginalized posterior (prob_array)
  //     -> symmetric uncertainty
  // 2. 68% confidence interval of marginalized posterior
  //     -> lower and upper uncertainty
  //
  
  double *sig_std, *sig_upper, *sig_lower, sqdelta, delta ;
  double val, val_min, val_max, val_mean, val_step ;
  double val_sum,  probsum, cdf_find, cdf_max ;
  double *val_array, *prob_array, *cdf_array, STD_WARN ;
  int  i, n_steps;
  bool ISVAR_w = false ;
  char fnam[] = "wfit_uncertainty_fitpar";

  // ---------- BEGIN -----------

  if ( strcmp(varname,varname_omm) == 0 ) {
    sig_std   = &WORKSPACE.omm_sig_std ;
    sig_lower = &WORKSPACE.omm_sig_lower ;
    sig_upper = &WORKSPACE.omm_sig_upper ;
    n_steps  =  INPUTS.omm_steps;
    val_min  =  INPUTS.omm_min;
    val_max  =  INPUTS.omm_max;
    val_step =  INPUTS.omm_stepsize;
    val_mean =  WORKSPACE.omm_mean;
    probsum  =  WORKSPACE.omm_probsum;    

    val_array     =  WORKSPACE.omm_val ;
    prob_array    =  WORKSPACE.omm_prob;
    cdf_array     =  WORKSPACE.omm_cdf;
  }
  else if ( strcmp(varname,varname_w) == 0 ) {
    ISVAR_w = true ;
    sig_std   = &WORKSPACE.w0_sig_std ;
    sig_lower = &WORKSPACE.w0_sig_lower ;
    sig_upper = &WORKSPACE.w0_sig_upper ;
    n_steps  =  INPUTS.w0_steps;
    val_min  =  INPUTS.w0_min;
    val_max  =  INPUTS.w0_max;
    val_step =  INPUTS.w0_stepsize;
    val_mean =  WORKSPACE.w0_mean;
    probsum  =  WORKSPACE.w0_probsum;

    val_array      =  WORKSPACE.w0_val ;    
    prob_array     =  WORKSPACE.w0_prob;
    cdf_array      =  WORKSPACE.w0_cdf;
  }
  else if ( strcmp(varname,varname_wa) == 0 ) {
    ISVAR_w = true ;
    sig_std   = &WORKSPACE.wa_sig_std ;
    sig_lower = &WORKSPACE.wa_sig_lower ;
    sig_upper = &WORKSPACE.wa_sig_upper ;
    n_steps  =  INPUTS.wa_steps;
    val_min  =  INPUTS.wa_min;
    val_max  =  INPUTS.wa_max;
    val_step =  INPUTS.wa_stepsize;
    val_mean =  WORKSPACE.wa_mean;
    probsum  =  WORKSPACE.wa_probsum;

    val_array      =  WORKSPACE.wa_val ;        
    prob_array     =  WORKSPACE.wa_prob;
    cdf_array      =  WORKSPACE.wa_cdf;
  }
  else {
    sprintf(c1err,"Invalid varname = '%s'", varname);
    sprintf(c2err,"Valid varnames are %s  %s  %s",
	    varname_omm, varname_w, varname_wa);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);

  }
  
  printf("      %s: \n", varname); fflush(stdout);
  
  *sig_std   = 0.0 ;
  *sig_upper = 0.0 ;
  *sig_lower = 0.0 ;

  // Get std dev for the weighted mean.  This is a reasonable   
  //  measure of uncertainty in w if distribution is ~Gaussian
  sqdelta = 0.0;
  cdf_array[0] = 0.0 ;
  for(i=0; i < n_steps; i++){ 	
    val      = val_min + i*val_step;
    val_array[i]   = val;
    prob_array[i] /= probsum;  
    delta    = val - val_mean ;
    sqdelta += prob_array[i] * (delta*delta);   
    if ( i > 0 ) { cdf_array[i] = cdf_array[i-1] + prob_array[i]; }
  }

  cdf_max  = cdf_array[n_steps-1] ;
  if ( cdf_max < 0.8 ) {
    printf(" WARNING-%s: cdf_max = %le (should be 1.0); check fitpar bounds.",
	   varname, cdf_max);    fflush(stdout);    
    WORKSPACE.NWARN++ ;
    return;
  }
  
  *sig_std = sqrt(sqdelta/probsum);
  printf("\t STD %s_sig estimate = %.4f (use=%d)\n",
	 varname, *sig_std, INPUTS.use_sig_std );

  if ( INPUTS.dofit_w0wa  &&  *sig_std < 0.01 && ISVAR_w ) {
    WORKSPACE.NWARN++ ;
    printf(" WARNING-%s: STD too small, likley BBC problem with MUERR\n", varname);
  }

  
  // - - - -
  cdf_find = (0.5 - PROBSUM_1SIGMA/2.0);
  val  = interp_1DFUN(1, cdf_find, n_steps,  cdf_array, val_array, fnam);
  *sig_lower = val_mean - val;

  cdf_find = (0.5 + PROBSUM_1SIGMA/2.0);
  if (cdf_max < 0.9999 ) {
    printf(" WARNING-%s: cdf_max = %f (not 1); rescale prob\n", varname, cdf_max);
    cdf_find *= cdf_max ;
    WORKSPACE.NWARN++ ;
  }
  val  = interp_1DFUN(1, cdf_find, n_steps,  cdf_array, val_array, fnam);
  *sig_upper = val - val_mean ;
  
  printf("\t 68 percent %s-err estimate: lower = %f, upper = %f (use=%d)\n", 
	 varname, *sig_lower, *sig_upper, INPUTS.use_sig_68);  
  fflush(stdout);

  return;
  
} // end wfit_uncertainty_fitpar


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
  double  omm_sig = WORKSPACE.omm_sig_std;
  double  w0_sig  = WORKSPACE.w0_sig_std ;
  double  wa_sig  = WORKSPACE.wa_sig_std ;
  double cov_w0wa=0., cov_w0omm =0., rho_w0wa=0., rho_w0omm=0.;
  double probsum = 0., prob=0;
  double diff_w0, diff_wa, diff_omm, sig_product ;
  bool dofit_w0wa = INPUTS.dofit_w0wa;
  
  // ----------- BEGIN --------------
  
  cpar.mushift = 0.0 ;

  for( i=0; i < INPUTS.w0_steps; i++) {     
    cpar.w0 = INPUTS.w0_min + i*INPUTS.w0_stepsize;
    for( kk=0; kk < INPUTS.wa_steps; kk++) { 
      cpar.wa = (INPUTS.wa_min + kk*INPUTS.wa_stepsize); 
      for(j=0; j < INPUTS.omm_steps; j++) {
	cpar.omm = INPUTS.omm_min + j*INPUTS.omm_stepsize;    

	prob     = WORKSPACE.extprob3d[i][kk][j];
	probsum += prob;
	diff_w0  = cpar.w0  - WORKSPACE.w0_mean;
	diff_omm = cpar.omm - WORKSPACE.omm_mean;
	cov_w0omm +=  (diff_w0 * diff_omm * prob) ;

	if ( dofit_w0wa ) {
	  diff_wa   = cpar.wa - WORKSPACE.wa_mean;
	  cov_w0wa +=  (diff_w0 * diff_wa * prob) ;
	}

      } // end j
    } // end kk 
  } // end i

  // - - - - - -
  // Compute cov and reduced covariances
  sig_product = w0_sig * omm_sig ;  
  cov_w0omm /= probsum;
  if ( sig_product > 0.0 ) 
    { rho_w0omm = cov_w0omm / sig_product; }
  else {
    rho_w0omm = 0.0; 
    WORKSPACE.NWARN++ ; 
    printf(" WARNING: cannot compute rho_w0omm because w0_sig=omm_sig=0\n");
    fflush(stdout);
  }

  if( dofit_w0wa ) { 
    sig_product = w0_sig * wa_sig ;
    cov_w0wa  /= probsum ;
    if ( sig_product > 0.0 ) 
      { rho_w0wa = cov_w0wa / sig_product ; }
    else {
      rho_w0wa = 0.0 ; 
      WORKSPACE.NWARN++ ;
      printf(" WARNING: cannot compute rho_w0wa because w0_sig=wa_sig=0\n");
      fflush(stdout);
    }
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
  //
  // Apr 26 2024: set [var]_sig_final = average of lower & upper
  
  int    Ndof = WORKSPACE.Ndof ;
  Cosparam cpar;
  double sigmu_int=0.0, sigmu_int1=0.0, sigmu_tmp, sqmusig_tmp ;
  double dif, mindif, muoff_tmp, snchi_tmp, chi2_tmp ;
  double w0_final, wa_final, omm_final ;
  double w0_sig_final, wa_sig_final, omm_sig_final ;
  double w0_ran=0.0,  wa_ran=0.0, omm_ran=0.0 ;
  double chi2_final, sigint_binsize ;
  int i;
  char fnam[] = "wfit_final" ;

  // ----------- BEGIN -----------

  printf("# =================================== \n");
  printf(" Extract final cosmology quantities \n");
  fflush(stdout);

  // load final cosmology parameters based on marg or minimize option
  if ( INPUTS.use_marg ) {
    omm_final = WORKSPACE.omm_mean ;
    w0_final  = WORKSPACE.w0_mean ;  
    wa_final  = WORKSPACE.wa_mean;

    // Apr 2024: use the real 68% prob

    if ( INPUTS.use_sig_68 ) {
      omm_sig_final = 0.5*(WORKSPACE.omm_sig_lower + WORKSPACE.omm_sig_upper);
      w0_sig_final  = 0.5*(WORKSPACE.w0_sig_lower  + WORKSPACE.w0_sig_upper);
      wa_sig_final  = 0.5*(WORKSPACE.wa_sig_lower  + WORKSPACE.wa_sig_upper);
    }
    else if ( INPUTS.use_sig_std ) {
      omm_sig_final = WORKSPACE.omm_sig_std;
      w0_sig_final  = WORKSPACE.w0_sig_std;
      wa_sig_final  = WORKSPACE.wa_sig_std;    
    }
     
  }
  else {
    w0_final  = WORKSPACE.w0_atchimin ; 
    wa_final  = WORKSPACE.wa_atchimin ; 
    omm_final = WORKSPACE.omm_atchimin ;

    // use STD approximation
    omm_sig_final = WORKSPACE.omm_sig_std;
    w0_sig_final  = WORKSPACE.w0_sig_std;
    wa_sig_final  = WORKSPACE.wa_sig_std;    
  }

  // store final results before adding random offsets
  cpar.w0  = w0_final;
  cpar.wa  = wa_final;
  cpar.omm = omm_final;
  cpar.ome = 1. - cpar.omm;
  cpar.mushift = 0.0 ;

  // add unknown offset if blind option 
  if ( INPUTS.blind ) {
    int dofit_wcdm = INPUTS.dofit_wcdm ;
    int dofit_lcdm = INPUTS.dofit_lcdm;
    int dofit_w0wa = INPUTS.dofit_w0wa ;
    srand(INPUTS.blind_seed);
    if (dofit_wcdm || dofit_w0wa)
      { w0_ran     = floor(1e6*rand()/RAND_MAX);}
    if (dofit_w0wa)
      { wa_ran     = floor(1e6*rand()/RAND_MAX);}
    omm_ran    = floor(1e6*rand()/RAND_MAX);
    w0_final   += sin(w0_ran);
    wa_final   += sin(wa_ran);
    omm_final  += sin(omm_ran);
  }

  // store final results and random numbers
  WORKSPACE.w0_final  = w0_final ;
  WORKSPACE.wa_final  = wa_final ; 
  WORKSPACE.omm_final = omm_final ; 

  WORKSPACE.w0_sig_final  = w0_sig_final ;
  WORKSPACE.wa_sig_final  = wa_sig_final ; 
  WORKSPACE.omm_sig_final = omm_sig_final ; 
  
  WORKSPACE.w0_ran  = w0_ran ;
  WORKSPACE.wa_ran  = wa_ran ;
  WORKSPACE.omm_ran = omm_ran ;

  // get final chi2 and mu_obs_list from final parameters.
  malloc_HDarrays(+1, HD_LIST[0].NSN, &HD_FINAL);
  HD_FINAL.NSN = HD_LIST[0].NSN;
  get_chi2_fit ( cpar.w0, cpar.wa, cpar.omm, INPUTS.sqsnrms,  // inputs
		 HD_FINAL.z, HD_FINAL.mu, HD_FINAL.f_interp,
		 &HD_FINAL.muresid_avg, &snchi_tmp, &chi2_final );   // return args

  // correctd HD_FINAL.mu so that wgt avg resid is zero (Dec 2024)
  // .xyz
  
  WORKSPACE.chi2_final  = chi2_final;

    // summarize chi2 and chi2 for each prior (Apr 2024)
  double chi2_om, chi2_cmb, chi2_bao, chi2_rd, rd ;
  get_chi2_priors(&cpar, &chi2_om, &chi2_cmb, &chi2_bao, &chi2_rd);
  printf("\t chi2tot = %8.2f \n", chi2_final);
  printf("\t chi2_prior(OMM) = %8.2f \n" , chi2_om);
  printf("\t chi2_prior(CMB) = %8.2f \n" , chi2_cmb);

  if ( INPUTS.use_bao > 0 ) {
    rd = rd_bao_prior(OPT_RD_CALC, &cpar);
    printf("\t chi2_prior(BAO)  = %8.2f   (rd=%.2f)\n" , chi2_bao, rd);
    printf("\t chi2_prior(rd)   = %8.2f \n" , chi2_rd);
  }
  
  fflush(stdout);

  // - - - -  sigmu_int - - - - 
  WORKSPACE.sigmu_int = 0.0;
  if ( INPUTS.fitnumber == 1 ) { return; } // Nov 24 2021
  
  return;
} // end wfit_final

// ==================================
void wfit_FoM(void) {
  
  // Estimate FoM for w0wa model or wCDM model
  // Apr 2024: use  WORKSPACE.[var]_sig_final
  // May 2024: compute FoM for either w0m or w0wa
  
  int Ndof = WORKSPACE.Ndof;
  double extchi_min = WORKSPACE.extchi_min;
  int i, kk, j;
  double extchi, extchi_dif, chi_approx, snchi_tmp, extchi_tmp, muoff_tmp;;
  double sig_product, rho ;
  Cosparam cpar;
  char string_vars[20];
  char fnam[] = "wfit_FoM" ;
  // --------------BEGIN --------------
  
  WORKSPACE.FoM_final = 0.0 ;

  if ( INPUTS.dofit_w0wa ) {
    sig_product = (WORKSPACE.w0_sig_final * WORKSPACE.wa_sig_final);
    rho         = WORKSPACE.rho_w0wa;
    sprintf(string_vars,"w0wa");
  }
  else {
    // wom for wCDM model
    sig_product = (WORKSPACE.w0_sig_final * WORKSPACE.omm_sig_final);
    rho = WORKSPACE.rho_w0omm;
    sprintf(string_vars,"womm");
  }


  if( fabs(rho) > 1. ){
    sprintf(c1err,"Invalid rho(%s) = %f \n",  string_vars, rho);
    sprintf(c2err,"Check Covariance Calculation ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  sig_product *= sqrt(1.0 - rho*rho);

  
  if ( sig_product > 0. ) 
    { WORKSPACE.FoM_final = 1.0/sig_product; }
  else 
    { WORKSPACE.FoM_final = -9.0; }

  printf("# ====================================== \n");
  printf(" FOM(%s) = %.2f\n", string_vars, WORKSPACE.FoM_final);
  fflush(stdout);

  return ;
} // end wfit_FoM

  

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
void invert_mucovar(COVMAT_DEF *MUCOV, double sqmurms_add) {

  // Mar 2023 - refactor to pass COVMAT struct
  //
  int  NSN    = MUCOV->NDIM ;
  int  i;
  time_t t0, t1;
  bool check_inverse = (INPUTS.debug_flag == 1000);
  double *MUCOV_ORIG ;
  int LDMP_MUCOV = INPUTS.ndump_mucov > 0 ;
  char fnam[] = "invert_mucovar" ;
  // ---------------- BEGIN --------------

  if ( INPUTS.use_mucov == 2 ) {
    printf("\t COV already inverted, so nothing to invert here.\n");
    fflush(stdout);
    return;
  }

  // - - - - - -
  printf("\t Invert %d x %d mucov matrix with COV_DIAG += %f \n", 
	 NSN, NSN, sqmurms_add);
  fflush(stdout);
    
  if ( sqmurms_add != 0.0 ) 
    { printf("\t\t xxx WARNING: fix bug and add sqmurms ... \n"); }

  if ( LDMP_MUCOV ) { dump_MUCOV(MUCOV,"MUCOV"); }
    
  fflush(stdout);

  t0 = time(NULL);

  if ( check_inverse ) {
    int MEMD = NSN * NSN * sizeof(double);
    MUCOV_ORIG = (double*) malloc(MEMD);
    for(i=0; i < NSN*NSN; i++ ) { MUCOV_ORIG[i] = MUCOV->ARRAY1D[i]; }
  }
  
  invertMatrix( NSN, NSN, MUCOV->ARRAY1D ) ;

  print_elapsed_time(t0, "invert matrix", UNIT_TIME_SECOND);

  if ( check_inverse ) {
    check_invertMatrix(NSN,MUCOV_ORIG,MUCOV->ARRAY1D);
    free(MUCOV_ORIG);
  }

  char *ptr_outFile = INPUTS.outFile_mucovtot_inv;
  if ( strlen(ptr_outFile) > 0 ) {

    printf("\n %s: write covtot_inv to %s\n", fnam, ptr_outFile);
    FILE *fpout = fopen(ptr_outFile,"wt");
    fprintf(fpout,"%d\n", NSN);
    for(i=0; i < NSN*NSN; i++ )
      { fprintf(fpout,"%14.8le\n", MUCOV->ARRAY1D[i]); }
    fclose(fpout);
    debugexit(fnam);
  }
 
  if ( LDMP_MUCOV ) { dump_MUCOV(MUCOV,"MUCOV^-1"); }
  
  fflush(stdout);

  return ;

} // end of invert_mucovar


// =========================================
void check_invertMatrix(int N, double *COV, double *COVINV ) {
  // compute C*C_inverse and compute maximum deviations from 1 on diagonal 
  // and max deviations from 0 on off-diagonal
  // print deviations to stdout

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
void get_chi2_fit ( 
	      double w0             // (I) 
	      ,double wa           // (I)
	      ,double OM           // (I) 
	      ,double sqmurms_add  // (I) anomalous mu-error squared
	      ,double *z_obs_list  // (O) redshifts used in fit
	      ,double *mu_obs_list // (O) mu per event (for HDIBC)
	      ,double *f_interp_list // (O) interp frac per event (for HDIBC)
	      ,double *muresid_avg   // (O) wgt-avg mu-residual (B/C)
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
  int  NSN       = HD_LIST[0].NSN;
  int  Ndof      = WORKSPACE.Ndof ;
  double sig_chi2min_naive = WORKSPACE.sig_chi2min_naive;
  double nsig_chi2min_skip = WORKSPACE.nsig_chi2min_skip;  
  double chi_hat_naive     = (double)Ndof;

  double OE, rz, sqmusig, sqmusiginv, Bsum, Csum ;
  double nsig_chi2, chi_hat, chi_tmp ;
  double dmu, dmu0, dmu1, mu_cos, mu_obs  ;
    
  double  chi2_prior = 0.0 ;
  double *rz_list  = (double*) malloc(NSN * sizeof(double) );
  double *dmu_list = (double*) malloc(NSN * sizeof(double) );
  Cosparam cparloc;
  int k, k0, k1, N0, N1, k1min, n_count=0, LDMP=0 ;
  
  HD_DEF *HD0 = &HD_LIST[0];
  HD_DEF *HD1 = &HD_LIST[1];

  double mushift0 = HD0->cospar_biasCor.mushift ;
  double mushift1 = HD1->cospar_biasCor.mushift ;

    
  // rz-interp variables
  int n_logz, iz;
  double z ;
  double mu_obs0, mu_obs1, f_interp;

  char fnam[] = "get_chi2_fit";

  // --------- BEGIN --------

  OE = 1.0 - OM ;
  cparloc.omm = OM ;
  cparloc.ome = OE ;
  cparloc.w0  = w0 ;
  cparloc.wa  = wa ;
  cparloc.mushift = 0.0 ;

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

    z = HD0->z[k] ;

    if ( USE_SPEED_INTERP )  { 
      exec_rz_interp(k, &cparloc, &rz, &mu_cos); 
    }
    else { 
      // brute force calculation of theory distance
      rz     = codist(z, &cparloc);
      mu_cos = get_mu_cos(z, rz) ;
    }

    rz_list[k]  = rz ;

    if ( INPUTS.USE_HDIBC ) {

      // interpolate two HDs
      double mu_bcor0 = HD0->mu_cospar_biascor[k] + mushift0 ;
      double mu_bcor1 = HD1->mu_cospar_biascor[k] + mushift1 ;

      // LDMP = ( fabs(OM-0.3) < 0.1 && fabs(w0+1.0)<0.5 && fabs(wa)<0.5 );
      if ( LDMP ) {
	printf(" xxx --------------------------------------------- \n");
	printf(" xxx %s: z=%.3f  mu_cos=%.3f   mu_bcor = %.3f %.3f \n",
	       fnam, z, mu_cos, mu_bcor0, mu_bcor1 ); 
	printf(" xxx %s: cospar(fit) OM,OE,w0,wa = %.3f %.3f %.3f %.3f \n",
	       fnam, OM, OE, w0, wa );
	fflush(stdout);
	debugexit(fnam);
      }

      mu_obs0 = HD0->mu[k] ;
      mu_obs1 = HD1->mu[k] ;
      f_interp = (mu_cos - mu_bcor0) / (mu_bcor1 - mu_bcor0) ;

      // avoid too much extrapolation
      if ( f_interp < -0.5 ) { f_interp = -0.5; }
      if ( f_interp >  1.5 ) { f_interp =  1.5; }

      mu_obs = mu_obs0 + f_interp * ( mu_obs1 - mu_obs0 );
    }
    else {
      // conventional : just one HD from one biasCor sim
      mu_obs = HD0->mu[k];
      f_interp = 0.0 ;
    }

    z_obs_list[k]    = z;      // Oct 23 2024 RK    
    mu_obs_list[k]   = mu_obs; // Oct 23 2024 RK
    f_interp_list[k] = f_interp;
    
    dmu_list[k] = mu_obs - mu_cos; 

    n_count++ ;
    if ( use_mucov ) {
      sqmusiginv = WORKSPACE.MUCOV_FINAL.ARRAY1D[k*(NSN+1)]; 
    }
    else {
      sqmusig     = HD0->mu_sqsig[k] + sqmurms_add ; // xxx INTERP?? xx
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
	sqmusiginv = WORKSPACE.MUCOV_FINAL.ARRAY1D[k]; // Inverse of the matrix 

	dmu0     = dmu_list[k0];
	dmu1     = dmu_list[k1];

	chi_tmp  = (sqmusiginv * dmu0 * dmu1 );

	Bsum    += sqmusiginv*(dmu0+dmu1); // Eq. A.11 of Goliath 2001  
	Csum    += (2.0*sqmusiginv);       // Eq. A.12 of Goliath 2001
	chi_hat += (2.0*chi_tmp);
	n_count++ ;

      } // end k1
    } // end k0
  }  // end do_offdiag

  //  - - - - - -

  /* Analytic marginalization over H0.  
     See appendices in Goliath et al, A&A, 2001 ,
     and add contribution from PRIOR(H0) = Gaussian instead of flat.
  */
  
  Csum    += 1./SQSIG_MUOFF ;  // H0 prior term
  *chi2sn  =  chi_hat - Bsum*Bsum/Csum ;


  *muresid_avg = Bsum/Csum; // Dec 7 2024
  
  // Ignore constant term:
  //  logCsum  = log(Csum/(2.*M_PI));
  //  *chi2sn += logCsum;

  *chi2tot = *chi2sn ;     // load intermediate function output

  double chi2_om, chi2_cmb, chi2_bao, chi2_rd;
  get_chi2_priors(&cparloc, &chi2_om, &chi2_cmb, &chi2_bao, &chi2_rd);
  
  *chi2tot += (chi2_om + chi2_cmb + chi2_bao + chi2_rd);

  free(rz_list);
  free(dmu_list);

  return ;

}  // end of get_chi2_fit

// =======================
void get_chi2_priors(Cosparam *cpar, double *chi2_om, double *chi2_cmb,
		     double *chi2_bao, double *chi2_rd) {

  // Creaetd Apr 2024
  // Compute each prior term separately to enable printing
  // each prior chi2 at end of job

  double rd;
  char fnam[] = "get_chi2_priors";

  // ---------- BEGIN ---------

  *chi2_om = *chi2_cmb = *chi2_bao = *chi2_rd = 0.0 ;
  
    /* Compute likelihood with external prior: Omega_m or BAO */
  if ( INPUTS.use_bao ) {
    *chi2_bao = chi2_bao_prior(cpar, &rd);

    // tack on optional rd prior
    double rd_pull = (rd - INPUTS.rd_prior) / INPUTS.rd_prior_sig ;
    *chi2_rd = rd_pull * rd_pull ;

  } 
  else {
    // Gaussian Omega_m prior
    double nsig;
    double OM = cpar->omm ;
    nsig = (OM-INPUTS.omm_prior)/INPUTS.omm_prior_sig ;
    *chi2_om = nsig*nsig;
  }

  // CMB prior
  if ( INPUTS.use_cmb ) {
    *chi2_cmb = chi2_cmb_prior(cpar);
  }

  return;
  
} // end get_chi2_priors

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
double chi2_bao_prior(Cosparam *cpar, double *rd_prior) {

  // Created Oct 2021
  // Refactored Apr 2024 to process more general input from $SNDATA_ROOT/models/BAO.
  // Return chi2 for bao prior.

  double chi2 = 0.0, covinv, chi2_tmp, rd_inverse ;
  int    NDATA = BAO_PRIOR.NDATA;
  int    i, i2, kk, IVAR;
  double rd, z, mean_meas, mean_model;
  double dif_vector[MXDATA_BAO_PRIOR], mean_model_vector[MXDATA_BAO_PRIOR];
  char   fnam[] = "chi2_bao_prior" ;

  // ------------ BEGIN -------------

  rd = rd_bao_prior(OPT_RD_CALC, cpar);
  *rd_prior = rd;  // return arg
  rd_inverse = 1.0/rd; // multipilcation below is faster than division?
  
  for(i=0; i < NDATA; i++ ) {
    z         = BAO_PRIOR.REDSHIFT[i];
    mean_meas = BAO_PRIOR.MEAN[i];
    IVAR      = BAO_PRIOR.IVAR[i];

    if ( IVAR == IVAR_BAO_DM_over_rd )
      { mean_model = DM_bao_prior(z, cpar) * rd_inverse ; }
    else if ( IVAR == IVAR_BAO_DH_over_rd )
      { mean_model = DH_bao_prior(z, cpar) * rd_inverse ; }
    else if ( IVAR == IVAR_BAO_DV_over_rd )
      { mean_model = DV_bao_prior(z, cpar) * rd_inverse ; }
    else {
      sprintf(c1err,"invalid IVAR=%d for %s", IVAR, BAO_PRIOR.VARNAME[i]);
      sprintf(c2err,"Valid IVAR = %d %d %d",
	      IVAR_BAO_DM_over_rd, IVAR_BAO_DH_over_rd, IVAR_BAO_DV_over_rd);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);          
    }

    mean_model_vector[i]  = mean_model ;
    dif_vector[i]         = mean_meas - mean_model ;

  } // end i loop over BAO data points

  // - - - - - - -
  // compute chi2, including full cov matrix       
  for(i=0; i < NDATA; i++ ) {
    for(i2=i; i2 < NDATA; i2++ ) {
      
      kk       = i*NDATA + i2;
      covinv   = BAO_PRIOR.COVINV1D[kk];
      
      chi2_tmp = dif_vector[i] * dif_vector[i2] * covinv;
      if ( i2 != i ) { chi2_tmp *= 2.0 ; }
      chi2 += chi2_tmp ;      
    }
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
  
  double H0     = H0_SALT2 ;
  double ld_cos, mu_cos ;
  
  ld_cos      = (1.0 + z) *  rz * c_light / H0;
  mu_cos      = 5.0 * log10(ld_cos) + 25. ;
  return mu_cos ;

}  // end get_mu_cos


// ==============================================
double get_minimized_cospar( double *w0_atchimin, double *wa_atchimin, 
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
    ,snchi_tmp, extchi_tmp, snchi_min, extchi_min, mures_tmp
    ;

  
  double w0step_tmp  = INPUTS.w0_stepsize / 10.0; 
  double wastep_tmp  = INPUTS.wa_stepsize / 10.0;
  double omstep_tmp  = INPUTS.omm_stepsize/ 10.0; 
  int    nbw0, nbwa, nbm,i, j, kk, imin=0, kmin =0, jmin=0 ;
  double nb_factor = 3.0;

  int  LDMP = 0;
  char fnam[] = "get_minimized_cospar";

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
  printf("   Minimize with refined OMstep=%6.4f (%d bins, omm=%6.4f) \n", 
	 omstep_tmp, 2*nbm, *omm_atchimin );
  printf("---------------------------------------\n");

  fflush(stdout);
  
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

	get_chi2_fit(w0_tmp, wa_tmp, om_tmp, INPUTS.sqsnrms, 
		     temp0_list, temp1_list, temp2_list,
		     &mures_tmp, &snchi_tmp, &extchi_tmp );

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

} // end of get_minimized_cospar

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
    {      rv = 1;    }
  else if (  *((double *) a)  >  *((double *) b)   )
    {      rv = -1;    }
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

double cs_over_E_integral(double amin, double amax, Cosparam *cptr) {
  // Created May 2 2024 by R.Kessler
  // return integral of c_s/E.
  return simpint(cs_over_EofA, amin, amax, cptr);
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

double cs_over_EofA(double a, Cosparam *cptr) {

  // Created May 2 2024 by R.Kessler
  //  [not yet used?]
  double Einv, Jacobian, c_over_E;
  double z = 1.0/a - 1.0 ;
  double c_s;

  z        = 1.0/a - 1.0;
  c_s      = c_sound(2, z, H0_Planck);
  
  Einv     = 1./EofA(a,cptr);
  Jacobian = 1.0/(a*a);  // dz = da/a^2, Jac=1/a^2
  c_over_E = c_s * Einv * Jacobian ;
  
  return c_over_E ;
}

double EofZ(double z, Cosparam *cptr){

  // Oct 31 2021 R.Kessler
  //  Implement omr (radiation) term such that setting omr=0.0 reduce
  //  to previous code.
  //  omr includes photons+neutrinos for early-universe integrals.
  //  This is inaccurate at late times when neutrinos become non-relativistic,
  //  but the late time omr contribution is so small that the approximation
  //  is insignificant.

  double E;
  double arg, arg_omm, arg_ome, arg_omk, arg_r, omm, omk, ome, w0, wa;
  double z1       = 1.0 + z;
  double z1_pow2  = z1*z1;
  double z1_pow3  = z1_pow2 * z1;
  double z1_pow4;
  double omr      = 0.9E-4 ; // photon+neutrino, Oct 31 2021, RK

  // ----------- BEGIN ----------
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
    arg_r     = omr * z1_pow4;
  }

  arg_omm = omm * z1_pow3 ;
  arg_omk = omk * z1_pow2 ;
  arg_ome = ome * pow(z1, 3.*(1+w0+wa))*exp(-3.0*wa*(z/z1)) ;

  arg = arg_omm + arg_omk + arg_ome ;
  if ( omr > 0.0 )  { arg += arg_r; }

  /* xxx
  if ( arg < 0.0 ) 
    { printf(" xxx z=%.3f, arg_omm=%.3f  arg_ome=%.3f  w0=%.3f\n", 
	     z, arg_omm, arg_ome, w0) ; } // xxx REMOVE
  xxx */

  // protect arg in case of om < 0 (Sep 13 2022)
  if ( arg < 0.01 ) { 
    arg = 0.01; 
    // printf(" xxx arg=%f\n", arg);
  }

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
  double ost=0.0, os=0.0, s_list[50];
  bool   is_converge, is_zero;

  for(j=1; j<=JMAX; j++){
    st = trapezoid(func,x1,x2,j,ost,vptr);
    s = (4.0*st-ost)/3.0;
    s_list[j] = s;
    if (j > 5) { 		/* Avoid spurious early convergence */     
      is_converge = (fabs(s-os) < EPS_CONVERGE_POSomm*fabs(os) ) ;
      is_zero     = (s==0.0 && os==0.0);
      if ( is_converge || is_zero )  { return s; }
    }
    os = s;
    ost = st;
  }

  // If we got to here, the integral did not converge in JMAX iterations 
  // e.g.,  simpint(one_over_EofZ, zero, z, cptr);
  printf("ERROR: simpint did not converge\n");
  double s_frac = fabs((s_list[j-1] - s_list[j-2])/s_list[j-1]);
  printf("\t [range=%.3f,%.3f om=%.3f w=%.3f s_frac=%le] \n",
	 x1, x2, vptr->omm, vptr->w0, s_frac);
  fflush(stdout);
  
  return s;
}

double trapint(double (*func)(double, Cosparam *), double x1, 
	       double x2, Cosparam *vptr)
/* based on NR 'qtrap' */
{
  int j;
  double s,olds=0.0;
  bool   is_converge, is_zero;

  for(j=1; j<=JMAX; j++){
    s = trapezoid(func,x1,x2,j,olds,vptr);
    if (j > 5) { 		/* Avoid spurious early convergence */
      is_converge = fabs(s-olds) < EPS_CONVERGE_POSomm*fabs(olds) ;
      is_zero     = (s==0.0 && olds==0.0);
      if ( is_converge || is_zero )  { return s; }
    }
    olds = s;
  }

  /* If we got to here, the integral did not converge in JMAX iterations */
  printf("ERROR: trapint did not converge\n");
  return 0.0;
}


double trapezoid(double (*func)(double, Cosparam *), double x1, double x2, 
		 int n, double s, Cosparam *vptr)
// Based on NR 'trapzd'.  Here, s must be fed back in from previous 
//   iterations, rather than using a static variable as in the NR version.
{
  double x,tnm,sum,del;
  int it,j;

  if (n == 1) {
    return 0.5*(x2-x1)*( (*func)(x1,vptr) + (*func)(x2,vptr) );
  } 
  else {
    for(it=1,j=1; j<n-1; j++) { it <<= 1; }
    tnm = it; 
    del = (x2-x1)/tnm;
    x = x1+0.5*del;
    for (sum=0.0,j=1; j<=it; j++,x+=del) 
      { sum += (*func)(x,vptr); }

    s = 0.5*(s+(x2-x1)*sum/tnm);
    return s;
  }
}


void test_c_sound(void) {

  // Created May 2024
  Cosparam cpar;
  double z, cs1, cs2, H0 = H0_Planck; 
  int iz, nz_list = 7;
  double z_list[7] = {
    1060.0, 1200.0, 1500.0, 2000.0, 5000.0, 2.0E4, 5.0E4
  } ;

  char fnam[] = "test_c_sound";

  // --------- BEGIN -----------
  
  cpar.omm  =  0.31;
  cpar.ome  =  0.69;
  cpar.w0   = -1.0;
  cpar.wa   =  0.0;
  cpar.mushift = 0.0 ;

  for(iz=0; iz < nz_list; iz++ ) {
    z = z_list[iz];
    cs1 = c_sound(1, z, H0);
    cs2 = c_sound(2, z, H0);

    printf(" xxx %s: z=%8.1f  c_s(1,2) = %.3f  %.3f \n",
	   fnam, z, cs1, cs2);
    
  }
  
  debugexit(fnam);
  return ;
} // end test_c_sound

void test_codist(void) {

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
  cpar.mushift = 0.0 ;

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
   return ;
} // end of getname

// =================================
int cidindex(char *cid) {

  int i;
  for ( i=0; i < HD_LIST[0].NSN; i++ ) {
    if ( strcmp(HD_LIST[0].cid[i],cid) == 0 ) return i;
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

  // Print residuals
  write_output_resid();
 
  // write chi2 & prob maps to fits files;
  write_output_chi2grid();

  return;
} // end WRITE_OUTPUT_DRIVER

// ********************************
void write_output_cospar(void) {

  // Created Aug 15 2020
  // Refactor Aug 2021 to avoid too much redundant code.
  //
  // format_cospar = 1 : legacy csv format
  // format_cospar = 2 : YAML format
  //
  // Sep 27 2022: write rho_wom

  int  use_marg  = INPUTS.use_marg;
  int  dofit_w0wa= INPUTS.dofit_w0wa ;
  int  format    = INPUTS.format_cospar;
  bool blind     = INPUTS.blind ;
  char *outFile  = INPUTS.outFile_cospar;
  double dt_fit  = (double)(t_end_fit  - t_start)/ 60.0 ;

#define MXVAR_WRITE 20
  int    ivar, NVAR = 0;
  FILE *fp ;
  char   VALUES_LIST[MXVAR_WRITE][40];
  char   VARNAMES_LIST[MXVAR_WRITE][20], LINE_STRING[200] ;
  char   ckey[40], cval[40], vv[20];
  char sep[] = " " ;
  char fnam[] = "write_output_cospar" ;

  // ----------- BEGIN -------------

  if ( strlen(outFile) == 0){
    sprintf(outFile, "%s.cospar", INPUTS.HD_infile_list[0]);
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
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.w0_sig_final ) ;
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
      sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.wa_sig_final ) ;
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
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.omm_sig_final ) ;
    NVAR++ ;
  }
  else {
    // WARNING: we don't have omm_sig_upper/lower, so just
    // write omm_sig_marg
    sprintf(VARNAMES_LIST[NVAR],"%ssig_marg", varname_omm); 
    sprintf(VALUES_LIST[NVAR], "%8.5f",  WORKSPACE.omm_sig_final ) ;
    NVAR++ ;
  }

  sprintf(vv,"%s%s", varname_w, varname_omm);
  sprintf(VARNAMES_LIST[NVAR],"rho_%s", vv );
  sprintf(VALUES_LIST[NVAR], "%6.3f",  WORKSPACE.rho_w0omm ) ;
  NVAR++;

  if ( dofit_w0wa ) {
    sprintf(vv,"%s%s", varname_w, varname_wa);
    sprintf(VARNAMES_LIST[NVAR],"rho_%s", vv );
    sprintf(VALUES_LIST[NVAR], "%6.3f",  WORKSPACE.rho_w0wa ) ;
    NVAR++;
  }

  // May 2024: always write FoM, even for wCDM(w,omm)
  sprintf(VARNAMES_LIST[NVAR],"FoM" );
  sprintf(VALUES_LIST[NVAR], "%6.1f",  WORKSPACE.FoM_final ) ;
  NVAR++ ;  

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
void write_output_resid(void) {

  // Jan 31 2023 R.Kessler
  // resurrect to write MU-resids and chi2 per SN
  // Oct 2024: Use HD_FINAL struct to write info.
  
  Cosparam cpar;
  FILE *fpresid;
  char TEMP_STRING[MXPATHLEN], HDIBC_STRING[200] ;
  char *outFile = INPUTS.outFile_resid ;
  char fnam[] = "write_output_resid" ;

  // ----------- BEGIN -------------

  if ( IGNOREFILE(outFile) ) { return; }

  printf("  Write chi2 and MU-residuals to %s  \n", outFile);
  fpresid = fopen(outFile, "wt");
          
  cpar.omm = WORKSPACE.omm_final ;
  cpar.ome = 1.0 - cpar.omm ;
  cpar.w0  = WORKSPACE.w0_final ;
  cpar.wa  = WORKSPACE.wa_final ;
  cpar.mushift = 0.0 ;

  int i;
  char   *cid ;
  double z, rz, ld_cos, mu_obs, mu_model, mu_sim, mu_dif, mu_sig;
  double chi2;
  double H0 = H0_SALT2;

  fprintf(fpresid,"# mu_res   = mu_obs - mu_model \n");
  fprintf(fpresid,"# mu_model = distance computed from best-fit cosmology model.\n");
  fprintf(fpresid,"# mu_sim   = distance computed from { OM, ODE, w0, wa } = "
	  "{ %.3f %.3f %.3f %.3f } \n",
	  COSPAR_SIM.omm, COSPAR_SIM.ome, COSPAR_SIM.w0, COSPAR_SIM.wa);
  fprintf(fpresid,"# chi2_diag = mu_res^2 / mu_sig^2 \n");
  fprintf(fpresid,"\n");

  sprintf(TEMP_STRING,"VARNAMES: "
	  "CID   zHD  mu_obs  mu_model  mu_sim  mu_res  mu_sig  chi2_diag");
  if ( INPUTS.USE_HDIBC ) {
    strcat(TEMP_STRING,"  mu_obs0_hdibc mu_obs1_hdibc mu_bcor0_hdibc mu_bcor1_hdibc "
	   "z_dif_hdibc f_interp_hdibc");
  }
  
  fprintf(fpresid,"%s\n", TEMP_STRING);
  
  for (i=0; i < HD_FINAL.NSN; i++) {
    cid          = HD_LIST[0].cid[i] ;    
    z            = HD_FINAL.z[i] ;  // can differ from HD_LIST[0] for HDIBC
    mu_obs       = HD_FINAL.mu[i];  // idem
    mu_obs      -= HD_FINAL.muresid_avg; // Dec 7 2024
    mu_sig       = HD_LIST[0].mu_sig[i];
    mu_sim       = HD_LIST[0].mu_sim[i];
    
    rz           =  codist( z, &cpar) ;    
    mu_model     =  get_mu_cos(z, rz) ; // best fit model mu    
    mu_dif       =  mu_obs - mu_model;  // muresid
    chi2         =  (mu_dif*mu_dif) / ( mu_sig*mu_sig );

    sprintf(TEMP_STRING,"SN: "
	    "%-12s %7.4f %7.4f %7.4f %7.4f "
	    " %7.4f %7.4f %7.1f  ", 
	    cid, z, mu_obs, mu_model, mu_sim,
	    mu_dif, mu_sig, chi2 );
    
    if ( INPUTS.USE_HDIBC )  {
      double z_dif, f_interp, mu_obs0, mu_obs1, mu_bcor0, mu_bcor1;
      double mushift0 = HD_LIST[0].cospar_biasCor.mushift;
      double mushift1 = HD_LIST[1].cospar_biasCor.mushift;      
      
      z_dif    = HD_LIST[1].z_orig[i] - HD_LIST[0].z_orig[i];
      f_interp = HD_FINAL.f_interp[i];
      mu_obs0  = HD_LIST[0].mu[i];
      mu_obs1  = HD_LIST[1].mu[i];
      mu_bcor0 = HD_LIST[0].mu_cospar_biascor[i] + mushift0 ;
      mu_bcor1 = HD_LIST[1].mu_cospar_biascor[i] + mushift1 ; 	
	
      sprintf(HDIBC_STRING,"  %7.4f %7.4f %7.4f %7.4f "
	      "%7.4f %.3f",
	      mu_obs0, mu_obs1, mu_bcor0, mu_bcor1,
	      z_dif, f_interp);
      strcat(TEMP_STRING,HDIBC_STRING);
    }
    fprintf(fpresid,"%s\n", TEMP_STRING);	    
    fflush(fpresid);

  }
  fclose(fpresid);

  return ;

} // end write_output_resid


// =======================================
void write_output_chi2grid(void) {

  // Created July 26 2022 by R.Kessler
  // Write chi2grid vs. fit params to text file.
  // E.g., for wCDM fit, the file output is
  // VARNAMES: w OM  chi2_sn chi2_tot
  //
  // where chi2_tot includes prior.
  //
  // Jan 31 2023: add blind offsets
  //
  char *outFile  = INPUTS.outFile_chi2grid;
  bool float_w0  = (INPUTS.w0_steps  > 1);
  bool float_wa  = (INPUTS.wa_steps  > 1);
  bool float_omm = (INPUTS.omm_steps > 1);

  FILE *fp;
  char VARLIST[100];
  char fnam[] = "write_output_chi2grid" ;

  // ------------ BEGIN -------------

  if ( IGNOREFILE(outFile) ) { return; }

  printf("  Write chi2 grid to %s \n", outFile);
  fp = fopen(outFile, "wt");
  if (fp == NULL){
    sprintf(c1err,"Could not open output chi2grid file:");
    sprintf(c2err,"%s", outFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }

  VARLIST[0] = 0 ;
  if ( float_w0  )  { strcat(VARLIST,"w0 "); }
  if ( float_wa  )  { strcat(VARLIST,"wa "); }
  if ( float_omm )  { strcat(VARLIST,"omm "); }
  strcat(VARLIST,"  chi2_sn chi2_tot ");

  fprintf(fp,"VARNAMES: ROW %s\n", VARLIST);

  // - - - -
  int i, kk, j, rownum=0;
  double snchi, extchi, w0, wa, omm ;
  char   line[100], cval[3][40];
  for(i=0; i < INPUTS.w0_steps; i++){
    for(kk=0; kk < INPUTS.wa_steps;kk++) {
      for(j=0; j < INPUTS.omm_steps; j++) {
	snchi =  WORKSPACE.snchi3d[i][kk][j]  - WORKSPACE.snchi_min;
	extchi = WORKSPACE.extchi3d[i][kk][j] - WORKSPACE.extchi_min;

	w0      = INPUTS.w0_min  + i  * INPUTS.w0_stepsize;
	wa      = INPUTS.wa_min  + kk * INPUTS.wa_stepsize ;
	omm     = INPUTS.omm_min + j  * INPUTS.omm_stepsize;

	if ( INPUTS.blind ) {
	  w0  += WORKSPACE.w0_ran ;
	  wa  += WORKSPACE.wa_ran ;
	  omm += WORKSPACE.omm_ran ;
	}

	cval[0][0] = cval[1][0] = cval[2][0] = 0;
	if ( float_w0  )  { sprintf(cval[0],"%.4f", w0); }
	if ( float_wa  )  { sprintf(cval[1],"%.4f", wa); }
	if ( float_omm )  { sprintf(cval[2],"%.4f", omm); }

	rownum++ ;
	sprintf(line,"ROW: %4d   "
		"%s %s %s   %12.5le %12.5le", 
		rownum, 
		cval[0], cval[1], cval[2],
		snchi, extchi );	
	fprintf(fp,"%s\n", line);
	if ( rownum % 10 == 0 ) { fflush(fp); }
      }
    }
  }
 

  fclose(fp);

  return;
} // end write_output_chi2grid

// ===========================
void test_cospar(void) {

  // Compare computed wfit distances against sim distances for 
  // three cosmologies:
  // w0,wa = {-1,0}  {-0.95, +0.15}  {-0.95,-0.15}

#define NCOSPAR_TEST 3
  double w0_test[NCOSPAR_TEST] = { -1.0, -0.95, -0.95 };
  double wa_test[NCOSPAR_TEST] = {  0.0, +0.15, -0.15 };

  Cosparam cpar;
  HzFUN_INFO_DEF HzFUN;
  ANISOTROPY_INFO_DEF ANISOTROPY;
  double H0 = H0_SALT2;

  FILE *fp;
  double z, rz, mu_wfit[NCOSPAR_TEST], mu_sim[NCOSPAR_TEST];
  double vpec=0.0 ;
  int  icos, irow=0;  
  char line[200];
  char outFile[] = "cospar_comparison.txt";   
  char fnam[] = "test_cospar";

  // -------------- BEGIN -----------

  printf("\n Open %s\n", outFile);
  fp = fopen(outFile,"wt");
  
  for(icos=0; icos < NCOSPAR_TEST; icos++ ) {
    fprintf(fp,"# om,w0,wa = %.4f %.4f %.4f for MU_fit%d and MU_sim%d \n",
	    OMEGA_MATTER_DEFAULT, w0_test[icos], wa_test[icos],
	    icos, icos);
  }

  fprintf(fp,"\n");

  fprintf(fp,"VARNAMES: ROW zCMB   MU_wfit0 MU_wfit1 MU_wfit2   "
	  "MU_sim0 MU_sim1 MU_sim2 \n");

  ANISOTROPY.USE_FLAG = false;

  for(z=0.01; z<2.0; z+=0.01) {

    for(icos=0; icos < NCOSPAR_TEST; icos++ ) {

      cpar.omm = OMEGA_MATTER_DEFAULT;
      cpar.ome = 1.0 - cpar.omm ;
      cpar.w0  = w0_test[icos];
      cpar.wa  = wa_test[icos];
      cpar.mushift = 0.0 ;
      rz            =  codist( z, &cpar) ;
      mu_wfit[icos] =  get_mu_cos(z,rz) ; 

      set_HzFUN_for_wfit(H0, cpar.omm, cpar.ome, cpar.w0, cpar.wa, &HzFUN);
      mu_sim[icos] = dLmag(z,z, vpec, &HzFUN, &ANISOTROPY);
    }

    irow++ ;
    sprintf(line,"ROW: %3d %.4f   "
	    "%.5f %.5f %.5f   %.5f %.5f %.5f ",
	    irow, z, 
	    mu_wfit[0], mu_wfit[1], mu_wfit[2],
	    mu_sim[0],  mu_sim[1],  mu_sim[2] );
    fprintf(fp, "%s\n", line);

  } // end z

  fclose(fp);

  debugexit(fnam);
  return;

} // end test_cospar

// ==========================
void CPU_summary(void) {

  // ----------- BEGIN -----------

  printf("# =================================== \n");
  printf(" CPU Summary \n");
  print_cputime(t_end_init, STRING_CPUTIME_PROC_ALL, UNIT_TIME_SECOND, 0);
  return;
} // end CPU_summary
