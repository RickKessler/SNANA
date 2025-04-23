#include "gsl/gsl_linalg.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

int init_genmag_BAYESN(char *MODEL_VERSION, char *MODEL_EXTRAP, int optmask);
int init_genmag_bayesn__( char *model_version, char *model_extrap, int *optmask);
void read_BAYESN_inputs(char *filename);

void genmag_BAYESN(
		  int OPTMASK     // (I) bit-mask of options (LSB=0)
		  ,int ifilt_obs  // (I) absolute filter index
		  ,double *parList_SN   // DLMAG, THETA, AV, RV
		  ,double mwebv   // (I) Galactic extinction: E(B-V)
		  ,double z       // (I) Supernova redshift
		  ,int    Nobs         // (I) number of epochs
		  ,double *Tobs_list   // (I) list of Tobs (w.r.t peakMJD)
		  ,double *magobs_list  // (O) observer-frame model mag values
		  ,double *magerr_list  // (O) observer-frame model mag errors
        );

void genmag_bayesn__(int *OPTMASK, int *ifilt_obs, double *parlist_SN,
	       	double *mwebv, double *z, int *Nobs,
	       	double *Tobs_list, double *magobs_list,
	       	double *magerr_list);

void dump_SED_element(FILE * file, double wave, double value);

gsl_matrix *invKD_irr(int Nk, double *xk);
gsl_matrix *spline_coeffs_irr(int N, int Nk, double *x, double *xk, gsl_matrix *invKD);
gsl_vector *sample_nu(int n_lam_knots, int n_tau_knots);
gsl_matrix *sample_epsilon(int n_lam_knots, int n_tau_knots, gsl_vector * nu, gsl_matrix * L_Sigma_epsilon);

void genSCATTER_BAYESN(); // renamed by ST Sep 14 2024

char BAYESN_MODELPATH[MXPATHLEN];
int VERBOSE_BAYESN;
#define OPTMASK_BAYESN_VERBOSE 128

// ST (14 Sep 2024)
int ENABLE_SCATTER_BAYESN;
#define OPTMASK_BAYESN_EPSILON 2 // enable only the non-gray (EPSILON) component of intrinsic scatter
#define OPTMASK_BAYESN_DELTAM  4 // enable only the gray (DELTAM) component of intrinsic scatter
#define OPTMASK_BAYESN_SCATTER_ALL 6 // sum of all non-default scatter bits
#define OPTMASK_BAYESN_SCATTER_DEFAULT 1 // enable all scatter components (default behaviour)

struct {
   int    n_lam_knots;
   int    n_tau_knots;
   int    n_sig_knots;

   // variables read from BayeSN model directory
   // variables that control the range of validity of the SED
   double l_filter_cen_min;
   double l_filter_cen_max;
   
   
   // variables that define the SED 
   // scalars
   double M0;
   double SIGMA0;
   double RV;
   double TAUA;

   // vectors
   double *lam_knots;
   double *tau_knots;

   // matrices
   gsl_matrix *L_Sigma_epsilon;
   gsl_matrix *W0;
   gsl_matrix *W1;

   // running epsilon variable (yikes)
   gsl_matrix *EPSILON;
   // running DELTAM variable (double yikes)
   double DELTAM;

   // for the base SED - typically Hsiao
   SEDMODEL_FLUX_DEF S0; 

   // computed quantities
   gsl_matrix *KD_tau;
   gsl_matrix *KD_lam;
   gsl_matrix *J_lam;

} BAYESN_MODEL_INFO;
