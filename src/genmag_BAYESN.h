#include "gsl/gsl_linalg.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

// init
int  init_genmag_BAYESN(
			char *MODEL_VERSION
			,char *MODEL_EXTRAP
			,int    optmask
			,char *NAMES_HOSTPAR 
			);

int  init_genmag_bayesn__(
			  char *model_version
			  ,char *model_extrap
			  ,int  *optmask
			  ,char *names_hostpar
			  );


void read_BAYESN_inputs( char * filename );

void init_HOSTPAR_BAYESN(int OPTMASK, char *NAMES_HOSTPAR) ;

// genmag
void genmag_BAYESN(
		   int      OPTMASK      // (I) bit-mask of options (LSB=0)
		   ,int      ifilt_obs    // (I) absolute filter index
		   ,double * parList_SN   // (I) DLMAG, THETA, AV, RV
		   ,double * parList_HOST // (I) see init_HOSTPAR_BAYESN (Jul 2025)
		   ,double   mwebv        // (I) Galactic extinction: E(B-V)
		   ,double   z            // (I) Supernova redshift
		   ,int      Nobs         // (I) number of epochs
		   ,double * Tobs_list    // (I) list of Tobs (w.r.t peakMJD)
		   ,double * magobs_list  // (O) observer-frame model mag values
		   ,double * magerr_list  // (O) observer-frame model mag errors
		   );

void genmag_bayesn__(
     int    * OPTMASK
    ,int    * ifilt_obs
    ,double * parlist_SN
    ,double * parlist_HOST
    ,double * mwebv
    ,double * z 
    ,int    * Nobs
    ,double * Tobs_list
    ,double * magobs_list
    ,double * magerr_list
);


void genSCATTER_BAYESN(); // renamed by ST Sep 14 2024


// ------------------------------------------
// -------------- globals -------------------
// ------------------------------------------

// splines
gsl_matrix * invKD_irr(
     int          Nk
    ,double     * xk
);
gsl_matrix * spline_coeffs_irr(
     int          N
    ,int          Nk
    ,double     * x
    ,double     * xk
    ,gsl_matrix * invKD
);
// intrinsic scatter
gsl_vector * sample_nu(
     int          n_lam_knots
    ,int          n_tau_knots
);
gsl_matrix * sample_epsilon(
     int          n_lam_knots
    ,int          n_tau_knots
    ,gsl_vector * nu
    ,gsl_matrix * L_Sigma_epsilon
);



// model path
char    BAYESN_MODELPATH[MXPATHLEN];

// optmask bits
int     VERBOSE_BAYESN;
#define OPTMASK_BAYESN_VERBOSE         128
int     ENABLE_SCATTER_BAYESN;             // ST (14 Sep 2024)
#define OPTMASK_BAYESN_EPSILON         2   // enable only the non-gray (EPSILON) scatter
#define OPTMASK_BAYESN_DELTAM          4   // enable only the gray (DELTAM) scatter
#define OPTMASK_BAYESN_SCATTER_ALL     6   // sum of all non-default scatter bits
#define OPTMASK_BAYESN_SCATTER_DEFAULT 1   // enable all scatter components (default)

#define MXHOSTPAR_BAYESN 20                // max number of host params (Jul 2025)

// model structure
struct {
   int          n_lam_knots;
   int          n_tau_knots;
   int          n_sig_knots;
   // variables read from BayeSN model directory
   // variables that control the range of validity of the SED
   double       l_filter_cen_min;
   double       l_filter_cen_max;
   // variables that define the SED 
   // scalars
   double       M0;
   double       SIGMA0;
   double       RV;
   double       TAUA;
   // vectors
   double     * lam_knots;
   double     * tau_knots;
   // matrices
   gsl_matrix * L_Sigma_epsilon;
   gsl_matrix * W0;
   gsl_matrix * W1;
   // running epsilon variable (yikes)
   gsl_matrix * EPSILON;
   // running DELTAM variable (double yikes)
   double       DELTAM;
   // pre-computed spline matrices
   gsl_matrix * KD_tau;
   gsl_matrix * KD_lam;
   gsl_matrix * J_lam;
   // for the base SED - typically Hsiao
   SEDMODEL_FLUX_DEF S0; 

  // Jul 28 2025 RK - tack on hostpar info, same as in genmag_PySEDMODEL
  int  N_HOSTPAR;
  char *NAME_ARRAY_HOSTPAR[MXHOSTPAR_BAYESN] ;
  int  IPAR_LOGMASS; 

} BAYESN_MODEL_INFO;
