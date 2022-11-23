#include "gsl/gsl_linalg.h"

// define pre-processor command to use python interface
#define USE_BAYESNxxx      

int init_genmag_BAYESN(char *version, int optmask);
int init_genmag_bayesn__( char *version, int *optmask);
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

gsl_matrix *invKD_irr(int Nk, double *xk);
gsl_matrix *spline_coeffs_irr(int N, int Nk, double *x, double *xk, gsl_matrix *invKD);

struct {
   int    n_lam_knots;
   int    n_tau_knots;
   int    n_sig_knots;

   // variables read from BayeSN model director
   // scalars
   double M0;
   double sigma0;
   double RV;
   double tauA;

   // vectors
   double *lam_knots;
   double *tau_knots;

   // matrices
   gsl_matrix *L_Sigma_epsilon;
   gsl_matrix *W0;
   gsl_matrix *W1;

   //double **L_Sigma_epsilon;
   //double **W0; 
   //double **W1; 

   // for the base SED - typically Hsiao
   SEDMODEL_FLUX_DEF S0; 

   // computed quantities
   gsl_matrix *KD_tau;
   gsl_matrix *KD_lam;
   gsl_matrix *J_lam;

   //double **KD_tau;
   //double **KD_lam;
   //double **J_lam;

} BAYESN_MODEL_INFO;
