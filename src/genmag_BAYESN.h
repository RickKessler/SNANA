
int init_genmag_BAYESN(char *version, int optmask);
int init_genmag_bayesn__( char *version, int *optmask);

void genmag_BAYESN(
		  int OPTMASK     // (I) bit-mask of options (LSB=0)
		  ,int ifilt_obs  // (I) absolute filter index
		  ,double *parList_SN   // DLMAG, THETA, AV, RV
		  ,double mwebv   // (I) Galactic extinction: E(B-V)
		  ,double z       // (I) Supernova redshift
		  ,int    Nobs         // (I) number of epochs
		  ,double *Tobs_list   // (I) list of Tobs (w.r.t peakMJD)
		  ,double *magobs_list  // (O) observed mag values
		  ,double *magerr_list  // (O) model mag errors
        );

void genmag_bayesn__(int *OPTMASK, int *ifilt_obs, double *parlist_SN,
	       	double *mwebv, double *z, int *Nobs,
	       	double *Tobs_list, double *magobs_list,
	       	double *magerr_list);

#define MAX_KNOTS_BAYESN 20 //max number of knots for BayeSN


struct {
   int    n_lam_knots;
   int    n_tau_knots;

   // variables read from BayeSN model director
   // scalars
   double M0;
   double sigma0;
   double RV;
   double tauA;

   // vectors
   double lam_knots[MAX_KNOTS_BAYESN];
   double tau_knots[MAX_KNOTS_BAYESN];

   // matrices
   double **L_Sigma_epsilon;
   double **W0; 
   double **W1; 

   // for the base SED - typically Hsiao
   SEDMODEL_FLUX_DEF S0; 

   // computed quantities
   double **KD_tau;
   double **KD_lam;
   double **J_lam;
   double **J_tau;

} BAYESN_MODEL_INFO;
