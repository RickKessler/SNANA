/***************************
  Created Aug 2013 by R.Kessler

  
**************************/
 
double RVMW_S11DM15 ;

char S11DM15_INFO_FILE[MXPATHLEN];
char S11DM15_MODELPATH[MXPATHLEN];

SEDMODEL_FLUX_DEF  S11DM15_SEDMODEL ; // define structure

struct S11DM15_INFO {
  double MBPEAK ;
  double DM15_REF;
  double a_LAMPAR[4];
  double b_PAR ;
  double TAUPOLY_STRETCH[4];
  double RESTLAMBDA_RANGE[2];
} S11DM15_INFO ;


double Trest_Bmax_S11DM15; // Trest at max B flux

// define flux-scale lookup tables vs. wavelength:
double **S11DM15_TABLE_FLUXSCALE;     // for 10**[-0.4(a*x + b*x^2)]
double **S11DM15_TABLE_XTHOST_FRAC ;  // for host extinction 


struct S11DM15_LAST {
  double z, dm15, RV, AV;  // to speed integration
} S11DM15_LAST ;


// ====== function prototypes ============

int  init_genmag_S11DM15(char *version, int OPTMASK );
void malloc_S11DM15_SEDMODEL(void);

void read_S11DM15_INFO_FILE(void) ;

void genmag_S11DM15(int ifilt_obs, double dm15, 
		    double AV, double RV, double mwebv,
		    double z, double mu, int Nobs, double *Tobs_list,
		    double *magobs_list, double *magerr_list );

double get_Tmax_S11DM15(void);

double get_FluxScale_S11DM15(double lambda, double dm15) ;

double get_tau_S11DM15(double dm15) ;

double get_magerr_S11DM15(double Trest);

void fill_TABLE_FLUXSCALE_S11DM15(double dm15, double RV, double AV, double z);

double INTEG_SEDFLUX_S11DM15(int ifilt_obs, double z, double Tobs);
double INTEG_BFLUX_S11DM15(double Trest) ;

void  dumpStuff_S11DM15(void) ;

// ========== END =============

