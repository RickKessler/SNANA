
// Mar 28 2024: Define pre-proc flag for using STEFFEN interp option.
//    Adequate GSL version available only on Perlmutter;
//    not on FNAL or RCC.

#define GSL_INTERP_STEFFEN   
#define METHOD_SPLINE_LINEAR "LINEAR"
#define	METHOD_SPLINE_CUBIC "CUBIC"
#define	METHOD_SPLINE_STEFFEN "STEFFEN"


struct{
  char method_spline[40];
  gsl_interp_accel *acc;
  gsl_spline       *spline;
  double zmin, zmax, dz; // dz used for derivative calculation
  double pdf_max ;
} zPDF_spline ;


void init_zPDF_spline(int N_Q, double* percentile_list, double* zphot_q_list, 
		      char *cid, char *method_spline, int verbose, double *mean, double *std_dev, int *error_flag );
double eval_zPDF_spline(double z);

void dump_zPDF(char *method_spline, int N_Q, double* percentile_list, double* zphot_q_list,
                      char *cid);

// Mangled fortran functions for snlc_fit

void init_zpdf_spline__( int *N_Q, double* percentile_list, double* zphot_q_list, 
			 char *cid, char *method_spline, int *verbose, double *mean, double *std_dev, int *error_flag );
double eval_zpdf_spline__(double *z);

