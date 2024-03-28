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
		      char *cid, char *method_spline, int verbose, double *mean, double *std_dev );
double eval_zPDF_spline(double z);


void dump_zPDF(char *method_spline, int N_Q, double* percentile_list, double* zphot_q_list,
                      char *cid);
// Mangled fortran functions


void init_zpdf_spline__( int *N_Q, double* percentile_list, double* zphot_q_list, 
			 char *cid, char *method_spline, int *verbose, double *mean, double *std_dev );
double eval_zpdf_spline__(double *z);

