
struct{
  gsl_interp_accel *acc;
  gsl_spline       *spline;
  double zmin, zmax;
} zPDF_spline ;


void init_zPDF_spline(int N_Q, double* percentile_list, double* zphot_q_list);
double eval_zPDF_spline(double z);

void init_zpdf_spline__(int *N_Q, double* percentile_list, double* zphot_q_list);
double eval_zpdf_spline__(double *z);
