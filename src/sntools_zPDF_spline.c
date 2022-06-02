// Created Jun 2022
// Tools to interpolate zPDF quantiles and return probability.

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 

#include "sntools.h"
#include "sntools_zPDF_spline.h"

#include <sys/types.h>
#include <sys/stat.h>

void init_zPDF_spline(int N_Q, double* percentile_list, double* zphot_q_list) {
  // created Jun 2022 R. Chen
  // Initialize spline interpolate for photo-z quantiles
  // allocate interpolation accelerator and gsl_spline object
  zPDF_spline.acc = gsl_interp_accel_alloc();
  zPDF_spline.spline = gsl_spline_alloc(gsl_interp_cspline, N_Q);
  
  // intialize spline
  gsl_spline_init(zPDF_spline.spline, zphot_q_list, percentile_list, N_Q);
}

double eval_zPDF_spline(double z) {
  // created Jun 2022 R. Chen
  // Returns PDF probability at redshift z
  // Must call init_zPDF_spline before calling this function
  double val;
  val = gsl_spline_eval(zPDF_spline.spline, z, zPDF_spline.acc);
  return val;
}

void init_zpdf_spline__(int *N_Q, double* percentile_list, 
			double* zphot_q_list){
  init_zPDF_spline(*N_Q, percentile_list, zphot_q_list);
}
double eval_zpdf_spline__(double *z)
{  eval_zPDF_spline(*z); }
