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
  char fnam[] = "init_zPDF_spline";
  // BEGIN
  zCDF_spline.acc = gsl_interp_accel_alloc();
  zCDF_spline.spline = gsl_spline_alloc(gsl_interp_cspline, N_Q);
  zCDF_spline.zmin = zphot_q_list[0];
  zCDF_spline.zmax = zphot_q_list[N_Q-1];
  zCDF_spline.dz = (zCDF_spline.zmax - zCDF_spline.zmin)/20.0; // Warning -- hack
  
  // check that percentile list covers 0 and 1.0
  double P0 = percentile_list[0];
  if (P0 > 1.0e-4) {
    sprintf(c1err,"First percentile prob =%f", P0);
    sprintf(c2err,"Expecting prob ~ 0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  double P1 = percentile_list[N_Q-1];
  if (P1 < 0.999) {
    sprintf(c1err,"Last percentile prob =%f", P1);
    sprintf(c2err,"Expecting prob ~ 1.0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  
  // intialize spline
  gsl_spline_init(zCDF_spline.spline, zphot_q_list, percentile_list, N_Q);
} // END OF init_zPDF_spline

double eval_zPDF_spline(double z) {
  // created Jun 2022 R. Chen
  // Returns PDF probability at redshift z
  // Must call init_zPDF_spline before calling this function
  double pdf,z0,z1,cdf0,cdf1;
  double dz = zCDF_spline.dz;
  char fnam[] = "eval_zPDF_spline";
  // BEGIN
  if (z < zCDF_spline.zmin) 
  	{pdf=0.;}
  else if (z > zCDF_spline.zmax)
        {pdf=0.;}
  else {
        z0=z1=z;
        if (z+dz < zCDF_spline.zmax) {z1 = z+dz;}
        if (z-dz > zCDF_spline.zmin) {z0 = z-dz;}
  	cdf0 = gsl_spline_eval(zCDF_spline.spline, z0, zCDF_spline.acc);
        cdf1 = gsl_spline_eval(zCDF_spline.spline, z1, zCDF_spline.acc);
        pdf = (cdf1-cdf0)/(z1-z0);
  }
  return pdf;
} // END OF eval_zPDF_spline

void init_zpdf_spline__(int *N_Q, double* percentile_list, 
			double* zphot_q_list){
  init_zPDF_spline(*N_Q, percentile_list, zphot_q_list);
}
double eval_zpdf_spline__(double *z)
{  eval_zPDF_spline(*z); }
