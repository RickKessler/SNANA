// Created Jun 2022
// Tools to interpolate zPDF quantiles and return probability.
//
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 

#include "sntools.h"
#include "sntools_zPDF_spline.h"

#include <sys/types.h>
#include <sys/stat.h>

void init_zPDF_spline(int N_Q, double* percentile_list, double* zphot_q_list, 
		      char *cid) {
  // created Jun 2022 R. Chen
  // Initialize spline interpolate for photo-z quantiles
  // allocate interpolation accelerator and gsl_spline object
  //
  // Jan 9, 2024: pass SNID to use for messaging.

  char fnam[] = "init_zPDF_spline";

  int i,i2;

  // ------ BEGIN ---------

  zCDF_spline.acc    = gsl_interp_accel_alloc();
  zCDF_spline.spline = gsl_spline_alloc(gsl_interp_cspline, N_Q);
  zCDF_spline.zmin   = zphot_q_list[0];
  zCDF_spline.zmax   = zphot_q_list[N_Q-1];
  zCDF_spline.dz     = (zCDF_spline.zmax - zCDF_spline.zmin)/20.0; // Warning -- hack
  
  // check that percentile list covers 0 and 1.0
  double P0 = percentile_list[0];
  if (P0 > 1.0e-4) {
    sprintf(c1err,"First percentile prob = %le for CID=%s", P0, cid);
    sprintf(c2err,"Expecting prob ~ 0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  double P1 = percentile_list[N_Q-1];
  if (P1 < 0.999) {
    sprintf(c1err,"Last percentile prob = %le for CID=%s", P1, cid);
    sprintf(c2err,"Expecting prob ~ 1.0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  for(i = 1; i<N_Q; i++){
    bool check_1 = percentile_list[i]>percentile_list[i-1];
    bool check_2 = zphot_q_list[i]>zphot_q_list[i-1];
    if(! check_1 || ! check_2){
      print_preAbort_banner(fnam);
      dump_zPDF(N_Q, percentile_list, zphot_q_list, cid);
      sprintf(c1err,"Quantile information is not monotinically incerasing for CID=%s", cid);
      sprintf(c2err,"Check both z and percentile in datafile");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  } 

  
  
  // intialize spline
  gsl_spline_init(zCDF_spline.spline, zphot_q_list, percentile_list, N_Q);

  // - - - - 
  double zmin, zmax, dz, z, pdf, pdf_max = 0.0;
  zmin = zphot_q_list[0] ;
  zmax = zphot_q_list[N_Q-1] ;
  dz   = zCDF_spline.dz ;
  for( z = zmin; z <= zmax; z += dz ) {
    pdf = gsl_spline_eval_deriv(zCDF_spline.spline, z, zCDF_spline.acc);
    if ( pdf > pdf_max ) { pdf_max = pdf; }
  }
  zCDF_spline.pdf_max = pdf_max;
  printf(" zPhot-quantile pdf(max) = %.2f  for CID = %s\n", pdf_max, cid);
  fflush(stdout);

} // END OF init_zPDF_spline

double eval_zPDF_spline(double z) {
  // created Jun 2022 R. Chen
  // Returns PDF probability at redshift z
  // Must call init_zPDF_spline before calling this function
  //
  // Jan 9 2024 RK : pdf /= pdf_max

  double pdf,z0,z1,cdf0,cdf1;
  double dz = zCDF_spline.dz;
  char fnam[] = "eval_zPDF_spline";
  // BEGIN
  if (z < zCDF_spline.zmin || z > zCDF_spline.zmax )  {
    pdf=0.0 ;
  }
  else {
    z0 = z1 = z;
    if ( z+dz < zCDF_spline.zmax) { z1 = z+dz ; }
    if ( z-dz > zCDF_spline.zmin) { z0 = z-dz ; }
    pdf = gsl_spline_eval_deriv(zCDF_spline.spline, z, zCDF_spline.acc);

    pdf /= zCDF_spline.pdf_max ; // RK Jan 9 2024: make sure that PDF(MAX)=1.0
    //printf("Derivative PDF value is %lf \n", pdf);
  }
  return pdf ;
} // END OF eval_zPDF_spline

void init_zpdf_spline__(int *N_Q, double* percentile_list, 
			double* zphot_q_list, char *cid){
  init_zPDF_spline(*N_Q, percentile_list, zphot_q_list, cid);
}
double eval_zpdf_spline__(double *z) {
  return eval_zPDF_spline(*z) ; 

}




void dump_zPDF(int N_Q, double* percentile_list, double* zphot_q_list,
	       char *cid){
  int i;
  for(i = 0; i<N_Q; i++){
    printf("Quantile %2d : z = %3f percentile = %3f\n",
	  i,zphot_q_list[i],percentile_list[i]);
    fflush(stdout);
  }
}




