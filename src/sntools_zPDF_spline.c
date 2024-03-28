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
		      char *cid, char *method_spline, int verbose, double *mean, double *std_dev ) {
  // created Jun 2022 R. Chen
  // Initialize spline interpolate for photo-z quantiles
  // allocate interpolation accelerator and gsl_spline object
  //
  // Jan 9, 2024: pass SNID to use for messaging.
  // Mar 25,2024: return mean and RMS

  dump_zPDF(method_spline, N_Q, percentile_list, zphot_q_list, cid);

  char fnam[] = "init_zPDF_spline";
  double sum    = 0. ;
  double sum_pdf= 0. ;
  double sum_sq = 0. ;
 

  int i,i2;
  int LDMP = 0;

  // ------ BEGIN ---------
 
  
  zPDF_spline.acc    = gsl_interp_accel_alloc();

  if (strcmp(method_spline,METHOD_SPLINE_LINEAR)==0 ){
  zPDF_spline.spline = gsl_spline_alloc(gsl_interp_linear, N_Q);     // Linear
  }
  else if (strcmp(method_spline,METHOD_SPLINE_CUBIC)==0){
    zPDF_spline.spline = gsl_spline_alloc(gsl_interp_cspline, N_Q); // cubic
  }
  else if  (strcmp(method_spline,METHOD_SPLINE_STEFFEN)==0){
#ifdef GSL_INTERP_STEFFEN
    zPDF_spline.spline = gsl_spline_alloc(gsl_interp_steffen, N_Q);    // Steffan
#else
    sprintf(c1err,"Compilation does not include %s method for GSL",
	    method_spline);
    sprintf(c2err,"Needs GSL 2.7 or higher");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
#endif
  }
  else {
    sprintf(c1err,"Invalid method_spline=%s for CID=%s", method_spline, cid );
    sprintf(c2err,"Valid methods are: %s  %s  %s", 
	    METHOD_SPLINE_LINEAR, METHOD_SPLINE_CUBIC, METHOD_SPLINE_STEFFEN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
    
  sprintf(zPDF_spline.method_spline, "%s", method_spline);
  zPDF_spline.zmin   = zphot_q_list[0];
  zPDF_spline.zmax   = zphot_q_list[N_Q-1];
  int NBIN_SPLINE    = 20; // Warning -- hack 
  zPDF_spline.dz     = (zPDF_spline.zmax - zPDF_spline.zmin)/(double)NBIN_SPLINE ;
  
  // check that percentile list covers 0 and 1.0
  double P0 = percentile_list[0];
  if (P0 > 1.0e-4) {
    print_preAbort_banner(fnam);
    dump_zPDF(method_spline, N_Q, percentile_list, zphot_q_list, cid);
    sprintf(c1err,"First percentile prob = %le for CID=%s", P0, cid);
    sprintf(c2err,"Expecting prob ~ 0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  double P1 = percentile_list[N_Q-1];
  if (P1 < 0.999) {
    print_preAbort_banner(fnam);
    dump_zPDF(method_spline, N_Q, percentile_list, zphot_q_list, cid);
    sprintf(c1err,"Last percentile prob = %le for CID=%s", P1, cid);
    sprintf(c2err,"Expecting prob ~ 1.0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  for(i = 1; i<N_Q; i++){
    bool check_1 = percentile_list[i] > percentile_list[i-1];
    bool check_2 = zphot_q_list[i]    > zphot_q_list[i-1];
    //bool force_dump = true ; // Internal debugging only ??
    if(  !(check_1 && check_2)) {
      print_preAbort_banner(fnam);
      dump_zPDF(method_spline, N_Q, percentile_list, zphot_q_list, cid);
      
      sprintf(c1err,"Quantile information is not monotonically increasing for CID=%s", cid);
      sprintf(c2err,"Check both z and percentile in datafile");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  } 

  
  
  // intialize spline
  gsl_spline_init(zPDF_spline.spline, zphot_q_list, percentile_list, N_Q);

  // - - - - 
  double zmin, zmax, dz, z, pdf, pdf_max = 0.0;
  int iz = 0; 
  double *pdf_store = (double*)malloc((NBIN_SPLINE+1)*sizeof(double));
  //double pdf_store[40];
  zmin = zphot_q_list[0] ;
  zmax = zphot_q_list[N_Q-1] ;
  dz   = zPDF_spline.dz ;
  for( z = zmin; z <= zmax; z += dz ) {
    pdf = gsl_spline_eval_deriv(zPDF_spline.spline, z, zPDF_spline.acc);
    if (pdf < 0.) {pdf = 0.0 ;} // avoid unphysical negative probability
    pdf_store[iz] =  pdf; iz++; 
    if ( pdf > pdf_max ) { pdf_max = pdf; }
    if(LDMP) {
      printf("XXX %s iz = %d, z = %le, pdf = %le \n",fnam,iz,z, pdf);
      //printf("ZPDF %le %le  \n",z, pdf); // XXX
    }
    sum     += z*pdf;
    sum_pdf += pdf;
  }


  *mean = sum/sum_pdf ;
  iz = 0; 
  for( z = zmin; z <= zmax; z += dz){
    pdf = pdf_store[iz]; iz++;
    sum_sq += (z - *mean)*(z - *mean)*pdf;
  }
  if(LDMP) {
    printf("XXX %s sum = %le, sum_pdf = %le, sum_sq = %le \n",fnam,sum,sum_pdf,sum_sq);
  }
  if(sum_sq > 0  && sum_pdf > 0){
  *std_dev  = sqrt(sum_sq/sum_pdf);
  }
  else {*std_dev = 0.;}

  if(LDMP) {
  printf("XXX %s std = %le \n",fnam,*std_dev);
  }
  //free(pdf_store);
    
  zPDF_spline.pdf_max = pdf_max;

  if ( verbose ) {
    printf("\t zPhot-quantile pdf(max) = %.2f mean = %.3f  std = %.3f for CID = %s\n", pdf_max,*mean, *std_dev,  cid);
    fflush(stdout);
  }


} // END OF init_zPDF_spline

double eval_zPDF_spline(double z) {
  // created Jun 2022 R. Chen
  // Returns PDF probability at redshift z
  // Must call init_zPDF_spline before calling this function
  //
  // Jan 9 2024 RK : pdf /= pdf_max

  double pdf,z0,z1,cdf0,cdf1;
  double dz = zPDF_spline.dz;
  char fnam[] = "eval_zPDF_spline";
  // BEGIN
  if (z < zPDF_spline.zmin || z > zPDF_spline.zmax )  {
    pdf=0.0 ;
  }
  else {
    z0 = z1 = z;
    if ( z+dz < zPDF_spline.zmax) { z1 = z+dz ; }
    if ( z-dz > zPDF_spline.zmin) { z0 = z-dz ; }
    pdf = gsl_spline_eval_deriv(zPDF_spline.spline, z, zPDF_spline.acc);

    pdf /= zPDF_spline.pdf_max ; // RK Jan 9 2024: make sure that PDF(MAX)=1.0
    //printf("Derivative PDF value is %lf \n", pdf);
  }
  return pdf ;
} // END OF eval_zPDF_spline

void init_zpdf_spline__(int *N_Q, double* percentile_list, 
			double* zphot_q_list, char *cid, char *method_spline, int *verbose, double *mean, double *rms) {
  //printf("xxx init_zpdf_spline__ method_spline = %s CID = %s \n", method_spline, cid);
  init_zPDF_spline(*N_Q, percentile_list, zphot_q_list, cid, method_spline, *verbose, mean, rms);
}
double eval_zpdf_spline__(double *z) {
  return eval_zPDF_spline(*z) ; 

}




void dump_zPDF(char *method_spline, int N_Q, double* percentile_list, double* zphot_q_list,
	       char *cid){
  int i;
  char fnam[] = "dump_zPDF";

  //printf(" xxx %s for cid = %s method_spline = %s \n", fnam, cid, method_spline);
  for(i = 0; i<N_Q; i++){
    //printf(" xxx %s: Quantile %2d : z = %3f percentile = %3f\n",
    //fnam, i, zphot_q_list[i], percentile_list[i]);
    fflush(stdout);
  }
} // end dump_zPDF





