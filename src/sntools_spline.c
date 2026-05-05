// generalized spline code
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 

#include "sntools.h"
#include "sntools_zPDF_spline.h"

#include <sys/types.h>
#include <sys/stat.h>


void init_spline(int *index, int NVAL, double* x_list, double* y_list, 
		 char *cid, char *method_spline, char *name, int verbose,
                      double *mean, double *std_dev, int *error_flag){
  char fnam[] = "init_spline()";
  printf("XXX Hello from %s, INDEX = %d, N = %d, name = %s\n",fnam, *index, NVAL, name);
}


void init_spline__(int *index, int *NVAL, double* x_list, double* y_list, 
		   char *cid, char *method_spline, char *name, int *verbose,
                      double *mean, double *std_dev, int *error_flag){
  init_spline(index, *NVAL,x_list, y_list, cid, method_spline, name, *verbose, mean, std_dev, error_flag);
}

double eval_spline(double x) {
  char fnam[] = "eval_spline";
  double y = 1.0;
  return y;
}


double eval_spline__(double *x) {
  return eval_spline(*x);
}





