// generalized spline code
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 

#include "sntools.h"
#include "sntools_zPDF_spline.h"

#include <sys/types.h>
#include <sys/stat.h>


void init_spline(int *index, int NVAL, double* x_list, double* y_list, 
                      char *cid, char *method_spline, int verbose,
                      double *mean, double *std_dev, int *error_flag){
  char fnam[] = "init_spline()";
  printf("XXX Hello from %s, INDEX = %d, N = %d\n",fnam, *index, NVAL);
}


void init_spline__(int *index, int *NVAL, double* x_list, double* y_list, 
                      char *cid, char *method_spline, int *verbose,
                      double *mean, double *std_dev, int *error_flag){
  init_spline(index, *NVAL,x_list, y_list, cid, method_spline, *verbose, mean, std_dev, error_flag);
}
