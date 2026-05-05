//generalized spline code


void init_spline(int *index, int NVAL, double* x_list, double* y_list,
                      char *cid, char *method_spline, int verbose,
                      double *mean, double *std_dev, int *error_flag);
void init_spline__(int *index, int *NVAL, double* x_list, double* y_list,
                      char *cid, char *method_spline, int *verbose,
                      double *mean, double *std_dev, int *error_flag);
