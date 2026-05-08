// sntools_spline_gen.h
//
// Generalized multi-instance spline interpolation.
//
// Core capability: given any arrays x[N] and y[N], fit a spline and then
// evaluate it at:
//   - a single query point x            --> eval_spline_gen(index, x)
//   - an array of query points x[M]     --> eval_spline_gen_array(index, M, x, y_out)
//
// Additional capability for photo-z marginalization of a y(z) quantity:
//   - Gaussian-weighted mean/std of y(x) --> eval_spline_gen_integral(...)
//   Used for: given logmass(z_grid), return <logmass> marginalized over
//   a Gaussian photo-z PDF G(z; z_ph, z_ph_err).
//
// Supports two evaluation modes, inferred automatically from the 'name' argument:
//
//   SPLINE_GEN_MODE_DIRECT  (default; e.g. name = "LOGMASS_ZGRID"):
//     Spline stores y(x) directly.
//     eval returns the interpolated y value.
//
//   SPLINE_GEN_MODE_DERIV   (triggered when name contains "QUANTILE"):
//     Spline stores a CDF: x = redshift, y = cumulative probability [0,1].
//     eval returns  d[CDF]/dx / val_max  -- i.e. the normalized PDF P(x).
//     This matches and replaces legacy eval_zPDF_spline().
//
// Multiple spline slots are maintained simultaneously (MXSPLINE_GEN slots).
// Each is addressed by an integer index, consistent with snana.F90 constants:
//   INDEX_SPLINE_QUANTILE_ZPHOT = 1   (DERIV  mode)
//   INDEX_SPLINE_LOGMASS_ZGRID  = 2   (DIRECT mode)
//
// Interpolation methods supported: LINEAR, CUBIC, STEFFEN.
// STEFFEN requires GSL >= 2.7 (available on Perlmutter; guarded by #ifdef).
//
// Design principles:
//   - No hidden global "current index" state for eval; index is always explicit.
//   - Style and error handling follows sntools_zPDF_spline.c (R.Chen Jun 2022).
//   - Fortran (gfortran) wrappers provided with trailing __ convention.
//
// Authors : A. Mitra  (May 2026) -- initial implementation
// Based on : sntools_zPDF_spline.c  by R. Chen (Jun 2022), R. Kessler (updates)


#ifndef SNTOOLS_SPLINE_GEN_H
#define SNTOOLS_SPLINE_GEN_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// - - - - - - - - - - - - - - - - - - - - - -
// Method name constants (guarded: also defined in sntools_zPDF_spline.h)
#ifndef METHOD_SPLINE_LINEAR
#define METHOD_SPLINE_LINEAR   "LINEAR"
#define METHOD_SPLINE_CUBIC    "CUBIC"
#define METHOD_SPLINE_STEFFEN  "STEFFEN"
#endif

// GSL STEFFEN availability flag (guarded)
#ifndef GSL_INTERP_STEFFEN
#define GSL_INTERP_STEFFEN
#endif

// - - - - - - - - - - - - - - - - - - - - - -
// Capacity and mode constants

#define MXSPLINE_GEN          10    // max simultaneous spline instances
#define SPLINE_GEN_MODE_DIRECT  0   // eval returns interpolated y(x)
#define SPLINE_GEN_MODE_DERIV   1   // eval returns d[CDF]/dx / val_max

#define NBIN_SPLINE_GEN        100  // internal integration grid for stats


// - - - - - - - - - - - - - - - - - - - - - -
// Spline instance struct (exposed for transparency and debugging)

typedef struct {
    char   name[40];         // user-supplied name, e.g. "LOGMASS_ZGRID"
    char   method[40];       // interpolation method: LINEAR / CUBIC / STEFFEN
    int    mode;             // SPLINE_GEN_MODE_DIRECT or SPLINE_GEN_MODE_DERIV
    int    NVAL;             // number of input (x, y) knots
    double xmin, xmax;       // range of input x array
    double dx;               // step size on internal integration grid
    double mean;             // mean of distribution (see below for definition)
    double std_dev;          // std dev of distribution
    double val_max;          // normalization factor (for DERIV: max PDF value)
    int    initialized;      // 1 if ready to eval, 0 otherwise
    gsl_interp_accel *acc;   // GSL acceleration object
    gsl_spline       *spline; // GSL spline object
} SPLINE_GEN_DEF;

static int SPLINE_GEN_INIT_DONE = 0;  // one-time zeroing flag   

// Exposed array (read-only from caller; do not write directly)
extern SPLINE_GEN_DEF SPLINE_GEN_ARRAY[MXSPLINE_GEN];

// - - - - - - - - - - - - - - - - - - - - - -
// C function prototypes

// Initialize spline in slot 'index'.
// Inputs:
//   index        : slot index in SPLINE_GEN_ARRAY [0 .. MXSPLINE_GEN-1]
//   NVAL         : number of knots (>= 2)
//   x_list[NVAL] : independent variable; must be strictly increasing
//   y_list[NVAL] : dependent variable
//   cid          : event ID string for error messages (null-terminated)
//   method_spline: "LINEAR", "CUBIC", or "STEFFEN"
//   name         : descriptive label; controls eval mode (null-terminated)
//   verbose      : > 0 prints summary
// Outputs:
//   mean         : mean of the quantity (PDF-weighted z for DERIV; mean y for DIRECT)
//   std_dev      : corresponding standard deviation
//   error_flag   :  0  --> ok
//                  -1  --> x_list contains negative values (bad redshifts)
//                  -2  --> x_list is not strictly increasing
void init_spline(int index, int NVAL,
                     double *x_list, double *y_list,
                     char *cid, char *method_spline, char *name, int verbose,
                     double *mean, double *std_dev, int *error_flag);

// Evaluate spline at a single x.
// DIRECT mode: returns interpolated y(x).
// DERIV  mode: returns d[CDF]/dx normalized to val_max -- the PDF P(x).
// Returns 0.0 if x is outside [xmin, xmax].
double eval_spline_gen(int index, double x);

// Evaluate spline at an array of query points.
// Fills y_out[N_query] with eval_spline_gen(index, x_query[i]) for each i.
// Useful for computing a PDF or logmass curve over a redshift grid.
void eval_spline_gen_array(int index, int N_query,
                            double *x_query, double *y_out);

// Gaussian-weighted integral of y(x) over the spline domain.
// For DIRECT mode: marginalizes y(x) over G(x; x_center, x_sigma).
//   mean_out = int[ y(x) * G(x) dx ] / int[ G(x) dx ]
//   std_out  = sqrt( int[ (y(x)-mean_out)^2 * G(x) dx ] / int[G(x) dx] )
// Primary use: logmass(z) marginalized over photo-z uncertainty:
//   eval_spline_gen_integral(INDEX_LOGMASS, z_ph, z_ph_err, &lm_mean, &lm_std)
// Also defined for DERIV mode (convolution of PDF with Gaussian).
// If x_sigma <= 0, returns point estimate at x_center and std_out = 0.
void eval_spline_gen_integral(int index,
                               double x_center, double x_sigma,
                               double *mean_out, double *std_out);

// Print a diagnostic summary of spline slot 'index' to stdout.
void dump_spline_gen(int index);

// - - - - - - - - - - - - - - - - - - - - - -
// Fortran (gfortran) wrappers  [ trailing __ naming convention ]
//
// Note on character arguments:
//   gfortran passes hidden character lengths as extra integer arguments
//   appended after all explicit args, in the order the char args appear.
//   We accept len_cid and len_method explicitly; the length for the
//   literal 'name' string is appended automatically by gfortran and ignored.

void   init_spline__(int *index, int *NVAL,
                          double *x_list, double *y_list,
                          char *cid, char *method_spline, char *name,
                          int *verbose,
                          double *mean, double *std_dev, int *error_flag,
                          int len_cid, int len_method);

double eval_spline_gen__(int *index, double *x);

void   eval_spline_gen_array__(int *index, int *N_query,
                                double *x_query, double *y_out);

void   eval_spline_gen_integral__(int *index,
                                   double *x_center, double *x_sigma,
                                   double *mean_out, double *std_out);

void   dump_spline_gen__(int *index);

#endif  // SNTOOLS_SPLINE_GEN_H
