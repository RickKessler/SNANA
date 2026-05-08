// sntools_spline_gen.c 
//
// Generalized multi-instance spline interpolation.
// See sntools_spline_gen.h for full API documentation.
//
// Authors : A. Mitra  (May 2026) -- initial implementation
// Based on : sntools_zPDF_spline.c  by R. Chen (Jun 2022), R. Kessler 

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "sntools.h"
#include "sntools_spline.h"

// - - - - - - - - - - - - - - - - - - - - - -
// Module-level storage

SPLINE_GEN_DEF SPLINE_GEN_ARRAY[MXSPLINE_GEN];

// Private helpers
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// First-call zero-initialization of the entire slot array.
static void _spline_gen_global_init(void) {
    int i;
    for (i = 0; i < MXSPLINE_GEN; i++) {
        SPLINE_GEN_ARRAY[i].initialized = 0;
        SPLINE_GEN_ARRAY[i].acc         = NULL;
        SPLINE_GEN_ARRAY[i].spline      = NULL;
        SPLINE_GEN_ARRAY[i].NVAL        = 0;
        SPLINE_GEN_ARRAY[i].val_max     = 1.0;
    }
    SPLINE_GEN_INIT_DONE = 1;
}

// Infer evaluation mode from the slot name.
// Names containing "QUANTILE" use DERIV mode (PDF from CDF derivative).
// All others use DIRECT mode (plain interpolation).
static int _get_mode_from_name(const char *name) {
    if (strstr(name, "QUANTILE") != NULL) { return SPLINE_GEN_MODE_DERIV; }
    return SPLINE_GEN_MODE_DIRECT;
}

// Free GSL objects for a slot if they were previously allocated.
static void _free_spline_slot(SPLINE_GEN_DEF *sp) {
    if (sp->initialized) {
        if (sp->spline) { gsl_spline_free(sp->spline);       sp->spline = NULL; }
        if (sp->acc)    { gsl_interp_accel_free(sp->acc);    sp->acc    = NULL; }
        sp->initialized = 0;
    }
}

// Allocate GSL spline and accelerator for the requested method.
// Returns 0 on success, -1 on unknown method (caller should abort).
static int _alloc_gsl_spline(SPLINE_GEN_DEF *sp, int NVAL,
                              const char *method_spline, const char *cid) {
    char fnam[] = "_alloc_gsl_spline";

    sp->acc = gsl_interp_accel_alloc();

    if (strcmp(method_spline, METHOD_SPLINE_LINEAR) == 0) {
        sp->spline = gsl_spline_alloc(gsl_interp_linear, NVAL);

    } else if (strcmp(method_spline, METHOD_SPLINE_CUBIC) == 0) {
        sp->spline = gsl_spline_alloc(gsl_interp_cspline, NVAL);

    } else if (strcmp(method_spline, METHOD_SPLINE_STEFFEN) == 0) {
#ifdef GSL_INTERP_STEFFEN
        sp->spline = gsl_spline_alloc(gsl_interp_steffen, NVAL);
#else
        sprintf(c1err, "Compilation does not include STEFFEN method");
        sprintf(c2err, "Needs GSL >= 2.7  (CID=%s)", cid);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
#endif
    } else {
        sprintf(c1err, "Unknown method_spline = '%s' for CID=%s", method_spline, cid);
        sprintf(c2err, "Valid options: %s  %s  %s",
                METHOD_SPLINE_LINEAR, METHOD_SPLINE_CUBIC, METHOD_SPLINE_STEFFEN);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
    return 0;
}

// Compute mean and std_dev on a uniform grid over [xmin, xmax].
//
// DERIV mode (CDF spline, x=z, y=percentile):
//   Uses the derivative  P(x) = d[CDF]/dx  as a probability weight.
//   mean   = integral[ x * P(x) dx ] / integral[ P(x) dx ]   (PDF-weighted mean of x)
//   std    = sqrt( integral[ (x-mean)^2 * P(x) dx ] / integral[P(x) dx] )
//   val_max is set to max(P(x)) for later normalization in eval.
//
// DIRECT mode (e.g. logmass(z)):
//   Computes simple mean and std of y over the x grid.
//   mean   = mean of y(x_i) over the uniform grid
//   std    = std  of y(x_i) over the uniform grid
//   val_max is set to max(|y|) (informational; not used in eval).
//
static void _compute_stats(SPLINE_GEN_DEF *sp) {
    int    ix, NBIN = NBIN_SPLINE_GEN;
    double x, val, sum_wval, sum_w, sum_sq, val_max;
    double xmin = sp->xmin, xmax = sp->xmax, dx = sp->dx;

    double *store = (double *)malloc((NBIN + 2) * sizeof(double));

    sum_wval = sum_w = val_max = 0.0;
    ix = 0;

    for (x = xmin; x <= xmax + 0.5 * dx; x += dx) {
        if (x > xmax) x = xmax;   // clamp last step

        if (sp->mode == SPLINE_GEN_MODE_DERIV) {
            // weight = PDF(x) = d[CDF]/dx; floor at 0
            val = gsl_spline_eval_deriv(sp->spline, x, sp->acc);
            if (val < 0.0) val = 0.0;
            sum_wval += x * val;   // PDF-weighted sum of x
            sum_w    += val;
            if (val > val_max) val_max = val;

        } else {
            // direct: accumulate y values uniformly
            val       = gsl_spline_eval(sp->spline, x, sp->acc);
            sum_wval += val;
            sum_w    += 1.0;
            if (fabs(val) > val_max) val_max = fabs(val);
        }
        store[ix++] = val;
    }

    sp->val_max = (val_max > 0.0) ? val_max : 1.0;
    sp->mean    = (sum_w > 0.0)   ? (sum_wval / sum_w) : 0.0;

    // second pass: standard deviation
    sum_sq = 0.0;
    ix = 0;
    for (x = xmin; x <= xmax + 0.5 * dx; x += dx) {
        if (x > xmax) x = xmax;
        double wt = store[ix];
        double q  = (sp->mode == SPLINE_GEN_MODE_DERIV) ? x : wt;
        double dq = (sp->mode == SPLINE_GEN_MODE_DERIV) ? wt : 1.0;
        sum_sq += (q - sp->mean) * (q - sp->mean) * dq;
        ix++;
    }
    sp->std_dev = (sum_w > 0.0 && sum_sq > 0.0) ? sqrt(sum_sq / sum_w) : 0.0;

    free(store);
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Public functions
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// ============================================================
void init_spline(int index, int NVAL,
                     double *x_list, double *y_list,
                     char *cid, char *method_spline, char *name, int verbose,
                     double *mean, double *std_dev, int *error_flag) {
    // See header for full documentation.

    char fnam[] = "init_spline";
    int  i;

    // one-time global initialization
    if (!SPLINE_GEN_INIT_DONE) { _spline_gen_global_init(); }

    // ---- validate index ----
    if (index < 0 || index >= MXSPLINE_GEN) {
        sprintf(c1err, "index=%d out of range [0, %d]", index, MXSPLINE_GEN - 1);
        sprintf(c2err, "CID=%s  name=%s", cid, name);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    // ---- validate NVAL ----
    if (NVAL < 2) {
        sprintf(c1err, "NVAL=%d is too small (need >= 2)", NVAL);
        sprintf(c2err, "CID=%s  name=%s", cid, name);
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    // ---- init outputs ----
    *error_flag = 0;
    *mean       = 0.0;
    *std_dev    = 0.0;

    // ---- check x is strictly increasing ----
    for (i = 1; i < NVAL; i++) {
        if (x_list[i] <= x_list[i - 1]) {
            if (x_list[i] < 0.0) { *error_flag = -1; }
            else                  { *error_flag = -2; }
            if (verbose) {
                printf(" WARNING %s: x_list not strictly increasing at i=%d"
                       " (x[%d]=%g, x[%d]=%g) for CID=%s name=%s\n",
                       fnam, i, i-1, x_list[i-1], i, x_list[i], cid, name);
                fflush(stdout);
            }
            return;
        }
    }

    // ---- free previous allocation for this slot ----
    _free_spline_slot(&SPLINE_GEN_ARRAY[index]);

    SPLINE_GEN_DEF *sp = &SPLINE_GEN_ARRAY[index];

    // ---- store metadata ----
    snprintf(sp->name,   sizeof(sp->name),   "%s", name);
    snprintf(sp->method, sizeof(sp->method), "%s", method_spline);
    sp->NVAL  = NVAL;
    sp->mode  = _get_mode_from_name(name);
    sp->xmin  = x_list[0];
    sp->xmax  = x_list[NVAL - 1];
    sp->dx    = (sp->xmax - sp->xmin) / (double)NBIN_SPLINE_GEN;

    // ---- allocate and fit GSL spline ----
    _alloc_gsl_spline(sp, NVAL, method_spline, cid);
    gsl_spline_init(sp->spline, x_list, y_list, NVAL);
    sp->initialized = 1;

    // ---- compute mean / std / val_max on uniform grid ----
    _compute_stats(sp);

    *mean    = sp->mean;
    *std_dev = sp->std_dev;

    if (verbose) {
        printf("\t init_spline[%d] '%s'  method=%s  mode=%s  N=%d"
               "  x=[%.4f, %.4f]\n",
               index, sp->name, sp->method,
               (sp->mode == SPLINE_GEN_MODE_DERIV) ? "DERIV" : "DIRECT",
               NVAL, sp->xmin, sp->xmax);
        printf("\t   mean=%.4f  std=%.4f  val_max=%.4f  CID=%s\n",
               *mean, *std_dev, sp->val_max, cid);
        fflush(stdout);
    }

} // end init_spline


// ============================================================
double eval_spline_gen(int index, double x) {
    // Evaluate spline slot 'index' at a single query point x.
    // DIRECT: returns interpolated y(x).
    // DERIV:  returns normalized PDF P(x) = d[CDF]/dx / val_max.
    // Returns 0.0 for x outside [xmin, xmax].

    char fnam[] = "eval_spline_gen";
    double y;
    SPLINE_GEN_DEF *sp = &SPLINE_GEN_ARRAY[index];

    if (!sp->initialized) {
        sprintf(c1err, "Spline slot %d ('%s') not initialized", index, sp->name);
        sprintf(c2err, "Call init_spline before eval_spline_gen");
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    if (x < sp->xmin || x > sp->xmax) { return 0.0; }

    if (sp->mode == SPLINE_GEN_MODE_DERIV) {
        y = gsl_spline_eval_deriv(sp->spline, x, sp->acc);
        if (y < 0.0)  y = 0.0;                // no unphysical negative PDF
        y /= sp->val_max;                      // normalize so peak = 1
    } else {
        y = gsl_spline_eval(sp->spline, x, sp->acc);
    }
    return y;

} // end eval_spline_gen


// ============================================================
void eval_spline_gen_array(int index, int N_query,
                            double *x_query, double *y_out) {
    // Evaluate spline slot 'index' at each of the N_query points in x_query[].
    // Writes results into y_out[N_query].
    // Points outside [xmin, xmax] are set to 0.0 (same as eval_spline_gen).
    //
    // Example uses:
    //   -- compute PDF on a fine z grid for integration
    //   -- compute logmass(z) on a redshift grid for plotting/debugging

    int i;
    for (i = 0; i < N_query; i++) {
        y_out[i] = eval_spline_gen(index, x_query[i]);
    }

} // end eval_spline_gen_array


// ============================================================
void eval_spline_gen_integral(int index,
                               double x_center, double x_sigma,
                               double *mean_out, double *std_out) {
    // Gaussian-weighted integral of y(x) (or P(x) for DERIV mode).
    //
    // Computes:
    //   mean_out = int[ y(x) * G(x) dx ] / int[ G(x) dx ]
    //   std_out  = sqrt( int[ (y(x) - mean_out)^2 * G(x) dx ] / int[G(x) dx] )
    // where G(x) = exp( -0.5 * (x - x_center)^2 / x_sigma^2 ).
    //
    // Integration is performed over [max(xmin, x_center - 4*sigma),
    //                                min(xmax, x_center + 4*sigma)]
    // on the same grid spacing dx as the stats grid.
    //
    // Degenerate case (x_sigma <= 0): returns point estimate at x_center.
    //
    // Primary use: logmass(z) marginalized over photo-z Gaussian uncertainty.
    //   eval_spline_gen_integral(INDEX_LOGMASS, z_ph, z_ph_err,
    //                            &lm_mean, &lm_std)

    char fnam[] = "eval_spline_gen_integral";
    SPLINE_GEN_DEF *sp = &SPLINE_GEN_ARRAY[index];

    if (!sp->initialized) {
        sprintf(c1err, "Spline slot %d ('%s') not initialized", index, sp->name);
        sprintf(c2err, "Call init_spline before eval_spline_gen_integral");
        errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    // degenerate case: no uncertainty -- just return point estimate
    if (x_sigma <= 0.0) {
        *mean_out = eval_spline_gen(index, x_center);
        *std_out  = 0.0;
        return;
    }

    double sigma2  = x_sigma * x_sigma;
    double dx      = sp->dx;

    // clamp integration range to ± 4 sigma and spline domain
    double x_lo = fmax(sp->xmin, x_center - 4.0 * x_sigma);
    double x_hi = fmin(sp->xmax, x_center + 4.0 * x_sigma);

    if (x_lo >= x_hi) {
        // Gaussian window entirely outside spline domain
        *mean_out = -9.0;
        *std_out  =  0.0;
        return;
    }

    double sum_yw = 0.0, sum_w = 0.0, sum_sq = 0.0;
    double x, y, g;

    // first pass: weighted mean
    for (x = x_lo; x <= x_hi + 0.5 * dx; x += dx) {
        if (x > x_hi) x = x_hi;

        g = exp(-0.5 * (x - x_center) * (x - x_center) / sigma2);

        if (sp->mode == SPLINE_GEN_MODE_DERIV) {
            y = gsl_spline_eval_deriv(sp->spline, x, sp->acc);
            if (y < 0.0) y = 0.0;
        } else {
            y = gsl_spline_eval(sp->spline, x, sp->acc);
        }

        sum_yw += y * g;
        sum_w  += g;
    }

    *mean_out = (sum_w > 0.0) ? (sum_yw / sum_w) : -9.0;

    // second pass: weighted std dev
    for (x = x_lo; x <= x_hi + 0.5 * dx; x += dx) {
        if (x > x_hi) x = x_hi;

        g = exp(-0.5 * (x - x_center) * (x - x_center) / sigma2);

        if (sp->mode == SPLINE_GEN_MODE_DERIV) {
            y = gsl_spline_eval_deriv(sp->spline, x, sp->acc);
            if (y < 0.0) y = 0.0;
        } else {
            y = gsl_spline_eval(sp->spline, x, sp->acc);
        }

        sum_sq += (y - *mean_out) * (y - *mean_out) * g;
    }

    *std_out = (sum_w > 0.0 && sum_sq > 0.0) ? sqrt(sum_sq / sum_w) : 0.0;

} // end eval_spline_gen_integral


// ============================================================
void dump_spline_gen(int index) {
    // Print diagnostic summary of slot 'index' to stdout.

    if (!SPLINE_GEN_INIT_DONE) {
        printf(" xxx dump_spline_gen: not yet initialized\n");
        return;
    }
    SPLINE_GEN_DEF *sp = &SPLINE_GEN_ARRAY[index];
    printf(" xxx dump_spline_gen[%d]:\n", index);
    printf(" xxx   name='%s'  method=%s  mode=%s  initialized=%d\n",
           sp->name, sp->method,
           (sp->mode == SPLINE_GEN_MODE_DERIV) ? "DERIV" : "DIRECT",
           sp->initialized);
    printf(" xxx   NVAL=%d  xmin=%.4f  xmax=%.4f  dx=%.4f\n",
           sp->NVAL, sp->xmin, sp->xmax, sp->dx);
    printf(" xxx   mean=%.4f  std_dev=%.4f  val_max=%.4f\n",
           sp->mean, sp->std_dev, sp->val_max);
    fflush(stdout);

} // end dump_spline_gen


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Fortran wrappers
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// gfortran appends hidden character lengths after all explicit arguments,
// in the order the character arguments appear in the call.
// We accept len_cid and len_method explicitly (matching the Fortran call site);
// the hidden length for any literal string (e.g. "LOGMASS_ZGRID") is appended
// by gfortran and lands after len_method -- it is unused here since all strings
// are null-terminated with // char(0) at the call site.

void init_spline__(int *index, int *NVAL,
                        double *x_list, double *y_list,
                        char *cid, char *method_spline, char *name,
                        int *verbose,
                        double *mean, double *std_dev, int *error_flag,
                        int len_cid, int len_method) {
    (void)len_cid;     // strings are null-terminated; lengths not needed
    (void)len_method;
    init_spline(*index, *NVAL, x_list, y_list,
                    cid, method_spline, name, *verbose,
                    mean, std_dev, error_flag);
}

double eval_spline_gen__(int *index, double *x) {
    return eval_spline_gen(*index, *x);
}

void eval_spline_gen_array__(int *index, int *N_query,
                              double *x_query, double *y_out) {
    eval_spline_gen_array(*index, *N_query, x_query, y_out);
}

void eval_spline_gen_integral__(int *index,
                                 double *x_center, double *x_sigma,
                                 double *mean_out, double *std_out) {
    eval_spline_gen_integral(*index, *x_center, *x_sigma, mean_out, std_out);
}

void dump_spline_gen__(int *index) {
    dump_spline_gen(*index);
}
