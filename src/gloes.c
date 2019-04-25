/* This is the GLoEs algorithm (Gaussian Local Estimation).  Original 1D
 * algorightm by Barry Madore (OCIW).  Further expanded to incorporate
 * an adaptive weight algorithm and work in 2D by Chris Burns (OCIW).  
 * Ported to C by Chris Burns*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_sort_vector.h>
/* #include "debugs.h" */

#define ORDER 2

/* This defines a maximum distance to consider for SVD.  Data points farther than
 * MAXDIST-sigma  will be rejected, allowing the matrix to be smaller */
#define MAXDIST 20

int nk(int k) {
   if(k==2){
      return(3);
   } else {
      return(k*k - nk(k-1));
   }
}

int fit2Dpoly(gsl_vector *x, gsl_vector *y, gsl_vector *z, gsl_vector *w, double x0, 
      double y0, int n, int k, gsl_vector *c, gsl_vector *cerror, double *chisq) {
   /* Fit a 2D polynomial of degree k to the data (x,y,z, weight) (length n)
    * around the point x0,y0.  The coefficients are returned as 2k+1 length double
    * *c and are in the sense: y(x) = c[0] + c[1]*(x-xo) + c[2]*(x-x0)^2 + ...
    * + c[k]*(x-x0)^k.  The errors in c[] are retured in cerror[] and chisq has
    * the minimum chisq of the solution.  The return value of 0 indicates
    * success.  The calling function is responsible for ensuring that memory
    * has been allocated for c and  cerror (each length k+1) */

   int i, j, l, m;
   int count;
   gsl_matrix *X, *cov;
   gsl_vector_int *perm;
   gsl_multifit_linear_workspace *work;
   gsl_vector *xsub, *ysub, *zsub, *wsub;

   perm = gsl_vector_int_alloc(n);
   count = 0;
   for (i = 0 ; i < n ; ++i ) {
      if (gsl_vector_get(w, i) > 0) {
         gsl_vector_int_set(perm, count, i);
         count += 1;
      }
   }
   xsub = gsl_vector_alloc(count);
   ysub = gsl_vector_alloc(count);
   zsub = gsl_vector_alloc(count);
   wsub = gsl_vector_alloc(count);
   for (i = 0 ; i < count ; ++i ) {
      j = gsl_vector_int_get(perm, i);
      gsl_vector_set(xsub, i, gsl_vector_get(x, j));
      gsl_vector_set(ysub, i, gsl_vector_get(y, j));
      gsl_vector_set(zsub, i, gsl_vector_get(z, j));
      gsl_vector_set(wsub, i, gsl_vector_get(w, j));
   }

   X = gsl_matrix_alloc(count, nk(k+1));
   cov = gsl_matrix_alloc(nk(k+1),nk(k+1));

   /* loop over data points */
   for (i = 0; i < count ; ++i) {
      m = 0;   /* a running index of parameters */
      /* outer loop over powers of y */
      for (j = 0 ; j < k+1 ; ++j){
         /* inner loop over powers of x */
         for (l = 0 ; l < k + 1 ; ++l) {
            /* fprintf(stderr, "j=%d, l=%d, m=%d\n", j, l, m); */
            if (j + l > k) continue;
            if (j == 0 && l == 0) {
               gsl_matrix_set(X, i, m, 1.0);
            }
            else if (j==1 && l == 0) {
               gsl_matrix_set(X, i, m, gsl_vector_get(xsub,i)-x0);
            }
            else if (j == 0 && l == 1) {
               gsl_matrix_set(X, i, m, gsl_vector_get(ysub,i)-y0);
            }
            else if (j == 1 && l == 1) {
               gsl_matrix_set(X, i, m, (gsl_vector_get(xsub,i)-x0)*(gsl_vector_get(ysub,i)-y0));
            }
            else {
               gsl_matrix_set(X, i, m, pow(gsl_vector_get(xsub,i)-x0, j)*pow(gsl_vector_get(ysub,i)-y0, l));
            }
            m = m + 1;
         }
      }
   }
   work = gsl_multifit_linear_alloc(count, nk(k+1));
   gsl_multifit_wlinear(X, wsub, zsub, c, cov, chisq, work);
   gsl_multifit_linear_free(work);

   for (i = 0 ; i < nk(k+1) ; ++i) {
      gsl_vector_set(cerror, i, sqrt(gsl_matrix_get(cov, i, i)));
   }
   gsl_matrix_free(X);
   gsl_matrix_free(cov);
   gsl_vector_free(xsub);
   gsl_vector_free(ysub);
   gsl_vector_free(zsub);
   gsl_vector_free(wsub);
   gsl_vector_int_free(perm);
   return 0;
}


int fitpoly(gsl_vector *x, gsl_vector *y, gsl_vector *w, 
	    double x0, int n, int k, 
	    gsl_vector *c, gsl_vector *cerror, double *chisq) {

   /* Fit a 1D polynomial of degree k to the data (x,y,weight) (length n)
    * around the point x0.  The coefficients are returned as (k+1)^2 
    * length double
    * c and are in the sense: 
    *    y(x) = \sum c[j*(k+1)+i] (x - x0)^i*(y-y0)^j
    * The errors in c[] are retured in cerror[] and chisq has
    * the minimum chisq of the solution.  The return value of 0 indicates
    * success.  The calling function is responsible for ensuring that memory
    * has been allocated for c and  cerror (each length (k+1)^2) */

   int i, j;
   gsl_matrix *X, *cov;
   gsl_multifit_linear_workspace *work;

   X = gsl_matrix_alloc(n, k+1);
   cov = gsl_matrix_alloc(k+1,k+1);

   for (i = 0; i < n ; ++i) {
      for (j = 0 ; j < k+1 ; ++j){
         if (j == 0) {
            gsl_matrix_set(X, i, j, 1.0);
         }
         else if (j==1) {
            gsl_matrix_set(X, i, j, gsl_vector_get(x,i)-x0);
         }
         else {
            gsl_matrix_set(X, i, j, pow(gsl_vector_get(x,i)-x0, j));
         }
      }
   }
   work = gsl_multifit_linear_alloc(n, k+1);
   gsl_multifit_wlinear(X, w, y, c, cov, chisq, work);
   gsl_multifit_linear_free(work);

   for (i = 0 ; i < k+1 ; ++i) {
      gsl_vector_set(cerror, i, sqrt(gsl_matrix_get(cov, i, i)));
   }
   gsl_matrix_free(X);
   gsl_matrix_free(cov);
   return 0;
}

int gloes_one_sigma( double *xp, double *yp, double *wp, int np, double *x, int n, 
      double sigma, double *y, double *ey) {
   /* Use the GLoEs algorithm to interpolate a given set of data (xp, yp, wp) (length np)
    * at the points in vector x (length n) using a window function with constant width
    * sigma.  Results are returned as an array y and errors in ey, assumed to be length n.
    * CAlling function is responsible for allocating memory for these variables.*/

   gsl_vector *weight;
   gsl_vector *gxp, *gyp;
   double max_w = 0;
   double w;
   int i, j;
   gsl_vector *c, *cerror;
   double dists2;
   double sigma2=sigma*sigma;
   int res;
   double chisq;

   /* Allocate memory */
   weight = gsl_vector_alloc(np);
   gxp = gsl_vector_alloc(np);
   gyp = gsl_vector_alloc(np);
   c = gsl_vector_alloc(3);
   cerror = gsl_vector_alloc(3);

   /* setup the anchor points */
   for (i = 0 ; i < np ; ++i) {
      gsl_vector_set(gxp, i, xp[i]);
      gsl_vector_set(gyp, i, yp[i]);
   }
    
   /* Now cycle through the evalutaion points */
   for (i = 0 ; i < n ; ++i) {
      /* cycle through the anchor points */
      max_w = 0;
      for (j = 0 ; j < np  ; ++j) {
         dists2 = (xp[j] - x[i])*(xp[j] - x[i]);
         w = exp(-0.5*dists2/sigma2);
         if (w > max_w) max_w = w;
         gsl_vector_set(weight, j, w);
      }
      /* now normalize by max_w */
      gsl_vector_scale(weight, 1.0/max_w);

      /* apply the data weights */
      for (j = 0 ; j < np ; ++j) {
         gsl_vector_set(weight, j, gsl_vector_get(weight, j)*wp[j]);
      }
      /* gsl_vector_fprintf(stdout, weight, "%g"); */

      /* Finally, fit the polynomial */
      res = fitpoly(gxp, gyp, weight, x[i], np, 2, c, cerror, &chisq);
      y[i] = gsl_vector_get(c,0);
      ey[i] = gsl_vector_get(cerror, 0);
   }
   return(0);
   gsl_vector_free(weight);
   gsl_vector_free(gxp);
   gsl_vector_free(gyp);
   gsl_vector_free(c);
   gsl_vector_free(cerror);
}

int gloes_n_sigma( double *xp, double *yp, double *wp, int np, double *x, int n, 
      double *sigma, double *y, double *ey) {
   /* Use the GLoEs algorithm to interpolate a given set of data (xp, yp, wp) (length np)
    * at the points in vector x (length n) using a window function with variable width
    * sigma[].  Results are returned as an array y and errors in ey, assumed to be length n.
    * CAlling function is responsible for allocating memory for these variables.*/

   gsl_vector *weight;
   gsl_vector *gxp, *gyp;
   double max_w = 0;
   double w;
   int i, j;
   gsl_vector *c, *cerror;
   double dists2;
   int res;
   double chisq;

   /* Allocate memory */
   weight = gsl_vector_alloc(np);
   gxp = gsl_vector_alloc(np);
   gyp = gsl_vector_alloc(np);
   c = gsl_vector_alloc(3);
   cerror = gsl_vector_alloc(3);

   /* setup the anchor points */
   for (i = 0 ; i < np ; ++i) {
      gsl_vector_set(gxp, i, xp[i]);
      gsl_vector_set(gyp, i, yp[i]);
   }
    
   /* Now cycle through the evalutaion points */
   for (i = 0 ; i < n ; ++i) {
      /* cycle through the anchor points */
      max_w = 0;
      for (j = 0 ; j < np  ; ++j) {
         dists2 = (xp[j] - x[i])*(xp[j] - x[i]);
         w = exp(-0.5*dists2/(sigma[i]*sigma[i]));
         if (w > max_w) max_w = w;
         gsl_vector_set(weight, j, w);
      }
      /* now normalize by max_w */
      gsl_vector_scale(weight, 1.0/max_w);

      /* apply the data weights */
      for (j = 0 ; j < np ; ++j) {
         gsl_vector_set(weight, j, gsl_vector_get(weight, j)*wp[j]);
      }
      /* gsl_vector_fprintf(stdout, weight, "%g"); */

      /* Finally, fit the polynomial */
      res = fitpoly(gxp, gyp, weight, x[i], np, 2, c, cerror, &chisq);
      y[i] = gsl_vector_get(c,0);
      ey[i] = gsl_vector_get(cerror, 0);
   }
   gsl_vector_free(weight);
   gsl_vector_free(gxp);
   gsl_vector_free(gyp);
   gsl_vector_free(c);
   gsl_vector_free(cerror);
   return(0);
}

int gloes_n( double *xp, double *yp, double *wp, int np, double *x, int n, 
      int N, double *y, double *ey) {
   /* Use the GLoEs algorithm to interpolate a given set of data (xp, yp, wp)
    * (length np) at the points in vector x (length n) using a window function
    * with variable width determined by keeping N data points within the FWHM
    * of the gaussian window.  Results are returned as an array y and errors in
    * ey, assumed to be length n.  CAlling function is responsible for
    * allocating memory for these variables.*/

   gsl_vector *weight;
   gsl_vector *gxp, *gyp;
   double max_w = 0;
   double w;
   int i, j, k;
   gsl_vector *c, *cerror;
   gsl_vector *dists2;
   int res;
   double sigma;
   gsl_permutation *index;
   double chisq;

   /* Allocate memory */
   weight = gsl_vector_alloc(np);
   dists2 = gsl_vector_alloc(np);
   index = gsl_permutation_alloc(np);
   gxp = gsl_vector_alloc(np);
   gyp = gsl_vector_alloc(np);
   c = gsl_vector_alloc(3);
   cerror = gsl_vector_alloc(3);

   /* setup the anchor points */
   for (i = 0 ; i < np ; ++i) {
      gsl_vector_set(gxp, i, xp[i]);
      gsl_vector_set(gyp, i, yp[i]);
   }
    
   /* Now cycle through the evalutaion points */
   for (i = 0 ; i < n ; ++i) {
      /* cycle through the anchor points */
      for (j = 0 ; j < np  ; ++j) {
         gsl_vector_set(dists2, j, (xp[j] - x[i])*(xp[j] - x[i]));
      }
      /* now find the permuation index that sorts the dists2 */
      gsl_sort_vector_index(index, dists2);

      /* We are interested in the index of the Nth distance */
      k = gsl_permutation_get(index, N);
      sigma = gsl_vector_get(dists2, k)/1.665;
      max_w = 0;
      for (j = 0 ; j < np ; ++j) {
         w = exp(-0.5*gsl_vector_get(dists2,j)/(sigma*sigma));
         if (w > max_w) max_w = w;
         gsl_vector_set(weight, j, w);
      }
      /* now normalize by max_w */
      gsl_vector_scale(weight, 1.0/max_w);

      /* apply the data weights */
      for (j = 0 ; j < np ; ++j) {
         gsl_vector_set(weight, j, gsl_vector_get(weight, j)*wp[j]);
      }
      /* gsl_vector_fprintf(stdout, weight, "%g"); */

      /* Finally, fit the polynomial */
      res = fitpoly(gxp, gyp, weight, x[i], np, 2, c, cerror, &chisq);
      y[i] = gsl_vector_get(c,0);
      ey[i] = gsl_vector_get(cerror, 0);
   }

   gsl_vector_free(weight);
   gsl_vector_free(dists2);
   gsl_vector_free(gxp);
   gsl_vector_free(gyp);
   gsl_vector_free(c);
   gsl_vector_free(cerror);
   return(0);
}

int gloes2D_n_sigma( double *xp, double *yp, double *zp, double *wp, int np, double *x, 
      double *y, int n, double *sigmax, double *sigmay, double *z, double *ez) {
   /* Use the GLoEs algorithm in 2D to interpolate a given set of data (xp, yp, zp,
    * wp) (length np) at the points in vectors x,y (length n) using a window
    * function with variable width defined by sigmax and sigmay.  Results are
    * returned as an array z and errors in ez, assumed to be length n.  CAlling
    * function is responsible for allocating memory for these variables.*/

   gsl_vector *weight;
   gsl_vector *gxp, *gyp, *gzp;
   double sum_w = 0;
   double w;
   int i, j;
   gsl_vector *c, *cerror;
   gsl_vector *dists2;
   int res;
   double chisq;
   int w_count = 0;

   /* Allocate memory */
   weight = gsl_vector_alloc(np);
   dists2 = gsl_vector_alloc(np);
   gxp = gsl_vector_alloc(np);
   gyp = gsl_vector_alloc(np);
   gzp = gsl_vector_alloc(np);
   c = gsl_vector_alloc(nk(ORDER+1));
   cerror = gsl_vector_alloc(nk(ORDER+1));

   /* setup the anchor points */
   for (i = 0 ; i < np ; ++i) {
      gsl_vector_set(gxp, i, xp[i]);
      gsl_vector_set(gyp, i, yp[i]);
      gsl_vector_set(gzp, i, zp[i]);
   }
    
   /* Now cycle through the evalutaion points */
   for (i = 0 ; i < n ; ++i) {
      /* cycle through the anchor points */
      for (j = 0 ; j < np  ; ++j) {
         gsl_vector_set(dists2, j, (xp[j] - x[i])*(xp[j] - x[i])/(sigmax[i]*sigmax[i])+
               (yp[j] - y[i])*(yp[j] - y[i])/(sigmay[i]*sigmay[i]));
      }
      sum_w = 0;
      w_count = 0;
      for (j = 0 ; j < np ; ++j) {
         if (gsl_vector_get(dists2,j) > MAXDIST) {
            w = 0;
         } else {
            w = exp(-0.5*gsl_vector_get(dists2,j));
            w_count += 1;
         }
         sum_w += w;
         gsl_vector_set(weight, j, w);
      }
      if (w_count < 7) {
         /* the evaluation point is too far from the anchor points */
         z[i] = 0.0;
         ez[i] = -1.0;
      } else {
         /* now normalize by sum_w */
         gsl_vector_scale(weight, 1.0/sum_w);
 
         /* apply the data weights */
         for (j = 0 ; j < np ; ++j) {
            gsl_vector_set(weight, j, gsl_vector_get(weight, j)*wp[j]);
         }
 
         /* Finally, fit the polynomial */
         res = fit2Dpoly(gxp, gyp, gzp, weight, x[i], y[i], np, ORDER, c, cerror, &chisq);
         z[i] = gsl_vector_get(c,0);
         ez[i] = gsl_vector_get(cerror, 0)*sqrt(chisq/(w_count - ORDER));
      }
   }

   gsl_vector_free(weight);
   gsl_vector_free(dists2);
   gsl_vector_free(gxp);
   gsl_vector_free(gyp);
   gsl_vector_free(gzp);
   gsl_vector_free(c);
   gsl_vector_free(cerror);
   return(0);
}

/* int main (int argc, char **argv) {

   double xp[100], yp[100], zp[100], wp[100];
   int i, res;
   double sigma;
   double x[1], y[1], ey[1];
   double chisq;

   x[0] = atof(argv[1]);
   sigma = atof(argv[2]);

   for (i = 0 ; i < 10 ; ++i) { 
      xp[i] = i*2*3.14159/10.0;
      yp[i] = sin(xp[i]);
      wp[i] = 0.001;
   }
   res  = gloes_n(xp, yp, wp, 10, x, 1, 3, y, ey);
   printf("Interpolation:  f(%f) = %f +/- %f\n", x[0], y[0], ey[0]);
   printf("Should be:  f(x) = sin(%f) = %f\n", x[0], sin(x[0]));
} */

