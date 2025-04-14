// Created Apr 2025
// Utilities to read/write python-style npz files 
// Initial use is for wfit.c to read npz-formatted covariance file 
// from create_covariance.py.
// Use package from https://github.com/rogersce/cnpy

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#include "sntools_npz.h"
#include "cnpy.h"



int read_npz_covmat(char *npz_file, double *array1d) {

  // Created Apr 14 2024
  // Read cov matrix from npz file of the form
  //    nsn = nsn
  //    cov = cov(utria).float
  //
  // The cov in npz file is float (4 byte) and upper triangular half
  // to save disk space. 
  // However, the returned cov is the full N x N in double precision.

  int k0, k1, j1du=0, j1d_tmp, NSN;
  double cov;
  bool CHECK_COV_LOADING = false;   int *n_cov_load;
  char fnam[] = "read_npz_covmat" ;

  // ------------ BEGIN --------------

  // load the entire npz file
  cnpy::npz_t  my_npz = cnpy::npz_load(npz_file);

  // pick off NSN
  cnpy::NpyArray nsn_npz = my_npz["nsn"];
  NSN = *nsn_npz.data<int>();

  // read upper-triangle part of cov
  cnpy::NpyArray cov_utri_npz = my_npz["cov"];
  float *cov_utri = cov_utri_npz.data<float>() ;

  // - - - - - - -
  // load full NSN x NSN output array and convert float to double

  if ( CHECK_COV_LOADING ) { 
    n_cov_load = (int*) malloc( NSN*NSN * sizeof(int) );
    for(j1d_tmp = 0; j1d_tmp < NSN*NSN; j1d_tmp++ ) { n_cov_load[j1d_tmp] = false; }
  }

  for (k0=0; k0 < NSN; k0++ ) {
    for (k1=k0; k1 < NSN; k1++ ) {
      cov    = cov_utri[j1du];
       
      j1d_tmp = k0*NSN + k1;
      array1d[j1d_tmp] = (double)cov ;
      if ( CHECK_COV_LOADING ) { n_cov_load[j1d_tmp]++ ; }

      if ( k0 != k1 ) {
	j1d_tmp = k1*NSN + k0;
	array1d[j1d_tmp] = (double)cov ;
	if ( CHECK_COV_LOADING ) { n_cov_load[j1d_tmp]++ ; }
      }

      j1du+=1 ;
    }
  }


  // - - - - - - - -

  if ( CHECK_COV_LOADING ) {    
    for(j1d_tmp = 0; j1d_tmp < NSN*NSN; j1d_tmp++ ) {
      int n = n_cov_load[j1d_tmp] ;
      if ( n != 1 ) {
	printf(" xxx %s ERROR: j1d=%d  n_cov_load = %d \n",
	       fnam, j1d_tmp, n ); fflush(stdout);
      }
    }
  }

  // - - - -
  return NSN;

} // read_npz_array

