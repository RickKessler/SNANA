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

int read_npz_array(char *npz_file, double *array1d) {

  int size_array = 0;
  char fnam[] = "read_npz_array" ;

  // load the entire npz file
  cnpy::npz_t  my_npz = cnpy::npz_load(npz_file);

  debugexit(fnam);
  return size_array;

} // read_npz_array

