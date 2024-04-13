
// Created July 2021 [extracted from sntools.h]

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#define MXGENVAR_SIMEFFMAP   10

struct SIMEFFMAP {
  int NGENVAR ; // number of variables to specify map
  char   VARNAME[MXGENVAR_SIMEFFMAP][20] ;  // variable name
  char   VARSCALE[MXGENVAR_SIMEFFMAP][8] ;  // LIN, LOG or INV
  int    IFLAGSCALE[MXGENVAR_SIMEFFMAP];    // 1  , 10,    -1
  int    NBINTOT;
  int    NBIN[MXGENVAR_SIMEFFMAP];
  double VARMIN[MXGENVAR_SIMEFFMAP] ;
  double VARMAX[MXGENVAR_SIMEFFMAP] ;
  double EFFMAX ; // max efficiency in map 

  double **TMPVAL ; // temp memory to read effic. map
  double  *TMPEFF ; // idem
} SIMEFFMAP ;

GRIDMAP_DEF  SIMEFF_GRIDMAP ;



// simeff utilities (include mangled functions for fortran)
int    init_SIMEFFMAP(char *file, char *varnamesList);
double get_SIMEFFMAP(int OPTMASK, int NVAR, double *GRIDVALS);
void   malloc_SIMEFFMAP(int flag);

int    init_simeffmap__(char *file, char *varnamesList);
double get_simeff__(int *OPTMASK, int *NVAR, double *GRIDVALS);

// === END ===
