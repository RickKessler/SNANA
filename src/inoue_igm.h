#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>


//// Maximum number of lines in template file
#define NTEMPLMAX (long)1.e5 

char LAF_FILE[1024];
char DLA_FILE[1024];

//// IGM parametrization from Inoue et al. 2014
double *lam1, *ALAF1, *ALAF2, *ALAF3, *ADLA1, *ADLA2;
int NA;
void read_Inoue_coeffs();
double tLSLAF();
double tLCLAF();
double tLSDLA();
double tLCDLA();

