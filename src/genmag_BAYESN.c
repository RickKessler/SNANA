/**********************************************

  July 1 2022 G. Narayan

********************************************/

#include "fitsio.h"
#include "sntools.h"
#include  "genmag_BAYESN.h"
// #include "sntools_modelgrid.h" 
// #include "sntools_genSmear.h" // Aug 30 2019



// ============ MANGLED FORTRAN FUNCTIONS =============
// GN - I think we will likely need these - copied from bayesn
int init_genmag_bayesn__( char *version, int *optmask) {
  int istat;
  istat = init_genmag_bayesn ( version, *optmask) ;
  return istat;
}


/* int genmag_bayesn__(int *ifilt, double *stretch, int *nobs, 
		  double *Trest, double *rest_mags, double *rest_magerrs ) {
  int istat;
  istat = genmag_bayesn(*ifilt, *stretch, *nobs, 
			Trest, rest_mags, rest_magerrs ) ;
  return istat;
} */

// void get_lamrange_bayesn__(double *lammin, double *lammax) {
//  get_LAMRANGE_bayesn(lammin,lammax);
// }

int init_genmag_bayesn(char *version, int optmask){
    printf("Hello from fortran hell");

}
