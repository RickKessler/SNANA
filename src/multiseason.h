
#define MXSEASON 100

struct {
  int   OPTMASK ;
  int   NREJECT_OUTLIER ;
  float TGAP ;
} INPUTS ;

// ---------------

void INIT_MULTISEASON(float *parList) ;
void init_multiseason__(float *parList);

void GET_MULTISEASON(char *CCID, int NOBS, double *MJD, 
		     double *FLUX, double *FLUXERR,
		     int *REJECT, int *NSEASON, double *CHI2RED, 
		     double *MJDMIN, double *MJDMAX, double *AVGFLUX);

void get_multiseason__(char *CCID, int *NOBS, double *MJD, 
		       double *FLUX, double *FLUXERR,
		       int *REJECT, int *NSEASON, double *CHI2RED,
		       double *MJDMIN, double *MJDMAX, double *AVGFLUX ) ;
