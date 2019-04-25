
#define MXSEASON 100

struct {
  int   OPTMASK ;
  int   NREJECT_OUTLIER ;
  float TGAP ;
} INPUTS ;

// ---------------

void INIT_MULTISEASON(float *parList) ;
void init_multiseason__(float *parList);

void GET_MULTISEASON(char *CCID, int NOBS, float *MJD, 
		     float *FLUX, float *FLUXERR,
		     int *REJECT, int *NSEASON, float *CHI2RED, 
		     float *MJDMIN, float *MJDMAX, float *AVGFLUX);

void get_multiseason__(char *CCID, int *NOBS, float *MJD, 
		       float *FLUX, float *FLUXERR,
		       int *REJECT, int *NSEASON, float *CHI2RED,
		       float *MJDMIN, float *MJDMAX, float *AVGFLUX ) ;
