// genmag_mlcs2k2.h

#define MAXDAYS_MLCS2k2 120        // should be the same as MXEPCOV in sndata.h
#define MAXFILT_MLCS2k2 9          // 012345678 <==> UBVRIYJHK
#define NFILT_MLCS2k2_REQUIRED 5   // require UBVRI templates

#define H0_MLCS2k2_TRAIN 65.0  // H0 where MLCS was trained

// define max lambda vs. IR filter
#define MINLAM_MLCS2k2_U   3200. ;
#define MAXLAM_MLCS2k2_I   9500. ;
#define MAXLAM_MLCS2k2_Y  10300. ;
#define MAXLAM_MLCS2k2_J  14000. ;
#define MAXLAM_MLCS2k2_H  18000. ;
#define MAXLAM_MLCS2k2_K  24500. ;

/* Global variables */

struct TEMPLATE_MAGS_MLCS2k2 {
  double MAGOFF;  // rest mag = MAGOFF
  double D1, D2;  //  + D1*DELTA + D2*(DELTA**2) 
  double MAGERR;  // error on mag
} TEMPLATE_MAGS_MLCS2k2[MAXDAYS_MLCS2k2][MAXFILT_MLCS2k2];



double TEMPLATE_DAYS_MLCS2k2[MAXDAYS_MLCS2k2];
int    NDAY_MLCS2k2;
int    NFILT_MLCS2k2; // either 5 or 6-9
double TMINDAY_MLCS2k2, TMAXDAY_MLCS2k2 ;
char   FILTSTRING_MLCS2k2[20]  ;
char   PATHMODEL_MLCS2k2[200] ;
int    IDAYPEAK_MLCS2k2 ;  // day-index for peakmag
int    USE_PEAKERR_ONLY ;  // 1=yes(test only);  0=no

float MLCS2k2_COVAR[MAXFILT_MLCS2k2][MAXDAYS_MLCS2k2][MAXFILT_MLCS2k2][MAXDAYS_MLCS2k2] ;

int NONZERO_COVAR;  // number of non-zero covariance elements

double LAMRANGE_MLCS2k2[2]; // lambda-range of MLCS2k2 model

int    QWGT_FLAG;
#define QWGT_TRESTMIN  (double)25.0
#define QWGT_DELTAMIN  (double)0.3 
#define QWGT_TRESTSIG  (double)5.0  // half-Gaussian sigma for QWGT
double QWGT_MLCS2k2 ;
double MOFF_MLCS2k2 ;

// ----------------------------------------
//    function declarations
// ----------------------------------------

int init_genmag_mlcs2k2(char *version, char *covFile
		     ,float scale_covar
		     ,char *filtlist );  // output

int genmag_mlcs2k2(int ifilt, double delta, int nobs, 
		double *rest_dates, double *rest_magval, double *rest_magerr);

int gencovar_mlcs2k2(int matsize, int *ifilt, double *rest_epoch, 
		  double *covar );

void rd_mlcs2k2_cov(char *covFile, float scale_covar);  // read cov matrix
void dmp_mlcs2k2_cov(void) ;
void dmp_mlcs2k2_magerr(void);

int  ifilt_mlcs2k2(char *cfilt) ;

void mag_extrap_mlcs2k2(int ifilt, double delta, double Trest, 
		     double *mag, double *magerr );

void JDATES ( double Trest, int *j1, int *j2 );
int  mlcs2k2_Tmin(void);
int  mlcs2k2_Tmax(void);

void set_LAMRANGE_mlcs2k2(void);
void get_LAMRANGE_mlcs2k2(double *lammin, double *lammax);
void get_MPQ_mlcs2k2(int j, int ifilt, double *M, double *P, double *Q);
double get_QWGT_mlcs2k2(double Trest, double delta);

// END
