// genmag_snoopy.h


#define MXEP_SNOOPY  2000  // max epochs per filter, summed over templates
#define MXEXCLUDE_SNOOPY 20
#define MXERR_SNOOPY     10000  // max size of error map (per filter)

#define FIX_MODELERR 0.0    // 0 => use model errors from Burns' map

#define MXFILT_SNOOPY    12  // max number of rest-frame filters to define
#define IFILT_SNOOPY_MIN  1
#define IFILT_SNOOPY_MAX  9
#define IFILT_SNOOPY_B   1
#define IFILT_SNOOPY_V   2
#define IFILT_SNOOPY_u   3
#define IFILT_SNOOPY_g   4
#define IFILT_SNOOPY_r   5
#define IFILT_SNOOPY_i   6
#define IFILT_SNOOPY_Y   7
#define IFILT_SNOOPY_J   8
#define IFILT_SNOOPY_H   9 
#define IFILT_SNOOPY_K  -999  // undefined

#define DM15REF_SNOOPY 1.1 //  SN calib ~ b*(dm15 - DM15REF)

#define MINLAM_SNOOPY   2900.0  // blue edge of u
#define MAXLAM_SNOOPY  20000.0  // H band

char SNOOPYNAME[20] ;
char FILTLIST_SNOOPY[20] ; 
char SNOOPY_MODELPATH[200];

SNGRID_DEF SNGRID_SNOOPY  ;

// =============================================
//   global snoopy arrays to contain the data 
// =============================================

struct SNOOPY_TEMPLATE {
  int  NSN ;             // number of SNe per filter
  int  NEPOCH ;          // number of epochs summed over SNe
  char FILTERNAME[2];    // 1 char name of each filter

  // define arrays used for light curve interpolation
  double x[MXEP_SNOOPY];    // Trest
  double y[MXEP_SNOOPY];    // dm15
  double z[MXEP_SNOOPY];    // flux
  double wz[MXEP_SNOOPY];   // weight

} SNOOPY_TEMPLATE[MXFILT_SNOOPY] ;  // defined by R.K


struct SNOOPY_MODELINFO {
  double TREST_RANGE[2];
  double DM15_RANGE[2];
  double PEAKMAG_REST[MXFILT_SNOOPY];
  double DM15SLOPE[MXFILT_SNOOPY];
  char   GRIDFILE[100]; // snana-generated GRID file
  int    OPT_GRIDFILE;  // >0 => use GRID file for speed

  char ERRORFILE[100];  // interpolate errors from this grid

  // define exclude ranges
  int    NRANGE_EXCLUDE ;
  int    EXCLUDE_IFILT[MXEXCLUDE_SNOOPY] ;
  double EXCLUDE_DM15[MXEXCLUDE_SNOOPY][2];
  double EXCLUDE_TREST[MXEXCLUDE_SNOOPY][2];
} SNOOPY_MODELINFO ;



struct SNOOPY_MODELERR {
  int      SIZE ;
  double **temp_MAP ;
  double  *temp_ERR ;
} SNOOPY_MODELERR[MXFILT_SNOOPY];  
struct GRIDMAP SNOOPY_MODELERR_GRIDMAP[MXFILT_SNOOPY] ;


int  OPT_RELFLUX_SNOOPY ;



// change these later since there are local variable with same name *#$*($#$*(
double sigx0   ;
double sigy0   ;
double xscale  ;
double sigxmax ;

// =============================================
//   function prototypes
// =============================================

int init_genmag_snoopy( char *version, int optmask, double *snoopyArgs, 
			char *filtlist);

int genmag_snoopy(int ifilt, double dm15, int nobs, double *rest_dates, 
		  double *rest_mags, double *rest_magerrs );

void parse_snoopy_modelinfo(char *modelPath);
void read_snoopy_errors(char *modelPath);
int  load_data_snoopy(char *);

void init_spline_snoopy(int ifilt);
void gridinterp_snoopy(int ifilt, double dm15, 
		       int nobs, double *rest_dates, 
		       double *relFlux, double *relFlux_err );

double Textrap_snoopy(int ifilt, double dm15, double Trest );

double snoopy_template_fun(double t, void *p) ;

int dm15temp(int, double, double *, int, double *, double *, 
	     int, double, double, double, double);

void dm15temp_errors(int ifilt, double dm15, double *rest_dates, 
			  int nobs, double *relFlux_err);


int FILTINDX_SNOOPY(char *cfilt);

void MADABORT(char *msg1, char *msg2);

void get_LAMRANGE_snoopy(double *lammin, double *lammax);

// END

