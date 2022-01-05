// genmag_snoopy.h

#define MXFILT_SNOOPY    12  // max number of rest-frame filters to define
#define IFILT_SNOOPY_MIN  0  //  xxx1
#define IFILT_SNOOPY_MAX  8  //  xxx9

#define MINLAM_SNOOPY   2900.0  // blue edge of u
#define MAXLAM_SNOOPY  20000.0  // H band

char SNOOPYNAME[20] ;
char FILTLIST_SNOOPY[MXFILT_SNOOPY] ; 
char SNOOPY_MODELPATH[200];

SNGRID_DEF SNGRID_SNOOPY  ;

// =============================================
//   global snoopy arrays to contain the data 
// =============================================


struct SNOOPY_MODELINFO {

  char   GRIDFILE[100]; // snana-generated GRID file
  int    OPT_GRIDFILE;  // >0 => use GRID file for speed

} SNOOPY_MODELINFO ;


// =============================================
//   function prototypes
// =============================================

int init_genmag_snoopy( char *version, int optmask, char *filtlist);

int genmag_snoopy(int ifilt, double stretch, int nobs, double *Trest, 
		  double *mag_list, double *magerr_list );

void parse_snoopy_modelinfo(char *modelPath);

void gridinterp_snoopy(int ifilt, double dm15, 
		       int nobs, double *Trest_list, 
		       double *mag_list, double *magerr_list );

void get_LAMRANGE_snoopy(double *lammin, double *lammax);

// END

