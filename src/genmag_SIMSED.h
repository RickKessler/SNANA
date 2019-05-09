// genmag_SIMSED.h


#define OPTMASK_SIMSED_PARAM    1 // continuous interp
#define OPTMASK_SIMSED_param    2 // baggage parameter
#define OPTMASK_SIMSED_GRIDONLY 4 // random gen snapped to GRID only

#define INFO_SIMSED_FILENAME  "SED.INFO" 
char INFO_SIMSED_FILENAME_FULL[MXPATHLEN] ;

// useful numbers

#define INTERP_SIMSED_MAX_DIM    8
#define INTERP_SIMSED_MAX_BAGGAGE_PARS 20
#define INTERP_SIMSED_INVALID_CORNER -1
#define INTERP_SIMSED_START_DIFF  1.0E8
#define INTERP_SIMSED_DELTA       1.0E-8

#define BINARYFLAG_KCORFILENAME 1  // 1 => read/write/check kcor filename

// define OPTMASK bits for init_genmag_SIMSED (LSB=0)
#define OPTMASK_SIMSED_BINARY    1  // --> make binary file
#define OPTMASK_SIMSED_TESTMODE  64 // used by SIMSED_check program

#define WRVERSION_SIMSED_BINARY  2  // July 30 2017:
int     IVERSION_SIMSED_BINARY ;     // actual version

#define LOGZBIN_SIMSED_DEFAULT 0.02

double Trange_SIMSED[2] ; // used for rd_sedflux
double Lrange_SIMSED[2] ;

int ISIMSED_SEQUENTIAL ;

/**********************************************
  Init Information
***********************************************/

char SIMSED_PATHMODEL[MXPATHLEN];
char SIMSED_KCORFILE[MXPATHLEN];
char SIMSED_PATHBINARY[MXPATHLEN]; // July 30 2017

/**********************************************
   Function Declarations
**********************************************/

int  init_genmag_SIMSED(char *version, char *PATH_BINARY, char *SURVEY,
			char *kcorFile, int OPTMASK );

int read_SIMSED_INFO(char *PATHMODEL );
int read_simsed_info__(char *PATHMODEL );
int count_SIMSED_INFO(char *PATHMODEL);

void set_SIMSED_MXDAY(char *PATHMODEL, FILE *fpbin, 
		      int RDFLAG_BINARY, int WRFLAG_BINARY );
void set_SIMSED_LOGZBIN(void);

void dump_SIMSED_INFO(void);

void open_SEDBINARY(char *fileName, FILE **fpbin, int *RDFLAG, int *WRFLAG);
void open_TABBINARY(char *fileName, FILE **fpbin, int *RDFLAG, int *WRFLAG);

void read_SIMSED_TABBINARY(FILE *fp, char *binFile);

void genmag_SIMSED(int OPTMASK, int ifilt, double x0, 
		   int NLUMIPAR, int *iflagpar, int *iparmap, double *lumipar,
		   double RV_host, double AV_host, double mwebv, double z, 
		   int nobs, double *Tobs_list, 
		   double *magobs_list, double *magerr_list, int *index_sed );

double interp_flux_SIMSED(int *iflagpar, int *iparmap, double *lumipar, 
			  int ifilt_obs, double z, double Trest );

double nextgrid_flux_SIMSED(int *iflagpar, int *iparmap, double *lumipar, 
			    int ifilt_obs, double z, double Trest );

double interp1D_flux_SIMSED(int *iflagpar, double *lumipar, int ifilt_obs,
			    double z, double Trest );

void checkBinary_SIMSED(char *binaryFile); // abort if earlier than SED.INFO

void read_SIMSED_flux(char *sedFile, char *sedComment) ;

// END
