// genmag_SIMSED.h


// define OPTMASK bits for generating params
#define OPTMASK_GEN_SIMSED_PARAM    1   // continuous interp
#define OPTMASK_GEN_SIMSED_param    2   // baggage parameter
#define OPTMASK_GEN_SIMSED_GRIDONLY 4   // random gen snapped to GRID only
#define OPTMASK_GEN_SIMSED_WGT  8  // select based on WGT column or external WGTMAP
#define OPTMASK_GEN_SIMSED_SEQ  16 // for "GRIDONLY: SEQUENTIAL" option

// define OPTMASK for init_genmag_SIMSED (LSB=0)
#define OPTMASK_INIT_SIMSED_BINARY    1  // make binary file(s) if not there
#define OPTMASK_INIT_SIMSED_BINARY1   2  // force creation of SED.BINARY
#define OPTMASK_INIT_SIMSED_BINARY2   4  // force create flux-table binary
#define OPTMASK_INIT_SIMSED_TESTMODE  64 // used by SIMSED_check program
#define OPTMASK_INIT_SIMSED_BATCH    128 // batch mode -> abort on stale binary


#define SIMSED_INFO_FILENAME    "SED.INFO" 
#define SIMSED_BINARY_FILENAME  "SED.BINARY" 
char SIMSED_INFO_FILENAME_FULL[MXPATHLEN] ;

// useful numbers

#define INTERP_SIMSED_MAX_DIM    8
#define INTERP_SIMSED_MAX_BAGGAGE_PARS 20
#define INTERP_SIMSED_INVALID_CORNER -1
#define INTERP_SIMSED_START_DIFF  1.0E8
#define INTERP_SIMSED_DELTA       1.0E-8

#define BINARYFLAG_KCORFILENAME 1  // 1 => read/write/check kcor filename


//#define WRVERSION_SIMSED_BINARY  2  // July 30 2017:
//#define WRVERSION_SIMSED_BINARY  3  // Dec 14 2021
#define WRVERSION_SIMSED_BINARY  4  // June 12 2022
int     IVERSION_SIMSED_BINARY ;     // actual version

#define LOGZBIN_SIMSED_DEFAULT 0.02

double Trange_SIMSED[2] ; // used for rd_sedflux
double Lrange_SIMSED[2] ;

int ISIMSED_SELECT;   // either SEQUENTIAL or WGT option

bool ISBATCH_SIMSED;   // T => running in batch mode
bool ISWGTMAP_SIMSED;  // True = user-provided SIMSED_WGTMAP file

/**********************************************
  Init Information
***********************************************/

char SIMSED_PATHMODEL[MXPATHLEN];
char SIMSED_KCORFILE[MXPATHLEN];
// xxx mark delete Dec 2021 char SIMSED_PATHBINARY[MXPATHLEN]; 

struct { 
  bool USE ;   // flat that binary table is used (read or write mode)
  char PATH[MXPATHLEN];  // path where binaries are read or written

  // read and write flag for both binary tables
  bool WRFLAG_SED;  // flag to write binary for SEDs
  bool RDFLAG_SED;  // flag to read existing binary
  bool WRFLAG_FLUX; // flag to write binary for flux-integral table
  bool RDFLAG_FLUX; // flag to read

  // force-create options 
  bool FORCE_CREATE_SED;       // set if SIMSED_USE_BINARY += 2
  bool FORCE_CREATE_FLUX;      // set if SIMSED_USE_BINARY += 4
 
} SIMSED_BINARY_INFO ;

/**********************************************
   Function Declarations
**********************************************/

int  init_genmag_SIMSED(char *version, char *PATH_BINARY, char *SURVEY, char *kcorFile, char *WGTMAP_FILE, int OPTMASK );

int read_SIMSED_INFO(char *PATHMODEL );
int read_simsed_info__(char *PATHMODEL );
int count_SIMSED_INFO(char *PATHMODEL);

void set_SIMSED_MXDAY(char *PATHMODEL, FILE *fpbin, 
		      bool RDFLAG_BINARY, bool WRFLAG_BINARY );
void set_SIMSED_LOGZBIN(void);

void set_SIMSED_WGT_SUM(char *WGTMAP_FILE);
int  pick_SIMSED_BY_WGT(void);
void dump_SIMSED_WGT(int NSED, char *msg);

void dump_SIMSED_INFO(void);

void open_SEDBINARY(char *fileName, bool force_create, 
		    FILE **fpbin, bool *RDFLAG, bool *WRFLAG);
void open_TABBINARY(char *fileName, bool force_create, 
		    FILE **fpbin, bool *RDFLAG, bool *WRFLAG);

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

int IS_INDEX_SIMSED(char *parName) ;

int IS_WGT_SIMSED(char *parName) ;

// END
