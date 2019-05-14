
#define MXFIELDGROUP_FLUXERRMAP 10
#define MXMAP_FLUXERRMAP        50
#define MXVAR_FLUXERRMAP        10
#define MXROW_FLUXERRMAP        10000

#define IPAR_FLUXERRMAP_MJD     0
#define IPAR_FLUXERRMAP_PSF     1  // FWHM, arcsec
#define IPAR_FLUXERRMAP_SKYSIG  2  // ADU/pixel
#define IPAR_FLUXERRMAP_ZP      3  // ADU
#define IPAR_FLUXERRMAP_LOGSNR  4
#define IPAR_FLUXERRMAP_SBMAG   5  // FLUXCAL per arcSec^2
#define IPAR_FLUXERRMAP_GALMAG  6  // 
#define IPAR_FLUXERRMAP_SNSEP   7  // SN-host sep, arcsec.
#define NPAR_FLUXERRMAP_REQUIRE 8  // number passed to get_FLUXERRMAP

#define IPAR_FLUXERRMAP_ERRSCALE   8  // ERRSCALE option
#define IPAR_FLUXERRMAP_ERRADD     9  // ERRADD   option
#define MXPAR_FLUXERRMAP          10
char VARNAMES_FLUXERRMAP[MXPAR_FLUXERRMAP][20];

#define MASK_APPLY_SIM_FLUXERRMAP  1
#define MASK_APPLY_DATA_FLUXERRMAP 2
#define MASK_DUMP_FLUXERRMAP       256


char FILENAME_FLUXERRMAP[MXPATHLEN];

int  NMAP_FLUXERRMODEL ;
struct {
  char NAME[40];  
  char FIELDLIST_ORIG[40]; // field list of pointer to DEFINE_FIELDGROUP
  char FIELDLIST[40];
  char BANDLIST[MXFILTINDX];
  int  NVAR ; // NDIM + NFUN
  char VARNAMES[MXVAR_FLUXERRMAP][20];
  // xxx  char VARLIST[200] ; // comma-separate list (for screen dump)
  int  IVARLIST[MXVAR_FLUXERRMAP] ; // list of IVAR from full list
  int  MASK_APPLY ;   // bit0 for sim, bit1 for data
  int  INDEX_SPARSE; 
  struct GRIDMAP  MAP ;
} FLUXERRMAP[MXMAP_FLUXERRMAP] ;
  

struct {
  // store info after DEFINE_FIELDGROUP key
  int  NDEFINE ;
  char NAME[MXFIELDGROUP_FLUXERRMAP][40];
  char FIELDLIST[MXFIELDGROUP_FLUXERRMAP][80];
} FLUXERR_FIELDGROUP ;


int NINDEX_SPARSE_FLUXERRMAP ;

double **TMP_ROWDATA_FLUXERRMAP ;

struct {
  int FLAG;
  int NOBS_TOT[MXMAP_FLUXERRMAP] ;
  //  int NOBS_EXTRAP[MXMAP_FLUXERRMAP] ;
  int NOBS_EXTRAP_HI[MXMAP_FLUXERRMAP] ;
  int NOBS_EXTRAP_LO[MXMAP_FLUXERRMAP] ;
} FLUXERRMAP_EXTRAP;

// xxx mark delete int EXTRAPFLAG_FLUXERRMAP;  // 1==> edge extrapolation
int NROW_DUMP_FLUXERRMAP;   

// ======== functions ==============
void  INIT_FLUXERRMODEL(int optmask, char *fileName, 
			char *mapList_ignore_dataErr);
void  init_fluxerrmodel__(int *optmask, char *fileName, 
			  char *mapList_ignore_dataErr);

void  DUMP_FLUXERRMAP(int imap);

void  END_FLUXERRMODEL(void);
void  end_fluxerrmodel__(void);

int   index_sparse_FLUXERRMAP(int NMAP, char *MAPNAME);
int   IVARLIST_FLUXERRMAP(char *varName) ;
void  parse_IGNORE_FLUXERRMAP(char *MAPLIST_IGNORE_DATAERR) ;
void  printSummary_FLUXERRMAP(void);
// xxx delete void  malloc_ROWDATA_FLUXERRMAP(int OPT, int NVAR);

void  get_FLUXERRMODEL(int OPT, double FLUXERR_IN, char *BAND, char *FIELD, 
		       int NPAR, double *PARLIST, 
		       double *FLUXERR_GEN, double *FLUXERR_DATA );
void  get_fluxerrmodel__(int *OPT, double *FLUXERR_IN, char *BAND, char *FIELD, 
			 int *NPAR, double *PARLIST, 
			 double *FLUXERR_GEN, double *FLUXERR_DATA );

void load_parList_FLUXERRMAP(int imap, double *PARLIST, double *parList) ;

double apply_FLUXERRMODEL(int imap, double errModelVal, double FLUXERR_IN);

