/*************************************************
  May 2011, R.Kessler


  Parameters and declarations for snfitsio.c .

    HISTORY

  Feb 15, 2013: 
      MXPAR_SNFITSIO  -> 200  (was 100) to handle 56 JPAS filters

  May 16, 2013
      MXPAR_SNFITSIO  -> 400  (was 200) to handle 56 JPAS filters
      with HOSTLIB

  Feb 2014:
   - define  SNFITSIO_CODE_IVERSION
   - SNFITSIO_PATH    -> SNFITSIO_DATA_PATH
   - SNFITSIO_VERSION -> SNFITSIO_PHOT_VERSION


  Apr 30 2014: MXFILE_SNFITSIO -> 50 (was 20)

  Aug 6 2014: SNFITSIO_SIMFLAG -> SNFITSIO_SIMFLAG_SNANA 
              and add SNFITSIO_SIMFLAG_MAGOBS

  Aug 2 2016: 
   + define SNFITSIO_SIMFLAG_SPECTROGRAPH
   + define MXTYPE_SNFITSIO and all dimention [2] -> [MXTYPE_SNFITSIO]

  July 25 2017: MXFILE_SNFITSIO -> 100 (was 50) 
  Mar  23 2019: MXFILE_SNFITSIO -> 300 (was 100) 

**************************************************/

// ==================================
// global variables


#define ITYPE_SNFITSIO_HEAD      0
#define ITYPE_SNFITSIO_PHOT      1
#define ITYPE_SNFITSIO_SPEC      2  // SPEC summary & fluxes tables
#define ITYPE_SNFITSIO_SPECTMP   3  // temp flux table
#define MXTYPE_SNFITSIO          4

#define MXFILE_SNFITSIO     300  // max number of fits-files to read
#define MXPAR_SNFITSIO      400  // max number of header variables

#define SNFITSIO_EOE_MARKER  -777.0  // from era of 9-track tapes


fitsfile  *fp_snfitsFile[MXTYPE_SNFITSIO] ;
#define  snfitsType  (char*[MXTYPE_SNFITSIO]) { "HEAD", "PHOT", "SPEC", "SPECTMP"  }

// name of FITS file without/with path
char  snfitsFile[MXFILE_SNFITSIO][MXTYPE_SNFITSIO][MXPATHLEN];   
char  snfitsFile_plusPath[MXFILE_SNFITSIO][MXTYPE_SNFITSIO][MXPATHLEN]; 

int   SNFITSIO_CODE_IVERSION  ; // internal: for back-compatibility
char  SNFITSIO_SNANA_VERSION[12]   ; // e.g, v10_79b (Oct 2020)

char  SNFITSIO_DATA_PATH[MXPATHLEN];
char  SNFITSIO_PHOT_VERSION[MXPATHLEN];
char  SNFITSIO_LISTFILE[MXPATHLEN];   // for read-back only
char  SNFITSIO_READMEFILE[MXPATHLEN]; // for read-back only


int NSNLC_SNFITSIO_TOT ;    // total number of SNe over all files.
int NSNLC_SNFITSIO[MXFILE_SNFITSIO];     // Number of SNe per file
int NSNLC_SNFITSIO_SUM[MXFILE_SNFITSIO]; // cumulative number


int NFILE_SNFITSIO ;     // number of fits files
int IFILE_SNFITSIO ;     // current fits-file index
int ISNFIRST_SNFITSIO ;  // first ISN in file
int NPAR_SNFITSIO[MXTYPE_SNFITSIO] ;  // Npar (i.e., columns) for each fits-file type
int MXOBS_SNFITSIO ;     // max NOBS among SNe (to allocate memory)

int MALLOC_LEN_SNFITSIO[MXTYPE_SNFITSIO] ; // malloc length per file type

bool  SNFITSIO_DATAFLAG ;       // true -> real data (not sim, not fakes)
bool  SNFITSIO_SIMFLAG_SNANA ;  // SNANA sim
bool  SNFITSIO_SIMFLAG_MAGOBS ; // data-like with SIM_MAGOBS
bool  SNFITSIO_SIMFLAG_SPECTROGRAPH ;  // simulated spectra (Aug 2016)
bool  SNFITSIO_SIMFLAG_SNRMON      ;   // SNR(MAGMONITOR)
bool  SNFITSIO_SIMFLAG_MODELPAR    ;   // model params for SIMSED, LCLIB
bool  SNFITSIO_SIMFLAG_NBR_LIST    ;   // HOSTLIB has NBR_LIST (Feb 2020)
bool  SNFITSIO_COMPACT_FLAG ;    // Jan 2018

char SNFITSIO_VARNAME_SNRMON[40] ;    // includes int mag value in name

int  SNFITSIO_SUBSURVEY_FLAG ;  // indicates subSurvey column
int  SNFITSIO_NSUBSAMPLE_MARK ; // indicates how many marked sub-samples
char SNFITSIO_DATATYPE[20] ;


struct TABLEDEF {
  // name of each header paramater (SNID, REDSHIFT, etc ...)
  char  name[MXPAR_SNFITSIO][40] ;
  char *ptrName[MXPAR_SNFITSIO] ;

  // cast of each header parameter (float, double, int ...)
  char  form[MXPAR_SNFITSIO][4] ;
  char *ptrForm[MXPAR_SNFITSIO] ;

  // these are all set to blank
  char *ptrUnit[MXPAR_SNFITSIO] ;

  int iform[MXPAR_SNFITSIO] ;

} SNFITSIO_TABLEDEF[MXTYPE_SNFITSIO] ;  // index is itype


struct WR_TABLEVAL_DEF {
  // temp values to fill Header table
  char          *value_A  ;  // ascii/text
  float          value_1E ;  // 4-byte float
  double         value_1D ;  // 8-byte double
  int            value_1J ;  // 4-byte signed int
  short int      value_1I ;  // 2-byte signed int
  unsigned short value_1U ;  // 2-byte unsigned int
  unsigned int   value_1V ;  // 4-byte unsigned int
  long long      value_1K ;  // 8 bytte long long int (May 2013)
  
  int NROW ; // increment numver of rows written

  // index used to speed search for header-param column during update
  // Note that the index is local and NOT NWR_SNFITSIO
  int COLNUM_LOOKUP[MXPAR_SNFITSIO] ;

} WR_SNFITSIO_TABLEVAL[MXTYPE_SNFITSIO] ;  // index is itype


#define IFORM_A   1
#define IFORM_1J  2
#define IFORM_1I  3
#define IFORM_1E  4
#define IFORM_1D  5
#define IFORM_1K  6
#define MXFORM_SNFITSIO 10

#define NULL_A      "NULL"
#define NULL_1J     (int)-9
#define NULL_1I     (short)-9
#define NULL_1E     (float)-9.0
#define NULL_1D     (double)-9.0
#define NULL_1K     (long long)-9


struct RD_TABLEVAL_DEF {
  int    NPAR[MXFORM_SNFITSIO] ;
  int    IPAR[MXFORM_SNFITSIO][MXPAR_SNFITSIO] ;    // absolute IPAR/column
  int    IPARINV[MXFORM_SNFITSIO][MXPAR_SNFITSIO] ; // inverse of above
} RD_SNFITSIO_TABLEVAL[MXTYPE_SNFITSIO] ;  // index is itype


// index is type (0=HEAD or 1=PHOT)
char     ***RD_SNFITSIO_TABLEVAL_A[MXTYPE_SNFITSIO] ; 
int       **RD_SNFITSIO_TABLEVAL_1J[MXTYPE_SNFITSIO] ;
short     **RD_SNFITSIO_TABLEVAL_1I[MXTYPE_SNFITSIO] ;
float     **RD_SNFITSIO_TABLEVAL_1E[MXTYPE_SNFITSIO] ;
double    **RD_SNFITSIO_TABLEVAL_1D[MXTYPE_SNFITSIO] ;
long long **RD_SNFITSIO_TABLEVAL_1K[MXTYPE_SNFITSIO] ;

// define absolute head-par indices for required elements
int IPAR_SNFITSIO_SNID ;
int IPAR_SNFITSIO_FAKE ;
int IPAR_SNFITSIO_NOBS ;
int IPAR_SNFITSIO_PTROBS_MIN ;
int IPAR_SNFITSIO_PTROBS_MAX ;

#define  stringBlank " " ;

// maks of epochs to store in RD_SNFITSIO_PARVAL function
// Default is all epochs unless SET_RDMASK_SNFITSIO is called
// with some epochs masked out.
int NEP_RDMASK_SNFITSIO_PARVAL ;
int RDMASK_SNFITSIO_PARVAL[MXEPOCH] ;


// define RDSPEC structures to read back SPEC info (April 2019) 
struct {
  int     NLAMBIN;
  double *LAMMIN_LIST, *LAMMAX_LIST ;
} RDSPEC_SNFITSIO_LAMINDEX ;

struct {  
  int     NROW  ;
  char   **SNID ;
  int    *NLAMBIN ;
  double *MJD;
  float  *TEXPOSE;
  int    *PTRSPEC_MIN, *PTRSPEC_MAX ;
} RDSPEC_SNFITSIO_HEADER ;

// ================================================
// ============= FUNCTION PROTOTYPES ==============
// ================================================


// note that WR_SNFITSIO_XXX are called externally;
// wr_snfitsio_xxx are called internally.

void WR_SNFITSIO_INIT(char *path, char *version, char *prefix,
		      int simFlag, int Nsubsample_mark, char *headFile);

int  is_fits(char *file);
void wr_snfitsio_create(int itype);
void wr_snfitsio_init_head(void);
void wr_snfitsio_init_phot(void);
void wr_snfitsio_init_spec(void);
void wr_snfitsio_addCol(char *tform, char *name, int  itype);

void WR_SNFITSIO_UPDATE(void);
void wr_snfitsio_update_head(void);
void wr_snfitsio_update_phot(int ep);
void wr_snfitsio_update_spec(int imjd);
void wr_snfitsio_fillTable(int *COLNUM, char *parName, int itype );

void WR_SNFITSIO_END(void);

void snfitsio_close(int ifile, int itype);
void snfitsio_errorCheck(char *comment, int status);
int  IPAR_SNFITSIO(int OPT, char *parName, int itype );
int  IPARFORM_SNFITSIO(int OPT, int iform, char *parName, int itype);

// Now the readback routines
int   RD_SNFITSIO_INIT(int MSKOPT, char *PATH, char *version);
int   RD_SNFITSIO_GLOBAL(char *parName, char *parString);
void  RD_SNFITSIO_CLOSE(char *version);
void  GET_SNFITSIO_INFO(char *VERSION, char *FILENAME_HEAD, 
			char *FILENAME_PHOT, int *IFILE );

int   rd_snfitsio_list(void);
void  rd_snfitsio_open(int ifile, int photflag_open, int vbose ); 
void  rd_snfitsio_file(int ifile);          // open and read everything
void  rd_snfitsio_simkeys(void);            // read optional SIMSED pars
void  rd_snfitsio_private(void);            // read optional PRIVATE vars
void  rd_snfitsio_free(int ifile, int itype);
void  rd_snfitsio_malloc(int ifile, int itype, int LEN);
void  rd_snfitsio_head(int ifile);
void  rd_snfitsio_tblpar(int ifile, int itype);
void  rd_snfitsio_tblcol(int itype, int icol, int firstRow, int lastRow);

void  rd_snfitsio_specFile(int ifile); 
void  rd_snfitsio_mallocSpec(int opt);

int RD_SNFITSIO_PARVAL(int isn, char *parName, 
		      double *parLIST, char *parString, int *iptr);

int RD_SNFITSIO_STR(int isn, char *parName, char *parString, int *ipar);
int RD_SNFITSIO_INT(int isn, char *parName, int    *parList, int *ipar);
int RD_SNFITSIO_SHT(int isn, char *parName, short int *parList, int *ipar);
int RD_SNFITSIO_FLT(int isn, char *parName, float  *parList, int *ipar);
int RD_SNFITSIO_DBL(int isn, char *parName, double *parList, int *ipar);
void RD_SNFITSIO_SPECROWS(char *SNID, int *ROWMIN, int *ROWMAX);
void RD_SNFITSIO_SPECDATA(int irow, double *METADATA, int *NLAMBIN,
			  double *LAMMIN, double *LAMMAX, 
			  double *FLAM, double *FLAMERR);

void  check_required_headkeys(void) ;
int   formIndex_snfitsio(char *form) ;

void SET_RDMASK_SNFITSIO(int N, int *mask) ;

// ------- mangled RD functions for snana/fortran -----------

int  rd_snfitsio_init__(int *MSKOPT, char *PATH,  char *version) ;
int  rd_snfitsio_global__(char *parName, char *parString) ;
void rd_snfitsio_close__(char *version) ;

void get_snfitsio_info__(char *VERSION, char *FILENAME_HEAD, 
			char *FILENAME_PHOT, int *IFILE );

int rd_snfitsio_parval__(int *isn, char *parName, 
			 double *parLIST, char *parString, int *iptr);

int rd_snfitsio_str__(int *isn, char *parName, char *parString, int *iptr) ;
int rd_snfitsio_int__(int *isn, char *parName, int *parLIST, int *iptr) ;
int rd_snfitsio_sht__(int *isn, char *parName, short int *parLIST, int *iptr);
int rd_snfitsio_flt__(int *isn,  char *parName, float *parLIST, int *iptr) ;
int rd_snfitsio_dbl__(int *isn,  char *parName, double *parLIST, int *iptr) ;

void set_rdmask_snfitsio__(int *N, int *mask) ;
void rd_snfitsio_specrows__(char *SNID, int *ROWMIN, int *ROWMAX );
void rd_snfitsio_specdata(int *irow, double *METADATA, int *NLAMBIN,
			  double *LAMMIN, double *LAMMAX, 
			  double *FLAM, double *FLAMERR);

// mangled write funs

void wr_snfitsio_update__(void) ;

void wr_snfitsio_init__(char *path, char *version, char *prefix, 
			int *simFlag, int *Nsubsample_mark, char *headFile ) ;

void wr_snfitsio_end__(void) ;


// ========= END ============

