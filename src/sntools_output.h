/*******************************************
  Created Feb 2013 by R.Kessler
  Generic functions for filling tables and packing lightcurves.
  

 Mar 28 2016: 
   + defined int NVAR_ADDCOL_TOT
   + MXVAR_TABLE -> 240 (was 200)

 Apr 13 2016: MXVAR_TABLE -> 400 (for Hounsell)

 Sep 12 2016: allow up to 2 pointers per variable to read a column 
              in 2 different arrays.

 Jan 6, 2017: 
    + MXCHAR_VARLIST->2000 (was 1000)
    + define MXFILE_AUTOSTORE & NFILE_AUTOSTORE 
   
 Feb 7 2017:
    + MXCHAR_FILENAME -> 240  (was 200)
 
 Apr 11 2017: add MXCHAR_MODELNAME

 May 11 2017: declare IVAR_READTABLE_POINTER

 Apr 4 2019: preproc flags HBOOK,ROOT,TEXT -> USE_[HBOOK,ROOT,TEXT]
 
 Dec 11 2019: replace lots of [40] with [MXCHAR_VARNAME]

 May 02 2020: add spectra format for MARZ (FITS with particular extensions)
 May 30 2020: MXSPEC_SPECPAK -> MXSPECTRA from sndata.h

*******************************************/


// define flags for software packages
#define USE_HBOOKxxx             
#define USE_ROOT     
#define USE_TEXT  // always leave this on; same logic as for HBOOK,ROOT, ...
#define USE_MARZ  // always leave this on

// ---------------------------------------
// flags to identify TABLEFILE_TYPE 

#define IFILETYPE_NULL   0   // undefined file type
#define IFILETYPE_HBOOK  1
#define IFILETYPE_ROOT   2
#define IFILETYPE_TEXT   3
#define IFILETYPE_MARZ   4
#define MXTABLEFILETYPE  5

#define MXCHAR_FILENAME  240
#define MXCHAR_VARLIST   2000  
#define MXCHAR_VARNAME   60
#define MXCHAR_CCID      20  // should be same as MXCHAR_CCID in snana.car
#define MXCHAR_MODELNAME 32  // max length of model name (e.g., SALT2.Guy10)

#define MXVAR_TABLE      400  // max number of variables in 1 table
#define MXLINE_TABLECOMMENT  20   // max number of user-passed comments

#define OPENFLAG_NULL   0
#define OPENFLAG_NEW    1
#define OPENFLAG_READ   2
#define MXOPENFLAG      3

// ------ define standard table name/id used in snana programs ------
#define TABLENAME_FITRES  "FITRES"
#define TABLENAME_SNANA   "SNANA"
#define TABLEID_FITRES     7788
#define TABLEID_SNANA      7100
#define TABLEID_MARZ       8100

char STRING_TABLEFILE_TYPE[MXTABLEFILETYPE][12] ;
  // =  { "NULL", "HBOOK", "ROOT", "TEXT", "MARZ" } ;

char STRING_TABLEFILE_OPENFLAG[MXOPENFLAG][12] ;
 //  = { "NULL", "NEW", "READ" } ;

char STRING_IDTABLE_SNANA[MXTABLEFILETYPE][12] ;
char STRING_IDTABLE_FITRES[MXTABLEFILETYPE][12] ;
char STRING_IDTABLE_OUTLIER[MXTABLEFILETYPE][12] ; // Mar 2021

int  CALLED_TABLEFILE_INIT ;  // flag to ensure call to TABLEFILE_INIT

// define TABLE-FILE name for each open-flag and file-type.
char NAME_TABLEFILE[MXOPENFLAG][MXTABLEFILETYPE][MXCHAR_FILENAME] ;

// define USE-flag for each open-flag and file-type to make
// sure that we don't open another.
int  USE_TABLEFILE[MXOPENFLAG][MXTABLEFILETYPE]  ; 
int  NOPEN_TABLEFILE ; // number of table files (should be >0)

// logical flag for first call to root,hbook...
// e.g., to call hlimit only once for hbook option.
int  FIRSTCALL_TABLEFILE_OPEN[MXTABLEFILETYPE];


int  NLINE_TABLECOMMENT ;
char LINE_TABLECOMMENT[MXLINE_TABLECOMMENT][MXCHAR_FILENAME];

// -------------------------------------
// define a few things from sntools.h so that we don't have to
// include all of sntools.h

#define SEV_FATAL  4      // must match value in sntools.h, for errmsg call
#define MXSPEC_SPECPAK MXSPECTRA

#define ICAST_L   16  // long long int (64 bits)
#define ICAST_D    8  // double
#define ICAST_F    4  // float
#define ICAST_S    3  // 2-byte short integer 
#define ICAST_I    2  // 4-byte integer 
#define ICAST_C    1  // char
#define MXCAST    20
char    CCAST_TABLEVAR[MXCAST] ;  // see TABLEVAR_INIT().

char  FILTERSTRING[MXFILTINDX] ;   // list of all filters


// -------------------
// define ADDCOL struct for creating table
#define MXVAR_ADDCOL 50 // max NVAR per call to SNTABLE_ADDCOL
int  NVAR_ADDCOL_TOT ;         // Mar 2016: sum of NVAR
typedef struct {
  int  NVAR ;
  char VARLIST_ORIG[MXCHAR_VARLIST] ;  // original VARLIST before parsing
  char VARNAME[MXVAR_ADDCOL][MXCHAR_VARNAME] ;   // each VARNAME
  int  ICAST[MXVAR_ADDCOL] ;    // integer cast
  char CCAST[MXVAR_ADDCOL][4] ;    // char cast: 'C', 'I', 'F', 'D'
  int  VECTOR_FLAG[MXVAR_ADDCOL] ; // 0=scaler, 1 = vector size, 2=vector
  int  ISIZE[MXVAR_ADDCOL];        // char size, or vector size
} SNTABLE_ADDCOL_VARDEF ; 

char ADDCOL_VARLIST_LAST[MXCHAR_VARLIST];


// Oct 14 2014: define struct for storing pointers to read table
struct READTABLE_POINTERS {

  int   IFILETYPE ;      // store file type (HBOOK,ROOT, TEXT ...)
  char  TABLENAME[100];  // store name of table to read

  // table-header info
  int    NVAR_TOT ; // total number of  variables
  int    NVAR_READ ; // total number to read
  char   VARNAME[MXVAR_TABLE][MXCHAR_VARNAME];   // index is NVARTOT
  int    PTRINDEX[MXVAR_TABLE];     // index is NVAR_READ

  // store pointer for each variable to read so that user can fill any cast
  int     NPTR[MXVAR_TABLE];          // number of pointers per var (9/2016)
  double     *PTRVAL_D[2][MXVAR_TABLE];      // index is NVARTOT
  float      *PTRVAL_F[2][MXVAR_TABLE];      // index is NVARTOT
  int        *PTRVAL_I[2][MXVAR_TABLE];      
  short int  *PTRVAL_S[2][MXVAR_TABLE];      
  long long int *PTRVAL_L[2][MXVAR_TABLE]; 
  char   **PTRVAL_C[2][MXVAR_TABLE];  // Aug 2013
  
  // ... or store file pointer for ascii dump
  FILE *FP_DUMP ;
  char LINEKEY_DUMP[40]; // key for reach dumped row  (e.g., 'SN: ')
  char SEPKEY_DUMP[2];   // var-separator string

  // store cast for each variable
  int    ICAST_READ[MXVAR_TABLE];   // cast read from file 
  int    ICAST_STORE[MXVAR_TABLE];  // cast to store after read

  int    MXLEN ;     // max size of PTRVAL_X arrays (for internal check)
  int    NROW;       // number of rows read from table

} READTABLE_POINTERS ;


// -------------------------------------------------
// SNLCPAK global declarations for light curves

#define SNLCPAK_INITDONE  12345
#define SPECPAK_INITDONE  98765
#define MXFIT_PER_SN      10

bool SNLCPAK_USE_HBOOK ; 
bool SNLCPAK_USE_ROOT  ;
bool SNLCPAK_USE_HDF5  ;
bool SNLCPAK_USE_TEXT  ;
int NCALL_SNLCPAK_FILL ; // number of calls to SNLCPAK_FILL

bool SPECPAK_USE_HBOOK ; 
bool SPECPAK_USE_ROOT  ;
bool SPECPAK_USE_HDF5  ;
bool SPECPAK_USE_TEXT  ;
bool SPECPAK_USE_MARZ  ;
int NCALL_SPECPAK_FILL ; // number of calls to SPECPAK_FILL

// define flags for epoch info
#define SNLCPAK_EPFLAG_FLUXDATA    1 // data flux 
#define SNLCPAK_EPFLAG_REJECT      2 // for data only
#define SNLCPAK_EPFLAG_CHI2        3 // for fitted data only
#define SNLCPAK_EPFLAG_FITFUN      4
#define SNLCPAK_EPFLAG_FLUXSIM     5  // simulated flux (May 20 2016)
#define SNLCPAK_EPFLAG_FLUXREST    6  // k-corrected rest-frame fluxes
#define SNLCPAK_EPFLAG_KCOR        7  // k-cor (Jan 2020)
#define SNLCPAK_EPFLAG_AVWARP      8  // AVwarp (Jan 2020)
#define SNLCPAK_EPFLAG_SIMFLUXREST 9  // simulated rest-flux
#define MXFLAG_SNLCPAK_EPOCH      10  

// define flags for filter info (one entry per filter)
#define SNLCPAK_BANDFLAG_NDOF      100
#define SNLCPAK_BANDFLAG_PKFLUX    101
#define SNLCPAK_BANDFLAG_PKMJD     102
#define SNLCPAK_BANDFLAG_CHI2      103
#define MXFLAG_SNLCPAK_FILTER      4

#define MXTEXT_SNLCPAK 5  // max number of displayText strings


// define structure containing information for 1 light curve.
struct  SNLCPAK_OUTPUT {
  int NCALL ;
  int NFIT_PER_SN ;  // number of fits per SN

  // define user-passed global variables set just once from SNLCPAK_INIT
  char SURVEY[60] ;
  char VERSION_PHOTOMETRY[200] ;
  char VERSION_SNANA[40];
  char SURVEY_FILTERS[MXFILTINDX] ;

  char TEXT_FORMAT[20]; // added Sep 2014
  int  OPT_TEXT_FORMAT ;  // integer index for above (-9 -> no text output)

  // set user-passed variables that change with each SN
  char CCID[40];

  int  NTEXT ;         // number of text strings
  int  MXTEXT ;        // to init text in tree
  char DISPLAYTEXT[MXTEXT_SNLCPAK][80] ;  // display strings

  int      NOBS_MAX ;      // max of all NOBS (for malloc)
  int      NOBS[MXFLAG_SNLCPAK_EPOCH] ;  // NOBS for each SNLCPAK_FLAG
  double  *MJD[MXFLAG_SNLCPAK_EPOCH] ;   // pointer to dynam allocate
  double  *TOBS[MXFLAG_SNLCPAK_EPOCH] ;  // pointer to dynam allocate
  double  *EPDATA[MXFLAG_SNLCPAK_EPOCH] ;
  double  *EPDATA_ERR[MXFLAG_SNLCPAK_EPOCH] ;
  int     *IFILTOBS[MXFLAG_SNLCPAK_EPOCH] ;

  double  FILTDATA[MXFLAG_SNLCPAK_FILTER][MXFILTINDX];  
  double  FILTDATA_ERR[MXFLAG_SNLCPAK_FILTER][MXFILTINDX];  

  int  NLCPAK_TOT ;  // expected number of LC to PACK 
  int  NLCPAK ;      // incremental counter
  int  SIMFLAG ;     // =1 for simulated fluxes (May 2016)

  // internally defined variable
  int  INITDONE  ;
  int  NFILTDEF_SURVEY ;
  int  NFILTDEF_PLOT ;             // Number of filters to plot  
  char FILTLIST_PLOT[MXFILTINDX] ; // char-list of plotted filters
  int  IFILT_MIN, IFILT_MAX;       // min,max sparse filter index


  int  USEFILT[MXFLAG_SNLCPAK_EPOCH][MXFILTINDX];   // vs. sparse IFILT indx
  int  NOBS_FILT[MXFLAG_SNLCPAK_EPOCH][MXFILTINDX] ;  // idem per filter
  int  *IFILT[MXFLAG_SNLCPAK_EPOCH] ;               // sparse IFILT

  char  TEXTFILE_PREFIX[200];

} SNLCPAK_OUTPUT ;


// Apr 2019: create struture for multiple spectra per event (SN & HOST)
struct SPECPAK_OUTPUT {
  int  NSPEC;  // per event
  char SURVEY[60] ;
  char VERSION_PHOTOMETRY[200] ;
  int  INITDONE ;

  char TEXT_FORMAT[20]; 
  int  OPT_TEXT_FORMAT ;  // integer index for above (-9 -> no text output)

  // set user-passed variables that change with each SN
  char    CCID[40];
  int     NLAMBIN_TOT;    // sum of lam bins for each event
  int     ID_LIST[MXSPEC_SPECPAK];
  int     NLAMBIN_LIST[MXSPEC_SPECPAK];
  double  MJD_LIST[MXSPEC_SPECPAK];
  double  TOBS_LIST[MXSPEC_SPECPAK]; 
  double  TEXPOSE_LIST[MXSPEC_SPECPAK];

  int    *ID ;
  double *LAMMIN, *LAMMAX, *FLAM, *FLAMERR ;

} SPECPAK_OUTPUT ;



// special struct for outlier option (8.03.2014)
#define  INDX_OUTLIER_NPTFIT  0   // sparse index for NPTFIT in FITRES+RESID
#define  INDX_OUTLIER_NOBS    0   // oidem for SNANA+EPOCHS table
#define  INDX_OUTLIER_IFILT   1   // sparse index for IFILTOBS
#define  INDX_OUTLIER_CHI2    2   // sparse index for CHI2FLUX
#define  NVAR_OUTLIER_DECODE  3

#define  OUTLIER_VARNAME_NPTFIT    "NPTFIT"
#define  OUTLIER_VARNAME_NOBS      "NOBS"
#define  OUTLIER_VARNAME_IFILT     "IFILTOBS"  
#define  OUTLIER_VARNAME_CHI2      "CHI2FLUX" 
#define  VARNAME_CUTFLAG_SNANA     "CUTFLAG_SNANA" 

// read-flags
#define OPT_SNTABLE_READ_forARRAY      1   // read to fill array in memory
#define OPT_SNTABLE_READ_forDUMP       2   // read to dump to ascii file
#define OPT_SNTABLE_READ_forOUTLIERS   3   // read to dump outlier to ascii


struct OUTLIER_INFO {
  int   USEFLAG ;
  int   IVAR[NVAR_OUTLIER_DECODE] ;
  char  VARNAME[NVAR_OUTLIER_DECODE][MXCHAR_VARNAME];

  float CUTWIN_NSIGMA[2];
  float CUTWIN_CHI2FLUX[2];

  // define stats vs. filter index; 0--> all
  int   NEP_TOT[MXFILTINDX];    // all epochs
  int   NEP_SELECT[MXFILTINDX]; // selected outliers

} OUTLIER_INFO ;


#define MXFILE_AUTOSTORE 10   // max files to autoStore (Jan 2017)
int NFILE_AUTOSTORE ;
int NREAD_AUTOSTORE ;
struct SNTABLE_AUTOSTORE {
  int     NVAR ;
  int     IVARMAP[MXVAR_TABLE]; 
  char    VARNAME[MXVAR_TABLE][MXCHAR_VARNAME];
  int     EXIST[MXVAR_TABLE];
  int     ICAST_READ[MXVAR_TABLE] ; // cast read from file
  int     IFILETYPE ;   // indicates ROOT,HBOOK,TEXT, etc ...

  int     NROW ;
  int     *LENCCID; // string len for each CCID (for faster lookup)
  char    **CCID ;
  double  **DVAL ;
  char    ***CVAL ;

} SNTABLE_AUTOSTORE[MXFILE_AUTOSTORE] ;


// define LASTREAD structure to speed up AUTOSTORE lookup
// when CCID is repeated.
struct LASTREAD_AUTOSTORE  {
  int  IFILE, IROW;
  char CCID[MXCHAR_CCID];
} LASTREAD_AUTOSTORE ;

// generic strings for errmsg 
char MSGERR1[200], MSGERR2[200] ;
char SNANA_VERSION[20] ;

// -------------------------------------------------
//                   FUNCTIONS
// -------------------------------------------------

#ifdef __cplusplus
extern"C" {
#endif

  void  get_SNANA_VERSION(char *SNANA_VERSION);  // defined in sntools.c
  void  set_FILTERSTRING(char *FILTERSTRING); // defined in sntools.c
  void  errmsg ( int isev, int iprompt, char *fnam, char *msg1, char *msg2 );
  void  print_banner ( const char *banner ) ;
  void  print_preAbort_banner(char *fnam);
  void  trim_blank_spaces(char *string) ;
  int   strcmp_ignoreCase(char *str1, char *str2) ;
  void  debugexit(char *string);
  void catVarList_with_comma(char *varList, char *addVarName);

  void  checkval_I(char *varname,int nval,int   *iptr, int imin, int imax );
  void  checkval_F(char *varname,int nval,float *fptr,float fmin,float fmax);
  void  checkval_D(char *varname, int nval, 
		   double *dptr, double dmin, double dmax);

  int  IGNOREFILE(char *fileName);

  // ------------------------------
  // functions added 4/26/2104 

  void TABLEFILE_INIT(void);
  void tablefile_init__(void);
  void TABLEFILE_INIT_VERIFY(char *FUNNAME, char *MSG) ;

  int  TABLEFILE_OPEN(char *FILENAME, char *STRINGOPT); // return file type
  int  tablefile_open__(char *FILENAME, char *STRINGOPT);

  void TABLEFILE_CLOSE(char *FILENAME) ;
  void tablefile_close__(char *FILENAME) ;

  void STORE_TABLEFILE_COMMENT(char *comment) ;
  void store_tablefile_comment__(char *comment) ;

  int  get_TABLEFILE_TYPE(char *FILENAME) ;

  void TABLEFILE_noFile_ABORT(char *FUNNAME, char *FILENAME) ;
  void TABLEFILE_notOpen_ABORT(char *FUNNAME, char *comment) ;
  void TABLEFILE_notCompiled_ABORT(char *FILENAME, char *FORMAT, char *ENV);

  void SNTABLE_LIST(char *FILENAME);

  void SNTABLE_DUMP_VARNAMES(char *FILENAME, char *TABLENAME); 

  int  SNTABLE_DUMP_VALUES(char *FILENAME, char *TABLENAME, 
			   int NVAR, char **VARLIST, int IVAR_NPT, 
			   FILE *FP_OUTFILE,
			   char *LINEKEY_DUMP, char *SEPKEY_DUMP );


  int  SNTABLE_DUMP_OUTLIERS(char *FILENAME, char *TABLENAME, 
			     int NVAR, char **VARLIST, int IVAR_NPT,
			     float *OUTLIER_NSIGMA, FILE *FP_OUTFILE, 
			     char *LINEKEY_DUMP, char *SEPKEY_DUMP );

  void SNTABLE_SUMMARY_OUTLIERS(void);
  bool ISTABLEVAR_IFILT(char *VARNAME);

  int  SNTABLE_NEVT  (char *FILENAME, char *TABLENAME); 
  int  sntable_nevt__(char *FILENAME, char *TABLENAME);

  void SNTABLE_DEBUG_DUMP  (char *fnam, int idump);
  void sntable_debug_dump__(char *fnam, int *idump);

  // ---------------------
  // TABLE-write wrappers
  void SNTABLE_CREATE  (int  TableID, char *TableName, char *TextFormat) ;
  void sntable_create__(int *TableID, char *TableName, char *TextFormat) ;

  void SNTABLE_FILL  (int  TableID) ;
  void sntable_fill__(int *TableID ) ;


  void SNTABLE_ADDCOL  (int  TableID, char *BLOCK, void * PTRVAR, 
			char *VARLIST, int USE4TEXT );
  void sntable_addcol__(int *TableID, char *BLOCK, void* PTRVAR, 
			char *VARLIST, int *USE4TEXT );

  void parse_ADDCOL_VARLIST(char *VARLIST,
			    SNTABLE_ADDCOL_VARDEF *ADDCOL_VARDEF); 
  void parse_TABLEVAR(char *varName_with_cast, char *varName, 
		      int *ICAST, int *VECFLAG, int *ISIZE) ;

  // prep table for reading
  int  SNTABLE_READPREP(int IFILETYPE, char *TABLENAME) ;

  int SNTABLE_READPREP_VARDEF(char *varList, void *ptr, 
			      int mxlen, int optMask );
  int sntable_readprep_vardef1(char *VARNAME_withCast, void *ptr, 
			       int mxlen, int vboseflag, char *varName_noCast);
  int SNTABLE_READ_EXEC(void);

  int  IVAR_READTABLE_POINTER(char *varName) ;
  void load_READTABLE_POINTER(int IROW, int IVAR, double DVAL, char *CVAL) ;
  void load_DUMPLINE(int OPT, char *LINE, double DVAL) ;
  void load_DUMPLINE_STR(char *LINE, char *STRING) ;

  int  get_ICAST_READTBLE_POINTER(char *varName);
  int  ICAST_for_textVar(char *varName) ;

  // - - - autoStore utility functions - - - - 
  void SNTABLE_AUTOSTORE_RESET(void);
  void sntable_autostore_reset__(void);

  int SNTABLE_AUTOSTORE_INIT(char *fileName, char *tableName, 
			     char *varList, int optMask);
  int sntable_autostore_init__(char *fileName, char *tableName, 
			       char *varList,int *optMask);

  void SNTABLE_AUTOSTORE_READ(char *CCID, char *varName, int *ISTAT, 
			      double *DVAL, char *CVAL );    // output value
  void sntable_autostore_read__(char *CCID, char *varName, int *ISTAT,
				double *DVAL, char *CVAL);  // output value
  
  void fetch_autostore_ccid(int ifile, int isn, char *ccid);
  void fetch_autostore_ccid__(int *ifile, int *isn, char *ccid);

  void   SNTABLE_AUTOSTORE_malloc(int OPT, int IFILE, int IVAR);
  

  int IVAR_VARNAME_AUTOSTORE(char *varName);
  int EXIST_VARNAME_AUTOSTORE(char *varName);
  int exist_varname_autostore__(char *varName);

  int  UNIQUE_AUTOSTORE_VARNAME(int IFILE, char *VARNAME);

  // histogram wrappers
  void SNHIST_INIT(int NDIM, int ID, char *TITLE, 
		   int *NBIN, double *XMIN, double *XMAX );
  void snhist_init__(int *NDIM, int *ID, char *TITLE,
		     int *NBIN, double *XMIN, double *XMAX );

  void SNHIST_FILL(int NDIM, int ID, double *VALUE, double WGT);
  void snhist_fill__(int *NDIM, int *ID, double *VALUE, double *WGT);

  // hist-read
  void SNHIST_RDBINS(int NDIM, int ID, char *CTIT, int *NBIN,
		     double *XMIN, double *XMAX);
  void snhist_rdbins__(int *NDIM, int *ID, char *CTIT, int *NBIN,
		       double *XMIN, double *XMAX);

  void SNHIST_RDCONT(int NDIM, int ID, int NRDBIN, double *CONTENTS);
  void snhist_rdcont__(int *NDIM, int *ID, int *NRDBIN, double *CONTENTS);

  // output-directory wrappers
  
  void MAKEDIR_OUTPUT(char *CCID, int CID ) ;
  void makedir_output__(char *CCID, int *CID ) ;
  
  void CDTOPDIR_OUTPUT(void) ;
  void cdtopdir_output__ (void) ;


  // SNLCPAK functions

  void SNLCPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *VERSION_SNANA,
		    char *SURVEY_FILTERS, int SIMFLAG, 
		    int NFIT, char *TEXT_FORMAT );

  void snlcpak_init__(char *SURVEY, char *VERSION_PHOT, char *VERSION_SNANA,
		      char *SURVEY_FILTERS, int *SIMFLAG, 
		      int *NFIT, char *TEXTFMT ); 

  void SNLCPAK_NFIT_PER_SN  (int  NFIT_PER_SN);
  void snlcpak_nfit_per_sn__(int *NFIT_PER_SN);

  void SNLCPAK_SURVEY(void);  // called after output file is opened
  void snlcpak_survey__(void);

  void SNLCPAK_DISPLAYTEXT(char *CCID, char *DISPLAYTEXT) ;
  void snlcpak_displaytext__(char *CCID, char *DISPLAYTEXT) ;

  void SNLCPAK_NFIT(int NFIT) ;
  void snlcpak_nfit__(int *NFIT) ;
  
  void SNLCPAK_DATA(char *CCID, int NOBS, double *MJD, double *TOBS, 
		    double *DATA, double *DATA_ERR, int *IFILTOBS, int FLAG);
  void snlcpak_data__(char *CCID, int *NOBS, double *MJD, double *TOBS, 
		      double *DATA, double *DATA_ERR, int *IFILTOBS,int *FLAG);
  
  void SNLCPAK_FILL(char *CCID) ;
  void snlcpak_fill__(char *CCID) ;
  
  void SNLCPAK_CLEAR_SN(void); 
  void SNLCPAK_CLEAR_PLOT(void); 
  void snlcpak_clear_plot__(void);
  
  // internal wrappers that do NOT need mangled fortran fun.
  void SNLCPAK_FILL_PREP(void);
  void SNLCPAK_CHECK(char *CCID, char *comment);
  char *replace_str(char *st, const char *orig, const char *repl) ;
  
  // SPECPAK functions (Apr 2019)
  void SPECPAK_INIT(char *SURVEY, char *VERSION_PHOT, char *TEXT_FORMAT );
  void specpak_init__(char *SURVEY, char *VERSION_PHOT, char *TEXTFMT ); 
  void SPECPAK_CLEAR_PLOT(void); 
  void specpak_clear_plot__(void);

  void SPECPAK_DATA(char *CCID, int IDSPEC, double MJD, double Tobs, 
		    double Texpose,int NLAMBIN,double *LAMMIN, double *LAMMAX,
		    double *FLAM,double *FLAMERR);
  void specpak_data__(char *CCID,int *IDSPEC, double *MJD, double *Tobs,
		      double *Texpose, 
		      int *NLAMBIN, double *LAMMIN, double *LAMMAX,
		      double *FLAM,double *FLAMERR);

  void SPECPAK_FILL(char *CCID) ;
  void specpak_fill__(char *CCID) ;

  // --- ISFILE_xxx functions
  int ISFILE_HBOOK(char *fileName);
  int ISFILE_ROOT(char *fileName);
  int ISFILE_TEXT(char *fileName);
  int ISFILE_MARZ(char *fileName);
				 
#ifdef __cplusplus
}          
#endif


// END
