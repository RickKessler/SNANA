/*******************************************************
     Created Jan 2005 by R.Kessler

     Defines:
      - Generic tools to read input from file
      - error handling routines
      - physical constants.


  Aug 22 2013:  define MODEL_S11DM15

  Aug 28 2013:  define struct SNTOOLS_FITRES to hold FITRES-variables;
                avoids accidental conflicts with other programs.

  Mar 20 2017: remove special chars from  FILTERSTRING_DEFAULT, except
               for '&' to add default SYN_SPECTROGRAPH filter if
               user does not define them.

  Apr 6 2017: define EXEC_CIDMASK()

  Jul 5 2017: add read_SURVEY

  Nov 2017: add XXX_PATH_SNDATA_SIM tools to keep track of
            user-defined output directories for simulations.

  Dec 10 2017: declare SNANA_VERSION here instad of snana.car
               so that it is more accessible. Also define function
               get_SNANA_VERSION.

  May 30 2020:
    MXWORDFILE_PARSE_WORDS -> 1 million (was 500k) to handle
    data files with lots of spectra.

  Jul 30 2020:
    + optional pre-proc flag ONE_RANDOM_STREAM to disable 2nd
      stream for Mac compilation. Default is still 2 streams.

  Jul 31 2020:
    + MXCHARWORD_PARSE_WORDS -> MXPATHLEN=300  (was 60) to handle file names

  Nov 12 2020:
    + Adding glob utility

  Feb 4 2021: add python-like dictionary utility (e.g., for FIELD-dependence)
  Jun 2 2021: MXWORDFILE_PARSE_WORDS -> 2M (was 1 million)
  
********************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <glob.h>

#include "sndata.h"

#include "sntools_gridmap.h"
#include "sntools_genGauss_asym.h"
#include "sntools_genExpHalfGauss.h"

// --------------------------------------------------
#define  SNANA_VERSION_CURRENT  "v11_04l"        
//#define  ONE_RANDOM_STREAM  // enable this for Mac (D.Jones, July 2020)
//#define  MACOS              // another MAC OS option, D.Jones, Sep 2020

#define KEYNAME_DOCANA_REQUIRED   "DOCUMENTATION:"
#define KEYNAME2_DOCANA_REQUIRED  "DOCUMENTATION_END:"
#define OPENMASK_VERBOSE        1  // see snana_openTextFile
#define OPENMASK_REQUIRE_DOCANA 2  // see snana_openTextFile
#define OPENMASK_IGNORE_DOCANA  4  // see snana_openTextFile
#define MXLINE_DOCANA         500  // max number of DOCANA lines

// default cosmo params from Planck 2018 (https://arxiv.org/abs/1807.06209)
#define OMEGA_MATTER_DEFAULT   0.315
#define OMEGA_LAMBDA_DEFAULT   0.685
#define w0_DEFAULT            -1.0
#define wa_DEFAULT             0.0
#define H0_SALT2            70.0    // km/s/Mpc : tied to SALT2 training
#define H0_MLCS             65.0    // km/s/Mpc : tied to MLCS training
#define H0_Planck          67.4   // 1807.06209 (Planck 2018)
#define H0_SH0ES           74.03  // 1903.07603 (Riess 2019)

#define LIGHT_km  2.99792458e5      // speed of light (km/s)
#define LIGHT_A   2.99792458e18     // speed of light (A/s)
#define PLANCK    6.6260755e-27     // Planck constant (erg s)
#define hc        LIGHT_A * PLANCK
#define PC_km     3.085678e13       // parsec (km)
#define TWOPI     6.28318530718
#define RADIAN    TWOPI / 360.0     // added Oct 2010
#define ZAT10PC    2.335e-9         // redshift at 10pc (H0=70)
#define ZMAX_SNANA   4.0        // max snana redshift, Dec 26 2016
#define PSFMAX_SNANA 5.0        // max allowed PSF, FWHM, arcsec (Mar 2021)
#define COMMA      ","              // to split comma-sep strings
#define COLON      ":"              // to split colon-sep strings
#define PERCENT    "%"              // idem for %-sep strings
#define PLUS       "+"
#define STAR       "*"
#define DOT        "."

// from Planck 2018 (installed June 8 2020)
#define  CMBapex_l  (double)264.031    // deg (galactic coords !!!)
#define  CMBapex_b  (double)48.253     // deg
#define  CMBapex_v  (double)369.82    // km/sec

#define FWHM_SIGMA_RATIO  2.3548    // FWHM/sigma for Gaussian
#define TEN        (double)10.0
#define LNTEN      (double)2.30259  // ln(10)

#define REJECT_FLAG    -1
#define ACCEPT_FLAG    +1
#define ERROR          -1  // use these as return args for int functions
#define SUCCESS        +1
#define NULLINT        -9
#define NULLFLOAT      -9.0
#define NULLDOUBLE     (double)-9.0
#define NULLSTRING     "NULL"
#define NODOUBLE       (double)1.7777E14
#define NOFLOAT        (float)1.7777E14
#define NOINT          555444333
#define BLANK_STRING   ""
#define NOTSET_STRING  "NOTSET" // init str pointers to avoid compile error
#define NULLTYPE    0   // SN type with no TYPE in data base.
#define FIELD_NONAME  "VOID"  // avoid using 'NULL' to avoid python problem.

#define INDEX_NOTSATURATE 0
#define INDEX_SATURATE    1
#define FLUXCALERR_SATURATE 1.0E8
#define MAG_SATURATE   -7.0   // saturated mag
#define MAG_ZEROFLUX   99.0   // expect flux=0 (e.g., pre-explosion)
#define MAG_NEGFLUX   128.0   // fluctuation results in flux < 0
#define MAG_UNDEFINED 666.0   // model mag wave-range doesn't cover filter
#define MAGERR_UNDEFINED 9.0
#define FLUX_UNDEFINED  -9.0  // undefined model flux
#define HOSTLIB_SNPAR_UNDEFINED    -9999.0 
#define HOSTLIB_IGAL_UNDEFINED     -9999
#define HOSTLIB_PROPERTY_UNDEFINED -9999.0    // Feb 10 2022

#define MODEL_STRETCH  1 // single-stretch (use rise/fall fudges for 2-stretch)
#define MODEL_MLCS2k2  3  // MLCS2k2 (Riess, Suarab, Hubert)
#define MODEL_SNOOPY   4  // dm15 model from C.Burns et al, 2011
#define MODEL_S11DM15  5  // hybrid-DM15 model for psnid, Sako et al, 2011
#define MODEL_SALT2    6  // SALT-II model
#define MODEL_SALT3    9  // reserve index for SALT3 if needed
#define MODEL_SIMSED   7  // set of SED sequences (e.g., sim-explosion models)
#define MODEL_BYOSED   8  // build-your-own SED model
#define MODEL_SNEMO    9  // SNEMO from SNFactory (Sep 2020)
#define MODEL_BAYESN    13  // BayeSN (Nov 2021)
#define MODEL_NON1ASED   10  // obs-frame NONIA from SED
#define MODEL_NON1AGRID  11  // obs-frame NONIA from mag-grid (Mar 2016)
#define MODEL_LCLIB      12  // light curve library (July 2017)
#define MODEL_FIXMAG     20  // constant/fixed mag (default is mag=20)
#define MODEL_RANMAG     21  // constant/fixed mag per LC, random each event
#define MODEL_SIMLIB     22  // use MAGOBS column of SIMLIB
#define MXMODEL_INDEX    40
#define PATH_SNDATA_SIM_LIST  "PATH_SNDATA_SIM.LIST"


#define KEYSOURCE_FILE 1 // moved from snlc_sim.h on June 24th 2021
#define KEYSOURCE_ARG  2

// keep this in sync with the fortran FILTDEF_STRING
// Oct 22 2015: add 27 special chars to hack an IFU
#define FILTERSTRING_DEFAULT  " ugrizYJHK UBVRIXy0123456789 abcdef ACDEFGLMNOPQSTWZ hjklmnopqstvwx &"

#define PRIVATE_MODELPATH_NAME "SNANA_MODELPATH"  // name of optional env

#define PARNAME_EBV "EBV"
#define PARNAME_AV  "AV"
#define PARNAME_RV  "RV"

char FILTERSTRING[MXFILTINDX] ;

// define variables for random number list
#define MXLIST_RAN      4  // max number of lists for stream0
#define MXSTORE_RAN  1000  // size of each RANLIST (for each event)
#define MXSTREAM_RAN    2  // max number of independent streams
#define BUFSIZE_RAN   256

struct {
  int     NSTREAM ; // number of srandom streams (legacy is 1)
  double  RANSTORE[MXLIST_RAN+1][MXSTORE_RAN] ;
  int     NLIST_RAN;         // Number of lists
  int     NSTORE_RAN[MXLIST_RAN+1] ;
  double  RANFIRST[MXLIST_RAN+1], RANLAST[MXLIST_RAN+1]; // for syncing.

  // for multi-stream randoms
  // xxx random_data_def ranStream[MXSTREAM_RAN];

#ifdef MACOS
  struct random_data { int someInteger; };  // from D.Jones, Sep 2020
#endif

  struct random_data  ranStream[MXSTREAM_RAN];
  char stateBuf[MXSTREAM_RAN][BUFSIZE_RAN];

  // wrap-around stats for how often each random is re-used.
  int    NCALL_fill_RANSTATs;
  double NWRAP_MIN[MXLIST_RAN+1] ;
  double NWRAP_MAX[MXLIST_RAN+1] ;
  double NWRAP_AVG[MXLIST_RAN+1] ;
  double NWRAP_RMS[MXLIST_RAN+1] ;
  double NWRAP[MXLIST_RAN+1] ;
  double NWRAP_SUM[MXLIST_RAN+1] ;
  double NWRAP_SUMSQ[MXLIST_RAN+1] ;

} GENRAN_INFO ;


// errmsg parameters
char c1err[200];   // for kcorerr utility
char c2err[200];   // for kcorerr utility
char BANNER[200];
int  EXIT_ERRCODE;  // program error code set by program (Jan 2019)

#define MXARGV 200  // 100->200 Apr 18 2022
int  NARGV_LIST;
char *ARGV_LIST[MXARGV];
int  USE_ARGV_LIST[MXARGV];  // 1 => line arg was used, 0=> not used.

// ---------
// generic DUMP_STRING for any function
#define MXCID_DUMP 10
#define DUMPFLAG_noABORT    1
#define DUMPFLAG_withABORT  2
struct {
  char  FUNNAME[80] ; // name of funtion to dump
  int   NCID ;        // number of CIDs to dump
  char  CCIDLIST[MXCID_DUMP][40];
  int   NFILT;
  char  FILTLIST[MXFILTINDX];
  int   ABORTFLAG ;  // 1 --> abort after dumping all CIDLIST

  double MJDRANGE[2] ;

  // number of DONE filters per CID
  int NFILT_DONE[MXCID_DUMP] ;

} DUMP_STRING_INFO ;

// Store GaussIntegrals from 0 to x
struct {
  int    INIT_FLAG;
  int    NBIN_XMAX ;
  double XMAX;
  double XMAX_LIST[MXSTORE_RAN];
  double GINT_LIST[MXSTORE_RAN];
} GAUSS_INTEGRAL_STORAGE ;

// variables for Landolt transformations (UBVRI,BX)
#define NFILT_LANDOLT 6
double LANDOLT_MAGPRIMARY[6];     // UBVRI,X synthetic mag for primary ref
double LANDOLT_COLOR_VALUE[6];   // Marriner's k0 - k4
double LANDOLT_COLOR_ERROR[6];   // errors on above


// offsets for quickly reducing n-dim indices to 1d-index
#define MXMAP_1DINDEX 200 // max number of 1d index maps to store
#define MXDIM_1DINDEX 10  // max number of dimenstions


int OFFSET_1DINDEX[MXMAP_1DINDEX][MXDIM_1DINDEX] ;
int NPT_PERDIM_1DINDEX[MXMAP_1DINDEX][MXDIM_1DINDEX] ;


// Feb 2021; define structure to enable python-like dictionary
typedef struct STRING_DICT_DEF {
  char STRING[100];     // original string to parse
  char NAME[60] ;       // name of dictionary
  int  N_ITEM;          // number of items to store
  int  MAX_ITEM;        // max number of items
  char **STRING_LIST ;  // list of string-item names
  double *VALUE_LIST ;  // double value for each string

  // if string is same as LAST_STRING, return LAST_VALUE;
  // this avoids looping thru lists of string-matching.
  char   LAST_STRING[60] ;
  double LAST_VALUE ;

} STRING_DICT_DEF ;


// define table info corresponding to hash storage
struct {
  int  NVAR, *IVAR_TABLE;
  char   **VARNAME_LIST; // table var name vs. ivar
  double **VAL_LIST;     // table value vs. ivar and isn
} HASH_STORAGE;

// Mar 2019: define user-input polynomial typedef with arbitrary order.
#define MXORDER_GENPOLY 20
typedef struct  {
  int ORDER;       // 2 -> 2nd order -> a + b*x + c*x^2
  double COEFF_RANGE[MXORDER_GENPOLY][2]; // range for each coeff.
  char   STRING[200]; // string that was parsed to get COEFF_RANGEs
  char   VARNAME[40]; // optional variable name
} GENPOLY_DEF ;



// Feb 2020: structure to handle correlated randoms using Cholesky decomp
typedef struct {
  int    MATSIZE;      // matrix size along each dimension
  double *COVMAT1D ;   // user-input COV matrix
  double **CHOLESKY2D; // Cholesky decomp matrix used to get correlated ran
  //  gsl_matrix_view chk; // internal matrix
} CHOLESKY_DECOMP_DEF ;


#define MXFILT_REMAP 20
struct  {
  int  NMAP ;
  char MAPSTRING[MXFILT_REMAP][20] ;
  int  IFILTOBS_MAP[MXFILTINDX] ;  // ifiltobs_orig -> ifiltobs_remap
  int  IFILT_MAP[MXFILTINDX] ;     // ifiltobs_orig -> ifilt_remap
} FILTER_REMAP  ;


// structure to hold warnings for old/obsolete inputs that should
// be replaced.
struct {
  int  NWARN ;
  char VARNAME_OLD[20][80];  // [iwarn][lenvarname]
  char VARNAME_NEW[20][80];  // [iwarn][lenvarname]
} OLD_INPUTS;


#define ADDBUF_PARSE_WORDS 10000
#define MXCHARWORD_PARSE_WORDS MXPATHLEN // MXCHAR per word
#define MXCHARLINE_PARSE_WORDS 2000      // max chars per line
#define MXWORDLINE_PARSE_WORDS  700      // max words per line
#define MXWORDFILE_PARSE_WORDS 2000000   // max words to parse in a file

#define MXWORDLINE_FLUX       10  // max words per line in SED file
#define MXCHARLINE_FLUX      120  // max char per line to read from SED

#define MSKOPT_PARSE_WORDS_FILE    1   // parse words in a file
#define MSKOPT_PARSE_WORDS_STRING  2   // parse string
#define MSKOPT_PARSE_WORDS_IGNORECOMMA    4  // parse blank space; ignore comma
#define MSKOPT_PARSE_WORDS_IGNORECOMMENT  8  // ignore after comment char
#define MSKOPT_PARSE_WORDS_FIRSTLINE     16  // read only 1st line only
#define LANGFLAG_PARSE_WORDS_C  0  // PARSE_WORDS language flag for C

struct {
  char  FILENAME[MXPATHLEN];
  int   BUFSIZE ;
  int   NWD;
  char **WDLIST;
} PARSE_WORDS ;


struct {
  int NFILE;
  char **FILE_NAMES;
} GLOB_LIST;

#define MXLIST_STRING_UNIQUE  200
#define MXLIST_KEY_UNIQUE    1000  // for key-dump only
#define STRINGMATCH_INIT      "INIT"
#define STRINGMATCH_KEY_DUMP  "KEY_DUMP"
struct {
  int  NLIST;
  char SOURCE_of_STRING[200]; // for messages
  char STRING[MXLIST_STRING_UNIQUE][60];  // list of uniqe strings
  int  NFOUND_STRING[MXLIST_STRING_UNIQUE];  // number of times each str found


  bool DUMPKEY_FLAG ;
  int  NKEY; // number of unique keys stored for dump
  char *KEY[MXLIST_KEY_UNIQUE] ;

} STRING_UNIQUE ;

struct {
  int     LAST_NOBS ;  // to know when a new malloc is needed.
  double *TLIST_SORTED, *MAGLIST_SORTED, *FLUXLIST_SORTED ;
  int    *INDEX_SORT;
} LCWIDTH ;


#define OPTMASK_SETPKMJD_FLUXMAX   8 // naive max flux
#define OPTMASK_SETPKMJD_FLUXMAX2 16 // Fmax-clump method, wgt=1 per obs
#define OPTMASK_SETPKMJD_FLUXMAX3 32 // idem, with wgt=log(SNR) per obs
#define OPTMASK_SETPKMJD_TRIGGER  64 // return MJD_TRIGGER

struct {
  int    OPTMASK;
  double SNRCUT, SNRCUT_BACKUP;
  double MJDWIN ;
} INPUTS_OBS_atFLUXMAX ;


#define MXFILE_ENVreplace 100
struct {
  int   NFILE;
  char  FILENAME_ORIG[MXFILE_ENVreplace][MXPATHLEN];
  char  FILENAME_ENVreplace[MXFILE_ENVreplace][MXPATHLEN];
}  ENVreplace_store;


#define MXVAR_WRITE_EPOCH_LIST 8
#define CUTMODE_REJECT -1
#define CUTMODE_ACCEPT +1
#define CUTTYPE_BITMASK 1
#define CUTTYPE_WINDOW  2
struct {
  FILE   *FP_OUT;
  int    NVAR;  // number of variables to apply cut
  char   VARNAME[MXVAR_WRITE_EPOCH_LIST][40] ;
  char   VARNAMES_LIST[200]; //
  double CUTWIN[MXVAR_WRITE_EPOCH_LIST][2] ;
  int    CUTMODE[MXVAR_WRITE_EPOCH_LIST] ;   // accept or reject
  int    CUTTYPE[MXVAR_WRITE_EPOCH_LIST] ;   // window or bit-mask

  char ROWKEY[20]; // REJECT: or ACCEPT:
  int  CUTMASK_ALL ; // cutmask if all cuts are satisfied

  int  NEPOCH_ALL ;    // number of epochs tested
  int  NEPOCH_WRITE ;  // number of epochs written
  int  NEPOCH_CUTFAIL[MXVAR_WRITE_EPOCH_LIST]; // NFAIL per cut
  int  NEPOCH_CUTFAIL_ONLY[MXVAR_WRITE_EPOCH_LIST]; // idem, only cut

} WRITE_EPOCH_LIST ;



// ##############################################################
//
//     functions
//
// ##############################################################

void write_epoch_list_init(char *outFile);
void write_epoch_list_addvar(char *varName, double *CUTWIN, char *CUTMODE);
void write_epoch_list_exec(char *CID,double MJD,char *BAND, double *VALUES);
void write_epoch_list_summary(void);

void write_epoch_list_init__(char *outFile);
void write_epoch_list_addvar__(char *varName, double *CUTWIN, char *CUTMODE);
void write_epoch_list_exec__(char *CID,double *MJD,char *BAND,double *VARLIST);
void write_epoch_list_summary__(void);

void catVarList_with_comma(char *varList, char *addVarName);

void init_Cholesky(int OPT, CHOLESKY_DECOMP_DEF *DECOMP ) ;
void getRan_GaussCorr(CHOLESKY_DECOMP_DEF *DECOMP,
		      double *RanList_noCorr, double *RanList_corr);

// xxx void GaussRanCorr(CHOLESKY_DECOMP_DEF *DECOMP,
// xxx		  double *RanList_noCorr, double *RanList_corr);

void INIT_SNANA_DUMP(char *STRING);
int  CHECK_SNANA_DUMP(char *FUNNAME, char *CCID, char *BAND, double MJD );
void ABORT_SNANA_DUMP(void) ;

void init_snana_dump__(char *STRING);
int  check_snana_dump__(char *FUNNAME, char *CCID, char *BAND, double *MJD);
void abort_snana_dump__(void);

// ---- utility to REMAP filters for tables -----
void FILTER_REMAP_INIT(char *remapString, char *VALID_FILTERLIST,
		       int *NFILT_REMAP,
		       int *IFILTLIST_REMAP, char *FILTLIST_REMAP);
void FILTER_REMAP_FETCH(int IFILTOBS_ORIG,
			int *IFILTOBS_REMAP, int *IFILT_REMAP);

void filter_remap_init__(char *remapString, char *VALID_FILTERLIST,
			 int *NFILT_REMAP,
			 int *IFILTLIST_REMAP, char *FILTLIST_REMAP);

void filter_remap_fetch(int *IFILTOBS_ORIG,
			int *IFILTOBS_REMAP, int *IFILT_REMAP);

int intrac_() ;  // needed by fortran minuit

void get_SNANA_VERSION(char *snana_version); // Dec 2017
void get_snana_version__(char *snana_version);

float get_SNANA_VERSION_FLOAT(char *snana_version); // Oct 2020
float get_snana_version_float__(char *snana_version);

bool correct_sign_vpec_data(char *snana_version_data);
bool correct_sign_vpec_data__(char *snana_version_data);

void print_full_command(FILE *fp, int argc, char** argv);
void print_KEYwarning(int ISEV, char *key_old, char *key_new);

void set_FILTERSTRING(char *FILTERSTRING) ;
void set_EXIT_ERRCODE(int ERRCODE);
void set_exit_errcode__(int *ERRCODE);

/* xxxxxxxx mark delete Dec 10 2021 xxxxxxxxx
void copy_SNDATA_GLOBAL(int copyFlag, char *key,
			int NVAL, char *stringVal, double *parVal);
void copy_SNDATA_HEAD(int copyFlag, char *key,
		      int NVAL, char *stringVal, double *parVal);
void copy_SNDATA_OBS(int copyFlag, char *key,
		     int NVAL,char *stringVal, double *parVal);
int select_MJD_SNDATA(double *CUTWIN_MJD);

void copy_GENSPEC(int copyFlag, char *key, int ispec, double *parVal);

void copy_int(int copyFlag, double *DVAL0, int    *IVAL1) ;
void copy_lli(int copyFlag, double *DVAL0, long long  *IVAL1) ;
void copy_flt(int copyFlag, double *DVAL0, float  *FVAL1) ;
void copy_dbl(int copyFlag, double *DVAL0, double *DVAL1) ;
void copy_str(int copyFlag, char   *STR0,  char   *STR1 );
//void copy_void(int copyFlag, double *DVAL0, void   *VAL1) ;
xxxxxxxxxxxxxx end mark xxxxxxxxxxx */

int  IGNOREFILE(char *fileName);
int  ignorefile_(char *fileName);

int strcmp_ignoreCase(char *str1, char *str2) ;


// data functions
void clr_VERSION ( char *version, int prompt );  // remove old *version files
int init_VERSION ( char *version);  // init VERSION_INFO struct
int init_SNPATH(void);
int init_SNDATA_EVENT(void) ;  // init SNDATA struct (for each event)
int init_SNDATA_GLOBAL(void); // init SNDATA globals (one-time init)
void set_SNDATA_FILTER(char *filter_list);

void init_GENSPEC_GLOBAL(void) ;
void init_GENSPEC_EVENT(int ISPEC, int NBLAM);

void ld_null(float *ptr, float value);

int rd_SNDATA(void);
int rd_filtband_int   ( FILE *fp, int isn, int i_epoch, int   *iptr ) ;
int rd_filtband_float ( FILE *fp, int isn, int i_epoch, float *fptr ) ;

int rd_sedFlux( char *sedFile, char *sedcomment,
		double DAYrange[2], double LAMrange[2],
		int MXDAY, int MXLAM, int OPTMASK,
		int *NDAY, double *DAY_LIST, double *DAY_STEP,
		int *NLAM, double *LAM_LIST, double *LAM_STEP,
		double *FLUX_LIST, double *FLUXERR_LIST );

int  PARSE_FILTLIST(char *filtlist_string, int *filtlist_array );

void update_covMatrix(char *name, int OPTMASK, int MATSIZE,
		      double (*covMat)[MATSIZE], double EIGMIN, int *istat_cov);
void update_covmatrix__(char *name, int *OPTMASK, int *MATSIZE,
			double (*covMat)[*MATSIZE], double *EIGMIN,
			int *istat_cov ) ;

int  store_PARSE_WORDS(int OPT, char *FILENAME);
void malloc_PARSE_WORDS(void);
void get_PARSE_WORD(int langFlag, int iwd, char *word);
void get_PARSE_WORD_INT(int langFlag, int iwd, int   *i_val);
void get_PARSE_WORD_FLT(int langFlag, int iwd, float *f_val);
void get_PARSE_WORD_NFLT(int langFlag, int NFLT, int iwd, float *f_val);
void get_PARSE_WORD_NFILTDEF(int langFlag, int iwd, float *f_val);

void get_PARSE_WORD_DBL(int langFlag, int iwd, double *d_val);
void get_parse_word__(int *langFlag, int *iwd, char  *word );
void get_parse_word_int__(int *langFlag, int *iwd, int   *i_val);
void get_parse_word_flt__(int *langFlag, int *iwd, float *f_val);
void get_parse_word_dbl__(int *langFlag, int *iwd, double *d_val);

int  match_cidlist_init(char *fileName, int *OPTMASK, char *varList_store);
int  match_cidlist_init__(char *fileName, int *OPTMASK, char *varList_store);

int  match_cidlist_exec(char *cid);
int  match_cidlist_exec__(char *cid);

double  match_cidlist_parval(int isn_match, char *varName, int abort_flag);
double  match_cidlist_parval__(int *isn_match, char *varName, int *abort_flag);


void   init_GENPOLY(GENPOLY_DEF *GENPOLY);
void   parse_GENPOLY(char *stringPoly, char *varName,
		     GENPOLY_DEF *GENPOLY, char *callFun );
double eval_GENPOLY(double VAL, GENPOLY_DEF *GENPOLY, char *callFun);
void   copy_GENPOLY(GENPOLY_DEF *GENPOLY_IN, GENPOLY_DEF *GENPOLY_OUT);

void parse_multiplier(char *inString, char *key, double *multiplier);
void check_uniform_bins(int NBIN, double *VAL, char *comment_forAbort);
void check_argv(void);
void check_arg_len(char *keyName, char *arg, int MXLEN) ;
void check_EOF(FILE *fp, char *file_name, char *fun_name, int nline_read) ;

double NoiseEquivAperture(double PSFSIG1, double PSFSIG2, double PSFratio);
double noiseequivaperture_(double *PSFSIG1, double *PSFSIG2, double *PSFratio);

double modelmag_extrap(double T, double Tref,
		       double magref, double magslope, int LDMP);
double modelflux_extrap(double T, double Tref,
			double fluxref, double fluxslope, int LDMP);

int fluxcal_SNDATA ( int iepoch, char *magfun, int opt ) ;
double asinhinv(double mag, int ifilt);


float effective_aperture ( float PSF_sigma, int VBOSE ) ;

/* xxxxxxxx legacy write function to hopefull delete soon (Feb 2021)
int   wr_SNDATA(int IFLAG_WR, int IFLAG_DBUG);
void  wr_HOSTGAL(FILE *fp);
void  wr_SIMKCOR(FILE *fp, int EPMIN, int EPMAX ) ;
int  wr_filtband_int   ( FILE *fp, char *keyword,
			 int NINT, int *iptr, char *comment, int opt );
int  wr_filtband_float ( FILE *fp, char *keyword,
			 int NFLOAT, float *fptr, char *comment, int idec );
int  header_merge(FILE *fp, char *auxheader_file);
int  sort_epochs_bymjd(void);
int   WRSTAT ( int wrflag, float value ) ;
xxxxxxxxxxxxxxxxxxxxx */

int   Landolt_ini(int opt, float *mag, float *kshift);
int   landolt_ini__(int *opt, float *mag, float *kshift);
int   Landolt_convert(int opt, double *mag_in, double *mag_out);
int   landolt_convert__(int *opt, double *mag_in, double *mag_out);

// errmsg" utilities
void  tabs_ABORT(int NTAB, char *fileName, char *callFun);
void  missingKey_ABORT(char *key, char *file, char *callFun) ;
void  legacyKey_abort(char *callFun,  char *legacyKey, char *newKey) ;

void  errlog ( FILE *fp, int  isev, char *fnam, char *msg1, char *msg2 );
void  errmsg ( int  isev, int iprompt, char *fnam, char *msg1, char *msg2 );
void  errmsg_( int *isev,int *iprompt, char *fnam, char *msg1, char *msg2 );

void  prompt(char *msg) ;
void  madend(FILE *fp, int flag);   // indicates bad end of program
void  happyend(void) ;    // happy end of program
void  parse_err ( char *infile, int NEWMJD, char *keyword );
void  print_preAbort_banner(char *fnam);
void  print_preabort_banner__(char *fnam);

void  warn_oldInputs(char *varName_old, char *varName_new);
void  warn_oldinputs__(char *varName_old, char* varName_new) ;

void  checkval_I(char *varname, int nval,
		 int *iptr, int imin, int imax );
void  checkval_F(char *varname, int nval,
		 float *fptr, float fmin, float fmax );
void  checkval_D(char *varname, int nval,
		 double *dptr, double dmin, double dmax );

void  checkval_i__(char *varname, int *nval,
		   int *iptr, int *imin, int *imax );


void  checkArrayBound(int i, int MIN, int MAX,
		      char *varName, char *comment, char *funName);
void  checkArrayBound_(int *i, int *MIN, int *MAX,
		       char *varName, char *comment, char *funName);

void  check_magUndefined(double mag, char *varName, char *callFun );
void  checkStringUnique(int MAX, char *string, char *msgSource, char *callFun);
bool  NstringMatch(int MAX, char *string, char *key);
bool  uniqueMatch(char *string, char *key);
int   uniqueOverlap(char *string, char *key);
int   keyMatch(char *string, char *key, char *keySuffix_optional);
bool  keyMatchSim(int MXKEY, char *KEY, char *WORD, int keySource);
void  dumpUniqueKey(char *key) ;

int   ivar_matchList(char *varName, int NVAR, char **varList);
int   match_cid_hash(char *cid, int ilist, int isn);
int   match_cid_hash__(char *cid, int *ilist, int *isn);

void read_VARNAMES_KEYS(FILE *fp, int MXVAR, int NVAR_SKIP, char *callFun,
			int *NVAR, int *NKEY, int *UNIQUE, char **VARNAMES );

unsigned int *CIDMASK_LIST;  int  MXCIDMASK, NCIDMASK_LIST ;
int  exec_cidmask(int mode, int CID);
int  exec_cidmask__(int *mode, int *CID);
void test_cidmask(void) ;

// parameters for errmsg utility
#define SEV_INFO   1  // severity flag => give info
#define SEV_WARN   2  // severity flag => give warning
#define SEV_ERROR  3  // severity flag => error
#define SEV_FATAL  4  // severity flag => abort program
#define EXIT_ERRCODE_KCOR           10
#define EXIT_ERRCODE_SIM            11
#define EXIT_ERRCODE_SALT2mu        12
#define EXIT_ERRCODE_combine_fitres 13
#define EXIT_ERRCODE_sntable_dump   14
#define EXIT_ERRCODE_wfit           15
#define EXIT_ERRCODE_merge_root     16
#define EXIT_ERRCODE_merge_hbook    17
#define EXIT_ERRCODE_UNKNOWN        99

// define old useful functions for reading/parsing input file
void  readint(FILE *fp, int nint, int *list) ;
void  readlong(FILE *fp, int nint, long long *list) ;
void  readfloat(FILE *fp, int nint, float *list);
void  readdouble(FILE *fp, int nint, double *list);
void  readchar(FILE *fp, char *clist) ;

int read_genpoly(char *KEYNAME, char **WORDS, int order_legacy,
                 GENPOLY_DEF *POLY) ;

void  read_SURVEYDEF(void);
int   get_IDSURVEY(char *SURVEY);

void  read_redshift(FILE *fp, float *redshift, float *redshift_err );
int   gtchars(char *string, char **argv);
void  debugexit(char *string);

int    INTFILTER ( char *cfilt );  // convert filter name into index


// miscellaneous

double MAGLIMIT_calculator(double ZPT, double PSF, double SKYMAG, double SNR);
double SNR_calculator(double ZPT, double PSF, double SKYMAG, double MAG,
		      double *FLUX_and_ERR ) ;

int getInfo_PHOTOMETRY_VERSION(char *VERSION,  char *DATADIR,
			       char *LISTFILE, char *READMEFILE);
int getinfo_photometry_version__(char *VERSION,  char *DATADIR,
				 char *LISTFILE, char *READMEFILE);

int file_timeDif(char *file1, char *file2);

void warn_NVAR_KEY(char *fileName);

void fillbins(int OPT, char *name, int NBIN, float *RANGE,
	      float *BINSIZE, float *GRIDVAL);

int commentchar(char *str);
int rd2columnFile(char *file, int MXROW, int *Nrow,
		   double *column1, double *column2 );
int rd2columnfile_(char *file, int *MXROW, int *Nrow,
		   double *column1, double *column2 );

int nrow_read(char *file, char *callFun) ;

// -----
double interp_SINFUN(double VAL, double *VALREF, double *FUNREF,
		     char *ABORT_COMMENT ) ;

#define OPT_INTERP_LINEAR    1
#define OPT_INTERP_QUADRATIC 2
double interp_1DFUN(int opt, double val, int NBIN,
		    double *VAL_LIST, double *FUN_LIST,
		    char *abort_comment );

double interp_1dfun__(int *opt, double *val, int *NBIN,
		      double *VAL_LIST, double *FUN_LIST,
		      char *abort_comment );

int quickBinSearch(double VAL, int NBIN, double *VAL_LIST,
		   char *abort_comment);
double quadInterp ( double VAL, double VAL_LIST[3], double FUN_LIST[3],
		    char *abort_comment );

double polyEval(int N, double *coef, double x);

void arrayStat(int N, double *array, double *AVG, double *RMS, double *MEDIAN);
void arraystat_(int *N, double *array, double *AVG, double *RMS, double *MEDIAN);
void test_arrayStat(void);
double STD_from_SUMS(int N, double SUM, double SQSUM);
double sigint_muresid_list(int N, double *MURES_LIST, double *MUCOV_LIST,
			   int OPTMASK, char *callFun );

void trim_blank_spaces(char *string) ;
void remove_string_termination(char *STRING, int LEN) ;

void splitString(char *string, char *sep, int MXsplit,
		 int *Nsplit, char **ptrSplit);
void splitString2(char *string, char *sep, int MXsplit,
		  int *Nsplit, char **ptrSplit) ;
void split2floats(char *string, char *sep, float *fval) ;

void remove_quote(char *string);
void extractStringOpt ( char *string, char *stringOpt) ;
void extractstringopt_( char *string, char *stringOpt) ;
void extract_MODELNAME(char *STRING, char *MODELPATH, char *MODELNAME);
void extract_modelname__(char *STRING, char *MODELPATH, char *MODELNAME);

void parse_commaSepList(char *item_name, char *item, int MAX_ITEM, int MXCHAR,
			int *n_item, char ***arrayList );
int  index_charString(char *c, char *s);

double PROB_Chi2Ndof(double chi2, int Ndof); // replace CERNLIB's PROB function
double prob_chi2ndof__(double *chi2, int *Ndof);

// pixel-distance to nearest edge
float edgedist(float XPIX, float YPIX,  int NXPIX, int NYPIX);

void print_banner ( const char *banner ) ;
void fprint_banner (FILE *FP, const char *banner ) ;

// shells to open text file
void find_pathfile(char *fileName, char *PATH_LIST, char *FILENAME, char *callFun);

FILE *open_TEXTgz(char *FILENAME, const char *mode,int *GZIPFLAG) ;
FILE *snana_openTextFile (int OPTMASK, char *PATH_LIST, char *fileName,
			  char *fullName, int *gzipFlag );
void snana_rewind(FILE *fp, char *FILENAME, int GZIPFLAG);
void abort_openTextFile(char *keyName, char *PATH_LIST,
			char *fileName, char *funCall);
bool check_openFile_docana(bool REQUIRE_DOCANA, FILE *fp, char *fileName); // check file is open
void check_file_docana(int optmask, char *fileName);   // open file and check
void check_file_docana__(int *optmask, char *fileName);

bool key_in_file(char *file, char *key, int nwd_check);
int  colnum_in_table(char *fileName, char *varName);

void react_missing_docana(bool FOUND_DOCANA, char *fileName);
void react_missing_docana__(bool *FOUND_DOCANA, char *fileName);
void abort_docana_tooLong(char *file, char *callFun);

void abort_bad_input(char *key,  char *word, int iArg, char *callFun);

int  ENVreplace(char *fileName, char *callFun, int ABORTFLAG);
void ENVrestore(char *fileName_noENV, char *fileName_orig);

void  init_string_dict(STRING_DICT_DEF *DICT, char *NAME, int MAXITEM);
void  load_string_dict(STRING_DICT_DEF *DICT, char *string, double val);
double get_string_dict(int OPT, char *string, STRING_DICT_DEF *DICT);
void  dump_string_dict(STRING_DICT_DEF *DICT);
void  parse_string_prescales(char *STRING, STRING_DICT_DEF *DICT);

// SLALIB functions translated by D. Cinabro
void slaEqgal ( double dr, double dd, double *dl, double *db );
void slaDcs2c ( double a, double b, double v[3] );
void slaDmxv ( double dm[3][3], double va[3], double vb[3] );
void slaDcc2s ( double v[3], double *a, double *b );
double slaDrange ( double angle );
double slaDranrm ( double angle );

// - - - - -

double angSep( double RA1,double DEC1,
	       double RA2,double DEC2, double  scale);

// random-number generators.
// May 2014: snran1 -> Flatran1,  float rangen -> double FlatRan
void   init_random_seed(int ISEED, int NSTREAM);
void   fill_RANLISTs(void);
void   sumstat_RANLISTs(int FLAG);
// xxx double unix_random(int istream) ;
// xxx double unix_GaussRan(int istream);

double unix_getRan_Flat1(int istream) ;
double unix_getRan_Gauss(int istream);

double getRan_Flat(int ilist, double *range);  //return rnmd on range[0-1]
double getRan_Flat1(int ilist);          // return 0 < random  < 1
double getRan_Gauss(int ilist);   // return Gauss randon (sigma=1)
double getRan_GaussClip(int ilist, double ranGmin, double ranGmax);
double getRan_GaussAsym(double siglo, double sighi, double peakinterval);
int    getRan_Poisson(double mean);

// mangled functions for fortran
double unix_random__(int *istream) ;
double getran_flat1__(int *ilist) ;          // for fortran
double getran_gauss__(int *ilist);         // for fortran

double biGaussRan_LEGACY(double siglo, double sighi);

double getRan_skewGauss(double rmin, double rmax, double siglo, double sighi,
			double skewlo, double skewhi);

double funVal_skewGauss(double x, double siglo,double sighi,
			double skewlo, double skewhi) ;

void   init_GaussIntegral(void);
double GaussIntegral(double nsig1, double nsig2);


// ------ sorting --------

void sortDouble(int NSORT, double *ARRAY, int ORDER,
		int *INDEX_SORT) ;
void sortFloat(int NSORT, float *ARRAY, int ORDER,
	       int *INDEX_SORT) ;
void sortInt(int NSORT, int *ARRAY, int ORDER,
	     int *INDEX_SORT) ;
void sortLong(int NSORT, long long *ARRAY, int ORDER,
	     int *INDEX_SORT) ;
void reverse_INDEX_SORT(int NSORT, int *INDEX_SORT) ;


// mangled fortran functions
void sortdouble_(int *NSORT, double *ARRAY, int *ORDER,
		 int *INDEX_SORT);
void sortfloat_(int *NSORT, float *ARRAY, int *ORDER,
		int *INDEX_SORT);
void sortint_(int *NSORT, int *ARRAY, int *ORDER,
	      int *INDEX_SORT);

// invert matrix to replace CERNLIB functions
void invertMatrix (int  N, int  n, double *Matrix ) ;
void invertmatrix_(int *N, int *n, double *Matrix ) ;


// functions for user-define PATH_SNDATA_SIM
#define MXPATH_SNDATA_SIM 20
void add_PATH_SNDATA_SIM(char *PATH);
int  getList_PATH_SNDATA_SIM(char **pathList);
FILE *openFile_PATH_SNDATA_SIM(char *mode);

// generic LC width util (for simulation)
void   init_lightCurveWidth(void);
double get_lightCurveWidth(int OPTMASK, int NOBS, double *TLIST,
			   double *MAGLIST, int *ERRFLAG, char *FUNCALL ) ;

void   init_lightcurvewidth__(void);
double get_lightcurvewidth__(int *OPTMASK, int *NOBS, double *TLIST,
			   double *MAGLIST, int *ERRFLAG, char *FUNCALL ) ;

void init_obs_atFLUXMAX(int OPTMASK, double *PARLIST, int VBOSE);
void get_obs_atFLUXMAX(char *CCID, int NOBS, float *FLUX, float *FLUXERR,
		       double *MJD, int *IFILTOBS, int *EP_atFLUXMAX);

void init_obs_atfluxmax__(int *OPTMASK, double *PARLIST, int *VBOSE);

void get_obs_atfluxmax__(char *CCID, int *NOBS, float *FLUX, float *FLUXERR,
			 double *MJD, int *IFILTOBS, int *EP_atFLUXMAX);

// glob functions
int  glob_file_list(char *wildcard, char ***file_list); // underlying util
int  store_glob_file_list(char *wildcard) ; 
int  store_glob_file_list__(char *wildcard) ;
void fetch_glob_file(int lang_flag, int ifile, char *file_name);
void fetch_glob_file__(int *lang_flag, int *ifile, char *file_name); 
void reset_glob_file_list(void);
void reset_glob_file_list__(void);

// multi-D malloc functions (copied from SALT2mu.c, July 2021)

void  print_debug_malloc(int opt, char *comment);
float malloc_double2D(int opt, int LEN1, int LEN2, double ***array2D );
float malloc_double3D(int opt, int LEN1, int LEN2, int LEN3,
                      double ****array3D );
float malloc_double4D(int opt, int LEN1, int LEN2, int LEN3, int LEN4,
                      double *****array4D );
float malloc_float3D(int opt, int LEN1, int LEN2, int LEN3,
                     float ****array3D );

float malloc_shortint2D(int opt, int LEN1, int LEN2,
			short int ***array2D );
float malloc_shortint4D(int opt, int LEN1, int LEN2, int LEN3, int LEN4,
			short int *****array4D );

// ============== END OF FILE =============
