/**************************************************
   sntools_trigger.h


   functions and Variables for maps related to
     + Detection efficiency
     + PHOTPROB/ML map
     + SPEC-confirmed efficiency
     + zHOST efficiency

  Jul 13 2020: MXMAP_SEARCHEFF_zHOST -> 20 (was 10)

  Feb 05 2021: define FIELDLIST_OVP and NFIELD_OVP

 **************************************************/


#define  IVARTYPE_SPECEFF_PEAKMAG  1  // flag val that var is a peakmag
#define  IVARTYPE_SPECEFF_COLOR    2  // flag val that var is a color
#define  IVARTYPE_SPECEFF_REDSHIFT 3  // flag val that var is a redshift
#define  IVARTYPE_SPECEFF_PEAKMJD  4  // flag val that var is a peakmjd
#define  IVARTYPE_SPECEFF_DTPEAK   5  // flag val that var is a dTpeak_min
#define  IVARTYPE_SPECEFF_HOSTMAG  6  // flag val that var is a HOSTmag
#define  IVARTYPE_SPECEFF_SBMAG    7  // flag val that var is a SBmag
#define  IVARTYPE_SPECEFF_SALT2mB   8  // flag val that var is a dTpeak_min
#define  IVARTYPE_SPECEFF          9  // flag val that var is SPECEFF
#define  MXVAR_SEARCHEFF_SPEC      10  // max number of SPEC-eff VARNAMES

#define  MXMAP_SEARCHEFF_DETECT   50  
#define  MXROW_SEARCHEFF_DETECT   10000

#define  MXMAP_SEARCHEFF_PHOTPROB     10
#define  MXROW_SEARCHEFF_PHOTPROB  10000
#define  MXVAR_SEARCHEFF_PHOTPROB     20
#define  MXOBS_PHOTPROB              400 // max obs per event with PHOTPROB

#define  MXMAP_SEARCHEFF_SPEC   20   // max number of SPEC-maps
#define  MXROW_SEARCHEFF_SPEC 30000  // temp max size of SPEC-eff map

#define  MXMAP_SEARCHEFF_zHOST   20    // max number of zHOST maps
#define  MXROW_SEARCHEFF_zHOST  5000   // max size of each map
#define  MXVAR_SEARCHEFF_zHOST   10    // max number of zHOST VARNAMES
#define  IVERSION_zHOST_LEGACY    1    // legacy z-only map
#define  IVERSION_zHOST_MULTID    2    // multi-D map of HoSTLIB properties

#define FLAG_EFFSNR_DETECT    1  // flag for EFF vs, SNR
#define FLAG_EFFABSSNR_DETECT 2  // flag for EFF vs, ABS(SNR)
#define FLAG_EFFMAG_DETECT    3  // flat for EFF vs. MAG
#define SEARCHEFF_PARNAME (char*[4]){ "", "SNR", "ABS(SNR)", "MAG" }

#define APPLYMASK_SEARCHEFF_PIPELINE    1  // pipeline detection
#define APPLYMASK_SEARCHEFF_SPEC        2  // spec confirmed
#define APPLYMASK_SEARCHEFF_zHOST       4  // zSpec from host

#define DETECT_MASK_SNR         1  // detect mask for SNR or MAG
#define DETECT_MASK_MJD_TRIGGER 2  // identify obs where trigger passes
#define DETECT_MASK_PHOTPROB    4  // detect mask for PHOTPROB

#define MXOBS_TRIGGER 10*MXEPOCH 

// define allowed variables to use in PHOTPROB map
#define IVARABS_PHOTPROB_SNR       0
#define IVARABS_PHOTPROB_LOGSNR    1
#define IVARABS_PHOTPROB_SBMAG     2
#define IVARABS_PHOTPROB_REDSHIFT  3  // Feb 14 2020
#define IVARABS_PHOTPROB_GALMAG    4  // Feb 16 2020
#define MXDEF_VARNAMES_PHOTPROB    5

char COMMENT_README_SEARCHEFF[200];

char VARDEF_SEARCHEFF_PHOTPROB[MXDEF_VARNAMES_PHOTPROB][20] ;

int  SEARCHEFF_FLAG;   // FLAG_EFFSNR or FLAG_EFFMAG

char PATH_SEARCHEFF[2*MXPATHLEN]; // path for PIPELINE_FILE, SPEC_FILE

int NONZERO_SEARCHEFF_SPEC  ;
int NONZERO_SEARCHEFF_zHOST ;
double **SEARCHEFF_TMPMAP2D ; // generic array to read/store map

struct  {
  char   USER_PIPELINE_LOGIC_FILE[MXPATHLEN]; // added Dec 2015
  char   USER_PIPELINE_EFF_FILE[MXPATHLEN]; // EFF vs. SNR or vs. MAG
  char   USER_SPEC_FILE[MXPATHLEN];     // SPEC eff vs. peakmag,z,etc ...
  char   USER_zHOST_FILE[MXPATHLEN] ;   // eff vs. z of getting zHOST

  // repeat input file names with full path 
  char   PIPELINE_LOGIC_FILE[MXPATHLEN]; 
  char   PIPELINE_EFF_FILE[MXPATHLEN]; 
  char   SPEC_FILE[MXPATHLEN];     
  char   zHOST_FILE[MXPATHLEN] ;   

  int    IFLAG_zHOST_EFFZERO;      // flag to set EFF_zHOST=0
  int    IFLAG_SPEC_EFFZERO;       // flag to set EFF_SPEC=0
  int    IVERSION_zHOST;           // 1=legacy, 2=multi-D

  double USER_SPECEFF_SCALE; // default=1.0
  double CUTWIN_SNRMAX_zHOST[2]; // extra requirement for zHOST

  // time-window (days) in which all detections count as 1 detection
  // e.g.,  0.3 -> 1 roughly night, 0.007 -> 10 minutes, or SDSS ugriz
  double TIME_SINGLE_DETECT ;  // days

  // global peakmag shift for GRIDMAP lookup
  double MAGSHIFT_SPECEFF ; 

  // min number of observations to evaluate trigger
  int    MINOBS ; 

    // define number of maps, which is also a logical flag.
  int NMAP_DETECT, NMAP_PHOTPROB, NMAP_SPEC, NMAP_zHOST ;

  int NREDUCED_CORR_PHOTPROB ;  // Nmap(PHOTPROB) with non-zero REDUCED_CORR
  int NPHOTPROB_DUMP;           // number of PHOTPROB>0 to dump per event

  int LEGACYMAP_DETECT; // no MAPNAME_DETECT of ENDMAP (MAr 2018)

  int  PHOTFLAG_DETECT ;  // optional PHOTFLAG bit for each detection
  int  PHOTFLAG_TRIGGER ; // idem, for the one epoch making trigger

  // - - - - - -  DEBUG STUFF - - - - - - - 
  // allow for fixing pipeline efficiency for debugging
  double FIX_EFF_PIPELINE ; 

  // option to hack z-dependent efficiency function (for debugging)
  int    FUNEFF_DEBUG ;  

  int  OPTMASK_OPENFILE ;

} INPUTS_SEARCHEFF ;


// Define structure for serach eff vs. SNR for each filter
int    MAPVERSION_SEARCHEFF_DETECT ; // allows legacy or new map style
struct SEARCHEFF_PIPELINE {
  char   MAPNAME[40] ;
  int    NBIN ;
  char   FILTERLIST[MXFILTINDX] ;
  double *VAL, *EFF ;
  int    NLINE_README;
  char   README[20][MXPATHLEN];
} SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT+1] ;


int MAPVERSION_SEARCHEFF_PHOTPROB ;

struct {
  char NAME[40] ;
  char FIELDLIST[60] ;
  char FILTERLIST[MXFILTINDX] ;
  char VARNAMES[MXVAR_SEARCHEFF_PHOTPROB][20] ;
  int  IVARABS[MXVAR_SEARCHEFF_PHOTPROB] ;
  double  VALMIN[MXVAR_SEARCHEFF_PHOTPROB] ;
  double  VALMAX[MXVAR_SEARCHEFF_PHOTPROB] ;
  double  BINSIZE; // photprob bin size
  int  NVAR_TOT, NVAR_MAP, NFUN_CDF ;
  int  NROW ;
  int  REQUIRE_DETECTION; // flag read from SEARCH-input file
  double CUTVAL ;         // reject epoch if PHOTPROB < CUTVAL (Feb 2020)
  double REDUCED_CORR ;   // reduced correlation (-1:1) read from input file
 
  double PHOTPROB_CDFBINS[MXVAR_SEARCHEFF_PHOTPROB] ;
  struct GRIDMAP GRIDMAP; 

  int NLINE_README;
  char README[4][MXPATHLEN];

} SEARCHEFF_PHOTPROB[MXMAP_SEARCHEFF_PHOTPROB];


// store obs-list for PHOTPROB; needed to compute correlated randoms.
struct {
  int    NSTORE;
  int    IMAP_LIST[MXOBS_PHOTPROB]; // IMAP for PHOTPROB
  int    OBS_LIST[MXOBS_PHOTPROB];  // obs index for SEARCHEFF_DATA
  double RAN_LIST[MXOBS_PHOTPROB];  // list of ran[0,1]
  int    OBSINV_LIST[MXOBS_TRIGGER];  // local index vs [obs]
} OBS_PHOTPROB;


#define MXMASK_SEARCHEFF_LOGIC 10 // max number of logic conditions
struct SEARCHEFF_LOGIC {
  int  NMJD;     // number of MJDs to have a detection
  int  NMASK;    // number of ORs
  int  IFILTDEF_MASK[MXMASK_SEARCHEFF_LOGIC];  // AND-mask vs. NMASK
  char INPUT_STRING[80];  // user-input logic-string
} SEARCHEFF_LOGIC ;


struct {
  char VARNAMES[MXVAR_SEARCHEFF_SPEC][40] ;
  //  int NROW, NVAR; 
  int IVAR, IVAR_REDSHIFT, IVAR_PEAKMJD, IVAR_DTPEAK, IVAR_SALT2mB;

  char FIELDLIST[60] ;
  int IVARTYPE[MXVAR_SEARCHEFF_SPEC] ;

  // ifilt_obs (mag and color) vs. IVAR index
  int NFILTLIST_PEAKMAG[MXVAR_SEARCHEFF_SPEC];
  int IFILTLIST_PEAKMAG[MXVAR_SEARCHEFF_SPEC][MXFILTINDX];  
  int IFILTOBS_PEAKCOLOR[MXVAR_SEARCHEFF_SPEC][2]; 
  int IFILTOBS_HOSTMAG[MXVAR_SEARCHEFF_SPEC]; 
  int IFILTOBS_SBMAG[MXVAR_SEARCHEFF_SPEC]; 

  struct GRIDMAP GRIDMAP ;

} SEARCHEFF_SPEC[MXMAP_SEARCHEFF_SPEC] ;


struct {
  int   IVARTYPE_MASK ;
  int   FLAG_PEAKMAG_ONLY ;
  int   BOOLEAN_OR, BOOLEAN_AND;
  int   NLINE_README;
  char  README[40][MXPATHLEN];
} SEARCHEFF_SPEC_INFO ;



struct {
  int  NVAR ;   // used by init_HOSTLIB
  char FIELDLIST[100];  // fieldList for each map
  double PEAKMJD_RANGE[2]; // PEAKMJD range for each map (7/2020)
  char VARNAMES_HOSTLIB[MXVAR_SEARCHEFF_zHOST][40] ; 
  int  IVAR_HOSTLIB[MXVAR_SEARCHEFF_zHOST] ; // points to HOSTLIB ivar
  struct GRIDMAP GRIDMAP ;
} SEARCHEFF_zHOST[MXMAP_SEARCHEFF_zHOST] ;


struct {
  char   FIELDLIST[100]; // fieldList for each map
  int    NROW ;
  double *REDSHIFT, *EFF;    // malloc for each new map
} SEARCHEFF_zHOST_LEGACY[MXMAP_SEARCHEFF_zHOST] ;


// Oct 2021 - definine MJDs associated with pipeline detections
typedef struct {
  double TRIGGER, FIRST, LAST;
} MJD_DETECT_DEF ;

// ------ data needed to evaluate trigger -------

struct {
  int    NOBS ;

  // scalars
  int    CID ;
  double REDSHIFT, PEAKMJD, DTPEAK_MIN, SALT2mB, SNRMAX, MWEBV ;
  int    SIMLIB_ID;


  char FIELDNAME[60]; // e.g., X3 or X1+X3 for overlap
  char FIELDLIST_OVP[MXFIELD_OVP][20]; //specify each ovp field 
  int  NFIELD_OVP;   // number of overlap fields

  // filter-dependent peak-mags
  double PEAKMAG[MXFILTINDX] ;
  double HOSTMAG[MXFILTINDX] ;
  double SBMAG[MXFILTINDX] ;

  // obs-dependent quantities
  double MJD[MXOBS_TRIGGER];
  double MAG[MXOBS_TRIGGER];  // mag for each obs
  double SNR[MXOBS_TRIGGER];  // signal-to-noise for each obs
  double FLUX[MXOBS_TRIGGER]; // flux in ADU, or errmsg only
  double FLUXERR[MXOBS_TRIGGER];
  int    IFILTOBS[MXOBS_TRIGGER];  // absolute filter index. each obs
  int    NPE_SAT[MXOBS_TRIGGER];   // Npe above sat (negative --> ok)
  int detectFlag[MXOBS_TRIGGER]; // detection flag for each obs (not trigger)
  double PHOTPROB[MXOBS_TRIGGER];  // Mar 13 2018
  int    NEXPOSE[MXOBS_TRIGGER]; // Number of Exposures in CO-ADD. 08/06/2022 
} SEARCHEFF_DATA ;


  // randoms
struct {
  double FLAT_PIPELINE[MXOBS_TRIGGER];     // flat ran for each obs [0,1]
  double FLAT_PHOTPROB[MXOBS_TRIGGER];     // flat ran for each detection
  double GAUSS_PHOTPROB[MXOBS_TRIGGER];     // Gauss ran for each detection
  double GAUSSCORR_PHOTPROB[MXOBS_TRIGGER] ;  // correlated Gauss Ran
  double FLAT_SPEC[MXFILTINDX+1] ;            // for each filter
} SEARCHEFF_RANDOMS ;


// ============== FUNCTION PROTOTYPES ===============

void   init_SEARCHEFF(char *SURVEY, int APPLYMASK_SEARCHEFF );
int    init_SEARCHEFF_PIPELINE(char *survey);
void   init_SEARCHEFF_LOGIC(char *survey) ;
void   init_SEARCHEFF_SPEC(char *survey)  ;
void   init_SEARCHEFF_zHOST(char *survey) ;
FILE   *open_zHOST_FILE(int OPT);
void   read_VARNAMES_zHOST(FILE *fp);
void   read_zHOST_FILE_LEGACY(FILE *fp);
void   read_zHOST_FILE(FILE *fp);
int    parse_VARNAMES_zHOST(FILE *fp, int *ivar_HOSTLIB, 
			    char **varName_HOSTLIB, char *varNameList );

int  readMap_SEARCHEFF_DETECT  (FILE *fp,  char *key);
int  readMap_SEARCHEFF_PHOTPROB(FILE *fp,  char *key);
int  malloc_NEXTMAP_SEARCHEFF_DETECT(void);

void   check_APPLYMASK_SEARCHEFF(char *SURVEY, int APPLYMASK_SEARCHEFF);

int    gen_SEARCHEFF(int ID, double *EFF_SPEC, double *EFF_zHOST, 
		     MJD_DETECT_DEF *MJD_DETECT );
int    gen_SEARCHEFF_PIPELINE(int ID, MJD_DETECT_DEF *MJD_DETECT );
int    gen_SEARCHEFF_SPEC(int ID, double *EFF_SPEC );
int    gen_SEARCHEFF_zHOST(int ID, double *EFF_zHOST );
int    gen_SEARCHEFF_DEBUG(char *what, double RAN, double *EFF);
double interp_SEARCHEFF_zHOST_LEGACY(void);
double interp_SEARCHEFF_zHOST(void);

void   check_SEARCHEFF_DETECT(int imap );
void   check_SEARCHEFF_PHOTPROB(int imap );
double LOAD_SPECEFF_VAR(int imap, int ivar);
void   LOAD_PHOTPROB_CDF(int NVAR_CDF, double *WGTLIST );
double LOAD_PHOTPROB_VAR(int OBS, int IMAP, int IVAR) ;
double GETEFF_PIPELINE_DETECT(int obs);

void   setObs_for_PHOTPROB(int DETECT_FLAG, int obs);
void   setRan_for_PHOTPROB(void) ;
double get_PIPELINE_PHOTPROB(int obs);
double get_PIPELINE_PHOTPROB_Obsolete(int DETECT_FLAG, int obs);
void   dumpLine_PIPELINE_PHOTPROB(void);

void   assign_SPECEFF(int imap, int ivar, char *VARNAME) ;
void   parse_search_eff_logic(char *survey, int NMJD, char *logic);
int    IFILTOBS_SPECEFF_VAR(char *VARNAME, char *PREFIX) ;
int    IVARABS_SEARCHEFF_PHOTPROB(char *VARNAME);

bool   MATCH_SEARCHEFF_FIELD(char *field_map);

// ============= END ===============
