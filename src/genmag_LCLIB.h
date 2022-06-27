
// ROW refers to row in LCLIB file that starts with 'T:'  or 'S:'
// OBS refers to real observation 

#define MXPAR_LCLIB 40            // max number of baggage params 
#define MXCUT_LCLIB 10
#define MXROW_LCLIB 1000000       // max number of DAY in LCLIB event.
#define MXOBS_TEMPLATE_LCLIB 20   // max number of template epochs to average
#define I2FLOAT_LCLIB 1000.0  // I2MAG = 1000 * MAG
#define MAGBRIGHT_LCLIB   5.0    // 
#define MAGFAINT_LCLIB   99.01   // allow mag=99 for zero flux
#define IFLAG_TEMPLATE_LCLIB 1
#define IFLAG_SEARCH_LCLIB   2

#define IFLAG_RECUR_PERIODIC   1  // recurring, but periodic (e.g., RRlyr)
#define IFLAG_RECUR_NONPERIODIC 2  // recurring, non-periodoc (e.g., AGN)
#define IFLAG_RECUR_NONRECUR   3  // non-recurring (e.g., SN)

#define OPTMASK_LCLIB_IGNORE_ANGLEMATCH 1 // option to ignore ANGLEMATCH cut
#define OPTMASK_LCLIB_useRADEC          8 // use RA,DEC from LCLIB
#define OPTMASK_LCLIB_DEBUG           512 // debug/refactor

#define DAYBACK_TEMPLATE_LCLIB 30.0 // used in forceTemplateRows
#define MODEL_RANMAG_LCLIB  "RANMAG" 

#define PARNAME_REDSHIFT_LCLIB  "REDSHIFT"
#define PARNAME_MWEBV_LCLIB     "MWEBV"

int LDUMP_EVENT_LCLIB ;

// info from global header in LCLIB_INFO struct
struct {
  char FILENAME[MXPATHLEN];
  FILE *FP;
  int  GZIPFLAG;

  char SURVEY[40];
  char FILTERS[MXFILTINDX];
  int  NFILTERS;

  int  NGENTOT ;      // total number of sim events to generate
  int  NEVENT ;       // expected number of LCLIB events in file
  int  NREPEAT ;      // Nrepeat each LCLIB event (for faster sim)
  int  OPTMASK ;      // user options

  bool   DO_ANGLEMATCH ;
  double TOBS_RANGE_MAX ; // max TOBS_MAX-TOBS_MIN (passed from main)

  int    IFLAG_RECUR_CLASS ;
  char   STRING_RECUR_CLASS[40];

  char   NAME_MODEL[60];
  int    NPAR_MODEL;       // number of model parameters
  int    NPAR_MODEL_STORE; // includes EVENTID and template mags

  char   PARNAME_MODEL[MXPAR_LCLIB][40] ; // list of modelPar names
  int    IPAR_REDSHIFT;    // non-zero for PARNAME=REDSHIFT
  int    IPAR_MWEBV ;      // non-zero of MWEBV is a model param
  double ZPHOTZ1ERR ;        // if REDSHIFT, Err[zphot/(1+z)] = 0.05

  // variables needed for setting photo-z
  double REDSHIFT_RANGE[2]; // read from LCLIB header
  char HOSTLIB_MSKOPT;      // store to reset after init_genmag_LCLIB()

  // template epoch info
  int    NEP_TEMPLATE;
  double EPLIST_TEMPLATE[MXOBS_TEMPLATE_LCLIB];
  double EPRANGE_TEMPLATE ;

  // debug info option
  int    DEBUGFLAG_RANMAG ;
  double GENRANGE_RANMAG[2];
  double GENRANGE_DIFMAG[2];

} LCLIB_INFO ;


// user cuts on PARVAL to reject LCLIB events (Oct 31 2017)
struct {
  int    NCUTWIN;
  char   PARNAME[MXCUT_LCLIB][40];
  double CUTWIN[MXCUT_LCLIB][2];
  int    ICUT_REDSHIFT;
} LCLIB_CUTS ;


// user inputs 
struct {
  double DAYSCALE ;             // DAY *= DAYSCALE (for debug only)
  double TOBS_OFFSET_RANGE[2] ; // if not zero, fix TOBS_OFFSET to this value
  int    ZERO_TEMPLATE_FLUX ;   // zero template flux if set
  int    FORCE_NREPEAT ;        // if=1, equal wgt per event (ignore dT_event)
} LCLIB_DEBUG ;


int ifiltmap_LCLIB[MXFILTINDX];  // ifilt = ifiltmap_LCLIB[ifilt_obs]

// store light curve for one evnt
struct {

  int NEVENT_READ;   // total number of events read from LCLIB
  int NREPEAT ;      // counter repeating each event.

  int NROW;             // combined number of S: and T: rows in event
  int NROW_T, NROW_S;
  int FIRSTROW_S, FIRSTROW_T ;
  int LASTROW_S,  LASTROW_T ;

  double FIRSTMAG[MXFILTINDX];
  double LASTMAG[MXFILTINDX];

  long long  ID ;      // ID after START_EVENT: key
  double RA, DEC, GLAT, GLON;  // RA, DEC, l, b
  double MWEBV;  // Feb 2021
  double PARVAL_MODEL[MXPAR_LCLIB] ;    // list of model par values

  double REDSHIFT;  // only if one of the parval names is 'REDSHIFT'
  double ZPHOT, ZPHOTERR ; // associated photo-z

  // declare light curve pointers with emphasis on using
  // minimal memory in case of very large event
  double     *DAY ;
  short int **I2MAG ;        // mag[ifilt][row]
  char      *STRING_FLAG ;   // S or T vs. row

  double DAYRANGE_S[2];  // min,max DAY for search
  double DAYRANGE_T[2];  // min,max DAY for template
  double DAYCOVER_S;     // DAYMAX-DAYMIN  
  double DAYCOVER_T;     // DAYMAX-DAYMIN  
  double DAYCOVER_ALL;   // idem for S & T
  double DAY_RANDOM ;    // random day to start in event

  int NROW_MAX, MEMLC_MAX; // keep track of max event size to inform user.

  int NforceTemplateRows ; 
  int LAST_EXTERNAL_ID ;

  double ANGLEMATCH, ANGLEMATCH_b;

  // computed quantities
  double PEAKMAG_S ;              // brightest mag
  double PEAKDAY_S ;              // LCLIB day at brightest mag
  double magTemplate[MXFILTINDX]; // mag in template

  double TGENRANGE_NONRECUR;     // range to generate non-recur events
  int    NROWADD_NONRECUR ;      // internally add S epochs with T mag.

  double TOBS_MIN, TOBS_MAX;      // passed to function set_TobsRange_LCLIB
  double TOBS_RANGE ;
  double TOBS_OFFSET ;            // shift Tobs overlap LCLIB DAY-range

  int NCHAR_ROW; // used by fseek to skip rows (Dec 2021)

} LCLIB_EVENT ;


// for Poisson generator (non-recurring)
//const gsl_rng_type *T_LCLIB;
//gsl_rng *r_LCLIB;


// ========= FUNCTION PROTOTYPES ================

// SNANA-sim wrappers
void init_genmag_LCLIB(char *INFILE, char *STRING_TEMPLATE_EPOCHS,
		       char *SURVEY, char *FILTERLIST, 
		       double TOBS_RANGE_MAX, int NGENTOT_EXPECT, 
		       int IFLAG_RANSTART, int OPTMASK );
void open_LCLIB(char *INFILE);
void read_GLOBAL_HEADER_LCLIB(void) ;
void parse_TEMPLATE_EPOCHS_LCLIB(char *STRING_TEMPLATE_EPOCHS);
void check_LCLIB(char *SURVEY, char *FILTERS) ;
void set_randomStart_LCLIB(void);

void parse_PARNAMES_LCLIB(char *parNameString);
void read_PARVAL_LCLIB(char *LINE);
void coord_translate_LCLIB(double *RA, double *DEC);
int  keep_PARVAL_LCLIB(void);
int  keep_ANGLEMATCH_LCLIB(double GalLat, double GalLong);

void addTemplateRows_LCLIB(void);
void addTemplateRows_PERIODIC(void) ;
void addTemplateRows_NONRECUR(void);
void addTemplateReset_NONRECUR(void);
void forceTemplateRows_LCLIB(void);
void forceTemplateReset_LCLIB(void);
void checkMag_PERIODIC_LCLIB(int ifilt);

void ranPhase_PERIODIC_LCLIB(void) ;

void dumpEvent_LCLIB(void);

int  set_IFLAG_RECUR_LCLIB(char *string_recur_class);
void set_TobsRange_LCLIB(double *TobsRange);
void set_NREPEAT_LCLIB(void);

void  set_REDSHIFT_LCLIB(void);

void genmag_LCLIB (int EXTERNAL_ID, int ifilt_obs, 
		   double *RA, double *DEC, double *mwebv, double lamFiltAvg,
		   int Nobs,  double *Tobs_list, double *magList_S,
		   double *mag_T, double *TobsPeak );

void readNext_LCLIB(double *RA, double *DEC);
void read_ROW_LCLIB(int IROW, char *KEY, char **ptrWDLIST); // read one row from LCLIB
void malloc_LCLIB_EVENT(int OPT);

void   set_TOBS_OFFSET_LCLIB(void) ;
void   get_TobsMinMax(int NOBS, double *TobsList, 
		      double *Tobs_min, double *Tobs_max);

double magSearch_LCLIB(int ifilt, double Tobs) ;
double magTemplate_LCLIB(int EXTERNAL_ID, int ifilt);
void   store_magTemplate_LCLIB(int EXTERNAL_ID, int ifilt, double XT_MW) ;

double magInterp_LCLIB(double T, int NROW, double *DAYLIST, short int *I2MAG);

// from snlc_sim:
double gen_MWEBV(double RA, double DEC);

