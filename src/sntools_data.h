/******************************************
  Created Dec 10 2021 by R.Kessler

  Declarations for Data-handling tools that don't depend on format.

 ******************************************/

int REFAC_DATA_FLAG ;  // Apr 3 2026: refactor to store LOGMASS[_ERR] in z-bins

#define MXFILE_OVERRIDE   10
#define IVARMAX_OVERRIDE  200  // max IVAR;  can be large even with only 1 override var

struct {
  bool USE;
  int NFILE; // number of override files
  int NVAR;  // number of override variables
  int NTOT_PER_VAR[IVARMAX_OVERRIDE] ; // tot overrides per IVAR
  int NRD_PER_VAR[IVARMAX_OVERRIDE];   // NRD per IVAR per event

  int NMATCH_by_CID, NMATCH_by_GALID;
  char VARNAME_MATCH[40];  // e.g., CID, SNID, GALID ...

  // logicals to decide if zCMB or zHEL needs to be recomputed.
  int  IVAR_zCMB, IVAR_zHEL;
  int  IVAR_HOSTGAL_ZPHOT[MXHOSTGAL], IVAR_HOSTGAL_ZPHOT_ERR[MXHOSTGAL];
  int  IVAR_HOSTGALz_ZPHOT_QUANTILE[MXHOSTGAL] ;

  int  NZPHOT_Q ; // number of ZPHOT_Q[nnn] quantiles (May 2023) LEGACY
  char **VARLIST_ZPHOT_Q; // LEGACY

  int IVAR_NAME_IAUC, IVAR_NAME_TRANSIENT; 

  // define variables to monitor what variable(s) are actually over-written.
  // If no vars are over-written, abort with error message 
  int  NEVT;     // number of events with unique CID or GALID
  int  NVAR_USE; // number of override variables used (Aug 2025)

  char ID_LAST[40]; // last CID or GALID passed; used to count NUSE_OVERRIDE

  int   N_HOSTGAL_PHOTOZ_REPLACE;

  // for special REDSHIFT_FINAL update
  float ORIG_HOSTGAL_PHOTOZ[MXHOSTGAL], ORIG_HOSTGAL_PHOTOZ_ERR[MXHOSTGAL]; 
  float ORIG_HOSTGALz_ZPHOT_QUANTILE[MXHOSTGAL][MXBIN_HOSTGALz_QUANTILE];

  bool IS_SIM;
} RD_OVERRIDE;


#define FORMAT_SNDATA_FITS 32
#define FORMAT_SNDATA_TEXT  2

int FORMAT_SNDATA_READ ;
int FORMAT_SNDATA_WRITE ;

// ======== function prototypes =============

void SET_REFAC_DATA_FLAG(int refac_data_flag);

void copy_SNDATA_GLOBAL(int copyFlag, char *key,
                        int NVAL, char *stringVal, double *parVal);
void copy_SNDATA_HEAD(int copyFlag, char *key,
                      int NVAL, char *stringVal, double *parVal);
void copy_SNDATA_OBS(int copyFlag, char *key,
                     int NVAL,char *stringVal, double *parVal);
int  select_MJD_SNDATA(double *CUTWIN_MJD);
void host_property_list_sndata(char *HOST_PROPERTY_LIST);
void LOAD_VARNAME_ZPHOT_Q_LEGACY(char *PREFIX, int PCT, char *VARNAME) ; 

void copy_GENSPEC(int copyFlag, char *key, int ispec, double *parVal);

void copy_int(int copyFlag, double *DVAL0, int    *IVAL1) ;
void copy_lli(int copyFlag, double *DVAL0, long long  *IVAL1) ;
void copy_flt(int copyFlag, double *DVAL0, float  *FVAL1) ;
void copy_dbl(int copyFlag, double *DVAL0, double *DVAL1) ;
void copy_str(int copyFlag, char   *STR0,  char   *STR1 );
void copy_HOSTGALz(int copyFlag, char *PREFIX, HOSTGALz_DEF *HOSTGALz);

bool IS_SIMKEY_SNDATA(char *key);

void RD_OVERRIDE_INIT(char *OVERRIDE_FILE, int REQUIRE_DOCANA);
int  RD_OVERRIDE_FETCH(char *CCID, long long int GALID, char *VARNAME, double *DVAL, char *STRVAL);
void RD_OVERRIDE_STORE_ORIG(int IVAR);
void RD_OVERRIDE_POSTPROC(void); // special updates for redshift variables

void get_override_file_list(char *OVERRIDE_PATH, char *OVERRIDE_FILE_LIST);
void rd_override_append(void);
void rd_override_zspec(void);
void rd_override_zphot(int igal);
void rd_override_zphot_legacy(void);
void rd_override_name(void);

void rd_override_check_mistake(char *varname_mistake, char *varname_correct);

void RD_PRIVATE_INIT(char *PRIVATE_VARNAME_LIST); 

// mangled functions for fortran
void set_refac_data_flag__(int *refac_data_flag);

void copy_sndata_global__(int *copyFlag, char *key,
                          int *NVAL, char *stringVal, double *parVal);
void copy_sndata_head__(int *copyFlag, char *key,
                        int *NVAL, char *stringVal, double *parVal);
void copy_sndata_obs__(int *copyFlag, char *key,
                       int *NVAL,char *stringVal, double *parVal);
int  select_mjd_sndata__(double *MJD_WINDOW);
void host_property_list_sndata__(char *HOST_PROPERTY_LIST);

void copy_genspec__(int *copyFlag, char *key, int *ispec, double *parVal ) ;

void rd_override_init__(char *OVERRIDE_FILE, int *REQUIRE_DOCANA);


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@  LEGACY_QUANTILE_FUNCTIONS @@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void rd_override_zphot_legacy(void);
void rd_override_zphot_q_legacy(int OPT);
void LOAD_VARNAME_ZPHOT_Q_LEGACY(char *PREFIX, int PCT, char *VARNAME) ; 

// end:
