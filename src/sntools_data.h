/******************************************
  Created Dec 10 2021 by R.Kessler

  Declarations for Data-handling tools that don't depend on format.

 ******************************************/

#define MXFILE_OVERRIDE 10
#define MXVAR_OVERRIDE  20

struct {
  bool USE;
  int NFILE; // number of override files
  int NVAR;  // number of override variables
  int N_PER_VAR[MXVAR_OVERRIDE] ;

  // logicals to decide if zCMB or zHEL needs to be recomputed.
  bool FOUND_zCMB, FOUND_zHEL;

} RD_OVERRIDE;


#define FORMAT_SNDATA_FITS 32
#define FORMAT_SNDATA_TEXT  2

int FORMAT_SNDATA_READ ;
int FORMAT_SNDATA_WRITE ;

// ======== function prototypes =============

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

bool IS_SIMKEY_SNDATA(char *key);

void RD_OVERRIDE_INIT(char *OVERRIDE_FILE);
int  RD_OVERRIDE_FETCH(char *CCID, char *VARNAME, double *DVAL);
void RD_OVERRIDE_POSTPROC(void); 
void rd_override_append(void);
void rd_override_zcalc(void);

// mangled functions for fortran
void copy_sndata_global__(int *copyFlag, char *key,
                          int *NVAL, char *stringVal,double *parVal);
void copy_sndata_head__(int *copyFlag, char *key,
                        int *NVAL, char *stringVal,double *parVal);
void copy_sndata_obs__(int *copyFlag, char *key,
                       int *NVAL,char *stringVal,double *parVal);
int  select_mjd_sndata__(double *MJD_WINDOW);

void copy_genspec__(int *copyFlag, char *key, int *ispec, double *parVal ) ;

void rd_override_init__(char *OVERRIDE_FILE);

// end 
