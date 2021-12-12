/******************************************
  Created Dec 10 2021 by R.Kessler

  Declarations for Data-handling tools that don't depend on format.

 ******************************************/

#define MXFILE_OVERRIDE 10
struct {
  bool USE;
  int N_FILE; // number of override files
  
} RD_OVERRIDE;

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


void RD_OVERRIDE_INIT(char *OVERRIDE_FILE);
int  RD_OVERRIDE_FETCH(char *CCID, char *VARNAME, double *DVAL);

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
