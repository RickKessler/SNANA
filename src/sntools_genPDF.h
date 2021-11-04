/********************************
  Created Jun 12 2020 by R.Kessler

  Oct 22 2020: MXITER_GENPDF -> 1000 (was 200)

 *******************************/

#define MXMAP_GENPDF    8  // max number of maps per file
#define MXVAR_GENPDF   10  // max varnames per mape
#define MXROW_GENPDF 15000   // max number of rows per map
#define MXITER_GENPDF  1000   // max number of iterations for genPDF
#define PROBMAX_REJECT_GENPDF 1.0E-3  // ignore range where P > this value

#define  OPTMASK_GENPDF_EXTRAP        1
#define  OPTMASK_GENPDF_SLOW          2  // use full val range
#define  OPTMASK_GENPDF_KEYSOURCE_ARG 4  // arg is from command line
#define  OPTMASK_GENPDF_EXTERNAL_FP   8

int      NMAP_GENPDF;
int      NCALL_GENPDF ;
int      OPT_EXTRAP_GENPDF;
int      OPTMASK_GENPDF ;
int      KEYSOURCE_GENPDF; // 1=file, 2=arg; used for prioritization

struct {
  char     MAPNAME[40];
  char     *VARNAMES[MXVAR_GENPDF];
  int      NVAR;              // May 26 2021
  GRIDMAP  GRIDMAP ;
  int      IVAR_HOSTLIB[MXVAR_GENPDF];
  
  // track stats on number of iterations to find value
  int N_CALL ;
  int N_ITER_SUM,  N_ITER_MAX ;
  
} GENPDF[MXMAP_GENPDF] ;

float TMPSTORE_PROB_REF_GENPDF[MXITER_GENPDF];
float TMPSTORE_PROB_GENPDF[MXITER_GENPDF];
float TMPSTORE_RAN_GENPDF[MXITER_GENPDF];

// ----------------------------------
//    FUNCTION PROTOTYPES
// ----------------------------------

void   init_genPDF(int OPTMASK, FILE *fp, char *fileName, char *ignore);
void   init_genPDF_from_GenGauss(int IMAP, GENGAUSS_ASYM_DEF *GENGAUSS);
void   assign_VARNAME_GENPDF(int imap, int ivar, char *varName) ;
void   checkAbort_VARNAME_GENPDF(char *varName);
double get_random_genPDF(char *parName, GENGAUSS_ASYM_DEF *GENGAUSS);
void   get_VAL_RANGE_genPDF(int IDMAP, double *val_inputs, double *VAL_RANGE, int dumpFlag);
void   free_memory_genPDF(void); // release memory of all genPDF maps

int  IDMAP_GENPDF(char *parName, bool *LOGPARAM);
void iter_summary_genPDF(void);
bool matchVar_GENPDF_GENGAUSS(char *varName_GENPDF, char *varName_GENGAUSS);

// END
