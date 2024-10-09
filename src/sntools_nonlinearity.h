

// -------------------------------
#define MXBIN_NONLIN 100 
#define OPTMASK_NONLIN_COUNT_TOT  1  // nonlin from total counts (e-)
#define OPTMASK_NONLIN_COUNT_RATE 2  // nonlin from count rate (e-/sec)

char MODELNAME_NONLIN[100];
int  NMAP_NONLIN ;
int  OPTMASK_NONLIN;
int  DUMPFLAG_NONLIN ;

typedef struct {
  char   FILTERS[MXFILTINDX];
  int    MAPSIZE ;
  double MAPVAL[2][MXBIN_NONLIN];  // definition depends on MODEL
} NONLIN_DEF ;


NONLIN_DEF  *NONLIN_MAP ;

struct {
  int   NLINE;
  char  LINE[10][MXPATHLEN];
} NONLIN_README ;


// --------------------------------------
void   INIT_NONLIN(char *inFile) ;
void   init_nonlin__(char *inFile);

double GET_NONLIN(char *cfilt, double Fpe_source, double Fpe_sky, 
		  double genmag) ;

double get_nonlin__(char *cfilt, double *Fpe_source, double *Fpe_sky, 
		     double *genmag);
