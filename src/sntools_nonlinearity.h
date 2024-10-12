

// -------------------------------
#define MXBIN_NONLIN 100 
#define OPTMASK_NONLIN_COUNT_TOT  1  // nonlin from total counts (e-)
#define OPTMASK_NONLIN_COUNT_RATE 2  // nonlin from count rate (e-/sec)
#define OPTMASK_NONLIN_PER_NEA  4  // nonlin map provided for total flux 
#define OPTMASK_NONLIN_PER_PIX  8  // nonlin map provided for each pixel (will convert to total flux)
#define OPTMASK_NONLIN_DUMPFLAG   1024  // dump flag
#define OPTMASK_NONLIN_DEBUGFLAG  2048  // internal test/debug flag


char MODELNAME_NONLIN[100];
int  NMAP_NONLIN ;
int  OPTMASK_NONLIN;
int  DUMPFLAG_NONLIN ;
int  DEBUGFLAG_NONLIN ;

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

double GET_NONLIN(char *CCID, char *cfilt, double Texpose, double NEA, double *Fpe_list, 
		  double genmag) ;

double get_nonlin__(char *CCID, char *cfilt, double *Texpose, double *NEA, double *Fpe_list,
		    double *genmag);

void check_OPTMASK_NONLIN(void) ;
double get_flux_scale_NONLIN(char *cfilt, double flux);

