/********************************
  Created Jun 12 2020 by R.Kessler

 *******************************/

#define MXMAP_GENPDF    8  // max number of maps per file
#define MXVAR_GENPDF    6  // max varnames per mape
#define MXROW_GENPDF 15000   // max number of rows per map
#define MXITER_GENPDF  200   // max number of iterations for genPDF

#define  OPTMASK_GENPDF_EXTRAP       1
#define  OPTMASK_GENPDF_EXTERNAL_FP  8

int      NMAP_GENPDF;
int      NCALL_GENPDF ;
int      OPT_EXTRAP_GENPDF;

struct {
  char     MAPNAME[40];
  char     *VARNAMES[MXVAR_GENPDF];
  GRIDMAP  GRIDMAP ;
  int      IVAR_HOSTLIB[MXVAR_GENPDF];
} GENPDF[MXMAP_GENPDF] ;

float TMPSTORE_PROB_REF_GENPDF[MXITER_GENPDF];
float TMPSTORE_PROB_GENPDF[MXITER_GENPDF];
float TMPSTORE_RAN_GENPDF[MXITER_GENPDF];

// ----------------------------------
//    FUNCTION PROTOTYPES
// ----------------------------------

void   init_genPDF(int OPTMASK, FILE *fp, char *fileName, char *ignore);
void   assign_VARNAME_GENPDF(int imap, int ivar, char *varName) ;
double get_random_genPDF(char *parName, GENGAUSS_ASYM_DEF *GENGAUSS);
int IDMAP_GENPDF(char *parName);

// END
