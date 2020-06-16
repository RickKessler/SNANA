/********************************
  Created Jun 12 2020 by R.Kessler

 *******************************/

#define MXMAP_GENPDF    8  // max number of maps per file
#define MXVAR_GENPDF    6  // max varnames per mape
#define MXROW_GENPDF 5000   // max number of rows per map

#define GENPDF_OPTMASK_EXTRAP 1

int      NMAP_GENPDF;
int      NCALL_GENPDF ;
int      OPT_EXTRAP_GENPDF;

struct {
  char     MAPNAME[40];
  char     *VARNAMES[MXVAR_GENPDF];
  GRIDMAP  GRIDMAP ;
  int      IVAR_HOSTLIB[MXVAR_GENPDF];
} GENPDF[MXMAP_GENPDF] ;



// ----------------------------------
//    FUNCTION PROTOTYPES
// ----------------------------------

void   init_genPDF(int OPTMASK, char *fileName, char *ignore);
void   assign_VARNAME_GENPDF(int imap, int ivar, char *varName) ;
double get_random_genPDF(char *parName, GENGAUSS_ASYM_DEF *GENGAUSS);
int IDMAP_GENPDF(char *parName);

// END
