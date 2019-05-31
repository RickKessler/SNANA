// Created Sep 2018

// define pre-processor command to use python interface
#define USE_PYTHONxxx 

// ===========================================
// global variables

#define MXLAM_BYOSED  10000 // max wave bins to define SED
#define MXPAR_BYOSED  10     // max number of params to describe SED

struct {
  int NLAM;           // number of bins in SED
  double *LAM, *SED;  // SED
  int LAST_EXTERNAL_ID ;

  int    NPAR ;
  char   **PARNAME;    // par names set during init stage
  double *PARVAL;      // par values update for each SED
} Event_BYOSED ;


// ===========================================
// function declarations
void init_genmag_BYOSED(char *PATH_VERSION, int OPTMASK, char *ARGLIST, 
			char *NAMES_HOSTPAR ) ;

void genmag_BYOSED(int EXTERNAL_ID, double zHEL, double MU, 
		   double MWEBV, int NHOSTPAR, double *HOSTPAR_LIST,
		   int IFILT, int NOBS, double *TOBS_list, 
		   double *MAGOBS_list, double *MAGERR_list );

int  fetchParNames_BYOSED(char **parNameList);
void fetchParVal_BYOSED(double *parVal);
void fetchSED_BYOSED(int EXTERNAL_ID, int NEWEVT_FLAG, double Tobs, int MXLAM, 
		     double *HOSTPAR_LIST, int *NLAM, double *LAM, double *FLUX);

void INTEG_zSED_BYOSED(int OPT_SPEC, int IFILT_OBS, double Tobs, 
		       double zHEL, double x0,
		       double RV, double AV,
		       int NLAM, double *LAM, double *SED,
		       double *Finteg, double *Fspec );

// return spectrum for spectrograph
void genSpec_BYOSED(double Tobs, double z, double MU, double MWEBV,
		    double RV_host, double AV_host, // (I)		     
		    double *GENFLUX_LIST,         // (O) fluxGen per bin 
		    double *GENMAG_LIST );        // (O) magGen per bin

void read_SALT2_template0(void); // for debug only (no python)

