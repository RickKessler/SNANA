/* =======================================
 Created July 24, 2012 by R.Kessler
========================================= */

// -----------------------
// Function Declarations
// -----------------------

void PSNID_INIT_VAR(void);
void PSNID_INIT_FILTERS(char *filtlist_fit);
void PSNID_READ_TEMPLATES(int TYPEINDX, char *file);
void PSNID_INIT_TEMPLATES(char *templates_nonIa_list, char *templates_nonIa_ignore) ;
void PSNID_INIT_XTMW(int OPT_MWCOLORLAW);
void PSNID_USER_INPUT(int NVAR, double *input_array, char *input_string );

void  psnid_init_nonIa_ignore(char *templates_nonIa_list, char *templates_nonIa_ignore) ;
int   psnid_match_NONIA_NAME(char *name);

void get_template_lc(int TYPEINDX, int iz, int ic1, int ic2, 
		     int ishape, int ifiltobs, int optDump,
		     int *NEPOCH, double *Trest,
		     double *mag, double *magerr) ;

void  psnid_check_filtlist_fit(int TYPEINDX) ;
void  psnid_init_nonIa_wgts(void) ;
void  psnid_dump_nonIa_templates(void) ;

void  psnid_dumpBins_SNGRID(void); 
void  psnid_dumpLine(int TYPEINDX, int IPAR, int ibin);
void  psnid_dumpLC_SNGRID(void); 
void  psnid_binCheck(int TYPEINDX, int IPAR, int ibin) ;


void psnid_store_data(char *CCID, int NOBS, int *IFILT, double *MJD,
		      double *FLUX, double *FLUXERR, double *FLUXSIM,
		      double *REDSHIFT, double *REDSHIFT_ERR,
		      double MWEBV, double MWEBVERR, int SIM_NON1A_INDEX );
void SNLCPAK_PSNID_DATA(double PKMJD);  // pass data to plot interface

void psnid_dumpInput_data(char *CCID, int NOBS, int *IFILTOBS, 
			  double *MJD, double *FLUX, double *FLUXERR,
			  double *REDSHIFT, double *REDSHIFT_ERR,
			  double MWEBV, double MWEBVERR);

void psnid_pogson2mJy(double mag, double mage,
		      double *flux, double *fluxe) ;

void psnid_pogson2fluxcal(double mag, double mage,
			  double *flux, double *fluxe) ;

void psnid_mktable_pogson2fluxcal(void);

// --------------------------------
// Global variable Declarations
// --------------------------------

#define PSNID_PI             3.14159265358979
#define MXTYPEINDX_PSNID     2
#define TYPEINDX_SNIA_PSNID  1  // for template struct
#define TYPEINDX_NONIA_PSNID 2  // for template struct
#define MXEP_PSNID      MXEPOCH // max epochs summed over filters


#define  PSNID_GOODMAG_LO       5.00     // changed from -8.99 (RK)
#define  PSNID_GOODMAG_HI      31.99
#define  PSNID_GOODMAGERR_LO    0.00 
#define  PSNID_GOODMAGERR_HI    4.99

#define  MXITER_PSNID 3  // should match MXITER_PNSNID in psnid.cat

// define SN grid structures using typedef SNGRID_DEF in sngridtools.h
SNGRID_DEF  SNGRID_PSNID[MXTYPEINDX_PSNID+1] ; // NULL, IA, NONIA


int  USEFLAG_TEMPLATES_PSNID[MXTYPEINDX_PSNID+1] ;   // logical flag
char TEMPLATES_FILE_PSNID[MXTYPEINDX_PSNID+1][200] ;    // vs. TYPEINDX
char TEMPLATETYPE_PSNID[MXTYPEINDX_PSNID+1][8] ;   // vs. TYPEINDX

char PATH_TEMPLATES_PSNID[200];

// define structure of inputs that are passed to psnid during init stage
struct  PSNID_INPUTS  {

  // start with variables passed directly from input file

  // &SNLCINP variables
  double H0, OMAT, OLAM, W0 ;

  // &PSNIDINP variables
  char CFILTLIST[80];  // char-list of filters; i.e, 'griz'
  char CFILTLIST_PEAKMAG_STORE[80];
  char OUTFILE[200];   // text-output for fit results

  char NMLFILE[200];    // fortran namelist file to allow reading

  char MODELNAME_MAGERR[60] ; // 

  double AV_TAU;       // exp(-av/AV_TAU) for av>0
  double AV_SMEAR ;    // exp(av/AV_SMEAR) for av<0
  double AV_PRIOR_STR; // strength of AV prior

  double WGT_ZPRIOR ;   // weight of external redshift prior

  int OPT_ZPRIOR ;         // bits 0,1,2 -> FLAT, ZSPEC, Zphot(host)
  double CUTWIN_ZERR[2] ;  // select redshift errro
  int MCMC_NSTEP ;         // number of MCMC steps

  double COLOR_MIN, COLOR_MAX;   // minimum and maximum color parameter
  int NCOLOR;                    // number of color param bins in grid search

  double DMU_MIN, DMU_MAX;       // minimum and maximum dmu (Ia and NON1A)
  double DMU_NON1A_MIN, DMU_NON1A_MAX; // NON1A only
  int NDMU;                      // number of dmu bins in grid search

  double TEMPLERR_SCALE_SNIA, TEMPLERR_SCALE_NONIA; // template error scale

  // bayesian prob-cuts for typing
  double PBAYES_IA_CUT, PBAYES_IBC_CUT, PBAYES_II_CUT, PBAYES_PEC1A_CUT,
    PBAYES_MODEL1_CUT, PBAYES_MODEL2_CUT,
    PBAYES_MODEL3_CUT, PBAYES_MODEL4_CUT ;
  
  // fitprob cuts for typing
  double FITPROB_IA_CUT, FITPROB_IBC_CUT, FITPROB_II_CUT, FITPROB_PEC1A_CUT,
    FITPROB_MODEL1_CUT, FITPROB_MODEL2_CUT,
    FITPROB_MODEL3_CUT, FITPROB_MODEL4_CUT ;
  
  int    OPT_RATEPRIOR ;
  double ZRATEPRIOR_SNIA[3] ;
  double ZRATEPRIOR_NONIA[3] ;

  int OPT_SIMCHEAT;        // 3/10/2017: 

  int NREJECT_OUTLIER;       // max number of outlier points to reject
  double CHISQMIN_OUTLIER;   // min chi2 for outlier rejection

  double TMAX_START[MXITER_PSNID];
  double TMAX_STOP[MXITER_PSNID];
  double TMAX_STEP[MXITER_PSNID];

  // quantities below are computed from the raw input above

  int  NFILT;                  // number of filters to fit
  int  USEFILT[MXFILTINDX];    // yes/no
  int  USEFILT_PEAKMAG_STORE[MXFILTINDX] ;
  int  USEFILT_GRIDOFF[MXTYPEINDX_PSNID+1][MXFILTINDX]; // SNGRID  offset
  int  IFILTLIST[80];      // integer-list of abs. filter indices
  int  IFILTLIST_INV[80];  // inverse map of above

  int  WRSTAT_OUTFILE_LEGACY;   // legacy flag to write outfile
  int  WRSTAT_TABLE;            // flag to make fitres table
  int  LSIM ;                   // =1 for simulation, 0 for data

  char SNANA_VERSION[20] ; // snana version
  // xxx mark delete Dec 10 2017  char VERSION_PHOTOMETRY[200];

  double XTMW_at_AV1[MXFILTINDX] ;  // Galactic extinctoin at AV=1

} PSNID_INPUTS ;




#define  PSNID_MXVAR_FITRES 100  // max number of output variables
struct PSNID_FITRES {
  int   NVAR ;
  char  VARNAMES[PSNID_MXVAR_FITRES][40] ;
  int   USE4NN[PSNID_MXVAR_FITRES] ;  // T => needed for NN analysis

  FILE *FP_OUTFILE ;  // pointer to legacy outfile option
} PSNID_FITRES ;


struct PSNID_NEARNBR {
  int NVAR ;
} PSNID_NEARNBR ;


// dump variables
#define MXLCDUMP_PSNID 10  // avoid huge accidental dumps
int NLCDUMP_PSNID ;


int FLAG_PSNID_INIT_VAR ; // to ensure that PSNID_INIT_VAR is called

// define data structure filled by psnid_store_data
struct DATA_PSNID_DOFIT  {

  // set by psnid_store_data
  char   CCID[40];
  int     NOBS ;  
  double  *MJD, *FLUXDATA, *FLUXERR, *FLUXSIM, *DUMERR0 ;
  int     *IFILT ; // sparse IFILT
  int     *IFILTOBS ; // absolute index
  double  REDSHIFT[4], REDSHIFT_ERR[4] ; // z and z_host
  double  MWEBV,    MWEBV_ERR ;
  int     SIM_NON1A_INDEX ;
} DATA_PSNID_DOFIT ;


// Oct 2013: define structure to store fit-resids.
//    Note that NOBS and obs index are the same as in DATA_PSNID_DOFIT
typedef struct RESIDS_PSNID_DOFIT_DEF {

  // epoch-dependent variables
  double *TOBS ;    // vs [obs]
  double *REJECT ;  // 1 -> epoch excluded from fit
  double *CHI2 ;
  double *FITFLUX ;
  double *FITFLUX_ERR ;

  // filter-dependent variables
  int     NFILT_USE ;
  int     *USE_FILT ;
  double  *NDOF_FILT ;  // Nepoch per filter
  double  *CHI2_FILT ;
  double  *PKMJD_FILT ;
  
} RESIDS_PSNID_DOFIT_DEF ;


// repeat with float to save memory in table var
// that stores resids for every event.
typedef struct F_RESIDS_PSNID_DOFIT_DEF {
  int    NFIT ;
  int    *IFILTOBS ; // absolute filter index
  double *MJD ;      
  float  *TOBS ;      // MJD - PKMJD, days
  float  *TREST ;     // Tobs/(1+z)
  float  *CHI2 ;      // PULL^2 
  float  *PULL ;      // DIF/ERR
  float  *FITFLUX ;      // best-fit model flux 
  float  *FITFLUX_ERR ;  // model error
  float  *DATAFLUX ;     // measured flux
  float  *DATAFLUX_ERR ; // measured error
} F_RESIDS_PSNID_DOFIT_DEF ; // to save memory in big ntuple


// define smooth fit-function (2-day bins) to overlay on plot
struct FITFUN_PSNID_DOFIT  {
  int     NOBS ;
  double *MJD, *TOBS, *FLUX, *FLUX_ERR ;
  int    *IFILT, *IFILTOBS ;
} FITFUN_PSNID_DOFIT ;


#define PSNID_MAGSCALE_for_LOOKUP 10000.0 // mag=24.234 -> index=242340
double *PSNID_FLUXCAL_LOOKUP ;
int     SIZEOF_PSNID_FLUXCAL_LOOKUP ;   // size of table, bytes

// ----------------------------------------------
// Declarations for numerical recipes (RK)
// Start with those needed for BEST method ; 
// add more declarations as needd for other methods.


float *vector(long nl, long nh) ;
int *ivector(long nl, long nh) ;
double *dvector(long nl, long nh) ;
int **imatrix(long nrl, long nrh, long ncl, long nch) ;
double **dmatrix(long nrl, long nrh, long ncl, long nch) ;
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) ;
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) ;
double ****d4tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh, long nel, long neh) ;


// include the free_xxx functions for c++ compiler
void free_dvector(double *v, long nl, long nh) ;
void free_ivector(int *v, long nl, long nh) ;
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh) ;
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) ;
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, 
		   long nch, long ndl, long ndh) ;

void free_d4tensor(double ****t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh, long nel, long neh) ;

void hunt(double xx[], int n, double x, int *jlo );
void indexx(int n, double arrin[], int indx[]);

float ran1(long *idum) ;
float ran2(long *idum) ;
float gasdev(long *idum) ;

/* moved to sntools.h (May 14 2014)
extern void get_snana_versions__(char *snana_version, char *version_photom, 
				 int len1, int len2);
*/
