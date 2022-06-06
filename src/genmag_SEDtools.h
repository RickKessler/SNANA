/*******************************************
 Created Oct 14, 2009 by R.Kessler

 Global variables used for SED-based models
 such as SALT2 or SN-explosion models. 

 Jun 24, 2010: MXLAMPOW -> MXLAMPOW_SEDMODEL = 2 (instead of 0) 
               and define new variable NLAMPOW_SEDMODEL

               Add 'ilampower' argument to interp_flux_SEDMODEL

 Jul 02, 2010: add SIMSED.PARVAL_STRING[ised]
 Jul 16, 2010: add FLUX_ERRFLAG
 Jul 19, 2010: increase MXFILT_SEDMODEL from 10 to 12
 Aug 31, 2010: increase MXFILT_SEDMODEL from 10 to MXFILTINDX
 Sep 02, 2010: declare zero_flux_SEDMODEL

 Apr 04, 2011: increase MXZBIN_SEDMODEL -> 500 (was 300)

 May 15, 2011: increase MXBIN_LAMSED_SEDMODEL from 2100 -> 5000

 Jun 09, 2011: define SEDMODEL.LAMMIN_ALL and SEDMODEL.MAXLAM_ALL

 Jun 20, 2011  increase MXPAR_SEDMODEL  from 10 -> 100

 Nov 1, 2011: MXBIN_DAYSED_SEDMODEL -> 200 (was 130)

 Nov 11, 2011: MXBIN_LAMFILT_SEDMODEL -> 800 (was 600)
 Feb 20, 2013: MXBIN_LAMFILT_SEDMODEL -> 900 (was 800)
 Mar 16, 2014: MXSEDMODEL -> 8000 (was 5000)

 Apr 26, 2014: MXFILT_SEDMODEL[20] -> MXFILT_SEDMODEL[MXFILT_SEDMODEL]
                 (fixes crash with jpas)

 Aug 05, 2015: MXBIN_LAMFILT_SEDMODEL -> 1500 (was 900) for WFIRST W band

 Aug 11, 2015: define void  shiftPeakDay_SEDMODEL();

 July 24 2016: define stuff for SPECTROGRAPH

 Jul 26, 2016: MXBIN_LAMFILT_SEDMODEL -> 2000 (was 1500)

 Jul 30 2016: define SPECTROGRAPH_SEDMODEL struct and a few 
              SPECTROGRAPH functions.

 Mar 6 2017: declare ISED_SEDMODEL

 Apr 17 2018: MXBIN_DAYSED_SEDMODEL -> 400 (was 250)

 May 21 2018: define struct INPUTS_SEDMODEL

 Jan 3 2019: 
   + new function fetch_parVal_SEDMODEL
   + new fortran-mangles for fetch_parinfo_sedmodel__ & fetch_parval_sedmodel__

 Aug 23 2019: MXBIN_LAMFILT_SEDMODEL -> 2400 (was 2000)

********************************************/

// define bounds for filter and SED arrays
#define MXSEDMODEL          8000   // max number of SED surfaces
#define MXFILT_SEDMODEL     MXFILTINDX     // max internal filter index
#define MXBIN_LAMFILT_SEDMODEL 2400   // length of largest filter file
#define MXBIN_LAMSED_SEDMODEL  5000   // max # lambda bins for SED
#define MXBIN_DAYSED_SEDMODEL   400 // max # epoch (day) bins for SED
#define MXBIN_SED_SEDMODEL     MXBIN_LAMSED_SEDMODEL*MXBIN_DAYSED_SEDMODEL 
#define MXBIN_PRIMARY_SEDMODEL  5000
#define MXZBIN_SEDMODEL         500     // max redshift bins
#define MXPAR_SEDMODEL         100      // max NPAR describing SEDs
#define MXLAMPOW_SEDMODEL 4   // max power of lambda for pre-computed integrals


#define LAMMIN_SEDMODEL 100     // global min Lambda for reading SED model
#define LAMMAX_SEDMODEL 30000   // global max Lambda


// define redshift bins for SEDMODEL_TABLE
#define NZBIN_SEDMODEL_DEFAULT   270
#define ZMIN_SEDMODEL_DEFAULT    0.005 
#define ZMAX_SEDMODEL_DEFAULT    2.500

#define OPTMASK_FLUXERR_SEDMODEL 1  // mask to read flux errros
#define OPTMASK_DAYLIST_SEDMODEL 2  // mask to allow non-uniform DAY array
#define OPTMASK_T0SHIFT_PEAKMAG  4  // shift T=0 to peak

typedef struct REDSHIFT_SEDMODEL_TYPE {
  double  ZMIN,    ZMAX;  
  double  LOGZMIN, LOGZMAX, LOGZBIN ;  // log10
  int     NZBIN;
  double  ZTABLE[MXZBIN_SEDMODEL];
  double  LOGZTABLE[MXZBIN_SEDMODEL];
} REDSHIFT_SEDMODEL_TYPE ;

struct REDSHIFT_SEDMODEL_TYPE   REDSHIFT_SEDMODEL ;


int NLAMPOW_SEDMODEL ;      // highest lambda power of pre-computed integrals
int IFILTMAP_SEDMODEL[100]; // translates absolute index to sparse index

int  NFILT_SEDMODEL ;
char FILTLIST_SEDMODEL[MXFILT_SEDMODEL] ;

int ISED_SEDMODEL;  // "ISED" at corner; else -9 for interpolation

struct FILTER_SEDMODEL {
  int     NLAM ;             // number of lambda bins for filter
  int     ifilt_obs ;        // absolute filter index passed by user
  double  magprimary;        // mag of primary reference (vega, AB ...)
  double  lam[MXBIN_LAMFILT_SEDMODEL];   // wavelengths
  double  transSN[MXBIN_LAMFILT_SEDMODEL]; // transmission for SN 
  double  transREF[MXBIN_LAMFILT_SEDMODEL]; // transmission for ref
  double  transSN_MAX ;      // Max transSN (added 4/30/2011)
  double  transREF_MAX ;     // Max transREF
  double  lamshift ;         // user-shift (A)
  double  lamstep;           // step, A
  double  lammin, lammax ;   // min/max wavelength
  double  mean;              // mean wavelength of filter
  double  ZP;                // reference zeropoint
  char    name[40];          // full name of filter
  char    survey[40];        // name of survey (Nov 2020)
} FILTER_SEDMODEL[MXFILT_SEDMODEL] ;


struct MWXT_SEDMODEL {
  int    OPT_COLORLAW ; // 89, 94, 99 ...
  double RV ;
} MWXT_SEDMODEL ; // Galactic extinction


struct PRIMARY_SEDMODEL {

  // info filled before filters are defined
  int     NLAM;
  double  lam[MXBIN_PRIMARY_SEDMODEL];   // wavelengths
  double  flux[MXBIN_PRIMARY_SEDMODEL];  // flux
  double  lamstep;                    // step, A
  double  lammin, lammax ;   // min/max wavelength
  char    name[20];

  // information filled after filters are defined
  double  MAG[MXFILTINDX] ;
  double  MAG_SORTED[MXFILTINDX] ;
  double  LAM_SORTED[MXFILTINDX] ;
} PRIMARY_SEDMODEL ;


// all info except for the SED flux
struct SEDMODEL {
  int    NSURFACE ;   // 2=default for SALT2; arbitrary for SIMSED model
  double FLUXSCALE ;  // scale all fluxes by this amount (default=1)
  double MAGERR_FIX ;  // fix mag-error if > 0
  double LOGZBIN ;     // Flux-table binsize for log10(redshift); default=.02

  int    NPAR ;        // Npar describing SEDs
  char   PARNAMES[MXPAR_SEDMODEL][40];
  double PARVAL[MXSEDMODEL][MXPAR_SEDMODEL];
  char   PARVAL_STRING[MXSEDMODEL][MXPAR_SEDMODEL][20]; // parval string to preserve format
  int    NBIN_PARVAL[MXPAR_SEDMODEL] ; // Number of bins per par val
  double PARVAL_MIN[MXPAR_SEDMODEL] ;
  double PARVAL_MAX[MXPAR_SEDMODEL] ;
  double PARVAL_BIN[MXPAR_SEDMODEL] ; // bin size
  int    IPAR_NON1A_INDEX ;   // ipar to use to fill SIM_xNON1A_INDEX

  char   FILENAME[MXSEDMODEL][80]; // NSURFACE of them

  double RESTLAMMIN_FILTERCEN;  // min-lambda for <lamfilt>/(1+z) 
  double RESTLAMMAX_FILTERCEN;  // max-lambda for <lamfilt>/(1+z)

  double LAMMIN_ALL ;  // min LAM covered by all SEDs
  double LAMMAX_ALL ;  // max LAM covered by all SEDs

  // define DAY and LAM binning for each SED
  int    NDAY[MXSEDMODEL],   NLAM[MXSEDMODEL] ;
  double LAMMIN[MXSEDMODEL], LAMMAX[MXSEDMODEL], LAMSTEP[MXSEDMODEL];
  double DAYMIN[MXSEDMODEL], DAYMAX[MXSEDMODEL], DAYSTEP[MXSEDMODEL];
  double *DAY[MXSEDMODEL];  // Aug 2017 - allows non-uniform DAY bins

  double DAYMIN_ALL, DAYMAX_ALL; // min & max day among all SEDs

  int MXDAY ; // max number of epochs to store => used to malloc (Nov 2011)

  int OPTMASK;  // 2 ==> use DAY array instead of uniform DAYSTEP

} SEDMODEL ;


// Aug 2013: move MWXT table from genmag_SALT2.h to here
double **SEDMODEL_TABLE_MWXT_FRAC ;
double SEDMODEL_MWEBV_LAST ;

// July 2016: define array for host extinction
double **SEDMODEL_TABLE_HOSTXT_FRAC ;
struct{double AV, z, RV; } SEDMODEL_HOSTXT_LAST ;


// define TEMP structure that gets over-written for each SED.
// This is mainly to avoid wasting memory storing each SED
// because we only need to store the flux-integrals 
typedef struct SEDMODEL_FLUX_DEF {
  int    NDAY, NLAM ;
  double LAMMIN, LAMMAX, LAMSTEP;  // based on SED
  double DAYMIN, DAYMAX, DAYSTEP;

  double *LAM ;       // lambda array
  double *DAY ;       // epoch array
  double *FLUX ;      // template flux array
  double *FLUXERR ;   

  double TSHIFT; // shift applied so that T=0 at peak (July 10 2018)

  // define fine-binned arrays so that linear interp
  // can be used in the integrations. The fine-bins
  // are only in the interpolated LAMDA dimension
  // (epoch is not interpolated in init_flux_SEDMODEL)
  int     N_FINEBIN ;
  double *FINEBIN_FLUX ;
  double *FINEBIN_LAM ;
  double  FINEBIN_LAMSTEP;

} SEDMODEL_FLUX_DEF ; // define structure with SED fluxes


SEDMODEL_FLUX_DEF  TEMP_SEDMODEL ;
SEDMODEL_FLUX_DEF *SEDMODEL_STORE ; // used with SPECTROGRAPH option

int NVAR_FIRSTBIN_SEDMODEL ;
struct FIRSTBIN_SEDMODEL  {
  int    NSED;
  char   VARNAME[20];  // variable name (DAY or LAM)
  int    NBIN;
  double VAL0; // first element of array
  double BINSIZE;
} FIRSTBIN_SEDMODEL[10] ; // store binning info about 1st SED


//#define   IVERSION_SEDBINARY 2          // Apr 12 2018
#define   IVERSION_SEDBINARY 3          // Apr 30 2018
#define   PADWORD_SEDBINARY  77777.     // pad-word after header (4.12.2018)
#define   ZEROWORD_SEDBINARY 909090909. // indicates list of zeros to follow
#define   MINZEROLIST_SEDBINARY 10      // at least this many to compress
int   NSEDBINARY;  // lenth of SEDBINARY array
float SEDBINARY[MXBIN_SED_SEDMODEL]; 

// Pre-calculate integrals: indices are
//  - ifilt
//  - log10(z)
//  - additional power of lambda inside integral
//  - epoch
//  - template-SED id



// Nov 24,2009: define array of pointers for dynamic allocation.
// Note that SEDMODEL.NSURFACE is the total number of surfaces

#define NDIM_SEDMODEL_FLUXTABLE  5    // number of dimensions to store 
#define IDIM_SEDMODEL_FILTER   1
#define IDIM_SEDMODEL_REDSHIFT 2
#define IDIM_SEDMODEL_LAMPOW   3
#define IDIM_SEDMODEL_DAY      4
#define IDIM_SEDMODEL_SED      5

float    *PTR_SEDMODEL_FLUXTABLE ;  // pointer array
long int  ISIZE_SEDMODEL_FLUXTABLE;  // total size
long int  NBTOT_SEDMODEL_FLUXTABLE;  // total number of fluxtable bins
int       NBIN_SEDMODEL_FLUXTABLE[NDIM_SEDMODEL_FLUXTABLE+1];
long int  N1DBINOFF_SEDMODEL_FLUXTABLE[NDIM_SEDMODEL_FLUXTABLE+1];
char      VARNAME_SEDMODEL_FLUXTABLE[NDIM_SEDMODEL_FLUXTABLE+1][12] ;


// SPECTROGRAPH STUFF (July 2016)
#define JFILT_SPECTROGRAPH 0  // to use some existing arrays
struct {
  int    NBLAM_TOT ;
  double *LAMMIN_LIST, *LAMMAX_LIST; // min,max lambda per bin
  double *LAMAVG_LIST;           // avg lambda per bin (e..g, for Extinction)
  double *ZP_LIST ;              // ZP per bin
  double LAMMIN, LAMMAX ;        // global min,max

} SPECTROGRAPH_SEDMODEL ;


// few inputs read from sim-input file
struct {
  int    OPTMASK_T0SHIFT_EXPLODE; // option to define T=0 at explosion
  double UVLAM_EXTRAPFLUX;       // extrapolate SED down to UV region
  double MINSLOPE_EXTRAPMAG_LATE;   // min mag/day slope for extrapolation
} INPUTS_SEDMODEL;

// ==============================================
// function declarations

int reset_SEDMODEL(void);

int init_primary_SEDMODEL(char *refname, int NLAM, 
			  double *LAMLIST, double *FLUXLIST ) ;

int init_filter_SEDMODEL(int ifilt_obs, char *filter_name, char *survey_name,
			 double magprimary, int NLAM,  double *LAMLIST, 
			 double *TRANSSNLIST, double *TRANSREFLIST, 
			 double LAMSHIFT ) ;

void filtdump_SEDMODEL(void);   // one-line dump per filter.

void init_redshift_SEDMODEL(int NZbin, double Zmin,  double Zmax);

void init_MWXT_SEDMODEL(int OPT_COLORLAW, double RV) ;

double interp_primaryFlux_SEDMODEL(double lam) ;
double interp_primaryMag_SEDMODEL(double lam) ;  // for SPECTROGRAPH bins

void malloc_FLUXTABLE_SEDMODEL(int NFILT, int NZBIN, int NLAMPOW, 
			       int NDAY, int NSED);
void malloc_SEDFLUX_SEDMODEL (SEDMODEL_FLUX_DEF *SEDMODEL, 
			     int NBIN_DAY, int NBIN_LAM, int NBIN_SED);

int IFILTSTAT_SEDMODEL(int ifilt_obs, double z) ;

void init_FINEBIN_SEDMODEL(int ised);
void init_flux_SEDMODEL( int ifilt_obs, int ised );
void zero_flux_SEDMODEL(void);

double interp_flux_SEDMODEL(int ISED, int ilampower, int ifilt_obs, 
			    double z, double Trest );
double get_flux_SEDMODEL(int ISED, int ilampow, int ifilt_obs,
			 double z, double Trest) ;
double get_magerr_SEDMODEL(int ISED, int ifilt_obs,
			   double z, double Trest) ;

double getFluxLam_SEDMODEL(int ISED, int IEP, double TOBS, double LAMOBS,
                           double z, char *funCall );

void get_LAMTRANS_SEDMODEL(int ifilt, int ilam, double *LAM, double *TRANS);

void get_LAMRANGE_SEDMODEL(int opt, double *lammin, double *lammax);
void checkLamRange_SEDMODEL(int ifilt, double z, char *callFun) ;
void get_DAYBIN_SEDMODEL(int ISED, double DAY, int *IDAY, double *FRAC);

void pack_SEDBINARY(int OPT); // OP+1 => PACK;  OPT=-1 => UNPACK

long int INDEX_SEDMODEL_FLUXTABLE(int ifilt, int iz, 
			     int ilampow, int iep, int  ised);

int get_SEDMODEL_INDICES( int IPAR, double LUMIPAR, 
			 int *ILOSED, int *IHISED);

void check_sedflux_bins(int ised, char *VARNAME, 
			int NBIN, double VAL0, double BINSIZE);
void check_surveyDefined_SEDMODEL(void);

double gridval_SIMSED(int ipar, int ibin);
double nearest_gridval_SIMSED (int ipar, double lumipar );

int NPAR_SEDMODEL(void);
int NSED_SEDMODEL(void);

int  fetch_parInfo_SEDMODEL(int ipar, char *parname,int *NBIN,double *range); 
void fetch_parVal_SEDMODEL(int ISED, int IPAR, double *PARVAL);

int IPAR_SEDMODEL ( char *parName); // return parName index

void fill_TABLE_MWXT_SEDMODEL(double RV, double mwebv) ;
void fill_TABLE_HOSTXT_SEDMODEL(double RV, double AV, double z) ;

double filterTrans_BessB(double lam) ;


void T0shiftPeak_SEDMODEL(SEDMODEL_FLUX_DEF *SEDFLUX, int vboseFlag);
void T0shiftExplode_SEDMODEL(int OPTMASK, SEDMODEL_FLUX_DEF *SEDFLUX, 
			     int vboseFlag);

void UVLAM_EXTRAPFLUX_SEDMODEL(double UVLAM, SEDMODEL_FLUX_DEF *SEDFLUX);
void set_UVLAM_EXTRAPFLUX_SEDMODEL(float UVLAM_MIN);

void print_ranges_SEDMODEL(SEDMODEL_FLUX_DEF *SEDFLUX);

void FLUX_SCALE_SEDMODEL(double SCALE, SEDMODEL_FLUX_DEF *SEDFLUX);

bool found_fluxerr_SEDMODEL(char *sedFile);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// SPECTROGRAPH functions

void INIT_SPECTROGRAPH_SEDMODEL(char *MODEL_NAME, int NBLAM, 
				double *LAMMIN_LIST, double *LAMMAX_LIST) ;

double getZP_SPECTROGRAPH_SEDMODEL(double LAMMIN, double LAMMAX, int DUMPFLAG);

void getSpec_SEDMODEL(int ised, 
		      double MWEBV, double RV_host, double AV_host,
		      double z, double MU, double Tobs, double MAGOFF,
		      double *FLUXGEN_LIST, double *GENMAG_LIST ) ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// define mangled functions with underscore (for fortran)

int reset_SEDMODEL__(void);

int init_primary_sedmodel__(char *refname, int *NLAM, 
			    double *LAMLIST, double *FLUXLIST );

int init_filter_sedmodel__(int *ifilt_obs, char *filter_name, 
			   char *survey_name, 
			   double *magprimary,
			   int *NLAM,  double *LAMLIST, 
			   double *TRANSSNLIST, 
			   double *TRANSREFLIST, 
			   double *LAMSHIFT) ;

void init_redshift_sedmodel__(int *NZbin, double *Zmin,  double *Zmax) ;

void init_mwxt_sedmodel__(int *OPT_COLORLAW, double *RV) ;

void get_lamrange_sedmodel__(int *opt, double *lammin, double *lammax);

int fetch_parinfo_sedmodel__(int *ipar,char *parname,int *NBIN,double *range); 
void fetch_parval_sedmodel__(int *ISED, int *IPAR, double *PARVAL);

void set_uvlam_extrapflux_sedmodel__(float *UVLAM);


// end of file.
