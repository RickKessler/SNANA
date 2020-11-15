/********************************************************
     Created Dec 2004 by R.Kessler 
   
     For K corrections.
     Note that "(I)" indicates from user input file.


  HISTORY

  Oct 22 2015: MXFILTDEF -> 100 (was 60)

  Jan 13 2017: define MXCHAR_FILENAME
  Jan 17 2017: MXLAM_SN -> 4000 (was 3000)

  Mar 27 2017: add LEAKAGE_CUT to INPUTS.

  Apr 10 2020: MXKCOR -> 200 (was 100)

  Nov 15 2020: IVERSION_KCOR -> 4 (was 3) for reading SURVEY key

********************************************************/

bool REQUIRE_SURVEY_KCOR = false ; // flip to require SURVEY in kcor-input 

//#define VERSION_KCOR 3    // internal version
#define VERSION_KCOR 4    // internal version
#define MXFILTDEF   100    // max number of defined filters 
#define MXLAM_SN    4000  // max number of lambda bins 
#define MXLAM_PRIMARY  5000  // max number of lambda bins 
#define MXLAM_FILT  5000  // max number of lambda bins for filters
#define MXEP        150   // 110
#define MXREDSHIFT  80    // max number of redshift bins 
#define MXKCOR     200    // max number of Kcors to make (Kgg,Kgr, etc..)
#define MXAV        50    // max size of AV grid
#define MXMWEBV      4    // max number of MW E(B-V) bins
#define MXPRIMARY    6    // max number of primary standards
#define MXCHAR_FILENAME 200

#define MXSED  MXLAM_SN*MXEP

// define interpolation options: 1=linear, 2=quadratic 
#define OPT_INTERP_FILTER   1    // for filter transmission
#define OPT_INTERP_VEGAFLUX 1    // for Vega flux
#define OPT_INTERP_SNFLUX   1    // for SN flux
#define NLAMBIN_INTERP      3    // number of lambda bins to interpolate

#define FILTSYSTEM_COUNT        1
#define FILTSYSTEM_ENERGY       2

#define FNU_AB  3.631E-20    // flat Fnu for AB, erg/cm^2*s*Hz

#define KCORMIN_VALID -31.0  // set kcor to NULLVAL if outside these limits
#define KCORMAX_VALID  31.0

#define MSKFRAME_REST 1
#define MSKFRAME_OBS  2

#define MXOUTFILE     4
#define FORMAT_FITS   1
#define FORMAT_HBOOK  2

double NULLVAL = 999.0 ;

char subdir_snsed[20]      = "snsed"   ;  // default location in SNDATA_ROOT
char subdir_standards[20]  = "standards" ;
char subdir_filter[20]     = "filters" ;
char PATH_SNDATA_FILTER[200];  // $SNDATA_ROOT/filters


struct INPUTS {

  char inFile_input[200];        // name of user input file
  char inFile_snsed[200];        // (I) name of file with SN templates 
  char inFile_snsed_fudge[200];  // (I) name of file with SED fudges
  char inFile_PRIMARY[MXPRIMARY][200];      // file name of primary ref
  
  char inFile_spectrograph[200];     // name of spectrograph table (July 2016)
  char stringOpt_spectrograph[40];   // options in () after fileName

  double snsed_powerlaw ;             // override inFile_snsed (4/2019)
  char name_PRIMARY[MXPRIMARY][60];  // user-defined name of primary ref
  int  NPRIMARY ;                    // number of defined primary refs

  int  N_OUTFILE;                  // number of outfiles
  char OUTFILE[MXOUTFILE][200];    // (I) xxx.his and or xxx.fits
  char FORMAT_FLAG[MXOUTFILE];     // internally set; format flag(s)

  char FILTPATH[200];             // (I) path for filter trans files 
  char FILTPATH_ORIG[200];        // (I) idem before ENVreplace (7/2020)
  char FILTPATH_replace1[200];   
  char FILTPATH_replace2[200];   

  char MAGSYSTEM_REPLACE1[40] ;
  char MAGSYSTEM_REPLACE2[40] ;
  char MAGSYSTEM_IGNORE[40]   ;

  // filter lam shifts ... entered via command-line override only
  int    NFILTER_LAMSHIFT ; // number of non-zero LAMSHIFTs
  double FILTER_LAMSHIFT[MXFILTDEF];

  // duplicate filter set with global LAMSHIFT for systematics (May 2017)
  double LAMSHIFT_GLOBAL ;

  double LEAKAGE_CUT ;  // drop out-of-band leakage below this

  double AV_MIN ;   // (I) min AV for grid
  double AV_MAX ;   
  double AV_BINSIZE;
  int    AV_OPTION;   // (I) 1=> avg lambda;  2=> fully convolved
  int    NBIN_AV ;    // (I) number of AV bins to store KCOR grid

  int    NBIN_REDSHIFT ;    // (I) No. redshift bins in KCOR grid 
  double REDSHIFT_MIN ;     // (I) min redshift 
  double REDSHIFT_MAX ;     // (I)  
  double REDSHIFT_BINSIZE;  // (I) 

  int IRD_ZPOFF ;  // 1 => read ZPOFF file.

  double RV_MWCOLORLAW ;    // (I) A(V)/E(B-V) in LMC
  int   OPT_MWCOLORLAW ;    // color law option for Galactic extinction
  char  STR_MWCOLORLAW[60]; // definition text string

  // inputs for SNMAG text dump
  int    DUMP_SNMAG;     // text-dump option of SN mag vs. epoch and filter
  char   SN_TYPE[20];  // (I) optional descriptor (Ia, IIP, etc ...)
  double SN_MAGOFF ;   // (I) global mag-offset 

  int FASTDEBUG ;    // 1 => skip calculations to run thru quickly
                     //      Mainly to debug hbook problems.

  int SKIPKCOR;      // (I) 1=> skip K-corrections (just do mags & zeropoints)
  int FLUXERR_FLAG;  // (I) flux errors are present => read & ignore
  int PLOTFLAG_SN ;   // (I) 0=> no SN SED plot to save space
  int PLOTFLAG_PRIMARY ;

  double TREF_EXPLODE ;  // defines rest-frame explosion date (days)

  int NLAMBIN_FT; // Number of Fourier Transform bins (must be power of 2)

} INPUTS ;


// for most cases, NAME = NAME_INPUT. However, there is an option
// where NAME_INPUT->NAME (e.g, VEGA->AB) means transform the
// input MAGSYS with 'NAME_INPUT' to MAGSYS with 'NAME.
typedef struct {
  char   NAME_INPUT[60], NAME[60];
  int    INDX_INPUT, INDX ;
  double OFFSET_INPUT, OFFSET ;
  int    DO_TRANSFORM ;  // T -> NAME_INPUT != NAME
  char   SURVEY_NAMES[200]; // optional assoication with survey names
} MAGSYSTEM_DEF ;

typedef struct {
  int INDX ;
  char NAME[60];
} FILTSYSTEM_DEF ;

typedef struct {
  // broad band filters
  char   FILTNAME[60];
  char   FILENAME[200];
  double MAGREF ;         // defines filter system

  // check if it's a synthetic filter from the spectrograph
  int    IFLAG_SYN ;       // 0 or 1
  double LAMRANGE_SYN[2] ; // min:max wavelength range of SYN_FILTER
} INPUT_FILTER_DEF ;

typedef struct {
  double LAMRANGE[2];    // wave range to apply OOB
  double TRANS_MAX_RATIO;
} OOB_DEF; // out-of-band transmission


int   IEPOCH_SNPEAK[MXFILTDEF] ;   // epoch index at SN peak (t=0)
int   IEPOCH_SN15DAY[MXFILTDEF]  ;

// option to use UBVRI Vega offsets from Astier06
int    VEGAOFF_FLAG_ASTIER06 ; 
double VEGAOFF_ASTIER06[7] = { -99.0, 0.02, 0.03, 0.03, 0.03, 0.024,.03 }; 

// define array of MW E(B-V) values to check linear relation
double MWEBV_LIST[MXMWEBV+1] = { 0.0, 0.1, 0.20, 0.40, 0.80 };


// define SN template
struct SNSED {

// start with quantities read from SN template file
  int    NEPOCH;                    // number of template epochs
  double EPOCH[MXEP];               //  rest-frame day since max 
  double LAMBDA[MXEP][MXLAM_SN];    //  central lambda value 
  double FLUX_WAVE[MXEP][MXLAM_SN];    // dE/dlam at lam 
  double FLUX_NU[MXEP][MXLAM_SN];      // dE/dnu at lam 


  // inferred from SN template file
  double LAMBDA_MIN;      // (I) from user input file 
  double LAMBDA_MAX;      // (I)  "   
  double LAMBDA_BINSIZE;  // determined from SN templates 
  int    NBIN_LAMBDA;      

  double TREST_MIN ;  // (I) min Trest read from SED
  double TREST_MAX ;  // (I) max ..
  double TREST_BIN ;  // binsize of SED (days)

  // computed from above
  double FLUXSUM[MXEP][MXFILTDEF];    // flux sum 

  // use float for MAGOBS  to save memory

  double MAG_RST[MXEP][MXFILTDEF];     // rest-frame mag (with AV=0)
  float *****R4MAG_OBS ;

  // define change in obs mag w.r.t. MW E(B-V)
  double ****MW_dXT_dEBV ;

} SNSED ;



// define SNSED FUDGE
struct SNSED_FUDGE {

  double  LAMBDA[MXLAM_SN];    // central lambda value 
  double  SCALE[MXLAM_SN];     // flux-fudge at each lambda

  double LAMBDA_MIN;      // (I) from user input file 
  double LAMBDA_MAX;      // (I)  "   
  double LAMBDA_BINSIZE;  // determined from SN templates 
  int   NBIN_LAMBDA;      

} SNSED_FUDGE ;



// define "primary" SED [VEGA, BD17, etc ...]

struct PRIMARYSED {   
  int USE ;    // 1-> this primary is actually used.
  int IS_AB;   // 1-> AB system

  double ZP[MXFILTDEF];          // primary zero point
  double FLUXSUM[MXFILTDEF];   // Flam dlam or Fnu dnu/nu
  double SYNMAG[MXFILTDEF];  // synthetic mags

  // define raw values read from file (note that labmda bins are not uniform)

  int     NBIN_LAMBDA_RAW ;              // raw number of lambda bins
  double  LAMBDA_MIN_RAW ;
  double  LAMBDA_MAX_RAW ;
  double  LAMBDA_RAW[MXLAM_PRIMARY];       // raw lambda bins
  double  FLUX_WAVE_RAW[MXLAM_PRIMARY];    // raw energy flux (erg/s/cm^2/A)

  // define variables re-binned for uniform labmda bins
  int     NBIN_LAMBDA ;              // raw number of lambda bins
  double  LAMBDA_MIN ;
  double  LAMBDA_MAX ;
  double  LAMBDA_BINSIZE;
  double  LAMBDA[MXLAM_PRIMARY];       // raw lambda bins
  double  FLUX_WAVE[MXLAM_PRIMARY];    // raw energy flux (erg/s/cm^2/A)
  double  FLUX_NU[MXLAM_PRIMARY];

  // add global info 4/25/2009
  double MAGSYSTEM_OFFSET   ;
  char   MAGSYSTEM_NAME[40] ;
  char   MAGSYSTEM_SEDFILE[200];


} PRIMARYSED[MXPRIMARY];




// define FILTER structure
int NFILTDEF;   // number of defined filtes
int NFILTPATH;  // number of filters paths (SDSS, SNLS ...)
int NFILTDEF_SYN; // number of synthetic filters from SPECTROGRAPH

struct  FILTER 
{
  char name[60] ;  // (I) full name of each filter 
  char band[4];    // (I) single-char representation of each filter/band
  char file[200] ;  // (I) data text file of each filter response 
      
  int IFLAG_SYN;    // (I) >0 if synthetic filter from spectrograph

  // define input filt transmssion 
  double TRANS[MXLAM_FILT]  ;    // efficiency per lambda bin 
  double LAMBDA[MXLAM_FILT] ;    // wavelen, A 
  double LAMBDA_MIN ;       // min lambda with non-zero trans
  double LAMBDA_MAX ;       // max lam with non-zero trans 

  // define lambda range with leakage cut
  double LAMBDA_MIN_LEAKAGE_CUT;
  double LAMBDA_MAX_LEAKAGE_CUT;
  double LEAKAGE_FRAC;  // out-of-band leakage for flat dN/dlam spectrum

  double LAMAVG ;           // <lam> in each filter
  int    NBIN_LAMBDA ;      // No. lambda bins 
  int  IFLAG_DUPLICATE ; // T -> duplicate band with LAMSHIFT (May 2017)

  // SSUM = int S(lam,nu) dlam,dnu/nu 
  double SSUM_PRIM ; // in primary ref lambda bins
  double SSUM_SN ;   // in SN lambda bins
  double SUMTOT;     // normal integral of T(lam) dlam

  char   MAGSYSTEM_NAME[40];     // (I) mag system (Vega, ABnu, ABlam ...)
  int    MAGSYSTEM_INDX;        // index corresponding to NAME
  int    MAGSYSTEM_INDX_INPUT; 
  double MAGSYSTEM_OFFSET;      // global system mag offset(same for each filt)
  double MAGFILTER_REF;         // (I) mag of primary ref 
  double MAGFILTER_ZPOFF;       // (I) AB offset from optional ZPOFF.DAT file
  double MAGFILTER_ZP;          // filter-dependent zero point

  int    IPATH ;                   // integer index to PATH
  char   PATH[200] ;              // subdir in $SNDATA_ROOT/filters/
  char   PATH_ORIG[200];          // idem, before ENVreplace
  char   FILTSYSTEM_NAME[20] ;    // (I) energy or count
  int    FILTSYSTEM_INDX ;        // energy or count
  char   SURVEY_NAMES[200];       // comma-sep list of surveys (Nov 2020)

  double SNMAG0_CHECK_VALUE ;     // external reference for check 
  char   SNMAG0_CHECK_NAME[40] ;  // comment on origin of refernece

  int MASKFRAME;  // bits 1,2 => rest,obs

  OOB_DEF OOB; // out-of-band info

} FILTER[MXFILTDEF] ;


double FILTER_LAMBDA_MAX; // max lambda among all filters
double FILTER_LAMBDA_MIN; // min lambda ...


struct R4KCOR_GRID {
  float ****VALUE ;
  float ****REDSHIFT ;         // z in same bins 
  float ****EPOCH ;            // rest frame day in same bins 
} R4KCOR_GRID ;


// K cor list applies to the grid and to the lightcurve list
int  NKCOR;                    // (I) No. K correction matrices to make 
char KCORLIST[MXKCOR][2][40];  // list of K cor filters (40 char/filter) 
char KCORSYM[MXKCOR][40];      // list of user-defined K cor symbols

int NKCOR_EXTRA;  // used to get synthetic mags for rest-frame filters


int  NZPOFF; // size of ZPOFF list
struct ZPOFF {
  int    IFILTPATH[MXFILTDEF];     // filter path index
  char   FILTLIST[MXFILTDEF] ;     // list of single-char filters
  double FILTVALUES[MXFILTDEF] ;   // ZPOFF values
} ZPOFF_STORE ;

// large temp arrays must be global to avoid valgrind errors
double  DUMMYERR[MXSED] ;  
double  DUMMYFLUX[MXSED] ;


// variables for fits to avoid redefining them in each function.
struct STRFITS {
  char  F4[4], U4[4], C20[4], D8[4], blank[8] ;
  char *tName[MXKCOR] ;
  char *tForm[MXKCOR] ;
  char *tUnit[MXKCOR] ;
} STRFITS ;

// =====================================================
//
//     !!!!  DECLARE FUNCTIONS   !!!!!
//
// =====================================================

int   rd_input(void) ;
void  get_NZBIN(void);
void  get_NAVBIN(void);
void  kcor_input_override(int OPT);
void  primary_override(char *primName, char *fileName);
int   kcor_ini(void) ;
void  set_kcorFile_format(void) ;
int   malloc_ini(void);
int   kcor_out(void) ;
int   kcor_grid(void) ;
void  primarymag_zp(int iprim);  // integrated fluxes, mags, and zero points/
void  primarymag_zp2(int iprim);
void  primarymag_summary(int iprim); 
int   snmag(void);    // determine integrated fluxes and magnitudes

void  storeFilterInfo(INPUT_FILTER_DEF *INPUT_FILTER,
		      MAGSYSTEM_DEF *MAGSYSTEM, FILTSYSTEM_DEF *FILTSYSTEM) ;

int   rd_filter( int ifilt ) ;  // read filter transmissions from files 
void  parse_MAGREF(char *FILTNAME, char *TXT_MAGREF, double *MAGREF ) ;
void  parse_MAGSYSTEM(char *MAGSYSTEM_TMP, MAGSYSTEM_DEF *MAGSYSTEM) ;
void  parse_OOB(char *bandList, double *LAMRANGE, double RATIO);
void  addOOBTrans_filter(int ifilt);

void  rd_ZPOFF(char *sdir);     // read optional ZPOFF.DAT file
int   rd_snsed(void);           // read SN template spectra from file 
void  hardWire_snsed_bins(void) ; // in case SN SED is not defined.
void check_snsed_bins(void) ;

void  ADDFILTERS_LAMSHIFT_GLOBAL(void);

double GET_SNSED_FUDGE(double lambda);
double get_ZPOFF(char *filtname, int IPATH);  // return ZPOFF for *filtname
int   FT_snsed(void);  // comput Fourier Transforms
void  parse_FILTER_LAMSHIFT(int *indx_ARGV) ;

int   rd_primary(int indx, char *subdir ); 
int   index_primary(char *name);

void filter_filename(FILE *fp, char *filter_path, char *filename );

void magflux_info( int iprim, int ifilt, int ilam, int epoch
		  ,double *lam, double *flux, double *w1, double *w2 ) ;


// return filter transmission at this "lam" 
double filter_trans8 ( double lam, int ifilt, int idump );  

// return SN flux at "epoch" and for "lambda/(1+redshift)"
double snflux8 ( double epoch, double lambda, double redshift, double av );

// return flux of primary standard
double primaryflux8 ( int iprim, double lambda );

// Compute K correction as in Nugent 2002 
// Note opt = OPT_KCOR_EPHOT, OPT_KCOR_NPHOT 
void kcor_eval ( int iopt
		 ,double av, double redshift, double epoch
		 ,int ifilt_rest, int ifilt_obs
		 ,int FLAG_MAGOBS
		 ,double *kcor_value, double *kcor_error
		 ,double *overlap, double *flux_obs
                        ) ;

// convert  epoch (days) to integer index
int  index_epoch ( double epoch );

// returns rest and observer filters for  user "ikcor" index
void index_filter ( int ikcor, int *ifilt_rest, int *ifilt_obs );

// return ifilt for FILTER[ifilt] struct
int INTFILTER_kcor(char *filterName) ;

// convert dE/dlam (erg/s/cm^2/A) to dE/dnu and dN/dnu
void flux_converter ( double lambda, double flux_wave,
			double *flux_nu, double *flux_count ); 

// rebin *flux_in array to *flux_out with same lambda bins as SN templates
void rebin_primary ( int nblam_in, double  *lam_in, double *flux_in, 
		     int *nblam_out, double *lam_out, double *flux_out );


void wr_hbook(char *outFile);  // store KCOR with hbook 

void wr_fits(char *outFile);   // store in FITS file
void wr_fits_HEAD(fitsfile *fp);
void wr_fits_PRIMARY(fitsfile *fp);
void wr_fits_SNSED(fitsfile *fp);
void wr_fits_ZPT(fitsfile *fp);
void wr_fits_FilterTrans(fitsfile *fp);
void wr_fits_KCOR(fitsfile *fp);
void wr_fits_MAGS(fitsfile *fp);
void wr_fits_SPECTROGRAPH(fitsfile *fp);

void wr_fits_errorCheck(char *comment, int status) ;

void wrsnmag_text(void);

void dmp_filtersummary(void);

// ======== END OF FILE ===========
