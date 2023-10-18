
// Created July 2016 by R.Kessler

#define MXSPEC                 MXSPECTRA  // max Nspec per event   
#define MXTEXPOSE_SPECTROGRAPH 50   // max size of TEXPOSE grid
#define MXLAM_SPECTROGRAPH     10000 // ->10k on May 27 2020 (was 2400)
#define MXLAM_SPECTROGRAPH_EXTEND 20 // max number of extra bins for resolution
#define MXLAMSMEAR_SPECTROGRAPH 40  // max number of smeared lambda bins
#define FITSTABLE_NAME_SPECTROGRAPH  "SPECTROGRAPH" 
#define NCOL_noSNR 3    // Ncolumns before SNR values
#define TEXPOSE_INFINITE_SPECTROGRAPH 1.0E8 // flag to ignore noise
#define ILIST_RANDOM_SPECTROGRAPH 3   // separate list for ran Trest
#define ISTREAM_RANDOM_SPECTROGRAPH 1 // independent random stream.
#define ISPEC_PEAK        3*MXSPECTRA // imjd=ISPEC_PEAK -> fetch peak spec

#define WRITE_MASK_SPEC_DEFAULT  2  // lammin lammax Flam FlamERR SIM_FLAM
#define WRITE_MASK_SPEC_SED_TRUE 4  // <lam> and SIM_FLAM only 

int  SPECTROGRAPH_USEFLAG ;
int  NERR_SNR_SPECTROGRAPH ;
int  NERR_BADSNR_SPECTROGRAPH ;

struct {

  char   INFILE_NAME[MXPATHLEN]; // input table with SPECBIN keys
  char   INSTRUMENT_NAME[60];
  double MAGREF_LIST[2];  // two MAGREFs to compute ZP and SIGSKY
  int    NBIN_LAM, NBIN_LAM_noREBIN;
  int    NBIN_TEXPOSE ;
  int    NREBIN_LAM;      // rebin by this integer value (Sep 2 2016)
  double SNR_POISSON_RATIO_ABORT_vsMAGREF ;
  double SNR_POISSON_RATIO_ABORT_vsTEXPOSE ;
  double MAGSNR_TOLERANCE_ABORT ; // tolerance on (m1-m0) - 2.5log(SNR0/SNR1)

  // define spectro table read from input file: 
  // needs to be malloced with NBIN_xxx
  double LAM_MIN,     LAM_MAX; // global min & max lambda
  double TEXPOSE_MIN, TEXPOSE_MAX; // global min & max TEXPOSE

  double *LAMMIN_LIST, *LAMMAX_LIST, *LAMAVG_LIST ;     // per spectro bin
  double *LAMBIN_LIST, *LAMSIGMA_LIST ; // binSize & resolution vs. lambda
  double TEXPOSE_LIST[MXTEXPOSE_SPECTROGRAPH] ;    // per expTime bin
  double LOGTEXPOSE_LIST[MXTEXPOSE_SPECTROGRAPH];  // used for ZP-interp
  double **SNR0, **SNR1  ;  // per spectro bin & expTime
  bool   *ISLAM_EXTEND_LIST; // True -> extended lam bin for lam-resulution

  int WRITE_MASK ; 

  // computed from table inputs
  double **ZP, **SQSIGSKY ; // per spectro & expTime bin

  // syn-filter info from kcor-iput file; filled when reading kcor file
  int    NSYN_FILTER;
  int    SYN_IFILTDEF_LIST[MXFILTINDX];
  int    SYN_IFILTINV_LIST[MXFILTINDX]; // map ifiltdef -> sparse ifilt
  bool   IS_SYN_FILTER[MXFILTINDX];     // vs. ifiltdef
  char   SYN_FILTERLIST_BAND[MXFILTINDX] ;
  char   SYN_FILTERLIST_NAME[MXFILTINDX][20] ;
  double SYN_FILTERLIST_LAMMIN[MXFILTINDX] ;
  double SYN_FILTERLIST_LAMMAX[MXFILTINDX] ;

} INPUTS_SPECTRO ;

#define MXVALUES_SPECBIN  10+2*MXTEXPOSE_SPECTROGRAPH
double  VALUES_SPECBIN[MXVALUES_SPECBIN];


// ------ GENERATED SPECTRA ------                                                        
struct {
  int    NMJD_TOT, NMJD_PROC ;
  int    NBLAM_TOT[MXSPEC];    // total number of wavelength bins
  int    NBLAM_VALID[MXSPEC] ; // number of valid wavelength bins per epoch
  int    NBLAM_READ[MXSPEC];   // used by reader to check SPECTRUM_NLAM key
  double LAMRANGE_VALID[MXSPEC][2];  // used for print only

  int    IMJD_HOST, IS_HOST[MXSPEC];
  bool   SKIP[MXSPEC];        // outside Trest range of model

  double MJDREF_LIST[MXSPEC]; // same for each event
  double MJD_LIST[MXSPEC], TOBS_LIST[MXSPEC], TREST_LIST[MXSPEC];
  int    IMJD_NEARPEAK;   // imjd for SN spectrum closest to peak

  int    OPT_TEXPOSE_LIST[MXSPEC] ; // =1(user-fixe), =2(compute from SNR)
  double TEXPOSE_LIST[MXSPEC];
  double TEXPOSE_TEMPLATE ;            // passed from SIMLIB or computed
  float  SNR_REQUEST_LIST[MXSPEC] ;    // from user input
  float  SNR_COMPUTE_LIST[MXSPEC] ;    // actually computed
  double LAMOBS_SNR_LIST[MXSPEC][2];  // LAMOBS range for each SNR

  int    INDEX_TAKE_SPECTRUM[MXSPEC] ; // refers to INPUTS.TAKE_SPECTRUM
  int    INV_TAKE_SPECTRUM[MXSPEC]  ; // inverse of above

  // everything below gets malloc'ed by NBIN_SPECTRO

  // true flux vs [NMJD][ILAM], and smeared flux
  double *GENMAG_LIST[MXSPEC] ;
  double *GENSNR_LIST[MXSPEC] ;
  double *GENFLUX_LIST[MXSPEC] ;    // true flux with no noise or lam-smear
  double *GENFLUX_LAMSMEAR_LIST[MXSPEC]; // lam-smeared flux, no Poisson noise
  double *GENFLUX_PEAK ; // GENFLUX at PEAKMJD; needed for HOSTSNFRAC option
  double *GENMAG_PEAK ;  
  double  SCALE_FLAM_HOST_CONTAM[MXSPEC]; // fraction of host spec included in SN spec

  // observed (noisy) flux vs [NMJD][ILAM] 
  double  *OBSFLUX_LIST[MXSPEC] ;     // obs flux with noise
  double  *OBSFLUXERR_LIST[MXSPEC] ;  // measured error
  double  *OBSFLUXERRSQ_LIST[MXSPEC]; // square of error

  double  *GENFLAM_LIST[MXSPEC] ; // true dF/dlam [ispec][ilam]
  double  *FLAM_LIST[MXSPEC];     // measured dF/dlam
  double  *FLAMERR_LIST[MXSPEC];  // error on above
  double  *FLAMWARP_LIST[MXSPEC]; // warp applied to FLAM
  int      USE_WARP;

  double  GENMAG_SYNFILT[MXSPEC][MXFILTINDX]; // true synthetic mag per filter

  // items below are used for read utils (not used for sim)
  int     ID_LIST[MXSPEC] ;
  double *LAMMIN_LIST[MXSPEC], *LAMMAX_LIST[MXSPEC], *LAMAVG_LIST[MXSPEC] ; 

  double FLATRAN_LIST[100]; // 0-1 randoms, e.g., for pre-scales

  // define array of Gaussian randoms for noise
  double *RANGauss_NOISE_TEMPLATE ; 

  bool IS_MALLOC[MXSPEC] ;

} GENSPEC ;

// =============================================
// =========== FUNCTION PROTOTYPES =============
// =============================================

void init_spectrograph(char *inFile, char *stringOpt ) ;
void parse_spectrograph_options(char *stringOpt) ;
void read_spectrograph_text(char *inFile) ;
void read_spectrograph_fits(char *inFile) ;
void init_INPUTS_SPECTRO(void) ;
void extend_spectrograph_lambins(void);
void copy_INPUTS_SPECTRO(int ilam0, int ilam1);
void dump_INPUTS_SPECTRO(int nbin_dump, char *comment);
int  read_TEXPOSE_LIST(FILE *fp); 
int  read_SPECBIN_spectrograph(FILE *fp);
void reset_VALUES_SPECBIN(void) ;

void malloc_spectrograph(int OPT, int NBIN_LAM, int NBIN_TEXPOSE) ;
void solve_spectrograph(void) ;

void get_FILTERtrans_spectrograph(double *LMIN, double *LMAX, int MXLAM, 
				  int *NLAMBIN,
				  double *LAM_ARRAY, double *TRANS_ARRAY );

double getSNR_spectrograph(int ilam, double Texpose_S, double Texpose_T, 
			   bool ALLOW_TEXTRAP,double genMag,double *ERRFRAC_T);

void check_SNR_SPECTROGRAPH(int l, int t);

int  IMJD_GENSPEC(double MJD); // return IMJD index such that MJD_LIST[IMJD] = MJD

void create_ideal_spectrograph(double lammin, double lammax, double lambin );

// === END === 

