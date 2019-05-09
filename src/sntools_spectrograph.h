
// Created July 2016 by R.Kessler

#define MXSPEC                 50   // max Nspec per event   
#define MXTEXPOSE_SPECTROGRAPH 50   // max size of TEXPOSE grid
#define MXLAM_SPECTROGRAPH    2000  // same as MXBIN_LAMFILT_SEDMODEL
#define MXLAMSMEAR_SPECTROGRAPH 20  // max number of smeared lambda bins
#define FITSTABLE_NAME_SPECTROGRAPH  "SPECTROGRAPH" 
#define NCOL_noSNR 3    // Ncolumns before SNR values
#define TEXPOSE_INFINITE_SPECTROGRAPH 1.0E8 // flag to ignore noise
#define ILIST_RANDOM_SPECTROGRAPH 3  // Jan 2018

int  SPECTROGRAPH_USEFLAG ;
int  NERR_SNR_SPECTROGRAPH ;

struct {

  char   INFILE_NAME[MXPATHLEN]; // input table with SPECBIN keys
  char   INSTRUMENT_NAME[60];
  double MAGREF_LIST[2];  // two MAGREFs to compute ZP and SIGSKY
  int    NBIN_LAM, NBIN_LAM_noREBIN;
  int    NBIN_TEXPOSE ;
  int    NREBIN_LAM;      // rebin by this integer value (Sep 2 2016)

  // define spectro table read from input file: 
  // needs to be malloced with NBIN_xxx
  double LAM_MIN,     LAM_MAX; // global min & max lambda
  double TEXPOSE_MIN, TEXPOSE_MAX; // global min & max TEXPOSE

  double *LAMMIN_LIST, *LAMMAX_LIST, *LAMAVG_LIST ;     // per spectro bin
  double *LAMBIN_LIST, *LAMSIGMA_LIST ; // binSize & resolution vs. lambda
  double TEXPOSE_LIST[MXTEXPOSE_SPECTROGRAPH] ;    // per expTime bin
  double **SNR0, **SNR1  ;  // per spectro bin & expTime

  int FORMAT_MASK ;  // 1=>LAMCEN only, 2=LAMMIN and LAMMAX

  // computed from table inputs
  double **ZP, **SQSIGSKY ; // per spectro & expTime bin

  // syn-filter info from kcor-iput file; filled when reading kcor file
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
  int    NBLAM_TOT;            // total number of wavelength bins
  int    NBLAM_VALID[MXSPEC] ; // number of valid wavelength bins per epoch

  double MJDREF_LIST[MXSPEC]; // same for each event
  double MJD_LIST[MXSPEC], TOBS_LIST[MXSPEC], TREST_LIST[MXSPEC];

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

  // observed (noisy) flux vs [NMJD][ILAM] 
  double  *OBSFLUX_LIST[MXSPEC] ;     // obs flux with noise
  double  *OBSFLUXERR_LIST[MXSPEC] ;  // measured error
  double  *OBSFLUXERRSQ_LIST[MXSPEC]; // square of error

  double  *GENFLAM_LIST[MXSPEC] ; // true dF/dlam
  double  *FLAM_LIST[MXSPEC];     // measured dF/dlam
  double  *FLAMERR_LIST[MXSPEC];  // error on above
  double  *FLAMWARP_LIST[MXSPEC]; // warp applied to FLAM
  int      USE_WARP;

  // define array of Gaussian randoms for noise
  double *RANGauss_NOISE_TEMPLATE[MXLAMSMEAR_SPECTROGRAPH] ; 

} GENSPEC ;

// =============================================
// =========== FUNCTION PROTOTYPES =============
// =============================================

void init_spectrograph(char *inFile, char *stringOpt ) ;
void parse_spectrograph_options(char *stringOpt) ;
void read_spectrograph_text(char *inFile) ;
void read_spectrograph_fits(char *inFile) ;
int  read_TEXPOSE_LIST(FILE *fp); 
int  read_SPECBIN_spectrograph(FILE *fp);
void reset_VALUES_SPECBIN(void) ;

void malloc_spectrograph(int OPT, int NBIN_LAM, int NBIN_TEXPOSE) ;
void solve_spectrograph(void) ;

void get_FILTERtrans_spectrograph(double *LMIN, double *LMAX, int MXLAM, 
				  int *NLAMBIN,
				  double *LAM_ARRAY, double *TRANS_ARRAY );

double getSNR_spectrograph(int ilam, double Texpose_S, double Texpose_T, 
			   double genMag,  double *ERRFRAC_T );

int check_SNR_SPECTROGRAPH(int l, int t);

// === END === 

