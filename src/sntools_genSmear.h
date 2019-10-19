// sntools_genSmear.h
//
// Mar 30 2018: MXLAM_GENSMEAR_SALT2 --> 4000 (was 1000)
// Oct 18 2019: add COV model

void  init_genSmear_FLAGS(double SCALE); // set flags to zero
int   istat_genSmear(void) ;

void  init_genSmear_USRFUN(int NPAR, double *parList, double *LAMRANGE ) ;

void   init_genSmear_SALT2(char *version, char *dispFile, double SIGCOH, 
			   double *GENRANGE_REDSHIFT);
void   read_genSmear_SALT2disp(char *smearFile) ;
void   read_genSmear_SALT2sigcoh(char *versionSALT2, GRIDMAP1D *SIGCOH_LAM ) ;
void   parse_SIGCOH_SALT2(char *KEYNAME, char *KEYARG, GRIDMAP1D *SIGCOH_LAM);
void   getFileName_SALT2colorDisp(char *fileName) ; // added Jan 2017

void  init_genSmear_Chotard11(int OPT_farUV) ;
void  init_genSmear_VCR(char *VCR_version, int index_SNmodel);
void  init_genSmear_CCM89(double *LAMRANGE) ;
void  init_genSmear_COH(void) ;
void  init_genSmear_biModalUV(void) ;
void  init_genSmear_OIR(void);
void  init_genSmear_COV(int OPTMASK, char *covFileName);

void  init_genSmear_private(void) ;

int   nval_genSmear_override(char *inputKey, char *parName);
void  store_genSmear_override(char *parName, int NVAL, double *tmpList);
int   exec_genSmear_override(int ipar, char *keyName, double *val) ;

void  read_VCR_VSI(void) ;     // read v_Si distribution
void  read_VCR_INFO(void) ;    // read misc. info/parameters
void  prep_VCR_forSALT2(void); // inits specific to SALT2
void  sort_VCR_BANDS(void);
void  prep_VCR_COVAR(void) ;
void  get_VCR_colorOffset(int ic, double *OFFSET, double *RMS );
void  parse_VCR_colorString(int ic);

void  get_NRAN_genSmear(int *NRANGauss, int *NRANFlat); // returns NRANxxx

void  SETRANGauss_genSmear(int NRAN, double *ranList );
void  SETRANFlat_genSmear(int NRAN, double *ranList );
void  SETSNPAR_genSmear(double shape, double color, double redshift) ;

void get_genSmear(double Trest, int NLam, double *Lam, double *magSmear) ;

void  get_genSmear_USRFUN(double Trest, int NLam, double *Lam, 
			  double *magSmear ) ;

void  get_genSmear_SALT2(double Trest, int NLam, double *Lam, 
			 double *magSmear ) ;

void  get_genSmear_Chotard11(double Trest, int NLam, double *Lam, 
			     double *magSmear ) ;

void  get_genSmear_VCR(double Trest, int NLam, double *Lam, 
		       double *magSmear ) ;


double get_CCM89_RV(void) ;
void  get_genSmear_CCM89(double Trest, int NLam, double *Lam, 
			 double *magSmear ) ;

void get_genSmear_COH(double Trest, int NLam, double *Lam, 
		      double *magSmear) ;

void  get_genSmear_biModalUV(double Trest, int NLam, double *Lam, 
			    double *magSmear ) ;

void  get_genSmear_OIR(double Trest, int NLam, double *Lam, 
		       double *magSmear ) ;

void  get_genSmear_COV(double Trest, int NLam, double *Lam, 
		       double *magSmear ) ;

void  get_genSmear_private(double Trest, int NLam, double *Lam, 
			   double *magSmear ) ;

int INODE_LAMBDA(double LAM, int NNODE, double *LAM_NODES);

void extraFilters_4genSmear(char *modelName,
			    int *NFILT_extra, int *IFILTOBS_extra) ;

#define MXBAND_GENSMEAR 20
void genSmear_nodes(int NBAND, double *LAMCEN, 
		    double CholMatrix[MXBAND_GENSMEAR][MXBAND_GENSMEAR],
		    int NLam, double *Lam, double *magSmear ) ;

// ============= GLOBALS ================

#define MXRAN_GENSMEAR  500  // 100->500 Oct 18 2019

// define central wavelength for commonly used Bessell filters 
#define LAMCEN_FILT_v   2500.0
#define LAMCEN_FILT_U   3560.0
#define LAMCEN_FILT_B   4390.0
#define LAMCEN_FILT_V   5490.0
#define LAMCEN_FILT_R   6545.0
#define LAMCEN_FILT_I   8045.0

double MAGSMEAR_COH;   // coherent part of scatter only, for SNTABLE


struct GENSMEAR {
  int    NGEN_RANGauss ;
  int    NGEN_RANFlat ;
  int    NSET_RANGauss ;
  int    NSET_RANFlat ;
  double RANGauss_LIST[MXRAN_GENSMEAR] ;
  double RANFlat_LIST[MXRAN_GENSMEAR] ;

  double SHAPE, COLOR ; // allows more complex magSmear models
  double REDSHIFT ;       // allows for redshift evolution (Jan 2014)
  double ZPOW[10];  // store powers of redshift for faster calculations.

  double SCALE ; // global scale to all mag-smearing (Oct 9 2018)

} GENSMEAR ;

// ---------- USRFUN struct -----------------
#define NPAR_GENSMEAR_USRFUN_REQUIRED   8
struct GENSMEAR_USRFUN {
  int     USE ;
  double  SIGCOH ; // sigma for coherent scatter  
  double  A5500 ;
  double  LAMINTERVAL;
  double  LAMPHASE ;
  double  TAU_LAM ; 
  double  TAU_EPOCH ;
  double  CORR_LAM ;     // wavelength correlation, Angstroms
  double  CORR_EPOCH ;   // epoch correlation, days

  double LAMREF ;  // Ampl defined here
  double LAMRANGE[2];

  // define temp variables in case we need to debug
  int    NNODE;
  double LAM_NODE[MXRAN_GENSMEAR]  ; // wavelength at each node
  double AMPL_NODE[MXRAN_GENSMEAR] ; // Amplitude at each node

}  GENSMEAR_USRFUN ;


// ---------- G10 struct -----------------
#define MXLAM_GENSMEAR_SALT2 4000

struct GENSMEAR_SALT2 {
  int    USE ;
  char   FILE[200] ;  // mag smear vs. wavelength
  int    NLAM ;       // no. lambda bins defining color smearing
  double LAM[MXLAM_GENSMEAR_SALT2];   // lambda at each bin
  double MINLAM, MAXLAM;              // min,max wavelenght
  double SIGMA[MXLAM_GENSMEAR_SALT2]; // sigma-smear value, mags
  double SIGMA_SCALE ;

  // SIGMA_INT from SALT2.INFO file or from user input
  GRIDMAP1D SIGCOH_LAM ;   // May 30 2018 : sigma_coh vs. wavelenth
  double    SIGCOH ;       // traditional coherent term 

  double LAMSEP_NODE ;
  double LAM_NODE[MXRAN_GENSMEAR];
  double SIG_NODE[MXRAN_GENSMEAR];
  int    NNODE ;

  double FUDGE_dSmear_dLAM[2] ; // Fudge G10 curve to match broadband smear
  double FUDGE_LAMSWITCH  ;

} GENSMEAR_SALT2 ;


// ---------- Chotard11 struct -----------------

#define NBAND_C11 6  // UBVRI correlations
struct GENSMEAR_C11 {
  int USE ;
  double Cholesky[NBAND_C11][NBAND_C11] ;
  int OPT_farUV;  // see sub-models C11_0, C11_1, C11_2
} GENSMEAR_C11 ;


struct GENSMEAR_OIR {
  int USE ;
} GENSMEAR_OIR;

struct GENSMEAR_COV {
  int USE ;
  int OPTMASK;
} GENSMEAR_COV;


// ------------- CCM89 struct ---------------
struct GENSMEAR_CCM89 {
  int USE ;
  double SIGCOH ;      // coherent term
  double RV_REF;       // reference RV; RV=RV_REF -> no fluctuation
  double RV_sigma ;    // sigma for random RV
  double RV_range[2] ; // generate withing this range

  int     NLAM_MAP ;
  double *LAM_MAP ;
  double *XTDIF_MAP;  // XT-XT0 vs. LAM at current RV

  double RV ;          // randome RV for each SN
} GENSMEAR_CCM89;


#define MXCOLOR_VCR  10

struct GENSMEAR_VCR {
  int      USE ;

  // path and filenames
  char MODELPATH[MXPATHLEN] ;
  char INFO_FILE[MXPATHLEN] ;
  char VSI_FILE[MXPATHLEN];

  // inputs from VCR.INFO file
  int      NCOLOR ;
  char     COLOR_STRING[MXCOLOR_VCR][8];  
  double   COLOR_SLOPE[MXCOLOR_VCR] ;
  double   SIGMACOH_MB ;
  double   COLOR_SIGMA_SCALE ; // multiplies all the COLOR_SIGMA values.
  double   COLOR_SIGMA[MXCOLOR_VCR] ;
  double   COLOR_CORMAT[MXCOLOR_VCR][MXCOLOR_VCR] ;

  double   LAMCEN_BAND[MXCOLOR_VCR*2];  // define warp nodes
  double   SPECSHIFT_SCALE;  // fudge-scale for spectral shifts

  int      USE_ZVAR ;  // logical flag
  double   ZVARIATION_POLY_VSI[4]; // VSI += sum [z^i * coef[i] ]

  // VSI distribution read from file
  int      NBIN_VSI; // number of bins for v_Si distribution
  double  *VSI, *PROB ;  // v_Si prob distribution
  double  *SUMPROB ;     // cumulative prov (for efficiency MC generation)
  double   VSI_MEAN ;

  // misc.
  double    COLOR_OFFSET[MXCOLOR_VCR] ;
  double    COLOR_RMS[MXCOLOR_VCR] ;
  char      COLOR_BAND[MXCOLOR_VCR][2][2];    // each band
  int       COLOR_IFILTDEF[MXCOLOR_VCR][2] ;  // absolute filter indices
  double    COLOR_COVAR2[MXCOLOR_VCR][MXCOLOR_VCR] ; // covariances
  double    COLOR_COVAR1[MXCOLOR_VCR*MXCOLOR_VCR] ;  // covariances
  double    Cholesky[MXCOLOR_VCR][MXCOLOR_VCR] ; 
  int       COLOR_IFILT[MXCOLOR_VCR] ;  // points to ORDERED_IFILTDEF

  // define list of filters ordered by wavelength
  int       NFILTDEF ;  // number of unique filters = size of ordered list
  int       ORDERED_IFILTDEF[MXCOLOR_VCR*2];  // sparse/ordered list
  double    ORDERED_LAMCEN[MXCOLOR_VCR*2];    // sparse/ordered list
  int       LAMCEN[MXFILTINDX] ;  // <LAM> vs. IFILTDEF
  int       IFILT_B ;

  // generated quantieis for SIMGEN-DUMP file
  double   GENRAN_VSI ;  // generated VSI
  double   GENRAN_COLORSHIFT[MXCOLOR_VCR];  // broadband color shift

} GENSMEAR_VCR ;



struct GENSMEAR_COH {
  int    USE ;
  double MAGSIGMA ;
} GENSMEAR_COH ;



struct GENSMEAR_BIMODAL_UV {
  int     USE ;
  double  SIGMACOH ; 
  double  LAMU_MIN, LAMU_CEN, LAMU_MAX;
  double  MAGU_SPLIT, MAGU_SIGMA ;
} GENSMEAR_BIMODAL_UV ;

// --------- private struct for testing --------------
struct GENSMEAR_PRIVATE {
  int USE ;
} GENSMEAR_PRIVATE ;


int NUSE_GENSMEAR ;
int NCALL_GENSMEAR ;


// --------------------------------------------------------------------
// Override params passed from simulation to allow command-line overrides
// See sim-input key GENMAG_SMEARPAR_OVERRIDE (Mar 22 2014)

#define MXSMEARPAR_OVERRIDE 20
int NSMEARPAR_OVERRIDE ;
struct {
  int    NVAL ;    // number of values passed per SMEARPAR (default=1)
  char   NAME[60]; // name of parameter to override
  double VALUE[MXSMEARPAR_OVERRIDE] ; // list of value(s)
} GENMAG_SMEARPAR_OVERRIDE[MXSMEARPAR_OVERRIDE] ;

