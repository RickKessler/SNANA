/***********
 Created Nov 2010 by R.Kessler

 
 parameters and global variables used to write/read 
 the GRID in fits format.

 Mar 16, 2016: add IFILTOBS[ifilt]
 Aug 27, 2017:  NON1A_NAME -> 200 chars instead of 40
 Sep 15, 2027:  MXGRIDGEN-> 500 (was 200)

**************/

// define internal version for backward-compatibility
//#define  IVERSION_GRID_WRITE 2    // Apr 3 2013
//#define  IVERSION_GRID_WRITE 3    // Aug 30 2013: add WGT and MAGOFF
#define  IVERSION_GRID_WRITE 4    // Jan 08 2017: add ITYPE_USER

#define MXGRIDGEN       500       // max number of grid point in any dimension
#define MAXMAG_GRIDGEN  32.0         // max mag to store with 16 bits
#define MAGPACK_GRIDGEN  1000.0       // I2MAG = MAG * MAGPACKSCALE
#define FLUXPACK_GRIDGEN 20000.      // for snoopy fluxes
#define NPADWD_LCBEGIN  2            // # pad words for each light curve
#define NPADWD_LCEND    2            // # pad words for each light curve
#define MARK_GRIDGEN_LCBEGIN -1111    // begin LC marker
#define MARK_GRIDGEN_LCEND   -9999    // end of LC marker

#define IPAR_GRIDGEN_LOGZ      1   // log10(redshift)
#define IPAR_GRIDGEN_COLORPAR  2   // AV or SALT2c
#define IPAR_GRIDGEN_COLORLAW  3   // RV or beta
#define IPAR_GRIDGEN_SHAPEPAR  4   // x1, Delta, dm15, nonIa-index ..
#define IPAR_GRIDGEN_FILTER    5   // filter index
#define IPAR_GRIDGEN_TREST     6   // rest-frame epoch
#define NPAR_GRIDGEN           6

#define EXTNAME_SNPAR_INFO   "SNPAR-INFO"
#define EXTNAME_LOGZ         "LOGZ-GRID" 
#define EXTNAME_COLORPAR     "COLOR-GRID" 
#define EXTNAME_COLORLAW     "RV/BETA-GRID" 
#define EXTNAME_SHAPEPAR     "SHAPE-GRID" 
#define EXTNAME_FILTER       "FILTER-GRID" 
#define EXTNAME_TREST        "TREST-GRID" 
#define EXTNAME_NONIa_INFO   "NONIa-INFO" 
#define EXTNAME_PTRI2LCMAG   "PTR_I2LCMAG" 
#define EXTNAME_I2LCMAG      "I2LCMAG" 
char    EXTNAME_GRIDGEN[NPAR_GRIDGEN+1][40] ;

#define SNTYPE_GRIDGEN_Ia     1
#define SNTYPE_GRIDGEN_NONIa  2

#define OPT_GRIDGEN_FORMAT_TEXT 1
#define OPT_GRIDGEN_FORMAT_FITS 2
int     OPT_GRIDGEN_FORMAT ;

FILE     *fp_GRIDGEN_TEXT ;
fitsfile *fp_GRIDGEN_FITS ;

#define NCOMMENT_GRIDGEN  8
char COMMENT_GRIDGEN[NCOMMENT_GRIDGEN][80];

// xxx mark delete char  GRIDGEN_FILE[200];   // copied from AUX_SIMFILE.GRIDGEN

float GRIDGEN_I2LCPACK;   // stored I*2 value = mag(or flux) x this value


// define user inputs
struct GRIDGEN_INPUTS {
  int  NBIN[NPAR_GRIDGEN+1] ;
  char FORMAT[20];  // from sim-input file: "TEXT" or "FITS"
} GRIDGEN_INPUTS ;



typedef struct {

  int   IVERSION ;          // internal version
  char  UNIQUE_KEY[100];

  char  SURVEY[60];
  char  MODEL[MXPATHLEN];
  char  FILTERS[MXFILTINDX] ;
  int   IFILTOBS[MXFILTINDX] ;  // Mar 2016, computed on readback only

  // start with info vs. grid parameter
  int   NBIN[NPAR_GRIDGEN+1];                // Number of bins
  float BINSIZE[NPAR_GRIDGEN+1] ;            // bin size
  float VALUE[NPAR_GRIDGEN+1][MXGRIDGEN] ;   // value at each grid point
  char  NAME[NPAR_GRIDGEN+1][20];        // name of param (i.e, DELTA, AV)
  float VALMIN[NPAR_GRIDGEN+1] ;   // same as VALUE[1]
  float VALMAX[NPAR_GRIDGEN+1] ;   // same as VALUE[NBIN]

  float FILTER_LAMAVG[MXGRIDGEN];  // mean wavelength for each filter

  // ilc = 1 + sum ILCOFF * (indx[IPAR] -1)
  int ILCOFF[NPAR_GRIDGEN+1] ; // to get abs LC index from iz,ic,ilum, etc ...

  int   *PTR_VALUE[NPAR_GRIDGEN+1];    // VALUE[ilc] = value for this 'ilc'
  float CURRENT_VALUE[NPAR_GRIDGEN+1];   // current par values for each LC

  // now the LC info
  int NGRIDGEN_LC  ;    // total number of lightcurves to generate
  int NGRIDGEN_PER_LC ; // number of measures per LC

  long int SIZEOF_GRIDGEN ; // total size of GRIDGEN (excluding header)
  int      NWD_I2GRIDGEN;

  short *I2GRIDGEN_LCMAG ;
  short *I2GRIDGEN_LCERR ;
  int   *PTR_GRIDGEN_LC ;

  // extra NONIA-INFO from simgen-input file
  int   NON1A_INDEX[MXGRIDGEN];     // SNANA index vs. sparse nonIa index
  char  NON1A_NAME[MXGRIDGEN][200];  // full SN name vs. idem
  char  NON1A_CTYPE[MXGRIDGEN][20]; // string Type(Ib,II..) vs. idem
  int   NON1A_ITYPE_AUTO[MXGRIDGEN];  // 1=Ib,Ic,Ibc, etc ... 2=II, IIP, IIL..
  int   NON1A_ITYPE_USER[MXGRIDGEN];  // SNTAG = user type (Jan 2017)
  float NON1A_WGT[MXGRIDGEN];       // relative rate (Aug 30 2013)
  float NON1A_WGTSUM[MXGRIDGEN];    // sumulative sum
  float NON1A_MAGOFF[MXGRIDGEN];    // mag offset
  float NON1A_MAGSMEAR[MXGRIDGEN] ; // mag smear in sim (not used to make grid)

  // Aug 2016: define info for PEC1A to allow different z-dependent rate
  int   ISPEC1A[MXGRIDGEN];         // flag peculiar Ia (Aug 2016)
  float FRAC_PEC1A ;

} SNGRID_DEF ;

SNGRID_DEF  SNGRID_WRITE ; // used by sim to write GRID

int OPT_SNOOPY_FLUXPACK ;

int NROW_WRITE_TOT ;

// ==========================


// ===================================
//     GRID  function declarations
// ===================================

// grid-write utils
void   init_GRIDsource(int opt);     // init stuff for GRID
void   init0_GRIDsource(void);   
void   init1_GRIDsource(void);    // init after genmodel is init'ed

void   gen_GRIDevent(int ilc);
void   update_GRIDarrays(void);

void   wr_GRIDfile(int OPT, char *GRIDfile);
void   wrhead_GRIDfile_text(void);
void   wrhead_GRIDfile_fits(void);
void   append_GRIDfile_text(void);
void   append_GRIDfile_fits(void);

void   end_GRIDfile(void);
void   get_GRIDKEY(void) ;

int    SNTYPE_GRIDGEN(void); 
void   load_EXTNAME_GRIDGEN(void);  // load EXTNAME_GRIDGEN array

// read back utils
void   load_EXTNAME_GRIDREAD(int IVERSION); // idem for different ifdef 
void   fits_read_SNGRID(int OPTMASK, char *sngridFile, 
			SNGRID_DEF *SNGRID); 

void   check_fitserror(char *comment, int status);

int  INDEX_GRIDGEN(int ipar, double parval, SNGRID_DEF *SNGRID);
int  get_NON1A_ITYPE_SNGRID(char *NON1A_CTYPE );

void renorm_wgts_SNGRID(SNGRID_DEF *SNGRID ) ;
void sort_by_PEC1A_SNGRID(SNGRID_DEF *SNGRID ) ;
void copy_NON1A_SNGRID(int IROW1, int IROW2, 
		       SNGRID_DEF *SNGRID1, SNGRID_DEF *SNGRID2 );

void dump_SNGRID(SNGRID_DEF *SNGRID ) ;

// ============= END ==========
