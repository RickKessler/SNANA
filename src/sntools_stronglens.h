
// ?? #define MXMAP_STRONGLENS 100000

#define zMIN_STRONGLENS  0.1 // avoid problems searching library

struct {
  int USE_FLAG; // logical flag for simulation
  int NCALL;

  // contents from model file/library
  int   NLENS;     // number of lenses in library
  int   *IDLENS;  // ID for each lens
  float *ZLENS;  // lens redshift
  float *LOGMASS_LENS, *LOGMASS_ERR_LENS;
  float *ZSRC; //source redshift
  int   *NIMG; // Number of images per lens
  float **XIMG_SRC, **YIMG_SRC ; // X and Y offsets of SN, arcsec
  float **XGAL_SRC, **YGAL_SRC ; // idem for host center
  float **DELAY ; // time delay of each image (days)
  float **MAGNIF; //magnification of each image
  
  char VARNAME_LENSID[40];
  char VARNAME_ZSRC[40];
  char VARNAME_ZLENS[40];
  char VARNAME_LOGMASS_LENS[40];
  char VARNAME_LOGMASS_ERR_LENS[40];
  char VARNAME_NIMG[40];
  char VARNAME_XIMG_SRC[40];
  char VARNAME_YIMG_SRC[40];
  char VARNAME_XGAL_SRC[40];
  char VARNAME_YGAL_SRC[40];
  char VARNAME_MAGNIF[40];
  char VARNAME_DELAY[40];

  int ICOL_LENSID, ICOL_ZSRC, ICOL_ZLENS;
  int ICOL_LOGMASS_LENS, ICOL_LOGMASS_ERR_LENS; // Jun 30 2022
  int ICOL_NIMG;
  int ICOL_XIMG_SRC, ICOL_YIMG_SRC;
  int ICOL_XGAL_SRC, ICOL_YGAL_SRC;
  int ICOL_MAGNIF, ICOL_DELAY;

} INPUTS_STRONGLENS;

typedef struct {                                                                  
  long long int IDLENS; 
  int    NIMG; 
  double zLENS, zSRC_MATCH ;  
  double LOGMASS, LOGMASS_ERR ;

  // store lists of image-dependent quantities 
  double XIMG_SRC_LIST[MXIMG_STRONGLENS], YIMG_SRC_LIST[MXIMG_STRONGLENS];  
  double XGAL_SRC_LIST[MXIMG_STRONGLENS], YGAL_SRC_LIST[MXIMG_STRONGLENS]; 
  double DELAY_LIST[MXIMG_STRONGLENS];    
  double MAGNIF_LIST[MXIMG_STRONGLENS], MAGSHIFT_LIST[MXIMG_STRONGLENS] ;  
} EVENT_STRONGLENS_DEF ;  

// function declarations
void init_stronglens(char *MODEL_FILE);
void malloc_stronglens(int NLENS);

void get_stronglens(double zSN, double *hostpar, int LDMP,
		    EVENT_STRONGLENS_DEF *EVENT_STRONGLENS) ;

double prob_stronglens(double z);

// end:
