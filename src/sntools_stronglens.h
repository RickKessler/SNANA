
#define MXMAP_STRONGLENS 100000
#define MXIMG_STRONGLENS 8   // max number of images per lens


struct {
  int USE_FLAG; // logical flag for simulation
  int NCALL;

  // Justin to do: Nimage => NIMG, Ximg,Yimg => XIMG, YIMG 

  // contents from model file/library
  int NLENS;     // number of lenses in library
  int    *IDLENS;  // ID for each lens
  float  *zLENS;  // lens redshift
  float  *zSRC; //source redshift
  int    *Nimage; // Number of images per lens
  float **Ximg, **Yimg ; // X and Y offsets, arcsec
  float **tdelay ; // time delay of each image (days)
  float **mu; //magnification of each image
  

  char VARNAME_LENSID[40];
  char VARNAME_zSRC[40];
  char VARNAME_zLENS[40];
  char VARNAME_NIMG[40];
  char VARNAME_XIMG[40];
  char VARNAME_YIMG[40];
  char VARNAME_MAG[40];
  char VARNAME_DELAY[40];

  int ICOL_LENSID, ICOL_zSRC, ICOL_zLENS, ICOL_NIMG, ICOL_XIMG, ICOL_YIMG;
  int ICOL_MAG, ICOL_DELAY;

} INPUTS_STRONGLENS;


// function declarations

void init_stronglens(char *MODEL_FILE);
void malloc_stronglens(int NLENS);

void get_stronglens(double zSN, double *hostpar, int *IDLENS, double *zLENS, 
		    int *blend_flag, int *Nimage,
		    double *tdelay, double *mu, double *Ximg, double *Yimg);

// end:
