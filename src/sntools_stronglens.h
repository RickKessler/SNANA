
#define MXMAP_STRONGLENS 100000
#define MXIMG_STRONGLENS 8   // max number of images per lens

struct {
  int USE_FLAG;
  int NCALL;

  // contents from model file/library
  int NLENS;     // number of lenses in library
  int    *IDLENS;  // ID for each lens
  float  *zLENS;  // lens redshift
  int    *Nimage; // Number of images per lens
  float **angSep, **phi ; // angle sep (arcsec), and phi (degrees)

} INPUTS_STRONGLENS;


// function declarations

void init_stronglens(char *MODEL_FILE);
void malloc_stronglens(int NLENS);

int get_stronglens(double zSN, double *hostpar, double *zLENS, 
		   int *blend_flag, int *Nimage,
		   double *mu, double *angSep, double *phi);

// end:
