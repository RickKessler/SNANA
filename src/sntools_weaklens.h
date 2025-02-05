

#define MXBIN_LENSING_z   100
#define MXBIN_LENSING_dmu 6000

struct {
  char   PROBMAP_FILE[MXPATHLEN];
  float  DMUSCALE;            // scale width of DMU profile 
  float  DMUERR_FRAC;         // frac error on lensdmu (Feb 2025)
  float  DSIGMADZ ;           // symmetric Gaussian model (not recommended)
} INPUTS_WEAKLENS ;


struct {
  int USEFLAG;
  int NBIN_z, NBIN_dmu;
  double *z_LIST, *dmu_LIST ;
  double **PROB;      // PROB[iz][dmu]

  double **FUNPROB ;  // DMU(SUMprob) in iz,imu bins
  double **FUNDMU ;

  double zMIN, zMAX, dmuMIN, dmuMAX;
} LENSING_PROBMAP ;


// xxx mark void   init_lensDMU(char *mapFileName, float dsigma_dz) ;
void   init_lensDMU(void) ;
double gen_lensDMU(double z, double ran1, int DUMP_FLAG);
double gen_lensDMU_smear(double lensDMU);
