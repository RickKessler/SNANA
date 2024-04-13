
#define MXBIN_LENSING_z   100
#define MXBIN_LENSING_dmu 6000

struct {
  int USEFLAG;
  int NBIN_z, NBIN_dmu;
  double *z_LIST, *dmu_LIST ;
  double **PROB;      // PROB[iz][dmu]

  double **FUNPROB ;  // DMU(SUMprob) in iz,imu bins
  double **FUNDMU ;

  double zMIN, zMAX, dmuMIN, dmuMAX;
} LENSING_PROBMAP ;

void   init_lensDMU(char *mapFileName, float dsigma_dz) ;
double gen_lensDMU(double z, double ran1, int DUMP_FLAG);
