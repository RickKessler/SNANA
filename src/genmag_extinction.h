
// Globals and functions for host-galaxy extinction (XT)

#define OPT_SNXT_CCM89  1  // use exact CCM89 model to apply host extinc  
#define OPT_SNXT_SJPAR  2  // use Saurabh alpha,beta,zeta paramitrization

// hard wire RVinv range and bin size for GRIDMAP
#define RVinv_MIN_XTMAG  0.25
#define RVinv_MAX_XTMAG  2.00
#define RVinv_BIN_XTMAG  0.25

#define OPT_MWCOLORLAW_XTMAG 94 // ODonnel 94 update to CCM89

double **TEMP_XTMAG_ARRAY;

struct {
  int OPT_SNXT;

  // use binInfo typedef from sntools_calib.h
  KCOR_BININFO_DEF BININFO_Trest;  
  KCOR_BININFO_DEF BININFO_RVinv;  

  struct GRIDMAP GRIDMAP3D ;

} XTMAG_INFO ;

void init_genmag_extinction(int OPT_SNXT);
void init_genmag_extinction__(int *OPT_SNXT);

void fill_binInfo_XTMAG(void);
void fill_GRIDMAP3D_XTMAG(void);
double eval_XTMAG_AV1(int ifilt, double Trest, double RVinv) ;
void dump_XTMAG(double Trest) ;

double genmag_extinction(int ifiltdef, double Trest, double RV, double AV );
double genmag_extinction__(int *ifiltdef, double *Trest, 
			   double *RV, double *AV );

