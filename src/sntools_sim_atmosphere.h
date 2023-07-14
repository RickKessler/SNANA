
// Created Jun 2023

#define COORD_SHIFT_NULL_ARCSEC 99.0
#define COORD_SHIFT_NULL_DEG    99.0/3600.0 // 99 arcsec  

// DCR effects in Table 1 of arXiv:2304.01858                                                                  
#define ATMOSPHERE_OPTMASK_DCR_COORD         1
#define ATMOSPHERE_OPTMASK_DCR_PSFSHAPE      2
#define ATMOSPHERE_OPTMASK_SIMGEN_DUMP_DCR   512 // write DCR SIMGEN DUMPfile                                 

#define KEYNAME_ATMOSPHERE_DCR_COORDRES_POLY   "ATMOSPHERE_DCR_COORDRES_POLY"
#define KEYNAME_ATMOSPHERE_DCR_MAGSHIFT_POLY   "ATMOSPHERE_DCR_MAGSHIFT_POLY"


// define inputs read from sim-input file
struct {
  int OPTMASK;
  bool DO_DCR_COORD ;
  bool DO_DCR_PSFSHAPE ;

  char SEDSTAR_FILE[MXPATHLEN]; // stellar SED to compute <lam> in each passband
  GENPOLY_DEF  DCR_COORDRES_POLY; // poly fun for astrometry resolution vs. PSF/SNR
  GENPOLY_DEF  DCR_MAGSHIFT_POLY; // poly fun for mag shift vs. PSF-fraction shift

  // define Gaussian variations in T,BH,PWV that impact index of refraction
  double SIGMA_SITE_TEMP; // Gauss sigma for site temperaure variation (Celsius)
  double SIGMA_SITE_BP;   // Gauss sigma for site barometri pressure variation (mm Hg)
  double SIGMA_SITE_PWV;  // Gauss sigma for site PWV variation (mm Hg)
  bool   APPLY_SIGMA_SITE; // true of 1 or more sigmas is non-zero

} INPUTS_ATMOSPHERE ;



typedef struct {
  double AVG;
  double SUM, WGTSUM ;
  double AVG_BAND[MXFILTINDX];
  double SUM_BAND[MXFILTINDX], WGTSUM_BAND[MXFILTINDX];
} COORD_AVG_DEF ;

struct {

  int    NBINLAM_CALSTAR;
  double *LAM_ARRAY_CALSTAR, *FLUX_ARRAY_CALSTAR ; // full SED
  double  LAMAVG_CALSTAR[MXFILTINDX];  // <lam> per band for calib stars 
  double *FLUX_CALSTAR[MXFILTINDX];    // flux-vs-lam on filter-lam grid
  double  n_CALSTAR_AVG[MXFILTINDX];   // index of refrac per band, calib stars

  double PRESSURE_AVG, TEMPERATURE_AVG, PWV_AVG;  // at telescope site location


  double SNRMIN; // min SNR to include in RA/DEC avg

  // info per event
  // GENLC struct in snlc_sim.h stores RA/DEC per obs;
  // here keep track of avg coord per band and over all obs.
  COORD_AVG_DEF COORD_RA;
  COORD_AVG_DEF COORD_DEC;
  COORD_AVG_DEF COORD_SIM_RA;
  COORD_AVG_DEF COORD_SIM_DEC;

  double COORDRES[MXEPSIM]; // computed coord resolution (asec) from measured SNR

} ATMOS_INFO ;

// ------- function prototypes ----------

void  GEN_ATMOSPHERE_DRIVER(void) ;
void  INIT_ATMOSPHERE(void);
void  INIT_EVENT_ATMOSPHERE(void) ;

void   read_stellar_sed_atmos(void) ;
void   init_stellar_sed_atmos(int ifilt_obs) ;

void reset_COORD_AVG(COORD_AVG_DEF *COORD);
void sum_COORD_AVG(COORD_AVG_DEF *COORD,
                   double RA_OBS, double WGT, int IFILT_OBS);

void  gen_airmass(int ep);
void  genSmear_coords(int ep);
void  gen_dcr_coordShift(int ep);

void  gen_dcr_magShift(int ep);

double compute_DCR_angle_approx(double LAM, double tan_ZENITH, 
				int IFILT_OBS, int DUMPFLAG); // xxx mark delete

double compute_DCR_angle(int ep, int DUMPFLAG);

double compute_index_refrac_atmos(double LAM, int DUMPFLAG) ;

void   test_compute_dcr(void);

