
// Created Jun 2023

#define COORD_SHIFT_NULL_ARCSEC 99.0
#define COORD_SHIFT_NULL_DEG    99.0/3600.0 // 99 arcsec  

typedef struct {
  double AVG;
  double SUM, WGTSUM ;
  double AVG_BAND[MXFILTINDX];
  double SUM_BAND[MXFILTINDX], WGTSUM_BAND[MXFILTINDX];
} COORD_AVG_DEF ;

struct {
  // fixed info that doesn't change
  double PRESSURE, TEMPERATURE, PWV;  // at telescope site location
  double LAMAVG_CALSTAR[MXFILTINDX]; // <lam> per band for calib stars 
  double n_CALSTAR[MXFILTINDX];      // index of refrac per band, calib stars

  double SNRMIN; // min SNR to include in RA/DEC avg

  // info per event

  // GENLC struct in snlc_sim.h stores RA/DEC per obs;
  // here keep track of avg coord per band and over all obs.
  COORD_AVG_DEF COORD_RA;
  COORD_AVG_DEF COORD_DEC;
  COORD_AVG_DEF COORD_SIM_RA;
  COORD_AVG_DEF COORD_SIM_DEC;

} ATMOS_INFO ;

// ------- function prototypes ----------

void  GEN_ATMOSPHERE_DRIVER(void) ;
void  INIT_ATMOSPHERE(void);
void  INIT_EVENT_ATMOSPHERE(void) ;

void reset_COORD_AVG(COORD_AVG_DEF *COORD);
void sum_COORD_AVG(COORD_AVG_DEF *COORD,
                   double RA_OBS, double WGT, int IFILT_OBS);

void  gen_airmass(int ep);
void  genSmear_coords(int ep);
void  gen_dcr_coordShift(int ep);

void  gen_dcr_magShift(int ep);
double gen_wave_sed_wgted(int ep);

double compute_DCR_angle(double LAM, double tan_ZENITH, 
			 int IFILT_OBS, int DUMPFLAG);
double compute_index_refrac_atmos(double LAM, int DUMPFLAG) ;
void   test_compute_dcr(void);

