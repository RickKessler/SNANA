

void  GEN_ATMOSPHERE_DRIVER(void) ;
void  gen_airmass(int ep);
void  genSmear_coords(int ep);
void  gen_dcr_coordShift(int ep);

void  gen_dcr_magShift(int ep);
double gen_wave_sed_wgted(int ep);

double compute_DCR_angle(double LAM, double TAN_ZENITH, int DUMPFLAG);
double compute_index_refrac_atmos(double LAM, int DUMPFLAG) ;
void   test_compute_dcr(void);
