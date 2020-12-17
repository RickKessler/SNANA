/*******************************************************
     Created Oct 2020 by R.Kessler 

     cosmology theory functions: H(z), MU(z) ...

********************************************************/

#define ICOSPAR_HzFUN_H0  0
#define ICOSPAR_HzFUN_OM  1
#define ICOSPAR_HzFUN_OL  2
#define ICOSPAR_HzFUN_w0  3
#define ICOSPAR_HzFUN_wa  4
#define NCOSPAR_HzFUN     5

#define MXMAP_HzFUN 5000  

typedef struct {
  double COSPAR_LIST[NCOSPAR_HzFUN];
  
  // optional 2-column map to define theory H(z)
  bool   USE_MAP ;
  char   *FILENAME ;
  int    Nzbin_MAP;
  double *zCMB_MAP, *HzFUN_MAP ;

} HzFUN_INFO_DEF ;


// ========= function prototypes =========

void init_HzFUN_INFO(int VBOSE, double *cosPar, char *fileName, 
		     HzFUN_INFO_DEF *HzFUN_INFO); 
void write_HzFUN_FILE(HzFUN_INFO_DEF *HzFUN_INFO);

double SFR_integral(double z, HzFUN_INFO_DEF *HzFUN_INFO);
double SFRfun_BG03(double z,  double H0 ) ;
double SFRfun_MD14(double z,  double *params);

double dVdz_integral(int OPT, double zmax, HzFUN_INFO_DEF *HzFUN_INFO);
		     
double dvdz_integral__(int *OPT, double *zmax, double *COSPAR);

double dVdz ( double z, HzFUN_INFO_DEF *HzFUN_INFO);
double Hzinv_integral(double zmin, double zmax, HzFUN_INFO_DEF *HzFUN_INFO); 

double Hainv_integral(double amin, double amax, HzFUN_INFO_DEF *HzFUN_INFO); 

double Hzfun ( double z, HzFUN_INFO_DEF *HzFUN_INFO); 
double Hzfun_wCDM ( double z, HzFUN_INFO_DEF *HzFUN_INFO); 
double Hzfun_interp ( double z, HzFUN_INFO_DEF *HzFUN_INFO); 
double dLmag ( double zCMB, double zHEL, HzFUN_INFO_DEF *HzFUN_INFO); 

double dlmag_fortc__(double *zCMB, double *zHEL, double *H0,
                     double *OM, double *OL, double *w0, double *wa);

double zcmb_dLmag_invert(double MU, HzFUN_INFO_DEF *HzFUN_INFO); 

double zhelio_zcmb_translator(double z_input, double RA, double DECL, 
			      char *coordSys, int OPT ) ;
double zhelio_zcmb_translator__(double *z_input, double *RA, double *DECL, 
				char *coordSys, int *OPT ) ;

// ============== END OF FILE =============
