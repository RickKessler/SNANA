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

typedef struct {
  double COSPAR_LIST[NCOSPAR_HzFUN];
  
  // optional 2-column map to define theory H(z)
  bool   USE_MAP ;
  int    Nzbin_MAP;
  double *zCMB_MAP, *HzFUN_MAP ;

} HzFUN_INFO_DEF ;


// ========= function prototypes =========

void init_HzFUN_INFO(double *cosPar, char *fileName, HzFUN_INFO_DEF *HzFUN);

double SFR_integral(double H0, double OM, double OL, double W, double Z );
double SFRfun_BG03(double H0, double z) ;
double SFRfun_MD14(double z, double *params);

double dVdz_integral(double H0, double OM, double OL, double W, 
		     double Zmax, int wgtopt ) ;

double dvdz_integral__(double *H0, double *OM, double *OL, double *W, 
		       double *Zmax, int *wgtopt ) ;

double dVdz ( double H0, double OM, double OL, double W, double Z ) ;
double Hzinv_integral(double H0, double OM, double OL, double W, 
		      double Zmin, double Zmax ) ;

double Hainv_integral(double H0, double OM, double OL, double W, 
		      double amin, double amax ) ;

double Hzfun ( double H0, double OM, double OL, double W, double Z ) ;
double dLmag ( double H0, double OM, double OL, double W, 
	       double zCMB, double zHEL ) ;

double zcmb_dLmag_invert(double H0, double OM, double OL, double W, double MU);

double zhelio_zcmb_translator(double z_input, double RA, double DECL, 
			      char *coordSys, int OPT ) ;
double zhelio_zcmb_translator__(double *z_input, double *RA, double *DECL, 
				char *coordSys, int *OPT ) ;

// ============== END OF FILE =============
