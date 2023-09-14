/******************************************

 Created Sep 2023 by R.Kessler
 
 *****************************************/

// define model parameters for late-time SNIa mag-extrapolation;
// based on data from K.Mandel.
#define MXLAMBIN_EXTRAP_LATETIME_Ia 20
#define MXPAR_EXTRAP_LATETIME_Ia    10 // includes computer params
#define NPAR_EXTRAP_LATETIME_Ia  4     // used in F(t) function
#define IPAR_EXTRAP_LAM_Ia       0
#define IPAR_EXTRAP_TAU1_Ia      1
#define IPAR_EXTRAP_TAU2_Ia      2
#define IPAR_EXTRAP_EXPRATIO_Ia  3
#define IPAR_EXTRAP_MAGSLOPE1_Ia 4  // computed magslope for tau1
#define IPAR_EXTRAP_MAGSLOPE2_Ia 5  // computed magslope for tau2
#define IPAR_EXTRAP_DAYPIVOT_Ia  6  // day when each flux is the same

struct {
  char   FILENAME[MXPATHLEN];
  double DAYMIN ; 
  int NLAMBIN;  // number of lambda bins to defie F(t)

  // Define parameters vs. ilambin.
  // index order allows easy interpolation vs. lambda
  //                 ipar                    ilambin
  double PARLIST[MXPAR_EXTRAP_LATETIME_Ia][MXLAMBIN_EXTRAP_LATETIME_Ia] ;

} INPUT_EXTRAP_LATETIME_Ia ;


// ========= function prototypes ============

// SNIa-specific
void   init_extrap_latetime_Ia(void);
double genmag_extrap_latetime_Ia(double mag_daymin, double day, double lam);
double FLUXFUN_EXTRAP_LATETIME_Ia(double t, double tau1, double tau2, double ratio);

// generic for any model
double modelflux_extrap(double T, double Tref,
			double fluxref, double fluxslope, int LDMP);

