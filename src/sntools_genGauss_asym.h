// Created Sep 2 2016 [yanked out of sntools.h]
// Tools to generate randoms from asymmetric Gaussian distributions

// ========= GLOBAL DECLARATIONS =================

#define MXGENGAUSS 100
int NFUN_GENGAUSS_ASYM ; // used in conjunction with FUNINDEX       


typedef struct  {
  bool   USE;       // T -> values are set (Jun 11 2020)  
  char   NAME[80];  // name of variable  
  double PEAK ;     // peak prob
  double PEAKRANGE[2]; // Range where PEAKPROB = 1 
  double SIGMA[2] ; // asymmetric Gaussian sigmas
  double SKEW[2] ;  // hack-skew; TrueSigma = SIGMA + SKEW*|x-PEAK|
  double RANGE[2] ; // allows truncation  
  int    NGRID ;      // if non-zero, snap to grid 

  // Mar 29 2017; add 2nd peak params (for low-z x1) 
  double PROB2;     // prob of generating 2nd peak (default=0.0) 
  double PEAK2;     // location of 2nd peak
  double SIGMA2[2]; // asym Gaussian sigmas of 2nd peak 
  int  FUNINDEX;    // = NFUN_GENGUASS_ASYM = unique index
  double RMS;  // RMS of asym Gaussian 

  int INDEX; // Generic index for internal use (not part of function)

} GENGAUSS_ASYM_DEF ;

GENGAUSS_ASYM_DEF  GENGAUSS_ASYM_LIST[MXGENGAUSS] ; 


// ========= FUNCTION PROTOTYPES ==================

void   init_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss, double VAL );
double getRan_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss);
double funVal_GENGAUSS_ASYM(double x, GENGAUSS_ASYM_DEF *genGauss);
void   copy_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss1,
                          GENGAUSS_ASYM_DEF *genGauss2) ;
void   dump_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss);

void   set_GENGAUSS_ASYM(double peak, double *sigma, double *range,
			 GENGAUSS_ASYM_DEF *genGauss);

void prepIndex_GENGAUSS(char *varName, GENGAUSS_ASYM_DEF *genGauss );

// function prototypes for skewNormal.c (moved here Jan 11 2017)
void  init_skewNormal(int seed); // one-time init
double skewNormalRan(int seed, double loc, double scale, double skew) ;


void checkVal_GENGAUSS(char *varName, double *val, char *fromFun );
int  parse_input_GENGAUSS(char *VARNAME, char **WORDS, int keySource,
                            GENGAUSS_ASYM_DEF *genGauss );

// == END ==
