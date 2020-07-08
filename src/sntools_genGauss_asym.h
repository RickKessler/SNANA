// Created Sep 2 2016 [yanked out of sntools.h]
// Tools to generate randoms from asymmetric Gaussian distributions

// ========= GLOBAL DECLARATIONS =================

#define MXGENGAUSS 100
int NFUN_GENGAUSS_ASYM ; // used in conjunction with FUNINDEX       

GENGAUSS_ASYM_DEF  GENGAUSS_ASYM_LIST[MXGENGAUSS] ; 


// ========= FUNCTION PROTOTYPES ==================

void   init_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss, double VAL );
double exec_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss);
void   copy_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss1,
                          GENGAUSS_ASYM_DEF *genGauss2) ;
void   dump_GENGAUSS_ASYM(GENGAUSS_ASYM_DEF *genGauss);

void   set_GENGAUSS_ASYM(double peak, double *sigma, double *range,
			 GENGAUSS_ASYM_DEF *genGauss);

void prepIndex_GENGAUSS(char *varName, GENGAUSS_ASYM_DEF *genGauss );

// function prototypes for skewNormal.c (moved here Jan 11 2017)
void  init_skewNormal(int seed); // one-time init
double skewNormalRan(int seed, double loc, double scale, double skew) ;

// == END ==
