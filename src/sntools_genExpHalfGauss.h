
#define ILIST_GETRAN_GENEXP 1  // use this random list

// March 20 2020: Generic struct for exponential and half gaussian. 
typedef struct {
  bool   USE;          // T => values are set 
  char   NAME[80];     // name of variable
  double EXP_TAU ;     // exponential component: exp(-x/EXP_TAU)
  double PEAK, SIGMA ; // peak & sigma of half gaussian component
  double RATIO ;       // Gauss(0)/Expon(0)
  double RANGE[2] ;    // generate random value in this RANGE

  double PROB_EXPON_REWGT, SQRT_PROB_EXPON_REWGT ;

  int INDEX; // Generic index for internal use (not part of function)
  int KEYSOURCE ;  // 1=FILE, 2=command line ; used to set priority

} GEN_EXP_HALFGAUSS_DEF ;


// ========= FUNCTION PROTOTYPES ==================

void init_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS, 
			    double VAL );

void setUseFlag_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS, 
				  char *name );

void set_GEN_EXPON(double tau, double *range, 
		   GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS);

void set_GEN_EXPON_REWGT(double expon_rewgt,
			 GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS);

double getRan_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS);
double funVal_GEN_EXP_HALFGAUSS(double x, GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS);

void copy_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *inp_EXP_HALFGAUSS, 
			    GEN_EXP_HALFGAUSS_DEF *out_EXP_HALFGAUSS);

void dump_GEN_EXP_HALFGAUSS(GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS);

int parse_input_EXP_HALFGAUSS(char *VARNAME, char **WORDS, int keySource,
			      GEN_EXP_HALFGAUSS_DEF *gen_EXP_HALFGAUSS );


