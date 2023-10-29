// ==============================
//  genmag_NON1ASED.h
//  
//  Mar 14 2016: add prep_NON1ASED & pick_NON1ASED
// ==============================

void  prep_NON1ASED(INPUTS_NON1ASED_DEF *INP_NON1ASED, 
		    GENLC_NON1ASED_DEF *GEN_NON1ASED);

bool getName_SED_FILE_NON1ASED(char *PATH, char *inpName, char *outName);

void read_NON1A_LIST(INPUTS_NON1ASED_DEF *INP_NON1ASED );
int  count_NON1A_LIST(char *PATHMODEL);
int  ISKEY_NON1A_LIST(char *string);

void sort_NON1ASED(INPUTS_NON1ASED_DEF *INP_NON1ASED);

void copy_NON1ASED(int i1, int i2,
                   INPUTS_NON1ASED_DEF *NON1ASED1,
                   INPUTS_NON1ASED_DEF *NON1ASED2 ) ;
			
void init_genmag_NON1ASED(int isparse, INPUTS_NON1ASED_DEF *INP_NON1ASED,
			  int OPTMASK );

void genmag_NON1ASED ( int index, int ifilt_obs, 
		       double mwebv,  double z, double x0,
		       double RVhost, double AVhost,
		       int    Nobs,   double *ptr_epoch,
		       double *ptr_genmag, double *ptr_generr );

void shift_NON1A_DAY(void);
		   
// global var
#define  ISED_NON1A  1
#define  MODELNAME_NON1ASED  "NON1ASED" 
