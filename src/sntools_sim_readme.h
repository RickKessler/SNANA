/***********************************
   Created Dec 22 2021

**********************************/

#define DOCANA_OVERVIEW         "OVERVIEW"
#define DOCANA_INPUT_KEYS       "INPUT_KEYS"
#define DOCANA_INPUT_NOTES      "INPUT_NOTES"
#define DOCANA_OUTPUT_SUMMARY   "OUTPUT_SUMMARY"

char  ORIG_FILE_README[MXPATHLEN]; // temp space to restore original filenames

// define KEY+ARG strings for keys that have variable number of 
// arguments, or that can have duplicate keys.
typedef struct {
  int NKEY;
  bool MALLOC1;
  char **KEY_LIST, **ARG_LIST;
} README_KEYPLUSARGS_DEF ;

README_KEYPLUSARGS_DEF README_KEYS_COSMO ;
README_KEYPLUSARGS_DEF README_KEYS_GENMODEL ;
README_KEYPLUSARGS_DEF README_KEYS_CID ;
README_KEYPLUSARGS_DEF README_KEYS_SIMLIB ; 
README_KEYPLUSARGS_DEF README_KEYS_HOSTLIB ; 
README_KEYPLUSARGS_DEF README_KEYS_RATEMODEL ; 
README_KEYPLUSARGS_DEF README_KEYS_LENS ; 
README_KEYPLUSARGS_DEF README_KEYS_SKY ;  // solid angl, RA, DEC ...
README_KEYPLUSARGS_DEF README_KEYS_MWEBV ; 
README_KEYPLUSARGS_DEF README_KEYS_NON1ASED ;
README_KEYPLUSARGS_DEF README_KEYS_SIMSED ;
README_KEYPLUSARGS_DEF README_KEYS_LCLIB ;
README_KEYPLUSARGS_DEF README_KEYS_FILTER ;
README_KEYPLUSARGS_DEF README_KEYS_FLUXERRMODEL ;
README_KEYPLUSARGS_DEF README_KEYS_GENMAG_OFF ;
README_KEYPLUSARGS_DEF README_KEYS_GENMAG_SMEAR ;
README_KEYPLUSARGS_DEF README_KEYS_TAKE_SPECTRUM ;
README_KEYPLUSARGS_DEF README_KEYS_RANSYSTPAR ;
README_KEYPLUSARGS_DEF README_KEYS_ZVARIATION ;
README_KEYPLUSARGS_DEF README_KEYS_GRIDGEN ; 
README_KEYPLUSARGS_DEF README_KEYS_CUTWIN ;
README_KEYPLUSARGS_DEF README_KEYS_COVMAT_SCATTER ;
README_KEYPLUSARGS_DEF README_KEYS_SIMGEN_DUMP ;


// -------- function prototypes -----------

void  README_DOCANA_DRIVER(int iflag_readme);
void  README_DOCANA_OVERVIEW(int *iline);
void  README_DOCANA_INPUT_KEYS(int *iline);
void  README_DOCANA_FILTERS(int *iline);
void  README_DOCANA_INPUT_NOTES(int *iline);
void  README_DOCANA_OUTPUT_SUMMARY(int *iline);
void  README_DOCANA_SED_TRUE(int *iline);
void  README_DOCANA_GENTYPE_MAP(int *iline);

void  readme_docana_output(int *iline, char *pad); 
void  readme_docana_genmodel(int *iline, char *pad);
void  readme_docana_cid(int *iline, char *pad);
void  readme_docana_redshift(int *iline, char *pad);
void  readme_docana_epoch(int *iline, char *pad); // MJD, Trest ...
void  readme_docana_instr(int *iline, char *pad);
void  readme_docana_hostlib(int *iline, char *pad);
void  readme_docana_hostmatch(int *iline, char *pad);
void  readme_docana_modelPar(int *iline, char *pad) ;
void  readme_docana_rate(int *iline, char *pad) ;
void  readme_docana_cutwin(int *iline, char *pad) ;
void  readme_docana_searcheff(int *iline, char *pad);
void  readme_docana_mwebv(int *iline, char *pad);
void  readme_docana_misc(int *iline, char *pad);

void onoff_readme_docana(int FLAG, char *onoff);

void  readme_docana_comment(int *iline, char *comment); // write # <comment> 

void get_TIME_START_readme_docana(char *TIME_START);

void VERSION_INFO_load(int *iline, char *pad, char *keyName, char *comment,
                       int lenkey, bool isint,
                       int nval, double *val, double valmin, double valmax,
                       double val_noprint ) ;

void README_KEYPLUSARGS_load(int MXKEY, int NWD, char **WORDS, int keySource,
			     README_KEYPLUSARGS_DEF *README_KEYS, 
			     char *callFun);
void README_KEYPLUSARGS_init(README_KEYPLUSARGS_DEF *README_KEYS) ;

void  readme_docana_load_list(int *iline, char *pad,
			      README_KEYPLUSARGS_DEF *README_KEYS);

void readme_docana_load_asymGauss(int *iline, char *pad, 
				  GENGAUSS_ASYM_DEF *GENGAUASS);

void readme_docana_load_expHalfGauss(int *iline, char *pad, 
				     GEN_EXP_HALFGAUSS_DEF *EXP_HALFGAUASS);

// === END ===
