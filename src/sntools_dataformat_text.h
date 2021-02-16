
#define OPTMASK_TEXT_HEAD 2  // matches snana.car
#define OPTMASK_TEXT_OBS  4
#define OPTMASK_TEXT_SPEC 8

#define MSKOPT_PARSE_TEXT_FILE  MSKOPT_PARSE_WORDS_FILE + MSKOPT_PARSE_WORDS_IGNORECOMMENT

#define MXVAROBS_TEXT 20
struct {
  int MJD, BAND, FIELD, FLUXCAL, FLUXCALERR ;
  int ZPFLUX, PSF, SKYSIG, SKYSIG_T, SKYSIG_GAIN ;
  int GAIN, PHOTFLAG, PHOTPROB, SIM_MAGOBS ;
} IVAROBS_TEXT ;


struct {
  int  NVERSION ;
  char DATA_PATH[MXPATHLEN];   
  char PHOT_VERSION[MXPATHLEN];   
  char LIST_FILE[MXPATHLEN] ;
  char README_FILE[MXPATHLEN];

  int  NFILE;
  char **DATA_FILE_LIST ;

} TEXT_VERSION_INFO ;


struct {
  int IPTR_READ ;  // pointer to current word in file
  int NWD_TOT ;   // total number of words in file
  int NVAROBS ;   // number of variables following each OBS key
  char VARNAME_LIST[MXVAROBS_TEXT][32] ;
  char STRING_LIST[MXVAROBS_TEXT][20] ;

  int NOBS_READ ;
} TEXT_FILE_INFO ;


void WR_DATAFILE_TEXT(char *OUTFILE);
void wr_datafile_text__(char *OUTFILE);

void wr_dataformat_text_HEADER(FILE *fp ) ;
void wr_dataformat_text_HOSTGAL(FILE *fp) ;
void wr_dataformat_text_SIMPAR(FILE *fp ) ;
void wr_dataformat_text_SNPHOT(FILE *fp ) ;
void wr_dataformat_text_SNSPEC(FILE *fp ) ;

void RD_TEXT_INIT(void); // one-time init
void rd_text_init__(void);

int RD_TEXT_PREP(int MSKOPT, char *PATH, char *VERSION);
int rd_text_prep__(int *MSKOPT, char *PATH, char *VERSION);

int  rd_text_list(void);
void rd_text_malloc_list(int OPT, int NFILE) ;
void rd_text_global(void);
void rd_text_varlist(int *iwd_file);

void RD_TEXT_EVENT(int OPTMASK, int ifile);
void rd_text_event__(int *OPTMASK, int *ifile);
bool parse_TEXT_HEAD(int *iwd);
bool parse_TEXT_OBS(int *iwd);
bool parse_TEXT_SPEC(int *iwd);

void check_plusminus_TEXT(int *iwd_file, float *PTR_ERR);
void parse_plusminus_TEXT(char *word, char *key, int *iwd_file, 
			  float *PTR_VAL, float *PTR_ERR) ;
