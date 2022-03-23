/*******************************************
   Created March 2022
    [separated from simlib_tools.c]


********************************************/

// define index params for INFO_HEAD array

#define NPAR_HEAD  7
#define IPAR_RA       0
#define IPAR_DEC      1
#define IPAR_MWEBV    2 
#define IPAR_PIXSIZE  3
#define IPAR_Z        4
#define IPAR_PKMJD    5 
#define IPAR_NEA_UNIT 6  // Dec 2021, RK

// define index parameters for.INFO_OBS array
#define NPAR_OBS     10
#define IPAR_MJD      0
#define IPAR_CCDGAIN  1
#define IPAR_CCDNOISE 2
#define IPAR_SKYSIG   3
#define IPAR_PSF0     4
#define IPAR_ZPT0     7
#define IPAR_MAG      9

// - - - - - -
void simlib_open_write(char *filename, char *surveyname, char *filters, 
		       char *telescope, char *comment, char *headFile );

void simlib_add_header(int optflag, int IDLIB, int NOBS, 
		       char *FIELD, float *INFO );

void simlib_add_mjd(int opt, double *INFO, char * STRINGID, char *FILTNAME);

void simlib_close_write(void);


int  CHECK_SIMLIB_VAL(char *varname, float value, float varmin, float varmax);

void parse_SIMLIB_IDplusNEXPOSE(char *inString, int *IDEXPT, int *NEXPOSE) ;

void PRINT_SIMLIB_ERROR(char *msgerr );
