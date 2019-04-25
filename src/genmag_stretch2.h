// genmag_stretch2.h

#define MXEPOCH_STRTEMPL 200
#define MXFILT_STRTEMPL 6

char   STRETCH2_MODELPATH[MXPATHLEN] ;
char   STRTEMPL_FILE[MXPATHLEN] ;
double STRTEMPL_MAGERR[MXFILT_STRTEMPL] ;
char   STRTEMPL_FILTERS[MXFILT_STRTEMPL] ;
char   STRTEMPL_SNTYPE[20];

// define array of rest-frame days to store template
double STRTEMPL_TREST[MXEPOCH_STRTEMPL] ;
double STRTEMPL_MAGREST[MXEPOCH_STRTEMPL][MXFILT_STRTEMPL];

double TMINDEF, TMAXDEF;
int    NEPOCH_STRTEMPL;
int    NFILT_STRTEMPL;
int    IFILTDEF_STRTEMPL[MXFILT_STRTEMPL];
char   PATH_STRETCH[200];


// legacy values for backward compatibility
#define NFILT_STRTEMPL_LEGACY 5   
char    STRTEMPL_FILTERS_LEGACY[MXFILT_STRTEMPL] ;
double  STRTEMPL_MAGERR_LEGACY[MXFILT_STRTEMPL] ;
int     LEGACYFLAG ; 

// prototypes

int  init_genmag_stretch2(char *version, char *cfilt_rest);
void init_stretch2_LEGACY(void);

int genmag_stretch2 (
		     double str1
		    ,double str2
		    ,int ifilt_rest
		    ,int nepoch
		    ,double *Trest
		    ,double *genmag 
		    ,double *genmag_err
		    );


void JDATES2 ( double Trest, int *j1, int *j2 );

// END
