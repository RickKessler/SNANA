/*******************************
  Created Feb 12 2018

  legacy fluxErr models from  
    sntools_fluxErrModels.h  sntools_imageNoiseModels.h

  
********************************/

// ==================================================
// ============ sntools_fluxErrModels.h =============
// ==================================================

struct {

  char SIMLIB_FILENAME[400];
  int  OPTFLAG;           // = INPUTS.FUDGEOPT_FLUXERR ;

  char   BAND[4], FIELD[40];
  int    IBAND;
  double MJD, ZP, SKYSIG, PSF ;

} FLUXERRMODEL ;


// function prototypes

void init_fluxErrModel_legacy(char *SIMLIB_FILENAME, int OPTFLAG);

double scale_fluxErrModel_legacy(char *BAND, char *FIELD, double MJD, 
				 double ZP, double SKYSIG, double PSF);
double scale_fluxErr_DES_DIFFIMG(void);
double scale_fluxErr_DES_SMP(void);




// ==================================================
// ======== sntools_imageNoiseModels.h ==============
// ==================================================


#define MXMAP_NOISEMODEL_HOST 200  // max number of maps (filterx field)
#define MXPAR_NOISEMODEL       10  // max number of params per HOSTMAG
#define MXBIN_HOSTMAG          20

char NOISEMODEL_FILE[MXPATHLEN]; // name of file
char NOISEMODEL_NAME[MXPATHLEN]; // name of noise model
int  NPAR_NOISEMODEL ;
int  NMAP_NOISEMODEL_HOST ;

struct {
  char BANDLIST[MXPATHLEN];  // e.g., g+r or r or i
  char FIELDLIST[MXPATHLEN]; // e.g., C3+X3

  int    LIBID ;
  int    NBIN_HOSTMAG ;
  double HOSTMAG[MXBIN_HOSTMAG] ;
  double PARAM[MXBIN_HOSTMAG][MXPAR_NOISEMODEL] ;
} NOISEMODEL_HOST[MXMAP_NOISEMODEL_HOST];


// ========== FUNCTIONS =============

void   INIT_NOISEMODEL_HOST_LEGACY(char *HOSTNOISE_FILE) ;
void   init_noisemodel_host_legacy__(char *HOSTNOISE_FILE);

void GEN_NOISEMODEL_HOST_LEGACY(char *BAND,  char *FIELD, int GALID, 
				double GALMAG, double SBMAG, 
				double SNSEP, double *noisePar);
void gen_noisemodel_host_legacy__(char *BAND, char *FIELD, int *GALID, 
				  double *GALMAG, double *SBMAG, 
				  double *SNSEP, double *noisePar);


