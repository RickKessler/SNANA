
#define MXPRIMARY_KCOR 6
#define MXFILT_KCOR    MXFILTINDX
#define MXTABLE_KCOR   100   // max numver of kcor tables

#define MXLAMBIN_FILT  4000  // max wave bins for filter trans
#define MXLAMBIN_PRIM  4000  // max wave bins for primary

#define MXTBIN_KCOR  150
#define MXZBIN_KCOR  100
#define MXAVBIN_KCOR 100

int KCOR_VERBOSE_FLAG;

typedef struct {
  char VARNAME[40];
  int NBIN;
  double BINSIZE, RANGE[2] ;  // binsize and min/max range
  double *GRIDVAL ;           // value at each bin
} KCOR_BININFO_DEF ;

struct KCOR_INFO {

  char FILENAME[MXPATHLEN] ;
  fitsfile *FP ;

  char FILTERS_SURVEY[MXFILT_KCOR]; // filter list read from SIMLIB file

  // header info
  int   VERSION, NPRIMARY, NFILTDEF, NFLAMSHIFT, NKCOR ;

  char  *PRIMARY_NAME[MXPRIMARY_KCOR] ;
  int    PRIMARY_INDX[MXFILT_KCOR];
  double PRIMARY_MAG[MXFILT_KCOR];    // native mag
  double PRIMARY_ZPOFF[MXFILT_KCOR];  // mag(native) - mag(synth)
  double SNPHOT_ZPOFF[MXFILT_KCOR];   // apply to photometry (from ZPOFF.DAT)

  char *FILTER_NAME[MXFILT_KCOR] ;
  int   IFILTDEF[MXFILT_KCOR] ; 
  int   MASK_FRAME_FILTER[MXFILT_KCOR]; // bits 0,1 --> rest, obs

  char *KCOR_STRING[MXTABLE_KCOR] ;
  bool  EXIST_KCOR[MXFILT_KCOR][MXFILT_KCOR]; // for each [rest][obs] combo
  bool  ISLAMSHIFT[MXFILT_KCOR] ;
  double RVMW;
  int    OPT_MWCOLORLAW ;
  
  KCOR_BININFO_DEF BININFO_LAM;
  KCOR_BININFO_DEF BININFO_T;
  KCOR_BININFO_DEF BININFO_z;
  KCOR_BININFO_DEF BININFO_AV ;

  double zRANGE_LOOKUP[2] ;

  float *FLUX_SNSED; // use float to save memory

  char SPECTROGRAPH_INSTRUMENT[60];
  char SPECTROGRAPH_FILTERLIST[MXFILT_KCOR];
  int  NFILTDEF_SPECTROGRAPH ;
  int  IFILTDEF_SPECTROGRAPH[MXFILT_KCOR] ;

  // misc init info
  int NCALL_READ ;
  int NKCOR_STORE ;
  bool STANDALONE;

} KCOR_INFO ;


// declare functions

void READ_KCOR_DRIVER(char *kcorFile, char *FILTERS_SURVEY );
void read_kcor_init(void);
void read_kcor_open(void);
void read_kcor_head(void);
void read_kcor_zpoff(void);
void read_kcor_snsed(void);
void read_kcor_tables(void);
void read_kcor_mags(void);
void read_kcor_filters(void);
void read_kcor_primarysed(void);

void read_kcor_binInfo(char *VARNAME, char *VARSYM, int MXBIN,
		       KCOR_BININFO_DEF *BININFO) ;

void parse_KCOR_STRING(char *STRING, 
		       char *strKcor, char *cfilt_rest, char *cfilt_obs);
