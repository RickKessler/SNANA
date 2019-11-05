
// Created Oct 2019 by R.Kessler

#define MXPRIMARY_KCOR 6

#define MXFILT_KCOR       MXFILTINDX
#define MXFILT_REST_KCOR  20
#define MXFILT_OBS_KCOR   62

#define MXTABLE_KCOR   100   // max numver of kcor tables

#define MXLAMBIN_FILT  4000  // max wave bins for filter trans
#define MXLAMBIN_PRIM  4000  // max wave bins for primary

#define MXTBIN_KCOR   150
#define MXZBIN_KCOR   100
#define MXAVBIN_KCOR  100
#define MXCBIN_AVWARP 100  // max color index for AVWARP table

#define MASK_FRAME_REST 1
#define MASK_FRAME_OBS  2

#define OPT_FRAME_REST 0
#define OPT_FRAME_OBS  1

#define KDIM_T        0
#define KDIM_z        1  
#define KDIM_AV       2
#define KDIM_IFILTr   3
#define KDIM_IFILTo   4
#define NKDIM_KCOR    5
#define N4DIM_KCOR    4

#define IDMAP_KCOR_TABLE   15
#define IDMAP_KCOR_AVWARP  16
#define IDMAP_KCOR_LCMAG   17
#define IDMAP_KCOR_MWXT    18

int KCOR_VERBOSE_FLAG;
int IFILTDEF_BESS_BX;

typedef struct {
  char VARNAME[40];
  int NBIN;
  double BINSIZE, RANGE[2] ;  // binsize and min/max range
  double *GRIDVAL ;           // value at each bin
} KCOR_BININFO_DEF ;

typedef struct {
  char NAME[40];
  int IDMAP ;
  int NBINTOT ;
  int NBIN[NKDIM_KCOR];
  int NDIM ;
} KCOR_MAPINFO_DEF ;

typedef struct {
  int  NFILTDEF ;
  int  IFILTDEF[MXFILT_KCOR];    // vs. sparse index
  int  IFILTDEF_INV[MXFILT_KCOR]; // vs. absolute index
  char FILTERSTRING[MXFILT_KCOR];
} KCOR_FILTERMAP_DEF ;


struct KCOR_INFO {

  char FILENAME[MXPATHLEN] ;
  fitsfile *FP ;

  char FILTERS_SURVEY[MXFILT_KCOR]; // filter list read from SIMLIB file
  int  NFILTDEF_SURVEY ;

  // header info
  int   VERSION, NFLAMSHIFT;

  int   NPRIMARY;
  char  *PRIMARY_NAME[MXPRIMARY_KCOR] ;
  int    PRIMARY_INDX[MXFILT_KCOR];
  double PRIMARY_MAG[MXFILT_KCOR];    // native mag
  double PRIMARY_ZPOFF[MXFILT_KCOR];  // mag(native) - mag(synth)
  double SNPHOT_ZPOFF[MXFILT_KCOR];   // apply to photometry (from ZPOFF.DAT)

  int   NFILTDEF;
  char *FILTER_NAME[MXFILT_KCOR] ;
  int   IFILTDEF[MXFILT_KCOR] ; 
  int   MASK_FRAME_FILTER[MXFILT_KCOR]; // bits 0,1 --> rest, obs
  int   MASK_EXIST_BXFILT; // logical for BX existing for rest/obs frame

  KCOR_FILTERMAP_DEF FILTERMAP_REST ;
  KCOR_FILTERMAP_DEF FILTERMAP_OBS  ;

  int   NKCOR, NKCOR_STORE;
  char *KCOR_STRING[MXTABLE_KCOR] ;
  bool  EXIST_KCOR[MXFILT_KCOR][MXFILT_KCOR]; // for each [rest][obs] combo
  int   IFILTMAP_KCOR[2][MXTABLE_KCOR];
  int   k_index[MXTABLE_KCOR];

  bool  ISLAMSHIFT[MXFILT_KCOR] ;
  double RVMW;
  int    OPT_MWCOLORLAW ;
  
  KCOR_BININFO_DEF BININFO_LAM;
  KCOR_BININFO_DEF BININFO_T;
  KCOR_BININFO_DEF BININFO_z;
  KCOR_BININFO_DEF BININFO_AV ;

  KCOR_MAPINFO_DEF MAPINFO_KCOR ;
  KCOR_MAPINFO_DEF MAPINFO_AVWARP ;
  KCOR_MAPINFO_DEF MAPINFO_LCMAG ;
  KCOR_MAPINFO_DEF MAPINFO_MWXT ;

  double zRANGE_LOOKUP[2] ;

  float *FLUX_SNSED; // use float to save memory

  char SPECTROGRAPH_INSTRUMENT[60];
  char SPECTROGRAPH_FILTERLIST[MXFILT_KCOR];
  int  NFILTDEF_SPECTROGRAPH ;
  int  IFILTDEF_SPECTROGRAPH[MXFILT_KCOR] ;

  // misc init info
  int NCALL_READ ;
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
int ISBXFILT_KCOR(char *cfilt);
void addFilter_kcor(int ifiltdef, KCOR_FILTERMAP_DEF *MAP);
void init_kcor_indices(void);
void get_MAPINFO_KCOR(char *what, KCOR_MAPINFO_DEF *MAPINFO); 

// end

