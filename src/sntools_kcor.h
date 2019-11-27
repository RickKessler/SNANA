
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

  // filled by addFilter_kcor
  int  NFILTDEF;
  int  IFILTDEF[MXFILT_KCOR];     // vs. sparse index
  int  IFILTDEF_INV[MXFILT_KCOR]; // vs. absolute index
  int  NDEFINE[MXFILT_KCOR];      // Number of times each IFILTDEF is defined
  int  NFILT_DUPLICATE;           // number of bands with duplicates
  char FILTERSTRING[MXFILT_KCOR]; // list of single-char bands
  char *FILTER_NAME[MXFILT_KCOR]; // full name of each filter, vs. sparse indx

  double PRIMARY_MAG[MXFILT_KCOR];
  double PRIMARY_ZPOFF[MXFILT_KCOR];

  // below is filled by loadFilterTrans_kcor
  int    NBIN_LAM[MXFILT_KCOR];
  double TRANS_MAX[MXFILT_KCOR];
  double LAMMEAN[MXFILT_KCOR], LAMRMS[MXFILT_KCOR];
  float *LAM[MXFILT_KCOR], *TRANS[MXFILT_KCOR];

} KCOR_FILTERMAP_DEF ;


struct KCOR_INFO {

  // info passed to driver
  char FILENAME[MXPATHLEN] ;
  fitsfile *FP ;

  char FILTERS_SURVEY[MXFILT_KCOR]; // filter list read from SIMLIB file
  int  NFILTDEF_SURVEY ;
  
  double MAGREST_SHIFT_PRIMARY[MXFILT_KCOR];
  double MAGOBS_SHIFT_PRIMARY[MXFILT_KCOR];

  // - - - - info read from FITS file - - - - - - 

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
  bool  IS_SURVEY_FILTER[MXFILT_KCOR];

  KCOR_FILTERMAP_DEF FILTERMAP_REST ;
  KCOR_FILTERMAP_DEF FILTERMAP_OBS  ;

  int   NKCOR, NKCOR_STORE;
  char *STRING_KCORLINE[MXTABLE_KCOR] ; // full FITS line with KCOR info
  char *STRING_KCORSYM[MXTABLE_KCOR];   // just K_xy part
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

  float *KCORTABLE1D_F;     // use float to save memory
  float *LCMAG_TABLE1D_F ;
  float *AVWARP_TABLE1D_F ;
  float *MWXT_TABLE1D_F ;
  float *FLUX_SNSED_F;      // use float to save memory

  char SPECTROGRAPH_INSTRUMENT[60];
  char SPECTROGRAPH_FILTERLIST[MXFILT_KCOR];
  int  NFILTDEF_SPECTROGRAPH ;
  int  IFILTDEF_SPECTROGRAPH[MXFILT_KCOR] ;

  // misc init info
  int NCALL_READ ;
  bool STANDALONE;


} KCOR_INFO ;


// declare functions

void READ_KCOR_DRIVER(char *kcorFile, char *FILTERS_SURVEY, 
		      double *MAGREST_SHIFT_PRIMARY, 
		      double *MAGPBS_SHIFT_PRIMARY );
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
int  ISBXFILT_KCOR(char *cfilt);
void addFilter_kcor(int ifiltdef, char *NAME, KCOR_FILTERMAP_DEF *MAP);
void init_kcor_indices(void);
void get_MAPINFO_KCOR(char *what, KCOR_MAPINFO_DEF *MAPINFO); 
void filter_match_kcor(char *NAME, int *IFILT_REST, int *IFILT_OBS);
void check_duplicate_filter(char *FRAME, int IFILTDEF, char *FILTER_NAME );

void loadFilterTrans_kcor(int IFILTDEF, int NBL, 
			  float *ARRAY_LAM, float *ARRAY_TRANS,
			  KCOR_FILTERMAP_DEF *MAP);
// end

