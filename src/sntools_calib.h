
// Created Oct 2019 by R.Kessler

#define MXPRIMARY_CALIB  6

#define MXFILT_CALIB       MXFILTINDX
#define MXFILT_REST_CALIB  20
#define MXFILT_OBS_CALIB   62

#define MXTABLE_KCOR   100   // max numver of kcor tables

#define MXLAMBIN_FILT_CALIB  4000  // max wave bins for filter trans
// xxx mark #define MXLAMBIN_PRIM  4000  // max wave bins for primary

#define MXTBIN_KCOR   150
#define MXZBIN_KCOR   100
#define MXAVBIN_KCOR  100
#define MXCBIN_AVWARP 100  // max color index for AVWARP table

#define KDIM_T        0
#define KDIM_z        1  
#define KDIM_AV       2
#define KDIM_IFILTr   3
#define KDIM_IFILTo   4
#define NKDIM_KCOR    5
#define N4DIM_KCOR    4
#define OPT_EXTRAP_KCOR 1  // extrap option for GRIDMAP utli

#define IDMAP_KCOR_TABLE   15
#define IDMAP_KCOR_AVWARP  16
#define IDMAP_KCOR_LCMAG   17
#define IDMAP_KCOR_MWXT    18

#define KCOR_CRAZYVAL 6.0 

int KCOR_VERBOSE_FLAG;
int IFILTDEF_BESS_BX;
int NERR_KCOR_AVWARP;

struct {
  bool USE_AVWARPTABLE;  // speed up AVwarp calculation
  bool DUMP_AVWARP, DUMP_MAG, DUMP_KCOR; // stdout dumps
} CALIB_OPTIONS ;


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
  int  OPT_FRAME; // indicates REST or OBS

  // filled by addFilter_kcor
  int  NFILTDEF;
  int  IFILTDEF[MXFILT_CALIB];     // vs. sparse index
  int  IFILTDEF_INV[MXFILT_CALIB]; // vs. absolute index
  int  NDEFINE[MXFILT_CALIB];      // Number of times each IFILTDEF is defined
  int  NFILT_DUPLICATE;           // number of bands with duplicates
  char FILTERSTRING[MXFILT_CALIB]; // list of single-char bands
  char *FILTER_NAME[MXFILT_CALIB]; // full name of each filter, vs. sparse indx
  char *SURVEY_NAME[MXFILT_CALIB];

  double PRIMARY_MAG[MXFILT_CALIB];
  double PRIMARY_ZPOFF_SYN[MXFILT_CALIB]; // from synthetic vs. native
  double PRIMARY_ZPOFF_FILE[MXFILT_CALIB]; // from ZPOFF.DAT file
  int    PRIMARY_KINDX[MXFILT_CALIB];  // index to CALIB_INFO.PRIMARY_XXX[]
  int    NBIN_LAM_PRIMARY ;
  float *PRIMARY_LAM, *PRIMARY_FLUX;

  // below is filled by loadFilterTrans_kcor
  int    NBIN_LAM[MXFILT_CALIB];
  double TRANS_MAX[MXFILT_CALIB];
  double LAMMEAN[MXFILT_CALIB], LAMRMS[MXFILT_CALIB];
  float *LAM[MXFILT_CALIB], *TRANS[MXFILT_CALIB];
  double LAMRANGE[MXFILT_CALIB][2]; // range of TRANS > 0

  // define rest-frame LAMRANGE for kcor lookup
  double LAMRANGE_KCOR[MXFILT_CALIB][2]; 

} FILTERCAL_DEF ;


struct CALIB_INFO {

  // info passed to driver
  char FILENAME[MXPATHLEN] ;
  fitsfile *FP ;

  char FILTERS_SURVEY[MXFILT_CALIB]; // filter list read from SIMLIB file
  int  NFILTDEF_SURVEY ;
  
  double MAGREST_SHIFT_PRIMARY[MXFILT_CALIB];
  double MAGOBS_SHIFT_PRIMARY[MXFILT_CALIB];

  // - - - - info read from FITS file - - - - - - 

  // header info
  int   VERSION, NFLAMSHIFT;

  int   NPRIMARY;
  char  *PRIMARY_NAME[MXPRIMARY_CALIB] ;
  int    PRIMARY_INDX[MXFILT_CALIB];
  double PRIMARY_MAG[MXFILT_CALIB];          // native mag
  double PRIMARY_ZPOFF_SYN[MXFILT_CALIB];    // mag(native) - mag(synth)
  double PRIMARY_ZPOFF_FILE[MXFILT_CALIB];   // from ZPOFF.DAT

  int   NFILTDEF;
  char *FILTER_NAME[MXFILT_CALIB] ;
  char *SURVEY_NAME[MXFILT_CALIB];
  int   IFILTDEF[MXFILT_CALIB] ; 
  int   MASK_FRAME_FILTER[MXFILT_CALIB]; // bits 0,1 --> rest, obs
  int   MASK_EXIST_BXFILT; // logical for BX existing for rest/obs frame
  bool  IS_SURVEY_FILTER[MXFILT_CALIB];

  FILTERCAL_DEF FILTERCAL_REST ;
  FILTERCAL_DEF FILTERCAL_OBS  ;

  int   NKCOR, NKCOR_STORE;
  char *STRING_KCORLINE[MXTABLE_KCOR] ; // full FITS line with KCOR info
  char *STRING_KCORSYM[MXTABLE_KCOR];   // just K_xy part
  bool  EXIST_KCOR[MXFILT_CALIB][MXFILT_CALIB]; // for each [rest][obs] combo
  int   IFILTMAP_KCOR[2][MXTABLE_KCOR];
  int   k_index[MXTABLE_KCOR];

  bool  ISLAMSHIFT[MXFILT_CALIB] ;

  double RVMW;
  int    OPT_MWCOLORLAW ;
  
  KCOR_BININFO_DEF BININFO_LAM;  // ??
  KCOR_BININFO_DEF BININFO_T;    // Trest
  KCOR_BININFO_DEF BININFO_z;    //
  KCOR_BININFO_DEF BININFO_AV ;  // AVwarp param to match model color to data
  KCOR_BININFO_DEF BININFO_C ;   // observed color (not SALT2c)

  KCOR_MAPINFO_DEF MAPINFO_KCOR ;
  KCOR_MAPINFO_DEF MAPINFO_AVWARP ;
  KCOR_MAPINFO_DEF MAPINFO_LCMAG ;
  KCOR_MAPINFO_DEF MAPINFO_MWXT ;

  double zRANGE_LOOKUP[2] ;

  int NBIN_KCOR_TABLE, NBIN_LCMAG_TABLE, NBIN_AVWARP_TABLE, NBIN_MWXT_TABLE;
  float *KCORTABLE1D_F;     // use float to save memory
  float *LCMAG_TABLE1D_F ;
  float *AVWARP_TABLE1D_F ;
  float *MWXT_TABLE1D_F ;
  float *FLUX_SNSED_F;      // use float to save memory

  // misc init info
  int NCALL_READ ;
  bool STANDALONE;

} CALIB_INFO ;


double **TEMP_KCOR_ARRAY;

struct {
  struct GRIDMAP GRIDMAP_LCMAG;
  struct GRIDMAP GRIDMAP_MWXT;
  struct GRIDMAP GRIDMAP_AVWARP ;
  struct GRIDMAP GRIDMAP_KCOR;
} KCOR_TABLE ;


  // ============================== 
// declare functions

void READ_CALIB_DRIVER(char *kcorFile, char *FILTERS_SURVEY, 
		       double *MAGREST_SHIFT_PRIMARY, 
		       double *MAGOBS_SHIFT_PRIMARY );

void read_calib_driver__(char *kcorFile, char *FILTERS_SURVEY,
			double *MAGREST_SHIFT_PRIMARY,
			double *MAGOBS_SHIFT_PRIMARY );

void init_calib_options(void);

void read_calib_init(void);
void read_calib_open(void);
void read_calib_head(void);
void read_calib_zpoff(void);
void read_calib_snsed(void);
void read_calib_filters(void);
void read_calib_primarysed(void);

void read_kcor_mags(void);
void read_kcor_tables(void);
void read_kcor_binInfo(char *VARNAME, char *VARSYM, int MXBIN,
		       KCOR_BININFO_DEF *BININFO) ;
void fill_kcor_binInfo_C(void);

void print_calib_summary(void);

void parse_KCOR_STRING(char *STRING, 
		       char *strKcor, char *cfilt_rest, char *cfilt_obs);
int  ISBXFILT_KCOR(char *cfilt);
void addFilter_kcor(int ifiltdef, char *NAME, FILTERCAL_DEF *MAP);
void init_kcor_indices(void);
void get_MAPINFO_KCOR(char *what, KCOR_MAPINFO_DEF *MAPINFO); 
void filter_match_kcor(char *NAME, int *IFILT_REST, int *IFILT_OBS);
void check_duplicate_filter(char *FRAME, int IFILTDEF, char *FILTER_NAME );

void loadFilterTrans_kcor(int IFILTDEF, int NBL, 
			  float *ARRAY_LAM, float *ARRAY_TRANS,
			  FILTERCAL_DEF *MAP);

void set_lamrest_range_KCOR(int ifilt);
void set_lamrest_range_UBVRI(int ifilt);

void get_calib_primary(char *primary_name, int *NBLAM, 
		      double *lam, double *flux);
void get_calib_primary__(char *primary_name, int *NBLAM, 
			double *lam, double *flux);

void get_calib_filterTrans(int OPT_FRAME, int ifilt_obs, char *surveyName, 
			  char *filterName, double *magprim, 
			  int *nblam, double *lam, 
			  double *transSN, double *transREF);
void get_calib_filtertrans__(int *OPT_FRAME, int *ifilt_obs, char *surveyName, 
			    char *filterName, double *magprim,
			    int *nblam, double *lam, 
			    double *transSN, double *transREF);

void get_calib_filtlam_stats(int opt_frame, int ifilt_obs,  
			    double *lamavg, double *lamrms,
			    double *lammin, double *lammax);

void get_calib_filtlam_stats__(int *opt_frame, int *ifilt_obs,  
			      double *lamavg, double *lamrms,
			      double *lammin, double *lammax);

void get_calib_filtindex_map(int OPT_FRAME, int *NFILTDEF, int *IFILTDEF_MAP,
			    int *IFILTDEF_INVMAP);
void get_calib_filtindex_map__(int *OPT_FRAME, int *NFILTDEF, int *IFILTDEF_MAP,
			      int *IFILTDEF_INVMAP);

double get_calib_zpoff_file(int ifiltdef);
double get_calib_zpoff_file__(int *ifiltdef);

// K-cor functions for MLCS, snoopy ...

void PREPARE_KCOR_TABLES(void);
void prepare_kcor_tables__(void);
void prepare_kcor_table_LCMAG(void);
void prepare_kcor_table_MWXT(void);
void prepare_kcor_table_AVWARP(void);
double fit_AVWARP(int ifiltdef_a, int ifiltdef_b, double T, double C);
void prepare_kcor_table_KCOR(void);

int nearest_ifiltdef_rest( int opt, int ifiltdef, int rank, double z, char *callFun,
			   double *lamdif_min );
int nearest_ifiltdef_rest__(int *opt, int *ifiltdef, int *rank, double *z, char *callFun,
			    double *lamdif_min );

// ?? void get_KCOR_FILTERCAL(int OPT_FRAME, char *fnam, FILTERCAL_DEF *MAP );
double GET_KCOR_DRIVER(int IFILT_OBS, int *IFILT_REST_LIST, 
		      double *MAG_REST_LIST, double *LAMDIF_LIST,
		      double Trest, double z, double *AVwarp);

double get_kcor_driver__(int *IFILT_OBS, int *IFILT_REST_LIST, 
			 double *MAG_REST_LIST, double *LAMDIF_LIST,
			 double *Trest, double *z, double *AVwarp);

double eval_kcor_table_LCMAG(int ifiltdef_rest, double Trest, double z, double AVwarp);

double eval_kcor_table_MWXT(int ifiltdef_obs, double Trest, double z, double AVwarp,
		     double MWEBV, double RV, int OPT_MWCOLORLAW);

double eval_kcor_table_AVWARP(int ifiltdef_a, int ifiltdef_b, double mag_a,double mag_b,
		       double Trest, int *istat);

double eval_kcor_table_value(int ifiltdef_rest, int ifiltdef_obs, double Trest, 
		      double z, double AVwarp);

// END

