/*************************
  Created Apr 15 2013 by R.Kessler

  Functions in caps are for user-interface.
  lower-case functions are internal.

************************/

#define MXTRAINFILE_NEARNBR  10 // max number of training files to catenate
#define MXVAR_NEARNBR        5  // max number of variables to use
#define MXBIN_SEPMAX_NEARNBR 100000 
#define MXTRUETYPE       1000  // max TRUETYPE value
#define NTRUETYPE_MAX    100   // max number of different true TYPES
#define HIDOFF_NEARNBR   800   // histogram ID offset

#define ID1D_CELLMAP_NEARNBR       10   // for translating to 1D index
#define BUFFSIZE_CELLMAP_NEARNBR  200   // realloc buf size

// define stupid params because HBOOK title limit is 80 chars
// to hold name of training file
#define MXCHAR_TITLE 80
#define NSPLIT_TITLE  3  // --> allows up to 3x80=240 chars


void NEARNBR_INIT(void);
void nearnbr_init__(void);

void NEARNBR_SET_TRAINPATH(char *path) ;
void nearnbr_set_trainpath__(char *path);

void NEARNBR_SET_TRAINFILE(char *file, float SCALE_NON1A);
void nearnbr_set_trainfile__(char *file, float *SCALE_NON1A);

void NEARNBR_SET_TRUETYPE(char *varName, int truetype_SNIa) ;
void nearnbr_set_truetype__(char *varName, int *truetype_SNIa);

void NEARNBR_SET_SEPMAX(char *varName, double *SEPMAX) ;
void nearnbr_set_sepmax__(char *varName, double *SEPMAX) ;

void NEARNBR_CELLMAP_INIT(int NCHOP);
void getInfo_CELLMAP(int OPT, float *VAL_ARRAY, int *ICELL_1D );
void realloc_NEARNBR_CELLMAP(int ICELL_1D);

void NEARNBR_SET_ODDEVEN(void);
void nearnbr_set_oddeven__(void);

void NEARNBR_INIT2(int ISPLIT);
void nearnbr_init2__(int *ISPLIT);

void NEARNBR_LOADVAL  (char *ccid, char *varName, double d_val) ;
void nearnbr_loadval__(char *ccid, char *varName, double *d_val) ;

void NEARNBR_APPLY(char *CCID);
void nearnbr_apply__(char *CCID);

void NEARNBR_GETRESULTS(char *CCID, int *ITYPE, int *NTYPE, int *ITYPE_LIST, 
			int *NCELL_TRAIN_LIST );
void nearnbr_getresults__(char *CCID, int *ITYPE, int *NTYPE, int *ITYPE_LIST,
			  int *NCELL_TRAIN_LIST ) ;

void  nearnbr_check_inputs(void) ;
void  nearnbr_read_trainLib(int ifile);
int   nearnbr_storeODD_trainLIB(int NROWTOT_NEW);
void  nearnbr_apply_trainLib();
void  nearnbr_init_SEPMAX(void) ;
void  nearnbr_storeTrueTypes(void);
void  nearnbr_preAnal_verify(void) ;
void  nearnbr_reset(void) ;
float nearnbr_SQDIST(int isep, int itrain) ;
int   nearnbr_whichType(int NTYPE, int *NCUTDIST, int *TYPE_CUTPROB ) ;
void  nearnbr_init_SUBSET(void) ;
void  nearnbr_fill_SUBSET_TRAIN(char *CCID); 
void  nearnbr_fill_SUBSET_APPLY(char *CCID); 
void  nearnbr_makeHist(int ISPLIT) ;
void  nearnbr_makeHist_once(void) ;
void  nearnbr_makeHist_allJobs(void) ;
void  nearnbr_fillHist(int ISEPMAX, int ISPARSE_TYPE) ;

void  nearnbr_TRAIN_FILENAME(int ifile, char *TRAIN_FILENAME);

// ================================

int NN_TRAINFLAG ;  // 1=> multiple SEPMAX bins for training; 0=> analysis
int NN_APPLYFLAG ;  // 1=> apply mode; 2=internal apply mode to training

struct {
  int DOFLAG;
  int NCHOP_PER_VAR;    // chop each dimention into this many pieces
  int NVAR, NCELL_TOT ; // NCELL_TOT = NCHOP^NDIM

  float VAL_MIN[MXVAR_NEARNBR];
  float VAL_MAX[MXVAR_NEARNBR];
  float VAL_BIN[MXVAR_NEARNBR];

  int *NLIST;          // list size for each cell
  int **ITRAIN_LIST ;      // train sample CID[icell][i=0,NLIST-1]
  int MEMTOT;

} NEARNBR_CELLMAP ;

struct NEARNBR_INPUTS {

  char   TRAINFILE_PATH[200];
  char   TRAINFILE_LIST[MXTRAINFILE_NEARNBR][200] ;
  int    NTRAINFILE ;
  char   VARNAME_TRUETYPE[60];
  int    TRUETYPE_SNIa ; // used with SCALE_NON1A
  int    TRAIN_ODDEVEN;  // if True, odd SNID for ref, even for TRAIN

  int    NVAR ;
  char   VARNAMES[MXVAR_NEARNBR][40] ;
  float  SEPMAX_MIN[MXVAR_NEARNBR]   ; 
  float  SEPMAX_MAX[MXVAR_NEARNBR]   ; 
  float  SEPMAX_BIN[MXVAR_NEARNBR]   ; 

  float  CUTPROB ;      // fraction of R<1 events to set tag
  float  NSIGMA_PROB ;  // required significance (sigma) of PROB > CUTPROB
  float  SCALE_NON1A ;  // extra sim scale to enhance CC sample;
                        // default=1 for physical CC/Ia rate-ratio

  int FILLHIST ;  // internally set to TRUE if multiple SEPMAX bins
} NEARNBR_INPUTS ;



struct NEARNBR_TRAINLIB {

  int    *CID_VALUES;   // note this int is for sims only
  float  *FITRES_VALUES[MXVAR_NEARNBR];
  int    *TRUETYPE ;  // true TYPE for each entry in the training set
  // xxx mark delete  float  **P_TRAIN; // train-prob for type & each event
  int     NTOT ; // total number of entries stored in FITRES_FILE
  int     NTOT_USE ; // number with valid TRUETYPE

  int MEMTOT ;       // total bytes allocated to memory for TRAINFILE(s)

  int NTRUETYPE ;   // Number of different true types in train set
  int TRUETYPE_LIST[NTRUETYPE_MAX] ; // sparse list of true types
  int TRUETYPE_MAP[MXTRUETYPE] ;     // sparse index vs. TRUE TYPE


  // number of each TRUETYPE
  int NSN[NTRUETYPE_MAX] ;

  // raw input purity without any training.
  double UNTRAINED_PURITY[NTRUETYPE_MAX];


  // define training sparse list for subset inside max k-sphere.
  // These variable are re-determined for each SN in ref sample.
  int NSUBSET ;
  int *ITRAIN ; // .ITRAIN[isubset] = itrain in full training sample

} NEARNBR_TRAINLIB ;



int   NBINTOT_SEPMAX_NEARNBR ;
float **NEARNBR_LIST_SEPMAX ;   // SEPMAX   vs. [ivar][isep]
float **NEARNBR_LIST_SQSEPMAX ; // SQSEPMAX vs. [ivar][isep]
int   **NEARNBR_LIST_NTYPE ;    // NTYPE  vs. [itype][isep]



struct NEARNBR_STORE {
  int   NVAL_LOAD ; // number of values loaded into VALUES
  float VALUE_LOAD[MXVAR_NEARNBR] ;   // loaded VALUE for each variables
  int   NLOAD[MXVAR_NEARNBR] ;   // should be 1 per variable

  int   TRUETYPE_LOAD ;          // actual TYPE of current SN
  int   TRUETYPE_SPARSE ;        // sparse index for above.

  float *SQSEP[MXVAR_NEARNBR] ;  // (VAL_SN - VAL_i)**2 for each trainSN

} NEARNBR_STORE ;


// RESULTS structure is filled only if there is 1 SEPMAX bin
typedef struct NEARNBR_RESULTS_DEF {
  int    ITYPE ;    // integer type (-1 => no type; -9 => not evaluated)
  int    NCELL_TOT ; // total number of reference SN in cell with SQDIST<1
  int    NCELL[MXTRUETYPE] ;    // NCELL vs. truetype sparse index
} NEARNBR_RESULTS_DEF ;

// define struct with training cuts, and again with training+NN cuts.
// ITYPE is the same for both, but the NCELL can be different
NEARNBR_RESULTS_DEF NEARNBR_RESULTS_TRAIN ; // with training cuts
NEARNBR_RESULTS_DEF NEARNBR_RESULTS_FINAL ; // with training+NN cuts
