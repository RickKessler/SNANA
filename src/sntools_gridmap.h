// Created July 2021 [moved from sntools.h]
// Apr 2024: change typedef GRIDMAP to GRIDMAP_DEF (follow SNANA convention)


// define prototype for multi-dimensionl grid; used for interpolation
#define MXDIM_GRIDMAP 20
typedef struct GRIDMAP_DEF {
  int     ID;        //    
  int     NDIM;     // Number of dimensions      
  int    *NBIN;     // Number of bins in each dimension   
  double *VALMIN;   // min value in each dimension          
  double *VALMAX;   // max value in each dimension     
  double *VALBIN;   // binsize in each dimension         
  double *RANGE;    // max-min in each dimension              

  int    NFUN ;        // Number of stored functions         
  double **FUNVAL ;    // function value at each grid point
  double  *FUNMIN ;    // min fun-val per function 
  double  *FUNMAX ;    // max fun-val per function 
  int    *INVMAP;      // covert multi-D indices into 1D index
  int  NROW;          // number or rows read from file 
  int  OPT_EXTRAP;    // 1=>snap outside values to edge 
  char VARLIST[80];   // comma-sep list of variables (optional to fill)   
  char *VARNAMES[MXDIM_GRIDMAP];    // array of variable names (internally computed)

  float MEMORY; // alloated memory, MB

} GRIDMAP_DEF ;


typedef struct GRIDMAP1D_DEF {
  int NBIN ;
  double  *XVAL, *YVAL ;
} GRIDMAP1D_DEF ;


#define IDGRIDMAP_SIMEFFMAP              8  // for MLCS-AV prior 
#define IDGRIDMAP_HOSTLIB_WGTMAP        20  // HOSTLIB weight map 
#define IDGRIDMAP_KCOR_LCMAG            21  // rest-frame mags for K-cor
#define IDGRIDMAP_KCOR_MWXT             22
#define IDGRIDMAP_KCOR_AVWARP           23
#define IDGRIDMAP_KCOR_VAL              24  // K-cor values
#define IDGRIDMAP_XTMAG                 25  // host-gal extinction map
#define IDGRIDMAP_SPECEFF_OFFSET        30  // id = OFFSET + imap 
#define IDGRIDMAP_zHOST_OFFSET          40  // id = OFFSET + imap 
#define IDGRIDMAP_PHOTPROB_OFFSET       50  // id = OFFSET + imap 
#define IDGRIDMAP_GENPDF                60  // populations
#define IDGRIDMAP_FLUXERRMODEL_OFFSET  100  // id = OFFSET + imap 


// ------ index mapping        
void clear_1DINDEX(int ID);
void  init_1DINDEX(int ID, int NDIM, int *NPT_PERDIM );
int    get_1DINDEX(int ID, int NDIM, int *indx );

// mangled functions for fortran                                                
void clear_1dindex__(int *ID);
void init_1dindex__(int *ID, int *NDIM, int *NPT_PERDIM );
int  get_1dindex__ (int *ID, int *NDIM, int *indx );


void malloc_GRIDMAP(int OPT, GRIDMAP_DEF *gridmap, int NFUN, int NDIM, int MAPSIZE);

void init_interp_GRIDMAP(int ID, char *MAPNAME, int MAPSIZE, int NDIM, int NFUN,
                         int OPT_EXTRAP, 
                         double **GRIDMAP_INPUT, double **GRIDFUN_INPUT, 
                         GRIDMAP_DEF *gridmap ); 
                                                                           
int  interp_GRIDMAP(GRIDMAP_DEF *gridmap, double *data, double *interpFun );
                                                                
void read_GRIDMAP(FILE *fp, char *MAPNAME, char *KEY_ROW, char *KEY_STOP,
                  int IDMAP, int NDIM, int NFUN, int OPT_EXTRAP, int MXROW,
                  char *callFun, GRIDMAP_DEF *GRIDMAP_LOAD ); 

// === END ===
