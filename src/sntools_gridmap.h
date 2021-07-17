// Created July 2021 [moved from sntools.h]

// define prototype for multi-dimensionl grid; used for interpolation
typedef struct GRIDMAP {
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
  int  OPT_EXTRAP;   // 1=>snap outside values to edge 
  char VARLIST[80]; // comma-sep list of variables (optional to fill)   
} GRIDMAP ;


typedef struct GRIDMAP1D {
  int NBIN ;
  double  *XVAL, *YVAL ;
} GRIDMAP1D ;

// xxx mark delete struct GRIDMAP  SIMEFF_GRIDMAP ;

#define IDGRIDMAP_SIMEFFMAP              8  // for MLCS-AV prior 
#define IDGRIDMAP_HOSTLIB_WGTMAP        20  // HOSTLIB weight map 
#define IDGRIDMAP_SPECEFF_OFFSET        30  // id = OFFSET + imap 
#define IDGRIDMAP_zHOST_OFFSET          40  // id = OFFSET + imap 
#define IDGRIDMAP_PHOTPROB_OFFSET       50  // id = OFFSET + imap 
#define IDGRIDMAP_GENPDF                60  // Jun 2020          
#define IDGRIDMAP_FLUXERRMODEL_OFFSET  100  // id = OFFSET + imap 


// ------ index mapping                                                         
void clear_1DINDEX(int ID);
void  init_1DINDEX(int ID, int NDIM, int *NPT_PERDIM );
int    get_1DINDEX(int ID, int NDIM, int *indx );

// mangled functions for fortran                                                
void clear_1dindex__(int *ID);
void init_1dindex__(int *ID, int *NDIM, int *NPT_PERDIM );
int  get_1dindex__ (int *ID, int *NDIM, int *indx );


void malloc_GRIDMAP(int OPT, GRIDMAP *gridmap, int NFUN, int NDIM, int MAPSIZE);

void init_interp_GRIDMAP(int ID, char *MAPNAME, int MAPSIZE, int NDIM, int NFUN,
                         int OPT_EXTRAP, 
                         double **GRIDMAP_INPUT, double **GRIDFUN_INPUT, 
                         GRIDMAP *gridmap ); 
                                                                           
int  interp_GRIDMAP(GRIDMAP *gridmap, double *data, double *interpFun );
                                                                
void read_GRIDMAP(FILE *fp, char *MAPNAME, char *KEY_ROW, char *KEY_STOP,
                  int IDMAP, int NDIM, int NFUN, int OPT_EXTRAP, int MXROW,
                  char *callFun, GRIDMAP *GRIDMAP_LOAD ); 

// === END ===
