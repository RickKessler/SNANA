

#define ORDER_PROB_WRONGHOST_POLY 4
#define DZTRUEMAX_WRONGHOST 0.01   // get zmatch when |ZTRUE-zSN| < this

struct {
  double USER_zMIN, USER_zMAX ;

  // polynom to give PROB_WRONGHOST(zSN)
  double PROB_POLY[ORDER_PROB_WRONGHOST_POLY] ;
  GENPOLY_DEF ZPOLY_PROB ;

  // wrongHost list of ZTRUE/ZMATCH pairs
  int     NLIST ;
  double *ZTRUE_LIST  ;
  double *ZMATCH_LIST ;
  
} WRONGHOST  ;

// ------------------------------------------------------------
// prototype functions
void    INIT_WRONGHOST(char *inFile, double ZMIN, double ZMAX);
double  gen_zHOST(int CID, double zSN, int *istat ) ;

