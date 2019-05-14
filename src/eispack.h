int bakvec ( int n, double t[], double e[], int m, double z[] );
void balanc ( int n, double a[], int *low, int *igh, double scale[] );
void balbak ( int n, int low, int igh, double scale[], int m, double z[] );
void bandr ( int n, int mb, double a[], double d[], double e[], double e2[], 
  bool matz, double z[] );
int bandv ( int n, int mbw, double a[], double e21, int m, double w[], 
  double z[] );
int bisect ( int n, double *eps1, double d[], double e[], double e2[], double lb, 
  double ub, int mm, int *m, double w[], int ind[] );
int bqr ( int nm, int n, int mb, double a[], double *t, double *r );
void cbabk2 ( int n, int low, int igh, double scale[], int m, double zr[], 
  double zi[] );
void cbal ( int n, double ar[], double ai[], int *low, int *igh, double scale[] );
void cdiv ( double ar, double ai, double br, double bi, double *cr, double *ci );
int cg_lr ( int n, double ar[], double ai[], double wr[], double wi[], bool matz, 
  double zr[], double zi[] );
int cg_qr ( int n, double ar[], double ai[], double wr[], double wi[], bool matz, 
  double zr[], double zi[] );
int ch ( int n, double ar[], double ai[], double w[], bool matz, double zr[], 
  double zi[] );
int ch3 ( int n, double a[], double w[], bool matz, double zr[], double zi[] );
int cinvit ( int n, double ar[], double ai[], double wr[], double wi[], 
  bool select[], int mm, int *m, double zr[], double zi[] );
void combak ( int n, int low, int igh, double ar[], double ai[], int inter[], 
  int m, double zr[], double zi[] );
void comhes ( int n, int low, int igh, double ar[], double ai[], int inter[] );
int comlr ( int n, int low, int igh, double hr[], double hi[], double wr[], 
  double wi[] );
int comlr2 ( int n, int low, int igh, int inter[], double hr[], double hi[], 
  double wr[], double wi[], double zr[], double zi[] );
int comqr ( int n, int low, int igh, double hr[], double hi[], double wr[], 
  double wi[] );
int comqr2 ( int n, int low, int igh, double ortr[], double orti[], 
  double hr[], double hi[], double wr[], double wi[], double zr[], 
  double zi[] );
void cortb ( int n, int low, int igh, double ar[], double ai[], double ortr[],
  double orti[], int m, double zr[], double zi[] );
void corth ( int n, int low, int igh, double ar[], double ai[], double ortr[], 
  double orti[] );
void csroot ( double xr, double xi, double *yr, double *yi );
void elmbak ( int n, int low, int igh, double a[], int ind[], int m, double z[] );
void elmhes ( int n, int low, int igh, double a[], int ind[] );
void eltran ( int n, int low, int igh, double a[], int ind[], double z[] );
int figi ( int n, double t[], double d[], double e[], double e2[] );
int figi2 ( int n, double t[], double d[], double e[], double z[] );
int hqr ( int n, int low, int igh, double h[], double wr[], double wi[] );
int hqr2 ( int n, int low, int igh, double h[], double wr[], double wi[], 
  double z[] );
void htrib3 ( int n, double a[], double tau[], int m, double zr[], double zi[] );
void htribk ( int n, double ar[], double ai[], double tau[], int m, double zr[], 
  double zi[] );
void htrid3 ( int n, double a[], double d[], double e[], double e2[], 
  double tau[] );
void htridi ( int n, double ar[], double ai[], double d[], double e[], 
  double e2[], double tau[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_print ( int n, int a[], char *title );
int imtql1 ( int n, double d[], double e[] );
int imtql2 ( int n, double d[], double e[], double z[] );
int imtqlv ( int n, double d[], double e[], double e2[], double w[], int ind[] );
int invit ( int n, double a[], double wr[], double wi[], bool select[], 
  int mm, int *m, double z[] );
int minfit ( int nm, int m, int n, double a[], double w[], int ip, double b[] );
void ortbak ( int n, int low, int igh, double a[], double ort[], int m, 
  double z[] );
void orthes ( int n, int low, int igh, double a[], double ort[] );
void ortran ( int n, int low, int igh, double a[], double ort[], double z[] );
double pythag ( double a, double b );
void qzhes ( int n, double a[], double b[], bool matz, double z[] );
int qzit ( int n, double a[], double b[], double eps1, bool matz, double z[] );
void qzval ( int n, double a[], double b[], double alfr[], double alfi[], 
  double beta[], bool matz, double z[] );
void qzvec ( int n, double a[], double b[], double alfr[], double alfi[], 
  double beta[], double z[] );
double r8_epsilon ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
void r8mat_identity  ( int n, double a[] );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8mat_zeros ( int m, int n, double a[] );
void r8vec_print ( int n, double a[], char *title );
void r8vec2_print ( int n, double a1[], double a2[], char *title );
int ratqr ( int n, double eps1, double d[], double e[], double e2[], int m, 
  double w[], int ind[], double bd[], bool type, int idef );
void rebak ( int n, double b[], double dl[], int m, double z[] );
void rebakb ( int n, double b[], double dl[], int m, double z[] );
int reduc ( int n, double a[], double b[], double dl[] );
int reduc2 ( int n, double a[], double b[], double dl[] );
int rg_elm ( int n, double a[], double wr[], double wi[], bool matz, 
  double z[] );
int rg_ort ( int n, double a[], double wr[], double wi[], bool matz, 
  double z[] );
int rgg ( int n, double a[], double b[], double alfr[], double alfi[], 
  double beta[], bool matz, double z[] );
int rs ( int n, double a[], double w[], bool matz, double z[] );
int rsb ( int n, int mb, double a[], double w[], bool matz, double z[] );
int rsg ( int n, double a[], double b[], double w[], bool matz, double z[] );
int rsgab ( int n, double a[], double b[], double w[], bool matz, double z[] );
int rsgba ( int n, double a[], double b[], double w[], bool matz, double z[] );
int rsm ( int n, double a[], double w[], int m, double z[] );
int rsp ( int n, int nv, double a[], double w[], bool matz, double z[] );
int rspp ( int n, int nv, double a[], double w[], bool matz, double z[], int m, 
  bool type );
int rst ( int n, double w[], double e[], bool matz, double z[] );
int rt ( int n, double a[], double w[], bool matz, double z[] );
int sturm_sequence ( double d[], double e[], double e2[], int n, int p, int q, 
  double x1 );
int svd ( int m, int n, double a[], double w[], bool matu, double u[], 
  bool matv, double v[] );
void timestamp ( );
int tinvit ( int n, double d[], double e[], double e2[], int m, double w[], 
  int ind[], double z[] );
int tql1 ( int n, double d[], double e[] );
int tql2 ( int n, double d[], double e[], double z[] );
int tqlrat ( int n, double w[], double fv2[] );
void trbak1 ( int n, double a[], double e[], int m, double z[] );
void trbak3 ( int n, int nv, double a[], int m, double z[] );
void tred1 ( int n, double a[], double w[], double fv1[], double fv2[] );
void tred2 ( int n, double a[], double d[], double e[], double z[] );
void tred3 ( int n, int nv, double a[], double d[], double e[], double e2[] );
int tridib ( int n, double *eps1, double d[], double e[], double e2[], 
  double *lb, double *ub, int m11, int m, double w[], int ind[] );
int tsturm ( int n, double *eps1, double d[], double e[], double e2[], double lb, 
  double ub, int mm, int *m, double w[], double z[] );
