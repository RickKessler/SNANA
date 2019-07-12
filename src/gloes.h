// gloes.h

extern int gloes_one_sigma( double *xp, double *yp, double *wp, int np, 
			    double *x, int n,  double sigma, double *y, double *ey);

extern int gloes_n_sigma( double *xp, double *yp, double *wp, int np, 
			  double *x, int n, double *sigma, double *y, double *ey);

extern int gloes_n( double *xp, double *yp, double *wp, int np, 
		    double *x, int n, int N, double *y, double *ey);

extern int gloes2D_n_sigma( double *xp, double *yp, double *zp, double *wp, 
			    int np, double *x, double *y, int n, 
			    double *sigmax, double *sigmay, double *z, double *ez);
