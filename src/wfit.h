//  from codist.h

typedef struct  {
  /* Cosmological parameter structure used by codist.c */
  double omm;			/* Omega_matter */
  double ome;			/* Omega_energy */
  double w0;			/* equation of state of Dark Energy: */
  double wa;			/*    w(z) = w0 + wa(1-a)  */
} Cosparam;

double EofZ(double z, Cosparam *cptr);
double one_over_EofZ(double z, Cosparam *cptr);
double codist(double z, Cosparam *cptr);
void   test_codist(void);




// from simpint.h

#define EPS 1.0e-6
#define JMAX 20

double simpint(double (*func)(double, Cosparam *), 
	       double x1, double x2, Cosparam *vptr);

double trapint(double (*func)(double, Cosparam *), 
	       double x1, double x2, Cosparam *vptr);

double trapezoid(double (*func)(double, Cosparam *), 
		 double x1, double x2, int n, double s, Cosparam *vptr);
