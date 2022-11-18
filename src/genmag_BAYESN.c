/**********************************************

  July 1 2022 G. Narayan

********************************************/

#include "math.h"
#include "gsl/gsl_linalg.h"
#include "fitsio.h"
#include "sntools.h"
#include "genmag_SEDtools.h"
#include  "genmag_BAYESN.h"
// #include "sntools_modelgrid.h" 
// #include "sntools_genSmear.h" // Aug 30 2019



// ============ MANGLED FORTRAN FUNCTIONS =============
// GN - I think we will likely need these - copied from bayesn
int init_genmag_bayesn__( char *version, int *optmask) {
  int istat;
  istat = init_genmag_BAYESN ( version, *optmask) ;
  return istat;
}


/* int genmag_bayesn__(int *ifilt, double *stretch, int *nobs, 
		  double *Trest, double *rest_mags, double *rest_magerrs ) {
  int istat;
  istat = genmag_bayesn(*ifilt, *stretch, *nobs, 
			Trest, rest_mags, rest_magerrs ) ;
  return istat;
} */

// void get_lamrange_bayesn__(double *lammin, double *lammax) {
//  get_LAMRANGE_bayesn(lammin,lammax);
// }

int init_genmag_BAYESN(char *version, int optmask){

    int  ised;
    int  retval = 0   ;
    int  ABORT_on_LAMRANGE_ERROR = 0;
    int  ABORT_on_BADVALUE_ERROR = 1;
    //char BANNER[120], tmpFile[200], sedcomment[40], version[60]  ;

    char fnam[] = "init_genmag_BAYESN";
    // -------------- BEGIN --------------

    // extrac OPTMASK options


    sprintf(BANNER, "%s : Initialize %s", fnam, version );
    print_banner(BANNER);
 
    // HACK HACK HACK 
    BAYESN_MODEL_INFO.n_lam_knots = 20;
    malloc_double2D(1, 20, 20, &BAYESN_MODEL_INFO.W0 ); 
    char SED_filepath[] = "/global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SNDATA_ROOT/snsed/Hsiao07.dat";
    int istat;
    SEDMODEL_FLUX_DEF *S0 = &BAYESN_MODEL_INFO.S0;
    malloc_SEDFLUX_SEDMODEL(S0,0,0,0);
    double Trange[2] = {-20.0, 90.0};
    double Lrange[2] = {1500.0, 16000.0};
    
    istat = rd_sedFlux(SED_filepath, "blah test rd_sedFlux", Trange, Lrange
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0
	       ,&S0->NDAY, S0->DAY, &S0->DAYSTEP
	       ,&S0->NLAM, S0->LAM, &S0->LAMSTEP
	       ,S0->FLUX,  S0->FLUXERR );
    printf("XXX istat %d\n", istat);

    if ( NFILT_SEDMODEL == 0 ) {
      sprintf(c1err,"No filters defined ?!?!?!? " );
      sprintf(c2err,"Need to call init_filter_SEDMODEL");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    //ABORT_on_LAMRANGE_ERROR = ( OPTMASK & OPTMASK_BAYESN_ABORT_LAMRANGE ) ;
    
    filtdump_SEDMODEL();

    // Hack wavelength ranges (R.Kessler) ... these should be read from model   
    SEDMODEL.LAMMIN_ALL =  2000.0 ;  // rest-frame SED range                    
    SEDMODEL.LAMMAX_ALL = 15000.0 ;
    SEDMODEL.RESTLAMMIN_FILTERCEN =  3000.0 ; // rest-frame central wavelength range
    SEDMODEL.RESTLAMMAX_FILTERCEN = 14000.0 ;
    
    printf("XXXX %s Hello from fortran hell\n", fnam);
    debugexit(fnam);

    return 0;

} // end init_genmag_BAYESN

// =====================================================
void genmag_BAYESN(
		  int OPTMASK     // (I) bit-mask of options (LSB=0)
		  ,int ifilt_obs  // (I) absolute filter index
		  ,double *parList_SN   // DLMAG, THETA, AV, RV
		  ,double mwebv   // (I) Galactic extinction: E(B-V)
		  ,double z       // (I) Supernova redshift
		  ,int    Nobs         // (I) number of epochs
		  ,double *Tobs_list   // (I) list of Tobs (w.r.t peakMJD) 
		  ,double *magobs_list  // (O) observed mag values
		  ,double *magerr_list  // (O) model mag errors
		  ) {
    int o;
    double mag;
    char fnam[] = "genmag_BAYESN";

    // ------- BEGIN -----------
    double zdum = 2.5*log10(1.0+z);
    for (o = 0; o < Nobs; o++) {
      mag   = fabs(.2*Tobs_list[o]) + 20.0 + zdum ;
      magobs_list[o] = mag;
      magerr_list[o] = 0.1;
    }

    return;

} //End of genmag_BAYESN
 
void genmag_bayesn__(int *OPTMASK, int *ifilt_obs, double *parlist_SN,
	       	double *mwebv, double *z, int *Nobs,
	       	double *Tobs_list, double *magobs_list,
	       	double *magerr_list) {
	genmag_BAYESN(*OPTMASK, *ifilt_obs, parlist_SN, *mwebv, *z,
		       	*Nobs, Tobs_list, magobs_list, magerr_list);
	return;
}

gsl_matrix *invKD_irr(int Nk, double *xk) {
	//FIXME I'm pretty sure this whole thing could be done better by using 
	//the GSL tridiagonal solver. This is just a direct port of my python.

	gsl_matrix * K = gsl_matrix_alloc(Nk-2, Nk-2);
	gsl_matrix * D = gsl_matrix_alloc(Nk-2, Nk);
	gsl_matrix * M = gsl_matrix_alloc(Nk, Nk);

	gsl_matrix_set_zero(K);
	gsl_matrix_set_zero(D);
	gsl_matrix_set_zero(M);

	gsl_matrix_set(K, 0, 0, (xk[2] - xk[0])/3.0);
	gsl_matrix_set(K, 0, 1, (xk[2] - xk[1])/6.0);
	gsl_matrix_set(K, Nk-3, Nk-4, (xk[Nk-2] - xk[Nk-3])/6.0);
	gsl_matrix_set(K, Nk-3, Nk-3, (xk[Nk-1] - xk[Nk-3])/3.0);

	int j, r;
	for (j=2; j<Nk-2; j++) {
		r = j - 1;
		gsl_matrix_set(K, r, r-1, (xk[j] - xk[j-1])/6.0);
		gsl_matrix_set(K, r, r, (xk[j+1] - xk[j-1])/3.0);
		gsl_matrix_set(K, r, r+1, (xk[j+1] - xk[j])/6.0);
	}
	for (j=1; j<Nk-1; j++) {
		r = j - 1;
		gsl_matrix_set(D, r, r, 1.0/(xk[j] - xk[j-1]));
		gsl_matrix_set(D, r, r+1, -1.0/(xk[j+1] - xk[j]) - 1.0/(xk[j] - xk[j-1]));
		gsl_matrix_set(D, r, r+2, 1.0/(xk[j+1] - xk[j]));
	}

	int c, s;

	gsl_permutation * p = gsl_permutation_alloc(Nk-2);
	gsl_linalg_LU_decomp(K, p, &s);

	for (c=0; c<Nk; c++) {
		//Solve for 1st to (Nk-2)th elements of cth column of M
		gsl_vector_view d = gsl_matrix_column(D, c);
		gsl_vector_view m = gsl_matrix_subcolumn(M, c, 1, Nk-2);
		gsl_linalg_LU_solve(K, p, &d.vector, &m.vector);
	}

	return M;
}

gsl_matrix *spline_coeffs_irr(int N, int Nk, double *x, double *xk, gsl_matrix *invKD) {
	gsl_matrix * J = gsl_matrix_alloc(N, Nk);
	gsl_matrix_set_zero(J);

	int i, j, q;
	double h, a, b, c, d, f;
	for (i=0; i<N; i++) {
		if (x[i] > xk[Nk-1]) {
			h = xk[Nk-1] - xk[Nk-2];
			a = (xk[Nk-1] - x[i])/h;
			b = 1.0 - a;
			f = (x[i] - xk[Nk-1])*h/6.0;

			gsl_matrix_set(J, i, Nk-2, a);
			gsl_matrix_set(J, i, Nk-1, b);
			//FIXME I can probably do this by accessing a row and modifying
			for (j=0; j<Nk; j++) {
				gsl_matrix_set(J, i, j, gsl_matrix_get(J, i, j) + f*gsl_matrix_get(invKD, Nk-2, j));
			}
		}
		else if (x[i] < xk[0]) {
			h = xk[1] - xk[0];
			b = (x[i] - xk[0])/h;
			a = 1.0 - b;
			f = (x[i] - xk[0])*h/6.0;

			gsl_matrix_set(J, i, 0, a);
			gsl_matrix_set(J, i, 1, b);
			for (j=0; j<Nk; j++) {
				gsl_matrix_set(J, i, j, gsl_matrix_get(J, i, j) - f*gsl_matrix_get(invKD, 1, j));
			}
		}
		else {
			q = 0;
			while (q < Nk && xk[q+1] <= x[i]) { q++; }
			h = xk[q+1] - xk[q];
			a = (xk[q+1] - x[i])/h;
			b = 1.0 - a;
			c = ((pow(a, 3) - a)/6.0)*h*h;
			d = ((pow(b, 3) - b)/6.0)*h*h;

			gsl_matrix_set(J, i, q, a);
			gsl_matrix_set(J, i, q+1, b);
			for (j=0; j<Nk; j++) {
				gsl_matrix_set(J, i, j, gsl_matrix_get(J, i, j) + c*gsl_matrix_get(invKD, q, j) + d*gsl_matrix_get(invKD, q+1, j));
			}
		}
	}

	return J;
}
