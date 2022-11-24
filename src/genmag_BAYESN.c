/**********************************************

  July 1 2022 G. Narayan

********************************************/

#include "math.h"
#include "gsl/gsl_linalg.h"
#include "fitsio.h"
#include "sntools.h"
#include "genmag_SEDtools.h"
#include  "genmag_BAYESN.h"


#ifdef USE_BAYESN
#include "yaml.h"
#endif


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


void read_BAYESN_inputs(char *filename)
{
    char fnam[] = "read_BAYESN_inputs";

    // -------------- BEGIN -------------
#ifdef USE_BAYESN
    FILE *fh = fopen(filename, "r");

    // declare YAML parser and event instances
    yaml_parser_t parser;
    yaml_event_t  event;  
    yaml_event_t last_event;

    // this is just a variable to specify the kind of data we are currently 
    // reading from the YAML file 
    int datatype = 0; // 0: we don't know 1: scalar 2: vector 3: matrix
                      // in principle we can support strings and other types
    // if we are reading a matrix from the YAML file, we need to keep track of row/col
    int col=0;
    int row=0;
    int rowsize = 0; // and rowsize 
                     //
    // default init for the sizes is negative to deliberately force an error
    // if we don't populate from the YAML
    // note that they are read as double because it is easier to assume 
    // that everything in the YAML is a float and just fix it later
    double N_LAM =  -1.0;
    double N_TAU =  -1.0;
    double N_SIG =  -1.0;

    double *L_Sigma_epsilon;
    double *W0;
    double *W1;
    
    // we need something to store the current scalar value from the YAML file
    double this_scalar = 0.0;
    // and something to point to the current BayeSN variable being populated
    // this only works if datatype is in 1-3 (assumed double)
    double *bayesn_var_dptr = &this_scalar;


    sprintf(BANNER, "%s : Begin reading BAYESN model components from %s", fnam, filename);
    print_banner(BANNER);

    /* Initialize parser */
    if(!yaml_parser_initialize(&parser))
      fputs("Failed to initialize parser!\n", stderr);
    if(fh == NULL)
      fputs("Failed to open file!\n", stderr);

    /* Set input file */
    yaml_parser_set_input_file(&parser, fh);
    
    /* Start parsing the input file */
    do {
      // everything action (open/read etc) is an event
      if (!yaml_parser_parse(&parser, &event)) {
         printf("Parser error %d\n", parser.error);
         exit(EXIT_FAILURE);
      }
    
      // We check what kind of event we get 
      switch(event.type)
      {
      // blank line
      case YAML_NO_EVENT: puts("No event!"); break;
      case YAML_STREAM_START_EVENT: puts("Reading BayeSN YAML file"); break;
      case YAML_STREAM_END_EVENT:   puts("Done loading from BayeSN YAML file");   break;
    
      /* Block delimeters - I don't actually need to do anything with these events */
      /*
      case YAML_DOCUMENT_START_EVENT: puts("<b>Start Document</b>"); break;
      case YAML_DOCUMENT_END_EVENT:   puts("<b>End Document</b>");   break;
      case YAML_MAPPING_START_EVENT:  puts("<b>Start Mapping</b>");  break;
      case YAML_MAPPING_END_EVENT:    puts("<b>End Mapping</b>");    break;
      case YAML_ALIAS_EVENT:  
          printf("Got alias (anchor %s)\n", event.data.alias.anchor); 
          break;
      */
    
      // the events we care about are all "Scalar" events
      // even if the data being read is a vector or a matrix, it is parsed element-by-element
      case YAML_SCALAR_EVENT: 
          // we have to decide how to handle each event 
          // we need to require three events to define the sizes of the rest of the arrays
          // these HAVE to be the first three events in the YAML file
          // but within those three, they can be in any order
          if (strcmp(event.data.scalar.value, "N_LAM")==0)
          {
              datatype = 1;
              bayesn_var_dptr = &N_LAM;
              break;
          }
          if (strcmp(event.data.scalar.value, "N_TAU")==0)
          {
              datatype = 1;
              bayesn_var_dptr = &N_TAU;
              break;
          }
          if (strcmp(event.data.scalar.value, "N_SIG")==0)
          {
              datatype = 1; 
              bayesn_var_dptr = &N_SIG;
              break;
          }
    
          // next we'll define how to handle the scalars
          if (strcmp(event.data.scalar.value, "M0")==0)
          {
              datatype = 1;
              bayesn_var_dptr = &BAYESN_MODEL_INFO.M0;
              break;
          }
          if (strcmp(event.data.scalar.value, "SIGMA0")==0)
          {
              datatype = 1;
              bayesn_var_dptr = &BAYESN_MODEL_INFO.sigma0;
              break;
          }
          if (strcmp(event.data.scalar.value, "RV")==0)
          {
              datatype = 1;
              bayesn_var_dptr = &BAYESN_MODEL_INFO.RV;
              break;
          }
          if (strcmp(event.data.scalar.value, "TAUA")==0)
          {
              datatype = 1;
              bayesn_var_dptr = &BAYESN_MODEL_INFO.tauA;
              break;
          }
    
          // next we'll define the vectors
          if (strcmp(event.data.scalar.value, "L_KNOTS")==0)
          {
              datatype = 2;
              BAYESN_MODEL_INFO.lam_knots = malloc(sizeof(double)*(int)N_LAM);
              for(int i=0; i<(int)N_LAM; i++)
              {
                  BAYESN_MODEL_INFO.lam_knots[i] = 0.0;
              }
              bayesn_var_dptr = &BAYESN_MODEL_INFO.lam_knots[col];
              break;
          }
          if (strcmp(event.data.scalar.value, "TAU_KNOTS")==0)
          {
              datatype = 2;
              BAYESN_MODEL_INFO.tau_knots = malloc(sizeof(double)*(int)N_TAU);
              for(int i=0; i<(int)N_TAU; i++)
              {
                  BAYESN_MODEL_INFO.tau_knots[i] = 0.0;
              }
              bayesn_var_dptr = &BAYESN_MODEL_INFO.tau_knots[col];
              break;
          }
    
          // finally parse the 2D matrices
          if (strcmp(event.data.scalar.value, "W0")==0)
          {
              datatype = 3;
              rowsize = (int) N_TAU;
              W0 = malloc(sizeof(double)*(int)N_LAM*(int)N_TAU);
              row = 0;
              col = 0;
              bayesn_var_dptr = &W0[col];
              break;
          }
          if (strcmp(event.data.scalar.value, "W1")==0)
          {
              datatype = 3;
              rowsize = (int) N_TAU;
              W1 = malloc(sizeof(double)*(int)N_LAM*(int)N_TAU);
              row = 0;
              col = 0;
              bayesn_var_dptr = &W1[col];
              break;
          }
          if (strcmp(event.data.scalar.value, "L_SIGMA_EPSILON")==0)
          {
              datatype = 3;
              rowsize = (int) N_SIG;
              L_Sigma_epsilon = malloc(sizeof(double)*(int)N_SIG*(int)N_SIG);
              row = 0;
              col = 0;
              bayesn_var_dptr = &L_Sigma_epsilon[col];
              break;
          }
    
          if (datatype !=0 )
          {
              this_scalar = atof(event.data.scalar.value);
              //printf("Got scalar (value %f) %d %d\n", this_scalar, row, col); 
    
              // if we read a scalar we're done with one read
              if (datatype == 1)
              { 
                  // update the value of the variable that the pointer points to 
                  *bayesn_var_dptr = this_scalar;
                  datatype = 0; 
                  double this_scalar = 0.0;
                  bayesn_var_dptr = &this_scalar;
              }
              else
              {
                  *(bayesn_var_dptr + row*rowsize + col) = this_scalar;
                  this_scalar = 0.0;
                  col = col + 1;
              }
              break;
          }
      case YAML_SEQUENCE_START_EVENT: 
          // looks like you get a start sequence for either scalars or arrays
          // but you only get an end sequence for arrays
          // you can also get nested start sequences 
          break;
      
      case YAML_SEQUENCE_END_EVENT:   
          if (datatype == 2)
          {
              // if we have a vector, the first end sequence tells us we are done 
              col = 0;
              datatype = 0;
              double this_scalar = 0.0;
              bayesn_var_dptr = &this_scalar;
          }
          if (datatype == 3)
          {
              // if we have a matrix, each row ends with an end sequence
              // so move to the next row 
              row = row + 1;
              col = 0;
              // after the last row we finally get two end sequences
              if (last_event.type == YAML_SEQUENCE_END_EVENT)
              {
                  row = 0;
                  datatype = 0;
                  rowsize = 0;
                  double this_scalar = 0.0;
                  bayesn_var_dptr = &this_scalar;
              }
          }
          break;
      }
      // store the last event because we need to know when we have two end-sequence events in a row
      // to know when we are done reading a matrix
      last_event = event;
    
      if(event.type != YAML_STREAM_END_EVENT)
        yaml_event_delete(&event);
    } while(event.type != YAML_STREAM_END_EVENT);
    /* END new code */

    /* Cleanup */
    yaml_event_delete(&event);
    yaml_event_delete(&last_event);
    yaml_parser_delete(&parser);
    fclose(fh);

    // in principle we should just set this at read from YAML
    // but it's easier to read from YAML as double and fix here
    BAYESN_MODEL_INFO.n_lam_knots = (int) N_LAM;
    BAYESN_MODEL_INFO.n_tau_knots = (int) N_TAU;
    BAYESN_MODEL_INFO.n_sig_knots = (int) N_SIG;
    BAYESN_MODEL_INFO.W0 = gsl_matrix_alloc(BAYESN_MODEL_INFO.n_lam_knots, 
                                            BAYESN_MODEL_INFO.n_tau_knots);
    BAYESN_MODEL_INFO.W1 = gsl_matrix_alloc(BAYESN_MODEL_INFO.n_lam_knots, 
                                            BAYESN_MODEL_INFO.n_tau_knots);
    BAYESN_MODEL_INFO.L_Sigma_epsilon = gsl_matrix_alloc(BAYESN_MODEL_INFO.n_sig_knots, 
                                                         BAYESN_MODEL_INFO.n_sig_knots);

    // finally initalize the GSL matrices 
    int k = 0;
    for(int i=0; i< BAYESN_MODEL_INFO.n_lam_knots; i++)
    {
        for(int j=0; j< BAYESN_MODEL_INFO.n_tau_knots; j++)
        {
            k = i*BAYESN_MODEL_INFO.n_tau_knots + j;
            gsl_matrix_set(BAYESN_MODEL_INFO.W0, i, j, W0[k]);
            gsl_matrix_set(BAYESN_MODEL_INFO.W1, i, j, W1[k]);
        }
    }
    for(int i=0; i< BAYESN_MODEL_INFO.n_sig_knots; i++)
    {
        for(int j=0; j< BAYESN_MODEL_INFO.n_sig_knots; j++)
        {

            k = i*BAYESN_MODEL_INFO.n_sig_knots + j;
            gsl_matrix_set(BAYESN_MODEL_INFO.L_Sigma_epsilon, i, j, L_Sigma_epsilon[k]);
        }
    }
#endif

#ifndef USE_BAYESN
      sprintf(c1err,"genmag_BAYESN.o compiled without libyaml." );
      sprintf(c2err,"Install libyaml, set env YAML_DIR to a non-null string, make clean; make; try again.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
#endif
}

int init_genmag_BAYESN(char *version, int optmask){

    int  ised;
    int  retval = 0   ;
    int  ABORT_on_LAMRANGE_ERROR = 0;
    int  ABORT_on_BADVALUE_ERROR = 1;
    //char BANNER[120], tmpFile[200], sedcomment[40], version[60]  ;
    //
    char fnam[] = "init_genmag_BAYESN";

    // -------------- BEGIN --------------

    // extrac OPTMASK options

    // this loads all the BAYESN model components into the BAYESN_MODEL_INFO struct
    char *filename = "/global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SNDATA_ROOT/models/bayesn/BAYESN.M20/BAYESN.YAML";
    read_BAYESN_inputs(filename);

    /*
    // some code to print the W0/W1 matrices and make sure they are read correctly
    for(int i=0; i< BAYESN_MODEL_INFO.n_lam_knots; i++)
    {
        for(int j=0; j< BAYESN_MODEL_INFO.n_tau_knots; j++)
        {
            printf("(%d %d) %+.3f ",i, j, gsl_matrix_get(BAYESN_MODEL_INFO.W1, i, j));
        }
        printf("\n");
    }
    printf("\n");
    */

    char SED_filepath[] = "/global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SNDATA_ROOT/snsed/Hsiao07.dat";
    int istat;
    SEDMODEL_FLUX_DEF *S0 = &BAYESN_MODEL_INFO.S0;
    malloc_SEDFLUX_SEDMODEL(S0,0,0,0);
    double Trange[2] = {BAYESN_MODEL_INFO.tau_knots[0], 
                        BAYESN_MODEL_INFO.tau_knots[BAYESN_MODEL_INFO.n_tau_knots-1]};
    double Lrange[2] = {BAYESN_MODEL_INFO.lam_knots[0], 
                        BAYESN_MODEL_INFO.lam_knots[BAYESN_MODEL_INFO.n_lam_knots-1]};
    
    istat = rd_sedFlux(SED_filepath, "blah test rd_sedFlux", Trange, Lrange
	       ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL, 0
	       ,&S0->NDAY, S0->DAY, &S0->DAYSTEP
	       ,&S0->NLAM, S0->LAM, &S0->LAMSTEP
	       ,S0->FLUX,  S0->FLUXERR );

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

    gsl_matrix_free(K);
    gsl_matrix_free(D);

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
