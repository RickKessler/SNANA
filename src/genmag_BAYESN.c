/**********************************************

  July 1 2022 G. Narayan

********************************************/

#include "stdio.h"
#include "stddef.h"
#include "math.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_matrix.h"
#include "fitsio.h"
#include "sntools.h"
#include "genmag_SEDtools.h"
#include  "genmag_BAYESN.h"
#include "MWgaldust.h"

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
    int i;

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
              for(i=0; i<(int)N_LAM; i++)
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
              for(i=0; i<(int)N_TAU; i++)
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
    
              // if we read a scalar we're done with one read
              if (datatype == 1)
              { 
                  if (VERBOSE_BAYESN > 0)
                  {
                    printf("DEBUG: Read this scalar from yaml %.2f\n", this_scalar);
                  }
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
    int k = 0, j;
    for(i=0; i< BAYESN_MODEL_INFO.n_lam_knots; i++)
    {
        for(j=0; j< BAYESN_MODEL_INFO.n_tau_knots; j++)
        {
            k = i*BAYESN_MODEL_INFO.n_tau_knots + j;
            gsl_matrix_set(BAYESN_MODEL_INFO.W0, i, j, W0[k]);
            gsl_matrix_set(BAYESN_MODEL_INFO.W1, i, j, W1[k]);
        }
    }
    for(i=0; i< BAYESN_MODEL_INFO.n_sig_knots; i++)
    {
        for(j=0; j< BAYESN_MODEL_INFO.n_sig_knots; j++)
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

int init_genmag_BAYESN(char *MODEL_VERSION, int optmask){

    int  ised;
    int  retval = 0   ;
    int  ABORT_on_LAMRANGE_ERROR = 0;
    int  ABORT_on_BADVALUE_ERROR = 1;
    //char BANNER[120], tmpFile[200], sedcomment[40], version[60]  ;
    char fnam[] = "init_genmag_BAYESN";

    // -------------- BEGIN --------------

    // extrac OPTMASK options
    VERBOSE_BAYESN = optmask & OPTMASK_BAYESN_VERBOSE;
    if (VERBOSE_BAYESN)
    {
        printf("DEBUG: VERBOSE_BAYESN flag is set\n");
        printf("DEBUG: %s MODEL_VERSION=%s", fnam, MODEL_VERSION);
    }


    // GN - there is one of these further below but unsure yet how RK wants us to use these

    // this loads all the BAYESN model components into the BAYESN_MODEL_INFO struct
    // HACK HACK HACK
    char version[60];
    extract_MODELNAME(MODEL_VERSION, BAYESN_MODELPATH, version);
    char yaml_file[MXPATHLEN];
    sprintf(yaml_file, "%s/BAYESN.YAML", BAYESN_MODELPATH);
    read_BAYESN_inputs(yaml_file);

    // HACK HACK HACK 
    char SED_filepath[] = "/global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SNDATA_ROOT/snsed/Hsiao07.dat";
    int istat;
    SEDMODEL_FLUX_DEF *S0 = &BAYESN_MODEL_INFO.S0;
    malloc_SEDFLUX_SEDMODEL(S0,0,0,0);
    double Trange[2] = {BAYESN_MODEL_INFO.tau_knots[0], 
                        BAYESN_MODEL_INFO.tau_knots[BAYESN_MODEL_INFO.n_tau_knots-1]};
    double Lrange[2] = {BAYESN_MODEL_INFO.lam_knots[0], 
                        BAYESN_MODEL_INFO.lam_knots[BAYESN_MODEL_INFO.n_lam_knots-1]};
    
    istat = rd_sedFlux(SED_filepath, "Hsiao Template", Trange, Lrange
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

    SEDMODEL.LAMMIN_ALL = BAYESN_MODEL_INFO.lam_knots[0] ;  // rest-frame SED range
    SEDMODEL.LAMMAX_ALL = BAYESN_MODEL_INFO.lam_knots[BAYESN_MODEL_INFO.n_lam_knots-1] ;
    SEDMODEL.RESTLAMMIN_FILTERCEN =  SEDMODEL.LAMMIN_ALL + 1200.0 ; // rest-frame central wavelength range
    SEDMODEL.RESTLAMMAX_FILTERCEN =  SEDMODEL.LAMMAX_ALL - 1200.0 ;

    if (VERBOSE_BAYESN > 0)
    {
        printf("DEBUG: LIMITS OF FILTER CENTRAL WAVELENGTH: %.1f, %.1f\n", SEDMODEL.RESTLAMMIN_FILTERCEN
            , SEDMODEL.RESTLAMMAX_FILTERCEN);
        printf("DEBUG: LIMITS OF SED WAVELENGTH: %.1f, %.1f\n", SEDMODEL.LAMMIN_ALL, SEDMODEL.LAMMAX_ALL);
    }

    //compute the inverse KD matrices and J_lam (SHOULD THIS BE DONE HERE??)
    BAYESN_MODEL_INFO.KD_tau = invKD_irr(BAYESN_MODEL_INFO.n_tau_knots,
            BAYESN_MODEL_INFO.tau_knots);
    BAYESN_MODEL_INFO.KD_lam = invKD_irr(BAYESN_MODEL_INFO.n_lam_knots,
            BAYESN_MODEL_INFO.lam_knots);
    BAYESN_MODEL_INFO.J_lam = spline_coeffs_irr(BAYESN_MODEL_INFO.S0.NLAM,
            BAYESN_MODEL_INFO.n_lam_knots, BAYESN_MODEL_INFO.S0.LAM,
            BAYESN_MODEL_INFO.lam_knots, BAYESN_MODEL_INFO.KD_lam);

    //debugexit(fnam);
    fflush(stdout);
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

    //bool dumpsed = false;
    bool enable_scatter = false;

    /*
    FILE * sedfile;
    if (dumpsed) {
        sedfile = fopen("sed_dump.txt", "w");
    }*/
    
    double DLMAG = parList_SN[0];
    double THETA = parList_SN[1];
    double AV    = parList_SN[2];
    double RV    = parList_SN[3];
    double z1, meanlam_obs,  meanlam_rest, ZP; 
    char *cfilt;
    int ifilt = 0, i;
    
    //SHOULD I BE DECLARING THESE HERE??
    gsl_matrix * J_tau; // for time interpolation
    gsl_matrix * W = gsl_matrix_alloc(BAYESN_MODEL_INFO.n_lam_knots,
            BAYESN_MODEL_INFO.n_tau_knots); // for W0 + THETA*W1
    gsl_matrix * WJ_tau = gsl_matrix_alloc(BAYESN_MODEL_INFO.n_lam_knots,
            Nobs); //to store matrix product W * J_tau
    gsl_vector_view j_lam; //to store a row of J_lam
    gsl_vector * jWJ = gsl_vector_alloc(Nobs); 

    double *lam_filt;
    double *trans_filt;
    double *lam_model;

    double mag;
    char fnam[] = "genmag_BAYESN";

    // ------- BEGIN -----------
    // translate absolute filter index into sparse index
    ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
    z1    = 1. + z ;


    // HACK HACK HACK - GN - why do the bloody phases not match by 1/1+z??? 20230210
    double *Trest_list   = malloc(sizeof(double)*Nobs);
    for(i=0;i<Nobs;i++)
    {
        Trest_list[i] = Tobs_list[i]/z1;
    }

    // filter info for this "ifilt"
    meanlam_obs  = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
    ZP           = FILTER_SEDMODEL[ifilt].ZP ;
    cfilt        = FILTER_SEDMODEL[ifilt].name ;
    meanlam_rest = meanlam_obs/z1 ;
    // make sure filter-lambda range is valid
    checkLamRange_SEDMODEL(ifilt,z,fnam);

    // get the filter wavelengths
    int nlam_filter = FILTER_SEDMODEL[ifilt].NLAM;
    lam_filt   = FILTER_SEDMODEL[ifilt].lam;
    trans_filt = FILTER_SEDMODEL[ifilt].transSN;

    // get the hsiao wavelengths
    int nlam_model = BAYESN_MODEL_INFO.S0.NLAM;
    lam_model      = BAYESN_MODEL_INFO.S0.LAM;
    double d_lam   = lam_model[1] - lam_model[0];

    // project the model into the observer frame 
    // get the rest-frame model wavelengths that overlap with the filter 
    int ilam_blue = 0;
    while (z1*lam_model[ilam_blue] <= lam_filt[0]) {
        ilam_blue++;
    }
    int ilam_red = nlam_model - 1;
    while (z1*lam_model[ilam_red] >= lam_filt[FILTER_SEDMODEL[ifilt].NLAM-1]) {
        ilam_red--;
    }

    // compute the matrix for time interpolation
    J_tau = spline_coeffs_irr(Nobs, BAYESN_MODEL_INFO.n_tau_knots,
            Trest_list, BAYESN_MODEL_INFO.tau_knots, BAYESN_MODEL_INFO.KD_tau);

    // compute W0 + theta*W1 (SHOULD THIS BE DONE HERE??)
    // also, this seems utterly unhinged as a way of computing this
    gsl_matrix_set_zero(W);
    gsl_matrix_add(W, BAYESN_MODEL_INFO.W1);
    gsl_matrix_scale(W, THETA);
    gsl_matrix_add(W, BAYESN_MODEL_INFO.W0);

    /* int wx, wy;
    printf("DEBUG: W matrix\n");
    for(wx=0; wx< 9; wx++)
    {
        for(wy=0; wy < 6; wy++)
        {
            printf("%.3f  ",gsl_matrix_get(W, wx, wy));
        }
        printf("\n");
    }
    printf("-----------\n"); */


    // compute W * J_tau^T
    gsl_matrix_set_zero(WJ_tau);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, W, J_tau, 0.0, WJ_tau);

    // interpolate the filter wavelengths on to the model in the observer frame
    // usually this is OK because the filters are more coarsely defined than the model
    // that may not be the case with future surveys and we should revisit
    int    this_nlam = ilam_red - ilam_blue + 1;
    int    o, q;
    double this_lam;
    double this_trans;
    double eA_lam_MW, eA_lam_host; //To store MW and host dust law evaluated at current wl
    double eW, S0_lam; //To store other SED  bits
    for (o = 0; o < Nobs; o++) { magobs_list[o] = 0.0; } //Set magnitudes to 0
    for(q=ilam_blue; q<ilam_red; q++)
      {
        this_lam   = lam_model[q]*z1;
        this_trans = interp_1DFUN(2, this_lam, nlam_filter, 
				       lam_filt, trans_filt, "DIE");

        // super weird computation
        // this finds a vector of length Nobs, giving the SED at the
        // current wavelength for all observations
        // basically j_lam * W * J_tau^T
        //
        // GN - 20230203 - Why the setting jWJ to zero here - we set it earlier
        // See enable_scatter
        //
        gsl_vector_set_zero(jWJ);
        j_lam = gsl_matrix_row(BAYESN_MODEL_INFO.J_lam, q);
        
        // GSN - 20230602 - J_lam matches - compare in restframe
        //printf("DEBUG lam_model: %.2f     lam filt: %.2f     j_lam: %.5f\n",lam_model[q], this_lam, gsl_vector_get(&j_lam.vector, 0));
        gsl_blas_dgemv(CblasTrans, 1.0, WJ_tau, &j_lam.vector, 0.0, jWJ);

        //GSN - 20230617 - get the right extinction
        eA_lam_MW = pow(10.0, -0.4*GALextinct(3.1, 3.1*mwebv, this_lam, 99));
        eA_lam_host = pow(10.0, -0.4*GALextinct(RV, AV, lam_model[q], 99));
        /*if (VERBOSE_BAYESN > 0)
        {
            printf("DEBUG: eA_lam_MW: %.3f     eA_lam_host: %.3f\n", eA_lam_MW, eA_lam_host);
        }*/

        int q_hsiao;
        for (o = 0; o < Nobs; o++) {
            /*if (o < 20){
                printf("DEBUG: Trest: %.2f    JWJ:   %.5f\n", Trest_list[o], gsl_vector_get(jWJ, o));
            }*/
            eW = pow(10.0, -0.4*gsl_vector_get(jWJ, o));
            // Seek the first Hsiao timestep above the current obs time
            q_hsiao = 0;
            while (BAYESN_MODEL_INFO.S0.DAY[q_hsiao] <= Trest_list[o]) { q_hsiao++; }
            if (q_hsiao < 0 || q_hsiao >= BAYESN_MODEL_INFO.S0.NDAY) {
                sprintf(c1err,"Time index outside of template range." );
                sprintf(c2err,"Invalid q_hsiao index %d. Valid range is [%d, %d]", q_hsiao, 0, BAYESN_MODEL_INFO.S0.NDAY);
                errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
            }
            double t1 = BAYESN_MODEL_INFO.S0.DAY[q_hsiao];
            double t0 = BAYESN_MODEL_INFO.S0.DAY[q_hsiao-1];
            double f1 = BAYESN_MODEL_INFO.S0.FLUX[nlam_model*q_hsiao + q];
            double f0 = BAYESN_MODEL_INFO.S0.FLUX[nlam_model*(q_hsiao-1) + q];
            // HACK HACK HACK throw Trest at this instead or Tobs
            S0_lam = (f0*(t1 - Trest_list[o]) + f1*(Trest_list[o] - t0))/(t1 - t0);
            /*if (o == 0 && q == ilam_blue + 10) {
                printf("XXX DEBUG Trest %.3f; q %d; this_trans %.6f; this_lam %.3f; d_lam %.3f; eA_lam_MW %.3f; eA_lam_host %.3f; eW %.3f; S0_lam %le\n", Trest_list[o], q, this_trans, this_lam, d_lam, eA_lam_MW, eA_lam_host, eW, S0_lam);
            }*/
            magobs_list[o] += this_trans*this_lam*d_lam*eA_lam_MW*eA_lam_host*eW*S0_lam; //Increment flux with contribution from this wl

            /*if (o == 0) {
                dump_sed_element(sedfile, lam_model[q], eA_lam_MW*eA_lam_host*eW*S0_lam);
            }*/
        }
    }

    /*if (dumpsed) {
        fclose(sedfile);
    }*/

    // GSN - 20230602 - J_tau matrix matches but JWJ does not - J_lam matches (look at rest-wavelengths)
    /*if (VERBOSE_BAYESN > 0)
    {
        printf("DEBUG: Printing J_tau matrix\n");
	    int crap = print_matrix(stdout, J_tau);
        printf("-----\n\n\n");
    }*/


    gsl_matrix_free(J_tau);
    gsl_matrix_free(W);
    gsl_matrix_free(WJ_tau);
    gsl_vector_free(jWJ);

    double zdum = 2.5*log10(1.0+z);
    if (VERBOSE_BAYESN > 0)
    {
        printf("DEBUG: BAYESN_MODEL_INFO.M0: %.2f   DLMAG: %.2f   ZP: %.2f   THETA: %.2f   AV: %.2f\n", BAYESN_MODEL_INFO.M0, DLMAG, ZP, THETA, AV);
    }
    double hc_local = hc;
    for (o = 0; o < Nobs; o++) {
      magobs_list[o] = BAYESN_MODEL_INFO.M0 + DLMAG 
	-2.5*log10(magobs_list[o]/hc_local) + ZP; //WE SHOULD CHECK THE ZERO POINTS

      magerr_list[o] = 0.1;
    }
    if (VERBOSE_BAYESN > 0)
    {
        printf("DEBUG: Printing lightcurve\n");
        for (o = 0; o < Nobs; o++) 
        { 
            printf("%.2f  ",magobs_list[o]);
        }
        printf("\n-----\n\n\n");
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

void dump_SED_element(FILE * file, double wave, double value) {
   fprintf(file, "%.2f %e", wave, value);
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
			while (q < Nk-2 && xk[q+1] <= x[i]) { q++; }
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

int print_matrix(FILE *f, const gsl_matrix *m)
{
        int status, n = 0;
        size_t i=0;
        size_t j=0;

        for (i = 0; i < m->size1; i++) {
                for (j = 0; j < m->size2; j++) {
                        if ((status = fprintf(f, "%g ", gsl_matrix_get(m, i, j))) < 0)
                                return -1;
                        n += status;
                }

                if ((status = fprintf(f, "\n")) < 0)
                        return -1;
                n += status;
        }

        return n;
}

