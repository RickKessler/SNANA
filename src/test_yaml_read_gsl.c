#include <stdio.h>
#include <yaml.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
//#include "genmag_SEDtools.h"
//#include "genmag_BAYESN.h"



float malloc_double2D(int opt, int LEN1, int LEN2, double ***array2D ) {
  // Created Jun 11 2019
  // Malloc array2D[LEN1][LEN2]  (intended for LEN1=NSN, LEN2=NCLPAR)
  float f_MEMTOT = 0.0 ;
  long long MEMTOT=0, i1 ;
  int MEM1 = LEN1 * sizeof(double *); 
  int MEM2 = LEN2 * sizeof(double);
  // ----------- BEGIN -------------

  if ( opt > 0 ) {

    *array2D = (double**) malloc(MEM1) ; MEMTOT += MEM1;
    for(i1=0; i1< LEN1; i1++ ) {
      (*array2D)[i1] = (double*) malloc(MEM2) ; MEMTOT += MEM2;
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1 < LEN1; i1++ ) { free((*array2D)[i1]); }
    free(array2D[i1]) ;    
  }
  return(f_MEMTOT);
}


int main(void)
{
  FILE *fh = fopen("/global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SNDATA_ROOT/models/bayesn/BAYESN.M20/BAYESN.YAML", "r");
  yaml_parser_t parser;
  yaml_event_t  event;  
  yaml_event_t last_event;

  int datatype = 0; // 0: we don't know 1: scalar 2: vector 3: matrix
                    // in principle we can support strings and other types
  int col=0;
  int row=0;
  int rowsize = 0;
  // default init for the sizes is negative to deliberately force an error
  // we should populate these values from the YAML file
  // they are also doubles because it is easier to assume everything is a double when parsing
  // and cast to int when needed
  double N_LAM =  -1.0;
  double N_TAU =  -1.0;
  double N_SIG =  -1.0;
  double RV    =   3.1;
  double M0    = -19.0;
  double SIGMA0=   0.1;
  double TAUA  =   0.4;
  double *LAM_KNOTS;
  double *TAU_KNOTS;
  gsl_matrix *L_SIGMA_EPSILON;
  gsl_matrix *W0;
  gsl_matrix *W1;
  gsl_matrix *bayesn_var_matrix;

  //double **LAM_SIGMA_EPSILON;
  //double **W0;
  //double **W1;

  // we need something to store the current scalar
  double this_scalar = 0.0;
  // and something to point to the current BayeSN variable being populated
  // this only works if datatype is in 1-3 (assumed double)
  double *bayesn_var_dptr = &this_scalar;

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
    case YAML_STREAM_START_EVENT: puts("STREAM START"); break;
    case YAML_STREAM_END_EVENT:   puts("STREAM END");   break;

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
            bayesn_var_dptr = &M0;
            break;
        }
        if (strcmp(event.data.scalar.value, "SIGMA0")==0)
        {
            datatype = 1;
            bayesn_var_dptr = &SIGMA0;
            break;
        }
        if (strcmp(event.data.scalar.value, "RV")==0)
        {
            datatype = 1;
            bayesn_var_dptr = &RV;
            break;
        }
        if (strcmp(event.data.scalar.value, "TAUA")==0)
        {
            datatype = 1;
            bayesn_var_dptr = &TAUA;
            break;
        }

        // next we'll define the vectors
        if (strcmp(event.data.scalar.value, "L_KNOTS")==0)
        {
            datatype = 2;
            LAM_KNOTS = malloc(sizeof(double)*(int)N_LAM);
            for(int i=0; i<(int)N_LAM; i++)
            {
                LAM_KNOTS[i] = 0.0;
            }
            bayesn_var_dptr = &LAM_KNOTS[col];
            break;
        }
        if (strcmp(event.data.scalar.value, "TAU_KNOTS")==0)
        {
            datatype = 2;
            TAU_KNOTS = malloc(sizeof(double)*(int)N_TAU);
            for(int i=0; i<(int)N_TAU; i++)
            {
                TAU_KNOTS[i] = 0.0;
            }
            bayesn_var_dptr = &TAU_KNOTS[col];
            break;
        }

        // finally parse the 2D matrices
        if (strcmp(event.data.scalar.value, "W0")==0)
        {
            datatype = 3;
            rowsize = (int) N_TAU;
            // malloc_double2D(1, (int) N_LAM, (int) N_TAU, &W0);

            row = 0;
            col = 0;
            bayesn_var_matrix = W0;
            //bayesn_var_dptr = &W0[row][col];
            gsl_matrix_set_zero(W0);
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
                if (datatype == 2)
                {
                    *(bayesn_var_dptr + col) = this_scalar;
                }
                else
                {
                    gsl_matrix_set(bayesn_var_matrix, row, col, this_scalar);
                }
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
  yaml_event_delete(&event);
  /* END new code */

  printf("Vars NLAM %d NTAU %d NSIG %d M0 %f SIGMA0 %f RV %f TAUA %f\n",
            (int)N_LAM, (int)N_TAU, (int)N_SIG, 
            M0, SIGMA0, RV, TAUA);
  printf("LAM_KNOTS:\n");
  for(int i=0; i<(int)N_LAM; i++)
  {
    printf("%f ",LAM_KNOTS[i]);
  }
  printf("\n");
  printf("TAU_KNOTS:\n");
  for(int i=0; i<(int)N_TAU; i++)
  {
    printf("%f ",TAU_KNOTS[i]);
  }
  printf("\n");
  
  printf("W0:\n");
  rowsize=(int) N_TAU;
  for(int i=0; i<(int)N_LAM; i++)
  {
    for(int j=0; j<(int)N_TAU; j++)
    {
        printf("(%d %d) %.2f",i,j, gsl_matrix_get(W0, i, j));
    }
    printf("\n");
  }
  printf("\n");

  /* Cleanup */
  yaml_parser_delete(&parser);
  fclose(fh);
  return 0;
}
