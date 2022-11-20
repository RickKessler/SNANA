#include <stdio.h>
#include <yaml.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
//#include "genmag_SEDtools.h"
//#include "genmag_BAYESN.h"

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


        if (strcmp(event.data.scalar.value, "W0")==0)
        {
            datatype = 3;
            break;
        }

        if (strcmp(event.data.scalar.value, "RV")==0)
        {
            datatype = 1;
            bayesn_var_dptr = &RV;
            break;
        }
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
        if (strcmp(event.data.scalar.value, "TAUA")==0)
        {
            datatype = 1;
            bayesn_var_dptr = &TAUA;
            break;
        }

        if (datatype !=0 )
        {
            this_scalar = atof(event.data.scalar.value);
             printf("Got scalar (value %f) %d %d\n", this_scalar, row, col); 

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


  /* Cleanup */
  yaml_parser_delete(&parser);
  fclose(fh);
  return 0;
}
