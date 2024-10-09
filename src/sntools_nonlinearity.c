/* =====================================================

  Package to read and apply non-linearity from map(s) in file.
  [initial use for WFIRST sims]

  Sim-input key  NONLINEARITY_FILE:  <inFile>
  where inFile is of the form:

  MODELNAME:  WHATEVER

  OPTMASK_NONLIN:  1  # 1=count tot, 2=count rate

  START_MAP:
  FILTERS: abcdef
  NONLIN:   4.0E0  0.982  #   Ftot(pe)  and Flux-scale
  NONLIN:   4.0E1  0.986
  NONLIN:   4.0E2  0.988 
  etc .
  END_MAP:

  Beware that Ftot should be the sum of Fsource + Fsky + Fgal.

  Repeat as many maps as needed in case different filters have 
  different non-linearities.

======================================================= */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_nonlinearity.h"


// =====================================
void INIT_NONLIN(char *inFile) {

  // Oct 2024:
  // Refactor to only require "FILTERS:" and "NONLIN:".
  // No longer need START_MAP and END_MAP keys.
  
  int langC = LANGFLAG_PARSE_WORDS_C;
  FILE *fp ;
  int  imap, nmap_read, MAPSIZE, NLINE, NWD ;
  double tmpD[4];
  bool  RDFLAG, END_OF_MAP, ISKEY_FILTERS, ISKEY_NONLIN ;
  char c_get[MXPATHLEN], LINE[MXPATHLEN], KEY[40], MSG[100] ;
  char fnam[] = "INIT_NONLIN" ;
  
  // ------------ BEGIN --------------

  NMAP_NONLIN = 0 ;
  MODELNAME_NONLIN[0] = 0 ;
  NONLIN_README.NLINE = NLINE = 0 ;
  OPTMASK_NONLIN = 0 ;
  
  if ( IGNOREFILE(inFile) ) { return ; }

  fp = fopen(inFile,"rt");
  if ( !fp ) {
    sprintf(c1err,"Cannot open non-linearity file");
    sprintf(c2err,"%s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sprintf(BANNER,"%s : prepare non-linearity map(s) in \n\t %s",
	  fnam, inFile);
  print_banner(BANNER);

  // first pass read to count how many tables for malloc
  // and to determine the model
  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"MODELNAME:"    ) == 0 ) { 
      readchar(fp,MODELNAME_NONLIN); 
      NONLIN_README.LINE[NLINE][0] = 0 ; NLINE++ ;
      sprintf(NONLIN_README.LINE[NLINE],
	      "  Apply NONLINEARITY MODEL '%s'", MODELNAME_NONLIN);
      NLINE++ ;
    }

    if ( strcmp(c_get,"OPTMASK_NONLIN:") == 0 ) 
      {  readint(fp, 1, &OPTMASK_NONLIN); }
    
    if ( strcmp(c_get,"FILTERS:") == 0 ) { NMAP_NONLIN++ ; }    
  }

  rewind(fp);

  if ( strlen(MODELNAME_NONLIN) == 0 ) {
    sprintf(c1err,"Must specify  MODELNAME: <model>" );
    sprintf(c2err,"in %s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  int NOPT_REQUIRE = 0;

  if ( (OPTMASK_NONLIN & OPTMASK_NONLIN_COUNT_TOT ) > 0 ) { NOPT_REQUIRE++; }
  if ( (OPTMASK_NONLIN & OPTMASK_NONLIN_COUNT_RATE) > 0 ) { NOPT_REQUIRE++; }  
  if ( NOPT_REQUIRE != 1 ) {
    sprintf(c1err,"Invalid NOPT_REQUIRE=%d; must require EITHER",
	    NOPT_REQUIRE);
    sprintf(c2err,"OPTMASK_NONLIN+=1(COUNT-TOTAL)  or +=2(COUNT-RATE) " );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // - - - -- 
  NONLIN_MAP = (NONLIN_DEF*) malloc( NMAP_NONLIN * sizeof(NONLIN_DEF));

  // init all maps
  for(imap=0 ; imap < NMAP_NONLIN; imap++ ) {
    NONLIN_MAP[imap].FILTERS[0] = 0 ;
    NONLIN_MAP[imap].MAPSIZE  = 0 ;
  }
  DUMPFLAG_NONLIN = 0 ;

  // - - - - - - -
  // read again and store each map
  nmap_read = MAPSIZE = RDFLAG = END_OF_MAP = 0 ;
  
  while ( fgets(LINE, MXPATHLEN, fp) != NULL ) {

    KEY[0] = 0 ;
    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING, LINE, fnam);

    if ( NWD >= 2 ) { get_PARSE_WORD(langC, 0, KEY); }
	    
    ISKEY_FILTERS = ( strcmp(KEY,"FILTERS:" ) == 0 );
    ISKEY_NONLIN  = ( strcmp(KEY,"NONLIN:"  ) == 0 );

    // check for end of map
    END_OF_MAP = RDFLAG && (NWD==0 || ISKEY_FILTERS );
    if ( END_OF_MAP ) {
      if ( MAPSIZE <= 0 ) {
	sprintf(c1err,"Found map with %d entries", MAPSIZE);
	sprintf(c2err,"Check %s", inFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      sprintf(MSG,"    Read NONLIN MAP%2.2d with %2d rows for %s",
	     imap, MAPSIZE, NONLIN_MAP[imap].FILTERS); 
      printf("%s\n", MSG);      fflush(stdout);
      sprintf(NONLIN_README.LINE[NLINE],"%s", MSG);
      NLINE++ ;
      RDFLAG = 0 ;
    }

    
    if ( ISKEY_FILTERS ) {
      MAPSIZE=0 ;  nmap_read++ ;  imap=nmap_read-1;  RDFLAG=1;
      get_PARSE_WORD(langC, 1, NONLIN_MAP[imap].FILTERS );
    }

    if ( ISKEY_NONLIN ) {
      get_PARSE_WORD_DBL(langC, 1, &tmpD[0] );
      get_PARSE_WORD_DBL(langC, 2, &tmpD[1] );      
      NONLIN_MAP[imap].MAPVAL[0][MAPSIZE]  = tmpD[0];
      NONLIN_MAP[imap].MAPVAL[1][MAPSIZE]  = tmpD[1];
      MAPSIZE++ ;
      NONLIN_MAP[imap].MAPSIZE = MAPSIZE ;
    }


  } // end while

  fclose(fp);

  printf("\n"); fflush(stdout);

  NONLIN_README.LINE[NLINE][0] = 0 ; NLINE++ ;
  NONLIN_README.NLINE = NLINE ;

  debugexit(fnam); // xxx REMOVE
  return ;

} // end INIT_NONLIN

void   init_nonlin__(char *inFile) { INIT_NONLIN(inFile); }


// ====================================
double GET_NONLIN(char *cfilt, double Fpe_source, double Fpe_sky, 
		  double mag ) {

  // Inputs
  //   + cfilt      = 1-char filter band
  //   + Fpe_source = the total source flux in photo-electrons.
  //   + Fpe_sky    = sky flux inside effective [NEA] aperture
  //   + mag        = true magnitude of object (not used yet)
  //
  // Returns F_meas/F_true .
  //
  // BEWARE: if nonLinearity is not defined for *cfilt, 
  //         function returns 1.0000
  //

  int    OPT_INTERP = 1;  // 1=linear interp
  int    imap ;
  double F_scale,  Fpe_tot;
  char  *ptrFilters, msg[100] ;
  char   fnam[] = "GET_NONLIN" ;

  // --------------- BEGIN ----------------

  F_scale = 1.000 ; // default is no non-linearity

  /*
  // xxxxxxxxxxxxx
  printf(" xxx NMAP=%d band=%s Fpe=%le  F_scale=%f\n",
	 NMAP_NONLIN, cfilt, Fpe, F_scale); fflush(stdout);
  // xxxxxxxxxxxxx
  */

  // bail if there are no maps.
  if ( NMAP_NONLIN == 0    ) { return(F_scale); }
  if ( Fpe_sky    <   0.0  ) { return(F_scale); }

  Fpe_tot = Fpe_source + Fpe_sky ;

  // find which map contains *cfilt.
  for(imap=0; imap < NMAP_NONLIN; imap++ ) {
    ptrFilters = NONLIN_MAP[imap].FILTERS ;

    if ( strstr(ptrFilters,cfilt) != NULL )  {
      sprintf(msg,"%s: band=%s imap=%d Fpe_tot=%le", 
	      fnam, cfilt, imap, Fpe_tot);


      F_scale = interp_1DFUN(OPT_INTERP, Fpe_tot,
			     NONLIN_MAP[imap].MAPSIZE,
			     NONLIN_MAP[imap].MAPVAL[0], // Ftot(pe)
			     NONLIN_MAP[imap].MAPVAL[1], // F_scale
			     msg ) ;      
    }
  }  // end imap loop


  int LDMP = DUMPFLAG_NONLIN ;
  if ( LDMP ) {
    printf(" xxx %s-mag=%6.3f  Fpe(sky,src)=%10.5le,%10.5le "
	   "--> F_scale=%6.4f \n",
	   cfilt, mag, Fpe_sky, Fpe_source, F_scale); 
    fflush(stdout);
  }

  return(F_scale);

} // end GET_NONLIN


double get_nonlin__(char *cfilt, double *Fpe_source, double *Fpe_sky,
		    double *genmag) {
  
  double F_scale = GET_NONLIN(cfilt, *Fpe_source, *Fpe_sky, *genmag);
  return(F_scale);
}



