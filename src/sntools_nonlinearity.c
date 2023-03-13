/* =====================================================

  Package to compute non-linearity.
  [initial use for WFIRST sims]

  Sim-input key  NONLINEARITY_FILE:  <inFile>
  read file of the form

  MODEL:  READOUT_RATE  ! obsolete xxxx
  MODEL:  FLUX_pe       ! May 27 2016

  START_MAP:
  FILTERS: abcdef
  NONLIN:   4.0  0.982  # 
  NONLIN:  40.0  0.986  # 
  NONLIN: 400.0  0.988  # 
  etc .
  END_MAP:

  Then repeat as many maps as needed in case
  different filters have different non-linearities.

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

  FILE *fp ;
  int  imap, RDFLAG, MAPSIZE, NLINE ;
  double tmpD[4];
  char c_get[100], MSG[100] ;
  char fnam[] = "INIT_NONLIN" ;
  
  // ------------ BEGIN --------------

  NMAP_NONLIN = 0 ;
  MODEL_NONLIN[0] = 0 ;
  NONLIN_README.NLINE = NLINE = 0 ;

  if ( IGNOREFILE(inFile) ) { return ; }

  fp = fopen(inFile,"rt");
  if ( !fp ) {
    sprintf(c1err,"Cannot open non-linearity file");
    sprintf(c2err,"%s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sprintf(BANNER,"%s : prepare non-linearity map \n", fnam);
  print_banner(BANNER);

  // first pass read to count how many tables for malloc
  // and to determine the model
  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"MODEL:"    ) == 0 ) { 
      readchar(fp,MODEL_NONLIN); 

      NONLIN_README.LINE[NLINE][0] = 0 ; NLINE++ ;
      sprintf(NONLIN_README.LINE[NLINE],
	      "  Apply NONLINEARITY MODEL '%s'", MODEL_NONLIN);
      NLINE++ ;
    }

    if ( strcmp(c_get,"START_MAP:") == 0 ) { NMAP_NONLIN++ ; }
  }

  rewind(fp);

  if ( strlen(MODEL_NONLIN) == 0 ) {
    sprintf(c1err,"Must specify  MODEL: <model>" );
    sprintf(c2err,"in %s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  NONLIN_MAP = (NONLIN_DEF*) malloc( NMAP_NONLIN * sizeof(NONLIN_DEF));

  // init all maps
  for(imap=0 ; imap < NMAP_NONLIN; imap++ ) {
    NONLIN_MAP[imap].FILTERS[0] = 0 ;
    NONLIN_MAP[imap].MAPSIZE  = 0 ;
  }
  DUMPFLAG_NONLIN = 0 ;

  // read again and store each map
  imap = RDFLAG = MAPSIZE = 0 ;
  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"START_MAP:") == 0 ) { RDFLAG = 1;  MAPSIZE=0 ; }
    
    if ( strcmp(c_get,"FILTERS:") == 0 ) 
      { readchar(fp, NONLIN_MAP[imap].FILTERS ); }

    if ( strcmp(c_get,"NONLIN:") == 0 ) {
      readdouble(fp, 2, tmpD );
      NONLIN_MAP[imap].MAPVAL[0][MAPSIZE]  = tmpD[0];
      NONLIN_MAP[imap].MAPVAL[1][MAPSIZE]  = tmpD[1];
      MAPSIZE++ ;
      NONLIN_MAP[imap].MAPSIZE = MAPSIZE ;
    }

    if ( strcmp(c_get,"DUMPFLAG:") == 0 ) 
      { readint(fp, 1, &DUMPFLAG_NONLIN ); }

    if ( strcmp(c_get,"END_MAP:") == 0 ) {       

      if ( MAPSIZE <= 0 ) {
	sprintf(c1err,"Found map with %d entries", MAPSIZE);
	sprintf(c2err,"Check %s", inFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      sprintf(MSG,"    Read NONLIN MAP%2.2d with %2d rows for %s",
	     imap, MAPSIZE, NONLIN_MAP[imap].FILTERS); 
      printf("%s\n", MSG);      fflush(stdout);
      imap++ ; RDFLAG=0;       

      sprintf(NONLIN_README.LINE[NLINE],"%s", MSG);
      NLINE++ ;
    }

  } // end while


  fclose(fp);

  printf("\n"); fflush(stdout);

  NONLIN_README.LINE[NLINE][0] = 0 ; NLINE++ ;
  NONLIN_README.NLINE = NLINE ;

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



