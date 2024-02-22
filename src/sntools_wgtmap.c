#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
//#include "sntools_cosmology.h"
//#include "snlc_sim.h"
//#include "sntools_host.h"
#include "sntools_wgtmap.h"
#include "sntools_output.h"

// ============================================
int read_VARNAMES_WGTMAP(char *WGTMAP_FILE, char *VARLIST_WGTMAP) {

  // ATG Reminder: Need to pass the fnam in as argument.
  // ATG Reminder: Refactor to call read_WGTMAP with flag to stop after 
  //               reading VARNAMES.
  // July 14 2020
  // pre-HOSTLIB-read utility to fetch & return comma-sep list of
  // VARNAMES appearing in WGTMAP so that WGTMAP variables can
  // be automatically stored in data files.
  // Check external WGTMAP file first; if no external WGTMAP file,
  // then check if WGTMAP is embedded in HOSTLIB.
  // Function returns number of WGTMAP variables found.
  //
  // Jan 15 2021: abort if EOF is reached.
  // Nov 18 2021: 
  //   + abort if WGT is not found
  //   + set FOUNDVAR_SNMAGSHIFT 
  //

  int NVAR   = 0 ;
  int MXCHAR = MXPATHLEN;
  int NWD, ivar, gzipFlag, NTMP=0 ;
  FILE *fp ;
  bool IS_SNVAR, FOUNDVAR_WGT=false, FOUNDVAR_SNMAGSHIFT=false;
  char FILENAME_FULL[MXPATHLEN], LINE[MXPATHLEN];
  char c_get[60], VARNAME[60];
  char KEY_VARNAMES[] = "VARNAMES_WGTMAP:" ;
  char KEY_STOP[]     = "GAL:" ; // stop reading when this key is found
  char fnam[]         = "read_VARNAMES_WGTMAP" ;
  int  LDMP = 0 ;
  // ------------- BEGIN ------------

  VARLIST_WGTMAP[0] = 0 ;

  // open file containing WGTMAP
  fp = open_TEXTgz(WGTMAP_FILE, "rt", &gzipFlag);
 
  bool STOP_READ = false;
  while( !STOP_READ ) { 

    if ( fscanf(fp, "%s", c_get) == EOF ) {
      sprintf(c1err,"Reached EOF before finding WGTMAP or GAL key.");
      sprintf(c2err,"Check format for HOSTLIB_FILE.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

    NTMP++;
    //    if ( NTMP < 200 ) 
    //{  printf(" xxx %s: c_get(%3d)='%s' \n", fnam, NTMP, c_get);  }

    // avoid reading the entire HOSTLIB file if there is no WGTMAP
    if ( strcmp(c_get,KEY_STOP) == 0 ) { STOP_READ = true; }

    if ( strcmp(c_get,KEY_VARNAMES) == 0 ) {
      STOP_READ = true ;
      fgets(LINE, MXCHAR, fp);
      NWD  = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE,fnam);
     
      for(ivar=0; ivar < NWD; ivar++ ) {
	get_PARSE_WORD(0,ivar,VARNAME);

	if ( strcmp(VARNAME,"WGT") == 0 )
	  { FOUNDVAR_WGT = true; continue; }

	//if ( strcmp(VARNAME,HOSTLIB_VARNAME_SNMAGSHIFT) == 0 )
	//  { FOUNDVAR_SNMAGSHIFT = true; continue; }

	NVAR++ ;
	/*
	IS_SNVAR = checkSNvar_HOSTLIB_WGTMAP(VARNAME);
	if ( !IS_SNVAR ) 
	  { catVarList_with_comma(VARLIST_WGTMAP,VARNAME); }
	if ( LDMP ) 
	  { printf(" xxx %s wgtmap var '%s'  (storeFlag=%d)\n", 
		   fnam, VARNAME, !IS_SNVAR ); } 
	*/
      }

    } // end reading VARNAMES_WGTMAP line
  } // end STOP_READ

  // - - - - -  
  /* ATG Might need this later?
  if ( NVAR > 0 && !FOUNDVAR_WGT ) {
    sprintf(c1err,"Missing required WGT column in WGTMAP");
    sprintf(c2err,"Check WGTMAP");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  HOSTLIB_WGTMAP.FOUNDVAR_SNMAGSHIFT = FOUNDVAR_SNMAGSHIFT;
  */

  if ( gzipFlag ){ pclose(fp); }     else { fclose(fp); }

  return(NVAR) ;

} // end read_VARNAMES_WGTMAP_LEGACY

