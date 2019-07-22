/*****************************************************

 Created Feb 12 2018

  Begin new version of fluxErrModels to treat all
  maps the same way. 

  This code will eventually replace the legacy codes in 
  sntools_fluxErrModels.c[h].

  Mar 16 2019: 
    refactor INIT_FLUXERRMODEL to use read_GRIDMAP, and to read OPT_EXTRAP.

****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_fluxErrModels.h"


// ===========================================================
void INIT_FLUXERRMODEL(int OPTMASK, char *fileName, char *MAPLIST_IGNORE_DATAERR) {

  // Created Feb 2018
  // Read and store maps to correct flux-uncertainties.
  // Corrections can be additive or scale, and depend on
  // and combination of 
  //    BAND, FIELD, MJD, SKYSIG, PSF, ZP, SNR, SB
  //
  // Inputs:
  //   OPTMASK  : bit options
  //   fileName : file containing maps to read
  //   MAPLIST_IGNORE_DATAERR : 
  //      optional [comma-separated] list of maps to use only for simulation,
  //      but not include in the reported data error.
  //                            
  //
  // Mar 16 2019:
  //  + refactor to use read_GRIDMAP
  //  + NVAR key is obsolete (see warn_MNVAR_KEY)
  //  + read OPT_EXTRAP key from map file, and pass to read_GRIDMAP().
  //

  FILE *fp;
  int gzipFlag, FOUNDMAP, NTMP, NVAR, NDIM, NFUN, ivar, igroup;
  int IDMAP, NMAP=0, imap, OPT_EXTRAP=0 ; 
  char PATH[MXPATHLEN], c_get[80];  
  char *fullName = FILENAME_FLUXERRMAP ;
  char *name, *fieldList, TMP_STRING[80], LINE[100];
  char MSGERR_FILE[200];
  char fnam[] = "INIT_FLUXERRMODEL" ;

  // ------------ BEGIN --------------

  NMAP_FLUXERRMODEL          = 0 ;
  FLUXERR_FIELDGROUP.NDEFINE = 0 ;  
  NINDEX_SPARSE_FLUXERRMAP   = 0 ;
  if ( IGNOREFILE(fileName) ) { return ; }

  sprintf(BANNER,"%s:", fnam);
  print_banner(BANNER);

  sprintf(PATH, "%s/simlib", PATH_SNDATA_ROOT);
  fp = snana_openTextFile(1,PATH, fileName, fullName, &gzipFlag);

  sprintf(MSGERR_FILE,"check FLUXERRMODEL_FILE: '%s'", fullName);
  FOUNDMAP=0;

  // function inputs
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_MJD],     "MJD");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_PSF],     "PSF");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_SKYSIG],  "SKYSIG");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_ZP],      "ZP");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_LOGSNR],  "LOGSNR");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_SBMAG],   "SBMAG");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_GALMAG],  "GALMAG");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_SNSEP],   "SNSEP");

  // function values
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_ERRSCALE], "ERRSCALE");
  sprintf(VARNAMES_FLUXERRMAP[IPAR_FLUXERRMAP_ERRADD],   "ERRADD");


  // --- init extrapolation flag ---
  FLUXERRMAP_EXTRAP.FLAG   = 1 ;  // default is edge extrapolation
  for(imap=0; imap < MXMAP_FLUXERRMAP ; imap++ ) {
    FLUXERRMAP_EXTRAP.NOBS_TOT[imap]    = 0 ;
    FLUXERRMAP_EXTRAP.NOBS_EXTRAP_LO[imap] = 0 ;
    FLUXERRMAP_EXTRAP.NOBS_EXTRAP_HI[imap] = 0 ;
  }

  // start reading file
  while ( fscanf(fp, "%s", c_get) != EOF ) {

    if ( strcmp(c_get,"OPT_EXTRAP:")==0 ) { readint(fp,1,&OPT_EXTRAP); }

    if ( strcmp(c_get,"DEFINE_FIELDGROUP:")==0 ) {
      NTMP      = FLUXERR_FIELDGROUP.NDEFINE ;
      name      = FLUXERR_FIELDGROUP.NAME[NTMP] ;
      fieldList = FLUXERR_FIELDGROUP.FIELDLIST[NTMP] ;
      readchar(fp, name );
      readchar(fp, fieldList );
      printf("\t FIELD-assign %s --> %s \n", name, fieldList );
      FLUXERR_FIELDGROUP.NDEFINE++ ;
    }

    if ( strcmp(c_get,"MAPNAME:")==0 ) {
      FOUNDMAP = 1 ;
      NMAP = NMAP_FLUXERRMODEL; 
      name = FLUXERRMAP[NMAP].NAME ;
      readchar(fp, name ) ;
      FLUXERRMAP[NMAP].NVAR = NVAR  = 0 ;
      FLUXERRMAP[NMAP].MAP.VARLIST[0]   = 0 ;
      FLUXERRMAP[NMAP].MASK_APPLY   = 3;    // SIM & DATA by default 
      
      sprintf(FLUXERRMAP[NMAP].BANDLIST,"ALL"); // default is all bands
      sprintf(FLUXERRMAP[NMAP].FIELDLIST,"ALL"); // default is all fields
      FLUXERRMAP[NMAP].INDEX_SPARSE = index_sparse_FLUXERRMAP(NMAP,name) ;

      NMAP_FLUXERRMODEL++ ; 
    }

    if ( strcmp(c_get,"BAND:")==0 || strcmp(c_get,"FILTER:")==0  ) {
      if ( FOUNDMAP == 0 ) {
	sprintf(c1err,"%s key not allowed outside MAP", c_get);
	errmsg(SEV_FATAL, 0, fnam, c1err, MSGERR_FILE ); 
      }
      readchar(fp, FLUXERRMAP[NMAP].BANDLIST );
    } // end BAND

    if ( strcmp(c_get,"FIELD:")==0 ) {
      if ( FOUNDMAP == 0 ) {
	sprintf(c1err,"%s key not allowed outside MAP", c_get);
	errmsg(SEV_FATAL, 0, fnam, c1err, MSGERR_FILE ); 
      }
      readchar(fp, TMP_STRING);
      sprintf(FLUXERRMAP[NMAP].FIELDLIST,      "%s", TMP_STRING);
      sprintf(FLUXERRMAP[NMAP].FIELDLIST_ORIG, "%s", TMP_STRING);
      for (igroup=0; igroup < FLUXERR_FIELDGROUP.NDEFINE; igroup++ ) {
	name      = FLUXERR_FIELDGROUP.NAME[igroup] ;
	fieldList = FLUXERR_FIELDGROUP.FIELDLIST[igroup] ;
	if ( strcmp(TMP_STRING,name) == 0 ) 
	  { sprintf(FLUXERRMAP[NMAP].FIELDLIST, "%s", fieldList); }
      }
    }  // end FIELD

    if ( strcmp(c_get,"NVAR:")==0 ) { warn_NVAR_KEY(fullName); }

    if ( strcmp(c_get,"VARNAMES:")==0 ) {
      fgets(LINE,100,fp);
      NVAR = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
      FLUXERRMAP[NMAP].NVAR = NVAR ;
      //      malloc_ROWDATA_FLUXERRMAP(+1,NVAR);

      for(ivar=0; ivar < NVAR; ivar++ ) { 
	get_PARSE_WORD(0,ivar,TMP_STRING);
	FLUXERRMAP[NMAP].IVARLIST[ivar] = IVARLIST_FLUXERRMAP(TMP_STRING);
	sprintf(FLUXERRMAP[NMAP].VARNAMES[ivar],"%s", TMP_STRING);
	if ( ivar>0 && ivar < NVAR-1) 
	  { strcat(FLUXERRMAP[NMAP].MAP.VARLIST,","); }
	if ( ivar < NVAR-1 ) 
	  { strcat(FLUXERRMAP[NMAP].MAP.VARLIST,TMP_STRING); }
      }

      IDMAP = IDGRIDMAP_FLUXERRMODEL_OFFSET + NMAP ;
      NDIM  = NVAR-1;  NFUN=1;
      read_GRIDMAP(fp,"ROW:", "ENDMAP:", IDMAP, NDIM, NFUN, OPT_EXTRAP,
		   MXROW_FLUXERRMAP, fnam, 
		   &FLUXERRMAP[NMAP].MAP );  // <== returned		   

      if ( (OPTMASK & MASK_DUMP_FLUXERRMAP)>0 )
	{ DUMP_FLUXERRMAP(NMAP); }

      FOUNDMAP=0 ;

    }  // end VARNAMES

    /* xxxxxxxxxxxxxxx mark delete Mar 16 2019 xxxxxxxx
    if ( strcmp(c_get,"ROW:")==0 ) {
      irow = FLUXERRMAP[NMAP].NROW;
      readdouble( fp, NVAR, TMPVAL );
      if ( irow < MXROW_FLUXERRMAP ) {
	for(ivar=0; ivar < NVAR; ivar++ ) 
	  { TMP_ROWDATA_FLUXERRMAP[ivar][irow] = TMPVAL[ivar]; }
      }
      FLUXERRMAP[NMAP].NROW++ ;
    }
    xxxxxxxxxx end mark xxxxxxxx*/


    /* xxxxxxxxxxxxxx mark delete xxxxxxxxxxxxxxxxxx
    if ( strcmp(c_get,"ENDMAP:")==0 ) {

      NROW = FLUXERRMAP[NMAP].NROW ;
      if ( NROW >= MXROW_FLUXERRMAP ) {
	sprintf(c1err,"NROW=%d exceeds bound (MXROW=%d)",
		NROW, MXROW_FLUXERRMAP);
	sprintf(c2err,"Check MAPNAME='%s'  BAND='%s'  FIELD='%s' "
		,FLUXERRMAP[NMAP].NAME 
		,FLUXERRMAP[NMAP].BANDLIST
		,FLUXERRMAP[NMAP].FIELDLIST );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      sprintf(TMP_STRING,"%s(%s-%s)",
	      FLUXERRMAP[NMAP].NAME,
	      FLUXERRMAP[NMAP].FIELDLIST,
	      FLUXERRMAP[NMAP].BANDLIST );

      IDMAP = IDGRIDMAP_FLUXERRMODEL_OFFSET + NMAP ;
      NROW  = FLUXERRMAP[NMAP].NROW ;
      
      init_interp_GRIDMAP(IDMAP, TMP_STRING, NROW, NVAR-1, NFUN, 0,
			  TMP_ROWDATA_FLUXERRMAP, 
			  &TMP_ROWDATA_FLUXERRMAP[NVAR-1], 
			  &FLUXERRMAP[NMAP].MAP );  // <== output
      			  
      if ( (OPTMASK & MASK_DUMP_FLUXERRMAP)>0 )
	{ DUMP_FLUXERRMAP(NMAP); }

      FOUNDMAP=0 ;
      malloc_ROWDATA_FLUXERRMAP(-1,NVAR);
    }
    xxxxxxxxxxx end mark xxxxxxxxxxxxxxx */

  }
  // done reading

  // --- check for maps to ignore in reported data error
  parse_IGNORE_FLUXERRMAP(MAPLIST_IGNORE_DATAERR);

  if ( gzipFlag == 0 ) 
    { fclose(fp); }
  else
    { pclose(fp); }

  // print summary of maps
  printSummary_FLUXERRMAP();

  return ;

  //  debugexit(fnam) ;

} // end INIT_FLUXERRMODEL

void  init_fluxerrmodel__(int *optmask, char *fileName, 
			  char *mapList_ignore_dataErr) 
{  INIT_FLUXERRMODEL(*optmask, fileName, mapList_ignore_dataErr); }


// =============================================
void  DUMP_FLUXERRMAP(int IMAP) {

  // Created April 2018
  // Re-write FLUXERRMODEL_FILE contents in fitres format
  // so that it is easier easier to make plots or analyze the maps.

  char *MAPNAME   = FLUXERRMAP[IMAP].NAME ;
  char *FIELD     = FLUXERRMAP[IMAP].FIELDLIST_ORIG ;
  char *BAND      = FLUXERRMAP[IMAP].BANDLIST ;
  int   NVAR_MAP  = FLUXERRMAP[IMAP].NVAR ;
  int   NROW      = FLUXERRMAP[IMAP].MAP.NROW ;

  char *LAST_MAPNAME;
  char dumpFileName[MXPATHLEN];
  FILE *fp;
  int  NVAR_DUMP=0, ivar, irow, OPENFLAG=0 ;
  char VARLIST[100][40];
  //  char fnam[] = "DUMP_FLUXERRMAP";

  // ---------------- BEGIN ---------------
  
  // check if file is opened new or append
  if ( IMAP == 0 ) { OPENFLAG=1; }
  if ( IMAP > 0 ) {
    LAST_MAPNAME = FLUXERRMAP[IMAP-1].NAME ;
    OPENFLAG     = ( strcmp(MAPNAME,LAST_MAPNAME) != 0 ) ;
  }

  sprintf(dumpFileName, "FLUXERRMAP-%s.TEXT", MAPNAME  );

  if ( OPENFLAG ) { 
    fp = fopen(dumpFileName,"wt"); // new file
    printf("\t Dump %s \n", dumpFileName );
    NROW_DUMP_FLUXERRMAP = 0 ;
  } 
  else
    { fp = fopen(dumpFileName,"at"); } // append mode


  if ( OPENFLAG ) {
    // write header when file is opened

    fprintf(fp,"# FLUXERRMAP from file:\n#\t %s\n", FILENAME_FLUXERRMAP);
    fprintf(fp,"#\n");

    sprintf( VARLIST[NVAR_DUMP], "ROW");      NVAR_DUMP++; 
    sprintf( VARLIST[NVAR_DUMP], "IMAP");     NVAR_DUMP++; 
    sprintf( VARLIST[NVAR_DUMP], "FIELD");    NVAR_DUMP++; 
    sprintf( VARLIST[NVAR_DUMP], "BAND");     NVAR_DUMP++; 
    for(ivar=0; ivar<NVAR_MAP; ivar++ ) {
      sprintf(VARLIST[NVAR_DUMP], "%s", FLUXERRMAP[IMAP].VARNAMES[ivar] );
      NVAR_DUMP++ ;
    }   
    fprintf(fp,"NVAR: %d \n", NVAR_DUMP);    
    fprintf(fp,"VARNAMES: ");
    for(ivar=0; ivar < NVAR_DUMP; ivar++ )
      { fprintf(fp,"%s ", VARLIST[ivar]); }
    fprintf(fp,"\n");
  } // end OPENFLAG


  for(irow=0; irow<NROW; irow++ ) {

    NROW_DUMP_FLUXERRMAP++ ;
    fprintf(fp,"ROW: %3d  ", NROW_DUMP_FLUXERRMAP ); 
    fprintf(fp,"%d %s %s ", IMAP, FIELD, BAND );

    for(ivar=0; ivar<NVAR_MAP; ivar++ ) {
      fprintf(fp," %7.3f",  TMP_ROWDATA_FLUXERRMAP[ivar][irow] );
    }

    fprintf(fp,"\n");
  }
  fclose(fp);

  return ;

} // end DUMP_FLUXERRMAP


// =======================================
int index_sparse_FLUXERRMAP(int NMAP, char *MAPNAME) {

  // if MAPNAME already used, return its index.
  // Otherwise increment index.

  int  imap, index_sparse = -1 ;
  int  INDEX_SPARSE, INDEX_MAX = 0 ;
  //  char fnam[] = "index_sparse_FLUXERRMAP" ;

  // ------------- BEGIN --------------

  if ( NMAP == 0 ) { NINDEX_SPARSE_FLUXERRMAP++ ; return(0); }

  for(imap=0; imap < NMAP; imap++ ) {
    INDEX_SPARSE = FLUXERRMAP[imap].INDEX_SPARSE ;

    /*
    printf(" xxx imap=%d of %d :  INDEX_SPARSE=%d \n",
	   imap, NMAP, INDEX_SPARSE );
    */

    if ( INDEX_SPARSE > INDEX_MAX ) 
      { INDEX_MAX = INDEX_SPARSE; }

    if ( strcmp(MAPNAME,FLUXERRMAP[imap].NAME)==0 ) 
      { return(INDEX_SPARSE); }
  }

  if ( index_sparse < 0 ) 
    { index_sparse = INDEX_MAX+1;  NINDEX_SPARSE_FLUXERRMAP++ ; }


  return(index_sparse);

} // end  index_sparse_FLUXERRMAP


// =========================================
void  parse_IGNORE_FLUXERRMAP(char *MAPLIST_IGNORE_DATAERR) {

  // parse input string for comma-separated list of MAP-NAMES
  // to ignore in the data errors. Note that all maps are used 
  // to simulate Poisson noise, so this function only affects
  // the flux errors used in the analysis.

  int MXMAP = MXMAP_FLUXERRMAP;
  int NMAPNAME_IGNORE, NMAP_IGNORE, imap, IMAP, MASK ;
  char **ptrMap ;
  char comma[4]  = "," ;
  char fnam[] = "parse_IGNORE_FLUXERRMAP" ;

  // ------------- BEGIN ------------

  if ( IGNOREFILE(MAPLIST_IGNORE_DATAERR) ) { return ; }

  ptrMap = (char**) malloc ( MXMAP * sizeof(char*) ) ;
  for(imap=0; imap < MXMAP; imap++ ) 
    { ptrMap[imap] = (char*) malloc ( 40 * sizeof(char) ) ;  }

  splitString(MAPLIST_IGNORE_DATAERR, comma, MXMAP_FLUXERRMAP,
	      &NMAPNAME_IGNORE, ptrMap );

  printf("\t Found %d FLUXERRMAPs to ignore in data errors: \n",
	 NMAPNAME_IGNORE );
  for(imap=0; imap < NMAPNAME_IGNORE; imap++ ) {
    NMAP_IGNORE=0;
    for(IMAP=0; IMAP < NMAP_FLUXERRMODEL; IMAP++ ) {
      if ( strcmp(ptrMap[imap],FLUXERRMAP[IMAP].NAME)== 0 ) {
	NMAP_IGNORE++ ;
	MASK = FLUXERRMAP[IMAP].MASK_APPLY ;
	if ( (MASK & MASK_APPLY_DATA_FLUXERRMAP)>0 ) 
	  { MASK -= MASK_APPLY_DATA_FLUXERRMAP; }
	FLUXERRMAP[IMAP].MASK_APPLY = MASK ;
      }
    }
    printf("\t\t %s (%d maps)\n", ptrMap[imap], NMAP_IGNORE );

    if ( NMAP_IGNORE == 0 ) {
      sprintf(c1err,"Invalid MAP-NAME='%s' to ignore data errors:", 
	      ptrMap[imap]);
      sprintf(c2err,"Check user input '%s'", 
	      MAPLIST_IGNORE_DATAERR); 
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  for(imap=0; imap < MXMAP; imap++ ) { free(ptrMap[imap]); }
  free(ptrMap);

  return ;
} // end parse_IGNORE_FLUXERRMAP


// =========================================
void printSummary_FLUXERRMAP(void) {

  int imap ;
  char NAME[60];
  //  char fnam[] = "printSummary_FLUXERRMAP" ;
  char dashLine[] = 
    "------------------------------------------------"
    "---------------------------" ;

  // ----------- BEGIN ------------

  printf("\n");
  printf("       MAPNAME       BANDS         FIELDLIST        "
	 "VARNAMES     SIZE \n");
  printf("# %s \n", dashLine);

  for(imap=0; imap < NMAP_FLUXERRMODEL; imap++ ) {

    sprintf(NAME, "%s:%d"  
	    , FLUXERRMAP[imap].NAME
	    , FLUXERRMAP[imap].INDEX_SPARSE );

    printf(" %-16.16s  %8s %-28.28s %-12.12s  %d\n"
	   , NAME
	   , FLUXERRMAP[imap].BANDLIST
	   , FLUXERRMAP[imap].FIELDLIST
	   , FLUXERRMAP[imap].MAP.VARLIST
	   , FLUXERRMAP[imap].MAP.NROW
	   );
    fflush(stdout);
  }
  printf("# %s \n", dashLine);
  printf(" NINDEX_SPARSE_FLUXERRMAP = %d \n", NINDEX_SPARSE_FLUXERRMAP );
  fflush(stdout);

} // end printSummary_FLUXERRMAP


// =======================================================
int IVARLIST_FLUXERRMAP(char *varName) {
  // Return absolute IVAR index for this variable name
  int  IVAR = -1 ; 
  char fnam[] = "IVAR_FLUXERRMAP" ;
  // ------------ BEGIN -----------
  for(IVAR=0; IVAR < MXVAR_FLUXERRMAP; IVAR++ ) {
    if ( strcmp(varName,VARNAMES_FLUXERRMAP[IVAR])==0 ) 
      { return(IVAR); }
  }

  printf("\n PRE-ABORT DUMP: Valid variables for FLUXERRMAP: \n");
  for(IVAR=0; IVAR < MXVAR_FLUXERRMAP; IVAR++ ) 
    { printf("\t %s \n", VARNAMES_FLUXERRMAP[IVAR] );  }
  
  sprintf(c1err,"Invalid FLUXERRMAP variable: '%s' ", varName);
  sprintf(c2err,"Check FLUXERRMAP_FILE");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 

  return(-1);

} // end IVAR_FLUXERRMAP

// =======================================================
void get_FLUXERRMODEL(int OPT, double FLUXERR_IN, char *BAND, char *FIELD, 
		      int NPAR, double *PARLIST,
		      double *FLUXERR_SIM, double *FLUXERR_DATA ) {

  // Created Feb 2018 by R.Kessler
  // For input flux-error (FLUXERR_IN), function returns corrected
  // flux error using model maps. Note that multiple maps can be
  // applied (with different MAPNAMEs), but only one map per 
  // MAPNAME can be used.
  //
  // Inputs:
  //   OPT         : option mask (not used)
  //   FLUXERR_IN  : input uncertainty, FLUXCAL units
  //   BAND        : single char band (e.g., g or r or i
  //   FIELD       : name of field, or NULL
  //   NPAR        : number of parameters in PARLIST
  //   PARLIST     : list of parameters for maps: PSF, SBMAG, etc ...
  //
  // Ouput:
  //   FLUXERR_SIM : generated flux error
  //   FLUXERR_DATA: reported flux-error in data files
  //
  // For nominal sims, FLUXERR_DATA = FLUXERR_SIM.
  // For systematic studies, however, it may be useful to check 
  // what happens if some error corrections are not accounted for 
  // in the analysis.
  //

  int NMAP      = NMAP_FLUXERRMODEL; 
  int NSPARSE[MXMAP_FLUXERRMAP];
  int IDMAP, istat, isp, imap, MATCH_BAND, MATCH_FIELD ;
  int NVAR, IVAR, ivar, MASK_APPLY ;
  int LDMP = 0 ;
  double FLUXERR_TMP, errModelVal, parList[MXPAR_FLUXERRMAP] ;
  char *tmpString ;
  char fnam[] = "get_FLUXERRMODEL";

  // ----------- BEGIN -------------

  //  LDMP = ( strcmp(BAND,"g") == 0 && fabs(FLUXERR_IN-12.8)<1.0 ); 

  *FLUXERR_SIM  = FLUXERR_IN ;
  *FLUXERR_DATA = FLUXERR_IN ;
  if ( NMAP_FLUXERRMODEL == 0 ) { return; }

  if ( LDMP ) {
    printf(" xxx ------------------------------------ \n");
    printf(" xxx %s DUMP: FIELD=%s-%s \n", fnam, FIELD, BAND );
  }


  if ( NPAR != NPAR_FLUXERRMAP_REQUIRE ) {
    sprintf(c1err,"NPAR=%d but expected %d", NPAR, NPAR_FLUXERRMAP_REQUIRE );
    sprintf(c2err,"grep IPAR_FLUXERRMAP "
	    "$SNANA_DIR/src/sntools_fluxErrModels.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  for(isp=0; isp<NINDEX_SPARSE_FLUXERRMAP; isp++ ) 
    { NSPARSE[isp] = 0 ; }

  for(imap=0; imap < NMAP; imap++ ) {
    MATCH_BAND = MATCH_FIELD = 0 ;

    // check BAND match
    tmpString = FLUXERRMAP[imap].BANDLIST ; 
    if ( strcmp(tmpString,"ALL") == 0 ) 
      { MATCH_BAND = 1; }
    else
      { if( strstr(tmpString,BAND)!=NULL)  { MATCH_BAND=1;}  }
    
    if ( MATCH_BAND == 0 ) { continue ; }

    // check FIELD match
    tmpString = FLUXERRMAP[imap].FIELDLIST ; 
    if ( strcmp(tmpString,"ALL") == 0 ) 
      { MATCH_FIELD = 1; }
    else
      { if( strstr(tmpString,FIELD)!=NULL)  { MATCH_FIELD=1;}  }

    if ( MATCH_FIELD == 0 ) { continue ; }

    // have valid map; increment number of times this MAPNAME is used.
    isp = FLUXERRMAP[imap].INDEX_SPARSE ; 
    NSPARSE[isp]++ ;

    if ( NSPARSE[isp] > 1 ) {
      sprintf(c1err,"%d FLUXERRMODEL maps for %s (BAND=%s, FIELD=%s)", 
	      NSPARSE[isp], FLUXERRMAP[imap].NAME, BAND, FIELD );
      sprintf(c2err,"Only 0 or 1 allowed per MAPNAME.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // load correct variables for this map
    IDMAP = IDGRIDMAP_FLUXERRMODEL_OFFSET + imap ;
    load_parList_FLUXERRMAP(imap, PARLIST, parList);
    istat = interp_GRIDMAP( &FLUXERRMAP[imap].MAP, parList, &errModelVal);

    if ( LDMP || istat<0 ) {
      char cparList[100] ;  cparList[0] = 0 ;
      NVAR      = FLUXERRMAP[imap].NVAR ;
      for(ivar=0; ivar < NVAR-1; ivar++ ) { 
	IVAR      = FLUXERRMAP[imap].IVARLIST[ivar] ;
	tmpString = FLUXERRMAP[imap].VARNAMES[ivar] ;
	sprintf(cparList,"%s %s=%.3f", cparList, tmpString, parList[ivar] ); 
      }
      printf(" xxx imap=%2d  %s(%s-%s)  MJD=%.3f  FLUXERR_IN=%.3f\n", 
	     imap, FLUXERRMAP[imap].NAME, FIELD, BAND,
	     PARLIST[IPAR_FLUXERRMAP_MJD],  FLUXERR_IN) ;
      printf(" xxx     %s  :  errModelVal=%.3f\n", 
	     cparList, errModelVal);
      fflush(stdout) ;
    }
   
    if ( istat < 0 ) {
      sprintf(c1err,"Cannot interpolate FLUXERRMAP");
      sprintf(c2err,"Need to extend range of map.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    MASK_APPLY = FLUXERRMAP[imap].MASK_APPLY ;
    if ( ( MASK_APPLY & MASK_APPLY_SIM_FLUXERRMAP)> 0 ) {
      FLUXERR_TMP  = *FLUXERR_SIM ;
      *FLUXERR_SIM = apply_FLUXERRMODEL(imap, errModelVal, FLUXERR_TMP);
    }
    if ( ( MASK_APPLY & MASK_APPLY_DATA_FLUXERRMAP)> 0 ) {
      FLUXERR_TMP   = *FLUXERR_DATA ;
      *FLUXERR_DATA = apply_FLUXERRMODEL(imap, errModelVal, FLUXERR_TMP);
    }


  } // end imap loop


  if ( LDMP ) {
    printf(" xxx FLUXERR[IN,SIM,DATA] = %.3f, %.3f, %.3f \n",
	   FLUXERR_IN, *FLUXERR_SIM, *FLUXERR_DATA );
    //  debugexit(fnam); 

  }


  return ;

} // end get_FLUXERRMODEL


// fortran wrapper
void  get_fluxerrmodel__(int *OPT, double *FLUXERR_IN, char *BAND, char *FIELD, 
			 int *NPAR, double *PARLIST, 
			 double *FLUXERR_GEN, double *FLUXERR_DATA ) {
  get_FLUXERRMODEL(*OPT, *FLUXERR_IN, BAND, FIELD, *NPAR, PARLIST,
		   FLUXERR_GEN, FLUXERR_DATA);
}


// =========================================================
double apply_FLUXERRMODEL(int imap, double errModelVal, double fluxErr) {

  
  int  NVAR     = FLUXERRMAP[imap].NVAR;
  char *VARNAME = FLUXERRMAP[imap].VARNAMES[NVAR-1];
  char *NAME    = FLUXERRMAP[imap].NAME ;
  int  IPAR1    = IPAR_FLUXERRMAP_ERRSCALE ;
  int  IPAR2    = IPAR_FLUXERRMAP_ERRADD ;
  double FLUXERR_OUT=0.0 ;
  char   fnam[] = "apply_FLUXERRMODEL" ;

  // ------------ BEGIN ------------
  

  if ( strcmp(VARNAME,VARNAMES_FLUXERRMAP[IPAR1])==0 ) {
    // scale error
    FLUXERR_OUT = (fluxErr * errModelVal) ;
  }
  else if ( strcmp(VARNAME,VARNAMES_FLUXERRMAP[IPAR2])==0 ) {
    // add err in quadrature, FLUXCAL units
    FLUXERR_OUT = sqrt(fluxErr*fluxErr + errModelVal*errModelVal) ;
  }
  else {
    sprintf(c1err,"Invalid last column '%s' in %s map", VARNAME, NAME );
    sprintf(c2err,"Last column must be ERRSCALE or ERRADD");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return(FLUXERR_OUT);

} // end apply_FLUXERRMODEL

// =========================================================
void load_parList_FLUXERRMAP(int imap, double *PARLIST, double *parList) {

  // Return parList for map 'imap'.
  // PARLIST is the full list of all parameters.


  int NVAR = FLUXERRMAP[imap].NVAR;
  int IVAR, ivar, EXTRAP_LO=0, EXTRAP_HI=0;
  double PARMIN, PARMAX;
  //  char fnam[] = "load_parList_FLUXERRMAP" ;

  // -------------- BEGIN -------------
  
  for(ivar=0; ivar < NVAR-1; ivar++ ) {
    IVAR          = FLUXERRMAP[imap].IVARLIST[ivar];
    parList[ivar] = PARLIST[IVAR];

    if ( FLUXERRMAP_EXTRAP.FLAG ) {
      // do edge-extrapolation
      PARMIN  = FLUXERRMAP[imap].MAP.VALMIN[ivar] ;
      PARMAX  = FLUXERRMAP[imap].MAP.VALMAX[ivar] ;      
      if ( parList[ivar] < PARMIN ) { parList[ivar] = PARMIN; EXTRAP_LO++; }
      if ( parList[ivar] > PARMAX ) { parList[ivar] = PARMAX; EXTRAP_HI++; }
    }
  }

  FLUXERRMAP_EXTRAP.NOBS_TOT[imap]++ ;
  if ( EXTRAP_LO ) { FLUXERRMAP_EXTRAP.NOBS_EXTRAP_LO[imap]++ ; }
  if ( EXTRAP_HI ) { FLUXERRMAP_EXTRAP.NOBS_EXTRAP_HI[imap]++ ; }

  return ;

} // end load_parList_FLUXERRMAP

// =================================
void END_FLUXERRMODEL(void) {

  int NMAP      = NMAP_FLUXERRMODEL; 
  int imap, N0, N1, NLO, NHI;
  double frac;
  char *NAME, *FIELD, *BAND, TMPNAME[80] ;
  //  char fnam[] = "END_FLUXERRMODEL" ;

  // ------------- BEGIN ------------

  if ( NMAP_FLUXERRMODEL == 0 ) { return; }

  if ( FLUXERRMAP_EXTRAP.FLAG ) {
    printf("\n");
    printf("   FLUXERRMAP                EXTRAP-FRACTION   "
	   "(NLO,NHI) \n");
    printf("   ------------------------------------------------------ \n");
    for(imap=0; imap < NMAP ; imap++ ) {
      N0    = FLUXERRMAP_EXTRAP.NOBS_TOT[imap];
      NLO   = FLUXERRMAP_EXTRAP.NOBS_EXTRAP_LO[imap];
      NHI   = FLUXERRMAP_EXTRAP.NOBS_EXTRAP_HI[imap];
      N1    = NLO + NHI ;
      NAME  = FLUXERRMAP[imap].NAME ;
      FIELD = FLUXERRMAP[imap].FIELDLIST ;
      BAND  = FLUXERRMAP[imap].BANDLIST ;

      if ( N0 > 0 )
	{ frac = (double)(NLO+NHI)/(double)N0 ; }
      else
	{ frac = 0.0 ; }

      sprintf(TMPNAME,"%s(%s-%s)", NAME, FIELD, BAND);

      printf("   %-24.24s  %d/%d = %6.4f   (%d,%d)\n",
	     TMPNAME,   N1,N0,frac,  NLO,NHI); fflush(stdout);
    }
  }

  return;
} // end END_FLUXERRMODEL

void end_fluxerrmodel__(void) { END_FLUXERRMODEL(); }
