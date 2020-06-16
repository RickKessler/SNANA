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
void INIT_FLUXERRMODEL(int OPTMASK, char *fileName, 
		       char *STRING_REDCOV, char *MAPLIST_IGNORE_DATAERR) {

  // Created Feb 2018
  // Read and store maps to correct flux-uncertainties.
  // Corrections can be additive or scale, and depend on
  // and combination of 
  //    BAND, FIELD, MJD, SKYSIG, PSF, ZP, SNR, SB
  //
  // Inputs:
  //   OPTMASK  : bit options
  //
  //   fileName : file containing maps to read
  //
  //   STRING_REDCOV : optional override of REDCOV argument(s).
  //      key val key val etc ...
  //      e.g., 
  //     FLUXERRMODEL_REDCOV(DEEP) griz:0.6  FLUXERRMODEL_REDCOV(SHAL) g:0.2
  //
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
  // Dec 9 2019: abort if map file cannot be opened.
  // Jan 10 2020: parse optional REDCOV 
  // Jan 16 2020: pass redcovString override.

  FILE *fp;
  int  gzipFlag, FOUNDMAP, NTMP, NVAR, NDIM, NFUN, ivar ;
  int  IDMAP, NMAP=0, imap, OPT_EXTRAP=0 ;
  int  USE_REDCOV=1, USE_REDCOV_OVERRIDE=0;
  bool HAS_COLON;
  char PATH[2*MXPATHLEN], c_get[80];  
  char *fullName = FILENAME_FLUXERRMAP ;
  char *name, *fieldList, TMP_STRING[80], LINE[100];
  char MSGERR_FILE[200];
  char colon[] = ":"  ;
  char fnam[] = "INIT_FLUXERRMODEL" ;

  // ------------ BEGIN --------------

  NMAP_FLUXERRMODEL          = 0 ;
  FLUXERR_FIELDGROUP.NDEFINE = 0 ;  
  NINDEX_SPARSE_FLUXERRMAP   = 0 ;
  NREDCOV_FLUXERRMODEL       = 0 ;

  if ( IGNOREFILE(fileName) ) { return ; }

  sprintf(BANNER,"%s:", fnam);
  print_banner(BANNER);

  sprintf(PATH, "%s %s/simlib", PATH_USER_INPUT, PATH_SNDATA_ROOT);
  fp = snana_openTextFile(1,PATH, fileName, fullName, &gzipFlag);

  if ( !fp ) {
    abort_openTextFile("FLUXERRMODEL_FILE", PATH, fileName, fnam);

    /* xxxxxxxx mark delete Feb 1 2020 xxxxxxxxxxxxx
    sprintf(c1err,"Cannot open FLUXERRMODEL_FILE");
    sprintf(c2err,"'%s'", fullName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    xxxxxxxxxxxxx */
  }

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

  // check user-option to disable REDCOV
  if ( IGNOREFILE(STRING_REDCOV)  ) { USE_REDCOV = 0 ; }
  if ( strlen(STRING_REDCOV) == 0 ) { USE_REDCOV = 1 ; }
 
  // check user-option to override REDCOV in FLUXERRMODEL_FILE
  if ( strlen(STRING_REDCOV) > 0 ) { USE_REDCOV_OVERRIDE =  1; }
  
  // - - - - - - - - - - - - - - - - - 
  // start reading file
  while ( fscanf(fp, "%s", c_get) != EOF ) {

    // on any comment char, scoop up rest of line and ignore it
    if ( commentchar(c_get) ) 
      { fgets(TMP_STRING, 80, fp) ; continue ; }
    
    HAS_COLON = ( strstr(c_get,colon) && c_get[0] != '#' );

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

    if ( USE_REDCOV && !USE_REDCOV_OVERRIDE && HAS_COLON ) {
      
      if ( strstr(c_get,"REDCOV") || strstr(c_get,"REDCOR") ) {
	readchar(fp, TMP_STRING);
	strcat(STRING_REDCOV,c_get);
	strcat(STRING_REDCOV," ");
	strcat(STRING_REDCOV,TMP_STRING);
	strcat(STRING_REDCOV," ");
      }
    }

    if ( strcmp(c_get,"MAPNAME:")==0 ) {
      FOUNDMAP = 1 ;
      NMAP = NMAP_FLUXERRMODEL; 
      name = FLUXERRMAP[NMAP].NAME ;
      readchar(fp, name ) ;
      FLUXERRMAP[NMAP].NVAR = NVAR  = 0 ;
      FLUXERRMAP[NMAP].MAP.VARLIST[0]   = 0 ;
      FLUXERRMAP[NMAP].MASK_APPLY   = 3;    // SIM & DATA by default 
      FLUXERRMAP[NMAP].SCALE_FLUXERR_TRUE = 1.0 ;
      FLUXERRMAP[NMAP].SCALE_FLUXERR_DATA = 1.0 ;

      // set defaults to all bands and fields
      sprintf(FLUXERRMAP[NMAP].BANDLIST,  "%s",  ALL_STRING); 
      sprintf(FLUXERRMAP[NMAP].FIELDLIST, "%s",  ALL_STRING );
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

      // if this field is a group; substitute field list.
      set_FIELDLIST_FLUXERRMODEL(TMP_STRING,FLUXERRMAP[NMAP].FIELDLIST);

    }  // end FIELD

    if ( strcmp(c_get,"SCALE_FLUXERR_DATA:")==0 ) {
      if ( FOUNDMAP == 0 ) {
	sprintf(c1err,"%s key not allowed outside MAP", c_get);
	errmsg(SEV_FATAL, 0, fnam, c1err, MSGERR_FILE ); 
      }
      readdouble(fp, 1, &FLUXERRMAP[NMAP].SCALE_FLUXERR_DATA );
    } 
    if ( strcmp(c_get,"SCALE_FLUXERR_TRUE:")==0 ) {
      if ( FOUNDMAP == 0 ) {
	sprintf(c1err,"%s key not allowed outside MAP", c_get);
	errmsg(SEV_FATAL, 0, fnam, c1err, MSGERR_FILE ); 
      }
      readdouble(fp, 1, &FLUXERRMAP[NMAP].SCALE_FLUXERR_TRUE );
    } 

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
      read_GRIDMAP(fp, "FLUXERR", "ROW:", "ENDMAP:", 
		   IDMAP, NDIM, NFUN, OPT_EXTRAP,
		   MXROW_FLUXERRMAP, fnam, 
		   &FLUXERRMAP[NMAP].MAP );  // <== returned		   

      if ( (OPTMASK & MASK_DUMP_MAP_FLUXERRMODEL)>0 )
	{ DUMP_MAP_FLUXERRMODEL(NMAP); }

      FOUNDMAP=0 ;

    }  // end VARNAMES

  }
  // done reading

  // --- check for maps to ignore in reported data error
  parse_IGNORE_FLUXERRMAP(MAPLIST_IGNORE_DATAERR);

  // check for reduced flux covarianve
  parse_REDCOV_FLUXERRMODEL(STRING_REDCOV);

  // close map file.
  if ( gzipFlag == 0 ) 
    { fclose(fp); }
  else
    { pclose(fp); }

  // print summary of maps
  printSummary_FLUXERRMODEL();

  //  if (NREDCOV_FLUXERRMODEL) { debugexit(fnam) ; }


  return ;

} // end INIT_FLUXERRMODEL

void  init_fluxerrmodel__(int *optmask, char *fileName, char *redcorString,
			  char *mapList_ignore_dataErr) 
{  INIT_FLUXERRMODEL(*optmask, fileName, redcorString,
		     mapList_ignore_dataErr); }


// =============================================
void  DUMP_MAP_FLUXERRMODEL(int IMAP) {

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
  //  char fnam[] = "DUMP_MAP_FLUXERRMODEL";

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

} // end DUMP_MAP_FLUXERRMODEL


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


// ===========================================
void parse_REDCOV_FLUXERRMODEL(char *STRING) {

  // Created Jan 10 2010
  //
  // Input STRING is of the form
  //   key val key val etc ...
  //   The key string has the form
  //      FLUXERRMODEL_REDCOV or REDCOV or
  //      FLUXERRMODEL_REDCOV(FIELDNAME) or REDCOV(FIELDNAME)
  //   and val has the form
  //      val = band1:rho1,band2:rho2,band3:rho3
  //   where band and REDCOV are colon separated, and comma
  //   separates different bands. 
  //
  // This syntax allows multiple REDCOV keys in the map file,
  // and also multiple FLUXERRMAP_REDCOV (override) keys in the 
  // sim-input file or command-line.
  //
  //
  // Split input STRING and load contents into struct REDCOV_FLUXERRMAP.
  // Each band can be defined no more than once; otherwise abort.

  int MXRED   = MXREDCOV_FLUXERRMAP;
  int MEMC    = MXCHAR_STRING_REDCOV * sizeof(char) ;

  int   NKEYVAL, NKEY, NITEM, i, ikey;
  char *ptrSplit0[MXRED], *ptrSplit1[MXRED];
  char FIELD[40];
  char space[] = " ", comma[]=  "," ;
  char fnam[] = "parse_REDCOV_FLUXERRMODEL" ;
  int  LDMP = 0 ;

  // ---------- BEGIN -----------

  if ( LDMP ) {  printf(" xxx %s: STRING='%s' \n", fnam, STRING);  }

  if ( IGNOREFILE(STRING) ) { return; }

  if ( strlen(STRING) >= MXCHAR_STRING_REDCOV ) {
    print_preAbort_banner(fnam);
    printf("  STRING = '%s' \n", STRING);
    sprintf(c1err,"STRING length exceeds bound of MXCHAR_STRING_REDCOV=%d", 
	    MXCHAR_STRING_REDCOV );
    sprintf(c2err,"Try fewer REDCOV keys are extend bound.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // allocate string memory for splitting strings
  for(i=0; i < MXRED; i++ ) {
    ptrSplit0[i] = (char*) malloc(MEMC);
    ptrSplit1[i] = (char*) malloc(MEMC);
  }

  // - - - - - - - -
  // split by blank space to get key & val
  splitString(STRING, space, MXRED, &NKEYVAL, ptrSplit0);
  NKEY = NKEYVAL/2; 

  int i2key=0;
  for(ikey=0; ikey < NKEY; ikey++ ) {
    i2key = 2*ikey;

    // get optional field
    extractStringOpt(ptrSplit0[i2key],FIELD);
    if ( IGNOREFILE(FIELD) ) { sprintf(FIELD,"%s",ALL_STRING); }
    
    if ( LDMP ) { 
      printf(" xxx ======================================== \n"); 
      printf(" xxx KEY[%d] = '%s'  -> FIELD='%s'\n", 
	     i2key, ptrSplit0[i2key], FIELD );
      printf(" xxx VAL[%d] = '%s' \n", 
	     i2key, ptrSplit0[i2key+1] );
    }
    
    // split by comma
    splitString(ptrSplit0[i2key+1], comma, MXRED, &NITEM, ptrSplit1);
    for(i=0; i < NITEM; i++ ) 
      {  load_REDCOV_FLUXERRMODEL(ptrSplit1[i],FIELD); }    

  } // end ikey

  // free local memory
  for(i=0; i < MXRED; i++ ) { free(ptrSplit0[i]); free(ptrSplit1[i]);  }

  return ;

} // end parse_REDCOV_FLUXERRMODEL


// ========================================
void load_REDCOV_FLUXERRMODEL(char *ITEM_REDCOV, char *FIELD) {

  // Input is argument item of REDCOV or FLUXERRMODEL_REDCOV
  // If argument is g:0.2,r,0.3,i:0.4
  // then ITEM_REDCOV is either g:0.2 or r:0.3 or i:0.4.
   
  int  NREDCOV = NREDCOV_FLUXERRMODEL ;
  int  N2, ifilt_obs, NBAND_TMP, iband, INDEX_CHECK ;
  double REDCOV ;
  char *ptr_BANDSTRING, *ptr_BANDLIST, *ptrSplit2[2];
  char *ptr_FIELDGRP, *ptr_FIELDLIST, band[2] ;
  char colon[] = ":" ;
  char fnam[]= "load_REDCOV_FLUXERRMODEL";

  // ------------- BEGIN ------------

  ptrSplit2[0] = (char*)malloc( 40*sizeof(char) );
  ptrSplit2[1] = (char*)malloc( 40*sizeof(char) );

  ptr_BANDSTRING = COVINFO_FLUXERRMODEL[NREDCOV].BANDSTRING ;
  ptr_BANDLIST   = COVINFO_FLUXERRMODEL[NREDCOV].BANDLIST ;
  ptr_FIELDGRP   = COVINFO_FLUXERRMODEL[NREDCOV].FIELDGROUP ;
  ptr_FIELDLIST  = COVINFO_FLUXERRMODEL[NREDCOV].FIELDLIST ;

  splitString(ITEM_REDCOV, colon, 2, &N2, ptrSplit2);    
  sprintf(ptr_BANDSTRING,  "%s", ITEM_REDCOV );
  sprintf(ptr_BANDLIST,    "%s", ptrSplit2[0] );
  sprintf(ptr_FIELDGRP,    "%s", FIELD);
  sscanf(ptrSplit2[1], "%le", &REDCOV );

  // set FIELDLIST based on FIELDGRP; 
  sprintf(ptr_FIELDLIST, "%s", ptr_FIELDGRP);
  set_FIELDLIST_FLUXERRMODEL(ptr_FIELDGRP,ptr_FIELDLIST);

  COVINFO_FLUXERRMODEL[NREDCOV].REDCOV = REDCOV ;
  COVINFO_FLUXERRMODEL[NREDCOV].ALL_FIELD = 
    ( strcmp(ptr_FIELDGRP,ALL_STRING) == 0 );

  // keep track of which bands are used, and abort if any
  // band is used more than once.
  NBAND_TMP = strlen(ptr_BANDLIST) ;
  for(iband=0; iband < NBAND_TMP; iband++ ) {
    sprintf(band, "%c", ptr_BANDLIST[iband] );
    ifilt_obs = INTFILTER(band);

    INDEX_CHECK = INDEX_REDCOV_FLUXERRMODEL(band,ptr_FIELDGRP,1,fnam);
    if ( INDEX_CHECK >= 0 ) {       
      sprintf(c1err,"Cannot define %s-%s band more than once.", 
	      band, ptr_FIELDGRP);
      sprintf(c2err,"Each band/FIELD can be defined only once in REDCOV keys.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

  } // end iband
  

  free(ptrSplit2[0]);      free(ptrSplit2[1]);
  NREDCOV_FLUXERRMODEL++ ;

  return ;

} // end load_REDCOV_FLUXERRMODEL


// ==========================================================
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
void printSummary_FLUXERRMODEL(void) {

  int imap ;
  char NAME[60];
  //  char fnam[] = "printSummary_FLUXERRMODEL" ;
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

  // print reduced covariances (Jan 2020)
  int NREDCOV = NREDCOV_FLUXERRMODEL ;
  printf(" NREDCOV = %d  (number of reduced flux-cov)\n", NREDCOV);
  for(imap=0; imap < NREDCOV; imap++ ) {
    printf("\t Excess scatter %d: REDCOV(%s) = %6.3f  (FIELD=%s)\n", imap,
	   COVINFO_FLUXERRMODEL[imap].BANDLIST, 
	   COVINFO_FLUXERRMODEL[imap].REDCOV,
	   COVINFO_FLUXERRMODEL[imap].FIELDGROUP    );
  }
  fflush(stdout);

} // end printSummary_FLUXERRMODEL


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

  print_preAbort_banner(fnam);
  printf("  Valid variables for FLUXERRMAP: \n");
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
		      double *FLUXERR_TRUE, double *FLUXERR_DATA ) {

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
  //   FLUXERR_TRUE: true generated flux error (for sim)
  //   FLUXERR_DATA: reported flux-error in data files
  //
  // For nominal sims, FLUXERR_DATA = FLUXERR_TRUE.
  // For systematic studies, however, it may be useful to check 
  // what happens if some error corrections are not accounted for 
  // in the analysis.
  //
  // Jan 22 2020: refactor to use INDEX_MAP_FLUXERRMODEL.
  // Feb 08 2020: scale errors with SCALE_FLUXERR_DATA[TRUE]

  //  int NMAP      = NMAP_FLUXERRMODEL; 
  int NSPARSE[MXMAP_FLUXERRMAP];
  int IDMAP, istat, isp, imap, NVAR, IVAR, ivar, MASK_APPLY ;
  int LDMP = 0 ;
  double FLUXERR_TMP, errModelVal, parList[MXPAR_FLUXERRMAP] ;
  char *tmpString ;
  char fnam[] = "get_FLUXERRMODEL";

  // ----------- BEGIN -------------

  //  LDMP = ( strcmp(BAND,"g") == 0 && fabs(FLUXERR_IN-12.8)<1.0 ); 

  *FLUXERR_TRUE = FLUXERR_IN ;
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

  imap = INDEX_MAP_FLUXERRMODEL(BAND, FIELD, fnam);
  if ( imap < 0 ) { return ; }
    
  // have valid map; increment number of times this MAPNAME is used.
  isp = FLUXERRMAP[imap].INDEX_SPARSE ; 
  NSPARSE[isp]++ ;

  if ( NSPARSE[isp] > 1 ) {
    sprintf(c1err,"%d FLUXERRMODEL maps for %s (BAND=%s, FIELD=%s)", 
	    NSPARSE[isp], FLUXERRMAP[imap].NAME, BAND, FIELD );
    sprintf(c2err,"Only 0 or 1 allowed per MAPNAME.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // load correct variables for this map;
  // errModelVal is the error scale from map
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
    FLUXERR_TMP  = *FLUXERR_TRUE ;
    *FLUXERR_TRUE = apply_FLUXERRMODEL(imap, errModelVal, FLUXERR_TMP);
    *FLUXERR_TRUE *= FLUXERRMAP[imap].SCALE_FLUXERR_TRUE;
  }

  if ( ( MASK_APPLY & MASK_APPLY_DATA_FLUXERRMAP)> 0 ) {
    FLUXERR_TMP    = *FLUXERR_DATA ;
    *FLUXERR_DATA  = apply_FLUXERRMODEL(imap, errModelVal, FLUXERR_TMP);
    *FLUXERR_DATA *= FLUXERRMAP[imap].SCALE_FLUXERR_DATA;
  }
  

  if ( LDMP ) {
    printf(" xxx FLUXERR[IN,TRUE,DATA] = %.3f, %.3f, %.3f \n",
	   FLUXERR_IN, *FLUXERR_TRUE, *FLUXERR_DATA );
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
void set_FIELDLIST_FLUXERRMODEL(char *FIELDGROUP, char *FIELDLIST) {

  // Created Jan 22 2020
  // For input FIELDGROUP, set ouptut FIELDLIST
  // Note that FIELDLIST is not initialized here, so if there
  // is no FIELDGROUP match, then input FIELDLIST is not modified.
  //
  int igroup ;
  char *tmp_FIELDGROUP, *tmp_FIELDLIST;
  

  for (igroup=0; igroup < FLUXERR_FIELDGROUP.NDEFINE; igroup++ ) {
    tmp_FIELDGROUP  = FLUXERR_FIELDGROUP.NAME[igroup] ;
    tmp_FIELDLIST   = FLUXERR_FIELDGROUP.FIELDLIST[igroup] ;
    if ( strcmp(tmp_FIELDGROUP,FIELDGROUP) == 0 ) 
      { sprintf(FIELDLIST, "%s", tmp_FIELDLIST); }
  }

} // end set_FIELDLIST_FLUXERRMODEL

// =====================================================
int INDEX_MAP_FLUXERRMODEL(char *BAND, char *FIELD, char *FUNCALL) {

  // Created Jan 22, 2020
  // For input BAND and FIELD, return index of map for FLUXERRMODEL.
  // FUNCALL is calling function, and used only for error message.

  int NMAP = NMAP_FLUXERRMODEL; 
  int  imap, IMAP=-9, NMATCH=0 ;
  bool MATCH_BAND, MATCH_FIELD;
  char *tmpString ;
  char fnam[] = "INDEX_MAP_FLUXERRMODEL" ;

  // ------------ BEGIN ---------

  for(imap=0; imap < NMAP; imap++ ) {
    MATCH_BAND = MATCH_FIELD = false ;

    // check BAND match
    tmpString = FLUXERRMAP[imap].BANDLIST ; 
    if ( strcmp(tmpString,ALL_STRING) == 0 ) 
      { MATCH_BAND = true ; }
    else
      { if( strstr(tmpString,BAND)!=NULL)  { MATCH_BAND=true;}  }
    
    if ( !MATCH_BAND ) { continue ; }

    // check FIELD match
    tmpString = FLUXERRMAP[imap].FIELDLIST ; 
    if ( strcmp(tmpString,ALL_STRING) == 0 ) 
      { MATCH_FIELD = true; }
    else
      { if( strstr(tmpString,FIELD)!=NULL)  { MATCH_FIELD=true;}  }

    if ( !MATCH_FIELD ) { continue ; }

    NMATCH++; IMAP=imap;
  }

  if ( NMATCH > 1 ) {
    sprintf(c1err,"Invalid NMATCH=%d for BAND=%s and FIELD=%s",
	    NMATCH, BAND, FIELD);
    sprintf(c2err,"Calling function is %s", FUNCALL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  
  return(IMAP) ;

} // end INDEX_MAP_FLUXERRMODEL

// ==========================================================
int INDEX_REDCOV_FLUXERRMODEL(char *BAND, char *FIELD, int OPT_FIELD,
			      char *FUNCALL) {

  // Return index of REDCOV map for input BAND and FIELD.
  // OPT_FIELD=1 -> check FIELDGROUP (e.g., DEEP, SHALLOW, etc...)
  // OPT_FIELD=2 -> check FIELDLIST (e.g., C3+X3, S1+S2, etc ...)
  //
  // FUNCALL is the calling function, and used only for error message.

  int NREDCOV = NREDCOV_FLUXERRMODEL ;
  int INDEX   = -9, NMATCH=0, i ;
  bool MATCH_BAND, MATCH_FIELD, ALL_FIELD ;
  char *tmp_FIELD, *tmp_BANDLIST ;
  char fnam[] = "INDEX_REDCOV_FLUXERRMODEL" ;

  // --------------- BEGIN ----------------

  for(i=0; i < NREDCOV; i++ ) {

    if ( OPT_FIELD == 1 ) 
      { tmp_FIELD    = COVINFO_FLUXERRMODEL[i].FIELDGROUP ; }
    else
      { tmp_FIELD    = COVINFO_FLUXERRMODEL[i].FIELDLIST ; }

    tmp_BANDLIST = COVINFO_FLUXERRMODEL[i].BANDLIST ;
    ALL_FIELD    = COVINFO_FLUXERRMODEL[i].ALL_FIELD ;
    MATCH_BAND   = ( strstr(tmp_BANDLIST,BAND) != NULL );
    MATCH_FIELD  = ( strstr(tmp_FIELD,FIELD)   != NULL || ALL_FIELD );
    if ( MATCH_BAND && MATCH_FIELD )  { INDEX = i; NMATCH++ ; }
  }

  if ( NMATCH > 1 ) {
    sprintf(c1err,"Invalid NMATCH=%d for BAND=%s and FIELD=%s",
	    NMATCH, BAND, FIELD);
    sprintf(c2err,"Calling function is %s", FUNCALL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return(INDEX) ;

} // end of INDEX_REDCOV_FLUXERRMODEL


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
