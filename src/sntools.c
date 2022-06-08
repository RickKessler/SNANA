// sntools.c

#include "sntools.h"
#include "sntools_spectrograph.h" // Feb 2021
#include "sntools_data.h"
#include "sntools_output.h"

#include <sys/types.h>
#include <sys/stat.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

#include <gsl/gsl_rng.h>      // for Poisson generator
#include <gsl/gsl_randist.h>  // idem

#include "eispack.h"
#include "eispack.c"

/*********************************************************
**********************************************************

               Utilities for SN Analysis
            (see declarations in sntools.h)

**********************************************************
**********************************************************/


// =====================================
int match_cidlist_init(char *fileName, int *OPTMASK, char *varList_store) {

  // Created June 2021
  //
  // if fileName == "", init hash table and return.
  //
  // OPTMASK += 1: use CID_IDSURVEY for matching, else CID
  //    WARNING: returned OPTMASK value is changed if 
  //             IDSURVEY doesnt exist.
  //
  // OPTMASK += 8 -> this is first file, so reset AUTOSTORE
  //
  // varList_store : optional comma-sep list of variables to store,
  //                 to be retreived later for any SNID.
  //
  // If fileName has a dot, read CID list from file;
  // else read comma or space sep list of CIDs from string.
  // File can be FITRES format with VARNAMES key, or a plain list.
  // Use hash table for fast access.
  // After calling this function, do matching with
  //   match = match_cidlist_exec(cid);
  //
  // Function returns number of stored CIDs.
  //

  bool IS_FILE = ( strstr(fileName,DOT) != NULL );
  bool USE_IDSURVEY       = ( *OPTMASK & 1 );
  bool FIRST_FILE         = ( *OPTMASK & 8 );
  bool REFAC        = ( *OPTMASK & 64 );
  bool LEGACY       = !REFAC ;

  bool FORMAT_TABLE = false ; // FITRES table format
  bool FORMAT_NONE  = false ; // cid list with no format
  int  colnum_idsurvey;
  int  NCID, NWD, isn, iwd, MSKOPT = -9 ;
  int  langC = LANGFLAG_PARSE_WORDS_C ;
  int  ILIST = 0, LDMP=0, OPT_AUTOSTORE ;
  double DVAL;
  char CCID[40], STRINGID[40], CVAL[12] ;
  char VARNAME_IDSURVEY[] = "IDSURVEY";
  char fnam[] = "match_cidlist_init";

  // ------------- BEGIN ------------

  if ( LDMP ) {
    printf(" xxx %s: fileName = '%s' \n", fnam, fileName);
    fflush(stdout);
  }

  // init hash table
  if ( strlen(fileName) == 0 )  { 
    match_cid_hash("",-1,0);  
    HASH_STORAGE.NVAR = 0;
    SNTABLE_AUTOSTORE_RESET();  // May 2022
    return 0; 
  }


  ENVreplace(fileName,fnam,1);

  // for input file, figure out format. If FITRES format,
  // check for IDSURVEY column to match by CID_IDSURVEY.
  if ( IS_FILE ) { 
    // ERROR codes: 
    //   colnum = -1 => file does not exist
    //   colnum = -2 => VARNAMES key does not exist
    //   colnum = -3 => VARNAMES key exists, but *varname not found.

    colnum_idsurvey = colnum_in_table(fileName, VARNAME_IDSURVEY);
    
    // if there is no IDSURVEY column, disable user's request 
    // to use IDSURVEY.
    if ( USE_IDSURVEY && colnum_idsurvey < 0 ) 
      { *OPTMASK -= 1;  USE_IDSURVEY = false; }

    if ( colnum_idsurvey == -1 ){ 
      sprintf(c1err,"CID TABLE DOES NOT EXIST");
      sprintf(c2err,"Check file %s",fileName);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    } 

    FORMAT_TABLE = ( colnum_idsurvey == -3 || colnum_idsurvey >= 0); 
    FORMAT_NONE  = !FORMAT_TABLE ;

  }
  else {
    // parse string
    MSKOPT  = MSKOPT_PARSE_WORDS_STRING + MSKOPT_PARSE_WORDS_IGNORECOMMA ;
    NCID    = store_PARSE_WORDS(MSKOPT,fileName);
    for(iwd = 0; iwd < NCID; iwd++ ) {
      get_PARSE_WORD(langC, iwd, CCID);
      match_cid_hash(CCID, ILIST, iwd);
    }
    return NCID ;
  }

  if ( LDMP ) {
    printf(" xxx %s: FORMAT_TABLE = %d  COLNUM_IDSURVEY = %d\n", 
	   fnam, FORMAT_TABLE, colnum_idsurvey );
  }

  // - - - - - - - -
  // if we get here, read file with appropriate format

  int  GZIPFLAG,  MXCHAR_LINE = 400, IDSURVEY, IFILE=0, ISTAT;
  bool IS_ROWKEY, LOAD_CID ;
  char tmpLine[MXCHAR_LINE], key[60], tmpWord[60] ;
  FILE *fp;
  NCID = 0;
  MSKOPT  = MSKOPT_PARSE_WORDS_STRING ;


  // if unformatted, do brute-force read of each CID
  if ( FORMAT_NONE ) {
    fp  = open_TEXTgz(fileName, "rt", &GZIPFLAG);
    while ( fgets(tmpLine,MXCHAR_LINE,fp) ) {
      if ( tmpLine[0] == '#' ) { continue ; }
      // parse words on this line  
      NWD = store_PARSE_WORDS(MSKOPT,tmpLine);
      if ( NWD == 0 ) { continue ; }

      // loop over words on this line     
      for ( iwd = 0; iwd < NWD; iwd++ ) {
	get_PARSE_WORD(langC, iwd, CCID);
	// xxx mark delete May 26 2022 match_cid_hash(STRINGID, ILIST, NCID);
	match_cid_hash(CCID, ILIST, NCID);
	NCID++ ;
        if ( strstr(CCID,COMMA) != NULL || strstr(CCID,COLON) != NULL ||
             strstr(CCID,"=")   != NULL )   {
          sprintf(c1err,"Invalid cid string = '%s'", CCID);
          sprintf(c2err,"Check cid_select_file %s",fileName);
          errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
        }

      } // end loop over CIDs on line
    } // end loop over lines in file

    fclose(fp);
  } // end FORMAT_NONE


  // - - - - - - - - - - - 
  if ( FORMAT_TABLE ) {
    OPT_AUTOSTORE = 1+4; // 1=print each var; 4=append next file
    if ( FIRST_FILE ) { SNTABLE_AUTOSTORE_RESET(); }

    NCID = SNTABLE_AUTOSTORE_INIT(fileName,"CIDLIST", "ALL", OPT_AUTOSTORE);

    if ( IFILE == 0 ) {
      printf("\n %s: store hash table to match CID list.\n", fnam);
      fflush(stdout);
    }

    int ifile, IVAR_IDSURVEY=-9, ISNOFF = 0, ivar, NVAR=0, MEMD, IVAR_TABLE;
    char *ptr_varname;

    // get isn offset to allow for multiple cid_select files
    for(ifile=0; ifile < NFILE_AUTOSTORE-1; ifile++ ) 
      { ISNOFF += SNTABLE_AUTOSTORE[ifile].NROW; }

    // current IFILE file index
    IFILE = NFILE_AUTOSTORE-1 ;

    // get column index for IDSURVEY
    if ( USE_IDSURVEY ) 
      { IVAR_IDSURVEY = IVAR_VARNAME_AUTOSTORE(VARNAME_IDSURVEY); }

    // check for additional columns to store (Apr 29 2022)
    if ( !IGNOREFILE(varList_store) ) {
      if ( IFILE == 0 ) {
	parse_commaSepList(fnam, varList_store, 10,60, 
			   &HASH_STORAGE.NVAR, &HASH_STORAGE.VARNAME_LIST);
	NVAR = HASH_STORAGE.NVAR ;
	MEMD = NCID * sizeof(double);
	HASH_STORAGE.VAL_LIST   = (double**)malloc(NVAR*sizeof(double*) );
	HASH_STORAGE.IVAR_TABLE = (int   * )malloc(NVAR*sizeof(int));
	for(ivar=0; ivar < NVAR; ivar++ ) {
	  ptr_varname = HASH_STORAGE.VARNAME_LIST[ivar] ;
	  IVAR_TABLE  = IVAR_VARNAME_AUTOSTORE(ptr_varname);
	  HASH_STORAGE.IVAR_TABLE[ivar] = IVAR_TABLE ;
	  if ( IVAR_TABLE >= 0 ) {
	    printf("\t store %12s with CID hash table\n", ptr_varname);
	    HASH_STORAGE.VAL_LIST[ivar] = (double*)malloc(MEMD);
	  }
	  else {
	    printf("\t Could not find %12s for CID hash table "
		   " (IVAR_TABLE=%d)\n", ptr_varname, IVAR_TABLE);
	  }
	  fflush(stdout);
	}
      } // end of 1st file init
      else {
	NVAR = HASH_STORAGE.NVAR ;
	MEMD = (ISNOFF+NCID) * sizeof(double);
	// beware: realloc is not tested ...
        for(ivar=0; ivar < NVAR; ivar++ ) {
	  IVAR_TABLE  = IVAR_VARNAME_AUTOSTORE(ptr_varname);
	  if ( IVAR_TABLE >= 0 ) {
	    HASH_STORAGE.VAL_LIST[ivar] = 
	      (double*)realloc(HASH_STORAGE.VAL_LIST[ivar],MEMD);
	  }
	}	
      }

    } // end varList_store

    for(isn=0; isn < NCID; isn++ ) {
      sprintf(CCID,"%s", SNTABLE_AUTOSTORE[IFILE].CCID[isn]);
      if ( USE_IDSURVEY ) {
	DVAL     = SNTABLE_AUTOSTORE[IFILE].DVAL[IVAR_IDSURVEY][isn];
	IDSURVEY = (int)DVAL ;
	sprintf(STRINGID,"%s_%d", CCID, IDSURVEY); 
      }
      else {
	sprintf(STRINGID,"%s", CCID); 
      }

      match_cid_hash(STRINGID, ILIST, ISNOFF+isn);
      
      // check option to store extra columns of info
      for(ivar=0; ivar < NVAR; ivar++ ) {
	IVAR_TABLE = HASH_STORAGE.IVAR_TABLE[ivar];
	if ( IVAR_TABLE >= 0 ) {
	  DVAL  = SNTABLE_AUTOSTORE[IFILE].DVAL[IVAR_TABLE][isn];
	  HASH_STORAGE.VAL_LIST[ivar][ISNOFF+isn] = DVAL ;
	}
      }

    } // end isn loop over sn

    //    printf(" xxx %s: NCID = %d \n", fnam, NCID);
    // debugexit(fnam); // xxxxxx

  } // end FORMAT_TABLE


  if ( LDMP ) {
    printf(" xxx %s: IS_FILE=%d  NCID=%d \n", 
	   fnam, IS_FILE, NCID );
    fflush(stdout);
  }

  return NCID ;

} // end match_cidlist_init


int match_cidlist_exec(char *cid) {
  // Created June 2021
  // Return isn0 index if input cid is on ILIST=0 that was
  // read in match_cidlist_init().
  // Apr 03 2022: replace function bool with int to return isn0
  //        instead of returning bool match=(isn0>0);
  
  int  isn0, ILIST = 1;
  char fnam[] = "match_cidlist_exec";
  // ------------- BEGIN --------------
  isn0  = match_cid_hash(cid, ILIST, -1);
  // printf(" xxx %s: cid=%s -> isn0 = %d \n", fnam, cid, isn0 );
  return isn0 ;
} // end match_cidlist_exec


double  match_cidlist_parval(int isn_match, char *varName, int abort_flag) {
  // Created Apr 29 2022 
  // return table value corresponding to isn_match (from hash table)
  // and *varName column name.
  //
  // abort_flag = 1 -> abort on missing varName
  // abort_flag = 0 -> return NULLVAL in missing varName

  int  IVAR_TABLE, ivar, NVAR = HASH_STORAGE.NVAR;
  char *ptr_varName;
  double DVAL, DVAL_MISSING = -9999.0 ;
  char fnam[] = "match_cidlist_parval";

  // ---------- BEGIN --------------

  for(ivar=0; ivar < NVAR; ivar++ ) {
    IVAR_TABLE  = HASH_STORAGE.IVAR_TABLE[ivar];
    if ( IVAR_TABLE >= 0 ) {
      ptr_varName = HASH_STORAGE.VARNAME_LIST[ivar];
      if ( strcmp(varName,ptr_varName) == 0 ) {
	DVAL = HASH_STORAGE.VAL_LIST[ivar][isn_match];
	return DVAL;
      }
    }
  }

  if ( abort_flag == 0 ) {
    return DVAL_MISSING ;
  }
  else {
    // if we get here, abort on error.
    print_preAbort_banner(fnam);
    for(ivar=0; ivar<NVAR; ivar++ ) {
      ptr_varName = HASH_STORAGE.VARNAME_LIST[ivar];
      printf("\t Valid HASH_STORAGE varName: %s \n", ptr_varName);
    }
    sprintf(c1err,"Invalid HASH_STORAGE varName = %s", varName);
    sprintf(c2err,"in CID match table.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return DVAL_MISSING;

} // end match_cidlist_parval


int match_cidlist_init__(char *fileName, int *OPTMASK, char *varList_store) 
{ return match_cidlist_init(fileName, OPTMASK, varList_store); }

int match_cidlist_exec__(char *cid) 
{ return match_cidlist_exec(cid); }

double  match_cidlist_parval__(int *isn_match, char *varName, int *abort_flag) 
{ return match_cidlist_parval(*isn_match, varName, *abort_flag); }

// ******************************************
#include "uthash.h"
// Jun 2021: define stuff for hash table; used to match CID lists.
struct hash_table_def {
  int id;               // key 
  char name[20];        // array size is max length of CID
  UT_hash_handle hh;    // makes this structure hashable 
} ;
struct hash_table_def *hash_table_users = NULL; 

int match_cid_hash(char *ccid, int ilist, int isn) {

  // Created Jun 2021
  // Use hash table to speed cid matching.             
  // Inputs:           
  //   ilist = 0 for original list, 1,2,3 for additional lists to match
  //   isn   = current SN index                 
  //                                                      
  // Function returns isn index for ilist=0
  // If there is no match, return -9.

  int isn0 = -9;
  struct hash_table_def *s, *tmp;
  char fnam[] = "match_cid_hash" ;

  // ---------------- BEGIN ---------------

  if ( ilist < 0 ) {
    /* free the hash table contents */
    HASH_ITER(hh, hash_table_users, s, tmp) {
      HASH_DEL(hash_table_users, s);
      free(s);
    }
    return(-1);
  }

  if ( ilist == 0 ) {
    // create hash table        
    s     = malloc(sizeof(struct hash_table_def));
    s->id = isn;

    strcpy(s->name, ccid);
    HASH_ADD_STR( hash_table_users, name, s );
    return(isn) ;
  }
  
  // if we get here, match input ccid to ilist=0               
  HASH_FIND_STR( hash_table_users, ccid, s);
  if ( s ) {  isn0 = s->id; }
  
  return isn0;

} // end match_cid_hash

int match_cid_hash__(char *cid, int *ilist, int *isn) {
  int isn0 = match_cid_hash(cid, *ilist, *isn);
  return(isn0);
}

// =================================================
// Feb 2021: functions to mimic python dictionary; func(string) = val
void init_string_dict(STRING_DICT_DEF *DICT, char *NAME, int MAXITEM) {

  // init dictionary with a few things.
  int i;

  sprintf(DICT->NAME, "%s", NAME);
  DICT->N_ITEM     = 0;
  DICT->MAX_ITEM   = MAXITEM ;
  DICT->VALUE_LIST = (double*) malloc( MAXITEM * sizeof(double) );

  DICT->STRING_LIST = (char**) malloc( MAXITEM * sizeof(char*) );
  for(i=0; i < MAXITEM; i++ ) {
    DICT->STRING_LIST[i] = (char*) malloc( 60 * sizeof(char) ); 
    DICT->STRING_LIST[i][0] = 0 ;
  }

  // init "last" values so that get_string_dict is a little faster
  DICT->LAST_STRING[0] = 0;
  DICT->LAST_VALUE     = -999.0 ;

  return ;
} // end init_string_dict

void  load_string_dict(STRING_DICT_DEF *DICT, char *string, double val) {

  // Load string and value for a dictionary element.
  char *NAME    = DICT->NAME ;
  int  MAX_ITEM = DICT->MAX_ITEM;
  int  N_ITEM   = DICT->N_ITEM ;
  int  i;
  bool LDMP = false ;
  char fnam[]   = "load_string_dict" ;

  //-------- BEGIN ---------
  
  if ( MAX_ITEM <=0 || MAX_ITEM > 1000 ) {
    sprintf(c1err,"Invalid MAX_ITEM = %d \n", MAX_ITEM);
    sprintf(c2err,"Probably forgot to call init_string_dict");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(DICT->STRING_LIST[N_ITEM], "%s", string );
  DICT->VALUE_LIST[N_ITEM] = val ;

  if ( LDMP ) {
    printf(" xxx %s %s[%s] = %f \n", fnam, NAME, string, val);
    fflush(stdout);
  }

  N_ITEM++ ;
  DICT->N_ITEM = N_ITEM ;

  if ( N_ITEM >= MAX_ITEM ) {
    print_preAbort_banner(fnam);
    for(i=0; i < N_ITEM-1; i++ ) {
      printf("\t DICT[%s] = %f \n", 
	     DICT->STRING_LIST[i], DICT->VALUE_LIST[i] );
    }
    sprintf(c1err,"N_ITEM=%d exceeds MAX_ITEM for DICT=%s",
	    N_ITEM, NAME );
    sprintf(c2err,"Increase MAX_ITEM or reduce number of load calls");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return;
} // end load_string_dict

void  dump_string_dict(STRING_DICT_DEF *DICT) {

  int   N_ITEM      = DICT->N_ITEM ;
  char *NAME        = DICT->NAME ;
  char *STRING_ORIG = DICT->STRING_LIST[N_ITEM] ;
  int  i;
  char fnam[]   = "dump_string_dict" ;

  //-------- BEGIN ---------
  
  sprintf(BANNER,"%s: dump %s Dictionary", fnam, NAME);
  print_banner(BANNER);
  printf("\t Original Parsed string: %s \n", STRING_ORIG);
  fflush(stdout);

  for (i=0; i < N_ITEM; i++ ) {
    printf("\t %s: %f \n", DICT->STRING_LIST[i], DICT->VALUE_LIST[i]); 
  }

  fflush(stdout);

} // end dump_string_dict

double get_string_dict(int OPT, char *string, STRING_DICT_DEF *DICT) {

  // Inputs
  //   OPT  = 0  -> return -999 if no string element
  //   OPT += 1  -> partial match with strstr instead of full match wit strcmp
  //   OPT += 16 -> abort if no string element
  //
  //  *string = item to return DICT[string] value
  
  int OPTMASK_PARTIAL_MATCH = 1;
  int OPTMASK_ABORT         = 16;

  double VAL  = -999.0 ;
  int i;
  int N_ITEM  = DICT->N_ITEM ;

  bool DO_PARTIAL_MATCH = (OPT & OPTMASK_PARTIAL_MATCH);
  bool MATCH ;
  char *NAME  = DICT->NAME ;
  char *STR ;
  char fnam[] = "get_string_dict" ;

  // ---------- BEGIN ------------

  // quick check if string is same as last call
  if ( strcmp(string,DICT->LAST_STRING) == 0 ) 
    { VAL = DICT->LAST_VALUE; return(VAL); }

  // brute force check over all strings
  for(i=0; i < N_ITEM; i++ ) {
    STR = DICT->STRING_LIST[i] ;
    if ( DO_PARTIAL_MATCH ) 
      { MATCH = ( strstr(string,STR) != NULL ) ; }
    else
      { MATCH = ( strcmp(string,STR) == 0 ) ; }

    if ( MATCH ) {
      VAL = DICT->VALUE_LIST[i];
      sprintf(DICT->LAST_STRING, "%s", string);
      DICT->LAST_VALUE = VAL ;

      return(VAL);
    }
  }

  // if we get here, no string match was found.
  // Check OPT for instructions to abort or return -999

  if ( OPT & OPTMASK_ABORT ) {
    print_preAbort_banner(fnam);
    for(i=0; i < N_ITEM; i++ ) {
      STR = DICT->STRING_LIST[i] ;
      VAL = DICT->VALUE_LIST[i];
      printf("\t DICT %s[%s] = %f \n", NAME, STR, VAL);
    }

    sprintf(c1err,"Invalid string = %s for DICTIONARY=%s",
	    string, DICT);
    sprintf(c2err,"Check valid string items above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return(VAL);

} // end get_string_dict


void parse_string_prescales(char *STRING, STRING_DICT_DEF *DICT) {

  // parse input STRING and return dictionary *DICT.
  // E.g., STRING = RALPH/4.3+SARAH/2.2+PEGGY/9.9
  // --> return pyhton-like dictionary of names and pre-scales;
  //       RALPH: 4.3
  //       SARAH: 2.2
  //       PEGGY: 9.9
  // Both plus and comma are valid separators.

  int  MAXITEM = DICT->MAX_ITEM;
  int    i, NLIST, MEMC = 40*sizeof(char);
  char **ptr_ITEMLIST;
  char sepKey[] = "+" ; // default FIELD separator
  char fnam[] = "parse_string_prescales" ;

  // -------------- BEGIN -------------

  sprintf(DICT->STRING,"%s", STRING);

  // although '+' is default separator, allow commas.
  if ( strstr(STRING,COMMA) != NULL ) 
    { sprintf(sepKey, "%s", COMMA); }

  // allocate memory for each item in FIELDLIST
  ptr_ITEMLIST = (char**)malloc( MAXITEM*sizeof(char*));
  for(i=0; i < MAXITEM; i++ )
    { ptr_ITEMLIST[i] = (char*)malloc(MEMC); }

  splitString(STRING, sepKey, MAXITEM,         // inputs               
	      &NLIST, ptr_ITEMLIST );         // outputs             

  // load dictionary
  int jslash;
  char *tmp, *tmp_ITEM, ITEM[40];
  double dval ;
  for(i=0; i < NLIST; i++ ) {
    tmp_ITEM = ptr_ITEMLIST[i] ; // may include slash
    dval = 1.0;  // default preScale
    sprintf(ITEM, "%s", tmp_ITEM); // 
    if ( strchr(tmp_ITEM,'/') != NULL ) {
      tmp    = strchr(tmp_ITEM, '/');
      jslash = (int)(tmp - tmp_ITEM);
      sscanf( &tmp_ITEM[jslash+1] , "%le", &dval );
      strncmp(ITEM, tmp_ITEM, jslash-1);
      ITEM[jslash] = '\0' ;
    }

    load_string_dict(DICT, ITEM, dval );
		     
  }

  return;

} // end parse_string_prescales 


// =================================================
int store_glob_file_list(char *wildcard) {

  // Created Mar 8 2022
  // Utility to store files based on wildcard;
  // use fetch_glob_file_list to retreive 1 at a time.
  // Includes fortran interface.
  char fnam[] = "store_glob_file_list";
  // --------- BEGIN -----------
  GLOB_LIST.NFILE = glob_file_list(wildcard, &GLOB_LIST.FILE_NAMES);
  return GLOB_LIST.NFILE;

} // end store_glob_file_list

int  store_glob_file_list__(char *wildcard) 
{ return store_glob_file_list(wildcard); }

void get_glob_file(int langFlag, int ifile, char *file_name) {

  // Created Mar 2022
  // Inputs:
  //   langFlag=0 ==> called by C code  ==> do NOT leave pad space
  //   langFlag=1 ==> called by fortran ==> leave pad space
  //   ifile      ==> file index to retrieve
  // Output:
  //   file_name 

  int NFILE = GLOB_LIST.NFILE;
  char fnam[] = "get_glob_file";

  // ----------- BEGIN -------------

  if ( ifile >= NFILE ) { 
    sprintf(c1err,"ifile=%d too large; NFILE=%d", ifile, NFILE);
    sprintf(c2err,"ifile must be < %d", NFILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(file_name, "%s", GLOB_LIST.FILE_NAMES[ifile]);
  if ( langFlag==0 ) 
    {     ; }
  else
    { strcat(file_name," "); }     // extra space for fortran
  return;
} // end get_glob_file

void get_glob_file__(int *langFlag, int *ifile, char *file_name) 
{ get_glob_file(*langFlag, *ifile,file_name); }

void reset_glob_file_list(void) {
  int ifile;
  for(ifile=0; ifile < GLOB_LIST.NFILE; ifile++ ) 
    { free(GLOB_LIST.FILE_NAMES[ifile]) ; }
  free(GLOB_LIST.FILE_NAMES);

  GLOB_LIST.NFILE = 0 ;

} // reset_glob_file_list

void reset_glob_file_list__(void) 
{ reset_glob_file_list(); }

// =================================================
int glob_file_list(char *wildcard, char ***file_list) {

  // Created by P.Armstrong and R.Kessler, 2020
  //  + abort if n_file > MXFILE_LIST

  int    i, n_file = 0; 
  glob_t glob_temp;
  char   fnam[] = "glob_file_list";
  // ----------- BEGIN --------------

  // run the glob
  glob(wildcard, 0, NULL, &glob_temp);
  n_file = glob_temp.gl_pathc;

  // allocate memory for each file name
  *file_list = (char**) malloc(n_file * sizeof(char*));
  for (i=0; i<n_file; i++) {
    (*file_list)[i] = (char*) malloc(MXPATHLEN * sizeof(char));
  }

  // load file names into return arg
  for (i=0; i<n_file; i++) {
    sprintf((*file_list)[i], "%s", glob_temp.gl_pathv[i]);
  }
  return n_file; 
}   // end glob_file_list


// ===============================================
void write_epoch_list_init(char *outFile) {

  // July 11 2020
  // Init utility to write epochs that pass (or fail) cuts
  // defined by calls to write_epoch_list_addvar
  //
  // Initial use is to make epoch-ignore list from snana.exe,
  // then pass this list to classifiers.

  int ivar;
  FILE *FP_OUT;
  char fnam[] = "write_epoch_list_init" ;

  // ------------- BEGIN -----------------
  WRITE_EPOCH_LIST.NVAR = 0;
  WRITE_EPOCH_LIST.NEPOCH_ALL   = 0 ;
  WRITE_EPOCH_LIST.NEPOCH_WRITE = 0 ;
  WRITE_EPOCH_LIST.CUTMASK_ALL  = 0 ;

  for(ivar=0; ivar < MXVAR_WRITE_EPOCH_LIST; ivar++ ) {
    WRITE_EPOCH_LIST.VARNAME[ivar][0]     = 0 ; 
    WRITE_EPOCH_LIST.CUTMODE[ivar]        = 0 ; 
    WRITE_EPOCH_LIST.NEPOCH_CUTFAIL[ivar]      = 0 ;
    WRITE_EPOCH_LIST.NEPOCH_CUTFAIL_ONLY[ivar] = 0 ;
  }

  WRITE_EPOCH_LIST.FP_OUT = fopen(outFile,"wt");
  if( !WRITE_EPOCH_LIST.FP_OUT ) {
    sprintf(c1err,"Could not open outFile");
    sprintf(c2err,"%s", outFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(WRITE_EPOCH_LIST.VARNAMES_LIST,"CID MJD BAND ");
  print_banner(fnam);
  printf("    Open outFile: %s\n", outFile);

  FP_OUT = WRITE_EPOCH_LIST.FP_OUT;
  fprintf(FP_OUT,"#  CUTVAR     CUTMASK    COMMENT \n");
  fflush(FP_OUT);

  return;

} // end write_epoch_list_init


// ==========================================================
void write_epoch_list_addvar(char *VARNAME, double *CUTWIN, 
			     char *CUTMODE_STRING) {

  // July 2020
  // Store info about this epoch-cut variable
  int  NVAR        = WRITE_EPOCH_LIST.NVAR ;
  int  CUTMASK_ADD = (1 << NVAR);
  int  CUTTYPE, CUTMODE ;
  char CUTTYPE_STRING[3][8] = { "NULL", "BITMASK", "WINDOW" };
  char *ptrCUTTYPE;
  char fnam[] = "write_epoch_list_addvar" ;

  // ------------- BEGIN -----------------

  sprintf(WRITE_EPOCH_LIST.VARNAME[NVAR], "%s", VARNAME);
  WRITE_EPOCH_LIST.CUTWIN[NVAR][0] = CUTWIN[0];
  WRITE_EPOCH_LIST.CUTWIN[NVAR][1] = CUTWIN[1];
  WRITE_EPOCH_LIST.CUTMASK_ALL |= CUTMASK_ADD ;



  //  if ( strstr(VARNAME,"PHOTFLAG") != NULL ) 
  if ( CUTWIN[1] < 0.0 && CUTWIN[0] >= 0.0 ) 
    { CUTTYPE = CUTTYPE_BITMASK; }
  else
    { CUTTYPE = CUTTYPE_WINDOW; }
  ptrCUTTYPE = CUTTYPE_STRING[CUTTYPE] ;

  if ( strcmp(CUTMODE_STRING,"REJECT") == 0 ) {
    CUTMODE = CUTMODE_REJECT;
  }
  else if ( strcmp(CUTMODE_STRING,"IGNORE") == 0 ) {
    CUTMODE = CUTMODE_REJECT;
  }
  else if ( strcmp(CUTMODE_STRING,"ACCEPT") == 0 ) {
    CUTMODE = CUTMODE_ACCEPT ;
  }

  printf("   Add CUTVAR %s (CUTMODE=%s, CUTTYPE=%s) \n",
	 VARNAME, CUTMODE_STRING, ptrCUTTYPE ) ;
  fflush(stdout);

  // update comments in outFile
  FILE *FP_OUT = WRITE_EPOCH_LIST.FP_OUT;
  char COMMENT[80];
  if ( CUTTYPE == CUTTYPE_BITMASK ) {
    sprintf(COMMENT,"%s %s %d", 
	    CUTMODE_STRING, ptrCUTTYPE, (int)CUTWIN[0] );
  }
  else {
    sprintf(COMMENT,"%s %s %.2f to %.2f", 
	    CUTMODE_STRING, ptrCUTTYPE, CUTWIN[0], CUTWIN[1] );
  }

  fprintf(FP_OUT,"# %-12.12s   %3d  %s\n", VARNAME, CUTMASK_ADD, COMMENT);
  fflush(FP_OUT);

  // update globals
  strcat(WRITE_EPOCH_LIST.VARNAMES_LIST,VARNAME);
  strcat(WRITE_EPOCH_LIST.VARNAMES_LIST," ");
  WRITE_EPOCH_LIST.CUTTYPE[NVAR] = CUTTYPE ;
  WRITE_EPOCH_LIST.CUTMODE[NVAR] = CUTMODE ;
  sprintf(WRITE_EPOCH_LIST.ROWKEY,"%s", CUTMODE_STRING);
  WRITE_EPOCH_LIST.NVAR++ ;

  return;

} // end write_epoch_list_addvar

// ==========================================================
void write_epoch_list_exec(char *CID,double MJD, char *BAND,double *VALUES) {

  bool DO_WRITE = false ;
  FILE *FP_OUT     = WRITE_EPOCH_LIST.FP_OUT;
  int  NVAR        = WRITE_EPOCH_LIST.NVAR ;
  int  CUTMASK_ALL = WRITE_EPOCH_LIST.CUTMASK_ALL ;

  int  ivar, MASK, CUTTYPE, CUTMODE, MASK_SELECT, CUTMASK_EPOCH ;
  double VAL, *CUTWIN ;
  bool PASSCUT;
  int  LDMP = 0; // (VALUES[0] > 1.01 );
  char STRING_LINE[200], STRING_VAL[20], *VARNAME ;
  char fnam[] = "write_epoch_list_exec" ;

  // ------------- BEGIN -----------------

  // on 1st event, write VARNAMES list (set in addvar)
  if ( WRITE_EPOCH_LIST.NEPOCH_ALL == 0 ) {
    fprintf(FP_OUT,"#\n");
    //  fprintf(FP_OUT,"CUTMASK_ALL: %d\n", WRITE_EPOCH_LIST.CUTMASK_ALL);
    fprintf(FP_OUT,"VARNAMES: %s CUTMASK\n", WRITE_EPOCH_LIST.VARNAMES_LIST);
    fflush(FP_OUT);
  }

  WRITE_EPOCH_LIST.NEPOCH_ALL++ ;

  if ( LDMP )
    { printf(" xxx %s -------- CID=%s --------------- \n", fnam,CID); }

  // apply cut to each variable
  CUTMASK_EPOCH = 0 ;
  for(ivar=0; ivar < NVAR; ivar++ ) {

    VAL     = VALUES[ivar] ;
    CUTTYPE = WRITE_EPOCH_LIST.CUTTYPE[ivar];
    CUTWIN  = WRITE_EPOCH_LIST.CUTWIN[ivar];
    CUTMODE = WRITE_EPOCH_LIST.CUTMODE[ivar];
    VARNAME = WRITE_EPOCH_LIST.VARNAME[ivar];

    if ( CUTTYPE == CUTTYPE_WINDOW ) {
      PASSCUT  = (VAL >= CUTWIN[0] && VAL <= CUTWIN[1] );
    }
    else {
      MASK_SELECT = (int)CUTWIN[0];
      MASK        = (int)VALUES[ivar] ;
      PASSCUT     = (MASK_SELECT & MASK) == 0 ;
    }
    
    // increment CUTMASK for this epoch
    if ( PASSCUT ) { 
      //      CUTMASK_EPOCH |= (1 << ivar) ; 
    }
    else {
      CUTMASK_EPOCH |= (1 << ivar) ;  // failed cuts
      WRITE_EPOCH_LIST.NEPOCH_CUTFAIL[ivar]++ ;
    }

    if ( LDMP ) {
      printf(" xxx %s   %s = %.2f  PASSCUT=%d \n",
	     fnam, VARNAME, VAL, PASSCUT); fflush(stdout);
    }
      
  } // end ivar loop


  if ( CUTMODE == CUTMODE_ACCEPT ) // all cuts must pass
    { DO_WRITE = ( CUTMASK_EPOCH == 0 ) ; }
  else 
    { DO_WRITE = ( CUTMASK_EPOCH > 0  ) ; } // any cut fails

  
  if ( LDMP ) {
    printf(" xxx %s: CUTMASK_EPOCH = %d, DO_WRITE=%d \n", 
	   fnam, CUTMASK_EPOCH, DO_WRITE); fflush(stdout);
  }

  // write to file of DO_WRITE = true
  if ( DO_WRITE ) {
    sprintf(STRING_LINE,"%s:  %10s %10.4f %s ", 
	    WRITE_EPOCH_LIST.ROWKEY, CID, MJD, BAND );

    // tack on variables
    for(ivar=0; ivar < NVAR; ivar++ ) {
      CUTTYPE = WRITE_EPOCH_LIST.CUTTYPE[ivar];
      if ( CUTTYPE == CUTTYPE_BITMASK ) 
	{ sprintf(STRING_VAL, "%6d ", (int)VALUES[ivar]) ; }
      else
	{ sprintf(STRING_VAL, "%.3f ", VALUES[ivar]) ; }
      strcat(STRING_LINE,STRING_VAL);
    }

    // finally, tack on CUTMASK 
    sprintf(STRING_VAL, "%3d ", CUTMASK_EPOCH );
    strcat(STRING_LINE,STRING_VAL);

    // print LINE to outFile
    fprintf(FP_OUT,"%s\n", STRING_LINE);
    fflush(FP_OUT);

    WRITE_EPOCH_LIST.NEPOCH_WRITE++ ;
  }


  return ;

} // end  write_epoch_list_exec

// ==========================================================
void write_epoch_list_summary(void) {

  // Write stat summary of epochs ACCEPTED or REJECTED.

  FILE *FP_OUT       = WRITE_EPOCH_LIST.FP_OUT;
  int  NVAR          = WRITE_EPOCH_LIST.NVAR ;
  int  NEPOCH_WRITE  = WRITE_EPOCH_LIST.NEPOCH_WRITE ;
  int  NEPOCH_ALL    = WRITE_EPOCH_LIST.NEPOCH_ALL ;
  char *ROWKEY       = WRITE_EPOCH_LIST.ROWKEY ; 
  int  ivar, NCUTFAIL;
  double frac;
  char *VARNAME ;
  char fnam[] = "write_epoch_list_summary" ;

  // ------------- BEGIN -----------------

  fprintf(FP_OUT,"\n");
  fprintf(FP_OUT,"# Wrote %d %s epochs out of %d total epochs\n",
	  NEPOCH_WRITE,  ROWKEY, NEPOCH_ALL );

  frac = (double)NEPOCH_WRITE / (double)NEPOCH_ALL ;
  fprintf(FP_OUT,"# %s fraction: %8.5f \n", ROWKEY, frac);
  fflush(FP_OUT); 

  // write cut stats per variable
  for(ivar=0 ; ivar < NVAR; ivar++ ) {
    NCUTFAIL    = WRITE_EPOCH_LIST.NEPOCH_CUTFAIL[ivar];
    VARNAME     = WRITE_EPOCH_LIST.VARNAME[ivar];
    frac        = (double)NCUTFAIL / (double)NEPOCH_ALL ;
    fprintf(FP_OUT,"#    NEPOCH_CUTFAIL[%10s] =  %6d  (frac = %7.5f)\n", 
	    VARNAME, NCUTFAIL, frac);
  }

  fprintf(FP_OUT, "#\n# Done.\n\n");

  fflush(FP_OUT); fclose(FP_OUT);
  return ;  

} // end write_epoch_list_summary


// mangled functions for fortran (snana.car)
void write_epoch_list_init__(char *outFile)
{ write_epoch_list_init(outFile); }
void write_epoch_list_addvar__(char *varName, double *CUTWIN, char *CUTMODE)
{write_epoch_list_addvar(varName,CUTWIN,CUTMODE); }
void write_epoch_list_exec__(char *CID,double *MJD,char *BAND,double *VALUES)
{ write_epoch_list_exec(CID,*MJD,BAND,VALUES); }
void write_epoch_list_summary__(void)
{ write_epoch_list_summary(); }


// ==========================================================
void catVarList_with_comma(char *varList, char *addVarName) {
  char comma[] = "," ;
  if ( strlen(varList) > 0 ) { strcat(varList,comma); }
  strcat(varList,addVarName);
} 

// ==========================================================
int ivar_matchList(char *varName, int NVAR, char **varList) {

  // Created May 2020
  // return ivar such that varList[ivar] that matches varName; 
  // else return -9
  int ivar;
  for(ivar=0; ivar < NVAR; ivar++ ) {
    if ( strcmp(varName,varList[ivar]) == 0 ) { return(ivar); }
  }

  return(-9);
} //

// =========================================================
void init_Cholesky(int OPT, CHOLESKY_DECOMP_DEF *DECOMP) {

  // Feb 2020
  // Inputs are 
  //   + OPT > 0 -> malloc COVMAT2D and load it
  //   + OPT < 0 -> free COVMAT2D and return
  //   + DECOMP->MATSIZE
  //   + DECOMP->COVMAT1D
  //
  // Outputs are
  //  + DECOMP->CHOLESKY2D

  int MATSIZE = DECOMP->MATSIZE ;
  int MEMD0   = MATSIZE * sizeof(double ) ;
  int MEMD1   = MATSIZE * sizeof(double*) ;
  int irow0, irow1;
  gsl_matrix_view chk; 

  // ------------ BEGIN -----------

  if ( OPT < 0 ) {
    for(irow0 = 0; irow0 < MATSIZE; irow0++ )
      { free(DECOMP->CHOLESKY2D[irow0]); }
    free(DECOMP->CHOLESKY2D);
    return ;
  }

  
  DECOMP->CHOLESKY2D = (double**) malloc ( MEMD1 );
  for(irow0 = 0; irow0 < MATSIZE; irow0++ )
    { DECOMP->CHOLESKY2D[irow0] = (double*) malloc ( MEMD0 ); }

  chk = gsl_matrix_view_array ( DECOMP->COVMAT1D, MATSIZE, MATSIZE);
  gsl_linalg_cholesky_decomp ( &chk.matrix)  ;    
  for (irow0=0; irow0 < MATSIZE ; irow0++){
    for (irow1 = 0; irow1 < MATSIZE ; irow1++) {    
      if ( irow0 <= irow1 ) {
	DECOMP->CHOLESKY2D[irow0][irow1] = 
	  gsl_matrix_get(&chk.matrix,irow0,irow1) ;
      }
      else
	{ DECOMP->CHOLESKY2D[irow0][irow1] = 0.0; }
    }    
  }

  return ;

} // end init_Cholesky


void getRan_GaussCorr(CHOLESKY_DECOMP_DEF *DECOMP,
		      double *RanList_noCorr, double *RanList_Corr) {

  // Feb 2020
  // For input list of MATSIZE Gaussian randoms in RanList_noCorr,  
  // return correlated randoms in RanList_Corr

  int MATSIZE   = DECOMP->MATSIZE;
  double GAURAN, tmpMat, tmpRan ;
  int irow0, irow1;
  char fnam[] = "getRan_GaussCorr" ;

  // ------------- BEGIN ------------

  for(irow0=0; irow0 < MATSIZE; irow0++ ) {
    GAURAN = 0.0 ;
    for(irow1 = 0; irow1 < MATSIZE ; irow1++ ) {
      tmpMat = DECOMP->CHOLESKY2D[irow1][irow0] ;
      tmpRan = RanList_noCorr[irow1];
      GAURAN += ( tmpMat * tmpRan) ;
    }
    RanList_Corr[irow0] = GAURAN;
  }

  return ;

} // end getRan_GaussCorr

// ==========================================================
void init_obs_atFLUXMAX(int OPTMASK, double *PARLIST, int VBOSE) {

  // May 24 2019
  // Initialize inputs to estimate PEAKMJD for each event,
  // using only flux and fluxerr.
  //

  char cmethod[100], cwgt[20] ;
  char fnam[] = "init_obs_atFLUXMAX" ;
  bool DO_FMAX       = ( (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX)  > 0 );
  bool DO_FMAXCLUMP2 = ( (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX2) > 0 );
  bool DO_FMAXCLUMP3 = ( (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX3) > 0 );
  // ------------- BEGIN ------------

  if ( OPTMASK == 0 ) { return ; }

  INPUTS_OBS_atFLUXMAX.OPTMASK        = OPTMASK ;
  INPUTS_OBS_atFLUXMAX.MJDWIN         = PARLIST[0];
  INPUTS_OBS_atFLUXMAX.SNRCUT         = PARLIST[1];
  INPUTS_OBS_atFLUXMAX.SNRCUT_BACKUP  = PARLIST[2];

  cmethod[0] = 0;

  if ( VBOSE ) {
    if ( DO_FMAX )  { 
      sprintf(cmethod,"max-flux");  
    }
    else if ( DO_FMAXCLUMP2 || DO_FMAXCLUMP3 ) {
      if ( DO_FMAXCLUMP2 ) { sprintf(cwgt,"wgt=1"); }
      else                 { sprintf(cwgt,"wgt=log(SNR)"); }

      sprintf(cmethod, "Fmax-clump [SNRCUT=%.1f(%.1f), MJDWIN=%.1f, %s]"
	      ,INPUTS_OBS_atFLUXMAX.SNRCUT        
	      ,INPUTS_OBS_atFLUXMAX.SNRCUT_BACKUP
	      ,INPUTS_OBS_atFLUXMAX.MJDWIN, cwgt   );    
    }
    else if ( (OPTMASK & OPTMASK_SETPKMJD_TRIGGER) > 0 )  { 

    }

    printf("\n %s: %s \n", fnam, cmethod);    fflush(stdout);

  } // end VBOSE

  return ;

} // end init_obs_atFLUXMAX


void get_obs_atFLUXMAX(char *CCID, int NOBS, 
		       float *FLUX_LIST, float *FLUXERR_LIST,
		       double *MJD_LIST, int *IFILTOBS_LIST, 
		       int *OBS_atFLUXMAX) {


  // Created May 2019
  // Find observation index (o) at max flux.
  // Must call init_obs_atFLUXMAX first.
  // Always do brute-force search on first iteration.
  // If FLUXMAX2 option is set in OPTMASK (Fmax in clump),
  // do 2nd iteration  over a restricted MJDWIN region 
  // that has the most observations with SNR > SNRCUT_USER. 
  // The goal here is to reject crazy-outler observations that can
  // introduce false PEAKMJDs. If there are no obs passing
  // SNR > SNRCUT_USER, then the process is repeated with
  // an SNR-cut of SNRCUT_BACJUP ... this is to ensure that we 
  // always get a PEAKMJD estimate.
  //
  // Be careful with two different MJD windows:
  // 1) Local MJDSTEP_SNRCUT=10:
  //    NSNRCUT (number of SNR detection) is evaluated in each
  //    10 day window
  // 2) User input MJDWIN_USER (several x 10 days)
  //    This larger window is used to find max number of SNR-detections
  //    from which Fmax is evaluated. SNRCUT_SUM is the sum over
  //    NWIN_COMBINE 10-day windows.
  //
  //  While only the 2nd window (MJDWIN_USER) is used in the end, 
  //  the smaller 10-day window is to improve speed of calculation.
  //
  //  OPTMASK = 8  --> return naive Fmax over all observatins
  //  OPTMASK = 16 --> use Fmax-clump method, wgt=1
  //  OPTMASK = 32 --> use Fmax-clump method, wgt=log10(SNR)
  //
  // Inputs:
  //   CCID         : SN id, for error message
  //   NOBS         : number of observations
  //   FLUX_LIST    : list of NOBS fluxes
  //   FLUXERR_LIST : list of NOBS flux-uncertainties
  //   MJD_LIST     : list of NOBS MJDs
  //   IFILTOBS_LIST: list of NOBS absolute filter indices
  //
  // OUTPUT:
  //   OBS_atFLUXMAX:  obs-index for max flux in each filter band
  //
  // Jul 2 2019: fix bug computing MJDMIN/MAX to work with unsorted MJD_LIST.
  // Jun 7 2020: sort by MJD to handle unsorted data files (e.g., LSST DC2)
  //
  // Sep 15 2020: 
  //  + replace int NSNRCUT with double WGT_SNRCUT to allow playing
  //    with different weights per event. 
  //  + Replace equal wieght per point with log10(SNR) to fix 
  //    problems near season boundary
  //
  // Mar 23 2021: bail if NOBS < 3
  //

  int  OPTMASK         = INPUTS_OBS_atFLUXMAX.OPTMASK ;
  if ( OPTMASK == 0 ) { return ; }
  if ( NOBS    <  3 ) { return ; }

  double MJDWIN_USER     = INPUTS_OBS_atFLUXMAX.MJDWIN ;
  double SNRCUT_USER     = INPUTS_OBS_atFLUXMAX.SNRCUT ;
  double SNRCUT_BACKUP   = INPUTS_OBS_atFLUXMAX.SNRCUT_BACKUP;
  double MJDSTEP_SNRCUT  = 10.0 ; // hard wired param

  // naive max flux
  int USE_MJDatFLUXMAX  = (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX ); 
  // fmax-clump method with wgt=1 per obs
  int USE_MJDatFLUXMAX2 = (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX2); 
  // fmax-clump method with wgt=log(SNR) per obs
  int USE_MJDatFLUXMAX3 = (OPTMASK & OPTMASK_SETPKMJD_FLUXMAX3); 

  int USE_MJDatFMAXCLUMP = (USE_MJDatFLUXMAX2 || USE_MJDatFLUXMAX3) ;

  int USE_BACKUP_SNRCUT, ITER, NITER, IMJD, IMJDMAX=0 ;
  int     NOBS_SNRCUT=0, NSNRCUT_MAXSUM=0, *NSNRCUT;
  double  WGT_SNRCUT_MAXSUM, WGT_SNRCUT_SUM, *WGT_SNRCUT ; 
  int IFILTOBS, o, omin, omax, omin2, omax2, o_sort, NOTHING ;
  int MALLOC=0 ;
  double SNR, SNRCUT=0.0, SNRMAX=0.0, FLUXMAX[MXFILTINDX] ;
  double MJD, MJDMIN, MJDMAX, FLUX, FLUXERR;

  int   *oMIN_SNRCUT=NULL, *oMAX_SNRCUT=NULL ;
  int    NWIN_COMBINE, MXWIN_SNRCUT, MEMI, MEMD ;
  int    LDMP = 0 ; // t(strcmp(CCID,"3530")==0 ) ;
  char fnam[] = "get_obs_atFLUXMAX" ;

  // ------------ BEGIN -------------

  NITER  = 1 ;         // always do max flux on 1st iter
  omin2  = omax2    = -9 ;
  NSNRCUT_MAXSUM    = 0 ; // ??
  WGT_SNRCUT_MAXSUM = 0.0;
  USE_BACKUP_SNRCUT = 0 ;

  // sort by MJD (needed for FmaxClump method)
  MEMI = NOBS*sizeof(int) ;
  int  ORDER_SORT = +1;
  int *INDEX_SORT   = (int*) malloc(MEMI) ;
  sortDouble(NOBS, MJD_LIST, ORDER_SORT, INDEX_SORT );

  // find MJDMIN,MAX (obs may not be time-ordered)
  MJDMIN = +999999.0 ;
  MJDMAX = -999999.0 ;
  for(o=0; o < NOBS; o++ ) {
    MJD = MJD_LIST[o];
    if ( MJD < MJDMIN ) { MJDMIN = MJD; }
    if ( MJD > MJDMAX ) { MJDMAX = MJD; }
  }
  
  NWIN_COMBINE = (int)(MJDWIN_USER/MJDSTEP_SNRCUT + 0.01) ;
  MXWIN_SNRCUT = (int)((MJDMAX-MJDMIN)/MJDSTEP_SNRCUT)+1 ;


  if ( MXWIN_SNRCUT < 0 ) {
    sprintf(c1err,"Crazy MXWIN_SNRCUT = %d",  MXWIN_SNRCUT);
    sprintf(c2err,"MJDMIN/MAX=%.2f/%.2f  MJDSTEP_SNRCUT=%.2f  NOBS=%d",
	    MJDMIN, MJDMAX, MJDSTEP_SNRCUT, NOBS);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  MEMI  = sizeof(int)    * MXWIN_SNRCUT ;
  MEMD  = sizeof(double) * MXWIN_SNRCUT ;

 START:

  if ( LDMP ) {
    printf(" xxx ---------  START CID=%s -------------- \n",
	   CCID ); fflush(stdout);
  }


  if ( USE_MJDatFLUXMAX ) {
    // for naive fluxmax, start with lower SNRCUT_BACKUP
    SNRCUT  = SNRCUT_BACKUP;
    USE_BACKUP_SNRCUT = 1;
  }
  else if ( USE_MJDatFMAXCLUMP ) {
    // fmax-clump method
    NITER   = 2 ;
    IMJDMAX = 0;
    SNRCUT  = SNRCUT_USER;
    if ( USE_BACKUP_SNRCUT ) { SNRCUT = SNRCUT_BACKUP; }

    if ( MALLOC == 0 ) {
      NSNRCUT       = (int*)malloc(MEMI);
      WGT_SNRCUT    = (double*)malloc(MEMD);
      oMIN_SNRCUT   = (int*)malloc(MEMI);
      oMAX_SNRCUT   = (int*)malloc(MEMI);     
      MALLOC        = 1 ; 
    }
    // initialize quantities in each 10-day bin
    for(o=0; o < MXWIN_SNRCUT; o++ ) {
      NSNRCUT[o]       =  0 ;
      WGT_SNRCUT[o]    = 0.0 ;
      oMIN_SNRCUT[o]   = oMAX_SNRCUT[o] = -9 ;
    }

  } // end if-block over USE_MJDatFMAXCLUMP


  //  - - - - - - - - - - - - -- - - - - 
  
  for(ITER=1; ITER <= NITER; ITER++ ) {

    for(IFILTOBS=0; IFILTOBS < MXFILTINDX; IFILTOBS++ ) {
      FLUXMAX[IFILTOBS]       = -9.0 ;
      OBS_atFLUXMAX[IFILTOBS] = -9 ;
    }

    if ( ITER == 1 )
      { omin = 0;  omax = NOBS-1; }
    else 
      { omin = omin2; omax=omax2; }


    if ( omin<0 || omax<0 || omin>= NOBS || omax>= NOBS ) {
      print_preAbort_banner(fnam);

      printf("\t WGT_SNRCUT_MAXSUM = %d \n", WGT_SNRCUT_MAXSUM);
      printf("\t WGT_SNRCUT[]  = %.2f, %.2f, %.2f, %.2f, %.2f ... \n",
	     WGT_SNRCUT[0], WGT_SNRCUT[1], WGT_SNRCUT[2],
	     WGT_SNRCUT[3], WGT_SNRCUT[4]);

      printf("\t NSNRCUT_MAXSUM    = %d \n", NSNRCUT_MAXSUM);
      printf("\t NSNRCUT[]  = %d, %d, %d, %d, %d ... \n",
	     NSNRCUT[0],NSNRCUT[1],NSNRCUT[2],NSNRCUT[3],NSNRCUT[4]);

      printf("\t NOBS_SNRCUT    = %d \n", NOBS_SNRCUT );
      printf("\t SNRMAX         = %.2f \n", SNRMAX );
      printf("\t IMJDMAX        = %d  \n",  IMJDMAX);      
      sprintf(c1err,"omin, omax = %d,%d   ITER=%d", omin, omax, ITER);
      sprintf(c2err,"CID=%s  NOBS=%d", CCID, NOBS);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( LDMP ) {
      printf(" xxx ITER=%d : omin,omax=%3d-%3d   MJDWIN=%.1f-%.1f"
	     " SNRCUT=%.1f \n", 
	     ITER,omin,omax, MJD_LIST[omin], MJD_LIST[omax], SNRCUT ); 
      fflush(stdout);
    }

    NOBS_SNRCUT=0;   SNRMAX = 0.0 ;
    for(o = omin; o <= omax; o++ ) {   // sorted index

      o_sort = INDEX_SORT[o];
      if ( o_sort < 0 || o_sort >= NOBS ) {
	sprintf(c1err,"Invalid o_sort=%d for o=%d \n", o_sort, o);
	sprintf(c2err,"problem sorting my MJD");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      MJD      = MJD_LIST[o_sort];
      FLUX     = (double)FLUX_LIST[o_sort];
      FLUXERR  = (double)FLUXERR_LIST[o_sort];
      IFILTOBS = IFILTOBS_LIST[o_sort];

      if ( IFILTOBS < 1 || IFILTOBS >= MXFILTINDX ) {
	sprintf(c1err,"Invalid IFILTOBS=%d for FLUX[%d]=%f", 
		IFILTOBS, o, FLUX );
	sprintf(c2err,"NOBS=%d", NOBS);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      IMJD = (int)MJD;    
      if ( IMJD    < 40000   ) { continue; }
      if ( FLUXERR < 1.0E-9  ) { continue ; }
      SNR = FLUX/FLUXERR ;
      
      if ( SNR > SNRMAX ) { SNRMAX = SNR; } // diagnostic 
      if ( SNR < SNRCUT ) { continue ; }

      if ( FLUX > FLUXMAX[IFILTOBS] ) { // max flux in each filter
	FLUXMAX[IFILTOBS] = FLUX ;
	OBS_atFLUXMAX[IFILTOBS] = o_sort;
      }

      if ( FLUX > FLUXMAX[0] ) {  // global FLUXMAX
	FLUXMAX[0] = FLUX ;
	OBS_atFLUXMAX[0] = o_sort;
      }


      // count number of SNR>cut epochs in sliding 10-day windows
      if ( USE_MJDatFMAXCLUMP && ITER==1 ) {
      
	IMJD = (int)((MJD - MJDMIN)/MJDSTEP_SNRCUT) ;
	IMJDMAX = IMJD ;
	   
	if ( IMJD < 0 || IMJD >= MXWIN_SNRCUT ) {
	  sprintf(c1err,"Invalid IMJD=%d (must be 0 to %d)", 
		  IMJD, MXWIN_SNRCUT-1);
	  sprintf(c2err,"Something is really messed up.");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
	}							  
	NSNRCUT[IMJD]++ ;
	WGT_SNRCUT[IMJD] += 1.0; // original from May 2019
	if(USE_MJDatFLUXMAX3) {WGT_SNRCUT[IMJD] += log10(SNR);} // Sep 2020

	if ( oMIN_SNRCUT[IMJD] < 0 ) { oMIN_SNRCUT[IMJD] = o; }
	if ( oMAX_SNRCUT[IMJD] < o ) { oMAX_SNRCUT[IMJD] = o; }
      } // end FmaxClump 

      NOBS_SNRCUT++ ;
      
    } // end o-loop over observations

    // =======================

    if ( LDMP ) {
      printf(" xxx ITER=%d : SNRMAX=%.2f NOBS_SNRCUT=%d \n",
	     ITER, SNRMAX, NOBS_SNRCUT ); fflush(stdout);
    }

    // if no obs pass SNRCUT on ITER=1, try again with SNRCUT_BACKUP         
    NOTHING = ( ITER==1 && NOBS_SNRCUT==0 ) ;
    if ( NOTHING ) {
      if ( USE_BACKUP_SNRCUT ) 
	{ goto FREE; }
      else
	{ USE_BACKUP_SNRCUT = 1; goto START; }
    }


    int iwin, iwin_shift, iwin_max, iwin_start;
    int oMIN_TMP, oMAX_TMP, NSNRCUT_SUM ;

    if ( USE_MJDatFMAXCLUMP && ITER==1 ) {
      //   check sliding combined windows for max NSNRCUT
      NSNRCUT_MAXSUM = 0 ; // ??
      WGT_SNRCUT_MAXSUM = 0.0;
      iwin_max = IMJDMAX - (NWIN_COMBINE-1) ;
      if ( iwin_max < 0 ) { iwin_max = 0 ; }
    
      for( iwin_start = 0; iwin_start <= iwin_max; iwin_start++ ) {

	// combine multiple 10-day windows to get a MJDWIN-day window
	NSNRCUT_SUM = 0; oMIN_TMP = oMAX_TMP = -9;
	WGT_SNRCUT_SUM = 0.0 ;

	for(iwin_shift = 0; iwin_shift<NWIN_COMBINE; iwin_shift++ ) {
	  iwin = iwin_start + iwin_shift ;
	  if ( iwin >= MXWIN_SNRCUT ) { continue; }
	  if ( NSNRCUT[iwin] > 0 ) {
	    NSNRCUT_SUM += NSNRCUT[iwin] ;
	    WGT_SNRCUT_SUM += WGT_SNRCUT[iwin];
	    if ( oMIN_TMP < 0 ) { oMIN_TMP = oMIN_SNRCUT[iwin]; }
	    oMAX_TMP = oMAX_SNRCUT[iwin];
	  }
	}

	//	if ( NSNRCUT_SUM > NSNRCUT_MAXSUM ) {
	if ( WGT_SNRCUT_SUM > WGT_SNRCUT_MAXSUM ) {
	  WGT_SNRCUT_MAXSUM = WGT_SNRCUT_SUM ;
	  NSNRCUT_MAXSUM = NSNRCUT_SUM ; // ??
	  omin2 = oMIN_TMP ;
	  omax2 = oMAX_TMP ;

	  if ( LDMP ) {
	    printf("\t xxx NSNRCUT_SUM=%3d for MJDWIN= %.1f to %.1f \n",
		   NSNRCUT_SUM, MJD_LIST[oMIN_TMP], MJD_LIST[oMAX_TMP] );
	    fflush(stdout);
	  }
	}

      } // end loop over iwin_start
    } 

  } // end ITER loop

 FREE:
  if ( MALLOC ) { 
    free(NSNRCUT);     free(WGT_SNRCUT);
    free(oMIN_SNRCUT); free(oMAX_SNRCUT);  
  }
  free(INDEX_SORT);  
  
  return;

} // end get_obs_atFLUXMAX


void init_obs_atfluxmax__(int *OPTMASK, double *PARLIST, int *VBOSE)
{ init_obs_atFLUXMAX(*OPTMASK, PARLIST, *VBOSE); }

void get_obs_atfluxmax__(char *CCID, int *NOBS, float *FLUX, float *FLUXERR,
			 double *MJD, int *IFILTOBS, int *EP_atFLUXMAX) 
{
  get_obs_atFLUXMAX(CCID,*NOBS,FLUX,FLUXERR,MJD,IFILTOBS,EP_atFLUXMAX);
}


// ==============================================
int keyMatch(char *string,char *key, char *keySuffix_optional ) {
  if ( strcmp(string,key)==0 )   
    { return(1); }

  else if ( strlen(keySuffix_optional) > 0 ) {
    // e.g. of keySuffix_optional == ':' then check for key:
    int ISTAT = 0;
    int MEMC  = (strlen(key)+10) * sizeof(char);
    char *KEY = (char*) malloc (MEMC);
    sprintf(KEY,"%s%s", key, keySuffix_optional ) ;
    if ( strcmp(string,KEY) == 0 ) { ISTAT = 1; }
    free(KEY); return(ISTAT);
  }

  // if we get here, there is no match.
  return(0);

} // end keyMatch

// ==============================================
bool keyMatchSim(int MXKEY, char *KEY, char *WORD, int keySource) {

  // Jul 20 2020                                     
  // special key-match for sim-inputs based on keySource:    
  //  FILE -> match with colon, and only MXKEY keys allowed                    
  //  ARG  -> match with or without colon, no limit to repeat keys 
  //    
  // if KEY has multiple space separated values, test them all.  
  // E.g., KEY = "GENSMEAR GEN_SMEAR" is equivalent to two
  // calls with KEY = GENSMEAR, and again with KEY = GEN_SMEAR

  bool IS_FILE = (keySource == KEYSOURCE_FILE);
  bool IS_ARG  = (keySource == KEYSOURCE_ARG );
  bool match = false ;
  int  MSKOPT = MSKOPT_PARSE_WORDS_STRING + MSKOPT_PARSE_WORDS_IGNORECOMMA;
  int  NKEY, ikey;
  char KEY_PLUS_COLON[MXPATHLEN], tmpKey[60];
  char fnam[] = "keyMatchSim";

  // ------------ BEGIN --------------                                          
  NKEY = store_PARSE_WORDS(MSKOPT,KEY);
  for(ikey=0; ikey < NKEY; ikey++ ) {
    get_PARSE_WORD(0, ikey, tmpKey);
    if ( IS_FILE ) {
      // read from file; key must have colon
      sprintf(KEY_PLUS_COLON, "%s%s", tmpKey, COLON);
      if ( NstringMatch( MXKEY, KEY_PLUS_COLON, WORD ) )
        { return(true); }
    }
    else {
      // read from command line arg, colon is optional              
      if ( keyMatch(WORD, tmpKey, COLON ) )  // COLON is optional suffix
        { return(true); }
    }
  }

  return(match);
} // end keyMatchSim                                                                                       

bool NstringMatch(int MAX, char *string, char *key) {

  // Created July 18 2020
  // Check for *string match to input *key, 
  // allowing up to MAX matches. MAX=1 -> uniqueMatch.
  //
  // Examples:
  // Initialize with
  //    NstringMatch(0, "INIT", "sim input file")
  //   
  // allow 1 and only 1 GENVERSION: key with
  //    NstringMatch(1, string, "GENVERSION:"); 
  //
  // allow up to 20 NON1A keys with
  //    NstringMatch(20, string, "NONIA:") ; 
  //
  // abort on OBSOLETE: key with
  //    NstringMatch(0, string, "OBSOLETE:") ; 
  //

  char *msgSource = STRING_UNIQUE.SOURCE_of_STRING ;
  char fnam[] = "NstringMatch";

  // ------------ BEGIN -------------

  if ( strcmp(string,STRINGMATCH_INIT) == 0 ) {
    // interpet *key  as "source of string" to store
    printf("  Initialize %s for %s\n", fnam, key); fflush(stdout);
    sprintf(STRING_UNIQUE.SOURCE_of_STRING, "%s", key);

    // init stuff for optional key dump
    STRING_UNIQUE.NLIST = STRING_UNIQUE.NKEY = 0 ;
    STRING_UNIQUE.DUMPKEY_FLAG = false ;
    if ( MAX < 0 ) { STRING_UNIQUE.DUMPKEY_FLAG = true; } 
        return(0);
  }

  // check option to print EVERY key
  if (  STRING_UNIQUE.DUMPKEY_FLAG ) { dumpUniqueKey(key);  }

  if ( strcmp(string,key) == 0 )
    { checkStringUnique(MAX,string,msgSource,fnam);  return(true); }
  else
    { return(false); }

} // end NstringMatch


int uniqueOverlap (char *string,char *key ) {

  // check Overlap string match up to length of key.
  // Example: 
  //   *key    = 'file='
  //   *string = 'file=anything' will return true.
  // 
  //  The overlap match must start at beginning of *string.
  //  Note that this function uses uniqueMatch util.

  int lenkey = strlen(key);
  int match ;
  char tmpString[1000];
  char fnam[] = "uniqueOverlap" ;
  // ---------- BEGIN -----------

  if ( strcmp(string,STRINGMATCH_INIT) == 0 )  
    {  NstringMatch( 0, STRINGMATCH_INIT, key); return(0); }

  if ( strcmp(string,STRINGMATCH_KEY_DUMP) == 0 )  
    {  NstringMatch(-1, STRINGMATCH_INIT, key); return(0); }

  strncpy(tmpString,string,lenkey); tmpString[lenkey]='\0';
  match = NstringMatch(1,tmpString,key);
  return(match);

} // end uniqueOverlap
 
void  checkStringUnique(int MAX, char *string, char *msgSource, char *callFun) {

  // Apr 2019
  // Utility to store strings and check if input *string is unique. 
  // If *string used more than MAX times,  abort with error message.
  //
  // Inputs:
  //   MAX       : max number of times string allowed to repeat
  //  *string    : string to check
  //  *msgSrouce : message about source; e.g., "sim-input file"
  //  *callFun   : name of calling function
  //
  // The latter two args are used only for error message.
  //
  // string = "INIT" -> initialize arrays.
  // 
  // Jul 17 2020: pass MAX arg.

  int i, NLIST, NFOUND=0 ;
  char *tmpString;
  bool LDMP ;
  char dumpString[] = "XXXZZZQQQ" ;  // "CLEAR" ;
  char fnam[] = "checkStringUnique" ;
  // --------------- BEGIN ------------

  if ( strcmp(string,STRINGMATCH_INIT) == 0 ) { 
    printf("  Initialize %s for %s\n", fnam, msgSource); fflush(stdout);
    STRING_UNIQUE.NLIST = 0;
    return ; 
  }

  NLIST = STRING_UNIQUE.NLIST ;

  if ( NLIST >= MXLIST_STRING_UNIQUE ) {
    sprintf(c1err,"NLIST=%d exceeds bound for dump", NLIST);
    sprintf(c2err,"Try increasing MXLIST_STRING_UNIQUE in sntools.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  for(i=0; i < NLIST; i++ ) {
    NFOUND = 0 ;
    tmpString = STRING_UNIQUE.STRING[i];
    if ( strcmp(tmpString,string) == 0 ) { 
      STRING_UNIQUE.NFOUND_STRING[i]++ ; 
      NFOUND = STRING_UNIQUE.NFOUND_STRING[i] ;
    }

    LDMP = ( strstr(string,dumpString) != NULL );
    if ( LDMP ) {
      printf(" xxx %s: string='%s'  tmpString[%d]=%s  NFOUND=%d  MAX=%d \n", 
	     fnam, string, i, tmpString, NFOUND, MAX); 
    }

    if ( NFOUND > MAX ) {
      sprintf(c1err,"%d  '%s' keys exceeds limit=%d in %s.", 
	      NFOUND, string, MAX, msgSource);
      sprintf(c2err,"Calling function is %s .", callFun);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // even if NFOUND > 0 , return to avoid storing same string 
    // more than onece.
    if ( NFOUND > 0 ) { return ; }
  }

  // - - - - - 
  // if we get here, NFOUND = 0, so add to list

  /*
  printf(" xxx %s: load UNIQUE.STRING[%d] = '%s' \n",
	 fnam, NLIST, string); fflush(stdout);
  */

  if ( MAX == 0 ) {
    sprintf(c1err,"Must remove obsolete '%s' key in %s", 
	    string, msgSource);
    sprintf(c2err,"Calling function is %s .", callFun);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(STRING_UNIQUE.STRING[NLIST],"%s", string);
  STRING_UNIQUE.NFOUND_STRING[NLIST] = 1; // 7.17.2020
  STRING_UNIQUE.NLIST++ ;

  if ( STRING_UNIQUE.NLIST >= MXLIST_STRING_UNIQUE ) {
    sprintf(c1err,"%d unique %s strings exceeds bound.", 
	    STRING_UNIQUE.NLIST, msgSource );
    sprintf(c2err,"Check MXLIST_STRING_UNIQUE=%d", MXLIST_STRING_UNIQUE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  return;

} // end checkStringUnique

void  dumpUniqueKey(char *key) {
  int i, NKEY = STRING_UNIQUE.NKEY ;
  bool EXIST = false ;
  char fnam[] = "dumpUniqueKey";

  // ----------- BEGIN -------------

  //  printf(" xxx %s: key = '%s' \n", fnam, key );

  if ( NKEY >= MXLIST_KEY_UNIQUE ) {
    sprintf(c1err,"NKEY=%d exceeds bound for dump", NKEY);
    sprintf(c2err,"Try increasing MXLIST_KEY_UNIQUE in sntools.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  //printf("\n# DUMP UNIQUE KEYS for %s\n", STRING_UNIQUE.SOURCE_of_STRING );
  for(i=0; i < NKEY; i++ ) {    
    if ( strcmp(key,STRING_UNIQUE.KEY[i]) == 0 ) { EXIST = true; }
  }
  
  if ( !EXIST ) {
    printf(" VALID_KEYNAME:  %4d   %s  \n", NKEY, key ); fflush(stdout);    
    STRING_UNIQUE.KEY[NKEY] = (char*) malloc( 60*sizeof(char) );
    sprintf(STRING_UNIQUE.KEY[NKEY], "%s", key);   
    STRING_UNIQUE.NKEY++ ;
  }

} // end void  dumpStringUniqueList


// ==========================================
void init_lightCurveWidth(void) {
  // ------------ BEGIN ---------------
  // mallac arrays so that they can be realloc later as NOBS increases
  LCWIDTH.LAST_NOBS = 10 ;
  int NTMP = LCWIDTH.LAST_NOBS ;
  LCWIDTH.TLIST_SORTED    = (double*)malloc( NTMP*sizeof(double) );
  LCWIDTH.MAGLIST_SORTED  = (double*)malloc( NTMP*sizeof(double) );
  LCWIDTH.FLUXLIST_SORTED = (double*)malloc( NTMP*sizeof(double) );
  LCWIDTH.INDEX_SORT      = (int   *)malloc( NTMP*sizeof(int)    );
  return ;

} // end init_lightCurveWidth

void   init_lightcurvewidth__(void) { init_lightCurveWidth(); }

double get_lightcurvewidth__(int *OPTMASK, int *NOBS, double *TLIST,
			     double *MAGLIST, int *ERRFLAG, char *FUNCALL ) {
  double width;
  width = get_lightCurveWidth(*OPTMASK,*NOBS,TLIST,MAGLIST,
			      ERRFLAG, FUNCALL);
  return(width);
}


double get_lightCurveWidth(int OPTMASK_LCWIDTH, int NOBS, 
			   double *TLIST, double *MAGLIST,
			   int *ERRFLAG, char *FUNCALL ) {

  // Created Aug 2017
  // For input light curve, return width (days) based on OPTMASK.
  // Inputs:
  //   OPTMASK_LCWIDTH  : bit-mask of options
  //   NOBS     : number of observations
  //   TLIST    : time (days) for each Obs
  //   MAGLIST  : mag for each Obs
  //   FUNCALL  : name of calling function (for error msg only)
  //
  //   Output
  //      ERRFLAG  : user warnings, error, etc.  ERRFLAG=0 for success.
  //      Function returns width in days. 
  //             If ERRFLAG != 0, Width -> -9.
  //
  // Feb 1 2019: add OPT_RISE (4-bit)
  //

  int  OPT_FWHM     = ( (OPTMASK_LCWIDTH & 1) > 0 ) ; // FWHM 
  int  OPT_MINMAX   = ( (OPTMASK_LCWIDTH & 2) > 0 ) ; // TMAX-TMIN
  int  OPT_RISE     = ( (OPTMASK_LCWIDTH & 4) > 0 ) ; // 10-100% rise time

  int  ERRFLAG_LOCAL = 0 ;
  int  obs, isort, obsFmax, FOUND ;
  int  LDUMP_DEBUG =  0 ; // brute-force prints for debug, then abort
  int  ORDER_SORT  = +1 ; // sort with increasing epochs.

  double Width = -9.0, MAG, FLUX, ARG, TMAX, FLUXMAX, T, F, T0,T1, F0,F1; 
  double THALF_lo, THALF_hi, T_FWHM, FHALF, FMIN  ;
  char fnam[] = "get_lightCurveWidth" ;

  // ------------------- BEGIN -------------------

  // error checks
  if ( OPTMASK_LCWIDTH == 0 ) {
    sprintf(c1err,"Invalid OPTMASK=%d passed by %s", 
	    OPTMASK_LCWIDTH, FUNCALL);
    sprintf(c2err,"Check valid options in function %s", fnam);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // extend local 'sort' arrays with realloc when NOBS increases
  if ( NOBS > LCWIDTH.LAST_NOBS  ) {
    int MEMD = NOBS * sizeof(double);
    int MEMI = NOBS * sizeof(int);
    LCWIDTH.TLIST_SORTED    = (double*)realloc(LCWIDTH.TLIST_SORTED,   MEMD );
    LCWIDTH.MAGLIST_SORTED  = (double*)realloc(LCWIDTH.MAGLIST_SORTED, MEMD );
    LCWIDTH.FLUXLIST_SORTED = (double*)realloc(LCWIDTH.FLUXLIST_SORTED,MEMD );
    LCWIDTH.INDEX_SORT      = (int   *)realloc(LCWIDTH.INDEX_SORT,     MEMI );
  }

  // sort TLIST (since T=0 is usually at the end of the list)
  sortDouble(NOBS, TLIST, ORDER_SORT, LCWIDTH.INDEX_SORT );

  if ( LDUMP_DEBUG ) { printf("\n DEBUG DUMP of SORTED Tobs,mag: \n");  }

  // create local time-ordered lists for easier computation,
  // and store fluxes and epoch of peak flux
  TMAX=-999.0; FLUXMAX = 0.0 ; obsFmax=-9;
  for(obs=0; obs<NOBS; obs++ ) {
      isort = LCWIDTH.INDEX_SORT[obs];      
      T     = TLIST[isort];
      MAG   = MAGLIST[isort] ;
      LCWIDTH.TLIST_SORTED[obs]   = T; 
      LCWIDTH.MAGLIST_SORTED[obs] = MAG ;

      ARG = 0.4*(ZEROPOINT_FLUXCAL_DEFAULT - MAG);
      FLUX = pow(TEN,ARG);
      LCWIDTH.FLUXLIST_SORTED[obs] = FLUX ;
      if ( FLUX > FLUXMAX ) { TMAX=T; FLUXMAX=FLUX; obsFmax=obs; }

      if ( LDUMP_DEBUG && MAG < 40.0 ) {
	printf("\t obs=%3d isort=%3d  T=%7.2f  MAG=%.2f  F=%10.6f\n",
	       obs, isort, T, MAG, FLUX ) ;
	fflush(stdout);
      }
  }


  // ----------- start -----------

  if ( OPT_MINMAX ) {
    // Tmin - Tmax option is for illustration only
    if ( NOBS > 1 ) 
      { Width = LCWIDTH.TLIST_SORTED[NOBS-1] - LCWIDTH.TLIST_SORTED[0]; }

  }
  else if ( OPT_FWHM ) {
     
    THALF_lo=THALF_hi = -9999.0; T_FWHM = -9.0 ;  
    FHALF=FLUXMAX/2.0 ; FMIN=0.1*FLUXMAX ;

    // find half max before peak
    FOUND=0;
    for (obs=obsFmax; obs>=0; obs-- ) {
      T=LCWIDTH.TLIST_SORTED[obs] ; F=LCWIDTH.FLUXLIST_SORTED[obs] ; 
      if ( F > FHALF ) { FOUND=0; THALF_lo=-9999.0 ; }
      if ( F < FHALF && F>FMIN && FOUND==0 ) {
	FOUND=1;
	T0 = T; T1=LCWIDTH.TLIST_SORTED[obs+1] ; 
	F0 = F; F1=LCWIDTH.FLUXLIST_SORTED[obs+1] ; 
	if ( LDUMP_DEBUG ) {
	  printf("\t xxx BEFORE PEAK T0=%.2f T1=%.2f  F0=%.3f F1=%.3f \n",
		 T0, T1, F0, F1 ); 
	}
	if ( F1 != F0 ) { THALF_lo = T0 + (T1-T0)*(FHALF-F0)/(F1-F0); }
      }
    }
    // next find half max after peak ... in case of 2nd bump,
    // identify last half-max epoch.
    FOUND=0;
    for (obs=obsFmax; obs<NOBS; obs++ ) {
      T=LCWIDTH.TLIST_SORTED[obs] ; F=LCWIDTH.FLUXLIST_SORTED[obs] ;
      if ( F > FHALF ) { FOUND=0; THALF_hi=-9.0 ; }
      if ( F < FHALF && F>FMIN && FOUND==0 ) {
	FOUND=1;
	T1 = T; T0=LCWIDTH.TLIST_SORTED[obs-1] ; 
	F1 = F; F0=LCWIDTH.FLUXLIST_SORTED[obs-1] ;
	if ( LDUMP_DEBUG ) {
	  printf("\t xxx AFTER PEAK T0=%.2f T1=%.2f  F0=%.3f F1=%.3f \n",
		 T0, T1, F0, F1 ); 
	}
	if ( F1 != F0 ) { THALF_hi = T0 + (T1-T0)*(FHALF-F0)/(F1-F0); }
      }
    }

    if ( THALF_lo > -9998. && THALF_hi > 0.0 ) 
      { T_FWHM = THALF_hi - THALF_lo; }

    Width = T_FWHM ;
    
  } // end OPT_FWHM
  
  else if ( OPT_RISE ) {
    
    // use existing 'HALF' variables, even though it's 10% of maxFlux.
    THALF_lo = -9999.0 ;  
    FHALF = 0.10*FLUXMAX ; FMIN=0.02*FLUXMAX ;

    FOUND=0;
    for (obs=obsFmax; obs>=0; obs-- ) {
      T=LCWIDTH.TLIST_SORTED[obs] ; F=LCWIDTH.FLUXLIST_SORTED[obs] ; 
      if ( F > FHALF ) { FOUND=0; THALF_lo=-9999.0 ; }
      if ( F < FHALF && F>FMIN && FOUND==0 ) {
	FOUND=1;
	T0 = T; T1=LCWIDTH.TLIST_SORTED[obs+1] ; 
	F0 = F; F1=LCWIDTH.FLUXLIST_SORTED[obs+1] ; 
	if ( LDUMP_DEBUG ) {
	  printf("\t xxx BEFORE PEAK T0=%.2f T1=%.2f  F0=%.3f F1=%.3f \n",
		 T0, T1, F0, F1 ); 
	}
	if ( F1 != F0 ) { THALF_lo = T0 + (T1-T0)*(FHALF-F0)/(F1-F0); }
      }
    }

    // Width = TMAX - T(10% of peakFlux)
    if ( THALF_lo > -9998. ) { Width = (TMAX - THALF_lo); }

  } // end OPT_RISE


  *ERRFLAG = ERRFLAG_LOCAL ; // load output arg
  LCWIDTH.LAST_NOBS = NOBS ;   // update last NOBS

  if ( LDUMP_DEBUG ) {
    printf("  Inputs:  OPTMASK=%d  FUNCALL=%s \n", OPTMASK_LCWIDTH, FUNCALL);
    printf("  Output:  WIDTH=%f  ERRFLAG=%d \n", Width, ERRFLAG_LOCAL);
    if ( OPT_FWHM ) {
      printf(" TMAX=%.3f : THALF(lo,hi,sum) = %.3f, %.3f, %.3f \n",
	     TMAX, THALF_lo, THALF_hi, T_FWHM );
    }

    debugexit(fnam); // friendly abort
  }

  return(Width);

} // end get_lightCurveWidth




// ==================================================
void  update_covMatrix(char *name, int OPTMASK, int MATSIZE,
		       double (*covMat)[MATSIZE], double EIGMIN, 
		       int *istat_cov ) {


  // May 2016
  // Part of re-factor to move code out of read_data so that
  // it can be used when reading simdata_biasCor also.
  // 
  // If covMat is not invertible, fix it.
  // Note that input covMat can change !!! 

  // Inputs:
  //  *name = name of SN, used for error msg only.
  //   OPTMASK: 
  //        += 1 -> abort on bad matrix
  //        += 2 -> allow some diag cov entries to be zero
  //        += 4 -> dump info for bad cov           
  //  
  //   MATSIZE = matrix size
  //   covMat  = matrix, MATSIZE x MATSIZE
  //
  // Outputs
  // *covMat    = fixed covMatrix 
  // *istat_cov =  0 if no change; 
  //            = -1 if covMat has changed.
  //            = -9 if any covMat element is -9 ==> undefined 
  //
  // Nov 9 2018:  
  //  + replace fortran rs_ with C version of rs using eispack.c[h] 
  //  + move function from SALT2mu.c into sntools.c
  //


  int nm = MATSIZE ;
  bool matz  ;
  int l, k, m, ierr, ipar, NBAD_EIGVAL ;
  int LDMP = 0; // (strcmp(name,"4249392") == 0);
  int ABORT_ON_BADCOV, ALLOW_ZERODIAG ;


  double eigval[MATSIZE];
  double eigvec[MATSIZE][MATSIZE];
  double eigvalOrig[MATSIZE], covOrig[MATSIZE][MATSIZE];
  double diag[3], EIGEN_MIN ;

  char fnam[] = "update_covMatrix" ;

  // ----------------- BEGIN -----------------

  *istat_cov = 0;

  // check for -9 entries --> undefined.
  // If found, reset cov=0 and return istat_cov=-9
  for (l=0; l<MATSIZE ; ++l)  {         
    for (k=0; k<MATSIZE; ++k)  { 
      if ( fabs(covMat[l][k]+9.0)<1.0E-6 ) {
	*istat_cov = -9 ;
	covMat[l][k] = 0.0 ;
      }
    }
  }
  if ( *istat_cov == -9 ) { return ; }

  // ------------------------
  // check OPTMASK bits
  
  ABORT_ON_BADCOV = ( (OPTMASK & 1) > 0 ) ;
  ALLOW_ZERODIAG  = ( (OPTMASK & 2) > 0 ) ;
  LDMP            = ( (OPTMASK & 4) > 0 ) ;

  // Find eigenvalues and eigenvectors, convention below
  // err[j][i]*eigvec[0][j] = eigval[0]*eigvec[0][i]

  if(LDMP){ printf("\t 1. xxx %s \n", fnam); fflush(stdout); }

  matz = true;
  ierr = rs(nm, &covMat[0][0], eigval, matz, &eigvec[0][0] );

  if(LDMP){ printf("\t 2. xxx %s \n", fnam); fflush(stdout); }

  EIGEN_MIN = 0.0 ;
  if ( ALLOW_ZERODIAG ) { EIGEN_MIN = -1.0E-6 ; }

  NBAD_EIGVAL = 0 ;
  for(ipar=0; ipar < MATSIZE; ipar++ )  { 
    if (eigval[ipar] <= EIGEN_MIN ){NBAD_EIGVAL += 1;}  
  }

  // Check for illegal error matrix
  if (ierr == 0 && NBAD_EIGVAL == 0 ) { return ; }

  // check option to abort
  if ( ABORT_ON_BADCOV ) {
    print_preAbort_banner(fnam);
    for(l=0; l < MATSIZE; l++ ) {
      printf("\t par%d: Eigval=%f  EigVec=%f\n", 
	     l, eigval[l], eigvec[l][l] );
    }
    sprintf(c1err,"%s matrix not positive", name);
    sprintf(c2err,"ierr(rs)=%d  NBAD_EIGVAL=%d", ierr, NBAD_EIGVAL);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  

  // ------------ FIX BAD COV -------------

  *istat_cov = -1 ;

  for (l=0;l<MATSIZE;++l)  { diag[l] = sqrt( covMat[l][l]) ; }
        
  for (k=0;k<MATSIZE;++k)  { 
    eigvalOrig[k] = eigval[k] ;
    if (eigval[k]<EIGMIN) {eigval[k]=EIGMIN;}  
  }
    
  for (l=0;l<MATSIZE;++l)  {
    for (m=0;m<MATSIZE;++m)  {
      covOrig[l][m] = covMat[l][m] ;
      covMat[l][m]=0.0;
      for (k=0;k<MATSIZE;++k) {
	covMat[l][m] += eigvec[k][l]*eigval[k]*eigvec[k][m];
      }
    }
  }


  // --------------------------------------
  // check option to dump info for cov fix
  if ( LDMP ) {
    printf("xxx -------- %s DUMP for SN=%s ------------- \n",  fnam, name);
    printf("xxx rsError=%i  Eigenval(mB,x1,c) = %f %f %f \n",
	   ierr, eigvalOrig[0], eigvalOrig[1], eigvalOrig[2] );
    fflush(stdout);

    printf("xxx  COV_ORIG(mB,x1,c)       -->      COV_FIXED(mB,x1,c) \n");
    for(l=0; l < MATSIZE; l++ ) {
      printf("xxx  ");
      for (m=0;m<MATSIZE;++m) { printf("%8.5f ", covOrig[l][m] ); }
      printf("       ");
      for (m=0;m<MATSIZE;++m) { printf("%8.5f ", covMat[l][m] ); }
      printf("\n");
    }

    printf("\n"); fflush(stdout);

  } // end LDMP


  return ;

} // end update_covMat

void update_covmatrix__(char *name, int *OPTMASK, int *MATSIZE,
			 double (*covMat)[*MATSIZE], double *EIGMIN, 
			 int *istat_cov ) {  
  update_covMatrix(name, *OPTMASK, *MATSIZE, covMat,
		   *EIGMIN, istat_cov); 
} 



// *******************************************************
int store_PARSE_WORDS(int OPT, char *FILENAME) {

  // Read FILENAME (ascii) and store each word to be fetched later
  // with next_PARSE_WORDS.
  // 
  // OPT  = -1 --> one-time init
  // OPT +=  1 --> parse file FILE
  // OPT +=  2 --> FILENAME is a string to parse, parse by space or comma
  // OPT +=  4 --> ignore comma in parsing string: space-sep only
  // OPT +=  8 --> ignore after comment char
  //                
  // Function returns number of stored words separated by either
  // space or comma.
  //
  // Jun 26 2018: free(tmpLine) --> fix awful memory leak
  // May 09 2019: refactor to allow option for ignoring comma in strings.
  // Jul 20 2020: new option to ignore comment char and anything after
  // Jul 31 2020: add abort trap on too-long string length
  // Aug 26 2020: new FIRSTLINE option to read only 1st line of file.
  // Feb 26 2021: for FIRSTLINE, read 5 lines for safety.
  // Feb 18 2022: read 2 lines for FIRSTLINE

  bool DO_STRING       = ( (OPT & MSKOPT_PARSE_WORDS_STRING) > 0 );
  bool DO_FILE         = ( (OPT & MSKOPT_PARSE_WORDS_FILE)   > 0 );
  bool CHECK_COMMA     = ( (OPT & MSKOPT_PARSE_WORDS_IGNORECOMMA) == 0 );
  bool IGNORE_COMMENTS = ( (OPT & MSKOPT_PARSE_WORDS_IGNORECOMMENT) > 0 );
  bool FIRSTLINE       = ( (OPT & MSKOPT_PARSE_WORDS_FIRSTLINE) > 0 );
  int LENF = strlen(FILENAME);
  int NWD, MXWD, iwdStart=0, GZIPFLAG, iwd, nline ;
  char LINE[MXCHARLINE_PARSE_WORDS], *pos, sepKey[4] = " ";
  FILE *fp;
  char fnam[] = "store_PARSE_WORDS" ;
  int LDMP =  0 ;  
  // ------------- BEGIN --------------------

  if ( LENF == 0  ) { PARSE_WORDS.NWD = 0 ; return(0); }

  // if input file (or string) is the same as before,
  // just return since everything is still stored.

  if ( LENF > 0 && strcmp(PARSE_WORDS.FILENAME,FILENAME)==0 ) 
    { return(PARSE_WORDS.NWD); }

  if ( LDMP ) {
    printf(" xxx %s: -----------------------------------------------\n",
	   fnam );
    printf(" xxx %s: OPT=%2d  BUFSIZE=%d    FILENAME='%s'\n", 
	   fnam, OPT, PARSE_WORDS.BUFSIZE, FILENAME ); fflush(stdout);
  }

  if ( OPT < 0 ) {
    PARSE_WORDS.BUFSIZE = PARSE_WORDS.NWD = 0 ;
    PARSE_WORDS.FILENAME[0] = 0 ;
    NWD = 0 ;
    return(NWD);
  }
  else if ( DO_STRING ) {
    // FILENAME is a string to parse
   
    char *tmpLine = (char*) malloc( (LENF+10)*sizeof(char) );
    sprintf(tmpLine, "%s", FILENAME); // Mar 13 2019
    if ( CHECK_COMMA && strchr(tmpLine,COMMA[0]) != NULL ) 
      { sprintf(sepKey,"%s", COMMA); }

    malloc_PARSE_WORDS() ;    MXWD = PARSE_WORDS.BUFSIZE ;

    splitString2(tmpLine, sepKey, MXWD,
		 &NWD, &PARSE_WORDS.WDLIST[0] ); // <== returned
    PARSE_WORDS.NWD = NWD ;
    
    // check that string lengths don't overwrite bounds (July 2020)
    for(iwd=0; iwd < NWD; iwd++ ) {
      int lwd = strlen(PARSE_WORDS.WDLIST[iwd]);
      if ( lwd >= MXCHARWORD_PARSE_WORDS ) {
	print_preAbort_banner(fnam); 
	printf("  split tmpLine = \n '%s' \n", tmpLine);
	sprintf(c1err,"split string len = %d for iwd=%d " 
		"(see tmpLine above)",	 lwd, iwd);
	sprintf(c2err,"Check MXCHARWORD_PARSE_WORDS=%d", 
		MXCHARWORD_PARSE_WORDS);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
    }
    
    if ( NWD >= MXWD ) {
      sprintf(c1err,"NWD=%d exceeds bound.", NWD);
      sprintf(c2err,"Check PARSE_WORDS.BUFSIZE=%d", MXWD);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
    free(tmpLine);
  }
  else if ( DO_FILE ) {
    // read text file
    fp = open_TEXTgz(FILENAME,"rt", &GZIPFLAG );
    if ( !fp ) {
      sprintf(c1err,"Could not open text file");
      sprintf(c2err,"%s", FILENAME);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    NWD = PARSE_WORDS.NWD = nline = LINE[0] = 0 ;
    while( fgets(LINE, MXCHARLINE_PARSE_WORDS, fp)  != NULL ) {
      if ( strlen(LINE) == 0 ) { continue; }
      nline++ ;
      malloc_PARSE_WORDS();
      if ( (pos=strchr(LINE,'\n') ) != NULL )  { *pos = '\0' ; }
      if ( PARSE_WORDS.NWD < MXWORDFILE_PARSE_WORDS ) 
	{ iwdStart = PARSE_WORDS.NWD; }
      splitString2(LINE, sepKey, MXWORDLINE_PARSE_WORDS, 
		   &NWD, &PARSE_WORDS.WDLIST[iwdStart] ); // <== returned

      if ( IGNORE_COMMENTS ) { // 7.2020
	int NWD_TMP = 0 ; bool FOUND_COMMENT=false;  char *ptrWD ;
	for(iwd = iwdStart; iwd < (iwdStart + NWD); iwd++ ) {
	  ptrWD = PARSE_WORDS.WDLIST[iwd] ; 
	  if ( commentchar(ptrWD) ) { FOUND_COMMENT = true ; }
	  if ( !FOUND_COMMENT ) { NWD_TMP++; }
	}
	NWD = NWD_TMP; // reset NWD to ignore comments
      }
      PARSE_WORDS.NWD += NWD;
      if ( FIRSTLINE && nline > 2 ) { break; }
    } // end while
    NWD = PARSE_WORDS.NWD ;

    check_EOF(fp, FILENAME, fnam, nline); 

    fclose(fp);
  }
  else {
    sprintf(c1err,"Invalid OPT=%d with FILENAME='%s'", OPT, FILENAME);
    sprintf(c2err,"grep MSKOPT_PARSE $SNANA_DIR/src/sntools.c");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  if ( NWD >= MXWORDFILE_PARSE_WORDS ) {
    sprintf(c1err,"NWD=%d exceeds bound, MXWORDFILE_PARSE_WORDS=%d",
	    NWD, MXWORDFILE_PARSE_WORDS );
    sprintf(c2err,"Check '%s' ", FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // abort on tabs
  int NTAB=0;
  for(iwd=0; iwd < NWD; iwd++ ) 
    {  if ( PARSE_WORDS.WDLIST[iwd][0] == '\t' ) { NTAB++ ; }   }
  if ( NTAB > 0 ) {
    sprintf(c1err,"Found %d invalid tabs.", NTAB);
    sprintf(c2err,"Check '%s' ", FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);    
  }


  if ( LENF < MXPATHLEN ) 
    { sprintf(PARSE_WORDS.FILENAME, "%s", FILENAME); }
  else
    { PARSE_WORDS.FILENAME[0] = 0 ; }


  if ( LDMP ) {
    printf(" xxx %s: NWD_STORE = %d \n", fnam, NWD);
    fflush(stdout);
  }

  return(NWD);;

} // end store_PARSE_WORDS

int store_parse_words__(int *OPT, char *FILENAME) 
{ return store_PARSE_WORDS(*OPT, FILENAME); }


void malloc_PARSE_WORDS(void) {
  int ADDBUF    = ADDBUF_PARSE_WORDS ;
  int MXCHARWD  = MXCHARWORD_PARSE_WORDS ;
  int NWD       = PARSE_WORDS.NWD ;
  int iwd, WD0, WD1, BUFSIZE, IFLAG=0 ;
  //  char fnam[] = "malloc_PARSE_WORDS" ;

  // ------------- BEGIN ----------------

  /*
  printf(" xxx %s: BUFSIZE = %d  (NWD=%d) \n", 
  fnam, PARSE_WORDS.BUFSIZE, NWD ); */

  if ( PARSE_WORDS.BUFSIZE == 0 ) {
    PARSE_WORDS.WDLIST    = (char**) malloc( ADDBUF*sizeof(char*) );
    PARSE_WORDS.BUFSIZE  += ADDBUF ;
    WD0 = 0;  WD1 = PARSE_WORDS.BUFSIZE ;
    IFLAG=1;
  }
  else if ( NWD > PARSE_WORDS.BUFSIZE - MXWORDLINE_PARSE_WORDS ) {
    WD0                   = PARSE_WORDS.BUFSIZE ;
    PARSE_WORDS.BUFSIZE  += ADDBUF ;
    WD1                   = PARSE_WORDS.BUFSIZE ;

    BUFSIZE               = PARSE_WORDS.BUFSIZE ;
    PARSE_WORDS.WDLIST    = 
      (char**) realloc(PARSE_WORDS.WDLIST, BUFSIZE*sizeof(char*) );
    IFLAG=2;    
  }

  else {
    IFLAG=3;
    return ;
  }

  
  for(iwd=WD0; iwd < WD1; iwd++ ) {
    PARSE_WORDS.WDLIST[iwd]  = (char*) malloc( MXCHARWD*sizeof(char) );
    PARSE_WORDS.WDLIST[iwd][0] = 0 ;
  }

  return ;

} // end malloc_PARSE_WORDS


void get_PARSE_WORD(int langFlag, int iwd, char *word) {

  // Code-language flag:
  // langFlag=0 ==> called by C code  ==> do NOT leave pad space
  // langFlag=1 ==> called by fortran ==> leave pad space

  int NWD = PARSE_WORDS.NWD ;
  char fnam[] = "get_PARSE_WORD" ;

  if ( iwd >= NWD ) {
    sprintf(c1err,"iwd=%d exceeds NWD_STORE=%d", iwd, NWD);
    sprintf(c2err,"Check '%s' ", PARSE_WORDS.FILENAME);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }
  
  sprintf(word, "%s", PARSE_WORDS.WDLIST[iwd] );
  if ( langFlag==0 ) 
    { trim_blank_spaces(word); }  // remove <CR>
  else
    { strcat(word," "); }     // extra space for fortran
  
} // end get_PARSE_WORD

void get_PARSE_WORD_INT(int langFlag, int iwd, int *i_val) {
  char word[100];   get_PARSE_WORD(langFlag, iwd, word);
  sscanf(word, "%d", i_val);
}
void get_PARSE_WORD_FLT(int langFlag, int iwd, float *f_val) {
  char word[100];   get_PARSE_WORD(langFlag, iwd, word);
  sscanf(word, "%f", f_val);
}
void get_PARSE_WORD_NFLT(int langFlag, int NFLT, int iwd, float *f_val) {
  // Created Dec 10 2021
  // return NFLT floats from store_PARSE_WORDS starting at iwd
  char word[100];   int i;
  for(i=0; i < NFLT; i++ ) { 
    get_PARSE_WORD_FLT(langFlag, iwd+i, &f_val[i]);
  }
}  // end get_PARSE_WORD_NFLT

void get_PARSE_WORD_NFILTDEF(int langFlag, int iwd, float *f_val) {
  // Created Dec 10 2021
  // load NFILT words and store into f_val[ifilt_obs]
  int ifilt, ifilt_obs, NFILT = SNDATA_FILTER.NDEF;
  for ( ifilt=0; ifilt < NFILT; ifilt++ ) {
    ifilt_obs  = SNDATA_FILTER.MAP[ifilt];
    get_PARSE_WORD_FLT(langFlag, iwd+ifilt, &f_val[ifilt_obs] ); 
  }

} // end get_PARSE_WORD_NFILTDEF

void get_PARSE_WORD_DBL(int langFlag, int iwd, double *d_val) {
  char word[100];   get_PARSE_WORD(langFlag, iwd, word);
  sscanf(word, "%le", d_val);
}

void get_parse_word__(int *langFlag, int *iwd, char *word) 
{ get_PARSE_WORD(*langFlag, *iwd, word); }

void get_parse_word_int__(int *langFlag, int *iwd, int *i_val) 
{ get_PARSE_WORD_INT(*langFlag, *iwd, i_val); }

void get_parse_word_flt__(int *langFlag, int *iwd, float *f_val) 
{ get_PARSE_WORD_FLT(*langFlag, *iwd, f_val); }

void get_parse_word_dbl__(int *langFlag, int *iwd, double *d_val) 
{ get_PARSE_WORD_DBL(*langFlag, *iwd, d_val); }


// ******************************************
void parse_multiplier(char *inString, char *key, double *multiplier) {

  // Mar 2019
  //
  // Examples:
  //  inString(I)   key(I)   multiplier(O)
  // ----------------------------------------------
  //   CC_S15        CC_S15       1.0
  //   .3*CC_S15     CC_S15       0.3
  //   CC_S15*.3     CC_S15       0.3
  //   BLABLA        CC_S15       0.0

  int    MEMC  = 40*sizeof(char);
  int    LDMP = 0;
  int    NSPLIT ;
  char  *splitValues[2], *cnum = NULL ;
  char star[] = "*" ;
  char fnam[] = "parse_multiplier" ;

  // ----------------- BEGIN ----------------

  if ( strstr(inString,key) == NULL ) { 
    *multiplier = 0.0;
  }
  else if ( strcmp(inString,key) == 0 ) {
    *multiplier = 1.0;
  }
  else {
    // do the parsing; start by splitting string around *
    splitValues[0] = (char*)malloc(MEMC);  
    splitValues[1] = (char*)malloc(MEMC);
    splitString(inString, star, 3,      // inputs
		&NSPLIT, splitValues );      // outputs         
    
    if ( strcmp(splitValues[0],key) == 0 ) 
      { cnum = splitValues[1]; }
    else if ( strcmp(splitValues[1],key) == 0 ) 
      { cnum = splitValues[0]; }
    else {
      sprintf(c1err,"Cannot find key='%s' after splitting string='%s'",
	      key, inString);
      sprintf(c2err,"Something is messed up.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }
    sscanf(cnum, "%le", multiplier ); // load output arg

    free(splitValues[0]);  free(splitValues[1]);
  }

  if ( LDMP ) {
    printf(" xxx ----------------------------- \n");
    printf(" xxx %s DUMP: \n", fnam);
    printf(" xxx inString = '%s'  key='%s' \n", inString, key);
    printf(" xxx multiplier = %le \n", *multiplier );
    debugexit(fnam);
  }

  return ;

} // end parse_multiplier

// ===========================================
void init_GENPOLY(GENPOLY_DEF *GENPOLY) {
  int o;
  GENPOLY->ORDER = -9;
  for(o=0; o < MXORDER_GENPOLY; o++ ) {
    GENPOLY->COEFF_RANGE[o][0] = 0.0 ;
    GENPOLY->COEFF_RANGE[o][1] = 0.0 ;
  }
} // end init_GENPOLY

void parse_GENPOLY(char *stringPoly, char *varName, 
		   GENPOLY_DEF *GENPOLY, char *callFun) {

  // Mar 23 2019
  // parse input stringPoly and load GENPOLY structure with 
  // polynomial of arbitrary order (up to 20).
  //
  // Examples:
  //  stringPoly = "1.0,.034,1.0E-4" 
  //     --> load 2nd order polynominal
  //  stringPoly = "0.9:1.1,.034,1.0E-4,2.2E-8" 
  //     --> load 3rd order polynominal, and a0 coefficient
  //         range is 0.9 to 1.1 for random selection.
  //
  // Input *callFun is for error or debug msg only.
  //
  // Be carefull to parse both commas and colons
  //
  // Mar 22 2020: pass varName and load it.
  //

  int MEMC    = sizeof(char);
  int MXSPLIT = MXORDER_GENPOLY;
  int LDMP    = 0 ;

  int o, NSPLIT, NRANGE, ORDER ;
  double DVAL0, DVAL1 ;
  char *splitValue[MXORDER_GENPOLY], *splitRange[2];
  char *tmpVal, colon[] = ":", comma[] = "," ;
  char fnam[] = "parse_GENPOLY" ;

  // ----------- BEGIN ------------

  sprintf(GENPOLY->STRING, "%s", stringPoly );
  sprintf(GENPOLY->VARNAME,"%s", varName ); // Mar 22 2020

  for(o=0; o < MXSPLIT; o++ ) 
    { splitValue[o] = (char*)malloc( 40*MEMC); }

  splitRange[0] = (char*)malloc( 40*MEMC);
  splitRange[1] = (char*)malloc( 40*MEMC);

  splitString(stringPoly, comma, MXSPLIT,      // inputs
              &NSPLIT, splitValue );      // outputs         

  ORDER = NSPLIT-1;
  GENPOLY->ORDER = ORDER ;

  // loop over each poly coeff and check for colon separating range
  for(o=0; o <= ORDER ; o++ ) {
    tmpVal = splitValue[o] ;
    if ( strstr(tmpVal,colon) == NULL ) {
      sscanf(tmpVal, "%le", &DVAL0 ); DVAL1=DVAL0;
    }
    else {
      splitString(tmpVal, colon, 4,
		  &NRANGE, splitRange );  // outputs
      if ( NRANGE != 2 ) {
	sprintf(c1err,"NRANGE=%d for order=%d. Expect 2 args"
		" around one colon.", 	NRANGE, o);
	sprintf(c2err,"Check %s element of %s (callFun=%s)", 
		tmpVal, stringPoly, callFun );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }
      sscanf(splitRange[0], "%le", &DVAL0 ); 
      sscanf(splitRange[1], "%le", &DVAL1 ); 
    }
    
    if ( DVAL1 < DVAL0 ) {
      sprintf(c1err,"Invalid range for poly-order %d (callFun=%s)", 
	      o, callFun);
      sprintf(c2err,"%s has 2nd value smaller than 1st.", tmpVal);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
    }

    GENPOLY->COEFF_RANGE[o][0] = DVAL0;
    GENPOLY->COEFF_RANGE[o][1] = DVAL1;

  } // end loop over poly terms

  for(o=0; o < MXSPLIT; o++ )  { free(splitValue[o]); }
  free(splitRange[0]);
  free(splitRange[1]);

  if ( LDMP ) {
    printf("\n xxx ======== %s DUMP =========== \n", fnam );
    printf(" xxx calling function: %s \n", callFun);
    printf(" xxx input stringPoly = '%s' \n", stringPoly);
    printf(" xxx NORDER = %d \n", ORDER );
    for(o=0; o <= ORDER; o++ ) {
      DVAL0 = GENPOLY->COEFF_RANGE[o][0];
      DVAL1 = GENPOLY->COEFF_RANGE[o][1];
      if ( DVAL1 > DVAL0 ) 
	{ printf(" xxx A%d : %le to %le \n", o, DVAL0, DVAL1); }
      else
	{ printf(" xxx A%d : %le \n", o, DVAL0 ); }
      fflush(stdout);
    }
    
    //    debugexit(fnam);
  } // end LDMP

  return ;

} // end parse_GENPOLY

void copy_GENPOLY(GENPOLY_DEF *GENPOLY_IN, GENPOLY_DEF *GENPOLY_OUT) {

  // Created April 5 2021
  // For input GENPOLY_IN, copy it to GENPOLY_OUT.

  int ORDER = GENPOLY_IN->ORDER ;
  int o;
  // ------- BEGIN ---------

  sprintf(GENPOLY_OUT->STRING,  "%s", GENPOLY_IN->STRING);
  sprintf(GENPOLY_OUT->VARNAME, "%s", GENPOLY_IN->VARNAME);
  GENPOLY_OUT->ORDER = GENPOLY_IN->ORDER ;

  for(o=0; o <= ORDER; o++ ) {
    GENPOLY_OUT->COEFF_RANGE[o][0] = GENPOLY_IN->COEFF_RANGE[o][0];
    GENPOLY_OUT->COEFF_RANGE[o][1] = GENPOLY_IN->COEFF_RANGE[o][1];
  }

  return;
} // end copy_GENPOLY


double eval_GENPOLY(double VAL, GENPOLY_DEF *GENPOLY, char *callFun) {

  // Mar 2019
  // evaluate polynominal for input value 'val'.
  // Note that random numbers are used for ranges.
  // here we avoid using slow 'pow' function.

  int o, ORDER = GENPOLY->ORDER;
  double VALPOLY = 0.0 ;
  double VALPOW  = 1.0 ;
  double COEFF_RANGE[2], RANCOEFF;
  char fnam[] = "eval_GENPOLY";
  int LDMP = 0 ;
  // ---------- begin ------------

  for(o=0; o <= ORDER; o++ ) {
    COEFF_RANGE[0] = GENPOLY->COEFF_RANGE[o][0];
    COEFF_RANGE[1] = GENPOLY->COEFF_RANGE[o][1];

    if ( COEFF_RANGE[0] < COEFF_RANGE[1] ) 
      { RANCOEFF = getRan_Flat ( 2, COEFF_RANGE ) ; }
    else
      { RANCOEFF = COEFF_RANGE[0]; }

    VALPOLY += RANCOEFF * VALPOW ;
    VALPOW = VALPOW * VAL;
  }

  if ( ORDER < 0 ) {
    sprintf(c1err,"Undefined POLYNOMIAL passed from fun=%s", callFun);
    sprintf(c2err,"VAL=%le  POLYSTRING='%s' ", VAL, GENPOLY->STRING);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  if ( LDMP ) {
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx %s DUMP for VAL=%le  and POLY='%s\n", 
	   fnam, VAL, GENPOLY->STRING );
    printf(" xxx ORDER=%d   VALPOLY = %le \n", ORDER, VALPOLY );
  }

  return(VALPOLY);

} // end eval_GENPOLY

// ====================================
const gsl_rng_type *T_Poisson ;
gsl_rng *r_Poisson ;

int getRan_Poisson(double mean){

  // Created DEc 27 2017
  // Return random Poisson integer from input mean.


  unsigned int k=0 ;

  // change random seed with gsl_rgn_set ??

  if ( mean < 1.0E-12 ) {
    gsl_rng_env_setup();
    T_Poisson = gsl_rng_default;
    r_Poisson = gsl_rng_alloc(T_Poisson);   
    return(0);
  }

  k = gsl_ran_poisson (r_Poisson,mean);
  
  //  gsl_rng_free(r);

  return( (int)k );

} // end getRan_Poisson

void get_SNANA_VERSION(char *snana_version) // pass global declaration
{ sprintf(snana_version, "%s", SNANA_VERSION_CURRENT); } 
void get_snana_version__(char *snana_version) 
{  get_SNANA_VERSION(snana_version); }


float get_SNANA_VERSION_FLOAT(char *snana_version) {
  // Oct 26 2020
  // Convert *snana_version string to float.
  // e.g. *snana_version = v10_78c -> return 10.78
  double dval0 = atof(&snana_version[1]) ;
  double dval1 = atof(&snana_version[4]) ;
  float  fval = (float)(dval0 + dval1/100.0) ;
  return(fval);
} // end get_SNANA_VERSION_FLOAT

float get_snana_version_float__(char *snana_version)
{ return get_SNANA_VERSION_FLOAT(snana_version); }


bool correct_sign_vpec_data(char *snana_version_data) {

  // Jan 2021: if SNANA_VERSION key is not known, assume vpec sign is correct
  if ( strcmp(snana_version_data,"UNKNOWN") == 0 ) { return true; }


  float version_f = get_SNANA_VERSION_FLOAT(snana_version_data);
  if ( version_f < 11.02 )
    { return false; } // VPEC in FITS data has wrong sign convention
  else
    { return true; }  // has correct convention

} // end correct_sign_vpec_data

bool correct_sign_vpec_data__(char *snana_version_data) 
{ return correct_sign_vpec_data(snana_version_data); }


void INIT_SNANA_DUMP(char *STRING) {

  // Created Sep 21 2017
  // Parse input *STRING and load DUMP_STRING_INFO structure.
  // This string can be set in the analysis programs via 
  // &SNLCINP input  DUMP_STRING = 'BLA BLA'. For simulation ... (not yet).
  //
  // Example:
  // DUMP_STRING = 
  //  'FUN FLUXERRCALC  FILTERS iz  CIDLIST 5001,5002
  //      MJDRANGE 53662 53663 ABORT'
  //
  // Will trigger the dump in subroutine FLUXERRCALC for filters
  // i & z, CIDs 5001 & 5002, and within MJDRANGE given above.
  // CIDLIST is interpreated as a list of comma-separated strings.
  // The ABORT flag triggers an abort after all CIDs and FILTERS
  // have been processed.  The job will finish to completion if:
  // a) no ABORT flag, or b) CIDLIST includes an invalid CID, or
  // c) FILTERS includes an invalid filter. For completed jobs,
  // search the stdout for dump info.
  //
  // Default values:
  //   + CIDLIST   --> none
  //   + FILTERS   --> none
  //   + MJDRANGE  --> wide open (0 to 9999999)
  //
  //
  // List of functions integrated with this utility:
  //      file          function         date
  //     ----------------------------------------------
  //    snana.car     FLUXERRCALC     2017-09-21
  //        ... more later ...
  //

  int i;
  char fnam[] = "INIT_SNANA_DUMP" ;

  // --------------- BEGIN --------------

  DUMP_STRING_INFO.NCID      = 0 ;
  DUMP_STRING_INFO.NFILT     = 0 ;
  DUMP_STRING_INFO.ABORTFLAG = 0 ;
  for(i=0; i < MXCID_DUMP; i++ ) 
    { DUMP_STRING_INFO.NFILT_DONE[i] = 0 ; }

  DUMP_STRING_INFO.FUNNAME[0]  = 0 ;
  DUMP_STRING_INFO.FILTLIST[0] = 0 ;

  DUMP_STRING_INFO.MJDRANGE[0] = 0.0 ;
  DUMP_STRING_INFO.MJDRANGE[1] = 9999999. ;

  // bail if required key is not present
  if ( strstr(STRING,"FUN") == NULL ) { return ; }

  // - - - - - - - - - - - - - - - - 
#define MXWORD_DUMP_STRING 20
  char BLANK[] = " " ;
  int MXSPLIT = MXWORD_DUMP_STRING ;
  int Nsplit, iwd, NCID=0 ;
  char *ptrSplit[MXWORD_DUMP_STRING];
  char *wordList[MXWORD_DUMP_STRING];
  int  MEMC = 60 * sizeof(char);
  // allocate wordList and assign ptrSplit
  for(i=0; i < MXWORD_DUMP_STRING; i++) {
    wordList[i] = (char*) malloc ( MEMC ) ;
    ptrSplit[i] = wordList[i];
  }
  splitString(STRING, BLANK, MXSPLIT,
	      &Nsplit, ptrSplit) ;   // returned

  char *cwd0=wordList[0], *cwd1=wordList[0], *cwd2=wordList[0];
  for(iwd=0; iwd < Nsplit; iwd++ ) {
    cwd0 = wordList[iwd] ;
    if ( iwd < (Nsplit-1) ) { cwd1 = wordList[iwd+1] ; }
    if ( iwd < (Nsplit-2) ) { cwd2 = wordList[iwd+2] ; }

    if ( strcmp(cwd0,"FUN") == 0 ) 
      { sscanf(cwd1, "%s", DUMP_STRING_INFO.FUNNAME );  }

    if ( strcmp(cwd0,"FILTLIST") == 0 || strcmp(cwd0,"FILTERS")==0 )  { 
      sscanf(cwd1, "%s", DUMP_STRING_INFO.FILTLIST );  
      DUMP_STRING_INFO.NFILT = strlen(DUMP_STRING_INFO.FILTLIST);
    }

    if ( strcmp(cwd0,"MJDRANGE") ==0 ) {
      sscanf(cwd1, "%le", &DUMP_STRING_INFO.MJDRANGE[0] );  
      sscanf(cwd2, "%le", &DUMP_STRING_INFO.MJDRANGE[1] );  
    }

    if ( strcmp(cwd0,"ABORT") ==0 ) 
      { DUMP_STRING_INFO.ABORTFLAG=1; }  // no argument

    if ( strcmp(cwd0,"CID") ==0 || strcmp(cwd0,"CIDLIST")==0 ) { 
      char *ptrCCID[MXCID_DUMP],  comma[] = "," ;
      for(i=0;i<MXCID_DUMP;i++) {ptrCCID[i] = DUMP_STRING_INFO.CCIDLIST[i];}
      splitString(cwd1, comma, MXCID_DUMP, &NCID, ptrCCID) ;  
      DUMP_STRING_INFO.NCID = NCID;
    }

  } // end iwd loop

  // - - - - - -  print summary - - - - - - - - 

  char BANNER[100];
  sprintf(BANNER,"%s: dump instructions", fnam );
  print_banner(BANNER);
  printf("\t Function to Dump: '%s' \n", DUMP_STRING_INFO.FUNNAME  );
  printf("\t Filters  to Dump: '%s' \n", DUMP_STRING_INFO.FILTLIST );

  printf("\t CID list to Dump: ");
  for(i=0; i < NCID; i++ ) { printf("%s ", DUMP_STRING_INFO.CCIDLIST[i] ); }
  printf("\n");

  printf("\t Abort after Dump: %d \n",   DUMP_STRING_INFO.ABORTFLAG );
  printf("\n");

  fflush(stdout);
  //  debugexit(fnam); // xxx REMOVE

  return ;

} // end INIT_SNANA_DUMP

// ==========================================================
int CHECK_SNANA_DUMP(char *FUNNAME, char *CCID, char *BAND, double MJD) {

  int IFLAG = 0 ;
  int NCID  = DUMP_STRING_INFO.NCID ;
  int NFILT = DUMP_STRING_INFO.NFILT ;
  //  char fnam[] = "CHECK_SNANA_DUMP" ;

  // ------------ BEGIN ----------------

  if ( NCID == 0 ) { return(IFLAG); }

  // check function name
  if ( strcmp( DUMP_STRING_INFO.FUNNAME,FUNNAME) != 0 ) 
    { return(IFLAG); }

  // check valid band
  if ( strstr(DUMP_STRING_INFO.FILTLIST,BAND) == NULL ) 
    { return(IFLAG); } 

  // check valid CCID
  int i, ICID=-9;
  for(i=0; i < NCID; i++ ) {
    if ( strcmp(DUMP_STRING_INFO.CCIDLIST[i],CCID)==0 ) { ICID=i; }
  }
  if ( ICID < 0 ) { return(IFLAG); }

  // check valid MJD
  if ( MJD < DUMP_STRING_INFO.MJDRANGE[0] ) { return(IFLAG); }
  if ( MJD > DUMP_STRING_INFO.MJDRANGE[1] ) { return(IFLAG); }

  // - - - - - - - - - - -
  DUMP_STRING_INFO.NFILT_DONE[ICID]++  ;

  // if we get here, set IFLAG > 0
  IFLAG = DUMPFLAG_noABORT ; // default

  // if ABORTFLAG is set, abort only when all CID and FILTERS have
  // been dumped.
  if ( DUMP_STRING_INFO.ABORTFLAG ) {
    int NFILT_DONE, NOTDONE = 0 ;
    for(i=0; i < NCID; i++ ) {
      NFILT_DONE = DUMP_STRING_INFO.NFILT_DONE[i];
      if ( NFILT_DONE < NFILT ) { NOTDONE=1; }
    }
    if ( NOTDONE==0 ) { IFLAG = DUMPFLAG_withABORT; }
  }

  return(IFLAG);

} // end CHECK_SNANA_DUMP

void ABORT_SNANA_DUMP(void) {

  // user should call this routine if CHECK_SNANA_DUMP returns 2,
  // but call this after doing the last dump.
  // The dump here is a friendly dump, not Evil-Abort face dump.

  printf("\n Done with all user-requested DUMPs for: \n");
  printf("\t %d CIDs, FILTERS=%s, MJDRANGE=%.3f-%.3f\n",
	 DUMP_STRING_INFO.NCID, DUMP_STRING_INFO.FILTLIST,
	 DUMP_STRING_INFO.MJDRANGE[0], DUMP_STRING_INFO.MJDRANGE[1] 
	 );
  printf(" User requests ABORT after all DUMPs. Bye Bye. \n");
  fflush(stdout);
  exit(66);

} // end ABORT_SNANA_DUMP


// - - - - - - - -
void init_snana_dump__(char *STRING) { INIT_SNANA_DUMP(STRING); }

int  check_snana_dump__(char *FUNNAME, char *CCID, char *BAND, double *MJD) 
{  return ( CHECK_SNANA_DUMP(FUNNAME,CCID,BAND,*MJD) ) ; } 

void abort_snana_dump__(void) { ABORT_SNANA_DUMP(); }


// ==========================================
void read_SURVEYDEF(void) {

  // Created Aug 19 2016 [moved frmo SALT2mu.c on July 2017]
  // Read SURVEY.DEF file and store each defined survey name in
  // SURVEY_INFO.SURVEYDEF_LIST[IDSURVEY].
  //

  char SURVEYDEF_FILE[MXPATHLEN], c_get[60], nameTmp[40];
  int  idTmp ;
  FILE *fp ;
  char fnam[] = "read_SURVEYDEF" ;

  // ------------ BEGIN -------------

  sprintf(SURVEYDEF_FILE,"%s/SURVEY.DEF", PATH_SNDATA_ROOT);
  fp = fopen(SURVEYDEF_FILE,"rt");
  if ( !fp ) {
    sprintf(c1err,"Could not open SURVEY.DEF file");
    sprintf(c2err,"%s", SURVEYDEF_FILE);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // init each survey name to NULL
  for(idTmp=0; idTmp < MXIDSURVEY; idTmp++ ) { 
    sprintf(SURVEY_INFO.SURVEYDEF_LIST[idTmp],"NULL"); 
  }

  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"SURVEY:") == 0 ) {
      readchar(fp, nameTmp ); readint(fp, 1, &idTmp );
      if ( idTmp >= MXIDSURVEY ) {
	sprintf(c1err,"IDSURVEY=%d(%s) exceeds MXIDSURVEY=%d", 
		idTmp, nameTmp, MXIDSURVEY);
	sprintf(c2err,"check %s", SURVEYDEF_FILE);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
      }

      sprintf(SURVEY_INFO.SURVEYDEF_LIST[idTmp],"%s", nameTmp);
    }
  }

  fclose(fp);
  return ;

} // end read_SURVEYDEF

// ==========================================                                   
int get_IDSURVEY(char *SURVEY) {
  // return integer IDSURVEY for input string *SURVEY
  int ID;
  //  char fnam[] = "get_IDSURVEY" ;     
  // -------------- BEGIN -------------  
  for(ID=0; ID < MXIDSURVEY; ID++ ) {
    if ( strcmp(SURVEY_INFO.SURVEYDEF_LIST[ID],SURVEY) == 0 )
      { return ID; }
  }

  return(-9);
} // end get_IDSURVEY 

// ============================================
int  exec_cidmask(int mode, int CID) {

  // Created Apri 2017 (code transfered from snana.car)
  // utility to store integer CIDs in bit mask (to save memory),
  // and to check if any CID-bit is set.  Used by simulation
  // and analysis to flag duplicates.

  // Inputs:
  //   mode =  0 --> allocate memory with MXCID=CID
  //   mode =  1 --> set CID'th bit
  //   mode =  2 --> return 1 if CID'th bit is set; 0 otherwise 
  //   mode = -1 --> free memory
  //
  // Apr 25 2017: bug fix, fmodf -> fmod for double

  int J, JBIT, NINT ;
  double dCID = (double)CID;
  char fnam[] = "exec_cidmask" ;

  // --------------- BEGIN -----------------

  if ( mode == 0 ) {
    NINT = (CID+1)/32 + 2 ;
    float xMB  = (float)(NINT*4)/1.0E6 ;
    printf("   Allocate %.2f MB for CIDMASK array "
	   "(to check duplicates)\n", xMB);  fflush(stdout);
    CIDMASK_LIST = (unsigned int*) malloc ( NINT * sizeof(unsigned int) ) ;
    MXCIDMASK = CID;
    NCIDMASK_LIST = NINT;
    for(J=0; J < NINT; J++ ) { CIDMASK_LIST[J]=0; }
    return(1);
  }
  else if ( mode == 1 ) {
    if ( CID > MXCIDMASK ) { return(0); }
    J = (CID/32) + 1 ;    JBIT = (int)fmod(dCID,32.0) ;
    CIDMASK_LIST[J] |= (1 << JBIT) ;
    /*
    printf(" xxx %s: CID=%2d --> J=%d JBIT=%2d  CIDMASK=%u,%u \n",
	   fnam, CID, J, JBIT, CIDMASK_LIST[1], CIDMASK_LIST[2] ); */
    return(1);
  }
  else if ( mode == 2 ) {
    J = (CID/32) + 1 ;    JBIT = (int)fmod(dCID,32.0) ;
    if ( (CIDMASK_LIST[J] & (1 << JBIT)) > 0 )
      { return(1); }
    else
      { return(0); }
  }
  else if ( mode == -1 ) {
    free(CIDMASK_LIST);
    return(1);
  }
  else {
    sprintf(c1err,"Invalid mode=%d for CID=%d", mode, CID);
    c2err[0] = 0 ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  return(1);

} // end exec_cidmask

int  exec_cidmask__(int *mode, int *CID) {  
  int istat = exec_cidmask(*mode,*CID); 
  return(istat);
}

void test_cidmask(void) {

  // test exec_cidmask function
#define NLIST 5
  int CIDLIST[NLIST] = { 4, 22, 23, 38, 78 } ;
  int i, CID, LSET ;
  int MXCID = 10000;
  char fnam[] = "test_cidmask";

  // ---------- BEGIN ----------
  
  exec_cidmask(0,MXCID);
  for(i=0; i < NLIST; i++ ) {
    CID = CIDLIST[i];
    exec_cidmask(1,CID);
  }

  for(CID=0; CID<100 ; CID++ ) {
    LSET = exec_cidmask(2,CID);
    printf(" xxx CID = %3d --> LSET=%d \n", CID, LSET);
    
  }

  fflush(stdout);
  debugexit(fnam);

  return ;

} // end test_cidmask


void extract_MODELNAME(char *STRING, char *MODELPATH, char *MODELNAME) {

  // Created Feb 22 2016
  // For input STRING = BLA
  //    return MODELPATH='' and MODELNAME=BLA
  // For input STRING = /project/models/SALT2/SALT2.ABC
  //    return MODELPATH=STRING and MODELNAME=SALT2.ABC
  //


  int i,i2, lastSlash, LENSTR, ENVstat ;
  char fnam[] = "extract_MODELNAME" ;

  // ------------ BEGIN -----------

  MODELPATH[0] = MODELNAME[0] = 0 ;

  ENVstat   = ENVreplace(STRING,fnam,1);    // check for ENV
  LENSTR    = strlen(STRING);
  lastSlash = -9;
  for(i=0; i < LENSTR; i++ ) {
    if ( STRING[i] == '/' ) { lastSlash=i; }
  }

  if ( lastSlash < 0 ) {
    sprintf(MODELNAME, "%s", STRING);
  }
  else {
    sprintf(MODELPATH, "%s", STRING);

    for(i=lastSlash+1; i < LENSTR; i++ ) {
      i2 = i - lastSlash - 1; 
      sprintf(&MODELNAME[i2], "%c", STRING[i] ); 
    } // end i loop over string chars

  } // end lastSlah if

  return ;

} // end extract_MODELNAME

void extract_modelname__(char *STRING, char *MODELPATH, char *MODELNAME) {
  extract_MODELNAME(STRING, MODELPATH, MODELNAME);
}

// ************************************
void parse_commaSepList(char *item_name, char *item, int MAX_ITEM, int MXCHAR,
			int *n_item, char ***arrayList ) {

  // Created Oct 2020
  // Utility to split comma-sep string *item and return number of 
  // items (n_item) and char arrayList[i]
  //  Inputs:
  //    *item_name : descriptor (for error message)
  //    *item      : comma-sep list
  //    MAX_ITEM   : abort if more than this many items
  //    MXCHAR     : allocate MXCHAR length per item
  //
  //  Output
  //    *n_item   : number of items in *item
  //    arrayList : array of individual items
  //
  // - - - -
  int i;
  int MEMC = MXCHAR * sizeof(char);
  char fnam[] = "parse_commaSepList" ;

  // ---------- BEGIN --------------
  // first allocate memory for file names 
  *arrayList = (char**)malloc( MAX_ITEM*sizeof(char*));
  for(i=0; i < MAX_ITEM; i++ ) 
    { (*arrayList)[i] = (char*)malloc(MEMC); }

  // bail on blank string
  if ( strlen(item) == 0 ) { *n_item=0; return; }

  // split item string
  splitString(item, COMMA, MAX_ITEM,    // inputs
	      n_item, *arrayList );      // outputs 
  
  char *f0 = *arrayList[0];
  if ( IGNOREFILE(f0) ) { *n_item = 0 ; }

} // end parse_commaSepList


// *****************************************************
void FILTER_REMAP_INIT(char *remapString, char *VALID_FILTERLIST,
		       int *NFILT_REMAP,
                       int *IFILTLIST_REMAP, char *FILTLIST_REMAP) {

  // Created Feb 2017
  // Example:
  //   Input remapString = 'abc->U def->B gh->V jklm->R'
  //   
  // Outputs
  //   NFILT_REMAP = 4
  //   IFILTLIST_REMAP = array of 4 indices corresponding to UBVR
  //   FILTLIST_REMAP  = 'UBVR'
  //
  // Inputs  VALID_FILTERLIST is a string of valid input filters
  // to check that the remapString does not contain invalid filters.
  //

  int  NFILT_ORIG, NFILT, ifilt, lenStr;
  int  USEFILT_ORIG[MXFILTINDX], USEFILT_REMAP[MXFILTINDX];
  char *localString, *ptrtok, ctmp[80] ;
  char fnam[] = "FILTER_REMAP_INIT" ;

  // ------------ BEGIN ------------

  // init output function args
  *NFILT_REMAP =  FILTLIST_REMAP[0]  = 0 ;

  // init global
  FILTER_REMAP.NMAP  = 0 ;

  lenStr = strlen(remapString);
  if ( lenStr == 0 ) { return ; }

  // - - - - -
  set_FILTERSTRING(FILTERSTRING);

  localString = (char*) malloc( lenStr* sizeof(char) );
  for(ifilt=0; ifilt<MXFILTINDX; ifilt++) { 
    USEFILT_ORIG[ifilt]  = 0; 
    USEFILT_REMAP[ifilt]  = 0; 

    FILTER_REMAP.IFILTOBS_MAP[ifilt] = -9 ;
    FILTER_REMAP.IFILT_MAP[ifilt]    = -9 ;
  }

  printf("\n %s: \n", fnam);
  printf(" remapString = \n  '%s' \n", remapString);


  // split remapString into pieces, where each piece is a map
  sprintf(localString, "%s", remapString); 
  ptrtok = strtok(localString," ") ; // split string
  NFILT  = NFILT_ORIG = 0;
  while ( ptrtok != NULL ) {
      sprintf(ctmp, "%s", ptrtok );
      if ( NFILT < MXFILT_REMAP )
	{  sprintf(FILTER_REMAP.MAPSTRING[NFILT], "%s", ctmp); }
      NFILT++ ;	
      ptrtok = strtok(NULL, " ");
  }
  free(localString) ;
  sprintf(c2err,"See remapString"); // generic abort comment.

  if ( NFILT >= MXFILT_REMAP ) {
    sprintf(c1err,"NFILT_REMAP=%d exeeds bound MXFILT_REMAP=%d.",
		NFILT, MXFILT_REMAP);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  // analyze each mapstring and load IFILTOBS_MAP
  char *ptrMap, BAND[2], band[2] ;  
  int lenMap, i, IFILTOBS, IFILT_REMAP, ifiltobs;
  for(IFILT_REMAP=0; IFILT_REMAP < NFILT; IFILT_REMAP++ ) {  
    ptrMap = FILTER_REMAP.MAPSTRING[IFILT_REMAP] ;
    lenMap = strlen(ptrMap);

    // start with last char, which is the mapped ifiltobs
    sprintf(BAND, "%c", ptrMap[lenMap-1] ) ;
    IFILTOBS     = INTFILTER(BAND);  // mapped filter
    FILTLIST_REMAP         = strcat(FILTLIST_REMAP,BAND);
    IFILTLIST_REMAP[IFILT_REMAP] = IFILTOBS ;
    // abort if remapped filter is used more than once
    if ( USEFILT_REMAP[IFILTOBS] ) {
      sprintf(c1err,"Remapped filter = %s used more than once.", BAND);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }
    USEFILT_REMAP[IFILTOBS] += 1;

    for(i=0; i < lenMap-3 ; i++ ) {
      sprintf(band, "%c", ptrMap[i] );
      ifiltobs = INTFILTER(band); // original filter

      // abort if original filter is used more than once
      if ( USEFILT_ORIG[ifiltobs] ) {
	sprintf(c1err,"original filter = %s used more than once.", band);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
      }
      // abort if original filter is not in the VALID_FILTERLIST      
      if ( strstr(VALID_FILTERLIST,band) == NULL ) {
	print_preAbort_banner(fnam);
	printf("  Defined bands: '%s'\n", VALID_FILTERLIST );
	sprintf(c1err,"Invalid filter=%s is not defined ", band);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err);  
      }

      USEFILT_ORIG[ifiltobs] += 1;
      FILTER_REMAP.IFILTOBS_MAP[ifiltobs] = IFILTOBS;
      FILTER_REMAP.IFILT_MAP[ifiltobs]    = IFILT_REMAP ;

      NFILT_ORIG++ ;
      printf("\t %s: %s->%s (%2.2d->%2.2d)  IFILT_REMAP=%d\n",
	     ptrMap, band, BAND, ifiltobs, IFILTOBS, IFILT_REMAP);
      fflush(stdout);

    }  // i loop
  } // ifilt
 

  *NFILT_REMAP = NFILT ;
  printf("  Done mapping %d original filters to %d mapped filters.\n",
	 NFILT_ORIG, NFILT);
  fflush(stdout);

  return;

} // end FILTER_REMAP_INIT

void filter_remap_init__(char *remapString, char *VALID_FILTERLIST,
			 int *NFILT_REMAP,
                         int *IFILTLIST_REMAP, char *FILTLIST_REMAP) {
  FILTER_REMAP_INIT(remapString, VALID_FILTERLIST,
		    NFILT_REMAP, IFILTLIST_REMAP, FILTLIST_REMAP);
}

void FILTER_REMAP_FETCH(int IFILTOBS_ORIG,
			int *IFILTOBS_REMAP, int *IFILT_REMAP) {

  // Created Feb 2017
  // For input IFILTOBS_ORIG, returns 
  //   + absolute IFILTOBS_REMAP,
  //   + sparse   IFILT_REMAP indices.

  //  char fnam[] = "FILTER_REMAP_FETCH" ;

  // ------------- BEGIN -----------
  
  *IFILT_REMAP    = FILTER_REMAP.IFILT_MAP[IFILTOBS_ORIG] ;
  *IFILTOBS_REMAP = FILTER_REMAP.IFILTOBS_MAP[IFILTOBS_ORIG]  ;

  return ;

} // end FILTER_REMAP_FETCH

void filter_remap_fetch__(int *IFILTOBS_ORIG,
			  int *IFILTOBS_REMAP, int *IFILT_REMAP) {
  FILTER_REMAP_FETCH(*IFILTOBS_ORIG, IFILTOBS_REMAP, IFILT_REMAP) ;
}



// =============================================
void init_GaussIntegral(void) {

  int ix, NX;
  double xmin = 0.0;
  double xmax, xbin, xi;
  double XMAX = 10.0 ;  // store integrals up to 10 sigma
  GAUSS_INTEGRAL_STORAGE.XMAX = XMAX ;

  char fnam[] = "init_GaussIntegral";

  // --------- BEGIN -----------

  GAUSS_INTEGRAL_STORAGE.INIT_FLAG = 1;

  
  NX = 500 ;
  GAUSS_INTEGRAL_STORAGE.NBIN_XMAX = NX ;

  xbin = XMAX / (double)NX ;

  if ( NX <= 0 || NX >= MXSTORE_RAN ) {
    sprintf(c1err,"Invalid NBIN_XMAX=%d", NX);
    sprintf(c2err,"Valid range is < %d", MXSTORE_RAN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }

  for(ix=0; ix < NX; ix++ ) {
    xi   = (double)ix;
    xmax = xmin + xbin * (xi - 0.5) ;
    GAUSS_INTEGRAL_STORAGE.XMAX_LIST[ix] = xmax ;
    GAUSS_INTEGRAL_STORAGE.GINT_LIST[ix] = GaussIntegral(xmin,xmax);   
  }

  // done with init
  GAUSS_INTEGRAL_STORAGE.INIT_FLAG = 0;


#define NDUMP 5
  int LDMP=0, idump;
  if ( LDMP ) {
    double xmin_dump[NDUMP] = { 0.0, -1.0, -1.0, -2.0, 1.0 } ;
    double xmax_dump[NDUMP] = { 1.0,  1.0,  0.0, -1.0, 2.0 } ;
    for(idump=0 ; idump < NDUMP ; idump++ ) {
      xmin= xmin_dump[idump];  xmax=xmax_dump[idump];
      printf(" xxx DUMP %s(%4.1f,%4.1f) = %f \n",  
	     fnam, xmin, xmax, GaussIntegral(xmin,xmax) );
    }
    debugexit(fnam);
  }

  return ;

} // end init_GaussIntegral

// =============================================
double GaussIntegral(double nsig1, double nsig2) {

  // Created Oct 2016
  // return Gaussian integral between nsig1 and nsig2
  //
  // Mar 19 2019: if nsig[1,2] > 10, set GINT=0 to avoid abort.

  int    NBIN, ibin ;
  double SUM, DENOM, SIGBIN, xsig, xi, arg ;
  char fnam[] = "GaussIntegral" ;

  // --------- BEGIN ----------

  SUM    = 0.0 ;

  if ( nsig1 == nsig2 ) { return SUM; }

  if ( GAUSS_INTEGRAL_STORAGE.INIT_FLAG  ) {
    DENOM = sqrt(TWOPI);
    // brute-force integral during init stage
    NBIN  = 1000*(int)(nsig2-nsig1);
    if ( NBIN < 5 ) { NBIN=5; }
    SIGBIN = (nsig2-nsig1)/(double)NBIN;
    for(ibin=0; ibin < NBIN; ibin++ ) {
      xi     = (double)ibin ;
      xsig   = nsig1 + SIGBIN*(xi+0.5);
      arg    = (xsig*xsig)/2.0 ;
      SUM   += exp(-arg);
    }
    SUM *= SIGBIN/DENOM ;
  }

  else {
    // interpolate stored integrals from 0 to xmax
    int    OPT_INTERP=1;
    double GINT1, GINT2, x1, x2, xsign1=1.0, xsign2=1.0 ;
    double *ptr_xmax = GAUSS_INTEGRAL_STORAGE.XMAX_LIST ;
    double *ptr_gint = GAUSS_INTEGRAL_STORAGE.GINT_LIST ;
    double  NSIGMAX  = 0.995*GAUSS_INTEGRAL_STORAGE.XMAX ;

    x1 = fabs(nsig1);  if(nsig1!=0.0) { xsign1 = x1/nsig1; }
    x2 = fabs(nsig2);  if(nsig2!=0.0) { xsign2 = x2/nsig2; }
    NBIN = GAUSS_INTEGRAL_STORAGE.NBIN_XMAX ;

    if ( NBIN <= 0 || NBIN >= MXSTORE_RAN ) {
      sprintf(c1err,"Invalid NBIN_XMAX=%d", NBIN);
      sprintf(c2err,"Check if init_GaussIntegral was called");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
    }


    if (x1 > NSIGMAX ) { x1 = NSIGMAX; }
    if (x2 > NSIGMAX ) { x2 = NSIGMAX; }

    GINT1 = interp_1DFUN(OPT_INTERP,x1,NBIN, ptr_xmax, ptr_gint, fnam);
    GINT2 = interp_1DFUN(OPT_INTERP,x2,NBIN, ptr_xmax, ptr_gint, fnam);

    SUM   = xsign2*GINT2 - xsign1*GINT1;
  }

  return(SUM) ;

} // end GaussIntegral

// =============================================
double angSep( double RA1,double DEC1, 
	       double RA2,double DEC2, double  scale) {

  // Copied from DIFFIMG on Nov 16 2015
  //
  // Oct 9, 2012 R. Kessler
  // for input coords of point 1 (RA1,DEC1) and point 2 (RA2,DEC2),
  // return angular separation. Inputs are in degrees and output
  // is in degrees x scale ->
  // * scale = 1    -> output is in degrees
  // * scale = 60   -> output is in arcmin
  // * scale = 3600 -> output is in arcsec

  double X1,Y1,Z1, X2, Y2, Z2, DOTPROD, sep ;
  double RAD = RADIAN ;

  // ------------- BEGIN ------------------

  X1 = cos(RA1*RAD) * cos(DEC1*RAD);
  Y1 = sin(RA1*RAD) * cos(DEC1*RAD);
  Z1 = sin(DEC1*RAD);

  X2 = cos(RA2*RAD) * cos(DEC2*RAD);
  Y2 = sin(RA2*RAD) * cos(DEC2*RAD);
  Z2 = sin(DEC2*RAD);

  DOTPROD = (1.0-1.0E-15)*(X1*X2 + Y1*Y2 + Z1*Z2);
  sep = acos(DOTPROD)/RAD ; // angular sep, degrees

  return (sep * scale) ;

} // end of angSep


// ==============================================================
// Altered from the fortran SLALIB by David Cinabro, June 2006.
// Translates equatorial coordinats (RA,DEC) to galactic
// longitude and latitude.  All in degrees and double precision.
// All the subroutines needed are included below.
// Usage: 
//    slaEqgal ( double RA, double DEC, double *GalLat, double *GalLong );
// ==============================================================

void slaEqgal ( double dr, double dd, double *dl, double *db )
/*
**  - - - - - - - - -
**   s l a E q g a l
**  - - - - - - - - -
**
**  Transformation from J2000.0 equatorial coordinates to
**  IAU 1958 Galactic coordinates.
**
**  (double precision)
**
**  Given:
**     dr,dd       double       J2000.0 RA,Dec
**
**  Returned:
**     *dl,*db     double       Galactic longitude and latitude l2,b2
**
**  (all arguments were radians, but translation from and to degrees done below)
**
**  Called:
**     slaDcs2c, slaDmxv, slaDcc2s, slaDranrm, slaDrange
**
**  Note:
**     The equatorial coordinates are J2000.0.  Use the routine
**     slaEg50 if conversion from B1950.0 'FK4' coordinates is
**     required.
**
**  Reference:
**     Blaauw et al, Mon.Not.R.astron.Soc.,121,123 (1960)
**
**  Last revision:   21 September 1998
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double v1[3], v2[3];
   double drr, ddr;
   double DPI = 3.1415926535897932384626433832795028841971693993751;

/*
**  l2,b2 system of Galactic coordinates
**
**  p = 192.25       RA of Galactic north pole (mean B1950.0)
**  q =  62.6        inclination of Galactic to mean B1950.0 equator
**  r =  33          longitude of ascending node
**
**  p,q,r are degrees
**
**  Equatorial to Galactic rotation matrix (J2000.0), obtained by
**  applying the standard FK4 to FK5 transformation, for zero proper
**  motion in FK5, to the columns of the B1950 equatorial to
**  Galactic rotation matrix:
*/
   static double rmat[3][3];

   rmat[0][0] = -0.054875539726;
   rmat[0][1] = -0.873437108010;
   rmat[0][2] = -0.483834985808;
   rmat[1][0] =  0.494109453312;
   rmat[1][1] = -0.444829589425;
   rmat[1][2] =  0.746982251810;
   rmat[2][0] = -0.867666135858;
   rmat[2][1] = -0.198076386122;
   rmat[2][2] =  0.455983795705;

   // Translate to radians
   drr = dr*DPI/180.0;
   ddr = dd*DPI/180.0;

/* Spherical to Cartesian */
   slaDcs2c ( drr, ddr, v1 );

/* Equatorial to Galactic */
   slaDmxv ( rmat, v1, v2 );

/* Cartesian to spherical */
   slaDcc2s ( v2, dl, db );

/* Express in conventional ranges */
   *dl = slaDranrm ( *dl );
   *db = slaDrange ( *db );
   // Translate back to degrees
   *dl = *dl*180.0/DPI;
   *db = *db*180.0/DPI;
   return ;
}

void slaDcs2c ( double a, double b, double v[3] )
/*
**  - - - - - - - - -
**   s l a D c s 2 c
**  - - - - - - - - -
**
**  Spherical coordinates to direction cosines (double precision)
**
**  Given:
**     a,b       double      spherical coordinates in radians
**                           (RA,Dec), (long,lat) etc
**
**  Returned:
**     v         double[3]   x,y,z unit vector
**
**  The spherical coordinates are longitude (+ve anticlockwise looking
**  from the +ve latitude pole) and latitude.  The Cartesian coordinates
**  are right handed, with the x axis at zero longitude and latitude,
**  and the z axis at the +ve latitude pole.
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double cosb;
   cosb = cos ( b );
   v[0] = cos ( a ) * cosb;
   v[1] = sin ( a ) * cosb;
   v[2] = sin ( b );
   return ;
}

void slaDmxv ( double dm[3][3], double va[3], double vb[3] )
/*
**  - - - - - - - -
**   s l a D m x v
**  - - - - - - - -
**
**  Performs the 3-d forward unitary transformation:
**     vector vb = matrix dm * vector va
**
**  (double precision)
**
**  Given:
**     dm       double[3][3]    matrix
**     va       double[3]       vector
**
**  Returned:
**     vb       double[3]       result vector
**
**  Note:  va and vb may be the same array.
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   int i, j;
   double w, vw[3];


/* Matrix dm * vector va -> vector vw. */
   for ( j = 0; j < 3; j++ ) {
      w = 0.0;
      for ( i = 0; i < 3; i++ ) {
         w += dm[j][i] * va[i];
      }
      vw[j] = w;
   }

/* Vector vw -> vector vb. */
   for ( j = 0; j < 3; j++ ) {
      vb[j] = vw[j];
   }

   return ;
}

void slaDcc2s ( double v[3], double *a, double *b )
/*
**  - - - - - - - - -
**   s l a D c c 2 s
**  - - - - - - - - -
**
**  Cartesian to spherical coordinates.
**
**  (double precision)
**
**  Given:
**     v       double[3]   x,y,z vector
**
**  Returned:
**     *a,*b   double      spherical coordinates in radians
**
**  The spherical coordinates are longitude (+ve anticlockwise looking
**  from the +ve latitude pole) and latitude.  The Cartesian coordinates
**  are right handed, with the x axis at zero longitude and latitude,
**  and the z axis at the +ve latitude pole.
**
**  If v is null, zero a and b are returned.  At either pole, zero a is
**  returned.
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double x, y, z, r;

   x = v[0];
   y = v[1];
   z = v[2];
   r = sqrt ( x * x + y * y );

   *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
   *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
   return ;
}

double slaDranrm ( double angle )
/*
**  - - - - - - - - - -
**   s l a D r a n r m
**  - - - - - - - - - -
**
**  Normalize angle into range 0-2 pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the range 0-2 pi (double).
**
**  Defined in slamac.h:  D2PI, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double w;
   double D2PI = 6.2831853071795864769252867665590057683943387987502;
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))

   w = dmod ( angle, D2PI );
   return ( w >= 0.0 ) ? w : w + D2PI;
}

double slaDrange ( double angle )
/*
**  - - - - - - - - - -
**   s l a D r a n g e
**  - - - - - - - - - -
**
**  Normalize angle into range +/- pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the range +/- pi.
**
**  Defined in slamac.h:  DPI, D2PI, dmod
**
**  Last revision:   22 July 2004
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
  double w;
  double DPI = 3.1415926535897932384626433832795028841971693993751;
  double D2PI = 6.2831853071795864769252867665590057683943387987502;
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))
/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-(A):(A))

  w = dmod ( angle, D2PI );
  return ( fabs ( w ) < DPI ) ? w : w - dsign ( D2PI, angle );
}

// =======================================
//      end of SLALIB functions
// =======================================

int ENVreplace(char *fileName, char *callFun, int ABORTFLAG) {

  // Feb 2015 [major overhaul Mar 30 2019]
  // if input *fileName starts with $XXX/yyy 
  // then replace $XXX with getenv("XXX").
  // Note that input fileName is modified.
  //
  // Inputs:
  //   *fileName   : file to check for ENV
  //   *callFun    : name of calling function to print on error
  //    ABORTFLAG  : non-zero -> abort on error; else return null fileName
  //
  //  Return SUCCESS or ERROR
  //

  int LL, i, FOUNDSLASH, SEV, NFILE; ;
  char firstChar[2], c1[2], ENVname[MXPATHLEN], suffix[MXPATHLEN] ;
  char fnam[] = "ENVreplace" ;  
  
  // ------------- BEGIN -------------

  if ( strcmp(fileName,"init") == 0 || strcmp(fileName,"INIT") ==0 ) 
    { ENVreplace_store.NFILE = 0 ; return(SUCCESS);  }

  sprintf(firstChar,"%c", fileName[0] );
  suffix[0]=0;
  ENVname[0]=0;

  NFILE = ENVreplace_store.NFILE;
  if ( NFILE < MXFILE_ENVreplace-1 ) 
    { sprintf(ENVreplace_store.FILENAME_ORIG[NFILE],"%s", fileName); }

  if ( *firstChar == '$' ) { 

    FOUNDSLASH = 0 ;
    LL = strlen(fileName);
    for(i=1; i < LL ; i++ ) {
      sprintf(c1,"%c", fileName[i] );
      if( *c1 == '/' ) {  FOUNDSLASH++ ;  }

      if ( FOUNDSLASH > 0  ) 
	{ strcat(suffix,c1); }
      else
	{ sprintf(ENVname, "%s%s", ENVname, c1) ; }
    }

    
    /*
    printf(" xxx --------------------------------------------- \n");
    printf(" xxx fileName = '%s' \n", fileName);
    printf(" xxx ENVname  = '%s' \n", ENVname );
    printf(" xxx suffix   = '%s' \n", suffix );
    */

    if ( getenv(ENVname) == NULL ) {
      if ( ABORTFLAG ) { SEV = SEV_FATAL; } else { SEV = SEV_WARN; }
      sprintf(c1err,"getenv(%s) failed (callFun=%s)", ENVname, callFun);
      sprintf(c2err,"check fileName='%s'", fileName);
      errmsg(SEV, 0, fnam, c1err, c2err);     
      return(ERROR);
    }

    sprintf(fileName, "%s%s", getenv(ENVname), suffix);
  } 

  if ( NFILE < MXFILE_ENVreplace-1 ) 
    { sprintf(ENVreplace_store.FILENAME_ENVreplace[NFILE],"%s", fileName); }
  ENVreplace_store.NFILE++ ;

  return(SUCCESS) ; 

} // end of ENVreplace

// void envreplace_(char *fileName) { ENVreplace(fileName); }

void ENVrestore(char *fileName_ENVreplace, char *fileName_orig) {

  // Created July 21, 2019
  // return fileName_orig before ENVreplace was called.

  int NFILE = ENVreplace_store.NFILE ;
  int ifile ;
  char *FILENAME_ENVreplace, *FILENAME_ORIG ;
  // ------------- BEGIN -----------

  sprintf(fileName_orig, "%s", fileName_ENVreplace); 

  for(ifile=0; ifile < NFILE; ifile++ ) {
    FILENAME_ENVreplace = ENVreplace_store.FILENAME_ENVreplace[ifile] ;
    FILENAME_ORIG       = ENVreplace_store.FILENAME_ORIG[ifile];
    if ( strcmp(fileName_ENVreplace,FILENAME_ENVreplace) == 0 ) {
      sprintf(fileName_orig,"%s", FILENAME_ORIG );
    }
  }

  return;
} // end ENVrestore


int intrac_() {  return ((int) isatty(0)); }  // needed by fortran minuit (from intrac.c)

void warn_oldInputs(char *varName_old, char *varName_new) {

  int NWARN, i;

  NWARN = OLD_INPUTS.NWARN ;
  if ( strcmp(varName_old,"init") == 0 ) 
    { OLD_INPUTS.NWARN = 0 ; return ; }

  if ( strcmp(varName_old,"list") == 0 ) {
    if ( NWARN == 0 ) { return ; }

    return ; // remove this when new system goes live
    printf("\n");
    for(i=0; i < NWARN; i++ ) {
      printf("  WARNING: REPLACE OLD INPUT  %s  with %s \n"
	     ,OLD_INPUTS.VARNAME_OLD[i]
	     ,OLD_INPUTS.VARNAME_NEW[i] );
      fflush(stdout);
    }
    return ;
  } // end list loop


  sprintf(OLD_INPUTS.VARNAME_OLD[NWARN], "%s", varName_old );
  sprintf(OLD_INPUTS.VARNAME_NEW[NWARN], "%s", varName_new );
  OLD_INPUTS.NWARN++ ;

} // end of warn_oldInputs

void warn_oldinputs__(char *varName_old, char* varName_new) 
{ warn_oldInputs(varName_old, varName_new); }

// ===================================
void set_FILTERSTRING(char *FILTERSTRING) {
  // Feb 2013: return list with every filter defined.
  sprintf(FILTERSTRING,"%s", FILTERSTRING_DEFAULT );
}


void set_EXIT_ERRCODE(int ERRCODE) {   EXIT_ERRCODE = ERRCODE; }
void set_exit_errcode__(int *ERRCODE) { set_EXIT_ERRCODE(*ERRCODE); }

// ====================================================
int IGNOREFILE(char *fileName) {

  // May 2014
  // Return 1 if fileName is "NULL" or "NONE".
  // These are valid names to override an existing fileName 
  // with a command to ignore it.
  // Return 0 otherwise.
  //
  // Dec 18 2017: include ' '

  if ( strcmp_ignoreCase(fileName,"null")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"none")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"blank")  == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"NULL")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"NONE")   == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName,"BLANK")  == 0 ) { return 1 ; }
  if ( strcmp_ignoreCase(fileName," ")      == 0 ) { return 1 ; }
  if ( strlen(fileName) == 0 ) { return 1; } // Sep 2016

  return  0 ;

}  // end of IGNOREFILE

int ignorefile_(char *fileName) { return IGNOREFILE(fileName); }

// ==================================================
int strcmp_ignoreCase(char *str1, char *str2) {

  // Feb 2014: analog of strcmp, but convert to lower case
  //           before testing so that result is case-insensitive.
  //
  // return strcmp on lower-case comparison.

  int len1, len2, j ;
  char str1_lc[100], str2_lc[100];

  // --------- BEGIN ----------
  len1 = strlen(str1) ;
  len2 = strlen(str2) ;
  if ( len1 != len2 ) { return -1 ; }

  for(j=0; j<len1; j++ ) { str1_lc[j] = tolower ( str1[j] ) ;  }
  for(j=0; j<len2; j++ ) { str2_lc[j] = tolower ( str2[j] ) ;  }

  str1_lc[len1] = '\0' ;
  str2_lc[len2] = '\0' ;

  return strcmp(str1_lc,str2_lc) ;
  
} // end of strcmp_ignoreCase

// ==================================================
void invertMatrix(int N, int n, double *Matrix ) {

  // Jan 2013
  // *Matrix is a 2D matrix of dimension NxN
  // but only the n x n subset is to be inverted.
  //
  // input *Matrix is overwritten with its inverse;
  // if you want to save the original matrix, save
  // it before calling this function.

  int s ;
  int i1, i2, J ;
  
  // -------------- BEGIN ---------------
  // Define all the used matrices
  gsl_matrix * m         = gsl_matrix_alloc (n, n);
  gsl_matrix * inverse   = gsl_matrix_alloc (n, n);
  gsl_permutation * perm = gsl_permutation_alloc (n);

  // Fill the matrix m

  J = 0;
  for ( i1=0; i1 < N; i1++ ) {
    for ( i2=0; i2 < N; i2++ ) {
      
      if ( i1 < n && i2 < n )  { 
	gsl_matrix_set(m,i1,i2, Matrix[J] );
      }

      J++ ;
    }
  }

  // Make LU decomposition of matrix m
  gsl_linalg_LU_decomp (m, perm, &s);

  // Invert the matrix m
  gsl_linalg_LU_invert (m, perm, inverse);

  // load inverse into Matrix
  J = 0;
  for ( i1=0; i1 < N; i1++ ) {
    for ( i2=0; i2 < N; i2++ ) {
      
      if ( i1 < n && i2 < n )  { 
	Matrix[J] = gsl_matrix_get(inverse,i1,i2) ;
      }

      J++ ;
    }
  }

  gsl_matrix_free(m) ;
  gsl_matrix_free(inverse) ;
  gsl_permutation_free(perm) ;

}  // end of invertMatrix

void invertmatrix_(int *N, int *n, double *Matrix ) {
  invertMatrix(*N, *n, Matrix ) ;
}


void sortDouble(int NSORT, double *ARRAY, int ORDER, 
		int *INDEX_SORT) {

  // Created Jan 12 2013 by R.Kessler
  // Wrapper to sort NSORT elements of ARRARY (replaces CERNLIB's SORTZV)
  // ORDER = -1/+1  -> decreasing/increasing order
  // Output is *INDEX_SORT (input ARRAY is NOT changed)

  size_t stride_t = 1;
  size_t n_t, *index_t ;
  int i;

  // ---------------- BEGIN ----------------

  n_t = NSORT ;
  index_t = (size_t*)malloc(n_t * sizeof(size_t) );

  // native function is double
  gsl_sort_index(index_t, ARRAY, stride_t, n_t);

  // load output arg.
  for(i=0; i < NSORT; i++ ) { INDEX_SORT[i] = index_t[i]; }

  // check to flip from increasing to decreasing order.
  if ( ORDER < 0 ) 
    { reverse_INDEX_SORT(NSORT, INDEX_SORT); }
  
  free(index_t);

} // end of sortDouble

// mangled sort function for fortran
void sortDouble_(int *NSORT, double *ARRAY, int *ORDER, 
		 int *INDEX_SORT) {
  sortDouble(*NSORT, ARRAY, *ORDER, INDEX_SORT) ;

  int i;
  // add 1 to INDEX_SORT so that index starts at 1 instead of 0
  for ( i=0; i < *NSORT; i++ ) { INDEX_SORT[i] += 1; }
}


void sortFloat(int NSORT, float *ARRAY, int ORDER, 
	       int *INDEX_SORT) {

  // Created Jan 12 2013 by R.Kessler
  // Wrapper to sort NSORT elements of ARRARY (replaces CERNLIB's SORTZV)
  // ORDER = -1/+1  -> decreasing/increasing order
  // Output is *INDEX_SORT (input ARRAY is NOT changed)

  size_t stride_t = 1;
  size_t n_t, *index_t ;
  int i;

  // ---------------- BEGIN ----------------

  n_t = NSORT ;
  index_t = (size_t*)malloc(n_t * sizeof(size_t) );

  gsl_sort_float_index(index_t, ARRAY, stride_t, n_t );

  for(i=0; i < NSORT; i++ ) { INDEX_SORT[i] = index_t[i]; }
  
  // check to flip from increasing to decreasing order.
  if ( ORDER < 0 )  { 
    reverse_INDEX_SORT(NSORT, INDEX_SORT); 
  }
  
  free(index_t);

} // end of sortFloat

// mangled sort function for fortran
void sortfloat_(int *NSORT, float *ARRAY, int *ORDER, 
		int *INDEX_SORT) {
  sortFloat(*NSORT, ARRAY, *ORDER, INDEX_SORT) ;

  int i;
  // add 1 to INDEX_SORT so that index starts at 1 instead of 0
  for ( i=0; i < *NSORT; i++ ) { INDEX_SORT[i] += 1; }
}



void sortInt(int NSORT, int *ARRAY, int ORDER, int *INDEX_SORT) {

  // Created Jan 12 2013 by R.Kessler
  // Wrapper to sort NSORT elements of ARRARY (replaces CERNLIB's SORTZV)
  // ORDER = -1/+1  -> decreasing/increasing order
  // Output is *INDEX_SORT (input ARRAY is NOT changed)

  size_t stride_t = 1;
  size_t n_t, *index_t ;
  int i;

  // ---------------- BEGIN ----------------

  n_t = NSORT ;
  index_t = (size_t*)malloc(n_t * sizeof(size_t) );

  gsl_sort_int_index(index_t, ARRAY, stride_t, n_t);
  
  for(i=0; i < NSORT; i++ ) { INDEX_SORT[i] = index_t[i]; }

  // check to flip from increasing to decreasing order.
  if ( ORDER < 0 ) 
    { reverse_INDEX_SORT(NSORT, INDEX_SORT); }
  

  free(index_t);

} // end of sortInt

// mangled sort function for fortran
void sortint_(int *NSORT, int *ARRAY, int *ORDER, 
	      int *INDEX_SORT) {
  sortInt(*NSORT, ARRAY, *ORDER, INDEX_SORT) ;

  int i;
  // add 1 to INDEX_SORT so that index starts at 1 instead of 0
  for ( i=0; i < *NSORT; i++ ) { INDEX_SORT[i] += 1; }
}


void reverse_INDEX_SORT(int NSORT, int *INDEX_SORT) {

  // revserse the order of INDEX_SORT.
  // Note that input INDEX_SORT array is changed.

  int *INDEX_TMP, i ;


  INDEX_TMP = (int*)malloc( NSORT * sizeof(int) );
    
  for ( i=0; i < NSORT; i++ ) 
    {  INDEX_TMP[i] = INDEX_SORT[i];  }

  for ( i=0; i < NSORT; i++ )  
    {  INDEX_SORT[i] = INDEX_TMP[NSORT-i-1]; }


  free(INDEX_TMP);
    
}  // end of revserse_INDEX_SORT

// ======================================================
void print_KEYwarning(int ISEV, char *key_old, char *key_new) {

  // print warning message to use the new key instead of the old key
  // ISEV = SEV_WARN  -> give warning but do not abort.
  // ISEV = SEV_FATAL -> abort

  char fnam[] = "print_KEYwarning";
  char line[] = "%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%" ;

  if ( ISEV == SEV_FATAL ) {
    sprintf(c1err,"key = '%s' is no longer valid.", key_old);
    sprintf(c2err,"Must use new key = '%s'", key_new);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);     
  }
  else {
    printf("\n");
    printf("%s%s\n", line, line);
    printf("\n");
    printf(" WARNING: use new '%s' key instead of old '%s' key' \n", 
	   key_new, key_old);
    printf("          Old '%s' key will soon become obsolete.\n", 
	   key_old);
    printf("\n");
    printf("%s%s\n", line, line);
    printf("\n");
  }

} // end of print_KEYwarning


// ***************************************************
double SNR_calculator(double ZPT, double PSF, double SKYMAG, double MAG,
		      double *FLUX_and_ERR ) {

  // Created May 2014 by R. Kessler
  // Return signal-to-noise (SNR) for input MAG and 
  // input observing conditions,
  //   ZPT:    Npe = 10**[-0.4*(mag-ZPT)]
  //   PSF:    FWHM, arcsec
  //   SKYMAG: mag/arcsec^2
  //
  // Jun 1 2014: return FLUX and its error in FLUX_and_ERR array.

  double SQSNR, SNR, tA, F, FSKY, arg ;
  double OMEGA = 1.51; // Pararam to define effective sky-noise area
                       //  A = ( OMEGA * PSF )^2
  // ------------- BEGIN -----------

  tA    = ( OMEGA * PSF ) * ( OMEGA * PSF );  // area
  arg   = -0.4*(MAG    - ZPT);   F    = pow(10.0,arg);
  arg   = -0.4*(SKYMAG - ZPT);   FSKY = pow(10.0,arg);
  SQSNR = (F*F)/( F + tA*FSKY );
  SNR   = sqrt(SQSNR);

  // load flux and its error to allow external code to get coadded SNR
  FLUX_and_ERR[0] = F ;
  FLUX_and_ERR[1] = sqrt( F + tA*FSKY );

  return SNR ;

} // end of SNRMAG

// ************************************************************
double MAGLIMIT_calculator(double ZPT, double PSF, double SKYMAG, double SNR){

  // October 2010 by P. Belov & S.Glazov
  //
  // compute and return point-source limiting mag corresponding to SNR = S/N
  // ZPT:    Npe = 10**[-0.4*(mag-ZPT)]
  // PSF:    FWHM, arcsec
  // SKYMAG: mag/arcsec^2
  // SNR   : S/N for limiting mag
  //
  // Aug 7, 2013 RK - fix bugs computing MLIM
  // May 17 2014 RK - moved from simulation code

  double MLIM;
  double tMLIM;  // Temp value for MLIM
  double EA    = 1.00; // The fraction of source-flux
  double OMEGA = 1.51; // Pararam to define effective sky-noise area
                       //  A = ( OMEGA * PSF )^2
  double eps; // epsilon
  double acc; // accuracy
  double arg;

  int cnt, cnt_orig; // counter in order to avoid infinite loop
  double tA, t1, t2; // Temp variables

  char fnam[] = "MAGLIMIT" ;

  // -------------- BEGIN --------------

  // check for crazy args
  if ( ZPT <= 1.01 || PSF <= 0.001 || SKYMAG <= 1.0 ) {
    MLIM = 0.0 ;
    return MLIM;
  }

  eps = 10000.0  ;
  acc =     0.001;  // required accuracy for convergence
  cnt = cnt_orig =  100 ;  
  t2  =  5.0 * log10(SNR/EA);
  tA  = ( OMEGA * PSF ) * ( OMEGA * PSF );  // area
  t1  = 2.5 * log10( tA );
  MLIM = 20.0 ; // very rough guess.

  // here we get the value for mlimit by iterative procedure
  while ( ( eps > acc ) && ( cnt > 0 ) ) {

    if ( cnt == cnt_orig ) 
      // initial guess is with no signal-noise
      {  tMLIM  = 0.5 * ( ZPT + SKYMAG - t1 - t2 ); }
    else
      {  tMLIM = MLIM ; } // previous iteration

    arg    = 0.4*(SKYMAG-tMLIM) ;  // RK added this Feb 2011


    MLIM = 0.5*(tMLIM + ZPT - t2) - 1.25*log10(1.0 +  tA*pow(10.0,-arg) );
    eps    = fabs( MLIM - tMLIM );

    /*
    printf(" xxx cnt=%3d: MLIM=%f  tMLIM=%f  arg=%f \n",
    cnt, MLIM, tMLIM, arg); */
    cnt   -= 1 ;
  }

  if ( cnt == 0 ) {    
    sprintf(c1err,"MLIM=%.3f value has not converged: eps=%f", 
	    MLIM, eps );
    sprintf(c2err,"ZP=%.3f  PSF=%.3f  SKYMAG=%.2f SNR=%.1f ",
	    ZPT, PSF, SKYMAG, SNR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  /*
  double SNR_CHECK = SNR_calculator(ZPT,PSF,SKYMAG, MLIM);
  printf(" xxxx SNR_CHECK/SNR = %.4f/%.4f = %f \n",
	SNR_CHECK, SNR, SNR_CHECK/SNR ); fflush(stdout);
  */


  return MLIM;

} // end of SIMLIB_maglimit


// ************************************************************
double NoiseEquivAperture(double PSFSIG1, double PSFSIG2, 
			  double PSFratio) {

/* ===========================================================
   Returns effective sky aperture, or noise-equivalent area
   
   = 1/ integral[PSF^2 rdr dtheta]
   
   for double-gaussian PSF where

   SIG1  : sigma of inner gaussian       (psf_sigma1 in psField)
   SIG2  : sigma of outer gaussian       (psf_sigma2 in psField)
   ratio : ratio of gaussians at origin  (pdf_b      in psField)

   PSF(r) is defined locally as
   
              exp(-r^2/2\sig1^2)            exp(-r^2/2\sig2^2)
   PSF = A1 * -------------------   +  A2 * ------------------
               2 * PI * sig1^2               2 * PI * sig1^2 
 
  where A1 + A2 = 1.00 for normalized PSF.
  With this definition,
 
            A2/sig2^2                    1
   ratio = -----------  ,  A1 = ------------------------
            A1/sig1^2            1 + ratio*sig2^2/sig1^2
 
 
  and effective aperture  = 
 
            4 * PI * ( sig1^2 + sig2^2 )
    =  ------------------------------------------
         1 + [ A1 * sig2/sig1 + A2*sig1/sig2]^2
 
  Special cases/checks:
   * sig2 = A2 = 0 : Aperture -> 4 * PI * sig1^2
   * sig1 = sig2   : Aperture -> 4 * PI * sig1^2 
   * sig2 = 2*sig1 : Aperture -> 4 * PI * sig1^2 / ( 1 - A2*6/5 )


  Note:  moved here from snana.car (SKY_APERTURE) on Feb 27, 2012,
         and REAL*4 -> double.
 
   ===================== */

  // local var


  //  double PI = 3.14159265 ;
  double PI = TWOPI/2.0 ;
  double A, A1, A2, SQSIG1, SQSIG2, SQSUM, TMP ;

  // ---------- BEGIN -----------

  A = 0.0 ;

  SQSIG1 = PSFSIG1 * PSFSIG1 ;
  SQSIG2 = PSFSIG2 * PSFSIG2 ;
  SQSUM  = SQSIG1 + SQSIG2 ;

  // check for single-gaussian PSF

  if ( PSFratio < 1.0E-5 || PSFSIG2 < 0.0001 ) {
    A = 4.0 * PI * SQSIG1 ;
    return A ;
  }

  TMP  = PSFratio * SQSIG2 / SQSIG1  ;
  A1   = 1.0 / ( 1.0 + TMP ) ;
  A2   = 1.0 - A1 ;
  TMP  = A1 * PSFSIG2/PSFSIG1  +  A2 * PSFSIG1/PSFSIG2 ;
  A    = 4.0 * PI * SQSUM / ( 1.0 + TMP*TMP ) ;
  return A ;

}  // end of   NoiseEquivAperture 

double noiseequivaperture_(double *PSFSIG1, double *PSFSIG2, 
			   double *PSFratio){
  return NoiseEquivAperture(*PSFSIG1, *PSFSIG2, *PSFratio) ;
}

// =================================================
double PROB_Chi2Ndof(double chi2, int Ndof ) {

  // Nov 17, 2011, R.Biswas
  // Return probability that chi2 > chi2 with Ndof
  // Note that this function should replaces CERNLIB's PROB function.
  //
  // For ideal Gaussian noise and correct uncertainties, FITPROB
  // distribution should be a flat between 0 and 1.
  // FITPROB pile-up near 1 suggests that errors are too large;
  // FITPROB pile-up near 0 suggests contamination, or errors too small.
  //
  // Feb 3, 2012: protect againsta Ndof < 1 -> returns P=0
  //
  
  double N, chi2red, P;
  N = (double)Ndof/2.0 ;
  chi2red = chi2/2.0 ;

  if ( Ndof < 1 ) 
    { P = 0.0 ; }   
  else if ( chi2 < 0.0 ) 
    { P = 1.0 ; }
  else
    { P = gsl_sf_gamma_inc_Q(N, chi2red); }

  return P ;

} // end of PROB_Chi2Ndof

double prob_chi2ndof__(double *chi2, int *Ndof) {
  return PROB_Chi2Ndof(*chi2, *Ndof);
} 

// ====================================
int getInfo_PHOTOMETRY_VERSION(char *VERSION      // (I) photometry version 
			       ,char *DATADIR     // (I/O) dir with data files 
			       ,char *LISTFILE    // (O) name of list file
			       ,char *READMEFILE  // (O) name of readme file
			       ) {

  // Created Mar 4, 2011 by R.Kessler
  // For input photometry VERSION, function returns 
  //   DATADIR:     SNANA data directory 
  //   LISTFILE:    file with list of data files in DATADIR
  //   READMEFILE:  file describing data sample
  //
  // If DATADIR is non-null, then it is passed from PRIVATE_DATA_PATH,
  // so use this directory to  determine LISTFILE and READMEFILE.
  // Do NOT change DATADIR in this case. If DATADIR=='', then set
  // to default $SNDATA_ROOT/lcmerge or path in PATH_SNDATA_SIM.DAT
  //
  // Function return argument is SUCCESS or ERROR
  // 
  // Dec 2 2012: rename function, 
  //   getInfo_SNANA_VERSION -> getInfo_PHOTOMETRY_VERSION
  //   to avoid confusion with snana version (i.e, v10_19)
  //
  // Nov 11 2014: allow data to be in a folder under lcmerge/[VERSION]
  // Feb 10, 2015: same fix for DATADIR and DATADIR/[VERSION]
  // Feb 26, 2015: for datadir, also check SNDATA_ROOT/SIM
  //
  // Nov 18 2017: 
  //  + check for user-define sim-paths to be compatible with
  //    sim-input option PATH_SNDATA_SIM.. See PATHLIST below.
  //    
  // Sep 19 2018: 
  //  remove local & obsolete MXDIR_CHECK and instead use global
  //  MXPATH_SNDATA_SIM. Abort if NDIR_CHECK >= MXPATH_SNDATA_SIM .
  //
  // Sep 12 2019: 
  //  + abort if DATADIR corresponds to $SNDATA_ROOT/SIM or any SIM path
  // 
  // Apr 20 2021
  //  To auto-find sub-folders under $SNDATA_ROOT/lmerge,
  //  if VERSION = PREFIX_SUFFIX, also check lcmerge/PREFIX/PREFIX_SUFFIX
  //

  char 
    SNDATA_ENV[20] = "SNDATA_ROOT"
    ,SNDATA_ROOT[MXPATHLEN]
    ,tmpDir[MXPATH_SNDATA_SIM][MXPATHLEN]
    ,tmpFile[MXPATH_SNDATA_SIM][MXPATHLEN]
    ,prefix[MXPATHLEN]
    ,fnam[] = "getInfo_PHOTOMETRY_VERSION"
    ;

  int idir, ifound, NFOUND, NDIR_CHECK, j_ ;
  int idirFOUND[MXPATH_SNDATA_SIM];
  int LDMP = 0;
  FILE *fp ;

  // ---------- BEGIN -----------

  /*
  printf(" xxx %s: check VER='%s' in PATH='%s' \n",
	 fnam, VERSION, DATADIR );
  */
  // init outputs to NULLSTRING value

  sprintf(LISTFILE,   "%s", NULLSTRING );
  sprintf(READMEFILE, "%s", NULLSTRING );

  // first make sure that required env is set.
  if ( getenv(SNDATA_ENV) == NULL ) {
    sprintf(c1err,"env variable '$%s' is not set.",  SNDATA_ENV );
    sprintf(c2err,"%s", "      ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  sprintf(SNDATA_ROOT, "%s", getenv(SNDATA_ENV) ) ;


  // define list of directories to check for data
  idir=0;

  if ( !IGNOREFILE(DATADIR) ) { 
    // private user dir
    sprintf(tmpDir[idir], "%s" ,          DATADIR ); 
    sprintf(tmpFile[idir],"%s/%s.LIST",   tmpDir[idir], VERSION  );
    idir++ ;

    sprintf(tmpDir[idir], "%s/%s",        DATADIR, VERSION);
    sprintf(tmpFile[idir],"%s/%s.LIST",   tmpDir[idir], VERSION  );
    idir++ ;
    
  }
  else {
    // default locations under SNDATA_ROOT
    sprintf(tmpDir[idir], "%s/lcmerge/%s",  SNDATA_ROOT, VERSION);
    sprintf(tmpFile[idir],"%s/%s.LIST",     tmpDir[idir], VERSION  );
    idir++ ;

    j_ = index_charString("_", VERSION);
    if ( j_ > 1 ) {
      strncpy(prefix,VERSION,j_);   prefix[j_] = '\0' ;

      /* 
      printf(" xxx j_ = %d  prefix = '%s' for VERSION = %s \n",
      j_, prefix, VERSION); */

      sprintf(tmpDir[idir], "%s/lcmerge/%s/%s",  
	      SNDATA_ROOT, prefix, VERSION);
      sprintf(tmpFile[idir],"%s/%s.LIST",  tmpDir[idir], VERSION  );
      idir++ ;
    }
    
  }

  // - - - - - - - - - - - - - - - - - - - - - 
  // Nov 18 2017: check for user-defined SIM-output dirs 
  //              based on PATH_SNDATA_SIM
  int ipath, NPATH;
  char **PATHLIST;

  PATHLIST = (char**) malloc ( MXPATH_SNDATA_SIM * sizeof(char*) );
  for(ipath=0; ipath < MXPATH_SNDATA_SIM; ipath++ ) 
    { PATHLIST[ipath] = (char*) malloc ( MXPATHLEN * sizeof(char) ); }
  NPATH = getList_PATH_SNDATA_SIM(PATHLIST);

  // tack on default SIM dir (Sep 2019)
  int IPATH_SIM_DEFAULT = NPATH;
  sprintf(PATHLIST[NPATH], "%s/SIM", SNDATA_ROOT ); NPATH++ ;

  if ( LDMP ) 
    { printf(" xxx DATADIR = '%s' \n", DATADIR); fflush(stdout); }

  for(ipath = 0 ; ipath < NPATH; ipath++ ) {

    if ( idir < MXPATH_SNDATA_SIM ) {
      sprintf(tmpDir[idir],  "%s/%s" , PATHLIST[ipath],  VERSION );
      sprintf(tmpFile[idir], "%s/%s.LIST", tmpDir[idir], VERSION  );

      if ( LDMP) {
	printf(" xxx check PATHLIST[%d] = '%s' \n", ipath, PATHLIST[ipath] );
	fflush(stdout); 
      }

      // Sep 12 2019: abort if DATADIR corresponds to any SIM path
      if ( strcmp(DATADIR,PATHLIST[ipath])== 0 ) {
	print_preAbort_banner(fnam);    
	printf("\t PRIVATE_DATA_PATH = '%s' \n", DATADIR);

	if ( ipath == IPATH_SIM_DEFAULT ) {
	  sprintf(c1err,"PRIVATE_DATA_PATH cannot be the same as");
	  sprintf(c2err,"$SNDATA_ROOT/SIM");
	}
	else {
	  sprintf(c1err,"PRIVATE_DATA_PATH cannot match any path in");
	  sprintf(c2err,"$SNDATA_ROOT/SIM/%s", PATH_SNDATA_SIM_LIST );
	}
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

    }
    idir++ ;
  }

  for(ipath=0; ipath < MXPATH_SNDATA_SIM; ipath++ ) 
    { free(PATHLIST[ipath]) ; }
  free(PATHLIST);

  // -------
  NDIR_CHECK = idir;

  if ( NDIR_CHECK >= MXPATH_SNDATA_SIM ) {
    sprintf(c1err,"NDIR_CHECK=%d exceeds bound", NDIR_CHECK);
    sprintf(c2err,"of MXPATH_SNDATA_SIM=%d", MXPATH_SNDATA_SIM);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // set DATADIR to directory with LIST file
  NFOUND = 0 ;

  for(idir=0; idir < NDIR_CHECK; idir++ ) {
    if ( (fp = fopen(tmpFile[idir], "rt")) != NULL )  { 
      sprintf(DATADIR, "%s", tmpDir[idir]); 
      idirFOUND[NFOUND] = idir;
      fclose(fp);
      NFOUND++ ;
    }
  } // end of idir


  if ( NFOUND == 0 ) {
    print_preAbort_banner(fnam);    
    for(idir=0; idir < NDIR_CHECK; idir++ ) {
      printf("   data not in '%s' \n", tmpDir[idir] );
    }	   
    sprintf(c1err,"Could not find SNANA version '%s'",  VERSION );
    sprintf(c2err,"Check directories above");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( NFOUND > 1 ) {
    print_preAbort_banner(fnam);    
    for(ifound=0; ifound < NFOUND; ifound++ ) {
      idir = idirFOUND[ifound];
      printf("   Found %s \n", tmpDir[idir] );
    }
    sprintf(c1err,"Found %d VERSION='%s'", NFOUND, VERSION);
    sprintf(c2err,"Only one unique version allowed.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


  // always use DATADIR to construct name of list-file and readme file.
  sprintf(LISTFILE,   "%s/%s.LIST",   DATADIR, VERSION );
  sprintf(READMEFILE, "%s/%s.README", DATADIR, VERSION );

  return SUCCESS ;

}  // end of function


int getinfo_photometry_version__(char *VERSION, char *DATADIR, 
				 char *LISTFILE, char *READMEFILE) {
  return getInfo_PHOTOMETRY_VERSION(VERSION,DATADIR,LISTFILE,READMEFILE);
}


// ==========================================
FILE *openFile_PATH_SNDATA_SIM(char *mode) {

  // Open file for 
  //  + reading (mode='read')
  //  + append  (mode='append')
  //
  // modeArg[2] -> modeArg[4] (fix Mac issue)
  //
  // Feb 2021: abort if open fails (e.g. file is write-protected)

  char fileName[MXPATHLEN], SNDATA_ROOT[MXPATHLEN] ;
  char modeArg[4];
  FILE *fp ;
  char fnam[] = "openFile_PATH_SNDATA_SIM" ;

  // ------------- BEGIN --------------

  sprintf(SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") ) ;

  // hard-wire name of file with list of alternate PATH_SNDATA_SIM
  sprintf(fileName, "%s/SIM/%s", SNDATA_ROOT, PATH_SNDATA_SIM_LIST );
  sprintf(modeArg, "%ct", mode[0] );
  fp = fopen(fileName,modeArg);

  if ( !fp ) {
    sprintf(c1err,"Cannot open PATH_SNDATA_SIM file in %s mode:", mode);
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  //  printf("\n Open %s in %s-mode (%s)\n", fileName, mode, modeArg);
  return(fp) ;

} // end openFile_PATH_SNDATA_SIM


void add_PATH_SNDATA_SIM(char *PATH) {

  // Created Nov 18 2017
  // Add PATH to file in $SNDATA_ROOT/SIM so that analysis codes 
  // know where else to check for simulated data files.  If PATH 
  // name is already written, don't write anything.

  FILE *fp;
  int  EXIST=0;
  int  NPATH_DEJA=2; // accounts for default /lcmerge and /SIM 
  char PATH_ENVreplace[MXPATHLEN];
  char ftmp[MXPATHLEN];
  char fnam[] = "add_PATH_SNDATA_SIM" ;

  // -------------- BEGIN ---------------

  // first read file to see if PATH already exists
  fp = openFile_PATH_SNDATA_SIM("read");
  sprintf(PATH_ENVreplace,"%s", PATH);  
  ENVreplace(PATH_ENVreplace,fnam,1);

  if ( fp ) {
    while( (fscanf(fp, "%s", ftmp)) != EOF) {
      ENVreplace(ftmp,fnam,1);
      NPATH_DEJA++ ;
      if ( strcmp(PATH_ENVreplace,ftmp) == 0  ) { EXIST=1; }
    }
    fclose(fp);
  }

  // if PATH does not already exist, append it to list file
  if ( EXIST == 0 ) {
    if ( NPATH_DEJA >= MXPATH_SNDATA_SIM ) {
      print_preAbort_banner(fnam);    
      printf("   NPATH_DEJA = %d (includes /lcmerge and /SIM) \n", NPATH_DEJA);
      printf("   MXPATH_SNDATA_SIM = %d \n", MXPATH_SNDATA_SIM);

      sprintf(c1err,"Too many paths defined in ");
      sprintf(c2err,"$SNDATA_ROOT/SIM/%s", PATH_SNDATA_SIM_LIST);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }
    fp = openFile_PATH_SNDATA_SIM("append");
    fprintf(fp, "%s\n", PATH);
    fclose(fp);
  }

  return;

} // end add_PATH_SNDATA_SIM


int  getList_PATH_SNDATA_SIM(char **pathList) {

  // Created Nov 2017
  // Return list (pathList) of user-defined PATH_SNDATA_SIM .
  // Function returns NPATH, number on list.
  // This list is used by analysis codes to look for
  // simulated data files so that user doesn't have
  // to worry about writing sim files to alternate
  // directories
  //
  // Mar 30 2019: refactor so that undefined paths are skipped
  //              without aborting.
  //
  int ENVstat, NPATH=0;
  FILE *fp ;
  char path[MXPATHLEN] ;
  char fnam[] = "getList_PATH_SNDATA_SIM" ;

  // ------------ BEGIN ---------------- 

  // open list file in read mode
  fp = openFile_PATH_SNDATA_SIM("read");
  if ( !fp ) { return(0); }   // return if it does not exist.
 
  // scoop up each PATH in list file.
  while( (fscanf(fp, "%s", path)) != EOF) {

    // do NOT abort on invalid ENV
    ENVstat = ENVreplace(path,fnam,0); 
    if ( ENVstat != SUCCESS ) { continue ; }

    if (NPATH<MXPATH_SNDATA_SIM) { sprintf(pathList[NPATH], "%s", path); }
    NPATH++ ;
  }
  fclose(fp);

  // abort if too many paths.
  if ( NPATH >= MXPATH_SNDATA_SIM ) {
    sprintf(c1err,"NPATH=%d exceeds bound of MXPATH_SNDATA_SIM=%d",
	    NPATH, MXPATH_SNDATA_SIM );
    sprintf(c2err,"Check opened file above." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  return(NPATH) ;

} // end getList_PATH_SNDATA_SIM


// =============================================================
void arrayStat(int N, double *array, double *AVG, double *STD, double *MEDIAN) {

  // For input *array return *AVG and *RMS
  // Jun 2 2020: include MEDIAN in output
  // Jul 21 2021: rename RMS -> STD to avoid confusion.
  // Sep 30 2021: fix median calc based on odd or even N

  int i;
  double XN, val, avg, sum, sqsum, std, median ;

  // ----------- BEGIN ------------

  *AVG = *STD = *MEDIAN = 0.0 ;
  if ( N <= 0 ) { return ; }

  avg  = std = sqsum = sum = median = 0.0 ;
  XN   = (double)N ;

  for ( i=0; i < N; i++ ) 
    { val=array[i] ; sum += val; sqsum += (val*val); }

  avg = sum/XN ; 
  std = STD_from_SUMS(N, sum, sqsum);

  // for median, sort list and them median is middle element
  int ORDER_SORT  = +1;
  int *INDEX_SORT = (int*) malloc( N * sizeof(int) ) ;
  int imed0, imed1, iHalf       = N/2;
  sortDouble(N, array, ORDER_SORT, INDEX_SORT );

  // xxx mark delete 9.30 2021: imedian = INDEX_SORT[iHalf];
  // xxx median  = array[imedian];

  // xxx  printf(" xxx N=%d ihalf=%d \n", N, iHalf);

  if ( N%2 == 1 ) { 
    // odd number of elements -> use middle value
    imed0   = INDEX_SORT[iHalf];
    median  = array[imed0]; 
  }
  else {
    // even number of elements; average middle two values
    imed0 = INDEX_SORT[iHalf-1];
    imed1 = INDEX_SORT[iHalf];
    median  = 0.5 * ( array[imed0] + array[imed1] );
  }
  

  // load output array.
  *AVG    = avg ;
  *STD    = std ;
  *MEDIAN = median;

  return ;
} // end of arrayStat

void arraystat_(int *N, double *array, double *AVG, double *RMS, 
		double *MEDIAN) 
{ arrayStat(*N, array, AVG, RMS, MEDIAN); }

void test_arrayStat(void) {
#define NTEST_arrayStat 9
  double AVG, RMS, MEDIAN, array[NTEST_arrayStat];
  int N, i;
  char fnam[] = "test_arrayStat" ;

  // ----------- begin -----------
  // load array values 1 ... N
  for(i=0; i < NTEST_arrayStat; i++ ) { array[i] = (double)(i+1); }

  for(N=NTEST_arrayStat; N > NTEST_arrayStat-2; N-- ) {
    arrayStat(N, array, &AVG, &RMS, &MEDIAN);
    printf(" xxx %s: N=%d -> AVG=%.3f, RMS=%.3f, MEDIAN=%.3f \n",
	   fnam, N, AVG, RMS, MEDIAN ); fflush(stdout);
  }
  return;

} // end test_arrayStat

// ========================================================
double STD_from_SUMS(int N, double SUM, double SQSUM) {

  // Created Aug 2017
  // Compute standard deviation (STD) from sums.

  double STD = 0.0 ;
  double XN  = (double)N;
  if ( N == 0 ) { return(STD); }

  double ARG = SQSUM/XN - pow((SUM/XN),2.0) ;
  if ( ARG > 0.0 ) { STD = sqrt(ARG); }

  return(STD);

} // end STD_from_SUMS

// ======================================================
double sigint_muresid_list(int N, double *MURES_LIST, double *MUCOV_LIST,
			   int OPTMASK, char *callFun ) { 


  // Created July 24 2021 by R.Kessler [extracted from SALT2mu code]
  // Unility to compute sigint from N mu-residuals
  //   MURES_LIST   : mu - mu(true,fit)
  //   MUCOV_LIST   : covariance per mu
  //
  // Strategy is to make first pass over LIST and compute sigint_approx. 
  // Then make another pass over LIST and compute RMS_PULL on a grid of 
  // sigint in small steps around sigint_approx. Finally, interpolate 
  // sigint vs. RMS_PULL at RMS_PULL=1.0
  //
  // OPTMASK & 1 --> do not abort on negative sigint
  //    will return negative sqrt(abs(arg)) as flag/quantitative info
  //
  // OPTMASK & 32 --> implement test feature
  // OPTMASK & 64 --> print debug dump
  //
  // *callFun is for error messages.
  //
  // Sep 27 2021
  //   + implement debug dump with OPTMASK & 64
  //   + reduce covtotfloor from 0.1^2 to 0.02^2
  //

  bool LABORT = (OPTMASK & 1) == 0 ;
  bool LTEST  = (OPTMASK & 32) > 0 ;
  bool LDMP   = (OPTMASK & 64) > 0 ;
  
  int    OPT_INTERP  = 1 ;
  double sigint_bin  = 0.01 ;
  double sigint_min = -0.3 ;
  double covtotfloor = 0.02*0.02 ; // protection for negative covtot
  if ( LTEST ) { covtotfloor = 0.1*0.1; } // see -927 flag in SALT2mu

  int    nbin_lo     = 30 ; // prep this many bins below sigint_approx
  int    nbin_hi     = 30 ; // idem above sigint_approx
  double XN          = (double)N;

#define MXSTORE_PULL 100

  int i ;
  double STD_MURES_ORIG, SQMURES, MURES, MUCOV ;
  double SUM_MUCOV = 0.0, SUM_MURES = 0.0, SUM_SQMURES=0.0 ;
  double sigint = 0.0, sigint_approx, tmp;
  double AVG_MUCOV, AVG_MUERR, AVG_MURES ;
  char fnam[] = "sigint_muresid_list";
  
  // ---------------- BEGIN -------------
 

  for ( i=0; i < N ; i ++ ) {
    MURES    = MURES_LIST[i];
    MUCOV    = MUCOV_LIST[i];
    SQMURES  = MURES * MURES ;
    if ( MUCOV <= 0.0 ) {
      print_preAbort_banner(fnam);
      printf("  %s called from %s\n", fnam, callFun);
      sprintf(c1err,"Invalid MUCOV = %le ", MUCOV);
      sprintf(c2err,"i=%d  MURES=%le", i, MURES);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
    }

    SUM_MUCOV   += MUCOV ;
    SUM_MURES   += MURES ;
    SUM_SQMURES += SQMURES ;
  }

  AVG_MURES = SUM_MURES / XN ;
  AVG_MUCOV = SUM_MUCOV / XN;
  AVG_MUERR = sqrt(AVG_MUCOV);
  STD_MURES_ORIG = STD_from_SUMS(N, SUM_MURES, SUM_SQMURES);

  tmp = STD_MURES_ORIG*STD_MURES_ORIG - AVG_MUCOV ;

  
  bool INVALID_SIGINT_APPROX = (STD_MURES_ORIG < AVG_MUERR);
  if  ( INVALID_SIGINT_APPROX && LABORT ) {
      print_preAbort_banner(fnam);
      printf("  %s called from %s\n", fnam, callFun);
      sprintf(c1err,"Cannot compute sigint because RMS < AVG_MUERR ??");
      sprintf(c2err,"RMS=%le, sqrt(AVG_COV)=%le  N=%d",
	      STD_MURES_ORIG, sqrt(AVG_MUCOV), N );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;       
  } // end INVALID_SIGINT && ABORT
  

  if (tmp<0) {
    sigint_approx = 0.;
  } else {
    sigint_approx = sqrt(tmp);
  }

  // - - - - - - - 
  // prepare interp-grid of sigint vs. RMS around sigint_approx.
  
  int    NBIN_SIGINT = 0 ;
  double sigTmp_lo   = sigint_approx - (nbin_lo*sigint_bin) - 1.0E-7 ;
  double sigTmp_hi   = sigint_approx + (nbin_hi*sigint_bin) ;
  double sigTmp, covTmp, covtot, sum_dif, sum_sqdif;
  double pull, sum_pull, sum_sqpull ;
  double sigTmp_store[MXSTORE_PULL], rmsPull_store[MXSTORE_PULL], rmsPull ;
  double ONE = 1.0 ;

  bool BOUND_ONE = false;
  
  if ( LDMP ) {
    printf(" xxx - - - - - - - - - - - - \n");
    printf(" xxx %s: debug dump for %s\n", fnam, callFun);
    printf(" xxx %s: AVG[MURES,MUCOV,MUERR] = %.3f, %.5f, %.3f \n",
	   fnam, AVG_MURES, AVG_MUCOV, AVG_MUERR);
    printf(" xxx %s: STD(MURES_ORIG) = %.3f \n", 
	   fnam, STD_MURES_ORIG);
    printf(" xxx %s: sigint[approx, (lo-hi)] = %.3f,  (%.3f to %.3f) \n",
	   fnam, sigint_approx, sigTmp_lo, sigTmp_hi);
    fflush(stdout);
  }

  // start with largest sigInt and decrease so that RMS is increasing     
  // for the interp function below                                        
  //for(sigTmp = sigTmp_hi; sigTmp >= sigTmp_lo; sigTmp -= sigint_bin ) {
  //printf("xxx RMS_MURES_ORIG=%f sqrt(AVG_MUCOV)=%f sigTmp_hi=%f\n",
  //	 RMS_MURES_ORIG,sqrt(AVG_MUCOV),sigTmp_hi);
  sigTmp = sigTmp_hi;
  while (!BOUND_ONE){
    
    if ( sigTmp < sigint_min ) {
      if ( LABORT ) {
	print_preAbort_banner(fnam);
	printf("  %s called from %s\n", fnam, callFun);
	sprintf(c1err,"Cannot compute sigint because sig trial < %f ??",
		sigint_min );
	sprintf(c2err,"STD=%le, sqrt(AVG_COV)=%le  N=%d",
		STD_MURES_ORIG, sqrt(AVG_MUCOV), N );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;       
      } 
      else { 
	return sigint_min; 
      }
    }
    
    sum_dif = sum_sqdif = sum_pull = sum_sqpull = 0.0 ;
    covTmp = sigTmp * fabs(sigTmp) ;

    for(i=0; i < N; i++ ) {
      covtot = MUCOV_LIST[i] + covTmp;
      if ( covTmp < 0  &&  covtot < covtotfloor ) {
	covtot = covtotfloor;
      }
      pull        = (MURES_LIST[i] - AVG_MURES) / sqrt(covtot);
      sum_pull   += pull ;
      sum_sqpull += ( pull * pull);
    }
    rmsPull = STD_from_SUMS(N, sum_pull, sum_sqpull);
    if ( rmsPull == 0.0 ){
      debugexit("xxx rmsPull = 0");
    }

    if (NBIN_SIGINT < MXSTORE_PULL) {
       rmsPull_store[NBIN_SIGINT] = rmsPull;
       sigTmp_store[NBIN_SIGINT]  = sigTmp ;
    }

    NBIN_SIGINT++ ;

    if (NBIN_SIGINT >= MXSTORE_PULL) {
      print_preAbort_banner(fnam);
      sprintf(c1err,"NBIN_SIGINT=%d exceeds bound MXSTORE_PULL=%d",
	      NBIN_SIGINT,MXSTORE_PULL);
      sprintf(c2err,"Increase bound or check array input");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;             
    }
    
    sigTmp -= sigint_bin;


    if (rmsPull>1.0){ BOUND_ONE = true; }

    //if (sigTmp<sigint_min) { 
    //  sprintf(c1err,"rmsPull > 1 for sigTmp=%f ???",sigint_min);
    //  sprintf(c2err,"called from %s",callFun);
    //  errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
    //}

  } // end sigTmnp loop

  bool ONE_TEST = ONE >= rmsPull_store[0] && ONE <= rmsPull_store[NBIN_SIGINT-1];
  if (!ONE_TEST){ 
    print_preAbort_banner(fnam);
    printf("  called from: '%s' \n", callFun);
    printf("  sigTmp_store range is %f to %f \n", 
	   sigTmp_store[0],sigTmp_store[NBIN_SIGINT-1]);
    printf("  NBIN_SIGINT=%d N_EVT=%d sigint_approx=%f\n", 
	   NBIN_SIGINT, N, sigint_approx );
    printf("  STD_MURES_ORIG=%f sqrt(AVG_MUCOV)=%f\n", 
	   STD_MURES_ORIG, sqrt(AVG_MUCOV) );
    sprintf(c1err,"ONE NOT CONTAINED by rmsPull_store array" );
    sprintf(c2err,"rmsPull_store range is %f to %f",
	    rmsPull_store[0],rmsPull_store[NBIN_SIGINT-1]);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ;
  }
  // interpolate sigInt vs. rmsPull at rmsPull=1                          
  sigint = interp_1DFUN(OPT_INTERP, ONE, NBIN_SIGINT,
			rmsPull_store, sigTmp_store, fnam);

  if ( LDMP ) {
    printf(" xxx %s: sigint(approx->final) = "
	   "%.4f -> %.4f  (<MURES>=%f)\n",
	   fnam, sigint_approx, sigint, AVG_MURES);
    printf(" xxx \n");
    fflush(stdout) ;
  }


  return sigint;

} // end sigint_muresid_list

// =============================================
void remove_quote(char *string) {

  // remove quote(s) from string
  // Input *string is returned without quotes.

  char q[]=  "'" ;
  char qq[] = "\"" ;
  if ( strstr(string,q) == NULL && strstr(string,qq) == NULL ) 
    { return ; }


  int i, i2, lens = strlen(string);
  char *string_orig = (char*) malloc( (lens+1) * sizeof(char) );  
  char ctmp[2] ;
  i2=0;  sprintf(string_orig,"%s", string);  string[0]=0; 
  for(i=0; i<lens; i++ ) {
    sprintf(ctmp, "%c", string_orig[i] ) ;
    if ( strcmp(ctmp,q)  == 0 ) { continue ; }
    if ( strcmp(ctmp,qq) == 0 ) { continue ; }
    sprintf(&string[i2], "%c", string_orig[i] );
    i2++ ;
  }

  free(string_orig);

  return ;

} // end remove_quote

void extractStringOpt(char *string, char *stringOpt) {

  // Created Aug 23 2016                          
  //                                                         
  // For input *string = 'blabla(stringOpt)moreBlaBla' 
  // return 
  //   *string    = 'blablamoreBlaBla' 
  //   *stringOpt = 'stringOpt'
  //  
  // i.e., remove () from *string and return contents of () 
  // in *stringOpt.  If there are no (), then *string is
  // returned unchanged and *stringOpt="".   

  int  i,lens = strlen(string);
  int  L=0, R=0 ;                // Left & Right parentheses logicals 
  char *stringLocal, ctmp[2] ;
  //  char fnam[] = "extractStringOpt";

  // --------------- BEGIN --------------    

  stringOpt[0]=0;
  if ( strstr(string,"(") == NULL ) { return ; }

  stringLocal = (char*) malloc ( (lens+2)*sizeof(char) );
  sprintf(stringLocal, "%s", string);  string[0]=0; ;

  for(i=0; i < lens; i++ ) {
    sprintf(ctmp, "%c", stringLocal[i] ) ;
    if ( *ctmp == '(' ) { L=1; R=0; continue ; }
    if ( *ctmp == ')' ) { R=1; L=0; continue ; }

    if ( L==1 && R==0 )
      { strcat(stringOpt,ctmp); }
    else
      { strcat(string,ctmp); }

  }

  free(stringLocal);
  return ;

} // end extractStringOpt

void extractstringopt_(char *string, char *stringOpt) 
{ extractStringOpt(string,stringOpt); }

// ===============================================
void  remove_string_termination(char *STRING, int LENTOT) {

  // Created Jan 2014 by R.Kessler
  // remove string termination and replace with padding
  // up to length LENTOT. Note that inputs *STRING is modified.
  int ic, LEN;
  // -------------- BEGIN -------------
  LEN = strlen(STRING);
  for(ic=LEN; ic < LENTOT; ic++ )  { STRING[ic] = ' ' ;  }
} 

// ====================================
void trim_blank_spaces(char *string) {

  // April 2013
  // return string without blank spaces.
  // Assume that first blank space after char is end of string
  // Examples:
  //   "BLA1   BLA2   " -> "BLA1".
  //   "   BLA1   BLA2" -> "BLA1".
  // Note that input string is overwritten !
  // 
  // strlen() is NOT used in case there is no '\0' termination.
  //
  // Jan 10 2017: check for termination char '\0' so that it now
  //              works for properly terminated strings with no
  //              extra blank spaces.
  //
  // Mar 13 2019:
  //  Check \r and \n for <CR> or line-feed. 
  //  

  int MXchar, i, FOUNDCHAR, ISCHAR, ISBLANK, ISTERM ;
  char *tmpString, c1[2] ;
  //  char fnam[] = "trim_blank_spaces" ;

  // -------------- BEGIN ---------

  MXchar  = 1000 ;

  // alloate temporary/huge array since we don't know the length
  // and cannot use strlen.
  tmpString = (char*) malloc( sizeof(char) * MXchar ) ;
  tmpString[0]=0;

  if ( strlen(string) == 0 ) { return; } // Aug 2014

  //  printf(" xxx trim string='%s' \n",string); fflush(stdout);

  // transfer non-null chars in string to tmpString
  FOUNDCHAR = 0 ;

  for ( i=0; i < MXchar-1; i++ ) {
    sprintf(c1, "%c", string[i] );
    ISBLANK = ( string[i] == ' '  )  ;
    ISCHAR  = ( ISBLANK == 0  );
    ISTERM  = ( string[i] == '\0' || string[i] == '\n' || string[i]=='\r' ) ;

    /*
    printf(" xxx ----------------- \n");
    printf(" xxx %s: '%s'(%d): ISCHAR=%d  ISBLANK=%d FOUNDCH=%d  ISTERM=%d \n",
	   fnam, string, i, ISCHAR, ISBLANK, FOUNDCHAR, ISTERM);
    printf(" xxx %s: tmpString = '%s'  will; add c1='%s' \n",  
	   fnam, tmpString, c1);
    */

    if ( ISCHAR ) { FOUNDCHAR++ ; } // once set, stays set

    if ( ISBLANK && FOUNDCHAR ) { goto DONE ; }
    if ( ISTERM               ) { goto DONE ; }

    if ( ISCHAR )  { strcat(tmpString,c1); } 
  }


  // overwrite input arg.
 DONE:

  sprintf(string, "%s", tmpString);
  free(tmpString);
  
} // end of trim_blank_spaces

int index_charString(char *c, char *s) {
  // Nov 9 2020
  // Return index of character *c in string *s
  int indx  = -9;
  char *e   = strstr(s,c);
  if ( e != NULL ) { indx = (int)(e-s); }

  return(indx);

} // index_charString

void splitString(char *string, char *sep, int MXsplit,
		 int *Nsplit, char **ptrSplit) {

  // Created July 2016
  //
  // Inputs:
  //    *string  : string to split (preserved)
  //    *sep     : separator, e.g., ',' or ' ' or '+'
  //    MXsplit  : abort if Nsplit >= MXsplit
  //
  // Output :
  //  *Nsplit     : number of split elements return in ptrSplit
  //  **ptrSplit  : array of pointers to split elements
  //
  // Aug 18 2021
  //   remove termination char in case extra blank spaces + <CR> are included.
  // ---------------

  bool ISTERM;
  int LEN, N;
  char *localString, *ptrtok, lastc[2] ;
  char fnam[] = "splitString" ;

  // ------------ BEGIN ---------------
  LEN         = strlen(string);
  localString = (char*) malloc( (LEN+10) * sizeof(char) );
  sprintf(localString, "%s", string);

  // remove termination char
  sprintf(lastc,"%c", localString[LEN-1]);
  ISTERM  = ( lastc[0] == '\0' || lastc[0] == '\n' || lastc[0] == '\r' );
  if ( ISTERM ) { localString[LEN-1] = 0; }

  ptrtok      = strtok(localString,sep) ; // split string
  N=0;

  while ( ptrtok != NULL  ) {
    if ( N < MXsplit ) 
      { sscanf(ptrtok,"%s", ptrSplit[N] );  }
    ptrtok = strtok(NULL, sep);
    N++ ;
  }

  if ( N > MXsplit ) {
    print_preAbort_banner(fnam);  
    printf("  string to split: '%s' \n", string);
    printf("  split separator: '%s' \n", sep);
    sprintf(c1err,"Nsplit = %d ", N );
    sprintf(c2err,"Exceeds bound MXsplit=%d", MXsplit);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }


  *Nsplit = N ; // load output arg
  free(localString);

  return ;

} // end splitString

void splitString2(char *string, char *sep, int MXsplit,
		  int *Nsplit, char **ptrSplit) {

  // ----
  // Use strtok_r ... much faster than strtok, 
  // but input string is destroyed.
  // ----
  //
  // Created July 2016                                
  //                                                                  
  // Inputs:                                                           
  //    *string  : string to split (destroyed)
  //    *sep     : separator, e.g., ',' or ' ' or '+'            
  //    MXsplit  : abort if Nsplit >= MXsplit                    
  //                                   
  // Output : 
  //  *Nsplit     : number of split elements return in ptrSplit 
  //  **ptrSplit  : array of pointers to split elements    
  //
  // Dec 27 2017: avoid <CR> in case of fgets scooping up extra
  //              blank spaces.
  // ---------------                                             

  int   N;
  char *localString, *token ;
  //  char fnam[] = "splitString2" ;

  // ------------ BEGIN ---------------

  localString = string ;
  N=0;

  while (( token = strtok_r(localString, sep, &localString ))) {

    /*
    printf(" xxx %s: token = '%s' L=%d \n", 
    fnam, token, strlen(token));  */

    if ( token[0] != '\0'  && token[0] != '\n' ) {
      if ( N < MXsplit ) { sprintf(ptrSplit[N],"%s", token ); }
      N++ ;
    }
  }
  *Nsplit = N ; // load output arg       

  return ;

}  // end of splitString2

void split2floats(char *string, char *sep, float *fval) {

  // Created Jun 26 2019
  // for *string = 'xxx[sep]yyy,
  // returns fval[0]=xxx and fval[1]=yyy.
  // Example:
  //   Input string   = 1.3,4.6
  //   Output fval[0] = 1.3
  //   Output fval[1] = 4.6
  //
  // Example:
  //   Input string   =  1.3
  //   Output fval[0] =  1.3
  //   Output fval[1] =  1.3
  //
  int Nsplit ;
  char cnum[2][40], *cptr[2];  cptr[0]=cnum[0]; cptr[1]=cnum[1];
  char fnam[] = "split2floats" ;
  // ---------------- BEGIN --------------------

  fval[0] = fval[1] = -9.0 ;

  if ( strstr(string,sep) == NULL ) 
    { sscanf(string, "%f", &fval[0] ); fval[1]=fval[0];  return ;  }
  

  // split the string by the sep input
  splitString(string, sep, 2, &Nsplit, cptr);
  if ( Nsplit != 2 ) {
    sprintf(c1err,"Invalid Nsplit=%d (expected 2)", Nsplit);
    sprintf(c2err,"Input string='%s'  sep='%s' ", string, sep);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err) ; 
  }

  sscanf(cnum[0], "%f", &fval[0] );
  sscanf(cnum[1], "%f", &fval[1] );

  return;
} // end split2floats


// =====================================
void warn_NVAR_KEY(char *fileName) {
  printf("   WARNING: Should remove obsolete NVAR key from %s\n", fileName);
} 

// =====================================
int file_timeDif(char *file1, char *file2) {

  // Mar  29, 2011
  // Return time-stamp difference between file1 and file2.
  // Dif > 0 if  file1 is  more  recent
  // Dif < 0 if  file1 is  older.
  //
  // Aug 18 2020: Allow for .gz in either file.

  char fnam[] = "file_timeDif" ;

  char file1gz[MXPATHLEN], file2gz[MXPATHLEN];
  struct stat statbuf1, statbuf1gz;
  struct stat statbuf2, statbuf2gz;

  int j1, j2, j1gz, j2gz, t1,  t2 ;

  // ------------ BEGIN ----------------

  j1 = stat(file1, &statbuf1); // returns 0 if file exists
  j2 = stat(file2, &statbuf2);

  sprintf(file1gz, "%s.gz", file1);
  sprintf(file2gz, "%s.gz", file2);
  j1gz = stat( file1gz, &statbuf1gz); 
  j2gz = stat( file2gz, &statbuf2gz);

  if ( j1 !=0 && j1gz != 0 ) {
    sprintf(c1err,"istat returned %d  for", j1);
    sprintf(c2err,"file1='%s' ", file1);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( j2 !=0  && j2gz != 0 ) {
    sprintf(c1err,"istat returned %d  for", j2);
    sprintf(c2err,"file2='%s' ", file2);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( j1 == 0 )
    { t1 = statbuf1.st_mtime ; }
  else
    { t1 = statbuf1gz.st_mtime ; }

  if ( j2 == 0 ) 
    { t2 = statbuf2.st_mtime ; }
  else
    { t2 = statbuf2gz.st_mtime ; }

  return t1 - t2 ;

} // end of  file_timeDif

// =====================================
int nrow_read(char *file, char *callFun) {

  // Created Jan 19 2016 by R.Kessler
  // read and return number of rows in ascii *file.
  // *calFun is the name of the calling function to 
  // print in case of error.
  // Useful to estimate malloc size.
  //
  // Dec 29 2017: use open_TEXTgz to read gzipped files.

  int NROW = 0, GZIPFLAG ;
  FILE *fp;
  char line[1000];
  char fnam[] = "nrow_read" ;

  fp = open_TEXTgz(file, "rt", &GZIPFLAG );

  if ( fp == NULL ) {
    sprintf(c1err,"Could not open file '%s'", file);
    sprintf(c2err,"Called from function %s", callFun );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  while ( fgets (line, 1000, fp) !=NULL  ) {  NROW++ ; }
  fclose(fp);

  return(NROW);

} // end nrow_read

// =====================================
int rd2columnFile(
		   char *file   // (I) file to read
		   , int MXROW  // (I) max size of column arrays
		   , int *Nrow  // (O) number of rows read from file
		   , double *column1, double *column2 // (O) column data
		   ) {

  // open  *file and read/return 2 data columns.
  // Abort if Nrow exceeds MXROW or if file does not exist.
  // Returns SUCCESS flag upon completion
  // Skips rows that have comments starting with #, !, %
  //
  // Dec 29 2017: use open_TEXTgz() to allow for gzipped file.
  // Oct 17 2020: skip optional DOCUMENTATION lines

  FILE *fp ;
  int  n, GZIPFLAG, IS_DOCANA=0 ;
  double tmp1, tmp2;  
  char line[MXPATHLEN], tmpline[MXPATHLEN], s1[40], s2[40] ;
  char *ptrtok ;
  char fnam[] = "rd2columnFile" ;

  // -------------- BEGIN ---------

  fp = open_TEXTgz(file, "rt", &GZIPFLAG );
  if ( fp == NULL ) {
    sprintf(c1err,"Could not open file");
    sprintf(c2err,"%s", file);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  *Nrow =  n = 0;

  while ( fgets (line, 180, fp) !=NULL  ) {

    // skip blank lines
    if ( strlen(line) <= 2 ) { continue ; }
    
    // break of tmpline into blank-separated strings
    sprintf(tmpline,"%s", line);
    ptrtok = strtok(tmpline," ");
    sscanf ( ptrtok, "%s", s1 );
    ptrtok = strtok(NULL, " ");

    // skip comment-lines
    if ( commentchar(s1) ) { continue ; }

    // skip DOCUMENTATION lines (Oct 2020)
    if (strcmp(s1,KEYNAME_DOCANA_REQUIRED)  == 0 ) 
      { IS_DOCANA = 1; }
    if (strcmp(s1,KEYNAME2_DOCANA_REQUIRED) == 0 ) 
      { IS_DOCANA = 0; continue; }

    if (IS_DOCANA) { continue; }

    // get 2nd string 
    sscanf ( ptrtok, "%s", s2 );
    ptrtok = strtok(NULL, " ");

    // strip off the numerical values
    sscanf(s1, "%le" , &tmp1 ) ;
    sscanf(s2, "%le" , &tmp2 ) ;

    column1[n] = tmp1 ;
    column2[n] = tmp2 ;

    n++ ;
    if ( n >= MXROW ) {
      sprintf(c1err,"Nrow exceeds MXROW=%d for file", MXROW );
      sprintf(c2err,"%s", file);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }

  } // end of while loop


  *Nrow = n ;
  fclose(fp);

  return SUCCESS ;

}  // rd2columnFile

// define mangle routine for fortran
int rd2columnfile_(char *file, int *MXROW, int *Nrow,
		   double *col1, double *col2) {
  return rd2columnFile(file,*MXROW, Nrow, col1, col2);
}


// ==========================================
int commentchar(char *str) {

  // returns 1 if input string s is a comment character
  // such as #  %  ! 
  // Nov 17 2021: Return true on blank string.

  char c1[2];
  if ( strlen(str) == 0 ) return 1;
  sprintf(c1,"%c", str[0]);

  if ( strcmp(c1,"#") == 0 ) return  1 ;
  if ( strcmp(c1,"!") == 0 ) return  1 ;
  if ( strcmp(c1,"%") == 0 ) return  1 ;
  if ( strcmp(c1,"@") == 0 ) return  1 ;
  return  0;
}

// ============================================
void fillbins(int OPT, char *name, int NBIN, float *RANGE, 
	      float *BINSIZE, float *GRIDVAL ) {

  // Created Octg 16, 2010 by RSK
  // For input NBIN and *RANGE (min,max),
  // return *BINSIZE  and *GRIDVAL array
  // OPT=1 : grid goes from MIN to MAX
  // OPT=2 : grid goes from MIN+BINSIZE/2 to MAX-BINSIZE/2
  //         (i.e, bin-centers like HBOOK histograms)
  //
  // Initial use is for sngridtools.c
  // if *name is NOT null, then print each gridval


  float VALMIN, VALMAX, DIF, xi, VAL, BIN, BINOFF ;
  int i, LDMP ;

  char fnam[] = "fillbins" ;

  // --------------- BEGIN --------------

  if ( NBIN <= 0 ) {
    sprintf(c1err,"invalid NBIN=%d for '%s' ", NBIN, name );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  }

  VALMIN = *(RANGE+0); 
  VALMAX = *(RANGE+1); 
  DIF    =  VALMAX - VALMIN ;

  if ( OPT == 2 || NBIN == 1 ) {
    BIN    =  DIF / (float)(NBIN);
    BINOFF = 0.5;
  }
  else {
    BIN    =  DIF / (float)(NBIN-1);
    BINOFF = 0.0 ;
  }

  *BINSIZE = BIN ;

  LDMP = strlen(name);
  if ( LDMP ) printf(" %2d GRID-BINS for %8s : ", NBIN, name );


  // now fill grid value at each grid center

  for ( i = 0; i < NBIN; i++ ) {
    xi = (float)i;
    VAL = VALMIN + BIN * (xi + BINOFF) ;
    *(GRIDVAL+i) = VAL;

    if ( LDMP ) printf("%6.2f ", VAL);
  }

  if ( LDMP ) printf("\n");


} // end of fillbins



// ***********************************
double getRan_Flat(int ilist, double *range ) {

  // return random number between range[0] and range[1]
  // 'ilist' selects which random list.
  // Jul 12 2021: rename FlatRan -> getRan_Flat

  double xran, dif, rantmp;
  double x0, x1;

  xran     = getRan_Flat1(ilist);
  x0       = range[0] ;
  x1       = range[1] ;
  dif      = x1 - x0;
  rantmp   = x0  + xran * dif ;

  return rantmp;

}  // end of rangen


// ****************************************
double getRan_GaussAsym(double siglo, double sighi, double peakinterval ) {

  // Return random number from bifurcate gaussian
  // with sigmas = "siglo" and "sighi" and peak = 0.0
  // Peak Interval extends from 0 to 0 + peakinterval
  // siglo extends below 0 
  // sighi extends above peakinterval
  //
  // Jan 2012: always pick random number to keep randoms synced.
  // July 6 2021 added peakinterval argument 

  double rr, rg, psum, biran, p[3]   ;
  double BIGAUSSNORMCON = 1.25331413732 ;
    // ---------- BEGIN ------------


  // pick random number to decide which half of the gaussian we are on
  rr = getRan_Flat1(1) ;  // pick random number between 0 and 1

  biran = 0.;

  if ( siglo == 0.0 && sighi == 0.0 ) { return biran;  } 

  psum = (siglo + sighi) * BIGAUSSNORMCON + peakinterval ;
  p[0] = (siglo * BIGAUSSNORMCON) / psum ;
  p[1] = p[0] + peakinterval / psum ; 
  p[2] = p[1] + (sighi * BIGAUSSNORMCON) / psum ; 

  if (rr < p[0]) {
    //useGauss_lo = true ; 
    rg = getRan_Gauss(1) ;
    biran = -fabs(rg) * siglo ; 
  }
  else if (rr < p[1]) {
    //use_peakinterval = true ;
    rg = getRan_Flat1(1) ; 
    biran = (peakinterval * rg) ; 
  }
  else {
    //useGauss_hi = true ; 
    rg = getRan_Gauss(1) ; 
    biran = +fabs(rg) * sighi + peakinterval ; 
  }

  return biran ;

} // end of getRan_GaussAsym

// **********************************************
double biGaussRan_LEGACY(double siglo, double sighi ) { 
 
  // !!!! HEY YOU! Mark Legacy July 2 2021
  // Return random number from bifurcate gaussian
  // with sigmas = "siglo" and "sighi" and peak = 0.0
  //
  // Jan 2012: always pick random number to keep randoms synced.

  double rr, rg, sigsum, plo, biran    ;

    // ---------- BEGIN ------------


  // pick random number to decide which half of the gaussian we are on
  rr = getRan_Flat1(1) ;  // pick random number between 0 and 1

  biran = 0.;

  if ( siglo == 0.0 && sighi == 0.0 ) { return biran;  } 

  sigsum = siglo + sighi ;
  plo    = siglo / sigsum ;    // prob of picking lo-side gaussian

  if ( sigsum <= 0.0 ) { return biran;  }


  // pick random gaussian number; force it positive
  rg = getRan_Gauss(1);
  if ( rg < 0.0 ) { rg = -1.0 * rg ; }  // force positive random

  if ( rr < plo ) 
    { biran = -rg * siglo;  }
  else
    { biran = +rg * sighi;  }


  return biran ;

} // end of biguassran_LEGACY

// **********************************************
double getRan_skewGauss(double xmin, double xmax, 
			double siglo, double sighi, 
			double skewlo, double skewhi) {

  //  Mar 16 2014
  //
  // Select random number from skewed Gaussian defined by
  //
  //      SG(x) = exp [ x^2 / (2 * SIG^2) ]
  // where
  //   SIGLO = siglo + |x|*skewlo
  //   SIGHI = sighi + |x|*skewhi
  //
  //
  // Notes:
  //   *  'skew' is not the usual definition of a skewed Gaussian.
  //   *  xmin, xmax are bounds such that peak of distribution is 
  //      xpeak=0; therefore pass xmin = XMIN-XPEAK and xmax = XMAX-XPEAK.  
  //   * xmin must be negative, and xmax must be positive.
  //   * siglo and sighi must both be positive
  //   * skew can be positive or negative, but if SIGLO or SIGHI
  //     become negative, code aborts.
  //
  // Method:
  // Since we cannot do exact solution, select bi_Gaussian random 'x' 
  // from Bounding BiGaussian Function [BBGF] defined with the skewed
  // sigma = |xmax| or |xmin|,  depending on the sign of skew.
  // Then assign probability Prob = SG(x)/BBGF(x) and verify that
  // 0 <= P <= 1. Keep selecting from BBGF until random(0-1) < Prob.
  //

  double SIGHI_BBGF, SIGLO_BBGF, x, P_BBGF, P_SG, Prob, ran1 ;
  char fnam[] = "getRan_skewGauss" ;
  int NTRY;
  int MXTRY = 100 ;

  // ------------ BEGIN -----------

  // define BBGF = Bounding Bi-Gaussian Function
  if ( skewlo > 0.0 ) 
    { SIGLO_BBGF = fabs(xmin) ; }
  else 
    { SIGLO_BBGF = siglo ; }


  if ( skewhi > 0.0 ) 
    { SIGHI_BBGF = fabs(xmax) ; }
  else 
    { SIGHI_BBGF = sighi ; }


  sprintf(c2err,"BND=%6.3f,%6.3f  sig=%6.3f,%6.3f  skew=%.3f,%.3f", 
	  xmin,xmax, siglo,sighi, skewlo,skewhi );
  


  NTRY = 0 ;
  P_SG = P_BBGF = Prob = ran1 = -9.0 ;

 PICKRAN:
  NTRY++ ;

  /* xxx
  printf(" xxx ------------------------------- \n");
  printf(" xxx NTRY = %d (xmin,xmax=%.3f,%.3f) \n", NTRY, xmin, xmax );
  printf(" xxx x=%.3f  Prob= %.4le/%.4le = %.4le   (ran1=%f) \n",
	 x, P_SG, P_BBGF, Prob, ran1); 
  fflush(stdout);
  xxx */


  if ( NTRY > MXTRY ) {
    sprintf(c1err,"Could not find random from skewed Gaussian after %d tries", 
	    NTRY );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  x = getRan_GaussAsym(SIGLO_BBGF, SIGHI_BBGF, 0.);
  if ( x < xmin ) { goto PICKRAN ; }
  if ( x > xmax ) { goto PICKRAN ; }

  P_BBGF = funVal_skewGauss(x, SIGLO_BBGF, SIGHI_BBGF, skewlo, skewhi );
  P_SG   = funVal_skewGauss(x, siglo,      sighi,      skewlo, skewhi );

  if ( P_BBGF <= 1.0E-9 ) {
    sprintf(c1err,"Invalid P_BBGF = %le for x=%.3f", P_BBGF, x);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  Prob = P_SG / P_BBGF ;
    
  if ( Prob < 0.0 || Prob > 1.0 ) {
    sprintf(c1err,"Invalid P_SG/P_BBGF = %.4le/%.4le = %.4f for x=%.3f", 
	    P_SG, P_BBGF, Prob, x);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // apply weight for this x value

  ran1 = getRan_Flat1(2); // 2-> 2nd list of randoms
  if ( ran1 < Prob ) 
    { return x ; }
  else
    { goto PICKRAN ; }


} // end of getRan_skewGauss

double funVal_skewGauss(double x, double siglo,double sighi, 
			double skewlo, double skewhi ) {

  // March 2014
  // See function definition at top of skewGaussRan().
  // Note that function peak value (mode) must correspond to x=0.

  double sqx, sig, sqsig, arg, funval ;
  char fnam[] = "funVal_skewGauss";

  // ---------------  BEGIN ----------------

  if ( x < 0.0 ) 
    { sig = siglo + skewlo * fabs(x); }
  else
    { sig = sighi + skewhi * x; }

  if ( sig < 0.0 ) {
    sprintf(c1err, "sig = %f < 0 for x=%f", sig, x);
    sprintf(c2err, "sig(lo,hi)=%.3f,%.3f  skew(lo,hi)=%.3f,%.3f",
	    siglo, sighi, skewlo, skewhi);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sqx    =   x * x ;
  sqsig  = sig * sig ;
  arg    = 0.5 * sqx/sqsig ;
  funval = exp(-arg);

  return funval ;

} // end of funVal_skewGauss


// ************************************
void init_random_seed(int ISEED, int NSTREAM) {

  // Create Jun 4 2020 by R.Kessler
  //  [moved init_simRandoms from snlc_sim.c to here, and re-named it]
  //
  // Init random seed(s) 
  // NSTREAM = 1 -> one random stream and regular init with srandom()
  // NSTREAM = 2 -> two independent streams, use srandom_r

  GENRAN_INFO.NSTREAM = NSTREAM ;
  int i ;
  int ISEED2 = ISEED * 7 + 137; // for 2nd stream, if requested
  int ISEED_LIST[MXSTREAM_RAN] = { ISEED, ISEED2 } ;
  char fnam[] = "init_random_seed" ;

  // ----------- BEGIN ----------------

  if ( NSTREAM == 1 ) 
    {   srandom(ISEED); }
  else {

#ifdef ONE_RANDOM_STREAM
    sprintf(c1err,"Invalid NSTREAM_RAN=%d because", NSTREAM);
    sprintf(c2err,"ONE_RANDOM_STREAM pre-proc flag is set in sntools.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
#else
    for(i=0; i < NSTREAM; i++ ) {
      memset( &GENRAN_INFO.ranStream[i], 0,  
	      sizeof(GENRAN_INFO.ranStream[i]) ) ;
      initstate_r(ISEED_LIST[i], GENRAN_INFO.stateBuf[i], BUFSIZE_RAN, 
		  &GENRAN_INFO.ranStream[i] ); 
      srandom_r( ISEED_LIST[i], &GENRAN_INFO.ranStream[i] ); 
    }
#endif
  }

  fill_RANLISTs(); 
  for ( i=1; i <= GENRAN_INFO.NLIST_RAN; i++ )  { 
    GENRAN_INFO.RANFIRST[i]    = getRan_Flat1(i); 
    GENRAN_INFO.NWRAP_MIN[i]   = 99999.0 ;
    GENRAN_INFO.NWRAP_MAX[i]   = 0.0; 
    GENRAN_INFO.NWRAP_SUM[i]   = 0.0; 
    GENRAN_INFO.NWRAP_SUMSQ[i] = 0.0; 
  }

  GENRAN_INFO.NCALL_fill_RANSTATs = 0;
  
  // ---------------- skewNormal stuff -------------------
  //
  // xxx  init_skewNormal(ISEED);  // one-time init, to set seed in python

  return ;

} // end init_random_seed


// **********************************************
void fill_RANLISTs(void) {

  // Dec 1, 2006 RSK
  // Load RANSTORE array with random numbers (uniform from 0-1)
  // from stream 0 using unix_getRan_Flat1(0)
  //
  // Jun 9 2018: use unix_getRan_Flat1() call.
  // Jun 4 2020: change function name from init_RANLIST -> fill_RANLISTs

  int ilist, istore, NLIST_RAN;
  char fnam[] = "fill_RANLISTs" ;

  // ---------------- BEGIN ----------------

  NLIST_RAN = 0 ;
  NLIST_RAN++ ;   // main generation
  NLIST_RAN++ ;   // genSmear
  NLIST_RAN++ ;   // spectrograph
  GENRAN_INFO.NLIST_RAN = NLIST_RAN ;

  if ( NLIST_RAN > MXLIST_RAN ) {
    sprintf(c1err,"NLIST_RAN=%d exceeds bound.", NLIST_RAN);
    sprintf(c2err,"Check MXLIST_RAN = %d", MXLIST_RAN);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sumstat_RANLISTs(0);

  for (ilist = 1; ilist <= NLIST_RAN; ilist++ ) {
    // fill new list of randoms
    GENRAN_INFO.NSTORE_RAN[ilist] = 0 ;
    for ( istore=0; istore < MXSTORE_RAN; istore++ ) {
      GENRAN_INFO.RANSTORE[ilist][istore] = unix_getRan_Flat1(0);
    }
  }

  return ;

}  // end of fill_RANLISTs

// **********************************
void sumstat_RANLISTs(int FLAG) {

  // Created Jun 4 2020
  // called by fill_RANLISTs
  // FLAG = 0 -> increment stats
  // FLAG = 2 -> print final symmary stats

  int NLIST_RAN = GENRAN_INFO.NLIST_RAN ;
  int ilist, NCALL ;
  double XNWRAP, XN1, XN0 = (double)MXSTORE_RAN ;

  // --------- BEGIN ----------

  if ( FLAG == 0 ) {
    if ( GENRAN_INFO.NSTORE_RAN[1] == 0 ) { return; }   
    GENRAN_INFO.NCALL_fill_RANSTATs++ ;
    for (ilist = 1; ilist <= NLIST_RAN; ilist++ ) {      
      // increment stats for previous wrap-usage
      XN1 = GENRAN_INFO.NSTORE_RAN[ilist];
      GENRAN_INFO.NWRAP[ilist]       += (XN1/XN0) ;
      XNWRAP = GENRAN_INFO.NWRAP[ilist];     
      GENRAN_INFO.NWRAP_SUM[ilist]   += XNWRAP;
      GENRAN_INFO.NWRAP_SUMSQ[ilist] += (XNWRAP*XNWRAP);      
      GENRAN_INFO.NWRAP[ilist]        = 0.0 ;
    } // end ilist
  }
  else {
    // compute final AVG and RMS
    double SUM, SUMSQ ;
    for (ilist = 1; ilist <= NLIST_RAN; ilist++ ) { 
      SUM   = GENRAN_INFO.NWRAP_SUM[ilist] ;
      SUMSQ = GENRAN_INFO.NWRAP_SUMSQ[ilist] ;
      NCALL = GENRAN_INFO.NCALL_fill_RANSTATs ;
      GENRAN_INFO.NWRAP_AVG[ilist] = SUM/(double)NCALL ; 
      GENRAN_INFO.NWRAP_RMS[ilist] = STD_from_SUMS(NCALL, SUM, SUMSQ);
    }
  }

  return;

} // end sumstat_RANLISTs

// **********************************
double unix_getRan_Flat1(int istream) {
  // Created Jun 9 2018
  // Input istream is the random stream: 0 or 1
  // Return random between 0 and 1.
  //
  // Jul 30 2020: check pre-proc flag ONE_RANDOM_STREAM

  int NSTREAM = GENRAN_INFO.NSTREAM ;
  int JRAN ;
  char fnam[] = "unix_getRan_Flat1";
  // ------------ BEGIN ----------------
  if ( NSTREAM == 1 )  { 
    JRAN = random(); 
  }
  else {
#ifdef ONE_RANDOM_STREAM
    sprintf(c1err,"Cannot use 2nd random stream because");
    sprintf(c1err,"ONE_RANDOM_STREAM pre-proc flag is set in sntools.h");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
#else
    random_r(&GENRAN_INFO.ranStream[istream], &JRAN); 
#endif
  }

  double r8 = (double)JRAN / (double)RAND_MAX ;  // 0 < r8 < 1 
  return(r8);
}

double unix_getRan_Flat1__(int *istream) 
{ return( unix_getRan_Flat1(*istream) ); }

// ***********************************
double getRan_Gauss(int ilist) {
  // return Gaussian random number using randoms from "ilist",
  // which uses stream 0.
  double R,  V1, V2, FAC, G ;
  // --------------- BEGIN ----------------
 BEGIN:
  V1 = 2.0 * getRan_Flat1(ilist) - 1.0;
  V2 = 2.0 * getRan_Flat1(ilist) - 1.0;
  R  = V1*V1 + V2*V2 ;
  if ( R >= 1.0 ) { goto BEGIN ; }
  FAC = sqrt(-2.*log(R)/R) ;
  G = V2 * FAC ;

  return G ;
}  // end of getRan_Gauss


double unix_getRan_Gauss(int istream) {
  // Created Jun 4 2020
  // pick random Gaussian directly from unix_getRan_Flat1 using 
  // independent "istream" input.
  double R,  V1, V2, FAC, G ;
  int    NSTREAM = GENRAN_INFO.NSTREAM ;
  char fnam[] = "unix_getRan_Gauss" ;

  // --------------- BEGIN ----------------
 BEGIN:
  if ( istream >= NSTREAM ) {
    sprintf(c1err,"Invalid istream = %d (NSTREAM=%d)", istream, NSTREAM);
    sprintf(c2err,"Check call to init_random_seed." );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  V1 = 2.0 * unix_getRan_Flat1(istream) - 1.0;
  V2 = 2.0 * unix_getRan_Flat1(istream) - 1.0;
  R  = V1*V1 + V2*V2 ;
  if ( R >= 1.0 ) { goto BEGIN ; }
  FAC = sqrt(-2.*log(R)/R) ;
  G = V2 * FAC ;
  return G ;
} // end unix_getRan_Gauss

double getRan_GaussClip(int ilist, double ranGmin, double ranGmax ) {
  // Created Aug 2016
  double ranG ;
 PICK_RANGAUSS:
  ranG = getRan_Gauss(ilist);
  if ( ranG < ranGmin || ranG > ranGmax ) { goto PICK_RANGAUSS; }
  return(ranG);

} // end getRan_GaussClip

// *********************************
double getRan_Flat1(int ilist) {

  int  N ;
  double   x8;
  int  NLIST_RAN = GENRAN_INFO.NLIST_RAN ;
  char fnam[] = "getRan_Flat1" ;

  // return random number between 0 and 1
  // Feb 2013: pass argument 'ilist' to pick random list.

  if ( ilist < 1 || ilist > NLIST_RAN ) {
    sprintf(c1err,"Invalid ilist = %d", ilist);
    sprintf(c2err,"Valid ilist is 1 to %d", NLIST_RAN );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // check to wrap around with random list.
  if ( GENRAN_INFO.NSTORE_RAN[ilist] >= MXSTORE_RAN ) { 
    GENRAN_INFO.NSTORE_RAN[ilist] = 0;  
    GENRAN_INFO.NWRAP[ilist] += 1.0 ;
  }

  // use current random in list
  N  = GENRAN_INFO.NSTORE_RAN[ilist] ;
  x8 = GENRAN_INFO.RANSTORE[ilist][N] ;

  // increment for next usage.
  GENRAN_INFO.NSTORE_RAN[ilist]++;  

  return x8;

}  // end of getRan_Flat1


double getran_gauss__(int *ilist) { return getRan_Gauss(*ilist); }
double getran_flat1__(int *ilist) { return getRan_Flat1(*ilist); }


// ********************************************************
double interp_SINFUN(double VAL, double *VALREF, double *FUNREF,
		     char *ABORT_COMMENT) {

  // Created April 11, 2012
  // sin-function interpolation.
  // Interpolate function(val) where VALREF[0] <= val <= VALREF[1]
  // and FUNREF[0-1] are the function values at the edges.
  // Use sin function to interpolate so that derivative=0
  // at the edges.
  // Initial usage is to pass wavelength for 'val' for
  // interpolating functions of intrinsic SN variations.
  //

  double VAL_MEAN, VAL_DIF, FUN_MEAN, FUN_DIF, FUN, ARG, S ;
  double PI = TWOPI/2.0 ;
  char fnam[] = "interp_SINFUN" ;

  // ----------------- BEGIN ----------------

  if ( VAL < VALREF[0] || VAL > VALREF[1] ) {
    sprintf(c1err,"Invalid VAL=%f passed from '%s'", VAL, ABORT_COMMENT);
    sprintf(c2err,"VALREF = %f to %f", VALREF[0], VALREF[1]) ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  VAL_MEAN = ( VALREF[1] + VALREF[0] ) * 0.5 ;
  VAL_DIF  = ( VALREF[1] - VALREF[0] );
  FUN_MEAN = ( FUNREF[1] + FUNREF[0] ) * 0.5 ;
  FUN_DIF  = ( FUNREF[1] - FUNREF[0] );
  
  ARG = PI * (VAL - VAL_MEAN)/VAL_DIF ;
  S   = sin(ARG);
  FUN = FUN_MEAN + ( 0.5 * FUN_DIF * S ) ;

  return FUN ;

} // end interp_SINFUN


// **************************************************************
double interp_1DFUN(
		    int OPT         // (I) 1=linear, 2=quadratic interp
		    ,double val     // (I) interp at this point (e.g., lambda)
		    ,int NBIN       // (I) number of function bins  passed
		    ,double *VAL_LIST  // (I) list of grid values
		    ,double *FUN_LIST  // (I) list of function values
		    ,char *abort_comment // (I) comment in case of abort
		    ) {

  // Created April 2011 by R.Kessler
  // 1D interp-utility.
  // Interpolate FUN_LIST(VAL_LIST) at the point *val.
  // Note that the input function need NOT be specified
  // on a uniform grid.
  // 
  // NBIN=3   => do the old interp8   functionality
  // NBIN > 3 => do the old lamInterp functionality
  //
  // Dec 13 2019: return min val immediately if NBIN==1

  int  IBIN, ibin0, ibin2 ;
  double 
    val0, val1, frac
    ,fun0, fun1, fun
    ,dif0, dif1
    ,*ptrVAL, *ptrFUN
    ;

  char fnam[] = "interp_1DFUN" ;

  // ------------- BEGIN -----------------

  if ( NBIN==1 ) { return(VAL_LIST[0]) ; }

  // do binary search to quickly find which bin contains 'val'
  IBIN = quickBinSearch(val, NBIN,VAL_LIST, abort_comment );

  if ( IBIN < 0 || IBIN >= NBIN-1 ) {
    sprintf(c1err,"quickBinSearch returned invalid IBIN=%d (NBIN=%d)", 
	    IBIN, NBIN );
    sprintf(c1err,"Check '%s'", abort_comment);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  if ( OPT == OPT_INTERP_LINEAR ) {
    val0 = VAL_LIST[IBIN] ;
    val1 = VAL_LIST[IBIN + 1] ;
    fun0 = FUN_LIST[IBIN ] ;
    fun1 = FUN_LIST[IBIN + 1] ;
    frac = (val - val0)/(val1-val0) ;
    fun  = fun0 + frac*(fun1-fun0);
    return fun ;
  }
  else if ( NBIN < 3 ) {
    sprintf(c1err,"Cannot do quadratic interp with %d bins.", NBIN);
    sprintf(c1err,"Check '%s'", abort_comment);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  else if ( OPT == OPT_INTERP_QUADRATIC ) {
    // quadratic interp
    dif0   = val - VAL_LIST[IBIN];
    dif1   = VAL_LIST[IBIN+1] - val;
    
    if ( dif0 <= dif1 && IBIN > 0 )  
      { ibin0 = IBIN - 1; }
    else 
      { ibin0 = IBIN ; }

    // if there are only three bins, then there is no choice
    // for which bins to use.
    if ( NBIN == 3 ) { ibin0 = 0; } 

    // make sure that all three bins are defined
    ibin2 = ibin0+2;
    if ( ibin0 < 0 || ibin2 > NBIN-1 ) {
      sprintf(c1err,"Invalid quadInterp bins: %d to %d", ibin0, ibin2);
      sprintf(c2err,"Must be within defined bins 0 to %d (val=%f)", 
	      NBIN-1, val );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    ptrVAL = &VAL_LIST[ibin0];
    ptrFUN = &FUN_LIST[ibin0];
    fun = quadInterp ( val, ptrVAL, ptrFUN, abort_comment);
    return fun;

  }
  else { 
    sprintf(c1err,"Invalid OPT=%d for '%s'", OPT, abort_comment);
    sprintf(c2err,"%s", "=> Cannot interpolate.") ;
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    return(-999.);
  }

  return(-999.);

} // end of interp_1DFUN


// **************************************************************
// mangle interp function for fortran.
double interp_1dfun__(int *OPT, double *val, int *NBIN
		    ,double *VAL_LIST, double *FUN_LIST
		    ,char *abort_comment ) {
  return interp_1DFUN(*OPT, *val, *NBIN, VAL_LIST, FUN_LIST, abort_comment);
}

// **************************************************************
double quadInterp ( double VAL, double VAL_LIST[3], double FUN_LIST[3],
		    char *abort_comment ) {
 
/***

  Created April 2011 by R.Kessler, 

  Do quadratic interp using the 3 function values.
  It is recommended that
       a_lam[0] <= lambda  <= a_lam[2]

  or you might get pathological solutions
  
  We solve 

     FUN_LIST[0] = A + B ( VAL_LIST[0] - LBAR )^2
     FUN_LIST[1] = A + B ( VAL_LIST[1] - LBAR )^2
     FUN_LIST[2] = A + B ( VAL_LIST[2] - LBAR )^2

 for A, B and LBAR,
 and then convert to

       quadInterp = C0 + C1*lambda + C2*(lambda**2)

 to avoid divide-by-zero when the 2nd order C2=0

 May 5, 2011:
   no need evaluate C[3]. For a little extra speed,
   compute   FUN = A + B(VAL-LBAR)^2

 Dec 30 2013: fix check on pure linear relation by allowing
              slope of E-12

 ****/

   double 
     A, B, LBAR   // defines quadratic function 
     , v0, v1, v2   // three values of input VAL_LIST
     , f0, f1, f2   // three values of input function 
     , sqv0, sqv1, sqv2
     , dprod, fprod  // internal
     , sqdval10, sqdval21, sqdval02
     , dval10,   dval21,   dval02
     , top, dval, dum
     , slope
     , dif
     , fun            // fun(lambda) => return value 
     , C[3]
        ;


   int LDMP;

   /* ------------------------ BEGIN -------------------- */

   f0 = FUN_LIST[0];
   f1 = FUN_LIST[1];
   f2 = FUN_LIST[2];

   v0 = VAL_LIST[0];
   v1 = VAL_LIST[1];
   v2 = VAL_LIST[2];

   // check trivial case of flat function
   if ( f0 == f1 && f1 == f2 ) {  fun = f1; return fun ; }

   // check for pure linear relation 
   if ( f1 != f0 ) {
     dum = (f2-f1)/(f1-f0) ;
     if ( fabs(dum-1.0) < 1.0E-12 ) {  // f2-f1 \simeq f1 - f0
       dval  = v1 - v0 ;
       slope = (f1-f0) / dval ;
       dval  = VAL - v0 ; 
       fun   = f0 + slope * dval ;       
       return fun;     
     }
   }

   // if we get here, do quadratic interpolation

   sqv0 = v0 * v0 ;
   sqv1 = v1 * v1 ;
   sqv2 = v2 * v2 ;

   sqdval10 = sqv1 - sqv0 ;
   sqdval21 = sqv2 - sqv1 ;
   sqdval02 = sqv0 - sqv2 ;

   dval10   = v1 - v0 ; 
   dval21   = v2 - v1 ; 
   dval02   = v0 - v2 ; 
   dprod = dval21 * dval10 * dval02;


   fprod = f0*dval21 + f1*dval02 + f2*dval10;

   if ( dprod == 0.0 ) 
      B = 0;
   else
      B = -fprod / dprod;


   top = f0 * sqdval21 + f1*sqdval02 + f2*sqdval10;
   if ( fprod == 0.0 )
      LBAR = 0.0 ;
   else
      LBAR = 0.5 * top / fprod;
  
   dval = v0 - LBAR;
   A = f0 - B*dval*dval;

   // compute function value

   dif  = VAL - LBAR ;
   fun  = A + B * (dif*dif) ; 

   /*
   C[0] = A + B * LBAR * LBAR;
   C[1] = top / dprod;
   C[2] = B;
   fun  = C[0] + C[1]*VAL + C[2]*(VAL*VAL) ;
   */

   LDMP = 0 ; // ( VAL == 4040.0 ) ; 
   if ( LDMP ) {
     printf(" xxxx -------- NEW interp_1DFUN DUMP ------------------ \n") ;
     printf(" xxxx v0,v1,v2 = %f %f %f \n", v0, v1, v2 );
     printf(" xxxx f0,f1,f2 = %f %f %f \n", f0, f1, f2 );
     printf(" xxxx quadInterp = %f  at val=%f \n", fun, VAL );
     printf(" xxxx C0, C1, C2 = %f %f %f \n", 
	    C[0], C[1], C[2] );
     printf(" xxx A=%le  B=%le  LBAR=%f \n", A, B, LBAR);
     //     debugexit("quad interp");
   }

   return fun ;

} // end of quadInterp


// ===================================================
int quickBinSearch(double VAL, int NBIN, double *VAL_LIST,
		   char *abort_comment) {

  // April 2011.
  // Return integer bin [0 < IBIN < NBIN-1] such that 
  // *(VAL_LIST+IBIN) contains VAL.
  // Use binary search to quickly find IBIN when NBIN is very large.
  //
  // Dec 13 2019: return(0) immediately of NBIN=1

  char fnam[] = "quickBinSearch" ;
  int  LDMP, NITER, ibin_min, ibin_max, ibin, ibin1, ibin2, ISTEP ;
  double    MINVAL, MAXVAL, VAL1, VAL2 ;

  // ------------- BEGIN --------------

  LDMP = 0 ; // ( fabs(VAL-7000.) < 0.01 );

  MINVAL = VAL_LIST[0] ;
  MAXVAL = VAL_LIST[NBIN-1];

  if ( VAL < MINVAL || VAL > MAXVAL )  {
    sprintf(c1err,"VAL = %le outside '%s' range",  
	    VAL, abort_comment );
    sprintf(c2err,"defined for %le < VAL < %le", 
	    MINVAL, MAXVAL );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  NITER    = 0;  // for efficiency testing only
  ibin_min = 0; 
  ibin_max = NBIN-1 ;
  
  if ( NBIN ==  1 ) { return(0); } 

 NEXTSTEP:

  ISTEP    = (ibin_max-ibin_min)/2 ;  // start with just two lambda bins

  // if we are down to fewer than 6 bins to check,
  // skip binary search and just do brute force search
  // one bin at a time.
  if ( ISTEP < 6 ) { ISTEP = 1; }

  /*
  printf("\t xxxx ibin = %3d to %3d  STEPS of %d \n", 
	 ibin_min, ibin_max, ISTEP);
  */

  for ( ibin=ibin_min; ibin <= ibin_max; ibin+=ISTEP ) {

    if ( ibin >= NBIN - 1 ) { continue ; }

    ibin1 = ibin ;
    ibin2 = ibin + ISTEP;
    if ( ibin2 >= NBIN-1 ) { ibin2 = NBIN-1 ; }

    VAL1 = VAL_LIST[ibin1] ;
    VAL2 = VAL_LIST[ibin2] ;
    NITER++ ;

    // abort if NITER gets larger than NBIN
    if ( NITER > NBIN ) {
      print_preAbort_banner(fnam);
      printf("\t ibin1=%d     ibin2=%d  (ISTEP=%d) \n", 
	     ibin1,  ibin2, ISTEP );
      printf("\t VAL1/VAL2 =%6.0f/%6.0f A\n", VAL1,  VAL2);
      printf("\t MINVAL/MAXVAL(entire grid) = %6.0f/%6.0f \n", MINVAL, MAXVAL);

      sprintf(c1err,"Could not find '%s' bin after NITER=%d (NBIN=%d)", 
	      abort_comment, NITER, NBIN);
      sprintf(c2err,"Something is wrong here.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    // just do linear interpolation when we find the right bin
    if ( VAL >= VAL1 && VAL <= VAL2 ) {
      if ( ISTEP == 1 ) 
	{ return ibin1 ; }
      else {
	ibin_min = ibin1 ;
	ibin_max = ibin2 ;
	goto NEXTSTEP ;
      }
    
    }  // VAL if-block
  }  // ibin loop


  // if we get here, then something is really wrong.
  print_preAbort_banner(fnam);
  printf("\t MINVAL=%f  MAXVAL=%f  NBIN=%d\n", MINVAL, MAXVAL, NBIN);
  sprintf(c1err,"Something is REALLY messed up:");
  sprintf(c2err,"Could not find '%s' bin for VAL=%le", abort_comment, VAL );
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  return(-9) ;

} // end of quickBinSearch


// ************************************************
double polyEval(int N, double *coef, double x) {
  // evaluate polynomial function sum_1^N coef[i] * x^i
  // and avoid using the slow pow function.
  // Note that poly order is N-1, not N. N is total number of terms.
  double F=0.0, xpow=1.0 ;
  int i;
  for(i=0; i<N; i++ ) { F += coef[i] * xpow ;  xpow *= x ; }
  return F ;
} // end of polyEval


// ***************************************************
double modelmag_extrap(
		       double T          // Time to extrapolate to
		       ,double Tref      // nearest time with model mag
		       ,double magref   // mag at Tref
		       ,double magslope  // dmag/dT  at Tref
		       ,int LDMP         // print to screen if true
		       ) {

  /****
    Created Apr 12, 2009 by R.Kessler

    Extrapolates model magnitide from Tref to time T
    (obs or rest-frame), where T is outside range of model.
    Note that T < Tref assumes on the rise, and T > Tref 
    assumes on the fall.
    
    Aug 7, 2014: fix bug and return extrapolated mag 
  ****/

  double mag, arg ;

  // -------- BEGIN ----------

  if ( T > Tref ) {

    // linear extrap at late times
    mag = magref + (T-Tref)*magslope ;
  }
  else {

    // exponential fall for early times
    arg = (T - Tref) * magslope/magref ;
    mag = magref * exp(arg);
  }

  return mag ; // Aug 7 2014

} // end of modelmag_extrap


// ***************************************************
double modelflux_extrap(
		       double T          // (I) Time to extrapolate to
		       ,double Tref      // (I) nearest time with model flux
		       ,double fluxref   // (I) flux at Tref
		       ,double fluxslope // (I) dflux/dT  at Tref
		       ,int LDMP         // (I) print to screen if true
		       ) {

  /****
    Created Apr 12, 2009 by R.Kessler

    Extrapolates model flux from Tref to time T
    (rest-frame), where T is outside range of model.
    T and Tref are in days, defined so that T=0 at peak brightness.

    T < Tref => 
    Rises like  F0*(T-T0)^2 x (T-T1) where the quadratic term 
    is a simple fireball model, and the 2nd (linear) term allows 
    matching both the function and its derivative at T=Tref.

    T > Tref => on the tail; use linear mag-extrap, which
                is a power-law flux-extrap.  If fluxref < 0,
                then just hard-wire flux=1E-30 since model
                is unphysical and extrapolation makes no sense.
  
  ****/

  double flux, slope, arg, F0, FP, T1 ;
  double DTref0, SQDTref0, CUBEDTref0 ,DT0, DT1 ;
  double T0  = -20.0 ; // time of explosion, wher T=0 at peak

  //  char fnam[] = "modelflux_extrap" ;

  // -------- BEGIN ----------

  // always use same sign to ensure that flux-> 0 at T->infinity
  slope = fabs(fluxslope) ;

  if ( T > Tref ) {

    // linear mag-extrap at late times => power-law extrap for flux
    if ( fluxref > 0.0 ) {
      arg  = -(T-Tref)*slope / (LNTEN * fluxref) ;
      flux = fluxref * pow(TEN,arg) ;
    }
    else { 
      flux = 1.0E-30 ; 
    }
  }
  else {

    if ( T < T0 ) { flux = 0.0 ; return flux ; }

    // solve for F0 and T1 to describe pre-max rise
    FP         = fluxslope ;
    DTref0     = Tref - T0;
    SQDTref0   = DTref0*DTref0;
    CUBEDTref0 = SQDTref0*DTref0;

    F0  = (FP * DTref0 - 2.*fluxref)/CUBEDTref0 ;
    T1  = Tref - fluxref / (F0*SQDTref0) ;

    DT0 = T - T0 ;
    DT1 = T - T1 ;
    flux = F0 * DT0*DT0 * DT1 ;
  }

  if ( LDMP > 0 ) {
    printf("FLUX_EXTRAP: F(T=%5.1f)=%10.3le", T, flux );
    printf(" [Fref(Tref=%3.0f)=%10.3le & dF/dT=%10.3le] \n",
	   Tref, fluxref, fluxslope );
  }

  return flux ;

} // end of modelflux_extrap



// ====================================================
int rd_sedFlux(
	     char *sedFile         // (I) name of SED file to read
	     ,char *sedcomment     // (I) comment to print in front of SED
	     ,double DAYrange[2]   // (I) rest-frame range (days) to accept
	     ,double LAMrange[2]   // (I) lambda-range to accept (A)
	     ,int    MXDAY         // (I) bound on *DAY_LIST
	     ,int    MXLAM         // (I) bound on *LAM_LIST
	     ,int    OPTMASK       // (I) see OPTMASK bits below
	     ,int    *NDAY         // (O) # of epochs within Trange
	     ,double *DAY_LIST     // (O) list of epochs (days)
	     ,double *DAY_STEP     // (O) day-step
	     ,int    *NLAM         // (O) # of lambda bins within LAMrange
	     ,double *LAM_LIST     // (O) list of lambdas
	     ,double *LAM_STEP     // (O) lambda-step
	     ,double *FLUX_LIST    // (O) flux list
	     ,double *FLUXERR_LIST // (O) flux error list
	    ) {

  /***************
   Created Apr 8, 2009 by R.Kessler
   General function to read SED flux with format

      DAY-0   WAVELENGTH-0   FLUX-0
      DAY-0   WAVELENGTH-1   FLUX-1
      DAY-0   WAVELENGTH-2   FLUX-2
      etc ...
 
   Binning in DAY and WAVELENGTH must be uniform, 
   although there is no explicit check here.

   OPTMASK += 1 --> read FLUXERR (4th column of SEDFILE)
   OPTMASK += 2 --> allow non-uniform DAY bins 

   The return-arg lengths are
     - DAY_LIST has length NDAY
     - LAM_LIST has length NLAM and 
     - FLUX_LIST has length NDAY * NLAM

    To extract
    for ( jflux=0; jflux < NDAY * NLAM ; jflux++ ) 
       iep  = jflux/NLAM
       ilam = jflux - NLAM * iep

       Flux   = *(FLUX_LIST+jflux);
       epoch  = *(DAY_LIST+iep);
       lambda = *(LAM_LIST+ilam);


     HISTORY

  Jan 19 2017: if just one bin, return DAY_STEP=0 or LAM_STEP=0.
               Fix divide-by-zero bug.

  Aug 7 2017: little refactoring and abort if DAYSTEP changes.
  Aug 11 2017: 
     + use splitString2 to simplify reading & parsing
     + read 80 chars per line instead of 180 ... hopefully faster
     + pass OPTMASK arg instead of IERRFLAG. See mask def above
  
  Dec 29 2017: use open_TEXTgz() to read gzipped file.

  May 29 2018: 
   + prepare to suppress DAYs in which FLUX=0 at all wavelenghts.
     See NONZEROFLUX counter. Not implemented

  Jan 4 2019:
   +  move error checking earlier, right after loop.
   + add error check for NDAY=NLAM=0 (maybe catch tabs)
        
  Feb 4 2020:
    + read 120 chars per line instead of 80 (and define MXCHAR_RDFLUX)
    + abort if line length is too long (to avoid corruption)
    + abort if fewer than 3 words are read


  **********/

  FILE *fpsed;

  char txterr[20], line[200], lastLine[200] ;
  //  char *ptrtok, s1[60], s2[60], s3[60], s4[60], tmpline[200] ;
  char *ptrStringVal[MXWORDLINE_FLUX], StringVal[MXWORDLINE_FLUX][40];
  char space[] = " ";
  char fnam[]  = "rd_sedFlux" ;

  double day, lam, day_last, lam_last, lam_expect, flux, fluxerr, XN ;
  double daystep_last, daystep, daystep_dif ;
  double lamstep_last, lamstep, lamstep_dif ;
  int iep, ilam, iflux, ival, LAMFILLED, OKBOUND_LAM, OKBOUND_DAY, NBIN ;
  int  NRDLINE, NRDWORD, GZIPFLAG, FIRST_NONZEROFLUX, NONZEROFLUX, LEN ;
  bool OPT_READ_FLUXERR, OPT_FIX_DAYSTEP;

  // define tolerances for binning uniformity (Aug 2017)
  double DAYSTEP_TOL = 0.5E-3; // tolerance on DAYSTEP uniformity
  double LAMSTEP_TOL = 0.01;   // tolerance on LAMSTEP uniformity

  // ------------- BEGIN -------------

  // init return args
  *NDAY = *NLAM = 0 ;

  // set flags that *NDAY and *NLAM are within bounds
  OKBOUND_DAY = OKBOUND_LAM = 1; 

  // open SED file.

  fpsed = open_TEXTgz(sedFile, "rt", &GZIPFLAG );
  if ( !fpsed  ) {
    sprintf(c1err,"Cannot open SED file: " );
    sprintf(c2err,"  '%s' ", sedFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // check OPTMASK args
  OPT_READ_FLUXERR = false ;      // default: ignore fluxerr column
  OPT_FIX_DAYSTEP  = true;     // default => require fixed bin size

  if ( (OPTMASK & 1) > 0 ) { OPT_READ_FLUXERR = true ; }
  if ( (OPTMASK & 2) > 0 ) { OPT_FIX_DAYSTEP  = false ; }

  if ( OPT_READ_FLUXERR  ) 
    { sprintf(txterr, "and errors"); }
  else
    { txterr[0]=0; }

  printf("  Read  %s  SED %s from : \n    %s \n", 
	 sedcomment, txterr, sedFile );

  fflush(stdout);

  iflux = iep = ilam = 0 ;

  LAMFILLED = NRDLINE = NONZEROFLUX = FIRST_NONZEROFLUX = 0;
  daystep_last = daystep=-9.0;  day_last = -999999. ;
  lamstep_last = lamstep=-9.0;  lam_last = -999999. ;

  for(ival=0; ival < MXWORDLINE_FLUX; ival++ ) 
    { ptrStringVal[ival] = StringVal[ival];   }

  line[0] = lastLine[0] = 0 ;

  while ( fgets (line, MXCHARLINE_FLUX+10, fpsed ) != NULL  ) {
    
    LEN = strlen(line);
    if ( LEN > MXCHARLINE_FLUX ) {
      print_preAbort_banner(fnam);
      printf("   sedFile: %s \n", sedFile);
      printf("   current  line: '%s' \n", line);
      sprintf(c1err,"Line length=%d exceeds bound of %d",
	      LEN, MXCHARLINE_FLUX );
      sprintf(c2err,"Either increase MXCHARLINE_FLUX, or reduce line len");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
    }

    NRDLINE++ ;   

    // skip blank lines
    if ( strlen(line) <= 2 ) { continue ; }
    if ( commentchar(line) ) { continue ; }

    //    splitString2(line, space, MXWORD_RDFLUX,  // input line is destroyed
    //		 &NRDWORD, ptrStringVal ) ;  // returned

    splitString(line, space, MXWORDLINE_FLUX, 
		&NRDWORD, ptrStringVal ) ;  // returned
   
    if ( NRDWORD < 3 ) {
      print_preAbort_banner(fnam);
      printf("   sedFile: %s \n", sedFile);
      printf("   previous line: '%s' \n", lastLine);
      printf("   current  line: '%s' \n", line);
      sprintf(c1err,"NRDWORD = %d, but expected at least 3", NRDWORD);
      sprintf(c2err,"Check dumped lines from file above.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
    }

    sscanf(StringVal[0], "%le" , &day  ) ;
    sscanf(StringVal[1], "%le" , &lam  ) ;
    sscanf(StringVal[2], "%le" , &flux ) ;

    sprintf(lastLine, "%s", line);

    if ( OPT_READ_FLUXERR ) { 
      if ( NRDWORD < 4 ) {
	printf("   sedFile: %s \n", sedFile);
	printf("   previous line: '%s' \n", lastLine);
	printf("   current  line: '%s' \n", line);
	sprintf(c1err,"NRDWORD = %d, but expected at least 4", NRDWORD);
	sprintf(c2err,"to read FLUXERR from 4th coloumn.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
      }
      sscanf(StringVal[3], "%le" , &fluxerr ) ;  
    }

    if ( day < DAYrange[0] ) { continue ; }
    if ( day > DAYrange[1] ) { continue ; }

    if ( lam < LAMrange[0] ) { continue ; }
    if ( lam > LAMrange[1] ) { continue ; }

    if ( day > day_last ) {
      if ( iep >= MXDAY ) 
	{ OKBOUND_DAY = 0 ; }
      else
	{ DAY_LIST[iep] = day; }
      
      // check that each DAYSTEP is the same (Aug 2017)
      if ( iep>0 && OKBOUND_DAY && OPT_FIX_DAYSTEP ) {
	daystep      = DAY_LIST[iep] - DAY_LIST[iep-1] ;
	daystep_dif  = fabs(daystep - daystep_last);
	if ( iep>1 && daystep_dif > DAYSTEP_TOL ) {
	  sprintf(c1err,"daystep[iep=%d] = (%.4f) - (%.4f) = %le",
		  iep, DAY_LIST[iep], DAY_LIST[iep-1], daystep);
	  sprintf(c2err,"daystep_last   = (%.4f) - (%.4f) = %le",
		  DAY_LIST[iep-1], DAY_LIST[iep-2], daystep_last);
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
	}
	daystep_last = daystep ;
      }

      iep++ ; *NDAY = iep ;  NONZEROFLUX = 0 ; 

    /* xxx not now
      // do NOT increment this day if all fluxes are zero.
      // Once a valid flux is found, alway store fluxes, 
      // even if non-zero. Thus only UV region is suppressed if zero.
      if ( NONZEROFLUX || FIRST_NONZEROFLUX ) 
	{ iep++ ; *NDAY = iep ;  NONZEROFLUX = 0 ; }
      else
	{ iflux -= *NLAM; }
    */

    }
    else if ( day < day_last ) {
      sprintf(c1err,"day_last = %f   but day=%f ", day_last, day );
      sprintf(c2err,"New day must increment." ) ;
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
    }

    if ( lam < lam_last )  { ilam = 0; LAMFILLED = 1; }
 
    if ( LAMFILLED == 0 ) {

      // increment NLAM and LAM_LIST on first epoch
      if ( ilam >= MXLAM ) 
	{ OKBOUND_LAM = 0 ; }
      else
	{ LAM_LIST[ilam] = lam; }

      *NLAM = ilam+1 ;

      // check that each LAMSTEP is the same (Aug 2017)
      if ( ilam>0 && OKBOUND_LAM ) {
	lamstep      = LAM_LIST[ilam] - LAM_LIST[ilam-1] ;
	lamstep_dif  = fabs(lamstep - lamstep_last);
	if ( ilam>1 && lamstep_dif > LAMSTEP_TOL ) {
	  sprintf(c1err,"lamstep[ilam=%d] = (%.4f) - (%.4f) = %le",
		  ilam, LAM_LIST[ilam], LAM_LIST[ilam-1], lamstep);
	  sprintf(c2err,"lamstep_last   = (%.4f) - (%.4f) = %le",
		  LAM_LIST[ilam-1], LAM_LIST[ilam-2], lamstep_last);
	  //	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
	  errmsg(SEV_WARN, 0, fnam, c1err, c2err ) ;
	}
	lamstep_last = lamstep ;
      }

    } else {

      // make sure that lambdas repeat exactly
      lam_expect = LAM_LIST[ilam] ;
      if ( lam != lam_expect && OKBOUND_LAM ) {
	print_preAbort_banner(fnam);
	printf("\t DAYrange = %7.2f to %7.2f  (NBIN=%d) \n", 
	       DAYrange[0], DAYrange[1], *NDAY  );
	printf("\t LAMrange = %7.1f to %7.1f  (NBIN=%d) \n", 
	       LAMrange[0], LAMrange[1], *NLAM );
	printf("\t iep=%d  ilam=%d \n", iep, ilam);

	sprintf(c1err,"Found LAM=%6.1f at day = %f ", lam, day);
	sprintf(c2err,"But expected LAM=%7.2f (ilam=%d)", lam_expect, ilam );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
    }

    ilam++ ; 

    // load flux after OKBOUND arrays are set
    if( OKBOUND_DAY && OKBOUND_LAM )  { 
      FLUX_LIST[iflux] = flux ;  
      if ( flux > 1.0E-12 ) { NONZEROFLUX++ ; FIRST_NONZEROFLUX = 1 ;}
      if ( OPT_READ_FLUXERR ) { FLUXERR_LIST[iflux] = fluxerr ;  }
      iflux++ ;  
    }

    day_last = day;
    lam_last = lam;

  }  // end of read-while loop


  // - - - - - - - - - 
  // error checking
  if ( OKBOUND_DAY == 0 ) {
    sprintf(c1err,"NDAY=%d exceeds bound of MXDAY=%d", *NDAY, MXDAY );
    sprintf(c2err,"%s", "Increase EPOCH  binsize or increase MXDAY");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  if ( OKBOUND_LAM == 0 ) {
    sprintf(c1err,"NLAM=%d exceeds bound of MXLAM=%d", *NLAM, MXLAM );
    sprintf(c2err,"%s", "Increase LAMBDA  binsize or increase MXLAM");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // more error checking related to tabs in SED file:
  if ( *NDAY == 0 || *NLAM == 0 ) {
    sprintf(c1err, "Invalid NDAY=%d and NLAM=%d", *NDAY, *NLAM);
    sprintf(c2err, "Check sed file (make sure there are no tabs)");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // - - - - - - - - - 
  NBIN = *NDAY;
  if ( NBIN > 1 ) {
    XN   = (double)(NBIN-1);
    *DAY_STEP = ( DAY_LIST[NBIN-1] - DAY_LIST[0] ) / XN ;
  }
  else 
    { *DAY_STEP = 0.0 ; }


  NBIN = *NLAM ;
  if ( NBIN > 1 ) {
    XN   = (double)(NBIN-1);
    *LAM_STEP = ( LAM_LIST[NBIN-1] - LAM_LIST[0] ) / XN ;
  }
  else
    { *LAM_STEP = 0.0 ; }


  if ( GZIPFLAG==0 ) 
    { fclose(fpsed); } // normal file stream
  else
    { pclose(fpsed); } // for gzip file

  return SUCCESS ;

} // end of rd_sedFlux


// ****************************************************
float effective_aperture ( float PSF_sigma, int VBOSE ) {

  // Returns effective aperture to determine background
  // for PSF defined as Gaussian with sigma = PSF_sigma
  // in pixels.
  // Aperture area = 1 / sum PSF^2
  // PSF = exp(-R^2/2*sigma^2)

  int NPIX, ix, iy ;
  double 
    x, y
    ,arg
    ,PSF, SQPSF
    ,SQR, SQSIG
    ,SUM_PSF, SUM_SQPSF
    ,area
    ,factor
    ;

  // -------------- BEGIN -------------



  // get number of pixels to sum
  NPIX = (int)(20.0 * PSF_sigma) ;
 
  SQSIG = (double)(PSF_sigma * PSF_sigma);

  SUM_PSF       = 0.0 ;
  SUM_SQPSF     = 0.0 ;

  for ( ix = -NPIX; ix <= NPIX; ix++ ) {
    x = (double)ix ;
    for ( iy= -NPIX; iy <= NPIX; iy++ ) {
      y          = (double)iy ;
      SQR        = x*x + y*y ;
      arg        = -0.5 * SQR / SQSIG;
      PSF        = exp(arg) ;
      SQPSF      = PSF * PSF ;
      SUM_PSF   += PSF ;
      SUM_SQPSF += SQPSF ;
    }
  }

  area = SUM_PSF * SUM_PSF / SUM_SQPSF ;

  if ( VBOSE == 1 ) {
    factor = area / (SQSIG * 3.14159265);
    printf(" APERTURE(PSF_sigma=%6.3f) = %10.2f = %8.3f * PI * PSFSIG^2 \n",
	   PSF_sigma, area, factor );
  }

  return area ;

} // end of effective_aperture



// ***************************************************
void clr_VERSION ( char *version, int prompt ) {

  /*****************
   remove files from 

      $SNDATA_ROOT/photometry/version*
      $SNDATA_ROOT/lcmerge/version*
      $SNDATA_ROOT/SIM/version*

     if prompt = 1 then prompt user first before removing.

   Aug 16, 2007: double char-string memory from 100 to 200.

   Oct 7, 2007: include new _SIM path as well.

   Nov 12, 2007: replace VERSION.LIST with VERSION.* to remove
                 .README file as well as .LIST file

  May 29, 2008: remove entire SIM/VERSION directory to avoid
                'arglist too long' for rm * command.

 Sep 30, 2010: replace gets() with scanf() => no more compile warnings !

 Jan 11 2017: remove user-query to remove old version; just do it.

  ***/

  char fnam[] = "clr_VERSION" ;
  char cmd[200];
  char listFile[200];
  char vprefix[200];
  int  RMFILE_PHOTOMETRY, RMFILE_LCMERGE, RMFILE_SIM, DO_RMFILE, isys ;
  FILE *fp;

  // ----------- BEGIN ------------

  print_banner(fnam);

  // check if this version exists by checking listFile

  DO_RMFILE = 0;  // init remove flag to false
  RMFILE_PHOTOMETRY = 0 ;
  RMFILE_LCMERGE = 0 ;
  RMFILE_SIM     = 0;

  sprintf(listFile, "%s/%s.LIST", PATH_SNDATA_LCMERGE, version);
  if ( (fp = fopen(listFile, "rt"))==NULL ) 
    printf("\t LCMERGE Version %s does not exist. \n", version );
  else { 
    printf("\t LCMERGE Version %s exists. \n", version );
    RMFILE_LCMERGE = 1;  DO_RMFILE = 1 ; fclose(fp); 
  }


  sprintf(listFile, "%s/%s.LIST", PATH_SNDATA_SIM, version);
  if ( (fp = fopen(listFile, "rt"))==NULL ) 
    printf("\t SIM Version %s does not exist. \n", version );
  else {
    printf("\t SIM Version %s exists. \n", version );
    RMFILE_SIM = 1;  DO_RMFILE = 1 ; fclose(fp); 
  }


  sprintf(listFile, "%s/%s.LIST", PATH_SNDATA_PHOTOMETRY, version);
  if ( (fp = fopen(listFile, "rt"))==NULL ) 
    printf("\t PHOTOMETRY Version %s does not exist. \n", version );
  else {
    printf("\t PHOTOMETRY Version %s exists. \n", version );
    RMFILE_PHOTOMETRY = 1; DO_RMFILE = 1 ; fclose(fp); 
  }



  if ( DO_RMFILE == 0 ) return ;


  if ( DO_RMFILE == 1 ) {

    printf("\t Removing version %s files ... ", version );

    // make sure to remove ONLY "version" files and not
    // "version*" files in case different versions have
    // similar names.
    // This means that the .LIST and _XXX files are removed separately.

    if ( RMFILE_PHOTOMETRY == 1 ) {
      sprintf(vprefix,"%s/%s", PATH_SNDATA_PHOTOMETRY, version );
      sprintf(cmd,"rm %s.*  %s_*", vprefix, vprefix );
      fflush(stdout);  isys = system ( cmd );
    }

    if ( RMFILE_LCMERGE == 1 ) {
      sprintf(vprefix,"%s/%s", PATH_SNDATA_LCMERGE, version );
      sprintf(cmd,"rm %s.*  %s_*", vprefix, vprefix );
      fflush(stdout);    
      isys = system ( cmd );
    }

    if ( RMFILE_SIM == 1 ) {
      sprintf(cmd,"rm -rf  %s/ ", PATH_SNDATA_SIM );
      printf("\n execute command: %s \n", cmd);
      fflush(stdout);    
      isys = system ( cmd );
    }

    printf(" Done. \n" );
  }

}  // end of clr_VERSION


// ***************************************************
int  init_VERSION ( char *version ) {

  // init VERSION_INFO structure

  // -------------- BEGIN ----------------

  sprintf(VERSION_INFO.NAME, "%s" , version );

  VERSION_INFO.N_SNLC     = 0;
  VERSION_INFO.N_SNFILE   = 0;
  VERSION_INFO.NEPOCH_TOT = 0;

  return SUCCESS ;

}  // end of init_VERSION



// ********************************************
void ld_null (float *ptr, float value) {

  // load *ptr with "value" if *ptr = NULLFLOAT ;
  // If *ptr is not null, then leave it with current value.

  if ( *ptr == NULLFLOAT ) *ptr = value ;

}




// ********************************************
int  init_SNPATH(void) {

  // Feb 19 2015: 
  //   checkg getenv() == NULL rather than if it returns "(null)".
  //   Set SURVEYNAME to "" instead of SDSS.
  //
  // Jan 31 2020: init PATH_USER_INPUT = ""

  char fnam[] = "init_SNPATH" ;

  // ------------ BEGIN ----------


  if ( getenv("SNDATA_ROOT") == NULL ) {
    sprintf(c1err,"Environment variable SNDATA_ROOT not defined.");
    sprintf(c2err,"See SNANA installation instructions.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
  printf (" SNDATA_ROOT = %s \n", PATH_SNDATA_ROOT);


  // ENV is defined, but now check that the directory is really there.
  // Abort if $SNDATA_ROOT/ does not exist
  int istat;
  struct stat statbuf ;
  istat = stat( PATH_SNDATA_ROOT, &statbuf);
  if ( istat != 0 ) {
    sprintf(c1err,"$SNDATA_ROOT = '%s'", PATH_SNDATA_ROOT);
    sprintf(c2err,"does not exist ?!?!? Check your $SNDATA_ROOT .");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  if ( getenv("SNANA_DIR") == NULL ) {
    sprintf(c1err,"Environment variable SNANA_DIR not defined.");
    sprintf(c2err,"See SNANA installation instructions.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
  sprintf( PATH_SNANA_DIR, "%s", getenv("SNANA_DIR") );
  printf (" SNANA_DIR   = %s \n", PATH_SNANA_DIR);

  // - - - - -

  sprintf(PATH_SNDATA_PHOTOMETRY, "%s/photometry", PATH_SNDATA_ROOT);
  sprintf(PATH_SNDATA_LCMERGE,    "%s/lcmerge",    PATH_SNDATA_ROOT);
  sprintf(PATH_SNDATA_SIM,        "%s/SIM",        PATH_SNDATA_ROOT);
  SNDATA.SURVEY_NAME[0]=0;
  SNDATA.SUBSURVEY_NAME[0] = 0 ;
  SNDATA.SUBSURVEY_LIST[0] = 0 ;

  PATH_USER_INPUT[0] = 0 ; 

  fflush(stdout);
  return(SUCCESS);

}   // end of init_SNPATH


// ================================
int init_SNDATA_GLOBAL(void) {

  int ifilt, ep, j ;
  char fnam[] = "init_SNDATA_GLOBAL" ;

  // ---------------- BEGIN -------------

  printf("\n  %s: \n", fnam); fflush(stdout);

  FORMAT_SNDATA_READ  = 0; 
  FORMAT_SNDATA_WRITE = 0;

  SNDATA.SURVEY_NAME[0]    =  0 ;
  SNDATA.MASK_FLUXCOR      =  0 ;
  SNDATA.VARNAME_SNRMON[0] =  0 ;
  SNDATA.DATATYPE[0]       =  0 ;

  SNDATA_FILTER.NDEF       =  0 ;
  SNDATA_FILTER.LIST[0]    =  0 ;
  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    SNDATA_FILTER.MAP[ifilt] = 0;
  }

  SNDATA.NVAR_PRIVATE      = 0 ;  // for data only
  SNDATA.NPAR_SIMSED       = 0 ;    
  SNDATA.NPAR_LCLIB        = 0 ;
  SNDATA.NPAR_PySEDMODEL   = 0 ;
  SNDATA.NPAR_SIM_HOSTLIB  = 0 ;
  
  SNDATA.SIMOPT_MWCOLORLAW = NULLINT ;
  SNDATA.SIMOPT_MWEBV      = NULLINT ;

  SNDATA.SIM_SL_FLAG    = 0 ;
  SNDATA.SIMLIB_FILE[0] = 0 ;
  SNDATA.SIMLIB_MSKOPT  = 0 ;

  SNDATA.APPLYFLAG_MWEBV = 0 ;

  SNDATA.SIM_BIASCOR_MASK = 0 ;
  
  for(ep=0; ep < MXEPOCH; ep++ ) {
   SNDATA.FILTCHAR[ep]  = (char*)malloc( 2  * sizeof(char) );
   SNDATA.FIELDNAME[ep] = (char*)malloc( 20 * sizeof(char) );
  }

  SNDATA.HOSTGAL_NFILT_MAGOBS = 0;
  SNDATA.HOSTGAL_USEMASK      = 0;
  SNDATA.HOSTGAL_NZPHOT_Q    = 0;
  for(j=0; j < MXBIN_ZPHOT_Q; j++)
    { SNDATA.HOSTGAL_PERCENTILE_ZPHOT_Q[j]  = -9.0;  }

  return(SUCCESS);

} // end init_SNDATA_GLOBAL

// *******************************************
int init_SNDATA_EVENT(void) {

  // initialize SNDATA.xxxx elements
  // Note that only one SN index is initialized per call.
  //
  // Mar 13 2021: ZEROPT_ERR[SIG] = 0 instead of -9 in case they are 
  //              not in data files.
  //
  int i_epoch, ifilt, i, igal, j ;
  char fnam[] = "init_SNDATA_EVENT" ;
  // --------- BEGIN -----------------

  /*
  printf(" xxx %s: init SNDATA struct  (SIM_RA=%f)\n",  
	 fnam, SNDATA.SIM_RA ); fflush(stdout);
  */


  sprintf(FLUXUNIT, "ADU");

  sprintf(SNDATA.IAUC_NAME,      "UNKNOWN" );
  sprintf(SNDATA.AUXHEADER_FILE, "UNKNOWN" );

  SNDATA.NLINES_AUXHEADER = 0;
  for ( i=0; i < 20; i++ ) 
    { sprintf(SNDATA.AUXHEADER_LINES[i],"GARBAGE AUXHEADER_LINE"); }

  SNDATA.RA     = NULLFLOAT ;
  SNDATA.DEC    = NULLFLOAT ;
  SNDATA.FAKE   = NULLINT ;
  SNDATA.MWEBV  = NULLFLOAT ;
  SNDATA.WRFLAG_BLINDTEST     = false ; 
  SNDATA.WRFLAG_PHOTPROB      = false ;
  SNDATA.SNTYPE = 0 ;

  SNDATA.FILTCHAR_1D[0] = 0 ;
  SNDATA.FIELDNAME_1D[0] = 0 ;
  SNDATA.NEPOCH = 0;
  SNDATA.NEWMJD = 0;
  SNDATA.MJD_TRIGGER      = 1.0E6 ;
  SNDATA.MJD_DETECT_FIRST = 1.0E6 ;
  SNDATA.MJD_DETECT_LAST  = 1.0E6 ;

  // default mag settings (1/25/2007)
  sprintf(SNDATA.MAGTYPE, "LOG10");
  sprintf(SNDATA.MAGREF,  "AB");

  // init redshift info

  SNDATA.REDSHIFT_HELIO        = NULLFLOAT ;
  SNDATA.REDSHIFT_HELIO_ERR    = NULLFLOAT ;
  SNDATA.REDSHIFT_FINAL        = NULLFLOAT ;
  SNDATA.REDSHIFT_FINAL_ERR    = NULLFLOAT ;
  SNDATA.VPEC = SNDATA.VPEC_ERR = 0.0 ;
  SNDATA.REDSHIFT_QUALITYFLAG  = 0;

  // init HOSTGAL info
  SNDATA.HOSTGAL_NMATCH[0] = 0;
  SNDATA.HOSTGAL_NMATCH[1] = 0;
  SNDATA.HOSTGAL_CONFUSION = -99.0;

  for(igal=0; igal<MXHOSTGAL; igal++ ) {  
    SNDATA.HOSTGAL_OBJID[igal]        =  0 ;
    SNDATA.HOSTGAL_FLAG[igal]         =  0 ;
    SNDATA.HOSTGAL_SPECZ[igal]        = -9.0 ;
    SNDATA.HOSTGAL_SPECZ_ERR[igal]    = -9.0 ;
    SNDATA.HOSTGAL_PHOTOZ[igal]       = -9.0 ;
    SNDATA.HOSTGAL_PHOTOZ_ERR[igal]   = -9.0 ;

    SNDATA.HOSTGAL_SNSEP[igal]        = -9.0 ;
    SNDATA.HOSTGAL_RA[igal]           = -999.0 ;
    SNDATA.HOSTGAL_DEC[igal]          = -999.0 ;
    SNDATA.HOSTGAL_DDLR[igal]         =  -9.0 ;

    SNDATA.HOSTGAL_LOGMASS_TRUE[igal] =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGMASS_OBS[igal]  =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGMASS_ERR[igal]  =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGSFR_TRUE[igal]  =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGSFR_OBS[igal]   =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGSFR_ERR[igal]   =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGsSFR_TRUE[igal] =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGsSFR_OBS[igal]  =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_LOGsSFR_ERR[igal]  =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_COLOR_TRUE[igal]   =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_COLOR_OBS[igal]    =  HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_COLOR_ERR[igal]    =  HOSTLIB_PROPERTY_UNDEFINED ;

    SNDATA.HOSTGAL_SQRADIUS[igal]     = HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_ELLIPTICITY[igal]  = HOSTLIB_PROPERTY_UNDEFINED ;
    SNDATA.HOSTGAL_OBJID2[igal]       = 0 ;
    SNDATA.HOSTGAL_OBJID_UNIQUE[igal] = 0 ;
    for(j=0; j<SNDATA.HOSTGAL_NZPHOT_Q; j++)
      { SNDATA.HOSTGAL_ZPHOT_Q[igal][j] = -9.0; }
  }


  // init SEARCH parameters
  SNDATA.SEARCH_TYPE      = NULLINT ;
  SNDATA.SEARCH_PEAKMJD   = NULLFLOAT ;
  SNDATA.MJD_TRIGGER      = NULLFLOAT ;
  SNDATA.MJD_DETECT_FIRST = NULLFLOAT ;
  SNDATA.MJD_DETECT_LAST  = NULLFLOAT ;
  
  // init sim parameters (used for simulation only)
  sprintf(SNDATA.SIM_MODEL_NAME, "NULL" );
  sprintf(SNDATA.SIM_COMMENT,    "NULL" );
  SNDATA.SIM_MODEL_INDEX    = NULLINT ; // model class; e.g., SIMSED
  SNDATA.SIM_TEMPLATE_INDEX = NULLINT ; // specific template or SED

  SNDATA.SIM_NOBS_UNDEFINED = 0 ;
  SNDATA.SIM_SEARCHEFF_MASK = 0 ;
  SNDATA.SIM_LIBID      = -9 ;
  SNDATA.SIM_NGEN_LIBID =  0 ;


  SNDATA.SIM_REDSHIFT_HELIO = NULLFLOAT ;
  SNDATA.SIM_REDSHIFT_CMB   = NULLFLOAT ;
  SNDATA.SIM_REDSHIFT_HOST  = NULLFLOAT ;
  SNDATA.SIM_REDSHIFT_FLAG        = -9 ;
  SNDATA.SIM_REDSHIFT_COMMENT[0]  =  0 ;

  SNDATA.SIM_DLMU     = NULLFLOAT ;
  SNDATA.SIM_RA       = NULLFLOAT ;
  SNDATA.SIM_DEC      = NULLFLOAT ;
  SNDATA.SIM_PEAKMJD  = NULLFLOAT ;
  SNDATA.SIM_AVTAU    = NULLFLOAT ;
  SNDATA.SIM_AV       = NULLFLOAT ;
  SNDATA.SIM_RV       = NULLFLOAT ;
  SNDATA.SIM_STRETCH  = NULLFLOAT ;
  SNDATA.SIM_DM15     = NULLFLOAT ;
  SNDATA.SIM_DELTA    = NULLFLOAT ;

  SNDATA.SIM_MWRV    = NULLFLOAT ;
  SNDATA.SIM_MWEBV   = NULLFLOAT ;


  SNDATA.SIM_SALT2alpha = NULLFLOAT ;
  SNDATA.SIM_SALT2beta  = NULLFLOAT ;
  SNDATA.SIM_SALT2x0  = NULLFLOAT ;
  SNDATA.SIM_SALT2x1  = NULLFLOAT ;
  SNDATA.SIM_SALT2c   = NULLFLOAT ;
  SNDATA.SIM_SALT2mB  = NULLFLOAT ;

  SNDATA.SIM_MAGSMEAR_COH = 0.0 ;

  SNDATA.SIM_RISETIME_SHIFT = 0.0 ;
  SNDATA.SIM_FALLTIME_SHIFT = 0.0 ;

  SNDATA.SIM_TRESTMIN = 0.0 ;
  SNDATA.SIM_TRESTMAX = 0.0 ;
  SNDATA.SIMFLAG_COVMAT_SCATTER = 0 ;

  SNDATA.SIM_HOSTLIB_GALID = -9 ;

  for ( ifilt=0; ifilt < MXFILTINDX; ifilt++ ) {
    SNDATA.SIM_PEAKMAG[ifilt]      = NULLFLOAT ;
    SNDATA.SIM_TEMPLATEMAG[ifilt]  = NULLFLOAT ;
    SNDATA.SIM_GALFRAC[ifilt]      = NULLFLOAT ;
    SNDATA.NPRESN[ifilt]                = NULLINT ;
    SNDATA.HOSTGAL_SB_FLUXCAL[ifilt]    = -999.0 ;
    SNDATA.HOSTGAL_SB_FLUXCALERR[ifilt] = -999.0 ;
    SNDATA.HOSTGAL_SB_MAG[ifilt]        = 99.0 ;

    for(igal=0; igal<MXHOSTGAL; igal++ ) {
      SNDATA.HOSTGAL_MAG[igal][ifilt]        = 999.0     ;
      SNDATA.HOSTGAL_MAGERR[igal][ifilt]     = 999.0     ;
    }
    SNDATA.SIM_EXPOSURE_TIME[ifilt]     = 1.0  ;
  }

  SNDATA.PIXSIZE     = NULLFLOAT ;
  SNDATA.NXPIX       = -9 ;  
  SNDATA.NYPIX       = -9 ;  

  SNDATA.SUBSAMPLE_INDEX = -9 ;


  //  -------------------------------------------
  // epoch info

  for ( i_epoch = 0; i_epoch < MXEPOCH; i_epoch++ ) {

    SNDATA.SEARCH_RUN[i_epoch]   = NULLINT ;
    SNDATA.TEMPLATE_RUN[i_epoch] = NULLINT ;

    SNDATA.MJD[i_epoch]          = (double)NULLFLOAT ;

    SNDATA.CCDNUM[i_epoch]   = NULLINT ; // Mar 15 2021
    SNDATA.IMGNUM[i_epoch]   = NULLINT ; // Oct 13 2021 

    // Mar 28 2021: replace 'NULL' with 'VOID' because pandas 
    //  gets confused with NULL
    sprintf ( SNDATA.FIELDNAME[i_epoch], FIELD_NONAME );

    SNDATA.IDTEL[i_epoch] = NULLINT ;
    sprintf(SNDATA.TELESCOPE[i_epoch], "BLANK" );


    SNDATA.FILTINDX[i_epoch]       = NULLINT ;
    SNDATA.SEARCH_FIELD[i_epoch]   = NULLINT ;
    SNDATA.TEMPLATE_FIELD[i_epoch] = NULLINT ;

    SNDATA.GAIN[i_epoch]           = NULLFLOAT ;
    SNDATA.READNOISE[i_epoch]      = NULLFLOAT ;

    SNDATA.XPIX[i_epoch]         = NULLFLOAT ;
    SNDATA.YPIX[i_epoch]         = NULLFLOAT ;
    SNDATA.EDGEDIST[i_epoch]     = NULLFLOAT ;

    SNDATA.SKY_SIG[i_epoch]      = NULLFLOAT ;
    SNDATA.SKY_SIG_T[i_epoch]    = 0.0       ; // Mar 26 2021
    SNDATA.PSF_SIG1[i_epoch]     = NULLFLOAT ;
    SNDATA.PSF_SIG2[i_epoch]     = 0.0 ;
    SNDATA.PSF_RATIO[i_epoch]    = 0.0 ;
    SNDATA.PSF_NEA[i_epoch]      = NULLFLOAT ;

    SNDATA.FLUXCAL[i_epoch]         = NULLFLOAT ;
    SNDATA.FLUXCAL_ERRTOT[i_epoch]  = NULLFLOAT ;

    SNDATA.MAG[i_epoch]           = NULLFLOAT ;
    SNDATA.MAG_ERR[i_epoch]       = NULLFLOAT ;

    SNDATA.SKYSUB_ERR[i_epoch]    = NULLFLOAT ;
    SNDATA.GALSUB_ERR[i_epoch]    = NULLFLOAT ;

    SNDATA.ZEROPT[i_epoch]         = NULLFLOAT ;
    SNDATA.ZEROPT_ERR[i_epoch]     = 0.0; 
    SNDATA.ZEROPT_SIG[i_epoch]     = 0.0; 

    SNDATA.PHOTFLAG[i_epoch]       = 0   ;
    SNDATA.PHOTPROB[i_epoch]       = 0.0 ;

    SNDATA.SIMEPOCH_MAG[i_epoch] = 99.0 ;

    SNDATA.SIMEPOCH_WARPCOLNAM[i_epoch][0] = 0 ;
    SNDATA.SIMEPOCH_KCORNAM[i_epoch][0]    = 0 ;
    SNDATA.SIMEPOCH_MAGSMEAR[i_epoch]      = 0.0 ;

  }  //  end i_epoch init loop

  return SUCCESS ;

}   // end of init_SNDATA_EVENT

// =================================================
void init_GENSPEC_GLOBAL(void) {
  int ispec;
  GENSPEC.NMJD_PROC = 0 ;
  GENSPEC.USE_WARP  = 0 ;
  for(ispec=0; ispec < MXSPECTRA; ispec++ )  { 
    GENSPEC.NBLAM_VALID[ispec] = 0;
    GENSPEC.SKIP[ispec]        = false;
  } 
  return;
} // end init_GENSPEC_GLOBAL

void init_GENSPEC_EVENT(int ispec, int NBLAM) {

  
  char fnam[] = "init_GENSPEC_EVENT";

  //  printf(" xxx %s: ispec=%d NB=%d  last NB=%d\n",
  //	 fnam, ispec, NBLAM, GENSPEC.NBLAM_VALID[ispec]  );

  if  ( ispec < 0 ) { init_GENSPEC_GLOBAL(); return; }
  
  if ( GENSPEC.NBLAM_VALID[ispec] > 0 ) {
    free(GENSPEC.LAMMIN_LIST[ispec])  ;
    free(GENSPEC.LAMMAX_LIST[ispec])  ;
    free(GENSPEC.LAMAVG_LIST[ispec])  ;
    free(GENSPEC.FLAM_LIST[ispec])    ;
    free(GENSPEC.FLAMERR_LIST[ispec]) ;
    free(GENSPEC.FLAMWARP_LIST[ispec]) ; // Apr 2 2021
    free(GENSPEC.GENFLAM_LIST[ispec]) ;
    free(GENSPEC.GENMAG_LIST[ispec])  ;
  }

  GENSPEC.NBLAM_VALID[ispec] = NBLAM;
  int MEMD = NBLAM * sizeof(double) ;
  GENSPEC.LAMMIN_LIST[ispec]   = (double*) malloc(MEMD);
  GENSPEC.LAMMAX_LIST[ispec]   = (double*) malloc(MEMD);
  GENSPEC.LAMAVG_LIST[ispec]   = (double*) malloc(MEMD);
  GENSPEC.FLAM_LIST[ispec]     = (double*) malloc(MEMD);
  GENSPEC.FLAMERR_LIST[ispec]  = (double*) malloc(MEMD);
  GENSPEC.FLAMWARP_LIST[ispec] = (double*) malloc(MEMD);
  GENSPEC.GENFLAM_LIST[ispec]  = (double*) malloc(MEMD);
  GENSPEC.GENMAG_LIST[ispec]   = (double*) malloc(MEMD);
  
  // Jul 1 2021: init wave-dependent arrays in case they aren't filled
  int ilam;
  for(ilam=0; ilam < NBLAM; ilam++ ) {
    GENSPEC.FLAM_LIST[ispec][ilam]     = -9.0 ;
    GENSPEC.FLAMERR_LIST[ispec][ilam]  = -9.0 ;
    GENSPEC.FLAMWARP_LIST[ispec][ilam] =  1.0 ;
    GENSPEC.GENFLAM_LIST[ispec][ilam]  = -9.0 ;
    GENSPEC.GENMAG_LIST[ispec][ilam]   = -9.0 ;
  }

} // end init_GENSPEC_EVENT


// =====================================
void set_SNDATA_FILTER(char *filter_list) {
  // Created Feb 15 2021
  // Restore SNDATA_FILTER struct using input *filter_list (e..g, 'ugriz')
  // This function should be called after reading FILTERS arg
  // from data file.
  int NFILT = strlen(filter_list);
  int ifilt, ifilt_obs;
  char cfilt[2];

  set_FILTERSTRING(FILTERSTRING);
  SNDATA_FILTER.NDEF = NFILT;
  sprintf(SNDATA_FILTER.LIST, "%s", filter_list);
  for(ifilt=0; ifilt < NFILT; ifilt++ ) {
    sprintf(cfilt, "%c", filter_list[ifilt] );
    ifilt_obs = INTFILTER(cfilt);
    SNDATA_FILTER.MAP[ifilt] = ifilt_obs; 
  }
  return ;
} // end set_SNDATA_FILTER

// ******************************************************
int  fluxcal_SNDATA ( int iepoch, char *magfun, int opt ) {


  /*********
    fill SNDATA.FLUXCAL[iepoch][ifilt] 
               and
         SNDATA.FLUXCAL_ERRTOT[iepoch][ifilt] 


   = FLUX * 10**{-0.4*ZP}    if  magfun = "log10"

   = SINH (  )               if  magfun = "asinh"


  Inputs:
     iepoch : epoch index for SNDATA.XXX arrays
     magfun : log10 or asinh
     opt    : 0 or 1 -> add ZP_sig uncertainty
                   2 -> do not add ZP_sig term (error added earlier)

           HISTORY 
        ~~~~~~~~~~~~
  Sep 5 2016: magfun = asinh is now obsolete

  Jan 3 2018: check for saturation
  Jan 23 2020: 
    pass opt argument to control use of ZP_sig in flux uncertainty.
   
  *********/

  double mag, mag_err, mag_tmp, ZP, ZP_err, ZP_scale, ZP_sig;
  double flux, flux_err, fluxcal, fluxcal_err, ferrp, ferrm;
  double tmperr, arg, sqerrtmp, relerr;
  int VALID_MAGFUN, IFILT, LTMP;
  char fnam[] = "fluxcal_SNDATA" ;

  // ------------- BEGIN ----------------

  VALID_MAGFUN = 0;

  ZP        = SNDATA.ZEROPT[iepoch];
  ZP_err    = SNDATA.ZEROPT_ERR[iepoch];
  ZP_sig    = SNDATA.ZEROPT_SIG[iepoch];
  flux      = SNDATA.FLUX[iepoch];
  flux_err  = SNDATA.FLUX_ERRTOT[iepoch];
  IFILT     = SNDATA.FILTINDX[iepoch]; // absolute obs filter indx

  fluxcal     = NULLDOUBLE ;
  fluxcal_err = NULLDOUBLE ;

  mag = 99.0 ; mag_err=99.0 ;

  // make some idiot checks

  if ( flux_err <   0.  ) { return ERROR; }

  // Jan 2018 : for saturated epoch, set fluxcal = flux
  if ( SNDATA.NPE_ABOVE_SAT[iepoch] > 0 ) {
    SNDATA.FLUXCAL[iepoch]        = flux ;
    SNDATA.FLUXCAL_ERRTOT[iepoch] = flux_err ;
    return(SUCCESS) ;
  }

  // for asinh mags, get calibrated flux from mag

  if ( strcmp(magfun,"asinh") == 0 ) {

    VALID_MAGFUN = 1 ;

    fluxcal = asinhinv ( mag, IFILT ) ;

    mag_tmp = mag - mag_err ;
    ferrp  = asinhinv ( mag_tmp, IFILT) - fluxcal;

    mag_tmp    = mag+mag_err ;
    ferrm      = fluxcal - asinhinv ( mag_tmp, IFILT ) ;

    fluxcal_err = (ferrp+ferrm)/2.0;

    if ( fluxcal_err > 10000. ) { fluxcal_err = 9999. ; }

    // xxxxxxxxxxxxxx dump + and - errors for test
    double mjd ;
    mjd = SNDATA.MJD[iepoch] ;
    LTMP = 0;  if ( fabs(mjd - 53626.37 ) < 0.1 ) { LTMP = 1; }
    if ( LTMP == -9 ) {
      printf(" MJD=%9.3f  magerr=%5.2f  "
	     "fluxerr(%c) = avg(+%5.1f -%5.1f) = %5.1f\n",
	     mjd, mag_err, FILTERSTRING[IFILT],  ferrp, ferrm, fluxcal_err );
    }
    // xxxxxxxxxxxx

  }

  // for log10 mag, get calibrated flux from FLUXCAL

  if ( strcmp(magfun,"log10") == 0 ) {
    VALID_MAGFUN = 1 ;
    arg      = -0.4 * (ZP - ZEROPOINT_FLUXCAL_DEFAULT) ;
    ZP_scale = pow(TEN,arg) ;
 
    if ( flux_err >= 0.0 && ZP > 10.0 ) {

      fluxcal     = flux     * ZP_scale ;
      fluxcal_err = flux_err * ZP_scale ;

      if (opt <= 1 ) {
	// add ZP_sig uncertainty here for FLUXCAL; 
	// note that flux in ADU does not have this ZP error.
	relerr      = powf(TEN,0.4*ZP_sig) - 1.0 ;
	tmperr      = fluxcal * relerr ;
	sqerrtmp    = fluxcal_err*fluxcal_err + tmperr*tmperr ;
	fluxcal_err = sqrt(sqerrtmp) ;
      }

    }
  }

  if ( VALID_MAGFUN == 0 ) {
    sprintf(c1err,"Invalid mag function: %s ", magfun );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
  }


  // load structure

  SNDATA.FLUXCAL[iepoch]        = fluxcal ;
  SNDATA.FLUXCAL_ERRTOT[iepoch] = fluxcal_err ;

  if ( IFILT == -666 ) 
    printf("Ep=%d flt=%d, ZP = %f +- %f, fluxcal=%f +- %f \n",
	   iepoch, IFILT, ZP, ZP_sig, fluxcal, fluxcal_err );


  return(SUCCESS) ;

}  // end of  fluxcal_SNDATA



// ******************************************* 
double asinhinv(double mag, int ifilt) {

  // Invert SDSS mag to get calibrated flux.

  // define SDSS softening parameter vs. filter
  double bb[6]     = {0.0, 1.4e-10, 0.9e-10, 1.2e-10, 1.8e-10, 7.4e-10};
  double arg, b, fluxCal, magoff, fluxScale;

    // ------------ BEGIN ------------

  magoff  = -2.5/log(10.0);
  b       = bb[ifilt];
  arg     = mag/magoff - log(b);

  fluxScale = pow(TEN,0.4*ZEROPOINT_FLUXCAL_DEFAULT);
  fluxCal   = fluxScale * (2*b) * sinh(arg);

  // xxx delete Mar 2013:  fluxcal = FLUXCAL_SCALE * 2*b * sinh(arg);
 
  return fluxCal;

} // end of asinhinv



// *****************************************
int PARSE_FILTLIST (char *filtlist_string, int *filtlist_array ) {

  // June 27 2008
  // parse input string 'filtlist_string' and return integer-array
  // "filtlist_array" of absolute filter indices.
  // Function arg is the number of filters.
  // 
  // Note that filters can be comma-separated:
  //  e.g.,   ugri and u,g,r,i return same integer list.
  //
  // Jun 22 2021: allow comma sep list ; refactor to use INTFILTER func.
  //

  int  LENLIST = strlen(FILTERSTRING);
  char fnam[] = "PARSE_FILTLIST";

  int NF_USER ;
  int ifilt_user, ifilt_list, ifilt_match;
  char cfilt_user[2], cfilt_tmp[2], **cfilt_list ;
  char filtlist_local[MXFILTINDX];

  //------------- BEGIN ----------

  filtlist_local[0] = 0 ;

  if ( strstr(filtlist_string,COMMA) != NULL )  {
    // e.g., g,r,i,z
    parse_commaSepList(fnam, filtlist_string, MXFILTINDX, 10,
		       &NF_USER, &cfilt_list); // <== returned
    for(ifilt_user=0; ifilt_user<NF_USER; ifilt_user++ ) {
      if ( strlen(cfilt_list[ifilt_user]) > 1 ) {
	sprintf(c1err,"Invalid band = '%s' in comma-sep list '%s' ", 
		cfilt_list[ifilt_user], filtlist_string );
	sprintf(c2err,"Must use single char band between commas.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      strcat(filtlist_local,cfilt_list[ifilt_user]); 
      free(cfilt_list[ifilt_user]);
    }
    free(cfilt_list);
  }
  else {
    // e.g., griz
    NF_USER = strlen(filtlist_string);
    sprintf(filtlist_local,"%s", filtlist_string);
  }
  // - - - - 

  if ( NF_USER >= LENLIST ) {
    sprintf(c1err,"%d filters is too many (> %d)", NF_USER, LENLIST);
    sprintf(c2err,"Check filter list  '%s' ", filtlist_string);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // loop over all user-defined filters to make sure
  // that they are all defined in FILTERSTRING.

  for ( ifilt_user = 0 ; ifilt_user < NF_USER; ifilt_user++ ) {
    sprintf(cfilt_user, "%c", filtlist_local[ifilt_user] ) ;
    ifilt_match = INTFILTER(cfilt_user);

    if ( ifilt_match < 0 ) {
      sprintf(c1err,"User-defined filter '%s' is not defined.", cfilt_user);
      sprintf(c2err,"Check filter list  %s' ", filtlist_string);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    filtlist_array[ifilt_user] = ifilt_match ;
  } // end ifilt_user

  return NF_USER ;

}  // end of function PARSE_FILTLIST


// ************************************************
void checkArrayBound(int i, int MIN, int MAX, 
		     char *varName, char *comment, char *funName) {

  // One-line call to check array bound.
  // Abort if index 'i' is outside index range MIN to MAX.
  // *varName, *comment and *fun are used to construct 
  // error message. Note that *fun is the name of the calling function.
  //

  char fnam[] = "checkArrayBound" ;

  // ------------------ BEGIN ---------------

  if ( (i >= MIN) && (i <= MAX) ) { return ; }

  sprintf(c1err,"%s = %d is outside valid range %d - %d (fun=%s) \n", 
	  varName ,i, MIN, MAX, funName);
  sprintf(c2err,"%s", comment);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  

} // end of checkArrayBound

void  checkArrayBound_(int *i, int *MIN, int *MAX, 
		       char *varName, char *comment, char *funName) {
  checkArrayBound(*i,*MIN,*MAX, varName, comment, funName) ;
}


// ******************************************************
void read_VARNAMES_KEYS(FILE *fp, int MXVAR, int NVAR_SKIP, char *callFun, 
			int *NVAR, int *NKEY,
			int *UNIQUE, char **VARNAMES ) {

  // Mar 2019
  // Read file (*fp) and return information about "VARNAMES: <varList>"
  //
  // Inputs
  //   fp      : file pointer to read
  //    MXVAR  : max number of variables to return after VARNAMES keys.
  //    NVAR_SKIP : number of variables to skip:
  //          postive  -> remove from end of list
  //          negative -> remove from start of list
  //  *callFun : name of calling function; for error message only
  //
  // Output:
  //   *NVAR     : total number of variables after all VARNAMES keys
  //   *NKEY     : total number of VARNAMES keys found
  //   *UNIQUE   : for each variable, 1=> unique, 0=> duplicate from previous
  //  **VARNAMES : list of all variables (0 to *NVAR-1)
  //
  // Jun 12 2020: new option for NVAR_SKIP < 0 -> skip first var(s).
  // Jul 13 2020: avoid storing duplicate var names

  int  NVAR_STORE = 0, NKEY_LOCAL = 0 ;
  int  FOUND_VARNAMES, ivar, ivar_start, ivar_end, ivar2, NVAR_TMP ;
  int  IVAR_EXIST, LDMP = 0 ;
  char c_get[60], LINE[100], tmpName[60] ;
  char fnam[] = "read_VARNAMES_KEYS" ;

  // -------------- BEGIN ------------

  if ( LDMP ) { printf(" xxx %s: ---------- DUMP ------------ \n", fnam);  }

  while( (fscanf(fp, "%s", c_get )) != EOF) {
    FOUND_VARNAMES = ( strcmp(c_get,"VARNAMES:")==0 ) ;
    if ( FOUND_VARNAMES ) {
      NKEY_LOCAL++ ;
      fgets(LINE, 100, fp ); // scoop up varnames
      NVAR_TMP  = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,LINE);
      
      ivar_start = 0; ivar_end = NVAR_TMP;
      if ( NVAR_SKIP < 0 ) { ivar_start -= NVAR_SKIP; }
      if ( NVAR_SKIP > 0 ) { ivar_end   -= NVAR_SKIP; }

      if(LDMP) 	{ 
	printf(" xxx %s: NVAR_TMP=%d ivar[start,end]=%d,%d \n", 
	       fnam, NVAR_TMP, ivar_start, ivar_end);  fflush(stdout);
      }
		       
      for ( ivar=ivar_start; ivar < ivar_end; ivar++ ) {
	get_PARSE_WORD(0, ivar, tmpName);
	IVAR_EXIST = ivar_matchList(tmpName, NVAR_STORE, VARNAMES );
	if ( LDMP ) {
	  printf(" xxx %s: ivar=%d(%s) EXIST=%d  \n",
		 fnam, ivar, tmpName, IVAR_EXIST );
	}
	if ( NVAR_STORE < MXVAR  &&   IVAR_EXIST < 0 ) {
	  if(LDMP) {
	    printf(" xxx %s: \t LOAD VARNAMES[%d] = %s \n",
		   fnam, NVAR_STORE, tmpName); fflush(stdout);
	  }
	  sprintf(VARNAMES[NVAR_STORE], "%s", tmpName);
	  NVAR_STORE++ ; // summed overal all VARNAMES
	}

      }
    } // end FOUND_VARNAMES
  } // end while    



  if ( NVAR_STORE >= MXVAR ) {
    sprintf(c1err,"NVAR_STORE=%d exceeds MXVAR=%d", NVAR_STORE, MXVAR);
    sprintf(c2err,"called by %s, VARNAMES[0]=%s", callFun, VARNAMES[0] );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // ------------------------------------
  // load output args
  *NVAR = NVAR_STORE ;
  *NKEY = NKEY_LOCAL ;

  // check which VARNAMES are unique
  char *NAME, *NAME2 ;
  for(ivar=0; ivar < NVAR_STORE; ivar++ ) {
    UNIQUE[ivar] = 1;
    NAME = VARNAMES[ivar] ;

    // loop over previous VARNAMES to check for duplicate
    for(ivar2=0; ivar2 < ivar; ivar2++ ) {
      NAME2 = VARNAMES[ivar2];
      if ( strcmp(NAME,NAME2) == 0 ) { UNIQUE[ivar]=0; } // duplicate
    }
  }

  
  return ;

} // end read_VARNAMES_KEYS

// **********************
void check_uniform_bins(int NBIN,double *VAL_ARRAY, char *comment_forAbort) {

  // July 2016: abort on non-uniform bin in *VAL array

  int i, j ;
  double VAL, VAL_LAST, DIF, DIF_LAST ;
  char fnam[] = "check_uniform_bins" ;

  // ------------ BEGIN ------------
  
  VAL_LAST = DIF_LAST = 0.0 ;
  for(i=0; i < NBIN; i++ ) {

    VAL = VAL_ARRAY[i];
    DIF = VAL - VAL_LAST ;

    if ( i > 1 && DIF != DIF_LAST ) {
      print_preAbort_banner(fnam);
      printf(" abort comment: %s \n", comment_forAbort);
      for(j=i-2; j < i+3; j++ ) {
	if ( j<0 || j >= NBIN ) { continue ; }
	printf("\t VAL_ARRAY[%d] = %f  (DIF=%f) \n", 
	       j, VAL_ARRAY[j], VAL_ARRAY[j]-VAL_ARRAY[j-1] );
      }
      sprintf(c1err,"Detected non-uniform grid (%s)", comment_forAbort);
      sprintf(c2err,"DIF=%f but DIF_LAST=%f", DIF, DIF_LAST);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    VAL_LAST = VAL;
    DIF_LAST = DIF ;
  }
  
  return ;

} // end check_uniform_bins

// ******************************************************
void  check_magUndefined(double mag, char *varName, char *callFun) {

  // Created Jun 23 2016
  char fnam[] = "check_magUndefined" ;

  if ( mag == MAG_UNDEFINED ) {
    sprintf(c1err,"Undefined %s = %f", varName, mag);
    sprintf(c2err,"Found in function %s", callFun );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );   
  }
}  // end check_magUndefined


// ******************************************************
void checkval_I(char *varname, int nval, int *iptr, int imin, int imax){

  // check that all iptr values are between imin and imax;
  // if not then abort. 
  // Jan 11 2017: isnan(ival) -> isnan ( float)ival )

  int i;
  int  ival ;
  char fnam[] = "checkval_I" ;

  // ----------------- BEGIN -----------------

  for ( i = 0; i < nval;  i++ ) {
    ival = *(iptr+i);
    
    if ( isnan( (float)ival ) ) {
      sprintf(c1err,"%s(item %d) = nan", varname, i );
      sprintf(c2err,"Expected  %i < %s < %i", imin, varname, imax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( ival < imin || ival > imax  ) {
      sprintf(c1err,"%s(item %d) = %d fails bounds-check", varname, i, ival );
      sprintf(c2err,"Expected  %d < %s < %d", imin, varname, imax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

}  // end of checkval_I


void checkval_i__(char *varname, int *nval, int *iptr, int *imin, int *imax){
  checkval_I(varname, *nval, iptr, *imin, *imax) ;
}


// ******************************************************
void checkval_F(char *varname, int nval, float *fptr, float fmin, float fmax){

  // check that all fptr values are between fmin and fmax;
  // if not then abort. 
  // Nov 2013: check for nan and rename function checkval -> checkval_F
  // May 27 2014: isnanf -> isnan (recommended by S. Rodney)

  int i;
  float val ;
  char fnam[] = "checkval_F" ;

  // ----------------- BEGIN -----------------

  for ( i = 0; i < nval;  i++ ) {
    val = *(fptr+i);
    
    if ( isnan(val) ) {
      sprintf(c1err,"%s(item %d) = nan", varname, i );
      sprintf(c2err,"Expected  %f < %s < %f", fmin, varname, fmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( val < fmin || val > fmax  ) {
      sprintf(c1err,"%s(item %d) = %f fails bounds-check", varname, i, val );
      sprintf(c2err,"Expected  %f < %s < %f", fmin, varname, fmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

}  // end of checkval_F

// ******************************************************
void checkval_D(char *varname, int nval, 
		double *dptr, double dmin, double dmax){

  // check that all fptr values are between fmin and fmax;
  // if not then abort. 
  // Nov 2013: check for nan and rename function checkval -> checkval_F

  int i;
  double val ;
  char fnam[] = "checkval_D" ;

  // ----------------- BEGIN -----------------

  for ( i = 0; i < nval;  i++ ) {
    val = *(dptr+i);
    
    if ( isnan(val) ) {
      sprintf(c1err,"%s(item %d) = nan", varname, i );
      sprintf(c2err,"Expected  %f < %s < %f", dmin, varname, dmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( val < dmin || val > dmax  ) {
      sprintf(c1err,"%s(item %d) = %f fails bounds-check", varname, i, val );
      sprintf(c2err,"Expected  %f < %s < %f", dmin, varname, dmax );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

}  // end of checkval_D

// ********************************
int landolt_ini__(int *opt, float *mag, float *kshift) {  
  // fortran wrapper
  int istat;
  istat = Landolt_ini(*opt, mag, kshift);
  return istat;
}

// ********************************
int Landolt_ini(
		int opt                // (I) BD17 or Vega
		,float *mag            // (I) primary mag  to store
		,float *kshift         // (I) shift k's for systematics
		) {

  // Read Landolt transformation parameters from file
  // and store primary mags in global var.
  // Landolt UBVRI,BX <=> 012345
  //
  // Jan 2, 2009: allow k_i values up to 1.0 instead of 0.1
  //              to allow for filter-adjustment tests that
  //              have larger color-transformations

  int ifilt, k ;

  float kval, kerr,  magtmp ;

  char 
    fnam[] = "Landolt_ini" 
    ,c_get[40]
    ,c_tmp[60]
    ,c_k[6]
    ,kfile[40]
    ,kfile_full[120]
    ;

  FILE *fp;

  // --------- BEGIN -------------


  print_banner("INIT  BESSELL <=> LANDOLT  TRANSFORMATIONS" );


  // init color terms to crazy value.

  for ( k=0; k <= 4 ; k++ ) {
    LANDOLT_COLOR_VALUE[k] = 0.0 ;
    LANDOLT_COLOR_ERROR[k] = 0.0 ;
  }

  // store mag in global array

  printf("   UBVRI,BX offsets: ");
  for ( ifilt=0; ifilt < NFILT_LANDOLT; ifilt++ ) {
    magtmp = *(mag + ifilt) ;
    LANDOLT_MAGPRIMARY[ifilt] = (double)magtmp ;
    printf(" %7.3f", magtmp );
  }
  printf("\n\n");

  if ( opt == 0 ) 
    goto PRINT_COLOR_TERMS ;
  else if ( opt < 4 ) 
    sprintf(kfile, "LANDOLT_COLOR_TERMS_BD17.DAT" );
  else
    sprintf(kfile, "LANDOLT_COLOR_TERMS_VEGA.DAT" );


  // read color terms from file
  // First try opening file in user dir ;
  // if not there than open official file in $SNDATA_ROOT.

  sprintf(c_tmp,"Ready to read Landolt color terms from");

  if ( (fp = fopen(kfile, "rt")) != NULL ) {
    printf("   %s user k-file: \n   %s \n", c_tmp, kfile );
    goto READ_KFILE ;
  }

  sprintf( PATH_SNDATA_ROOT, "%s", getenv("SNDATA_ROOT") );
  sprintf(kfile_full,"%s/filters/Landolt/%s", 
	  PATH_SNDATA_ROOT, kfile );

  if ( (fp = fopen(kfile_full, "rt"))==NULL ) {
    sprintf(c1err,"Cannot open: %s", kfile_full );
    errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
    fclose(fp); 
  }

  printf("   %s default k-file: \n   %s \n", c_tmp, kfile_full );

 READ_KFILE:

  while( fscanf(fp,"%s", c_get )!= EOF ) {

    for ( k=0; k <= 4; k++ ) {
      sprintf(c_k, "k%d:", k );
      if ( strcmp(c_get,c_k) == 0 ) {

	readfloat(fp, 1, &kval );
	LANDOLT_COLOR_VALUE[k] = (double)kval + *(kshift+k);

	readchar(fp, c_tmp);  // skip over '+-' symbol

	readfloat(fp, 1, &kerr );
	LANDOLT_COLOR_ERROR[k] = (double)kerr ;

      }
    }
  } // end of fscanf loop


  fclose(fp);


  // print color terms k0 - k4

 PRINT_COLOR_TERMS:

  for ( k=0; k <=4 ; k++ ) {
    kval = LANDOLT_COLOR_VALUE[k] ;
    kerr = LANDOLT_COLOR_ERROR[k] ;

    printf("\t   k%d = %6.3f +- %6.3f \n", k, kval, kerr ) ;

    if ( fabsf(kval) > 1.0 ) {
      sprintf(c1err,"k%d has invalid value above", k);
      errmsg(SEV_FATAL, 0, fnam, c1err, BLANK_STRING );
    }
  }

  printf("\t   U(reported) -> [UX - BX + B]_(synthetic) \n");


  return SUCCESS;

}  // end 




/**********************************************
  SALT-II color correction formula
**********************************************/


// ********************************
// fortran wrapper
int landolt_convert__(int *opt, double *mag_in, double *mag_out) {
  int istat ;
  istat = Landolt_convert(*opt, mag_in, mag_out);
  return istat ;
}

// ********************************
int Landolt_convert(int opt, double *mag_in, double *mag_out) {

  /*****
   Created Mar 11, 2007 by R.Kessler

   opt = +1 : *mag_in  = synthetic Bessell UBVRI,BX
              *mag_out = reported  Landolt UBVRI 

   opt = -1 : *mag_in  = reported  Landolt UBVRI, BX-B
              *mag_out = synthetic Bessell UBVRI, BX

         (but note that reported U is (UX - BX + B)_synth


  ******/

  int ifilt;
  double k0, k1, k2, k3, k4 ;

  int off_U  = 0 ;
  int off_B  = 1 ;
  int off_V  = 2 ;
  int off_R  = 3 ;
  int off_I  = 4 ;
  int off_BX = 5 ;

  double del, delref;
  double DEL_V, DEL_BV, DEL_UBX, DEL_VR, DEL_RI;
  double Vout, Utmp, DEL_BXB ;

  // ------------ BEGIN ----------------

  // init *mag_out
  for ( ifilt=0; ifilt < NFILT_LANDOLT ; ifilt++ ) {
    *(mag_out+ifilt) = -99.0 ; 
  }


  k0 = LANDOLT_COLOR_VALUE[0];
  k1 = LANDOLT_COLOR_VALUE[1];
  k2 = LANDOLT_COLOR_VALUE[2];
  k3 = LANDOLT_COLOR_VALUE[3];
  k4 = LANDOLT_COLOR_VALUE[4];

  // apply magdif array

  if ( opt > 0 ) {  // convert Bessell -> Landolt


    del    = *(mag_in+off_B)       - *(mag_in+off_V);
    delref = LANDOLT_MAGPRIMARY[off_B] - LANDOLT_MAGPRIMARY[off_V] ;
    DEL_V  = k0*(del-delref);
    DEL_BV = k1 * (del-delref);

    del     = *(mag_in+off_U)       - *(mag_in+off_BX);
    delref  = LANDOLT_MAGPRIMARY[off_U] - LANDOLT_MAGPRIMARY[off_BX] ;
    DEL_UBX = k2 * (del-delref);

    del     = *(mag_in+off_V)       - *(mag_in+off_R);
    delref  = LANDOLT_MAGPRIMARY[off_V] - LANDOLT_MAGPRIMARY[off_R] ;
    DEL_VR  = k3 * (del-delref);

    del     = *(mag_in+off_R)       - *(mag_in+off_I);
    delref  = LANDOLT_MAGPRIMARY[off_R] - LANDOLT_MAGPRIMARY[off_I] ;
    DEL_RI  = k4 * (del-delref);

    *(mag_out + off_V) = *(mag_in + off_V) + DEL_V ;
    *(mag_out + off_B) = *(mag_in + off_B) + DEL_V + DEL_BV ;
    *(mag_out + off_R) = *(mag_in + off_R) + DEL_V - DEL_VR ;
    *(mag_out + off_I) = *(mag_in + off_I) + DEL_V - DEL_VR - DEL_RI ;

    *(mag_out + off_U) = *(mag_in+off_U) - *(mag_in+off_BX) + *(mag_in+off_B)
     + DEL_V + DEL_BV + DEL_UBX ;

  } 

  else if ( opt < 0 ) {  // convert Landolt -> Bessell

    del    = *(mag_in+off_B)       - *(mag_in+off_V) ;
    delref = LANDOLT_MAGPRIMARY[off_B] - LANDOLT_MAGPRIMARY[off_V] ;
    DEL_BV = (del + k1*delref) / ( 1. + k1 ) ;  // (B-V)_Bess
    DEL_V  = k0 * (delref - DEL_BV );  // V_Bess - V_Land

    del    = *(mag_in+off_V)       - *(mag_in+off_R) ;
    delref = LANDOLT_MAGPRIMARY[off_V] - LANDOLT_MAGPRIMARY[off_R] ;
    DEL_VR = (del + k3*delref) / ( 1. + k3 ) ;  // (V-R)_Bess

    del    = *(mag_in+off_R)       - *(mag_in+off_I) ;
    delref = LANDOLT_MAGPRIMARY[off_R] - LANDOLT_MAGPRIMARY[off_I] ;
    DEL_RI = (del + k4*delref) / ( 1. + k4 ) ;  // (R-I)_Bess

    del    = *(mag_in+off_U)       - *(mag_in+off_B) ;
    delref = LANDOLT_MAGPRIMARY[off_U] - LANDOLT_MAGPRIMARY[off_B] ;
    DEL_UBX = (del + k2*delref) / ( 1. + k2 );   // (UX-BX)_Bess


    Vout = *(mag_in  + off_V) + DEL_V ;
    *(mag_out+off_V) = Vout ;
    *(mag_out+off_B) = Vout + DEL_BV ;
    *(mag_out+off_R) = Vout - DEL_VR ;
    *(mag_out+off_I) = *(mag_out+off_R) - DEL_RI ;

    Utmp = DEL_UBX + *(mag_out+off_B) ;  // reported U = UX-BX+B

    // to get synthetic U, add synthetic BX-B

    DEL_BXB = *(mag_in+off_BX) ;
    *(mag_out+off_U) = Utmp + DEL_BXB ;

    // BX = (BX-B)_in + B_out
    *(mag_out+off_BX) = DEL_BXB + *(mag_out+off_B);

  }

  return SUCCESS ; // add Aug 7 2014 to avoid compile warning.

}  // end of function Landolt_convert



// *************************************
float edgedist ( float X, float Y, int NXPIX, int NYPIX ) {

  float XX, YY, DX, DY, DMIN;

  XX = (float)NXPIX - X;
  YY = (float)NYPIX - Y;      

  // take min by brute force since fminf does not work ???
  if ( X < XX ) DX = X; else DX = XX;
  if ( Y < YY ) DY = Y; else DY = YY;

  if ( DX < DY ) DMIN = DX; else DMIN = DY;

  return DMIN ;

} // end of edgedist


// ==================================
void find_pathfile(char *fileName, char *PATH_LIST, char *FILENAME, char *funCall){

  // Created Sep 28 2020
  // For input *fileName, check if it exists in CWD, or in any
  // path in (space-separated) PATH_LIST.
  // Return FILENAME that includes full path (if needed).
  // Abort if file not found.

  // ------------ BEGIN ----------------

#define MXPATH_CHECK 4

  struct stat statbuf, statbuf_gz;
  bool FOUNDIT = false;
  int  jstat, jstat_gz, ipath, NPATH ;
  char *path, *PATH[MXPATH_CHECK], sepKey[] = " ";
  char  tmpName[MXPATHLEN], tmpName_gz[MXPATHLEN] ;
  char fnam[] = "find_pathfile";

  // ------------- BEGIN -----------

  for(ipath=0; ipath < MXPATH_CHECK; ipath++ )
    { PATH[ipath] = (char*) malloc(MXPATHLEN*sizeof(char) ); }

  splitString(PATH_LIST, sepKey, MXPATH_CHECK,
	       &NPATH, &PATH[1] ); // <== returned

  NPATH++; sprintf(PATH[0],"");

  for ( ipath=0; ipath < NPATH; ipath++ ) {
    path = PATH[ipath] ;
    if ( strlen(path) == 0 ) 
      { sprintf(tmpName,"%s", fileName); }
    else
      { sprintf(tmpName,"%s/%s", path, fileName); }

    sprintf(tmpName_gz, "%s.gz", tmpName);
    jstat    = stat(tmpName,    &statbuf);    // returns 0 if file exists
    jstat_gz = stat(tmpName_gz, &statbuf_gz); // returns 0 if file exists

    if ( jstat==0 || jstat_gz==0 ) 
      { FOUNDIT = true ;   sprintf(FILENAME, "%s", tmpName);  break; }

  } // end ipath
  
  // - - - - - - - - - - 
  if ( !FOUNDIT ) {
    print_preAbort_banner(fnam);
    printf("   Called by function %s \n", funCall);
    printf("   jstat=%d jstat_gz=%d \n", jstat, jstat_gz);
    for ( ipath=1; ipath < NPATH; ipath++ ) 
      { printf("   File not in path: %s \n", PATH[ipath] ); }
    sprintf(c1err,"Could not find file in paths listed above:");
    sprintf(c2err,"fileName = '%s'", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);   
  }

  //  debugexit(fnam);

  for(ipath=0; ipath < MXPATH_CHECK; ipath++ )  { free(PATH[ipath]); }


  return;
} // end find_pathfile

// ==================================
FILE *open_TEXTgz(char *FILENAME, const char *mode, int *GZIPFLAG ) {

  // Dec 1 2017:
  // Shell to call fopen.
  // If mode is read, then check for gzipped file.
  // If both gzip and unzip file exist, wait 5 sec and
  // try again ... if both still exist then ABORT.
  // Return GZIPFLAG=1 if gzip file is opened
  //
  // Mar 6 2019: replace zcat with 'gunzip -c' so that it works on Mac.
  // Dec 17 2021: add logic to wait 5 sec if file and file.gz both exist.
  //
  FILE *fp ;
  struct stat statbuf ;
  int istat_gzip, istat_unzip, LEN, N_ITER=0;
  char gzipFile[MXPATHLEN], unzipFile[MXPATHLEN];
  char cmd_zcat[MXPATHLEN] ;
  char fnam[]=  "open_TEXTgz" ;

  // -------------- BEGIN ------------

 START:

  N_ITER++ ;
  *GZIPFLAG = 0 ;

  // for reading, check for gz file
  if ( strstr(mode,"r") != NULL ) {

    if ( strstr(FILENAME,".gz") != NULL )  { 
      // .gz file name passed as argument
      LEN = strlen(FILENAME);
      sprintf(gzipFile,   "%s", FILENAME); 
      sprintf(unzipFile,  "%s", FILENAME);
      unzipFile[LEN-3] = 0 ; // remove .gz extension
    } 
    else {
      // check for .gz file if not given
      sprintf(gzipFile,  "%s.gz", FILENAME);
      sprintf(unzipFile, "%s",    FILENAME);
    } 

    istat_gzip  = stat(gzipFile,  &statbuf) ;
    istat_unzip = stat(unzipFile, &statbuf) ;

    //    printf(" xxx ------------------------------- \n");
    //    printf(" xxx istat=%3d for '%s' \n", istat_gzip,  gzipFile);
    //    printf(" xxx istat=%3d for '%s' \n", istat_unzip, unzipFile);

    bool FOUND_2FILES = ( istat_gzip==0 && istat_unzip==0 );
    if ( FOUND_2FILES && N_ITER==1 ) { 
      printf("\t found gzip and unzip file ... try again in 5 sec... \n");
      fflush(stdout);
      sleep(5.0); 
      goto START;    
    }

    if ( FOUND_2FILES ) {
      print_preAbort_banner(fnam);
      printf("  Found %s \n", gzipFile );
      printf("  Found %s \n", unzipFile );
      sprintf(c1err, "Found gzipped and unzipped file."); 
      sprintf(c2err, "Cannot open both files.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    if ( istat_gzip == 0 ) {
      sprintf(cmd_zcat, "gunzip -c %s", gzipFile);
      fp = popen(cmd_zcat,"r");
      *GZIPFLAG = 1 ;
      return(fp);
    }
  }

  // if we get here, do regular text open
  fp = fopen(FILENAME,mode);
  return(fp);

} // end open_TEXTgz


// =====================================
void snana_rewind(FILE *fp, char *FILENAME, int GZIPFLAG) {

  int gzipFlag ;
  char fnam[] = "snana_rewind" ;

  // --------------- BEGIN ----------------

  if ( GZIPFLAG == 0 ) {
    rewind(fp);
  }
  else {
    // cannot rewind popen for gz file; must close and re-open
    int istat = pclose(fp);
    if ( istat == -1 ) {
      sprintf(c1err,"pclose Error = -1 for file=");
      sprintf(c2err,"%s", FILENAME);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    }
    fp = open_TEXTgz(FILENAME, "rt", &gzipFlag);
  }

} // end snana_rewind


// *************************************************
FILE *snana_openTextFile (int OPTMASK, char *PATH_LIST, char *fileName, 
			  char *fullName, int *gzipFlag ) {

  /* ----------------------------------------------
    Shell to open text file for reading.
    First search local directory; if file not found
    then search paths in [space-separated] PATH_LIST.
    The search priority is
      + current pwd
      + 1st path in PATH_LIST
      + 2nd path in PATH_LIST
    There is no warning or error if file exists in 
    multiple directories.

    This allows user to easily over-ride default
    with private version.

   Inputs
     + OPTMASK :
         += 1 -> verbose mode
         += 2 -> abort if DOCUMENTATION is not first string in file
         += 4 -> don't check for DOCUMENTATION

     + PATH_LIST   : optional space-separated list of paths to check
     + fileName    : name of file to option
   Outputs
     + fullName : full name of file, including path
     + gzipFlag : 1 if gzipped, 0 otherwise

   Function returns file pointer.

   Dec 29 2017: use open_TEXTgz to allow reading gzipped files.
   Jan 11 2018: add gzipFile output arg
   Mar 20 2019: padd vboseFlag to print comment to stdout
   Feb 01 2020: SNPATH -> PATH_LIST (space separated)
   Aug 26 2020: change vboseflag to more generic OPTMASK
   May 20 2021: OPTMASK +=4 -> don't check for DOCANA

  ----------------------------------------------- */

#define TEXTMODE_read "rt"
#define MXPATH_CHECK 4

  bool VBOSE          = ( (OPTMASK & OPENMASK_VERBOSE)        > 0 ) ;
  bool REQUIRE_DOCANA = ( (OPTMASK & OPENMASK_REQUIRE_DOCANA) > 0 ) ;
  bool IGNORE_DOCANA  = ( (OPTMASK & OPENMASK_IGNORE_DOCANA)  > 0 ) ;
  int ipath, NPATH ;
  bool IS_OPEN = false ;
  char *PATH[MXPATH_CHECK], sepKey[]= " " ; 
  FILE *fp ;
  char fnam[] = "snana_openTextFile" ;

  // --------------- BEGIN ----------------

  // First try current working directory

  sprintf(fullName, "%s", fileName );

  fp = open_TEXTgz(fullName,TEXTMODE_read, gzipFlag );
  if ( fp != NULL ) {       
    if ( VBOSE )  { printf("\t Opened : %s \n", fullName ); }
    goto DONE ;
  } 

  // if we get here, try paths in PATH_LIST

  for(ipath=0; ipath < MXPATH_CHECK; ipath++ )
    { PATH[ipath] = (char*) malloc(MXPATHLEN*sizeof(char) ); }

  splitString(PATH_LIST, sepKey, MXPATH_CHECK,
	       &NPATH, PATH ); // <== returned

  for(ipath=0; ipath < NPATH; ipath++ ) {
    if ( IS_OPEN ) { continue ; }

    sprintf(fullName, "%s/%s", PATH[ipath],  fileName );
    fp = open_TEXTgz(fullName,TEXTMODE_read, gzipFlag );

    if ( fp != NULL ) {
      IS_OPEN = true ;
      if ( VBOSE ) { printf("\t Opened : %s \n", fullName ); }
    }

  } // end ipath

  // free memory
  for(ipath=0; ipath < MXPATH_CHECK; ipath++ )   { free(PATH[ipath]); }

 DONE:

  if  ( fp != NULL  && !IGNORE_DOCANA ) { 
    bool FOUND_DOCANA ;
    FOUND_DOCANA = check_openFile_docana(REQUIRE_DOCANA,fp,fullName); 

    // if we didn't find DOCANA key (and didn't abort), then need to
    // rewind in case that first alread-read key is needed later.
    // Note that rewind doesn't work for gzip files (?)
    if ( !FOUND_DOCANA && *gzipFlag == 0 ) 
      { snana_rewind(fp, fullName, *gzipFlag); }
  }

  // return pointer regardless of status
  return fp ;

}  // end of snana_openTextFile

// ******************************************************
bool key_in_file(char *fileName, char *key, int nwd_check) {

  // Created Jun 2021
  // Return true of *key exists in *fileName, 
  // and it is among fist nwd_check words.

  int gzipFlag, nwd_read=0;
  bool found_key = false ;
  char c_get[MXPATHLEN];
  FILE *fp;
  // ------------- BEGIN -----------

  fp  = open_TEXTgz(fileName, "rt", &gzipFlag);

  while( fscanf(fp,"%s", c_get )!= EOF ) {
    if ( strcmp(c_get,key) == 0 ) { found_key=true; break; }
    nwd_read++ ;
    if ( nwd_read > nwd_check ) { break; }
  }

  if ( gzipFlag ) { pclose(fp); }    else  { fclose(fp); }

  return found_key ;

} // end key_in_file

// ===============================================
int colnum_in_table(char *fileName, char *varName) {

  // Created Jun 25 2021
  // if "VARNAMES:" key exists among first nwd_check words,
  // search for *varname and return column number colnum, 
  // where colnum=0 for CID.
  //
  // ERROR codes:
  //   colnum = -1 => file does not exist
  //   colnun = -2 => VARNAMES key does not exist
  //   colnum = -3 => VARNAMES key exists, but *varname not found.
  //

  int colnum = -1;
  int nwd_check = 1000; // stop checking after this many words
  int gzipFlag, nvar, ivar, nwd_read=0;
  char KEY[] = "VARNAMES:" ;
  char c_get[MXPATHLEN], VARNAME_STRING[MXPATHLEN];
  FILE *fp;

  char fnam[] = "colnum_in_file" ;

  // ------------ BEGIN --------------

  fp  = open_TEXTgz(fileName, "rt", &gzipFlag);

  colnum = -1;
  if ( !fp ) { return(colnum); }

  colnum = -2;
  while( fscanf(fp,"%s", c_get )!= EOF ) {

    if ( strcmp(c_get,KEY) == 0 ) { 
      colnum = -3;
      fgets(VARNAME_STRING, MXPATHLEN, fp);
      nvar = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING,VARNAME_STRING);
      for(ivar=0; ivar < nvar; ivar++ ) {
	get_PARSE_WORD(0, ivar, c_get);
	if ( strcmp(varName,c_get) == 0 ) { colnum = ivar; }
      }
      break; 
    }

    nwd_read++ ;
    if ( nwd_read > nwd_check ) { break; }
  }

  if ( gzipFlag ) { pclose(fp); }  else  { fclose(fp); }

  return colnum ;

} // end colnum_in_file



// *************************************
bool check_openFile_docana(bool REQUIRE_DOCANA, FILE *fp, char *fileName) {
  // Created Aug 26 2020
  // For already open file (fp), abort if first word in file is not 
  // a DOCANA key.
  // A separate abort function is available to both C and fortran.
  //
  // Return FOUND_DOCANA logical
  //
  // Nov 18 2021; use fgets instead of fscanf
  // Nov 23 2021: back to simpler fscanf 
  
  char line[MXPATHLEN], key[60], *ptrtok, *pos;
  bool FOUND_DOCANA ;
  char fnam[] = "check_openFile_docana";
  // ------------- BEGIN --------

  fscanf(fp,"%s",key); 
  FOUND_DOCANA = (strcmp(key,KEYNAME_DOCANA_REQUIRED) == 0 );

  if ( !FOUND_DOCANA ) 
    { react_missing_docana(REQUIRE_DOCANA,fileName); }

  return FOUND_DOCANA ;

} // end check_openFile_docana

// *************************************
void check_file_docana(int optmask, char *fileName) {
  // Created Aug 26 2020
  // Open and read first line of fileName; if not DOCANA key 
  // then abort or give warning based on input optmask.
  // A separate abort function is available to both C and fortran.
  //
  // Input 
  //    optmask & 1 -> abort if no DOCANA; else give warning.

  int MSKOPT   = MSKOPT_PARSE_WORDS_FILE + MSKOPT_PARSE_WORDS_FIRSTLINE ;
  bool REQUIRE = ( (optmask & 1) > 0 ) ;
  bool FOUND_DOCANA = false ;
  int  langFlag=0, iwd, NWD;
  char key[60];
  char fnam[] = "check_file_docana";
  // ------------- BEGIN --------

  NWD = store_PARSE_WORDS(MSKOPT, fileName);

  for(iwd=0; iwd < NWD; iwd++ ) {
    get_PARSE_WORD(langFlag, iwd, key);
    if ( strcmp(key,KEYNAME_DOCANA_REQUIRED)==0 ) 
      { FOUND_DOCANA = true; }
  }

  if ( !FOUND_DOCANA )  { react_missing_docana(REQUIRE,fileName); }

  return;
} // end check_file_docana

void abort_bad_input(char *key, char *word, int iArg, char *callFun) {
  sprintf(c1err,"Could not read arg[%d] = '%s'", iArg, word);
  sprintf(c2err,"for KEY = %s (check N_arg after input key)", key );
  errmsg(SEV_FATAL, 0, callFun, c1err, c2err ) ;
} // end abort_bad_input

void react_missing_docana(bool REQUIRE_DOCANA, char *fileName) {
  // Calling function has already determined that DOCANA is missing;
  // here, either abort or give warning based on REQUIRE_DOCANA.
  // 
  bool GIVE_WARNING = !REQUIRE_DOCANA ;
  char fnam[] = "react_missing_docana" ;

  if ( REQUIRE_DOCANA ) { print_preAbort_banner(fnam); }

  // always print this stuff on missing DOCANA
  printf("\n");

  if ( GIVE_WARNING ) 
    { printf("\t ***** DOCANA WARNING ******* \n"); }
  printf("  Missing required '%s'  key in \n", KEYNAME_DOCANA_REQUIRED);
  printf("    %s \n", fileName);
  printf("  See DOCANA examples with linux command: \n");
  printf("    grep -R DOCUMENTATION_END $SNDATA_ROOT \n") ;
  printf("  File must begin with 'DOCUMENTATION:' key\n");
  printf("\n");

  if( REQUIRE_DOCANA ) {
    sprintf(c1err,"See DOCANA error above. Must add DOCUMENTATION block to");
    sprintf(c2err,"%s", fileName );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ) ;
  }

}   // end react_missing_docana

void react_missing_docana__(bool *REQUIRE_DOCANA, char *fileName) 
{ react_missing_docana(*REQUIRE_DOCANA,fileName); }
void check_file_docana__(int *optmask, char *fileName) 
{ check_file_docana(*optmask, fileName); }


// **********************************
void abort_docana_tooLong(char *fileName, char *callFun) {
  sprintf(c1err,"Could not find %s",  KEYNAME2_DOCANA_REQUIRED);
  sprintf(c2err,"in %s", fileName );
  errmsg(SEV_FATAL, 0, callFun, c1err, c2err ) ;
} // end abort_docana_tooLong

// *****************************************************
void abort_openTextFile(char *keyName, char *PATH_LIST,
                        char *fileName, char *funCall) {

  // if *snana_openTextFile returns NULL pointer, call this
  // function to print error info and abort. Main thing is
  // to print out all directories that were checked from
  // current dir and PATH_LIST.
  
  //#define MXPATH_CHECK 4
  int NPATH, ipath;
  char *PATH[MXPATH_CHECK], sepKey[] = " " ;

  // -------------- BEGIN ----------------
  print_preAbort_banner(funCall);

  // append path(s) from PATH_LIST
  for(ipath=0; ipath < MXPATH_CHECK; ipath++ )
    { PATH[ipath] = (char*) malloc( MXPATHLEN * sizeof(char) ); }

  // first path is current working dir
  getcwd(PATH[0],MXPATHLEN);

  // split string to get other paths
  splitString(PATH_LIST, sepKey, MXPATH_CHECK,
	      &NPATH, &PATH[1] ); // <== returned
  
  printf("\n  The following paths were checked for input file read from"
	 "\n  %s: %s \n", keyName, fileName);
  for(ipath=0; ipath <= NPATH; ipath++ ) 
    { printf("     %s \n", PATH[ipath] );   }

  sprintf(c1err,"Could not find '%s' input file:", keyName );
  sprintf(c2err,"%s  (see preAbort info above)", fileName);
  errmsg(SEV_FATAL, 0, funCall, c1err, c2err );
  
  return ;

} // end abort_openTextFile

// ***************************
int INTFILTER ( char *cfilt ) {

  // returns absolute filter index  for string *cfilt
  // Oct 29 2019: use last char of cfilt to work with arbitrary string

  int len = strlen(cfilt);
  int ifilt, itmp;
  char ctmp[2], cfilt1[2];
  //---------- BEGIN ----------------

  sprintf(cfilt1, "%c", cfilt[len-1]);
  ifilt = 0 ;
  for ( itmp=0; itmp < MXFILTINDX; itmp++ ) {
    sprintf(ctmp,"%c", FILTERSTRING[itmp] );
    if (strcmp(ctmp,cfilt1) == 0 ) { ifilt = itmp; return ifilt ; }
  }

  return ifilt;
} // end of INTFILTER



// *****************************************************
int wr_filtband_int ( 
		     FILE *fp        // write to this  *fp
		     ,char *keyword  //keyword before integers
		     ,int NINT       // number of integers to write
		     ,int *iptr      // pointer to first integer
		     ,char *comment  // comment after integers
		     ,int opt        // 0=> header format; 1=> BAND format
		     ) {

  int i;

  // ------------- BEGIN -----------

  if ( opt == 1 ) {
    fprintf(fp,"  %16s  ", keyword);

    for ( i = 0; i < NINT ; i++ ) 
      fprintf(fp, "%8d ", *(iptr+i) );
  }
  else {
    fprintf(fp,"%s  ", keyword);

    for ( i = 0; i < NINT ; i++ ) 
      fprintf(fp, "%4d ", *(iptr+i) );
  }

  fprintf(fp," %s \n" , comment ) ;

  return(0); // added Aug 7 2014 to avoid compile warning

} // end of wr_filtband_int


// *****************************************************
int wr_filtband_float(
		      FILE *fp 
		      ,char *keyword  // keyword
		      ,int NFLOAT     // number of floats to write
		      ,float *fptr    // point to 1st float
		      ,char *comment  // comment after floats
		      ,int idec       // number of digits after decimal
     ) {

  int i;

  // ---------------- BEGIN --------------


  if ( idec > 0 ) 
    fprintf(fp,"  %16s  ", keyword  ) ;  // for epoch
  else
    fprintf(fp,"%20s", keyword  ) ;  // for header

  for ( i=0; i<NFLOAT; i++ ) {
    if ( idec == 2 ) 
      fprintf(fp,"%8.2f ", *(fptr+i) );
    else if ( idec == 3 ) 
      fprintf(fp,"%8.3f ", *(fptr+i) );
    else if ( idec == 4 ) 
      fprintf(fp,"%8.4f ", *(fptr+i) );
    else if ( idec == 5 ) 
      fprintf(fp,"%8.5f ", *(fptr+i) );

    else if ( idec < 0  ) 
      fprintf(fp,"%6.1f ", *(fptr+i) );
  }

  fprintf(fp," %s \n" , comment ) ;

  return(0); // add Aug 7 2014 to avoid compile warnings.

} // end of wr_filtband_float

// *********************************************************
void check_EOF(FILE *fp, char *file_name, char *fun_name, int nline_read) {

  // Created Dec 3 2021 by R.Kessler
  // if feof(fp) is True, return with no comment.
  // if !feof(fp), abort on error that EOF was not reached.
  // Initial motivation is to better identify file-read
  // errors on cori.

  int ieof = feof(fp);
  int ierr = ferror(fp);
  char fnam[] = "check_EOF" ;

  // ---------- BEGIN -----------

  if ( ierr == 0 ) { return; }

  print_preAbort_banner(fnam);
  printf("   EOF failed for: \n");
  printf("    file name: %s\n", file_name);
  printf("    from function: %s\n", fun_name);
  printf("    after reading %d lines. \n", nline_read);

  printf("    feof(fp)   = %d \n", ieof );
  printf("    ferror(fp) = %d \n", ierr );

  sprintf(c1err,"File reading ended without reaching EOF;");
  sprintf(c2err,"Check file-system or if file was altered during read.");
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 

  return;

} // end check_EOF


// ***********************************
void check_argv(void) {

  // ABORT if there are un-used command-line args
  // (to catch mis-typed commands)

  int NBAD, i ;
  char fnam[] = "check_argv";

  // ----------- BEGIN ---------

  NBAD = 0;
  for ( i = 1 ; i < NARGV_LIST; i++ ) {
    if ( USE_ARGV_LIST[i] == 0 ) {
      NBAD++ ;
      if ( NBAD == 1 ) printf("  CHECK_ARGV ERRORS: \n" );

      printf(" \t ERROR: detected un-used command-line arg: %s \n", 
	     ARGV_LIST[i]);
    }
  }

  if ( NBAD > 0 ) {
    sprintf(c1err,"%d invalid/unknown command line arg(s)", NBAD);
    sprintf(c2err,"Check all %d command-line args", NARGV_LIST );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

} // end check_argv

void check_arg_len(char *keyName, char *arg, int MXLEN) {
  char fnam[] = "check_arg_len";
  int LEN = strlen(arg);
  if ( LEN >= MXLEN ) {
    sprintf(c1err,"len(%s) = %d exceeds bound of %d", keyName, LEN, MXLEN);
    sprintf(c2err,"Shorten arg, or increase array bound.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
} // end check_argLen
 
// *******************************************************
void parse_err ( char *inFile, int NEWMJD, char *keyword ) {
  // print standard error message for parsing.
  // errmsg called with ABORT flag !
  char fnam[] = "parse_err" ;
  sprintf(c1err,"Problem parsing %s", inFile );
  sprintf(c2err,"Check NEWMJD-Epoch %d, keyword '%s' ", NEWMJD, keyword);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
}

// *********************************************************
void  legacyKey_abort(char *callFun, char *legacyKey, char *newKey) {

  if ( strlen(newKey) > 0 ) {
    sprintf(c1err,"Legacy key '%s' is no longer valid.", legacyKey);
    sprintf(c2err,"Use %s key instead .", newKey );
    errmsg(SEV_FATAL, 0, callFun, c1err, c2err );
  }
  else {
    sprintf(c1err,"Legacy key '%s' is no longer valid.", legacyKey);
    sprintf(c2err,"Remove key from file ."  );
    errmsg(SEV_FATAL, 0, callFun, c1err, c2err );
  }

} // end legacyKey_abort

// *********************************************************
void missingKey_ABORT(char *key, char *file, char *callFun) {
  // created Jan 2014
  char fnam[] = "missingKey_ABORT";
  sprintf(c1err,"%s could not find KEY='%s'", callFun, key);
  sprintf(c2err,"in file = %s", file);
  errmsg(SEV_FATAL, 0, fnam, c1err, c2err );

} // end of missingKey_ABORT

void tabs_ABORT(int NTAB, char *fileName, char *callFun) {
  char fnam[] = "tabs_ABORT";
  if ( NTAB > 0 ) {
    sprintf(c1err,"%s found %d invalid tabs in file ", 
	    callFun, NTAB);
    sprintf(c2err,"%s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }
} // end tabs_ABORT

// ************************************************
void errmsg(
       int isev            /* (I) severity flag */
      ,int iprompt         /* (I) 1=>prompt user to continue */
      ,char *fnam          /* (I) name of function calling errmsg */
      ,char *msg1          /* (I) message to print           */
      ,char *msg2          /* (I) 2nd msg to print ("" => no 2nd msg) */
          ) 
/*************************************
   Created Jan 30, 2004  by R.S. Kessler.

   ABSTRACT:
      Print error message(s), and count no. of errors.
      Abort program if isevere == SEV_ABORT
        
   Dec 17 2019: re-structure a bit to prep for fortran use

*************************************/
{
  char c_severe[40];  // char string for serverity 
        /* ------- begin function execution ------- */
   switch ( isev )
     {
        case SEV_INFO  : sprintf(c_severe,"INFO");     break;
        case SEV_WARN  : sprintf(c_severe,"WARNING");  break;
        case SEV_ERROR : sprintf(c_severe,"ERROR");    break;
        case SEV_FATAL : sprintf(c_severe,"FATAL ERROR ABORT");  break;
     }
   if ( isev == SEV_FATAL ) {  madend(stdout,0); }
   // print severity and name of function
   printf(" %s called by %s\n", c_severe, fnam);
   // print each message
   if( strlen(msg1) > 0 ) { printf("   %s\n", msg1 ); }
   if( strlen(msg2) > 0 ) { printf("   %s\n", msg2 ); }
   printf("\n");   fflush(stdout);
   if( isev == SEV_FATAL ) { exit(EXIT_ERRCODE); }
   return ;

}   /* end of function "errmsg" */

void errmsg_( int *isev,int *iprompt, char *fnam, char *msg1, char *msg2 ) 
{  errmsg(*isev, *iprompt, fnam, msg1, msg2); }


// ************************************************
void errlog(
	    FILE *fp      // FP or stdout
	    ,int isev     // (I) severity flag
	    ,char *fnam   // (I) name of function calling errmsg
	    ,char *msg1   // (I) message to print 
	    ,char *msg2   // (I) 2nd msg to print
          ) 
/*************************************
   Created Jun 2021

   Same ass errmsg, except:
    + here fp is an argument to enable sending message to log file
    + remove obolete iprompt argument.

*************************************/
{
  char c_severe[40];  // char string for serverity 
  // ------- begin function execution ------- 
   switch ( isev )
     {
        case SEV_INFO  : sprintf(c_severe,"INFO");     break;
        case SEV_WARN  : sprintf(c_severe,"WARNING");  break;
        case SEV_ERROR : sprintf(c_severe,"ERROR");    break;
        case SEV_FATAL : sprintf(c_severe,"FATAL ERROR ABORT");  break;
     }
   if ( isev == SEV_FATAL ) {  madend(fp,0); }
   // print severity and name of function
   fprintf(fp, " %s called by %s\n", c_severe, fnam);
   // print each message
   if( strlen(msg1) > 0 ) { fprintf(fp,"   %s\n", msg1 ); }
   if( strlen(msg2) > 0 ) { fprintf(fp,"   %s\n", msg2 ); }
   fprintf(fp,"\n");   fflush(stdout);
   if( isev == SEV_FATAL ) { exit(EXIT_ERRCODE); }
   return ;

}   // end of function "errlog"

// ************************************************
void madend(FILE *fp, int flag) {

   char cmsg[40] = { "ABORT program on Fatal Error." };

   //   printf("\n");   printf("\n");
   fprintf(fp, "\n   `|```````|`    ");
   fprintf(fp, "\n   <| o\\ /o |>    ");
   fprintf(fp, "\n    | ' ; ' |     ");
   fprintf(fp, "\n    |  ___  |     %s ", cmsg);
   fprintf(fp, "\n    | |' '| |     ");
   fprintf(fp, "\n    | `---' |     ");
   fprintf(fp, "\n    \\_______/    ");
   fprintf(fp, "\n");

   fprintf(fp,"\n");   
   fflush(fp);

   if ( flag == 1 ) { exit(EXIT_ERRCODE); }

}    //  end of madend  


// ************************************************
void happyend(void) {
   fflush(stdout);
   printf("\n Program stopping gracefully. Bye. \n");
   exit(0);
}


void  print_preAbort_banner(char *fnam) {  
  printf("\n\n");
  printf("# ================================================ \n");
  printf("  PRE-ABORT DUMP from function %s : \n", fnam); 
  fflush(stdout);
}

void  print_preabort_banner__(char *fnam) 
{  print_preAbort_banner(fnam); }
 

// *************************************************
int read_genpoly(char *KEYNAME, char **WORDS, int order_legacy, 
		 GENPOLY_DEF *POLY) {

  // Created Sep 30 2020
  // Parse polynomial as either comma-sep or space sep; e.g,
  //
  //   1.01,0.3,0.04   # 2nd order, command sep
  //      or
  //   1.01 0.3 0.04   # same, but space sparated -> legacy
  //
  //   Comma sep can have any order, while space-sep is
  //   restricted to order_legacy passed as argument.
  //
  // Dec 21 2021: increase char arrays from len=60 to 200
  //
  int  MEMD, nread, j, N=0;
  double *zpoly_legacy ;
  char first_word[200], cpoly[200] ;
  char fnam[] = "read_genpoly" ;

  // -------------- BEGIN -----------

  /*
  printf(" xxx %s: WORDS = %s %s %s   order=%d \n",
	 fnam, WORDS[0], WORDS[1], WORDS[2], order_legacy );
  */

  sscanf(WORDS[N], "%s", first_word); N++ ;

  // if there is a comma, read comma-sep poly coefficients in one read
  if ( strstr(first_word,COMMA) ) {
    sprintf(cpoly, "%s", first_word);
  }
  else {
    // read zpoly_legacy as space-separated, and check that each
    // each coeff is indeed float and not a string
    cpoly[0] = 0;
    MEMD = (order_legacy+2)*sizeof(double);
    zpoly_legacy = (double*) malloc(MEMD);

    nread = sscanf(first_word, "%le", &zpoly_legacy[0]);
    catVarList_with_comma(cpoly,first_word); 

    for(j=1; j <= order_legacy; j++ ) { 
      nread = sscanf(WORDS[N], "%le", &zpoly_legacy[j] ); 
      if ( nread != 1 ) { abort_bad_input(KEYNAME,WORDS[N],j,fnam); }
      catVarList_with_comma(cpoly,WORDS[N]); 
      N++ ;
    }

    free(zpoly_legacy);
  }

  // parse POLY struct
  parse_GENPOLY(cpoly, KEYNAME, POLY, fnam);

  return N ;

} // end read_genpoly

// ************************************************
void readint(FILE *fp, int nint, int *list)   
/*****
   Created 7-feb-95 by R.Kessler

   ABSTRACT: 
     Read "nint" integers from configuration file, 
     and put in array "list." 

  Apr 2019: abort on reading string
*****/
{
    char c_get[80];
    int  i, itmp, istat, scanStat ;
    char fnam[] = "readint";

    for ( i=0; i<nint; i++) {
      scanStat = fscanf(fp,"%s",c_get);         // read next string 
      itmp     = NOINT  ;
      istat    = sscanf ( c_get, "%d", &itmp );

      if ( itmp == NOINT ) {
	sprintf(c1err,"Could not read int from string='%s' ; ", c_get);
	sprintf(c2err,"reading item %d of %d items.", i+1, nint);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      list[i] = itmp; // xxxx atoi(c_get);       // convert to int   
     }

}  // end of function "readint" 


// ************************************************
void readlong(FILE *fp, int nint, long long *list)   
/*****
   Created Feb 2015 by R.Kessler

   ABSTRACT: 
     Read "nint" long long integers from configuration file, 
     and put in array "list." 
*****/
{
    int  i, scanStat ;

    for ( i=0; i<nint; i++)
      {  scanStat = fscanf(fp,"%lld", &list[i] );  }

}  // end of function "readlong" 




// ************************************************
void readfloat(FILE *fp, int nint, float *list)   
/************* 
   Created 22-JAN-96 by A. Roodman
   From READINT

   ABSTRACT: 
     Read "nint" floats from configuration file, 
     and put in array "list." 

   April 2019: abort reading string by mistake.

***************/
{
    char c_get[80];
    int  i, scanStat, fstat ;
    float ftmp;
    char fnam[] = "readfloat";

    for ( i=0; i<nint; i++) {
      scanStat = fscanf(fp,"%s",c_get);         // read next string 
      ftmp     = NOFLOAT  ;
      fstat    = sscanf ( c_get, "%f", &ftmp );
      if ( ftmp == NOFLOAT ) {
	sprintf(c1err,"Could not read float from string='%s' ; ", c_get);
	sprintf(c2err,"reading item %d of %d items.", i+1, nint);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      list[i] = ftmp;
    }

}  // end of function "readfloat"


// ************************************************
void readdouble(FILE *fp, int nint, double *list)   
/************* 
   Feb 2009 R.Kessler
   ABSTRACT: 
     Read "nint" doubles from configuration file, 
     and put in array "list." 

  Apr 2019: abort reading string by mistake.

***************/
{
    char c_get[80];
    int  i,scanStat, dstat ;
    double dtmp;
    char fnam[] = "readdouble";

    for ( i=0; i<nint; i++) {
      scanStat = fscanf(fp,"%s",c_get);         // read next string 
      dtmp     = NODOUBLE ;
      dstat    = sscanf ( c_get, "%le", &dtmp );

      if ( dtmp == NODOUBLE ) {
	sprintf(c1err,"Could not read double from string='%s' ; ", c_get);
	sprintf(c2err,"reading item %d of %d items.", i+1, nint);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }

      list[i] = dtmp;
    }

}  // end of function "readdouble"




// ************************************************
void readchar(FILE *fp, char *clist)    
/*****
   Created 7-feb-95 by R.Kessler
   ABSTRACT: 
     Read next char. string from configuration file, 
     and put in the array "clist". 

******/
{
  int scanStat ;
  char fnam[] = "readchar";
  scanStat = fscanf(fp,"%s",clist);   // read next string 
  if ( clist[0] == '\t' ) {
    sprintf(c1err,"Invalid tab read in string '%s'", clist);
    sprintf(c2err,"Check file being read.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

}  // end of function "readchar"


// ******************************************************
void print_full_command(FILE *fp, int argc, char** argv) {

  int i;
  fprintf(fp,"\n Full command: ");
  for ( i=0; i < argc; i++ ) {
    fprintf(fp,"%s ", argv[i] );
  }
  fprintf(fp,"\n\n"); fflush(fp);
}

// ************************************************
void print_banner (const char *banner ) {
  printf("\n ********************************************"
	 "********************** \n");
  printf("   %s  \n" , banner );   fflush(stdout);
}
void fprint_banner (FILE *FP, const char *banner ) {
  fprintf(FP, "\n ********************************************"
	 "********************** \n");
  fprintf(FP, "   %s  \n" , banner );   fflush(FP);
}

void debugexit(char *string) {
  printf("\n xxx DEBUG EXIT: %s \n", string);
  fflush(stdout);
  exit(1);
}


// ***************************************
void print_debug_malloc(int opt, char *fnam) {

  // print debug statement if |opt| > 1
  char what[12];

  if ( abs(opt) <= 1 ) { return; }
  if ( opt > 0 )
    { sprintf(what,"alloc"); }
  else
    { sprintf(what,"free"); }
  
  printf(" DEBUG_MALLOC-%s for %s\n", what, fnam);
  fflush(stdout);

}  // end print_debug_malloc  


float malloc_shortint2D(int opt, int LEN1, int LEN2, short int ***array2D ) {
  // Created Sep 2021
  // Malloc array2D[LEN1][LEN2]  (intended for LEN1=NSN, LEN2=NCLPAR)
  float f_MEMTOT = 0.0 ;
  long long MEMTOT=0, i1 ;
  int MEM1 = LEN1 * sizeof(short int*); 
  int MEM2 = LEN2 * sizeof(short int);
  char fnam[] = "malloc_shortint2D";
  // ----------- BEGIN -------------

  print_debug_malloc(opt,fnam);

  if ( opt > 0 ) {

    *array2D = (short int**) malloc(MEM1) ; MEMTOT += MEM1;
    for(i1=0; i1< LEN1; i1++ ) {
      (*array2D)[i1] = (short int*) malloc(MEM2) ; MEMTOT += MEM2;
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1 < LEN1; i1++ ) { free((*array2D)[i1]); }
    free(array2D[i1]) ;    
  }

  return(f_MEMTOT);
} // end malloc_shortint2D

float malloc_double2D(int opt, int LEN1, int LEN2, double ***array2D ) {
  // Created Jun 11 2019
  // Malloc array2D[LEN1][LEN2]  (intended for LEN1=NSN, LEN2=NCLPAR)
  float f_MEMTOT = 0.0 ;
  long long MEMTOT=0, i1 ;
  int MEM1 = LEN1 * sizeof(double*); 
  int MEM2 = LEN2 * sizeof(double);
  char fnam[] = "malloc_double2D";
  // ----------- BEGIN -------------

  print_debug_malloc(opt,fnam);

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

float malloc_double3D(int opt, int LEN1, int LEN2, int LEN3, 
		      double ****array3D ) {
  // Created Jun 11 2019
  // Malloc array3D[LEN1][LEN2][LEN3] 
  //   (intended for LEN1=NSN, LEN2=MXa, LEN3=MXb)

  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1, i2 ;
  int MEM1 = LEN1 * sizeof(double**); 
  int MEM2 = LEN2 * sizeof(double*);
  int MEM3 = LEN3 * sizeof(double);
  char fnam[] = "malloc_double3D";
  // ----------- BEGIN -------------

  print_debug_malloc(opt,fnam);

  if ( opt > 0 ) {

    *array3D = (double***) malloc(MEM1) ; MEMTOT+=MEM1;
    for(i1=0; i1<LEN1; i1++ ) {
      (*array3D)[i1] = (double**) malloc(MEM2) ; MEMTOT+=MEM2;
      for(i2=0; i2<LEN2; i2++ ) {
	(*array3D)[i1][i2] = (double*) malloc(MEM3); MEMTOT+=MEM3;
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1<MEM1; i1++ ) {
      for(i2=0; i2<LEN2; i2++ ) {  free( (*array3D)[i1][i2]); }
      free( (*array3D)[i1]) ;
    }
    free(array3D);
  }


  return(f_MEMTOT);

} // end malloc_double3D


float malloc_float3D(int opt, int LEN1, int LEN2, int LEN3, 
		      float ****array3D ) {
  // Created Jan 2021
  // Malloc array3D[LEN1][LEN2][LEN3] 

  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1, i2 ;
  int MEM1 = LEN1 * sizeof(float**); 
  int MEM2 = LEN2 * sizeof(float*);
  int MEM3 = LEN3 * sizeof(float);
  char fnam[] = "malloc_float3D";
  // ----------- BEGIN -------------

  print_debug_malloc(opt,fnam);

  if ( opt > 0 ) {

    *array3D = (float***) malloc(MEM1) ; MEMTOT+=MEM1;
    for(i1=0; i1<LEN1; i1++ ) {
      (*array3D)[i1] = (float**) malloc(MEM2) ; MEMTOT+=MEM2;
      for(i2=0; i2<LEN2; i2++ ) {
	(*array3D)[i1][i2] = (float*) malloc(MEM3); MEMTOT+=MEM3;
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1<MEM1; i1++ ) {
      for(i2=0; i2<LEN2; i2++ ) {  free( (*array3D)[i1][i2]); }
      free( (*array3D)[i1]) ;
    }
    free(array3D);
  }


  return(f_MEMTOT);

} // end malloc_float3D


float malloc_double4D(int opt, int LEN1, int LEN2, int LEN3, int LEN4,
		      double *****array4D ) {
  // Created July 2019
  // Malloc array3D[LEN1][LEN2][LEN3][LEN4] 
  //   (intended for LEN1=NSN, LEN2=MXa, LEN3=MXb, LEN4=MXg)

  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1, i2, i3 ;
  int MEM1 = LEN1 * sizeof(double***); 
  int MEM2 = LEN2 * sizeof(double**);
  int MEM3 = LEN3 * sizeof(double*);
  int MEM4 = LEN4 * sizeof(double);
  char fnam[] = "malloc_double4D";

  // ----------- BEGIN -------------

  print_debug_malloc(opt,fnam);

  if ( opt > 0 ) {

    *array4D = (double****) malloc(MEM1) ; MEMTOT+=MEM1;
    for(i1=0; i1<LEN1; i1++ ) {
      (*array4D)[i1] = (double***) malloc(MEM2) ; MEMTOT+=MEM2;
      for(i2=0; i2<LEN2; i2++ ) {
	(*array4D)[i1][i2] = (double**) malloc(MEM3); MEMTOT+=MEM3;
	for(i3=0; i3<LEN3; i3++ ) {
	  (*array4D)[i1][i2][i3] = (double*) malloc(MEM4); MEMTOT+=MEM4;
	}
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1<MEM1; i1++ ) {
      for(i2=0; i2<LEN2; i2++ ) {
	for(i3=0; i3<LEN3; i3++ ) {
	  free( (*array4D)[i1][i2][i3] ); 
	}
	free( (*array4D)[i1][i2]); 
      }
      free( (*array4D)[i1]) ;
    }
    free(array4D);
  }


  return(f_MEMTOT);

}   // end malloc_double4D


float malloc_shortint4D(int opt, int LEN1, int LEN2, int LEN3, int LEN4,
			short int *****array4D ) {
  // Created Sep 2021
  // Malloc array3D[LEN1][LEN2][LEN3][LEN4] 
  //   (intended for LEN1=NSN, LEN2=MXa, LEN3=MXb, LEN4=MXg)

  float f_MEMTOT = 0.0 ;
  int MEMTOT=0, i1, i2, i3 ;
  int MEM1 = LEN1 * sizeof(short int***); 
  int MEM2 = LEN2 * sizeof(short int**);
  int MEM3 = LEN3 * sizeof(short int*);
  int MEM4 = LEN4 * sizeof(short int);
  char fnam[] = "malloc_shortint4D";

  // ----------- BEGIN -------------

  print_debug_malloc(opt,fnam);

  if ( opt > 0 ) {

    *array4D = (short int****) malloc(MEM1) ; MEMTOT+=MEM1;
    for(i1=0; i1<LEN1; i1++ ) {
      (*array4D)[i1] = (short int***) malloc(MEM2) ; MEMTOT+=MEM2;
      for(i2=0; i2<LEN2; i2++ ) {
	(*array4D)[i1][i2] = (short int**) malloc(MEM3); MEMTOT+=MEM3;
	for(i3=0; i3<LEN3; i3++ ) {
	  (*array4D)[i1][i2][i3] = (short int*) malloc(MEM4); MEMTOT+=MEM4;
	}
      }
    }

    f_MEMTOT = (float)(MEMTOT)/1.0E6;
    return(f_MEMTOT);
  } 
  else {  
    for(i1=0; i1<MEM1; i1++ ) {
      for(i2=0; i2<LEN2; i2++ ) {
	for(i3=0; i3<LEN3; i3++ ) {
	  free( (*array4D)[i1][i2][i3] ); 
	}
	free( (*array4D)[i1][i2]); 
      }
      free( (*array4D)[i1]) ;
    }
    free(array4D);
  }


  return(f_MEMTOT);

}   // end malloc_shortint4D
