/************************************
  Created April 2013 by R.Kessler

  Program to merge ROOT files produces by 
  an snana job. THe merging is mainly needed for the
  'split_and_fit' output in which each job is split
  into many  sub-jobs and thus the output files need
  to be re-combined, or merged.



  Usage: 
     merge_root.exe <outFileList>  <mergeFile>

  If <mergeFile> already exists then abort.
  File Extensions .root/.ROOT are assumed to be ROOT.

  Nov 27 2020:
    + abort if number of input files exceeds MXFILE_MERGE
    + MXFILE_MERGE -> 200 (was 130)

**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#include "sntools.h"
#include "sntools_output.h"

#ifdef USE_ROOT
void MERGE_ROOT(int NFILE, char **INFILES, char *OUTFILE);
#endif


#define MXFILE_MERGE 200 

char msgerr1[80], msgerr2[80];

struct INPUTS {
  int  NFILE_IN ;
  char INFILES[MXFILE_MERGE][200] ;
  char OUTFILE[200];  // final merged file
  int  IFILETYPE ;
  char FILETYPE_NAME[20];
} INPUTS ;



void parse_args(int argc, char **argv) ;
void checkFiles(void);


// ==================================
int main(int argc, char **argv) {

  int NF, i ;
  char *inF[MXFILE_MERGE], *outF ;
  char  fnam[] = "merge_root" ;

  set_EXIT_ERRCODE(EXIT_ERRCODE_merge_root);

#ifndef USE_ROOT
  sprintf(msgerr1,
	  "Cannot run %s because it's not linked to root.", fnam);
  sprintf(msgerr2,
	  "Set USE_ROOT in Makefile and #define ROOT in sntools_output.h");
  errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 ); 
#endif


  printf(" Begin %s \n", fnam ); fflush(stdout);
  parse_args(argc,argv);

  // do lots of sanity checks
  checkFiles();

  NF   = INPUTS.NFILE_IN ;
  outF = INPUTS.OUTFILE ;
  for(i=0; i<NF; i++ ) { inF[i] = INPUTS.INFILES[i] ; }

#ifdef USE_ROOT
  MERGE_ROOT(NF,inF,outF); 
#endif

  return(0);

} // end of main



// ====================
void parse_args(int NARG, char **argv) {

  int NIN, i ;
  char fnam[] = "parse_args" ;

  // --------- BEGIN -----------

  INPUTS.NFILE_IN = NIN = 0;
  INPUTS.OUTFILE[0] = 0 ;

  for(i=1; i < NARG; i++ ) {


    if ( i < NARG-1 )  { 
      if ( NIN < MXFILE_MERGE ) 
	{ sprintf(INPUTS.INFILES[NIN],"%s", argv[i]) ; }
      NIN++ ; 
    }
    else
      { sprintf(INPUTS.OUTFILE,"%s", argv[i]) ; }
  }

  // abort if too many files (Nov 2020)
  if ( NIN >= MXFILE_MERGE ) {
    sprintf(msgerr1,"%d input files exceeds bound of MXFILE_MERGE=%d",
	    NIN, MXFILE_MERGE);
    sprintf(msgerr2,"Reduce number of files, or increase MXFILE_MERGE");
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  INPUTS.NFILE_IN = NIN ;
  
} // end of parse_args



// =========================
void checkFiles(void) {

  // sanity checks in the inputs
  // Make sure that
  // - at least 1 input files
  // - input files exist
  // - output file does not exist
  // - all file types are the same

  int NIN, jstat, ifile, IFILETYPE, IFILETYPE_FIRST ;
  char *inFile, *outFile, FTYPE[20] ;
  char fnam[] = "checkFiles" ;

  struct stat statbuf ;

  // ------------ BEGIN --------

  // Make sure we have at least 2 files to mege

  IFILETYPE_FIRST = -9;

  NIN = INPUTS.NFILE_IN ;
  if ( NIN < 1 ) {
    sprintf(msgerr1,"Only %d input files to merge.", NIN);
    sprintf(msgerr2,"At least 2 are required.");
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }

  // make sure that each input file exists using stat function

  for(ifile=0; ifile < NIN; ifile++ ) {
    inFile  = INPUTS.INFILES[ifile] ;
    jstat   = stat(inFile, &statbuf );
    IFILETYPE = -9 ;
    sprintf(FTYPE, "Unknown");

    if ( jstat == 0 )  { 

#ifdef USE_ROOT
      if ( ISFILE_ROOT(inFile)  ) { 
	IFILETYPE = IFILETYPE_ROOT  ; 
	sprintf(FTYPE,"ROOT");
      }
#endif

      printf(" Found %-8.8s File %3d: '%s' \n", FTYPE, ifile+1, inFile ); 
      fflush(stdout);

      if ( ifile == 0 ) { 
	IFILETYPE_FIRST  = IFILETYPE ; 
	INPUTS.IFILETYPE = IFILETYPE; 
	sprintf(INPUTS.FILETYPE_NAME, "%s", FTYPE);

	if ( IFILETYPE < 0 ) {
	  sprintf(msgerr1,"Invalid file type for '%s'", inFile);
	  sprintf(msgerr2,"Must be a ROOT file.");
	  errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );	
	}
      }

      if ( IFILETYPE != IFILETYPE_FIRST ) {
	sprintf(msgerr1,"'%s' is bad file type", FTYPE);
	sprintf(msgerr2,"for in-file = '%s'", inFile);
	errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );	
      }

    }
    else {
      sprintf(msgerr1,"Cannot find requested input file:") ;
      sprintf(msgerr2,"%s", inFile);
      errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
    }
  }


  // make sure that output file does NOT exist
  outFile  = INPUTS.OUTFILE ;
  jstat    = stat(outFile, &statbuf );
  printf(" merged output --> '%s' \n", outFile) ; fflush(stdout);
  if ( jstat == 0  ){
    sprintf(msgerr1,"Output file already exists ?!?!?!");
    sprintf(msgerr2,"See '%s' ", outFile);
    errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
  }


  IFILETYPE = -9 ;
#ifdef USE_HBOOK
      if ( ISFILE_HBOOK(outFile) ) { IFILETYPE = IFILETYPE_HBOOK ;  }
#endif

#ifdef USE_ROOT
      if ( ISFILE_ROOT(outFile)  ) { IFILETYPE = IFILETYPE_ROOT  ; }
#endif

      if ( IFILETYPE != INPUTS.IFILETYPE ) {
	sprintf(msgerr1,"Output file type does not match input.");
	sprintf(msgerr2,"Output file type must be '%s'",
		INPUTS.FILETYPE_NAME);
	errmsg(SEV_FATAL, 0, fnam, msgerr1, msgerr2 );
      }

} // end of checkFiles

