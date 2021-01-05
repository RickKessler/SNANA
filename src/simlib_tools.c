/*******************************************
  simlib_tools.c:  Created April 4, 2007 by R.Kessler


  Public tools to write simulation library in standard
  format for snlc_sim program.
  Does lots of error checking and aborts on any
  hint of trouble !   You are not allowed to pass 
  nonsense values to this these functions.
  When you pass "NOBS" ovservations, it will then count
  how many times you call simlib_add_mjd to make sure
  that you add the correct number of epochs.

  Here is how to use these four functions.
  Make sure to read comments at the top of each function
  so that you know what arguments to pass.

  1. call simlib_open() to open file and write global info.

  2. for each RA,DECL position where a supernova could have been
     observed, call function
          simlib_add_header(...)

  3. For RA,DECL above, loop over observed MJDs & filters and call
        simlib_add_mjd(...) 

  3b. After last enry, call 
          simlib_add_header(-1, ...)
      to leave  END marker.

  4. call simlib_close() at the end of job.

  
      HISTORY
     ~~~~~~~~~~

  May 20,. 2008: add char *telescope argument to simlib_open

  May 27, 2009: define new global CFILT_LAST so that both MHJ & filter
                are compared to determine new MJD.  Handles case where
                MJD for different filters is the same to given precision.

  Jun 8, 2009: add FIELD as argument to simlib_add_header(..).

  Feb 25, 2012: add HEADERFILE argument to simlib_open to paste in
                arbitrary header info from another file.

  Aug 31, 2012: "NOBS" check 1-1400 changed to 1-3000

  Jun 20 2017: add MJD arg to simlib_add_mjd() so that MJD is double
               rather than float.

  Apr 3 2018:
   + long int IDEXPT -> char STRINGID[20] to allow ID*NEXPOSE
     e.g., 123438*2

 Jan 5 2021: allow PSF max to 50 (was 10) [for LSST]

********************************************/


void simlib_open(char *filename, char *surveyname, char *filters, 
		 char *telescope, char *comment, char *headFile );

void simlib_add_header(int optflag, int IDLIB, int NOBS, 
		       char *FIELD, float *INFO );

void simlib_add_mjd(int opt, double MJD, char * STRINGID, char *FILTNAME, float *INFO);


void simlib_close(void);


int  CHECK_LIBVAL(char *varname, float value, float varmin, float varmax);
void MADABORT(void);
void PRINT_ERROR(char *msgerr );

FILE *FPLIB;      // global pointer to library file.

int   NSIMLIB        ;  // number of simlib entries          
int   NMJD_FOUND     ;  // number of add_mjd  calls; should equal NOBS
int   NMJD_EXPECT = 0;  // = NOBS
double    MJD_LAST;
// xxx mark delete long int IDEXPT_LAST;
char     STRINGID_LAST[20];
char     CFILT_LAST[2];

int OPT_CHECKVAL ;


// *************************************
void simlib_open(
		 char *filename     // full name of file to open
		 ,char *surveyname  // name of survey; i.e, "SDSS"
		 ,char *filters     // filter list; i.e, "ugriz"
		 ,char *telescope   // name of telescope
		 ,char *comment     // user comment
		 ,char *headFile    // optional file with header contents
		 ) {

  // Created Apr 6, 2007 by R.Kessler
  // Open library file and write global information 
  // May 20, 2008: add telescope as argument
  // Feb 25, 2012: add *headFile argument

  time_t tstart;

  FILE *fp_head ;
  char cline[200];
  char fnam[20] = "simlib_open";
  
  // -------------- BEGIN --------------


  NSIMLIB = 0;

  if ( (FPLIB = fopen(filename, "wt"))==NULL ) {       
    printf( "FATAL ERROR(%s): \n", fnam );
    printf( "Cannot open libfile :\n" );
    printf( " '%s' \n", filename);
    MADABORT();
  }

  printf("\n Opened sim-library output file : \n\t %s \n", filename);

  // write global info at top of lib-file

  fprintf(FPLIB,"SURVEY: %s     FILTERS: %s   TELESCOPE:  %s \n", 
	  surveyname, filters, telescope );
  fprintf(FPLIB,"USER: %s     HOST: %s \n", 
	  getenv("USER"), getenv("HOST") );
  fprintf(FPLIB,"COMMENT: '%s'  \n", comment);

  // Check or optional header file

  if ( strlen(headFile) > 0 ) {
    if ( (fp_head = fopen(headFile, "rt"))==NULL ) {       
      printf( "FATAL ERROR(%s): \n", fnam );
      printf( "Cannot open headFile: headFile :\n" );
      printf( " '%s' \n", headFile);
      MADABORT();
    }

    printf("\t Extract header contents from : %s \n", headFile);

    // copy each line of header file into simlib header.
    fprintf(FPLIB,"COMMENT: Header below extracted from %s \n", headFile);
    while ( fgets (cline, 100, fp_head) !=NULL  ) 
      {      fprintf(FPLIB,"%s", cline);    }

    fclose(fp_head);

  }  // end of headFile 

  time(&tstart);
  fprintf(FPLIB,"\nBEGIN LIBGEN  %s \n", ctime(&tstart) ) ;

  fprintf(FPLIB,"\n");

  fflush(FPLIB);


}  // end of function simlib_open
 

// *****************************************
void simlib_add_header(
		       int optflag     // (I) option for debug (0=nominal)
		       ,int   IDLIB    // (I) incremental lib id (1,2,3 ... )
		       ,int   NOBS     // (I) # obs  to follow (must be >=3)
		       ,char *FIELD    // (I) name of field
		       ,float *INFO    // (I) header info; see unpacking below
		       ) {
  /*****
   Created April 2007 by R.Kessler
   Write header for new SN event at new RA,DECL position.
   Unpacking comments below specify required or optional for each 
   *INFO element.

   "optional" means that only positive values are used.
   If optional Z or PEAKMJD  are positive, these values will be used
   in the simulation instead of picking them randomly.

   optflag = 0 => nominal (recommended)

   optflag = 1 => skip tests on variables to allow garbage written to file.
                     (for debugging)


   optflag = -1  leave END_LIBID marker after all epochs 
                      (all other args ignored) 

  Jun  8, 2009: pass FIELDNAME to include in header

  Sep 22 2013: write MWEBV ad .4f instead of .3f to have better
               precision on verifying that it's from SFD.

  *****/


  char fnam[20] = "simlib_add_header";
  char msgerr[100];
  int istat, NOPT;

  float RA, DECL, MWEBV, PIXSIZE, Z, PEAKMJD ;

  // ------------ BEGIN ------------

  OPT_CHECKVAL = 1;  // default: make checks on variable values

  if ( optflag == 1 ) 
    OPT_CHECKVAL = 0; // turn off sanity checks (for debugging)



  // check number of MJDs for last lib entry
  //  if ( NMJD_EXPECT > 0 &&  NMJD_FOUND != NMJD_EXPECT )

  if ( optflag == -1 ) {

    fprintf(FPLIB,"END_LIBID: %d \n", IDLIB );

    if ( NMJD_FOUND != NMJD_EXPECT ) {

      PRINT_ERROR("\n");

      sprintf(msgerr,"  FATAL ERROR in %s \n", fnam);
      PRINT_ERROR(msgerr);

      sprintf(msgerr,"\t Found %d MJDs, but expect %d from user NOBS arg. \n", 
	      NMJD_FOUND, NMJD_EXPECT);
      PRINT_ERROR(msgerr);

      sprintf(msgerr,"\t Check last entry in your library file. \n");
      PRINT_ERROR(msgerr);
      
      MADABORT();
    }
    return;
  }


  NMJD_EXPECT = NOBS;
  NMJD_FOUND  = 0;
  MJD_LAST    = 0.0 ;
  // xxx mark delete   IDEXPT_LAST = 0;
  STRINGID_LAST[0] = 0 ;
  sprintf(CFILT_LAST,"?");

  // unpack header info

  RA      = *(INFO+0) ;  // required: right-ascension, degrees
  DECL    = *(INFO+1) ;  // required: declination, degrees
  MWEBV   = *(INFO+2) ;  // required: MilkyWay extinction = AV/RV
  PIXSIZE = *(INFO+3) ;  // required: pixel size, arcseconds
  Z       = *(INFO+4) ;  // optional: redshift
  PEAKMJD = *(INFO+5) ;  // optional: MJD at maximum B-band luminosity

  // start by writing required info.

  fprintf(FPLIB, "# -------------------------------------------- \n");
  fprintf(FPLIB, "LIBID: %d \n", IDLIB);

  fprintf(FPLIB, "RA: %f   DECL: %f   NOBS: %d    MWEBV: %.4f   PIXSIZE: %5.3f \n", 
	  RA, DECL,  NOBS, MWEBV, PIXSIZE );

  if ( strlen(FIELD) > 0 ) 
    fprintf(FPLIB,"FIELD: %s \n", FIELD );

  // now check for optional info
  NOPT=0;
  if ( Z       >= 0.0 ) {NOPT++; fprintf(FPLIB,"REDSHIFT: %6.4f    ", Z); }
  if ( PEAKMJD >= 0.0 ) {NOPT++; fprintf(FPLIB,"PEAKMJD:  %9.3f    ", PEAKMJD); }
  if ( NOPT > 0 ) fprintf(FPLIB,"\n");


  if ( optflag == 1 && NOBS == 0 ) {
    fprintf(FPLIB, "\t WARNING: ZERO OBSERVATIONS !!! \n");
    fflush(FPLIB);
    return ;
  }


  // make table header to that file is self-documented.

  fprintf(FPLIB, "\n");
  fprintf(FPLIB,"#                           CCD  CCD         PSF1 PSF2 PSF2/1                    \n");
  fprintf(FPLIB,"#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG \n");



  NSIMLIB++ ;

  // make sanity checks on input values.

  if ( OPT_CHECKVAL == 1 ) {
    istat = CHECK_LIBVAL("IDLIB",   (float)IDLIB,  0.0, 100000. );
    istat = CHECK_LIBVAL("NOBS",    (float)NOBS,   1.0, 3000.   );
    istat = CHECK_LIBVAL("RA",      RA,           -200., 400.0  );
    istat = CHECK_LIBVAL("DECL",    DECL,         -200., 400.0  );
    istat = CHECK_LIBVAL("MWEBV",   MWEBV,         0.0, 2.0     );
    istat = CHECK_LIBVAL("PIXSIZE", PIXSIZE,       0.0, 2.0     );
  }

} // end of simlib_add_header



// *********************************************
void simlib_add_mjd(
		    int opt            // 1=>search info;  2=> template info
		    ,double MJD        // MJD
		    ,char *STRINGID    // exposure ID/RUN/FIELD/whatever
		    ,char *FILTNAME    // filter name, abbreviation
		    ,float *INFO       // see unpacking comments below
		    ) {

  /*****
   Created April 2007 by R.Kessler
   Write info related to observed MJD; can then be used by snlc_sim.
   Here are a few comments:

      - OPT =  1 for search image is required
      - OPT =  2 for template image is optional.

      - SKYSIG (in measured ADU) is for search or template image.
        However, if you do not give template info, then SKYSIG(search)
        shuold include your best estimate of template noise when doing
        galaxy-subtraction.

      - If you only want single Guassian PSF, then set PSF = { PSF, 0.0, 0.0 }

      - to ignore MAGOBS, set it to -99.0; otherwise this magnitude is used
        in the simulation rather than generating a mag from a model.
        This feature allows you to over-ride the simulated mags with your 
        own mags.

      - Note that the units are measured ADU and pixels.
        Since you also provide CCDGAIN (electrons/ADU) and pixel size in
        arcseconds, the simulation will make the necessary conversions.

   May 27, 2009:  ISNWMJD requires different MJD and different filter.
                  See new global declaration CFILT_LAST[2].

   Oct 25, 2010: %10d -> %10ld  for IDEXPT

   Jan 6, 2011: increase SKYSIG CHECK_LIBVAL-limit from 400 to 600 
                so that it works for LSST coadd.

  Jun 20 2017: 
    + SKYSIG ABORT = 1000 -> 2000 (for LSST-DDF)
    + add MJD argument

  *****/


  char key[2];
  char fnam[20] = " simlib_add_mjd" ;
  int istat;
  int ISNEWMJD, ISNEWID ;

  // define local variables for unpacking

  float PSF[3], ZPT[2], SKYSIG, CCDGAIN, CCDNOISE, MAGOBS;

  // --------------- BEGIN --------------

  // unpack *INFO array
  // xxx mark delete   MJD      = *(INFO+0) ;  // 53000 in year 2006
  CCDGAIN  = *(INFO+1) ;  // electrons per ADU
  CCDNOISE = *(INFO+2) ;  // CCD read noise in electrons
  SKYSIG   = *(INFO+3) ;  // skynoise in ADU per pixel
  PSF[0]   = *(INFO+4) ;  // PSF-sigma (pixels) for inner Gaussian
  PSF[1]   = *(INFO+5) ;  // PSF-sigma (pixels) for outer Gaussian
  PSF[2]   = *(INFO+6) ;  // PSF(outer)/PSF(inner) ratio at origin
  ZPT[0]   = *(INFO+7) ;  // zero point
  ZPT[1]   = *(INFO+8) ;  // zero point error, or spread
  MAGOBS   = *(INFO+9) ;  // model of observed mag (optional)


  if ( opt == 1 ) {
    sprintf(key, "S") ;

    ISNEWMJD = ISNEWID = 0 ;
    if ( MJD != MJD_LAST || strcmp(FILTNAME,CFILT_LAST) != 0 ) 
      { ISNEWMJD = 1 ; }

    if ( strcmp(STRINGID,STRINGID_LAST) != 0 ) 
      { ISNEWID = 1; }

    if ( ISNEWMJD == 1  || ISNEWID )  {
      NMJD_FOUND++ ;
    }

    MJD_LAST = MJD ;
    sprintf(CFILT_LAST,"%s", FILTNAME);
    // xxx mark delete     IDEXPT_LAST = IDEXPT;
    sprintf(STRINGID_LAST, "%s", STRINGID);
  }
  else if ( opt == 2 ) 
    sprintf(key, "T") ;
  else {
    printf("ERROR: opt=%d is invalid for function %s \n", opt, fnam);
    MADABORT();
  }


  fprintf(FPLIB,"%s: "
	  "%9.4f %10.10s %s %5.2f %5.2f %6.2f "
	  "%4.2f %4.2f %5.3f "
	  "%6.2f %6.3f"
	  , key
	  , MJD, STRINGID, FILTNAME, CCDGAIN, CCDNOISE, SKYSIG
	  , PSF[0], PSF[1], PSF[2]
	  , ZPT[0], ZPT[1]
	  );


  // print MAGOBS for "Search-image" option only
  if ( opt == 1 ) fprintf(FPLIB," %7.3f", MAGOBS );

  fflush(FPLIB);

  // add <CR>
  fprintf(FPLIB, "\n");

  // make sanity checks after adding entry to file

  if ( OPT_CHECKVAL == 1 ) {
    istat = CHECK_LIBVAL("MJD",         MJD,     20000., 80000. );
    istat = CHECK_LIBVAL("CCDGAIN",     CCDGAIN, 0.0, 100. );
    istat = CHECK_LIBVAL("CCDNOISE",    CCDNOISE,0.0, 100. );
    istat = CHECK_LIBVAL("SKYSIG",      SKYSIG,  0.0,2000. );
    istat = CHECK_LIBVAL("PSF(inner)",  PSF[0], 0.0, 10. );
    istat = CHECK_LIBVAL("PSF(outer)",  PSF[1], 0.0, 10. );
    istat = CHECK_LIBVAL("PSF-ratio",   PSF[2], 0.0, 10. );
    
    // xxx    istat = CHECK_LIBVAL("ZeroPoint",       ZPT[0], 10.0, 40. );
    istat = CHECK_LIBVAL("ZeroPoint",       ZPT[0], 50.0, 40. );
    istat = CHECK_LIBVAL("ZeroPoint-sigma", ZPT[1],  0.0, 4.0 );
  }

}  // end of simlib_add_mjd


// ************************************************
void simlib_close(void) {

  fprintf(FPLIB, "\n");
  fprintf(FPLIB, "END_OF_SIMLIB: %d ENTRIES \n", NSIMLIB );

  fclose(FPLIB) ;
}  

// ************************************************
int CHECK_LIBVAL(char *varname, 
		 float value, float varmin, float varmax) {

  // check that 'value' is between varmin and varmax;
  // if not, then print message and MADABORT

  char msgerr[10][80];
  int i;

  // ---------------- BEGIN --------------

  if ( value < varmin || value > varmax ) {
    sprintf(msgerr[0],"\n ERROR WRITING LIBRARY: \n");
    sprintf(msgerr[1],"\t %s = %f is invalid. \n", varname, value );
    sprintf(msgerr[2],"\t %s valid range is %f to %f \n", 
	    varname, varmin, varmax );
    sprintf(msgerr[3],"\t Check last entry in your library file. \n");

    for ( i=0; i<=3; i++ ) {
      PRINT_ERROR ( msgerr[i] );
    }

    MADABORT();
  }

  return 0 ;


} // end of CHECK_LIBVAL


// ***************************************
void PRINT_ERROR(char *msgerr ) {

  printf("%s", msgerr );
  fprintf(FPLIB, "%s", msgerr );

}  // end of PRINT_ERROR

// ************************************************
void MADABORT(void) {
   char cmsg[40] = { "ABORT program on Fatal Error." };

   printf("\n");
   printf("\n   `|```````|`    ");
   printf("\n   <| o\\ /o |>    ");
   printf("\n    | ' ; ' |     ");
   printf("\n    |  ___  |     %s ", cmsg);
   printf("\n    | |' '| |     ");
   printf("\n    | `---' |     ");
   printf("\n    \\_______/    ");
   printf("\n");

   printf("\n");   exit(1);

}    //  end of "MADABORT"  
