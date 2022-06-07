/**************************************************************
 Created Oct 14, 2009 by R.Kessler

 Generic SED model based on simulations of the SN explosion.
 Takes advantage of the SALT2 infrastructure for initializing
 filters, primary reference (i.e, BD17), and the flux-integrals
 for each passbands, redshift, epoch and SED surface.
 All fluxes are computed in the observer-frame.

  HISTORY       

 Mar 1 2017: pass RV_host and AV_host to genmag_SIMSED() to allow
             analytcally computed host-galaxy extinction.

 Jul 28 2017: remove obsolete binary with warning instead of abort.
              Makes SIMSED easier to use.

 July 30 2017: 
   + refactor with new functions open_SEDBINARY and open_TABBINARY.
     Needed to add IVERSION and MXDAY into binary file so that
     for binary read mode there is no need to read ASCII-SED files
     to get MXDAY.  Older binaries are no longer supported.
     New binaries will have IVERSION so that future changes can 
     be made backward-compatible.

   + set default LOGZINB=.02 and parse optional LOGZBIN: <xxx>
     from SED.INFO file.

 July 20 2018:
  +  replace interp_flux_SEDMODEL calls with get_flux_SEDMODEL;
     latter takes care of extrapolation for Trest outside range of model.
  + in genmag_SIMSED, (OPTMAX & 8)>0 -->  debug dump.

  Jul 27 2018
   + if UVLAM extrap is set, set SEDMODEL.MINLAMFILT= UVLAM

  Jul 30 2018
    + add INDEX_SED output argument

  Aug 26 2018: read MINSLOPE_EXTRAP_LATE from SED.INFO file

  Dec 20 2018: refactor using loops from 0 to SEDMODEL.NPAR-1
                (instead of 1 to SEDMODEL.NPAR)

  Apr 28 2019: set Lrange_SIMSED when reading "RESTLAMBDA_RANGE:" key.
               --> affects value of T0shiftPeak, bolometric flux,
                   and memory used to read.

 Dec 15 2020:
   + checkBinary_SIMSED aborts if in batch mode and SED.BINARY needs
     to be remade. Avoids conflict among multiple batch jobs.

 Dec 14 2021: check new OPTMASK options to force creation of binaries.
              See FORCE_SEDBINARY and FORCE_TABBINARY

 Mar 02 2022: fix bug so that UVLAM_EXTRAP works when reading binary file
              or original text files.

*************************************/

#include  <stdio.h> 
#include  <math.h>     
#include  <stdlib.h>   
#include  <sys/stat.h>

#include  "sntools.h"           // SNANA community tools
#include  "genmag_SEDtools.h"
#include  "genmag_SIMSED.h"
#include  "MWgaldust.h"

const char dual_bits_SIMSED[INTERP_SIMSED_MAX_DIM] =
  {1, 2, 4, 8, 16, 32, 64, 128};

// =======================================================
// define mangled functions with underscore (for fortran)


int init_genmag_simsed__(char *version, char *PATH_BINARY, char *SURVEY,
			 char *kcorFile, int *OPTMASK) {
  int istat;
  istat = init_genmag_SIMSED(version, PATH_BINARY, SURVEY, kcorFile, *OPTMASK);
  return istat;
} 


void genmag_simsed__(int *OPTMASK, int *ifilt, double *x0,
		     int *NLUMIPAR,int *iflagpar,int *iparmap, double *lumipar,
		     double *RV_host, double *AV_host, double *mwebv, 
		     double *z, int *nobs, double *Tobs_list, 
		     double *magobs_list, double *magerr_list,
		     int *index_sed ) {

  genmag_SIMSED(*OPTMASK, *ifilt, *x0, *NLUMIPAR, iflagpar, iparmap, lumipar,
		*RV_host, *AV_host,
		*mwebv, *z, *nobs, Tobs_list, magobs_list, magerr_list,
		index_sed );
}



/****************************************************************
  init_genmag_SIMSED:
    o reads in the filters from FilterFiles
    o calculates filter mean and AB zeropoint
    o reads the templates from TemplateFiles
    o returns true if successful

   Note: must call init_filter_SEDMODEL() and init_primary_SEDMODEL()
          before calling this function.


  Feb 22 2017: allow input *VERSION to incude user-defined path.

  Aug 9 2017: 
    + Lrange_SIMSED[1] -> 30000 (was 20000)
    + Trange_SIMSED -> -100 to +500 (was -20 to +200)

 Apr 28 2018: 
   + always call shiftPeakDay_SEDMODEL(), not just when T0>0.
   + for SED.BINARY, fix buggy fwrite under if(WRFLAG_SEDBINARY).
       sizeof(int*)   -> sizeof(int) 
       sizeof(float*) -> sizeof(float) 

 May 17 2018:
   Shift T=0 to be at peak or at time of explosion (see OPTMASK_T0SHIFT).
   Beware: latter option tested for NON1ASED, but not for SIMSED.

 May 27 2018:
   +  Trange_SIMSED -> -500 to +500 (was -100 to +500)

 July 11 2018: check OPTMASK_T0SHIFT_PEAKMAG

 July 26 2018: Trange_SIMSED -> 1500 (was 500)

****************************************************************/
int init_genmag_SIMSED(char *VERSION      // SIMSED version
		       ,char *PATH_BINARY // directory to write/read binaries
		       ,char *SURVEY      // name of survey  
		       ,char *kcorFile    // kcor filename
		       ,int OPTMASK       // bit-mask of options
		       ) {   

  // OPTMASK +=  1 --> create binary file if it doesn't exist
  // OPTMASK +=  2 --> force creation of SED.BINARY
  // OPTMASK +=  4 --> force creaton of flux-table binary
  // OPTMASK += 64 --> test mode only, no binary, no time-stamp checks
  // OPTMASK += 128 -> batch mode, thus abort on stale binary
  //
  // Aug 18 2020: check kcor file with .gz
  // Dec 14 2021: new OPTMASK 2 and 4
  // Mar 02 2022: check UVLAM_EXTRAP

  int NZBIN, IZSIZE, ifilt, ifilt_obs, ised, istat;
  int retval = SUCCESS ;

  bool USE_BINARY=false, USE_TESTMODE=false;
  bool FORCE_SEDBINARY=false, FORCE_TABBINARY=false;

  char
    BANNER[120]
    ,tmpFile[MXPATHLEN]
    ,sedFile[MXPATHLEN]
    ,bin1File[MXPATHLEN]  // SED binary file
    ,bin2File[MXPATHLEN]  // flux-table binary file
    ,sedcomment[40], version[60], kcorFile_gz[MXPATHLEN]
    ;

  FILE *fpbin1, *fpbin2 ;

  struct stat statbuf ; // to check if BINARY dir exists
  char fnam[] = "init_genmag_SIMSED" ;

  // -------------- BEGIN --------------

  sprintf(BANNER, "%s : Initialize SIMSED Light Curve Model", fnam );
  print_banner(BANNER);

  // get local logicals for bit-mask  options
  USE_BINARY      = ( OPTMASK &  OPTMASK_INIT_SIMSED_BINARY   ) > 0 ;
  USE_TESTMODE    = ( OPTMASK &  OPTMASK_INIT_SIMSED_TESTMODE ) > 0 ;
  ISBATCH_SIMSED  = ( OPTMASK &  OPTMASK_INIT_SIMSED_BATCH    ) > 0 ;

  if ( (OPTMASK & OPTMASK_INIT_SIMSED_BINARY1)> 0 )
    { FORCE_SEDBINARY = true; USE_BINARY = true;  }
  if ( (OPTMASK & OPTMASK_INIT_SIMSED_BINARY2)> 0 )
    { FORCE_TABBINARY = true; USE_BINARY = true;  }

  if ( NFILT_SEDMODEL == 0  && !USE_TESTMODE ) {
    sprintf(c1err,"No filters defined ?!?!?!? " );
    sprintf(c2err,"Need to call init_filter_SEDMODEL");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  // get full path of kcor file ... as passed or in $SNDATA_ROOT/kcor
  if ( BINARYFLAG_KCORFILENAME && !USE_TESTMODE ) {
    sprintf(kcorFile_gz, "%s.gz", kcorFile);
    bool ACCESS_KCOR = ( access(kcorFile,   F_OK) == 0 || 
			 access(kcorFile_gz,F_OK) == 0 );

    if ( ACCESS_KCOR ) 
      { sprintf(SIMSED_KCORFILE,"%s", kcorFile); }
    else {
      sprintf(SIMSED_KCORFILE,"%s/kcor/%s", 
	      PATH_SNDATA_ROOT, kcorFile); 
    }
  }


  // summarize filter info
  if ( !USE_TESTMODE  ) {  filtdump_SEDMODEL(); }

  // ==========================================
  // construct path to SIMSED surfaces

  extract_MODELNAME(VERSION,
		    SIMSED_PATHMODEL, version); // returned

  // if no path, assume default path 
  if ( strlen(SIMSED_PATHMODEL) == 0 ) {
    sprintf( SIMSED_PATHMODEL, "%s/models/SIMSED/%s", 
	     getenv("SNDATA_ROOT"), version );
  }

  printf("  Read SIMSED model parameters from \n     %s\n",
	 SIMSED_PATHMODEL );

  printf("\n");

  sprintf(SIMSED_INFO_FILENAME_FULL, "%s/%s", 
	  SIMSED_PATHMODEL, SIMSED_INFO_FILENAME );
  
  // set defaults

  ISIMSED_SEQUENTIAL = 0 ;  // for sequential GRIDONLY option 

  SEDMODEL.NSURFACE   = 0 ;
  SEDMODEL.FLUXSCALE  = 1.0 ;
  SEDMODEL.MAGERR_FIX = 0.1 ;
  SEDMODEL.LOGZBIN    = LOGZBIN_SIMSED_DEFAULT ; 
  SEDMODEL.NPAR       = -1 ;
  SEDMODEL.OPTMASK    = 
    OPTMASK_DAYLIST_SEDMODEL +   // allow non-uniform day bins
    OPTMASK_T0SHIFT_PEAKMAG      // shift T=0 to peak (July 2018)
    ;


  // set extreme ranges to read anything
  Trange_SIMSED[0] = -500. ; 
  Trange_SIMSED[1] = 1500. ;
  Lrange_SIMSED[0] =  LAMMIN_SEDMODEL ;
  Lrange_SIMSED[1] =  LAMMAX_SEDMODEL ;

  // -------------------------------------------------
  // check for binary file to read SEDs much quicker.
  // If binary does not exist, then set flag to create binary for next time.

  SIMSED_BINARY_INFO.WRFLAG_SED  = false;
  SIMSED_BINARY_INFO.RDFLAG_SED  = false;
  SIMSED_BINARY_INFO.WRFLAG_FLUX = false;
  SIMSED_BINARY_INFO.RDFLAG_FLUX = false;

  if ( USE_BINARY ) {

    IVERSION_SIMSED_BINARY = WRVERSION_SIMSED_BINARY ;

    // make sure that PATH_BINARY exists, and copy it to global
    sprintf(SIMSED_BINARY_INFO.PATH, "%s", PATH_BINARY);
    istat = stat(PATH_BINARY, &statbuf);
    if ( istat != 0 ) {
      sprintf(c1err,"PATH_BINARY='%s'", PATH_BINARY);
      sprintf(c2err,"does not exist.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    sprintf(bin1File, "%s/%s", SIMSED_PATHMODEL, SIMSED_BINARY_FILENAME );
    sprintf(bin2File,"%s/%s_%s-%s.BINARY", 
	    PATH_BINARY, version, SURVEY, FILTLIST_SEDMODEL );

    open_SEDBINARY(bin1File, FORCE_SEDBINARY,&fpbin1, 
		   &SIMSED_BINARY_INFO.RDFLAG_SED, 
		   &SIMSED_BINARY_INFO.WRFLAG_SED);

    open_TABBINARY(bin2File, FORCE_TABBINARY, &fpbin2,
		   &SIMSED_BINARY_INFO.RDFLAG_FLUX, 
		   &SIMSED_BINARY_INFO.WRFLAG_FLUX);

  }

  // -------------------------------------- 
  malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL,0,0,0);
  read_SIMSED_INFO(SIMSED_PATHMODEL);
  dump_SIMSED_INFO();

  // - - - - - - - - - - - - - - - - - - - - 
  if ( USE_TESTMODE ) { return(retval); } // July 28 2018
  // - - - - - - - - - - - - - - - - - - - - 


  // determine SEDMODEL.MXDAY from ASCII or binary file
  set_SIMSED_MXDAY(SIMSED_PATHMODEL,fpbin1, 
		   SIMSED_BINARY_INFO.RDFLAG_SED, 
		   SIMSED_BINARY_INFO.WRFLAG_SED); 

  // check to change default logz binning
  set_SIMSED_LOGZBIN();

  // =======================================================
  // allocate memory for storing flux-integral tables
  NZBIN  = REDSHIFT_SEDMODEL.NZBIN ;
  NLAMPOW_SEDMODEL = 0 ;
  malloc_FLUXTABLE_SEDMODEL ( NFILT_SEDMODEL, NZBIN, NLAMPOW_SEDMODEL, 
			      SEDMODEL.MXDAY, SEDMODEL.NSURFACE );
  fflush(stdout);

  // ------- Now read the spectral templates -----------

  for ( ised = 1 ; ised <= SEDMODEL.NSURFACE ; ised++ ) {
    
    sprintf(tmpFile, "%s/%s", SIMSED_PATHMODEL, SEDMODEL.FILENAME[ised] );
    sprintf(sedcomment,"(ised=%d/%d)", ised, SEDMODEL.NSURFACE );

    // check whether to read from binary or from text files.
    // Text files are slow (1 sec per SED)

    if ( SIMSED_BINARY_INFO.RDFLAG_SED ) {
      // read from binary file
      fread(sedFile, sizeof(sedFile), 1, fpbin1 );
      if ( strcmp(tmpFile,sedFile) != 0 ) {
	printf("\n\n");
	printf("BINARY   SED File: '%s' \n", sedFile );
	printf("EXPECTED SED File: '%s' \n", tmpFile );
	sprintf(c1err,"binary SED file does not match expected file.");
	sprintf(c2err,"Try deleting %s file.", SIMSED_BINARY_FILENAME);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      printf("  Read %s SED surface from binary file : \n", sedcomment);
      fflush(stdout);

      fread( &NSEDBINARY, sizeof(int   ),   1,        fpbin1 );
      fread( SEDBINARY,   sizeof(float ), NSEDBINARY, fpbin1 );
      pack_SEDBINARY(-1);  // transfer SEDBINARY to TEMP_SEDMODEL struct

    } else {      
      // read from text file
      read_SIMSED_flux(tmpFile, sedcomment) ;

      // check day-shift for MJD_EXPLODE or PEAKMAG
      int OPTMASK_EXPLODE = INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE ;
      int OPTMASK_PEAKMAG = (SEDMODEL.OPTMASK & OPTMASK_T0SHIFT_PEAKMAG);
      TEMP_SEDMODEL.TSHIFT = 0.0 ;
      if ( OPTMASK_EXPLODE >= 0 ) { 
	T0shiftExplode_SEDMODEL(OPTMASK_EXPLODE, &TEMP_SEDMODEL, 1); 
      }
      else if ( OPTMASK_PEAKMAG > 0 ) {
	T0shiftPeak_SEDMODEL(&TEMP_SEDMODEL,1); 
      }
    }


    double UVLAM = INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX;
    if ( UVLAM > 0.0 ) { UVLAM_EXTRAPFLUX_SEDMODEL(UVLAM, &TEMP_SEDMODEL); }

    // check array bounds
    if ( TEMP_SEDMODEL.NDAY > SEDMODEL.MXDAY ) {
      sprintf(c1err,"NDAY=%d exceeds expected bound of %d",
	      TEMP_SEDMODEL.NDAY, SEDMODEL.MXDAY ) ;
      sprintf(c2err,"Check DAY bins for SED.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( TEMP_SEDMODEL.NLAM > MXBIN_LAMSED_SEDMODEL ) {
      sprintf(c1err,"NLAM=%d exceeds bound of %d",
	      TEMP_SEDMODEL.NLAM, MXBIN_LAMSED_SEDMODEL ) ;
      sprintf(c2err,"Check LAMBDA bins for SED.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }


    if ( !SIMSED_BINARY_INFO.RDFLAG_FLUX  ) {

      // make fine lambda bins for faster integration
      init_FINEBIN_SEDMODEL(ised); 

      for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
	ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs ;
	init_flux_SEDMODEL(ifilt_obs,ised); 
      }  // ifilt
      init_FINEBIN_SEDMODEL(-1);  // free FINEBIN memory alloc.

    } else {
      // read from binary table
      init_flux_SEDMODEL(0,ised);  // set a few things, then skip integrals
    }

    if ( SIMSED_BINARY_INFO.WRFLAG_SED ) {

      pack_SEDBINARY(+1);   // transfer SEDMODEL to SEDBINARY array
      fwrite( tmpFile,     sizeof(tmpFile), 1,          fpbin1 ) ;
      fwrite( &NSEDBINARY, sizeof(int  ),   1,          fpbin1 ) ;
      fwrite( SEDBINARY,   sizeof(float),  NSEDBINARY, fpbin1 ) ;
    }

    print_ranges_SEDMODEL(&TEMP_SEDMODEL);

  }    //  end loop over ised templates


  fflush(stdout);

  if ( SIMSED_BINARY_INFO.WRFLAG_SED || SIMSED_BINARY_INFO.RDFLAG_SED ) 
    {  fclose(fpbin1);  }


  // write binary integral-flux table to current directory;
  // saves lots of init-time when reading this back

  if ( SIMSED_BINARY_INFO.WRFLAG_FLUX ) {
    IZSIZE = sizeof(REDSHIFT_SEDMODEL) ;
    fwrite(NBIN_SEDMODEL_FLUXTABLE, sizeof(NBIN_SEDMODEL_FLUXTABLE),1,fpbin2);
    fwrite(&IZSIZE, sizeof(IZSIZE),    1, fpbin2); // size of REDSHIFT struct
    fwrite(&REDSHIFT_SEDMODEL, IZSIZE, 1, fpbin2); 

    if ( BINARYFLAG_KCORFILENAME ) 
      { fwrite(SIMSED_KCORFILE,  MXPATHLEN, 1, fpbin2 ); } // Apr 2011

    fwrite(PTR_SEDMODEL_FLUXTABLE, ISIZE_SEDMODEL_FLUXTABLE,1,fpbin2);
    fclose(fpbin2);
    printf("\n  Write filter-integral flux-table to binary file: \n");
    printf("\t %s \n\n", bin2File);
    fflush(stdout);
  }
  else if ( SIMSED_BINARY_INFO.RDFLAG_FLUX ) {
    read_SIMSED_TABBINARY(fpbin2,bin2File);
    fclose(fpbin2);
  }


  // ===========================================
  printf("\n");
  printf(" Global DAY range: %.2f to %.2f \n",
	 SEDMODEL.DAYMIN_ALL, SEDMODEL.DAYMAX_ALL);

  printf("  %s : Done. \n", fnam );
  fflush(stdout);

  return(retval);

} // end of function init_genmag_SIMSED


// ******************************************************
void open_SEDBINARY(char *binFile, bool force_create, 
		    FILE **fpbin, bool *RDFLAG, bool *WRFLAG) {

  // Created July 30 2017
  // Open SED binary in read or write mode.
  // Write/read IVERSION_SIMSED_BINARY.
  //
  // Inputs:  
  //    binFile = name of binary file.
  //    force_create -> force creation of binary of true;
  //                   else create only if file does not exist.
  // Outputs:
  //   fpbin = file pointer
  //   *RDFLAG & WRFLAG to indicate which mode.
  //
  // Dec 14 2021: add force_create arg.
  //
  // ----------- BEGIN -----------

  *RDFLAG = *WRFLAG = false ; *fpbin = NULL;

  checkBinary_SIMSED(binFile); // remove obsolete binary

  *fpbin = fopen(binFile,"rb") ;

  if ( *fpbin == NULL || force_create ) {
    *WRFLAG = true ;
    // re-open in write mode
    *fpbin = fopen(binFile,"wb") ;
    printf("\n Create SED-BINARY file for quicker initialization: \n");
    printf("  %s\n\n", binFile ); fflush(stdout);
    fwrite(&IVERSION_SIMSED_BINARY, sizeof(int *), 1, *fpbin ) ;
  }
  else {
    *RDFLAG = true ;
    printf("\n Read SED-BINARY file for quicker initialization: \n");
    printf("  %s\n\n", binFile );
    IVERSION_SIMSED_BINARY = -9 ;
    fread(&IVERSION_SIMSED_BINARY, sizeof(int *),  1, *fpbin);
  }

  fflush(stdout);

  return ;

} // end open_SEDBINARY


void open_TABBINARY(char *binFile, bool force_create, 
		    FILE **fpbin, bool *RDFLAG, bool *WRFLAG) {

  // Created July 30 2017
  // Open Flux-table binary file in read or write mode.
  // This file depends on survey filters.
  //
  // Inputs:
  //   binFile = name of binary file
  //    force_create -> force creation of binary of true;
  //                   else create only if file does not exist.
  // Outputs:
  //   fpbin = file pointer
  //   *RDFLAG & WRFLAG to indicate which mode.

  // ----------- BEGIN -----------


  *RDFLAG = *WRFLAG = false ; *fpbin = NULL;

  checkBinary_SIMSED(binFile); // remove obsolete binary

  *fpbin = fopen(binFile,"rb") ;

  if ( *fpbin == NULL || force_create ) {
    *WRFLAG = true ;
    *fpbin = fopen(binFile,"wb") ;
  }
  else {
    *RDFLAG = true ;
  }

  fflush(stdout);

  return ;

} // end open_TABBINARY

// ******************************************************
void read_SIMSED_flux(char *sedFile, char *sedComment) {

  // Nov 02, 2011
  // Shell to call rd_sedFlux to read flux and errors from text file.
  // Fill global struct TEMP_SEDMODEL.
  //
  // Aug 12 2017: replace SEDMODEL.FLUX_ERRFLAG argument with OPTMASK
  //
  char fnam[] = "read_SIMSED_flux";

  // -------------- BEGIN -----------------

  printf("# - - - - - - - - - - - - - - - - - - - - - - - - - \n");
  rd_sedFlux(sedFile, sedComment
	     ,Trange_SIMSED, Lrange_SIMSED
	     ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL
	     ,SEDMODEL.OPTMASK
	     ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
	     ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
	     ,TEMP_SEDMODEL.FLUX,  TEMP_SEDMODEL.FLUXERR );  

  int NDAY = TEMP_SEDMODEL.NDAY;
  int NLAM = TEMP_SEDMODEL.NLAM;
  TEMP_SEDMODEL.DAYMIN = TEMP_SEDMODEL.DAY[0];
  TEMP_SEDMODEL.DAYMAX = TEMP_SEDMODEL.DAY[NDAY-1];
  TEMP_SEDMODEL.LAMMIN = TEMP_SEDMODEL.LAM[0];
  TEMP_SEDMODEL.LAMMAX = TEMP_SEDMODEL.LAM[NLAM-1];

} // end of read_SIMSED_flux

// ****************************************************************
void read_SIMSED_TABBINARY(FILE *fp, char *binFile) {

  // Created Dec 16, 2010 by R.Kessler
  // Before reading flux-integral table from binary file,
  // read REDSHIFT information. If table has wider redshift
  // range, then re-allocate memory to use the wider range
  // (intead of aborting due to different NZBIN).
  // This allows creating one binary table with a wide z-range,
  // and using the same table for smaller z-ranges.
  //
  // Mar 24 2021: improve error messaging with CTAG.
  //

  int NERR, idim, IZSIZE_RD, IZSIZE_ACTUAL;
  bool LZSAME, LZOK, LZBAD, LZMIN_OK, LZMAX_OK, LNZBIN_OK ;
  int NBINTMP[NDIM_SEDMODEL_FLUXTABLE+1]  ;

  struct REDSHIFT_SEDMODEL_TYPE  ZTMP ;

  char kcorFile_tmp[MXPATHLEN], CTAG[20] ;
  char fnam[] = "read_SIMSED_TABBINARY" ;

  // ------------ BEGIN ----------

  NERR = 0 ;

  printf("\n  Read filter-integral flux-table from binary file: \n");
  printf("\t %s \n\n", binFile);
  fflush(stdout);

  ZTMP.NZBIN = REDSHIFT_SEDMODEL.NZBIN ;
  ZTMP.ZMIN  = REDSHIFT_SEDMODEL.ZMIN ; // user-requested ZMIN
  ZTMP.ZMAX  = REDSHIFT_SEDMODEL.ZMAX ; // user-requested ZMAX
  IZSIZE_ACTUAL = sizeof(REDSHIFT_SEDMODEL);

  // read header info
  fread( NBINTMP,           sizeof(NBINTMP),    1, fp);
  fread(&IZSIZE_RD,         sizeof(IZSIZE_RD),  1, fp);
  fread(&REDSHIFT_SEDMODEL, IZSIZE_ACTUAL,      1, fp);

  if ( IZSIZE_RD != IZSIZE_ACTUAL ) {
    print_preAbort_banner(fnam);
    printf("\t sizeof(REDSHIFT_SEDMODEL) struct in binary file: %d bytes \n",
	    IZSIZE_RD);
    printf("\t sizeof(REDSHIFT_SEDMODEL) struct now : %d bytes \n",
	    IZSIZE_ACTUAL );

    sprintf(c1err,"%s", "Must remove binary file");
    sprintf(c2err,"%s", binFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
 
  // check if table-redshift range is the same, wider(OK), or smaller(bad)
  LZSAME = ( REDSHIFT_SEDMODEL.ZMIN  == ZTMP.ZMIN &&
	     REDSHIFT_SEDMODEL.ZMAX  == ZTMP.ZMAX &&
	     REDSHIFT_SEDMODEL.NZBIN == ZTMP.NZBIN );

  LZMIN_OK  = ( REDSHIFT_SEDMODEL.ZMIN  <= ZTMP.ZMIN ) ;
  LZMAX_OK  = ( REDSHIFT_SEDMODEL.ZMAX  >= ZTMP.ZMAX ) ;
  LNZBIN_OK = ( REDSHIFT_SEDMODEL.NZBIN >= ZTMP.NZBIN ) ; 
  LZOK      = ( LZMIN_OK && LZMAX_OK && LNZBIN_OK ) ;
  LZBAD     = !LZOK ;

  if ( LZBAD ) {
    NERR++ ;
    sprintf(CTAG,"INFO ");  if ( !LNZBIN_OK ) { sprintf(CTAG,"ERROR"); }
    printf(" %s: NZBIN(request,table) = %d , %d \n", 
	   CTAG, ZTMP.NZBIN,  REDSHIFT_SEDMODEL.NZBIN );

    sprintf(CTAG,"INFO ");  if ( !LZMIN_OK ) { sprintf(CTAG,"ERROR"); }
    printf(" %s: ZMIN(request,table) = %6.4f , %6.4f \n", 
	   CTAG, ZTMP.ZMIN,  REDSHIFT_SEDMODEL.ZMIN );

    sprintf(CTAG,"INFO ");  if ( !LZMAX_OK ) { sprintf(CTAG,"ERROR"); }
    printf(" %s: ZMAX(request,table) = %6.4f , %6.4f \n", 
	   CTAG, ZTMP.ZMAX,  REDSHIFT_SEDMODEL.ZMAX );

    sprintf(c1err,"GENRANGE_REDSHIFT is not compatible with binary table.");
    sprintf(c2err,"Restrict GENRANGE_REDSHIFT or re-make binary table.");
 
    fflush(stdout);
  }

  // if binary table redshift range contains current ZTMP, use larger range
  if ( LZOK  &&  LZSAME == 0 ) {
    printf("  Re-allocate memory with larger redshift range from table. \n");
    fflush(stdout);
    free(PTR_SEDMODEL_FLUXTABLE) ;
    malloc_FLUXTABLE_SEDMODEL ( NFILT_SEDMODEL, REDSHIFT_SEDMODEL.NZBIN,
				NLAMPOW_SEDMODEL, SEDMODEL.MXDAY, 
				SEDMODEL.NSURFACE );
  }


  // check that NBIN matches for each dimension
  if ( !LZBAD ) {
    for ( idim=1; idim <= NDIM_SEDMODEL_FLUXTABLE; idim++ ) {
      if  ( NBINTMP[idim] != NBIN_SEDMODEL_FLUXTABLE[idim] ) {
	printf(" ERROR: NBIN(%s)=%d from binary table, but request is %d \n"
	       , VARNAME_SEDMODEL_FLUXTABLE[idim]
	       , NBINTMP[idim]
	       , NBIN_SEDMODEL_FLUXTABLE[idim] );
	NERR++ ;
	sprintf(c1err,"Binary table mis-match => " ) ;
	sprintf(c2err,"Try deleting %s", binFile );
      }
    }  // idim
  }

  if ( NERR > 0 ) errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 


  // ---------------------------------
  // read name of kcor file form binary, and check for match
  if ( BINARYFLAG_KCORFILENAME ) {
    fread(kcorFile_tmp, MXPATHLEN, 1, fp );
    if ( strcmp(SIMSED_KCORFILE,kcorFile_tmp) != 0 ) {
      sprintf(c1err,"Binary file KCOR_FILE: '%s' ", kcorFile_tmp);
      sprintf(c2err,"but current KCOR_FILE: '%s' ", SIMSED_KCORFILE);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }
  }

  // ------------
  // read entire flux table
  printf("\t Read entire flux table ... "); fflush(stdout);
  fread(PTR_SEDMODEL_FLUXTABLE, ISIZE_SEDMODEL_FLUXTABLE, 1, fp);
  printf("Done reading. \n"); fflush(stdout);

  return ;


} // end of read_SIMSED_TABBINARY



// ****************************************************************
int read_SIMSED_INFO(char *PATHMODEL) {

  // July 2, 2010 RK
  // Read SED.INFO file and store values in SEDMODEL structures
  // Code moved from init_genmag_SIMSED so that this function
  // can be shared with other programs.
  //
  // Jul 16, 2010 RK - read FLUX_ERRFLAG 
  // Nov 02, 2011 RK - read first SED to fill MXDAY
  // Jul 28, 2017 RK - check first and last SED for NDAY in case
  //                   there are more epochs in later SEDs.
  //                   Beware that logic is still fragile.
  //
  // Jul 10 2018 - read option OPT_DAYSHIFT
  // Jul 26 2018 - if UVLAM extrap is set, MINLAMFILT = UVLAM
  // Jul 30 2018 - set SEDMODEL.IPAR_NON1A_INDEX 
  // Aug 26 2018 - read MINSLOPE_EXTRAP_LATE
  // Jan 03 2019 - SEDMODEL.NSURFACE=0 to allow multiple calls.
  //             - Return function arg = NSED (for fortran calls)
  //
  // Apr 28 2019: set Lrange_SIMSED when reading "RESTLAMBDA_RANGE:" key.
  //
  // Jun 6 2022: if PARVAL range > 1E6, write %le format

  char *ptrFile, c_get[80], *ptr_parval, tmpName[60], c_parval[80] ;
  double PARLIM[2], DIF, XN;
  int NPAR, ipar, NSED, NBPAR, ERRFLAG, OPTFLAG ;

  FILE *fp;
  char fnam[] = "read_SIMSED_INFO" ;

  // --------- BEGIN --------

  // open & read info file

  ptrFile = SIMSED_INFO_FILENAME_FULL ;
  sprintf(ptrFile, "%s/%s",PATHMODEL, SIMSED_INFO_FILENAME);
  SEDMODEL.IPAR_NON1A_INDEX = -9;

  if (( fp = fopen(ptrFile,"rt")) == NULL ) {
    sprintf(c1err,"Could not open info file:");
    sprintf(c2err,"%s", ptrFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  printf("\n Read list of SEDs and parameters from \n  %s \n", ptrFile);

  NPAR = NSED = SEDMODEL.NSURFACE = 0 ;

  while( (fscanf(fp, "%s", c_get )) != EOF) {

    if ( strcmp(c_get,"FLUX_SCALE:") == 0 ) 
      { readdouble(fp, 1, &SEDMODEL.FLUXSCALE ); }

    if ( strcmp(c_get,"FLUX_ERRFLAG:") == 0 ) { 
      readint (fp, 1, &ERRFLAG ); 
      if ( ERRFLAG > 0 ) { SEDMODEL.OPTMASK += OPTMASK_FLUXERR_SEDMODEL; }
    }

    if ( strcmp(c_get,"OPTFLAG_T0SHIFT_PEAKMAG:") == 0 ) { 
      readint (fp, 1, &OPTFLAG ); 
      int OVP = ( SEDMODEL.OPTMASK & OPTMASK_T0SHIFT_PEAKMAG);
      // default OPTMASK_T0SHIFT_PEAKMAG = 1; check option to turn it off
      if ( OPTFLAG == 0 ) { SEDMODEL.OPTMASK -= OVP; } 
    }

    if ( strcmp(c_get,"MAG_ERR:") == 0 ) 
      { readdouble(fp, 1, &SEDMODEL.MAGERR_FIX ); }


    // read MINSLOPE_EXTRAPMAG from SED.INFO file only if this
    // parameter was not already set in the sim-input file.
    if ( INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE == 0.0 ) {
      if ( strcmp(c_get,"MINSLOPE_EXTRAPMAG_LATE:") == 0 ) 
	{ readdouble(fp, 1, &INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE ); }
    }

    // read UVLAM_EXTRAP from SED.INFO file only if this
    // parameter was not already set in the sim-input file.
    if ( INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX < 0.0 ) {
      if ( strcmp(c_get,"UVLAM_EXTRAPFLUX:") == 0 ) 
	{ readdouble(fp, 1, &INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX ); }
    }

    if ( strcmp(c_get,"LOGZBIN:") == 0 ) 
      { readdouble(fp, 1, &SEDMODEL.LOGZBIN ); }
    

    if ( strcmp(c_get,"RESTLAMBDA_RANGE:") == 0 ) {
      readdouble(fp, 1, &SEDMODEL.RESTLAMMIN_FILTERCEN );
      readdouble(fp, 1, &SEDMODEL.RESTLAMMAX_FILTERCEN );
      
      Lrange_SIMSED[0] =  SEDMODEL.RESTLAMMIN_FILTERCEN ;
      Lrange_SIMSED[1] =  SEDMODEL.RESTLAMMAX_FILTERCEN ;

      double UVLAM = INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX;
      if (UVLAM > 0.0 ) { SEDMODEL.RESTLAMMIN_FILTERCEN = UVLAM; }

    }

    if ( strcmp(c_get,"NPAR:") == 0 ) 
      {  readint(fp, 1, &SEDMODEL.NPAR ); }

    if ( strcmp(c_get,"PARNAMES:") == 0 ) {
      if ( SEDMODEL.NPAR < 0 ) {
	sprintf(c1err,"PARNAMES key specified before NPAR key");
	sprintf(c2err,"Check %s", ptrFile );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
	readchar(fp, tmpName);
	sprintf(SEDMODEL.PARNAMES[ipar],"%s", tmpName);
	if ( IS_INDEX_SIMSED(tmpName) ) { SEDMODEL.IPAR_NON1A_INDEX=ipar; }
      }
    }

    // ------------

    if ( strcmp(c_get,"SED:") == 0 ) {

      SEDMODEL.NSURFACE++ ;  NSED = SEDMODEL.NSURFACE ;
      readchar(fp, SEDMODEL.FILENAME[NSED] );
      NPAR = SEDMODEL.NPAR ;

      // read parameters into char array to also store them
      // in PARVAL_STRING to preserve format
      for ( ipar = 0; ipar < NPAR; ipar++ ) {
	ptr_parval = SEDMODEL.PARVAL_STRING[NSED][ipar] ;
	readchar(fp, ptr_parval );
	sscanf(ptr_parval, "%le", &SEDMODEL.PARVAL[NSED][ipar] );
      }
    }

  } // end of reading info file


  fclose(fp);

  // for each SIMSED parameter, get NBIN, MIN and MAX,
  // and print summary info.
  
  for ( ipar =0; ipar < NPAR; ipar++ ) {
    fetch_parInfo_SEDMODEL(ipar, tmpName, &NBPAR, PARLIM );
    SEDMODEL.NBIN_PARVAL[ipar] = NBPAR ;
    SEDMODEL.PARVAL_MIN[ipar]  = PARLIM[0] ;
    SEDMODEL.PARVAL_MAX[ipar]  = PARLIM[1] ;

    if ( NBPAR > 1 ) { 
      DIF = PARLIM[1] - PARLIM[0] ;
      XN  = (double)(NBPAR-1) ;
      SEDMODEL.PARVAL_BIN[ipar]  = DIF / XN ;
    }
    else
      { SEDMODEL.PARVAL_BIN[ipar]  = 0.0 ; }


    if ( PARLIM[1] < 1.0E6 ) 
      { sprintf(c_parval,"%8.3f to %8.3f", PARLIM[0], PARLIM[1]); }
    else
      { sprintf(c_parval,"%10.3le to %10.3le", PARLIM[0], PARLIM[1]); }
    printf("    Found '%16s' with %2d bins from %s\n",
	   tmpName, NBPAR, c_parval );

  }

  // -------

  printf("\n Finished reading parameters for %d SEDs \n", NSED );
  fflush(stdout);

  return(SEDMODEL.NSURFACE) ;

} // end of read_SIMSED_INFO


// allow fortran call
int read_simsed_info__(char *PATHMODEL ) { 
  int NSED = read_SIMSED_INFO(PATHMODEL);
  return(NSED);
}

// ===============================
int IS_INDEX_SIMSED(char *parName) {

  int IS_INDEX = 0 ;
  if ( strstr(parName,"INDEX") != NULL ) { IS_INDEX = 1; }
  if ( strstr(parName,"INDX" ) != NULL ) { IS_INDEX = 1; }
  if ( strstr(parName,"index") != NULL ) { IS_INDEX = 1; }
  if ( strstr(parName,"indx" ) != NULL ) { IS_INDEX = 1; }
  return(IS_INDEX);
} // end IS_INDEX_SIMSED

// ===================================
int count_SIMSED_INFO(char *PATHMODEL ) {

  // Apr 28 2019
  // return nubmer of "SED:" keys in SED.INFO file

  int NSED=0;
  FILE *fp;
  char FILENAME[200], c_get[60] ;
  char fnam[] = "count_SIMSED_INFO";

  // ---------------- BEGIN --------------

  sprintf(FILENAME, "%s/%s",PATHMODEL, SIMSED_INFO_FILENAME);

  if (( fp = fopen(FILENAME,"rt")) == NULL ) {
    sprintf(c1err,"Could not open info file:");
    sprintf(c2err,"%s", FILENAME );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }
  printf("\n Read number of SEDs  from \n  %s \n", FILENAME );


  while( (fscanf(fp, "%s", c_get )) != EOF) {
    if ( strcmp(c_get,"SED:")==0 ) { NSED++ ; }
  }
  
  fclose(fp);

  return(NSED);

} // end count_SIMSED_INFO


// **********************************************
void set_SIMSED_MXDAY(char *PATHMODEL, FILE *fpbin, 
		      bool RDFLAG_BINARY, bool WRFLAG_BINARY ) {

  // Created July 30 2017
  // Determine and fill SEDMODEL.MXDAY 
  // If RDFLAG_BINARY = true --> 
  //    read SEDMODEL.MXDAY from fpbin
  // If RDFLAG_BINARY == false --> 
  //   find & read largest ASCII-SED file; then read NDAY
  //                  
  // Aug 10 2017: SEDMODEL.MXDAY = NDAY  and not NDAY+10
  // Dec 29 2017: check sedFile and gzipped file too.
  // Jun 06 2022: MXDAY += 5 in case of very close file sizes.

  int NSED = SEDMODEL.NSURFACE ;
  int ised, istat, size, MXsize, ised_MXsize, NDAY ;
  char sedFile[MXPATHLEN], sedFile_gz[MXPATHLEN], comment[60] ;
  struct stat statbuf ; 
  char fnam[] = "set_SIMSED_MXDAY";

  // ---------------- BEGIN -----------------

  if ( !RDFLAG_BINARY  ) {
    // find and read largest ASCII-SED file
    MXsize = 0 ;  ised_MXsize=-9;
    for(ised=1; ised <= NSED; ised++ ) {
      sprintf(sedFile,    "%s/%s", PATHMODEL, SEDMODEL.FILENAME[ised] );
      sprintf(sedFile_gz, "%s.gz", sedFile);
      istat = stat(sedFile, &statbuf);
      if ( istat != 0 ) { istat = stat(sedFile_gz, &statbuf); }
      size  = statbuf.st_size ;
      if ( size > MXsize ) { MXsize=size; ised_MXsize=ised; }
    }
    
    // read largest file to get MXDAY
    sprintf(sedFile, "%s/%s", PATHMODEL, SEDMODEL.FILENAME[ised_MXsize] );
    sprintf(comment, "NDAY from largest file");
    read_SIMSED_flux(sedFile, comment);  
    NDAY           = TEMP_SEDMODEL.NDAY ;
  }


  /*
  printf("\t NDAY(largest file)=%d  => allocate %d epochs in SEDMODEL \n", 
  	 NDAY, SEDMODEL.MXDAY );
  */

  if ( RDFLAG_BINARY ) {
    fread(&SEDMODEL.MXDAY, sizeof(int*), 1, fpbin);
  }
  else {
    SEDMODEL.MXDAY = NDAY + 5 ; // leave a little slop in file sizes
    if ( WRFLAG_BINARY ) 
      { fwrite(&SEDMODEL.MXDAY, sizeof(int*), 1, fpbin); }
  }

  // allocate DAY array for each SED (Aug 2017)
  int MEM =  SEDMODEL.MXDAY * sizeof(double) ;
  for(ised=1; ised <= NSED; ised++ )  { 
    SEDMODEL.DAY[ised] = (double*) malloc(MEM);  
  }

  return ;

} // end set_SIMSED_MXDAY

// **********************************************
void  set_SIMSED_LOGZBIN(void) {

  // Created July 31 2017
  // Adjust REDSHIFT_SEDMODEL.NZBIN based on SEDMODEL.LOGZBIN.

  double LOGZBIN =  SEDMODEL.LOGZBIN ;
  double ZMIN    =  REDSHIFT_SEDMODEL.ZMIN ;
  double ZMAX    =  REDSHIFT_SEDMODEL.ZMAX ;

  double ZDIF = log10(ZMAX) - log10(ZMIN) ;
  int    NZBIN ;
  // ----------- BEGIN ------------

  NZBIN = (int)( ZDIF/LOGZBIN ) + 1 ;
  init_redshift_SEDMODEL(NZBIN, ZMIN, ZMAX);
  return ;

}  // end set_SIMSED_LOGZBIN

// **********************************************
void dump_SIMSED_INFO(void) {

  printf("   SED FLUX_SCALE = %le \n",    SEDMODEL.FLUXSCALE );
  printf("   SED MAG_ERR    = %5.2f \n",  SEDMODEL.MAGERR_FIX );

  if ( INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE > 0.0 ) {
    printf("   Force late-time extrap >= %.2f mag/day \n",
	   INPUTS_SEDMODEL.MINSLOPE_EXTRAPMAG_LATE );
  }

  // - - - - - 
  printf("   Read Flux errors : " );
  if ( ( SEDMODEL.OPTMASK & OPTMASK_FLUXERR_SEDMODEL) == 0 ) 
    { printf(" NO \n"); }
  else
    { printf(" YES \n"); }

  // - - - - - 
  printf("   Shift T=0 to time of peak  : " );
  if ( (SEDMODEL.OPTMASK & OPTMASK_T0SHIFT_PEAKMAG)==0 ) 
    { printf(" NO \n"); }
  else
    { printf(" YES \n"); }


  /*
  printf("   Shift T=0 to time of explosion : " );
  if ( (SEDMODEL.OPTMASK & OPTMASK_EXPLODE)==0 ) 
    { printf(" NO \n"); }
  else
    { printf(" YES \n"); }
  */


  printf("\n");

  fflush(stdout);

  return ;

}  // end of dump_SIMSED_INFO


// *********************************************
void checkBinary_SIMSED(char *binaryFile) {

  // Mar 29, 2011
  // Remove binaryFile if it has a time-stamp that is earlier than
  // the time-stamp of the SED.INFO file.
  // This check prevents updating the SIMSED model, and then using
  // a stale binary-table file of fluxes.
  //
  // July 2017: replace ABORT with warning and remove
  //            stale binary file. 
  //
  // Dec 15 2020: 
  //  + if ISBATCH_SIMSED and stale binary, then abort.
  //  + check files in loop.

  int tdif_sec, ifile ;
  double tdif_day ;
  char file_type[2][20] = { "SED.INFO file", "KCOR_FILE" } ;
  char *ptrFile, *ptrType ;
  char rm[400];
  char line[] = "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-";
  char fnam[] = "checkBinary_SIMSED" ;
  FILE *fp;
  
  // ------------ BEGIN --------------

  // if binaryFile does not exist, return.
  if (( fp = fopen(binaryFile,"r")) == NULL ) { return ; }
  fclose(fp);

  for ( ifile=0; ifile < 2; ifile++ ) {

    ptrType = file_type[ifile]; // SED.INFO, KCOR

    if ( ifile == 0 ) 
      { ptrFile = SIMSED_INFO_FILENAME_FULL; }
    else {
      ptrFile = SIMSED_KCORFILE ; 
      if ( BINARYFLAG_KCORFILENAME == 0 ) { continue; }
    }

    tdif_sec = file_timeDif(binaryFile, ptrFile );

    // xxx    if ( ifile ==1 ) { tdif_sec = -54.0 ; } // xxx REMOVE

    if ( tdif_sec < 0 ) {
      tdif_day = -(double)tdif_sec / 86400. ;

      printf("\n%s\n WARNING INFO: \n", line);
      printf(" %s: \n   %s\n",  ptrType, ptrFile);
      printf(" Stale binary file: \n   %s \n", binaryFile );
      printf(" is %.3f days older than %s", tdif_day, ptrType );

      if ( ISBATCH_SIMSED ) {
	sprintf(c1err,"Cannot remake binary file in batch mode.");
	sprintf(c2err,"Must run sim interactively to re-make binary.");
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }
      
      sprintf(c1err,"Removing stale binary file above,");
      sprintf(c2err,"and re-making new binary file.");
      errmsg(SEV_WARN, 0, fnam, c1err, c2err ); 
      printf("%s\n\n", line);
      
      sprintf(rm,"rm %s", binaryFile);
      system(rm); return ;
    }
  } // end ifile loop

  return ;

} // end of checkBinary_SIMSED



// ****************************************************************
void genmag_SIMSED(
		  int OPTMASK       // (I) bit-mask of options
		  ,int ifilt_obs    // (I) absolute filter index
		  ,double x0        // (I) SIMSED flux-scale
		  ,int  NLUMIPAR    // (I) Number of lumi-params 
		  ,int    *iflagpar // (I) parameter flags
		  ,int    *iparmap  // (I) ipar index map
		  ,double *lumipar  // (I/O) shape/lumi parameter(s)
		  ,double RV_host   // (I) RV of host
		  ,double AV_host   // (I) AV of host
		  ,double mwebv     // (I) Galactic extinction: E(B-V)
		  ,double z         // (I) Supernova redshift
		  ,int    Nobs      // (I) number of epochs
		  ,double *Tobs_list   // (I) list of obs times (since mB max) 
		  ,double *magobs_list  // (O) observed mag values
		  ,double *magerr_list  // (O) model mag errors
		  ,int *index_sed       // (O) SED index (if defined)
		  ) {


  /****
  Return observer frame mag in absolute filter index "ifilt_obs" 
  for input SIMSED parameters.

   OPTMASK+=1 (bit0) => return flux instead of mag ; magerr still in mag.
                       (to avoid discontinuity for negative flux)

   OPTMASK+=2 (bit1) => print warning message when model flux < 0

   OPTMASK+=8 (bit3) => dump flag

        HISTORY


  Mar 1 2017: add RV_host & AV_host arguments

  Apri 30 2018: 
    + Sinterp=0 of Trest is outside SEDMODEL.MINDAY_ALL-SEDMODEL.MAXDAY_ALL

  May 7 2018:  
   + for undefined flux, magerr_list[epobs] = MAGERR_UNDEFINED ;
   + extrapolate for Trest > MAXDAY

  July 20 2018:  OPTMASK & 8 --> debug DUMP flag

  Jul 30 2018: add output arg *index_sed

  ***/

  double  meanlam_obs, meanlam_rest, ZP, z1, Tobs, Trest, flux, arg, Sinterp  ;
  int ifilt, epobs, OPT_COLORLAW    ;
  int  LDMP_BADFLUX, LDMP_DEBUG, LRETURN_MAG, LRETURN_FLUX, LSEDSEQ ;
  double AV, XT_MW, XT_HOST ;
  double magobs, magerr, tmpPar;
  char *cfilt ;

  char fnam[] = "genmag_SIMSED" ;

  // ----------------- BEGIN -----------------


  // set default optopn and then parse OPTMASK

  LDMP_BADFLUX = LRETURN_FLUX = LDMP_DEBUG = 0 ;
  LRETURN_MAG  = 1;
  LSEDSEQ      = 0;

  if ( (OPTMASK & 1) > 0 ) { LRETURN_MAG    = 0; LRETURN_FLUX = 1; }
  if ( (OPTMASK & 2) > 0 ) { LDMP_BADFLUX   = 1; }
  if ( (OPTMASK & 4) > 0 ) { LSEDSEQ        = 1; }
  if ( (OPTMASK & 8) > 0 ) { LDMP_DEBUG     = 1; }
  
  // make sure that user NLUMIPAR matches expected number 
  // of parameters.

  if ( NLUMIPAR != SEDMODEL.NPAR ) {
    sprintf(c1err, "You passed NLUMIPAR=%d", NLUMIPAR);
    sprintf(c2err, "but expected SEDMODEL.NPAR=%d", SEDMODEL.NPAR );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }


  // translate absolute filter index into sparse index
  ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
  z1    = 1. + z;

  // filter info for this "ifilt"
  meanlam_obs  = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda  
  meanlam_rest = meanlam_obs/z1;
  ZP           = FILTER_SEDMODEL[ifilt].ZP ;
  cfilt        = FILTER_SEDMODEL[ifilt].name ;

  // get approx Galactic extinction using central wavelength of filter
  AV   = RV_MWDUST * mwebv ;
  OPT_COLORLAW = MWXT_SEDMODEL.OPT_COLORLAW ;
  XT_MW = GALextinct(RV_MWDUST,AV,meanlam_obs,OPT_COLORLAW); 

  // get approx extinction from host in rest-frame (Mar 1 2017)
  if ( AV_host > 1.0E-9 ) 
    { XT_HOST = GALextinct ( RV_host, AV_host, meanlam_rest, 94 ); }
  else
    { XT_HOST = 0.0 ; }

  // - - - - - - -  - - 
  //determine integer times which sandwich the times in Tobs

  if ( LDMP_DEBUG ) {
    printf(" xxx ------- %s DEBUG DUMP ------------ \n", fnam );
    printf(" xxx Nobs=%d  ifilt_obs=%d z=%.3f  DAYMAX_ALL=%.1f  lumipar=%.1f\n",
	   Nobs, ifilt_obs, z, SEDMODEL.DAYMAX_ALL, *lumipar );
    fflush(stdout);
  }
  
  for ( epobs=0; epobs < Nobs; epobs++ ) {
    
    Tobs  = Tobs_list[epobs];
    Trest = Tobs / z1;

    if ( LSEDSEQ ) {
      // get flux from next SED on the grid; fill all of the lumipar
      Sinterp = nextgrid_flux_SIMSED(iflagpar, iparmap, lumipar,
				     ifilt_obs, z, Trest);
    }
       
    else {
      // interpolate flux for these lumipar. Note that
      // early/late extrap is included in underlying get_flux_SEDMODEL.
      Sinterp = interp_flux_SIMSED(iflagpar,  iparmap,  lumipar,
				   ifilt_obs, z, Trest);
    }

    
    // - - - - - - - - - -
    if ( LDMP_DEBUG ) {
      printf(" xxx Trest=%8.2f  x0=%9.3le   Sinterp=%9.3le  "
	     " ISED=%d  DAYMAX=%.1f \n",
	     Trest, x0,  Sinterp,
	     ISED_SEDMODEL, SEDMODEL.DAYMAX[ISED_SEDMODEL]  );
      fflush(stdout);
      
    }
    // - - - - - - - - - -

    flux = x0 * Sinterp ;  

    if ( Sinterp == FLUX_UNDEFINED || isnan(flux) ) {
      magobs = MAG_UNDEFINED ;
      magerr = MAGERR_UNDEFINED ;
    }
    else if ( flux <= 1.0E-30 ) {
      if ( LDMP_BADFLUX > 0 ) {
	printf("  genmag_SIMSED Warning:");
	printf(" Flux(%s)<0 at Trest = %6.2f => return mag=99 \n",
	     cfilt, Trest );
      }
      magobs = MAG_ZEROFLUX ;
      magerr = MAGERR_UNDEFINED ;
    }
    else{
      magobs = (ZP + XT_MW + XT_HOST) - 2.5*log10(flux);
      if ( magobs > MAG_ZEROFLUX ) { magobs = MAG_ZEROFLUX; }
      magerr = SEDMODEL.MAGERR_FIX ;
    }


    magobs_list[epobs] = magobs ;
    magerr_list[epobs] = magerr ;
          
    if ( fabs(magobs) < -400.  && ifilt_obs==3) {
      printf(" XXX ------------------------------- \n");
      printf(" xxx %s  OPTMASK=%d  LRETURN_FLUX=%d \n", 
	     fnam, OPTMASK,LRETURN_FLUX);
      printf(" xxx %s  Tobs=%6.3f  x0=%8.3le  Sint=%8.3le  "
	     "flux=%8.3le mag=%.2f \n",
	     cfilt, Tobs, x0, Sinterp, flux, magobs_list[epobs] );
      
      printf(" xxx \t iflagpar=%d,%d   iparmap=%d,%d    lumipar=%f,%f\n",
	     iflagpar[0], iflagpar[1], iparmap[0], iparmap[1],
	     lumipar[0], lumipar[1] );
      //      debugexit(fnam); // xxxxxxxxxx
    }
    
    
    // check option to return flux intead of mag;
    // preserves continutity for negative model fluxes.

    if ( LRETURN_FLUX ) {
      arg     = -0.4 * magobs_list[epobs] ;
      magobs_list[epobs] = pow(TEN,arg);
    }


  } // end epobs loop over epochs


  // if any SIMSED parameter looks or smells like an index,
  // return this value as *index_sed.
  int IPAR = SEDMODEL.IPAR_NON1A_INDEX;
  if ( IPAR >= 0 ) {
    tmpPar = SEDMODEL.PARVAL[ISED_SEDMODEL][IPAR] + 0.01 ;
    *index_sed = (int)(tmpPar);
  }
  else
    { *index_sed = -9; }

} // end of genmag_SIMSED


// ****************************************************
double nextgrid_flux_SIMSED (
			  int *iflag,        // (I) flag params
			  int *iparmap,      // (I) ipar index map
			  double *lumipar,   // (I) SIMSED lumipars
			  int ifilt_obs,     // (I) obs filter
			  double z,          // (I) redshift
			  double Trest       // (I) Trest (days)
			  ) 
{

  // Increment ISIMSED_SEQUENTIAL and return integrated flux for this SED.
  // Also return the *lumipar array.
  //
  // Mar 6 2017: set global ISED_SEDMODEL
  //
  // Dec 20 2018: ISIMSED_SEQUENTIAL is set from snlc_sim, not passed
  //              as lumipar argument.

  double Sinterp ;
  int  ipar, ISED ;

  // ----------- BEGIN -----------

  // get SED index for GRIDONLY-sequential

  ISED = ISIMSED_SEQUENTIAL ;

  // extract SED index for GRIDONLY-sequential option
  if ( ISED  > SEDMODEL.NSURFACE ) { 
    sprintf(c1err,"ISED =%d exceeds number of SEDS(%d)",
	    ISED , SEDMODEL.NSURFACE);
  }


  Sinterp = get_flux_SEDMODEL( ISED, 0, ifilt_obs, z, Trest);

  ISED_SEDMODEL = ISED; // set globa, Mar 6 2017

  // load *lumipar array
  for ( ipar=0; ipar < SEDMODEL.NPAR ; ipar++ ) {
    lumipar[ipar] = SEDMODEL.PARVAL[ISED][ipar];
  }


  return(Sinterp) ;

} // end of nextgrid_flux_SIMSED


// ****************************************************
double interp_flux_SIMSED(
			  int *iflag,        // (I) flag params
			  int *iparmap,      // (I) ipar index map
			  double *lumipar,   // (I) SIMSED lumipars
			  int ifilt_obs,     // (I) obs filter
			  double z,          // (I) redshift
			  double Trest       // (I) Trest (days)
)
{

  /* ------------------------------------------------
    Created May 03, 2010 by B. Diemer
    Interpolate in multi-dimensional space of SIMSED parameters
    to get flux-integral.

    Jun 24, 2010 (RK) - add ilampow=0 argument to interp_flux_SEDMODEL.
                        No effect for SIMSED models, but the extra arg
                        is needed for the SALT2  model.

    Jun 29, 2010 (RK) 
        - replace case test on iflag with bit-mask test
        - fix tabs

    May 16, 2011 (RK) skip interpolation if there is just 1 SED.

    Nov 18, 2011 (RK) pass *iparmap so that user can pick SIMSED params
                      in any order.  See ipar_model below.

   Nov 28,  2011:  bugfix: fill *lumipar array if there is just one SED.

   Aug 17 2015: if all params are GRIDONLY, then find ISED at grid-node
                instead of interpolating . See NGRIDONLY & NMATCH.

   Mar 6 2017: set global ISED_SEDMODEL

   Dec 26 2018:
     fix bug from v10_63g (July 2018). For GRIDONLY option, 
     make sure to load  *lumipar.

   Jun 3 2022: 
     + fix index bug computing range
     + fix index bug loading *lumipar ... before it worked only if
       selected model params were in same order is in SED.INFO file

  -------------------------------------------------- */

  /*
    Verbose can have values greater than 1; those create a lot of
    debug output.
  */
  int verbose = 0 ;
  int pars[INTERP_SIMSED_MAX_DIM];
  int pars_baggage[INTERP_SIMSED_MAX_BAGGAGE_PARS];
  int i, j, k, ISED, num_dims, num_pars_baggage;
  int found_corner, index, ipar_model, NPAR, ipar, ipar_user ;
  int flag, NGRIDONLY, NMATCH=0 ;

  double Sinterp, left_min_diff, right_min_diff;
  double diff, diff0, diff1, parval, range, term;

  char fnam[] = "interp_flux_SIMSED";

  // --------------- BEGIN --------------

  //  if ( fabsf(Trest) < 1.0  ) { verbose = 1; } // xxxxxxxxxxxxxx

  num_dims = 0;
  num_pars_baggage = 0;
  ISED_SEDMODEL = -9;
  
  if ( verbose ) { 
    printf("\n DDDDDDDDump %s  at Trest=%6.3f DDDDDDDDDD \n", fnam, Trest ); 
  }

  // if just one SED surface, then there is no need to 
  // interpolate in the SED-parameter space;
  // just interpolate in the space of redshift and Trest.
  if ( SEDMODEL.NSURFACE == 1 ) {
    ISED = 1;
    ISED_SEDMODEL = ISED; // set globa, Mar 6 2017
    Sinterp = get_flux_SEDMODEL( ISED, 0, ifilt_obs, z, Trest);

    // load *lumipar array
    for ( ipar=0; ipar < SEDMODEL.NPAR ; ipar++ ) 
      { lumipar[ipar] = SEDMODEL.PARVAL[ISED][ipar];  }

    return(Sinterp) ;
  }


  NGRIDONLY = 0 ;

  for (i = 0; i < SEDMODEL.NPAR; i++) {

    flag       = iflag[i];
    ipar_model = iparmap[i];
    pars_baggage[i] = 0 ;
		 
    if ( flag & OPTMASK_GEN_SIMSED_GRIDONLY )
      { NGRIDONLY++ ; } // Aug 17 2015 RK

    if ( (flag & OPTMASK_GEN_SIMSED_PARAM ) > 0 )  {
      pars[num_dims] = i;

      if(verbose) {
	printf("\t => Found par=%d, name=%s, value = %f (num_dims=%d) \n", 
	       i, SEDMODEL.PARNAMES[ipar_model], lumipar[pars[i]], num_dims);
      }
      num_dims++;

    }
    else if ( (flag & OPTMASK_GEN_SIMSED_param ) > 0 ) {
      pars_baggage[i] = 1 ;
      num_pars_baggage++;
    }
    else {
      // need to abort here
    }
  } // i  loop


  // Aug 17 2015
  // if all params are GRID-ONLY, just find SED and skip interpolation
  // to avoid pathologies.
  
  NPAR = SEDMODEL.NPAR ;
  if( NGRIDONLY == NPAR-num_pars_baggage ) {
    for ( ISED = 1; ISED <= SEDMODEL.NSURFACE ; ISED++ ) {
      
      NMATCH = 0 ;
      for(j=0; j < NPAR; j++ ) {
	if ( pars_baggage[j] ) { continue ; } // skip baggage	     
	// xxx mark delete range  = SEDMODEL.PARVAL_MAX[j] - SEDMODEL.PARVAL_MIN[j] ;
	ipar_model = iparmap[j];
        range      = SEDMODEL.PARVAL_MAX[ipar_model] - SEDMODEL.PARVAL_MIN[ipar_model] ;
	parval     = SEDMODEL.PARVAL[ISED][ipar_model];
	diff       = (parval - lumipar[j]) / range ;
	
	/*
	printf(" xxx ISED=%d ipar_model=%d  parval=%6.2f "
	       "lumipar[%d]=%6.2f  range=%.2f  diff=%.2f\n",
	       ISED, ipar_model, parval, j,lumipar[j], range, diff ); 
	*/
	
	if ( fabs(diff) < 0.0001 ) { NMATCH++ ; }
      }
      if ( NMATCH == NGRIDONLY ) { 
	Sinterp = get_flux_SEDMODEL(ISED, 0, ifilt_obs, z, Trest);	
	ISED_SEDMODEL = ISED; // set globa, Mar 6 2017

	// load *lumipar array
	for ( ipar=0; ipar < SEDMODEL.NPAR ; ipar++ ) {
	  ipar_model = iparmap[ipar];
	  lumipar[ipar] = SEDMODEL.PARVAL[ISED][ipar_model];
	  // xxx mark delete lumipar[ipar] = SEDMODEL.PARVAL[ISED][ipar];  
	}

	return(Sinterp) ;
      }
    } // end ISED loop


    // if we get here abort on error.
    sprintf(c1err, "Could not find GRIDONLY match for %d params", NGRIDONLY);
    sprintf(c2err, "z=%f Trest=%f ifilt_obs=%d  NMATCH=%d",
	    z, Trest, ifilt_obs, NMATCH );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);

  } // end of NGRIDONLY


  if(verbose)
    {
      printf("\t => Found %d interp params and %d baggage params.\n", 
	     num_dims, num_pars_baggage);
    }
  
  if (num_dims > INTERP_SIMSED_MAX_DIM) {
      sprintf(c1err, "Cannot interpolate in %d dimensions.", num_dims);
      sprintf(c2err, " ");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  /*
   * Define arrays which depend upon the number of dimensions
   */
  int sheet_size[num_dims][2];
  double hypercube[num_dims][2];
  int bits[num_dims];
  double pars_baggage_interp[num_pars_baggage];
  double pars_baggage_term[num_pars_baggage];

  /*
   Find indexes bracketing values; for each dimension, 
   there is an (n-1)dim sheet of data points below and 
   a similar sheet above the value we wish to interpolate, 
   if the values lie on the required uniform grid.   
   The sheet should have at least 2^(n - 1) data points in it, 
   otherwise interpolation becomes impossible (this usually 
   happens when the data point lies outside the allowed region).
   */

  int num_corners = (int)pow(2.0, (double)num_dims);
  int corners[num_corners];
  
  for (i = 0; i < num_corners; i++ )
    {   corners[i] = INTERP_SIMSED_INVALID_CORNER;    }
  
  for (i = 0; i < num_dims; i++)
    {
      sheet_size[i][0] = 0;
      sheet_size[i][1] = 0;
      left_min_diff = -INTERP_SIMSED_START_DIFF;
      right_min_diff = INTERP_SIMSED_START_DIFF;

      ipar_model = iparmap[i];

      /*
       * Find out how many values there are in each sheet closest to the
       * point of interpolation in both directions.
       */
      for (j = 0; j < SEDMODEL.NSURFACE; j++)
	{

	  // BD orig   diff = SEDMODEL.PARVAL[j + 1][pars[i] + 1] - lumipar[pars[i]];
	  diff = SEDMODEL.PARVAL[j + 1][ipar_model] - lumipar[pars[i]]; // RK
	  
	  if (diff <= 0)
	    {
	      if (fabs(diff - left_min_diff) < INTERP_SIMSED_DELTA)
		{
		  sheet_size[i][0]++;
		} else if (diff > left_min_diff)
		{
		  left_min_diff = diff;
		  sheet_size[i][0] = 1;
		}
	    } else
	    {
	      if (fabs(diff - right_min_diff) < INTERP_SIMSED_DELTA)
		{
		  sheet_size[i][1]++;
		} else if (diff < right_min_diff)
		{
		  right_min_diff = diff;
		  sheet_size[i][1] = 1;
		}
	    }
	}
            
      /*
       * The hypercube stores the coordinates of the sheets in all dimensions.
       */

      hypercube[i][0] = lumipar[pars[i]] + left_min_diff;
      hypercube[i][1] = lumipar[pars[i]] + right_min_diff;

      if (verbose)	{
	printf("Dim %d has %d data points in left sheet.\n", 
	       i, sheet_size[i][0]);
	printf("Dim %d has %d data points in right sheet.\n", 
	       i, sheet_size[i][1]);
	printf("Hypercube[%d] (left) = %f.\n",  i, hypercube[i][0]);
	printf("Hypercube[%d] (right) = %f.\n", i, hypercube[i][1]);
      }

      if (sheet_size[i][0] < pow(2, num_dims - 1))
	{	  	  
	  sprintf(c1err, "Dimension %d does not have enough data points in lower sheet (%d).", i, sheet_size[i][0]);
	  sprintf(c2err, "%s = %f",
		  SEDMODEL.PARNAMES[ipar_model], lumipar[pars[i]] );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	}
      
      if (sheet_size[i][1] < pow(2, num_dims - 1))
	{
	  sprintf(c1err, "Dimension %d does not have enough data points in higher sheet (%d).", i, sheet_size[i][1]);
	  sprintf(c2err, "%s = %f (sheet_size=%d,%d  num_dims=%d)",
		  SEDMODEL.PARNAMES[ipar_model], lumipar[pars[i]],
		  sheet_size[i][0], sheet_size[i][1],  num_dims );
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	}


    }

  /*
   * Find corner points from known coordinates of the hypercube: a corner point will
   * have coordinates lying on the hypercube in all dimensions.
   */
  for (j = 0; j < SEDMODEL.NSURFACE; j++)
    {
      found_corner = 1;
      
      for (i = 0; i < num_dims; i++)
	{
	  ipar_model = iparmap[i];
	  diff0 = fabs(SEDMODEL.PARVAL[j + 1][ipar_model] - hypercube[i][0]);
	  diff1 = fabs(SEDMODEL.PARVAL[j + 1][ipar_model] - hypercube[i][1]);

	  if ( diff0 < INTERP_SIMSED_DELTA)
	    {
	      if(verbose > 1)
		{
		  printf("Successful comparison between %f and %f, delta %f.\n", 
			 SEDMODEL.PARVAL[j + 1][pars[i] + 1], 
			 hypercube[i][0], SEDMODEL.PARVAL[j + 1][ipar_model] - hypercube[i][0]);
		}
	      bits[i] = 0;
	    } else if ( diff1 < INTERP_SIMSED_DELTA)
	    {
	      if(verbose > 1)
		{
		  printf("Successful comparison between %f and %f, delta %f.\n", 
			 SEDMODEL.PARVAL[j + 1][pars[i] + 1], 
			 hypercube[i][1], SEDMODEL.PARVAL[j + 1][ipar_model] - hypercube[i][1]);
		}
	      bits[i] = 1;
	    } else
	    {
	      if(verbose > 1)
		{
		  printf("Unsuccessful comparison between %f, %f and %f.\n", 
			 SEDMODEL.PARVAL[j + 1][ipar_model], 
			 hypercube[i][0], hypercube[i][1]);
		}
	      found_corner = 0;
	      continue;
	    }
	}
      
      /*
        Store the found corner in the appropriate bin. 
	The binning is done using binary number, 
	e.g. in 3D it would be low-low-low, high-low-low, low-high-low  etc.
       */

      if (found_corner)
	{
	  
	  index = 0;
	  for (k = 0; k < num_dims; k++)
	    {
	      index = index + bits[k] * dual_bits_SIMSED[k];
	    }
	  if (corners[index] != INTERP_SIMSED_INVALID_CORNER)
	    {
	      sprintf(c1err, "Corners[%d] already filled; two data points have the same coordinates.", index);
	      sprintf(c2err, " ");
	      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	    } else
	    {
	      corners[index] = j;
	      if (verbose)
		{
		  printf("Found index %d, j = %d.\n", index, j);
		}
	    }
	}
      
      // TODO abort if all corners found? Would avoid check on double filled corners
    }
  
  /*
   * Check whether all corners have been filled. If not, there are data points missing
   * from the grid.
   */
  for (i = 0; i < num_corners; i++)
    {
      if (corners[i] == INTERP_SIMSED_INVALID_CORNER)
	{
	  sprintf(c1err, "Corners[%d] could not be found.", i);
	  sprintf(c2err, " ");
	  errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
	}
    }

  /*
   * Interpolation; each corner represents a distance-weighted term.
   */
  Sinterp = 0;
  
  for(k = 0; k < num_pars_baggage; k++)
    {
      pars_baggage_interp[k] = 0;
    }

  for (i = 0; i < num_corners; i++)
    {
      /*
       * Set 0/1 bit for each dimension, corresponding to the lower and upper
       * index on the cube in coordinate space.
       */
      for (j = 0; j < num_dims; j++)
	{
	  bits[j] = ((i & dual_bits_SIMSED[j]) != 0);
	}
      
      /*
       * Set term to data point in question; 
       * multiply term by distance weightings for each
       * dimension.
       */
      term = get_flux_SEDMODEL(corners[i] + 1, 0, ifilt_obs, z, Trest);
      ISED_SEDMODEL = corners[0]+1;
      
      for(k = 0; k < num_pars_baggage; k++)
	{
	  ipar_model = iparmap[num_dims+k]; // RK
	  pars_baggage_term[k] = 
	    SEDMODEL.PARVAL[corners[i] + 1][ipar_model]; // RK
	}
      
      for (j = 0; j < num_dims; j++)
	{
	  term = term * fabs(lumipar[pars[j]] - hypercube[j][1 - bits[j]]);
	  
	  for(k = 0; k < num_pars_baggage; k++)
	    {
	      pars_baggage_term[k] = pars_baggage_term[k] * fabs(lumipar[pars[j]] - hypercube[j][1 - bits[j]]);
	    }
	}
      
      Sinterp = Sinterp + term;

      for(k = 0; k < num_pars_baggage; k++)
	{
	  pars_baggage_interp[k] = pars_baggage_interp[k] + pars_baggage_term[k];
	}
    }
  
  /*
   * Divide by the volume of the hypercube in coordinate space.
   */
  for (i = 0; i < num_dims; i++)
    {
      Sinterp = Sinterp / (hypercube[i][1] - hypercube[i][0]);

      for(k = 0; k < num_pars_baggage; k++)
	{
	  pars_baggage_interp[k] = pars_baggage_interp[k] / (hypercube[i][1] - hypercube[i][0]);
	}
    }


  // fill baggage parameters in lumipar array

  for(k = 0; k < num_pars_baggage; k++) {
    ipar_user = num_dims + k ;
    lumipar[ipar_user] = pars_baggage_interp[k]; // RK
  }


  return(Sinterp) ;
  
} // end of interp_flux_SIMSED


// ****************************************************
double interp1D_flux_SIMSED(
			    int *iflag,        // (I) flag params
			    double *lumipar,  // (I) SIMSED lumipars
			    int ifilt_obs,    // (I) obs filter
			    double z,         // (I) redshift
			    double Trest      // (I) Trest (days)
			    ) {

  /*
   Created Dec 21, 2009 by R.Kessler
   When NPAR=1, interpolate in 1D space of SIMSED parameter
   to get flux-integral. This is just for initial testing,
   and to check the limiting case of the more general
   interp_flux_SIMSED function.

   Feb 19, 2010: allow for LUMIPAR grid to be in any order,
                 rather than assume monotonically increasing order. 
  */

  int  ipar, IPAR, I0SED, I1SED, NPAR, NPAR_RETURN, istat, ilampow ;
  double parval, parval0, parval1, LUMIPAR, frac, S0int, S1int, Sinterp  ;
  char fnam[] = "interp1D_flux_SIMSED";

  // ------------------ BEGIN ----------------

  Sinterp = 0.0;

  NPAR = IPAR = NPAR_RETURN = 0;
  LUMIPAR = -999.;
  for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {
    if ( ( iflag[ipar] & 1) > 0 ) { 
      NPAR++ ; 
      IPAR = ipar; 
      LUMIPAR = lumipar[ipar];
    } 
    else NPAR_RETURN++;
  }

  if ( NPAR > 1 ) { return Sinterp ; }

  istat = get_SEDMODEL_INDICES( IPAR, LUMIPAR, &I0SED, &I1SED ); 

  if ( istat != SUCCESS ) {
    sprintf(c1err,"Invalid I0SED=%d  I1SED=%d for LUMIPAR=%f", 
	    I0SED, I1SED, LUMIPAR ) ;
    sprintf(c2err,"iflag=%d  ifilt_obs=%d  z=%le  Trest=%6.2f",
	    *iflag, ifilt_obs, z, Trest );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );  
  }


  parval0 = SEDMODEL.PARVAL[I0SED][IPAR];
  parval1 = SEDMODEL.PARVAL[I1SED][IPAR];
  frac    = (LUMIPAR - parval0)/(parval1-parval0);

  ilampow = 0;

  S0int = get_flux_SEDMODEL(I0SED, ilampow, ifilt_obs, z, Trest );
  S1int = get_flux_SEDMODEL(I1SED, ilampow, ifilt_obs, z, Trest );
  
  Sinterp  = S0int + (S1int-S0int)*frac;

  if ( NPAR_RETURN == 0 ) return Sinterp ;

  // compute interpolated values for return (baggage) parameters

  for ( ipar=0; ipar < SEDMODEL.NPAR; ipar++ ) {

    if ( ( iflag[ipar-1] & 2) > 0 ) { 
      parval0 = SEDMODEL.PARVAL[I0SED][ipar] ;
      parval1 = SEDMODEL.PARVAL[I1SED][ipar] ;
      parval  = parval0 + (parval1-parval0)*frac;
      lumipar[ipar] = parval ;
    }
  }


  if ( ifilt_obs == -2 && fabs(Trest) < 2.0 && LUMIPAR > 160.) {
    printf(" xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
    printf(" xxx lumipar = %f \n", LUMIPAR );
    printf(" xxx I0SED=%d  I1SED=%d \n", I0SED, I1SED);
    printf(" xxx parval[0,1] = %f %f   frac=%f \n", parval0, parval1, frac );
    printf(" xxx S0int=%le   S1int=%le   Sint=%le \n", 
	   S0int, S1int, Sinterp);
    debugexit("genmag_SIMSED Sint");
  }

  return Sinterp ;

} // interp1D_flux_SIMSED



