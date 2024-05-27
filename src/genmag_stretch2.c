/***************************************
 Created Aug 19, 2008 by R.Kessler

 double-stretch model where str1 & str2 are the
 stretches for the rise & fall, respectively.
 Single-stretch model is obtained by passing 
 args with str2 = str1.

 Model is rest-frame only, filter system is UBVRI.


 Nov 10, 2010: change init_genmag_stretch2() to input model_rootdir
               instead of template filename.
  
****************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sntools.h"
#include "genmag_stretch2.h"

//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>


// =======================================================
// define mangled function to be called by fortran

int init_genmag_stretch2__ ( char *version, char *cfilt_rest ) {
  int istat;
  istat = init_genmag_stretch2(version,cfilt_rest);
  return istat;
}


int genmag_stretch2__ 
(
  double  *str1         // (I) stretch
 ,double  *str2         // (I) stretch
 ,int     *ifilt_rest   // (I) absolute filter index
 ,int     *nepoch       // (I) Number of epochs to compute mag
 ,double  *Trest        // (I) list of T-Tpeak, rest-frame days
 ,double  *genmag       // (O) list of rest mags
 ,double  *genmag_err   // (O) list of rest mags errors
 ) 
{

  int istat;


  istat = genmag_stretch2(*str1, *str2, *ifilt_rest, *nepoch, 
			  Trest, genmag, genmag_err );
  return istat;

}
  


// =============== Regular C functions ==========================



// **********************************************************************
int init_genmag_stretch2 ( 
			  char *VERSION    // (I) template version
			  ,char *cfilt_rest   // (O) rest-frame filters
			  ) {

  /*****
  Read stretch template from
   1. official area ($SNDATA_ROOT/models/stretch)
   2. non1a area    ($SNDATA_ROOT/snsed/non1a)
   3. current directory

   Return string of rest-frame filters; 
   i.e., "UBVRI", "abcde", etc ...

   Nov 10, 2010: change first argument to modelPath and hard-wire
                 name of template file in *modelPath.

   Jan 18, 2011: change modelPath arg to version , and construct
                 STRETCH2_MODELPATH locally

   Feb 22 2017: allow input *VERSION to include path

  ****/

  int iep, ifilt, ifiltdef, NWD    ;
  char  c_get[40], cfilt[2], version[60] ;
  double Trest, *ptrmag;
  FILE* fp ;

  char  fnam[] = "init_genmag_stretch2" ;

  // ------- BEGIN ------


  init_stretch2_LEGACY();

  printf("\n init_genmag_stretch2: INITIALIZE DOUBLE-STRETCH MODEL \n");

  sprintf(FILTERSTRING,"%s", FILTERSTRING_DEFAULT );

  extract_MODELNAME(VERSION,
		    STRETCH2_MODELPATH, version); // returned

  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    sprintf( STRETCH2_MODELPATH, "%s/%s", 
	     getenv(PRIVATE_MODELPATH_NAME), version );
  }
  else if ( strlen(STRETCH2_MODELPATH) > 0 ) {
    // do nothing for user-specified path/model
  }
  else {
    sprintf( STRETCH2_MODELPATH, "%s/models/stretch2/%s", 
             getenv("SNDATA_ROOT"), version );
  }

  printf("  Read stretch2 model parameters from \n\t  %s\n",
         STRETCH2_MODELPATH );




  sprintf(STRTEMPL_FILE, "%s/%s", 
	  STRETCH2_MODELPATH, "stretch_UBVRI_template.dat" );
  if ( (fp = fopen(STRTEMPL_FILE, "rt")) == NULL ) {
    sprintf(c1err,"Could not find stretch-template file:");
    sprintf(c2err,"%s", STRTEMPL_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  // -------------------

  printf("   Read stretch-template file from: \n   %s \n", STRTEMPL_FILE);


  // init a few things ...

  STRTEMPL_FILTERS[0] = 0 ;
  NFILT_STRTEMPL = NWD = 0;

  // set default mag-errors ... need to allow user input somehow ??
  for ( ifilt=0; ifilt < MXFILT_STRTEMPL ; ifilt++ )
    STRTEMPL_MAGERR[ifilt] = 0.1 ;


  // ===================
  // first try to read header; if nothing found in first ten words,
  // then assume it's a legacy file and hard-wird old defaults.

  while( (fscanf(fp, "%s", c_get)) != EOF) {

    NWD++ ;

    if ( strcmp(c_get,"FILTER:") == 0 ) {
      readchar(fp,cfilt);
      NFILT_STRTEMPL++ ;
      sprintf(STRTEMPL_FILTERS,"%s%s", STRTEMPL_FILTERS, cfilt);
    }


    if ( strcmp(c_get,"SNTYPE:") == 0 )
      readchar(fp,STRTEMPL_SNTYPE);

    if ( strcmp(c_get,"EPOCH:") == 0 ) { 
      printf("\t Found filters '%s' for SN type %s \n", 
	     STRTEMPL_FILTERS, STRTEMPL_SNTYPE );
      goto RDEPOCH ;
    }


    if ( NWD >= 50 && strlen(STRTEMPL_FILTERS) == 0 ) {
      NFILT_STRTEMPL = NFILT_STRTEMPL_LEGACY ;
      sprintf(STRTEMPL_FILTERS,"%s", STRTEMPL_FILTERS_LEGACY );
      printf("\t No header => assume legacy file (<= v7_07) with filters %s\n",
	     STRTEMPL_FILTERS );
      LEGACYFLAG = 1 ;
      for ( ifilt=0; ifilt < NFILT_STRTEMPL; ifilt++ )
	STRTEMPL_MAGERR[ifilt] = STRTEMPL_MAGERR_LEGACY[ifilt] ;

      goto RDEPOCH ;
    }

  }  // end of fscanf loop


 RDEPOCH:

  rewind(fp);

  iep = 0;

  if ( LEGACYFLAG > 0 ) {

    while( fscanf(fp, "%lf", &Trest ) > 0 ) {
      iep++ ;
      STRTEMPL_TREST[iep] = Trest;
      ptrmag = STRTEMPL_MAGREST[iep];
      readdouble(fp, NFILT_STRTEMPL, ptrmag);
    }


  } else {

    while( (fscanf(fp, "%s", c_get)) != EOF) {
      if ( strcmp(c_get,"EPOCH:") == 0 ) {
	iep++ ;
	ptrmag = STRTEMPL_MAGREST[iep];
	readdouble(fp, 1, &Trest);
	readdouble(fp, NFILT_STRTEMPL, ptrmag);
	STRTEMPL_TREST[iep] = Trest;
      }

    }

  } // end of LEGACYFLAG if-block


  if ( iep <= 0 ) {
    sprintf(c1err,"%d epochs in stretch template.", iep);
    sprintf(c2err,"Check template file listed above.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  NEPOCH_STRTEMPL = iep ;
  printf("\t Read %d epochs from template \n", NEPOCH_STRTEMPL );


  sprintf(cfilt_rest, "%s", STRTEMPL_FILTERS ); // load return arg

  // convert filter string into integer indices
  for ( ifilt=0; ifilt < NFILT_STRTEMPL; ifilt++ ) {

    sprintf(cfilt,"%c",  STRTEMPL_FILTERS[ifilt] ) ;
    ifiltdef = INTFILTER(cfilt);
    IFILTDEF_STRTEMPL[ifilt] = ifiltdef ;

    printf("\t Filter %s : INDEX = %2.2d  MAGERR = %5.2f \n",
	   cfilt, ifiltdef, STRTEMPL_MAGERR[ifilt] );
  }


  // get min,max Trest

  TMINDEF =  9999. ;   TMAXDEF = -9999. ;
  for( iep=1; iep <= NEPOCH_STRTEMPL; iep++ ) {
    Trest = STRTEMPL_TREST[iep];
    if ( Trest < TMINDEF ) TMINDEF = Trest;
    if ( Trest > TMAXDEF ) TMAXDEF = Trest;
  } // end of iep loop


  fclose(fp);

  /**
  printf("\n xxxxx GRACEFUL EXIT from init_genmag_stretch xxxxx \n");
  exit(1) ; //  xxxxxxx
  */

  return SUCCESS ;
  
} // end init_genmag_stretch2



void init_stretch2_LEGACY(void) {

  // Moved here, Feb 21, 2013

  LEGACYFLAG = 0; 
  sprintf(STRTEMPL_FILTERS_LEGACY,"UBVRI");
  STRTEMPL_MAGERR_LEGACY[0] = 2.0 ;  // U
  STRTEMPL_MAGERR_LEGACY[1] = 0.1 ;  // B
  STRTEMPL_MAGERR_LEGACY[2] = 0.1 ;  // V
  STRTEMPL_MAGERR_LEGACY[3] = 0.1 ;
  STRTEMPL_MAGERR_LEGACY[4] = 0.5 ;  // I
} 

// *************************************************************
int genmag_stretch2 (
		    double  str1         // (I) stretch
		    ,double  str2         // (I) stretch
		    ,int     ifilt_rest   // (I) absolute filter index
		    ,int     nepoch       // (I) Nepochs to compute mag
		    ,double  *Trest_list  // (I) list of T-Tpeak, rest days
		    ,double  *genmag_list // (O) list of rest mags
		    ,double  *genmag_err_list // (O) list of errors
		    ){

  char fnam[] = "genmag_stretch2" ;
  int iep, j1, j2, ifilt, i ;

  double Trest, Tstr, mag, w1, w2, tempmag1, tempmag2 ;
  double magoff_lumi ;

  // -------------------- START -------------------

  // translate absolute "ifilt_rest" into sparse "ifilt"
  ifilt = -9;
  for ( i = 0; i < NFILT_STRTEMPL; i++ ) {
    if ( IFILTDEF_STRTEMPL[i] == ifilt_rest ) ifilt = i;
  }


  if ( ifilt < 0 ) {
    sprintf(c1err,"Invalid ifilt_rest=%d (%c). Check filter-defs in", 
	    ifilt_rest, FILTERSTRING[ifilt_rest]);
    sprintf(c2err, "%s", STRTEMPL_FILE );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }


  for ( iep = 0; iep < nepoch; iep++ ) {

    Trest = *(Trest_list + iep);

    // determine stretch-corrected time in rest-frame
    if ( Trest < 0.0 ) 
      Tstr = Trest / str1 ;
    else
      Tstr = Trest / str2 ;


    if ( Tstr < TMINDEF ) {
      mag = STRTEMPL_MAGREST[1][ifilt] ;
      goto LOADMAG;
    }

    if ( Tstr > TMAXDEF ) {
      mag = STRTEMPL_MAGREST[NEPOCH_STRTEMPL][ifilt] ;
      goto LOADMAG;
    }

    // get j1,j2 date-indices that bracket Tstr

    JDATES2 ( Tstr, &j1, &j2 ) ; 

    w1 = fabs(STRTEMPL_TREST[j2] - Tstr );
    w2 = fabs(STRTEMPL_TREST[j1] - Tstr );

    tempmag1  = STRTEMPL_MAGREST[j1][ifilt] ;
    tempmag2  = STRTEMPL_MAGREST[j2][ifilt] ;
 
    mag = (w1*tempmag1 + w2*tempmag2) / (w1 + w2);

  LOADMAG:

    magoff_lumi = -0.78 * (str2-1.0); // added 6.11.2009

    *(genmag_list + iep)     = mag + magoff_lumi ;
    *(genmag_err_list + iep) = STRTEMPL_MAGERR[ifilt] ;

    /**
    printf(" %d-%d xxx Trest = %6.1f, Tstr = %6.1f, mag[%d] = %f \n",
	   j1,j2, Trest, Tstr, ifilt, mag );
    **/

  }


  return SUCCESS ;
}


// ************************************************
void JDATES2 ( double Trest, int *j1, int *j2 ) {

  // return integer indices j1 and j2 such that
  // template_dates[j1,j2] bracket Trest

  int j;

  // ------------- BEGIN -----------

  *j1 = 0 ;  // init
  *j2 = 1 ;

  for(j=1; j <= NEPOCH_STRTEMPL ; j++ ) {
    if( STRTEMPL_TREST[j] > Trest ){
      *j1 = j-1 ;
      *j2 = j ;
      break ;
    }
  }

  if ( *j1 < 0 ) { *j1=0; *j1=1 ; }

}  // end of JDATES
