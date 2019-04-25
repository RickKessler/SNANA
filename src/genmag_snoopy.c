/********** 
  dm15temp.c:  driver program for GLoEs (Gaussial Local Estimation) to 
  smooth the 2-d surface of time-dm15-flux  lightcurve data.  
 
  Author:  Chris Burns (cburns@ociw.edu) 
  The problem:  you've got a set of low-z lightcurves in certain filters.  
  You have analyzed each lightcurve to determine the time of maximum, 
  the maximum magnitude, and the decline rate (dm15).  Now, each 
  observation can be placed on `a grid (x -> epoch, y-> dm15) and the 
  flux (magnitude) defines a surface on this grid.  We now wish to 
  sample this surface (defined by heterogenously spaced points) at any 
  arbitrary (epoch,dm15) point.  Typically, this is so that you can 
  fit another lightcurve (high-z or without pre-max data) using least 
  squares.
  The solution:  Use GLoEs (Gaussian Local Estimation) to smooth the  
  2D surface and interpolate to any given point you wish.  It does 
  this by fitting a 2nd degree polynomial centered at the point of 
  interest.  The data points are weighted by both their internal 
  errors as well as a sliding 2D elliptical Gaussian of variable 
  width (hence the "Gaussian" in the name).  In this way, only points 
  close to the interpolating point are used (hence the "local" in the 
  name).  See gloes.c for more details about the algorithm.
 
  This program expects to find a file named "templates.dat" in which 
  you have listed the locations of the template lightcurve data.  
  The format of templates.dat should be:
    SN-name ; filter ; redshift ; dm15 ; T(Bmax) ; Mmax ; filename
  (ie, semi-colon separated).  SN-name really isn't used, it's just 
  for your own sanity.  Here's an example:
 
  SN2004gf ; B ; 0.01 ; 1.1 ; 452.5 ; 16.5 ; SN2004gf_B.dat
  SN2004gf ; V ; 0.01 ; 1.1 ; 452.5 ; 15.9 ; SN2004gf_V.dat
  SN2004gf ; g ; 0.01 ; 1.1 ; 452.5 ; 15.8 ; SN2004gf_g_prime.dat
  SN2004gf ; r ; 0.01 ; 1.1 ; 452.5 ; 15.4 ; SN2004gf_r_prime.dat
  ....
 
  Then, each datafile listed in this file has the following format:
  time magnitude magnitude-error
 
  The time should be in days.  Tmax listed in templates.dat will be 
  subtracted to get the epoch.  The Mmax from templates.dat will be 
  subtracted from the magnitudes so that Mmax = 0 for each template 
  (or M(Tmax) = 0 if you choose Mmax accordingly).
  
  This program will read in templates.dat and then the files listed 
  therein, filling in the arrays needed to define the surface.  It 
  reads the value of dm15 from the command line and then samples the 
  surface on the interval [-10,70] in one-day increments.  The 
  interpolated values and errors are printed to STDOUT. 
 
  The width of the window function is a function of position on the grid.
  Currently, a fixed formula is used (see sigmax() and sigmay()), since
  automatic algorithms in 2D are non-trivial.  These formulas are purely
  trial and error:  they seem to work well.  As data are added, they will
  likely need to be updated to improve efficiency (or experiment for
  yourself).  But generically, one can say that there is very little
  curvature in the dm15-direction (hence I use a fixed width), and high
  curvature near max in epoch-direction, which decreases away from max,
  hence I use a linearly increasing window.
 
  You need the GSL libraries and include files.  Compile using:
 
  gcc -o dm15temp dm15temp.c gloes.c -lm -lgsl -lgslcblas -I{Idir} -L{Ldir}
 
  where {Idir} and {Ldir} should be replaced with the location of 
  the gsl include files and libraries, respectfully, if they are 
  not in the standard search path.  Then, simply call the program like:
 
  dm15temp dm15 dm15 dm15 > data
 
  You can list any number of values for dm15, they will be output
  sequencially to STDOUT.  The format is described in the output 
  header.  Note that the lightcurves are in flux units and are 
  normalized to 1.0 at the maximum of the _individual lightcurves_.  
  So, at t(Bmax), the SN will not, in general, have 0 color.  
  It will, however, have all pseudo-colors at max (e.g., Bmax-Vmax) = 0
  so it's up to you to impose any colors yourself 
  (e.g., Phillips et al, 1999).

  
  @@@@@ Installation modifications from R.Kessler (Aug 2010) @@@@

  - integrate init_genmag_snoopy and genmag_snoopy with the dm15 
    functions so that they are all in one file.

  - screen dump info during the init stage

  - Original code had ifilt = 0-8 in the init, but 1-9 later.
    Modified to use 1-9 scheme everywhere to avoid confusion.
    Also define parameters IFILT_SNOOPY_B[VugriYJH] to avoid 
    hard-wired filter indices.

  - put global variables inside global structure 
    SNOOPY_TEMPLATE[ifilt]



       HISTORY
    ~~~~~~~~~~~~~
 Nov 1, 2010 RK - leave out blank space in *filtlist return arg from 
                  init_genmag_snoopy().

 Jan 18, 2011 RK - change modelPath arg to *version (for code uniformity),
                   and allow PRIVATE_MODELPATH_NAME override

 Jul 27, 2011 RK 
   - read EXCLUDE_RANGE keys from SNOOPY.INFO file
   - extrapolate early times to have continuous model; see Textrap

 Aug 12, 2011 RK: read ERROR map from text file in snoopy area.
                  See ERRORFILE key in the SNooPy.INFO file.

 Feb 13, 2012 RK - add get_LAMRANGE_snoopy() to define rest-frame
                    wavelength range.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include "gloes.h"

#define  SNGRIDREAD
#include "fitsio.h"
#include "sntools.h"
#include "sntools_grid.h" 

#include "genmag_snoopy.h" 


// ============ MANGLED FORTRAN FUNCTIONS =============
// (Nov 1,2010, R.Kessler)

int init_genmag_snoopy__( char *version, int *optmask, double *snoopyArgs, 
			  char *filtlist) {
  int istat;
  istat = init_genmag_snoopy ( version, *optmask, snoopyArgs, filtlist) ;
  return istat;
}


int genmag_snoopy__(int *ifilt, double *dm15, int *nobs, 
		  double *rest_dates, double *rest_mags, double *rest_magerrs ) {
  int istat;
  istat = genmag_snoopy(*ifilt, *dm15, *nobs, 
			rest_dates, rest_mags, rest_magerrs ) ;
  return istat;
}

void get_lamrange_snoopy__(double *lammin, double *lammax) {
  get_LAMRANGE_snoopy(lammin,lammax);
}


// external spline function

extern void in2dex_(int *ispline, int *N2D,
		    double *XX, double *YY, double *ZZ,
		    double *XLIM, double *YLIM, double *SS, int *IERR );

extern double ge2dex_ ( int *IND, double *Trest, double *Lrest, int *IERR ) ;


// =========== CODE ================

/* Initiallize the SNooPy template generator 
 *    Inputs:  version  (where templates.dat and templates/ folder reside)
 *             optmask   bitmask of options
 *                       bit1 => return relative flux instead of mag
 *             snoopyArgs (NULL pointer.  We don't need any arguments)
 *    Outputs: filtlist:  "BVugriYJH"
 *    return value:  -3:  memory allocation error
 *                   -2:  templates.dat or one of templates/ files not found
 *                   -1:  A path name was too long.  Need to increase array
 *                        sizes in the code
 *                    0:  success
 *                    1:  Bad parsing of the templates.dat file
 *  
 * Feb 22 2017: allow input *VERSION to include model path
 */
int init_genmag_snoopy( char *VERSION, int optmask, 
			double *snoopyArgs, char *filtlist) {
  /* Initialization function for the SNooPy templates */
  int ifilt, res, opt_fits_read ;
  char cfilt[2], msg[100], tmpFile[400], version[60];
  double mag;

  // ----------- BEGIN -------------

  // hard-wire a few strings
  sprintf(SNOOPYNAME,"SNooPy");
  sprintf(FILTLIST_SNOOPY, "BVugriYJH" );

  res = 0;
  OPT_RELFLUX_SNOOPY = (optmask & 1); 

  if ( OPT_RELFLUX_SNOOPY == 1 ) 
    { sprintf(msg, "Initialize %s to return relative flux", SNOOPYNAME); }
  else
    { sprintf(msg, "Initialize %s to return mag", SNOOPYNAME); }


  print_banner(msg);
  printf("\t (Based on Burns et al., arXiv:1010:4040) \n\n" );


  // construct SNOOPY_MODELPATH
  extract_MODELNAME(VERSION,
		    SNOOPY_MODELPATH, version); // returned

  if ( getenv(PRIVATE_MODELPATH_NAME) != NULL ) {
    // ENV
    sprintf( SNOOPY_MODELPATH, "%s/%s", 
	     getenv(PRIVATE_MODELPATH_NAME), version );
  }
  else if ( strlen(SNOOPY_MODELPATH) > 0 ) {
    // do nothing for user-defined path/model
  }
  else {
    // default location under $SNDATA_ROOT
    sprintf( SNOOPY_MODELPATH, "%s/models/snoopy/%s", 
             getenv("SNDATA_ROOT"), version );
  }

  printf("  Read SNooPy model parameters from \n\t  %s\n",
         SNOOPY_MODELPATH );

  // always read templates to get filters
  res = load_data_snoopy(SNOOPY_MODELPATH);

  // model info must be read after load_data_snoopy so that we have filters
  parse_snoopy_modelinfo(SNOOPY_MODELPATH);

  // read snoopy model errors
  read_snoopy_errors(SNOOPY_MODELPATH);

  // fill filter-string output arg and check that each peak mag is filled.
  filtlist[0] = 0 ;
  for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) {
    
    sprintf(cfilt, "%s", SNOOPY_TEMPLATE[ifilt].FILTERNAME );  
    sprintf(filtlist,"%s%s", filtlist, cfilt );

    mag = SNOOPY_MODELINFO.PEAKMAG_REST[ifilt] ;
    if ( mag < -900.0 ) {
      sprintf(c1err,"Undefined %s  PEAKMAG_REST[%s]", SNOOPYNAME, cfilt );
      sprintf(c2err,"   " );
      MADABORT(c1err,c2err);
    }
  }

#ifdef SNGRIDREAD
  if ( SNOOPY_MODELINFO.OPT_GRIDFILE == 1 ) {

    sprintf(msg,"Ignore templates above and use pre-calculated snoopy grid.");
    printf("\n  ==> %s \n", msg);

    sprintf(tmpFile,"%s/%s", SNOOPY_MODELPATH, SNOOPY_MODELINFO.GRIDFILE ) ;
    opt_fits_read = 1;  // bit1 => verbose 
    fits_read_SNGRID(opt_fits_read,tmpFile, &SNGRID_SNOOPY );


    if ( strcmp(filtlist,SNGRID_SNOOPY.FILTERS) != 0 ) {
      sprintf(c1err,"Expected filters = '%s' from snoopy templates.", 
	      filtlist );
      sprintf(c2err,"But found filters '%s' in snoopy grid.", 
	      SNGRID_SNOOPY.FILTERS );
      MADABORT(c1err,c2err);
    }


    /*
    printf("\n");
    for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) 
      { init_spline_snoopy(ifilt); }

    debugexit("init_spline");
    */

  }
#endif


  return(res);

} // end of init_genmag_snoopy



/* Generate SNooPy templates and return them
 *
 *  Inputs:   ifilt      (1 <= ifilt <= 9)  the filter you want.
 *            dm15       The decline rate parameter
 *            nobs       Number of rest-frame epochs
 *            rest_dates The rest-frame epochs at which to generate mags.
 *  Outputs:  rest_mags  Normalized rest-frame magnitudes:  m=0 at peak
 *            rest_magerrs mag errors (argument added by R.K)
 *  return value:  -2:  memory allocation error
 *                 -1:  invalid filter
 *                  0:  success 
 *                  5:  internal warning:  might need to check output.  
*/
int genmag_snoopy(int ifilt, double dm15, int nobs, 
		  double *rest_dates, double *rest_mags, double *rest_magerrs ) {

   /* These are the default window function paramters.  */
   double 
     sigx0    = 3.0
     ,scalex  = 0.1
     ,maxsigx = 10.0
     ,sigy0   = 0.3
     ,*emag
     ,relFlux[500]         // max relFlux = 1
     ,relFlux_exact[500]   // used for testing only
     ,relFlux_err[500]
     ,t, tmp
     ,nullflux = 1.0E-50
     ,NULLFLUX = 2.0E-50
     ,magtmp, magerr
     ,Textrap[2]
     ,dTextrap = 1.0 
     ,relFlux_ref[2], relFlux_slope
     ;

   int normalize = 1;
   int res, i, LTEST_INTERP, LDMP ;
   char cfilt[2];

   // ----------- BEGIN -------------

   LTEST_INTERP = 0 * SNOOPY_MODELINFO.OPT_GRIDFILE  ;

   sprintf(cfilt, "%s", SNOOPY_TEMPLATE[ifilt].FILTERNAME );  // for messages

   if (! (emag = (double *) malloc(nobs*sizeof(double))) ) {
      return(-2);
   }
   /* check inputs */
   if (ifilt < IFILT_SNOOPY_MIN || ifilt > IFILT_SNOOPY_MAX ) {
     sprintf(c1err,"Invalid ifilt=%d passed to genmag_snoopy", ifilt );
     sprintf(c2err,"Valid ifilt range is %d - %d ", 
	     IFILT_SNOOPY_MIN,  IFILT_SNOOPY_MAX ) ;
     MADABORT(c1err,c2err);
   }

   for  ( i=0; i < nobs; i++ ) { 
     relFlux[i]      = nullflux; 
     relFlux_err[i]  = nullflux; 
   }


   res = -9.0 ;

   if ( SNOOPY_MODELINFO.OPT_GRIDFILE == 1 ) {
#ifdef SNGRIDREAD
     gridinterp_snoopy(ifilt, dm15, nobs, rest_dates, 
		       relFlux, relFlux_err );  // <== returned

     if  ( LTEST_INTERP ) {
       res = dm15temp(ifilt, dm15, rest_dates, nobs, relFlux_exact, emag, 
		      normalize, sigx0, scalex, maxsigx, sigy0);
     }
#endif
   } 
   else {
     res = dm15temp(ifilt, dm15, rest_dates, nobs, relFlux, emag, 
		    normalize, sigx0, scalex, maxsigx, sigy0);

     // return array of relFlux_err
     dm15temp_errors(ifilt, dm15, rest_dates, nobs, relFlux_err);
   }

   // Now we need to apply some masks, and convert relFlux -> mag

   if ( LTEST_INTERP ) 
     { printf("  xxx Test interp for DM15 = %f \n", dm15 ); }

   for (i = 0 ; i < nobs ; ++i ) {
     
      if ( LTEST_INTERP ) {
	tmp = relFlux[i] / relFlux_exact[i] ;
	printf("   xxx %s-relFlux @ t = %6.2f : interp/exact=%6.4f/%6.4f = %6.4f \n",
	       cfilt, t, relFlux[i], relFlux_exact[i],tmp );
      }

      // global cuts on the grid

     t = rest_dates[i] ;

     if ( dm15 < SNOOPY_MODELINFO.DM15_RANGE[0] ) 
       { relFlux[i] = relFlux_err[i] = NULLFLUX ; }
     if ( dm15 > SNOOPY_MODELINFO.DM15_RANGE[1] )
       { relFlux[i] = relFlux_err[i] = NULLFLUX ; }

     // avoid mag of negative  flux
     if ( relFlux[i] < nullflux )
       { relFlux[i] = relFlux_err[i] = NULLFLUX ; }

     // Filter-specific restrictions 
     if ( relFlux[i] > 2*NULLFLUX ) 
       {  Textrap[0] = Textrap_snoopy(ifilt,dm15,t); }
     else
       {  Textrap[1] = 0.0 ; }

     /* 
     if ( abs(t) < 3.0 ) {  // DDDDDDDDDDDDDD
       printf(" xxx ifilt=%d  dm15=%5.2f  Trest=%6.1f : relFlux=%6.3f  Textrap = %f %f \n", 
	      ifilt, dm15, t, relFlux[i], Textrap[0], Textrap[1] );
     }
     */

     // extrapolate when reading model from GRID
     // extrapolate model from Textrap to 't'
     if ( Textrap[0] < 0.0 && SNOOPY_MODELINFO.OPT_GRIDFILE > 0 ) {
       LDMP = 0 ;
       Textrap[1] = Textrap[0] + dTextrap ;
       gridinterp_snoopy(ifilt, dm15, 2, Textrap, relFlux_ref, relFlux_err );
       relFlux_slope = (relFlux[1] - relFlux[0])/dTextrap ;
       relFlux[i] = 
	 modelflux_extrap(t, Textrap[0], 
			  relFlux_ref[0], relFlux_slope, LDMP );
     }
     
     
     // convert relFlux to mag
     magtmp = -2.5 * log10(relFlux[i])  
       + SNOOPY_MODELINFO.PEAKMAG_REST[ifilt] 
       + SNOOPY_MODELINFO.DM15SLOPE[ifilt] * ( dm15 - DM15REF_SNOOPY );


     // conver relFlux_err to mag-error
     //  rest_magerrs[i] = 2.5*log10(1.0+relFlux_err[i]/relFlux[i]) ;
     rest_magerrs[i] = relFlux_err[i] / relFlux[i]  ;


     if ( OPT_RELFLUX_SNOOPY ) {
       magerr = rest_magerrs[i];
       rest_mags[i]    = relFlux[i] ;
       rest_magerrs[i] = relFlux[i] * magerr ;
     }
     else {
       rest_mags[i] = magtmp ;
     }
     
     /*
      if ( fabs(t) < 1.0 ) {
	printf(" xxxx rest_mag[%s,T=%6.2f, dm15=%5.2f]= %6.3f (relFlux=%6.4f) \n", 
	     cfilt, rest_dates[i], rest_mags[i], dm15, relFlux[i] );
      }
      */

   }

   return(res);

} // end of genmag_snoopy



// =============================
double Textrap_snoopy(int ifilt, double dm15, double Trest ) {


  // Created July 27, 2011
  // If input 'ifilt, dm15, Trest' is excluded by the model,
  // return min/max value of valid Trest for which to 
  // extrapolate. If input values are OK then return 0.0

  int i, NRANGE;
  double Tret, TT;

  // ----------------- BEGIN --------------------
  Tret = 0.0 ;

  NRANGE = SNOOPY_MODELINFO.NRANGE_EXCLUDE ;
  if ( NRANGE <= 0 ) { return Tret ; }

  TT = SNOOPY_MODELINFO.TREST_RANGE[0] ;
  if ( Trest < TT ) { return TT; }

  TT = SNOOPY_MODELINFO.TREST_RANGE[1] ;
  if ( Trest < TT ) { return TT; }


  for ( i=1; i <= NRANGE ; i++ ) {
    if ( ifilt != SNOOPY_MODELINFO.EXCLUDE_IFILT[i]   ) { continue ; }
    if ( dm15  < SNOOPY_MODELINFO.EXCLUDE_DM15[i][0]  ) { continue ; }
    if ( dm15  > SNOOPY_MODELINFO.EXCLUDE_DM15[i][1]  ) { continue ; }
    if ( Trest < SNOOPY_MODELINFO.EXCLUDE_TREST[i][0] ) { continue ; }
    if ( Trest > SNOOPY_MODELINFO.EXCLUDE_TREST[i][1] ) { continue ; }
    
    TT = SNOOPY_MODELINFO.EXCLUDE_TREST[i][1] ;
    if ( Trest < -1.0 && Trest < TT ) { return TT ; }

    TT = SNOOPY_MODELINFO.EXCLUDE_TREST[i][0] ;
    if ( Trest > +1.0 && Trest > TT ) { return TT ; }
  }

  return Tret;

} // end of Textrap_snoopy



/* function for returning the width of the window function in epoch 
   as a function of epoch and dm15.  Right now, just a linearly 
   increasing window with a cap 
*/
double  sigmax(double x, double y) {
   double sigma = sigx0 + xscale*fabs(x);
   if (sigma > sigxmax) sigma = sigxmax;
   return(sigma);
}

/* function for returning the width of the window function in dm15 
   as a function of epoch and dm15.  Right now, just a constant value 
*/
double sigmay(double x, double y) {
   double sigma = sigy0;
   return(sigma);
}

int load_data_snoopy(char *path) {
   FILE *fp1, *fp2;
   int i, j, nchars, NTMP, nptmp ;

   char 
     filename[1024]
     ,filename2[256]
     ,line[256]
     ,SN[80]
     ,tmpname[MXFILT_SNOOPY][MXEP_SNOOPY][80]
     ,filter[15]
     ;

   //     char fnam[] = "load_data_snoopy" ;

   double dm15, Tmax, Mmax, zed;
   double epoch, mag, emag, arg, flux, Trest, wtmp ;

   // ----------- BEGIN -----------

   nchars = snprintf(filename, 1024, "%s/%s", path, "templates.dat");
   if (nchars > 1024) {
     sprintf(c1err,"Path name for %s templates.dat is too long", SNOOPYNAME );
     sprintf(c2err,"Char-len of path is %d ", nchars ) ;
     MADABORT(c1err,c2err);
   }

   printf("  Read SNooPy templates from: \n");
   printf("  %s \n", filename );

   if (! (fp1 = fopen(filename, "r"))) { 
     sprintf(c1err,"Could not find file with list of %s templates:", 
	     SNOOPYNAME );
     sprintf(c2err," '%s'", filename ) ;
     MADABORT(c1err,c2err);
   }

   /* Initialize arrays */
   for (i = 0 ; i < MXFILT_SNOOPY ; ++i) {

     SNOOPY_TEMPLATE[i].NSN = 0 ;
     SNOOPY_TEMPLATE[i].NEPOCH = 0 ;

     for(j = 0 ; j < MXEP_SNOOPY ; ++j) {
       SNOOPY_TEMPLATE[i].x[j]  = 0.0 ;
       SNOOPY_TEMPLATE[i].y[j]  = 0.0 ;
       SNOOPY_TEMPLATE[i].z[j]  = 0.0 ;
       SNOOPY_TEMPLATE[i].wz[j] = 0.0 ;
     }
   } // end of i-loop over SNoopy filters


   /* Read in template information */
   while(fgets(line, 256, fp1) != NULL) {
      if(line[0] == '#') continue;
      line[255] = '\0';
      strcpy(SN, strtok(line, ";"));
      sscanf(strtok(NULL, ";"), "%s", filter);
      zed = atof(strtok(NULL, ";"));
      dm15 = atof(strtok(NULL, ";"));
      Tmax = atof(strtok(NULL, ";"));
      Mmax = atof(strtok(NULL, ";"));

      sscanf(strtok(NULL, ";"), "%s", filename2);
      nchars = snprintf(filename, 1024, "%s/%s", path, filename2);
      if (nchars > 1024) {
	sprintf(c1err,"Path name for %s template is too long", SNOOPYNAME );
	sprintf(c2err,"Char-len of path is %d ", nchars ) ;
	MADABORT(c1err,c2err);
      }

      if ( ! ( fp2 = fopen(filename, "r") ) ) {
	sprintf(c1err,"Could not find %s template file", SNOOPYNAME );
	sprintf(c2err," '%s'", filename ) ;
	MADABORT(c1err,c2err);
      }

      if (strcmp(filter, "B") == 0) {
	i = IFILT_SNOOPY_B ;
      } else if (strcmp(filter, "V") == 0) {
	i = IFILT_SNOOPY_V ;
      } else if (strcmp(filter, "u") == 0) {
	i = IFILT_SNOOPY_u ;
      } else if (strcmp(filter, "g") == 0) {
	i = IFILT_SNOOPY_g ;
      } else if (strcmp(filter, "r") == 0) {
	i = IFILT_SNOOPY_r ;
      } else if (strcmp(filter, "i") == 0) {
	i = IFILT_SNOOPY_i ;
      } else if (strcmp(filter, "Y") == 0) {
	i = IFILT_SNOOPY_Y ;
      } else if (strcmp(filter, "J") == 0) {
	i = IFILT_SNOOPY_J ;
      } else if (strcmp(filter, "H") == 0) {
	i = IFILT_SNOOPY_H ;
      } else if (strcmp(filter, "K") == 0) {
	i = IFILT_SNOOPY_K ;
      } else {
	sprintf(c1err,"%s filter = '%s' is unknwown.", SNOOPYNAME, filter );
	sprintf(c2err,"  " ) ;
	MADABORT(c1err,c2err);
      }

      SNOOPY_TEMPLATE[i].NSN++ ;
      sprintf(SNOOPY_TEMPLATE[i].FILTERNAME, "%s", filter );

      while( fscanf(fp2, "%lf %lf %lf", &epoch, &mag, &emag) != EOF) {
	nptmp = SNOOPY_TEMPLATE[i].NEPOCH ;

	if ( nptmp >= MXEP_SNOOPY ) {
	  sprintf(c1err,"%s nepoch exceeds array bound of %d", 
		  SNOOPYNAME, MXEP_SNOOPY );
	  sprintf(c2err,"Check parameter MXEP_SNOOPY" ) ;
	  MADABORT(c1err,c2err);
	}

	arg   =  -0.4*(mag - Mmax) ;
	flux  = pow(10.0,arg);
	Trest = (epoch - Tmax)/(1.0+zed);
	wtmp  = 1.0857/(flux*emag) ;      // weight

	SNOOPY_TEMPLATE[i].x[nptmp]   = Trest ;
	SNOOPY_TEMPLATE[i].y[nptmp]   = dm15 ;
	SNOOPY_TEMPLATE[i].z[nptmp]   = flux ;

	// check for 1/0 
	SNOOPY_TEMPLATE[i].wz[nptmp]  = emag > 0 ? flux : 0.0; 

	/* xxxxxxxx Aug 29, 2010 - mark for deletion ...
	x[i][nptmp] = (epoch - Tmax)/(1.0+zed);
	y[i][nptmp] = dm15;
	z[i][nptmp] = pow(10, -0.4*(mag - Mmax));
	wz[i][nptmp] = emag > 0 ? 1.0857/(z[i][nptmp]*emag) : 0.0; 
	xxxxxxxxxx */

	nchars = snprintf( tmpname[i][nptmp], 80, "%s", SN);
	SNOOPY_TEMPLATE[i].NEPOCH++ ;
      }
      fclose(fp2);
   }
   fclose(fp1);

   // print summary of templates (RK)
   for  ( i=0 ; i < MXFILT_SNOOPY; i++ ) {
     NTMP = SNOOPY_TEMPLATE[i].NSN ;
     if ( NTMP > 0 ) 
       printf("\t Found %2d templates for filter='%s' (%4d total epochs) \n"
	      ,SNOOPY_TEMPLATE[i].NSN
	      ,SNOOPY_TEMPLATE[i].FILTERNAME
	      ,SNOOPY_TEMPLATE[i].NEPOCH
	      );

   }

   return(0);
}


struct my_params { double dm15; int filter; } ;


double snoopy_template_fun(double t, void *p) {
   struct my_params *params = (struct my_params * ) p;
   double sigx, sigy;
   double dm15 = (params->dm15);
   int filter = (params->filter);
   int np, res;
   double mag, emag;
   sigx = sigmax(t, dm15);
   sigy = sigmay(t, dm15);
   np = SNOOPY_TEMPLATE[filter].NEPOCH ;

   res = gloes2D_n_sigma(
			 SNOOPY_TEMPLATE[filter].x
			 ,SNOOPY_TEMPLATE[filter].y
			 ,SNOOPY_TEMPLATE[filter].z
			 ,SNOOPY_TEMPLATE[filter].wz
			 ,np, &t, &dm15, 1, &sigx, &sigy, &mag, &emag);

   return (-1.0*mag);
}


int dm15temp(int filter, double dm15, double *t, int size, 
	     double *mag, double *emag, int normalize, 
	     double arg1, double arg2, double arg3, double arg4) {

   int i;
   int np, res;
   double *sigx, *sigy;
   struct my_params p;
   double *dm15s;
   const gsl_min_fminimizer_type *T;
   gsl_min_fminimizer *s;
   int iter=0, maxiter=100;
   int status;
   double m, a, b, max;
   gsl_function func;

   sigx0   = arg1;
   sigy0   = arg4;
   xscale  = arg2;
   sigxmax = arg3;

   sigx = (double *) malloc(size*sizeof(double));
   sigy = (double *) malloc(size*sizeof(double));
   dm15s = (double *) malloc(size*sizeof(double));
   for (i = 0 ; i < size ; ++i ) {
      sigx[i] = sigmax(t[i], dm15);
      sigy[i] = sigmay(t[i], dm15);
      dm15s[i] = dm15;

      mag[i] = MAG_UNDEFINED ; emag[i] = MAGERR_UNDEFINED ; // RK
   }

   np  = SNOOPY_TEMPLATE[filter].NEPOCH;
 
   
   res = gloes2D_n_sigma(
			 SNOOPY_TEMPLATE[filter].x
			 ,SNOOPY_TEMPLATE[filter].y
			 ,SNOOPY_TEMPLATE[filter].z
			 ,SNOOPY_TEMPLATE[filter].wz
			 ,np, t, dm15s, size, sigx, sigy, mag, emag);



   if (normalize) {
      gsl_set_error_handler_off();
      /* Normalize the lightcurves ...  find the maximum */
      func.function = &snoopy_template_fun;
      func.params = &p;
      T = gsl_min_fminimizer_brent;
      s = gsl_min_fminimizer_alloc(T);
      p.dm15 = dm15;
      p.filter = filter;
      status = gsl_min_fminimizer_set(s, &func, 0.0, -10.0, 20.0);
      if (status == GSL_EINVAL) {
         /* fprintf(stderr, "Normalization Warning:  gsl_min_fminimizer failed.\n"); */
         gsl_min_fminimizer_free(s);
         free(sigx);
         free(sigy);
         free(dm15s);
         return(5);
      }

      do {
         iter++;
         status = gsl_min_fminimizer_iterate(s);
         m = gsl_min_fminimizer_x_minimum(s);
         max = -1.0*gsl_min_fminimizer_f_minimum(s);
         a = gsl_min_fminimizer_x_lower(s);
         b = gsl_min_fminimizer_x_upper(s);
         status = gsl_min_test_interval(a, b, 0.0001, 0.0);
      } while (status == GSL_CONTINUE && iter < maxiter);
      for (i = 0 ; i < size ; ++i ) {
         mag[i] = mag[i]/max;
         emag[i] = emag[i]/max;
      }
      gsl_min_fminimizer_free(s);
   }
   free(sigx);
   free(sigy);
   free(dm15s);
   return(res);
}

// ************************************************
void parse_snoopy_modelinfo(char *modelPath ) {

  char infoFile[256];
  char c_get[80], cfilt[2]; 
  FILE *fp;
  int ifilt, N ;
  double mag, slope;

  // ------- BEGIN ----------

  sprintf(infoFile,"%s/%s.INFO", modelPath, SNOOPYNAME );

   if (! (fp = fopen(infoFile, "rt"))) { 
     sprintf(c1err,"Could not find info file:" ); 
     sprintf(c2err," '%s'", infoFile ) ;
     MADABORT(c1err,c2err);
   }
   
   SNOOPY_MODELINFO.TREST_RANGE[0] = -999. ;
   SNOOPY_MODELINFO.TREST_RANGE[1] = -999. ;

   SNOOPY_MODELINFO.DM15_RANGE[0] = -999. ;
   SNOOPY_MODELINFO.DM15_RANGE[1] = -999. ;

   SNOOPY_MODELINFO.NRANGE_EXCLUDE = 0 ;

   for ( ifilt=0; ifilt < MXFILT_SNOOPY; ifilt++ ) 
     {  SNOOPY_MODELINFO.PEAKMAG_REST[ifilt] = -999.0 ; }

   sprintf(SNOOPY_MODELINFO.ERRORFILE,"%s", "BLANK");
   sprintf(SNOOPY_MODELINFO.GRIDFILE,"%s", "BLANK");
   SNOOPY_MODELINFO.OPT_GRIDFILE = 0;

   while( fscanf(fp, "%s", c_get) != EOF) {

     if ( strcmp(c_get,"GRIDFILE:") == 0 ) {
       SNOOPY_MODELINFO.OPT_GRIDFILE = 1 ;
       fscanf(fp, "%s", SNOOPY_MODELINFO.GRIDFILE );
     }


     if ( strcmp(c_get,"ERRORFILE:") == 0 ) {
       fscanf(fp, "%s", SNOOPY_MODELINFO.ERRORFILE );
     }

     if ( strcmp(c_get,"TREST_RANGE:") == 0 ) {
       fscanf(fp, "%le", &SNOOPY_MODELINFO.TREST_RANGE[0] );
       fscanf(fp, "%le", &SNOOPY_MODELINFO.TREST_RANGE[1] );
       SNOOPY_MODELINFO.TREST_RANGE[0] -= 1.0E-6 ;
       SNOOPY_MODELINFO.TREST_RANGE[1] += 1.0E-6 ;
     }

     if ( strcmp(c_get,"DM15_RANGE:") == 0 ) {
       fscanf(fp, "%le", &SNOOPY_MODELINFO.DM15_RANGE[0] );
       fscanf(fp, "%le", &SNOOPY_MODELINFO.DM15_RANGE[1] );
       SNOOPY_MODELINFO.DM15_RANGE[0] -= 1.0E-6 ;
       SNOOPY_MODELINFO.DM15_RANGE[1] += 1.0E-6 ;
     }

     if ( strcmp(c_get,"SNMAG:") == 0 ) {
       fscanf(fp, "%s",  cfilt );
       fscanf(fp, "%le", &mag );
       fscanf(fp, "%le", &slope );

       ifilt = FILTINDX_SNOOPY(cfilt);
       if ( ifilt < 0 ) {
	 sprintf(c1err,"Undefined filter '%s' after SNMAG key.",  cfilt); 
	 sprintf(c2err,"Check '%s'", infoFile ) ;
	 MADABORT(c1err,c2err);
       }
       SNOOPY_MODELINFO.PEAKMAG_REST[ifilt] = mag;
       SNOOPY_MODELINFO.DM15SLOPE[ifilt]    = slope;
     }


     if ( strcmp(c_get,"EXCLUDE_RANGE:") == 0 ) {

       SNOOPY_MODELINFO.NRANGE_EXCLUDE++ ;
       N = SNOOPY_MODELINFO.NRANGE_EXCLUDE ;

       if ( N >= MXEXCLUDE_SNOOPY ) {
	 sprintf(c1err,"Number of EXCLUDE_RANGE keys exceeds bound of %d",
		 MXEXCLUDE_SNOOPY ) ;
	 sprintf(c2err,"Check '%s'", infoFile ) ;
	 MADABORT(c1err,c2err);
       }

       fscanf(fp, "%s", cfilt );       
       ifilt = FILTINDX_SNOOPY(cfilt);
       if ( ifilt < 0 ) {
	 sprintf(c1err,"Undefined filter '%s' after EXCLUDE_RANGE:", cfilt); 
	 sprintf(c2err,"Check '%s'", infoFile ) ;
	 MADABORT(c1err,c2err);
       }

       SNOOPY_MODELINFO.EXCLUDE_IFILT[N] = ifilt; 
       fscanf(fp, "%le", &SNOOPY_MODELINFO.EXCLUDE_DM15[N][0]  );       
       fscanf(fp, "%le", &SNOOPY_MODELINFO.EXCLUDE_DM15[N][1]  ); 
       fscanf(fp, "%le", &SNOOPY_MODELINFO.EXCLUDE_TREST[N][0] );       
       fscanf(fp, "%le", &SNOOPY_MODELINFO.EXCLUDE_TREST[N][1] ); 

       /*
       printf(" xxx EXCLUDE %s  %4.2f < DM15 < %4.2f && %6.2f < T < %6.2f\n"
	      ,	  cfilt
	      ,SNOOPY_MODELINFO.EXCLUDE_DM15[N][0]
	      ,SNOOPY_MODELINFO.EXCLUDE_DM15[N][1]
	      ,SNOOPY_MODELINFO.EXCLUDE_TREST[N][0]
	      ,SNOOPY_MODELINFO.EXCLUDE_TREST[N][1] );
       */

     }

   }
   fclose(fp);

   printf("\t %s TREST_RANGE: %6.1f  to %6.1f days \n"
	  ,SNOOPYNAME
	  ,SNOOPY_MODELINFO.TREST_RANGE[0] 
	  ,SNOOPY_MODELINFO.TREST_RANGE[1] );

   printf("\t %s DM15_RANGE:  %6.1f  to %6.1f \n"
	  ,SNOOPYNAME
	  ,SNOOPY_MODELINFO.DM15_RANGE[0] 
	  ,SNOOPY_MODELINFO.DM15_RANGE[1] );

} // end of parse_snoopy_modelinfo


// ************************************************
void read_snoopy_errors(char *modelPath) {

  // Created Aug 2011 by R.Kessler
  // read snoopy errors from text file and initialize interpolation map.
  // Errors skipped if using GRIDMAP option because the GRIDMAP
  // includes the errors.

  int ifilt, ISIZE, NVAR, NFUN, IDMAP, I8, I8p, I4, gzipFlag ;
  double TREST, DM15, ERR;
  FILE *fp;
  char fullName[200], cfilt[12];
  //  char fnam[] = "read_snoopy_errors" ;
  

  // -------------- BEGIN --------------

  // don't bother if using GRIDFILE
  if ( SNOOPY_MODELINFO.OPT_GRIDFILE ) { return ; }

  fp = snana_openTextFile(0,modelPath, SNOOPY_MODELINFO.ERRORFILE, 
			  fullName, &gzipFlag );

  if ( fp == NULL ) {
    sprintf(c1err,"Could not find snoopy error file");
    sprintf(c1err,"'%s'", SNOOPY_MODELINFO.ERRORFILE );
    MADABORT(c1err,c2err);
  }

  printf("\n");
  printf(" Read SNooPy errors from \n   %s\n", fullName);

  NVAR = 2; // Number of variables to define ERRORs (TREST and DM15)
  NFUN = 1; // just one function to store (ERR)

  I8  = sizeof(double);
  I8p = sizeof(double*);
  I4  = sizeof(int);

  // allocate temp memory for each filter to store error map
  for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) {

    SNOOPY_MODELERR[ifilt].temp_MAP    = 
      (double**)malloc(I8p*NVAR);

    SNOOPY_MODELERR[ifilt].temp_MAP[0] = 
      (double*)malloc(I8*MXERR_SNOOPY); // TREST

    SNOOPY_MODELERR[ifilt].temp_MAP[1] = 
      (double*)malloc(I8*MXERR_SNOOPY); // DM15

    SNOOPY_MODELERR[ifilt].temp_ERR    = 
      (double*)malloc(I8*MXERR_SNOOPY);

    SNOOPY_MODELERR[ifilt].SIZE = 0;
  }

  while( fscanf(fp, "%s", cfilt) != EOF) {

    ifilt = FILTINDX_SNOOPY(cfilt);
    readdouble(fp, 1, &TREST);
    readdouble(fp, 1, &DM15);
    readdouble(fp, 1, &ERR);

    if ( ifilt < 0 ) {
      sprintf(c1err,"Invalid filter '%s' in error map.", cfilt);
      sprintf(c2err,"See TREST=%f   DM15=%f  and ERR=%f", TREST, DM15, ERR);
      MADABORT(c1err,c2err);
    }

    SNOOPY_MODELERR[ifilt].SIZE++ ;
    ISIZE = SNOOPY_MODELERR[ifilt].SIZE ;
    SNOOPY_MODELERR[ifilt].temp_MAP[0][ISIZE-1] = TREST ;
    SNOOPY_MODELERR[ifilt].temp_MAP[1][ISIZE-1] = DM15  ;    
    SNOOPY_MODELERR[ifilt].temp_ERR[ISIZE-1]    = ERR   ; 

    if ( ISIZE >= MXERR_SNOOPY ) {
      sprintf(c1err,"MAP-SIZE(%s) exceeds bound of MXERR_SNOOPY=%d",
	      cfilt, MXERR_SNOOPY);
      sprintf(c2err,"Either increase MXERR_SNOOPY or shrink error map.");
      MADABORT(c1err,c2err);
    }
  }  // fscanf

  fclose(fp);
  // print map size for each filter
  printf(" Error Map-size(%s): \n   ", FILTLIST_SNOOPY );
  for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) 
    { printf("%d ", SNOOPY_MODELERR[ifilt].SIZE );  }
  printf("\n" );

  
  // init interp 
  // Start by allocating GRIDMAP memory
  for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) {
    ISIZE = SNOOPY_MODELERR[ifilt].SIZE ;
    SNOOPY_MODELERR_GRIDMAP[ifilt].VALMIN = (double*)malloc(I8*NVAR+I8);
    SNOOPY_MODELERR_GRIDMAP[ifilt].VALMAX = (double*)malloc(I8*NVAR+I8);
    SNOOPY_MODELERR_GRIDMAP[ifilt].VALBIN = (double*)malloc(I8*NVAR+I8);
    SNOOPY_MODELERR_GRIDMAP[ifilt].RANGE  = (double*)malloc(I8*NVAR+I8);
    SNOOPY_MODELERR_GRIDMAP[ifilt].NBIN   = (int   *)malloc(4*NVAR+4);
    SNOOPY_MODELERR_GRIDMAP[ifilt].INVMAP = (int   *)malloc(4*ISIZE+4);
    SNOOPY_MODELERR_GRIDMAP[ifilt].FUNVAL    = (double**)malloc(I8p);
    SNOOPY_MODELERR_GRIDMAP[ifilt].FUNVAL[0] = (double *)malloc(I8*ISIZE);

    IDMAP = 30 + ifilt;
    init_interp_GRIDMAP(IDMAP, "SNOOPY", ISIZE, NVAR, NFUN,0,
                        SNOOPY_MODELERR[ifilt].temp_MAP, 
                        &SNOOPY_MODELERR[ifilt].temp_ERR, 
                        &SNOOPY_MODELERR_GRIDMAP[ifilt] ); // <== returns this
  }


  // free temp memory after interp_GRIDMAP is  initialized.
  for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) {
    free(SNOOPY_MODELERR[ifilt].temp_MAP[0]);
    free(SNOOPY_MODELERR[ifilt].temp_MAP[1]);
    free(SNOOPY_MODELERR[ifilt].temp_MAP) ;
    free(SNOOPY_MODELERR[ifilt].temp_ERR) ;
  }


} // end of read_snoopy_errors



// ************************************************
void dm15temp_errors(int ifilt, double dm15, double *rest_dates, 
		       int nobs, double *relFlux_err) {


  int i, istat;
  double xval[2], ERR ;
  //  char fnam[] =  "dm15temp_errors" ;

  // ---------- BEGIN ------------

  for ( i=0; i < nobs; i++ ) {
    xval[0] = rest_dates[i];   
    xval[1] = dm15;

    if ( FIX_MODELERR > 0.0 ) 
      { ERR = FIX_MODELERR ; }
    else
      { istat = interp_GRIDMAP(&SNOOPY_MODELERR_GRIDMAP[ifilt], xval, &ERR );}

    relFlux_err[i] = ERR ;
  }

} // end of dm15temp_errors

// ************************************************
int FILTINDX_SNOOPY(char *cfilt) {

  int ifilt;

  // ------------ BEGIN -----------

  for ( ifilt=IFILT_SNOOPY_MIN; ifilt <= IFILT_SNOOPY_MAX; ifilt++ ) {
    if ( strcmp(cfilt,SNOOPY_TEMPLATE[ifilt].FILTERNAME) == 0 ) {
      return ifilt;
    }
  }

  return -9;

} // end of FILTINDX_SNOOPY


void get_LAMRANGE_snoopy(double *lammin, double *lammax) {
  *lammin = MINLAM_SNOOPY ;
  *lammax = MAXLAM_SNOOPY ;
} // end of get_LAMRANGE_snoopy

#ifdef SNGRIDREAD
// *****************************************************************
void gridinterp_snoopy(int ifilt, double dm15, 
		       int nobs, double *rest_dates, 
		       double *relFlux, double *relFlux_err ) {


  // Nov 2010 R.Kessler
  // Interpolate relative Flux on snoopy grid 
  // (made by snana using snoopy model)
  // and return relFlux[i] for each epoch 'rest_date[i]'
  // Note that fits_read_SNGRID() must be called before
  // calling this function.

  int 
    i, ioff, im
    ,index_dm15, index_Trest
    ,ILC, IPTRLC_OFF[2], IPTRLC
    ,NBIN_TREST, NBIN_DM15
    ;

  double  
    Trest
    ,ratio_Trest
    ,ratio_dm15
    ,gridFlux2[2][2]     // [Trest][dm15]
    ,gridFlux1[2]
    ,gridErr2[2][2]     // [Trest][dm15]
    ,gridErr1[2]
    ,dif
    ,TREST_MIN,  DM15_MIN
    ,TREST_MAX,  DM15_MAX
    ,TREST_BIN,  DM15_BIN
    ;

  short   I2TMP, I2ERR ;


  // ----------- BEGIN -------------

  NBIN_TREST = SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_TREST] ;
  NBIN_DM15  = SNGRID_SNOOPY.NBIN[IPAR_GRIDGEN_SHAPEPAR] ;

  TREST_BIN  = SNGRID_SNOOPY.BINSIZE[IPAR_GRIDGEN_TREST] ;
  DM15_BIN   = SNGRID_SNOOPY.BINSIZE[IPAR_GRIDGEN_SHAPEPAR] ;

  // get dm15 grid-index

  index_dm15 = INDEX_GRIDGEN(IPAR_GRIDGEN_SHAPEPAR, dm15, &SNGRID_SNOOPY );
  if ( index_dm15 >= NBIN_DM15 ) 
    { index_dm15-- ; }

  ILC           = index_dm15; 
  IPTRLC_OFF[0] = SNGRID_SNOOPY.PTR_GRIDGEN_LC[ILC+0] ; 
  IPTRLC_OFF[1] = SNGRID_SNOOPY.PTR_GRIDGEN_LC[ILC+1] ;

  DM15_MIN   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][index_dm15] ;
  DM15_MAX   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_SHAPEPAR][index_dm15+1] ;
  ratio_dm15 = (dm15 -  DM15_MIN)/DM15_BIN ;

  /*
  if ( ratio_dm15 < -1.0E-4 || ratio_dm15 > 1.000000001 ) {
    sprintf(c1err, "Invalid ratio_dm15 = %f for dm15=%6.3f (index_dm15=%d)", 
	    ratio_dm15, dm15, index_dm15 );
    sprintf(c2err,"DM15(MIN,MAX,BIN) = %6.3f,%6.3f,%6.3f", 
	    DM15_MIN, DM15_MAX, DM15_BIN );
    MADABORT(c1err,c2err);
  }
  */


  // make sure that 1st word in BEGIN-LC marker
  I2TMP = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC_OFF[0]];
  if ( I2TMP != MARK_GRIDGEN_LCBEGIN ){
    sprintf(c1err,"First I*2 word of ILC=%d is %d .", ILC, I2TMP );
    sprintf(c2err,"But expected %d", MARK_GRIDGEN_LCBEGIN );
    MADABORT(c1err,c2err);
  }

  // make sure that 2nd word is 8-bits of ILC
  I2TMP = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC_OFF[0]+1];
  if ( I2TMP != ( ILC & 127 ) ){
    sprintf(c1err,"2nd I*2 word of ILC=%d is %d .", ILC, I2TMP );
    sprintf(c2err,"But expected %d", (ILC & 127) );
    MADABORT(c1err,c2err);
  }

  ioff       = (ifilt-1)*NBIN_TREST + NPADWD_LCBEGIN-1 ;

  for ( i=0; i < nobs; i++ ) {

    Trest       = rest_dates[i];
    index_Trest = INDEX_GRIDGEN(IPAR_GRIDGEN_TREST,Trest, &SNGRID_SNOOPY );
    if ( index_Trest >= NBIN_TREST ) 
      { index_Trest-- ; }

    TREST_MIN   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][index_Trest] ;
    TREST_MAX   = SNGRID_SNOOPY.VALUE[IPAR_GRIDGEN_TREST][index_Trest+1] ;
    ratio_Trest = (Trest - TREST_MIN)/TREST_BIN ;

    /*
    if ( ratio_Trest < 0.0 || ratio_Trest > 1.000000001 ) {
      sprintf(c1err, "Bad ratio_Trest = %f at Trest(%s)=%6.2f  dm15=%6.2f", 
	      ratio_Trest, SNOOPY_TEMPLATE[ifilt].FILTERNAME, Trest, dm15 );
      sprintf(c2err,"index_Trest=%d  TREST(MIN,MAX,BIN)=%6.2f,%6.2f,%6.2f", 
	      index_Trest, TREST_MIN, TREST_MAX, TREST_BIN );
      MADABORT(c1err,c2err);
    }
    */


    // get gridFlux2 at the 4 corners bounding Trest,dm15
    for ( im=0; im <=1; im++ ) {    // dm15 edges
      IPTRLC      = IPTRLC_OFF[im] + ioff + index_Trest;
      I2TMP       = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC] ;
      I2ERR       = SNGRID_SNOOPY.I2GRIDGEN_LCERR[IPTRLC] ;
      gridFlux2[0][im]  = (double)I2TMP / GRIDGEN_I2LCPACK ;
      gridErr2[0][im]   = (double)I2ERR / GRIDGEN_I2LCPACK ;

      IPTRLC      = IPTRLC_OFF[im] + ioff + (index_Trest+1) ;
      I2TMP       = SNGRID_SNOOPY.I2GRIDGEN_LCMAG[IPTRLC] ;
      I2ERR       = SNGRID_SNOOPY.I2GRIDGEN_LCERR[IPTRLC] ;
      gridFlux2[1][im]  = (double)I2TMP / GRIDGEN_I2LCPACK ;
      gridErr2[1][im]   = (double)I2ERR / GRIDGEN_I2LCPACK ;

      dif = gridFlux2[1][im] - gridFlux2[0][im] ;
      gridFlux1[im] = gridFlux2[0][im] + (dif * ratio_Trest);

      dif = gridErr2[1][im] - gridErr2[0][im] ;
      gridErr1[im] = gridErr2[0][im] + (dif * ratio_Trest);
    }

    dif = gridFlux1[1] - gridFlux1[0] ;
    relFlux[i]  = gridFlux1[0] + (dif * ratio_dm15);

    dif = gridErr1[1] - gridErr1[0] ;
    relFlux_err[i]  = gridErr1[0] + (dif * ratio_dm15);

  } // i-loop over observations
  


} // end of gridinterp_snoopy

#endif


// ************************************************
void MADABORT(char *msg1, char *msg2) {

   char cmsg[40] = { "ABORT program on Fatal Error." };

   printf("\n");
   printf("\n");
   printf("\n   `|```````|`    ");
   printf("\n   <| o\\ /o |>    ");
   printf("\n    | ' ; ' |     ");
   printf("\n    |  ___  |     %s ", cmsg);
   printf("\n    | |' '| |     ");
   printf("\n    | `---' |     ");
   printf("\n    \\_______/    ");
   printf("\n");
   printf("\n");   
   printf(" %s \n",msg1 );   
   printf(" %s \n",msg2 );   

   exit(1);

}    //  end of "MADABORT"  
