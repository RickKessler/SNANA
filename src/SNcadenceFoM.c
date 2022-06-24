
/*************************
 Tools used to evaluate cadence Figure of Merit.
 Jan 2011 by P.Belov and S.Glazov


   History
  ~~~~~~~~~~~

 Jun 7, 2011 RK  (for snana v9_44)
  - re-organize inputs with IPAR_cadence_XXX parameters

  - define DEFAULT_cadence_XXX parameters

  - call new function arraystat(...) to comput AVG and RMS for M5SIG

  - pass new parameter MJDGAP_IGNORE to ignore large seasonal gaps;
    for now this parameter is not used, but it's there for Sasha
    and Pavel.
 
 Jun 16, 2011 P.Belov (for snana v9_45)
  - some erros are patched
  - sort procedure is added for m5sig_array
 
 Jun 22, 2011 P.Belov (for snana v9_45)
  - the subroutine SNcadenceFoM and dependencies are modified in order to
    take into account large gaps between epochs.

    RK: if just 1 epoch then return .001 instead of -9


 Jan 12, 2013  MAXNPIECES_BETWEEN_GAPS -> 20  (was 5)
               to allow for full LSST survey period.

 Mar 17 2013: add local GaussRan_FOM() that does not use a list.
              Removes chance that some random arrays are the same.

 Jun 1 2014  MXRAN_MJDARRAY -> 2000 (was 1000)
************************/

#define MAXRAN_MJDARRAY 2000   // max length of internal MJD arrays
#define MINOBS_SNCADENCFOM  2  // min Nobs to compute FoM
#define NPAR_SNCADENCEFOM   7
#define MAXNPIECES_BETWEEN_GAPS 20  // max number of pieces of epochs in array
#define MAXNRANDEPOCHS 10000 // max number of randomly generated epochs


#define MJDGAP_IGNORE_DEFAULT 100. // ignore MJD gaps longer than this

// define parameter indices if user passes parameters
#define IPAR_cadence_NRAN_MJDARRAY  0 // number of random MJD arrays
#define IPAR_cadence_M5SIG_AVG      1 // <M5SIG> (5-sigma limiting mag)
#define IPAR_cadence_M5SIG_RMS      2 // rms spread on above
#define IPAR_cadence_QUALITY        3 // 
#define IPAR_cadence_CUTOFF_DIFMJD  4 // max MJD-dif  with zero added to chi2
#define IPAR_cadence_CUTOFF_FRAC    5 // above = frac * avg MJD-dif in cadence
#define IPAR_cadence_MJDGAP_IGNORE  6 // ignore MJDGAP longer than this

#define DEFAULT_cadence_NRAN_MJDARRAY   500
#define DEFAULT_cadence_M5SIG_AVG      -9.0  // ==> internally calculated
#define DEFAULT_cadence_M5SIG_RMS      -9.0  // ==> internally calculated
#define DEFAULT_cadence_QUALITY         3.0
#define DEFAULT_cadence_CUTOFF_DIFMJD  -9.0  // ==> internally calculated
#define DEFAULT_cadence_CUTOFF_FRAC     0.5 
#define DEFAULT_cadence_MJDGAP_IGNORE   100.


//char MSGERR1[80];
//char MSGERR2[80];

double cadence_randMJD (double MJDMIN, double MJDMAX);

double cadence_eqMJD (int k, int N, double MJDMIN, double MJDMAX);

double cadence_chi2time (int N, double *MJDMIN, double MJDMAX);

double cadence_m5sigRAND (double M5AVER, double M5SIG3, double QUALITY);

double cadence_chi2total ( int FirstItem, int N, double *MJDLIST, 
			   double *M5SIGLIST, 
			   double M5AVER, double M5SIG3, double QUALITY, 
			   double CUTOFF, char *SIMLIB_IDENTIFIER);

int cadence_gen_RANDMJD (int FirstItem, int N, double *MJDRAND, double MJDMIN, double MJDMAX);

int cadence_gen_RANDM5SIG (int FirstItem, int N, double *M5SRAND, double M5AVER, 
			   double M5SIG3, double QUALITY) ;

double SNcadenceFoM (int OPTMASK, 
		     int Nobs, double *MJDLIST, double *M5SIGLIST, 
		     double *parList, char *SIMLIB_IDENTIFIER );


double GaussRan_FOM(void) ; // does NOT use RANLIST



// **************************************************
//
//                BEGIN FUNCTIONS
//
// **************************************************



// ===========================================
double GaussRan_FOM(void) {

  // Return Gaussian random number
  // Since we don't care about syncing randoms,
  // the RANLIST is NOT used.

  double R1, R2 ;
  double R,  V1, V2, FAC, G ;
  long int i8 ;

  // --------------- BEGIN ----------------
 BEGIN:

  i8 = random();
  R1 = (double)i8 / (double)RAND_MAX ;

  i8 = random();
  R2 = (double)i8 / (double)RAND_MAX ;

  V1 = 2.0 * R1 - 1.0;
  V2 = 2.0 * R2 - 1.0;
  R  = V1*V1 + V2*V2 ;
  if ( R >= 1.0 ) { goto BEGIN ; }
  FAC = sqrt(-2.*log(R)/R) ;
  G = V2 * FAC ;

  return G ;
}  // end of Gaussran

// ============================================
double cadence_randMJD (double MJDMIN, double MJDMAX) {
  // function returns random MJD: MJDMIN<=MJD<=MJDMAX
  // 1) use the predefined list of pseudorandom numbers
  //    return MJDMIN + rand * (MJDMAX-MJDMIN); 
  // 2) use the standart pseudorandom generator in order to
  //    get the pseudorandom numbers in the range [0,1] with accuracy 0.0001
  return MJDMIN + ((double)(rand() % 10000) / 10000.0) * (MJDMAX-MJDMIN);
}

// ========================================================
double cadence_eqMJD (int k, int N, double MJDMIN, double MJDMAX) {
  // function returns equidistant i-th MJD i=0,...,N
  //   0             N
  //  |0 1 2 3 4 5 6 7|
  // full amount of MJDs = N+1; N in the denominator
  return MJDMIN + (double)k * (MJDMAX-MJDMIN) / (double)N;
}

// ==============================================================
double cadence_m5sigRAND (double M5AVER, double M5SIG3, double QUALITY) {
  // function generates the gauss distributed random value for the magnitude
  //  in the range (-QUALITY, QUALITY) with the mean = 0
  // Feb 21 2013 RK : abs -> fabs
  double rand_value;

  rand_value = GaussRan_FOM();
  
  while ( fabs(rand_value) > QUALITY ) {
    rand_value = GaussRan_FOM();
  }

  return  M5AVER + (rand_value * M5SIG3);
}

// =====================================================
double cadence_chi2time (int N, double *MJDLIST, double CUTOFF) {
  // returns $\chi^2$ of the given array MJDLIST[N]
  //  some contributions in the $\chi^2$ are rejected ( cut-off parameter )
  int i ;
  double temp
   ,MJDMIN
   ,MJDMAX;
  double tchi2 = 0.0;
  double mjd_array[MAXRAN_MJDARRAY];

  // Fill in the arrays for epochs by the values of the epochs
  for ( i=0; i < N; i++ ) {
    mjd_array[i] = *(MJDLIST+i);
  }

  // Save the first and last time knots
  MJDMIN = mjd_array[0];
  MJDMAX = mjd_array[N-1];

  // Chi2 calculation
  for ( i=0; i < N; i++ ) {
    //temp = mjd_array[i] - eqMJD(i,N-1,MJDMIN,MJDMAX);
    temp = mjd_array[i] - ( MJDMIN + ( MJDMAX - MJDMIN ) * (double) i / (double) (N-1) );
    if ( fabs(temp) > CUTOFF ) {
        tchi2 += (temp - CUTOFF) * (temp - CUTOFF) / (CUTOFF * CUTOFF);
    }
  }

  return tchi2;
} // end of function 'chi2time'


// ====================================================================
// ====================================================================
double cadence_chi2total( int FirstItem, int N, double *MJDLIST, 
			  double *M5SIGLIST, double M5AVER, 
			  double M5SIG3, double QUALITY, double CUTOFF,
			  char *SIMLIB_IDENTIFIER ) {

  // function calculates FoM according to the explanation in the paper
  // 
  // Jul 22, 2011 (RK) define deltaChi2Mag_uni and deltaChi2Mag_grid
  //                   to test these quantities separately.
  //                   Pass SIMLIB_IDENTIFIER for special dumps.
  //

  double 
    MJDMIN, MJDMAX
    ,chi2Tot  // total chi2 for given epochs (mag contribution is included)
    ,deltaChi2Mag      // magnitud contribution
    ,deltaChi2Mag_uni  // from non-uniform grid
    ,deltaChi2Mag_grid // from grid spacing
    ,mjd_array[MAXRAN_MJDARRAY]
    ,m5sig_array[MAXRAN_MJDARRAY]
    ,mjd_copy[MAXRAN_MJDARRAY]
    ,temp, temp0, temp1, temp2, arg
    ;

  int ind[MAXRAN_MJDARRAY];   // array for indexes of sorted elements
  int i,j,it;

  // Copy the values into local arrays
  for ( i=0; i < N; i++ ) {
    mjd_array[i]   = *(MJDLIST+FirstItem+i);
    m5sig_array[i] = *(M5SIGLIST+FirstItem+i);
    ind[i] = i;
  }

  // Sort the arrays and the indexes according to the sorted MJDs
  for ( i=0; i < N-1; i++ ) {
    for ( j=0; j < N-1; j++ ) {
      if ( mjd_array[j] > mjd_array[j+1] ) {
        temp = mjd_array[j+1]; 
	mjd_array[j+1] = mjd_array[j]; 
	mjd_array[j]   = temp;
        it        = ind[j+1]; 
	ind[j+1] = ind[j]; 
	ind[j]   = it;
      }
    }
  }

  // Define min and max values of MJDs
  MJDMIN = mjd_array[0];
  MJDMAX = mjd_array[N-1];

  // Calculate chi2Tot ================================================

  // chi2 of the MJD array with N epochs
  temp0 = cadence_chi2time( N, mjd_array, CUTOFF );
  chi2Tot = temp0;

  // Loop over all epochs
  for ( i=0; i < N; i++ ) {

    if ( m5sig_array[ind[i]] < M5AVER ) {
    // add contribution from one epoch removed from the cadence.
    //  otherwise: no extra penalty for total chi2.

      // Create the array without one epoch
      for ( j=0;   j < i; j++ ) {
        mjd_copy[j] = mjd_array[j];
      }
      for ( j=i+1; j < N; j++ ) {
        mjd_copy[j-1] = mjd_array[j];
      }

      // chi2 of the MJD array with N-1 epochs
      temp1 = cadence_chi2time( N-1, mjd_copy, CUTOFF );

      // calulate $\Delta \chi^2_{m^i}$ from the note:
      deltaChi2Mag      = 0.0;
      deltaChi2Mag_uni  = 0.0;
      deltaChi2Mag_grid = 0.0;

      if ( temp1 > temp0 ) {
      // Contribution due to worse uniformity:
        deltaChi2Mag_uni += ( temp1 - temp0 );
	deltaChi2Mag     += deltaChi2Mag_uni ;
      }

      // Now get contribution due to increased optimal grid spacing:
      temp2 = ( MJDMAX - MJDMIN ) / (double) ( N - 1 );
      if ( temp2 > CUTOFF ) {
	arg = ( temp2 / CUTOFF - 1.0 );
	//        deltaChi2Mag_grid += pow(arg,2.0);
        deltaChi2Mag_grid += (arg*arg) ;
	deltaChi2Mag      += deltaChi2Mag_grid ;

	/*
	if ( deltaChi2Mag_grid > 20.0  ){
	  printf(" xxxx deltaChi2Mag_grid = %f for MJD=%9.3f (%s) \n",
		 deltaChi2Mag_grid, mjd_array[i], SIMLIB_IDENTIFIER );
	}
	*/
      }

      
      // now time for extra penalty to the total chi2:
      if ( m5sig_array[ind[i]]  <= M5AVER - QUALITY * M5SIG3) {
           chi2Tot += deltaChi2Mag;
      }
      else {
           chi2Tot += deltaChi2Mag
            *( M5AVER - m5sig_array[ind[i]] ) / ( QUALITY * M5SIG3 )
            *( M5AVER - m5sig_array[ind[i]] ) / ( QUALITY * M5SIG3 );
      }

    } // end of if condition

  } // end of loop over all epochs

  return chi2Tot;
} // end of the function 'chi2total'

// ===========================================================================
int cadence_gen_RANDMJD (int FirstItem, int N, double *MJDRAND, 
			 double MJDMIN, double MJDMAX) {
  // function generates the random MJD array

  int i,j;
  double arrayMJD[MAXRAN_MJDARRAY];
  double temp;

  // Generate MJD array
  for ( i=0; i < N; i++ ) {
    arrayMJD[i] = cadence_randMJD(MJDMIN,MJDMAX);
  }
  
  // Sort
  for ( i=0; i < N-1; i++ ) {
    for ( j=0; j < N-1; j++ ) {
      if ( arrayMJD[j] > arrayMJD[j+1] ) {
        temp = arrayMJD[j+1]; arrayMJD[j+1] = arrayMJD[j]; arrayMJD[j] = temp;
      }
    }
  }

  // Write new array over previous
  for ( i=0; i < N; i++ ) {
    *(MJDRAND+FirstItem+i) = arrayMJD[i];
  }

  return 1;
}

// =======================================================================
int cadence_gen_RANDM5SIG (int FirstItem, int N, double *M5SRAND, double M5AVER, 
			     double M5SIG3, double QUALITY) {
  // function generates the random magnitude array according to
  //  the gauss distribution with mean = M5AVER and $\sigma$ = M5SIG3
  //  QUALITY defines the range (QUALITY=3 means $3\sigma$)

  int i;
  double arrayM5SIG[MAXRAN_MJDARRAY];

  // Generate MJD array
  for ( i=0; i < N; i++ ) {
    arrayM5SIG[i] = cadence_m5sigRAND(M5AVER, M5SIG3, QUALITY);
  }

  // Write new array over previous
  for ( i=0; i < N; i++ ) {
    *(M5SRAND+FirstItem+i) = arrayM5SIG[i];
  }

  return 1;
}


// ====================================================================
// ====================================================================
double SNcadenceFoM( int OPTMASK, int Nobs, double *MJDLIST, double *M5SIGLIST,
		     double *parList, char *SIMLIB_IDENTIFIER ) {

  // Template function Created Sep 27, 2010 by R.Kessler
  // Return figure of merit for observered cadence specified by
  // list of 'Nobs' MJDs and 5-sigma limiting magnitudes.
  // Assumes that all observations are in the same filter,
  // so this function should be called separately for each filter.

/***
 *
 * Figure of Merit Evaluator
 * A.Glazov & P.Belov (DESY), December 2010
 *
 * Description:
 *  Function returns the Figure of Merit for each couple of arrays:
 *   array of epochs (MJDLIST) and array of M5SIG luminosities (M5SIGLIST).
 *   The MJDLIST array is obtained from SIMLIB file for each filter and
 *   the M5SIGLIST array is calculated from the data in SIMLIB file.
 *
 * A brief strategy is as follows. A number of gaps is found and 
 * this gaps have impact on the parameters, for example CUTOFF. 
 * The main loop over all epochs in an array selects pieces between 
 * gaps and calculate FoM for each piece of epochs. Randoms arrays 
 * are generated for this piece and FoMs for this arrays are also 
 * calculated. The weighted average of all pieces are evaluated for 
 * real and randomly generated FoMs.
 *
 * Input Function args:
 *  
 * OPTMASK buts (LSB=1):
 *   Bit 1 => use *parList instead of default values. These params are 
 *      passed as follows:
 *        snlc_sim.exe myinput.txt  SIMLIB_CADENCEFOM_PARLIST <par1> ... <parN>
 *      where the par-definitions are given below where they are stripped.
 *      Can also put "SIMLIB_CADENCEFOM_PARLIST: <par1> ... <parN>"
 *      inside your sim-input file.
 *
 *   Bit 2 => Dump info to screen.
 *
 *  Nobs      : Number of observations passed
 * *MJDLIST   : list of MJD dates for each obs
 * *M5SIGLIST : list of M5sigma limiting mags for each obs
 * *parList   : see IPAR_cadence_XXX parameters at top of file
 *
 * SIMLIB_IDENTIFIER: text string to identify SIMLIB entry for
 *                    dump or error message.
 *
 *****************
 * Modifications:
 *****************
 * P. Belov
 * Jun 16, 2011. Small errors were patched. Thanks to Rick!
 *               The sort procedure was added for m5sig_array.
 *
 * P. Belov
 * Jun 22, 2011. The subroutine was modified in order to take into account
 *               large gaps between epochs.
 *
 *
 ***/


  int 
    OPT_PARAM, OPT_DUMP 
    ,i,j,k,it
    ,belowBins    // Number of FoMs of random arrays below the given FoM
    ,ind[MAXRAN_MJDARRAY]    // array for indexes of sorted elements
    ,NRAN_t      // temp for NRAN_MJDARRAY
    ,LDMP_LOCAL
    ,NMJDGAP_IGNORE = 0 // P.Belov: Number of MJD gaps which we will ignore
    ,MJDMIN_PIECE // a number of the first observation in this piece of MJDs
    ,MJDMAX_PIECE // a number of the last observation in this piece of MJDs
    ,MJD_PIECE // a number of epochs in this piece of MJDs
    ;

  double 
    MJDMIN, MJDMAX , MJDRANGE, DIFMJD_AVG
    ,CUTOFF_DIFMJD, CUTOFF_FRAC
    ,FoM[MAXNPIECES_BETWEEN_GAPS][MAXNRANDEPOCHS]   // FoM of random generated array of epochs
    ,FoM_0[MAXNPIECES_BETWEEN_GAPS]    // FoM of the given array of epochs
    ,ResultFoM
    ,ResultFoM_0
    ,RandFoM[MAXNRANDEPOCHS]
    ,temp, XNobs
    ,CUTOFF_t        // Temp variable for CUTOFF
    ,M5SIG_AVG_t     // Temp variable for M5SIG_AVG
    ,M5SIG_RMS_t     // Temp variable for M5SIG_RMS
    ,QUALITY_t       // Temp variable for QUALITY
    ,mjd_array[MAXRAN_MJDARRAY]
    ,mjd_rand_array[MAXRAN_MJDARRAY]
    ,m5sig_array[MAXRAN_MJDARRAY]
    ,m5sig_rand_array[MAXRAN_MJDARRAY]
    ,MEDIAN
    ;

  double result;    // save fraction of FoM bins above the current FoM

  char BANNER[100];
  char fnam[] = "SNcadenceFoM" ;

  // -------------- BEGIN --------------

  // check OPTMASK bits
  OPT_PARAM = (OPTMASK & 1) ;
  OPT_DUMP  = (OPTMASK & 2) ;


  // define a few quick things from the inputs
  XNobs  = (double)Nobs ;        // convert to double for later use

  // bail if there are not enough observations
  if ( Nobs <= MINOBS_SNCADENCFOM ) {
    result = 0.0;
    return result;
  }


  // Fill in the arrays for epochs by the values of the epochs
  for ( i=0; i < Nobs; i++ ) {
    mjd_array[i]   = MJDLIST[i] ;
    m5sig_array[i] = M5SIGLIST[i] ;
    ind[i] = i;  // the indexer
  }

  // Sort
  for ( i=0; i < Nobs-1; i++ ) {
    for ( j=0; j < Nobs-1; j++ ) {
      if ( mjd_array[j] > mjd_array[j+1] ) {
        temp           = mjd_array[j+1]; 
	mjd_array[j+1] = mjd_array[j]; 
	mjd_array[j]   = temp;
	temp = m5sig_array[j+1];
	m5sig_array[j+1] = m5sig_array[j];
	m5sig_array[j] = temp;
        it       = ind[j+1]; 
	ind[j+1] = ind[j]; 
	ind[j]   = it;
      }
    }
  }

  MJDMIN     = mjd_array[0] ;
  MJDMAX     = mjd_array[Nobs-1] ;
  
  MJDRANGE   = MJDMAX - MJDMIN ;


  // get the number of gaps
  for ( i=1; i < Nobs; i++ ) {
    temp = mjd_array[i] - mjd_array[i-1] ;
    if ( temp > MJDGAP_IGNORE_DEFAULT ) {
        NMJDGAP_IGNORE += 1 ;
        MJDRANGE -= temp ; // subtract a magnitude of a gap
    }
  }
  
  DIFMJD_AVG = MJDRANGE /(double)(Nobs-1-NMJDGAP_IGNORE) ;

  if ( OPT_PARAM  ) {  // user-passed values

    NRAN_t          = (int)*(parList+IPAR_cadence_NRAN_MJDARRAY) ;
    M5SIG_AVG_t     = *(parList+IPAR_cadence_M5SIG_AVG ) ;
    M5SIG_RMS_t     = *(parList+IPAR_cadence_M5SIG_RMS ) ;
    QUALITY_t       = *(parList+IPAR_cadence_QUALITY   ) ;

    CUTOFF_DIFMJD = *(parList+IPAR_cadence_CUTOFF_DIFMJD ) ; 
    CUTOFF_FRAC   = *(parList+IPAR_cadence_CUTOFF_FRAC ) ;
    if ( CUTOFF_DIFMJD > 0.0 ) 
      {  CUTOFF_t = CUTOFF_DIFMJD ; }    
    else 
      { CUTOFF_t = CUTOFF_FRAC * DIFMJD_AVG ;  }
  }

  else {    // default parameters

    NRAN_t          = DEFAULT_cadence_NRAN_MJDARRAY ;
    QUALITY_t       = DEFAULT_cadence_QUALITY ;
    CUTOFF_FRAC     = DEFAULT_cadence_CUTOFF_FRAC ;
    CUTOFF_t        = CUTOFF_FRAC * DIFMJD_AVG ;
    CUTOFF_DIFMJD   = DEFAULT_cadence_CUTOFF_DIFMJD ; // not used
  }

  arrayStat(Nobs, m5sig_array,               // <== input
	    &M5SIG_AVG_t, &M5SIG_RMS_t, &MEDIAN );  // <== returned

  LDMP_LOCAL = 0 ; // set to whatever needed for debugging

  // dump here after all values are set (RK Jan 26, 2011)
  if  ( OPT_DUMP || LDMP_LOCAL ) {
    sprintf(BANNER,"SNcadenceFoM Parameters for %s [OPTMASK=%d]", 
	    SIMLIB_IDENTIFIER, OPTMASK);

    printf("\n ************************************************** \n");
    printf("   %s  \n" , BANNER );
    fflush(stdout);

    printf("\t Number of SIMLIB observations to examine: %d \n", Nobs );
    printf("\t Number of random MJD arrays to generate: %d \n", NRAN_t);
    printf("\t <M5SIG> = %6.3f (rms=%5.3f) \n", M5SIG_AVG_t, M5SIG_RMS_t );
    printf("\t QUALITY = %6.3f  \n", QUALITY_t );
    printf("\t CUTOFF=%6.3f days (CUTOFF_DIFMJD=%6.3f, CUTOFF_FRAC=%5.3f) \n", 
	   CUTOFF_t, CUTOFF_DIFMJD, CUTOFF_FRAC );
    //    debugexit("bla");
    fflush(stdout);
  }


  // error check arguments after dump-option above
  if ( CUTOFF_DIFMJD > 0.0 && CUTOFF_FRAC > 0.0 ) {
    sprintf(MSGERR1,"CUTOFF_DIFMJD=%6.2f and CUTOFF_FRAC=%6.3f",
	    CUTOFF_DIFMJD, CUTOFF_FRAC );
    sprintf(MSGERR2,"but both cannot be > 0. Error at %s", SIMLIB_IDENTIFIER );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
  }

  if ( NRAN_t >= MAXRAN_MJDARRAY ) {
    sprintf(MSGERR1,"NRAN_t=%d exceeds bound of MAXRAN_MJDARRAY=%d", 
	    NRAN_t, MAXRAN_MJDARRAY);
    sprintf(MSGERR2, "Error at %s", SIMLIB_IDENTIFIER );
    errmsg(SEV_FATAL, 0, fnam, MSGERR1, MSGERR2 );
  }
  // !!! error message; SHOULD BE FATAL !!!
  if ( (NMJDGAP_IGNORE+1) >= MAXNPIECES_BETWEEN_GAPS ) { 
    // number of MJD pieces between gaps bigger than defined constant in header
    printf(" Error: (NMJDGAP_IGNORE+1) >= MAXNPIECES_BETWEEN_GAPS  !!! \n");
    printf("         NMJDGAP_IGNORE = %d\n", NMJDGAP_IGNORE );
    printf("         MAXNPIECES_BETWEEN_GAPS = %d\n", MAXNPIECES_BETWEEN_GAPS);
    return 0.0;
  }

  // define initial values for variables
  MJDMIN_PIECE = 0 ; // denote the first epoch in a piece
  MJDMAX_PIECE = 0 ; // denote the last epoch in a piece
  k = 0; // counter for number of pieces between gaps

  // main loop over epochs; the gaps are rejected
  for ( i=1; i < Nobs; i++ ) {

    // reject a gap
    if ( ((mjd_array[i]-mjd_array[i-1]) > MJDGAP_IGNORE_DEFAULT)
     || (i==(Nobs-1)) ) {

      // !!!  error message; should be FATAL !!!
      if ( (mjd_array[i]-mjd_array[i-1]) > MJDGAP_IGNORE_DEFAULT ) {
        if ( (i==(Nobs-1)) || (i==1) ) {
	  /*
          printf(" Error: Only 1 epoch in the piece !!! (%s)\n", 
	  SIMLIB_IDENTIFIER ); */
          return 0.001;
        }
        else {
          if ( (mjd_array[i-1]-mjd_array[i-2]) > MJDGAP_IGNORE_DEFAULT ) {
	    /*
            printf(" Error: Only 1 epoch in the piece !!! (%s) \n",
	    SIMLIB_IDENTIFIER ); */
            return 0.001;
          }
        }
      }

      // define a number of last epoch in a piece before the gap
      if (i==(Nobs-1)) {
        MJDMAX_PIECE = i ;  // in case of no gaps in array
      } else {
        MJDMAX_PIECE = i-1 ;
      }
      
      // define a length of the piece
      MJD_PIECE = MJDMAX_PIECE - MJDMIN_PIECE + 1 ;

        // debug --->
        //if ( (MJD_PIECE < 3) && !(i==(Nobs-1)) ) { printf(" Warning: Less than 3 epochs in the piece %d \n", MJD_PIECE) ; }
        // <--- debug
      
      // calculate fom for the piece

      FoM_0[k] = ((double) MJD_PIECE) * 
	cadence_chi2total ( MJDMIN_PIECE,
			    MJD_PIECE,
			    mjd_array,
			    m5sig_array,
			    M5SIG_AVG_t,
			    M5SIG_RMS_t,
			    QUALITY_t,
			    CUTOFF_t,
			    SIMLIB_IDENTIFIER  // RK Jul 22, 2011
			    ) ;

      if  ( OPT_DUMP || LDMP_LOCAL ) {  
	printf("\t\t k=%d MJD_PIECE=%d  <M5SIG>=%.2f +- %.2f  FoM_0=%f \n",
	       k, MJD_PIECE, M5SIG_AVG_t, M5SIG_RMS_t, FoM_0[k] );
      }

      // calculate random FoMs
      for ( j=0; j < NRAN_t; j++ ) {

          it = cadence_gen_RANDMJD ( MJDMIN_PIECE,
          			     MJD_PIECE,
			             mjd_rand_array,
			             mjd_array[MJDMIN_PIECE],
			             mjd_array[MJDMIN_PIECE+MJD_PIECE-1]
			             ) ;

          it = cadence_gen_RANDM5SIG ( MJDMIN_PIECE,
          			       MJD_PIECE,
				       m5sig_rand_array,
				       M5SIG_AVG_t,
				       M5SIG_RMS_t,
				       QUALITY_t
				       ) ;

          FoM[k][j] = ((double) MJD_PIECE) * 
	    cadence_chi2total ( MJDMIN_PIECE,
				MJD_PIECE,
				mjd_rand_array,
				m5sig_rand_array,
				M5SIG_AVG_t,
				M5SIG_RMS_t,
				QUALITY_t,
				CUTOFF_t,
				SIMLIB_IDENTIFIER // RK Jul 22 2011
				) ;
      }
      
      // move the start indexer to the first epoch after gap
      MJDMIN_PIECE = i ;
      
      // increase counter of the pieces of epochs between gaps
      k += 1;
      
    } // if MJDGAP_IGNORE_DEFAULT
  } // end of the for loop
  
  //printf (" k = %d\n",k);

  // calculate FoMs from real surveys
  ResultFoM_0 = 0.0;
  for ( i=0; i < k; i++ ) {  // here it is possible to use NMJDGAP_IGNORE+1
    ResultFoM_0 += FoM_0[i];
  }
  ResultFoM_0 = ResultFoM_0 / XNobs ;

  // collect FoMs from all pieces with randomly generated epochs
  // the collecting goes over all pieces i for each random generation j
  for ( j=0; j < NRAN_t; j++ ) {
    ResultFoM = 0.0;
    for ( i=0; i < k; i++ ) {  // here it is possible to use NMJDGAP_IGNORE+1
      ResultFoM += FoM[i][j];
    }
    RandFoM[j] = ResultFoM / XNobs ;
  }

  belowBins = 0 ;
  // calculate a number of random FoMs below given FoM
  for ( i=0; i < NRAN_t; i++ ) {
    if ( RandFoM[i] < ResultFoM_0 ) { belowBins++ ; }
  }

  // get the fraction of FoM bins above the current FoM
  result = 1.0 - (double) belowBins / (double) NRAN_t ;

  if  ( OPT_DUMP || LDMP_LOCAL ) {  
    printf("\t ResultFoM_0 = %f    belowBins=%d   k=%d\n", 
	   ResultFoM_0, belowBins, k);
  }
  // Output to the standart out stream
  //printf(" === SNcadence output after processing === \n");
  //printf("     Nobs = %d \n", Nobs);
  //printf("     FoM_0= %f \n", ResultFoM_0);
  //printf(" === Bins above the FoM: %f \n", result);

  return result;

} // end of SNcadenceFoM

