/* =====================================================

  Package to read and apply non-linearity from map(s) in file.
  [initial use for WFIRST sims]

  Sim-input key  NONLINEARITY_FILE:  <inFile>
  where inFile is of the form:

  MODELNAME:  WHATEVER

  OPTMASK_NONLIN:  1  # 1=count tot, 2=count rate

  FILTERS: abcdef
  NONLIN:   4.0E0  0.982  #   Ftot(pe)  and Flux-scale
  NONLIN:   4.0E1  0.986
  NONLIN:   4.0E2  0.988 
  etc .

  Notes:
    + interpolation is done in log10(Ftot) space
    + Ftot should be the sum of Fsource + Fsky + Fgal.
    + Repeat as many maps as needed in case different filters have 
      different non-linearities.
    + for OPTMASK_NONLIN += 2 (count-rate), first NONLIN arg is
      Ftot_pe/TEXPOSE, where TEXPOSE per observation is read from
      the PHOT data stream. For sims, add TEXPOSE([FIELD]) keys in 
      global header to automatically fill TEXPOSE in data.

======================================================= */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "sntools.h"
#include "sntools_nonlinearity.h"


// =====================================
void INIT_NONLIN(char *inFile) {

  // Oct 2024:
  // Refactor to only require "FILTERS:" and "NONLIN:".
  // No longer need START_MAP and END_MAP keys.
  
  int langC = LANGFLAG_PARSE_WORDS_C;
  FILE *fp ;
  int  imap, nmap_read, MAPSIZE, NLINE, NWD ;
  double tmpD[4];
  bool  RDFLAG, END_OF_MAP, ISKEY_FILTERS, ISKEY_NONLIN ;
  char c_get[MXPATHLEN], LINE[MXPATHLEN], KEY[40], MSG[100] ;
  char fnam[] = "INIT_NONLIN" ;
  
  // ------------ BEGIN --------------

  NMAP_NONLIN = 0 ;
  MODELNAME_NONLIN[0] = 0 ;
  NONLIN_README.NLINE = NLINE = 0 ;
  OPTMASK_NONLIN = 0 ;
  
  if ( IGNOREFILE(inFile) ) { return ; }

  fp = fopen(inFile,"rt");
  if ( !fp ) {
    sprintf(c1err,"Cannot open non-linearity file");
    sprintf(c2err,"%s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  sprintf(BANNER,"%s : prepare non-linearity map(s) in \n\t %s",
	  fnam, inFile);
  print_banner(BANNER);

  // first pass read to count how many tables for malloc
  // and to determine the model
  while( (fscanf(fp, "%s", c_get)) != EOF) {
    if ( strcmp(c_get,"MODELNAME:"    ) == 0 ) { 
      readchar(fp,MODELNAME_NONLIN); 
      NONLIN_README.LINE[NLINE][0] = 0 ; NLINE++ ;
      sprintf(NONLIN_README.LINE[NLINE],
	      "  Apply NONLINEARITY MODEL '%s'", MODELNAME_NONLIN);
      NLINE++ ;
    }

    if ( strcmp(c_get,"OPTMASK_NONLIN:") == 0 ) 
      {  readint(fp, 1, &OPTMASK_NONLIN); }
    
    if ( strcmp(c_get,"FILTERS:") == 0 ) { NMAP_NONLIN++ ; }    
  }

  rewind(fp);

  if ( strlen(MODELNAME_NONLIN) == 0 ) {
    sprintf(c1err,"Must specify  MODELNAME: <model>" );
    sprintf(c2err,"in %s", inFile);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  check_OPTMASK_NONLIN();

  // - - - -- 
  NONLIN_MAP = (NONLIN_DEF*) malloc( NMAP_NONLIN * sizeof(NONLIN_DEF));
  
  // init all maps
  for(imap=0 ; imap < NMAP_NONLIN; imap++ ) {
    NONLIN_MAP[imap].FILTERS[0] = 0 ;
    NONLIN_MAP[imap].MAPSIZE  = 0 ;
  }
  DUMPFLAG_NONLIN  = (OPTMASK_NONLIN & OPTMASK_NONLIN_DUMPFLAG);
  DEBUGFLAG_NONLIN = (OPTMASK_NONLIN & OPTMASK_NONLIN_DEBUGFLAG);
	 
  // - - - - - - -
  // read again and store each map
  nmap_read = MAPSIZE = RDFLAG = END_OF_MAP = 0 ;
  
  while ( fgets(LINE, MXPATHLEN, fp) != NULL ) {

    KEY[0] = 0 ;
    NWD = store_PARSE_WORDS(MSKOPT_PARSE_WORDS_STRING, LINE, fnam);

    if ( NWD >= 2 ) { get_PARSE_WORD(langC, 0, KEY, fnam); }
	    
    ISKEY_FILTERS = ( strcmp(KEY,"FILTERS:" ) == 0 );
    ISKEY_NONLIN  = ( strcmp(KEY,"NONLIN:"  ) == 0 );

    // check for end of map
    END_OF_MAP = RDFLAG && (NWD==0 || ISKEY_FILTERS );
    if ( END_OF_MAP ) {
      if ( MAPSIZE <= 0 ) {
	sprintf(c1err,"Found map with %d entries", MAPSIZE);
	sprintf(c2err,"Check %s", inFile);
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
      }
      sprintf(MSG,"\t Read NONLIN MAP%2.2d with %2d rows for %s",
	     imap, MAPSIZE, NONLIN_MAP[imap].FILTERS); 
      printf("%s\n", MSG);      fflush(stdout);
      sprintf(NONLIN_README.LINE[NLINE],"%s", MSG);
      NLINE++ ;
      RDFLAG = 0 ;
    }

    
    if ( ISKEY_FILTERS ) {
      MAPSIZE=0 ;  nmap_read++ ;  imap=nmap_read-1;  RDFLAG=1;
      get_PARSE_WORD(langC, 1, NONLIN_MAP[imap].FILTERS, fnam );
    }

    if ( ISKEY_NONLIN ) {
      get_PARSE_WORD_DBL(langC, 1, &tmpD[0], fnam );
      get_PARSE_WORD_DBL(langC, 2, &tmpD[1], fnam );      
      NONLIN_MAP[imap].MAPVAL[0][MAPSIZE]  = log10(tmpD[0]); //map is in flux, but stored in log flux!
      NONLIN_MAP[imap].MAPVAL[1][MAPSIZE]  = tmpD[1] ;
      MAPSIZE++ ;
      NONLIN_MAP[imap].MAPSIZE = MAPSIZE ;
    }


  } // end while

  fclose(fp);

  printf("\n"); fflush(stdout);

  NONLIN_README.LINE[NLINE][0] = 0 ; NLINE++ ;
  NONLIN_README.NLINE = NLINE ;

  //  debugexit(fnam); // xxx REMOVE
  return ;

} // end INIT_NONLIN


void   init_nonlin__(char *inFile) { INIT_NONLIN(inFile); }


void check_OPTMASK_NONLIN(){
    char fnam[] = "check_OPTMASK_NONLIN" ;
    
  int NOPT_REQUIRE = 0;

  if ( (OPTMASK_NONLIN & OPTMASK_NONLIN_COUNT_TOT ) > 0 ) {
    printf("\t Compute NONLIN from COUNT-TOTAL\n");
    NOPT_REQUIRE++;
  }
  if ( (OPTMASK_NONLIN & OPTMASK_NONLIN_COUNT_RATE) > 0 ) {
    printf("\t Compute NONLIN from COUNT/Texpose\n");    
    NOPT_REQUIRE++;
  }
  
  if ( NOPT_REQUIRE != 1 ) {
    sprintf(c1err,"Invalid NOPT_REQUIRE=%d; must require EITHER",
	    NOPT_REQUIRE);
    sprintf(c2err,"OPTMASK_NONLIN+=1(COUNT-TOTAL)  or +=2(COUNT-RATE) " );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  // xxx

  NOPT_REQUIRE = 0;

  if ( (OPTMASK_NONLIN & OPTMASK_NONLIN_PER_NEA ) > 0 ) {
    printf("\t Compute NONLIN for total flux in NEA\n");
    NOPT_REQUIRE++;
  }
  if ( (OPTMASK_NONLIN & OPTMASK_NONLIN_PER_PIX) > 0 ) {
    printf("\t Convert NONLIN per pixel to NONLIN for flux in NEA\n");    
    NOPT_REQUIRE++;
  }
  
  if ( NOPT_REQUIRE != 1 ) {
    sprintf(c1err,"Invalid NOPT_REQUIRE=%d; must require EITHER",
	    NOPT_REQUIRE);
    sprintf(c2err,"OPTMASK_NONLIN+=4(PER-NEA)  or +=8(PER-PIX) " );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }


;}


// ====================================
double GET_NONLIN(char *CCID, char *cfilt, double Texpose, double NEA, double *Fpe_list, 
		  double mag ) {

  // Inputs
  //   + CCID        = char CCID; for abort and dump messages
  //   + cfilt       = 1-char filter band
  //   + Texpose     = exposure time (sec); for count-rate nonlin option
  //   + NEA         = noise equiv area, in pixels
  //                    WARNING: do not pass NEA in arcsec^2
  //   + Fpe_list[0] = Fpe(source)
  //   + Fpe_list[1] = Fpe(sky) within NEA aperture
  //   + Fpe_list[2] = Fpe(galaxy)  
  //   + mag         = true magnitude of object (not used yet)
  //
  // Returns F_meas/F_true .
  //
  // BEWARE: if nonLinearity is not defined for *cfilt, 
  //         function returns 1.0000
  //

  double Fpe_source      = Fpe_list[0];
  double Fpe_sky         = Fpe_list[1];
  double Fpe_galaxy      = Fpe_list[2];

  // if ( DEBUGFLAG_NONLIN ) { Fpe_sky = 0.0; } // test
  
  bool NONLIN_COUNT_TOT  = (OPTMASK_NONLIN & OPTMASK_NONLIN_COUNT_TOT  ) > 0 ;
  bool NONLIN_COUNT_RATE = (OPTMASK_NONLIN & OPTMASK_NONLIN_COUNT_RATE ) > 0 ;
  bool NONLIN_PER_NEA    = (OPTMASK_NONLIN & OPTMASK_NONLIN_PER_NEA  ) > 0 ;
  bool NONLIN_PER_PIX    = (OPTMASK_NONLIN & OPTMASK_NONLIN_PER_PIX ) > 0 ;
  int LDMP = DUMPFLAG_NONLIN ;
  
  int    imap ;
  double scale_nonlin,  Fpe_tot, flux_scale_count;
  char  *ptrFilters, msg[100] ;
  char   fnam[] = "GET_NONLIN" ;

  // --------------- BEGIN ----------------

  
  if ( LDMP ) {
    printf(" xxx ----------------------------------------------- \n");
    printf(" xxx %s dump for CCID=%s  BAND=%s  Texpose=%d  NEA=%.3f  mag=%6.3f \n",
	   fnam, CCID, cfilt, (int)Texpose, NEA, mag);
    printf(" xxx \t Fpe_total(src,sky,gal) = %10.3le  %10.3le  %10.3le  (e-)\n",	 
	   Fpe_source, Fpe_sky, Fpe_galaxy);
    fflush(stdout);
  }
  

  /* xxxxxxxx mark delete Oct 15 2024 xxxxxxxxxxx
  if ( DO_SPEED_TEST_PSF ) {
    double rsq, PSF_DUMMY, sigsq = 0.334 ;
    int ncalc = 0.0 ;
    for(rsq = 0.01; rsq< 1.0; rsq += 0.01 ) {
      PSF_DUMMY = exp(-rsq/sigsq);
      ncalc++ ;
    }
    printf(" xxx %s: finished %s SPEED_TEST_PSF with %d exp calcs.\n",
	   fnam, cfilt, ncalc); fflush(stdout);
  } // end SPEED_TEST
  xxxxxxxxxxx end mark xxxxxxxxx */
  
  scale_nonlin = 1.000 ; // default is no non-linearity

  // bail if there are no maps.
  if ( NMAP_NONLIN == 0    ) { return(scale_nonlin); }
  if ( Fpe_sky    <   0.0  ) { return(scale_nonlin); }


  // - - - - - - - -

  // convert Fpe-count to count-rate if count-rate option is selected  
  if ( NONLIN_COUNT_RATE ) {
    if ( Texpose < 2 ) {
      sprintf(c1err,"Invalid Texpose = %.3f for band=%s", Texpose, cfilt);
      sprintf(c2err,"Cannot compute count-rate NONLIN");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
    }
    flux_scale_count = 1.0/Texpose ; 
    Fpe_source *= flux_scale_count ;
    Fpe_sky    *= flux_scale_count ;
    Fpe_galaxy *= flux_scale_count ;

    if ( LDMP ) {
      printf(" xxx \t Fpe_rate(src,sky,gal)  = %10.3le  %10.3le  %10.3le  "
	     "(e-/sec)\n",	 
	     Fpe_source, Fpe_sky, Fpe_galaxy);
    }
  }
  
  // - - - - - - - -

  if (NONLIN_PER_NEA) { 
    Fpe_tot = (Fpe_source + Fpe_sky + Fpe_galaxy) ;
    scale_nonlin = get_flux_scale_NONLIN(cfilt, Fpe_tot);
  }
  else if ( NONLIN_PER_PIX ) {
    // convert nonlin per pix to nea
    double NEA_over_PI = NEA/3.14159 ;
    double  Fpix_sky   = (Fpe_sky/NEA) ; //sky per pixel
    int    MXpix       = (int)4.*NEA_over_PI + 2; // max Npix for NEA
    double *PSF_grid   = (double*) malloc(MXpix * sizeof(double));
    double sum_PSF     = 0.0, invsum_PSF, PSF ;
    double sigsq_PSF   = 0.25*NEA_over_PI; // sigma^2 for effective Gauss PSF
    double RSQ_BORDER  = NEA_over_PI ;   // effective RSQ at border of NEA
    int    npix_psf = 0, i ;

    // area of square containing NEA is 4*NEA/pi ... half size of
    // side is therefore 0.5*sqrt(4*NEA/PI)
    double x, x_half  = (int)(0.5*sqrt(4.*NEA_over_PI)) + 1.0 ; 
    double y, y_half  = x_half ;
    double xoff=0.0, yoff=0.0; // perhaps later, select these randomly
    double rsq, rsq_min=9999999.0, PSFmax=0.0  ;
    double Fpix_sky_nonlin;
    
    for(x=-x_half; x < x_half; x+=1.0 ) {
      for(y=-y_half; y < y_half; y+=1.0 ) {
	rsq = (x-xoff)*(x-xoff) + (y-yoff)*(y-yoff);
	if ( rsq < RSQ_BORDER ) {
	  if ( npix_psf < MXpix) {
	    PSF                = exp(-0.5*rsq/sigsq_PSF);
	    if ( DEBUGFLAG_NONLIN && rsq != 0.0 ) { PSF=0.0; } // test
	    PSF_grid[npix_psf] = PSF;
	    sum_PSF += PSF ;
	  }
	  if ( rsq < rsq_min ) { rsq_min=rsq; PSFmax = PSF; }
	  npix_psf++ ;
	}
      } // end y
    } // end x
    
    if ( npix_psf >= MXpix ) {
      sprintf(c1err,"Invalid npix_psf=%d  exceeds  MXpix=%d", npix_psf, MXpix);
      sprintf(c2err,"BAND=%s  NEA=%.3f  RSQ_BORDER=%.3f",
	      cfilt, NEA, RSQ_BORDER );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
    }

    if ( sum_PSF <= 0.0 ) {
      sprintf(c1err,"Invalid sum_PSF = %le (npix_psf=%d)", sum_PSF, npix_psf);
      sprintf(c2err,"BAND=%s  NEA=%.3f  RSQ_BORDER=%.3f",
	      cfilt, NEA, RSQ_BORDER );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err );	
    }


    Fpix_sky_nonlin = Fpix_sky * get_flux_scale_NONLIN(cfilt,Fpix_sky);
    
    if ( LDMP ) {
      printf(" xxx \t Fpix_sky = %.4f -> %.4f(nonlin)\n",
	     Fpix_sky, Fpix_sky_nonlin );
      printf(" xxx \t npix_psf=%d  MXpix=%d  sum_PSF=%.4f  PSFmax/sum_PSF = %.3f\n",
	     npix_psf, MXpix, sum_PSF, PSFmax/sum_PSF ); fflush(stdout);
    }
    
    // loop over PSF and compute weighted sums
    double Fpix_tot, Fpix_tot_nonlin ;
    double Fpix_wgtsum_src = 0.0, Fpix_wgtsum_nonlin=0.0 ;
    invsum_PSF = 1.0/sum_PSF;

    
    for(i=0; i < npix_psf; i++ ) {
      PSF_grid[i] *= invsum_PSF;
      PSF          = PSF_grid[i] ;
      Fpix_wgtsum_src += ( Fpe_source * PSF * PSF); // denominator
      
      Fpix_tot        = (Fpe_source * PSF + Fpix_sky) ;
      Fpix_tot_nonlin = Fpix_tot * get_flux_scale_NONLIN(cfilt,Fpix_tot) ;

      Fpix_wgtsum_nonlin += (Fpix_tot_nonlin - Fpix_sky_nonlin)* PSF;
    }

    
    scale_nonlin = Fpix_wgtsum_nonlin / Fpix_wgtsum_src ;
    
    free(PSF_grid);

  } // NONLIN_PER_PIX
  
  // - - - - - - 
  if ( LDMP ) { 
    printf(" xxx \t --> scale_nonlin=%6.4f \n",  scale_nonlin); 
    fflush(stdout);
  }

  return(scale_nonlin);

} // end GET_NONLIN

double get_nonlin__(char *CCID, char *cfilt, double *Texpose, double *NEA, double *Fpe_list, 
		    double *genmag) {
  
  double F_scale = GET_NONLIN(CCID, cfilt, *Texpose, *NEA, Fpe_list,  *genmag);
  return(F_scale);
}

// =============================
double get_flux_scale_NONLIN(char *cfilt, double flux) {

  // Return F(+nonlin) / F(perfect linearity)

  if ( flux <= 0.0 ) { return 1.000; }
  
  double log10_flux = log10(flux);
  double f_scale = 0.0 ;
  int    OPT_INTERP = 1;  // 1=linear interp
  int    imap ;
  char  *ptrFilters, msg[100] ;
  char fnam[] = "get_flux_scale_NONLIN";

  // ---------- BEGIN ---------

  // find which map contains *cfilt.
  for(imap=0; imap < NMAP_NONLIN; imap++ ) {
    ptrFilters = NONLIN_MAP[imap].FILTERS ;

    if ( strstr(ptrFilters,cfilt) != NULL )  {
      sprintf(msg,"%s: band=%s imap=%d Flux = %le", 
	      fnam, cfilt, imap, flux);

      f_scale = interp_1DFUN(OPT_INTERP, log10_flux,
			     NONLIN_MAP[imap].MAPSIZE,
			     NONLIN_MAP[imap].MAPVAL[0], // log10(flux)
			     NONLIN_MAP[imap].MAPVAL[1], // F_scale
			     msg ) ;
      return f_scale;
    }
  }  // end imap loop

  return f_scale ;
  
} // end get_flux_scale_NONLIN
