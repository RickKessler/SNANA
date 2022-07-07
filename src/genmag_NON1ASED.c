/***************************************************
  Created Aug 31, 2010 by R.Kessler

  
  Use non-Ia SED to generate observer-frame mags
  (SALT2-like model) without using K-corrections.

    HISTORY
   ~~~~~~~~~

  Mar 14 2016: include prep_NON1ASED and pick_NON1ASED here,
               along with include snlc_sim.h .
               Need more refactoring to define NON1ASED structure
               in genmag_NON1ASED.h so that snlc_sim.h reference
               can be removed.

  Aug 28 2017: new function getName_SED_FILE_NON1ASED()

  Jul 20 2018: 
    + replace interp_flux_SEDMODEL with get_flux_SEDMODEL to take
      care of late-time extrapolation. Expect no change here,
      but uses same utility as genmag_SIMSED.

  May 4 2019
    + read and apply index-dependent FLUXSCALE_NON1ASED (default SCALE=1)
      in NON1A.LIST file.

****************************************************/

#include <sys/stat.h>

#include "sntools.h" 
#include "sntools_cosmology.h" 
#include "fitsio.h"
#include "sntools_modelgrid.h"
#include "sntools_stronglens.h"
#include "snlc_sim.h"      // need the typedefs only

#include "MWgaldust.h"
#include "genmag_NON1ASED.h"
#include "genmag_SEDtools.h"

#define MODELGRID_GEN

// *************************************
void init_genmag_NON1ASED(int isparse, INPUTS_NON1ASED_DEF *INP_NON1ASED) {

  // Init SED for NON1ASED model.
  // Note that prep_NON1SED is called earlier.
  //
  // Inputs:
  //   isparse      = sparse index for template
  //   INP_NON1ASED = structure with all input information
  // 
  // May 15 2018: 
  //  + check option to call T0shiftExplode_SEDMODEL
  //  + new default is to call T0shiftPeak_SEDMODEL, so that t=0 at peak
  //
  // May 06 2019: refactor to pass INPUTS_NON1ASED_DEF
  //
  // Feb 06, 2020: 
  //   + for GENGRID option, call found_fluxerr_SEDMODEL()
  //     to check for FLUXERR column in SED file.

  double UVLAM     = INPUTS_SEDMODEL.UVLAM_EXTRAPFLUX ;
  int   DO_GENGRID = ( INP_NON1ASED->IFLAG_GEN == IFLAG_GENGRID ) ;
  int ifilt, ifilt_obs, NZBIN, NON1A_INDEX ;
  double Trange[2], Lrange[2] ;
  char sedcomment[40], *sedFile ;
  //  char fnam[] = "init_genmag_NON1ASED"  ;

  // ------------- BEGIN -------------

  if ( isparse < 0 ) {          // one-time init

    // summarize filter info
    filtdump_SEDMODEL();

    NLAMPOW_SEDMODEL      = 0;
    NZBIN                 = REDSHIFT_SEDMODEL.NZBIN ;
    SEDMODEL.NSURFACE     = 1 ;  // process 1 NONIA sed at a time
    SEDMODEL.RESTLAMMIN_FILTERCEN  = INP_NON1ASED->RESTLAMBDA_RANGE[0];
    SEDMODEL.RESTLAMMAX_FILTERCEN  = INP_NON1ASED->RESTLAMBDA_RANGE[1];
    SEDMODEL_MWEBV_LAST   = -999.  ;

    SEDMODEL.OPTMASK      = 
      OPTMASK_DAYLIST_SEDMODEL  +  // allow non-uniform day bins
      OPTMASK_T0SHIFT_PEAKMAG      // shift T=0 to be at peakmag
      ;

    // if error column exists in first SED file, 
    // add read-err flag to OPTMASK 
    if ( DO_GENGRID ) {
      sedFile     = INP_NON1ASED->SED_FILE[1] ;
      if ( found_fluxerr_SEDMODEL(sedFile) )  { 
	SEDMODEL.OPTMASK += OPTMASK_FLUXERR_SEDMODEL ; 
	printf("\t Read FLUXERR from SED files.\n");
      }
      else {
	printf("\t Do NOT read FLUXERR from SED files.\n");
      }
    }


    malloc_SEDFLUX_SEDMODEL(&TEMP_SEDMODEL,0,0,0);

    malloc_FLUXTABLE_SEDMODEL ( NFILT_SEDMODEL, NZBIN, NLAMPOW_SEDMODEL, 
				MXBIN_DAYSED_SEDMODEL, SEDMODEL.NSURFACE );

    // malloc DAY array (Aug 2017)
    int MEM = MXBIN_DAYSED_SEDMODEL * sizeof(double) ;
    SEDMODEL.DAY[ISED_NON1A] = (double*) malloc(MEM) ;

    return;
  }

  NON1A_INDEX = INP_NON1ASED->INDEX[isparse];
  sedFile     = INP_NON1ASED->SED_FILE[isparse] ;

  Trange[0] =  -150. ;  // widen Trange Apr 2 2018 
  Trange[1] =   500. ;  
  Lrange[0] = SEDMODEL.RESTLAMMIN_FILTERCEN ;
  Lrange[1] = SEDMODEL.RESTLAMMAX_FILTERCEN ;

  sprintf(sedcomment,"NON1A-%3.3d", NON1A_INDEX );

  rd_sedFlux(sedFile, sedcomment, Trange, Lrange
	     ,MXBIN_DAYSED_SEDMODEL, MXBIN_LAMSED_SEDMODEL
	     ,SEDMODEL.OPTMASK
	     ,&TEMP_SEDMODEL.NDAY, TEMP_SEDMODEL.DAY, &TEMP_SEDMODEL.DAYSTEP
	     ,&TEMP_SEDMODEL.NLAM, TEMP_SEDMODEL.LAM, &TEMP_SEDMODEL.LAMSTEP
	     ,TEMP_SEDMODEL.FLUX,  TEMP_SEDMODEL.FLUXERR );

  if ( UVLAM > 0.0 ) { UVLAM_EXTRAPFLUX_SEDMODEL(UVLAM, &TEMP_SEDMODEL); } 

  SEDMODEL.LAMSTEP[0]          = TEMP_SEDMODEL.LAMSTEP ; // Nov 17 2016  
  SEDMODEL.LAMSTEP[ISED_NON1A] = TEMP_SEDMODEL.LAMSTEP ; // Nov 17 2016  
  SEDMODEL.FLUXSCALE           = INP_NON1ASED->FLUXSCALE[NON1A_INDEX] ;

  // make sure that DAY=0 at peak (Oct 30 2014)

  int OPTMASK_EXPLODE = INPUTS_SEDMODEL.OPTMASK_T0SHIFT_EXPLODE ;
  int OPTMASK_PEAKMAG = ( SEDMODEL.OPTMASK & OPTMASK_T0SHIFT_PEAKMAG ) ;  
  if ( OPTMASK_EXPLODE >= 0 ) {
    T0shiftExplode_SEDMODEL(OPTMASK_EXPLODE, &TEMP_SEDMODEL, 1); 
  }
  else if ( OPTMASK_PEAKMAG > 0 ) {
    T0shiftPeak_SEDMODEL(&TEMP_SEDMODEL,1); 
  }


  init_FINEBIN_SEDMODEL(-1);           // free FINEBIN memory alloc.
  init_FINEBIN_SEDMODEL(ISED_NON1A);   // allocate memory

  zero_flux_SEDMODEL();
  for(ifilt=1; ifilt <= NFILT_SEDMODEL; ifilt++) {
    ifilt_obs = FILTER_SEDMODEL[ifilt].ifilt_obs ;
    init_flux_SEDMODEL(ifilt_obs,ISED_NON1A);
  }  

  init_flux_SEDMODEL(0,0); // flag to print SED-summary

  return;

} // end of init_genmag_NON1ASED


// *****************************************
void genmag_NON1ASED (
		   int index             // (I) NON1A class
		   ,int ifilt_obs        // (I) obs filter index
		   ,double mwebv         // (I) Galactic extinction: E(B-V)
		   ,double z             // (I) Supernova redshift
		   ,double mu            // (I) distance modulus
		   ,double RV_host       // (I) RV of host galaxy
		   ,double AV_host       // (I) AV of host galaxy
		   ,int     Nobs         // (I) number of epochs
		   ,double *Tobs_list    // (I) Tobs-Tpeak (days)
		   ,double *magobs_list   // (O) obs mags
		   ,double *magerr_list   // (O) obs mag-errs
		   ) {

  // Oct 23 2018: replace x0 argument with mu

  int  ifilt, epobs, ILAMPOW = 0 ;
  double  z1, ZP, meanlam_obs, meanlam_rest, Tobs, Trest ;
  double  AV_MW, XT_MW, XT_HOST, flux, FLUX, magerr, magobs;
  char *cfilt;
  //  char fnam[] = "genmag_NON1ASED" ;

  // --------- BEGIN ----------

  ifilt = IFILTMAP_SEDMODEL[ifilt_obs] ;
  z1    = 1. + z;

  // filter info for this "ifilt"
  meanlam_obs  = FILTER_SEDMODEL[ifilt].mean ;  // mean lambda
  meanlam_rest = meanlam_obs/z1 ;
  ZP           = FILTER_SEDMODEL[ifilt].ZP ;
  cfilt        = FILTER_SEDMODEL[ifilt].name ;
  
  // get approx Galactic extinction using central wavelength of filter
  // No need to be as precise as SALT2.
  AV_MW = RV_MWDUST * mwebv ; 
  XT_MW = GALextinct ( RV_MWDUST, AV_MW, meanlam_obs, 94 );

  // get extinction from host in rest-frame
  if ( AV_host > 1.0E-9 ) 
    { XT_HOST = GALextinct ( RV_host, AV_host, meanlam_rest, 94 ); }
  else
    { XT_HOST = 0.0 ; }

  // Nov 16 2016: load MWEBV tables for each band and spectrum
  fill_TABLE_MWXT_SEDMODEL(RV_MWDUST, mwebv); 

  // - - - - - - - 
  for ( epobs=0; epobs < Nobs; epobs++ ) {

    Tobs = Tobs_list[epobs];
    Trest = Tobs / z1;

    magobs_list[epobs] = MAG_ZEROFLUX ;
    magerr_list[epobs] = MAGERR_UNDEFINED ;

    if ( Trest < TEMP_SEDMODEL.DAYMIN ) { continue ; }
    
    // call function to handle Trest extrap if needed
    flux    = get_flux_SEDMODEL(ISED_NON1A, ILAMPOW, ifilt_obs, z, Trest) ;
    magerr  = get_magerr_SEDMODEL(ISED_NON1A, ifilt_obs, z, Trest);
    FLUX    = flux * pow(TEN,-0.4*mu);

    // if(epobs==10) { printf(" xxx flux=%le magerr=%f \n", flux, magerr); }

    if ( flux == FLUX_UNDEFINED ) {
      FLUX = FLUX_UNDEFINED ;
      magobs_list[epobs] = MAG_UNDEFINED  ;
      magerr_list[epobs] = MAGERR_UNDEFINED ;
    }
    else if ( FLUX > 1.0E-30 ) {
      magobs =  (ZP + XT_MW + XT_HOST) - 2.5*log10(FLUX);
      magobs_list[epobs] = magobs ; 
      magerr_list[epobs] = magerr ;
    }
    else {
      magobs_list[epobs] = MAG_ZEROFLUX ;
      magerr_list[epobs] = MAGERR_UNDEFINED ;
    }
    
    if ( fabs(Trest) < 5 && ifilt_obs < -8 ) 
      printf(" xxxxx Trest(%s)=%6.1f  ep=%2d  FLUX=%9.3le  mag=%6.3f  \n", 
	     cfilt, Trest, epobs, FLUX, magobs_list[epobs] );
    
  } // epobs

  return ;

} // end of genmag_NON1ASED


// ********************************************
void prep_NON1ASED(INPUTS_NON1ASED_DEF *INP_NON1ASED, 
		   GENLC_NON1ASED_DEF *GEN_NON1ASED) {

  // Feb 11, 2009 R.Kessler
  // One-time call to prepare picking a NON1ASED template
  // If NGENTOT > NSED, then distribute according to weight.
  // If NGENTOT < NSED, then generate 1 each for the first
  // NGENTOT SEDs, then none for the rest.
  //
  // Inputs:
  //  *INP_NON1ASED is the input structure (from reading input files)
  //  *GEN_NON1ASED is used in generation
  //
  // Mar 14 2016: refactor with all _non1a stuff removed.
  // Mar 25 2016: refactor again using INP_NON1ASED and GEN_NON1ASED
  // Aug 15 2016: call read_NON1A_LIST() 
  // Feb 01 2017: check ONEFLAG to handle NGEN < N(SED)
  // Jul 29 2017: allow sum(NON1A wgt)=0 as long as sum(PEC1A wgt)>0
  // May 16 2019: set NINDEX after read_NON1A_LIST in case it changes

  float FRAC_PEC1A = GEN_NON1ASED->FRAC_PEC1A ; // Ngen(pecIa)/NgenTot
  int   DO_GENGRID = ( GEN_NON1ASED->IFLAG_GEN == IFLAG_GENGRID ) ;

  int index, isp, NINDEX, ISPEC1A ;
  int NGENTMP, NGENSUM, NGENTOT, NSED_SKIP, ONEFLAG ;

  double  wgt, WGTSUM_TOT, WGTSUM[2], XNGEN[2], XN ;

  char  name[MXPATHLEN],  *ptrSed ;
  char  fnam[] = "prep_NON1ASED" ;
    
  // ------------ BEGIN ------------

  INP_NON1ASED->RESTLAMBDA_RANGE[0] = LAMMIN_SEDMODEL ;
  INP_NON1ASED->RESTLAMBDA_RANGE[1] = LAMMAX_SEDMODEL ;

  // first check path to NON1ASEDs
  if ( strlen(INP_NON1ASED->PATH) == 0 ) 
    { sprintf(INP_NON1ASED->PATH,"%s/snsed/NON1A", PATH_SNDATA_ROOT);  }

  NINDEX = INP_NON1ASED->NINDEX ;
  if ( NINDEX == 0 ) {
    sprintf(c1err,"User specified 'GENMODEL: %s'", MODELNAME_NON1ASED );
    sprintf(c2err,"but found no 'NON1A: <index> <wgt>' ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  sprintf(BANNER, " PREPARE %s TEMPLATE SELECTION ", MODELNAME_NON1ASED );
  print_banner ( BANNER );
  printf("    MODELPATH = %s \n", INP_NON1ASED->PATH );

  // init NON1A variables
  for ( isp=0; isp < MXNON1A_TYPE; isp++ ) {
    INP_NON1ASED->FLUXSCALE[isp]   = 1.0 ;
    INP_NON1ASED->INDEXVALID[isp]  = 0;
    GEN_NON1ASED->NGENWR[isp]      = 0; 
    GEN_NON1ASED->NGENTOT[isp]     = 0; 
    GEN_NON1ASED->CIDRANGE[isp][0] = INP_NON1ASED->CIDOFF+1; 
    GEN_NON1ASED->CIDRANGE[isp][1] = -9; 
    sprintf(GEN_NON1ASED->TYPE[isp],  "NULL" );
  }


  // -------------------------------------
  // read master list of available NON1A
  read_NON1A_LIST(INP_NON1ASED ); // load NON1ASED->LIST_XXX

  // ---------------------------
  // Aug 15 2016:
  // sort user-input NON1A & PEC1A so that PEC1A are at end of list
  sort_NON1ASED(INP_NON1ASED);

  NINDEX = INP_NON1ASED->NINDEX ; // may have changed after read_NON1A_LIST
			 
  // ---------------------------
  // get sum of wgts to renormalize to wgt-sum = 1
  WGTSUM[0] = WGTSUM[1] = WGTSUM_TOT = 0.0 ;
  for ( isp=1; isp <= NINDEX; isp++ ) {   
    wgt     = INP_NON1ASED->WGT[isp];
    ISPEC1A = INP_NON1ASED->ISPEC1A[isp];
    WGTSUM[ISPEC1A] += wgt ;
    WGTSUM_TOT      += wgt ;
  }

  if ( WGTSUM_TOT <= 0.0 ) {
    sprintf(c1err,"Sum of NON1A+PEC1A weights = %f", WGTSUM_TOT );
    sprintf(c2err,"Check NON1A keywords in sim-input file.");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
  }

  // determine first NGEN for each NON1A class
  NGENSUM = NSED_SKIP = 0;
  NGENTOT = INP_NON1ASED->NGENTOT ; 

  XN       = (double)NGENTOT ;
  XNGEN[0] = XN * ( 1.0 - FRAC_PEC1A );
  XNGEN[1] = XN * FRAC_PEC1A ;

  // check flag to process only 1 per SED (Feb 2017)
  ONEFLAG = ( NGENTOT < NINDEX ) ;

  //loop over sparse index
  for ( isp=1; isp <= NINDEX; isp++ ) {
    index   = INP_NON1ASED->INDEX[isp] ;
    wgt     = INP_NON1ASED->WGT[isp];
    ISPEC1A = INP_NON1ASED->ISPEC1A[isp];
    XN      = ( wgt / WGTSUM[ISPEC1A] ) * XNGEN[ISPEC1A] ;

    if ( !DO_GENGRID && (ISPEC1A && FRAC_PEC1A == 0.0)  ) {
      print_preAbort_banner(fnam);
      printf("\t Check that DNDZ_PEC1A key is defined.\n");
      printf("\t Check that GENRANGE_REDSHIFT is not a delta function.\n");
      printf("\t Check that GENRANGE_PEAKMJD  is not a delta function.\n");
      printf("\t Check that SOLID_ANGLE > 0 \n");
      sprintf(c1err,"Generating %s (%d)",
	      INP_NON1ASED->LIST_TYPE[index], index);
      sprintf(c2err,"But PEC1A rate=0. Check keys above.");
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    
    /*
    printf(" xxx isp=%2d ISPEC1A=%d  wgt=%.3f WGTSUM=%.3f XN=%f \n",
	   isp, ISPEC1A, wgt, WGTSUM[ISPEC1A], XNGEN[ISPEC1A] ); 
	   xxxxxxxxx */

    sprintf(name,"%s", INP_NON1ASED->LIST_NAME[index] ) ;
    NGENTMP = (int)(XN + 0.5);

    // check one-per-SED option when NGEN is <= Number of SEDs
    if ( ONEFLAG ) {
      if ( isp <= NGENTOT ) { NGENTMP = 1; }
      else                  { NGENTMP = 0; }
    }

    // adjust last class so that NGENTMP sums to NGEN
    if ( isp == NINDEX  || NGENTMP > NGENTOT ) 
      { NGENTMP = NGENTOT - NGENSUM ;  }

#ifdef MODELGRID_GEN
    if ( DO_GENGRID ) { 
      NGENTMP  = 1;
      NGENTMP *= SNGRID_WRITE.NBIN[IPAR_GRIDGEN_LOGZ] ;
      NGENTMP *= SNGRID_WRITE.NBIN[IPAR_GRIDGEN_COLORPAR] ;
      NGENTMP *= SNGRID_WRITE.NBIN[IPAR_GRIDGEN_COLORLAW] ;
    }
#endif

    if ( INP_NON1ASED->INDEXVALID[index] == 0 ) {
      sprintf(c1err,"Invalid NON1A index = %d in sim-input file.",index );
      sprintf(c2err,"Compare with %s", INP_NON1ASED->LISTFILE );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
    }

    if ( NGENTMP < 0 ) { NGENTMP = 0; } // Jan 9 2012

    NGENSUM += NGENTMP ;
    INP_NON1ASED->NGEN[isp]  = NGENTMP ;
    INP_NON1ASED->MXGEN[isp] = NGENSUM ;

    ptrSed  = INP_NON1ASED->SED_FILE[isp] ;
    getName_SED_FILE_NON1ASED(INP_NON1ASED->PATH, name, 
			      ptrSed); // return ptrSed

    if ( NGENTMP > 0 ) {
      printf("     Will generate %4d %s (MXGEN=%4d) with %s\n", 
	     NGENTMP, MODELNAME_NON1ASED, NGENSUM, name );
      fflush(stdout);
    }
    else {
      NSED_SKIP++ ; 
    }


    if ( NGENTMP < 0 ) {
      sprintf(c1err,"Cannot generate %d %s events.", 
	      NGENTMP, MODELNAME_NON1ASED );
      sprintf(c2err,"for %s", ptrSed );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);       
    }

    // check that NONIA and PECIA are not mixed up (Aug 2016)
    int ISPEC1A_USER = INP_NON1ASED->ISPEC1A[isp] ;
    int ISPEC1A_LIST = INP_NON1ASED->LIST_ISPEC1A[index] ;
    if ( ISPEC1A_USER != ISPEC1A_LIST ) {
      sprintf(c1err, "ISPEC1A(USER)=%d but ISPEC1A(NON1A_LIST)=%d",
	      ISPEC1A_USER, ISPEC1A_LIST );
      sprintf(c2err,"for %s (index=%d)", name, index );
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);         
    }

  } // isp


  if ( NSED_SKIP ) {
    printf("     WARNING: skipped %d of %d SEDs \n",
	   NSED_SKIP, NINDEX); 
  }

  fflush(stdout);

  return ;

} // end of prep_NON1ASED


// =========================================================
void  getName_SED_FILE_NON1ASED(char *PATH, char *inpName, char *outName) {

  // Created Aug 28 2017   
  // Inputs
  //   *PATH     is the path containing NON1A.LIST file
  //   *inpName  is a name in the NON1A.LIST file
  //
  // Output 
  //   *outName is the full file name, including path.

  struct stat statbuff;
  int jstat ;
  //  char fnam[] = "getName_SED_FILE_NON1ASED" ;

  // -------------- BEGIN --------------

  // First check legacy option where outName = PATH/inpName.SED
  sprintf(outName, "%s/%s.SED", PATH, inpName ); 
  jstat = stat(outName, &statbuff); // returns 0 if file exists
  if ( jstat == 0 ) { return ; }

  // if we get here, then inpName is assumed to be a complete
  // file name. If it has a slash (/), then it included a path.
  // If no slash, then append PATH.

  
  if ( strstr(inpName, "/") ) { 
    // user gave full path and name
    sprintf(outName, "%s", inpName );
  }
  else {
    // user gives file name without PATH, so assume same path as NON1A.LIST
    sprintf(outName,"%s/%s", PATH, inpName ); 
  } 

  return;

} // end getName_SED_FILE_NON1ASED


// =========================================================
int count_NON1A_LIST(char *PATHMODEL) {

  // count number of keys (i.e., templates) defined in NON1A.LIST file.
  // Called by init_SNgrid().

  int  Ncount = 0 ;
  char listFile[MXPATHLEN], cget[200];
  FILE *fp;
  char fnam[] = "count_NON1A_LIST" ;

  // --------------- BEGIN ----------------

  sprintf(listFile, "%s/NON1A.LIST",  PATHMODEL );
  if ( (fp = fopen(listFile, "rt"))==NULL ) { 
    sprintf ( c1err, "Cannot open file :" );
    sprintf ( c2err," '%s' ", listFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  while( (fscanf(fp, "%s", cget)) != EOF) {
    if ( ISKEY_NON1A_LIST(cget) ) { Ncount++ ; }
  } 

  fclose(fp);

  return(Ncount);

} // end count_NON1A_LIST

// =========================================================
void read_NON1A_LIST(INPUTS_NON1ASED_DEF *INP_NON1ASED ) {
                        
  // Created Aug 15 2016
  // Read NON1A.LIST file.
  // Syntax:
  //
  // FLUX_SCALE: 1.0                # optional
  // RESTLAMBDA_RANGE: 2000 12000   # optional
  // NON1A: <ID1> <TYPE1> <fileName1>
  // NON1A: <ID2> <TYPE2> <fileName2>
  // NON1A: <ID3> <TYPE3> <fileName3>
  // etc ...
  // where ID=integer id for TYPE, and TYPE is a string
  // (e.g., Ia, II, Ibc, etc)
  //
  // Jan 19 2017: allow NON1ASED key as well as NON1A.
  // Jan 30 2017: skip comment lines
  // Apr 28 2019: 
  //   + read optional FLUX_SCALE and RESTLAMBDA_RANGE keys
  //     FLUX_SCALE can repeat within NON1A.LIST file so that different
  //     model groups have different FLUX_SCALE.
  //
  FILE *fp;
  int NLIST, L_NON1A, L_PEC1A, index, NINDEX ;
  int FOUND_PEC1A=0;
  double SCALE = 1.0 ;
  int    ALLNON1A_FLAG=0, SNTAG_ALLNON1A=0 ;
  double MAGOFF_ALLNON1A=0.0, MAGSMEAR_ALLNON1A=0.0;
  char cget[80], type[20], fname[MXPATHLEN], listFile[MXPATHLEN] ;
  char tmpLine[100], *ptrTmp = tmpLine ;
  char fnam[] = "read_NON1A_LIST" ;

  // ------------- BEGIN -----------

  sprintf(listFile, "%s/NON1A.LIST",  INP_NON1ASED->PATH );
  sprintf(INP_NON1ASED->LISTFILE, "%s", listFile);

  if ( (fp = fopen(listFile, "rt"))==NULL ) { 
    sprintf ( c1err, "Cannot open file :" );
    sprintf ( c2err," '%s' ", listFile );
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }


  printf("    Read all NON1A from '%s' \n", listFile);

  NLIST = INP_NON1ASED->NLIST = 0 ;

  // If the first NON1A class is zero, then use all NON1A with equal prob.
  if ( INP_NON1ASED->INDEX[1] == 0 ) {  
    ALLNON1A_FLAG = 1;   INP_NON1ASED->NINDEX = 0;  
    MAGOFF_ALLNON1A   = INP_NON1ASED->MAGOFF[1];
    MAGSMEAR_ALLNON1A = INP_NON1ASED->MAGSMEAR[1][0];
    SNTAG_ALLNON1A    = INP_NON1ASED->SNTAG[1];
  }
  else { 
    ALLNON1A_FLAG = 0; 
  }

  while( (fscanf(fp, "%s", cget)) != EOF) {

    // if comment key is found, read remainder of line into dummy string 
    // so that anything after comment key is ignored (even a valid key) 
    if ( cget[0] == '#' || cget[0] == '!' || cget[0] == '%' )
      { ptrTmp = fgets(tmpLine, 80, fp) ; continue ; }

    // 4.2019: check keys that are the same as in SED.INFO file ...
    // later should just read SED.INFO file if it's there.
    if ( strcmp(cget,"FLUX_SCALE:") == 0 ) // 
      { readdouble(fp, 1, &SCALE ); }  // applied below
    if ( strcmp(cget,"RESTLAMBDA_RANGE:") == 0 ) 
      { readdouble(fp, 2, INP_NON1ASED->RESTLAMBDA_RANGE); }


    L_NON1A = ISKEY_NON1A_LIST(cget) ;
    L_PEC1A = ( strcmp(cget,"PEC1A:") == 0 || strcmp(cget,"PECIA:") == 0 );

    if ( L_NON1A || L_PEC1A ) {

      readint(fp,  1, &index );
      readchar(fp, type );
      readchar(fp, fname );  ENVreplace(fname,fnam,1);

      if ( index <= 0 || index >= MXNON1A_TYPE ) {
	sprintf(c1err,"Invalid %s index = %d read from ", 
		MODELNAME_NON1ASED, index );
	sprintf(c2err,"%s", INP_NON1ASED->LISTFILE );
	errmsg(SEV_FATAL, 0, fnam, c1err, c2err ); 
      }

      NLIST++ ;
      INP_NON1ASED->LIST_INDEX[NLIST] = index;

      // store inputs bases on index rather than sparse-index
      INP_NON1ASED->INDEXVALID[index] = 1; 
      INP_NON1ASED->FLUXSCALE[index]  = SCALE ;
      sprintf(INP_NON1ASED->LIST_TYPE[index], "%s", type );
      sprintf(INP_NON1ASED->LIST_NAME[index], "%s", fname );

      // check for peculiar SN1A (Aug 2016)
      INP_NON1ASED->LIST_ISPEC1A[index] = 0 ;
      if ( L_PEC1A ) { 
	INP_NON1ASED->LIST_ISPEC1A[index] = 1; 
	sprintf(INP_NON1ASED->LIST_TYPE[index], "PEC1A-%s", type );
	FOUND_PEC1A=1;
      }

      // check option to use all NON1A with equal prob.
      if ( ALLNON1A_FLAG ) {
	INP_NON1ASED->NINDEX++ ; NINDEX = INP_NON1ASED->NINDEX;
	INP_NON1ASED->INDEX[NINDEX]        = index ; 
	INP_NON1ASED->WGT[NINDEX]          = 1.0 ; 
	INP_NON1ASED->KEYVAL[NINDEX][2]    = 1.0 ; // 2nd key is WGT
	INP_NON1ASED->MAGOFF[NINDEX]       = MAGOFF_ALLNON1A ;
	INP_NON1ASED->MAGSMEAR[NINDEX][0]  = MAGSMEAR_ALLNON1A ;
	INP_NON1ASED->SNTAG[NINDEX]        = SNTAG_ALLNON1A ;
	sprintf(INP_NON1ASED->SED_FILE[NINDEX],"%s", fname ) ;
      }

    }
  } // end of fscanf loop
  INP_NON1ASED->NLIST = NLIST ;

  fclose(fp);

  return ;

} // end read_NON1A_LIST

int ISKEY_NON1A_LIST(char *string) {
  if ( strcmp(string,"NON1A:")    == 0  ) { return(1); }
  if ( strcmp(string,"NONIA:")    == 0  ) { return(1); }
  if ( strcmp(string,"NON1ASED:") == 0  ) { return(1); }
  if ( strcmp(string,"NONIASED:") == 0  ) { return(1); }
  return(0);
}

// ======================================================
void sort_NON1ASED(INPUTS_NON1ASED_DEF *INP_NON1ASED) {

  // Created Aug 2016
  // sort user NON1ASED list  so that all PEC1A are at
  // the end of the list. This sorting doesn't matter
  // for the NON1ASED model, but it does matter for
  // writing NON1A to a grid and then using NON1AGRID.
  // The sorting for NON1AGRID allows generating redshift
  // with the correct (NON1A vs. PEC1A) rate-vs-z model.

  INPUTS_NON1ASED_DEF *TMP_NON1ASED ;
  int isp, inew, ISPEC1A, NNON1A=0, NPEC1A=0 ;
  int  NINDEX = INP_NON1ASED->NINDEX ;
  //  char fnam[] = "sort_NON1ASED" ;

  // --------------- BEGIN ----------------

  // check number of PEC1A keys
  for ( isp=1; isp <= NINDEX; isp++ )
    { if ( INP_NON1ASED->ISPEC1A[isp] ) { NPEC1A++ ; } }

  if ( NPEC1A == 0 ) { return ; }

  // -------------------------------
  // allocate temp structure
  TMP_NON1ASED = (INPUTS_NON1ASED_DEF*) malloc( sizeof(INPUTS_NON1ASED_DEF));

  // copy NON1ASED arrays into TMP array
  for ( isp=1; isp <= NINDEX; isp++ )
    { copy_NON1ASED(isp,isp, INP_NON1ASED, TMP_NON1ASED) ;  }


  int ISPOFF_PEC1A = NINDEX-NPEC1A ;
  NNON1A = NPEC1A = 0;
  for ( isp=1; isp <= NINDEX; isp++ ) {
    
    ISPEC1A = INP_NON1ASED->ISPEC1A[isp] ;
    if ( ISPEC1A ) 
      { NPEC1A++ ; inew = ISPOFF_PEC1A + NPEC1A; }
    else
      { NNON1A++ ; inew = NNON1A; }

    copy_NON1ASED(isp,inew, TMP_NON1ASED, INP_NON1ASED) ; 
  }
  
  free(TMP_NON1ASED);

  return;

} // end sort_NON1ASED

// =================================================
void copy_NON1ASED(int i1, int i2, 
		   INPUTS_NON1ASED_DEF *NON1ASED1, 
		   INPUTS_NON1ASED_DEF *NON1ASED2 ) {

  //  char fnam[] = "copy_NON1ASED" ;

  // ----------- BEGIN -------------

  NON1ASED2->INDEX[i2]        = NON1ASED1->INDEX[i1] ;
  NON1ASED2->WGT[i2]          = NON1ASED1->WGT[i1] ;
  NON1ASED2->MAGOFF[i2]       = NON1ASED1->MAGOFF[i1] ;
  NON1ASED2->MAGSMEAR[i2][0]  = NON1ASED1->MAGSMEAR[i1][0] ;
  NON1ASED2->MAGSMEAR[i2][1]  = NON1ASED1->MAGSMEAR[i1][1] ;
  NON1ASED2->SNTAG[i2]        = NON1ASED1->SNTAG[i1] ;
  NON1ASED2->ISPEC1A[i2]      = NON1ASED1->ISPEC1A[i1] ;
  NON1ASED2->FLUXSCALE[i2]    = NON1ASED1->FLUXSCALE[i1] ;

  return ;

} // end copy_NON1ASED

