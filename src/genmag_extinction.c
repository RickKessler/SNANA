/*********************************************************
  Created Nov 2022

  Compute/store/fetch host-galaxy extinction (XTMAG) in rest-frame bands
  for SNIa models that require explicit K-corrections: mlcs2k2, snoopy.
  XTMAG is stored on a GRIDMAP vs. {1/RV, Trest, IFILT_REF}
  for fast lookup in sim and LCfit programs. 
  
  The SNSED flux (read from kcor/calib file) is used to perform full
  integral over each band to get a more accurate XTMAG compared
  to using only central wavelength of the band.

  Translated from fortran -> C code as part of K-cor refactor to
  translate all fortran utilities into C.

 ********************************************************/

#include "fitsio.h"
#include "MWgaldust.h"
#include "sntools.h"
#include "sntools_calib.h"
#include "genmag_extinction.h"

// ================================
void init_genmag_extinction(int OPT_SNXT) { 

  // Initialize XTMAG on 3D grid of 1/RV, Trest, IFILTDEF_REST.
  // Roughly translated from original fortran subroutine INIT_XTHOST.

  double Trest ;
  char fnam[] = "init_genmag_extinction" ;

  // ----------- BEGIN ------------

  print_banner(fnam);

  XTMAG_INFO.OPT_SNXT = OPT_SNXT;


  if ( OPT_SNXT == OPT_SNXT_CCM89 ) {
    fill_binInfo_XTMAG();
    fill_GRIDMAP3D_XTMAG();
  }
  else if ( OPT_SNXT == OPT_SNXT_SJPAR ) {
    sprintf(c1err,"OPT_SNXT_SJPAR not yet available in refactored C code.");
    sprintf(c2err,"Need to translate fortran RDXTPAR");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
    
  }
  else {
    sprintf(c1err,"Invalid OPT_SNXT=%d", OPT_SNXT);
    sprintf(c2err,"Valid options are %d(CCM89) or %d(SJPAR)",
	    OPT_SNXT_CCM89, OPT_SNXT_SJPAR);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

    // print human-readable tables to stdout
    Trest =  0.0 ; dump_XTMAG(Trest); 
    Trest = 20.0 ; dump_XTMAG(Trest); 

  return;

} // end init_genmag_extinction

void init_genmag_extinction__(int *OPT_SNXT) 
{ init_genmag_extinction(*OPT_SNXT) ; }


// =====================================
void fill_binInfo_XTMAG(void) {

  // fill binInfo structure for Trest and RVinv
  
  KCOR_BININFO_DEF  *BININFO;
  double RVinv;
  int i, NBIN;
  char fnam[] = "fill_binInfo_XTMAG" ;

  // ------------- BEGIN -----------

  // load RVinv grid using hard-wired params from .h file
  BININFO = &XTMAG_INFO.BININFO_RVinv ;
  NBIN    = (int)((RVinv_MAX_XTMAG - RVinv_MIN_XTMAG) / RVinv_BIN_XTMAG );
  BININFO->GRIDVAL = (double*) malloc( (NBIN+1)*sizeof(double) );
  NBIN=0;
  for(RVinv=RVinv_MIN_XTMAG; RVinv <= RVinv_MAX_XTMAG; RVinv+=RVinv_BIN_XTMAG )
    { BININFO->GRIDVAL[NBIN] = RVinv ;    NBIN++ ;   }
  BININFO->NBIN     = NBIN ;
  BININFO->RANGE[0] = RVinv_MIN_XTMAG;
  BININFO->RANGE[1] = RVinv_MAX_XTMAG;

  // copy Trest GRID from KCOR structure
  BININFO = &XTMAG_INFO.BININFO_Trest ;
  NBIN    = CALIB_INFO.BININFO_T.NBIN ;
  BININFO->GRIDVAL = (double*) malloc( NBIN * sizeof(double) ) ;
  BININFO->NBIN = NBIN;
  BININFO->RANGE[0] = CALIB_INFO.BININFO_T.RANGE[0];
  BININFO->RANGE[1] = CALIB_INFO.BININFO_T.RANGE[1];
  for(i=0; i < NBIN; i++ ) 
    { BININFO->GRIDVAL[i] = CALIB_INFO.BININFO_T.GRIDVAL[i] ;    }

  return;

} // end fill_binInfo_XTMAG


// ==================================
void fill_GRIDMAP3D_XTMAG(void) {

  FILTERCAL_DEF *FILTERCAL_REST = &CALIB_INFO.FILTERCAL_REST;
  int NFILTDEF_REST             = FILTERCAL_REST->NFILTDEF ;
  int NBIN_Trest                = XTMAG_INFO.BININFO_Trest.NBIN;
  int NBIN_RVinv                = XTMAG_INFO.BININFO_RVinv.NBIN;
  int ifilt, J1D, iT, iRv ;
  double Trest, RVinv, XTMAG ;
  float temp_mem;
  char MAPNAME[] = "HostExtinction" ;
  char fnam[] = "fill_GRIDMAP3D_XTMAG" ;

  // ------------ BEGIN ------------

  if ( NFILTDEF_REST <= 0 ) {
    sprintf(c1err,"Invalid NFILTDEF_REST = %d ", NFILTDEF_REST);
    sprintf(c2err,"but rest-frame bands are expected");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );
  }

  printf("\t %s(ifilt_rest, Trest, RVinv) \n", fnam);

  int NBIN_TOT = NBIN_RVinv * NBIN_Trest * NFILTDEF_REST;
  int NDIM_INP = 3;
  int NDIM_FUN = 1;
  temp_mem = malloc_double2D(+1,NDIM_INP+NDIM_FUN,NBIN_TOT,&TEMP_XTMAG_ARRAY);
  J1D = 0 ;
  // start the 3D loop
  
  for(ifilt=0; ifilt < NFILTDEF_REST; ifilt++ ) {
    for(iT=0; iT < NBIN_Trest; iT++ ) {
      for(iRv=0; iRv < NBIN_RVinv; iRv++ ) {     

	RVinv    = XTMAG_INFO.BININFO_RVinv.GRIDVAL[iRv];
	Trest    = XTMAG_INFO.BININFO_Trest.GRIDVAL[iT];

	XTMAG = eval_XTMAG_AV1(ifilt, Trest, RVinv);

	TEMP_XTMAG_ARRAY[0][J1D]  = RVinv;
	TEMP_XTMAG_ARRAY[1][J1D]  = Trest;
	TEMP_XTMAG_ARRAY[2][J1D]  = (double)ifilt;
	TEMP_XTMAG_ARRAY[3][J1D]  = XTMAG;

	J1D++ ;

      } // end iRVinv loop
    }  // end iTrest loop   
  } // end ifilt

  init_interp_GRIDMAP(IDGRIDMAP_XTMAG, MAPNAME, 
		      NBIN_TOT, NDIM_INP, NDIM_FUN, OPT_EXTRAP_KCOR, 
		      TEMP_XTMAG_ARRAY, &TEMP_XTMAG_ARRAY[NDIM_INP],
		      &XTMAG_INFO.GRIDMAP3D); // <== returned

  printf("\t Allocate %.1f/%.1f MB of GRIDMAP/temp memory for %s\n", 
	 XTMAG_INFO.GRIDMAP3D.MEMORY, temp_mem, MAPNAME);
  fflush(stdout);
 
  // free TEMP_XTMAG_ARRAY
  malloc_double2D(-1, NDIM_INP+NDIM_FUN, NBIN_TOT, &TEMP_XTMAG_ARRAY);

  return;

} // end fill_GRIDMAP3D_XTMAG



// ================================================
double eval_XTMAG_AV1(int ifilt, double Trest, double RVinv) {

  // Return host extinction for sparse band index ifilt,
  // rest-frame phase "Trest" and RVinv = 1/RV/
  // Computation is for AV=1.

  FILTERCAL_DEF *FILTERCAL_REST = &CALIB_INFO.FILTERCAL_REST ;
  
  int    ifiltdef  = FILTERCAL_REST->IFILTDEF[ifilt];
  int    NBL_FILT  = FILTERCAL_REST->NBIN_LAM[ifilt] ;

  int    NBL_SED   = CALIB_INFO.BININFO_LAM.NBIN;
  int    NBT_SED   = CALIB_INFO.BININFO_T.NBIN;
  double Trest_MIN = CALIB_INFO.BININFO_T.RANGE[0];
  double Trest_MAX = CALIB_INFO.BININFO_T.RANGE[1];
  double Trest_BIN = CALIB_INFO.BININFO_T.BINSIZE ;

  double AV1    = 1.0 ;
  double RV     = 1.0/RVinv ;

  int  ilam_filt, ilam_sed, it, jflux ;
  double XTMAG_AV1, XTMAG_TMP, LAM, TRANS, PARDUM=0.0 ;
  double FLUX, FT, arg, FLUX_RATIO, XTFRAC, FLUXSUM=0.0, FLUXSUMXT=0.0 ;
  char fnam[] = "eval_XTMAG_AV1" ;

  // ---------- BEGIN ----------

  XTMAG_AV1 = 0.0 ; // init output

  // loop over wave bins for filter
  for(ilam_filt=0; ilam_filt < NBL_FILT; ilam_filt++ ) {
    LAM       = FILTERCAL_REST->LAM[ifilt][ilam_filt];
    TRANS     = FILTERCAL_REST->TRANS[ifilt][ilam_filt];

    XTMAG_TMP = GALextinct ( RV, AV1, LAM, OPT_MWCOLORLAW_XTMAG, &PARDUM );

    // get index for SN flux
    ilam_sed  = FILTERCAL_REST->ILAM_SED[ifilt][ilam_filt];
    it        = (int)((Trest - Trest_MIN)/Trest_BIN);
    jflux     = NBL_SED*it + ilam_sed ;
    FLUX      = (double)CALIB_INFO.FLUX_SNSED_F[jflux] ;

    FT        = FLUX * TRANS ;
    arg       = -0.4*XTMAG_TMP;
    XTFRAC    = pow(10.0,arg);

    FLUXSUM   += FT ;
    FLUXSUMXT += (FT * XTFRAC) ;
  }

  if ( FLUXSUM < 1.0E-40 || FLUXSUMXT < 1.0E-40 ) {
    sprintf(c1err,"Invalid FLUXSUM=%le / FLUXSUMXT=%le \n",
	    FLUXSUM, FLUXSUMXT );
    sprintf(c2err,"FLUXSUM[XT] should be > 0");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err );

  }

  FLUX_RATIO   = FLUXSUMXT / FLUXSUM ;
  XTMAG_AV1    = -2.5*log10(FLUX_RATIO) ;

  return XTMAG_AV1 ;

} // end eval_XTMAG_AV1


// ==================================
void dump_XTMAG(double Trest) {

  FILTERCAL_DEF *FILTERCAL_REST = &CALIB_INFO.FILTERCAL_REST;
  int NFILTDEF_REST             = FILTERCAL_REST->NFILTDEF ;

#define NBIN_DUMP_RVinv 5
  double RVinv_LIST[NBIN_DUMP_RVinv] = { 0.25, 0.50, 0.75, 1.0, 2.0 } ;

  int i, ifilt, ifiltdef, lenf ;
  double RVinv, XTMAG_AV1, RV, AV1=1.0 ;
  char   *name, band[4] ;
  char fnam[] = "dump_XTMAG" ;

  // ------------- BEGIN -------------

  printf("\n\t\t A_X/AV at Trest = %.1f with 1/RV= \n", Trest);
  printf("\t BAND  ");
  for(i=0; i < NBIN_DUMP_RVinv; i++ ) {
    RVinv = RVinv_LIST[i];
    printf("%5.3f  ", RVinv); 
  }
  printf("\n");
  printf("\t--------------------------------------------- \n");

  for(ifilt=0; ifilt < NFILTDEF_REST; ifilt++ ) {
    ifiltdef = FILTERCAL_REST->IFILTDEF[ifilt];
    name     = FILTERCAL_REST->FILTER_NAME[ifilt];
    lenf     = strlen(name);
    sprintf(band,"%c", name[lenf-1] );

    printf("\t  %s    ", band);
    for(i=0; i < NBIN_DUMP_RVinv; i++ ) {
      RVinv = RVinv_LIST[i];  RV = 1.0/RVinv ;

      XTMAG_AV1 = genmag_extinction(ifiltdef, Trest, RV, AV1 ) ;
      //XTMAG_AV1 =  eval_XTMAG_AV1(ifilt, Trest, RVinv);
      printf("%6.3f ", XTMAG_AV1);
    }
    printf("\n"); fflush(stdout);
  }

  printf("\t--------------------------------------------- \n");

  fflush(stdout);


  //  if ( Trest > 10.0 ) {  debugexit(fnam); }

  return ;

} // end dump_XTMAG


double genmag_extinction(int ifiltdef, double Trest, double RV, double AV ) {

  // Return XTMAG for inputs Trest, AV, RV.

  int    ifilt = CALIB_INFO.FILTERCAL_REST.IFILTDEF_INV[ifiltdef];
  double XTMAG=0.0, XTMAG_AV1, GRIDVAL_LIST[3];
  char fnam[] = "genmag_extinction" ;

  // ----------- BEGIN ----------

  // check for crazy values

  if ( fabs(AV) > 100.0 ) { return 10.0; }
  if ( RV      < -10.0  ) { return XTMAG; }

  GRIDVAL_LIST[0] = 1.0/RV ;
  GRIDVAL_LIST[1] = Trest ;
  GRIDVAL_LIST[2] = (double)ifilt;

  interp_GRIDMAP( &XTMAG_INFO.GRIDMAP3D, GRIDVAL_LIST, &XTMAG_AV1); 
  XTMAG = AV * XTMAG_AV1 ;

  return XTMAG;

} // end genmag_extinction

double genmag_extinction__(int *ifiltdef, double *Trest, 
			   double *RV, double *AV ) {
  return genmag_extinction(*ifiltdef, *Trest, *RV, *AV); 
}
