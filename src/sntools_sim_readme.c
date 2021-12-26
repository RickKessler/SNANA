/******************************************************
  
   Created Dec 22 2021 [code moved out of snlc_sim]
   Write DOCUMENTATION block to [VERSION].README

******************************************************/

#include "sntools.h"
#include "sntools_cosmology.h"
#include "snlc_sim.h"
#include "sntools_host.h"
#include "sntools_wronghost.h"
#include "sntools_nonlinearity.h"
#include "sntools_genSmear.h"
#include "sntools_trigger.h" 
#include "sntools_genPDF.h"
#include "sntools_sim_readme.h"


// *************************************************
void README_DOCANA_DRIVER(int iflag_readme) {

  // Created Dec 22 2021 by R.Kessler
  //
  // Complete refactor for output README.
  // Replace casual comments with yaml-compliant format
  // under DOCUMENTATION block in several parts:
  // 1. OVERVIEW       # survey, model, host, user ...
  // 2. INPUTS_KEYS    # summarize all user input keys
  // 3. INPUTS_NOTES   # compute rates, etc ...
  // 4. OUTPUT_NOTES   # stats, cpu time ...
  //
  // Input arg:
  //   iflag_readme = 0 => one-time init 
  //   iflag_readme = 1 => write first 3 parts
  //   iflag_readme = 2 => write OUTPUT_NOTES
  //
  // The iflag_readme scheme enables viewing most of the README
  // while job is running; e.g., can check INPUT_KEYS before 
  // long job finishes.

  int i;
  char fnam[] = "README_DOCANA_DRIVER";

  // ----------- BEGIN -----------

  if ( iflag_readme == 0 ) {
    VERSION_INFO.NLINE_README      = 0;
    VERSION_INFO.NLINE_README_INIT = 0;

    README_KEYS_COSMO.NKEY             = 0 ;
    README_KEYS_GENMODEL.NKEY          = 0 ;
    README_KEYS_SIMLIB.NKEY            = 0 ;
    README_KEYS_HOSTLIB.NKEY           = 0 ;
    README_KEYS_RATEMODEL.NKEY         = 0 ;
    README_KEYS_LENS.NKEY              = 0 ;
    README_KEYS_SKY.NKEY               = 0 ;
    README_KEYS_MWEBV.NKEY             = 0 ;
    README_KEYS_NON1ASED.NKEY          = 0 ;
    README_KEYS_SIMSED.NKEY            = 0 ;
    README_KEYS_LCLIB.NKEY             = 0 ;
    README_KEYS_FILTER.NKEY            = 0 ; // keys with _FILTER
    README_KEYS_FLUXERRMODEL.NKEY      = 0 ;

    README_KEYS_GENMAG_OFF.NKEY        = 0 ;
    README_KEYS_GENMAG_SMEAR.NKEY      = 0 ;

    README_KEYS_TAKE_SPECTRUM.NKEY     = 0 ;
    README_KEYS_RANSYSTPAR.NKEY        = 0 ;
    README_KEYS_ZVARIATION.NKEY        = 0 ;
    README_KEYS_GRIDGEN.NKEY           = 0 ;
    README_KEYS_CUTWIN.NKEY            = 0 ;
    README_KEYS_COVMAT_SCATTER.NKEY    = 0 ;
    README_KEYS_SIMGEN_DUMP.NKEY       = 0 ;
    return;
  }


  i = VERSION_INFO.NLINE_README ;

  sprintf(BANNER,"%s: Prepare README content (iflag=%d)", 
	  fnam, iflag_readme);
  print_banner(BANNER);

  if ( iflag_readme == 1 ) { 
    i++; 
    sprintf(VERSION_INFO.README_DOC[i] ,"%s", KEYNAME_DOCANA_REQUIRED ); 

    README_DOCANA_OVERVIEW(&i);
    README_DOCANA_INPUT_KEYS(&i);
  }
  else {
    README_DOCANA_NOTES(&i);

    i++; 
    sprintf(VERSION_INFO.README_DOC[i],"%s",KEYNAME2_DOCANA_REQUIRED); 
  }

  VERSION_INFO.NLINE_README = i;  
  if ( iflag_readme == 1 ) 
    { VERSION_INFO.NLINE_README_INIT = i; }

  return;

} // end README_DOCANA_DRIVER


void README_DOCANA_OVERVIEW(int *iline) {
  int i = *iline;
  char pad[] = "    ", *cptr, cwd[MXPATHLEN] ;

  // ----------- BEGIN ------------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  %s:", DOCANA_OVERVIEW ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSURVEY:       %s",  pad, GENLC.SURVEY_NAME);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sGENMODEL:     %s", pad, INPUTS.GENMODEL);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sHOST_MACHINE: %s", pad, getenv("HOST") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sUSERNAME:     %s",  pad, getenv("USER") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSNDATA_ROOT:  %s", pad, PATH_SNDATA_ROOT );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSNANA_DIR:     %s", pad, PATH_SNANA_DIR );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sSNANA_VERSION: %s", pad, SNANA_VERSION_CURRENT );

  // write current directory (Sep 5 2013)
  if ( getcwd(cwd,MXPATHLEN) != NULL ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sCWD:   %s", pad, cwd );
  }

  int ifile;
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%sINPUT_FILE:", pad );
  for(ifile=0; ifile < INPUTS.NREAD_INPUT_FILE; ifile++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    ENVrestore(INPUTS.INPUT_FILE_LIST[ifile], ORIG_FILE_README);
    sprintf(cptr,"%s- %s", pad, ORIG_FILE_README );
  }

  *iline = i;
  return ;
} // end README_DOCANA_OVERVIEW


void  README_DOCANA_INPUT_KEYS(int *iline) {
  int i = *iline;
  char *cptr, pad[] = "    ", noComment[]="" ;
  int nval1=1, nval2=2, lenkey=20;
  double *dptr, dval, dval_list[10];
  char fnam[] = "README_DOCANA_INPUT_KEYS";

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  INPUT_KEYS:");

  // output genversion, format, nevt ...
  readme_docana_output(&i, pad);

  // source model
  readme_docana_genmodel(&i, pad);

  // instrument: filter, kcor, simlib, noise ...
  readme_docana_instr(&i, pad);

  // hostlib 
  readme_docana_hostlib(&i,pad);

  // search eff maps
  readme_docana_searcheff(&i, pad);

  // redshift and vpec 
  readme_docana_redshift(&i,pad);

  // epoch info (MJD range, Trest range .. )
  readme_docana_epoch(&i,pad);

  // - - - - -  MWEBV - - - - - -
  readme_docana_mwebv(&i,pad);

  // population params for color & stretch
  readme_docana_modelPar(&i,pad);

  // rate model (DNDZ, DNDB ...)
  readme_docana_rate(&i, pad) ;

  // misc
  readme_docana_misc(&i,pad);

  // CUTWIN cut-windows
  readme_docana_cutwin(&i, pad) ;

  *iline = i;
  return;
} // end README_DOCANA_INPUT_KEYS

void README_DOCANA_NOTES(int *iline) {
  int i = *iline;
  // ----------- BEGIN ------------


  *iline = i;
  return;
} // end README_DOCANA_NOTES


void readme_docana_comment(int *iline, char *comment) {
  int i = *iline;
  i++; sprintf(VERSION_INFO.README_DOC[i],"# %s", comment);
  *iline = i;
} // end readme_docana_comment


void readme_docana_genmodel(int *iline, char *pad) {
  
  // store GENMODEL, intrinsic scatter model etc ...
  int i = *iline;
  char *cptr, noComment[]="" ;
  int nval1=1, nval2=2, lenkey=24, o ;  
  double dval, *dptr;

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Source model");

  readme_docana_load_list(&i, pad, &README_KEYS_GENMODEL);

  readme_docana_load_list(&i, pad, &README_KEYS_NON1ASED);

  readme_docana_load_list(&i, pad, &README_KEYS_SIMSED);

  readme_docana_load_list(&i, pad, &README_KEYS_LCLIB);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s ",  pad,lenkey,
	  "GENMAG_SMEAR_MODELNAME:", INPUTS.GENMAG_SMEAR_MODELNAME);

  dval = (double)INPUTS.GENMAG_OFF_GLOBAL ;
  VERSION_INFO_load(&i, pad, "GENMAG_OFF_GLOBAL:", noComment,
		    lenkey, false, nval1, &dval, -20.0,20.0, 0.0); 

  dval = (double)INPUTS.GENMAG_SMEAR[0] ;
  VERSION_INFO_load(&i, pad, "GENMAG_SMEAR:", "coherent scatter (mag)",
		    lenkey, false, nval1, &dval, -20.0,20.0, 0.0); 

  dval = (double)INPUTS.GENMODEL_ERRSCALE ;
  VERSION_INFO_load(&i, pad, "GENMODEL_ERRSCALE:", noComment,
		    lenkey, false, nval1, &dval, 0.0, 20.0, 0.0); 

  readme_docana_load_list(&i, pad, &README_KEYS_GENMAG_OFF);
  readme_docana_load_list(&i, pad, &README_KEYS_GENMAG_SMEAR);

  // - - - - 
  if ( README_KEYS_COSMO.NKEY>0 || README_KEYS_LENS.NKEY > 0 ) {
    readme_docana_comment(&i, "Cosmology inputs");
    readme_docana_load_list(&i, pad, &README_KEYS_COSMO);
    readme_docana_load_list(&i, pad, &README_KEYS_LENS);
  }

  *iline = i;
  return;

} // end readme_docana_genmodel


// ========================================
void readme_docana_instr(int *iline, char *pad) {

  // instrumental: filters, kcor/calib file, simlib, noise
  int i = *iline;
  int nval1=1, nval2=2, lenkey = 20;
  double dval;
  char noComment[] = "";
  char *cptr;
  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Instrumental inputs");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s ",
	  pad, lenkey, "GENFILTERS:", INPUTS.GENFILTERS);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS.KCOR_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s ",
	  pad, lenkey, "KCOR_FILE:", ORIG_FILE_README);

  readme_docana_load_list(&i, pad, &README_KEYS_SIMLIB);

  /* xxxx mark delete 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS.SIMLIB_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s", 
	  pad, lenkey, "SIMLIB_FILE:", ORIG_FILE_README);


  dval = (double)INPUTS.SIMLIB_MSKOPT ;
  VERSION_INFO_load(&i, pad, "SIMLIB_MSKOPT:", noComment,
		    lenkey, true, nval1, &dval, 0.0,1.E12, -1.0); 

  dval = (double)INPUTS.SIMLIB_NREPEAT ;
  VERSION_INFO_load(&i, pad, "SIMLIB_NREPEAT:", noComment,
		    lenkey, true, nval1, &dval, 0.5,100.0, -9.0); 

  if ( strcmp(INPUTS.SIMLIB_FIELDLIST,"ALL") != 0  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s%-*s %s", 
	    pad, lenkey, "SIMLIB_FIELDLIST:", INPUTS.SIMLIB_FIELDLIST);
  }

  dval = (double)INPUTS.SIMLIB_IDSTART ;
  VERSION_INFO_load(&i, pad, "SIMLIB_IDSTART:", noComment,
		    lenkey, true, nval1, &dval, 0.5,1.0E9, -9.0); 

  dval = (double)INPUTS.SIMLIB_IDLOCK ;
  VERSION_INFO_load(&i, pad, "SIMLIB_IDLOCK:", noComment,
		    lenkey, true, nval1, &dval, 0.5,1.0E9, -9.0); 

  dval = (double)INPUTS.SIMLIB_MAXRANSTART ;
  VERSION_INFO_load(&i, pad, "SIMLIB_MAXRANSTART:", noComment,
		    lenkey, true, nval1, &dval, 0.5,1.0E9, -9.0); 

  dval = (double)INPUTS.SIMLIB_MINOBS ;
  VERSION_INFO_load(&i, pad, "SIMLIB_MINOBS:", noComment,
		    lenkey, true, nval1, &dval, 0.5,100.0, -9.0); 

  dval = (double)INPUTS.SIMLIB_MINSEASON ;
  VERSION_INFO_load(&i, pad, "SIMLIB_MINSEASON:", noComment,
		    lenkey, true, nval1, &dval, 0.5,1.0E9, -9.0); 

  dval = (double)INPUTS.SIMLIB_NSKIPMJD ;
  VERSION_INFO_load(&i, pad, "SIMLIB_NSKIPMJD:", noComment,
		    lenkey, true, nval1, &dval, 0.5,1000.0, -9.0); 
  
  dval = (double)INPUTS.SIMLIB_DUMP ;
  VERSION_INFO_load(&i, pad, "SIMLIB_DUMP:", noComment,
		    lenkey, true, nval1, &dval, 0.5,1000.0, -9.0); 

  dval = (double)INPUTS.USE_SIMLIB_REDSHIFT ;
  VERSION_INFO_load(&i, pad, "USE_SIMLIB_REDSHIFT:", noComment,
		    lenkey, true, nval1, &dval, 0.0,9.0, 0.0); 
  dval = (double)INPUTS.USE_SIMLIB_PEAKMJD ;
  VERSION_INFO_load(&i, pad, "USE_SIMLIB_PEAKMJD:", noComment,
		    lenkey, true, nval1, &dval, 0.0,9.0, 0.0); 
  dval = (double)INPUTS.USE_SIMLIB_DISTANCE ;
  VERSION_INFO_load(&i, pad, "USE_SIMLIB_DISTANCE:", noComment,
		    lenkey, true, nval1, &dval, 0.0,9.0, 0.0); 
  dval = (double)INPUTS.USE_SIMLIB_MAGOBS ;
  VERSION_INFO_load(&i, pad, "USE_SIMLIB_MAGOBS:", noComment,
		    lenkey, true, nval1, &dval, 0.0,9.0, 0.0); 
  dval = (double)INPUTS.USE_SIMLIB_SPECTRA ;
  VERSION_INFO_load(&i, pad, "USE_SIMLIB_SPECTRA:", noComment,
		    lenkey, true, nval1, &dval, 0.0,9.0, 0.0); 
  dval = (double)INPUTS.USE_SIMLIB_SALT2 ;
  VERSION_INFO_load(&i, pad, "USE_SIMLIB_SALT2:", noComment,
		    lenkey, true, nval1, &dval, 0.0,9.0, 0.0); 

  xxxxxxxx end mark xxxxxxxx */

  dval = (double)INPUTS.SMEARFLAG_FLUX ;
  VERSION_INFO_load(&i, pad, "SMEARFLAG_FLUX:", "1->add Poisson noise",
		    lenkey, true, nval1, &dval, 0.0,100.0, -9.0); 

  dval = (double)INPUTS.SMEARFLAG_ZEROPT ;
  VERSION_INFO_load(&i, pad, "SMEARFLAG_ZEROPT:", 
		    "+=1->apply scatter, +=2->add to FLUXERRCAL", 
		    lenkey, true, nval1, &dval, 0.0,100.0, -9.0); 

  dval = (double)INPUTS.FUDGE_SNRMAX ;
  VERSION_INFO_load(&i, pad, "FUDGE_SNRMAX:", noComment,
		    lenkey, true, nval1, &dval, 0.0,1.0E5, -9.0); 

  readme_docana_load_list(&i, pad, &README_KEYS_FILTER);

  readme_docana_load_list(&i, pad, &README_KEYS_FLUXERRMODEL);

  if ( README_KEYS_TAKE_SPECTRUM.NKEY > 0 ) {
    readme_docana_comment(&i, "Spectrograph inputs");
    readme_docana_load_list(&i, pad, &README_KEYS_TAKE_SPECTRUM);
  }

  *iline = i;
  return;

} // end readme_docana_instr


void readme_docana_hostlib(int *iline, char *pad) {

  // Created Dec 22 2021
  // Write HOSTLIB_XXX keys and values to DOCUMENTATION block
  //
  // TO DO:
  //  HOSTLIB_FIXSERSIC ??

  int i = *iline;
  int lenkey = 24; // to format printing of key
  double dval, dval_list[20] ;
  char *cptr, *sval, *key;
  int n, nval;
  char  noComment[] = "" ;
  char  fnam[] = "readme_docana_hostlib" ;

  // ----------- BEGIN ------------

  if ( !IGNOREFILE(INPUTS.WRONGHOST_FILE) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    ENVrestore(INPUTS.WRONGHOST_FILE, ORIG_FILE_README);
    sprintf(cptr,"%s%-*s  %s ", 
	    pad, lenkey, "WRONGHOST_FILE:", ORIG_FILE_README);
    *iline = i;
  }

  if ( IGNOREFILE(INPUTS.HOSTLIB_FILE) ) { return; }

  readme_docana_comment(&i, "HOSTLIB inputs");

  readme_docana_load_list(&i, pad, &README_KEYS_HOSTLIB);


  // - - - - -
  *iline = i;
  return;

} // end readme_docana_hostlib


void readme_docana_modelPar(int *iline, char *pad) {

  // load population params for stretch, color, alpha, beta ...

  int i = *iline;
  char *cptr, *s;

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Population and rate-model params");

  s = INPUTS.GENPDF_FILE;
  if ( !IGNOREFILE(s) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%sGENPDF_FILE:  %s \n", pad, s);
  }

  if ( INDEX_GENMODEL == MODEL_SALT2 ) {
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2x1);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2c);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2ALPHA);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_SALT2BETA);
  }
  else {
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_DM15);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_DELTA);
    readme_docana_load_asymGauss(&i, pad, &INPUTS.GENGAUSS_STRETCH);
  }

  readme_docana_load_expHalfGauss(&i, pad, &INPUTS.GENPROFILE_AV);
  readme_docana_load_expHalfGauss(&i, pad, &INPUTS.GENPROFILE_EBV_HOST);
  readme_docana_load_asymGauss   (&i, pad, &INPUTS.GENGAUSS_RV );


  *iline = i;
  return;

} // end readme_docana_modelPar


void readme_docana_rate(int *iline, char *pad) {

  int i = *iline;
  int nval1=1, nval2=2, lenkey=24, m, NMODEL ;
  char *cptr, *s;
  double *dptr, dval;
  // ----------- BEGIN ------------

  readme_docana_load_list(&i, pad, &README_KEYS_RATEMODEL);

  if ( INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT.ORDER > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    s = INPUTS.RATEPAR.DNDZ_ZPOLY_REWGT.STRING;
    sprintf(cptr,"%s%-*s %s # %s",
	    pad, lenkey, "DNDZ_ZPOLY_REWGT:", s, 
	    "dNdz *= polyFun(z)" );
  }

  dptr = &INPUTS.RATEPAR.DNDZ_ZEXP_REWGT ;
  VERSION_INFO_load(&i, pad, "DNDZ_ZEXP_REWGT:", "dN/dz *= z^REWGT",
		    lenkey, false, nval1, dptr, -10.0,10.0, 0.0); 

  dptr = &INPUTS.RATEPAR.DNDZ_ALLSCALE ;
  VERSION_INFO_load(&i, pad, "DNDZ_ALLSCALE:", "dN/dz *= ALLSCALE",
		    lenkey, false, nval1, dptr, 0.0,1.0E4, 1.0); 

  // for DNDZ_SCALE, logic in VERSION_INFO_load doesn't work
  // so need to make explicit check here.
  dptr = INPUTS.RATEPAR.DNDZ_SCALE ;
  if ( dptr[0] != 1.0 || dptr[1] != 1.0 ) {
    VERSION_INFO_load(&i, pad, "DNDZ_SCALE:", 
		      "dN/dz(SNIa,NON1A) *= SCALE[0,1]",
		      lenkey, false, nval2, dptr, 0.0,1.0E4, -9.0 );
  }

  *iline = i;
  return;

} // end readme_docana_rate


void readme_docana_cutwin(int *iline, char *pad) {

  int i = *iline;
  int icut, NCUTWIN = INPUTS.NCUTWIN_TOT;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval;
  char *KEY, *ARG, KEY_LAST[100]="" ;

  // ----------- BEGIN ------------

  if ( NCUTWIN == 0 ) { return; }

  readme_docana_comment(&i,"CUTWIN inputs");

  dval = (double)INPUTS.APPLY_CUTWIN_OPT ;
  VERSION_INFO_load(&i, pad, "APPLY_CUTWIN_OPT:", noComment,
                      lenkey, true, nval1, &dval, 0.0, 100, 0.0 );
  
  readme_docana_load_list(&i, pad, &README_KEYS_CUTWIN);

  *iline = i;
  return;

} // end readme_docana_cutwin


void readme_docana_redshift(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval;

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Redshift inputs");

  dptr = INPUTS.GENRANGE_REDSHIFT;
  VERSION_INFO_load(&i, pad, "GENRANGE_REDSHIFT:", noComment, 
		    lenkey, false, nval2, dptr, 0.0,10.0, -1.0); 

  dval = (double)INPUTS.GENSIGMA_REDSHIFT;
  VERSION_INFO_load(&i, pad, "GENSIGMA_REDSHIFT:", noComment, 
		    lenkey, false, nval1, &dval, 0.0,10.0, -1.0); 

  dval = (double)INPUTS.GENSIGMA_VPEC ;
  VERSION_INFO_load(&i, pad, "GENSIGMA_VPEC:", "true vpec scatter (km/sec)",
		    lenkey, false, nval1, &dval, 0.0,9000.0, -1.0); 

  dval = (double)INPUTS.VPEC_ERR ;
  VERSION_INFO_load(&i, pad, "VPEC_ERR:", 
		    "vpec scatter after correction (km/sec)",
		    lenkey, false, nval1, &dval, 0.0,9000.0, -1.0); 

  dval = (double)INPUTS.VEL_CMBAPEX ;
  VERSION_INFO_load(&i, pad, "VEL_CMBAPEX:", "km/sec",
		    lenkey, false, nval1, &dval, 0.0,400.0, CMBapex_v); 
  

  *iline = i;
  return;

} // end readme_docana_redshift

void readme_docana_epoch(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Epoch & sky coverage inputs");

  dptr = INPUTS.GENRANGE_MJD;
  VERSION_INFO_load(&i, pad, "GENRANGE_MJD:", noComment, 
		    lenkey, false, nval2, dptr, 21000.0,79000.0, -1.0); 

  dptr = INPUTS.GENRANGE_PEAKMJD;
  VERSION_INFO_load(&i, pad, "GENRANGE_PEAKMJD:", noComment, 
		    lenkey, false, nval2, dptr, 1.0E3,1.0E5, -1.0); 

  dval = (double)INPUTS.GENSIGMA_PEAKMJD;
  VERSION_INFO_load(&i, pad, "GENSIGMA_PEAKMJD:", noComment, 
		    lenkey, false, nval1, &dval, 0.0,10.0, 0.0); 

  dval_list[0] = (double)INPUTS.GENRANGE_TREST[0];
  dval_list[1] = (double)INPUTS.GENRANGE_TREST[1];
  VERSION_INFO_load(&i, pad, "GENRANGE_TREST:", noComment, 
		    lenkey, false, nval2, dval_list, -1.0E3,1.0E4, 0.111); 

  dptr = INPUTS.GENRANGE_RA ;
  VERSION_INFO_load(&i, pad, "GENRANGE_RA:", noComment, 
		    lenkey, false, nval2, dptr, -359.0, 360.0, -999.0); 

  dptr = INPUTS.GENRANGE_DEC ;
  VERSION_INFO_load(&i, pad, "GENRANGE_DEC:", noComment, 
		    lenkey, false, nval2, dptr, -359.0, 360.0, -999.0); 

  dval = (double)INPUTS.SOLID_ANGLE ;
  VERSION_INFO_load(&i, pad, "SOLID_ANGLE:", noComment, 
		    lenkey, false, nval1, &dval, 0.0,20.0, 0.0); 

  *iline = i;
  return;

} // end readme_docana_epoch

void readme_docana_misc(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Misc inputs");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s", pad, lenkey, "GENSOURCE:", INPUTS.GENSOURCE);

  dval = (double)INPUTS.ISEED ;
  VERSION_INFO_load(&i, pad, "RANSEED:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  dval = (double)INPUTS.DEBUG_FLAG ;
  VERSION_INFO_load(&i, pad, "DEBUG_FLAG:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  readme_docana_load_list(&i, pad, &README_KEYS_RANSYSTPAR);

  readme_docana_load_list(&i, pad, &README_KEYS_ZVARIATION);

  readme_docana_load_list(&i, pad, &README_KEYS_GRIDGEN);

  readme_docana_load_list(&i, pad, &README_KEYS_SIMGEN_DUMP);

  *iline = i;
  return;

} // end readme_docana_misc


void readme_docana_mwebv(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Galactic extinction");

  readme_docana_load_list(&i, pad, &README_KEYS_MWEBV);

  *iline = i;
  return;

} // end readme_docana_mwebv

void readme_docana_searcheff(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];
  // ----------- BEGIN ------------

  readme_docana_comment(&i, "SEARCHEFF/detections");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s  %s", 
	  pad, lenkey, "SEARCHEFF_PIPELINE_LOGIC_FILE:", ORIG_FILE_README);

  dval = (double)INPUTS.NEWMJD_DIF ;
  VERSION_INFO_load(&i, pad, "NEWMJD_DIF:", 
		    "day-sep if > 1 detections required",
		    lenkey, false, nval1, &dval, 0.0,2000.0, -1.0); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s  %s", 
	  pad, lenkey, "SEARCHEFF_PIPELINE_EFF_FILE:", ORIG_FILE_README);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_SPEC_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s", 
	  pad, lenkey, "SEARCHEFF_SPEC_FILE:", ORIG_FILE_README);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  ENVrestore(INPUTS_SEARCHEFF.USER_zHOST_FILE, ORIG_FILE_README);
  sprintf(cptr,"%s%-*s %s", 
	  pad, lenkey, "SEARCHEFF_zHOST_FILE:", ORIG_FILE_README);

  dval = (double)INPUTS.APPLY_SEARCHEFF_OPT ;
  VERSION_INFO_load(&i, pad, "APPLY_SEARCHEFF_OPT:", 
		    "+= 1,2,4 => pipe,spec,zhost", 
		    lenkey, true, nval1, &dval, 0.0,2000.0, -1.0); 

  *iline = i;
  return;

} // end readme_docana_searcheff


void readme_docana_output(int *iline, char *pad) {

  // Store GENVERSTION, FORMAT_MASK, etc ...

  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------

  readme_docana_comment(&i, "Output data");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s %s", pad, lenkey, "GENVERSION:", INPUTS.GENVERSION);

  dval = (double)INPUTS.NGENTOT_LC ;
  VERSION_INFO_load(&i, pad, "NGENTOT_LC:", noComment,
		    lenkey, true, nval1, &dval, 1.0,1.0E8, 0.0); 

  dval = (double)INPUTS.NGEN_SEASON ;
  VERSION_INFO_load(&i, pad, "NGEN_SEASON:", noComment,
		    lenkey, true, nval1, &dval, 0.0,1.0E8, 0.0); 

  dval = (double)INPUTS.FORMAT_MASK ;
  VERSION_INFO_load(&i, pad, "FORMAT_MASK:", 
		    " += 2,32,16 -> TEXT, FITS, randomCID", 
		    lenkey, true, nval1, &dval, 0.0,2000.0, -1.0); 

  // - - - -- types - - - - -
  dval = (double)GENLC.SIMTYPE ;
  VERSION_INFO_load(&i, pad, "GENTYPE:", "true type", 
		    lenkey, true, nval1, &dval, 0.0,2000.0, -1.0); 

  dval_list[0] = (double)INPUTS.SNTYPE_Ia_SPEC;
  dval_list[1] = (double)INPUTS.SNTYPE_Ia_PHOT;
  VERSION_INFO_load(&i, pad, "SNTYPE:", "spec Type, photID type", 
		    lenkey, true, nval2, dval_list, 0.0,2000.0, -1.0); 

  dval = (double)INPUTS.CIDOFF ;
  VERSION_INFO_load(&i, pad, "CIDOFF:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  dval = (double)INPUTS.CIDRAN_MIN ;
  VERSION_INFO_load(&i, pad, "CIDRAN_MIN:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 

  dval = (double)INPUTS.CIDRAN_MAX ;
  VERSION_INFO_load(&i, pad, "CIDRAN_MAX:", noComment, 
		    lenkey, true, nval1, &dval, 0.0,1.0E9, -1.0); 
 
  *iline = i;
  return;

} // end readme_docana_output

void readme_docana_template(int *iline, char *pad) {
  int i = *iline;
  int nval1=1, nval2=2, lenkey=24 ;
  char *cptr, noComment[]="" ;
  double *dptr, dval, dval_list[10];

  // ----------- BEGIN ------------
  // copy this and change 'template' to something else
  *iline = i;
  return;

} // end readme_docana_template

void  readme_docana_load_list(int *iline, char *pad,
                              README_KEYPLUSARGS_DEF *README_KEYS) {

  // Created Dec 2021
  // Load list of NKEY keys. If KEY is unique, write it as
  //   KEY ARG
  // If there are duplicate keys, write them as
  //   KEY
  //   - ARG0
  //   - ARG1
  //   etc ...
  // Each KEY_LIST item is assumed to include a colon.
  // Beware that duplidates are assumed to be sequential in the list.
  // Non-sequential duplicates will both be written as "KEY ARG".
  // 

  int i = *iline;
  int NKEY = README_KEYS->NKEY;
  int k, lenkey=24 ;
  bool NEW, UNIQUE;
  char *KEY, *ARG, *KEY_NEXT, KEY_LAST[100] = "", *cptr ;

  // ------------- BEGIN -------------

  for(k=0; k < NKEY; k++ ) {
    KEY = README_KEYS->KEY_LIST[k]; 
    ARG = README_KEYS->ARG_LIST[k];
    NEW = (strcmp(KEY,KEY_LAST) != 0) ; // it's a new key

    if ( k < NKEY-1 ) {
      KEY_NEXT = README_KEYS->KEY_LIST[k+1];
      UNIQUE = NEW && (strcmp(KEY,KEY_NEXT) != 0 ); 
    }
    else
      { UNIQUE = NEW; } 

    if ( NEW && UNIQUE ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s%-*s %s", pad, lenkey, KEY, ARG); // write "KEY ARG"
    }
    else if ( NEW ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s%-*s", pad, lenkey, KEY); // write "KEY"
    }

    if ( !UNIQUE ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s- %s", pad, ARG); // write "- ARG"
    }

    sprintf(KEY_LAST,"%s", KEY);
  }

  *iline = i;
  return;

} // end readme_docana_load_list


void readme_docana_load_asymGauss(int *iline, char *pad,
				 GENGAUSS_ASYM_DEF *GENGAUSS) {

  // Utility to load asymGauss params for DOCANA
  int i = *iline ;
  int nval1=1, nval2=2, lenkey=24 ;
  bool USE      = GENGAUSS->USE ;
  char *VARNAME = GENGAUSS->NAME;
  char noComment[] = "" ;

  char keyName[60];
  double *dptr, dval ;
  // ------------ BEGIN -------------
 
  if ( !USE ) { return; }

  // start with required elements

  sprintf(keyName,"GENPEAK_%s:", VARNAME);
  dptr = &GENGAUSS->PEAK ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, -1.E4,1.E4, -9.0); // val,min,max noprint

  sprintf(keyName,"GENSIGMA_%s:", VARNAME);
  dptr = GENGAUSS->SIGMA ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, -1.E4,1.E4, -9.0); // val,min,max noprint

  sprintf(keyName,"GENRANGE_%s:", VARNAME);
  dptr = GENGAUSS->RANGE ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, -1.E4,1.E4, -9.0); // val,min,max noprint

  // optional elements

  sprintf(keyName,"GENGRID_%s:", VARNAME);
  dval = (double)GENGAUSS->NGRID ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, true,
		    nval1, &dval, 0,5, 0.0); // val,min,max noprint

  sprintf(keyName,"GENPEAK2_%s:", VARNAME);
  dptr = &GENGAUSS->PEAK2 ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint

  sprintf(keyName,"GENPROB2_%s:", VARNAME);
  dptr = &GENGAUSS->PROB2 ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint

  sprintf(keyName,"GENSIGMA2_%s:", VARNAME);
  dptr = GENGAUSS->SIGMA2 ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint

  if(  GENGAUSS->PEAKRANGE[1] >  GENGAUSS->PEAKRANGE[0] ) {
    sprintf(keyName,"PEAKRANGE_%s:", VARNAME);
    dptr = GENGAUSS->PEAKRANGE ;
    VERSION_INFO_load(&i, pad, keyName,
		      noComment, lenkey, false,
		      nval2, dptr, -1.E4,1.E4, 0.0); // val,min,max noprint
  }

  *iline = i;
  return;
} // end readme_docana_load_asymGauss


void readme_docana_load_expHalfGauss(int *iline, char *pad,
                                     GEN_EXP_HALFGAUSS_DEF *EXP_HALFGAUASS) {

  // load params for Exp + halfGauss function
  int i = *iline ;
  int nval1=1, nval2=2, lenkey=24 ;
  bool USE      = EXP_HALFGAUASS->USE;
  char *VARNAME = EXP_HALFGAUASS->NAME;
  char noComment[] = "" ;

  char keyName[60];
  double *dptr, dval;

  // ------------- BEGIN -----------
	 
  if ( !USE ) { return; }

  sprintf(keyName,"GENTAU_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->EXP_TAU ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval1, dptr, 0.0, 1.0E5, -9.0); // val,min,max noprint

  sprintf(keyName,"GENRANGE_%s:", VARNAME);
  dptr = EXP_HALFGAUASS->RANGE ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    noComment, lenkey, false, 
		    nval2, dptr, 0.0,100.0, -9.0); // val,min,max noprint
  
  // optonal half-Gauss params
  if ( EXP_HALFGAUASS->RATIO == 0.0 ) { *iline=i; return; }

  sprintf(keyName,"GENGAUPEAK_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->PEAK ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    "peak of half-Gauss", lenkey, false, 
		    nval1, dptr, 0.0,100.0, -9.0); // val,min,max noprint

  sprintf(keyName,"GENSIGMA_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->SIGMA ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    "sigma of half-Gauss", lenkey, false, 
		    nval1, dptr, 0.0,100.0, -9.0); // val,min,max noprint

  sprintf(keyName,"GENRATIO_%s:", VARNAME);
  dptr = &EXP_HALFGAUASS->RATIO ;
  VERSION_INFO_load(&i, pad, keyName, 
  		    "Gauss(0)/Expon(0)", lenkey, false, 
		    nval1, dptr, 0.0,100.0, -9.0); // val,min,max noprint

  *iline = i;
} // end readme_docana_load_expHalfGauss


// *********************************************
void VERSION_INFO_load(int *iline, char *pad, char *keyName,  char *comment,
		       int lenkey, bool isint, int nval,
		       double *val_list, double valmin, double valmax, 
		       double val_noprint ) {

  // created Dec 22 2021
  // Load info to VERSION_INFO.README_DOC if:
  //   valmin <= val[0] <= valmax && val[0] != val_noprint
  // "nval" is the number of values following keyname.
  // include *comment if not null.

  int  ival, i = *iline ;
  double val0   = val_list[0];
  double val ;
  bool noprint  = fabs(val0 - val_noprint) < 1.0E-8;
  bool passcut  = (val0 >= valmin   &&   val0 <= valmax );
  bool do_print = passcut && !noprint ;
  char *cptr, cval[40] ;
  char fnam[] = "VERSION_INFO_load" ;

  // ------------ BEGIN ------------

  if ( !do_print ) { return; }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s%-*s ",  pad, lenkey, keyName );

  for(ival=0; ival < nval; ival++ ) {
    val = val_list[ival];
    if ( isint || val == 0.0 ) 
      { sprintf(cval,"%d  ", (int)val ); }
    else if ( fabs(val) > 0.01 ) 
      { sprintf(cval,"%.4f  ", val ); }
    else
      { sprintf(cval,"%.4le  ", val ); }

    strcat(cptr,cval);
  }
    
  if ( comment[0] != 0 ) { strcat(cptr,"# "); strcat(cptr,comment); }

  *iline = i;
  return ;
} // end VERSION_INFO_load

// =============================================================
void README_KEYPLUSARGS_load(int MXKEY, int NWD, char **WORDS,
			     README_KEYPLUSARGS_DEF *README_KEYS,
			     char *callFun ) {

  // Store NWD WORDS in README_KEYS strut.
  // Inputs:
  //   MXKEY : max number of keys to store
  //   NWD   : number of words past 
  //   WORDS : list of words;  KEY=WORDS[0], ARGS=WORDS[1:N]
  //   callFun: calling function, for error message
  //
  // Output:
  //   README_KEYS : structure to load
  //

  int NKEY = README_KEYS->NKEY;
  int MEMC1 = MXKEY * sizeof(char*);
  int MEMC0 = 100   * sizeof(char);
  int iwd;
  char *KEY, *ARG;
  // ------------ BEGIN ----------

  if ( NKEY == 0 ) {
    // allocate pointer for all MXKEY possible keys
    README_KEYS->KEY_LIST = (char**) malloc(MEMC1);
    README_KEYS->ARG_LIST = (char**) malloc(MEMC1);
  }

  // allodate way more memory for SIMGEN_DUMP
  if ( strstr(WORDS[0],"SIMGEN_DUMP") != NULL ) 
    { MEMC0 = 1000   * sizeof(char);  }

  // allocate 100 chars for this key
  README_KEYS->KEY_LIST[NKEY] = (char*) malloc(MEMC0);
  README_KEYS->ARG_LIST[NKEY] = (char*) malloc(MEMC0);

  KEY = README_KEYS->KEY_LIST[NKEY];  KEY[0]=0;
  ARG = README_KEYS->ARG_LIST[NKEY];  ARG[0]=0;

  sprintf(KEY, "%s", WORDS[0]);

  // load args
  for(iwd=1; iwd<=NWD; iwd++ ) 
    { strcat(ARG,WORDS[iwd]); strcat(ARG," ");  }

  // increment number of stored keys in this structure.
  README_KEYS->NKEY++ ;

  return;

} // README_KEYPLUSARGS_load



// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
// xxxxxxxxxx LEGACY functions below xxxxxxxxxxxx
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

// *********************************************
void readme_doc_legacy(int iflag_readme) {

  // Created 2007 by R.Kessler
  // fill VERSION_INFO.README_DOC structure with 
  // "VERSION_INFO.NLINE_README" lines of README content.
  // Note that the contents are printed from elsewhere
  //
  // iflag_readme = 1 => fill init-part
  // iflag_readme = 2 => fill post-sim part
  //

  char ctmp[MXPATHLEN], cfilt[2], cwd[MXPATHLEN] ;
  char *cptr;
  char conoff[2][4] = { "OFF" , "ON" } ;

  int i, j, j2, ifilt_obs, ifilt_rest, ifilt, itmp, imap, iopt, ipar ;
  int NLINE, NOV, NON1A_non1a, OVP1, OVP2 ;

  double XN, XNERR ;
  float xt, xtprod, val, ZMIN, ZMAX, shift[2];

  char  fnam[] = "readme_doc_legacy" ;

  // ------------ BEGIN readme_doc() ---------------

  i = VERSION_INFO.NLINE_README ;

  print_banner ( " Fill legacy README " );

  if ( iflag_readme == 2 ) goto AFTERSIM ;

  // -----------------------------

  //--- brief description

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s \n", KEYNAME_DOCANA_REQUIRED ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SURVEY:       %s\n",  GENLC.SURVEY_NAME);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  GENMODEL:     %s \n", INPUTS.MODELNAME);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  HOST_MACHINE: %s \n", getenv("HOST") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf( cptr, "  USERNAME:  %s \n", getenv("USER") );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SNDATA_ROOT:  %s \n", PATH_SNDATA_ROOT );
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SNANA_DIR:     %s \n", PATH_SNANA_DIR );
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SNANA_VERSION: %s \n", SNANA_VERSION_CURRENT );

  // write current directory (Sep 5 2013)
  if ( getcwd(cwd,MXPATHLEN) != NULL ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"  CWD:   %s \n", cwd );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s \n", KEYNAME2_DOCANA_REQUIRED ); 

  // indicate changes if "GENPERFECT" is requested.
  readme_doc_GENPERFECT(&i);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n LEGACY_DESCRIPTION: \n" ); 

  readme_doc_SIMLIB(&i);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generation VERSION: %s \n", INPUTS.GENVERSION ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generation source : %s \n", INPUTS.GENSOURCE ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generation FILTERS: %s \n", INPUTS.GENFILTERS ); 

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Rest-frame FILTERS: %s  for  model=%s\n", 
	    GENLC.FILTLIST_REST, INPUTS.MODELNAME );

  }
  // ------------
  // print info for skipped filters

  readme_doc_filterWarn(&i);

  NON1A_non1a = ( INDEX_GENMODEL == MODEL_NON1ASED ) ;
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t KCOR lookup tables: %s \n", INPUTS.KCOR_FILE ); 
    

  if ( GENFRAME_OPT == GENFRAME_REST ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t KCOR lookup  time : Trest" );
    strcat ( cptr, "\n" );

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t AVwarp option for SED : color from " );
    if ( INPUTS.KCORFLAG_COLOR == 1 ) 
      { strcat ( cptr, "closest passbands \n" ); }
    else if ( INPUTS.KCORFLAG_COLOR == 2 ) 
      { strcat ( cptr, "Jha table \n" ); }
    else
      { strcat ( cptr, "UNKNOWN\n" ); }
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.SMEARFLAG_FLUX>0 ) { j=1; } else { j=0; }
  sprintf(cptr, "\t Flux-smearing is %s \n", conoff[j] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Reported flux-uncertainty includes " ); 
  j2 = (INPUTS.SMEARFLAG_FLUX & 2);
  if ( j2 == 0 ) 
    { strcat(cptr,"SKY+GALAXY+SOURCE\n"); }
  else
    { strcat(cptr,"SKY only\n"); } // SMP-like errors


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if( INPUTS.SMEARFLAG_ZEROPT > 0 ) { j=1; } else { j=0; }
  sprintf(cptr, "\t Zeropt-smearing is %s \n", conoff[j] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;  j=0;
  OVP1 = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_PHOT ; 
  if ( OVP1 ) { j = 1; }
  sprintf(cptr, "\t Host-galaxy shot-noise  is %s \n", conoff[j] );


  i++; cptr = VERSION_INFO.README_DOC[i] ;  j=0;
  OVP2 = INPUTS.SMEARFLAG_HOSTGAL & SMEARMASK_HOSTGAL_IMAGE ; 
  if ( OVP2 ) { j = 1; }
  sprintf(cptr, "\t Host-galaxy image-noise  is %s \n", conoff[j] );


  readme_doc_MWXT(&i);


  // ------- RATE INFO ---------

  for ( j=0; j < NLINE_RATE_INFO; j++ ) {
     i++; cptr = VERSION_INFO.README_DOC[i] ;
     sprintf(cptr,"%s\n", LINE_RATE_INFO[j] );
  }

  // ------- generation --------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  GENERATION RANGES: \n");

  if ( INPUTS.USE_SIMLIB_GENOPT ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t ==> Override gen-values with SIMLIB header values.\n");
  }

  if ( INPUTS.USE_SIMLIB_REDSHIFT == 0 ) 
    { sprintf(ctmp,"%s distribution.", INPUTS.RATEPAR.NAME ); }
  else
    { sprintf(ctmp,"SIMLIB values." ); }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,
	  "\t Generate Redshift : %6.3f to %6.3f  using  %s \n"
	  ,INPUTS.GENRANGE_REDSHIFT[0]
	  ,INPUTS.GENRANGE_REDSHIFT[1] 
	  ,ctmp	  );


  if ( INPUTS.GENBIAS_REDSHIFT != 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Generate redshift bias : %f \n", 
	    INPUTS.GENBIAS_REDSHIFT );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.GENSIGMA_REDSHIFT >= 0.0 ) {
    sprintf(cptr,"\t REDSHIFT_FINAL is ZCMB_GEN smeared by : %8.5f \n", 
	    INPUTS.GENSIGMA_REDSHIFT);
  }
  else
    { sprintf(cptr,"\t REDSHIFT_FINAL is host-photoZ \n"); }

  
  if ( WRONGHOST.NLIST > 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t WRONGHOST model frac: %.3f + %.3f*z + %.5f*z^2 \n"
	    ,WRONGHOST.PROB_POLY[0]
	    ,WRONGHOST.PROB_POLY[1]
	    ,WRONGHOST.PROB_POLY[2]     );
  }

  if ( INPUTS.VEL_CMBAPEX == 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,
	    "\t v(cmb)=0 -> REDSHIFT_HELIO=REDSHIFT_CMB (legacy option)\n");
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  bool  USE_HOSTLIB_VPEC = (INPUTS.HOSTLIB_MSKOPT & HOSTLIB_MSKOPT_USEVPEC );
  if ( USE_HOSTLIB_VPEC  ) {
    sprintf(cptr,"\t Pec. Velocity (VPEC) HOSTLIB RMS : %.1f km/sec\n", 
	    HOSTLIB.VPEC_RMS );
  }
  else {
    sprintf(cptr,"\t Pec. Velocity (VPEC) sigma(gen,correct): "
	    "%.1f, %.1f km/sec\n", 
	    INPUTS.GENSIGMA_VPEC, INPUTS.VPEC_ERR );
  }
  if ( INPUTS.RESTORE_WRONG_VPEC ) {
    i++ ; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t\t (RESTORE wrong VPEC sign convention)\n") ;
  }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %s   ZP   offsets : ", INPUTS.GENFILTERS);
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    sprintf(ctmp,"%6.3f ", INPUTS.GENMAG_OFF_ZP[ifilt_obs] ); 
    strcat(cptr,ctmp); 
  }
  strcat(cptr,"\n"); 


  float magoff;
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %s  MODEL offsets : ", INPUTS.GENFILTERS);
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    magoff = INPUTS.GENMAG_OFF_MODEL[ifilt_obs] + INPUTS.GENMAG_OFF_GLOBAL ;
    sprintf(ctmp,"%6.3f ", magoff ); 
    strcat(cptr,ctmp); 
  }
  strcat(cptr,"\n"); 


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %s  exposure times: ", INPUTS.GENFILTERS);
  xtprod = 1.0 ;
  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
    xt = INPUTS.EXPOSURE_TIME_FILTER[ifilt_obs] ;
    xtprod *= xt ;

    if ( xt > 10.0 )   sprintf(ctmp,"%6.1f ", xt ); 
    else               sprintf(ctmp,"%6.4f ", xt ); 

    strcat(cptr,ctmp); 
  } // end if ifilt
  strcat(cptr,"\n"); 

  // indicate which quantities were scaled by EXPOSURE_TIME  
  if ( xtprod != 1.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %s  exposure MSKOPT=%d => ", 
	    INPUTS.GENFILTERS, INPUTS.EXPOSURE_TIME_MSKOPT );
    
    if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 0) )
      { strcat(cptr, "ZPT  "); }
    if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 1) )
      { strcat(cptr, "SKYSIG  "); }
    if ( INPUTS.EXPOSURE_TIME_MSKOPT  & (1 << 2) )
      { strcat(cptr, "READNOISE  "); }

    strcat(cptr,"\n");
  }

  // - - - - - - -  -
  // optional MJD_TEMPLAT (Sep 2017)
  if ( INPUTS.USE_MJD_TEMPLATE ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t MJD_TEMPLATE: " );
    
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
      ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];
      val       = INPUTS.MJD_TEMPLATE_FILTER[ifilt_obs] ;
      if ( val > 1.0 ) {
	sprintf(ctmp,"%6.1f(%c) ", val, FILTERSTRING[ifilt_obs] ); 
	strcat(cptr,ctmp); 
      }
    } // end if ifilt
    strcat(cptr,"\n"); 
  }

  // - - - - - - -
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t RA       : %6.2f to %6.2f  deg\n", 
	  INPUTS.GENRANGE_RA[0], INPUTS.GENRANGE_RA[1] );

  // - - - -  PEAKMJD stuff - - - - - 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t PEAKMJD  : %8.1f to %8.1f   \n", 
	  INPUTS.GENRANGE_PEAKMJD[0], INPUTS.GENRANGE_PEAKMJD[1] );

  if ( INPUTS.GENSIGMA_PEAKMJD > 1.0E-9 ) {
    sprintf(ctmp,"Gauss smear, sigma=%5.2f days", INPUTS.GENSIGMA_PEAKMJD);
  }
  else if ( (INPUTS.OPT_SETPKMJD & OPTMASK_SETPKMJD_FLUXMAX2)>0 ) {
    sprintf(ctmp,"Fmax-clump, MJDWIN=%.1f, SNRCUT>%.1f(3.0)",
	    INPUTS.MJDWIN_SETPKMJD, INPUTS.SNRCUT_SETPKMJD );
  }
  else if ( (INPUTS.OPT_SETPKMJD & OPTMASK_SETPKMJD_FLUXMAX)>0 ) {
    sprintf(ctmp,"naive max flux.");
  }
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t PEAKMJD-estimate  : %s\n", ctmp);

  // - - - - - - - - - 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Trest    : %8.2f to %8.2f  days \n", 
	  INPUTS.GENRANGE_TREST[0], INPUTS.GENRANGE_TREST[1] );

  
  if ( INPUTS.TGRIDSTEP_MODEL_INTERP > 0.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t TGRIDSTEP  : %.1f days (for model-mag interp) \n",
	    INPUTS.TGRIDSTEP_MODEL_INTERP );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t RISETIME-SHIFT(days) SIGMA(lo,hi) : %3.1f , %3.1f  (Mean= %3.1f) \n"
	  ,INPUTS.GENGAUSS_RISETIME_SHIFT.SIGMA[0]
	  ,INPUTS.GENGAUSS_RISETIME_SHIFT.SIGMA[1] 
	  ,INPUTS.GENGAUSS_RISETIME_SHIFT.PEAK
	  );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t FALLTIME-SHIFT(days) SIGMA(lo,hi) : %3.1f , %3.1f  (Mean= %3.1f) \n"
	  ,INPUTS.GENGAUSS_FALLTIME_SHIFT.SIGMA[0]
	  ,INPUTS.GENGAUSS_FALLTIME_SHIFT.SIGMA[1] 
	  ,INPUTS.GENGAUSS_FALLTIME_SHIFT.PEAK
	  );

  // print shape-par info for all models except SIMSED
  if ( INPUTS.NPAR_SIMSED == 0  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(ctmp,"\t Shape-par(%s)", GENLC.SHAPEPAR_NAME );
    sprintf_GENGAUSS(cptr, ctmp, &INPUTS.GENGAUSS_SHAPEPAR);
  }

  // SIMSED parameters
  readme_doc_SIMSED(&i);

  // ---- SALT2 params
  readme_doc_SALT2params(&i);

  // ---- GENPDF populations (Aug 2021)
  readme_doc_GENPDF(&i);

  // ---- FIXMAG params
  readme_doc_FIXMAG(&i);


  // --------------------------------
  // dump host extinction params

  if ( INPUTS.GENPROFILE_AV.USE )
    { readme_doc_hostxt(&i, &INPUTS.GENPROFILE_AV); }
  else if ( INPUTS.GENPROFILE_EBV_HOST.USE ) 
    { readme_doc_hostxt(&i, &INPUTS.GENPROFILE_EBV_HOST); }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  Host Extinction Parameters: NONE  (AV=0) \n");    
  }


  // ====================================
  // Z-dependent parameters (if requested)

  if ( NPAR_ZVAR_USR > 0 ) {
    
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\n  Z-dependent SN-parameter shifts: \n");

    for ( ipar=0; ipar < NPAR_ZVAR_USR; ipar++ ) {
      sprintf(ctmp,"%s", INPUT_ZVARIATION[ipar].PARNAME );
      ZMIN = INPUTS.GENRANGE_REDSHIFT[0] ; 
      ZMAX = INPUTS.GENRANGE_REDSHIFT[1] ; 
      shift[0] = get_zvariation(ZMIN,ctmp) ;
      shift[1] = get_zvariation(ZMAX,ctmp) ;

      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr," %20s-shift = %6.3f(ZMIN=%5.3f), %6.3f(ZMAX=%5.3f) \n",
	      ctmp, shift[0],ZMIN, shift[1], ZMAX );
    }  // end of NPAR_ZVAR_USR loop
  } else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\n  Z-dependent SN-parameter shifts:  None. \n");
  }


  // ===========================================

  readme_doc_magSmear(&i);

  readme_doc_nonLin(&i); // May 27 2016

  // ============================================
  // AVWARP-overflows vs. rest-frame filter

  i++; cptr = VERSION_INFO.README_DOC[i] ;

  sprintf(ctmp, "\n  AVWARP_OVERFLOWS: ");
  if ( NAVWARP_OVERFLOW[0] == 0 ) {
    strcat(ctmp," NONE. ");
  }
  else {
    
    for ( ifilt=0; ifilt < GENLC.NFILTDEF_REST; ifilt++ ) {
      ifilt_rest = GENLC.IFILTMAP_REST[ifilt] ;
      sprintf(cfilt,"%c", FILTERSTRING[ifilt_rest] ); 
      NOV = NAVWARP_OVERFLOW[ifilt_rest] ;
      if ( NOV > 0 ) sprintf(ctmp, "%s %d(%s)", ctmp, NOV, cfilt);
    } 
  }
  strcat(ctmp,"\n");
  strcat(cptr,ctmp);
  sprintf(WARNING_AVWARP_OVERFLOW,"\n  WARNING: %s", ctmp); 

  readme_doc_TAKE_SPECTRUM(&i);

  // ----- cosmology parameters

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Cosmology Parameters: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t H0 = %6.2f km/s per MPc \n", INPUTS.H0 );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Omega_{M,L} = %6.3f, %6.3f   w0,wa = %5.2f,%5.3f  \n",
	  INPUTS.OMEGA_MATTER, INPUTS.OMEGA_LAMBDA, 
	  INPUTS.w0_LAMBDA, INPUTS.wa_LAMBDA );



  // ------ software Search efficiency

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n --------------------------------------------------- \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Software-Pipeline Search Efficiency (MINOBS=%d) from \n", 
	  INPUTS_SEARCHEFF.MINOBS );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[0] ) ;
  for ( imap=0; imap < INPUTS_SEARCHEFF.NMAP_DETECT; imap++ ) {
    NLINE = SEARCHEFF_DETECT[imap].NLINE_README ; 
    for ( itmp = 0; itmp < NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr, "%s", SEARCHEFF_DETECT[imap].README[itmp] ) ;
    }    
  }

  for ( imap=0; imap < INPUTS_SEARCHEFF.NMAP_PHOTPROB; imap++ ) {
    NLINE = SEARCHEFF_PHOTPROB[imap].NLINE_README ; 
    for ( itmp = 0; itmp < NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr, "%s", SEARCHEFF_PHOTPROB[imap].README[itmp] ) ;
    }    
  }

  // print detection logic
  NLINE = SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].NLINE_README ; 
  for ( itmp = 0; itmp < NLINE; itmp++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", SEARCHEFF_DETECT[MXMAP_SEARCHEFF_DETECT].README[itmp] ) ;
  }    
  

  // ------ Spec Search efficiency

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", "\n  Spectroscopic Efficiency : \n" );
  for ( iopt=0; iopt < SEARCHEFF_SPEC_INFO.NLINE_README; iopt++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr, "%s \n", SEARCHEFF_SPEC_INFO.README[iopt] ) ;    
  }
  

  if ( INPUTS_SEARCHEFF.NMAP_zHOST ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", "\n  Unconfirmed zHOST Efficiency map from \n" );
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %s \n", INPUTS_SEARCHEFF.zHOST_FILE );
  }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", "\n  Unconfirmed zHOST Efficiency : 100% \n" );
  }
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  %s \n", COMMENT_README_TRIGGER);

  // print SNTYPE values for SPEC and PHOT Ia-subsets 
  // For NONIA there is only one type specified with the NONIA keys.
  if ( LGEN_SNIA ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"  SNTYPE(Ia) = %d(SPEC)  and %d(PHOT) \n", 
	    INPUTS.SNTYPE_Ia_SPEC, INPUTS.SNTYPE_Ia_PHOT);    
  }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr," --------------------------------------------------- \n");


  // ------ software cuts -----------
  if ( INPUTS.APPLY_CUTWIN_OPT  ) {  readme_doc_CUTWIN(&i) ; }


  // ---- HOST-GALAXY 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  HOSTLIB Summary: " );

  if ( INPUTS.HOSTLIB_USE == 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"None. \n" );
  }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n" );
    NLINE = HOSTLIB.NLINE_COMMENT ;
    for ( itmp=0; itmp < NLINE; itmp++ ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t %s \n", HOSTLIB.COMMENT[itmp] );
    }
  }


  // -----  FUDGES on observing conditions ------
  readme_doc_FUDGES(&i);

  // write list of MAPs so that pipelines can check time-stamps, etc ...
  readme_doc_mapFileList(&i);

  // ======================================
  VERSION_INFO.NLINE_README_INIT = i;   // can dump to here after init
  if ( iflag_readme == 1 ) return ;
  // ======================================

 AFTERSIM:
  i = VERSION_INFO.NLINE_README_INIT ;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n ============ END OF SIMULATION SUMMARY ============== \n");

  // ---- random seed and first/last random

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Random Number Sync: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t RANDOM SEED: %d   (RANLIST_START_GENSMEAR: %d)\n", 
	  INPUTS.ISEED, INPUTS.RANLIST_START_GENSMEAR );


  int ilist;
  sumstat_RANLISTs(2);
  for ( ilist=1; ilist <= GENRAN_INFO.NLIST_RAN; ilist++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"   FIRST/LAST Random Number (List=%d): %f %f  "
	    "AVG(wrap) = %.1f +_ %.1f \n", ilist, 
	    GENRAN_INFO.RANFIRST[ilist], GENRAN_INFO.RANLAST[ilist],
	    GENRAN_INFO.NWRAP_AVG[ilist], GENRAN_INFO.NWRAP_RMS[ilist]	);
  }

  // ---- statistics

  double t_gen   = (TIMERS.t_end - TIMERS.t_end_init); // total time after init
  double R_gen   = (double)NGENLC_TOT / t_gen ;  // NGEN/sec
  double R_write = (double)NGENLC_WRITE/t_gen ;  // NWRITE/sec

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Generation Statistics (gen CPU=%.1f minutes): \n", 
	  t_gen/60.);

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Generated %5d simulated light curves "
	  "(%.f/sec) \n",  NGENLC_TOT, R_gen );
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Wrote     %5d simulated light curves to SNDATA files "
	  "(%.f/sec). \n",  NGENLC_WRITE, R_write );

  if ( NGENSPEC_WRITE > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Wrote     %5d simulated spectra to SNDATA files \n",
	    NGENSPEC_WRITE );
  }

  // spectroscopic tags

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Spectroscopic-type: %d -> %d (before -> after cuts)\n",
	  GENLC.NTYPE_SPEC, GENLC.NTYPE_SPEC_CUTS);
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Photometric-type:   %d -> %d (before -> after cuts)\n",
	  GENLC.NTYPE_PHOT, GENLC.NTYPE_PHOT_CUTS);


  if ( !IGNOREFILE(INPUTS.WRONGHOST_FILE) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    double frac = 0.0 ;
    int N0      = NGENLC_WRITE ;
    int N1      = GENLC.NTYPE_PHOT_WRONGHOST ;
    if ( GENLC.NTYPE_PHOT_CUTS > 0 ) { frac = (double)N1 / (double)N0 ; }
    sprintf(cptr,"  Wrong-host fraction: %d/%d = %.4f\n", N1,N0,frac);
  }

  // ----------------

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  Rejection Statistics: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by NEPOCH<%d \n",
	  NGEN_REJECT.NEPOCH, (int)INPUTS.CUTWIN_NEPOCH[0] ) ;


  double MAGMIN = INPUTS.GENRANGE_PEAKMAG[0];
  double MAGMAX = INPUTS.GENRANGE_PEAKMAG[1];
  if ( MAGMIN > 0.0 && MAGMAX < 9999. ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t %5d rejected by GENRANGE_PEAKMAG(%.1f to %.1f) \n",  
	    NGEN_REJECT.GENMAG, MAGMIN, MAGMAX );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by GENRANGEs \n",  
	  NGEN_REJECT.GENRANGE );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by SEARCH-TRIGGER \n",  
	  NGEN_REJECT.SEARCHEFF );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t %5d rejected by CUTWIN-SELECTION \n",  
	  NGEN_REJECT.CUTWIN );

  // ---

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"  SEARCH+CUTS Efficiency: %7.4f +- %7.4f \n", 
	  GENLC.GENEFF, GENLC.GENEFFERR);

  // give warning if generation stops early
  if ( GENLC.STOPGEN_FLAG  ) {

    bool QUIT_NOREWIND=((INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_QUIT_NOREWIND)>0);
  
    if ( QUIT_NOREWIND ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,
	      "\n  WARNING: GENERATION STOPPED AFTER ONE PASS THRU SIMLIB\n");
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t  AS REQUESTED BY SIM-INPUT SIMLIB_MSKOPT += %d\n",
	      SIMLIB_MSKOPT_QUIT_NOREWIND );
    }
    else {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\n  WARNING: GENERATION STOPPED WHEN ERROR\n");
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t   on ERR(EFF) <= EFFERR_STOPGEN(=%f)\n",  
	      INPUTS.EFFERR_STOPGEN );
      
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t  => YOU HAVE %d FEWER LIGHT CURVES THAN REQUESTED\n",
	      INPUTS.NGEN_LC - NGENLC_WRITE );
    }
  }


  // give SN-stats with cuts
  if ( NLINE_RATE_INFO > 0 ) {
    XN   = (INPUTS.RATEPAR.SEASON_COUNT  + INPUTS.RATEPAR_PEC1A.SEASON_COUNT) ;
    XN  *= GENLC.GENEFF ; // multiply by cut-efficiency
    if ( NGENLC_WRITE > 0 ) 
      { XNERR = XN/sqrt((double)NGENLC_WRITE); }
    else
      { XNERR = 0.0 ; }
    
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\n  Number of SNe per season AFTER CUTS : "
	    "%6.0f +- %5.0f \n", XN, XNERR );
  }


  // NON1ASED
  if ( INPUTS.NON1ASED.NINDEX > 0 ) { readme_doc_NON1ASED(&i); }

  // end marker

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n\t ===== END OF README FILE ====== \n");

  // ------
  if ( i >= MXDOCLINE ) {
    sprintf ( c1err, "%d README.DOC lines exceeds array bound of %d",
	      i, MXDOCLINE );
    sprintf ( c2err,"Increase parameter MXDOCLINE and re-make ");
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err); 
  }

  VERSION_INFO.NLINE_README  = i;

  return;

}  // end of readme_doc_legacy


// ********************************************
void readme_doc_CUTWIN(int *iline) {
  
  int i, icut ;
  char *cptr, ctmp[80] ;
  // ---------- BEGIN ---------

  i = *iline ;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  SOFTWARE CUTS: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t EPOCH CUT: %5.0f Lambda(rest) < %5.0f A \n",
	  INPUTS.EPCUTWIN_LAMREST[0], INPUTS.EPCUTWIN_LAMREST[1] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t EPOCH CUT: SNR >= %2.0f  \n",
	  INPUTS.EPCUTWIN_SNRMIN[0] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t TrestMIN < %5.1f  &&  TrestMAX > %4.1f days \n",
	  INPUTS.CUTWIN_TRESTMIN[1], INPUTS.CUTWIN_TRESTMAX[0] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Max TGAP(rest) <= %5.1f  days \n",
	  INPUTS.CUTWIN_TGAPMAX[1] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Max T0GAP(rest) <= %5.1f  days \n",
	  INPUTS.CUTWIN_T0GAPMAX[1] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t NOBS(MJDDIF > %3.1f) >= %2d \n",
	  INPUTS.CUTWIN_MJDDIF[0], INPUTS.CUTWIN_NOBSDIF[0] );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t NEPOCH(SNR > %4.1f) >= %2.0f \n",
	  INPUTS.CUTWIN_NEPOCH[1], INPUTS.CUTWIN_NEPOCH[0] ) ;

  for ( icut=1; icut <= INPUTS.NCUTWIN_SNRMAX; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SNRMAX > %4.1f for %d of the '%s' filters "
	    "(%5.1f < Trest < %5.1f) \n"
	    , INPUTS.CUTWIN_SNRMAX[icut][0]
	    , INPUTS.CUTWIN_SNRMAX_NFILT[icut]
	    , INPUTS.CUTWIN_SNRMAX_FILTERS[icut]
	    , INPUTS.CUTWIN_SNRMAX_TREST[icut][0]
	    , INPUTS.CUTWIN_SNRMAX_TREST[icut][1]
	    );
  }


  for ( icut=0; icut < INPUTS.NCUTWIN_SATURATE; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t NOBS_SATURATE = %d to %d for each %s filter. \n"
	    , INPUTS.CUTWIN_SATURATE_NOBS[icut][0]
	    , INPUTS.CUTWIN_SATURATE_NOBS[icut][1]
	    , INPUTS.CUTWIN_SATURATE_FILTERS[icut]	    );
  }
  for ( icut=0; icut < INPUTS.NCUTWIN_NOSATURATE; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t NOBS_NOSATURATE = %d to %d for each %s filter. \n"
	    , INPUTS.CUTWIN_NOSATURATE_NOBS[icut][0]
	    , INPUTS.CUTWIN_NOSATURATE_NOBS[icut][1]
	    , INPUTS.CUTWIN_NOSATURATE_FILTERS[icut]	    );
  }


  /* need to add N_DUMP args ....
  if ( INPUTS.NVAR_SIMGEN_DUMP > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(ctmp,"\t SIMGEN_DUMP file includes" );
    if ( INPUTS.APPLY_CUTWIN_OPT  == 1 ) 
      { sprintf(cptr,"%s SNe passing software cuts. \n", ctmp ); }
    else 
      { sprintf(cptr,"%s ALL SNe; use CUTMASK for software cuts.\n", ctmp); }
  }
  xxxx */

  if ( INPUTS.CUTWIN_PEAKMAG[1] < 998.0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t PEAKMAG(All bands) :  %.1f - %.1f\n", 
	    INPUTS.CUTWIN_PEAKMAG_ALL[0], INPUTS.CUTWIN_PEAKMAG_ALL[1] ) ;
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t PEAKMAG(any filter) < %3.1f \n", 
	  INPUTS.CUTWIN_PEAKMAG[1] ) ;
  
  for ( icut=1; icut <= INPUTS.NCUTWIN_PEAKMAG_BYFIELD; icut++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t PEAKMAG(%s): %.1f to %1.f \n", 
	    INPUTS.CUTWIN_BYFIELDLIST[icut],
	    INPUTS.CUTWIN_PEAKMAG_BYFIELD[icut][0],
	    INPUTS.CUTWIN_PEAKMAG_BYFIELD[icut][1] );
  }
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t MWEBV <= %3.1f \n", INPUTS.CUTWIN_MWEBV[1] ) ;
  
  // optional cut on time above SNRMIN
  if ( INPUTS.CUTWIN_EPOCHS_NFILT > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Time < %.1f days for SNR(%s) > %.1f \n",
	    INPUTS.CUTWIN_EPOCHS_TRANGE[1],
	    INPUTS.CUTWIN_EPOCHS_FILTERS,
	    INPUTS.CUTWIN_EPOCHS_SNRMIN );
  }

  // ------------
  *iline  = i ;

  return;

}  // end of readme_doc_CUTWIN


// ********************************************
void  readme_doc_TAKE_SPECTRUM(int *iline) {

  // Mar 14 2017: 
  // write TAKE_SPECTRUM info to readme file.
  //
  // Mar 23 2019: use GENPOLY typedefs and write WARP info

  int N = NPEREVT_TAKE_SPECTRUM ;
  int i, j , OPT_FRAME_EPOCH, OPT_TEXPOSE, IS_HOST ;
  float *ptrEP, *ptrLAM ;
  char *cptr, Tname[20], name2[20], zpolyString[100], warpString[100] ;
  char *ptrFIELD, fieldString[60] ;
  GENPOLY_DEF *GENZPOLY_SNR, *GENZPOLY_TEXPOSE, *GENLAMPOLY_WARP; 

  // ------- BEGIN -------

  if ( N == 0 ) { return ; }

  i = *iline ;
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  TAKE_SPECTRUM: \n");

  
  for(j=0; j < N; j++ ) {

    ptrEP    = INPUTS.TAKE_SPECTRUM[j].EPOCH_RANGE; 
    ptrLAM   = INPUTS.TAKE_SPECTRUM[j].SNR_LAMRANGE ;
    ptrFIELD = INPUTS.TAKE_SPECTRUM[j].FIELD ;

    GENZPOLY_SNR     = &INPUTS.TAKE_SPECTRUM[j].GENZPOLY_SNR ;
    GENZPOLY_TEXPOSE = &INPUTS.TAKE_SPECTRUM[j].GENZPOLY_TEXPOSE ;
    GENLAMPOLY_WARP  = &INPUTS.TAKE_SPECTRUM[j].GENLAMPOLY_WARP ;

    OPT_FRAME_EPOCH = INPUTS.TAKE_SPECTRUM[j].OPT_FRAME_EPOCH ;
    OPT_TEXPOSE = INPUTS.TAKE_SPECTRUM[j].OPT_TEXPOSE ;
    Tname[0] = name2[0] = zpolyString[0] = warpString[0] = IS_HOST = 0 ;
 
    if ( OPT_FRAME_EPOCH == GENFRAME_REST ) 
      { sprintf(Tname,"TREST"); }
    else if ( OPT_FRAME_EPOCH == GENFRAME_OBS )
      { sprintf(Tname,"TOBS"); }
    else if ( OPT_FRAME_EPOCH == GENFRAME_HOST )
      { sprintf(Tname,"HOST");  IS_HOST=1; }

    if ( OPT_TEXPOSE == 1 ) { 
      sprintf(name2, "TEXPOSE" ); 
      sprintf(zpolyString,"%s", GENZPOLY_TEXPOSE->STRING) ;
    }
    else { 
      // sprintf(name2, "SNR(%.0f:%.0f)", ptrLAM[0], ptrLAM[1] ) ; 
      sprintf(name2, "SNR" );
      sprintf(zpolyString,"%s", GENZPOLY_SNR->STRING) ;
    }

    if ( GENLAMPOLY_WARP->ORDER > 0 ) 
      { sprintf(warpString,"WARP:%s", GENLAMPOLY_WARP->STRING ); }

    if ( strlen(ptrFIELD) > 0 ) 
      { sprintf(fieldString,"FIELD=%s", ptrFIELD); }
    else
      { fieldString[0] = 0 ; }

    // - - - - 
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    if ( IS_HOST ) {
      sprintf(cptr,"    %s                      %s-zPOLY:%s  %s  %s\n",
	      Tname, name2, zpolyString, warpString, fieldString );
    }
    else {
      sprintf(cptr,"    %s = %5.1f to %5.1f  "
	      "  %s-zPOLY:%s  %s  %s\n",
	      Tname, ptrEP[0], ptrEP[1],
	      name2, zpolyString, warpString, fieldString );
    }

  } // end j loop over spectra

  // print optional prescale string
  STRING_DICT_DEF *DICT = &INPUTS.DICT_SPECTRUM_FIELDLIST_PRESCALE ;
  if ( DICT->N_ITEM > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"    SPECTRUM_PRESCALE: %s\n", DICT->STRING);
  }

  *iline = i;
  return ;

} // end readme_doc_TAKE_SPECTRUM


// ********************************************
void readme_doc_FUDGES(int *iline) {

  // Created May 2014 (pulled from readme_doc)
  
  int i, NLINE_FUDGE, j ;
  char *cptr ;
  char fudgeLine[10][100] ;

  // ------- BEGIN -------
  i = *iline ;

  NLINE_FUDGE = 0 ;

  // create local comment lines based on used fudges.

  if ( INPUTS.FUDGESCALE_NOISE_SKY != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB NOISE(SKY) : %5.2f ", 
	  INPUTS.FUDGESCALE_NOISE_SKY );
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_NOISE_READ != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB NOISE(CCD-read) : %5.2f ", 
	  INPUTS.FUDGESCALE_NOISE_READ );
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_PSF != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for SIMLIB PSF      : %5.2f ", 
	  INPUTS.FUDGESCALE_PSF);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESHIFT_ZPT != 0.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-shift for SIMLIB ZPT      : %5.2f ", 
	  INPUTS.FUDGESHIFT_ZPT);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGE_MAGERR != 0.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge MAGERR (add in quad)      : %6.4f ", 
	  INPUTS.FUDGE_MAGERR);
    NLINE_FUDGE++ ;
  }
  if ( INPUTS.FUDGE_ZPTERR != 0.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge ZPTERR in SIMLIB          : %6.4f ", 
	  INPUTS.FUDGE_ZPTERR);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.FUDGESCALE_FLUXERR != 1.0 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t Fudge-scale for FLUX-ERROR      : %5.2f ", 
	  INPUTS.FUDGESCALE_FLUXERR);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.OPT_FUDGE_SNRMAX == 1 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t FUDGE_SNRMAX    : %s  (adjust exposure time)",  
	    INPUTS.STRING_FUDGE_SNRMAX);
    NLINE_FUDGE++ ;
  }

  if ( INPUTS.OPT_FUDGE_SNRMAX == 2 ) {
    cptr = fudgeLine[NLINE_FUDGE] ;
    sprintf(cptr,"\t FUDGE_SNRMAX    : %s  (adjust sky noise)",  
	    INPUTS.STRING_FUDGE_SNRMAX);
    NLINE_FUDGE++ ;
  }

  // ----------------------
  // make README comment based on whether any of the fudges are used.

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if( NLINE_FUDGE == 0 ) 
    {  sprintf(cptr,"\n  Fudges on SIMLIB Seeing Conditions: NONE. \n"); }
  else
    {  sprintf(cptr,"\n  Fudges on SIMLIB Seeing Conditions: \n"); }


  // print summary of USED fudges only.
  for(j=0; j < NLINE_FUDGE; j++ ) {
    i++ ; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "%s\n", fudgeLine[j]);
  }


  *iline = i;

  return ;

} // end readme_doc_FUDGES

// ******************************************
void readme_doc_mapFileList(int *iline) {

  // write each map file so that higher level pipelines
  // can grep out these files and compare time stamps
  // against the sim data time stamp.

  int i;
  // -------- BEGIN ---------

  i = *iline;

  i++; 
  sprintf(VERSION_INFO.README_DOC[i], "\n");

  readme_doc_mapFile(&i, "SIMLIB_FILE:", INPUTS.SIMLIB_FILE);
  readme_doc_mapFile(&i, "KCOR_FILE:",   INPUTS.KCOR_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_WGTMAP_FILE:", 
		     INPUTS.HOSTLIB_WGTMAP_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_ZPHOTEFF_FILE:", 
		     INPUTS.HOSTLIB_ZPHOTEFF_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_SPECBASIS_FILE:", 
		     INPUTS.HOSTLIB_SPECBASIS_FILE);
  readme_doc_mapFile(&i, "HOSTLIB_SPECDATA_FILE:", 
		     INPUTS.WRONGHOST_FILE);
  readme_doc_mapFile(&i, "WRONGHOST_FILE:", 
		     INPUTS.WRONGHOST_FILE);
  readme_doc_mapFile(&i, "FLUXERRMODEL_FILE:",
		     INPUTS.FLUXERRMODEL_FILE);
  readme_doc_mapFile(&i, "NONLINEARITY_FILE:",
		     INPUTS.NONLINEARITY_FILE);
  readme_doc_mapFile(&i, "ZVARIATION_FILE:",
		     INPUT_ZVARIATION_FILE );
  readme_doc_mapFile(&i, "WEAKLENS_PROBMAP_FILE" ,
		     INPUTS.WEAKLENS_PROBMAP_FILE );
  readme_doc_mapFile(&i, "SEARCHEFF_PIPELINE_LOGIC_FILE:" ,
		     INPUTS_SEARCHEFF.USER_PIPELINE_LOGIC_FILE);
  readme_doc_mapFile(&i, "SEARCHEFF_PIPELINE_EFF_FILE:" ,
		     INPUTS_SEARCHEFF.USER_PIPELINE_EFF_FILE);
  readme_doc_mapFile(&i, "SEARCHEFF_SPEC_FILE:" ,
		     INPUTS_SEARCHEFF.USER_SPEC_FILE);
  readme_doc_mapFile(&i, "SEARCHEFF_zHOST_FILE:" ,
		     INPUTS_SEARCHEFF.USER_zHOST_FILE);

  *iline = i;

  return ;

}  // end readme_doc_mapFile_list

void readme_doc_mapFile(int *iline, char *KEY, char *FILENAME) {

  int i;
  char *cptr  ;
  char KEY_MAP[] = "MAP:" ;

  // -------------- BEGIN ------------
  i = *iline ;
  
  if ( !IGNOREFILE(FILENAME) ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s %s %s\n", KEY_MAP, KEY, FILENAME);
  }


  *iline = i;
  return ;

} // end readme_doc_mapFile

// ********************************************
void readme_doc_GENPERFECT(int *iline) {

  // April 2014

  int i, itmp ;
  char *cptr, *PARNAME ;
  double DPARVAL_ORIG, DPARVAL_USED ;
  int    IPARVAL_ORIG, IPARVAL_USED ;
  // --------- BEGIN ---------

  i = *iline ;

  if ( GENPERFECT.NVAR <= 0 ) { return ; }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n PERFECT LIGHTCURVES REQUESTED: \n");
  
  for ( itmp=1; itmp <= GENPERFECT.NVAR; itmp++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    
    PARNAME      = GENPERFECT.parnam[itmp] ;
    DPARVAL_ORIG = GENPERFECT.parval[itmp][0] ;
    DPARVAL_USED = GENPERFECT.parval[itmp][1] ;

    if ( GENPERFECT.partype[itmp] == 1 ) {
      IPARVAL_ORIG = (int)DPARVAL_ORIG ;
      IPARVAL_USED = (int)DPARVAL_USED ;
      sprintf(cptr,"\t %-28.28s : %d -> %d \n",
	      PARNAME, IPARVAL_ORIG, IPARVAL_USED );
    }
    else {
      sprintf(cptr,"\t %-28.28s : %.5f -> %.5f \n" ,
	      PARNAME, DPARVAL_ORIG, DPARVAL_USED );
    }
  }

  *iline = i;

}  // end of readme_doc_GENPERFECT

// ********************************************
void readme_doc_NON1ASED(int *iline) {

  // Created Feb 1 2017    [code moved out of readme_doc()]
  //
  // Print SEDs which have NGENTOT>0.

  int i, j, isp, index, NINDEX ;
  int NGENTOT, NGENWR, NGENUSR, CID0, CID1 ;
  float eff;
  char *cptr,  *ptrtype, ctmp[100], cline[200] ;
  char fnam[] = "readme_doc_NON1ASED" ;

  // --------------- BEGIN --------------

  i = *iline ;

  NINDEX = INPUTS.NON1ASED.NINDEX ;

  sprintf(cline,
	  " ---------------------------------------------------"
	  "-------------------------\n");

  // summarize user inputs
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n\n  NON1A USER INPUTS vs. INDEX: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr," INDEX(TYPE) " );
  for ( j = 2; j <= INPUTS.NON1ASED.NKEY; j++ )
    { sprintf(cptr,"%s %8s", cptr, INPUTS.NON1ASED.KEYLIST[j] ); }
  strcat(cptr,"\n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"%s", cline);

  for ( isp=1; isp <= NINDEX; isp++ ) {
    index   = INPUTS.NON1ASED.INDEX[isp];
    ptrtype = GENLC.NON1ASED.TYPE[isp];
    NGENTOT = GENLC.NON1ASED.NGENTOT[isp];
    if ( NGENTOT == 0 ) { continue ; }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr," %3.3d (%5s) ", index, ptrtype );
    
    for ( j = 2; j <= INPUTS.NON1ASED.NKEY; j++ )
      { sprintf(cptr,"%s %8.3f", cptr, INPUTS.NON1ASED.KEYVAL[isp][j] ); }
    
    sprintf(cptr,"%s  %s\n", cptr,  INPUTS.NON1ASED.LIST_NAME[index] );
  } // end isp loop
  
  // ----------------------------------------
  // now summarize stats, effic, CID-range, peakmags ...
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  NON1A-SUMMARY vs. INDEX: \n");
  
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"                  NGEN    NGEN   SEARCH   \n");
  
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr," INDEX(TYPE)      written total  Effic    CID-range \n" );

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s", cline);

    for ( isp=1; isp <= NINDEX; isp++ ) {
      index     = INPUTS.NON1ASED.INDEX[isp];
      NGENTOT   = GENLC.NON1ASED.NGENTOT[isp];
      NGENWR    = GENLC.NON1ASED.NGENWR[isp];
      NGENUSR   = INPUTS.NON1ASED.NGEN[isp];
      ptrtype   = GENLC.NON1ASED.TYPE[isp];
      CID0      = GENLC.NON1ASED.CIDRANGE[isp][0];
      CID1      = GENLC.NON1ASED.CIDRANGE[isp][1];

      if ( NGENTOT > 0 ) { eff = (float)NGENWR/(float)NGENTOT ; }
      else               { eff = 0.0 ; continue ; }
      
      // glue together rest-frame peakmags

      ctmp[0] = 0 ;
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr," %3.3d (%10s)  %4d  %5d   %5.3f %6d -%6d\n", 
	      index, ptrtype, NGENWR, NGENTOT, eff, CID0, CID1 );

      if ( NGENUSR != NGENWR  && INPUTS.NGEN_LC > 0 ) {
	sprintf(c1err,"NGEN[NONIA index=%d]=%d, but expected NGEN=%d",
		index, NGENWR, NGENUSR );
	errmsg(SEV_WARN, 0, fnam, c1err, ""); 
      }

    }  // end of isp loop

  *iline = i;

  return ;

} // end readme_doc_NON1ASED


// ********************************************
void readme_doc_SIMLIB(int *iline) {

  // add SIMLIB info to readme file.
  int i, j, NSKIP ;
  char *cptr, ctmp[100];

  // ------------- BEGIN ---------------

  i = *iline ;

  // - - - - - - 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( GENLC.SIMLIB_IDLOCK >= 0 ) 
    { sprintf(ctmp,"Lock LIBID=%d", GENLC.SIMLIB_IDLOCK ); }
  else
    { sprintf(ctmp,"start at LIBID=%d", INPUTS.SIMLIB_IDSTART ); }

  if ( INPUTS.SIMLIB_NREPEAT > 1 ) // Apr 26 2017
    { sprintf(ctmp,"NREPEAT=%d", INPUTS.SIMLIB_NREPEAT ); } 

  sprintf(cptr,"\t SIMLIB filename  : %s (%s) \n", 
	  INPUTS.SIMLIB_FILE, ctmp ); 

  // - - - - - - 
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t SIMLIB SURVEY    : %s  (TELESCOPE=%s, MINOBS=%d) \n", 
	  GENLC.SURVEY_NAME, GENLC.TELESCOPE[0], INPUTS.SIMLIB_MINOBS  ); 

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t SIMLIB UNITS     : %s for PSF,  %s for SKYSIG \n", 
	  SIMLIB_GLOBAL_HEADER.PSF_UNIT, SIMLIB_GLOBAL_HEADER.SKYSIG_UNIT  ); 

  if ( INPUTS.SIMLIB_MSKOPT != 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB MSKOPT   : %d \n", INPUTS.SIMLIB_MSKOPT );
  }

  if ( strcmp(INPUTS.SIMLIB_FIELDLIST,"ALL") != 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\t SIMLIB FIELDLIST : %s\n", INPUTS.SIMLIB_FIELDLIST); 
  }
  


  double MINSEASON = INPUTS.SIMLIB_MINSEASON;
  if ( MINSEASON > 0.01 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB min season len : %.0f days", MINSEASON);
  }
 

  int NPE_SAT  = SIMLIB_GLOBAL_HEADER.NPE_PIXEL_SATURATE ;
  int PHOTFLAG = SIMLIB_GLOBAL_HEADER.PHOTFLAG_SATURATE ;
  if ( NPE_SAT < 999999999 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB Saturation : %d pe in central pixel "
	    "-> PHOTFLAG=%d\n",  NPE_SAT, PHOTFLAG);
  }

  if ( INPUTS.FUDGEOPT_FLUXERR ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t FUDGEOPT_FLUXERR : %d \n",
	    INPUTS.FUDGEOPT_FLUXERR );
  }

  if ( INPUTS.SIMLIB_IDLOCK > 1 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB ID LOCKED to LIBID : %d \n", 
	    GENLC.SIMLIB_IDLOCK ); 
  }
  if ( INPUTS.SIMLIB_IDLOCK == 1 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB ID LOCKED to first accepted LIBID : %d \n", 
	    GENLC.SIMLIB_IDLOCK ); 
  }

  NSKIP = INPUTS.NSKIP_SIMLIB ;
  if ( NSKIP > 0 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t SIMLIB IDs skipped : " );
    for ( j=0; j < NSKIP ; j++ ) {
      sprintf(ctmp," %d", INPUTS.SIMLIB_IDSKIP[j] );
      strcat(cptr, ctmp);
    }
    strcat ( cptr, "\n" );
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t NEWMJD_DIF       : "
	  "%5.2f  minutes (defines trigger epoch)\n", 
	  24.*60.*INPUTS.NEWMJD_DIF ); 

  *iline = i;


  //  int    KEEP_ENTIRE_SEASON = 
  //    (INPUTS.SIMLIB_MSKOPT & SIMLIB_MSKOPT_ENTIRE_SEASON );

} // end of readme_doc_SIMLIB


// *******************************
void  readme_doc_magSmear(int *iline) {

  int i, j, ifilt, onoff;
  char *cptr, ctmp[80] ;
  char conoff[2][4] = { "OFF" , "ON" } ;
  //  char fnam[] = "readme_doc_magSmear" ;

  // ----------------- BEGIN -----------------
  i = *iline ;

  // intrinsic smearing params

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Intrinsic MAG-smearing models "
	  "(sigma clip %4.1f to %4.1f) : \n"
	  ,INPUTS.SIGMACLIP_MAGSMEAR[0]
	  ,INPUTS.SIGMACLIP_MAGSMEAR[1]
	  );

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 1: Coherent MAG-smearing (GENMAG_SMEAR) : %6.3f  \n", 
	  INPUTS.GENMAG_SMEAR[0] );


  onoff = 0;
  if ( INPUTS.GENMODEL_ERRSCALE > 0.0 ) { onoff=1; }
  if ( INPUTS.NFILT_SMEAR       > 0   ) { onoff=1; }


  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 2: passband MAG-smearing is %s \n", conoff[onoff] );

  if ( onoff > 0 ) {
    if ( INPUTS.GENMODEL_ERRSCALE_OPT == 1 )  
      { sprintf(ctmp,"PEAK"); }
    else           
      { sprintf(ctmp,"EPOCH-DEPENDENT"); }


    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Scale %s model errors (GENMODEL_ERRSCALE) : %6.3f  \n", 
	    ctmp, INPUTS.GENMODEL_ERRSCALE );

    if ( INPUTS.NFILT_SMEAR > 0 ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t Fixed ");
      for ( j=1; j <= INPUTS.NFILT_SMEAR; j++ ) {
	ifilt = INPUTS.IFILT_SMEAR[j]; 
	sprintf(cptr, "%s%c", cptr, FILTERSTRING[ifilt] ) ;
      }
      strcat(cptr," mag-smear:");

      for ( j=1; j <= INPUTS.NFILT_SMEAR; j++ ) {
	ifilt = INPUTS.IFILT_SMEAR[j]; 
	sprintf(cptr, "%s %4.2f", cptr, INPUTS.GENMAG_SMEAR_FILTER[ifilt] ) ;
      }
      strcat(cptr,"\n");

    }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t Correlation between models 1 & 2 : %5.2f \n",
	    INPUTS.GENMODEL_ERRSCALE_CORRELATION );

  } // NFILT_SMEAR

  // --------------------------------
  // restlambda-dependent smearing using  GENMAG_SMEAR_MODELNAME

  onoff = 0;   ctmp[0] = 0 ;
  if ( istat_genSmear() > 0 ) 
    { onoff=1;     sprintf(ctmp,"%s", INPUTS.GENMAG_SMEAR_MODELNAME); }
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 3: %s model-smear is %s  \n", 
	  ctmp, conoff[onoff] );

  // print override parameter names xxxgenmag
  for(j=0; j < NSMEARPAR_OVERRIDE; j++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\t    %s-Override :  %d values for %s\n" 
	    , ctmp
	    , GENMAG_SMEARPAR_OVERRIDE[j].NVAL
	    , GENMAG_SMEARPAR_OVERRIDE[j].NAME );
  }

  // -------------------------
  // model 4: intrinsic scatter matrix

  if ( INPUTS.NCOVMAT_SCATTER > 0 ) { onoff=1;} else { onoff=0;}

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 4: intrinsic scatter matrix is %s \n", 
	  conoff[onoff] );

  /* xxx mark delete 
  char *ptrScat; 
  for ( j=0; j < 8; j++ ) {
    ptrScat = README_KEYPLUSARGS.COVMAT_SCATTER[j] ;
    if ( strlen(ptrScat) > 0 ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"%s\n", ptrScat ); 
    }
  }
  xxxxxxxx */
  
  // model 5: SMEAR_USRFUN
  if ( INPUTS.NPAR_GENSMEAR_USRFUN > 0 )  { onoff=1;} else { onoff=0;}
  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"   Model 5: GENMAG_SMEAR_USRFUN is  %s \n", 
	  conoff[onoff] );
  
  if ( onoff == 1 ) {
    i++ ; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"\t FUNPAR[1-%d] = ", INPUTS.NPAR_GENSMEAR_USRFUN );
    for ( j=0; j < INPUTS.NPAR_GENSMEAR_USRFUN; j++ )  { 
      sprintf(ctmp,"%6.2f ", INPUTS.GENMAG_SMEAR_USRFUN[j] );  
      cptr = strcat(cptr,ctmp);
    }
    cptr = strcat(cptr,"\n");
  }



  *iline = i;
  return ;

} // end of readme_doc_magSmear

// ***************************************
void  readme_doc_nonLin(int *iline) {

  int i,j;
  char *cptr ;
  // ------------- BEGIN ----------------

  if ( NONLIN_README.NLINE == 0 ) { return ; }

  i = *iline ;

  for(j=0; j < NONLIN_README.NLINE; j++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr,"%s\n", NONLIN_README.LINE[j] );
  }

  *iline = i;

  return ;

} // end readme_doc_nonLin

// ***************************************
void  readme_doc_SIMSED(int *iline) {

  int i, ipar, iflag ;
  char *cptr, ctmp[80], ctmp2[80] ;

  i = *iline ;

  // xxx mark? if ( INPUTS.OPTMASK_SIMSED == OPTMASK_GEN_SIMSED_GRIDONLY ) {
  if ( INPUTS.OPTMASK_SIMSED == 4 ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;  
    sprintf(cptr,"\n   %s model option is GRIDONLY & SEQUENTIAL \n", 
	    INPUTS.MODELNAME );
  } 
  else if ( INPUTS.NPAR_SIMSED > 0  ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;  
    sprintf(cptr,"\n   %s model parameters: \n", 
	    INPUTS.MODELNAME );
  }


  for ( ipar = 0; ipar < INPUTS.NPAR_SIMSED; ipar++ ) {

    iflag = INPUTS.GENFLAG_SIMSED[ipar];

    if ( (iflag & 1) == 0 ) 
      { continue ; } // baggage

    if ( (iflag & 4) > 0 ) 
      { sprintf(ctmp, "GRIDONLY" ); }
    else
      { sprintf(ctmp, "contin." ); }

    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(ctmp2,"    %s(%s)", INPUTS.PARNAME_SIMSED[ipar], ctmp );
    sprintf_GENGAUSS(cptr, ctmp2, &INPUTS.GENGAUSS_SIMSED[ipar] );

  }


  *iline = i;

} // end of  readme_doc_SIMSED

// *******************************
void  readme_doc_MWXT(int *iline) {

  int i ;
  char *cptr ;

  i = *iline ;
   
  if ( INPUTS.MWEBV_FLAG ) {

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t MilkyWay extinction  is ON  \n" );

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t    Color law: %s  (OPT_MWCOLORLAW=%d) \n",  
	    INPUTS.STR_MWCOLORLAW, INPUTS.OPT_MWCOLORLAW );


    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    if ( INPUTS.GENRANGE_MWEBV[0] >= 0.0 ) {
      sprintf(cptr, "\t    E(B-V) randomly picked between %.3f and %.3f\n",  
	      INPUTS.GENRANGE_MWEBV[0], INPUTS.GENRANGE_MWEBV[1] );
    }
    else {
      sprintf(cptr, "\t    E(B-V): %s   (OPT_MWEBV=%d)\n",  
	      INPUTS.STR_MWEBV, INPUTS.OPT_MWEBV );
    }

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t    sigma(MWEBV) = %.2f*MWEBV + %.2f  \n", 
	    INPUTS.MWEBV_SIGRATIO, INPUTS.MWEBV_SIG );

    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t    shift(MWEBV) = %5.3f mag \n", 
	    INPUTS.MWEBV_SHIFT );

    if ( INPUTS.APPLYFLAG_MWEBV ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ; 
      sprintf(cptr, "\t    Correct dataFile FLUXCAL for MWEBV \n");
    }
    
  }
  else {
    i++; cptr = VERSION_INFO.README_DOC[i] ; 
    sprintf(cptr, "\t MilkyWay extinction  is OFF \n" );
  }

  *iline = i ;

} // end of readme_doc_MWXT

// *******************************
void readme_doc_filterWarn(int *iline) {

  double ZMIN, ZMAX;
  int i, ifilt, ifilt_obs, NSKIP;
  char *cptr, cfilt[2];

  i = *iline ;

  for ( ifilt=0; ifilt < GENLC.NFILTDEF_OBS; ifilt++ ) {
    ifilt_obs = GENLC.IFILTMAP_OBS[ifilt];

    sprintf(cfilt,"%c", FILTERSTRING[ifilt_obs] ); 
    NSKIP = NSKIP_FILTER[ifilt_obs] ;
    ZMIN  = ZVALID_FILTER[0][ifilt_obs] ;
    ZMAX  = ZVALID_FILTER[1][ifilt_obs] ;

    if ( NSKIP > 0 && ZMIN > 100. ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t FILTER-WARNING: %s not generated for any z-range. \n", 
	      cfilt);
    }

    if ( NSKIP > 0 && ZMIN < 10. ) {
      i++; cptr = VERSION_INFO.README_DOC[i] ;
      sprintf(cptr,"\t FILTER-WARNING: %s generated only for %5.3f < z < %5.3f\n",
	      cfilt, ZMIN, ZMAX );       
    }

  }

  *iline = i;

} // end of readme_doc_filterWarn


// *******************************
void readme_doc_hostxt(int *iline, GEN_EXP_HALFGAUSS_DEF *GENPROFILE) {

  // add host extinction info to README 
  // Apr 4 2020: refactor and update to allow EBV or AV spec.

  char *VARNAME = GENPROFILE->NAME;
  int i;
  char *cptr;
  char cEXPON[40], cGAUSS[40] ;

  // ------------ BEGIN --------------

  i = *iline ;

  // write host extinct info

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\n  Host Extinction Parameters: \n");

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf_GENGAUSS(cptr, "\t RV ", &INPUTS.GENGAUSS_RV);
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr,"\t Gen-Range for %s  : %4.2f to %4.2f  (model=%s) \n", 
	  VARNAME,GENPROFILE->RANGE[0], GENPROFILE->RANGE[1], INPUTS.GENSNXT );
  
 
  
  double TAU   = GENPROFILE->EXP_TAU;
  double RATIO = GENPROFILE->RATIO;
  double SIG   = GENPROFILE->SIGMA;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( TAU < 50.0 ) {
    
    sprintf(cEXPON,"ExpNONE");
    sprintf(cGAUSS,"GaussNONE");

    if ( TAU > 1.0E-9 ) 
      { sprintf(cEXPON,"exp(-%s/%4.2f)", VARNAME, TAU );   }
    if ( INPUTS.GENGAUSIG_AV > 1.0E-9 ) 
      { sprintf(cGAUSS,"%4.2f x Gauss(%s,sig=%4.2f)", RATIO, VARNAME,SIG );  }
    
    sprintf(cptr,"\t dN/d%s = %s + %s \n", VARNAME, cEXPON, cGAUSS ); 
  }
  else  { 
    sprintf(cptr,"\t dN/d%s = flat \n", VARNAME ); 
  }
  

  *iline = i ;

} // end of readme_doc_hostxt


// **************************************
void readme_doc_FIXMAG(int *iline ) {

  int i ;
  char *cptr ;

  if ( INDEX_GENMODEL != MODEL_FIXMAG ) { return ; }

  i = *iline ;

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(cptr, "\t %s Range: %.2f to %.2f",
	  INPUTS.GENMODEL, INPUTS.FIXMAG[0], INPUTS.FIXMAG[1] );

  *iline = i;

} // end readme_doc_FIXMAG


// **************************************
void readme_doc_SALT2params(int *iline ) {

  // Aug 11 2021: write SALT2 only if asymGauss fun is used

  int i ;
  char *cptr, string[40] ;
  char star[2];

  // ------------ BEGIN ---------

  if ( INDEX_GENMODEL != MODEL_SALT2 ) { return ; }

  i = *iline ;

  if ( INPUTS.GENGAUSS_SALT2c.USE ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf_GENGAUSS(cptr, "\t SALT2c", &INPUTS.GENGAUSS_SALT2c);
  }

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  sprintf(string, "\t Alpha" );
  sprintf_GENGAUSS(cptr, string, &INPUTS.GENGAUSS_SALT2ALPHA );
  

  i++; cptr = VERSION_INFO.README_DOC[i] ;
  if ( INPUTS.SALT2BETA_cPOLY.ORDER >= 0 ) { 
    sprintf(cptr,"\t Beta(c) = %s \n", INPUTS.SALT2BETA_cPOLY.STRING);
  }
  else {
    sprintf(string, "\t Beta ");
    sprintf_GENGAUSS(cptr, string, &INPUTS.GENGAUSS_SALT2BETA);
  }


  *iline = i ;

  return; 

} // end of readme_doc_SALT2params
 

// **************************************
void readme_doc_GENPDF(int *iline ) {

  // Aug 11 2021: write SALT2 only if asymGauss fun is used

  int i, imap ;
  char *cptr, string[40] ;
  char star[2];

  // ------------ BEGIN ---------

  if ( INDEX_GENMODEL != MODEL_SALT2 ) { return ; }

  i = *iline ;

  for (imap=0; imap < NMAP_GENPDF; imap++ ) {
    i++; cptr = VERSION_INFO.README_DOC[i] ;
    sprintf(cptr, "\t GENPDF: %s(%s)\n", 
	    GENPDF[imap].MAPNAME,  GENPDF[imap].GRIDMAP.VARLIST );
  }

  *iline = i ;

  return; 

} // end of readme_doc_GENPDF 
