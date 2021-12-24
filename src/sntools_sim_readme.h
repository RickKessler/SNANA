/***********************************
   Created Dec 22 2021

**********************************/

#define DOCANA_OVERVIEW      "OVERVIEW"
#define DOCANA_INPUT_KEYS    "INPUT_KEYS"
#define DOCANA_INPUT_NOTES   "INPUT_NOTES"
#define DOCANA_OUTPUT_NOTES  "OUTPUT_NOTES"

char  ORIG_FILE_README[MXPATHLEN]; // temp space to restore original filenames

// define KEY+ARG strings for keys that have variable number of 
// arguments, or that can have duplicate keys.
struct {

  char GENMODEL[MXPATHLEN];                // GENMODEL 

  char *RATEMODEL[MXRATEPAR_ZRANGE+1];          // for DNDZ and DNDB
  char *RATEMODEL_PEC1A[MXRATEPAR_ZRANGE+1];    // for DNDZ_PEC1A

  // store each CUTWIN; store key and arg separately to track duplicates
  char *KEY_CUTWIN[MXCUTWIN];  
  char *ARG_CUTWIN[MXCUTWIN];

  // store filter-dependent keys that can appear multiple times.
  int NKEY_FILTER;     
  char *KEY_FILTER[MXFILTINDX]; // store each key to track duplicdates
  char *ARG_FILTER[MXFILTINDX]; // store each arg separately

  char COVMAT_SCATTER[10][100]; // very old ... need to check??

} README_KEYPLUSARGS ;

// -------- function prototypes -----------

void  README_DOCANA_DRIVER(int iflag_readme);
void  README_DOCANA_OVERVIEW(int *iline);
void  README_DOCANA_INPUT_KEYS(int *iline);
void  README_DOCANA_NOTES(int *iline);

void  readme_docana_genmodel(int *iline, char *pad);
void  readme_docana_redshift(int *iline, char *pad);
void  readme_docana_epoch(int *iline, char *pad); // MJD, Trest ...
void  readme_docana_instr(int *iline, char *pad);
void  readme_docana_hostlib(int *iline, char *pad);
void  readme_docana_modelPar(int *iline, char *pad) ;
void  readme_docana_rate(int *iline, char *pad) ;
void  readme_docana_cutwin(int *iline, char *pad) ;
void  readme_docana_searcheff(int *iline, char *pad);
void  readme_docana_mwebv(int *iline, char *pad);
void  readme_docana_misc(int *iline, char *pad);

void  readme_docana_comment(int *iline, char *comment); // write # <comment> 

void VERSION_INFO_load(int *iline, char *pad, char *keyName, char *comment,
                       int lenkey, bool isint,
                       int nval, double *val, double valmin, double valmax,
                       double val_noprint ) ;

void readme_docana_load_asymGauss(int *iline, char *pad, 
				  GENGAUSS_ASYM_DEF *GENGAUASS);

void readme_docana_load_expHalfGauss(int *iline, char *pad, 
				     GEN_EXP_HALFGAUSS_DEF *EXP_HALFGAUASS);

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//    legacy readme functions below
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

void   readme_doc_legacy(int iflag_readme);
void   readme_doc_SIMLIB(int *iline) ;
void   readme_doc_filterWarn(int *iline);
void   readme_doc_hostxt(int *iline, GEN_EXP_HALFGAUSS_DEF *GENPROFILE) ;
void   readme_doc_MWXT(int *iline);
void   readme_doc_NON1ASED(int *iline);
void   readme_doc_SIMSED(int *iline);
void   readme_doc_magSmear(int *iline);
void   readme_doc_nonLin(int *iline);
void   readme_doc_SALT2params(int *iline ) ;
void   readme_doc_GENPDF(int *iline ) ;
void   readme_doc_FIXMAG(int *iline ) ;
void   readme_doc_GENPERFECT(int *iline ) ;
void   readme_doc_FUDGES(int *iline) ;
void   readme_doc_mapFileList(int *iline) ;
void   readme_doc_mapFile(int *iline, char *KEY, char *FILENAME) ;
void   readme_doc_CUTWIN(int *iline) ;
void   readme_doc_TAKE_SPECTRUM(int *iline);
