/*********************************************
  Created March 2016

 *******************************************/

// globals
SNGRID_DEF NON1AGRID ;    // read from file

double LOGZ_NON1AGRID ;
int    ILOGZ_NON1AGRID ;
int    INDEX_NON1AGRID ;

// -------- prototype functions ------------

void init_genmag_NON1AGRID(char *GRIDFILE, double FRAC_PEC1A );

void init_interp_NON1AGRIDMAP(int ifilt, SNGRID_DEF *SNGRID, GRIDMAP *GRIDMAP );

void genmag_NON1AGRID (int ifilt_obs, double mwebv, double z,
		       double RVhost, double AVhost,
		       double ranWgt, double ranSmear,
		       int NOBS, double *TobsList,
		       double *magList, double *magerrList, double *magSmear);

double  magInterp_NON1AGRID(int ifilt, int NON1A_INDEX, double z, double Trest) ;

double magNode_NON1AGRID(int ifilt, int NON1A_INDEX, int iz, int ep) ;

double fetchInfo_NON1AGRID(char*what) ;

void  checkRange_NON1AGRID(int IPAR, double VAL) ;

// end
