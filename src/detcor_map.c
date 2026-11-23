/* detcor_map.c -- see detcor_map.h
 *
 *            
 * MAIN PURPOSE
 *   A detector efficiency correction (detcor) that depends only on 
 *   wavelngth, but NOT on filter. Thus eff(detcor) in filter-overlap 
 *   regions is the same for each overlap band.
 *
 *   An external code tabulates the multiplicative correction factor
 *   vs. wavelength and vs. a detector property (p);
 *       R(lambda; p) = T(p) / T(p0)
 *   on a grid of wavelengths and detcor property (long-format rows). 
 *   The function here load that grid into memory once (read_detcor_map), 
 *   and during program execution get_detcor_ratio() interpolates R(lambda)
 *   [linear interp] for any p-value.
 *
 *   This correction can be used by SNANA sim to model the effect,
 *   and also by LC fit to correct for it (if survey team does not
 *   apply such a correction).
 *
 *   Initial motivation is to correct Roman detector efficiency (and thus
 *   filter transmission) for p = ice thickness (in sim), where the correction 
 *   grid is produced by handshaker. 
 * 
 * Out-of-range: p is clamped to the map range; R=1 outside the WAVE range.
 *
 * Units: wavelength in Angstrom
 *
 * user interface functions
 *   read_detcor_map          read map file into an DETCOR_MAP_DEF (call once):
 *   get_detcor_ratio         R at every map wavelength for one property value
 *
 * local utility functions:
 *   free_detcor_map          release the memory read_detcor_map allocated
 *   weights_detcor_map       property bin + interpolation weight
 *   interp_detcor_ratio      R at one wavelength and property value
 *   apply_detcor_to_filter   detector-corrected transmission on any filter grid
 */
#include "sntools.h"
#include "sntools_output.h"
#include "detcor_map.h"

// ======================================================
static int cmp_double_detcor(const void *a, const void *b) {
  double da = *(const double*)a, db = *(const double*)b;
  return (da > db) - (da < db);
}

// ======================================================
static int unique_sorted_detcor(int N, double *vals, double *uniq) {
  // Copy vals[N] into uniq, sort ascending, drop duplicates.
  // Returns number of unique values; uniq must hold N elements.
  int i, NU = 0 ;
  for(i=0; i < N; i++ ) { uniq[i] = vals[i]; }
  qsort(uniq, N, sizeof(double), cmp_double_detcor);
  for(i=0; i < N; i++ ) {
    if ( NU == 0 || uniq[i] != uniq[NU-1] ) { uniq[NU++] = uniq[i]; }
  }
  return NU ;
}

// ======================================================
static int exact_index_detcor(double val, int N, double *list) {
  // Binary search for val in ascending list[N]; return index or -1.
  int lo = 0, hi = N-1, mid ;
  while ( lo <= hi ) {
    mid = (lo + hi) / 2 ;
    if      ( list[mid] < val ) { lo = mid + 1; }
    else if ( list[mid] > val ) { hi = mid - 1; }
    else                        { return mid ;  }
  }
  return -1 ;
}

// ======================================================
void read_detcor_map(char *fileName, DETCOR_MAP_DEF *map,
		     char *column_name, char *ratio_name) {

  // PURPOSE:
  //   Load a LONG-format detcor map from disk into memory, so every later
  //   call (interp_detcor_ratio, get_detcor_ratio, apply_detcor_to_filter)
  //   works from memory without touching the file again. Call once per map.
  //
  // INPUTS:
  //   fileName    : path to the map, an SNANA text table, e.g.
  //                    VARNAMES: ROW WAVE ICE_THICKNESS RATIO_MODEL1 ...
  //                    ROW: 1 4800.00 0.0 1.000000 ...
  //                 One row per (WAVE, property) pair; rows in any order,
  //                 but together they must fill a complete grid.
  //   *map        : structure to load with contents of fileName
  //   column_name : name of property column; e.g., "ICE_THICKNESS"
  //   ratio_name  : name of ratio column to use; e.g., "RATIO_MODEL1"
  //
  // OUTPUTS (filled into *map; arrays are malloc'ed here):
  //   map->FILENAME       : copy of fileName
  //   map->COLUMN_NAME    : copy of column_name
  //   map->RATIO_NAME     : copy of ratio_name
  //   map->NWAVE, WAVE[]  : unique wavelengths [A], ascending
  //   map->NPROP, PROP[]  : unique property values, ascending
  //   map->RATIO[ip][iw]  : R(WAVE[iw]; PROP[ip])
  //
  // WHAT IT DOES:
  //   1. Reads WAVE, <column_name>, <ratio_name> by name (abort if missing);
  //      other columns (ROW, other RATIO models) are ignored.
  //   2. Builds sorted unique WAVE and property lists.
  //   3. Checks NWAVE*NPROP == NROW, then places each row in the grid and
  //      aborts on a duplicate (WAVE,property) pair -> grid is complete.
  //
  // Release memory with free_detcor_map(map).

  int    NROW, NROW_READ, NWAVE, NPROP, IFILETYPE, irow, iw, ip ;
  int    OPT_VARDEF = 2 ;  // abort if column is missing
  double *wave_row, *prop_row, *ratio_row, *uniq ;
  char   **filled ;
  char   TBLNAME[] = "DETCOR" ;
  char   fnam[]    = "read_detcor_map" ;

  // ------------ BEGIN -------------

  memset(map, 0, sizeof(*map));
  snprintf(map->FILENAME, sizeof(map->FILENAME), "%s", fileName);

  map->COLUMN_NAME = (char*) malloc( strlen(column_name) + 1 );
  sprintf(map->COLUMN_NAME, "%s", column_name);
  map->RATIO_NAME  = (char*) malloc( strlen(ratio_name) + 1 );
  sprintf(map->RATIO_NAME, "%s", ratio_name);

  TABLEFILE_INIT();
  NROW      = SNTABLE_NEVT(fileName, TBLNAME);
  IFILETYPE = TABLEFILE_OPEN(fileName, "read", fnam);
  SNTABLE_READPREP(IFILETYPE, TBLNAME);

  if ( NROW < 2 ) {
    sprintf(c1err,"Found NROW=%d rows; need at least 2", NROW);
    sprintf(c2err,"Check DETCOR file %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // read the 3 needed columns, row by row as stored in the file
  wave_row  = (double*) malloc(NROW * sizeof(double));
  prop_row  = (double*) malloc(NROW * sizeof(double));
  ratio_row = (double*) malloc(NROW * sizeof(double));
  SNTABLE_READPREP_VARDEF("WAVE",      wave_row,  NROW, OPT_VARDEF);
  SNTABLE_READPREP_VARDEF(column_name, prop_row,  NROW, OPT_VARDEF);
  SNTABLE_READPREP_VARDEF(ratio_name,  ratio_row, NROW, OPT_VARDEF);

  NROW_READ = SNTABLE_READ_EXEC(0);  // also closes the file
  if ( NROW_READ != NROW ) {
    sprintf(c1err,"Read %d rows, but SNTABLE_NEVT counted %d", NROW_READ, NROW);
    sprintf(c2err,"Check DETCOR file %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // unique sorted WAVE and property grids
  uniq  = (double*) malloc(NROW * sizeof(double));

  NWAVE = unique_sorted_detcor(NROW, wave_row, uniq);
  map->WAVE = (double*) malloc(NWAVE * sizeof(double));
  for(iw=0; iw < NWAVE; iw++ ) { map->WAVE[iw] = uniq[iw]; }

  NPROP = unique_sorted_detcor(NROW, prop_row, uniq);
  if ( NPROP > MXPROP_DETCOR ) {
    sprintf(c1err,"Found %d unique %s values; MXPROP_DETCOR=%d",
	    NPROP, column_name, MXPROP_DETCOR);
    sprintf(c2err,"in %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  for(ip=0; ip < NPROP; ip++ ) { map->PROP[ip] = uniq[ip]; }
  free(uniq);

  if ( NWAVE < 2 ) {
    sprintf(c1err,"Found NWAVE=%d unique wavelengths; need at least 2", NWAVE);
    sprintf(c2err,"Check DETCOR file %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  if ( NWAVE*NPROP != NROW ) {
    sprintf(c1err,"NWAVE(%d) x N%s(%d) = %d, but NROW=%d -> incomplete grid",
	    NWAVE, column_name, NPROP, NWAVE*NPROP, NROW);
    sprintf(c2err,"Check DETCOR file %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // place each row in RATIO[ip][iw]; flag duplicates
  filled = (char**) malloc(NPROP * sizeof(char*));
  for(ip=0; ip < NPROP; ip++ ) {
    map->RATIO[ip] = (double*) malloc(NWAVE * sizeof(double));
    filled[ip]     = (char*)   calloc(NWAVE, sizeof(char));
  }

  for(irow=0; irow < NROW; irow++ ) {
    iw = exact_index_detcor(wave_row[irow], NWAVE, map->WAVE);
    ip = exact_index_detcor(prop_row[irow], NPROP, map->PROP);
    if ( filled[ip][iw] ) {
      sprintf(c1err,"Duplicate row for WAVE=%.2f %s=%g",
	      wave_row[irow], column_name, prop_row[irow]);
      sprintf(c2err,"Check DETCOR file %s", fileName);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
    map->RATIO[ip][iw] = ratio_row[irow];
    filled[ip][iw]     = 1 ;
  }

  map->NWAVE = NWAVE ;
  map->NPROP = NPROP ;

  for(ip=0; ip < NPROP; ip++ ) { free(filled[ip]); }
  free(filled);
  free(wave_row); free(prop_row); free(ratio_row);

  printf("\t Loaded DETCOR map %s: NWAVE=%d (%.1f-%.1f A), "
	 "N%s=%d (%g-%g), ratio column %s\n",
	 fileName, NWAVE, map->WAVE[0], map->WAVE[NWAVE-1],
	 column_name, NPROP, map->PROP[0], map->PROP[NPROP-1], ratio_name);
  fflush(stdout);

  return ;

} // end read_detcor_map


// ======================================================
void free_detcor_map(DETCOR_MAP_DEF *map) {

  // PURPOSE:
  //   Release the memory read_detcor_map allocated, so a map can be
  //   re-read or discarded without leaking.
  //
  // INPUTS:
  //   map : a map previously filled by read_detcor_map
  //
  // OUTPUTS:
  //   map : WAVE, every RATIO[ip] array and the name strings freed; the
  //         whole struct is then zeroed (NWAVE = NPROP = 0, pointers NULL),
  //         so it is safe to pass to read_detcor_map again.

  int ip;
  free(map->WAVE);
  for(ip=0; ip < map->NPROP; ip++ ) { free(map->RATIO[ip]); }
  free(map->COLUMN_NAME);
  free(map->RATIO_NAME);
  memset(map, 0, sizeof(*map));
}


// ===================================================================
static void weights_detcor_map(DETCOR_MAP_DEF *map, double prop,
			       int *iprop, double *w) {

  // PURPOSE (internal helper, not in detcor_map.h):
  //   Find which two property columns bracket the requested property,
  //   and how far between them it is, so callers can interpolate linearly:
  //       R = (1-w) * RATIO[iprop] + w * RATIO[iprop+1]
  //
  // INPUTS:
  //   map   : map filled by read_detcor_map
  //   prop : requested property value, e.g. ice thickness [nm]. Values below PROP[0] or above
  //           PROP[NPROP-1] are clamped to that end of the map (no
  //           extrapolation).
  //
  // OUTPUTS:
  //   *iprop : index of the lower bracketing column, 0 .. NPROP-2
  //             (0 for a single-property map)
  //   *w      : weight of the upper column, 0 .. 1
  //             (0 for a single-property map)
  //
  // WHAT IT DOES:
  //   Clamps prop, finds the bin with SNANA's quickBinSearch, and keeps the
  //   index at most NPROP-2 so that iprop+1 is always a valid column
  //   (needed when prop sits exactly on the last column).

  int    NPROP = map->NPROP ;
  double *PROP = map->PROP ;
  int    it ;
  char   fnam[] = "weights_detcor_map" ;

  if ( NPROP == 1 ) { *iprop = 0; *w = 0.0; return; }

  if ( prop < PROP[0]        ) { prop = PROP[0];        }
  if ( prop > PROP[NPROP-1] ) { prop = PROP[NPROP-1]; }

  it = quickBinSearch(prop, NPROP, PROP, "DETCOR property", fnam);
  if ( it > NPROP-2 ) { it = NPROP-2; }

  *iprop = it ;
  *w      = (prop - PROP[it]) / (PROP[it+1] - PROP[it]) ;

} // end weights_detcor_map


// ======================================================
double interp_detcor_ratio(DETCOR_MAP_DEF *map, double prop, double wave) {

  // PURPOSE:
  //   Return the detcor multiplier R at one arbitrary wavelength and
  //   property. This is the basic lookup the other functions build on.
  //
  // INPUTS:
  //   map   : map filled by read_detcor_map
  //   prop : property value, e.g. ice thickness [nm] (clamped to the map's property range)
  //   wave  : wavelength [A]; does not need to be on the map's WAVE grid
  //
  // OUTPUT (return value):
  //   R(wave; prop) = T(prop)/T(p0) at that wavelength.
  //   Returns exactly 1.0 (no detcor effect) if wave is outside the map's
  //   WAVE range.
  //
  // WHAT IT DOES:
  //   1. Gets the two bracketing property columns and weight w from
  //      weights_detcor_map.
  //   2. Interpolates each column linearly in wavelength with SNANA's
  //      interp_1DFUN -> R0 (lower column), R1 (upper column).
  //   3. Returns (1-w)*R0 + w*R1. A single-property map returns R0.

  int    NWAVE = map->NWAVE, it ;
  double w, R0, R1 ;
  char   fnam[] = "interp_detcor_ratio" ;

  if ( wave < map->WAVE[0] || wave > map->WAVE[NWAVE-1] ) { return 1.0 ; }

  weights_detcor_map(map, prop, &it, &w);
  R0 = interp_1DFUN(OPT_INTERP_LINEAR, wave, NWAVE, map->WAVE,
		    map->RATIO[it], fnam);
  if ( map->NPROP == 1 ) { return R0 ; }
  R1 = interp_1DFUN(OPT_INTERP_LINEAR, wave, NWAVE, map->WAVE,
		    map->RATIO[it+1], fnam);
  return (1.0-w)*R0 + w*R1 ;

} // end interp_detcor_ratio


// ======================================================
void get_detcor_ratio(DETCOR_MAP_DEF *map, double prop, double *ratio_out) {

  // PURPOSE:
  //   Return the whole multiplier curve R(lambda) for one property, on
  //   the map's own wavelength grid. It is the fast path (no wavelength
  //   interpolation) and is what the test compares against handshaker.
  //
  // INPUTS:
  //   map       : map filled by read_detcor_map
  //   prop     : property value, e.g. ice thickness [nm] (clamped to the map's property range)
  //
  // OUTPUTS:
  //   ratio_out : caller-allocated array of length map->NWAVE;
  //               ratio_out[iw] = R(map->WAVE[iw]; prop)
  //
  // WHAT IT DOES:
  //   Gets the bracketing columns and weight w once from
  //   weights_detcor_map, then for every wavelength row computes
  //   (1-w)*RATIO[it][iw] + w*RATIO[it+1][iw]. A single-property map
  //   uses the same column for both.

  int    it, it1, iw ;
  double w ;

  weights_detcor_map(map, prop, &it, &w);
  it1 = ( map->NPROP > 1 ) ? it+1 : it ;
  for(iw=0; iw < map->NWAVE; iw++ )
    { ratio_out[iw] = (1.0-w)*map->RATIO[it][iw] + w*map->RATIO[it1][iw]; }

} // end get_detcor_ratio


// ======================================================
void apply_detcor_to_filter(DETCOR_MAP_DEF *map, double prop, int NLAM,
			    double *lam, double *trans_nom, double *trans_detcor) {

  // PURPOSE:
  //   Make the detector-corrected version of a filter. This is the end product:
  //   SNANA would use trans_detcor in place of the nominal transmission for an
  //   observation with this detcor property.
  //
  // INPUTS:
  //   map       : map filled by read_detcor_map
  //   prop     : property value, e.g. ice thickness [nm] (clamped to the map's property range)
  //   NLAM      : number of points in the filter
  //   lam       : [NLAM] filter wavelengths [A]; any grid, it does not need
  //               to match the map's WAVE grid
  //   trans_nom : [NLAM] nominal (uncorrected) filter transmission
  //
  // OUTPUTS:
  //   trans_detcor : [NLAM] caller-allocated;
  //               trans_detcor[ilam] = trans_nom[ilam] * R(lam[ilam]; prop)
  //               Points outside the map's WAVE range keep trans_nom (R = 1).
  //
  // WHAT IT DOES:
  //   Calls interp_detcor_ratio at every filter wavelength and multiplies.
  //   trans_detcor may be the same array as trans_nom (in-place update).

  int ilam;
  for(ilam=0; ilam < NLAM; ilam++ )
    { trans_detcor[ilam] = trans_nom[ilam] * interp_detcor_ratio(map, prop, lam[ilam]); }
}
