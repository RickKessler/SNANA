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
 *   on a grid of wavelengths (rows) and detcor property (columns). 
 *   The function here load that grid into memory once (read_detctor_map), 
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
 * WARNING: if p-value or wavelength is outside grid range, program aborts !!!
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
void read_detcor_map(char *fileName, DETCOR_MAP_DEF *map, char *comment_string) {

  // PURPOSE:
  //   Load the handshaker detcor map from disk into memory, so every later
  //   call (interp_detcor_ratio, get_detcor_ratio, apply_detcor_to_filter) works
  //   from memory without touching the file again. Call once per map.
  //
  // INPUTS:
  //   fileName : path to the map, an SNANA text table written by
  //              handshaker-icemap / icecor.write_ice_map:
  //                 DOCUMENTATION: ... DOCUMENTATION_END:   (skipped)
  //                 VARNAMES: ROW WAVE R_<t0>nm R_<t1>nm ...
  //                 ROW: 1 4800.00 1.000000 1.016542 ...
  //              Any number of R_<thick>nm columns is allowed (including 1),
  //              but they must be in ascending thickness order.
  //
  //   *map:  structure to load with contents of fileName
  //
  //   *comment_string : human-readable text to identify map; e.g., "ICE_THICKNESS"
  //
  // OUTPUTS (filled into *map; arrays are malloc'ed here):
  //   map->FILENAME      : copy of fileName
  //   map->NWAVE         : number of wavelength rows
  //   map->WAVE[iw]      : wavelength of row iw [A], checked ascending
  //   map->NTHICK        : number of R_<thick>nm columns
  //   map->THICK[it]     : thickness of column it [nm], parsed from its name
  //   map->RATIO[it][iw] : R(WAVE[iw]; THICK[it])
  //   map->COMMENT_STRING: comment_string
  //
  // WHAT IT DOES:
  //   1. Uses SNANA's table reader: SNTABLE_NEVT counts rows, TABLEFILE_OPEN
  //      opens the file, SNTABLE_READPREP reads the VARNAMES list.
  //   2. For each column: WAVE -> map->WAVE; R_<thick>nm -> new map->RATIO
  //      array (thickness taken from the column name). Other columns
  //      (e.g. ROW) are ignored.
  //   3. SNTABLE_READ_EXEC reads all rows into those arrays, then closes
  //      the file.
  //   4. Sanity checks; any failure aborts through SNANA's errmsg:
  //      NWAVE >= 2, exactly 1 WAVE column, >= 1 R column,
  //      NTHICK <= MXTHICK_DETCOR, thickness and WAVE both ascending,
  //      rows read == rows counted.
  //
  // Release memory with free_detcor_map(map).

  int    NWAVE, NVAR, NROW, IFILETYPE, ivar, iw, NCOL_WAVE = 0 ;
  int    OPT_VARDEF = 0 ;
  double thick ;
  char  *varName ;
  char   TBLNAME[] = "DETCOR" ;
  char   fnam[]    = "read_detcor_map" ;

  // ------------ BEGIN -------------

  memset(map, 0, sizeof(*map));
  snprintf(map->FILENAME, sizeof(map->FILENAME), "%s", fileName);
  
  map->COMMENT_STRING = (char*) malloc( sizeof(char) * sizeof(comment_string) );
  sprintf(map->COMMENT_STRING, "%s", comment_string);

  TABLEFILE_INIT();
  NWAVE     = SNTABLE_NEVT(fileName, TBLNAME);
  IFILETYPE = TABLEFILE_OPEN(fileName, "read", fnam);
  NVAR      = SNTABLE_READPREP(IFILETYPE, TBLNAME);

  if ( NWAVE < 2 ) {
    sprintf(c1err,"Found NWAVE=%d rows; need at least 2", NWAVE);
    sprintf(c2err,"Check DETCOR file %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // map each VARNAME to an array: WAVE, or R_<thick>nm
  map->WAVE = (double*) malloc(NWAVE * sizeof(double));
  for(ivar=0; ivar < NVAR; ivar++ ) {
    varName = READTABLE_POINTERS.VARNAME[ivar];

    if ( strcmp(varName,"WAVE") == 0 ) {
      SNTABLE_READPREP_VARDEF(varName, map->WAVE, NWAVE, OPT_VARDEF);
      NCOL_WAVE++ ;
      continue ;
    }

    if ( sscanf(varName, "R_%lfnm", &thick) != 1 ) { continue; }

    if ( map->NTHICK >= MXTHICK_DETCOR ) {
      sprintf(c1err,"More than MXTHICK_DETCOR=%d thickness columns",
	      MXTHICK_DETCOR);
      sprintf(c2err,"in %s", fileName);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
    if ( map->NTHICK > 0 && thick <= map->THICK[map->NTHICK-1] ) {
      sprintf(c1err,"Thickness columns not ascending: %s after R_%gnm",
	      varName, map->THICK[map->NTHICK-1]);
      sprintf(c2err,"in %s", fileName);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }

    map->THICK[map->NTHICK] = thick ;
    map->RATIO[map->NTHICK] = (double*) malloc(NWAVE * sizeof(double));
    SNTABLE_READPREP_VARDEF(varName, map->RATIO[map->NTHICK], NWAVE, OPT_VARDEF);
    map->NTHICK++ ;
  }

  if ( NCOL_WAVE != 1 || map->NTHICK == 0 ) {
    sprintf(c1err,"Found %d WAVE columns and %d R_<thick>nm columns",
	    NCOL_WAVE, map->NTHICK);
    sprintf(c2err,"Need 1 WAVE and >=1 R_<thick>nm; check %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }

  // read the whole table (also closes the file)
  NROW = SNTABLE_READ_EXEC(0);
  if ( NROW != NWAVE ) {
    sprintf(c1err,"Read %d rows, but SNTABLE_NEVT counted %d", NROW, NWAVE);
    sprintf(c2err,"Check DETCOR file %s", fileName);
    errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
  }
  map->NWAVE = NWAVE ;

  for(iw=1; iw < NWAVE; iw++ ) {
    if ( map->WAVE[iw] <= map->WAVE[iw-1] ) {
      sprintf(c1err,"WAVE not ascending at row %d (%.2f after %.2f)",
	      iw+1, map->WAVE[iw], map->WAVE[iw-1]);
      sprintf(c2err,"Check DETCOR file %s", fileName);
      errmsg(SEV_FATAL, 0, fnam, c1err, c2err);
    }
  }

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
  //   map : WAVE and every RATIO[it] array freed; the whole struct is
  //         then zeroed (NWAVE = NTHICK = 0, pointers NULL), so it is safe
  //         to pass to read_detcor_map again.

  int it;
  free(map->WAVE);
  for(it=0; it < map->NTHICK; it++ ) { free(map->RATIO[it]); }
  memset(map, 0, sizeof(*map));
}


// ===================================================================
static void weights_detcor_map(DETCOR_MAP_DEF *map, double thick,
			       int *ithick, double *w) {

  // PURPOSE (internal helper, not in detcor_map.h):
  //   Find which two thickness columns bracket the requested thickness,
  //   and how far between them it is, so callers can interpolate linearly:
  //       R = (1-w) * RATIO[ithick] + w * RATIO[ithick+1]
  //
  // INPUTS:
  //   map   : map filled by read_detcor_map
  //   thick : requested detcor thickness [nm]. Values below THICK[0] or above
  //           THICK[NTHICK-1] are clamped to that end of the map (no
  //           extrapolation).
  //
  // OUTPUTS:
  //   *ithick : index of the lower bracketing column, 0 .. NTHICK-2
  //             (0 for a single-column map)
  //   *w      : weight of the upper column, 0 .. 1
  //             (0 for a single-column map)
  //
  // WHAT IT DOES:
  //   Clamps thick, finds the bin with SNANA's quickBinSearch, and keeps the
  //   index at most NTHICK-2 so that ithick+1 is always a valid column
  //   (needed when thick sits exactly on the last column).

  int    NTHICK = map->NTHICK ;
  double *THICK = map->THICK ;
  int    it ;
  char   fnam[] = "weights_detcor_map" ;

  if ( NTHICK == 1 ) { *ithick = 0; *w = 0.0; return; }

  if ( thick < THICK[0]        ) { thick = THICK[0];        }
  if ( thick > THICK[NTHICK-1] ) { thick = THICK[NTHICK-1]; }

  it = quickBinSearch(thick, NTHICK, THICK, "DETCOR thickness", fnam);
  if ( it > NTHICK-2 ) { it = NTHICK-2; }

  *ithick = it ;
  *w      = (thick - THICK[it]) / (THICK[it+1] - THICK[it]) ;

} // end weights_detcor_map


// ======================================================
double interp_detcor_ratio(DETCOR_MAP_DEF *map, double thick, double wave) {

  // PURPOSE:
  //   Return the detcor multiplier R at one arbitrary wavelength and
  //   thickness. This is the basic lookup the other functions build on.
  //
  // INPUTS:
  //   map   : map filled by read_detcor_map
  //   thick : detcor thickness [nm] (clamped to the map's thickness range)
  //   wave  : wavelength [A]; does not need to be on the map's WAVE grid
  //
  // OUTPUT (return value):
  //   R(wave; thick) = T(thick)/T(0) at that wavelength.
  //   Returns exactly 1.0 (no detcor effect) if wave is outside the map's
  //   WAVE range.
  //
  // WHAT IT DOES:
  //   1. Gets the two bracketing thickness columns and weight w from
  //      weights_detcor_map.
  //   2. Interpolates each column linearly in wavelength with SNANA's
  //      interp_1DFUN -> R0 (lower column), R1 (upper column).
  //   3. Returns (1-w)*R0 + w*R1. A single-column map returns R0.

  int    NWAVE = map->NWAVE, it ;
  double w, R0, R1 ;
  char   fnam[] = "interp_detcor_ratio" ;

  if ( wave < map->WAVE[0] || wave > map->WAVE[NWAVE-1] ) { return 1.0 ; }

  weights_detcor_map(map, thick, &it, &w);
  R0 = interp_1DFUN(OPT_INTERP_LINEAR, wave, NWAVE, map->WAVE,
		    map->RATIO[it], fnam);
  if ( map->NTHICK == 1 ) { return R0 ; }
  R1 = interp_1DFUN(OPT_INTERP_LINEAR, wave, NWAVE, map->WAVE,
		    map->RATIO[it+1], fnam);
  return (1.0-w)*R0 + w*R1 ;

} // end interp_detcor_ratio


// ======================================================
void get_detcor_ratio(DETCOR_MAP_DEF *map, double thick, double *ratio_out) {

  // PURPOSE:
  //   Return the whole multiplier curve R(lambda) for one thickness, on
  //   the map's own wavelength grid. It is the fast path (no wavelength
  //   interpolation) and is what the test compares against handshaker.
  //
  // INPUTS:
  //   map       : map filled by read_detcor_map
  //   thick     : detcor thickness [nm] (clamped to the map's thickness range)
  //
  // OUTPUTS:
  //   ratio_out : caller-allocated array of length map->NWAVE;
  //               ratio_out[iw] = R(map->WAVE[iw]; thick)
  //
  // WHAT IT DOES:
  //   Gets the bracketing columns and weight w once from
  //   weights_detcor_map, then for every wavelength row computes
  //   (1-w)*RATIO[it][iw] + w*RATIO[it+1][iw]. A single-column map
  //   uses the same column for both.

  int    it, it1, iw ;
  double w ;

  weights_detcor_map(map, thick, &it, &w);
  it1 = ( map->NTHICK > 1 ) ? it+1 : it ;
  for(iw=0; iw < map->NWAVE; iw++ )
    { ratio_out[iw] = (1.0-w)*map->RATIO[it][iw] + w*map->RATIO[it1][iw]; }

} // end get_detcor_ratio


// ======================================================
void apply_detcor_to_filter(DETCOR_MAP_DEF *map, double thick, int NLAM,
			 double *lam, double *trans_nom, double *trans_detcor) {

  // PURPOSE:
  //   Make the detector-corrected version of a filter. This is the end product:
  //   SNANA would use trans_detcor in place of the nominal transmission for an
  //   observation with this detcor thickness.
  //
  // INPUTS:
  //   map       : map filled by read_detcor_map
  //   thick     : detcor thickness [nm] (clamped to the map's thickness range)
  //   NLAM      : number of points in the filter
  //   lam       : [NLAM] filter wavelengths [A]; any grid, it does not need
  //               to match the map's WAVE grid
  //   trans_nom : [NLAM] nominal (uncorrected) filter transmission
  //
  // OUTPUTS:
  //   trans_detcor : [NLAM] caller-allocated;
  //               trans_detcor[ilam] = trans_nom[ilam] * R(lam[ilam]; thick)
  //               Points outside the map's WAVE range keep trans_nom (R = 1).
  //
  // WHAT IT DOES:
  //   Calls interp_detcor_ratio at every filter wavelength and multiplies.
  //   trans_detcor may be the same array as trans_nom (in-place update).

  int ilam;
  for(ilam=0; ilam < NLAM; ilam++ )
    { trans_detcor[ilam] = trans_nom[ilam] * interp_detcor_ratio(map, thick, lam[ilam]); }
}
