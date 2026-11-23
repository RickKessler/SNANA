/* detcor_map.h -- Roman detector-correction (detcor) transmission-multiplier map, built on SNANA tools.
 *
 * Developed OUTSIDE SNANA (handshaker/snana/detcor_map). It calls SNANA's
 * table reader (SNTABLE_*), errmsg, quickBinSearch and interp_1DFUN, so it is
 * compiled with -I$SNANA_DIR/src and linked against SNANA's sntools.o and
 * sntools_output.o. SNANA programs can call these functions directly later.
 *
 * Map file = handshaker's `handshaker-icemap` / icecor.write_ice_map output,
 * a standard SNANA text table:
 *
 *   DOCUMENTATION: ... DOCUMENTATION_END:
 *   VARNAMES: ROW WAVE R_0nm R_10nm ... R_300nm
 *   ROW: 1 4800.00 1.000000 1.016542 ...
 *
 * R(lambda; t) = T(t)/T(0). The detector-corrected filter is
 *   T_detcor(lambda) = T_nom(lambda) * R(lambda; ICE_THICK),
 * with R linearly interpolated in thickness (clamped to the map's range) and in
 * wavelength (R = 1 outside the map's WAVE range) -- the same convention as
 * handshaker's ReleasedIceModel.
 */
#ifndef DETCOR_MAP_H
#define DETCOR_MAP_H

#define MXTHICK_DETCOR 200

typedef struct {
  char    FILENAME[1024];
  int     NWAVE, NTHICK;
  double *WAVE;                     /* [NWAVE]  Angstrom, ascending        */
  double  THICK[MXTHICK_DETCOR];    /* [NTHICK] nm, ascending              */
  double *RATIO[MXTHICK_DETCOR];    /* RATIO[ithick][iwave]                */
} DETCOR_MAP_DEF;

/* Read map file into *map (allocates). Aborts via errmsg on any problem. */
void   read_detcor_map(char *fileName, DETCOR_MAP_DEF *map);
void   free_detcor_map(DETCOR_MAP_DEF *map);

/* R(wave; thick) at one wavelength [A]; 1 outside the map's WAVE range. */
double interp_detcor_ratio(DETCOR_MAP_DEF *map, double thick_nm, double wave_A);

/* R(WAVE[iw]; thick) for all map wavelengths -> ratio_out[NWAVE]. */
void   get_detcor_ratio(DETCOR_MAP_DEF *map, double thick_nm, double *ratio_out);

/* trans_detcor[i] = trans_nom[i] * R(lam[i]; thick), for any filter grid. */
void   apply_detcor_to_filter(DETCOR_MAP_DEF *map, double thick_nm, int NLAM,
                           double *lam_A, double *trans_nom, double *trans_detcor);

#endif
