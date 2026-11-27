/* detcor_map.h -- detector-correction (detcor) transmission-multiplier map,
 * built on SNANA tools.
 *
 * Map file = SNANA text table in LONG format (one row per wavelength x
 * property value), e.g. written by make_detcor_map.py:
 *
 *   DOCUMENTATION: ... DOCUMENTATION_END:
 *   VARNAMES: ROW WAVE ICE_THICKNESS RATIO_MODEL1 [RATIO_MODEL2 ...]
 *   ROW: 1  4800.00    0.0  1.000000
 *   ROW: 2  4800.00   10.0  1.016542
 *
 * The property column (e.g. ICE_THICKNESS) and the ratio column
 * (e.g. RATIO_MODEL1) are selected by name from the sim-input key
 *   DETCOR_MJD:  <COLUMN_NAME>  <RATIO_NAME>  <MJD_FILE>  <DETCOR_FILE>
 * so several models can live in one file. Rows may be in any order, but
 * must fill a complete (WAVE x property) grid.
 *
 * R(lambda; p) = T(p)/T(p0). The detector-corrected filter is
 *   T_detcor(lambda) = T_nom(lambda) * R(lambda; p),
 * with R linearly interpolated in p (clamped to the map's range) and in
 * wavelength (R = 1 outside the map's WAVE range) -- the same convention as
 * handshaker's ReleasedIceModel.
 */

#ifndef DETCOR_MAP_H
#define DETCOR_MAP_H

#define MXPROP_DETCOR 200

typedef struct {
  char    FILENAME[1024];
  int     NWAVE, NPROP;
  double *WAVE;                    // [NWAVE] Angstrom, ascending
  double  PROP[MXPROP_DETCOR];     // [NPROP] property values (e.g. ice thickness, nm), ascending
  double *RATIO[MXPROP_DETCOR];    // RATIO[iprop][iwave]
  char   *COLUMN_NAME ;            // property column name in map file; e.g. "ICE_THICKNESS"
  char   *RATIO_NAME ;             // ratio (model) column name;        e.g. "RATIO_THUSHARA"
} DETCOR_MAP_DEF;

/* Property value vs. MJD (e.g., ice thickness history), read from MJD_FILE:
 *
 *   VARNAMES: ROW MJD ICE_THICKNESS
 *   ROW: 1  61000.00000    0.0
 *   ROW: 2  61001.00000   12.0
 *
 * The property column is selected by the same COLUMN_NAME as the map.
 * MJD must be strictly ascending; the property is linearly interpolated
 * in MJD, so a sawtooth (e.g. decon cycle) is exact if the file lists its
 * breakpoints. An MJD outside the table range aborts.
 */
typedef struct {
  char    FILENAME[1024];
  int     NMJD;
  double *MJD;                     // [NMJD] strictly ascending
  double *PROP;                    // [NMJD] property value at each MJD
  char   *COLUMN_NAME ;            // property column name; e.g. "ICE_THICKNESS"
} DETCOR_MJD_DEF;

DETCOR_MAP_DEF DETCOR_SEDMODEL ; // Nov 24 2026 RK 
DETCOR_MJD_DEF DETCOR_MJD_TABLE ; // property vs. MJD from MJD_FILE (Nov 2026)

/* Read long-format map file into *map (allocates). Aborts via errmsg on any problem. */
void   READ_DETCOR_MAP(char *fileName, DETCOR_MAP_DEF *map,
                       char *column_name, char *ratio_name);
void   free_detcor_map(DETCOR_MAP_DEF *map);

/* Read property-vs-MJD file into *tab (allocates); free with free_detcor_mjd. */
void   READ_DETCOR_MJD(char *fileName, DETCOR_MJD_DEF *tab, char *column_name);
void   free_detcor_mjd(DETCOR_MJD_DEF *tab);

/* Property value at mjd (linear interp); aborts if mjd outside table range. */
double GET_DETCOR_PROP(DETCOR_MJD_DEF *tab, double mjd);

/* R(wave; prop) at one wavelength [A]; 1 outside the map's WAVE range. */
double interp_detcor_ratio(DETCOR_MAP_DEF *map, double prop, double wave_A);

/* R(WAVE[iw]; prop) for all map wavelengths -> ratio_out[NWAVE]. */
void   GET_DETCOR_RATIO(DETCOR_MAP_DEF *map, double prop, double *ratio_out);

/* trans_detcor[i] = trans_nom[i] * R(lam[i]; prop), for any filter grid. */
void   apply_detcor_to_filter(DETCOR_MAP_DEF *map, double prop, int NLAM,
                              double *lam_A, double *trans_nom, double *trans_detcor);

#endif
