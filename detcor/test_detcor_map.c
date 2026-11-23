/* test_detcor_map.c -- exercise detcor_map.c.
 *
 * Usage:
 *   test_detcor_map <map.dat> [--mode timing|check] [--ntest N] [--seed S]
 *                [--tmax T] [--filter FILE]... [--dump OUT] [--python "CMD"]
 *
 * --mode timing (default)
 *     time read_detcor_map(), then N random thicknesses in [0, tmax] through
 *     get_detcor_ratio() (and apply_detcor_to_filter() for every --filter given).
 *     Dumps NDUMP random-thickness multipliers to OUT, and for each --filter the
 *     detector-corrected filter T_nom*R at the same thicknesses to OUT_<filterbase>.
 * --mode check
 *     for N random thicknesses, system-call handshaker for that single thickness
 *     (CMD --thickness T -o tmp), read its output, and compare to the C
 *     interpolation of the grid. PASS if max|dR| <= 2e-6 (both are %.6f text).
 *     Prints the python wall time next to the C time for the same thickness.
 * --mode profile
 *     look for bottlenecks (see do_profile): scaling with number of calls,
 *     memory growth, per-thickness sweep (written to profile_thickness.dat),
 *     cost of each SNANA tool inside one lookup vs a hand-coded reference,
 *     and repeated map reloads. --ntest sets the largest scaling N (default 1e6).
 *
 * CMD defaults to $HANDSHAKER_DETCOR, else "handshaker-icemap".
 * --filter FILE: 2 columns (wave[A] trans), read with SNANA's rd2columnFile.
 *                Repeat --filter for several bands (max MXFILT).
 *
 * Built against SNANA objects (sntools.o, sntools_output.o): see Makefile.
 */
#include "sntools.h"
#include "detcor_map.h"
#include <sys/resource.h>

#define NDUMP  5
#define TOL    2.0e-6
#define MXFILT 32
#define MXLAM_FILT 20000

typedef struct {
  char       *file;
  char        base[256];   /* filename without directory or extension */
  int         NLAM;
  double     *lam, *trans;
} FILTER_DEF;

static double now_ms(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return 1.0e3 * ts.tv_sec + 1.0e-6 * ts.tv_nsec;
}

static double urand(double lo, double hi) { return lo + (hi - lo) * (rand() / (double)RAND_MAX); }

static void read_filter(char *file, FILTER_DEF *f)
{
  char *p;
  f->file = file;
  p = strrchr(file, '/');
  snprintf(f->base, sizeof(f->base), "%s", p ? p + 1 : file);
  if ((p = strrchr(f->base, '.')) != NULL) *p = '\0';
  f->lam   = malloc(sizeof(double) * MXLAM_FILT);
  f->trans = malloc(sizeof(double) * MXLAM_FILT);
  rd2columnFile(file, MXLAM_FILT, &f->NLAM, f->lam, f->trans, 0);  /* SNANA; aborts on error */
  if (f->NLAM < 2) {
    sprintf(c1err, "filter has %d points", f->NLAM);
    sprintf(c2err, "check %s", file);
    errmsg(SEV_FATAL, 0, "read_filter", c1err, c2err);
  }
}

static int dump_filter(DETCOR_MAP_DEF *map, const FILTER_DEF *f, const double *thk, const char *out)
{
  int i, iw;
  FILE *fp;
  double *cols = malloc(sizeof(double) * NDUMP * f->NLAM);
  if ((fp = fopen(out, "w")) == NULL) { fprintf(stderr, "ERROR: cannot write %s\n", out); free(cols); return 1; }
  for (i = 0; i < NDUMP; i++) apply_detcor_to_filter(map, thk[i], f->NLAM, f->lam, f->trans, cols + i * f->NLAM);
  fprintf(fp, "# detector-corrected filter T_nom*R(lambda;t)\n# filter: %s\n# map:    %s\n", f->file, map->FILENAME);
  fprintf(fp, "# WAVE[A] T_nom");
  for (i = 0; i < NDUMP; i++) fprintf(fp, " t=%.2fnm", thk[i]);
  fprintf(fp, "\n");
  for (iw = 0; iw < f->NLAM; iw++) {
    fprintf(fp, "%.2f %.6f", f->lam[iw], f->trans[iw]);
    for (i = 0; i < NDUMP; i++) fprintf(fp, " %.6f", cols[i * f->NLAM + iw]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  free(cols);
  return 0;
}

static int do_timing(DETCOR_MAP_DEF *map, int ntest, double tmax, FILTER_DEF *filt, int nfilt, const char *dump)
{
  int i, iw, k, rc = 0;
  double *R = malloc(sizeof(double) * map->NWAVE);
  double thk[NDUMP], t0, dt, sink = 0.0;
  FILE *fp;

  t0 = now_ms();
  for (i = 0; i < ntest; i++) { get_detcor_ratio(map, urand(0, tmax), R); sink += R[map->NWAVE / 2]; }
  dt = now_ms() - t0;
  printf("  get_detcor_ratio      : %d calls, %.3f ms total, %.3f us/call (NWAVE=%d)\n",
         ntest, dt, 1e3 * dt / ntest, map->NWAVE);

  for (k = 0; k < nfilt; k++) {
    double *tout = malloc(sizeof(double) * filt[k].NLAM);
    t0 = now_ms();
    for (i = 0; i < ntest; i++) {
      apply_detcor_to_filter(map, urand(0, tmax), filt[k].NLAM, filt[k].lam, filt[k].trans, tout);
      sink += tout[filt[k].NLAM / 2];
    }
    dt = now_ms() - t0;
    printf("  apply_detcor_to_filter: %d calls, %.3f ms total, %.3f us/call (NLAM=%4d, %s)\n",
           ntest, dt, 1e3 * dt / ntest, filt[k].NLAM, filt[k].base);
    free(tout);
  }

  /* dump NDUMP random thicknesses: multiplier R, then each detector-corrected filter */
  for (i = 0; i < NDUMP; i++) thk[i] = urand(0, tmax);
  if ((fp = fopen(dump, "w")) == NULL) { fprintf(stderr, "ERROR: cannot write %s\n", dump); free(R); return 1; }
  fprintf(fp, "# detcor multiplier R(lambda;t) from %s\n# WAVE[A]", map->FILENAME);
  for (i = 0; i < NDUMP; i++) fprintf(fp, " t=%.2fnm", thk[i]);
  fprintf(fp, "\n");
  {
    double *cols = malloc(sizeof(double) * NDUMP * map->NWAVE);
    for (i = 0; i < NDUMP; i++) get_detcor_ratio(map, thk[i], cols + i * map->NWAVE);
    for (iw = 0; iw < map->NWAVE; iw++) {
      fprintf(fp, "%.2f", map->WAVE[iw]);
      for (i = 0; i < NDUMP; i++) fprintf(fp, " %.6f", cols[i * map->NWAVE + iw]);
      fprintf(fp, "\n");
    }
    free(cols);
  }
  fclose(fp);
  printf("  dumped R at %d random thicknesses -> %s\n", NDUMP, dump);

  for (k = 0; k < nfilt; k++) {
    char out[1024];
    snprintf(out, sizeof(out), "%s_%s", dump, filt[k].base);
    rc |= dump_filter(map, &filt[k], thk, out);
    printf("  dumped %-10s at the same thicknesses -> %s\n", filt[k].base, out);
  }
  printf("  (checksum %.6f)\n", sink);
  free(R);
  return rc;
}

static int do_check(DETCOR_MAP_DEF *map, int ntest, double tmax, const char *pycmd)
{
  int i, iw, nfail = 0;
  /* .dat suffix so SNANA's table reader recognizes a text table */
  char cmd[2048], tmp[] = "/tmp/detcor_map_check_XXXXXX.dat";
  double *R = malloc(sizeof(double) * map->NWAVE), sum_py = 0.0, sum_c = 0.0;
  int fd = mkstemps(tmp, 4);
  if (fd < 0) { perror("mkstemps"); free(R); return 1; }
  close(fd);

  printf("  python: %s\n", pycmd);
  for (i = 0; i < ntest; i++) {
    DETCOR_MAP_DEF py;
    double t = urand(0, tmax), dmax = 0.0, t0, dt_py, dt_c;
    snprintf(cmd, sizeof(cmd), "%s --thickness %.4f -o %s > /dev/null", pycmd, t, tmp);
    t0 = now_ms();
    if (system(cmd) != 0) { fprintf(stderr, "ERROR: '%s' failed\n", cmd); nfail++; continue; }
    dt_py = now_ms() - t0;
    read_detcor_map(tmp, &py);          /* SNANA reader; aborts on a bad file */
    if (py.NWAVE != map->NWAVE || py.NTHICK != 1 || fabs(py.THICK[0] - t) > 1e-3) {
      fprintf(stderr, "ERROR: handshaker output for t=%.4f has NWAVE=%d NTHICK=%d THICK=%g\n",
              t, py.NWAVE, py.NTHICK, py.NTHICK ? py.THICK[0] : -1.0);
      free_detcor_map(&py); nfail++; continue;
    }
    t0 = now_ms();
    get_detcor_ratio(map, t, R);
    dt_c = now_ms() - t0;
    for (iw = 0; iw < map->NWAVE; iw++) {
      double d = fabs(R[iw] - py.RATIO[0][iw]);
      if (fabs(map->WAVE[iw] - py.WAVE[iw]) > 0.01) d = 1e9;   /* wave grids must match */
      if (d > dmax) dmax = d;
    }
    printf("  t=%8.4f nm  max|R_C - R_handshaker| = %.2e  %s   python %7.1f ms | C %.4f ms\n",
           t, dmax, dmax <= TOL ? "PASS" : "FAIL", dt_py, dt_c);
    sum_py += dt_py; sum_c += dt_c;
    if (dmax > TOL) nfail++;
    free_detcor_map(&py);
  }
  remove(tmp);
  free(R);
  if (ntest > nfail)
    printf("  mean time per thickness: python %.1f ms, C %.4f ms\n",
           sum_py / (ntest - nfail), sum_c / (ntest - nfail));
  printf("  check: %d/%d PASS\n", ntest - nfail, ntest);
  return nfail ? 1 : 0;
}

/* ===================== profile mode ===================== */

static double maxrss_MB(void)
{
  /* peak resident memory of this process so far [MB] */
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
  return ru.ru_maxrss / (1024.0 * 1024.0);   /* bytes on macOS */
#else
  return ru.ru_maxrss / 1024.0;              /* kB on Linux    */
#endif
}

static int cmp_dbl(const void *a, const void *b)
{
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static double median(double *v, int n) { qsort(v, n, sizeof(double), cmp_dbl); return v[n / 2]; }

static int bsearch_lo(const double *x, int n, double v)
{
  /* largest i in [0, n-2] with x[i] <= v (x ascending) */
  int lo = 0, hi = n - 1;
  if (v <= x[0])     return 0;
  if (v >= x[n - 1]) return n - 2;
  while (hi - lo > 1) { int m = (lo + hi) / 2; if (x[m] <= v) lo = m; else hi = m; }
  return lo;
}

static double ref_ratio(DETCOR_MAP_DEF *map, double thick, double wave)
{
  /* hand-coded reference: same math as interp_detcor_ratio, no SNANA calls.
   * One wavelength search shared by both thickness columns. */
  int n = map->NWAVE, it, iw, it1;
  double w, u, a, b;
  if (wave < map->WAVE[0] || wave > map->WAVE[n - 1]) return 1.0;
  if (map->NTHICK == 1) { it = 0; w = 0.0; }
  else {
    double t = thick;
    if (t < map->THICK[0]) t = map->THICK[0];
    if (t > map->THICK[map->NTHICK - 1]) t = map->THICK[map->NTHICK - 1];
    it = bsearch_lo(map->THICK, map->NTHICK, t);
    w  = (t - map->THICK[it]) / (map->THICK[it + 1] - map->THICK[it]);
  }
  it1 = (map->NTHICK > 1) ? it + 1 : it;
  iw = bsearch_lo(map->WAVE, n, wave);
  u  = (wave - map->WAVE[iw]) / (map->WAVE[iw + 1] - map->WAVE[iw]);
  a  = (1.0 - u) * map->RATIO[it][iw]  + u * map->RATIO[it][iw + 1];
  b  = (1.0 - u) * map->RATIO[it1][iw] + u * map->RATIO[it1][iw + 1];
  return (1.0 - w) * a + w * b;
}

static int do_profile(DETCOR_MAP_DEF *map, int nmax, double tmax, FILTER_DEF *filt, int nfilt)
{
  /* Look for bottlenecks / slow cases. Sections:
   *  P1 timer overhead
   *  P2 scaling: cost per call vs number of calls (should be flat) + memory
   *  P3 thickness sweep: every thickness -30 .. tmax+30 nm in 1 nm steps,
   *     median cost per call at each (slow spots? clamping? exact grid nodes?)
   *  P4 breakdown of one lookup into SNANA tools vs hand-coded reference
   *  P5 repeated read_detcor_map/free_detcor_map: time + memory growth (leaks)
   */
  FILTER_DEF syn, *F;
  int i, k, N, nsweep, M = 101, nbad = 0;
  double t0, dt, sink = 0.0, rss0, *R, *tout, *tm, dmax = 0.0;
  double WMIN = map->WAVE[0], WMAX = map->WAVE[map->NWAVE - 1];
  char callFun[] = "do_profile", fnam_concat[200];
  FILE *fp;

  /* filter for apply_detcor_to_filter: widest given, else synthetic 200 pts */
  F = NULL;
  for (k = 0; k < nfilt; k++) if (!F || filt[k].NLAM > F->NLAM) F = &filt[k];
  if (!F) {
    syn.file = "synthetic"; snprintf(syn.base, sizeof(syn.base), "synthetic200");
    syn.NLAM = 200; syn.lam = malloc(200 * sizeof(double)); syn.trans = malloc(200 * sizeof(double));
    for (i = 0; i < 200; i++) { syn.lam[i] = 9000.0 + 60.0 * i; syn.trans[i] = 1.0; }
    F = &syn;
  }
  R    = malloc(sizeof(double) * map->NWAVE);
  tout = malloc(sizeof(double) * F->NLAM);
  printf("  profile filter: %s (NLAM=%d)\n", F->base, F->NLAM);

  /* ---- P1 ---- */
  t0 = now_ms();
  for (i = 0; i < 1000000; i++) sink += now_ms();
  dt = now_ms() - t0;
  printf("\n  P1 timer overhead: %.3f us per now_ms() call\n", dt / 1000.0);

  /* ---- P2 ---- */
  printf("\n  P2 scaling (random thickness 0-%g nm); cost/call should stay flat\n", tmax);
  printf("     %-20s %9s %12s %12s %10s\n", "function", "N", "total[ms]", "us/call", "maxRSS[MB]");
  for (N = 1000; N <= nmax; N *= 10) {
    t0 = now_ms();
    for (i = 0; i < N; i++) { get_detcor_ratio(map, urand(0, tmax), R); sink += R[i % map->NWAVE]; }
    dt = now_ms() - t0;
    printf("     %-20s %9d %12.2f %12.4f %10.1f\n", "get_detcor_ratio", N, dt, 1e3 * dt / N, maxrss_MB());
  }
  for (N = 100; N <= nmax / 10; N *= 10) {
    t0 = now_ms();
    for (i = 0; i < N; i++) { apply_detcor_to_filter(map, urand(0, tmax), F->NLAM, F->lam, F->trans, tout); sink += tout[i % F->NLAM]; }
    dt = now_ms() - t0;
    printf("     %-20s %9d %12.2f %12.4f %10.1f\n", "apply_detcor_to_filter", N, dt, 1e3 * dt / N, maxrss_MB());
  }

  /* ---- P3 ---- */
  nsweep = (int)(tmax + 60.0) + 1;               /* -30 .. tmax+30 in 1 nm */
  tm = malloc(sizeof(double) * M);
  {
    double *med_apply = malloc(sizeof(double) * nsweep), *med_ratio = malloc(sizeof(double) * nsweep);
    double *all = malloc(sizeof(double) * nsweep), gmed, worst = 0.0, tworst = 0.0;
    int j;
    fp = fopen("profile_thickness.dat", "w");
    if (fp) fprintf(fp, "# thick[nm] median_us_get_detcor_ratio median_us_apply_detcor_to_filter(%s) (%d calls each)\n", F->base, M);
    for (j = 0; j < nsweep; j++) {
      double t = -30.0 + j;
      for (i = 0; i < M; i++) { t0 = now_ms(); get_detcor_ratio(map, t, R); tm[i] = now_ms() - t0; sink += R[0]; }
      med_ratio[j] = 1e3 * median(tm, M);
      for (i = 0; i < M; i++) { t0 = now_ms(); apply_detcor_to_filter(map, t, F->NLAM, F->lam, F->trans, tout); tm[i] = now_ms() - t0; sink += tout[0]; }
      med_apply[j] = 1e3 * median(tm, M);
      if (fp) fprintf(fp, "%7.1f %10.4f %10.4f\n", t, med_ratio[j], med_apply[j]);
    }
    if (fp) fclose(fp);
    for (j = 0; j < nsweep; j++) all[j] = med_apply[j];
    gmed = median(all, nsweep);
    for (j = 0; j < nsweep; j++) {
      if (med_apply[j] > worst) { worst = med_apply[j]; tworst = -30.0 + j; }
      if (med_apply[j] > 1.5 * gmed) nbad++;
    }
    printf("\n  P3 thickness sweep: %d thicknesses (-30 to %g nm, 1 nm steps), %d calls each\n", nsweep, tmax + 30, M);
    printf("     apply_detcor_to_filter median over sweep: %.2f us; slowest %.2f us at t=%.0f nm\n", gmed, worst, tworst);
    printf("     thicknesses > 1.5x median: %d of %d   (full table: profile_thickness.dat)\n", nbad, nsweep);
    free(med_apply); free(med_ratio); free(all);
  }

  /* ---- P4 ---- */
  {
    int NP = 1000000;
    double *tv = malloc(NP * sizeof(double)), *wv = malloc(NP * sizeof(double));
    double d_full, d_qbs, d_i1d, d_cat, d_ref;
    for (i = 0; i < NP; i++) { tv[i] = urand(0, tmax); wv[i] = urand(WMIN, WMAX); }

    t0 = now_ms(); for (i = 0; i < NP; i++) sink += interp_detcor_ratio(map, tv[i], wv[i]); d_full = now_ms() - t0;
    t0 = now_ms(); for (i = 0; i < NP; i++) sink += quickBinSearch(tv[i], map->NTHICK, map->THICK, "prof", callFun); d_qbs = now_ms() - t0;
    t0 = now_ms(); for (i = 0; i < NP; i++) sink += interp_1DFUN(OPT_INTERP_LINEAR, wv[i], map->NWAVE, map->WAVE, map->RATIO[0], callFun); d_i1d = now_ms() - t0;
    t0 = now_ms(); for (i = 0; i < NP; i++) { concat_callfun_plus_fnam(callFun, "quickBinSearch", fnam_concat); sink += fnam_concat[i % 8]; } d_cat = now_ms() - t0;
    t0 = now_ms(); for (i = 0; i < NP; i++) sink += ref_ratio(map, tv[i], wv[i]); d_ref = now_ms() - t0;
    for (i = 0; i < NP; i++) { double d = fabs(ref_ratio(map, tv[i], wv[i]) - interp_detcor_ratio(map, tv[i], wv[i])); if (d > dmax) dmax = d; }

    printf("\n  P4 one lookup R(wave;thick), %d random points (ns per call)\n", NP);
    printf("     interp_detcor_ratio (SNANA tools)         : %8.1f ns\n", 1e6 * d_full / NP);
    printf("       = 1 x quickBinSearch on THICK[%d]    : %8.1f ns\n", map->NTHICK, 1e6 * d_qbs / NP);
    printf("       + 2 x interp_1DFUN on WAVE[%d]     : %8.1f ns  (%.1f each)\n", map->NWAVE, 2e6 * d_i1d / NP, 1e6 * d_i1d / NP);
    printf("         each quickBinSearch does sprintf : %8.1f ns  (concat_callfun_plus_fnam; 3 per lookup)\n", 1e6 * d_cat / NP);
    printf("     hand-coded reference (no SNANA calls)  : %8.1f ns  -> %.1fx faster; max|diff| = %.1e\n",
           1e6 * d_ref / NP, d_full / d_ref, dmax);
    free(tv); free(wv);
  }

  /* ---- P5 ---- */
  {
    int NRELOAD = 5;
    DETCOR_MAP_DEF m2;
    printf("\n  P5 reload map %d times (SNANA reader prints below)\n", NRELOAD);
    rss0 = maxrss_MB();
    for (k = 0; k < NRELOAD; k++) {
      t0 = now_ms(); read_detcor_map(map->FILENAME, &m2); dt = now_ms() - t0;
      tm[k] = dt; free_detcor_map(&m2);
    }
    printf("  P5 read_detcor_map times [ms]:");
    for (k = 0; k < NRELOAD; k++) printf(" %.1f", tm[k]);
    printf("   maxRSS %.1f -> %.1f MB\n", rss0, maxrss_MB());
  }

  printf("\n  (checksum %.6g)\n", sink);
  free(R); free(tout); free(tm);
  if (F == &syn) { free(syn.lam); free(syn.trans); }
  return 0;
}

int main(int argc, char **argv)
{
  DETCOR_MAP_DEF map;
  FILTER_DEF filt[MXFILT];
  char *mapfile = NULL;
  const char *mode = "timing", *dump = "detcor_map_test_dump.dat";
  const char *pycmd = getenv("HANDSHAKER_DETCOR") ? getenv("HANDSHAKER_DETCOR") : "handshaker-icemap";
  char *filtfile[MXFILT];
  int i, ntest = -1, seed = (int)time(NULL), rc, nfilt = 0;
  double tmax = -1.0, t0, dt;

  for (i = 1; i < argc; i++) {
    if      (!strcmp(argv[i], "--mode")   && i + 1 < argc) mode   = argv[++i];
    else if (!strcmp(argv[i], "--ntest")  && i + 1 < argc) ntest  = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--seed")   && i + 1 < argc) seed   = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--tmax")   && i + 1 < argc) tmax   = atof(argv[++i]);
    else if (!strcmp(argv[i], "--dump")   && i + 1 < argc) dump   = argv[++i];
    else if (!strcmp(argv[i], "--python") && i + 1 < argc) pycmd  = argv[++i];
    else if (!strcmp(argv[i], "--filter") && i + 1 < argc) {
      if (nfilt == MXFILT) { fprintf(stderr, "too many --filter (max %d)\n", MXFILT); return 2; }
      filtfile[nfilt++] = argv[++i];
    }
    else if (argv[i][0] != '-' && !mapfile) mapfile = argv[i];
    else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 2; }
  }
  if (!mapfile) {
    fprintf(stderr, "usage: %s <map.dat> [--mode timing|check] [--ntest N] [--seed S] "
            "[--tmax T] [--filter FILE]... [--dump OUT] [--python CMD]\n", argv[0]);
    return 2;
  }
  srand(seed);

  t0 = now_ms();
  read_detcor_map(mapfile, &map);       /* SNANA reader; aborts on error */
  dt = now_ms() - t0;
  if (tmax < 0) tmax = map.THICK[map.NTHICK - 1];
  printf("read_detcor_map: %s\n  NWAVE=%d (%.1f-%.1f A)  NTHICK=%d (%g-%g nm)  read in %.2f ms\n",
         mapfile, map.NWAVE, map.WAVE[0], map.WAVE[map.NWAVE - 1],
         map.NTHICK, map.THICK[0], map.THICK[map.NTHICK - 1], dt);
  printf("mode=%s seed=%d\n", mode, seed);

  if (!strcmp(mode, "check"))
    rc = do_check(&map, ntest > 0 ? ntest : 3, tmax, pycmd);
  else if (!strcmp(mode, "timing")) {
    for (i = 0; i < nfilt; i++) read_filter(filtfile[i], &filt[i]);
    rc = do_timing(&map, ntest > 0 ? ntest : 10000, tmax, filt, nfilt, dump);
    for (i = 0; i < nfilt; i++) { free(filt[i].lam); free(filt[i].trans); }
  }
  else if (!strcmp(mode, "profile")) {
    for (i = 0; i < nfilt; i++) read_filter(filtfile[i], &filt[i]);
    rc = do_profile(&map, ntest > 0 ? ntest : 1000000, tmax, filt, nfilt);
    for (i = 0; i < nfilt; i++) { free(filt[i].lam); free(filt[i].trans); }
  }
  else { fprintf(stderr, "unknown --mode %s (timing|check|profile)\n", mode); rc = 2; }

  free_detcor_map(&map);
  return rc;
}
