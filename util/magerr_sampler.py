"""
Efficient sampling of magnitude errors from the Diffsky magerr model.

Strategy
--------
The YAML gives, per band and per mag bin, a 2-component mixture:
    normal Gaussian  (amp, mu, sigma)
    skew Gaussian    (amp, a, loc, scale)

Instead of interpolating the 7 fit parameters between mag bins (fragile: the
amplitudes are not on a comparable normalization, and individual fits can be
pathological), we precompute the inverse CDF (quantile function) of each
fitted mixture on a fixed grid of probabilities.  Sampling is then

    u  ~ U(0,1)                                   one draw per galaxy
    t  = fractional mag-bin coordinate
    err = bilinear_interp(Q[band], t, u)

which is quantile ("Wasserstein") interpolation between adjacent mag bins:
a genuine distribution that correctly morphs both location and width, and is
immune to component swapping or bad individual fits.

Cost at runtime: one RNG draw + one map_coordinates call per band. No loops,
no rejection sampling, no per-object branching.
"""

import sys
import warnings

import numpy as np
import yaml
from scipy.ndimage import map_coordinates
from scipy.special import ndtr  # standard normal CDF, vectorized
from scipy.stats import norm, skewnorm

# How the fitted amplitudes should be read.  The fitter records this in the yaml
# as FIT_QUALITY.AMP_CONVENTION, and that is what we use -- reading a file with
# the wrong convention produces a wrong-but-plausible distribution, so it must
# never be guessed.
#   "pdf"  : model = A1 * skewnormpdf(x; a, loc, scale) + A2 * normpdf(x; mu, sigma)
#            -> both amplitudes are areas and are directly comparable
#   "bare" : model = A1 * skewnormpdf(...) + A2 * exp(-(x-mu)^2/(2 sigma^2))
#            -> A2 is a peak HEIGHT, area = A2*sigma*sqrt(2pi)
#
# Only used when the yaml does not say: files written before 2026-09-16 have no
# AMP_CONVENTION key and are all "bare".
AMP_CONVENTION_DEFAULT = "bare"

N_QUANTILES = 2048   # resolution of the stored inverse CDF
N_XGRID = 4096       # resolution of the grid used to build each CDF

# Fit quality is decided by the FITTER, not here: diffsky_magerr_model.classify_fit
# writes a FITFLAG / FITFLAG_REASON column per band into the yaml, along with a
# top-level FIT_QUALITY block recording the thresholds it used.  This sampler
# only reads those.  The fitter can check things we cannot -- in particular the
# four parameters whose upper bound is the per-bin FIT_MAX, and the dof -- so
# duplicating its criteria here would be both redundant and less capable.
#
# LEGACY_CHI2_MAX below is used *only* for the fallback path that handles yaml
# files written before FITFLAG existed.
LEGACY_CHI2_MAX = 20.0

FLAG_GOOD, FLAG_WARN, FLAG_BAD = "GOOD", "WARN", "BAD"


class MagErrSampler:
    def __init__(self, yaml_path, bands=("g", "r", "i", "z"),
                 field="NONE", amp_convention=None,
                 nq=N_QUANTILES, clip_negative=True, bin_label="center",
                 chi2_max=None, repair_bad_bins=False,
                 bad_severities=(FLAG_BAD,), verbose=True):
        with open(yaml_path) as f:
            model = yaml.safe_load(f)

        # diffsky_magerr_model.fit_magerr_func builds
        #     magbins = np.arange(lo, hi + 0.001, d)          -> 20, 20.5, ... 27.0
        # and fits one row per bin [magbins[i], magbins[i+1]] over magbins[:-1].
        # So np.arange(lo, hi, d) gives each row's LEFT EDGE (20, 20.5, ... 26.5),
        # and the magnitude representative of the galaxies inside a row is the bin
        # CENTER, half a bin higher (20.25, 20.75, ... 26.75).  Use bin_label="left"
        # to anchor on the edges instead.
        lo, hi, d = model["MAG_BINS"]
        if bin_label not in ("center", "left"):
            raise ValueError(f"bin_label must be 'center' or 'left', got {bin_label!r}")
        edges = np.arange(lo, hi, d)
        nbin = len(edges)
        self.bin_edges = edges
        self.mag_centers = edges + (0.5 * d if bin_label == "center" else 0.0)
        self.dmag = d
        self.bands = list(bands)
        self.nq = nq
        self.clip_negative = clip_negative
        # Resolve the amplitude convention: explicit argument wins, then what the
        # fitter recorded in the yaml, then the pre-2026-09-16 default.
        self.fit_quality = model.get("FIT_QUALITY", {})
        if amp_convention is None:
            amp_convention = self.fit_quality.get("AMP_CONVENTION",
                                                  AMP_CONVENTION_DEFAULT)
        if amp_convention not in ("bare", "pdf"):
            raise ValueError(f"amp_convention must be 'bare' or 'pdf', got {amp_convention!r}")
        self.amp_convention = amp_convention

        self.gauss_par = {}
        self.skew_par = {}
        self.median = {}
        self.data_mean = {}
        self.data_std = {}
        self.fitchi2 = {}
        self.fitflag = {}
        self.fitreason = {}
        self.fitdof = {}
        self.ndata = {}
        self.fit_max = {}
        for b in self.bands:
            blk = model[field][b]
            g = np.asarray(blk["ERR_FUNC"]["GAUSS_PAR"], float)
            s = np.asarray(blk["ERR_FUNC"]["SKEWGAUSS_PAR"], float)
            assert g.shape == (nbin, 3) and s.shape == (nbin, 4), (b, g.shape, s.shape)
            self.gauss_par[b] = g
            self.skew_par[b] = s
            self.median[b] = np.asarray(blk["MEDIAN"], float)
            # Summary stats of the magerr data per bin, used by
            # sample(check_bin_means=...).  Absent in older yaml.
            self.data_mean[b] = np.asarray(blk["MEAN"], float) \
                if "MEAN" in blk else None
            self.data_std[b] = np.asarray(blk["STDDEV"], float) \
                if "STDDEV" in blk else None
            # All of these are absent in yaml written before the fitter learned
            # to record fit quality.
            self.fitchi2[b] = np.asarray(blk["FITCHI2"], float) \
                if "FITCHI2" in blk else None
            self.fitflag[b] = np.asarray(blk["FITFLAG"]) \
                if "FITFLAG" in blk else None
            self.fitreason[b] = list(blk.get("FITFLAG_REASON", [""] * nbin))
            self.fitdof[b] = np.asarray(blk["DOF"], float) if "DOF" in blk else None
            self.ndata[b] = np.asarray(blk["NDATA"], float) if "NDATA" in blk else None
            self.fit_max[b] = np.asarray(blk["FIT_MAX"], float) \
                if "FIT_MAX" in blk else None

        # Cross-band correlation -> Cholesky factor for the Gaussian copula.
        self.corr = np.asarray(model.get("PEARSON_COEFS"), float) \
            if "PEARSON_COEFS" in model else None
        self.chol = np.linalg.cholesky(self.corr) if self.corr is not None else None

        # Fit quality: read the fitter's verdict, or fall back for old files.
        self.chi2_max = chi2_max
        self.bad_severities = tuple(bad_severities)
        if all(self.fitflag[b] is not None for b in self.bands):
            self.flag_source = "FITFLAG"
            self.fit_flags = {b: self._read_fit_flags(b) for b in self.bands}
        elif any(self.fitchi2[b] is not None for b in self.bands):
            self.flag_source = "legacy-chi2"
            warnings.warn(
                "yaml has no FITFLAG column, so it predates fit-quality flagging "
                "in diffsky_magerr_model.py. Falling back to a chi2-only check; "
                "regenerate the yaml for the full set of flags.")
            self.fit_flags = {b: self._legacy_flag_bins(b) for b in self.bands}
        else:
            self.flag_source = "none"
            warnings.warn("yaml has neither FITFLAG nor FITCHI2: no fit-quality "
                          "information at all, treating every fit as GOOD.")
            self.fit_flags = {b: {} for b in self.bands}

        # Q[band] has shape (nbin, nq): the quantile function of each mag bin.
        self.p_grid = (np.arange(nq) + 0.5) / nq
        self.Q = {b: self._build_quantile_table(b) for b in self.bands}

        # Rebuild untrustworthy rows from their nearest trustworthy neighbours.
        # xxx I am not sure if this repair_bad_bins feature works properly, so
        # keep it set to False
        self.repaired = {b: (self._repair_band(b) if repair_bad_bins else {})
                         for b in self.bands}

        # Handing out NaN magerrs is bad, but handing out silently plausible
        # ones would be worse -- so leave the NaNs and say so clearly.
        nan_rows = {b: sorted(np.where(~np.all(np.isfinite(self.Q[b]), axis=1))[0].tolist())
                    for b in self.bands}
        nan_rows = {b: v for b, v in nan_rows.items() if v}
        if nan_rows:
            warnings.warn(
                f"quantile table still has non-finite rows {nan_rows} "
                f"(the fitter could not fit those bins). Sampling at those "
                f"magnitudes will return NaN. "
                + ("Pass repair_bad_bins=True to rebuild them from neighbours."
                   if not repair_bad_bins else
                   "Every neighbouring bin is flagged too, so there was nothing "
                   "to rebuild them from."))

        if verbose:
            nflag = sum(len(v) for v in self.fit_flags.values())
            nfix = sum(len(v) for v in self.repaired.values())
            if nflag:
                nbad = sum(1 for v in self.fit_flags.values()
                           for sev, _ in v.values() if sev == FLAG_BAD)
                print(f"MagErrSampler: {nflag} of {nbin * len(self.bands)} fits "
                      f"flagged ({nbad} {FLAG_BAD}), {nfix} rebuilt from "
                      f"neighbours [flags: {self.flag_source}]. "
                      "Call .report() for details.")

    # ------------------------------------------------------------ fit quality

    def _read_fit_flags(self, band):
        """Read the fitter's FITFLAG / FITFLAG_REASON columns.

        Returns {bin_index: (severity, [reasons])} for every non-GOOD bin, the
        same shape the rest of this class already expects.
        """
        flag = self.fitflag[band]
        reason = self.fitreason[band]
        chi2 = self.fitchi2[band]
        nbin = len(self.mag_centers)
        assert flag.shape == (nbin,), (band, flag.shape)

        # Forward compatibility: a newer fitter may invent severities we do not
        # know.  Treat anything unrecognised as WARN rather than crashing, and
        # never silently drop it.
        known = (FLAG_GOOD, FLAG_WARN, FLAG_BAD)
        unknown = sorted(set(np.asarray(flag).tolist()) - set(known))
        if unknown:
            warnings.warn(f"band {band}: unrecognised FITFLAG value(s) {unknown}; "
                          f"treating as {FLAG_WARN}")

        flags = {}
        for k in range(nbin):
            sev = str(flag[k])
            if sev == FLAG_GOOD:
                continue
            if sev not in known:
                sev = FLAG_WARN
            why = [r for r in str(reason[k]).split("; ") if r]
            if not why:
                why = ["flagged by the fitter, no reason recorded"]
            flags[k] = (sev, why)

        # An explicit chi2_max overrides the file, but only to be STRICTER --
        # it can promote a bin to BAD, never acquit one the fitter condemned.
        if self.chi2_max is not None and chi2 is not None:
            for k in range(nbin):
                if np.isfinite(chi2[k]) and chi2[k] > self.chi2_max:
                    sev, why = flags.get(k, (FLAG_GOOD, []))
                    extra = (f"chi2/dof={chi2[k]:.1f} > {self.chi2_max:g} "
                             "(caller override)")
                    flags[k] = (FLAG_BAD, why + [extra])
        return flags

    def _legacy_flag_bins(self, band):
        """chi2-only fallback for yaml written before FITFLAG existed.

        Deliberately does NOT try to reproduce the fitter's bound-pin checks.
        Those only ever produced WARN, which is reported but never acted on, so
        dropping them changes no sampling result -- and it lets the duplicated
        copies of curve_fit's bounds disappear from this file for good.
        """
        chi2 = self.fitchi2[band]
        if chi2 is None:
            return {}
        chi2_max = self.chi2_max if self.chi2_max is not None else LEGACY_CHI2_MAX
        flags = {}
        for k in range(len(self.mag_centers)):
            if not np.isfinite(chi2[k]) or chi2[k] <= 0:
                flags[k] = (FLAG_BAD, [f"chi2/dof={chi2[k]:.3g} is <= 0 or not "
                                       "finite (degenerate fit)"])
            elif chi2[k] > chi2_max:
                flags[k] = (FLAG_BAD, [f"chi2/dof={chi2[k]:.1f} > {chi2_max:g}"])
        return flags

    def _repair_band(self, band):
        """Replace each "bad" row of Q with a blend of its nearest good rows.

        This is the same quantile interpolation sample() uses between bins, just
        applied to overwrite a row rather than to fill the gap between two.
        Returns {bin_index: (k_lo, k_hi, weight)} describing what was done.
        """
        bad = [k for k, (sev, _) in self.fit_flags[band].items()
               if sev in self.bad_severities]
        # A row whose parameters are non-finite (the fitter could not fit it at
        # all) must be rebuilt even if it somehow escaped flagging, or it would
        # hand out NaN magerrs.
        for k in range(self.Q[band].shape[0]):
            if k not in bad and not np.all(np.isfinite(self.Q[band][k])):
                bad.append(k)
        bad = sorted(bad)
        if not bad:
            return {}
        Q = self.Q[band]
        good = [k for k in range(Q.shape[0]) if k not in bad]
        if not good:
            warnings.warn(f"band {band}: every fit is flagged, leaving all rows as fit")
            return {}

        Q_orig = Q.copy()          # so repairs never chain off other repairs
        done = {}
        for k in bad:
            below = [g for g in good if g < k]
            above = [g for g in good if g > k]
            if below and above:
                k_lo, k_hi = below[-1], above[0]
                w = (k - k_lo) / (k_hi - k_lo)
            else:                   # flagged bin at an end: copy nearest good row
                k_lo = k_hi = below[-1] if below else above[0]
                w = 0.0
            Q[k] = (1.0 - w) * Q_orig[k_lo] + w * Q_orig[k_hi]
            done[k] = (k_lo, k_hi, w)
        return done

    def report(self, stream=None):
        """Print a per-band fit-quality table and what was rebuilt."""
        out = stream if stream is not None else sys.stdout
        def w(line=""):
            out.write(line + "\n")

        lo, d = self.bin_edges[0], self.dmag

        # Provenance first.  Without this you cannot tell a freshly flagged yaml
        # from a stale one, or know which rule produced a given verdict.
        if self.flag_source == "FITFLAG":
            w("Fit-quality flags: read from the yaml (written by "
              "diffsky_magerr_model.classify_fit)")
            fq = self.fit_quality
            if fq:
                shown = ", ".join(f"{k}={fq[k]}" for k in
                                  ("CHI2_MAX", "NDATA_MIN", "NDATA_WARN",
                                   "DOF_MIN", "PIN_RTOL", "VESTIGIAL_FRAC")
                                  if k in fq)
                w(f"  fitter thresholds: {shown}")
        elif self.flag_source == "legacy-chi2":
            cm = self.chi2_max if self.chi2_max is not None else LEGACY_CHI2_MAX
            w(f"Fit-quality flags: RECOMPUTED here, chi2-only (chi2_max={cm:g}) "
              "-- this yaml predates FITFLAG")
        else:
            w("Fit-quality flags: NONE AVAILABLE in this yaml")
        if self.chi2_max is not None and self.flag_source == "FITFLAG":
            w(f"  caller override: chi2_max={self.chi2_max:g} (stricter only)")
        w(f"  bins treated as replaceable: {', '.join(self.bad_severities)}")
        w(f"  amplitude convention: {self.amp_convention}"
          + ("" if "AMP_CONVENTION" in self.fit_quality
             else "  (not recorded in the yaml; assumed)"))
        w()

        w(f"{'mag bin':>13} " + " ".join(f"{b:>9}" for b in self.bands))
        w("-" * (13 + 10 * len(self.bands)))
        for k in range(len(self.mag_centers)):
            cells = []
            for b in self.bands:
                c = self.fitchi2[b]
                txt = "    n/a" if c is None else f"{c[k]:7.2f}"
                mark = "*" if k in self.repaired[b] else \
                       ("!" if k in self.fit_flags[b] else " ")
                cells.append(f"{txt}{mark} ")
            w(f"{lo + d * k:6.1f}-{lo + d * (k + 1):5.1f} " + " ".join(cells))
        w()
        w("  (* rebuilt from neighbours, ! flagged but kept)")

        for b in self.bands:
            if not self.fit_flags[b]:
                continue
            w()
            w(f"{b} band:")
            for k, (sev, why) in sorted(self.fit_flags[b].items()):
                fix = self.repaired[b].get(k)
                if fix is None:
                    action = "kept as fit"
                elif fix[0] == fix[1]:
                    action = f"copied bin {fix[0]}"
                else:
                    # Spanning more than one bin is a much weaker statement than
                    # spanning one, so say how wide the gap actually was.
                    span = (fix[1] - fix[0]) * d
                    action = (f"rebuilt from bins {fix[0]} and {fix[1]} "
                              f"(w={fix[2]:.2f}, spanning {span:.1f} mag)")
                extra = ""
                if self.ndata[b] is not None and np.isfinite(self.ndata[b][k]):
                    extra += f", ndata={int(self.ndata[b][k])}"
                if self.fitdof[b] is not None and np.isfinite(self.fitdof[b][k]):
                    extra += f", dof={int(self.fitdof[b][k])}"
                w(f"  bin {k:2d} ({lo + d * k:.1f}-{lo + d * (k + 1):.1f}) "
                  f"{sev}{extra}: {action}")
                for r in why:
                    w(f"           - {r}")

    # ---------------------------------------------------------------- tables

    def _mixture_weights(self, amp_g, sigma, amp_s):
        """Integrated (unnormalized) mass of each component."""
        if self.amp_convention == "bare":
            m_g = amp_g * sigma * np.sqrt(2.0 * np.pi)
        else:
            m_g = amp_g
        return m_g, amp_s

    def _build_quantile_table(self, band):
        gpar, spar = self.gauss_par[band], self.skew_par[band]
        nbin = gpar.shape[0]
        Q = np.empty((nbin, self.nq))

        for k in range(nbin):
            # The fitter writes NaN parameters when curve_fit failed outright.
            # np.linspace(nan, nan) would silently poison the whole row, so mark
            # it NaN here and let _repair_band rebuild it from its neighbours.
            if not (np.all(np.isfinite(gpar[k])) and np.all(np.isfinite(spar[k]))):
                Q[k] = np.nan
                continue

            amp_g, mu, sig = gpar[k]
            amp_s, a, loc, scale = spar[k]
            sig, scale = abs(sig), abs(scale)

            m_g, m_s = self._mixture_weights(amp_g, sig, amp_s)
            tot = m_g + m_s
            w_g, w_s = m_g / tot, m_s / tot

            # Grid wide enough to hold both components in either tail.
            xlo = min(mu - 8 * sig, loc - 8 * scale)
            xhi = max(mu + 8 * sig, loc + 8 * scale + abs(a) * scale)
            if self.clip_negative:
                xlo = max(xlo, 0.0)
            x = np.linspace(xlo, xhi, N_XGRID)

            pdf = w_g * norm.pdf(x, mu, sig) + w_s * skewnorm.pdf(x, a, loc, scale)

            # Trapezoid CDF, normalized over the (possibly truncated) support.
            cdf = np.concatenate(([0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(x))))
            cdf /= cdf[-1]

            # Invert.  Keep only strictly increasing samples so np.interp is well posed.
            keep = np.concatenate(([True], np.diff(cdf) > 0))
            Q[k] = np.interp(self.p_grid, cdf[keep], x[keep])

        return Q

    # --------------------------------------------------------------- sampling

    def _mag_coord(self, mag):
        """Fractional bin coordinate, clipped (flat extrapolation off the ends)."""
        t = (np.asarray(mag, float) - self.mag_centers[0]) / self.dmag
        return np.clip(t, 0.0, len(self.mag_centers) - 1.0)

    def _check_bin_means(self, band, mag, err, nsig):
        """Raise if, in any mag bin, |mean(err) - MEAN| / STDDEV > nsig.

        MEAN and STDDEV are the yaml's summary stats of the magerr data in that
        bin.  The offset is measured in units of STDDEV (the scatter of the
        data), NOT the standard error STDDEV/sqrt(n): the fitted model is
        slightly off the data mean in most bins, so a standard-error test would
        fail on any large enough sample even when the sampler works correctly.
        Mags outside the binned range, and bins with no usable MEAN/STDDEV, are
        skipped.
        """
        mu, sd = self.data_mean[band], self.data_std[band]
        if mu is None or sd is None:
            raise ValueError(f"check_bin_means needs MEAN and STDDEV for band "
                             f"{band}, but the yaml does not have them")
        nbin = len(self.bin_edges)
        k = np.floor((np.ravel(mag) - self.bin_edges[0]) / self.dmag).astype(int)
        err = np.ravel(err)
        inside = (k >= 0) & (k < nbin) & np.isfinite(err)
        n = np.bincount(k[inside], minlength=nbin)
        s = np.bincount(k[inside], weights=err[inside], minlength=nbin)
        for j in np.where(n > 0)[0]:
            if not (np.isfinite(mu[j]) and np.isfinite(sd[j]) and sd[j] > 0):
                continue
            z = (s[j] / n[j] - mu[j]) / sd[j]
            if abs(z) > nsig:
                lo = self.bin_edges[j]
                raise ValueError(
                    f"Mean of sampled errors is not consistent with error model "
                    f"mean in {band} band, mag bin {lo:.2f}-{lo + self.dmag:.2f}: "
                    f"sampled mean={s[j] / n[j]:.5g} (n={n[j]}), "
                    f"model MEAN={mu[j]:.5g}, STDDEV={sd[j]:.5g}, "
                    f"|diff|/STDDEV={abs(z):.3g} > {nsig:g}")

    def sample(self, band, mag, rng=None, u=None, check_bin_means=None):
        """Draw one magerr per entry of `mag` for a single band.
        Pass `u` (uniforms in (0,1)) to drive the draw yourself -- that is how
        cross-band correlation is injected.
        check_bin_means: if a number, after drawing, compare the mean sampled
        error in each mag bin with the yaml's MEAN for that bin and raise
        ValueError if |diff|/STDDEV exceeds it.  None (default) skips the check.
        t and u are the coordinates used to locate a magerr in the Q table
        t=8.7 means take a mag between the 8th and 9th bin center, weighted 70%
        towards the 9th
        u is the percentile, selected randomly, in the magerr distribution to draw
        the magerr
        """
        rng = rng or np.random.default_rng()
        # map_coordinates cannot take 0-d coordinates, so a bare scalar mag
        # becomes a length-1 array.
        mag = np.atleast_1d(mag)
        t = self._mag_coord(mag)
        if u is None:
            u = rng.random(t.shape)
        u = np.atleast_1d(u)
        # Bilinear: linear in mag between bins, linear in u within the quantile
        # table -> exactly quantile interpolation between adjacent mag bins.
        coords = np.stack([t, np.asarray(u) * (self.nq - 1)])
        err = map_coordinates(self.Q[band], coords, order=1, mode="nearest")
        if self.clip_negative:
            err = np.maximum(err, 0.0)
        if check_bin_means is not None:
            self._check_bin_means(band, mag, err, check_bin_means)
        return err
