#!/usr/bin/env python3
#
# Created Mar 2026 by A. Mitra and CLAUDE CODE [AI]
#
# Config-driven utility to generate SNANA HOSTLIB weight maps (WGTMAP).
# Reads a YAML config file specifying the hostlib, grid variables,
# and one or more WGTMAP entries (each with a MODE, rate parameters,
# and output file path).
#
# Supported modes (MODE key in WGTMAP_ENTRIES):
#   RATEMOD   [implemented]  -- pure analytic rate 10^(beta*LOGMASS) on uniform grid
#   MODEL1    [implemented]  -- same rate formula, zeroed on empty HOSTLIB bins
#   MODEL2    [implemented]  -- beta fitted from observed SN/galaxy rate ratio
#
# Run with -H for full config format and examples.
#
# ================================================

import argparse
import gzip
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HELP_MSG = """
make_hostlib_weightmap.py  --  Config-driven SNANA HOSTLIB weight map generator

Usage:
    python make_hostlib_weightmap.py <config.yml>
    python make_hostlib_weightmap.py <config.yml> --outdir /path/to/output
    python make_hostlib_weightmap.py -H

Supported modes (MODE key in each WGTMAP_ENTRIES block):

  RATEMOD    [implemented]
    Pure analytic rate formula on a uniform grid — no HOSTLIB needed.
    W(M, S) = 10^(beta * LOGMASS) for cells passing GALAXY_SELECTION, else 0.
    P_host ∝ W × n_HOSTLIB = rate × n_HOSTLIB  (no double-counting).
    Reproduces the DES default WGTMAP logic.

  MODEL1     [implemented]
    Data-masked analytic rate. Same W = 10^(beta*LOGMASS) formula as RATEMOD,
    but cells with zero HOSTLIB galaxies are forced to W=0.
    Requires HOSTLIB_FILE in config.

  MODEL2     [implemented]
    Empirical rate-fitted beta. Beta is derived from the observed SN rate per
    galaxy as a function of host mass:
        rate(M) = N_SN(M) / N_gal(M)   [per LOGMASS bin]
    log10(rate(M)) is fit vs LOGMASS to extract beta_fit.
    N_SN comes from an observed FITRES file (HOST_LOGMASS column).
    N_gal comes from the HOSTLIB (after galaxy selection).
    W = 10^(beta_fit*LOGMASS) is then applied on the uniform grid.
    Requires HOSTLIB_FILE and FITRES_FILE in config. BETA is overridden.

Config file format (YAML):

    # Grid variables. Format: VARNAME(min:max:step)
    # step is optional; uses DVAR if omitted.
    VARLIST:  LOGMASS(7.0:13.0:0.1)  LOGSFR(-12.5:3.0:0.1)

    DVAR:         0.1          # default bin step when step not given in VARLIST
    SSFR_THRESH:  -11.5       # log(sSFR/yr^-1) passive/SF boundary
    OUTPUT_DIR:   /path/to/output/dir  # optional; prepended to relative OUTFILE paths
    HOSTLIB_FILE:        /path/to/hostlib.HOSTLIB  # required for MODEL1 / MODEL2
    FITRES_FILE:         /path/to/data.FITRES      # required for MODEL2
    FITRES_LOGMASS_COL:  HOST_LOGMASS              # optional; default HOST_LOGMASS

    WGTMAP_ENTRIES:
      - LABEL:             SNIA
        MODE:              RATEMOD
        BETA:              0.63
        GALAXY_SELECTION:  ALL          # ALL | PASSIVE | STAR_FORMING
        OUTFILE:           SNIA.WGTMAP.gz

      - LABEL:             SNII
        MODE:              RATEMOD
        BETA:              0.16
        GALAXY_SELECTION:  STAR_FORMING
        OUTFILE:           SNII.WGTMAP.gz

Key SNANA requirements enforced here:
  - WGTMAP grid must be COMPLETE and RECTANGULAR (all cells written, including zeros).
  - VARNAMES_WGTMAP: line must list all axis variables before WGT.
  - Output is gzip-compressed when OUTFILE ends in .gz.
"""

def expand_path(p):
    """Expand $ENV_VAR, ${ENV_VAR}, and ~ in a path string, return Path."""
    return Path(os.path.expanduser(os.path.expandvars(str(p))))


# ── global defaults ────────────────────────────────────────────────────────────
DEFAULT_DVAR        =   0.1
DEFAULT_SSFR_THRESH = -11.5

# Axis name aliases used to locate LOGMASS / LOGSFR in VARLIST
LOGMASS_CANDIDATES = ["LOGMASS", "LOGMASS_TRUE"]
LOGSFR_CANDIDATES  = ["LOGSFR", "LOG_SFR", "LOGSFR_TRUE"]


# ══════════════════════════════════════════════════════════════════════════════
# VARLIST parsing
# ══════════════════════════════════════════════════════════════════════════════

def parse_varlist(varlist_str, default_step=DEFAULT_DVAR):
    """
    Parse a VARLIST string into axis specifications.

    Supported token formats (space-separated):
        VARNAME(min:max:step)   -- explicit bin step
        VARNAME(min:max)        -- uses default_step

    Returns list of dicts:
        [{"name": str, "min": float, "max": float, "step": float,
          "edges": ndarray, "centers": ndarray}, ...]
    """
    axes   = []
    pat    = re.compile(r"(\w+)\(\s*([^:)]+)\s*:\s*([^:)]+)\s*(?::\s*([^)]+))?\s*\)")
    tokens = varlist_str.split()
    if not tokens:
        sys.exit("ERROR: VARLIST is empty.")

    for token in tokens:
        m = pat.match(token.strip())
        if not m:
            sys.exit(
                f"ERROR: cannot parse VARLIST token '{token}'.\n" \
                f"  Expected: VARNAME(min:max) or VARNAME(min:max:step)"
            )
        name  = m.group(1)
        vmin  = float(m.group(2))
        vmax  = float(m.group(3))
        step  = float(m.group(4)) if m.group(4) is not None else default_step

        # Build edges and centres; round centres to avoid float drift
        edges   = np.arange(vmin, vmax + step * 0.5001, step)
        centers = np.round(0.5 * (edges[:-1] + edges[1:]), 6)

        axes.append({
            "name": name, "min": vmin, "max": vmax,
            "step": step, "edges": edges, "centers": centers,
        })
        print(f"    {name}: [{vmin}, {vmax}], step={step}, {len(centers)} bins")

    return axes


# ══════════════════════════════════════════════════════════════════════════════
# WGTMAP writer
# ══════════════════════════════════════════════════════════════════════════════

def write_wgtmap(out_path, W_nd, axes, label, notes, config):
    """
    Write a SNANA-format WGTMAP to out_path (gzip if .gz extension).

    W_nd  : N-D weight array with shape (n_ax0, n_ax1, ..., n_axN-1).
            MUST be a complete rectangular grid (no sparse writing).
    axes  : list of axis specs from parse_varlist (same order as W_nd dims).
    label : string label for documentation block.
    notes : list of note strings written to DOCUMENTATION block.
    config: top-level config dict (used for provenance only).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_total = int(W_nd.size)
    n_nz    = int(np.count_nonzero(W_nd))
    wgt_sum = float(W_nd.sum())

    # Build VARNAMES_WGTMAP line: all axes then WGT then SNMAGSHIFT
    varnames_str = "  ".join(ax["name"] for ax in axes)

    doc_lines = [
        "DOCUMENTATION:",
        f"    LABEL:         {label}",
        f"    GENERATED_BY:  make_hostlib_weightmap.py",
        f"    NOTES:",
    ]
    for note in notes:
        doc_lines.append(f"    - {note}")
    doc_lines += [
        f"    GRID_TOTAL_BINS:  {n_total}",
        f"    GRID_NONZERO:     {n_nz}",
        f"    WGT_SUM:          {wgt_sum:.6f}",
        "DOCUMENTATION_END:",
        f"VARNAMES_WGTMAP:  {varnames_str}  WGT  SNMAGSHIFT",
    ]
    header = "\n".join(doc_lines) + "\n"

    # Flatten N-D grid: meshgrid gives coordinates matching W_nd.ravel() order
    grids   = np.meshgrid(*[ax["centers"] for ax in axes], indexing="ij")
    coords  = [g.ravel() for g in grids]   # one array per axis
    weights = W_nd.ravel()

    opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(out_path, "wt") as fh:
        fh.write(header)
        for i, w in enumerate(weights):
            coord_str = "  ".join(f"{coords[k][i]:9.4f}" for k in range(len(axes)))
            fh.write(f"WGT:  {coord_str}   {w:.6e}  0.000\n")

    size_kb = out_path.stat().st_size / 1024
    print(f"  Wrote: {out_path}  "
          f"({n_nz} non-zero / {n_total} total,  WGT sum={wgt_sum:.6f},  "
          f"{size_kb:.1f} KB)")


# ══════════════════════════════════════════════════════════════════════════════
# HOSTLIB reader  (shared by MODEL1 and MODEL2)
# ══════════════════════════════════════════════════════════════════════════════

def read_hostlib(hostlib_path):
    """
    Read a SNANA HOSTLIB file (plain or .gz) into a pandas DataFrame.

    Parses lines starting with GAL: using the VARNAMES: header.
    Returns a DataFrame with numeric columns where possible.
    """
    path = expand_path(hostlib_path)
    if not path.exists():
        sys.exit(f"ERROR: HOSTLIB_FILE not found: {path}")

    opener = gzip.open if str(path).endswith(".gz") else open
    varnames = None
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("VARNAMES:"):
                varnames = line.split()[1:]
            elif line.startswith("GAL:"):
                rows.append(line.split()[1:])

    if varnames is None:
        sys.exit(f"ERROR: no VARNAMES: line found in HOSTLIB {path}")
    if not rows:
        sys.exit(f"ERROR: no GAL: rows found in HOSTLIB {path}")

    df = pd.DataFrame(rows, columns=varnames)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(df[col])
    print(f"  Read {len(df):,} galaxies  ({len(varnames)} columns)  from {path.name}")
    return df


def read_fitres(fitres_path):
    """
    Read a SNANA FITRES file (plain or .gz) into a pandas DataFrame.

    Parses SN: rows using the VARNAMES: header, line by line.
    Skips all metadata/comment lines — robust to SNANA's variable-width headers.
    """
    path = expand_path(fitres_path)
    if not path.exists():
        sys.exit(f"ERROR: FITRES_FILE not found: {path}")

    opener = gzip.open if str(path).endswith(".gz") else open
    varnames = None
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("VARNAMES:"):
                varnames = line.split()[1:]
            elif line.startswith("SN:"):
                rows.append(line.split()[1:])

    if varnames is None:
        sys.exit(f"ERROR: no VARNAMES: line found in FITRES {path}")
    if not rows:
        sys.exit(f"ERROR: no SN: rows found in FITRES {path}")

    df = pd.DataFrame(rows, columns=varnames[:len(rows[0])])
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(df[col])

    print(f"  Read {len(df):,} SNe  ({len(varnames)} columns)  from {path.name}")
    return df


def _resolve_hostlib_col(df, ax_name, label):
    """
    Map a VARLIST axis name to the matching column in df.

    Checks the axis name directly, then falls back to LOGMASS_CANDIDATES /
    LOGSFR_CANDIDATES alias lists.  Exits with an error if not found.
    """
    if ax_name in df.columns:
        return ax_name
    for candidates in (LOGMASS_CANDIDATES, LOGSFR_CANDIDATES):
        if ax_name in candidates:
            for c in candidates:
                if c in df.columns:
                    return c
            sys.exit(
                f"ERROR [{label}]: axis '{ax_name}' not found in HOSTLIB columns "
                f"{list(df.columns)}"
            )
    sys.exit(
        f"ERROR [{label}]: VARLIST axis '{ax_name}' not found in HOSTLIB columns "
        f"{list(df.columns)}"
    )


def _apply_selection_mask_hostlib(df, selection, ssfr_thresh, label):
    """
    Filter a HOSTLIB DataFrame by GALAXY_SELECTION.

    Returns (df_selected, logm_col, logsfr_col_or_None).
    """
    logm_col = None
    for c in LOGMASS_CANDIDATES:
        if c in df.columns:
            logm_col = c
            break
    if logm_col is None:
        sys.exit(f"ERROR [{label}]: HOSTLIB has no LOGMASS column.")

    if selection == "ALL":
        return df, logm_col, None

    logsfr_col = None
    for c in LOGSFR_CANDIDATES:
        if c in df.columns:
            logsfr_col = c
            break
    if logsfr_col is None:
        sys.exit(
            f"ERROR [{label}]: GALAXY_SELECTION={selection} requires LOGSFR "
            f"column in HOSTLIB."
        )

    ssfr = df[logsfr_col] - df[logm_col]
    if selection == "PASSIVE":
        return df[ssfr < ssfr_thresh].copy(), logm_col, logsfr_col
    else:  # STAR_FORMING
        return df[ssfr >= ssfr_thresh].copy(), logm_col, logsfr_col


# ══════════════════════════════════════════════════════════════════════════════
# MODE: MODEL1  — analytic rate masked to HOSTLIB-populated bins
# ══════════════════════════════════════════════════════════════════════════════

def build_wgtmap_model1(entry, config, df_hostlib, axes):
    """
    Build a WGTMAP using the same rate formula as RATEMOD but restricted to
    grid cells that actually contain HOSTLIB galaxies.

    W(M, S) = 10^(beta * LOGMASS)   for cells where:
                                        n_HOSTLIB(M,S) > 0
                                        AND galaxy selection passes
    W(M, S) = 0                      otherwise

    Compared to RATEMOD, this zeros out physically empty (M, S) regions of
    parameter space, so SNANA will never draw a host from a combination that
    does not exist in the actual galaxy catalog.

    Required config key:
        HOSTLIB_FILE    str  -- path to SNANA HOSTLIB file (plain or .gz)

    Required entry keys (same as RATEMOD):
        BETA             float  -- mass-rate exponent
        GALAXY_SELECTION str   -- ALL | PASSIVE | STAR_FORMING
    """
    label       = entry.get("LABEL", "UNKNOWN")
    beta        = float(entry.get("BETA", 0.63))
    selection   = str(entry.get("GALAXY_SELECTION", "ALL")).upper()
    ssfr_thresh = float(config.get("SSFR_THRESH", DEFAULT_SSFR_THRESH))

    if df_hostlib is None:
        sys.exit(f"ERROR [{label}] MODEL1: HOSTLIB_FILE must be set in config.")

    ax_names = [ax["name"] for ax in axes]

    # ── locate axes ────────────────────────────────────────────────────────────
    logm_ax_idx = next(
        (i for i, n in enumerate(ax_names) if n in LOGMASS_CANDIDATES), None
    )
    if logm_ax_idx is None:
        sys.exit(f"ERROR [{label}] MODEL1: VARLIST must contain LOGMASS (or alias).")

    if selection != "ALL":
        logsfr_ax_idx = next(
            (i for i, n in enumerate(ax_names) if n in LOGSFR_CANDIDATES), None
        )
        if logsfr_ax_idx is None:
            sys.exit(
                f"ERROR [{label}] MODEL1: GALAXY_SELECTION={selection} requires "
                f"LOGSFR in VARLIST."
            )
    else:
        logsfr_ax_idx = None

    # ── bin HOSTLIB galaxies onto the VARLIST grid ─────────────────────────────
    hostlib_vals = [
        df_hostlib[_resolve_hostlib_col(df_hostlib, ax["name"], label)].values
        for ax in axes
    ]
    bins    = [ax["edges"] for ax in axes]
    n_data, _ = np.histogramdd(np.column_stack(hostlib_vals), bins=bins)
    # n_data shape: (n_ax0, n_ax1, ...) — count of HOSTLIB galaxies per cell

    # ── build meshgrid for rate formula ───────────────────────────────────────
    grids   = np.meshgrid(*[ax["centers"] for ax in axes], indexing="ij")
    logm_nd = grids[logm_ax_idx]

    # ── galaxy selection mask on grid ─────────────────────────────────────────
    if selection == "ALL":
        sel_mask = np.ones(logm_nd.shape, dtype=bool)
    else:
        logsfr_nd = grids[logsfr_ax_idx]
        ssfr_nd   = logsfr_nd - logm_nd
        sel_mask  = (
            (ssfr_nd < ssfr_thresh) if selection == "PASSIVE"
            else (ssfr_nd >= ssfr_thresh)
        )

    data_mask = n_data > 0
    mask      = sel_mask & data_mask

    n_total     = int(mask.size)
    n_sel       = int(sel_mask.sum())
    n_populated = int(data_mask.sum())
    n_used      = int(mask.sum())
    n_zeroed    = n_sel - n_used

    print(f"  GALAXY_SELECTION={selection}: {n_sel:,} / {n_total:,} grid cells pass cut")
    print(f"  HOSTLIB-populated: {n_populated:,} / {n_total:,} cells have galaxies")
    print(f"  Combined (selection AND data): {n_used:,} / {n_total:,} cells active")
    print(f"  Zeroed vs RATEMOD: {n_zeroed:,} cells set to 0 (no HOSTLIB galaxies)")

    # ── rate on grid ──────────────────────────────────────────────────────────
    W = np.where(mask, np.power(10.0, beta * logm_nd), 0.0)
    total = W.sum()
    if total > 0:
        W /= total
    else:
        print(f"  WARNING [{label}]: all-zero WGTMAP — check SSFR threshold / HOSTLIB coverage.")

    notes = [
        "MODE: MODEL1  (W = 10^(beta*LOGMASS), restricted to HOSTLIB-populated bins)",
        f"BETA = {beta}",
        f"GALAXY_SELECTION = {selection}",
        f"SSFR_THRESH = {ssfr_thresh}",
        f"HOSTLIB galaxies = {len(df_hostlib):,}",
        f"Grid cells active = {n_used:,}  of {n_total:,}",
        f"Grid cells zeroed (empty HOSTLIB bins) = {n_zeroed:,}  of {n_sel:,} selected",
    ]
    return W, notes


# ══════════════════════════════════════════════════════════════════════════════
# MODE: MODEL2  — beta fitted from observed SN rate per galaxy (N_SN / N_gal)
# ══════════════════════════════════════════════════════════════════════════════

def build_wgtmap_model2(entry, config, df_hostlib, df_fitres, axes):
    """
    Build a WGTMAP with beta derived from the empirical SN rate per galaxy:

        rate(M) = N_SN(M) / N_gal(M)   [per LOGMASS bin]

    N_SN(M) comes from an observed FITRES file (HOST_LOGMASS column).
    N_gal(M) comes from the HOSTLIB (after applying galaxy selection).

    Fit:  log10(rate(M)) = beta_fit * M + const
          weighted least squares, weights = sqrt(N_SN) for Poisson noise,
          over bins where both N_SN > 0 and N_gal > 0.

    W = 10^(beta_fit * LOGMASS) is then applied on the full uniform grid
    (same coverage as RATEMOD).

    Required config keys:
        HOSTLIB_FILE          str  -- path to SNANA HOSTLIB (plain or .gz)
        FITRES_FILE           str  -- path to observed SN FITRES (plain or .gz)
        FITRES_LOGMASS_COL    str  -- FITRES column for host logmass
                                      (optional; default HOST_LOGMASS)

    Required entry key:
        GALAXY_SELECTION str  -- ALL | PASSIVE | STAR_FORMING
        (BETA key is ignored; overridden by the fit)
    """
    label       = entry.get("LABEL", "UNKNOWN")
    selection   = str(entry.get("GALAXY_SELECTION", "ALL")).upper()
    ssfr_thresh = float(config.get("SSFR_THRESH", DEFAULT_SSFR_THRESH))
    logm_sn_col = str(config.get("FITRES_LOGMASS_COL", "HOST_LOGMASS"))

    if df_hostlib is None:
        sys.exit(f"ERROR [{label}] MODEL2: HOSTLIB_FILE must be set in config.")
    if df_fitres is None:
        sys.exit(f"ERROR [{label}] MODEL2: FITRES_FILE must be set in config.")
    if logm_sn_col not in df_fitres.columns:
        sys.exit(
            f"ERROR [{label}] MODEL2: column '{logm_sn_col}' not found in FITRES. "
            f"Available: {list(df_fitres.columns)}"
        )

    ax_names = [ax["name"] for ax in axes]

    # ── locate LOGMASS axis (and LOGSFR if needed) ────────────────────────────
    logm_ax_idx = next(
        (i for i, n in enumerate(ax_names) if n in LOGMASS_CANDIDATES), None
    )
    if logm_ax_idx is None:
        sys.exit(f"ERROR [{label}] MODEL2: VARLIST must contain LOGMASS (or alias).")

    if selection != "ALL":
        logsfr_ax_idx = next(
            (i for i, n in enumerate(ax_names) if n in LOGSFR_CANDIDATES), None
        )
        if logsfr_ax_idx is None:
            sys.exit(
                f"ERROR [{label}] MODEL2: GALAXY_SELECTION={selection} requires "
                f"LOGSFR in VARLIST."
            )
    else:
        logsfr_ax_idx = None

    logm_ax_spec = axes[logm_ax_idx]
    logm_edges   = logm_ax_spec["edges"]
    logm_centers = logm_ax_spec["centers"]

    # ── N_gal(M): HOSTLIB galaxies per LOGMASS bin (after galaxy selection) ───
    df_gal, logm_col, _ = _apply_selection_mask_hostlib(
        df_hostlib, selection, ssfr_thresh, label
    )
    print(
        f"  GALAXY_SELECTION={selection}: "
        f"{len(df_gal):,} / {len(df_hostlib):,} HOSTLIB galaxies pass cut"
    )
    if len(df_gal) == 0:
        sys.exit(f"ERROR [{label}] MODEL2: no HOSTLIB galaxies pass galaxy selection.")

    n_gal, _ = np.histogram(df_gal[logm_col].values, bins=logm_edges)

    # ── N_SN(M): observed SNe per LOGMASS bin (after dropping sentinels) ──────
    logm_sn = pd.to_numeric(df_fitres[logm_sn_col], errors="coerce")
    logm_sn = logm_sn[logm_sn > -9.0].dropna()   # drop -99 sentinels / NaN
    n_sn_used = len(logm_sn)
    print(f"  SNe with valid {logm_sn_col}: {n_sn_used:,}  of {len(df_fitres):,}")
    if n_sn_used == 0:
        sys.exit(f"ERROR [{label}] MODEL2: no valid SN host masses in FITRES.")

    n_sn, _ = np.histogram(logm_sn.values, bins=logm_edges)

    # ── rate(M) = N_SN / N_gal per bin, fit log10(rate) ~ beta * M ───────────
    fit_mask = (n_sn > 0) & (n_gal > 0)
    n_fit_bins = int(fit_mask.sum())
    print(f"  Bins with N_SN>0 AND N_gal>0: {n_fit_bins}  of {len(logm_centers)}")
    if n_fit_bins < 2:
        sys.exit(
            f"ERROR [{label}] MODEL2: fewer than 2 usable LOGMASS bins ({n_fit_bins}) "
            f"— cannot fit beta. Check FITRES host mass range vs VARLIST grid."
        )

    M_fit    = logm_centers[fit_mask]
    rate_fit = n_sn[fit_mask].astype(float) / n_gal[fit_mask].astype(float)
    logr_fit = np.log10(rate_fit)
    wts      = np.sqrt(n_sn[fit_mask].astype(float))   # Poisson weights on N_SN

    A = np.column_stack([M_fit, np.ones_like(M_fit)])
    coeffs, _, _, _ = np.linalg.lstsq(A * wts[:, None], logr_fit * wts, rcond=None)
    beta_fit, intercept_fit = coeffs

    print(f"  Fitted beta = {beta_fit:.4f}  (intercept = {intercept_fit:.4f})")
    print(f"  [reference beta from YAML: {entry.get('BETA', 'N/A')}]")

    # ── apply fitted beta on uniform grid (same footprint as RATEMOD) ─────────
    grids   = np.meshgrid(*[ax["centers"] for ax in axes], indexing="ij")
    logm_nd = grids[logm_ax_idx]

    if selection == "ALL":
        sel_mask = np.ones(logm_nd.shape, dtype=bool)
    else:
        logsfr_nd = grids[logsfr_ax_idx]
        ssfr_nd   = logsfr_nd - logm_nd
        sel_mask  = (
            (ssfr_nd < ssfr_thresh) if selection == "PASSIVE"
            else (ssfr_nd >= ssfr_thresh)
        )

    W = np.where(sel_mask, np.power(10.0, beta_fit * logm_nd), 0.0)
    total = W.sum()
    if total > 0:
        W /= total
    else:
        print(f"  WARNING [{label}]: all-zero WGTMAP.")

    notes = [
        "MODE: MODEL2  (W = 10^(beta_fit*LOGMASS); beta from empirical N_SN/N_gal rate fit)",
        f"BETA_FIT = {beta_fit:.6f}",
        f"INTERCEPT_FIT = {intercept_fit:.6f}",
        f"BETA_REFERENCE = {entry.get('BETA', 'N/A')}  (from YAML, not used)",
        f"GALAXY_SELECTION = {selection}",
        f"SSFR_THRESH = {ssfr_thresh}",
        f"N_SN used in fit: {n_sn_used:,}  (FITRES column: {logm_sn_col})",
        f"N_gal used in fit: {len(df_gal):,}  (HOSTLIB, after selection)",
        f"LOGMASS bins used in fit: {n_fit_bins}",
    ]
    return W, notes


# ══════════════════════════════════════════════════════════════════════════════
# MODE: RATEMOD  — pure analytic rate on uniform grid (no HOSTLIB binning)
# ══════════════════════════════════════════════════════════════════════════════

def build_wgtmap_ratemod(entry, config, df_hostlib, axes):
    """
    Build a WGTMAP using a pure analytic rate formula on a uniform grid.

    Reproduces the logic of make_analytic_wgtmap.py (DES default WGTMAPs):
        W(M, S) = 10^(beta * LOGMASS)   for cells passing GALAXY_SELECTION
        W(M, S) = 0                      otherwise

    No HOSTLIB is used — the grid is evaluated analytically at every cell.
    This avoids double-counting when the same HOSTLIB drives the SNANA sim:
        P(host=galaxy_i) ∝ W(M_i,S_i) × n_HOSTLIB(M_i,S_i)
                         = rate(M_i)   × n_HOSTLIB   (correct)

    Required entry keys:
        BETA             float  -- mass-rate exponent (e.g. 0.63 for SNe Ia)
        GALAXY_SELECTION str   -- ALL | PASSIVE | STAR_FORMING

    Returns (W_nd, notes_list).
    """
    label       = entry.get("LABEL", "UNKNOWN")
    beta        = float(entry.get("BETA", 0.63))
    selection   = str(entry.get("GALAXY_SELECTION", "ALL")).upper()
    ssfr_thresh = float(config.get("SSFR_THRESH", DEFAULT_SSFR_THRESH))

    valid_selections = {"ALL", "PASSIVE", "STAR_FORMING"}
    if selection not in valid_selections:
        sys.exit(
            f"ERROR [{label}]: unknown GALAXY_SELECTION '{selection}'.\n"
            f"  Valid options: {sorted(valid_selections)}"
        )

    ax_names = [ax["name"] for ax in axes]

    # ── locate LOGMASS axis (required for rate) ───────────────────────────────
    logm_ax = next(
        (i for i, n in enumerate(ax_names) if n in LOGMASS_CANDIDATES), None
    )
    if logm_ax is None:
        sys.exit(f"ERROR [{label}] RATEMOD: VARLIST must contain LOGMASS (or alias).")

    # ── build N-D meshgrid of bin centres ─────────────────────────────────────
    grids    = np.meshgrid(*[ax["centers"] for ax in axes], indexing="ij")
    logm_nd  = grids[logm_ax]   # shape (n_ax0, n_ax1, ...)

    # ── galaxy selection mask via sSFR = LOGSFR - LOGMASS ─────────────────────
    if selection == "ALL":
        mask = np.ones(logm_nd.shape, dtype=bool)
    else:
        logsfr_ax = next(
            (i for i, n in enumerate(ax_names) if n in LOGSFR_CANDIDATES), None
        )
        if logsfr_ax is None:
            sys.exit(
                f"ERROR [{label}] RATEMOD: GALAXY_SELECTION={selection} requires "
                f"LOGSFR (or alias) in VARLIST."
            )
        ssfr_nd = grids[logsfr_ax] - logm_nd
        mask = (ssfr_nd < ssfr_thresh) if selection == "PASSIVE" else (ssfr_nd >= ssfr_thresh)

    n_pass  = int(mask.sum())
    n_total = int(mask.size)
    print(f"  GALAXY_SELECTION={selection}: "
          f"{n_pass:,} / {n_total:,} grid cells pass cut  "
          f"(sSFR threshold={ssfr_thresh})")

    # ── rate on grid ──────────────────────────────────────────────────────────
    W = np.where(mask, np.power(10.0, beta * logm_nd), 0.0)
    total = W.sum()
    if total > 0:
        W /= total
    else:
        print(f"  WARNING [{label}]: all-zero WGTMAP — check SSFR threshold.")

    notes = [
        "MODE: RATEMOD  (W = 10^(beta*LOGMASS) on uniform grid, no HOSTLIB binning)",
        f"BETA = {beta}",
        f"GALAXY_SELECTION = {selection}",
        f"SSFR_THRESH = {ssfr_thresh}  (log[SFR/M*], only used for PASSIVE/STAR_FORMING)",
        f"Grid cells used = {n_pass:,}  of {n_total:,}",
    ]
    return W, notes


# ── dispatch table: register new modes here ───────────────────────────────────
MODE_DISPATCH = {
    "RATEMOD": build_wgtmap_ratemod,
    "MODEL1":  build_wgtmap_model1,
    "MODEL2":  build_wgtmap_model2,
}


# ══════════════════════════════════════════════════════════════════════════════
# CLI entry point
# ══════════════════════════════════════════════════════════════════════════════

def get_args():
    ap = argparse.ArgumentParser(
        description="Config-driven SNANA HOSTLIB weight map generator.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "config_file", nargs="?", default=None,
        help="Path to YAML config file.",
    )
    ap.add_argument(
        "-H", "--HELP", action="store_true",
        help="Print extended help with config file format and examples, then exit.",
    )
    ap.add_argument(
        "--outdir", default=None,
        help="Output directory. Prepended to any relative OUTFILE path in the config. "
             "Overrides OUTPUT_DIR if also set in the config. "
             "Absolute OUTFILE paths are always used as-is.",
    )
    return ap.parse_args()


def main():
    args = get_args()

    if args.HELP:
        print(HELP_MSG)
        sys.exit(0)

    if args.config_file is None:
        sys.exit("ERROR: missing config_file argument.\n"
                 "  Usage: make_hostlib_weightmap.py <config.yml>\n"
                 "  Run with -H for full help.")

    cfg_path = Path(args.config_file)
    if not cfg_path.exists():
        sys.exit(f"ERROR: config file not found: {cfg_path}")

    # ── resolve output directory (CLI > config > none) ──────────────────────
    # Priority: --outdir CLI arg > OUTPUT_DIR in config > use OUTFILE as-is
    # Only applied to relative OUTFILE paths; absolute paths are never modified.
    outdir_override = None
    if args.outdir:
        outdir_override = expand_path(args.outdir)
    # (config OUTPUT_DIR is resolved later, after config is loaded)

    print(f"\n{'='*64}")
    print(f"  make_hostlib_weightmap.py")
    print(f"  config: {cfg_path}")
    print(f"{'='*64}\n")

    with open(cfg_path) as fh:
        config = yaml.safe_load(fh)
    #sys.exit(f"XXX config = {config}")    

    # ── resolve output directory ─────────────────────────────────────────────
    # CLI --outdir beats config OUTPUT_DIR; both only affect relative OUTFILE paths.
    if outdir_override is None and config.get("OUTPUT_DIR"):
        outdir_override = expand_path(str(config["OUTPUT_DIR"]).strip())
    if outdir_override is not None:
        print(f"Output directory: {outdir_override}")

    # ── parse global VARLIST ─────────────────────────────────────────────────
    varlist_str = str(config.get("VARLIST", "")).strip()
    if not varlist_str:
        sys.exit("ERROR: VARLIST not specified in config.")
    default_dvar = float(config.get("DVAR", DEFAULT_DVAR))

    print("Parsing VARLIST ...")
    axes = parse_varlist(varlist_str, default_step=default_dvar)

    # ── validate entries ─────────────────────────────────────────────────────
    entries = config.get("WGTMAP_ENTRIES", [])
    if not entries:
        sys.exit("ERROR: WGTMAP_ENTRIES is empty or missing in config.")

    known_modes = list(MODE_DISPATCH.keys())
    for i, entry in enumerate(entries):
        label = entry.get("LABEL", f"entry_{i}")
        mode  = str(entry.get("MODE", "")).upper()
        print(f"Processing label {label} ")
        if not mode:
            sys.exit(f"ERROR [{label}]: MODE not specified.")
        if mode not in MODE_DISPATCH:
            sys.exit(
                f"ERROR [{label}]: unknown MODE '{mode}'.\n"
                f"  Known modes: {known_modes}"
            )
        if not entry.get("OUTFILE"):
            sys.exit(f"ERROR [{label}]: OUTFILE not specified.")

    # ── load HOSTLIB once if any entry needs it ───────────────────────────────
    DATA_MODES   = {"MODEL1", "MODEL2"}
    MODEL2_MODES = {"MODEL2"}
    entry_modes  = {str(e.get("MODE", "")).upper() for e in entries}

    df_hostlib = None
    if entry_modes & DATA_MODES:
        hostlib_path = config.get("HOSTLIB_FILE")
        if not hostlib_path:
            sys.exit(
                "ERROR: one or more entries use MODEL1/MODEL2 but HOSTLIB_FILE "
                "is not set in the config."
            )
        print(f"\nLoading HOSTLIB: {hostlib_path}")
        df_hostlib = read_hostlib(hostlib_path)

    # ── load FITRES once if any entry needs MODEL2 ────────────────────────────
    df_fitres = None
    if entry_modes & MODEL2_MODES:
        fitres_path = config.get("FITRES_FILE")
        if not fitres_path:
            sys.exit(
                "ERROR: one or more entries use MODEL2 but FITRES_FILE "
                "is not set in the config."
            )
        print(f"\nLoading FITRES: {fitres_path}")
        df_fitres = read_fitres(fitres_path)

    # ── process each WGTMAP entry ────────────────────────────────────────────
    n_ok = 0
    for entry in entries:
        label   = entry.get("LABEL", "UNKNOWN")
        mode    = str(entry.get("MODE", "")).upper()
        outfile = Path(entry.get("OUTFILE"))
        if outdir_override is not None and not outfile.is_absolute():
            outfile = outdir_override / outfile

        print(f"\n{'─'*64}")
        print(f"  LABEL={label}   MODE={mode}")
        print(f"{'─'*64}")

        builder = MODE_DISPATCH[mode]
        try:
            if mode == "MODEL2":
                W, notes = builder(entry, config, df_hostlib, df_fitres, axes)
            else:
                W, notes = builder(entry, config, df_hostlib, axes)
        except NotImplementedError as exc:
            sys.exit(f"ERROR: {exc}")

        write_wgtmap(outfile, W, axes, label=label, notes=notes, config=config)
        n_ok += 1

    print(f"\n{'='*64}")
    print(f"  Done. {n_ok}/{len(entries)} WGTMAP(s) written successfully.")
    print(f"{'='*64}\n")


if __name__ == "__main__":
    main()
