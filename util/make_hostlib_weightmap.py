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
#   MODEL1    [placeholder]  -- future WGTMAP approach #1
#   MODEL2    [placeholder]  -- future WGTMAP approach #2
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

  MODEL1     [placeholder -- not yet implemented]
    Future WGTMAP approach #1 (to be defined by developer).

  MODEL2     [placeholder -- not yet implemented]
    Future WGTMAP approach #2 (to be defined by developer).

Config file format (YAML):

    # Grid variables. Format: VARNAME(min:max:step)
    # step is optional; uses DVAR if omitted.
    VARLIST:  LOGMASS(7.0:13.0:0.1)  LOGSFR(-12.5:3.0:0.1)

    DVAR:        0.1         # default bin step when step not given in VARLIST
    SSFR_THRESH: -11.5       # log(sSFR/yr^-1) passive/SF boundary
    OUTPUT_DIR:  /path/to/output/dir  # optional; prepended to relative OUTFILE paths

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
# Placeholder modes  (not yet implemented)
# ══════════════════════════════════════════════════════════════════════════════

def build_wgtmap_model1(entry, config, df_hostlib, axes):
    """
    [PLACEHOLDER]  Future WGTMAP approach #1.

    Replace this docstring and body with the new implementation when ready.
    Suggested signature is unchanged: (entry, config, df_hostlib, axes) → (W_nd, notes).
    """
    label = entry.get("LABEL", "UNKNOWN")
    raise NotImplementedError(
        f"[{label}] MODE: MODEL1 is not yet implemented."
    )


def build_wgtmap_model2(entry, config, df_hostlib, axes):
    """
    [PLACEHOLDER]  Future WGTMAP approach #2.

    Replace this docstring and body with the new implementation when ready.
    Suggested signature is unchanged: (entry, config, df_hostlib, axes) → (W_nd, notes).
    """
    label = entry.get("LABEL", "UNKNOWN")
    raise NotImplementedError(
        f"[{label}] MODE: MODEL2 is not yet implemented."
    )


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
            W, notes = builder(entry, config, None, axes)
        except NotImplementedError as exc:
            sys.exit(f"ERROR: {exc}")

        write_wgtmap(outfile, W, axes, label=label, notes=notes, config=config)
        n_ok += 1

    print(f"\n{'='*64}")
    print(f"  Done. {n_ok}/{len(entries)} WGTMAP(s) written successfully.")
    print(f"{'='*64}\n")


if __name__ == "__main__":
    main()
