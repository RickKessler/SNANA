#!/usr/bin/env python
#
# Created May 2026 by Jonah Medoff
# Translate tables (csv or fitres) into cigale-readable input tables;
# Translate cigale output into usable table for SNANA.
#
#

import argparse
import csv
import gzip
import numpy as np
import yaml
import sys
from astropy.io import fits


### Functions used by all modes ###

def get_args():
    parser = argparse.ArgumentParser(description="Wrapper around cigale to make it compatible with SNANA")
    parser.add_argument("config", help="Path to yml config file")
    parser.add_argument("--input_table_file", help="Input SNANA table file; e.g., HOSTLIB or FITRES")
    parser.add_argument("--input_cigale_results", help="Cigale results file (results.fits)")
    parser.add_argument("--output_cigale_file", help="Outputted cigale input file")
    parser.add_argument("--output_galid_map", help="Map between cigale ids and SNANA GALIDs")
    parser.add_argument("--output_snana_file", help="Outputted SNANA file with host properties included")
    parser.add_argument("--mode", required=True, help="snana_to_cigale or cigale_to_snana")
    parser.add_argument("--zgrid", nargs=3, help="Redshift grid")
    parser.add_argument("--varname_z", help="Redshift column in input file")
    args = parser.parse_args()
    return args


def parse_config(config_path, mode):
    """mode = 'SNANA_TO_CIGALE' or 'CIGALE_TO_SNANA' """
    with open(config_path) as f:
        raw = yaml.safe_load(f)
    return raw[mode]


def resolve(arg_value, config, config_key, cli_flag, mode):
    """Return arg_value if given, else config[config_key]; raise if neither."""
    if arg_value is not None:
        return arg_value
    if config_key in config:
        return config[config_key]
    raise ValueError(
        f"Missing required value for {mode}: pass {cli_flag} on the command line "
        f"or set {config_key} in config file."
    )


### SNANA_TO_CIGALE-specific functions ###

def mag_err(mags, depth):
    """
    Calculate mag_err if given a depth
    """
    SNR = 5.0 * 10 ** (0.4 * (depth - mags))
    return 2.5 / (np.log(10) * SNR)


def mag_to_mjy(mag_ab):
    """
    Convert AB magnitude to flux in milliJansky
    """
    return 3.631e6 * 10 ** (-0.4 * mag_ab)


def mag_err_to_mjy_err(flux_mjy, mag_err_val):
    """
    Propagate mag error to flux error in mJy: sigma_F = 0.4 * ln(10) * F * sigma_m
    Notice that if you expand this formula by plugging in the formulas for F and sigma_m,
    you find that sigma_F does not depend at all on the mags. This means that, for constant
    depths, I will get the same flux error on all galaxies
    """
    return 0.4 * np.log(10) * flux_mjy * mag_err_val


REQUIRED_COLUMN_KEYS = {"id"}
REQUIRED_MAG_KEYS = {"ctio.des.g", "ctio.des.r", "ctio.des.i", "ctio.des.z"}


def parse_column_map(column_map, varname_z = None):
    """
    Parse CIGALE_COLUMN_MAP: simple 1-to-1 column mappings (no errors)
    """
    keys = set(column_map.keys())
    if keys != REQUIRED_COLUMN_KEYS:
        raise ValueError(
            f"CIGALE_COLUMN_MAP must contain exactly these entries: {sorted(REQUIRED_COLUMN_KEYS)}. "
            f"Got: {sorted(keys)}"
        )
    columns = []
    for cigale_col, snana_col in column_map.items():
        columns.append({
            "cigale_col": cigale_col,
            "snana_col": snana_col,
        })
    if varname_z:
        columns.append({
            "cigale_col": 'redshift',
            "snana_col": varname_z,
        })
    return columns


def parse_mag_map(mag_map):
    """
    Parse CIGALE_MAG_MAP: each entry must have an SNANA column and an error source
    (either a column name or a depth value)
    """
    keys = set(mag_map.keys())
    if keys != REQUIRED_MAG_KEYS:
        raise ValueError(
            f"CIGALE_MAG_MAP must contain exactly these entries: {sorted(REQUIRED_MAG_KEYS)}. "
            f"Got: {sorted(keys)}"
        )
    columns = []
    for cigale_col, value_str in mag_map.items():
        tokens = str(value_str).split()
        if len(tokens) != 2:
            raise ValueError(
                f"CIGALE_MAG_MAP entry '{cigale_col}' must have exactly 2 values (snana_col + error source), got: '{value_str}'"
            )
        snana_col = tokens[0]
        try:
            depth = float(tokens[1])
            err_source = ("depth", depth)
        except ValueError:
            err_source = ("column", tokens[1])
        columns.append({
            "cigale_col": cigale_col,
            "snana_col": snana_col,
            "err_source": err_source,
        })
    return columns


def parse_redshift_grid(zgrid_map):
    """
    Parse REDSHIFT_GRID: zmin zmax numbins
    """
    #parts = zgrid_map.split()
    if len(zgrid_map) != 3:
        raise ValueError(
            "REDSHIFT_GRID must contain exactly 3 values: "
            "'zmin zmax numbins'"
        )
    zmin_str, zmax_str, numbins_str = zgrid_map
    try:
        zmin = float(zmin_str)
        zmax = float(zmax_str)
    except ValueError:
        raise ValueError(
            f"zmin and zmax must be numeric values. "
            f"Got: {zmin_str}, {zmax_str}"
        )
    try:
        numbins = int(numbins_str)
    except ValueError:
        raise ValueError(
            f"numbins must be an integer. Got: {numbins_str}"
        )
    if numbins <= 0:
        raise ValueError("numbins must be positive")
    return np.linspace(zmin, zmax, numbins)


def read_snana_input(filepath):
    """
    Reads snana input file and saves it as a table (dict)
    I'm assuming that all SNANA input files have columns labeled by VARNAMES and rows
    labeled by either GAL (hostlib) or SN (fitres)
    """
    opener = gzip.open if filepath.endswith(".gz") else open
    varnames = None
    rows = []
    with opener(filepath, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("VARNAMES:"):
                varnames = line.split()[1:]
                continue
            if line.startswith("GAL:") or line.startswith("SN:"):
                vals = line.split()[1:]
                rows.append(vals)
    if varnames is None:
        raise ValueError(f"No VARNAMES line found in {filepath}")
    data = {}
    for i, name in enumerate(varnames):
        col_vals = [row[i] for row in rows]
        if i == 0:
            try:
                data[name] = np.array(col_vals, dtype=int)
            except (ValueError, TypeError):
                # some columns have non-numeric entries (e.g., strings, lists)
                data[name] = np.array(col_vals, dtype=object)
        else:
            try:
                data[name] = np.array(col_vals, dtype=float)
            except (ValueError, TypeError):
                # some columns have non-numeric entries (e.g., strings, lists)
                data[name] = np.array(col_vals, dtype=object)
    return data


def create_output_header_and_cols(col_entries, mag_entries, data, zgrid):
    """
    Use parsed maps to create the header and columns for the cigale output file.
    REQUIRES REDSHIFT GRID: To deal with a redshift grid, this function now manually
    adds a redshift column and repeats all of the rows for each redshift in the grid.
    Note:
    np.repeat([1,2,3], 3) = [1,1,1,2,2,2,3,3,3]
    np.tile([1,2,3], 3) = [1,2,3,1,2,3,1,2,3]
    """
    header_parts = []
    output_cols = []
    snana_galids_tiled = None

    for col in col_entries:
        header_parts.append(col["cigale_col"])
        snana_vals = data[col["snana_col"]]
        if col["cigale_col"] == "id":
            n = len(snana_vals)
            cigale_ids = np.arange(n * len(zgrid), dtype=np.int64)
            snana_galids_tiled = np.tile(snana_vals.astype(np.int64), len(zgrid))
            output_cols.append(cigale_ids)
        else:
            output_cols.append(np.tile(snana_vals, len(zgrid)))

    # Redshift column
    header_parts.append('redshift')
    output_cols.append(np.repeat(zgrid, len(data[col["snana_col"]])))

    for col in mag_entries:
        mags = data[col["snana_col"]]
        flux_mjy = mag_to_mjy(mags)

        kind, val = col["err_source"]
        if kind == "depth":
            sigma_m = mag_err(mags, val)
        else:
            sigma_m = data[val]
        flux_err_mjy = mag_err_to_mjy_err(flux_mjy, sigma_m)

        header_parts.append(col["cigale_col"])
        output_cols.append(np.tile(flux_mjy, len(zgrid)))
        header_parts.append(col["cigale_col"] + "_err")
        output_cols.append(np.tile(flux_err_mjy, len(zgrid)))
    return header_parts, output_cols, snana_galids_tiled


def create_output_header_and_cols_nozgrid(col_entries, mag_entries, data):
    """
    Use parsed maps to create the header and columns for the cigale output file
    REQUIRES VARNAME_Z (i.e., NO REDSHIFT GRID) 
    """
    header_parts = []
    output_cols = []

    for col in col_entries:
        header_parts.append(col["cigale_col"])
        output_cols.append(data[col["snana_col"]])

    for col in mag_entries:
        mags = data[col["snana_col"]]
        flux_mjy = mag_to_mjy(mags)

        kind, val = col["err_source"]
        if kind == "depth":
            sigma_m = mag_err(mags, val)
        else:
            sigma_m = data[val]
        flux_err_mjy = mag_err_to_mjy_err(flux_mjy, sigma_m)

        header_parts.append(col["cigale_col"])
        output_cols.append(flux_mjy)
        header_parts.append(col["cigale_col"] + "_err")
        output_cols.append(flux_err_mjy)
    return header_parts, output_cols


def write_galid_map(galid_map_path, output_cols, snana_galids_tiled):
    """
    Write GALID map. Maps cigale ids to SNANA GALIDs. If using a redshift grid,
    the SNANA GALIDs will have duplicates, whereas cigale ids cannot have duplicates,
    hence the distinction (and therefore mapping) between these two sets of ids.
    """
    cigale_id_col = output_cols[0]
    with open(galid_map_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cigale_id", "snana_galid"])
        for cid, sgid in zip(cigale_id_col, snana_galids_tiled):
            w.writerow([int(cid), int(sgid)])
    print(f"Wrote galid map ({len(cigale_id_col)} rows) to {galid_map_path}")


def write_cigale_output(cigale_output_path, header_parts, output_cols):
    """
    Writes cigale file in the required format of a cigale input.
    """
    header_parts[0] = "#" + header_parts[0]
    nrows = len(output_cols[0])
    with open(cigale_output_path, "w") as f:
        f.write(" ".join(header_parts) + "\n")
        for i in range(nrows):
            vals = []
            for col in output_cols:
                if np.issubdtype(col.dtype, np.integer):
                    vals.append(f"{int(col[i]):d}")
                else:
                    vals.append(f"{col[i]:.6f}")
            f.write(" ".join(vals) + "\n")
    print(f"Wrote {nrows} rows to {cigale_output_path}")


def snana_to_cigale(args, config, mode):
    """
    Converts SNANA input file (e.g., HOSTLIB or FITRES) to a cigale input file.
    """
    # Required params for SNANA_TO_CIGALE
    snana_input_path = resolve(args.input_table_file, config, 'INPUT_TABLE_FILE', '--input_table_file', mode)
    cigale_output_path = resolve(args.output_cigale_file, config, 'OUTPUT_CIGALE_FILE', '--output_cigale_file', mode)

    # Other params
    galid_map_path = args.output_galid_map if args.output_galid_map is not None else config.get("OUTPUT_GALID_MAP")
    zgrid = args.zgrid if args.zgrid is not None else config.get("REDSHIFT_GRID").split()
    varname_z = args.varname_z if args.varname_z is not None else config.get("REDSHIFT_COL") 

    col_entries = parse_column_map(config["CIGALE_COLUMN_MAP"], varname_z = varname_z)
    mag_entries = parse_mag_map(config["CIGALE_MAG_MAP"])

    data = read_snana_input(snana_input_path)

    if zgrid:
        redshift_grid = parse_redshift_grid(zgrid)

        header_parts, output_cols, snana_galids_tiled = create_output_header_and_cols(
            col_entries, mag_entries, data, redshift_grid
        )

        if galid_map_path:
            write_galid_map(galid_map_path, output_cols, snana_galids_tiled)
        else:
            raise ValueError("galid_map_path must be provided if using zgrid")

    elif varname_z:
        header_parts, output_cols = create_output_header_and_cols_nozgrid(
            col_entries, mag_entries, data
        )

    else:
        raise ValueError(f"Either zgrid or varname_z must be provided for mode '{mode}'")

    write_cigale_output(cigale_output_path, header_parts, output_cols)


### CIGALE_TO_SNANA-specific functions ###

REQUIRED_COLUMN_KEYS_LOGMASS_OVERRIDE = {
    "GALID",
    "HOSTGALz_LOGMASS_ZGRID",
    "HOSTGALz_LOGMASS_VALGRID",
    "HOSTGALz_LOGMASS_ERRGRID",
}


def build_doc_header(snana_cols):
    varnames_line = "VARNAMES:  " + "  ".join(snana_cols)
    return (
        "DOCUMENTATION:\n"
        "\n"
        "DOCUMENTATION_END:\n"
        "\n"
        "\n"
        f"{varnames_line}\n"
        "\n"
        "\n"
    )


def validate_column_map(column_map, snana_format):
    """
    Validate COLUMN_MAP against the required keys for the given SNANA_FORMAT.
    For now only LOGMASS_OVERRIDE is supported.
    """
    if snana_format != "LOGMASS_OVERRIDE":
        raise ValueError(
            f"SNANA_FORMAT '{snana_format}' is not supported. "
            f"Currently only 'LOGMASS_OVERRIDE' is implemented."
        )
    keys = set(column_map.keys())
    if keys != REQUIRED_COLUMN_KEYS_LOGMASS_OVERRIDE:
        raise ValueError(
            f"COLUMN_MAP for SNANA_FORMAT=LOGMASS_OVERRIDE must contain exactly "
            f"these entries: {sorted(REQUIRED_COLUMN_KEYS_LOGMASS_OVERRIDE)}. "
            f"Got: {sorted(keys)}"
        )
    return dict(column_map)


def load_galid_map(path):
    """
    Load the sidecar CSV written by SNANA_TO_CIGALE.
    Returns a dict {cigale_id (int): snana_galid (int)}.
    """
    mapping = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mapping[int(row["cigale_id"])] = int(row["snana_galid"])
    return mapping


def read_cigale_results(filepath, cigale_cols):
    """
    Read needed columns from a Cigale results.fits file (BinTableHDU at HDU 1).
    Returns a dict mapping cigale column name -> numpy array.
    """
    with fits.open(filepath) as hdul:
        table = hdul[1].data
        out = {}
        for col in cigale_cols:
            out[col] = np.asarray(table[col])
    return out


def linear_mass_to_logmass(mass, mass_err):
    """
    Convert linear stellar mass (solMass) and its error to log10(M) and the
    propagated error sigma_logM = sigma_M / (M * ln 10).
    Non-finite or non-positive masses produce NaN; those get replaced with -9
    later when the output is written.
    """
    mass = np.asarray(mass, dtype=float)
    mass_err = np.asarray(mass_err, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        logmass = np.where(mass > 0, np.log10(mass), np.nan)
        logmass_err = np.where(mass > 0, mass_err / (mass * np.log(10)), np.nan)
    return logmass, logmass_err


def fmt_val(x):
    """Format a numeric value for the output file, using -9 for non-finite."""
    if not np.isfinite(x):
        return "-9"
    return f"{x:.6f}"


def write_snana_output(output_path, snana_cols, galids, redshifts, logmass, logmass_err, num_pad_rows=0):
    """
    Write the SNANA LOGMASS_OVERRIDE-format file. Galaxies appear in order of
    first appearance in `galids`; within each galaxy, rows are sorted by
    ascending redshift, followed by `num_pad_rows` of all-'-9' padding.
    """
    galids_int = galids.astype(np.int64)

    seen_order = list(dict.fromkeys(galids_int.tolist()))

    with open(output_path, "w") as f:
        f.write(build_doc_header(snana_cols))
        for i, gid in enumerate(seen_order):
            mask = galids_int == gid
            z_g = redshifts[mask]
            m_g = logmass[mask]
            e_g = logmass_err[mask]

            order = np.argsort(z_g)
            z_g = z_g[order]
            m_g = m_g[order]
            e_g = e_g[order]

            for z, m, e in zip(z_g, m_g, e_g):
                f.write(f"GAL: {gid}  {fmt_val(z)}  {fmt_val(m)}  {fmt_val(e)}\n")

            for _ in range(num_pad_rows):
                f.write(f"GAL: {gid}  -9  -9  -9\n")

            if i != len(seen_order) - 1:
                f.write("\n")

    total_real_rows = int(np.sum(np.isfinite(redshifts)))
    return len(seen_order), total_real_rows


def cigale_to_snana(args, config, mode):
    """
    Converts cigale output file (results.fits) to a SNANA LOGMASS_OVERRIDE file.
    """
    # Required params for CIGALE_TO_SNANA
    cigale_results_path = resolve(args.input_cigale_results, config, 'INPUT_CIGALE_RESULTS', '--input_cigale_results', mode)
    output_path = resolve(args.output_snana_file, config, 'OUTPUT_SNANA_FILE', '--output_snana_file', mode)

    # Other params
    galid_map_path = args.output_galid_map if args.output_galid_map is not None else config.get("INPUT_GALID_MAP")

    # Currently, snana_format needs to be given in config
    snana_format = config["SNANA_FORMAT"]
    column_map = validate_column_map(config["COLUMN_MAP"], snana_format)

    cigale_id_col = column_map["GALID"]
    cigale_z_col = column_map["HOSTGALz_LOGMASS_ZGRID"]
    cigale_mass_col = column_map["HOSTGALz_LOGMASS_VALGRID"]
    cigale_mass_err_col = column_map["HOSTGALz_LOGMASS_ERRGRID"]

    data = read_cigale_results(
        cigale_results_path,
        [cigale_id_col, cigale_z_col, cigale_mass_col, cigale_mass_err_col],
    )

    if galid_map_path:
        galid_map = load_galid_map(galid_map_path)
        cigale_ids = data[cigale_id_col].astype(np.int64)
        missing = [int(c) for c in cigale_ids if int(c) not in galid_map]
        if missing:
            raise ValueError(
                f"{len(missing)} cigale ids from {cigale_results_path} are not in "
                f"{galid_map_path} (e.g., {missing[:5]}). Make sure the map matches the run."
            )
        galids = np.array([galid_map[int(c)] for c in cigale_ids], dtype=np.int64)
    else:
        # If no galid map path is provided, that means the cigale ids are
        # already equal to the galids
        galids = data[cigale_id_col].astype(np.int64)

    redshifts = data[cigale_z_col].astype(float)
    logmass, logmass_err = linear_mass_to_logmass(
        data[cigale_mass_col], data[cigale_mass_err_col]
    )

    snana_cols = list(column_map.keys())
    n_gal, n_rows = write_snana_output(
        output_path, snana_cols, galids, redshifts, logmass, logmass_err
    )
    print(f"Wrote {n_gal} galaxies ({n_rows} grid rows) to {output_path}")


if __name__ == "__main__":

    # Read command-line arguments
    args = get_args()

    # Read 'mode' (i.e., SNANA_TO_CIGALE or CIGALE_TO_SNANA)
    # and config file
    mode = args.mode
    config = parse_config(args.config, mode)

    if mode == 'SNANA_TO_CIGALE':
        #Run SNANA_TO_CIGALE
        snana_to_cigale(args, config, mode)

    elif mode == 'CIGALE_TO_SNANA':
        # Run CIGALE_TO_SNANA
        cigale_to_snana(args, config, mode)

    else:
        sys.exit(f"mode must be 'SNANA_TO_CIGALE' or 'CIGALE_TO_SNANA', got {mode!r}")


    print(f"DONE.")

    # === END: ====
