#!/usr/bin/env python3
""" Created Summer 2026 by  A.Mitra
    Installed into SNANA, Sep 2026

    Prepare a synthetic photometry override table for SNANA (LastJourney layout).

    Usage:
      python diffsky_synthetic_photometry.py config.yml [-w WILDCARD]

    Example config.yml:
      CATALOG_DIR:  /path/to/diffsky/catalog_release
      MODEL_DIR:    /path/to/matching/model_files   # optional; default = CATALOG_DIR
      MOCK_VERSION: catalog_release_name            # optional; default = basename(CATALOG_DIR)
      OUTPUT_DIR:   /path/to/output                 # new directory; never overwritten
      Z_MIN: 0.05
      Z_MAX: 1.8
      GRID_SIZE:   200          # optional, default 200
      BATCH_SIZE:  10000        # optional, default 10000
      SCATTER_POLICY: catalog   # optional, default catalog; choices: catalog, zero

      CUTWIN:     # optional; selection cuts on diffsky-catalog columns, applied before
                  # photometry synthesis so galaxies that would be discarded downstream
                  # are never photometered (mirrors make_hostlib_diffsky.py's CUTWIN)
      - COLUMN: ra
        MIN: 240.0
        MAX: 245.0
      - COLUMN: dec
        MIN: 53.0
        MAX: 58.0
      - COLUMN: logsm_obs
        MIN: 8.0             # MIN and/or MAX; at least one is required per entry

      FILTERS:    # sedpy filter id -> output column name, in output column order
      - SEDPY_ID: decam_g
        COLUMN:   des_g
      - SEDPY_ID: decam_r
        COLUMN:   des_r
      # ...or, instead of an inline FILTERS list:
      # FILTER_FILE: /path/to/filters.yml   # a YAML file holding its own FILTERS list
"""

from __future__ import annotations

import argparse
from collections import namedtuple
from datetime import datetime, timezone
import glob
from importlib.metadata import version, PackageNotFoundError
import json
from pathlib import Path
import sqlite3
import tempfile
from types import SimpleNamespace

import numpy as np
import yaml


def read_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def load_filters(config):
    """Return (sedpy_ids, column_names) from the config's FILTERS list or FILTER_FILE."""
    if 'FILTERS' in config:
        entries = config['FILTERS']
    elif 'FILTER_FILE' in config:
        entries = read_yaml(config['FILTER_FILE'])['FILTERS']
    else:
        raise ValueError('config_file must define FILTERS (inline list) or '
                          'FILTER_FILE (path to a YAML file with a FILTERS list)')
    if not entries:
        raise ValueError('FILTERS must list at least one filter')
    sedpy_ids = tuple(str(e['SEDPY_ID']) for e in entries)
    column_names = tuple(str(e['COLUMN']) for e in entries)
    if len(set(column_names)) != len(column_names):
        raise ValueError('FILTERS COLUMN names must be unique')
    return sedpy_ids, column_names


def load_cutwin(config):
    """Return ((column, min_or_None, max_or_None), ...) from the config's optional CUTWIN list.

    Applied before photometry synthesis so galaxies that would be cut downstream are
    never photometered; mirrors make_hostlib_diffsky.py's CUTWIN.
    """
    cuts = []
    for entry in config.get('CUTWIN', []):
        column = str(entry['COLUMN'])
        lo = entry.get('MIN')
        hi = entry.get('MAX')
        if lo is None and hi is None:
            raise ValueError(f"CUTWIN entry for column {column!r} must set MIN and/or MAX")
        cuts.append((column, None if lo is None else float(lo), None if hi is None else float(hi)))
    return tuple(cuts)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('config_file', type=Path,
                    help='YAML config with catalog paths, redshift range, and filter definitions')
    p.add_argument('--wildcard', '-w', type=str, default=None,
                    help='Optional substring to select a subset of lc_cores-*.diffsky_gals.hdf5 '
                         'files (e.g. for a quick test), matched as lc_cores-*<wildcard>*.diffsky_gals.hdf5')
    return p.parse_args(argv)


def build_run_args(config, wildcard):
    """Validate the YAML config and merge it with CLI overrides into a run-argument namespace."""
    missing = [key for key in ('CATALOG_DIR', 'OUTPUT_DIR', 'Z_MIN', 'Z_MAX') if key not in config]
    if missing:
        raise ValueError(f'config_file is missing required key(s): {", ".join(missing)}')
    z_min, z_max = float(config['Z_MIN']), float(config['Z_MAX'])
    if not (0 < z_min < z_max < np.inf):
        raise ValueError('Require finite 0 < Z_MIN < Z_MAX')
    grid_size = int(config.get('GRID_SIZE', 200))
    batch_size = int(config.get('BATCH_SIZE', 10000))
    if grid_size < 2 or batch_size < 1:
        raise ValueError('GRID_SIZE must be >= 2 and BATCH_SIZE >= 1')
    scatter_policy = config.get('SCATTER_POLICY', 'catalog')
    if scatter_policy not in ('catalog', 'zero'):
        raise ValueError("SCATTER_POLICY must be 'catalog' or 'zero'")
    sedpy_ids, column_names = load_filters(config)
    cutwin = load_cutwin(config)
    return SimpleNamespace(
        catalog_dir=Path(config['CATALOG_DIR']),
        model_dir=Path(config['MODEL_DIR']) if config.get('MODEL_DIR') else None,
        mock_version=str(config['MOCK_VERSION']) if config.get('MOCK_VERSION') else None,
        output_dir=Path(config['OUTPUT_DIR']),
        z_min=z_min, z_max=z_max, grid_size=grid_size, batch_size=batch_size,
        scatter_policy=scatter_policy, wildcard=wildcard, cutwin=cutwin,
        sedpy_ids=sedpy_ids, column_names=column_names,
    )


def compute_synthetic_photometry(
        diffsky_data, *, ssp_data, param_collection, sim_info, z_phot_table, tcurves,
        column_names, scatter_policy="catalog", n_ssp_cols=None, ):

    """Return row-aligned synthetic magnitude columns using the Diffsky engine.

    Model objects and input fields follow compute_phot_from_diffsky_mock.
    Redshifts must be positive and bracketed by the supplied interpolation grid.
    ``catalog`` requires an existing finite (N, N_ssp) scatter array.
    ``zero`` explicitly reproduces the local production approximation and
    requires n_ssp_cols from the matching catalog/template metadata. Zero
    scatter is not asserted to equal mean flux or mean magnitude.

    Does not modify the input, reorder rows, assign IDs, add observational
    noise, or apply an additional extinction or magnification correction.
    Physical conventions are inherited from the supplied Diffsky engine.
    """
    data = dict(diffsky_data)
    z = np.asarray(data["redshift_true"])
    grid = np.asarray(z_phot_table)
    if z.ndim != 1 or not len(z) or not np.all(np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("redshift_true must be a nonempty finite positive 1D array")
    if (grid.ndim != 1 or len(grid) < 2 or not np.all(np.isfinite(grid))
            or np.any(grid <= 0) or np.any(np.diff(grid) <= 0)):
        raise ValueError("z_phot_table must be finite, positive, and strictly increasing")
    if z.min() < grid[0] or z.max() > grid[-1]:
        raise ValueError("z_phot_table must bracket all galaxy redshifts")
    if scatter_policy == "zero":
        if not isinstance(n_ssp_cols, int) or isinstance(n_ssp_cols, bool) or n_ssp_cols < 1:
            raise ValueError("zero scatter requires a positive integer n_ssp_cols")
        data["delta_mag_ssp_scatter"] = np.zeros((len(z), n_ssp_cols), dtype=np.float32)
    elif scatter_policy != "catalog":
        raise ValueError("scatter_policy must be 'catalog' or 'zero'")
    if "delta_mag_ssp_scatter" not in data:
        raise ValueError("catalog scatter is required; zero scatter must be explicitly requested")
    scatter = np.asarray(data["delta_mag_ssp_scatter"])
    if (scatter.ndim != 2 or scatter.shape[0] != len(z) or scatter.shape[1] < 1
            or not np.all(np.isfinite(scatter))):
        raise ValueError("delta_mag_ssp_scatter must be finite with shape (N_gal, N_ssp)")

    from diffsky.data_loaders.hacc_utils import load_lc_mock
    engine = getattr(load_lc_mock, "compute_phot_from_diffsky_mock", None)
    if not callable(engine):
        raise RuntimeError(
            "This adapter requires Diffsky's load_lc_mock.compute_phot_from_diffsky_mock API. "
            "Use the environment matching the catalog production release; "
            "newer Diffsky APIs are not interchangeable."
        )
    result = engine(
        diffsky_data=data, ssp_data=ssp_data, param_collection=param_collection,
        sim_info=sim_info, z_phot_table=grid, tcurves=tcurves,
    )
    mags = np.asarray(result["obs_mags"])
    if mags.shape != (len(z), len(column_names)) or not np.all(np.isfinite(mags)):
        raise ValueError(f"Diffsky returned invalid magnitudes: expected finite (N_gal, {len(column_names)})")
    return {name: mags[:, i].copy() for i, name in enumerate(column_names)}



def read_patch(path, scatter_policy, extra_columns=()):
    # Register compression plugins before opening catalog data.
    import hdf5plugin  # noqa: F401
    import opencosmo as oc
    from diffmah import DEFAULT_MAH_PARAMS
    from diffstar import DEFAULT_DIFFSTAR_PARAMS
    columns = list(dict.fromkeys([
        *DEFAULT_MAH_PARAMS._fields, *DEFAULT_DIFFSTAR_PARAMS._fields,
        'redshift_true', 'mc_sfh_type', 'uran_av', 'uran_delta', 'uran_funo',
        'uran_pburst', 'gal_id', *extra_columns,
    ]))
    if scatter_policy == 'catalog':
        columns.append('delta_mag_ssp_scatter')
    catalog = oc.open(path, synth_cores=True)
    try:
        table = catalog.select(columns).get_data()
        # Preserve multidimensional scatter; pandas cannot represent it directly.
        return {name: np.asarray(table[name]) for name in columns}
    finally:
        close = getattr(catalog, 'close', None)
        if close is not None:
            close()


def load_model(args, patch):
    import hdf5plugin  # noqa: F401
    import h5py
    from diffsky.data_loaders.hacc_utils import lc_mock, load_lc_mock
    if not callable(getattr(load_lc_mock, 'compute_phot_from_diffsky_mock', None)):
        raise RuntimeError('Installed Diffsky lacks compute_phot_from_diffsky_mock; '
                           'use the environment matching the catalog production release.')
    from dsps.cosmology.flat_wcdm import CosmoParams
    with h5py.File(patch, 'r') as hdf:
        attrs = hdf['header/simulation/cosmology'].attrs
        cosmology = {key: float(attrs[key]) for key in ('omega_b', 'omega_m', 'h', 'w_0', 'w_a')}
        n_ssp = int(hdf['cores/data/delta_mag_ssp_scatter'].shape[1])
    cp = CosmoParams(Om0=cosmology['omega_m'], w0=cosmology['w_0'],
                     wa=cosmology['w_a'], h=cosmology['h'])
    sim = namedtuple('SimInfo', ['cosmo_params', 'fb'])(cp, cosmology['omega_b'] / cosmology['omega_m'])
    model_dir = str(args.model_dir or args.catalog_dir)
    return dict(ssp_data=lc_mock.load_diffsky_ssp_data(model_dir, args.mock_version),
                param_collection=lc_mock.load_diffsky_param_collection(model_dir, args.mock_version),
                sim_info=sim), cosmology, n_ssp


def check_patch_metadata(path, cosmology, n_ssp):
    import h5py
    with h5py.File(path, 'r') as hdf:
        attrs = hdf['header/simulation/cosmology'].attrs
        if any(float(attrs[k]) != v for k, v in cosmology.items()):
            raise ValueError(f'Inconsistent cosmology: {path}')
        if hdf['cores/data/delta_mag_ssp_scatter'].shape[1] != n_ssp:
            raise ValueError(f'Inconsistent SSP scatter dimension: {path}')


def register_ids(db, ids):
    ids = np.asarray(ids)
    if ids.ndim != 1 or ids.dtype.kind not in 'iu':
        raise ValueError('gal_id must be a one-dimensional integer array; float IDs are unsafe')
    if ids.size and (int(ids.min()) < -(2**63) or int(ids.max()) >= 2**63):
        raise ValueError('gal_id does not fit signed int64')
    try:
        db.executemany('INSERT INTO ids VALUES (?)', ((int(i),) for i in ids))
        db.commit()
    except sqlite3.IntegrityError as exc:
        raise ValueError('Duplicate gal_id: cannot safely use the SNANA serial_tag join. '
                         'Agree on a composite key with the catalog/SNANA maintainers.') from exc
    return ids.astype(np.int64, copy=False)


def select_rows(data, z_min, z_max, cutwin):
    """Row indices passing the redshift range and all CUTWIN cuts."""
    z = np.asarray(data['redshift_true'])
    mask = (z >= z_min) & (z <= z_max)
    for column, lo, hi in cutwin:
        values = np.asarray(data[column])
        if lo is not None:
            mask &= values >= lo
        if hi is not None:
            mask &= values <= hi
    return np.flatnonzero(mask)


def write_batches(patches, args, work, model, curves, grid, n_ssp, column_names, reader=read_patch,
                  compute=compute_synthetic_photometry):
    import gc
    import jax
    import pyarrow as pa
    import pyarrow.parquet as pq
    schema = pa.schema([('serial_tag', pa.int64()), *[(name, pa.float64()) for name in column_names]])
    # CUTWIN columns are only needed to build the selection mask; they are not
    # inputs to the photometry engine, so drop them from each batch before compute().
    cutwin_only_columns = tuple(column for column, _, _ in args.cutwin)
    counts = []
    with sqlite3.connect(work / 'ids.sqlite') as db, pq.ParquetWriter(work / 'photometry.parquet', schema) as writer:
        db.execute('CREATE TABLE ids (id INTEGER PRIMARY KEY)')
        for path in patches:
            data = reader(path, args.scatter_policy, cutwin_only_columns)
            z = np.asarray(data['redshift_true'])
            if z.ndim != 1 or not np.all(np.isfinite(z)):
                raise ValueError(f'Invalid catalog redshifts: {path}')
            if any(np.asarray(value).ndim == 0 or len(value) != len(z) for value in data.values()):
                raise ValueError(f'Input columns are not row aligned: {path}')
            selected = select_rows(data, args.z_min, args.z_max, args.cutwin)
            data = {k: v for k, v in data.items() if k not in cutwin_only_columns}
            for start in range(0, len(selected), args.batch_size):
                take = selected[start:start + args.batch_size]
                batch = {k: v[take] for k, v in data.items()}
                ids = register_ids(db, batch.pop('gal_id'))
                mags = compute(batch, **model, tcurves=curves, z_phot_table=grid,
                               column_names=column_names,
                               scatter_policy=args.scatter_policy, n_ssp_cols=n_ssp)
                writer.write_table(pa.Table.from_pydict({'serial_tag': ids, **mags}, schema=schema))
            counts.append({'path': str(path.resolve()), 'selected_rows': len(selected),
                           'size_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            print(f'{path.name}: {len(selected):,} selected galaxies', flush=True)
            # Each patch has a different N_gal, so JAX JIT-compiles a new kernel per
            # patch; clear the cache to prevent unbounded growth over many patches.
            jax.clear_caches()
            gc.collect()
    (work / 'ids.sqlite').unlink()
    if not sum(item['selected_rows'] for item in counts):
        raise ValueError('No galaxies selected; no output published')
    return counts



"""Synthetic-photometry adapter for Diffsky's existing galaxy photometry engine.

No catalog paths, selection cuts, geometry prescriptions, or SNANA dependency.
The caller supplies row-aligned Diffsky arrays and matching model objects.
"""



def build_photometry_curves(sedpy_ids, column_names):
    """Load sedpy throughput curves for the configured filters; wavelengths are in Angstroms."""
    from sedpy.observate import load_filters
    from dsps.data_loaders.defaults import TransmissionCurve

    Curves = namedtuple('Curves', column_names)
    return Curves(*[
        TransmissionCurve(
            wave=np.asarray(f.wavelength, dtype=np.float64),
            transmission=np.asarray(f.transmission, dtype=np.float64),
        ) for f in load_filters(list(sedpy_ids))
    ])



def main(argv=None):
    cli = parse_args(argv)
    config = read_yaml(cli.config_file)
    args = build_run_args(config, cli.wildcard)
    args.catalog_dir = args.catalog_dir.resolve()
    if not args.catalog_dir.is_dir():
        raise ValueError('CATALOG_DIR must be an existing release directory')
    if args.mock_version is None:
        args.mock_version = args.catalog_dir.name
    if args.output_dir.exists():
        raise FileExistsError(f'Output already exists: {args.output_dir}')
    # Only direct children of the explicitly supplied catalog directory.
    if args.wildcard:
        pattern = f'lc_cores-*{args.wildcard}*.diffsky_gals.hdf5'
    else:
        pattern = 'lc_cores-*.diffsky_gals.hdf5'
    patches = sorted(Path(p) for p in glob.glob(str(args.catalog_dir / pattern)))
    if not patches:
        raise FileNotFoundError(f'No {pattern} files in CATALOG_DIR')
    model, cosmology, n_ssp = load_model(args, patches[0])
    for patch in patches:
        check_patch_metadata(patch, cosmology, n_ssp)
    curves = build_photometry_curves(args.sedpy_ids, args.column_names)
    grid = np.linspace(args.z_min, args.z_max, args.grid_size)
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.photometry-', dir=args.output_dir.parent) as temp:
        work = Path(temp) / 'result'
        work.mkdir()
        counts = write_batches(patches, args, work, model, curves, grid, n_ssp, args.column_names)
        versions = {}
        for package in ('numpy', 'sedpy', 'dsps', 'diffsky', 'diffmah', 'diffstar',
                        'jax', 'jaxlib', 'opencosmo', 'h5py', 'hdf5plugin', 'pyarrow'):
            try:
                versions[package] = version(package)
            except PackageNotFoundError:
                versions[package] = 'not recorded by package metadata'
        arrays = {'z_phot_table': grid}
        for name, curve in zip(args.column_names, curves):
            arrays[name + '_wave'] = curve.wave
            arrays[name + '_transmission'] = curve.transmission
        np.savez(work / 'filters_and_grid.npz', **arrays)
        metadata = {'created_utc': datetime.now(timezone.utc).isoformat(),
                    'config_file': str(cli.config_file.resolve()),
                    'arguments': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    'model_dir': str((args.model_dir or args.catalog_dir).resolve()),
                    'cosmology': cosmology, 'n_ssp_cols': n_ssp, 'versions': versions,
                    'patches': counts, 'rows': sum(p['selected_rows'] for p in counts),
                    'synthetic_cores': True, 'join': 'gal_id -> serial_tag (int64; unique in selected rows)',
                    'scientific_validation': 'Not certified by this driver'}
        (work / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
        if args.output_dir.exists():
            raise FileExistsError(args.output_dir)
        work.rename(args.output_dir)
    print(f'Created {args.output_dir}/photometry.parquet')

    return


if __name__ == '__main__':
    main()
    
