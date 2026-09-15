#!/usr/bin/env python3
""" Created Summer 2026 by  A.Mitra
    Installed into SNANA, Sep 2026

    Prepare a DECam photometry override table for SNANA (LastJourney layout).

   TO-DO list for SNANA:
     * change reader to glob.glob for lc-core files using optional wildcard;
     * add --wildcard (-w) areg to select small subset of lc-cores for quick test

     * replace command-line inputs with config_file that includes filter definitions/file
     * remove FILTER_NAMES and COLUMN_NAMED 
     * determin --mock-version from catalog basename, instead of separate input
     * remove reference to DES/DECam; use more generic language

"""

from __future__ import annotations

import argparse
from collections import namedtuple
from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError
import json
from pathlib import Path
import sqlite3
import tempfile

import numpy as np


FILTER_NAMES = ("decam_g", "decam_r", "decam_i", "decam_z", "decam_Y")
COLUMN_NAMES = ("des_g", "des_r", "des_i", "des_z", "des_Y")
DECamCurves = namedtuple("DECamCurves", COLUMN_NAMES)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog-dir', type=Path, required=True)
    p.add_argument('--model-dir', type=Path, help='Default: catalog directory')
    p.add_argument('--mock-version', required=True)
    p.add_argument('--output-dir', type=Path, required=True, help='New directory; never overwritten')
    p.add_argument('--z-min', type=float, required=True)
    p.add_argument('--z-max', type=float, required=True)
    p.add_argument('--grid-size', type=int, default=200)
    p.add_argument('--batch-size', type=int, default=10000)
    p.add_argument('--scatter-policy', choices=['catalog', 'zero'], default='catalog')
    args = p.parse_args(argv)
    if not (0 < args.z_min < args.z_max < np.inf):
        p.error('Require finite 0 < z-min < z-max')
    if args.grid_size < 2 or args.batch_size < 1:
        p.error('grid-size must be >= 2 and batch-size >= 1')
    return args


def compute_decam_photometry(
        diffsky_data, *, ssp_data, param_collection, sim_info, z_phot_table,
        tcurves=None, scatter_policy="catalog", n_ssp_cols=None, ):
    
    """Return row-aligned DES magnitude columns using the Diffsky engine.

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
        sim_info=sim_info, z_phot_table=grid,
        tcurves=build_decam_curves() if tcurves is None else tcurves,
    )
    mags = np.asarray(result["obs_mags"])
    if mags.shape != (len(z), len(COLUMN_NAMES)) or not np.all(np.isfinite(mags)):
        raise ValueError("Diffsky returned invalid magnitudes: expected finite (N_gal, 5)")
    return {name: mags[:, i].copy() for i, name in enumerate(COLUMN_NAMES)}



def read_patch(path, scatter_policy):
    # Register compression plugins before opening catalog data.
    import hdf5plugin  # noqa: F401
    import opencosmo as oc
    from diffmah import DEFAULT_MAH_PARAMS
    from diffstar import DEFAULT_DIFFSTAR_PARAMS
    columns = list(dict.fromkeys([
        *DEFAULT_MAH_PARAMS._fields, *DEFAULT_DIFFSTAR_PARAMS._fields,
        'redshift_true', 'mc_sfh_type', 'uran_av', 'uran_delta', 'uran_funo',
        'uran_pburst', 'gal_id',
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


def write_batches(patches, args, work, model, curves, grid, n_ssp, reader=read_patch,
                  compute=compute_decam_photometry):
    import pyarrow as pa
    import pyarrow.parquet as pq
    schema = pa.schema([('serial_tag', pa.int64()), *[(name, pa.float64()) for name in COLUMN_NAMES]])
    counts = []
    with sqlite3.connect(work / 'ids.sqlite') as db, pq.ParquetWriter(work / 'photometry.parquet', schema) as writer:
        db.execute('CREATE TABLE ids (id INTEGER PRIMARY KEY)')
        for path in patches:
            data = reader(path, args.scatter_policy)
            z = np.asarray(data['redshift_true'])
            if z.ndim != 1 or not np.all(np.isfinite(z)):
                raise ValueError(f'Invalid catalog redshifts: {path}')
            if any(np.asarray(value).ndim == 0 or len(value) != len(z) for value in data.values()):
                raise ValueError(f'Input columns are not row aligned: {path}')
            selected = np.flatnonzero((z >= args.z_min) & (z <= args.z_max))
            for start in range(0, len(selected), args.batch_size):
                take = selected[start:start + args.batch_size]
                batch = {k: v[take] for k, v in data.items()}
                ids = register_ids(db, batch.pop('gal_id'))
                mags = compute(batch, **model, tcurves=curves, z_phot_table=grid,
                               scatter_policy=args.scatter_policy, n_ssp_cols=n_ssp)
                writer.write_table(pa.Table.from_pydict({'serial_tag': ids, **mags}, schema=schema))
            counts.append({'path': str(path.resolve()), 'selected_rows': len(selected),
                           'size_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            print(f'{path.name}: {len(selected):,} selected galaxies', flush=True)
    (work / 'ids.sqlite').unlink()
    if not sum(item['selected_rows'] for item in counts):
        raise ValueError('No galaxies selected; no output published')
    return counts



"""DECam photometry adapter for Diffsky's existing galaxy photometry engine.

No catalog paths, selection cuts, geometry prescriptions, or SNANA dependency.
The caller supplies row-aligned Diffsky arrays and matching model objects.
"""



def build_decam_curves():
    """Load sedpy DECam throughput curves; wavelengths are in Angstroms."""
    from sedpy.observate import load_filters
    from dsps.data_loaders.defaults import TransmissionCurve

    return DECamCurves(*[
        TransmissionCurve(
            wave=np.asarray(f.wavelength, dtype=np.float64),
            transmission=np.asarray(f.transmission, dtype=np.float64),
        ) for f in load_filters(list(FILTER_NAMES))
    ])



def main(argv=None):
    args = parse_args(argv)
    args.catalog_dir = args.catalog_dir.resolve()
    if not args.catalog_dir.is_dir():
        raise ValueError('catalog-dir must be an existing release directory')
    if args.output_dir.exists():
        raise FileExistsError(f'Output already exists: {args.output_dir}')
    # Only direct children of the explicitly supplied catalog directory.
    patches = sorted(args.catalog_dir.glob('lc_cores-*.diffsky_gals.hdf5'))
    if not patches:
        raise FileNotFoundError('No lc_cores-*.diffsky_gals.hdf5 files in catalog-dir')
    model, cosmology, n_ssp = load_model(args, patches[0])
    for patch in patches:
        check_patch_metadata(patch, cosmology, n_ssp)
    curves = build_decam_curves()
    grid = np.linspace(args.z_min, args.z_max, args.grid_size)
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.decam-', dir=args.output_dir.parent) as temp:
        work = Path(temp) / 'result'
        work.mkdir()
        counts = write_batches(patches, args, work, model, curves, grid, n_ssp)
        versions = {}
        for package in ('numpy', 'sedpy', 'dsps', 'diffsky', 'diffmah', 'diffstar',
                        'jax', 'jaxlib', 'opencosmo', 'h5py', 'hdf5plugin', 'pyarrow'):
            try:
                versions[package] = version(package)
            except PackageNotFoundError:
                versions[package] = 'not recorded by package metadata'
        arrays = {'z_phot_table': grid}
        for name, curve in zip(FILTER_NAMES, curves):
            arrays[name + '_wave'] = curve.wave
            arrays[name + '_transmission'] = curve.transmission
        np.savez(work / 'filters_and_grid.npz', **arrays)
        metadata = {'created_utc': datetime.now(timezone.utc).isoformat(),
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
    
