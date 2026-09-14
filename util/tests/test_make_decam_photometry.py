"""Driver regression tests with real Parquet output and a substitute photometry engine."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from types import SimpleNamespace
import sqlite3

import numpy as np
import pyarrow.parquet as pq
import pytest

from des_photometry import COLUMN_NAMES
from make_decam_photometry import parse_args, register_ids, write_batches


def test_integer_ids_are_exact_and_duplicates_rejected():
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE ids (id INTEGER PRIMARY KEY)')
        ids = np.array([2**53 + 1, 2**53 + 2], dtype=np.int64)
        np.testing.assert_array_equal(register_ids(db, ids), ids)
        with pytest.raises(ValueError, match='Duplicate gal_id'):
            register_ids(db, ids[:1])
        with pytest.raises(ValueError, match='float IDs'):
            register_ids(db, np.array([1.0]))


def test_batches_keep_scatter_and_large_ids_aligned(tmp_path):
    patch = tmp_path / 'patch.hdf5'
    patch.touch()
    args = SimpleNamespace(z_min=0.1, z_max=1.0, batch_size=2, scatter_policy='catalog')
    data = {'gal_id': np.array([2**53+1, 2**53+2, 2**53+3, 2**53+4], dtype=np.int64),
            'redshift_true': np.array([0.5, 0.01, 0.8, 1.0]),
            'delta_mag_ssp_scatter': np.arange(12).reshape(4, 3)}
    seen = []

    def compute(batch, **kwargs):
        seen.append(batch['delta_mag_ssp_scatter'].copy())
        return {name: batch['redshift_true'] + i + 20 for i, name in enumerate(COLUMN_NAMES)}

    counts = write_batches([patch], args, tmp_path, {}, None, None, 3,
                          reader=lambda *a: data, compute=compute)
    out = pq.read_table(tmp_path / 'photometry.parquet').to_pandas()
    assert counts[0]['selected_rows'] == 3
    np.testing.assert_array_equal(out.serial_tag, data['gal_id'][[0, 2, 3]])
    np.testing.assert_array_equal(np.concatenate(seen), data['delta_mag_ssp_scatter'][[0, 2, 3]])
    np.testing.assert_allclose(out.des_g, [20.5, 20.8, 21])
    assert list(out) == ['serial_tag', *COLUMN_NAMES]


def test_duplicate_across_patches_fails(tmp_path):
    paths = [tmp_path / name for name in ['a', 'b']]
    for path in paths:
        path.touch()
    args = SimpleNamespace(z_min=0.1, z_max=1, batch_size=1, scatter_policy='zero')
    with pytest.raises(ValueError, match='Duplicate gal_id'):
        write_batches(paths, args, tmp_path, {}, None, None, 3,
                      reader=lambda *a: {'gal_id': np.array([7]), 'redshift_true': np.array([0.5])},
                      compute=lambda batch, **kw: {name: np.array([22.]) for name in COLUMN_NAMES})


def test_empty_selection_fails(tmp_path):
    path = tmp_path / 'empty'
    path.touch()
    args = SimpleNamespace(z_min=0.1, z_max=1, batch_size=1, scatter_policy='catalog')
    with pytest.raises(ValueError, match='No galaxies selected'):
        write_batches([path], args, tmp_path, {}, None, None, 3,
                      reader=lambda *a: {'gal_id': np.array([7]), 'redshift_true': np.array([2.])})


@pytest.mark.parametrize('extra', [['--batch-size', '0'], ['--grid-size', '1'], ['--z-min', 'nan']])
def test_invalid_cli_settings(extra):
    with pytest.raises(SystemExit):
        parse_args(['--catalog-dir', '.', '--mock-version', 'test', '--output-dir', 'out',
                    '--z-min', '0.05', '--z-max', '1.8', *extra])


@pytest.mark.parametrize('duplicate', [False, True])
def test_main_publishes_complete_output_only(tmp_path, monkeypatch, duplicate):
    import make_decam_photometry as driver
    catalog = tmp_path / 'catalog'
    catalog.mkdir()
    (catalog / 'lc_cores-0.diffsky_gals.hdf5').touch()
    output = tmp_path / 'output'
    monkeypatch.setattr(driver, 'load_model', lambda *a: ({}, {'h': 0.7}, 3))
    monkeypatch.setattr(driver, 'check_patch_metadata', lambda *a: None)
    monkeypatch.setattr(driver, 'build_decam_curves', lambda: [
        SimpleNamespace(wave=np.array([4000., 5000.]), transmission=np.array([0., 1.]))
    ] * 5)
    original = driver.write_batches
    data = {'gal_id': np.array([1, 1 if duplicate else 2]), 'redshift_true': np.array([0.2, 0.5])}

    def write(*args):
        return original(*args, reader=lambda *a: data,
                        compute=lambda batch, **kw: {name: batch['redshift_true'] + 20 for name in COLUMN_NAMES})

    monkeypatch.setattr(driver, 'write_batches', write)
    argv = ['--catalog-dir', str(catalog), '--mock-version', 'fixture',
            '--output-dir', str(output), '--z-min', '0.1', '--z-max', '1.0']
    if duplicate:
        with pytest.raises(ValueError, match='Duplicate gal_id'):
            driver.main(argv)
        assert not output.exists()
        assert not list(tmp_path.glob('.decam-*'))
    else:
        driver.main(argv)
        import json
        assert json.loads((output / 'metadata.json').read_text())['rows'] == 2
        assert (output / 'filters_and_grid.npz').exists()
        assert pq.read_table(output / 'photometry.parquet').num_rows == 2
        with pytest.raises(FileExistsError):
            driver.main(argv)
