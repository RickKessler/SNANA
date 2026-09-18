"""Regression tests for diffsky_synthetic_photometry.py: config parsing, the
generic filter interface, and the batch-writing driver. No real catalog
access or JAX computation -- the Diffsky engine and filter curves are mocked.
"""

import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from types import ModuleType, SimpleNamespace
import sqlite3

import numpy as np
import pyarrow.parquet as pq
import pytest
import yaml

import diffsky_synthetic_photometry as driver
from diffsky_synthetic_photometry import (
    build_run_args, compute_synthetic_photometry, load_cutwin, load_filters,
    parse_args, read_yaml, register_ids, select_rows, write_batches,
)


COLUMN_NAMES = ('band_a', 'band_b', 'band_c', 'band_d', 'band_e')


# ---- config / filter parsing --------------------------------------------

def write_config(tmp_path, name='config.yml', **overrides):
    base = {
        'CATALOG_DIR': str(tmp_path / 'catalog'),
        'OUTPUT_DIR': str(tmp_path / 'output'),
        'Z_MIN': 0.1, 'Z_MAX': 1.0,
        'FILTERS': [{'SEDPY_ID': 'survey_g', 'COLUMN': 'band_g'},
                    {'SEDPY_ID': 'survey_r', 'COLUMN': 'band_r'}],
    }
    base.update(overrides)
    path = tmp_path / name
    path.write_text(yaml.safe_dump(base))
    return path


def test_load_filters_inline():
    config = {'FILTERS': [{'SEDPY_ID': 'survey_g', 'COLUMN': 'band_g'},
                          {'SEDPY_ID': 'survey_r', 'COLUMN': 'band_r'}]}
    sedpy_ids, column_names = load_filters(config)
    assert sedpy_ids == ('survey_g', 'survey_r')
    assert column_names == ('band_g', 'band_r')


def test_load_filters_from_file(tmp_path):
    filter_file = tmp_path / 'filters.yml'
    filter_file.write_text('FILTERS:\n- SEDPY_ID: survey_i\n  COLUMN: band_i\n')
    sedpy_ids, column_names = load_filters({'FILTER_FILE': str(filter_file)})
    assert sedpy_ids == ('survey_i',)
    assert column_names == ('band_i',)


def test_load_filters_requires_one_source():
    with pytest.raises(ValueError, match='FILTERS.*FILTER_FILE'):
        load_filters({})


def test_load_filters_rejects_duplicate_columns():
    config = {'FILTERS': [{'SEDPY_ID': 'a', 'COLUMN': 'x'}, {'SEDPY_ID': 'b', 'COLUMN': 'x'}]}
    with pytest.raises(ValueError, match='unique'):
        load_filters(config)


def test_build_run_args_defaults_and_overrides(tmp_path):
    config = read_yaml(write_config(tmp_path))
    args = build_run_args(config, wildcard='007')
    assert args.catalog_dir == Path(config['CATALOG_DIR'])
    assert args.model_dir is None
    assert args.mock_version is None
    assert args.grid_size == 200 and args.batch_size == 10000
    assert args.scatter_policy == 'catalog'
    assert args.wildcard == '007'
    assert args.sedpy_ids == ('survey_g', 'survey_r')
    assert args.column_names == ('band_g', 'band_r')
    assert args.cutwin == ()


def test_build_run_args_wires_cutwin(tmp_path):
    config = read_yaml(write_config(tmp_path, CUTWIN=[{'COLUMN': 'ra', 'MIN': 240.0, 'MAX': 245.0}]))
    args = build_run_args(config, wildcard=None)
    assert args.cutwin == (('ra', 240.0, 245.0),)


@pytest.mark.parametrize('missing_key', ['CATALOG_DIR', 'OUTPUT_DIR', 'Z_MIN', 'Z_MAX'])
def test_build_run_args_missing_required_key(tmp_path, missing_key):
    config = read_yaml(write_config(tmp_path))
    del config[missing_key]
    with pytest.raises(ValueError, match='missing required key'):
        build_run_args(config, wildcard=None)


@pytest.mark.parametrize('overrides', [
    {'Z_MIN': 1.0, 'Z_MAX': 0.5},
    {'GRID_SIZE': 1},
    {'BATCH_SIZE': 0},
    {'SCATTER_POLICY': 'bogus'},
])
def test_build_run_args_invalid_values(tmp_path, overrides):
    config = read_yaml(write_config(tmp_path, **overrides))
    with pytest.raises(ValueError):
        build_run_args(config, wildcard=None)


def test_parse_args_requires_config_file():
    with pytest.raises(SystemExit):
        parse_args([])


def test_parse_args_wildcard_flag():
    args = parse_args(['some_config.yml', '-w', '14'])
    assert args.config_file == Path('some_config.yml')
    assert args.wildcard == '14'


# ---- compute_synthetic_photometry (engine interface) ---------------------

@pytest.fixture
def engine(monkeypatch):
    calls = []

    def compute(**kwargs):
        calls.append(kwargs)
        return {"obs_mags": np.arange(10).reshape(2, 5) + 20.0}

    package = ModuleType("diffsky.data_loaders.hacc_utils")
    package.load_lc_mock = SimpleNamespace(compute_phot_from_diffsky_mock=compute)
    monkeypatch.setitem(sys.modules, "diffsky.data_loaders.hacc_utils", package)
    return calls


def inputs():
    return {
        "redshift_true": np.array([0.8, 0.2]),
        "delta_mag_ssp_scatter": np.ones((2, 3)),
    }


def run(data, **kwargs):
    return compute_synthetic_photometry(
        data, ssp_data=None, param_collection=None, sim_info=None,
        z_phot_table=np.array([0.1, 0.9]), tcurves="preloaded",
        column_names=COLUMN_NAMES, **kwargs,
    )


def test_column_and_row_order_and_scatter_preserved(engine):
    data = inputs()
    result = run(data)
    assert tuple(result) == COLUMN_NAMES
    np.testing.assert_array_equal(result["band_a"], [20, 25])
    np.testing.assert_array_equal(result["band_e"], [24, 29])
    assert engine[0]["diffsky_data"]["delta_mag_ssp_scatter"] is data["delta_mag_ssp_scatter"]


def test_explicit_zero_does_not_mutate_input(engine):
    data = inputs()
    run(data, scatter_policy="zero", n_ssp_cols=3)
    assert np.all(data["delta_mag_ssp_scatter"] == 1)
    assert np.all(engine[0]["diffsky_data"]["delta_mag_ssp_scatter"] == 0)


def test_missing_scatter_is_not_silently_replaced(engine):
    with pytest.raises(ValueError, match="catalog scatter is required"):
        run({"redshift_true": np.array([0.8, 0.2])})
    assert not engine


def test_nonfinite_photometry_rejected(engine, monkeypatch):
    package = sys.modules["diffsky.data_loaders.hacc_utils"]
    monkeypatch.setattr(package.load_lc_mock, "compute_phot_from_diffsky_mock",
                        lambda **kw: {"obs_mags": np.full((2, 5), np.nan)})
    with pytest.raises(ValueError, match="invalid magnitudes"):
        run(inputs())


def test_wrong_column_count_rejected(engine):
    with pytest.raises(ValueError, match=r"N_gal, 2"):
        compute_synthetic_photometry(
            inputs(), ssp_data=None, param_collection=None, sim_info=None,
            z_phot_table=np.array([0.1, 0.9]), tcurves="preloaded",
            column_names=("band_a", "band_b"),
        )


# ---- write_batches (unchanged logic, now column_names is explicit) ------

def test_integer_ids_are_exact_and_duplicates_rejected():
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE ids (id INTEGER PRIMARY KEY)')
        ids = np.array([2**53 + 1, 2**53 + 2], dtype=np.int64)
        np.testing.assert_array_equal(register_ids(db, ids), ids)
        with pytest.raises(ValueError, match='Duplicate gal_id'):
            register_ids(db, ids[:1])


def test_batches_keep_scatter_and_large_ids_aligned(tmp_path):
    patch = tmp_path / 'patch.hdf5'
    patch.touch()
    args = SimpleNamespace(z_min=0.1, z_max=1.0, batch_size=2, scatter_policy='catalog', cutwin=())
    data = {'gal_id': np.array([2**53+1, 2**53+2, 2**53+3, 2**53+4], dtype=np.int64),
            'redshift_true': np.array([0.5, 0.01, 0.8, 1.0]),
            'delta_mag_ssp_scatter': np.arange(12).reshape(4, 3)}
    seen = []

    def compute(batch, **kwargs):
        seen.append(batch['delta_mag_ssp_scatter'].copy())
        return {name: batch['redshift_true'] + i + 20 for i, name in enumerate(COLUMN_NAMES)}

    counts = write_batches([patch], args, tmp_path, {}, None, None, 3, COLUMN_NAMES,
                          reader=lambda *a: data, compute=compute)
    out = pq.read_table(tmp_path / 'photometry.parquet').to_pandas()
    assert counts[0]['selected_rows'] == 3
    np.testing.assert_array_equal(out.serial_tag, data['gal_id'][[0, 2, 3]])
    np.testing.assert_array_equal(np.concatenate(seen), data['delta_mag_ssp_scatter'][[0, 2, 3]])
    np.testing.assert_allclose(out.band_a, [20.5, 20.8, 21])
    assert list(out) == ['serial_tag', *COLUMN_NAMES]


def test_empty_selection_fails(tmp_path):
    path = tmp_path / 'empty'
    path.touch()
    args = SimpleNamespace(z_min=0.1, z_max=1, batch_size=1, scatter_policy='catalog', cutwin=())
    with pytest.raises(ValueError, match='No galaxies selected'):
        write_batches([path], args, tmp_path, {}, None, None, 3, COLUMN_NAMES,
                      reader=lambda *a: {'gal_id': np.array([7]), 'redshift_true': np.array([2.])})


# ---- CUTWIN pre-selection (speed: avoid photometering discarded galaxies) ----

def test_load_cutwin_empty_by_default():
    assert load_cutwin({}) == ()


def test_load_cutwin_parses_min_and_max():
    config = {'CUTWIN': [{'COLUMN': 'ra', 'MIN': 240.0, 'MAX': 245.0},
                         {'COLUMN': 'logsm_obs', 'MIN': 8.0}]}
    cuts = load_cutwin(config)
    assert cuts == (('ra', 240.0, 245.0), ('logsm_obs', 8.0, None))


def test_load_cutwin_requires_min_or_max():
    with pytest.raises(ValueError, match='MIN and/or MAX'):
        load_cutwin({'CUTWIN': [{'COLUMN': 'ra'}]})


def test_select_rows_combines_zrange_and_cutwin():
    data = {'redshift_true': np.array([0.5, 0.5, 0.5, 1.5]),
            'ra': np.array([241.0, 250.0, 242.0, 241.0]),
            'logsm_obs': np.array([9.0, 9.0, 7.0, 9.0])}
    cutwin = (('ra', 240.0, 245.0), ('logsm_obs', 8.0, None))
    np.testing.assert_array_equal(select_rows(data, 0.1, 1.0, cutwin), [0])


def test_write_batches_applies_cutwin_and_drops_cut_columns_before_compute(tmp_path):
    patch = tmp_path / 'patch.hdf5'
    patch.touch()
    args = SimpleNamespace(z_min=0.1, z_max=1.0, batch_size=10, scatter_policy='catalog',
                          cutwin=(('ra', 240.0, 245.0),))
    data = {'gal_id': np.array([1, 2, 3]),
            'redshift_true': np.array([0.5, 0.5, 0.5]),
            'ra': np.array([241.0, 250.0, 242.0]),
            'delta_mag_ssp_scatter': np.arange(9).reshape(3, 3)}
    seen_batches = []

    def compute(batch, **kwargs):
        seen_batches.append(batch)
        return {name: np.full(len(batch['redshift_true']), 20.0) for name in COLUMN_NAMES}

    counts = write_batches([patch], args, tmp_path, {}, None, None, 3, COLUMN_NAMES,
                          reader=lambda *a: dict(data), compute=compute)
    assert counts[0]['selected_rows'] == 2
    assert 'ra' not in seen_batches[0]


# ---- main(): full config-driven run, engine/curves mocked ---------------

def write_full_config(tmp_path, catalog, output, **overrides):
    config = {
        'CATALOG_DIR': str(catalog), 'OUTPUT_DIR': str(output),
        'Z_MIN': 0.1, 'Z_MAX': 1.0,
        'FILTERS': [{'SEDPY_ID': f'survey_{c}', 'COLUMN': c} for c in COLUMN_NAMES],
    }
    config.update(overrides)
    path = tmp_path / 'config.yml'
    path.write_text(yaml.safe_dump(config))
    return path


@pytest.mark.parametrize('duplicate', [False, True])
def test_main_publishes_complete_output_only(tmp_path, monkeypatch, duplicate):
    catalog = tmp_path / 'catalog'
    catalog.mkdir()
    (catalog / 'lc_cores-0.diffsky_gals.hdf5').touch()
    output = tmp_path / 'output'
    monkeypatch.setattr(driver, 'load_model', lambda *a: ({}, {'h': 0.7}, 3))
    monkeypatch.setattr(driver, 'check_patch_metadata', lambda *a: None)
    monkeypatch.setattr(driver, 'build_photometry_curves', lambda *a: [
        SimpleNamespace(wave=np.array([4000., 5000.]), transmission=np.array([0., 1.]))
    ] * len(COLUMN_NAMES))
    original = driver.write_batches
    data = {'gal_id': np.array([1, 1 if duplicate else 2]), 'redshift_true': np.array([0.2, 0.5])}

    def write(*args):
        return original(*args, reader=lambda *a: data,
                        compute=lambda batch, **kw: {name: batch['redshift_true'] + 20 for name in COLUMN_NAMES})

    monkeypatch.setattr(driver, 'write_batches', write)
    config_path = write_full_config(tmp_path, catalog, output)
    argv = [str(config_path)]
    if duplicate:
        with pytest.raises(ValueError, match='Duplicate gal_id'):
            driver.main(argv)
        assert not output.exists()
        assert not list(tmp_path.glob('.photometry-*'))
    else:
        driver.main(argv)
        assert json.loads((output / 'metadata.json').read_text())['rows'] == 2
        assert (output / 'filters_and_grid.npz').exists()
        assert pq.read_table(output / 'photometry.parquet').num_rows == 2
        with pytest.raises(FileExistsError):
            driver.main(argv)


def test_main_mock_version_defaults_to_catalog_dir_basename(tmp_path, monkeypatch):
    catalog = tmp_path / 'my_release_v1'
    catalog.mkdir()
    (catalog / 'lc_cores-0.diffsky_gals.hdf5').touch()
    output = tmp_path / 'output'
    seen_mock_version = []
    monkeypatch.setattr(driver, 'load_model', lambda args, patch: (
        seen_mock_version.append(args.mock_version), ({}, {'h': 0.7}, 3))[1])
    monkeypatch.setattr(driver, 'check_patch_metadata', lambda *a: None)
    monkeypatch.setattr(driver, 'build_photometry_curves', lambda *a: [
        SimpleNamespace(wave=np.array([4000., 5000.]), transmission=np.array([0., 1.]))
    ] * len(COLUMN_NAMES))
    monkeypatch.setattr(driver, 'write_batches', lambda *a, **kw: (_ for _ in ()).throw(
        FileNotFoundError('stop after mock_version check')))
    config_path = write_full_config(tmp_path, catalog, output)
    with pytest.raises(FileNotFoundError, match='stop after mock_version check'):
        driver.main([str(config_path)])
    assert seen_mock_version == ['my_release_v1']


def test_main_wildcard_selects_subset(tmp_path, monkeypatch):
    catalog = tmp_path / 'catalog'
    catalog.mkdir()
    (catalog / 'lc_cores-014.diffsky_gals.hdf5').touch()
    (catalog / 'lc_cores-999.diffsky_gals.hdf5').touch()
    output = tmp_path / 'output'
    seen_patches = []
    monkeypatch.setattr(driver, 'load_model', lambda args, patch: (
        seen_patches.append(patch.name), ({}, {'h': 0.7}, 3))[1])
    monkeypatch.setattr(driver, 'check_patch_metadata', lambda *a: None)
    monkeypatch.setattr(driver, 'build_photometry_curves', lambda *a: [
        SimpleNamespace(wave=np.array([4000., 5000.]), transmission=np.array([0., 1.]))
    ] * len(COLUMN_NAMES))
    monkeypatch.setattr(driver, 'write_batches', lambda *a, **kw: (_ for _ in ()).throw(
        FileNotFoundError('stop after patch selection')))
    config_path = write_full_config(tmp_path, catalog, output)
    with pytest.raises(FileNotFoundError, match='stop after patch selection'):
        driver.main([str(config_path), '-w', '014'])
    assert seen_patches == ['lc_cores-014.diffsky_gals.hdf5']
