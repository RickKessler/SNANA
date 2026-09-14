"""Interface regression tests; no catalog access or JAX computation."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from des_photometry import COLUMN_NAMES, compute_decam_photometry


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
    return compute_decam_photometry(
        data, ssp_data=None, param_collection=None, sim_info=None,
        z_phot_table=np.array([0.1, 0.9]), tcurves="preloaded", **kwargs,
    )


def test_column_and_row_order_and_scatter_preserved(engine):
    data = inputs()
    result = run(data)
    assert tuple(result) == COLUMN_NAMES
    np.testing.assert_array_equal(result["des_g"], [20, 25])
    np.testing.assert_array_equal(result["des_Y"], [24, 29])
    assert engine[0]["diffsky_data"]["delta_mag_ssp_scatter"] is data["delta_mag_ssp_scatter"]
    np.testing.assert_array_equal(engine[0]["diffsky_data"]["redshift_true"], [0.8, 0.2])


def test_explicit_zero_does_not_mutate_input(engine):
    data = inputs()
    run(data, scatter_policy="zero", n_ssp_cols=3)
    assert np.all(data["delta_mag_ssp_scatter"] == 1)
    assert np.all(engine[0]["diffsky_data"]["delta_mag_ssp_scatter"] == 0)


def test_missing_scatter_is_not_silently_replaced(engine):
    with pytest.raises(ValueError, match="catalog scatter is required"):
        run({"redshift_true": np.array([0.8, 0.2])})
    assert not engine


@pytest.mark.parametrize("redshift", [[0.0, 0.2], [np.nan, 0.2], [0.95, 0.2]])
def test_invalid_redshift_rejected_before_engine(engine, redshift):
    data = inputs()
    data["redshift_true"] = redshift
    with pytest.raises(ValueError):
        run(data)
    assert not engine


def test_nonfinite_photometry_rejected(engine, monkeypatch):
    package = sys.modules["diffsky.data_loaders.hacc_utils"]
    monkeypatch.setattr(package.load_lc_mock, "compute_phot_from_diffsky_mock",
                        lambda **kw: {"obs_mags": np.full((2, 5), np.nan)})
    with pytest.raises(ValueError, match="invalid magnitudes"):
        run(inputs())


def test_incompatible_diffsky_reports_required_api(engine, monkeypatch):
    package = sys.modules['diffsky.data_loaders.hacc_utils']
    monkeypatch.delattr(package.load_lc_mock, 'compute_phot_from_diffsky_mock')
    with pytest.raises(RuntimeError, match='environment matching the catalog'):
        run(inputs())
