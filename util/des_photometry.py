"""DECam photometry adapter for Diffsky's existing galaxy photometry engine.

No catalog paths, selection cuts, geometry prescriptions, or SNANA dependency.
The caller supplies row-aligned Diffsky arrays and matching model objects.
"""
from collections import namedtuple
import numpy as np

FILTER_NAMES = ("decam_g", "decam_r", "decam_i", "decam_z", "decam_Y")
COLUMN_NAMES = ("des_g", "des_r", "des_i", "des_z", "des_Y")
DECamCurves = namedtuple("DECamCurves", COLUMN_NAMES)


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


def compute_decam_photometry(
    diffsky_data, *, ssp_data, param_collection, sim_info, z_phot_table,
    tcurves=None, scatter_policy="catalog", n_ssp_cols=None,
):
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
