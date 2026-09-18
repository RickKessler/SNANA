# Diffsky DECam photometry for SNANA

Compute synthetic DECam **grizY** magnitudes from a Diffsky LastJourney catalog,
then supply them to SNANA's `make_hostlib_diffsky.py` to build a DES host library.
The calculation uses Diffsky's existing galaxy photometry engine with sedpy's
DECam throughput curves.

```text
Diffsky catalog + matching SSP templates and model parameters
    → make_decam_photometry.py
    → photometry.parquet
    → SNANA make_hostlib_diffsky.py + your configuration
    → final HOSTLIB
```

This guide describes the optional preprocessing utility and an example SNANA
configuration. Running the preprocessor does not modify the SNANA installation. The driver is adapted from
our LastJourney production workflow; its real-catalog photometry and complete
SNANA integration still require validation in the target environment.

## Requirements

Use an environment compatible with the **catalog release that you are reading**:

- Python 3.10 or newer and NumPy
- Diffsky, Diffmah, Diffstar, DSPS, and JAX
- Astronomy sedpy, providing `sedpy.observate` and the `decam_*` filters
- OpenCosmo, h5py, and hdf5plugin, including support for the catalog's compression
- PyArrow for the intermediate Parquet table
- pytest for the regression tests

The adapter requires
`diffsky.data_loaders.hacc_utils.load_lc_mock.compute_phot_from_diffsky_mock`.
The historical environment contains Diffsky 0.3.4, Diffmah 0.7.3, Diffstar 1.0.1,
OpenCosmo 1.2.4, and JAX/JAXlib 0.9.0. Its DSPS metadata reports 0.0.0, so these
version strings alone cannot reconstruct that environment. The newer environment
inspected during development no longer provides the required Diffsky API; the
driver rejects that mismatch explicitly.

The original production workflow used OpenCosmo 1.2.4. This is historical
context, not a tested version lock for this driver. Obtain compatible versions
and matching SSP/model files from the Diffsky maintainers; installing unrelated
latest versions is not a reproducibility specification. The driver records
installed package versions in its output metadata.

For the final conversion, your SNANA checkout must provide
`util/make_hostlib_diffsky.py` with **both `OVERRIDE_FILE` and `OVERRIDE_COLUMNS`
support**, reading `gal_id` from OpenCosmo and matching it to `serial_tag` in the
external table. This interface was inspected in our local SNANA checkout;
compatibility with every upstream SNANA version is not claimed.

## 1. Generate the photometry table

Set `SNANA_DIR` to your SNANA checkout, then run in your Diffsky environment:

```bash
python "$SNANA_DIR/util/make_decam_photometry.py" \
  --catalog-dir /path/to/diffsky/catalog_release \
  --mock-version catalog_release_name \
  --output-dir /path/to/decam_output \
  --z-min 0.05 \
  --z-max 1.8 \
  --batch-size 10000 \
  --grid-size 200
```

Replace every `/path/to/...` and `catalog_release_name` with your actual inputs.
`--mock-version` is the name expected by Diffsky's model loaders. SSP templates
and model parameters are loaded from `--catalog-dir` unless you supply
`--model-dir /path/to/matching/model_files`.

### Concrete input locations used on Perlmutter

Our existing production workflow reads the following catalog release. The SSP
and parameter loaders use this same directory as their model directory:

```bash
export SNANA_DIR="$HOME/SNANA"
export DECAM_CATALOG=/global/cfs/cdirs/hacc/OpenCosmo/LastJourney/synthetic_galaxies/hltds_cosmos_260215_02_17_2026
export DECAM_MODEL_DIR="$DECAM_CATALOG"
export DECAM_MOCK_VERSION=hltds_cosmos_260215_02_17_2026
```

These are NERSC filesystem locations, not public download URLs; access to the
shared HACC catalog is required. Outside NERSC, obtain the same release and
matching SSP/model files from the Diffsky/OpenCosmo maintainers and substitute
your local paths. Do not substitute a different release's model files.

The original local scripts and configuration are in:

```text
/global/cfs/cdirs/desc-sn/SNANA/SURVEYS/LSST/USERS/ayanmitr/HOSTLIB_DIFFSKY/DIFFSKY_MAG
```

The historical photometry Python executable is:

```bash
export DECAM_PHOT_PYTHON=/global/common/software/lsst/install/td_env/2026-04-07-37-02/py/envs/td_env/bin/python
```

This identifies the environment used by the earlier production scripts; the new
end-to-end workflow still needs real-catalog validation. Use your current
SNANA-compatible Python for the conversion stage; the runner accepts separate
Python executables for the two stages.

The supported input layout is the LastJourney layout used by the original
workflow: direct children named `lc_cores-*.diffsky_gals.hdf5`, cosmology attributes
under `header/simulation/cosmology`, and the SSP scatter dimension under
`cores/data/delta_mag_ssp_scatter`. Other catalog layouts need a reader adapter.
All patches must belong to the same release and matching model realization.

The driver:

- Reads resolved and synthetic cores through `opencosmo.open(..., synth_cores=True)`.
- Selects the inclusive redshift interval; applies no sky or mass cuts.
- Reads one patch into memory, then computes photometry in batches. Batch size
  controls photometry memory, **not** the size of the catalog read.
- Preserves the integer galaxy IDs and rejects duplicates across all selected
  rows, including duplicates across patches.
- Stops on read or computation errors. It does not skip failed patches.
- Publishes a new output directory only after successful completion. Existing
  output directories are never overwritten; a failed run must be restarted.

Run substantial catalogs on a compute node using your site's scheduler. The
repository does not submit jobs or assume a particular scheduler. The grid
size of 200 is a starting setting, not an established accuracy guarantee.

The output directory contains:

| File | Contents |
| --- | --- |
| `photometry.parquet` | `serial_tag` as signed int64, followed by `des_g`, `des_r`, `des_i`, `des_z`, `des_Y` |
| `metadata.json` | Arguments, versions, cosmology, input patch sizes/timestamps, row counts, and scatter policy |
| `filters_and_grid.npz` | Actual throughput arrays and photometry redshift grid |

`serial_tag` is copied from the catalog's `gal_id`; it is not an arbitrary row
number. This table is an **intermediate photometry table**, not a complete HOSTLIB.
Retain the matching model files and environment specification alongside these
outputs; metadata alone does not archive the model inputs.

### SSP scatter

The default `--scatter-policy catalog` requires the catalog's actual
`delta_mag_ssp_scatter` array, including synthetic-core values. If your
OpenCosmo/compression setup cannot expose that multidimensional column, the run
fails. Fix the reader/environment with the catalog maintainers before nominal
production.

For an explicitly approximate comparison with the old production workflow:

```bash
python "$SNANA_DIR/util/make_decam_photometry.py" \
  --catalog-dir /path/to/diffsky/catalog_release \
  --mock-version catalog_release_name \
  --output-dir /path/to/decam_zero_scatter \
  --z-min 0.05 --z-max 1.8 \
  --scatter-policy zero
```

Zero scatter is not guaranteed to equal ensemble-mean flux or magnitude. The
choice is recorded in `metadata.json`; there is no automatic fallback to zero.

## 2. Configure the SNANA HOSTLIB builder

Copy [examples/snana_des.config](examples/snana_des.config) into your own run
directory and edit its paths, sky/mass/redshift selection, output filename, and
survey depths. The example follows the local LastJourney HOSTLIB configuration
and includes disk/bulge morphology. Check that your release provides those fields.

The connection between the two steps is:

```yaml
CAT_DIR: /path/to/diffsky/catalog_release
OVERRIDE_FILE: /path/to/decam_output/photometry.parquet
OVERRIDE_COLUMNS: des_g des_r des_i des_z
```

The mapping in that same configuration includes:

```yaml
HOSTLIB_VARNAMES_MAP:
- serial_tag GALID 20.0f
- des_g g_obs 6.3f
- des_r r_obs 6.3f
- des_i i_obs 6.3f
- des_z z_obs 6.3f
```

This excerpt only illustrates the magnitude mapping; use the full example for
host properties. SNANA's output `GALID` may be generated by its own serial-number
logic. It is separate from the intermediate `serial_tag` used for the input join.

The driver computes all five bands; the example uses DES griz. To include Y,
add `des_Y` to `OVERRIDE_COLUMNS`, add `des_Y Y_obs 6.3f` to
`HOSTLIB_VARNAMES_MAP`, and, if using depth-based errors, add a `des_Y` depth and
`des_Y_err Y_obs_err 6.3f` mapping. Uppercase `Y` is intentional.

`MAGERR_SNR5` specifies the survey depths from which SNANA's utility derives
magnitude uncertainties. The example's values of 25 are placeholders. The
photometry adapter itself adds no measurement noise.

Use the **same catalog release and input patch set** in both steps. Keep the
SNANA redshift cuts inside the photometry interval. Start with cuts on native
catalog quantities; verify your SNANA version's handling before adding cuts on
injected DES magnitudes. Sky and mass cuts can be narrower in SNANA because the
photometry driver does not apply those cuts.

### Galaxy identity and missing matches

Never replace `gal_id` with `core_tag`: synthetic cores can share `core_tag=-1`.
The driver rejects duplicate selected `gal_id` values because the local SNANA
join otherwise retains the first occurrence. If this fails for light-cone
replicas, a composite-key change must be agreed with the catalog and SNANA
maintainers. Renumbering only the photometry table would break the join.

Uniqueness within this run does not establish identity across other releases or
patch selections. Review the SNANA conversion log and require **zero missing
photometry matches and zero duplicate-ID warnings** for the intended sample.
The local SNANA utility can warn and fill missing values rather than fail; this
repository does not change that behavior.

## 3. Run SNANA's conversion

After preparing your configuration, run this yourself in the SNANA environment:

```bash
python /path/to/SNANA/util/make_hostlib_diffsky.py /path/to/snana_des.config
```

The resulting HOSTLIB path comes from `HOSTLIB_FILE` in the configuration.
SNANA reads the other galaxy properties from Diffsky, injects the DES magnitudes,
and derives the configured host quantities and magnitude errors.

Running the preprocessor does not edit your SNANA checkout, existing production
configuration, or previous HOSTLIBs. The command above is an instruction for a later run.

## Runnable example: configuration plus both stages

[examples/run_decam_example.py](examples/run_decam_example.py) reads the example
configuration, fills in the input/output paths, and runs the two commands in
order. It uses the configuration's redshift bounds for photometry, retains its
sky/mass cuts and magnitude depths, and writes one final HOSTLIB. You can supply
an edited configuration with `--config /path/to/my_snana_des.config`.

After setting the Perlmutter variables above, first prepare and inspect the run:

```bash
python "$SNANA_DIR/doc/diffsky_decam/examples/run_decam_example.py" \
  --catalog-dir "$DECAM_CATALOG" \
  --model-dir "$DECAM_MODEL_DIR" \
  --mock-version "$DECAM_MOCK_VERSION" \
  --photometry-python "$DECAM_PHOT_PYTHON" \
  --snana-python "$(command -v python)" \
  --run-dir "$PWD/decam_prepare_example" \
  --prepare-only
```

This writes `decam_prepare_example/snana_des.config` and prints the exact
commands without reading galaxy data or running either stage. The Python used
to launch the runner needs PyYAML. Paths and the selected config are checked.

To execute both stages **on a compute node**, use a new output directory:

```bash
python "$SNANA_DIR/doc/diffsky_decam/examples/run_decam_example.py" \
  --catalog-dir "$DECAM_CATALOG" \
  --model-dir "$DECAM_MODEL_DIR" \
  --mock-version "$DECAM_MOCK_VERSION" \
  --photometry-python "$DECAM_PHOT_PYTHON" \
  --snana-python "$(command -v python)" \
  --run-dir "$PWD/decam_run_example"
```

This is a workflow example, not a cheap two-galaxy smoke test: it processes all
patches in the supplied catalog directory over the config's redshift interval.
Use a small compatible catalog fixture for a short scientific test. The runner
does not allocate compute resources or submit a job.

The resulting directory contains:

```text
decam_run_example/
  snana_des.config              # resolved configuration passed to SNANA
  photometry/                  # Parquet magnitudes, metadata, filters and grid
  stage1.log                   # photometry output
  stage2.log                   # SNANA conversion output
  DES_DIFFSKY.HOSTLIB           # final host library
```

A failed first stage prevents the SNANA stage from starting. Existing run
directories are not overwritten. Review the logs before using the HOSTLIB;
subprocess success alone does not establish scientific correctness or complete
photometry matching. For an explicitly approximate reproduction of the old
zero-scatter setting, add `--scatter-policy zero`; the default requires catalog
scatter. Replace the example magnitude depths before production.

## Validation before production

Run the local regression tests:

```bash
python -m pytest -q "$SNANA_DIR/util/tests/test_des_photometry.py" \
  "$SNANA_DIR/util/tests/test_make_decam_photometry.py" \
  "$SNANA_DIR/util/tests/test_decam_snana_interface.py"
```

The tests cover adapter behavior, batching, scatter/row alignment, exact integer
IDs above 2^53, duplicate rejection, actual Parquet serialization, successful
output publication, and cleanup after a failed run. An interface test also uses
SNANA's actual configuration parser, override join, and depth-based error calculation. They use
a substitute photometry engine and do **not** validate the astrophysics.

Before a full production run:

1. Use a small catalog fixture containing resolved and synthetic cores, with
   matching model files. Run the driver and compare against a direct Diffsky
   call with identical inputs and filters.
2. Reproduce native catalog photometry with its native filters to check the
   model/scatter setup. Compare DECam results to a trusted independent integration.
3. Refine the redshift grid and agree on acceptable magnitude differences,
   including low-redshift, dusty, quenched, and high-redshift galaxies.
4. Convert that same sample with SNANA, check match counts and magnitude columns,
   and run a minimal SNANA simulation that reads the final HOSTLIB.

Confirm DES throughput/atmosphere definitions, observed-frame AB conventions,
dust, lensing, and Milky Way extinction with the maintainers. The adapter inherits
Diffsky's conventions and applies no additional extinction or magnification.
Morphology and band-dependent component fluxes require separate validation for
host surface-brightness applications.

## Using the adapter in another pipeline

If your Diffsky pipeline already prepares the model objects and row-aligned
arrays, use the adapter directly:

```python
# Put SNANA/util on your Python import path.
from des_photometry import build_decam_curves, compute_decam_photometry

curves = build_decam_curves()  # reuse across batches
columns = compute_decam_photometry(
    diffsky_data,
    ssp_data=ssp_data,
    param_collection=param_collection,
    sim_info=sim_info,
    z_phot_table=z_phot_table,
    tcurves=curves,
)
# Attach columns to the same galaxies in the same order.
```

The caller supplies the MAH/SFH and stochastic parameters required by its
Diffsky engine. Empty batches must be skipped. Native DES catalog columns could
ultimately remove the external override step entirely.

See [DISCUSSION_QUESTIONS.md](DISCUSSION_QUESTIONS.md) for the outstanding
scientific and upstream-integration questions.
