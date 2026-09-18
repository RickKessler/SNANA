"""Exercise the existing SNANA parser and override join with a real Parquet file."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
import yaml

import make_hostlib_diffsky as snana


def test_example_config_and_snana_join(tmp_path):
    config_path = Path(__file__).resolve().parents[2] / 'doc/diffsky_decam/examples/snana_des.config'
    config = yaml.safe_load(config_path.read_text())
    ids = np.array([2**53 + 1, 2**53 + 2], dtype=np.int64)
    mags = {'serial_tag': ids, **{f'des_{band}': np.array([25., 24.]) for band in 'grizY'}}
    override = tmp_path / 'photometry.parquet'
    pd.DataFrame(mags).to_parquet(override, index=False)
    config['OVERRIDE_FILE'] = str(override)
    config['HOSTLIB_FILE'] = [f'{tmp_path}/unused.HOSTLIB 1.0']
    config = snana.parse_config_driver(config)
    assert config['band_magerr_list_diffsky'] == ['des_g', 'des_r', 'des_i', 'des_z']
    assert config['hostlib_varname_dict']['des_g'] == 'g_obs'
    # Reverse catalog order: the join must use exact IDs rather than row position.
    catalog = pd.DataFrame({'gal_id': ids[::-1], 'row_marker': ['second', 'first']})
    result = snana.inject_override_columns(catalog, config)
    assert result.row_marker.tolist() == ['second', 'first']
    np.testing.assert_array_equal(result.des_g, [24., 25.])
    assert 'gal_id' not in result
    result = snana.add_col_magerr_snr5(result, config)
    np.testing.assert_allclose(result.des_g_err, 1.086 / 5 * 10**(0.4 * (result.des_g - 25)))
