#!/usr/bin/env python3
"""Prepare a DES config and run Diffsky photometry followed by SNANA conversion."""
import argparse
from pathlib import Path
import shlex
import subprocess
import sys

import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog-dir', type=Path, required=True)
    parser.add_argument('--mock-version', required=True)
    parser.add_argument('--model-dir', type=Path)
    parser.add_argument('--run-dir', type=Path, required=True, help='New output directory')
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('snana_des.config'))
    parser.add_argument('--snana-dir', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--photometry-python', default=sys.executable)
    parser.add_argument('--snana-python', default=sys.executable)
    parser.add_argument('--scatter-policy', choices=['catalog', 'zero'], default='catalog')
    parser.add_argument('--prepare-only', action='store_true', help='Write config and print commands; do not run either stage')
    args = parser.parse_args()
    catalog = args.catalog_dir.resolve()
    model = (args.model_dir or catalog).resolve()
    run = args.run_dir.resolve()
    snana = args.snana_dir.resolve()
    if not catalog.is_dir() or not model.is_dir():
        parser.error('Catalog and model directories must exist')
    if any(c.isspace() for c in str(run)):
        parser.error('Use a run directory without whitespace (SNANA HOSTLIB_FILE syntax)')
    config = yaml.safe_load(args.config.read_text())
    cuts = [row.split() for row in config['CUTWIN'] if row.split()[0] == 'redshift_true']
    if len(cuts) != 1 or len(cuts[0]) != 3:
        parser.error('Config must contain one redshift_true CUTWIN with two bounds')
    z_min, z_max = map(float, cuts[0][1:])
    if not 0 < z_min < z_max < float('inf'):
        parser.error('Require finite 0 < redshift lower bound < upper bound')
    for filename in ['make_decam_photometry.py', 'make_hostlib_diffsky.py']:
        if not (snana / 'util' / filename).is_file():
            parser.error(f'Missing {snana / "util" / filename}')
    config['CAT_DIR'] = str(catalog)
    config['OVERRIDE_FILE'] = str(run / 'photometry' / 'photometry.parquet')
    config['HOSTLIB_FILE'] = [f'{run}/DES_DIFFSKY.HOSTLIB 1.0']
    if config.get('OVERRIDE_COLUMNS', '').split() != ['des_g', 'des_r', 'des_i', 'des_z']:
        parser.error('This example expects OVERRIDE_COLUMNS: des_g des_r des_i des_z')
    run.mkdir(parents=True, exist_ok=False)
    config_path = run / 'snana_des.config'
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    commands = [
        [args.photometry_python, str(snana / 'util/make_decam_photometry.py'),
         '--catalog-dir', str(catalog), '--model-dir', str(model),
         '--mock-version', args.mock_version, '--output-dir', str(run / 'photometry'),
         '--z-min', str(z_min), '--z-max', str(z_max), '--scatter-policy', args.scatter_policy],
        [args.snana_python, str(snana / 'util/make_hostlib_diffsky.py'), str(config_path)],
    ]
    print(f'Configuration: {config_path}', flush=True)
    for index, command in enumerate(commands, start=1):
        print(f'Stage {index}: {shlex.join(command)}', flush=True)
        if not args.prepare_only:
            log = run / f'stage{index}.log'
            print(f'Writing output to {log}', flush=True)
            with log.open('x') as stream:
                subprocess.run(command, cwd=run, stdout=stream, stderr=subprocess.STDOUT, check=True)
    if args.prepare_only:
        print('Prepared only. Run the printed commands from the run directory, or rerun this example with a NEW run directory and without --prepare-only.')
    else:
        print(f'Created {run}/DES_DIFFSKY.HOSTLIB; review stage2.log for missing matches and duplicate warnings.')


if __name__ == '__main__':
    main()
