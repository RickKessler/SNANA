#!/usr/bin/env python
import sys
import os
import argparse
import numpy as np


def make_simblib_ideal(survey='LSST', filters='ugrizY', zpt=35.,\
                    mjd_start=53000., mjd_end=53400., mjd_step=2.,\
                    simlib_file='IDEAL.SIMLIB', clobber=False):

    mjds = np.arange(mjd_start, mjd_end+mjd_step, mjd_step)
    nobs = len(mjds)*len(filters)

    with open(simlib_file, 'w') as outf:
        outf.write(f'SURVEY: {survey}   FILTERS: {filters} \n')
        outf.write('BEGIN LIBGEN\n\n')
        outf.write('LIBID: 1\n')
        outf.write(f'RA: 0.0   DECL: 0.0   NOBS: {nobs:d}   MWEBV: 0.0   PIXSIZE: 0.2\n\n')
        outf.write('#                           CCD  CCD         PSF1 PSF2 PSF2/1 \n')
        outf.write('#     MJD      IDEXPT  FLT GAIN NOISE SKYSIG (pixels)  RATIO  ZPTAVG ZPTERR  MAG \n')

        idexpt = 10000 +1
        for mjd in mjds:
            for pb in filters:
                outf.write(f'S: {mjd:.3f}    {idexpt:d}   {pb}   1.00  1.00  40.0  1.00 0.00 0.000  35.00  0.00  99\n')
                idexpt +=1
            outf.write('\n')
        outf.write('END_LIBID: 1 \n')
        outf.write('END_OF_SIMLIB: \n')
    print(f'Wrote {nobs} exposures, check: {idexpt-10000-1}')
    return


def get_options(argv):
    parser = argparse.ArgumentParser(description='Generate a SNANA IDEAL SIMLIB')
    parser.add_argument('--survey', '-s', default='LSST',\
                        help='Survey Name')
    parser.add_argument('--filters', '-f', default='ugrizY',\
                        help='Filters (one letter per, no delimiter)')
    parser.add_argument('--zpt', '-z', default=35., type=float,\
                        help='Specify the survey zeropoint')
    parser.add_argument('--mjd_start','-x', default=53000., type=float,\
                        help='Specify MJD for SIMLIB start')
    parser.add_argument('--mjd_end','-y', default=53400., type=float,\
                        help='Specify MJD for SIMLIB end')
    parser.add_argument('--mjd_step','-d', default=2., type=float,\
                        help='Specify MJD step for SIMLIB')
    parser.add_argument('--simlib_file', '-o', default='IDEAL.SIMLIB',\
                        help='Specify the output file for the SIMLIB')
    parser.add_argument('--clobber', action='store_true',\
                        help='Clobber the output if it exists')
    args = parser.parse_args(args=argv)

    if args.mjd_start >= args.mjd_end:
        message = f'mjd_start ({args.mjd_start}) must be greater than mjd_end ({args.mjd_end})'
        raise ValueError(message)

    if args.zpt <= 0 or args.mjd_step <= 0:
        message = f'zpt ({args.zpt}) and mjd_step ({args.mjd_step}) must be > 0'
        raise ValueError(message)

    if os.path.isfile(args.simlib_file) and not args.clobber:
        message = f'SIMLIB file {args.simlib_file} exists. Not overwriting without --clobber'
        raise OSError(message)

    return(args)


def main(argv=None):
    args = get_options(argv=argv)
    survey      = args.survey
    filters     = args.filters
    zpt         = args.zpt
    mjd_start   = args.mjd_start
    mjd_end     = args.mjd_end
    mjd_step    = args.mjd_step
    simlib_file = args.simlib_file
    clobber     = args.clobber

    return make_simblib_ideal(survey=survey, filters=filters, zpt=zpt,\
                    mjd_start=mjd_start, mjd_end=mjd_end, mjd_step=mjd_step,\
                    simlib_file=simlib_file, clobber=clobber)


if __name__=='__main__':
    sys.exit(main(argv=sys.argv[1:]))
