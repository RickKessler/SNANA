#!/usr/bin/env python
#
# Created Mar 19 2021 by R.Kessler
#
# Mini-pipeline to create fluxError corrections for data and sim.
# The following stages are run here:
#  + create SIMLIB from fakes overlaid on images
#  + run simulation using SIMLIB to have same epochs and mags as fakes
#  + create tables with every observation
#  + make fluxError maps
#
# Aug 11 2021: if SBMAG-dependence is set, keep epochs with SBMAG<50
#             (initial use is LSST-DC2)
# Aug 12 2021: if too few fakes to compute cor, take value from 
#              nearest bin.
#
# Nov 28 2022: fix bug computing nrow_per_filter in make_fluxerr_model_map()
#               (divide by NFIELD_GROUP)
#              At stage07, rename DESfakes to {survey}_fakes
#
# ========================

import os, sys, argparse, glob, yaml, math
import numpy as np
from   argparse import Namespace
import pandas as pd

#JOBNAME_SNANA = "/home/rkessler/SNANA/bin/snana.exe"
JOBNAME_SNANA = "snana.exe"
JOBNAME_SIM   = "snlc_sim.exe"

TABLE_NAME    = "OUTLIER"

STRING_FAKE   = "FAKE"
STRING_SIM    = "SIM"
STRING_FIELDS = "FIELDS"

FLUXERRMODEL_FILENAME_FAKE = f"FLUXERRMODEL_{STRING_FAKE}.DAT"
FLUXERRMODEL_FILENAME_SIM  = f"FLUXERRMODEL_{STRING_SIM}.DAT"
MAP_NAME = "FLUXERR_SCALE"
 
FORCE_ERRCORR1 = True  # force Error corr>=1 to avoid sim nan

USERNAME      = os.environ['USER']
USERNAME4     = os.environ['USER'][0:4]
HOSTNAME      = os.uname()[1].split('.')[0]

COLNAME_BIN1D    = "BIN1D"
COLNAME_IFILTOBS = "IFILTOBS"
COLNAME_BAND     = "BAND"
COLNAME_IFIELD   = "IFIELD"

IFILTOBS_MAX = 80

ISTAGE_MAKEMAP = 5

# list of reduced flux correlations to try in sim 
REDCOV_LIST = [ 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 ]

NMLKEY_DATA_PATH   = 'PRIVATE_DATA_PATH'
NMLKEY_VERSION     = 'VERSION_PHOTOMETRY'
NMLKEY_KCOR_FILE   = 'KCOR_FILE'
NMLKEY_SNTABLE     = 'SNTABLE_LIST'
NMLKEY_SIMLIB_FILE = 'SIMLIB_OUTFILE'
NMLKEY_OPT_SIMLIB  = 'OPT_SIMLIB_OUT'
NMLKEY_LPROB_TRUEFLUX    = 'LPROB_TRUEFLUX'
NMLKEY_MODEL_FILE        = 'FLUXERRMODEL_FILE'
NMLKEY_SIM_MODEL_FILE    = 'SIM_FLUXERRMODEL_FILE' # for redcov_test
NMLKEY_TEXTFILE_PREFIX   = 'TEXTFILE_PREFIX'
NMLKEY_HFILE_OUT         = 'HFILE_OUT'

NMLKEY_LIST = [NMLKEY_DATA_PATH, NMLKEY_VERSION, NMLKEY_KCOR_FILE, 
               NMLKEY_SNTABLE, NMLKEY_TEXTFILE_PREFIX, NMLKEY_HFILE_OUT,
               NMLKEY_SIMLIB_FILE, NMLKEY_OPT_SIMLIB, NMLKEY_LPROB_TRUEFLUX,
               NMLKEY_MODEL_FILE, NMLKEY_SIM_MODEL_FILE ]

KEY_ALLBANDS        = "ALL"
VARNAME_PROB_PREFIX = "PROB_TRUEFLUX"  # varname in SNANA table

TABLE_SUFFIX_SNANA   = "SNANA.TEXT"
TABLE_SUFFIX_OUTLIER = "OUTLIER.TEXT"

REDCOV_SUMMARY_FILE = "REDCOV.SUMMARY"

# - - - - - - - - - - - - - 
HELP_CONFIG = """
# keys for input config file

OUTDIR: [OUTDIR]   # name of output directory

VERSION_FAKES:    [full path to data folder with fakes that include true mag]
HOSTLIB_FILE:     [full path to HOSTLIB]
KCOR_FILE:        [full path to KCOR/calib file]
PATH_SNDATA_SIM:  [optional path for sim data output]

OPT_SIMLIB:  2  # 1=all epochs, 2=only epochs with Ftrue>0

# optional, reject extreme NSIG outliers
CUTWIN_NSIG:  0.0  [NSIG_MAX]

# additional/optional &SNLCINP arguments to select events/epochs
EXTRA_SNLCINP_ARGS:
  - PHOTFLAG_MSKREJ = [MSK]  # reject events with these PHOTFLAG bits
  - OPT_SETPKMJD    = -9     # don't wast time computing PEAKMJD
  - CUTWIN_MJD      = [MJDMIN], [MJDMAX]  # select season(s)
  - CUTWIN_NFIELD   = 1, 1                # reject field overlaps
  - CUTWIN_ERRTEST  = 0.5, 2.0            # reject crazy ERR_CALC/ERR_DATA
  - SIMVAR_CUTWIN_STRING = 'NEP_SIM_MAGOBS 4 9999' # at least 4 Ftrue>0

#  - MXEVT_PROCESS = 500  # quick test with small subset of fakes

# Optional map-computation in independent groups of fields. 
# Here the maps are computed for SHALLOW and DEEP field groups.
FIELDS: 
  SHALLOW:  S1 S2 C1 C2 X1 X2 E1 E2
  DEEP:     X3 C3

# Define multi-dimensional map bins using variable(s) from OUTLIER table.
# Values outside map-bin range are pulled into first/last map bin so 
# that all obs are used.

FLUXERRMAP_BINS: 
  - FILTER                 # auto compute bins from filters in data file
  - SBMAG   8  20   28     # nbin min max (histogram bins)
  - PSF     3  1.0  4.0    # idem



# ==========================================================
# what to do with the output (argument of OUTDIR in config file)
#  [see snana manual section 4.14.1 FLUXERRMODEL Tables]

These are the three output files to use/examine:
   FLUXERRMODEL_SIM.DAT  FLUXERRMODEL_FAKE.DAT  REDCOV.SUMMARY

1. In the simulation, add input
   FLUXERRMODEL_FILE:  $PATH/FLUXERRMODEL_SIM.DAT

2. If real data errors have not been inflated for SBMAG, 
   add this input
   &SNLCINP
      FLUXERRMODEL_FILE = '$PATH/FLUXERRMODEL_FAKE.DAT'

  Beware that LCFIT implements FLUXERRMODEL_FILE only for real data;
  it is ignored for sims because sim errors should already be inflated.
  This auto-ignore feature allow using the same LCFIT input file for
  data and sim.

3. The final step is to select reduced correlations from REDCOV.SUMMARY.
   See snana manual section  4.14.2 Modeling Flux Correlations.
   This feature describes correlations of the EXCESS SCATTER among 
   observations of the same event. The Poisson fluctuations are
   independent.
   Visually examing REDCOV.SUMMARY and select rho-value (reduced correlation)
   that roughly minimizes reduced chi2 between fake-data and sim.
   Following the manual syntax, manually enter the REDCOV key(s) in either 
   the sim-input or FLUXERRMODEL_FILE.

"""

# ==================================

def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP on input config file"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of input file"
    parser.add_argument("input_file", help=msg, nargs="?", default=None)

    msg = "clobber everything and start over"
    parser.add_argument("--clobber", help=msg, action="store_true")

    msg = f"start stage number (default=1, {ISTAGE_MAKEMAP}=make maps)"
    parser.add_argument("-s", "--start_stage", help=msg, nargs='?', 
                        type=int, default=1)

    msg = f"Jump to make map stage {ISTAGE_MAKEMAP} " \
          f"(previous stages already run)" 
    parser.add_argument("-m", "--makemap", help=msg, action="store_true")

    msg = "verify on fakes using map; output scales should be 1"
    parser.add_argument("--verify", help=msg, action="store_true")

    msg = "skip reduced cov check"
    parser.add_argument("--redcov_skip", help=msg, action="store_true")

    msg = "reduced cov for sim test using  FLUXERRMODEL_SIM map "
    parser.add_argument("--redcov_test", help=msg, nargs='?', 
                        type=float, default=-9.0 )

    msg = "include HBOOK file for each output table"
    parser.add_argument("--hbook", help=msg, action="store_true")

    args = parser.parse_args()

    if args.makemap : args.start_stage = ISTAGE_MAKEMAP

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return args

    # end get_args

def read_yaml(input_file):
    input_lines = []
    with open(input_file, 'r') as f :
        for line in f:
            input_lines.append(line)

    config_yaml = yaml.safe_load("\n".join(input_lines))
    return config_yaml
    # end read_yaml

def read_input(input_file):

    input_yaml = read_yaml(input_file)

    # parse VERSION_FAKES into path and version for &SNLCINP inputs
    VERSION_FAKES     = input_yaml['VERSION_FAKES']
    PRIVATE_DATA_PATH = os.path.dirname(VERSION_FAKES)
    if len(PRIVATE_DATA_PATH) < 2: PRIVATE_DATA_PATH = None

    VERSION           = os.path.basename(VERSION_FAKES)
    input_yaml['VERSION']           = VERSION
    input_yaml['PRIVATE_DATA_PATH'] = PRIVATE_DATA_PATH

    # expand FIELDs into group names and list per group
    FIELD_GROUP_NAMES = []
    FIELD_GROUP_LISTS = []
    if STRING_FIELDS in input_yaml:
        FIELDS = input_yaml[STRING_FIELDS]
        for field_group_name in FIELDS :
            field_list = FIELDS[field_group_name].split()
            FIELD_GROUP_NAMES.append(field_group_name)
            FIELD_GROUP_LISTS.append(field_list)
            #print(f" xxx field group {field_group_name} = {field_list}")

    input_yaml['NFIELD_GROUP']      = len(FIELD_GROUP_NAMES)
    input_yaml['FIELD_GROUP_NAMES'] = FIELD_GROUP_NAMES
    input_yaml['FIELD_GROUP_LISTS'] = FIELD_GROUP_LISTS

    # check for SBMAG depdendence in mag (Aug 15 2021)
    FLUXERRMAP_BINS = input_yaml['FLUXERRMAP_BINS']
    USE_SBMAG = False
    for row in FLUXERRMAP_BINS :
        if row.split()[0] == 'SBMAG': USE_SBMAG = True
    input_yaml['USE_SBMAG'] = USE_SBMAG
    

    return input_yaml
    # end read_input

def prep_outdir(ISTAGE,config):

    prefix     = stage_prefix(ISTAGE)
    input_yaml = config.input_yaml
    input_file = config.args.input_file
    args       = config.args

    key = 'OUTDIR'
    if key not in input_yaml:
        sys.exit(f"\n ERROR: missing required OUTDIR key in {input_file}\n")

    OUTDIR = input_yaml[key]

    if args.verify : 
        OUTDIR_ORIG = OUTDIR
        OUTDIR     += "_VERIFY"
        input_yaml['OUTDIR']      = OUTDIR
        input_yaml['OUTDIR_ORIG'] = OUTDIR_ORIG

    do_mkdir = False

    if os.path.exists(OUTDIR) :
        if args.clobber : 
            do_mkdir = True
            cmd = f"rm -r {OUTDIR}"
            os.system(cmd)
    else :
        do_mkdir = True

    
    if do_mkdir:  
        print(f"{prefix}: Create OUTDIR  /{OUTDIR}")
        os.mkdir(OUTDIR)
    else:
        print(f"{prefix}: Skip creating existing OUTDIR  /{OUTDIR}")

    sys.stdout.flush()

    # end prep_outdir

def get_survey_info(config):

    # run snana.
    VERSION           = config.input_yaml['VERSION']
    PRIVATE_DATA_PATH = config.input_yaml['PRIVATE_DATA_PATH']

    TEXTFILE_PREFIX   = "TEMP_GET_SURVEY" # prefix for YAML output
    yaml_file         = f"{TEXTFILE_PREFIX}.YAML"
    log_file          = f"{TEXTFILE_PREFIX}.LOG"

    if PRIVATE_DATA_PATH is None :
        arg_path = ""
    else:
        arg_path = f"PRIVATE_DATA_PATH {PRIVATE_DATA_PATH} "

    print(f" Extract SURVEY-FILTER info from {VERSION} :")

    cmd = f"{JOBNAME_SNANA} NOFILE " \
          f"VERSION_PHOTOMETRY {VERSION} " \
          f"{arg_path} " \
          f"TEXTFILE_PREFIX {TEXTFILE_PREFIX} " \
          f"MXEVT_PROCESS 0 OPT_YAML 1 > {log_file} "

    os.system(cmd)

    snana_yaml = read_yaml(yaml_file)
    survey  = snana_yaml['SURVEY']
    filters = snana_yaml['FILTERS']
    print(f"\t -> Found SURVEY-FILTERS = {survey}-{filters} ")

    cmd_rm = f"rm {yaml_file} {log_file}"
    os.system(cmd_rm)

    return survey, filters

    # read yaml file to get survey

    # end get_survey_info

def stage_prefix(ISTAGE):
    prefix = f"STAGE{ISTAGE:02d}"
    return prefix

def create_simdata(ISTAGE,config):

    # use files in nominal OUTDIR to create a sim data set that
    # uses FLUXERRMODEL_FILE(SIM) and user-input redcov.
    # Also overwrite YAML inputs as if this sim data version
    # was input instead of the fakes. Goal here is to test
    # if the minimum chi2red in output corresponds to true redcov.

    OUTDIR_ORIG   = config.input_yaml['OUTDIR']
    VERSION_ORIG  = os.path.basename(config.input_yaml['VERSION_FAKES'])
    redcov        = config.args.redcov_test
    Jrho          = int(100*redcov)

    if not os.path.exists(OUTDIR_ORIG) :
        msgerr = "\n"
        msgerr += f"ERROR: /{OUTDIR_ORIG} does not exist. \n"
        msgerr += f"Process real fakes before SIMTEST with --redcov_test .\n"
        sys.exit(msgerr)

    # overwrite selected user inputs    
    OUTDIR        = f"{OUTDIR_ORIG}_REDCOV{Jrho:03d}"
    VERSION       = f"SIMTEST_{VERSION_ORIG}"
    config.input_yaml['OUTDIR']            = OUTDIR
    config.input_yaml['VERSION_FAKES']     = VERSION
    config.input_yaml['VERSION']           = VERSION
    config.input_yaml['PRIVATE_DATA_PATH'] = None

    prefix        = stage_prefix(ISTAGE)

    print(f"{prefix}: create sim data to replace fakes")
    if config.args.start_stage > 1 :
        print(f"\t Already done --> SKIP")
        return
    
    # find simgen-input file in OUTDIR_ORIG; easier than re-creating it
    wildcard    = f"*simgen*input"
    input_files = glob.glob1(OUTDIR_ORIG,wildcard)
    n_file      = len(input_files)

    if n_file != 1 :
        msgerr = "\n"
        msgerr += f"ERROR: Found {n_file} files searching for\n"
        msgerr += f"   {wildcard}\n"
        msgerr += f" Expecting 1 and only 1 file.\n"
        msgerr += f" Files found are:\n"
        msgerr += f"   {input_files} \n"
        sys.exit(msgerr)

    simarg_redcov = prep_simarg_redcov(config,redcov)
    log_file = "out_create_simtest.log"

    cmd  = f"cd {OUTDIR_ORIG}; {JOBNAME_SIM} {input_files[0]} "
    cmd += f"GENVERSION {VERSION} "
    cmd += f"FLUXERRMODEL_FILE {FLUXERRMODEL_FILENAME_SIM} "
    cmd += f"FLUXERRMODEL_REDCOV {simarg_redcov} "
    cmd += f"FLUXERRMAP_IGNORE_DATAERR {MAP_NAME} "
    cmd += f" > {log_file} "

    #print(f"\n xxx run \n{cmd}\n")
    os.system(cmd)


    # end create_simdata

def create_fake_simlib(ISTAGE,config):

    # Run snana.exe job on fakes with option to create SIMLIB
    # where MAG column is true mag.

    OUTDIR   = config.input_yaml['OUTDIR']
    prefix   = stage_prefix(ISTAGE)

    # load name of SIMLIB file here before ISTAGE check
    SIMLIB_OUTFILE = f"{prefix}_FAKES.SIMLIB"
    config.SIMLIB_FILE = SIMLIB_OUTFILE

    print(f"{prefix}: create SIMLIB/cadence file from fakes")
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return

    KCOR_FILE         = config.input_yaml['KCOR_FILE']
    VERSION           = config.input_yaml['VERSION']
    PRIVATE_DATA_PATH = config.input_yaml['PRIVATE_DATA_PATH']
    OPT_SIMLIB        = config.input_yaml['OPT_SIMLIB']
    
    nmlarg_dict = init_nmlargs()
    nmlarg_dict[NMLKEY_DATA_PATH]   =  PRIVATE_DATA_PATH
    nmlarg_dict[NMLKEY_VERSION]     =  VERSION
    nmlarg_dict[NMLKEY_KCOR_FILE]   =  KCOR_FILE
    nmlarg_dict[NMLKEY_SNTABLE]     =  ' '
    nmlarg_dict[NMLKEY_SIMLIB_FILE] = SIMLIB_OUTFILE
    nmlarg_dict[NMLKEY_OPT_SIMLIB]  = OPT_SIMLIB

    nml_prefix = f"{prefix}_make_simlib"
    nml_file, NML_FILE = create_nml_file(config, nmlarg_dict, nml_prefix)
            
    # - - - - - - 
    run_snana_job(config, nml_file, "")
    sys.stdout.flush()

    # end create_fake_simlib

def init_nmlargs():
    nmlarg_dict = {}
    for nmlkey in NMLKEY_LIST:  nmlarg_dict[nmlkey] = None
    return nmlarg_dict
    # end init_nmlargs

def create_nml_file(config, nmlarg_dict, nml_prefix):

    input_args = config.args
    input_yaml = config.input_yaml

    OUTDIR   = input_yaml['OUTDIR']
    nml_lines_auto = []
    nml_lines_user = []

    # Create input name list file for snana.exe
    for nmlkey in NMLKEY_LIST:
        arg  = nmlarg_dict[nmlkey]
        if arg is not None:
            if isinstance(arg,str) :  arg = f"'{arg}'"
            nml_lines_auto.append(f"   {nmlkey:<20} = {arg} ")

    # - - - -
    if 'EXTRA_SNLCINP_ARGS' in config.input_yaml:
        for arg in config.input_yaml['EXTRA_SNLCINP_ARGS']:
            nml_lines_user.append(f"   {arg}")

    # - - - --  -
    # write nml_lines to nml file.

    nml_file = f"{nml_prefix}.nml"
    NML_FILE = f"{OUTDIR}/{nml_file}"
    n_lines  = len(nml_lines_auto) + len(nml_lines_user)

    print(f"\t Create {nml_file} with {n_lines} keys.")

    with open(NML_FILE,"wt") as f:
        f.write(f" &SNLCINP\n")

        f.write(f"\n ! auto-generated input\n")
        for line in nml_lines_auto: f.write(f"{line}\n")

        if input_args.verify:
            OUTDIR_ORIG = input_yaml['OUTDIR_ORIG']
            orig_file   = f"../{OUTDIR_ORIG}/{FLUXERRMODEL_FILENAME_FAKE}"
            f.write(f"   FLUXERRMODEL_FILE = '{orig_file}' \n")

        f.write(f"\n ! user-generated input\n")
        for line in nml_lines_user: f.write(f"{line}\n")
        f.write(f" &END\n") 
    
    return nml_file, NML_FILE

    # end create_nml_file

def run_snana_job(config, nml_file, args_command_line):

    args       = config.args
    input_yaml = config.input_yaml
    OUTDIR     = input_yaml['OUTDIR']

    # use prefix to construct name of log file
    prefix   = nml_file.split(".")[0]
    log_file = f"{prefix}.log"

    print(f"\t Run {JOBNAME_SNANA} on {nml_file} ")
    sys.stdout.flush()

    # run it ...
    cmd  = f"cd {OUTDIR}; {JOBNAME_SNANA} {nml_file} "
    cmd += f"{args_command_line} "
    cmd += f" > {log_file} "

    #sys.exit(f"\n xxx cmd(snana) = \n{cmd}")
    os.system(cmd)

    # end run_snana_job


def simgen_nocorr(ISTAGE,config):

    args           = config.args
    OUTDIR         = config.input_yaml['OUTDIR']
    filters        = config.filters
    SIMLIB_FILE    = config.SIMLIB_FILE
    KCOR_FILE      = config.input_yaml['KCOR_FILE']
    USE_SBMAG      = config.input_yaml['USE_SBMAG']

    # init optional keys
    HOSTLIB_FILE         = "NONE"
    HOSTLIB_MSKOPT       = None
    HOSTLIB_SBRADIUS     = None 
    HOSTLIB_SERSIC_SCALE = None
    PATH_SNDATA_SIM      = None 

    SIMKEY_REQUIRE_SBMAG = [ 'HOSTLIB_FILE', 'HOSTLIB_MSKOPT', 
                             'HOSTLIB_SBRADIUS', 
                             'HOSTLIB_SCALE_SERSIC_SIZE' ]
    HOSTLIB_DICT = {}

    if USE_SBMAG :
        for simkey in SIMKEY_REQUIRE_SBMAG :
            if simkey in config.input_yaml:
                HOSTLIB_DICT[simkey] = config.input_yaml[simkey]
            else :
                sys.exit(f"\n ERROR: SBMAG depdendence requires missing " \
                         f"{simkey} key \n\t in {args.input_file}")      

            #sys.exit(f"\n xxx HOSTLIB_DICT = {HOSTLIB_DICT}\n")

    simkey = 'PATH_SNDATA_SIM'
    if simkey in config.input_yaml :
        PATH_SNDATA_SIM = config.input_yaml[simkey]
    
    prefix         = stage_prefix(ISTAGE)
    sim_input_file = f"{prefix}_simgen_fakes.input"
    sim_log_file   = f"{prefix}_simgen_fakes.log"

    SIM_INPUT_FILE = f"{OUTDIR}/{sim_input_file}"
    SIM_LOG_FILE   = f"{OUTDIR}/{sim_log_file}"
    GENVERSION     = f"{prefix}_simgen_fakes_{USERNAME4}"
    config.SIM_GENVERSION = GENVERSION
    config.sim_input_file = sim_input_file

    ranseed = 12345 
    sim_input_lines = []

    print(f"{prefix}: run SNANA simulation using SIMLIB model")
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return

    if args.verify :
        OUTDIR_ORIG = config.input_yaml['OUTDIR_ORIG']
        orig_file   = f"../{OUTDIR_ORIG}/{FLUXERRMODEL_FILENAME_SIM}"
        sim_input_lines.append(f"FLUXERRMODEL_FILE:  {orig_file}")
        sim_input_lines.append(f" ")

    if PATH_SNDATA_SIM is not None :
        sim_input_lines.append(f"PATH_SNDATA_SIM:        {PATH_SNDATA_SIM}")

    sim_input_lines.append(f"GENVERSION:        {GENVERSION}")
    sim_input_lines.append(f"SIMLIB_FILE:       {SIMLIB_FILE}")
    sim_input_lines.append(f"SIMLIB_MSKOPT:     4      " \
                           f" # stop at end of SIMLIB file")
    sim_input_lines.append(f"NGENTOT_LC:        1000000  # any large number")
    sim_input_lines.append(f"GENSOURCE:         RANDOM")
    sim_input_lines.append(f"GENMODEL:          SIMLIB")
    sim_input_lines.append(f"GENFILTERS:        {filters}")
    sim_input_lines.append(f"KCOR_FILE:         {KCOR_FILE}")    

    if USE_SBMAG:
        sim_input_lines.append("\n# Include HOSTLIB info for SBMAG dependene");
        for simkey in HOSTLIB_DICT:
            arg = HOSTLIB_DICT[simkey]
            sim_input_lines.append(f"{simkey}:      {arg}")
        sim_input_lines.append("")
        
    sim_input_lines.append(f"RANSEED:           {ranseed} ")
    sim_input_lines.append(f"FORMAT_MASK:       32  # 2=TEXT  32=FITS ")
    sim_input_lines.append(f"SMEARFLAG_FLUX:    1   # Poisson noise from sky+source")
    sim_input_lines.append(f"SMEARFLAG_ZEROPT:  0")
    sim_input_lines.append(f"OPT_MWEBV:         1  # write MWEBV to data; not applied")
    sim_input_lines.append(f"GENSIGMA_MWEBV_RATIO: 0")
    sim_input_lines.append(f" ")
    sim_input_lines.append(f"USE_SIMLIB_PEAKMJD:   1 ")
    sim_input_lines.append(f"USE_SIMLIB_REDSHIFT:  1 ")
    sim_input_lines.append(f"GENRANGE_PEAKMJD:     40000  90000 ")
    sim_input_lines.append(f"GENRANGE_REDSHIFT:    0.012  0.98 ")
    sim_input_lines.append(f"GENRANGE_TREST:      -100 100 ")

    with open(SIM_INPUT_FILE,"wt") as f:
        for line in sim_input_lines:
            f.write(f"{line}\n")

    print(f"\t Run {JOBNAME_SIM} to generate {GENVERSION} ")
    sys.stdout.flush()

    cmd = f"cd {OUTDIR} ; {JOBNAME_SIM} {sim_input_file} > {sim_log_file}"
    os.system(cmd)

    # check for FATAL error
    f = open(f"{SIM_LOG_FILE}",  "r")
    if 'FATAL' in f.read():
        sys.exit(f"\n FATAL ERROR: check {SIM_LOG_FILE} \n")

    return
    # end simgen_nocorr

def make_outlier_table(ISTAGE,config,what):

    # run snana.exe with OUTLIER(nsig:0) to create flux table
    # for all observations.
    # Input what = FAKE or SIM

    OUTDIR = config.input_yaml['OUTDIR']
    prefix = stage_prefix(ISTAGE)
    
    nml_prefix   = f"{prefix}_fluxTable_{what}"
    table_file   = f"{nml_prefix}.{TABLE_SUFFIX_OUTLIER}"

    print(f"{prefix}: make {TABLE_NAME} table for {what}")
    sys.stdout.flush()
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return table_file

    input_yaml        = config.input_yaml
    KCOR_FILE         = input_yaml['KCOR_FILE']
    PRIVATE_DATA_PATH = input_yaml['PRIVATE_DATA_PATH']
    USE_SBMAG         = input_yaml['USE_SBMAG']

    # - - - - -
    # set arg(s) for OUTLIER table.
    # if SBMAG-dependence is set, require SBMAG<50.
    arg_outlier = 'nsig:0.0'  # default arg for OUTLIER table
    if USE_SBMAG : arg_outlier += ',sbmag:50.0'
    
    # - - - - -
    if what == STRING_FAKE :
        VERSION  = input_yaml['VERSION']
    else:
        # for SIM
        VERSION = config.SIM_GENVERSION
        
    # - - - -
    nmlarg_dict = init_nmlargs()
    nmlarg_dict[NMLKEY_DATA_PATH]   =  PRIVATE_DATA_PATH
    nmlarg_dict[NMLKEY_VERSION]     =  VERSION
    nmlarg_dict[NMLKEY_KCOR_FILE]   =  KCOR_FILE
    nmlarg_dict[NMLKEY_SNTABLE]     =  f"SNANA OUTLIER({arg_outlier})"
    nmlarg_dict[NMLKEY_TEXTFILE_PREFIX] = nml_prefix

    nml_file, NML_FILE = create_nml_file(config, nmlarg_dict, nml_prefix)
    # - - - - - - 
    run_snana_job(config, nml_file, "")

    
    # compress large TEXT tables
    print(f"\t gzip TEXT tables from {JOBNAME_SNANA} ... ")
    cmd = f"cd {OUTDIR}; gzip STAGE*.TEXT"
    os.system(cmd)

    return table_file 

    # end make_outlier_table

def parse_map_bins(config):

    # add map_bins dictionary to config
    # Each input row includes; varname  Nbin min max

    input_yaml      = config.input_yaml
    FLUXERRMAP_BINS = input_yaml['FLUXERRMAP_BINS']
    NFIELD_GROUP    = input_yaml['NFIELD_GROUP']

    # if FIELD dependent, insert FIELD as first element in map
    ivar_field = -9 
    if NFIELD_GROUP > 0 :  
        nbin=NFIELD_GROUP;  valmin=-0.5; valmax = nbin - 0.5
        row = f"{COLNAME_IFIELD}  {nbin}  {valmin}  {valmax}"
        FLUXERRMAP_BINS.insert(0,row)
        ivar_field = 0

    # if FILTER var is given, hard-wire nbin and range 
    NDIM = 0 ; ivar_filter = -9
    for row in FLUXERRMAP_BINS:
        NDIM += 1 ;   row=row.split();   varname = row[0]
        if varname == 'FILTER' or varname == 'BAND' :
            nbin = IFILTOBS_MAX; valmin = -0.5; valmax=nbin-0.5
            row = f"{COLNAME_IFILTOBS} {nbin} {valmin} {valmax}"
            FLUXERRMAP_BINS[NDIM-1] = row
            ivar_filter = NDIM - 1

    #print(f" xxx FLUXERRMAP_BINS = {FLUXERRMAP_BINS} ")

    # - - - - - - - - - - - - - - 
    map_bin_dict = {}
    varname_list=[] ;  nbin_list=[] ;  valmin_list=[]; valmax_list=[]
    bin_edge_list=[]

    #sys.exit(f"\n xxx FLUXERRMAP_BINS = \n{FLUXERRMAP_BINS}")

    NDIM   = 0 
    NBIN1D = 1
    for row in FLUXERRMAP_BINS:
        NDIM   += 1
        row     = row.split()        
        varname = row[0]
        nbin=int(row[1]);  valmin=float(row[2]); valmax=float(row[3])
        NBIN1D *= nbin
        bins    = np.linspace(valmin,valmax,nbin+1)
        print(f"    Store {nbin:2d} {varname} bins from {valmin} to {valmax}")
        varname_list.append(varname)
        nbin_list.append(nbin)
        valmin_list.append(valmin)
        valmax_list.append(valmax)
        bin_edge_list.append(bins)

    print(f"    Store total of {NBIN1D} 1D bins")
    sys.stdout.flush()

    # make list for header without FIELD or IFILTOBS 
    varname_header_list = varname_list.copy() 
    if COLNAME_IFIELD in varname_header_list:
        varname_header_list.remove(COLNAME_IFIELD)
    if COLNAME_IFILTOBS in varname_header_list :
        varname_header_list.remove(COLNAME_IFILTOBS)

    # - - - - -
    id_1d = np.arange(NBIN1D)  # 0,1,2 ... NBIN1D-1
    id_nd = []                 # array of indices for each dimension

    if NDIM == 1 :
        id0 = np.unravel_index(id_1d,(nbin_list[0]))
        id_nd.append(id0)
        indexing_array = np.arange(NBIN1D).reshape((nbin_list[0]))
    elif NDIM == 2 :
        id0,id1 = np.unravel_index(id_1d,(nbin_list[0],nbin_list[1]))
        id_nd.append(id0)
        id_nd.append(id1)
        indexing_array = np.arange(NBIN1D).reshape((nbin_list[0],nbin_list[1]))
    elif NDIM == 3 :
        id0,id1,id2 = \
            np.unravel_index(id_1d,(nbin_list[0],nbin_list[1],nbin_list[2]))
        id_nd.append(id0)
        id_nd.append(id1)
        id_nd.append(id2)
        indexing_array = \
            np.arange(NBIN1D).reshape((nbin_list[0],nbin_list[1],nbin_list[2]))
    elif NDIM == 4 :
        id0,id1,id2,id3 = \
            np.unravel_index(id_1d,(nbin_list[0],nbin_list[1],nbin_list[2],nbin_list[3]))
        id_nd.append(id0)
        id_nd.append(id1)
        id_nd.append(id2)
        id_nd.append(id3)
        indexing_array = \
            np.arange(NBIN1D).reshape((nbin_list[0],nbin_list[1],nbin_list[2],nbin_list[3]))

    else :
        sys.exit("\n ERROR: cannot process NDIM={NDIM}\n")

    #print(f" xxx nbin_list = {nbin_list} ")
    #print(f" xxx id_nd = \n{id_nd}")
    #print(f" xxx indexing_array=\n{indexing_array} ")
    #print(f" xxx  idem(1,2) = {indexing_array[1,2]} ")

    map_bin_dict = {
        'NDIM'          : NDIM,       # number of map dimensions
        'NBIN1D'        : NBIN1D,      # total number of multiD bins
        'NVAR'          : len(varname_list),
        'varname_list'  : varname_list ,
        'varname_header_list'  : varname_header_list ,
        'nbin_list'     : nbin_list ,
        'valmin_list'   : valmin_list ,
        'valmax_list'   : valmax_list ,
        'bin_edge_list' : bin_edge_list ,   # histogram bin edges
        'id_1d'         : id_1d,
        'id_nd'         : id_nd,
        'indexing_array': indexing_array,
        'ivar_field'    : ivar_field,
        'ivar_filter'   : ivar_filter,   # flag to make filter-dependent maps
    }

    config.map_bin_dict = map_bin_dict
    sys.stdout.flush()

    # end parse_map_bins

def make_fluxerr_model_map(ISTAGE,config):

    #
    # Driver to read OUTLIER tables, compute fluxerr maps,
    # and write maps for SNANA simulation. 
    # Output is one map from fakes (for real data) and
    # one map for SNANA sim.
    #  

    OUTDIR       = config.input_yaml['OUTDIR']
    map_bin_dict = config.map_bin_dict
    prefix       = stage_prefix(ISTAGE)
    nfilters     = len(config.filters)
        
    fluxerrmodel_file_fake = f"{OUTDIR}/{FLUXERRMODEL_FILENAME_FAKE}"
    fluxerrmodel_file_sim  = f"{OUTDIR}/{FLUXERRMODEL_FILENAME_SIM}"

    print(f"{prefix}: create FLUXERRMODEL maps. ")
    sys.stdout.flush()
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return

    flux_table_fake = f"{OUTDIR}/{config.flux_table_fake}"
    flux_table_sim  = f"{OUTDIR}/{config.flux_table_sim}"
    if not os.path.exists(flux_table_fake):  flux_table_fake += '.gz'
    if not os.path.exists(flux_table_sim):   flux_table_sim  += '.gz'

    # read each table
    df_fake = store_flux_table(flux_table_fake, map_bin_dict)
    df_sim  = store_flux_table(flux_table_sim,  map_bin_dict)

    # load list of unique ifiltobs & band into 
    # map_bin_dict.ifiltobs_set, band_set
    get_filter_list(df_fake, map_bin_dict)

    # add index columns, force bounds, apply cuts ...
    df_fake, df_sim = modify_tables(df_fake, df_sim, config)

    # - - - - -
    #print(f"\n xxx df_fake = \n{df_fake[['LOGSNR','PSF', 'i_LOGSNR']]}")
    #print(f"\n xxx df_fake = \n{df_fake}")
    # - - - - - 

    # open output map files
    f_fake = open(fluxerrmodel_file_fake,"wt")
    f_sim  = open(fluxerrmodel_file_sim, "wt")

    write_map_global_header(f_fake, STRING_FAKE, config)
    write_map_global_header(f_sim,  STRING_SIM,  config)

    # get a few things for 1D loop over bins
    NBIN1D        = map_bin_dict['NBIN1D']
    id_nd         = map_bin_dict['id_nd'] 
    ivar_filter   = map_bin_dict['ivar_filter']
    ivar_field    = map_bin_dict['ivar_field']
    ifiltobs_list = map_bin_dict['ifiltobs_list']
    ifiltobs_last = -9
    ifield_last   = -9

    # if there is no filter dependence, write one header 
    # with BAND specifying all bands.
    if ivar_filter < 0 :
        write_map_header(f_fake, -9, config)
        write_map_header(f_sim,  -9, config)

    # start loop over 1D bins (which loops over all dimensions of map)

    print(f" Begin loop over {NBIN1D} 1D map bins ... ")
    sys.stdout.flush()
    errscale_dict = \
        { 'bin1d':[], 'n_fake':[], 'n_sim':[], 'cor_fake':[], 'cor_sim':[] }
    nrow = 0
    if nfilters > 0 : 
        NFIELD_GROUP    = config.input_yaml['NFIELD_GROUP']
        nrow_per_filter = int(NBIN1D/IFILTOBS_MAX/NFIELD_GROUP)
    
    for BIN1D in range(0,NBIN1D):

        ifield = -9
        if ivar_field >= 0: ifield = id_nd[ivar_field][BIN1D]

        # check for start of filter-dependent map
        use_filter = True
        if ivar_filter >= 0:
            ifiltobs   = id_nd[ivar_filter][BIN1D]
            use_filter = ifiltobs in ifiltobs_list
            if use_filter and ifiltobs != ifiltobs_last :
                write_map_header(f_fake, ifield, ifiltobs, config)
                write_map_header(f_sim,  ifield, ifiltobs, config)
                nrow = 0
                errscale_dict = \
                    { 'bin1d':[], 'n_fake':[], 'n_sim':[], 'cor_fake':[], 'cor_sim':[] }
            ifiltobs_last = ifiltobs

        if not use_filter : continue  # skip the pad zeros in ifiltobs_list

        # select sample in this BIN1D (this multi-D cell)
        pull_fake = \
            df_fake['PULL'].to_numpy()[np.where(df_fake[COLNAME_BIN1D]==BIN1D)]
        pull_sim = \
            df_sim['PULL'].to_numpy()[np.where(df_sim[COLNAME_BIN1D]==BIN1D)]
        ratio_fake = \
            df_fake['ERR_RATIO'].to_numpy()[np.where(df_fake[COLNAME_BIN1D]==BIN1D)]

        # compute errScale correction for fake and sim
        n_fake, n_sim, cor_fake, cor_sim = \
            compute_errscale_cor ( pull_fake, pull_sim, ratio_fake )

        nrow += 1
        errscale_dict['bin1d'].append(BIN1D)
        errscale_dict['n_fake'].append(n_fake)
        errscale_dict['n_sim'].append(n_sim)
        errscale_dict['cor_fake'].append(cor_fake)
        errscale_dict['cor_sim'].append(cor_sim)

        # set write flag on last filter or last bin
        if ivar_filter >=0 :
            write_flag = (nrow == nrow_per_filter)
        else:
            write_flag = (nrow == NBIN1D)

        #print(f" xxx nrow={nrow}  ifiltobs={ifiltobs}  write_flag={write_flag}")
        
        if write_flag :
            errscale_extrap(errscale_dict)
            for i in range(0,nrow):
                bin1d    = errscale_dict['bin1d'][i]
                n_fake   = errscale_dict['n_fake'][i]
                n_sim    = errscale_dict['n_sim'][i]
                cor_fake = errscale_dict['cor_fake'][i]
                cor_sim  = errscale_dict['cor_sim'][i]                
                write_map_row(f_fake, config, bin1d, cor_fake, n_fake, -9)
                write_map_row(f_sim,  config, bin1d, cor_sim,  n_fake, n_sim )
        
        # xxx mark delete
        # update map files.
        #write_map_row(f_fake, config, BIN1D, cor_fake, n_fake, -9)
        #write_map_row(f_sim,  config, BIN1D, cor_sim,  n_fake, n_sim )
        # xxxx
        
    # - - - 
    print("\n")
    print(f" Done creating {fluxerrmodel_file_fake} ")
    print(f" Done creating {fluxerrmodel_file_sim} ")
    sys.stdout.flush()

    return
    # end make_fluxerr_model_map


def  errscale_extrap(errscale_dict):
    # if cor = 0.0, this means there are not enough fake in this bin
    # to compute a correction. In this case, set cor = cor(nearest bin).
    # This functions returns modified errscale_dict with cor=0 replaced
    # with cor=cor(nearest).

    cor_fake_orig = errscale_dict['cor_fake']
    cor_sim_orig  = errscale_dict['cor_sim']

    # init extrap correction 
    cor_fake_extrap = []
    cor_sim_extrap  = []

    npcor_fake_orig = np.array(cor_fake_orig)
    npcor_sim_orig  = np.array(cor_sim_orig)
    non_zeros       = np.nonzero(npcor_fake_orig)[0]

    i = 0
    for cor_fake,cor_sim in zip(cor_fake_orig,cor_sim_orig) :
        if cor_fake == 0.0 :
            distances   = np.abs(non_zeros - i)
            closest_idx = np.min(np.where(distances == np.min(distances)))
            cor_fake_near  = npcor_fake_orig[non_zeros[closest_idx]]
            cor_sim_near   = npcor_sim_orig[non_zeros[closest_idx]]
            cor_fake       = cor_fake_near
            cor_sim        = cor_sim_near
            
        cor_fake_extrap.append(cor_fake)
        cor_sim_extrap.append(cor_sim)
        i += 1

    errscale_dict['cor_fake'] = cor_fake_extrap
    errscale_dict['cor_sim']  = cor_sim_extrap

    # end errscale_extrap
    
def  modify_tables(df_fake, df_sim, config):

    # Modify tables by
    #  + applying cuts to reject some rows
    #  + force map variable values to lie within map bounds
    #  + add bin-index columns


    map_bin_dict   = config.map_bin_dict
    args           = config.args
    input_yaml     = config.input_yaml

    NDIM           = map_bin_dict['NDIM']
    NBIN1D         = map_bin_dict['NBIN1D']
    varname_list   = map_bin_dict['varname_list']
    nbin_list      = map_bin_dict['nbin_list']
    valmin_list    = map_bin_dict['valmin_list']
    valmax_list    = map_bin_dict['valmax_list']
    bin_edge_list  = map_bin_dict['bin_edge_list']
    id_1d          = map_bin_dict['id_1d'] 
    id_nd          = map_bin_dict['id_nd'] 
    ivar_filter    = map_bin_dict['ivar_filter']

    # apply optional cuts on NSIG and ERR_RATIO

    nrow_orig_fake = len(df_fake)
    nrow_orig_sim  = len(df_sim)

    key = 'CUTWIN_NSIG'
    if key in input_yaml :
        cutwin_nsig = input_yaml[key].split()
        nsig_max    = float(cutwin_nsig[1])
        df_fake     = df_fake.loc[ df_fake['NSIG'] < nsig_max ]
        df_sim      = df_sim.loc[  df_sim['NSIG']  < nsig_max ]

        # reset indices
        df_fake     = df_fake.reset_index()
        df_sim      = df_sim.reset_index()

    #sys.exit(f"\n xxx df_fake = \n{df_fake}\n")

    nrow_fake   = len(df_fake)
    nrow_sim    = len(df_sim)
    config.nrow_fake  = nrow_fake
    config.nrow_sim   = nrow_sim

    print(f" Nrow(fake) = {nrow_orig_fake} -> {nrow_fake} after cuts.")
    print(f" Nrow(sim)  = {nrow_orig_sim} -> {nrow_sim} after cuts.")
    sys.stdout.flush()

    # - - - - 
    # force variables in map to lie within map ranges e.g., 
    # if LOGSNR has 5 bins from 0.3 to 2.3, then LOGSNR<0.3 -> 0.30001
    force_bounds(df_fake, config)
    force_bounds(df_sim,  config)

    # assign integer IDFIELD = 0, 1, 2, ... NFIELD_GROUP-1 to each row
    NFIELD_GROUP = input_yaml['NFIELD_GROUP']
    if NFIELD_GROUP > 0 :
        FIELD_LISTS = input_yaml['FIELD_GROUP_LISTS']
        print(f"   Add {COLNAME_IFIELD} column to tables ...")
        df_fake[COLNAME_IFIELD] = \
            df_fake.apply(lambda row: apply_field(row,FIELD_LISTS), axis=1)

        df_sim[COLNAME_IFIELD] = \
            df_sim.apply(lambda row: apply_field(row,FIELD_LISTS), axis=1)
        
    #sys.exit(f"\n xxx BYE BYE df_fake=\n{df_fake}\n")

    # - - - - - 
    # digitize the variables used in map; i.e., compute multi-D
    # index for each map-variable and each row in table.
    ivar=0
    ibins_fake = [ ];   ibins_sim = [ ]  # init arrays of multi-D indices
    for varname, bins in zip(varname_list, bin_edge_list):
        dfcol_fake = df_fake[varname]
        dfcol_sim  = df_sim[varname]
        ibins_fake.append(np.digitize(dfcol_fake, bins))
        ibins_sim.append(np.digitize( dfcol_sim,  bins))
        ivarname = f"i_{varname}"  
        df_fake[ivarname] = ibins_fake[ivar] - 1
        df_sim[ivarname]  = ibins_sim[ivar]  - 1
        ivar += 1

    # add 1d index colum to each flux table to enable easy selection
    # of multi-D cells from 1D index.
    print(f"   Add {COLNAME_BIN1D} column to tables ...")
    df_fake[COLNAME_BIN1D] = \
        df_fake.apply(lambda row: apply_id_1d(row,map_bin_dict), axis=1)

    df_sim[COLNAME_BIN1D] = \
        df_sim.apply(lambda row: apply_id_1d(row,map_bin_dict), axis=1)

    return df_fake, df_sim

    # end modify_tables

def force_bounds(df,config):

    map_bin_dict = config.map_bin_dict
    varname_list = map_bin_dict['varname_list']
    valmin_list  = map_bin_dict['valmin_list']
    valmax_list  = map_bin_dict['valmax_list']

    for varname,valmin,valmax in zip(varname_list,valmin_list,valmax_list):
        if varname == COLNAME_IFILTOBS : continue
        if varname == COLNAME_IFIELD   : continue
        print(f"\t Force {valmin} < {varname} < {valmax}")
        df[varname] = df[varname].where(df[varname] > valmin, valmin+0.0001)
        df[varname] = df[varname].where(df[varname] < valmax, valmax-0.0001)
        
    sys.stdout.flush()
    # end force_bounds


def compute_errscale_cor(pull_fake, pull_sim, ratio_fake):
    
    # for this 1D bin, compute ERRSCALE correction for FAKE(data) and SIM.
    #  + pull_fake a list of PULL = (F-Ftrue)/ERR_CALC
    #  + pull_sim  is the same for sime
    #  + ratio_fake is list of ERR_DATA/ERR_CALC [fakes]
    #
    #  From  Sec 6.4 of https://arxiv.org/pdf/1811.02379.pdf 
    #
    #  Eq 15 for FAKE (intended for data)
    #     scale = RMS[(F-Ftrue)/ERRCALC]_fake / <ERR_F/ERR_CALC>_fake
    #
    #  Eq 14 for SIM (intended for sim)
    #     scale = RMS[(F-Ftrue)/ERRCALC]_fake / RMS[(F-Ftrue)/ERRCALC]_sim

    n_fake = len(pull_fake)
    n_sim  = len(pull_sim)
    if n_fake > 5 and n_sim > 5 :

        # shift pulls so that median/avg is zero
        avg_pull_fake  = np.median(pull_fake)
        avg_pull_sim   = np.median(pull_sim)
        pull_fake      = pull_fake - avg_pull_fake
        pull_sim       = pull_sim  - avg_pull_sim

        # for RMS, compute 1.48*median|pull| to reduce sensitivity to outliers
        rms_pull_fake  = 1.48 * np.median(np.absolute(pull_fake))
        rms_pull_sim   = 1.48 * np.median(np.absolute(pull_sim))

        avg_ratio      = np.median(ratio_fake)  # ERR_DATA/ERR_CALC

        # finally, the map corrections
        cor_fake       = rms_pull_fake / avg_ratio     # correct fake & data
        cor_sim        = rms_pull_fake / rms_pull_sim  # correct sims

        if FORCE_ERRCORR1 :
            if cor_fake < 1.0 : cor_fake = 1.001
            if cor_sim  < 1.0 : cor_sim  = 1.001

        #print(f"\t xxx cor_sim = {rms_pull_fake:.3f} / {rms_pull_sim:.3f}" \
        #      f" = {cor_fake:.3f}  " \
        #      f" (avgPull={avg_pull_fake:.3f},{avg_pull_sim:0.3f}) " )

    else:
        # set cor=0 here and later it will be extrapolated to nearest bin
        cor_fake = 0.0 ; cor_sim = 0.0


    return n_fake, n_sim, cor_fake, cor_sim

    # end compute_errscale_cor

def write_map_global_header(f, what, config):
    
    # write global comments and one-time keys.
    # what = FAKE or SIM
    # If map depends on fields, define each field group.

    input_yaml        = config.input_yaml
    NFIELD_GROUP      = input_yaml['NFIELD_GROUP']
    FIELD_GROUP_NAMES = input_yaml['FIELD_GROUP_NAMES']
    FIELD_GROUP_LISTS = input_yaml['FIELD_GROUP_LISTS']

    map_bin_dict      = config.map_bin_dict
    NDIM              = map_bin_dict['NDIM']
    varname_list      = map_bin_dict['varname_list']

    nrow_fake = config.nrow_fake
    nrow_sim  = config.nrow_sim

    # for filter and field, replace internal index names with
    # more reasonable FIELD and BAND
    varname_string    = ' '.join(varname_list)
    varname_string    = varname_string.replace(COLNAME_IFIELD,"FIELD")
    varname_string    = varname_string.replace(COLNAME_IFILTOBS,"BAND")

    if what == STRING_SIM:
        usage_code  = "snana.exe and snlc_sim.exe" ; 
        usage_key   = "FLUXERRMODEL_FILE"
        item        = what
        string_nrow = f"Nobs(FAKE,SIM) = {nrow_fake}, {nrow_sim}"
    else:
        usage_code  = "snlc_fit.exe" 
        usage_key   = "FLUXERRMODEL_FILE in &SNLCINP"
        item        = f"FAKES and DATA"
        string_nrow = f"Nobs(FAKE) = {nrow_fake} "

    # - - - 
    f.write(f"DOCUMENTATION:")
    f.write(f"""
  PURPOSE: correct flux uncertainty for {item}
  REF:
  - AUTHOR: Kessler et al, 2019 (DES3YR sims, see Sec 6.4)
    ADS:    https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.1171K
  INTENT:  Test
  USAGE_CODE:  {usage_code}
  USAGE_KEY:   {usage_key}
  NOTES:
  - map dependence is {varname_string}
  - {string_nrow}
  - map-create command =  {sys.argv[0]} {sys.argv[1]} 
  - created by user={USERNAME} on HOST={HOSTNAME}  
    """)
    f.write(f"\nDOCUMENTATION_END:\n")
    f.write(f"\n\n")

    for field_name, field_list in zip(FIELD_GROUP_NAMES,FIELD_GROUP_LISTS):
        snana_field_list = '+'.join(field_list)
        f.write(f"DEFINE_FIELDGROUP: {field_name}  {snana_field_list}\n")

    f.flush()
    # end write_map_define_fields

def write_map_header(f, ifield, ifiltobs, config):

    map_bin_dict  = config.map_bin_dict
    input_yaml    = config.input_yaml

    FIELD_GROUP_NAMES   = input_yaml['FIELD_GROUP_NAMES']
    band_map            = map_bin_dict['band_map']
    band_list           = map_bin_dict['band_list']
    varname_header_list = map_bin_dict['varname_header_list']
    varlist             = ' '.join(varname_header_list)

    if ifield >=0 : 
        field_arg = FIELD_GROUP_NAMES[ifield]

    if ifiltobs < 0 :
        band_arg = ''.join(band_list)
    else:
        band_arg = band_map[ifiltobs]

    f.write("\n")
    f.write(f"MAPNAME: {MAP_NAME} \n")
    if ifield >= 0 : f.write(f"FIELD: {field_arg} \n")
    f.write(f"BAND: {band_arg} \n")
    f.write(f"VARNAMES: {varlist}   ERRSCALE\n")

    f.flush()
    # end write_map_header

def write_map_row(f, config, BIN1D, cor, n_fake, n_sim ):

    # write map row corresponding to 1D index BIN1D
    # cor is the ERRSCALE correction 

    map_bin_dict    = config.map_bin_dict
    bin_edge_list   = map_bin_dict['bin_edge_list']
    NVAR            = map_bin_dict['NVAR']
    id_nd           = map_bin_dict['id_nd']
    varname_list    = map_bin_dict['varname_list']
    ivar_field      = map_bin_dict['ivar_field']
    ivar_filter     = map_bin_dict['ivar_filter']
    nbin_list       = map_bin_dict['nbin_list']

    #print(f" bin_list = {bin_list}")
    # convert to multi-D indices and extract grid values
    row_line = 'ROW: '
    last_row = True

    for ivar in range(0,NVAR):
        if ivar == ivar_field  : continue
        if ivar == ivar_filter : continue
        itmp       = id_nd[ivar][BIN1D]
        lo_edge    = bin_edge_list[ivar][itmp]
        hi_edge    = bin_edge_list[ivar][itmp+1]
        bin_center = 0.5*(lo_edge + hi_edge)
        row_line  += f"{bin_center:8.4f}  "
        if itmp < nbin_list[ivar]-1 : last_row = False

    row_line += f"{cor:8.3f}"

    # prepare comment with stats
    if n_sim > 0:
        comment = f"N_[FAKE,SIM] = {n_fake} , {n_sim}"
    else :
        comment = f"N_FAKE = {n_fake} "

    f.write(f"{row_line}    # {comment} \n")

    if last_row:   
        f.write("ENDMAP:\n\n")
        f.flush()

    # end write_map_row

def store_flux_table(flux_table, map_bin_dict):

    df = pd.read_csv(flux_table, comment="#", delim_whitespace=True)
    nrow = len(df)
    print(f"    Read/store {flux_table} with {nrow} rows.")

    STR_F        = 'FLUXCAL_DATA' ; 
    STR_FTRUE    = 'FLUXCAL_TRUE' ; 
    STR_ERR      = 'FLUXCAL_ERR_DATA'
    STR_ERR_CALC = 'FLUXCAL_ERR_CALC'

    # compute modified PULL with ERR -> ERR_CALC
    pull = (df[STR_F]-df[STR_FTRUE])/df[STR_ERR_CALC] 
    df['PULL'] = pull.values

    err_ratio  = df[STR_ERR] / df[STR_ERR_CALC]
    df['ERR_RATIO'] = err_ratio.values

    return df

    # end store_flux_table

def get_filter_list(df, map_bin_dict):

    # get unique list of IFILTOBS and BAND 

    ifiltobs_list = sorted(list(set(df[COLNAME_IFILTOBS])))
    band_list     = []
    band_map      = [ -9 ]*100

    for ifiltobs in ifiltobs_list :
        indx = df[df[COLNAME_IFILTOBS] == ifiltobs].index[0]
        band = df[COLNAME_BAND][indx]
        band_list.append(band)
        band_map[ifiltobs] = band

    band_string = ''.join(band_list)
    print(f"    Found IFILTOBS set = {ifiltobs_list}")
    print(f"       -->    BAND set = {band_list}  ({band_string})")

    map_bin_dict['ifiltobs_list']  = ifiltobs_list
    map_bin_dict['band_list']      = band_list
    map_bin_dict['band_map']       = band_map
    map_bin_dict['band_string']    = band_string
    return
    # end get_filter_list


def apply_field(row,FIELD_LISTS):

    # return IFIELD index for this row.
    # Input FIELD_LISTS is a list of lists; e..g,
    # [ ['X3','C3'] , ['S1', 'S2', 'X1', 'X2'] ]

    FIELD = row['FIELD']
    ifield = 0
    for field_list in FIELD_LISTS:
        if FIELD in field_list: return ifield
        ifield += 1

    # if we get here, abort
    sys.exit(f"\n ERROR: FIELD={FIELD} is not in \n\t {FIELD_LISTS}. " \
             f"\n\t See FIELDS arg in input file.")

    return ifield
    # end apply_field

def apply_id_1d(row, map_bin_dict):

    # return 1D index for this row.
    
    id_1d          = map_bin_dict['id_1d']
    NDIM           = map_bin_dict['NDIM']
    varname_list   = map_bin_dict['varname_list']
    nbin_list      = map_bin_dict['nbin_list']
    indexing_array = map_bin_dict['indexing_array']
    ib_list       = []

    for varname, nbin in zip(varname_list, nbin_list) : 
        ivarname = f"i_{varname}" 
        ib       = row[ivarname]
        ib_list.append(ib)

    # - - - 
    # beware: this NDIM logic is goofy & fragile :(

    SMART_WAY = False

    if SMART_WAY:
        # couldn't get this to work; need somebody even smarter ??
        id_1d = indexing_array[np.ix_([[ib_list[i]] for i in ib_list])]
    else:
        # dumb way with NDIM if-block
        if NDIM == 1:
            id_1d  = indexing_array[ib_list[0]]
        elif NDIM == 2 :
            id_1d  = indexing_array[ib_list[0],ib_list[1]]
        elif NDIM == 3 :
            id_1d  = indexing_array[ib_list[0],ib_list[1],ib_list[2]]
        elif NDIM == 4 :
            id_1d = indexing_array[ib_list[0],ib_list[1],ib_list[2],ib_list[3]]
    return id_1d
    # end apply_id_1d


def create_nml_redcov(ISTAGE,config):

    # create nml file to analyze for reduced cov with SNANA table.

    prefix = stage_prefix(ISTAGE)

    nml_prefix = f"{prefix}_redcov"
    print(f"{prefix}: create NML file for analyzing redcov. ")
    sys.stdout.flush()
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return f"{nml_prefix}.nml"

    KCOR_FILE         = config.input_yaml['KCOR_FILE']
    VERSION           = config.input_yaml['VERSION']
    PRIVATE_DATA_PATH = config.input_yaml['PRIVATE_DATA_PATH']
    map_file          = FLUXERRMODEL_FILENAME_FAKE

    nmlarg_dict = init_nmlargs()
    nmlarg_dict[NMLKEY_DATA_PATH]   =  PRIVATE_DATA_PATH
    nmlarg_dict[NMLKEY_VERSION]     =  VERSION
    nmlarg_dict[NMLKEY_KCOR_FILE]   =  KCOR_FILE
    nmlarg_dict[NMLKEY_SNTABLE]     =  'SNANA'
    nmlarg_dict[NMLKEY_TEXTFILE_PREFIX]   = nml_prefix
    nmlarg_dict[NMLKEY_LPROB_TRUEFLUX]    = True

    if config.args.redcov_test < 0.0 :
        nmlarg_dict[NMLKEY_MODEL_FILE]  = map_file

    nml_file, NML_FILE = create_nml_file(config, nmlarg_dict, nml_prefix)
 
    args_cmd_line = ""

    # for redcov test using sims, use special key to force sim to
    # use flux-error map just like real data
    if config.args.redcov_test > -0.001 :
        args_cmd_line += f"{NMLKEY_SIM_MODEL_FILE} {map_file} "

    if config.args.hbook :
        args_cmd_line += f"{NMLKEY_HFILE_OUT} {nml_prefix}.HBOOK "

    # run it for data; analyze sims with overrides in next stage
    run_snana_job(config, nml_file, args_cmd_line )

    return nml_file

    # end create_nml_redcov

def redcov_simgen_plus_snana(ISTAGE,config,rho):

    # + simulate with reduced cov "rho" using sim-input file from STAGE02.
    # + Run SNANA job using nml file created in previous stage

    survey       = config.survey
    prefix       = stage_prefix(ISTAGE)
    Jrho         = int(100*rho)  # used for file names
    GENVERSION   = f"{prefix}_{survey}_fakes_REDCOV{Jrho:03d}"

    print(f"{prefix}: generate and process sim with rho = {rho}")
    sys.stdout.flush()
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return GENVERSION

    sim_input_file = config.sim_input_file
    nml_file       = config.nml_file_redcov    
    OUTDIR         = config.input_yaml['OUTDIR']

    # construct special sim arg for reduced cov;
    simarg_redcov = prep_simarg_redcov(config,rho)

    log_file_sim   = f"{GENVERSION}_SIM.LOG"
    log_file_snana = f"{GENVERSION}_SNANA.LOG"
  
    # - - - - - - - 
    # create and run sim command
    cmd  = f"cd {OUTDIR}; {JOBNAME_SIM} {sim_input_file} "
    cmd += f"GENVERSION {GENVERSION} "
    cmd += f"FLUXERRMODEL_FILE {FLUXERRMODEL_FILENAME_SIM} "
    if rho > 0.0 :  cmd += f"FLUXERRMODEL_REDCOV {simarg_redcov}"
    cmd += f" > {log_file_sim}"

    print(f"\t Generage {GENVERSION}")
    sys.stdout.flush()

    #sys.exit(f"\n xxx cmd=\n{cmd}")
    os.system(cmd)

    # - - - - - -
    # create and run snana command:
    cmd  = f"cd {OUTDIR}; {JOBNAME_SNANA} {nml_file} "
    cmd += f"VERSION_PHOTOMETRY {GENVERSION} "
    cmd += f"TEXTFILE_PREFIX {GENVERSION} "

    # xxxx mark delete xxxx
    #if config.args.redcov_test > -0.01 :
    #    cmd += f"{NMLKEY_SIM_MODEL_FILE} NONE "

    if config.args.hbook :
        cmd += f"{NMLKEY_HFILE_OUT} {GENVERSION}.HBOOK "

    cmd += f" > {log_file_snana}"
    print(f"\t Run {JOBNAME_SNANA} on {GENVERSION} to produce SNANA table")
    os.system(cmd)

    return GENVERSION

    # end redcov_simgen_plus_snana

def prep_simarg_redcov(config,rho):

    # construct sim argument for REDCOV key
    # e.g, filters= gri and rho = 0.4
    # simarg_redcov = g:0.4,r:0.4,i:0.4

    rho_arg = rho
    if rho == 1.0 : rho_arg = 0.99

    filters        = config.filters
    filter_list    = list(filters)

    simarg_redcov_list = []
    for band in filter_list: 
        simarg_redcov_list.append(f"{band}:{rho_arg}")

    simarg_redcov = ','.join(simarg_redcov_list)

    return simarg_redcov
    # end prep_simarg_redcov

def redcov_analyze(ISTAGE,config, rho, prefix_sim, f_summary ):
    
    # + Analyze SNANA tables (data & sim) and compute FoM
    # + Update summary table using passed pointer f_summary

    prefix       = stage_prefix(ISTAGE)
    print(f"{prefix}: Compute FoM per band for rho = {rho}")
    sys.stdout.flush()
    if ISTAGE < config.args.start_stage :
        print(f"\t Already done --> SKIP")
        return
    
    OUTDIR           = config.input_yaml['OUTDIR']
    nml_file         = config.nml_file_redcov
    filters          = config.filters
    filter_list      = list(filters)
    filter_list.append(KEY_ALLBANDS)

    prefix_data      = nml_file.split('.')[0]
    table_file_data  = f"{prefix_data}.{TABLE_SUFFIX_SNANA}"
    table_file_sim   = f"{prefix_sim}.{TABLE_SUFFIX_SNANA}"
    print(f"\t Compare data to {prefix_sim}")

    # fom_list is FoM per band
    chi2red_list = compute_chi2red(config, table_file_data, table_file_sim )

    # update summay file
    chi2red_sum = 0.0
    for band,chi2red in zip(filter_list,chi2red_list):

        if band == KEY_ALLBANDS: 
            comment = '# include ALL bands for chi2red'
        else:
            chi2red_sum += chi2red
            comment  = ''

        f_summary.write(f"  - {band:<3} {rho:4.2f}  {chi2red:7.2f}" \
                        f"   {comment}\n")

    # - - - - 
    f_summary.write(f"  - SUM {rho:4.2f}  {chi2red_sum:7.2f} " \
                    f"  # sum of chi2red over bands\n\n")

    # end redcov_analyze

def compute_chi2red(config, table_file_data, table_file_sim):

    # compute chi2red(data/sim) for each band;
    # return chi2red_list vs. band.

    OUTDIR         = config.input_yaml['OUTDIR']
    filters        = config.filters
    filter_list    = list(filters)
    filter_list.append(KEY_ALLBANDS)
    chi2red_list   = []

    TABLE_FILE_DATA = f"{OUTDIR}/{table_file_data}"
    TABLE_FILE_SIM  = f"{OUTDIR}/{table_file_sim}"
    df_data = pd.read_csv(TABLE_FILE_DATA, comment="#", delim_whitespace=True)
    df_sim  = pd.read_csv(TABLE_FILE_SIM,  comment="#", delim_whitespace=True)

    nrow_data = len(df_data)
    nrow_sim  = len(df_sim)

    prob_min  = 0.10
    prob_max  = 1.00
    nbin_prob = 9   # number of histogram bins
    prob_bins = np.linspace(prob_min, prob_max, nbin_prob+1)

    print(f"\t nrow(data,sim) = {nrow_data} , {nrow_sim} ")
    dump_flag = False
    if dump_flag : print(" xxx ----------------------------------- ")

    
    for band in filter_list:

        if band == KEY_ALLBANDS :
            varname_prob     = f"{VARNAME_PROB_PREFIX}"
        else:
            varname_prob     = f"{VARNAME_PROB_PREFIX}_{band}"

        # make sure all PROB values are < 1.00
        df_data[varname_prob] = \
            df_data[varname_prob].where(df_data[varname_prob] < 1.0, 0.9999)
        df_sim[varname_prob] = \
            df_sim[varname_prob].where(df_sim[varname_prob] < 1.0, 0.9999)

        prob_data        = df_data[varname_prob]
        prob_sim         = df_sim[varname_prob]
        prob_data_binned = np.digitize(prob_data, prob_bins)
        prob_sim_binned  = np.digitize(prob_sim,  prob_bins)

        # bin data and sim; exclude underflow bin 0
        bin_data       = np.bincount(prob_data_binned)[1:]
        bin_sim_raw    = np.bincount(prob_sim_binned)[1:]
        sim_scale      = np.sum(bin_data)/np.sum(bin_sim_raw)

        if dump_flag:
            print(f" xxx {band}: bin_data = {bin_data}")
            print(f" xxx {band}: bin_sim  = {bin_sim_raw} x {sim_scale:.3f}")
            sys.stdout.flush()

        # compute chi2 (there must be a better way; this is awful ?!?!?)

        bin_sim_norm   = sim_scale * bin_sim_raw
        bin_sim_norm2  = sim_scale * sim_scale * bin_sim_raw
        bin_dif        = np.subtract(bin_data,bin_sim_norm)
        bin_variance   = np.add(bin_data,bin_sim_norm2)
        bin_chi2       = np.divide(bin_dif*bin_dif,bin_variance)
        chi2           = np.sum(bin_chi2)
        chi2red        = chi2 / float(nbin_prob-1)

        if dump_flag :
            print(f" xxx chi2red = {chi2:.2f}/{nbin_prob} = {chi2red:.2f}" )
            sys.stdout.flush()

        chi2red_list.append(chi2red)

    return chi2red_list

    # end compute_chi2red

def create_redcov_summary_file(ISTAGE,config):
    OUTDIR       = config.input_yaml['OUTDIR']
    summary_file = f"{OUTDIR}/{REDCOV_SUMMARY_FILE}"

    prefix       = stage_prefix(ISTAGE)
    print(f"{prefix}: Open {REDCOV_SUMMARY_FILE}")

    f_summary    = open(summary_file,"wt")
    f_summary.write("# Data/sim chi2red vs. Reduced Flux Correlation \n\n")
    f_summary.write("REDCOV:   # band  rho  chi2red(data/sim)\n")
    return f_summary
    # end create_redcov_summary_file

# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    config      = Namespace()
    config.args = get_args()

    # option for long HELP menus
    if config.args.HELP :
        see_me = (f" !!! ************************************************ !!!")
        print(f"\n{see_me}\n{see_me}\n{see_me}")
        print(f"{HELP_CONFIG}")
        sys.exit(' Scroll up to see full HELP menu.\n Done: exiting Main.')

    config.input_yaml = read_input(config.args.input_file)

    # run snana job to extract name of SURVEY from fake data
    config.survey, config.filters = get_survey_info(config)

    ISTAGE = 0
    print()

    # check option to replace fake data with sim data that includes
    # fluxerr mape and user-input redcov.
    if config.args.redcov_test >= 0.0 :
        create_simdata(ISTAGE,config)

    # - - - - - - - - - - - 
    ISTAGE += 1
    prep_outdir(ISTAGE,config)

    # - - - - - -
    # create simlib from fakes
    ISTAGE += 1
    create_fake_simlib(ISTAGE,config)

    # simulate fakes with snana
    ISTAGE += 1
    simgen_nocorr(ISTAGE,config)

    # run snana on fakes and sim; create OUTLIER table with nsig>=0
    # to catch all flux observations
    ISTAGE += 1
    config.flux_table_fake = make_outlier_table(ISTAGE,config,STRING_FAKE)
    config.flux_table_sim  = make_outlier_table(ISTAGE,config,STRING_SIM)

    # create fluxerrmodel map files
    ISTAGE += 1
    parse_map_bins(config)
    make_fluxerr_model_map(ISTAGE,config)

    if config.args.verify:   sys.exit(0)

    if config.args.redcov_skip:  sys.exit(0)

    # - - - - - - - -
    print("")
    # stages below examine correlations in excess scatter
    # by comparing data to sims with different RHO

    ISTAGE += 1
    config.nml_file_redcov = create_nml_redcov(ISTAGE,config)

    ISTAGE += 1
    prefix_sim_list = []
    for rho in REDCOV_LIST:
        prefix_sim = redcov_simgen_plus_snana(ISTAGE,config,rho)
        prefix_sim_list.append(prefix_sim)

    ISTAGE += 1
    f_summary = create_redcov_summary_file(ISTAGE,config)
    for rho,prefix_sim in zip(REDCOV_LIST,prefix_sim_list):
        redcov_analyze(ISTAGE,config, rho, prefix_sim, f_summary)
        
# === END ===

