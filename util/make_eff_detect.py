#!/usr/bin/env python
#
# Created Dec 2022 by R.Kessler and B.Sanchez
#
# Create map of EFF(DIA detection) vs. SNR_calc, where SNR_calc is 
# calculated SNR(ZP,PSF,SKY), and is the same as SNR compuation in 
# the simulation. The resulting map is intended as sim-input via key 
# SEARCHEFF_PIPELINE_EFF_FILE
#
# The Eff-vs-SNR is a function of band and FIELD;
# bands are determined by the data and FIELD groups are
# specified in the input config file for this script.
#
# 
# =================

import os, sys, argparse, glob, yaml, math
import numpy as np
from   argparse import Namespace
import pandas as pd

# --------------------------------
# define a few handy globals

JOBNAME_SNANA = "snana.exe"

USERNAME      = os.environ['USER']
USERNAME4     = os.environ['USER'][0:4]
HOSTNAME      = os.uname()[1].split('.')[0]


# We use the snana 'OUTLIER' table that was designed to write
# obervations for Nsigma outliers. Here we don't care about outliers
# and use this feature with Nsig >=0 to get everything.
TABLE_NAME    = "OUTLIER"

STRING_FIELDS = "FIELDS"

COLNAME_IFILTOBS = "IFILTOBS"
COLNAME_BAND     = "BAND"
COLNAME_IFIELD   = "IFIELD"

NMLKEY_DATA_PATH   = 'PRIVATE_DATA_PATH'
NMLKEY_VERSION     = 'VERSION_PHOTOMETRY'
NMLKEY_KCOR_FILE   = 'KCOR_FILE'
NMLKEY_SNTABLE     = 'SNTABLE_LIST'
NMLKEY_TEXTFILE_PREFIX   = 'TEXTFILE_PREFIX'

NMLKEY_LIST = [ NMLKEY_DATA_PATH, NMLKEY_VERSION, NMLKEY_KCOR_FILE,
                NMLKEY_SNTABLE, NMLKEY_TEXTFILE_PREFIX ]

TABLE_SUFFIX_SNANA   = "SNANA.TEXT"
TABLE_SUFFIX_OUTLIER = "OUTLIER.TEXT"


def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP on input config file"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of input config file"
    parser.add_argument("input_file", 
                        help=msg, nargs="?", default=None)

    msg = "clobber everything and start over"
    parser.add_argument("--clobber", 
                        help=msg, action="store_true")

    msg = "skip making flux table and make map from existing table"
    parser.add_argument("-s", "--skip_fluxtable",
                        help=msg, action="store_true")

    args = parser.parse_args()

    # ?? if args.makemap : args.start_stage = ISTAGE_MAKEMAP

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

    # read and parse input config file

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

    do_mkdir = False

    if os.path.exists(OUTDIR) :
        if args.clobber : 
            do_mkdir = True
            cmd      = f"rm -r {OUTDIR}"
            os.system(cmd)
    else :
        do_mkdir = True

    
    if do_mkdir:  
        print(f"{prefix}: Create OUTDIR:    /{OUTDIR}")
        os.mkdir(OUTDIR)
    else:
        print(f"{prefix}: Skip creating existing OUTDIR  /{OUTDIR}")

    sys.stdout.flush()
    return

    # end prep_outdir

def get_survey_info(config):

    # later replace hard-wired values with call to quick_commands
    # using YAML output
    survey = 'DES'
    filters = 'griz'

    return survey, filters
    # end get_survey_info

def stage_prefix(ISTAGE):
    prefix = f"STAGE{ISTAGE:02d}"
    return prefix

def init_nmlargs():
    nmlarg_dict = {}
    for nmlkey in NMLKEY_LIST:  nmlarg_dict[nmlkey] = None
    return nmlarg_dict
    # end init_nmlargs                                                                 

def create_nml_file(config, nmlarg_dict, nml_prefix):

    # To do: convert this info quick_command utility ??

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

        f.write(f"\n ! user-generated input\n")
        for line in nml_lines_user: f.write(f"{line}\n")

        f.write(f" &END\n") 
    
    return nml_file, NML_FILE

    # end create_nml_file

def run_snana_job(config, nml_file, args_command_line):

    # to do: convert to quick_command util ?

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

    os.system(cmd)
    # end run_snana_job         


def make_flux_table(ISTAGE,config):

    # run snana.exe with OUTLIER(nsig:0) to create flux table
    # for all observations.

    OUTDIR = config.input_yaml['OUTDIR']
    prefix = stage_prefix(ISTAGE)
    
    nml_prefix   = f"{prefix}_fluxTable"
    table_file   = f"{nml_prefix}.{TABLE_SUFFIX_OUTLIER}"

    print(f"{prefix}: make {TABLE_NAME} table")
    sys.stdout.flush()

    if config.args.skip_fluxtable:
        print(f"\t SKIP this stage and use existing {table_file}")
        return table_file

    input_yaml        = config.input_yaml
    KCOR_FILE         = input_yaml['KCOR_FILE']
    PRIVATE_DATA_PATH = input_yaml['PRIVATE_DATA_PATH']
    VERSION           = input_yaml['VERSION']

    # - - - - -
    # set arg(s) for OUTLIER table.
    arg_outlier = 'nsig:0.0'  # default arg for OUTLIER table
            
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

    # end make_flux_table

def read_table(ISTAGE, config):

    input_yaml = config.input_yaml

    OUTDIR     = input_yaml['OUTDIR']
    prefix     = stage_prefix(ISTAGE)
    table_file = config.table_file    

    print(f"{prefix}: read {table_file} ")
    sys.stdout.flush()

    table_file_path = f"{OUTDIR}/{table_file}.gz"
    df = pd.read_csv(table_file_path, comment="#", delim_whitespace=True)

    # assign integer IDFIELD = 0, 1, 2, ... NFIELD_GROUP-1 to each row         
    NFIELD_GROUP = input_yaml['NFIELD_GROUP']
    if NFIELD_GROUP > 0 :
        FIELD_LISTS = input_yaml['FIELD_GROUP_LISTS']
        print(f"   Add {COLNAME_IFIELD} column to table ...")
        df[COLNAME_IFIELD] = \
            df.apply(lambda row: get_IFIELD(row['FIELD'],config), axis=1)

    print(f"\n xxx df = \n{df} \n")

    return df
    # end read_table

def get_IFIELD(FIELD,config):

    # return IFIELD index for this FIELD string input..
    # Input FIELD_LISTS is a list of lists; e..g,
    # [ ['X3','C3'] , ['S1', 'S2', 'X1', 'X2'] ]
    #
    # Note that overlap fields (e.g., S1+S2) will return -9

    input_yaml     = config.input_yaml
    FIELD_LISTS    = input_yaml['FIELD_GROUP_LISTS']

    ifield = 0
    for field_list in FIELD_LISTS:
        if FIELD in field_list: return ifield
        ifield += 1

    return -9
    # end get_IFIELD

def effsnr_binned(ISTAGE, config):

    prefix     = stage_prefix(ISTAGE)
    print(f"{prefix}: compute EFF(detect) in  SNR_calc bins");
    sys.stdout.flush()

    effsnr_binned_dict = {}

    return effsnr_binned_dict
    # end effsnr_binned

def effsnr_fit(ISTAGE, config):

    prefix   = stage_prefix(ISTAGE)
    print(f"{prefix}: fit sigmoid to smooth EFF(detect) vs SNR_calc");
    sys.stdout.flush()

    effsnr_fit_dict = {}

    return effsnr_fit_dict

    # end effsnr_fit

def write_map(ISTAGE, config):
    
    map_file = f"SEARCHEFF_PIPELINE_EFF_{config.survey}.DAT"

    prefix   = stage_prefix(ISTAGE)
    print(f"{prefix}: write map to {map_file}")
    sys.stdout.flush()

    # end write_map

# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    config      = Namespace()
    config.args = get_args()

    # option for long HELP menus
    # ...

    config.input_yaml = read_input(config.args.input_file)

    # run snana job to extract name of SURVEY from fake data
    config.survey, config.filters = get_survey_info(config)

    ISTAGE = 0
    print()

    # - - - - - - - - - - - 
    ISTAGE += 1
    prep_outdir(ISTAGE,config)

    # run snana on fakes (or sim); create OUTLIER table with nsig>=0
    # to include all flux observations
    ISTAGE += 1
    config.table_file = make_flux_table(ISTAGE,config)

    # =================================
    #  xxxxxxx BELOW IS FOR BRUNO
    # =================================

    # read table and apply cuts
    ISTAGE += 1    
    config.df_table = read_table(ISTAGE, config)

    # xxxxxx need strategy to loop over FIELD group and band xxxxxx
    
    # construct Efficiency in SNR_calc bins
    ISTAGE += 1    
    config.effsnr_binned_dict = effsnr_binned(ISTAGE, config)

    # smooth eff-vs-SNR with fit to sigmoid
    ISTAGE += 1    
    config.effsnr_smooth_dict = effsnr_fit(ISTAGE, config)

    # finally, write the map for the simulation
    ISTAGE += 1
    write_map(ISTAGE, config)

    # END:


