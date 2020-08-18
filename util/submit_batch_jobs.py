#!/usr/bin/env python
#
# Created July 2020 by R.Kessler & S. Hinton
#
# TO-DO LIST for
#  BASE/util: 
#   - more elegant HELP menu per program?
#   - implement ssh (on FNAL cluster)
#   - translate legacy input file (w/Justin)
#  SIM:
#   - for sim, leave symbolic links for redundant sim job
#   - problem reading SIMGEN-input file when SIMGEN_DUMP breaks
#      to another line that is not YAML compatible
#
#  FIT:
#   - track down why NEVT(HBOOK) sometimes fails
#  BBC
#
# - - - - - - - - - -

#import os
import sys
import yaml
import argparse
import logging
import submit_batch.submit_util as util

from   submit_batch.submit_params   import *
from   submit_batch.submit_prog_sim import Simulation
from   submit_batch.submit_prog_fit import LightCurveFit
from   submit_batch.submit_prog_bbc import BBC
from   argparse import Namespace

# =====================================
def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "name of input file"
    parser.add_argument("input_file", help=msg, nargs="?", default=None)

    # misc user args
    msg = "increase output verbosity"
    parser.add_argument("-v", "--verbose", help=msg, action="store_true")

    msg = "KILL jobs"
    parser.add_argument("-k", "--kill", help=msg, action="store_true")

    msg = "Create & init outdir, but do NOT submit jobs"
    parser.add_argument("-n", "--nosubmit", help=msg, action="store_true")

    msg = "process x10 fewer events for sim,fit,bbc (applies only to sim data)"
    parser.add_argument("--fast", help=msg, action="store_true")

    msg = "debug mode: submit jobs, but skip merge process"
    parser.add_argument("--nomerge", help=msg, action="store_true")

    msg = (f"debug mode: reset merge process (nevt->0, rm merged files ...)"\
           f" to allow interactive -m")
    parser.add_argument("--merge_reset", help=msg, action="store_true")

    # args passed internally from command files
    msg = "INTERNAL: launch merge process"
    parser.add_argument("-m", "--merge", help=msg, action="store_true")

    msg = "INTERNAL: launch merge when all done files exist"
    parser.add_argument("-M", "--MERGE", help=msg, action="store_true")

    msg = "INTERNAL: time stamp (Nsec since midnight) to verify" \
           " merge process examines correct output_dir"
    parser.add_argument('-t', nargs='+', help=msg, type=int )

    msg = "INTERNAL: cpu number "
    parser.add_argument('--cpunum', nargs='+', help=msg, type=int )

    msg = "HELP with input file config(s)"
    parser.add_argument("-H", "--HELP", help=msg, default=None, type=str, \
            choices=["SIM", "FIT", "BBC"])
    args = parser.parse_args()

    return parser.parse_args()

def which_program_class(config):

    # check YAML/input_file keys to determine which program class
    # (sim,fit,bbc) will run

    program_class = None 
    input_file    = config['args'].input_file
    merge_flag    = config_yaml['args'].merge_flag

    if "GENVERSION_LIST" in config :
        program_class = Simulation
    elif "VERSION" in config['CONFIG'] :
        program_class = LightCurveFit
    else:
        # BBC does not have a unique required batch key,
        # so instead check for unique input key for code
        with open(input_file, 'r') as f :    
            if 'u1=' in f.read():
                program_class = BBC

    # keep quiet for merge process
    if not merge_flag :
        logging.info(f" Program class : {program_class.__name__} \n" )

    return program_class

def set_merge_flag(config):
    merge_flag = config['args'].merge  or \
                 config['args'].MERGE  or \
                 config['args'].merge_reset
    return merge_flag

def check_legacy_input_file(input_file):
    with open(input_file,"r") as f:
        flat_word_list=[word for line in f for word in line.split()]
        #f_read = f.read()

    if 'CONFIG:' in flat_word_list :
        # check for obsolete keys
        for item in flat_word_list :
            key = item.rstrip(':')
            if key in OBSOLETE_CONFIG_KEYS :
                msgerr = []
                msgerr.append(f" Obsolete key '{key}' no longer valid.")
                util.log_assert(False,msgerr)

        return  # file ok, do nothing.
    else:
        # define fake config_yaml, and pass to program driver
        config_yaml = \
            { 'args' : Namespace(input_file=input_file, legacy_input=True) }

    if  'GENVERSION:' in flat_word_list :
        program = Simulation(config_yaml)  

    elif 'VERSION:' in flat_word_list :
        program = LightCurveFit(config_yaml) 

    elif 'u1=' in str(flat_word_list) :  # check for u1= substring
        program = BBC(config_yaml) 
    else:
        print(f" xxx word_list = {flat_word_list}")
        msgerr = ['Unrecognized legacy input file:', input_file ]
        util.log_assert(False,msgerr)

    #sys.exit("\n xxx DEBUG DIE check legacy xxx ")
    # end check_legacy_input_file

def print_submit_messages(config_yaml):

    # print final info to screen for user
    CONFIG = config_yaml['CONFIG']

    print(f" Done launching jobs. Sit back and relax.")

    if 'OUTDIR' in CONFIG :
        OUTDIR = CONFIG['OUTDIR']
        MERGE_LOG = (f"{OUTDIR}/{MERGE_LOG_FILE}")
        print(f" Check status in {MERGE_LOG} ")

    if config_yaml['args'].nomerge :
        print(f" REMEMBER: you disabled the merge process.")

    if config_yaml['args'].fast :
        print(f" REMEMBER: fast option will process 1/{FASTFAC} of request.")

    # end print_submit_messages

def print_nosubmit_messages(config_yaml):
    # print final info to screen for user
    CONFIG = config_yaml['CONFIG']
    if 'OUTDIR' in CONFIG :
        print(f"\n Check job preparation in {CONFIG['OUTDIR']}/ ")

    print(f" Jobs NOT sumbitted. Bye Bye.")

    # end print_nosubmit_messages

# =============================================
if __name__ == "__main__":
    args  = get_args()
    store = util.setup_logging(args)

    if args.HELP :
        print(f"{HELP_CONFIG[args.HELP]}")
        sys.exit(' Done')

    # check for legacy input; if so, translate and quit
    check_legacy_input_file(args.input_file)

    # Here we know it's got a CONFIG block, so read the YAML input
    config_yaml = util.extract_yaml(args.input_file)
    config_yaml['args'] = args

    #sys.exit(f" xxx config_yaml = {config_yaml} ")

    # set merge flag before running program_class
    config_yaml['args'].merge_flag   = set_merge_flag(config_yaml)
    config_yaml['args'].legacy_input = False

    logging.debug(config_yaml)

    # - - - - - -
    # determine which program class (sim, fit, bbc)
    program_class = which_program_class(config_yaml)

    # - - - - - -
    # run the class
    program = program_class(config_yaml)  # calls __init only

    # check merge options
    if config_yaml['args'].merge_flag :
        program.merge_driver()
        print('  Done with merge process in Main -> exit.')
        exit(0)

    # check option to kill jobs 
    if config_yaml['args'].kill :
        kill_jobs(config_prep)

    # create output dir
    program.create_output_dir()

    # prepare files, lists, program args
    program.submit_prepare_driver() 

    # write .BATCH and .CMD scripts
    program.write_script_driver()
    
    # Create MERGE.LOG file with all jobs in WAIT state.          
    # This file gets updated later by merge process.
    program.create_merge_file()

    # create SUBMIT.INFO file for merge process ... 
    # unlike MERGE.LOG, this file never changes.
    program.create_info_file()

    if args.nosubmit :
        print_nosubmit_messages(config_yaml)
        exit(0)

    # submit jobs via batch or ssh
    program.launch_jobs()

    # Print any warnings or errors that we captured at the end to make
    # sure they arent missed
    store.print_warnings()
    store.print_errors()
    
    print_submit_messages(config_yaml)

    exit(0)

# === END ===
