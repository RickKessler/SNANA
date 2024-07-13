#!/usr/bin/env python
#
# Created May 2024 by R.Kessler
#
# Wrapper to run launch submit_batch_jobs.sh multiple times.
# Read list file containing list of submit_batch input files to 
# launch with submit_batch_jobs.sh. Submit each set, waiting for
# ALL.DONE before launching next set.
# Multiple input files on a line are launched together as a set;
# each successive line is launched when the ALL.DONEs exist with
# SUCCESS.
#
# July 9 2024: remove scancel --me ... later, need to scancel specific pids
#
# ==================================================

import os, sys, datetime, shutil, subprocess, time, glob, yaml, argparse
import getpass, logging

SNANA_DIR        = os.getenv('SNANA_DIR')
USERNAME         = getpass.getuser()
HOSTNAME         = os.uname()[1]
CWD              = os.getcwd()
SUBMIT_LOG_DIR   = CWD + '/submit_logs'

KEY_SUBMIT_LIST = 'SUBMIT_LIST'
KEY_SUBMIT_DIR  = 'SUBMIT_DIR'

SUBMIT_JOB_NAME     = "submit_batch_jobs.sh"


MERGE_LOG_FILE     = "MERGE.LOG"
SUBMIT_INFO_FILE   = "SUBMIT.INFO"
ALL_DONE_FILE      = "ALL.DONE"

STRING_SUCCESS = "SUCCESS"

WAIT_TIME_CHECK_DONE    =  20.0   # time (sec) between checking for done files
WAIT_TIME_STDOUT_UPDATE = 600.0   # time (sec) betweeing printing status to stdout 

# =======================================

def parse_args():

    parser = argparse.ArgumentParser()

    msg = "print help menu for config file"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "Submit batch job name"
    parser.add_argument("--jobname", help=msg, default=SUBMIT_JOB_NAME, 
                        type=str)

    msg = f"private snana devel directory (replaces {SNANA_DIR}"
    parser.add_argument("--snana_dir", help=msg, default=None, 
                        type=str)

    msg = f"Number of rows to skip in job-list file"
    parser.add_argument("--nrowskip", help=msg, default=0, 
                        type=int)

    msg = f"file with list of submit_batch input files."
    parser.add_argument("-l", "--list_file", help=msg, default=None, 
                        type=str)

    args = parser.parse_args()

    if args.snana_dir is None:
        args.snana_dir = SNANA_DIR

    if args.snana_dir is not None :
        args.jobname = f"{args.snana_dir}/util/{SUBMIT_JOB_NAME}"

    if args.HELP : 
        print_help_menu()
                    
    logging.info(f" Inputs: ")
    logging.info(f"   SNANA_DIR       = {args.snana_dir} ")
    logging.info(f"   SUBMIT_JOB_NAME = {args.jobname}")
    logging.info(f" Submit-HOSTNAME: {HOSTNAME} ")

    if args.list_file is None:
        sys.exit(f"\n ERROR: missing required --list_file arg.")

    return args

# ===========================
def print_help_menu():

    help_menu = f"""

# example input config file for testing {SUBMIT_JOB_NAME}

# Optional location of submit_batch input files;
# if location is not given, cwd is assumed.
SUBMIT_DIR: $SNANA_TESTS/inputs_submit_batch 

# Each row is a list of submit_batch inputs to launch in parallel.
# Next row waits for all jobs in previous row to finish.

SUBMIT_LIST:
# first set of jobs
  - SIMGEN_MASTER_DES.INPUT  SIMGEN_MASTER_LOWZ.INPUT
  - LCFIT_DES_DATA+SIM.NML   LCFIT_LOWZ_DATA+SIM.NML
  - BBC_DATA+SIM.INPUT       BBC_DATA+SIM_0+0.INPUT  BBC_DATA+SIM_w0wa.INPUT
  - CREATE_COV_DATA+SIM.INPUT
  - WFIT_DATA+SIM.INPUT


# 2nd set of jobs starts with LCFIT applied to above sims
  - LCFIT_DES_1REPEAT.NML  LCFIT_LOWZ_1REPEAT.NML
  - BBC_1REPEAT.INPUT      BBC_1REPEAT_NSPLITRAN.INPUT

# 3rd set of jobs
  - SIMGEN_MASTER_DES_1CHANGE.INPUT  SIMGEN_MASTER_LOWZ_1CHANGE.INPUT
  - LCFIT_DES_1CHANGE.NML  LCFIT_LOWZ_1CHANGE.NML
  - BBC_1CHANGE.INPUT
  - CREATE_COV_1CHANGE.INPUT
  - WFIT_1CHANGE.INPUT

# 4th set of jobs
  - LCFIT_DES_FITOPT_SPARSE.NML  LCFIT_LOWZ_FITOPT_SPARSE.NML
  - BBC_FITOPT_MAP.INPUT
"""
    sys.exit(f"{help_menu}")

    return

# ============================
def setup_logging():

    #logging.basicConfig(level=logging.DEBUG,
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)6s] %(message)s")
    #format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)

    return

def read_list_file(input_submit_file):

    if not os.path.exists(input_submit_file):
        sys.exit(f"\n ERROR: cannot find {input_submit_file}")

    config     = extract_yaml(input_submit_file)
    submit_dir = config.setdefault(KEY_SUBMIT_DIR,CWD)
    config[KEY_SUBMIT_DIR] = os.path.expandvars(submit_dir)

    return config

def scancel_all_jobs(submit_info_file):

    if os.path.exists(submit_info_file):
        # scancel pid list read from submit_info_file
        logging.info(f"scancel all slurm jobs in {submit_info_file}")
        submit_info_contents = extract_yaml(submit_info_file)
        SBATCH_LIST          = submit_info_contents['SBATCH_LIST']
        for tmp in SBATCH_LIST :
            cpu     = tmp[0]    # cpu id = 0,1,2, etc ...
            pid     = tmp[1]
            jobname = tmp[2]
            logging.info(f"\t scancel {pid} for {jobname}")
            
            cmd_cancel = f"scancel {pid}"
            os.system(cmd_cancel)
            
        logging.info(f"EXIT after scancel on all pids")
    else:
        logging.info(f"WARNING: cannot find {submit_info_file}")
        logging.info(f"EXIT without scancel on pids")

    sys.exit()
    return

def check_file_exists(file_name):
    # abort if file does not exist
    exist = os.path.isfile(file_name)
    if not exist:
        msgerr = f"ERROR: Cannot find required file: {file_name}"
        logging.info(f"{msgerr}")
        scancel_all_jobs("dum.log")
        sys.exit(f"\n ABORT {sys.argv[0]}")

    # end check_file_exists

def check_done_file(done_file):
    # Return true if done_file exists and has SUCCESS in it.
    # Return false if done_file does not exist.
    # If done_file exists with FAIL, cancel jobs and abort.

    dirname = os.path.dirname(done_file)
    submit_info_file = f"{dirname}/{SUBMIT_INFO_FILE}" # in case of error
    
    if os.path.isfile(done_file):
        with open(done_file,"rt") as f:
            line = f.read() ;  status = line.rstrip("\n")
            if status != STRING_SUCCESS :
                scancel_all_jobs(submit_info_file)
                msgerr = f"\n ERROR: found {status} in \n\t {done_file}\n"
                sys.exit(msgerr)
        return True
    else:
        return False

    # end check_done_file
            
def extract_yaml(input_file):

    # after expandvars on input_file, read and return yaml contents
    
    input_file = os.path.expandvars(input_file)

    exist  = os.path.isfile(input_file)

    lines = []
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#END_YAML"):
                break
            lines.append(line)
    config = yaml.safe_load("\n".join(lines))
    return config
 
def get_outdir_list(config):

    # for each submit-input file, read & store output dir.
    # Under CONFIG yaml block, read either OUTDIR or LOGDIR key.

    outdir_list = []
    submit_dir  = config[KEY_SUBMIT_DIR]

    for submit in config[KEY_SUBMIT_LIST]:
        set_list = ''
        for infile in submit.split() :
            INFILE = (f"{submit_dir}/{infile}")
            input_yaml = extract_yaml(INFILE)

            if 'CONFIG' in input_yaml:
                CONFIG     = input_yaml['CONFIG']
            else:
                msgerr = f"\n ERROR: cannot find CONFIG in {INFILE}"
                sys.exit(msgerr)

            outdir = None 
            if 'OUTDIR' in CONFIG :
                outdir = CONFIG['OUTDIR']
            if 'LOGDIR' in CONFIG :
                outdir = CONFIG['LOGDIR']

            set_list += f"{outdir} "

        outdir_list.append(set_list)

    #sys.exit(f" xxx outdir_list = \n{outdir_list}")
    return outdir_list
    # end parse_outdir


def run_submit(infile_list, outdir_list, SUBMIT_INFO ):

    INPUTS = SUBMIT_INFO['INPUTS']  # command line inputs
    config = SUBMIT_INFO['config']  # contents of config input file

    # infile_list is space-separated list of input files
    # to launch with submit_batch_jobs
    SUBMIT_JOB_NAME = os.path.expandvars(INPUTS.jobname)
    snana_dir       = INPUTS.snana_dir
    submit_dir      = config[KEY_SUBMIT_DIR]

    if not os.path.exists(SUBMIT_LOG_DIR):
        os.mkdir(SUBMIT_LOG_DIR)

    done_file_list = []
    info_file_list = []
    for infile, outdir in zip(infile_list, outdir_list) :
        merge_file = f"{submit_dir}/{outdir}/{MERGE_LOG_FILE}"
        info_file  = f"{submit_dir}/{outdir}/{SUBMIT_INFO_FILE}"
        done_file  = f"{submit_dir}/{outdir}/{ALL_DONE_FILE}"      
        info_file_list.append(info_file)
        done_file_list.append(done_file)
        logging.info(f" submit {infile}  -> {outdir}")

        arg_list   = [ infile ]
        arg_string = infile
        if snana_dir != SNANA_DIR:  
            arg_list   = [ infile, f"--snana_dir", f"{snana_dir}" ]
            arg_string = f"{infile} --snana_dir {snana_dir}"

        prefix          = infile.split('.')[0]
        submit_log_file = f"{SUBMIT_LOG_DIR}/submit_{prefix}.log"

        cmd_list   = [ SUBMIT_JOB_NAME ] + arg_list
        cmd_string = f"{SUBMIT_JOB_NAME} {arg_string} >& {submit_log_file} "

        ret = subprocess.run( [cmd_string], cwd=submit_dir,
                              shell=True, capture_output=True, text=True )

        time.sleep(1)
        check_file_exists(merge_file)
        check_file_exists(info_file)
        time.sleep(3)

    # - - - - - - - 
    # wait for done files
    NDONE_EXPECT = len(done_file_list)
    NDONE_FIND   = 0
    NDONE_LAST   = -9
    t_proc_sec   = 0.0
    t_proc_upd   = 0.0
    
    while NDONE_FIND < NDONE_EXPECT :
        NDONE_FIND = 0
        time.sleep(WAIT_TIME_CHECK_DONE)
        t_proc_sec += WAIT_TIME_CHECK_DONE
        t_proc_upd += WAIT_TIME_CHECK_DONE
        t_proc_min  = t_proc_sec/60.0

        for done_file in done_file_list :
            found = check_done_file(done_file)
            if found : NDONE_FIND += 1


        DO_STDOUT_UPD = (NDONE_FIND > NDONE_LAST) or \
                        (t_proc_upd >= WAIT_TIME_STDOUT_UPDATE)
        if DO_STDOUT_UPD:
            msg = f"\t found {NDONE_FIND} of {NDONE_EXPECT} {ALL_DONE_FILE} files " \
                  f"({t_proc_min:.1f} minutes)"
            logging.info(msg)
            t_proc_upd = 0.0

        NDONE_LAST = NDONE_FIND

    # - - - - - -
    logging.info(f"\t Finished set with {STRING_SUCCESS}" )
    return

    # end run_submit

# ===================================
# ============== MAIN ===============
# ===================================

if __name__ == "__main__":

    setup_logging()

    job_string = ' '.join(sys.argv)
    logging.info(f"Begin {job_string}")

    # parse input arguments
    INPUTS = parse_args()

    SUBMIT_INFO = { 'INPUTS' : INPUTS }

    # read list of input files 
    logging.info(f" Read submit-batch input files from  {INPUTS.list_file}")
    config = read_list_file(INPUTS.list_file)

    SUBMIT_INFO.update( {'config' : config} )
    infile_submit_list = config[KEY_SUBMIT_LIST]
    submit_dir         = config[KEY_SUBMIT_DIR]
    
    nset = 0
    for infile_set in infile_submit_list:
        nset += 1
        if nset <= INPUTS.nrowskip : continue
        logging.info(f"   inFile set:  {infile_set}")
       
    outdir_submit_list = get_outdir_list(config)

    SUBMIT_INFO.update( { 'outdir_submit_list' : outdir_submit_list } )
    for outdir_set in outdir_submit_list:
        logging.info(f"   outdir set: {outdir_set}")


    logging.info(f"\n# - - - - - - - - - - - - - - - - - -"); 
    logging.info(f" Output subDirs under {submit_dir}\n")

    nset =  0
    for infile_set, outdir_set in zip(infile_submit_list,outdir_submit_list):
        nset += 1
        if nset <= INPUTS.nrowskip : continue
        infile_list  = infile_set.split()
        outdir_list  = outdir_set.split()

        #submit, wait for ALL.DONE
        run_submit(infile_list, outdir_list, SUBMIT_INFO )  

    # - - - -
    logging.info(f"")
    logging.info(f"Done.")
    logging.info(f"All jobs reported {STRING_SUCCESS} in {ALL_DONE_FILE}")
    sys.exit(0)

    # === END ===

    
