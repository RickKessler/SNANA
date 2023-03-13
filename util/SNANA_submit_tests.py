#!/usr/bin/env python
#
# Read list of batch jobs to submit, submit each set,
# wait for ALL.DONE to launch next set.
# Multiple input files on a line are launched together as a set;
# each successive line is launched when ALL.DONEs exist with
# SUCCESS.
#
# Aug 09 2021: add --snana_dir and --nrowskip args
# Jan 20 2022: run create_coariance.py if CREATE_COV is in input file name
# Mar 08 2023: add -l arg to switch list file
#
# ==================================================

import os, sys, datetime, shutil, subprocess, time, glob, yaml, argparse
import getpass

USERNAME         = getpass.getuser()
CWD              = os.getcwd()
SNANA_DIR        = os.environ['SNANA_DIR']
SNANA_TESTS_DIR  = os.environ['SNANA_TESTS'] + '/inputs_submit_batch'

SUBMIT_LIST_FILE_DEFAULT = f"{SNANA_TESTS_DIR}/SNANA_submit_tests.LIST"

KEY_SUBMIT_LIST = 'SUBMIT_LIST'

SUBMIT_JOB_NAME     = "submit_batch_jobs.sh"
CREATE_COV_JOB_NAME = "create_covariance.py"

MERGE_LOG_FILE     = "MERGE.LOG"
SUBMIT_INFO_FILE   = "SUBMIT.INFO"
ALL_DONE_FILE      = "ALL.DONE"

STRING_SUCCESS = "SUCCESS"

PREFIX_CREATE_COV = "CREATE_COV" # for input file

# =======================================

def parse_args():

    parser = argparse.ArgumentParser()

    msg = "Submit batch job name"
    parser.add_argument("--jobname", help=msg, default=SUBMIT_JOB_NAME, 
                        type=str)

    msg = f"private snana devel directory (replaces {SNANA_DIR}"
    parser.add_argument("--snana_dir", help=msg, default=None, 
                        type=str)

    msg = f"Number of rows to skip in job-list file"
    parser.add_argument("--nrowskip", help=msg, default=0, 
                        type=int)

    msg = f"alternate job-list file"
    parser.add_argument("-l", "--list_file", help=msg, 
                        default=SUBMIT_LIST_FILE_DEFAULT, 
                        type=str)

    args = parser.parse_args()

    if args.snana_dir is None:
        args.snana_dir = SNANA_DIR

    if args.snana_dir is not None :
        args.jobname = f"{args.snana_dir}/util/{SUBMIT_JOB_NAME}"

    print(f" Inputs: ")
    print(f"   SNANA_DIR       = {args.snana_dir} ")
    print(f"   SUBMIT_JOB_NAME = {args.jobname}")

    return args

def scancel_all_jobs():
    cmd_cancel = (f"scancel --user={USERNAME}")
    print(f"{cmd_cancel}"); sys.stdout.flush()
    os.system(cmd_cancel)

def check_file_exists(file_name):
    # abort if file does not exist
    exist = os.path.isfile(file_name)
    if not exist:
        msgerr = (f"\n ERROR: Cannot find required file: \n\t {file_name}\n")
        print(f"{msgerr}")
        scancel_all_jobs()
        sys.exit(f"\n ABORT {sys.argv[0]}")

    # end check_file_exists

def check_done_file(done_file):
    # Return true if done_file exists and has SUCCESS in it.
    # Return false if done_file does not exist.
    # If done_file exists with FAIL, cancel jobs and abort.

    if os.path.isfile(done_file):
        with open(done_file,"rt") as f:
            line = f.read() ;  status = line.rstrip("\n")
            if status != STRING_SUCCESS :
                scancel_all_jobs()
                msgerr = (f"\n ERROR: found {status} in \n\t {done_file}\n")
                sys.exit(msgerr)
        return True
    else:
        return False

    # end check_done_file
            
def extract_yaml(input_file):
    msgerr = [(f"Cannot find the input yaml file:\n   {input_file}")]
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

    # for each submit-input file, read & store output.
    # Under CONFIG yaml block, read either OUTDIR or LOGDIR key.

    outdir_list = []
    
    for submit in config[KEY_SUBMIT_LIST]:
        #print(f" xxx submit  = {submit} ")
        set_list = ''
        for infile in submit.split() :
            INFILE = (f"{SNANA_TESTS_DIR}/{infile}")
            #print(f" xxx infile = {infile} ")
            input_yaml = extract_yaml(INFILE)

            if 'CONFIG' in input_yaml:
                CONFIG     = input_yaml['CONFIG']
            else:
                CONFIG     = input_yaml  # for create_cov input

            outdir = None 
            if 'OUTDIR' in CONFIG :
                outdir = CONFIG['OUTDIR']
            if 'LOGDIR' in CONFIG :
                outdir = CONFIG['LOGDIR']

            #sys.exit(f" xxx infile={infile} -> outdir={outdir}")
            set_list += (f"{outdir} ")

        outdir_list.append(set_list)
    return outdir_list
    # end parse_outdir

def run_create_cov(infile_list, outdir_list, INPUTS):

    # Created Jan 20 2022
    # Run create_covariance job(s) interactively.
    # Applies to any input file with CREATE_COV in the name.

    for infile,outdir in zip(infile_list,outdir_list) :
        print(f" Create cov with {infile}  -> {outdir}")
        sys.stdout.flush()

        if INPUTS.snana_dir :
            jobname = f"{INPUTS.snana_dir}/util/{CREATE_COV_JOB_NAME}"
        else:
            jobname = CREATE_COV_JOB_NAME

        cmd_list = [ jobname, infile ]
        ret = subprocess.run( cmd_list, 
                              cwd=SNANA_TESTS_DIR,
                              capture_output=True, text=True )
    # end run_create_cov

def run_submit(infile_list, outdir_list, INPUTS):

    # infile_list is space-separated list of input files
    # to launch with submit_batch_jobs
    SUBMIT_JOB_NAME = os.path.expandvars(INPUTS.jobname)
    snana_dir       = INPUTS.snana_dir

    done_file_list = []
    for infile,outdir in zip(infile_list,outdir_list) :
        merge_file = f"{SNANA_TESTS_DIR}/{outdir}/{MERGE_LOG_FILE}"
        info_file  = f"{SNANA_TESTS_DIR}/{outdir}/{SUBMIT_INFO_FILE}"
        done_file  = f"{SNANA_TESTS_DIR}/{outdir}/{ALL_DONE_FILE}"      
        done_file_list.append(done_file)
        print(f" submit {infile}  -> {outdir}")
        sys.stdout.flush()

        arg_list = [ infile ]
        if snana_dir != SNANA_DIR:  
            arg_list = [ infile, f"--snana_dir", f"{snana_dir}" ]

        cmd_list = [ SUBMIT_JOB_NAME] + arg_list
        ret = subprocess.run( cmd_list, 
                              cwd=SNANA_TESTS_DIR,
                              capture_output=True, text=True )
        
        check_file_exists(merge_file)
        check_file_exists(info_file)
        time.sleep(3)

    # - - - - - - - 
    # wait for done files
    NDONE_EXPECT = len(done_file_list)
    NDONE_FIND   = 0
    while NDONE_FIND < NDONE_EXPECT :
        NDONE_FIND = 0
        time.sleep(20)
        for done_file in done_file_list :
            found = check_done_file(done_file)
            if found : NDONE_FIND += 1

        print(f"\t found {NDONE_FIND} of {NDONE_EXPECT} {ALL_DONE_FILE} files")
        sys.stdout.flush()

    # check for SUCCESS or FAIL

    print(f"\t Finished set with {STRING_SUCCESS}" )
    sys.stdout.flush()

    # end run_submit

# ===================================
# ============== MAIN ===============
# ===================================

if __name__ == "__main__":

    # parse input arguments
    INPUTS = parse_args()

    SUBMIT_INFO = {}

    # read list of input files 
    print(f" Read input files from :\n\t {INPUTS.list_file}")
    config = extract_yaml(INPUTS.list_file)
    SUBMIT_INFO.update( {'config' : config} )
    infile_submit_list = config[KEY_SUBMIT_LIST]
    nset = 0
    for infile_set in infile_submit_list:
        nset += 1
        if nset <= INPUTS.nrowskip : continue
        print(f"   inFile set:  {infile_set}")
        sys.stdout.flush() 
       
    outdir_submit_list = get_outdir_list(config)
    SUBMIT_INFO.update( { 'outdir_submit_list' : outdir_submit_list } )
    for outdir_set in outdir_submit_list:
        print(f"   outdir set: {outdir_set}")
        sys.stdout.flush()

    print(f"\n# - - - - - - - - - - - - - - - - - - \n"); 
    print(f" Output subDirs under \n  {SNANA_TESTS_DIR}\n")
    sys.stdout.flush()

    nset =  0
    for infile_set,outdir_set in zip(infile_submit_list,outdir_submit_list) :
        nset += 1
        if nset <= INPUTS.nrowskip : continue
        infile_list  = infile_set.split()
        outdir_list  = outdir_set.split()
        if PREFIX_CREATE_COV in infile_list[0]:
            run_create_cov(infile_list, outdir_list, INPUTS)
        else:
            #submit, wait for ALL.DONE
            run_submit(infile_list, outdir_list, INPUTS)  

    # - - - - 
    msg = (f"\n Done. " \
           f"All jobs reported {STRING_SUCCESS} in {ALL_DONE_FILE}\n")
    sys.exit(msg)

    # === END ===

    
