#!/usr/bin/env python
#
# Created July 2020 by R.Kessler & S. Hinton  
#
# Oct 29 2020: add SALT2train framework
# Nov 24 2020: add --ncore arg
# Dec 17 2020: purge now works on train_SALT2 outputs
# Jan 22 2021: garbage above CONFIG is ignored.
# Jan 23 2021: begin adding train_SALT3
# May 24 2021: call submit_iter2()
# Aug 09 2021: add --snana_dir arg
# Oct 04 2021; add wfit class; maybe later can generalize to cosmofit?
# Oct 12 2021: implement --check_abort for lcfit
# Oct 20 2021: add makeDataFiles
# Dec 04 2021: new input --merge_background
# Feb 02 2022: add --faster arg to prescale by 100 (e.g., for WFD sim)
# Apr 08 2022: add --merge_force arg for sync_evt option
# Sep 30 2022: begin new create_covmat class (stat+syst covar matrix) 
# Oct 18 2022: rename wfit class to cosmifit (more general name)
# Feb 27 2023: undo hack that allowed missing f90nml for makeDataFiles
# May 19 2023: remove obsolete translate class
# May 22 2023: add train_BAYESN class
# Jan 17 2023: new option --cpu_sum
# May 07 2024: add heatmaps to purge list.
#
# - - - - - - - - - -

#import os
import sys, yaml, argparse, subprocess, logging
import submit_util      as util


from   submit_params      import *
from   submit_prog_sim    import Simulation
from   submit_prog_lcfit  import LightCurveFit


from   submit_prog_bbc      import BBC
from   submit_prog_covmat   import create_covmat
from   submit_prog_cosmofit import cosmofit
from   submit_train_SALT2   import train_SALT2
from   submit_train_SALT3   import train_SALT3
from   submit_train_BAYESN  import train_BAYESN
from   submit_makeDataFiles import MakeDataFiles
from   argparse import Namespace

# =====================================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP with input file config(s); then exit"
    parser.add_argument("-H", "--HELP", help=msg, default=None, type=str,
                        choices = ["SIM", "LCFIT", "BBC", "COVMAT", "COSMOFIT",
                                   "TRAIN_SALT2", "TRAIN_SALT3", "TRAIN_BAYESN",
                                   "MERGE", "AIZ" ])
    msg = "name of input file"
    parser.add_argument("input_file", help=msg, nargs="?", default=None)

    # misc user args

    msg = "Create & init outdir, but do NOT submit jobs"
    parser.add_argument("-n", "--nosubmit", help=msg, action="store_true")

    # - - - - -
    # change number of cores
    msg = "number of cores"
    parser.add_argument('--ncore', help=msg, type=int, default=None )
    # xxx mark parser.add_argument('--ncore', nargs='+', help=msg, type=int )

    msg = "override OUTDIR in config file"
    parser.add_argument('--outdir', help=msg, type=str, default=None )
    #parser.add_argument('--outdir', nargs='+', help=msg, type=str )

    # reduce processing
    msg = "process x10 fewer events for sim,fit,bbc (applies only to sim data)"
    parser.add_argument("--fast", help=msg, action="store_true")

    msg = "process x100 fewer events"
    parser.add_argument("--faster", help=msg, action="store_true")

    msg = "ignore FITOPT (LC & BBC fits)"
    parser.add_argument("--ignore_fitopt", help=msg, action="store_true")

    msg = "ignore MUOPT for BBC fit"
    parser.add_argument("--ignore_muopt", help=msg, action="store_true")

    # - - - -
    # purge files
    msg = "Use 'find' to locate and remove non-essential output."
    parser.add_argument("--purge", help=msg, action="store_true")

    msg = f"Diagnostic: Sum all CPU under current dir"
    parser.add_argument("--cpu_sum", help=msg, action="store_true")

    # - - -
    msg = "increase output verbosity (default=True)"
    parser.add_argument("-v", "--verbose", help=msg, action="store_true")

    # - - - -
    msg = "kill current jobs (requires input file as 1st arg)"
    parser.add_argument("-k", "--kill", help=msg, action="store_true")

    msg = "kill jobs if FAIL is detected"
    parser.add_argument("--kill_on_fail", help=msg, action="store_true")

    msg = "check for abort using interactive job for 300 events"
    parser.add_argument("--check_abort", help=msg, action="store_true")

    msg = "abort on missing DOCANA keys in maps & libraries"
    parser.add_argument("--require_docana", help=msg, action="store_true")

    msg = "run merge as background process instead of via batch"
    parser.add_argument("--merge_background", help=msg, action="store_true")

    # - - - - - DEBUG options - - - - 
    msg = "DEBUG MODE: submit jobs, but skip merge process"
    parser.add_argument("--nomerge", help=msg, action="store_true")
    msg = f"DEBUG MODE: reset (undo) merge process "
    parser.add_argument("--merge_reset", help=msg, action="store_true")
    msg = f"DEBUG MODE: developer flag to avoid conflicts. "
    parser.add_argument("--devel_flag", help=msg, type=int, default=0 )
    msg = f"DEBUG MODE: force crash in batch-prep "
    parser.add_argument("--force_crash_prep", help=msg, action="store_true")
    msg = f"DEBUG MODE: force crash in merge "
    parser.add_argument("--force_crash_merge", help=msg, action="store_true")
    msg = f"DEBUG MODE: force abort in merge "
    parser.add_argument("--force_abort_merge", help=msg, action="store_true")
    msg = f"DEBUG MODE: force garbage argument for this job id"  # not implemented yet ...
    parser.add_argument("--jobid_force_bad_arg", help=msg, type=int, default=None )

    msg = f"DEBUG MODE: run codes from private snana_dir "
    parser.add_argument("--snana_dir", help=msg, type=str, default=None )

    # args passed internally from command files
    msg = "INTERNAL:  merge process (if no BUSY file)"
    parser.add_argument("-m", "--merge", help=msg, action="store_true")

    msg = "INTERNAL: last merge process when all done files exist"
    parser.add_argument("-M", "--MERGE_LAST", help=msg, action="store_true")

    msg = "INTERNAL:  force merge process (wait for BUSY files to disappear)"
    parser.add_argument("--merge_force", help=msg, action="store_true")

    msg = "INTERNAL: time stamp (Nsec since midnight) to verify" \
           " merge process examines correct output_dir"
    parser.add_argument('-t', nargs='+', help=msg, type=int )

    msg = "INTERNAL: cpu number for BUSY-merge file-name"
    parser.add_argument('--cpunum', nargs='+', help=msg, type=int )

    msg = "INTERNAL:  2nd iteration"
    parser.add_argument("--iter2", help=msg, action="store_true")

    args = parser.parse_args()

    # internally set prescale arg
    args.prescale = 1
    if args.fast   : args.prescale = FASTFAC
    if args.faster : args.prescale = FASTFAC2

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return args

    # end get_args

def which_program_class(config):

    # check YAML/input_file keys to determine which program class
    # (sim,fit,bbc) will run

    program_class = None
    input_file    = config['args'].input_file
    merge_flag    = config_yaml['args'].merge_flag
    CONFIG        = config['CONFIG']
    if "GENVERSION_LIST" in config :
        program_class = Simulation

    elif "VERSION" in CONFIG :
        program_class = LightCurveFit # SALT2 LC fits

    elif "INPDIR+" in CONFIG :
        program_class = BBC          # Beams with Bias Corr (KS17)

    elif "WFITOPT" in CONFIG :
        program_class = cosmofit    # wfit ...

    elif "FIRECROWN_INPUT_FILE" in CONFIG :
        program_class = cosmofit    # firecrown/Cosmosis ...   
        
    elif "PATH_INPUT_TRAIN" in CONFIG :
        program_class = train_SALT2  # original snpca from J.Guy

    elif "JACOBIAN_MATRIX" in CONFIG :
        program_class = train_SALT2 # Patrick Armstrong - 17 Mar 22

    elif "SALT3_CONFIG_FILE" in CONFIG :
        program_class = train_SALT3  # saltshaker from D'Arcy & David

    elif "BAYESN_CONFIG_FILE" in CONFIG :
        program_class = train_BAYESN  # from Grayling, Thorp, Narayan, Mandel

    elif "MAKEDATAFILE_SOURCE" in CONFIG:
        program_class = MakeDataFiles

    elif "INPUT_COVMAT_FILE" in CONFIG:  # create_covariance, Sep 30 2022
        program_class = create_covmat

    else :
        sys.exit("\nERROR: Could not determine program_class")

    # keep quiet for merge process
    if not merge_flag :
        logging.info(f" Program class : {program_class.__name__} \n" )

    return program_class

def set_merge_flag(config):

    args = config['args']
    merge_flag = args.merge  or \
                 args.MERGE_LAST   or \
                 args.merge_force  or \
                 args.merge_reset

    set_cpunum0 = args.merge_reset or args.MERGE_LAST
    cpunum      = args.cpunum
    if cpunum is None and set_cpunum0 :
        args.cpunum = [ 0 ]

    if args.merge_background:
        args.nomerge = True

    return merge_flag
    # end set_merge_flag

def check_input_file_name(args):

    input_file    = args.input_file

    # abort if path is included in the input file name.
    if '/' in input_file :
        msgerr = []
        msgerr.append(f"Invalid input file: {input_file}")
        msgerr.append(f"because path not allowed as part of name.")
        msgerr.append(f"Must submit in same dir as input_file.")
        util.log_assert(False,msgerr)


    #end check_input_file_name

def print_submit_messages(config_yaml):

    # print final info to screen for user
    CONFIG = config_yaml['CONFIG']
    args   = config_yaml['args']

    logging.info(f" Done launching jobs. Sit back and relax.")

    if 'OUTDIR' in CONFIG :
        OUTDIR = CONFIG['OUTDIR']
        MERGE_LOG = (f"{OUTDIR}/{MERGE_LOG_FILE}")
        logging.info(f" Check status in {MERGE_LOG} ")

    if args.merge_background :
        logging.info(f" REMEMBER: merge is background process (view with ps -f).")
    elif args.nomerge :
        logging.info(f" REMEMBER: you disabled the merge process.")

    if args.prescale > 1 :
        logging.info(f" REMEMBER: fast option will process " \
                     f"1/{args.prescale} of request.")

    if args.force_crash_merge :
        logging.info(f" REMEMBER: there is a forced crash in MERGE process.")

    if args.force_abort_merge :
        logging.info(f" REMEMBER: there is a forced abort in MERGE process.")

    return

    # end print_submit_messages

def print_nosubmit_messages(config_yaml):
    # print final info to screen for user
    CONFIG = config_yaml['CONFIG']
    if 'OUTDIR' in CONFIG :
        logging.info(f"\n Check job preparation in {CONFIG['OUTDIR']}/ ")

    logging.info(f" Jobs NOT sumbitted. Bye Bye.")

    # end print_nosubmit_messages

def purge_old_submit_output():

    # LC fitting
    util.find_and_remove(f"{SUBDIR_SCRIPTS_LCFIT}*")
    util.find_and_remove(f"FITOPT*.LCPLOT*")
    util.find_and_remove(f"FITOPT*.HBOOK*")
    util.find_and_remove(f"FITOPT*.ROOT*")

    # Scone heatmaps (May 2024)
    util.find_and_remove(f"heatmaps*.tfrecord")
        
    # BBC
    util.find_and_remove(f"{SUBDIR_SCRIPTS_BBC}*")

    # SALT2 training
    util.find_and_remove(f"{SUBDIR_SCRIPTS_TRAIN}*")
    util.find_and_remove(f"{SUBDIR_CALIB_TRAIN}")
    util.find_and_remove(f"{SUBDIR_OUTPUT_TRAIN}")

    # SALT3 training
    util.find_and_remove(f"{SUBDIR_MISC}*")

    # end purge_old_submit_output

def print_cpu_sum():
    
    # Created Jan 17 2024 (Perlmutter down for maintenance)
    # Find all MERGE.LOG files in current pwd;
    # sum CPU by class; print CPU_SUM by class, and total CPU.
    # This task is strictly diaganostic to help track CPU needs.
    # This task does not submit or process jobs.

    SUMMARY_FILE_LIST = [ MERGE_LOG_FILE, MERGE_wfit_LOG_FILE, 'SCONE_SUMMARY.LOG' ] 

    logging.info(f"\n Sum CPU in {SUMMARY_FILE_LIST} files ... ")

    summary_log_list = []

    for summ_file in SUMMARY_FILE_LIST:
        cmd_find  = f"find . -name {summ_file} " + "-exec du -mc {} +"

        OPT_FIND = 'run'  # or 'check_output'
        if OPT_FIND == 'check_output' :
            find_list    = subprocess.check_output(cmd_find, shell=True)
            find_list    = (find_list.rstrip()).decode('utf-8')
        elif OPT_FIND == 'run' :
            find_list  = subprocess.run( cmd_find.split(), cwd='./' ,
                                         capture_output=True, text=True ).stdout

        if len(find_list) == 0 : continue

        # every other elment is file or dir   .xyz
        find_list  = find_list.split()[1::2]  
        find_list  = find_list[:-1]

        # remove last element for 'total'
        summary_log_list += find_list
        #print(f"\n xxx summary_log_list = {summary_log_list}")

    # - - - - 
    # define expected keys in MERGE.LOG
    KEY_PROGRAM_CLASS     = 'PROGRAM_CLASS'
    KEY_CPU_SUM           = 'CPU_SUM'

    # PROGRAM_CLASS key does not exist before about mid-Jan 2024,
    # so substitute "UnknownClass" in this case
    PROGRAM_CLASS_UNKNOWN = "UnknownClass"

    cpu_dict = {}
    cpu_sum_total = 0.0 
    NMISSING_KEY_CPU_SUM = 0
    N_SKIP               = 0

    SKIP_LIST = [ '5_MERGE', 'denied' ] # skip these duplicates

    for merge_log in summary_log_list:

        do_skip = False
        for str_skip in SKIP_LIST :
            if str_skip in merge_log: do_skip = True
        if do_skip:
            N_SKIP += 1
            continue

        yaml_contents ,comment_lines = util.read_merge_file(merge_log)

        logging.info(f" Read {merge_log}")
        if KEY_PROGRAM_CLASS in yaml_contents:
            program_class = yaml_contents[KEY_PROGRAM_CLASS]
        else:
            program_class = PROGRAM_CLASS_UNKNOWN

        if KEY_CPU_SUM in yaml_contents:
            cpu  =  float(yaml_contents[KEY_CPU_SUM])
        else:
            cpu = 0.0
            NMISSING_KEY_CPU_SUM += 1
            logging.info(f" **** WARNING **** Missing required {KEY_CPU_SUM} key.")

        cpu_sum_total += cpu
        logging.info(f"\t CPU = {cpu} hr for program class {program_class}")
        if program_class not in cpu_dict:
            cpu_dict[program_class] = 0.0
        cpu_dict[program_class] += cpu

    # sort by CPU in decreasing order
    cpu_dict_sorted = dict(sorted(cpu_dict.items(), key=lambda x:x[1], reverse=True ))
        
    CLASS_TOTAL = "*** TOTAL ***"
    cpu_dict_sorted[CLASS_TOTAL] = cpu_sum_total

    logging.info("")
    for program_class, cpu in cpu_dict_sorted.items():
        percent = 100.0 * (cpu / cpu_sum_total)
        logging.info(f"  CPU-sum of {program_class:16s} :  {cpu:8.2f} hr   " \
                     f"({percent:6.1f} %)")

    if NMISSING_KEY_CPU_SUM > 0 :
        logging.info(f"\n WARNING: {NMISSING_KEY_CPU_SUM} " \
                     f"{MERGE_LOG_FILE} files are missing {KEY_CPU_SUM} key.")

    if N_SKIP > 0 :
        logging.info(f"\n WARNING: {N_SKIP} skipped files containing {SKIP_LIST}")
        
    logging.info("")

    return
    # end print_cpu_sum

def print_HELP():
    see_me = (f" !!! ************************************************ !!!")
    print(f"\n{see_me}\n{see_me}\n{see_me}")
    HELP_upcase = (args.HELP).upper()
    print(f"{HELP_MENU[HELP_upcase]}")
    sys.exit(' Scroll up to see full HELP menu.\n Done: exiting Main.')

def run_merge_driver(program,args):
    try:
        program.merge_driver()
        if args.check_abort:
            exit(0)
        else:
            logging.info('  Done with merge process -> exit Main.')
            exit(0)
    except Exception as e:
        logging.exception(e, exc_info=True)
        cpunum   = args.cpunum[0]
        cpu_file = f"CPU{cpunum:04d}*.LOG"
        print(f"{e}")
        msg      = [ f"Check {cpu_file} for merge crash" ]
        program.log_assert(False, msg )
    
    # end run_merge_driver

# =============================================
if __name__ == "__main__":

    args  = get_args()
    store = util.setup_logging(args)

    # check option for long HELP menus
    if args.HELP : 
        print_HELP()

    # check option to "purge" un-needed files with linux find and rm;
    # removes tarred script-dirs, root & hbook files, etc ...
    if args.purge :
        purge_old_submit_output()
        sys.exit(' Done with purge: exiting Main.')

    if args.cpu_sum :
        print_cpu_sum()
        sys.exit(' Done with summing CPU times.')

    # check input file: does it have a path?
    check_input_file_name(args)

    # Here we know there's a CONFIG block, so read the YAML input
    config_yaml = util.extract_yaml(args.input_file, "CONFIG:", KEY_END_YAML)
    config_yaml['args'] = args  # store args here for convenience

    # set logical merge flag before running program_class
    config_yaml['args'].merge_flag   = set_merge_flag(config_yaml)

    logging.debug(config_yaml)  # ???

    # - - - - - -
    # determine which program class (sim, fit, bbc, train ...)
    program_class = which_program_class(config_yaml)
    args.program_class = f"{program_class.__name__}"

    # run the class
    program = program_class(config_yaml)  # calls __init only

    # - - - - - - - -
    # check merge options
    if args.merge_flag:
        run_merge_driver(program,args)

    # - - - - - -
    # check option to kill jobs
    if args.kill :
        program.kill_jobs()
        logging.info('  Done killing jobs -> exit Main.')
        exit(0) 

    # - - - - - -
    try:
        # create output dir
        program.create_output_dir()

        # prepare files, lists, program args
        program.submit_prepare_driver()

        # write .BATCH and .CMD scripts
        program.write_script_driver()

        # Create MERGE.LOG file with all jobs in WAIT state.
        # This file gets updated later by merge process.
        program.create_merge_file()

        # create static SUBMIT.INFO file for merge process
        # (unlike MERGE.LOG, SUBMIT.INFO never changes)
        program.create_info_file()

    except Exception as e:
        logging.exception(e, exc_info=True)
        msg    = [ e, "Crashed while preparing batch jobs.",
                   "Check Traceback" ]
        program.log_assert(False, msg )

    # - - - - - -
    if args.nosubmit :
        print_nosubmit_messages(config_yaml);
        program.submit_iter2()
        exit(0)

    # - - - - - - - - -
    program.launch_jobs() # submit via batch or ssh

    # Print warnings and errors to make sure they aren't missed
    store.print_warnings()
    store.print_errors()
    print_submit_messages(config_yaml) # final stuff for user to REMEMBER

    # check for iterative submit (e.g., sync events in BBC)
    program.submit_iter2()

    exit(0)    # bye bye.

# === END ===
