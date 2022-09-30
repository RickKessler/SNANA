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
#
# - - - - - - - - - -

#import os
import sys, yaml, argparse, subprocess, logging
import submit_util      as util
import submit_translate as tr

from   submit_params      import *
from   submit_prog_sim    import Simulation

try:  # hack allows missing f90nml for makeDataFiles env (Oct 2021)
    from   submit_prog_lcfit  import LightCurveFit
except Exception as e:
    print(f" WARNING: could not import LightCurveFit")
    pass

from   submit_prog_bbc    import BBC
from   submit_prog_covmat import create_covmat
from   submit_prog_wfit   import wFit
from   submit_train_SALT2 import train_SALT2
from   submit_train_SALT3 import train_SALT3
from   submit_makeDataFiles import MakeDataFiles
from   argparse import Namespace

# =====================================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP with input file config(s); then exit"
    parser.add_argument("-H", "--HELP", help=msg, default=None, type=str,
                        choices = ["SIM", "LCFIT", "BBC", "WFIT",
                                   "TRAIN_SALT2", "TRAIN_SALT3",
                                   "TRANSLATE", "MERGE", "AIZ" ])
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

    msg = "+=1 -> new input file has REFAC_ prefix; " + \
          "+=2 -> old input file has LEGACY_ prefix ; " + \
          "+=4 -> continue submit with new file. "
    parser.add_argument('--opt_translate' , help=msg, type=int, default=1 )

    msg = "abort on missing DOCANA keys in maps & libraries"
    parser.add_argument("--require_docana", help=msg, action="store_true")

    msg = "run merge as background process instead of via batch"
    parser.add_argument("--merge_background", help=msg, action="store_true")

    msg = "DEBUG MODE: submit jobs, but skip merge process"
    parser.add_argument("--nomerge", help=msg, action="store_true")

    msg = (f"DEBUG MODE: reset (undo) merge process ")
    parser.add_argument("--merge_reset", help=msg, action="store_true")

    msg = (f"DEBUG MODE: developer flag to avoid conflicts. ")
    parser.add_argument("--devel_flag", help=msg, type=int, default=0 )

    msg = (f"DEBUG MODE: force crash in batch-prep ")
    parser.add_argument("--force_crash_prep", help=msg, action="store_true")
    msg = (f"DEBUG MODE: force crash in merge ")
    parser.add_argument("--force_crash_merge", help=msg, action="store_true")
    msg = (f"DEBUG MODE: force abort in merge ")
    parser.add_argument("--force_abort_merge", help=msg, action="store_true")

    msg = (f"DEBUG MODE: run codes from private snana_dir ")
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
        program_class = wFit         # fast cosmology fitter "wfit"

    elif "PATH_INPUT_TRAIN" in CONFIG :
        program_class = train_SALT2  # original snpca from J.Guy

    elif "JACOBIAN_MATRIX" in CONFIG :
        program_class = train_SALT2 # Patrick Armstrong - 17 Mar 22

    elif "SALT3_CONFIG_FILE" in CONFIG :
        program_class = train_SALT3  # saltshaker from D'Arcy & David

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
    opt_translate = args.opt_translate

    # abort if path is included in the input file name.
    if '/' in input_file :
        msgerr = []
        msgerr.append(f"Invalid input file: {input_file}")
        msgerr.append(f"because path not allowed as part of name.")
        msgerr.append(f"Must submit in same dir as input_file.")
        util.log_assert(False,msgerr)

    # check to translate legacy input
    args.input_file = check_legacy_input_file(input_file, opt_translate)

    #end check_input_file_name

def check_legacy_input_file(input_file, opt_translate):

    # if there is no 'CONFIG:' key, this is a legacy input file ;
    # translate using file-name convention based on user input
    # --opt_translate. See opt_translate details with
    #   submit_batch_jobs.sh -H TRANSLATE
    #
    # Function returns name of input file ... original of already
    # in correct YAML format, or translated.

    exit_always = (opt_translate & 8 ) > 0 # exit for legacy or refac file

    msgerr = []
    with open(input_file,"r") as f:
        flat_word_list=[word for line in f for word in line.split()]
        #f_read = f.read()

    if 'CONFIG:' in flat_word_list :
        # check for obsolete keys that are not translated
        for item in flat_word_list :
            key = item.rstrip(':')
            if key in OBSOLETE_CONFIG_KEYS :
                comment = OBSOLETE_CONFIG_KEYS[key]
                msgerr.append(f" Obsolete key '{key}' no longer valid.")
                msgerr.append(f" Comment: {comment}")
                util.log_assert(False,msgerr)

        if exit_always :
            sys.exit("\n Input file already translated; exit anyway.")

        return input_file    # file ok, do nothing.

    # - - - -  -

    #if opt_translate is None:  opt_translate = 1

    # prepare options
    rename_refac_file    = (opt_translate & 1 ) > 0
    rename_legacy_file   = (opt_translate & 2 ) > 0
    exit_after_translate = (opt_translate & 4 ) == 0 # default is to exit

    if '/' in input_file:
        msgerr.append(f"Will not translate input file in another directory.")
        msgerr.append(f"Recommend")
        msgerr.append(f"  cd {os.path.dirname(input_file)}")
        msgerr.append(f"  {os.path.basename(sys.argv[0])} " \
                      f"{os.path.basename(input_file)}")
        util.log_assert(False,msgerr)

    if rename_refac_file :
        legacy_input_file = input_file
        refac_input_file  = (f"REFAC_{input_file}")
    elif rename_legacy_file :
        if input_file[0:7] == 'LEGACY_' :  # don't add another LEGACY prefix
            legacy_input_file = input_file
            refac_input_file  = input_file[7:]
        else :
            legacy_input_file = (f"LEGACY_{input_file}")
            refac_input_file  = input_file
            cmd_mv = (f"mv {input_file} {legacy_input_file}")
            print(f" Save {input_file} as {legacy_input_file}")
            os.system(cmd_mv)
    else :
        msgerr.append(f" Must invalid opt_transate = {opt_translate} ")
        msgerr.append(f" Must have either ")
        msgerr.append(f"     opt_translate & 1 (rename refac file) or ")
        msgerr.append(f"     opt_translate & 2 (rename legacy file) ")
        util.log_assert(False,msgerr)

    msg_translate = (f"\n TRANSLATE LEGACY INPUT file for ")
    print(f" opt_translate = {opt_translate}")

    IS_SIM = False;   IS_FIT = False;  IS_BBC = False

    if  'GENVERSION:' in flat_word_list :  IS_SIM = True
    if  'VERSION:'    in flat_word_list :  IS_FIT = True
    if  '&SNLCINP'    in flat_word_list :  IS_FIT = True
    if  'u1='    in str(flat_word_list) :  IS_BBC = True

    if  IS_SIM :
        logging.info(f"{msg_translate} sim_SNmix.pl :")
        tr.SIM_legacy_to_refac( legacy_input_file, refac_input_file )

    elif IS_FIT :
        logging.info(f"{msg_translate} split_and_fit.pl :")
        tr.FIT_legacy_to_refac( legacy_input_file, refac_input_file )

    elif IS_BBC :
        logging.info(f"{msg_translate} SALT2mu_fit.pl: ")
        tr.BBC_legacy_to_refac( legacy_input_file, refac_input_file )
    #    program = BBC(config_yaml)
    else:
        msgerr = ['Unrecognized legacy input file:', input_file ]
        util.log_assert(False,msgerr)

    if exit_after_translate :
        sys.exit("\n Exit after input file translation.")


    return refac_input_file

    # end check_legacy_input_file

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

    # BBC
    util.find_and_remove(f"{SUBDIR_SCRIPTS_BBC}*")

    # SALT2 training
    util.find_and_remove(f"{SUBDIR_SCRIPTS_TRAIN}*")
    util.find_and_remove(f"{SUBDIR_CALIB_TRAIN}")
    util.find_and_remove(f"{SUBDIR_OUTPUT_TRAIN}")

    # SALT3 training
    util.find_and_remove(f"{SUBDIR_MISC}*")

    # end purge_old_submit_output

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

    # check input file: does it have a path? Does it need to be translated?
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
