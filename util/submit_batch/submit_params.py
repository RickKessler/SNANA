# ==============================================
# Created July 2020 by R.Kessler & S. Hinton
#
#  Constant parameters & names for submit script.
#
#  Default memory -> 4GB (was 2GB)
#
#  Mar 3 2023: improve -H TRAIN_SALT2/3
#
# ==============================================

import os, sys
import datetime
import time
import getpass

#from colorama import Fore, Style

# start with flags that should be switched to command-line args
NCPU_MERGE_DISTRIBUTE  = 10000  # default: use all CPUs to merge
#NCPU_MERGE_DISTRIBUTE  = 0  # 0 -> merge only with CPU=0 (no conflict issues)

# debug feature: copy each MERGE.LOG to MERGE.LOG_{Nsec}
KEEP_EVERY_MERGELOG = False

ECHO_TIME = f"`date +%Y-%m-%d` `date +%H:%M:%S`"

# --fast option prescales by this factor
FASTFAC   = 10    # for --fast
FASTFAC2 = 100    # for --faster

# - - - - - -
SHELL            = os.environ['SHELL']

SNANA_DIR        = os.environ['SNANA_DIR']  # must have SNANA_DIR
try:
    SNDATA_ROOT  = os.environ['SNDATA_ROOT']
else:
    pass  # allow classes to run without SNDATA_ROOT (Aug 28 2025)


# generic program types to control batch flow
PROGRAM_TYPE_SIM    = "SIM"    # simulation
PROGRAM_TYPE_LCFIT  = "FIT"    # light curve fit (e.g., SALT2, PSNID, ...)
PROGRAM_TYPE_BBC    = "BBC"    # BEAMS with bias corrections
PROGRAM_TYPE_COSMOFIT  = "COSMOFIT"   # cosmology fitter

# default program names ... can be changed by user
PROGRAM_NAME_SIM            =  "snlc_sim.exe"
PROGRAM_NAME_LCFIT          =  "snlc_fit.exe" 
PROGRAM_NAME_BBC            =  "SALT2mu.exe"
PROGRAM_NAME_COVMAT         =  "create_covariance.py"
PROGRAM_NAME_WFIT           =  "wfit.exe"
PROGRAM_NAME_COMBINE_FITRES = "combine_fitres.exe"
PROGRAM_NAME_SMP            = "campari.py"
PROGRAM_NAME_FIRECROWN      =  None  
PROGRAM_NAME_MKDATA         =  "makeDataFiles.sh"
PROGRAM_NAME_UNKNOWN        =  "UNKNOWN"          # must be specified by JOBNAME key

# utilities
PROGRAM_NAME_CAT           = "sntable_cat.py"    # utility used by BBC class
PROGRAM_NAME_WAIT_FOR_FILE = "wait_for_file.py"  # util to wait for file and optionally check for string inside

SUBMIT_MODE_BATCH = "BATCH"
SUBMIT_MODE_SSH   = "SSH"
SBATCH_COMMAND    = 'sbatch'

# define subDir for batch scripts
SUBDIR_SCRIPTS_SIM    = ""
SUBDIR_SCRIPTS_LCFIT  = "SPLIT_JOBS_LCFIT"
SUBDIR_SCRIPTS_BBC    = "SCRIPTS_BBCFIT"
SUBDIR_SCRIPTS_COVMAT = "SCRIPTS_COVMAT"
SUBDIR_SCRIPTS_COSMOFIT   = "SCRIPTS_COSMOFIT"
SUBDIR_SCRIPTS_TRAIN  = "SCRIPTS_TRAIN"
SUBDIR_SCRIPTS_MKDATA = "SCRIPTS_MKDATA"
SUBDIR_SCRIPTS_COMBINE = ""
SUBDIR_OUTPUT_TRAIN   = "OUTPUT_TRAIN"
SUBDIR_CALIB_TRAIN    = "CALIB_TRAIN"
SUBDIR_SCRIPTS_SMP    = "SCRIPTS_SMP"    # July 8 2025
SUBDIR_MISC           = "misc"           # used by train_SALT3

MODEL_SNIa  = "SNIa"
MODEL_NONIa = "NONIa"

HOSTNAME = os.uname()[1].split('.')[0]

BATCH_MEM_DEFAULT      = "4GB"       # default mem is 4GB (9.24.2022) 
BATCH_WALLTIME_DEFAULT = '24:00:00'  # default wall time is 24hr
BATCH_MAXJOB_DEFAULT   = 500         # max number of jobs allowed in queue
BATCH_NTHREADS_DEFAULT = 1           # number of threads per job 08/apr/2022

USERNAME = getpass.getuser()
USER4    = USERNAME[0:4]

CWD      = os.getcwd()

today = datetime.date.today()

seconds_since_midnight = int(time.time() - time.mktime(today.timetuple()))

# define current time (e.g, 2022-03-21 04:54:12) and use rsplit
# to remove decimal places for seconds
time_submit_start = str(datetime.datetime.now()).rsplit('.',1)[0]

SUFFIX_FITRES = "FITRES"
SUFFIX_M0DIF  = "M0DIF"
SUFFIX_COV    = "COV"

# define monitor files
MERGE_LOG_FILE         = "MERGE.LOG"
MERGE_wfit_LOG_FILE    = "MERGE_wfit.LOG"  # only for BBC's wfit afterburner
SUBMIT_INFO_FILE  = "SUBMIT.INFO"
TABLE_SPLIT       = "SPLIT"  # yaml table in MERGE.LOG
TABLE_MERGE       = "MERGE"  # yaml table in MERGE.LOG
TABLE_EXTRA       = "EXTRA"  # do global overwrite on this name

MERGE_MODE_DEFAULT     = "DEFAULT"     # default
MERGE_MODE_SKIP        = "SKIP"        # if --nomerge
MERGE_MODE_BACKGROUND  = "BACKGROUND"  # if --merge_background

# True -> uses 'set -e' in each CPU*.CMD script to stop all
# all future processing upon any merge failures.
# BEWARE: doesn't work, and not clear we want this feature.
STOP_ALL_ON_MERGE_ERROR = False

SNANA_ABORT_STRING    = "FATAL ERROR ABORT"
FAIL_SUMMARY_FILE     = "FAIL_SUMMARY.LOG"

# FAIL_INDEX_XXX are used to select FAIL_MODE and FAIL_COMMENT (above)
FAIL_INDEX_ABORT    = 1  # explicit abort from code
FAIL_INDEX_ZERO     = 2  # zero events in output
FAIL_INDEX_BAD      = 3  # bad output detected (e.g.. crazy fit params)
FAIL_INDEX_UNKNOWN  = 4  # undetectable error such as crash 

FAIL_MODE_STRINGS  = [ 'TOTAL', 'ABORT', 'ZERO_EVENTS', 'BAD_OUTPUT', 'UNKNOWN' ]

FAIL_MODE_COMMENTS = [ '', \
                       'tail LOG file for more info', \
                       'check cuts, genRanges, etc...', \
                       'band output', \
                       'segFault? diskQuota? wallTime?'   ]


if len(FAIL_MODE_STRINGS) != len(FAIL_MODE_COMMENTS) :
    sys.exit(f"\n ERROR: sizes of FAIL_MODE_STRINGS & FAIL_MODE_COMMENTS are different. \n" \
             f"\t Fix it in submit_params.py !!!")

CMD_wildcard        = "CPU*.CMD"  # need to read all CMD files
JOB_SUFFIX_TAR_LIST  = [ 'YAML', 'DONE', 'LOG'  ]

# strings for DONE files
STRING_SUCCESS = "SUCCESS"
STRING_FAIL    = "FAIL"
STRING_STOP    = "STOP"

TMAX_EXE_WAIT_ABORT = 600  # max wait for binary exe to exist (in case of make)

KEY_END_YAML = "#END_YAML"

# either of these keys allowed in CONFIG
CONFIG_KEYLIST_DONE_FILE = [ 'DONE_STAMP', 'DONE_STAMP_FILE' ]
DEFAULT_DONE_FILE = "ALL.DONE"  # default if DONE_STAMP not in CONFIG

CONFIG_KEYNAME_ENV_REQUIRE  = "ENV_REQUIRE"  # name of key with required ENV
ENV_SNANA_SETUP_COMMAND     = "SNANA_SETUP_COMMAND"
ENV_SNANA_IMAGE_DOCKER      = "SNANA_IMAGE_DOCKER"

CONFIG_KEYLIST_SNANA_LOGIN_SETUP = \
        [ 'SNANA_LOGIN_SETUP', 'SNANA_SSHLOGIN_SETUP' ]

# lok file for merge process
BUSY_FILE_PREFIX = "BUSY_MERGE_CPU"
BUSY_FILE_SUFFIX = "LOCK"
BACKUP_PREFIX    = "BACKUP"

# define processing states
COLNUM_MERGE_STATE = 0  # first colmun of any MERGE table must be STATE
SUBMIT_STATE_WAIT = "WAIT"
SUBMIT_STATE_RUN  = "RUN "
SUBMIT_STATE_DONE = "DONE"
SUBMIT_STATE_FAIL = "FAIL"
SUBMIT_STATE_BUSY = "BUSY"

# xxx SUBMIT_STATE_WAIT = Fore.YELLOW + Style.BRIGHT + "WAIT" + Style.RESET_ALL

arg_check_abort  = "check_abort"
arg_kill_on_fail = "kill_on_fail"
 
KEY_DOCANA_START = 'DOCUMENTATION' 
KEY_DOCANA_END   = 'DOCUMENTATION_END'

# column ids for FITOPT_LIST written by fit job and read by BBC
COLNUM_FITOPT_NUM   = 0   # e.g., FITOPT001
COLNUM_FITOPT_LABEL = 1   # optional user label
COLNUM_FITOPT_ARG   = 2   # command-line args for fit job

# abort on any of these obsolete CONFIG keys (list includes sim, fit, bbc)
# Obsolete Dictionary includes comment printed to screen.
COMMENT_NOMORE_SALT2mu = "No more SALT2mu from FIT-input; use BBC input."
COMMENT_MAYBE_LATER    = "Might add this feature later."
COMMENT_NOT_NEEDED     = "No longer needed or relevant."
OBSOLETE_CONFIG_KEYS = \
{
    'SALT2mu_INFILE'            : COMMENT_NOMORE_SALT2mu ,
    'SALT2mu_SIMVERSION_INPUT'  : COMMENT_NOMORE_SALT2mu ,
    'SALT2mu_BIASCOR_PATH'      : COMMENT_NOMORE_SALT2mu ,
    'SALT2mu_CCPRIOR_PATH'      : COMMENT_NOMORE_SALT2mu ,
    'DO_FITOPT000'              : COMMENT_MAYBE_LATER ,
    'DOSKIP_DUPLICATE_SIMJOBS'  : COMMENT_MAYBE_LATER ,
    'VERSION_AFTERBURNER'       : COMMENT_MAYBE_LATER ,
    'VERSION+FITOPT'            : COMMENT_MAYBE_LATER ,
    'PLOTOPT'                   : COMMENT_MAYBE_LATER ,
    'CONVERT_SIMGEN_DUMP'       : COMMENT_NOT_NEEDED ,
    'DELAY_SUBMIT'              : COMMENT_NOT_NEEDED ,
    'H2ROOT_FLAG'               : COMMENT_NOT_NEEDED ,
    'MIN_SNANA_VERSION'         : COMMENT_NOT_NEEDED ,
    'TOPDIR_OVERRIDE'           : COMMENT_NOT_NEEDED ,
    'TOTAL_WAIT_ABORT'          : COMMENT_NOT_NEEDED ,
    'OUTDIR_OVERRIDE'           : "Use OUTDIR key instead (same key as in FIT input)" ,
    'GZIP_FLAG'                 : "gzip automatic; see CLEANUP_FLAG to NOT gzip",
    'APPEND_FITRES'             : "see APPEND_TABLE_VARLIST with -H FIT" ,
    'APPEND_TABLE_TEXT'         : "see APPEND_TABLE_VARLIST with -H FIT" ,
    'FITRES_COMBINE_FILE'       : "see APPEND_TABLE_TEXTFILE with -H FIT" ,
    'DUMMY'                     : "no comma here"
}


# === END ===
