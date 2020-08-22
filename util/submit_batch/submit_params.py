# ==============================================
# Created July 2020 by R.Kessler & S. Hinton
#
#  Constant parameters & names for submit script.
#  Giant HELP_CONFIG per task are at the bottom.
#
# ==============================================

import os
import datetime
import time
import getpass

# start with flags that should be switched to command-line args
#NCPU_MERGE_DISTRIBUTE  = 10000  # default: use all CPUs to merge
NCPU_MERGE_DISTRIBUTE  = 0  # 0 -> merge only with CPU=0 (no conflict issues)

# debug feature: copy each MERGE.LOG to MERGE.LOG_{Nsec}
KEEP_EVERY_MERGELOG = False  

# --fast option prescales by this factor
FASTFAC = 10

# - - - - - - 
SNANA_DIR        = os.environ['SNANA_DIR']
SNDATA_ROOT      = os.environ['SNDATA_ROOT']
SHELL            = os.environ['SHELL']

# generic program types to control batch flow
PROGRAM_TYPE_SIM  = "SIM"  # simulation
PROGRAM_TYPE_FIT  = "FIT"  # light curve fit (e.g., SALT2, PSNID, ...)
PROGRAM_TYPE_BBC  = "BBC"  # BEAMS with bias corrections

# default program names ... can be changed by user
PROGRAM_NAME_SIM  =  "snlc_sim.exe"
PROGRAM_NAME_FIT  =  "snlc_fit.exe"
PROGRAM_NAME_BBC  =  "SALT2mu.exe"

SUBMIT_MODE_BATCH = "BATCH"
SUBMIT_MODE_SSH   = "SSH"

# define subDir for batch scripts
SUBDIR_SCRIPTS_SIM = "" 
SUBDIR_SCRIPTS_FIT = "SPLIT_JOBS_LCFIT"
SUBDIR_SCRIPTS_BBC = "SCRIPTS_BBCFIT"

MODEL_SNIa  = "SNIa"
MODEL_NONIa = "NONIa"

HOSTNAME = os.uname()[1].split('.')[0]

MEMORY_DEFAULT = "2000"  # default memory request for batch jobs

USERNAME = getpass.getuser()
USER4    = USERNAME[0:4]

CWD      = os.getcwd()

today = datetime.date.today()

seconds_since_midnight = int(time.time() - time.mktime(today.timetuple()))

SUFFIX_FITRES = "FITRES"
SUFFIX_M0DIF  = "M0DIF"

# define monitor files
MERGE_LOG_FILE    = "MERGE.LOG"                                      
SUBMIT_INFO_FILE  = "SUBMIT.INFO"
TABLE_SPLIT       = "SPLIT"  # yaml table in MERGE.LOG
TABLE_MERGE       = "MERGE"  # yaml table in MERGE.LOG

# True -> uses 'set -e' in each CPU*.CMD script to stop all
# all future processing upon any merge failures.
# BEWARE: doesn't work, and not clear we want this feature.
STOP_ALL_ON_MERGE_ERROR = False  

SNANA_ABORT_STRING    = "FATAL ERROR ABORT"
FAIL_SUMMARY_FILE     = "FAIL_SUMMARY.LOG"

FAIL_MODE_STRINGS  = [ 'TOTAL', 'ABORT', 'ZERO_EVENTS', 'UNKNOWN' ]
FAIL_MODE_COMMENTS = [ '', \
                       'tail LOG file for more info', \
                       'check cuts, genRanges, etc...', \
                       'segFault? diskQuota? wallTime?'   ]

# FAIL_INDEX_XXX are used to select FAIL_MODE and FAIL_COMMENT (above)
FAIL_INDEX_ABORT    = 1
FAIL_INDEX_ZERO     = 2
FAIL_INDEX_UNKNOWN  = 3
CMD_wildcard        = "CPU*.CMD"  # need to read all CMD files

# strings for DONE files
STRING_SUCCESS = "SUCCESS"
STRING_FAIL    = "FAIL"
STRING_STOP    = "STOP"

# either of these keys allowed in CONFIG
CONFIG_KEYLIST_DONE_FILE = [ 'DONE_STAMP', 'DONE_STAMP_FILE' ] 
DEFAULT_DONE_FILE = "ALL.DONE"  # default if DONE_STAMP not in CONFIG

# lok file for merge process
BUSY_FILE_PREFIX = "BUSY_MERGE_CPU"
BUSY_FILE_SUFFIX = "LOCK"

# define processing states
COLNUM_MERGE_STATE = 0  # first colmun of any MERGE table must be STATE
SUBMIT_STATE_WAIT = "WAIT"
SUBMIT_STATE_RUN  = "RUN "
SUBMIT_STATE_DONE = "DONE"
SUBMIT_STATE_FAIL = "FAIL"
SUBMIT_STATE_BUSY = "BUSY"  # ?? 

# column ids for FITOPT_LIST written by fit job and read by BBC
COLNUM_FITOPT_NUM   = 0   # e.g., FITOPT001
COLNUM_FITOPT_LABEL = 1   # optional user label
COLNUM_FITOPT_ARG   = 2   # command-line args for fit job

# abort on any of these obsolete CONFIG keys 
#  (both sim and fit included in same list)
OBSOLETE_CONFIG_KEYS = \
[ 'TOPDIR_OVERRIDE', 'DOSKIP_DUPLICATE_SIMJOBS', 'CONVERT_SIMGEN_DUMP',
  'SALT2mu_INFILE', 'SALT2mu_SIMVERSION_INPUT', 'SALT2mu_BIASCOR_PATH',
  'SALT2mu_CCPRIOR_PATH', 'DO_FITOPT000', 'DELAY_SUBMIT', 'GZIP_FLAG',
  'H2ROOT_FLAG', 'APPEND_FITRES', 'APPEND_TABLE_TEXT', 'FITRES_COMBINE_FILE', 
  'MIN_SNANA_VERSION', 'VERSION_AFTERBURNER', 'PLOTOPT' ]
  
# ================================================
#   HELP_CONFIG


HELP_CONFIG_GENERIC = f"""
CONFIG:
  BATCH_INFO: sbatch [batch_template_file]  [n_core] 
     or 
  NODELIST: [node1] [node2] ...  # for ssh 
  
  # optional switch from using default $SNANA_DIR/bin/snlc_sim.exe
  JOBNAME: $MY_PATH/snlc_sim.exe 

  # optional memory request (default is 2 GB)
  BATCH_MEM: 4000     # 4GB (e..g, extra mem for big SIMSED models)

  # optional switch from using default $SNANA_DIR/bin/snlc_fit.exe
  JOBNAME: $MY_PATH/snlc_fit.exe 

  # default ALL.DONE is created under OUTDIR; here can specify
  # an optional/additionl done file anywhere
  DONE_STAMP_FILE: $MYPATH/PIPE_STAGE4.DONE
"""

HELP_CONFIG_SIM =  f"""
  ***** HELP/MENU for Simulation YAML Input ***** 

  """ +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  # option to re-route data files (default is $SNDATA_ROOT/SIM)
  PATH_SNDATA_SIM:  $SCRATCH_SIMDIR 

  RANSEED_REPEAT: 4 34212  # split into 4 jobs, then merge
       or
  RANSEED_CHANGE: 4 34212  # split into 4 jobs, do NOT merge

  FORMAT_MASK: 48    # +32=FITS, +16=random CID (TEXT not allowed)
  RESET_CIDOFF: 2    # unique CID among all GENVERSIONs
  NGEN_UNIT:   0.5   # 0.5 x NGENTOT_LC computed from
  RANGE(z,PKMJD) \n\t\t\t and SOLID_ANGLE
    (if no NGEN_UNIT, use NGENTOT_LC from sim-input or from GENOPT)
  GENPREFIX:    DES  # out_file name prefix (please keep it short)
                     # and suffix to default SIMLOGS_[GENPREFIX] 
  CLEANUP_FLAG:  0   # turn off default cleanup (for debug)

  SIMGEN_INFILE_SNIa:  # default SNIa input file(s) for all GENVERSIONs
  - SIMGEN_SNIa-SALT2+G10.input
  - etc ...            # can add more Ia models (e.g., C11, BS20..)
  SIMGEN_INFILE_NONIa: # default NONIa input files for all GENVERSIONs
  - SIMGEN_SNIbc-Templates.INPUT
  - SIMGEN_SNII-NMF.INPUT
  - etc ...            # can add more NONIa models

# NOTES:
#   If no NGEN_UNIT, use NGENTOT_LC from sim-input or from GENOPT
#   Valid SNIa  input file keys: SIMGEN_INFILE_[SNIa,SN1a,Ia,1a]
#   Valid NONIa input file keys: SIMGEN_INFILE_[NONIa,NON1a]

GENVERSION_LIST:
  - GENVERSION: DES_TEST0     # use SIMGEN_INFILE with no changes

  - GENVERSION: DES_TEST1     # apply GENOPT changes below
    GENOPT:
      GENPEAK_SALT2x1:   0.84
      GENSIGMA_SALT2x1:  1.4  0.2 
    SIMGEN_INFILE_NONIa:
    - SIMGEN_SNII-NMF_OVERRIDE.INPUT    # override default NON1a INFILEs above

  - GENVERSION: DES_TEST2
    GENOPT(SNIa):                  # apply only to SNIa jobs 
      GENPEAK_SALT2c:   -0.054
      GENSIGMA_SALT2c:   0.042  0.105
      GENTAU_AV:  0.18
    GENOPT(NONIa):                 # apply only to NONIa jobs
      GENTAU_AV:  0.35
    GENOPT(Iax):                   # apply only to SIMGEN_SNIax.INPUT
      GENTAU_AV:  0.55
    SIMGEN_INFILE_Ia:
    - SIMGEN_SNIa-SALT2+C11.input      # override default SNIa INFILE above
    SIMGEN_INFILE_NONIa:               # override NONIa INFILEs with these 3
    - SIMGEN_SNIbc-Templates_OVERRIDE.INPUT
    - SIMGEN_SNII-NMF_OVERRIDE.INPUT
    - SIMGEN_SNIax.INPUT

GENOPT_GLOBAL:   # OPTIONAL commands applied to all GENVERSIONs
   KEY1: ARG1
   KEY2: ARG2
   etc ...
# NOTE: INPUT_INCLUDE_FILE in GENOPT_GLOBAL is used for NGEN_UNIT and to
#       verify required keys ... but INPUT_INCLUDE_FILE inside GENOPT
#       is not used for these init functions.

"""


HELP_CONFIG_FIT = f"""    
   ***** HELP/MENU for LightCurveFit YAML Input *****

    All YAML input must go above &SNLCINP nameList block.

  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  OUTDIR:  [outdir]              # all output goes here
  VERSION:
  - MY_DATA    # in $SNDATA_ROOT/lcmerge, or PRIVATE_DATA_PATH 
  - MY_SIMDATA 
  - MY_SIMBIASCOR_*     # wildcard allowed
  - etc ... 
  FITOPT:
  - /ZPCAL/    MAGOBS_SHIFT_ZP g .01  # optional ZPCAL label for other codes
  - /ZPCAL/    MAGOBS_SHIFT_ZP r .01
  - /ZPCAL/    MAGOBS_SHIFT_ZP i .01 
  - /ZPCAL/    MAGOBS_SHIFT_ZP z .01 
  - /RETRAIN/  FITMODEL_NAME  SALT2.retrain1 
  - /RETRAIN/  FITMODEL_NAME  SALT2.retrain2 
  - /RETRAIN/  FITMODEL_NAME  SALT2.retrain3 
  - CUTWIN_SNRMAX 6 999  \t\t\t# no label needed
  - CUTWIN_SNRMAX 6 999  CUTWIN_SNRMAX2 4 999  # multiple options allowed   
  - USE_MINOS
  - FITOPT000     # no fit; ln -s FITOPT000.[SUFFIX] FITOPT011.[SUFFIX]
  - FITOPT000     # another synLink for FITOPT012
  - etc ... 

# Sym Link Notes for FITOPT000: this feature is useful for systematics
# with multiple surveys. For example above, FITOPT011 and FITOP012 could 
# be calibration variatios for a different survey, so here the sym link
# uses the default LCFIT (FITOPT000) without wasting CPU.

  # optional append variables from root/hbook into FITRES-TEXT table.
  #  (in old split_and_fit script, this key was APPEND_TABLE_TEXT)
  # To see full list of varables to append,
  #     stable_dump.pl <hbook_or_root_file> FITRES
   APPEND_TABLE_VARLIST: SNRMAX_g SNRMAX_r SNRMAX_i SNRMAX_z

  # optional append variables from external file into FITRES-TEXT table.
  #  (in old split_and_fit script, this key was FITRES_COMBINE_FILE)
  # E.g., supplement list of host properties, v_pec, etc ...
  APPEND_TABLE_TEXTFILE:  APPEND_THIS_FILE.FITRES

  # debug options to force failure in table-merge:
  FORCE_MERGE_TABLE_MISSING(HBOOK):  force missing HBOOK merge
  - DES_TEST1_FITOPT001
  - DES_TEST2_FITOPT001
  FORCE_MERGE_TABLE_CORRUPT(ROOT): # force corrupt ROOT file
  - DES_TEST3_FITOPT002

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# * Namelists below are used by snana.exe, psnid.exe, snlc_fit.exe.
# * Supplemental input files (e.g., KCOR_FILE) are copied to
#     {SUBDIR_SCRIPTS_FIT} if path is not included in file name.
# * HFILE_OUT and ROOTFILE_OUT are logical flags for batch script.

  &SNLCINP
   ! input for snana.exe, psnid.exe, snlc_fit.exe
    HFILE_OUT    = 'XYZ.HBOOK' # any name is logical flag for batch job
    ROOTFILE_OUT = 'XYZ.ROOT'  # any name is logical flag for batch job
   ! TEXT format is used regardless of whether TEXTFILE_PREFIX is defined
  &END

  &FITINP
     ! input for snlc_fit.exe
  &END
"""

HELP_CONFIG_BBC = f"""
    ***** HELP/MENU for BBC YAML Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  INPDIR+: 
  - dirSurvey1  # LCFIT OUTDIR for Survey1
  - dirSurvey2  # LCFIT OUTDIR for Survey2
  - dirSurvey3  # LCFIT OUTDIR for Survey3
    etc ...

  OUTDIR:   [outdir]   # all output goes here

  # if there are multiple versions per INPDIR, method is to "IGNORE"
  # substrings and then match version names. In example below, 
  # TEST1_DES is combined with TEST1_LOWZ,  TEST2_DES is matched to 
  # TEST2_LOWZ, etc... If there is 1 and only 1 version per INPDIR, 
  # there is no need for this string-match key.
  STRINGMATCH_IGNORE:   _DES  _LOWZ 
    
  # BBC variations
  MUOPT: 
  - p1=0.2 p2=3.3 
  - redchi2_tol=.02
  - sig1=0.14
  - simfile_biascor=[something_else]

  # Option to run wfit (fast, but ancient) as merge process. Useful
  # for quick cosmology cross-checks. To get help on argList, 
  # run "wfit.exe" with no args. Output cosmology fit params are 
  # in YAML format.
  WFITMUDIF_OPT: <argList>

  # process M independent random sum-samples; useful to compare RMS vs. errors.
  # Be careful that every VERSION+FITOPT+MUOPT is divided into NSPLITRAN jobs.
  NSPLITRAN: <nsplitran>
#END_YAML

"""

HELP_CONFIG = { 
    'SIM' : HELP_CONFIG_SIM,
    'FIT' : HELP_CONFIG_FIT,
    'BBC' : HELP_CONFIG_BBC
}

# === END ===
