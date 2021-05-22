# ==============================================
# Created July 2020 by R.Kessler & S. Hinton
#
#  Constant parameters & names for submit script.
#  Giant HELP_CONFIG per task are at the bottom.
#
#
# ==============================================

import os
import datetime
import time
import getpass

# start with flags that should be switched to command-line args
NCPU_MERGE_DISTRIBUTE  = 10000  # default: use all CPUs to merge
#NCPU_MERGE_DISTRIBUTE  = 0  # 0 -> merge only with CPU=0 (no conflict issues)

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
PROGRAM_NAME_SIM   =  "snlc_sim.exe"
PROGRAM_NAME_FIT   =  "snlc_fit.exe"
PROGRAM_NAME_BBC   =  "SALT2mu.exe"
PROGRAM_NAME_UNKNOWN =  "UNKNOWN"     # must be specified by JOBNAME key

SUBMIT_MODE_BATCH = "BATCH"
SUBMIT_MODE_SSH   = "SSH"

# define subDir for batch scripts
SUBDIR_SCRIPTS_SIM   = "" 
SUBDIR_SCRIPTS_FIT   = "SPLIT_JOBS_LCFIT"
SUBDIR_SCRIPTS_BBC   = "SCRIPTS_BBCFIT"
SUBDIR_SCRIPTS_TRAIN = "SCRIPTS_TRAIN"
SUBDIR_OUTPUT_TRAIN  = "OUTPUT_TRAIN"
SUBDIR_CALIB_TRAIN   = "CALIB_TRAIN"

MODEL_SNIa  = "SNIa"
MODEL_NONIa = "NONIa"

HOSTNAME = os.uname()[1].split('.')[0]

BATCH_MEM_DEFAULT      = "2000"      # default memory request is 2GB
BATCH_WALLTIME_DEFAULT = '24:00:00'  # default wall time is 24hr
BATCH_MAXJOB_DEFAULT   = 500         # max number of jobs allowed in queue

USERNAME = getpass.getuser()
USER4    = USERNAME[0:4]

CWD      = os.getcwd()

today = datetime.date.today()

seconds_since_midnight = int(time.time() - time.mktime(today.timetuple()))

SUFFIX_FITRES = "FITRES"
SUFFIX_M0DIF  = "M0DIF"
SUFFIX_COV    = "COV"

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

KEY_END_YAML = "#END_YAML"

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
SUBMIT_STATE_BUSY = "BUSY"

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



# ================================================
#   HELP_CONFIG


HELP_CONFIG_GENERIC = f"""
CONFIG:
  BATCH_INFO: sbatch [batch_template_file]  [n_core] 
     or 
  NODELIST: [node1] [node2] ...  # for ssh 

  # optional memory request (default is 2 GB)
  BATCH_MEM: 4000     # 4GB (e.g., extra mem for big SIMSED models)

  # optional max walltime request (default is 24hr)
  BATCH_WALLTIME: '1:00:00'  # 1hr max wall time

  # default ALL.DONE is created under OUTDIR; here can specify
  # an optional/additionl done file anywhere
  DONE_STAMP_FILE: $MYPATH/PIPE_STAGE4.DONE
"""


HELP_TRANSLATE = f"""
          TRANSLATING LEGACY INPUT FILES 

The 'LEGACY' input files for [sim_SNmix, split_and_fit, SALT2mu_fit] will 
not work with submit_batch_jobs.py, and therefore submit_batch_jobs includes 
an automatic translation of the input file. Command line option
     --opt_translate <opt>
controls the file-name convention, and also whether to exit or continue after 
translation. Note that opt_translate is a bit mask. LEGACY input files are 
automatically detected by the lack of a 'CONFIG:' key. If a CONFIG key exists,
opt_translate is ignored.

  opt_translate +=1 ->
    This default behavior produces a translated input file with name
    REFAC_[input_file]. The original input file is not modified.

  opt_translate +=2 -> 
    The original input file is saved as LEGACY_[input_file]; the translated 
    input file has the original name. If the original input file already has 
    a 'LEGACY_' prefix, the file name is not modified and the translated input
    file has the 'LEGACY_'  prefix removed.
    Example 1: input_file = abc.input is saved as LEGACY_abc.input;
               translated input file is abc.input
    Example 2: input_file = LEGACY_abc.input is not modified;
               translated input file is abc.input.

  opt_translate +=4 ->
    continue running submit_batch_jobs using translated input file.

  opt_translate +=8 ->
    always exit, regardless of whether input file is legacy or not.

Setting opt_translate to 1 or 2 results in translation followed by exiting 
submit_batch_jobs. This option enables visual inspection of translated input 
file before launching batch jobs. 

Setting opt_tranlate = 5 (1+4) or 6 (2+4) results in translation that is 
immediately followed by executation of the batch script (no questions asked). 
This option enables pipelines to run without interruption.

If the input file is already in the correct YAML format, opt_translate is 
ignored; therefore it is safe to always include an opt_translate argument.

  """


HELP_CONFIG_SIM =  f"""
  ***** HELP/MENU for Simulation YAML Input ***** 

  """ +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  # optional switch from using default $SNANA_DIR/bin/snlc_sim.exe
  JOBNAME: $MY_PATH/snlc_sim.exe 

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
                     # and default log dir is SIMLOGS_[GENPREFIX] 
  LOGDIR:   MY_LOGS  # override default SIMLOGS_[GENPREFIX]
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
  # optional switch from using default $SNANA_DIR/bin/snlc_fit.exe
  JOBNAME: $MY_PATH/snlc_fit.exe 

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
  - /NOREJECT_misc/  CUTWIN_MJD 56550 57040    # not in 2nd-iter BBC reject
  - CUTWIN_SNRMAX 6 999                        # no label needed
  - CUTWIN_SNRMAX 6 999  CUTWIN_SNRMAX2 4 999  # multiple options allowed   
  - USE_MINOS
  - FITOPT000     # no fit; ln -s FITOPT000.[SUFFIX] FITOPT011.[SUFFIX]
  - FITOPT000     # another synLink for FITOPT012
  - etc ... 

  # Sym Link Notes for FITOPT000: this feature is useful for systematics
  # with multiple surveys. For example above, FITOPT011 and FITOP012 could 
  # be calibration variatios for a different survey, so here the sym link
  # uses the default LCFIT (FITOPT000) without wasting CPU.
  # Labels with NOREJECT are used by BBC to exlcude these tests in 
  # determining reject list in 2nd BBC fit iteration.

  # Option to use events from FITOPT000 in all FITOPTs ...
  # except those with NOREJECT in label.
  OPT_SNCID_LIST: 1
     or
  OPT_SNCID_LIST: 3  # same, but also use PKMJDINI from FITOPT000
  
  # optional append variables from root/hbook into FITRES-TEXT table.
  #  (in old split_and_fit script, this key was APPEND_TABLE_TEXT)
  # To see full list of varables to append,
  #     stable_dump.pl <hbook_or_root_file> FITRES
   APPEND_TABLE_VARLIST: SNRMAX_g SNRMAX_r SNRMAX_i SNRMAX_z

  # optional append variables from external file(s) into FITRES-TEXT table.
  #  (in old split_and_fit script, this key was FITRES_COMBINE_FILE)
  # E.g., supplement list of host properties, v_pec, etc ...
  APPEND_TABLE_TEXTFILE:  APPEND_THIS_FILE.FITRES
  APPEND_TABLE_TEXTFILE:  APPEND1.FITRES,APPEND2.FITRES,APPEND3.FITRES
  APPEND_TABLE_TEXTFILE:  APPEND1.FITRES APPEND2.FITRES APPEND3.FITRES
  #  (accepts comma-sep or space sep list of files)

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

# If the merge process fails for silly reason (e.g, insufficiency memory,
# mis-typed variable name in APPEND_TABLE_VARLIST), or you want to add
# more variables to APPEND_TABLE_VARLIST, the merge process can be
# repeated interactively without re-doing the LC fits:
   0. visually check that MERGE.LOG is yaml compliant
   1. submit_batch_jobs.sh  <inputFile>  --merge_reset
   2. submit_batch_jobs.sh <inputFile> -M
      (be patient, and beware that merge repeat may be buggy)

"""

HELP_CONFIG_BBC = f"""
    ***** HELP/MENU for BBC YAML Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  # optional switch from using default $SNANA_DIR/bin/SALT2mu.exe
  JOBNAME: $MY_PATH/SALT2mu.exe 

  INPDIR+: 
  - dirSurvey1  # LCFIT OUTDIR for Survey1
  - dirSurvey2  # LCFIT OUTDIR for Survey2
  - dirSurvey3  # LCFIT OUTDIR for Survey3
    etc ...

       or
  INPDIR+: None  # flag to ignore INPDIR+; instead, use datafile= argument

  OUTDIR:   [outdir]   # all output goes here

  # if there are multiple versions per INPDIR, method is to "IGNORE"
  # substrings and then match version names. In example below, 
  # TEST1_DES is combined with TEST1_LOWZ,  TEST2_DES is matched to 
  # TEST2_LOWZ, etc... If there is 1 and only 1 version per INPDIR, 
  # there is no need for this string-match key.
  STRINGMATCH_IGNORE:   _DES  _LOWZ 
    
  STRING_VERSION_IGNORE:  ABC DEF  # ignore versions with these strings

  # BBC variations (for each VERSION and each FITOPT). Note that the
  # NOREJECT label excludes this MUOPT from defining reject.list.
  MUOPT: 
  - /SYST_ab/   p1=0.2 p2=3.3 
  - /SYST_tol/  redchi2_tol=.02
  - /SYST_sig1/ sig1=0.14
  - /NOREJECT/  simfile_biascor=[something_else]
  
  # Option to run wfit (fast, but ancient) as merge process. Useful
  # for quick cosmology cross-checks. To get help on argList, 
  # run "wfit.exe" with no args. Output cosmology fit params are 
  # in YAML format.
  WFITMUDIF_OPT: <argList>

  # Process subset of FITOPT x MUOPT matrix. Examples are
  FITOPTxMUOPT: 0+0   # process only FITOPT=000 or  MUOPT=000
  FITOPTxMUOPT: 2+3   # process only FITOPT=002 or  MUOPT=003
  FITOPTxMUOPT: 0&0   # process only FITOPT=000 and MUOPT=000
  FITOPTxMUOPT: 2&3   # process only FITOPT=002 and MUOPT=003

  # and to specify multiple options,
  FITOPTxMUOPT: 
  - 0+0   # process FITOPT=000 or MUOPT=000 ...
  - 2&4   # and include FITOPT002 x MUOPT004
  - DUMP  # optional flag to dump explicit matrix for FITOPTxMUOPT

  # For first two examples, the number of BBC jobs per version is
  # NFITOPT + NMUOPT + 1. For next 2 examples, just 1 BBC job per version.
  # This key cannot be configured to mimic command line args 
  # --ignore_fitopt or --ignore_muopt.  However, "FITOPTxMUOPT: 0&0" is 
  # equivalent to setting both with "--ignore_fitopt --ignore_muopt"

  # process independent random sum-samples; useful to compare RMS vs. errors.
  # Be careful that every VERSION+FITOPT+MUOPT is divided into NSPLITRAN jobs.
  NSPLITRAN: <nsplitran>

  # Optional FITOPT map of SURVEY x FITOPT(LCFIT) for each output FITOPT(BBC).
  # For each FITOPT(BBC), this map gives instructions for which FITRES file 
  # (from LCFIT task) to catenate from each INPDIR. This map is designed to
  # be created automatically by Pippin, but there may be cases where users 
  # create this map manually, or with a different script. The order of
  # LCFIT-FITOPTs follows the INPDIR order. With this map feature, links to 
  # FITOPT000 are NOT needed in LCFIT task. Batch init process for BBC writes
  # FITOPT_OUT_LIST (yaml block) to SUBMIT.INFO file with summary of labels 
  # and args for each FITOPT(BBC). Beware that LCFIT FITOPT numbers do not 
  # align with BBC's FITOPT_OUT_LIST numbers.
  FITOPT_MAP:
    SURVEY_LIST:  LOWZ  DES  PS1              # human-readable table header
    FITOPT000: FITOPT000 FITOPT000 FITOPT000  # global default
    FITOPT001: FITOPT001 FITOPT001 FITOPT001  # global (e.g., MWEBV_SCALE)
    FITOPT002: FITOPT002 FITOPT000 FITOPT000  # change only LOWZ (e.g., calib)
    FITOPT003: FITOPT003 FITOPT000 FITOPT000  # change only LOWZ
    FITOPT004: FITOPT000 FITOPT002 FITOPT000  # change only DES
    FITOPT005: FITOPT000 FITOPT003 FITOPT000  # change only DES
    FITOPT006: FITOPT000 FITOPT000 FITOPT002  # change only PS1
    FITOPT007: FITOPT000 FITOPT000 FITOPT003  # change only PS1

#END_YAML

"""


HELP_CONFIG_TRAIN_SALT2 = f"""
    ***** HELP/MENU for TRAIN_SALT2 YAML Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  # Must specify name of training script (because it's outside SNANA)
  JOBNAME: [train_script_name]

  # input data files and config for snpca. Includes survey.yaml, which
  # maps snpca surveys/instruments into SNANA survey/bandpasses.
  PATH_INPUT_TRAIN: [path]  

  # input Instrument and MagSys (aka SALTPATH)
  PATH_INPUT_CALIB: [path] 

  # TRAINOPT args specify calibration systematics per band or group of bands.
  # An independent training is done for each TRAINOPT argument.
  # SHIFTLIST_FILE is a file containing a list of MAGSHIFT and WAVESHIFT
  # keys; <CR> are stripped so that contents can be distributed among 
  # multiple lines for human readability. The explicit MAGSHIFT and
  # WAVESHIFT keys are intended for linear perturbations to measure
  # derivatives for systematics; the SHIFTLIST_FILE feature is intended
  # for a random calibration offset in every band.
  # PATH_INPUT_CALIB key specifies a different calibration directory.

  TRAINOPT:
  - MAGSHIFT  SDSS  g 0.01
  - MAGSHIFT  SDSS  g,z 0.01,-0.01    MAGSHIFT CfA2 B 0.01
  - WAVESHIFT CfA3  r,i 10,10         MAGSHIFT CfA3 U .01
  - SHIFTLIST_FILE  shifts_01.dat
  - SHIFTLIST_FILE  shifts_02.dat
  - SHIFTLIST_FILE  shifts_03.dat
  - PATH_INPUT_CALIB  $PATH/calib_different

  OUTDIR:   [outdir]   # all output goes here

 # The TRAINOPT-calibration shifts in the training are propagated to 
 # SNANA's light curve fitting via MAGSHIFT and WAVESHIFT keys written
 # to the SALT2.INFO file for each SALT2.MODELnnn directory. The mapping
 # "snpca survey/instrument -> SNANA survey/passbands" is contained in 
 # [PATH_INPUT_TRAIN]/survey.yaml.

"""


HELP_CONFIG_TRAIN_SALT3 = f"""
    ***** HELP/MENU for TRAIN_SALT3 YAML Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  # Must specify name of training code (because it's outside SNANA)
  JOBNAME: [train_script_name]

  CONDA_DEFAULT_ENV:  saltshaker  # abort if $CONDA_DEFAULT_ENV != saltshaker

  # top-level config input files
  SALT3_CONFIG_FILE:  <fileName>  # top-level input for SALTshaker code

  # global command-line options for all TRAINOPTs below
  TRAINOPT_GLOBAL: <list of options>
     e.g.,
  TRAINOPT_GLOBAL: --resume_from_outputdir $SNTRAIN_ROOT/SALT3/SALT3.K21


  # TRAINOPT args specify calibration systematics per band or group of bands.
  # An independent training is done for each TRAINOPT argument.
  # SHIFTLIST_FILE is a file containing a list of MAGSHIFT and WAVESHIFT
  # keys; <CR> are stripped so that contents can be distributed among 
  # multiple lines for human readability. The explicit MAGSHIFT and
  # WAVESHIFT keys are intended for linear perturbations to measure
  # derivatives for systematics; the SHIFTLIST_FILE feature is intended
  # for a random calibration offset in every band.
  # PATH_INPUT_CALIB key specifies a different calibration directory.

  TRAINOPT:
  - MAGSHIFT  SDSS  g 0.01
  - MAGSHIFT  SDSS  g,z 0.01,-0.01    MAGSHIFT CfA2 B 0.01
  - LAMSHIFT  CfA3  r,i 10,10         MAGSHIFT CfA3 U .01
  - SHIFTLIST_FILE  shifts_01.dat
  - SHIFTLIST_FILE  shifts_02.dat
  - SHIFTLIST_FILE  shifts_03.dat

  OUTDIR:   [outdir]   # all output goes here

 # The TRAINOPT-calibration shifts in the training are propagated to 
 # SNANA's light curve fitting via MAGSHIFT and LAMSHIFT keys written
 # to the SALT3.INFO file in each SALT3.MODELnnn directory. 

"""


HELP_MERGE = f"""
          MERGE LOGIC

While there are no merge options, this help section may be useful in case
debugging is needed for a merge process that doesn't finish properly.
'MERGE' refers to tasks run after a science job (SciJob). For SIM, the merge
process combines sim data files from multiple split jobs into a single
data version. For FIT, the merge process combines tables from the split
jobs into a single table (per FITOPT). For BBC, there is no merging since
a BBC job cannot be split among multiple cores. Merge tasks also include
organization of science files (e.g. FITRES tables), such as moving them 
to a more appropriate location outside of the messy script-directory where 
jobs had run.

The general merge strategy is that each CPU launches a merge process
after each SciJob has finished. Thus each merge process must wake
up, figure out what (if any) action is needed, and take action. The
advantages of this strategy are 
  1:) distribute merge task load among multiple CPUs.
  2:) merge tasks are done in batch job, not on login node.
  3:) if batch jobs are killed, no lingering background jobs on login node.
Difficulties are 
  1:( collect adequate information.
  2:( avoid conflict when mutliple merge tasks are launched simultaneously.
  3:( ensure merge process runs when all is finished.

In the CPU*.CMD files, each SciJob is followed by a merge
task as follows
   python submit_batch_jobs.py <inputFile> -m -t 21014 --cpunum 0 

where -m tells the batch script to function as a merge task (instead of 
submit task), the -t argument is a time stamp (Nsec since midnight)
to verify against time stamp in OUTDIR, and --cpunum is used to label
BUSY files. 

The difficulties are addressed as follows:

1:(
Each merge task collects information from two required files under OUTDIR
that are created before batch jobs are submitted: 
     SUBMIT.INFO    # fixed info, never changes
     MERGE.LOG      # state of each SciJob: WAIT, RUN, DONE, FAIL
 
The merge process analyzes these two files and decides what action is
necessary. After taking action, MERGE.LOG is updated so that the next 
merge process knows not to repeat already finished merge tasks. When all 
SciJobs and merge tasks have finished, a final "cleanup" task is run to 
do things like compress files and create a summary file.

2:(
Merge conflicts are avoided using busy files named
   BUSY_MERGE_CPU[cpunum].LOCK
A merge process exits immediately if a BUSY*.LOCK file exists; otherwise 
it creates a BUSY file to lock out other merge processes. The BUSY*LOCK
file is removed after the merge process has finished and updated MERGE.LOG.
It is possible that multiple BUSY*LOCK files are created simultaneously. 
To handle this situation, after a BUSY*LOCK file is created the merge 
process waits a few seconds and then checks again for a list of BUSY*LOCK 
files. If more than 1 exists, only the first in the sorted list remains 
active, while the others exit. For example, suppose
       BUSY_MERGE_CPU0002.LOCK  
       BUSY_MERGE_CPU0006.LOCK  
both exist; CPU0006 merge process will exit while CPU0002 merge process 
remains to carry out its merge tasks. To track the occurance of multiple 
BUSY*LOCK files, "grep simultaneous CPU*.LOG"


3:(
The last issue is that all merge tasks can exit before finishing. For example,
consider 30 jobs on 30 CPUs, and suppose that CPU-4 SciJob finishes first. 
Next, suppose that the merge process performs a few tasks, and during these 
tasks all other CPUs finish and exit because of the BUSY_MERGE_CPU0004.LOCK 
file. In summary, a single LOCK file can result in no remaiming merge task 
when all SciJobs finish.

To ensure a final merge process after all SciJobs finish, the merge task 
after the last SciJob gets a -M argument instead of -m. The -M argument is 
an instruction to wait for all expected DONE files, and to wait for any 
remaining BUSY*LOCK files to clear. Beware that the last SciJob does not 
always run last due to different wait times in the batch queue. For example, 
consider 3 SciJobs submitted to 2 CPUs (CPU000, CPU001). CPU000 has SciJob 
1 & 3, while CPU001 has SciJob 2. If CPU000 runs before CPU001, the last 
SciJob (3) can finish long before SciJob 2. In this scenario, here is the
expected sequence of events:
  + CPU000 runs, while CPU001 waits in the queue.
  + SciJob 1 finishes on CPU000.
  + merge proc 1 (-m) runs after SciJob 1
  + SciJob 3 finishes on CPU000.
  + merge proc 3 (-M) sees only 2 DONE files, so does nothing while waiting 
    for 3rd DONE file.
  + CPU001 finally runs
  + SciJob 2 finishes on CPU001.
  + merge proc 2 (-m) runs arter SciJob 3 and runs merge tasks for everything,
    including unfinished tasks from CPU000. It skips cleanup because only
    the last merge proc can do cleanup.
  + waiting merge proc 3 (CPU000) sees all 3 DONE files, and also sees that
    all merge tasks are done. It runs only the cleanup to compress and 
    create summary file(s).


"""

HELP_AIZ = f"""
    LOGIC for "ABORT IF ZERO" (AIZ) 

There are many mistakes which can result in zero output events. To trap such
mistakes (or code bugs), the YAML output for each science job includes 
ABORT_IF_ZERO (AIZ), which is the number of events processed; AIZ=0 is a 
failure flag for the merge process.

However, care is needed to avoid aborting on low-stat jobs where Poisson 
fluctuations result in zero events for a subset of the split jobs. For example, 
consider a Kilonova simulation where there are 50 total events (after trigger), 
and the sim job is split over 50 cores. The average number per core is 1, 
meaning that about 1/3 of the split jobs will have zero events. In this case, 
there is no failure, but a naive check on AIZ=0 would result in a false failure.
On the other hand, if 50 SNIa jobs result in AIZ>10000 on 49 cores and AIZ=0
on 1 core, this should result in a failure.

The AIZ abort logic is as follows in the merge process. If any AIZ >= 30, the 
naive logic applies where any AIZ=0 flags a failure. If all AIZ < 30 (grep
for aiz_thresh), the low-stat Poisson regime is assumed and a subset of AIZ=0 
is allowed without flagging a failure. If all AIZ=0, failure is flagged.

Determining aiz_thresh:
A false failure occurs if the true <AIZ> is below aiz_thresh, at least one job
has AIZ=0 (with probabilithy P_0), and at least one job has AIZ >= aiz_thresh
(with prob P_thresh). For 'Nsplit' jobs, the false failure probability is
       P_FF = [1-(1-P_0)**Nsplit)] x [1-(1-P_thresh)**Nsplit)]

With Nsplit = 100, the table below shows P_FF values where the mean aiz is 
taken to be <aiz> = aiz_thresh/2:

    aiz_thresh <aiz>    P_0     P_thresh    P_FF
     -----------------------------------------------
       10        5     6.74E03  3.18E-2   0.47
       20       10     4.54E-5  3.45E-3   1.32E-3
       24       12     6.14E-6  1.47E-3   8.42E-5
       30       15     3.06E-7  4.18E-4   1.25E-6
     -----------------------------------------------

A reasonable choice is aiz_thresh=30 so that P_FF ~ E-6

"""

# - - - - - - - 
HELP_MENU = { 
    'SIM'         : HELP_CONFIG_SIM,
    'FIT'         : HELP_CONFIG_FIT,
    'BBC'         : HELP_CONFIG_BBC,
    'TRAIN_SALT2' : HELP_CONFIG_TRAIN_SALT2,
    'TRAIN_SALT3' : HELP_CONFIG_TRAIN_SALT3,
    'TRANSLATE'   : HELP_TRANSLATE,
    'MERGE'       : HELP_MERGE,
    'AIZ'         : HELP_AIZ     # ABORT_IF_ZERO
}

# === END ===
