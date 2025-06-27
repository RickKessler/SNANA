
from submit_params    import *

# ================================================
#   HELP_CONFIG per task. Each HELP_CONFIG is a mini-manual that should
#   help users create/modify a config file for ther task.
#   HELP_CONFIG_GENERIC is common info for all tasks.
#   For new help menu, make sure to append HELP_MENU dictionary at bottom 
#   of this file.



HELP_CONFIG_GENERIC = f"""
CONFIG:   # generic info for any task

  BATCH_INFO: sbatch [batch_template_file]  [n_core]
     or
  NODELIST: [node1] [node2] ...  # for ssh

  # optional memory request (default is 2 GB)
  BATCH_MEM: 8GB     # e.g., extra mem for big SIMSED models

  # optional max walltime request (default is 24hr)
  BATCH_WALLTIME: '1:00:00'  # 1hr max wall time

  # option to force all jobs on single node
  BATCH_SINGLE_NODE: True

  # optional list of required ENVs (aborts if any ENV is not defined)
  ENV_REQUIRE: SNANA_LSST_SIM  LSST_STACK_VERSION

  # snana setup for SSH-login (not for batch)
  SNANA_SSHLOGIN_SETUP: <SNANA setup command>

  # default ALL.DONE is created under OUTDIR; here can specify
  # an optional/additionl done file anywhere
  DONE_STAMP_FILE: $MYPATH/PIPE_STAGE4.DONE

  # define code if not in path, or optional code switch from 
  # default code in path (e.g., for development)
  JOBNAME: $MY_PATH/MYCODE

  # ---------------------------------------------------------------
  #    Below are task-specific options 
  # ---------------------------------------------------------------
"""



HELP_CONFIG_SIM =  f"""
  ***** HELP/MENU for Simulation CONFIG-yaml Input *****

  """ +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""

  # option to re-route data files (default is $SNDATA_ROOT/SIM)
  PATH_SNDATA_SIM:  $SCRATCH_SIMDIR

  RANSEED_REPEAT: 4 34212  # split into 4 jobs, then merge
       or
  RANSEED_CHANGE: 4 34212  # split into 4 jobs, do NOT merge

  FORMAT_MASK: 48        # +32=FITS, +16=random CID (TEXT not allowed)
  RESET_CIDOFF: 1        # unique CID within each GENVERSION
  RESET_CIDOFF: 2        # unique CID among all GENVERSIONs
  CIDOFF_GLOBAL:  100000 # add this offset to all CIDOFF

  NGEN_UNIT:   0.5   # 0.5 x NGENTOT_LC computed from
  RANGE(z,PKMJD) \n\t\t\t and SOLID_ANGLE
    (if no NGEN_UNIT, use NGENTOT_LC from sim-input or from GENOPT)
  GENPREFIX:    DES  # out_file name prefix (please keep it short)
                     # and default log dir is SIMLOGS_[GENPREFIX]
  LOGDIR:   MY_LOGS  # override default SIMLOGS_[GENPREFIX]
  CLEANUP_FLAG:  0       # turn off default SIMLOGS cleanup (for debug)
  NOGZIP_DATA_FILES:     # turn off default gzip of output FITS data files
                         # (only if there is a very good reason)

  SIMGEN_INFILE_SNIa:  # default SNIa input file(s) for all GENVERSIONs
  - SIMGEN_SNIa-SALT2+G10.input
  - etc ...            # can add more Ia models (e.g., C11, BS20..)
  SIMGEN_INFILE_NONIa: # default NONIa input files for all GENVERSIONs
  - SIMGEN_SNIbc-Templates.INPUT
  - SIMGEN_SNII-NMF.INPUT
  - etc ...            # can add more NONIa models

  # for RANSEED_CHANGE, optional list of input-include files where each
  # include files has systematic shifts for one sim.
  INPUT_INCLUDE_FILE: <$PATH/myList>

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

# special option to use different sim code version; enables running on 
# multiple old code versons to identify large change.
  - GENVERSION: DES_TEST_SNANA_v11_04e
      GENOPT:  JOBNAME  /products/SNANA_v11_04e/bin/snlc_sim.exe
   
GENOPT_GLOBAL:   # OPTIONAL commands applied to all GENVERSIONs
   KEY1: ARG1
   KEY2: ARG2
   etc ...
# NOTE: INPUT_INCLUDE_FILE in GENOPT_GLOBAL is used for NGEN_UNIT and to
#       verify required keys ... but INPUT_INCLUDE_FILE inside GENOPT
#       is not used for these init functions.

"""


HELP_CONFIG_LCFIT = f"""
   ***** HELP/MENU for LightCurveFit CONFIG-yaml Input *****

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
  - /NOREJECT_misc/  CUTWIN_MJD 56550 57040    # not in 2nd-iter BBC reject
  - CUTWIN_SNRMAX 6 999                        # no label needed
  - CUTWIN_SNRMAX 6 999  CUTWIN_SNRMAX2 4 999  # multiple options allowed
  - USE_MINOS
  - FITOPT000     # no fit; ln -s FITOPT000.[SUFFIX] FITOPT011.[SUFFIX]
  - FITOPT000     # another synLink for FITOPT012
  - etc ...

  FITOPT_GLOBAL: <command-line arg list>  # global options for all FITOPTs

  # Version-dependent fit-options can be specified with
  FITOPT_GLOBAL(MY_DATA): <command-line data args>
  FITOPT_GLOBAL(SIM):     <command-line sim  args>
  
  # where the first set of fit-options applies only to VERSION MY_DATA 
  # (see VERSION: key above), and the 2nd set applies to VERSIONs containing
  # "SIM" in the name:  MY_SIMDATA and MY_SIMBIASCOR_* . The options in 
  # "FITOPT_GLOBAL:" are always included, regardless of whether or not 
  # version-dependent FITOPTs are provided.

  # Sym Link Notes for FITOPT000: this feature is useful for systematics
  # with multiple surveys. For example above, FITOPT011 and FITOP012 could
  # be calibration variatios for a different survey, so here the sym link
  # uses the default LCFIT (FITOPT000) without wasting CPU.
  # Labels with NOREJECT are used by BBC to exlcude these tests in
  # determining reject list in 2nd BBC fit iteration.

  # To repeat fits with multiple code versions (e.g., to identify code version 
  # causing a big change);
  FITOPT:
  - /v11_04e/  JOBNAME /products/SNANA_v11_04e/bin/snlc_fit.exe
  - /v11_04f/  JOBNAME /products/SNANA_v11_04f/bin/snlc_fit.exe
  - /v11_04g/  JOBNAME /products/SNANA_v11_04g/bin/snlc_fit.exe
  etc ...
  # will use the specified JOBNAME and also suppress the JOBNAME key and arg
  # from the argument list.

  # Option to use events from FITOPT000 in all FITOPTs ...
  # except those with NOREJECT in label.
  OPT_SNCID_LIST: 1         # use CIDs from FITOPT000 (original SNANA key)
  OPT_SNCID_LIST: 3         # use CIDs and use FITPAR as INIVAL
  OPT_SNCID_LIST: 7         # use CIDs and use FITPAR as INIVAL and Gauss prior
     or
  FLAG_USE_SAME_EVENTS: 1   # alternate key used by Pippin
  FLAG_USE_SAME_EVENTS: 3   # idem
  FLAG_USE_SAME_EVENTS: 7   # idem

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
# * Namelists below are used by snana.exe, psnid.exe, {PROGRAM_NAME_LCFIT}.
# * Supplemental input files (e.g., KCOR_FILE) are copied to
#     {SUBDIR_SCRIPTS_LCFIT} if path is not included in file name.
# * HFILE_OUT and ROOTFILE_OUT are logical flags for batch script.

  &SNLCINP
   ! input for snana.exe, psnid.exe, {PROGRAM_NAME_LCFIT}
    HFILE_OUT    = 'XYZ.HBOOK' # any name is logical flag for batch job
    ROOTFILE_OUT = 'XYZ.ROOT'  # any name is logical flag for batch job
   ! TEXT format is used regardless of whether TEXTFILE_PREFIX is defined
  &END

  &FITINP
     ! input for {PROGRAM_NAME_LCFIT}
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

HELP_CONFIG_COSMOFIT = f"""
   ***** HELP/MENU for wfit or Firecrown CONFIG-yaml Input *****

  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""

  # input directories for wfit are output by create_covariance.
  # Each dir contains hubble diagram file, covsys_* files, and INFO.YML
  INPDIR:
  - <dir0>
  - <dir1>
  - <dir2>
  - /global/mydir/output_test1*    # wildcards allowed
  - /global/mydir/output_test2*    # multiple times
  - etc ...

  # List wfit command-line options; each row is for a separate wfit job.
  # Beware that there is no default setting, so WFITOPT starts here at 0.
  # The /label/ are optional, and are included in the wfit-summary table
  # for human- and machine-readability.
  WFITOPT:
  - /OMpri/     -ompri 0.31 -dompri 0.01
  - /CMBpri/    -cmb_sim -sigma_Rcmb 0.007
  - /CMB+BAO/   -cmb_sim -sigma_Rcmb 0.007 -bao_sim
  - /w0wa+CMB/  -wa -wasteps 51 -w0steps 51 -cmb_sim -sigma_Rcmb 0.007
  - /w0wa+CMB/  -wa -wasteps 51 -w0steps 51 -outfile_chi2grid -outfile_resid 
#      (for outfile_xxx options, output files are stored in separate subdir)
      or
  FCOPT:
  - /label1/ <firecrown options> 
  - /label2/ <firecrown options> 
 
  COVOPT:  ALL NOSYS  # select subset of cov systematics options

# optional global wfit options appended to each WFITOPT above
  WFITOPT_GLOBAL:  -hsteps 61 -wsteps 101 -omsteps 81
       or
  FCOPT_GLOBAL:  <firecrown options>

# Default blind flag is set for data, but not for sim.
# These defaults can be changed with
  BLIND_DATA: False   # no blinding -> show results (check with your collaborators first!!!)
  BLIND_SIM:  True    # blind the sim results
# The -blind flag cannot be included in WFITOPT above,
# otherwise submit_batch will abort.

  OUTDIR:  [outdir]              # all output goes here

# Optional avg calculation of cosmology parameters or differences;
# error is computed as stddev/sqrt(Nsamples).
# X-Y means compute the average differences for w_X-w_Y and same for 
# other fitted cosmology parameters. Note that each entry is a sub-string 
# of the directory under 7_CREATE_COV,and sometimes /X is needed to ensure 
# uniqueness.
# 
  WFITAVG:  # also works with WEIGHT_AVG and FCAVG key
  - /BIN5YR_5YR_OnlyIa - /BIN5YR_5YR_IaCC    # w(wa) diff-avg
  - UNBIN5YR_5YR_OnlyIa - UNBIN5YR_5YR_IaCC  # w(wa) diff-avg
  - UNBIN5YR_5YR_IaCC                        # w(wa) avg

  USE_COVSYS_INV: False  # disable default use of already-inverted COVSYS
                         # and let cosmo-fit code invert COVSYS

"""

HELP_CONFIG_BBC = f"""
    ***** HELP/MENU for BBC  CONFIG-yaml  Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""

  INPDIR+:
  - dirSurvey1  # LCFIT OUTDIR for Survey1
  - dirSurvey2  # LCFIT OUTDIR for Survey2
  - dirSurvey3  # LCFIT OUTDIR for Survey3
    etc ...

       or
  INPDIR+: None  # flag to ignore INPDIR+; instead, use datafile= argument

  # to include external FITRES files (e.g., from a separate analysis)
  INPFILE+:
  - file1   # include this explicit FITRES file; no subdir structure
  - file2
    etc ...

  OUTDIR:   [outdir]   # all output goes here

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # if there are multiple versions per INPDIR, method is to "IGNORE"
  # substrings and then match version names. In example below,
  # TEST1_DES is combined with TEST1_LOWZ,  TEST2_DES is matched to
  # TEST2_LOWZ, etc... If there is 1 and only 1 version per INPDIR,
  # there is no need for this string-match key.
  STRINGMATCH_IGNORE:   _DES  _LOWZ

  STRING_VERSION_IGNORE:  ABC DEF  # ignore versions with these strings

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BBC variations (for each VERSION and each FITOPT). Note that the
  # NOREJECT label excludes this MUOPT from defining reject.list.
  # MUOPT000 corresponds to input config file with no alterations
  # and hence it is not specified here (i.e., removing MUOPT key
  # here results in MUOPT000 output)
  MUOPT:
  - /SYST_ab/   p1=0.2 p2=3.3                     # MUOPT001
  - /SYST_tol/  redchi2_tol=.02                   # MUOPT002
  - /SYST_sig1/ sig1=0.14                         # MUOPT003
  - /NOREJECT/  simfile_biascor=[something_else]  # MUOPT004

  # To run with older code versions, use JOBNAME key in MUOPT; e.g.
  MUOPT:
  - /v11_04e/  JOBNAME  /products/SNANA_v11_04e/bin/SALT2mu.exe
  - /v11_04f/  JOBNAME  /products/SNANA_v11_04f/bin/SALT2mu.exe
  etc ...

  # default label for MUOPT000 is None; to output a more informative label
  # (e.g., for downstream codes to parse);
  MUOPT000_LABEL:  MY_NOMINAL_LABEL

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Option to run wfit (fast, but ancient) as merge process. Useful
  # for quick cosmology cross-checks. To get help on argList,
  # run "wfit.exe" with no args. Output cosmology fit params are
  # in YAML format. Since there is no config file with default params,
  # WFITMUDIF_OPTs start at 0, not 1 as for MUOPTs.
  #
  WFITMUDIF_OPT: <argList>
  #   or
  WFITMUDIF_OPT:   # labels below are optional
  - /label0/  <argList0>         # wfit0
  - /lable1/  <argList1>         # wfit1
  - /label2/  <argList2>         # wfit2
  - etc ...

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 #       SAME-EVENT LOGIC for SYSTEMATICS
 # The BBC class reads SUBMIT.INFO from the LCFIT output and checks for
 #       OPT_SNCID_LIST > 0   or   FLAG_USE_SAME_EVENTS > 0
 # which is a flag to analyze common events from FITOPT000. In this case,
 # submit_batch automatically applies similar logic to BBC by submiting
 # two sequential iterations:
 #   iter1: analyze list of events passing cuts in FITOPT000
 #            (BiasCor rejects might still be different per FITOPT)
 #   iter2: analyze common subset found in all FITOPTs
 #            (accounts for biasCor rejects)
 # This process results in two sets of OUTDIRs: [OUTDIR]_ITER1 and [OUTDIR],
 # where the latter is for the 2nd (final) iteration. Beware that 2 x NCORE
 # jobs are launched; the 2nd set of jobs waits until the ALL.DONE file
 # appears from 1st set. FITOPT labels with NOREJECT (e.g., NOREJECT_TEST)
 # are analyzed with cuts and not used for same-event logic. To use this
 # feature with Pippin, put "FLAG_USE_SAME_EVENTS: 1" in the FITOPT.yml
 # file (same file where systematic variations are defined).

 FLAG_USE_SAME_EVENTS: 0  # disable sync-event flag passed from 2_LCFIT stage

 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 # Flag to monitor rejected events due to common-event requirement;
 # create self-documented BBC_REJECT_MONITOR.FITRES and add reject stats to BBC_FITPAR_SUMMARY
 WRITE_REJECT_MONITOR: True  

#END_YAML

"""

HELP_CONFIG_COVMAT = f"""
    ***** HELP/MENU for COVMAT CONFIG-yaml Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""

  # input file for {PROGRAM_NAME_COVMAT}
  INPUT_COVMAT_FILE:   CREATE_COVMAT_TEMPLATE.INPUT

  # define BBC outdirs; these override INPUT_DIR in the above input file.
  # Labels are required and are part of outdir folder names; submit_batch_jobs
  # aborts if a label is missing.
  BBC_OUTDIR:  
  - /5YR_IaCC/   $PIPPIN_OUTPUT/MYJOBS/6_BIASCOR/5YR_IaCC/output
  - /5YR_OnlyIa/ $PIPPIN_OUTPUT/MYJOBS/6_BIASCOR/5YR_OnlyIa/output

  # repeat job with these COVMATOPT options. While LCFIT and BBC have default 
  # FITOPT000, there is no default here and thus an option with no arguments is 
  # default. Each option must have label to use in outdir folder name; 
  # submit_batch_jobs aborts if label is missing.
  COVMATOPT:
  - /BINNED/                  # default with no args, but must have label
  - /UNBIN/      --unbinned
  - /REBIN_2x4/  --nbin_x1 2 --nbin_c 4
  - /REBIN_4x8/  --nbin_x1 4 --nbin_c 8

  # To test with fewer jobs, this optional key selects strings to match the
  # version-subdir names. In this example, 11 version-subdirs are selected;
  # 0003 and 0030-0039.
  VERSION_SUBSET:
  - OUTPUT_BBCFIT-0003  # select one version-subdir
  - OUTPUT_BBCFIT-003   # select OUTPUT_BBCFIT-0030 thru OUTPUT_BBCFIT-0039

  # direct all output here
  OUTDIR: 7_CREATE_COV

"""

HELP_TRAINOPT_GENERIC = f"""
  # The nominal training with salt-input files is output as TRAINOPT000.
  # For additional trainings with variations (e.g, calibration shifts,
  # training options), use the TRAINOPT key as illustrated by the example
  # below that produces TRAINOPT001 thru TRAINOPT007.  Remove TRAINOPT key
  # to produce only TRAINOP000.
  # SHIFTLIST_FILE is a file containing a list of MAGSHIFT and WAVESHIFT
  # keys; <CR> are stripped so that contents can be distributed among
  # multiple lines for human readability. The explicit MAGSHIFT and
  # WAVESHIFT keys are intended for linear perturbations to measure
  # derivatives for systematics; the SHIFTLIST_FILE feature is intended 
  # for a random calibration offset in every band.

  TRAINOPT: # survey  band shift
  - MAGSHIFT  SDSS  g 0.01
  - MAGSHIFT  SDSS  g,z 0.01,-0.01    MAGSHIFT CfA2 B 0.01
  - WAVESHIFT CfA3  r,i 10,10         MAGSHIFT CfA3 U .01
  - SHIFTLIST_FILE  shifts_01.dat 
  - SHIFTLIST_FILE  shifts_02.dat
  - SHIFTLIST_FILE  shifts_03.dat
  """


HELP_CONFIG_SMP = f"""
    ***** HELP/MENU for SMP YAML Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""
  SMP-specific info ... TO DO LATER ...
  etc ...
  """

HELP_CONFIG_TRAIN_SALT2 = f"""
    ***** HELP/MENU for TRAIN_SALT2 YAML Input *****
  """  +  (f"{HELP_CONFIG_GENERIC}") +  \
  f"""

  # input data files and config for snpca. Includes survey.yaml, which
  # maps snpca surveys/instruments into SNANA survey/bandpasses.
  PATH_INPUT_TRAIN: [path]

  # input Instrument and MagSys (aka SALTPATH)
  PATH_INPUT_CALIB: [path]

  {HELP_TRAINOPT_GENERIC}
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
  """  +  f"{HELP_CONFIG_GENERIC}" +  \
  f"""

  CONDA_DEFAULT_ENV:  saltshaker  # abort if $CONDA_DEFAULT_ENV != saltshaker

  # top-level config input files
  SALT3_CONFIG_FILE:  <fileName>  # top-level input for SALTshaker code

  # global command-line options for all TRAINOPTs below
  TRAINOPT_GLOBAL: <list of options>
     e.g.,
  TRAINOPT_GLOBAL: --resume_from_outputdir $SNTRAIN_ROOT/SALT3/SALT3.K21

  {HELP_TRAINOPT_GENERIC}

  OUTDIR:   [outdir]   # all output goes here

 # The TRAINOPT-calibration shifts in the training are propagated to
 # SNANA's light curve fitting via MAGSHIFT and LAMSHIFT keys written
 # to the SALT3.INFO file in each SALT3.MODELnnn directory.

"""


HELP_CONFIG_TRAIN_BAYESN = f"""
    ***** HELP/MENU for TRAIN_BAYESN YAML Input *****
  """  +  f"{HELP_CONFIG_GENERIC}" +  \
  f"""

  # top-level config input files
  BAYESN_CONFIG_FILE:  <fileName> 

  # global command-line options for all TRAINOPTs below
  TRAINOPT_GLOBAL: <list of options>

  {HELP_TRAINOPT_GENERIC}

  OUTDIR:   [outdir]   # all output goes here

 # The TRAINOPT-calibration shifts in the training are propagated to
 # SNANA's light curve fitting via MAGSHIFT and LAMSHIFT keys written
 # to the BAYESN.INFO file in each BAYESN.MODELnnn directory.

"""

HELP_CONFIG_MAKEDATAFILES = f"""
    ***** HELP/MENU for MAKEDATAFILES Input *****
  """  +  f"{HELP_CONFIG_GENERIC}" +  \
  f"""

  MAKEDATAFILE_SOURCE:  SNANA_FOLDER #  or DES or SIRAH_FOLDER or LSST_FASTDB

  MAKEDATAFILE_INPUTS:
  - <input1>  # data folder or db name
  - <input2>  

  MAKEDATAFILE_ARGS:  <optional args for makeDataFiles.sh>

  NEVT: 500  # process this many from each input (for test only)
  
  FIELD:         WFD
  NSPLITRAN:     5
  
  OUTDIR_SNANA: <outDir>

"""



HELP_MERGE = f"""
          MERGE LOGIC

While there are no explicit merge options, this help section may be useful for 
adding a new class, or debugging a merge process that doesn't finish properly.
'MERGE' refers to tasks run after a science job (SciJob). For SIM, the merge
process combines sim data files from multiple split jobs into a single
data version. For FIT, the merge process combines tables from the split
jobs into a single table (per FITOPT). For BBC, there is no combine since
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
   submit_batch_jobs.sh <inputFile> -m -t 21014 --cpunum 0

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
    'SIM'           : HELP_CONFIG_SIM,
    'LCFIT'         : HELP_CONFIG_LCFIT,
    'BBC'           : HELP_CONFIG_BBC,
    'COVMAT'        : HELP_CONFIG_COVMAT,
    'COSMOFIT'      : HELP_CONFIG_COSMOFIT,
    'SMP'           : HELP_CONFIG_SMP,
    'TRAIN_SALT2'   : HELP_CONFIG_TRAIN_SALT2,
    'TRAIN_SALT3'   : HELP_CONFIG_TRAIN_SALT3,
    'TRAIN_BAYESN'  : HELP_CONFIG_TRAIN_BAYESN,
    'MAKEDATAFILES' : HELP_CONFIG_MAKEDATAFILES,
    'MERGE'         : HELP_MERGE,
    'AIZ'           : HELP_AIZ     # ABORT_IF_ZERO
}


# === END: ===
