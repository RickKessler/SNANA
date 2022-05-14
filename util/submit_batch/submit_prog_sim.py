# ======================================
# Created July 2020 by R. Kessler
# 
# Batch control for class Simulation; replaces lagacy sim_SNmix.pl
#
# Upgrades over legacy sim_SNmix.pl :
#    + can run multiple SNIa jobs/models, analogous to SNCC
#
#    + SIMLOGS_XXX/MERGE.LOG is updated regularly with current state
#      (WAIT, RUN, DONE) and statistics. YAML format parsable by Pippin.
#      
#    + each job+args written to SIMLOGS_XXX/SCRIPTS/CPU*.CMD to 
#      leave explicit record of jobs
#
#    + FAIL-REPEAT scripts for failed jobs
#
#    + CPU diagnostics added to MERGE.LOG   
#
# MAYBE, TO-DO ...
#    - if SIMLOGS exists with no done file, give warning 
#
# ////////////////////////
#
#      HISTORY
#  ~~~~~~~~~~~~~~~
# Nov 16 2020: protect GENOPT args with () -> \(\)
# Dec 03 2020: replace legacy n_file_max=6 with global MAX_SIMGEN_INFILE=12
# Dec 09 2020: if NGENTOT_LC is in GENOPT, remove it after parsing it.
# Dec 10 2020: extract NGENTOT_LC frmo GENOPT_GLOBAL
# Jan 04 2021: fix setting ranseed_list when FORMAT_MASK=32
# Jan 06 2021: cidadd safety margin -> 1000 (was 10) to reduce chance
#               of running out of random CIDs
# Jan 14 2021: add MERGE.LOG column for NSPEC_WRITE
# Aug 18 2021: implement --snana_dir for SIMnorm jobs
# Sep 09 2021: abort of NGEN_UNIT is mixed with NGENTOT_LC in GENOPT
# Nov 09 2021: strip comment fields before reading sim-input files
#                (see is_comment_line)
# Dec 08 2021: write NLC_GEN[WRITE]_SUM to MERGE.LOG via get_misc_merge_info()
#
# Mar 07 2022: 
#  copy sim-input files of not merge_flag ... hopefully fixes rare problem
#  of sim jobs reading zero words from sim-input file and then aborting.
#
# Apr 7 2022
#   leave SUBMIT.INFO outside SIMLOGS.tar, and create BACKUP_*.tar.gz files
#   before SIMLOGS.tar (better organized) 
#
# May 3 2022: protect arg wildcards by adding quotes; see util.protect_wildcard
#
# May 14 2022: if RANSYSTPAR_RANSEED_SYST is set in genopt, sync the 
#               ranseed_syst args among isplit.
#
# ==========================================

import os,sys,glob,yaml,shutil
import logging, coloredlogs

import submit_util  as  util
from   submit_params    import *
from   submit_prog_base import Program

# define columns of MERGE.LOG 
COLNUM_SIM_MERGE_STATE        = 0     # current state; e.g., WAIT, RUN, DONE
COLNUM_SIM_MERGE_IVER         = 1     # GENVERSION index
COLNUM_SIM_MERGE_GENVERSION   = 2
COLNUM_SIM_MERGE_NLC_GEN      = 3
COLNUM_SIM_MERGE_NLC_WRITE    = 4
COLNUM_SIM_MERGE_NSPEC_WRITE  = 5     # added Jan 14 2021
COLNUM_SIM_MERGE_CPU          = 6     # CPU minutes
COLNUM_SIM_MERGE_NSPLIT       = 7     

# define keys for sim master-input ...
SIMGEN_MASTERFILE_KEYLIST_SNIa = [ 
    'SIMGEN_INFILE_Ia', 'SIMGEN_INFILE_SNIa', 
    'SIMGEN_INFILE_1a', 'SIMGEN_INFILE_SN1a' ]
SIMGEN_MASTERFILE_KEYLIST_NONIa = [ 
    'SIMGEN_INFILE_NONIa', 'SIMGEN_INFILE_NON1a' ]

SIMGEN_INPUT_LISTFILE = "INPUT_FILE.LIST" # contains list of input files

MAX_SIMGEN_INFILE = 12

RANSEED_KEYLIST = [ 'RANSEED_REPEAT', 'RANSEED_CHANGE' ]

KEY_INPUT_INCLUDE_FILE ='INPUT_INCLUDE_FILE' # Optional key

# Define keys to read from each underlying sim-input file
SIMGEN_INFILE_KEYCHECK = { # Narg  Required     sync-Verify
    "GENMODEL"          :     [ 1,    True,      False   ],
    "NGENTOT_LC"        :     [ 1,    False,     False   ],
    "TAKE_SPECTRUM"     :     [ 1,    False,     False   ],
    "FORMAT_MASK"       :     [ 1,    False,     True    ],
    "GENFILTERS"        :     [ 1,    True,      True    ],
    "PATH_USER_INPUT"   :     [ 1,    False,     True    ],
    "GENRANGE_REDSHIFT" :     [ 2,    True,      True    ],
    "GENRANGE_PEAKMJD"  :     [ 2,    True,      True    ],
    "SOLID_ANGLE"       :     [ 1,    True,      True    ]    }

# define GENOPT_GLOBAL key subStrings to ignore in the SIMnorm process;
# makes no difference in result, but in case of debug there is no need
# to sift thru so many unused arguments.
GENOPT_GLOBAL_IGNORE_SIMnorm = [
    "SIMGEN_DUMP", "HOSTLIB_", "SEARCHEFF_" ,  "GENMODEL_EXTRAP" ]

# format mask options for output data files
FORMAT_MASK_TEXT   = 2
FORMAT_MASK_FITS   = 32
FORMAT_MASK_CIDRAN = 16

FORMAT_TEXT = "TEXT"
FORMAT_FITS = "FITS"

# define max ranseed to avoid exceeding 4-byte limit of snlc_sim storage
RANSEED_MAX = 1000000000   # 1 billion

KEY_NGENTOT    = "NGENTOT_LC"

NGENTOT_CHECK_ABORT = 300

# - - - - - - - - - - - - - - - - - - -     -
class Simulation(Program):
    def __init__(self, config_yaml) :
        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_SIM
        super().__init__(config_yaml, config_prep)

    def set_output_dir_name(self):
        # check user-option LOGDIR; else default SIMLOGS_[GENPREFIX]
        # Since there is no explicit OUTDIR key in CONFIG, 
        # set CONFIG['OUTDIR'] here.

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file
        msgerr       = []
        if 'LOGDIR' in CONFIG :
            OUTDIR = CONFIG['LOGDIR']
            output_dir_name = os.path.expandvars(OUTDIR)
        else :
        # set default name based on GENPREFIX
            if 'GENPREFIX' not in CONFIG :
                msgerr.append(f"Must define GENPREFIX: <genprefix>")
                msgerr.append(f"Inside yaml-CONFIG")
                util.log_assert(False, msgerr)

            GENPREFIX = CONFIG['GENPREFIX']
            MAXLEN = 50 ; LEN = len(GENPREFIX)
            if LEN > MAXLEN :
                msgerr.append(f"GENPREFIX = {GENPREFIX}")
                msgerr.append(f"has stringLen = {LEN} ; ")
                msgerr.append(f"Len(GENPREFIX) must be below {MAXLEN} " 
                              f"due to FITS file limitation in header")
                msgerr.append(f"Check {input_file}")
                util.log_assert(False,msgerr)

            OUTDIR = f"SIMLOGS_{GENPREFIX}"
            output_dir_name = f"{CWD}/{OUTDIR}"

        self.config_yaml['CONFIG']['OUTDIR'] = OUTDIR
        return output_dir_name,SUBDIR_SCRIPTS_SIM
        # set_output_dir_name
 
    def submit_prepare_driver(self):
        # July 2020
        # prepare simulation arguments for each GENVERSION

        print("")

        # check where output data files are written
        self.sim_prep_PATH_SNDATA_SIM()

        # collect GENOPT_GLOBAL to append onto each version-specific GENOPT
        self.sim_prep_GENOPT_GLOBAL()

        # parse GENVERSION_LIST block
        self.sim_prep_GENVERSION_LIST()

        # extract list of SIMGEN-input files
        self.sim_prep_SIMGEN_INFILE()

        self.sim_prep_GENOPT() 

        # fetch format and do_cidran 
        self.sim_prep_FORMAT_MASK()

        # parse random seed options (REPEAT or CHANGE)
        self.sim_prep_RANSEED()

        # prepare 1D index lists with length n_job_tot
        self.sim_prep_index_lists()
        
        # determine now many to generate per job, and CIDRAN range          
        self.sim_prep_NGENTOT_LC()

        # determine CIDRAN parameters to ensure unique randoms
        self.sim_prep_CIDRAN()

        # abort on conflicts
        self.sim_check_conflicts()

        # end submit_prepare_driver (for sim)

    def sim_prep_index_lists(self):

        # Construct 1D sparse lists of
        #    version      index iver    0 to n_job_tot
        #    infile      index ifile    0 to n_job_tot
        #    split-job index isplit    0 to n_job_tot
        #
        # These 1D arrays can be used in 1D loops in range(0,n_job_tot)
        # instead of 3D for blocks.

        n_genversion  = self.config_prep['n_genversion']
        infile_list2d = self.config_prep['infile_list2d']
        n_job_split   = self.config_prep['n_job_split']
        n_job_tot     = 0
        iver_list=[];  ifile_list=[];  isplit_list=[]

        for iver in range(0,n_genversion) :
            n_file = len(infile_list2d[iver])
            for ifile in range(0,n_file):
                for isplit in range(0,n_job_split):
                    n_job_tot += 1
                    iver_list.append(iver)
                    ifile_list.append(ifile)
                    isplit_list.append(isplit)

        self.config_prep['n_job_tot']    = n_job_tot
        self.config_prep['n_done_tot']   = n_job_tot # same as n_job_tot
        self.config_prep['iver_list']    = iver_list
        self.config_prep['ifile_list']   = ifile_list
        self.config_prep['isplit_list']  = isplit_list

        # sim_prep_index_lists

    def sim_prep_GENOPT_GLOBAL(self):

        # there is a separate GENOPT_GLOBAL for SIMnorm jobs
        # to simplify debugging if needed.
        # M. Vincenzi Febr 2022: added protect_parenthesis in keys and values

        IGNORE_SIMnorm = GENOPT_GLOBAL_IGNORE_SIMnorm

        GENOPT_GLOBAL_STRING  = ""  # nominal
        GENOPT_GLOBAL_SIMnorm = ""  # for SIMnorm

        if 'GENOPT_GLOBAL' in self.config_yaml :
            GENOPT_GLOBAL  = self.config_yaml['GENOPT_GLOBAL']
            for key,value in GENOPT_GLOBAL.items():
                key_protect   = util.protect_parentheses(key)
                value_protect = util.protect_parentheses(value)
                value_protect = util.protect_wildcard(value_protect)

                GENOPT_GLOBAL_STRING += f"{key_protect} {value_protect}  "
                
                SKIP_SIMnorm = \
                    any(substring in key for substring in IGNORE_SIMnorm)
                if not SKIP_SIMnorm :
                    GENOPT_GLOBAL_SIMnorm += f"{key} {value}  "

        # check for ngentot_lc in GENOPT_GLOBAL (Dec 2020)
        genopt_replace, strval = \
                self.extract_value_from_genopt(GENOPT_GLOBAL_STRING,KEY_NGENTOT)
        if strval is not None : 
            ngentot_global       = int(strval)
            GENOPT_GLOBAL_STRING = genopt_replace
        else:
            ngentot_global = -9

        # load config_prep
        self.config_prep['genopt_global']          = GENOPT_GLOBAL_STRING
        self.config_prep['genopt_global_SIMnorm']  = GENOPT_GLOBAL_SIMnorm
        self.config_prep['ngentot_global']         = ngentot_global

        # end sim_prep_genopt_global

    def sim_prep_GENOPT(self):    

        # parge GENOPT and GENOPT(ARG) keys where ARG can be
        # Ia, SNIa, NONIa ... or ARG is a subtring to match against 
        # infile names
        # If '(ARG)' is not given, the associated GENOPT is applied
        # to all sim-input infiles for the specified GENVERSION.

        GENVERSION_LIST   = self.config_yaml['GENVERSION_LIST']
        merge_flag        = self.config_yaml['args'].merge_flag
        path_sndata_sim   = self.config_prep['path_sndata_sim']  
        n_genversion      = self.config_prep['n_genversion'] 
        infile_list2d     = self.config_prep['infile_list2d'] 
        model_list2d      = self.config_prep['model_list2d'] 

        genopt_list2d = []    # init array to load below
        verbose          = (not merge_flag)

        for iver in range(0,n_genversion) :
            GENV         = GENVERSION_LIST[iver]
            GENVERSION     = GENV['GENVERSION']
            n_file         = len(infile_list2d[iver])

            if verbose :
                print(f" sim_prep_GENOPT({GENVERSION}) ")

            # catenate the GENOPTs into one long GENOPT_FINAL string
            genopt_list = []  # list for this GENVERSION
            genarg_list = []
            genopt_list2d.append([''] * n_file)     # init ifile index

            for key,value in GENV.items():
                if 'GENOPT' in key :
                    genarg = util.extract_arg(key) # optional arg from ()
                    #print(f"\t xxx {key} -> arg = {genarg}")
                    for key2,value2 in value.items():
                        # if key includes (), replace with \( \)
                        key2_protect = util.protect_parentheses(key2)
                        value2_protect = util.protect_parentheses(value2)
                        value2_protect = util.protect_wildcard(value2_protect)
                        genopt = f"{key2_protect} {value2_protect}     "
                        #print(f" xxx genopt = {genopt}")
                        genopt_list.append(genopt)
                        genarg_list.append(genarg)

            # set file-specific GENOPT 
            n_arg = len(genopt_list)

            for ifile in range(0,n_file):
                infile = infile_list2d[iver][ifile]
                model  = model_list2d[iver][ifile]
                for i_arg in range(0,n_arg):
                    genarg = genarg_list[i_arg]
                    genopt = genopt_list[i_arg]
                    match  = self.genopt_arg_match(infile,model,genarg)
                    if match :
                        genopt_list2d[iver][ifile] += genopt

                if verbose:
                    print(f"\t GENOPT({infile}): {genopt_list2d[iver][ifile]}")

            iver += 1

        # store genopt 
        self.config_prep['genopt_list2d']  = genopt_list2d

        #sys.exit("\n xxx DEBUG DIE xxx \n")

        # end sim_prep_GENOPT

    def genopt_arg_match(self, infile, model, genarg):
        # Returns True if GENOPT is applied to param inputs :
        #  infile   = name of sim-input file
        #  model    = SNIa or NONIa
        #  genarg   = argument in GENOPT(genarg)
        #

        match_SNIa        = False
        match_NONIa       = False
        match_infile      = False
        genarg_SNIa_list  = [ 'Ia', '1a', 'SNIa', 'SN1a' ] 
        genarg_NONIa_list = [ 'NONIa', 'NON1a' ]

        # return True immediatly if genopt_arg is blank
        if len(genarg) == 0 :
            return True

        # check for explicit SNIa tag
        if model == MODEL_SNIa :
            for string in genarg_SNIa_list :
                if genarg == string :
                    match_SNIa = True

        # check for explicit NONIa tag
        if model == MODEL_NONIa :
            for string in genarg_NONIa_list :
                if genarg == string :
                    match_NONIa = True

        # check for substring of input file name
        if genarg in infile :
            match_infile = True

        # - - - - - - - - -
        found_match = (match_SNIa or match_NONIa or match_infile)

        dump_flag = False
        if dump_flag :
            print(f" xxx -------------------------------------- ")
            print(f" xxx infile = {infile} ")
            print(f" xxx genarg = {genarg} ")
            print(f" xxx model  = {model}  ")
            print(f" xxx match[SNIa,NON1a,infile] = " \
                  f"{match_SNIa},{match_NONIa},{match_infile}")
            print(f" xxx found_match = {found_match}")

        # - - - - -
        if found_match :
            return True
        else:
            return False
            #msgerr = [] 
            #msgerr.append(f"Invalid GENOPT({genarg}):")
            #msgerr.append(f"Valid GENOPT args in () are:")
            #msgerr.append(f"   Any element of {genarg_SNIa_list}")
            #msgerr.append(f"   Any element of {genarg_NONIa_list}")
            #msgerr.append(f"   substring of '{infile}'")
            #self.log_assert(False,msgerr)            

        # end genopt_arg_match

    def sim_prep_GENVERSION_LIST(self):
        # read and store each genversion, and remove pre-existing
        # output directory

        CONFIG              = self.config_yaml['CONFIG']
        GENVERSION_LIST     = self.config_yaml['GENVERSION_LIST']
        nosubmit            = self.config_yaml['args'].nosubmit
        merge_flag          = self.config_yaml['args'].merge_flag
        path_sndata_sim     = self.config_prep['path_sndata_sim']  
        genversion_list     = []
        IS_RANSEED_CHANGE   = 'RANSEED_CHANGE' in CONFIG

        for GENV in GENVERSION_LIST :
            GENVERSION     = GENV['GENVERSION']
            genversion_list.append(GENVERSION)

            # if GENVERSION exists on disk, remove it ... 
            # unless nosubmit or merge_flag flag is set
            if nosubmit is False and merge_flag is False :
                if IS_RANSEED_CHANGE :
                    wildcard_genversion = f"{GENVERSION}-*"
                else:
                    wildcard_genversion = f"{GENVERSION}"

                v_list = glob.glob1(path_sndata_sim,wildcard_genversion)
                for v in v_list:
                    shutil.rmtree(f"{path_sndata_sim}/{v}")

        self.config_prep['genversion_list']  = genversion_list
        self.config_prep['n_genversion']     = len(genversion_list)

    # end sim_prep_GENVERSION_LIST

    def sim_prep_RANSEED(self):
        # Parse randon seed(s) and determine number of split jobs.
        # Note that RANSEED_REPEAT splits a sim into sub-jobs, 
        # then re-combines into one effective sim job.
        # RANSEED_CHANGE splits into sub-jobs; does NOT re-combine.
        #
        # Oct 12 2021: if --check_abort, n_job_split -> 1

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file     # for msgerr
        check_abort  = self.config_yaml['args'].check_abort
        do_cidran    = self.config_prep['do_cidran']
        nkey_found   = 0  # local: number of valid RANSEED_XXX keys
        ranseed_list = [] # array vs. split index
        msgerr       = []

        if 'RANSEED' in CONFIG:
            msgerr.append(f"Invalid RANSEED: key")
            msgerr.append(f"Use RANSEED_REPEAT: or RANSEED_CHANGE:")
            self.log_assert(False,msgerr)

        RANSEED_KEY = ""
        for key in RANSEED_KEYLIST:
            if key in CONFIG:
                nkey_found    += 1                
                RANSEED_KEY     = key
                RANSEED_LIST    = CONFIG[key]
                n_job_split     = int(RANSEED_LIST.split()[0])
                ranseed         = int(RANSEED_LIST.split()[1])
                if check_abort: n_job_split=1

                if ranseed > RANSEED_MAX :
                    msgerr.append(f"ranseed = {ranseed} is too big.")
                    msgerr.append(f"ranseed must be under {RANSEED_MAX} " \
                                  f"to store as 4 byte integer.")
                    msgerr.append(f"Check {key} in {input_file}")
                    self.log_assert(False,msgerr)                    

                if key == 'RANSEED_REPEAT' and do_cidran :
                    # same seed per core since snlc_sim.exe will change
                    # each seed internally as it finds new CID list
                    ranseed_list = [ranseed] * n_job_split
                else:
                    # change seed each core 
                    self.prep_input_file_include(n_job_split) # Check for optional input systematics
                    for job in range(0,n_job_split):
                        ranseed_tmp = ranseed + 10000*job + job*job + 13
                        ranseed_list.append(ranseed_tmp)

        # - - - - - - -
        # require 1 and only 1 RANSEED_XXX key; otherwise abort.
        if nkey_found != 1 :
            msgerr.append(f"Found {nkey_found} RANSEED keys -> INVALID.")
            msgerr.append(f"Must specify either: "
                          "RANSEED_REPEAT or RANSEED_CHANGE")
            msgerr.append(f"in YAML CONFIG.")
            msgerr.append(f" (and note that RANSEED: is not valid)")
            self.log_assert(False,msgerr)

        # - - - - - - - 
        # expand list of versions based on ranseed key
        genv_list = self.config_prep['genversion_list']
        genversion_list_all,igenver_list_all = \
            self.genversion_expand_list(genv_list,RANSEED_KEY,n_job_split) 

        self.config_prep['n_job_split']         = n_job_split 
        self.config_prep['ranseed_list']        = ranseed_list
        self.config_prep['ranseed_key']         = RANSEED_KEY
        self.config_prep['genversion_list_all'] = genversion_list_all
        self.config_prep['igenver_list_all']    = igenver_list_all

        # end sim_prep_RANSEED

    def prep_input_file_include(self, n_job_split):
       
        CONFIG        = self.config_yaml['CONFIG']
        key           = KEY_INPUT_INCLUDE_FILE
        msgerr        = []
        self.config_prep['include_list'] = []
        self.config_prep['n_include_files'] = 0
        
        if key not in CONFIG:
            return 
        include_list_file = os.path.expandvars(CONFIG[key]) # File containing list of include files

        if not os.path.isfile(include_list_file):
            msgerr.append(f"include_list_file = {include_list_file}")
            msgerr.append(f"{include_list_file} does not exist!")
            self.log_assert(False, msgerr)
        
        with open(include_list_file, "rt") as f:
            open_file = f.read()
        if open_file == "":
            msgerr.append(f"Contents of {include_list_file} = '{open_file}'")
            msgerr.append(f"is empty!")
            self.log_assert(False, msgerr)

        include_list = [f for f in open_file.split("\n") if f != ""] # Assumed each line is ONLY a file name.
        n_include_list = len(include_list)
        print(f"Found {n_include_list} include files from {include_list_file}")

        if n_include_list < n_job_split:
            msgerr.append(f"n_job_split = {n_job_split}, n_include_list = {n_include_list}")
            msgerr.append(f"n_include_list must be >= n_job_split")
            self.log_assert(False, msgerr)

        # Get the full path to each file
        parent_path = os.path.dirname(include_list_file)
        full_path_include_list = []
        for filename in include_list:
            if '/' in filename:
                full_filename = filename
            else:
                full_filename = os.path.join(parent_path, filename)
            if not os.path.isfile(full_filename):
                msgerr.append(f"full_filename = {full_filename}, derived from '{filename}'")
                msgerr.append(f"does not exist!")
                msgerr.append(f"Check contents of {include_list_file}")
                self.log_assert(False, msgerr)
            full_path_include_list.append(full_filename)

        self.config_prep['include_list'] = full_path_include_list
        self.config_prep['n_include_files'] = n_include_list

        # end prep_input_file_include

    def sim_prep_NGENTOT_LC(self):

        CONFIG        = self.config_yaml['CONFIG']
        prescale      = self.config_yaml['args'].prescale
        n_genversion  = self.config_prep['n_genversion']
        infile_list2d = self.config_prep['infile_list2d']
        INFILE_KEYS   = self.config_prep['INFILE_KEYS']
        n_core        = self.config_prep['n_core']
        n_job_split   = self.config_prep['n_job_split']
        do_cidran     = self.config_prep['do_cidran']
        key_ngen_unit = "NGEN_UNIT"
        ngentot_sum   = 0 

        # check for NGEN_UNIT
        if key_ngen_unit in CONFIG:
            ngen_unit = float(CONFIG[key_ngen_unit])
        else:
            ngen_unit = -1.0
            
        # print(f" xxx ngen_unit = {ngen_unit}" )
        ngentot_list2d    = []
        for iver in range(0,n_genversion):
            n_file = len(infile_list2d[iver])
            ngentot_list = []
            for ifile in range(0,n_file):
                if ngen_unit < 0.0 :
                    ngentot = self.get_ngentot_from_input(iver,ifile,+1)
                else:
                    dummy   = self.get_ngentot_from_input(iver,ifile,-1)
                    ngentmp = self.get_ngentot_from_rate(iver,ifile) 
                    ngentot = int(ngen_unit * ngentmp)

                # finally, check for fast option to divide by 10 or 100
                if prescale > 1 :  
                    ngentot = int(ngentot/prescale)

                ngentot_list.append(ngentot) # append ifile dimension
                ngentot_sum += ngentot

            ngentot_list2d.append(ngentot_list)     # fill iver dimension

        self.config_prep['ngentot_list2d'] = ngentot_list2d
        self.config_prep['ngen_unit']      = ngen_unit
        self.config_prep['ngentot_sum']    = ngentot_sum

        # Jan 4 2021: move this test here from sim_prep_FORMAT_MASK
        if do_cidran and ngentot_sum == 0 :
            msgerr = []
            msgerr.append(f"Invalid NGENTOT_LC=0  with CIDRAN option " \
                          f" (FORMAT_MASK += {FORMAT_MASK_CIDRAN}) ")
            msgerr.append(f"Try one of the following:")
            msgerr.append(f"  - FORMAT_MASK -= {FORMAT_MASK_CIDRAN}:")
            msgerr.append(f"  - specify NGENTOT_LC")
            msgerr.append(f"  - specify NGEN_UNIT")
            self.log_assert(False,msgerr)

        # end sim_prep_NGENTOT_LC

    def get_ngentot_from_input(self,iver,ifile,opt):

        # Sep 9 2021: add opt input:
        # opt > 0 -> return ngentot from input file or GENOPT key
        # opt < 0 -> abort of NGENTOT_LC is in GENOPT key 
        #             (to avoid conflict with NGEN_UNIT)

        INFILE_KEYS      = self.config_prep['INFILE_KEYS']
        ngentot_global   = self.config_prep['ngentot_global']
        genopt_list2d    = self.config_prep['genopt_list2d']
        genopt           = genopt_list2d[iver][ifile]
        found_GENOPT_override = False 
        
        ngentot = 0

        # Start with default NGENTOT_LC is from the sim-input file
        if KEY_NGENTOT in INFILE_KEYS[iver][ifile]:
            ngentot  = INFILE_KEYS[iver][ifile][KEY_NGENTOT]

        # check for GENOPT override
        genopt_replace, strval = \
                self.extract_value_from_genopt(genopt,KEY_NGENTOT)
        if strval is not None : 
            ngentot = int(strval)
            genopt_list2d[iver][ifile]        = genopt_replace
            self.config_prep['genopt_list2d'] = genopt_list2d
            found_GENOPT_override = True

        # check for GENOPT_GLOBAL override (Dec 10 2020)
        if ngentot_global > 0 :
            ngentot = ngentot_global
            found_GENOPT_override = True

        # 9.09.2021: check option to abort on NGENTOT key
        if found_GENOPT_override and opt < 0 :
            msgerr = []
            input_file   = self.config_yaml['args'].input_file 
            msgerr.append(f"Cannot mix NGEN_UNIT and NGENTOT_LC in")
            msgerr.append(f"{input_file} .")
            msgerr.append(f"Either remove NGENTOT_LC from GENOPT[_GLOBAL]")
            msgerr.append(f"or remove NGEN_UNIT from CONFIG.")
            self.log_assert(False,msgerr)

        return ngentot

    def extract_value_from_genopt(self,genopt,key):
        # Created Dec 10 2020
        # Extract value after key from genopt;
        # return genopt with "key value=" removed, and return value.
        # Example:
        #   input genopt = "KEY1 46  KEY2 33  KEY3 1.66" and key= KEY2
        #   Function returns
        #     genopt_replace = "KEY1 46  KEY3 1.66
        #     value          = 33
        # If key not in genopt, return genopt_replace=genopt and value=None

        genopt_replace = genopt  # default return value
        value          = None    # default return value
        genopt_words   = genopt.split()
        if key in genopt_words :
            jindx       = genopt_words.index(key)
            value       = genopt_words[jindx+1]
            tmp = genopt_words[:jindx] + genopt_words[jindx+2:]
            genopt_replace = "  ".join(tmp)

        return genopt_replace, value
        # end extract_value_from_genopt

    def get_ngentot_from_rate(self,iver,ifile):

        # run sim with INIT_ONLY flag to use sim as rate-calculator
        # that quits immediately without generating events. Basic
        # idea is that one "UNIT" of generated events is computed from
        # GENRANGE_PEAKMJD, GENRANGE_REDSHIFT and SOLID_ANGLE.
        # User-input NGEN_UNIT (master-input) multiples NGENTOT_LC
        # for one UNIT to get final NGENTOT_LC. Example: if user ranges
        # result in NGENTOT_LC = 4000 for one UNIT, and user input
        # has "NGEN_UNIT: 5", then NGENTOT_LC -> 4000x5 = 20000.
        # (see get_normalization in sim_SNmix.pl)
        # Seens that sim_SNmix did not inlcude GENOPTs for SIMnorm
        # step, so we'll include those here. Should rarely, if ever,
        # make a difference, unless GENOPT change REDSHIFT, PEAKMJD,
        # or SOLID_ANGLE.
        
        msgerr        = []
        key_ngentot   = "NGENTOT_RATECALC:"
        genversion    = self.config_prep['genversion_list'][iver]
        model         = self.config_prep['model_list2d'][iver][ifile] # SNIa or NONIa
        infile        = self.config_prep['infile_list2d'][iver][ifile]
        program       = self.config_prep['program']
        output_dir    = self.config_prep['output_dir']
        genopt_global = self.config_prep['genopt_global_SIMnorm']
        genopt        = self.config_prep['genopt_list2d'][iver][ifile]
        snana_dir     = self.config_yaml['args'].snana_dir

        cddir         = f"cd {output_dir}"
        ngentot       = 0
        
        model_string  = self.model_string_suffix(model,ifile)
        prefix        = f"SIMnorm_{genversion}_{model_string}"
        log_file      = f"{prefix}.LOG"
        LOG_FILE      = f"{output_dir}/{log_file}"

        arg_list  = ""
        arg_list += f"INIT_ONLY 1 "
        arg_list += f"{genopt} "
        arg_list += f"{genopt_global} "

        # 1.31.2021: check for PATH_USER_INPU
        INFILE_KEYS  = self.config_prep['INFILE_KEYS']
        key = 'PATH_USER_INPUT'
        if key not in INFILE_KEYS[iver][ifile] :
            arg_list += f"{key} {CWD} "

        # contruct two sets of strings.
        # cmd_string is passed to os.system
        # msgerr is error message in case of problem.

        cmd_array = []     # for screen dump
        cmd_array.append(f"{cddir} ; ")

        if snana_dir is None :
            cmd_array.append(f"{program} {infile} ")
        else:
            cmd_array.append(f"{snana_dir}/bin/{program} {infile} ")

        cmd_array.append(f"{arg_list} ")
        cmd_array.append(f"> {log_file}")

        cmd_stdout = []
        cmd_string = "" # for os.system
        for cmd in cmd_array:
            cmd_string += cmd
            cmd_stdout.append(f"  {cmd} \\")  # allows cut-and-paste

        #print(f"{cmd_string}")
        os.system(cmd_string)

        # check if LOG_FILE is really there
        msgerr = []
        msgerr.append(f"Sim normalization commands :")
        msgerr += cmd_stdout
        util.check_file_exists(LOG_FILE,msgerr)

        # read value after NGENTOT_RATECALC key
        found_key = False
        with open(LOG_FILE, 'r') as f :
            for line in f:
                if len(line.strip()) > 1 :
                    words = line.split()
                    if words[0] == key_ngentot :
                        ngentot = int(words[1])
                        found_key = True

        # - - - - - -
        msgerr = []
        if not found_key:
            msgerr.append(f"Unable to find {key_ngentot} key in {log_file} ;")
            msgerr.append(f"LOG created from sim normalization commands :")
            msgerr += cmd_stdout
            self.log_assert(False,msgerr)

        if ngentot == 0 :
            msgerr.append(f"ngentot=0 in {log_file} ;")
            msgerr.append(f"LOG created from sim normalization commands :")
            msgerr += cmd_stdout
            self.log_assert(False,msgerr)
            
        print(f"  Compute NGENTOT={ngentot:6d} for {prefix}")

        return ngentot
        # end get_ngentot_from_rate

    def sim_prep_FORMAT_MASK(self):

        # Determine format_mask from input.
        # format_mask can be either in GENOPT_GLOBAL, or as CONFIG key;
        # require 1 and only 1; abort on 0 or 2 keys
        # Also parse format_mask to determine if TEXT or FITS,
        # and option for random CIDs 

        CONFIG          = self.config_yaml['CONFIG']
        input_file      = self.config_yaml['args'].input_file     # for msgerr
        genopt_global   = self.config_prep['genopt_global'].split()
        key             = 'FORMAT_MASK'
        format_mask     = -1 # init value
        msgerr = []

        nkey = 0  
        if key in CONFIG :
            format_mask = CONFIG[key]
            nkey += 1

        if key in genopt_global :
            jindx       = genopt_global.index(key)
            format_mask = int(genopt_global[jindx+1])
            nkey += 1

        if nkey != 1 :
            msgerr.append(f"Found {nkey} {key} keys;")
            msgerr.append(f"1 and only 1 {key} key must be specified.")
            msgerr.append(f"Check CONFIG and GENOPT_GLOBAL in {input_file}")
            self.log_assert( False , msgerr)

        if format_mask & FORMAT_MASK_TEXT > 0 :
            format_type = FORMAT_TEXT
            msgerr.append(f"TEXT format not yet supported;")
            msgerr.append(f"Only FITS format is supported " \
                          f"(FORMAT_MASK += {FORMAT_MASK_FITS})")
            self.log_assert( False , msgerr)

        elif format_mask & FORMAT_MASK_FITS > 0 :
            format_type = FORMAT_FITS
        else:
            strtmp = f"{FORMAT_MASK_TEXT}(TEXT) or {FORMAT_MASK_FITS}(FITS)"
            msgerr.append(f"Invalid FORMAT_MASK = {format_mask}")
            msgerr.append(f"FORMAT_MASK must include either {strtmp}")
            self.log_assert( False , msgerr)

        print(f"  FORMAT_MASK = {format_mask} ({format}) ")

        # check option for random CIDs
        do_cidran  = (format_mask & FORMAT_MASK_CIDRAN) > 0

        self.config_prep['format_mask'] = format_mask
        self.config_prep['format']      = format_type        # TEXT or FITS
        self.config_prep['do_cidran']   = do_cidran

        # end sim_prep_FORMAT_MASK

    def sim_prep_CIDRAN(self):

        # Compute CIDOFF for each version/model/splitJob
        # based on logic of RANSEED_REPEAT[CHANGE] and RESET_CIDOFF.
        # This logic is quite tricky, so be careful.
        #
        # A key user-input flag is      RESET_CIDOFF: <value>
        #    0 -> use CIDOFF from input files
        #          beware that duplicate CIDs may exist 
        #
        #    1 -> start CIDOFF=0, and +gentot for each split job and version.
        #         Ensures unique CID among all events with each GENVERSION,
        #         but repeat CID may exist in different GENVERSIONs
        #
        #    2 -> start CIDOFF=0, and +gentot (never resets)
        #         Ensures unique CID among all models and GENVERSIONs
        #
        # Note that (FORMAT_MASK & 16)>0 (random CIDs) will automatically
        # set RESET_CIDOFF=1 if not set by user.
        #
        # ------------

        CONFIG         = self.config_yaml['CONFIG']
        INFILE_KEYS    = self.config_prep['INFILE_KEYS']
        genopt_global  = self.config_prep['genopt_global'].split()
        format_mask    = self.config_prep['format_mask']
        do_cidran      = self.config_prep['do_cidran']
        ngentot_list2d = self.config_prep['ngentot_list2d'] # per split job
        infile_list2d  = self.config_prep['infile_list2d']

        iver_list      = self.config_prep['iver_list']
        ifile_list     = self.config_prep['ifile_list']
        isplit_list    = self.config_prep['isplit_list']

        n_genversion   = self.config_prep['n_genversion']
        n_job_tot      = self.config_prep['n_job_tot']
        n_job_split    = self.config_prep['n_job_split']
        n_file_max     = MAX_SIMGEN_INFILE
        cidoff_list3d  = [[[0 for k in range(0,n_job_split)] for j in range(0,n_file_max)] for i in range(0,n_genversion)]
        cidran_max_list = [0] * n_genversion
        cidran_min      =    0
        cidran_max      =    0
        reset_cidoff    =    0

        # check for optional min CIDOFF = CIDRAN_MIN
        key = 'CIDRAN_MIN'
        if key in CONFIG :
            cidran_min = CONFIG[key]

        # - - - - - 
        key = 'RESET_CIDOFF'
        if key in CONFIG :
            reset_cidoff = CONFIG[key]

        # check for auto-settings of reset_cidoff based on other user options

        # check if user forgot to set RESET_CIDOFF with random CIDs
        if do_cidran and reset_cidoff == 0 :
            reset_cidoff = 1  

        # no need for 2-option if split jobs are NOT combined
        if reset_cidoff == 2 and 'RANSEED_CHANGE' in CONFIG :
             reset_cidoff = 1

        # - - - - - 
        # xxx ngentot_sum = 0
        cidoff        = cidran_min
        for iver,ifile,isplit in zip(iver_list,ifile_list,isplit_list):

            new_version = (ifile == 0 and isplit==0)
            if reset_cidoff < 2 and new_version :
                cidoff = 0      # reset CIDOFF for new version
            if reset_cidoff < 2 and isplit == 0 :  # next file in version 
                cidoff = util.roundup_first_digit(cidoff) # 
            # note: for reset_cidoff=2, cidoff never resets

            cidoff_list3d[iver][ifile][isplit] = cidoff

            ngentot      = ngentot_list2d[iver][ifile] # per split job
            # xxx ngentot_sum += ngentot    # increment total number generated
            if reset_cidoff > 0 :
                cidadd       = int(ngentot*1.1)+1000   # leave safety margin
                cidoff      += cidadd        # for random CIDs in snlc_sim
                cidran_max   = cidoff
            else:
                cidoff += ngentot
                
            if isplit == n_job_split-1 :
                cidran_max_list[iver] = cidran_max

        # for unique CIDs everywhere, cidran_max must be the
        # same for all jobs
        if reset_cidoff == 2 :
            cidran_max_list = [cidran_max] * n_genversion

        # store info
        self.config_prep['reset_cidoff']      = reset_cidoff
        self.config_prep['cidran_min']        = cidran_min
        self.config_prep['cidran_max_list']   = cidran_max_list
        self.config_prep['cidoff_list3d']     = cidoff_list3d
        self.config_prep['n_job_tot']         = n_job_tot

        self.sim_prep_dump_cidoff()

        # end sim_prep_CIDRAN

    def sim_check_conflicts(self):

        # misc conflict checks/aborts.
        
        INFILE_KEYS    = self.config_prep['INFILE_KEYS']
        ngentot_sum    = self.config_prep['ngentot_sum']
        n_genversion   = self.config_prep['n_genversion']
        infile_list2d  = self.config_prep['infile_list2d']
        TAKE_SPECTRUM  = False 
        msgerr         = []

        # - - - - - - -
        # avoid generating spectra for very large jobs such as biasCor.
        for iver in range(0,n_genversion):
            n_infile   = len(infile_list2d[iver])
            for ifile in range(0,n_infile):
                if 'TAKE_SPECTRUM' in INFILE_KEYS[iver][ifile]:
                    TAKE_SPECTRUM = True

        MXGENTOT_TAKE_SPECTRUM = 200000
        if TAKE_SPECTRUM and ngentot_sum > MXGENTOT_TAKE_SPECTRUM :
            msgerr.append(f"ngentot_sum = {ngentot_sum} is too large " \
                          f"with TAKE_SPECTRUM keys.")
            msgerr.append(f"MXGENTOT_TAKE_SPECTRUM = {MXGENTOT_TAKE_SPECTRUM}")
            msgerr.append(f"Are spectra really needed for such " \
                          f"a large sample?")
            self.log_assert( False , msgerr)

        # - - - - - 

        # end sim_check_conflicts

    def sim_prep_dump_cidoff(self):
        # debug dump of cidoff for each genversion and model/file
        n_genversion      = self.config_prep['n_genversion']
        genversion_list   = self.config_prep['genversion_list']
        infile_list2d     = self.config_prep['infile_list2d']
        model_list2d      = self.config_prep['model_list2d']
        cidoff_list3d     = self.config_prep['cidoff_list3d']
        reset_cidoff      = self.config_prep['reset_cidoff']
        cidran_min        = self.config_prep['cidran_min']
        cidran_max_list   = self.config_prep['cidran_max_list']
        n_job_tot         = self.config_prep['n_job_tot']
        n_job_split       = self.config_prep['n_job_split']
        ngentot_sum       = self.config_prep['ngentot_sum']
        
        msgerr = []

        print(f"")
        #print(f" DUMP CIDOFF vs. GENVERSION and MODEL/INFILE")
        print(f"     GENVERSION        MODEL        CIDOFF(GENVERSION)")
        print(f"# -------------------------------------------------- ")

        for iver in range(0,n_genversion) :
            genv = genversion_list[iver]
            n_file = len(infile_list2d[iver])
            if n_file > MAX_SIMGEN_INFILE :
                msgerr.append(f"n_file[iver={iver}] = {n_file}")
                msgerr.append(f"exceeds bound of MAX_SIMGEN_INFILE " \
                              f"= {MAX_SIMGEN_INFILE}")
                self.log_assert(False,msgerr)

            for ifile in range(0,n_file):
                str_genv   = f"{genv:20.20}"
                str_model  = str(model_list2d[iver][ifile]) + '-' + str(ifile)
                str_model  = f"{str_model:8.8}"
                cidoff_list = cidoff_list3d[iver][ifile]
                len_list    = len(cidoff_list)
                if len_list < 6 :
                    str_cidoff = f"{cidoff_list}"
                else:
                    str_cidoff = f"{cidoff_list[0:2]} ... {cidoff_list[len_list-2:len_list]}"

                print(f" {str_genv} {str_model} {str_cidoff}")

        print(f"# -------------------------------------------------- ")
        print(f"  RESET_CIDOFF        = {reset_cidoff} ")
        print(f"  CIDRAN_MIN        = {cidran_min}")
        print(f"  CIDRAN_MAX(genver)= {cidran_max_list}")
        print(f"  Sum of NGENTOT_LC = {ngentot_sum} x {n_job_split} " \
              f"(distribute among {n_job_tot} jobs)")
        print(f"")

        # end sim_prep_dump_cidoff

    def sim_prep_SIMGEN_INFILE(self):

        # get default list of sim-input infiles
        infile_list2d_SNIa_default, infile_list2d_NONIa_default = \
                self.sim_prep_SIMGEN_INFILE_defaults()

        # check for genv-dependent overrides
        infile_list2d_SNIa_genv, infile_list2d_NONIa_genv = \
                self.sim_prep_SIMGEN_INFILE_genversion()

        # Prepare final list of sim-input files using defaults
        # and genversion-overrides.
        # Logic here is tricky: allow overriding either Ia or NONIa or
        # both Ia+NONIa sim-input files.

        CONFIG                = self.config_yaml['CONFIG']
        n_genversion          = self.config_prep['n_genversion']
        infile_list2d_SNIa    = infile_list2d_SNIa_default
        infile_list2d_NONIa   = infile_list2d_NONIa_default
        infile_list2d         = [] * n_genversion
        model_list2d          = [] * n_genversion
        for iver in range(0,n_genversion):              

            tmp_list = infile_list2d_SNIa_genv[iver]
            if len(tmp_list) > 0 :
                infile_list2d_SNIa[iver] = tmp_list

            tmp_list = infile_list2d_NONIa_genv[iver]
            if len(tmp_list) > 0 :
                infile_list2d_NONIa[iver] = tmp_list

            n_SNIa  = len(infile_list2d_SNIa[iver])
            n_NONIa = len(infile_list2d_NONIa[iver])

            infile_list = infile_list2d_SNIa[iver] + infile_list2d_NONIa[iver]
            model_list  = [MODEL_SNIa]*n_SNIa + [MODEL_NONIa] * n_NONIa
            infile_list2d.append(infile_list)
            model_list2d.append(model_list)

        # store final 2D infile lists
        self.config_prep['infile_list2d']        =  infile_list2d  
        self.config_prep['infile_list2d_SNIa']   =  infile_list2d_SNIa 
        self.config_prep['infile_list2d_NONIa']  =  infile_list2d_NONIa
        self.config_prep['model_list2d']         =  model_list2d 
                
        # read a few keys from each INFILE, including INCLUDE file(s).
        # ?? for duplicate infile, copy key_dict instead of re-reading ??
        INFILE_KEYS              = []
        include_file_list_unique = [] # needed later to copy files
        for iver in range(0,n_genversion):
            keyval_dict_list = []
            n_file = len(infile_list2d[iver])
            for ifile in range(0,n_file):
                infile   = infile_list2d[iver][ifile]
                key_dict = {} ; include_files = []
                key_dict,include_file_list = \
                    self.sim_prep_SIMGEN_INFILE_read(infile)
                keyval_dict_list.append(key_dict)

                # store unique list of include files, and don't bother
                # keep track of which genversion or infile they are in
                for infile in include_file_list:
                    if infile not in include_file_list_unique :
                        include_file_list_unique.append(infile)

            INFILE_KEYS.append(keyval_dict_list)

        # update config_prep with 2D arrays: [iver][ifile]
        self.config_prep['INFILE_KEYS']              = INFILE_KEYS 
        self.config_prep['include_file_list_unique'] = include_file_list_unique

        # After loading config_prep, verify some sim-input keys in each
        # GENVERSION.
        # INFILE_PAIRS_VERIFY_ERROR is used to avoid duplicate printing
        # of errors for repeated infile pairs.
        msgerr = [] ; nerr=0
        self.config_prep['INFILE_PAIRS_VERIFY_ERROR'] = [] 
        for iver in range(0,n_genversion):
            for keycheck in SIMGEN_INFILE_KEYCHECK :
                nerr+=self.sim_prep_SIMGEN_INFILE_verify(iver,keycheck,msgerr)
                        
        msgerr.append(f"{nerr} verify errors")
        msgerr.append(f"Check sim-input files and INCLUDE files")
        self.log_assert(nerr==0,msgerr)

        # copy input files to outdir/simlogs directory
        merge_flag        = self.config_yaml['args'].merge_flag
        if not merge_flag:
            self.sim_prep_SIMGEN_INFILE_copy()

        # end sim_prep_SIMGEN_INFILE


    def sim_prep_SIMGEN_INFILE_defaults(self):

        # read default sim-input files from CONFIG yaml block
        # Functions returns two 2D lists
        #  1. list of SNIa input files per version
        #  2. list of NONIa input files per version
        #
        # The defaults here have the same infile lists for each version.

        CONFIG              = self.config_yaml['CONFIG']
        # read default infile_list

        KEYLIST = SIMGEN_MASTERFILE_KEYLIST_SNIa
        infile_list_SNIa = util.get_YAML_key_values(CONFIG,KEYLIST)
                                 
        KEYLIST = SIMGEN_MASTERFILE_KEYLIST_NONIa
        infile_list_NONIa = util.get_YAML_key_values(CONFIG,KEYLIST)

        # fix partial paths by adding CWD if needed.
        infile_list_SNIa  = util.fix_partial_path(infile_list_SNIa)
        infile_list_NONIa = util.fix_partial_path(infile_list_NONIa)

        n_infile_SNIa  = len(infile_list_SNIa)
        n_infile_NONIa = len(infile_list_NONIa)
        n_infile       = n_infile_SNIa + n_infile_NONIa

        if n_infile == 0 :
            msgerr.append(f"Found no SIMGEN-input file keys.")
            msgerr.append(f"Check YAML keys SIMGEN_INFILE_Ia & "
                          "SIMGEN_INFILE_NONIa.")
            self.log_assert(False,msgerr)

        # load the same default infile lists for each genversion
        infile_list2d_SNIa_default  = []
        infile_list2d_NONIa_default = []
        n_genversion  = self.config_prep['n_genversion']
        for iver in range(0,n_genversion):
            infile_list2d_SNIa_default.append(infile_list_SNIa)
            infile_list2d_NONIa_default.append(infile_list_NONIa)

        return infile_list2d_SNIa_default, infile_list2d_NONIa_default

        #end sim_prep_SIMGEN_INFILE_defaults

    
    def sim_prep_SIMGEN_INFILE_genversion(self):

        # load 2D list of infile[version][ifile] for
        # overrides under GENOPT keys

        GENVERSION_LIST       = self.config_yaml['GENVERSION_LIST']
        infile_list2d_SNIa    = [] # init output 2D array of infiles
        infile_list2d_NONIa   = [] # idem for NONIa

        for GENV in GENVERSION_LIST :

            # fetch optional sim-input files here to override CONFIG defaults
            infile_list_SNIa  = []
            infile_list_NONIa = []

            KEYLIST = SIMGEN_MASTERFILE_KEYLIST_SNIa
            infile_list_SNIa = util.get_YAML_key_values(GENV,KEYLIST)
                                     
            KEYLIST = SIMGEN_MASTERFILE_KEYLIST_NONIa
            infile_list_NONIa = util.get_YAML_key_values(GENV,KEYLIST)
        
            # fix partial paths by adding CWD if needed.
            infile_list_SNIa  = util.fix_partial_path(infile_list_SNIa)
            infile_list_NONIa = util.fix_partial_path(infile_list_NONIa)
                             
            # store list of lists (list for each GENVERSION)
            infile_list2d_SNIa.append(infile_list_SNIa)
            infile_list2d_NONIa.append(infile_list_NONIa)

        return infile_list2d_SNIa, infile_list2d_NONIa
        
        # end sim_prep_SIMGEN_INFILE_genversion

    def sim_prep_SIMGEN_INFILE_copy(self):

        # Copy sim-input files to SIMLOGS where jobs run.
        # Pass list to generic 'copy' function so that it can
        # create a standard INPUT_FILE.LIST for future reference

        input_file         = self.config_yaml['args'].input_file 
        infile_list2d      = self.config_prep['infile_list2d']
        output_dir         = self.config_prep['output_dir']
        include_file_list_unique=self.config_prep['include_file_list_unique']
        verbose            = False

        infile_copy_list = [ input_file ] 
        for infile_list in infile_list2d:
            for infile in infile_list:
                if infile not in infile_copy_list :
                    infile_copy_list.append(infile)

        # and now the INCLUDE files; list is already trimmed to be
        # unique and avoid duplicates
        for infile in include_file_list_unique :
            infile_copy_list.append(infile)

        n_copy = len(infile_copy_list)
        base_dir = os.path.basename(output_dir)
        logging.info(f"\t Copy {n_copy} sim-input files to {base_dir}")

        util.copy_input_files(infile_copy_list, output_dir,
                              SIMGEN_INPUT_LISTFILE )

        #print(f" xxx copy_list = {infile_copy_list} ")

        # end sim_prep_SIMGEN_INFILE_copy

    def sim_prep_SIMGEN_INFILE_read(self,infile):

        # read infile (and any INCLUDE files therein) and load a
        # dictionary of keys that are needed later to
        #  + check consistency of key-values among SNIa & NONIa
        #  + check FORMAT_MASK
        #  + check user-input path
        #  + etc ...
        #
        # The key dictionary is globally defined by SIMGEN_INFILE_KEYCHECK.
        # Functions returns
        #    + input_dict dictionary of key+values
        #    + list of include files that were found and read

        key_list_include = ["INPUT_INCLUDE_FILE", "INPUT_FILE_INCLUDE"]
        key_list         = SIMGEN_INFILE_KEYCHECK # keys to read
        input_lines      = [] ;
        input_word_list  = []
        do_dump          = False

        # first make sure that infile exists
        msgerr = [ f"Check SIMGEN_INFILE_SNIa[NONIa] keys" ]
        util.check_file_exists(infile,msgerr)

        # read everything as YAML (take advantage of KEY: [VALUE] syntax)
        with open(infile, 'rt') as f :
            for line in f: 
                if util.is_comment_line(line) : continue
                if '#' in line: line = line.split('#')[0] 
                input_lines.append(line)
                input_word_list += line.split()
                #flat_word_list = [word for line in f for word in line.split()]

        #sys.exit(f"\n xxx input_lines({infile}) = \n{input_lines}\n")
        input_yaml = yaml.safe_load("\n".join(input_lines))

        # search input_word_list for include file keys the old-fashion 
        # way because this key can appear multiple times and thus can 
        # fail YAML read.
        inc_file       = ''
        inc_file_list  = []
        index_inc_list = []
        for key in key_list_include :
            key_raw = key + ':'
            index_inc_list += \
                [i for i, x in enumerate(input_word_list) if x == key_raw ]
        for indx in index_inc_list :
            inc_file = os.path.expandvars(input_word_list[indx+1])
            if inc_file not in inc_file_list: 
                inc_file_list.append(inc_file)

        if do_dump:
            print(f" 1.xxx ------------------------- ")
            print(f" 1.xxx read keys from {infile}")
            print(f" 1.xxx nlines(infile) = {len(input_lines)}")
            print(f" 1.xxx INCLUDE key indices = {index_inc_list} ")
            print(f" 1.xxx inc_file_list = {inc_file_list} ")

        # check for INCLUDE file in GENOPT
        GENOPT_GLOBAL = self.config_prep['genopt_global'].split()
        for key in key_list_include :  # check contents of GENOPT_GLOBAL
            if key in GENOPT_GLOBAL :
                j = GENOPT_GLOBAL.index(key)
                inc_file = os.path.expandvars(GENOPT_GLOBAL[j+1])
                if inc_file not in inc_file_list :
                    inc_file_list.append(inc_file)

        #print(f" xxx {infile} : include_file_list = {include_file_list} ")
        for inc_file in inc_file_list :
            util.check_file_exists(inc_file,
                                   [(f"Check INCLUDE files in {infile}")] )
            with open(inc_file, 'rt') as finc :
                for line in finc:
                    if util.is_comment_line(line) : continue
                    if '#' in line: line = line.split('#')[0]
                    input_lines.append(line)

        # read YAML again, but with all include file lines 
        input2_yaml = yaml.safe_load("\n".join(input_lines))

        # store input_dict with keys from key_list
        # If key does not exist, then it's not defined in input_dict.
        input_dict = {}
        for key in key_list:
            nkey = key_list[key][0]
            if key in input2_yaml:
                input_dict[key] = input2_yaml[key]

        if do_dump:
            print(f" 2.xxx include file = {inc_file_list}" )
            print(f" 2.xxx nlines(infile+include) = {len(input_lines)}")
            sys.exit("\n xxx DEBUG DIE xxx \n")

        return input_dict, inc_file_list

        # end sim_prep_SIMGEN_INFILE_read

    def sim_prep_SIMGEN_INFILE_verify(self,iver,keycheck,msgerr):
        # For input genversion, verify keycheck among all sim-input
        # files and abort if particular keys are not the same;
        #      e.g., REDSHIFT, PEAKMJD, SOLID_ANGLE.
        # This check avoids unphysical sims in which SNIa and SNCC 
        # have different generation ranges.
        #
        # This function determines
        #    + if keycheck should be verified
        #    + how many arguments to check,
        #    + if the key is required to exist,
        #    + if they key must be the same in each sim-input file
        #
        # Inputs :
        #    iver     = version index
        #    keycheck = element of SIMGEN_INFILE_KEYCHECK dictionary
        #
        # Output: 
        #    msgerr.append(another error message for abort)
        #    Function returns number of verify errors
        #
        # Note that ERROR messages are NOT printed here, but instead
        # they are all sent to log_assert to ensure that these error
        # message are not hidden by the python Traceback dump.
        #

        # strip off list of input files, and keys inside each file
        infile_list2d  = self.config_prep['infile_list2d']
        model_list2d   = self.config_prep['model_list2d']
        INFILE_KEYS    = self.config_prep['INFILE_KEYS']
        INFILE_PAIRS   = self.config_prep['INFILE_PAIRS_VERIFY_ERROR']
        n_infile       = len(infile_list2d[iver])
        narg           = SIMGEN_INFILE_KEYCHECK[keycheck][0]
        do_require     = SIMGEN_INFILE_KEYCHECK[keycheck][1]
        do_verify      = SIMGEN_INFILE_KEYCHECK[keycheck][2]

        nerr   = 0 
        nfound = 0

        for ifile in range(0,n_infile):
            infile_ref       = infile_list2d[iver][0] 
            infile_tmp       = infile_list2d[iver][ifile] 
            key_exists       = False
            infile_pair      = f"{infile_ref}+{infile_tmp}+{keycheck}"
            unique_pair      = infile_pair not in INFILE_PAIRS

            # first make check on required keys
            if keycheck in INFILE_KEYS[iver][ifile] :
                nfound      += 1
                key_exists   = True 
                    
            elif do_require :
                nerr      += 1
                key_exists = False 
                msg=(f"ERROR: required key {keycheck} missing in {infile_ref}")
                msgerr.append(msg)

            #print(f" xxx key={keycheck} exist={key_exists} req={do_require}")

            # next, check that key values are the same in all files
            if key_exists and do_verify and unique_pair :
                key_value_ref = str(INFILE_KEYS[iver][0][keycheck])
                key_value_tmp = str(INFILE_KEYS[iver][ifile][keycheck])
                verify        = True
                for i in range(0,narg):
                    val_ref = key_value_ref.split()[i]
                    val_tmp = key_value_tmp.split()[i]
                    if val_tmp != val_ref  :
                        verify = False

                if not verify :
                    nerr += 1 ; 
                    msg_ref = f"{keycheck}: {key_value_ref} in {infile_ref}"
                    msg_tmp = f"{keycheck}: {key_value_tmp} in {infile_tmp}"
                    msg=(f"KEY-VERIFY ERROR:\n\t {msg_ref}\n\t {msg_tmp}")
                    msgerr.append(msg)
                    INFILE_PAIRS.append(infile_pair)

#         print(f" xxx check '{keycheck}' for {iver}: nfound={nfound}" )
        return nerr

        # end sim_prep_SIMGEN_INFILE_verify

    def sim_prep_PATH_SNDATA_SIM(self):
        # figure out where sim data files will be written
        CONFIG = self.config_yaml['CONFIG']
        key       = 'PATH_SNDATA_SIM'
        if key in CONFIG:
            path_sndata_sim = os.path.expandvars(CONFIG[key])
            if os.path.exists(path_sndata_sim) == False :
                msgerr = []
                msgerr.append(f"PATH_SNDATA_SIM does not exist:")
                msgerr.append(f"  {path_sndata_sim} ")
                self.log_assert(False,msgerr)
            flag = True
        else:
            path_sndata_sim = f"{SNDATA_ROOT}/SIM"
            flag = False

        self.config_prep['path_sndata_sim']      = path_sndata_sim
        self.config_prep['user_path_sndata_sim'] = flag

        # end sim_prep_PATH_SNDATA_SIM

    def write_command_file(self, icpu, f):
        # For this icpu, write full set of sim commands to 
        # already-opened command file with pointer f.
        # Function returns number of jobs for this cpu

        n_job_tot      = self.config_prep['n_job_tot']
        n_job_split    = self.config_prep['n_job_split']
        n_core         = self.config_prep['n_core']

        # open CMD file for this icpu
        # xxx mark delete  f = open(COMMAND_FILE, 'a') 

        n_job_cpu    = 0     # number of jobs for this CPU
        iver_list    = self.config_prep['iver_list']
        ifile_list   = self.config_prep['ifile_list']
        isplit_list  = self.config_prep['isplit_list']

        # keep track of TMP version names for MERGE process
        if icpu == 0 :
            n_genv      = self.config_prep['n_genversion']
            list2d = [['' for j in range(0,10)] for i in range(0,n_genv)]
            self.config_prep['TMP_genversion'] = list2d
    
        TMP_list2d = self.config_prep['TMP_genversion']  # init array

        # loop over ALL jobs, and pick out the ones for this ICPU
        n_job_local = 0
        for iver,ifile,isplit in zip(iver_list,ifile_list,isplit_list):

            index_dict = {
                'iver':iver, 'ifile':ifile, 'isplit':isplit, 'icpu':icpu
            }  
            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :
                n_job_cpu += 1

                # define sim job and merge job; then glue together
                job_info_sim   = self.prep_JOB_INFO_sim(index_dict)
                util.write_job_info(f, job_info_sim, icpu)

                job_info_merge = \
                    self.prep_JOB_INFO_merge(icpu,n_job_local,False) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

                # store TMP_VERSION for later
                TMP_list2d[iver][ifile] = job_info_sim['tmp_genversion']

        # store TMP version strings needed later in MERGE.LOG file
        self.config_prep['TMP_genversion_list2d'] = TMP_list2d

        # xxx mark delete xxx  f.close()

        if n_job_local != n_job_tot :
            msgerr = []
            msgerr.append(f"Expected {n_job_tot} total jobs;")
            msgerr.append(f"but found {n_job_local} jobs.")
            self.log_assert(False,msgerr)

        return n_job_cpu
        # end write_command_file for sim
        

    def prep_JOB_INFO_sim(self, job_index_dict):

        # Return JOB_INFO dictionary with 
        #    cd job_dir
        #    program.exe arg_list  > log_file
        #    touch TMP_[xxx].DONE
        #
        # Inputs
        #    job_index_dict = dictionary of indices for this job
        #

        # strip off indices from input dictionary
        iver   = job_index_dict['iver']
        ifile  = job_index_dict['ifile']
        isplit = job_index_dict['isplit']
        icpu   = job_index_dict['icpu']

        # pick off a few globals
        CONFIG            = self.config_yaml['CONFIG']
        GENPREFIX         = CONFIG['GENPREFIX']

        args              = self.config_yaml['args']
        no_merge          = args.nomerge
        kill_on_fail      = args.kill_on_fail
        check_abort       = args.check_abort

        program           = self.config_prep['program'] 
        n_job_split       = self.config_prep['n_job_split']
        output_dir        = self.config_prep['output_dir']
        infile_list2d     = self.config_prep['infile_list2d']
        model_list2d      = self.config_prep['model_list2d']
        INFILE_KEYS       = self.config_prep['INFILE_KEYS']
        n_genversion      = self.config_prep['n_genversion']
        genversion_list   = self.config_prep['genversion_list']
        genopt_list2d     = self.config_prep['genopt_list2d']
        ngentot_list2d    = self.config_prep['ngentot_list2d']
        ranseed_list      = self.config_prep['ranseed_list']
        ranseed_key       = self.config_prep['ranseed_key']
        IS_CHANGE         = 'CHANGE' in ranseed_key
        n_include_files   = self.config_prep.get('n_include_files', 0)
        genopt_global     = self.config_prep['genopt_global']
        user_path_sndata  = self.config_prep['user_path_sndata_sim']
        path_sndata       = self.config_prep['path_sndata_sim']
        format_mask       = self.config_prep['format_mask']
        
        Nsec  = seconds_since_midnight

        reset_cidoff     = self.config_prep['reset_cidoff']
        cidran_min       = self.config_prep['cidran_min']
        cidran_max_list  = self.config_prep['cidran_max_list']
        cidoff_list3d    = self.config_prep['cidoff_list3d']

        # init JOB_INFO dictionary. Note that sim job runs in same
        # dir where simgen-master file resides; this avoids copying
        # lots of iinputs & nclude files and potential infile-clobber
        # issues.

        JOB_INFO = {}
        JOB_INFO['job_dir']   = output_dir  # where to run job
        JOB_INFO['program']   = program

        isplit1      = isplit+1               # for TMP-genversion names 
        genversion   = genversion_list[iver]
        genopt       = genopt_list2d[iver][ifile]
        ranseed      = ranseed_list[isplit]
        infile       = infile_list2d[iver][ifile]
        model        = model_list2d[iver][ifile]
        ngentot      = ngentot_list2d[iver][ifile]
        if check_abort: ngentot = NGENTOT_CHECK_ABORT
        Nsec         = seconds_since_midnight

        split_string = f"{isplit1:04d}"         # e.g., 0010
        model_string = self.model_string_suffix(model,ifile)

        tmp1       = f"TMP_{USER4}_{genversion}_{model_string}" 
        tmp2       = self.genversion_split_suffix(isplit1,Nsec)

        tmp_ver    = f"{tmp1}-{tmp2}"        # temp GENVERSION
        log_file   = f"{tmp_ver}.LOG"
        done_file  = f"{tmp_ver}.DONE"
        genprefix  = f"{GENPREFIX}_{model_string}-{split_string}"
        
        arg_list    = []
        arg_list.append(f"    GENVERSION {tmp_ver}")
        arg_list.append(f"    GENPREFIX  {genprefix}")
        arg_list.append(f"    NGENTOT_LC    {ngentot}    NGEN_LC 0")
        arg_list.append(f"    RANSEED {ranseed}")
        arg_list.append(f"    FORMAT_MASK {format_mask}")
        

        if reset_cidoff > 0 :
            cidoff = cidoff_list3d[iver][ifile][isplit]
            str1   = f"CIDRAN_MIN {cidran_min}"
            str2   = f"CIDRAN_MAX {cidran_max_list[iver]}"
            str3   = f"CIDOFF {cidoff}"
            arg_list.append(f"    {str1}    {str2}    {str3}")

        arg_list.append(f"    JOBID {isplit1}     NJOBTOT {n_job_split}")
        arg_list.append(f"    WRFLAG_MODELPAR  0") # disable model-par output
        arg_list.append(f"    WRFLAG_YAML_FILE 1") # enable YAML output

        # check for user-option to require DOCANA
        if self.config_yaml['args'].require_docana :
            arg_list.append(f"    REQUIRE_DOCANA 1")

        if user_path_sndata :                 
            arg_list.append(f"    PATH_SNDATA_SIM {path_sndata}")
        
        
        # suppress this if PATH_USER_INPUT already defined by user ??
        key = 'PATH_USER_INPUT'
        if key not in INFILE_KEYS[iver][ifile] :
            arg_list.append(f"    {key} {CWD}")

        # May 2022: check option to sync ransystpar random seeds
        #   Note that genopt may be modified.
        genopt = self.sync_ransystpar_ranseed_syst(genopt,isplit1)
        #sys.exit(f"\n genopt= {genopt}")

        arg_list.append(f"{genopt}")         # user args by version
        arg_list.append(f"{genopt_global}")     # user global args
        if IS_CHANGE and n_include_files > 0:
            input_include_file = self.config_prep['include_list'][isplit]
            arg_list.append(f"INPUT_INCLUDE_FILE {input_include_file}")

        JOB_INFO['input_file']    = infile
        JOB_INFO['log_file']      = log_file
        JOB_INFO['done_file']     = done_file
        JOB_INFO['arg_list']      = arg_list
        JOB_INFO['tmp_genversion_split']  = tmp_ver
        JOB_INFO['tmp_genversion']        = tmp1    # combined genv        
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['check_abort']   = check_abort

        return JOB_INFO

        # end prep_JOB_INFO_sim

    def model_string_suffix(self,model,ifile):
        model_string = f"{model}MODEL{ifile}"     # e.g., SNIaMODEL0
        return model_string
        # end model_string_suffix

    def genversion_split_suffix(self,isplit,Nsec):
        suffix = f"{isplit:04d}_{Nsec}"
        return suffix

    def sync_ransystpar_ranseed_syst(self,genopt,isplit1):

        # if genopt contains "RANSYSTPAR_RANSEED_SYST <ranseed>",
        # then replace ranseed with value that depends on isplit1
        # so that all ISPLIT=1 have same seed, all ISPLIT=2 have same
        # seed, etc ... Goal is to sync correlated systematics among
        # samples; e.g., galactic extinction, cosmology params ...
        # Note that argument genopt is modified.

        key = "RANSYSTPAR_RANSEED_SYST"
        genopt_out = genopt  
        if key in genopt:
            genopt_list   = genopt.split()
            j             = genopt_list.index(key)
            ranseed_orig  = int(genopt_list[j+1])
            ranseed_split = ranseed_orig*7 + isplit1*41 

            str_orig      = str(ranseed_orig)
            str_split     = str(ranseed_split)
            genopt_out    = genopt.replace(str_orig,str_split)
            print(f"\t ransystpar_ranseed_syst -> {ranseed_split}" \
                  f" for ISPLIT={isplit1}")
            
        return genopt_out
        
        # end sync_ransystpar_ranseed_syst

    def append_info_file(self,f):

        # append sim-specific information to SUBMIT.INFO file

        CONFIG            = self.config_yaml['CONFIG']
        simlog_dir        = self.config_prep['output_dir']
        script_dir        = self.config_prep['script_dir']
        path_sndata_sim   = self.config_prep['path_sndata_sim']
        n_genversion      = self.config_prep['n_genversion']
        ngentot_sum       = self.config_prep['ngentot_sum']
        format_mask       = self.config_prep['format_mask']
        ranseed_key       = self.config_prep['ranseed_key'] 
        ngen_unit         = self.config_prep['ngen_unit']
        n_job_split       = self.config_prep['n_job_split']

        # - - - - - - - 
        f.write("\n# Original user input \n")

        comment = "1 unit per GENRANGE_(PEAKMJD,REDSHIFT) + SOLID_ANGLE"
        f.write(f"NGEN_UNIT:            {ngen_unit}      # ({comment})\n")

        f.write(f"PATH_SNDATA_SIM:      {path_sndata_sim} \n")

        f.write(f"SIMLOG_DIR:           {simlog_dir} \n")

        f.write(f"GENPREFIX:            {CONFIG['GENPREFIX']} \n")

        f.write(f"RANSEED_KEY:          {ranseed_key}  \n")

        f.write(f"FORMAT_MASK:          {format_mask} \n")

        # - - - - - - - 
        f.write("\n# Computed from original input \n")

        f.write(f"JOBFILE_WILDCARD:     'TMP_{USER4}*' \n")

        f.write(f"N_GENVERSION:         {n_genversion} \n")

        comment = "NGEN/split summed over all GENVERSIONs"
        f.write(f"NGENTOT_SUM:        {ngentot_sum}  # ({comment})\n")

        # write grand total over all cores (July 26 2021)
        ngentot_sumsplit = ngentot_sum * n_job_split
        comment = f"NGEN x {n_job_split} split-jobs " \
                  f"summed over all GENVERSIONs"
        f.write(f"NGENTOT_SUMSPLIT:   {ngentot_sumsplit}  # ({comment})\n")
    
        #end append_info_file

    def create_merge_table(self,f):

        # Required element of submit process. Before submitting jobs,
        # create initial merge file with all WAIT states.
        # This file is read and updated frequently by merge
        # process invoked by -m argument to submit_batch_jobs.sh .
        # A locally defined MERGE_INFO structure is passed to 
        # a generic write_MERGE_INFO function to create MERGE.LOG/
        # Uses YAML format, and for human-readability there is a 
        # one line commented header before each YAML table.

        n_job_tot            = self.config_prep['n_job_tot']
        n_job_split          = self.config_prep['n_job_split']
        simlog_dir           = self.config_prep['output_dir']
        path_sndata_sim      = self.config_prep['path_sndata_sim']
        n_genversion         = self.config_prep['n_genversion']
        genversion_list      = self.config_prep['genversion_list']
        genversion_list_all  = self.config_prep['genversion_list_all']
        igenver_list_all     = self.config_prep['igenver_list_all']
        infile_list2d        = self.config_prep['infile_list2d']
        model_list2d         = self.config_prep['model_list2d']
        ranseed_key          = self.config_prep['ranseed_key']
        IS_REPEAT            = 'REPEAT' in ranseed_key
        IS_CHANGE            = 'CHANGE' in ranseed_key

        n_all = len(genversion_list_all)

        TMP_genversion_list2d  = self.config_prep['TMP_genversion_list2d']

        # create SPLIT table
        # write TMP_ genversions per SNIa/NONIa model
        header_line = " STATE  IVER     GENVERSION     " \
                      "NLC_GEN NLC_WRITE NSPEC_WRITE  CPU  NSPLIT"
        MERGE_INFO = { 
            'primary_key' : TABLE_SPLIT,
            'header_line' : header_line,
            'row_list'      : []
        }
        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for iver in range(0,n_genversion):
            n_file     = len(infile_list2d[iver]) 
            genversion = genversion_list[iver]
            for ifile in range(0,n_file):
                TMP_genv = TMP_genversion_list2d[iver][ifile]
                # define ROW here is fragile in case columns are changed
                ROW = [ STATE, iver, TMP_genv, 0, 0, 0, 0, n_job_split ]
                MERGE_INFO['row_list'].append(ROW)    
        util.write_merge_file(f, MERGE_INFO, [] ) 

        # - - - - - - 
        # finally, the combined versions (remove NSPLIT column)
        header_line = " STATE     IVER  GENVERSION    " \
                      " NLC_GEN NLC_WRITE NSPEC_WRITE CPU"
        MERGE_INFO = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line,
            'row_list'      : []
        }

        # for combined versions, use genversion_list_all to account
        # RANSEED_REPEAT or RANSEED_CHANGE
        for iver_all in range(0,n_all):
            genversion = genversion_list_all[iver_all]
            iver       = igenver_list_all[iver_all]
            ROW = [ STATE, iver, genversion, 0, 0, 0, 0 ]
            MERGE_INFO['row_list'].append(ROW)    
        util.write_merge_file(f, MERGE_INFO, [] )

        #printf(" xxx force crash with C-like printf xxx \n")

        # end create_merge_table

    def genversion_expand_list(self,genversion_list,ranseed_key,n_job_split):

        # Define genversion_list_all to include split jobs for RANSEED_CHANGE.
        # Example with GENVERSIONS test1 and test2, and 3 split jobs with
        # RANSEED_CHANGE:
        #  Inputs is
        #     genversion_list = ['test1', test2' ] 
        #
        #    output is
        #      genversion_list_all =            
        #       [ 'test1-0001' test1-0002', 'test1-0003',
        #         'test2-0001' test2-0002', 'test2-0003' ] (determined here)
        #
        #      iver_list_all = 0,0,0, 1,1,1
        #
        genversion_list_all  = []
        iver_list_all        = []
        n_ver                = len(genversion_list)

        if 'CHANGE' in ranseed_key :
            for iver in range(0,n_ver) :
                genv = genversion_list[iver]
                for isplit in range(1,n_job_split+1) : # split starts at 1, not 0
                    genv_split = self.genversion_split_name(genv,isplit)
                    genversion_list_all.append(genv_split)
                    iver_list_all.append(iver)
        else:
            genversion_list_all = genversion_list # no change for RANSEED_REPEAT
            iver_list_all        = list(range(n_ver))
        
        return genversion_list_all,iver_list_all

        # end genversion_expand_list

    def genversion_split_name(self,genversion,isplit):
        # if genversion = ABC and isplit=4, return ABC-0004
        name = f"{genversion}-{isplit:04d}"
        return name

# ==============================================
#        SIM-FILE MERGE UTILS
# ==============================================

    def merge_config_prep(self,output_dir):

        # sim-specific settings to config_prep that are needed later.
        submit_info_yaml = self.config_prep['submit_info_yaml'] 

        self.config_prep['path_sndata_sim']     = \
                        submit_info_yaml['PATH_SNDATA_SIM']
        self.config_prep['output_dir']     = output_dir 

        self.sim_prep_GENOPT_GLOBAL()
        self.sim_prep_GENVERSION_LIST()
        self.sim_prep_SIMGEN_INFILE()

        # end merge_config_prep

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # Check for genversion(s) where all DONE files exist;
        # in this case, update both the split and merge tables,
        # and also move data files. Watch logic dependence on
        # RANSEED_REPEAT vs. RANSEED_CHANGE.
        #
        # Inputs:
        #    MERGE_INFO_CONTENTS : contents of MERGE log file
        #

        row_list_split   = MERGE_INFO_CONTENTS[TABLE_SPLIT]
        row_list_merge   = MERGE_INFO_CONTENTS[TABLE_MERGE]

        submit_info_yaml = self.config_prep['submit_info_yaml']
        simlog_dir       = submit_info_yaml['SIMLOG_DIR']
        Nsec_time_stamp  = submit_info_yaml['TIME_STAMP_NSEC']
        ranseed_key      = submit_info_yaml['RANSEED_KEY']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        n_genversion     = submit_info_yaml['N_GENVERSION']
        msgerr           = []
        COLNUM_STATE     = COLNUM_SIM_MERGE_STATE
        COLNUM_IVER      = COLNUM_SIM_MERGE_IVER
        COLNUM_GENV      = COLNUM_SIM_MERGE_GENVERSION
        COLNUM_NSPLIT    = COLNUM_SIM_MERGE_NSPLIT
        COLNUM_NGEN      = COLNUM_SIM_MERGE_NLC_GEN
        COLNUM_NLC       = COLNUM_SIM_MERGE_NLC_WRITE
        COLNUM_NSPEC     = COLNUM_SIM_MERGE_NSPEC_WRITE
        COLNUM_CPU       = COLNUM_SIM_MERGE_CPU

        # define keys to read and sum from YAML file produced by science job
        key_nlc_gen, key_nlc_gen_sum, key_nlc_gen_list = \
                    self.keynames_for_job_stats('NGENLC_TOT')

        key_nlc_write, key_nlc_write_sum, key_nlc_write_list = \
                    self.keynames_for_job_stats('NGENLC_WRITE')

        key_nspec_write, key_nspec_write_sum, key_nspec_write_list = \
                    self.keynames_for_job_stats('NGENSPEC_WRITE')

        key_cpu, key_cpu_sum, key_cpu_list = \
                    self.keynames_for_job_stats('CPU_MINUTES')

        KEY_YAML_LIST   = [ key_nlc_gen, key_nlc_write, key_nspec_write, 
                            key_cpu ]
        # xxxx 'SURVEY', 'IDSURVEY' ]

        # init outputs of function
        n_state_change    = 0
        row_split_new     = []
        row_merge_new     = []

        # init new combined rows to be the input ... will change
        # below if DONE files are found
        for row in row_list_merge:
            row_merge_new.append(row) # default output is same as input

        irow_split = 0
        for row in row_list_split:
            row_split_new.append(row)  # default output is same as input
            STATE        = row[COLNUM_STATE]
            IVER         = row[COLNUM_IVER]
            TMP_GENV     = row[COLNUM_GENV] + "*" + str(Nsec_time_stamp)
            NSPLIT       = row[COLNUM_NSPLIT]

            # if not done, update the STATE
            Finished = (STATE == SUBMIT_STATE_DONE) or \
                       (STATE == SUBMIT_STATE_FAIL)

            if not Finished :

                TMP_LOG_LIST, TMP_DONE_LIST, TMP_YAML_LIST = \
                    util.get_file_lists_wildcard(simlog_dir,TMP_GENV)

                # DONE and YAML lists are forced to have same length 
                # as LOG list, 
                # careful to sum only the files that are NOT None 
                NLOG   = sum(x is not None for x in TMP_LOG_LIST)
                NDONE  = sum(x is not None for x in TMP_DONE_LIST)
                NYAML  = sum(x is not None for x in TMP_YAML_LIST)
                NEW_STATE      = STATE

                if NLOG > 0:
                    NEW_STATE = SUBMIT_STATE_RUN
                if NDONE == NSPLIT :
                    NEW_STATE = SUBMIT_STATE_DONE

                    # update sim stats and check for errors
                    job_stats = self.get_job_stats(simlog_dir, 
                                                   TMP_LOG_LIST, TMP_YAML_LIST,
                                                   KEY_YAML_LIST)

                    #print(f"\n xxx KEY_YAML_LIST = {KEY_YAML_LIST} ")
                    #print(f" xxx job_stats = \n {job_stats} \n")

                    # check for failures
                    nfail = job_stats['nfail']
                    if nfail > 0 : NEW_STATE = SUBMIT_STATE_FAIL

                    # update stats for SPLIT table; same for REPEAT or CHANGE
                    row[COLNUM_NGEN]   = job_stats[key_nlc_gen_sum]
                    row[COLNUM_NLC]    = job_stats[key_nlc_write_sum]
                    row[COLNUM_NSPEC]  = job_stats[key_nspec_write_sum]
                    row[COLNUM_CPU]    = job_stats[key_cpu_sum]
                    row_split_new[irow_split] = row     # update split stats

                    # move data files and updat combine table depending
                    # on RANSEED_REPEAT or RANSEED_CHANGE
                    if 'REPEAT' in ranseed_key :
                        row_merge  = row_list_merge[IVER]
                        GENV_MERGE = row_merge[COLNUM_GENV]
                        self.move_sim_data_files(TMP_GENV,GENV_MERGE, nfail)
                        row_merge[COLNUM_NGEN] += job_stats[key_nlc_gen_sum]
                        row_merge[COLNUM_NLC]  += job_stats[key_nlc_write_sum]
                        row_merge[COLNUM_NSPEC]+= job_stats[key_nspec_write_sum]
                        row_merge[COLNUM_CPU]  += job_stats[key_cpu_sum]
                        row_merge_new[IVER] = row_merge
                    else:
                        for isplit in range(0,n_job_split):
                            Nsec   = Nsec_time_stamp
                            suffix = self.genversion_split_suffix(isplit+1,Nsec)
                            iver_all    = isplit + IVER*n_job_split
                            row_merge   = row_list_merge[iver_all]
                            tmp_genv    = f"{row[COLNUM_GENV]}-{suffix}"
                            genv_merge  = row_merge[COLNUM_GENV]

                            self.move_sim_data_files(tmp_genv,genv_merge,nfail)

                            row_merge[COLNUM_NGEN] += \
                                    job_stats[key_nlc_gen_list][isplit]
                            row_merge[COLNUM_NLC] += \
                                    job_stats[key_nlc_write_list][isplit]
                            row_merge[COLNUM_NSPEC] += \
                                    job_stats[key_nspec_write_list][isplit]
                            row_merge[COLNUM_CPU] += \
                                    job_stats[key_cpu_list][isplit]
                            row_merge_new[iver_all] = row_merge

                if NEW_STATE != STATE :
                    row[COLNUM_SIM_MERGE_STATE] = NEW_STATE
                    n_state_change += 1

            irow_split += 1

        # - - - - - - - - - - - - - - - - - - - - - - - -
        # Update DONE states for MERGE    table.
        # This is tricky because all split-job MODELS must be done
        # before declaring MERGE job to be done.

        # Check which split jobs are done for all models,
        # and also which jobs are running for any model
        iver_done_flag    = [ True ] * n_genversion
        iver_run_flag     = [ False] * n_genversion
        iver_fail_flag    = [ False] * n_genversion
        for row in row_split_new :
            cpu             = float(row[COLNUM_CPU])
            row[COLNUM_CPU] = float(f"{cpu:.1f}")    # e.g., 45.2

            STATE     = row[COLNUM_STATE]
            iver      = row[COLNUM_IVER]
            if STATE != SUBMIT_STATE_DONE :
                iver_done_flag[iver] = False # at least 1 model NOT done
            if STATE == SUBMIT_STATE_RUN :
                iver_run_flag[iver]     = True     # at least 1 model running
            if STATE == SUBMIT_STATE_FAIL :
                iver_fail_flag[iver]  = True  # any split-failure fails merge

        # check which merged genversions are done
        iver_all = 0
        for row in row_merge_new:
            cpu             = float(row[COLNUM_CPU])
            row[COLNUM_CPU] = float(f"{cpu:.1f}")    # e.g., 45.2

            iver = row[COLNUM_IVER] 
            if iver_done_flag[iver] :  # merged genversion is ready
                row_merge_new[iver_all][COLNUM_STATE] = SUBMIT_STATE_DONE
            if iver_run_flag[iver] :  # at least 1 model is running
                row_merge_new[iver_all][COLNUM_STATE] = SUBMIT_STATE_RUN
            if iver_fail_flag[iver] : # something failed
                row_merge_new[iver_all][COLNUM_STATE] = SUBMIT_STATE_FAIL
            iver_all += 1

        row_list_dict = {
            'row_split_list' : row_split_new,
            'row_merge_list' : row_merge_new,
            'row_extra_list' : []
        }
        return row_list_dict, n_state_change

        # xxx mark return row_split_new, row_merge_new, n_state_change

        # end merge_update_state

    def move_sim_data_files(self,genversion_split, genversion_combine, nfail):

        # Move sim data files from genversion_split to genversion_combine,
        # and gzip FITS files. Genversion_split may include wildcard.
        # Also update auxilliary files.
        # If genversion_combine directory does not exist, then create it.
        # Note that submit task must ensure that prevous/lingering
        # genversion_combine is removed so that here we know to 
        # create a new directory if it doesn't exist.
        #
        # If nfail > 0, create genversion_combine and write FAIL 
        # in README file; the quit. This allows downstream analysis codes 
        # to quickly check for FAIL.
        #
        # Dec 20 2021: gzip FITS if not already gzipped.

        args             = self.config_yaml['args']
        if args.check_abort: return

        submit_info_yaml = self.config_prep['submit_info_yaml']
        simlog_dir       = submit_info_yaml['SIMLOG_DIR']
        path_sndata_sim  = submit_info_yaml['PATH_SNDATA_SIM']    
        Nsec_time_stamp  = submit_info_yaml['TIME_STAMP_NSEC']

        from_dir   = f"{path_sndata_sim}/{genversion_split}"  
        target_dir = f"{path_sndata_sim}/{genversion_combine}"

        msg = f"  move {genversion_split} files to {genversion_combine}"
        logging.info(msg)

        dump_split_list     = glob.glob(f"{from_dir}/TMP*.DUMP")

        # defin aux files for combined version
        dump_file     = f"{target_dir}/{genversion_combine}.DUMP"
        readme_file   = f"{target_dir}/{genversion_combine}.README"
        list_file     = f"{target_dir}/{genversion_combine}.LIST"

        # if target dir does NOT exist, create target dir along
        # with aux files.
        if os.path.exists(target_dir) == False :
            os.mkdir(target_dir)

            # create blank IGNORE file
            # xxx mark delete with open (ignore_file,"w") as f :   pass

            # create blank README file
            with open (readme_file,"w") as f :
                if nfail > 0 : 
                    f.write("FAIL\n")  # leave message for snlc_fit
                pass

            # create combined DUMP file with VARNAMES & comments from
            # first DUMP file. Protect against job failure.
            if len(dump_split_list) > 0 :
                dump_file_template = dump_split_list[0]
                self.create_simgen_dump_file(dump_file_template,dump_file)

        # if there were failures, return
        if nfail > 0 : 
            return

        # - - - - - - - - - - 
        # Move the FITS files, update [VERSION].LIST file, gzip FITS files.
        # Make sure to list the SNIa first, then NONIa, so that analysis 
        # init is based on SNIa. 2>/dev/null suppresses linux error
        # message when SNIa or NONIa FITS files are not there.
        # 

        cd_dir    = f"cd {target_dir}"
        mv_FITS   = f"mv {from_dir}/*.FITS* ."
        
        tmp_FITS  = f"*{MODEL_SNIa}MODEL*HEAD.FITS*"
        ls_SNIa   = f"ls {tmp_FITS}  >     {list_file} 2>/dev/null"

        tmp_FITS  = f"*{MODEL_NONIa}MODEL*HEAD.FITS*"
        ls_NONIa  = f"ls {tmp_FITS} >> {list_file} 2>/dev/null"

        ls_LIST    = f"{ls_SNIa};  {ls_NONIa}"

        #remove .gz extensions in LIST file
        rm_gz      = f"sed -i 's/FITS.gz/FITS/g' {list_file}"

        cmd     = f"{cd_dir}; {mv_FITS}; {ls_LIST}; {rm_gz} "

        # gzip FITS files if not already gzipped. 
        fits_list = glob.glob1(target_dir,"*.FITS")
        if len(fits_list) > 0 :  cmd += "; gzip *.FITS"

        os.system(cmd)

        # loop over TMP_*DUMP files and append combined DUMP file
        for dump_split_file in dump_split_list :
            self.append_merge_dump_file(dump_split_file,dump_file)

        # end move_sim_data_files


    def append_merge_dump_file(self,dump_split_file,dump_file):

        # read input dump_split_file and save the lines
        # with "SN:" key. Append these SN lines to already
        # existing dump_file.

        # make sure both dump files exist
        msgerr = []
        msgerr.append(f"Check TMP-DUMP files")
        util.check_file_exists(dump_split_file,msgerr)

        msgerr.append(f"Combine DUMP file not created??")
        util.check_file_exists(dump_file,msgerr)

        # - - - - 
        lines_SN = []
        with open (dump_split_file,"r") as f :
            for line in f:
                word_list = line.split()
                if len(word_list) > 0 :
                    if word_list[0] == 'SN:' :
                        lines_SN.append(line)

#         nline_SN = len(lines_SN)
#         print(f" xxx grab {nline_SN} lines from {dump_split_file} ")

        with open (dump_file,"a") as f :
            f.write("".join(lines_SN) )

        # end update_merge_dump_file

    def create_simgen_dump_file(self,dump_file_template,dump_file):

        # dump_file_template is an existing dump file from which
        # strip out comment lines and varnames. This info is 
        # written to the top of the new "dump_file".

        dump_comment_lines = []
        dump_varnames_line = ""
        nline_read           = 0
        with open (dump_file_template,"r") as f :
            for line in f:
                if len(line.strip()) > 1 :
                    nline_read += 1
                    if line[0] == '#' :
                        dump_comment_lines.append(line.rstrip("\n"))
                    if line.split()[0] == 'VARNAMES:' :
                        dump_varnames_line = line
                    if line.split()[0] == 'SN:' :
                        break

        with open (dump_file,"w") as f :
            for line in dump_comment_lines :
                f.write(f"{line}\n")
            f.write(f"\n{dump_varnames_line}\n")

        # end create_simgen_dump_file

    def merge_job_wrapup(self,iver_all,MERGE_INFO_CONTENTS):

        # GENVERSION index iver_all is done, so perform 
        # wrap up tasks:
        #    + write to README
        #    + create misc/ subdir; move & copy files to misc/
        #    + remove TMP* GENVERSIONs
        #
        # Input 'iver_all' is a reminder that this function works
        # for both RANSEED_REPEAT and RANSEED_CHANGE

        input_file          = self.config_yaml['args'].input_file
        submit_info_yaml    = self.config_prep['submit_info_yaml']
        path_sndata_sim     = submit_info_yaml['PATH_SNDATA_SIM'] 
        simlog_dir          = submit_info_yaml['SIMLOG_DIR'] 
        ranseed_key         = submit_info_yaml['RANSEED_KEY'] 
        cleanup_flag        = submit_info_yaml['CLEANUP_FLAG']
        Nsec_time_stamp     = submit_info_yaml['TIME_STAMP_NSEC']          
        Nsec_now  = seconds_since_midnight # current time since midnight

        row_list_merge  = MERGE_INFO_CONTENTS[TABLE_MERGE]
        row_list_split  = MERGE_INFO_CONTENTS[TABLE_SPLIT]
        genversion      = row_list_merge[iver_all][COLNUM_SIM_MERGE_GENVERSION]
        iver            = row_list_merge[iver_all][COLNUM_SIM_MERGE_IVER]
        path_genv       = f"{path_sndata_sim}/{genversion}"

        IS_REPEAT = 'REPEAT' in ranseed_key  # do 2nd stage merge
        IS_CHANGE = 'CHANGE' in ranseed_key  # do NOT do 2nd stage

        misc_subdir   = "misc"
        misc_dir      = f"{path_genv}/{misc_subdir}"

        # - - - - - - - -
        msg = f"    Perform wrap-up tasks for {genversion} ({Nsec_now})"
        logging.info(msg)

        # - - - - - - - - 
        # create and write README file
        msg = f"\t Create README for {genversion}"
        logging.info(msg)
        readme_file   = f"{path_genv}/{genversion}.README"
        with open(readme_file,"w") as f : 
            self.merge_write_readme(f, iver_all, MERGE_INFO_CONTENTS)

        # move README files to misc/ for REPEAT only
        os.mkdir(misc_dir)
        TMP_GENV_LIST = []
        for row in row_list_split :
            if iver == row[COLNUM_SIM_MERGE_IVER] :
                if IS_CHANGE: 
                    suffix = f"-{genversion[-4:]}" # strip split-job number
                else:
                    suffix = ""

                TMP_GENV    = row[COLNUM_SIM_MERGE_GENVERSION] + suffix
                TMP_GENV   += f"*{Nsec_time_stamp}"
                TMP_GENV_LIST.append(TMP_GENV)
                readme_list = f"{path_sndata_sim}/{TMP_GENV}/TMP*README"
                mv_readme   = f"mv {readme_list} {misc_dir}/"
                os.system(mv_readme)

        # - - - - - - - - - - 
        # copy input files to misc/
        if IS_REPEAT:
            infile_list_2copy = [ f"{simlog_dir}/{input_file}" ]
            infile_list2d     = self.config_prep['infile_list2d']
            n_file            = len(infile_list2d[iver])
            for ifile in range(0,n_file) :
                infile      = infile_list2d[iver][ifile]
                infile_base = os.path.basename(infile)
                infile_list_2copy.append(f"{simlog_dir}/{infile_base}")

            util.copy_input_files(infile_list_2copy,misc_dir,"")

        # - - - - - - - - - - 
        # remove TMP versions for this iver ... if CLEANUP_FLAG is set
        # Make sure to include Nsec stamp in dir name to avoid clobbering
        # other jobs in progress.
        if cleanup_flag :
            # first tar up misc dir
            tar      = f"tar -cf {misc_subdir}.tar {misc_subdir} "
            gzip     = f"gzip {misc_subdir}.tar"
            rm       = f"rm -r {misc_subdir}"
            clean    = f"cd {path_genv}; {tar}; {gzip}; {rm}"
            os.system(clean)

            msg = f"\t Remove TMP_ GENVERSIONs"
            logging.info(msg)
            for TMP_GENV in TMP_GENV_LIST:
                rm_TMP      = f"rm -rf {path_sndata_sim}/{TMP_GENV} "
                os.system(rm_TMP)

        #if iver_all == 1 :   sys.exit("\n xxx DEBUG DIE xxx \n")

        # end merge_job_wrapup

    def merge_write_readme(self, f, iver_all, MERGE_INFO_CONTENTS):

        input_file          = self.config_yaml['args'].input_file
        submit_info_yaml    = self.config_prep['submit_info_yaml']
        infile_list2d       = self.config_prep['infile_list2d']
        model_list2d        = self.config_prep['model_list2d']
        path_sndata_sim     = submit_info_yaml['PATH_SNDATA_SIM'] 
        simlog_dir          = submit_info_yaml['SIMLOG_DIR'] 
        ranseed_key         = submit_info_yaml['RANSEED_KEY'] 
        cleanup_flag        = submit_info_yaml['CLEANUP_FLAG']
        Nsec_time_stamp     = submit_info_yaml['TIME_STAMP_NSEC'] 

        Nsec_now  = seconds_since_midnight # current time since midnight

        row_list_merge  = MERGE_INFO_CONTENTS[TABLE_MERGE]
        row_list_split  = MERGE_INFO_CONTENTS[TABLE_SPLIT]
        genversion      = row_list_merge[iver_all][COLNUM_SIM_MERGE_GENVERSION]
        iver            = row_list_merge[iver_all][COLNUM_SIM_MERGE_IVER]

        nlc_gen   = row_list_merge[iver_all][COLNUM_SIM_MERGE_NLC_GEN]
        nlc_write = row_list_merge[iver_all][COLNUM_SIM_MERGE_NLC_WRITE]
        nspec_write = row_list_merge[iver_all][COLNUM_SIM_MERGE_NSPEC_WRITE]
        cpu    = row_list_merge[iver_all][COLNUM_SIM_MERGE_CPU]

        IS_REPEAT =  'REPEAT' in ranseed_key
        IS_CHANGE =  'CHANGE' in ranseed_key

        # - - - - - -
        f.write(f"DOCUMENTATION:\n")
        f.write(f"  GENVERSION: {genversion} \n")  
        f.write(f"\n#                   NLC_GEN   NLC_WRITE " \
                f" NSPEC_WRITE  CPU(minutes)\n")
        f.write(f"  STAT_SUMMARY: \n")  
        f.write(f"  - {'TOTAL':<12}   {nlc_gen:8}   {nlc_write:8} " \
                f" {nspec_write:8}        {cpu}\n")

        # write out same info for each model ... only for RANSEED_REPEAT
        if IS_REPEAT :
            g0 = len("TMP_" + USER4) + len(genversion) + 2
            for row in row_list_split :
                if iver == row[COLNUM_SIM_MERGE_IVER] :
                    TMP_GENV    = row[COLNUM_SIM_MERGE_GENVERSION]
                    g1          = len(TMP_GENV)
                    genv   = f"{TMP_GENV[g0:g1]}" # e.g. SNIaMODEL0
                    nlc_gen   = row[COLNUM_SIM_MERGE_NLC_GEN]
                    nlc_write = row[COLNUM_SIM_MERGE_NLC_WRITE]
                    nspec_write = row[COLNUM_SIM_MERGE_NSPEC_WRITE]
                    cpu       = row[COLNUM_SIM_MERGE_CPU]
                    f.write(f"  - {genv:<12}   " \
                            f"{nlc_gen:8}   {nlc_write:8}  {nspec_write:8} " \
                            f"       {cpu}\n")
        # - - - -
        f.write(f"\n")
        f.write(f"  INPUT_FILES:\n")
        f.write(f"  - {'MASTER':<8}  {input_file} \n")

        # print list of model-input files for this version
        n_file = len(infile_list2d[iver])
        for ifile in range(0,n_file) :
            infile = infile_list2d[iver][ifile]
            model  = model_list2d[iver][ifile] 
            f.write(f"  - {model:<8}  {infile} \n")

        f.write("\n")
        f.write(f"  SUBMIT_DIR:  {CWD}\n")
        f.write("DOCUMENTATION_END:\n")

        # end merge_write_readme

    def merge_cleanup_final(self) :

        # Called after all tasks have complete (see -M arg);
        # tar & gzip most of the contents of SIMLOGS;
        # leave the following files outside tar file so that 
        # they are always visible:
        #    MERGE.LOG  SUBMIT.INFO  ALL.DONE
        # Everything that gets tarred is also removed; therefore
        # specify each item in tar_list and be careful wild cards.
        #
        # Apr 7 2022:
        #  + leave SUBMIT.INFO out of SIMLOGS.tar
        #  + construct several BACKUP*.tar.gz files before SIMLOGS.tar

        submit_info_yaml   = self.config_prep['submit_info_yaml']
        ngen_unit          = submit_info_yaml['NGEN_UNIT']
        simlog_dir         = submit_info_yaml['SIMLOG_DIR']
        
        msg = "\n SIM Clean up SIMLOGS (tar+gzip)"
        logging.info(msg)

        # - - - -
        input_list = f"{SIMGEN_INPUT_LISTFILE} "
        # read list of sim-input files fom list file
        list_file = f"{simlog_dir}/{SIMGEN_INPUT_LISTFILE}"
        with open(list_file,"r") as f:
            for infile in f.read().split():
                infile.strip("\n")
                input_list += f"{infile} "
        
        # - - - -
        util.compress_files(+1, simlog_dir, "TMP_*.LOG",  "TMP-LOG",   "" )
        util.compress_files(+1, simlog_dir, "TMP_*.DONE", "TMP-DONE",  "" )
        util.compress_files(+1, simlog_dir, "TMP_*.YAML", "TMP-YAML",  "" )     
        util.compress_files(+1, simlog_dir, "CPU*",       "CPU",   "" ) 
        util.compress_files(+1, simlog_dir, input_list,   "INPUTS","" ) 
        if ngen_unit > 0 :
            util.compress_files(+1, simlog_dir, "SIMnorm*",  "SIMnorm",  "" ) 

        # combine all of the BACKUP*.tar.gz files into one tar file
        tar_file   = "SIMLOGS.tar"
        tar_list   = "BACKUP_*.tar.gz"
        cd_log     = f"cd {simlog_dir}"
        cmd_tar    = f"tar -cf {tar_file} {tar_list}"
        cmd_rm     = f"rm -rf {tar_list}"
        CMD        = f"{cd_log}; {cmd_tar}; {cmd_rm} "
        os.system(CMD)

        return

    # end merge_cleanup_final

    def merge_cleanup_final_legacy(self) :

        # Called after all tasks have complete (see -M arg);
        # tar & gzip most of the contents of SIMLOGS;
        # leave the following files outside tar file so that 
        # they are always visible:
        #    MERGE.LOG  SUBMIT.INFO  ALL.DONE
        # Everything that gets tarred is also removed; therefore
        # specify each item in tar_list and be careful wild cards.

        submit_info_yaml   = self.config_prep['submit_info_yaml']
        ngen_unit          = submit_info_yaml['NGEN_UNIT']
        simlog_dir         = submit_info_yaml['SIMLOG_DIR']
        
        msg = "\n SIM Clean up SIMLOGS (tar+gzip)"
        logging.info(msg)

        tar_list  = ""
        tar_list += "TMP_* "
        tar_list += "CPU* "
        if ngen_unit > 0 : tar_list += "SIMnorm* "

        tar_list += f"{SIMGEN_INPUT_LISTFILE} "
        # xxx mark tar_list += f"{SUBMIT_INFO_FILE} "

        if KEEP_EVERY_MERGELOG :
            tar_list += f"{MERGE_LOG_FILE}_* "

        # read list of sim-input files fom list file
        list_file = f"{simlog_dir}/{SIMGEN_INPUT_LISTFILE}"
        with open(list_file,"r") as f:
            for infile in f.read().split():
                infile.strip("\n")
                tar_list += f"{infile} "
     
        tar_file   = "SIMLOGS.tar"
        cd_log     = f"cd {simlog_dir}"
        cmd_tar    = f"tar -cf {tar_file} {tar_list}"
        cmd_gzip   = f"gzip {tar_file}"
        cmd_rm     = f"rm -rf {tar_list} {tar_file}"
        CMD        = f"{cd_log}; {cmd_tar}; {cmd_gzip}; {cmd_rm} "
        os.system(CMD)

        return

    # end merge_cleanup_final_legacy

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file.
        # Each info line must be of the form
        #  KEYNAME:  VALUE

        output_dir      = self.config_prep['output_dir']
        survey,idsurvey = util.get_survey_info(output_dir)

        info_lines   = []
        info_lines.append(f"SURVEY:         {survey}")
        info_lines.append(f"IDSURVEY:       {idsurvey}")

        # report sum of NGEN and NWRITE (Dec 2021)
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        MERGE_INFO_CONTENTS,comment_lines = \
                util.read_merge_file(MERGE_LOG_PATHFILE)
        NLC_GEN_SUM = 0
        NLC_WRITE_SUM = 0
        NSPEC_WRITE_SUM = 0
        for row in MERGE_INFO_CONTENTS[TABLE_MERGE] :
            NLC_GEN_SUM     += row[COLNUM_SIM_MERGE_NLC_GEN]
            NLC_WRITE_SUM   += row[COLNUM_SIM_MERGE_NLC_WRITE]
            NSPEC_WRITE_SUM += row[COLNUM_SIM_MERGE_NSPEC_WRITE]
        
        info_lines.append(f"NLC_GEN_SUM:      {NLC_GEN_SUM}")
        info_lines.append(f"NLC_WRITE_SUM:    {NLC_WRITE_SUM}")
        info_lines.append(f"NSPEC_WRITE_SUM:  {NSPEC_WRITE_SUM}")

        return info_lines
        # end get_misc_merge_info

    def get_merge_COLNUM_CPU(self):
        return COLNUM_SIM_MERGE_CPU

# ======= END ======
