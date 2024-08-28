# Created July 2020 by R.Kessler & S. Hinton
#
# Upgrades w.r.t. split_and_fit.pl:
#
#  + NEVT summed from YAML output -> more robust
#  + each merged table examined to ensure NEVT in merged table.
#  + error results in FAIL_SUMMARY.LOG and FAIL_REPEAT scripts
#  + error results in killing all other jobs
#  + options to force failure ... to check how Pippin reacts
#  + CPU diagnostics added to MERGE.LOG
#
# Jan 8 2021:  
#  add OPT_SNCID_LIST option to use events from FITOPT000. To ensure
#  that all of the FITOPT000 are processed first, the iver,fitopt loop
#  was switched to fitopt,iver.
#
# Mar 17 2021: replace sntable_dump.pl with sntable_dump.py (perl->python)
# Apr 07 2021: better error message when appending TEXT table fails.
# Apr 17 2021: write PRIVATE_DATA_PATH to SUBMIT.INFO 
# Apr 25 2021: pass command to nevt_table_check and print command on failure.
# May 13 2021: fix bug setting NOREJECT and OPT_SNCID_LIST option
#
# May 19 2021: if OPT_SNCID_LIST>0 then
#     + set n_job_split = n_core to quickly process FITOPT000
#     + write NEVT_COMMON  to MERGE.LOG
#
# Apr 04 2022: wait for merge_root.exe/merge_hbook.exe to exist in case
#              SNANA build is in progress.
#
# Aug 31 2022 RK - set abort trap for VARNAMES mis-match among split jobs.
#
# Oct 15 2022 RK
#   +  set varnames_ref to be first non-EMPTY varnames in
#          case first varnames is empty.
#   + replace remove_locf_messages.pl with remove_locf_messages.py
#
# Nov 15 2022 RK  implement NMAX_STATE_CHANGE = 10 to help distribute
#                 merge tasks among more cores.
#
# Jul 15 2023 RK - add print_elapse_time calls in merge process to help
#                   identify slow post-process steps.
#
# Oct 12 2023 RK - begin refactor/update to allow stand-alone BayeSN-python LCFIT code.
#                  See dependence on LCFIT_SUBCLASS.
#
# Feb 14 2024 RK 
#    + rename NEVT_SNANA_CUTS to NEVT_LC_CUTS to have a more generic
#                  name for SALT3 and BayeSN.
#    + minor tweaks to post-process for subclass LCFIT_BAYESN.
#    + auto-set kill_on_fail flag if sync-event flag is set to avoid infinite 
#      wait-for-file if FITOPT000.FITRES is never created.
#
# Aug 28 2024 RK
#    + allow new kind of FITOPT to give JOBNAME <jobname> to easily run tests
#      over many old code versions. See -H LCFIT for help.
# - - - - - - - - - -

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time, subprocess
import f90nml
import submit_util as util
import pandas as pd

from   submit_params import *
from   submit_prog_base import Program

# =======================================

LCFIT_SNANA  = 'SNANA'    # default uses SNANA's snlc_fit.exe
LCFIT_BAYESN = 'BAYESN'   # begin integratoin, Oct 2023
LCFIT_SALT3  = 'SALT3'    # placeholder for future ??
LCFIT_SUBCLASS  = None    # see prep_subclass to determine this from progam name

# - - - - - 
# define output table-format info
ITABLE_TEXT  = 0
ITABLE_HBOOK = 1
ITABLE_ROOT  = 2
TABLE_NAME_LIST    = [SUFFIX_FITRES, 'SNANA', 'LCPLOT'] # for TEXT format only

FORMAT_TEXT  = 'TEXT'
FORMAT_HBOOK = 'HBOOK'
FORMAT_ROOT  = 'ROOT'


# these globals will be updated in prep_subclass
TABLE_INPKEY_LIST  = []   # input key per format
TABLE_FORMAT_LIST  = []   # format is used as suffix for file names
NTABLE_FORMAT      = 0 

KEY_SUBCLASS_DICT = {
    'version_photometry' :  {  # key for folder name for real or sim data 
        LCFIT_SNANA  : 'VERSION_PHOTOMETRY',
        LCFIT_BAYESN : '--version_photometry',
        LCFIT_SALT3  : None
    },

    'jobsplit' :  {   # used to internal select independent subsets
        LCFIT_SNANA  : 'JOBSPLIT',
        LCFIT_BAYESN : '--jobsplit',
        LCFIT_SALT3  : None
    },    
        
    'sim_prescale' : {   # optional prescale based on --fast or --faster command line input
        LCFIT_SNANA  : 'SIM_PRESCALE',
        LCFIT_BAYESN : '--sim_prescale',
        LCFIT_SALT3  : None
    },

    # outfile specifications
    'outfile_prefix' :  {
        LCFIT_SNANA  : 'TEXTFILE_PREFIX',
        LCFIT_BAYESN : '--outfile_prefix',
        LCFIT_SALT3  : None
    },
    'outfile_hbook' :  {
        LCFIT_SNANA  : 'HFILE_OUT',
        LCFIT_BAYESN : None,
        LCFIT_SALT3  : None
    },
    'outfile_root' :  {
        LCFIT_SNANA  : 'ROOTFILE_OUT',
        LCFIT_BAYESN : None,
        LCFIT_SALT3  : None
    },

    'dummy_nocomma' : None
}

# - - - - - - - 
# here are a few things for SNANA &SNLCINP namelist,
# but NOT for other sub-classes.
# list of supplemental file inputs to copy IF no path is specified.
COPY_SNLCINP_FILES = \
    [ "KCOR_FILE", "FLUXERRMODEL_FILE", "HEADER_OVERRIDE_FILE", \
      "MAGCOR_FILE", "SIM_MAGCOR_FILE", "SNCID_LIST_FILE",  \
      "NONLINEARITY_FILE", "FUDGE_HOSTNOISE_FILE", "USERTAGS_FILE" ]

# abort if any of these &SNLCINP inputs are specified 
ABORT_ON_SNLCINP_INPUTS =  \
  [ "OUT_EPOCH_IGNORE_FILE",  "SNMJD_LIST_FILE", "SNMJD_OUT_FILE", \
    "SIMLIB_OUT", "SIMLIB_OUTFILE", "MARZFILE_OUT" ]

# - - - - - - -
# define columns in merge file
COLNUM_FIT_MERGE_STATE           = 0  # STATE required in col=0
COLNUM_FIT_MERGE_VERSION         = 1
COLNUM_FIT_MERGE_FITOPT          = 2
COLNUM_FIT_MERGE_NEVT_ALL        = 3  # NEVT0
COLNUM_FIT_MERGE_NEVT_LC_CUTS    = 4  # NEVT1
COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS = 5  # NEVT2
COLNUM_FIT_MERGE_CPU             = 6

NMAX_STATE_CHANGE = 10  # max number of state changes per merge task (Nov 15 2022)

# hard-code option to validate each version; later pass as CONFIG arg
# 1->weak check on README file, 2->strong check, run snana.exe job
OPT_VALIDATE_VERSION = 1

# all merged files names have MERGE prefix before moving them
# to VERSION subdir and removing MERGE prefix.
PREFIX_MERGE = "MERGE"  

# CONFIG key options to append output FITRES file
KEY_APPEND_TABLE_VARLIST  = "APPEND_TABLE_VARLIST"
KEY_APPEND_TABLE_TEXTFILE = "APPEND_TABLE_TEXTFILE"

# prefix for short snana jobs to validate each version
PREFIX_TEMP_SNANA = "TEMP_SNANA"

# define script to dump number of events in table (for hbook & root)
SCRIPT_SNTABLE_DUMP    = "sntable_dump.py"   # Mar 17 2021
Cprogram_SNTABLE_DUMP  = "sntable_dump.exe"  

# define program to merge text-fitres files
PROGRAM_COMBINE_FITRES = "combine_fitres.exe"
NULLVAL_COMBINE_FITRES = -9191   # value for missing CID in extern file

# flags for debug utility to force table-merge failure
FLAG_FORCE_MERGE_TABLE_MISSING = 1
FLAG_FORCE_MERGE_TABLE_CORRUPT = 2

FITOPT_STRING = "FITOPT"

# optional part of FITOPT label to exclude from BBC reject list
# and to exclude from using CID list from FITOPT000
FITOPT_STRING_NOREJECT = "NOREJECT" 

# key to use FITOPT000 sample for all other FITOPTs; for CONFIG and snlc_fit
#                        pippin/submit key       snlc_fit key name
KEY_OPT_SNCID_LIST  = [ 'FLAG_USE_SAME_EVENTS', 'OPT_SNCID_LIST' ]

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================

class LightCurveFit(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_LCFIT
        super().__init__(config_yaml, config_prep)

    def set_output_dir_name(self):
        CONFIG     = self.config_yaml['CONFIG']
        input_file = self.config_yaml['args'].input_file  # for msgerr
        msgerr     = []
        if 'OUTDIR' in CONFIG :
            output_dir_name = os.path.expandvars(CONFIG['OUTDIR'])
        else:
            msgerr.append(f"OUTDIR key missing in yaml-CONFIG")
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp

        return output_dir_name,SUBDIR_SCRIPTS_LCFIT
        # end set_output_dir_name

    def submit_prepare_driver(self):
        
        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file 

        # check which subclass/code-type (SNANA, BayeSN, SALT3
        self.fit_prep_subclass()

        # read/store &SNLCINP namelist (for SNANA subclass only
        if LCFIT_SUBCLASS == LCFIT_SNANA:
            nml = f90nml.read(input_file)
            self.config_prep['snlcinp'] = nml['snlcinp']
            self.fit_prep_check_SNLCINP()
        else:
            self.config_prep['snlcinp'] = None


        # assemble list all possible paths where data/sim could be
        self.fit_prep_path_list()

        # get list of versions to fit, along with path to each version
        self.fit_prep_VERSION()

        # get list of FITOPT-options
        self.fit_prep_FITOPT()

        # compute and store indexing for split jobs
        self.fit_prep_index_lists()

        # check output table options
        self.fit_prep_table_options()

        # copy supplemental input files that don't have a path
        self.fit_prep_copy_files()

        # check OPT_SNCID_LIST option to use same SNe as in FITOPT000
        self.fit_prep_same_sncid()

        # end submit_prepare_driver

        return
    # end submit_prepare_driver

    def fit_prep_subclass(self):
        # Created Oct 2023
        # Default LCFIT_SUBCLASS = LCFIT_SNANA;
        # here check program name for bayesn and salt3 subclass.

        program_name = self.config_prep['program'].lower()
        global LCFIT_SUBCLASS, TABLE_FORMAT_LIST, NTABLE_FORMAT, TABLE_INPKEY_LIST
        
        if 'bayesn' in program_name :
            LCFIT_SUBCLASS = LCFIT_BAYESN  # begin integration Oct 2023
            TABLE_FORMAT_LIST = [ FORMAT_TEXT ] # remove HBOOK and ROOT options
            
        elif 'salt3' in program_name :
            LCFIT_SUBCLASS = LCFIT_SALT3  # placeholder for future
            TABLE_FORMAT_LIST = [ FORMAT_TEXT ] 
            
        else :
            # default is SNANA lcfit
            LCFIT_SUBCLASS = LCFIT_SNANA
            TABLE_FORMAT_LIST = [ FORMAT_TEXT,  FORMAT_HBOOK,  FORMAT_ROOT ] 
        
        NTABLE_FORMAT     = len(TABLE_FORMAT_LIST)

        # load input key(s) TABLE_INPKEY_LIST for each format
        key_dict_list = [ 'outfile_prefix', 'outfile_hbook', 'outfile_root' ] # internal only
        TABLE_INPKEY_LIST = []
        for itable in range(0,NTABLE_FORMAT):
            key_dict = key_dict_list[itable]  # internal key for dictionary
            key_input_file = KEY_SUBCLASS_DICT[key_dict][LCFIT_SUBCLASS] # key in fit-input file
            TABLE_INPKEY_LIST.append(key_input_file) # key per defined out-table format

        msg = f"\n\t !!!! LCFIT_SUBCLASS = {LCFIT_SUBCLASS} !!!! \n"
        logging.info(msg)
        return
        # end fit_prep_subclass
        # - - - - - - - - - 

    def fit_prep_same_sncid(self):

        # Created Jan 8 2021
        # if user sets OPT_SNCID_LIST in CONFIG, then use FITOPT000.FITRES
        # as event mask for all successive FITOPTs so that all 
        # re-analyses use the same events. However, if NOREJECT
        # is part of FITOPT label, ignore OPT_SNCID_LIST.

        CONFIG   = self.config_yaml['CONFIG']
        opt_sncid_list = self.config_prep['opt_sncid_list']  

        KEY_OPT  = KEY_OPT_SNCID_LIST[1]   # key for snlc_fit
        KEY_FILE = 'SNCID_LIST_FILE'       # key for snlc_fit
        argdict_same_sncid = {}

        if opt_sncid_list > 0 :         
            # make list of reference FITOPT000.FITRES file for each 
            # data version
            print(f"\n  PREPARE {KEY_OPT} to USE SAME EVENTS " \
                  f"FOR ALL FITOPTs\n")
            arg_opt       = f"{KEY_OPT} {opt_sncid_list}"
            arg_file_list = []
            output_dir    = self.config_prep['output_dir']
            version_list  = self.config_prep['version_list']   
            for version in version_list :
                #sncid_file = f"{output_dir}/{version}/FITOPT000.FITRES"
                sncid_file = f"../{version}/FITOPT000.FITRES"
                arg_file_list.append(f"{KEY_FILE} {sncid_file}")

            argdict_same_sncid['arg_opt']        = arg_opt
            argdict_same_sncid['arg_file_list']  = arg_file_list

        #print(f" xxx opt_sncid_list = {opt_sncid_list} ")
        #print(f" xxx argdict_same_sncid = {argdict_same_sncid} ")

        # load the goodies
        self.config_prep['argdict_same_sncid'] = argdict_same_sncid

        # end fit_prep_same_scnid


    def fit_prep_copy_files(self):

        # if supplemental input file (passed from &SNLCINP, &FIT, ... ) 
        # has no path (i., not '/'),
        #  + make sure that it exists (abort if not)
        #  + copy to SPLIT_JOBS_LCFIT where split-jobs run
        #
        # Ignore null strings so that things like
        #    MY_WHATEVER_FILE = ''
        # are ignored and doesn't abort.
        #
        # Full paths should be included for all input files defined
        # inside &SNLCINP, in which case this function does nothing.
        # For testing, however, it may be convenient to work with
        # local files; in this case, such files are copied to 
        # SPLIT_JOBS_LCFIT where the split-jobs are run.
        #
        # Oct 13 2023: update to select copy files based on subclass

        script_dir  = self.config_prep['script_dir']
        input_file  = self.config_yaml['args'].input_file 

        # always copy primary input file
        shutil.copy(input_file,script_dir)

        # fetch list of files to copy based on sub-class

        if LCFIT_SUBCLASS == LCFIT_SNANA:
            copy_list = self.get_copy_file_list_SNANA()

        elif LCFIT_SUBCLASS == LCFIT_BAYESN:
            copy_list = self.get_copy_file_list_BAYESN()

        elif LCFIT_SUBCLASS == LCFIT_SALT3:
            copy_list = self.get_copy_file_list_SALT3()


        if len(copy_list) > 0 :
            copy_list_string = " ".join(copy_list)
            msg = f" Copy these input files to {SUBDIR_SCRIPTS_LCFIT}:\n" \
                  f"   {copy_list_string} \n"
            logging.info(msg)
            os.system(f"cp {copy_list_string} {script_dir}/")

        # fit_prep_copy_files

    def get_copy_file_list_SNANA(self):

        # return list of files to copy into script_dir.
        # These are files with no path; files with a full path
        # are not copied (because they can be large)

        snlcinp     = self.config_prep['snlcinp']
        msgerr      = []
        copy_list   = []
        for key_infile in COPY_SNLCINP_FILES :
            if key_infile in snlcinp :
                infile     = snlcinp[key_infile]
                no_path    = "/" not in infile
                if no_path and len(infile) > 0 : 
                    if not os.path.isfile(infile):
                        msgerr.append(f" Missing input file for "\
                                      f"{key_infile} = '{infile}'")
                        msgerr.append(f"Check &SNLCINP in {input_file}")
                        self.log_assert(False,msgerr)
                    copy_list.append(infile)
        return copy_list

        # end get_copy_file_list_SNANA

    def get_copy_file_list_BAYESN(self):
        return []
        # end get_copy_file_list_BAYESN

    def get_copy_file_list_SALT3(self):
        return []
        # end get_copy_file_list_SALT3

        
    def fit_prep_check_SNLCINP(self):

        # Abort on particulare SNLCINP inputs that are not allowed
        # to run in batch mode.

        snlcinp    = self.config_prep['snlcinp']
        key_list_invalid = ""
        for key in ABORT_ON_SNLCINP_INPUTS :
            if key in snlcinp :
                key_list_invalid += f"{key} "
        
        if len(key_list_invalid) > 0 :
            msgerr = []
            msgerr.append(f"Found &SNLCINP keys not allowed in batch mode:")
            msgerr.append(f"  {key_list_invalid} ")
            self.log_assert(False,msgerr)

        # end fit_prep_check_SNLCINP

    def fit_prep_path_list(self):

        # find all possible paths where data or sim could be
        # Store array path_check_list
        #
        # Oct 12 2023: 
        #   refactor to find private_data_path with aribtrary SUBCLASS.

        path_sim_default    = f"{SNDATA_ROOT}/SIM"
        path_data_default   = f"{SNDATA_ROOT}/lcmerge"
        path_sim_list_file  = f"{path_sim_default}/PATH_SNDATA_SIM.LIST"
        path_check_list = []
        
        path_check_list.append(path_data_default)
        path_check_list.append(path_sim_default)

        # check private_data_path, and beware that input is based on LCFIT_SUBCLASS        
        key      = 'private_data_path'
        private_data_path       = None

        if LCFIT_SUBCLASS == LCFIT_SNANA:
            snlcinp  = self.config_prep['snlcinp']
            if key in snlcinp :
                private_data_path = snlcinp[key] # keep env in path name

        elif LCFIT_SUBCLASS == LCFIT_BAYESN :
            pass  # MG: finish this
        
        elif LCFIT_SUBCLASS == LCFIT_SALT3 :
            pass        

        # add private path to list if not accidentally re-defining default path
        # (i.e., ignore mistake to specify PRIVATE_DATA_PATH = SNDATA_ROOT/lcmerge)
        self.config_prep[key] = private_data_path
        if private_data_path is not None:
            path_expand = os.path.expandvars(private_data_path)
            if path_expand != path_data_default: 
                path_check_list.append(path_expand)

        # - - - - 
        with open(path_sim_list_file,"r") as f :
            for line in f:
                path_expand = os.path.expandvars(line.rstrip("\n"))
                path_check_list.append(path_expand)
                
        #print(f" xxx path_check_list = {path_check_list} ")

        self.config_prep['path_check_list']   = path_check_list
        # end fit_prep_path_list

    def fit_prep_VERSION(self):        

        # read/store each VERSION and path where it is locaated; 
        # if wildcard, go fishing in data areas and sim areas 
        # specified in $SNDATA_ROOT/SIM/PATH_SNDATA_SIM.LIST
        # Finally, create output_dir/[VERSION] for each version ...
        # this is where merged table files end up.

        input_file      = self.config_yaml['args'].input_file 
        path_check_list = self.config_prep['path_check_list']
        output_dir      = self.config_prep['output_dir']
        CONFIG          = self.config_yaml['CONFIG']
        KEYLIST         = [ 'VERSION' ]
        version_list_tmp = util.get_YAML_key_values(CONFIG,KEYLIST)

        # for each version, check all paths ... expand
        # version wildcards and make sure that each version
        # appears once and only once.

        version_list      = []
        path_version_list = []
        msg_list          = []
        version_notFound_list = [] # for err msg only
        msgerr  = []
        nerr_validate = 0
        opt_validate  = OPT_VALIDATE_VERSION

        logging.info("\n Search for VERSION wildcards and paths: ")
        for v_tmp in version_list_tmp :  # may include wild cards
            found = False
            for path in path_check_list :
                v_list  = sorted(glob.glob(f"{path}/{v_tmp}"))
                for v in v_list:
                    found   = True
                    #j_slash = v.rindex('/');    version = v[j_slash+1:]
                    version = os.path.basename(v)
                    # avoid tar files and gz files
                    if '.tar' in v : continue
                    if '.gz'  in v : continue

                    validate,msg_status = \
                        self.fit_validate_VERSION(opt_validate,path,version)
                    msg = f"   Found VERSION {version}    {msg_status}"
                    logging.info(f"{msg}")
                    if validate is False :
                        nerr_validate += 1  # store all errors before abort
                    version_list.append(version)
                    path_version_list.append(path)
                    msg_list.append(msg)  # used below in case of log_assert

            if not found :
                logging.info(f"   ERROR: could not find {v_tmp} ")
                version_notFound_list.append(v_tmp)

        if len(version_notFound_list) > 0 :
            msgerr.append("Could not find the following VERSION(s)")
            for v in version_notFound_list :
                msgerr.append(f"   {v} ")
            self.log_assert(False,msgerr)

        # check for duplicates 
        n_dupl_list,dupl_list = util.find_duplicates(version_list)
        n_list = len(n_dupl_list)
        if n_list > 0 :
            msgerr.append("Invalid multiple VERSIONs")
            for i in range(0,n_list):
                dupl   = dupl_list[i]
                n_dupl = n_dupl_list[i]
                msgerr.append(f"\t Found {n_dupl} {dupl} directories")
            msgerr.append(f"Only one VERSION directory allowed.")
            msgerr.append(f"Check these directories:")
            for path in set(path_version_list):
                msgerr.append(f"\t {path}")
            self.log_assert(False,msgerr)

        # abort on validation errors, otherwise remove TEMP_SNANA files
        # Note that errors above are stored so that here ALL errors
        # are printed in abort message.
        if nerr_validate > 0 :
            msgerr.append(f"Validation error on {nerr_validate} versions; ")
            msgerr += msg_list
            self.log_assert(False,msgerr)
        else :
            cmd = f"cd {output_dir}; rm {PREFIX_TEMP_SNANA}* 2>/dev/null"
            os.system(cmd)
            #print(f" xxxx cmd to remove TEMP_SNANA, \n {cmd} ")

        # finally, create subdir for each version
        for v in version_list :
            v_dir = f"{output_dir}/{v}"
            os.mkdir(v_dir)

        # load lists for version and path_version
        self.config_prep['n_version']         = len(version_list)
        self.config_prep['version_list']      = version_list
        self.config_prep['path_version_list'] = path_version_list

        logging.info("")

        # end fit_prep_VERSION

    def fit_validate_VERSION(self, opt, path, version):
        # Driver to select version-validaton based on input opt.
        validate=True ;   msg= ''
        if opt == 1 :
            validate,msg = self.fit_validate1_VERSION(path,version)
        if opt == 2 :
            validate,msg = self.fit_validate2_VERSION(path,version)

        return validate,msg
        # end fit_validate_VERSION

    def fit_validate1_VERSION(self, path, version):

        # validate based on README file in path.
        # Validation fails if README does not exist,
        # or if FAIL appears on first line of README.

        validate      = True  # default is validation success
        string_status = ""    # no string needed for success
        readme_file = f"{path}/{version}/{version}.README"
        
        # if no README file, it's a no-brainer failure
        if os.path.isfile(readme_file) is False :
            validate = False
            string_status = "Validate ERROR:  missing README"

        # check if first word in readme file is FAIL 
        # left by sim-submission.
        with open(readme_file,"r") as r :
            line = r.readline().rstrip("\n")
            if 'FAIL' in line :
                validate = False
                string_status = "Validate ERROR : 'FAIL' found in README"

        return validate, string_status
        # end fit_validate1_VERSION

    def fit_validate2_VERSION(self, path, version):

        # Strong validation method:
        # run short snana.exe job on input  version and create
        # YAML file; if NEVT_TOT=0, or YAML file is not produced,
        # return False. If NEVT_TOT > 0, return True
        # Also return comment string; SUCCESS for True,
        # and error msg for False.

        snana_dir       = self.config_yaml['args'].snana_dir
        output_dir      = self.config_prep['output_dir']
        cddir           = f"cd {output_dir}"
        key_nevt        = 'NEVT_TOT'
        nevt_proc       = 10  # process this many to validate
        validate        = True 

        textfile_prefix = f"{PREFIX_TEMP_SNANA}_{version}"
        log_file        = f"{textfile_prefix}.LOG"
        yaml_file       = f"{textfile_prefix}.YAML"

        if snana_dir is None :
            cmd_snana = f"snana.exe NOFILE  "
        else:
            cmd_snana = f"{snana_dir}/bin/snana.exe NOFILE  "

        cmd_snana += f"VERSION_PHOTOMETRY {version}  "
        cmd_snana += f"MXEVT_PROCESS {nevt_proc}  "
        cmd_snana += f"JOBSPLIT 1 1  "
        cmd_snana += f"TEXTFILE_PREFIX {textfile_prefix}  " \
                     f"SNTABLE_LIST '' " 

        os.system(f"{cddir} ; {cmd_snana} > {log_file} ")

        # read and parse yaml file
        YAML_FILE = f"{output_dir}/{yaml_file}"
        string_status = "Validate SUCCESS"  # default status

        if os.path.isfile(YAML_FILE) :
            snana_yaml = util.extract_yaml(YAML_FILE, None, None )
            nevt       = snana_yaml[key_nevt]
            if nevt == 0 :
                validate      = False
                string_status = f"Validate ERROR: NEVT=0"
        else:
            validate      = False
            string_status = f"Validate ERROR: cannot analyze with snana"

        return validate, string_status

        # end fit_validate2_VERSION

    def fit_prep_FITOPT(self) :

        # read/store list of FITOPT options, along with 'fitnums'
        # FITOPT000, FITOPT001, etc ... These fitnums are used 
        # to create file names.
        # If there is a lable in /LABEL/, strip it out and store it
        # in SUBMIT_INFO, FITOPT.README along with fitopt and fitnum,
        #
        # Jan 24 2021: refactor using prep_jobopt_list utility

        output_dir        = self.config_prep['output_dir']
        CONFIG            = self.config_yaml['CONFIG']
        args              = self.config_yaml['args']
        ignore_fitopt     = args.ignore_fitopt

        if ignore_fitopt: 
            # user option to ignore FITOPTs
            fitopt_rows = []
        else:
            # default: read FITOPT info
            KEYLIST       = [ FITOPT_STRING ]    # key under CONFIG
            fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)

        # check for OPT_SNCID_LIST ... just store it here for later
        # Check multiple key-options
        opt_sncid_list = 0
        for key in KEY_OPT_SNCID_LIST:
            if key in CONFIG :  
                opt_sncid_list = CONFIG[key]
                if opt_sncid_list>0 and args.nomerge:
                    msgerr = []
                    msgerr.append(f"nomerge option does not work " \
                                  f" with OPT_SNCID_LIST")
                    self.log_assert(False,msgerr)
            
        # - - - - - -
        fitopt_dict = util.prep_jobopt_list(fitopt_rows,FITOPT_STRING,1,None)

        fitopt_arg_list   = fitopt_dict['jobopt_arg_list']
        fitopt_num_list   = fitopt_dict['jobopt_num_list']
        fitopt_label_list = fitopt_dict['jobopt_label_list']
        n_fitopt          = fitopt_dict['n_jobopt']

        # - - - - - - - - - 
        # update list for symbolic links to FITOPT000 [DEFAULT]
        link_FITOPT000_list = []
        for arg,num in zip(fitopt_arg_list,fitopt_num_list) :
            if self.is_sym_link(arg) :
                link_FITOPT000_list.append(num)

        # - - - - - - - -
        logging.info(f"  Found {n_fitopt-1} FITOPT variations.")
        logging.info(f"  link_FITOPT000_list: {link_FITOPT000_list}")

        self.config_prep['n_fitopt']            = n_fitopt
        self.config_prep['fitopt_num_list']     = fitopt_num_list
        self.config_prep['fitopt_arg_list']     = fitopt_arg_list
        self.config_prep['fitopt_label_list']   = fitopt_label_list
        self.config_prep['link_FITOPT000_list'] = link_FITOPT000_list
        self.config_prep['n_fitopt_link']       = len(link_FITOPT000_list)
        self.config_prep['opt_sncid_list']      = opt_sncid_list

        # end fit_prep_FITOPT

    def fit_prep_index_lists(self):

        # Construct 1D sparse lists of                   
        #   version   index iver    0 to n_job_tot 
        #   fitopt    index iopt    0 to n_job_tot 
        #   split-job index isplit  0 to n_job_tot 
        #                                                  
        # These 1D arrays can be used in 1D loops in range(0,n_job_tot)
        # instead of 3D for blocks.
        #
        # Notation warning: 
        #   n_job_tot is number of real jobs that are NOT sym links.
        #   n_fitopt  is total number of FITOPS that includes sym links
        #   n_fitopt_link is number of FITOPTs that are sym links
        #
        n_version        = self.config_prep['n_version']
        n_fitopt_tot     = self.config_prep['n_fitopt']  # all FITOPT
        n_fitopt_link    = self.config_prep['n_fitopt_link'] # links to FITOPT000
        version_list     = self.config_prep['version_list']
        fitopt_arg_list  = self.config_prep['fitopt_arg_list']
        n_core           = self.config_prep['n_core']
        opt_sncid_list   = self.config_prep['opt_sncid_list']
        do_dump = True

        # first figure out how many split jobs
        n_fitopt_tmp = n_fitopt_tot - n_fitopt_link # number of FITOPTs to process
        n_job_tmp    = n_version * n_fitopt_tmp  # N_job if no splitting
        n_job_split  = int(n_core/n_job_tmp)

        # - - - - - - -
        # check special cases to alter n_job_split

        # require at least 1 job
        if n_job_split == 0 : n_job_split = 1

        # if waiting for FITOPT000, distribute over all cores to avoid
        # long wait for FITOPT000. But no more than 100 splits
        if opt_sncid_list > 0 :
            n_job_split = n_core
            if n_job_split > 100: n_job_split = 100

        # note that n_job_tot does NOT include symbolic links
        n_job_tot = n_job_split * n_job_tmp

        logging.info(f"  Determine number of split jobs that are not sym links:")
        logging.info(f"    n_version x n_fitopt = {n_version} x {n_fitopt_tmp} "\
                     f"= {n_job_tmp}   # excludes sym links")
        logging.info(f"    n_job[tot,split] = {n_job_tot},{n_job_split}" \
                     f" for {n_core} cores   # excludes sym links" )
        logging.info("")

        # create sparse index lists that include sym links
        iver_list=[] ;  iopt_list=[];   isplit_list=[]
        size_sparse_list = 0

        # Loop over iopt (fitopt) first to ensure that all FITOPT000
        # are processed first ... matters when using OPT_SNCID_LIST.
        for iopt in range(0,n_fitopt_tot):
            for iver in range(0,n_version):
                for isplit in range(0,n_job_split):
                    iver_list.append(iver)
                    iopt_list.append(iopt)
                    isplit_list.append(isplit)
                    size_sparse_list += 1  # n_job(proc+links)

        self.config_prep['size_sparse_list'] = size_sparse_list
        self.config_prep['n_job_tot']     = n_job_tot  # excluded symLinks
        self.config_prep['n_done_tot']    = size_sparse_list # proc+links
        self.config_prep['n_job_split']   = n_job_split
        self.config_prep['iver_list']     = iver_list
        self.config_prep['iopt_list']     = iopt_list
        self.config_prep['isplit_list']   = isplit_list

        # store number of jobs that are simply symbolic links
        n_job_link = n_version * n_fitopt_link * n_job_split
        self.config_prep['n_job_link']    = n_job_link

        # end fit_prep_index_lists

    def fit_prep_table_options(self):
        
        # check table options in &SNLCINP namelist, and set
        # logical flag for each table format: TEXT, HBOOK, ROOT
        # In &SNLCINP, existance of HFILE_OUT is a logical flag
        # to create an output HBOOK file; ROOTFILE_OUT is a flag
        # to create a root file. TEXT format is always output
        # for back-compatibility.

        msgerr = []
        snlcinp     = self.config_prep['snlcinp']        
        script_dir  = self.config_prep['script_dir']

        self.config_prep['use_table_format'] = [ False ] * NTABLE_FORMAT
        use_table_format = self.config_prep['use_table_format']
        # xxx mark use_table_format[ITABLE_TEXT] = True  # always force this format

        use_table_format = []
        for key, fmt in zip(TABLE_INPKEY_LIST, TABLE_FORMAT_LIST):
            if fmt == FORMAT_TEXT:
                use = True  # always force this to be true
            else:
                use = False

            if LCFIT_SUBCLASS == LCFIT_SNANA :
                if key in snlcinp :  # snana selects among TEXT, ROOT, HBOOK   
                    use = True

            use_table_format.append(use)
            logging.info(f"  Write {fmt:6s} output table : {use}")

        self.config_prep['use_table_format'] = use_table_format

        # get list of tables in SNTABLE_LIST

        if LCFIT_SUBCLASS == LCFIT_SNANA:
            # SNANA tables include SNANA, FITRES, LCPLOT, OUTLIER,  ...
            sntable_string = snlcinp['sntable_list'] 
        else:
            # only one kind of table for non-SNANA fitters?
            sntable_string = 'FITRES'

        sntable_list = sntable_string.split()

        logging.info(f"  SNTABLE_LIST = {sntable_list} ")

        if LCFIT_SUBCLASS == LCFIT_SNANA:
            self.check_options_SNANA_table()

        logging.info("")

        return
        # end fit_prep_table_options

    def check_options_SNANA_table(self):

        # Created Oct 2023
        # sanity checks for SNANA tables

        script_dir       = self.config_prep['script_dir']
        use_table_format = self.config_prep['use_table_format']

        #sys.exit(f"\n xxx use_table_format = {use_table_format}  ITABLE_ROOT={ITABLE_ROOT}\n")

        msgerr = []
        # if appending variables to FITRES file, make sure
        # that either HBOOK or ROOT is specified
        CONFIG    = self.config_yaml['CONFIG']
        key       = KEY_APPEND_TABLE_VARLIST
        if key in CONFIG :
            require = False
            if use_table_format[ITABLE_ROOT]:
                require = True
            if use_table_format[ITABLE_HBOOK]:
                require = True
            if not require:
                msgerr.append(f" {key} found in CONFIG input")
                msgerr.append(f" but could not find {FORMAT_HBOOK} or {FORMAT_ROOT}. ")
                self.log_assert(False,msgerr)
            
        # for LCFIT_SNANA
        # if APPEND_TABLE_TEXTFILE is defined, make sure that it exists
        # and that it has a full path.
        key = KEY_APPEND_TABLE_TEXTFILE
        if key in CONFIG :
            delim = None
            if ',' in CONFIG[key]: delim=','
            for text_file in CONFIG[key].split(delim):
                text_file = os.path.expandvars(text_file)
                msg = f"  Every output FITRES file will be appended with: \n"\
                      f"    {text_file} "
                logging.info(msg)

                if not os.path.isfile(text_file):
                    msgerr.append(f"{text_file} does not exist.")
                    msgerr.append(f"Check {key} argument under CONFIG block.")
                    self.log_assert(False,msgerr)

                if '/' not in text_file :
                    shutil.copy(text_file,script_dir)


        return
        # end check_options_SNANA_table

    def write_command_file(self,icpu,f):
        # For this icpu, write full set of sim commands to
        # already-opened command file with pointer f. 
        # Function returns number of jobs for this cpu 

        # loop over version, fitopt
        size_sparse_list = self.config_prep['size_sparse_list'] 
        n_job_tot        = self.config_prep['n_job_tot']  # does NOT include symlinks
        iver_list        = self.config_prep['iver_list']
        iopt_list        = self.config_prep['iopt_list']
        isplit_list      = self.config_prep['isplit_list']

        n_version        = self.config_prep['n_version']
        n_fitopt         = self.config_prep['n_fitopt']
        version_list     = self.config_prep['version_list']
        fitopt_arg_list  = self.config_prep['fitopt_arg_list']
        fitopt_num_list  = self.config_prep['fitopt_num_list']
        n_core           = self.config_prep['n_core']
        n_job_cpu        = 0

        n_job_local = 0 ;   n_job_real=0 

        for iver,iopt,isplit in zip(iver_list,iopt_list,isplit_list):

            index_dict = {
                'iver':iver, 'iopt':iopt, 'isplit':isplit, 'icpu':icpu
            }  
            n_job_local += 1
            if self.is_sym_link(fitopt_arg_list[iopt]) : continue
            n_job_real += 1  # use this to skip links

            #if ( (n_job_local-1) % n_core ) == icpu :
            if ( (n_job_real-1) % n_core ) == icpu :

                n_job_cpu += 1

                job_info_fit   = self.prep_JOB_INFO_fit(index_dict)
                util.write_job_info(f, job_info_fit, icpu)

                job_info_merge = \
                    self.prep_JOB_INFO_merge(icpu,n_job_real,False) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

        # - - - - 

        if n_job_real != n_job_tot :
            msgerr = []
            msgerr.append(f"Expected {n_job_tot} total jobs;")
            msgerr.append(f"but found {n_job_local} jobs.")
            self.log_assert(False,msgerr)

        return n_job_cpu

        # end write_command_file

    def is_sym_link(self,fitopt_arg):
        # for input fitopt argument, return True if it means
        # symbolic link to FITOPT000.
        if fitopt_arg == f"{FITOPT_STRING}000"    : return True
        if fitopt_arg == 'DEFAULT'                : return True
        return False

    def prep_JOB_INFO_fit(self,index_dict):
        # Return JOB_INFO dictionary with 
        #   cd job_dir
        #   program arg_list  > log_file
        #   touch TMP_[xxx].DONE
        #
        # Inputs
        #   index_dict = dictionary of indices for this job
        #
        # May 13 2021: fix bug setting NOREJECT and OPT_SNCID_LIST option

        # strip off indices from input dictionary
        iver   = index_dict['iver']  # version index for data or sim
        iopt   = index_dict['iopt']  # FITOPT index
        isplit = index_dict['isplit']+1  # fortran like index for file names
        icpu   = index_dict['icpu']      # cpu index

        CONFIG        = self.config_yaml['CONFIG']
        args          = self.config_yaml['args']
        input_file    = args.input_file 
        kill_on_fail  = args.kill_on_fail
        check_abort   = args.check_abort
        program       = self.config_prep['program']
        output_dir    = self.config_prep['output_dir']
        script_dir    = self.config_prep['script_dir']
        version       = self.config_prep['version_list'][iver]
        fitopt_arg    = self.config_prep['fitopt_arg_list'][iopt]
        fitopt_num    = self.config_prep['fitopt_num_list'][iopt]
        fitopt_label  = self.config_prep['fitopt_label_list'][iopt]
        #fitopt_global = CONFIG.setdefault('FITOPT_GLOBAL',None)
        fitopt_global = self.get_fitopt_global(version)

        # Aug 28 2024: check FITOPT that is actually a different program name,
        #              such as an older code version; e.g., 
        #    FITOPT:
        #    - /v11_04e/  JOBNAME  /products/SNANA_v11_04e/bin/snlc_fit.exe
        #
        if 'JOBNAME' in fitopt_arg:
            program = os.path.expandvars(fitopt_arg.split()[1])
            fitopt_arg = ''
            
        use_table_format = self.config_prep['use_table_format']
        n_job_split   = self.config_prep['n_job_split']
        split_num     = f"SPLIT{isplit:03d}"
        prefix        = f"{version}_{fitopt_num}_{split_num}"
        done_file     = f"{prefix}.DONE"
        log_file      = f"{prefix}.LOG"
        yaml_file     = f"{prefix}.YAML"
        arg_list      = []
        JOB_INFO      = {}

        JOB_INFO['job_dir']     = script_dir  # where to run job
        JOB_INFO['program']     = program
        JOB_INFO['input_file']  = input_file
        JOB_INFO['log_file']    = log_file
        JOB_INFO['done_file']   = done_file
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['check_abort']   = check_abort

        # set command line arguments

        key = KEY_SUBCLASS_DICT['version_photometry'][LCFIT_SUBCLASS]
        arg_list.append(f"  {key} {version}")

        key = KEY_SUBCLASS_DICT['jobsplit'][LCFIT_SUBCLASS]
        arg_list.append(f"  {key} {isplit} {n_job_split}")
            
        # check fast option to prescale sims by 10 (data never pre-scaled)
        if args.prescale > 1 :
            key = KEY_SUBCLASS_DICT['sim_prescale'][LCFIT_SUBCLASS]
            arg_list.append(f"  {key} {args.prescale}")
            
        if args.require_docana :
            arg_list.append(f"  REQUIRE_DOCANA 1")

        # tack on outFile for each table format. For TEXT, do NOT
        # include suffix in TEXTFILE_PREFIX argument
        for itab in range(0,NTABLE_FORMAT) :
            if use_table_format[itab] :
                key    = TABLE_INPKEY_LIST[itab] 
                fmt    = TABLE_FORMAT_LIST[itab]
                arg    = f"{key:<16} {prefix}"
                if fmt != TABLE_FORMAT_LIST[ITABLE_TEXT] :
                    arg += f".{fmt}"  # format is suffix
                arg_list.append(f"  {arg}")

        
        # user-define FITOPT options. 
        if fitopt_global is not None:
            arg_list.append(f"  {fitopt_global}")

            
        arg_list.append(f"{fitopt_arg}")

        # Jan 8, 2021: option to use CID list from FITOPT000
        opt_sncid_list = self.config_prep['opt_sncid_list']

        if args.check_abort: 
            arg_list.append("MXEVT_CUTS 1")
            opt_sncid_list = 0  # disable event-sync feature

        if fitopt_label is None:
            NOREJECT = None
        else:
            NOREJECT = FITOPT_STRING_NOREJECT in fitopt_label

        # for sync-event feature, each FITOPT must wait for FITOPT000
        # to get its cid list.
        if iopt > 0 and opt_sncid_list > 0  and NOREJECT is False :
            argdict_same_sncid = self.config_prep['argdict_same_sncid']
            arg_opt   = argdict_same_sncid['arg_opt']             # KEY OPT
            arg_file  = argdict_same_sncid['arg_file_list'][iver] # KEY FILE
            wait_file = arg_file.split()[1]  # just the FILE name
            arg_list.append(f"{arg_file}")
            arg_list.append(f"{arg_opt}")
            JOB_INFO['wait_file']     = wait_file
            JOB_INFO['kill_on_fail']  = True
            args.kill_on_fail         = True   # avoid infinite wait on FITOPT000 failure

        # - - - 
        JOB_INFO['arg_list']   = arg_list

        # - - - - - - - - 
        # On FITOPT000, check to create sym links from other FITOPTs
        # The sym links are created AFTER FITOPT000 finishes to 
        # ensure time-ordered file creation for merge process.
        link_FITOPT000_list = self.config_prep['link_FITOPT000_list']
        n_link              = len(link_FITOPT000_list)
        if iopt == 0 and n_link > 0 :
            sym_link_dict = { 
                'version':version, 'prefix_file':prefix, 
                'fitopt_num':fitopt_num }
            JOB_INFO['sym_link_list'] = self.get_sym_link_list(sym_link_dict)
        # - - - - - - - - 
            
        return JOB_INFO

        # end prep_JOB_INFO_fit

    def get_fitopt_global(self,version):
        # Created Aug 29 2024
        # Return fitopt_global st
        CONFIG        = self.config_yaml['CONFIG']
        key_fitopt_pairs = [(key,value) for key,value in \
                            CONFIG.items() if key.startswith("FITOPT_GLOBAL")]

        if len(key_fitopt_pairs) == 0:
            fitopt_global = None
        else:
            fitopt_global = ''            
            for key, fitopt  in key_fitopt_pairs:       
                version_pattern = util.extract_arg(key)
                if version_pattern == '' or version_pattern in version:
                    fitopt_global += f"{fitopt}  "

        LDMP = False
        if LDMP:
            print(f" xxx ------------------------------------------- ")
            print(f" xxx get_fitopt_global DUMP for version = {version} :")
            print(f" xxx key_fitopt_pairs = \n{key_fitopt_pairs}")        
            print(f" xxx fitopt_global = \n{fitopt_global}")
                
        
        return fitopt_global
        # end get_fitopt_global
        
    def get_sym_link_list(self,sym_link_dict):
        # prepare commands for symbolic links from FITOPTnnn to FITOPT000,
        # where nnn > 0
        # This function is called when user input includes example such as
        # FITOPT:
        #  - bla bla
        #  - bla2 bla2
        #  - FITOPT000  
        # Here the sym link commands are created for DONE, LOG, YAML 
        # outputs, and for FITOP003 -> FITOPT000

        link_FITOPT000_list = self.config_prep['link_FITOPT000_list']
        version             = sym_link_dict['version']
        fitopt_num          = sym_link_dict['fitopt_num']
        prefix_file         = sym_link_dict['prefix_file']
        suffix_link_list    = [ 'LOG ', 'DONE', 'YAML' ]
        sym_link_list   = []  # init return arg list

        for fitopt_link in link_FITOPT000_list : # e.g., FITOPT003
            prefix_link  = prefix_file.replace(fitopt_num,fitopt_link)
            for suffix in suffix_link_list :
                orig_file   = f"{prefix_file}.{suffix}" # with FITOPT000
                link_file   = f"{prefix_link}.{suffix}" # with FITOPTnnn
                sym_link    = f"ln -s {orig_file} {link_file}"
                sym_link_list.append(sym_link)

        return sym_link_list
        # end get_sym_link_list


    def append_info_file(self,f):

        # Create SUBMIT.INFO file for merge task and also for
        # downstream scripts.

        CONFIG            = self.config_yaml['CONFIG']
        n_job_link        = self.config_prep['n_job_link']
        output_dir        = self.config_prep['output_dir']
        n_fitopt          = self.config_prep['n_fitopt']
        version_list      = self.config_prep['version_list']
        fitopt_arg_list   = self.config_prep['fitopt_arg_list']
        fitopt_num_list   = self.config_prep['fitopt_num_list']
        fitopt_label_list = self.config_prep['fitopt_label_list']
        link_FITOPT000_list = self.config_prep['link_FITOPT000_list']
        use_table_format  = self.config_prep['use_table_format']
        ignore_fitopt     = self.config_yaml['args'].ignore_fitopt
        private_data_path = self.config_prep['private_data_path']
        opt_sncid_list    = self.config_prep['opt_sncid_list']

        f.write(f"\n# Fit info\n")
        f.write(f"N_JOB_LINK:          {n_job_link}   " \
                f"# Njob with link to FITOPT000\n")
        f.write(f"JOBFILE_WILDCARD:    '*SPLIT*' \n")
        f.write(f"TABLE_FORMAT_LIST:   {TABLE_FORMAT_LIST} \n")
        f.write(f"USE_TABLE_FORMAT:    {use_table_format} \n")
        f.write(f"IGNORE_FITOPT:       {ignore_fitopt}\n")
        f.write(f"PRIVATE_DATA_PATH:   {private_data_path} \n")
        f.write(f"{KEY_OPT_SNCID_LIST[0]}:   {opt_sncid_list}   " \
                "# >0 -> FITOPT>0 uses events from FITOPT000\n")

        key_misc_list = [ KEY_APPEND_TABLE_VARLIST, KEY_APPEND_TABLE_TEXTFILE]
        for key in key_misc_list :
            key_yaml = key + ':'
            if key in CONFIG :
                f.write(f"{key_yaml:<22}  {CONFIG[key]} \n")

        f.write("\n")
        f.write("VERSION_LIST: \n")
        for v in version_list :
            f.write(f"  - {v}\n")

        f.write("\n")
        f.write(f"FITOPT_LIST: \n")
        row = [ 0 ] * 3
        for i in range(0,n_fitopt):
            num   = fitopt_num_list[i]
            label = fitopt_label_list[i] 
            arg   = fitopt_arg_list[i]
            row[COLNUM_FITOPT_NUM]   =  num
            row[COLNUM_FITOPT_LABEL] =  label
            row[COLNUM_FITOPT_ARG]   =  arg
            #row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")

        f.write(f"LINK_FITOPT000_LIST: {link_FITOPT000_list}\n")
        f.write("\n")

        # end  append_info_file

    def create_merge_table(self,f):

        # Required element of submit process. Before submitting jobs,
        # create initial merge file with all WAIT states.
        # This file is read and updated frequently by merge
        # process invoked by -m argument to submit_batch_jobs.sh
        # A locally defined MERGE_INFO structure is passed to 
        # a generic write_MERGE_INFO function to create MERGE.LOG/
        # Uses YAML format, and for human-readability there is a 
        # one line commented header before each YAML table.

        n_version       = self.config_prep['n_version']        
        version_list    = self.config_prep['version_list']
        output_dir      = self.config_prep['output_dir']
        fitopt_num_list = self.config_prep['fitopt_num_list']
        n_job_split     = self.config_prep['n_job_split']

        # create only MERGE table ... no need for SPLIT table
        header_line_merge = \
            " STATE   VERSION  FITOPT  " \
            "NEVT_ALL  NEVT_LC_CUTS NEVT_FIT_CUTS  CPU"

        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for version in version_list :
            for num in fitopt_num_list :
                # define ROW here is fragile in case columns are changed
                ROW_MERGE = [ STATE, version, num, 0, 0, 0, 0.0 ]
                INFO_MERGE['row_list'].append(ROW_MERGE)  
        util.write_merge_file(f, INFO_MERGE, [] ) 

        # end create_merge_file

    def merge_config_prep(self,output_dir):

        self.fit_prep_subclass()  # Feb 15 2024

        # fit-specific settings to config_prep that are needed later.
        submit_info_yaml = self.config_prep['submit_info_yaml'] 

        # restore list of possible table formats for this subclass (Oct 2023)
        global TABLE_FORMAT_LIST, NTABLE_FORMAT
        TABLE_FORMAT_LIST = submit_info_yaml['TABLE_FORMAT_LIST']
        NTABLE_FORMAT     = len(TABLE_FORMAT_LIST)

        ## self.config_prep['output_dir']   = output_dir 

        # end merge_config_prep

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.
        # To continuously update statistics during RUN state,
        # n_state_change is always set > 0 if in RUN state,
        # even if there is not STATE change.

        MERGE_LAST          = self.config_yaml['args'].MERGE_LAST
        submit_info_yaml    = self.config_prep['submit_info_yaml']
        snana_version       = self.config_prep['snana_version']
        script_dir          = submit_info_yaml['SCRIPT_DIR']
        n_job_split         = submit_info_yaml['N_JOB_SPLIT']
        link_FITOPT000_list = submit_info_yaml['LINK_FITOPT000_LIST']
        COLNUM_STATE   = COLNUM_FIT_MERGE_STATE 
        COLNUM_VERS    = COLNUM_FIT_MERGE_VERSION 
        COLNUM_FITOPT  = COLNUM_FIT_MERGE_FITOPT
        COLNUM_NEVT0   = COLNUM_FIT_MERGE_NEVT_ALL 
        COLNUM_NEVT1   = COLNUM_FIT_MERGE_NEVT_LC_CUTS
        COLNUM_NEVT2   = COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS
        COLNUM_CPU     = COLNUM_FIT_MERGE_CPU
        
        key_tot, key_tot_sum, key_list = \
                self.keynames_for_job_stats('NEVT_TOT')

        # use old key NEVT_SNANA_CUTS if snana_version < v11_05k
        # Beware that this fails for reading snlc_fit.exe from restored git clone
        KEY_YAML = 'NEVT_LC_CUTS'
        if snana_version < 'v11_05k' : KEY_YAML = 'NEVT_SNANA_CUTS'  # .xyz
        key_snana, key_snana_sum, key_snana_list = \
                self.keynames_for_job_stats(KEY_YAML)
        
        key_lcfit, key_lcfit_sum, key_lcfit_list = \
                self.keynames_for_job_stats('NEVT_LCFIT_CUTS')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')
        key_list  = [ key_tot, key_snana, key_lcfit, key_cpu ]

        row_list_merge   = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []

        if MERGE_LAST:
            nmax_state_change = 9999999
        else:
            nmax_state_change = NMAX_STATE_CHANGE

        #  - - - - -
        irow = 0
        version_last  = "bla"
        Finished_last = False
 
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input

            # strip off row info
            STATE       = row[COLNUM_STATE]
            version     = row[COLNUM_VERS]
            fitopt_num  = row[COLNUM_FITOPT]
            search_wildcard = f"{version}_{fitopt_num}_SPLIT*"

            # check if DONE or FAIL ; i.e., if Finished
            Finished = (STATE == SUBMIT_STATE_DONE) or \
                       (STATE == SUBMIT_STATE_FAIL)

            do_check_state = (not Finished  and n_state_change < nmax_state_change)

            # avoid checking FITOPT files that have no chance when last
            # FITOPT didn't finish. Make sure to check again when version updates.
            if version == version_last and not Finished_last and not MERGE_LAST:
                do_check_state = False

            # - - - - - - - - - - - -  -
            if do_check_state :
                NEW_STATE = STATE

                # get list of LOG, DONE, and YAML files 
                log_list, done_list, yaml_list = \
                    util.get_file_lists_wildcard(script_dir,search_wildcard)

                # careful to sum only the files that are NOT None
                NLOG   = sum(x is not None for x in log_list)  
                NDONE  = sum(x is not None for x in done_list)  
                NYAML  = sum(x is not None for x in yaml_list)  

                if NLOG > 0:
                    NEW_STATE = SUBMIT_STATE_RUN
                if NDONE == n_job_split :
                    NEW_STATE = SUBMIT_STATE_DONE
                    
                    job_stats = self.get_job_stats(script_dir, 
                                                   log_list, yaml_list, key_list)
                    # check for failures in snlc_fit jobs.
                    nfail = job_stats['nfail']
                    if nfail > 0 :  NEW_STATE = SUBMIT_STATE_FAIL
                
                    row[COLNUM_STATE]  = NEW_STATE
                    row[COLNUM_NEVT0]  = job_stats[key_tot_sum]
                    row[COLNUM_NEVT1]  = job_stats[key_snana_sum]
                    row[COLNUM_NEVT2]  = job_stats[key_lcfit_sum]
                    row[COLNUM_CPU]    = job_stats[key_cpu_sum]

                    if fitopt_num in link_FITOPT000_list :
                        row[COLNUM_CPU] = 0.0 # zero CPU for sym links

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             # assume nevt changes

            irow += 1
            version_last  = version
            Finished_last = Finished

        # - - - - - - - 
        # first return arg (row_split) is null since there is 
        # no need for a SPLIT table
        row_list_dict = {
            'row_split_list' : [],
            'row_merge_list' : row_list_merge_new,
            'row_extra_list' : []
        }
        return row_list_dict, n_state_change

        # end merge_update_state

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        # irow is the row to wrapup in MERGE_INFO_CONTENTS
        # One row corresonds to one VERSION and one FITOPT;
        # combine the SPLIT output tables into one merged table file.
        # Merge separately for each table format; HBOOK, ROOT, TEXT ...
        # Also check for symLink option to replace table-merge
        # with sym link to FITOPT000.

        # init name of merged table file for each format
        self.config_prep['merge_table_file_list'] = [''] * NTABLE_FORMAT

        row         = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        state       = row[COLNUM_FIT_MERGE_STATE]
        version     = row[COLNUM_FIT_MERGE_VERSION]
        fitopt_num  = row[COLNUM_FIT_MERGE_FITOPT]  # e.g., FITOPT001
        nevt_expect = row[COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS]

        submit_info_yaml = self.config_prep['submit_info_yaml']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_table_format = submit_info_yaml['USE_TABLE_FORMAT']
        msgerr           = []

        version_fitopt = f"{version}_{fitopt_num}"

        # start with idiot check to make sure that state is DONE.
        if state != SUBMIT_STATE_DONE :
            msgerr.append(f"Unexpected STATE = {state}; " \
                          f"expected {SUBMIT_STATE_DONE} state.")
            msgerr.append(f"Check {version_fitopt} files.")
            self.log_assert(False,msgerr)

        logging.info(f" Begin table merge for {version_fitopt} ")

        # check for sym link first
        if self.create_sym_link_tables(version,fitopt_num) :
            return

        # - - - - -
        # make sure we have correct number of files to merge
        for itable in range(0,NTABLE_FORMAT):
            fmt         =  TABLE_FORMAT_LIST[itable]
            if itable == ITABLE_TEXT :
                suffix = f"{SUFFIX_FITRES}.{fmt}"  # special case for TEXT
            else:
                suffix = fmt
            use_format     =  use_table_format[itable] 
            f_wildcard     =  f"{script_dir}/{version_fitopt}*{suffix}"
            if use_format :
                util.check_file_count(n_job_split, f_wildcard )

        # - - - - -
        # create dictionary to pass to merge_table_XXX functions below.
        version_fitopt_dict = {
            'version'         : version,
            'fitopt'          : fitopt_num,
            'version_fitopt'  : version_fitopt,
            'nevt_expect'     : nevt_expect 
        }

        # - - - - -
        # merge strategy is to merge split files into temp files 
        # named MERGE_{version_fitopt}.{suffix}; in case of abort,
        # these temp files remain. After NEVT validation,
        # each temp file is moved to {version}/{fitopt}.{suffix}

        # only SNANA has root and hbook
        if LCFIT_SUBCLASS == LCFIT_SNANA:
            if use_table_format[ITABLE_HBOOK] :
                self.merge_table_CERN(ITABLE_HBOOK, version_fitopt_dict)

            if use_table_format[ITABLE_ROOT] :
                self.merge_table_CERN(ITABLE_ROOT, version_fitopt_dict)

        # all subclass must have TEXT format.
        # Process TEXT format after ROOT & HBOOK to allow for append feature
        if use_table_format[ITABLE_TEXT] :
            for table_name in TABLE_NAME_LIST :
                self.merge_table_TEXT(table_name, version_fitopt_dict)

        # move MERGE files, and remove 'MERGE_' prefix
        self.move_merge_table_files(version_fitopt_dict)

        # Nov 15 2022
        # compress this version_fitopt here to reduce file count for big jobs.
        
        n_job_split  = submit_info_yaml['N_JOB_SPLIT']
        n_job_link   = submit_info_yaml['N_JOB_LINK']
        do_compress = n_job_split > 3  and n_job_link == 0

        if do_compress :
            suffix_tar_list = self.get_suffix_tar_list_lcfit()
            logging.info(f"  Compress {version_fitopt} for {suffix_tar_list}")
            for suffix in suffix_tar_list :
                wildcard = f"{version_fitopt}*SPLIT*.{suffix}"
                tar_file = f"{version_fitopt}_SPLITALL_{suffix}.tar"
                util.compress_files(+1, script_dir, wildcard, tar_file, "" )    
 
        logging.info("")

        if irow == 33333:
            sys.exit(f"\n\t xxxxx DEBUG DIE from fit wrapup ... xxxx ")

        return

        # end merge_job_wrapup 

    def get_suffix_tar_list_lcfit(self):

        submit_info_yaml      = self.config_prep['submit_info_yaml']
        use_table_format      = submit_info_yaml['USE_TABLE_FORMAT']

        suffix_tar_list = JOB_SUFFIX_TAR_LIST.copy()
        for use, suffix in zip(use_table_format,TABLE_FORMAT_LIST):
            if use:
                suffix_tar_list.append(suffix)

        return suffix_tar_list
        # end get_suffix_tar_list_lcfit

    def create_sym_link_tables(self, version, fitopt_num):

        # if this fitopt_num (e.g., FITOPT003) is on list of sym links,
        # create sym link for each merged table expected to appear 
        # in ../version.
        # Function returns true of sym links are create; false otherwise.

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']

        use_table_format    = submit_info_yaml['USE_TABLE_FORMAT']
        link_FITOPT000_list = submit_info_yaml['LINK_FITOPT000_LIST']
        cdv                 = f"cd {output_dir}/{version}"
        create_links        = False  # init return function arg

        if fitopt_num in link_FITOPT000_list :
            fitopt_ref = f"{FITOPT_STRING}000"
            cmd_link_all = ''

            for itab in range(0,NTABLE_FORMAT):
                if not use_table_format[itab] : continue
                suf   =  TABLE_FORMAT_LIST[itab]
                if itab == ITABLE_TEXT :   suf  = f"{SUFFIX_FITRES}"  # special case

                suf += '.gz'
                cmd_link = f"ln -s {fitopt_ref}.{suf} {fitopt_num}.{suf} ; "
                cmd_link_all += cmd_link
                logging.info(f"   create sym-link to {fitopt_ref}.{suf}")

            cmd = f"{cdv} ; {cmd_link_all}"
            create_links = True
            os.system(cmd)
        
        return create_links
        # end create_sym_link_tables

    def merge_table_TEXT(self, table_name, version_fitopt_dict):

        # Input table_name is either 'FITRES', 'SNANA', or 'LCPLOT'.
        # Problem with TEXT tables is that each table must be written
        # to separate set of files, so here each table name is processed
        # separately. Usually only the FITRES table is requested, but 
        # sometimes other tables are included.
        # HBOOK & ROOT don't have this issue because all tables
        # reside in one file.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        version          = version_fitopt_dict['version']
        fitopt           = version_fitopt_dict['fitopt']
        version_fitopt   = version_fitopt_dict['version_fitopt']
        nevt_expect      = version_fitopt_dict['nevt_expect']

        itable           = ITABLE_TEXT
        suffix           = TABLE_FORMAT_LIST[itable]
        prefix           = PREFIX_MERGE  
        table_wildcard  = f"{version_fitopt}*{table_name}.TEXT"

        # check debug option to force merge-table failure
        flag_force_fail = \
            self.flag_force_merge_table_fail(itable,version_fitopt)

        if flag_force_fail == FLAG_FORCE_MERGE_TABLE_CORRUPT :
            table_wildcard += f(" {table_wildcard}") # double output

        # if no tables exist, bail out
        table_list = sorted(glob.glob1(script_dir,table_wildcard))
        if len(table_list) == 0 :
            return

        # for FITRES table, check that all of the VARNAMES lists 
        # are the same; else abort
        if table_name == SUFFIX_FITRES :
            self.check_table_varnames_TEXT(table_list)

        # construct linux command to catenate TEXT files,
        # and then a special awk command to remove all VARNAMES
        # lines EXCEPT for the first VARNAMES line.
        # Note that output extension is table name, not TEXT
        cddir           = f"cd {script_dir}"
        out_table_file  = f"{prefix}_{version_fitopt}.{table_name}"
        out_table_file2 = f"{prefix}2_{version_fitopt}.{table_name}" # for awk

        cmd_cat  = f"cat {table_wildcard} > {out_table_file}"
        cmd_awk  = f"awk '!/^VARNAMES/ || ++n <= 1' {out_table_file} > " \
                   f"{out_table_file2} ; " \
                   f"mv {out_table_file2} {out_table_file}"

        msg = f"   merge {n_job_split} {suffix}-{table_name} table files."
        logging.info(msg)

        if flag_force_fail != FLAG_FORCE_MERGE_TABLE_MISSING :
            tref = datetime.datetime.now()
            cmd_all = f"{cmd_cat} ; {cmd_awk}"
            os.system(f"{cddir} ; {cmd_all}")
            util.print_elapse_time(tref,f"merge {n_job_split} table files")

        OUT_TABLE_FILE = f"{script_dir}/{out_table_file}"
        self.check_file_exists(OUT_TABLE_FILE,["Problem with table-merge"])

        # - - - - - -
        # for FITRES table only, do unitarity check to make sure that 
        # nevt_expect rows are really there. Also check options to 
        # append variables to table.

        if table_name == SUFFIX_FITRES :
            OUT_TABLE_FILE = f"{script_dir}/{out_table_file}"
            nevt_find, n_nan = util.nrow_table_TEXT(OUT_TABLE_FILE,"SN:")
            if n_nan > 0:
                msgerr = [f"found {n_nan} nan in {OUT_TABLE_FILE}"]
                self.log_assert(False,msgerr) 

            self.nevt_table_check(nevt_expect, nevt_find, out_table_file,
                                  cmd_all)
            self.config_prep['merge_table_file_list'][itable] = out_table_file

            # check options to append FITRES file
            # 1. APPEND_TABLE_VARLIST  -> extract vars from HBOOK or ROOT file
            # 2. APPEND_TABLE_TEXTFILE -> append vars from external file.
            self.append_table_varlist(version_fitopt_dict)  # optional
            self.append_table_textfile(version_fitopt_dict)  # optional

        return
        # end merge_table_TEXT

    def check_table_varnames_TEXT(self,table_list):

        # Created Aug 31 2022
        # read VARNAMES list for each file in table_list,
        # and abort if any VARNAMES list differs from first table.
        # Beware that this works only for key (FITRES) format;
        # does not work for CSV formatted tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        msgerr           = []
        n_list = len(table_list)
        if n_list == 1 : return

        KEY_VARNAMES   = 'VARNAMES:'
        VARNAMES_EMPTY = 'EMPTY'  # for empty file where 0 events pass

        varnames_list = []
        for tb_file in table_list:
            tb_file_full = f"{script_dir}/{tb_file}"
            with open(tb_file_full,"rt") as t:
                line_list = t.readlines()
                if len(line_list) == 0:
                    # flag empty FITRES file to avoid false error
                    varnames_list.append(VARNAMES_EMPTY)
                else:
                    for line in line_list:
                        if KEY_VARNAMES in line:
                            varnames = line.split(KEY_VARNAMES)[1]
                            varnames_list.append(varnames)
                            break

        # - - - - - - - 
        n2_list = len(varnames_list)
        if n2_list != n_list:
            msgerr.append(f"Found {n2_list} VARNAMES keys, " \
                          f"but expected {n_list}. ")
            msgerr.append(f"table_list = \n\t{table_list}")
            msgerr.append(f"varnames_list = \n\t{varnames_list}")
            self.log_assert(False,msgerr) 
            
        # find varnames_ref to be first varnames that isn't empty
        varnames_ref = None
        for varnames in varnames_list:
            if varnames != VARNAMES_EMPTY:
                varnames_ref = varnames
                break

        # xxx mark varnames_ref = varnames_list[0] # this is a comma-sep string
        nerr = 0
        for varnames,tb_file in zip(varnames_list, table_list):
            is_empty = (varnames == VARNAMES_EMPTY)
            if is_empty:
                is_match = True
            else:
                is_match = self.match_varnames_list(varnames_ref, varnames)

            if not (is_empty or is_match):
                logging.info(f"ERROR: VARNAMES mis-match for {tb_file}")
                nerr += 1

        # - - - - - 
        if nerr > 0 :
            msgerr.append(f"Found {nerr} VARNAMES mis-matches " \
                          f"among {n_list} split jobs.")
            msgerr.append(f"Reference FITRES file is: {table_list[0]}")
            self.log_assert(False,msgerr) 

        return
        # end check_table_varnames_TEXT

    def match_varnames_list(self,varnames_0, varnames_1):
        # Created Sep 2022 
        # Returns True of the two sets of varnames has the same length
        # and all elements match except those with SIMNULL in the name.
        # For example, a mix of SNIa-SALT2 + SNCC sim will have SIM_SALT2c
        # in the SALT2 fitres file, and SIMNULL2 (or SIM_AV) in the SNCC
        # fitres file. These don't match, but this particular mis-match 
        # should not cause a problem.

        is_match = True
        KEYSTR_SKIP_LIST = [ "SIMNULL", "SIM_SALT2", "SIM_AV", "SIM_RV" ]

        vname_list_0 = varnames_0.split()
        vname_list_1 = varnames_1.split()
        len0 = len(vname_list_0)
        len1 = len(vname_list_1)
        if len0 != len1: 
            logging.info(f" VARNAME Match ERROR: len0={len0} len1={len1}")
            logging.info(f" VARNAME list0 = {vname_list_0}")
            logging.info(f" VARNAME list1 = {vname_list_1}")
            return False

        for vname_0, vname_1 in zip(vname_list_0, vname_list_1):

            # check for keys to skip because they differ for SNIa-SALT2
            # and non-SNIA
            skip = False
            for key_skip in KEYSTR_SKIP_LIST:
                if key_skip in vname_0: skip = True
                if key_skip in vname_1: skip = True
            if skip : 
                continue

            if vname_0 != vname_1 :
                logging.info(f" VARNAME Match ERROR: {vname_0} != {vname_1}")
            
                is_match = False

        return is_match
        # end match_varnames_list

    def append_table_varlist(self,version_fitopt_dict) :

        # Check option to extract variables from HBOOK/ROOT file,
        # and append TEXT-FITRES file.
        # See CONFIG key APPEND_TABLE_VARLIST
        #

        CONFIG                = self.config_yaml['CONFIG']
        submit_info_yaml      = self.config_prep['submit_info_yaml']
        merge_table_file_list = self.config_prep['merge_table_file_list']
        script_dir            = submit_info_yaml['SCRIPT_DIR']
        use_table_format      = submit_info_yaml['USE_TABLE_FORMAT']
        version_fitopt        = version_fitopt_dict['version_fitopt']
        nevt_expect           = version_fitopt_dict['nevt_expect']
        
        # use CONFIG from input file (rather than SUBMIT.INFO) in case
        # user fixes APPEND_TABLE_VARLIST after jobs ran.
        key  = KEY_APPEND_TABLE_VARLIST
        if key in CONFIG:
             varlist_append = CONFIG[key]
        else:
            return

        #if key in submit_info_yaml :
        #    varlist_append = submit_info_yaml[key]
        #else:
        #    return

        msgerr = []

        text_table_file = merge_table_file_list[ITABLE_TEXT]  # file to append
        if use_table_format[ITABLE_ROOT] :
            full_table_file = merge_table_file_list[ITABLE_ROOT]
        elif use_table_format[ITABLE_HBOOK] :
            full_table_file = merge_table_file_list[ITABLE_HBOOK]
        else :
            msgerr.append("Cannot append TEXT table without full table")
            msgerr.append("Must define HBOOK or ROOT to append TEXT table")
            self.log_assert(False,msgerr) 

        append_log_file = f"sntable_append_{PREFIX_MERGE}_{version_fitopt}.log"
        append_out_file = f"sntable_append_{PREFIX_MERGE}_{version_fitopt}.text"

        tref = datetime.datetime.now()
        logging.info(f"   Append TEXT table from {full_table_file} ")

        cddir = f"cd {script_dir}"
        cmd_append = f"{SCRIPT_SNTABLE_DUMP} {full_table_file} FITRES " \
                     f"--VERBOSE " \
                     f"-v '{varlist_append}' " \
                     f"-a {text_table_file} > {append_log_file} 2>/dev/null"


        # abort if required C program exe is not found within 10 minutes        
        snana_dir      = self.config_yaml['args'].snana_dir
        Cprogram_path  = util.get_SNANA_program_path(snana_dir,Cprogram_SNTABLE_DUMP)
        found_Cprogram = util.program_exists(Cprogram_path, TMAX_EXE_WAIT_ABORT, True) 
        
        istat = os.system(f"{cddir} ; {cmd_append} ")
        util.print_elapse_time(tref,SCRIPT_SNTABLE_DUMP)

        if istat != 0 :
            msgerr.append(f"Failed to apppend variables with command")
            msgerr.append(f"   {cmd_append}")
            msgerr.append(f"See {append_log_file}")
            self.log_assert(False,msgerr) 

        # replace FITRES file with append file, 
        # and remove sntable* junk files

        tref = datetime.datetime.now()
        if os.path.isfile(f"{script_dir}/{append_out_file}") :
            cmd_mv = f"mv {append_out_file} {text_table_file}"
            cmd_rm = f"rm sntable_*"
            cmd    = f"{cddir}; {cmd_mv} ; {cmd_rm}"
            os.system(cmd)
        else : 
            msgerr.append(f"Append FITRES table failed with command")
            msgerr.append(f"{cmd_append}")
            msgerr.append(f"Check if these variables are all valid:")
            msgerr.append(f"    {varlist_append}")
            self.log_assert(False,msgerr) 

        # finally, make sure that the number of rows still matches nevt_expect
        OUT_TABLE_FILE = f"{script_dir}/{text_table_file}"
        nevt_find, n_nan = util.nrow_table_TEXT(OUT_TABLE_FILE,"SN:")

        if n_nan > 0 :
            msgerr = [f"found {n_nan} nan in {OUT_TABLE_FILE}"]
            self.log_assert(False,msgerr) 

        self.nevt_table_check(nevt_expect, nevt_find, text_table_file,
                              cmd_append )
                  
        util.print_elapse_time(tref,"append cleanup & validate")

        # end append_table_varlist

    def append_table_textfile(self,version_fitopt_dict) :

        # Check option to extract and append variables from 
        # external TEXT file. See CONFIG key APPEND_TABLE_TEXFILE
        #
        submit_info_yaml      = self.config_prep['submit_info_yaml']
        merge_table_file_list = self.config_prep['merge_table_file_list']
        script_dir            = submit_info_yaml['SCRIPT_DIR']
        use_table_format      = submit_info_yaml['USE_TABLE_FORMAT']
        version_fitopt        = version_fitopt_dict['version_fitopt']
        nevt_expect           = version_fitopt_dict['nevt_expect']
        
        key  = KEY_APPEND_TABLE_TEXTFILE
        if key in submit_info_yaml :
            external_file_list = submit_info_yaml[key]
        else:
            return

        orig_file = merge_table_file_list[ITABLE_TEXT]  # file to append
        out_file  = "TEMP_" + orig_file
        log_file  = "TEMP_COMBINE.LOG"

        cddir = f"cd {script_dir}"
        cmd1  = f"{PROGRAM_COMBINE_FITRES} {orig_file} {external_file_list} "\
                f"-outfile_text {out_file} " \
                f"t " \
                f"-nullval_float {NULLVAL_COMBINE_FITRES} " \
                f">& {log_file} " 
        cmd2  = f"mv {out_file} {orig_file}"
        cmd3  = f"rm {log_file}"
        cmd   = f"{cmd1} ; {cmd2} ; {cmd3}"
        os.system(f"{cddir} ; {cmd}" )

        # finally, make sure that the number of rows still matches nevt_expect
        ORIG_FILE = f"{script_dir}/{orig_file}"
        nevt_find, n_nan = util.nrow_table_TEXT(ORIG_FILE,"SN:")
        self.nevt_table_check(nevt_expect, nevt_find, orig_file, cmd)

        #sys.exit('\n xxx DEBUG DIE xxxx ')

    # end append_table_textfile

    def merge_table_CERN(self, itable, version_fitopt_dict):

        # call merge program for HBOOK or ROOT based on itable arg.
        # into one.
        version          = version_fitopt_dict['version']
        fitopt           = version_fitopt_dict['fitopt']
        version_fitopt   = version_fitopt_dict['version_fitopt']
        nevt_expect      = version_fitopt_dict['nevt_expect']

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        suffix           = TABLE_FORMAT_LIST[itable]
        prefix           = PREFIX_MERGE
        msgerr           = []
        tref = datetime.datetime.now()

        # check debug option to force merge-table failure
        flag_force_fail = self.flag_force_merge_table_fail(itable,version_fitopt)

        logging.info(f"   merge {n_job_split} {suffix} table files.")
        
        # get name of program
        program_merge = f"merge_{suffix.lower()}.exe"
        self.check_program_merge_table_CERN(program_merge)

        cddir           = f"cd {script_dir}"
        out_table_file  = f"{prefix}_{version_fitopt}.{suffix}"
        log_table_file  = f"{prefix}_{version_fitopt}_{suffix}.LOG"
        f_list          = f"{version_fitopt}*{suffix}"
        if flag_force_fail == FLAG_FORCE_MERGE_TABLE_CORRUPT :
            f_list += f" XXX_FORCE_CORRUPT.{suffix}"

        cmd_merge       = f"{program_merge} {f_list} {out_table_file} "
        cmd             = f"{cmd_merge} &> {log_table_file}  "

        if flag_force_fail != FLAG_FORCE_MERGE_TABLE_MISSING :
            os.system(f"{cddir}; {cmd}" )

        OUT_TABLE_FILE = f"{script_dir}/{out_table_file}"
        self.check_file_exists(OUT_TABLE_FILE, ["Problem with table-merge"])

        # for HBOOK, remove garbage from log file
        if itable == ITABLE_HBOOK :
            cmd_clean_log = \
                f"{cddir} ; remove_locf_messages.py {log_table_file} -q"
            os.system(cmd_clean_log)

        # - - - -
        util.print_elapse_time(tref,f"merge {n_job_split} table files")
        

        # ?? check log file for success message ??
        tref = datetime.datetime.now()

        # extract number of events in final HBOOK/ROOT file
        NTRY_MAX = 2; ntry=0; nevt_find = -9
        while nevt_find < 0 and ntry < NTRY_MAX :
            nevt_find = self.nrow_table_CERN(f"{OUT_TABLE_FILE}")
            ntry += 1

        if nevt_find < 0 :
            msgerr.append(f"Cannot get NEVT from merged {suffix} table " \
                          f"after {ntry} attempts; ")
            msgerr.append(f"   {out_table_file}")
            msgerr.append(f"---> table may be corrupted.")
            msgerr.append(f" ")
            msgerr.append(f"Full table-merge command was")
            msgerr.append(f"  {cmd}")
            msgerr.append(f" ")
            self.log_assert(False,msgerr)

        if ntry > 1 :
            msg=(f"\tWARNING: ntry={ntry} to get NEVT for {out_table_file}")
            logging.warning(msg)

        # abort if nevt_find != nevt_expect
        self.nevt_table_check(nevt_expect, nevt_find, out_table_file, cmd)

        # store merge table file for append_table_text()
        merge_table_file_list =  self.config_prep['merge_table_file_list']
        merge_table_file_list[itable] = out_table_file

        util.print_elapse_time(tref,f"validate NEVT in merge table(s)")

        # end merge_table_CERN

    def check_program_merge_table_CERN(self,program_merge):
        # Created Apr 2022 by R.Kessler
        # wait for program_merge to exist as executable (allow for make during processing)
        # Abort if wait is too long.
        
        snana_dir     = self.config_yaml['args'].snana_dir
        program_path  = util.get_SNANA_program_path(snana_dir,program_merge)
        found_program = util.program_exists(program_path, TMAX_EXE_WAIT_ABORT, True)
                
        return
        # end check_program_merge_table_CERN

    
    def nrow_table_CERN(self,table_file):
        # return number of table rows in CERN file that
        # has either HBOOK or ROOT extension
        # Return -9 on error.
        # Script prints "NEVT:  <nevt>", so parse the 2nd element.

        script   = SCRIPT_SNTABLE_DUMP
        arg_NEVT = "--NEVT"
        if '.pl' in script:  arg_NEVT = "NEVT" # legacy perl arg

        # abort if required C program exe is not found within 10 minutes        
        snana_dir      = self.config_yaml['args'].snana_dir
        Cprogram_path  = util.get_SNANA_program_path(snana_dir,Cprogram_SNTABLE_DUMP)
        found_Cprogram = util.program_exists(Cprogram_path, TMAX_EXE_WAIT_ABORT, True) 
        
        cmd_nevt = f"{script} {table_file} FITRES {arg_NEVT} " \
                   f" | grep 'NEVT:' "
        try: 
            result_line = subprocess.check_output(cmd_nevt, shell=True)
            result_line = (result_line.rstrip()).decode('utf-8')
            nevt_find   = int(result_line.split()[1])
            # nevt_find = 0 # xxx REMOVE
            return nevt_find
        except :
            return -9

        #end nrow_table_CERN

    def move_merge_table_files(self,version_fitopt_dict):
       
        # move files of the form  MERGE_FITOPTnnn.FORMAT
        # to VERSION/FITOPTnnn.FORMAT
        # (i.e., strip of the MERGE_ prefix that is for debugging)

        version          = version_fitopt_dict['version']
        fitopt           = version_fitopt_dict['fitopt']

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        table_wildcard   = f"{PREFIX_MERGE}_*"
        table_list       = glob.glob1(script_dir,table_wildcard)

        cmd = f"cd {script_dir} ; "
        for tfile in table_list :
            jdot   = tfile.index('.')
            suffix = tfile[jdot+1:]
            if suffix == 'LOG' :
                cmd += f"rm {tfile} ; "
            else:
                move_file = f"../{version}/{fitopt}.{suffix}"
                cmd      += f"mv {tfile} {move_file} ; "

        os.system(cmd)

        # end move_merge_table_files

    def nevt_table_check(self, nevt_expect, nevt_find, out_file, cmd):
        # abort if nevt_expect != nevt_find
        # Apr 25 2021: pass command {cmd} and print it on error.
        if nevt_expect != nevt_find :
            msgerr = []
            msgerr.append(f"Found {nevt_find} events in merged {out_file}")
            msgerr.append(f"but expected {nevt_expect} from summing YAML files.")            
            msgerr.append(f"Table merge has a problem. Check command")
            msgerr.append(f"  {cmd}")
            self.log_assert(False,msgerr)
    # end nevt_table_check

    def merge_cleanup_final(self):
        
        # every snlc_fit job succeeded, so here we simply compress output.
        # Before compressing, check that every expected merged table
        # file exists because previous error checking is for snlc_fit,
        # not for merging tables.
        #
        # If any merged table file is missing, 
        #  + print list of missing table files in MERGE.LOG (and to stdout)
        #  + abort and change SUCCESS state to FAIL in global DONE file
        #

        output_dir          = self.config_prep['output_dir']
        submit_info_yaml    = self.config_prep['submit_info_yaml']
        use_table_format    = submit_info_yaml['USE_TABLE_FORMAT']
        link_FITOPT000_list = submit_info_yaml['LINK_FITOPT000_LIST']
        script_dir          = submit_info_yaml['SCRIPT_DIR']
        subdir         = SUBDIR_SCRIPTS_LCFIT
        tar_file       = f"{subdir}.tar"

        logging.info(f" FIT cleanup: check if all merged tables exist.")

        # check that every expected merged table file exits;
        # read status from MERGE file to get list of version & fitopts
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        MERGE_INFO_CONTENTS,comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        row_list  = MERGE_INFO_CONTENTS[TABLE_MERGE]
        msgerr = ['ERROR merging table files: ' ]
        nerr   = 0
        for row in row_list :
            version = row[COLNUM_FIT_MERGE_VERSION]
            fitopt  = row[COLNUM_FIT_MERGE_FITOPT]  # e.g., FITOPT002
            if fitopt in link_FITOPT000_list : continue # skip sym links

            for itable in range(0,NTABLE_FORMAT):
                use    = use_table_format[itable]
                suffix = TABLE_FORMAT_LIST[itable]
                if itable == ITABLE_TEXT:
                    suffix = SUFFIX_FITRES  # special case here
                table_file = f"{version}/{fitopt}.{suffix}"
                TABLE_FILE = f"{output_dir}/{table_file}"
                if use and not os.path.isfile(TABLE_FILE):
                    nerr += 1
                    msgerr.append(f"    Missing expected {table_file}")
                
        if nerr > 0 :
            self.log_assert(False,msgerr) # abort and write FAIL to global DONE file.
        else :
            logging.info(f" FIT cleanup: all merged tables exist.")

        # if we get here, table-merging seems to have worked so tar and zip
        # use *SPLIT*[suffix] so that it works on *SPLIT*.[suffix]
        # and also on *SPLIT_[suffix].tar

        suffix_tar_list = self.get_suffix_tar_list_lcfit()
        logging.info(f" LCFIT cleanup: tar {suffix_tar_list} under {subdir}/")
        for suffix in suffix_tar_list :
            wildcard = f"*SPLIT*{suffix}*"
            util.compress_files(+1, script_dir, wildcard,  suffix, "" )    

        logging.info(f" FIT cleanup: gzip merged tables.")
        cmd_gzip = f"cd {output_dir} ; gzip */FITOPT* 2>/dev/null"
        os.system(cmd_gzip)

        self.merge_cleanup_script_dir() 

        # xxx nothing to write to? logging.info(f" FIT cleanup: Done.")

        # end merge_cleanup_final

    def merge_reset(self,output_dir):

        # remove all merge products to allow re-doing the merge with -m option.
        # Renders output_dir as if --nomerge option had been given with
        # original submit command.
        #  + set all STATE values to WAIT
        #  + set all NEVT to 0
        #  + set CPU to 0
        #  + unpack SPLIT_JOBS_LCFIT
        #  + remove merged table files: [VERSION]/FITOPT* files
        #  + remove all-done file
        #
        # Apr 23 2021: check that file/subdir exists before removing it.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        version_list     = submit_info_yaml['VERSION_LIST']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        fnam = "merge_reset"

        # read status from MERGE file          
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        colnum_zero_list    = [ COLNUM_FIT_MERGE_NEVT_ALL, 
                                COLNUM_FIT_MERGE_NEVT_LC_CUTS,
                                COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS,
                                COLNUM_FIT_MERGE_CPU]
        logging.info(f"  {fnam}: STATE->WAIT and NEVT->0 in {MERGE_LOG_FILE}")
        util.merge_table_reset(MERGE_LOG_PATHFILE, TABLE_MERGE,  \
                               COLNUM_MERGE_STATE, colnum_zero_list)

        # loop over each version and clean out FITOPT* files 
        logging.info(f"  {fnam}: remove FITOPT* from version subdirs")
        for version in version_list :
            v_dir  = f"{output_dir}/{version}"
            cmd_rm = f"rm {v_dir}/FITOPT* 2>/dev/null"
            os.system(cmd_rm)

        # if script_dir is tarred & gzipped, unpack it
        logging.info(f"  {fnam}: unapck {SUBDIR_SCRIPTS_LCFIT}")
        util.compress_subdir(-1, script_dir)

        # untar and unzip file inside SUBDIR_SCRIPTS_LCFIT
        util.untar_script_dir(script_dir)

        # remove lingering temp files in SPLIT_JOBS_LCFIT.
        # Also remove MERGE.LOG_{Nsec} backups, and DONE file

        # define list of wildcards to remove under script_dir and output_dir
        wildcard_script_dir = [ f"{PREFIX_MERGE}*", "sntable*" ]
        wildcard_output_dir = [ f"{MERGE_LOG_FILE}_*", "*.DONE", "BUSY*" ]

        logging.info(f"  {fnam}: remove misc junk files from " \
                     f"{SUBDIR_SCRIPTS_LCFIT}")
        for wildcard in wildcard_script_dir :
            if len(wildcard) < 2: continue  # avoid accidental rm *
            if len(glob.glob1(script_dir,wildcard)) > 0 :
                   cmd_rm = f"rm {script_dir}/{wildcard}"
                   print(f"\t Remove {wildcard} ")
                   os.system(cmd_rm)

        output_base = os.path.basename(output_dir)
        logging.info(f"  {fnam}: remove misc junk files from {output_base}")
        for wildcard in wildcard_output_dir :
            if len(wildcard) < 2: continue  # avoid accidental rm *
            if len(glob.glob1(output_dir,wildcard)) > 0 :
                cmd_rm = f"rm {output_dir}/{wildcard}"
                print(f"\t Remove {wildcard}  ")
                os.system(cmd_rm)

        # end merge_reset

    def flag_force_merge_table_fail(self, itable, version_fitopt):

        # Check CONFIG (user input) for 
        #   FORCE_MERGE_TABLE_MISSING({table_format})
        #   FORCE_MERGE_TABLE_CORRUPT({table_format})
        #
        # and return appropriate flag. This flag is used bu
        # table-merge function to force failure.
        # Goal is to debug error handling on table-merge failures.

        flag   = 0 # init output
        CONFIG = self.config_yaml['CONFIG']

        #if 'FORCE' not in CONFIG :
        #    return flag

        suffix = TABLE_FORMAT_LIST[itable]

        key_list = [ f"FORCE_MERGE_TABLE_MISSING({suffix}", 
                     f"FORCE_MERGE_TABLE_CORRUPT({suffix}" ]
        flag_list = [ FLAG_FORCE_MERGE_TABLE_MISSING , 
                      FLAG_FORCE_MERGE_TABLE_CORRUPT  ]

        nkey = len(key_list)

        for ikey in range(0,nkey):
            tmp_key  = key_list[ikey]
            tmp_flag = flag_list[ikey]
            if tmp_key in CONFIG :
                for item in CONFIG[tmp_key] :
                    if item == version_fitopt:
                        flag = tmp_flag
                        msg = f"force {suffix}-table fail for {version_fitopt}"
                        logging.info(msg)
        return flag
        # end force_merge_table_fail

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file .
        # Each line is of the form
        #   KEYNAM:  VALUE

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        opt_sncid_list = 0 
        for key in KEY_OPT_SNCID_LIST :
            if key in submit_info_yaml :
                opt_sncid_list = submit_info_yaml[key]

        survey,idsurvey  = util.get_survey_info(script_dir,"*SPLIT*.YAML")

        info_lines  = []
        info_lines.append(f"SURVEY:         {survey}")
        info_lines.append(f"IDSURVEY:       {idsurvey}")

        # May 19 2021: count number of common events for option to 
        # use FITOPT000 events in all FITOPTs
        if opt_sncid_list > 0 :
            version = submit_info_yaml['VERSION_LIST'][0]
            nevt_common = self.get_nevt_common(version)
            info_lines.append(f"NEVT_COMMON:    {nevt_common}     " \
                              f"# among FITOPTs in {version}")

        return info_lines
        # end get_misc_merge_info    

    def get_merge_COLNUM_CPU(self):
        return COLNUM_FIT_MERGE_CPU

    def get_nevt_common(self,version):
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        FITOPT_LIST      = submit_info_yaml['FITOPT_LIST']
        
        fitres_dir = f"{output_dir}/{version}"
        combined   = None

        for row in FITOPT_LIST:
            num   = row[0]  # e.g., FITOPT001
            label = row[1]  
            fitres_file = f"{num}.{SUFFIX_FITRES}"
            FITRES_FILE = f"{fitres_dir}/{fitres_file}"
            NOREJECT    = FITOPT_STRING_NOREJECT in label
            if NOREJECT : continue
            #print(f" xxx process {FITRES_FILE}")
            df = pd.read_csv(FITRES_FILE, delim_whitespace=True, comment="#")
            if combined is None:
                combined = df.index
            else:
                combined = combined.intersection(df.index)

        nevt_common = combined.shape[0]
        
        #sys.exit('\n xxx DEBUG DIE xxxx ')
        return nevt_common
        # end get_nevt_common


