# Created July 2020 by R.Kessler
#
# Improvements w.r.t. SALT2mu_fit.pl:
#
#  + if there is 1 and only 1 version per INPDIR, there is no
#    need to specify STRINGMATCH or STRINGMATCH_IGNORE ... works
#    for RANSEED_REPEAT and RANSEED_CHANGE.
#
#  + If no errors are detected, everything is gzipped by default.
#
#  + Added MERGE.LOG to monitor STATE(WAIT,RUN,DONE,FAIL) and to
#    monitor stats for DATA,BIASCOR,CCPRIOR. This MERGE.LOG
#    similar to that for sim and fit jobs.
#
#  + output from optional wfit is in YAML format for easier parsing.
#
#  + NSPLITRAN runs on all FITOPT & MUOPT ... beware of large N_JOB_TOT
#
#  + new FITOPTxMUOPT input to select subset of FITOPTxMUOPT matrix.
#     (see help with -H BBC arg)
#
#  + Automatically creates summary files:
#     BBC_REJECT_SUMMRY.LIST      -> SN rejected by some FITOPT/MUOPT
#     BBC_SUMMARY_wfit.FITRES     -> summary of wfit results
#     BBC_SUMMARY_SPLITRAN.FITRES -> NSPLITRAN summary of AVG, RMS, ...
#
#  + CPU diagnostics added to MERGE.LOG 
#
# - - - - 
# Potential Upgrades:
#   read list of BIASCOR outdirs from fit process, and automatically 
#   fill in simfile_biascor arg. Same for simfile_ccprior.
#
#
#        HISTORY
#
# Dec 02 2020: add SUFFIX_COV to list of files to move
# Dec 17 2020: update get_matrix_FITOPTxMUOPT to process multiple FITOPTxMUOPT
# Jan 12 2021: write BBC_ACCEPT_SUMMARY for CIDs in all FITOPT*.FITRES files.
# Mar 08 2021: if INPDIR+: None, use argument of datafile=
# Apr 23 2021: abort if any version dir does not exist.
# May 24 2021: check option to use events from FITOPT000
# May 27 2021: new def make_FITOPT_OUT_LIST 
#                 (append_fitopt_info_file is obsolete)
# Aug 09 2021: 
#     use def get_wfit_values() to handel legacy vs. refact yaml keys.
#     Beware that new w0,wa model in wfit is NOT handled here [yet]
# Sep 19 2021: write NEVT_bySAMPLE in BBC_SUMMARY_FITPAR.YAML
# Sep 28 2021: include wa if -wa arg is specified for wfit
# Oct 05 2021: move get_wfit_values to submit_util.py so that
#                dedicated wfit class can use it too.
# Nov 17 2021: add list protection in def make_fitpar_summary()
# Nov 24 2021: write OLAM_REF and w_REF to submit info
# Jan 18 2022: fix writing REJECT_FRAC_BIASCOR to yaml file.
# Mar 03 2022: add zPRIOR* to append_varname_missing 
# Mar 28 2022: write IZBIN to BBC_ACCEPT summary file
# Apr 08 2022: 
#   + fix missing-IZBIN bug for M11-style fit without biasCor.
#   + local sync_evt is no True/False instead of 1/0
#   + use merge_force logic on FITOPT000 if sync_evt is set
#
# Aug 18 2022: change 4D loop order so that iver is inside instead of outside;
#      --> ensures that FITOPT000 is always done first.
#
# Oct 01 2022 RK - minor refactor for merge_reset()
#
# ================================================================

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time
import submit_util as util
import numpy  as np
import pandas as pd

from submit_params    import *
from submit_prog_base import Program

USE_INPDIR = True

PREFIX_SALT2mu           = "SALT2mu"

# define FIT column to extract FIT VERSION for FIT-MERGE.LOG
COLNUM_FIT_MERGE_VERSION = 1  # same param in submit_prog_lcfit.py -> fragile

# define colums in BBC MERGE.LOG
COLNUM_BBC_MERGE_VERSION      = 1
COLNUM_BBC_MERGE_FITOPT       = 2
COLNUM_BBC_MERGE_MUOPT        = 3
COLNUM_BBC_MERGE_NEVT_DATA    = 4
COLNUM_BBC_MERGE_NEVT_BIASCOR = 5
COLNUM_BBC_MERGE_NEVT_CCPRIOR = 6
COLNUM_BBC_MERGE_SPLITRAN     = 7

# list used in wrapup, cleanup, and merge_reset
SUFFIX_MOVE_LIST = [ SUFFIX_FITRES, SUFFIX_M0DIF, SUFFIX_COV ]

# hard-wire output subDir name if there is 1-and-only-1 version
SUBDIR_OUTPUT_ONE_VERSION = "OUTPUT_BBCFIT"

# name of quick-and-dirty cosmology fitting program
PROGRAM_wfit = "wfit.exe"
PREFIX_wfit = "wfit"

FITPAR_SUMMARY_FILE   = "BBC_SUMMARY_FITPAR.YAML"   # Mar 28 2021
SPLITRAN_SUMMARY_FILE = "BBC_SUMMARY_SPLITRAN.FITRES"
WFIT_SUMMARY_FILE     = "BBC_SUMMARY_wfit.FITRES"

BBC_REJECT_SUMMARY_FILE  = "BBC_REJECT_SUMMARY.LIST"
BBC_ACCEPT_SUMMARY_FILE  = "BBC_ACCEPT_SUMMARY.LIST"

KEY_ROW               = "ROW:"

KEY_FITOPTxMUOPT      = 'FITOPTxMUOPT'
BLOCKNAME_FITOPT_MAP  = 'FITOPT_MAP'

#  Allow either of two keys
#                        pippin/submit key       key for snlc_fit
KEYLIST_SYNC_EVT     = [ 'FLAG_USE_SAME_EVENTS', 'OPT_SNCID_LIST' ]

MUOPT_STRING           = "MUOPT"
FITOPT_STRING_NOREJECT = "NOREJECT" # optional part of FITOPT label

OUTDIR_ITER1_SUFFIX    = "_ITER1"


KEYNAME_VARNAMES = "VARNAMES"
TABLE_VARNAME_CID      = "CID"
TABLE_VARNAME_IDSURVEY = "IDSURVEY"
TABLE_VARNAME_IZBIN    = "IZBIN"

# - - - - - - - - - - - - - - - - - - -  -
class BBC(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_BBC
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
            log_assert(False,msgerr)

        return output_dir_name,SUBDIR_SCRIPTS_BBC

    def submit_prepare_driver(self):
        logging.info("")

        # check for devel flag(s)
        devel_flag = self.config_yaml['args'].devel_flag

        # - - - - - - -
        # read C code inputs (not YAML block)
        self.bbc_read_input_file()

        # store list of BBC MUOPTs
        self.bbc_prep_muopt_list()

        self.bbc_prep_splitran()

        # read/store version and fitopt info
        self.bbc_prep_version_list()

        # figure out which versions to combine and create sorted lists
        self.bbc_prep_version_match() 

        # determine output FITOPTs based on FITOPT_MAP, or input FITOPTs
        self.bbc_prep_fitopt_outlist()

        # convet 2D and 4D nested loops into 1D loops
        self.bbc_prep_index_lists()

        # create output dir for each version or each version-splitran
        self.bbc_prep_mkdir()

        # copy & combine tables from INPDIR+ directories
        self.bbc_prep_combine_tables()

        self.bbc_prep_copy_files()

        # if sync-FITOPT000 option, change output for 1st iteration
        self.change_outdir_iter1()

        logging.info("")
        #sys.exit(f"\n xxx DEBUG DIE xxx \n")
        # end submit_prepare_driver

    def bbc_read_input_file(self):

        # read and store input file contents AFTER yaml block.
        # For each abc=xyz, store dictionary element.

        CONFIG     = self.config_yaml['CONFIG']
        input_file = self.config_yaml['args'].input_file
        FOUND_END_YAML = False
        contents_local = [] 
        input_file_dict = { } 
        with open(input_file,"r") as f:
            for line in f:
                word_list = (line.rstrip("\n")).split()
                if len(word_list) == 0 : continue
                word = word_list[0]
                if word    == '#END_YAML' : FOUND_END_YAML = True
                if word[0] == '#'  : continue
                if FOUND_END_YAML is False : continue 
                contents_local.append(word)
                if '=' in word :
                    jeq = word.index("=")
                    key = word[0:jeq] ; val = word[jeq+1:]
                    input_file_dict.update( {key:val} )
                    
        # - - - - 
        self.config_prep['input_file_dict'] = input_file_dict

        # end bbc_read_input_file

    def bbc_prep_version_list(self):
        # read/store list of versions for each INPDIR+.
        # Make sure it exists, along with MERGE.LOG and SUBMIT.INFO
        # Store list of versions, and list of FITOPTs for each version.

        CONFIG          = self.config_yaml['CONFIG']
        input_file      = self.config_yaml['args'].input_file 
        devel_flag      = self.config_yaml['args'].devel_flag
        IS_FITOPT_MAP   = BLOCKNAME_FITOPT_MAP in self.config_yaml
        msgerr = []
        
        key = 'STRING_VERSION_IGNORE'
        if key in CONFIG :
            string_version_ignore = CONFIG[key].split()
        else :
            string_version_ignore = []

        # - - - -
        key = 'INPDIR+'
        if key not in CONFIG :
            msgerr.append(f"Missing require key = {key} under CONFIG block")
            msgerr.append(f"Check {input_file}")
            self.log_assert(False,msgerr)

        CONFIG_INPDIR = CONFIG[key]
        n_inpdir      = len(CONFIG_INPDIR)

        inpdir_list         = [ ]
        survey_list         = [ ]    # for each inpdir
        inpdir_list_orig    = [ ]  # before expandvar
        version_list2d      = [ ] * n_inpdir  # vs. inpdir, iver
        fitopt_table_list2d = [ ] * n_inpdir
        # xxx fitopt_num_list     = [ ]
        n_fitopt_list       = [ ]
        n_version_list      = [ ]
        sync_evt_list       = [ ] 
        idir = 0; 

        config_inpdir_list, config_label_list = \
                    self.get_inpdir_list(CONFIG_INPDIR)

        if config_inpdir_list is None:
            global USE_INPDIR ; USE_INPDIR = False
            self.bbc_prep_noINPDIR()
            return

        # - - - - - 
        for path_orig in config_inpdir_list: 
            logging.info(f"  Prepare INPDIR {path_orig}")
            path_expand        = os.path.expandvars(path_orig)
            MERGE_LOG_PATHFILE = f"{path_expand}/{MERGE_LOG_FILE}"
            INFO_PATHFILE      = f"{path_expand}/{SUBMIT_INFO_FILE}"
            DONE_PATHFILE      = f"{path_expand}/{DEFAULT_DONE_FILE}"

            # check that required files exist
            msgerr = [f"Missing required {DEFAULT_DONE_FILE} file in", 
                      f"{path_orig}" ] 
            self.check_file_exists(DONE_PATHFILE,msgerr)

            msgerr = [f"Missing required {MERGE_LOG_FILE} file in", 
                      f"{path_orig}" ] 
            self.check_file_exists(MERGE_LOG_PATHFILE,msgerr)

            msgerr = [f"Missing required {SUBMIT_INFO_FILE} file in", 
                      f"{path_orig}" ] 
            self.check_file_exists(INFO_PATHFILE,msgerr)

            #  make sure DONE stamp exists with SUCCESS
            with open(DONE_PATHFILE,'r') as f :
                word = f.readlines();  word = word[0].rstrip("\n")
                if word != STRING_SUCCESS :
                    msgerr = []
                    msgerr.append(f"Expecting {STRING_SUCCESS} string written in ")
                    msgerr.append(f"   {DONE_PATHFILE}")
                    msgerr.append(f"but found {word} instead.")
                    msgerr.append(f"BBC cannot process FAILED LCFIT output.")
                    self.log_assert(False,msgerr)

            # read MERGE LOG file from LCFIT job
            MERGE_INFO,comment_lines = util.read_merge_file(MERGE_LOG_PATHFILE)
            row_list = MERGE_INFO[TABLE_MERGE]
            version_list = []
            msgerr       = [] ; nverr=0
            for row in row_list :
                version = row[COLNUM_FIT_MERGE_VERSION] 

                # check that verson dir exists 
                path_version = f"{path_expand}/{version}"
                if not os.path.exists(path_version):
                    msgerr.append(f" ERROR: found {version} in MERGE.LOG, ")
                    msgerr.append(f"\t but cannot find {path_version}")
                    nverr += 1

                ignore  = any(s in version for s in string_version_ignore)
                if ignore: continue
                if version in version_list : continue
                version_list.append(version)

            if nverr > 0:
                self.log_assert(False,msgerr)                

            survey = MERGE_INFO['SURVEY']
            survey_list.append(survey)

            # read LC FITOPT table from FIT job's submit info file
            fit_info_yaml  = util.extract_yaml(INFO_PATHFILE, None, None)
            fitopt_table   = fit_info_yaml['FITOPT_LIST']

            sync_evt       = False  # back-compatible if key isn't there
            KEY_SYNC_EVT = ""
            for key in KEYLIST_SYNC_EVT :  # check both key name options
                if key in fit_info_yaml:
                    sync_evt       = fit_info_yaml[key] > 0
                    KEY_SYNC_EVT   = key

            if devel_flag == -20: sync_evt = False # disable event sync

            # - - - 
            n_fitopt       = len(fitopt_table)

            # udpates lists vs. idir
            n_fitopt_list.append(n_fitopt)
            fitopt_table_list2d.append(fitopt_table) 
            inpdir_list.append(path_expand)
            inpdir_list_orig.append(path_orig)
            version_list2d.append(version_list)
            n_version_list.append(len(version_list))
            sync_evt_list.append(sync_evt)

            idir += 1
            #print(f" xxx ------------------------------------------")
            #print(f" xxx version_list = {version_list} \n xxx in {path} ") 
            #print(f" xxx fitopt_list({n_fitopt}) = {fitopt_table}")

        # - - - - -
        # abort if n_fitopt is different for any INPDIR
        SAME = len(set(n_fitopt_list)) == 1
        if not SAME and not IS_FITOPT_MAP :
            msgerr = []
            msgerr.append(f"Mis-match number of FITOPT; "\
                          f"n_fitopt = {n_fitopt_list} for") 
            for path in inpdir_list_orig :
                msgerr.append(f"\t {path}")
            self.log_assert(False,msgerr)

        # - - - - - -
        # abort if sync_evt is different for any INPDIR
        SAME = len(set(sync_evt_list)) == 1
        if not SAME :
            msgerr = []
            msgerr.append(f"Mis-match for {KEY_SYNC_EVT} flag; ")
            for path,sync in zip(inpdir_list_orig,sync_evt_list) :
                msgerr.append(f" {KEY_SYNC_EVT}={sync} for {path}")
            msgerr.append(f"All {KEY_SYNC_EVT} must be the same.")
            self.log_assert(False,msgerr)

        # - - - - - - - - - -
        # store the goodies
        self.config_prep['n_inpdir']          = n_inpdir
        self.config_prep['inpdir_list']       = inpdir_list
        self.config_prep['survey_list']       = survey_list
        self.config_prep['label_inpdir_list'] = config_label_list
        self.config_prep['version_list2d']    = version_list2d  # vs.idir,iver
        self.config_prep['n_version_list']      = n_version_list
        self.config_prep['n_fitopt_inplist']    = n_fitopt_list
        self.config_prep['fitopt_table_list2d'] = fitopt_table_list2d
        self.config_prep['sync_evt_list']       = sync_evt_list

        # Jun 27 2022; for sync_evt, force the kill_on_fail option
        #    since the output is useless if any task fails.
        sync_evt = sync_evt_list[0]
        if sync_evt:
            self.config_yaml['args'].kill_on_fail = True

        return;

        # end bbc_prep_version_list

    def bbc_prep_noINPDIR(self):

        logging.info(f"\n\t *** WARNING: INPDIR ignore --> " \
                     f"use datafile argument *** ")
        n_inpdir = 1
        self.config_prep['n_inpdir']            = n_inpdir
        self.config_prep['inpdir_list']         = None
        self.config_prep['survey_list']         = [ 'UNKNOWN' ]
        self.config_prep['version_list2d']      = [] * n_inpdir
        self.config_prep['fitopt_table_list2d'] = [] * n_inpdir
        self.config_prep['fitopt_num_outlist_map'] = []
        self.config_prep['sync_evt_list']          = [ False ]
        return;


    def get_inpdir_list(self,CONFIG_INPDIR):

        # for input CONFIG_INPDIR (arg of INPDIR+ key), return
        # inpdir_list, label_list.
        # CONFIG_INPDIR can either be an unlabelled list,
        #   INPDIR+:
        #     - path1
        #     - path2
        #
        # or a dictionary with labels,
        #   INPDIR+:
        #     - SURVEY1: path1
        #     - SURVEY2: path2
        
        is_list     = isinstance(CONFIG_INPDIR, list)
        n_inpdir    = len(CONFIG_INPDIR)

        # Mar 8 2021: check option to ignore INPDIR+ and use datafile= arg
        if CONFIG_INPDIR == 'None' or CONFIG_INPDIR == 'Ignore' :
            return None, None
            
        # - - - -
        if is_list :
            inpdir_list = CONFIG_INPDIR
            label_list  = [ None ] * n_inpdir
        else :
            inpdir_list = [ ]
            label_list  = [ ] 
            for label in CONFIG_INPDIR:
                path = CONFIG_INPDIR[label]
                inpdir_list.append(path)
                label_list.append(label)

        return inpdir_list, label_list

        # end get_inpdir_list

    def bbc_prep_version_match(self):
        # using input IGNORE_STRING to figure out which version in each
        # inpdir to combine with other INPDIRs. Do not assume versions
        # are in same order in each INPDIR because some INPDIRs may
        # have extra test versions that are not relevant.
        # Beware, nasty logic !

        if not USE_INPDIR : 
            self.config_prep['n_version_out']      = 1
            self.config_prep['version_out_list']   = [ 'datafile_arg' ]
            return 

        CONFIG           = self.config_yaml['CONFIG']
        n_inpdir         = self.config_prep['n_inpdir']
        inpdir_list      = self.config_prep['inpdir_list']
        n_version_list   = self.config_prep['n_version_list']
        version_list2d   = self.config_prep['version_list2d']
        msgerr = []
        key    = 'STRINGMATCH_IGNORE'
        
        # if STRINGMATCH is not defined, then there must be
        # 1 and only one version in each inpdir ... if not, abort.
        key    = 'STRINGMATCH_IGNORE'
        if key in CONFIG :
            stringmatch_ignore = CONFIG[key].split()
        else:
            stringmatch_ignore = [ 'IGNORE' ]

        # - - - - 
        if stringmatch_ignore[0] == 'IGNORE' :
            self.bbc_prep_1version_match()
            return

        # - - - - -
        logging.info(f" STRINGMATCH_IGNORE = {stringmatch_ignore} \n")

        # start by removing the stringmatch_ignore from every version
        # Versions that have nothing replaced are tossed.

        version_orig_list2d = [] * n_inpdir  # version_list minus unmatched versions
        version_out_list2d  = [] * n_inpdir 
        n_version_out_list  = []  
        isort_list2d        = [] * n_inpdir 

        for idir in range(0,n_inpdir) :
            n_version = n_version_list[idir]
            #print(f" xxx ------ idir={idir} ------------ ")
            version_out_list  = []
            version_orig_list = []
            for iver in range(0,n_version) :
                version_orig = version_list2d[idir][iver]

                version_out  = version_orig
                for str_ignore in stringmatch_ignore :
                    version_out = version_out.replace(str_ignore,"")

                if version_out != version_orig :
                    version_out_list.append(version_out)
                    version_orig_list.append(version_orig)
                    #print(f" xxx version = {version_orig} -> {version_out} ")

            version_orig_list2d.append(version_orig_list)
            version_out_list2d.append(sorted(version_out_list))
            n_version_out_list.append(len(version_out_list))

            # store sort-map in isort_list2d
            x_list = sorted((e,i) for i,e in enumerate(version_out_list))
            isort_list = []
            for v,isort in x_list :
                isort_list.append(isort)
                #print(f"\t xxx   v = {v}  isort = {isort} ")
            isort_list2d.append(isort_list)
            # end idir loop over INPDIR+

        # - - - - - 
        # make sure that number of string-replaced version_out is the same 
        # in each inpdir, otherwise it's hopeless -> abort.
        same = len(set(n_version_out_list)) == 1 
        if not same :
            msgerr = []
            msgerr.append(f"Problem applying with STRINGMATCH_IGNORE.")
            msgerr.append(f"n_version_out_list = {n_version_out_list} ")
            msgerr.append(f"has different number of versions in each INPDIR.")
            self.log_assert(False,msgerr)
            
        # create sorted list of versions to combine
        n_version_out      = n_version_out_list[0] # Nvers with string replace
        version_orig_sort_list2d = \
            [['' for i in range(n_version_out)] for j in range(n_inpdir)]
        version_out_sort_list2d = \
            [['' for i in range(n_version_out)] for j in range(n_inpdir)]
        version_out_list = [] * n_version_out  # 1D array of output versions

        for iver in range(0,n_version_out):
            v_out  = version_out_list2d[0][iver]
            version_out_list.append(v_out)
            #print(f" xxx --------------- iver={iver} ------------------ ")
            for idir in range(0,n_inpdir):
                isort  = isort_list2d[idir][iver]
                v_orig = version_orig_list2d[idir][isort]
                v_out  = version_out_list2d[idir][iver]
                version_orig_sort_list2d[idir][isort] = v_orig
                version_out_sort_list2d[idir][isort]  = v_out
                #print(f" xxx idir,iver={idir},{iver}; {v_orig}->{v_out}")

                valid_v_out = len(v_out)>0
                msgerr = []
                msgerr.append(f"idir,iver={idir},{iver}; {v_orig}->'{v_out}'")
                msgerr.append(f"Check STRINGMATCH_IGNORE")
                self.log_assert(valid_v_out, msgerr)

        self.config_prep['n_version_out']            = n_version_out
        self.config_prep['version_out_list']         = version_out_list
        self.config_prep['version_orig_sort_list2d'] = version_orig_sort_list2d
        self.config_prep['version_out_sort_list2d']  = version_out_sort_list2d

        # end bbc_prep_version_match

    def bbc_prep_1version_match(self):

        # STRINGMATCH_IGNORE is not set (i.e., ignored), so try auto-match.
        # Match is trivial for sims generated by RANSEED_REPEAT, but tricky
        # for RANSEED_CHANGE that has multiple version with -nnnn extensions.

        key              = 'STRINGMATCH_IGNORE'
        n_version_list   = self.config_prep['n_version_list']
        version_list2d   = self.config_prep['version_list2d'] # [idir][iver]
        inpdir_list      = self.config_prep['inpdir_list']
        n_inpdir         = self.config_prep['n_inpdir']

        logging.info(f"  {key} = IGNORE -> " \
                     f"Attempt auto-match version(s) in each INPDIR")

        # - - - - - - - - - - -
        # check easy case with 1 and only 1 version per INPDIR
        all_one = len(set(n_version_list)) == 1 and n_version_list[0] == 1
        if all_one :
            logging.info(f"  Auto-match success: " \
                         f"found one and only one version per INPDIR")
            self.config_prep['n_version_out']            = 1
            self.config_prep['version_orig_sort_list2d'] = version_list2d
            self.config_prep['version_out_sort_list2d']  = \
                            [[SUBDIR_OUTPUT_ONE_VERSION]]
            self.config_prep['version_out_list'] = \
                            [SUBDIR_OUTPUT_ONE_VERSION]
            return

        # Mar 10 2021: check case with only one INPDIR (with multiple versions)
        if n_inpdir == 1 :        
            logging.info(f"  Auto-match success: only one INPDIR")
            self.config_prep['n_version_out']            = n_version_list[0]
            self.config_prep['version_orig_sort_list2d'] = version_list2d
            self.config_prep['version_out_sort_list2d']  = version_list2d
            self.config_prep['version_out_list']         = version_list2d[0]
            return
            
        # - - - - - - - - 
        # Tricky case: there are multiple verions, so check for RANSEED_CHANGE 
        # sim that has suffix index per version. 
        # Example:  RANSEED_CHANGE: 3 123345 results in
        #   TEST_DES-001
        #   TEST_DES-002
        #   TEST_DES-003
        # and similarly with DES replaced by LOWZ. For this case, there is 
        # really just one simulated version with many "splitsim" versions 
        # that can be matched by the splitsim index. Note that splitsim is 
        # different than splitran used elsewhere.

        n_version        = n_version_list[0]
        all_one_splitsim = True
        version_out_sort_list2d = \
            [['' for i in range(n_version)] for j in range(n_inpdir)]
        version_out_list = [] 

        for iver in range(0,n_version):
            suffix_expect = self.suffix_splitran(n_version,iver+1) # e.g. -0001
            len_suffix    = len(suffix_expect)
            v_out         = SUBDIR_OUTPUT_ONE_VERSION + suffix_expect
            version_out_list.append(v_out)
            for idir in range(0,n_inpdir):
                version_out_sort_list2d[idir][iver] = v_out
                version = version_list2d[idir][iver]
                match   = version[-len_suffix:] == suffix_expect
                if not match:   all_one_splitsim = False

        if all_one_splitsim :
            logging.info(f"  Auto-match success: one split-sim version " \
                         f"(from RANSEED_CHANGE) per INPDIR")
            self.config_prep['n_version_out']            = n_version
            self.config_prep['version_orig_sort_list2d'] = version_list2d
            self.config_prep['version_out_sort_list2d']  = version_out_sort_list2d
            self.config_prep['version_out_list'] = version_out_list
            return
                
        # - - - - - 
        # if we get here, there is no way to auto-compute the matching
        # and therefore abort with message saying that STRINGMATCH_IGNORE 
        # is needed.
        CONFIG     = self.config_yaml['CONFIG']
        input_file = self.config_yaml['args'].input_file  # for msgerr
        msg        = []
        msg.append(f"Cannot auto-match data versions and therefore {key} key")
        msg.append(f"is needed in {input_file}")
        msg.append(f"\t(for details: submit_batch_jobs.sh -H BBC)" )
        msg.append(f"For auto-match, only 1 VERSION per INPDIR is allowed,")
        msg.append(f"or 1 sim version using RANSEED_CHANGE.")
        msg.append(f"n_version_list = {n_version_list} for ")
        msg += inpdir_list
        self.log_assert(False,msg)

        # end bbc_prep_1version_match

    def bbc_prep_fitopt_outlist(self):

        # Determine output fitopts in 1 of 2 ways:
        #  1) no FITOPT_MAP in the input file -->
        #     output FITOPT list = input FITOPT list
        #  2) use FITOPT_MAP (likely created by Pippin)
        #

        if not USE_INPDIR : 
            self.config_prep['n_fitopt']            = 1
            self.config_prep['fitopt_num_outlist']  = [ 'FITOPT000' ]
            self.config_prep['FITOPT_OUT_LIST']     = [ ]
            return 

        IS_FITOPT_MAP = BLOCKNAME_FITOPT_MAP in self.config_yaml
        n_inpdir            = self.config_prep['n_inpdir']
        n_fitopt_inplist    = self.config_prep['n_fitopt_inplist']  
        fitopt_table_list2d = self.config_prep['fitopt_table_list2d'] #idir,ifit
        survey_inplist      = self.config_prep['survey_list'] 

        fitopt_num_outlist     = []
        fitopt_num_outlist_map = []
        msgerr = []

        if IS_FITOPT_MAP :
            FITOPT_MAP = self.config_yaml[BLOCKNAME_FITOPT_MAP]

        else :
            # there is no input FITOPT_MAP, so create a trivial map where
            # each output FITOPT maps onto the input FITOPTs with the
            # same FITOPT num-index
            n_fitopt   = n_fitopt_inplist[0] # same for all INPDIR
            FITOPT_MAP = {}
            for ifit in range(0,n_fitopt) :
                fitopt_row   = fitopt_table_list2d[0][ifit]
                fitopt_num   = fitopt_row[COLNUM_FITOPT_NUM] # e.g. FITOP004
                fitopt_list  = f"{fitopt_num} " * n_inpdir
                FITOPT_MAP[fitopt_num] = fitopt_list
                
        # - - - - -
        #print(f" xxx USE ME: FITOPT_MAP = {FITOPT_MAP}")

        # if FITOPT_MAP includes SURVEY_LIST, check that order matches
        # that of INPDIRs
        if 'SURVEY_LIST' in FITOPT_MAP :
            survey_maplist = FITOPT_MAP['SURVEY_LIST'].split()
            if survey_maplist != survey_inplist :
                msgerr.append(f"SURVEY_LIST mis-match.")
                msgerr.append(f"  SURVEY_LIST(INPDIR+)    = {survey_inplist}")
                msgerr.append(f"  SURVEY_LIST(FITOPT_MAP) = {survey_maplist}")
                msgerr.append(f"SURVEY_LIST order in INPDIR+ and FITOPT_MAP "
                              f"must match.")
                self.log_assert(False,msgerr)

        # Create outlist arrays for submit logic.
        # Count n_fitopt explicitly in case FITOPT_MAP contains
        # other keys that are not FITOPT; e.g., SURVEY_LIST
        n_fitopt = 0
        for fitopt_num in FITOPT_MAP :
            if fitopt_num[0:6] != 'FITOPT' : continue
            fitopt_num_outlist.append(fitopt_num)
             
            # for each INPDIR, load map of original FITOPT needed to copy
            # original FITRES file
            fitopt_map = FITOPT_MAP[fitopt_num].split()
            fitopt_num_outlist_map.append(fitopt_map)
            n_fitopt    += 1
            

        #print(f"\n xxx map = {fitopt_num_outlist_map} \n")

        # - - - - - - - - - 
        print(f"\n Prepare output FITOPT list:")
        for survey,n in zip(survey_inplist,n_fitopt_inplist):
            logging.info(f"   Found {n:3d} FITOPTs for SURVEY = {survey}")

        if IS_FITOPT_MAP :
            logging.info(f"   Found {BLOCKNAME_FITOPT_MAP}  ")
        else :
            logging.info(f"   Did not find {BLOCKNAME_FITOPT_MAP}  ")

        logging.info(f"   --> {n_fitopt} output FITOPTs")            

        # store *outlist arrays
        self.config_prep['n_fitopt']               = n_fitopt
        self.config_prep['fitopt_num_outlist']     = fitopt_num_outlist
        self.config_prep['fitopt_num_outlist_map'] = fitopt_num_outlist_map

        # - - - - - - - - - - - - - 
        # prepare FITOPT_OUT_LIST table for SUBMIT.INFO file
        self.make_FITOPT_OUT_LIST()

        # end bbc_prep_fitopt_outlist(self)

    def make_FITOPT_OUT_LIST(self):

        # Created May 27 2021
        # Construct output FITOPT_OUT_LIST for SUBMIT.INFO file.
        # Each table row contains:
        #    'FITOPTNUM'  'SURVEY'  'user_label'   'user_args'
        #

        n_fitopt        = self.config_prep['n_fitopt']      
        n_inpdir        = self.config_prep['n_inpdir'] 
        survey_list     = self.config_prep['survey_list'] 
        fitopt_num_list = self.config_prep['fitopt_num_outlist']    
        fitopt_num_map  = self.config_prep['fitopt_num_outlist_map']
        fitopt_table_list2d = self.config_prep['fitopt_table_list2d'] #idir,ifit

        dump_flag      = False
        FITOPT_OUT_LIST = []

        if not USE_INPDIR : 
            item_list = [ 'GLOBAL', 'FITOPT000', None, None ]
            FITOPT_OUT_LIST.append(item_list)
            self.config_prep['FITOPT_OUT_LIST'] = FITOPT_OUT_LIST
            return

        ifit_out = 0
        for fitopt_num_out in fitopt_num_list:            

            if dump_flag :
                logging.info(" xxx ---------------------------------------- ")
            # check if this FITOPT is global, or specific to one survey
            fitopt_num_inplist  = fitopt_num_map[ifit_out][0:n_inpdir]

            n_arg_none = 0 ;  n_arg_FITOPT000 = 0; n_arg_define = 0
            survey_store = 'ERROR' ;  label_store = None; arg_store = None 
            for idir in range(0,n_inpdir):          
                survey         = survey_list[idir]   
                fitopt_num_inp = fitopt_num_inplist[idir]
                ifit_inp       = int(fitopt_num_inp[6:])
                row    = fitopt_table_list2d[idir][ifit_inp]
                num    = row[COLNUM_FITOPT_NUM]  # e.g., FITOPT003
                label  = row[COLNUM_FITOPT_LABEL]
                arg    = row[COLNUM_FITOPT_ARG]
                if arg == '' : # only for FITOPT000
                    n_arg_none  += 1
                elif arg == 'FITOPT000' : # sym link back to FITOPT000
                    n_arg_FITOPT000 += 1  
                else :                    # genuine LC fit arg list
                    survey_store = survey
                    arg_store    = arg
                    label_store  = label
                    n_arg_define += 1

                if dump_flag :
                    print(f" xxx {fitopt_num_out}: idir={idir} num={num} " \
                          f"label={label} arg='{arg}'")

            # - - - - - - - - - - - - - - - - - - - - - 
            # if all args are valid, set survey_store to GLOBAL
            if n_inpdir > 1 :
                if n_arg_define == n_inpdir or ifit_out == 0 :
                    survey_store = 'GLOBAL'
                else:
                    pass # leave survey_store as is
            else:
                # no need for GLOBAL if only one survey
                survey_store = survey

            ifit_out += 1

            # construct and write yaml-compliant info list
            item_list = []
            item_list.append(fitopt_num_out)
            item_list.append(survey_store)
            item_list.append(label_store)
            item_list.append(arg_store)
            FITOPT_OUT_LIST.append(item_list)

        self.config_prep['FITOPT_OUT_LIST'] = FITOPT_OUT_LIST
        return
        # end make_FITOPT_OUT_LIST

    def bbc_prep_index_lists(self):
        # construct sparse 1D lists to loop over version and FITOPT

        CONFIG        = self.config_yaml['CONFIG']
        n_version     = self.config_prep['n_version_out']  
        n_fitopt      = self.config_prep['n_fitopt']
        n_muopt       = self.config_prep['n_muopt']
        n_splitran    = self.config_prep['n_splitran']

        n2d_index = n_version * n_fitopt
        n4d_index = n_version * n_fitopt * n_muopt * n_splitran

        # fetch matrix for which ifit x imu to use
        n_use2d,use_matrix2d,use_fitopt = self.get_matrix_FITOPTxMUOPT()

        # ---------------------------
        # create 2D index lists used to prepare inputs
        # (create subdirs, catenate input FITRES files)
        iver_list2=[]; ifit_list2=[]; 
        for iver in range(0,n_version):
            for ifit in range(0,n_fitopt):
                if not use_fitopt[ifit]: continue
                iver_list2.append(iver)
                ifit_list2.append(ifit)

        # 4D lists are for prep_JOB_INFO
        iver_list4=[]; ifit_list4=[]; imu_list4=[]; isplitran_list4=[]
        for ifit in range(0,n_fitopt):
            for imu in range(0,n_muopt):
                for iver in range(0,n_version):
                    if not use_matrix2d[ifit][imu] : continue
                    for isplitran in range(0,n_splitran):
                        iver_list4.append(iver)
                        ifit_list4.append(ifit)
                        imu_list4.append(imu)
                        isplitran_list4.append(isplitran+1) # 1-n_splitran

        self.config_prep['n2d_index']  = n2d_index
        self.config_prep['iver_list2'] = iver_list2
        self.config_prep['ifit_list2'] = ifit_list2

        self.config_prep['n4d_index']  = n4d_index
        self.config_prep['iver_list4'] = iver_list4
        self.config_prep['ifit_list4'] = ifit_list4
        self.config_prep['imu_list4']  = imu_list4
        self.config_prep['isplitran_list4']  = isplitran_list4

        self.config_prep['use_matrix2d']    = use_matrix2d
        self.config_prep['use_fitopt']      = use_fitopt
        self.config_prep['n_use_matrix2d']  = n_use2d

        # end bbc_prep_index_lists

    def get_matrix_FITOPTxMUOPT(self):

        # return use_matrix2d (ifit x imu) for which matrix elements
        # are used (True) or ignored (False).
        # use_matrix1d is a 1D array of which ifit are used.
        #
        # See more info with submit_batch_jobs.sh -H BBC
        # Dec 17 2020: update to process list of FITOPTxMUOPT 

        n_fitopt        = self.config_prep['n_fitopt']
        n_muopt         = self.config_prep['n_muopt']
        CONFIG          = self.config_yaml['CONFIG']
        ignore_fitopt   = self.config_yaml['args'].ignore_fitopt
        ignore_muopt    = self.config_yaml['args'].ignore_muopt
        msgerr        = []
        ALL_FLAG      = False
        or_char       = '+'
        and_char      = '&'
        bool_logic_list  = []
        ifit_logic_list  = []
        imu_logic_list   = []
        dump_flag_matrix = False
        ifit_logic = -9;  imu_logic = -9

        # check if FITOPTxMUOPT is specified, and whether it's 
        # one value (str) or a list. If one value, convert to list.
        FITOPTxMUOPT_LIST = []
        if KEY_FITOPTxMUOPT in CONFIG :
            FITOPTxMUOPT_LIST = CONFIG[KEY_FITOPTxMUOPT]
            if isinstance(FITOPTxMUOPT_LIST,str): 
                FITOPTxMUOPT_LIST =  [ FITOPTxMUOPT_LIST ]

        if 'DUMP' in FITOPTxMUOPT_LIST :
            dump_flag_matrix  = True
            FITOPTxMUOPT_LIST.remove('DUMP')

        if len(FITOPTxMUOPT_LIST) == 0 :
            ALL_FLAG     = True
            CONFIG[KEY_FITOPTxMUOPT] = "ALL"

        #print(f" xxx FITOPTxMUOPT_LIST = {FITOPTxMUOPT_LIST}")
        #print(f" xxx ALL_FLAG = {ALL_FLAG} ")
        #print(f" xxx dump_flag_matrix = {dump_flag_matrix} ")
        # - - - - - - 

        for FITOPTxMUOPT in FITOPTxMUOPT_LIST :

            if ignore_fitopt:
                msgerr.append(f"Cannot mix {KEY_FITOPTxMUOPT} key with " \
                              f"command line arg --ignore_fitopt")
                self.log_assert(False,msgerr)

            if ignore_muopt:
                msgerr.append(f"Cannot mix {KEY_FITOPTxMUOPT} key with " \
                              f"command line arg --ignore_muopt")
                self.log_assert(False,msgerr)

            if or_char in FITOPTxMUOPT: 
                bool_logic = or_char ;  bool_string = "or"
            elif and_char in FITOPTxMUOPT:
                bool_logic = and_char ; bool_string = "and"
            else:
                msgerr.append(f"Invalid {KEY_FITOPTxMUOPT}: {FITOPTxMUOPT}")
                msgerr.append(f"Must have {or_char} or {and_char} symbol ")
                self.log_assert(False,msgerr)
                
            j_bool       = FITOPTxMUOPT.index(bool_logic)    
            ifit_logic   = int(FITOPTxMUOPT[0:j_bool])   # fitopt number
            imu_logic    = int(FITOPTxMUOPT[j_bool+1:])  # muopt number

            msg = f"  {KEY_FITOPTxMUOPT} logic: process " \
                  f"FITOPT={ifit_logic} {bool_string} MUOPT={imu_logic} "
            logging.info(msg)

            bool_logic_list.append(bool_logic)
            ifit_logic_list.append(ifit_logic)
            imu_logic_list.append(imu_logic)


        # - - - - - - - - - - 
        n_use2d = 0
        use_matrix2d = \
            [[ False for i in range(n_muopt)] for j in range(n_fitopt)]
        use_matrix1d = [ False ]  * n_fitopt
        
        for ifit in range(0,n_fitopt):
            for imu in range(0,n_muopt):

                if ALL_FLAG:
                    use2d = True
                else:
                    use2d = False

                for bool_logic, ifit_logic, imu_logic in \
                    zip(bool_logic_list, ifit_logic_list, imu_logic_list):

                    use_ifit = (ifit == ifit_logic )
                    use_imu  = (imu  == imu_logic  )
                    if bool_logic == or_char :
                        if use_ifit or use_imu: use2d = True
                    elif bool_logic == and_char :
                        if use_ifit and use_imu : use2d = True

                if ignore_fitopt : use2d = (ifit == 0)
                if ignore_muopt  : use2d = (use2d and imu==0)

                #print(f" xxx check ifit,imu = {ifit}({ifit_logic}, "
                #      f"{imu}({imu_logic})  " \
                #      f" use[fit,mu,2d]={use_ifit}, {use_imu}, {use2d}")

                if use2d :
                    use_matrix2d[ifit][imu] = True
                    use_matrix1d[ifit]      = True
                    n_use2d += 1

        #print(f"  xxx use_matrix2d = \n\t {use_matrix2d} ")
        #print(f"  xxx use_matrix1d = \n\t {use_matrix1d} ")
        #sys.exit(f"\n xxx DEBUG DIE xxx\n")
        
        if dump_flag_matrix :    
            self.dump_matrix_FITOPTxMUOPT(use_matrix2d)

        return n_use2d, use_matrix2d, use_matrix1d

        # end get_matrix_FITOPTxMUOPT

    def dump_matrix_FITOPTxMUOPT(self,matrix):
        
        n_fitopt        = self.config_prep['n_fitopt']
        n_muopt         = self.config_prep['n_muopt']

        logging.info(f"\n Dump FITOPTxMUOPT matrix: ")
        for ifit in range(0,n_fitopt):
            line = f"   FITOPT{ifit:03d}:  "
            for imu in range(0,n_muopt):
                USE = "F"
                if matrix[ifit][imu] : USE = "T"
                line += f"{USE} "
            logging.info(f"{line}")

        # end dump_matrix_FITOPTxMUOPT

    def bbc_prep_mkdir(self):

        # create version-output directory(s).
        # For NSPLITRAN, split outdir into outdir per splitran.

        output_dir       = self.config_prep['output_dir']  
        version_out_list = self.config_prep['version_out_list']
        n_splitran       = self.config_prep['n_splitran']
        USE_SPLITRAN     = n_splitran > 1

        logging.info("\n  Create BBC output directories:")
        for v_out in version_out_list :
            for i in range(0,n_splitran):
                isplitran = i + 1  # 1 to n_splitran
                v_dir     = v_out + self.suffix_splitran(n_splitran,isplitran)
                V_DIR     = f"{output_dir}/{v_dir}"
                logging.info(f"    Create BBC output dir {v_dir} ")
                os.mkdir(V_DIR)

        # end bbc_prep_outdirs

    def bbc_prep_combine_tables(self):

        # Catenate FITRES files from INPDIR+ so that each copied
        # FITRES file includes multiple surveys
        # For NSPLITRAN, copy only to first split-dir to avoid
        # duplicate copies of input FITRES files.

        if not USE_INPDIR: return

        output_dir         = self.config_prep['output_dir']  
        n_inpdir           = self.config_prep['n_inpdir']  
        v_out_list         = self.config_prep['version_out_sort_list2d']
        iver_list2         = self.config_prep['iver_list2'] 
        ifit_list2         = self.config_prep['ifit_list2']
        fitopt_num_outlist = self.config_prep['fitopt_num_outlist']
        n_splitran         = self.config_prep['n_splitran']
        USE_SPLITRAN       = n_splitran > 1

        cat_file_log   = f"{output_dir}/cat_FITRES_SALT2mu.LOG"

        logging.info("\n  Prepare input FITRES files")
        iver_last = -9

        for iver,ifit in zip(iver_list2, ifit_list2):

            idir0 = 0  # some things just need first INPDIR index

            # get output dir name
            v_dir   = v_out_list[idir0][iver]
            v_dir  += self.suffix_splitran(n_splitran,1)
            V_DIR   = f"{output_dir}/{v_dir}"

            cat_list   = self.make_cat_fitres_list(iver,ifit)

            #print(" xxx ---------------------------------- " )
            #print(f" xxx iver={iver} ifit={ifit} :")
            #print(f" xxx cat_list = {cat_list} ")

            # execute the FITRES catenation
            fitopt_num     = fitopt_num_outlist[ifit]
            ff             = f"{fitopt_num}.{SUFFIX_FITRES}"
            input_ff       = "INPUT_" + ff
            cat_file_out   = f"{V_DIR}/{input_ff}"
            nrow = self.exec_cat_fitres(cat_list, cat_file_out, cat_file_log)

            if iver != iver_last : logging.info(f"    {v_dir}: ")
            iver_last = iver

            logging.info(f"\t Catenate {n_inpdir} {ff} files"\
                         f" -> {nrow} events ")

        # - - - - - 
        logging.info("   gzip the catenated FITRES files.")
        cmd_gzip = f"cd {output_dir}; gzip */INPUT_FITOPT*.{SUFFIX_FITRES}"
        os.system(cmd_gzip)

        # remove cat log file
        rm_log = f"cd {output_dir}; rm {cat_file_log}"
        os.system(rm_log)

        # end bbc_prep_combine_tables
    
        
    def exec_cat_fitres(self,cat_list, cat_file_out, cat_file_log):

        # prepare & execute catenate command for this cat_list 
        # (comma-sep list of files) into cat_file_out.
        # Use the "SALT2mu.exe cat_only" option to handle different
        # columns in each FITRES file. While .gz extensions are
        # not included, SALT2mu handles files with or without
        # .gz extensions.
        #
        # function returns number of rows in catenated file

        snana_dir    = self.config_yaml['args'].snana_dir
        if snana_dir is None:
            code_name = PROGRAM_NAME_BBC
        else:
            code_name = f"{snana_dir}/bin/{PROGRAM_NAME_BBC}"

        cmd_cat = f"{code_name}  " \
                  f"cat_only  "    \
                  f"datafile={cat_list}  " \
                  f"append_varname_missing='PROB*,zPRIOR*'  " \
                  f"catfile_out={cat_file_out}  " \
                  f" &> {cat_file_log}"

        #print(f" xxx command to cat fitres files, \n {cmd_cat} \n")
        os.system(cmd_cat)

        # check number of rows
        nrow = util.nrow_table_TEXT(cat_file_out, "SN:")

        return nrow

        # end exec_cat_fitres

    def make_cat_fitres_list(self, iver, ifit ):
        
        # Use input indices for version (iver) and fitopt (ifit)
        # to construct comma-sep list of fitres files over INPDIRs.
        # This list is used to catenate FITRES files from different
        # input directories

        n_inpdir        = self.config_prep['n_inpdir']  
        v_orig_list     = self.config_prep['version_orig_sort_list2d']
        inpdir_list     = self.config_prep['inpdir_list']
        fitopt_num_list = self.config_prep['fitopt_num_outlist'] 
        fitopt_num_map  = self.config_prep['fitopt_num_outlist_map'] 

        cat_list = ''
        fitopt_num_out   = fitopt_num_list[ifit]
        nmissing         = 0 ; msgerr=[]

        for idir in range(0,n_inpdir):  
            inpdir       = inpdir_list[idir]
            v_orig       = v_orig_list[idir][iver]
            fitopt_num_inp  = fitopt_num_map[ifit][idir] 

            ff    = f"{v_orig}/{fitopt_num_inp}.{SUFFIX_FITRES}"
            FF    = f"{inpdir}/{ff}"
            FFgz  = f"{FF}.gz"
            cat_list += f"{FF},"

            EXIST_FF = os.path.isfile(FF) or os.path.isfile(FFgz)
            if not EXIST_FF:
                logging.info(f" ERROR: missing {FF}")
                nmissing += 1

        cat_list = cat_list.rstrip(",")  # remove trailing comma

        if nmissing > 0 :
            msgerr.append(f"missing {nmissing} input {SUFFIX_FITRES} files")
            self.log_assert(False,msgerr)

        return cat_list

        # end make_cat_fitres_list

    def bbc_prep_muopt_list(self):
    
        # Jan 24 2021: refactor using prep_jobopt_list util

        CONFIG   = self.config_yaml['CONFIG']
        key      = MUOPT_STRING
        if key in CONFIG  :
            muopt_rows = CONFIG[key]
        else:
            muopt_rows = []

        # - - - - - -
        muopt_dict = util.prep_jobopt_list(muopt_rows, MUOPT_STRING, 1, None)

        n_muopt          = muopt_dict['n_jobopt']
        muopt_arg_list   = muopt_dict['jobopt_arg_list']
        muopt_num_list   = muopt_dict['jobopt_num_list']
        muopt_label_list = muopt_dict['jobopt_label_list']
                
        logging.info(f" Store {n_muopt-1} BBC options from MUOPT keys")

        self.config_prep['n_muopt']          = n_muopt
        self.config_prep['muopt_arg_list']   = muopt_arg_list
        self.config_prep['muopt_num_list']   = muopt_num_list
        self.config_prep['muopt_label_list'] = muopt_label_list

        # end bbc_prep_muopt_list


    def bbc_prep_splitran(self) :

        CONFIG        = self.config_yaml['CONFIG']

        # check NSPLITRAN option to fit random sub-samples.
        # Running on sim data, this feature is useful to measure
        # rms on cosmo params, and compare with fitted error.
        key_nsplitran = 'NSPLITRAN'
        n_splitran    = 1
        if key_nsplitran in CONFIG : 
            # check YAML block input, NSPLITRAN:  <nsplitran> 
            n_splitran = CONFIG[key_nsplitran]
        else :
            # check alternate C-code input, NSPLITRAN=<nsplitran>
            input_file_dict =  self.config_prep['input_file_dict']
            if key_nsplitran in input_file_dict :
                n_splitran = int(input_file_dict[key_nsplitran])

        # - - - 
        self.config_prep['n_splitran']  = n_splitran

        # end bbc_prep_splitran


    def bbc_prep_copy_files(self):
        input_file    = self.config_yaml['args'].input_file 
        script_dir    = self.config_prep['script_dir']
        shutil.copy(input_file,script_dir)

        # end bbc_prep_copy_files

    def change_outdir_iter1(self):

        output_dir   = self.config_prep['output_dir']  
        sync_evt     = self.config_prep['sync_evt_list'][0]
        iter2        = self.config_yaml['args'].iter2
        iter1        = not iter2

        if not sync_evt : return

        # set submit_iter for global control
        if iter1: 
            self.config_prep['submit_iter'] = 1
        else:
            self.config_prep['submit_iter'] = 2 
            return

        # change config arguments
        override_output_dir = f"{output_dir}{OUTDIR_ITER1_SUFFIX}"
        override_script_dir = f"{override_output_dir}/{SUBDIR_SCRIPTS_BBC}"

        self.config_prep['script_dir']  = override_script_dir
        self.config_prep['output_dir']  = override_output_dir
        self.config_yaml['args'].outdir = override_output_dir

        # move current outdir to override_outdir
        output_subdir          = os.path.basename(output_dir)
        override_output_subdir = os.path.basename(override_output_dir)

        logging.info(f"\n\t move  {output_subdir} -> {override_output_subdir}")

        # remove target dir if it already  exists
        if  os.path.exists(override_output_dir) : 
            shutil.rmtree(override_output_dir)

        # do the move
        shutil.move(output_subdir,override_output_subdir)

        # end change_outdir_iter1

    # ===================================================
    def write_command_file(self, icpu, f):
        # For this icpu, write full set of sim commands to  
        # already-opened command file with pointer f.
        # Function returns number of jobs for this cpu 

        input_file      = self.config_yaml['args'].input_file 
        n_version       = self.config_prep['n_version_out']  
        n_fitopt        = self.config_prep['n_fitopt']
        n_muopt         = self.config_prep['n_muopt']
        n_splitran      = self.config_prep['n_splitran']
        muopt_arg_list  = self.config_prep['muopt_arg_list']
        n_core          = self.config_prep['n_core']

        n4d_index        = self.config_prep['n4d_index']
        iver_list4       = self.config_prep['iver_list4'] 
        ifit_list4       = self.config_prep['ifit_list4']
        imu_list4        = self.config_prep['imu_list4']
        isplitran_list4  = self.config_prep['isplitran_list4']

        n_use_matrix2d = self.config_prep['n_use_matrix2d'] # FITOPT x MUOPT

        CONFIG   = self.config_yaml['CONFIG']
        use_wfit = 'WFITMUDIF_OPT' in CONFIG  # check follow-up job after bbc

        #n_job_tot   = n_version * n_fitopt * n_muopt * n_splitran
        n_job_tot   = n_version * n_use_matrix2d * n_splitran
        n_job_split = 1     # cannot break up BBC job as with sim or fit
        n_job_cpu   = 0

        self.config_prep['n_job_split'] = n_job_split
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot
        self.config_prep['use_wfit']    = use_wfit

        n_job_local = 0

        for iver,ifit,imu,isplitran in \
            zip(iver_list4,ifit_list4,imu_list4,isplitran_list4):

            n_job_local += 1
            index_dict = \
                { 'iver':iver, 'ifit':ifit, 'imu':imu, 'icpu':icpu,
                  'isplitran': isplitran }

            if ( (n_job_local-1) % n_core ) == icpu :

                n_job_cpu += 1
                job_info_bbc   = self.prep_JOB_INFO_bbc(index_dict)
                util.write_job_info(f, job_info_bbc, icpu)

                if use_wfit :
                    job_info_wfit  = self.prep_JOB_INFO_wfit(index_dict)
                    util.write_job_info(f, job_info_wfit, icpu)
    
                merge_force = self.set_merge_force_bbc(ifit)
                job_info_merge = \
                    self.prep_JOB_INFO_merge(icpu,n_job_local,merge_force) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

        # - - - - 

        return n_job_cpu

        # end write_command_file

    def set_merge_force_bbc(self,ifit):

        # Created Apr 8 2022 by R.Kessler
        # If event-sync is set for first iteration,
        # then return merge_force=True on FITOPT000.
        # Goal is to avoid hang up when another
        # merge process sets a BUSY.

        sync_evt     = self.config_prep['sync_evt_list'][0]
        iter2        = self.config_yaml['args'].iter2
        iter1        = not iter2

        if not sync_evt: return False

        FITOPT_OUT_LIST = self.config_prep['FITOPT_OUT_LIST']
        label       = FITOPT_OUT_LIST[ifit][COLNUM_FITOPT_LABEL]
        skip_sync   = FITOPT_STRING_NOREJECT in label
        if skip_sync: return False

        if iter1 and ifit==0:
            return True
        else:
            return False

        # end set_merge_force_bbc

    def prep_JOB_INFO_bbc(self,index_dict):
        # Return JOB_INFO dictionary with 
        #   cd job_dir
        #   program.exe arg_list  > log_file
        #   touch TMP_[xxx].DONE
        #
        # Inputs
        #   index_dict = dictionary of indices for this job
        #
        # Beware that input FITRES files are not in script_dir,
        # but they are in ../[version]

        # strip off indices from input dictionary
        iver      = index_dict['iver']
        ifit      = index_dict['ifit']
        imu       = index_dict['imu'] 
        isplitran = index_dict['isplitran'] 

        #print(f" xxx iver={iver}, ifit={ifit}, imu={imu} ", \
            #flush=True)

        input_file    = self.config_yaml['args'].input_file 
        prescale      = self.config_yaml['args'].prescale
        kill_on_fail  = self.config_yaml['args'].kill_on_fail
        iter2         = self.config_yaml['args'].iter2
        iter1         = not iter2

        program      = self.config_prep['program']
        output_dir   = self.config_prep['output_dir']
        script_dir   = self.config_prep['script_dir']
        version      = self.config_prep['version_out_list'][iver]
        fitopt_num   = self.config_prep['fitopt_num_outlist'][ifit] #e.g FITOPT002
        muopt_num    = self.config_prep['muopt_num_list'][imu] # e.g MUOPT003
        muopt_arg    = self.config_prep['muopt_arg_list'][imu]
        n_splitran   = self.config_prep['n_splitran']
        USE_SPLITRAN = n_splitran > 1
        use_wfit     = self.config_prep['use_wfit']
        sync_evt     = self.config_prep['sync_evt_list'][0]

        # construct row mimicking MERGE.LOG
        
        row = [ None, version, fitopt_num, muopt_num, 0,0,0, isplitran ]
        prefix_orig, prefix_final = self.bbc_prefix("bbc", row)
        input_ff      = f"INPUT_{fitopt_num}.{SUFFIX_FITRES}"
        log_file      = f"{prefix_orig}.LOG"
        done_file     = f"{prefix_orig}.DONE"
        all_done_file = f"{output_dir}/{DEFAULT_DONE_FILE}"

        JOB_INFO = {}
        JOB_INFO['program']       = program
        JOB_INFO['input_file']    = input_file
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = f"{log_file}"
        JOB_INFO['done_file']     = f"{done_file}"
        JOB_INFO['all_done_file'] = f"{all_done_file}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
                
        # if wfit job will run, suppress DONE file here and wait for
        # wfit to finish before writing DONE files. This logic avoids
        # confusing the merge process where SALT2mu has finished but
        # wfit still runs.
        if use_wfit : JOB_INFO['done_file'] = ''

        arg_list = []
        arg_list.append(f"  prefix={prefix_orig}")

        # get data file from ../version, or for splitran,
        # get data file from first splitran directory
        if USE_INPDIR:
            version_datafile = version + self.suffix_splitran(n_splitran,1)
            arg_list.append(f"  datafile=../{version_datafile}/{input_ff}")

        arg_list.append(f"  write_yaml=1")

        if USE_SPLITRAN :
            # note that fortran-like isplitran index is used here
            arg = f"NSPLITRAN={n_splitran} JOBID_SPLITRAN={isplitran}"
            arg_list.append(f"  {arg}")

        arg_list.append(f"{muopt_arg}")     # user input

        # check command line input --fast option to prescale by 10
        # Only sim is pre-scaled; don't pre-scale data.
        if prescale > 1 :
            arg_list.append(f"prescale_simdata={prescale}")

        # - - - - - -
        # check option to use FITOPT000 events for all syst.

        if sync_evt : 
            FITOPT_OUT_LIST = self.config_prep['FITOPT_OUT_LIST']
            label       = FITOPT_OUT_LIST[ifit][COLNUM_FITOPT_LABEL]
            skip_sync   = FITOPT_STRING_NOREJECT in label
            wait_file   = None
            select_file = None

            # check logic for each iterations. 
            if skip_sync : 
                pass

            elif iter1 and ifit > 0 :
                # process events from FITOPT000_MUOPT000
                row = [ None, version, "FITOPT000", "MUOPT000", 0,0,0,  0 ]
                prefix_orig, prefix_final = self.bbc_prefix("bbc", row)
                ff_file     = f"../{version}/{prefix_final}.{SUFFIX_FITRES}"
                wait_file   = f"{ff_file}.gz"
                select_file = ff_file

            elif iter2 :
                # process ONLY the common events from output of iter1 FITOPTs
                outdir_iter1 = f"{output_dir}{OUTDIR_ITER1_SUFFIX}"
                wait_file    = f"{outdir_iter1}/{DEFAULT_DONE_FILE}"
                # xxx?? f" {STRING_SUCCESS}" # require SUCCESS in file
                select_file  = f"{outdir_iter1}/{version}/{BBC_ACCEPT_SUMMARY_FILE}"

            if wait_file is not None :
                JOB_INFO['wait_file'] = wait_file
            if select_file is not None :
                arg_list.append(f"cid_select_file={select_file}")
                arg_list.append(f"CUTWIN NONE")  # turn off CUTWIN cuts

        # - - - - - 
        JOB_INFO['arg_list'] = arg_list

        return JOB_INFO

        # end prep_JOB_INFO_bbc


    def prep_JOB_INFO_wfit(self,index_dict):
        # optional: run wfit cosmology fitter if WFITMUDIF_OPT is set in CONFIG

        iver      = index_dict['iver']
        ifit      = index_dict['ifit']
        imu       = index_dict['imu'] 
        isplitran = index_dict['isplitran'] 

        CONFIG        = self.config_yaml['CONFIG']
        output_dir    = self.config_prep['output_dir']
        script_dir    = self.config_prep['script_dir']
        version       = self.config_prep['version_out_list'][iver]
        fitopt_num    = self.config_prep['fitopt_num_outlist'][ifit] 
        muopt_num     = self.config_prep['muopt_num_list'][imu] # e.g MUOPT003

        row = [ None, version, fitopt_num, muopt_num, 0,0,0, isplitran ]
        prefix_bbc_orig,  prefix_bbc_final  = self.bbc_prefix("bbc",  row)
        prefix_wfit_orig, prefix_wfit_final = self.bbc_prefix("wfit", row)
    
        # note that the done file has the SALT2mu/BBC done stamp,
        # not a wfit done stamp.
        wfit_inp_file   = f"{prefix_bbc_orig}.{SUFFIX_M0DIF}"
        wfit_done_file  = f"{prefix_bbc_orig}.DONE"
        wfit_out_file   = f"{prefix_wfit_orig}.YAML"
        wfit_log_file   = f"{prefix_wfit_orig}.LOG"

        arg_list = []
        arg_list.append(f"-cospar_yaml {wfit_out_file} ") 
        arg_list.append(CONFIG['WFITMUDIF_OPT'])

        JOB_INFO = {}
        JOB_INFO['program']     = PROGRAM_wfit
        JOB_INFO['input_file']  = wfit_inp_file
        JOB_INFO['log_file']    = wfit_log_file
        JOB_INFO['done_file']   = wfit_done_file
        JOB_INFO['job_dir']     = ""  # same job dir as SALT2mu.exe job
        JOB_INFO['arg_list']    = arg_list

        return JOB_INFO

        # prep_JOB_INFO_wfit 

    def append_info_file(self,f):
        # append info to SUBMIT.INFO file

        vout_list         = self.config_prep['version_out_list']
        n_version         = self.config_prep['n_version_out']
        n_fitopt          = self.config_prep['n_fitopt']
        n_muopt           = self.config_prep['n_muopt']
        muopt_arg_list    = self.config_prep['muopt_arg_list']
        muopt_num_list    = self.config_prep['muopt_num_list']
        muopt_label_list  = self.config_prep['muopt_label_list']
        inpdir_list       = self.config_prep['inpdir_list']
        survey_list       = self.config_prep['survey_list']
        n_splitran        = self.config_prep['n_splitran']
        use_wfit          = self.config_prep['use_wfit']
        input_file_dict   = self.config_prep['input_file_dict']
        ignore_muopt      = self.config_yaml['args'].ignore_muopt
        ignore_fitopt     = self.config_yaml['args'].ignore_fitopt
        iter2             = self.config_yaml['args'].iter2
        sync_evt          = self.config_prep['sync_evt_list'][0]
        FITOPT_OUT_LIST   = self.config_prep['FITOPT_OUT_LIST']
        CONFIG            = self.config_yaml['CONFIG']
        FITOPTxMUOPT      = CONFIG[KEY_FITOPTxMUOPT]

        f.write(f"\n# BBC info\n")

        # beware that LOG,DONE,YAML files are not under script_dir,
        # but under ../[VERSION]
        f.write(f"JOBFILE_WILDCARD:  '*FITOPT*MUOPT*' \n")

        f.write(f"NVERSION:       {n_version}      " \
                f"# number of data VERSIONs\n")
        f.write(f"NFITOPT:        {n_fitopt}      " \
                f"# number of FITOPTs\n")
        f.write(f"NMUOPT:         {n_muopt}      " \
                f"# number of BBC options\n")
        f.write(f"NSPLITRAN:      {n_splitran}      " \
                f"# number of random sub-samples\n")

        olam_ref = input_file_dict['p9']
        w_ref    = input_file_dict['p11']
        f.write(f"OLAM_REF:  {olam_ref} \n")
        f.write(f"w_REF:     {w_ref}  \n")
        
        f.write(f"ITER2:        {iter2} # False, True -> ITER1, ITER2 \n")
        f.write(f"SYNC_EVT:     {sync_evt} # T -> use evnts from FITOP000\n")
        f.write(f"USE_WFIT:       {use_wfit}     " \
                f"# option to run wfit on BBC output\n")
        if use_wfit :
            f.write(f"OPT_WFIT:       {CONFIG['WFITMUDIF_OPT']}    " \
                    f"# wfit arguments\n")

        f.write(f"IGNORE_FITOPT:  {ignore_fitopt} \n")
        f.write(f"IGNORE_MUOPT:   {ignore_muopt} \n")
        f.write(f"{KEY_FITOPTxMUOPT}:   {FITOPTxMUOPT} \n")

        if USE_INPDIR :
            f.write("\n")
            f.write("INPDIR_LIST:\n")
            for inpdir in inpdir_list:  f.write(f"  - {inpdir}\n")

        f.write("\n")
        f.write("SURVEY_LIST:\n")
        for survey in survey_list:    f.write(f"  - {survey}\n")

        f.write("\n")
        f.write("VERSION_OUT_LIST:\n")
        for v in vout_list :   f.write(f"  - {v}\n")


        # - - - - - -
        f.write("\n")
        f.write("FITOPT_OUT_LIST:  # 'FITOPTNUM'  'SURVEY'  " \
                f"'user_label'   'user_args'\n")
        for row in FITOPT_OUT_LIST: 
            f.write(f"  - {row}\n")

        # - - - - -
        f.write("\n")
        f.write("MUOPT_OUT_LIST:  " \
                "# 'MUOPTNUM'  'user_label'  'user_args'\n")
        for num,arg,label in zip(muopt_num_list, muopt_arg_list,
                                 muopt_label_list):
            if arg == '' : arg = None  # 9.19.2021
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")

        # end append_info_file

    def create_merge_table(self,f):

        n4d_index       = self.config_prep['n4d_index']
        iver_list4      = self.config_prep['iver_list4'] 
        ifit_list4      = self.config_prep['ifit_list4']
        imu_list4       = self.config_prep['imu_list4']
        isplitran_list4 = self.config_prep['isplitran_list4']
        n_splitran      = self.config_prep['n_splitran']
        USE_SPLITRAN    = n_splitran > 1

        # create only MERGE table ... no need for SPLIT table

        header_line_merge = \
            f" STATE   VERSION  FITOPT  MUOPT " \
            f"NEVT_DATA  NEVT_BIASCOR  NEVT_CCPRIOR  SPLITRAN"

        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state

        for iver,ifit,imu,isplitran in \
            zip(iver_list4,ifit_list4,imu_list4,isplitran_list4):

            version    = self.config_prep['version_out_list'][iver]
            version   += self.suffix_splitran(n_splitran,isplitran)
            fitopt_num = f"FITOPT{ifit:03d}"
            muopt_num  = f"{MUOPT_STRING}{imu:03d}"

            # ROW here is fragile in case columns are changed
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(version)
            ROW_MERGE.append(fitopt_num)
            ROW_MERGE.append(muopt_num)
            ROW_MERGE.append(0)    # NEVT_DATA
            ROW_MERGE.append(0)    # NEVT_BIASCOR
            ROW_MERGE.append(0)    # NEVT_CCPRIOR
            if USE_SPLITRAN :
                ROW_MERGE.append(isplitran)
            else:
                ROW_MERGE.append(None)

            INFO_MERGE['row_list'].append(ROW_MERGE)  
        util.write_merge_file(f, INFO_MERGE, [] ) 

        # end create_merge_table

    def merge_config_prep(self,output_dir):

        submit_info_yaml = self.config_prep['submit_info_yaml']
        vout_list   = submit_info_yaml['VERSION_OUT_LIST']
        n_fitopt    = submit_info_yaml['NFITOPT']
        n_muopt     = submit_info_yaml['NMUOPT']
        n_splitran  = submit_info_yaml['NSPLITRAN']

        self.config_prep['version_out_list'] = vout_list
        self.config_prep['n_splitran']       = n_splitran

        # end merge_config_prep

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        
        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_VERSION   = COLNUM_BBC_MERGE_VERSION
        COLNUM_FITOPT    = COLNUM_BBC_MERGE_FITOPT  
        COLNUM_MUOPT     = COLNUM_BBC_MERGE_MUOPT 
        COLNUM_NDATA     = COLNUM_BBC_MERGE_NEVT_DATA
        COLNUM_NBIASCOR  = COLNUM_BBC_MERGE_NEVT_BIASCOR
        COLNUM_NCCPRIOR  = COLNUM_BBC_MERGE_NEVT_CCPRIOR
        NROW_DUMP   = 0

        # keynames_for_job_stats returns 3 keynames : 
        #   {base}, {base}_sum, {base}_list
        key_ndata, key_ndata_sum, key_ndata_list = \
                self.keynames_for_job_stats('NEVT_DATA')
        key_nbiascor, key_nbiascor_sum, key_nbiascor_list = \
                 self.keynames_for_job_stats('NEVT_BIASCOR')
        key_nccprior, key_nccprior_sum, key_nccprior_list = \
                 self.keynames_for_job_stats('NEVT_CCPRIOR')
        key_list = [ key_ndata, key_nbiascor, key_nccprior  ] 

        row_list_merge   = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []

        nrow_check = 0
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input
            nrow_check += 1
            irow        = nrow_check - 1 # row index

            # strip off row info
            STATE       = row[COLNUM_STATE]

            prefix_orig, prefix_final = self.bbc_prefix("bbc", row)            
            search_wildcard = f"{prefix_orig}*"

            if irow < NROW_DUMP  :
                print(f" xxx ------------------------ ") 
                print(f" xxx DUMP merge_update_state for irow = {irow} ")
                print(f" xxx   prefix_orig  = {prefix_orig}  ->")
                print(f" xxx   prefix_final = {prefix_final}   ")

            # check if DONE or FAIL ; i.e., if Finished
            Finished = (STATE == SUBMIT_STATE_DONE) or \
                       (STATE == SUBMIT_STATE_FAIL)

            if not Finished :
                NEW_STATE = STATE

                # get list of LOG, DONE, and YAML files 
                log_list, done_list, yaml_list = \
                    util.get_file_lists_wildcard(script_dir,search_wildcard)

                # careful to sum only the files that are NOT None
                NLOG   = sum(x is not None for x in log_list)  
                NDONE  = sum(x is not None for x in done_list)  
                NYAML  = sum(x is not None for x in yaml_list)  

                if irow < NROW_DUMP:
                    print(f" xxx N(LOG,DONE,YAML) = {NLOG} {NDONE} {NYAML} ")

                if NLOG > 0:
                    NEW_STATE = SUBMIT_STATE_RUN
                if NDONE == n_job_split :
                    NEW_STATE = SUBMIT_STATE_DONE

                    bbc_stats = self.get_job_stats(script_dir,
                                                   log_list, yaml_list, key_list)
                    
                    # check for failures in snlc_fit jobs.
                    nfail = bbc_stats['nfail']
                    if nfail > 0 :  NEW_STATE = SUBMIT_STATE_FAIL
                 
                    row[COLNUM_STATE]     = NEW_STATE
                    row[COLNUM_NDATA]     = bbc_stats[key_ndata_sum]
                    row[COLNUM_NBIASCOR]  = bbc_stats[key_nbiascor_sum]
                    row[COLNUM_NCCPRIOR]  = bbc_stats[key_nccprior_sum]
                    
                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             # assume nevt changes

        # - - - - - -  -
        # The first return arg (row_split) is null since there is 
        # no need for a SPLIT table
        row_list_dict = {
            'row_split_list' : [],
            'row_merge_list' : row_list_merge_new,
            'row_extra_list' : []
        }
        return row_list_dict, n_state_change

        # end merge_update_state
        
    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):

        # For irow of MERGE.LOG,
        # copy ouput FITRES and M0DIFF files to their final ../VERSION
        # location, and remove VERSION_ from the file names since
        # they are under VERSION/subdir.
        # Example:
        #  TEST_FITOPT000_MUOPT000.FITRES is moved to 
        #  ../TEST/FITOPT000_MUOPT000.FITRES
        # ff means fitres file.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wfit         = submit_info_yaml['USE_WFIT']

        row     = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        version = row[COLNUM_BBC_MERGE_VERSION]
        prefix_orig, prefix_final = self.bbc_prefix("bbc", row)

        cddir         = f"cd {script_dir}"
        cdv           = f"cd {output_dir}/{version}"

        logging.info(f"\t Move {prefix_orig} files to {version}/ ")
        for suffix_move in SUFFIX_MOVE_LIST :
            orig_file = f"{prefix_orig}.{suffix_move}"
            move_file = f"{prefix_final}.{suffix_move}"
            cmd_move  = f"{cddir}; mv {orig_file} ../{version}/{move_file}"
            cmd_gzip  = f"gzip ../{version}/{move_file}"
            cmd_all   = f"{cmd_move} ; {cmd_gzip}"
            #print(f" xxx cmd_all = {cmd_all}")
            os.system(cmd_all)

        # check to move wfit YAML file (don't bother gzipping)
        if use_wfit :
            prefix_orig, prefix_final = self.bbc_prefix("wfit", row)
            suffix_move = "YAML"
            orig_file = f"{prefix_orig}.{suffix_move}"

            EXIST_ORIG_FILE = os.path.isfile(f"{script_dir}/{orig_file}")
            msgerr = [ f"{prefix_orig} problem", "No YAML output" ] 
            self.log_assert( EXIST_ORIG_FILE, msgerr )

            move_file = f"{prefix_final}.{suffix_move}"
            cmd_move  = f"{cddir}; mv {orig_file} ../{version}/{move_file}"
            cmd_all   = f"{cmd_move}"
            os.system(cmd_all)

        if irow == 9999 :
            sys.exit("\n xxx DEBUG DIE xxx \n")

    #end merge_job_wrapup
    
        
    def merge_cleanup_final(self):

        # Every SALT2mu job succeeded, so here we simply compress output.
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        vout_list        = submit_info_yaml['VERSION_OUT_LIST']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_splitran       = submit_info_yaml['NSPLITRAN']
        use_wfit         = submit_info_yaml['USE_WFIT']
        script_subdir    = SUBDIR_SCRIPTS_BBC

        logging.info(f"  BBC cleanup: create {FITPAR_SUMMARY_FILE}") 
        self.make_fitpar_summary()
        
        if use_wfit :
            logging.info(f"  BBC cleanup: create {WFIT_SUMMARY_FILE}")
            self.make_wfit_summary()
            
        if n_splitran > 1 :
            logging.info(f"  BBC cleanup: create {SPLITRAN_SUMMARY_FILE}")
            self.make_splitran_summary()

        for vout in vout_list :
            self.make_reject_summary(vout)

        logging.info(f"  BBC cleanup: compress {JOB_SUFFIX_TAR_LIST}")
        for suffix in JOB_SUFFIX_TAR_LIST :
            wildcard = f"{jobfile_wildcard}*.{suffix}"
            util.compress_files(+1, script_dir, wildcard, suffix, "" )

        logging.info("")

        self.merge_cleanup_script_dir()

        # end merge_cleanup_final

    def make_reject_summary(self,vout):

        # get list of all FITRES files in /vout, then find number
        # of matches for each SN. Finally, write file with
        #   VARNAMES: CID  NJOB_REJECT
        # where NJOB_REJECT is the number of FITRES files where
        # CID was rejected.
        # Goal is to use this file in 2nd round of BBC and include
        # only events that pass in all FITOPT and MUOPTs.
        #
        # Jan 12 2021: write ACCEPT file as well. Return if n_splitran>1

        n_splitran    = self.config_prep['n_splitran']
        if n_splitran > 1 : return

        output_dir    = self.config_prep['output_dir']
        VOUT          = f"{output_dir}/{vout}"
        fitres_list   = self.get_fflist_reject_summary(VOUT)

        reject_file   = BBC_REJECT_SUMMARY_FILE
        REJECT_FILE   = f"{VOUT}/{reject_file}"

        accept_file   = BBC_ACCEPT_SUMMARY_FILE
        ACCEPT_FILE   = f"{VOUT}/{accept_file}"

        logging.info(f"  BBC cleanup: create {vout}/{reject_file}")
        logging.info(f"  BBC cleanup: create {vout}/{accept_file}")

        n_ff     = len(fitres_list) # number of FITRES files
        
        first_fitres_file = VOUT + "/" + fitres_list[0]

        # check for duplicates from 
        # same data light curve measured by multiple surveys, or
        # multiple sims (e.g., LOWZ + HIGHZ) with random overlap CIDs

        df_first  = pd.read_csv(first_fitres_file, 
                                comment="#", delim_whitespace=True)
        first_cids  = df_first[TABLE_VARNAME_CID]
        first_cids_unique, counts = np.unique(first_cids,return_counts=True)
        n_dupl   = len(counts[counts>1])
        has_dupl = n_dupl > 0

        dump_dupl = False 
        if dump_dupl:
            for cid,cnt in zip(first_cids_unique,counts) :
                if cnt > 1:
                    print(f" xxx duplicate cid={cid} has cnt={cnt}")
                
        # - - - - - - - - 
        if has_dupl :
            logging.info(f"\t {n_dupl} duplicates found in first fitres file.")
            cid_dict    = self.get_cid_list_duplicates(fitres_list, VOUT)
            unique_dict = cid_dict['unique_dict']
        else:
            cid_dict = self.get_cid_list(fitres_list, VOUT)

        cid_list        = cid_dict['cid_list']
        cid_unique      = cid_dict['cid_unique']
        n_count         = cid_dict['n_count']
        n_reject        = cid_dict['n_reject']
                    
        cid_all_pass    = cid_unique[n_count == n_ff]
        cid_some_fail   = cid_unique[n_count <  n_ff]
        n_all           = len(cid_unique)
        n_some_fail     = len(cid_some_fail)
        n_all_pass      = len(cid_all_pass)
        f_some_fail     = float(n_some_fail)/float(n_all)
        str_some_fail   = f"{f_some_fail:.4f}"

        KEYVAR = KEYNAME_VARNAMES

        # - - - - - - - -
        with open(REJECT_FILE,"wt") as f:
            f.write(f"# BBC-FF = BBC FITRES file.\n")
            f.write(f"# Total number of BBC-FF: " \
                    f"{n_ff} (FITOPT x MUOPT). \n")
            f.write(f"# {n_some_fail} of {n_all} CIDs ({str_some_fail}) "\
                    f"fail cuts in 1 or more BBC-FF\n")
            f.write(f"#  and also pass cuts in 1 or more BBC-FF.\n#\n")
            f.write(f"# These CIDs are rejected in {PROGRAM_NAME_BBC} with\n")
            f.write(f"#    reject_list_file={reject_file} \n")
            f.write(f"\n")
            if has_dupl :
                f.write(f"# Beware of Duplicate CIDs "
                        f"(each CID + IDSURVEY is unique) \n")
                f.write(f"{KEYVAR}: CID IDSURVEY NJOB_REJECT \n")
                for ucid,nrej in zip(cid_unique,n_reject) :
                    cid    = unique_dict[ucid][TABLE_VARNAME_CID]
                    idsurv = unique_dict[ucid][TABLE_VARNAME_IDSURVEY]
                    if nrej>0: 
                        f.write(f"SN:  {cid:<12} {idsurv}  {nrej:3d} \n")
            else:
                f.write(f"{KEYVAR}: CID NJOB_REJECT \n")
                for cid,nrej in zip(cid_unique,n_reject) :
                    if nrej>0: 
                        f.write(f"SN:  {cid:<12}   {nrej:3d} \n")
            f.write(f"\n")

        # - - - - 

        with open(ACCEPT_FILE,"wt") as f:
            f.write(f"# BBC-FF = BBC FITRES file.\n")
            f.write(f"# Total number of BBC-FF: " \
                    f"{n_ff} (FITOPT x MUOPT). \n")
            f.write(f"# {n_all_pass} CIDs " \
                    f"pass cuts in all BBC-FF\n")
            f.write(f"# These CIDs are selected in {PROGRAM_NAME_BBC} with\n")
            f.write(f"#    accept_list_file={accept_file} \n")
            f.write(f"\n")
            if has_dupl :
                f.write(f"# Beware of Duplicate CIDs " \
                        f"(each CID + IDSURVEY is unique) \n")
                f.write(f"{KEYVAR}: {TABLE_VARNAME_CID} " \
                        f"{TABLE_VARNAME_IDSURVEY} {TABLE_VARNAME_IZBIN}\n")
                for ucid,nrej in zip(cid_unique,n_reject) :
                    cid    = unique_dict[ucid][TABLE_VARNAME_CID]
                    idsurv = unique_dict[ucid][TABLE_VARNAME_IDSURVEY]
                    izbin  = unique_dict[ucid][TABLE_VARNAME_IZBIN]
                    if nrej==0: 
                        f.write(f"SN:  {cid:<12} {idsurv}  {izbin}\n")
            else:
                f.write(f"{KEYVAR}: {TABLE_VARNAME_CID} "\
                        f" {TABLE_VARNAME_IZBIN}\n")
                izbin_unique      = cid_dict['izbin_unique']

                for cid, izbin, nrej in \
                    zip(cid_unique,izbin_unique,n_reject) :
                    if nrej==0: 
                        f.write(f"SN:  {cid:<12}  {izbin}\n")

            f.write(f"\n")

        # - - - - -
        return
        # end make_reject_summary

    def get_cid_list(self,fitres_list,VOUT):
        # get cid_list of all CIDs in all files. If same events appear in 
        # each file, each CID appears n_ff times. If a CID appears less 
        # than n_ff times, it goes into reject list.
        n_ff       = len(fitres_list)
        cid_list   = []
        izbin_list = []
        found_first_file = False

        for ff in fitres_list:
            FF       = f"{VOUT}/{ff}"
            df       = pd.read_csv(FF, comment="#", delim_whitespace=True)
            cid_list = np.concatenate((cid_list, df.CID.astype(str)))

            if not found_first_file: 
                df0 = df.copy()
                df0[TABLE_VARNAME_CID] = df0[TABLE_VARNAME_CID].astype(str)
                found_first_file = True

        # - - - - - - - - - - - - -
        # get list of unique CIDs, and how many times each CID appears
        cid_unique, n_count = np.unique(cid_list, return_counts=True)

        # number of times each CID does not appear in a fitres file
        n_reject        = n_ff - n_count

        cid_dict = {}
        cid_dict['cid_list']   = cid_list
        cid_dict['cid_unique'] = cid_unique
        cid_dict['n_count']    = n_count
        cid_dict['n_reject']   = n_reject

        # Mar 28 2022: fetch list of izbin 
        if TABLE_VARNAME_IZBIN in df0:
            izbin_unique = []
            for cid in cid_unique:
                izbin_tmplist = df0.loc[df0[TABLE_VARNAME_CID]==cid][TABLE_VARNAME_IZBIN].values
                izbin = -9
                if len(izbin_tmplist) > 0 : izbin = izbin_tmplist[0]
                izbin_unique.append(izbin)
        else:
            ncid         = len(cid_unique)
            izbin_unique = [ -9 ] * ncid

        cid_dict['izbin_unique']  = izbin_unique

        return cid_dict
        # end of get_cid_list

    def get_cid_list_duplicates(self,fitres_list,VOUT):
        # get cid_list of all CIDs in all files. If same events appear in 
        # each file, each CID appears n_ff times. If a CID appears less 
        # than n_ff times, it goes into reject list.

        devel_flag = self.config_yaml['args'].devel_flag

        n_ff        = len(fitres_list)
        unique_dict = {}
        ucid_list   = []
        found_first_file = False

        for ff in fitres_list:
            FF       = f"{VOUT}/{ff}"
            df       = pd.read_csv(FF, comment="#", delim_whitespace=True)
            ucid_list = np.concatenate((ucid_list, 
                                        df.CID.astype(str)+"__"+df.IDSURVEY.astype(str)))

            if not found_first_file:
                df0 = df.copy()
                df0[TABLE_VARNAME_CID] = df0[TABLE_VARNAME_CID].astype(str)
                found_first_file = True

        # - - - - - - - - - - - - -
        # get list of unique CIDs, and how many times each CID appears

        cid_unique, n_count = np.unique(ucid_list, return_counts=True)
        HAS_IZBIN = TABLE_VARNAME_IZBIN in df0

        for ucid in cid_unique:
            unique_dict[ucid] = {}
            cid      = str(ucid.split("__")[0])
            idsurvey = int(ucid.split("__")[1])

            if HAS_IZBIN:
                izbin_list = df0.loc[(df0[TABLE_VARNAME_CID]==cid) & \
                                     (df0[TABLE_VARNAME_IDSURVEY]==idsurvey)][TABLE_VARNAME_IZBIN].values
            else:
                izbin_list = []

            if len(izbin_list) > 0 :
                izbin = int(izbin_list[0])
            else:
                izbin = -9

            unique_dict[ucid][TABLE_VARNAME_CID]      = cid
            unique_dict[ucid][TABLE_VARNAME_IDSURVEY] = idsurvey
            unique_dict[ucid][TABLE_VARNAME_IZBIN]    = izbin

        # number of times each CID does not appear in a fitres file
        n_reject        = n_ff - n_count

        cid_dict = {}
        cid_dict['cid_list']    = ucid_list
        cid_dict['cid_unique']  = cid_unique
        cid_dict['unique_dict'] = unique_dict
        cid_dict['n_count']     = n_count
        cid_dict['n_reject']    = n_reject

        return cid_dict 

        # end of get_cid_list_duplicates

    def get_fflist_reject_summary(self,VOUT):

        # Oct 29 2020
        # get list of fitres_files needed to create reject summary.
        # Default is all files using wildcard. However, if NOREJECT
        # string is part of user label, then remove associated fitres
        # file(s) from list. This allows adding  FITOPT /NOREJECT_TEST/
        # to FITOPT /SYST_BLA/, and the NOREJECT tests won't affect the 
        # reject.list used in Pippin's 2nd BBC iteration.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        FITOPT_OUT_LIST  = submit_info_yaml['FITOPT_OUT_LIST']
        MUOPT_OUT_LIST   = submit_info_yaml['MUOPT_OUT_LIST']

        wildcard          = "FITOPT*MUOPT*.FITRES.gz"
        fitres_list_all   = sorted(glob.glob1(VOUT,wildcard))
        fitres_list       = []
        NOREJECT_list     = []

        for row in FITOPT_OUT_LIST:
            fitnum = row[0]
            label  = row[2]
            if FITOPT_STRING_NOREJECT in label : NOREJECT_list.append(fitnum)

        for row in MUOPT_OUT_LIST:
            munum = row[0]
            label = row[1]
            if FITOPT_STRING_NOREJECT in label : NOREJECT_list.append(munum)

        # loop thru fitres_list and remove anything on NOREJECT_list
        for ff in fitres_list_all :
            EXCLUDE = False
            for string in NOREJECT_list :
                if string in ff:  EXCLUDE = True
            if not EXCLUDE : fitres_list.append(ff)

        return fitres_list

        # end get_fflist_reject_summary

    def make_fitpar_summary(self):

        # Mar 28 2021: 
        # write summary info for each version/fitopt/muopt.
        # Output is YAML, but designed mainly for human readability.
        # When this function was written, there were no codes expected
        # to read this; only for human eyes.
        #
        # Nov 17 2021: 
        #     protect against FITOPT_LIST=None and MUOPT_LIST=None
        
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        iter2            = submit_info_yaml['ITER2']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_splitran       = submit_info_yaml['NSPLITRAN']
        FITOPT_LIST      = submit_info_yaml['FITOPT_OUT_LIST']
        MUOPT_LIST       = submit_info_yaml['MUOPT_OUT_LIST']
        
        if n_splitran > 1 : return

        # read the whole MERGE.LOG file to figure out where things are
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        MERGE_INFO_CONTENTS,comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        # - - - 
        SUMMARYF_FILE     = f"{output_dir}/{FITPAR_SUMMARY_FILE}"
        f = open(SUMMARYF_FILE,"wt") 
        version_last = "BLEH"

        for row in MERGE_INFO_CONTENTS[TABLE_MERGE]:
            version    = row[COLNUM_BBC_MERGE_VERSION] # sim data version
            fitopt_num = row[COLNUM_BBC_MERGE_FITOPT]  # e.g., FITOPT002
            muopt_num  = row[COLNUM_BBC_MERGE_MUOPT]   # e.g., MUOPT003
            
            # get indices for summary file
            ifit = int(f"{fitopt_num[6:]}")
            imu  = int(f"{muopt_num[5:]}")
            
            # figure out name of BBC-YAML file and read it 
            prefix_orig, prefix_final = self.bbc_prefix("bbc", row)
            YAML_FILE  = f"{script_dir}/{version}_{prefix_final}.YAML"
            #print(f"  YAML_FILE = {YAML_FILE}")
            bbc_yaml   = util.extract_yaml(YAML_FILE, None, None )
            BBCFIT_RESULTS = bbc_yaml['BBCFIT_RESULTS']

            NEVT_DATA            = bbc_yaml['NEVT_DATA']
            NEVT_BIASCOR         = bbc_yaml['NEVT_BIASCOR']
            NEVT_CCPRIOR         = bbc_yaml['NEVT_CCPRIOR']
            NEVT_REJECT_BIASCOR  = bbc_yaml['NEVT_REJECT_BIASCOR']
            frac_reject = float(NEVT_REJECT_BIASCOR)/float(NEVT_DATA)

            is_new_version  = version != version_last
            version_last = version

            if is_new_version: 
                f.write(f"\n# ==================================== \n")
                f.write(f"{version}: \n")

            f.write(f" \n")
            f.write(f"- {fitopt_num}_{muopt_num}: \n")

            # check list sizes to accomodate noINPDIR option
            n_list = 0
            if FITOPT_LIST : n_list = len(FITOPT_LIST)
            if n_list > 0 :
                f.write(f"    FITOPT: {FITOPT_LIST[ifit][3]} \n")

            n_list = 0
            if MUOPT_LIST : n_list = len(MUOPT_LIST)
            if n_list > 0 :
                f.write(f"    MUOPT:  {MUOPT_LIST[imu][2]} \n")

            f.write(f"    NEVT:   {NEVT_DATA}, {NEVT_BIASCOR}, {NEVT_CCPRIOR}"
                    f"        # DATA, BIASCOR, CCPRIOR\n")

            # - - - - - - - - -
            # check for NEVT by sample (9.19.2021)
            if 'SAMPLE_LIST' in bbc_yaml :
                # split string by commas and remove white space
                tmp         = bbc_yaml['SAMPLE_LIST']
                SAMPLE_LIST = [x.strip() for x in tmp.split(',')]

                tmp = bbc_yaml['NEVT_DATA_bySAMPLE']
                NEVT_DATA_bySAMPLE    = [x.strip() for x in tmp.split(',')]

                tmp = bbc_yaml['NEVT_BIASCOR_bySAMPLE']
                NEVT_BIASCOR_bySAMPLE = [x.strip() for x in tmp.split(',')]

                tmp = bbc_yaml['NEVT_CCPRIOR_bySAMPLE']
                NEVT_CCPRIOR_bySAMPLE = [x.strip() for x in tmp.split(',')]

                #print(f" xxx ---------------------------- ")
                #print(f" xxx {fitopt_num}_{muopt_num}")
                #print(f" xxx SAMPLE_LIST = {SAMPLE_LIST} ")
                #print(f" xxx NEVT_DATA_bySAMPLE = {NEVT_DATA_bySAMPLE}")
                
                f.write(f"    NEVT_bySAMPLE:"
                        f"                # DATA, BIASCOR, CCPRIOR\n")
                for sample,ndata,nbias,ncc in zip(SAMPLE_LIST,
                                                  NEVT_DATA_bySAMPLE,
                                                  NEVT_BIASCOR_bySAMPLE,
                                                  NEVT_CCPRIOR_bySAMPLE) :
                    key = f"{sample}:"
                    f.write(f"      {key:<20} {ndata:>5s}, {nbias:>7s}, {ncc:>4s}\n")

            # - - - - -
            nrej = NEVT_REJECT_BIASCOR
            if iter2:
                comment = f"{nrej} evts have no biasCor at ITER2; check ITER1"
            else:
                comment = f"{nrej} evts have no biasCor "
            f.write(f"    REJECT_FRAC_BIASCOR: {frac_reject:.4f}    # {comment}\n")

# - - - - - - - fix YAML output before doing frac-by-sample - - - - - 
#            f.write(f"    REJECT_FRAC_BIASCOR_bySAMPLE:\n"
# XYZ            for sample,ndata,nbias,ncc in zip(SAMPLE_LIST,
#                                                  NEVT_DATA_bySAMPLE,
#                                                  NEVT_BIASCOR_bySAMPLE,
#                                                  NEVT_CCPRIOR_bySAMPLE) :
#                    key = f"{sample}:"
#                    f.write(f"      {key:<20} {ndata:>5s}, {nbias:>7s}, {ncc:>4s}\n")
# - - - - - - - - - - - 

            # - - - - 
            for result in BBCFIT_RESULTS:
                #print(f" xxx result = {result}  key = {result.keys()} ")
                for key in result:
                    val = str(result[key].split()[0])
                    err = str(result[key].split()[1])
                    KEY = f"{key}:"
                    f.write(f"    {KEY:<12}  {val:>8} +_ {err:<8} \n")

            #f.write(f" \n")
            #f.write(f" \n")
        f.close()

        #sys.exit("\n xxx DEBUG STOP in make_fitpar_summary\n")

        # end make_fitpar_summary

    def make_wfit_summary(self):
        
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wfit         = submit_info_yaml['USE_WFIT']
        n_splitran       = submit_info_yaml['NSPLITRAN']
        opt_wfit         = submit_info_yaml['OPT_WFIT']
        use_wfit_w0wa    = '-wa'    in opt_wfit
        use_wfit_blind   = '-blind' in opt_wfit

        # - - - 
        SUMMARYF_FILE     = f"{output_dir}/{WFIT_SUMMARY_FILE}"
        f = open(SUMMARYF_FILE,"w") 

        varname_w   = "w"
        if use_wfit_w0wa : varname_w = "w0"
        varname_wa  = "wa"
        varname_omm = "omm"
        
        # prepare w-varnames for varnames
        varlist_w = f"{varname_w} {varname_w}sig"
        if use_wfit_w0wa :
            varlist_w += f" {varname_wa} {varname_wa}sig"  #w0waCDM


        varnames = f"VARNAMES: ROW VERSION FITOPT MUOPT  " \
                   f"{varlist_w}  {varname_omm} {varname_omm}_sig  "\
                   f"chi2 sigint   \n"

        # read the whole MERGE.LOG file to figure out where things are
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        MERGE_INFO_CONTENTS,comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        nrow = 0 
        for row in MERGE_INFO_CONTENTS[TABLE_MERGE]:
            nrow += 1
            version    = row[COLNUM_BBC_MERGE_VERSION] # sim data version
            fitopt_num = row[COLNUM_BBC_MERGE_FITOPT]  # e.g., FITOPT002
            muopt_num  = row[COLNUM_BBC_MERGE_MUOPT]   # e.g., MUOPT003
            isplitran  = row[COLNUM_BBC_MERGE_SPLITRAN]
            
            # get indices for summary file
            ifit = f"{fitopt_num[6:]}"
            imu  = f"{muopt_num[5:]}"
            
            # figure out name of wfit-YAML file and read it
            prefix_orig,prefix_final = self.bbc_prefix("wfit", row)
            YAML_FILE  = f"{output_dir}/{version}/{prefix_final}.YAML"
            wfit_yaml  = util.extract_yaml(YAML_FILE, None, None )

            # extract wfit values into local variables
            wfit_values_dict = util.get_wfit_values(wfit_yaml)

            w       = wfit_values_dict['w']  
            w_sig   = wfit_values_dict['w_sig']
            wa      = wfit_values_dict['wa']    
            wa_sig  = wfit_values_dict['wa_sig']
            omm     = wfit_values_dict['omm']  
            omm_sig = wfit_values_dict['omm_sig']
            chi2    = wfit_values_dict['chi2'] 
            sigint  = wfit_values_dict['sigint']
            w_ran   = int(wfit_values_dict['w_ran']) 
            wa_ran  = int(wfit_values_dict['wa_ran'])
            omm_ran = int(wfit_values_dict['omm_ran'])

            if use_wfit_w0wa :
                w0 = w ; w0_sig = w_sig
                w_values = f"{w0:7.4f} {w0_sig:6.4f} {wa:7.4f} {wa_sig:6.4f}"
            else:
                w_values = f"{w:7.4f} {w_sig:6.4f}"

            string_values = \
                f"{nrow:3d}  {version} {ifit} {imu} " \
                f"{w_values}  {omm:6.4f} {omm_sig:6.4f} " \
                f"{chi2:.1f} {sigint:.3f} "

            if nrow == 1 and use_wfit_blind: 
                f.write(f"# cosmology params blinded.\n")
                f.write(f"#   {varname_w:<3} includes sin({w_ran}) \n")
                if use_wfit_w0wa :
                    f.write(f"#   {varname_wa:<3} includes sin({wa_ran}) \n")
                f.write(f"#   {varname_omm:<3} includes sin({omm_ran}) \n")
                f.write(f"\n")

            if nrow==1 : 
                f.write(f"{varnames}\n")

            f.write(f"{KEY_ROW} {string_values}\n")

        f.close()

        # end make_wfit_summary

    def make_splitran_summary(self):

        # collect all BBC fit params, and optional w(wfit);
        # write them out into a FITRES-formatted text file.
        # Include column indices for FITOPT and MUOPT.

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        n_splitran       = submit_info_yaml['NSPLITRAN']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wfit         = submit_info_yaml['USE_WFIT']
        vout_list        = submit_info_yaml['VERSION_OUT_LIST']

        SUMMARYF_FILE     = f"{output_dir}/{SPLITRAN_SUMMARY_FILE}"
        f = open(SUMMARYF_FILE,"w") 

        self.write_splitran_comments(f)
        self.write_splitran_header(f)

        # read the whole MERGE.LOG file to figure out where things are
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        MERGE_INFO_CONTENTS,comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        nrow = 0 
        for row in MERGE_INFO_CONTENTS[TABLE_MERGE]:
            version    = row[COLNUM_BBC_MERGE_VERSION] # sim data version
            fitopt_num = row[COLNUM_BBC_MERGE_FITOPT]  # e.g., FITOPT002
            muopt_num  = row[COLNUM_BBC_MERGE_MUOPT]   # e.g., MUOPT003
            isplitran  = row[COLNUM_BBC_MERGE_SPLITRAN]
            
            # remove suffix from version to get base version
            suffix       = self.suffix_splitran(n_splitran,isplitran)
            len_base     = len(version) - len(suffix)
            version_base = version[0:len_base]
            sys.stdout.flush()

            # get indices for summary file
            iver = vout_list.index(version_base)
            ifit = f"{fitopt_num[6:]}"
            imu  = f"{muopt_num[5:]}"

            # process all splitran files upon reaching SPLITRAN=1
            # in MERGE.LOG file
            if isplitran > 1 : continue
            nrow += 1  # for row number in summary file

            # the ugly code is in get_splitran_values 
            varname_list, value_list2d, error_list2d = \
                    self.get_splitran_values(row)

            # for each list of values, get statistics, then print to table.
            n_var = len(varname_list)
            for ivar in range(0,n_var):
                varname    = varname_list[ivar]
                value_list = value_list2d[ivar]
                error_list = error_list2d[ivar]
                stat_dict  = util.get_stat_dict(value_list,error_list)
                AVG_VAL    = stat_dict['AVG_VAL'] # avg fit value
                AVG_ERR    = stat_dict['AVG_ERR'] # avg fit error
                ERR_AVG    = stat_dict['ERR_AVG'] # error on mean 
                RMS        = stat_dict['RMS']     # RMS on fit values
                ERR_RMS    = stat_dict['ERR_RMS'] # error on RMS

                string_values = \
                    f"{nrow:3d} {iver} {ifit} {imu} {varname:<10} "\
                    f"{AVG_VAL:8.4f} {ERR_AVG:8.4f} " \
                    f"{AVG_ERR:8.4f} {RMS:8.4f} {ERR_RMS:8.4f} "

                f.write(f"{KEY_ROW} {string_values}\n")

            f.write(f"\n")

            if nrow == 77777 : break  # debug only

        f.close()
        # end make_splitran_summary
    
    def get_splitran_values(self,row):

        # for input row from MERGE.LOG, return
        #   varnames_list (list of variables names with BBC results)
        #   values_list2d (list of values vs. splitran for each variable)

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wfit         = submit_info_yaml['USE_WFIT']
        n_splitran       = submit_info_yaml['NSPLITRAN']
        split_string     = self.suffix_splitran(n_splitran,1)

        version          = row[COLNUM_BBC_MERGE_VERSION]
        version_base     = version.rstrip(split_string) # remove "-0001"
        prefix_orig,prefix_final = self.bbc_prefix("bbc", row)

        # scoop up YAML files. Be careful that '-{isplitran} is the
        # part we need to exlude from prefix, but we don't want to 
        # remove other dashes in version name.

        prefix_search = prefix_orig.replace(split_string,'*')
        wildcard_yaml = f"{prefix_search}*.YAML"
        yaml_list     = glob.glob1(script_dir, wildcard_yaml)

        bbc_results_yaml   = []
        for yaml_file in yaml_list :
            YAML_FILE = f"{script_dir}/{yaml_file}"
            tmp_yaml  = util.extract_yaml(YAML_FILE, None, None )
            n_var     = len(tmp_yaml['BBCFIT_RESULTS'])
            bbc_results_yaml.append(tmp_yaml)

        if use_wfit :  n_var += 1 ;  

        # Make list of varnames[ivar] and value_list2d[ivar][isplitran]
        # Trick is to convert [isplitran][ivar] -> [ivar][isplitran]
        #   (I hate this code)
        
        #print(f" xxx -----------------------")
        #print(f" xxx bbc_results_yaml = {bbc_results_yaml} ")
        #print(f" xxx n_var = {n_var} ")

        varname_list = []
        value_list2d = [ 0.0 ] * n_var  # [ivar][isplitran]
        error_list2d = [ 0.0 ] * n_var

        for ivar in range(0,n_var): 
            value_list2d[ivar] = []
            error_list2d[ivar] = []
        isplitran    = 0
        
        for results in bbc_results_yaml:  # loop over splitran
            BBCFIT_RESULTS = results['BBCFIT_RESULTS']
            ivar = 0 
            for item in BBCFIT_RESULTS:  # loop over variables
                #print(f" xxx check item({ivar}) = {item}")
                for key,val in item.items() :
                    str_val = str(val).split()[0]
                    str_err = str(val).split()[1]
                    value_list2d[ivar].append(float(str_val))
                    error_list2d[ivar].append(float(str_err))
                    if isplitran == 0 : varname_list.append(key)
                ivar += 1
            isplitran += 1

        # - - - - - - - - 
        # check option to include w(wfit)
        # Note that wfit_*YAML files have already been moved to 
        # ../version-[splitran] and had version string removed from 
        # file name; therefore, use prefix_final instead of prefix_orig

        if use_wfit :
            ivar = n_var - 1
            w_list = [] ;  werr_list = []
            varname_list.append("w_wfit")
            prefix_orig,prefix_final = self.bbc_prefix("wfit", row)
            yaml_file     = f"{prefix_final}.YAML"
            v_wildcard    = f"{output_dir}/{version_base}*"
            yaml_list     = glob.glob(f"{v_wildcard}/{yaml_file}")

            #print(f" xxx wildcard = {v_wildcard}/{yaml_file}")
            #print(f" xxx yaml_list = {yaml_list} ") 

            for yaml_file in yaml_list :
                tmp_yaml  = util.extract_yaml(yaml_file, None, None )
                wfit_values_dict = util.get_wfit_values(tmp_yaml)
                w       = wfit_values_dict['w']  
                w_sig   = wfit_values_dict['w_sig']
                w_list.append(w)
                werr_list.append(w_sig)
            value_list2d[ivar] = w_list
            error_list2d[ivar] = werr_list
        
        return varname_list, value_list2d, error_list2d

        # end get_splitran_values

    def write_splitran_comments(self, f):

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wfit         = submit_info_yaml['USE_WFIT']
        vout_list        = submit_info_yaml['VERSION_OUT_LIST']
        muopt_list       = submit_info_yaml['MUOPT_OUT_LIST']
        n_splitran       = submit_info_yaml['NSPLITRAN']

        f.write(f"# ========================================= \n")
        f.write(f"# NSPLITRAN: {n_splitran} \n#\n")
        # write comments with map if IVER -> VERSION
        iver = 0
        for vout in vout_list :
            f.write(f"# IVER = {iver:2d} --> {vout} \n")
            iver += 1
        f.write(f"#\n")

        # element 0: FITOPTnnn
        # element 1: optional label
        # element 2: arg list
        for muopt in muopt_list :
            muopt_num = muopt[0]
            muopt_arg = muopt[2]
            f.write(f"# {muopt_num}: {muopt_arg} \n")

        f.write(f"#\n")
        f.write(f"# AVG_VAL = sum[VALUES]/NSNFIT   # avg fit value \n")
        f.write(f"# AVG_ERR = sum[ERRORS]/NSNFIT   # avg fit error \n")
        f.write(f"# RMS     = r.m.s of VALUES      # rms if fit values \n")
        f.write(f"# ERR_AVG = RMS/sqrt(NSNFIT)     # error on mean \n")
        f.write(f"# ERR_RMS = RMS/sqrt(2*NSNFIT)   # error on rms  \n")
        f.write(f"# ============================================== \n\n")

        #end write_splitran_comments

    def write_splitran_header(self, f):
        # write header using names under BBCFIT_RESULTS
        f.write("VARNAMES: ROW IVER FITOPT MUOPT  PARNAME  " \
                f"AVG_VAL  ERR_AVG  AVG_ERR  RMS  ERR_RMS \n")

    def merge_reset(self,output_dir):

        # unpack things in merge_cleanup_final, but in reverse order

        submit_info_yaml = self.config_prep['submit_info_yaml']
        vout_list        = submit_info_yaml['VERSION_OUT_LIST']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        script_subdir    = SUBDIR_SCRIPTS_BBC
        fnam = "merge_reset"

        logging.info(f"   {fnam}: reset STATE and NEVT in {MERGE_LOG_FILE}")
        MERGE_LOG_PATHFILE = f"{output_dir}/{MERGE_LOG_FILE}"
        colnum_zero_list = [ COLNUM_BBC_MERGE_NEVT_DATA, 
                             COLNUM_BBC_MERGE_NEVT_BIASCOR,
                             COLNUM_BBC_MERGE_NEVT_CCPRIOR ]
        util.merge_table_reset(MERGE_LOG_PATHFILE, TABLE_MERGE,  \
                               COLNUM_MERGE_STATE, colnum_zero_list)

        # xxx unpack_script_dir(scrippt_dir, jobfile_wildcard)
        util.untar_script_dir(script_dir)

        # xxxxxxxxx mark delete Oct 1 2022 xxxxxx
        #logging.info(f"  {fnam}: uncompress {script_subdir}/")
        #util.compress_subdir(-1, f"{output_dir}/{script_subdir}" )

        #logging.info(f"  {fnam}: uncompress {JOB_SUFFIX_TAR_LIST}")
        #for suffix in JOB_SUFFIX_TAR_LIST :
        #    wildcard = f"{jobfile_wildcard}*.{suffix}"
        #    util.compress_files(-1, script_dir, wildcard, suffix, "" )

        #logging.info(f"  {fnam}: uncompress CPU* files")
        #util.compress_files(-1, script_dir, "CPU*", "CPU", "" )
        # xxxxxxxxx end mark xxxxxxxx

        # - - - - - - - - - 
        logging.info(f"  {fnam}: restore {SUFFIX_MOVE_LIST} to {script_subdir}")
        for vout in vout_list : 
            vout_dir = f"{output_dir}/{vout}"
            cdv      = f"cd {vout_dir}"
            logging.info(f" \t\t restore {vout}")
            for suffix_move in SUFFIX_MOVE_LIST :
                #logging.info(f" \t\t restore {vout}/*{suffix_move}")
                wildcard  = f"FITOPT*.{suffix_move}"
                cmd_unzip = f"{cdv} ; gunzip {wildcard}.gz"
                os.system(cmd_unzip)

            # restore each file with version_ appended to name
            ff_list = sorted(glob.glob1(vout_dir,"FITOPT*"))
            #print(f"\t xxx ff_list = {ff_list} ")
            for ff in ff_list:
                ff_move   = f"{script_dir}/{vout}_{ff}"
                cmd_move  = f"mv {ff} {ff_move}"
                cmd_all   = f"{cdv} ; {cmd_move}"
                os.system(cmd_all)

            # DILLON June 24th 2021: REPEAT FOR wfit YAML      
            if len(glob.glob1(vout_dir,"wfit*")) > 0:
                wfit_list = sorted(glob.glob1(vout_dir,"wfit*"))
                for wff in wfit_list: # wfit filenames in OUTPUT directory               
                    wf = wff[len(PREFIX_wfit)+1:] # shortened temporary version 
                    wf_move   = f"{script_dir}/wfit_{vout}_{wf}" # for SCRIPTS directory 
                    cmd_move  = f"mv {wff} {wf_move}"
                    cmd_all   = f"{cdv} ; {cmd_move}"
                    os.system(cmd_all)

        # end merge_reset

    def bbc_prefix(self, program, row):

        # Input program can be
        #   SALT2mu or bbc -> nominal code
        #   wift           -> tack on extra 'wfit' to prefix
        #
        # Input row has the format of a row in MERGE.LOG
        #
        # Function returns 
        #    prefix_orig  (includes version_ )
        #    prefix_final (does not include version_

        version       = row[COLNUM_BBC_MERGE_VERSION]
        fitopt_num    = row[COLNUM_BBC_MERGE_FITOPT]
        muopt_num     = row[COLNUM_BBC_MERGE_MUOPT]

        n_splitran    = self.config_prep['n_splitran']
        isplitran     = 1  # default is no splitran
        if n_splitran > 1 : isplitran  = row[COLNUM_BBC_MERGE_SPLITRAN]
        suffix        = self.suffix_splitran(n_splitran,isplitran)

        # if suffix is already part of version, don't re-apply it.
        if len(suffix) > 3 and suffix in version :   suffix = ''

        prefix_orig   = f"{version}{suffix}_{fitopt_num}_{muopt_num}"
        prefix_final  = f"{fitopt_num}_{muopt_num}"

        # check for adding 'wfit' to prefix
        if program.lower() == 'wfit' :
            prefix_orig  = f"wfit_{prefix_orig}"
            prefix_final = f"wfit_{prefix_final}"

        return prefix_orig, prefix_final

        # end bbc_prefix

    def suffix_splitran(self,n_splitran,isplitran):
        suffix = ''
        if n_splitran > 1 : suffix = f"-{isplitran:04d}"
        return suffix

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file  
        submit_info_yaml = self.config_prep['submit_info_yaml']
        survey_list      = submit_info_yaml['SURVEY_LIST']

        info_lines = []
        info_lines.append(f"SURVEY_LIST:  {survey_list}")
        return info_lines
        # end get_misc_merge_info 

    def get_merge_COLNUM_CPU(self):
        return -9  # there is no CPU column


