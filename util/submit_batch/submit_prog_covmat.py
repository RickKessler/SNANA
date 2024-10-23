# Created Sep 30 2022 by R.Kessler
#
# Read output of BBC and run create_covariance.py on a 3D grid of
#  BBC_OUTDIR x BBC_SUBDIR x COVMATOPT
#
# Typically N BBC_SUBDIRs reflect N independent sim samples.
#

import os, sys, shutil, yaml, glob
import logging
#import coloredlogs
import datetime, time
import submit_util as util
import numpy as np

from submit_params    import *
from submit_prog_base import Program

# keys in config file
KEY_BBC_OUTDIR     = 'BBC_OUTDIR'
KEY_COVMATOPT      = 'COVMATOPT'
KEY_VERSION_SUBSET = 'VERSION_SUBSET'

# keys in native input file
KEY_SYS_SCALE_FILE = 'SYS_SCALE_FILE'

# define columns for MERGE.LOG;  column 0 is always for STATE
COLNUM_COVMAT_MERGE_COVMATOPT    = 1
COLNUM_COVMAT_MERGE_BBCDIR       = 2
COLNUM_COVMAT_MERGE_SUBDIR       = 3
COLNUM_COVMAT_MERGE_NCOVMAT      = 4 
COLNUM_COVMAT_MERGE_CPU          = 5

PREFIX_JOB_FILES = 'COVMAT'  # for LOG, DONE, YAML 

# - - - - - - - - - - - - - - - - - - -  -
class create_covmat(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_COVMAT
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

        return output_dir_name, SUBDIR_SCRIPTS_COVMAT
        # end set_output_dir_name

    def submit_prepare_driver(self):

        print('')
        # read yaml input from create_cov input file
        self.prep_covmat_read_inputs()

        # read and store BBC outdirs and subdirs
        self.prep_covmat_find_bbc_outdirs()

        # read and store covmat options
        self.prep_covmat_options()

        # copy input files to script-dir
        self.prep_covmat_copyFiles()

        # prepare internal 3D indices 
        self.prep_covmat_index_lists()

        print('')
        return
        # end submit_prepare_driver


    def prep_covmat_read_inputs(self):

        CONFIG     = self.config_yaml['CONFIG']
        input_covmat_file = os.path.expandvars(CONFIG['INPUT_COVMAT_FILE'])
        
        input_lines = []
        with open(input_covmat_file, 'rt') as f :
            for line in f: 
                if util.is_comment_line(line) : continue
                input_lines.append(line)
        input_covmat_yaml = yaml.safe_load("\n".join(input_lines))

        self.config_prep['input_covmat_file'] = input_covmat_file
        self.config_prep['input_covmat_yaml'] = input_covmat_yaml
        return
        # end prep_covmat_bbc_input_dirs

    def prep_covmat_find_bbc_outdirs(self):

        # Read list of BBC outdirs and then
        # read MERGE.LOG from bbc output to identify all of the
        # output directories that contain FITRES files.
        # These BBC outputs are inputs to create_covariance.

        CONFIG          = self.config_yaml['CONFIG']
        bbc_outdir_rows = CONFIG[KEY_BBC_OUTDIR]

        if KEY_VERSION_SUBSET in CONFIG:
            version_subset = CONFIG[KEY_VERSION_SUBSET]
        else:
            version_subset = []

        bbc_outdir_dict = \
            util.prep_jobopt_list(bbc_outdir_rows, KEY_BBC_OUTDIR, 0, None)
        
        # abort if any label is None
        util.require_jobopt_labels(bbc_outdir_dict)

        #print(f"\n xxx bbc_outdir_dict = \n{bbc_outdir_dict} \n")

        COLNUM_BBC_MERGE_VERSION   = 1
        bbc_rawdir_list  = bbc_outdir_dict['jobopt_arg_list']

        bbc_outdir_list  = []  # list of BBC outdirs
        bbc_subdir_list2 = []  # each BBC outdir can have subdir list
                               # e.g., 50 independent sim samples
        n_subdir_list = []

        for outdir in bbc_rawdir_list:
            outdir = os.path.expandvars(outdir)
            bbc_outdir_list.append(outdir)
            merge_log = f"{outdir}/{MERGE_LOG_FILE}"            
            merge_yaml, comment_lines = util.read_merge_file(merge_log) 
            
            bbc_subdir_list = []
            for row in merge_yaml['MERGE']:
                version = row[COLNUM_BBC_MERGE_VERSION]
                if version in bbc_subdir_list: continue
                keep = True
                if len(version_subset) > 0:
                    keep = False
                    for subset in version_subset:
                        if subset in version: keep = True
                if keep: 
                    bbc_subdir_list.append(version)

            bbc_subdir_list2.append(bbc_subdir_list)
            n_subdir_list.append(len(bbc_subdir_list))

        # - - - - - - -
        bbc_outdir_dict['outdir_list']  = bbc_outdir_list
        bbc_outdir_dict['subdir_list2'] = bbc_subdir_list2
        
        # print summary info
        n_bbc_outdir = len(bbc_rawdir_list)
        print(f"  Found {n_bbc_outdir} BBC outdirs")
        print(f"  Number of subdirs per outdir: {n_subdir_list}")
        #print(f"\n xxx bbc_outdir_dict = \n{bbc_outdir_dict}\n")
                                 
        self.config_prep['bbc_outdir_dict']  = bbc_outdir_dict

        return
        # end prep_covmat_find_bbc_outdirs

    def prep_covmat_options(self):

        CONFIG          = self.config_yaml['CONFIG']
        covmatopt_rows = CONFIG[KEY_COVMATOPT]

        covmatopt_dict = \
            util.prep_jobopt_list(covmatopt_rows, KEY_COVMATOPT, 0, None)

        # abort if any label is None
        util.require_jobopt_labels(covmatopt_dict)

        n_covmatopt = covmatopt_dict['n_jobopt']
        print(f"  Found {n_covmatopt} COVMAT options")

        self.config_prep['covmatopt_dict']  = covmatopt_dict
        #print(f"\n xxx covmatopt_dict = \n{covmatopt_dict}\n")

        return
        # end prep_covmat_options
              
    def prep_covmat_copyFiles(self):

        # copy input file to outdir
        input_master_file = self.config_yaml['args'].input_file
        input_covmat_file = self.config_prep['input_covmat_file']
        input_covmat_yaml = self.config_prep['input_covmat_yaml']
        script_dir        = self.config_prep['script_dir']

        shutil.copy(input_master_file, script_dir)
        shutil.copy(input_covmat_file, script_dir)

        if KEY_SYS_SCALE_FILE in input_covmat_yaml:
            sys_scale_file = input_covmat_yaml[KEY_SYS_SCALE_FILE]
            if isinstance(sys_scale_file,str):
                shutil.copy(sys_scale_file, script_dir)
        
        return
        # end prep_covmat_copyFiles

    def prep_covmat_index_lists(self):

        CONFIG             = self.config_yaml['CONFIG']
        bbc_outdir_dict    = self.config_prep['bbc_outdir_dict']
        covmatopt_dict     = self.config_prep['covmatopt_dict'] 

        # extract the 3 lists used for 3D loop
        # bbc_outdir and arg have fixed list size;
        # subdir_list2 might depend on bbc_outdir

        bbc_outdir_list    = bbc_outdir_dict['outdir_list']  
        bbc_subdir_list2   = bbc_outdir_dict['subdir_list2'] 
        arg_list           = covmatopt_dict['jobopt_arg_list']

        n_bbc_outdir = len(bbc_outdir_list)
        n_arg_covmat = len(arg_list)

        icovmat_list3 = []; idir_list3 = [];  isubdir_list3 = [];   
        n_job_tot = 0

        for icovmat in range(0,n_arg_covmat):
            for idir in range(0,n_bbc_outdir):
                bbc_subdir_list = bbc_subdir_list2[idir]
                n_bbc_subdir    = len(bbc_subdir_list)
                for sdir in range(0,n_bbc_subdir):
                    n_job_tot += 1
                    icovmat_list3.append(icovmat)
                    idir_list3.append(idir)
                    isubdir_list3.append(sdir)
                    
        # - - - - - -
        print(f"  Prep index lists for {n_job_tot} COVMAT jobs")

        self.config_prep['n_job_split'] = 1
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        self.config_prep['icovmat_list3']  = icovmat_list3
        self.config_prep['idir_list3']     = idir_list3
        self.config_prep['isubdir_list3']  = isubdir_list3

 
        return
        # end prep_covmat_index_lists

    def write_command_file(self, icpu, f):

        n_core             = self.config_prep['n_core']
        bbc_outdir_dict    = self.config_prep['bbc_outdir_dict']
        covmatopt_dict     = self.config_prep['covmatopt_dict'] 

        bbc_outdir_list    = bbc_outdir_dict['outdir_list']  
        bbc_subdir_list2   = bbc_outdir_dict['subdir_list2'] 
        arg_list           = covmatopt_dict['jobopt_arg_list']

        icovmat_list3 = self.config_prep['icovmat_list3']
        idir_list3    = self.config_prep['idir_list3'] 
        isubdir_list3 = self.config_prep['isubdir_list3']
        
        n_job_cpu   = 0
        n_job_local = 0
        
        for icovmat,idir,isubdir in \
            zip(icovmat_list3,idir_list3,isubdir_list3):

            n_job_local += 1
            index_dict = { 'icovmat':icovmat, 'idir':idir, 'isubdir':isubdir, 
                           'icpu':icpu }

            if ( (n_job_local-1) % n_core ) != icpu : continue

            n_job_cpu += 1
            job_info_covmat   = self.prep_JOB_INFO_covmat(index_dict)
            util.write_job_info(f, job_info_covmat, icpu)

            job_info_merge = \
                self.prep_JOB_INFO_merge(icpu,n_job_local,False) 
            util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu
        # end write_command_file

    def prep_JOB_INFO_covmat(self,index_dict):

        icovmat = index_dict['icovmat']
        idir    = index_dict['idir']
        isubdir = index_dict['isubdir']

        kill_on_fail = self.config_yaml['args'].kill_on_fail
        program      = self.config_prep['program']
        output_dir   = self.config_prep['output_dir']
        script_dir   = self.config_prep['script_dir']

        input_covmat_file  = self.config_prep['input_covmat_file']
        bbc_outdir_dict    = self.config_prep['bbc_outdir_dict']
        covmatopt_dict     = self.config_prep['covmatopt_dict'] 

        bbc_outdir         = bbc_outdir_dict['outdir_list'][idir]  
        bbc_subdir         = bbc_outdir_dict['subdir_list2'][idir][isubdir]
        bbc_outdir_label   = bbc_outdir_dict['jobopt_label_list'][idir]

        arg_string         = covmatopt_dict['jobopt_arg_list'][icovmat]
        arg_label          = covmatopt_dict['jobopt_label_list'][icovmat]

        job_label = f"{arg_label}_{bbc_outdir_label}_{bbc_subdir}"

        prefix        = f"{PREFIX_JOB_FILES}_{job_label}"
        log_file      = f"{prefix}.LOG" 
        done_file     = f"{prefix}.DONE"
        yaml_file     = f"{prefix}.YAML"
        all_done_file = f"{output_dir}/{DEFAULT_DONE_FILE}"

        covmat_outdir = f"{output_dir}/{job_label}"

        # define command-line arguments
        arg_list = []
        arg_list.append(f"--input_dir {bbc_outdir}")
        arg_list.append(f"--version   {bbc_subdir}")
        arg_list.append(f"--outdir    {covmat_outdir}")
        arg_list.append(f"--yaml_file {yaml_file}")
        arg_list.append(arg_string)

        JOB_INFO = {}
        JOB_INFO['program']       = program
        JOB_INFO['input_file']    = input_covmat_file
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = f"{log_file}"
        JOB_INFO['done_file']     = f"{done_file}"
        JOB_INFO['all_done_file'] = f"{all_done_file}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['arg_list']      = arg_list
  
        return JOB_INFO


    def append_info_file(self,f):
        CONFIG             = self.config_yaml['CONFIG']

        f.write(f"\n# covmat info\n")

        f.write(f"JOBFILE_WILDCARD:  '{PREFIX_JOB_FILES}*' \n\n")

        f.write(f"{KEY_BBC_OUTDIR}: \n")
        for row in CONFIG[KEY_BBC_OUTDIR] :
            f.write(f"- {row} \n")

        f.write("\n")

        f.write(f"{KEY_COVMATOPT}: \n")
        for row in CONFIG[KEY_COVMATOPT] :
            f.write(f"- {row} \n")

        return
        # end append_info_file

    def create_merge_table(self,f):

        # create merge table with NDOF=0; later this table
        # gets updated as jobs finish.
        icovmat_list3   = self.config_prep['icovmat_list3']
        idir_list3      = self.config_prep['idir_list3']
        isubdir_list3   = self.config_prep['isubdir_list3'] 

        bbc_outdir_dict    = self.config_prep['bbc_outdir_dict']
        covmatopt_dict     = self.config_prep['covmatopt_dict'] 
        
        bbc_label_list   = bbc_outdir_dict['jobopt_label_list']
        bbc_subdir_list2 = bbc_outdir_dict['subdir_list2']
        covmatopt_label_list = covmatopt_dict['jobopt_label_list']

        # create only MERGE table ... no need for SPLIT table
        header_line_merge = \
                f" STATE  COVMATOPT  BBCDIR  SUBDIR  NCOVMAT  CPU "

        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line_merge,
            'row_list'    : []   
        }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state

        for icovmat,idir,isubdir in \
            zip(icovmat_list3,idir_list3,isubdir_list3):

            covmatopt_label  = covmatopt_label_list[icovmat]
            bbc_label        = bbc_label_list[idir]
            bbc_subdir       = bbc_subdir_list2[idir][isubdir]

            # ROW here is fragile in case columns are changed
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(covmatopt_label)
            ROW_MERGE.append(bbc_label)
            ROW_MERGE.append(bbc_subdir)
            ROW_MERGE.append(0)       # N_COVMAT
            ROW_MERGE.append(0.0)     # CPU
            INFO_MERGE['row_list'].append(ROW_MERGE)  

        # - - - - -
        util.write_merge_file(f, INFO_MERGE, [] ) 

        return
        # end create_merge_table

    def merge_config_prep(self,output_dir):
        submit_info_yaml = self.config_prep['submit_info_yaml']
 

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_COVMATOPT = COLNUM_COVMAT_MERGE_COVMATOPT
        COLNUM_BBCDIR    = COLNUM_COVMAT_MERGE_BBCDIR
        COLNUM_SUBDIR    = COLNUM_COVMAT_MERGE_SUBDIR
        COLNUM_NCOVMAT   = COLNUM_COVMAT_MERGE_NCOVMAT
        COLNUM_CPU       = COLNUM_COVMAT_MERGE_CPU
        NROW_DUMP   = 0

        key_ncov, key_ncov_sum, key_ncov_list = \
                self.keynames_for_job_stats('N_COVMAT')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')

        key_list = [ key_ncov, key_cpu ] 

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
            prefix      = self.covmat_prefix(row) 
            search_wildcard = f"{PREFIX_JOB_FILES}_{prefix}*"

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

                if NLOG > 0:
                    NEW_STATE = SUBMIT_STATE_RUN
                if NDONE == n_job_split :
                    NEW_STATE = SUBMIT_STATE_DONE
     
                    covmat_stats = self.get_job_stats(script_dir,
                                                      log_list, 
                                                      yaml_list, 
                                                      key_list)
                    
                    # check for failures in snlc_fit jobs.
                    nfail = covmat_stats['nfail']
                    if nfail > 0 :  NEW_STATE = SUBMIT_STATE_FAIL
                 
                    row[COLNUM_STATE]     = NEW_STATE
                    row[COLNUM_NCOVMAT]   = covmat_stats[key_ncov_sum]
                    row[COLNUM_CPU]       = covmat_stats[key_cpu_sum]
                    
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
        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        row     = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        # end merge_job_wrapup

    def merge_cleanup_final(self):

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        script_subdir    = SUBDIR_SCRIPTS_COVMAT

        logging.info(f"  covmat cleanup: compress {JOB_SUFFIX_TAR_LIST}")
        for suffix in JOB_SUFFIX_TAR_LIST :
            wildcard = (f"{PREFIX_JOB_FILES}*.{suffix}") 
            util.compress_files(+1, script_dir, wildcard, suffix, "" )

        logging.info("")

        self.merge_cleanup_script_dir()

        # end merge_cleanup_final

    def covmat_prefix(self,row):
        # parse input row passed from MERGE.LOG and construct
        # prefix for output files
        covmatopt  = row[COLNUM_COVMAT_MERGE_COVMATOPT]  
        bbcdir     = row[COLNUM_COVMAT_MERGE_BBCDIR] 
        subdir     = row[COLNUM_COVMAT_MERGE_SUBDIR]
        prefix     = f"{covmatopt}_{bbcdir}_{subdir}"
        return prefix

        # end covmat_prefix

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file  
        submit_info_yaml = self.config_prep['submit_info_yaml']
 
        return []
        # end get_misc_merge_info

    def get_merge_COLNUM_CPU(self):
        return COLNUM_COVMAT_MERGE_CPU

    # .xyz END
