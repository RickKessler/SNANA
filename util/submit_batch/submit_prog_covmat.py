# Created Sep 30 2022 by R.Kessler
#
# Read output of BBC and run create_covariance.

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time
import submit_util as util
import numpy as np

from submit_params    import *
from submit_prog_base import Program

# keys in config file
KEY_BBC_OUTDIR = 'BBC_OUTDIR'
KEY_COVMATOPT  = 'COVMATOPT'

# keys in native input file
KEY_SYS_SCALE_FILE = 'SYS_SCALE_FILE'

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
                if version not in bbc_subdir_list:
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

        idir_list3 = [];  isubdir_list3 = [];   icovmat_list3 = []
        n_job_tot = 0

        for idir in range(0,n_bbc_outdir):
            bbc_subdir_list = bbc_subdir_list2[idir]
            n_bbc_subdir    = len(bbc_subdir_list)
            for sdir in range(0,n_bbc_subdir):
                for icovmat in range(0,n_arg_covmat):
                    n_job_tot += 1
                    idir_list3.append(idir)
                    isubdir_list3.append(sdir)
                    icovmat_list3.append(icovmat)
                    
        # - - - - - -
        print(f"  Prep index lists for {n_job_tot} COVMAT jobs")

        self.config_prep['n_job_split'] = 1
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        self.config_prep['idir_list3']     = idir_list3
        self.config_prep['isubdir_list3']  = isubdir_list3
        self.config_prep['icovmat_list3']  = icovmat_list3

 
        return
        # end prep_covmat_index_lists

    def write_command_file(self, icpu, f):

        n_core             = self.config_prep['n_core']
        bbc_outdir_dict    = self.config_prep['bbc_outdir_dict']
        covmatopt_dict     = self.config_prep['covmatopt_dict'] 

        bbc_outdir_list    = bbc_outdir_dict['outdir_list']  
        bbc_subdir_list2   = bbc_outdir_dict['subdir_list2'] 
        arg_list           = covmatopt_dict['jobopt_arg_list']

        idir_list3    = self.config_prep['idir_list3'] 
        isubdir_list3 = self.config_prep['isubdir_list3']
        icovmat_list3 = self.config_prep['icovmat_list3']
        
        n_job_cpu   = 0
        n_job_local = 0
        
        for idir,isubdir,icovmat in \
            zip(idir_list3,isubdir_list3,icovmat_list3):

            n_job_local += 1
            index_dict = { 'idir':idir, 'isubdir':isubdir, 
                           'icovmat':icovmat, 'icpu':icpu }

            if ( (n_job_local-1) % n_core ) != icpu : continue

            n_job_cpu += 1
            job_info_wfit   = self.prep_JOB_INFO_covmat(index_dict)
            util.write_job_info(f, job_info_wfit, icpu)

            #job_info_merge = \
            #    self.prep_JOB_INFO_merge(icpu,n_job_local,False) 
            #util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu
        # end write_command_file

    def prep_JOB_INFO_covmat(self,index_dict):

        idir    = index_dict['idir']
        isubdir = index_dict['isubdir']
        icovmat = index_dict['icovmat']

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

        prefix        = f"COVMAT_{job_label}"
        log_file      = f"{prefix}.LOG" 
        done_file     = f"{prefix}.DONE"
        all_done_file = f"{output_dir}/{DEFAULT_DONE_FILE}"

        covmat_outdir = f"{output_dir}/{job_label}"

        # start with user-defined args from WFITOPT[_GLOBAL] key
        arg_list = []
        arg_list.append(f"--input_dir {bbc_outdir}")
        arg_list.append(f"--version   {bbc_subdir}")
        arg_list.append(f"--outdir    {covmat_outdir}")
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


        # .xyz END
