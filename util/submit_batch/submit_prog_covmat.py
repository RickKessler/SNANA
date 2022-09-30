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

        # read yaml input from create_cov input file
        self.prep_covmat_read_inputs()

        self.prep_covmat_find_bbc_outdirs()

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
        bbc_rawdir_list = CONFIG['BBC_OUTDIR']

        # use util.prep_jobopt_list to strip labels and args

        COLNUM_BBC_MERGE_VERSION   = 1
        bbc_outdir_list  = []  # list of BBC outdirs
        bbc_subdir_list2 = []  # each BBC outdir can have subdir list
                               # e.g., 50 independent sim samples

        for outdir in bbc_rawdir_list:
            outdir = os.path.expandvars(outdir)
            bbc_outdir_list.append(outdir)
            merge_log = f"{outdir}/{MERGE_LOG_FILE}"            
            merge_yaml, comment_lines = util.read_merge_file(merge_log) 
            
            bbc_subdir_list = []
            for row in merge_yaml['MERGE']:
                version = row[COLNUM_BBC_MERGE_VERSION]
                if version not in bbc_outdir_list:
                    bbc_subdir_list.append(version)
            bbc_subdir_list2.append(bbc_subdir_list)

        print(f"\n xxx bbc_outdir_list = \n{bbc_outdir_list}\n")
        self.config_prep['bbc_outdir_list']  = bbc_outdir_list
        self.config_prep['bbc_subdir_list1'] = bbc_outdir_list2

        return
        # end prep_covmat_find_bbc_outdirs


        # .xyz END
