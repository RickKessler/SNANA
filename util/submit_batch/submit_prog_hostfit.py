import os, sys, shutil, yaml, glob, logging
import datetime, time, subprocess
import submit_util as util
import pandas as pd

from   submit_params import *
from   submit_prog_base import Program

SUBCLASS_HOSTFIT_CIGALE = 'cigale'

# Created May 2026
# run SED fit code on galaxies to get host properties.
#
class HostPropertyFit(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = None
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

        logging.info(f" xxx hello from new task HostPropertyFit: outdir = {output_dir_name} ")
        return output_dir_name, SUBDIR_SCRIPTS_HOSTFIT

        
    def submit_prepare_driver(self):
        input_file = self.config_yaml['args'].input_file
        SNANA_TO_CIGALE = util.extract_yaml(input_file, "SNANA_TO_CIGALE:", None)

        CIGALE_TO_SNANA = util.extract_yaml(input_file, "CIGALE_TO_SNANA:", None)

        if SNANA_TO_CIGALE:
            SUBCLASS = SUBCLASS_HOSTFIT_CIGALE
        else:
            sys.exit('Cannot find subclass for HOSTFIT')

        self.config_prep['n_job_tot'] = self.config_prep['n_core']
        self.config_prep['n_job_split'] = self.config_prep['n_core']
        self.config_prep['SUBCLASS'] = SUBCLASS
        self.config_prep['SNANA_TO_CIGALE'] = SNANA_TO_CIGALE
        self.config_prep['CIGALE_TO_SNANA'] = CIGALE_TO_SNANA
        logging.info(f'SUBCLASS = {SUBCLASS}')

        #CONFIG_HOSTFIT     = self.config_yaml['SNANA_TO_CIGALE']
        #FILTERS       = CONFIG.setdefault(KEY_FILTERS, None)
        #sys.exit(self.config_yaml)
        return

    def write_command_file(self, icpu, f):
        n_core = self.config_prep['n_core']
        print('xxx n_core = ', n_core)

        for ijob in range(n_core):
            job_info_hostfit   = self.prep_JOB_INFO_hostfit(ijob)
            util.write_job_info(f, job_info_hostfit, icpu)


        sys.exit('xxx bye')
        return

    def prep_JOB_INFO_hostfit(self, ijob):
        SNANA_TO_CIGALE = self.config_prep['SNANA_TO_CIGALE']
        program       = self.config_prep['program']
        output_dir    = self.config_prep['output_dir']
        script_dir    = self.config_prep['script_dir']
        input_file    = args.input_file
        n_job_split   = self.config_prep['n_job_split']
        SUBCLASS = self.config_prep['SUBCLASS']
        

        split_num     = f"SPLIT{isplit:03d}"
        prefix        = f"{SUBCLASS}_{split_num}"
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
