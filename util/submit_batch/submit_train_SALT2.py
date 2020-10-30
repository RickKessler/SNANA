# Created Oct 29 2020

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time, subprocess
import submit_util as util
from   submit_params    import *
from   submit_prog_base import Program


# ====================================================
#    BEGIN FUNCTIONS
# ====================================================


class train_SALT2(Program):
    def __init__(self, config_yaml) :
        config_prep = {}
        config_prep['program'] = "undefined"
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

        return output_dir_name, SUBDIR_SCRIPTS_TRAIN
        # end set_output_dir_name

    def submit_prepare_driver(self):

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file 

        self.train_prep_trainopt_list()

        # end submit_prepare_driver

    def train_prep_trainopt_list(self):
        CONFIG           = self.config_yaml['CONFIG']
        input_file       = self.config_yaml['args'].input_file 
        n_trainopt          = 1
        trainopt_arg_list   = [ '' ] 
        trainopt_num_list   = [ 'TRAINOPT000' ] 
        trainopt_label_list = [ None ]
        
        key = 'TRAINOPT'
        if key in CONFIG  :
            for trainopt_raw in CONFIG[key] : # might include label
                num = (f"TRAINOPT{n_trainopt:03d}")
                label, trainopt = util.separate_label_from_arg(trainopt_raw)
                trainopt_arg_list.append(trainopt)
                trainopt_num_list.append(num)
                trainopt_label_list.append(label)
                n_trainopt += 1
                
        logging.info(f" Store {n_trainopt-1} TRAIN-SALT2 options " \
                     "from TRAINOPT keys")

        self.config_prep['n_trainopt']          = n_trainopt
        self.config_prep['trainopt_arg_list']   = trainopt_arg_list
        self.config_prep['trainopt_num_list']   = trainopt_num_list
        self.config_prep['trainopt_label_list'] = trainopt_label_list

        # end train_prep_trainopt_list
        

    def write_command_file(self, icpu, COMMAND_FILE):
        pass

