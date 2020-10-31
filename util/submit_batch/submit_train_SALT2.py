# Created Oct 29 2020

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time, subprocess
import submit_util as util
from   submit_params    import *
from   submit_prog_base import Program

SALTPATH_STRING = "SALTPATH"
SALTPATH_FILES  = [ 'fitmodel.card',  'Instruments',  'MagSys' ]

SUBDIR_TRAIN = "TRAIN"
SUBDIR_MODEL = "MODEL"

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================


class train_SALT2(Program):
    def __init__(self, config_yaml) :
        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_UNKNOWN
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

        # scoop up TRAINOPT list from user CONFIG
        self.train_prep_trainopt_list()

        # foreach training, prepare path for SALTPATH and for training output
        self.train_prep_paths()

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
        
    def train_prep_paths(self):

        # for each TRAINOPT, create path for SALTPATH and for training output

        output_dir        = self.config_prep['output_dir']
        trainopt_num_list = self.config_prep['trainopt_num_list']
        trainopt_arg_list = self.config_prep['trainopt_arg_list']
        
        outdir_saltpath_list   = []
        outdir_train_list      = []

        topdir_saltpath = (f"{output_dir}/{SALTPATH_STRING}")
        os.mkdir(topdir_saltpath)

        for num,arg in zip(trainopt_num_list,trainopt_arg_list) :
            outdir_saltpath    = (f"{topdir_saltpath}/{num}")
            outdir_train       = (f"{output_dir}/{num}")
            outdir_saltpath_list.append(outdir_saltpath)
            outdir_train_list.append(outdir_train)
            os.mkdir(outdir_saltpath)
            os.mkdir(outdir_train)
            self.train_prep_SALTPATH(outdir_saltpath,num,arg)
            
        self.config_prep['outdir_saltpath_list']  =  outdir_saltpath_list
        self.config_prep['outdir_train_list']     =  outdir_train_list
        
            # .xyz
        # end train_prep_paths

    def train_prep_SALTPATH(self,outdir_saltpath,num,arg):

        CONFIG           = self.config_yaml['CONFIG']
        SALTPATH_BASE    = CONFIG['SALTPATH_BASE']

        logging.info(f"\t Prepare {SALTPATH_STRING}/{num}")
        for item in SALTPATH_FILES :
            cmd_rsync = (f"rsync -r {SALTPATH_BASE}/{item} {outdir_saltpath}")
            os.system(cmd_rsync)

        # insert code here from Georgie/Chris to modify outdir_saltpath
        # based on user {arg} passed from input file.


        # end train_prep_SALTPATH

    def write_command_file(self, icpu, COMMAND_FILE):

        n_core          = self.config_prep['n_core']
        n_trainopt      = self.config_prep['n_trainopt']            
        n_job_tot   = n_trainopt
        n_job_split = 1     # cannot break up train job
        n_job_local = 0

        # open CMD file for this icpu  
        f = open(COMMAND_FILE, 'a')

        for itrain in range(0,n_trainopt):
            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :

                job_info_train   = self.prep_JOB_INFO_train(itrain)
                util.write_job_info(f, job_info_train, icpu)
    
        f.close()
            
        self.config_prep['n_job_split'] = n_job_split
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        # end write_command_file

    def prep_JOB_INFO_train(self,itrain):

        input_file        = self.config_yaml['args'].input_file 
        program           = self.config_prep['program']
        script_dir        = self.config_prep['script_dir']
        kill_on_fail      = self.config_yaml['args'].kill_on_fail

        output_dir        = self.config_prep['output_dir']
        trainopt_num_list = self.config_prep['trainopt_num_list']
        trainopt_arg_list = self.config_prep['trainopt_arg_list']
        saltpath_list     = self.config_prep['outdir_saltpath_list']
        saltpath          = saltpath_list[itrain]
        prefix            = trainopt_num_list[itrain]
        arg_list          = [ ]

        trainDir_file  = (f"{prefix}.CONFIG")
        self.create_trainDir_file(itrain,trainDir_file)

        JOB_INFO = {}
        JOB_INFO['setenv']        = (f"export {SALTPATH_STRING}={saltpath}")
        JOB_INFO['program']       = program
        JOB_INFO['input_file']    = trainDir_file
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = (f"{prefix}.LOG")
        JOB_INFO['done_file']     = (f"{prefix}.DONE")
        JOB_INFO['all_done_file'] = (f"{output_dir}/{DEFAULT_DONE_FILE}")
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['arg_list']      = arg_list
        
        return JOB_INFO

        # end prep_JOB_INFO_train

    def create_trainDir_file(self,itrain,trainDir_file):

        # create input config file that is passed to shell training script

        CONFIG            = self.config_yaml['CONFIG']
        script_dir        = self.config_prep['script_dir']
        trainpath_list    = self.config_prep['outdir_train_list']

        trainpath         = trainpath_list[itrain]
        TRAINDIR_FILE     = (f"{script_dir}/{trainDir_file}")
        SALT2_MODEL_BASE  = CONFIG['SALT2_MODEL_BASE']

        with open(TRAINDIR_FILE,"wt") as f :
            f.write(f"initDir:     {SALT2_MODEL_BASE} \n")
            f.write(f"trainingDir: {trainpath}/{SUBDIR_TRAIN} \n")
            f.write(f"outputDir:   {trainpath}/{SUBDIR_MODEL} \n")

        # end create_trainDir_file

    def create_merge_table(self,f):

        n_trainopt        = self.config_prep['n_trainopt'] 
        trainopt_num_list = self.config_prep['trainopt_num_list']

        header_line_merge = (f" STATE   TRAINOPT  NEVT_DATA ")
        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for num in trainopt_num_list :
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(num)     # e..g, TRAINOPT002
            ROW_MERGE.append(0)       # NEVT_DATA

            INFO_MERGE['row_list'].append(ROW_MERGE)  
        # - - - -
        util.write_merge_file(f, INFO_MERGE, [] ) 

        # end create_merge_table

    def append_info_file(self,f):

        # add info to SUBMIT.INFO file; use passed file pointer f

        n_trainopt = self.config_prep['n_trainopt'] 
        num_list   = self.config_prep['trainopt_num_list']
        arg_list   = self.config_prep['trainopt_arg_list'] 
        label_list = self.config_prep['trainopt_label_list']

        f.write("\n")
        f.write("TRAINOPT_OUT_LIST:  " \
                "# 'TRAINOPTNUM'  'user_label'  'user_args'\n")
        for num,arg,label in zip(num_list,arg_list,label_list):
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")
        
        # end append_info_file

