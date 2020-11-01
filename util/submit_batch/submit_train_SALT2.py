# Created Oct 29 2020

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time, subprocess
import submit_util as util
from   submit_params    import *
from   submit_prog_base import Program

TRAINOPT_STRING = "TRAINOPT"
SALTPATH_STRING = "SALTPATH"
SALTPATH_FILES  = [ 'fitmodel.card',  'Instruments',  'MagSys' ]

KEY_MODEL_SUFFIX     = "MODEL_SUFFIX"
MODEL_SUFFIX_DEFAULT = "TRAIN"

SUBDIR_TRAIN = "TRAIN"
SUBDIR_MODEL = "MODEL"

COLNUM_TRAIN_MERGE_TRAINOPT       = 1
COLNUM_TRAIN_MERGE_NEVT           = 2
COLNUM_TRAIN_MERGE_CPU            = 3

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
        trainopt_num_list   = [ f"{TRAINOPT_STRING}000" ] 
        trainopt_label_list = [ None ]
                
        key = TRAINOPT_STRING
        if key in CONFIG  :
            for trainopt_raw in CONFIG[key] : # might include label
                num = (f"{TRAINOPT_STRING}{n_trainopt:03d}")
                label, trainopt = util.separate_label_from_arg(trainopt_raw)
                trainopt_arg_list.append(trainopt)
                trainopt_num_list.append(num)
                trainopt_label_list.append(label)
                n_trainopt += 1
            
        logging.info(f" Store {n_trainopt-1} TRAIN-SALT2 options " \
                     f"from {TRAINOPT_STRING} keys")

        self.config_prep['n_trainopt']          = n_trainopt
        self.config_prep['trainopt_arg_list']   = trainopt_arg_list
        self.config_prep['trainopt_num_list']   = trainopt_num_list
        self.config_prep['trainopt_label_list'] = trainopt_label_list

        # end train_prep_trainopt_list
        
    def get_path_trainopt(self,which,trainopt):
        output_dir    = self.config_prep['output_dir']
        model_suffix  = self.config_prep['model_suffix']
        nnn           = trainopt[-3:]
        path          = ""
        # print(f" xxx trainopt={trainopt}  nnn={nnn}")
        if which == SALTPATH_STRING :
            path = (f"{output_dir}/{SALTPATH_STRING}/{trainopt}")
        elif which == SUBDIR_MODEL :
            path = (f"{output_dir}/SALT2.{model_suffix}{nnn}")
        elif which == SUBDIR_TRAIN :
            path = (f"{output_dir}/{SUBDIR_TRAIN}/{trainopt}")
            
        return path
        # end get_path_trainopt

    def train_prep_paths(self):

        # for each TRAINOPT, create path for SALTPATH and for training output

        output_dir        = self.config_prep['output_dir']
        trainopt_num_list = self.config_prep['trainopt_num_list']
        trainopt_arg_list = self.config_prep['trainopt_arg_list']
        CONFIG            = self.config_yaml['CONFIG']

        if KEY_MODEL_SUFFIX in CONFIG : 
            model_suffix = CONFIG[KEY_MODEL_SUFFIX] 
        else:
            model_suffix = MODEL_SUFFIX_DEFAULT
        self.config_prep['model_suffix']  = model_suffix 
        
        outdir_saltpath_list   = []
        outdir_train_list      = []
        outdir_model_list      = []

        topdir_saltpath = (f"{output_dir}/{SALTPATH_STRING}")
        topdir_train    = (f"{output_dir}/{SUBDIR_TRAIN}")
        os.mkdir(topdir_saltpath)
        os.mkdir(topdir_train)

        for trainopt,arg in zip(trainopt_num_list,trainopt_arg_list) :
            outdir_saltpath = self.get_path_trainopt(SALTPATH_STRING,trainopt)
            outdir_train    = self.get_path_trainopt(SUBDIR_TRAIN,trainopt)
            outdir_model    = self.get_path_trainopt(SUBDIR_MODEL,trainopt)
            outdir_saltpath_list.append(outdir_saltpath)
            outdir_train_list.append(outdir_train)
            outdir_model_list.append(outdir_model)
            os.mkdir(outdir_saltpath)
            os.mkdir(outdir_train)
            os.mkdir(outdir_model)
            self.train_prep_SALTPATH(outdir_saltpath,trainopt,arg)
            
        self.config_prep['outdir_saltpath_list']  =  outdir_saltpath_list
        self.config_prep['outdir_train_list']     =  outdir_train_list
        self.config_prep['outdir_model_list']     =  outdir_model_list

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

        self.config_prep['n_job_split'] = n_job_split
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        # open CMD file for this icpu  
        f = open(COMMAND_FILE, 'a')

        for itrain in range(0,n_trainopt):
            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :

                job_info_train   = self.prep_JOB_INFO_train(itrain)
                util.write_job_info(f, job_info_train, icpu)
    
                job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_local) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

        f.close()            

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
        JOB_INFO['program']       = (f"python {program}")
        JOB_INFO['input_file']    = "" 
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = (f"{prefix}.LOG")
        JOB_INFO['done_file']     = (f"{prefix}.DONE")
        JOB_INFO['all_done_file'] = (f"{output_dir}/{DEFAULT_DONE_FILE}")
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['arg_list']      = arg_list
        
        arg_list.append(f"-c {trainDir_file}")

        return JOB_INFO

        # end prep_JOB_INFO_train

    def create_trainDir_file(self,itrain,trainDir_file):

        # create input config file (trainDir_file) that is passed 
        # to shell training script. Config file contains directory
        # names for output.
        
        CONFIG            = self.config_yaml['CONFIG']
        script_dir        = self.config_prep['script_dir']
        output_dir        = self.config_prep['output_dir']

        #outdir_saltpath_list = self.config_prep['outdir_saltpath_list'] 
        outdir_train_list    = self.config_prep['outdir_train_list']
        outdir_model_list    = self.config_prep['outdir_model_list']

        PATH_SALT2_BASE   = os.path.expandvars(CONFIG['PATH_SALT2_BASE'])
        outdir_train      = outdir_train_list[itrain]
        outdir_model      = outdir_model_list[itrain]

        # config file passed to pythin training shell
        TRAINDIR_FILE     = (f"{script_dir}/{trainDir_file}")

        with open(TRAINDIR_FILE, "wt") as f :
            f.write(f"initDir      {PATH_SALT2_BASE} \n")
            f.write(f"trainingDir  {outdir_train} \n")
            f.write(f"outputDir    {outdir_model} \n")

        # end create_trainDir_file

    def create_merge_table(self,f):

        n_trainopt        = self.config_prep['n_trainopt'] 
        trainopt_num_list = self.config_prep['trainopt_num_list']

        header_line_merge = (f"    STATE   {TRAINOPT_STRING}  NEVT  CPU")
        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for num in trainopt_num_list :
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(num)     # e..g, TRAINOPT002
            ROW_MERGE.append(0)       # NEVT
            ROW_MERGE.append(0.0)     # CPU

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

        f.write(f"# train_SALT2 info \n")
        f.write(f"JOBFILE_WILDCARD: {TRAINOPT_STRING}* \n")
        f.write(f"MODEL_SUFFIX: {model_suffix} \n")

        f.write(f"\n")
        f.write(f"TRAINOPT_OUT_LIST:  " \
                f"# 'TRAINOPTNUM'  'user_label'  'user_args'\n")
        for num,arg,label in zip(num_list,arg_list,label_list):
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")
        
        # end append_info_file

    def merge_config_prep(self,output_dir):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        model_suffix     = submit_info_yaml['MODEL_SUFFIX']
        self.config_prep['output_dir']     = output_dir 
        self.config_prep['script_dir']     = script_dir 
        self.config_prep['model_suffix']   = model_suffix
        # end merge_config_prep

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_TRAINOPT  = COLNUM_TRAIN_MERGE_TRAINOPT  
        COLNUM_NEVT      = COLNUM_TRAIN_MERGE_NEVT
        COLNUM_CPU       = COLNUM_TRAIN_MERGE_CPU

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []
        row_list_merge     = MERGE_INFO_CONTENTS[TABLE_MERGE]

        nrow_check = 0
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input
            nrow_check += 1
            irow        = nrow_check - 1 # row index
            trainopt    = row[COLNUM_TRAINOPT] # e.g., TRAINOPT001
            search_wildcard = (f"{trainopt}*")

            # strip off row info
            STATE       = row[COLNUM_STATE]

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

                    # since there is no YAML file to examine, we have a 
                    # kluge check on success
                    success,tproc = self.get_train_status(trainopt)
                    if not success : 
                        self.check_for_failure(log_list[0], -1, +1)
                        NEW_STATE = SUBMIT_STATE_FAIL

                    row[COLNUM_STATE]     = NEW_STATE
                    row[COLNUM_NEVT]      = 0  # ??? fill this later
                    row[COLNUM_CPU]       = tproc

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             

        # - - - - - -  -
        # The first return arg (row_split) is null since there is 
        # no need for a SPLIT table
        return [], row_list_merge_new, n_state_change

        # end merge_update_state

    def get_train_status(self,trainopt):

        # for input trainopt (e.g., TRAINOPT001), check if SALT2
        # training succeeded (returns True) or failed (returns False).
        # Initial check (Halloween, 2020) is existance of
        #   output_dir/trainopt/MODEL/salt2_template_0.dat
        #   output_dir/trainopt/MODEL/salt2_template_1.dat
        #
        # Function returns success(True or False),tproc (process time, minutes)

        #submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = self.config_prep['script_dir']

        model_dir  = self.get_path_trainopt(SUBDIR_MODEL,trainopt)
        train_dir  = self.get_path_trainopt(SUBDIR_TRAIN,trainopt)

        s0_file    = (f"{model_dir}/salt2_template_0.dat")
        s1_file    = (f"{model_dir}/salt2_template_1.dat")

        # get process time between original AAA_README and DONE file
        done_file  = (f"{script_dir}/{trainopt}.DONE")
        start_file = (f"{train_dir}/AAA_README")
        tdone      = os.path.getmtime(done_file)
        tstart     = os.path.getmtime(start_file)
        tproc      = int((tdone - tstart)/60.0)

        if not os.path.exists(s0_file) : return False,tproc
        if not os.path.exists(s1_file) : return False,tproc

        # make sure each file has something in it
        num_lines_s0 = sum(1 for line in open(s0_file))
        num_lines_s1 = sum(1 for line in open(s1_file))
        if num_lines_s0 < 100 : return False,tproc
        if num_lines_s1 < 100 : return False,tproc

        return True,tproc
        # end get_train_status

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):

        # gzip and tar for 'irow' training job.
        output_dir       = self.config_prep['output_dir']
        script_dir       = self.config_prep['script_dir']

        row      = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        trainopt = row[COLNUM_TRAIN_MERGE_TRAINOPT]

        model_dir = self.get_path_trainopt(SUBDIR_MODEL,trainopt)
        train_dir = self.get_path_trainopt(SUBDIR_TRAIN,trainopt)

        logging.info(f"    Compress output for {trainopt} :")

        # make tar file from SALTPATH
        logging.info(f"\t Compress {SALTPATH_STRING}/{trainopt}")
        cmd_cddir = (f"cd {output_dir}/{SALTPATH_STRING}")
        cmd_tar   = (f"tar -cf {trainopt}.tar {trainopt}; gzip {trainopt}.tar")
        cmd_rm    = (f"rm -r {trainopt}")
        cmd       = (f"{cmd_cddir}; {cmd_tar} ; {cmd_rm}")
        os.system(cmd)

        # Gzip contents of TRAIN and TRAIN -> TRAIN.tar
        logging.info(f"\t Compress {SUBDIR_TRAIN}/{trainopt}")
        cmd_rmfits = (f"cd {train_dir}; rm *.fits")
        cmd_gzip   = (f"cd {train_dir}; gzip *.dat *.list")
        cmd_tar    = (f"cd {train_dir}/../ ; " \
                      f"tar -cf {trainopt}.tar {trainopt} ; " \
                      f"rm -r {trainopt}" )
        os.system(cmd_rmfits)
        os.system(cmd_gzip)
        os.system(cmd_tar)

        # gzip contents of MODEL, leave directory
        logging.info(f"\t gzip contents of {model_dir}")
        cmd = (f"cd {model_dir}; gzip salt2*.dat")
        os.system(cmd)

        # end merge_job_wrapup

    def merge_cleanup_final(self):

        script_dir       = self.config_prep['script_dir']

        # start with script_dir; separate tar files for CPU* and TRAINOPT*
        prefix = (f"CPU")
        util.compress_files(+1, script_dir, f"{prefix}*", prefix, "" )

        prefix = (f"{TRAINOPT_STRING}")
        util.compress_files(+1, script_dir, f"{prefix}*", prefix, "" )

        # - - - - - -
        # maybe later, SALTPATH -> SALTPATH.tar and TRAIN -> TRAIN.tar

        # end merge_cleanup_final

    def get_misc_merge_info(self):
        # return misc info to append in MERGE.LOG file
        misc_info = []
        return misc_info

    def get_merge_COLNUM_CPU(self):
        return COLNUM_TRAIN_MERGE_CPU

# === END: ===
