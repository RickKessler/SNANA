# Jan 23 2021
# SALT3 training using saltshaker code from arxiv:xxxx
# 
# TODO:
#  + check that $CONDA_DEFAULT_ENV = salt3 before launching
#  + tack on --yamloutputfile YAMLOUTPUTFILE

import  os, sys, shutil, yaml, glob
import  logging, coloredlogs
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *
from    submit_prog_base import Program

# define yaml keys for input files
CONFIG_KEYLIST_INPUT_FILE = [ 'INPUT_TRAIN_FILE', 'INPUT_MODEL_FILE' ]

# define input file keys passed to trainsalt code
CODE_KEYLIST_INPUT_FILE   = [ '--configfile', '--trainingconfig' ]

TRAINOPT_STRING = "TRAINOPT"

# Define suffix for output model used by LC fitters: 
#    SALT2.[MODEL_SUFFIX][nnn]
# Default output dirs are SALT3.MODEL000, SALT3.MODEL001, ...
MODEL_SUFFIX_DEFAULT = "MODEL"

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================


class train_SALT3(Program):
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

        # get input config files
        self.train_prep_input_files()

        # scoop up TRAINOPT list from user CONFIG
        self.train_prep_trainopt_list()

        # foreach training, prepare output paths 
        self.train_prep_paths()

        # copy input file
        self.train_prep_copy_files()

        # end submit_prepare_driver

    def train_prep_input_files(self):

        # read & store names of input files for trainsalt

        input_master_file = self.config_yaml['args'].input_file 
        CONFIG            = self.config_yaml['CONFIG']
        
        input_file_list = []  # init input files for trainsalt
        msgerr = []

        for key in CONFIG_KEYLIST_INPUT_FILE:
            if key not in CONFIG:
                msgerr.append(f"Missing required yaml key '{key}' ")
                msgerr.append(f"in {input_master_file} ")
                self.log_assert(False,msgerr)
            else:
                input_file = CONFIG[key]
                input_file_list.append(input_file)
                
        self.config_prep['input_file_list'] = input_file_list

        # end train_prep_input_files

    def train_prep_trainopt_list(self):

        CONFIG   = self.config_yaml['CONFIG']
        key      = TRAINOPT_STRING
        if key in CONFIG :
            trainopt_rows = CONFIG[key]
        else:
            trainopt_rows = []

        # - - - - - 
        trainopt_dict = util.prep_jobopt_list(trainopt_rows, 
                                              TRAINOPT_STRING, 
                                              None )

        n_trainopt          = trainopt_dict['n_jobopt']
        trainopt_arg_list   = trainopt_dict['jobopt_arg_list']
        trainopt_ARG_list   = trainopt_dict['jobopt_ARG_list']
        trainopt_num_list   = trainopt_dict['jobopt_num_list']  
        trainopt_label_list = trainopt_dict['jobopt_label_list']
        trainopt_shift_file = trainopt_dict['jobopt_file_list']
        use_arg_file        = trainopt_dict['use_arg_file']

        logging.info(f" Store {n_trainopt-1} TRAIN-SALT3 options " \
                     f"from {TRAINOPT_STRING} keys")

        self.config_prep['n_trainopt']          = n_trainopt
        self.config_prep['trainopt_arg_list']   = trainopt_arg_list
        self.config_prep['trainopt_ARG_list']   = trainopt_ARG_list
        self.config_prep['trainopt_num_list']   = trainopt_num_list
        self.config_prep['trainopt_label_list'] = trainopt_label_list
        self.config_prep['trainopt_shift_file'] = trainopt_shift_file
        self.config_prep['use_arg_file']        = use_arg_file

        # end train_prep_trainopt_list

    def train_prep_paths(self):

        # for each TRAINOPT, create path for SALTPATH and for training output

        output_dir        = self.config_prep['output_dir']
        trainopt_num_list = self.config_prep['trainopt_num_list']
        model_suffix      = MODEL_SUFFIX_DEFAULT  
        outdir_model_list  = []

        for trainopt_num in trainopt_num_list:
            nnn           = trainopt_num[-3:]
            sdir_model    = f"SALT3.{model_suffix}{nnn}"
            outdir_model  = f"{output_dir}/{sdir_model}"
            print(f"\t Create {sdir_model}")
            outdir_model_list.append(outdir_model)
            os.mkdir(outdir_model)

        self.config_prep['outdir_model_list']  =  outdir_model_list
        # end train_prep_paths

    def train_prep_copy_files(self):

        input_file_list = self.config_prep['input_file_list']
        script_dir      = self.config_prep['script_dir']

        for input_file in input_file_list:
            os.system(f"cp {input_file} {script_dir}/")

        # end train_prep_copy_files

    def write_command_file(self, icpu, f):
        # For this icpu, write full set of sim commands to
        # already-opened command file with pointer f. 
        # Function returns number of jobs for this cpu

        n_core          = self.config_prep['n_core']
        n_trainopt      = self.config_prep['n_trainopt'] 
        n_job_tot   = n_trainopt
        n_job_split = 1     # cannot break up train job
        n_job_local = 0
        n_job_cpu   = 0

        self.config_prep['n_job_split'] = n_job_split
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        for itrain in range(0,n_trainopt):
            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :

                n_job_cpu += 1
                job_info_train   = self.prep_JOB_INFO_train(itrain)
                util.write_job_info(f, job_info_train, icpu)
    
                job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_local) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu

        # end write_command_file

    def prep_JOB_INFO_train(self,itrain):

        CONFIG            = self.config_yaml['CONFIG']
        program           = self.config_prep['program']
        input_file_list   = self.config_prep['input_file_list']
        script_dir        = self.config_prep['script_dir']
        kill_on_fail      = self.config_yaml['args'].kill_on_fail

        output_dir        = self.config_prep['output_dir']
        trainopt_num      = self.config_prep['trainopt_num_list'][itrain]
        trainopt_arg      = self.config_prep['trainopt_arg_list'][itrain]
        outdir_model      = self.config_prep['outdir_model_list'][itrain]

        prefix            = trainopt_num
        arg_list          = [ ]
        msgerr            = [ ]

        log_file   = f"{prefix}.LOG"
        done_file  = f"{prefix}.DONE"
        start_file = f"{prefix}.START"
        yaml_file  = f"{prefix}.YAML"

        for key,input_file in zip(CODE_KEYLIST_INPUT_FILE,input_file_list):
            arg_list.append(f"{key} {input_file}")

        arg_list.append(f"--outputdir {outdir_model}")
        arg_list.append(f"--yamloutputfile {yaml_file}")
        arg_list.append(f"{trainopt_arg}")

        JOB_INFO = {}
        JOB_INFO['program']       = f"{program}"
        JOB_INFO['input_file']    = ""  
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = log_file
        JOB_INFO['done_file']     = done_file
        JOB_INFO['start_file']    = start_file
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['arg_list']      = arg_list

        return JOB_INFO

        # end prep_JOB_INFO_train

    def create_merge_table(self,f):

        n_trainopt        = self.config_prep['n_trainopt'] 
        trainopt_num_list = self.config_prep['trainopt_num_list']

        header_line_merge = (f"    STATE   {TRAINOPT_STRING}  NLC NSPEC  CPU")
        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for num in trainopt_num_list :
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(num)     # e..g, TRAINOPT002
            ROW_MERGE.append(0)       # NLC
            ROW_MERGE.append(0)       # NSPEC
            ROW_MERGE.append(0.0)     # CPU

            INFO_MERGE['row_list'].append(ROW_MERGE)  
        # - - - -
        util.write_merge_file(f, INFO_MERGE, [] ) 

        # end create_merge_table

    def append_info_file(self,f):

        # append info to SUBMIT.INFO file; use passed file pointer f

        n_trainopt   = self.config_prep['n_trainopt'] 
        num_list     = self.config_prep['trainopt_num_list']
        arg_list     = self.config_prep['trainopt_arg_list'] 
        ARG_list     = self.config_prep['trainopt_ARG_list'] 
        label_list   = self.config_prep['trainopt_label_list']

        f.write(f"# train_SALT2 info \n")
        f.write(f"JOBFILE_WILDCARD: {TRAINOPT_STRING}* \n")

        f.write(f"\n")
        f.write(f"TRAINOPT_OUT_LIST:  " \
                f"# 'TRAINOPTNUM'  'user_label'  'user_args'\n")
        # use original ARG_list instead of arg_list; the latter may
        # include contents of shiftlist_file.
        for num, arg, label in zip(num_list, ARG_list, label_list):
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")
        
       
        # .xyz
