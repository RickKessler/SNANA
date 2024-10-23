# Created May 2023 by R.Kessler
# 
import  os, sys, shutil, yaml, configparser, glob
import  logging
#import  coloredlogs
import  datetime, time, subprocess
import  submit_util as util
import  submit_train_util  as train_util

from    submit_train_util  import *
from    submit_params    import *
from    submit_prog_base import Program

# define key for main trainsalt config file
KEY_CONFIG_FILE = 'BAYESN_CONFIG_FILE'

# define keys that point to other files to copy
KEY_LIST_FILE_COPY = []

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================


class train_BAYESN(Program):
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
        script_dir   = self.config_prep['script_dir']
        output_dir   = self.config_prep['output_dir']

        # scoop up and store TRAINOPT list from user CONFIG.
        # This is before prepping input files in case TRAINOPT
        # have additional input files.
        config_prep_trainopt = \
            train_util.train_prep_trainopt_list(METHOD_TRAIN_SALT3,
                                                CONFIG, script_dir)
        self.config_prep.update(config_prep_trainopt)
 

        # get input config files
        self.train_prep_input_files()

        # copy input files to script_dir
        input_file_list = self.train_prep_input_files()
        self.config_prep['input_file_list'] = input_file_list
        util.copy_input_files(input_file_list, script_dir, TRAIN_INPUT_FILENAMES)

        # foreach training, prepare output paths 
        trainopt_num_list = self.config_prep['trainopt_num_list']
        outdir_model_list = \
            train_util.prep_model_paths(METHOD_TRAIN_BAYESN, trainopt_num_list,
                        output_dir )
        self.config_prep['outdir_model_list'] = outdir_model_list

        sys.stdout.flush()

        # end submit_prepare_driver

    def train_prep_input_files(self):

        # read & store names of input files for trainsalt

        input_master_file = self.config_yaml['args'].input_file 
        CONFIG            = self.config_yaml['CONFIG']
        
        # start list of all input files for trainsalt
        input_file_list = [ ]  
        msgerr = []

        if KEY_CONFIG_FILE not in CONFIG:
            msgerr.append(f"Missig required key {KEY_CONFIG_FILE} ")
            self.log_assert(False,msgerr)

        config_file = CONFIG[KEY_CONFIG_FILE]
        input_file_list.append(config_file)

        for input_file in input_file_list:
            if not os.path.exists(input_file):
                msgerr.append(f"Input file {input_file}")
                msgerr.append(f"does not exist.")
                msgerr.append(f"Check '{key}' arg in {input_master_file}")
                util.log_assert(False,msgerr) # just abort, no done stamp
                
        
        input_file_list += [ input_master_file ]

        # return list of all input files
        return input_file_list

        # end train_prep_input_files

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
    
                job_info_merge = \
                    self.prep_JOB_INFO_merge(icpu,n_job_local,False) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu

        # end write_command_file

    def prep_JOB_INFO_train(self,itrain):

        CONFIG            = self.config_yaml['CONFIG']
        program           = self.config_prep['program']
        input_file_list   = self.config_prep['input_file_list']
        config_file       = self.config_prep['input_file_list'][0]
        script_dir        = self.config_prep['script_dir']
        kill_on_fail      = self.config_yaml['args'].kill_on_fail

        output_dir        = self.config_prep['output_dir']
        trainopt_num      = self.config_prep['trainopt_num_list'][itrain]
        trainopt_arg      = self.config_prep['trainopt_arg_list'][itrain]
        trainopt_global   = self.config_prep['trainopt_global']
        outdir_model      = self.config_prep['outdir_model_list'][itrain]
        do_fast           = self.config_yaml['args'].fast
        do_global         = len(trainopt_global) > 0

        prefix            = trainopt_num
        arg_list          = [ ]
        msgerr            = [ ]

        log_file   = f"{prefix}.LOG"
        done_file  = f"{prefix}.DONE"
        start_file = f"{prefix}.START"
        yaml_file  = f"{prefix}.YAML"

        arg_list.append(f"--input {config_file}")
        arg_list.append(f"--outputdir {outdir_model}")
        arg_list.append(f"--yamloutputfile {yaml_file}")
        arg_list.append(f"{trainopt_arg}")
        if do_global : arg_list.append(trainopt_global)
        if do_fast   : arg_list.append("--fast")

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

        append_info_dict = {
            'METHOD'      : METHOD_TRAIN_BAYESN,
            'CONFIG'      : self.config_yaml['CONFIG'],
            'n_trainopt'  : self.config_prep['n_trainopt'],
            'num_list'    : self.config_prep['trainopt_num_list'],
            'arg_list'    : self.config_prep['trainopt_arg_list'],
            'ARG_list'    : self.config_prep['trainopt_ARG_list'],
            'label_list'  : self.config_prep['trainopt_label_list'],
            'calib_shift_list'  : self.config_prep['calib_shift_list'],
            'outdir_model_list' : self.config_prep['outdir_model_list']
        }

        train_util.append_info_file(f, append_info_dict)
        return

        # end append_info_file

    def merge_config_prep(self,output_dir):
        pass

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.
        # Oct 10 2022: fix indentation bug under 'if NDONE == n_job_split'
        #               --> fixes phony failures.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_TRAINOPT  = COLNUM_TRAIN_MERGE_TRAINOPT  
        COLNUM_NLC       = COLNUM_TRAIN_MERGE_NLC
        COLNUM_NSPEC     = COLNUM_TRAIN_MERGE_NSPEC
        COLNUM_CPU       = COLNUM_TRAIN_MERGE_CPU

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []
        row_list_merge     = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # keynames_for_job_stats returns 3 keynames : 
        #   {base}, {base}_sum, {base}_list
        key_nlc, key_nlc_sum, key_nlc_list = \
                self.keynames_for_job_stats('NUM_SNE')
        key_nspec, key_nspec_sum, key_nspec_list = \
                 self.keynames_for_job_stats('NUM_SPECTRA')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')
        key_list = [ key_nlc, key_nspec, key_cpu ] 

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

                    job_stats = self.get_job_stats(script_dir, 
                                                   log_list, 
                                                   yaml_list, 
                                                   key_list)

                    row[COLNUM_STATE]     = NEW_STATE
                    row[COLNUM_NLC]       = job_stats[key_nlc_sum]
                    row[COLNUM_NSPEC]     = job_stats[key_nspec_sum]
                    row[COLNUM_CPU]       = job_stats[key_cpu_sum]

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             
                
        # - - - - - 
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

        # cleanup for 'irow' training job.

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        modeldir_list    = submit_info_yaml['MODELDIR_LIST']

        row      = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        trainopt = row[COLNUM_TRAIN_MERGE_TRAINOPT]
        model_dir = f"{output_dir}/{modeldir_list[irow]}"
        misc_dir  = f"{model_dir}/{SUBDIR_MISC}"
        print(f" Cleanup {model_dir}")

        # - - - - - -
        # check if SALTshaker already cleaned up
        tar_file = f"{SUBDIR_MISC}.tar"
        tmp_list = glob.glob1(misc_dir,f"{tar_file}*")
        done_misc = len(tmp_list) > 0

        if not done_misc :
            os.mkdir(misc_dir)
            mv_string = "*.png *.pdf *.npy gauss* *parameters* " \
                        "salt3train_snparams.txt"
            cmd = f"cd {model_dir} ; mv {mv_string} {SUBDIR_MISC}/"
            os.system(cmd)

            # tar and gzip misc sub dir
            cmd_tar = f"cd {model_dir} ; " \
                    f"tar -cf {tar_file} {SUBDIR_MISC} ; " \
                    f"gzip  {tar_file} ; " \
                    f"rm -r {SUBDIR_MISC} "
            os.system(cmd_tar)

        # gzip dat files
        cmd_gzip = f"cd {model_dir} ; gzip *.dat "
        os.system(cmd_gzip)
        
        # append SALT3.INFO file with calib_shift info
        self.append_SALT3_INFO_FILE(trainopt,model_dir)

        # end  merge_job_wrapup

    def append_SALT3_INFO_FILE(self,trainopt,model_dir):

        # append SALT2.INFO file with calib info;
        # used by SNANA light curve fitting code to use same shift
        # as in the training.

        submit_info_yaml = self.config_prep['submit_info_yaml']

        SURVEY_LIST_SAMEMAGSYS = submit_info_yaml['SURVEY_LIST_SAMEMAGSYS']
        SURVEY_LIST_SAMEFILTER = submit_info_yaml['SURVEY_LIST_SAMEFILTER']

        if trainopt == f"{TRAINOPT_STRING}000" : return

        # read SNANA_SALT3_INFO from SUBMIT.INFO; this includes
        # info for all trainopts.
        SNANA_INFO         = submit_info_yaml[KEY_SNANA_SALT3_INFO]

        SALT3_INFO_FILE = f"{model_dir}/SALT3.INFO"
        print(f"\t Append calib info to {SALT3_INFO_FILE}")

        f = open(SALT3_INFO_FILE,"at")
        f.write(f"\n\n# Calibration shifts used in SALTshaker Training\n")

        for row in SNANA_INFO:
            if row[0] == trainopt: 
                key    = row[1]
                survey = row[2]
                band   = row[3]
                shift  = row[4]
                f.write(f"{key}: {survey:<10} {band} {shift} \n")

                # write other surveys with same magsys (Jul 9 2021)
                if key == KEY_MAGSHIFT:
                    s_list  = SURVEY_LIST_SAMEMAGSYS
                else:
                    s_list  = SURVEY_LIST_SAMEFILTER

                if survey in s_list :
                    for s in filter(lambda s: s not in [survey], s_list):
                        comment = f"same {key} as {survey}"
                        f.write(f"{key}: {s:<10} {band} {shift}  " \
                                f"  # {comment}\n")
        f.close()

        # end append_SALT2_INFO_FILE

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file.
        # Each info line must be of the form
        #  KEYNAME:  VALUE

        return []
        # end get_misc_merge_info

    def merge_cleanup_final(self):
        # every snlc_fit job succeeded, so here we simply compress output.
  
        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        wildcard_list = [ 'TRAINOPT', 'CPU', 'CALIB_SHIFT' ]
        for w in wildcard_list :
            wstar = f"{w}*"
            tmp_list = glob.glob1(script_dir,wstar)
            if len(tmp_list) == 0 : continue
            print(f"\t Compress {wstar}")
            util.compress_files(+1, script_dir, wstar, w, "" )

        # - - - -
        # tar up entire script dir
        util.compress_subdir(+1, script_dir)

        # end merge_cleanup_final

    def get_merge_COLNUM_CPU(self):
        return COLNUM_TRAIN_MERGE_CPU

# =========== END: =======



