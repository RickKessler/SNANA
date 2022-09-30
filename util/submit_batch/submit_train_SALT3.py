# Jan 23 2021
# SALT3 training using saltshaker code from 
#    https://arxiv.org/abs/2104.07795
# 
# Jul 09 2021: 
#   new keys 'SURVEY_LIST_SAMEMAGSYS and 'SURVEY_LIST_SAMEFILTER'  
#   to record extra surveys in SALT2.INFO file.      


import  os, sys, shutil, yaml, configparser, glob
import  logging, coloredlogs
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *
from    submit_prog_base import Program

# define key for main trainsalt config file
KEY_CONFIG_FILE = 'SALT3_CONFIG_FILE'

TRAINOPT_STRING        = "TRAINOPT"
TRAINOPT_GLOBAL_STRING = "TRAINOPT_GLOBAL"

# Define suffix for output model used by LC fitters: 
#    SALT2.[MODEL_SUFFIX][nnn]
# Default output dirs are SALT3.MODEL000, SALT3.MODEL001, ...
MODEL_SUFFIX_DEFAULT = "MODEL"

# Define columns in MERGE.LOG. Column 0 is always the STATE.                   
COLNUM_TRAIN_MERGE_TRAINOPT    = 1
COLNUM_TRAIN_MERGE_NLC         = 2
COLNUM_TRAIN_MERGE_NSPEC       = 3
COLNUM_TRAIN_MERGE_CPU         = 4


# config keys for calibration shifts (same as for train_SALT2)
KEY_MAGSHIFT       = "MAGSHIFT"
KEY_WAVESHIFT      = "WAVESHIFT"
KEY_LAMSHIFT       = "LAMSHIFT"
KEY_SHIFTLIST_FILE = "SHIFTLIST_FILE"
KEY_CALIBSHIFT_LIST  = [ KEY_MAGSHIFT, KEY_WAVESHIFT, KEY_LAMSHIFT ]

KEYS_SURVEY_LIST_SAME = ['SURVEY_LIST_SAMEMAGSYS', 'SURVEY_LIST_SAMEFILTER']

# define prefix for files with calib shifts.
PREFIX_CALIB_SHIFT   = "CALIB_SHIFT"  

# define command-line override key to specify file with calibration shifts
KEY_SALTshaker_CALIBSHIFT_FILE = "--calibrationshiftfile"

# define key for SALT3.INFO read by SNANA codes.
KEY_SNANA_SALT3_INFO = "SNANA_SALT3_INFO"

# create list of config keys whose argument is a file that
# gets copied to script_dir
SECTION_FILE_COPY  = 'iodata'
KEY_LIST_FILE_COPY = [ 'trainingconfig', 'tmaxlist', 'snparlist' ]

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

        # scoop up and store TRAINOPT list from user CONFIG.
        # This is before prepping input files in case TRAINOPT
        # have additional input files.
        self.train_prep_trainopt_list()

        # get input config files
        self.train_prep_input_files()

        # copy input files to script_dir
        self.train_prep_copy_files()

        # scoop up TRAINOPT list from user CONFIG
        # xxxx  self.train_prep_trainopt_list()

        # foreach training, prepare output paths 
        self.train_prep_paths()

        sys.stdout.flush()

        # end submit_prepare_driver

    def train_prep_input_files(self):

        # read & store names of input files for trainsalt

        input_master_file = self.config_yaml['args'].input_file 
        CONFIG            = self.config_yaml['CONFIG']
        
        input_file_list = []  # init input files for trainsalt
        msgerr = []

        if KEY_CONFIG_FILE not in CONFIG:
            msgerr.append(f"Missig required key {KEY_CONFIG_FILE} ")
            self.log_assert(False,msgerr)

        config_file = CONFIG[KEY_CONFIG_FILE]
        input_file_list.append(config_file)

        print(f"\n Check other input files inside main config file: " \
              f"{config_file}")

        # parse config_file to get other potential input files
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read(config_file)

        for key in KEY_LIST_FILE_COPY :
            if key in config[SECTION_FILE_COPY]:
                input_file = config[SECTION_FILE_COPY][key]
                input_file_list.append(input_file)

        for input_file in input_file_list:
            if not os.path.exists(input_file):
                msgerr.append(f"Input file {input_file}")
                msgerr.append(f"does not exist.")
                msgerr.append(f"Check '{key}' arg in {input_master_file}")
                util.log_assert(False,msgerr) # just abort, no done stamp
                
        # - - - - - 
        # check for additional input files in the TRAINOPT
        
        trainopt_global   = self.config_prep['trainopt_global']
        trainopt_arg_list = self.config_prep['trainopt_arg_list']
        trainopt_all      = trainopt_arg_list + [ trainopt_global ]

        for item in trainopt_all:
            item_list = item.split()
            for key in KEY_LIST_FILE_COPY :
                key_override = f"--{key}"
                if key_override in item_list:
                    j = item_list.index(key_override)                    
                    input_file = item_list[j+1]
                    input_file_list.append(input_file)

        #sys.exit(f"\n xxx trainopt_all = \n{trainopt_all}")
        
        # store list of all input files
        self.config_prep['input_file_list'] = input_file_list

        # end train_prep_input_files

    def train_prep_trainopt_list(self):

        CONFIG   = self.config_yaml['CONFIG']

        trainopt_global = ""  # apply to each TRAINOPT
        trainopt_rows   = []  # TRAINOT-specified commands

        # start with global settings
        key = TRAINOPT_GLOBAL_STRING
        if key in CONFIG:  trainopt_global = CONFIG[key]

        # next, TRAINOPT per job
        key      = TRAINOPT_STRING
        if key in CONFIG :  trainopt_rows = CONFIG[key]

        # - - - - - 
        trainopt_dict = util.prep_jobopt_list(trainopt_rows, 
                                              TRAINOPT_STRING, 1,
                                              KEY_SHIFTLIST_FILE )

        n_trainopt          = trainopt_dict['n_jobopt']
        trainopt_arg_list   = trainopt_dict['jobopt_arg_list']
        trainopt_ARG_list   = trainopt_dict['jobopt_ARG_list']
        trainopt_num_list   = trainopt_dict['jobopt_num_list']  
        trainopt_label_list = trainopt_dict['jobopt_label_list']
        trainopt_shift_file = trainopt_dict['jobopt_file_list']
        use_arg_file        = trainopt_dict['use_arg_file']

        logging.info(f" Store {n_trainopt-1} TRAIN-SALT3 options " \
                     f"from {TRAINOPT_STRING} keys")

        # for MAGSHIFT and./or WAVESHIFT, create calibration file
        # in script_dir. The returned calib_shift_file has each
        # calib arg unpacked so that each row has only 1 shift.
        # arg_replace = arg, but calibration shifts are replaced
        # with command to read calib-shift file.
        arg_replace_list = []
        calib_shift_list = []
        for num,arg in zip(trainopt_num_list,trainopt_arg_list):
            arg_replace, calib_shift = self.make_calib_shift_file(num,arg)
            arg_replace_list.append(arg_replace)
            calib_shift_list.append(calib_shift)  # list of lists

        self.config_prep['n_trainopt']          = n_trainopt
        self.config_prep['trainopt_arg_list']   = arg_replace_list
        self.config_prep['trainopt_ARG_list']   = trainopt_ARG_list
        self.config_prep['trainopt_num_list']   = trainopt_num_list
        self.config_prep['trainopt_label_list'] = trainopt_label_list
        self.config_prep['trainopt_shift_file'] = trainopt_shift_file
        self.config_prep['trainopt_global']     = trainopt_global
        self.config_prep['use_arg_file']        = use_arg_file
        self.config_prep['calib_shift_list']    = calib_shift_list

        print('')

        # end train_prep_trainopt_list

    def make_calib_shift_file(self,num,arg):

        # if arg contains MAGSHIFT or WAVESHIFT key, write them out
        # calib-shift file.
        # Example:
        #  num = TRAINOPT003
        #  arg = WAVESHIFT CfA3  r,i 10,8     MAGSHIFT CfA3 U .01
        #
        # ==> write the following  to CALIB_SHIFT_TRAINOPT003.DAT:
        #
        # WAVESHIFT CfA3  r 10
        # WAVESHIFT CfA3  i  8
        # MAGSHFIT  cfA3  U 0.01
        #
        # and function returns arg_replace = 
        #    "calibrationshiftsfile = CALIB_SHIFT_TRAINOPT003.DAT"
        # and also returns calib_shift for SUBMIT.INFO file
        #
        # Make sure that arg_replace retains non-calib optoins
        # to allow mixing calib and non-calib arguments.

        script_dir      = self.config_prep['script_dir']

        # first check if any calib shift key is in this arg
        found_calshift = False
        for key in KEY_CALIBSHIFT_LIST :
            if key in arg: found_calshift = True
        if not found_calshift : return arg, []

        # if we get here, there is at least one valid calib-shift key,
        # so unpack arg and write calib-shift file for SALTshaker.

        calib_shift_file = f"{PREFIX_CALIB_SHIFT}_{num}.DAT"
        CALIB_SHIFT_FILE = f"{script_dir}/{calib_shift_file}"
        f = open(CALIB_SHIFT_FILE,"wt")

        print(f"\t Create {calib_shift_file}")
        arg_list = arg.split()
        arg_replace = "" 
        n_arg       = len(arg_list)
        use_list    = [ False ] * n_arg
        calib_shift_list = []

        for iarg in range(0,n_arg):            
            item = arg_list[iarg]
            if item not in KEY_CALIBSHIFT_LIST :  continue 
            key       = arg_list[iarg+0]
            survey    = arg_list[iarg+1]
            band_list = arg_list[iarg+2].split(',')
            val_list  = arg_list[iarg+3].split(',')
            use_list[iarg:iarg+4] = [True] * 4

            for band,val in zip(band_list,val_list):
                line = f"{key} {survey} {band} {val}"
                f.write(f"{line}\n")
                calib_shift_list.append(line)

        f.close()

        # - - - - 
        # store un-used arguments in arg_replace
        for item,use in zip(arg_list,use_list):
            if not use : arg_replace += f"{item} " 

        # tack on calibshift file
        arg_replace += f"{KEY_SALTshaker_CALIBSHIFT_FILE} {calib_shift_file} "

        return arg_replace, calib_shift_list

        # end make_calshift_file

    def train_prep_paths(self):

        # for each TRAINOPT, create path for SALTPATH and for training output

        output_dir        = self.config_prep['output_dir']
        trainopt_num_list = self.config_prep['trainopt_num_list']
        n_trainopt        = self.config_prep['n_trainopt']
        model_suffix      = MODEL_SUFFIX_DEFAULT  
        outdir_model_list  = []
        outdir_model_list_base  = []

        print(f" Create {n_trainopt} output model directories:")
        for trainopt_num in trainopt_num_list:
            nnn           = trainopt_num[-3:]
            sdir_model    = f"SALT3.{model_suffix}{nnn}"
            outdir_model  = f"{output_dir}/{sdir_model}"
            print(f"\t Create {sdir_model}")
            outdir_model_list.append(outdir_model)
            outdir_model_list_base.append(sdir_model)
            os.mkdir(outdir_model)

        sys.stdout.flush()
        self.config_prep['outdir_model_list']  =  outdir_model_list
        self.config_prep['outdir_model_list_base'] = outdir_model_list_base
        # end train_prep_paths

    def train_prep_copy_files(self):

        # copy input files to script_dir, but only those which don't
        # have full path ... to avoid clobbering same input tile names
        # in different directories

        input_file_list = self.config_prep['input_file_list']
        script_dir      = self.config_prep['script_dir']

        for input_file in input_file_list:
            if '/' not in input_file:
                print(f"\t Copy input file: {input_file}")
                os.system(f"cp {input_file} {script_dir}/")

        print('')
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

        # xxxx mark delete 
        #for key,input_file in zip(CODE_KEYLIST_INPUT_FILE,input_file_list):
        #    arg_list.append(f"{key} {input_file}")
        # xxxx

        arg_list.append(f"--configfile {config_file}")
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

        CONFIG       = self.config_yaml['CONFIG']
        n_trainopt   = self.config_prep['n_trainopt'] 
        num_list     = self.config_prep['trainopt_num_list']
        arg_list     = self.config_prep['trainopt_arg_list'] 
        ARG_list     = self.config_prep['trainopt_ARG_list'] 
        label_list   = self.config_prep['trainopt_label_list']
        calib_shift_list = self.config_prep['calib_shift_list']
        outdir_model_list_base = self.config_prep['outdir_model_list_base']

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

        f.write("MODELDIR_LIST:\n")
        for model_dir in outdir_model_list_base :
            f.write(f"  - {model_dir}\n")
        f.write("\n")

        for key in KEYS_SURVEY_LIST_SAME:
            if key in CONFIG :
                f.write(f"{key}:  {CONFIG[key]} \n")
            else:
                f.write(f"{key}:  [ ] \n")

        # write keys for SALT2.INFO to be read by SNANA code 
        # each row is
        #  [ 'TRAINOPTnnn', KEY, SURVEY, SHIFT_VAL ]
        f.write(f"{KEY_SNANA_SALT3_INFO}: \n")
        for num, item_list in zip(num_list,calib_shift_list) :
            for item in item_list:
                row = [ num ] + item.split()
                f.write(f"  - {row} \n")
        f.write("\n")

        # end append_info_file

    def merge_config_prep(self,output_dir):
        pass

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

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
                                               log_list, yaml_list, key_list)

                row[COLNUM_STATE]     = NEW_STATE
                row[COLNUM_NLC]       = job_stats[key_nlc_sum]
                row[COLNUM_NSPEC]     = job_stats[key_nspec_sum]
                row[COLNUM_CPU]       = job_stats[key_cpu_sum]

                row_list_merge_new[irow] = row  # update new row
                n_state_change += 1             

        # first return arg (row_split) is null since there is 
        # no need for a SPLIT table
        row_list_dict = {
            'row_split_list' : [],
            'row_merge_list' : row_list_merge_new,
            'row_extra_list' : []
        }
        return row_list_dict, n_state_change
        # xxx mark return [], row_list_merge_new, n_state_change

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



