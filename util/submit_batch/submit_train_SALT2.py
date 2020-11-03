# Created Oct 31 2020
#
# Remaining issues:
#   * how to specify and implement TRAINOPT in CONFIG block
#   * how are MAG_OFFSET and SIGMA_INT determined ??
#   * can train_SALT2_37_mw.py  produce YAML output for submit_batch ?
#       NEVT:           348
#       ABORT_IF_ZERO:  348
#       CPU_MINUTES:    137.3
#       ANYTHING_ELSE:  ??
#   * how is SUCCESS determined ? Existence of salt2_template_[01].dat ?
#   * fast option ?
# 
import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time, subprocess
import submit_util as util
from   submit_params    import *
from   submit_prog_base import Program

TRAINOPT_STRING = "TRAINOPT"
SALTPATH_STRING = "SALTPATH"  # env name of calib path for snpca
SALTPATH_FILES  = [ 'fitmodel.card',  'Instruments',  'MagSys' ]

# define subdirs to contain outputs of snpca that are not used by LC fitters
SUBDIR_CALIB    = "CALIB_TRAIN"    # a.k.a SALTPATH
SUBDIR_TRAIN    = "OUTPUT_TRAIN"

# define subdirs under SUBDIR_CALIB ($SALTPATH)
SUBDIR_MAGSYS   = "MagSys"
SUBDIR_INSTR    = "Instruments"

# define keys for CONFIG inputs
KEY_PATH_INPUT_TRAIN = "PATH_INPUT_TRAIN"  # required (data,conf files)
KEY_PATH_INPUT_CALIB = "PATH_INPUT_CALIB"  # required (a.k.a, SALTPATH)
KEY_MODEL_SUFFIX     = "MODEL_SUFFIX"      # optional: change MODEL suffix

# Define suffix for output model used by LC fitters.
# Default output dirs are SALT2.MODEL000, SALT2.MODEL001, ...
MODEL_SUFFIX_DEFAULT = "MODEL"        

# Define list of trained SALT2 model files to check for existence;
# their existence defines SUCCESS/FAIL of training job.
COLORLAW_FILE   = "salt2_color_correction.dat"
CHECK_FILE_LIST = [ "salt2_template_0.dat", "salt2_template_1.dat", 
                    COLORLAW_FILE ]

# define content for SALT2_INFO file
SALT2_INFO_FILE    = "SALT2.INFO"    # created for SNANA programs
SALT2_MAG_OFFSET   = 0.27
SALT2_SIGMA_INT    = 0.106 
SALT2_COLOR_OFFSET = 0.0

SALT2_INFO_INCLUDE = f"""
SEDFLUX_INTERP_OPT: 2  # 1=>linear,    2=>spline
ERRMAP_INTERP_OPT:  1  # 0=snake off;  1=>linear  2=>spline
ERRMAP_KCOR_OPT:    1  # 1/0 => on/off

MAGERR_FLOOR:   0.005            # model-error floor
MAGERR_LAMOBS:  0.0  2000  4000  # magerr minlam maxlam
MAGERR_LAMREST: 0.1   100   200  # magerr minlam maxlam
"""

# Define columns in MERGE.LOG. Column 0 is always the STATE.
COLNUM_TRAIN_MERGE_TRAINOPT    = 1
COLNUM_TRAIN_MERGE_NEVT        = 2
COLNUM_TRAIN_MERGE_CPU         = 3

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
        path          = ""
        # print(f" xxx trainopt={trainopt}  nnn={nnn}")
        if which == SUBDIR_CALIB :
            path = (f"{output_dir}/{SUBDIR_CALIB}/{trainopt}")
        elif which == SUBDIR_TRAIN :
            path = (f"{output_dir}/{SUBDIR_TRAIN}/{trainopt}")
        elif which == "MODEL" :
            model_suffix  = self.config_prep['model_suffix']
            nnn           = trainopt[-3:]
            path = (f"{output_dir}/SALT2.{model_suffix}{nnn}")
            
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
        
        outdir_calib_list  = []
        outdir_train_list  = []
        outdir_model_list  = []

        topdir_calib    = (f"{output_dir}/{SUBDIR_CALIB}")
        topdir_train    = (f"{output_dir}/{SUBDIR_TRAIN}")
        os.mkdir(topdir_calib)
        os.mkdir(topdir_train)

        for trainopt,arg in zip(trainopt_num_list,trainopt_arg_list) :
            outdir_calib   = self.get_path_trainopt(SUBDIR_CALIB,trainopt)
            outdir_train   = self.get_path_trainopt(SUBDIR_TRAIN,trainopt)
            outdir_model   = self.get_path_trainopt("MODEL",trainopt)
            outdir_calib_list.append(outdir_calib)
            outdir_train_list.append(outdir_train)
            outdir_model_list.append(outdir_model)
            os.mkdir(outdir_calib)
            os.mkdir(outdir_train)
            os.mkdir(outdir_model)
            self.train_prep_SALTPATH(outdir_calib,trainopt,arg)
            
        self.config_prep['outdir_calib_list']  =  outdir_calib_list
        self.config_prep['outdir_train_list']  =  outdir_train_list
        self.config_prep['outdir_model_list']  =  outdir_model_list

        # end train_prep_paths

    def train_prep_SALTPATH(self,outdir_calib,num,arg):

        # prepare calibration file(s) pointed to by ENV $SALTPATH.
        # First copy base calibration files, then check input arg
        # (from TRAINOPT input) for instructions on modifying 
        # calibration file(s).

        CONFIG           = self.config_yaml['CONFIG']
        PATH_INPUT_CALIB = CONFIG[KEY_PATH_INPUT_CALIB] # aka SALTPATH
        msgerr = []

        logging.info(f"\t Prepare {SUBDIR_CALIB}/{num}")
        for item in SALTPATH_FILES :
            cmd_rsync = (f"rsync -r {PATH_INPUT_CALIB}/{item} {outdir_calib}")
            os.system(cmd_rsync)

        if arg == '' : return

        #return # xxxx REMOVE

        if 'MAGSHIFT' in arg :
            mag_file, mag_arg_list, mag_shift = \
                    self.parse_MAGSHIFT(outdir_calib,arg)
            nchange = \
                self.update_file_shift_mag(mag_file,mag_arg_list,mag_shift)
            if nchange != 1:
                msgerr.append(f"Applied {nchange} mag changes, but expect 1.")
                msgerr.append(f"Check user TRAINOPT '{arg}' ")
                msgerr.append(f"and check mag file")
                msgerr.append(f"   {mag_file}")
                self.log_assert(False,msgerr)

        elif 'WAVESHIFT' in arg :
            filter_file, wave_shift = \
                    self.parse_WAVESHIFT(outdir_calib,arg)
            status = self.update_file_shift_wave(filter_file, wave_shift)
            if status != 0 :
                msgerr.append(f"Unable to open filter-transmissino file")
                msgerr.append(f"  {filter_file}")
                self.log_assert(False,msgerr)            
        else :
            pass

        # end train_prep_SALTPATH

    def parse_MAGSHIFT(self,outdir_calib,arg):
        # parse TRAINOPT arg and return 
        #  1) mag file name to modify
        #  2) list of keys to identify row to modify
        #  3) mag shift to add

        arg_split    = arg.split()
        survey       = arg_split[1]
        band         = arg_split[2]
        mag_shift    = float(arg_split[3])
        mag_arg_list = [ survey, band ]

        # use map to convert survey/band into file .xyz
        mag_file = (f"{outdir_calib}/{SUBDIR_MAGSYS}/SDSS-AB-off.dat")

        print(f"\n XXX TEST CHANGE IN MAG-FILE {mag_file}")
        return mag_file, mag_arg_list, mag_shift
        # end parse_MAGSHIFT

    def parse_WAVESHIFT(self,outdir_calib,arg):
        # parse TRAINOPT arg and return 
        #  1) filter-transmission file name to modify
        #  2) wave shift.

        arg_split    = arg.split()
        survey       = arg_split[1]
        band         = arg_split[2]
        wave_shift   = float(arg_split[3])

        # use map to convert survey/band into file name .xyz
        filter_file  = outdir_calib 
        filter_file += (f"/{SUBDIR_INSTR}/Keplercam/{band}_Keplercam.txt")

        print(f"\n XXX TEST CHANGE IN FILTER-FILE {filter_file}")
        return filter_file, wave_shift
        # end parse_WAVESHIFT

    def update_file_shift_wave(self,filter_file, wave_shift):
        # modify filter_file to have wavelength shift of wave_shift. 
        if not os.path.exists(filter_file): return -1
        msgerr     = []
        line_list_out = []
        line_list_inp = []
        line_list_out.append(f"# wave_shift = {wave_shift} A has been " \
                             "applied by\n")
        line_list_out.append(f"#    {sys.argv[0]}\n")

        # scoop up lines in file
        with open(filter_file,"rt") as f:
            for line in f:
                line_list_inp.append(line)

        # shift value in first column
        for line in line_list_inp:
            word_list = line.split()
            line_list_out.append(line) # default new line = old line

            if line[0] == '#' : continue
            word_list = (line.rstrip("\n")).split()
            if len(word_list) < 2 : continue
            wave_orig = float(word_list[0])
            wave_out  = wave_orig + wave_shift
            new_line  = str(wave_out) + "  " + " ".join(word_list[1:])
            line_list_out[-1:] = (f" {new_line} \n")

        # - - - - - - - -
        with open(filter_file,"wt") as f:
            f.write("".join(line_list_out))

        return 0
        # end update_file_shift_wave

    def update_file_shift_mag(self, file_name, key_list, mag_shift ):

        # for input file_name, search for key_list sequence and 
        # shift value by mag_shift.
        # Example: key_list = [ SDSS, u ] and mag_shift = 0.01
        # If file_name includes a line with
        #    SDSS u 0.022
        # then modify the file line to be
        #    SDSS u 0.032  #  0.022 + 0.01
        #
        # F
        # If file_name does not exist, return -1.
        #
        # Perhaps move to submit_util?

        if not os.path.exists(file_name): return -1

        nchange    = 0
        nkey       = len(key_list)  # number of key strings to match
        key_string = " ".join(key_list)
        msgerr     = []
        line_list_out = []
        line_list_inp = []
        line_list_out.append(f"# This file has been modified by\n")
        line_list_out.append(f"#    {sys.argv[0]}\n\n")

        # scoop up lines in file
        with open(file_name,"rt") as f:
            for line in f:
                line_list_inp.append(line)

        # look for line(s) to modify
        for line in line_list_inp:
            line_list_out.append(line) # default new line = old line
            nmatch_key = 0
            if line[0] == '#' : continue
            word_list = (line.rstrip("\n")).split()
            if len(word_list) < nkey + 1 : continue

            for i in range(0,nkey):
                if key_list[i] == word_list[i] : nmatch_key += 1
            if nmatch_key == nkey:
                strval = word_list[nkey] 
                try :
                    value = float(strval)
                except:
                    msgerr.append(f"Expected float after '{key_string}'" )
                    msgerr.append(f"but found {strval}")
                    msgerr.append(f"Check file: {file_name}")
                    self.log_assert(False,msgerr)

                new_value  = value + mag_shift
                new_line   = (f"{key_string} {new_value:.6f}" \
                              f"   # {value} + {mag_shift} \n")
                line_list_out[-1:] = new_line  # overwrite last line_out
                nchange += 1

        # modify file with line_list_out
        if nchange > 0 :
            with open(file_name,"wt") as f:
                f.write("".join(line_list_out))

        return nchange
        # end shift_float_value
    
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
        path_calib_list   = self.config_prep['outdir_calib_list']
        path_calib        = path_calib_list[itrain]
        prefix            = trainopt_num_list[itrain]
        arg_list          = [ ]

        trainDir_file  = (f"{prefix}.CONFIG")
        self.create_trainDir_file(itrain,trainDir_file)

        JOB_INFO = {}
        JOB_INFO['setenv']        = (f"export {SALTPATH_STRING}={path_calib}")
        JOB_INFO['program']       = (f"python {program}")
        JOB_INFO['input_file']    = "" 
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = (f"{prefix}.LOG")
        JOB_INFO['done_file']     = (f"{prefix}.DONE")
        JOB_INFO['start_file']    = (f"{prefix}.START")
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

        PATH_INPUT_TRAIN  = os.path.expandvars(CONFIG[KEY_PATH_INPUT_TRAIN])
        outdir_train      = outdir_train_list[itrain]
        outdir_model      = outdir_model_list[itrain]

        # config file passed to pythin training shell
        TRAINDIR_FILE     = (f"{script_dir}/{trainDir_file}")

        with open(TRAINDIR_FILE, "wt") as f :
            f.write(f"initDir      {PATH_INPUT_TRAIN} \n")
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

        # append info to SUBMIT.INFO file; use passed file pointer f

        n_trainopt   = self.config_prep['n_trainopt'] 
        num_list     = self.config_prep['trainopt_num_list']
        arg_list     = self.config_prep['trainopt_arg_list'] 
        label_list   = self.config_prep['trainopt_label_list']
        model_suffix = self.config_prep['model_suffix']

        f.write(f"# train_SALT2 info \n")
        f.write(f"JOBFILE_WILDCARD: {TRAINOPT_STRING}* \n")
        f.write(f"MODEL_SUFFIX: {model_suffix}   " \
                f"# -> create SALT2.{model_suffix}nnn/ \n")

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

        # For input trainopt (e.g., TRAINOPT001), check if SALT2
        # training succeeds and return status = True or False.
        # Initial check (Halloween, 2020) is existance of files in
        # global CHECK_FILE_LIST array.
        #
        # Also get CPU process time (tproc) based on START and DONE files.
        # Function returns 
        #     status,tproc

        output_dir   = self.config_prep['output_dir']
        script_dir   = self.config_prep['script_dir']
        model_dir    = self.get_path_trainopt("MODEL",trainopt)
        train_dir    = self.get_path_trainopt(SUBDIR_TRAIN,trainopt)

        # get process time between START and DONE files
        done_file  = (f"{script_dir}/{trainopt}.DONE")
        start_file = (f"{script_dir}/{trainopt}.START")
        tdone      = os.path.getmtime(done_file)
        tstart     = os.path.getmtime(start_file)
        tproc      = int((tdone - tstart)/60.0)

        # check for existence of SALT2 model files
        nerr = 0
        for check_file in CHECK_FILE_LIST:
            CHECK_FILE = (f"{model_dir}/{check_file}")
            if os.path.exists(CHECK_FILE) : 
                # make sure each file has something in it
                num_lines = sum(1 for line in open(CHECK_FILE))
                if num_lines < 3 : 
                    nerr += 1
                    logging.info(f" ERROR: only {num_lines} in {CHECK_FILE}")
            else:
                nerr += 1
                logging.info(f" ERROR: missing {CHECK_FILE}")

        # - - - - - - - - -
        status = (nerr == 0)
        return status,tproc
        # end get_train_status

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):

        # gzip and tar for 'irow' training job.
        output_dir       = self.config_prep['output_dir']
        script_dir       = self.config_prep['script_dir']

        row      = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        trainopt = row[COLNUM_TRAIN_MERGE_TRAINOPT]

        self.merge_create_SALT2_INFO_file(trainopt)

        model_dir  = self.get_path_trainopt("MODEL",trainopt)
        train_dir  = self.get_path_trainopt(SUBDIR_TRAIN,trainopt)
        calib_dir  = self.get_path_trainopt(SUBDIR_CALIB,trainopt)

        logging.info(f"    Compress output for {trainopt} :")

        # make tar file from CALIB/TRAINOPTnnn  (aka SALTPATH)
        logging.info(f"\t Compress {SUBDIR_CALIB}/{trainopt}")
        util.compress_subdir(+1,calib_dir)

        # Gzip contents of TRAINOPT, then  TRAINOPTnnn -> TRAINOPTnnn.tar.gz
        logging.info(f"\t Compress {SUBDIR_TRAIN}/{trainopt}")
        cmd_clean = (f"cd {train_dir}; rm *.fits; gzip *.dat *.list")
        os.system(cmd_clean)
        util.compress_subdir(+1,train_dir)

        # gzip contents of MODEL, leave directory intact for LC fitter
        logging.info(f"\t gzip contents of {model_dir}")
        cmd = (f"cd {model_dir}; gzip salt2*.dat")
        os.system(cmd)

        # end merge_job_wrapup

    def merge_create_SALT2_INFO_file(self,trainopt):
        # create SALT2.INFO file read by SNANA programs
        # (snlc_sim.exe & snlc_fit.exe)

        model_dir    = self.get_path_trainopt("MODEL",trainopt)
        info_file    = (f"{model_dir}/{SALT2_INFO_FILE}")

        color_law_dict = self.get_color_law(model_dir)
        min_lambda   = color_law_dict['min_lambda']
        max_lambda   = color_law_dict['max_lambda']
        range_lambda = (f"{min_lambda} {max_lambda}")
        npar         = color_law_dict['npar']
        par_list     = " ".join(color_law_dict['par_list'])

        print(f"    Create {SALT2_INFO_FILE}")
        with open(info_file,"wt") as f:

            f.write(f"RESTLAMBDA_RANGE: {range_lambda}\n")
            f.write(f"COLORLAW_VERSION: {color_law_dict['version']} \n")
            f.write(f"COLORCOR_PARAMS:  {range_lambda} {npar} {par_list}\n")
            f.write(f"\n")

            f.write(f"MAG_OFFSET:    {SALT2_MAG_OFFSET} \n")
            f.write(f"SIGMA_INT:     {SALT2_SIGMA_INT}  \n")
            f.write(f"COLOR_OFFSET:  {SALT2_COLOR_OFFSET} \n")

            f.write(f"{SALT2_INFO_INCLUDE}\n")

            f.write(f"# {trainopt} \n")
            f.write(f"# ?? calib option keys for SNANA ... ?? \n")
            
        # #end merge_create_SALT2_INFO_file
        
    def get_color_law(self,model_dir):
        file_name = (f"{model_dir}/{COLORLAW_FILE}")
        color_law_dict = {}

        nline = 0; npar_tot = -1; npar_read=0; cpar_list = []

        with open(file_name,"r") as f:
            for line in f:
                nline += 1
                word_list = (line.rstrip("\n")).split()
                #print(f" xxx colorlaw line = {word_list} ")
                if nline == 1:
                    npar_tot = int(word_list[0])
                elif npar_read < npar_tot :
                    npar_read += 1
                    cpar_list.append(word_list[0])
                elif 'version' in word_list[0] :
                    version = int(word_list[1])
                elif 'min_lambda' in word_list[0] :
                    min_lambda = int(word_list[1])
                elif 'max_lambda' in word_list[0] :
                    max_lambda = int(word_list[1])
        # - - - - - 
        color_law_dict['npar']        = npar_tot
        color_law_dict['par_list']    = cpar_list
        color_law_dict['version']     = version
        color_law_dict['min_lambda']  = min_lambda
        color_law_dict['max_lambda']  = max_lambda
        #print(f" xxx color_law_dict = {color_law_dict} ")

        return color_law_dict
        # end get_color_law

    def merge_cleanup_final(self):

        script_dir       = self.config_prep['script_dir']

        # start with script_dir; separate tar files for CPU* and TRAINOPT*
        prefix = (f"CPU")
        util.compress_files(+1, script_dir, f"{prefix}*", prefix, "" )

        prefix = (f"{TRAINOPT_STRING}")
        util.compress_files(+1, script_dir, f"{prefix}*", prefix, "" )

        # - - - - - -
        # maybe later, CALIB -> CALIB.tar and TRAIN -> TRAIN.tar

        # end merge_cleanup_final

    def get_misc_merge_info(self):
        # return misc info to append in MERGE.LOG file
        misc_info = []
        return misc_info

    def get_merge_COLNUM_CPU(self):
        return COLNUM_TRAIN_MERGE_CPU

# === END: ===
