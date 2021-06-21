# Created Oct 31 2020
# 
# Nov 28 2020: 
#   + check for multiple filter-trans files; e.g., SNLS radial dependence.
#
# Jun 20 2021:
#   new keys 'SURVEY_LIST_SAMEMAGSYS and 'SURVEY_LIST_SAMEFILTER'
#   to record extra surveys in SALT2.INFO file.
#

import  os, sys, shutil, yaml, glob
import  logging, coloredlogs
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *
from    submit_prog_base import Program

TRAINOPT_STRING  = "TRAINOPT"
SALTPATH_STRING  = "SALTPATH"  # env name of calib path for snpca
SALTPATH_FILES   = [ 'fitmodel.card',  'Instruments',  'MagSys' ]
FILTERWHEEL_FILE = "FilterWheel"

# define subdirs under SUBDIR_CALIB ($SALTPATH)
SUBDIR_MAGSYS   = "MagSys"
SUBDIR_INSTR    = "Instruments"

# define keys for CONFIG inputs
KEY_PATH_INPUT_TRAIN = "PATH_INPUT_TRAIN"  # required (data,conf files)
KEY_PATH_INPUT_CALIB = "PATH_INPUT_CALIB"  # required (a.k.a, SALTPATH)
KEY_MODEL_SUFFIX     = "MODEL_SUFFIX"      # optional: change MODEL suffix

KEY_MAGSHIFT       = "MAGSHIFT"
KEY_WAVESHIFT      = "WAVESHIFT"
KEY_SHIFTLIST_FILE = "SHIFTLIST_FILE"

# Define suffix for output model used by LC fitters.
# Default output dirs are SALT2.MODEL000, SALT2.MODEL001, ...
MODEL_SUFFIX_DEFAULT = "MODEL"        

# Define list of trained SALT2 model files to check for existence;
# their existence defines SUCCESS/FAIL of training job.
COLORLAW_FILE   = "salt2_color_correction.dat"
CHECK_FILE_LIST = [ "salt2_template_0.dat", "salt2_template_1.dat", 
                    COLORLAW_FILE ]

ICOL_INFO_TRAINOPT    = 0
ICOL_INFO_FILE        = 1
ICOL_INFO_KEY         = 2
ICOL_INFO_SURVEY      = 3
ICOL_INFO_INSTR       = 4
ICOL_INFO_BAND        = 5
ICOL_INFO_SHIFT       = 6
ICOL_INFO_COMMENT     = 7

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

KEYS_SURVEY_LIST_SAME = ['SURVEY_LIST_SAMEMAGSYS', 'SURVEY_LIST_SAMEFILTER']

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

        # read survey map to magSys and Instruments
        self.train_prep_survey_map()

        # scoop up TRAINOPT list from user CONFIG
        self.train_prep_trainopt_list()

        # foreach training, prepare output paths
        self.train_prep_paths()

        # check for duplicate shifts and flag erorrs/warnings
        self.train_prep_error_checks()

        # end submit_prepare_driver

    def train_prep_survey_map(self):

        # read survey.yaml map so that user-input calibration shifts
        # can be converted into a specific file & key-value to modify.
        # This map is NOT used for nominal training; it is used only
        # for systematic variations in the calibration.

        CONFIG           = self.config_yaml['CONFIG']
        PATH_INPUT_CALIB = CONFIG[KEY_PATH_INPUT_CALIB] # aka SALTPATH
        PATH_EXPAND      = os.path.expandvars(PATH_INPUT_CALIB)
        msgerr = []

        if not os.path.exists(PATH_EXPAND):
            msgerr.append(f"Cannot find path for")
            msgerr.append(f"  {PATH_INPUT_CALIB}")
            msgerr.append(f"Check arg of {KEY_PATH_INPUT_CALIB}.")
            self.log_assert(False,msgerr)

        survey_map_file  = (f"{PATH_EXPAND}/survey.yaml")
        survey_yaml      = util.extract_yaml(survey_map_file, None, None )
        
        self.config_prep['survey_yaml']      = survey_yaml
        self.config_prep['survey_map_file']  = survey_map_file

        # end train_prep_survey_map

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
                                              KEY_SHIFTLIST_FILE )

        n_trainopt          = trainopt_dict['n_jobopt']
        trainopt_arg_list   = trainopt_dict['jobopt_arg_list']
        trainopt_ARG_list   = trainopt_dict['jobopt_ARG_list']
        trainopt_num_list   = trainopt_dict['jobopt_num_list']  
        trainopt_label_list = trainopt_dict['jobopt_label_list']
        trainopt_shift_file = trainopt_dict['jobopt_file_list']
        use_shift_file      = trainopt_dict['use_arg_file']

        logging.info(f" Store {n_trainopt-1} TRAIN-SALT2 options " \
                     f"from {TRAINOPT_STRING} keys")

        self.config_prep['n_trainopt']          = n_trainopt
        self.config_prep['trainopt_arg_list']   = trainopt_arg_list
        self.config_prep['trainopt_ARG_list']   = trainopt_ARG_list
        self.config_prep['trainopt_num_list']   = trainopt_num_list
        self.config_prep['trainopt_label_list'] = trainopt_label_list
        self.config_prep['trainopt_shift_file'] = trainopt_shift_file
        self.config_prep['use_shift_file']      = use_shift_file

        # end train_prep_trainopt_list
    
    def get_path_trainopt(self,which,trainopt):
        output_dir    = self.config_prep['output_dir']
        path          = ""
        if which == SUBDIR_CALIB_TRAIN :
            path = (f"{output_dir}/{SUBDIR_CALIB_TRAIN}/{trainopt}")
        elif which == SUBDIR_OUTPUT_TRAIN :
            path = (f"{output_dir}/{SUBDIR_OUTPUT_TRAIN}/{trainopt}")
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
        update_calib_info  = []

        topdir_calib    = (f"{output_dir}/{SUBDIR_CALIB_TRAIN}")
        topdir_train    = (f"{output_dir}/{SUBDIR_OUTPUT_TRAIN}")
        os.mkdir(topdir_calib)
        os.mkdir(topdir_train)

        for trainopt,arg in zip(trainopt_num_list,trainopt_arg_list) :
            outdir_calib = self.get_path_trainopt(SUBDIR_CALIB_TRAIN,trainopt)
            outdir_train = self.get_path_trainopt(SUBDIR_OUTPUT_TRAIN,trainopt)
            outdir_model = self.get_path_trainopt("MODEL",trainopt)
            outdir_calib_list.append(outdir_calib)
            outdir_train_list.append(outdir_train)
            outdir_model_list.append(outdir_model)
            os.mkdir(outdir_calib)
            os.mkdir(outdir_train)
            os.mkdir(outdir_model)
            info = self.train_prep_SALTPATH(outdir_calib,trainopt,arg)
            if info is not None : update_calib_info += info

        self.config_prep['outdir_calib_list']  =  outdir_calib_list
        self.config_prep['outdir_train_list']  =  outdir_train_list
        self.config_prep['outdir_model_list']  =  outdir_model_list
        self.config_prep['update_calib_info' ] =  update_calib_info

        logging.info("")

        #sys.exit(f"\n xxx DEBUG EXIT xxxx\n")

        # end train_prep_paths

    def train_prep_SALTPATH(self, outdir_calib, trainopt, arg):

        # prepare calibration file(s) pointed to by ENV $SALTPATH.
        # First copy base calibration files, then examine input arg
        # (passed from TRAINOPT input) and modify calibratoin files
        # (MagSys and Instrument filters) accordingly.
        #
        # Function returns calib_updates, a diagnostic info array for 
        # the SUBMIT_INFO file to enable humans to trace changes.
        # Note that calib_updates is strictly for diagnostics and
        # not used here internally.

        CONFIG           = self.config_yaml['CONFIG']
        PATH_INPUT_CALIB = CONFIG[KEY_PATH_INPUT_CALIB] # aka SALTPATH
        msgerr           = []
        calib_updates    = []
        replace_path_calib = False

        logging.info(f"   Prepare {SUBDIR_CALIB_TRAIN}/{trainopt}")

        # check for alternative calib directory (Nov 25 2020)
        arg_list = arg.split()
        if len(arg_list) > 1 :
            if arg_list[0] == KEY_PATH_INPUT_CALIB :
                PATH_INPUT_CALIB = os.path.expandvars(arg_list[1])
                logging.info(f"\t Replace {KEY_PATH_INPUT_CALIB}")
                replace_path_calib = True
                
        # - - - - - 
        for item in SALTPATH_FILES :
            cmd_rsync = f"rsync -r {PATH_INPUT_CALIB}/{item} {outdir_calib}"
            os.system(cmd_rsync)

        if arg == '' or replace_path_calib : return []

        # search each item in arg_list to allow for multiple sets
        # of MAGSHIFT or WAVESHIFT in a given TRAINOPT
        iarg = 0
        for item in arg_list :
            if item == KEY_MAGSHIFT or item == KEY_WAVESHIFT :
                arg_sublist = arg_list[iarg:iarg+4]
                key         = item
                survey      = arg_sublist[1]
                band_list   = arg_sublist[2].split(",")
                shift_list  = arg_sublist[3].split(",")

                if key == KEY_MAGSHIFT :
                    magsys_file,Instr_list = \
                        self.get_magsys_file(outdir_calib, survey)

                    magsys_dict = {
                        'magsys_file' :  magsys_file,
                        'survey'      :  survey,
                        'Instr_list'  :  Instr_list,
                        'band_list'   :  band_list,
                        'shift_list'  :  shift_list
                    }
                    info = self.update_magsys_file(magsys_dict)
                    calib_updates += info 
                    
                elif key == KEY_WAVESHIFT :
                    n_update    = 0
                    for band,shift in zip(band_list,shift_list):
                        shift = float(shift)
                        filter_file_list,Instr_list = \
                            self.get_filter_file(outdir_calib, survey, band)

                        for filt_file, Instr in zip(filter_file_list,Instr_list):
                            filter_dict = {
                                'filter_file' :  filt_file,
                                'survey'      :  survey,
                                'Instr'       :  Instr,
                                'band'        :  band,
                                'shift'       :  shift
                            }
                            info = self.update_filter_file(filter_dict)

                            # for SNLS, only append once to represent
                            # the radial dependence
                            n_update += 1
                            add_info = True
                            if survey == 'SNLS' and n_update>1 : add_info=False

                            if add_info : calib_updates += info
                else:
                    pass # probably should abort 

            # - - - -
            iarg += 1

        # ---------------------
        # prepend trainopt to each calib_updates item
        calib_updates_final = []
        for item in calib_updates:
            calib_updates_final.append( [trainopt] + item)

        return calib_updates_final
        # end train_prep_SALTPATH

    def get_magsys_file(self,outdir_calib,survey):

        # for input outdir_calib and survey, return full name of
        # magsys_file, and also return Instrument list.

        survey_yaml = self.config_prep['survey_yaml']
        magsys_file = "null"
        Instr       = "Null"
        msgerr      = []

        if survey in survey_yaml :
            Instr_list       = survey_yaml[survey]['Instrument']
            magsys_file_base = survey_yaml[survey]['MagSys_filename']
            magsys_dir  = (f"{outdir_calib}/{SUBDIR_MAGSYS}")
            magsys_file = (f"{magsys_dir}/{magsys_file_base}")
            if not os.path.exists(magsys_file):
                msgerr.append(f"Cannot find MagSys file:")
                msgerr.append(f"  {magsys_file} ")
                msgerr.append(f"for survey = {survey} and " \
                              f"Instrument = {Instr_list} ")
                self.log_assert(False,msgerr)
        else:
            input_file  = self.config_yaml['args'].input_file 
            msgerr.append(f"Invalid TRAINOPT survey: {survey}")
            msgerr.append(f"Check TRAINOPT in {input_file}")
            self.log_assert(False,msgerr)
    
        return magsys_file, Instr_list
        # end get_magsys_file

    def update_magsys_file(self, magsys_dict):

        # Edit/update magsys_file : search for row containing 
        # instrument in Intr_list and band in band_list;
        # shift mag value by shift_list.
        #
        # E.g., 
        #   Instr_list = [STANDARD, LANDOLT]
        #   band_list  = U,I
        #   shift_list = -0.01, 0.01
        #
        # Mag values are shifted by -0.01 for U band, and +0.01 for I-band.
        # These shifts are done for Instrument names STANDARD & LANDOLT.
        #
        # Function return human-readable "info" array destined for 
        # SUBMIT.INFO file.

        magsys_file = magsys_dict['magsys_file']
        survey      = magsys_dict['survey']
        Instr_list  = magsys_dict['Instr_list']
        band_list   = magsys_dict['band_list']
        shift_list  = magsys_dict['shift_list']

        info = []
        if not os.path.exists(magsys_file): return info

        magsys_base = os.path.basename(magsys_file)

        msgerr     = []
        nchange    = 0
        line_list_out = []
        line_list_inp = []
        line_list_out.append(f"# This file has been modified by\n")
        line_list_out.append(f"#    {sys.argv[0]}\n\n")

        # scoop up lines in file
        with open(magsys_file,"rt") as f:
            for line in f:
                line_list_inp.append(line)

        # look for line(s) to modify
        for line in line_list_inp:
            line_list_out.append(line) # default new line = old line

            # xxx mark if line[0] == '#' : continue
            if util.is_comment_line(line) : continue

            word_list = (line.rstrip("\n")).split()
            if len(word_list) < 3 : continue

            for band, shift in zip(band_list,shift_list):
                Instr = word_list[0]
                match = (Instr in Instr_list and word_list[1] == band)
                if not match : continue 
                strval = word_list[2] 
                try :
                    value = float(strval)
                except:
                    msgerr.append(f"Expected float after '{key_string}'" )
                    msgerr.append(f"but found {strval}")
                    msgerr.append(f"Check file: {file_name}")
                    self.log_assert(False,msgerr)

                new_value  = value + float(shift)
                new_line   = (f"{Instr} {band} {new_value:.6f}" \
                              f"   # {value} + {shift} \n")
                line_list_out[-1:] = new_line  # overwrite last line_out
                nchange += 1
                comment  = (f"{Instr}-{band} = {value} -> {new_value:.6f}")
                logging.info(f"\t Update {magsys_base}: {comment}" )
                info.append( [magsys_base, KEY_MAGSHIFT, 
                              survey, Instr, band, str(shift), comment] )

        # modify file with line_list_out
        if nchange > 0 :
            with open(magsys_file,"wt") as f:
                f.write("".join(line_list_out))

        return info

        # end update_magsys_file

    def get_filter_file(self,outdir_calib, survey, band ):
        
        # Return list if filter-transmission files corresponding to
        # inputs
        #   + outdir_calib  : full dir name within SALTPATH
        #   + survey        : name of survey (e..g, CfA1)
        #   + band          : band char
        #
        # Note that more than one filter-trans file can be returned;
        # hence a list. Instrument list is also returned.

        survey_yaml = self.config_prep['survey_yaml']
        FILTER_FILE = "NULL"
        filter_file = "null"
        Instr       = "Null"
        Instr_list          = []
        Instr_outlist       = []
        filter_file_outlist = []

        
        if survey in survey_yaml :
            subdir_list = survey_yaml[survey]['Instrument_subdir']
            Instr_list  = survey_yaml[survey]['Instrument']

            for subdir,Instr in zip(subdir_list, Instr_list):
                filter_dir  = (f"{outdir_calib}/{SUBDIR_INSTR}/{subdir}")
                FilterWheel = (f"{filter_dir}/{FILTERWHEEL_FILE}")
                flist       = self.parse_FilterWheel(FilterWheel,band)

                #print(f" xxx - - - - - - - -")
                #print(f" xxx FilterWheel = {FilterWheel}")
                #print(f" xxx flist = {flist} ")
                for filter_file in flist :
                    filter_file_outlist.append(f"{filter_dir}/{filter_file}")
                    Instr_outlist.append(Instr)
        else:
            msgerr      = []
            input_file  = self.config_yaml['args'].input_file 
            msgerr.append(f"Invalid TRAINOPT survey: {survey}")
            msgerr.append(f"Check TRAINOPT in {input_file}")
            self.log_assert(False,msgerr)
    
        return filter_file_outlist, Instr_outlist

        # end get_filter_file

    def parse_FilterWheel(self,FilterWheel,band):

        # parse FilterWheel file and for band return name(s) of
        # filter transmission file. FilterWheel file format is
        #   band filter_name fileName  
        # Note that filter_name is ignore here.
        # Note that SNLS has radial dependence, hence a list of
        # filter file names can be returned.
        #
        # Beware that FilterWheel can have 2 or 3 columns,
        # so filter_file is last column.

        filter_file_list = [ ]
        with open(FilterWheel,"rt") as f:
            for line in f:
                word_list = line.split()
                if word_list[0] == band:
                    filter_file_list.append(word_list[-1])
    
        return filter_file_list
        # end parse_FilterWheel

    def update_filter_file(self, filter_dict):
        # Modify filter_file to have wavelength shift of wave_shift. 
        # Function returns human-readable "info" array destined for 
        # SUBMIT.INFO file.

        filter_file = filter_dict['filter_file']
        survey      = filter_dict['survey']
        Instr       = filter_dict['Instr']
        band        = filter_dict['band']
        shift       = filter_dict['shift']  # wave shift, A

        if not os.path.exists(filter_file): 
            msgerr = []
            msgerr.append(f"Cannot update non-existent filter_file:")
            msgerr.append(f"{filter_file}")
            self.log_assert(False,msgerr)

        filter_base = os.path.basename(filter_file)
        comment     = (f"{shift} A shift")
        logging.info(f"\t Update {filter_base} with {comment}.")

        info = [ [ filter_base, KEY_WAVESHIFT, 
                   survey, Instr, band, str(shift), comment ] ]

        msgerr     = []
        line_list_out = []
        line_list_inp = []
        line_list_out.append(f"# wave_shift = {shift} A has been " \
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

            if util.is_comment_line(line) : continue
            # xxxx mark if line[0] == '#' : continue

            word_list = (line.rstrip("\n")).split()
            if len(word_list) < 2 : continue
            wave_orig = float(word_list[0])
            wave_out  = wave_orig + shift
            new_line  = str(wave_out) + "  " + " ".join(word_list[1:])
            line_list_out[-1:] = (f" {new_line} \n")

        # - - - - - - - -
        with open(filter_file,"wt") as f:
            f.write("".join(line_list_out))
            
        return info
        # end update_filter_file

    def train_prep_error_checks(self):

        update_calib_info = self.config_prep['update_calib_info']
        use_shift_file    = self.config_prep['use_shift_file']
        n_update = len(update_calib_info)
        nerr = 0 ; nwarn = 0; msgerr=[] 
        txt_error = f"ERROR: DUPLICATE CALIB-SHIFT in same TRAINOPT:"
        txt_warn  = f"WARNING: DUPLICATE CALIB-SHIFT in different TRAINOPTs:"

        for i in range(0,n_update):
            for j in range(i,n_update):

                update_i = update_calib_info[i]
                update_j = update_calib_info[j]

                #survey_snana_i, band_snana_i = self.get_SNANA_INFO(update_i)
                #survey_snana_j, band_snana_j = self.get_SNANA_INFO(update_j)
                #print(f" xxx i={i},j={j} survey_snana = " \
                #      f"{survey_snana_i} {survey_snana_j}")

                if i == j : continue

                TRAINOPT_i = update_i[ICOL_INFO_TRAINOPT]
                TRAINOPT_j = update_j[ICOL_INFO_TRAINOPT]

                FILE_i = update_i[ICOL_INFO_FILE]
                FILE_j = update_j[ICOL_INFO_FILE]

                KEY_i = update_i[ICOL_INFO_KEY]  # MAGSHIFT or WAVESHIFT
                KEY_j = update_j[ICOL_INFO_KEY]

                INSTR_i = update_i[ICOL_INFO_INSTR]
                INSTR_j = update_j[ICOL_INFO_INSTR]

                BAND_i = update_i[ICOL_INFO_BAND]
                BAND_j = update_j[ICOL_INFO_BAND]

                match_file  = (FILE_i     == FILE_j)
                match_instr = (INSTR_i    == INSTR_j) 
                match_band  = (BAND_i     == BAND_j)         
                match_opt   = (TRAINOPT_i == TRAINOPT_j)
                
                if match_file and match_instr and match_band :
                    txt_update = (f"   {update_i[0:7]}\n   {update_j[0:7]}")
                    if match_opt: 
                        logging.info(f"{txt_error}\n" \
                                     f"{txt_update}" )
                        nerr += 1
                    elif not use_shift_file :
                        logging.info(f"{txt_warn}\n" \
                                     f"{txt_update}")
                        nwarn += 1

        if nerr > 0 :
            msgerr.append(f"{nerr} duplicate calib shifts for same TRAINOPT")
            msgerr.append(f"Check ERROR messages above.")
            self.log_assert(False,msgerr)

        ###sys.exit("\n xxx DEBUG EXIST xxx ")
        # end train_prep_error_checks

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

        # open CMD file for this icpu  
        # xxxx mark delete f = open(COMMAND_FILE, 'a')

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
        msgerr = []

        #outdir_saltpath_list = self.config_prep['outdir_saltpath_list'] 
        outdir_train_list    = self.config_prep['outdir_train_list']
        outdir_model_list    = self.config_prep['outdir_model_list']

        PATH_INPUT_TRAIN  = os.path.expandvars(CONFIG[KEY_PATH_INPUT_TRAIN])
        if not os.path.exists(PATH_INPUT_TRAIN):
            msgerr.append(f"Cannot find path for")
            msgerr.append(f"  {PATH_INPUT_TRAIN}")
            msgerr.append(f"Check arg of {KEY_PATH_INPUT_TRAIN}.")
            self.log_assert(False,msgerr)

        outdir_train      = outdir_train_list[itrain]
        outdir_model      = outdir_model_list[itrain]

        # config file passed to pythin training shell
        TRAINDIR_FILE     = (f"{script_dir}/{trainDir_file}")

        # Nov 18 2020: add colons for yaml compatibility
        with open(TRAINDIR_FILE, "wt") as f :
            f.write(f"initDir:      {PATH_INPUT_TRAIN} \n")
            f.write(f"trainingDir:  {outdir_train} \n")
            f.write(f"outputDir:    {outdir_model} \n")

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

        CONFIG       = self.config_yaml['CONFIG']
        n_trainopt   = self.config_prep['n_trainopt'] 
        num_list     = self.config_prep['trainopt_num_list']
        arg_list     = self.config_prep['trainopt_arg_list'] 
        ARG_list     = self.config_prep['trainopt_ARG_list'] 
        label_list   = self.config_prep['trainopt_label_list']
        model_suffix = self.config_prep['model_suffix']
        survey_map_file = self.config_prep['survey_map_file'] 

        f.write(f"# train_SALT2 info \n")
        f.write(f"JOBFILE_WILDCARD: {TRAINOPT_STRING}* \n")
        f.write(f"MODEL_SUFFIX: {model_suffix}   " \
                f"# -> create SALT2.{model_suffix}nnn/ \n")

        f.write(f"SURVEY_MAP_FILE:  {survey_map_file} \n")

        f.write(f"\n")

        for key in KEYS_SURVEY_LIST_SAME:
            if key in CONFIG :
                f.write(f"{key}:  {CONFIG[key]} \n")
            else:
                f.write(f"{key}:  [ ] \n")

        f.write(f"\n")
        f.write(f"TRAINOPT_OUT_LIST:  " \
                f"# 'TRAINOPTNUM'  'user_label'  'user_args'\n")
        # use original ARG_list instead of arg_list; the latter may
        # include contents of shiftlist_file.
        for num, arg, label in zip(num_list, ARG_list, label_list):
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")
        
        # write which calib files were modified for systematics 
        # (human readable)
        update_calib_info = self.config_prep['update_calib_info' ]
        f.write(f"CALIB_UPDATES: # TRAINOPT    file          comment\n")
        for item_full in update_calib_info:
            item = [ item_full[ICOL_INFO_TRAINOPT], 
                     item_full[ICOL_INFO_FILE],
                     item_full[ICOL_INFO_COMMENT] ]
            f.write(f"  - {item}\n")
        f.write("\n")

        # write info for SNANA's SALT2.INFO file
        n_item = 0
        for item_full in update_calib_info:
            #print(f" xxx item_full = {item_full} ")
            survey_snana, band_snana = self.get_SNANA_INFO(item_full)
            if survey_snana is None : continue 

            trainopt  = item_full[ICOL_INFO_TRAINOPT]
            key_snana = item_full[ICOL_INFO_KEY]
            shift     = item_full[ICOL_INFO_SHIFT]
            arg_snana = f"{survey_snana} {band_snana} {shift}"
            item      = [ trainopt, key_snana, arg_snana ]

            n_item += 1
            if n_item == 1 : f.write(f"SNANA_SALT2_INFO: \n")
            f.write(f"  - {item} \n")
        f.write("\n")

        # end append_info_file

    def get_SNANA_INFO(self,info_list):

        survey_yaml  = self.config_prep['survey_yaml']
        survey_train = info_list[ICOL_INFO_SURVEY]
        instr_train  = info_list[ICOL_INFO_INSTR]
        band_train   = info_list[ICOL_INFO_BAND]
        Instr_list   = survey_yaml[survey_train]['Instrument']

        survey_snana = None 
        band_snana   = None         

        key_snana = 'SNANA_INSTR'
        if key_snana in survey_yaml[survey_train]:
            temp_list    = survey_yaml[survey_train][key_snana]
            j            = Instr_list.index(instr_train)
            survey_snana = temp_list[j]
            #survey_snana = ",".join(temp_list)

        key_snana = 'SNANA_BANDS'
        if key_snana in survey_yaml[survey_train]:
            band_list_train = survey_yaml[survey_train]['Bands']
            band_list_snana = survey_yaml[survey_train][key_snana]
            j               = band_list_train.index(band_train)
            band_snana      = band_list_snana[j]

        return survey_snana, band_snana
        # end get_SNANA_INFO

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
        train_dir    = self.get_path_trainopt(SUBDIR_OUTPUT_TRAIN,trainopt)

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
        train_dir  = self.get_path_trainopt(SUBDIR_OUTPUT_TRAIN,trainopt)
        calib_dir  = self.get_path_trainopt(SUBDIR_CALIB_TRAIN,trainopt)

        logging.info(f"    Compress output for {trainopt} :")

        # make tar file from CALIB/TRAINOPTnnn  (aka SALTPATH)
        logging.info(f"\t Compress {SUBDIR_CALIB_TRAIN}/{trainopt}")
        util.compress_subdir(+1,calib_dir)

        # Gzip contents of TRAINOPT, then  TRAINOPTnnn -> TRAINOPTnnn.tar.gz
        logging.info(f"\t Compress {SUBDIR_OUTPUT_TRAIN}/{trainopt}")
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
        info_file    = f"{model_dir}/{SALT2_INFO_FILE}"

        color_law_dict = self.get_color_law(model_dir)
        min_lambda   = color_law_dict['min_lambda']
        max_lambda   = color_law_dict['max_lambda']
        range_lambda = f"{min_lambda} {max_lambda}"
        npar         = color_law_dict['npar']
        par_list     = " ".join(color_law_dict['par_list'])

        logging.info(f"    Create {SALT2_INFO_FILE}")
        with open(info_file,"wt") as f:

            f.write(f"RESTLAMBDA_RANGE: {range_lambda}\n")
            f.write(f"COLORLAW_VERSION: {color_law_dict['version']} \n")
            f.write(f"COLORCOR_PARAMS:  {range_lambda} {npar} {par_list}\n")
            f.write(f"\n")

            f.write(f"MAG_OFFSET:    {SALT2_MAG_OFFSET} \n")
            f.write(f"SIGMA_INT:     {SALT2_SIGMA_INT}  \n")
            f.write(f"COLOR_OFFSET:  {SALT2_COLOR_OFFSET} \n")

            f.write(f"{SALT2_INFO_INCLUDE}\n")

            self.append_SALT2_INFO_TRAINOPT(f,trainopt)  

        #if trainopt == "TRAINOPT003" :
        #    sys.exit(f"\n xxx DEBUG STOP at merge_create_SALT2_INFO_file\n")

        # #end merge_create_SALT2_INFO_file
        
    def append_SALT2_INFO_TRAINOPT(self,f,trainopt):

        # Created Jun 2021
        # write MAGSHIFT and WAVESHIFT keys for SNANA to use.

        submit_info_yaml       = self.config_prep['submit_info_yaml']

        if 'SNANA_SALT2_INFO' not in submit_info_yaml : return

        SNANA_SALT2_INFO       = submit_info_yaml['SNANA_SALT2_INFO']
        SURVEY_LIST_SAMEMAGSYS = submit_info_yaml['SURVEY_LIST_SAMEMAGSYS']
        SURVEY_LIST_SAMEFILTER = submit_info_yaml['SURVEY_LIST_SAMEFILTER']

        f.write(f"# {trainopt} \n")
        for item in SNANA_SALT2_INFO :
            if trainopt == item[0] :
                key      = item[1]
                arg_list = item[2].split()
                survey   = arg_list[0]      # extract survey for check below
                arg      = " ".join(arg_list[1:])  # band value
                f.write(f"{key}: {survey:<10} {arg} \n")

                # write other surveys with same magsys
                if key == KEY_MAGSHIFT:
                    s_list  = SURVEY_LIST_SAMEMAGSYS
                else:
                    s_list  = SURVEY_LIST_SAMEFILTER

                if survey in s_list :
                    for s in filter(lambda s: s not in [survey], s_list):
                        comment = f"same {key} as {survey}"
                        f.write(f"{key}: {s:<10} {arg}    # {comment}\n")

        # end append_SALT2_INFO_TRAINOPT

    def get_color_law(self,model_dir):
        file_name = (f"{model_dir}/{COLORLAW_FILE}")
        color_law_dict = {}

        nline = 0; npar_tot = -1; npar_read=0; cpar_list = []

        with open(file_name,"r") as f:
            for line in f:
                nline += 1
                word_list = (line.rstrip("\n")).split()
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
