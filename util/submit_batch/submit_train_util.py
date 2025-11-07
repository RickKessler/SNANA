# Created Jun 2023 by R.Kessler
# training utilities to be shared by SALT3 and BayeSN
# 

import  os, sys, shutil, yaml, glob, logging
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *
import  numpy as np

METHOD_TRAIN_SALT3  = "SALT3"  # also known as SALTshaker
METHOD_TRAIN_BAYESN = "BAYESN"

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

# define calib key for SUBMIT.INFO file that passes to snana 
# LC fit via model output
KEY_SNANA_CALIB_INFO  = "SNANA_CALIB_INFO" 

# define prefix for files with calib shifts.
PREFIX_CALIB_SHIFT   = "CALIB_SHIFT"  

KEYDICT_CALIBSHIFT_FILE = {
    METHOD_TRAIN_SALT3  : '--calibrationshiftfile' ,
    METHOD_TRAIN_BAYESN : '--calibrationshiftfile' 
}


# name of file with list of input files
TRAIN_INPUT_FILENAMES = "INPUT_FILE.LIST" 

#KEY_SALTshaker_CALIBSHIFT_FILE

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================


def train_prep_trainopt_list(METHOD,CONFIG,script_dir):

    # Created Jun 15 2023 [Moved from submit_train_SALT3.py]
    # Read TRAINOPT keys, extract calib shifts, create input
    # calib-shift file(s), and return dictionary of trainopt info.
    #
    # Inputs:
    #  METHOD = METHOD_TRAIN_SALT3 or METHOD_TRAIN_BAYESN
    #  CONFIG = CONFIG block from submit-config-input file
    #  script_dir = location for writing calib-shift files (if needed)
    #


    trainopt_global = ""  # apply to each TRAINOPT
    trainopt_rows   = []  # TRAINOPT-specified commands

    # start with global settings
    key = TRAINOPT_GLOBAL_STRING
    if key in CONFIG:  
        trainopt_global = CONFIG[key]

    # next, TRAINOPT per job
    key      = TRAINOPT_STRING
    if key in CONFIG :  
        trainopt_rows = CONFIG[key]

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

    logging.info(f" Store {n_trainopt-1} TRAIN-{METHOD} options " \
                 f"from {TRAINOPT_STRING} keys")

    # for MAGSHIFT and./or WAVESHIFT, create calibration file
    # in script_dir. The returned calib_shift_file has each
    # calib arg unpacked so that each row has only 1 shift.
    # arg_replace = arg, but calibration shifts are replaced
    # with command to read calib-shift file.
    arg_replace_list = []
    calib_shift_list = []
    for num,arg in zip(trainopt_num_list,trainopt_arg_list):
        arg_replace, calib_shift = \
                make_calib_shift_file(num,arg,METHOD,script_dir)
        arg_replace_list.append(arg_replace)
        calib_shift_list.append(calib_shift)  # list of lists

    config_prep_local = {
        'n_trainopt'          : n_trainopt ,
        'trainopt_arg_list'   : arg_replace_list,
        'trainopt_ARG_list'   : trainopt_ARG_list,
        'trainopt_num_list'   : trainopt_num_list,
        'trainopt_label_list' : trainopt_label_list,
        'trainopt_shift_file' : trainopt_shift_file,
        'trainopt_global'     : trainopt_global,
        'use_arg_file'        : use_arg_file,
        'calib_shift_list'    : calib_shift_list
    }

    # xxx self.config_prep.update(config_prep_local)

    print('')

    return config_prep_local

    # end train_prep_trainopt_list

def make_calib_shift_file(num, arg, METHOD, script_dir):
    
    # if arg contains MAGSHIFT or WAVESHIFT key, write them out
    # in a calib-shift file.
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
    # Make sure that arg_replace retains non-calib options
    # to allow mixing calib and non-calib arguments.
    
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

    # tack on key & calibshift file
    KEY = KEYDICT_CALIBSHIFT_FILE[METHOD]
    arg_replace += f"{KEY} {calib_shift_file} "

    return arg_replace, calib_shift_list

    # end make_calshift_file


def prep_model_paths(METHOD, trainopt_num_list, output_dir ):

    # Created Jun 16 2023
    # for each TRAINOPT, create path for SALTPATH and for training output
    # Inputs:
    #   METHOD ;  TRAIN_SALT3 or TRAIN_BAYESN or ....
    #   trainopt_num_list: e.g., 000, 001, 002 ...
    #   outpout_dir: create [METHOD].[SUFFIX][nnn] directories here
    #     e.g., SALT3.MODEL001, SALT3.MODEL002, etc ...
    #
    #  Return list of created directories 

    model_suffix      = MODEL_SUFFIX_DEFAULT  
    outdir_model_list  = []
    n_trainopt = len(trainopt_num_list)

    print(f" Create {n_trainopt} output model directories:")
    for trainopt_num in trainopt_num_list:
        nnn           = trainopt_num[-3:]
        sdir_model    = f"{METHOD}.{model_suffix}{nnn}"
        outdir_model  = f"{output_dir}/{sdir_model}"
        print(f"\t Create {sdir_model}")
        outdir_model_list.append(outdir_model)
        os.mkdir(outdir_model)

    sys.stdout.flush()
    return outdir_model_list

    # end prep_model_paths


def append_info_file(f, append_info_dict):

    # append training info to SUBMIT.INFO file with pointer f.
    
    METHOD             = append_info_dict['METHOD']  # SALT3 or BAYESN ...
    CONFIG             = append_info_dict['CONFIG']
    n_trainopt         = append_info_dict['n_trainopt'] 
    trainopt_num_list  = append_info_dict['num_list']
    arg_list           = append_info_dict['arg_list'] 
    ARG_list           = append_info_dict['ARG_list'] 
    label_list         = append_info_dict['label_list']
    calib_shift_list   = append_info_dict['calib_shift_list']
    outdir_model_list  = append_info_dict['outdir_model_list']
    output_dir         = append_info_dict['output_dir']

    # - - - - - - 
    f.write(f"# train_{METHOD} info \n")
    f.write(f"JOBFILE_WILDCARD: {TRAINOPT_STRING}* \n")

    f.write(f"\n")
    f.write(f"TRAINOPT_OUT_LIST:  " \
            f"# 'TRAINOPTNUM'  'user_label'  'user_args'\n")
    # use original ARG_list instead of arg_list; the latter may
    # include contents of shiftlist_file.
    for trainopt, arg, label in zip(trainopt_num_list, ARG_list, label_list):
        row   = [ trainopt, label, arg ]
        f.write(f"  - {row} \n")
    f.write("\n")

    f.write("MODELDIR_LIST:\n")
    for model_dir in outdir_model_list :
        model_dir_base = os.path.basename(model_dir)
        f.write(f"  - {model_dir_base}\n")
    f.write("\n")
    
    for key in KEYS_SURVEY_LIST_SAME:
        if key in CONFIG :
            f.write(f"{key}:  {CONFIG[key]} \n")
        else:
            f.write(f"{key}:  [ ] \n")

    # - - - - - 
    # write keys for SALT3.INFO to be read by SNANA code 
    # each row is
    #  [ 'TRAINOPTnnn', KEY, SURVEY, SHIFT_VAL ]
    f.write(f"{KEY_SNANA_CALIB_INFO}: \n")
    for trainopt, item_list in zip(trainopt_num_list, calib_shift_list) :
        for item in item_list:
            row = [ trainopt ] + item.split()
            f.write(f"  - {row} \n")
    f.write("\n")

    # - - - - - - - - - 
    # Oct 28 2025: write csv table of shifts to enable diagnostic plots
    rownum = 0
    calib_table =  f"{output_dir}/CALIB_SHIFTS_TABLE.DAT"

    magshift_dict = {}
    waveshift_dict = {}
    table_lines = []  # store table lines to allow pre-pending stats at stop of file
    table_lines.append(f"ROW  MODELNUM  SURVEY  BAND  VARNAME_SHIFT  SHIFT")

    for trainopt, item_list in zip(trainopt_num_list, calib_shift_list) :
        trainopt_num = int(trainopt.split('TRAINOPT')[1])
        for item in item_list:
            wdlist = item.split()
            var_shift = wdlist[0]  # MAGSHFIT or WAVESHIFT
            survey    = wdlist[1]
            band      = wdlist[2]
            shift     = wdlist[3]  # shift value (mag or Angstroms)
            rownum += 1
            
            line = f"{rownum:4d}  {trainopt_num:2d}  {survey:<12}  {band:<8}  " \
                   f"{var_shift:<10} {shift} "

            table_lines.append(line)

            # store list of shift values by survey and var to get stats
            if 'MAG' in var_shift:
                if survey not in magshift_dict:  magshift_dict[survey]  = []
                magshift_dict[survey].append(float(shift))
            if 'WAVE' in var_shift:
                if survey not in waveshift_dict: waveshift_dict[survey] = []
                waveshift_dict[survey].append(float(shift))

    # - - - - - -
    t = open(calib_table,"wt")
    t.write(f"# ======================================================================== \n")
    print_calib_shift_stats(t, "MAGSHIFT",  magshift_dict)
    print_calib_shift_stats(t, "WAVESHIFT", waveshift_dict)
    t.write(f"# ======================================================================== \n")
    t.write('\n')
    for line in table_lines:
        t.write(f"{line}\n")
    t.close()

    return

    # end append_info_file
    
def print_calib_shift_stats(f, shift_type, shift_dict):
    # Created Oct 29 2025
    # Inputs
    #  f = file pointer to write to
    #  shift_type = MAGSHIFT or WAVESHIFT
    #  shift_dict[survey] = [ list of shift values ]
    #
    # Print mean, stddev, rms to file pointer f

    for survey, shift_list in shift_dict.items():
        shift_list_np = np.array(shift_list)
        key     = f"{survey}_{shift_type}:"
        nval    = len(shift_list)
        mean    = np.mean(shift_list_np)
        std_dev = np.std(shift_list_np)
        rms     = np.sqrt(np.mean(shift_list_np**2))
        f.write(f"# {key:<20}  n={nval:3d}  mean={mean:7.4f}  std_dev={std_dev:7.4f}  rms={rms:7.4f}\n")

    return

# =========== END: =======



