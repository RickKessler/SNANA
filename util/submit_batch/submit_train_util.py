# Created Jun 2023 by R.Kessler
# training utilities to be shared by SALT3 and BayeSN
# 

import  os, sys, shutil, yaml, glob, logging
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *

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

# define prefix for files with calib shifts.
PREFIX_CALIB_SHIFT   = "CALIB_SHIFT"  


KEYDICT_CALIBSHIFT_FILE = {
    METHOD_TRAIN_SALT3  : '--calibrationshiftfile' ,
    METHOD_TRAIN_BAYESN : '--calibrationshiftfile' 
}


#KEY_SALTshaker_CALIBSHIFT_FILE

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================


def train_prep_trainopt_list(METHOD,CONFIG,script_dir):

    # return dictionary of trainopt info

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


# =========== END: =======



