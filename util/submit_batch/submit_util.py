# ==============================================
# Created July 2020 by R.Kessler & S. Hinton
#
# generic utilites for submit  script
#
# Mar 7 2022: use -x instead of -f to check for program ...
#             hope to avoid running .exe during make.
#
# ==============================================

import os, sys, yaml, shutil, glob, math, ntpath
import logging, coloredlogs, subprocess
import pandas as pd
from   submit_params import *

# =================================================

def combine_csv_files(wildcard, combined_file, remove_flag=False):

    # Created Nov 18 2021
    # combine csv files from wildcard -> 
    # output combined_file includes contents from all wildcard files.
    # Nov 30 2021: return if there are no files with wildcard.

    csv_file_list = sorted(glob.glob(wildcard))

    n_csv = len(csv_file_list)
    logging.info(f"  Combine {n_csv} csv truth tables for\n\t {wildcard}")
    if n_csv ==0 : return

    # read table contents as strings to avoid modifying float format 
    # in the combined csv.            
    combined_csv = pd.concat([pd.read_csv(f,dtype=str) \
                              for f in csv_file_list ] )

    # write it all out in one combined file                                 
    combined_csv.to_csv(combined_file, index=False)

    # remove original csv files 
    if remove_flag :
        cmd_rm = f"rm {wildcard}"
        os.system(cmd_rm)

    # end combine_csv_files

def get_wfit_values(wfit_yaml):

    # Created Aug 9 2021
    # parse yaml for wfit values, allowing for legacy and 
    # refactored (Aug 2021) wfit. 
    # Also check for wsig_marg vs. wsig_lo/wsig_hi
    # Sep 28 2021: check for wa and its uncertainty
    # Oct 23 2021: check for Rho
    # Feb 22 2022: check for NWARNINGS

    key_list = [ 'w', 'w0' ]
    for key in key_list:
        if  key in wfit_yaml:
            w  = wfit_yaml[key]  

    key_list = [ 'w_sig', 'wsig_marg', 'wsig_lo', 
                 'w0sig_marg', 'w0sig_lo' ]
    w_sig    = -9.0
    for key in key_list:
        if key in wfit_yaml: 
            w_sig = wfit_yaml[key]
            if key == 'wsig_lo' :
                w_sig_lo = wfit_yaml['wsig_lo'] 
                w_sig_hi = wfit_yaml['wsig_hi'] 
                w_sig    = 0.5*(w_sig_lo + w_sig_hi)
            if key == 'w0sig_lo' :
                w_sig_lo = wfit_yaml['w0sig_lo'] 
                w_sig_hi = wfit_yaml['w0sig_hi'] 
                w_sig    = 0.5*(w_sig_lo + w_sig_hi)


    # - - - repeat for optoinal wa
    key_list = [ 'wa' ]
    wa       = 0.0
    for key in key_list:
        if  key in wfit_yaml:
            wa  = wfit_yaml[key]  

    key_list = [ 'wasig_marg', 'wasig_lo' ]
    wa_sig   = 0.0
    for key in key_list:
        if key in wfit_yaml: 
            wa_sig = wfit_yaml[key]
            if key == 'wasig_lo' :
                wa_sig_lo = wfit_yaml['wasig_lo'] 
                wa_sig_hi = wfit_yaml['wasig_hi'] 
                wa_sig    = 0.5*(wa_sig_lo + wa_sig_hi)

    # - - - OM - - - -
    key_list = [ 'omm', 'OM' ]
    OM = -9.0
    for key in key_list:
        if  key in wfit_yaml:
            omm  = wfit_yaml[key]  

    key_list = [ 'omm_sig', 'OMsig', 'OMsig_marg' ]
    omm_sig = -9.0
    for key in key_list:
        if  key in wfit_yaml:
            omm_sig  = wfit_yaml[key]  

    # - - - repeat for FoM (for w0wa fit)
    key_list = [ 'FoM' ]
    FoM       = 0.0
    for key in key_list:
        if  key in wfit_yaml:
            FoM  = wfit_yaml[key]  

    # - - - repeat for reduced cov Rho(w0wa) - Oct 23 2021
    key_list = [ 'Rho' ]
    Rho       = 0.0
    for key in key_list:
        if  key in wfit_yaml:
            Rho  = wfit_yaml[key]  

    # - - - misc - - - - 
    chi2    = wfit_yaml['chi2'] 
    sigint  = wfit_yaml['sigint']

    key_list = [ 'wrand', 'wran', 'w0ran' ]
    w_ran = -9.0 
    for key in key_list:
        if key in wfit_yaml:
            w_ran   = wfit_yaml[key]

    key_list = [ 'warand', 'waran' ]
    wa_ran   = 0
    for key in key_list:
        if key in wfit_yaml:
            wa_ran   = wfit_yaml[key]

    key_list = [ 'ommrand', 'ommran', 'OMran' ]
    omm_ran = -9.0
    for key in key_list:
        if key in wfit_yaml:
            omm_ran   = wfit_yaml[key]

    key_list = [ 'BLIND', 'blind' ]
    blind    = 0
    for key in key_list:
        if key in wfit_yaml:
            blind = wfit_yaml[key]

    key_list = [ 'NWARNINGS' ]
    nwarn    = 0
    for key in key_list:
        if key in wfit_yaml:
            nwarn = wfit_yaml[key]

    wfit_values_dict = {
        'w'        : w ,
        'w_sig'    : w_sig ,
        'omm'      : omm ,
        'omm_sig'  : omm_sig ,
        'chi2'     : chi2 ,
        'sigint'   : sigint ,
        'w_ran'    : w_ran,
        'omm_ran'  : omm_ran,
        'blind'    : blind ,
        'nwarn'    : nwarn ,
        # optional below
        'wa'       : wa,
        'wa_sig'   : wa_sig,
        'wa_ran'   : wa_ran,
        'FoM'      : FoM,
        'Rho'      : Rho
    }

    return wfit_values_dict

    # end get_wfit_values

def prep_jobopt_list(config_rows, string_jobopt, key_arg_file):

    # Created Jan 23 2021
    # Generic utility to strip args from config_rows and return
    # dictionary of arg_list, num_list, label_list.
    # The args are from FITOPT key for FIT, MUOPT key for BBC,
    # TRAINOPT key for training ... Each jobopt can include
    # a label surrounding by slashes; e.g., /LABEL_THIS/
    #
    # Inputs:
    #   config_rows = CONFIG[string_jobopt] passed from user input file
    #   string_jobopt = 'FITOPT' or'MUOPT'or 'TRAINOPT' ...
    #   key_arg_file = optional key for file name containing args

    jobopt_dict = {} # init output dictionary
    n_jobopt          = 1
    jobopt_ARG_list   = [ '' ] # original user arg
    jobopt_arg_list   = [ '' ] # expanded args used by script
    jobopt_num_list   = [ f"{string_jobopt}000" ] 
    jobopt_label_list = [ None ]
    jobopt_file_list  = [ None ]   
    use_arg_file      = False

    if 'WFIT' in string_jobopt :  # no default 000 job for wfit 
        n_jobopt = 0
        jobopt_ARG_list = [] 
        jobopt_arg_list = []
        jobopt_num_list = []
        jobopt_label_list = []
        jobopt_file_list  = []

    for jobopt_raw in config_rows :    # might include label
        num = f"{string_jobopt}{n_jobopt:03d}"

        # separate  label and ARG in
        #    jobopt_string: /label/ ARG
        label, ARG = separate_label_from_arg(jobopt_raw)

        # if arg points to an arg_file, read and return args from file
        arg, arg_file  = read_arg_file(ARG,key_arg_file)
        if arg_file is not None:  use_arg_file = True

        # uopdate jobopt lists
        jobopt_arg_list.append(arg)
        jobopt_ARG_list.append(ARG)
        jobopt_num_list.append(num)
        jobopt_label_list.append(label)
        jobopt_file_list.append(arg_file)
        n_jobopt += 1
     
    # - - - -
    jobopt_dict = {
        'jobopt_arg_list'   : jobopt_arg_list ,
        'jobopt_ARG_list'   : jobopt_ARG_list ,
        'jobopt_num_list'   : jobopt_num_list ,
        'jobopt_label_list' : jobopt_label_list ,
        'jobopt_file_list'  : jobopt_file_list ,
        'n_jobopt'          : n_jobopt,
        'use_arg_file'      : use_arg_file
    }

    return jobopt_dict

    # end prep_jobopt_list

def read_arg_file(ARG, KEY_ARG_FILE):

    # if ARG starts with KEY_ARG_FILE, the return arg is equal
    # to the contents of arg_file; otherwise return arg = ARG.
    # Motivation is that user can build long list of arguments
    # (e.g., random calib variations for training)
    # and store each set of args in a separate file.

    arg      = ARG   # init return arg to input ARG
    arg_file = None  # init return arg_file

    if KEY_ARG_FILE is None :
        return arg, arg_file
    
    # - - - - -
    word_list = ""
    arg_list = ARG.split()

    if arg_list[0] == KEY_ARG_FILE :
        arg_file = arg_list[1]
        with open(arg_file,"rt") as f:
            for line in f:
                if is_comment_line(line) : continue
                word_list += line.replace("\n"," ")
            arg = word_list
    # - - - -
    return arg, arg_file
    # end read_arg_file


def protect_wildcard(arg):
    # Created May 2022 by R.Kessler
    # if arg = abc*, return 'abc*'
    # if arg = 'abc*', return 'abc*' [no change]

    arg_protect = arg
    if isinstance(arg,str):
        has_quote = '\'' in arg  or  '"' in arg
        if '*' in arg and not has_quote:
            arg_protect = f"'{arg}'"

    return arg_protect
    # end protect_wildcard

def protect_parentheses(arg):
    # Created Dec 10 2020
    # if arg = abc(option), returns abc\(option\).
    # If arg = abc, returns abc (no change)
    # This protection is needed to read GENOPT, FITOPT, MUOPT  args .
    # M. Vincenzi Febr 2022: only protects strings to avoid error on int/float
    if isinstance(arg,str):
        arg_protect = arg.replace('(','\(').replace(')','\)')
    else: 
        arg_protect = arg
    return arg_protect
    # end protect_parentheses

def is_comment_line(line):
    if len(line) == 0  : return True
    if line[0] == '#'  : return True
    if line[0] == '@'  : return True
    if line[0] == '%'  : return True
    if line[0] == '!'  : return True
    if line[0] == '\n' : return True

    return False

def fix_partial_path(file_list):

    # if any file in file_list has a partial path, tack on CWD
    # e.g., file = sim/abc.input, then out_file_list has
    # [CWD]/sim/abc.input
    
    out_file_list = [] 
    for f0 in file_list :
        f1 = os.path.expandvars(f0)
        if '/' in f1 and f1[0] != '/' :  f1 = f"{CWD}/{f0}"
        out_file_list.append(f1)

    return out_file_list
    # end fix_partial_path
        
def separate_label_from_arg(input_arg_string):

    # If input_arg_string = /LABEL/ x1min=-2.0 nzbin=20
    # then return
    #   label = LABEL
    #   arg_string = x1min=-2.0 nzbin=20
    #
    #  If there is no label, return label = None

    # init output for no label
    label = None ;   arg_string = input_arg_string
    
    if len(input_arg_string) > 0 :
        word_list = input_arg_string.split()
        has_label = word_list[0][0] == '/'
        if has_label :
            label        = word_list[0].strip('/') 
            arg_string   = " ".join(word_list[1:])
            arg_string   = protect_parentheses(arg_string)

    return label, arg_string
    # end separate_label_from_arg

def standardise_path(path,cwd):
    if "$" in path:
        path = os.path.expandvars(path)
        if "$" in path:
            msgerr = []
            msgerr.append(f"Unable to resolve ENV in {path}")
            msgerr.append(f"Check how your ENV is set.")
            log_assert(False,msgerr)

    if not path.startswith("/"):
        path = os.path.join(cwd, path)

    return path

def print_debug_line(line):
    print(f"\n DEBUG_DUMP: \n {line} \n DEBUG_DUMP_END: \n")

def find_and_remove(find_arg):

    # Called as part of purge option:
    # + search for 'find_arg' files using linux find command
    # + print each file, along with size (MB)
    # + ask user to remove (y/[n])
    #
    # Make sure that input find_arg includes appropriate wildcards;
    # e.g, SPLIT_JOBS_LCFIT* to include tar files.
    
    remove_list = []
    remove_size = []
    cmd_find     = f"find . -name {find_arg}" + " -exec du -mc {} +"
    find_list    = subprocess.check_output(cmd_find, shell=True)
    find_list    = (find_list.rstrip()).decode('utf-8')
    find_list    = find_list.split()
    remove_size += find_list[0::2]  # every other element is size (MB)
    remove_list += find_list[1::2]  # every other elment if file or dir

    # - - - - - - 
    print(f"# ======================================================== ")

    # print summary of files/directories to remove (but don't remove, yet)
    n_file = len(remove_list) - 1  # leave out the 'total' line
    for i in range(0,n_file+1):    # include 'total' line for summary
        remove_file = remove_list[i]
        size        = remove_size[i]
        print(f"   Found {remove_file}  ({size} MB)")

    if n_file < 0 :
        print(f"  No files found to remove for {find_arg} ")
        return

    response = input(f"\n   Remove {n_file} {find_arg} files above y/[n] ? ")
    if response == 'y' :
        print(f"\t Removing {n_file} files ... ")
        for i in range(0,n_file):
            remove_file = remove_list[i]
            cmd_rm = f"rm -r {remove_file}"
            os.system(cmd_rm)
    else:
        print(f"\t Do not remove {find_arg} files")

    # end find_and_remove

def get_stat_dict(value_list,error_list):
    # For input list of values_list, return dictionary of
    # 'AVG_VAL', 'ERR_AVG',  "AVG_ERR',  'RMS', 'ERR_RMS'

    n_val    = len(value_list)
    
    if n_val > 0 :
        sumval   = 0.0 ; sqsumval = 0.0;  sumerr=0.0
        for val,err in zip(value_list,error_list) :
            sumval   += val
            sqsumval += val*val
            sumerr   += err
            
            AVG_VAL = sumval / n_val
            RMS     = math.sqrt( sqsumval/n_val - AVG_VAL*AVG_VAL )
            AVG_ERR = sumerr / n_val

            ERR_AVG = RMS/math.sqrt(n_val)
            ERR_RMS = ERR_AVG / 1.414   # sigma/sqrt(2*n)
            
    else:
        AVG_VAL = 0.0; AVG_ERR=0.0; ERR_AVG = 0.0; RMS=0.0; ERR_RMS=0.0

    stat_dict = { 'AVG_VAL': AVG_VAL,  'AVG_ERR': AVG_ERR,
                  'ERR_AVG': ERR_AVG,  'RMS': RMS, 'ERR_RMS': ERR_RMS }

    return stat_dict

    # end get_stat_dict

def merge_table_reset(merge_file, table_name, colnum_state, colnum_zero_list):

    # read {table_name} from {merge_file} and reset table as follows:
    #   STATE -> WAIT
    #   value -> zero for columns in colnum_zero_list
    # Then re-write MERGE.LOG

    INFO_CONTENTS,comment_lines = read_merge_file(merge_file)
    row_list       = INFO_CONTENTS[table_name]
    row_reset_list = []
    
    for row in row_list :
        row_reset = row
        row[colnum_state] = SUBMIT_STATE_WAIT
        for colnum in colnum_zero_list :
            row[colnum] = 0
        row_reset_list.append(row_reset)

    # re-write MERGE.LOG with everything zero'ed out
    INFO_MERGE = { 
        'primary_key' : table_name, 
        'header_line' : comment_lines[0],
        'row_list'    : row_reset_list   
    }

    f = open(merge_file, 'w') 
    write_merge_file(f, INFO_MERGE, [] ) 
    f.close()

    # end merge_table_reset

def compress_files(flag, dir_name, wildcard, name_backup, wildcard_keep ):

    # name of tar file is BACKUP_{name_backup}.tar
    # Inputs
    #  flag > 0 -> compress
    #  flag < 0 -> uncompress
    #  dir_name -> cd to this directory
    #  wildcard -> include these files in tar file
    #  name_backup -> tar file name is BACKUP_{name_backup}.tar
    #  wildcard_keep -> do NOT remove these files
    #

    tar_file   = f"BACKUP_{name_backup}.tar"
    targz_file = f"{tar_file}.gz"
    cddir      = f"cd {dir_name}"

    # be careful if wildcard string is too short; don't want rm *
    if len(wildcard) < 3 :
        msgerr = []
        msgerr = f"wildcard = '{wildcard}' is dangerously short string"
        msgerr = f"that could result in removing too much."
        msgerr = f"Provide longer wildcard string."
        log_assert(False,msgerr)

    if flag > 0 :
        cmd_tar  = f"tar -cf {tar_file} {wildcard} "
        cmd_gzip = f"gzip {tar_file}"

        if len(wildcard_keep) == 0 :
            cmd_rm   = f"rm {wildcard}"
        else:
            # remove all wildcard files EXCEPT for wildcard_keep
            cmd_rm = f"find {wildcard} ! -name '{wildcard_keep}' " + \
                     "-type f -exec rm {} +"

        cmd_all  = f"{cddir} ; {cmd_tar} ; {cmd_gzip} ; {cmd_rm} "
    else:
        cmd_unpack = f"tar -xzf {targz_file}"
        cmd_rm     = f"rm {targz_file}"
        cmd_all    = f"{cddir} ; {cmd_unpack} ; {cmd_rm} "

    #logging.info(f"\n xxx cmd_all = {cmd_all}\n")
    os.system(cmd_all)

    # end compress_files

def compress_subdir(flag,dir_name):

    # flag  > 0 --> tar and gzip dir_name
    # flag  < 0 --> unzip and un-tar
    # Example:
    #   flag=1 and dir_name = $MY_PATH/LOTS_OF_JUNK
    #   --> create $MY_PATH/LOTS_OF_JUNK.tar.gz and remove LOTS_OF_JUNK
    #
    #  flag = -1 restores $MY_PATH/LOTS_OF_JUNK  and removes tar file.
    #
    # Initial use is for cleanup_job_files(flag=1) and 
    # merge_reset(flag=-1)


    topdir_name = os.path.dirname(dir_name)
    subdir_name = os.path.basename(dir_name)

    #logging.info(f" xxx topdir_name = {topdir_name}")
    #logging.info(f" xxx subdir_name = {subdir_name}")

    cddir        = f"cd {topdir_name}"
    tar_file     = f"{subdir_name}.tar"
    targz_file   = f"{tar_file}.gz"

    if flag > 0:  # compress
        cmd_tar    = f"tar -cf {tar_file} {subdir_name}"
        cmd_gzip   = f"gzip {tar_file}"
        cmd_rmdir  = f"rm -rf {subdir_name}"
        cmd_all    = f"{cddir} ; {cmd_tar}; {cmd_gzip} ; {cmd_rmdir}"
        os.system(cmd_all)
    else:  # uncompress if tar file exists
        exist_tar = os.path.exists(f"{topdir_name}/{targz_file}")
        if exist_tar :
            cmd_unpack = f"tar -xzf {targz_file}"
            cmd_rmgz   = f"rm {targz_file}"
            cmd_all    = f"{cddir} ; {cmd_unpack}; {cmd_rmgz} "
            os.system(cmd_all)

    # end compress_subdir

def get_file_lists_wildcard(search_dir, search_wildcard):

    # for input search_wildcard, search for the following file lists:
    #    {search_dir}/{search_wildcard}.LOG
    #    {search_dir}/{search_wildcard}.DONE
    #    {search_dir}/{search_wildcard}.YAML
    # 
    # This function returns these three lists.
    # 
    # To keep associations betwee all three lists, DONE and YAML lists 
    # are forced to have the same length as LOG list. If DONE or YAML
    # don't exist, they are filled as None.
    # Example, there are 4 LOG files; 1st and 3rd have associated
    # DONE file, but only 3rd has YAML file. Output is
    #
    #  log_list  = [ F1.LOG,  F2.LOG, F3.LOG,  F4.LOG ]
    #  done_list = [ F1.DONE, None,   F3.DONE, None ]
    #  yaml_list = [ None,    None,   F3.YAML, None ]
    #
    # Mar 25 2021: find last dot (with .rindex) instead of first dot
    #              to allow for dot in version name.
    #

    # search .LOG first to define list.
    search_log = f"{search_wildcard}.LOG"
    log_list   = sorted(glob.glob1(search_dir, f"{search_log}") )
    done_list = []
    yaml_list = []

    # for each log file, check for DONE and .YAML                          
    for log_file in log_list :
        # xxx mark delete jdot      = log_file.index(".")
        jdot      = log_file.rindex(".")
        prefix    = log_file[0:jdot]
        done_file = f"{prefix}.DONE"
        DONE_FILE = f"{search_dir}/{done_file}"

        yaml_file = f"{prefix}.YAML"
        YAML_FILE = f"{search_dir}/{yaml_file}"
        if not os.path.isfile(DONE_FILE) :
            done_file = None
        if not os.path.isfile(YAML_FILE) :
            yaml_file = None
                
        done_list.append(done_file)
        yaml_list.append(yaml_file)

    return  log_list, done_list, yaml_list
    
    # end get_file_lists_wildcard


def nrow_table_TEXT(table_file, row_key):
    # For input TEXT file, return number rows with 'row_key'
    nrow        = 0
    cmd_grep    = f"grep '{row_key}' {table_file} | wc "
    result_line = subprocess.check_output(cmd_grep,shell=True).rstrip()
    result_line = result_line.decode('utf-8')
    nrow        = int(result_line.split()[0])
    return nrow
    # end nrow_table_TEXT                                                         

def extract_arg(key):
    # If key is of the form  KEY(ARG), function returns ARG.
    # If not (), function returns ''
    arg = ''
    if '(' in key and ')' in key :
        j0  = key.index('(') + 1
        j1  = key.index(')')
        arg = key[j0:j1]
    return arg
    # extract_arg

def parse_done_stamp(output_dir,CONFIG):
    # if DONE_STAMP key is in the CONFIG yaml block, use it;
    # otherwise  return ''
    done_stamp_file = ''
    found = False
    for key in CONFIG_KEYLIST_DONE_FILE :
        if key in CONFIG :
            done_stamp_file = os.path.expandvars(CONFIG[key])
            found = True

    return done_stamp_file,found
    # end parse_done_stamp

# ========================================
def roundup_pow10(ngen):
    # return integer rounded up to nearest power of 10  
    # Examples   
    #   ngen  Nreturn  
    #      0   0
    #     74   100    
    #    423   1000
    #   7842   10000
    
    if ngen == 0 :
        return(0)
    xlog = math.log10(float(ngen))
    nlog = int(xlog) + 1
    x    = math.pow(10.0,float(nlog))
    return int(x)
    # end round_pow10  

# ========================================
def roundup_first_digit(ngen):
    # return integer rounded up based on first digit
    # Examples   
    #   ngen  Nreturn  
    #      0   0
    #     74   80
    #    423   500
    #   7842   8000
    
    if ngen == 0 :
        return(0)
    # comment using example with ngen = 4432
    xgen  = float(ngen)                    # 4432.0
    xlog  = math.log10(xgen)               # = 3.6466
    jlog  = int(xlog)                      # = 3
    xbase = math.pow(10.0,float(jlog))     # = 1000
    x     = xbase*int(xgen/xbase) + xbase  # 5000

    return int(x)
    # end round_first_digit

# ======================================
def copy_input_files(infile_copy_list,output_dir,list_file):

    # copy each file in 'infile_copy_list' to 'output_dir'
    # If list_file != '', write list of files to output_dir/list_file
    # This function avoids duplicate copies so that input list 
    # can contain duplicates.
    #
    # 
    # BEWARE of duplicate input file names in different directories.
    #   If $DIR1/sim_include.input and $DIR2/sim_include.input are
    #   both on the copy list, they would both be copied and clobber
    #   each other. To avoid clobbering, this function aborts if any
    #   base file name appears more than once.
    #   Note that abort also occurs if include file name is 
    #   'sim_include.input' in on job, and $DIRNAME/sim_include.input 
    #   in another job ... even if DIRNAME is CWD.
    #

    done_copy_list        = []
    done_copy_list_nopath = []
    msgerr                = []
    bad_file_name         = False

    for infile in infile_copy_list:
        if infile not in done_copy_list :
            shutil.copy(infile,output_dir)
            done_copy_list.append(infile)

            infile_nopath = os.path.basename(infile)
            # xxx infile_nopath = ntpath.split(infile)[1]
            if infile_nopath in done_copy_list_nopath:
                j = done_copy_list_nopath.index(infile_nopath)
                msgerr.append(f"Cannot define duplicate input/include file names")
                msgerr.append(f"in different paths; see")
                msgerr.append(f"   {infile} ")
                msgerr.append(f"   {infile_copy_list[j]} ")
                log_assert(False,msgerr)

            done_copy_list_nopath.append(infile_nopath)

    # check option to write all input file names to a list file
    # make sure to write only the base name without path, otherwise
    # the merge process will delete the input file from $PATH.

    if list_file != '' :
        LIST_FILE = f"{output_dir}/{list_file}"
        with open(LIST_FILE, 'w') as f : 
            for infile in done_copy_list:
                infile_base = os.path.basename(infile) # exclude path           
                f.write(f"{infile_base}\n")

    # end copy_input_files

def find_duplicates(string_array):
    # returns two arrays
    #   + array of strings that appear more than once
    #   + array if integers, Nappear

    unique_array = set(string_array)
    dupl_list    = []
    n_dupl_list  = []
    for string in unique_array :
        n = string_array.count(string)        
        #print(f" xxx n={n} for string={string}")
        if n > 1:
            dupl_list.append(string)
            n_dupl_list.append(n)

    return n_dupl_list, dupl_list

    # end find_duplicates

def check_file_count(n_expect, file_wildcard):
    # abort if ls {file_wildcard} does not yeidl n_expect files.
    f_list  = glob.glob(f"{file_wildcard}")
    n_list  = len(f_list)
    msgerr  = []
    if n_list != n_expect :
        msgerr.append(f"Found {n_list} file with")
        msgerr.append(f"ls {file_wildcard} ")
        msgerr.append(f"but expected to find {n_expect}")
        log_assert(False,msgerr)
    # end check_file_count

def check_file_exists(file_name,msg_user):
    # abort if file does not exist
    exist = os.path.isfile(file_name)
    if not exist:
        msgerr = [ f"Cannot find file:", f"\t {file_name}" ]
        for msg in msg_user:
            msgerr.append(msg)
        log_assert(exist, msgerr)
    # end check_file_exists

def read_merge_file(merge_file) :
    comment_lines = []
    input_lines   = []
    with open(merge_file, 'r') as f :
        for line in f:
            input_lines.append(line)
            if line[0] == '#' :
                comment_lines.append(line[1:].strip("\n"))

    input_yaml = yaml.safe_load("\n".join(input_lines))
    return input_yaml,comment_lines

    # end read_merge_file

def write_done_stamp(output_dir,done_stamp_list,string):

    # if string is SUCCESS, write it only if stamp file doesn't exist.
    # output_dir is used if there is no slash in done_file.
    for done_file in done_stamp_list :
        if '/' not in done_file :
            DONE_FILE = f"{output_dir}/{done_file}"
        else:
            DONE_FILE = done_file

        if string == STRING_SUCCESS and os.path.isfile(DONE_FILE):
            # if file has FAIL, don't overwrite it. 
            # if file has SUCCESS, no point in re-writing same value.
            pass
        else :
            msg = f"\n Write {string} to done stamp file: \n\t {done_file}"
            logging.info(msg)
            with open(DONE_FILE,"w") as f :
                f.write(f"{string}\n") 

def write_merge_file(f, MERGE_INFO, comment_lines ):
    # write MERGE_INFO to file pointer f
    primary_key = MERGE_INFO['primary_key']
    header_line = MERGE_INFO['header_line']
    row_list    = MERGE_INFO['row_list']    

    f.write(f"#{header_line} \n")
    f.write(f"{primary_key}: \n")
    for row in row_list : 
        f.write(f"  - {row}\n")

    f.write("\n")
    for comment in comment_lines :  f.write(f"#{comment}\n")

    # end write_merge_file

def backup_merge_file(merge_file):
    # Util to debug sequence of updating MERGE.LOG file
    # Input merge_file should include full path.
    Nsec = seconds_since_midnight  # current Nsec, not from submit info
    merge_file_save = f"{merge_file}_{Nsec}"
    shutil.copyfile(merge_file, merge_file_save )
    # end backup_merge_file

def wait_for_files(n_file_wait, wait_dir, wait_files):
    # go to sleep until all wait_files exist
    # Inputs:
    #  n_file_wait     = number of files to wait for
    #  wait_dir        = directory to search for wait_files
    #  wait_files      = file specifier with wildcare; e.g, TMP*.DONE

    T_SLEEP = 20  # sleep time until next file-exist check
 
    logging.info(f"  Wait for {n_file_wait} {wait_files} files")
    n_file_exist = 0; n_try=0
    while  n_file_exist < n_file_wait :
        n_try += 1
        if n_try > 1:  time.sleep(T_SLEEP) 
        wait_file_list  = glob.glob1(wait_dir,wait_files)
        n_file_exist    = len(wait_file_list)
        time_now        = datetime.datetime.now()
        tstr            = time_now.strftime("%Y-%m-%d %H:%M:%S") 
        msg = f"\t Found {n_file_exist} of {n_file_wait} files ({tstr})"
        logging.info(msg)

    # end wait_for_file

def write_job_info(f,JOB_INFO,icpu):

    # write job program plus arguemnts to file pointer f.
    # All job-info params are passed via JOB_INFO dictionary.
    # Jan 8 2021: check optional wait_file

    job_dir      = JOB_INFO['job_dir']    # cd here; where job runs
    program      = JOB_INFO['program']    # name of program; e.g, snlc_sim.exe
    input_file   = JOB_INFO['input_file'] # input file name
    log_file     = JOB_INFO['log_file']   # pipe stdout here
    done_file    = JOB_INFO['done_file']  # DONE stamp for monitor tasks
    arg_list     = JOB_INFO['arg_list']   # argumets of program
    msgerr       = []
    check_abort = JOB_INFO.get(arg_check_abort,False)

    if len(job_dir) > 1 :
        f.write(f"# ---------------------------------------------------- \n")
        f.write(f"cd {job_dir} \n\n")

    CHECK_CODE_EXISTS = '.exe' in program and not check_abort

    CHECK_ALL_DONE    = 'all_done_file' in JOB_INFO  and \
                        'kill_on_fail'  in JOB_INFO  and \
                        not check_abort

    CHECK_WAIT_FILE   = 'wait_file' in JOB_INFO

    if CHECK_ALL_DONE :
        # if ALL.DONE file exists, something else failed ... 
        # so no point in continuing.
        all_done_file = JOB_INFO['all_done_file']  # exists only on failure
        kill_on_fail  = JOB_INFO['kill_on_fail']
        f.write(f"if [ -f {all_done_file} ] ; then \n")
        msg_echo = f"Found unexpected {DEFAULT_DONE_FILE} -> something FAILED."
        f.write(f"  echo '  {msg_echo}' \n")
        if kill_on_fail :
            msg_echo = f"Kill all remaining jobs."
            cmd_kill = f"  cd {CWD}\n"\
                       f"  {sys.argv[0]} \\\n" \
                       f"     {sys.argv[1]} --cpunum {icpu} -k "
                        # f"  exit")
        else:
            msg_echo = "Continue with next job."

        f.write(f"  echo '  {msg_echo}' \n")
        if kill_on_fail: f.write(f"{cmd_kill} \n")
        f.write(f"fi \n\n")

    if CHECK_CODE_EXISTS :
        # wait for program to appear in case SNANA make is in progress

        program_plus_path = find_program(program)

        wait_for_code = f"while [ ! -x {program_plus_path} ]; " \
                        f"do sleep 5; done"
        f.write(f"echo 'Wait for {program} if SNANA make is in progress'\n")
        f.write(f"{wait_for_code}\n")
        f.write(f"echo '{program} exists -> continue' \n\n")

    if CHECK_WAIT_FILE:
        wait_file     = JOB_INFO['wait_file']  # wait for this file to exist
        wait_for_file = f"while [ ! -f {wait_file} ]; " \
                        f"do sleep 10; done"
        f.write(f"echo 'Wait for {wait_file}'\n")
        f.write(f"{wait_for_file}\n")
        f.write(f"echo '{wait_file} exists -> continue' \n\n")

    if check_abort:  # leave human readable marker for each job
        f.write("echo '# ======================================= '\n")

    # check optional ENV to set before running program
    if 'setenv' in JOB_INFO :
        f.write(f"{JOB_INFO['setenv']} \n")

    # check optional start-file stamp (alternate way to get CPU time)
    if 'start_file' in JOB_INFO :
        f.write(f"touch {JOB_INFO['start_file']} \n")
    # - - - - - - - - - 
    f.write(f"{program} {input_file} \\\n")
    # write each arg on separte line for easier viewing
    for arg in arg_list :  
        if arg != '' :
            f.write(f"   {arg} \\\n")
        
    f.write(f"  &>  {log_file} \n" )   # write to stdout and stderr

    # Apr 2022: for lcfit program, remove minuit stdout from log file
    if PROGRAM_NAME_LCFIT in program:
        f.write(f"remove_minuit_stdout.py {log_file}\n")

    if len(done_file) > 4 :
        f.write(f"touch {done_file} \n")
        f.write(f"echo 'Finished {program} -> create {done_file}' \n")

    f.write(f"\n")

    # check for symbolic link(s) from other FITOPTs
    if 'sym_link_list' in JOB_INFO :
        sym_link_list = JOB_INFO['sym_link_list']
        for link in sym_link_list :  
            f.write(f"{link} \n")
        f.write(f"\n")

    # end write_job_info

def find_program(program):

    # Created Oct 11 2021 by R.Kessler
    # use unix "which" to find program. If not found, keep searching
    # every 5 seconds in case code update/make is in progress.
    # If program is not found after very long time -> abort.
    # Function returns full name of program including path.

    t_next        = 5    # check for program this often (sec)
    t_abort       = 300  # abort after this many seconds of searching
    t_sum         = 0    # total time searching for program (sec)

    found_program = False
    while not found_program :
        program_plus_path = shutil.which(program)
        if program_plus_path is None: 
            time_now        = datetime.datetime.now()
            tstr            = time_now.strftime("%Y-%m-%d %H:%M:%S") 
            print(f" Cannot find program {program} at {tstr}; " \
                  f"will try again in {t_next} sec.")
            t_sum += t_next
            if t_sum > t_abort :
                msgerr.append(f"Cannot find program {program}")
                msgerr.append(f"after {t_sum} seconds.");  
                log_assert(False, msgerr) 
            
            time.sleep(t_next)
        else:
            found_program = True                

    return program_plus_path
    # end find_program

def write_jobmerge_info(f,JOB_INFO,icpu):
    # write merge task 
    merge_input_file = JOB_INFO['merge_input_file']
    merge_arg_list   = JOB_INFO['merge_arg_list']
    check_abort      = JOB_INFO.get(arg_check_abort,False)
    match_cpu    = icpu <= NCPU_MERGE_DISTRIBUTE
    do_merge     = len(merge_input_file) > 1  # undefined file -> no merge
    
    if match_cpu and do_merge :
        merge_task = f"{sys.argv[0]} {merge_input_file} {merge_arg_list}"
        f.write(f"cd {CWD} \n")
        f.write(f"echo Run merge_driver monitor task. \n")
        f.write(f"{merge_task} \n")
        if not check_abort: 
            f.write(f"echo $?")
        f.write(f"\n")

    # end write_jobmerge_info

def get_YAML_key_values(YAML_BLOCK, KEYLIST):

    # return keys from CONFIG as follow:
    # YAML_BLOCK:
    #   KEY_BLA:
    #   - bla1
    #   - bla2
    #   - bla3
    #
    # if input KEYLIST = ['KEY_BLA'], 
    # then value_list is returned as ['bla1', 'bla2', 'bla3']
    # If KEYLIST has multiple values, check each value;
    # this feature allow different key-values; 
    #  e.g, KEYLIST = ['INFILE_Ia', 'INFILE_SNIa']
    #
    
    value_list = []
    for key in KEYLIST :
        if key in YAML_BLOCK :
            for item in YAML_BLOCK[key]:
                value_list.append(item)

    return value_list
    # end get_YAML_key_values

def get_survey_info(yaml_path):
    # Read SURVEY (string) and IDSURVEY (int) from YAML file,
    # and return these quantities.
    # If yaml_path is a directory, read first file in glob list;
    # if yaml_path is a file, read this particular file.

    if  os.path.isfile(yaml_path) :
        yaml_file = yaml_path
    else :
        # it's a directory
        yaml_list = glob.glob(f"{yaml_path}/*.YAML")
        yaml_file = yaml_list[0]

    yaml_info = extract_yaml(yaml_file, None, None )    
    return yaml_info['SURVEY'], yaml_info['IDSURVEY']
    # end get_survey_info


# ---------------------------------------------------
# MESSAGING
# ------------------------

class MessageStore(logging.Handler):
    store = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.store = {}

    def emit(self, record):
        l = record.levelname
        if l not in self.store:
            self.store[l] = []
        self.store[l].append(record)

    def get_warnings(self):
        return self.store.get("WARNING", []) + []

    def get_errors(self):
        return self.store.get("CRITICAL", []) + self.store.get("ERROR", []) + self.store.get("EXCEPTION", []) + []

    def print_warnings(self):
        items = self.get_warnings()
        if not items:
            logging.notice("No warnings")
        else:
            logging.warning(f"{len(items)} Warnings:")
            for w in items:
                logging.warning(f"\t{w.msg}")

    def print_errors(self):
        items = self.get_errors()
        if not items:
            logging.notice("No errors")
        else:
            logging.error(f"{len(items)} Errors:")
            for w in items:
                logging.error(f"\t{w.msg}")


def setup_logging(args):
    level = logging.DEBUG if args.verbose else logging.INFO

    # Adding notice level for more important lines that arent warning
    message_store = MessageStore()
    message_store.setLevel(logging.WARNING)
    notice_level = 25
    logging.addLevelName(notice_level, "NOTICE")

    def notice(self, message, *args, **kws):
        if self.isEnabledFor(notice_level):
            self._log(notice_level, message, args, **kws)

    def notice_root(message, *args, **kwargs):
        logging.log(notice_level, message, *args, **kwargs)

    logging.Logger.notice = notice
    logging.notice = notice_root
    fmt = "[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s" if args.verbose else "%(message)s"
    handlers = [logging.StreamHandler(), message_store]
    handlers[0].setLevel(level)
    logging.basicConfig(level=level, format=fmt, handlers=handlers)
    coloredlogs.install(level=level, fmt=fmt, reconfigure=True, level_styles=coloredlogs.parse_encoded_styles("debug=8;notice=green;warning=yellow;error=red,bold;critical=red,inverse"),)
    return message_store


def log_assert(condition, message):
    if not condition:

        msg_abort_face = (
            f"\n\n"
            f"\n   `|```````|`    "
            f"\n   <| o\\ /o |>   "
            f"\n    | ' ; ' |     "
            f"\n    |  ___  |     ABORT submit on Fatal Error. "
            f"\n    | |' '| |     "         
            f"\n    | `---' |     "
            f"\n    \\_______/    " 
            f"\n"
            f"\n{SNANA_ABORT_STRING} : "
        )

        for item in message :
            msg_abort_face += f"\n   {item}"

        logging.exception(message)
        assert condition, msg_abort_face
        # end log_assert

def extract_yaml(input_file, key_start, key_end):

    # Jan 2021: ignore everything before key_start and after key_end.
    #           if either key is None, ignore the key

    msgerr = [(f"Cannot find the input yaml file:\n   {input_file}")]
    exist  = os.path.isfile(input_file)
    log_assert(exist,msgerr)
               
    line_list = []
    FLAG_START = key_start is None
    FLAG_END   = False

    with open(input_file, "r") as f:
        for line in f:
            if key_start is not None:
                if line.startswith(key_start) : FLAG_START = True

            if key_end is not None:
                if line.startswith(key_end) : FLAG_END = True

            #print(f"\t xxx FLAG(START,END) = {FLAG_START} {FLAG_END} " \
            #      f" for line={line}")

            if not FLAG_START : continue
            if FLAG_END: break

            # xxx mark delete if line.startswith("#END_YAML"): break
            line_list.append(line)

    config = yaml.safe_load("\n".join(line_list))

    #logging.info(f" YAML config loaded successfully from {input_file}")
    return config
    # end extract_yaml

