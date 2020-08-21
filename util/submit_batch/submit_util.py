# ==============================================
# Created July 2020 by R.Kessler & S. Hinton
#
# generic utilites for submit  script
# ==============================================

import os, sys, yaml, shutil, glob, math, ntpath
import logging, coloredlogs, subprocess
from   submit_params import *

# =================================================

def get_stat_dict(value_list):
    # For input list of values_list, return dictionary of
    # 'AVG, 'ERR_AVG', 'RMS', 'ERR_RMS'

    n_val    = len(value_list)
    
    if n_val > 0 :
        sumval   = 0.0 ; sqsumval = 0.0
        for val in value_list :
            sumval   += val
            sqsumval += val*val
            
            AVG = sumval / n_val
            RMS = math.sqrt( sqsumval/n_val - AVG*AVG )
            
            ERR_AVG = RMS/math.sqrt(n_val)
            ERR_RMS = ERR_AVG / 1.414   # sigma/sqrt(2*n)
            
    else:
        AVG = 0.0; ERR_AVG = 0.0; RMS=0.0; ERR_RMS=0.0

    stat_dict = { 'AVG':AVG, 'ERR_AVG':ERR_AVG, 'RMS':RMS, 'ERR_RMS': ERR_RMS }
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

def compress_files(flag, dir_name, wildcard, name_backup ):

    # name of tar file is BACKUP_{name_backup}.tar

    #
    tar_file   = (f"BACKUP_{name_backup}.tar")
    targz_file = (f"{tar_file}.gz")
    cddir      = (f"cd {dir_name}")

    # be careful if wildcard string is too short; don't want rm *
    if len(wildcard) < 3 :
        msgerr = []
        msgerr = (f"wildcard = '{wildcard}' is dangerously short string")
        msgerr = (f"that could result in removing too much.")
        msgerr = (f"Provide longer wildcard string.")
        log_assert(False,msgerr)

    if flag > 0 :
        cmd_tar  = (f"tar -cf {tar_file} {wildcard} ")
        cmd_gzip = (f"gzip {tar_file}")
        cmd_rm   = (f"rm {wildcard}")
        cmd_all  = (f"{cddir} ; {cmd_tar} ; {cmd_gzip} ; {cmd_rm} ")
    else:
        cmd_unpack = (f"tar -xzf {targz_file}")
        cmd_rm     = (f"rm {targz_file}")
        cmd_all  = (f"{cddir} ; {cmd_unpack} ; {cmd_rm} ")

    #print(f" xxx cmd_all = {cmd_all}")
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

    jdot        = dir_name.rindex("/")
    topdir_name = dir_name[0:jdot]
    subdir_name = dir_name[jdot+1:]

    #print(f" xxx topdir_name = {topdir_name}")
    #print(f" xxx subdir_name = {subdir_name}")

    cddir        = (f"cd {topdir_name}")
    tar_file     = (f"{subdir_name}.tar")
    targz_file   = (f"{tar_file}.gz")

    if flag > 0:  # compress
        cmd_tar    = (f"tar -cf {tar_file} {subdir_name}")
        cmd_gzip   = (f"gzip {tar_file}")
        cmd_rmdir  = (f"rm -r {subdir_name}")
        cmd_all    = (f"{cddir} ; {cmd_tar}; {cmd_gzip} ; {cmd_rmdir}")
    else:  # uncompress
        cmd_unpack = (f"tar -xzf {targz_file}")
        cmd_rmgz   = (f"rm {targz_file}")
        cmd_all    = (f"{cddir} ; {cmd_unpack}; {cmd_rmgz} ")

    #print(f" xxxx {cmd_all}")
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

    # search .LOG first to define list.
    search_log = (f"{search_wildcard}.LOG")
    log_list   = sorted(glob.glob1(search_dir, f"{search_log}") )
    done_list = []
    yaml_list = []

    # for each log file, check for DONE and .YAML                          
    for log_file in log_list :
        jdot      = log_file.index(".")
        prefix    = log_file[0:jdot]
        done_file = (f"{prefix}.DONE")
        DONE_FILE = (f"{search_dir}/{done_file}")

        yaml_file = (f"{prefix}.YAML")
        YAML_FILE = (f"{search_dir}/{yaml_file}")
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
    cmd_grep    = (f"grep '{row_key}' {table_file} | wc " )
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

            infile_nopath = ntpath.split(infile)[1]
            if infile_nopath in done_copy_list_nopath:
                j = done_copy_list_nopath.index(infile_nopath)
                msgerr.append(f"Cannot define duplicate input/include file names")
                msgerr.append(f"in different paths; see")
                msgerr.append(f"   {infile} ")
                msgerr.append(f"   {infile_copy_list[j]} ")
                log_assert(False,msgerr)

            done_copy_list_nopath.append(infile_nopath)

    # check option to write all input file names to a list file
    if list_file != '' :
        LIST_FILE = (f"{output_dir}/{list_file}")
        with open(LIST_FILE, 'w') as f : 
            for infile in done_copy_list:
                f.write(f"{infile}\n")

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
        msgerr = [ (f"Cannot find file:"), (f"\t {file_name}") ]
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
            DONE_FILE = (f"{output_dir}/{done_file}" )
        else:
            DONE_FILE = done_file

        if string == STRING_SUCCESS and os.path.isfile(DONE_FILE):
            # if file has FAIL, don't overwrite it. 
            # if file has SUCCESS, no point in re-writing same value.
            pass
        else :
            msg = (f"\n Write {string} to done stamp file: \n\t {done_file}")
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

    for comment in comment_lines :
        f.write(f"#{comment}\n")

    # end write_merge_file

def backup_merge_file(merge_file):
    # Util to debug sequence of updating MERGE.LOG file
    # Input merge_file should include full path.
    Nsec = seconds_since_midnight  # current Nsec, not from submit info
    merge_file_save = (f"{merge_file}_{Nsec}")
    shutil.copyfile(merge_file, merge_file_save )
    # end backup_merge_file

def wait_for_files(n_file_wait, wait_dir, wait_files):
    # go to sleep until all wait_files exist
    # Inputs:
    #  n_file_wait     = number of files to wait for
    #  wait_dir        = directory to search for wait_files
    #  wait_files      = file specifier with wildcare; e.g, TMP*.DONE

    logging.info(f"  Wait for {n_file_wait} {wait_files} files")
    n_file_exist = 0
    while  n_file_exist < n_file_wait :
        time.sleep(5)  # sleep time should be param, or passed as arg??
        wait_file_list  = glob.glob1(wait_dir,wait_files)
        n_file_exist    = len(wait_file_list)
        time_now        = datetime.datetime.now()
        tstr            = time_now.strftime("%Y-%m-%d %H:%M:%S") 
        msg = (f"\t Found {n_file_exist} of {n_file_wait} files ({tstr})")
        logging.info(msg)

    # end wait_for_file

def write_job_info(f,JOB_INFO,icpu):

    # write job program plus arguemnts to file pointer f.
    # All job-info are passed via JOB_INFO.

    job_dir    = JOB_INFO['job_dir']    # cd here; where job runs
    program    = JOB_INFO['program']    # name of program; e.g, snlc_sim.exe
    input_file = JOB_INFO['input_file'] # input file name
    log_file   = JOB_INFO['log_file']   # pipe stdout here
    done_file  = JOB_INFO['done_file']  # DONE stamp for monitor tasks
    arg_list   = JOB_INFO['arg_list']   # argumets of program

    if len(job_dir) > 1 :
        f.write(f"# ---------------------------------------------------- \n")
        f.write(f"cd {job_dir} \n\n")

    # for bash, wait for program to appear if SNANA make is in progress.
    # Not sure how to do this in csh.
    if 'bash' in SHELL :
        program_plus_path = shutil.which(program)
        wait_for_code = (f"while [ ! -f {program_plus_path} ]; " \
                         f"do sleep 5; done" )
        f.write(f"echo 'Wait for {program} if SNANA make is in progress'\n")
        f.write(f"{wait_for_code}\n")
        f.write(f"echo {program} exists. \n\n")

    f.write(f"{program} {input_file} \\\n")
    # write each arg on separte line for easier viewing
    for arg in arg_list :  
        if arg != '' :
            f.write(f"   {arg} \\\n")

    f.write(f"   >  {log_file} \n" )

    if len(done_file) > 4 :
        f.write(f"touch {done_file} \n")

    f.write(f"\n")

    # check for symbolic link(s) from other FITOPTs
    if 'sym_link_list' in JOB_INFO :
        sym_link_list = JOB_INFO['sym_link_list']
        for link in sym_link_list :  
            f.write(f"{link} \n")
        f.write(f"\n")

    # end write_job_info

def write_jobmerge_info(f,JOB_INFO,icpu):
    # write merge task 
    merge_input_file = JOB_INFO['merge_input_file']
    merge_arg_list   = JOB_INFO['merge_arg_list']
    match_cpu    = icpu <= NCPU_MERGE_DISTRIBUTE
    do_merge     = len(merge_input_file) > 1  # undefined file -> no merge
    if match_cpu and do_merge :
        merge_task = (f"{sys.argv[0]} {merge_input_file} {merge_arg_list}")
        f.write(f"cd {CWD} \n")
        f.write(f"python {merge_task} \n")
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


def kill_jobs(config_prep):

    # kill jobs and exit 
    submit_mode = config_prep['submit_mode']
    n_core      = config_prep['n_core']

    if ( submit_mode == SUBMIT_MODE_SSH ) :
        # BEWARE: need to check for unique SSH nodes ???
        nodelist = config_prep['nodelist']
        for node in nodelist.split() :
            cmd_kill = ("ssh -x %s 'kill -KILL -1' " % node )
            print('\t %s ' % cmd_kill )
            #            os.system(cmd_kill)

    elif (submit_mode == SUBMIT_MODE_BATCH ):
        pass

    msg = (f"\n Done killing {submit_mode} jobs.")
    sys.exit(msg)
    return

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
            f"\n    |  ___  |     ABORT program on Fatal Error. "
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
#        assert condition, message


def extract_yaml(input_file):
    msgerr = [(f"Cannot find the input yaml file:\n   {input_file}")]
    exist  = os.path.isfile(input_file)
    # xxx mark delete log_assert(os.path.exists(input_file), 
    log_assert(exist,msgerr)
               
    lines = []
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#END_YAML"):
                break
            lines.append(line)
    config = yaml.safe_load("\n".join(lines))

    #logging.info(f" YAML config loaded successfully from {input_file}")
    return config

