# ============================================
# Created July 2020 by R.Kessler & S. Hinton
#
# Base class Program
# ============================================

#import argparse
import os, sys, shutil, yaml
import logging, coloredlogs
import datetime, time, subprocess
import getpass, ntpath, glob

#from   datetime import datetime
from   abc import ABC, abstractmethod
from   submit_params import *

import submit_util as util


# ======================================
class Program:
    def __init__(self, config_yaml, config_prep):

# config_yaml     (user input)
# config_prep     (prepared from yaml input; independent of config_yaml)
# config          (concatenation of the above)

        self.config_yaml   = config_yaml
        self.config_prep   = config_prep
        self.config        = None

        CONFIG = config_yaml['CONFIG']
        if 'JOBNAME' in CONFIG :
            config_prep['program'] = CONFIG['JOBNAME']

        if config_yaml['args'].merge_flag :  # bail for merge process
            return

        logging.info(f" Program name: {config_prep['program']}")

        # unpack BATCH_INFO or NODELIST keys; update config_prep
        self.parse_batch_info(config_yaml,config_prep)

        # end Program __init__

    @abstractmethod       
    def set_output_dir_name(self):
        raise NotImplementedError()

    @abstractmethod       
    def submit_prepare_driver(self):
        raise NotImplementedError()

    @abstractmethod       
    def write_command_file(self, icpu, COMMAND_FILE):
        raise NotImplementedError()

    @abstractmethod       
    def create_merge_table(self,f):
        raise NotImplementedError()

    @abstractmethod       
    def append_info_file(self,f):
        raise NotImplementedError()

    @abstractmethod       
    def merge_job_wrapup(self,job_merge,MERGE_INFO_CONTENTS) :
        raise NotImplementedError()

    @abstractmethod       
    def merge_update_state(self,MERGE_INFO_CONTENTS):
        raise NotImplementedError()

    @abstractmethod       
    def merge_cleanup_final(self):
        raise NotImplementedError()

    @abstractmethod       
    def merge_reset(self,output_dir):
        raise NotImplementedError()

    @abstractmethod       
    def get_merge_COLNUM_CPU(self):
        raise NotImplementedError()

    @abstractmethod
    def get_misc_merge_info(self):
        raise NotImplementedError()

    def parse_batch_info(self,config_yaml,config_prep):
    
        # check of SSH or BATCH, and parse relevant strings

        NODELIST      = ''
        n_core        = 0
        submit_mode   = "NULL"
        node_list     = []
        memory        = MEMORY_DEFAULT
        kill_flag     = config_yaml['args'].kill

        config_prep['nodelist']       = ''
        config_prep['batch_command']  = ''

        CONFIG = config_yaml['CONFIG']

        if 'NODELIST' in CONFIG  :
            NODELIST    = CONFIG['NODELIST']
            n_core      = len(NODELIST.split())
            submit_mode = SUBMIT_MODE_SSH
            config_prep['nodelist'] = NODELIST

            if not kill_flag :
                logging.info(f"\t ssh to {n_core} nodes: {NODELIST} " )
            
            for node in NODELIST.split():
                node_list.append(node)

        elif  'BATCH_INFO' in CONFIG  :
            BATCH_INFO  = CONFIG['BATCH_INFO'].split()
            command     = BATCH_INFO[0]
            template    = os.path.expandvars(BATCH_INFO[1])
            n_core      = int(BATCH_INFO[2])
            node_list   = [HOSTNAME] * n_core  # used for script file name
            submit_mode = SUBMIT_MODE_BATCH
            config_prep['batch_command']  = command
            config_prep['BATCH_TEMPLATE'] = template  # all caps -> full path
            logging.info(f"\t Batch command:  {command}" )
            logging.info(f"\t Batch template: {template}" )
            logging.info(f"\t Batch n_core:   {n_core}" )
        else :
            log_assert(False, f"Could not find BATCH_INFO or NODELIST")

        # check optional memory spec
        if 'BATCH_MEM' in CONFIG :
            memory = CONFIG['BATCH_MEM']
        
        sys.stdout.flush()

        config_prep['n_core']      = n_core 
        config_prep['submit_mode'] = submit_mode
        config_prep['node_list']   = node_list
        config_prep['memory']      = memory
    
    # end parse_batch_info

    def kill_jobs(self):

        # create_output_dir sets dir names, but won't create anything
        # because kill flag is set
        self.create_output_dir()

        output_dir    = self.config_prep['output_dir']
        done_list     = self.config_prep['done_stamp_list']
        submit_mode   = self.config_prep['submit_mode'] 
        batch_command = self.config_prep['batch_command'] 
        IS_SSH        = submit_mode == SUBMIT_MODE_SSH
        IS_BATCH      = submit_mode == SUBMIT_MODE_BATCH

        if IS_SSH :
            self.kill_ssh_jobs()

        elif IS_BATCH and batch_command == 'sbatch' :
            self.kill_sbatch_jobs()

        else:
            msgerr = []
            msgerr.append(f"Unable to kill jobs for:")
            msgerr.append(f"  submit_mode   = '{submit_mode} ")
            msgerr.append(f"  batch_command = '{batch_command}")
            util.log_assert(False,msgerr)

        # - - - - - - - - 
        # if we get here, leave notice in both MERGE.LOG and ALL.DONE files
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        time_now            = datetime.datetime.now()
        with open(MERGE_LOG_PATHFILE, 'a') as f:
            f.write(f"\n !!! JOBS KILLED at {time_now} !!! \n")

        util.write_done_stamp(output_dir, done_list, STRING_STOP)

        # end kill_jobs

    def kill_ssh_jobs(self):
        node_list   = self.config_prep['node_list']
        for node in node_list :
            print(f" Kill jobs on {node}")
            cmd_kill  = "kill -KILL -1"                
            ret       = subprocess.call(["ssh", node, cmd_kill] )

        # end kill_ssh_jobs


    def kill_sbatch_jobs(self):
        
        # kill sbatch jobs in slurm.
        # Read sbatch info from SUBMIT.INFO:
        #   list of CPUNUM, PID, JOB_NAME
        # Then loop over list and execute 'scancel --name=JOBNAME'
        # If cpunum is an argument, kill this cpu last.

        output_dir       = self.config_prep['output_dir']
        INFO_PATHFILE    = (f"{output_dir}/{SUBMIT_INFO_FILE}")
        submit_info_yaml = util.extract_yaml(INFO_PATHFILE)

        # if cpunum is an argument, this cpu is kill last.
        cpunum_last = -9;  job_name_last = ''
        if self.config_yaml['args'].cpunum is not None :
            cpunum_last = self.config_yaml['args'].cpunum[0]

        SBATCH_LIST  = submit_info_yaml['SBATCH_LIST'] 
        njob_kill    = len(SBATCH_LIST)
        for item in SBATCH_LIST :
            cpunum     = item[0]
            pid        = item[1]
            job_name   = item[2]
            if cpunum == cpunum_last : job_name_last = job_name; continue 
            cmd_kill   = (f"scancel --name={job_name}")
            logging.info(f"\t {cmd_kill}")
            os.system(cmd_kill)
            
        # check to kill job on cpunum_last
        if cpunum_last >= 0 :
            cmd_kill   = (f"scancel --name={job_name_last}")
            logging.info(f"\t {cmd_kill}")
            os.system(cmd_kill)

        logging.info(f" Done cancelling {njob_kill} batch jobs.")
        logging.info(f" Check remaining {USERNAME} jobs ... ")
        time.sleep(3)
        os.system(f"squeue -u {USERNAME}")
        logging.info(f"\n Beware that batch jobs might remain from other submits.")

        # end kill_sbatch_jobs

    def write_script_driver(self):

        # For each CPU, create batch script (CPU*.BATCH) which sources
        # command file (CPU*.CMD) that has list of native SNANA commands.
        # For ssh mode, these .CMD files are sourced upon login, and
        # thus .BATCH files are not needed for ssh mode.
        # BATCH files are written here and do not depend on task.
        # CMD files are program-dependent via call to
        #       self.write_command_file(icpu,COMMAND_FILE)
        #
        
        CONFIG      = self.config_yaml['CONFIG']
        input_file  = self.config_yaml['args'].input_file 
        output_dir  = self.config_prep['output_dir']
        script_dir  = self.config_prep['script_dir']
        n_core      = self.config_prep['n_core']
        submit_mode = self.config_prep['submit_mode']
        node_list   = self.config_prep['node_list']
        program     = self.config_prep['program']

        # for each cpu, store name of script and batch file
        command_file_list = []  # command file name, no path
        batch_file_list   = []   # batch file name, no patch
        COMMAND_FILE_LIST = []   # includes full path
        BATCH_FILE_LIST   = []   # idem
        cmdlog_file_list  = []
        job_name_list     = []

        # loop over each core and create CPU[nnn]_JOBS.CMD that
        # are used for either batch or ssh. For batch, also 
        # create batch_file using BATCH_TEMPLATE
        for icpu in range(0,n_core) :
            node          = node_list[icpu]
            cpu_name      = (f"CPU{icpu:04d}")
            job_name      = (f"{input_file}-{cpu_name}") 
            prefix        = (f"CPU{icpu:04d}_JOBLIST_{node}")
            command_file  = (f"{prefix}.CMD")
            log_file      = (f"{prefix}.LOG")
            COMMAND_FILE  = (f"{script_dir}/{command_file}")
            print(f"\t Create {command_file}")

            command_file_list.append(command_file)
            cmdlog_file_list.append(log_file)
            COMMAND_FILE_LIST.append(COMMAND_FILE)
            job_name_list.append(job_name)

            # write few global things to COMMAND_FILE
            with open(COMMAND_FILE, 'w') as f :

                # write linux command to echo start time in python-like
                # format so that python merge process can read it back
                # and measure time pending in batch queue.
                f.write(f"#!/usr/bin/env bash \n")
                f.write(f"echo TIME_START: " \
                        f"`date +%Y-%m-%d` `date +%H:%M:%S` \n")
                f.write(f"echo 'Begin {command_file}' \n\n")

                if STOP_ALL_ON_MERGE_ERROR :
                    f.write(f"set -e \n") 
                #if 'SNANA_LOGIN_SETUP' in CONFIG:
                #    f.write(f"{CONFIG['SNANA_LOGIN_SETUP']} \n")

            # write program-specific content
            self.write_command_file(icpu,COMMAND_FILE)

            # write extra batch file for batch mode
            if ( submit_mode == SUBMIT_MODE_BATCH ):
                batch_file = (f"{prefix}.BATCH")
                BATCH_FILE = (f"{script_dir}/{batch_file}")
                batch_file_list.append(batch_file)
                BATCH_FILE_LIST.append(BATCH_FILE)
                self.write_batch_file(batch_file, log_file,
                                      command_file, job_name)

        # store few thigs for later
        self.config_prep['cmdlog_file_list']  = cmdlog_file_list
        self.config_prep['command_file_list'] = command_file_list
        self.config_prep['COMMAND_FILE_LIST'] = COMMAND_FILE_LIST
        self.config_prep['job_name_list']     = job_name_list
        self.config_prep['batch_file_list']   = batch_file_list
        self.config_prep['BATCH_FILE_LIST']   = BATCH_FILE_LIST

        # make all CMD files group-executable with g+x.
        # Python os.chmod is tricky because it may only apply 
        # permission to user, or wipe out group privs. 
        # To make things easier here, we just use os.system().
        cmd_chmod = (f"chmod g+x {script_dir}/CPU*.CMD")
        os.system(cmd_chmod);

        n_job_tot   = self.config_prep['n_job_tot']
        print(f" BATCH DRIVER JOB COUNT SUMMARY:",
              f"{n_job_tot} {program} jobs on {n_core} cores")

        
        # check option to force crash (to test higher level pipelines)
        if self.config_yaml['args'].force_crash_prep :
            printf(" xxx force batch-prep crash with C-like printf xxx \n")

        # end write_script_driver

    def write_batch_file(self, batch_file, log_file, command_file, job_name):

        # Create batch_file that executes "source command_file"
        # BATCH_TEMPLATE file is read, and lines are modified using
        # string replacements for
        #   REPLACE_NAME, REPLACE_LOGFILE, REPLACE_MEM, REPLACE_JOB
        # with job-specific info.
        # Note that lower-case xxx_file has no path; 
        # upper case XXX_FILE includes full path

        BATCH_TEMPLATE  = self.config_prep['BATCH_TEMPLATE'] 
        script_dir      = self.config_prep['script_dir']
        replace_memory  = self.config_prep['memory']
        debug_batch     = self.config_yaml['args'].debug_batch

        BATCH_FILE      = (f"{script_dir}/{batch_file}")

        with open(BATCH_TEMPLATE,"r") as f:
            template_batch_lines = f.readlines()

        #print(f" xxx template_batch_lines = {template_batch_lines} ")

        # get strings to replace 
        replace_job_name   = job_name

        # nothing to change for log file
        replace_log_file   = log_file  

        #replace_job_cmd = (f"cd {script_dir} \nsource {command_file}")
        replace_job_cmd = (f"cd {script_dir} \nsh {command_file}")

        # - - - define list of strings to replace - - - - 
        batch_lines = []            
        REPLACE_KEY_LIST = [ 'REPLACE_NAME', 'REPLACE_MEM',
                             'REPLACE_LOGFILE', 'REPLACE_JOB' ]
        replace_string_list = [ replace_job_name, replace_memory,
                                replace_log_file, replace_job_cmd ]
        NKEY = len(REPLACE_KEY_LIST)

        # replace strings in batch_lines
        for line in template_batch_lines:
            for ikey in range(0,NKEY):
                REPLACE_KEY    = REPLACE_KEY_LIST[ikey]
                replace_string = replace_string_list[ikey] 
                line           = line.replace(REPLACE_KEY,replace_string)
            batch_lines.append(line)

        with open(BATCH_FILE,"w") as f:
            f.write("".join(batch_lines))

        # end write_batch_file

    def prep_JOB_INFO_merge(self,icpu,ijob):
        # Return JOB_INFO dictionary of strings to run merge process.
        # Inputs:
        #   icpu = 0 to n_core-1
        #   ijob = 1 to n_job_tot
        #
        # Merge task must is the form
        #   python <thisScript> <inputFile> arg_list
        # arg_list includes
        #  -m -> merge and quit or
        #  -M -> wait for all DONE files to appear, then merge it all.
        #  -t <Nsec>   time stamp to verify merge and submit jobs
        #  --cpunum <cpunum>  in case specific CPU needs to be identified    

        input_file     = self.config_yaml['args'].input_file
        no_merge       = self.config_yaml['args'].nomerge
        n_core         = self.config_prep['n_core']
        n_job_tot      = self.config_prep['n_job_tot']
        Nsec           = seconds_since_midnight

        JOB_INFO = {}

        if no_merge :
            JOB_INFO['merge_input_file'] = ""
            JOB_INFO['merge_arg_list']   = ""
            return JOB_INFO

        # - - - - 
        # determine if this is last job for this cpu
        last_job_cpu  = (n_job_tot - ijob) < n_core
 
        # check where to flag last merge process with -M
        if NCPU_MERGE_DISTRIBUTE == 0 :
            # cpu=0 has last merge process
            last_merge =  last_job_cpu and icpu == 0
        else:
            # last merge is after last job, regardless of cpunum
            last_merge =  ijob == n_job_tot

        m_arg = "-m"
        if last_merge :  m_arg = "-M" 

        arg_list = (f"{m_arg} -t {Nsec} --cpunum {icpu}")
        JOB_INFO['merge_input_file']  = input_file
        JOB_INFO['merge_arg_list']    = arg_list
        
        return JOB_INFO

        # prep_JOB_INFO_merge


    def create_info_file(self):

        # create file with permanent information that can be used
        # in MERGE process, or by downstream scrupts.
        # Note that MERGE.LOG changes during processing, but this
        # SUBMIT.INFO file never changes.
        #
        # This function creates file and writes required info.
        # Then it calls append_info_file() to get program-specific info.

        CONFIG          = self.config_yaml['CONFIG']
        n_job_tot       = self.config_prep['n_job_tot']
        n_done_tot      = self.config_prep['n_done_tot']
        n_job_split     = self.config_prep['n_job_split']
        n_core          = self.config_prep['n_core']
        output_dir      = self.config_prep['output_dir']
        script_dir      = self.config_prep['script_dir']
        done_stamp_list = self.config_prep['done_stamp_list']
        Nsec            = seconds_since_midnight
        time_now        = datetime.datetime.now()


        cleanup_flag = 1     # default
        if 'CLEANUP_FLAG' in CONFIG :
            cleanup_flag = CONFIG['CLEANUP_FLAG']  # override deault

        print(f"  Create {SUBMIT_INFO_FILE}")
        INFO_PATHFILE  = (f"{output_dir}/{SUBMIT_INFO_FILE}")
        f = open(INFO_PATHFILE, 'w') 

        # required info for all tasks
        f.write("\n# Required info \n")

        comment = "submit time; Nsec since midnight"
        f.write(f"TIME_STAMP_NSEC:   {Nsec}    # {comment}\n")
        f.write(f"TIME_STAMP_SUBMIT: {time_now}    \n")

        f.write(f"SCRIPT_DIR:       {script_dir} \n")
        f.write(f"DONE_STAMP_LIST:  {done_stamp_list} \n")
        f.write(f"CLEANUP_FLAG:     {cleanup_flag} \n")

        comment = "total number of jobs (excludes sym links)"
        f.write(f"N_JOB_TOT:        {n_job_tot}     # {comment}\n")

        comment = "total number of done files (includes sym links)"
        f.write(f"N_DONE_TOT:       {n_done_tot}     # {comment}\n")

        comment = "njob to merge per task after processing"
        f.write(f"N_JOB_SPLIT:      {n_job_split}     # {comment}\n")

        comment = "number of cores"
        f.write(f"N_CORE:           {n_core}     # {comment} \n")

        force_crash_prep  = self.config_yaml['args'].force_crash_prep
        force_crash_merge = self.config_yaml['args'].force_crash_merge
        f.write(f"FORCE_CRASH_PREP:    {force_crash_prep} \n")
        f.write(f"FORCE_CRASH_MERGE:   {force_crash_merge}\n")
  
        # append program-specific information
        f.write("\n")
        self.append_info_file(f)

        f.close()

        # end create_info_file

    def create_output_dir(self):

        # if there is no slash in output_dir name, tack on cwd
        # to make sure full path is included. Then create output_dir.
        # if output_dir already exists, clobber it.

        kill_flag = self.config_yaml['args'].kill
        output_dir_temp,script_subdir = self.set_output_dir_name()
        
        if '/' not in output_dir_temp :
            output_dir = (f"{CWD}/{output_dir_temp}")
        else:
            output_dir = output_dir_temp

        # define script_dir based on script_subdir
        if ( len(script_subdir) > 2 ):
            script_dir = (f"{output_dir}/{script_subdir}")
        else:
            script_dir = output_dir 

        # store dir names 
        self.config_prep['output_dir'] = output_dir
        self.config_prep['script_dir'] = script_dir

        # fetch & store done stamp file(s)
        CONFIG          = self.config_yaml['CONFIG']
        done_stamp_list = [ DEFAULT_DONE_FILE ]  # always write default ALL.DONE
        done_stamp_file,Found = util.parse_done_stamp(output_dir,CONFIG)
        if Found : done_stamp_list.append(done_stamp_file)  # check for another
        self.config_prep['done_stamp_list'] = done_stamp_list

        # - - - - - - - - - 
        if kill_flag : return
        # - - - - - - - - - 

        logging.info(f" Create output dir:\n   {output_dir}")
        if  os.path.exists(output_dir) :
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)

        # next create subdir for scripts, unless subDir is ./
        if script_dir != output_dir :
            os.mkdir(script_dir)
            
        # end create_output_dir


    def create_merge_file(self):
        output_dir      = self.config_prep['output_dir'] 
        logging.info(f"  Create {MERGE_LOG_FILE}")
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        f = open(MERGE_LOG_PATHFILE, 'w') 

        # create/write the actual tabls(s) ; this is program specific
        self.create_merge_table(f) 

        f.close()

        if KEEP_EVERY_MERGELOG :
            util.backup_merge_file(MERGE_LOG_PATHFILE)
        # end create_merge_file

    def launch_jobs(self):
        
        # submit all the jobs; either batch or ssh
        submit_mode    = self.config_prep['submit_mode']
        script_dir     = self.config_prep['script_dir']
        cddir          = (f"cd {script_dir}")

        if submit_mode == SUBMIT_MODE_BATCH :
            batch_command    = self.config_prep['batch_command'] 
            batch_file_list  = self.config_prep['batch_file_list']
            for batch_file in batch_file_list :
                cmd = (f"{cddir} ; {batch_command} {batch_file}")
                #os.system(cmd)

                ret = subprocess.run( [ batch_command, batch_file], 
                                      cwd=script_dir,
                                      capture_output=True, text=True )
                #print(f" xxx launch {batch_file} -> '{ret.stdout}' ")

            self.fetch_slurm_pid_list()

        else:
            n_core         = self.config_prep['n_core']
            node_list      = self.config_prep['node_list']
            command_file_list = self.config_prep['command_file_list']
            cmdlog_file_list  = self.config_prep['cmdlog_file_list'] 
            CONFIG            = self.config_yaml['CONFIG']
            if 'SNANA_LOGIN_SETUP' in CONFIG:
                login_setup = (f"{CONFIG['SNANA_LOGIN_SETUP']}")
            else:
                login_setup = ""

            logging.info(f"  login_setup (for ssh):  {login_setup} ")
            qq = '"'
            for inode in range(0,n_core):
                node       = node_list[inode]
                log_file   = cmdlog_file_list[inode]
                cmd_file   = command_file_list[inode]

                #cmd_source = (f"{cddir} ; source {cmd_file} >& {log_file} &")
                cmd_source = (f"{login_setup}; {cddir} ; " \
                              f"sh {cmd_file} >& {log_file} &")

                #ret = subprocess.call(["ssh", node, cmd_source] )
                logging.info(f"  Submit jobs via ssh -x {node}")
                ret = subprocess.Popen(["ssh", "-x", node, cmd_source ],
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE)
                #print(f" xxx {node} ret = {ret}")

        # end launch_jobs

    def fetch_slurm_pid_list(self):

        # for sbatch, fetch process id for each CPU; otherwise do nothing.
        # 
        batch_command    = self.config_prep['batch_command']
        output_dir       = self.config_prep['output_dir']
        script_dir       = self.config_prep['script_dir']
        job_name_list    = self.config_prep['job_name_list']
        msgerr = []
        if batch_command != 'sbatch' : return

        # prep squeue command with format: i=pid, j=jobname            
        cmd = (f"squeue -u {USERNAME} -h -o '%i %j' ")
        ret = subprocess.run( [cmd], shell=True, 
                              capture_output=True, text=True )
        pid_all = ret.stdout.split()
        
        INFO_PATHFILE  = (f"{output_dir}/{SUBMIT_INFO_FILE}")
        f = open(INFO_PATHFILE, 'a') 
        f.write(f"\nSBATCH_LIST:  # [CPU,PID,JOB_NAME] \n")

        for job_name in job_name_list :
            j_job    = pid_all.index(job_name)
            if j_job <= 0 :
                msgerr.append(f"Could not find pid for job = {job_name}")
                msgerr.append(f"Check for sbatch problem")
                self.log_assert(False, msgerr)

            pid      = pid_all[j_job-1]
            cpunum   = int(job_name[-4:])
            print(f"\t pid = {pid} for {job_name}")
            f.write(f"  - [ {cpunum:3d}, {pid}, {job_name} ] \n")
        f.close()
        
        # end fetch_slurm_pid_list

    def merge_driver(self):

        # Called when -m or -M argument is passed to submit_batch_jobs.
        # Read & write MERGE.LOG file with SPLIT and MERGE tables
        # in YAML format.
        # Program-specific functions modify the SPLIT and MERGE
        # tables. With each change in processing status, the
        # MERGE.LOG file is modified ... so be careful.
        # 
        # Each SPLIT job and MERGE job has a state. The state
        # sequence is
        #     WAIT -> RUN -> DONE [-> FAIL]
        # where DONE indicates success, and FAIL is described
        # in the failure_* functions.
        #
        # After RUN has finished, DONE or FAIL is determined by
        # program-specific analysis of the stdout (log file).
        # The key function is 'merge_update_state', which figures
        # out which jobs have changed state, and also takes 
        # appropriate action such as moving and removing files.

        # - - - - - -
        # collect time/date info to clearly mark updates in
        # CPU*.LOG file(s)
        Nsec     = seconds_since_midnight  # current Nsec, not from submit info
        time_now = datetime.datetime.now()
        tstr     = time_now.strftime("%Y-%m-%d %H:%M:%S") 
        fnam = "merge_driver"
        MERGE_LAST  = self.config_yaml['args'].MERGE_LAST
        cpunum      = self.config_yaml['args'].cpunum[0]

        logging.info(f"\n")
        logging.info(f"# ================================================== ")
        logging.info(f"# {fnam}: Begin at {tstr} ({Nsec})")

        # need to re-compute output_dir to find submit info file
        output_dir,script_subdir = self.set_output_dir_name()
        self.config_prep['output_dir']  = output_dir

        # read SUBMIT.INFO passed from original submit job... 
        # this info never changes
        logging.info(f"# {fnam}: read {SUBMIT_INFO_FILE}")
        INFO_PATHFILE    = (f"{output_dir}/{SUBMIT_INFO_FILE}")
        submit_info_yaml = util.extract_yaml(INFO_PATHFILE)
        self.config_prep['submit_info_yaml'] = submit_info_yaml

        # check option to reset merge process 
        if self.config_yaml['args'].merge_reset :
            logging.info(f"# {fnam}: Execute merge_reset debug utility")
            self.merge_reset(output_dir)
            sys.exit("\n Done with merge_reset. " \
                     f"Try merge again with -m option.")

        # make sure DONE_STAMP is stored in both config_prep and info_yaml
        # in case a prep function is called again during merge process.
        self.config_prep['done_stamp_list'] = submit_info_yaml['DONE_STAMP_LIST']
        self.config_prep['cleanup_flag']    = submit_info_yaml['CLEANUP_FLAG']

        # make sure time stamps are consistent
        self.merge_check_time_stamp(output_dir)

        # if last merge call (-M), then must wait for all of the done
        # files since there will be no more chances to merge.
        if MERGE_LAST : self.merge_last_wait()

        # set busy lock file to prevent a simultaneous  merge task
        self.set_merge_busy_lock(+1)

        # read status from MERGE file. There is one comment line above
        # each table to provide a human-readable header. The remaining
        # comments after the tables are probably merge-abort messages; 
        # these post-table comments are saved in comment_lines, and 
        # re-written below in write_merge_file function.

        logging.info(f"# {fnam}: examine {MERGE_LOG_FILE}")
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        MERGE_INFO_CONTENTS,comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        self.merge_config_prep(output_dir)  # restore config_prep

        # make boolean list of which MERGED processes have already finished;
        # below will check which processes are newly DONE.
        done_list_start = self.get_merge_done_list(3,MERGE_INFO_CONTENTS)

        # check for changes in state for SPLIT and MERGE tables.
        # function returns updated set of SPLIT and MERGE table rows.
        row_list_split, row_list_merge, n_change = \
                   self.merge_update_state(MERGE_INFO_CONTENTS)
        
        # check option to force crash (to test higher level pipelines)
        if submit_info_yaml['FORCE_CRASH_MERGE'] and cpunum == 0 :
            printf(" xxx force merge crash with C-like printf xxx \n")

        use_split = len(row_list_split) > 0
        use_merge = len(row_list_merge) > 0

        # Modify MERGE.LOG if there is a change in the processing STATE
        if n_change > 0 :
            itable = 0 
            # redefine SPLIT table
            if use_split :
                INFO_STATE_SPLIT = {
                    'header_line' : comment_lines[itable],
                    'primary_key' : TABLE_SPLIT,
                    'row_list'    : row_list_split }
                itable += 1

            if use_merge :
                INFO_STATE_MERGE = {
                    'header_line' : comment_lines[itable],
                    'primary_key' : TABLE_MERGE, 
                    'row_list'    : row_list_merge }
                itable += 1

            # re-write MERGE.LOG file with new set of STATEs
            # Note that SPLIT table is optional; MERGE table is required.
            # Any comment_lines after the tables are re-written so that
            # we don't lose information or merge-abort messages.
            with open(MERGE_LOG_PATHFILE, 'w') as f :
                if use_split :
                    util.write_merge_file(f, INFO_STATE_SPLIT, [] )
                if use_merge :
                    util.write_merge_file(f, INFO_STATE_MERGE, \
                                          comment_lines[itable:] )

            msg_update = (f"Finished {n_change} STATE updates ({Nsec}).")
        else:
            msg_update = (f"No merge updates -> do nothing. ")

        logging.info(f"# {fnam}: {msg_update}")

        # -----------------------------------
        # for debug, keep each MERGE.LOG as MERGE.LOG_{Nsec}
        if KEEP_EVERY_MERGELOG :
            util.backup_merge_file(MERGE_LOG_PATHFILE)

        # call wrap-up function for each newly MERGED process in MERGE table
        done_list_end = self.get_merge_done_list(3, MERGE_INFO_CONTENTS)
        fail_list_end = self.get_merge_done_list(2, MERGE_INFO_CONTENTS)
        n_job_merge   = len(MERGE_INFO_CONTENTS[TABLE_MERGE]) # entire list
        n_done        = 0  
        n_wrapup      = 0
        for job_merge in range(0,n_job_merge):
            done_now     = done_list_end[job_merge]  # DONE or FAIL
            fail_now     = fail_list_end[job_merge]  # FAILED
            done_before  = done_list_start[job_merge] 
            done_new     = done_now and not done_before
            if done_now :
                n_done += 1
            if done_new and not fail_now :
                n_wrapup += 1
                self.merge_job_wrapup(job_merge,MERGE_INFO_CONTENTS)

        logging.info(f"# {fnam}: finished {n_wrapup} wrapup tasks ")

        # Only last merge process does cleanup tasks and DONE stamps
        if MERGE_LAST and n_done == n_job_merge :
            nfail_tot = self.failure_summary()

            #self.merge_write_misc_info()     # write task-specific info
            misc_info     = self.get_misc_merge_info()
            proctime_info = self.get_proctime_info() 
            self.append_merge_file(misc_info+proctime_info)

            if nfail_tot == 0 :
                logging.info(f"\n# {fnam}: ALL JOBS DONE -> " \
                             f"BEGIN FINAL CLEANUP ")
                cleanup_flag  = submit_info_yaml['CLEANUP_FLAG']
                if cleanup_flag:  self.merge_cleanup_final()  
                STRING_STATUS = STRING_SUCCESS
            else:
                STRING_STATUS = STRING_FAIL

            done_stamp_list = submit_info_yaml['DONE_STAMP_LIST']
            util.write_done_stamp(output_dir, done_stamp_list, STRING_STATUS)

            logging.info(f"\n# {fnam}: finished with {STRING_STATUS}. " \
                         f"Bye Bye" )

        # remove busy lock file
        self.set_merge_busy_lock(-1)
        
        # end merge_driver

    def merge_last_wait(self):

        # called after last science job (jobid = n_job_tot), but beware 
        # that other jobs may run later due to pending in queue. 
        # Here we wait for
        #  + all DONE files to exist
        #  + no BUSY files from other merge process
        # 

        submit_info_yaml = self.config_prep['submit_info_yaml'] 
        n_job_tot        = submit_info_yaml['N_JOB_TOT']
        n_done_tot       = submit_info_yaml['N_DONE_TOT']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_dir       = submit_info_yaml['SCRIPT_DIR'] 
        done_wildcard    = (f"{jobfile_wildcard}.DONE")
        util.wait_for_files(n_done_tot, script_dir, done_wildcard) 

        time.sleep(1)

        # sleep until there are no more busy files.
        n_busy,busy_list = self.get_busy_list()
        while n_busy > 0 :
            logging.info(f"\t Wait for {busy_list} to clear")
            time.sleep(5)
            n_busy,busy_list = self.get_busy_list()

        # end merge_last_wait

    def merge_cleanup_script_dir(self):
        # Tar of CPU* files, then tar+gzip script_dir
        # Beware that after this task, nothing more can be written
        # to the CPU*LOG files.

        cpunum            = self.config_yaml['args'].cpunum[0]
        submit_info_yaml  = self.config_prep['submit_info_yaml']
        script_dir        = submit_info_yaml['SCRIPT_DIR'] 

        log_file_keep  = (f"CPU{cpunum:04d}*.LOG")
        logging.info(f"  Standard cleanup: compress CPU* files " \
                     f"except for {log_file_keep}")
        util.compress_files(+1, script_dir, "CPU*", "CPU", log_file_keep )

        # tar and zip the script dir
        logging.info(f"  Standard cleanup: compress {script_dir}")
        util.compress_subdir(+1, f"{script_dir}" )

        # merge_cleanup_script_dir

    def set_merge_busy_lock(self,flag):
        # flag > 0 --> create BUSY LOCK file
        # flag < 0 --> remove it

        Nsec  = seconds_since_midnight  # current Nsec, not from submit info
        t_msg = (f"T_midnight={Nsec}")

        output_dir     = self.config_prep['output_dir']

        if self.config_yaml['args'].cpunum :
            cpunum = self.config_yaml['args'].cpunum[0] # passed from script
        else:
            return 
            #cpunum = 0  # interactive debug

        busy_file     = (f"{BUSY_FILE_PREFIX}{cpunum:04d}.{BUSY_FILE_SUFFIX}")
        BUSY_FILE     = (f"{output_dir}/{busy_file}")
        # xxx busy_wildcard = (f"{BUSY_FILE_PREFIX}*.{BUSY_FILE_SUFFIX}")

        if flag > 0 :
            # check for other busy files to avoid conflict
            n_busy,busy_list = self.get_busy_list()
            if n_busy > 0 :
                msg = (f"\n Found existing {busy_list[0]} --> "\
                       f"exit merge process.")
                sys.exit(msg)  
            else: 
                logging.info(f"\t Create {busy_file} for {t_msg}")
                with open(BUSY_FILE,"w") as f:
                    f.write(f"{Nsec}\n")  # maybe useful for debug situation

                # check for simultaneously created busy file(s);
                # --> avoid conflict by keeping only first one in sort list
                time.sleep(1)
                n_busy,busy_list = self.get_busy_list()
                if n_busy > 1 and busy_file != busy_list[0] :
                    cmd_rm = (f"rm {BUSY_FILE}")
                    os.system(cmd_rm)
                    msg = (f"\n Found simultaneous {busy_list[0]} --> "\
                           f" exit merge process.")
                    sys.exit(msg)  

        elif len(BUSY_FILE) > 10 :  # avoid disaster with rm
            logging.info(f"\t Remove {busy_file} for {t_msg}")
            cmd_rm = (f"rm {BUSY_FILE}")
            os.system(cmd_rm)

        # end set_merge_busy_lock

    def get_busy_list(self):
        # return numbrer of busy files,  and sorted list of busy files.
        output_dir     = self.config_prep['output_dir']
        busy_wildcard  = (f"{BUSY_FILE_PREFIX}*.{BUSY_FILE_SUFFIX}")
        busy_list  = sorted(glob.glob1(output_dir,busy_wildcard))
        n_busy     = len(busy_list)
        return n_busy, busy_list
        # end get_busy_list

    def get_merge_done_list(self, mask_check, MERGE_INFO_CONTENTS):

        # return list of logical flags per MERGED GENVERSION, 
        # where flag=True if Done, flag=false if not done. 
        # Note that this operates on final GENVERSION and  
        # NOT on the TMP-GENVERSIONs. 
        #
        # mask_check = 1 : return True if DONE
        # mask_check = 2 : return True if FAIL
        # mask_check = 3 : return True if DONE or FAIL

        row_list         = MERGE_INFO_CONTENTS[TABLE_MERGE]
        n_genversion     = len(row_list)
        iver_done_flag   = [ False ] * n_genversion
        iver_all         = 0

        check_DONE = (mask_check & 1) > 0
        check_FAIL = (mask_check & 2) > 0

        for row in row_list :
            STATE     = row[COLNUM_MERGE_STATE]
            is_done   = (STATE == SUBMIT_STATE_DONE )
            is_fail   = (STATE == SUBMIT_STATE_FAIL )
            if check_DONE and is_done :
                iver_done_flag[iver_all] = True
            if check_FAIL and is_fail :
                iver_done_flag[iver_all] = True

            iver_all += 1

        return iver_done_flag
        # end get_merge_done_list   
    
    def merge_check_time_stamp(self,output_dir):
        # compare time stamp in SUBMIT.INFO file against time
        # stamp passed to this merge process; if they don't match,
        #   + create STOP-MERGE_{Nsec_now} file (with comments inside)
        #   + create DONE_STAMP with STOP written inside
        #   + exit this script with error code
        #
        # Most likely explanation for time-stamp mis match is 
        # re-submitting without killing previous jobs
        #
        # if 'set -e' is in the CPU*.CMD script, then stoppig
        # here will also stop the entire CPU*.CMD script so that
        # the next jobs are not run.

        submit_info_yaml  = self.config_prep['submit_info_yaml']
        #output_dir        = self.config_prep['output_dir']
        time_stamp_submit = submit_info_yaml['TIME_STAMP_NSEC']
        done_stamp_list   = submit_info_yaml['DONE_STAMP_LIST']
        Nsec_now          = seconds_since_midnight

        # allow interactive debug without -t so that we don't have
        # to manually extract the new time stamp from SUBMIT.INFO 
        # file to debug merge process.
        # This allows debugging with
        #   submit_batch_jobs.py <inFile> -m
        if not self.config_yaml['args'].t :
            return

        # time_stamp and cpunum args are given only in the CPU*CMD files
        time_stamp_merge  = self.config_yaml['args'].t[0]
        cpunum            = self.config_yaml['args'].cpunum[0]

        if time_stamp_submit != time_stamp_merge :
            stop_file = (f"{output_dir}/STOP_CPU{cpunum:04d}_{Nsec_now}.INFO")
            with open(stop_file,"w") as f :
                f.write(f"Time stamp mis-match: \n")
                f.write(f"  T_midnight(submit) = {time_stamp_submit} sec \n")
                f.write(f"  T_midnight(merge)  = {time_stamp_merge} sec\n")
                f.write(f"\n")
                f.write(f"STOP MERGE and force CPU script to exit.\n")
                f.write(f"\n")
                f.write(f" Likely explanation: job re-submitted without " \
                        f"killing previous batch jobs.\n")

            for done_file in done_stamp_list:
                with open(done_file,"w") as f :
                    f.write(STRING_STOP) 

            msgerr = []
            msgerr.append(f"Time stamp mis-match in merge process")           
            msgerr.append(f"See info in file {stop_file}")
            util.log_assert(False,msgerr)
            #exit(EXIT_CODE_TIME_STAMP)
            
        # end merge_check_time_stamp    
    
    def append_merge_file(self,info_list):
        # append extra keys to merge log file
        # Each info line is of the form
        #   "KEY: VALUE"

        output_dir          = self.config_prep['output_dir']
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")   
        
        print(f" xxx MERGE_LOG_PATHFILE = {MERGE_LOG_PATHFILE} ")

        #  append to bottom of MERGE.LOG
        with open(MERGE_LOG_PATHFILE,"a") as f:
            for line in info_list:  f.write(f"{line} \n")

        # end append_merge_file

    def get_proctime_info(self):
        # return proc time info for MERGE.LOG file including
        #  + wall time
        #  + EFF(CPU) = average CPU/core divided by wall-time 
        #
        # Note that time_xxx are datetime objects;
        # t_xxx is a computed time (float)

        submit_info_yaml  = self.config_prep['submit_info_yaml']
        script_dir        = submit_info_yaml['SCRIPT_DIR'] 
        n_core            = submit_info_yaml['N_CORE']
        time_submit       = submit_info_yaml['TIME_STAMP_SUBMIT']
        time_now          = datetime.datetime.now()
        time_dif          = time_now - time_submit

        output_dir          = self.config_prep['output_dir']
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")   
        
        t_seconds = time_dif.total_seconds()
        if t_seconds < 3000.0 :
            t_unit = 60.0;      unit = "minutes"
        else:
            t_unit = 3600.0 ;    unit = "hours"

        t_wall   = t_seconds/t_unit
        msg_time = [ ]
        msg_time.append(f"UNIT_TIME:      {unit} ")
        msg_time.append(f"WALL_TIME:      {t_wall:.2f}  ")

        # - - - - - - - - 
        # Read TIME_START value from each CPU*LOG file, and measure
        # min & max time pending in batch queue. First or 2nd line of 
        # CPU*LOG should be of the form
        #         TIME_START: YYYY-MM-DD HH:MM:SS
        # Ideally this key would be read as YAML file, but unfortunately
        # the batch systems write non-YAML output at the top of the 
        # CPU*LOG files. The TIME_START argument has to be read as a
        # string and converted to a datetime object ... hence the ugly code.
        cpu_wildcard  = "CPU*.LOG"
        cpu_log_list  = sorted(glob.glob1(script_dir,cpu_wildcard))
        t_pend_list   = []
        for cpu_log_file in cpu_log_list :            
            LOG_FILE = (f"{script_dir}/{cpu_log_file}")
            found_time = False
            #print(f" xxx -------------------------------------- ")
            #print(f" xxx examine {cpu_log_file}")
            with open(LOG_FILE,'r') as f :
                for line in f :
                    word_list = (line.rstrip("\n")).split()
                    if len(word_list) == 0 : continue 
                    key       = word_list[0]
                    if key == "TIME_START:" :
                        found_time       = True
                        t_str            = (f"{word_list[1]} {word_list[2]}")
                        time_start_cpu   = \
                            datetime.datetime.strptime(t_str, '%Y-%m-%d %H:%M:%S')  
                        t_sec = (time_start_cpu - time_submit).total_seconds()
                        t_sec         += 1.0  # avoid violating causality
                        t_pend         = t_sec/t_unit
                        t_pend_list.append(t_pend)
                    if found_time : break

        if len(t_pend_list) > 0 :
            t_pend_min = min(t_pend_list)
            t_pend_max = max(t_pend_list)
            msg_time.append(f"TMIN_PENDING:   {t_pend_min:.02f} ")
            msg_time.append(f"TMAX_PENDING:   {t_pend_max:.02f} ")

        # - - - - - - - - 
        # if there is a CPU column, compute total CPU and avg CPU/core/T_wal
        colnum_cpu = self.get_merge_COLNUM_CPU()
        if colnum_cpu > 0 :
            cpu_sum    = 0.0
            MERGE_INFO_CONTENTS,comment_lines = \
                util.read_merge_file(MERGE_LOG_PATHFILE)
            for row in MERGE_INFO_CONTENTS[TABLE_MERGE] :
                cpu = row[colnum_cpu]
                cpu_sum += cpu  # units are minutes, not seconds

            cpu_sum *= 60.0/t_unit
            cpu_avg  = cpu_sum / n_core
            eff_cpu  = cpu_avg / t_wall

            msg_time.append(f"CPU_SUM:        {cpu_sum:.3f} ")
            msg_time.append(f"EFFIC_CPU:      {eff_cpu:.3f}   # CPU/core/T_wall")


        return msg_time

        # end get_proctime_info

    def check_file_exists(self,file_name,msg_user):
        # abort if file does not exist and use self.log_assert
        # that writes error message to MERGE.LOG and FAIL to done stamp
        exist = os.path.isfile(file_name)
        if not exist:
            msgerr = [ (f"Cannot find file:"), (f"\t {file_name}") ]
            for msg in msg_user:
                msgerr.append(msg)
            self.log_assert(exist, msgerr)
    # end check_file_exists

    def check_for_failure(self, log_file, nevt, isplit):

        # Top-level function to check for failures. External
        # program must determine nevt.
        # nevt <=0 means job failed, so call faulure_update
        # function that will update FAILURES.LOG
        # This function returns true if failure is identified.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        cleanup_flag     = submit_info_yaml['CLEANUP_FLAG']
        MXSPLIT_FAIL_REPEAT = 2       # max isplit to make FAIL_REPEAT script
        
        fail_no_output   = (nevt <  0)  # no output
        fail_zero_evt    = (nevt == 0)  # zero events 
        found_fail       = (fail_no_output or fail_zero_evt)
        create           = (isplit <= MXSPLIT_FAIL_REPEAT )
        
        if found_fail :
            cleanup_flag   = 0 # STOP all cleanup activities
            self.failure_update(log_file, create, fail_zero_evt )

        return found_fail

        # end check_for_failure
    
    def failure_update(self, job_log_file, create_script, found_zero):

        # Generic diagnostic for jobs that do not finish successfully.
        # Inspect job_log_file in log_dir to determine if
        #  + job aborted
        #  + job crashed (no valid end keys)
        #  + zero events  (if inoput found_zero=True)
        #
        # If there is a failure, update FAILURES.LOG with yaml block;
        # Each YAML block key is
        #   FAILURE-{num_fail}
        #
        # Inputs:
        #   job_log_file  : file to search for ABORT message
        #   create_script : True to create script to repeate failure
        #   found_zero    : T -> job finished ok, but zero events processed.
        #
        # Beware that fail-stats are not stored in memory because the 
        # merge/monitor tasks exit and re-start. Therefore, failure 
        # stats are re-read from FAILURES.LOG each time this function
        # is called. Similarly, all of the CMD command lines are re-read.
        # The file-reading is excessive, but only when there are lots 
        # of failures. 

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR'] 
        msgerr = []
        NMSG_ABORT_SHOW  = 3  # show 3 msg lines from log file

        FAIL_SUMMARY_PATHFILE  = (f"{script_dir}/{FAIL_SUMMARY_FILE}")
        JOB_LOG_PATHFILE       = (f"{script_dir}/{job_log_file}")

        # on 1st failure, create fail log
        if not os.path.isfile(FAIL_SUMMARY_PATHFILE):
            msg = (f"\n FIRST JOB FAILURE -> Create {FAIL_SUMMARY_FILE}")
            logging.info(msg)
            with open(FAIL_SUMMARY_PATHFILE,"w") as f:
                pass  # create fail-log, but don't write 

        # - - - - - - - -
        # for each new merge process, re-read from info from disk.
        if 'n_job_fail_list' not in self.config_prep :
            n_job_fail_list = self.read_failure_stats()
            command_lines   = self.read_command_lines()
            self.config_prep['n_job_fail_list'] = n_job_fail_list
            self.config_prep['command_lines']   = command_lines

        # load local pointers to info from files
        n_job_fail_list = self.config_prep['n_job_fail_list']
        command_lines   = self.config_prep['command_lines']

        # - - - - - - - - -
        # construct name of FAIL-REPEAT script.
        # Prepend FAIL_REPEAT and replace .LOG with .CMD.
        j_dot       = job_log_file.index(".")
        j_start     = 0

        # special case for SIM jobs: remove TMP_[USER4] prefix
        sim_prefix = (f"TMP_{USER4}_")
        if sim_prefix in job_log_file:
            j_start = len(sim_prefix)
        
        fail_repeat_script = (f"FAIL_REPEAT_{job_log_file[j_start:j_dot]}.CMD")
        FAIL_REPEAT_SCRIPT = (f"{script_dir}/{fail_repeat_script}")

        # - - - - - - - - - - - - - -  
        # Check log file for abort.
        msg_abort  = []
        nmsg_abort = -9
        with open(JOB_LOG_PATHFILE) as f_log :
            #for i, line in enumerate(f_log.read().splitlines()):
            for line in f_log:
                if SNANA_ABORT_STRING in line:
                    nmsg_abort = 0
                if nmsg_abort >= 0 :
                    msg_abort.append(line.rstrip("\n"))
                    nmsg_abort += 1

        # check logicals to find out what kind of failure
        FAIL_INDEX = -9 
        if nmsg_abort > 0 :
            FAIL_INDEX = FAIL_INDEX_ABORT
        elif found_zero :
            FAIL_INDEX = FAIL_INDEX_ZERO
        else :
            FAIL_INDEX = FAIL_INDEX_UNKNOWN

        FAIL_MODE  = FAIL_MODE_STRINGS[FAIL_INDEX]

        # increment total number of failures
        n_job_fail_list[0] += 1
        n_job_fail          = n_job_fail_list[0]
        msg = (f"    *** FAIL {n_job_fail:3d} for " \
               f"{job_log_file} ")
        logging.info(msg)

        # - - - - - - - - - - - - - - -
        # open FAILURES LOG and append this failure with new YAML block
        with open(FAIL_SUMMARY_PATHFILE,"a") as f:
            f.write(f"\n")
            f.write(f"FAILURE-{n_job_fail:04d}:\n")
            f.write(f"  JOB_LOG_FILE:       {job_log_file} \n")

            script_arg = fail_repeat_script
            if create_script is False :
                script_arg = "NOT_CREATED"
            f.write(f"  FAIL_REPEAT_SCRIPT: {script_arg} \n")

            f.write(f"  FAIL_MODE:          {FAIL_MODE} \n")

            if len(msg_abort) > 0 :
                f.write(f"  ABORT_MESSAGES: \n")
                for i in range(0,NMSG_ABORT_SHOW):
                    f.write(f"  - {msg_abort[i].lstrip()} \n")

        # - - - - - - - - - - - - - - - -
        # search full command in original CMD files needed to
        # repeat failure
        program       = self.config_prep['program']  # e.g., snlc_sim.exe
        lnum_program  = -1  # last line number containing program name
        lnum_log      = -1  # line number containing log file name
        lnum = 0  # keep track of line numberr
        for line in command_lines:
            if program in line :
                lnum_program = lnum
            if job_log_file in line :
                lnum_log = lnum
                break
            lnum += 1

        #print(f"\t xxx lnum[prog,log] = {lnum_program} {lnum_log}")

        msgerr.append(f"Could not find log file name in " \
                      f"{CMD_wildcard} scripts")
        msgerr.append(f"Log file is {job_log_file}")
        msgerr.append(f"{CMD_wildcard} scripts are in\n\t {script_dir}")
        util.log_assert(lnum_log > 0, msgerr)

        # write full repeat-failure command to bash script FAIL_xxx.CMD
        if create_script :
            with open(FAIL_REPEAT_SCRIPT,"w") as f :
                f.write("\n".join(command_lines[lnum_program:lnum_log]) )
                os.system(f"chmod +x {FAIL_REPEAT_SCRIPT}")

    # end  failure_update

    def read_failure_stats(self):

        # read FAILURES.LOG file and count total number of
        # failures, and also now many failures per fail-mode.
        # Read lines for the form
        #   FAIL_MODE:  [MODE]
        #
        # where [MODE] = ABORT, ZERO_EVENTS, or UNKNOWN
        # This read function allows independent merge/monitor 
        # tasks to catch up on the stats from previous tasks.

        n_fail_mode = len(FAIL_MODE_STRINGS)  # includes TOTAL
        n_job_fail_list = [ 0 ] * n_fail_mode

        #output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR'] 
        FAIL_SUMMARY_PATHFILE = (f"{script_dir}/{FAIL_SUMMARY_FILE}")

        # ?? is there a better way to parse this ???
        with open(FAIL_SUMMARY_PATHFILE,"r") as f :
            for line in f:
                if 'FAIL_MODE:' in line:
                    FAIL_MODE  = line.split()[1]
                    # convert FAIL_MODE string to fail mode index
                    itmp = 0 ;  IFAIL = -9
                    for string in FAIL_MODE_STRINGS:
                        if FAIL_MODE == string:
                            IFAIL = itmp
                        itmp += 1
                    if IFAIL >= 0 :
                        n_job_fail_list[0]     += 1
                        n_job_fail_list[IFAIL] += 1

        return n_job_fail_list

        # end read_failure_stats

    def read_command_lines(self):

        # read all CPU*CMD files and store lines, excluding
        # lines that are blank or comments. This allows parsing
        # memory for commands to reproduce failures.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR'] 

        command_file_list  = glob.glob1(script_dir,CMD_wildcard)
        command_lines      = []
        for cmd_file in command_file_list :
            CMD_FILE = (f"{script_dir}/{cmd_file}")
            with open(CMD_FILE,"r") as f_cmd :
                for line in f_cmd :
                    if len(line) > 0 and line[0] != '#' :
                        command_lines.append(line.rstrip("\n"))

        nlines = len(command_lines)
        msg = (f"  Stored {nlines} command lines from " \
               f"{CMD_wildcard} files")
        #logging.info(msg)

        return command_lines
        # end read_command_lines

    def failure_summary(self):
        # pre-pend summary info at top of failure-log file,
        # and also at top of MERGE.LOG file.
        # Use linux sed tricks to add lines at top of existing file.
        # If there are no failures, MERGE.LOG is still updated
        # but non-existing FAILURES.LOG is ignored.

        #submit_info_yaml = self.config_prep['submit_info_yaml']

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR'] 
        FAIL_SUMMARY_PATHFILE = (f"{script_dir}/{FAIL_SUMMARY_FILE}")
        MERGE_LOG_PATHFILE    = (f"{output_dir}/{MERGE_LOG_FILE}")

        n_list = len(FAIL_MODE_STRINGS) # TOT, ABORT, ZERO_EVT, UNKNOWN

        if os.path.isfile(FAIL_SUMMARY_PATHFILE):
            n_job_fail_list = self.read_failure_stats()
            found_failures  = True 
        else:
            # no failures, so load all zeros for n_fail per mode
            n_job_fail_list = [ 0 ] * n_list
            found_failures  = False 
            n_list = 1  # print only total failures

        # prepare string for sed command
        string_insert = "FAIL_SUMMARY:\\\n"
        for i in range(0,n_list) :
            mode    = FAIL_MODE_STRINGS[i]
            comment = FAIL_MODE_COMMENTS[i]
            nfail   = n_job_fail_list[i]
            key     = (f"NJOB_FAIL({mode}):")
            string  = (f"   {key:<25}  {nfail:2d}  # {comment}")
            string_insert += (f"{string}\\\n")  # for sed command
            
        # construct sed command to pre-pend summary at top of file
        # ?? can python do this ??
        cmd_sed = (f"sed -i '1i{string_insert}'")

        if found_failures :
            CMD = (f"{cmd_sed} {FAIL_SUMMARY_PATHFILE}")
            os.system(CMD)

            # if script_dir is not output_dir, copy FAILURES.LOG to
            # output_dir so that it's easier to find
            if script_dir != output_dir :
                shutil.copy(FAIL_SUMMARY_PATHFILE,output_dir)

        # always update merge log file
        CMD = (f"{cmd_sed} {MERGE_LOG_PATHFILE}")
        os.system(CMD)

        # return total number of failures
        return n_job_fail_list[0]

        # end failure_summary

    def get_job_stats(self, script_dir, log_file_list, yaml_file_list, 
                      yaml_key_list):
        
        # Called when all DONE files exist.
        # Loop over input log_file_list and for each file,
        #    + check if yaml file exists
        #    + if yaml exists, read & extract stats for each
        #       yaml_key_list
        #
        # Store both stat sums and a stat-list overy YAML files.
        #
        # Only yaml files are parsed here; log_file is passed to
        # check_for_failures so that ABORT message can be
        # extracted elsewhere .
        #
        
        n_key              = len(yaml_key_list)
        n_split            = len(log_file_list)
        key_AIZ            = 'ABORT_IF_ZERO'
        aiz_list           = [ -9 ] * n_split
        aiz_max            = 0
        n_aiz_zero         = 0
        
        # init output dictionary
        job_stats = { 'nfail' : 0 }
        for key in yaml_key_list :
            key_sum  = (f"{key}_sum")   # store sum
            key_list = (f"{key}_list")  # store stats for each split job
            job_stats[key_list] = [ 0 ] * n_split
            job_stats[key_sum]  =   0
    
        for isplit in range(0,n_split):
            log_file  = log_file_list[isplit]
            yaml_file = yaml_file_list[isplit]
            LOG_FILE  = (f"{script_dir}/{log_file}")
            YAML_FILE = (f"{script_dir}/{yaml_file}")

            if os.path.isfile(YAML_FILE) :
                stats_yaml       = util.extract_yaml(YAML_FILE)
            
                aiz              = stats_yaml[key_AIZ]
                aiz_list[isplit] = aiz
                if aiz > aiz_max : aiz_max = aiz
                if aiz == 0 : n_aiz_zero += 1
            
                for item in yaml_key_list :
                    key, key_sum, key_list  = self.keynames_for_job_stats(item)
                    job_stats[key_list][isplit] = stats_yaml[key]
                    job_stats[key_sum]         += stats_yaml[key]
                
                    # fix format for CPU
                    if 'CPU' in key :
                        cpu_sum = (f"{job_stats[key_sum]:.1f}")
                        job_stats[key_sum] = float(cpu_sum)

        # - - - - - - - - - - - - - - - - - - - - a
        # Examine failures. If there is no output YAML file, aiz=-9 flags
        # a clear failure -> abort.
        # If any aiz == 0, it's tricky because very low-stat jobs (e.g., KN)
        # can result in aiz=0 due to random Poisson fluctuations. Goal is to
        # abort only if aiz=0 cannot be due to low-stat random fluctuation.
        #
        # If any aiz > aiz_thresh, then abort if any aiz=0; otherwise abort 
        # only if all aiz=0. The logic is to intervene only if all 
        # aiz < aiz_thresh and a subset of aiz=0; in this case, aiz += 1
        # to suppress failure from aiz=0. Note that aiz = -9 -> -8,
        # and still fails since YAML output wasn't found.
        #

        aiz_thresh  = 30  # P_FF ~ E-6; submit_batch_jobs.sh -H AIZ 

        subset_aiz_zero = n_aiz_zero>0 and n_aiz_zero < n_split # some 0, not all
        if aiz_max < aiz_thresh and subset_aiz_zero :
            for i in range(0,n_split): aiz_list[i] += 1
            msg = (f"\t max({key_AIZ})={aiz_max} -> suppress abort from AIZ=0.")
            logging.info(msg)

        # loop again over split jobs and check for failures.
        for isplit in range(0,n_split):
            log_file   = log_file_list[isplit]
            aiz        = aiz_list[isplit]
            found_fail =  self.check_for_failure(log_file, aiz, isplit+1)
            if found_fail : job_stats['nfail'] += 1

        return job_stats

        # end get_job_stats

    def keynames_for_job_stats(self,keyname_base):
        keyname_sum  = (f"{keyname_base}_sum")
        keyname_list = (f"{keyname_base}_list")
        return keyname_base, keyname_sum, keyname_list
        # end keynames_for_job_stats

    def log_assert(self,condition, msgerr):
        # same as log_assert in util, except here it also
        # + writes FAIL to done_stamp_file
        # + appends msgerr to MERGE.LOG

        if condition is False :
            output_dir      = self.config_prep['output_dir']
            done_stamp_list = self.config_prep['done_stamp_list']
            MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")

            util.write_done_stamp(output_dir, done_stamp_list,STRING_FAIL)

            if os.path.isfile(MERGE_LOG_PATHFILE) :
                with open(MERGE_LOG_PATHFILE,"a") as f:
                    f.write(f"# {SNANA_ABORT_STRING}\n")
                    for msg in msgerr :
                        f.write(f"#   {msg}\n")
                    f.write(f"#\n")

            util.log_assert(condition, msgerr)        

# ======= END OF FILE =========

