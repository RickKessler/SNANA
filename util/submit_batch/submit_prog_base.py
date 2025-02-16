#
# Created July 2020 by R.Kessler & S. Hinton
#
# TO DO (JULY 2022); for single node option, remove mem-per-cpu 
#
# Base class Program
#
#     HISTORY
# Jan 2 2021: add small delay in each CPU* file to avoid jobs finishing
#             before all are submitted, resulting in pid-check failure.
#
# Jan 6 2021: 
#   in write_batch_file, add REPLACE_[WALLTIME,NTASK,CPU_PER_TASK].
#   These keys are hard-coded for submit_batch, but can be altered
#   by Pippin for non-SNANA jobs.
#
# Jan 21 2021: 
#   write_command_file return n_job_cpu; if n_job_cpu==0, sleep 10 sec
#   so that npid check doesn't fail.
#
# Apr 23 2021: fix --merge_reset and -M to hopefully allow interactive 
#              merge recover (-H FIT for instructions)
#
# May 24 2021: new function submit_iter2()
# Aug 09 2021: implement optional --snana_dir arg
# Oct 29 2021: add optional TABLE_EXTRA to MERGE.LOG file
# Oct 30 2021: new ENV_REQUIRE key for any task.
# Dec 08 2021: write wall time and CPU sum in hours (no option for minutes)
# Apr 12 2022: check for docker image in sbatch-file; if there, use shifter
# Jul 18 2022: check new CONFIG option "BATCH_SINGLE_NODE: True"
# Nov 15 2022: call merge_driver_exit() to write merge-process time
# Dec 02 2022: write SNANA_VERSION to SUBMIT.INFO file
# Jul 14 2023: fix pid-check logic to avoid jobs finishing before all
#               jobs are launched. See update_slurm_pid_list().
#
# Jan 23 2025: fix nomerge_last() to work properly when OUTDIR is not under CWD.
# Feb 06 2025: check for BAD_OUTPUT flag in output YAML files.
#              Initial use is to check for crazy errors in BBC fit.
#
# Feb 16 2025: fix kill-job logic to work properly when first iteration BBC job fails,
#              so that both iterations are stopped and produce STOP in ALL.DONE file.
#
# ============================================

#import argparse
import os, sys, shutil, yaml
import logging
#import coloredlogs
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
        args        = config_yaml['args']
        msgerr      = []
            
        CONFIG = config_yaml['CONFIG']
        if 'JOBNAME' in CONFIG :
            config_prep['program'] = CONFIG['JOBNAME']

        config_prep['snana_version'] = util.get_snana_version()
        
        # default is 1 submit, so iter is None
        # (bbc might change this to 1 and 2) 
        config_prep['submit_iter'] = None  # May 24 2021

        if args.merge_flag :  # bail for merge process
            return

        self.check_env_required(config_yaml)

        # - - - - -
        program = config_prep['program']
        logging.info(f" Submit host;  {HOSTNAME}")
        logging.info(f" Program name: {program}")

        if program == PROGRAM_NAME_UNKNOWN :
            input_file  = args.input_file 
            msgerr.append(f"Unknown program name.")
            msgerr.append(f"Must define program with JOBNAME key")
            msgerr.append(f"in {input_file}")
            util.log_assert(False,msgerr)

        # Oct 12 2021: option to run interactive jobs and check if job aborted
        if args.check_abort :
            self.prep_check_abort(config_yaml)

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
    def merge_config_prep(self,output_dir):
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

    @abstractmethod
    def check_abort(self):
        print(f"\n WARNING: not implemented.")

    def prep_check_abort(self,config_yaml):
        # Prepare to run interactive job with 1 event and check for abort.
        # set ncore=1

        CONFIG = config_yaml['CONFIG']

        if 'NODELIST' in CONFIG  :
            NODELIST    = CONFIG['NODELIST']
            FIRSTNODE   = NODELIST.split()[0]
            CONFIG['NODELIST'] = FIRSTNODE
        else :
            config_yaml['args'].ncore = 1

        # disable cleanup flag
        CONFIG['CLEANUP_FLAG'] = 0

        # end prep_check_abort

    def check_env_required(self, config_yaml):

        CONFIG = config_yaml['CONFIG']
        msgerr = []

        # check conda env for SALT3 (Feb 2021)
        ENV = 'CONDA_DEFAULT_ENV'
        if ENV in CONFIG :
            ENV_value_expect = CONFIG[ENV]
            ENV_value        = os.getenv(ENV)
            if ( ENV_value != ENV_value_expect ) :
                msgerr.append(f"Expected  ${ENV} = {ENV_value_expect} ; ")
                msgerr.append(f"but found ${ENV} = {ENV_value} ")
                util.log_assert(False,msgerr)

        # check optional ENV required to exist (doesn't matter what value)
        key = CONFIG_KEYNAME_ENV_REQUIRE
        if key in CONFIG:
            ENV_name_list  = CONFIG[key].split()
            missing_ENV_list = []
            for ENV_name in ENV_name_list:
                ENV_value = os.getenv(ENV_name,None)
                if ENV_value is None:
                    missing_ENV_list.append(ENV_name)
                    logging.info(f"  ERROR: missing required env ${ENV_name}")
                else:
                    logging.info(f"  Found required ${ENV_name} = {ENV_value}")

            n_missing = len(missing_ENV_list)
            if n_missing > 0:
                msgerr.append(f"This/these {n_missing} required envs " \
                              f"are not set:")
                msgerr.append(f"    {missing_ENV_list}")
                msgerr.append(f"Check {key} key in config input file.")
                util.log_assert(False,msgerr)
            else:
                print('')

        sys.stdout.flush()
        # end check_env_required

    def parse_batch_info(self,config_yaml,config_prep):
    
        # check of SSH or BATCH, and parse relevant strings


        NODELIST      = ''
        n_core        = 0
        submit_mode   = "NULL"
        node_list     = []

        memory        = BATCH_MEM_DEFAULT
        maxjob        = BATCH_MAXJOB_DEFAULT
        walltime      = BATCH_WALLTIME_DEFAULT
        nthreads      = BATCH_NTHREADS_DEFAULT
        batch_single_node = False

        kill_flag     = config_yaml['args'].kill
        n_core_arg    = config_yaml['args'].ncore
        msgerr        = []
        config_prep['nodelist']       = ''
        config_prep['batch_command']  = ''
        config_prep['command_docker']        = None
        config_prep[ENV_SNANA_SETUP_COMMAND] = None

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

            if n_core_arg is None :
                n_core = int(BATCH_INFO[2]) # from CONFIG input
            else:
                n_core = n_core_arg  # command-line override

            self.check_docker_image(template)
                
            node_list   = [HOSTNAME] * n_core  # used for script file name
            submit_mode = SUBMIT_MODE_BATCH
            config_prep['batch_command']  = command
            config_prep['BATCH_TEMPLATE'] = template  # all caps -> full path

            logging.info(f"\t Batch command:    {command}" )
            logging.info(f"\t Batch template:   {template}" )
            logging.info(f"\t Batch n_core:     {n_core}" )
        else :
            msgerr.append(f"Could not find BATCH_INFO or NODELIST.")
            msgerr.append(f"Check CONFIG block in the input file.")
            util.log_assert(False, msgerr)

        # check optional memory input
        if 'BATCH_MEM' in CONFIG :
            memory = str(CONFIG['BATCH_MEM'])

        # check optional walltime (Jan 6 2021)
        if 'BATCH_WALLTIME' in CONFIG :
            walltime = CONFIG['BATCH_WALLTIME']

        if 'BATCH_NTHREADS' in CONFIG :
            nthreads = CONFIG['BATCH_NTHREADS']

        # option to run everything on one node
        if 'BATCH_SINGLE_NODE' in CONFIG :
            batch_single_node = CONFIG['BATCH_SINGLE_NODE']
            config_prep['batch_command']  += f" -n {n_core}"

        sys.stdout.flush()

        config_prep['n_core']      = n_core 
        config_prep['submit_mode'] = submit_mode
        config_prep['node_list']   = node_list
        config_prep['memory']      = memory
        config_prep['walltime']    = walltime
        config_prep['maxjob']      = maxjob
        config_prep['nthreads']    = nthreads
        config_prep['batch_single_node'] = batch_single_node

        return
        
    # end parse_batch_info

    def check_docker_image(self,sbatch_template):
        
        # Created Apr 12 2022
        # check for docker image key "--image" in sbatch-template file

        use_docker_image    = False
        ENV_HOST            = None
        command_docker      = None
        SNANA_SETUP_COMMAND = None
        msgerr           = []

        batch_lines = open(sbatch_template,'r').read()
        if "--image" in batch_lines:
            use_docker_image = True

        # Beware hard-wired hack:
        # Since docker has security issues, labs develop their own implementation
        # of docker with a different name. The "docker_command_dict" below maps 
        # ENV_HOST into the docker command to launch jobs.
        # For example NERSC machines define ENV_HOST:  $NERSC_HOST = Cori,
        # and they user shifter to launch jobs.
        # Beware that shifter -> Podman at some point on Perlmutter
        docker_command_dict = {
            #  ENV_HOST      docker-command
            'NERSC_HOST'  : 'shifter --module=none ' ,
            'SLAC_HOST'   : 'singularity'  # ?? check if we ever use SLAC cluster
        }

        if use_docker_image:
            for ENV_tmp, cmd_tmp in docker_command_dict.items():
                ENV_VAL_tmp = os.getenv(ENV_tmp,None)
                if ENV_VAL_tmp is not None:
                    ENV_HOST       = ENV_tmp
                    ENV_HOST_VAL   = ENV_VAL_tmp
                    command_docker = cmd_tmp
                
            if command_docker is None:
                msgerr.append("Cannot use docker image on this machine.")
                msgerr.append("Known docker machines are")
                msgerr.append("  {docker_command_dict}")
                util.log_assert(False,msgerr)

            # also check ENV for command to setup snana
            SNANA_SETUP_COMMAND = os.getenv(ENV_SNANA_SETUP_COMMAND,None)
            if SNANA_SETUP_COMMAND is None:
                msgerr.append(f"Using docker image requires ENV")
                msgerr.append(f"  ${ENV_SNANA_SETUP_COMMAND}")
                msgerr.append(f"but it is not set.")
                util.log_assert(False,msgerr)
                
            logging.info(f"\t Found docker image: " \
                         f"use {command_docker} for ${ENV_HOST}={ENV_HOST_VAL}")

        # - - - - -
        # load output
        self.config_prep['command_docker']        = command_docker
        self.config_prep[ENV_SNANA_SETUP_COMMAND] = SNANA_SETUP_COMMAND


        if command_docker is None :  
            sh = 'sh'
        else:
            sh = f"{command_docker} sh"
        self.config_prep['command_sh'] = sh

        return
        #end check_docker_image

    def kill_jobs(self):

        # create_output_dir sets dir names, but won't create anything
        # because kill flag is set
        self.create_output_dir()
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep.setdefault('submit_info_yaml',None)

        if submit_info_yaml is None:
            INFO_PATHFILE    = f"{output_dir}/{SUBMIT_INFO_FILE}"
            submit_info_yaml = util.extract_yaml(INFO_PATHFILE, None, None )

        # xxxxxxxx
        #print(f"\n xxx kill_jobs: submit_info_yaml = \n{submit_info_yaml}\n") # .xyz
        #print(f"\n xxx kill_jobs: config_prep = \n{self.config_prep} \n\n")
        #sys.stdout.flush() # xxxx
        # xxxxxxx

        done_list     = submit_info_yaml['DONE_STAMP_LIST']
        submit_mode   = submit_info_yaml['SUBMIT_MODE'] 
        batch_command = SBATCH_COMMAND   # fragile alert warning

        # xxxx mark delete Feb 16 2025 xxxxxxx
        #output_dir    = self.config_prep['output_dir']
        #done_list     = self.config_prep['done_stamp_list']
        #submit_mode   = self.config_prep['submit_mode'] 
        #batch_command = self.config_prep['batch_command'] 
        # xxxxxxxxxxxx

        IS_SSH        = submit_mode == SUBMIT_MODE_SSH
        IS_BATCH      = submit_mode == SUBMIT_MODE_BATCH


        # write FAIL stamps before killing the jobs
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        time_now            = datetime.datetime.now()
        with open(MERGE_LOG_PATHFILE, 'a') as f:
            f.write(f"\n !!! JOBS KILLED at {time_now} !!! \n")

        util.write_done_stamp(output_dir, done_list, STRING_STOP)

        # - - - - - - - - - 
        if IS_SSH :
            self.kill_ssh_jobs()

        elif IS_BATCH and batch_command == SBATCH_COMMAND :
            self.kill_sbatch_jobs()

        else:
            msgerr = []
            msgerr.append(f"Unable to kill jobs for:")
            msgerr.append(f"  submit_mode   = '{submit_mode} ")
            msgerr.append(f"  batch_command = '{batch_command}")
            util.log_assert(False,msgerr)


        # xxxxxxxxxxx mark delete Feb 2025 xxxxxxxxxxxx
        # if we get here, leave notice in both MERGE.LOG and ALL.DONE files
        #MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        #time_now            = datetime.datetime.now()
        #with open(MERGE_LOG_PATHFILE, 'a') as f:
        #    f.write(f"\n !!! JOBS KILLED at {time_now} !!! \n")
        #        util.write_done_stamp(output_dir, done_list, STRING_STOP)
        # xxxxxxxxxxxxxxxxxx end mark xxxxxxxxxx


        # end kill_jobs

    def kill_ssh_jobs(self):
        node_list   = self.config_prep['node_list']
        for node in node_list :
            logging.info(f" Kill jobs on {node}")
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
        INFO_PATHFILE    = f"{output_dir}/{SUBMIT_INFO_FILE}"
        submit_info_yaml = util.extract_yaml(INFO_PATHFILE,None,None)

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
            cmd_kill   = f"scancel --name={job_name}"
            logging.info(f"\t {cmd_kill}")
            os.system(cmd_kill)
            
        # check to kill job on cpunum_last
        if cpunum_last >= 0 :
            cmd_kill   = f"scancel --name={job_name_last}"
            logging.info(f"\t {cmd_kill}")
            os.system(cmd_kill)

        logging.info(f" Done cancelling {njob_kill} batch jobs.")
        logging.info(f" Check remaining {USERNAME} jobs ... ")
        time.sleep(3)
        os.system(f"squeue -u {USERNAME}")
        logging.info(f"\n Beware that batch jobs might remain from other submits.")

        # end kill_sbatch_jobs

    def write_script_driver(self):

        # For each CPU, create batch script (CPU[nnn]*.BATCH and also
        # command script (CPU[nnn].CMD) that has list of native commands.
        # For ssh mode, these .CMD files are sourced upon login, and
        # thus .BATCH files are not needed for ssh mode.
        # BATCH files are written here and do not depend on task.
        # CMD files are program-dependent via call to
        #       self.write_command_file(icpu,COMMAND_FILE)
        #
        # Jan 21 2201: write_command_file returns n_job_cpu;
        #     if n_job_cpu==0, add extra delay to avoid npid error
        #
        # Dec 04 2021: write CPU*DONE file (for merge_background)

        CONFIG      = self.config_yaml['CONFIG']
        args        = self.config_yaml['args']
        input_file  = args.input_file 
        snana_dir   = args.snana_dir
        output_dir  = self.config_prep['output_dir']
        script_dir  = self.config_prep['script_dir']
        n_core      = self.config_prep['n_core']
        submit_mode = self.config_prep['submit_mode']
        command_docker = self.config_prep['command_docker']
        batch_single_node = self.config_prep['batch_single_node']
        node_list      = self.config_prep['node_list']
        program        = self.config_prep['program']
        submit_iter    = self.config_prep['submit_iter']

        # for each cpu, store name of script and batch file
        command_file_list = []  # command file name, no path
        batch_file_list   = []   # batch file name, no patch
        COMMAND_FILE_LIST = []   # includes full path
        BATCH_FILE_LIST   = []   # idem
        cmdlog_file_list  = []
        job_name_list     = []
        n_core_with_jobs  = 0

        # loop over each core and create CPU[nnn]_JOBS.CMD that
        # are used for either batch or ssh. For batch, also 
        # create batch_file using BATCH_TEMPLATE

        logging.info(f"  Create command files:")

        for icpu in range(0,n_core) :
            node          = node_list[icpu]
            cpu_name      = f"CPU{icpu:04d}"
            prefix        = f"CPU{icpu:04d}_JOBLIST_{node}"
            command_file  = f"{prefix}.CMD"
            log_file      = f"{prefix}.LOG"
            COMMAND_FILE  = f"{script_dir}/{command_file}"
            done_file     = f"{prefix}.DONE"
            DONE_FILE     = f"{script_dir}/{done_file}"

            logging.info(f"\t Create {command_file}")

            if submit_iter is None:
                # normal task is single submit, so ignore iter in job name
                job_name      = f"{input_file}-{cpu_name}"
            else:
                # job name depends on submit_iter
                ii        = f"iter{submit_iter}"
                job_name  = f"{input_file}_{ii}-{cpu_name}"

            # compute small delay per core to avoid first jobs
            # finishing before all are submitted, then failing
            # the pid-submit check. Delay is largest for core 0, 
            # then is reduced by 0.2 sec per core. For 100 cores,
            # first delay is 20 sec.
            delay = float(n_core - icpu)/5

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

                # print actual start time for this batch job
                f.write(f"echo TIME_START: " \
                        f"`date +%Y-%m-%d` `date +%H:%M:%S` \n")

                # set ENV for global/uniform start time that propagates
                # into all sim-output READMEs. 
                f.write(f"export SNANA_TIME_START='{time_submit_start}'\n")

                f.write(f"echo 'Sleep {delay} sec " \
                        f"(wait for remaining batch-submits)' \n")
                f.write(f"sleep {delay} \n")

                if command_docker is not None:
                    SNANA_SETUP_COMMAND = \
                        self.config_prep[ENV_SNANA_SETUP_COMMAND]
                    f.write(f"echo ' ' \n")
                    f.write(f"echo 'Setup SNANA for docker using command:' \n")
                    f.write(f"echo '   {SNANA_SETUP_COMMAND}' \n")  # <= print setup command
                    f.write(f"{SNANA_SETUP_COMMAND}\n")             # <= do snana setup
                    
                f.write(f"echo ' ' \n")

                f.write(f"echo 'Begin {command_file}' \n\n")

                if snana_dir is not None:
                    path_list = f"{snana_dir}/bin:{snana_dir}/util:$PATH"
                    f.write(f"export SNANA_DIR={snana_dir} \n")
                    f.write(f"export PATH={path_list}\n" )

                f.write(f"echo SNANA_DIR = $SNANA_DIR \n")
                f.write(f"python --version \n")  # Dec 2023

                if STOP_ALL_ON_MERGE_ERROR :
                    f.write(f"set -e \n") 

                # write program-specific content
                n_job_cpu = self.write_command_file(icpu,f)

                # if there are no jobs, sleep another 5 seconds so that
                # batch job does not immediately exit and fail npid check.
                if n_job_cpu == 0 :
                    f.write(f"echo 'No jobs -> sleep 10s to " \
                            f"ensure npid check is ok' \n")
                    f.write(f"sleep 10 \n")
                else:
                    n_core_with_jobs += 1

                # create CPU*DONE file to make clear that all tasks
                # for this core have run
                f.write(f"\ntouch {DONE_FILE}\n") 

            # - - - - - 
            # write extra batch file for batch mode
            if ( submit_mode == SUBMIT_MODE_BATCH ):
                batch_file = f"{prefix}.BATCH"
                BATCH_FILE = f"{script_dir}/{batch_file}"
                new_batch_file = True
                if batch_single_node and icpu > 0: new_batch_file = False

                if new_batch_file :
                    batch_file_list.append(batch_file)
                    BATCH_FILE_LIST.append(BATCH_FILE)
                    self.write_batch_file(batch_file, log_file,
                                          command_file, job_name)
                if batch_single_node:
                    batch_file0 = batch_file_list[0]
                    last_job    = (icpu == n_core-1)
                    self.append_batch_file(batch_file0, command_file, last_job)

        # store few thigs for later
        self.config_prep['cmdlog_file_list']  = cmdlog_file_list
        self.config_prep['command_file_list'] = command_file_list
        self.config_prep['COMMAND_FILE_LIST'] = COMMAND_FILE_LIST
        self.config_prep['job_name_list']     = job_name_list
        self.config_prep['batch_file_list']   = batch_file_list
        self.config_prep['BATCH_FILE_LIST']   = BATCH_FILE_LIST
        self.config_prep['n_core_with_jobs']  = n_core_with_jobs

        # check option to run merge process in the background
        if args.merge_background:
            self.write_script_merge_background()

        # make all CMD files group-executable with +x.
        # Python os.chmod is tricky because it may only apply 
        # permission to user, or wipe out group privs. 
        # To make things easier here, we just use os.system().
        cmd_chmod = f"chmod +x {script_dir}/CPU*.CMD"
        os.system(cmd_chmod);

        n_job_tot   = self.config_prep['n_job_tot']
        logging.info(f" BATCH DRIVER JOB COUNT SUMMARY: " \
                     f"{n_job_tot} total jobs distributed on {n_core} cores")

        
        # check option to force crash (to test higher level pipelines)
        if args.force_crash_prep :
            printf(" xxx force batch-prep crash with C-like printf xxx \n")

        # end write_script_driver

    def write_script_merge_background(self):

        # Created Dec 4 2021 by R.Kessler
        # Invoked by --merge_background option.
        # Write bash script to execute merge process until all of the
        # CPU*DONE files exist.
        # This bash script is launched in the background after the
        # nominal tasks have been submitted.

        n_core         = self.config_prep['n_core']
        script_dir     = self.config_prep['script_dir']
        args           = self.config_yaml['args']
        input_file     = args.input_file 
        t_stamp        = seconds_since_midnight
        cpunum         = 0

        if args.prescale > 1 :
            t_sleep = 20  # sleep time between checking merge
        else:
            t_sleep = 100

        base_name         = "CPU_MERGE_BACKGROUND"
        cpu_merge_script  =  f"{script_dir}/{base_name}.CMD"
        cpu_merge_log     =  f"{script_dir}/{base_name}.LOG"
        
        merge_script     = sys.argv[0]  # $path/submit_batch_jobs

        arg_cpu = f"--cpunum {cpunum}"
        arg_t   = f"-t {t_stamp}"
        arg_m   = f"--merge"
        arg_M   = f"--MERGE_LAST"

        merge_args       = f"{input_file} {arg_m} {arg_t} {arg_cpu}"
        merge_args_final = f"{input_file} {arg_M} {arg_t} {arg_cpu}"
        wildcard = "CPU*DONE"
        wildcard_echo = "CPU\*DONE"

        # store for when it's launched
        self.config_prep['cpu_merge_script'] = cpu_merge_script
        self.config_prep['cpu_merge_log']    = cpu_merge_log

        with open(cpu_merge_script,"wt") as f:
            f.write("#!/usr/bin/env bash \n")
            f.write(f"echo HOST = {HOSTNAME} \n")
            f.write(f"echo \n")
            f.write(f"n_done=0\n\n")

            f.write(f"while [ $n_done -lt {n_core} ] \n")
            f.write(f"do\n")
            f.write(f"  sleep {t_sleep} \n")
            f.write(f"  echo '# ======== MERGE_BACKGROUND CHECK ==========' \n")

            f.write(f"  echo Found $n_done of {n_core} {wildcard_echo} files.\n")
            f.write(f"  echo Run merge at "
                    f"`date +%Y-%m-%d` `date +%H:%M:%S` \n")

            f.write(f"  cd {CWD}\n")
            f.write(f"  {merge_script} {merge_args} \n")
            f.write(f"  \n")
            f.write(f"  cd {script_dir}\n")
            f.write(f"  n_done=`ls {wildcard} 2> /dev/null | wc -l` \n")
            f.write(f"  echo \n")
            f.write(f"done\n")

            f.write("\n")
            f.write(f"echo Found all $n_done of {n_core} {wildcard_echo} files.\n")
            f.write(f"echo Run final merge   \n")
            f.write(f"  cd {CWD}\n")
            f.write(f"  {merge_script} {merge_args_final} \n")
            f.write(f"echo Done with merge_background.\n")
        return
        # end write_script_merge_background

    def write_batch_file(self, batch_file, log_file, command_file, job_name):

        # Create batch_file that executes "source command_file"
        # BATCH_TEMPLATE file is read, and lines are modified using
        # string replacements for
        #   REPLACE_NAME, REPLACE_LOGFILE, REPLACE_MEM, REPLACE_JOB
        # with job-specific info.
        # Note that lower-case xxx_file has no path; 
        # upper case XXX_FILE includes full path
        #
        # Apr 12 2022: check for docker command (e.g., 'shifter')
        
        BATCH_TEMPLATE   = self.config_prep['BATCH_TEMPLATE'] 
        script_dir       = self.config_prep['script_dir']
        command_docker   = self.config_prep['command_docker']
        command_sh       = self.config_prep['command_sh']
        replace_memory   = self.config_prep['memory']
        replace_walltime = self.config_prep['walltime']
        replace_cpus_per_task = self.config_prep['nthreads'] # 08/apr/2022
        batch_single_node = self.config_prep['batch_single_node']

        BATCH_FILE      = f"{script_dir}/{batch_file}"

        # get strings to replace 
        replace_job_name   = job_name

        # nothing to change for log file
        replace_log_file   = log_file  

        use_docker = (command_docker is not None)
        
        # xxxx replace_job_cmd=f"cd {script_dir} \n{command_sh} {command_file}"
        replace_job_cmd  = f"cd {script_dir}"
        if not batch_single_node:
            replace_job_cmd += f"\n{command_sh} {command_file}"            

        # Jan 6 2021: add few more replace keys that can be modified
        # by non-SNANA tasks (e.g., classifiers, CosmoMC ...)
        replace_ntask        = 1
        #replace_cpu_per_task = 1

        # - - - define list of strings to replace - - - - 

        REPLACE_KEY_DICT = { 
            'REPLACE_NAME'          : replace_job_name,
            'REPLACE_MEM'           : replace_memory,
            'REPLACE_LOGFILE'       : replace_log_file,
            'REPLACE_JOB'           : replace_job_cmd,
            'REPLACE_WALLTIME'      : replace_walltime,
            'REPLACE_NTASK'         : replace_ntask,
            'REPLACE_CPUS_PER_TASK'  : replace_cpus_per_task,
        }

        if use_docker: 
            ENV_name = ENV_SNANA_IMAGE_DOCKER
            replace_image_docker = os.getenv(ENV_name,None)
            if replace_image_docker is None:
                msgerr = []
                msgerr.append(f"Missing required ${ENV_name}")
                msgerr.append(f"for --image argument in batch file.")
                self.log_assert(False, msgerr)
            REPLACE_KEY_DICT['REPLACE_IMAGE_DOCKER'] = replace_image_docker

        # - - - - 
        batch_line_list = open(BATCH_TEMPLATE,'r').readlines()
        b = open(BATCH_FILE,"w")
        for line in batch_line_list:
            if batch_single_node and 'REPLACE_MEM' in line:
                continue
            for KEY,VALUE in REPLACE_KEY_DICT.items():
                if KEY in line:
                    line = line.replace(KEY,str(VALUE))
            b.write(f"{line}")  # line includes \n
        b.close()

        return
        # end write_batch_file

    def append_batch_file(self, batch_file, command_file, last_job):
        script_dir       = self.config_prep['script_dir']
        command_sh       = self.config_prep['command_sh']
        BATCH_FILE       = f"{script_dir}/{batch_file}"
        ampsand          = '&'
        # xxx put back after testing if last_job:  ampsand = '' 

        with open(BATCH_FILE,"at") as f :
            f.write(f"{command_sh} {command_file} {ampsand} \n")
            
            # on last job, wait for done file to exit script,
            # otherwise remaining batch jobs/merge is killed.
            if last_job :
                output_dir  = self.config_prep['output_dir']
                f.write(f"\n# \n") 
                f.write(f"echo ' ' \n")
                f.write(f"echo Wait for {DEFAULT_DONE_FILE} "
                        f"file before exiting. \n")
                done_file = f"{output_dir}/{DEFAULT_DONE_FILE}"
                f.write(f"while [ ! -f {done_file} ] ; " \
                        f"do sleep 60; done\n" )
        return
        #end append_batch_file

    def prep_JOB_INFO_merge(self,icpu,ijob,merge_force):
        # Return JOB_INFO dictionary of strings to run merge process.
        # Inputs:
        #   icpu = 0 to n_core-1
        #   ijob = 1 to n_job_tot
        #   merge_force = True => pass --merge_force instead of default --merge
        #
        # Merge task must is the form
        #   python <thisScript> <inputFile> arg_list
        # arg_list includes
        #  --merge      -> merge and quit or
        #  --MERGE_LAST -> wait for all DONE files to appear, 
        #                   then merge it all.
        #  -t <Nsec>   time stamp to verify merge and submit jobs
        #  --cpunum <cpunum>  in case specific CPU needs to be identified    
        #
        #  May 24 2021: check outdir override from command line
        #  Apr 08 2022: pass merge_force arg

        args                = self.config_yaml['args']
        input_file          = args.input_file
        nomerge             = args.nomerge
        output_dir_override = args.outdir 
        devel_flag          = args.devel_flag
        check_abort         = args.check_abort

        n_core         = self.config_prep['n_core']
        n_job_tot      = self.config_prep['n_job_tot']
        Nsec           = seconds_since_midnight

        JOB_INFO = {}

        # determine if this is last job for this cpu
        last_job_cpu  = (n_job_tot - ijob) < n_core

        # if check_abort, skip merge except for last job
        # xxx skip_merge = check_abort and not last_job_cpu

        # - - - - 
 
        # check where to flag last merge process with -M
        if NCPU_MERGE_DISTRIBUTE == 0 :
            # cpu=0 has last merge process 
            last_merge =  last_job_cpu and icpu == 0
        else :
            # last merge is after last job, regardless of cpunum (default)
            last_merge =  ijob == n_job_tot

        if nomerge and not last_merge :
            JOB_INFO['merge_input_file'] = ""
            JOB_INFO['merge_arg_list']   = ""
            return JOB_INFO

        m_arg = "--merge"
        if merge_force:  m_arg = "--merge_force"
        if last_merge :  m_arg = "--MERGE_LAST" 

        arg_list = f"{m_arg} -t {Nsec} --cpunum {icpu}"
        
        # check for outdir override (May 24 2021)
        if output_dir_override is not None:
            arg_list += f"  --outdir {output_dir_override}"

        # check for check_abort (Oct 12 2021)
        if check_abort :
            arg_list += f" --{arg_check_abort}"

        # check for devel flag
        if devel_flag != 0 :
            arg_list += f" --devel_flag {devel_flag}"

        # check for nomerge (Mar 28 2022)
        if nomerge :
            arg_list += f" --nomerge"

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
        #
        # Write FAST if set.

        args             = self.config_yaml['args']
        prescale         = self.config_yaml['args'].prescale
        CONFIG           = self.config_yaml['CONFIG']
        n_job_tot        = self.config_prep['n_job_tot']
        n_done_tot       = self.config_prep['n_done_tot']
        n_job_split      = self.config_prep['n_job_split']
        n_core           = self.config_prep['n_core']
        n_core_with_jobs = self.config_prep['n_core_with_jobs']
        output_dir       = self.config_prep['output_dir']
        script_dir       = self.config_prep['script_dir']
        done_stamp_list  = self.config_prep['done_stamp_list']
        submit_mode      = self.config_prep['submit_mode']
        Nsec             = seconds_since_midnight
        time_now         = time_submit_start

        cleanup_flag = 1     # default
        if 'CLEANUP_FLAG' in CONFIG :
            cleanup_flag = CONFIG['CLEANUP_FLAG']  # override default

        logging.info(f"  Create {SUBMIT_INFO_FILE}")
        INFO_PATHFILE  = f"{output_dir}/{SUBMIT_INFO_FILE}"
        f = open(INFO_PATHFILE, 'w') 

        # required info for all tasks
        f.write("\n# Required info \n")

        f.write(f"CWD:         {CWD}\n")

        arg_string = " ".join(sys.argv[1:])
        f.write(f"ARG_LIST:    {arg_string}\n")

        comment = "submit time; Nsec since midnight"
        f.write(f"TIME_STAMP_NSEC:   {Nsec}    # {comment}\n")
        f.write(f"TIME_STAMP_SUBMIT: {time_now}    \n")
        f.write(f"SUBMIT_MODE:       {submit_mode} \n")

        comment = "submit from this HOST"
        f.write(f"SUBMIT_HOST:       {HOSTNAME}   # {comment} \n")

        if args.merge_background :
            merge_mode = MERGE_MODE_BACKGROUND
        elif args.nomerge : 
            merge_mode = MERGE_MODE_SKIP
        else:
            merge_mode = MERGE_MODE_DEFAULT            
        f.write(f"MERGE_MODE:       {merge_mode} \n")
        
        f.write(f"SCRIPT_DIR:       {script_dir} \n")
        f.write(f"DONE_STAMP_LIST:  {done_stamp_list} \n")
        f.write(f"CLEANUP_FLAG:     {cleanup_flag}   # flag to compress run-job dir\n")
        f.write(f"KILL_ON_FAIL:     {args.kill_on_fail} \n")

        comment = "total number of jobs (excludes sym links)"
        f.write(f"N_JOB_TOT:        {n_job_tot}     # {comment}\n")

        comment = "total number of done files (includes sym links)"
        f.write(f"N_DONE_TOT:       {n_done_tot}     # {comment}\n")

        comment = "njob to merge per task after processing"
        f.write(f"N_JOB_SPLIT:      {n_job_split}     # {comment}\n")

        comment = "number of cores"
        f.write(f"N_CORE:           {n_core}     # {comment} \n")

        comment = "n_core with jobs"
        f.write(f"N_CORE_WITH_JOBS: {n_core_with_jobs}   # {comment} \n")

        if prescale > 1 :
            f.write(f"FAST:             {prescale}    " \
                    f"# process 1/{prescale} of request\n")

        force_crash_prep  = self.config_yaml['args'].force_crash_prep
        force_crash_merge = self.config_yaml['args'].force_crash_merge
        force_abort_merge = self.config_yaml['args'].force_abort_merge
        f.write(f"FORCE_CRASH_PREP:    {force_crash_prep} \n")
        f.write(f"FORCE_CRASH_MERGE:   {force_crash_merge}\n")
        f.write(f"FORCE_ABORT_MERGE:   {force_abort_merge}\n")

        snana_version = util.get_snana_version()
        f.write(f"SNANA_VERSION:    {snana_version}\n")
  
        # append program-specific information
        f.write("\n")
        self.append_info_file(f)

        f.close()

        # end create_info_file

    def create_output_dir(self):

        # if there is no slash in output_dir name, tack on cwd
        # to make sure full path is included. Then create output_dir.
        # if output_dir already exists, clobber it.

        kill_flag       = self.config_yaml['args'].kill

        output_dir_temp,script_subdir = self.set_output_dir_name()  

        # check command line option to override outdir (May 2021)
        output_dir_temp = self.override_output_dir_name(output_dir_temp)

        if '/' not in output_dir_temp :
            output_dir = f"{CWD}/{output_dir_temp}"
        else:
            output_dir = output_dir_temp

        # define script_dir based on script_subdir
        if ( len(script_subdir) > 2 ):
            script_dir = f"{output_dir}/{script_subdir}"
        else:
            script_dir = output_dir 

        # store dir names 
        self.config_prep['output_dir'] = output_dir
        self.config_prep['script_dir'] = script_dir

        # fetch & store done stamp file(s)
        CONFIG          = self.config_yaml['CONFIG']
        done_stamp_list = [ DEFAULT_DONE_FILE ]  # always write default ALL.DON
        done_stamp_file,Found = util.parse_done_stamp(output_dir,CONFIG)
        if Found : done_stamp_list.append(done_stamp_file)  # check for another
        self.config_prep['done_stamp_list'] = done_stamp_list

        # - - - - - - - - - 
        if kill_flag : return
        # - - - - - - - - - 

        logging.info(f" Create output dir:\n   {output_dir}")
        if  os.path.exists(output_dir) : shutil.rmtree(output_dir)
        os.mkdir(output_dir)

        # next create subdir for scripts, unless subDir is ./
        if script_dir != output_dir :
            os.mkdir(script_dir)
            
        # end create_output_dir

    def get_output_dir_name(self):
        return  self.config_prep['output_dir']
        # end

    def override_output_dir_name(self,output_dir):
    
        output_dir_override = self.config_yaml['args'].outdir

        if output_dir_override is None: 
            return output_dir
        else:        
            return output_dir_override

        # end override_output_dir_name

    def create_merge_file(self):
        output_dir      = self.config_prep['output_dir'] 
        logging.info(f"  Create {MERGE_LOG_FILE}")
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        f = open(MERGE_LOG_PATHFILE, 'w') 

        # create/write the actual tabls(s) ; this is program specific
        self.create_merge_table(f) 

        f.close()

        if KEEP_EVERY_MERGELOG :
            util.backup_merge_file(MERGE_LOG_PATHFILE)
        # end create_merge_file

    def launch_jobs(self):
        
        args = self.config_yaml['args']
        check_abort = args.check_abort

        # submit all the jobs; either batch or ssh
        submit_mode    = self.config_prep['submit_mode']
        script_dir     = self.config_prep['script_dir']
        cddir          = f"cd {script_dir}"
        
        if check_abort :
            command_file_list = self.config_prep['command_file_list']
            command_file      = './' + command_file_list[0]
            ret = subprocess.run( [ command_file ], 
                                  cwd=script_dir,
                                  capture_output=False, text=True )
            sys.exit("\n Done with abort check.")

        
        elif submit_mode == SUBMIT_MODE_BATCH :
            batch_command    = self.config_prep['batch_command'] 
            batch_file_list  = self.config_prep['batch_file_list']
            nslurm_submit    = 0
            slurm_pid_list = [] ; slurm_job_name_list = []
            nbatch  = len(batch_file_list)

            for batch_file in batch_file_list :
                batch_command_list = batch_command.split()
                batch_command_list.append(batch_file)
                logging.info(f"\t Launch {batch_file} ... ")
                #ret = subprocess.run( [ batch_command, batch_file], 
                ret = subprocess.run( batch_command_list, 
                                      cwd=script_dir,
                                      capture_output=True, text=True )

                # update slurm pid list every 10 submits in case some
                # tasks finish before all jobs are launched.
                nslurm_submit += 1
                if (nslurm_submit % 10 == 0 or nslurm_submit == nbatch ) :
                    slurm_pid_list, slurm_job_name_list = \
                        self.update_slurm_pid_list(slurm_pid_list, 
                                                   slurm_job_name_list)
                    len_pid = len(slurm_pid_list)
                    logging.info(f"\t   [ len(slurm_pid_list) = {len_pid} ]")

            # - - - -
            self.examine_slurm_pid_list(slurm_pid_list, slurm_job_name_list)
            #self.fetch_slurm_pid_list_obsolete()

        elif submit_mode == SUBMIT_MODE_SSH :
            # SSH
            n_core            = self.config_prep['n_core']
            node_list         = self.config_prep['node_list']
            command_file_list = self.config_prep['command_file_list']
            cmdlog_file_list  = self.config_prep['cmdlog_file_list'] 
            CONFIG            = self.config_yaml['CONFIG']
            
            login_setup = ''
            for key in CONFIG_KEYLIST_SNANA_LOGIN_SETUP :
                if key in CONFIG:
                    login_setup = f"{CONFIG[key]}"

            logging.info(f"  login_setup (for ssh):  {login_setup} ")
            qq = '"'
            for inode in range(0,n_core):
                node       = node_list[inode]
                log_file   = cmdlog_file_list[inode]
                cmd_file   = command_file_list[inode]
                cmd_source = f"{login_setup} 2> /dev/null; {cddir} ; " \
                             f"sh {cmd_file} >& {log_file} &"

                #ret = subprocess.call(["ssh", node, cmd_source] )
                logging.info(f"  Submit jobs via ssh -x {node}")
                ret = subprocess.Popen(["ssh", "-x", node, cmd_source ],
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE)

        # check to launch background merge process (Dec 2021)
        if args.merge_background :
            self.launch_merge_background()

        return
        # end launch_jobs

    def launch_merge_background(self):

        # Created Dec 4 2021 by R.Kessler
        # Launch merge-background script into background on login node.

        script_dir       = self.config_prep['script_dir']
        cpu_merge_script = self.config_prep['cpu_merge_script']
        cpu_merge_log    = self.config_prep['cpu_merge_log']
        
        script_base_name = os.path.basename(cpu_merge_script)
        log_base_name    = os.path.basename(cpu_merge_log)
        logging.info(f"\n\t Launch {script_base_name}\n")

        command = f"sh ./{script_base_name} >& {log_base_name} & "
        ret = subprocess.run( [ command ], cwd=script_dir, 
                              shell=True, capture_output=False, text=True )

        # end launch_merge_background

    def submit_iter2(self):

        # if this is the first of two submit jobs, then launch 
        # 2nd (iter) submit. Note that typical submit has only
        # one iteration in which submit_iter = None and this
        # function does nothing. Action here is only when submit_iter==1.

        args        = self.config_yaml['args']
        submit_iter = self.config_prep['submit_iter']
        if submit_iter != 1 : return

        line = "# ============================================="
        logging.info(f"")
        logging.info(f"{line}")
        logging.info(f"{line}")
        logging.info(f"{line}")

        for t in range(3,0,-1):
            logging.info(f"\t Will auto-submit 2nd iteration in {t} seconds ...")
            time.sleep(1)


        arg_list   = sys.argv 
        if not args.kill_on_fail: 
            arg_list.append('--kill_on_fail')  # Jun 27 2022

        arg_list.append('--iter2')

        if args.refac_file_check:
            arg_list.append('--refac_file_check')  # temporary, Feb 14 2025

        arg_string = " ".join(arg_list) 
        logging.info(f"\n submit_iter2 with \n  {arg_string}\n")

        ret  = subprocess.call( arg_list )

        # end launch_jobs_iter2

    def update_slurm_pid_list(self, slurm_pid_list, slurm_job_name_list):

        # for sbatch, fetch process id for each CPU; otherwise do nothing.

        batch_command     = self.config_prep['batch_command']
        script_dir        = self.config_prep['script_dir']

        msgerr = []
        if batch_command != SBATCH_COMMAND : return

        # prep squeue command with format: i=pid, j=jobname            
        cmd = f"squeue -u {USERNAME} -h -o '%i %j' "
        ret = subprocess.run( [cmd], shell=True, 
                              capture_output=True, text=True )
        pid_all = ret.stdout.split()  # list of [ pid, jobname, pid, jobname, etc ...]

        n = len(pid_all)
        for i in range(0,n,2):
            pid      = pid_all[i]
            job_name = pid_all[i+1]
            #logging.info(f" xxx found pid={pid}  jobname={jobname}")
            if job_name not in slurm_job_name_list:
                slurm_pid_list.append(pid)
                slurm_job_name_list.append(job_name)

        return slurm_pid_list, slurm_job_name_list

        # end update_slurm_pid_list

    def examine_slurm_pid_list(self, slurm_pid_list, slurm_job_name_list):

        # Created July 14 2023
        # examine pid_list and job_name_list and abort if any jobs
        # were/are not in the slurm queue.

        output_dir       = self.config_prep['output_dir']
        job_name_list    = self.config_prep['job_name_list']

        INFO_PATHFILE  = f"{output_dir}/{SUBMIT_INFO_FILE}"
        f = open(INFO_PATHFILE, 'a') 
        f.write(f"\nSBATCH_LIST:  # [CPU, PID, JOB_NAME] \n")

        npid_fail = 0 ; njob_tot = len(job_name_list)
        for job_name in job_name_list :
            if job_name in slurm_job_name_list:
                j_job    = slurm_job_name_list.index(job_name)
                pid      = slurm_pid_list[j_job]
                cpunum   = int(job_name[-4:])
                logging.info(f"\t pid = {pid} for {job_name}")
                f.write(f"  - [ {cpunum:3d}, {pid}, {job_name} ] \n")
            else:
                npid_fail += 1
                logging.info(f" ERROR: cannot find pid for job = {job_name}")
                continue

        f.close()
        
        if npid_fail > 0 :
            msgerr = []
            msgerr.append(f"{npid_fail} of {njob_tot} jobs NOT in queue.")
            msgerr.append(f"Check for sbatch problem; e.g., njob limit.")
            self.log_assert(False, msgerr)

        return

        return

        # end examine_slurm_pid_list

    def fetch_slurm_pid_list_obsolete(self):

        # for sbatch, fetch process id for each CPU; otherwise do nothing.
        # Nov 25 2020: list all pid failures before aborting.

        batch_command    = self.config_prep['batch_command']
        output_dir       = self.config_prep['output_dir']
        script_dir       = self.config_prep['script_dir']
        job_name_list    = self.config_prep['job_name_list']
        batch_single_node = self.config_prep['batch_single_node']

        msgerr = []
        if batch_command != SBATCH_COMMAND : return

        # for single-node option, there is only one pid to check
        if batch_single_node:
            job_name0     = job_name_list[0]
            job_name_list = [ job_name0 ]

        # prep squeue command with format: i=pid, j=jobname            
        cmd = f"squeue -u {USERNAME} -h -o '%i %j' "
        ret = subprocess.run( [cmd], shell=True, 
                              capture_output=True, text=True )
        pid_all = ret.stdout.split()

        INFO_PATHFILE  = f"{output_dir}/{SUBMIT_INFO_FILE}"
        f = open(INFO_PATHFILE, 'a') 
        f.write(f"\nSBATCH_LIST:  # [CPU,PID,JOB_NAME] \n")

        npid_fail = 0 ; njob_tot = len(job_name_list)
        for job_name in job_name_list :
            if job_name in pid_all:
                j_job    = pid_all.index(job_name)
                pid      = pid_all[j_job-1]
                cpunum   = int(job_name[-4:])
                logging.info(f"\t pid = {pid} for {job_name}")
                f.write(f"  - [ {cpunum:3d}, {pid}, {job_name} ] \n")
            else:
                npid_fail += 1
                logging.info(f" ERROR: cannot find pid for job = {job_name}")
                continue

        f.close()
        
        if npid_fail > 0 :
            msgerr.append(f"{npid_fail} of {njob_tot} jobs NOT in queue.")
            msgerr.append(f"Check for sbatch problem; e.g., njob limit.")
            self.log_assert(False, msgerr)

        return

        # end fetch_slurm_pid_list_obsolete

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
        #
        # Nov 15 2022: add timer information
        #
        # - - - - - -
        # collect time/date info to clearly mark updates in
        # CPU*.LOG file(s)

        t_merge_start = time.time()

        Nsec     = seconds_since_midnight  # current Nsec, not from submit info
        time_now = datetime.datetime.now()
        tstr     = time_now.strftime("%Y-%m-%d %H:%M:%S") 
        fnam     = "merge_driver"

        args             = self.config_yaml['args']
        MERGE_LAST       = args.MERGE_LAST
        nomerge          = args.nomerge
        merge_background = args.merge_background
        cpunum           = args.cpunum[0]
        check_abort      = args.check_abort 
        verbose_flag     = not check_abort

        if verbose_flag :
            logging.info(f"\n")
            logging.info(f"# =========================================== ")
            logging.info(f"# {fnam}: Begin at {tstr} ({Nsec})")            
            if MERGE_LAST:
                logging.info(f"# {fnam}: MERGE_LAST = {MERGE_LAST}")

        # need to re-compute output_dir to find submit info file
        output_dir,script_subdir = self.set_output_dir_name()
        output_dir               = self.override_output_dir_name(output_dir)
        self.config_prep['output_dir']  = output_dir

        # read SUBMIT.INFO passed from original submit job... 
        # this info never changes
        if verbose_flag:
            logging.info(f"# {fnam}: read {SUBMIT_INFO_FILE}")

        INFO_PATHFILE    = f"{output_dir}/{SUBMIT_INFO_FILE}"
        submit_info_yaml = util.extract_yaml(INFO_PATHFILE, None, None )
        self.config_prep['submit_info_yaml'] = submit_info_yaml

        # check option to reset merge process 
        if self.config_yaml['args'].merge_reset :
            self.merge_reset_driver()

        # make sure DONE_STAMP is stored in both config_prep and info_yaml
        # in case a prep function is called again during merge process.
        self.config_prep['done_stamp_list'] = submit_info_yaml['DONE_STAMP_LIST']
        self.config_prep['cleanup_flag']    = submit_info_yaml['CLEANUP_FLAG']

        # make sure time stamps are consistent
        self.merge_check_time_stamp(output_dir)

        # if last merge call (-M), then must wait for all of the done
        # files since there will be no more chances to merge.
        
        if MERGE_LAST : 
            self.merge_last_wait()
            if nomerge and not merge_background :
                self.nomerge_last()
                exit(0)

        # set busy lock file to prevent a simultaneous  merge task
        self.set_merge_busy_lock(+1,t_merge_start)

        # read status from MERGE file. There is one comment line above
        # each table to provide a human-readable header. The remaining
        # comments after the tables are probably merge-abort messages; 
        # these post-table comments are saved in comment_lines, and 
        # re-written below in write_merge_file function.

        if verbose_flag:
            logging.info(f"# {fnam}: examine {MERGE_LOG_FILE}")

        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        MERGE_INFO_CONTENTS, comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        self.merge_config_prep(output_dir)  # restore config_prep

        # make boolean list of which MERGED processes have already finished;
        # below will check which processes are newly DONE.
        done_list_start = self.get_merge_done_list(3,MERGE_INFO_CONTENTS)

        # check for changes in state for SPLIT and MERGE tables.
        # function returns updated set of SPLIT and MERGE table rows.

        row_list_dict, n_change = \
            self.merge_update_state(MERGE_INFO_CONTENTS)

        row_split_list = row_list_dict['row_split_list'] # optional
        row_merge_list = row_list_dict['row_merge_list'] # required
        row_extra_list = row_list_dict['row_extra_list'] # optional
        
        if not MERGE_LAST: 
            self.force_merge_failure(submit_info_yaml)

        use_split = len(row_split_list) > 0
        use_merge = len(row_merge_list) > 0
        use_extra = len(row_extra_list) > 0

        # Modify MERGE.LOG if there is a change in the processing STATE
        if n_change > 0 :
            itable = 0
            # redefine SPLIT table
            if use_split :
                INFO_STATE_SPLIT = {
                    'header_line' : comment_lines[itable],
                    'primary_key' : TABLE_SPLIT,
                    'row_list'    : row_split_list }
                itable += 1

            if use_extra :  # Oct 29 2021
                # table_names[0,1] are always SPLIT and MERGE
                # table_names[2] for extra table can have any name.
                table_name = row_list_dict['table_names'][2] # fragile warning!!
                INFO_STATE_EXTRA = {
                    'header_line' : comment_lines[itable],
                    'primary_key' : table_name,
                    'row_list'    : row_extra_list }
                itable += 1
            
            if use_merge :  # MERGE table must always be last
                INFO_STATE_MERGE = {
                    'header_line' : comment_lines[itable],
                    'primary_key' : TABLE_MERGE, 
                    'row_list'    : row_merge_list }
                itable += 1
                
            # re-write MERGE.LOG file with new set of STATEs. Note that
            # SPLIT & EXTRA tables are optional; MERGE table is required.
            # Any comment_lines after the tables are re-written so that
            # we don't lose information or merge-abort messages.
            with open(MERGE_LOG_PATHFILE, 'w') as f :
                if use_split :
                    util.write_merge_file(f, INFO_STATE_SPLIT, [] )

                if use_extra :
                    util.write_merge_file(f, INFO_STATE_EXTRA, [] )

                # note that merge table must be last because it includes
                # the comment lines for after the tables.
                if use_merge :
                    util.write_merge_file(f, INFO_STATE_MERGE, \
                                          comment_lines[itable:] )

            msg_update = f"Finished {n_change} STATE updates ({Nsec})."
        else:
            msg_update = f"No merge updates -> do nothing. "

        if verbose_flag:
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
        n_fail        = 0 
        n_wrapup      = 0
        for job_merge in range(0,n_job_merge):
            done_now     = done_list_end[job_merge]  # DONE or FAIL
            fail_now     = fail_list_end[job_merge]  # FAILED
            done_before  = done_list_start[job_merge] 
            done_new     = done_now and not done_before
            if done_now :
                n_done += 1
            if done_new and not fail_now and not check_abort :
                n_wrapup += 1
                self.merge_job_wrapup(job_merge,MERGE_INFO_CONTENTS)
            if fail_now:
                n_fail += 1

        if verbose_flag:
            logging.info(f"# {fnam}: finished {n_wrapup} wrapup tasks ")

        # Feb 2024: check kill_on_fail here in addition to check in CPU*.CMD
        if n_fail > 0 and submit_info_yaml['KILL_ON_FAIL'] :
            args.kill = True   # set flag as if -k were comman-line arg
            self.kill_jobs()

        # Only last merge process does cleanup tasks and DONE stamps
        if MERGE_LAST and n_done == n_job_merge :
            nfail_tot = self.failure_summary()

            if not check_abort :
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
        self.set_merge_busy_lock(-1,t_merge_start)
    
        self.merge_driver_exit(t_merge_start,False)
        return
        # end merge_driver

    def merge_driver_exit(self, t_merge_start, exit_flag):

        # Created Nov 15 2022
        # Write exit statement with timing info to enable 
        # checking unusually long merge time.
        #
        # Inputs:
        #   t_merge_start : time merge process started
        #   exit_flag     : if True, cal exit(0)
        #

        fnam = 'merge_driver_exit'
        t_merge_end = time.time()
        
        if t_merge_start is None:
            t_merge = 0
        else:
            t_merge     = t_merge_end - t_merge_start

        if exit_flag :
            action = 'skipped'
        else:
            action = 'completed'

        msg = f"{action} merge process after {t_merge:.2f} sec"
        logging.info(f"# {fnam}: {msg}")

        if exit_flag:  
            exit(0)

        return
        # end merge_graceful_exit

    def merge_reset_driver(self):

        # Created Apr 24 2021

        # parse batch info to enable kill_jobs
        self.parse_batch_info(self.config_yaml,self.config_prep)


        output_dir    = self.config_prep['output_dir']
        submit_mode   = self.config_prep['submit_mode'] 

        IS_SSH        = submit_mode == SUBMIT_MODE_SSH
        IS_BATCH      = submit_mode == SUBMIT_MODE_BATCH

        logging.info(f"# merge_rest: Execute merge_reset debug utility")

        # kill lingering jobs to avoid waiting monitor task to
        # mess things up.
        if IS_SSH :
            self.kill_ssh_jobs()

        elif IS_BATCH : 
            self.kill_sbatch_jobs()

        # call class-specific merge_reset function
        self.merge_reset(output_dir)
        sys.exit(f"\n Done with merge_reset. " \
                 f"Try merge again with -M option.")

        return

        # end merge_reset_driver

    def force_merge_failure(self,submit_info_yaml):

        # Apr 23 2021
        # check user option to force crash or force abort in merge process.

        # check option to force crash (to test higher level pipelines)
        if submit_info_yaml['FORCE_CRASH_MERGE']  :
            logging.info(f" xxx force merge crash with C-like printf xxx \n")
            printf(" xxx force merge crash with C-like printf xxx \n")

        # check option to force abort (to test higher level pipelines)
        if submit_info_yaml['FORCE_ABORT_MERGE']  :
            msgerr = []
            msgerr.append(f" Force ABORT in merge process;")
            msgerr.append(f" see --force_abort_merge argument")
            util.log_assert(False,msgerr)

        return

        # end force_merge_failure

    def merge_last_wait(self):

        # called after last science job (jobid = n_job_tot), but beware 
        # that other jobs may run later due to pending in queue. 
        # Here we wait for
        #  + all DONE files to exist
        #  + no BUSY files from other merge process
        # 
        # Nov 15 2022:
        #  Adjust logic to count DONE files and DONE files contained
        #  in already-compressed tar files.

        submit_info_yaml = self.config_prep['submit_info_yaml'] 
        n_job_tot        = submit_info_yaml['N_JOB_TOT']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        n_done_tot       = submit_info_yaml['N_DONE_TOT']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_dir       = submit_info_yaml['SCRIPT_DIR'] 
        done_file_wildcard    = f"{jobfile_wildcard}.DONE"
        done_tar_wildcard     = f"{jobfile_wildcard}DONE.tar.gz"

        #util.wait_for_files(n_done_tot, script_dir, done_wildcard) 
        wait_for_all_done = True
        ntry = 0
        while wait_for_all_done :
            done_file_list  = glob.glob1(script_dir, done_file_wildcard)
            done_tar_list   = glob.glob1(script_dir, done_tar_wildcard)
            n_done_file     = len(done_file_list)
            n_done_tar      = len(done_tar_list) * n_job_split
            n_done          = n_done_file + n_done_tar
            time_now        = datetime.datetime.now()
            tstr            = time_now.strftime("%Y-%m-%d %H:%M:%S") 
            msg = f"\t Found {n_done} of {n_done_tot} DONE files ({tstr})"
            logging.info(msg)

            ntry += 1
            if n_done != n_done_tot :
                if ntry > 1: time.sleep(20)
            else:
                msg = f"\t    ({n_done_file} from DONE files, " \
                      f"{n_done_tar} from tar files)"
                logging.info(msg)
                wait_for_all_done = False
                
        # - - - - - - - - - -  -
        time.sleep(1)

        # sleep until there are no more busy files.
        n_busy,busy_list = self.get_busy_list()
        while n_busy > 0 :
            logging.info(f"\t Wait for {busy_list} to clear")
            time.sleep(5)
            n_busy,busy_list = self.get_busy_list()

        logging.info("")
        return

        # end merge_last_wait

    def nomerge_last(self):

        # Created Mar 28 2022
        # for --nomerge option, create RUN_MERGE_[inputFile] script and
        # create tar file of outdir. The output RUN_MERGE script is a
        # debug tool to quickly run the merge process interactively.
        #
        # Jan 23 2025: fix to work if outdir is not in same place in submit-input file.
        # Feb 04 2025: fix Jan 23 bug when output_dir has no slashes

        output_dir          = self.config_prep['output_dir'] 
        args                = self.config_yaml['args']
        submit_info_yaml    = self.config_prep['submit_info_yaml'] 
        input_file          = args.input_file

        #xxx mark output_dir_path = os.path.dirname(output_dir)
        output_dir_path = output_dir + '/..'
        output_dir_base = os.path.basename(output_dir)
        input_file_base = os.path.basename(input_file.split('.')[0])

        merge_script   = "RUN_MERGE_" + input_file_base
        program_submit = sys.argv[0]  # $path/submit_batch_jobs

        merge_mode = submit_info_yaml['MERGE_MODE']
        if merge_mode == MERGE_MODE_BACKGROUND:
            return

        logging.info(f"\n Create {merge_script} ")
        with open(merge_script,"wt") as s:
            s.write(f"# test merge process interactively. \n")
            s.write(f"echo 'Update OUTDIR before merge process ...' \n")
            s.write(f"cd {output_dir_path} \n")
            s.write(f"rm -r   {output_dir_base}\n")
            s.write(f"tar -xf {output_dir_base}.tar\n")
            s.write(f"\n")
            s.write(f"cd {CWD}\n")            
            s.write(f"{program_submit} \\\n")
            s.write(f"  {input_file}   \\\n")
            s.write(f"  -M \n")

        cmd = f"chmod +x {merge_script}"
        os.system(cmd)

        # - - - - - - - - - - - - 
        # create tar file of output dir
        tar_file      = f"{output_dir_base}.tar"
        tar_file_path = f"{output_dir_path}/{output_dir_base}.tar"        

        # remove old tar file if it exists
        if os.path.exists(tar_file_path):  os.remove(tar_file_path)

        cmd = f"cd {output_dir_path} ;  tar -cf {tar_file} {output_dir_base}"
        logging.info(f" Create {tar_file}")
        logging.info(f" with command: {cmd}")
        os.system(cmd)        

        logging.info(f" Done.")

        return
        # end nomerge_last

    def merge_cleanup_script_dir(self):
        # Tar of CPU* files, then tar+gzip script_dir
        # Beware that after this task, nothing more can be written
        # to the CPU*LOG files.

        cpunum            = self.config_yaml['args'].cpunum[0]
        submit_info_yaml  = self.config_prep['submit_info_yaml']
        script_dir        = submit_info_yaml['SCRIPT_DIR'] 

        log_file_keep  = f"CPU{cpunum:04d}*.LOG"
        logging.info(f"  Standard cleanup: compress CPU* files " \
                     f"except for {log_file_keep}")
        util.compress_files(+1, script_dir, "CPU*", "CPU", log_file_keep )

        # tar and zip the script dir
        logging.info(f"  Standard cleanup: compress {script_dir}")
        util.compress_subdir(+1, f"{script_dir}" )

        # merge_cleanup_script_dir

    def set_merge_busy_lock(self,flag, t_merge_start):

        # Inputs:
        #   flag > 0 --> create BUSY LOCK file
        #   flag < 0 --> remove it
        #
        #   t_merge_start : used to print process time upon exit

        Nsec  = seconds_since_midnight  # current Nsec, not from submit info
        t_msg = f"T_midnight={Nsec}"
        output_dir   = self.config_prep['output_dir']
        args         = self.config_yaml['args']
        check_abort  = args.check_abort
        merge_force  = args.merge_force  # wait for BUSY to clear

        merge_normal = not merge_force # exit if BUSY elsewhere
        verbose_flag = not args.check_abort
        fnam = "merge_driver"

        if self.config_yaml['args'].cpunum :
            cpunum = self.config_yaml['args'].cpunum[0] # passed from script
        else:
            return 

        busy_file     = f"{BUSY_FILE_PREFIX}{cpunum:04d}.{BUSY_FILE_SUFFIX}"
        BUSY_FILE     = f"{output_dir}/{busy_file}"
        n_busy,busy_list = self.get_busy_list()

        if flag > 0 :
            # check for other busy files to avoid conflict
            if merge_normal and n_busy > 0 :
                msg = f"\n# {fnam}: Found existing " \
                      f"{busy_list[0]} --> skip merge process."
                logging.info(msg)
                self.merge_driver_exit(t_merge_start, True)
                # xxx mark sys.exit(msg)  
            else: 
                if merge_force :
                    while n_busy > 0:
                        t_now = datetime.datetime.now()
                        msg = f"# {fnam}: \t merge_force -> " \
                              f"wait for BUSY to clear ({t_now})."
                        logging.info(msg)
                        time.sleep(1)
                        n_busy, busy_list = self.get_busy_list()

                msg = f"# {fnam}: \t Create {busy_file} for {t_msg}"
                if verbose_flag: logging.info(msg)
                with open(BUSY_FILE,"w") as f:
                    f.write(f"{Nsec}\n")  # maybe useful for debug situation

                # check for simultaneously created busy file(s);
                # --> avoid conflict by keeping only first one in sort list
                time.sleep(1)
                n_busy,busy_list = self.get_busy_list()
                if n_busy > 1 and busy_file != busy_list[0] :
                    cmd_rm = f"rm {BUSY_FILE}"
                    os.system(cmd_rm)
                    msg = f"\n# {fnam}: Found simultaneous " \
                          f"{busy_list[0]} --> skip merge process."
                    logging.info(msg)
                    self.merge_driver_exit(t_merge_start, True)
                    # xxx mark sys.exit(msg)  

        elif len(BUSY_FILE)>5 and os.path.exists(BUSY_FILE):  # avoid rm *
            if verbose_flag:
                logging.info(f"# {fnam}: " \
                             f"\t Remove {busy_file} for {t_msg}")
            cmd_rm = f"rm {BUSY_FILE}"
            os.system(cmd_rm)

        # end set_merge_busy_lock

    def get_busy_list(self):
        # return numbrer of busy files,  and sorted list of busy files.
        output_dir      = self.config_prep['output_dir']
        busy_wildcard   = f"{BUSY_FILE_PREFIX}*.{BUSY_FILE_SUFFIX}"
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
        #   submit_batch_jobs.sh <inFile> -m
        if not self.config_yaml['args'].t :
            return

        # time_stamp and cpunum args are given only in the CPU*CMD files
        time_stamp_merge  = self.config_yaml['args'].t[0]
        cpunum            = self.config_yaml['args'].cpunum[0]

        if time_stamp_submit != time_stamp_merge :
            stop_file = f"{output_dir}/STOP_CPU{cpunum:04d}_{Nsec_now}.INFO"
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

        args                = self.config_yaml['args']
        output_dir          = self.config_prep['output_dir']
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        
        #print(f" xxx MERGE_LOG_PATHFILE = {MERGE_LOG_PATHFILE} ")

        #  append to bottom of MERGE.LOG
        with open(MERGE_LOG_PATHFILE,"a") as f:
            # write class name (Jan 2024)
            f.write(f"PROGRAM_CLASS:  {args.program_class}\n")

            # write class-specific information
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
        n_core_with_jobs  = submit_info_yaml['N_CORE_WITH_JOBS']
        time_submit       = submit_info_yaml['TIME_STAMP_SUBMIT']
        time_now          = datetime.datetime.now()
        time_dif          = time_now - time_submit

        output_dir          = self.config_prep['output_dir']
        MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"
        
        t_seconds = time_dif.total_seconds()
        t_unit    = 3600.0 ;    unit = "hours"

        t_wall   = t_seconds/t_unit
        msg_time = [ ]
        msg_time.append(f"UNIT_TIME:      {unit}   ")
        msg_time.append(f"WALL_TIME:      {t_wall:.3f}  # {unit}  ")

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
            LOG_FILE = f"{script_dir}/{cpu_log_file}"
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
                        t_str            = f"{word_list[1]} {word_list[2]}"
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
            cpu_avg  = cpu_sum / n_core_with_jobs
            eff_cpu  = cpu_avg / t_wall

            msg_time.append(f"CPU_SUM:        {cpu_sum:.3f}   # {unit}")
            msg_time.append(f"EFFIC_CPU:      {eff_cpu:.3f}   " \
                            f"# CPU/core/T_wall")

        return msg_time

        # end get_proctime_info

    def check_file_exists(self,file_name,msg_user):
        # abort if file does not exist and use self.log_assert
        # that writes error message to MERGE.LOG and FAIL to done stamp
        exist = os.path.isfile(file_name)
        if not exist:
            msgerr = [ f"Cannot find file:", f"\t {file_name}" ]
            for msg in msg_user:
                msgerr.append(msg)
            self.log_assert(exist, msgerr)
    # end check_file_exists

    def check_for_failure(self, log_file, nevt, bad_output, isplit):

        # Top-level function to check for failures. External
        # program must determine nevt.
        # nevt <=0 means job failed, so call faulure_update
        # function that will update FAILURES.LOG
        # This function returns true if failure is identified.
        #
        # Mar 26 2021: create=True always. isplit<4 missing too much.
        # Feb 06 2025: pass bad_output logical

        submit_info_yaml = self.config_prep['submit_info_yaml']
        cleanup_flag     = submit_info_yaml['CLEANUP_FLAG']
        MXSPLIT_FAIL_REPEAT = 4       # max isplit to make FAIL_REPEAT script
        
        fail_no_output   = (nevt <  0)  # no output
        fail_zero_evt    = (nevt == 0)  # zero events 
        found_fail       = (fail_no_output or fail_zero_evt)
        create           = True   # Mar 26 2021

        if found_fail or bad_output :
            cleanup_flag   = 0 # STOP all cleanup activities
            self.failure_update(log_file, create, fail_zero_evt, bad_output )

        return found_fail

        # end check_for_failure
    
    def failure_update(self, job_log_file, create_script, found_zero, bad_output ):

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
        #   bad_output    : T -> bad output detected by science code
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

        FAIL_SUMMARY_PATHFILE  = f"{script_dir}/{FAIL_SUMMARY_FILE}"
        JOB_LOG_PATHFILE       = f"{script_dir}/{job_log_file}"

        # on 1st failure, create fail log
        if not os.path.isfile(FAIL_SUMMARY_PATHFILE):
            msg = f"\n FIRST JOB FAILURE -> Create {FAIL_SUMMARY_FILE}"
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
        sim_prefix = f"TMP_{USER4}_"
        if sim_prefix in job_log_file:
            j_start = len(sim_prefix)
        
        fail_repeat_script = f"FAIL_REPEAT_{job_log_file[j_start:j_dot]}.CMD"
        FAIL_REPEAT_SCRIPT = f"{script_dir}/{fail_repeat_script}"

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
        elif bad_output :
            FAIL_INDEX = FAIL_INDEX_BAD    # Feb 2025
        else :
            FAIL_INDEX = FAIL_INDEX_UNKNOWN

        FAIL_MODE  = FAIL_MODE_STRINGS[FAIL_INDEX]

        # increment total number of failures
        n_job_fail_list[0] += 1
        n_job_fail          = n_job_fail_list[0]
        msg = f"    *** FAIL {n_job_fail:3d} for " \
              f"{job_log_file}   FAIL_MODE={FAIL_MODE}  " \
              f"found_zero={found_zero}  bad_output={bad_output}"
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
        FAIL_SUMMARY_PATHFILE = f"{script_dir}/{FAIL_SUMMARY_FILE}"

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
            CMD_FILE = f"{script_dir}/{cmd_file}"
            with open(CMD_FILE,"r") as f_cmd :
                for line in f_cmd :
                    if util.is_comment_line(line) : continue
                    if len(line) == 0 : continue
                    #if len(line) > 0 and line[0] != '#' :
                    command_lines.append(line.rstrip("\n"))

        nlines = len(command_lines)
        msg = f"  Stored {nlines} command lines from " \
              f"{CMD_wildcard} files"
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
        FAIL_SUMMARY_PATHFILE = f"{script_dir}/{FAIL_SUMMARY_FILE}"
        MERGE_LOG_PATHFILE    = f"{output_dir}/{MERGE_LOG_FILE}"

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
            key     = f"NJOB_FAIL({mode}):"
            string  = f"   {key:<25}  {nfail:2d}  # {comment}"
            string_insert += f"{string}\\\n"  # for sed command
            
        # construct sed command to pre-pend summary at top of file
        # ?? can python do this ??
        cmd_sed = f"sed -i '1i{string_insert}'"

        if found_failures :
            CMD = f"{cmd_sed} {FAIL_SUMMARY_PATHFILE}"
            os.system(CMD)

            # if script_dir is not output_dir, copy FAILURES.LOG to
            # output_dir so that it's easier to find
            if script_dir != output_dir :
                shutil.copy(FAIL_SUMMARY_PATHFILE,output_dir)

        # always update merge log file
        CMD = f"{cmd_sed} {MERGE_LOG_PATHFILE}"
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
        # Feb 2025: check for BAD_OUTPUT key

        n_key              = len(yaml_key_list)
        n_split            = len(log_file_list)
        key_AIZ            = 'ABORT_IF_ZERO'
        key_BAD            = 'BAD_OUTPUT'   
        aiz_list           = [ -9 ]    * n_split
        bad_list           = [ False ] * n_split
        aiz_max            = 0
        n_aiz_zero         = 0
        n_bad              = 0

        # init output dictionary
        job_stats = { 'nfail' : 0 }
        for key in yaml_key_list :
            key_sum  = f"{key}_sum"   # store sum
            key_list = f"{key}_list"  # store stats for each split job
            job_stats[key_list] = [ 0 ] * n_split
            job_stats[key_sum]  =   0
    
        for isplit in range(0,n_split):
            log_file  = log_file_list[isplit]
            yaml_file = yaml_file_list[isplit]
            LOG_FILE  = f"{script_dir}/{log_file}"
            YAML_FILE = f"{script_dir}/{yaml_file}"

            if os.path.isfile(YAML_FILE) :
                stats_yaml       = util.extract_yaml(YAML_FILE, None, None )
                aiz              = stats_yaml[key_AIZ]                    # required key
                bad              = stats_yaml.setdefault(key_BAD,False)   # optional key (Feb 2025)

                aiz_list[isplit] = aiz
                bad_list[isplit] = bad
                if aiz > aiz_max : aiz_max = aiz
                if aiz == 0 : n_aiz_zero += 1
                if bad      : n_bad += 1
            
                for item in yaml_key_list :
                    key, key_sum, key_list  = self.keynames_for_job_stats(item)
                    self.log_assert(key in stats_yaml,
                                    [ f"Missing yaml key={key}", f"in {YAML_FILE}"] )

                    job_stats[key_list][isplit] = stats_yaml[key]
                    job_stats[key_sum]         += stats_yaml[key]
                
                    # fix format for CPU
                    if 'CPU' in key :
                        cpu_sum = f"{job_stats[key_sum]:.2f}"
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
            msg = f"\t max({key_AIZ})={aiz_max} -> suppress abort from AIZ=0."
            logging.info(msg)

        # loop again over split jobs and check for failures.
        for isplit in range(0,n_split):
            log_file   = log_file_list[isplit]
            aiz        = aiz_list[isplit]            
            bad        = bad_list[isplit]            
            found_fail =  self.check_for_failure(log_file, aiz, bad, isplit+1)
            if found_fail : job_stats['nfail'] += 1

        return job_stats

        # end get_job_stats

    def keynames_for_job_stats(self,keyname_base):
        keyname_sum  = f"{keyname_base}_sum"
        keyname_list = f"{keyname_base}_list"
        return keyname_base, keyname_sum, keyname_list
        # end keynames_for_job_stats

    def log_assert(self,condition, msgerr):
        # same as log_assert in util, except here it also
        # + writes FAIL to done_stamp_file
        # + appends msgerr to MERGE.LOG
        # + remove BUSY file

        if condition is False :
            output_dir      = self.config_prep['output_dir']
            done_stamp_list = self.config_prep['done_stamp_list']
            MERGE_LOG_PATHFILE  = f"{output_dir}/{MERGE_LOG_FILE}"

            util.write_done_stamp(output_dir, done_stamp_list,STRING_FAIL)

            if os.path.isfile(MERGE_LOG_PATHFILE) :
                with open(MERGE_LOG_PATHFILE,"a") as f:
                    f.write(f"# {SNANA_ABORT_STRING}\n")
                    for msg in msgerr :
                        f.write(f"#   {msg}\n")
                    f.write(f"#\n")

            self.set_merge_busy_lock(-1, None)  # remove busy file Apr 24 2021
            util.log_assert(condition, msgerr)        
            return

# ======= END OF FILE =========

