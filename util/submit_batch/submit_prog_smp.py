# Created June 2025 by R.Kessler and C. Meldorf
# Initial use is MultiSMP for Roman.
#

import os, sys, logging, datetime, time, subprocess, yaml
import submit_util as util
import pandas as pd
import numpy  as np

from   submit_params    import *
from   submit_prog_base import Program

# ===================================================


# not sure we need this subclass dictionary, but maybe later if there are multiple
# SMP implementations
KEY_SUBCLASS_DICT = {
    'campari'  : [] 
}


FITOPT_STRING = "FITOPT"

# define columns in merge file
COLNUM_SMP_MERGE_STATE           = 0  # STATE required in col=0
COLNUM_SMP_MERGE_FITOPT          = 1
COLNUM_SMP_MERGE_ISPLIT          = 2
COLNUM_SMP_MERGE_NDET            = 3  # number of detectors (not detections)
COLNUM_SMP_MERGE_NLC             = 4  # number of light curves
COLNUM_SMP_MERGE_NOBS            = 5  # number of observations for which flux is determined
COLNUM_SMP_MERGE_CPU             = 6


# ====================================================
#    BEGIN CLASS
# ====================================================


class SceneModelPhotometry(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_SMP  # default; override later with optional JOBNAME key
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

        return output_dir_name, SUBDIR_SCRIPTS_SMP

        
    def submit_prepare_driver(self):

        # init driver for smp
        self.prep_smp_subclass()
        self.prep_smp_fitopt()
        self.prep_smp_healpix()
        self.prep_smp_flatten()
        self.prep_smp_output_dirs()  

        # load misc globals used by base program

        n_fitopt = self.config_prep['fitopt_dict']['n_jobopt']
        self.config_prep['n_job_tot']   = n_fitopt
        self.config_prep['n_done_tot']  = n_fitopt
        
        return
        # end submit_prepare_driver

    def prep_smp_subclass(self):
        SMP_SUBCLASS = self.config_yaml['CONFIG']['SMP_METHOD']
        msg = f"\n\t !!!! SMP_SUBCLASS = {SMP_SUBCLASS} !!!! \n"
        logging.info(msg)

    def prep_smp_fitopt(self):

        CONFIG        = self.config_yaml['CONFIG']
        KEYLIST       = [ FITOPT_STRING ]    # optional key under CONFIG 
        fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)  # parse FITOPTs

        # convert FITOPT string to dictionary; element 000 (no options) is
        # automatically pre-pended as first FITOPT
        fitopt_dict   = util.prep_jobopt_list(fitopt_rows, FITOPT_STRING, 1, None)

        n_fitopt = fitopt_dict['n_jobopt']
        logging.info(f"\t Found {n_fitopt}  FITOPTs")
    
        self.config_prep['fitopt_dict'] = fitopt_dict
        return



    def prep_smp_output_dirs(self):

        # prepare output dir for each fitopt; passed to campari code
    
        output_dir        = self.config_prep['output_dir']
        iopt_list         = self.config_prep['iopt_list']  
        isplit_list       = self.config_prep['isplit_list']

        fitopt_dict       = self.config_prep['fitopt_dict']
        fitopt_num_list   = fitopt_dict['jobopt_num_list']  # e.g., "FITOPT001"=

        output_dir_smp_list = []
        logging.info('')
        logging.info(f"  prepare SMP outdir(s):")        
        for iopt, isplit in zip(iopt_list,isplit_list):
            fitopt_num = fitopt_num_list[iopt]
            prefix         = self.smp_prefix(fitopt_num,isplit) 
            output_dir_smp = f"{output_dir}/OUTPUT_{prefix}"
            output_dir_smp_list.append(output_dir_smp)
            logging.info(f"\t {output_dir_smp}")
            os.mkdir(output_dir_smp)

        logging.info('')
        self.config_prep['output_dir_smp_list'] = output_dir_smp_list
        
        return
    # end prep_smp_output_dirs

    def prep_smp_healpix(self):

        CONFIG        = self.config_yaml['CONFIG']
        n_fitopt      = self.config_prep['fitopt_dict']['n_jobopt']
        n_core        = self.config_prep['n_core']

        n_hpsubset    = int(n_core/n_fitopt)  # number of healpix subsets
        if n_hpsubset < 1: n_hpsubset = 1
        n_job_split   = n_hpsubset
        
        # revise n_core
        n_core = n_hpsubset * n_fitopt
        self.config_prep['n_core'] = n_core
        
        logging.info(f" Calculated nubmer of healpix subsets: {n_hpsubset}")
        logging.info(f" Revised n_core: {n_core} = " \
                     f"{n_fitopt}(FITOPT) x {n_hpsubset}(HEALPIX SUBSETS)")

        # read full list of healpixels
        key = 'HEALPIX_LIST_FILE' 
        if key not in CONFIG:
            sys.exit(f"\n ABORT: missing required key = {key}\n CONFIG = \n{CONFIG}") # use official abort function

        HEALPIX_LIST_FILE = CONFIG[key]
        with open(HEALPIX_LIST_FILE,"rt") as h:
            lines = h.readlines()
            
        contents = yaml.safe_load("\n".join(lines))
        nside        = contents['NSIDE']
        healpix_global_list = contents['HEALPIX']
        n_healpix_global = len(healpix_global_list)
        logging.info(f" Found {n_healpix_global} total healpixels to process in {HEALPIX_LIST_FILE}")
        # split healpix_list into n_job_split sub-lists
        # (welcome AI overlords for figuring this out)
        split_list = [list(arr) for arr in \
                      np.array_split(np.array(healpix_global_list), n_job_split)]

        # construct separate healpix list file for each task/list
        healpix_file_list = []
        for i_split, healpix_list in enumerate(split_list):
            hp_file = self.write_healpix_sublist(i_split, nside,healpix_list)
            healpix_file_list.append(hp_file)
            
        self.config_prep['n_job_split']       = n_job_split
        self.config_prep['healpix_file_list'] = healpix_file_list
        return
        # end prep_smp_healpix

    def prep_smp_flatten(self):
        n_job_split      = self.config_prep['n_job_split']
        fitopt_dict      = self.config_prep['fitopt_dict']
        n_fitopt         = fitopt_dict['n_jobopt']

        iopt_list   = []  # for FITOPT
        isplit_list = []

        for iopt in range(0,n_fitopt):
            for isplit in range(0,n_job_split):
                iopt_list.append(iopt)
                isplit_list.append(isplit)

        self.config_prep['iopt_list']   = iopt_list
        self.config_prep['isplit_list'] = isplit_list        

        return
    
    def write_healpix_sublist(self, i_split, nside, healpix_list):
                                            
        script_dir       = self.config_prep['script_dir']
        
        n_healpix = len(healpix_list)
        hp_file = f"HEALPIX_SPLIT{i_split:03d}.DAT"
        hp_path = f"{script_dir}/{hp_file}"
        logging.info(f"\t write {hp_file} with {n_healpix} healpixels")

        with open(hp_path,"wt") as h:
            h.write(f"Healpix sub-list for SPLIT job {i_split}\n\n")
            h.write(f"NSIDE: {nside} \n")
            h.write(f"\nHEALPIX: \n")
            for hp in healpix_list:
                h.write(f"- {hp}\n")
                
        return hp_path
    
        # end write_healpix_sublist
        
    def write_command_file(self, icpu, f):

        # For this icpu, write full set of sim commands to 
        # already-opened command file with pointer f.
        # Function returns number of jobs for this cpu
        
        n_job_tot        = self.config_prep['n_job_tot']
        n_core           = self.config_prep['n_core']
        fitopt_dict      = self.config_prep['fitopt_dict']
        n_job_split      = self.config_prep['n_job_split']

        iopt_list        = self.config_prep['iopt_list']
        isplit_list      = self.config_prep['isplit_list']
        
        n_fitopt         = fitopt_dict['n_jobopt']
        n_job_local       = 0
        
        for iopt, isplit in zip(iopt_list, isplit_list):
            
            index_dict = { 'iopt':iopt, 'isplit':isplit, 'icpu':icpu  }

            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :        
                job_info_exec  = self.prep_JOB_INFO_smp(index_dict)
                util.write_job_info(f, job_info_exec, icpu)

                #job_info_merge = self.prep_JOB_INFO_merge(icpu, n_job_local, False)
                #util.write_jobmerge_info(f, job_info_merge, icpu)

        # - - - - - -
        return        
        # end write_command_file
        

    def prep_JOB_INFO_smp(self, index_dict):

        iopt   = index_dict['iopt']  # this is the FITOPT
        isplit = index_dict['isplit']
        icpu   = index_dict['icpu']
        
        program             = self.config_prep['program']
        output_dir          = self.config_prep['output_dir']
        fitopt_dict         = self.config_prep['fitopt_dict']
        healpix_file_list   = self.config_prep['healpix_file_list']
        output_dir_smp_list = self.config_prep['output_dir_smp_list']
        args                = self.config_yaml['args']  # user args to submit_batch
        
        fitopt_num     = fitopt_dict['jobopt_num_list'][iopt]
        fitopt_arg     = fitopt_dict['jobopt_arg_list'][iopt]
        fitopt_label   = fitopt_dict['jobopt_label_list'][iopt]
        output_dir_smp = output_dir_smp_list[iopt]
        healpix_file   = healpix_file_list[isplit]

        if fitopt_label is None:
            fitopt_label = 'NOLABEL'
            
        prefix    = self.smp_prefix(fitopt_num, isplit)
        log_file  = f"{prefix}.LOG"
        done_file = f"{prefix}.DONE"
        yaml_file = f"{prefix}.YAML"
        
        JOB_INFO = {}
        JOB_INFO['job_dir']      = output_dir  # where to run job 
        JOB_INFO['program']      = program
        JOB_INFO['input_file']   = args.input_file
        
        arg_list = []
        arg_list.append(fitopt_arg)

        # append required yaml_file
        yaml_arg = f"--outfile_yaml {yaml_file}"
        arg_list.append(yaml_arg)

        # specify outdir for this FITOPT
        output_arg     = f"photometry-campari-paths-output_dir {output_dir_smp}"
        arg_list.append(output_arg)
        
        # specify file with healpix subset
        hp_subset_arg = f"photometry-campari-paths-healpix-file {healpix_file}"  # RK guess; ask Cole
        arg_list.append(hp_subset_arg)

        JOB_INFO['log_file']      = log_file
        JOB_INFO['done_file']     = done_file
        JOB_INFO['arg_list']      = arg_list
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['kill_on_fail']  = args.kill_on_fail
        JOB_INFO['check_abort']   = args.check_abort
        
        return JOB_INFO

	# end prep_JOB_INFO_smp        

    def smp_prefix(self, fitopt_num, isplit):
        prefix = f"SMP_{fitopt_num}_SPLIT{isplit:03d}"
        return prefix
    
    # ========== MERGE STUFF =============

    def append_info_file(self,f):
        # append info to SUBMIT.INFO file  

        fitopt_dict       = self.config_prep['fitopt_dict']
        n_fitopt          = fitopt_dict['n_jobopt']
        
        f.write(f"JOBFILE_WILDCARD:  'SMP*' \n")
        f.write(f"\n# {PROGRAM_NAME_SMP} info\n")
        
        f.write("\n")
        f.write(f"FITOPT_LIST: \n")
        row = [ 0 ] * 3
        for i in range(0,n_fitopt):
            row[COLNUM_FITOPT_NUM]   =  fitopt_dict['jobopt_num_list'][i]
            row[COLNUM_FITOPT_LABEL] =  fitopt_dict['jobopt_label_list'][i]
            row[COLNUM_FITOPT_ARG]   =  fitopt_dict['jobopt_arg_list'][i]
            f.write(f"  - {row} \n")
        f.write("\n")


        healpix_file_list = self.config_prep['healpix_file_list']
        f.write(f"\nHEALPIX_FILE_LIST: \n")
        for h in healpix_file_list:
            f.write(f"  - {h}\n")
            
        return
    # end append_info_file
               

    def create_merge_table(self,f):

        iopt_list         = self.config_prep['iopt_list']
        isplit_list       = self.config_prep['isplit_list']
        fitopt_dict       = self.config_prep['fitopt_dict']        
        n_fitopt          = fitopt_dict['n_jobopt']
        fitopt_num_list   = fitopt_dict['jobopt_num_list']  # e.g., "FITOPT001"    

        header_line = " STATE   FITOPT  ISPLIT  NDET  NLC  NOBS   CPU "
        MERGE_INFO = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line,
            'row_list'      : []
        }
        
        for iopt, isplit in zip(iopt_list, isplit_list):
            fitopt_num = fitopt_num_list[iopt]
            ROW = [ SUBMIT_STATE_WAIT, fitopt_num,  isplit, 0,0,0,0 ]
            MERGE_INFO['row_list'].append(ROW)

        util.write_merge_file(f, MERGE_INFO, [] )
        return
        # end create_merge_table

        
    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.                  
        # Return update row list MERGE tables.     

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']

        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        #header_line = " STATE   FITOPT  ISPLIT  NDET  NLC  NOBS   CPU "        

        key_ndet, key_ndet_sum, key_ndet_list = \
                self.keynames_for_job_stats('NDET')
        key_nlc, key_nlc_sum, key_nlc_list = \
                self.keynames_for_job_stats('NLC')
        key_nobs, key_nobs_sum, key_nobs_list = \
                self.keynames_for_job_stats('NOBS')   
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')

        key_list = [ key_ndet, key_nlc, key_nobs, key_cpu ]
        row_list_merge   = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # init outputs of function                  
        n_state_change     = 0
        row_list_merge_new = []

        nrow_check = 0
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input     
            nrow_check += 1
            irow        = nrow_check - 1 # row index
            # strip off row info                   
            STATE       = row[COLNUM_SMP_MERGE_STATE]
            FITOPT_NUM  = row[COLNUM_SMP_MERGE_FITOPT]
            ISPLIT      = row[COLNUM_SMP_MERGE_ISPLIT]
            prefix      = self.smp_prefix(FITOPT_NUM, ISPLIT ) 
            search_wildcard = f"{prefix}*"

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

                    fit_stats = self.get_job_stats(script_dir,
                                                   log_list, yaml_list, key_list)

                    # check for failures in snlc_fit jobs.       
                    nfail = fit_stats['nfail']
                    if nfail > 0 :  NEW_STATE = SUBMIT_STATE_FAIL

                    row[COLNUM_SMP_MERGE_STATE]     = NEW_STATE
                    row[COLNUM_SMP_MERGE_NDET]      = fit_stats[key_ndet_sum]
                    row[COLNUM_SMP_MERGE_NLC]       = fit_stats[key_nlc_sum]
                    row[COLNUM_SMP_MERGE_NOBS]      = fit_stats[key_nobs_sum]  
                    row[COLNUM_SMP_MERGE_CPU]       = fit_stats[key_cpu_sum]

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             # assume nevt changes
        # - - - - -
        row_list_dict = {
            'row_split_list' : [],   # no split table, except for sim
            'row_merge_list' : row_list_merge_new,
            'row_extra_list' : []    # no extras
        }

        return row_list_dict, n_state_change
        # end merge_update_state


    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        row     = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]

    def merge_cleanup_final(self):
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_subdir    = SUBDIR_SCRIPTS_COSMOFIT

        logging.info(f"  SMP cleanup for {JOB_SUFFIX_TAR_LIST}")
        for suffix in JOB_SUFFIX_TAR_LIST:
            wildcard = (f"SMP*.{suffix}")
            util.compress_files(+1, script_dir, wildcard, suffix, "" )

        util.compress_files(+1, script_dir, "CPU*", "CPU", "" )

        logging.info("")
        # end merge_cleanup_final                        

    def merge_config_prep(self,output_dir):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        # .xyz restore SUBMIT.INFO information needed in MERGE process
        return

   
    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file        
        submit_info_yaml = self.config_prep['submit_info_yaml']
        info_lines = []

        # .xyz any SMP-specific information to append to MERGE.LOG ?
        # example:
        info_lines.append(f"DUMMY_SMP_INFO:  replace-this-with-something-useful")
        
        return info_lines
        # end get_misc_merge_info

    def get_merge_COLNUM_CPU(self):
        # used by base program to sum CPU over all cores
        return COLNUM_SMP_MERGE_CPU
