# Created June 2025 by R.Kessler
# class to run combine_fitres.exe as part of pippin pipeline;
# in particular, to combine classifier probabilities (SCONE, SNN ....)
# with the LCFIT table output, so that it is ready for BBC.
# Pippin did this task interactively, and for DES-SN5YR it takes
# about 4-5 hr of wall time ... dividing this task among cores
# should greatly reduce the wall time.
#
# Feb 2 2026: fix bug; copy ALL.DONE when all is finished at end of merge_cleanup_final
#          ( not during prep stage)

import os, sys, shutil, yaml, glob
import logging
#import coloredlogs
import datetime, time
import submit_util      as util
from   submit_params    import *
from   submit_prog_base import Program

COLNUM_COMBINE_MERGE_JOBNUM   = 1
COLNUM_COMBINE_MERGE_JOBNAME  = 2
COLNUM_COMBINE_MERGE_NFILE    = 3
COLNUM_COMBINE_MERGE_NEVT     = 4
COLNUM_COMBINE_MERGE_CPU      = 5

# define files produced by submit_batch_jobs that are needed
# to mimic another output dir from submit_batch_jobs
# Make it look like real job by copying SUBMIT.INFO first,
# then copy other files when all is finished.

STAGE_PREP = "PREP"
STAGE_DONE = "DONE"
MIMIC_FILE_LIST      = [SUBMIT_INFO_FILE, MERGE_LOG_FILE, DEFAULT_DONE_FILE ]
MIMIC_STAGE_LIST     = [STAGE_PREP,       STAGE_DONE,     STAGE_DONE        ]

# =============================================================
class combine_fitres(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_COMBINE_FITRES
        super().__init__(config_yaml, config_prep)

    def set_output_dir_name(self):
        CONFIG     = self.config_yaml['CONFIG']
        input_file = self.config_yaml['args'].input_file  # for msgerr
        msgerr     = []

        if 'LOGDIR' in CONFIG :
            OUTDIR = CONFIG['LOGDIR']
            output_dir_name = os.path.expandvars(OUTDIR)
        else:
            msgerr.append(f"LOGDIR key missing in yaml-CONFIG")
            msgerr.append(f"Check {input_file}")
            log_assert(False,msgerr)

        return output_dir_name, SUBDIR_SCRIPTS_COMBINE
        # end set_output_dir_name

    def submit_prepare_driver(self):

        self.prep_combine_inputs()
        self.prep_combine_outdirs()
        self.prep_combine_mimic_outdirs(STAGE_PREP)
        self.prep_combine_symlinks()
        return

    def prep_combine_inputs(self):
    
        # check that all input files and directories exist;
        # abort on missing inputs. Make sure to print all missing
        # items before aborting.

        COMBINE_FITRES = self.config_yaml['COMBINE_FITRES']
        TASK_LIST    = COMBINE_FITRES['TASK_LIST']
        INPUT_TOPDIR = COMBINE_FITRES.setdefault('INPUT_TOPDIR',None)

        nerr = 0
        n_symlink = 0
        n_job_tot = 0
        n_inp_base_tot = 0
        n_inp_append_tot = 0
        combine_arg_list = []  # args for every combine_fitres  job
        key_unique_list  = []  # define unique key per job for LOG and DONE files
        combine_outfile_list = [] 
        mimic_outdir_inp_list = []
        mimic_outdir_out_list = []
        symlink_src_list      = []
        symlink_dest_list     = []

        for task_dict in TASK_LIST:
            key = key = list(task_dict.keys())[0]
            #sys.exit(f"\n xxx task_dict = \n{task_dict}")
            INPUT_BASE      = task_dict[key]['INPUT_BASE']
            INPUT_APPEND    = task_dict[key].setdefault('INPUT_APPEND', [])
            OUTDIR_COMBINE  = task_dict[key]['OUTDIR_COMBINE']
            MIMIC_OUTDIR    = task_dict[key].setdefault('MIMIC_OUTDIR_SUBMIT_BATCH', None)
            if MIMIC_OUTDIR:
                MIMIC_OUTDIR_INP = MIMIC_OUTDIR.split()[0]
                MIMIC_OUTDIR_OUT = MIMIC_OUTDIR.split()[1]
            else:
                MIMIC_OUTDIR_INP = None
                MIMIC_OUTDIR_OUT = None

            if INPUT_TOPDIR :  
                INPUT_BASE     = f"{INPUT_TOPDIR}/{INPUT_BASE}"
                INPUT_APPEND   = [f"{INPUT_TOPDIR}/{x}" for x in INPUT_APPEND ]
                OUTDIR_COMBINE = f"{INPUT_TOPDIR}/{OUTDIR_COMBINE}"
                task_dict[key]['OUTDIR_COMBINE'] = OUTDIR_COMBINE
                if MIMIC_OUTDIR:
                    MIMIC_OUTDIR_INP = f"{INPUT_TOPDIR}/{MIMIC_OUTDIR_INP}"
                    MIMIC_OUTDIR_OUT = f"{INPUT_TOPDIR}/{MIMIC_OUTDIR_OUT}"

            # - - - - 
            inp_base_list = sorted(glob.glob(INPUT_BASE))
            n_inp_base    = len(inp_base_list) 
            n_inp_base_tot += n_inp_base
            if n_inp_base == 0 :
                logging.info(f" ERROR: cannot find INPUT_BASE = {INPUT_BASE}")
                nerr += 1


            inp_append_string = ' '.join(INPUT_APPEND)
            n_inp_append      = len(INPUT_APPEND)                
            n_inp_append_tot += n_inp_append
            for inp in INPUT_APPEND:
                if not os.path.exists(inp):
                    logging.info(f" ERROR: cannot find INPUT_APPEND = {inp}") 
                    nerr += 1
            # - - - -            
            # construct args for combine_fitres 
            for inp_base in inp_base_list:                
                base          = os.path.basename(inp_base)
                outfile       = f"{OUTDIR_COMBINE}/{base}"

                if n_inp_append > 0:
                    combine_args  = f"{inp_base} {inp_append_string} --outfile_text {outfile} -T" 
                    key_unique    = f"{key}_{base}".split('.gz')[0]  # remove .gz
                    key_unique    = key_unique.split('.')[0]
                    n_job_tot    += 1
                else :
                    combine_args  = None
                    key_unique    = None
                    symlink_src_list.append(inp_base)
                    symlink_dest_list.append(outfile)
                    n_symlink += 1

                combine_arg_list.append(combine_args)
                key_unique_list.append(key_unique)
                combine_outfile_list.append(outfile)
                mimic_outdir_inp_list.append(MIMIC_OUTDIR_INP)
                mimic_outdir_out_list.append(MIMIC_OUTDIR_OUT)

        # - - - -

        logging.info(f" Found {n_inp_base_tot:5d} INPUT_BASE files.")
        logging.info(f" Found {n_inp_append_tot:5d} INPUT_APPEND files.")
        logging.info(f" Prepared {n_job_tot} {PROGRAM_NAME_COMBINE_FITRES} jobs.")
        logging.info(f" Prepared {n_symlink} symbolic links tasks with nothing to append.")

        self.config_prep['n_job_tot']              = n_job_tot
        self.config_prep['n_done_tot']             = n_job_tot
        self.config_prep['n_job_split'] = 1

        self.config_prep['combine_arg_list']       = combine_arg_list
        self.config_prep['key_unique_list']        = key_unique_list
        self.config_prep['combine_outfile_list']   = combine_outfile_list
        self.config_prep['mimic_outdir_inp_list']  = mimic_outdir_inp_list
        self.config_prep['mimic_outdir_out_list']  = mimic_outdir_out_list
        self.config_prep['symlink_src_list']       = symlink_src_list
        self.config_prep['symlink_dest_list']      = symlink_dest_list

        if nerr > 0:
            msgerr = []
            msgerr.append(f"{nerr} errors checking inputs; see above messages")
            util.log_assert(False, msgerr)

        return
        # end prep_combine_check_inputs

    def prep_combine_outdirs(self):
        # loop over every output file and make sure that dirname is created
        combine_outfile_list  = self.config_prep['combine_outfile_list']

        logging.info(f"") 
        logging.info(f" Create outdirs for {PROGRAM_NAME_COMBINE_FITRES} output: ")

        outdir_checked_list = []
        for outfile in combine_outfile_list:
            outdir = os.path.dirname(outfile)
            if outdir not in outdir_checked_list:
                logging.info(f"   + {outdir}")
                if os.path.exists(outdir):  shutil.rmtree(outdir)
                desired_mode = 0o2770  # group rw
                os.makedirs(outdir, mode=desired_mode)
                outdir_checked_list.append(outdir)

        return
        # end prep_combine_outdirs

    def prep_combine_mimic_outdirs(self, stage):

        mimic_outdir_inp_list = self.config_prep['mimic_outdir_inp_list']
        mimic_outdir_out_list = self.config_prep['mimic_outdir_out_list']
        use_mimic_list = []

        logging.info(f"") 
        logging.info(f" Mimic submit_batch outdirs at {stage}-STAGE for {MIMIC_FILE_LIST}: ")

        for mimic_outdir_inp, mimic_outdir_out in zip(mimic_outdir_inp_list, mimic_outdir_out_list):
            USE_INP   = mimic_outdir_inp and mimic_outdir_inp not in use_mimic_list
            if USE_INP:
                logging.info(f"   + {mimic_outdir_out}")
                self.mimic_outdir_submit_batch(stage, mimic_outdir_inp, mimic_outdir_out)
                use_mimic_list.append(mimic_outdir_inp)
                                           
        return

    def prep_combine_symlinks(self):
        symlink_src_list  = self.config_prep['symlink_src_list']  
        symlink_dest_list = self.config_prep['symlink_dest_list'] 

        n_symlink = len(symlink_src_list)
        logging.info(f"\n Create {n_symlink} symbolic links (tasks with nothing to append)")
        use_dest_list = []
        for src, dest in zip(symlink_src_list,symlink_dest_list):
            if dest not in use_dest_list:   # avoid duplicate symlink requests
                os.symlink(src,dest)
                use_dest_list.append(dest)

        return

    def mimic_outdir_submit_batch(self, stage, mimic_outdir_inp, mimic_outdir_out):
        
        #logging.info(f" Copy {MIMIC_FILE_LIST} to {mimic_outdir_out}")

        for mimic_stage, mimic_file in zip(MIMIC_STAGE_LIST, MIMIC_FILE_LIST):

            if mimic_stage != stage: continue

            source_file = f"{mimic_outdir_inp}/{mimic_file}"
            dest_file   = f"{mimic_outdir_out}/{mimic_file}"
            if not os.path.exists(source_file):
                msgerr = [ f"Cound not find {source_file} to copy" ]
                util.log_assert(False, msgerr)
            
            if 'DONE' in source_file:
                shutil.copy(source_file, dest_file)  # current time stamp for ALL.DONE
            else:
                shutil.copy2(source_file, dest_file) # preserve original time stamp

        return

    def write_command_file(self, icpu, f):

        # For this icpu, write full set of sim commands to 
        # already-opened command file with pointer f.
        # Function returns number of jobs for this cpu
        # Nov 2024: create 20 blank list2d entries instead of 10
        
        n_job_tot        = self.config_prep['n_job_tot']
        n_core           = self.config_prep['n_core']
        combine_arg_list = self.config_prep['combine_arg_list']
        key_unique_list  = self.config_prep['key_unique_list']

        n_job_local    = 0     
        ii = 0
        for i, combine_args in enumerate(combine_arg_list):
            if combine_args is not None :  

                n_job_local += 1
                if ii % n_core == icpu:                
                    key_unique = key_unique_list[i]

                    job_info_combine   = self.prep_JOB_INFO_combine(combine_args, key_unique)
                    util.write_job_info(f, job_info_combine, icpu)

                    job_info_merge = self.prep_JOB_INFO_merge(icpu, n_job_local, False)
                    util.write_jobmerge_info(f, job_info_merge, icpu)

                ii += 1

        return
        # end write_command_file

    def prep_JOB_INFO_combine(self, combine_args, key_unique):

        program           = self.config_prep['program']
        output_dir        = self.config_prep['output_dir']

        args              = self.config_yaml['args']  # user args to submit_batch

        prefix    = f"COMBINE_{key_unique}"
        log_file  = f"{prefix}.LOG"
        done_file = f"{prefix}.DONE"
        yaml_file = f"{prefix}.YAML"

        JOB_INFO = {}
        JOB_INFO['job_dir']   = output_dir  # where to run job    
        JOB_INFO['program']   = program        
        JOB_INFO['input_file']    = ""

        arg_list = combine_args.split()
        arg_list.append(f"--outfile_yaml {yaml_file}")

        JOB_INFO['log_file']      = log_file
        JOB_INFO['done_file']     = done_file
        JOB_INFO['arg_list']      = arg_list
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['kill_on_fail']  = args.kill_on_fail
        JOB_INFO['check_abort']   = args.check_abort

        return JOB_INFO

        # end prep_JOB_INFO_combine



    # ========== MERGE STUFF =============

    def append_info_file(self,f):
        # append info to SUBMIT.INFO file  

        key_unique_list  = self.config_prep['key_unique_list']

        f.write(f"JOBFILE_WILDCARD:  'COMBINE*' \n")

        f.write(f"\n# {PROGRAM_NAME_COMBINE_FITRES} info\n")
        
        f.write(f"JOBNAME_LIST:\n")
        for key in key_unique_list:
            if key is not None:
                f.write(f"  - {key}\n")

        return
        # end append_info_file
                   

    def create_merge_table(self,f):

        key_unique_list  = self.config_prep['key_unique_list']

        header_line = " STATE    JOBNUM  JOBNAME   NFILE_COMBINE  NEVT_COMMON   CPU "

        MERGE_INFO = { 
            'primary_key' : TABLE_MERGE, 
            'header_line' : header_line,
            'row_list'      : []
        }
        
        jobnum = 0
        nff    = 0
        nevt   = 0
        CPU    = 0

        for jobname in key_unique_list:
            if jobname:
                ROW = [ SUBMIT_STATE_WAIT, jobnum, jobname, nff, nevt, CPU ]
                MERGE_INFO['row_list'].append(ROW)
                jobnum += 1

        util.write_merge_file(f, MERGE_INFO, [] )
        return
        # end create_merge_table

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.                  
        # Return update row list MERGE tables.     

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        key_unique_list  = self.config_prep['key_unique_list'] 

        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        #header_line = " STATE    JOBNUM  JOBNAME   NFILE_COMBINE  NEVT_COMMON   CPU "
        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_JOBNUM    = COLNUM_COMBINE_MERGE_JOBNUM
        COLNUM_JOBNAME   = COLNUM_COMBINE_MERGE_JOBNAME
        COLNUM_NFILE     = COLNUM_COMBINE_MERGE_NFILE
        COLNUM_NEVT      = COLNUM_COMBINE_MERGE_NEVT
        COLNUM_CPU       = COLNUM_COMBINE_MERGE_CPU

        key_nff, key_nff_sum, key_nff_list = \
                self.keynames_for_job_stats('NFILE_COMBINE')
        key_nevt, key_nevt_sum, key_nevt_list = \
                self.keynames_for_job_stats('NEVT_COMMON')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')

        key_list = [ key_nff, key_nevt, key_cpu ]
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
            STATE       = row[COLNUM_STATE]
            JOBNAME     = row[COLNUM_JOBNAME]
            search_wildcard = f"COMBINE_{JOBNAME}*"

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

                    row[COLNUM_STATE]     = NEW_STATE
                    row[COLNUM_NFILE]     = fit_stats[key_nff_sum]
                    row[COLNUM_NEVT]      = fit_stats[key_nevt_sum]
                    row[COLNUM_CPU]       = fit_stats[key_cpu_sum]

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
        key_unique_list  = self.config_prep['key_unique_list'] 
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_subdir    = SUBDIR_SCRIPTS_COSMOFIT

        logging.info(f"  combine cleanup for {JOB_SUFFIX_TAR_LIST}")
        for suffix in JOB_SUFFIX_TAR_LIST:
            wildcard = (f"COMBINE*.{suffix}")
            util.compress_files(+1, script_dir, wildcard, suffix, "" )

        util.compress_files(+1, script_dir, "CPU*", "CPU", "" )
        
        logging.info("")

        # do final MIMIC copy here so that ALL.DONE appears after all is merged
        self.prep_combine_mimic_outdirs(STAGE_DONE)

        # end merge_cleanup_final                        

    def merge_config_prep(self,output_dir):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        self.config_prep['key_unique_list'] = submit_info_yaml['JOBNAME_LIST']
        return

   
    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file        
        submit_info_yaml = self.config_prep['submit_info_yaml']
        info_lines = []
        return info_lines
        # end get_misc_merge_info     

    def get_merge_COLNUM_CPU(self):
        return COLNUM_COMBINE_MERGE_CPU
