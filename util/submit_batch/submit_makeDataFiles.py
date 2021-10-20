# Oct 20 2021
# MakeDataFiles.py batch code
# - Gautham Narayan and Rick Kessler (LSST DESC)


import  os, sys, shutil, yaml, configparser, glob
import  logging, coloredlogs
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *
from    submit_prog_base import Program


# Define columns in MERGE.LOG. Column 0 is always the STATE.
COLNUM_MKDATA_MERGE_DATAUNIT    = 1
COLNUM_MKDATA_MERGE_NEVENTS     = 2
COLNUM_MKDATA_MERGE_NSPECZ      = 3
COLNUM_MKDATA_MERGE_CPU         = 4
OUTPUT_FORMAT_LSST_ALERTS       = 'lsst_avro'
OUTPUT_FORMAT_SNANA             = 'snana'



# ====================================================
#    BEGIN FUNCTIONS
# ====================================================

class MakeDataFiles(Program):
    def __init__(self, config_yaml) :
        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_MKDATA
        super().__init__(config_yaml, config_prep)

    def set_output_dir_name(self):
        CONFIG     = self.config_yaml['CONFIG']
        input_file = self.config_yaml['args'].input_file  # for msgerr
        msgerr     = []
        if 'OUTDIR_ALERTS' in CONFIG:
            output_format = OUTPUT_FORMAT_LSST_ALERTS
            output_dir_name = os.path.expandvars(CONFIG['OUTDIR_ALERTS'])
        elif 'OUTDIR_SNANA' in CONFIG:
            output_format = OUTPUT_FORMAT_SNANA
            output_dir_name = os.path.expandvars(CONFIG['OUTDIR_SNANA'])
        else:
            msgerr.append(f"OUTDIR key missing in yaml-CONFIG")
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp
        self.config_yaml['args'].output_format = output_format
        return output_dir_name, SUBDIR_SCRIPTS_MKDATA
        # end set_output_dir_name

    def submit_prepare_driver(self):

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file

        # end submit_prepare_driver

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

        for itrain in range(0,n_trainopt):
            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :

                n_job_cpu += 1
                job_info_train   = self.prep_JOB_INFO_mkdata(itrain)
                util.write_job_info(f, job_info_train, icpu)

                job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_local)
                util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu

        # end write_command_file

    def prep_JOB_INFO_mkdata(self,itrain):

        CONFIG            = self.config_yaml['CONFIG']
        program           = self.config_prep['program']
        input_file_list   = self.config_prep['input_file_list']
        config_file       = self.config_prep['input_file_list'][0]
        script_dir        = self.config_prep['script_dir']
        kill_on_fail      = self.config_yaml['args'].kill_on_fail

        output_dir        = self.config_prep['output_dir']
        trainopt_num      = self.config_prep['trainopt_num_list'][itrain]
        trainopt_arg      = self.config_prep['trainopt_arg_list'][itrain]
        trainopt_global   = self.config_prep['trainopt_global']
        outdir_model      = self.config_prep['outdir_model_list'][itrain]
        do_fast           = self.config_yaml['args'].fast
        do_global         = len(trainopt_global) > 0

        prefix            = trainopt_num
        arg_list          = [ ]
        msgerr            = [ ]

        log_file   = f"{prefix}.LOG"
        done_file  = f"{prefix}.DONE"
        start_file = f"{prefix}.START"
        yaml_file  = f"{prefix}.YAML"

        # xxxx mark delete
        #for key,input_file in zip(CODE_KEYLIST_INPUT_FILE,input_file_list):
        #    arg_list.append(f"{key} {input_file}")
        # xxxx

        arg_list.append(f"--configfile {config_file}")
        arg_list.append(f"--outputdir {outdir_model}")
        arg_list.append(f"--yamloutputfile {yaml_file}")
        arg_list.append(f"{trainopt_arg}")
        if do_global : arg_list.append(trainopt_global)
        if do_fast   : arg_list.append("--fast")

        JOB_INFO = {}
        JOB_INFO['program']       = f"{program}"
        JOB_INFO['input_file']    = ""
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = log_file
        JOB_INFO['done_file']     = done_file
        JOB_INFO['start_file']    = start_file
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['arg_list']      = arg_list

        return JOB_INFO

        # end prep_JOB_INFO_mkdata

    def create_merge_table(self,f):

        n_trainopt        = self.config_prep['n_trainopt']
        trainopt_num_list = self.config_prep['trainopt_num_list']

        header_line_merge = (f"    STATE   {TRAINOPT_STRING}  NLC NSPEC  CPU")
        INFO_MERGE = {
            'primary_key' : TABLE_MERGE,
            'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for num in trainopt_num_list :
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(num)     # e..g, TRAINOPT002
            ROW_MERGE.append(0)       # NLC
            ROW_MERGE.append(0)       # NSPEC
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
        calib_shift_list = self.config_prep['calib_shift_list']
        outdir_model_list_base = self.config_prep['outdir_model_list_base']

        f.write(f"# train_SALT2 info \n")
        f.write(f"JOBFILE_WILDCARD: {TRAINOPT_STRING}* \n")

        f.write(f"\n")
        f.write(f"TRAINOPT_OUT_LIST:  " \
                f"# 'TRAINOPTNUM'  'user_label'  'user_args'\n")
        # use original ARG_list instead of arg_list; the latter may
        # include contents of shiftlist_file.
        for num, arg, label in zip(num_list, ARG_list, label_list):
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")

        f.write("MODELDIR_LIST:\n")
        for model_dir in outdir_model_list_base :
            f.write(f"  - {model_dir}\n")
        f.write("\n")

        for key in KEYS_SURVEY_LIST_SAME:
            if key in CONFIG :
                f.write(f"{key}:  {CONFIG[key]} \n")
            else:
                f.write(f"{key}:  [ ] \n")

        # write keys for SALT2.INFO to be read by SNANA code
        # each row is
        #  [ 'TRAINOPTnnn', KEY, SURVEY, SHIFT_VAL ]
        f.write(f"{KEY_SNANA_SALT3_INFO}: \n")
        for num, item_list in zip(num_list,calib_shift_list) :
            for item in item_list:
                row = [ num ] + item.split()
                f.write(f"  - {row} \n")
        f.write("\n")

        # end append_info_file

    def merge_config_prep(self,output_dir):
        pass

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_TRAINOPT  = COLNUM_TRAIN_MERGE_TRAINOPT
        COLNUM_NLC       = COLNUM_TRAIN_MERGE_NLC
        COLNUM_NSPEC     = COLNUM_TRAIN_MERGE_NSPEC
        COLNUM_CPU       = COLNUM_TRAIN_MERGE_CPU

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []
        row_list_merge     = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # keynames_for_job_stats returns 3 keynames :
        #   {base}, {base}_sum, {base}_list
        key_nlc, key_nlc_sum, key_nlc_list = \
                self.keynames_for_job_stats('NUM_SNE')
        key_nspec, key_nspec_sum, key_nspec_list = \
                 self.keynames_for_job_stats('NUM_SPECTRA')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')
        key_list = [ key_nlc, key_nspec, key_cpu ]

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

                job_stats = self.get_job_stats(script_dir,
                                               log_list, yaml_list, key_list)

                row[COLNUM_STATE]     = NEW_STATE
                row[COLNUM_NLC]       = job_stats[key_nlc_sum]
                row[COLNUM_NSPEC]     = job_stats[key_nspec_sum]
                row[COLNUM_CPU]       = job_stats[key_cpu_sum]

                row_list_merge_new[irow] = row  # update new row
                n_state_change += 1

        # first return arg (row_split) is null since there is
        # no need for a SPLIT table
        return [], row_list_merge_new, n_state_change
        # end merge_update_state

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):

        # cleanup for 'irow' training job.

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        modeldir_list    = submit_info_yaml['MODELDIR_LIST']

        row      = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        trainopt = row[COLNUM_TRAIN_MERGE_TRAINOPT]
        model_dir = f"{output_dir}/{modeldir_list[irow]}"
        misc_dir  = f"{model_dir}/{SUBDIR_MISC}"
        print(f" Cleanup {model_dir}")

        # - - - - - -
        # check if SALTshaker already cleaned up
        tar_file = f"{SUBDIR_MISC}.tar"
        tmp_list = glob.glob1(misc_dir,f"{tar_file}*")
        done_misc = len(tmp_list) > 0

        if not done_misc :
            os.mkdir(misc_dir)
            mv_string = "*.png *.pdf *.npy gauss* *parameters* " \
                        "salt3train_snparams.txt"
            cmd = f"cd {model_dir} ; mv {mv_string} {SUBDIR_MISC}/"
            os.system(cmd)

            # tar and gzip misc sub dir
            cmd_tar = f"cd {model_dir} ; " \
                    f"tar -cf {tar_file} {SUBDIR_MISC} ; " \
                    f"gzip  {tar_file} ; " \
                    f"rm -r {SUBDIR_MISC} "
            os.system(cmd_tar)

        # gzip dat files
        cmd_gzip = f"cd {model_dir} ; gzip *.dat "
        os.system(cmd_gzip)

        # end  merge_job_wrapup


    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file.
        # Each info line must be of the form
        #  KEYNAME:  VALUE

        return []
        # end get_misc_merge_info

    def merge_cleanup_final(self):
        # every snlc_fit job succeeded, so here we simply compress output.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        wildcard_list = [ 'TRAINOPT', 'CPU', 'CALIB_SHIFT' ]
        for w in wildcard_list :
            wstar = f"{w}*"
            tmp_list = glob.glob1(script_dir,wstar)
            if len(tmp_list) == 0 : continue
            print(f"\t Compress {wstar}")
            util.compress_files(+1, script_dir, wstar, w, "" )

        # - - - -
        # tar up entire script dir
        util.compress_subdir(+1, script_dir)

        # end merge_cleanup_final

    def get_merge_COLNUM_CPU(self):
        return COLNUM_TRAIN_MERGE_CPU

# =========== END: =======



