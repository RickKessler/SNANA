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
COLNUM_MKDATA_MERGE_NEVT        = 2
COLNUM_MKDATA_MERGE_NEVT_SPECZ  = 3
COLNUM_MKDATA_MERGE_NEVT_PHOTOZ = 4
OUTPUT_FORMAT_LSST_ALERTS       = 'lsst_avro'
OUTPUT_FORMAT_SNANA             = 'snana'
KEYLIST_OUTPUT                  = ['OUTDIR_SNANA', 'OUTDIR_ALERTS']
KEYLIST_OUTPUT_OPTIONS          = ['--outdir_snana', '--outdir_alerts']
BASE_PREFIX                     = 'MAKEDATA'
DATA_UNIT_STR                   = 'DATA_UNIT'


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

    def prepare_make_data_units(self):
        CONFIG      = self.config_yaml['CONFIG']
        inputs_list = CONFIG.get('MAKEDATAFILE_INPUTS', None)
        input_file  = self.config_yaml['args'].input_file  # for msgerr
        nsplitran   = CONFIG.get('NSPLITRAN', 1)
        field       = CONFIG.get('FIELD', None)
        split_mjd_detect_in = CONFIG.get('SPLIT_MJD_DETECT',None)
        msgerr      = []

        output_args  = None
        noutkeys = 0
        for key, opt in zip(KEYLIST_OUTPUT, KEYLIST_OUTPUT_OPTIONS):
            if key in CONFIG:
                outdir = CONFIG[key]
                if '/' not in outdir: # checking to make sure that the outdir has a full path
                    outdir = os.getcwd() + '/' + outdir
                output_args = f'{opt} {outdir}'
                noutkeys += 1

        if noutkeys > 1:
            msgerr.append(f"Multiply defined key for output format in yaml-CONFIG")
            msgerr.append(f'Require EXACTLY one of {KEYLIST_OUTPUT}')
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp

        if output_args is None:
            msgerr.append(f"Missing key for output format in yaml-CONFIG")
            msgerr.append(f'Require one of {KEYLIST_OUTPUT}')
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp

        if inputs_list is None:
            msgerr.append(f"MAKEDATAFILE_INPUTS key missing in yaml-CONFIG")
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp

        split_mjd_detect = {}
        if split_mjd_detect_in is None:
            split_mjd_detect['nbin'] = 1
            split_mjd_detect['step'] = 0
        else:
            mjdmin, mjdmax, mjdbin  = split_mjd_detect_in.split()
            split_mjd_detect['min']  = int(mjdmin)
            split_mjd_detect['max']  = int(mjdmax)
            split_mjd_detect['nbin'] = int(mjdbin)
            split_mjd_detect['step'] = (split_mjd_detect['max'] - split_mjd_detect['min'])/split_mjd_detect['nbin']

        makeDataFiles_args_list = []
        prefix_output_list = []
        isplitmjd = 0 ### HACK HACK HACK - fix this when we add a look over MJD
        for input_src in inputs_list:  # e.g. folder or name of DB
            for isplitran in range(nsplitran):
                # add a third index for MJD here
                base_name = os.path.basename(input_src)
                prefix_output = self.get_prefix_output(isplitmjd, base_name, isplitran)
                prefix_output_list.append(prefix_output)

                args = f'--snana_folder {input_src} ' ### HACK HACK HACK - will need to generalize for other inputs
                args += f'{output_args} '
                args += f'--field {field} '
                if nsplitran > 1:
                    args += f'--nsplitran {nsplitran} '
                    args += f'--isplitran {isplitran+1} ' # note that argument for isplitran starts with 1
                makeDataFiles_args_list.append(args)

        self.config_prep['inputs_list'] = inputs_list
        self.config_prep['split_mjd_detect'] = split_mjd_detect
        self.config_prep['makeDataFiles_args_list'] = makeDataFiles_args_list
        self.config_prep['prefix_output_list'] = prefix_output_list
        # end prepare_make_data_units

    def get_prefix_output(self, isplitmjd, base_name, isplitran):

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file
        output_format = self.config_yaml['args'].output_format
        msgerr = []
        splitran_str = f'SPLITRAN{isplitran+1:03d}'
        if output_format == OUTPUT_FORMAT_LSST_ALERTS:
            # isplit
            prefix_output = f'{BASE_PREFIX}_SPLITMJD{isplitmjd:03d}_{base_name}_{splitran_str}'
        elif output_format == OUTPUT_FORMAT_SNANA:
            prefix_output = f'{BASE_PREFIX}_{base_name}_{splitran_str}'
        else:
            msgerr.append(f'Invalid Output Format {output_format}')
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp
        return prefix_output


    def submit_prepare_driver(self):

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file
        self.prepare_make_data_units()

        # end submit_prepare_driver

    def write_command_file(self, icpu, f):
        # For this icpu, write full set of sim commands to
        # already-opened command file with pointer f.
        # Function returns number of jobs for this cpu

        n_core          = self.config_prep['n_core']
        makeDataFiles_args_list = self.config_prep['makeDataFiles_args_list']
        n_job_tot   = len(makeDataFiles_args_list)
        n_job_split = 1     # cannot break up makeDataFiles job since already broken up
        n_job_local = 0
        n_job_cpu   = 0

        self.config_prep['n_job_split'] = n_job_split
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        for idata_unit in range(0,n_job_tot):
            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :

                n_job_cpu += 1
                job_info_data_unit   = self.prep_JOB_INFO_mkdata(idata_unit)
                util.write_job_info(f, job_info_data_unit, icpu)

                job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_local)
                util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu

        # end write_command_file

    def prep_JOB_INFO_mkdata(self,idata_unit):

        CONFIG            = self.config_yaml['CONFIG']
        makeDataFiles_arg = self.config_prep['makeDataFiles_args_list'][idata_unit]
        prefix            = self.config_prep['prefix_output_list'][idata_unit]
        program           = self.config_prep['program']
        script_dir        = self.config_prep['script_dir']
        kill_on_fail      = self.config_yaml['args'].kill_on_fail

        output_dir        = self.config_prep['output_dir']
        # do_fast           = self.config_yaml['args'].fast

        arg_list          = [makeDataFiles_arg,]
        msgerr            = [ ]

        log_file   = f"{prefix}.LOG"
        done_file  = f"{prefix}.DONE"
        start_file = f"{prefix}.START"
        yaml_file  = f"{prefix}.YAML"

        arg_list.append(f"--output_yaml_file {yaml_file}")
        # if do_fast   : arg_list.append("--fast")        # may need later

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

        prefix_output_list        = self.config_prep['prefix_output_list']

        header_line_merge = (f"    STATE   {DATA_UNIT_STR}  NEVT NEVT_SPECZ NEVT_PHOTOZ")
        INFO_MERGE = {
            'primary_key' : TABLE_MERGE,
            'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for prefix in prefix_output_list :
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(prefix)
            ROW_MERGE.append(0)       # NEVT
            ROW_MERGE.append(0)       # NEVT_SPECZ
            ROW_MERGE.append(0)       # NEVT_PHOTOZ

            INFO_MERGE['row_list'].append(ROW_MERGE)
        # - - - -
        util.write_merge_file(f, INFO_MERGE, [] )

        # end create_merge_table

    def append_info_file(self,f):

        # append info to SUBMIT.INFO file; use passed file pointer f

        CONFIG       = self.config_yaml['CONFIG']
        prefix_output_list   = self.config_prep['prefix_output_list']

        f.write(f"# makeDataFiles info \n")
        f.write(f"JOBFILE_WILDCARD: {BASE_PREFIX}* \n")

        f.write(f"\n")
        f.write(f"PREFIX_OUTPUT_LIST:  \n" )
        # use original ARG_list instead of arg_list; the latter may
        # include contents of shiftlist_file.
        for prefix in prefix_output_list:
            f.write(f"  - {prefix} \n")
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

        COLNUM_STATE       = COLNUM_MERGE_STATE
        COLNUM_DATAUNIT    = COLNUM_MKDATA_MERGE_DATAUNIT
        COLNUM_NEVT        = COLNUM_MKDATA_MERGE_NEVT
        COLNUM_NEVT_SPECZ  = COLNUM_MKDATA_MERGE_NEVT_SPECZ
        COLNUM_NEVT_PHOTOZ = COLNUM_MKDATA_MERGE_NEVT_PHOTOZ

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []
        row_list_merge     = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # keys from YAML file
        # NEVT_ALL:                118
        # NEVT_HOSTGAL_SPECZ:      60
        # NEVT_HOSTGAL_PHOTOZ:     118
        # NEVT_SPECTRA:            0

        # keynames_for_job_stats returns 3 keynames :
        #   {base}, {base}_sum, {base}_list
        key_nall, key_nall_sum, key_nall_list = \
                self.keynames_for_job_stats('NEVT_ALL')
        key_nspecz, key_nspecz_sum, key_nspecz_list = \
                 self.keynames_for_job_stats('NEVT_HOSTGAL_SPECZ')
        key_nphotz, key_nphotz_sum, key_nphotz_list = \
                 self.keynames_for_job_stats('NEVT_HOSTGAL_PHOTOZ')
        #key_cpu, key_cpu_sum, key_cpu_list = \
        #        self.keynames_for_job_stats('CPU_MINUTES')
        key_list = [ key_nall, key_nspecz, key_nphotz ]

        nrow_check = 0
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input
            nrow_check += 1
            irow        = nrow_check - 1 # row index
            data_unit    = row[COLNUM_DATAUNIT]
            search_wildcard = (f"{data_unit}*")

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

                if NLOG > 0 :
                    NEW_STATE = SUBMIT_STATE_RUN
                if NDONE == n_job_split :
                    NEW_STATE = SUBMIT_STATE_DONE

                    job_stats = self.get_job_stats(script_dir,
                                                   log_list, 
                                                   yaml_list, 
                                                   key_list)

                    row[COLNUM_STATE]       = NEW_STATE
                    row[COLNUM_NEVT]        = job_stats[key_nall_sum]
                    row[COLNUM_NEVT_SPECZ]  = job_stats[key_nspecz_sum]
                    row[COLNUM_NEVT_PHOTOZ] = job_stats[key_nphotz_sum]

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

        row      = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        print(f" Do nothing!")

        # Gautham - 20211026 - for alerts - will need to edit this and create gz files
        # gzip dat files
        # cmd_gzip = f"cd {model_dir} ; gzip *.dat "
        # os.system(cmd_gzip)

        # end  merge_job_wrapup


    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file.
        # Each info line must be of the form
        #  KEYNAME:  VALUE

        return []
        # end get_misc_merge_info

    def merge_cleanup_final(self):
        # every makeDataFiles succeeded, so here we simply compress output.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        cwd              = submit_info_yaml['CWD']

        command_list = ['makeDataFiles.sh', '--outdir_snana', output_dir, '--merge']
        ret = subprocess.run(  command_list,
                                 cwd=cwd,
                                 capture_output=False, text=True )
        print(ret,' debugging')

        wildcard_list = [ 'MAKEDATA', 'CPU',  ]
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
        return -9

# =========== END: =======



