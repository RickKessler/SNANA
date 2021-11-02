# Oct 20 2021
# MakeDataFiles.py batch code
# - Gautham Narayan and Rick Kessler (LSST DESC)


import  os, sys, shutil, yaml, configparser, glob
import  logging, coloredlogs, tarfile
import  datetime, time, subprocess
import  submit_util as util
from    submit_params    import *
from    submit_prog_base import Program
import numpy as np


# Define columns in MERGE.LOG. Column 0 is always the STATE.
COLNUM_MKDATA_MERGE_ISPLITMJD   = 1
COLNUM_MKDATA_MERGE_DATAUNIT    = 2
COLNUM_MKDATA_MERGE_NEVT        = 3
COLNUM_MKDATA_MERGE_NEVT_SPECZ  = 4
COLNUM_MKDATA_MERGE_NEVT_PHOTOZ = 5
COLNUM_MKDATA_MERGE_NOBS_ALERT  = 6  # for lsst alerts only
COLNUM_MKDATA_MERGE_RATE        = 6  # or add 1 for lsst_alerts

OUTPUT_FORMAT_LSST_ALERTS       = 'lsst_avro'
OUTPUT_FORMAT_SNANA             = 'snana'

KEYLIST_OUTPUT           = ['OUTDIR_SNANA',   'OUTDIR_LSST_ALERT']
KEYLIST_OUTPUT_OPTIONS   = ['--outdir_snana', '--outdir_lsst_alert']
OUTPUT_FORMAT            = [OUTPUT_FORMAT_SNANA, OUTPUT_FORMAT_LSST_ALERTS]

KEYLIST_SPLIT_MJD          = ['SPLIT_MJD_DETECT', 'SPLIT_PEAKMJD']
KEYLIST_SPLIT_MJD_OPTIONS  = ['--mjd_detect_range', '--peakmjd_range']

BASE_PREFIX          = 'MAKEDATA'   # base for log,yaml,done files
DATA_UNIT_STR        = 'DATA_UNIT'  # merge table comment

# params for lsst alerts
ALERT_SUBDIR        = "ALERTS"     # move mjd tar files here
ALERT_DAY_NAME      = "mjd"      # xxx later switch to nite 
TABLE_COMPRESS      = "COMPRESS" # name of extra table in MERGE.LOG file
COLNUM_COMPRESS_ISPLITMJD = 1  # index for split MJD range
COLNUM_COMPRESS_MJDRANGE  = 2
COLNUM_COMPRESS_NMJD_DIR  = 3  # number of day directories
COLNUM_COMPRESS_TIME      = 4  # Nsec
COLNUM_COMPRESS_RATE      = 5  # Ndir/sec


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
        output_format = None
        
        for key_output, format_str in zip(KEYLIST_OUTPUT, OUTPUT_FORMAT):
            
            if key_output in CONFIG:
                output_format = format_str
                output_dir_name = os.path.expandvars(CONFIG[key_output])

        if output_format is None:
            msgerr.append(f"OUTDIR key missing in yaml-CONFIG")
            msgerr.append(f"Must provide one of {KEYLIST_OUTPUT}")
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp
            
        self.config_yaml['args'].output_format = output_format

        
        return output_dir_name, SUBDIR_SCRIPTS_MKDATA
        # end set_output_dir_name

    def prepare_output_args(self):
        '''
        Prepare the output directory based on format options (SNANA/ALERTS/etc)
        '''
        CONFIG      = self.config_yaml['CONFIG']
        input_file  = self.config_yaml['args'].input_file  # for msgerr
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

        self.config_prep['output_args'] = output_args
        # end prepare_output_args


    def prepare_input_args(self):
        '''
        Prepare input arguments from config file
        '''
        CONFIG        = self.config_yaml['CONFIG']
        inputs_list   = CONFIG.get('MAKEDATAFILE_INPUTS', None)
        input_source  = CONFIG.get('MAKEDATAFILE_SOURCE', None)
        nevt          = CONFIG.get('NEVT', None)
                                   
        input_file    = self.config_yaml['args'].input_file  # for msgerr
        msgerr        = []

        if inputs_list is None:
            msgerr.append(f"MAKEDATAFILE_INPUTS key missing in yaml-CONFIG")
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp

        # select the SPLIT_MJD option
        # abort if more than one SPLIT_MJD option is specified
        n_mjd_split_opts   = 0
        split_mjd_in       = None
        split_mjd_key_name = None
        split_mjd_option   = None
        for key, opt in zip(KEYLIST_SPLIT_MJD, KEYLIST_SPLIT_MJD_OPTIONS):
            if key in CONFIG:
                n_mjd_split_opts  += 1
                split_mjd_key_name = key
                split_mjd_option   = opt
                split_mjd_in       = CONFIG[key]
        if n_mjd_split_opts > 1:
            msgerr.append(f"DEFINE ONLY ONE OF {KEYLIST_SPLIT_MJD}")
            msgerr.append(f"Check {input_file}")
            util.log_assert(False,msgerr) # just abort, no done stamp

        # parse the input SPLIT_MJD string into ranges for makeDataFiles
        split_mjd = {}
        if split_mjd_in is None:
            split_mjd['nbin'] = 1
            split_mjd['step'] = 0
        else:
            mjdmin, mjdmax, mjdbin  = split_mjd_in.split()
            imjdmin = int(mjdmin)
            imjdmax = int(mjdmax)
            nbin    = int(mjdbin)
            split_mjd['min']  = imjdmin
            split_mjd['max']  = imjdmax
            split_mjd['nbin'] = nbin
            split_mjd['step'] = (imjdmax - imjdmin) / nbin
            # use nbin + 1 to include the edges
            grid  = np.linspace(imjdmin, imjdmax, nbin+1)
            split_mjd['min_edge'] = grid[0:-1]
            split_mjd['max_edge'] = grid[1:]

        self.config_prep['inputs_list'] = inputs_list
        self.config_prep['split_mjd']   = split_mjd
        self.config_prep['split_mjd_key_name'] = split_mjd_key_name  # CONFIG YAML keyname
        self.config_prep['split_mjd_option'] = split_mjd_option #makeDataFiles.sh option

        self.config_prep['input_source'] = input_source
        self.config_prep['nevt']         = nevt

        # end prepare_input_args

    def prepare_data_units(self):
        CONFIG      = self.config_yaml['CONFIG']
        output_args = self.config_prep['output_args']
        inputs_list = self.config_prep['inputs_list']
        input_file  = self.config_yaml['args'].input_file  # for msgerr
        split_mjd   = self.config_prep['split_mjd']
        split_mjd_option = self.config_prep['split_mjd_option']
        n_splitmjd  = split_mjd['nbin']
        n_splitran  = CONFIG.get('NSPLITRAN', 1)
        field       = CONFIG.get('FIELD', None)
        msgerr      = []

        makeDataFiles_args_list = []
        prefix_output_list = []
        isplitmjd_list     = []

        if n_splitmjd > 1:
            isplitmjd_temp_list = zip(range(0, n_splitmjd),\
                                            split_mjd['min_edge'],\
                                            split_mjd['max_edge'])

        else:
            isplitmjd_temp_list = [(-1, -1, -1),]
        # When using zip in Python 3, you get an iterator, not a list,
        # and thus if you print the iterator or cause anything to loop 
        # over it e.g. len(), then the generator is exhausted and has 
        # to be recreated
        # i.e. do not mess with any variable that is created from a zip
        # DO NOT MESS WITH isplitmjd_temp_list

        for isplitmjd, min_edge, max_edge in isplitmjd_temp_list:
            for input_src in inputs_list:  # e.g. folder or name of DB
                args_list = []
                base_name = os.path.basename(input_src)

                # construct base prefix without isplitran
                prefix_output = self.get_prefix_output(min_edge, base_name, -1)
                prefix_output_list.append(prefix_output)
                isplitmjd_list.append(isplitmjd)

                args_list.append(f'--snana_folder {input_src}') ### HACK need to generalize for other inputs

                args_list.append(f'{output_args}')
                args_list.append(f'--field {field}')
                if n_splitmjd > 1:
                    args_list.append(f'{split_mjd_option} {min_edge} {max_edge}')
                if n_splitran > 1:
                    args_list.append(f'--nsplitran {n_splitran}')
                    # args_list.append(f'--isplitran {isplitran+1}') # note that argument for isplitran starts with 1
                makeDataFiles_args_list.append(args_list)

        n_job       = len(makeDataFiles_args_list)
        n_job_tot   = n_job*n_splitran
        n_job_split = n_splitran
        n_job_local = 0
        n_job_cpu   = 0

        idata_unit_list = []
        isplitran_list  = []

        for idata_unit in range(0, n_job):
            for isplitran in range(0, n_splitran):
                idata_unit_list.append(idata_unit)
                isplitran_list.append(isplitran)

        self.config_prep['n_job']       = n_job
        self.config_prep['n_job_split'] = n_job_split
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot
        self.config_prep['idata_unit_list'] = idata_unit_list
        self.config_prep['isplitran_list']  = isplitran_list

        self.config_prep['makeDataFiles_args_list'] = makeDataFiles_args_list
        self.config_prep['prefix_output_list'] = prefix_output_list
        self.config_prep['isplitmjd_list']     = isplitmjd_list

        # end prepare_data_units


    def get_prefix_output(self, mjd, base_name, isplitran):

        CONFIG        = self.config_yaml['CONFIG']
        input_file    = self.config_yaml['args'].input_file
        output_format = self.config_yaml['args'].output_format
        msgerr = []
        imjd = int(mjd)
        splitran_str  = f'SPLITRAN{isplitran+1:03d}'
        splitmjd_str  = f'SPLITMJD{imjd:05d}'
        prefix_output = f'{BASE_PREFIX}'

        if imjd >= 10000:
            prefix_output += f'_{splitmjd_str}'

        prefix_output += f'_{base_name}'

        if isplitran >= 0:
            prefix_output += f'_{splitran_str}'
        return prefix_output
        # end get_prefix_output
        
    def submit_prepare_driver(self):

        # called from base to prepare makeDataFile arguments for batch.

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file
        script_dir   = self.config_prep['script_dir']
        output_dir   = self.config_prep['output_dir']

        self.prepare_output_args()
        self.prepare_input_args()
        self.prepare_data_units()
        
        # copy input config file to script-dir
        shutil.copy(input_file,script_dir)

        # create ALERTS subdir for final mjd-tar files
        alerts_dir    = f"{output_dir}/{ALERT_SUBDIR}"
        self.config_prep['alerts_dir'] = alerts_dir
        os.mkdir(alerts_dir)

        # end submit_prepare_driver

    def write_command_file(self, icpu, f):

        # Called from base;
        # For this icpu, write full set of sim commands to
        # already-opened command file with pointer f.
        # Function returns number of jobs for this cpu

        n_core          = self.config_prep['n_core']
        CONFIG          = self.config_yaml['CONFIG']
        makeDataFiles_args_list = self.config_prep['makeDataFiles_args_list']

        n_job       = self.config_prep['n_job']
        n_job_split = self.config_prep['n_job_split']
        n_job_tot   = self.config_prep['n_job_tot']
        n_job_tot   = self.config_prep['n_done_tot']
        idata_unit_list = self.config_prep['idata_unit_list']
        isplitran_list  = self.config_prep['isplitran_list']
        n_job_local     = 0
        n_job_cpu       = 0

        index_dict = {}
        for idata_unit, isplitran in zip(idata_unit_list, isplitran_list):

            n_job_local += 1
            if ( (n_job_local-1) % n_core ) == icpu :

                n_job_cpu += 1
                index_dict['isplitran']   = isplitran
                index_dict['icpu']        = icpu
                index_dict['idata_unit']  = idata_unit
                index_dict['n_job_local'] = n_job_local

                job_info_data_unit   = self.prep_JOB_INFO_mkdata(index_dict)
                util.write_job_info(f, job_info_data_unit, icpu)

                job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_local)
                util.write_jobmerge_info(f, job_info_merge, icpu)

        #print(f" xxx n_core={n_core}   n_job_cpu = {n_job_cpu}  " \
        #      f" n_job_local = {n_job_local}")
        
        return n_job_cpu
        # end write_command_file


    def prep_JOB_INFO_mkdata(self,index_dict):

        CONFIG  = self.config_yaml['CONFIG']

        isplitran   = index_dict['isplitran']
        isplitarg   = isplitran + 1           # passed to makeDataFiles.py
        icpu        = index_dict['icpu']
        idata_unit  = index_dict['idata_unit']
        n_job_local = index_dict['n_job_local']

        makeDataFiles_arg = self.config_prep['makeDataFiles_args_list'][idata_unit]
        isplitmjd         = self.config_prep['isplitmjd_list'][idata_unit]
        prefix_base       = self.config_prep['prefix_output_list'][idata_unit]
        prefix            = f'{prefix_base}_SPLITRAN{isplitarg:03d}'
        program           = self.config_prep['program']
        script_dir        = self.config_prep['script_dir']
        output_dir        = self.config_prep['output_dir']
        nevt              = self.config_prep['nevt']

        kill_on_fail      = self.config_yaml['args'].kill_on_fail
        output_format     = self.config_yaml['args'].output_format
        no_merge          = self.config_yaml['args'].nomerge

        out_lsst_alert    = (output_format == OUTPUT_FORMAT_LSST_ALERTS)
        # do_fast           = self.config_yaml['args'].fast

        msgerr            = [ ]
        log_file   = f"{prefix}.LOG"
        done_file  = f"{prefix}.DONE"
        start_file = f"{prefix}.START"
        yaml_file  = f"{prefix}.YAML"

        arg_split         = f'--isplitran {isplitarg}'
        arg_list          = makeDataFiles_arg + [arg_split,]
        
        if out_lsst_alert :
            schema_file = CONFIG['LSST_ALERT_SCHEMA']
            arg_list.append(f"--lsst_alert_schema {schema_file}")

        if nevt is not None:
            arg_list.append(f"--nevt {nevt}")

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

        # for lsst alerts, wait for previous compress-ISPLITMJD to finish
        # to avoid piling up too many alert files. Remember that 
        # isplitmjd goes from 1 - nsplitmjd.
        set_wait_file = out_lsst_alert and (isplitmjd > 0) \
                        and (not no_merge)
        if set_wait_file :
            split_mjd   = self.config_prep['split_mjd']
            isplit_previous = isplitmjd-1 
            min_edge    = split_mjd['min_edge'][isplit_previous]
            max_edge    = split_mjd['max_edge'][isplit_previous]
            wait_file   = self.get_compress_done_file(min_edge,max_edge)
            JOB_INFO['wait_file'] = wait_file

        return JOB_INFO
        # end prep_JOB_INFO_mkdata

    def create_merge_table(self,f):

        # Called from base to create rows for table in  MERGE.LOG 
        # Always create required MERGE table.
        # For LSST alerts, also create supplemental "COMPRESS" table
        # to track mjd compression
        
        isplitmjd_list      = self.config_prep['isplitmjd_list']
        prefix_output_list  = self.config_prep['prefix_output_list']
        output_format       = self.config_yaml['args'].output_format 
        out_lsst_alert      = (output_format == OUTPUT_FORMAT_LSST_ALERTS)

        # 1. required MERGE table
        header_line_merge = f"    STATE   ISPLITMJD  {DATA_UNIT_STR}  " \
                            f"NEVT NEVT_SPECZ NEVT_PHOTOZ  "
        if out_lsst_alert : 
            header_line_merge += "NOBS_ALERT  ALERT/sec"
        else:
            header_line_merge += "NEVT/sec"
        
        INFO_MERGE = {
            'primary_key' : TABLE_MERGE,
            'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for prefix, isplitmjd in zip(prefix_output_list,isplitmjd_list):
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(isplitmjd) 
            ROW_MERGE.append(prefix)    # data unit name
            ROW_MERGE.append(0)         # NEVT
            ROW_MERGE.append(0)         # NEVT_SPECZ
            ROW_MERGE.append(0)         # NEVT_PHOTOZ
            if out_lsst_alert: ROW_MERGE.append(0)  # NOBS_ALERT
            ROW_MERGE.append(0.0)       # rate/sec
            
            INFO_MERGE['row_list'].append(ROW_MERGE)

        # call util to write tables to MERGE.LOG
        if out_lsst_alert: 
            self.create_compress_table(f)

        util.write_merge_file(f, INFO_MERGE, [] )
        
        # end create_merge_table

    def create_compress_table(self, f):
        # Called for lsst alert output

        split_mjd            = self.config_prep['split_mjd']
        nsplitmjd            = split_mjd['nbin']
        min_edge_list        = split_mjd['min_edge']
        max_edge_list        = split_mjd['max_edge']
        header_line_compress = \
            f"    STATE   ISPLITMJD MJD-RANGE  NDIR_{ALERT_DAY_NAME}  " \
            f"Nsec NDIR/sec"
                        
        INFO_COMPRESS = {
            'primary_key' : TABLE_COMPRESS,
            'header_line' : header_line_compress,
            'row_list'    : []   }
        
        STATE = SUBMIT_STATE_WAIT    # all start in WAIT state
        for isplitmjd in range(0,nsplitmjd):
            imin         = int(min_edge_list[isplitmjd])
            imax         = int(max_edge_list[isplitmjd])
            str_mjd_range = f"{imin}-{imax}"
                
            ROW_COMPRESS = []
            ROW_COMPRESS.append(STATE)
            ROW_COMPRESS.append(isplitmjd)      # index: 0,1,...
            ROW_COMPRESS.append(str_mjd_range)  # e.g., 59000-59200
            ROW_COMPRESS.append(0)              # init NDIR_MJD=0
            ROW_COMPRESS.append(0)              # init Nsec
            ROW_COMPRESS.append(0.0)            # init rate = NDIR/sec
            
            INFO_COMPRESS['row_list'].append(ROW_COMPRESS)
        util.write_merge_file(f, INFO_COMPRESS, [] )

        # end create_compress_table

    def append_info_file(self,f):

        # Called from base to
        # append info to SUBMIT.INFO file; use passed file pointer f

        CONFIG              = self.config_yaml['CONFIG']       
        output_format       = self.config_yaml['args'].output_format
        prefix_output_list  = self.config_prep['prefix_output_list']
        input_source        = self.config_prep['input_source']
        alerts_dir          = self.config_prep['alerts_dir']
        split_mjd_key_name  = self.config_prep['split_mjd_key_name']
        split_mjd           = self.config_prep['split_mjd']
        nsplitmjd           = split_mjd['nbin']

        f.write(f"# makeDataFiles info \n")
        f.write(f"JOBFILE_WILDCARD: {BASE_PREFIX}* \n")
        f.write(f"\n")

        f.write(f"MAKEDATAFILE_SOURCE: {input_source} \n")
        f.write(f"OUTPUT_FORMAT:   {output_format} \n")
        f.write(f"ALERTS_DIR:      {alerts_dir}\n")
        f.write(f"\n")

        f.write(f"KEYNAME_SPLITMJD:  {split_mjd_key_name}\n")
        f.write(f"NSPLITMJD: {nsplitmjd} \n");
        if nsplitmjd > 1:
            min_edge = list(split_mjd['min_edge'])
            max_edge = list(split_mjd['max_edge'])
            f.write(f"MIN_MJD_EDGE: {min_edge} \n")
            f.write(f"MAX_MJD_EDGE: {max_edge} \n")
        f.write(f"\n")

        # write out each job prefix
        f.write(f"PREFIX_OUTPUT_LIST:  \n" )
        for prefix in prefix_output_list:
            f.write(f"  - {prefix} \n")
        f.write("\n")
        # end append_info_file

    def merge_config_prep(self,output_dir):
        pass

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # Called from base to
        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.
        # For lsst alerts, also update COMPRESS table.
        
        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        output_format       = submit_info_yaml['OUTPUT_FORMAT']
        out_lsst_alert      = (output_format == OUTPUT_FORMAT_LSST_ALERTS)
        
        COLNUM_STATE       = COLNUM_MERGE_STATE
        COLNUM_DATAUNIT    = COLNUM_MKDATA_MERGE_DATAUNIT
        COLNUM_NEVT        = COLNUM_MKDATA_MERGE_NEVT
        COLNUM_NEVT_SPECZ  = COLNUM_MKDATA_MERGE_NEVT_SPECZ
        COLNUM_NEVT_PHOTOZ = COLNUM_MKDATA_MERGE_NEVT_PHOTOZ
        COLNUM_NOBS_ALERT  = COLNUM_MKDATA_MERGE_NOBS_ALERT
        COLNUM_RATE        = COLNUM_MKDATA_MERGE_RATE
        if out_lsst_alert: COLNUM_RATE += 1

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []
        row_list_merge     = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # keynames_for_job_stats returns 3 keynames :
        #   {base}, {base}_sum, {base}_list
        key_nall, key_nall_sum, key_nall_list = \
                self.keynames_for_job_stats('NEVT_ALL')
        key_nspecz, key_nspecz_sum, key_nspecz_list = \
                 self.keynames_for_job_stats('NEVT_HOSTGAL_SPECZ')
        key_nphotz, key_nphotz_sum, key_nphotz_list = \
                 self.keynames_for_job_stats('NEVT_HOSTGAL_PHOTOZ')        
        key_tproc, key_tproc_sum, key_tproc_list = \
                 self.keynames_for_job_stats('WALLTIME')  

        key_list = [ key_nall, key_nspecz, key_nphotz, key_tproc ]

        if out_lsst_alert:
            key_nalert, key_nalert_sum, key_nalert_list = \
                self.keynames_for_job_stats('NOBS_ALERT')
            key_list += [ key_nalert ]
                                                          
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
                                                   key_list )

                    row[COLNUM_STATE]       = NEW_STATE
                    row[COLNUM_NEVT]        = job_stats[key_nall_sum]
                    row[COLNUM_NEVT_SPECZ]  = job_stats[key_nspecz_sum]
                    row[COLNUM_NEVT_PHOTOZ] = job_stats[key_nphotz_sum]
                    
                    if out_lsst_alert:
                        row[COLNUM_NOBS_ALERT] = job_stats[key_nalert_sum]
                        n_tmp = row[COLNUM_NOBS_ALERT]
                    else:
                        n_tmp = row[COLNUM_NEVT]

                    # load N/sec instead of CPU time
                    t_proc = job_stats[key_tproc_sum]
                    rate   = 0.0
                    if t_proc > 0.0 : rate   = n_tmp / t_proc
                    row[COLNUM_RATE] = float(f"{rate:.1f}")

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1

        # - - - - - - - - - - -
        # check for optional extra table
        row_extra_list = []
        if out_lsst_alert:
            row_extra_list = self.compress_update_state(MERGE_INFO_CONTENTS)
            
        # first return arg (row_split) is null since there is
        # no need for a SPLIT table
        
        row_list_dict = {
            'row_split_list'   : [],
            'row_merge_list'   : row_list_merge_new,
            'row_extra_list'   : row_extra_list,
            'table_names'      : [ TABLE_SPLIT, TABLE_MERGE, 
                                   TABLE_COMPRESS ]
        }
        return row_list_dict, n_state_change
    
        # end merge_update_state

    def compress_update_state(self,MERGE_INFO_CONTENTS):
        
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        nsplitmjd        = submit_info_yaml['NSPLITMJD']
        
        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_ISPLITMJD = COLNUM_COMPRESS_ISPLITMJD
        COLNUM_NMJD_DIR  = COLNUM_COMPRESS_NMJD_DIR
        COLNUM_TIME      = COLNUM_COMPRESS_TIME  # Nsec
        COLNUM_RATE      = COLNUM_COMPRESS_RATE  # Ndir/sec
        
        row_merge_list        = MERGE_INFO_CONTENTS[TABLE_MERGE]
        row_compress_list     = MERGE_INFO_CONTENTS[TABLE_COMPRESS]
        row_compress_list_new = [] 
        nrow = 0
        
        for row in row_compress_list:
            # strip off row info
            row_compress_list_new.append(row)
            nrow += 1
            STATE       = row[COLNUM_STATE]
            ISPLITMJD   = row[COLNUM_ISPLITMJD]
            # check if DONE or FAIL ; i.e., if Finished
            Finished = (STATE==SUBMIT_STATE_DONE) or (STATE==SUBMIT_STATE_FAIL)

            if Finished:
                continue  # already compressed; try next

            # xxx if ISPLITMJD > 0 : continue  # temp xxx REMOVE
            
            # Check if makeDataFile tasks have finished for this MJD range        
            splitmjd_done_list = [True] * nsplitmjd
            for row in row_merge_list:
                state     = row[COLNUM_MERGE_STATE]
                isplitmjd = row[COLNUM_MKDATA_MERGE_ISPLITMJD]
                if state != SUBMIT_STATE_DONE:
                    splitmjd_done_list[isplitmjd] = False

            if not splitmjd_done_list[ISPLITMJD]:
                continue    # avro file creation tasks not done; bye bye

            # compress it !
            wildcard = "mjd*"
            mjd_dir_list = sorted(glob.glob1(output_dir,wildcard))

            if nsplitmjd > 1 :
                min_edge_list = submit_info_yaml['MIN_MJD_EDGE']
                max_edge_list = submit_info_yaml['MAX_MJD_EDGE']
            else:
                min_edge_list = [ 10000 ]
                max_edge_list = [ 99000 ]

            min_edge = min_edge_list[ISPLITMJD]
            max_edge = max_edge_list[ISPLITMJD]
            
            time_0     = datetime.datetime.now()
            n_compress = self.compress_mjd_dirs(mjd_dir_list, 
                                                min_edge, max_edge)
            time_1     = datetime.datetime.now()
            time_dif   = (time_1 - time_0).total_seconds()
            rate       = n_compress / time_dif
            rate_str   = f"{rate:.1f}"

            irow = nrow - 1
            row_compress_list_new[irow][COLNUM_STATE]    = SUBMIT_STATE_DONE
            row_compress_list_new[irow][COLNUM_NMJD_DIR] = n_compress
            row_compress_list_new[irow][COLNUM_TIME]     = int(time_dif)
            row_compress_list_new[irow][COLNUM_RATE]     = float(rate_str)
            
        return row_compress_list_new
    
        # end compress_update_state
        
    
    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        # nothing to do here (yet)
        pass
        # end  merge_job_wrapup

    def compress_mjd_dirs_slow(self, mjd_dir_list, min_edge, max_edge):

        # xxx uses python tar util, which is much slower than os.system

        # For mjd_dirs in mjd_dir_list, compress those within
        # min_edge and max_edge-1.
        # "Compress"  mjd[mjd] diretory -> mjd[mjd].tar.gz
        
        output_dir = self.config_prep['output_dir']

        compress_done_file = self.get_compress_done_file(min_edge,max_edge)
        n_compress = 0

        imin = int(min_edge); imax = int(max_edge)-1
        logging.info(f"  Begin compression for " \
                     f"{ALERT_DAY_NAME}{imin} to " \
                     f"{ALERT_DAY_NAME}{imax} ")

        # loop over each mjd within range and tar it up
        sys.stdout.flush()
        for mjd_dir in mjd_dir_list:
            mjd         = int(mjd_dir[3:])
            do_compress = mjd>= min_edge and mjd < max_edge
            if do_compress :
                n_compress += 1
                tar_file     = f"{output_dir}/{ALERT_SUBDIR}/{mjd_dir}.tar.gz"
                mjd_dir_full = f"{output_dir}/{mjd_dir}"
                with tarfile.open(tar_file,"w:gz") as tar:
                    tar.add(mjd_dir_full,
                            arcname=os.path.basename(mjd_dir_full))
                shutil.rmtree(mjd_dir_full)                
        
        # touch done file to flag that this MJD range is compressed
        cmd_done = f"touch {compress_done_file}"
        os.system(cmd_done)

        return n_compress
    
        # end compress_mjd_dirs_slow

    def compress_mjd_dirs(self, mjd_dir_list, min_edge, max_edge):

        #
        # For mjd_dirs in mjd_dir_list, compress those within
        # min_edge and max_edge-1.
        # "Compress"  mjd[mjd] diretory -> mjd[mjd].tar.gz
        # using os.system ... it's faster than python tar.
        
        output_dir = self.config_prep['output_dir']

        compress_done_file = self.get_compress_done_file(min_edge,max_edge)
        n_compress = 0

        imin = int(min_edge); imax = int(max_edge)-1
        logging.info(f"  Begin compression for " \
                     f"{ALERT_DAY_NAME}mjd{imin} to " \
                     f"{ALERT_DAY_NAME}mjd{imax} ")

        # construct big tar command
        sys.stdout.flush()
        cmd_tar = f"cd {output_dir} ; "
        for mjd_dir in mjd_dir_list:
            mjd         = int(mjd_dir[3:])
            do_compress = mjd>= min_edge and mjd < max_edge
            if do_compress:
                n_compress += 1
                tar_file = f"{mjd_dir}.tar.gz"
                cmd_tar += f"tar -czf {tar_file} {mjd_dir} --remove-files ; "

        # tack on command to move it all into /ALERTS,
        # and ten run it all with one os command
        if n_compress > 0:
            cmd_tar += f"mv {ALERT_DAY_NAME}*.tar.gz {ALERT_SUBDIR}"
            os.system(cmd_tar)

        # touch done file to flag that this MJD range is compressed
        cmd_done = f"touch {compress_done_file}"
        os.system(cmd_done)

        return n_compress
    
        # end compress_mjd_dirs
                    
    def get_compress_done_file(self,min_edge,max_edge):
        # return name of done file for compress mjd range defined by
        # min_edge to max_edge (for LSST alerts)
        output_dir       = self.config_prep['output_dir']
        imin = int(min_edge)
        imax = int(max_edge)
        mjd_range_str = f"{ALERT_DAY_NAME}{imin}-{imax}"
        alert_dir = f"{output_dir}/{ALERT_SUBDIR}"
        done_file = f"{alert_dir}/compress_{mjd_range_str}.done"
        return done_file

        # end get_done_file_compress

    def get_misc_merge_info(self):

        # Called at end of all jobs, return misc info lines 
        # (yaml format) to write at the end of the MERGE.LOG file.
        # Each info yaml-compatible line must be of the form
        #  KEYNAME:  VALUE

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_format    = submit_info_yaml['OUTPUT_FORMAT']
        isfmt_snana      = (output_format == OUTPUT_FORMAT_SNANA)
        isfmt_lsst_alert = (output_format == OUTPUT_FORMAT_LSST_ALERTS)

        # sum NEVT column in MERGE.LOG
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        MERGE_INFO_CONTENTS, comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)
        
        row_list  = MERGE_INFO_CONTENTS[TABLE_MERGE]
        NEVT_SUM = 0
        NOBS_SUM = 0
        for row in row_list:            
            NEVT_SUM += int(row[COLNUM_MKDATA_MERGE_NEVT]) 
            if isfmt_lsst_alert:
                NOBS_SUM += int(row[COLNUM_MKDATA_MERGE_NOBS_ALERT]) 

        info_lines = [ f"NEVT_SUM:        {NEVT_SUM}" ]

        if isfmt_lsst_alert:
            info_lines += [ f"NOBS_ALERT_SUM:  {NOBS_SUM}" ]

            # count mjd*.tar files in ALERT_SUBDIR
            alert_dir    = f"{output_dir}/{ALERT_SUBDIR}"
            wildcard     = f"{ALERT_DAY_NAME}*.tar.gz"
            mjd_tar_list = glob.glob1(alert_dir, wildcard)
            ndir         = len(mjd_tar_list)
            info_lines  += [ f"NDIR_MJD_SUM:    {ndir}" ]

            # sum alert-compress times
            nsec_sum = 0.0
            row_list = MERGE_INFO_CONTENTS[TABLE_COMPRESS]
            for row in row_list: nsec_sum += row[COLNUM_COMPRESS_TIME] 
            t_compress = nsec_sum/60.0
            t_compress = float(f"{t_compress:.2f}") 
            info_lines += [ f"TIME_COMPRESS_SUM:  {t_compress}  # minutes" ]
            
        # - - - - -
        return info_lines

        # end get_misc_merge_info

    def merge_cleanup_final(self):
        # every makeDataFiles succeeded, so here we simply compress output.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        cwd              = submit_info_yaml['CWD']
        output_format    = submit_info_yaml['OUTPUT_FORMAT']
        isfmt_snana      = (output_format == OUTPUT_FORMAT_SNANA)
        isfmt_lsst_alert = (output_format == OUTPUT_FORMAT_LSST_ALERTS)
        msgerr = []
        
        if isfmt_snana :
            command_list = ['makeDataFiles.sh',
                            '--outdir_snana', output_dir, '--merge']
            ret = subprocess.run(command_list, capture_output=False, text=True )
            
        elif isfmt_lsst_alert :
            pass
            
        else:
            msgerr.append(f"Unknown format '{output_format}" )
            util.log_assert(False,msgerr) # just abort, no done stamp
            
        #print(ret,' debugging')

        wildcard_list = [ 'MAKEDATA', 'CPU',  ]
        for w in wildcard_list :
            wstar = f"{w}*"
            tmp_list = glob.glob1(script_dir,wstar)
            if len(tmp_list) == 0 : continue
            print(f"\t Compress {wstar}")
            sys.stdout.flush()
            util.compress_files(+1, script_dir, wstar, w, "" )

        # - - - -
        # tar up entire script dir
        util.compress_subdir(+1, script_dir)

        # end merge_cleanup_final

    def get_merge_COLNUM_CPU(self):
        return -9

# =========== END: =======



