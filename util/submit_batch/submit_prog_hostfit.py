import os, sys, shutil, yaml, glob, logging
import datetime, time, subprocess
import submit_util as util
import pandas as pd
import re

from   submit_params import *
from   submit_prog_base import Program

SUBCLASS_HOSTFIT_CIGALE = 'CIGALE'
PROGRAM_CIGALE_TRANSLATOR = '/home/jmedoff/SNANA/util/cigale_translator.py'
#PROGRAM_CIGALE_TRANSLATOR = 'cigale_translator.py'
CIGALE_INPUT_SUBDIR = 'CIGALE_INPUT'
CIGALE_CSV_FILE = 'cigale_input.in'
GALID_MAP_FILE = 'galid_map.csv'
FITOPT_STRING = 'FITOPT'

# define columns for MERGE.LOG;  column 0 is always for STATE                  
COLNUM_HOSTFIT_MERGE_FITOPT      = 1
COLNUM_HOSTFIT_MERGE_NGAL        = 2
COLNUM_HOSTFIT_MERGE_CPU         = 3

# Created May 2026
# run SED fit code on galaxies to get host properties.
#
class HostPropertyFit(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = None
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

        
        return output_dir_name, SUBDIR_SCRIPTS_HOSTFIT

        
    def submit_prepare_driver(self):
        input_file = self.config_yaml['args'].input_file
        output_dir = self.config_prep['output_dir']

        #SNANA_TO_CIGALE = util.extract_yaml(input_file, "SNANA_TO_CIGALE:", None)

        #CIGALE_TO_SNANA = util.extract_yaml(input_file, "CIGALE_TO_SNANA:", None)

        SUBCLASS = SUBCLASS_HOSTFIT_CIGALE
        
        cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR
        self.prep_cigale_translator(cigale_input_dir)
        self.prep_cigale_fitopt()
        self.prep_cigale_symlinks()

        self.config_prep['cigale_input_dir'] = cigale_input_dir
        self.config_prep['n_job_tot'] = self.config_prep['n_core']
        self.config_prep['n_done_tot'] = self.config_prep['n_core']
        self.config_prep['n_job_split'] = 1
        self.config_prep['SUBCLASS'] = SUBCLASS
        #self.config_prep['SNANA_TO_CIGALE'] = SNANA_TO_CIGALE
        #self.config_prep['CIGALE_TO_SNANA'] = CIGALE_TO_SNANA
        logging.info(f'SUBCLASS = {SUBCLASS}')

        #CONFIG_HOSTFIT     = self.config_yaml['SNANA_TO_CIGALE']
        #FILTERS       = CONFIG.setdefault(KEY_FILTERS, None)
        #sys.exit(self.config_yaml)
        return

    def prep_cigale_translator(self, cigale_input_dir):
        CONFIG     = self.config_yaml['CONFIG']
        cigale_translator_file = CONFIG['CIGALE_TRANSLATOR_FILE']
        #galid_map_file = self.get_filepath(GALID_MAP_FILE, CIGALE_INPUT_SUBDIR)
        #cigale_csv_file = self.get_filepath(CIGALE_CSV_FILE, CIGALE_INPUT_SUBDIR)

        command_copy = f'cp {cigale_translator_file} {cigale_input_dir}/{cigale_translator_file}'

        command_exe = f'cd {cigale_input_dir}; {PROGRAM_CIGALE_TRANSLATOR} {cigale_translator_file} --mode SNANA_TO_CIGALE --output_cigale_file {CIGALE_CSV_FILE} --output_galid_map {GALID_MAP_FILE}'

        os.mkdir(cigale_input_dir)
        os.system(command_copy)
        # Execute via subprocess to get nrows output
        result = subprocess.run(command_exe, shell=True, stdout=subprocess.PIPE, text=True, check=True)
        nrows = int(result.stdout.strip())
        self.config_prep['cigale_input_nrows'] = nrows

        return


    def prep_cigale_fitopt(self):
        CONFIG     = self.config_yaml['CONFIG']
        KEYLIST       = [ FITOPT_STRING ]    # key under CONFIG
        fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)
        fitopt_dict = util.prep_jobopt_list(fitopt_rows,FITOPT_STRING,1,None)
        
        self.config_prep['fitopt_dict'] = fitopt_dict
        self.config_prep['n_fitopt'] = fitopt_dict['n_jobopt']
        
        # CREATE WORKING DIR FOR EACH FITOPT
        output_dir = self.config_prep['output_dir'] 
        fitopt_dir_list = []
        jobopt_num_list = fitopt_dict['jobopt_num_list']
        jobopt_arg_list = fitopt_dict['jobopt_arg_list']
        
        for fitopt_num, fitopt_arg in zip(jobopt_num_list, jobopt_arg_list):
            logging.info(f'prepare {fitopt_num}')
            fitopt_dir = f'{output_dir}/{fitopt_num}'
            fitopt_dir_list.append(fitopt_dir)
            os.mkdir(fitopt_dir)
            self.prep_pcigale_ini_file(fitopt_dir, fitopt_arg)

        self.config_prep['fitopt_dir_list'] = fitopt_dir_list

        return

    def prep_cigale_symlinks(self):
        fitopt_dict = self.config_prep['fitopt_dict']
        fitopt_dir_list = self.config_prep['fitopt_dir_list']
        script_dir    = self.config_prep['script_dir']

        jobopt_num_list = fitopt_dict['jobopt_num_list']
        prefix        =self.get_prefix_name()
        sym_link_log_list = []
        sym_link_done_list = []

        for fitopt_num, fitopt_dir in zip(jobopt_num_list, fitopt_dir_list):
            symlink_log_name = fitopt_num + '_' + f'{prefix}.LOG'
            log_file_orig = fitopt_dir + '/' +f'{prefix}.LOG' 
            symlink_log_command = f'cd {script_dir}; ln -s {log_file_orig} {symlink_log_name}'

            symlink_done_name = fitopt_num + '_' + f'{prefix}.DONE'
            done_file_orig = fitopt_dir + '/' +f'{prefix}.DONE'
            symlink_done_command = f'cd {script_dir}; ln -s {done_file_orig} {symlink_done_name}'
            
            sym_link_log_list.append(symlink_log_command)
            sym_link_done_list.append(symlink_done_command)

        self.config_prep['sym_link_log_list'] = sym_link_log_list
        self.config_prep['sym_link_done_list'] = sym_link_done_list
        return

    def fitopt_str_to_dict(self, s):
        pattern = re.compile(r'(\w+)\s*=\s*')
        matches = list(pattern.finditer(s))
        result = {}
        for i, m in enumerate(matches):
            key = m.group(1)
            start = m.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(s)
            result[key] = s[start:end].strip()
        return result

    def get_cigale_bands_str(self):
        CONFIG     = self.config_yaml['CONFIG']
        cigale_translator_file = CONFIG['CIGALE_TRANSLATOR_FILE']

        with open(cigale_translator_file, "r") as f:
            cigale_translator_config = yaml.safe_load(f)
    
        mag_map_keys = cigale_translator_config["SNANA_TO_CIGALE"]["CIGALE_MAG_MAP"].keys()
        cigale_bands_str = ''

        for mag_key in list(mag_map_keys):
            #cigale_bands_str += f'{mag_key}, {mag_key}_err, '
            cigale_bands_str += f'{mag_key}, '
        cigale_bands_str = cigale_bands_str[:-2]

        return cigale_bands_str

    def replace_keys_pcigale(self, pcigale_ini_file_target, fitopt_arg_dict):
        # SUBSTITUTE KEYS IN COPIED CIGALE INPUT FILE                         
        with open(pcigale_ini_file_target, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            stripped = line.lstrip()
            if '=' in stripped and not stripped.startswith('#'):
                key = stripped.split('=', 1)[0].strip()
                if key in fitopt_arg_dict:
                    indent = line[:len(line) - len(stripped)]
                    if key == 'bands' and indent == '': # Special case for first 'bands' parameter
                        bands_str = [x.strip() for x in fitopt_arg_dict[key].split(',')]
                        bands_str_err = ', '.join(item for band in bands_str for item in (band, f'{band}_err'))
                        line = f'{indent}{key} = {bands_str_err}\n'
                    else:
                        line = f'{indent}{key} = {fitopt_arg_dict[key]}\n'
            new_lines.append(line)

        with open(pcigale_ini_file_target, 'w') as f:
            f.writelines(new_lines)

        return

    def prep_pcigale_ini_file(self, fitopt_dir, fitopt_arg):
        CONFIG     = self.config_yaml['CONFIG']
        nthread    = self.config_prep['nthreads']
        output_dir = self.config_prep['output_dir']
        cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR

        pcigale_ini_file_orig = CONFIG['CIGALE_INPUT_FILE_LIST'].split()[0]
        pcigale_ini_spec_file_orig = CONFIG['CIGALE_INPUT_FILE_LIST'].split()[1]
        pcigale_ini_file_target = f'{fitopt_dir}/{pcigale_ini_file_orig}'
        pcigale_ini_spec_file_target = f'{fitopt_dir}/{pcigale_ini_spec_file_orig}'
        
        # MODIFY .ini FILE
        cigale_bands_str = self.get_cigale_bands_str()
        cigale_datafile_str = cigale_input_dir + '/' + CIGALE_CSV_FILE
        fitopt_arg_full = f'data_file = {cigale_datafile_str} cores = {nthread} bands = {cigale_bands_str} {fitopt_arg}'
        fitopt_arg_dict = self.fitopt_str_to_dict(fitopt_arg_full)
        

        # FUTURE MODIFICATION OF .ini.spec FILE?

        command_copy = f'cp {pcigale_ini_file_orig} {pcigale_ini_file_target}'
        os.system(command_copy)
        command_copy = f'cp {pcigale_ini_spec_file_orig} {pcigale_ini_spec_file_target}'
        os.system(command_copy)

        self.replace_keys_pcigale(pcigale_ini_file_target, fitopt_arg_dict)
        self.replace_keys_pcigale(pcigale_ini_spec_file_target, {})

        return


    def write_command_file(self, icpu, f):
        n_core = self.config_prep['n_core']
        n_fitopt = self.config_prep['n_fitopt']

        for ijob in range(n_fitopt):
            if ijob == icpu:
                job_info_hostfit   = self.prep_JOB_INFO_hostfit(ijob)
                util.write_job_info(f, job_info_hostfit, icpu)
                # NEED TO FIX WHEN ncore != njob

        #sys.exit('xxx bye')
        return        

    def get_prefix_name(self):
        prefix        = f"{SUBCLASS_HOSTFIT_CIGALE}"
        return prefix

    def prep_JOB_INFO_hostfit(self, ijob):
        #SNANA_TO_CIGALE = self.config_prep['SNANA_TO_CIGALE']
        program       = self.config_prep['program']
        output_dir    = self.config_prep['output_dir']
        script_dir    = self.config_prep['script_dir']
        args          = self.config_yaml['args']
        input_file    = args.input_file
        SUBCLASS = self.config_prep['SUBCLASS']
        fitopt_dict = self.config_prep['fitopt_dict']
        fitopt_dir = self.config_prep['fitopt_dir_list'][ijob]
        sym_log_link = self.config_prep['sym_link_log_list'][ijob]
        sym_done_link = self.config_prep['sym_link_done_list'][ijob]

        fitopt_num    = fitopt_dict['jobopt_num_list'][ijob] # e.g., "FITOPT000"
        prefix        = self.get_prefix_name()
        done_file     = f"{prefix}.DONE"
        log_file      = f"{prefix}.LOG"
        yaml_file     = f"{prefix}.YAML"
        arg_list      = []
        start_file    = f'{prefix}.START'
        
        JOB_INFO      = {}
        JOB_INFO['job_dir']     = fitopt_dir  # where to run job
        JOB_INFO['program']     = program
        JOB_INFO['input_file']  = 'run'
        JOB_INFO['log_file']    = log_file
        JOB_INFO['done_file']   = done_file
        JOB_INFO['all_done_file'] = f"{output_dir}/{DEFAULT_DONE_FILE}"
        JOB_INFO['start_file'] = start_file
        JOB_INFO['arg_list'] = arg_list
        JOB_INFO['sym_link_list'] = [sym_log_link, sym_done_link]

        return JOB_INFO

    def create_merge_table(self,f):
        n_fitopt = self.config_prep['n_fitopt']
        fitopt_dict = self.config_prep['fitopt_dict']

        header_line_merge = \
                        f" STATE  FITOPT  NGAL  CPU  "

        INFO_MERGE = {
            'primary_key' : TABLE_MERGE,
            'header_line' : header_line_merge,
            'row_list'    : []
        }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        NGAL = self.config_prep['cigale_input_nrows']

        for ijob in range(n_fitopt):
            fitopt_num    = fitopt_dict['jobopt_num_list'][ijob] # e.g., "FITOPT000"   
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(fitopt_num)
            ROW_MERGE.append(NGAL)
            ROW_MERGE.append(0.0) # CPU
            INFO_MERGE['row_list'].append(ROW_MERGE)
        # - - - - -                                                           
        util.write_merge_file(f, INFO_MERGE, [] )

        return

    def append_info_file(self,f):
        CONFIG             = self.config_yaml['CONFIG']
        f.write("# BEGIN HOSTFIT INFO \n")
        f.write(f"JOBFILE_WILDCARD:    '{FITOPT_STRING}*' \n")

        KEY = FITOPT_STRING
        f.write(f"{KEY}: \n")
        for row in CONFIG[KEY] :
            f.write(f"- {row} \n")

        return

    def merge_config_prep(self,output_dir):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        self.config_prep['output_dir']     = output_dir
        self.config_prep['script_dir']     = script_dir
        
        return

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []
        row_list_merge     = MERGE_INFO_CONTENTS[TABLE_MERGE]

        nrow_check = 0
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input
            nrow_check += 1
            irow        = nrow_check - 1 # row index   
            fitopt    = row[COLNUM_HOSTFIT_MERGE_FITOPT] # e.g., FITOPT001
            search_wildcard = (f"{fitopt}*")

            # strip off row info 
            STATE       = row[COLNUM_MERGE_STATE]

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

                    # since there is no YAML file to examine, we have a 
                    # kluge check on success
                    success,tproc = self.get_cigale_status(fitopt)
                    if not success :
                        self.check_for_failure(log_list[0], -1, +1)
                        NEW_STATE = SUBMIT_STATE_FAIL

                    row[COLNUM_MERGE_STATE]     = NEW_STATE
                    row[COLNUM_HOSTFIT_MERGE_CPU]       = tproc

                    #COLNUM_HOSTFIT_MERGE_FITOPT      = 1
                    #COLNUM_HOSTFIT_MERGE_NGAL        = 2
                    #COLNUM_HOSTFIT_MERGE_CPU         = 3

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1

        # - - - - - -  -     
        # The first return arg (row_split) is null since there is
        # no need for a SPLIT table
        row_list_dict = {
            'row_split_list' : [],
            'row_merge_list' : row_list_merge_new,
            'row_extra_list' : []
        }

        return row_list_dict, n_state_change
        # end merge_update_state 

    def get_cigale_status(self, fitopt):
        success = True
        tproc = 1.5 # minutes
        return success, tproc

    def get_merge_COLNUM_CPU(self):
        return COLNUM_HOSTFIT_MERGE_CPU

    def get_filepath(self, filename, subdir = ''):
        output_dir = self.config_prep['output_dir']
        filepath = output_dir + '/' + subdir + '/' + filename
        return filepath

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        CONFIG     = self.config_yaml['CONFIG']
        KEYLIST       = [ FITOPT_STRING ]    # key under CONFIG
        fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)
        fitopt_dict = util.prep_jobopt_list(fitopt_rows,FITOPT_STRING,1,None)
        fitopt_num = fitopt_dict['jobopt_num_list'][irow]
        cigale_results_subdir =fitopt_num + '/out'

        #output_dir = self.config_prep['output_dir']
        #cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR
        #cigale_translator_file = cigale_input_dir + '/' + CONFIG['CIGALE_TRANSLATOR_FILE']
        
        cigale_translator_file = self.get_filepath(CONFIG['CIGALE_TRANSLATOR_FILE'], CIGALE_INPUT_SUBDIR)
        galid_map_file = self.get_filepath(GALID_MAP_FILE, CIGALE_INPUT_SUBDIR)
        cigale_result_file = self.get_filepath('results.fits', cigale_results_subdir)
        output_snana_file = self.get_filepath('LOGMASS_GRID.DAT', fitopt_num)

        row  = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        fitopt_num = row[COLNUM_HOSTFIT_MERGE_FITOPT]
        
        command_exe = f'{PROGRAM_CIGALE_TRANSLATOR} {cigale_translator_file} --mode CIGALE_TO_SNANA --input_cigale_results {cigale_result_file} --output_galid_map {galid_map_file} --output_snana_file {output_snana_file}'

        os.system(command_exe)

        return
