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
        
        command_copy = f'cp {cigale_translator_file} {cigale_input_dir}/{cigale_translator_file}'

        command_exe = f'cd {cigale_input_dir}; {PROGRAM_CIGALE_TRANSLATOR} {cigale_translator_file} --mode SNANA_TO_CIGALE --output_cigale_file {CIGALE_CSV_FILE} --output_galid_map {GALID_MAP_FILE}'
        #command_exe += ' --zgrid 0.1 1.5 15' # temp hack

        os.mkdir(cigale_input_dir)
        os.system(command_copy)
        os.system(command_exe)

        #sys.exit(f'command_exe = {command_copy}')

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

        fitopt_num    = fitopt_dict['jobopt_num_list'][ijob] # e.g., "FITOPT000"
        prefix        = f"{SUBCLASS}"
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
        NGAL = 0 # fix this (have translator return ngal)

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

        KEY = FITOPT_STRING
        f.write(f"{KEY}: \n")
        for row in CONFIG[KEY] :
            f.write(f"- {row} \n")

        return
