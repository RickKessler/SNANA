import os, sys, shutil, yaml, glob, logging, gzip
import datetime, time, subprocess
import submit_util as util
#import pandas as pd
import re
from astropy.io import fits
import numpy as np

from   submit_params import *
from   submit_prog_base import Program

SUBCLASS_HOSTFIT_CIGALE = 'CIGALE'
SUBCLASS_HOSTFIT_CIGALE_LEGACY = 'CIGALE_LEGACY'
#PROGRAM_CIGALE_TRANSLATOR = '/home/jmedoff/SNANA/util/cigale_translator.py'
#PROGRAM_CIGALE_TRANSLATOR = 'cigale_translator.py'
CIGALE_INPUT_SUBDIR = 'CIGALE_INPUT'
CIGALE_CSV_FILE = 'cigale_input.in'
FITOPT_STRING = 'FITOPT'
MAX_GAL_PER_TASK_DEFAULT = 5000000

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


    def get_nrows_cigale(self):
        CONFIG     = self.config_yaml['CONFIG']
        cigale_translator_file = CONFIG['CIGALE_TRANSLATOR_FILE']
        
        with open(cigale_translator_file, "r") as f:
            cigale_translator_config = yaml.safe_load(f)
        input_table_file = cigale_translator_config["SNANA_TO_CIGALE"]["INPUT_TABLE_FILE"]
        redshift_grid = cigale_translator_config["SNANA_TO_CIGALE"].get("REDSHIFT_GRID")

        row_keys=("GAL", "SN") # Allows input table files with rows indicated by GAL: or SN: (matches with cigale_translator.py)

        keys = {k.rstrip(":") for k in row_keys}
        counts = {k: 0 for k in keys}

        opener = gzip.open if input_table_file.endswith(".gz") else open
        with opener(input_table_file, "rt") as f:
            for line in f:
                tokens = line.split(None, 1)
                if tokens:
                    tok = tokens[0].rstrip(":")
                    if tok in counts:
                        counts[tok] += 1

        present = {k: n for k, n in counts.items() if n > 0}
        if len(present) == 0:
            raise ValueError(f"{input_table_file}: no rows found for any of {sorted(keys)}")
        if len(present) > 1:
            raise ValueError(f"{input_table_file}: found multiple row types {present}; expected exactly one")
        
        nrows_input_table = next(iter(present.values()))
        nbins_zgrid = int(redshift_grid.split()[-1])
        return nrows_input_table*nbins_zgrid


    def submit_prepare_driver(self):
        args = self.config_yaml['args']
        SUBCLASS = SUBCLASS_HOSTFIT_CIGALE
        
        #if not args.devel_flag:
        #    self.submit_prepare_driver_legacy()
        #    return
        
        # Start devel here
        input_file = args.input_file
        output_dir = self.config_prep['output_dir']

        CONFIG     = self.config_yaml['CONFIG']
        MAX_GAL_PER_TASK = int(CONFIG.get('MAX_GAL_PER_TASK', MAX_GAL_PER_TASK_DEFAULT))

        nrows = self.get_nrows_cigale()
        nsplit = (nrows + MAX_GAL_PER_TASK - 1) // MAX_GAL_PER_TASK
        print(f'nrows = {nrows}, nsplit = {nsplit}')
        self.config_prep['cigale_input_nsplit'] = nsplit

        cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR
        self.prep_cigale_translator(cigale_input_dir)
        self.prep_cigale_fitopt()
        self.prep_cigale_symlinks()

        self.config_prep['cigale_input_dir'] = cigale_input_dir
        self.config_prep['n_job_tot'] = self.config_prep['n_core']
        self.config_prep['n_done_tot'] = self.config_prep['n_core']
        self.config_prep['n_job_split'] = 1
        self.config_prep['SUBCLASS'] = SUBCLASS

        logging.info(f'SUBCLASS = {SUBCLASS}')

        return

    def submit_prepare_driver_legacy(self):
        SUBCLASS = SUBCLASS_HOSTFIT_CIGALE_LEGACY
        args = self.config_yaml['args']
        input_file = args.input_file
        output_dir = self.config_prep['output_dir']

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
        args = self.config_yaml['args']
        #if not args.devel_flag:
        #    self.prep_cigale_translator_legacy(cigale_input_dir)
        #    return

        # Start devel here
        CONFIG     = self.config_yaml['CONFIG']
        cigale_translator_file = CONFIG['CIGALE_TRANSLATOR_FILE']
        program_cigale_translator = CONFIG['CIGALE_TRANSLATOR_SCRIPT']
        #cigale_csv_file = self.get_filepath(CIGALE_CSV_FILE, CIGALE_INPUT_SUBDIR)

        prescale_num = self.config_yaml['args'].prescale
        #print('xxx prescale = ', prescale_num)
        nsplit = self.config_prep['cigale_input_nsplit']

        command_copy = f'cp {cigale_translator_file} {cigale_input_dir}/{cigale_translator_file}'
        command_exe = f'cd {cigale_input_dir}; {program_cigale_translator} {cigale_translator_file} --mode SNANA_TO_CIGALE --output_cigale_file {CIGALE_CSV_FILE} --prescale {prescale_num} --nsplit {nsplit}'
        # Test nsplit
        #command_exe = f'cd {cigale_input_dir}; {program_cigale_translator} {cigale_translator_file} --mode SNANA_TO_CIGALE --output_cigale_file {CIGALE_CSV_FILE} --prescale {prescale_num} --nsplit 3'

        os.mkdir(cigale_input_dir)
        os.system(command_copy)
        # Execute via subprocess to get nrows output (this is now calculated above in get_nrows_cigale(), so I should maybe change this)
        result = subprocess.run(command_exe, shell=True, stdout=subprocess.PIPE, text=True, check=True)
        nrows = [int(x) for x in result.stdout.splitlines()] # list of ints (e.g., [50, 50, 50])
        self.config_prep['cigale_input_nrows'] = nrows

        return

    def prep_cigale_translator_legacy(self, cigale_input_dir):
        CONFIG     = self.config_yaml['CONFIG']
        cigale_translator_file = CONFIG['CIGALE_TRANSLATOR_FILE']
        program_cigale_translator = CONFIG['CIGALE_TRANSLATOR_SCRIPT']
        #cigale_csv_file = self.get_filepath(CIGALE_CSV_FILE, CIGALE_INPUT_SUBDIR)

        prescale_num = self.config_yaml['args'].prescale
        #print('xxx prescale = ', prescale_num)

        command_copy = f'cp {cigale_translator_file} {cigale_input_dir}/{cigale_translator_file}'

        command_exe = f'cd {cigale_input_dir}; {program_cigale_translator} {cigale_translator_file} --mode SNANA_TO_CIGALE --output_cigale_file {CIGALE_CSV_FILE} --prescale {prescale_num}'
        # Test nsplit
        #command_exe = f'cd {cigale_input_dir}; {program_cigale_translator} {cigale_translator_file} --mode SNANA_TO_CIGALE --output_cigale_file {CIGALE_CSV_FILE} --prescale {prescale_num} --nsplit 3'

        os.mkdir(cigale_input_dir)
        os.system(command_copy)
        # Execute via subprocess to get nrows output (this is now calculated above in get_nrows_cigale(), so I should change this)
        result = subprocess.run(command_exe, shell=True, stdout=subprocess.PIPE, text=True, check=True)
        nrows = int(result.stdout.strip())
        self.config_prep['cigale_input_nrows'] = nrows

        return

    def prep_cigale_fitopt(self):
        args = self.config_yaml['args']
        #if not args.devel_flag:
        #    self.prep_cigale_fitopt_legacy()
        #    return

        # Start devel here
        CONFIG     = self.config_yaml['CONFIG']
        KEYLIST       = [ FITOPT_STRING ]    # key under CONFIG
        fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)
        fitopt_dict = util.prep_jobopt_list(fitopt_rows,FITOPT_STRING,1,None)
        n_subset = self.config_prep['cigale_input_nsplit']

        self.config_prep['fitopt_dict'] = fitopt_dict
        self.config_prep['n_fitopt'] = fitopt_dict['n_jobopt']

        # CREATE WORKING DIR FOR EACH FITOPT (AND NOW ALSO SUBSET)
        output_dir = self.config_prep['output_dir']
        fitopt_dir_list = []
        subset_dir_list = []
        jobopt_num_list = fitopt_dict['jobopt_num_list']
        jobopt_arg_list = fitopt_dict['jobopt_arg_list']

        for fitopt_num, fitopt_arg in zip(jobopt_num_list, jobopt_arg_list):
            logging.info(f'prepare {fitopt_num}')
            fitopt_dir = f'{output_dir}/{fitopt_num}'
            fitopt_dir_list.append(fitopt_dir)
            os.mkdir(fitopt_dir)

            for s in range(n_subset):
                subset_dir = f'{fitopt_dir}/SUBSET{s:03d}'
                subset_dir_list.append(subset_dir)
                os.mkdir(subset_dir)   
                self.prep_pcigale_ini_file(subset_dir, fitopt_arg)

        self.config_prep['fitopt_dir_list'] = fitopt_dir_list
        self.config_prep['subset_dir_list'] = subset_dir_list

        return

    def prep_cigale_fitopt_legacy(self):
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
            self.prep_pcigale_ini_file_legacy(fitopt_dir, fitopt_arg)

        self.config_prep['fitopt_dir_list'] = fitopt_dir_list

        return

    def prep_cigale_symlinks(self):
        args = self.config_yaml['args']
        #if not args.devel_flag:
        #    self.prep_cigale_symlinks_legacy()
        #    return

        # Start devel here 
        fitopt_dict = self.config_prep['fitopt_dict']
        #fitopt_dir_list = self.config_prep['fitopt_dir_list']
        subset_dir_list = self.config_prep['subset_dir_list']
        #print(f'xxx subset_dir_list = {subset_dir_list}')
        script_dir    = self.config_prep['script_dir']
        n_subset = self.config_prep['cigale_input_nsplit']

        jobopt_num_list = fitopt_dict['jobopt_num_list']
        subset_num_list = [f'{jobopt}_SUBSET{s:03d}' for jobopt in jobopt_num_list for s in range(n_subset)]
        self.config_prep['subset_num_list'] = subset_num_list
        #print(f'xxx subset_num_list = {subset_num_list}')
        prefix        =self.get_prefix_name()
        sym_link_log_list = []
        sym_link_done_list = []

        #for fitopt_num, fitopt_dir in zip(jobopt_num_list, fitopt_dir_list):
        for subset_num, subset_dir in zip(subset_num_list, subset_dir_list):
            symlink_log_name = subset_num + '_' + f'{prefix}.LOG'
            log_file_orig = subset_dir + '/' +f'{prefix}.LOG'
            symlink_log_command = f'cd {script_dir}; ln -s {log_file_orig} {symlink_log_name}'

            symlink_done_name = subset_num + '_' + f'{prefix}.DONE'
            done_file_orig = subset_dir + '/' +f'{prefix}.DONE'
            symlink_done_command = f'cd {script_dir}; ln -s {done_file_orig} {symlink_done_name}'

            sym_link_log_list.append(symlink_log_command)
            sym_link_done_list.append(symlink_done_command)

        self.config_prep['sym_link_log_list'] = sym_link_log_list
        self.config_prep['sym_link_done_list'] = sym_link_done_list
        return

    def prep_cigale_symlinks_legacy(self):
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

    def prep_pcigale_ini_file(self, subset_dir, fitopt_arg):
        CONFIG     = self.config_yaml['CONFIG']
        nthread    = self.config_prep['nthreads']
        output_dir = self.config_prep['output_dir']
        cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR
        #n_subset = len(glob.glob(os.path.join(cigale_input_dir, "cigale_input_SUBSET*.in")))
        n_subset = self.config_prep['cigale_input_nsplit']

        pcigale_ini_file_orig = CONFIG['CIGALE_INPUT_FILE_LIST'].split()[0]
        pcigale_ini_spec_file_orig = CONFIG['CIGALE_INPUT_FILE_LIST'].split()[1]
        pcigale_ini_file_target = f'{subset_dir}/{pcigale_ini_file_orig}'
        pcigale_ini_spec_file_target = f'{subset_dir}/{pcigale_ini_spec_file_orig}'

        # MODIFY .ini FILE
        cigale_bands_str = self.get_cigale_bands_str()
        subset_num = os.path.basename(subset_dir) # e.g., SUBSET000
        root, ext = os.path.splitext(CIGALE_CSV_FILE)
        cigale_csv_subset_file = root + '_' + subset_num + ext
        cigale_datafile_str = cigale_input_dir + '/' + cigale_csv_subset_file
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

    def prep_pcigale_ini_file_legacy(self, fitopt_dir, fitopt_arg):
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
        n_core   = self.config_prep['n_core']
        n_fitopt = self.config_prep['n_fitopt']
        n_subset = self.config_prep['cigale_input_nsplit']
        n_job    = 0
        for ijob in range(n_fitopt*n_subset):
            n_job += 1  # track total number pof jobs; counter starts at 1, not 0
            if ijob % n_core == icpu:
            #if ijob == icpu:
                job_info_hostfit   = self.prep_JOB_INFO_hostfit(ijob)
                util.write_job_info(f, job_info_hostfit, icpu)
                # NEED TO FIX WHEN ncore != njob

                job_info_merge = self.prep_JOB_INFO_merge(icpu, n_job, False)
                util.write_jobmerge_info(f, job_info_merge, icpu)

        #sys.exit('xxx bye')
        return        

    def get_prefix_name(self):
        prefix        = f"{SUBCLASS_HOSTFIT_CIGALE}"
        return prefix

    def prep_JOB_INFO_hostfit(self, ijob):
        program       = self.config_prep['program']
        output_dir    = self.config_prep['output_dir']
        script_dir    = self.config_prep['script_dir']
        args          = self.config_yaml['args']
        input_file    = args.input_file
        SUBCLASS = self.config_prep['SUBCLASS']
        fitopt_dict = self.config_prep['fitopt_dict']
        #fitopt_dir = self.config_prep['fitopt_dir_list'][ijob]
        subset_dir = self.config_prep['subset_dir_list'][ijob]
        sym_log_link = self.config_prep['sym_link_log_list'][ijob]
        sym_done_link = self.config_prep['sym_link_done_list'][ijob]

        #fitopt_num    = fitopt_dict['jobopt_num_list'][ijob] # e.g., "FITOPT000"
        prefix        = self.get_prefix_name()
        done_file     = f"{prefix}.DONE"
        log_file      = f"{prefix}.LOG"
        yaml_file     = f"{prefix}.YAML"
        arg_list      = []
        start_file    = f'{prefix}.START'
        
        JOB_INFO      = {}
        #JOB_INFO['job_dir']     = fitopt_dir  # where to run job
        JOB_INFO['job_dir']     = subset_dir  # where to run job  
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
        n_subset = self.config_prep['cigale_input_nsplit']

        header_line_merge = \
                        f" STATE  FITOPT  NGAL  CPU  "

        INFO_MERGE = {
            'primary_key' : TABLE_MERGE,
            'header_line' : header_line_merge,
            'row_list'    : []
        }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        NGALs = self.config_prep['cigale_input_nrows'] # list of len n_subset

        for ijob in range(n_fitopt*n_subset):
            #fitopt_num    = fitopt_dict['jobopt_num_list'][ijob] # e.g., "FITOPT000"   
            subset_num = self.config_prep['subset_num_list'][ijob] # e.g., "FITOPT000_SUBSET000"
            _, subset_index = self.strip_fitopt_subset(subset_num, type_flag = int)
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(subset_num)
            ROW_MERGE.append(NGALs[subset_index])
            ROW_MERGE.append(0.0) # CPU
            INFO_MERGE['row_list'].append(ROW_MERGE)
        # - - - - -                                                           
        util.write_merge_file(f, INFO_MERGE, [] )

        return

    def append_info_file(self,f):
        CONFIG             = self.config_yaml['CONFIG']
        n_fitopt = self.config_prep['n_fitopt']
        #print('xxx n_fitopt ', n_fitopt)
        f.write("# BEGIN HOSTFIT INFO \n")
        f.write(f"JOBFILE_WILDCARD:    '{FITOPT_STRING}*' \n")

        KEY = FITOPT_STRING
        f.write(f"{KEY}: \n")
        # This loop should only run if there's another fitopt besides FITOPT000
        for ijob in range(n_fitopt-1):
            row = CONFIG[KEY][ijob]
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

    def get_cigale_status(self, fitopt_num):
        script_dir    = self.config_prep['script_dir']
        prefix = self.get_prefix_name()
        symlink_log_name = script_dir + '/' + fitopt_num + '_' + f'{prefix}.LOG'
        #sys.exit(f'xxx Log name = {symlink_log_name}')
        success = False
        tproc = None

        with open(symlink_log_name) as f:
            log = f.read()
        if "Run completed!" in log:
            success = True
        duration = re.search(r"Total duration:\s*(\d+):(\d+):(\d+)", log)
        if duration:
            h, mins, s = (int(x) for x in duration.groups())
            tproc = datetime.timedelta(hours=h, minutes=mins, seconds=s).total_seconds() / 60  # minutes
        return success, round(tproc, 4)

    def get_merge_COLNUM_CPU(self):
        return COLNUM_HOSTFIT_MERGE_CPU

    def get_filepath(self, filename, subdir = ''):
        output_dir = self.config_prep['output_dir']
        filepath = output_dir + '/' + subdir + '/' + filename
        return filepath

    def strip_fitopt_subset(self, fitopt_subset, type_flag = str):
        if type_flag is str:
            fitopt, subset = fitopt_subset.split('_')
        elif type_flag is int:
            fitopt_str, subset_str = fitopt_subset.split('_')
            fitopt = int(fitopt_str[-3:])
            subset = int(subset_str[-3:])
        else:
            raise ValueError(f"type_flag should be str or int. Given: {type_flag}")
        return fitopt, subset

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        row  = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        fitopt_subset_num = row[COLNUM_HOSTFIT_MERGE_FITOPT]
        fitopt_num, subset_num = self.strip_fitopt_subset(fitopt_subset_num)
        fitopt_subset_dir = fitopt_num + '/' + subset_num

        CONFIG     = self.config_yaml['CONFIG']
        #KEYLIST       = [ FITOPT_STRING ]    # key under CONFIG
        #fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)
        #fitopt_dict = util.prep_jobopt_list(fitopt_rows,FITOPT_STRING,1,None)
        #fitopt_num = fitopt_dict['jobopt_num_list'][irow]
        #cigale_results_subdir =fitopt_num + '/out'
        cigale_results_subdir =fitopt_subset_dir + '/out'
        program_cigale_translator = CONFIG['CIGALE_TRANSLATOR_SCRIPT']

        output_dir = self.config_prep['output_dir']
        cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR
        #cigale_translator_file = cigale_input_dir + '/' + CONFIG['CIGALE_TRANSLATOR_FILE']

        cigale_translator_file = self.get_filepath(CONFIG['CIGALE_TRANSLATOR_FILE'], CIGALE_INPUT_SUBDIR)
        cigale_result_file = self.get_filepath('results.fits', cigale_results_subdir)
        output_snana_file = self.get_filepath('LOGMASS_GRID.DAT.gz', fitopt_subset_dir)

        command_exe = f'{program_cigale_translator} {cigale_translator_file} --mode CIGALE_TO_SNANA --input_cigale_results {cigale_result_file} --output_snana_file {output_snana_file}'
        os.system(command_exe)

        num_invalid, tot_len = self.count_invalid_logmasses(cigale_result_file)
        print(f'{num_invalid} out of {tot_len} rows ({num_invalid/tot_len*100:.4}%) in cigale results are NaN.')

        with open(cigale_translator_file, "r") as f:
            cigale_translator_config = yaml.safe_load(f)
        input_table_file = cigale_translator_config["SNANA_TO_CIGALE"]["INPUT_TABLE_FILE"]
        self.update_output_documentation(output_snana_file, input_table_file)

        # Check if all subsets are done (count DONE files) - if so, combine them
        #n_done = len(glob.glob(f"{fitopt_num}/SUBSET*/*.DONE"))
        n_done = len(glob.glob(os.path.join(output_dir, f"{fitopt_num}/SUBSET*/LOGMASS_GRID.DAT.gz")))
        n_subset = len(glob.glob(os.path.join(cigale_input_dir, "cigale_input_SUBSET*.in")))
        #print(f'xxx n_done, n_subset = {n_done}, {n_subset}')
        if n_done == n_subset:
            subset_output_files = []
            for s in range(n_subset):
                output_snana_file = self.get_filepath('LOGMASS_GRID.DAT.gz', fitopt_num + f'/SUBSET{s:03d}')
                subset_output_files.append(output_snana_file)
            combined_output_file = self.get_filepath('LOGMASS_GRID.DAT.gz', fitopt_num)
            self.combine_logmass_grids(subset_output_files, combined_output_file)

        return

    def count_invalid_logmasses(self, cigale_result_file):
        with fits.open(cigale_result_file) as hdul:
            dat = hdul[1].data
        num_invalid = np.sum(np.isnan(dat['bayes.stellar.mass_total']))
        tot_len = len(dat)
        return num_invalid, tot_len

    def update_output_documentation(self, output_snana_file, input_file):
        # Add input host information to output file documentation
        # (reads input_file directly from pcigale_translator.input)
        opener = gzip.open if output_snana_file.endswith(".gz") else open
        with opener(output_snana_file, 'rt') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            stripped = line.lstrip()
            if stripped.startswith("INPUT_HOST_INFORMATION"):    
                indent = line[:len(line) - len(stripped)]
                line = f'{indent}INPUT_HOST_INFORMATION: {input_file}\n'
            new_lines.append(line)

        with opener(output_snana_file, 'wt') as f:
            f.writelines(new_lines)

        return

    def get_misc_merge_info(self):
        lines = []

        CONFIG     = self.config_yaml['CONFIG']
        n_core   = CONFIG['BATCH_NTHREADS']

        KEYLIST       = [ FITOPT_STRING ]    # key under CONFIG                                                                                                                                                
        fitopt_rows   = util.get_YAML_key_values(CONFIG,KEYLIST)
        fitopt_dict = util.prep_jobopt_list(fitopt_rows,FITOPT_STRING,1,None)

        output_dir = self.config_prep['output_dir']
        cigale_input_dir = output_dir + '/' + CIGALE_INPUT_SUBDIR
        n_subset = len(glob.glob(os.path.join(cigale_input_dir, "cigale_input_SUBSET*.in")))

        script_dir    = self.config_prep['script_dir']
        prefix = self.get_prefix_name()

        CPU_PER_GAL_LIST = []
        for fitopt_num in fitopt_dict['jobopt_num_list']:
            for s in range(n_subset):
                symlink_log_name = script_dir + '/' + f'{fitopt_num}_SUBSET{s:03d}_{prefix}.LOG'

                with open(symlink_log_name) as f:
                    log = f.read()

                duration = re.search(r"Total duration:\s*(\d+):(\d+):(\d+)", log)
                num_objects = re.search(r"Number of objects\s*[│|]\s*(\d+)", log)
        
                if duration and num_objects:
                    h, mins, s = (int(x) for x in duration.groups())
                    tproc = datetime.timedelta(hours=h, minutes=mins, seconds=s).total_seconds() / 60  # minutes
                    ngal = int(num_objects.group(1))
                    CPU_PER_GAL = tproc * n_core / ngal # minutes
                    CPU_PER_GAL_LIST.append(CPU_PER_GAL)
        MEAN_CPU_PER_GAL = round(np.mean(CPU_PER_GAL_LIST), 4)
        LINE = f'CPU_PER_GALAXY: {MEAN_CPU_PER_GAL}  # minutes'
        lines.append(LINE)

        # Hard code FITOPT000 since number of NaN logmasses should be the same for all fitopts
        cigale_results_subdir = 'FITOPT000/SUBSET000/out'
        cigale_result_file = self.get_filepath('results.fits', cigale_results_subdir)
        num_invalid, tot_len = self.count_invalid_logmasses(cigale_result_file)
        LINE_nan_num = f'NaN_LOGMASSES: {num_invalid}'
        LINE_nan_percent = f'%_NaN_LOGMASSES: {num_invalid/tot_len*100:.4}%'
        lines.append(LINE_nan_num)
        lines.append(LINE_nan_percent)

        return lines

    def merge_cleanup_final(self):
        return

    def combine_logmass_grids(self, input_files, output_file):
        """Concatenate LOGMASS_GRID files from different subsets"""
        input_files = [str(p) for p in input_files]
        if len(input_files) < 2:
            raise ValueError("need at least two input files to combine")

        header_lines = None
        ref_varnames = None
        grids = {}
        n_dupes = 0
        cigale_files = []

        for ifile, path in enumerate(input_files):
            header, varnames, n_read = [], None, 0
            seen_varnames = in_body = False

            opener = gzip.open if path.endswith(".gz") else open
            with opener(path, 'rt') as f:
                for lineno, raw in enumerate(f, start=1):
                    line = raw.rstrip("\n")
                    stripped = line.strip()
                    if not stripped.startswith("GAL:"):
                        if not in_body:
                            header.append(line)
                            if stripped.startswith("VARNAMES:"):
                                seen_varnames = True
                                varnames = line.split()[1:]
                        continue

                    in_body = True
                    n_read += 1
                    parts = stripped.split()
                    if len(parts) < 3:
                        raise ValueError(f"{path}:{lineno}: GAL: row has too few columns: {stripped!r}")
                    galid = parts[1]
                    try:
                        zgrid = float(parts[2])
                    except ValueError:
                        raise ValueError(f"{path}:{lineno}: cannot parse ZGRID {parts[2]!r} as float")

                    zmap = grids.setdefault(galid, {})
                    zkey = round(zgrid, 6)   # avoid float-repr mismatches between files
                    if zkey not in zmap:
                        zmap[zkey] = [line]
                        continue
                    n_dupes += 1

            if not seen_varnames:
                raise ValueError(f"{path}: no VARNAMES line found")
            if n_read == 0:
                raise ValueError(f"{path}: no GAL: rows found")

            for value in self.cigale_block(header)[0]:
                if value not in cigale_files:
                    cigale_files.append(value)

            if ifile == 0:
                header_lines, ref_varnames = header, varnames
            elif varnames != ref_varnames:
                print(f"WARNING: VARNAMES in {path} differ from {input_files[0]}\n"
                      f"         {input_files[0]}: {ref_varnames}\n"
                      f"         {path}: {varnames}",
                      file=sys.stderr)

        galids = list(grids)

        # List all the CIGALE_RESULTS_FILEs in the combined LOGMASS_GRID file
        if cigale_files:
            _, istart, iend, indent = self.cigale_block(header_lines)
            #print(f'xxx len(cigale_files) = {len(cigale_files)}, istart = {istart}')
            if len(cigale_files) == 1:
                header_lines[istart:iend] = [f"{indent}CIGALE_RESULTS_FILE: {cigale_files[0]}"]
            else:
                header_lines[istart:iend] = ([f"{indent}CIGALE_RESULTS_FILE:"]
                                         + [f"{indent}- {value}" for value in cigale_files])

        n_rows = 0
        opener = gzip.open if output_file.endswith(".gz") else open
        with opener(output_file, "wt") as out:
            for line in header_lines:
                out.write(line + "\n")
            for igal, galid in enumerate(galids):
                if igal:
                    out.write("\n")
                zmap = grids[galid]
                for zkey in sorted(zmap):
                    for line in zmap[zkey]:
                        out.write(line + "\n")
                        n_rows += 1

        print('Combined LOGMASS_GRID subsets:')
        print(f'n_files: {len(input_files)}')
        print(f'n_galid: {len(galids)}')
        print(f'n_rows: {n_rows}')
        print(f'n_duplicates: {n_dupes}')

        return

    def cigale_block(self, header):
        """Locate CIGALE_RESULTS_FILE in documentation block"""
        for i, line in enumerate(header):
            stripped = line.strip()
            if not stripped.startswith("CIGALE_RESULTS_FILE:"):
                continue
            indent = line[:len(line) - len(line.lstrip())]
            inline = stripped[len("CIGALE_RESULTS_FILE:"):].strip()
            if inline:
                return [inline], i, i + 1, indent
            # YAML-list form: consume the "- item" lines that follow
            values, j = [], i + 1
            while j < len(header) and header[j].strip().startswith("- "):
                values.append(header[j].strip()[2:].strip())
                j += 1
            return values, i, j, indent
        return [], None, None, ""

