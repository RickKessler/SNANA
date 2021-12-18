#!/usr/bin/env python
# Created Dec 2021 by Maria Vincenzi and Rick Kessler
# 
# Script to create SIMSED binaries as explained in
# snana manual section "9.5 SIMSED". 
# Syntax:
#    make_simsed_binaries.py <configFile>
#    make_simsed_binaries.py -H    # help on configFile
#
# A separate configFile is needed for each survey.
#
# ============================================

import os, argparse, logging, shutil, time
import re, yaml, sys, gzip, math, subprocess

# define a few globals
GENVERSION_PREFIX = 'SIMSED_BINARY'             # prefix for sim genversion
BASH_PREFIX       = f'MAKE_{GENVERSION_PREFIX}' # prefix for bash file
LOG_PREFIX        = BASH_PREFIX                 # prefix for log files
JOB_NAME          = 'snlc_sim.exe'              # sim executable name in snana

# 2=force new SED.BINARY, 4=force new fluxTable binary, 6=force both
USE_BINARY        = 4  

# define key names for input config file
KEYNAME_SIM_KEYS      = "SIM_KEYS"
KEYNANE_SIMSED_MODELS = "SIMSED_MODELS"
KEYNAME_NJOB_PARALLEL = "NJOB_PARALLEL"

# ============================
def setup_logging():

    #logging.basicConfig(level=logging.DEBUG,
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)

def read_yaml(path):
    path_expand = os.path.expandvars(path)
    logging.debug(f"Reading YAML from {path_expand}")
    with open(path_expand) as f:
        return yaml.safe_load(f.read())

def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP menu for config options"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of the yml config file to run"
    parser.add_argument("input_file", help=msg,
                        nargs="?", default=None)

    msg = "Create bash script but do not submit it "
    parser.add_argument("--nosubmit", "-n", help=msg,
                            action="store_true")
    # parse it
    args = parser.parse_args()

    #if len(sys.argv) == 1:
    #    parser.print_help()
    #    sys.exit()

    if args.HELP : 
        print_help_menu()

    return args
    # end get_args

def print_help_menu():

    help_menu = f"""

# example input config file for DES

# Provide a few sim-input keys 
{KEYNAME_SIM_KEYS}:
  GENFILTERS:         griz
  KCOR_FILE:          $SNDATA_ROOT/kcor/DES/DES-SN3YR/kcor_DECam.fits
  SIMLIB_FILE:        $DES_ROOT/simlibs/DES_DIFFIMG.SIMLIB
  SIMSED_PATH_BINARY: $PLASTICC_ROOT/model_libs/SIMSED_BINARIES

# provide models
{KEYNANE_SIMSED_MODELS}:    # model                         zmin  zmax
  - $PLASTICC_ROOT/model_libs/SIMSED.CART-MOSFIT  0.02  1.5
  - $PLASTICC_ROOT/model_libs/SIMSED.KN-K17       0.01  0.4
  - $PLASTICC_ROOT/model_libs/SIMSED.SNII-NMF     0.01  1.30
  - $PLASTICC_ROOT/model_libs/SIMSED.SNIa-91bg    0.01  0.5
  - etc ...

# optional parallel jobs (interactive; not batch)
{KEYNAME_NJOB_PARALLEL}: 3 # process in groups of 3

# after binaries are created, check processing status in 
#      MAKE_SIMSED_BINARY_SUMMARY.INFO
# Ideally zmax(SNR>5) should be about 0.1 to 0.2 below zmax.

"""
    
    sys.exit(f"{help_menu}")

    # end print_help_menu


def check_paths(info_dict):
    SIM_KEYS=info_dict['config'][KEYNAME_SIM_KEYS]
    model_path_list=info_dict['model_path_list']
    model_basename_list=info_dict['model_basename_list']
    path_dict = {'KCOR_FILE':   SIM_KEYS['KCOR_FILE'],
                 'SIMLIB_FILE': SIM_KEYS['SIMLIB_FILE'],
                 'SIMSED_PATH_BINARY': SIM_KEYS['SIMSED_PATH_BINARY']}
    for model, model_base in zip(model_path_list, model_basename_list):
        path_dict[model_base]=model
    found_error = False
    for key,path in path_dict.items():
        path_expand = os.path.expandvars(path)
        path_expand_gz = path_expand+'.gz'
        found_path = os.path.exists(path_expand) or os.path.exists(path_expand_gz) 
        if not found_path: 
            logging.warning(f'{key}: path not found (check {path})\n')
            found_error = True
    assert not found_error, "ABORT, fix missing file"
    
    return

def parse_models(info_dict):
    logging.info('Parse simsed models')
    SIMSED_MODELS=info_dict['config'][KEYNANE_SIMSED_MODELS]
    survey=info_dict['survey']
    model_path_list=[]
    model_zmin_list=[]
    model_zmax_list=[]
    model_basename_list=[]
    model_logfile_list=[]
    genversion_list=[]
    yaml_list=[]
    for line in SIMSED_MODELS:
        model=line.split()[0]
        model_base=os.path.basename(model)
        zmin=line.split()[1]
        zmax=line.split()[2]
        model_path_list.append(model)
        model_zmin_list.append(zmin)
       	model_zmax_list.append(zmax)
       	model_basename_list.append(model_base)
        version_prefix = f'{LOG_PREFIX}_{model_base}_{survey}'
        model_logfile_list.append(f'{version_prefix}.LOG')
        genversion_list.append(f'{version_prefix}')
        yaml_list.append(f'{version_prefix}.YAML')
    info_dict['model_path_list']=model_path_list
    info_dict['model_zmin_list']=model_zmin_list
    info_dict['model_zmax_list']=model_zmax_list
    info_dict['model_basename_list']=model_basename_list
    info_dict['model_logfile_list']=model_logfile_list
    info_dict['n_model']=len(model_path_list)
    info_dict['genversion_list']=genversion_list
    info_dict['yaml_list']=yaml_list
    return


def get_simlib_info(info_dict):
    logging.info(f'Extracting info from simlib')
    args=info_dict['args']
    SIM_KEYS=info_dict['config'][KEYNAME_SIM_KEYS]
    simlib_file_name = os.path.expandvars(SIM_KEYS['SIMLIB_FILE'])
     
    with open(simlib_file_name) as s:
        for line in s:
            if 'SURVEY:' in line:
                survey=line.split()[1]
            elif line[0:2]=='S:':
                start_mjd=float(line.split()[1])
                GENRANGE_PEAKMJD = f"{start_mjd+10.0}  {start_mjd+100.0}"
                break
    info_dict['survey']=survey
    info_dict['GENRANGE_PEAKMJD'] = GENRANGE_PEAKMJD
    logging.info(f'\t Found survey: {survey}')
    logging.info(f'\t Determined MJD range: {GENRANGE_PEAKMJD}')
    return


def create_simgen_file(info_dict):
    logging.info(f'Creating simulations input file')
    survey = info_dict['survey']
    args=info_dict['args']
    genversion=f'{GENVERSION_PREFIX}_{survey}'
    sim_input_file_name= f"SIM_{survey}.input"
    SIM_KEYS=info_dict['config'][KEYNAME_SIM_KEYS]
    GENRANGE_PEAKMJD = info_dict['GENRANGE_PEAKMJD']
    with open(sim_input_file_name, 'wt') as f:
        f.write(f"# Input keys passed from {args.input_file}\n")
        for k,v in SIM_KEYS.items():
            f.write(f'{k}: {v}\n')
        f.write(f'#\n# auto-generated keys with required values\n')
        f.write(f'GENVERSION: {genversion}\n')
        f.write(f'SIMSED_USE_BINARY: {USE_BINARY}\n')
        f.write(f'GENSOURCE:   RANDOM\n')
        f.write(f'GENRANGE_PEAKMJD: {GENRANGE_PEAKMJD}\n')
        f.write(f'#\n# auto-generated keys to prevent abort (values do not matter)\n')
        f.write(f'NGENTOT_LC: 300\n')
        f.write(f'RANSEED: 12945\n')
        f.write(f'GENRANGE_TREST:  -40  100\n') 
        f.write(f'DNDZ: FLAT\n')
        f.write(f'WRFLAG_YAML_FILE: 1\n')

    # load item in the dictionary
    info_dict['GENVERSION']=genversion
    info_dict['sim_input_file']=sim_input_file_name
    logging.info(f'Simgen input file created: {sim_input_file_name}')

    #end create_simgen_file

def create_bash_script(info_dict):
    survey = info_dict['survey']    
    sim_input_file      = info_dict['sim_input_file']
    model_path_list     = info_dict['model_path_list']
    model_zmin_list     = info_dict['model_zmin_list']
    model_zmax_list     = info_dict['model_zmax_list']
    model_basename_list = info_dict['model_basename_list']
    model_logfile_list  = info_dict['model_logfile_list']
    genversion_list     = info_dict['genversion_list']
    config              = info_dict['config']
    n_model             = info_dict['n_model']

    bash_name           = f'{BASH_PREFIX}_{survey}.sh'

    if KEYNAME_NJOB_PARALLEL in config:
        njobs = config[KEYNAME_NJOB_PARALLEL]
    else:
        njobs = 1

    with open(bash_name, 'wt') as b:
        for index, model, zmin, zmax, model_base, logfile, genversion\
            in zip(range(n_model), model_path_list,
                   model_zmin_list, model_zmax_list, 
                   model_basename_list, model_logfile_list, genversion_list):
            arg_string  = f'GENMODEL {model} '
            arg_string += f'GENRANGE_REDSHIFT {zmin} {zmax} '
            arg_string += f'GENVERSION {genversion} '

            output_line  = f'{JOB_NAME} {sim_input_file} {arg_string} '
       	    output_line	+= f' > {logfile} '
            if (index+1)%njobs!=0: 
                output_line +=' & '
            b.write(f'{output_line}\n')
    cmd = f'chmod +x {bash_name}'
    os.system(cmd)        
    info_dict['bash_name']=bash_name

def cleanup(info_dict):
    logging.info(f'Cleaning up old LOG files')
    model_logfile_list = info_dict['model_logfile_list']
    yaml_list = info_dict['yaml_list']
    for log,yml in zip(model_logfile_list, yaml_list):
        if os.path.exists(log):
            logging.info(f'\t Found old version of {log}, removing it.')
            os.remove(log)
        if os.path.exists(yml):
            logging.info(f'\t Found old version of {yml}, removing it.')
            os.remove(yml)
    return None

def submit_bash(info_dict):
    bash_name = info_dict['bash_name']
    cmd_list = [f'./{bash_name}']
    logging.info(f'Submitting bash job, {cmd_list}')

    command = f'sh ./{bash_name} & ' ## {log_base_name} &'
    ret = subprocess.run( [ command ], cwd=os.getcwd(),
               shell=True, capture_output=False, text=True )
    return

def monitor_sims(info_dict):
    logging.info('Monitoring sims:')
    sleep_time = 3
    n_model = info_dict['n_model']
    model_logfile_list = info_dict['model_logfile_list']
    yaml_list = info_dict['yaml_list']
    model_basename_list=info_dict['model_basename_list']

    n_done=0
    t_last_update_list = [0.0]*n_model # times
    done_list = [False]*n_model # times
    while n_done<n_model:
        time.sleep(sleep_time)
        for yml, log, index in zip(yaml_list, model_logfile_list, range(n_model)):
             if done_list[index]: 
                 continue
             if not os.path.exists(log): 
                 continue
             t_last_update = os.path.getmtime(log)
             if t_last_update_list[index]==t_last_update:
                 n_done += 1
                 done_list[index]=True
                 logging.info(f'\t Done LOG file {log}')
             else:
                 t_last_update_list[index]=t_last_update
    logging.info(f'\t {n_done} LOG files are done. Checking for output YAML')
    
    success_list=[False]*n_model
    all_success=True
    for yml,model_basename,log,index in zip(yaml_list, model_basename_list,\
                                            model_logfile_list, range(n_model)):
        if not os.path.exists(yml):
            logging.warning(f'\t Did not find YAML file for model {model_basename}.')
            logging.warning(f'\t Check {log} log file.')
            success_list[index]=False
            all_success=False
        else:
            success_list[index]=True
    if not all_success:
        logging.warning(f'One or more sims were unsuccessful.')
    else:
        logging.info(f'\t All YAMLs found.')
    info_dict['success_list']=success_list
    return 

def make_summary(info_dict):
    logging.info(r'Creating summary.')
    success_list = info_dict['success_list']
    n_model = info_dict['n_model']
    model_logfile_list = info_dict['model_logfile_list']
    yaml_list = info_dict['yaml_list']
    model_basename_list = info_dict['model_basename_list']
    summary_file = f'{BASH_PREFIX}_SUMMARY.INFO' 
    summary_header = f'# CREATION_DATE:  {time.strftime("%a, %d %b %Y %H:%M:%S ", time.gmtime())}\n'
    summary_header += f'# USERNAME:    {os.environ["USER"]}\n'
    summary_header += f'# HOST:    {os.environ["HOSTNAME"]}\n'
    summary_header += f'# CWD:     {os.getcwd()}\n# \n'
    summary_header += f'{"ROW:":<5}{"MODEL":<20} {"zmin":<7} {"zmax":<7} {"zmax(SNR>5)":<15} {"size(MB)":<10}\n'
    with open(summary_file, 'wt') as sum:
        sum.write(summary_header)
        for yml,model_basename,index,success in \
            zip(yaml_list, model_basename_list,\
                range(n_model), success_list):
            if success:
                yml_content = yaml.load(open(yml), Loader=yaml.FullLoader)
                [zmin,zmax] = yml_content['REDSHIFT_RANGE_FLUXTABLE'].split()
                zmax_snr = yml_content['REDSHIFT_MAX_SNR5']
                size = yml_content['SIZE_SIMSED_FLUXTABLE']
                sum.write(f'{index:<5}{model_basename:<20} {zmin:<7} {zmax:<7} {zmax_snr:<15} {size:<10}\n')
            else:
                sum.write(f'{index:<5}{model_basename:<20} {"FAIL":<7} {"FAIL":<7} {"FAIL":<15} {"FAIL":<10}\n')
    logging.info(f'Summary created in {summary_file}')
    return

# ===================================================
if __name__ == "__main__":
    try:
        setup_logging()
        logging.info("# ========== BEGIN make_simsed_binaries ===============")
        args   = get_args() 
        config = read_yaml(args.input_file)
        # print(f"xxxx config = {config}")        
        info_dict = {"args":args, "config":config}
        
        # extract survey name and approx MJD range from simlib
        get_simlib_info(info_dict)
        # parse list of models, paths and redshift ranges
        parse_models(info_dict) 
        
        # check paths of simsed models and sim input files exist
        check_paths(info_dict)

        # create minimal sim input file
        create_simgen_file(info_dict)

        # create bash script with snslc.exe commands
        create_bash_script(info_dict) 

        # cleaning up old LOG/YAML files. 
        # This is necessary for monitoring the status of the sims
        cleanup(info_dict)
        
        if args.nosubmit: sys.exit('Exiting without submitting bash script')
        # run bash script
        submit_bash(info_dict)        

        # monitor sim status
        monitor_sims(info_dict)
        
        # print out summary
        make_summary(info_dict)
        logging.info("# ========== END make_simsed_binaries ===============")

    except Exception as e:
        logging.exception(e)
        raise e



