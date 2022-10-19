# Created Oct 4 2021 by R.Kessler
# Point to directory created by "create_covariance.py"
# and run wfit on hubble_diagram with all covsys_* files
# and with all user-specified WFITOPT options.
#
# Jan 17 2022 M.Vincenzi - add option to make wgt averages
#                          See KEYNAME_WFITAVG
#
# Feb 22 2022 RK - write nwarn column (read from yaml file from wfit)
# Feb 23 2022 RK - allow WFITOPT or WFITOPTS key
# Feb 24 2022 RK,A.Mitra: new COVOPTS key to select subset
# Feb 26 2022 RK - minor refac for WFITAVG feature
# Aug 16 2022 RK - if -outfile_chi2grid is passed, replace its arg
#                  with a standard file name.
# Oct 02 2022 RK - add and test merge_reset() util [standard method]
#
# Oct 13 2022 RK
#   + refactor to read covsys file names and HD filename from INFO.YML;
#     allows more flexibility in file names.
#     See new method read_hd_info_file
#
# Oct 1 82022 RK - change class name from wFit to cosmofit (more general name)
# ====================================================================

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time
import submit_util as util
import numpy as np

from submit_params    import *
from submit_prog_base import Program

# ------------------------------------------------
# xxx mark delete PREFIX_wfit   = "wfit"
COSMOFIT_CODE_WFIT = 'WFIT'
COSMOFIT_CODE_FIRECROWN = 'FIRECROWN'

# define columns for MERGE.LOG;  column 0 is always for STATE
COLNUM_WFIT_MERGE_DIROPT       = 1
COLNUM_WFIT_MERGE_COVOPT       = 2
COLNUM_WFIT_MERGE_WFITOPT      = 3
COLNUM_WFIT_MERGE_NDOF         = 4 
COLNUM_WFIT_MERGE_CPU          = 5

# define CONFIG key names
KEYNAME_WFITOPT_LIST  = [ "WFITOPT", "WFITOPTS" ] # allow either key
KEYNAME_COVOPT_LIST   = ["COVOPT", "COVOPTS"] 
KEYNAME_BLIND_DATA    = "BLIND_DATA"
KEYNAME_BLIND_SIM     = "BLIND_SIM"
KEYNAME_WFITAVG_LIST  = [ "WFITAVG", "WEIGHT_AVG" ]

BLIND_DATA_DEFAULT = True
BLIND_SIM_DEFAULT  = False

WFIT_SUMMARY_FILE     = "WFIT_SUMMARY.FITRES"
WFIT_SUMMARY_AVG_FILE = "WFIT_SUMMARY_AVG.FITRES"
WFIT_AVGTYPE_SINGLE   = "AVG_SINGLE"
WFIT_AVGTYPE_DIFF     = "AVG_DIFF"

# - - - - - - - - - - - - - - - - - - -  -
class cosmofit(Program):
    def __init__(self, config_yaml):
        CONFIG     = config_yaml['CONFIG']
        config_prep = {}
        if 'WFITOPT' in CONFIG:
            config_prep['COSMOFIT_CODE'] = COSMOFIT_CODE_WFIT
            config_prep['program'] = PROGRAM_NAME_WFIT
        elif 'FCOPT' in CONFIG :
            config_prep['COSMOFIT_CODE'] = COSMOFIT_CODE_FIRECROWN
            config_prep['program'] = PROGRAM_NAME_FIRECROWN
        else :
            msgerr=[]
            msgerr.append(f"Cannot determine program name")
            msgerr.append(f"Expecting WFITOPT or FCOPT in CONFIG block")
            log_assert(False,msgerr)            
        super().__init__(config_yaml, config_prep)
        return

        
    def set_output_dir_name(self):
        CONFIG     = self.config_yaml['CONFIG']
        input_file = self.config_yaml['args'].input_file  # for msgerr
        COSMOFIT_CODE = self.config_prep['COSMOFIT_CODE']
        msgerr     = []
        if 'OUTDIR' in CONFIG :
            output_dir_name = os.path.expandvars(CONFIG['OUTDIR'])
        else:
            msgerr.append(f"OUTDIR key missing in yaml-CONFIG")
            msgerr.append(f"Check {input_file}")
            log_assert(False,msgerr)
            
        return output_dir_name,SUBDIR_SCRIPTS_COSMOFIT
       
    
    def submit_prepare_driver(self):

        # Continue Firecrown development here.
        
        print("")

        # store input directories, and covsys_*
        self.cosmofit_prep_input_list()

        # store wfit options under WFITOPT key
        self.cosmofit_prep_wfitopt_list()

        # prepare blind flag for each inpdir based on data or sim
        self.cosmofit_prep_blind()

        # prepare index list for inpdir/covsys/wfitopt to simplify
        # 3D loops
        self.cosmofit_prep_index_lists()

        # prepare WFITAVG: mean and std err on mean
        self.cosmofit_prep_wfitavg()

        # copy input file to outdir
        input_file    = self.config_yaml['args'].input_file
        script_dir    = self.config_prep['script_dir']
        shutil.copy(input_file,script_dir)

        # end submit_prepare_driver

    def cosmofit_prep_input_list(self):

        # store user list of input directories (INPDIR key),
        # where each inpdir is an output of create_covariance.py.
        # For eac INPDIR, read and store list of covsys_[nnn].txt
        # files so that wfit runs for each systematic test. 

        msgerr = []
        input_file      = self.config_yaml['args'].input_file 
        CONFIG          = self.config_yaml['CONFIG']

        key = 'INPDIR'
        if key not in CONFIG:
            msgerr.append(f"Missing required {key} key in CONFIG block")
            msgerr.append(f"Check {input_file}")
            self.log_assert(False, msgerr)
        else:
            inpdir_list_orig = CONFIG[key] # orig list may include ENVs

        if inpdir_list_orig is None :
            msgerr.append(f"{key} is empty;")
            msgerr.append(f"Check {input_file}")
            self.log_assert(False, msgerr)

        # expand inpdir for wildcards and ENVs 
        inpdir_list = []
        for inpdir_orig in inpdir_list_orig:
            inpdir = os.path.expandvars(inpdir_orig)
            if '*' in inpdir_orig:
                tmp_list = sorted(glob.glob(inpdir))
            else:
                tmp_list = [ inpdir ]
            inpdir_list += tmp_list

        # - - - - -
        isdata_list        = []
        hd_file_list       = []
        covsys_file_list2d = [] # file list per inpdir
        covsys_num_list2d  = [] # cov index list per inpdir
        covinfo_list       = [] # list of yaml info per inpdir

        covsys_select_list = None  # Select All cov by default
        for key in KEYNAME_COVOPT_LIST:
            if key in CONFIG:
                covsys_select_list = CONFIG[key].split()

        for inpdir in inpdir_list:

            if not os.path.exists(inpdir) :
                msgerr.append(f"Cannot find input directory")
                msgerr.append(f"  {inpdir}")       
                self.log_assert(False, msgerr)            

            yaml_info, dict_info = \
                self.read_hd_info_file(inpdir,covsys_select_list)

            hd_file           = dict_info['hd_file'] 
            covsys_file_list  = dict_info['covsys_file_list']
            covsys_num_list   = dict_info['covsys_num_list']
            n_covsys          = len(covsys_file_list)
            isdata            = dict_info['isdata']   

            hd_file_list.append(hd_file)
            covsys_file_list2d.append(covsys_file_list)
            covsys_num_list2d.append(covsys_num_list)
            isdata_list.append(isdata)
            covinfo_list.append(yaml_info)

            #print(f" xxx {covsys_file_list}")
            #print(f" xxx yaml info = {yaml_info}")

            print(f" Found {inpdir} \n" \
                  f" \t with {n_covsys} covsys files and " \
                  f"ISDATA_REAL={isdata} ")


        #sys.exit(f"\n xxx inpdir_list = \n{inpdir_list} \n")
        #print(f" xxx covsys_list = {covsys_list} ")
        # - - - - - -
        #self.config_prep['inpdir_list_orig']  = inpdir_list_orig
        self.config_prep['inpdir_list']        = inpdir_list
        self.config_prep['n_inpdir']           = len(inpdir_list)

        self.config_prep['hd_file_list']       = hd_file_list
        self.config_prep['covsys_file_list2d'] = covsys_file_list2d
        self.config_prep['covsys_num_list2d']  = covsys_num_list2d
        self.config_prep['isdata_list']        = isdata_list
        self.config_prep['covinfo_list']       = covinfo_list

        self.wfit_error_check_input_list()

        #print(f" isdata_list = {isdata_list}")

        # end cosmofit_prep_input_list

    def cosmofit_prep_wfitavg(self):

        # parse WFITAVG key in CONFIG block.
        # Only do error checking at this point, no computation.
        # The avg calculations are done at the MERGE stage.

        CONFIG      = self.config_yaml['CONFIG']
        inpdir_list = self.config_prep['inpdir_list']

        KEYNAME_WFITAVG = self.get_keyname_wfit(KEYNAME_WFITAVG_LIST)
        if KEYNAME_WFITAVG is None: 
            return

        # - - - - - - - 
        # check that each wildcard corresponds to at least 1 input directory
        for wfitavg in CONFIG[KEYNAME_WFITAVG]:
            wfitavg_dirs = wfitavg.replace(' ','').split('-')
            for wildcard in wfitavg_dirs:
                matches = [f for f in inpdir_list if wildcard in f]
                if len(matches)<1:
                    msgerr=[f'Found no matches for widlcard {wildcard}']
                    self.log_assert(False, msgerr)

            if len(wfitavg_dirs) == 2:
                # check that both dirs have the same structure
                wildcard1 = wfitavg_dirs[0]
                wildcard2 = wfitavg_dirs[1]
                suffixes1 = [f.replace(wildcard1,'*') \
                             for f in inpdir_list if wildcard1 in f]
                suffixes2 = [f.replace(wildcard2,'*') \
                             for f in inpdir_list if wildcard2 in f]
                text_wildcard = f"wildcards {wildcard1} and {wildcard2}"

                if suffixes1 == suffixes2:
                    print(f"Found matching number/names of dirs for " \
                          f"{text_wildcard}")
                else:
                    msgerr = []
                    len1 = len(suffixes1) ; len2=len(suffixes2)
                    msgerr.append(f"Mis-match for {text_wildcard}.")
                    msgerr.append(f"Found {len1} dirs for wildcard={wildcard1}")
                    msgerr.append(f"  {suffixes1}")
                    msgerr.append(f"Found {len2} dirs for wildcard={wildcard2}")
                    msgerr.append(f"  {suffixes2}")
                    self.log_assert(False, msgerr)

        # - - - - - 
        return
        # end  cosmofit_prep_wfitavg

    def wfit_error_check_input_list(self):

        # loop over each inpdir and abort on problems such as
        # non-existing inpdir, n_covsys=0, etc ...
        # Print all ERRORS before aborting.

        #inpdir_list_orig = self.config_prep['inpdir_list_orig']
        inpdir_list        = self.config_prep['inpdir_list'] 
        hd_file_list       = self.config_prep['hd_file_list']
        covsys_file_list2d = self.config_prep['covsys_file_list2d']
        nerr = 0
        msgerr = []

        for inpdir, hd_base, covsys_file_list in \
            zip(inpdir_list, hd_file_list, covsys_file_list2d):

            hd_file    = f"{inpdir}/{hd_base}"
            n_covsys   = len(covsys_file_list)
            if n_covsys == 0 :            
                nerr += 1
                msgerr.append(f"ERROR: cannot find covsys files in")
                msgerr.append(f"   {inpdir}")  

            if not os.path.exists(hd_file):
                nerr += 1
                msgerr.append(f"ERROR: cannot find expected HD file:")
                msgerr.append(f"   {hd_file}")       
        
        # - - - - - - -
        if nerr > 0 :
            msgerr.append(f"Found {nerr} errors with INPDIR list.")
            msgerr.append(f"See {nerr} ERROR messages above.")
            self.log_assert(False, msgerr)

        # end wfit_error_check_input_list

    def read_hd_info_file(self, inpdir, covsys_select_list):

        # Ceated Oct 13 2022 by RK
        # Read Hubble-diagram (HD) info file created by create_covariance.py.
        # Translate yaml contents into a more practical dictionary.
        #
        # Inputs:
        #  inpdir: directory containing INFO.YML
        #  covsys_select_list: list of covsys labels to select for 
        #       separate wfit task;
        #       e.g., covsys_select_list = [ 'ALL', 'ZP' ]
        #       If covsys_select_list is None, then each covsys is used.
        #
        # Output dictionary includes
        #  + base name of HD file
        #  + isdata flag (True for real data, False for sim)
        #  + list of covsys num-indices
        #  + list of covsys labels
        #  + list of covsys base file names
        #  
        
        # define yaml keys written by create_covariance.py
        INFO_FILENAME          = "INFO.YML" # read this from inpdir

        INFO_KEYNAME_HD        = "HD"
        INFO_KEYNAME_COVOPTS   = "COVOPTS"
        INFO_KEYNAME_ISDATA    = "ISDATA_REAL"     # key in cov info file
        HD_BASENAME_LEGACY     = "hubble_diagram.txt" 

        yaml_file = f"{inpdir}/{INFO_FILENAME}"
        yaml_info = util.extract_yaml(yaml_file, None, None )

        # - - - - - -
        # load separate info dictionary to store info in more convenient way

        # get name of Hubble diagram file
        key     = INFO_KEYNAME_HD
        if key in yaml_info:
            hd_file = yaml_info[key]  # refac read from INFO file
        else:
            hd_file = HD_BASENAME_LEGACY  # legacy hard wite
    

        # read flag indicating real data
        key     = INFO_KEYNAME_ISDATA
        if key in yaml_info:
            isdata = (yaml_info[key] > 0)
        else:
            isdata  = False  # default in dase key is missing
            
        # - - - -
        COVOPTS = yaml_info[INFO_KEYNAME_COVOPTS]
        covsys_num_list   = []
        covsys_label_list = []
        covsys_file_list  = []

        for covnum, covinfo in COVOPTS.items():
            covinfo_split = covinfo.split()

            covsys_label  = covinfo_split[0]
            if len(covinfo_split) > 1:
                covsys_file = covinfo_split[1] # refactored
            else:
                covsys_file = f"covsys_{covnum:03d}.txt.gz"  # legacy hard wire

            covsys_num_list.append(covnum)
            covsys_label_list.append(covsys_label)
            covsys_file_list.append(covsys_file)

        # - - - - -
        # check optional subset of covsys options to store
        if covsys_select_list is not None :
            COVOPTS_DICT             = yaml_info["COVOPTS"].copy()
            yaml_info_select         = yaml_info.copy()
            covsys_file_list_select  = []
            covsys_num_list_select   = []
            for covsys_num_tmp, covsys_file_tmp in \
                zip(COVOPTS_DICT, covsys_file_list):
                covsys_name_tmp = COVOPTS_DICT[covsys_num_tmp]

                if covsys_name_tmp in covsys_select_list:
                    covsys_file_list_select.append(covsys_file_tmp)
                    covsys_num_list_select.append(covsys_num_tmp)
                else:
                    yaml_info_select["COVOPTS"].pop(covsys_num_tmp)

            # update lists to include only the user-requsted subset
            covsys_file_list = covsys_file_list_select
            covsys_num_list  = covsys_num_list_select
            yaml_info        = yaml_info_select

        # - - - - - - - - - - 
        # store dictionary
        dict_info = {
            'hd_file'           : hd_file,
            'covsys_num_list'   : covsys_num_list,
            'covsys_label_list' : covsys_label_list,
            'covsys_file_list'  : covsys_file_list,
            'isdata'            : isdata            
        }
        
        #print(f"\n xxx dict_info = \n{dict_info}\n")

        return yaml_info, dict_info
        # end read_hd_info_file


    def cosmofit_prep_wfitopt_list(self):

        msgerr = []
        input_file      = self.config_yaml['args'].input_file 
        CONFIG          = self.config_yaml['CONFIG']
        output_dir      = self.config_prep['output_dir']
        wfitopt_rows    = None
        key_found       = None

        for key in KEYNAME_WFITOPT_LIST:
            if key in CONFIG:
                wfitopt_rows = CONFIG[key]
                key_found    = key

        if wfitopt_rows is None:
            msgerr.append(f"Missing required CONFIG key ")
            msgerr.append(f"   {KEYNAME_WFITOPT_LIST} ")
            msgerr.append(f"One of these keys must be in {input_file}")
            self.log_assert(False, msgerr)

        wfitopt_dict = util.prep_jobopt_list(wfitopt_rows, key_found, 0, None)

        n_wfitopt          = wfitopt_dict['n_jobopt']
        wfitopt_arg_list   = wfitopt_dict['jobopt_arg_list']
        wfitopt_num_list   = wfitopt_dict['jobopt_num_list']
        wfitopt_label_list = wfitopt_dict['jobopt_label_list']

        logging.info(f"\n Store {n_wfitopt} wfit options from " \
                     f"{key_found} keys" )

        self.config_prep['n_wfitopt']          = n_wfitopt
        self.config_prep['wfitopt_arg_list']   = wfitopt_arg_list
        self.config_prep['wfitopt_num_list']   = wfitopt_num_list
        self.config_prep['wfitopt_label_list'] = wfitopt_label_list

        # check for global wfitopt
        wfitopt_global = ""
        for key_base in KEYNAME_WFITOPT_LIST:
            key   = f"{key_base}_GLOBAL"            
            if key in CONFIG :
                wfitopt_global = CONFIG[key]
        
        self.config_prep['wfitopt_global'] = wfitopt_global

        # - - - - - 
        # check for wa in fit
        use_wa          = False
        outdir_chi2grid = None
        tmp_list = wfitopt_arg_list + [ wfitopt_global ]
        for tmp in tmp_list:
            if '-wa'               in tmp : 
                use_wa = True
            if '-outfile_chi2grid' in tmp : 
                outdir_chi2grid = f"{output_dir}/CHI2GRID"
                if not os.path.exists(outdir_chi2grid):
                    os.mkdir(outdir_chi2grid)

        self.config_prep['use_wa'] = use_wa
        self.config_prep['outdir_chi2grid'] = outdir_chi2grid

        return

        # end cosmofit_prep_wfitopt_list

    def cosmofit_prep_blind(self):

        CONFIG   = self.config_yaml['CONFIG']

        # if no user-override for blind optoin, set CONFIG to default

        # check user override for data
        key = KEYNAME_BLIND_DATA
        if key not in CONFIG:  CONFIG[key] = BLIND_DATA_DEFAULT
        blind_data = CONFIG[key]

        # ... and for sim
        key = KEYNAME_BLIND_SIM
        if key not in CONFIG:  CONFIG[key] = BLIND_SIM_DEFAULT
        blind_sim = CONFIG[key]

        logging.info(f" ")
        logging.info(f"\t BLIND DATA: {blind_data}")
        logging.info(f"\t BLIND SIM:  {blind_sim}")
        logging.info(f" ")
        
        # abort if any WFITOPT has -blind ... to avoid interference
        # with BLIND_DATA and BLIND_SIM yaml flags
        
        wfitopt_list     = self.config_prep['wfitopt_arg_list']
        wfitopt_global   = self.config_prep['wfitopt_global']
        tmp_list = wfitopt_list + [ wfitopt_global ]
        for wfitopt in tmp_list:
            if "-blind" in wfitopt :
                msgerr=[]
                msgerr.append(f"Cannot use -blind arg in WFITOPT.")
                msgerr.append(f"Control blinding in CONFIG with")
                msgerr.append(f"  {KEYNAME_BLIND_DATA} and {KEYNAME_BLIND_SIM} keys")
                self.log_assert(False, msgerr)

        # check isdata flag per inpdir
        inpdir_list      = self.config_prep['inpdir_list']
        isdata_list      = self.config_prep['isdata_list']
        arg_blind_list   = [] # set to either "-blind" or ""
        for inpdir,isdata in zip(inpdir_list,isdata_list):
            issim = not isdata
            arg_blind = "-blind"
            if isdata and not blind_data : arg_blind = ""
            if issim  and not blind_sim  : arg_blind = ""
            arg_blind_list.append(arg_blind)

        self.config_prep['arg_blind_list'] = arg_blind_list

        # end cosmofit_prep_blind

    def cosmofit_prep_index_lists(self):

        # prepare internal index lists for efficient looping
        CONFIG             = self.config_yaml['CONFIG']
        inpdir_list        = self.config_prep['inpdir_list']  
        covsys_file_list2d = self.config_prep['covsys_file_list2d']
        covsys_num_list2d  = self.config_prep['covsys_num_list2d']
        wfitopt_list       = self.config_prep['wfitopt_arg_list']

        n_inpdir         = self.config_prep['n_inpdir']
        n_wfitopt        = self.config_prep['n_wfitopt']        

        # count total number of jobs, and beware that number
        # of covsys can be different in each inpdir.
        n_job_tot = 0
        idir_list3 = [];  ifit_list3 = [];   icov_list3 = []

        for idir in range(0,n_inpdir):
            covsys_file_list = covsys_file_list2d[idir]
            n_covsys         = len(covsys_file_list)
            for icov in range(0,n_covsys):
                for ifit in range(0,n_wfitopt):
                    n_job_tot += 1
                    idir_list3.append(idir)
                    icov_list3.append(icov)
                    ifit_list3.append(ifit)

        self.config_prep['n_job_split'] = 1
        self.config_prep['n_job_tot']   = n_job_tot
        self.config_prep['n_done_tot']  = n_job_tot

        self.config_prep['idir_list3']  = idir_list3
        self.config_prep['icov_list3']  = icov_list3
        self.config_prep['ifit_list3']  = ifit_list3

        # end cosmofit_prep_index_lists

    def write_command_file(self, icpu, f):

        input_file         = self.config_yaml['args'].input_file 
        inpdir_list        = self.config_prep['inpdir_list']  
        covsys_file_list2d = self.config_prep['covsys_file_list2d']
        wfitopt_list       = self.config_prep['wfitopt_arg_list']
        n_core             = self.config_prep['n_core']

        idir_list3 = self.config_prep['idir_list3'] 
        icov_list3 = self.config_prep['icov_list3']
        ifit_list3 = self.config_prep['ifit_list3']
        
        n_job_cpu   = 0
        n_job_local = 0

        for idir,icov,ifit in zip(idir_list3,icov_list3,ifit_list3):

            n_job_local += 1
            index_dict = \
                { 'idir':idir, 'ifit':ifit, 'icov':icov, 'icpu':icpu }

            if ( (n_job_local-1) % n_core ) != icpu : continue

            n_job_cpu += 1
            job_info_wfit   = self.prep_JOB_INFO_wfit(index_dict)
            util.write_job_info(f, job_info_wfit, icpu)

            job_info_merge = \
                self.prep_JOB_INFO_merge(icpu,n_job_local,False) 
            util.write_jobmerge_info(f, job_info_merge, icpu)

        return n_job_cpu

        # end write_command_file

    def prep_JOB_INFO_wfit(self,index_dict):

        idir = index_dict['idir']
        ifit = index_dict['ifit']
        icov = index_dict['icov']

        kill_on_fail = self.config_yaml['args'].kill_on_fail
        program      = self.config_prep['program']
        output_dir   = self.config_prep['output_dir']
        script_dir   = self.config_prep['script_dir']

        inpdir       = self.config_prep['inpdir_list'][idir]
        arg_blind    = self.config_prep['arg_blind_list'][idir]
        arg_string   = self.config_prep['wfitopt_arg_list'][ifit]
        arg_global   = self.config_prep['wfitopt_global']
        cov_basename = self.config_prep['covsys_file_list2d'][idir][icov]
        hd_basename  = self.config_prep['hd_file_list'][idir]
        outdir_chi2grid = self.config_prep['outdir_chi2grid']

        prefix = self.wfit_num_string(idir,icov,ifit)

        covsys_file   = f"{inpdir}/{cov_basename}"
        hd_file       = f"{inpdir}/{hd_basename}"
        log_file      = f"{prefix}.LOG" 
        done_file     = f"{prefix}.DONE"
        all_done_file = f"{output_dir}/{DEFAULT_DONE_FILE}"


        # start with user-defined args from WFITOPT[_GLOBAL] key
        arg_list =  [ arg_string ]
        if len(arg_global) > 0: arg_list.append(arg_global)

        if outdir_chi2grid is not None :
            outfile  = f"{outdir_chi2grid}/{prefix}_CHI2GRID.DAT"
            arg_list = util.replace_arg(arg_list,"-outfile_chi2grid",outfile)
            #print(f" xxx arg_list -> {arg_list}")

        # define covsys file from create_cov
        arg_list.append(f"-mucov_file {covsys_file}")

        # tack on blind arg
        arg_list.append(arg_blind)

        # define output YAML file to be parsed by submit-merge process
        arg_list.append(f"-cospar_yaml {prefix}.YAML")

        JOB_INFO = {}
        JOB_INFO['program']       = program
        JOB_INFO['input_file']    = hd_file
        JOB_INFO['job_dir']       = script_dir
        JOB_INFO['log_file']      = f"{log_file}"
        JOB_INFO['done_file']     = f"{done_file}"
        JOB_INFO['all_done_file'] = f"{all_done_file}"
        JOB_INFO['kill_on_fail']  = kill_on_fail
        JOB_INFO['arg_list']      = arg_list
  
        return JOB_INFO

        # end prep_JOB_INFO_wfit

    def wfit_num_string(self,idir,icov,ifit):

        if idir >= 0 and icov < 0 and ifit < 0:
            string = f"DIROPT{idir:05d}"

        elif icov >= 0 and idir < 0 and ifit < 0 :
            string = f"COVOPT{icov:03d}"
        elif ifit >= 0 and idir < 0 and icov < 0 :
            string = f"WFITOPT{ifit:03d}"

        elif idir < 0 and ifit>=0 and icov>=0 :
            string = f"COVOPT{icov:03d}_WFITOPT{ifit:03d}"

        elif idir >= 0 and ifit>=0 and icov>=0 :
            string = f"DIROPT{idir:05d}_COVOPT{icov:03d}_WFITOPT{ifit:03d}"
        else:
            string = "ERROR"

        return string
        # end wfit_num_string

    def cosmofit_prefix(self,row):
        # parse input row passed from MERGE.LOG and construct
        # prefix for output files
        dirnum  = row[COLNUM_WFIT_MERGE_DIROPT]  # e.g, DIROPT00003
        covnum  = row[COLNUM_WFIT_MERGE_COVOPT]  # e.g, COVOPT002
        wfitnum = row[COLNUM_WFIT_MERGE_WFITOPT] # e.g  WFITOPT001
        prefix = f"{dirnum}_{covnum}_{wfitnum}"
        return prefix

        # end cosmofit_prefix

    def append_info_file(self,f):
        # append info to SUBMIT.INFO file

        CONFIG             = self.config_yaml['CONFIG']
        n_wfitopt          = self.config_prep['n_wfitopt']
        wfitopt_arg_list   = self.config_prep['wfitopt_arg_list'] 
        wfitopt_num_list   = self.config_prep['wfitopt_num_list']
        wfitopt_label_list = self.config_prep['wfitopt_label_list']
        wfitopt_global     = self.config_prep['wfitopt_global']
        inpdir_list        = self.config_prep['inpdir_list']
        covinfo_list       = self.config_prep['covinfo_list']
        use_wa             = self.config_prep['use_wa']

        blind_data   = CONFIG[KEYNAME_BLIND_DATA] # T or F
        blind_sim    = CONFIG[KEYNAME_BLIND_SIM] # T or F

        f.write(f"\n# wfit info\n")

        f.write(f"JOBFILE_WILDCARD:  'DIR*COVOPT*WFITOPT*' \n")

        idir=-1
        f.write("\n")
        f.write(f"INPDIR_LIST: \n")
        for inpdir, covinfo in zip(inpdir_list, covinfo_list):
            idir += 1
            diropt_num  = self.wfit_num_string(idir,-1,-1)
            row         = [ diropt_num, inpdir ]
            f.write(f"  {diropt_num}: {inpdir} \n")
            #f.write(f"  - {row} \n")

            COVOPTS      = covinfo['COVOPTS']
            COVOPTS_keys = list(COVOPTS.keys())

            n_covopt = len(COVOPTS)
            f.write(f"  COVOPTS({diropt_num}): \n")
            for icov in range(0,n_covopt):
                covopt_num  = self.wfit_num_string(-1,icov,-1) 
                covindx     = COVOPTS_keys[icov]  # original index
                covopt      = COVOPTS[covindx]    # label  covFile
                covopt_label = covopt.split()[0]  # just the label
                f.write(f"    {covopt_num}: {covopt_label} \n")

        f.write("\n")
        f.write(f"{KEYNAME_BLIND_DATA}:  {blind_data}\n")
        f.write(f"{KEYNAME_BLIND_SIM}:   {blind_sim}\n")

        f.write("\n")
        f.write(f"N_WFITOPT:         {n_wfitopt}      " \
                f"# number of wfit options\n")
        f.write(f"USE_wa:        {use_wa}   " \
                f"# T if any WFITOPT uses waw0CDM model\n")

        f.write("\n")
        f.write("WFITOPT_LIST:  " \
                "# 'WFITOPTNUM'  'user_label'  'user_args'\n")
        for num,arg,label in zip(wfitopt_num_list, wfitopt_arg_list,
                                 wfitopt_label_list):
            row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")
        f.write(f"WFITOPT_GLOBAL: {wfitopt_global} \n")

        # end append_info_file  

    def create_merge_table(self,f):

        # create merge table with NDOF=0; later this table
        # gets updated as jobs finish.
        idir_list3      = self.config_prep['idir_list3']
        icov_list3      = self.config_prep['icov_list3'] 
        ifit_list3      = self.config_prep['ifit_list3']
        
        # create only MERGE table ... no need for SPLIT table
        header_line_merge = \
                f" STATE  DIROPT  COVOPT  WFITOPT  NDOF CPU "
        # xxx f" SPLITRAN"

        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state

        for idir,icov,ifit in \
            zip(idir_list3,icov_list3,ifit_list3):

            diropt_num  = self.wfit_num_string(idir,-1,-1)
            covopt_num  = self.wfit_num_string(-1,icov,-1)
            wfitopt_num = self.wfit_num_string(-1,-1,ifit)

            # ROW here is fragile in case columns are changed
            ROW_MERGE = []
            ROW_MERGE.append(STATE)
            ROW_MERGE.append(diropt_num)
            ROW_MERGE.append(covopt_num)
            ROW_MERGE.append(wfitopt_num)
            ROW_MERGE.append(0)    # Ndof
            ROW_MERGE.append(0.0)  # CPU
            INFO_MERGE['row_list'].append(ROW_MERGE)  

        # - - - - -
        util.write_merge_file(f, INFO_MERGE, [] ) 
            
        # end create_merge_file

    def merge_config_prep(self,output_dir):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        
    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']

        COLNUM_STATE     = COLNUM_MERGE_STATE
        COLNUM_DIROPT    = COLNUM_WFIT_MERGE_DIROPT
        COLNUM_COVOPT    = COLNUM_WFIT_MERGE_COVOPT  
        COLNUM_WFITOPT   = COLNUM_WFIT_MERGE_WFITOPT  
        COLNUM_NDOF      = COLNUM_WFIT_MERGE_NDOF
        COLNUM_CPU       = COLNUM_WFIT_MERGE_CPU
        NROW_DUMP   = 0

        key_ndof, key_ndof_sum, key_ndof_list = \
                self.keynames_for_job_stats('Ndof')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')

        key_list = [ key_ndof, key_cpu ] 

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
            prefix      = self.cosmofit_prefix(row) 
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
                    
                    wfit_stats = self.get_job_stats(script_dir,
                                                    log_list, 
                                                    yaml_list, 
                                                    key_list)
                    
                    # check for failures in snlc_fit jobs.
                    nfail = wfit_stats['nfail']
                    if nfail > 0 :  NEW_STATE = SUBMIT_STATE_FAIL
                 
                    row[COLNUM_STATE]     = NEW_STATE
                    row[COLNUM_NDOF]      = wfit_stats[key_ndof_sum]
                    row[COLNUM_CPU]       = wfit_stats[key_cpu_sum]
                    
                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             # assume nevt changes

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

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        submit_info_yaml = self.config_prep['submit_info_yaml']
        output_dir       = self.config_prep['output_dir']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        row     = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]

        # end merge_job_wrapup

    def merge_cleanup_final(self):
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_subdir    = SUBDIR_SCRIPTS_WFIT

        self.make_wfit_summary()
        self.make_wfitavg_summary()

        logging.info(f"  wfit cleanup: compress {JOB_SUFFIX_TAR_LIST}")
        for suffix in JOB_SUFFIX_TAR_LIST :
            wildcard = (f"{jobfile_wildcard}*.{suffix}") 
            util.compress_files(+1, script_dir, wildcard, suffix, "" )

        logging.info("")

        self.merge_cleanup_script_dir()

        # end merge_cleanup_final

    def make_wfit_summary(self):

        CONFIG           = self.config_yaml['CONFIG']
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wa           = submit_info_yaml['USE_wa']
        INPDIR_LIST      = submit_info_yaml['INPDIR_LIST']
        WFITOPT_LIST     = submit_info_yaml['WFITOPT_LIST']

        SUMMARYF_FILE    = f"{output_dir}/{WFIT_SUMMARY_FILE}"
        logging.info(f"\t Writing wfit summary to {WFIT_SUMMARY_FILE}")
        out_lines_list = []

        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        MERGE_INFO_CONTENTS,comment_lines = \
                util.read_merge_file(MERGE_LOG_PATHFILE)

        dirnum_last = "xyz"
        wfit_summary_table = {}
        nrow = 0 ; nrow_warn = 0
        for row in MERGE_INFO_CONTENTS[TABLE_MERGE]:
            nrow += 1
            dirnum     = row[COLNUM_WFIT_MERGE_DIROPT][-5:] # e.g., DIROPT00000
            covnum     = row[COLNUM_WFIT_MERGE_COVOPT][-3:] # e.g., COVOPT001
            wfitnum    = row[COLNUM_WFIT_MERGE_WFITOPT][-3:] # idem
            prefix     = self.cosmofit_prefix(row)
            YAML_FILE  = f"{script_dir}/{prefix}.YAML"

            wfit_yaml        = util.extract_yaml(YAML_FILE, None, None )
            wfit_values_dict = util.get_wfit_values(wfit_yaml)

            w        = wfit_values_dict['w']  
            w_sig    = wfit_values_dict['w_sig']
            omm      = wfit_values_dict['omm']  
            omm_sig  = wfit_values_dict['omm_sig']
            rho_womm = wfit_values_dict['rho_womm']
            chi2     = wfit_values_dict['chi2'] 
            sigint   = wfit_values_dict['sigint']
            blind    = wfit_values_dict['blind']
            nwarn    = wfit_values_dict['nwarn']
            if nwarn > 0 : nrow_warn += 1

            # extract user labels for cov and wfit
            str_diropt    = 'DIROPT' + dirnum
            str_covopt    = 'COVOPT' + covnum
            dir_name      = INPDIR_LIST[str_diropt]  # create_cov dir
            covopt_dict   = INPDIR_LIST[f'COVOPTS({str_diropt})']
            covopt_label  = covopt_dict[str_covopt]
            wfitopt_label = WFITOPT_LIST[int(wfitnum)][1]
            
            if wfitopt_label == "None" : 
                wfitopt_label = "NoLabel"

            if use_wa:
                wa      = wfit_values_dict['wa']    
                wa_sig  = wfit_values_dict['wa_sig']
                FoM     = wfit_values_dict['FoM']
                if 'Rho' in wfit_values_dict:
                    rho_w0wa = wfit_values_dict['Rho']  # legacy name
                else:
                    rho_w0wa = wfit_values_dict['rho_w0wa']  # 9.27.2022      
            else:
                wa       = 0
                wa_sig   = 0
                FoM      = 0
                rho_w0wa = 0

            # load table for mean and std err on mean (table not used here)
            local_dict = {'dirnum': dirnum, 
                          'covnum': covnum, 
                          'wfitnum': wfitnum, 
                          'w':w, 'w_sig':w_sig, 
                          'omm':omm, 'omm_sig':omm_sig, 
                          'wa':wa, 'wa_sig': wa_sig,
                          'rho_w0omm' : rho_womm,
                          'rho_w0wa'  : rho_w0wa,
                          'FoM':FoM, 
                          'covopt_label':covopt_label,
                          'wfitopt_label':wfitopt_label,
                          'nwarn' : nwarn  # R.Kessler Feb 23 2022
            }

            unique_key = dir_name + '_' + covnum + '_' + wfitnum 
            wfit_summary_table[unique_key] = local_dict
            
            if nrow == 1:
                out_lines_list += \
                    self.write_wfit_summary_header(wfit_values_dict)
                
            if dirnum != dirnum_last:
                out_lines_list.append(f"#\n# {str_diropt}={dir_name}")

            str_nums    = f"{dirnum} {covnum} {wfitnum}"
            if use_wa : 
                str_results  = f"{w:.5f} {w_sig:.5f} "
                str_results += f"{wa:7.5f} {wa_sig:7.5f} "
                str_results += f"{FoM:5.1f} {rho_w0wa:6.3f} "
                str_results += f"{omm:.3f} {omm_sig:.3f}  "
            else:
                str_results  = f"{w:.6f} {w_sig:.6f}  "
                str_results += f"{omm:.5f} {omm_sig:.5f}  "
                str_results += f"{rho_womm:6.3f} "

            str_misc    = f"{chi2:6.1f} {blind} {nwarn} "
            str_labels  = f"{covopt_label:<10} {wfitopt_label}"
            line = f"ROW: {nrow:3d} {str_nums} {str_results}" \
                                  f"{str_misc} {str_labels}"
            out_lines_list.append(f"{line}")
            dirnum_last = dirnum

        # - - - - - - - -
        with open(SUMMARYF_FILE,"w") as f:
            if nrow_warn > 0:
                text_warn = f"WARNING: {nrow_warn} of {nrow} fits have warnings;"\
                            f" check nwarn column."
                f.write(f"# {text_warn}\n\n")
            for line in out_lines_list:  f.write(f"{line}\n")
        
        self.config_prep['wfit_summary_table'] = wfit_summary_table
        # end make_wfit_summary


    def write_wfit_summary_header(self,wfit_values_dict):

        # return lines with header info and VARNAMES for wfit-summary file

        submit_info_yaml = self.config_prep['submit_info_yaml']
        use_wa           = submit_info_yaml['USE_wa']

        varnames_w  = "w w_sig"
        varnames_om = "omm omm_sig rho_womm"
        if use_wa: 
            varnames_w  = "w0 w0_sig wa wa_sig FoM rho_w0wa"
            varnames_om = "omm omm_sig"

        VARNAMES_STRING = \
            f"ROW  iDIR iCOV iWFIT {varnames_w} "  \
            f"{varnames_om} chi2 blind nwarn COVOPT WFITOPT"

        w_ran   = int(wfit_values_dict['w_ran']) 
        wa_ran  = int(wfit_values_dict['wa_ran'])
        omm_ran = int(wfit_values_dict['omm_ran'])

        lines_list = []
        if w_ran > 0 : 
            lines_list.append(f"#  blind=1 ==> w0,wa,omm += " \
                              f"sin({w_ran},{wa_ran},{omm_ran}) ")
            lines_list.append("")
                             
        lines_list.append(f"VARNAMES: {VARNAMES_STRING} ")
        return lines_list

        # write_wfit_summary_header

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file  
        submit_info_yaml = self.config_prep['submit_info_yaml']
        info_lines = []
        return info_lines
        # end get_misc_merge_info 

    def get_merge_COLNUM_CPU(self):
        return -9  # there is no CPU column

    def merge_reset(self,output_dir):

        # unpack things in merge_cleanup_final, but in reverse order
        submit_info_yaml = self.config_prep['submit_info_yaml']
        jobfile_wildcard = submit_info_yaml['JOBFILE_WILDCARD']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        script_subdir    = SUBDIR_SCRIPTS_WFIT
        fnam = "merge_reset"

        logging.info(f"   {fnam}: reset STATE and NEVT in {MERGE_LOG_FILE}")
        MERGE_LOG_PATHFILE = f"{output_dir}/{MERGE_LOG_FILE}"
        colnum_zero_list = [ COLNUM_WFIT_MERGE_NDOF ]
        util.merge_table_reset(MERGE_LOG_PATHFILE, TABLE_MERGE,  \
                               COLNUM_MERGE_STATE, colnum_zero_list)

        util.untar_script_dir(script_dir)

        # xxx logging.info(f"  {fnam}: uncompress {script_subdir}/")
        # xxx util.compress_subdir(-1, f"{output_dir}/{script_subdir}" )

        return
        # end merge_reset

    def get_keyname_wfit(self, KEYNAME_LIST):

        # for list of possible key names in KEYNAME_LIST,
        # return the keyname that exists in CONFIG block.

        CONFIG  = self.config_yaml['CONFIG']
        keyname = None
        for key in KEYNAME_LIST:
            if key in CONFIG: 
                keyname = key
        return keyname
        # end get_keyname_wfit

    def make_wfitavg_lists(self):

        CONFIG           = self.config_yaml['CONFIG']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        INPDIR_LIST      = submit_info_yaml['INPDIR_LIST']

        inpdirs_full_paths = [INPDIR_LIST[k] \
                            for k in INPDIR_LIST.keys() if k.startswith('DIROPT')]

        KEYNAME_WFITAVG = self.get_keyname_wfit(KEYNAME_WFITAVG_LIST)

        wfitavg_list = {}
        for wfitavg in CONFIG[KEYNAME_WFITAVG]:
            wfitavg_dirs = wfitavg.replace(' ','').split('-')
            wildcard = wfitavg_dirs[0]
            dirslist_fullpath = [f for f in inpdirs_full_paths if wildcard in f]
            wfitavg_list[wfitavg] = {
                    'avg_type' : WFIT_AVGTYPE_SINGLE,
                    'wildcard1' : wildcard,
                    'wildcard2' : None,
                    'dirslist_fullpath1' : dirslist_fullpath,
                    'dirslist_fullpath2': None
            }
            if len(wfitavg_dirs)==2: 
                # check first that dirs match
                wildcard2 = wfitavg_dirs[1]
                dirslist_fullpath2 = [f for f in inpdirs_full_paths if wildcard2 in f]
                wfitavg_list[wfitavg]['avg_type'] = WFIT_AVGTYPE_DIFF
                wfitavg_list[wfitavg]['wildcard2'] = wildcard2
                wfitavg_list[wfitavg]['dirslist_fullpath2'] = dirslist_fullpath2

        self.config_prep['wfitavg_list'] = wfitavg_list

        # end make_wfitavg_lists                                                                                                                              

    def compute_average(self, fit_list):
        fit_array = np.array(fit_list)
        
        mean        = np.mean(fit_array)
        std_of_mean = np.std(fit_array)/np.sqrt(len(fit_list))

        if np.isnan(mean):   # RK
            mean = 0.0 
            std_of_mean = 0.0

        return mean, std_of_mean
        # end compute_average
        
    def make_wfitavg_summary(self):

        # Driver utility to compute means and std err on mean among directories
        # See MEAN_STDERRMEAN key in input CONFIG file

        CONFIG           = self.config_yaml['CONFIG']
        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_wa           = submit_info_yaml['USE_wa']
        INPDIR_LIST      = submit_info_yaml['INPDIR_LIST']
        WFITOPT_LIST     = submit_info_yaml['WFITOPT_LIST']
        wfit_summary_table  = self.config_prep['wfit_summary_table']

        KEYNAME_WFITAVG = self.get_keyname_wfit(KEYNAME_WFITAVG_LIST)
        if KEYNAME_WFITAVG is None: return

        # load lists of files needed to comput avgs
        self.make_wfitavg_lists()
        wfitavg_list  = self.config_prep['wfitavg_list']
        
        logging.info(f"\t Writing means summary to {WFIT_SUMMARY_AVG_FILE}")

        AVG_FILE    = f"{output_dir}/{WFIT_SUMMARY_AVG_FILE}"
        f = open(AVG_FILE,"w")
        nrow = 0
        
        if use_wa: VARNAME_FOM = "<w_sig> <w_sig>_sig FoM FoM_sig"
        else: VARNAME_FOM = "<w_sig> <w_sig>_sig"

        VARNAMES_STRING = \
                          f"ROW  iCOV iWFIT <w> <w>_sig   <wa> <wa>_sig  "  \
                          f"<omm> <omm>_sig {VARNAME_FOM} N_DIRs COVOPT WFITOPT"
        f.write(f"VARNAMES: {VARNAMES_STRING} \n")
        
        avg_comment_dict = {
            WFIT_AVGTYPE_SINGLE : ' (mean and std err on fitted values)',
            WFIT_AVGTYPE_DIFF   : ' (mean and std err on differences in fitted values)'
        }

        for wfitavg in wfitavg_list:
            # use the first directory in the list to find the set of 
            # unique covopts and unique wfitopts
            dir_0 = wfitavg_list[wfitavg]['dirslist_fullpath1'][0]
            unique_matching_covopts = np.unique([f.replace(dir_0,'')[1:4] \
                                                 for f in wfit_summary_table.keys()\
                                                 if f[:-8]==dir_0])
            unique_matching_wfitopts = np.unique([f.replace(dir_0,'')[5:] \
                                                  for f in wfit_summary_table.keys()\
                                                  if f[:-8]==dir_0])
            f.write(f"#\n# Mean and std err on mean for option: " \
                    f"{wfitavg} {avg_comment_dict[wfitavg_list[wfitavg]['avg_type']]}\n")

            for covnum in unique_matching_covopts:
                for wfitnum in unique_matching_wfitopts:
                    omm_list = []; w_list = []; wa_list = []; wsig_list = []; FoM_list = []
                    for dir_ in wfitavg_list[wfitavg]['dirslist_fullpath1']:
                        unique_key = dir_+'_%s_%s'%(covnum,wfitnum)
                        omm_list.append(wfit_summary_table[unique_key]['omm'])
                        w_list.append(wfit_summary_table[unique_key]['w'])
                        wa_list.append(wfit_summary_table[unique_key]['wa'])
                        wsig_list.append(wfit_summary_table[unique_key]['w_sig'])
                        FoM_list.append(wfit_summary_table[unique_key]['FoM'])
                    covopt_label  = wfit_summary_table[unique_key]['covopt_label']
                    wfitopt_label = wfit_summary_table[unique_key]['wfitopt_label']

                    if wfitavg_list[wfitavg]['dirslist_fullpath2'] is not None:
                        omm_list2 = []; w_list2 = []; wa_list2 = []; wsig_list2 = []; FoM_list2 = []
                        for dir2_ in wfitavg_list[wfitavg]['dirslist_fullpath2']:
                            unique_key = dir2_+'_%s_%s'%(covnum,wfitnum)
                            omm_list2.append(wfit_summary_table[unique_key]['omm'])
                            w_list2.append(wfit_summary_table[unique_key]['w'])
                            wa_list2.append(wfit_summary_table[unique_key]['wa'])
                            wsig_list2.append(wfit_summary_table[unique_key]['w_sig'])
                            FoM_list2.append(wfit_summary_table[unique_key]['FoM'])
                    else: # if theres no second set of dirs, it mean this is not a difference so just set x_list2 to zero 
                        zero_list2 = np.zeros(len(wfitavg_list[wfitavg]['dirslist_fullpath1']))
                        omm_list2 = zero_list2
                        w_list2 = zero_list2
                        wa_list2 = zero_list2
                        wsig_list2 = zero_list2
                        FoM_list2 = zero_list2
                    
                    ##compute mean and std err on mean
                    print(f"\t Compute averages for '{wfitopt_label}' " \
                          f"with COVOPT={covopt_label}")
                    sys.stdout.flush()

                    omm_avg, omm_avg_std = self.compute_average(np.array(omm_list)-np.array(omm_list2))
                    w_avg, w_avg_std     = self.compute_average(np.array(w_list)-np.array(w_list2))
                    wsig_avg, wsig_avg_std = self.compute_average(np.array(wsig_list)-np.array(wsig_list2))
                    if use_wa:
                        wa_avg, wa_avg_std = self.compute_average(np.array(wa_list)-np.array(wa_list2))
                        FoM_avg, FoM_avg_std = self.compute_average(np.array(FoM_list)-np.array(FoM_list2))
                    else:
                        wa_avg, wa_avg_std = 0.0, 0.0
                        FoM_avg, FoM_avg_std = 0.0, 0.0

                    str_nums    = f"{covnum} {wfitnum} "
                    str_results = f"{w_avg:7.4f} {w_avg_std:7.4f} "
                    str_results += f"{wa_avg:7.4f} {wa_avg_std:7.4f} "
                    str_results += f"{omm_avg:7.3f} {omm_avg_std:7.3f}  "
                    str_results += f"{wsig_avg:7.4f} {wsig_avg_std:7.4f}  "
                    
                    if use_wa:
                        str_results += f"{FoM_avg:5.0f} {FoM_avg_std:5.0f}  "

                    str_misc    = f"{len(w_list)}"
                    str_labels  = f"{covopt_label:<10} {wfitopt_label}"
                    nrow +=1
                    f.write(f"ROW: {nrow:3d} {str_nums} {str_results}  " \
                                f"{str_misc} {str_labels}\n")

        f.close()

    # end make_wfitavg_summary

    # === END: ===
