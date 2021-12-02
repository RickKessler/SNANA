# Created Oct 4 2021 by R.Kessler
# Point to directory created by "create_covariance.py"
# and run wfit on hubble_diagram with all covsys_* files
# and with all user-specified WFITOPT options.
 
import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time
import submit_util as util

from submit_params    import *
from submit_prog_base import Program

PREFIX_wfit   = "wfit"

# define names of files produced by create_covariance.py
PREFIX_covsys = "covsys"
HD_FILENAME   = "hubble_diagram.txt"  #
COV_INFO_FILENAME = "INFO.YML"        # contains ISREAL_DATA flag
KEYNAME_ISDATA    = "ISDATA_REAL"     # key in cov info file

# define columns for MERGE.LOG
# column 0 is always for STATE
COLNUM_WFIT_MERGE_DIROPT       = 1
COLNUM_WFIT_MERGE_COVOPT       = 2
COLNUM_WFIT_MERGE_WFITOPT      = 3
COLNUM_WFIT_MERGE_NDOF         = 4   # Ndof
COLNUM_WFIT_MERGE_CPU          = 5

KEYNAME_WFITOPT    = "WFITOPT"
KEYNAME_BLIND_DATA = "BLIND_DATA"
KEYNAME_BLIND_SIM  = "BLIND_SIM"

BLIND_DATA_DEFAULT = True
BLIND_SIM_DEFAULT  = False

JOB_SUFFIX_TAR_LIST  = [ 'YAML', 'DONE', 'LOG'  ]

WFIT_SUMMARY_FILE     = "WFIT_SUMMARY.FITRES"

# - - - - - - - - - - - - - - - - - - -  -
class wFit(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_WFIT
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
            log_assert(False,msgerr)

        return output_dir_name,SUBDIR_SCRIPTS_WFIT

    def submit_prepare_driver(self):

        print("")

        # store input directories, and covsys_*
        self.wfit_prep_input_list()

        # store wfit options under WFITOPT key
        self.wfit_prep_wfitopt_list()

        # prepare blind flag for each inpdir based on data or sim
        self.wfit_prep_blind()

        # prepare index list for inpdir/covsys/wfitopt to simplify
        # 3D loops
        self.wfit_prep_index_lists()

        # copy input file to outdir
        input_file    = self.config_yaml['args'].input_file
        script_dir    = self.config_prep['script_dir']
        shutil.copy(input_file,script_dir)

        # end submit_prepare_driver

    def wfit_prep_input_list(self):

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

            
        # expand inpdir for wildcards and ENVs .xyz
        inpdir_list = []
        for inpdir_orig in inpdir_list_orig:
            inpdir = os.path.expandvars(inpdir_orig)
            if '*' in inpdir_orig:
                tmp_list = sorted(glob.glob(inpdir))
            else:
                tmp_list = [ inpdir ]
            inpdir_list += tmp_list

        # - - - - -
        wildcard      = f"{PREFIX_covsys}*"
        isdata_list   = []
        covsys_list2d = [] # file list per inpdir
        covinfo_list  = [] # list of yaml info per inpdir

        for inpdir in inpdir_list:

            if not os.path.exists(inpdir) :
                msgerr.append(f"Cannot find input directory")
                msgerr.append(f"  {inpdir}")       
                #if inpdir != inpdir_orig : msgerr.append(f"  {inpdir}")
                self.log_assert(False, msgerr)
                
            covsys_list  = sorted(glob.glob1(inpdir,wildcard))
            n_covsys     = len(covsys_list)
            covsys_list2d.append(covsys_list)

            # read ISDATA_REAL flag from INFO.YML file
            isdata, yaml_info = self.read_isdata(inpdir)
            isdata_list.append(isdata)
            covinfo_list.append(yaml_info)

            print(f" Found {inpdir} \n" \
                  f" \t with {n_covsys} {PREFIX_covsys} files and " \
                  f"ISDATA_REAL={isdata} ")


        #sys.exit(f"\n xxx inpdir_list = \n{inpdir_list} \n")
        #print(f" xxx covsys_list = {covsys_list} ")
        # - - - - - -
        #self.config_prep['inpdir_list_orig']  = inpdir_list_orig
        self.config_prep['inpdir_list']       = inpdir_list
        self.config_prep['n_inpdir']          = len(inpdir_list)
        self.config_prep['covsys_list2d']     = covsys_list2d
        self.config_prep['isdata_list']       = isdata_list
        self.config_prep['covinfo_list']      = covinfo_list

        self.wfit_error_check_input_list()

        #print(f" isdata_list = {isdata_list}")

        # end wfit_prep_input_list

    def wfit_error_check_input_list(self):
        # loop over each inpdir and abort on problems such as
        # non-existing inpdir, n_covsys=0, etc ...
        # Print all ERRORS before aborting.

        #inpdir_list_orig = self.config_prep['inpdir_list_orig']
        inpdir_list      = self.config_prep['inpdir_list'] 
        covsys_list2d    = self.config_prep['covsys_list2d']
        nerr = 0
        msgerr = []

        for inpdir, covsys_list in \
            zip(inpdir_list, covsys_list2d):

            hd_file    = f"{inpdir}/{HD_FILENAME}"
            n_covsys   = len(covsys_list)
            if n_covsys == 0 :            
                nerr += 1
                msgerr.append(f"ERROR: cannot find {PREFIX_covsys}* files in")
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

    def read_isdata(self,inpdir):

        yaml_file = f"{inpdir}/{COV_INFO_FILENAME}"
        yaml_info = util.extract_yaml(yaml_file, None, None )
        isdata    = False  # arbitrary default in dase key is missing

        key = KEYNAME_ISDATA
        if key in yaml_info:
            isdata = (yaml_info[KEYNAME_ISDATA] > 0)

        return isdata, yaml_info

        # end read_isdata 
        
    def wfit_prep_wfitopt_list(self):

        msgerr = []
        input_file      = self.config_yaml['args'].input_file 
        CONFIG          = self.config_yaml['CONFIG']
        key             = KEYNAME_WFITOPT

        if key not in CONFIG:
            msgerr.append(f"Missing required {key} key in CONFIG block")
            msgerr.append(f"Check {input_file}")
            self.log_assert(False, msgerr)

        wfitopt_rows = CONFIG[key]
        wfitopt_dict = util.prep_jobopt_list(wfitopt_rows, key, None)

        n_wfitopt          = wfitopt_dict['n_jobopt']
        wfitopt_arg_list   = wfitopt_dict['jobopt_arg_list']
        wfitopt_num_list   = wfitopt_dict['jobopt_num_list']
        wfitopt_label_list = wfitopt_dict['jobopt_label_list']

        logging.info(f"\n Store {n_wfitopt} wfit options from " \
                     f"{key} keys" )

        self.config_prep['n_wfitopt']          = n_wfitopt
        self.config_prep['wfitopt_arg_list']   = wfitopt_arg_list
        self.config_prep['wfitopt_num_list']   = wfitopt_num_list
        self.config_prep['wfitopt_label_list'] = wfitopt_label_list

        # check for global wfitopt
        key   = f"{KEYNAME_WFITOPT}_GLOBAL"
        wfitopt_global = ""
        if key in CONFIG :
            wfitopt_global = CONFIG[key]
        
        self.config_prep['wfitopt_global'] = wfitopt_global

        # check for wa in fit
        use_wa = False
        tmp_list = wfitopt_arg_list + [ wfitopt_global ]
        for tmp in tmp_list:
            if '-wa' in tmp : use_wa = True

        self.config_prep['use_wa'] = use_wa

        # end wfit_prep_wfitopt_list

    def wfit_prep_blind(self):

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

        # end wfit_prep_blind

    def wfit_prep_index_lists(self):

        # prepare internal index lists for efficient looping
        CONFIG           = self.config_yaml['CONFIG']
        inpdir_list      = self.config_prep['inpdir_list']  
        covsys_list2d    = self.config_prep['covsys_list2d']
        wfitopt_list     = self.config_prep['wfitopt_arg_list']

        n_inpdir         = self.config_prep['n_inpdir']
        n_wfitopt        = self.config_prep['n_wfitopt']        

        # count total number of jobs, and beware that number
        # of covsys can be different in each inpdir.
        n_job_tot = 0
        idir_list3 = [];  ifit_list3 = [];   icov_list3 = []

        for idir in range(0,n_inpdir):
            covsys_list = covsys_list2d[idir]
            n_covsys    = len(covsys_list)
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

        # end wfit_prep_index_lists

    def write_command_file(self, icpu, f):

        input_file      = self.config_yaml['args'].input_file 
        inpdir_list     = self.config_prep['inpdir_list']  
        covsys_list2d   = self.config_prep['covsys_list2d']
        wfitopt_list    = self.config_prep['wfitopt_arg_list']
        n_core          = self.config_prep['n_core']

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

            job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_local) 
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

        inpdir      = self.config_prep['inpdir_list'][idir]
        arg_blind   = self.config_prep['arg_blind_list'][idir]
        arg_string  = self.config_prep['wfitopt_arg_list'][ifit]
        arg_global  = self.config_prep['wfitopt_global']
        tmpcov_file = self.config_prep['covsys_list2d'][idir][icov]

        prefix = self.wfit_num_string(idir,icov,ifit)

        covsys_file   = f"{inpdir}/{tmpcov_file}"
        hd_file       = f"{inpdir}/{HD_FILENAME}"
        log_file      = f"{prefix}.LOG" 
        done_file     = f"{prefix}.DONE"
        all_done_file = f"{output_dir}/{DEFAULT_DONE_FILE}"
        
        # start with user-defined args from WFITOPT[_GLOBAL] key
        arg_list =  [ arg_string ]
        if len(arg_global) > 0: arg_list.append(arg_global)

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
            string = f"DIROPT{idir:03d}"

        elif icov >= 0 and idir < 0 and ifit < 0 :
            string = f"COVOPT{icov:03d}"

        elif ifit >= 0 and idir < 0 and icov < 0 :
            string = f"WFITOPT{ifit:03d}"

        elif idir < 0 and ifit>=0 and icov>=0 :
            string = f"COVOPT{icov:03d}_WFITOPT{ifit:03d}"

        elif idir >= 0 and ifit>=0 and icov>=0 :
            string = f"DIROPT{idir:03d}_COVOPT{icov:03d}_WFITOPT{ifit:03d}"
        else:
            string = "ERROR"

        return string
        # end wfit_num_string

    def wfit_prefix(self,row):
        # parse input row passed from MERGE.LOG and construct
        # prefix for output files
        dirnum  = row[COLNUM_WFIT_MERGE_DIROPT]
        covnum  = row[COLNUM_WFIT_MERGE_COVOPT]
        wfitnum = row[COLNUM_WFIT_MERGE_WFITOPT]
        prefix = f"{dirnum}_{covnum}_{wfitnum}"
        return prefix

        # end wfit_prefix

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

            COVOPTS = covinfo['COVOPTS']
            
            n_covopt = len(COVOPTS)
            f.write(f"  COVOPTS({diropt_num}): \n")
            for icov in range(0,n_covopt):
                covopt_num  = self.wfit_num_string(-1,icov,-1)
                covopt      = COVOPTS[icov]
                f.write(f"    {covopt_num}: {covopt} \n")

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
            prefix      = self.wfit_prefix(row) 
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

        # xxx mark delete return [], row_list_merge_new, n_state_change

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
        f = open(SUMMARYF_FILE,"w") 

        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        MERGE_INFO_CONTENTS,comment_lines = \
                util.read_merge_file(MERGE_LOG_PATHFILE)

        dirnum_last = "xyz"

        nrow = 0 
        for row in MERGE_INFO_CONTENTS[TABLE_MERGE]:
            nrow += 1
            dirnum     = row[COLNUM_WFIT_MERGE_DIROPT][-3:] # e.g., 000
            covnum     = row[COLNUM_WFIT_MERGE_COVOPT][-3:] # e.g., 001
            wfitnum    = row[COLNUM_WFIT_MERGE_WFITOPT][-3:] # idem
            prefix     = self.wfit_prefix(row)
            YAML_FILE  = f"{script_dir}/{prefix}.YAML"

            wfit_yaml        = util.extract_yaml(YAML_FILE, None, None )
            wfit_values_dict = util.get_wfit_values(wfit_yaml)

            w       = wfit_values_dict['w']  
            w_sig   = wfit_values_dict['w_sig']
            omm     = wfit_values_dict['omm']  
            omm_sig = wfit_values_dict['omm_sig']
            chi2    = wfit_values_dict['chi2'] 
            sigint  = wfit_values_dict['sigint']
            blind   = wfit_values_dict['blind']
         
            # extract user labels for cov and wfit
            str_diropt  = 'DIROPT' + dirnum
            str_covopt  = 'COVOPT' + covnum
            dir_name    = INPDIR_LIST[str_diropt]  # create_cov dir
            covopt_dict = INPDIR_LIST[f'COVOPTS({str_diropt})']
            #print(f" xxx covopt_dict = {covopt_dict} ")
            covopt_label  = covopt_dict[str_covopt]
            wfitopt_label = WFITOPT_LIST[int(covnum)][1]

            if use_wa:
                wa      = wfit_values_dict['wa']    
                wa_sig  = wfit_values_dict['wa_sig']
                FoM     = wfit_values_dict['FoM']
                Rho     = wfit_values_dict['Rho']
                
            if nrow == 1:
                self.write_wfit_summary_header(f,wfit_values_dict)

            if dirnum != dirnum_last:
                f.write(f"#\n# {str_diropt}={dir_name}\n")
            str_nums    = f"{dirnum} {covnum} {wfitnum}"
            if use_wa : 
                str_results  = f"{w:.4f} {w_sig:.4f} "
                str_results += f"{wa:7.4f} {wa_sig:7.4f} "
                str_results += f"{FoM:5.1f} {Rho:6.3f} "
                str_results += f"{omm:.3f} {omm_sig:.3f}  "
            else:
                str_results  = f"{w:.5f} {w_sig:.5f}  "
                str_results += f"{omm:.5f} {omm_sig:.5f}  "

            str_misc    = f"{chi2:4.1f} {blind} "
            str_labels  = f"{covopt_label:<10} {wfitopt_label}"
            f.write(f"ROW: {nrow:3d} {str_nums} {str_results}  " \
                    f"{str_misc} {str_labels}\n")

            dirnum_last = dirnum

        f.close()
        # end make_wfit_summary

    def write_wfit_summary_header(self,f,wfit_values_dict):
        # write header info and VARNAMES for wfit-summary file

        submit_info_yaml = self.config_prep['submit_info_yaml']
        use_wa           = submit_info_yaml['USE_wa']

        varnames_w  = "w w_sig"
        if use_wa: varnames_w = "w0 w0_sig wa wa_sig FoM Rho"
        VARNAMES_STRING = \
            f"ROW  iDIR iCOV iWFIT {varnames_w} "  \
            f"omm omm_sig  chi2 blind  COVOPT WFITOPT"

        w_ran   = int(wfit_values_dict['w_ran']) 
        wa_ran  = int(wfit_values_dict['wa_ran'])
        omm_ran = int(wfit_values_dict['omm_ran'])

        if w_ran > 0 : 
            f.write(f"#  blind=1 ==> w0,wa,omm += " \
                    f"sin({w_ran},{wa_ran},{omm_ran}) \n")
            f.write(f"\n")
                             
        f.write(f"VARNAMES: {VARNAMES_STRING} \n")
        # write_wfit_summary_header

    def get_misc_merge_info(self):
        # return misc info lines to write into MERGE.LOG file  
        submit_info_yaml = self.config_prep['submit_info_yaml']
        info_lines = []
        return info_lines
        # end get_misc_merge_info 

    def get_merge_COLNUM_CPU(self):
        return -9  # there is no CPU column

    # === END: ===
