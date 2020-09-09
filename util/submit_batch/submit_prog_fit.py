# Created July 2020 by R.Kessler & S. Hinton
#
# - - - - - - - - - -

import os, sys, shutil, yaml, glob
import logging, coloredlogs
import datetime, time, subprocess
import f90nml
#import getpass
#import ntpath
import submit_util as util
import submit_translate as tr
from   submit_params import *
from   submit_prog_base import Program

# define output table-format info
ITABLE_TEXT  = 0
ITABLE_HBOOK = 1
ITABLE_ROOT  = 2
TABLE_NAME_LIST    = [SUFFIX_FITRES, 'SNANA', 'LCPLOT'] # for TEXT format only
TABLE_SUFFIX_LIST  = [ 'TEXT', 'HBOOK', 'ROOT' ]
TABLE_SNLCINP_LIST = [ 'TEXTFILE_PREFIX', 'HFILE_OUT', 'ROOTFILE_OUT' ]
NTABLE_FORMAT = len(TABLE_SUFFIX_LIST)


# list of supplemental file inputs to copy IF no path is specified.
COPY_SNLCINP_FILES = \
    [ "KCOR_FILE", "FLUXERRMODEL_FILE", "HEADER_OVERRIDE_FILE", \
      "MAGCOR_FILE", "SIM_MAGCOR_FILE", "SNCID_LIST_FILE",  \
      "NONLINEARITY_FILE", "FUDGE_HOSTNOISE_FILE", "USERTAGS_FILE" ]

# abort if any of these &SNLCINP inputs are specified 
ABORT_ON_SNLCINP_INPUTS =  \
  [ "OUT_EPOCH_IGNORE_FILE",  "SNMJD_LIST_FILE", "SNMJD_OUT_FILE", \
    "SIMLIB_OUT", "SIMLIB_OUTFILE", "MARZFILE_OUT" ]

# define columns in merge file
COLNUM_FIT_MERGE_STATE           = 0  # STATE required in col=0
COLNUM_FIT_MERGE_VERSION         = 1
COLNUM_FIT_MERGE_FITOPT          = 2
COLNUM_FIT_MERGE_NEVT_ALL        = 3  # NEVT0
COLNUM_FIT_MERGE_NEVT_SNANA_CUTS = 4  # NEVT1
COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS = 5  # NEVT2
COLNUM_FIT_MERGE_CPU             = 6

# hard-code option to validate each version; later pass as CONFIG arg
# 1->weak check on README file, 2->strong check, run snana.exe job
OPT_VALIDATE_VERSION = 1

# all merged files names have MERGE prefix before moving them
# to VERSION subdir and removing MERGE prefix.
PREFIX_MERGE = "MERGE"  

# CONFIG key options to append output FITRES file
KEY_APPEND_TABLE_VARLIST  = "APPEND_TABLE_VARLIST"
KEY_APPEND_TABLE_TEXTFILE = "APPEND_TABLE_TEXTFILE"

# prefix for short snana jobs to validate each version
PREFIX_TEMP_SNANA = "TEMP_SNANA"

# define script to dump number of events in table (for hbook & root)
SCRIPT_SNTABLE_DUMP    = "sntable_dump.pl"  # ?? convert to python ??

# define program to merge text-fitres files
PROGRAM_COMBINE_FITRES = "combine_fitres.exe"

# flags for debug utility to force table-merge failure
FLAG_FORCE_MERGE_TABLE_MISSING = 1
FLAG_FORCE_MERGE_TABLE_CORRUPT = 2

# ====================================================
#    BEGIN FUNCTIONS
# ====================================================

class LightCurveFit(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_FIT
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
            util.log_assert(False,msgerr) # just abort, no done stamp

        return output_dir_name,SUBDIR_SCRIPTS_FIT
        # end set_output_dir_name

    def translate_input_file(self, legacy_input_file, refac_input_file):
        logging.info(f"\n TRANSLATE LEGACY split_and_fit INPUT FILE: " \
                     f"{legacy_input_file}")
        tr.FIT_legacy_to_refac(legacy_input_file, refac_input_file)

    def submit_prepare_driver(self):

        CONFIG       = self.config_yaml['CONFIG']
        input_file   = self.config_yaml['args'].input_file 

        # read/store &SNLCINP namelist
        nml = f90nml.read(input_file)
        self.config_prep['snlcinp'] = nml['snlcinp']

        # abort immediately if particular snlcinp options are/aren't set
        self.fit_prep_check_SNLCINP()

        # assemble list all possible paths where data/sim could be
        self.fit_prep_path_list()

        # get list of versions to fit, along with path to each version
        self.fit_prep_VERSION()

        # get list of FITOPT-options
        self.fit_prep_FITOPT()

        # compute and store indexing for split jobs
        self.fit_prep_index_lists()

        # check output table options
        self.fit_prep_table_options()

        # copy supplemental input files that don't have a path
        self.fit_prep_copy_files()

        # end submit_prepare_driver

    # - - - - - - - - - 

    def fit_prep_copy_files(self):

        # if supplemental input file (passed from &SNLCINP, &FIT, ... ) 
        # has no path (i., not '/'),
        #  + make sure that it exists (abort if not)
        #  + copy to SPLIT_JOBS_LCFIT where split-jobs run
        #
        # Full paths should be included for all input files defined
        # inside &SNLCINP, in which case this function does nothing.
        # For testing, however, it may be convenient to work with
        # local files; in this case, such files are copied to 
        # SPLIT_JOBS_LCFIT where the split-jobs are run.

        script_dir  = self.config_prep['script_dir']
        snlcinp     = self.config_prep['snlcinp']
        input_file  = self.config_yaml['args'].input_file 

        # always copy primary input file
        shutil.copy(input_file,script_dir)

        msgerr=[]
        copy_list_string = ""  # arg of cp
        for key_infile in COPY_SNLCINP_FILES :
            if key_infile in snlcinp :
                infile     = snlcinp[key_infile]
                no_path    = "/" not in infile
                if no_path :
                    if not os.path.isfile(infile):
                        msgerr.append(f" Missing &SNLCINP input file {infile}")
                        self.log_assert(False,msgerr)
                    copy_list_string += (f"{infile} ")

        if len(copy_list_string) > 0 :
            msg = (f" Copy these input files to {SUBDIR_SCRIPTS_FIT}:\n" \
                   f"   {copy_list_string} \n")
            logging.info(msg)
            os.system(f"cp {copy_list_string} {script_dir}/")

        # fit_prep_copy_files

    def fit_prep_check_SNLCINP(self):

        # Abort on particulare SNLCINP inputs that are not allowed
        # to run in batch mode.

        snlcinp    = self.config_prep['snlcinp']
        key_list_invalid = ""
        for key in ABORT_ON_SNLCINP_INPUTS :
            if key in snlcinp :
                key_list_invalid += (f"{key} ")
        
        if len(key_list_invalid) > 0 :
            msgerr = []
            msgerr.append(f"Found &SNLCINP keys not allowed in batch mode:")
            msgerr.append(f"  {key_list_invalid} ")
            self.log_assert(False,msgerr)

        # end fit_prep_check_SNLCINP

    def fit_prep_path_list(self):

        # find all possible paths where data or sim could be
        # Store array path_check_list
        snlcinp             = self.config_prep['snlcinp']
        path_sim_default    = (f"{SNDATA_ROOT}/SIM") 
        path_data_default   = (f"{SNDATA_ROOT}/lcmerge") 
        path_sim_list_file  = (f"{path_sim_default}/PATH_SNDATA_SIM.LIST")
        path_check_list = []
        
        path_check_list.append(path_data_default)
        path_check_list.append(path_sim_default)

        key = 'private_data_path'
        if key in snlcinp :
            path = os.path.expandvars(snlcinp[key])
            path_check_list.append(path) 

        with open(path_sim_list_file,"r") as f :
            for line in f:
                path = os.path.expandvars(line.rstrip("\n"))
                path_check_list.append(path)
                
        #print(f" xxx path_check_list = {path_check_list} ")

        self.config_prep['path_check_list'] = path_check_list
        # end fit_prep_path_list

    def fit_prep_VERSION(self):        

        # read/store each VERSION and path where it is locaated; 
        # if wildcard, go fishing in data areas and sim areas 
        # specified in $SNDATA_ROOT/SIM/PATH_SNDATA_SIM.LIST
        # Finally, create output_dir/[VERSION] for each version ...
        # this is where merged table files end up.

        input_file      = self.config_yaml['args'].input_file 
        path_check_list = self.config_prep['path_check_list']
        output_dir      = self.config_prep['output_dir']
        CONFIG          = self.config_yaml['CONFIG']
        KEYLIST         = [ 'VERSION' ]
        version_list_tmp = util.get_YAML_key_values(CONFIG,KEYLIST)

        # for each version, check all paths ... expand
        # version wildcards and make sure that each version
        # appears once and only once.

        version_list      = []
        path_version_list = []
        msg_list          = []
        version_notFound_list = [] # for err msg only
        msgerr  = []
        nerr_validate = 0
        opt_validate  = OPT_VALIDATE_VERSION

        logging.info("\n Search for VERSION wildcards and paths: ")
        for v_tmp in version_list_tmp :  # may include wild cards
            found = False
            for path in path_check_list :
                v_list  = sorted(glob.glob(f"{path}/{v_tmp}"))
                for v in v_list:
                    found   = True
                    #j_slash = v.rindex('/');    version = v[j_slash+1:]
                    version = os.path.basename(v)
                    # avoid tar files and gz files
                    if '.tar' in v : continue
                    if '.gz'  in v : continue

                    validate,msg_status = \
                        self.fit_validate_VERSION(opt_validate,path,version)
                    msg = (f"   Found VERSION {version}    {msg_status}")
                    logging.info(f"{msg}")
                    if validate is False :
                        nerr_validate += 1  # store all errors before abort
                    version_list.append(version)
                    path_version_list.append(path)
                    msg_list.append(msg)  # used below in case of log_assert

            if not found :
                logging.info(f"   ERROR: could not find {v_tmp} ")
                version_notFound_list.append(v_tmp)

        if len(version_notFound_list) > 0 :
            msgerr.append("Could not find the following VERSION(s)")
            for v in version_notFound_list :
                msgerr.append(f"   {v} ")
            self.log_assert(False,msgerr)

        # check for duplicates 
        n_dupl_list,dupl_list = util.find_duplicates(version_list)
        n_list = len(n_dupl_list)
        if n_list > 0 :
            msgerr.append("Invalid multiple VERSIONs")
            for i in range(0,n_list):
                dupl   = dupl_list[i]
                n_dupl = n_dupl_list[i]
                msgerr.append(f"\t Found {n_dupl} {dupl} directories")
            msgerr.append(f"Only one VERSION directory allowed.")
            msgerr.append(f"Check these directories:")
            for path in set(path_version_list):
                msgerr.append(f"\t {path}")
            self.log_assert(False,msgerr)

        # abort on validation errors, otherwise remove TEMP_SNANA files
        # Note that errors above are stored so that here ALL errors
        # are printed in abort message.
        if nerr_validate > 0 :
            msgerr.append(f"Validation error on {nerr_validate} versions; ")
            msgerr += msg_list
            self.log_assert(False,msgerr)
        else :
            cmd = (f"cd {output_dir}; rm {PREFIX_TEMP_SNANA}* 2>/dev/null")
            os.system(cmd)
            #print(f" xxxx cmd to remove TEMP_SNANA, \n {cmd} ")

        # finally, create subdir for each version
        for v in version_list :
            v_dir = (f"{output_dir}/{v}")
            os.mkdir(v_dir)

        # load lists for version and path_version
        self.config_prep['n_version']         = len(version_list)
        self.config_prep['version_list']      = version_list
        self.config_prep['path_version_list'] = path_version_list

        logging.info("")

        # end fit_prep_VERSION

    def fit_validate_VERSION(self, opt, path, version):
        # Driver to select version-validaton based on input opt.
        validate=True ;   msg= ''
        if opt == 1 :
            validate,msg = self.fit_validate1_VERSION(path,version)
        if opt == 2 :
            validate,msg = self.fit_validate2_VERSION(path,version)

        return validate,msg
        # end fit_validate_VERSION

    def fit_validate1_VERSION(self, path, version):

        # validate based on README file in path.
        # Validation fails if README does not exist,
        # or if FAIL appears on first line of README.

        validate      = True  # default is validation success
        string_status = ""    # no string needed for success
        readme_file = (f"{path}/{version}/{version}.README")
        
        # if no README file, it's a no-brainer failure
        if os.path.isfile(readme_file) is False :
            validate = False
            string_status = "Validate ERROR:  missing README"

        # check if first word in readme file is FAIL 
        # left by sim-submission.
        with open(readme_file,"r") as r :
            line = r.readline().rstrip("\n")
            if 'FAIL' in line :
                validate = False
                string_status = "Validate ERROR : 'FAIL' found in README"

        return validate, string_status
        # end fit_validate1_VERSION

    def fit_validate2_VERSION(self, path, version):

        # Stron validation method:
        # run short snana.exe job on input  version and create
        # YAML file; if NEVT_TOT=0, or YAML file is not produced,
        # return False. If NEVT_TOT > 0, return True
        # Also return comment string; SUCCESS for True,
        # and error msg for False.

        output_dir      = self.config_prep['output_dir']
        cddir           = (f"cd {output_dir}")
        key_nevt        = 'NEVT_TOT'
        nevt_proc       = 10  # process this many to validate
        validate        = True 

        textfile_prefix = (f"{PREFIX_TEMP_SNANA}_{version}")
        log_file        = (f"{textfile_prefix}.LOG")
        yaml_file       = (f"{textfile_prefix}.YAML")
        cmd_snana = "snana.exe NOFILE  "
        cmd_snana += (f"VERSION_PHOTOMETRY {version}  ")
        cmd_snana += (f"MXEVT_PROCESS {nevt_proc}  ")
        cmd_snana += (f"JOBSPLIT 1 1  ")
        cmd_snana += (f"TEXTFILE_PREFIX {textfile_prefix}  " \
                          f"SNTABLE_LIST '' " )

        os.system(f"{cddir} ; {cmd_snana} > {log_file} ")

        # read and parse yaml file
        YAML_FILE = (f"{output_dir}/{yaml_file}")
        string_status = "Validate SUCCESS"  # default status

        if os.path.isfile(YAML_FILE) :
            snana_yaml = util.extract_yaml(YAML_FILE)
            nevt       = snana_yaml[key_nevt]
            if nevt == 0 :
                validate      = False
                string_status = (f"Validate ERROR: NEVT=0")
        else:
            validate      = False
            string_status = (f"Validate ERROR: cannot analyze with snana")

        return validate, string_status

        # end fit_validate2_VERSION

    def fit_prep_FITOPT(self) :

        # read/store list of FITOPT options, along with 'fitnums'
        # FITOPT000, FITOPT001, etc ... These fitnums are used 
        # to create file names.
        # If there is a lable in /LABEL/, strip it out and store it
        # in SUBMIT_INFO, FITOPT.README along with fitopt and fitnum,

        output_dir        = self.config_prep['output_dir']
        CONFIG            = self.config_yaml['CONFIG']
        KEYLIST           = [ 'FITOPT' ]    # key under CONFIG
        fitopt_arg_list   = [ '' ]  # FITOPT000 is always blank
        fitopt_arg_list  += (util.get_YAML_key_values(CONFIG,KEYLIST))
        n_fitopt          = len(fitopt_arg_list)
        fitopt_num_list   = [ '' ] * n_fitopt
        fitopt_label_list = [ 'DEFAULT' ] + ['']*(n_fitopt-1)
        link_FITOPT000_list = []

        # prepare fitnum FITOPT[nnn]
        for i in range(0,n_fitopt):
            fitopt_num_list[i] = (f"FITOPT{i:03d}")
            word_list          = fitopt_arg_list[i].split()
            if len(word_list) > 0 :
                has_label =  word_list[0][0] == '/'
                if has_label :
                    label = word_list[0].strip('/')                    
                    fitopt_label_list[i] = label
                    fitopt_arg_list[i]   = " ".join(word_list[1:])
                else :
                    fitopt_label_list[i] = 'None' 
                    fitopt_arg_list[i]   = " ".join(word_list)

            # update list for symbolic links to FITOPT000 [DEFAULT]
            if self.is_sym_link(fitopt_arg_list[i]) :
                link_FITOPT000_list.append(fitopt_num_list[i])

        # - - - 

        logging.info(f"  Found {n_fitopt-1} FITOPT variations.")
        logging.info(f"  link_FITOPT000_list: {link_FITOPT000_list}")

        self.config_prep['n_fitopt']            = n_fitopt
        self.config_prep['fitopt_num_list']     = fitopt_num_list
        self.config_prep['fitopt_arg_list']     = fitopt_arg_list
        self.config_prep['fitopt_label_list']   = fitopt_label_list
        self.config_prep['link_FITOPT000_list'] = link_FITOPT000_list
        self.config_prep['n_fitopt_link']       = len(link_FITOPT000_list)

        self.write_legacy_FITOPT_README()

        # end fit_prep_FITOPT


    def fit_prep_index_lists(self):

        # Construct 1D sparse lists of                   
        #   version   index iver    0 to n_job_tot 
        #   fitopt    index iopt    0 to n_job_tot 
        #   split-job index isplit  0 to n_job_tot 
        #                                                  
        # These 1D arrays can be used in 1D loops in range(0,n_job_tot)
        # instead of 3D for blocks.
        #
        # Notation warning: 
        #   n_job_tot is number of real jobs that are NOT sym links.
        #   n_fitopt  is total number of FITOPS that includes sym links
        #   n_fitopt_link is number of FITOPTs that are sym links
        #
        n_version        = self.config_prep['n_version']
        n_fitopt_tot     = self.config_prep['n_fitopt']  # all FITOPT
        n_fitopt_link    = self.config_prep['n_fitopt_link'] # links to FITOPT000
        version_list     = self.config_prep['version_list']
        fitopt_arg_list  = self.config_prep['fitopt_arg_list']
        n_core           = self.config_prep['n_core']
        do_dump = True

        # first figure out how many split jobs
        n_fitopt_tmp = n_fitopt_tot - n_fitopt_link # number of FITOPTs to process
        n_job_tmp    = n_version * n_fitopt_tmp  # N_job if no splitting
        n_job_split  = int(n_core/n_job_tmp)
        if n_job_split == 0 : n_job_split = 1

        # note that n_job_tot does NOT include symbolic links
        n_job_tot = n_job_split * n_job_tmp

        logging.info(f"  Determine number of split jobs that are not sym links:")
        logging.info(f"    n_version x n_fitopt = {n_version} x {n_fitopt_tmp} "\
                     f"= {n_job_tmp}   # excludes sym links")
        logging.info(f"    n_job[tot,split] = {n_job_tot},{n_job_split}" \
                     f" for {n_core} cores   # excludes sym links" )
        logging.info("")

        # create sparse index lists that include sym links
        iver_list=[] ;  iopt_list=[];   isplit_list=[]
        size_sparse_list = 0
        for iver in range(0,n_version):
            for iopt in range(0,n_fitopt_tot):
                for isplit in range(0,n_job_split):
                    iver_list.append(iver)
                    iopt_list.append(iopt)
                    isplit_list.append(isplit)
                    size_sparse_list += 1  # n_job(proc+links)

        self.config_prep['size_sparse_list'] = size_sparse_list
        self.config_prep['n_job_tot']     = n_job_tot  # does NOT include symLinks
        self.config_prep['n_done_tot']    = size_sparse_list # proc+links
        self.config_prep['n_job_split']   = n_job_split
        self.config_prep['iver_list']     = iver_list
        self.config_prep['iopt_list']     = iopt_list
        self.config_prep['isplit_list']   = isplit_list

        # store number of jobs that are simply symbolic links
        n_job_link = n_version * n_fitopt_link * n_job_split
        self.config_prep['n_job_link']    = n_job_link

        # end fit_prep_index_lists

    def fit_prep_table_options(self):
        
        # check table options in &SNLCINP namelist, and set
        # logical flag for each table format: TEXT, HBOOK, ROOT
        # In &SNLCINP, existance of HFILE_OUT is a logical flag
        # to create an output HBOOK file; ROOTFILE_OUT is a flag
        # to create a root file. TEXT format is always output
        # for back-compatibility.

        msgerr = []
        snlcinp     = self.config_prep['snlcinp']        
        script_dir  = self.config_prep['script_dir']

        self.config_prep['use_table_format'] = [ False ] * NTABLE_FORMAT
        use_table_format = self.config_prep['use_table_format']
        use_table_format[ITABLE_TEXT] = True  # always force this format

        for fmt in range(0,NTABLE_FORMAT) :
            key = TABLE_SNLCINP_LIST[fmt]
            suffix = TABLE_SUFFIX_LIST[fmt]
            if key in snlcinp :
                use_table_format[fmt] = True

            logging.info(f"  Write {suffix:6s} output table : " \
                         f"{use_table_format[fmt]}")

        # get list of tables in SNTABLE_LIST
        sntable_string = snlcinp['sntable_list']
        sntable_list   = []
        #print(f" xxx sntable_string = {sntable_string}")

        for string in sntable_string.split() :
            sntable_list.append(string)

        logging.info(f"  SNTABLE_LIST = {sntable_list} ")

        # if appending variables to FITRES file, make sure
        # that either HBOOK or ROOT is specified
        CONFIG    = self.config_yaml['CONFIG']
        key       = KEY_APPEND_TABLE_VARLIST
        if key in CONFIG :
            require = False
            if use_table_format[ITABLE_ROOT]:
                require = True
            if use_table_format[ITABLE_HBOOK]:
                require = True
            if not require:
                msgerr.append(f" {key} found in CONFIG input")
                msgerr.append(f" but could not find HBOOK or ROOT. ")
                self.log_assert(False,msgerr)
            
        # if APPEND_TABLE_TEXTFILE is defined, make sure that it exists
        # and that it has a full path.
        key = KEY_APPEND_TABLE_TEXTFILE
        if key in CONFIG :
            text_file = os.path.expandvars(CONFIG[key])
            msg = (f"  Every output FITRES file will be appended with: \n"\
                   f"    {text_file} ")
            logging.info(msg)

            if not os.path.isfile(text_file):
                msgerr.append(f"{text_file} does not exist.")
                msgerr.append(f"Check {key} argument under CONFIG block.")
                self.log_assert(False,msgerr)

            if '/' not in text_file :
                shutil.copy(text_file,script_dir)

        logging.info("")

        # end fit_prep_table_options


    def write_command_file(self,icpu,COMMAND_FILE):
        
        # loop over version, fitopt
        size_sparse_list = self.config_prep['size_sparse_list'] 
        n_job_tot        = self.config_prep['n_job_tot']  # does NOT include symlinks
        iver_list        = self.config_prep['iver_list']
        iopt_list        = self.config_prep['iopt_list']
        isplit_list      = self.config_prep['isplit_list']

        n_version        = self.config_prep['n_version']
        n_fitopt         = self.config_prep['n_fitopt']
        version_list     = self.config_prep['version_list']
        fitopt_arg_list  = self.config_prep['fitopt_arg_list']
        fitopt_num_list  = self.config_prep['fitopt_num_list']
        n_core           = self.config_prep['n_core']

        # open CMD file for this icpu
        f = open(COMMAND_FILE, 'a') 

        n_job_local = 0 ;   n_job_real=0 

        for job in range(0,size_sparse_list):  # n_job(proc+links)
            iver   = iver_list[job]
            iopt   = iopt_list[job]
            isplit = isplit_list[job]
            index_dict = {
                'iver':iver, 'iopt':iopt, 'isplit':isplit, 'icpu':icpu
            }  
            n_job_local += 1
            if self.is_sym_link(fitopt_arg_list[iopt]) : continue
            n_job_real += 1  # use this to skip links

            #if ( (n_job_local-1) % n_core ) == icpu :
            if ( (n_job_real-1) % n_core ) == icpu :

                job_info_fit   = self.prep_JOB_INFO_fit(index_dict)
                util.write_job_info(f, job_info_fit, icpu)

                job_info_merge = self.prep_JOB_INFO_merge(icpu,n_job_real) 
                util.write_jobmerge_info(f, job_info_merge, icpu)

        # - - - - 
        if n_job_real != n_job_tot :
            msgerr = []
            msgerr.append(f"Expected {n_job_tot} total jobs;")
            msgerr.append(f"but found {n_job_local} jobs.")
            self.log_assert(False,msgerr)

        # end write_command_file

    def is_sym_link(self,fitopt_arg):
        # for input fitopt argument, return True if it means
        # symbolic link to FITOPT000.
        if fitopt_arg == 'FITOPT000'    : return True
        if fitopt_arg == 'DEFAULT'      : return True
        return False

    def prep_JOB_INFO_fit(self,index_dict):
        # Return JOB_INFO dictionary with 
        #   cd job_dir
        #   program arg_list  > log_file
        #   touch TMP_[xxx].DONE
        #
        # Inputs
        #   index_dict = dictionary of indices for this job
        #

        # strip off indices from input dictionary
        iver   = index_dict['iver']
        iopt   = index_dict['iopt']
        isplit = index_dict['isplit']+1  # fortran like index for file names
        icpu   = index_dict['icpu']

        input_file    = self.config_yaml['args'].input_file 
        program       = self.config_prep['program']
        output_dir    = self.config_prep['output_dir']
        script_dir    = self.config_prep['script_dir']
        version       = self.config_prep['version_list'][iver]
        fitopt_arg    = self.config_prep['fitopt_arg_list'][iopt]
        fitopt_num    = self.config_prep['fitopt_num_list'][iopt]
        use_table_format = self.config_prep['use_table_format']
        #n_job_tot     = self.config_prep['n_job_tot']
        n_job_split   = self.config_prep['n_job_split']
        split_num     = (f"SPLIT{isplit:03d}")
        prefix        = (f"{version}_{fitopt_num}_{split_num}")
        done_file     = (f"{prefix}.DONE")
        log_file      = (f"{prefix}.LOG")
        yaml_file     = (f"{prefix}.YAML")
        arg_list      = []
        JOB_INFO      = {}

        JOB_INFO['job_dir']     = script_dir  # where to run job
        JOB_INFO['program']     = program
        JOB_INFO['input_file']  = input_file
        JOB_INFO['log_file']    = log_file
        JOB_INFO['done_file']   = done_file

        # set command line arguments
        arg_list.append(f"  VERSION_PHOTOMETRY {version}")
        arg_list.append(f"  JOBSPLIT {isplit} {n_job_split}")

        # check fast option to prescale sims by 10 (data never pre-scaled)
        if self.config_yaml['args'].fast :
            arg_list.append(f"  SIM_PRESCALE {FASTFAC}")

        if self.config_yaml['args'].require_docana :
            arg_list.append(f"  REQUIRE_DOCANA 1")

        # tack on outFile for each table format. For TEXT, do NOT
        # include suffix in TEXTFILE_PREFIX argument
        for itab in range(0,NTABLE_FORMAT) :
            if use_table_format[itab] :
                key    = TABLE_SNLCINP_LIST[itab].upper()
                suffix = TABLE_SUFFIX_LIST[itab]
                arg    = (f"{key:<16} {prefix}")
                if suffix != TABLE_SUFFIX_LIST[ITABLE_TEXT] :
                    arg += (f".{suffix}")
                arg_list.append(f"  {arg}")

        # finally, the user-define FITOPT options. 
        arg_list.append(f"{fitopt_arg}")

        JOB_INFO['arg_list']   = arg_list

        # - - - - - - - - 
        # On FITOPT000, check to create sym links from other FITOPTs
        # The sym links are created AFTER FITOPT000 finishes to 
        # ensure time-ordered file creation for merge process.
        link_FITOPT000_list = self.config_prep['link_FITOPT000_list']
        n_link              = len(link_FITOPT000_list)
        if iopt == 0 and n_link > 0 :
            sym_link_dict = { 
                'version':version, 'prefix_file':prefix, 
                'fitopt_num':fitopt_num }
            JOB_INFO['sym_link_list'] = self.get_sym_link_list(sym_link_dict)
        # - - - - - - - - 
            
        return JOB_INFO

        # end prep_JOB_INFO_fit
 
    def get_sym_link_list(self,sym_link_dict):
        # prepare commands for symbolic links from FITOPTnnn to FITOPT000,
        # where nnn > 0
        # This function is called when user input includes example such as
        # FITOPT:
        #  - bla bla
        #  - bla2 bla2
        #  - FITOPT000  
        # Here the sym link commands are created for DONE, LOG, YAML 
        # outputs, and for FITOP003 -> FITOPT000

        link_FITOPT000_list = self.config_prep['link_FITOPT000_list']
        version             = sym_link_dict['version']
        fitopt_num          = sym_link_dict['fitopt_num']
        prefix_file         = sym_link_dict['prefix_file']
        suffix_link_list    = [ 'LOG ', 'DONE', 'YAML' ]
        sym_link_list   = []  # init return arg list

        for fitopt_link in link_FITOPT000_list : # e.g., FITOPT003
            prefix_link  = prefix_file.replace(fitopt_num,fitopt_link)
            for suffix in suffix_link_list :
                orig_file   = (f"{prefix_file}.{suffix}") # with FITOPT000
                link_file   = (f"{prefix_link}.{suffix}") # with FITOPTnnn
                sym_link    = (f"ln -s {orig_file} {link_file}")
                sym_link_list.append(sym_link)

        return sym_link_list
        # end get_sym_link_list

    def write_legacy_FITOPT_README(self):

        # although SUBMIT_INFO contains the FITOPT info in YAML format,
        # here the legacy FITOPT.README file is also created for 
        # back-compatibility with dowstream scripts.

        output_dir        = self.config_prep['output_dir']
        n_fitopt          = self.config_prep['n_fitopt']
        fitopt_arg_list   = self.config_prep['fitopt_arg_list']
        fitopt_num_list   = self.config_prep['fitopt_num_list']
        fitopt_label_list = self.config_prep['fitopt_label_list']

        legacy_readme_file = (f"{output_dir}/FITOPT.README")
        with open(legacy_readme_file,"w") as f:
            f.write(f"Legacy file for back-compatibility; " \
                    f"please switch to using {SUBMIT_INFO_FILE}\n")
            for i in range(0,n_fitopt):
                num   = fitopt_num_list[i]    # FITOPnnn
                label = fitopt_label_list[i]  # optional input
                if label == 'None' :
                    legacy_label = ''
                else:
                    legacy_label = (f"[{label}]")
                arg   = fitopt_arg_list[i]
                f.write(f"FITOPT: {num[6:]} {legacy_label} {arg} \n")

        # end write_legacy_FITOPT_README


    def append_info_file(self,f):

        # Create SUBMIT.INFO file for merge task and also for
        # downstream scripts.

        CONFIG            = self.config_yaml['CONFIG']
        #n_job_tot         = self.config_prep['n_job_tot']
        #n_job_split       = self.config_prep['n_job_split']
        n_job_link        = self.config_prep['n_job_link']
        output_dir        = self.config_prep['output_dir']
        n_fitopt          = self.config_prep['n_fitopt']
        version_list      = self.config_prep['version_list']
        fitopt_arg_list   = self.config_prep['fitopt_arg_list']
        fitopt_num_list   = self.config_prep['fitopt_num_list']
        fitopt_label_list = self.config_prep['fitopt_label_list']
        link_FITOPT000_list = self.config_prep['link_FITOPT000_list']
        use_table_format  = self.config_prep['use_table_format']

        f.write(f"\n# Fit info\n")
        f.write(f"N_JOB_LINK:          {n_job_link}   " \
                f"# Njob with link to FITOPT000\n")
        f.write(f"JOBFILE_WILDCARD:    '*SPLIT*' \n")
        f.write(f"TABLE_FORMATS:       {TABLE_SUFFIX_LIST} \n")
        f.write(f"USE_TABLE_FORMAT:    {use_table_format} \n")

        key_misc_list = [ KEY_APPEND_TABLE_VARLIST, KEY_APPEND_TABLE_TEXTFILE]
        for key in key_misc_list :
            key_yaml = key + ':'
            if key in CONFIG :
                f.write(f"{key_yaml:<22}  {CONFIG[key]} \n")

        f.write("\n")
        f.write("VERSION_LIST: \n")
        for v in version_list :
            f.write(f"  - {v}\n")

        f.write("\n")
        f.write(f"FITOPT_LIST: \n")
        row = [ 0 ] * 3
        for i in range(0,n_fitopt):
            num   = fitopt_num_list[i]
            label = fitopt_label_list[i] 
            arg   = fitopt_arg_list[i]
            row[COLNUM_FITOPT_NUM]   =  num
            row[COLNUM_FITOPT_LABEL] =  label
            row[COLNUM_FITOPT_ARG]   =  arg
            #row   = [ num, label, arg ]
            f.write(f"  - {row} \n")
        f.write("\n")

        f.write(f"LINK_FITOPT000_LIST: {link_FITOPT000_list}\n")
        f.write("\n")

        # end  append_info_file

    def create_merge_table(self,f):

        # Required element of submit process. Before submitting jobs,
        # create initial merge file with all WAIT states.
        # This file is read and updated frequently by merge
        # process invoked by -m argument to submit_batch_jobs.py.
        # A locally defined MERGE_INFO structure is passed to 
        # a generic write_MERGE_INFO function to create MERGE.LOG/
        # Uses YAML format, and for human-readability there is a 
        # one line commented header before each YAML table.

        n_version       = self.config_prep['n_version']        
        version_list    = self.config_prep['version_list']
        output_dir      = self.config_prep['output_dir']
        fitopt_num_list = self.config_prep['fitopt_num_list']
        n_job_split     = self.config_prep['n_job_split']

        # create only MERGE table ... no need for SPLIT table
        header_line_merge = \
            " STATE   VERSION  FITOPT  " \
            "NEVT_ALL  NEVT_SNANACUT NEVT_FITCUT  CPU"

        INFO_MERGE = { 
            'primary_key' : TABLE_MERGE, 'header_line' : header_line_merge,
            'row_list'    : []   }

        STATE = SUBMIT_STATE_WAIT # all start in WAIT state
        for version in version_list :
            for num in fitopt_num_list :
                # define ROW here is fragile in case columns are changed
                ROW_MERGE = [ STATE, version, num, 0, 0, 0, 0.0 ]
                INFO_MERGE['row_list'].append(ROW_MERGE)  
        util.write_merge_file(f, INFO_MERGE, [] ) 

        # end create_merge_file

    def merge_config_prep(self,output_dir):

        # fit-specific settings to config_prep that are needed later.
        submit_info_yaml = self.config_prep['submit_info_yaml'] 
        ## self.config_prep['output_dir']   = output_dir 

        # end merge_config_prep

    def merge_update_state(self, MERGE_INFO_CONTENTS):

        # read MERGE.LOG, check LOG & DONE files.
        # Return update row list MERGE tables.
        # To continuously update statistics during RUN state,
        # n_state_change is always set > 0 if in RUN state,
        # even if there is not STATE change.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        link_FITOPT000_list = submit_info_yaml['LINK_FITOPT000_LIST']
        COLNUM_STATE   = COLNUM_FIT_MERGE_STATE 
        COLNUM_VERS    = COLNUM_FIT_MERGE_VERSION 
        COLNUM_FITOPT  = COLNUM_FIT_MERGE_FITOPT
        COLNUM_NEVT0   = COLNUM_FIT_MERGE_NEVT_ALL 
        COLNUM_NEVT1   = COLNUM_FIT_MERGE_NEVT_SNANA_CUTS
        COLNUM_NEVT2   = COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS
        COLNUM_CPU     = COLNUM_FIT_MERGE_CPU

        key_tot, key_tot_sum, key_list = \
                self.keynames_for_job_stats('NEVT_TOT')
        key_snana, key_snana_sum, key_snana_list = \
                self.keynames_for_job_stats('NEVT_SNANA_CUTS')
        key_lcfit, key_lcfit_sum, key_lcfit_list = \
                self.keynames_for_job_stats('NEVT_LCFIT_CUTS')
        key_cpu, key_cpu_sum, key_cpu_list = \
                self.keynames_for_job_stats('CPU_MINUTES')
        key_list  = [ key_tot, key_snana, key_lcfit, key_cpu ]

        row_list_merge   = MERGE_INFO_CONTENTS[TABLE_MERGE]

        # init outputs of function
        n_state_change     = 0
        row_list_merge_new = []

        #  - - - - -
        irow = 0
        for row in row_list_merge :
            row_list_merge_new.append(row) # default output is same as input

            # strip off row info
            STATE       = row[COLNUM_STATE]
            version     = row[COLNUM_VERS]
            fitopt_num  = row[COLNUM_FITOPT]
            search_wildcard = (f"{version}_{fitopt_num}_SPLIT*")

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
                    # check for failures in snlc_fit jobs.
                    nfail = job_stats['nfail']
                    if nfail > 0 :  NEW_STATE = SUBMIT_STATE_FAIL
                
                    # xxx sum_stats=self.split_sum_stats(False,log_list,yaml_list)
                    row[COLNUM_STATE]  = NEW_STATE
                    row[COLNUM_NEVT0]  = job_stats[key_tot_sum]
                    row[COLNUM_NEVT1]  = job_stats[key_snana_sum]
                    row[COLNUM_NEVT2]  = job_stats[key_lcfit_sum]
                    row[COLNUM_CPU]    = job_stats[key_cpu_sum]

                    if fitopt_num in link_FITOPT000_list :
                        row[COLNUM_CPU] = 0.0 # zero CPU for sym links

                    row_list_merge_new[irow] = row  # update new row
                    n_state_change += 1             # assume nevt changes

            irow += 1

        # first return arg (row_split) is null since there is 
        # no need for a SPLIT table
        return [], row_list_merge_new, n_state_change

        # end merge_update_state

    def split_sum_stats(self, search_failure_flag, log_list, yaml_list):

        # xxxxxxxxxx OBSOLETE MARK DELETE xxxxxxxxxxx

        # Return statistics sums for yaml_list files.
        # If search_failure_flag = Flase, then examine only the yaml_list
        # and do not check for failures.
        # When all DONE files exist, this function is called with
        # search_failure_flag=True so that failures are examined too.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_log_file       = len(log_list)
        split_stats = {
            'nevt_sum_tot'        : 0, 
            'nevt_sum_cut_snana'  : 0,
            'nevt_sum_cut_lcfit'  : 0,
            'cpu_sum'             : 0.0,
            'nfail_sum'           : 0
        }

        # xxxxxxxxxx OBSOLETE MARK DELETE xxxxxxxxxxx

        for isplit in range(0,n_log_file):            
            yaml_file = yaml_list[isplit]            
            nevt_test = -9  # used to search for failures
            if yaml_file :
                YAML_FILE = (f"{script_dir}/{yaml_file}")
                yaml_data = util.extract_yaml(YAML_FILE)
                split_stats['nevt_sum_tot']       += yaml_data['NEVT_TOT']
                split_stats['nevt_sum_cut_snana'] += yaml_data['NEVT_SNANA_CUTS']
                split_stats['nevt_sum_cut_lcfit'] += yaml_data['NEVT_LCFIT_CUTS']
                split_stats['cpu_sum']            += yaml_data['CPU_TIME']

                # fix cpu format
                cpu = (f"{split_stats['cpu_sum']:.1f}")
                split_stats['cpu_sum']  = float(cpu)

                # test value for failure testing below
                nevt_test = yaml_data['ABORT_IF_ZERO'] 

        # xxxxxxxxxx OBSOLETE MARK DELETE xxxxxxxxxxx

            # check flag to check for failure.        
            if search_failure_flag and nevt_test <= 0 :
                log_file   = log_list[isplit]
                found_fail = self.check_for_failure(log_file,nevt_test,isplit+1)
                if found_fail :
                    split_stats['nfail_sum'] += 1

        return split_stats

        # end split_sum_stats
        # xxxxxxxxxx END OBSOLETE MARK DELETE xxxxxxxxxxx

    def merge_job_wrapup(self, irow, MERGE_INFO_CONTENTS):
        # irow is the row to wrapup in MERGE_INFO_CONTENTS
        # One row corresonds to one VERSION and one FITOPT;
        # combine the SPLIT output tables into one merged table file.
        # Merge separately for each table format; HBOOK, ROOT, TEXT ...
        # Also check for symLink option to replace table-merge
        # with sym link to FITOPT000.

        # init name of merged table file for each format
        self.config_prep['merge_table_file_list'] = [''] * NTABLE_FORMAT

        row         = MERGE_INFO_CONTENTS[TABLE_MERGE][irow]
        state       = row[COLNUM_FIT_MERGE_STATE]
        version     = row[COLNUM_FIT_MERGE_VERSION]
        fitopt_num  = row[COLNUM_FIT_MERGE_FITOPT]  # e.g., FITOPT001
        nevt_expect = row[COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS]

        submit_info_yaml = self.config_prep['submit_info_yaml']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        use_table_format = submit_info_yaml['USE_TABLE_FORMAT']
        msgerr           = []

        version_fitopt = (f"{version}_{fitopt_num}")

        # start with idiot check to make sure that state is DONE.
        if state != SUBMIT_STATE_DONE :
            msgerr.append(f"Unexpected STATE = {state}; " \
                          f"expected {SUBMIT_STATE_DONE} state.")
            msgerr.append(f"Check {version_fitopt} files.")
            self.log_assert(False,msgerr)

        logging.info(f" Begin table merge for {version_fitopt} ")

        # check for sym link first
        if self.create_sym_link_tables(version,fitopt_num) :
            return

        # - - - - -
        # make sure we have correct number of files to merge
        for itable in range(0,NTABLE_FORMAT):
            suffix         =  TABLE_SUFFIX_LIST[itable]
            if itable == ITABLE_TEXT :
                suffix = (f"{SUFFIX_FITRES}.{suffix}")  # special case for TEXT
            use_format     =  use_table_format[itable] 
            f_wildcard     =  (f"{script_dir}/{version_fitopt}*{suffix}")
            if use_format :
                util.check_file_count(n_job_split, f_wildcard )

        # - - - - -
        # create dictionary to pass to merge_table_XXX functions below.
        version_fitopt_dict = {
            'version'         : version,
            'fitopt'          : fitopt_num,
            'version_fitopt'  : version_fitopt,
            'nevt_expect'     : nevt_expect 
        }

        # - - - - -
        # merge strategy is to merge split files into temp files 
        # named MERGE_{version_fitopt}.{suffix}; in case of abort,
        # these temp files remain. After NEVT validation,
        # each temp file is moved to {version}/{fitopt}.{suffix}

        if use_table_format[ITABLE_HBOOK] :
            self.merge_table_CERN(ITABLE_HBOOK, version_fitopt_dict)

        if use_table_format[ITABLE_ROOT] :
            self.merge_table_CERN(ITABLE_ROOT, version_fitopt_dict)

        # process TEXT format after ROOT & HBOOK to allow for append feature
        if use_table_format[ITABLE_TEXT] :
            for table_name in TABLE_NAME_LIST :
                self.merge_table_TEXT(table_name, version_fitopt_dict)

        # move MERGE files, and remove 'MERGE_' prefix
        self.move_merge_table_files(version_fitopt_dict)

        logging.info("")

        if irow == 33333:
            sys.exit(f"\n\t xxxxx DEBUG DIE from fit wrapup ... xxxx ")

        # end merge_job_wrapup 

    def create_sym_link_tables(self, version, fitopt_num):

        # if this fitopt_num (e.g., FITOPT003) is on list of sym links,
        # create sym link for each merged table expected to appear 
        # in ../version.
        # Function returns true of sym links are create; false otherwise.

        output_dir       = self.config_prep['output_dir']
        submit_info_yaml = self.config_prep['submit_info_yaml']

        use_table_format    = submit_info_yaml['USE_TABLE_FORMAT']
        link_FITOPT000_list = submit_info_yaml['LINK_FITOPT000_LIST']
        cdv                 = (f"cd {output_dir}/{version}")
        create_links        = False  # init return function arg

        if fitopt_num in link_FITOPT000_list :
            fitopt_ref = 'FITOPT000'
            cmd_link_all = ''

            for itab in range(0,NTABLE_FORMAT):
                if not use_table_format[itab] : continue

                suf   =  TABLE_SUFFIX_LIST[itab]
                if itab == ITABLE_TEXT :
                    suf  = (f"{SUFFIX_FITRES}")  # special case

                suf += '.gz'
                cmd_link = (f"ln -s {fitopt_ref}.{suf} {fitopt_num}.{suf} ; ")
                cmd_link_all += cmd_link
                logging.info(f"   create sym-link to {fitopt_ref}.{suf}")
            cmd = (f"{cdv} ; {cmd_link_all}")
            create_links = True
            os.system(cmd)
        
        return create_links
        # end create_sym_link_tables

    def merge_table_TEXT(self, table_name, version_fitopt_dict):

        # Input table_name is either 'FITRES', 'SNANA', or 'LCPLOT'.
        # Problem with TEXT tables is that each table must be written
        # to separate set of files, so here each table name is processed
        # separately. Usually only the FITRES table is requested, but 
        # sometimes other tables are included.
        # HBOOK & ROOT don't have this issue because all tables
        # reside in one file.

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        version          = version_fitopt_dict['version']
        fitopt           = version_fitopt_dict['fitopt']
        version_fitopt   = version_fitopt_dict['version_fitopt']
        nevt_expect      = version_fitopt_dict['nevt_expect']

        itable           = ITABLE_TEXT
        suffix           = TABLE_SUFFIX_LIST[itable]
        prefix           = PREFIX_MERGE  
        table_wildcard  = (f"{version_fitopt}*{table_name}.TEXT")

        # check debug option to force merge-table failure
        flag_force_fail = \
            self.flag_force_merge_table_fail(itable,version_fitopt)

        if flag_force_fail == FLAG_FORCE_MERGE_TABLE_CORRUPT :
            table_wildcard += f(" {table_wildcard}") # double output

        # if no tables exist, bail out
        table_list = glob.glob1(script_dir,table_wildcard)
        if len(table_list) == 0 :
            return

        # construct linux command to catenate TEXT files,
        # and then a special awk command to remove all VARNAMES
        # lines EXCEPT for the first VARNAMES line.
        # Note that output extension is table name, not TEXT
        cddir           = (f"cd {script_dir}")
        out_table_file  = (f"{prefix}_{version_fitopt}.{table_name}")
        out_table_file2 = (f"{prefix}2_{version_fitopt}.{table_name}")# for awk

        cmd_cat  = (f"cat {table_wildcard} > {out_table_file}")
        cmd_awk  = (f"awk '!/^VARNAMES/ || ++n <= 1' {out_table_file} > " \
                    f"{out_table_file2} ; " \
                    f"mv {out_table_file2} {out_table_file}")

        msg = (f"   merge {n_job_split} {suffix}-{table_name} table files.")
        logging.info(msg)

        if flag_force_fail != FLAG_FORCE_MERGE_TABLE_MISSING :
            os.system(f"{cddir}; {cmd_cat} ; {cmd_awk}")

        OUT_TABLE_FILE = (f"{script_dir}/{out_table_file}")
        self.check_file_exists(OUT_TABLE_FILE,["Problem with table-merge"])

        # - - - - - -
        # for FITRES table only, do unitarity check to make sure that 
        # nevt_expect rows are really there. Also check options to 
        # append variables to table.

        if table_name == SUFFIX_FITRES :
            OUT_TABLE_FILE = (f"{script_dir}/{out_table_file}")
            nevt_find = util.nrow_table_TEXT(OUT_TABLE_FILE,"SN:")
            self.nevt_table_check(nevt_expect, nevt_find, out_table_file)
            self.config_prep['merge_table_file_list'][itable] = out_table_file

            # check options to append FITRES file
            # 1. APPEND_TABLE_VARLIST  -> extract vars from HBOOK or ROOT file
            # 2. APPEND_TABLE_TEXTFILE -> append vars from external file.
            self.append_table_varlist(version_fitopt_dict)  # optional
            self.append_table_textfile(version_fitopt_dict)  # optional

        # end merge_table_TEXT

    def append_table_varlist(self,version_fitopt_dict) :

        # Check option to extract variables from HBOOK/ROOT file,
        # and append TEXT-FITRES file.
        # See CONFIG key APPEND_TABLE_VARLIST
        #
        submit_info_yaml      = self.config_prep['submit_info_yaml']
        merge_table_file_list = self.config_prep['merge_table_file_list']
        script_dir            = submit_info_yaml['SCRIPT_DIR']
        use_table_format      = submit_info_yaml['USE_TABLE_FORMAT']
        version_fitopt        = version_fitopt_dict['version_fitopt']
        nevt_expect           = version_fitopt_dict['nevt_expect']
        
        key  = KEY_APPEND_TABLE_VARLIST
        if key in submit_info_yaml :
            varlist_append = submit_info_yaml[key]
        else:
            return

        msgerr = []

        text_table_file = merge_table_file_list[ITABLE_TEXT]  # file to append
        if use_table_format[ITABLE_ROOT] :
            full_table_file = merge_table_file_list[ITABLE_ROOT]
        elif use_table_format[ITABLE_HBOOK] :
            full_table_file = merge_table_file_list[ITABLE_HBOOK]
        else :
            msgerr.append("Cannot append TEXT table without full table")
            msgerr.append("Must define HBOOK or ROOT to append TEXT table")
            self.log_assert(False,msgerr) 

        #FULL_TABLE_FILE = (f"{script_dir}/{full_table_file}")
        #TEXT_TABLE_FILE = (f"{script_dir}/{text_table_file}")

        append_log_file = (f"sntable_append_{PREFIX_MERGE}_{version_fitopt}.log")
        append_out_file = (f"sntable_append_{PREFIX_MERGE}_{version_fitopt}.text")

        logging.info(f"\t Append TEXT table from {full_table_file} ")

        cddir = (f"cd {script_dir}")
        cmd_append = (f"{SCRIPT_SNTABLE_DUMP} {full_table_file} FITRES " \
                      f"-v {varlist_append} " \
                      f"-a {text_table_file} > {append_log_file} 2>/dev/null")
        os.system(f"{cddir} ; {cmd_append} ")

        # replace FITRES file with append file, 
        # and remove sntable* junk files

        if os.path.isfile(f"{script_dir}/{append_out_file}") :
            cmd_mv = (f"mv {append_out_file} {text_table_file}")
            cmd_rm = (f"rm sntable_*")
            cmd    = (f"{cddir}; {cmd_mv} ; {cmd_rm}")
            os.system(cmd)
        else : 
            msgerr.append(f"Append FITRES table failed with command")
            msgerr.append(f"{cmd_append}")
            msgerr.append(f"Check if these variables are all valid:")
            msgerr.append(f"    {varlist_append}")
            self.log_assert(False,msgerr) 

        # finally, make sure that the number of rows still matches nevt_expect
        OUT_TABLE_FILE = (f"{script_dir}/{text_table_file}")
        nevt_find = util.nrow_table_TEXT(OUT_TABLE_FILE,"SN:")
        self.nevt_table_check(nevt_expect, nevt_find, text_table_file)
                  
        # end append_table_varlist

    def append_table_textfile(self,version_fitopt_dict) :

        # Check option to extract and append variables from 
        # external TEXT file. See CONFIG key APPEND_TABLE_TEXFILE
        #
        submit_info_yaml      = self.config_prep['submit_info_yaml']
        merge_table_file_list = self.config_prep['merge_table_file_list']
        script_dir            = submit_info_yaml['SCRIPT_DIR']
        use_table_format      = submit_info_yaml['USE_TABLE_FORMAT']
        version_fitopt        = version_fitopt_dict['version_fitopt']
        nevt_expect           = version_fitopt_dict['nevt_expect']
        
        key  = KEY_APPEND_TABLE_TEXTFILE
        if key in submit_info_yaml :
            external_file = submit_info_yaml[key]
        else:
            return

        orig_file = merge_table_file_list[ITABLE_TEXT]  # file to append
        out_file  = "TEMP_" + orig_file
        log_file  = "TEMP_COMBINE.LOG"

        cddir = (f"cd {script_dir}")
        cmd1  = (f"{PROGRAM_COMBINE_FITRES} {orig_file} {external_file} " \
                 f"-outfile_text {out_file} >& {log_file}" )
        cmd2  = (f"mv {out_file} {orig_file}")
        cmd3  = (f"rm {log_file}")
        cmd   = (f"{cddir} ; {cmd1} ; {cmd2} ; {cmd3}")
        os.system(cmd)

        # finally, make sure that the number of rows still matches nevt_expect
        ORIG_FILE = (f"{script_dir}/{orig_file}")
        nevt_find = util.nrow_table_TEXT(ORIG_FILE,"SN:")
        self.nevt_table_check(nevt_expect, nevt_find, orig_file)

        #sys.exit('\n xxx DEBUG DIE xxxx ')

    # end append_table_textfile

    def merge_table_CERN(self, itable, version_fitopt_dict):

        # call merge program for HBOOK or ROOT based on itable arg.
        # into one.
        version          = version_fitopt_dict['version']
        fitopt           = version_fitopt_dict['fitopt']
        version_fitopt   = version_fitopt_dict['version_fitopt']
        nevt_expect      = version_fitopt_dict['nevt_expect']

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']
        n_job_split      = submit_info_yaml['N_JOB_SPLIT']
        suffix           = TABLE_SUFFIX_LIST[itable]
        prefix           = PREFIX_MERGE
        msgerr           = []

        # check debug option to force merge-table failure
        flag_force_fail = self.flag_force_merge_table_fail(itable,version_fitopt)

        logging.info(f"   merge {n_job_split} {suffix} table files.")
        
        # get name of program
        program_merge = (f"merge_{suffix.lower()}.exe")

        cddir           = (f"cd {script_dir}")
        out_table_file  = (f"{prefix}_{version_fitopt}.{suffix}") 
        log_table_file  = (f"{prefix}_{version_fitopt}_{suffix}.LOG") 
        f_list          = (f"{version_fitopt}*{suffix}")
        if flag_force_fail == FLAG_FORCE_MERGE_TABLE_CORRUPT :
            f_list += (f" XXX_FORCE_CORRUPT.{suffix}")

        cmd_merge       = (f"{program_merge} {f_list} {out_table_file} ")
        cmd             = (f"{cddir} ; {cmd_merge} > {log_table_file}  ")

        if flag_force_fail != FLAG_FORCE_MERGE_TABLE_MISSING :
            os.system(cmd)

        OUT_TABLE_FILE = (f"{script_dir}/{out_table_file}")
        self.check_file_exists(OUT_TABLE_FILE, ["Problem with table-merge"])

        # for HBOOK, remove garbage from log file
        if itable == ITABLE_HBOOK :
            cmd_clean_log = \
                (f"{cddir} ; remove_locf_messages.pl {log_table_file} QUIET") 
            os.system(cmd_clean_log)

        # ?? check log file for success message ??
        
        # extract number of events in final HBOOK/ROOT file
        NTRY_MAX = 2; ntry=0; nevt_find = -9
        while nevt_find < 0 and ntry < NTRY_MAX :
            nevt_find = self.nrow_table_CERN(f"{OUT_TABLE_FILE}")
            ntry += 1

        if nevt_find < 0 :
            msgerr.append(f"Cannot get NEVT from merged {suffix} table " \
                          f"after {ntry} attempts; ")
            msgerr.append(f"   {out_table_file}")
            msgerr.append(f"---> table may be corrupted.")
            msgerr.append(f" ")
            msgerr.append(f"Full table-merge command was")
            msgerr.append(f"  {cmd}")
            msgerr.append(f" ")
            self.log_assert(False,msgerr)

        if ntry > 1 :
            msg=(f"\tWARNING: ntry={ntry} to get NEVT for {out_table_file}")
            logging.warning(msg)

        # abort if nevt_find != nevt_expect
        self.nevt_table_check(nevt_expect, nevt_find, out_table_file)

        # store merge table file for append_table_text()
        merge_table_file_list =  self.config_prep['merge_table_file_list']
        merge_table_file_list[itable] = out_table_file

        # end merge_table_CERN

    def nrow_table_CERN(self,table_file):
        # return number of table rows in CERN file that
        # has either HBOOK or ROOT extension
        # Return -9 on error.
        # Script prints "NEVT:  <nevt>", so parse the 2nd element.

        script   = SCRIPT_SNTABLE_DUMP
        cmd_nevt = (f"{script} {table_file} FITRES NEVT " \
                    f" | grep 'NEVT:' ")
        try: 
            result_line = subprocess.check_output(cmd_nevt, shell=True)
            result_line = (result_line.rstrip()).decode('utf-8')
            nevt_find   = int(result_line.split()[1])
            return nevt_find
        except :
            return -9

        #end nrow_table_CERN

    def move_merge_table_files(self,version_fitopt_dict):
       
        # move files of the form  MERGE_FITOPTnnn.FORMAT
        # to VERSION/FITOPTnnn.FORMAT
        # (i.e., strip of the MERGE_ prefix that is for debugging)

        version          = version_fitopt_dict['version']
        fitopt           = version_fitopt_dict['fitopt']

        submit_info_yaml = self.config_prep['submit_info_yaml']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        table_wildcard   = (f"{PREFIX_MERGE}_*")
        table_list       = glob.glob1(script_dir,table_wildcard)

        cmd = (f"cd {script_dir} ; ")
        for tfile in table_list :
            jdot   = tfile.index('.')
            suffix = tfile[jdot+1:]
            if suffix == 'LOG' :
                cmd += (f"rm {tfile} ; ")
            else:
                move_file = (f"../{version}/{fitopt}.{suffix}")
                cmd      += (f"mv {tfile} {move_file} ; ")

        os.system(cmd)

        # end move_merge_table_files

    def nevt_table_check(self, nevt_expect, nevt_find, out_file):
        # abort if nevt_expect != nevt_find
        if nevt_expect != nevt_find :
            msgerr = []
            msgerr.append(f"Found {nevt_find} events in merged {out_file}")
            msgerr.append(f"but expected {nevt_expect} from summing YAML files.")
            msgerr.append(f"Table merge has a problem.")
            self.log_assert(False,msgerr)
    # end nevt_table_check

    def merge_cleanup_final(self):
        
        # every snlc_fit job succeeded, so here we simply compress output.
        # Before compressing, check that every expected merged table
        # file exists because previous error checking is for snlc_fit,
        # not for merging tables.
        #
        # If any merged table file is missing, 
        #  + print list of missing table files in MERGE.LOG (and to stdout)
        #  + abort and change SUCCESS state to FAIL in global DONE file
        #

        output_dir          = self.config_prep['output_dir']
        submit_info_yaml    = self.config_prep['submit_info_yaml']
        use_table_format    = submit_info_yaml['USE_TABLE_FORMAT']
        link_FITOPT000_list = submit_info_yaml['LINK_FITOPT000_LIST']
        script_dir          = submit_info_yaml['SCRIPT_DIR']
        subdir         = SUBDIR_SCRIPTS_FIT
        tar_file       = (f"{subdir}.tar")

        logging.info(f" FIT cleanup: check if all merged tables exist.")

        # check that every expected merged table file exits;
        # read status from MERGE file to get list of version & fitopts
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        MERGE_INFO_CONTENTS,comment_lines = \
            util.read_merge_file(MERGE_LOG_PATHFILE)

        row_list  = MERGE_INFO_CONTENTS[TABLE_MERGE]
        msgerr = ['ERROR merging table files: ' ]
        nerr   = 0
        for row in row_list :
            version = row[COLNUM_FIT_MERGE_VERSION]
            fitopt  = row[COLNUM_FIT_MERGE_FITOPT]  # e.g., FITOPT002
            if fitopt in link_FITOPT000_list : continue # skip sym links

            for itable in range(0,NTABLE_FORMAT):
                use    = use_table_format[itable]
                suffix = TABLE_SUFFIX_LIST[itable]
                if itable == ITABLE_TEXT:
                    suffix = SUFFIX_FITRES  # special case here
                table_file = (f"{version}/{fitopt}.{suffix}")
                TABLE_FILE = (f"{output_dir}/{table_file}")
                if use and not os.path.isfile(TABLE_FILE):
                    nerr += 1
                    msgerr.append(f"    Missing expected {table_file}")
                
        if nerr > 0 :
            self.log_assert(False,msgerr) # abort and write FAIL to global DONE file.
        else :
            logging.info(f" FIT cleanup: all merged tables exist.")

        # if we get here, table-merging seems to have worked so tar and zip

        logging.info(f" FIT cleanup: tar up files under {subdir}/")
        util.compress_files(+1, script_dir, "*SPLIT*.LOG",  "LOG", "" )
        util.compress_files(+1, script_dir, "*SPLIT*.YAML", "YAML", "" )
        util.compress_files(+1, script_dir, "*SPLIT*.DONE", "DONE", "" )

        for itable in range(0,NTABLE_FORMAT):
            use    = use_table_format[itable]
            suffix = TABLE_SUFFIX_LIST[itable]
            if use :
                wildcard = (f"*SPLIT*.{suffix}")
                util.compress_files(+1, script_dir, wildcard, suffix, "" )
                # ?? at some point, should delete these since merged table is there ??

        logging.info(f" FIT cleanup: gzip merged tables.")
        cmd_gzip = (f"cd {output_dir} ; gzip */FITOPT* 2>/dev/null")
        os.system(cmd_gzip)

        self.merge_cleanup_script_dir() 

        # xxx nothing to write to? logging.info(f" FIT cleanup: Done.")

        # end merge_cleanup_final

    def merge_reset(self,output_dir):
        # remove all merge products to allow re-doing the merge with -m option.
        # Renders output_dir as if --nomerge option had been given with
        # original submit command.
        #  + set all STATE values to WAIT
        #  + set all NEVT to 0
        #  + set CPU to 0
        #  + unpack SPLIT_JOBS_LCFIT
        #  + remove merged table files: [VERSION]/FITOPT* files
        #  + remove all-done file
        
        submit_info_yaml = self.config_prep['submit_info_yaml']
        version_list     = submit_info_yaml['VERSION_LIST']
        script_dir       = submit_info_yaml['SCRIPT_DIR']

        fnam = "merge_reset"

        # read status from MERGE file          
        MERGE_LOG_PATHFILE  = (f"{output_dir}/{MERGE_LOG_FILE}")
        colnum_zero_list    = [ COLNUM_FIT_MERGE_NEVT_ALL, 
                                COLNUM_FIT_MERGE_NEVT_SNANA_CUTS,
                                COLNUM_FIT_MERGE_NEVT_LCFIT_CUTS,
                                COLNUM_FIT_MERGE_CPU]
        logging.info(f"  {fnam}: STATE->WAIT and NEVT->0 in {MERGE_LOG_FILE}")
        util.merge_table_reset(MERGE_LOG_PATHFILE, TABLE_MERGE,  \
                               COLNUM_MERGE_STATE, colnum_zero_list)

        # loop over each version and clean out FITOPT* files 
        logging.info(f"  {fnam}: remove FITOPT* from version subdirs")
        for version in version_list :
            v_dir = (f"{output_dir}/{version}")
            cmd_rm = (f"cd {v_dir} ; rm FITOPT* 2>/dev/null")
            os.system(cmd_rm)


        # if script_dir is tarred & gzipped, unpack it
        logging.info(f"  {fnam}: unapck {SUBDIR_SCRIPTS_FIT}")
        util.compress_subdir(-1, script_dir)

        # untar and unzip file inside SUBDIR_SCRIPTS_FIT
        cmd_unzip = (f"cd {script_dir}; cat BACKUP*.tar.gz | tar xzf - -i ; " \
                     f"rm BACKUP*.gz")
        os.system(cmd_unzip)

        # remove lingering temp files in SPLIT_JOBS_LCFIT.
        # Also remove MERGE.LOG_{Nsec} backups, and DONE file
        logging.info(f"  {fnam}: remove misc junk files from {SUBDIR_SCRIPTS_FIT}")
        cmd_rm1 = (f"cd {script_dir} ; rm {PREFIX_MERGE}* sntable* 2>/dev/null")
        cmd_rm2 = (f"cd {output_dir} ; rm {MERGE_LOG_FILE}_*  *.DONE 2>/dev/null")
        os.system(cmd_rm1)
        os.system(cmd_rm2)
        #print(f" cmd_rm1 = {cmd_rm1}")
        #print(f" cmd_rm2 = {cmd_rm2}")

        # end merge_reset

    def flag_force_merge_table_fail(self, itable, version_fitopt):

        # Check CONFIG (user input) for 
        #   FORCE_MERGE_TABLE_MISSING({table_format})
        #   FORCE_MERGE_TABLE_CORRUPT({table_format})
        #
        # and return appropriate flag. This flag is used bu
        # table-merge function to force failure.
        # Goal is to debug error handling on table-merge failures.

        flag   = 0 # init output
        CONFIG = self.config_yaml['CONFIG']

        #if 'FORCE' not in CONFIG :
        #    return flag

        suffix = TABLE_SUFFIX_LIST[itable]

        key_list = [ (f"FORCE_MERGE_TABLE_MISSING({suffix})"), 
                     (f"FORCE_MERGE_TABLE_CORRUPT({suffix})") ]
        flag_list = [ FLAG_FORCE_MERGE_TABLE_MISSING , 
                      FLAG_FORCE_MERGE_TABLE_CORRUPT  ]

        nkey = len(key_list)

        for ikey in range(0,nkey):
            tmp_key  = key_list[ikey]
            tmp_flag = flag_list[ikey]
            if tmp_key in CONFIG :
                for item in CONFIG[tmp_key] :
                    if item == version_fitopt:
                        flag = tmp_flag
                        msg = (f"force {suffix}-table fail for {version_fitopt}")
                        logging.info(msg)
        return flag
        # end force_merge_table_fail

    def get_merge_COLNUM_CPU(self):
        return COLNUM_FIT_MERGE_CPU
