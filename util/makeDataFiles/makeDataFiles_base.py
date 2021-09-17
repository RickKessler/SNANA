# Base class Program

import os, sys, shutil, yaml
import logging, coloredlogs
import datetime, time, subprocess
import getpass, ntpath, glob

import numpy as np
import makeDataFiles_util as util

from   makeDataFiles_params  import *
from   abc import ABC, abstractmethod

# ======================================
class Program:
    def __init__(self, config_inputs, config_data):

        self.config_inputs = config_inputs
        self.config_data   = config_data
        args = config_inputs['args']
        
        print("  Base init")
        sys.stdout.flush()

        self.config_data['t_start'] = datetime.datetime.now()
                    
        # init all possible data units
        self.init_data_unit(config_inputs, config_data)

        # create top-level outdir
        outdir_list = [ args.outdir_snana ]
        for outdir  in outdir_list :
            if not os.path.exists(outdir):
                logging.info(f" Create top-dir for light curves: {outdir}")
                sys.stdout.flush()
                os.mkdir(outdir)

        # store info for phot varnames
        self.store_varlist_obs(config_inputs, config_data)
        
        # end Program __init__
        
    def init_data_unit(self, config_inputs, config_data ):

        # define every possible data unit here and store them in list.
        # Only units with data will have a directory created.
        # The name in each list is a name that will be part of the
        # folder name.

        args      = self.config_inputs['args'] # user command line args
        nseason       = 10
        nsplit        = args.nsplitran
        isplit_select = args.isplitran # 1 to nsplit, or -1 for all
        iyear_select  = args.year      # 1-NYEAR, or -1 for all
        field         = args.field
        survey        = args.survey
        
        unit_name_list   = []
        unit_nevent_list = []
        msgerr    = []
        
        for iseason in range(0,nseason):
            iyear  = iseason + 1    # starts at 1
            if iyear_select > 0 and iyear != iyear_select :
                continue
            
            for isplit in range(0,nsplit):
                ISPLIT  = -9
                if nsplit > 1 : ISPLIT = isplit + 1

                if isplit_select > 0 and ISPLIT != isplit_select :
                    continue ;

                if isplit == 0 and iseason==0 :
                    unit_name = \
                        self.assign_data_unit_name(survey, field, -1, ISPLIT)
                    unit_name_list.append(unit_name)
                
                unit_name = \
                    self.assign_data_unit_name(survey, field, iyear, ISPLIT)
                unit_name_list.append(unit_name)

                
        #sys.exit(f"\n xxx unit_name_list = {unit_name_list}")
        
        # init 'exist' logical to false for each data unit
        n_data_unit      = len(unit_name_list)

        unit_nevent_list = [ 0 ] * n_data_unit
        
        config_data['data_folder_prefix']  = survey
        config_data['data_unit_name_list']     = unit_name_list
        config_data['data_unit_nevent_list']   = unit_nevent_list

        return
    
        # end init_data_unit

    def assign_data_unit_name(self, survey, field, iseason, iran):

        name = f"{survey}"

        if field != FIELD_VOID : name += f"_{field}"

        if iseason >= 0 : name += f"_{PREFIX_SEASON}{iseason:02d}"
            
        if iran > 0:  name += f"_{PREFIX_SPLIT}{iran:03d}"
            
        return name
        # end assign_data_unit_name
        
    def which_data_unit(self, data_dict):

	# use data header info to figure out which data unit.
        # If no data unit is matched, return None

        args          = self.config_inputs['args'] # user command line args
        nsplit        = args.nsplitran
        isplit_select = args.isplitran  # 1 to nsplit, or -1 for all
        iyear_select  = args.year       # 1 to nyear, or -1 for all
        field         = args.field
        survey        = args.survey
        
        data_unit_name = None
        d_raw    = data_dict['head_raw']
        d_calc   = data_dict['head_calc']

        SNID    = d_raw[DATAKEY_SNID]
        RA      = d_raw[DATAKEY_RA]
        DEC     = d_raw[DATAKEY_DEC]
        FIELD   = d_raw[DATAKEY_FIELD]
        MJD     = d_calc[DATAKEY_PEAKMJD] # later should use MJD_TRIGGER

        # create dictionary needed to determine iyear
        event_dict = { 'mjd':MJD,  'ra': RA, 'dec':DEC, 'field':field }
        YY = util.iyear_survey(survey, event_dict)
                
        # check field match
        match_field = False
        if field == FIELD       : match_field = True
        if field == FIELD_VOID  : match_field = True
        if not match_field : return None

        # check year/season match
        match_year = True
        if iyear_select > 0 : match_year = (YY == iyear_select)
        if not match_year : return None
        
        # check split job    
        match_split = True ; ISPLIT = -9
        if nsplit > 1 :
            iSNID  = int(SNID)
            isplit = iSNID % nsplit       # counter starts at 0
            ISPLIT = isplit + 1           # counter starts at 1
            if isplit_select > 0 :
                match_split = (ISPLIT == isplit_select)
        if not match_split : return None

        # - - - - - - - -        
        data_unit_name  =  self.assign_data_unit_name(survey,FIELD, YY, ISPLIT)
        data_unit_name_list = self.config_data['data_unit_name_list']
        
        if data_unit_name not in data_unit_name_list :
            msgerr = []
            msgerr.append(f"Invalid data_unit_name = {data_unit_name}")
            msgerr.append(f"for SNID = {SNID} .")
            msgerr.append(f"RA={RA}  DEC={DEC}  MJD={MJD}")
            msgerr.append(f"Valid data_unit_name_list = ")
            msgerr.append(f"    {data_unit_name_list}")
            util.log_assert(False,msgerr)

        # - - - - -
        return data_unit_name
	# end which_data_unit
    
    def output_data_folder_name(self, data_unit_name, ISTEXT):
        prefix              = self.config_data['data_folder_prefix']
        data_unit_name_list = self.config_data['data_unit_name_list']
        
        if data_unit_name not in data_unit_name_list:
            msgerr = []
            msgerr.append(f" Invalid data unit '{data_unit_name}")
            msgerr.append(f" Valid data units are : ")
            msgerr.append(f"   {data_unit_list}")
            util.log_assert(False,msgerr)
            
        folder = f"{data_unit_name}"
        if ISTEXT:
            folder = f"{FORMAT_TEXT}_{folder}"

        return folder
    
        # end output_data_folder_name

    def store_varlist_obs(self, config_inputs, config_data):

        # for fakes, tack on true mag to list of variables per obs
        fake      = config_inputs['args'].fake
        survey    = config_inputs['args'].survey
        varnames  = VARNAMES_OBS
        varfmt    = VARNAMES_FMT
        val_undef = VAL_UNDEFINED_LIST
        
        if fake :
            varnames   += f" {VARNAME_TRUEMAG}"
            varfmt     += f" 8.4f"
            val_undef  += VAL_NULL 
            
        # convert space-sep string list into python list
        varlist_obs   = varnames.split()
        varlist_fmt   = varfmt.split()    # format per var
        vallist_undef = val_undef         # already a list

        # check for unfilled PHOT columns to exclude
        for var in self.exclude_varlist_obs() :
            k          = varlist_obs.index(var)
            temp_obs   = varlist_obs.pop(k)
            temp_fmt   = varlist_fmt.pop(k)
            temp_undef = vallist_undef.pop(k)                
                
        nvar    = len(varlist_obs)
        config_data['vallist_undef']  = vallist_undef
        config_data['varlist_fmt']   = varlist_fmt
        config_data['varlist_obs']   = varlist_obs
        config_data['nvar_obs']      = nvar
                    
        #print(f" xxx nvar={nvar}  varlist = {varlist}")

        return
    
        # end store_varlist_obs

    def exclude_varlist_obs(self):
        # return optional list of phot columns to exclude from the
        # output text data files; default is exclude nothing and
        # write out all columns
        return []
    
    def compute_data_event(self, data_event_dict):
        
        # compute & append a few varaibles to
        #   data_event_dict['head_raw']
        #   data_event_dict['head_calc'] 

        msgerr   = []
        fake     = self.config_inputs['args'].fake
        survey   = self.config_inputs['args'].survey
        d_raw    = data_event_dict['head_raw']
        d_calc   = data_event_dict['head_calc']        
        
        snid     = d_raw[DATAKEY_SNID]
        zhel     = d_raw[DATAKEY_zHEL]
        zhel_err = d_raw[DATAKEY_zHEL_ERR]
        ra       = d_raw[DATAKEY_RA]
        dec      = d_raw[DATAKEY_DEC]

        if fake :
            snana_flag_fake = SNANA_FLAG_FAKE
        else:
            snana_flag_fake = SNANA_FLAG_DATA
            
        zcmb      = util.helio_to_cmb(zhel, ra, dec)

        # no urgency for loading MWEBV because TEXT->FITS translator
        # computes and stores MWEBV. However, if we want correct MWEBV
        # in the TEXT files, need to compute it here:

        if DATAKEY_MWEBV in d_calc:
            mwebv     = d_calc[DATAKEY_MWEBV]
            mwebv_err = d_calc[DATAKEY_MWEBV_ERR]
        else:
            mwebv     = -9.0
            mwebv_err = -9.0

        # - - - -
        dump_flag = False
        if dump_flag :
            print(f" xxx ------------------------------")
            print(f" xxx DUMP for compute_data_event")
            print(f" xxx SNID={snid}   RA={ra}  DEC={dec}  zhel={zhel:8.5f}")
            if DATAKEY_zCMB in d_calc:
                zcmb_deja = d_calc[DATAKEY_zCMB]
                print(f"\t already existing zcmb = {zcmb_deja:8.5f}")
            if DATAKEY_MWEBV in d_calc:
                mwebv_deja     = d_calc[DATAKEY_MWEBV]
                mwebv_deja_err = d_calc[DATAKEY_MWEBV_ERR]
                print(f"\t already existing mwebv = " \
                      f"{mwebv_deja:8.5f} +_ {mwebv_deja_err:8.5f} ") 
                
            print(f" xxx COMPUTE zcmb  = {zcmb:8.5f}")
            print(f" xxx COMPUTE mwebv = {mwebv:8.5f} +_ {mwebv:8.5f}")
            sys.stdout.flush()
            
        # - - - - - - -
        # load goodies
        d_raw[DATAKEY_SURVEY]      = survey
        d_raw[DATAKEY_FAKE]        = snana_flag_fake

        if survey not in SURVEY_INFO['FILTERS']:
            msgerr.append(f"{survey} filters not defined")
            msgerr.append(f"Check SURVEY_INFO dictionary in " \
                          f"makeDataFiles_params.py")
            util.log_assert(False,msgerr)
        else:
            d_raw[DATAKEY_FILTERS]     = SURVEY_INFO['FILTERS'][survey]

        if survey in SURVEY_INFO['CCD']:
            d_raw[DATAKEY_NXPIX]   = SURVEY_INFO['CCD'][survey][0]
            d_raw[DATAKEY_NYPIX]   = SURVEY_INFO['CCD'][survey][1]
            d_raw[DATAKEY_PIXSIZE] = SURVEY_INFO['CCD'][survey][2]
                
        d_calc[DATAKEY_zCMB]      = zcmb
        d_calc[DATAKEY_zCMB_ERR]  = zhel_err
        d_calc[DATAKEY_MWEBV]     = mwebv
        d_calc[DATAKEY_MWEBV_ERR] = mwebv_err

        # if there is no VPEC, tack on default
        if DATAKEY_VPEC not in d_calc :
            d_calc[DATAKEY_VPEC]     = VPEC_DEFAULT[0]
            d_calc[DATAKEY_VPEC_ERR] = VPEC_DEFAULT[1]
            
        return
        # compute_data_event

    def init_phot_dict(self,NOBS):
        
        # The read_event function for each class should call this
        # function before reading event photometry.
        # Init value = None for each column and observation.
        # If a PHOT variable isn't set, base code sees None value
        # and can take appropriate action (abort, set to -9, etc ...)

        phot_dict   = {}
        phot_dict['NOBS'] = NOBS

        varlist_obs = self.config_data['varlist_obs']
        for varname in varlist_obs:
            phot_dict[varname] = [ None ] * NOBS
        
        return phot_dict
        # end init_phot_dict

    def write_aux_files_snana(self):

        # write auxilary files (.LIST and .README) for each created folder
        # with TEXT formatted files.
        # Here we loop over folders created from this task.
        
        args          = self.config_inputs['args']
        outdir        = args.outdir_snana
        nevent_list   = self.config_data['data_unit_nevent_list']
        name_list     = self.config_data['data_unit_name_list']
        prefix        = self.config_data['data_folder_prefix']

        print(f"")
        sys.stdout.flush()
        
        for nevent, name in zip(nevent_list, name_list):

            if nevent == 0 : continue
            folder_out    = self.output_data_folder_name(name,True)

            msg = f" Create aux files for {folder_out} and gzip " \
                  f"{TEXTFILE_SUFFIX} files."
            logging.info(msg)
            sys.stdout.flush()
            
            data_dir      = f"{outdir}/{folder_out}"
            search_string = f"{prefix}*{TEXTFILE_SUFFIX}"

            data_file_list = glob.glob1(data_dir, f"{search_string}" )
            list_file      = f"{data_dir}/{folder_out}.LIST"
            readme_file    = f"{data_dir}/{folder_out}.README"        
        
            with open(list_file,"wt") as f:
                for data_file in data_file_list:
                    f.write(f"{data_file}\n")

            self.write_readme(folder_out, nevent, FORMAT_TEXT)
            
            # gzip data files
            cmd = f"cd {data_dir} ; gzip {search_string}"
            os.system(cmd)
        
        # end write_aux_files_snana
        
    def write_event_text_snana(self, data_event_dict, data_unit_name):

        # create one text file and write one event described by
        # dictionary data_event_dict.
        # Input data_unit_name is used to determine name of folder.
        
        ISTEXT = True
        outdir                = self.config_inputs['args'].outdir_snana
        data_unit_name_list   = self.config_data['data_unit_name_list']
        data_unit_nevent_list = self.config_data['data_unit_nevent_list'] 
        prefix                = self.config_data['data_folder_prefix']
        indx_unit      = data_unit_name_list.index(data_unit_name)  
        nevent         = data_unit_nevent_list[indx_unit]
        folder         = self.output_data_folder_name(data_unit_name,ISTEXT)
        data_dir       = f"{outdir}/{folder}"
        tar_file       = f"{outdir}/{folder}.tar.gz"
        exist_folder   = os.path.exists(data_dir)
        exist_tar_file = os.path.exists(tar_file)
        
        if nevent == 0 :
            # remove folder if it already exists
            if exist_folder :
                cmd_rm = f"rm -r {data_dir}"
                os.system(cmd_rm)

            if exist_tar_file :
                cmd_rm = f"rm -r {tar_file}"
                os.system(cmd_rm)
                
            # create folder
            logging.info(f"\t Create folder {folder}")
            sys.stdout.flush()
            os.mkdir(data_dir)
            
        head_raw  = data_event_dict['head_raw']
        head_calc = data_event_dict['head_calc']
        phot_raw  = data_event_dict['phot_raw']
        
        SNID      = head_raw[DATAKEY_SNID]
        NOBS      = phot_raw[DATAKEY_NOBS]

        str_SNID     = SNID
        if SNID.isdigit(): str_SNID = f"{int(SNID):010d}"
        
        data_file     = f"{data_dir}/{prefix}_{str_SNID}.DAT"
        nvar_obs      = self.config_data['nvar_obs']
        varlist_obs   = self.config_data['varlist_obs']
        varlist_fmt   = self.config_data['varlist_fmt']
        vallist_undef = self.config_data['vallist_undef'] 
        varstring_obs = ' '.join(varlist_obs)
        msgerr = []

        self.config_data['data_unit_nevent_list'][indx_unit] += 1
        
        with open(data_file, "wt") as f :

            # write header info
            self.write_header_snana(f,head_raw)
            
            f.write("\n# computed quantities \n")
            self.write_header_snana(f,head_calc)

            # write epoch/phot info
            f.write(f"\n# -------------------------------------- \n" \
                    f"# obs info\n")
            f.write(f"NOBS: {NOBS}\nNVAR: {nvar_obs} \n"
                    f"VARLIST: {varstring_obs}\n")
            
            for obs in range(0,NOBS):
                LINE = "OBS:"
                for varname,fmt,val_undef in \
                    zip(varlist_obs,varlist_fmt,vallist_undef):
                    val = phot_raw[varname][obs]
                    if val == None:
                        val = val_undef
                        if val == VAL_ABORT :
                            msgerr.append(f"Missing required PHOT column {varname}")
                            msgerr.append(f"Check SNID = {SNID}")
                            util.log_assert(False,msgerr)
                            
                    #print(f" xxx load varnamne={varname} val={val} fmt={fmt}")
                    LINE += f" {val:{fmt}}"
                f.write(f"{LINE}\n")

            f.write(f"END:\n")
        # end write_event_text_snana

    def write_header_snana(self, f, data_head):

        # write list of header key,val pairs in data_head.
        # If XXX and XXX_ERR both exist, write as
        #   KEYXXX:  XXX +_ XXX_ERR
        
        for key in data_head:
            #print(f" xxx header key = {key}")
            if '_ERR' in key: continue
            key_plus_err   = f"{key}_ERR"
            key_plus_colon = f"{key}:"
            val            = data_head[key]
            string_val     = f"{val}"
            if val == VAL_NULL : continue
            if key_plus_err in data_head:
                err  = data_head[key_plus_err]
                string_val = f"{val:9.6f} +- {err:9.6f}"
            f.write(f"{key_plus_colon:<20s}  {string_val} \n")
            f.flush()
            
        # and write_head_keys

    def convert2fits_snana(self):

        # loop over newly created TEXT file versions and convert
        # to fits format ... then tar up TEXT folder.

        args          = self.config_inputs['args']
        outdir        = args.outdir_snana
        text          = args.text 
        nevent_list   = self.config_data['data_unit_nevent_list']
        name_list     = self.config_data['data_unit_name_list']
        prefix        = self.config_data['data_folder_prefix']

        print(f"")
        sys.stdout.flush()
        
        NEVT_TOT  = 0
        NUNIT_TOT = 0
        for nevent,name in zip(nevent_list, name_list):
            if nevent == 0 : continue
            folder_text    = self.output_data_folder_name(name,True)
            folder_fits    = self.output_data_folder_name(name,False)
            log_file       = f"{folder_text}/convert2fits_{folder_fits}.log"
            yaml_file      = f"{outdir}/{folder_text}.YAML"  # expected output
            
            NEVT_TOT  += nevent
            NUNIT_TOT += 1
            
            msg = f"  Convert TEXT -> FITS for {folder_fits}" \
                  f" NEVT={nevent}"
            logging.info(msg)
            sys.stdout.flush()
            
            time_0 = datetime.datetime.now()        
            outdir_text  = f"{outdir}/{folder_text}"
            outdir_fits  = f"{outdir}/{folder_fits}"

            # rm fits folder if still there from previous job
            if os.path.exists(outdir_fits):
                cmd_rm = f"cd {outdir} ; rm -r {folder_fits}"
                os.system(cmd_rm)            

            cmd_snana   = f"{PROGRAM_SNANA} NOFILE " \
                          f"PRIVATE_DATA_PATH ./ " \
                          f"VERSION_PHOTOMETRY    {folder_text} " \
                          f"VERSION_REFORMAT_FITS {folder_fits} " \
                          f"{OPTIONS_TEXT2FITS_SNANA} "
            cmd = f"cd {outdir}; {cmd_snana} > {log_file}"
            os.system(cmd)

            # - - - - 
            # if YAML file doesn't exist, abort with message that
            # convert job probably aborted or crashed.
            if not os.path.exists(yaml_file):
                msgerr = []
                msgerr.append(f"Cannot find expected yaml file:")
                msgerr.append(f"    {yaml_file}")
                msgerr.append(f"TEXT->FITS convert job probably aborted or crashed;")
                msgerr.append(f"See convert-log file:")
                msgerr.append(f"    {outdir}/{log_file} ")
                util.log_assert(False,msgerr)

            # - - - - - 
            # clean up
            # gzip FITS files and make compressed tar file from TEXT dir
            cmd_gzip_fits = f"cd {outdir_fits} ; gzip *.FITS"

            tar_file = f"{folder_text}.tar"
            cmd_tar_text  = f"cd {outdir} ; " \
                            f"tar -cf {tar_file} {folder_text} ; " \
                            f"gzip {tar_file} ; " \
                            f"rm -r {folder_text} "
            
            os.system(cmd_gzip_fits)
            if not text:
                os.system(cmd_tar_text)

            # re-write FITS readme in data folder
            self.write_readme(folder_fits, nevent, FORMAT_FITS )
            
            time_1   = datetime.datetime.now()
            time_dif = (time_1 - time_0).total_seconds()
            rate     = int(float(nevent)/float(time_dif))
            logging.info(f"\t Rate(convert+cleanup): {rate}/sec ")
            sys.stdout.flush()

        # - - - - -
        self.config_data['t_end']     = datetime.datetime.now()
        self.config_data['NEVT_TOT']  = NEVT_TOT
        self.config_data['NUNIT_TOT'] = NUNIT_TOT

        return
    
        # end convert2fits_snana

    def write_readme(self, folder, nevent, data_format):

        # write readme file for SNANA formatted data folder.
        # data_format is either TEXT or FITS ... for FITS, include
        # extra info from YAML file created by snana.exe.
        
        args          = self.config_inputs['args']
        outdir        = args.outdir_snana
        outdir_data  = f"{outdir}/{folder}"
        readme_file  = f"{outdir_data}/{folder}.README"
        script_command = ' '.join(sys.argv)
        IS_FITS      = (data_format == FORMAT_FITS)
        
        # for FITS format, pick up extra info from YAML file created by snana.exe
        if IS_FITS :
            yaml_file       = f"{outdir}/{FORMAT_TEXT}_{folder}.YAML"  
            snana_yaml      = util.read_yaml(yaml_file)
            NEVT_HOST_ZSPEC = snana_yaml['NEVT_HOST_ZSPEC']
            NEVT_HOST_ZPHOT = snana_yaml['NEVT_HOST_ZPHOT']
                #print(f"\n xxx snana_yaml = \n{snana_yaml}")
            # remove YAML file
            cmd_rm = f"rm {yaml_file}"
            os.system(cmd_rm)
            
        with open(readme_file,"wt") as f:
            f.write(f"{DOCANA_KEY}: \n")
            f.write(f"  PURPOSE:  transient lightcurve data files " \
                        f"for analysis\n")

            if args.lsst_ap :
                f.write(f"  SOURCE_LSST_AP:   {args.lsst_ap} \n")
                
            if args.lsst_drp :
                f.write(f"  SOURCE_LSST_DRP:  {args.lsst_drp} \n")

            if args.sirah_folder is not None :
                f.write(f"  SOURCE_SIRAH_FOLDER:  {args.sirah_folder} \n")

            if args.snana_folder is not None:
                f.write(f"  SOURCE_SNANA_FOLDER:  {args.snana_folder} \n")
                
            f.write(f"  SURVEY:           {args.survey} \n")
            f.write(f"  FIELD:            {args.field} \n")
            f.write(f"  FORMAT:           {data_format} \n")
            f.write(f"  SCRIPT_COMMAND:   {script_command} \n")
            f.write(f"  USERNAME:         {USERNAME} \n")
            f.write(f"  HOSTNAME:         {HOSTNAME} \n")
            f.write(f"  NEVT_ALL:         {nevent} \n")

            if IS_FITS:
                f.write(f"  NEVT_HOST_ZSPEC:  {NEVT_HOST_ZSPEC} \n")
                f.write(f"  NEVT_HOST_ZPHOT:  {NEVT_HOST_ZPHOT} \n")                
                
            f.write(f"{DOCANA_KEY_END}: \n")
            
        # end write_readme
        
    def final_summary(self):
        
        # grand summary
        t_start = self.config_data['t_start']
        t_end   = self.config_data['t_end']
        t_dif_sec  = (t_end-t_start).total_seconds()

        if t_dif_sec < 2000.0:
            t_dif   = t_dif_sec/60.0
            t_unit  = "minutes"
        else:
            t_dif  = t_dif_sec/3600.0
            t_unit = "hours"

        NEVT_TOT  = self.config_data['NEVT_TOT']
        NUNIT_TOT = self.config_data['NUNIT_TOT']

        logging.info("\n FINAL MAKE-DATA-FILE SUMMARY: \n")
        logging.info(f" Total number of data units created: {NUNIT_TOT}")
        logging.info(f" Total number of events written:     {NEVT_TOT}")
        logging.info(f" Total processing time ({t_unit}):    {t_dif:.1f}" )
        sys.stdout.flush()
        
        # end final_summary

    def screen_update(self,evt,NEVT_TOT):

        if evt < 0 :
            time_0  = datetime.datetime.now()
            self.config_data['time_0'] =  time_0
            return

        if evt == 0 : return
        
        rmd = evt % NEVT_SCREEN_UPDATE
        if rmd == 0 or evt == NEVT_TOT-1 :
            time_0   = self.config_data['time_0']
            time_now = datetime.datetime.now()
            time_dif = (time_now - time_0).total_seconds()
            rate     = int(float(evt)/float(time_dif))
            logging.info(f"\t\t Process evt={evt:6d} of {NEVT_TOT} "
                         f" ({rate}/sec)")
            sys.stdout.flush()
            
        return
    
    # end screen_update

    def read_data_driver(self):

        args = self.config_inputs['args']  # command line args

        # one-time init
        self.init_read_data()

        NEVT_READ = 0
        NEVT_WRITE = 0
        nevent_subgroup = 1
        i_subgroup = 0
        
        while nevent_subgroup > 0 :
            
            nevent_subgroup = self.prep_read_data_subgroup(i_subgroup)
            if nevent_subgroup == 0 : break

            self.screen_update(-1,0)  # init clock for rate monitor

            for evt in range(0,nevent_subgroup):

                NEVT_READ += 1

                # call class-dependent function to read event
                data_event_dict = self.read_event(evt)

                # add computed variables; e.g., zCMB, MWEBV ...
                self.compute_data_event(data_event_dict)

                # figure out which data unit
                data_unit_name = self.which_data_unit(data_event_dict)
                if data_unit_name is None : continue

                NEVT_WRITE += 1

                self.write_event_text_snana(data_event_dict,data_unit_name)

                self.screen_update(evt,nevent_subgroup)
                
                if NEVT_READ >= args.nevt : break

            self.end_read_data_subgroup()
            if NEVT_READ >= args.nevt : break
            i_subgroup += 1
            
        # - - - - -
        self.end_read_data()

        # create LIST and README file 
        self.write_aux_files_snana()
            
        # end read_data_driver
    
    def reset_data_event_dict(self):

        # reset all data values to -9 to ensure that every
        # key gets written to data files, even if read_event
        # code fails to set a value.
        
        raw_dict  = {}
        calc_dict = {}
        
        for key in DATAKEY_LIST_RAW :
            raw_dict[key] = -9
                
        for key in DATAKEY_LIST_CALC :
            calc_dict[key] = -9            
            
        return raw_dict, calc_dict
    
    # ================================================
    #    MERGE PROCESS
    # ================================================

    def merge_snana_driver(self):

        args   = self.config_inputs['args']
        outdir = args.outdir_snana
        survey = args.survey
        
        if outdir is None:
            sys.exit(f"\n ERROR: must specify --outdir_snana\n")

        print(f"\n Merge {outdir}")

        # merge SPLIT folders; get all prefixes by scooping up all SPLIT001 jobs
        search_string = f"{survey}*{PREFIX_SPLIT}001"
        split_dir_list = sorted(glob.glob1(outdir, search_string ))
        for split_dir in split_dir_list:
            # if split_dir = LSST_WFDY01_SPLIT001, base_name=LSST_WFDY01
            merge_folder = split_dir.split(f"_{PREFIX_SPLIT}")[0]
            search_string = f"{merge_folder}_{PREFIX_SPLIT}*"
            self.merge_snana_folders(MODE_MERGE_MOVE,
                                     outdir, search_string, merge_folder)


        # merge Y## folders into folder with all seasons
        search_string = f"{survey}_*{PREFIX_SEASON}*"
        year_dir_list = sorted(glob.glob1(outdir, search_string ))
        merge_folder  = year_dir_list[0].split(f"_{PREFIX_SEASON}")[0]
        self.merge_snana_folders(MODE_MERGE_LINK,
                                 outdir, search_string, merge_folder)

        # archive TEXT versions 
        TEXT_archive_dir = "TEXT_archive"
        TEXT_list = glob.glob1(outdir, f"TEXT*" )
        if len(TEXT_list) > 0 :
            os.mkdir(f"{outdir}/{TEXT_archive_dir}")
            cmd_mv = f"cd {outdir} ; mv TEXT*.tar.gz {TEXT_archive_dir}"
            os.system(cmd_mv)

        return

    # end merge_snana_driver
    
    def merge_snana_folders(self, MODE, outdir, folder_list_string, merge_folder):
    
        # e.g., folder_list_string = LSST_WFD01_SPLIT*
        #       merge_folder       = LSST_WFD01
        #    ->
        #      combine data from all LSST_WFD01_SPLIT* folders into 
        #      single folder LSST_WFD01

        msgerr = []
        logging.info(f"\t Create merge-folder {merge_folder}")
        merge_folder_full = f"{outdir}/{merge_folder}"
        if os.path.exists(merge_folder_full):
            msgerr.append(f"{merge_folder_full} already exists ?!?!?!")
            msgerr.append(f"Something is wacky.")
            sys.exit(msgerr)

        os.mkdir(merge_folder_full)

        n_move = 0
        statsum_dict = {}
        for key in KEYLIST_README_STATS:   statsum_dict[key] = 0


        folder_list = glob.glob1(outdir, f"{folder_list_string}" )
        for folder in folder_list :
            n_move += 1
            FOLDER = f"{outdir}/{folder}"

            if MODE == MODE_MERGE_MOVE :
                mv_string = f"*.FITS.gz"
                cmd_mv = f"mv {mv_string} ../{merge_folder}"
                cmd    = f"cd {FOLDER}; {cmd_mv}"
            else:
                cmd = f"cd {outdir}/{merge_folder} ; "
                fits_file_list = glob.glob1(FOLDER, f"*.FITS.gz" )
                for fits_file in fits_file_list :
                    cmd += f"ln -s ../{folder}/{fits_file} {fits_file} ; "
            os.system(cmd)

            if n_move == 1 : 
                cmd_mv = f"mv *.IGNORE ../{merge_folder}/{merge_folder}.IGNORE"
                cmd    = f"cd {FOLDER}; {cmd_mv}" 
                os.system(cmd)

            # increment sum stats from readme
            README_file  = f"{FOLDER}/{folder}.README"
            README_yaml  = util.read_yaml(README_file)
            for key in KEYLIST_README_STATS:  
                NEVT  = README_yaml[DOCANA_KEY][key]
                statsum_dict[key] += NEVT

        # - - - - - - - -
        # update sum stats and re-write readme
        README_file = f"{merge_folder_full}/{merge_folder}.README"
        for key in KEYLIST_README_STATS:  
            NEVT = statsum_dict[key]
            README_yaml[DOCANA_KEY][key] = statsum_dict[key]
        self.write_merge_readme(README_file,README_yaml)

        # - - - - 
        # remove original folders
        if MODE == MODE_MERGE_MOVE:
            cmd_rm = f"rm -r {folder_list_string}"
            cmd    = f"cd {outdir}; {cmd_rm}"
            os.system(cmd)

        # create merged LIST file
        HEAD_list = glob.glob1(merge_folder_full, "*HEAD*.FITS.gz" )
        LIST_file = f"{merge_folder_full}/{merge_folder}.LIST"
        with open(LIST_file,"wt") as l :
            for item in HEAD_list:
                HEAD_file = item
                # strip off .gz extension
                if '.gz' in item: HEAD_file = item.split('.gz')[0]  
                l.write(f"{HEAD_file}\n")

        #print(f" xxx ---------------------------- ")    
        #print(f" xxx base_name = {base_name} ")
        #print(f" xxx split_dir_list = {split_dir_list} ")

        return
        # end merge_snana_folders

    def write_merge_readme(self, README_file, README_yaml):

        # write to README_file with contents README_yaml.
        with open(README_file,"wt") as r:
            yaml.dump(README_yaml,r, sort_keys=False)
            #r.write(f"{README_yaml}\n")

    # ====================================
    @abstractmethod
    def init_read_data(self):   # one-time init
        raise NotImplementedError()

    @abstractmethod
    def end_read_data(self):    # global end for reading data
        raise NotImplementedError()
    
    @abstractmethod
    def prep_read_data_subgroup(self): # prepare to read subgroup
        raise NotImplementedError()

    @abstractmethod
    def end_read_data_subgroup(self):  # end reading subgroup
        raise NotImplementedError()

    @abstractmethod
    def read_event(self):    # read one event; fill data_dict
        raise NotImplementedError()

    @abstractmethod
    def iyear_data(self,MJD,RA,DEC,FIELD):
        raise NotImplementedError()
    
