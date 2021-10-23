# Base class Program

import os, sys, shutil, yaml
import logging   # , coloredlogs
import datetime, time, subprocess
import getpass, ntpath, glob

import numpy as np
import makeDataFiles_util as util
import write_data_snana   as snana

try:
    import write_data_lsst_alert as lsst_alert
except ImportError:
    pass

from   makeDataFiles_params  import *

import numpy as np
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
        outdir_list = [ args.outdir_snana, args.outdir_lsst_alert ]
        for outdir in outdir_list :
            if outdir is None: continue
            if not os.path.exists(outdir):
                logging.info(f" Create top-dir for light curves: {outdir}")
                sys.stdout.flush()
                os.mkdir(outdir)

        # - - - - - - - - 
        self.extend_DATAKEY_LIST(config_inputs)

        # store info for phot varnames
        self.store_varlist_obs(config_inputs, config_data)
    
        # end Program __init__
        
    def init_data_unit(self, config_inputs, config_data ):

        # define every possible data unit here and store them in list.
        # Only units with data will have a directory created.
        # The name in each list is a name that will be part of the
        # folder name.

        args          = self.config_inputs['args'] # user command line args
        nsplit        = args.nsplitran
        isplit_select = args.isplitran # 1 to nsplit, or -1 for all
        iyear_select  = args.year      # 1-NYEAR, or -1 for all
        field         = args.field
        survey        = args.survey
        peakmjd_range    = args.peakmjd_range
        mjd_detect_range = args.mjd_detect_range

        n_season       = MXSEASON

        # for MJD-related cuts, set n_season=1 so that there is
        # no explicit season breakdown
        if peakmjd_range    is not None: n_season = 1 
        if mjd_detect_range is not None: n_season = 1

        unit_name_list   = []
        unit_nevent_list = []
        msgerr    = []
        
        if isplit_select == 0 or isplit_select > nsplit:
            msgerr = []
            msgerr.append(f"Invalid --isplitran {isplit_select}")
            msgerr.append(f"Valid --isplitran arg range is 1 to {nsplit}")
            util.log_assert(False,msgerr)

        for iseason in range(0,n_season):
            iyear  = iseason + 1    # starts at 1
            if iyear_select > 0 and iyear != iyear_select :
                continue
            
            for isplit in range(0,nsplit):
                ISPLIT  = -9
                if nsplit > 1 : ISPLIT = isplit + 1

                if isplit_select > 0 and ISPLIT != isplit_select :
                    continue ;

                # define unit name for all seasons combined
                if (isplit == 0 or isplit_select>0) and iseason==0 :
                    unit_name = \
                        self.assign_data_unit_name(survey, field, -1, ISPLIT)
                    unit_name_list.append(unit_name)
                
                # define unit name for this season/iyear
                unit_name = \
                    self.assign_data_unit_name(survey, field, iyear, ISPLIT)
                unit_name_list.append(unit_name)                
        
        # init 'exist' logical to false for each data unit
        n_data_unit      = len(unit_name_list)

        unit_nevent_list = [ 0 ] * n_data_unit
        
        config_data['data_folder_prefix']      = survey
        config_data['data_unit_name_list']     = unit_name_list
        config_data['data_unit_nevent_list']   = unit_nevent_list
        config_data['n_season']                = n_season

        readme_stats_list = []
        for i in range(0,n_data_unit):
            readme_stats_list.append(util.init_readme_stats())

        config_data['readme_stats_list'] = readme_stats_list
        config_data['NEVT_SPECTRA']      = 0
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
        n_season      = self.config_data['n_season']
        args          = self.config_inputs['args'] # user command line args
        nsplit        = args.nsplitran
        isplit_select = args.isplitran  # 1 to nsplit, or -1 for all
        iyear_select  = args.year       # 1 to nyear, or -1 for all
        field         = args.field
        survey        = args.survey
        peakmjd_range    = args.peakmjd_range
        mjd_detect_range = args.mjd_detect_range

        data_unit_name = None
        d_raw    = data_dict['head_raw']
        d_calc   = data_dict['head_calc']

        SNID       = d_raw[DATAKEY_SNID]
        RA         = d_raw[DATAKEY_RA]
        DEC        = d_raw[DATAKEY_DEC]
        FIELD      = d_raw[DATAKEY_FIELD]
        PEAKMJD    = d_calc[DATAKEY_PEAKMJD] 
        MJD_DETECT = d_calc[DATAKEY_MJD_DETECT] 

        # check field match
        match_field = False
        if field == FIELD       : match_field = True
        if field == FIELD_VOID  : match_field = True
        if not match_field : return None

        # create dictionary needed to determine iyear
        small_event_dict = { 'peakmjd': PEAKMJD,  'mjd_detect': MJD_DETECT,
                             'ra': RA,  'dec': DEC,  'field': field }

        if n_season > 1:
            YY = util.iyear_survey(survey, small_event_dict)
        else:
            YY = -1  # no explicit season dependence

        # check year/season match
        match_year = True
            
        if iyear_select > 0 : 
            match_year = (YY == iyear_select)

        if not match_year : 
            return None

        # - - - - - - - - - - - - - 
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
            msgerr.append(f"RA={RA}  DEC={DEC}  PEAKMJD={PEAKMJD}")
            msgerr.append(f"Valid data_unit_name_list = ")
            msgerr.append(f"    {data_unit_name_list}")
            util.log_assert(False,msgerr)

        # - - - - -
        return data_unit_name
        # end which_data_unit
    
    def extend_DATAKEY_LIST(self,config_inputs):

        # expand global DATAKEY_LIST_RAW to include filter-dependent
        # HOSTGAL keys.

        survey   = config_inputs['args'].survey
        filters  = list(SURVEY_INFO['FILTERS'][survey])

        global DATAKEY_LIST_RAW
        prefix_list = [ HOSTKEY_PREFIX_MAG, HOSTKEY_PREFIX_MAGERR, 
                        HOSTKEY_PREFIX_SB ]
        for prefix in prefix_list :
            for band in filters:
                datakey = f"{prefix}_{band}"
                DATAKEY_LIST_RAW.append(datakey)

        # end load_HOSTKEY_band

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
        # Also count how many spectra and append to data_event_dict

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
            
        # check if there are spectra
        if 'spec_raw' in data_event_dict :
            spec_raw    = data_event_dict['spec_raw']
            n_spectra   = len(spec_raw)
        else:
            # if read source ignores spectra, add empty dictionary
            # to avoid crash later
            data_event_dict['spec_raw'] = {}
            n_spectra = 0

        data_event_dict['n_spectra'] = n_spectra

        return
        # compute_data_event

    def pass_data_cuts(self,data_event_dict):

        # Apply optional user cuts from command line
        # Return pass_cuts = True of False.

        pass_cuts = True  # init output to True

        args             = self.config_inputs['args']
        peakmjd_range    = args.peakmjd_range
        mjd_detect_range = args.mjd_detect_range
        #d_raw         = data_event_dict['head_raw']
        d_calc        = data_event_dict['head_calc']

        if peakmjd_range is not None:
            PEAKMJD       = d_calc[DATAKEY_PEAKMJD] 
            if PEAKMJD < peakmjd_range[0] : pass_cuts = False
            if PEAKMJD > peakmjd_range[1] : pass_cuts = False

        if mjd_detect_range is not None:
            MJD_DETECT    = d_calc[DATAKEY_MJD_DETECT] 
            if MJD_DETECT < mjd_detect_range[0]: pass_cuts = False
            if MJD_DETECT > mjd_detect_range[1]: pass_cuts = False
            
        return pass_cuts
        # end pass_data_cuts

    def init_phot_dict(self,NOBS):
        
        # The read_event function for each source should call this
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

    def init_spec_dict(self,NSPEC):
        spec_dict = {}
        spec_dict['NSPEC'] = NSPEC
        return spec_dict
        
    def final_summary(self):
    
        # comput total number of events and number of data units created
        NEVT_TOT  = 0
        NUNIT_TOT = 0
        nevent_list   = self.config_data['data_unit_nevent_list']
        for nevent in nevent_list:
            if nevent == 0 : continue 
            NEVT_TOT  += nevent
            NUNIT_TOT += 1

        # get total process time
        t_start = self.config_data['t_start']
        t_end   = datetime.datetime.now()  # xxx self.config_data['t_end']
        t_dif_sec  = (t_end-t_start).total_seconds()

        if t_dif_sec < 2000.0:
            t_dif   = t_dif_sec/60.0
            t_unit  = "minutes"
        else:
            t_dif  = t_dif_sec/3600.0
            t_unit = "hours"


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
            logging.info(f"\t\t Process evt={evt+1:6d} of {NEVT_TOT} "
                         f" ({rate}/sec)")
            sys.stdout.flush()
            
        return
    
    # end screen_update

    def read_data_driver(self):

        args = self.config_inputs['args']  # command line args

        # one-time init
        self.init_read_data()

        NEVT_READ    = 0

        nevent_subgroup = 1  # anything > 0 to pass while block below
        i_subgroup      = 0  # start with this subbroup index
        
        data_unit_name_list   = self.config_data['data_unit_name_list']

        while nevent_subgroup > 0 :
            
            nevent_subgroup = self.prep_read_data_subgroup(i_subgroup)
            if nevent_subgroup == 0 : break

            self.screen_update(-1,0)  # init clock for rate monitor

            for evt in range(0,nevent_subgroup):

                NEVT_READ += 1

                # call read-source-dependent function to read event
                data_event_dict = self.read_event(evt)
                
                # add computed variables; e.g., zCMB, MWEBV ...
                self.compute_data_event(data_event_dict)

                # apply optional user-cuts (e.g., PEAKMJD cut, etc...)
                pass_cuts = self.pass_data_cuts(data_event_dict)
                if pass_cuts is False:
                    continue

                # figure out which data unit
                data_unit_name = self.which_data_unit(data_event_dict)
                if data_unit_name is None : 
                    continue

                # add more info to data event dictionary
                index_unit  = data_unit_name_list.index(data_unit_name)
                data_event_dict['data_unit_name'] = data_unit_name
                data_event_dict['index_unit']     = index_unit

                if args.outdir_snana is not None :
                    snana.write_event_text_snana(args, self.config_data,
                                                 data_event_dict) 
                if args.outdir_lsst_alert is not None:
                    lsst_alert.write_event_lsst_alert(args, self.config_data,
                                                      data_event_dict) 
                # increment number of events for this data unit

                self.config_data['data_unit_nevent_list'][index_unit] += 1

                self.update_readme_stats(data_event_dict)

                self.screen_update(evt,nevent_subgroup)
                
                if NEVT_READ >= args.nevt : break

            self.end_read_data_subgroup()
            if NEVT_READ >= args.nevt : break
            i_subgroup += 1
            
        # - - - - -
        # option end-read tasks; e.g., close data base connection
        self.end_read_data()
        
        # write auxilary files for each data unit (name); 
        # e.g.  LIST and README files.
        # Use aux_file util based on choice of output format.
        nevent_list   = self.config_data['data_unit_nevent_list']
        name_list     = self.config_data['data_unit_name_list']
        for nevent, name in zip(nevent_list, name_list):
            if nevent == 0 : continue
            index_unit   = data_unit_name_list.index(name)
            if args.output_yaml_file:
                self.write_yaml_file(index_unit)
            if args.outdir_snana:
                snana.write_aux_files_snana(name, args, self.config_data)
            elif args.outdir_lsst_alert:
                pass # ???

        # end read_data_driver
    
    def write_yaml_file(self, index_unit):
        # write yaml file to be parsed by pipeline.
        # This is the same file as README file in output directory.
        args         = self.config_inputs['args'] 
        readme_stats = self.config_data['readme_stats_list'][index_unit]
        readme_dict = {
            'readme_file'  : args.output_yaml_file,
            'readme_stats' : readme_stats,
            'data_format'  : FORMAT_TEXT
        }
        util.write_readme(args,readme_dict)

        # end write_yaml_file

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

        # end reset_data_event_dict

    def update_readme_stats(self, data_event_dict):
        
        head_raw   = data_event_dict['head_raw']
        head_calc  = data_event_dict['head_calc']
        n_spectra  = data_event_dict['n_spectra'] 
        index_unit = data_event_dict['index_unit']
        
        # update stats that will eventually written to README file
        specz = -9.0;  photoz = -9.0
        if HOSTKEY_SPECZ  in head_raw:   specz   = head_raw[HOSTKEY_SPECZ]
        if HOSTKEY_PHOTOZ in head_calc:  photoz  = head_calc[HOSTKEY_PHOTOZ]

        readme_stats = self.config_data['readme_stats_list'][index_unit]

        readme_stats['NEVT_ALL'] += 1

        if specz > 0.0 :
            readme_stats['NEVT_HOSTGAL_SPECZ'] += 1

        if photoz > 0.0 :
            readme_stats['NEVT_HOSTGAL_PHOTOZ'] += 1

        if n_spectra > 0 : 
            key = 'NEVT_SPECTRA'
            readme_stats[key]      += 1  # increment this data unit
            self.config_data[key]  += 1  # sum over all data units

        # end update_readme_stats

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
    
