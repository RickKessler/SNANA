#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  makeDataFiles_base.py
#
#  Copyright 2021 R. Kessler
#
# April 9 2022 RK
#   Continue on nevt=0 and stop on nevt < 0 ; allows small samples
#   in which some of the FITS files have zero events.
#
"""Base meta-class Program, built to include the basic infrastructre for more
general data file parsers.
"""

import os
import sys
import logging
import datetime

# import numpy as np
# import yaml

from abc import abstractmethod


import makeDataFiles_params as gpar
import makeDataFiles_util as util
import write_data_snana 

from write_data_csv import csvWriter


# =============================================================================
#  Main Meta-Class Program
# =============================================================================
class Program:
    """Structural base meta-Class.
    This is the bare bones logic for the makedatafiles framework
    """

    def __init__(self, config_inputs, config_data):

        args = config_inputs['args']
        self.config_inputs = config_inputs
        self.config_data   = config_data

        logging.info("  Base init")

        self.config_data['t_start'] = datetime.datetime.now()

        # init all possible data units
        self.init_data_unit()

        # create top-level outdir
        outdir_list = [
            args.outdir_snana,
            args.outdir_csv
        ]

        for outdir in outdir_list:
            if outdir is None:
                continue

            if not os.path.exists(outdir):
                logging.info(f" Create top-dir for light curves: {outdir}")
                sys.stdout.flush()
                os.mkdir(outdir)


        self.extend_DATAKEY_LIST()

        # store info for phot varnames
        self.store_varlist_obs()


        # end Program __init__

    def init_data_unit(self):

        # define every possible data unit here and store them in list.
        # Only units with data will have a directory created.
        # The name in each list is a name that will be part of the
        # folder name.

        args          = self.config_inputs['args']  # user command line args
        nsplit        = args.nsplitran
        isplit_select = args.isplitran  # 1 to nsplit, or -1 for all
        iyear_select  = args.year       # 1-NYEAR, or -1 for all
        field_select  = args.field
        survey        = args.survey
        peakmjd_range      = args.peakmjd_range
        nite_detect_range  = args.nite_detect_range
        outdir_csv         = args.outdir_csv
        n_season           = gpar.MXSEASON

        # for MJD-related cuts, set n_season=1 so that there is
        # no explicit season breakdown
        if peakmjd_range is not None:
            n_season = 1
        if nite_detect_range is not None:
            n_season = 1
        if outdir_csv :
            n_season = 1

        unit_name_list   = []
        unit_nevent_list = []
        msgerr           = []

        if isplit_select == 0 or isplit_select > nsplit:
            msgerr = []
            msgerr.append(f"Invalid --isplitran {isplit_select}")
            msgerr.append(f"Valid --isplitran arg range is 1 to {nsplit}")
            util.log_assert(False,msgerr)

        # - - - - - -
        for iseason in range(0, n_season):
            iyear  = iseason + 1    # starts at 1
            if iyear_select > 0 and iyear != iyear_select:
                continue

            for isplit in range(0,nsplit):
                ISPLIT  = -9
                if nsplit > 1:
                    ISPLIT = isplit + 1

                if isplit_select > 0 and ISPLIT != isplit_select:
                    continue

                do_all_seasons = (isplit == 0 or isplit_select>0) and iseason==0
                do_one_season  = True  # (not outdir_lsst_alert)
                if outdir_csv : 
                    do_all_seasons = True
                    do_one_season  = False

                # define unit name for all seasons combined
                if do_all_seasons:
                    unit_name = self.assign_data_unit_name(
                        survey,
                        field_select,
                        -1,
                        ISPLIT
                    )
                    unit_name_list.append(unit_name)

                # define unit name for this season/iyear
                if do_one_season:
                    unit_name = \
                        self.assign_data_unit_name(survey, field_select,
                                                   iyear, ISPLIT)
                    unit_name_list.append(unit_name)

        # - - - - - - - - - - -
        # init 'exist' logical to false for each data unit
        n_data_unit      = len(unit_name_list)

        unit_nevent_list = [ 0 ] * n_data_unit

        self.config_data['data_folder_prefix']      = survey
        self.config_data['data_unit_name_list']     = unit_name_list
        self.config_data['data_unit_nevent_list']   = unit_nevent_list
        self.config_data['n_season']                = n_season

        #print(f"\n xxx unit_name_list = \n{unit_name_list}")

        readme_stats_list = []
        readme_stats_sum = util.init_readme_stats()
        for i in range(0, n_data_unit):
            readme_stats_list.append(util.init_readme_stats())

        self.config_data['readme_stats_list'] = readme_stats_list
        self.config_data['readme_stats_sum']  = readme_stats_sum
        self.config_data['NEVT_SPECTRA']      = 0
        return
        # end init_data_unit

    def assign_data_unit_name(self, survey, field, iseason, iran):

        name = f"{survey}"

        if field != gpar.FIELD_VOID:
            name += f"_{field}"

        if iseason >= 0:
            name += f"_{gpar.PREFIX_SEASON}{iseason:02d}"

        if iran > 0:
            name += f"_{gpar.PREFIX_SPLIT}{iran:03d}"

        return name
        # end assign_data_unit_name

    def which_data_unit(self, data_dict):

        # use data header info to figure out which data unit.
        # If no data unit is matched, return None
        n_season      = self.config_data['n_season']
        args          = self.config_inputs['args']  # user command line args
        nsplit        = args.nsplitran
        isplit_select = args.isplitran  # 1 to nsplit, or -1 for all
        iyear_select  = args.year       # 1 to nyear, or -1 for all
        field_select  = args.field
        survey        = args.survey
        peakmjd_range    = args.peakmjd_range
        nite_detect_range = args.nite_detect_range

        data_unit_name = None
        d_raw    = data_dict['head_raw']
        d_calc   = data_dict['head_calc']

        SNID       = d_raw[gpar.DATAKEY_SNID]
        RA         = d_raw[gpar.DATAKEY_RA]
        DEC        = d_raw[gpar.DATAKEY_DEC]
        FIELD      = d_raw[gpar.DATAKEY_FIELD]
        PEAKMJD    = d_calc[gpar.DATAKEY_PEAKMJD]
        MJD_DETECT = d_calc[gpar.DATAKEY_MJD_DETECT_FIRST]

        # - - - - - - - - - - - - - - - - -
        # check match for field
        #
        ### Could this block be just??:
        ### define the match_field and break:
        ## match_field = (field_select == gpar.FIELD_VOID) or (field_select == FIELD)
        ## if not match_field: return None
        
        ### or just
        ## if field_select!=gpar.FIELD_VOID and field_select != FIELD:
        ##    return None

        if field_select == gpar.FIELD_VOID:
            match_field = True  # no --field arg
        else:
            match_field = False
            if field_select == FIELD:
                match_field = True

        if not match_field :
            return None

        # - - - - - - - - -
        # check match for season/year

        YY = -1  # no explicit season dependence
        if n_season > 1:
        # create dictionary needed to determine iyear
            small_event_dict = {
                'peakmjd':PEAKMJD,
                'mjd_detect':MJD_DETECT,
                'ra':RA,
                'dec':DEC,
                'field':FIELD
            }
            YY = util.iyear_survey(survey, small_event_dict)

        match_year = True
        if iyear_select > 0:
            match_year = (YY == iyear_select)

        if not match_year:
            return None

        # - - - - - - - - - - - - -
        # check match for split job
        match_split = True
        ISPLIT = -9
        if nsplit > 1:
            iSNID  = int(SNID)
            isplit = iSNID % nsplit       # counter starts at 0
            ISPLIT = isplit + 1           # counter starts at 1
            if isplit_select > 0:
                match_split = (ISPLIT == isplit_select)
        if not match_split:
            return None

        # - - - - - - - -
        data_unit_name  = \
            self.assign_data_unit_name(survey, field_select, YY, ISPLIT)
        data_unit_name_list = \
            self.config_data['data_unit_name_list']

        if data_unit_name not in data_unit_name_list:
            msgerr = []
            msgerr.append(f"Invalid data_unit_name = {data_unit_name}")
            msgerr.append(f"for SNID = {SNID} .")
            msgerr.append(f"RA={RA}  DEC={DEC}  PEAKMJD={PEAKMJD}")
            msgerr.append(f"Valid data_unit_name_list = ")
            msgerr.append(f"    {data_unit_name_list}")
            util.log_assert(False, msgerr)

        # - - - - -
        return data_unit_name
        # end which_data_unit

    def extend_DATAKEY_LIST(self):

        # expand global DATAKEY_LIST_RAW to include filter-dependent
        # HOSTGAL keys.

        survey   = self.config_inputs['args'].survey
        filters  = list(gpar.SURVEY_INFO['FILTERS'][survey])

        for prefix in gpar.HOSTKEY_PREFIX_LIST :
            for band in filters:
                datakey = f"{prefix}_{band}"
                gpar.DATAKEY_LIST_RAW.append(datakey)

        # end load_HOSTKEY_band

    def store_varlist_obs(self):

        # for fakes, tack on true mag to list of variables per obs
        survey    = self.config_inputs['args'].survey

        varlist_obs   = gpar.VARNAMES_OBS_LIST    # var names
        varlist_fmt   = gpar.VARNAMES_FMT_LIST    # format per var
        vallist_undef = gpar.VAL_UNDEFINED_LIST   # value for undefined

        # check for unfilled PHOT columns to exclude in output
        for var in self.exclude_varlist_obs() :
            k          = varlist_obs.index(var)
            temp_obs   = varlist_obs.pop(k)
            temp_fmt   = varlist_fmt.pop(k)
            temp_undef = vallist_undef.pop(k)

        nvar    = len(varlist_obs)
        self.config_data['vallist_undef']  = vallist_undef
        self.config_data['varlist_fmt']    = varlist_fmt
        self.config_data['varlist_obs']    = varlist_obs
        self.config_data['nvar_obs']       = nvar

        #print(f" xxx nvar={nvar}  varlist = {varlist}")

        return

        # end store_varlist_obs

    def append_truemag_obs(self):

        # add true mag variable to list of variables in PHOT table
        # (called by read function of true mag exists; e.g. ,sim or fakes)

        logging.info(f"\t Append {gpar.VARNAME_TRUEMAG} to PHOT table")

        self.config_data['vallist_undef']  += [gpar.VAL_NULL]
        self.config_data['varlist_fmt']    += ["8.4f"]
        self.config_data['varlist_obs']    += [gpar.VARNAME_TRUEMAG]
        self.config_data['nvar_obs']       += 1

        return

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
        survey   = self.config_inputs['args'].survey
        d_raw    = data_event_dict['head_raw']
        d_calc   = data_event_dict['head_calc']

        snid     = d_raw[gpar.DATAKEY_SNID]
        zhel     = d_raw[gpar.DATAKEY_zHEL]
        zhel_err = d_raw[gpar.DATAKEY_zHEL_ERR]
        ra       = d_raw[gpar.DATAKEY_RA]
        dec      = d_raw[gpar.DATAKEY_DEC]

        snana_flag_fake = gpar.SNANA_FLAG_DATA

        # is SIM_MAGOBS columm exists, label data type as FAKE
        # Note that SIMs get labeled as FAKE.
        d_phot = data_event_dict['phot_raw']
        if gpar.VARNAME_TRUEMAG in d_phot:
            snana_flag_fake = gpar.SNANA_FLAG_FAKE

        if zhel > 0.0:
            zcmb      = util.helio_to_cmb(zhel, ra, dec)
        else:
            zcmb = gpar.VAL_NULL

        # no urgency for loading MWEBV because TEXT->FITS translator
        # computes and stores MWEBV. However, if we want correct MWEBV
        # in the TEXT files, need to compute it here:

        if gpar.DATAKEY_MWEBV in d_calc:
            mwebv     = d_calc[gpar.DATAKEY_MWEBV]
            mwebv_err = d_calc[gpar.DATAKEY_MWEBV_ERR]
        else:
            mwebv     = -9.0
            mwebv_err = -9.0

        # - - - -
        dump_flag = False
        if dump_flag:
            print(f" xxx ------------------------------")
            print(f" xxx DUMP for compute_data_event")
            print(f" xxx SNID={snid}   RA={ra}  DEC={dec}  zhel={zhel:8.5f}")
            if gpar.DATAKEY_zCMB in d_calc:
                zcmb_deja = d_calc[gpar.DATAKEY_zCMB]
                print(f"\t already existing zcmb = {zcmb_deja:8.5f}")
            if gpar.DATAKEY_MWEBV in d_calc:
                mwebv_deja     = d_calc[gpar.DATAKEY_MWEBV]
                mwebv_deja_err = d_calc[gpar.DATAKEY_MWEBV_ERR]
                print(f"\t already existing mwebv = " \
                      f"{mwebv_deja:8.5f} +_ {mwebv_deja_err:8.5f} ")

            print(f" xxx COMPUTE zcmb  = {zcmb:8.5f}")
            print(f" xxx COMPUTE mwebv = {mwebv:8.5f} +_ {mwebv:8.5f}")
            sys.stdout.flush()

        # - - - - - - -
        # load goodies
        d_raw[gpar.DATAKEY_SURVEY] = survey
        d_raw[gpar.DATAKEY_FAKE] = snana_flag_fake

        if survey not in gpar.SURVEY_INFO['FILTERS']:
            msgerr.append(f"{survey} filters not defined")
            msgerr.append(f"Check SURVEY_INFO dictionary in " \
                          f"makeDataFiles_params.py")
            util.log_assert(False, msgerr)
        else:
            d_raw[gpar.DATAKEY_FILTERS] = gpar.SURVEY_INFO['FILTERS'][survey]

        if survey in gpar.SURVEY_INFO['CCD']:
            d_raw[gpar.DATAKEY_NXPIX]   = gpar.SURVEY_INFO['CCD'][survey][0]
            d_raw[gpar.DATAKEY_NYPIX]   = gpar.SURVEY_INFO['CCD'][survey][1]
            d_raw[gpar.DATAKEY_PIXSIZE] = gpar.SURVEY_INFO['CCD'][survey][2]

        d_calc[gpar.DATAKEY_zCMB]      = zcmb
        d_calc[gpar.DATAKEY_zCMB_ERR]  = zhel_err
        d_calc[gpar.DATAKEY_MWEBV]     = mwebv
        d_calc[gpar.DATAKEY_MWEBV_ERR] = mwebv_err

        # if there is no VPEC, tack on default
        if gpar.DATAKEY_VPEC not in d_calc :
            d_calc[gpar.DATAKEY_VPEC]     = gpar.VPEC_DEFAULT[0]
            d_calc[gpar.DATAKEY_VPEC_ERR] = gpar.VPEC_DEFAULT[1]

        # check if there are spectra
        if 'spec_raw' in data_event_dict:
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

    def select_subsample(self, data_event_dict):

        # Apply optional subsample selection from command line
        # Return sel = True of False.

        args          = self.config_inputs['args']
        d_raw         = data_event_dict['head_raw']
        d_calc        = data_event_dict['head_calc']

        SNID_raw = d_raw[gpar.DATAKEY_SNID]
        # if SNID_raw.isdigit():
        #     SNID = int(SNID_raw)
        # else:
        #     SNID = SNID_raw

        SNID = int(SNID_raw) if SNID_raw.isdigit() else SNID_raw

        var_dict = {
            gpar.DATAKEY_SNID             : SNID,
            gpar.DATAKEY_RA               : d_raw[gpar.DATAKEY_RA],
            gpar.DATAKEY_DEC              : d_raw[gpar.DATAKEY_DEC],
            gpar.DATAKEY_PEAKMJD          : d_calc[gpar.DATAKEY_PEAKMJD],
            gpar.DATAKEY_MJD_DETECT_FIRST : d_calc[gpar.DATAKEY_MJD_DETECT_FIRST]
        }
        sel = util.select_subsample(args, var_dict)

        return sel

        # end select_subsample

    def init_phot_dict(self, NOBS):

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

        if t_dif_sec < 20000.0:
            t_dif   = t_dif_sec/60.0
            t_unit  = "minutes"
        else:
            t_dif  = t_dif_sec/3600.0
            t_unit = "hours"


        logging.info("\n FINAL MAKE-DATA-FILE SUMMARY: \n")
        logging.info(f" Total number of data units created: {NUNIT_TOT}")
        logging.info(f" Total number of events written:     {NEVT_TOT}")
        logging.info(f" CPUTIME(total):   {t_dif:.1f} {t_unit}" )
        # xxx mark logging.info(f" Total processing time ({t_unit}):    {t_dif:.1f}" )        
        sys.stdout.flush()

        # end final_summary

    def screen_update(self, evt, NEVT_TOT):

        if evt == 0:
            return

        if evt < 0:
            #time_0  = datetime.datetime.now()
            #self.config_data['time_0'] =  time_0
            self.config_data['time_0'] = datetime.datetime.now()
            return

        rmd = evt % gpar.NEVT_SCREEN_UPDATE
        if rmd == 0 or evt == NEVT_TOT-1:
            time_0   = self.config_data['time_0']
            time_now = datetime.datetime.now()
            time_dif = (time_now - time_0).total_seconds()
            rate     = int(float(evt)/float(time_dif))
            logging.info(f"\t\t Read evt={evt+1:6d} of {NEVT_TOT} "
                         f" ({rate}/sec)")
            sys.stdout.flush()

        return

    def read_data_driver(self):

        args = self.config_inputs['args']  # command line args

        # one-time init
        self.init_read_data()

        NEVT_READ    = 0

        nevent_subgroup = 1  # anything > 0 to pass while block below
        i_subgroup      = 0  # start with this subbroup index

        data_unit_name_list   = self.config_data['data_unit_name_list']

        if args.outdir_csv is not None:   # R.Kessler, July 2022
            csv_writer = csvWriter( args, self.config_data )

        while nevent_subgroup >= 0:

            nevent_subgroup = self.prep_read_data_subgroup(i_subgroup)

            if nevent_subgroup <  0:    break

            self.screen_update(-1,0)  # init clock for rate monitor

            for evt in range(0, nevent_subgroup):

                NEVT_READ += 1

                # call source-dependent function to read event
                data_event_dict = self.read_event(evt)

                # check optional subsample selection defined by reader;
                # if selection is not evaluated by reader, evaluate here.
                if 'select' in data_event_dict:
                    sel = data_event_dict['select']
                else:
                    # read_event did not select subsample, so do it here.
                    sel = self.select_subsample(data_event_dict)

                if sel is False:
                    continue

                # figure out which data unit
                data_unit_name = self.which_data_unit(data_event_dict)
                if data_unit_name is None:
                    continue

                # add computed variables; e.g., zCMB, MWEBV ...
                self.compute_data_event(data_event_dict)

                # add more info to data event dictionary
                index_unit  = data_unit_name_list.index(data_unit_name)
                data_event_dict['data_unit_name'] = data_unit_name
                data_event_dict['index_unit']     = index_unit

                if args.outdir_snana :
                    write_data_snana.write_event_text_snana(args, self.config_data,
                                                            data_event_dict)

                elif args.outdir_csv :
                    csv_writer.write_event_csv(data_event_dict)

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

        self.walltime_read_data()

        # write auxilary files for each data unit (name);
        # e.g.  LIST and README files.
        # Use aux_file util based on choice of output format.
        nevent_list   = self.config_data['data_unit_nevent_list']
        name_list     = self.config_data['data_unit_name_list']
        for nevent, name in zip(nevent_list, name_list):

            index_unit   = data_unit_name_list.index(name)


            # force creating yaml file even if there are zero events
            # to avoid merge problem later. This gets overwritten below
            # for nevent > 0.
            if index_unit == 0 and args.output_yaml_file:
                self.write_yaml_file(index_unit)
            
            if nevent == 0 : continue

            if args.output_yaml_file:
                self.write_yaml_file(index_unit)
                
            if args.outdir_snana:
                write_data_snana.write_aux_files_snana(name, args, self.config_data)

            elif args.outdir_csv :
                csv_writer.end_write(index_unit)

        # end read_data_driver

    def write_yaml_file(self, index_unit):

        # write yaml file to be parsed by pipeline.
        # This is the same file as README file in output directory.
        args         = self.config_inputs['args']
        readme_stats = self.config_data['readme_stats_list'][index_unit]
        t_proc       = self.config_data['t_proc'] # seconds


        readme_dict = {
            'readme_file'  : args.output_yaml_file,
            'readme_stats' : readme_stats,
            'data_format'  : gpar.FORMAT_TEXT,
            'docana_flag'  : False       # no DOCUMENTATION block
        }
        util.write_readme(args,readme_dict, t_proc )

        # end write_yaml_file

    def walltime_read_data(self):
        # compute wall time (t_proc) since t_start
        t_start = self.config_data['t_start']
        t_end   = datetime.datetime.now()
        t_dif_sec  = (t_end - t_start).total_seconds()
        self.config_data['t_end']      = t_end
        self.config_data['t_proc']     = t_dif_sec
        return
        # end walltime_read_data

    def update_readme_stats(self, data_event_dict):
        """Evaluate stats for readme of datafile. Has no impact on actual 
        datafiles. 
        """
        head_raw   = data_event_dict['head_raw']
        head_calc  = data_event_dict['head_calc']
        n_spectra  = data_event_dict['n_spectra']
        index_unit = data_event_dict['index_unit']

        # update stats that will eventually written to README file
        specz = -9.0
        photoz = -9.0

        if gpar.HOSTKEY_SPECZ  in head_raw:
            specz  = head_raw[gpar.HOSTKEY_SPECZ]
        if gpar.HOSTKEY_PHOTOZ in head_calc:
            photoz = head_calc[gpar.HOSTKEY_PHOTOZ]

        readme_stats = self.config_data['readme_stats_list'][index_unit]
        readme_sum   = self.config_data['readme_stats_sum']

        key = 'NEVT_ALL'
        readme_stats[key] += 1
        readme_sum[key]   += 1

        if specz > 0.0:
            key = 'NEVT_HOSTGAL_SPECZ'
            readme_stats[key] += 1
            readme_sum[key]   += 1

        if photoz > 0.0:
            key = 'NEVT_HOSTGAL_PHOTOZ'
            readme_stats[key] += 1
            readme_sum[key] += 1

        if n_spectra > 0:
            key = 'NEVT_SPECTRA'
            readme_stats[key] += 1  # increment this data unit
            readme_sum[key]   += 1  # increment this data unit
            self.config_data[key] += 1  # sum over all data units

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
    def iyear_data(self, MJD, RA, DEC, FIELD):
        raise NotImplementedError()

