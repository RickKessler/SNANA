# Read data from already existing folder in SNANA-FITS format.
# Intended only for testing makeDataFiles; not for production.
# Nov 30 2021: fix bug to account for PTROBS starting at 1 instead of 0.
# Dec 20 2021: fix dumb bug to read last HEAD & PHOT file.

import glob
import logging  # , coloredlogs
import os
import shutil
import sys

import yaml
from astropy.io import fits

import makeDataFiles_params as gpar
import makeDataFiles_util as util
from makeDataFiles_base import Program


#- - - - - - - - - - - - - - - - - -     -
#TODO turn this into a utility. This is more general behavior.
class data_snana_folder(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        print(" Init data_snana_folder class.")
        super().__init__(config_inputs,config_data)

        args          = self.config_inputs['args']  # command line args
        refac, legacy = self.snana_refac_legacy()

        if refac:
            # run __init__ in snana-reader class
            snana_folder = args.snana_folder
            SNANA_READER = util.READ_SNANA_FOLDER(snana_folder) 
            config_data['SNANA_READER'] = SNANA_READER

    def init_read_data(self):

        refac, legacy = self.snana_refac_legacy()

        if refac: 
            return

        if legacy:
            self.init_read_data_legacy()

        return

        # end init_read_data

    def init_read_data_legacy(self):
        args         = self.config_inputs['args']  # command line args
        data_folder  = os.path.expandvars(args.snana_folder)
        version      = os.path.basename(data_folder)
        list_file    = f"{data_folder}/{version}.LIST"
        logging.info(f" Read data version = {version}")
        logging.info(f" from data_foler = {data_folder}")
        # scoop up contents of LIST file
        with open(list_file, 'r') as f:
            HEAD_file_list = f.read().replace('\n', ' ').split()
        self.config_data['HEAD_file_list'] = HEAD_file_list
        self.config_data['n_HEAD_file']    = len(HEAD_file_list)
        self.config_data['version']        = version
        self.config_data['data_folder']    = data_folder
        return
        # end init_read_data_legacy


    def prep_read_data_subgroup(self, i_subgroup):

        refac, legacy = self.snana_refac_legacy()

        if refac :
            SNANA_READER = self.config_data['SNANA_READER']
            nevt = SNANA_READER.exec_read(i_subgroup)

        if legacy:
            nevt = self.prep_read_data_legacy(i_subgroup)

        # - - - -
        return nevt

        #end prep_read_data_subgroup

    def prep_read_data_legacy(self, i_subgroup):

        n_HEAD_file      = self.config_data['n_HEAD_file']
        if i_subgroup == n_HEAD_file  :  return 0 # done reading

        data_folder      = self.config_data['data_folder']
        HEAD_file_base   = self.config_data['HEAD_file_list'][i_subgroup]

        HEAD_file       = f"{data_folder}/{HEAD_file_base}"
        nevt, hdu_head  = util.open_fits(HEAD_file)

        PHOT_file_base  = hdu_head[0].header['PHOTFILE']
        PHOT_file       = f"{data_folder}/{PHOT_file_base}"
        NROW, hdu_phot  = util.open_fits(PHOT_file)

        logging.info(f"   Read {nevt} events from {HEAD_file_base}")

        table_head = hdu_head[1].data
        table_phot = hdu_phot[1].data
        head_names = table_head.columns.names
        phot_names = table_phot.columns.names

        # on first subgroup, check for true mag in PHOT table
        # e.g., fakes overlaid on images or sim
        if i_subgroup == 0 and gpar.VARNAME_TRUEMAG in phot_names:
            self.append_truemag_obs()

        table_dict = {
            'head_file'  : HEAD_file_base,
            'table_head' : table_head,
            'table_phot' : table_phot,
            'head_names' : head_names,
            'phot_names' : phot_names
        }

        self.config_data['nevt_subgroup'] = nevt
        self.config_data['table_dict'] = table_dict
        self.config_data['hdu_head']   = hdu_head
        self.config_data['hdu_phot']   = hdu_phot

        # - - - -
        return nevt

        #end prep_read_data_legacy

    def end_read_data_subgroup(self):
        
        refac, legacy = self.snana_refac_legacy()

        if refac:
            SNANA_READER = self.config_data['SNANA_READER']
            SNANA_READER.end_read()

        if legacy:
            self.config_data['hdu_head'].close()
            self.config_data['hdu_phot'].close()

        # end end_read_data_subgroup

    def end_read_data(self):
        # global end for reading data
        pass

    def read_event(self,evt):

        args          = self.config_inputs['args']  # command line args
        refac, legacy = self.snana_refac_legacy()

        if refac: 
            SNANA_READER = self.config_data['SNANA_READER']
            data_dict = SNANA_READER.get_data_dict(args,evt)

        if legacy:
            data_dict = self.read_event_legacy(evt)

        return data_dict

        # end read_event


    def read_event_legacy(self,evt):

        msgerr     = []
        table_dict = self.config_data['table_dict']
        args       = self.config_inputs['args']  # command line args

        # read and store one event for row "evt" and return data_dict.
        varlist_obs = self.config_data['varlist_obs']

        # define local pointers to head and phot tables from FITS file
        table_head = table_dict['table_head']
        table_phot = table_dict['table_phot']
        head_names = table_dict['head_names']
        phot_names = table_dict['phot_names']

        # init output dictionaries
        head_raw, head_calc, head_sim = util.reset_data_event_dict()

        try:
            SNID = table_head.SNID[evt].decode('utf-8').replace(' ','')
        except:
            SNID = table_head.SNID[evt]
        head_raw[gpar.DATAKEY_SNID]    = SNID

        head_raw[gpar.DATAKEY_SNTYPE]  = table_head.SNTYPE[evt]
        
        head_raw[gpar.DATAKEY_RA]    = table_head.RA[evt]

        # check 'DEC' and legacy column name 'DECL'
        head_raw[gpar.DATAKEY_DEC] = \
            util.get_snana_table_value(['DEC','DECL'],evt,table_head)

        # lightcurve-MJD info. Note that MJD_DETECT_FIRST is optional
        head_calc[gpar.DATAKEY_PEAKMJD]   = int(table_head.PEAKMJD[evt])

        if gpar.DATAKEY_MJD_DETECT_FIRST in head_names:
            head_calc[gpar.DATAKEY_MJD_DETECT_FIRST] = \
                table_head.MJD_DETECT_FIRST[evt]
            head_calc[gpar.DATAKEY_MJD_DETECT_LAST] = \
                table_head.MJD_DETECT_LAST[evt]
        else:
            if args.nite_detect_range is not None:
                msgerr.append(f"Cannot implement args.nite_detect_range = " \
                              f"{args.nite_detect_range}")
                msgerr.append(f"Because {gpar.DATAKEY_MJD_DETECT_FIRST} is not in "\
                              f"data header")
                util.log_assert(False,msgerr)

        # - - - - - - -
        # check user sub-sample selection here to avoid reading
        # remainder of header and photometry for rejected events.
        apply_select = True
        if apply_select :
            var_dict = {
                gpar.DATAKEY_SNID       : int(SNID),
                gpar.DATAKEY_RA         : head_raw[gpar.DATAKEY_RA],
                gpar.DATAKEY_DEC        : head_raw[gpar.DATAKEY_DEC],
                gpar.DATAKEY_PEAKMJD    : head_calc[gpar.DATAKEY_PEAKMJD],
                gpar.DATAKEY_MJD_DETECT_FIRST : head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]
            }
            sel = util.select_subsample(args,var_dict)
            if sel is False :
                data_dict = {
                    'head_raw'  : head_raw,
                    'head_calc' : head_calc,
                    'select'    : False
                }
                return data_dict

        # - - - - - - -
        head_raw[gpar.DATAKEY_zHEL]       = table_head.REDSHIFT_HELIO[evt]
        head_raw[gpar.DATAKEY_zHEL_ERR]   = table_head.REDSHIFT_HELIO_ERR[evt]

        # strip off calculated values
        head_calc[gpar.DATAKEY_zCMB]          = table_head.REDSHIFT_FINAL[evt]
        head_calc[gpar.DATAKEY_zCMB_ERR]      = table_head.REDSHIFT_FINAL_ERR[evt]
        head_calc[gpar.DATAKEY_MWEBV]         = table_head.MWEBV[evt]
        head_calc[gpar.DATAKEY_MWEBV_ERR]     = table_head.MWEBV_ERR[evt]

        # - - - - - -
        # store HOSTGAL and HOSTGAL2 keys in head_raw[calc]
        util.store_snana_hostgal(gpar.DATAKEY_LIST_RAW,  evt, table_dict, 
                                 head_raw )
        util.store_snana_hostgal(gpar.DATAKEY_LIST_CALC, evt, table_dict, 
                                 head_calc)

        # check for true sim type (sim or fakes), Nov 14 2021
        if gpar.SIMKEY_TYPE_INDEX in head_names:
            head_sim[gpar.SIMKEY_TYPE_INDEX] = table_head[gpar.SIMKEY_TYPE_INDEX][evt]

        # - - - - - - - - - - -
        # get pointers to PHOT table.
        # Beware that PTROBS pointers start at 1 instead of 0,
        # so subtract 1 here to have python indexing.
        ROWMIN = table_head.PTROBS_MIN[evt] - 1
        ROWMAX = table_head.PTROBS_MAX[evt] - 1

        NOBS     = ROWMAX - ROWMIN + 1
        phot_raw = self.init_phot_dict(NOBS)

        table_column_names = table_phot.columns.names

        if 'FLT' in table_column_names:
            LEGACY_FLT = True  # legacy column name is FLT for band
        else:
            LEGACY_FLT = False

        for varname in varlist_obs:
            varname_table = varname
            if LEGACY_FLT:
                if varname == 'BAND' : varname_table = 'FLT'

            if varname_table in table_column_names :
                phot_raw[varname] = \
                    table_phot[varname_table][ROWMIN:ROWMAX+1].copy()

        # - - - - -
        # get field from from first observation,
        # Beware that event can overlap multiple fields.
        field = phot_raw[gpar.DATAKEY_FIELD][0]
        if args.survey == 'LSST' :
            field = util.field_plasticc_hack(field,table_dict['head_file'])

        head_raw[gpar.DATAKEY_FIELD] = field

        # - - - -
        spec_raw = {}

        # - - - - -
        # load output dictionary
        data_dict = {
            'head_raw'  : head_raw,
            'head_calc' : head_calc,
            'phot_raw'  : phot_raw,
            'spec_raw'  : spec_raw,
        }
        if len(head_sim) > 0:
            data_dict['head_sim'] =  head_sim

        if apply_select :
            data_dict['select'] = True

        return data_dict

        # end read_event_legacy

    def set_dump_flag(self, isn, data_event_dict):
        d_raw = data_event_dict['head_raw']
        zhel  = d_raw['REDSHIFT_HELIO']
        dump_flag = False # isn<200 and zhel < 0.25
        return dump_flag
        # set_dump_flag

    def snana_refac_legacy(self):

        args   = self.config_inputs['args']  # command line args

        # default is legacy unless refac=110
        #refac  = args.refac == gpar.REFAC_READ_SNANA_FOLDER 
        #legacy = not refac

        # default is refac unless refac=-110
        legacy  = args.refac == gpar.LEGACY_READ_SNANA_FOLDER 
        refac   = not legacy

        return refac, legacy
        # end refac_legacy
