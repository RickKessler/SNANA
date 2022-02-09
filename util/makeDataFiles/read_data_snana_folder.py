# Read data from already existing folder in SNANA-FITS format.
# Intended only for testing makeDataFiles; not for production.
# Nov 30 2021: fix bug to account for PTROBS starting at 1 instead of 0.
# Dec 20 2021: fix dumb bug to read last HEAD & PHOT file.
# Feb 03 2022: integrate refac code; remove legacy (snana-reader util)
#
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
        super().__init__(config_inputs, config_data)

        args          = self.config_inputs['args']  # command line args
        #refac, legacy = self.snana_refac_legacy()

        # run __init__ in snana-reader class
        snana_folder = args.snana_folder
        SNANA_READER = util.READ_SNANA_FOLDER(snana_folder)
        config_data['SNANA_READER'] = SNANA_READER

    def init_read_data(self):
        return
        # end init_read_data

    def prep_read_data_subgroup(self, i_subgroup):
        SNANA_READER = self.config_data['SNANA_READER']
        nevt = SNANA_READER.exec_read(i_subgroup)
        # - - - -
        return nevt

        #end prep_read_data_subgroup

    def end_read_data_subgroup(self):

        SNANA_READER = self.config_data['SNANA_READER']
        SNANA_READER.end_read()

        # end end_read_data_subgroup

    def end_read_data(self):
        # global end for reading data
        pass

    def read_event(self, evt):
        args          = self.config_inputs['args']  # command line args
        SNANA_READER = self.config_data['SNANA_READER']
        data_dict = SNANA_READER.get_data_dict(args,evt)
        return data_dict
        # end read_event


    def set_dump_flag(self, isn, data_event_dict):
        d_raw = data_event_dict['head_raw']
        zhel  = d_raw['REDSHIFT_HELIO']
        dump_flag = False # isn<200 and zhel < 0.25
        return dump_flag
        # set_dump_flag

    def snana_refac_legacy(self):
        args   = self.config_inputs['args']  # command line args
        legacy  = true
        refac   = not legacy
        return refac, legacy
        # end refac_legacy
