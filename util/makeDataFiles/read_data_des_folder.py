# Read DES(SMP) data from directory.

import glob
import logging
import os
import shutil
import sys

import yaml
from astropy.io import fits

import makeDataFiles_params as gpar
import makeDataFiles_util as util
from makeDataFiles_base import Program

# - - - - - - - - - - - - - - - - - - -     -
class DataFolderDES(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        logging.info(" Init DataFolderDES class.")
        super().__init__(config_inputs, config_data)

        args          = self.config_inputs['args']  # command line args

        des_folder    = args.des_folder
        SNANA_READER  = util.SNANA_FolderReader(des_folder)
        self.config_data['SNANA_READER'] = SNANA_READER

    def init_read_data(self):

        # init which private variables to keep, along with
        # comment for TEXT file
        private_dict = {
            'DES_numepochs_ml_Y1' : 'Ndetect in Y1 passing autoscan',
            'DES_numepochs_ml_Y2' : 'Ndetect in Y2 passing autoscan',
            'DES_numepochs_ml_Y3' : 'Ndetect in Y3 passing autoscan',
            'DES_numepochs_ml_Y4' : 'Ndetect in Y4 passing autoscan',
            'DES_numepochs_ml_Y5' : 'Ndetect in Y5 passing autoscan',
            'AGN_SCAN'            : 'reject on value = 2'
        }

        SNANA_READER = self.config_data['SNANA_READER']
        SNANA_READER.init_private_dict(private_dict)

        # .xyz initialize SMP; read master list ....
        # e.g., check $DES_SMP; else abort.
        # self.config_data['smp_master_list'] = master_list

        # end init_read_data

    def prep_read_data_subgroup(self, i_subgroup):

        SNANA_READER = self.config_data['SNANA_READER']
        nevt         = SNANA_READER.exec_read(i_subgroup)
        return nevt

        #end prep_read_data_subgroup

    def end_read_data_subgroup(self):
        SNANA_READER = self.config_data['SNANA_READER']
        SNANA_READER.end_read()
        # end end_read_data_subgroup

    def end_read_data(self):
        # global end for reading data
        pass

    def read_event(self, evt ):

        args         = self.config_inputs['args']  # command line args
        SNANA_READER = self.config_data['SNANA_READER']
        data_dict    = SNANA_READER.get_data_dict(args,evt)

        # MJD_trigger is a private variable, so move it to 
        # nominal SNANA variable. .xyz
        key_private_list = [ 'PRIVATE(DES_mjd_trigger)' ]
        key_head_list    = [ gpar.DATAKEY_MJD_DETECT_FIRST ]
        head_calc        = data_dict['head_calc']
        for key_private, key_head in zip(key_private_list,key_head_list):
            head_calc[key_head] = SNANA_READER.get_data_val(key_private,evt)

        # .xyz read SMP for this SNID and overwrite FLUXCAL[ERR]
        # tricky part: DIFFIMG has ~500 epochs spanning 5 years,
        # but SMP has ~100 epochs spanning 1 season. So need to
        # copy epoch meta-data (SKY,GAIN,ZP,PSF) from original
        # diffimg phot array to final smp array,


        return data_dict

        # end read_event

    def set_dump_flag(self, isn, data_event_dict):
        d_raw = data_event_dict['head_raw']
        zhel  = d_raw['REDSHIFT_HELIO']
        dump_flag = False # isn<200 and zhel < 0.25
        return dump_flag
        # set_dump_flag


