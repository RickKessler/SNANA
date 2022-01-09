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

#from makeDataFiles_params import *

# - - - - - - - - - - - - - - - - - - -     -
class data_des_folder(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        print(" Init data_des_folder class.")
        super().__init__(config_inputs,config_data)

        args          = self.config_inputs['args']  # command line args

        des_folder    = args.des_folder
        SNANA_READER  = util.READ_SNANA_FOLDER(des_folder) 
        config_data['SNANA_READER'] = SNANA_READER

    def init_read_data(self):
        pass
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

    def open_snana_fits(self,file_name):
       # check file_name and file_name.gz, and open the file that exists.
       # Function returns hdu pointer and number of rows in table.

        msgerr = []
        file_namegz = f"{file_name}.gz"
        if os.path.exists(file_namegz) :
            hdul = fits.open(file_namegz)
        elif os.path.exists(file_name):
            hdul = fits.open(file_name)
        else:
            msgerr.append(f"Cannot find fits file")
            msgerr.append(f" {file_name}   not")
            msgerr.append(f" {file_namegz} ")
            util.log_assert(False,msgerr)

        NROW = hdul[1].header['NAXIS2']
        return NROW, hdul

    # end open_snana_fits

    def read_event(self, evt ):

        args         = self.config_inputs['args']  # command line args
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


