import os,sys,glob,yaml,shutil
import logging, coloredlogs

import makeDataFiles_util  as    util
from   makeDataFiles_base import Program

# - - - - - - - - - - - - - - - - - - -     -
class data_lsst_ap(Program):
    def __init__(self, config_inputs, config_data) :
        config_data = {}
        print(" Init data_lsst_ap class.")
        super().__init__(config_inputs,config_data)
        

    def read_data_driver(self):
        pass
    
