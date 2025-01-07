# Created Jan 2025 by R.Kessler
# Read from F.A.S.T data base (for LSST-DESC) 

import os, sys, glob, yaml, shutil
import logging 
import numpy as np
import makeDataFiles_util  as    util
from   makeDataFiles_base    import Program
from   makeDataFiles_params  import *

class data_lsst_fastdb(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        print(" Init data_lsst_fastdb class.")
        super().__init__(config_inputs, config_data)

    def init_read_data(self):

        args = self.config_inputs['args']  # command line args

        
