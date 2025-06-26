# Created June 2025 by R.Kessler and C. Meldorf
# Initial use is MultiSMP for Roman.
#

import os, sys, logging, datetime, time, subprocess
import submit_util as util
import pandas as pd

from   submit_params    import *
from   submit_prog_base import Program

# ===================================================


# not sure we need this subclass dictionary, but maybe later if there are multiple
# SMP implementations
KEY_SUBCLASS_DICT = {
    'MultiSMP'  : [] 
}


# define columns in merge file
COLNUM_FIT_MERGE_STATE           = 0  # STATE required in col=0
COLNUM_FIT_MERGE_VERSION         = 1  # ??? naively copied from LCFIT class
COLNUM_FIT_MERGE_FITOPT          = 2  # ??? idem
COLNUM_FIT_MERGE_NCCD            = 3
COLNUM_FIT_MERGE_NLC             = 4
COLNUM_FIT_MERGE_CPU             = 5

# ====================================================
#    BEGIN CLASS
# ====================================================


class SceneModelPhotometry(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = PROGRAM_NAME_SMP
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
        return
        
    def submit_prepare_driver(self):
        
        # ...

        return
        # end submit_prepare_driver

