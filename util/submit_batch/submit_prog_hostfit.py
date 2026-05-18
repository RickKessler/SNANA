import os, sys, shutil, yaml, glob, logging
import datetime, time, subprocess
import submit_util as util
import pandas as pd

from   submit_params import *
from   submit_prog_base import Program


# Created May 2026
# run SED fit code on galaxies to get host properties.
#
class HostPropertyFit(Program):
    def __init__(self, config_yaml):

        config_prep = {}
        config_prep['program'] = None
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

        logging.info(f" xxx hello from new task HostPropertyFit: outdir = {output_dir_name} ")
        return output_dir_name, SUBDIR_SCRIPTS_HOSTFIT

        
    def submit_prepare_driver(self):

        pass
        return
