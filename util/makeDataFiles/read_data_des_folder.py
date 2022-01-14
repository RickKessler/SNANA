# Read DES(SMP) data from directory.

import glob
import logging
import os
import shutil
import sys
import tarfile

import yaml
import pandas as pd

from astropy.io import fits

import makeDataFiles_params as gpar
import makeDataFiles_util as util
from makeDataFiles_base import Program


# Hard coded constant until we have environment variable
# '/project2/rkessler/SURVEYS/DES/ROOT/smp_cori'
PATH_DES_SMP = os.getenv("DES_SMP")  

SMP_MASTERLIST_FILE ='masterlist.csv'

# Select private DES variables to keep
PRIVATE_VAR_DICT = {
    'DES_numepochs_ml_Y1' : 'Ndetect in Y1 passing autoscan',
    'DES_numepochs_ml_Y2' : 'Ndetect in Y2 passing autoscan',
    'DES_numepochs_ml_Y3' : 'Ndetect in Y3 passing autoscan',
    'DES_numepochs_ml_Y4' : 'Ndetect in Y4 passing autoscan',
    'DES_numepochs_ml_Y5' : 'Ndetect in Y5 passing autoscan',
    'AGN_SCAN'            : 'reject on value = 2'
}


# - - - - - - - - - - - - - - - - - - -     -
class data_des_folder(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        logging.info(" Init data_des_folder class.")
        super().__init__(config_inputs, config_data)

        args          = self.config_inputs['args']  # command line args

        des_folder    = args.des_folder
        SNANA_READER  = util.READ_SNANA_FOLDER(des_folder)
        self.config_data['SNANA_READER'] = SNANA_READER
        

    def init_read_data(self):

        # init which private variables to keep, along with
        # comment for TEXT file

        SNANA_READER = self.config_data['SNANA_READER']
        SNANA_READER.init_private_dict(PRIVATE_VAR_DICT)

        # .xyz initialize SMP; read master list ....
        
        # e.g., check $DES_SMP; else abort.
        #HARD CODE THE SMP location
        #self.PATH_SMP = PATH_SMP
        if not os.path.isdir(PATH_DES_SMP):
            msg_err = []
            msg_err.append(f"$DES_SMP path was not found")
            msg_err.append(f"Check {PATH_DES_SMP}")
            
            util.log_assert(False, msg_err)
        
        logging.info("Prepare DES SMP, setting up masterlist and tarballs")
        self.file_cache = {}
        self.file_cache['tarballs'] = {}

        self.masterlistpath = os.path.join(PATH_DES_SMP, SMP_MASTERLIST_FILE)
        self.smp_master_list = pd.read_csv(self.masterlistpath)

        self.n_smp_files = np.max(self.smp_master_list.tar_id)
        logging.info(f"Masterlist shows {self.n_smp_files} files")

        # end init_read_data

    def _get_phot_table(self, evt_id):

        if not str(evt_id) in self.file_cache.keys():
            phot_info = self.smp_master_list.query(f"cid=={evt_id}").iloc[0]

            self.file_cache[str(evt_id)] = dict(phot_info)
            tarball = self._get_tarball(phot_info.tar_id)

            lc_table_basename = os.path.join(
                os.path.dirname(phot_info.tar_path),
                f"csvfiles/phot_{evt_id}.csv"
            )
            lctab = pd.read_csv(tarball.extractfile(lc_table_basename))
            self.file_cache[str(evt_id)]['lctab'] = lctab

            return lctab

        return self.file_cache[str(evt_id)]['lctab']


    def _get_tarball(self, tar_id):
        
        if not str(tar_id) in self.file_cache['tarballs'].keys():
            phot_table_basename = f"phot_{str(tar_id).zfill(2)}.tar.gz"
            if tarfile.is_tarfile(phot_table_basename):
                tarball_file = tarfile.open(phot_table_basename, "r:gz")
            else:
                raise IOError(f"File {phot_table_basename} not found")

            self.file_cache['tarballs'][str(tar_id)] = tarball_file
        
        return self.file_cache['tarballs'][str(tar_id)]


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
        """Teardown method for closing open files and clearing cache."""

        # close open tarballs
        for tar_id, afile in self.file_cache['tarballs'].iteritems():
            afile.close()

        return

    def read_event(self, evt):

        args         = self.config_inputs['args']  # command line args
        SNANA_READER = self.config_data['SNANA_READER']
        data_dict    = SNANA_READER.get_data_dict(args, evt)

        # MJD_trigger is a private variable, so move it to 
        # nominal SNANA variable. .xyz
        key_private_list = [ 'PRIVATE(DES_mjd_trigger)' ]
        key_head_list    = [ gpar.DATAKEY_MJD_DETECT_FIRST ]
        head_calc        = data_dict['head_calc']
        for key_private, key_head in zip(key_private_list, key_head_list):
            head_calc[key_head] = SNANA_READER.get_data_val(key_private, evt)

        # .xyz read SMP for this SNID and overwrite FLUXCAL[ERR]
        # tricky part: DIFFIMG has ~500 epochs spanning 5 years,
        # but SMP has ~100 epochs spanning 1 season. So need to
        # copy epoch meta-data (SKY,GAIN,ZP,PSF) from original
        # diffimg phot array to final smp array,

        smp_lc = self._get_phot_table(evt)
        

        return data_dict

        # end read_event


    def set_dump_flag(self, isn, data_event_dict):
        d_raw = data_event_dict['head_raw']
        zhel  = d_raw['REDSHIFT_HELIO']
        dump_flag = False # isn<200 and zhel < 0.25
        return dump_flag
        # set_dump_flag


