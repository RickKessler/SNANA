import os, sys, glob, yaml, shutil, pickle
import logging  # , coloredlogs
import numpy as np
import makeDataFiles_util  as    util
from   makeDataFiles_base    import Program
from   makeDataFiles_params  import *

# FILTERMAP_PKL2SNANA = {
#     # CSV ->     SNANA
#     'p48g'     : 'ZTF-g' ,
#     'p48r'     : 'ZTF-r' ,
#     'p48i'     : 'ZTF-i' ,
# }

# - - - - - - - - - - - - - - - - - - -     -
class data_lsst_tom_db(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        print(" Init lsst data tom class.")
        super().__init__(config_inputs, config_data)

    def init_read_data(self):

        args = self.config_inputs['args']  # command line args
        db_connect_info = args.lsst_tom_db

        msgerr = [ "LSST TOM db reading code isn't written yet, no connection actually made" ]
        util.log_assert(False,msgerr)

        # Figure out number of events
        nevt = 1
        self.config_inputs['nevt']      = nevt

        # end read_data_driver

    def prep_read_data_subgroup(self, i_subgroup):
        # No subgroups defined yet
        nevt = self.config_inputs['nevt']
        if i_subgroup == 0:
            return nevt
        else:
            return -1
        # end prep_read_data_subgroup
        
    def end_read_data_subgroup(self):
        pass
    def end_read_data(self):
        # global end for reading data
        pass

    def	exclude_varlist_obs(self):
        # return list of PHOT columns to excude from output text files
        return []
    
    def read_event(self, evt ):

        varlist_obs = self.config_data['varlist_obs']

        # xxx pkl_phot, plk_spec = pickle.load(open(PKL_FILE,'rb'))
        # replace pkl read with astropy read ?
        
        SNID     = 'abc'   # pkl_phot.meta['NAME']
        zHEL     = 0.0132  # pkl_phot.meta['zHEL']
        zHEL_ERR = 0.001   # pkl_phot.meta['zHEL_ERR']            
        RA       = 24.343453  # pkl_phot.meta['RA']
        DEC      = 14.32234   # pkl_phot.meta['DEC']
        PEAKMJD  = 56908.32   # pkl_phot.meta['fit_t0']
        MW_EBV   = 0.023      # pkl_phot.meta['MW_EBV']

        # init output dictionaries
        head_raw, head_calc, head_sim = util.reset_data_event_dict()
        head_raw[DATAKEY_SNID]        = SNID
        head_raw[DATAKEY_RA]          = RA
        head_raw[DATAKEY_DEC]         = DEC
        head_raw[DATAKEY_zHEL]        = zHEL
        head_raw[DATAKEY_zHEL_ERR]    = zHEL_ERR
        head_raw[DATAKEY_FIELD]       = FIELD_VOID

        # calc quantities
        head_calc[DATAKEY_PEAKMJD]     = int(PEAKMJD)
        head_calc[DATAKEY_MWEBV]       = MW_EBV
        head_calc[DATAKEY_MWEBV_ERR]   = MW_EBV*0.16

        # copy photometry ...
        NOBS = 22  # len(pkl_phot)
        phot_raw = self.init_phot_dict(NOBS)

        phot_raw['NOBS']       = NOBS
        #phot_raw['MJD']        = pkl_phot['MJD']
        #phot_raw['ZPFLUX']     = pkl_phot['ZP']
        #phot_raw['FIELD']      = [ FIELD_VOID ] * NOBS

        # use filter string map to translate pkl filter names
        # to names for SNANA.
        #pkl_band_list = pkl_phot['Filter'] 
        #phot_raw['BAND'] = \
        #    [FILTERMAP_PKL2SNANA[pkl_band] for pkl_band in pkl_band_list ]

        # convert ZTF flux to FLUXCAL using ZP
        #phot_Flux     = np.array(pkl_phot['Flux'].astype(np.float))
        #phot_Fluxerr  = np.array(pkl_phot['Fluxerr'].astype(np.float))
        #phot_ZP       = np.array(pkl_phot['ZP'])
        #phot_scale    = np.power(10.0,(0.4*(SNANA_ZP-phot_ZP)))
        #phot_raw['FLUXCAL']    = phot_Flux    * phot_scale
        #phot_raw['FLUXCALERR'] = phot_Fluxerr * phot_scale
        
        #if evt == 0 :
        #    print(f"\n xxx NOBS = {NOBS} \n{pkl_phot.meta}")
                
        # - - - - -
        # load output dictionary
        data_dict = {
            'head_raw'  : head_raw,   # measured by ztf
            'head_calc' : head_calc,  # computed 
            'phot_raw'  : phot_raw   
        }
        
        return data_dict
  
        # end read_event
