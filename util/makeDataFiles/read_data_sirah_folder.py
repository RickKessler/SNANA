import os, sys, glob, yaml, shutil, pickle
import logging  # , coloredlogs
import numpy as np
import makeDataFiles_util  as    util
from   makeDataFiles_base    import Program
from   makeDataFiles_params  import *

# J.Pierel 10/14 -- updated filter names
#filt_swap = {'F098M':'WFC3_IR_F098M-L','F105W':'WFC3_IR_F105W-Y',
#             'F110W':'WFC3_IR_F110W-M','F125W':'WFC3_IR_F125W-J',
#             'F140W':'WFC3_IR_F140W-N','F160W':'WFC3_IR_F160W-H',
#             'ZTFg':'ZTF-G','ZTFr':'ZTF-R','ZTFr':'ZTF-I'
#             }

FILTERMAP_PKL2SNANA = {
    # PKL ->     SNANA
    'ZTFg'     : 'ZTF-G' ,
    'ZTFr'     : 'ZTF-R' ,
    'ATLAS-c'  : 'ATLAS-c' ,
    'ATLAS-o'  : 'ATLAS-o' ,
    'PS1-g'    : 'PS1-g' ,
    'PS1-r'    : 'PS1-r' ,
    'PS1-i'    : 'PS1-i' ,
    'PS1-z'    : 'PS1-z' ,
    'PS1-y'    : 'PS1-y' ,
    'NIRIY'    : 'NIRI-Y/A', # Gemini Photometry
    'NIRIJ'    : 'NIRI-J/B', # Gemini Photometry
    'NIRIH'    : 'NIRI-H/C',	# Gemini Photometry
    'F098M'    : 'WFC3_IR_F098M-L' ,  # <lam> ~  9,900 A
    'F105W'    : 'WFC3_IR_F105W-Y' ,  # <lam> ~ 10,600 A
    'F125W'    : 'WFC3_IR_F125W-J' ,  # <lam> ~ 12,500 
    'F140W'    : 'WFC3_IR_F140W-N' ,  # <lam> ~ 14,000 A
    'F160W'    : 'WFC3_IR_F160W-H'    # <lam> ~ 15,400 A
}

# - - - - - - - - - - - - - - - - - - -     -
class data_sirah_folder(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        print(" Init data_sirah_folder class.")
        super().__init__(config_inputs, config_data)

    def init_read_data(self):

        args = self.config_inputs['args']  # command line args
        sirah_folder = args.sirah_folder
        
        # get list of pkl files in sirah folder
        wildcard = "*.pkl"
        pkl_file_list = glob.glob1(sirah_folder, wildcard )
        nevt = len(pkl_file_list)
        
        if nevt == 0 :
            msgerr = []
            msgerr.append(f"Could not find any {wildcard} files")
            msgerr.append(f"in {sirah_folder}")
            util.log_assert(False,msgerr)
            
        self.config_inputs['pkl_file_list'] = pkl_file_list
        self.config_inputs['nevt']          = nevt

        # end read_data_driver

    def prep_read_data_subgroup(self, i_subgroup):
        # There is only one subgroup for SIRAH
        nevt = self.config_inputs['nevt']
        if i_subgroup == 0:
            return nevt
        else:
            return 0
        # end prep_read_data_subgroup
        
    def end_read_data_subgroup(self):
        pass
    def end_read_data(self):
        # global end for reading data
        pass

    def	exclude_varlist_obs(self):
        # return list of PHOT columns to excude from output text files
        return [ 'XPIX', 'YPIX', 'CCDNUM', 'PHOTFLAG',
                 'GAIN', 'NEA', 'SKYSIG' ]
    
    def read_event(self, evt ):

        args = self.config_inputs['args']  # command line args 
        sirah_folder = args.sirah_folder
        pkl_file     = self.config_inputs['pkl_file_list'][evt]
        PKL_FILE     = f"{sirah_folder}/{pkl_file}"

        varlist_obs = self.config_data['varlist_obs']

        pkl_phot, pkl_spec = pickle.load(open(PKL_FILE,'rb'))
        SNID     = pkl_phot.meta['NAME']
        zHEL     = pkl_phot.meta['zHEL']
        zHEL_ERR = pkl_phot.meta['zHEL_ERR']            
        RA       = pkl_phot.meta['RA']
        DEC      = pkl_phot.meta['DEC']
        PEAKMJD  = pkl_phot.meta['fit_t0']
        MW_EBV   = pkl_phot.meta['MW_EBV']

        # init output dictionaries
        head_raw, head_calc = self.reset_data_event_dict()
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
        NOBS = len(pkl_phot)
        phot_raw = self.init_phot_dict(NOBS)

        phot_raw['NOBS']       = NOBS
        phot_raw['MJD']        = pkl_phot['MJD']
        phot_raw['ZPFLUX']     = pkl_phot['ZP']
        phot_raw['FIELD']      = [ FIELD_VOID ] * NOBS

        # use filter string map to translate pkl filter names
        # to names for SNANA.
        pkl_band_list = pkl_phot['Filter'] 
        phot_raw['BAND'] = \
            [FILTERMAP_PKL2SNANA[pkl_band] for pkl_band in pkl_band_list ]

        # convert SIRAH flux to FLUXCAL using ZP
        phot_Flux     = np.array(pkl_phot['Flux'].astype(np.float))
        phot_Fluxerr  = np.array(pkl_phot['Fluxerr'].astype(np.float))
        phot_ZP       = np.array(pkl_phot['ZP'])
        phot_scale    = np.power(10.0,(0.4*(SNANA_ZP-phot_ZP)))
        phot_raw['FLUXCAL']    = phot_Flux    * phot_scale
        phot_raw['FLUXCALERR'] = phot_Fluxerr * phot_scale
        
        #if evt == 0 :
        #    print(f"\n xxx NOBS = {NOBS} \n{pkl_phot.meta}")
                
        # copy spectra
        NSPEC    = len(pkl_spec) 
        spec_raw = pkl_spec

        # - - - - -
        # load output dictionary
        data_dict = {
            'head_raw'  : head_raw,
            'head_calc' : head_calc,
            'phot_raw'  : phot_raw,
            'spec_raw'  : spec_raw   # added Oct 14 2021
        }
        
        print(f"\t Read {pkl_file:<14}    " \
              f"NOBS={NOBS:3d}  NSPECTRA={NSPEC:2d}")

        return data_dict
  
        # end read_event
