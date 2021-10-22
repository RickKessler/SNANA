# Created Oct 22, 2021
# write data in lsst-alert format for broker test.
# [R.Hlozek, R.Kessler ...]

import os, sys, yaml, shutil, glob, math
import logging, coloredlogs, subprocess

import numpy as np
from   makeDataFiles_params    import *
import makeDataFiles_util  as util

#import lsst.alert.packet
#from pathlib import Path
#from fastavro import writer, reader


# map dictionary(SNANA) varName to alert varName 
lc = "lc"  # instruction to take lower case of dict value
VARNAME_HEADER_MAP = {
    DATAKEY_SNID            : 'diaObjectId',
    DATAKEY_RA              : lc,
    DATAKEY_DEC             : 'decl',
    DATAKEY_MWEBV           : lc,
    DATAKEY_MWEBV_ERR       : lc,
    DATAKEY_zHEL            : 'z_final' ,
    DATAKEY_zHEL_ERR        : 'z_final_err',
    DATAKEY_NOBS            : lc  # in phot_raw, not header
}

VARNAME_OBS_MAP = {
    'MJD'        : 'midPointTai',
    'BAND'       : 'filterName',
    'FLUXCAL'    : 'apFlux',
    'FLUXCALERR' : 'apFluxErr'
}

# ====================================================

def write_event_lsst_alert(args, config_data, data_event_dict):

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']
    SNID      = head_raw[DATAKEY_SNID] # for error message
    NOBS      = phot_raw[DATAKEY_NOBS]
    diasrc    = {}

    translate_dict_alert(-1, data_event_dict, diasrc) # translate header

    # translate each obs
    for o in range(0,NOBS):
        translate_dict_alert(o, data_event_dict, diasrc)
        #print(f"\t xxx write obs {o:3d} of {NOBS} for {SNID} ")


    sys.exit(f"\n xxx NOBS={NOBS}  diasrc = \n{diasrc}")
    # end write_event_lsst_alert


def translate_dict_alert(obs, data_event_dict, diasrc):
    # obs = -1 -> set header info
    # obs >= 0 -> set info for obs 

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']

    if obs < 0 :
        for varName_inp in VARNAME_HEADER_MAP:
            varName_avro = VARNAME_HEADER_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()

            if varName_inp in head_raw:
                diasrc[varName_avro] = head_raw[varName_inp]
            elif varName_inp in head_calc :
                diasrc[varName_avro] = head_calc[varName_inp]
            else:
                diasrc[varName_avro] = phot_raw[varName_inp]
                
            print(f" xxx {varName_inp} -> {varName_avro} ")
    else:
        for varName_inp in VARNAME_OBS_MAP:
            varName_avro = VARNAME_OBS_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()
            diasrc[varName_avro] = phot_raw[varName_inp][obs]
        
# end translate_dict_alert
