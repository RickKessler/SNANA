# Created Oct 22, 2021
# write data in lsst-alert format for broker test.
# [R.Hlozek, R.Kessler ...]

import os, sys, yaml, shutil, glob, math
import logging, subprocess, json

import numpy as np
from   makeDataFiles_params    import *
import makeDataFiles_util  as util

from pathlib import Path
import lsst.alert.packet
from fastavro import writer, reader
from copy import copy


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
    DATAKEY_NOBS            : lc,       # in phot_raw, not header
    #
    HOSTKEY_SNSEP           : lc,
    HOSTKEY_SPECZ           : 'hostgal_z' ,
    HOSTKEY_SPECZ_ERR       : 'hostgal_z_err'
}

#HOSTKEY_OBJID         = "HOSTGAL_OBJID"
#HOSTKEY_PHOTOZ        = "HOSTGAL_PHOTOZ"
#HOSTKEY_PHOTOZ_ERR    = "HOSTGAL_PHOTOZ_ERR"
#HOSTKEY_LOGMASS       = "HOSTGAL_LOGMASS"

for prefix in [ 'HOSTGAL_MAG', 'HOSTGAL_MAGERR' ] :
    for band in list(SURVEY_INFO['FILTERS']['LSST']):
        key = f"{prefix}_{band}"
        VARNAME_HEADER_MAP[key] = lc

VARNAME_OBS_MAP = {
    'MJD'        : 'midPointTai',
    'BAND'       : 'filterName',
    'FLUXCAL'    : 'apFlux',
    'FLUXCALERR' : 'apFluxErr'
}

PHOTFLAG_DETECT = 4096  # should read this from global data header ??


# ===============================================================
def init_schema_lsst_alert(schema_file):

    schema     = lsst.alert.packet.Schema.from_file(filename=schema_file)
    schema_dir = os.path.dirname(schema_file)
    json_file  = f"{schema_dir}/sample_data/plasticc.json"  # too much hard coding

    print(f"\n Init alert schema based on\n\t schema_file={schema_file}\n" \
          f"\t jon_file={json_file}")
    
    # Load an example json alert, and clear the numberical input
    with open(json_file) as f:
        alert_data = json.load(f)
    
    return schema, alert_data

    # end prep_write_lsst_alert
    
def write_event_lsst_alert(args, config_data, data_event_dict):

    # Inputs:
    #   args : user command line inputs
    #   config_data       : info about data units and phot varnames
    #   data_event_dict   : current event: header, phot, spec

    head_raw = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']
    SNID      = head_raw[DATAKEY_SNID] # for error message
    NOBS      = phot_raw[DATAKEY_NOBS]

    # strip off number of processed events; init stuff on nevent=0
    data_unit_name    = data_event_dict['data_unit_name']
    index_unit        = data_event_dict['index_unit']
    data_unit_name_list   = config_data['data_unit_name_list']
    data_unit_nevent_list = config_data['data_unit_nevent_list']    
    nevent            = data_unit_nevent_list[index_unit]
    outdir            = args.outdir_lsst_alert

    if nevent == 0 :
        # later check for removing old folders ??
        schema, alert_data  = init_schema_lsst_alert(args.lsst_alert_schema)
        alert_data_orig = alert_data.copy()
        config_data['schema']          = schema
        config_data['alert_data_orig'] = alert_data_orig
        config_data['diaSourceId']     = 1000000
        
    schema          = config_data['schema'] 
    diaSourceId     = config_data['diaSourceId'] 
    alert_data_orig = config_data['alert_data_orig']
    alert           = copy(alert_data_orig)
    
    # copy structure of original sample alert to local diasrc
    prvDiaSources = alert_data_orig['prvDiaSources']
    diasrc = prvDiaSources[0] 

    # # print this out for testing
    #for key in diasrc.keys():
    #    print(key)
    #print('---------- diasrc keys ----------')
        
    alert['prvDiaSources'].clear() # clear out all the past histories

    # - - - - - -
    # translate snana header and create diasrc dictionary for lsst alert
    my_diasrc = diasrc #{}
    translate_dict_diasrc(-1, data_event_dict, my_diasrc)
    alert['diaSource'] = my_diasrc

    diaObjectId = my_diasrc['diaObjectId'] # same as SNID in snana sim file

    # translate each obs to diasrc dictionary 
    for o in range(0,NOBS)

        diaSourceId += 1
        my_diasrc['diaSourceId'] = diaSourceId

        translate_dict_diasrc(o, data_event_dict, my_diasrc) # update my_diasrc
    
        if o == 0 :
            # Save the my_diasrc info input the diaSource to be saved
            alert['diaSource'] = my_diasrc
            continue

        # check for detection
        photflag    = data_event_dict['phot_raw']['PHOTFLAG'][o]
        detect      = (photflag & PHOTFLAG_DETECT) > 0

        # construct name of avro file using mjd, objid, srcid
        outdir_mjd  = outdir  # xxx later tack on /[mjdint]
        mjd         = data_event_dict['phot_raw']['MJD'][o]
        mjd_file    = f"{outdir_mjd}/{mjd}_{diaObjectId}_{diaSourceId}.avro"

        # serialize the alert    
        avro_bytes = schema.serialize(alert)
        messg      = schema.deserialize(avro_bytes)
        
        with open(mjd_file,"wb") as f:
            schema.store_alerts(f, [alert])

        # now that you have written out this alert,
        # move the diasource info to the "past" for the next observation
        alert['prvDiaSources'].append(alert['diaSource'])
        
        #print(f" xxx o={o}  mjd={mjd}")
        
    # end write_event_lsst_alert

def translate_dict_diasrc(obs, data_event_dict, diasrc):

    # obs = -1 -> translate header info in data_event_dict to diasrc
    # obs >= 0 -> translate obs in data_event_dict to diasrc

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']
    
    if obs < 0 :
        for key in diasrc.keys():
            print(key)
        for varName_inp in VARNAME_HEADER_MAP:
            varName_avro = VARNAME_HEADER_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()

            if varName_inp in head_raw:
                if varName_inp == DATAKEY_SNID:
                    diasrc[varName_avro] = int(head_raw[varName_inp])
                else:
                    diasrc[varName_avro] = head_raw[varName_inp]
                    
            elif varName_inp in head_calc :
                diasrc[varName_avro] = head_calc[varName_inp]
            else:
                diasrc[varName_avro] = phot_raw[varName_inp]
                
            #print(f" xxx {varName_inp} -> {varName_avro} ")
    else:
        for varName_inp in VARNAME_OBS_MAP:
            varName_avro = VARNAME_OBS_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()
            diasrc[varName_avro] = phot_raw[varName_inp][obs]
        
    # end translate_dict_diasrc
