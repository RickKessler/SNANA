# Created Oct 22, 2021
# write data in lsst-alert format for broker test.
# [R.Hlozek, R.Kessler ...]

import os, sys, yaml, shutil, glob, math, datetime
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
TIMEBACK_FORCE  = 50    #how many days before 1st detect to include forced phot.
NOBS_ALERT_MAX  = 2000  # used to compute diaSource

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

    print('')
    sys.stdout.flush()

    return schema, alert_data

    # end prep_write_lsst_alert
    
def write_event_lsst_alert(args, config_data, data_event_dict):

    # Inputs:
    #   args : user command line inputs
    #   config_data       : info about data units and phot varnames
    #   data_event_dict   : current event: header, phot, spec

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']
    SNID      = int(head_raw[DATAKEY_SNID]) # to compare sourceID
    NOBS      = phot_raw[DATAKEY_NOBS]

    # strip off number of processed events; init stuff on nevent=0
    data_unit_name        = data_event_dict['data_unit_name']
    index_unit            = data_event_dict['index_unit']
    data_unit_name_list   = config_data['data_unit_name_list']
    data_unit_nevent_list = config_data['data_unit_nevent_list']    
    nevent        = data_unit_nevent_list[index_unit]
    outdir        = args.outdir_lsst_alert

    if nevent == 0 :
        # later check for removing old folders ??
        schema, alert_data  = init_schema_lsst_alert(args.lsst_alert_schema)
        alert_data_orig     = alert_data.copy()
        config_data['schema']          = schema
        config_data['alert_data_orig'] = alert_data_orig
        config_data['diaSourceId']     = 0
        config_data['n_alert_write']   = 0
        config_data['n_event_write']   = 0
        config_data['t_start_alert']   = datetime.datetime.now()
        
    schema          = config_data['schema'] 
    diaSourceId     = config_data['diaSourceId'] 
    alert_data_orig = config_data['alert_data_orig']
    alert           = copy(alert_data_orig)
    
    # copy structure of original sample alert to local diasrc
    prvDiaSources = alert_data_orig['prvDiaSources']
    diasrc = prvDiaSources[0] 

    # # print this out for testing
    #for key in diasrc.keys():
    #    print(f" xxx found original alert key: {key}")
        
    alert['prvDiaSources'].clear() # clear out all the past histories
        
    # - - - - - -
    # translate snana header and create diasrc dictionary for lsst alert
    my_diasrc = diasrc #{}
    translate_dict_diasrc(-1, data_event_dict, my_diasrc)
    alert['diaSource'] = my_diasrc

    diaObjectId = my_diasrc['diaObjectId'] # same as SNID in snana sim file

    config_data['n_event_write'] += 1
    FIRST_OBS = True
    MJD_REF=data_event_dict['head_calc'][DATAKEY_PEAKMJD] # later change to  DATAKEY_MJD_DETECT
    
    #translate each obs to diasrc dictionary 
    for o in range(0,NOBS):
        mjd         = data_event_dict['phot_raw']['MJD'][o]
        keep_force = (MJD_REF - mjd) < TIMEBACK_FORCE and \
                     (mjd - MJD_REF) < 100  # temp until we have last MJD_DETECT
        if not keep_force: continue
        
        # skip non-detections (maybe later, add force photo after 1st detect?)
        photflag    = data_event_dict['phot_raw']['PHOTFLAG'][o]
        detect      = (photflag & PHOTFLAG_DETECT) > 0

        # compute UNIQUE diaSource from already unique SNID
        diaSourceId = NOBS_ALERT_MAX*SNID + o
        my_diasrc['diaSourceId'] = diaSourceId

        # update my_diasrc with this obs
        translate_dict_diasrc(o, data_event_dict, my_diasrc) 
    
        if FIRST_OBS:
            # Save my_diasrc info on 1st observation
            alert['diaSource'] = my_diasrc
            FIRST_OBS = False

        # serialize the alert    
        avro_bytes = schema.serialize(alert)
        messg      = schema.deserialize(avro_bytes)

        # write alerts ONLY for detection.
        # problem: first alert includes previous force photometry
        #   which violates causality.
        if detect :
            # construct name of avro file using mjd, objid, srcid
            outdir_mjd  = make_outdir_mjd(outdir,mjd)
            mjd_file  = f"{outdir_mjd}/" \
                        f"alert_mjd{mjd:.4f}_obj{diaObjectId}_src{diaSourceId}.avro"
            with open(mjd_file,"wb") as f:
                schema.store_alerts(f, [alert])
                config_data['n_alert_write'] += 1
                print_alert_stats(config_data)            
                
        # now that you have written out this alert,
        # move the diasource info to the "past" for the next observation
        alert['prvDiaSources'].append(alert['diaSource'])
        
        #print(f" xxx o={o}  mjd={mjd}")

    # - - - - - -           
    # if first obs was never found, hack alert to avoid crash.
    # This error can happen using MJD_REF=PEAKMJD, but should not 
    # occur when using MJD_REF = MJD_DETECT_FIRST.
    if  FIRST_OBS :
        print(f"   No detections found for SNID={SNID}; " \
              f"hack alert to avoid crash")
        sys.stdout.flush()

        o = 0
        translate_dict_diasrc(o, data_event_dict, my_diasrc)
        alert['diaSource'] = my_diasrc
        alert['prvDiaSources'].append(alert['diaSource'])

    return

# end write_event_lsst_alert

def make_outdir_mjd(outdir,mjd):

    # + construct name of mjd-specific outdir for this outdir and mjd
    # + create outdir_mjd if it does not already exit

    mjdint = int(mjd)
    outdir_mjd = outdir + '/mjd' + str(mjdint)
    if not os.path.exists(outdir_mjd) :
        cmd = f"mkdir {outdir_mjd}"
        #print(f"\t Create {outdir_mjd}")
        os.system(cmd)

    return outdir_mjd

# end make_outdir_mjd

def translate_dict_diasrc(obs, data_event_dict, diasrc):

    # obs = -1 -> translate header info in data_event_dict to diasrc
    # obs >= 0 -> translate obs in data_event_dict to diasrc

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']
    
    if obs < 0 :
        # xxx for key in diasrc.keys():    print(key)
        for varName_inp in VARNAME_HEADER_MAP:
            varName_avro = VARNAME_HEADER_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()

            if varName_inp in head_raw:
                if varName_inp == DATAKEY_SNID: # convert str to int
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

def print_alert_stats(config_data):
    
    n_alert = config_data['n_alert_write']
    n_event = config_data['n_event_write']
            
    if n_alert % 50 == 0 :
        t_start_alert = config_data['t_start_alert']
        t_now         = datetime.datetime.now()
        t_dif_sec  = (t_now - t_start_alert).total_seconds()
        rate       = int(n_alert / t_dif_sec)
        print(f"\t Wrote {n_alert:8d} alerts ({rate}/sec) " \
              f"for {n_event:6d} events.")
        sys.stdout.flush()

        #self.config_data['t_start'] = datetime.datetime.now()
        #t_start = self.config_data['t_start']

        
# end update_alert_stats
    
