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

import gzip

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

LSST_ZP_nJy     = 31.4  # report calibrated flux in this unit

# FLXUCAL(ALERT) = FLUXCAL(SNANA)*SCALE_FLUXCAL
ARG_ZPDIF       = 0.4*(LSST_ZP_nJy-SNANA_ZP)
SCALE_FLUXCAL   = math.pow(10.0,ARG_ZPDIF) 
KEYNAME_SUBSTRING_FLUXCAL = 'FLUXCAL' # scale variables with this substring

PHOTFLAG_DETECT = 4096  # should read this from global data header ??

NOBS_ALERT_MAX  = 2000  # used to compute diaSource
NOBS_ALERT_UPDATE = 100 # std update after this many alerts

ALERT_DAY_NAME    = "NITE"
LSST_SITE_NAME    = "CTIO"

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

    # end init_schema_lsst_alert

def init_truth_dict(outfile):
    truth_dict = {
        'outfile'      : outfile,
        'nobs_store'   : 0,
        'diaSourceId'  : [],
        'snid'         : [],
        'mjd'          : [],
        'detect'       : [],
        'true_gentype' : [],
        'true_genmag'  : []
    }
    return truth_dict
    # end init_truth_dict
    
def write_event_lsst_alert(args, config_data, data_event_dict):

    # Inputs:
    #   args : user command line inputs
    #   config_data       : info about data units and phot varnames
    #   data_event_dict   : current event: header, phot, spec

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']

    # scale SNANA-sim flux[Err] to alert flux unit (nJy)
    for key in VARNAME_OBS_MAP:
        if KEYNAME_SUBSTRING_FLUXCAL in key:
            phot_raw[key] = [f*SCALE_FLUXCAL for f in phot_raw[key] ]

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
        schema, alert_data  = init_schema_lsst_alert(args.lsst_alert_schema)
        alert_data_orig     = alert_data.copy()
        config_data['schema']          = schema
        config_data['alert_data_orig'] = alert_data_orig
        config_data['diaSourceId']     = 0
        config_data['n_alert_write']   = 0
        config_data['n_event_write']   = 0
        config_data['t_start_alert']   = datetime.datetime.now()

        outfile = args.outfile_alert_truth
        if outfile :
            config_data['truth_dict'] = init_truth_dict(outfile)
        else:
            config_data['truth_dict'] = None
            
    # - - - - - 
    schema             = config_data['schema']
    diaSourceId        = config_data['diaSourceId']
    alert_data_orig    = config_data['alert_data_orig']
    alert              = copy(alert_data_orig)
    alert_first_detect = copy(alert_data_orig)

    # copy structure of original sample alert to local diasrc
    prvDiaSources = alert_data_orig['prvDiaSources']
    diasrc = prvDiaSources[0]

    # # print this out for testing
    #for key in diasrc.keys():
    #    print(f" xxx found original alert key: {key}")

    alert['prvDiaSources'].clear() # clear out all the past histories
    alert_first_detect['prvDiaSources'].clear()

    # - - - - - -
    # translate snana header and create diasrc dictionary for lsst alert
    my_diasrc = diasrc #{}
    translate_dict_diasrc(-1, data_event_dict, my_diasrc)
    true_gentype = data_event_dict['head_sim'][SIMKEY_TYPE_INDEX]
    
    alert['diaSource']              = my_diasrc
    alert_first_detect['diaSource'] = my_diasrc

    diaObjectId = my_diasrc['diaObjectId'] # same as SNID in snana sim file

    config_data['n_event_write'] += 1
    FIRST_OBS    = True
    nobs_detect  = 0
    nobs_keep    = 0 # nobs_detect + nobs_forcePhoto
    
    if args.nite_detect_range:
        MJD_REF  = head_calc[DATAKEY_MJD_DETECT_FIRST]
        MJD_LAST = head_calc[DATAKEY_MJD_DETECT_LAST]
        TIME_BACK_FORCE  = 30   #Ndays before MJD_REF to include forced phot.
        TIME_FORWARD_FORCE = MJD_LAST - MJD_REF + 0.1
    elif args.peakmjd_range:
        MJD_REF = head_calc[DATAKEY_PEAKMJD]
        TIME_BACK_FORCE  = 50   #Ndays before MJD_REF to include forced phot.
        TIME_FORWARD_FORCE = 100

    # compute things for each obs
    mjd_list         = data_event_dict['phot_raw']['MJD']
    photflag_list    = data_event_dict['phot_raw']['PHOTFLAG']
    true_genmag_list = data_event_dict['phot_raw'][VARNAME_TRUEMAG]

    # get boolean list of obs within TIME_BACK_FORCE from 1st detect
    # and up to last detection
    keep_force_list = [ (MJD_REF-mjd)<TIME_BACK_FORCE and \
                        (mjd-MJD_REF)<TIME_FORWARD_FORCE \
                        for mjd in mjd_list ]

    set_limits_o = True # False  # testing for optimal speed for obs loop
    if set_limits_o :
        # find min and max obs index where keep_force = True
        reversed_list = keep_force_list[::-1]
        o_start = keep_force_list.index(True)
        o_end   = NOBS - reversed_list.index(True)
    else:
        # brute-force loop over all obs
        o_start = 0
        o_end   = NOBS
        
    # get list of boolean detection flag per observation
    detect_list     = [ (photflag & PHOTFLAG_DETECT)>0 \
                        for photflag in photflag_list ]

    # - - - - - - - - - - -
    #translate each obs to diasrc dictionary
    for o in range(o_start,o_end):
        keep_force = keep_force_list[o]
        if not keep_force: continue
        
        nobs_keep += 1
        mjd       = data_event_dict['phot_raw']['MJD'][o]
        detect    = detect_list[o]
        
        # compute UNIQUE diaSource from already unique SNID
        diaSourceId = NOBS_ALERT_MAX*SNID + o
        my_diasrc['diaSourceId'] = diaSourceId

        # update my_diasrc with this obs
        translate_dict_diasrc(o, data_event_dict, my_diasrc)

        if FIRST_OBS  :
            # Save my_diasrc info on 1st observation
            alert['diaSource'] = my_diasrc
            FIRST_OBS = False

        # serialize the alert
        avro_bytes = schema.serialize(alert)
        messg      = schema.deserialize(avro_bytes)

        # write truth table for detections and forced photo
        if args.outfile_alert_truth:
            truth_dict = config_data['truth_dict']
            #print(f"\n xxx truth_dict = {truth_dict} \n")
            truth_dict['nobs_store'] += 1
            truth_dict['diaSourceId'].append(diaSourceId)
            truth_dict['snid'].append(SNID)
            truth_dict['mjd'].append(mjd)
            truth_dict['detect'].append(detect)
            truth_dict['true_gentype'].append(true_gentype)
            truth_dict['true_genmag'].append(true_genmag_list[o])
        
        # write alerts ONLY for detection.
        # problem: first alert includes previous force photometry
        #   which violates causality.
        if detect :
            nobs_detect += 1

            # construct name of avro file using mjd, objid, srcid
            outdir_nite  = make_outdir_nite(outdir,mjd)
            str_day = f"mjd{mjd:.4f}"
            str_obj = f"obj{diaObjectId}"
            str_src = f"src{diaSourceId}"
            mjd_file  = f"{outdir_nite}/" \
                        f"alert_{str_day}_{str_obj}_{str_src}.avro"

            gzip_mjd_file = mjd_file + '.gz'
            with gzip.GzipFile(filename=gzip_mjd_file,
                               mode='wb', compresslevel=9) as f:
                if nobs_detect == 1 :
                    # store only 1st detection; no force photo yet.
                    avro_bytes = schema.serialize(alert_first_detect)
                    messg      = schema.deserialize(avro_bytes)
                    schema.store_alerts(f, [alert_first_detect] )
                else:
                    #store alert with previous epochs
                    schema.store_alerts(f, [alert])

                config_data['n_alert_write'] += 1
                print_alert_stats(config_data)
                
        # after writing each aler, copy diasource info to the "past"
        # for the next observation. Make sure to use .copy() to
        # avoid storing a pointer.
        alert['prvDiaSources'].append(alert['diaSource'].copy())

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

def make_outdir_nite(outdir,mjd):

    # + construct name of mjd-specific outdir for this outdir and mjd
    # + create outdir_nite if it does not already exit

    # TODO - pass site name
    nite   = util.get_sunset_mjd(mjd, site=LSST_SITE_NAME)
    niteint = int(nite)
    outdir_nite = outdir + '/' + f"{ALERT_DAY_NAME}" + str(niteint)
    if not os.path.exists(outdir_nite) :
        cmd = f"mkdir {outdir_nite}"
        #print(f"\t Create {outdir_nite}")
        os.system(cmd)

    return outdir_nite
# end make_outdir_nite


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

def print_alert_stats(config_data, done_flag=False):

    n_alert = config_data['n_alert_write']
    n_event = config_data['n_event_write']

    if n_alert % NOBS_ALERT_UPDATE == 0 or done_flag :
        t_start_alert = config_data['t_start_alert']
        t_now         = datetime.datetime.now()
        t_dif_sec  = (t_now - t_start_alert).total_seconds()
        rate       = int(n_alert / t_dif_sec)
        print(f"\t Wrote {n_alert:8d} alerts ({rate}/sec) " \
              f"for {n_event:6d} events.")
        sys.stdout.flush()

    if done_flag:
        print(f"\t Finished writing {n_alert} LSST alerts ({rate}/sec).")
        sys.stdout.flush()

# end print_alert_stats


def write_summary_lsst_alert(name, config_data):
    # write final summary to stdout for "name" of data unit
    done_flag = True
    print_alert_stats(config_data, done_flag)

    truth_dict = config_data['truth_dict']
    if truth_dict :
        write_truth(truth_dict)
    
    # end write_summary_lsst_alert
    
def write_truth(truth_dict):
    
    outfile = truth_dict['outfile'] + '.gz'
    nobs    = truth_dict['nobs_store']
    
    logging.info(f"\n Write {nobs} obs to truth table in {outfile}")

    with gzip.open(outfile, mode='wt') as f:
        
        f.write(f"# SourceID, SNID, MJD, DETECT, " \
                f"TRUE_GENTYPE, TRUE_GENMAG\n")
        
        for o in range(0,nobs):
            diaSourceId  = truth_dict['diaSourceId'][o]
            snid         = truth_dict['snid'][o]
            mjd          = truth_dict['mjd'][o]
            detect       = truth_dict['detect'][o]
            idetect=0
            if detect: idetect = 1
            true_gentype = truth_dict['true_gentype'][o]
            true_genmag  = truth_dict['true_genmag'][o]        

            line = f"{diaSourceId}, {snid}, {mjd:.4f}, {idetect}, " \
                   f"{true_gentype}, {true_genmag:.3f}"
            f.write(f"{line}\n")
            
    # end write_truth
        

