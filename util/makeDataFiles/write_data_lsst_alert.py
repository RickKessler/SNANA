# Created Oct 22, 2021
# write data in lsst-alert format for broker test.
# [R.Kessler, G. Narayan, R.Hlozek, ]
#
# Jan 14 2022 G.Narayan - fix diaObject bug found by Rob K.
# Jan 15 2022 R.Kessler - minor cleanup; start working on reducing output
# Jan 31 2022 RK^2 fix bug setting alert_first_detect
# Feb 22 2022 RK - integreate ZPHOT_Q and check for 2nd hostgal
# Mar 30 2022 RK - fix setting alertId (not alertID) for all epochs
#                   (not just for FIRST_OBS)
# Mar 20 2022 RK - create & write simVersion string
# 
# Apr 09 2022 RK - few updates to run x10 faster:
#   + for sunset-MJD, replace astroplan call with grid search on
#       pre-computed list of sunset-MJDs at CTIO (x5 faster)
#   + use os.mkdir to create nite-dir instead os is.system (x2 faster)
#
# May 12 2022 RK - hostgal_z[_err] -> hostgal_zspec[_err] and add
#                  hostgal_zphot[_err]
#
# ==================================================

import os, sys, glob, gzip, math, yaml, json
import datetime, logging, shutil, subprocess
from copy import copy
from pathlib import Path

import lsst.alert.packet
import numpy as np
from fastavro import reader, writer

import makeDataFiles_params as gpar
import makeDataFiles_util as util
#from makeDataFiles_params import *

# ========================================================
# map dictionary(SNANA) varName to alert varName
lc = "lc"  # instruction to take lower case of dict value
VARNAME_DIASRC_MAP = {
    gpar.DATAKEY_SNID            : 'diaObjectId',
    gpar.DATAKEY_RA              : lc,
    gpar.DATAKEY_DEC             : 'decl',
    gpar.DATAKEY_zHEL            : 'z_final' ,
    gpar.DATAKEY_zHEL_ERR        : 'z_final_err',
    gpar.DATAKEY_NOBS            : lc       # in phot_raw, not header
    #
}

VARNAME_DIAOBJ_MAP = {
    gpar.DATAKEY_SNID            : 'diaObjectId',
    gpar.DATAKEY_RA              : lc,
    gpar.DATAKEY_DEC             : 'decl',
    gpar.DATAKEY_MWEBV           : lc,
    gpar.DATAKEY_MWEBV_ERR       : lc,
    gpar.HOSTKEY_SNSEP           : lc,
    gpar.HOSTKEY_SPECZ           : 'hostgal_zspec',
    gpar.HOSTKEY_SPECZ_ERR       : 'hostgal_zspec_err',
    gpar.HOSTKEY_PHOTOZ          : 'hostgal_zphot',
    gpar.HOSTKEY_PHOTOZ_ERR      : 'hostgal_zphot_err'
}


VARNAME_OBS_MAP = {
    'MJD'        : 'midPointTai',
    'BAND'       : 'filterName',
    'FLUXCAL'    : 'psFlux',
    'FLUXCALERR' : 'psFluxErr'
}

LSST_ZP_nJy     = 31.4   # report calibrated flux in this unit

# FLXUCAL(ALERT) = FLUXCAL(SNANA)*SCALE_FLUXCAL
ARG_ZPDIF       = 0.4*(LSST_ZP_nJy-gpar.SNANA_ZP)
SCALE_FLUXCAL   = math.pow(10.0,ARG_ZPDIF)
KEYNAME_SUBSTRING_FLUXCAL = 'FLUXCAL'  # scale variables with this substring

PHOTFLAG_DETECT = 4096   # should read this from global data header ??

NOBS_ALERT_MAX  = 2000   # used to compute diaSource
NOBS_ALERT_UPDATE = 100  # std update after this many alerts

ALERT_DAY_NAME    = "NITE"
LSST_SITE_NAME    = "CTIO"

#'alertId', 'l1dbId', 'diaSource', 'prvDiaSources', 'diaObject', 'ssObject']

TIME_WAIT_FORCEPHOTO = 0.9  # days,  wait this long for previous sources

# ===============================================================
def append_HOSTGAL_DIAOBJ_MAP():

    # Created Feb 11 2022 by R.Kessler
    # After reading columm names, add more hostgal variables to
    # global VARNAME_DIAOBJ_MAP
    
    for prefix in ['HOSTGAL_MAG', 'HOSTGAL_MAGERR'] :
        for band in list(gpar.SURVEY_INFO['FILTERS']['LSST']):
            key_snana = f"{prefix}_{band}"
            key_alert = f"{prefix.lower()}_{band}"
            VARNAME_DIAOBJ_MAP[key_snana] = key_alert

    for key in gpar.DATAKEY_LIST_ZPHOT_Q :
        if 'HOSTGAL2' in key: continue  # 2nd host treated elsewhere
        key_snana = key
        key_alert = lc
        VARNAME_DIAOBJ_MAP[key_snana] = key_alert

    # for each HOSTGAL_XXX key, add another key with HOSTGAL2_XXX
    TMP_MAP = VARNAME_DIAOBJ_MAP.copy()
    for key_snana, key_alert in TMP_MAP.items() :
        if gpar.HOSTKEY_BASE in key_snana :
            key2_snana = util.key_hostgal_nbr(key_snana,2)
            key2_alert = util.key_hostgal_nbr(key_alert,2)
            VARNAME_DIAOBJ_MAP[key2_snana] = key2_alert

    #sys.exit(f"\n VARNAME_DIAOBJ_MAP = \n{VARNAME_DIAOBJ_MAP}")
    
    # end append_HOSTGAL_DIAOBJ_MAP
    
def init_schema_lsst_alert(schema_file, snana_folder):

    schema     = lsst.alert.packet.Schema.from_file(filename=schema_file)
    schema_dir = os.path.dirname(schema_file)

    # beware: too much hard coding
    json_file  = f"{schema_dir}/sample_data/plasticc.json"  

    logging.info(f"\n Init alert schema based on\n" \
                 f"\t schema_file={schema_file}\n" \
                 f"\t jon_file={json_file} ")

    # Load an example json alert, and clear the numberical input
    with open(json_file,'r') as f:
        alert_data = json.load(f)
    
    # - - - - - - - - - - - - - 
    # Mar 30 2022: construct simVersion string
    key_list   = ['TIME_START', 'SNANA_VERSION', 'SIMLIB_FILE' ]
    value_dict = util.extract_sim_readme_info(snana_folder, key_list);
    
    TIME_START    = value_dict['TIME_START']
    SNANA_VERSION = value_dict['SNANA_VERSION']
    SIMLIB_FILE   = value_dict['SIMLIB_FILE']

    base          = os.path.basename(SIMLIB_FILE)
    cadence       = base.rsplit('.',1)[0]
    t_start       = str(TIME_START).split()[0] # just date; leave out hr:min:sec
    
    # get schema version
    with open(schema_file,'r') as s:
        content        = json.load(s)
        schema_version = content['namespace'].split('.')[-1]

    simVersion = f"date({t_start})_" \
                 f"snana({SNANA_VERSION})_" \
                 f"cadence({cadence})_" \
                 f"schema({schema_version})"
      
    logging.info(f" simVersion = {simVersion}\n")

    alert_data['diaObject']['simVersion'] = simVersion
    
    return schema, alert_data, simVersion

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
    head_sim  = data_event_dict['head_sim']
    phot_raw  = data_event_dict['phot_raw']

    # scale SNANA-sim flux[Err] to alert flux unit (nJy)
    for key in VARNAME_OBS_MAP:
        if KEYNAME_SUBSTRING_FLUXCAL in key:
            phot_raw[key] = [f*SCALE_FLUXCAL for f in phot_raw[key] ]

    # - - - -
    SNID         = int(head_raw[gpar.DATAKEY_SNID]) # to compare sourceID
    NOBS         = phot_raw[gpar.DATAKEY_NOBS]
    true_gentype  = head_sim[gpar.SIMKEY_TYPE_INDEX]

    # strip off number of processed events; init stuff on nevent=0
    data_unit_name        = data_event_dict['data_unit_name']
    index_unit            = data_event_dict['index_unit']
    data_unit_name_list   = config_data['data_unit_name_list']
    data_unit_nevent_list = config_data['data_unit_nevent_list']
    nevent        = data_unit_nevent_list[index_unit]
    outdir        = args.outdir_lsst_alert

    if nevent == 0 :
        append_HOSTGAL_DIAOBJ_MAP()
        
        schema, alert_data, simVersion  = \
            init_schema_lsst_alert(args.lsst_alert_schema,args.snana_folder)
        alert_data_orig     = alert_data.copy()
        config_data['schema']          = schema
        config_data['alert_data_orig'] = alert_data_orig
        config_data['diaSourceId']     = 0
        config_data['n_alert_write']   = 0
        config_data['n_event_write']   = 0
        config_data['t_start_alert']   = datetime.datetime.now()
        config_data['simVersion']      = simVersion
        
        outfile = args.outfile_alert_truth
        if outfile :
            config_data['truth_dict'] = init_truth_dict(outfile)
        else:
            config_data['truth_dict'] = None

        # check for file with pre-computed list of sunset MJDs (Apr 9 2022)
        mjd_sunset_file = args.mjd_sunset_file
        if mjd_sunset_file :
            mjd_sunset_dict = { 'mjd_file' : mjd_sunset_file }
        else:
            mjd_sunset_dict = {}
        config_data['mjd_sunset_dict'] = mjd_sunset_dict
        
    # - - - - -
    schema             = config_data['schema']
    diaSourceId        = config_data['diaSourceId']
    alert_data_orig    = config_data['alert_data_orig']
    alert              = copy(alert_data_orig)
    alert_first_detect = copy(alert_data_orig)

    # copy structure of original sample alert to local diaSource
    prvDiaSources = alert_data_orig['prvDiaSources']
    diaObject     = alert_data_orig['diaObject']
    diaSource     = prvDiaSources[0]

    # # print this out for testing
    #for key in diasrc.keys():
    #    print(f" xxx found original alert key: {key}")

    alert['prvDiaSources'].clear() # clear out all the past histories
    alert_first_detect['prvDiaSources'].clear() # clear out all the past histories

    # - - - - - -
    my_diaSource = diaSource    # not empty - has default quantities from schema
    my_diaObject = diaObject    # not empty - has default quantites from schema

    
    # translate snana header and create diaSource[Object] dictionaries for lsst alert
    translate_dict_alert(-1, data_event_dict,           # <== I
                         my_diaSource, my_diaObject)    # <== I/O
    
    alert['diaSource']              = my_diaSource
    alert['diaObject']              = my_diaObject

    # idiot check ...
    diaObjectId = my_diaObject['diaObjectId']
    check_ObjId = my_diaSource['diaObjectId']
    if diaObjectId != check_ObjId:
        msgerr = ''
        msgerr += f"ObjectId={diaObjectId} from my_diaObject, but\n"
        msgerr += f"ObjectId={check_ObjId} from my_diaSource\n"
        msgerr += f"-> Something is messed up."
        raise ValueError(msgerr)

    config_data['n_event_write'] += 1
    FIRST_OBS    = True
    nobs_detect  = 0
    nobs_keep    = 0 # nobs_detect + nobs_forcePhoto

    if args.nite_detect_range:
        MJD_REF  = head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]
        MJD_LAST = head_calc[gpar.DATAKEY_MJD_DETECT_LAST]
        TIME_BACK_FORCE  = 30   #Ndays before MJD_REF to include forced phot.
        TIME_FORWARD_FORCE = MJD_LAST - MJD_REF + 0.1
    elif args.peakmjd_range:
        MJD_REF = head_calc[gpar.DATAKEY_PEAKMJD]
        TIME_BACK_FORCE  = 50   #Ndays before MJD_REF to include forced phot.
        TIME_FORWARD_FORCE = 100

    # compute things for each obs
    mjd_list         = data_event_dict['phot_raw']['MJD']
    photflag_list    = data_event_dict['phot_raw']['PHOTFLAG']
    true_genmag_list = data_event_dict['phot_raw'][gpar.VARNAME_TRUEMAG]

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
    #translate each obs to diaSrc dictionary
    for o in range(o_start,o_end):
        keep_force = keep_force_list[o]
        if not keep_force: continue

        nobs_keep += 1
        mjd       = data_event_dict['phot_raw']['MJD'][o]
        detect    = detect_list[o]

        # compute UNIQUE diaSource from already unique SNID
        diaSourceId = NOBS_ALERT_MAX*SNID + o
        my_l1dbId   = 0                           # we don't know what this is
        my_diaSource['diaSourceId'] = diaSourceId

        # update Source Object with this obs
        translate_dict_alert(o, data_event_dict, my_diaSource, my_diaObject)

        if FIRST_OBS  :
            # Save my_diasrc info on 1st observation
            alert['diaSource'] = my_diaSource
            alert['diaObject'] = my_diaObject
            alert['l1dbId']    = my_l1dbId
            FIRST_OBS = False

        alert['alertId']   = diaSourceId
        
        # serialize the alert
        avro_bytes = schema.serialize(alert)
        messg      = schema.deserialize(avro_bytes)

        # write truth table for detections and forced photo
        if args.outfile_alert_truth:
            truth_dict = config_data['truth_dict']
            truth_dict['nobs_store'] += 1
            truth_dict['diaSourceId'].append(diaSourceId)
            truth_dict['snid'].append(SNID)
            truth_dict['mjd'].append(mjd)
            truth_dict['detect'].append(detect)
            truth_dict['true_gentype'].append(true_gentype)
            truth_dict['true_genmag'].append(true_genmag_list[o])

        # write alerts ONLY for detection.
        if detect :
            nobs_detect += 1

            # construct name of avro file using mjd, objid, srcid
            mjd_sunset_dict = config_data['mjd_sunset_dict']
            outdir_nite  = make_outdir_nite(outdir, mjd, mjd_sunset_dict)
            str_day = f"mjd{mjd:.4f}"
            str_obj = f"obj{diaObjectId}"
            str_src = f"src{diaSourceId}"
            mjd_file  = f"{outdir_nite}/" \
                        f"alert_{str_day}_{str_obj}_{str_src}.avro"

            delta_t = mjd - MJD_REF  # elapsed time since 1st detect
            gzip_mjd_file = mjd_file + '.gz'
            with gzip.GzipFile(filename=gzip_mjd_file,
                               mode='wb', compresslevel=9) as f:

                # xxx mark delete if nobs_detect == 1 :
                # ignore previous sources within ~1day of 1st detection
                # because it takes 24 hr to run forcePhoto

                if delta_t < TIME_WAIT_FORCEPHOTO :
                    # store only 1st detection; no force photo yet.
                    alert_first_detect['alertId']   = diaSourceId
                    alert_first_detect['diaSource'] = copy(my_diaSource)
                    alert_first_detect['diaObject'] = copy(my_diaObject)
                    alert_first_detect['prvDiaSources'].clear()

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
        logging.warn(f"No detections found for SNID={SNID}; " \
                     f"hack alert to avoid crash")
        o = 0
        translate_dict_alert(o, data_event_dict, my_diaSource, my_diaObject)
        alert['diaSource'] = my_diaSource
        alert['alertId']   = my_diaSource
        alert['prvDiaSources'].append(alert['diaSource'])

    return

# end write_event_lsst_alert

def make_outdir_nite(outdir,mjd, mjd_sunset_dict):

    # + construct name of mjd-specific outdir for this outdir and mjd
    # + create outdir_nite if it does not already exit
    
    nite   = util.get_sunset_mjd(mjd, LSST_SITE_NAME, mjd_sunset_dict )
    niteint = int(nite)
    outdir_nite = outdir + '/' + f"{ALERT_DAY_NAME}" + str(niteint)
    if not os.path.exists(outdir_nite) :
        os.mkdir(outdir_nite)

        # xxx mark cmd = f"mkdir {outdir_nite}"
        # xxx mark os.system(cmd)

    return outdir_nite
# end make_outdir_nite


def translate_dict_alert(obs, data_event_dict, diaSource, diaObject):

    # obs = -1 -> translate header info in data_event_dict to diaSource
    # obs >= 0 -> translate obs in data_event_dict to diaSource.
    #               Take care for speed optimization

    if obs < 0 : # load the header info

        # set flag to store host info only if there is host info;
        # e.g., host info should never be set for Galactic transients.
        head_raw       = data_event_dict['head_raw']
        HOSTGAL_NMATCH = head_raw[gpar.HOSTKEY_NMATCH]

        for varName_inp in VARNAME_DIASRC_MAP:  # loop over SNANA keys
            varName_avro = VARNAME_DIASRC_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()

            diaSource[varName_avro] = \
                get_data_alert_value(data_event_dict, varName_inp)

        for varName_inp in VARNAME_DIAOBJ_MAP:
            varName_avro = VARNAME_DIAOBJ_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()

            diaObject[varName_avro] = \
                get_data_alert_value(data_event_dict, varName_inp)

            #print(f" xxx {varName_inp} -> {varName_avro} ")

    else: # load the observations
        phot_raw  = data_event_dict['phot_raw']
        for varName_inp in VARNAME_OBS_MAP:
            varName_avro = VARNAME_OBS_MAP[varName_inp]
            if varName_avro == lc:  varName_avro = varName_inp.lower()
            diaSource[varName_avro] = phot_raw[varName_inp][obs]

    # end translate_dict_alert

def get_data_alert_value(data_event_dict, varName):
    # return data value for varName
    # Check head_raw, head_calc and phot_raw.

    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']
    phot_raw  = data_event_dict['phot_raw']
    value     = None

    if varName in head_raw:
        if varName == gpar.DATAKEY_SNID: # convert str to int
            value = int(head_raw[varName])
        else:
            value = head_raw[varName]

    elif varName in head_calc :
        value = head_calc[varName]
    else:
        value = phot_raw[varName]

    # - - - - - - 
    # format with fewer digits to reduce size of output files (Apr 6 2022)
    # Use SNANA key names (not avro keys)
    #if 'HOSTGAL_MAG' in varName or 'REDxSHIFT' in varName:
    #    ival   = int(1.0E4 * value + 0.5)
    #    value  = float(ival) * 1.0E-4  # doesn't work ?????
        
    return value
    # end get_data_alert_value

def print_alert_stats(config_data, done_flag=False):

    n_alert = config_data['n_alert_write']
    n_event = config_data['n_event_write']

    if n_alert % NOBS_ALERT_UPDATE == 0 or done_flag :
        t_start_alert = config_data['t_start_alert']
        t_now         = datetime.datetime.now()
        t_dif_sec  = (t_now - t_start_alert).total_seconds()
        rate       = int(n_alert / t_dif_sec)
        logging.info(f"\t Wrote {n_alert:8d} alerts ({rate}/sec) " \
                     f"for {n_event:6d} events.")
        sys.stdout.flush()

    if done_flag:
        logging.info(f"\t Finished writing {n_alert} LSST alerts ({rate}/sec).")
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


