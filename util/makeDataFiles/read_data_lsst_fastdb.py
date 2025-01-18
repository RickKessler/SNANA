# Created Jan 2025 by R.Kessler
# Read from F.A.S.T data base (for LSST-DESC) 

import os, sys, glob, yaml, shutil, time, logging, math
import numpy  as np
import pandas as pd
import makeDataFiles_util  as    util
import makeDataFiles_params as gpar

from   makeDataFiles_base    import Program
from   makeDataFiles_params  import *

try:
    path_fastdb_api=os.path.expandvars("$TD_SOFTWARE/tom_deployment/tom_desc_fastdb_dev/fastdb_api")
    sys.path.insert(1, path_fastdb_api)
    from fastdb_api import FASTDB
except:
    pass

# hard wire here for now, but should read from map SNANA <--> FASTDB

KEYMAP = { # SNANA -> fastdb

    # header
    gpar.DATAKEY_SNID  : "dia_object",
    gpar.DATAKEY_RA    : "ra",
    gpar.DATAKEY_DEC   : "decl" ,

    # PHOT
    gpar.DATAKEY_MJD         : 'mid_point_tai',
    gpar.DATAKEY_BAND        : 'filter_name',
    gpar.DATAKEY_FLUXCAL     : 'ps_flux',
    gpar.DATAKEY_FLUXCALERR  : 'ps_flux_err',    
    
    "DUMMY"       : "dummmy"  # no comma here
}

DATAKEY_HEAD_LIST = [ DATAKEY_RA, DATAKEY_DEC ]
DATAKEY_PHOT_LIST = [ DATAKEY_MJD, DATAKEY_BAND, DATAKEY_FLUXCAL, DATAKEY_FLUXCALERR ]

FASTDB_KEYNAME_SNID = KEYMAP[DATAKEY_SNID]  

FASTDB_ZP = 31.4  # nJy
FLUXSCALE_SNANA = math.pow(10.0, (SNANA_ZP-FASTDB_ZP)/2.5 )


TABLENAME_DIA_OBJECT = "dia_object"
TABLENAME_DIA_SOURCE = "dia_source"         # detections only
#TABLENAME_DIA_SOURCE = "dia_forced_source"  # forced photo, including detections


# ======================================================

class data_lsst_fastdb(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        logging.info("Init data_lsst_fastdb class.")
        gpar.PREFIX_SEASON = "SEASON"  # intermediate prefix on data files
        
        super().__init__(config_inputs, config_data)

    def init_read_data(self):

        args = self.config_inputs['args']  # command line args
        logging.info("Connect to FASTDB\n")
        self.data_access = FASTDB()
    
        return
        # end init_read_data

        
    def prep_read_data_subgroup(self, i_subgroup):
        #

        if i_subgroup != 0 :            return -1
        
        data_access = self.data_access
        args        = self.config_inputs['args'] 
        nsplit = args.nsplitran
        isplit = args.isplitran
        nevt   = args.nevt
        
        logging.info(f"# ----------------------------------------------- ")
        logging.info(f"Prepare isplit = {isplit} of {nsplit}")

        # break up query into small pieces that can be debugged more easily
        if args.season >=0 :
            q_season = f"season={args.season}"
        else:
            q_season = 'season>=0'
            
        q_nobs   = f"nobs>5"
        q_mod    = f"mod({FASTDB_KEYNAME_SNID},{nsplit})+1={isplit}"            
        q_list   = [ q_season, q_nobs, q_mod ]
        q_join   = " and ".join(q_list)
            
        query = f"select * from {TABLENAME_DIA_OBJECT} where {q_join}"
        logging.info(f"  Select dia objects with query = \n\t {query}")
        t0 = time.time()
        dia_object_all = data_access.submit_short_query(query)
        util.print_elapsed_time(t0, "Perform dia_object query")
        self.dia_object_all = dia_object_all  # store for later
            
        nobj       = len(dia_object_all)
        snid_list  = [ obj[FASTDB_KEYNAME_SNID] for obj in dia_object_all ]        
            
        logging.info(f"  Found {nobj} objects.")
            
        # read all of the light curves on single query
        if nevt < nobj:
            nobj = nevt
            logging.info(f"  Truncate number of objects to {nevt} (--nevt arg)")
            snid_list = snid_list[0:nevt]

        snid_first = snid_list[0]
        snid_last  = snid_list[-1]
        logging.info(f"  First/Last SNID = {snid_first} / {snid_last} ")                        
        logging.info('')
        
        query = f"select * FROM {TABLENAME_DIA_SOURCE}  WHERE dia_object IN %(objs)s"
        logging.info(f"  Select dia sources with query = \n\t {query}")
        t0 = time.time()
        dia_source_all = data_access.submit_short_query( query,
                                                         subdict={'objs': snid_list} )
        util.print_elapsed_time(t0, "Perform dia_source query" )

        self.dia_source_all = dia_source_all
        nsrc = len(dia_source_all)
        logging.info(f"  Found {nsrc} sources among {nobj} objects ")

        # finally, load pointer dictionary element for each snid to avoid duplicate
        # searching in read_event.
        t0 = time.time()
        uuid_snid_pointer = {}
        ptr_uuid = 0 
        for uuid in dia_source_all:
            snid = str(uuid[FASTDB_KEYNAME_SNID])
            if snid not in uuid_snid_pointer:
                ptr_uuid_min = ptr_uuid                    
            uuid_snid_pointer[snid] = [ ptr_uuid_min, ptr_uuid ]                    
            ptr_uuid += 1
        util.print_elapsed_time(t0, "Set pointers in returned dia_source list" )

        CHECK_NPTR = False
        if CHECK_NPTR:
            nptr = len(uuid_snid_pointer)
            if nptr != nobj :
                errmsg = [ f"Found {nobj} objects, but found {nptr} snid pointers",
                           f"Something is WRONG !!" ]                
                util.log_assert( False, errmsg )

        self.uuid_snid_pointer = uuid_snid_pointer
        logging.info('')
        return nobj

        # end prep_read_data_subgroup

    
    def read_event(self, evt ):

        args            = self.config_inputs['args']
        # init output dictionaries
        DEBUG_DUMP = True
        dia_object_all    = self.dia_object_all
        dia_source_all    = self.dia_source_all
        uuid_snid_pointer = self.uuid_snid_pointer
        dia_object_evt    = dia_object_all[evt]
        
        SNID = str(dia_object_evt[FASTDB_KEYNAME_SNID])

        if DEBUG_DUMP:
            logging.debug(f"Read {args.read_class} event {evt} : SNID = {SNID}")
        
        snana_head_raw, snana_head_calc, snana_head_sim = util.reset_data_event_dict()
        snana_head_raw[DATAKEY_SNID]  = str(SNID)

        for key in DATAKEY_HEAD_LIST:
            snana_head_raw[key]   = dia_object_evt[KEYMAP[key]]

        # - - - - - - --  -
        # load light curve info using pointers to extract list of dia-dictionary elements

        if SNID in uuid_snid_pointer :
            ptrmin    = uuid_snid_pointer[SNID][0]
            ptrmax    = uuid_snid_pointer[SNID][1]
            nobs      = ptrmax - ptrmin + 1
            
            snana_phot_raw          = self.init_phot_dict(nobs)
            snana_phot_raw['NOBS']  = nobs
            dia_source_evt          = dia_source_all[ptrmin:ptrmax+1]
            
            for key_snana in DATAKEY_PHOT_LIST:
                key_fastdb = KEYMAP[key_snana]
                scale      = 1.0                
                if 'FLUXCAL' in key_snana:
                    scale = FLUXSCALE_SNANA
                    snana_phot_raw[key_snana] = [ scale * tmp[key_fastdb] for \
                                                  tmp in dia_source_evt ]
                else:
                    snana_phot_raw[key_snana] = [ tmp[key_fastdb] for \
                                                  tmp in dia_source_evt ] 
        else :
            # event with no observations; perhaps with cuts such as PSF/ZP/PHOTFLAG ?
            snana_phot_raw = self.init_phot_dict(0)

        # - - - - - -
        nobs = snana_phot_raw['NOBS']
        
        if DEBUG_DUMP:
            str_warn = ''
            if nobs == 0: str_warn = "BEWARE NO OBSERVATIONS !!!"
            logging.debug(f"\t store {nobs} sources for SNID = {SNID}   {str_warn}")

        # construct SNANA dictionary
        snana_data_dict = {
            'head_raw'  : snana_head_raw,
            'head_calc' : snana_head_calc,
            'phot_raw'  : snana_phot_raw
        }


        # check option to fudge detections with SNR; this is a placeholder
        # until detections are passed from fastdb
        if args.snr_detect:
            self.force_snr_detections(snana_data_dict)

        # store first/second/last MJD_DETECT            
        self.store_mjd_detections(snana_data_dict)  
        
        # reject this event if there are no detections
        MJD_DETECT_FIRST = snana_data_dict['head_calc'][gpar.DATAKEY_MJD_DETECT_FIRST]
        if MJD_DETECT_FIRST < 1000.0 :
            snana_data_dict['select'] = False
            
        return snana_data_dict
    
    # end read_event

    def force_snr_detections(self, snana_data_dict):
        # if args.snr_detect is set, compute detection for each obs and set
        # args.photflag_detect mask of PHOTFLAG.
        #
        
        args            = self.config_inputs['args']
        snr_detect      = args.snr_detect
        photflag_detect = args.photflag_detect
        if snr_detect is None: return 
        
        head_raw   = snana_data_dict['head_raw']
        head_calc  = snana_data_dict['head_calc']
        phot_raw   = snana_data_dict['phot_raw']
        
        KEY_PHOTFLAG    = gpar.DATAKEY_PHOTFLAG
        KEY_FLUXCAL     = gpar.DATAKEY_FLUXCAL
        KEY_FLUXCALERR  = gpar.DATAKEY_FLUXCALERR

        fluxcal_list    = phot_raw[KEY_FLUXCAL]
        fluxcalerr_list = phot_raw[KEY_FLUXCALERR]
        snr_list        = [x / y for x, y in zip(fluxcal_list, fluxcalerr_list)]
        detect_list     = [x > snr_detect for x in snr_list ]
        #self.detect_list = detect_list
        
        if True in detect_list:
            photflag_list   = [ photflag_detect if x else 0 for x in detect_list ]
        else:
            photflag_list = [0] * len(snr_list)


        phot_raw[KEY_PHOTFLAG] = photflag_list            
        #print(f" xxx -----------")
        #print(f" xxx force_snr_detect: photflag_list = \n{photflag_list[0:20]}")
        
        return 

    def store_mjd_detections(self,snana_data_dict):

        # compute and store MJD_DETECT_[FIRST/SECOND/LAST]
        args            = self.config_inputs['args']        
        photflag_detect = args.photflag_detect
        
        head_calc    = snana_data_dict['head_calc']
        phot_raw     = snana_data_dict['phot_raw']
        KEY_MJD      = gpar.DATAKEY_MJD
        KEY_PHOTFLAG = gpar.DATAKEY_PHOTFLAG

        if KEY_PHOTFLAG not in phot_raw : return

        mjd_detect_first  = -9.0
        mjd_detect_second = -9.0
        mjd_detect_last   = -9.0

        j_first  = -9
        j_second = -9
        j_last   = -9
        
        # make list of detections per observation. This should work for
        # real detection and also for forced snr_detect.
        
        photflag_list = phot_raw[KEY_PHOTFLAG]
        if photflag_list[0] is None: return
        
        #print(f" xxx photflag_list = {photflag_list[0:20]}")
        detect_list   = [ (int(x) & photflag_detect)>0 for x in photflag_list ]
        
        # find MJD for first/second/last detection
        if True in detect_list:        
            j_first  = detect_list.index(True)
            j_last   = len(detect_list) - detect_list[::-1].index(True) - 1

            mjd_detect_first  = phot_raw[KEY_MJD][j_first]
            mjd_detect_last   = phot_raw[KEY_MJD][j_last]
            
            tmp_list = detect_list[j_first+1:]
            if True in tmp_list:
                j_second = tmp_list.index(True) + j_first+1
                mjd_detect_second = phot_raw[KEY_MJD][j_second]
                
            #print(f" xxx j = {j_first}  {j_second}  {j_last}")
        
        #print(f" xxx mjd_detect = {mjd_detect_first}  {mjd_detect_second}  {mjd_detect_last} ")
        # store output
        head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]  = mjd_detect_first
        head_calc[gpar.DATAKEY_MJD_DETECT_SECOND] = mjd_detect_second        
        head_calc[gpar.DATAKEY_MJD_DETECT_LAST]   = mjd_detect_last        

        return
    
    def end_read_data_subgroup(self):
        pass
    def end_read_data(self):
        # global end for reading data                           
        pass

    def exclude_varlist_obs(self):
        # return list of PHOT columns to excude from output text files
        # see VARNAMES_OBS in makeDataFiles_params.py
        return [ 'PSF_SIG1', 'ZEROPT',  'SKY_SIG', 'XPIX', 'YPIX', 'GAIN',
                 'CCDNUM', 'IMGNUM' ]
    
