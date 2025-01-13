# Created Jan 2025 by R.Kessler
# Read from F.A.S.T data base (for LSST-DESC) 

import os, sys, glob, yaml, shutil
import logging 
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

        logging.info(f"# ----------------------------------------------- ")
        logging.info(f"Prepare isplit = {isplit} of {nsplit}")

        # break up query into small pieces that can be debugged more easily
        q_season = f"season=1"
        q_nobs   = f"nobs>5"
        q_mod    = f"mod({FASTDB_KEYNAME_SNID},{nsplit})+1={isplit}"            
        q_list   = [ q_season, q_nobs, q_mod ]
        q_join   = " and ".join(q_list)
            
        query = f"select * from dia_object where {q_join}"
        logging.info(f"  Select dia objects with query = \n\t {query}")
        dia_object_all = data_access.submit_short_query(query)
        self.dia_object_all = dia_object_all  # store for later
            
        nobj = len(dia_object_all)
        snid_first = dia_object_all[0][FASTDB_KEYNAME_SNID]
        snid_last  = dia_object_all[nobj-1][FASTDB_KEYNAME_SNID]
            
        logging.info(f"  Found {nobj} objects; first/last SNID = " \
                     f"{snid_first} / {snid_last}")
        logging.info('')
            
        # read all of the light curves on single query
        snid_list = [ obj[FASTDB_KEYNAME_SNID] for obj in dia_object_all ]
        query = "select * FROM dia_source_current WHERE dia_object IN %(objs)s"
        logging.info(f"  Select dia sources with query = \n\t {query}")
        dia_source_all = data_access.submit_short_query( query,
                                                         subdict={'objs': snid_list} )
        self.dia_source_all = dia_source_all
        nsrc = len(dia_source_all)
        logging.info(f"  Found {nsrc} sources among {nobj} objects")

        # finally, load pointer dictionary element for each snid to avoid duplicate
        # searching in read_event.
        uuid_snid_pointer = {}
        ptr_uuid = 0 
        for uuid in dia_source_all:
            snid = str(uuid[FASTDB_KEYNAME_SNID])
            if snid not in uuid_snid_pointer:
                ptr_uuid_min = ptr_uuid                    
            uuid_snid_pointer[snid] = [ ptr_uuid_min, ptr_uuid ]                    
            ptr_uuid += 1

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

        # init output dictionaries
        DEBUG_DUMP = True
        dia_object_all    = self.dia_object_all
        dia_source_all    = self.dia_source_all
        uuid_snid_pointer = self.uuid_snid_pointer
        dia_object_evt    = dia_object_all[evt]
        
        SNID = str(dia_object_evt[FASTDB_KEYNAME_SNID])

        if DEBUG_DUMP:
            logging.debug(f"Read event {evt} : SNID = {SNID}")
        
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
                snana_phot_raw[key_snana] = [ tmp[key_fastdb] for \
                                              tmp in dia_source_evt ]
        else:
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

        # reject this event if there are no observations ...
        # might need a warning somewhere ?
        if nobs == 0:
            snana_data_dict['select'] = False
           
        return snana_data_dict
    
    # end read_event


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
    
