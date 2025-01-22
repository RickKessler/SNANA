# Created Jan 2025 by R.Kessler
# Read from F.A.S.T data base (for LSST-DESC) 

import os, sys, glob, yaml, shutil, time, logging, math, io, copy
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
    gpar.DATAKEY_PHOTFLAG    : 'valid_flag',   # temp until real photflag is there
    gpar.DATAKEY_BAND        : 'filter_name',
    gpar.DATAKEY_FLUXCAL     : 'ps_flux',
    gpar.DATAKEY_FLUXCALERR  : 'ps_flux_err',    
    
    "DUMMY"       : "dummmy"  # no comma here
}


DATAKEY_HEAD_LIST = [ DATAKEY_RA, DATAKEY_DEC ]
DATAKEY_PHOT_LIST = [ DATAKEY_MJD, gpar.DATAKEY_PHOTFLAG,
                      DATAKEY_BAND, DATAKEY_FLUXCAL, DATAKEY_FLUXCALERR ]

FASTDB_KEYNAME_SNID = KEYMAP[DATAKEY_SNID]  
FASTDB_KEYNAME_MJD  = KEYMAP[gpar.DATAKEY_MJD]

FASTDB_ZP = 31.4  # nJy
FLUXSCALE_SNANA = math.pow(10.0, (SNANA_ZP-FASTDB_ZP)/2.5 )

# hard-wired stuff
TABLENAME_DIA_OBJECT = "dia_object"
#TABLENAME_DIA_SOURCE = "dia_source"         # detections only
TABLENAME_DIA_SOURCE = "dia_forced_source"  # forced photo, including detections

QUERY_INTERFACE_LONG  = "LONG"   # recommended
QUERY_INTERFACE_SHORT = "SHORT"
QUERY_INTERFACE       = QUERY_INTERFACE_LONG

QMAXWAIT_OBJECT = 500   # max query wait (seconds) for dia_object query
QMAXWAIT_SOURCE = 3000  # max query wait (seconds) for dia_source query

SELECT_PARAMS_DICT = {

    # cuts to select events
    'NDETECT_MIN'        : 2,
    'MJDDIF_DETECT_MIN'  : 0.5,     # 2 detections must be separated by at least 0.5 days

    # - - - -
    # cuts to select observations
    
    # select obs 30 days before 1st detection, and up to 60 days after last detection
    'MJD_DETECT_SELECT_RANGE' : [ -30, 60 ]
    # PHOTFLAG mask ?? (TBD)
}

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
        
        # read fast data base for isplitran and store in data frame;
        # will parse it later for writing. Objects and sources are stored separately.

        if i_subgroup != 0 :            return -1
        
        data_access = self.data_access
        args        = self.config_inputs['args'] 
        
        logging.info(f"# ----------------------------------------------- ")
        logging.info(f"Prepare isplit = {args.isplitran} of {args.nsplitran}")

        query = self.construct_query_object()
        
        logging.info(f"  Select dia objects with {QUERY_INTERFACE} query = \n\t {query}")
        t0 = time.time()

        if QUERY_INTERFACE == QUERY_INTERFACE_LONG :
            # return csv, then convert to data frame
            csv_obj_tmp     = data_access.synchronous_long_query(query,
                                                             checkeach=5,   # check every 5 sec
                                                             maxwait=QMAXWAIT_OBJECT ) # abort after this time
            df_dia_object  = pd.read_csv(io.StringIO(csv_obj_tmp), sep=',')
        else:
            # return dictionary; then convert to data frame
            dict_obj_tmp      = data_access.submit_short_query(query)
            df_dia_object     = pd.DataFrame(dict_obj_tmp)

        
        util.print_elapsed_time(t0, "Perform dia_object query")
        self.df_dia_object = df_dia_object  # store for later
            
        nobj       = len(df_dia_object)
        snid_list  = list(df_dia_object['dia_object'])
        snid_list  = list(map(int,snid_list))
        logging.info(f"  Found {nobj} objects.")

        snid_first = snid_list[0]
        snid_last  = snid_list[-1]
        logging.info(f"  First/Last SNID = {snid_first} / {snid_last} ")
        logging.info(f"  dia_object columns: {df_dia_object.columns}")
        logging.info('')

        
        query = f"select * FROM {TABLENAME_DIA_SOURCE}  WHERE dia_object IN %(objs)s order by " \
            f"{FASTDB_KEYNAME_SNID}, {FASTDB_KEYNAME_MJD} "
        logging.info(f"  Select dia sources with {QUERY_INTERFACE} query = \n\t {query}")
        t0 = time.time()

        if QUERY_INTERFACE == QUERY_INTERFACE_LONG :
            csv_src_tmp = data_access.synchronous_long_query(query,
                                                           subdict={'objs': snid_list},
                                                           checkeach=5,             # check every 5 sec
                                                           maxwait=QMAXWAIT_SOURCE ) # abort after this time
            df_dia_source  = pd.read_csv(io.StringIO(csv_src_tmp),sep=',')            
        else:
            dict_src_tmp = data_access.submit_short_query( query,
                                                           subdict={'objs': snid_list} )
            df_dia_source = pd.DataFrame(dict_src_tmp)
            
        util.print_elapsed_time(t0, "Perform dia_source query" )

        self.df_dia_source = df_dia_source
        nsrc = len(df_dia_source)
        logging.info(f"  Found {nsrc} sources among {nobj} objects ")
        logging.info(f"  dia_source columns: {df_dia_source.columns}")    

        # - - - - -
            
        logging.info('')
        return nobj

        # end prep_read_data_subgroup
        
    def read_event(self, evt ):

        args            = self.config_inputs['args']
        # init output dictionaries
        DEBUG_DUMP = True
        df_dia_object        = self.df_dia_object
        df_dia_source        = self.df_dia_source

        #import pdb
        #pdb.set_trace()
        
        dia_object_evt    = df_dia_object.iloc[evt]
        
        SNID = str(dia_object_evt[FASTDB_KEYNAME_SNID])
        nobs     = dia_object_evt['nobs']
        
        if DEBUG_DUMP:
            logging.debug(f"Read {args.read_class} event {evt} : SNID = {SNID}")
        
        snana_head_raw, snana_head_calc, snana_head_sim = util.reset_data_event_dict()
        snana_head_raw[DATAKEY_SNID]  = str(SNID)

        for key in DATAKEY_HEAD_LIST:
            snana_head_raw[key]   = dia_object_evt[KEYMAP[key]]

        # - - - - - - --  -
        # load light curve info 
        snana_phot_raw          = self.init_phot_dict(nobs)
        snana_phot_raw['NOBS']  = nobs
        df_lightcurve    = df_dia_source[df_dia_source[FASTDB_KEYNAME_SNID] == int(SNID)]
        
        keylist_snana_phot_store = []
        for key_snana in DATAKEY_PHOT_LIST:
            key_fastdb = KEYMAP[key_snana]
            if key_fastdb in df_lightcurve.columns :
                val_list   = list(df_lightcurve[key_fastdb])
                if 'FLUXCAL' in key_snana:
                    val_list = [ x*FLUXSCALE_SNANA for x in val_list ]

                snana_phot_raw[key_snana] = val_list
                keylist_snana_phot_store.append(key_snana)

        # zero out PHOTFLAG for placeholder 'valid_flag'; detection bits are added below
        if KEYMAP[gpar.DATAKEY_PHOTFLAG] == 'valid_flag':
            snana_phot_raw[gpar.DATAKEY_PHOTFLAG] = [0] * nobs
            
        # - - - - - -
        
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
        
        # Apply event selection requirements
        select = self.select_event(snana_data_dict)
        if select:
            pass  # do NOT set 'select' to True
        else:
            snana_data_dict['select'] = False

        # select observations using MJD cuts
        if select :
            self.select_obs(snana_data_dict, keylist_snana_phot_store )
                            
        return snana_data_dict
    
    # end read_event

    def construct_query_object(self):
        args        = self.config_inputs['args']
        season      = args.season
        nsplit      = args.nsplitran
        isplit      = args.isplitran
        
        # break up query into small pieces that can be debugged more easily

        if season >=0 :
            q_season = f"season={season}"
        else:
            q_season = 'season>=0'
            
        q_nobs   = f"nobs>5"
        
        q_mod    = f"mod({FASTDB_KEYNAME_SNID},{nsplit})+1={isplit}"
        
        
        q_list   = [ q_season, q_nobs, q_mod ]
        q_join   = " and ".join(q_list)

        if args.nevt:
            q_join  += f" order by dia_object limit {args.nevt}"    
        
        query = f"select * from {TABLENAME_DIA_OBJECT} where {q_join} "
        
        return query
    
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
        ndetect = 0
        
        j_first  = -9
        j_second = -9
        j_last   = -9
        
        # make list of detections per observation. This should work for
        # real detection and also for forced snr_detect.
        
        photflag_list = phot_raw[KEY_PHOTFLAG]
        if photflag_list[0] is None: return
        
        #print(f" xxx photflag_list = {photflag_list[0:20]}")
        detect_list   = [ (int(x) & photflag_detect)>0 for x in photflag_list ]
        ndetect       = detect_list.count(True)
        
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
        head_calc[gpar.DATAKEY_NDETECT]           = ndetect
        return

    def select_event(self, snana_data_dict):
        
        # return true if this event passes selection cuts:
        #  * at least 2 detections separated by a night
        #    -> Avoids spurious detections and moving objects

        head_calc         = snana_data_dict['head_calc']        
        MJD_DETECT_FIRST  = head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]
        MJD_DETECT_LAST   = head_calc[gpar.DATAKEY_MJD_DETECT_LAST]
        NDETECT           = head_calc[gpar.DATAKEY_NDETECT]

        NDETECT_MIN       = SELECT_PARAMS_DICT['NDETECT_MIN']
        MJDDIF_DETECT_MIN = SELECT_PARAMS_DICT['MJDDIF_DETECT_MIN']
        
        if NDETECT < NDETECT_MIN :
            return False

        MJDDIF_DETECT = MJD_DETECT_LAST - MJD_DETECT_FIRST
        if MJDDIF_DETECT < MJDDIF_DETECT_MIN:
            return False

        return True
    
    def select_obs(self, snana_data_dict, keylist_snana_phot_store):

        # Truncate MJD-dependent arrays based on MJD cut
        # Inputs:
        #   snana_data_dict : data dictionary; modify phot_raw arrays
        #   keylist_snana_phot_store : list of phot keys to modify
        
        head_calc = snana_data_dict['head_calc']
        phot_raw  = snana_data_dict['phot_raw']   
        
        MJD_DETECT_FIRST = head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]
        MJD_DETECT_LAST  = head_calc[gpar.DATAKEY_MJD_DETECT_LAST]

        MJD_DETECT_SELECT_RANGE = SELECT_PARAMS_DICT['MJD_DETECT_SELECT_RANGE']
        mjdmin           = MJD_DETECT_FIRST + MJD_DETECT_SELECT_RANGE[0]
        mjdmax           = MJD_DETECT_LAST  + MJD_DETECT_SELECT_RANGE[1]

        # store original mjd_list in separate memory using deepcopy
        mjd_list         = copy.deepcopy(phot_raw[gpar.DATAKEY_MJD])
        
        for key_snana in keylist_snana_phot_store :
            val_list_orig = phot_raw[key_snana]
            val_list_trim = [x for x,t in zip(val_list_orig,mjd_list)  if mjdmin<t<mjdmax ]
            phot_raw[key_snana] = val_list_trim  # overwrite input
            
        # - - - - - 
        del mjd_list
        nobs = len(phot_raw[gpar.DATAKEY_MJD])  # updated nobs
        phot_raw['NOBS']  = nobs                      # overwrite input
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
    
