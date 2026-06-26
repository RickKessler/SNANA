# Created Jan 2025 by R.Kessler
# Read from F.A.S.T data base (for LSST-DESC) 
#
# Jun ?? 2026: begin complete overhaual with real data
# Jun 26 2026: minor overhaul after Rob-Refac on fastdb

import os, sys, glob, yaml, shutil, time, logging, math, io, copy
import numpy  as np
import pandas as pd

import makeDataFiles_util   as  util
import makeDataFiles_params as  gpar

from   makeDataFiles_base    import Program
from   makeDataFiles_params  import *

# use try/except for imports that might crash on other data sets
try:
    import pyarrow
    from  fastdb_client import FASTDBClient
    from  astropy.table import Table
    from  astropy.table import vstack
except:
    pass


# =======================================================================
# hard-wired stuff for fastdb access


MXOBJ_PER_FETCH = 1000 ; # max number of objects per source query

SURVEY_LSST     = 'LSST'
FILTERLIST_LSST = SURVEY_INFO['FILTERS'][SURVEY_LSST]
BAND_LIST_LSST  = list(FILTERLIST_LSST)

# - - - -
# hard wire key name map here for now, but should read from map SNANA <--> FASTDB
KEYMAP = { # SNANA -> fastdb

    # header
    gpar.DATAKEY_SNID           : 'rootid',
    gpar.DATAKEY_NAME_TRNS      : 'diaobjectid',
    gpar.DATAKEY_RA             : 'ra',
    gpar.DATAKEY_DEC            : 'dec' ,

    # PHOT
    gpar.DATAKEY_MJD         : 'mjd',
    gpar.DATAKEY_PHOTFLAG    : 'photflag',  # temp dummy  
    gpar.DATAKEY_BAND        : 'band',
    gpar.DATAKEY_FLUXCAL     : 'flux',
    gpar.DATAKEY_FLUXCALERR  : 'fluxerr',    
    
    "DUMMY"       : "dummmy"  # no comma here
}


PHOTFLAG_DETECT       = 1    # set this PHOTFLAG bit for isdet detection
PHOTFLAG_GARBAGE      = 64   # found garbage value

TOL_MJD_NIGHT     = 0.4  # coadd obs withing this tolerance, days

CUT_NOBS_COADD_RATIO = 0.4  # separates WFD and DDF

DATAKEY_HEAD_LIST = [ DATAKEY_RA, DATAKEY_DEC ]
DATAKEY_PHOT_LIST = [ DATAKEY_MJD, gpar.DATAKEY_PHOTFLAG,
                      DATAKEY_BAND, DATAKEY_FLUXCAL, DATAKEY_FLUXCALERR ]

FASTDB_KEYNAME_ROOTID  = KEYMAP[gpar.DATAKEY_SNID]
FASTDB_KEYNAME_OBJID   = KEYMAP[gpar.DATAKEY_NAME_TRNS]

#FASTDB_KEYNAME_MJD  = KEYMAP[gpar.DATAKEY_MJD]

FASTDB_ZP = 31.4  # nJy
FLUXSCALE_SNANA = math.pow(10.0, (SNANA_ZP-FASTDB_ZP)/2.5 )


# ======================================================

class data_lsst_fastdb(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        logging.info("Init data_lsst_fastdb class.")
        gpar.PREFIX_SEASON = "CHUNK"  # intermediate prefix on data files

        super().__init__(config_inputs, config_data)

    def init_read_data(self):
                
        args = self.config_inputs['args']  # command line args        

        if args.photflag_detect == 0:
            args.photflag_detect = PHOTFLAG_DETECT
        if args.photflag_garbage == 0 :
            args.photflag_garbage = PHOTFLAG_GARBAGE
                    
        logging.info('')
        logging.info("Begin init_read_data")        
        logging.info("Connect to FASTDBClient")
        fdb = FASTDBClient( "production" )

        # get total number of light curves
        res = fdb.post( "count/rootid/realtime" )
        nlc_tot = res['count']
        logging.info(f" Total number of light curves in fastdb:  {nlc_tot} ")
        if args.nevt < nlc_tot:      nlc_tot = args.nevt
        logging.info(f" Total number of light curves to extract:  {nlc_tot} ")

        n_subgroup = int( ( nlc_tot-1) / MXOBJ_PER_FETCH) + 1
        logging.info(f" Will read {n_subgroup} groups of {MXOBJ_PER_FETCH} objects")
        logging.info(f" coadd_by_nite: {args.coadd_by_nite} ")
        logging.info(f" PHOTFLAG[DETECT,GARBAGE] = " \
                     "{args.photflag_detect} {args.photflag_garbage}")
        
        self.fdb        = fdb
        self.nlc_tot    = nlc_tot
        self.n_subgroup = n_subgroup
        self.garbage_table_extra_columns = { 'ltcvs_offset' : 0, 'ltcsv_limit': 0 }

        self.nobs_coadd_ratio = []  # for diagnostic histogram to check WFD vs. DDF
        
        logging.info("Done with init_read_data")
        logging.info('')

        return
    # end init_read_data

        
    def prep_read_data_subgroup(self, i_subgroup):

        # read fast data base for this subgroup;
        # will parse it later for writing.
        # Objects and sources are stored separately.
        
        args          = self.config_inputs['args']
        fdb           = self.fdb
        nlc_tot       = self.nlc_tot
        n_subgroup    = self.n_subgroup

        nsplit        = args.nsplitran
        isplit_select = args.isplitran  # 1 to nsplit, or -1 for all                                               
        n_fetch    = min(MXOBJ_PER_FETCH,args.nevt)
        n_offset   = i_subgroup * MXOBJ_PER_FETCH

        if n_offset + n_fetch > args.nevt:
            n_fetch = args.nevt - n_offset
            
        if n_offset >= nlc_tot:
            return -1  # done reading

        # check split option to read this group (typically to split jobs among cores)
        if nsplit > 1:
            match_split = util.select_split(i_subgroup, args)            
        else:
            match_split = True

        # - - - -
        verb = "Begin" if match_split else "Skip "
        logging.info(f"# ---------------------------------------------------------------- ")
        logging.info(f" {verb} reading {n_fetch} lightcurves for " \
                     f"subgroup {i_subgroup:3d} of {n_subgroup} " \
                     f"(offset={n_offset})")

        if not match_split: return 0
        
        t0 = time.perf_counter()
        # here is the fastdb magic:
        manyltcvs = fdb.post( "ltcv/getmanyltcvs/realtime",                          
                              json={
                                  'limit': n_fetch,
                                  'offset' : n_offset,
                                  'include_source_positions':      True,
                                  'return_object_info':            True,
                                  'return_diaobject_positions':    True,
                                  'use_weighted_source_positions': True,
                                  'nonevalue': -999
                          } )

        # include_object_positions with return_diaobject_positions
        # xxx mark 'always_use_weighted_source_positions': True,

        if args.pdb:
            import pdb;  pdb.set_trace()   # so that Rob can take over my screen
        
        dict_objinfo = manyltcvs['objinfo']        
        dict_ltcvs   = manyltcvs['ltcvs']
        nevt         = len(dict_ltcvs);
    
        t1 = time.perf_counter()
        logging.info(f" Fetched {nevt} lightcurves in {t1-t0:.2f} sec.\n" )
        logging.info('')

        
        self.dict_objinfo = dict_objinfo
        self.dict_ltcvs   = dict_ltcvs
        self.ltcsv_offset = n_offset     # used only in GARBAGE table
        self.ltcsv_limit  = n_fetch      # used only in GARBAGE table
        self.garbage_table_extra_columns = { 'ltcvs_offset' : n_offset, 'ltcsv_limit': n_fetch }

        if i_subgroup == -9 :
            self.dump_dict_objinfo()
            
        return nevt

    # end prep_read_data_subgroup

    def dump_dict_objinfo(self):
        
        # dump structure of dict_objinfo to help debug this beast
        
        keys_objinfo = list(self.dict_objinfo.keys())

        logging.info("")
        logging.info(" dump_dict_objinfo: ")
        for key in keys_objinfo:
            ty = type(self.dict_objinfo[key])
            logging.info(f"\t dict_objinfo[{key}]  is type {ty}")

        logging.info("")            
        #sys.exit("\n xxx TEMP EXIT xxx")
        return
    
    def read_event(self, evt ):

        # Called by base code:
        # read event for event number "evt" in this subgroup.

        rootid          = self.get_rootid(evt)  # also checks rootid consistency        
        args            = self.config_inputs['args']
        
        # init output dictionaries
        dict_objinfo      = self.dict_objinfo 
        dict_ltcvs        = self.dict_ltcvs 
        lc_dict           = dict_ltcvs[evt]

        diaobjid   = dict_objinfo[FASTDB_KEYNAME_OBJID][evt][0]  # fragile; take first one on list
        SNID       = rootid
        nobs       = len(lc_dict['mjd'])
        self.nobs  = nobs
        
        lc_dict['photflag'] = [0] * nobs  # hack this in until we get more info from Rubin
        
        
        snana_head_raw, snana_head_calc, snana_head_sim = util.reset_data_event_dict()
        snana_head_raw[gpar.DATAKEY_SNID]        = str(SNID)
        snana_head_raw[gpar.DATAKEY_NAME_TRNS]   = diaobjid
        snana_head_raw[gpar.DATAKEY_SURVEY]   = SURVEY_LSST
        snana_head_raw[gpar.DATAKEY_FILTERS]  = FILTERLIST_LSST
        snana_head_raw[gpar.DATAKEY_FAKE]     = 0
        snana_head_raw[gpar.DATAKEY_SNTYPE]   = 0
        
        for key_snana in DATAKEY_HEAD_LIST:
            key_fastdb = KEYMAP[key_snana]
            if key_fastdb in dict_objinfo:                
                snana_head_raw[key_snana]  = dict_objinfo[key_fastdb][evt]

    
        nobs_before_coadd = nobs
            
        # - - - - - - --  -
        # load light curve info 
        snana_phot_raw          = self.init_phot_dict(nobs)
        snana_phot_raw['NOBS']  = nobs        

        self.set_garbage_list(evt, rootid, None, [])  # init garbage counter
        
        for key_snana in DATAKEY_PHOT_LIST:
            key_fastdb = KEYMAP.setdefault(key_snana,None)
            if key_fastdb in lc_dict :
                val_list   = lc_dict[key_fastdb]
                    
                if 'FLUXCAL' in key_snana:                    
                    self.set_garbage_list(evt, rootid, key_snana, val_list)
                    
                    if self.n_garbage_dict[gpar.GARBAGEKEY_FLUX_ALL] == 0:
                        val_list = [ x*FLUXSCALE_SNANA for x in val_list ]
                
                if key_snana == gpar.DATAKEY_PHOTFLAG:
                    val_list = [ args.photflag_detect * int(x) for x in lc_dict['isdet'] ]
                
                snana_phot_raw[key_snana] = val_list

        
        # tack on GARBAGE bit
        self.check_garbage(snana_head_raw, snana_phot_raw)
            
        # store first/second/last MJD_DETECT before coadd        
        self.store_mjd_detections(snana_head_calc, snana_phot_raw)
        
        # do nightly coadd on snana dictionary if (1) coadd is requested
        # as command line arg, and (2) there is no garbage
        # [coadding would hide the garbage]
        nobs_garbage = self.n_garbage_dict[gpar.GARBAGEKEY_FLUX_ALL]
        do_coadd = args.coadd_by_nite and nobs_garbage == 0

        nobs_after_coadd, snana_phot_coadd, nite_detect_dict = \
            self.coadd_by_nite(snana_phot_raw, do_coadd)
            
        # --- private LSST variables that are not standard SNANA vars
        ratio = nobs_after_coadd / nobs_before_coadd 
        lsst_head_private = {}
        lsst_head_private['PRIVATE(NOBS_BEFORE_COADD)']  =  nobs_before_coadd
        lsst_head_private['PRIVATE(NOBS_AFTER_COADD)']   =  nobs_after_coadd
        lsst_head_private['PRIVATE(RATIO_NOBS_COADD)']   =  ratio
        lsst_head_private['NOBS_GARBAGE']                =  nobs_garbage
        for b, n_nite in nite_detect_dict.items():
            key_private = f'PRIVATE(N_NITE_DETECT_{b})'
            lsst_head_private[key_private] = n_nite


        self.nobs_coadd_ratio.append(ratio)  # increment for ascii histogram

        # figure out DDF or WFD based on coadd/total NOBS ratio.
        # TO DO later: figure out DDF+WFD overlaps
        
        if ratio < CUT_NOBS_COADD_RATIO:
            field = FIELD_DDF
        else:
            field = FIELD_WFD
            
        snana_phot_coadd[gpar.DATAKEY_FIELD] = [field] * snana_phot_coadd[gpar.DATAKEY_NOBS]
        snana_head_raw[gpar.DATAKEY_FIELD]   = field
        
        # - - - - - -
        # construct SNANA dictionary            
        snana_data_dict = {
            'head_raw'     : snana_head_raw,
            'head_calc'    : snana_head_calc,
            'head_private' : lsst_head_private,
            'phot_raw'     : snana_phot_coadd            
        }

        #sys.exit(f"\n xxx BYE snana_data_dict = \n{snana_data_dict}")
        
        fix_now = False;
        if fix_now:
            # Apply event selection requirements
            select = self.select_event(snana_data_dict)
            if select:
                pass  # do NOT set 'select' to True
            else:
                snana_data_dict['select'] = False
                            
        return snana_data_dict
    
    # end read_event

    def get_rootid(self, evt):

        rootid_obj     = self.dict_objinfo[FASTDB_KEYNAME_ROOTID][evt]
        rootid_ltcvs   = self.dict_ltcvs[evt][FASTDB_KEYNAME_ROOTID] 

        if rootid_obj != rootid_ltcvs :
            msgerr = [ f"objinfo rootid={rootid_obj} but ",
                       f"ltcvs   rootid={rootid_ltcvs}",
                       f"evt={evt}  offset={self.ltcsv_offset}  limit={self.ltcsv_limit}" ]
            util.log_assert(False,msgerr)
            
        return rootid_obj
        
    # end get_rootid
    
    def fudge_garbabe_obs(self, evt, key_snana, val_list):

        # fudge garbage in observation to test garbage tracking infrastructure
        val_list_garbage = val_list
        match_key = (gpar.DATAKEY_FLUXCALERR == key_snana)
        match_evt = (evt % 23 == 0)
        if match_key and match_evt :
            val_list[0] = None 

        return val_list_garbage
        
    def check_garbage(self, snana_head_raw, snana_phot_raw):

        # if there is garbage, set garbage keys in snana_head_raw & snana_phot_raw

        args   = self.config_inputs['args']
        
        # check for garbage/missing coordinates
        RA  = snana_head_raw[gpar.DATAKEY_RA]
        DEC = snana_head_raw[gpar.DATAKEY_DEC]
        if not (isinstance(RA, float) and isinstance(DEC,float) ):
            snana_head_raw[gpar.GARBAGEKEY_RADEC] = True
        
        nobs_garbage = self.n_garbage_dict[gpar.GARBAGEKEY_FLUX_ALL] 
        snana_phot_raw[gpar.DATAKEY_NOBS_GARBAGE] = nobs_garbage

        if nobs_garbage > 0:
            photflag_before  = snana_phot_raw[gpar.DATAKEY_PHOTFLAG] 
            photflag_garbage = [ args.photflag_garbage * int(x) for x in self.garbage_list ]

            # set PHOTFLAG garbage bit for each obs that has garbage flux or fluxerr
            snana_phot_raw[gpar.DATAKEY_PHOTFLAG] = \
                [a+b for a,b in zip(photflag_before,photflag_garbage)]

            # for each garbage type, set header tag to print at top of each data file
            for key in list(self.n_garbage_dict.keys()):
                n = self.n_garbage_dict[key]  # number of obs with garbage flux
                if n > 0:
                    snana_head_raw[key]  = True

        # - - - - - - - -
        # if any garbage key is set, set the ALL flag
        for k in GARBAGEKEY_LIST:
            if snana_head_raw[k]:
                snana_head_raw[gpar.GARBAGEKEY_ALL] = True
                
        return
    
    def set_garbage_list(self, evt, rootid, key, val_list):
        
        # check for garbage in this val_list defined by
        #   + not a float
        #   + exact value is zero

        if len(val_list) == 0:
            # init garbage list and dictionary of counters
            self.garbage_list = [False] * self.nobs
            self.n_garbage_dict = {
                gpar.GARBAGEKEY_FLUX_ALL :     0,
                gpar.GARBAGEKEY_FLUX_ZERO:     0,
                gpar.GARBAGEKEY_FLUX_NOTFLOAT: 0
            }
            return

        # - - - - - - -
        FUDGE_GARBAGE_DEBUG = 0  # manual enable/disable
        if FUDGE_GARBAGE_DEBUG > 0 :
            val_list = self.fudge_garbabe_obs(evt, key, val_list)

        # - - - - - -
        old_garbage_list        = self.garbage_list
        new_n_garbage_dict, new_garbage_list = util.get_garbage_list(val_list)
        
        n_new = new_n_garbage_dict[gpar.GARBAGEKEY_FLUX_ALL]
        if n_new > 0:
            logging.info(f" WARNING: {n_new:3d} garbage {key} values for ROOTID={rootid}")
            self.garbage_list = [a or b for a, b in zip(old_garbage_list,new_garbage_list)]

            # increment garbage counter for each type of flt garbage
            for key in list(new_n_garbage_dict.keys()) :
                self.n_garbage_dict[key] += new_n_garbage_dict[key]
                
        return 0
    
    def coadd_by_nite(self, phot_raw, do_coadd):

        # for inpyut phot_raw dictionary of lists (table columns),
        # coadd each band grouped by nights and return coadd dictionary
        #
        # do_coadd = True -> return coadded photometry;
        #          = False -> return original photometr = other stuff
        #              computed from coadd; this is used for garbage data
        #            
        
        colnames = list(phot_raw)

        # drop columns with NOBS in the name
        nobs_keys = [key for key in phot_raw if 'NOBS' in key]
        for key in nobs_keys:
            colnames.remove(key)

        # convert to astropy table for easier manipulations
        phot_local_dict = {k: phot_raw[k] for k in colnames if k in phot_raw}
        t_phot = Table( phot_local_dict )

        t_coadd_list = []
        nite_detect_dict = {}
        for b in BAND_LIST_LSST:
            t_band     = t_phot[t_phot[gpar.DATAKEY_BAND] == b]  
            n_obs_band = len(t_band)
            nite_detect_dict[b] = 0
            if n_obs_band > 0:
                t_coadd_band = self.coadd_single_band(t_band, do_coadd)
                t_coadd_list.append(t_coadd_band)                

                photflag_list         = t_coadd_band[gpar.DATAKEY_PHOTFLAG].tolist()
                n_detect, detect_list = self.count_detect(photflag_list)                
                nite_detect_dict[b] = n_detect
                
        # combine bands into single table
        t_coadd = vstack(t_coadd_list)
        t_coadd.sort(gpar.DATAKEY_MJD)  # re-sort by MJD

        # covert coadd astropy table back to dictionary of lists 
        nobs_coadd = len(t_coadd)

        if do_coadd:
            phot_coadd_dict = { gpar.DATAKEY_NOBS :        nobs_coadd,
                                gpar.DATAKEY_NOBS_GARBAGE: 0 }
            for col in colnames:
                phot_coadd_dict[col] = t_coadd[col].tolist()
        else:
            phot_coadd_dict = phot_raw
            
        del phot_local_dict
        del t_coadd
        del t_phot
        
        return nobs_coadd, phot_coadd_dict, nite_detect_dict
    # end coadd_by_nite
        
    def coadd_single_band(self,t_band, do_coadd):

        # coadd table for this band, and also compute number of nites with detection
        
        from functools import reduce
        from operator import ior

        # sort by MJD
        t_band.sort(gpar.DATAKEY_MJD)

        # Identify where the difference exceeds tolerance
        diff           = np.diff( t_band[gpar.DATAKEY_MJD].data )
        new_group_mask = diff > TOL_MJD_NIGHT

        # 4. Generate group IDs using cumulative sum
        group_ids          = np.zeros(len(t_band), dtype=int)
        group_ids[1:]      = np.cumsum(new_group_mask)
        t_band['group_id'] = group_ids
    
        # 5. Group by new IDs
        grouped_table = t_band.group_by('group_id')        

        t_coadd_list = []
        colnames = t_band.colnames
        fscale   = FLUXSCALE_SNANA
        
        # Iterate through groups and take average
        for group in grouped_table.groups:
            t_coadd = Table()
            n_group = len(group)
            for col in colnames:
                val_list    = group[col].data  # list of column valoues over all rows

                if not do_coadd:
                    t_coadd[col] = None
                    continue ;
                
                if col == gpar.DATAKEY_MJD:
                    t_coadd[col]        = [ np.mean(val_list) ]
                
                elif col == gpar.DATAKEY_FLUXCAL:
                    # approximate assuming same ZP per exposure
                    t_coadd[col]        = [ fscale*np.mean(val_list) ]
                
                elif col == gpar.DATAKEY_FLUXCALERR:
                    errscale         = fscale/n_group
                    t_coadd[col]     = [ errscale*np.sqrt(np.sum(val_list**2)) ]

                elif col == gpar.DATAKEY_PHOTFLAG:
                    photflag = reduce(ior, val_list)
                    t_coadd[col] = photflag
                    
                elif 'group_id' not in col:
                    t_coadd[col]        = [ val_list[0] ]  # no change
                    
            t_coadd_list.append(t_coadd)

        # - - - - - - - 
        t_coadd = vstack(t_coadd_list)

        return t_coadd
    # end coadd_single_band
    
    def force_snr_detections(self, snana_data_dict):

        # if args.snr_detect is set, compute detection for each obs and set
        # args.photflag_detect mask of PHOTFLAG.
        #
        
        args            = self.config_inputs['args']
        snr_detect      = args.snr_detect
        photflag_detect = PHOTFLAG_DETECT
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

    def store_mjd_detections(self, head_calc, phot_raw):

        # compute and store MJD_DETECT_[FIRST/SECOND/LAST]
        args            = self.config_inputs['args']        
        photflag_detect = args.photflag_detect
        
        #head_calc    = snana_data_dict['head_calc']
        #phot_raw     = snana_data_dict['phot_raw']
        KEY_MJD      = gpar.DATAKEY_MJD
        KEY_PHOTFLAG = gpar.DATAKEY_PHOTFLAG

        if KEY_PHOTFLAG not in phot_raw : return

        mjd_detect_first  = -9.0
        mjd_detect_second = -9.0
        mjd_detect_last   = -9.0
        n_detect = 0
        
        j_first  = -9
        j_second = -9
        j_last   = -9
        
        # make list of detections per observation. This should work for
        # real detection and also for forced snr_detect.

        photflag_list = phot_raw[KEY_PHOTFLAG]
        if photflag_list[0] is None: return

        n_detect, detect_list = self.count_detect(photflag_list)

        #sys.exit(f"\n xxx photflag_list = {photflag_list} \n xxx detect_list = {detect_list}")
        
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
        
        #print(f" xxx mjd_detect = {mjd_detect_first}  {mjd_detect_second}  {mjd_detect_last} ")
        # store output

        head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]  = mjd_detect_first
        head_calc[gpar.DATAKEY_MJD_DETECT_SECOND] = mjd_detect_second        
        head_calc[gpar.DATAKEY_MJD_DETECT_LAST]   = mjd_detect_last        
        head_calc[gpar.DATAKEY_N_DETECT]          = n_detect
        
        return

    def count_detect(self,photflag_list):
        args          = self.config_inputs['args']        
        detect_list   = [ (int(x) & args.photflag_detect )>0 for x in photflag_list ]
        n_detect      = detect_list.count(True)
        return n_detect, detect_list
    
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

    # ----------------------------------------------
    def end_read_data_subgroup(self):
        pass
    
    def end_read_data(self):

        self.print_hist_coadd_ratio()
        # global end for reading data                           
        return

    def print_hist_coadd_ratio(self):
        # print ascii historgram of nobs_after_coadd / nobs_before_coadd
        # to check WFD/DDF separation


        ratio_list = self.nobs_coadd_ratio

        ratio_min = 0.0
        ratio_max = 1.0
        ratio_bin = 0.1
        bins = np.arange(ratio_min, ratio_max + ratio_bin, ratio_bin)
        binned_sums, _ = np.histogram(ratio_list, bins = bins )

        mx_bsum = max(binned_sums)
        nbin    = len(binned_sums)

        if mx_bsum == 0 : return

        logging.info("")
        logging.info(" Print distribution of coadd_ratio = nobs_after_coadd/nobs_before_coadd") 
        sym = '*'   # symbol for histogram
        mxsym = 70  # scale content axis so this max number of symbols displayed
        
        logging.info("")
        logging.info("   coadd_ratio   Nevt")
        for i in range(0,nbin):
            r0 = bins[i]
            r1 = bins[i+1]
            bsum       = binned_sums[i]
            nsym       = int(mxsym * (bsum/mx_bsum)+.5)
            string_sym = sym * nsym
            hline       = f"  {r0:3.1f} - {r1:3.1f} :  {string_sym} {bsum}"
            logging.info(f"{hline}")
                      
        logging.info("")
        
        #from  ascii_graph   import Pyasciigraph
        #graph = Pyasciigraph()
        #for line in graph.graph('Data Distribution', binnned_sum):
        #    print(line)
        
        return
    
    # end print_hist_coadd_ratio
    
    def exclude_varlist_obs(self):
        # return list of PHOT columns to excude from output text files
        # see VARNAMES_OBS in makeDataFiles_params.py
        return [ 'PSF_SIG1', 'ZEROPT',  'SKY_SIG', 'XPIX', 'YPIX', 'GAIN',
                 'CCDNUM', 'IMGNUM' ]
    

    
