import os, sys, glob, yaml, shutil, pickle, re
import requests
import json
import logging  # , coloredlogs
import numpy as np
import makeDataFiles_util  as    util
from   makeDataFiles_base    import Program
from   makeDataFiles_params  import *

# FILTERMAP_PKL2SNANA = {
#     # CSV ->     SNANA
#     'p48g'     : 'ZTF-g' ,
#     'p48r'     : 'ZTF-r' ,
#     'p48i'     : 'ZTF-i' ,
# }

# - - - - - - - - - - - - - - - - - - -     -
class data_lsst_tom_db(Program):
    def __init__(self, config_inputs) :
        config_data = {}
        print(" Init lsst data tom class.")
        super().__init__(config_inputs, config_data)

    def init_read_data(self):
        args = self.config_inputs['args']  # command line args
        match = re.search( "(.*):(.*)@(http.*)", args.lsst_tom_db )
        if match is None:
            errmsg = f"Failed to parse TOM connection string {args.lsst_tom_db}"
            util.log_assert( False, errmsg )
            raise ValueError( errmsg )
        user = match.group(1)
        password = match.group(2)
        self.url = match.group(3)

        self.rqs = requests.session()
        self.rqs.get( f'{self.url}/accounts/login' )
        res = self.rqs.post( f'{self.url}/accounts/login/',
                             data={ "username": user,
                                    "password": password,
                                    "csrfmiddlewaretoken": self.rqs.cookies['csrftoken'] } )
        if res.status_code != 200:
            raise RuntimeError( f"Failed to log in; http status: {res.status_code}" )
        if 'Please enter a correct' in res.text:
            # ...hoping that the parsing is still accurate...
            # import pdb; pdb.set_trace()
            raise RuntimeError( "Failed to log in.  I think.  Put in a debug break and look at res.text" )
        self.rqs.headers.update( { 'X-CSRFToken': self.rqs.cookies['csrftoken'] } )
        
        rows = self.send_query( 'SELECT COUNT("diaObjectId") AS count FROM elasticc_diaobject' )
        self.config_inputs['nevt'] = int( rows[0]['count'] )
        sys.stderr.write( f"Found {self.config_inputs['nevt']} events.\n" )
        
        self.cache_size = 100
        self.diaobjectcache = []
        self.diaobjectoffset = 0

        # end read_data_driver

    def send_query( self, query, subs={} ):
        res = self.rqs.post( f'{self.url}/db/runsqlquery/',
                             json={ 'query': query, 'subdict': subs } )
        if res.status_code != 200:
            raise RuntimeError( f'Got status {res.status_code} from {self.url}/db/runsqlquery/' )
        data = json.loads( res.text )
        if ( 'status' not in data ) or ( data['status'] != 'ok' ):
            if ( data['status'] == 'error' ) and ( 'error' in data ):
                raise RuntimeError( f"Error response from server: {data['error']}" )
            else:
                raise RuntimeError( f'Got unexpected response from {self.url}/db/runsqlquery/' )
        return data['rows']
        
    def prep_read_data_subgroup(self, i_subgroup):
        # No subgroups defined yet
        nevt = self.config_inputs['nevt']
        if i_subgroup == 0:
            return nevt
        else:
            return -1
        # end prep_read_data_subgroup
        
    def end_read_data_subgroup(self):
        pass
    def end_read_data(self):
        # global end for reading data
        pass

    def	exclude_varlist_obs(self):
        # return list of PHOT columns to excude from output text files
        return []
    
    def read_event(self, evt ):

        varlist_obs = self.config_data['varlist_obs']

        if ( evt < self.diaobjectoffset ) or ( evt >= self.diaobjectoffset + len(self.diaobjectcache) ):
            # Have to read a new bactch
            offset = self.cache_size * int( evt / self.cache_size )
            q = ( 'SELECT o."diaObjectId",o.ra,o."decl",o.mwebv,o.mwebv_err,t.zhelio,t.peakmjd '
                  'FROM elasticc_diaobject o '
                  'INNER JOIN elasticc_diaobjecttruth t ON o."diaObjectId"=t."diaObject_id" '
                  'ORDER BY o."diaObjectId" '
                  'LIMIT %(limit)s offset %(offset)s' )
            subs = { 'limit': self.cache_size,
                     'offset': offset }
            rows = self.send_query( q, subs )
            self.diaobjectoffset = offset
            self.diaobjectcache = rows

        diaobject = self.diaobjectcache[ evt - self.diaobjectoffset ]

        SNID = diaobject['diaObjectId']
        zHEL = diaobject['zhelio']
        zHEL_ERR = 0.001
        RA = diaobject['ra']
        DEC = diaobject['decl']
        PEAKMJD = diaobject['peakmjd']
        MW_EBV = diaobject['mwebv']
        MW_EBV_ERR = diaobject['mwebv_err']
        
        # init output dictionaries
        head_raw, head_calc, head_sim = util.reset_data_event_dict()
        head_raw[DATAKEY_SNID]        = str(SNID)    # I think this needs to be str?
        head_raw[DATAKEY_RA]          = RA
        head_raw[DATAKEY_DEC]         = DEC
        head_raw[DATAKEY_zHEL]        = zHEL
        head_raw[DATAKEY_zHEL_ERR]    = zHEL_ERR
        head_raw[DATAKEY_FIELD]       = FIELD_VOID

        # calc quantities
        head_calc[DATAKEY_PEAKMJD]          = int(PEAKMJD)
        head_calc[DATAKEY_MWEBV]            = MW_EBV
        head_calc[DATAKEY_MWEBV_ERR]        = MW_EBV_ERR
        head_calc[DATAKEY_MJD_DETECT_FIRST] = 1e32
        head_calc[DATAKEY_MJD_DETECT_LAST]  = -1e32
        
        # photometry
        # Photometry
        # Even though all the sources also show up in ForcedSources, we have to
        #  query both tables to figure out how to set PHOTFLAG

        rows = self.send_query( 'SELECT "diaSourceId","midPointTai","psFlux","psFluxErr","filterName" '
                                'FROM elasticc_diasource '
                                'WHERE "diaObject_id"=%(snid)s '
                                'ORDER BY "midPointTai"',
                                { 'snid': SNID } )
        sources = { row['diaSourceId']: row  for row in rows }

        for key in sources:
            sources[key]['PHOTFLAG'] = 4096
            if sources[key]['midPointTai'] < head_calc[DATAKEY_MJD_DETECT_FIRST]:
                head_calc[DATAKEY_MJD_DETECT_FIRST] = sources[key]['midPointTai']
            if sources[key]['midPointTai'] > head_calc[DATAKEY_MJD_DETECT_LAST]:
                head_calc[DATAKEY_MJD_DETECT_LAST] = sources[key]['midPointTai']

        sources[ list(sources.keys())[0] ]['PHOTFLAG'] = 6144
        
        rows = self.send_query( 'SELECT "diaForcedSourceId","midPointTai","psFlux","psFluxErr","filterName" '
                                'FROM elasticc_diaforcedsource '
                                'WHERE "diaObject_id"=%(snid)s '
                                'ORDER BY "midPointTai" ', { 'snid': SNID } )
        forcedsources = { row['diaForcedSourceId']: row  for row in rows }

        # I know from how I (rknop) write the alerts that if
        #  the diaForcedSourceId is the same as the diaSourceId,
        #  then these are the same observation
        forcedkeys = [ key for key in forcedsources if key not in sources ]
        for key in forcedkeys:
            sources[key] = forcedsources[key]
            sources[key]['PHOTFLAG'] = 0

        # Now that forced and unforced searches are merged, re-sort by midpointTai
        keys = list( sources.keys() )
        keys.sort( key=lambda key: sources[key]['midPointTai'] )

        phot_raw = self.init_phot_dict( len(keys) )
        phot_raw['NOBS'] = len(keys)
        phot_raw['MJD'] = np.array( [ sources[key]['midPointTai'] for key in keys ], dtype=np.float )
        phot_raw['ZPFLUX'] = np.array( [ 31.4 ] * len(keys), dtype=np.float )
        phot_raw['FIELD'] = [ 'WFD' ] * len(keys)
        phot_raw['BAND'] = [ sources[key]['filterName'][-1] for key in keys ]
        phot_raw['FLUXCAL'] = np.array( [ sources[key]['psFlux'] for key in keys ], dtype=np.float )
        phot_raw['FLUXCALERR'] = np.array( [ sources[key]['psFluxErr'] for key in keys ], dtype=np.float )
        phot_raw['PHOTFLAG'] = np.array( [ sources[key]['PHOTFLAG'] for key in keys ], dtype=np.int )

        # sys.stderr.write( f'Returning head_raw with {len(head_raw)} entries, '
        #                   f'head_calc with {len(head_calc)} entries, and'
        #                   f'phot_raw with {len(phot_raw)} entries; '
        #                   f'phot_raw["NOBS"]={phot_raw["NOBS"]} and '
        #                   f'len(phot_raw["MJD"]={len(phot_raw["MJD"])}\n' )
        # import pdb; pdb.set_trace()
        
        return {
            'head_raw': head_raw,
            'head_calc': head_calc,
            'phot_raw': phot_raw
        }
