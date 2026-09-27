# Created Sep 2026 by R.Kessler
#   [pulled out of read_fastdb_test.py to be independent uti]
#
# Generic import utility to pass photometry dictionary of lists,
# and return co-added phot dict, where co-add is nitely and in each band.
#
# Note that user-passed dictionary must include some required keys
# as indicated in the hard-coded KEYLIST_COADD_DICT lists below.
#

import os, sys,  logging
import numpy  as np
from  astropy.table import Table, vstack

TOL_MJD_ANY_BAND    = 0.6  # for any-band nite detection

# map local KEY_TYPE (left) to possible key names in passed phot dictionary (right)
KEYLIST_COADD_DICT = {
    'NOBS'      : [ 'nobs', 'nmjd', 'nep' ],
    'BAND'      : [ 'band', 'filter' ],
    'MJD'       : [ 'mjd' ],
    'DETECT'    : [ 'isdet', 'detect' ],
    'PHOTFLAG'  : [ 'photflag' ],
    'PHOTPROB'  : [ 'photprob', 'reliability', 'real_bogus' ],    
    'FLUX'      : [ 'flux', 'fluxcal' ],
    'FLUXERR'   : [ 'fluxerr', 'fluxcalerr' ],
    'DUMMY'     : [ 'dummy' ]
    }

# expand each list to includ upper-case equivalents
for key in list(KEYLIST_COADD_DICT):
    KEYLIST_COADD_DICT[key] += [s.upper() for s in KEYLIST_COADD_DICT[key] ]


# ====================================================
def coadd_by_nite(phot_dict, band_list_full, tol_mjd_coadd, do_coadd):

    # for input phot_dict dictionary of lists (table columns),
    # coadd each band grouped by nights and return coadd dictionary
    #
    # Inputs:
    #   phot_dict: dictionary of lists, mjd, band, flux, fluxerr, etc ...
    #              See KEYLIST_COADD_DICT for allowed keynames.
    #              Scalars in this dictionary (e.g., NOBS) are ignored.
    #
    # band_list_full = full list of bands for which nite_detect_dict is returned.
    #                  For missing bands in phot_dict, nite_detect_dict[band] = 0.
    #
    # tol_mjd_coadd: MJD tolerance for coadd; e.g, 0.2 -> coadd obs (same band) within 5.8hr
    #
    # do_coadd = True -> return coadded photometry;
    #          = False -> return original photometry;
    #            this option may be useful for data flagged as having garbage.
    #             
    # - - - - - - - 
    keynames  = list(phot_dict)

    # find needed keynames in phot_dict
    ABORT_FLAG  = 1 # -> abort on missing/required key
    key_nobs    = get_keyname_phot('NOBS',    keynames, ABORT_FLAG)
    key_band    = get_keyname_phot('BAND',    keynames, ABORT_FLAG)
    key_mjd     = get_keyname_phot('MJD',     keynames, ABORT_FLAG)
    key_detect  = get_keyname_phot('DETECT',  keynames, ABORT_FLAG)

    # get list of bands to coadd separately (below) in each band
    band_list = list(set(phot_dict[key_band]))

    DUMP_FLAG = 0
        
    # drop any key that is not a list since only a list can be coadded.
    # This avoids astropy confusion with scalar.
    keynames_removed = []
    keynames_for_table = []
    for key in keynames:
        if isinstance( phot_dict[key], list):
            keynames_for_table.append(key)
        else:
            keynames_removed.append(key)
            
            #keynames.remove(key)


    if DUMP_FLAG > 0:
        print(f" xxx coadd_by_nite DUMP : ")
        print(f" xxx keynames_orig = {keynames}")
        print(f" xxx keynames_for_table = {keynames_for_table}")
        print(f" xxx keynames_removed   = {keynames_removed}")                 
    
    # convert to astropy table for easier manipulations
    phot_local_dict = {k: phot_dict[k] for k in keynames_for_table if k in phot_dict}
    t_phot   = Table( phot_local_dict )
    colnames = keynames_for_table
    
    t_coadd_list = []

    # init all bands in nite_detect_dict, even those that don't exist in phot table
    nite_detect_dict = dict.fromkeys(band_list_full, 0)
        
    for b in band_list :
        t_band     = t_phot[t_phot[key_band] == b]  
        n_obs_band = len(t_band)
        if n_obs_band > 0:
            t_coadd_band        = coadd_single_band(t_band, colnames, tol_mjd_coadd, do_coadd)
            nite_detect_dict[b] = t_coadd_band[key_detect].sum()
            t_coadd_list.append(t_coadd_band) 

    # - - - - - -
    # combine bands into single table
    t_coadd = vstack(t_coadd_list)
    t_coadd.sort(key_mjd)       # re-sort by MJD
    
    # determine number of nites with detection, regardless of band : 'ANY_BAND'
    # need to coadd again using all bands together, and pass wide tolerance
    # to cover full nite.
    t_coadd_dummy = coadd_single_band(t_phot, colnames, TOL_MJD_ANY_BAND, do_coadd)  

    nite_detect_dict['ANY_BAND'] = t_coadd_dummy[key_detect].sum()
    
    # covert coadd astropy table back to dictionary of lists 
    nobs_coadd = len(t_coadd)

    # construct co-added (output) phot dictionary
    if do_coadd:
        phot_coadd_dict = { key_nobs: nobs_coadd }  # restore updated nobs key 
        for col in colnames:
            phot_coadd_dict[col] = t_coadd[col].tolist()
    else:
        # return original dictionary
        phot_coadd_dict = phot_dict

    # - - - - -
    if DUMP_FLAG > 0:
        sys.exit(f"\n xxx  nobs_coadd = {nobs_coadd}\n xxx t_coadd = \n{t_coadd}")


    del phot_local_dict
    del t_coadd
    del t_phot
    
    return nobs_coadd, phot_coadd_dict, nite_detect_dict
# end coadd_by_nite
        
def coadd_single_band(t_band, colnames, tol_mjd, do_coadd):

    # coadd table for this single band; return co-added table
    # Inputs:
    #   t_band:   astropy table with mjd, band, flux[err], etc ...
    #   colnames: list of column names to consider in coadd
    #   tol_mjd:  coadd within this time-window (days)
    #   do_coadd: do coadd if true; else return t_band unmodified
    #
    #
    
    from functools import reduce
    from operator import ior

    ABORT_FLAG  = 1  # abort on required key that is missing
    key_mjd        = get_keyname_phot('MJD',        colnames, ABORT_FLAG)
    key_flux       = get_keyname_phot('FLUX',       colnames, ABORT_FLAG)
    key_fluxerr    = get_keyname_phot('FLUXERR',    colnames, ABORT_FLAG)
    key_detect     = get_keyname_phot('DETECT',     colnames, ABORT_FLAG)
    
    ABORT_FLAG  = 0  # for optional key, return None on missing key
    key_photflag   = get_keyname_phot('PHOTFLAG',   colnames, ABORT_FLAG)
    key_photprob   = get_keyname_phot('PHOTPROB',   colnames, ABORT_FLAG)
    
    # sort by MJD
    t_band.sort(key_mjd)

    # Identify where the difference exceeds tolerance
    diff           = np.diff( t_band[key_mjd].data )
    new_group_mask = diff > tol_mjd

    # 4. Generate group IDs using cumulative sum
    group_ids          = np.zeros(len(t_band), dtype=int)
    group_ids[1:]      = np.cumsum(new_group_mask)
    t_band['group_id'] = group_ids
    
    # 5. Group by new IDs
    grouped_table = t_band.group_by('group_id')        

    t_coadd_list = []
    fscale   = 1.0
        
    # Iterate through groups and take average
    for group in grouped_table.groups:
        t_coadd = Table()
        n_group = len(group)
        for col in colnames:
            val_list    = group[col].data  # list of column valoues over all rows

            if not do_coadd:
                t_coadd[col] = None
                continue ;
                
            if col == key_mjd:
                t_coadd[col]        = [ np.mean(val_list) ]
                
            elif col == key_flux:
                # approximate assuming same ZP per exposure
                t_coadd[col]        = [ fscale*np.mean(val_list) ]
                
            elif col == key_fluxerr:
                errscale         = fscale/n_group
                t_coadd[col]     = [ errscale*np.sqrt(np.sum(val_list**2)) ]

            elif col == key_detect :
                detect = any(val_list)
                t_coadd[col] = detect

            elif col == key_photflag :
                photflag = reduce(ior, val_list)
                t_coadd[col] = photflag
                
            elif col == key_photprob:
                # average over observations; perhaps later can be wgted avg?
                t_coadd[col]        = [ fscale*np.mean(val_list) ]
                
            elif 'group_id' not in col:
                t_coadd[col]  = [ val_list[0] ]  # no change
                    
        t_coadd_list.append(t_coadd)

    # - - - - - - - 
    t_coadd = vstack(t_coadd_list)

    return t_coadd
# end coadd_single_band


def get_keyname_phot(KEY_TYPE, colnames, ABORT_FLAG):

    # for input KEY_TYPE and list of phot_dict colnames,
    # return relevant key in phot_dict.

    # ABORT_FLAG = 0 -> return None on missing key
    # ABORT_FLAG = 1 -> abort on missing key
    
    if KEY_TYPE not in KEYLIST_COADD_DICT:
        sys.exit(f"\n ERROR: Invalid KEY_TYPE = {KEY_TYPE} \n" \
                 f"\t Valid KEY_TYPES: {list(KEYLIST_COADD_DICT)} ")
        
    valid_keylist_phot = KEYLIST_COADD_DICT[KEY_TYPE]
    
    for key in valid_keylist_phot:
        if key in colnames:
            return key

    # if we get here, abort
    if ABORT_FLAG > 0:
        sys.exit(f"\n ERROR: could not find {KEY_TYPE}-type key in phot table\n" \
             f"\t Valid keys are {valid_keylist_phot}")
    else:
        return None
    # - - - - 



    
