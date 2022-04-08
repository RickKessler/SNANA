# Generic Utilities for makeDataFiles.
#
# Mar 30 2022: add extract_sim_readme_info(...)

import os, sys, glob, math, yaml, tarfile
import logging, shutil, subprocess

import numpy as np
from astropy.io import fits

#from astropy.table import Table

#from makeDataFiles_params import *
import makeDataFiles_params as gpar

ASTROPLAN_EXISTS = False
try:
    from astroplan import Observer
    ASTROPLAN_EXISTS = True
except ImportError as e:
    pass
from astropy.time import Time

# =============================
def extract_sim_readme_info(path_sim,key_list):

    # Created Mar 31 2022
    # parse README file from path_sim and read yaml keys in key_list. 
    # Return dictionary of values corresponding to key_list.
    # Checks for README inside misc.tar.gz from submit_batch_jobs, 
    # and also checks for [GENVERSION].README from interactive job.

    value_dict = {}
    genversion = os.path.basename(path_sim)
    
    misc_tar_file = "misc.tar.gz"
    readme_file   = f"{genversion}.README"

    path_sim_expand = os.path.expandvars(path_sim)
    misc_tar_path = f"{path_sim_expand}/{misc_tar_file}"
    readme_path   = f"{path_sim_expand}/{readme_file}"

    # check of path_sim contains misc.tar (from submit_batch_jobs)
    # or README file (from interactive job)
    found_misc   = False
    found_readme = False
    if os.path.exists(misc_tar_path):
        found_misc = True
    elif os.path.exists(readme_path):
        found_readme = True
    else:
        msgerr = []
        msgerr.append(f"Could not find {misc_tar_file} nor {readme_file}")
        msgerr.append(f"Something not right in path_sim")
        msgerr.append(f"  {path_sim}")
        log_assert(False,msgerr)

    # - - - - - 
    if found_misc :
        # get list of files/members inside tar file
        tar = tarfile.open(misc_tar_path)
        misc_file_list   = tar.getnames()
        misc_member_list = tar.getmembers()

        # find first README file
        for fnam, member in zip(misc_file_list, misc_member_list):
            if "README" in fnam:
                readme_file   = fnam
                readme_member = member
                break

        # read README inside tar file
        f = tar.extractfile(readme_member)
    else:
        f = open(readme_path,"rt")

    # - - - -
    line_list = []
    for line_tmp in f:
        try:
            line = line_tmp.decode('utf-8')
        except:
            line = line_tmp

        line_list.append(line)

    # - - - - - - - - -
    if found_misc:
        tar.close()
    else:
        f.close()

    # - - - 
    readme_yaml = yaml.safe_load("\n".join(line_list))

    # get list of DOC yaml blocks to check since input key_list could
    # be under any yaml block
    DOC              = readme_yaml['DOCUMENTATION']
    DOC_BLOCK_LIST   = DOC.keys()  # e.g., OVERVIEW, INPUT_KEYS, INPUT_NOTES

    for key in key_list:
        value = None
        for block in DOC_BLOCK_LIST:
            if key in DOC[block]:
                value = DOC[block][key]
        value_dict[key] = value

    return value_dict
    # end extract_sim_readme_info


# =======================================
def select_subsample(args, var_dict):

    # return True of variables in var_dict pass selection
    # defined in args. Return False if any selection fails.
    # If there is a window, require xmin <= x < xmax.
    # The selection here is intened to be a CPU sub-sample
    # rather than a traditional cut.

    SNID       = var_dict[gpar.DATAKEY_SNID]    # int
    PEAKMJD    = var_dict[gpar.DATAKEY_PEAKMJD] # float
    MJD_DETECT_FIRST = var_dict[gpar.DATAKEY_MJD_DETECT_FIRST] # float

    nsplitran         = args.nsplitran  # from command line input
    isplitran_select  = args.isplitran

    # random split cut. Be careful that isplitran = 1 to N
    # (not 0 to N-1)
    if nsplitran > 1 and isplitran_select>=0 :
        isplitran = (SNID % nsplitran) + 1
        if isplitran != isplitran_select:  return False

    # PEAKMJD
    if args.peakmjd_range :
        if PEAKMJD <  args.peakmjd_range[0]: return False
        if PEAKMJD >= args.peakmjd_range[1]: return False

    # MJD of first detection
    if args.nite_detect_range :
        # make coarse approximate cut using MJD for speed
        if MJD_DETECT_FIRST <  args.nite_detect_range[0]-1.: return False
        if MJD_DETECT_FIRST >= args.nite_detect_range[1]+1.: return False

        # TODO - pass site into select_subsample
        NITE = get_sunset_mjd(MJD_DETECT_FIRST, site='CTIO')

        # make exact cut for MJD using NITE
        if NITE <  args.nite_detect_range[0]: return False
        if NITE >= args.nite_detect_range[1]: return False

    return True
    # end select_subsample

def get_sunset_mjd(mjd, site='CTIO'):
    '''
    Returns an MJD of sunset prior to input mjd as float - not a Time Object
    '''
    if ASTROPLAN_EXISTS:
        ctio = Observer.at_site(site)
        detect_time = Time(mjd, format='mjd')
        sun_set = ctio.sun_set_time(detect_time, which='previous').mjd
        NITE = sun_set
    else:
        NITE = mjd
    return NITE

def init_readme_stats():
    readme_stats = {}
    for key in gpar.KEYLIST_README_STATS:
        readme_stats[key] = 0
    return readme_stats

def write_readme(args, readme_dict, walltime=-1.0):

    # input args are the user-command line args.
    # input readme_dict is prepared by write_data_xxx module.

    readme_file  = readme_dict['readme_file']
    readme_stats = readme_dict['readme_stats']
    data_format  = readme_dict['data_format']
    docana_flag  = readme_dict['docana_flag']

    script_command = ' '.join(sys.argv)
    indent_str = ''

    # store line of yanl lines without indentation
    line_list = []
    line_list.append(f"PURPOSE:  transient lightcurve data files " \
                f"for analysis")

    if args.lsst_ap :
        line_list.append(f"SOURCE_LSST_AP:   {args.lsst_ap} ")

    if args.lsst_drp :
        line_list.append(f"SOURCE_LSST_DRP:  {args.lsst_drp} ")

    if args.sirah_folder is not None :
        line_list.append(f"SOURCE_SIRAH_FOLDER:  {args.sirah_folder}")

    if args.snana_folder is not None:
        line_list.append(f"SOURCE_SNANA_FOLDER:  {args.snana_folder}")

    line_list.append(f"SURVEY:           {args.survey}")
    line_list.append(f"FIELD:            {args.field} ")
    line_list.append(f"FORMAT:           {data_format} ")
    line_list.append(f"SCRIPT_COMMAND:   {script_command} ")
    line_list.append(f"USERNAME:         {gpar.USERNAME} ")
    line_list.append(f"HOSTNAME:         {gpar.HOSTNAME}")

    n_list = []
    for key in gpar.KEYLIST_README_STATS:
        key_plus_colon = f"{key}:"
        n = readme_stats[key]
        n_list.append(n)
        line_list.append(f"{key_plus_colon:<22}   {n}")

    if walltime > 0.0 :
        key_plus_colon = "WALLTIME:"
        line_list.append(f"{key_plus_colon:<22}   {walltime:.2f}   # seconds")

    nevt_all = n_list[0]
    if nevt_all == 0 : nevt_all=1  # never allow abort on 0 events.
    line_list.append(f"ABORT_IF_ZERO:  {nevt_all}")

    # - - - - - - - -
    with open(readme_file,"wt") as f:
        if docana_flag:
            f.write(f"{gpar.DOCANA_KEY}: \n")
            indent_str = '  '

        for line in line_list:
            f.write(f"{indent_str}{line}\n")

        if docana_flag:
            f.write(f"{gpar.DOCANA_KEY_END}: \n")

    return
    # end write_readme

def write_yaml(file_name, yaml_contents):
    # write yaml_contents to file_name (e.g., for README)
    with open(file_name,"wt") as r:
        yaml.dump(yaml_contents, r, sort_keys=False)

def get_survey_snana(snana_folder):
    # for input snana_folder, run snana.exe GETINFO folder
    # and extract survey

    cmd = f"{gpar.PROGRAM_SNANA} GETINFO {snana_folder}"
    ret = subprocess.run( [cmd], shell=True,
                          capture_output=True, text=True )
    ret_stdout = ret.stdout.split()

    key_survey = "SURVEY:"
    if key_survey not in ret_stdout:
        msgerr = []
        msgerr.append(f"Cannot find {key_survey} key from command")
        msgerr.append(f"   {cmd}")
        log_assert(False,msgerr)

    k      = ret_stdout.index(key_survey)
    survey = ret_stdout[k+1]
    return survey
    # end get_survey_snana

def read_yaml(yaml_file):
    yaml_lines = []
    with open(yaml_file,"rt") as y:
        for line in y: yaml_lines.append(line)

    contents  = yaml.safe_load("\n".join(yaml_lines))
    return contents
    # end read_yaml

def iyear_survey(survey, event_dict):
    iyear = -1

    if survey == 'LSST' :
        iyear = iyear_LSST(event_dict)
    elif survey == 'DES' :
        iyear = iyear_DES(event_dict)
    elif survey == 'SIRAH' :
        iyear = -1
    else:
        msgerr = []
        msgerr.append(f"Cannot determine iyear for survey={survey}")
        log_assert(False,msgerr)

    return iyear
    #end iyear_survey

def iyear_LSST(event_dict):

    # need to read season map to replace fixed MJD range.
    mjd   = event_dict['peakmjd']
    ra    = event_dict['ra']
    dec   = event_dict['dec']
    field = event_dict['field']  # DDF or WFD
    iyear = -1
    if mjd < 59800:
        iyear = 1
    elif mjd < 60100 :
        iyear = 2
    else:
        iyear = 3

    return iyear

def iyear_DES(event_dict):

    # DES year/season starts at 1.
    # SV is technically year 0, but there is no SV data in
    # final DES-SN data.

    mjd       = event_dict['peakmjd']
    mjd_start = 56500.0  # July 27 2013 -> start of Y1

    iyear = int((mjd - mjd_start)/365.0) + 1
    return iyear

# ==============================================
#  Utilities to transform coords and redshifts from
#  https://github.com/sncosmo/sndatasets/blob/master/sndatasets/utils.py
#
# ================================================

def hms_to_deg(h, m, s):
    return 15. * (h + m / 60. + s / 3600.)

def sxhr_to_deg(s):
    """sexagesimal hours to degrees"""
    h, m, s = s.split(':')
    return hms_to_deg(int(h), int(m), float(s))

def sdms_to_deg(sign, d, m, s):
    sign = 1. - 2. * (sign == '-')
    return sign * (d + m / 60. + s / 3600.)

def sx_to_deg(s):
    """sexagesimal to degrees. Sign must be the first character"""
    sign = s[0]
    d, m, s = s.split(':')
    return sdms_to_deg(sign, abs(int(d)), int(m), float(s))

def jd_to_mjd(t):
    return t - 2400000.5

def mag_to_flux(m, me, zp):
    """Convert magnitude and magnitude error to flux, given a zeropoint."""

    f = 10.**(0.4 * (zp - m))
    fe = math.log(10.) * 0.4 * me * f

    return f, fe

def radec_to_xyz(ra, dec):
    x = math.cos(np.deg2rad(dec)) * math.cos(np.deg2rad(ra))
    y = math.cos(np.deg2rad(dec)) * math.sin(np.deg2rad(ra))
    z = math.sin(np.deg2rad(dec))

    return np.array([x, y, z], dtype=np.float64)

def cmb_dz(ra, dec):
    """See http://arxiv.org/pdf/astro-ph/9609034
     CMBcoordsRA = 167.98750000 # J2000 Lineweaver
     CMBcoordsDEC = -7.22000000
    """

    # J2000 coords from NED
    CMB_DZ  = 371000. / 299792458.
    CMB_RA  = 168.01190437
    CMB_DEC = -6.98296811
    CMB_XYZ = radec_to_xyz(CMB_RA, CMB_DEC)

    coords_xyz = radec_to_xyz(ra, dec)

    dz = CMB_DZ * np.dot(CMB_XYZ, coords_xyz)

    return dz

def helio_to_cmb(z, ra, dec):
    """Convert from heliocentric redshift to CMB-frame redshift.

    Parameters
    ----------
    z : float
        Heliocentric redshift.
    ra, dec: float
        RA and Declination in degrees (J2000).
    """

    dz = -cmb_dz(ra, dec)
    one_plus_z_pec = math.sqrt((1. + dz) / (1. - dz))
    one_plus_z_CMB = (1. + z) / one_plus_z_pec

    return one_plus_z_CMB - 1.

def cmb_to_helio(z, ra, dec):
    """Convert from CMB-frame redshift to heliocentric redshift.

    Parameters
    ----------
    z : float
        CMB-frame redshift.
    ra, dec: float
        RA and Declination in degrees (J2000).
    """

    dz = -cmb_dz(ra, dec)
    one_plus_z_pec = math.sqrt((1. + dz) / (1. - dz))
    one_plus_z_helio = (1. + z) * one_plus_z_pec

    return one_plus_z_helio - 1.


# ========================================
#  Jan 10 2022:
#  These snana-reader utilities are outside the READ_SNANA_FOLDER class
#  so that legacy option works in same code version as --refac 110.
#  After legacy code is removed, these snana-reader utilities should be
#  moved inside the READ_SNANA_FOLDER class.
#
def open_fits(file_name):
    # check file_name and file_name.gz, and open the file that exists.
    # Function returns hdu pointer and number of rows in table.

    msgerr = []
    file_namegz = f"{file_name}.gz"
    if os.path.exists(file_namegz) :
        hdul = fits.open(file_namegz)
    elif os.path.exists(file_name):
        hdul = fits.open(file_name)
    else:
        msgerr.append(f"Cannot find fits file")
        msgerr.append(f" {file_name}   not")
        msgerr.append(f" {file_namegz} ")
        log_assert(False, msgerr)

    NROW = hdul[1].header['NAXIS2']
    return NROW, hdul

    # end open_fits

def reset_data_event_dict():

    # Util for reading SNANA FITS format.
    # reset all data values to -9 to ensure that every
    # key gets written to data files, even if read_event
    # code fails to set a value.
    # Jan 8 2022: move this function out of base.

    raw_dict  = {}
    calc_dict = {}
    sim_dict  = {}

    for key in gpar.DATAKEY_LIST_RAW :
        raw_dict[key] = -9

    for key in gpar.DATAKEY_LIST_CALC :
        calc_dict[key] = -9

    for key in gpar.DATAKEY_LIST_SIM :
        sim_dict[key] = -9

    return raw_dict, calc_dict, sim_dict
    # end reset_data_event_dict

def get_snana_table_value(varlist, irow, table):

    # Util for reading SNANA FITS format.
    # Read "irow" of SNANA FITS table for varlist, and return value.
    # Varlist = ['NAME1', 'NAME2', etc] is a list of column names to check.
    # This allows reading older legacy  names;
    # E.g., varlist = ['DEC', 'DECL'] returns value if either key is present.

    value = None
    for varname in varlist:
        try:
            value = table[varname][irow]
            return value
        except:
            pass  # just try next varname
    return value

    # end get_snana_table_value

def store_snana_hostgal(datakey_list, evt, table_dict, head_store):

    # store hostgal values in head_store dictionary.
    # Note that input head_store is modified here.
    #
    # Inputs:
    #  datakey_list:   list of keys to load
    #  evt:            event number/row number of table
    #  table_dict      fits table info
    #
    # Output:
    #   head_store:   output dictionary
    #
    # If HOSTGAL_SPECZ[PHOTOZ] < 0, set values to VAL_NULL for clarity.

    table_head = table_dict['table_head']
    head_names = table_dict['head_names']

#    len_base = len(gpar.HOSTKEY_BASE)
    for key in datakey_list:
        if gpar.HOSTKEY_BASE not in key:
            continue

        is_z = gpar.HOSTKEY_SPECZ in key or gpar.HOSTKEY_PHOTOZ in key
        key2 = key_hostgal_nbr(key,2)
        key3 = key_hostgal_nbr(key,3)
# xxx mark key2 = gpar.HOSTKEY_BASE + '2' + key[len_base:] # neighbor host
# xxx mark key3 = gpar.HOSTKEY_BASE + '3' + key[len_base:]
        key_list = [ key, key2, key3]
        for k in key_list:
            if k in head_names :
                val = table_head[k][evt]
                if is_z and val < 0.0 : val = gpar.VAL_NULL
                head_store[k] = val

    # end store_snana_hostgal

def key_hostgal_nbr(key,n):
    # if key = HOSTGAL_XXX and n=2, return HOSTGAL2_XXX

    len_base = len(gpar.HOSTKEY_BASE)
    BASE     = gpar.HOSTKEY_BASE
    base     = BASE.lower()
    if BASE in key:
        key_nbr  = BASE + str(n) + key[len_base:]
    elif base in key:
        key_nbr  = base + str(n) + key[len_base:]        
    else:
        key_nbr  = key
        
    return key_nbr
    # end key_hostgal_nbr

def store_snana_private(datakey_list, evt, table_dict):

    # store private values in head_store dictionary.
    #
    # Inputs:
    #  datakey_list:   list of keys to load
    #  evt:            event number/row number of table
    #  table_dict      fits table info
    #
    # function output:
    #   head_store:   output dictionary
    #

    table_head = table_dict['table_head']
    head_names = table_dict['head_names']
    head_private = {}

    for key in datakey_list:
        head_private[key] = table_head[key][evt]

    return head_private

    # end store_snana_private

def field_plasticc_hack(field, head_file_name):

    # ugly/embarassing hack to get field (DDF or WFD) from filename
    # because original plasticc data didn't store field.
    # If input field is NULL or VOID, set it based on head_file_name.

    missing_field = (field == gpar.FIELD_NULL or field == gpar.FIELD_VOID )
    if not missing_field : return field

    if gpar.FIELD_DDF in head_file_name:
        field = gpar.FIELD_DDF
    elif gpar.FIELD_WFD in head_file_name:
        field = gpar.FIELD_WFD
    else:
        msgerr.append(f"Unable to determine FIELD for")
        msgerr.append(f"{head_file_name}")
        log_assert(False,msgerr)

    return field
    # end field_plasticc_hack

# ====================================================
# utility class to read SNANA folder in FITS format
# ====================================================

class READ_SNANA_FOLDER:

    """
    Created Jan 8 2022 by R.Kessler
    Tools to read SNANA data files in FITS format.

    __init__ " read LIST file and store list of HEAD-FITS files

    Loop over ifile until nevt=0:
       nevt = exec_read(ifile):
           read and store contents of HEAD and PHOT for file index ifile

       data_dict = get_data_dict(args, evt):
         Return data_dict for input args(cuts) and event/row 'evt'.
         This is the standard data_dict that is written to output.

       end_read(): close HEAD and PHOT fits file

    """

    def __init__(self, data_folder):

        # Created Jan 8 2022 by R.Kessler

        version      = os.path.basename(data_folder)
        list_file    = f"{data_folder}/{version}.LIST"

        # scoop up contents of LIST file
        with open(list_file, 'r') as f:
            HEAD_file_list = f.read().replace('\n', ' ').split()

        n_HEAD_file = len(HEAD_file_list)

        logging.info(f" Read data version = {version}")
        logging.info(f" from data_folder = {data_folder}")
        logging.info(f" Found {n_HEAD_file} FITS-HEAD files.")

        self.snana_folder_dict  = {}
        self.snana_folder_dict['data_folder']     = data_folder
        self.snana_folder_dict['version']         = version
        self.snana_folder_dict['HEAD_file_list']  = HEAD_file_list
        self.snana_folder_dict['n_HEAD_file']     = n_HEAD_file
        self.snana_folder_dict['private_dict']    = {}
        self.snana_folder_dict['zphot_q_dict']    = {} # photo-z quantiles


    def exec_read(self, ifile):

        n_HEAD_file      = self.snana_folder_dict['n_HEAD_file']
        if ifile >= n_HEAD_file  :
            return 0 # done reading

        data_folder      = self.snana_folder_dict['data_folder']
        HEAD_file_base   = self.snana_folder_dict['HEAD_file_list'][ifile]

        HEAD_file       = f"{data_folder}/{HEAD_file_base}"
        nevt, hdu_head  = open_fits(HEAD_file)

        PHOT_file_base  = hdu_head[0].header['PHOTFILE']
        PHOT_file       = f"{data_folder}/{PHOT_file_base}"
        NROW, hdu_phot  = open_fits(PHOT_file)

        logging.info(f"   Read {nevt} events from {HEAD_file_base}")

        table_head = hdu_head[1].data
        table_phot = hdu_phot[1].data

        head_names = table_head.columns.names
        phot_names = table_phot.columns.names

        table_dict = {
            'head_file'  : HEAD_file_base,
            'table_head' : table_head,
            'table_phot' : table_phot,
            'head_names' : head_names,
            'phot_names' : phot_names
        }

        self.snana_folder_dict['nevt_subgroup'] = nevt
        self.snana_folder_dict['table_dict'] = table_dict
        self.snana_folder_dict['hdu_head']   = hdu_head
        self.snana_folder_dict['hdu_phot']   = hdu_phot

        if ifile == 0:
            self.init_first_file()

        return nevt
        # end exec_read

    def init_first_file(self):

        # init a few things after opening first fits file

        table_dict   = self.snana_folder_dict['table_dict']
        private_dict = self.snana_folder_dict['private_dict']

        head_names = table_dict['head_names']
        phot_names = table_dict['phot_names']

        # check for true mag in PHOT table
        # e.g., fakes overlaid on images or sim
        if gpar.VARNAME_TRUEMAG in phot_names:
            gpar.VAL_UNDEFINED_LIST   += [gpar.VAL_NULL]
            gpar.VARNAMES_FMT_LIST    += ["8.4f"]
            gpar.VARNAMES_OBS_LIST    += [gpar.VARNAME_TRUEMAG]

        # if there are PRIVATE variables, add them to list to read/store
        # If no private variables are specified in private_dict,
        # then keep them all. Note that key = 'PRIVATE(subkey)'
        n_private_keep = len(private_dict)
        for key in head_names:
            if key[0:7] != "PRIVATE" : continue
            keep = True
            if n_private_keep > 0 :
                subkey = key.split('(')[1][:-1]
                keep   = subkey in private_dict
                #print(f" xxx key={key}  subkey={subkey}  keep={keep}")
            if keep :
                gpar.DATAKEY_LIST_PRIVATE  += [ key ]

        # Feb 11 2022 check for zphot quantiles in global header
        #    avoid catching NZPHOT_Q
        for key in head_names:
            if '_ZPHOT_Q' in key:   # fragile alert
                gpar.DATAKEY_LIST_ZPHOT_Q  += [ key ]

        if len(gpar.DATAKEY_LIST_ZPHOT_Q) > 0:
            gpar.DATAKEY_LIST_CALC += gpar.DATAKEY_LIST_ZPHOT_Q
        
        # end init_first_file

    def init_private_dict(self, private_dict):
        # store dictionary of private variables to keep
        # If this function isn't called, then all private variables
        # are kept by default.
        self.snana_folder_dict['private_dict']   = private_dict
        # end init_private_dict

    def get_data_val(self, key, evt):
        # return value for key and evt=row
        table_head = self.snana_folder_dict['table_dict']['table_head']
        return table_head[key][evt]

    def get_data_dict(self, args, evt):

        # Inputs:
        #   args:  input command line args to implement user cuts
        #   evt:   event/row number to get data_dict
        #
        # Outut:
        #   Return data_dict for this 'evt'

        msgerr     = []

        table_dict = self.snana_folder_dict['table_dict']

        # define local pointers to head and phot tables from FITS file
        table_head = table_dict['table_head']
        table_phot = table_dict['table_phot']
        head_names = table_dict['head_names']
        phot_names = table_dict['phot_names']

        # init output dictionaries
        head_raw, head_calc, head_sim = reset_data_event_dict()

        # read and store SNID
        try:
            SNID = table_head.SNID[evt].decode('utf-8').replace(' ','')
        except:
            SNID = table_head.SNID[evt]
        head_raw[gpar.DATAKEY_SNID]    = SNID

        # read and store SNTYPE (typically this is spectroscopic type or 0)
        head_raw[gpar.DATAKEY_SNTYPE]  = table_head.SNTYPE[evt]

        # read and store coords; allow DEC or legacy DECL name
        head_raw[gpar.DATAKEY_RA]    = table_head.RA[evt]
        head_raw[gpar.DATAKEY_DEC] = \
            get_snana_table_value(['DEC','DECL'],evt,table_head)


        # lightcurve-MJD info. Note that MJD_DETECT_FIRST is optional
        head_calc[gpar.DATAKEY_PEAKMJD]   = table_head.PEAKMJD[evt]

        KEY0 = gpar.DATAKEY_MJD_DETECT_FIRST
        KEY1 = gpar.DATAKEY_MJD_DETECT_LAST
        if gpar.DATAKEY_MJD_DETECT_FIRST in head_names:
            head_calc[KEY0] =  table_head.MJD_DETECT_FIRST[evt]
            head_calc[KEY1] =  table_head.MJD_DETECT_LAST[evt]
        else:
            # data does not have [MJD_DETECT_FIRST,LAST], so abort if
            # if nite-range selection is requsted
            nite_range = args.nite_detect_range
            if nite_range is not None:
                msgerr.append(f"Cannot implement nite_detect_range={nite_range}")
                msgerr.append(f"because {KEY0} is not in data header")
                log_assert(False,msgerr)

        # - - - - - - -
        # check user sub-sample selection here to avoid reading
        # remainder of header and photometry for rejected events.
        apply_select = True
        if apply_select :
            var_dict = {
                gpar.DATAKEY_SNID       : int(SNID),
                gpar.DATAKEY_RA         : head_raw[gpar.DATAKEY_RA],
                gpar.DATAKEY_DEC        : head_raw[gpar.DATAKEY_DEC],
                gpar.DATAKEY_PEAKMJD    : head_calc[gpar.DATAKEY_PEAKMJD],
                gpar.DATAKEY_MJD_DETECT_FIRST : \
                head_calc[gpar.DATAKEY_MJD_DETECT_FIRST]
            }
            sel = select_subsample(args,var_dict)
            # If this event is not selected, just return and skip reading phot
            if not sel:
                data_dict = {
                    'head_raw'  : head_raw,
                    'head_calc' : head_calc,
                    'select'    : False
                }
                return data_dict

        # - - - - - - -
        # load helio redshift info
        head_raw[gpar.DATAKEY_zHEL]      = table_head.REDSHIFT_HELIO[evt]
        head_raw[gpar.DATAKEY_zHEL_ERR]  = table_head.REDSHIFT_HELIO_ERR[evt]

        # strip off calculated zCMB redshfit
        head_calc[gpar.DATAKEY_zCMB]       = table_head.REDSHIFT_FINAL[evt]
        head_calc[gpar.DATAKEY_zCMB_ERR]   = table_head.REDSHIFT_FINAL_ERR[evt]


        # Galactic extinction
        head_calc[gpar.DATAKEY_MWEBV]      = table_head.MWEBV[evt]
        head_calc[gpar.DATAKEY_MWEBV_ERR]  = table_head.MWEBV_ERR[evt]

        # - - - - - -
        # store HOSTGAL and HOSTGAL2 keys in head_raw[calc].
        # head_[raw,calc] is both input and output of store_snana_hostgal
        store_snana_hostgal(gpar.DATAKEY_LIST_RAW,  evt, table_dict, head_raw)
        store_snana_hostgal(gpar.DATAKEY_LIST_CALC, evt, table_dict, head_calc)

        # - - - - -
        # store optional PRIVATE variables
        head_private = \
            store_snana_private(gpar.DATAKEY_LIST_PRIVATE, evt, table_dict)
        
        # check for true sim type (sim or fakes), Nov 14 2021
        KEY = gpar.SIMKEY_TYPE_INDEX
        if KEY in head_names:   head_sim[KEY] = table_head[KEY][evt]

        # - - - - - - - - - - -
        # get pointers to PHOT table.
        # Beware that PTROBS pointers start at 1 instead of 0,
        # so subtract 1 here to have python indexing.
        ROWMIN = table_head.PTROBS_MIN[evt] - 1
        ROWMAX = table_head.PTROBS_MAX[evt] - 1
        NOBS   = ROWMAX - ROWMIN + 1

        phot_raw   = {}
        phot_raw['NOBS'] = NOBS

        table_column_names = table_phot.columns.names

        # check for reading legacy FLT column name (instead of BAND).
        # Note that data_dict loads 'BAND' even if the input column is FLT.
        LEGACY_FLT = False
        if 'FLT' in table_column_names:  LEGACY_FLT = True

        for varname in gpar.VARNAMES_OBS_LIST:
            phot_raw[varname] = [ None ] * NOBS
            varname_table = varname
            if LEGACY_FLT:
                if varname == 'BAND' : varname_table = 'FLT'

            if varname_table in table_column_names :
                phot_raw[varname] = \
                    table_phot[varname_table][ROWMIN:ROWMAX+1].copy()

        # - - - - -
        # get field from from first observation,
        # Beware that light curve can overlap multiple fields.
        field = phot_raw[gpar.DATAKEY_FIELD][0]
        if args.survey == 'LSST' :
            field = field_plasticc_hack(field, table_dict['head_file'])
        head_raw[gpar.DATAKEY_FIELD] = field

        # - - - -
        # load blank dictionary for spectra ... to be filled later.
        spec_raw = {}

        # - - - - -
        # load output dictionary
        data_dict = {
            'head_raw'     : head_raw,
            'head_calc'    : head_calc,
            'head_private' : head_private,
            'phot_raw'     : phot_raw,
            'spec_raw'     : spec_raw,
            'select'       : True
        }

        # check optional dictionary items to append
        if len(head_sim) > 0:
            data_dict['head_sim'] = head_sim

        return data_dict

        # end get_data_dict

    def end_read(self):
        self.snana_folder_dict['hdu_head'].close()
        self.snana_folder_dict['hdu_phot'].close()


# -----------------------
# LOG MESSAGING
# ------------------------
class MessageStoreLogger(logging.Handler):
    store = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.store = {}

    def emit(self,record):
        l  =  record.levelname
        if l not in self.store:
            self.store[l] = []
        self.store[l].append(record)

    def get_warnings(self):
        return self.store.get("WARNING", []) + []

    def get_errors(self):
        return self.store.get("CRITICAL", []) + self.store.get("ERROR", []) + \
                self.store.get("EXCEPTION", []) + []

    def print_warnings(self):
        items = self.get_warnings()
        if not items:
            logging.notice("No warnings")
        else:
            logging.warning(f"{len(items)} Warnings:")
            for w in items:
                logging.warning(f"\t{w.msg}")

    def print_errors(self):
        items = self.get_errors()
        if not items:
            logging.notice("No errors")
        else:
            logging.error(f"{len(items)} Errors:")
            for w in items:
                logging.error(f"\t{w.msg}")


def setup_logging(args):
    level = logging.DEBUG if args.verbose else logging.INFO

    # Adding notice level for more important lines that arent warning
    message_store = MessageStoreLogger()
    message_store.setLevel(logging.WARNING)
    notice_level = 25
    logging.addLevelName(notice_level, "NOTICE")

    def notice(self, message, *args, **kws):
        if self.isEnabledFor(notice_level):
            self._log(notice_level, message, args, **kws)

    def notice_root(message, *args, **kwargs):
        logging.log(notice_level, message, *args, **kwargs)

    logging.Logger.notice = notice
    logging.notice = notice_root
    fmt = "[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s" \
          if args.verbose else "%(message)s"
    handlers = [logging.StreamHandler(), message_store]
    handlers[0].setLevel(level)
    logging.basicConfig(level=level, format=fmt, handlers=handlers)

    # coloredlogs.install(level=level, fmt=fmt, reconfigure=True,
    #                     level_styles=coloredlogs.parse_encoded_styles(
    # "debug=8#;notice=green;warning=yellow;error=red,bold;critical=red,inverse"
    # ),)
    return message_store


def log_assert(condition, message):
    if not condition:
        for item in message :
            gpar.ABORT_FACE_MSSG += f"\n   {item}"

        logging.exception(message)
        assert condition, gpar.ABORT_FACE_MSSG
        # end log_assert
