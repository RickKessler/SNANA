# Generic Utilities for makeDataFiles.
#

import glob
import logging  # ,coloredlogs
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import yaml
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

# - - - - -
def iyear_DES(event_dict):
    return -1

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
        log_assert(False,msgerr)

    NROW = hdul[1].header['NAXIS2']
    return NROW, hdul

    # end open_fits

def reset_data_event_dict():

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

# -----------------------------------
#   read SNANA folder ?? .xyz
# -----------------------------------
#  def init_read_snana_folder(folder, snana_folder_dict) 
#      store goodies in snana_folder_dict
#      return 
#
#  def exec_read_snana_folder(ifile,snana_folder_dict)
#      scoop up contents of HEAD and PHOT for ifile
#      perform actions of prep_read_data_subgroup
#
#  def end_read_snana_folder(snana_folder_dict)
#     close hdus
#
#  def open_snana_fits
#      return hdul
#
#  def get_event_snana_folder(evt, snana_folder_dict)
#     return data_dict for evt

class READ_SNANA_FOLDER:
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

        # on first subgroup, check for true mag in PHOT table
        # e.g., fakes overlaid on images or sim
        if ifile == 0 and gpar.VARNAME_TRUEMAG in phot_names:
            gpar.VAL_UNDEFINED_LIST   += [gpar.VAL_NULL]
            gpar.VARNAMES_FMT_LIST    += ["8.4f"]
            gpar.VARNAMES_OBS_LIST    += [gpar.VARNAME_TRUEMAG]

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

        return nevt
        # end exec_read

    def get_data_dict(self):
        msgerr     = []
        table_dict = self.snana_folder_dict['table_dict']

        # define local pointers to head and phot tables from FITS file
        table_head = table_dict['table_head']
        table_phot = table_dict['table_phot']
        head_names = table_dict['head_names']
        phot_names = table_dict['phot_names']

        # end get_data_dict

    def end_read(self):
        self.snana_folder_dict['hdu_head'].close()
        self.snana_folder_dict['hdu_phot'].close()


# -----------------------
# MESSAGING
# ------------------------
class MessageStore(logging.Handler):
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
        return self.store.get("CRITICAL", []) + self.store.get("ERROR", []) + self.store.get("EXCEPTION", []) + []

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
    message_store = MessageStore()
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

#    coloredlogs.install(level=level, fmt=fmt, reconfigure=True,
#                        level_styles=coloredlogs.parse_encoded_styles("debug=8#;notice=green;warning=yellow;error=red,bold;critical=red,inverse"),)
    return message_store


def log_assert(condition, message):
    if not condition:

        msg_abort_face = (
            f"\n\n"
            f"\n   `|```````|`    "
            f"\n   <| o\\ /o |>   "
            f"\n    | ' ; ' |     "
            f"\n    |  ___  |     ABORT makeDataFiles on Fatal Error. "
            f"\n    | |' '| |     "
            f"\n    | `---' |     "
            f"\n    \\_______/    "
            f"\n"
            f"\nFATAL ERROR ABORT : "
        )

        for item in message :
            msg_abort_face += f"\n   {item}"

        logging.exception(message)
        assert condition, msg_abort_face
        # end log_assert
