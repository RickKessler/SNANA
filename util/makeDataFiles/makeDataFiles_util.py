# Generic Utilities for makeDataFiles.
#


import os, sys, yaml, shutil, glob, math
import logging, subprocess  # ,coloredlogs

import numpy as np
from   astropy.table import Table
from   makeDataFiles_params    import *

# =======================================
def select_subsample(args, var_dict):

    # return True of variables in var_dict pass selection
    # defined in args. Return False if any selection fails.
    # If there is a window, require xmin <= x < xmax.
    # The selection here is intened to be a CPU sub-sample 
    # rather than a traditional cut. 

    SNID       = var_dict[DATAKEY_SNID]    # int
    PEAKMJD    = var_dict[DATAKEY_PEAKMJD] # float
    MJD_DETECT = var_dict[DATAKEY_MJD_DETECT] # float

    nsplitran         = args.nsplitran  # from command line input
    isplitran_select  = args.isplitran

    # random split cut. Be careful that isplitran = 1 to N
    # (not 0 to N-1)
    if nsplitran > 1 :
        isplitran = (SNID % nsplitran) + 1
        if isplitran != isplitran_select:  return False

    # PEAKMJD
    if args.peakmjd_range :
        if PEAKMJD <  args.peakmjd_range[0]: return False
        if PEAKMJD >= args.peakmjd_range[1]: return False

    # MJD of first detection
    if args.mjd_detect_range :
        if MJD_DETECT <  args.mjd_detect_range[0]: return False
        if MJD_DETECT >= args.mjd_detect_range[1]: return False
    
    return True
    # end select_subsample

def init_readme_stats():
    readme_stats = {}
    for key in KEYLIST_README_STATS:   readme_stats[key] = 0
    return readme_stats

def write_readme(args, readme_dict, walltime = -1.0):

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
    line_list.append(f"USERNAME:         {USERNAME} ")
    line_list.append(f"HOSTNAME:         {HOSTNAME}")

    n_list = []
    for key in KEYLIST_README_STATS:
        key_plus_colon = f"{key}:"
        n = readme_stats[key]
        n_list.append(n)
        line_list.append(f"{key_plus_colon:<22}   {n}")

    if walltime > 0.0 :
        key_plus_colon = "WALLTIME:"
        line_list.append(f"{key_plus_colon:<22}   {walltime:.2f}   # seconds")

    nevt_all = n_list[0]
    line_list.append(f"ABORT_IF_ZERO:  {nevt_all}")

    # - - - - - - - -
    with open(readme_file,"wt") as f:
        if docana_flag:
            f.write(f"{DOCANA_KEY}: \n")
            indent_str = '  '

        for line in line_list:
            f.write(f"{indent_str}{line}\n")

        if docana_flag:
            f.write(f"{DOCANA_KEY_END}: \n")

    return
    # end write_readme

def write_yaml(file_name, yaml_contents):
    # write yaml_contents to file_name (e.g., for README)
    with open(file_name,"wt") as r:
        yaml.dump(yaml_contents, r, sort_keys=False)

def get_survey_snana(snana_folder):
    # for input snana_folder, run snana.exe GETINFO folder
    # and extract survey
    cmd = f"{PROGRAM_SNANA} GETINFO {snana_folder}"
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

# ---------------------------------------------------
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
