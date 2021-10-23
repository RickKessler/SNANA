# ==================================
# Created July 2021 by R.Kessler
#
# Define program constants here ...
# maybe later these can be read from config file ?
# ================================

import os
import datetime
import time
import getpass

DOCANA_KEY     = "DOCUMENTATION"
DOCANA_KEY_END = "DOCUMENTATION_END"

FORMAT_TEXT = "TEXT"
FORMAT_FITS = "FITS"

USERNAME = getpass.getuser()
HOSTNAME = os.environ['HOSTNAME']

MXSEASON = 12 # max number of seasons

# - - - - - - - 
# define survey info that is fixed for all events;
# e.g., NXPIX = SURVEY_INFO['CCD'][survey][0]
SURVEY_INFO = {
    'FILTERS' : {           # mandatory 
        'LSST'  : "ugrizY",
        'SIRAH' : "GRcogrizyABCLYJNH",  # ztf/ATLAS/PS1/WFC3
        'DES'   : "griz" ,
        'PS1'   : "griz"
    },
    'CCD' : {   # optional
        'LSST' : [ 4072, 4000, 0.199598], # NXPIX, NYPIX, pixsize
        'DES'  : [ 2048, 4096, 0.263   ]
    }
} # end SURVEY_INFO


SNANA_FLAG_DATA = 0
SNANA_FLAG_FAKE = 1
SNANA_FLAG_SIM  = 2

SNANA_ZP = 27.5 # FLUXCAL = FLUX * 10^[ -0.4*(ZP-SNANA_ZP) ]
VPEC_DEFAULT = [ 0.0, 300.0 ]  # VPEC and error, km/sec

# define list of variable names for each observation;
VARNAMES_OBS = "MJD BAND FIELD PHOTFLAG  " \
               "XPIX YPIX CCDNUM GAIN " \
               "FLUXCAL FLUXCALERR ZPFLUX NEA SKYSIG"

# define text format for each VARNAMES_OBS
VARNAMES_FMT = "10.4f 2s   8s    4d      6.1f 6.1f  4d     6.3f "\
               "12.4e    12.4e    8.4f  6.3f 6.2f"

# define values for undefined variables ...
# value set to VAL_ABORT will trigger abort because it is required.
VAL_ABORT    = 666  
VAL_NULL     = -9
VAL_UNDEFINED_LIST = [
    VAL_ABORT, VAL_ABORT, "VOID",  0,      # for MJD BAND FIELD PHOTFLAG
    VAL_NULL,  VAL_NULL,  VAL_NULL, VAL_NULL,
    VAL_ABORT, VAL_ABORT, VAL_NULL, VAL_NULL, VAL_NULL ]

VARNAME_TRUEMAG = "SIM_MAGOBS"  # for fakes, add this to VARNAMES_OBS

FIELD_DDF      = "DDF"
FIELD_WFD      = "WFD"
FIELD_HST      = "HST"
FIELD_DEEP     = "DEEP"
FIELD_MEDIUM   = "MEDIUM"
FIELD_SHALLOW  = "SHALLOW"
FIELD_VOID     = "VOID"

PREFIX_SEASON   = "Y"       # e.g.. for Y03 in file name
PREFIX_SPLIT    = "SPLIT"   # e.g., for SPLIT012 in file name
TEXTFILE_SUFFIX = ".DAT"    # used for intermediate TEXT  file name

# snana program to convert TEXT to FITS
PROGRAM_SNANA = "snana.exe"

# set TEXT->FITS options to
#  + define MWEBV from SFD98
#  + estimate PEAKMJD from fmax-clump method

OPTIONS_TEXT2FITS_SNANA = \
        "OPT_YAML 1  OPT_MWEBV 2  OPT_SETPKMJD 20"
OPTION_TEXT2FITS_SPECTRA_SNANA =  \
        "OPT_REFORMAT_FITS 128"

# for writing events, update screen after this many
NEVT_SCREEN_UPDATE = 500

# define yaml keys to store statistics for README
KEYLIST_README_STATS = [ 'NEVT_ALL', 
                         'NEVT_HOSTGAL_SPECZ', 'NEVT_HOSTGAL_PHOTOZ',
                         'NEVT_SPECTRA' ]
  
# define key names for data_event_dict dictionary,
# and for TEXT-formatted data files (readable by snana codes)
# The global list is used to initial all values to -9

DATAKEY_SURVEY      = "SURVEY"
DATAKEY_FILTERS     = "FILTERS"
DATAKEY_NXPIX       = "NXPIX"
DATAKEY_NYPIX       = "NYPIX"
DATAKEY_PIXSIZE     = "PIXSIZE"

DATAKEY_SNID        = "SNID"
DATAKEY_FAKE        = "FAKE"
DATAKEY_RA          = "RA"
DATAKEY_DEC         = "DEC"
DATAKEY_zHEL        = "REDSHIFT_HELIO"
DATAKEY_zHEL_ERR    = "REDSHIFT_HELIO_ERR"
DATAKEY_zCMB        = "REDSHIFT_CMB"
DATAKEY_zCMB_ERR    = "REDSHIFT_CMB_ERR"
DATAKEY_VPEC        = "VPEC"
DATAKEY_VPEC_ERR    = "VPEC_ERR"

DATAKEY_MWEBV       = "MWEBV"
DATAKEY_MWEBV_ERR   = "MWEBV_ERR"
DATAKEY_PEAKMJD     = "PEAKMJD"
DATAKEY_MJD_DETECT  = "MJD_DETECT_FIRST"
DATAKEY_FIELD       = "FIELD"
DATAKEY_NOBS        = "NOBS"

HOSTKEY_OBJID         = "HOSTGAL_OBJID"
HOSTKEY_PHOTOZ        = "HOSTGAL_PHOTOZ"
HOSTKEY_PHOTOZ_ERR    = "HOSTGAL_PHOTOZ_ERR"
HOSTKEY_SPECZ         = "HOSTGAL_SPECZ"
HOSTKEY_SPECZ_ERR     = "HOSTGAL_SPECZ_ERR"
HOSTKEY_SNSEP         = "HOSTGAL_SNSEP"
HOSTKEY_LOGMASS       = "HOSTGAL_LOGMASS"
HOSTKEY_ELLIP         = "HOSTGAL_ELLIPTICITY"
HOSTKEY_SQRADIUS      = "HOSTGAL_SQRADIUS"
HOSTKEY_PREFIX_MAG    = "HOSTGAL_MAG"         # band-dependent
HOSTKEY_PREFIX_MAGERR = "HOSTGAL_MAGERR"      # idem
HOSTKEY_PREFIX_SB     = "HOSTGAL_SB_FLUXCAL"  # idem

HOSTKEY_PREFIX_LIST = [HOSTKEY_PREFIX_MAG, HOSTKEY_PREFIX_MAGERR, 
                       HOSTKEY_PREFIX_SB ]

DATAKEY_LIST_RAW = \
    [ DATAKEY_SURVEY, DATAKEY_SNID, DATAKEY_FAKE, DATAKEY_FILTERS,
      DATAKEY_NXPIX, DATAKEY_NYPIX, DATAKEY_PIXSIZE,
      DATAKEY_RA, DATAKEY_DEC,
      DATAKEY_zHEL, DATAKEY_zHEL_ERR, DATAKEY_FIELD,
      HOSTKEY_OBJID, HOSTKEY_SPECZ, HOSTKEY_SPECZ_ERR, HOSTKEY_SNSEP,
      HOSTKEY_ELLIP, HOSTKEY_SQRADIUS
    ]

DATAKEY_LIST_CALC = \
    [ DATAKEY_zCMB, DATAKEY_zCMB_ERR, DATAKEY_MWEBV, DATAKEY_MWEBV_ERR, 
      DATAKEY_PEAKMJD, DATAKEY_MJD_DETECT, 
      HOSTKEY_PHOTOZ, HOSTKEY_PHOTOZ_ERR, HOSTKEY_LOGMASS
    ]

# -------
MODE_MERGE_MOVE = "MERGE_MOVE" # move files, then remove original folder
MODE_MERGE_LINK = "MERGE_LINK" # merge with sym links; keep orig folder

# === END ===
