#!/usr/bin/env python
#
# Compute syst covariance matrices
#
# Jun 2017: Original script written by D.Scolnic for Pantheon
#
# Oct 2020: Re-written by S.Hinton for Pippin compatibility,
#            and unbinned option
#
# Dec 2020: R.Kessler, S.Hinton
#  Define JLA and BBC method. JLA is the legacy method with no changes.
#  BBC method updates:
#    + include optional COV_stat matrix
#    + set cov_syst=0 by ignoring key, or setting to null
#        (i.e., no longer need all-zero matrix)
#    + enable rebinning in c,x1,logmass space
#    + more clear table headers
#
# Feb 15 2021 RK
#    + new command line args to replace contents of input file.
#      The tricky one is -m --muopt to select only one MUOPT.
#
# May 14 2021 RK -
#   + fix bug in get_cov_from_covopt to allow pad spaces between [] terms.
#
# Aug 02 2021 Dillon
#    + update to handle option to subtract MUERR_VPEC
#
#
# Sep 22 2021 R.Kessler
#    + major overhaul to create a "standard" output under OUTDIR that is
#      intended for multiple cosmology fitting programs. 
#      If COSMOMC_TEMPALTES_PATH is defined, the legacy /cosmomc 
#      subDir is created as well.
#
# Oct 01 2021 R,Kessler
#    + create seprate write_HD_binned and write_HD_unbinned functions
#    + output HD is fitres format (same format as other SNANA codes)
#    + unbinned HD includes CID IDSURVEY MUERR_VPEC
#    + comment out df = df.set_index(["IDSURVEY", "CID"]) and pray
#
# Oct 6 2021: INFO.YML now includes ISDATA_REAL so that cosmology-fitting
#             programs can auto-detect real data and apply default blinding.
#             Only works on SALT2mu/BBC jobs run after Oct 6 2021.
#
# Nov 23 2021 
#   + for unbinned replace MUERR -> MUERR_RENORM
#   + rebin working, but <z> still not quite right.
#
# Jan 20 2022: fix to work if there is no SYS_SCALE_FILE input
# Apr 22 2022: write standard cov matrix gzipped.
#              Cov for cosmomc still unzipped until somebody shows
#              that cosmomc can read gzip file.
#
# Apr 27 2022: reduce cov labels to new row and diag only -> 
#               reduce output file size by ~50%
#
# Apr 30 2022: write MUERR_SYS = sqrt(COVSYS_DIAG) to hubble diagram
#              (diagnostic only; not used in cosmology fit)
#
# May 2 2022 A.Mitra: check for unitary matrix and pos-definite.
# Jul 19 2022 RK 
#   + define INFO_YML_FILENAME as global
#   + write ISDATA_REAL flag to HD header
#
# Aug 2022 P. Armstrong: create hubble_diagram for systematics
#
# Sep 28 2022 RK - 
#   + new optional inputs --sys_scale_file 
#   + read optional 3rd arg to scale errors in COVOPTS
# Sep 30 2022 RK
#   + optional --yaml_file arg to communicate with sumbit_batch_jobs
#
# Oct 10 2022 A.Mitra and R.K.
#   + New input arguments --label_cov_rows with default = False.
#
# Oct 13 2022 RK 
#   + for INFO.YAML, add HD file name and name of each covsys file.
#   + change cov format from 12.8f to 12.6e to get 7 digits of 
#     precision instead of only 3 or 4 digits using 12.8f
#
# Nov 9 2022 RK
#   + replace zCMB with zHD in hubble_diagram.txt
#   + for unbinned, write correct zHEL (no vpec correction)
#   + write zHD and zHEL comments at top of hubble_diagram.txt
#
# Dec 7 2022 RK 
#   + fix bug in prep_config(); affects pippin integration
#   + write SNANA_VERSION to INFO.YML file
#   + new --nosys arg
#
# Mar 08 2023 RK 
#     - read VERSION_PHOTOMETRY from hubble diagram table header
#     - read cospar_biascor from sim-README
#     - write above info to INFO.YML (to pass on to cosmology fitting codes)
#
# ===============================================

import os, argparse, logging, shutil, time, subprocess
import re, yaml, sys, gzip, math
import numpy  as np
import pandas as pd
from   pathlib import Path
from   functools import reduce
from   sklearn.linear_model import LinearRegression
import seaborn as sb
import matplotlib.pyplot as plt

import astropy.units as u
from   astropy.cosmology import Planck13, z_at_value
from   astropy.cosmology import FlatLambdaCDM

SUFFIX_M0DIF  = "M0DIF"
SUFFIX_FITRES = "FITRES"

PREFIX_COVSYS     = "covsys"
HD_FILENAME       = "hubble_diagram.txt"
INFO_YML_FILENAME = "INFO.YML"

VARNAME_CID          = "CID"  # for unbinned fitres files
VARNAME_ROW          = "ROW"  # for binned M0DIF files
VARNAME_IDSURVEY     = "IDSURVEY"
VARNAME_MU           = "MU"
VARNAME_M0DIF        = "M0DIF"
VARNAME_MUDIF        = "MUDIF"
VARNAME_MUDIFERR     = "MUDIFERR"
VARNAME_MUREF        = "MUREF"
VARNAME_MURES        = "MURES"
VARNAME_MUERR        = "MUERR"
VARNAME_MUERR_VPEC   = "MUERR_VPEC"
VARNAME_MUERR_RENORM = "MUERR_RENORM"
VARNAME_MUERR_SYS    = "MUERR_SYS"
VARNAME_zHD          = "zHD"
VARNAME_zHEL         = "zHEL"

VARNAME_iz     = "IZBIN"
# xxx mark delete VARNAME_z      = "z"  # note that zHD is internally renamed z
VARNAME_x1     = "x1"
VARNAME_c      = "c"

VARNAME_NEVT_BIN = 'NEVT'

SUBDIR_COSMOMC = "cosmomc"

# keys in header of FITRES and M0DIF tables output by SALT2mu
KEYNAME_ISDATA             = "ISDATA_REAL" 
KEYNAME_VERSION_PHOTOMETRY = "VERSION_PHOTOMETRY"

KEYNAME_SYS_SCALE_FILE = "SYS_SCALE_FILE"

m_REF = 0  # MUOPT reference number for cov
f_REF = 0  # FITOPT reference number for cov


# define list of sim-input keys for cosmology params; needed to recover biasCor cospar
KEY_DOCANA    = "DOCUMENTATION"
KEY_SIM_INPUT = "INPUT_KEYS_SNIaMODEL0"
KEYLIST_COSPAR_SIM = [ 'OMEGA_MATTER', 'OMEGA_LAMBDA', 'w0_LAMBDA', 'wa_LAMBDA',
                       'MUSHIFT' ]

# ============================
def setup_logging():

    #logging.basicConfig(level=logging.DEBUG, 
    logging.basicConfig(level=logging.INFO, 
        format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)


def read_yaml(path):
    path_expand = os.path.expandvars(path)
    logging.debug(f"Reading YAML from {path_expand}")
    with open(path_expand) as f:
        return yaml.safe_load(f.read())


def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "HELP menu for config options"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of the yml config file to run"
    parser.add_argument("input_file", help=msg, 
                        nargs="?", default=None)

    msg = "override INPUT_DIR (=BBC output)" 
    parser.add_argument("-i", "--input_dir", help=msg, 
                        nargs='?', type=str, default=None)
    
    msg = "override OUTDIR" 
    parser.add_argument("-o", "--outdir", help=msg, 
                        nargs='?', type=str, default=None)
    
    msg = "override data VERSION"
    parser.add_argument("-v", "--version", help=msg, 
                        nargs='?', type=str, default=None)
    
    msg = "override METHOD (e.g., JLA, BBC)" 
    parser.add_argument("--method", help=msg, 
                        nargs='?', type=str, default=None)

    msg = "use only this muopt number"
    parser.add_argument("-m", "--muopt", help=msg, 
                        nargs='?', type=int, default=-1 ) 

    msg = f"override {KEYNAME_SYS_SCALE_FILE} in the input file"
    parser.add_argument("--sys_scale_file", help=msg, 
                        nargs='?', type=str, default=None ) 

    msg = f"no systematics (stat only)"
    parser.add_argument("--nosys", help=msg, action="store_true")

    # xxx not yet ...
    #msg = "scale all systematics by this factor"
    #parser.add_argument("--sys_scale_global", help=msg, 
    #                    nargs='?', type=float, default=1.0 ) 
    
    msg = "Add labels in covariance matrix output (visual debug)"                                                                                      
    parser.add_argument("--label_cov_rows", help=msg,                                                                                                      
                        nargs='?', type=bool, default=False )

    msg = "Use each SN instead of BBC binning"
    parser.add_argument("-u", "--unbinned", help=msg, action="store_true")

    msg = "number of x1 bins (default=0)"
    parser.add_argument("--nbin_x1", help=msg, 
                        nargs='?', type=int, default=0 )

    msg = "number of c bins (default=0)"
    parser.add_argument("--nbin_c", help=msg, 
                        nargs='?', type=int, default=0 )

    #msg = "rebin args; 0 (unbinned) or e.g., c:5,x1:2 (5 c bins,2 x1 bins)"
    #parser.add_argument("--rebin", help=msg, nargs='+', type=str )

    msg = "Subtract MUERR(VPEC) from MUERR. Forces unbinned."
    parser.add_argument("-s", "--subtract_vpec", help=msg, 
                        action="store_true")

    msg = "Produce a hubble diagram for every systematic"
    parser.add_argument("--systematic_HD", help=msg, action="store_true")

    msg = "output yaml file (for submit_batch_jobs)"
    parser.add_argument("--yaml_file", help=msg, 
                        nargs='?', type=str, default=None )

    # parse it
    args = parser.parse_args()
    if args.subtract_vpec:
        args.unbinned = True

    if (args.nbin_x1>0 or args.nbin_c>0): 
        args.unbinned = False  # recommended by D.Brout, Aug 8 2022

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    if args.HELP : 
        print_help_menu()

    return args
    # end get_args

def print_help_menu():
    menu = f"""

   HELP MENU for create_covariance input file
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INPUT_DIR:  <BBC output dir from submit_batch_jobs or from Pippin>
VERSION:    <subDir under INPUT_DIR>  

{KEYNAME_SYS_SCALE_FILE}: <optional yaml file with systematic scales>
  [allows legacy key SYSFILE: ]

# Define extra cov matrices using covsys options:
#  [COVOPT-label]  [FITOPT-label,MUOPT-label,errSystScale]
#   and no pad spaces allowed inside []
COVOPTS:
- '[NOSYS]  [=DEFAULT,=DEFAULT]' # no syst 
- '[PS1CAL] [+PS1CAL,=DEFAULT]'  # PS1CAL syst only MUOPT=0
- '[ALLCAL] [+CAL,=DEFAULT]'     # All syst with 'CAL' in name, MUOPT=0
- '[SALT2]  [+SALT2,=DEFAULT]'   # SALT2 syst only, MUOPT=0
- '[NOCAL]  [-CAL,=DEFAULT]'     # all syst except for calib, MUOPT=0
- '[SCAT]   [=DEFAULT,+SCAT]     # FITOPT=0, MUOPTs with SCAT in label

# optional 3rd arg scales syst error
- '[SCAT]   [=DEFAULT,+SCAT,1.2] # errSys *= 2, and covSys*= 1.2^2
- '[ALL]    [,,1.4]              # scale ALL errSys by 1.4, covSys*= 1.4^2

MUOPT_SCALES:  # replace scales in {KEYNAME_SYS_SCALE_FILE}
  CLAS_SNIRF:  1.0
  SCATTER_C11: 1.0

OUTDIR:  <name of output directory>

# Some of the above keys have command line overrides (-h for options):
  create_covariance.py <configFile> \\
     -i INPUT_DIR \\
     -v VERSION   \\
     --sys_scale_file SYS_SCALE_FILE \\

# specialized keys for CosmoMC
COSMOMC_TEMPLATES_PATH: <out path of ini files for Cosmomc>
  [allows legacy key COSMOMC_TEMPLATES:]

COSMOMC_DATASET_FILE: <out file with cosmomc instructions>
  [allows legacy key DATASET_FILE: ]

    """
    print(f"{menu}")

    sys.exit()

    # end print_help_menu

def get_snana_version():
    # fetch snana version that includes tag + commit;
    # e.g., v11_05-4-gd033611.
    # Use same git command as in Makefile for C code
    SNANA_DIR        = os.environ['SNANA_DIR']
    cmd = f"cd {SNANA_DIR};  git describe --always --tags"
    ret = subprocess.run( [ cmd ], cwd=os.getcwd(),
                      shell=True, capture_output=True, text=True )
    snana_version = ret.stdout.replace('\n','')
    return snana_version

def read_header_info(hd_file):

    # Created oct 6 2021 by R.Kessler
    # read first few lines of hd_file and look for ISDATA_REAL  key
    # that is in a comment field ... hence read as yaml.
    # SALT2mu begain wriring ISDATA_REAL key on Oct 6 2021;
    # Function returns ISDATA_REAL = 0 or 1 if key is found;
    # returns -1 (unknown) if key not found.
    # Mar 2023: rename isdata_real to read_header_info

    with gzip.open(hd_file, 'r') as f:
        line_list = f.readlines()

    isdata_real = -1  # init to unknonw
    maxline_read = 10 # bail after this many lines
    nline_read   = 0
    isdata_real  = -1
    header_info = {}  # define output dictionary

    found_isdata   = False
    found_ver_phot = False

    for line in line_list:
        line  = line.rstrip()  # remove trailing space and linefeed  
        line  = line.decode('utf-8')
        wd_list = line.split()

        # xxx mark delete: if key_isdata in wd_list:
        if any(KEYNAME_ISDATA in s for s in wd_list):  
            key_isdata     = f"{KEYNAME_ISDATA}:"
            j              = wd_list.index(key_isdata)
            if j > 0 :
                isdata_real    = int(wd_list[j+1])            
                found_isdata   = True

        if any(KEYNAME_VERSION_PHOTOMETRY in s for s in wd_list):  
            key = wd_list[1].replace(':','')  # item 0 is hash, item 1 is key
            if len(wd_list) > 2 :
                arg = wd_list[2].split(',')
            else:
                arg = [ "UNKNOWN" ]
            header_info[key] = arg
            found_ver_phot   = True

        nline_read += 1
        if nline_read == maxline_read or isdata_real>=0 : break

    header_info[KEYNAME_ISDATA] = isdata_real
    logging.info(f"ISDATA_REAL = {isdata_real}")

    # - - - - - give warnings for key(s) not found - - - - - -
    found_list = [ found_isdata, found_ver_phot ]
    key_list   = [ KEYNAME_ISDATA, KEYNAME_VERSION_PHOTOMETRY ]

    for found,key in zip(found_list,key_list):
        if not found:
            logging.warning(f"did not find {key} in table header")

    return header_info

    # end read_header_info

def load_hubble_diagram(hd_file, args, config):

    # read single M0DIF or FITRES file from BBC output,
    # and return contents.

    if not os.path.exists(hd_file):
        raise ValueError(f"Cannot load Hubble diagram data from {hd_file}" \
                         f" - it doesnt exist")

    df = pd.read_csv(hd_file, delim_whitespace=True, comment="#")
    logging.debug(f"\tLoaded data with Nrow x Ncol {df.shape} from {hd_file}")

    #sys.exit("\n xxx DEBUG STOP xxx\n")

    # Do a bit of data cleaning: replace 999 with nan 
    # (beware that 999 in BBC output means no info, does not mean nan)
    df = df.replace(999.0, np.nan)

    # M0DIF doesnt have MU column, so add it back in
    # For FITRES file with all events, do nothing sinve MU exists
    if "MU" not in df.columns:
        df[VARNAME_MU] = df[VARNAME_MUREF] + df[VARNAME_MUDIF]
    if VARNAME_MUERR not in df.columns:
        df[VARNAME_MUERR] = df[VARNAME_MUDIFERR]

    # Sort by CID for unbinned; sort by z for M0DIF
    # --> ensure direct subtraction comparison
    if "CID" in df.columns:
        df["CID"] = df["CID"].astype(str)
        df = df.sort_values([ "CID"])
        df = df.rename(columns={"MUMODEL": VARNAME_MUREF})
        # xxx mark df = df.rename(columns={"zHD": VARNAME_z, "MUMODEL": VARNAME_MUREF})
        if args.subtract_vpec:
            msgerr = f"Cannot subtract VPEC because MUERR_VPEC " \
                     f"doesn't exist in {hd_file}"
            assert "MUERR_VPEC" in df.columns, msgerr

            df[VARNAME_MUERR] = np.sqrt(df[VARNAME_MUERR]**2 - \
                                        df[VARNAME_MUERR_VPEC]**2)  
            logging.debug(f"Subtracted {VARNAME_MUERR_VPEC} " \
                          f"from {VARNAME_MUERR}")

        df['CIDstr'] = df['CID'].astype(str) + "_" + \
                       df['IDSURVEY'].astype(str)

        # we need to make a new column to index the dataframe on 
        # unique rows for merging two different systematic dfs - Dillon
        df['CIDindex'] = df['CID'].astype(str) + "_" + \
                         df['IDSURVEY'].astype(str)
        df = df.set_index("CIDindex")

    elif 'z' in df.columns:  # for M0DIF file that has z instead of zHD
        # should not need z-sorting here for M0DIF, but what the heck
        df = df.rename(columns={"z": VARNAME_zHD})  # Nov 9 2022
        df = df.sort_values(VARNAME_zHD)
    elif VARNAME_zHD in df.columns:
        df = df.sort_values(VARNAME_zHD) # sort, but likely not needed
        pass

    # - - - -  -
    
    return df
    # end load_hubble_diagram


def get_hubble_diagrams(folder, args, config):

    # return table for each Hubble diagram 

    is_rebin    = config['nbin_x1'] > 0 and config['nbin_c'] > 0
    is_unbinned = args.unbinned 
    is_binned   = not (is_unbinned or is_rebin)

    folder_expand = Path(os.path.expandvars(folder))
    logging.debug(f"Loading all data files in {folder_expand}")
    HD_list     = {}
    infile_list = []
    label_list  = []
    first_load  = True

    for infile in sorted(os.listdir(folder_expand)):

        is_M0DIF    = f".{SUFFIX_M0DIF}"  in infile
        is_FITRES   = f".{SUFFIX_FITRES}" in infile and "MUOPT" in infile

        # Feb 15 2021 RK - check option to select only one of the MUOPT
        if args.muopt >= 0 :
            muopt_num = f"MUOPT{args.muopt:03d}"
            is_M0DIF  = is_M0DIF and muopt_num in infile
            
        do_load = False
        if is_binned   and is_M0DIF  : do_load = True 
        if is_unbinned and is_FITRES : do_load = True
        if is_rebin    and is_FITRES : do_load = True

        if do_load :
            replace_list = [".gz", f".{SUFFIX_M0DIF}", f".{SUFFIX_FITRES}" ]
            label        = infile

            for rep in replace_list:
                label = label.replace(rep,"") # remove suffix

            infile_list.append(infile)
            label_list.append(label)
            hd_file = folder_expand/infile
            # grab contents of every M0DIF(binned) or FITRES(unbinned) file 
            HD_list[label] = load_hubble_diagram(hd_file, args, config)
            if first_load:  
                header_info = read_header_info(hd_file)
            first_load = False

    #sys.exit(f"\n xxx\n result = {result} \n")

    config[KEYNAME_ISDATA] = header_info[KEYNAME_ISDATA]
    config['header_info']  = header_info

    # - - - - -  -
    # for unbinned or rebin, select SNe that are common to all 
    # FITOPTs and MUOPTs
    if is_unbinned or is_rebin:
        HD_list = get_common_set_of_sne(HD_list)
        
    # - - - - -    
    # df['e'] = e.values or  df1 = df1.assign(e=e.values)

    if is_rebin : 
        label0 = label_list[0]
        get_rebin_info(config,HD_list[label0])
        for label in label_list:
            print(f"   Rebin {label}")
            HD_list[label]['iHD'] = config['col_iHD']
            HD_list[label] = rebin_hubble_diagram(config,HD_list[label])

    if is_unbinned :
        HD_list = update_MUERR(HD_list)

    return HD_list
    # end get_hubble_diagrams

def get_common_set_of_sne(datadict):

    combined = None
    n_file = 0
    for label, df in datadict.items():
        if combined is None:
            combined = df.index
            missing = None
        else:
            missing = combined.difference(df.index)
            combined = combined.intersection(df.index)
        
        n_sn = combined.shape[0]
        n_file += 1
        logging.info(f"Common set from {label} has {n_sn} elements")
        logging.info(f"\n{missing}")
        assert combined.shape[0], "\t No common SNe ?!?!?"

    n_sn = combined.shape[0]
    logging.info(f"Common SN set from {n_file} files has {n_sn} events")
    for label, df in datadict.items():
        datadict[label] = df.loc[combined, :]

    return datadict

def update_MUERR(HDs):
    for label,df in HDs.items():
        if VARNAME_MUERR_RENORM in df.columns:
            HDs[label][VARNAME_MUERR] =  df[VARNAME_MUERR_RENORM]
    return HDs

def get_rebin_info(config,HD):

    # Dec 29 2020
    # Figure out x1 and c rebinning for input Hubble diagram table 'HD'.
    # Load bin info into config.

    nbin_x1 = config['nbin_x1']  # user input
    nbin_c  = config['nbin_c']   # user input
    nbin_z  = len(HD[VARNAME_iz].unique())  # from FITRES file

    epsilon = 1.0E-6

    zmin  = HD[VARNAME_zHD].min() - epsilon 
    zmax  = HD[VARNAME_zHD].max() + epsilon 
    x1min = HD[VARNAME_x1].min()  - epsilon 
    x1max = HD[VARNAME_x1].max()  + epsilon
    cmin  = HD[VARNAME_c].min()   - epsilon
    cmax  = HD[VARNAME_c].max()   + epsilon

    logging.info(f"\t z (min,max) = {zmin:.3f} {zmax:.3f}   (nbin={nbin_z}) ")
    logging.info(f"\t x1(min,max) = {x1min:.3f} {x1max:.3f}  (nbin={nbin_x1}) ")
    logging.info(f"\t c (min,max) = {cmin:.3f} {cmax:.3f}  (nbin={nbin_c}) ") 
    
    bins_x1 = np.linspace(x1min, x1max, num=nbin_x1+1, endpoint=True)
    bins_c  = np.linspace(cmin,   cmax, num=nbin_c+1,  endpoint=True)

    # get ix1 and ic index for each event.
    # Beware that digitize return fortran-like index starting at 1
    col_ix1 = np.digitize( HD[VARNAME_x1], bins_x1)
    col_ic  = np.digitize( HD[VARNAME_c],  bins_c)
    col_iz  = HD[VARNAME_iz]

    col_ix1 -= 1  # convert to C like index starting at 0
    col_ic  -= 1  # idem

    # convert three 1D indices into a single 1D index
    # i = x + (y * max_x) + (z * max_x * max_y)
    # x,y,z -> ic,ix1,iz
    col_iHD  = col_ic + (col_ix1*nbin_c)  + (col_iz*nbin_x1*nbin_c)

    logging.info(f"\t iHD(min,max) = {col_iHD.min()}  {col_iHD.max()}  ")
    dump_flag = False
    if dump_flag :
        print(f" xxx bins_x1 = {bins_x1} ")
        print(f" xxx bins_c  = {bins_c} ")
        #print(f"\n xxx col(ix) = \n{col_ix1}")
        #print(f"\n xxx col(ic) = \n{col_ic}")
        print(f"\n xxx col(iHD) = \n{col_iHD}")

    # - - - - -
    # load bins into config
    config['nbin_z']         = nbin_z
    config['bins_x1']        = bins_x1
    config['bins_c']         = bins_c
    config['nbin_HD']        = nbin_x1 * nbin_c * nbin_z
    config['col_iHD']        = col_iHD

    config['col_iz']         = col_iz
    config['col_ix1']        = col_ix1
    config['col_ic']         = col_ic

    #sys.exit(f"\n xxx DEBUG DIE from get_rebins xxx ")
    return

    # end get_rebin_info

def rebin_hubble_diagram(config, HD_unbinned):

    # Dec 29 2020
    # rebin unbinned results based on nbin_xxx user inputs
    #
    # 

    nbin_x1      = config['nbin_x1']
    nbin_c       = config['nbin_c']
    nbin_HD      = config['nbin_HD']
    col_iz       = config['col_iz']
    col_ix1      = config['col_ix1']
    col_ic       = config['col_ic']

    col_iHD      = HD_unbinned['iHD']
    col_c        = HD_unbinned[VARNAME_c]
    col_z        = HD_unbinned[VARNAME_zHD]
    col_mures    = HD_unbinned[VARNAME_MURES]
    col_mudif    = HD_unbinned[VARNAME_M0DIF]
    col_muref    = HD_unbinned[VARNAME_MUREF]
    col_mu       = HD_unbinned[VARNAME_MU]
    col_muerr    = HD_unbinned[VARNAME_MUERR]
    col_muerr_renorm = HD_unbinned[VARNAME_MUERR_RENORM]
    HD_rebin_dict = { VARNAME_ROW: [], 
                      VARNAME_zHD: [], VARNAME_MU: [] , VARNAME_MUERR: [],
                      VARNAME_MUREF: [], VARNAME_NEVT_BIN: [] }

    wgt0      = 1.0/(col_muerr*col_muerr)
    wgt1      = 1.0/(col_muerr_renorm*col_muerr_renorm)
    wgtmuref  = wgt0 * col_muref
    wgtmudif  = wgt1 * col_mudif
    wgtmu     = wgt1 * col_mu

    # - - - - -
    OM_ref = 0.30  # .xyz should read OM_ref from BBC output
    zcalc_grid, mucalc_grid = get_HDcalc(OM_ref)

    for i in range(0,nbin_HD):
        binmask       = (col_iHD == i)
        wgt0sum       = np.sum(wgt0[binmask])
        wgt1sum       = np.sum(wgt1[binmask])
        nevt          = np.sum(binmask)
        if nevt == 0 :   continue

        iz  = col_iz[binmask][0]
        ix1 = col_ix1[binmask][0]
        ic  = col_ic[binmask][0]

        row_name     = f"BIN{i:04d}_z{iz:02d}-x{ix1}-c{ic}"
        muerr_wgtavg = math.sqrt(1.0/wgt1sum)
        muref_wgtavg = np.sum(wgtmuref[binmask]) / wgt0sum
        #mudif_wgtavg = np.sum(wgtmudif[binmask]) / wgt1sum
        #mu_wgtavg    = muref_wgtavg + mudif_wgtavg
        mu_wgtavg    = np.sum(wgtmu[binmask]) / wgt1sum
        z_invert     = np.interp(muref_wgtavg, mucalc_grid.value, zcalc_grid)

        #print(f"  xxx rebinned mu[{i:2d}] = " \
        #      f"{mu_wgtavg:7.4f} +_ {muerr_wgtavg:7.4f}")

        HD_rebin_dict[VARNAME_ROW].append(row_name)
        HD_rebin_dict[VARNAME_zHD].append(z_invert)
        HD_rebin_dict[VARNAME_MU].append(mu_wgtavg)
        HD_rebin_dict[VARNAME_MUERR].append(muerr_wgtavg)
        HD_rebin_dict[VARNAME_MUREF].append(muref_wgtavg)
        HD_rebin_dict[VARNAME_NEVT_BIN].append(nevt)

    # - - - - -
    HD_rebin = pd.DataFrame(HD_rebin_dict)
    return HD_rebin

    # end rebin_hubble_diagram

def get_HDcalc(OM):

    # return calculated HD for input OM and flat LCDM.
    zmin    = 1.0E-04
    zmax    = 4.0
    z_grid  = np.geomspace(zmin, zmax, 500)
    cosmo   = FlatLambdaCDM(H0=70, Om0=OM)
    mu_grid = cosmo.distmod(z_grid)
    return z_grid, mu_grid

def get_fitopt_muopt_from_name(name):
    f = int(name.split("FITOPT")[1][:3])
    m = int(name.split("MUOPT")[1][:3])
    return f, m


def get_name_from_fitopt_muopt(f, m):
    return f"FITOPT{f:03d}_MUOPT{m:03d}"


def get_fitopt_scales(lcfit_info, sys_scales):
    """ Returns a dict mapping FITOPT numbers to (label, scale) """
    fitopt_list = lcfit_info["FITOPT_OUT_LIST"]

    if fitopt_list is None: return None

    result = {}
    for number, _, label, _ in fitopt_list:
        if label != "DEFAULT" and label is not None:
            if label not in sys_scales:
                pass
                #logging.warning(f"No FITOPT scale found for {label}")
        scale = sys_scales.get(label, 1.0)
        d = int(number.replace("FITOPT", ""))
        result[d] = (label, scale)
    return result

def get_cov_from_diff(df1, df2, scale):
    """ Returns both the covariance contribution and summary stats 
    (slope and mean abs diff) """

    len1 = df1[VARNAME_MU].size
    len2 = df2[VARNAME_MU].size
    errmsg =   f"Oh no, len1={len1} and len2={len2}; " \
               f"looks like you have a different number of bins/supernova " \
               f"for your systematic and this is not good."
    assert len1 == len2, errmsg


    diff = scale * ((df1[VARNAME_MU] - df1[VARNAME_MUREF]) - \
                    (df2[VARNAME_MU] - df2[VARNAME_MUREF])).to_numpy()
    diff[~np.isfinite(diff)] = 0
    cov = diff[:, None] @ diff[None, :]

    # Determine the gradient using simple linear regression
    reg = LinearRegression()
    weights = 1 / np.sqrt(0.003**2 + df1[VARNAME_MUERR].to_numpy()**2 \
                          + df2[VARNAME_MUERR].to_numpy()**2)
    mask = np.isfinite(weights)

    reg.fit(df1.loc[mask, [VARNAME_zHD]], diff[mask], 
            sample_weight=weights[mask])
    coef = reg.coef_[0]

    mean_abs_deviation = np.average(np.abs(diff), weights=weights)
    max_abs_deviation = np.max(np.abs(diff))
    return cov, (coef, mean_abs_deviation, max_abs_deviation)
    # end get_cov_from_diff

def get_cov_from_covfile(data, covfile, scale):
    covindf = pd.read_csv(covfile,float_precision='high',low_memory=False)
    covindf['CID1'] = covindf['CID1'].astype(str)+"_"+covindf['IDSURVEY1'].astype(str)
    covindf['CID2'] = covindf['CID2'].astype(str)+"_"+covindf['IDSURVEY2'].astype(str)
    covout = np.zeros((len(data),len(data)))
    for i,row in covindf.iterrows():
        if len(np.argwhere(data['CIDstr'].array == row['CID1'])) == 0:
            print(row['CID1'],'1 missing from output/cosmomc/data_wCID.txt')
            continue
        if len(np.argwhere(data['CIDstr'].array == row['CID2'])) == 0:
            print(row['CID2'],'2 missing from output/cosmomc/data_wCID.txt')
            continue
        ww1 = np.argwhere(data['CIDstr'].array == row['CID1'])[0][0]
        ww2 = np.argwhere(data['CIDstr'].array == row['CID2'])[0][0]
        #if ww1 == ww2:
        #    print('skipping, not doing same SN diagonals')
        #    continue
        covout[ww1,ww2] = row['MU_COV']

    return covout, (0, 0, 0)

def get_contributions(m0difs, fitopt_scales, muopt_labels, 
                      muopt_scales, extracovdict):
    """ Gets a dict mapping 'FITOPT_LABEL|MUOPT_LABEL' to covariance)"""
    result, slopes = {}, []

    for name, df in m0difs.items():
        f, m = get_fitopt_muopt_from_name(name)
        logging.debug(f"Determining contribution for FITOPT {f}, MUOPT {m}")

        # Get label and scale for FITOPTS and MUOPTS. Note 0 index is DEFAULT
        if f == 0:
            fitopt_label = "DEFAULT"
            fitopt_scale = 1.0
        else:
            fitopt_label, fitopt_scale = fitopt_scales[f]

        muopt_label = muopt_labels[m] if m else "DEFAULT"
        muopt_scale = muopt_scales.get(muopt_label, 1.0)

        logging.debug(f"FITOPT {f} has scale {fitopt_scale}, " \
                      f"MUOPT {m} has scale {muopt_scale}")

        # Depending on f and m, compute the contribution to the covariance matrix
        if f == f_REF and m == m_REF :
            # This is the base file, so don't return anything.
            # CosmoMC will add the diag terms itself.
            cov = np.zeros((df[VARNAME_MU].size, df[VARNAME_MU].size))
            summary = 0, 0, 0
        elif m != m_REF :
            # This is a muopt, to compare it against the MUOPT000 for the same FITOPT
            df_compare = m0difs[get_name_from_fitopt_muopt(f, 0)]
            cov,summary = get_cov_from_diff(df, df_compare, fitopt_scale*muopt_scale)
        else:
            # This is a fitopt with MUOPT000, compare to base file
            df_compare = m0difs[get_name_from_fitopt_muopt(0, 0)]
            cov, summary = get_cov_from_diff(df, df_compare, fitopt_scale)

        result[f"{fitopt_label}|{muopt_label}"] = cov
        slopes.append([name, fitopt_label, muopt_label, *summary])

    #now loop over extracov_labels
    for key,value in extracovdict.items():
        cov, summary = get_cov_from_covfile(df, os.path.expandvars(value), muopt_scales[key])
        fitopt_label = 'DEFAULT'
        muopt_label = key
        result[f"{fitopt_label}|{muopt_label}"] = cov
        slopes.append([name, fitopt_label, muopt_label, *summary])

    summary_df = pd.DataFrame(slopes, columns=["name", "fitopt_label", "muopt_label", "slope", "mean_abs_deviation", "max_abs_deviation"])
    summary_df = summary_df.sort_values(["slope", "mean_abs_deviation", "max_abs_deviation"], ascending=False)
    return result, summary_df


def apply_filter(string, pattern):
    """ Used for matching COVOPTs to FITOPTs and MUOPTs"""
    string, pattern = string.upper(), pattern.upper()
    if pattern.startswith("="):
        return string == pattern[1:]
    elif pattern.startswith("+"):
        return pattern[1:] in string
    elif pattern.startswith("-"):
        return pattern[1:] not in string
    elif pattern == "":
        return True
    else:
        raise ValueError(f"Unable to parse COVOPT matching pattern {pattern}")


def get_cov_from_covopt(covopt, contributions, base, calibrators):

    # Parse covopts (from input config file) that look like 
    #        "[cal] [+cal,=DEFAULT]"
    #
    # Split covopt into two terms so that extra pad spaces don't
    # break findall command (RK May 14 2021)
    #
    # 9.29.2022 RK - optional 3rd arg with sys scale ?
    #         "[cal] [+cal,=DEFAULT, SCALE=1.3]" 

    
    covopt_list = covopt.split() # break into two terms
    tmp0 = covopt_list[0]
    tmp1 = covopt_list[1]

    bracket_content0 = re.findall(r"\[(.*)\]",  tmp0)[0]
    bracket_content1 = re.findall(r"\[(.*)\]",  tmp1)[0]
    bracket_content1_list = bracket_content1.split(',')

    #  Sep 29 2022 RK - refactor parsing to allow 2 or 3 comma-sep
    #  items in 2nd bracket.
    label         = bracket_content0
    fitopt_filter = bracket_content1_list[0]
    muopt_filter  = bracket_content1_list[1]
    covopt_scale     = 1.0
    if len(bracket_content1_list) > 2 :  # err_scale (3rd item) is optional
        sig_scale    = float(bracket_content1_list[2])
        covopt_scale = sig_scale * sig_scale
        
    # generic message-content for debug or error
    msg_content1 =  \
        f"COV({label}): FITOPT/MUOPT filters = " \
        f"'{fitopt_filter}' / {muopt_filter} | " \
        f" covopt_scale={covopt_scale}"

    # xxx mark delete  print(f" xxx refac {msg_content1}")

    fitopt_filter = fitopt_filter.strip()
    muopt_filter  = muopt_filter.strip()
    logging.debug(f"Compute {msg_content1}")

    final_cov = None

    if calibrators: # Cepheid calibrators don't have z-syst
        mask_calib = base.reset_index()["CID"].isin(calibrators)

    for key, cov in contributions.items():

        fitopt_label, muopt_label = key.split("|")

        apply_fitopt = apply_filter(fitopt_label, fitopt_filter)
        apply_muopt  = apply_filter(muopt_label,  muopt_filter)
        apply_vpec   = apply_filter(fitopt_label, "+VPEC") or \
                       apply_filter(muopt_label,  "+VPEC")
        apply_zshift = apply_filter(fitopt_label, "+ZSHIFT") or \
                       apply_filter(muopt_label,  "+ZSHIFT")

        if apply_fitopt and apply_muopt :
            if final_cov is None:
                final_cov = cov.copy()*0
            if True:
                # If calibrators and  VPEC term, filter out calib
                if calibrators and (apply_vpec or apply_zshift):
                    print(f"FITOPT {fitopt_label} MUOPT {muopt_label} " \
                          f"ignored for calibrators...")
                    cov2 = cov.copy()
                    cov2[mask_calib, :] = 0
                    cov2[:, mask_calib] = 0
                    # xxx mark delete final_cov += cov2
                    final_cov += cov2 * covopt_scale
                else:
                    # xxx mark delete final_cov += cov
                    final_cov += cov * covopt_scale

    assert final_cov is not None,  f"No syst matches {msg_content1} " 

    # Validate that the final_cov is invertible
    try:
        # CosmoMC will add the diag terms, 
        # so lets do it here and make sure its all good
        effective_cov = final_cov + np.diag(base[VARNAME_MUERR] ** 2)

        # First just try and invert it to catch singular matrix errors
        precision = np.linalg.inv(effective_cov)

        # A.Mitra, May 2022
        # Check if matrix is unitary and pos-definite.
        pr = np.dot(effective_cov,precision)
        pr = np.round(pr,decimals=3)
        flag = is_unitary(np.round(pr,decimals=2))
        if flag :
            logging.info(f"{label} Matrix is UNITARY")
        else :
            logging.info(f"WARNING: {label} Matrix is not UNITARY")

        flag2 = is_pos_def(np.round(effective_cov,decimals=3))
        if flag2:
            logging.info(f"{label} Matrix is Positive-Definite")
        else :
            logging.warning(f"{label} Matrix is not Positive-Definite")
        
        # check that COV is well conditioned to deal with float precision
        epsilon = sys.float_info.epsilon
        cond = np.linalg.cond(effective_cov)
        assert cond < 1 / epsilon, "Cov matrix is ill-conditioned and cannot be inverted"
        logging.info(f"Covar condition for COVOPT {label} is {cond:.3f}")

    except np.linalg.LinAlgError as ex:
        logging.exception(f"Unable to invert covariance matrix for COVOPT {label}")
        raise ex

    return label, final_cov

def is_unitary(matrix: np.ndarray) -> bool:
    # Created May 2022 by A.Mitra
    # Return true if input matrix is unitary
    unitary = True
    n = len(matrix)
    error = np.linalg.norm(np.eye(n) - matrix.dot( matrix.transpose().conjugate()))
    if not(error < np.finfo(matrix.dtype).eps * 10.0 *n):
        unitary = False
    return unitary
    # end is_unitary

def is_pos_def(x):
    # Created May 2022 by A.Mitra
    # return True of input matrix x is Postive Definite
    return np.all(np.linalg.eigvals(x) > 0)
    # end is_pos_def

def write_standard_output(config, args, covs, data, label_list):

    # Created 9.22.2021 by R.Kessler
    # Write standard cov matrices and HD for cosmology fitting programs;
    # e.g., wfit, CosmoSIS, firecrown ...
    # Note that CosmoMC uses a more specialized output created
    # by write_cosmomc_output().

    # P. Armstrong 05 Aug 2022
    # Add option to create hubble_diagram.txt for each systematic as well

    logging.info("")
    logging.info("   OUTPUT  ")
    unbinned       = args.unbinned
    label_cov_rows = args.label_cov_rows
    outdir         = Path(config["OUTDIR"])
    os.makedirs(outdir, exist_ok=True)

    # Apr 30 2022: get array of muerr_sys(ALL) for output
    muerr_sys_list = get_muerr_sys(covs)
            
    # - - - -
    # P. Armstrong 05 Aug 2022: Write HD for each label
    for label in label_list:
        base_name = get_name_from_fitopt_muopt(f_REF, m_REF)

        # If creating HD for nominal, write to hubble_diagram.txt
        if label == base_name:
            data_file = outdir / HD_FILENAME
        else:
            data_file = outdir / f"{label}_{HD_FILENAME}"

        if unbinned :
            write_HD_unbinned(data_file, data[label], muerr_sys_list)
        else:
            write_HD_binned(data_file, data[label], muerr_sys_list)

    
    # write covariance matrices and datasets
    opt_cov = 0  
    if label_cov_rows: opt_cov+=1
    for i, (label, cov) in enumerate(covs):

        # xxx mark delete Oct 2022
        #base_file   = f"{PREFIX_COVSYS}_{i:03d}.txt" 
        #base_file  += '.gz'  # force gzip file, Apr 22 2022

        base_file   = get_covsys_filename(i)
        covsys_file = outdir / base_file
        write_covariance(covsys_file, cov, opt_cov)
        
    return

    # end write_standard_output

def get_covsys_filename(i):
    # Created Oct 13 2022 by R.Kessler: 
    # return name of covsys file for systematic index i
    covsys_file   = f"{PREFIX_COVSYS}_{i:03d}.txt.gz"
    return covsys_file
    # end get_covsys_filename

def get_muerr_sys(covs):

    # Created April 30 2022
    # return array of muerr_sys = sqrt(COVSYS_diag)
    # to indicate size of syst. error for each HD entry.

    muerr_sys_list = None

    covsys_all = None
    for i, (label, cov) in enumerate(covs):
        if label == "ALL":
            covsys_all = cov
        
    if covsys_all is not None :
        covdiag_list = covsys_all.diagonal()
        muerr_sys_list = np.sqrt(covdiag_list)

    return muerr_sys_list
    # end get_muerr_sys

# =====================================
# cosmomc utilities

def write_cosmomc_output(config, args, covs, base):

    # Copy & modify INI files. 
    # Create covariance matrices
    # Create dataset file
 
    logging.info("")
    logging.info("# - - - - - - - - - - - - - - - - - - - - - - - -")
    logging.info(f"       OUTPUT FOR COSMOMC-JLA ")

    OUTDIR           = config["OUTDIR"]
    cosmomc_path     = config["COSMOMC_TEMPLATES_PATH"]
    dataset_file     = config["COSMOMC_DATASET_FILE"]
    unbinned         = args.unbinned

    out              = Path(OUTDIR) / SUBDIR_COSMOMC
    dataset_template = Path(cosmomc_path) / dataset_file
    dataset_files    = []

    os.makedirs(out, exist_ok=True)

    # Create lcparam file
    
    data_file      = out / f"data.txt"
    data_file_wCID = out / f"data_wCID.txt"
    prefix_covsys  = "sys"
    write_cosmomc_HD(data_file, base, unbinned)
    write_cosmomc_HD(data_file_wCID, base, unbinned, cosmomc_format=False)

    # Create covariance matrices and datasets
    opt_cov = 0     # write cov with no comments
    for i, (label, cov) in enumerate(covs):
        dataset_file = out / f"dataset_{i}.txt"
        covsyst_file = out / f"{prefix_covsys}_{i}.txt" 

        write_covariance(covsyst_file, cov, opt_cov)
        write_cosmomc_dataset(dataset_file, data_file, 
                              covsyst_file, dataset_template)
        dataset_files.append(dataset_file)

    # Copy some INI files
    ini_files = [f for f in os.listdir(cosmomc_path) if f.endswith(".ini") or f.endswith(".yml") or f.endswith(".md")]
    for ini in ini_files:
        op = Path(cosmomc_path) / ini

        if ini in ["base.ini"]:
            # If its the base.ini, just copy it
            npath = out / ini
            shutil.copy(op, npath)
        else:
            # Else we need one of each ini per covopt
            for i, _ in enumerate(covs):
                # Copy with new index
                npath = out / ini.replace(".ini", f"_{i}.ini")
                shutil.copy(op, npath)

                basename = os.path.basename(npath).replace(".ini", "")

                # Append the dataset info
                with open(npath, "a+") as f:
                    f.write(f"\nfile_root={basename}\n")
                    f.write(f"jla_dataset={dataset_files[i]}\n")


    logging.info("# - - - - - - - - - - - - - - - - - - - - - - - -")
    logging.info("")

    # end write_cosmomc_output

def write_cosmomc_dataset(path, data_file, cov_file, template_path):
    with open(template_path) as f_in:
        with open(path, "w") as f_out:
            f_out.write(f_in.read().format(data_file=data_file, 
                                           cov_file=cov_file))
    # end write_cosmomc_dataset

def write_cosmomc_HD(path, base, unbinned, cosmomc_format=True):

    if not cosmomc_format:
        if unbinned:
            varlist = [VARNAME_CID, VARNAME_IDSURVEY, VARNAME_zHD, 
                       VARNAME_MU, VARNAME_MUERR]
        else:
            varlist = [VARNAME_zHD, VARNAME_MU, VARNAME_MUERR]

        base[varlist].to_csv(path, sep=" ", index=False, float_format="%.5f")
        return

    zs   = base[VARNAME_zHD].to_numpy()
    mu   = base[VARNAME_MU].to_numpy()
    mbs  = -19.36 + mu
    mbes = base[VARNAME_MUERR].to_numpy()

    # I am so sorry about this, but CosmoMC is very particular
    logging.info(f"Write HD to {path}")
    with open(path, "w") as f:
        f.write("#name zcmb    zhel    dz mb        dmb     x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor\n")
        for i, (z, mb, mbe) in enumerate(zip(zs, mbs, mbes)):
            f.write(f"{i:5d} {z:6.5f} {z:6.5f} 0  {mb:8.5f} {mbe:8.5f} " \
                    f"0 0 0 0 0 0 0 0 0 0 0 0\n")

    # end write_cosmomc_HD

# ========= end cosmomc utils ====================


def write_HD_binned(path, base, muerr_sys_list):

    # Dec 2020
    # Write standard binned HD format for BBC method
    # Sep 30 2021: replace csv format with SNANA fitres format
    # Apr 30 3022: check muerr_sys_list

    #if "CID" in df.columns:

    unbinned = False

    logging.info(f"Write binned HD to {path}")

    wrflag_nevt   = (VARNAME_NEVT_BIN in base)
    wrflag_syserr = (muerr_sys_list is not None)

    keyname_row = f"{VARNAME_ROW}:"
    # xxx mark varlist = f"{VARNAME_ROW} zCMB zHEL {VARNAME_MU} {VARNAME_MUERR}"
    varlist = f"{VARNAME_ROW} {VARNAME_zHD} {VARNAME_zHEL} " \
              f"{VARNAME_MU} {VARNAME_MUERR}"

    name_list   = base[VARNAME_ROW].to_numpy()
    zHD_list    = base[VARNAME_zHD].to_numpy()
    mu_list     = base[VARNAME_MU].to_numpy()
    muerr_list  = base[VARNAME_MUERR].to_numpy()

    if wrflag_nevt:
        varlist += f"  {VARNAME_NEVT_BIN}"
        nevt_list = base[VARNAME_NEVT_BIN].to_numpy()
    else:
        nevt_list = muerr_list  # anything to allow for loop with zip

    if wrflag_syserr:
        varlist += f" {VARNAME_MUERR_SYS}"
        syserr_list = muerr_sys_list
    else:
        syserr_list = muerr_list # anything to allow zip loop

    with open(path, "w") as f:
        write_HD_comments(f, unbinned, wrflag_syserr)
        f.write(f"VARNAMES: {varlist}\n")
        for (name, z, mu, muerr, nevt, syserr) in \
            zip(name_list, zHD_list, mu_list, muerr_list, 
                nevt_list, syserr_list):
            val_list = f"{name:<6}  {z:6.5f} {z:6.5f} {mu:8.5f} {muerr:8.5f} "
            if wrflag_nevt: val_list += f" {nevt} "
            if wrflag_syserr: val_list += f" {syserr:8.5f}"
            f.write(f"{keyname_row} {val_list}\n")
    return

    # end write_HD_binned

def write_HD_unbinned(path, base, muerr_sys_list):

    # Dec 2020
    # Write standard unbinned HD format for BBC method
    # Sep 30 2021: replace csv format with SNANA fitres format
    # Apr 30 2022: pass muerr_sys_list

    #if "CID" in df.columns:

    unbinned = True
    logging.info(f"Write unbinned HD to {path}")
    
    #print(f"\n xxx base=\n{base} \n")

    varname_row = "CID"
    keyname_row = "SN:"

    varlist = f"{varname_row} {VARNAME_IDSURVEY} " \
              f"{VARNAME_zHD} {VARNAME_zHEL} " \
              f"{VARNAME_MU} {VARNAME_MUERR}"
              # xxx mark delete Nov 9 2022 f"zCMB zHEL " \

    name_list   = base[VARNAME_CID].to_numpy()
    idsurv_list = base[VARNAME_IDSURVEY].to_numpy()
    zHD_list    = base[VARNAME_zHD].to_numpy()
    zHEL_list   = base[VARNAME_zHEL].to_numpy()
    mu_list     = base[VARNAME_MU].to_numpy()
    muerr_list  = base[VARNAME_MUERR].to_numpy()

    # check for optional quantities that may not exist in older files
    found_muerr_vpec = VARNAME_MUERR_VPEC in base
    found_muerr_sys  = muerr_sys_list is not None

    if found_muerr_vpec :   
        varlist += f" {VARNAME_MUERR_VPEC}"
        muerr2_list = base[VARNAME_MUERR_VPEC].to_numpy()
    else:
        muerr2_list = muerr_list # anything for zip command

    if found_muerr_sys:
        varlist += f" {VARNAME_MUERR_SYS}"
        syserr_list = muerr_sys_list
    else:
        syserr_list = muerr_list # anything for zip command

    # - - - - - - -
    with open(path, "w") as f:
        write_HD_comments(f, unbinned, found_muerr_sys)
        f.write(f"VARNAMES: {varlist}\n")
        for (name, idsurv, zHD, zHEL, mu, muerr, muerr2, syserr) in \
            zip(name_list, idsurv_list, zHD_list, zHEL_list,
                mu_list, muerr_list, muerr2_list, syserr_list ):
            val_list = f"{name:<10} {idsurv:3d} " \
                       f"{zHD:6.5f} {zHEL:6.5f} " \
                       f"{mu:8.5f} {muerr:8.5f}"
            if found_muerr_vpec: val_list += f" {muerr2:8.5f}"
            if found_muerr_sys:  val_list += f" {syserr:8.5f}"

            f.write(f"{keyname_row} {val_list}\n")
    return
    # end write_HD_unbinned

def write_HD_comments(f, unbinned, wrflag_syserr):

    f.write(f"# zHD       = redshift in CMB frame with VPEC correction\n")

    if unbinned:
        txt_zHEL = "helio redshift (beware: no VPEC corr)"
    else:
        txt_zHEL = "zHD"
    f.write(f"# zHEL      = {txt_zHEL} \n")

    f.write(f"# MU        = distance modulus corrected for bias and " \
            "contamination\n")
    f.write(f"# MUERR     = stat-uncertainty on MU \n")

    if wrflag_syserr:
        f.write(f"# MUERR_SYS = sqrt(COVSYS_DIAG) for 'ALL' sys " \
                "(diagnostic)\n")

    # write ISDATA flag as comment in HD 
    ISDATA = config[KEYNAME_ISDATA]
    f.write(f"# {KEYNAME_ISDATA}: {ISDATA}   "\
            "# flag for cosmology fitter to choose blind option\n")
    f.write(f"#\n")

    return
    # end write_HD_comments

def write_covariance(path, cov, opt_cov):
    
    add_labels     = (opt_cov == 1) # label some elements for human readability
    file_base      = os.path.basename(path)
    covdet         = np.linalg.det(cov)
    nrow           = cov.shape[0]

    logging.info("")
    logging.info(f"Write cov to {path}")

    # RK - write diagnostic to check if anything changes
    logging.info(f"    {file_base}: size={nrow}  |cov| = {covdet:.5e}")
    sys.stdout.flush() 

    # - - - - -
    # Write out the matrix
    nwr = 0 ; rownum = -1; colnum=-1

    if '.gz' in str(path):
        f = gzip.open(path,"wt")
    else:
        f = open(path,"wt") 

    f.write(f"{nrow}\n")
    for c in cov.flatten():
        nwr += 1
        is_new_row = False
        if (nwr-1) % nrow == 0 : 
            is_new_row = True ; rownum += 1 ; colnum = -1
        colnum += 1

        label = ""
        if add_labels:            
            if colnum == 0 or colnum == rownum : 
                label = f"# ({rownum},{colnum})"
        f.write(f"{c:13.6e}  {label}\n")
        # xxx mark delete f.write(f"{c:12.8f}  {label}\n")

    f.close()

    # end write_covariance

def write_summary_output(config, covariances, base):

    # write information to INFO.YAML that is intended to be 
    # picked up by cosmology fitting progam. Info includes
    # each covsys label & file, name of hubble diagram file.
    # and ISDATA_REAL flag.
    # Mar 2023: include VERSION_PHOTOMETRY and COSPAR_BIASCOR

    out  = Path(config["OUTDIR"])
    info = {} # init dictionary to dump to info file

    info['HD'] = HD_FILENAME

    cov_info = {}
    for i, (label, cov) in enumerate(covariances):
        covsys_file = get_covsys_filename(i)
        cov_info[i] = f"{label:<20} {covsys_file}"
        # xxx mark delete oct 13 2022 cov_info[i] = label

    info["COVOPTS"] = cov_info
        
    # xxx mark del Mar 2023:  info[KEYNAME_ISDATA] = config[KEYNAME_ISDATA]

    SNANA_VERSION = get_snana_version()
    info['SNANA_VERSION'] = SNANA_VERSION

    sim_version = None
    for key,item in config['header_info'].items() :
        info[key] = item     # store stuff from BBC table 
        if 'BIASCOR' in key :  # fetch biasCor sim version
            sim_version    = item[0]
 
    logging.info(f"Write {INFO_YML_FILENAME}")
    with open(out / INFO_YML_FILENAME, "w") as f:
        yaml.safe_dump(info, f )

    # - - - - - - - - - - - - - 
    # append cospar_biascor so that it's at the end, rather than 
    # at the beginning with default alphabetical ordering
    # Beware that this method of fetching cospar_biascor requries
    # original biasCor sim folder to still exist on disk; if this folder
    # is purged (e.g., to deal with quota crisis) then cospar will not
    # be found using this method.

    cospar_biascor = []
    if sim_version is not None:
        cospar_biascor = get_cospar_sim(sim_version)
        with open(out / INFO_YML_FILENAME, "at") as f:
            info_cospar = { 'COSPAR_BIASCOR': cospar_biascor }
            yaml.safe_dump(info_cospar, f )

    # end write_summary_output

def get_cospar_sim(sim_version):

    # for input snana_folder 'sim_version', run snana.exe GETINFO folder 
    # to extract name of README file, then parse README to get cosmo params

    cospar_sim = {}  # init output dictionary
    msgerr     = '\n'

    cmd = f"snana.exe GETINFO {sim_version}"
    
    ret = subprocess.run( [cmd], shell=True,
                          capture_output=True, text=True )
    ret_stdout = ret.stdout.split()

    key_readme = "README_FILE:"
    if 'FATAL' in ret_stdout:
        msgerr  = f"Cannot find biasCor sim-data folder {sim_version}\n"
        msgerr += f"\t May need to regenerate biasCor {sim_version}  \n"
        msgerr += f"\t Only need to recreate README, so try --faster with submit_batch_jobs.sh\n"
        logging.warning(msgerr)
        return cospar_sim
        #assert False, msgerr

    if key_readme not in ret_stdout:
        msgerr += f"  Cannot find {key_readme} key from command\n"
        msgerr += f"     {cmd}\n"
        assert False, msgerr

    k           = ret_stdout.index(key_readme)
    readme_file = ret_stdout[k+1]

    readme_contents = read_yaml(readme_file)
    docana          = readme_contents[KEY_DOCANA] # DOCUMENTATION block
    
    if KEY_SIM_INPUT in docana :
        sim_inputs      = docana[KEY_SIM_INPUT]
    else :
        sim_inputs = []  # allow for for older SNANA versions 
        logging.warning(f"did not find {KEY_SIM_INPUT} in biasCor sim README")

    for cospar in KEYLIST_COSPAR_SIM:
        if cospar in sim_inputs :
            cospar_sim[cospar] = sim_inputs[cospar]

    return cospar_sim

    # end get_cospar_sim

def write_correlation(path, label, base_cov, diag, base):
    logging.debug(f"\tWrite out cov for COVOPT {label}")

    zs = base[VARNAME_zHD].round(decimals=5)
    cov = diag + base_cov

    diag = np.sqrt(np.diag(cov))
    corr = cov / (diag[:, None] @ diag[None, :])

    np.fill_diagonal(corr, 1.0)
    np.savetxt(path, corr, fmt="%5.2f")
    np.savetxt(str(path).replace("corr", "cov"), cov, fmt="%9.6f")
    covdf = pd.DataFrame(cov, columns=zs, index=zs)
    base_covdf = pd.DataFrame(base_cov, columns=zs, index=zs)
    corr = pd.DataFrame((corr * 100).astype(int), columns=zs, index=zs)
    precision = pd.DataFrame(np.arcsinh(np.linalg.inv(cov)), columns=zs, index=zs)

    if corr.shape[0] < 1 :
        logging.debug("\tCreating precision and correlation matrix plots. Sit tight.")
        fig, axes = plt.subplots(figsize=(16, 14), ncols=2, nrows=2)
        axes = axes.flatten()
        annot = corr.shape[0] < 30

        args = {"shrink": 0.9}
        cm = np.nanmax(np.abs(cov))
        pm = np.nanmax(np.abs(precision.to_numpy()))
        bcm = np.nanmax(np.abs(base_covdf.to_numpy()))
        sb.heatmap(covdf.replace([np.inf, -np.inf], np.nan), annot=False, ax=axes[0], vmin=-cm, vmax=cm, cmap="PuOr", square=True, cbar_kws=args)
        sb.heatmap(precision, annot=False, ax=axes[1], cmap="PuOr", square=True, vmin=-pm, vmax=pm, cbar_kws=args)
        sb.heatmap(corr, annot=annot, fmt="d", ax=axes[2], cmap="RdBu", vmin=-100, vmax=100, square=True, cbar_kws=args)
        sb.heatmap(base_covdf.replace([np.inf, -np.inf], np.nan), annot=False, ax=axes[3], cmap="PuOr", vmin=-bcm, vmax=bcm, square=True, cbar_kws=args)
        axes[0].set_title("Covariance matrix")
        axes[1].set_title(f"Arcsinh(precision) ")
        axes[2].set_title("Correlation matrix (percent)")
        axes[3].set_title("Covariance matrix without diagonal contributions")
        plt.tight_layout()
        fig.savefig(path.with_suffix(".png"), bbox_inches="tight", dpi=300)
    else:
        logging.info("\tMatrix is large, skipping plotting")


def write_debug_output(config, covariances, base, summary):
    # Plot correlation matrix
    out = Path(config["OUTDIR"])

    # The slopes can be used to figure out what systematics have largest impact on cosmology
    logging.info("Write out summary.csv information")
    with open(out / "summary.csv", "w") as f:
        with pd.option_context("display.max_rows", 100000, "display.max_columns", 100, "display.width", 1000):
            f.write(summary.__repr__())

    logging.info("Showing correlation matrices:")
    diag = np.diag(base[VARNAME_MUERR] ** 2)
    for i, (label, cov) in enumerate(covariances):
        write_correlation(out / f"corr_{i}_{label}.txt", label, cov, diag, base)


def get_lcfit_info(submit_info):
    path = Path(submit_info["INPDIR_LIST"][0]) / "SUBMIT.INFO"
    logging.info(f"Loading LCFIT SUBMIT.INFO from {path}")
    return read_yaml(path)


def remove_nans(data):
    base_name = get_name_from_fitopt_muopt(f_REF, m_REF)
    base      = data[base_name]
    try:
        base = base.drop('FIELD', axis=1)
    except:
        #FIELD not in df  
        pass
    mask = ~(base.isnull().any(axis=1))
    for name, df in data.items():
        data[name] = df.loc[mask, :]
    return data, base.loc[mask, :]


def create_covariance(config, args):
    # Define all our pathing to be super obvious about where it all is
    input_dir       = Path(config["INPUT_DIR"])
    version         = config["VERSION"]
    data_dir        = input_dir / version
    extra_covs      = config.get("EXTRA_COVS",[])
    use_cosmomc     = config['use_cosmomc']

    # Read in all the needed data
    submit_info   = read_yaml(input_dir / "SUBMIT.INFO")

    # read optional sys scales 
    if KEYNAME_SYS_SCALE_FILE in config:
        sys_scale_file  = Path(config[KEYNAME_SYS_SCALE_FILE])
        sys_scale       = read_yaml(sys_scale_file)
    else:
        sys_scale = { 0: (None,1.0) }

    fitopt_scales = get_fitopt_scales(submit_info, sys_scale)

    # Also need to get the MUOPT labels from the original LCFIT directory
    muopt_labels = {int(x.replace("MUOPT", "")): l for x, l, _ in  \
                    submit_info.get("MUOPT_OUT_LIST", [])}

    #tack on extra covs as muopts
    extracovdict = {}
    for extra in extra_covs:
        muopt_labels[len(muopt_labels)+1] = extra.split()[0]
        extracovdict[extra.split()[0]] = extra.split()[2]


    # Feb 15 2021 RK - check option to selet only one of the muopt
    if args.muopt >= 0 :
        tmp_label    = muopt_labels[args.muopt]
        muopt_labels = { args.muopt : "DEFAULT" }     
    
    if 'MUOPT_SCALES' in config:
        muopt_scales = config["MUOPT_SCALES"]
    else:
        muopt_scales = { "DEFAULT" : 1.0 }

    for extra in extra_covs:
        muopt_scales[extra.split()[0]] = extra.split()[1]

    # Load in all the hubble diagrams
    logging.info(f"Read Hubble diagrams for version = {version}")
    data = get_hubble_diagrams(data_dir, args, config)

    # Filter data to remove rows with infinite error
    data, base = remove_nans(data)

    # Now that we have the data, figure out how each much each
    # FITOPT/MUOPT pair contributes to cov
    contributions, summary = get_contributions(data, fitopt_scales,
                                               muopt_labels, 
                                               muopt_scales, 
                                               extracovdict)

    # find contributions which match to construct covs for each COVOPT
    logging.info(f"Compute covariance for COVOPTS")

    # Add covopt to compute everything
    if args.nosys:
        covopts_default = ["[=DEFAULT] [,=DEFAULT]"] # Dec 2022 
    else:
        covopts_default = ["[ALL] [,]"]

    covopts = covopts_default + config.get("COVOPTS",[])  

    covariances = \
        [ get_cov_from_covopt(c, contributions, base, 
                              config.get("CALIBRATORS")) for c in covopts]
    
    # P. Armstrong 05 Aug 2022
    # Create hubble_diagram.txt for every systematic, not just nominal
    if args.systematic_HD:
        label_list = list(data.keys())
    # Only create hubble_diagram.txt for the nominal
    else:
        label_list = [get_name_from_fitopt_muopt(f_REF, m_REF)]

    # write standard output for cov(s) and hubble diagram (9.22.2021)
    write_standard_output(config, args, covariances, data, label_list)

    # write specialized output for cosmoMC sampler
    if use_cosmomc :
        write_cosmomc_output(config, args, covariances, base)

    write_summary_output(config, covariances, base)

    # Sep 30 2022 RK - submit_batch_jobs needs output yaml to communicate
    if args.yaml_file is not None:
        write_yaml(args, len(covariances) )

    # end create_covariance

def prep_config(config,args):

    # Dec 7 2022: fix bug setting override args at start of method instead 
    #   of at the end.

    # - - - - -
    # check override args (RK, Feb 15 2021)
    if args.input_dir is not None :
        config["INPUT_DIR"] = args.input_dir
        logging.info(f"OPTION: override INPUT_DIR with {args.input_dir}")

    if args.outdir is not None :
        config["OUTDIR"] = args.outdir
        logging.info(f"OPTION: override OUTDIR with {args.outdir}")        

    if args.sys_scale_file is not None :
        config[KEYNAME_SYS_SCALE_FILE] = args.sys_scale_file
        logging.info(f"OPTION: override {KEYNAME_SYS_SCALE_FILE} with " \
                     f"{args.sys_scale_file}")

    if args.version is not None:
        config["VERSION"] = args.version
        logging.info(f"OPTION: override VERSION with {args.version}")

    if args.muopt >= 0 :
        global m_REF ; m_REF = args.muopt  # RK, Feb 2021
        logging.info(f"OPTION: use only MUOPT{m_REF:03d}")        

    # - - -  - -    

    key_legacy_list = [ 'COSMOMC_TEMPLATES',      'DATASET_FILE', 
                        'SYSFILE' ]
    key_update_list = [ 'COSMOMC_TEMPLATES_PATH', 'COSMOMC_DATASET_FILE',
                        'SYS_SCALE_FILE' ]

    for key_legacy,key_update in zip(key_legacy_list,key_update_list):
        if key_legacy in config:
            msg = f"Replace legacy key '{key_legacy}' with {key_update}"
            logging.info(msg)
            config[key_update] = config[key_legacy] 

    path_list = [ 'INPUT_DIR', 'OUTDIR', 'SYS_SCALE_FILE', 
                  'COSMOMC_TEMPLATES_PATH' ]
    for path in path_list:
        if path in config:
            config[path] = os.path.expandvars(config[path]) ;   

    # check special/legacy features for cosmoMC/JLA
    config['use_cosmomc'] = False
    if 'COSMOMC_TEMPLATES_PATH' in config: 
        config['use_cosmomc'] = True

    # WARNING: later add option to read from input file
    config['nbin_x1'] = args.nbin_x1
    config['nbin_c']  = args.nbin_c

        
    # end prep_config

def write_yaml(args, n_cov):

    start_time = args.start_time
    end_time   = time.time()
    cpu_min    = (end_time - start_time)/60.0

    with open(args.yaml_file,"wt") as y:
        y.write(f"N_COVMAT:       {n_cov}\n")
        y.write(f"ABORT_IF_ZERO:  {n_cov}    # same as N_COVMAT\n")
        y.write(f"CPU_MINUTES:    {cpu_min:.3f} \n")

    return
    # end write_yaml

# ===================================================
if __name__ == "__main__":


    try:
        setup_logging()
        logging.info("# ========== BEGIN create_covariance ===============")

        command = " ".join(sys.argv)
        logging.info(f"# Command: {command}")
        
        args   = get_args()
        args.start_time = time.time()
        config = read_yaml(args.input_file)
        prep_config(config,args)  # expand vars, set defaults, etc ...
        create_covariance(config, args)

    except Exception as e:
        logging.exception(e)
        raise e
