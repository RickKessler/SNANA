#!/usr/bin/env python
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
#
# ===============================================

import os, argparse, logging, shutil
import re, yaml, sys, gzip, math
import numpy  as np
import pandas as pd
from pathlib import Path
from functools import reduce
from sklearn.linear_model import LinearRegression
import seaborn as sb
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.cosmology import Planck13, z_at_value
from astropy.cosmology import FlatLambdaCDM

SUFFIX_M0DIF  = "M0DIF"
SUFFIX_FITRES = "FITRES"

PREFIX_COVSYS  = "covsys"
HD_FILENAME    = "hubble_diagram.txt"

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

VARNAME_iz     = "IZBIN"
VARNAME_z      = "z"  # note that zHD is internally renamed z
VARNAME_x1     = "x1"
VARNAME_c      = "c"

VARNAME_NEVT_BIN = 'NEVT'

SUBDIR_COSMOMC = "cosmomc"
KEYNAME_ISDATA = 'ISDATA_REAL'   # key in fitres of M0DIF file from SALT2mu

KEYNAME_SYS_SCALE_FILE = "SYS_SCALE_FILE"

m_REF = 0  # MUOPT reference number for cov
f_REF = 0  # FITOPT reference number for cov

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

    # parse it
    args = parser.parse_args()
    if args.subtract_vpec: args.unbinned = True

    #if len(sys.argv) == 1:
    #    parser.print_help()
    #    sys.exit()

    if args.HELP : print_help_menu()

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

#  [CosmoMC-label]  [FITOPT-label, MUOPT-label]
COVOPTS:
- '[NOSYS]  [=DEFAULT,=DEFAULT]' # no syst (for JLA method only)
- '[PS1CAL] [+PS1CAL,=DEFAULT]'  # PS1CAL syst only MUOPT=0
- '[ALLCAL] [+CAL,=DEFAULT]'     # All syst with 'CAL' in name, MUOPT=0
- '[SALT2]  [+SALT2,=DEFAULT]'   # SALT2 syst only, MUOPT=0
- '[NOCAL]  [-CAL,=DEFAULT]'     # all syst except for calib, MUOPT=0
- '[SCAT]   [=DEFAULT,+SCAT]     # FITOPT=0, MUOPTs with SCAT in label

MUOPT_SCALES:  # replace scales in {KEYNAME_SYS_SCALE_FILE}
  CLAS_SNIRF:  1.0
  SCATTER_C11: 1.0

OUTDIR:  <name of output directory>

# specialized keys for CosmoMC
COSMOMC_TEMPLATES_PATH: <out path of ini files for Cosmomc>
  [allows legacy key COSMOMC_TEMPLATES:]

COSMOMC_DATASET_FILE: <out file with cosmomc instructions>
  [allows legacy key DATASET_FILE: ]

    """
    print(f"{menu}")

    sys.exit()

    # end print_help_menu


def check_isdata_real(hd_file):

    # Created oct 6 2021 by R.Kessler
    # read first few lines of hd_file and look for ISDATA_REAL  key
    # that is in a comment field ... hence read as yaml.
    # SALT2mu begain wriring ISDATA_REAL key on Oct 6 2021;
    # Function returns ISDATA_REAL = 0 or 1 of key is found;
    # returns -1 (unknown) if key not found.

    with gzip.open(hd_file, 'r') as f:
        line_list = f.readlines()

    isdata_real = -1  # init to unknonw
    maxline_read = 10 # bail after this many lines
    nline_read   = 0

    key = f"{KEYNAME_ISDATA}:" # include colon in brute force search

    for line in line_list:
        line  = line.rstrip()  # remove trailing space and linefeed  
        line  = line.decode('utf-8')
        wd_list = line.split()
        if key in wd_list:
            j = wd_list.index(key)
            isdata_real = int(wd_list[j+1])
                        
        nline_read += 1
        if nline_read == maxline_read or isdata_real>=0 : break

    logging.info(f"ISDATA_REAL = {isdata_real}")
    return isdata_real

    # end check_isdata_real

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
        df = df.sort_values(["zHD", "CID"])
        df = df.rename(columns={"zHD": VARNAME_z, "MUMODEL": VARNAME_MUREF})
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

    elif VARNAME_z in df.columns:
        # should not need z-sorting here for M0DIF, but what the heck
        df = df.sort_values(VARNAME_z)

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
    isdata_real = -1  # init to unknwon data or sim

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
            if first_load:  isdata_real = check_isdata_real(hd_file)
            first_load = False

    #sys.exit(f"\n xxx\n result = {result} \n")

    config[KEYNAME_ISDATA] = isdata_real # Oct 6 2021, R.Kessler

    # - - - - -  -
    # for unbinned or rebin, select SNe that are common to all 
    # FITOPTs and MUOPTs
    if is_unbinned or is_rebin:
        HD_list = get_common_set_of_sne(HD_list)
        
    # - - - - -    
    # df['e'] = e.values or  df1 = df1.assign(e=e.values)

    if is_rebin : 
        label0 = label_list[0]
        get_rebin_info(config,HD_list[label])
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

    zmin  = HD[VARNAME_z].min()   - epsilon 
    zmax  = HD[VARNAME_z].max()   + epsilon 
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
    col_z        = HD_unbinned[VARNAME_z]
    col_mures    = HD_unbinned[VARNAME_MURES]
    col_mudif    = HD_unbinned[VARNAME_M0DIF]
    col_muref    = HD_unbinned[VARNAME_MUREF]
    col_mu       = HD_unbinned[VARNAME_MU]
    col_muerr    = HD_unbinned[VARNAME_MUERR]
    col_muerr_renorm = HD_unbinned[VARNAME_MUERR_RENORM]
    HD_rebin_dict = { VARNAME_ROW: [], 
                      VARNAME_z: [], VARNAME_MU: [] , VARNAME_MUERR: [],
                      VARNAME_MUREF: [], VARNAME_NEVT_BIN: [] }

    wgt0      = 1.0/(col_muerr*col_muerr)
    wgt1      = 1.0/(col_muerr_renorm*col_muerr_renorm)
    wgtmuref  = wgt0 * col_muref
    wgtmudif  = wgt1 * col_mudif
    wgtmu     = wgt1 * col_mu

    # - - - - -
    OM_ref = 0.30  # .xyz should read OM_ref from BBC output
    zcalc_grid, mucalc_grid = get_HDcalc(OM_ref)
    # .xyz

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
        HD_rebin_dict[VARNAME_z].append(z_invert)
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
    assert df1[VARNAME_MU].size == df2[VARNAME_MU].size, \
        "Oh no, looks like you have a different number of bins/supernova " \
        "for your systematic and this is not good."

    diff = scale * ((df1[VARNAME_MU] - df1[VARNAME_MUREF]) - \
                    (df2[VARNAME_MU] - df2[VARNAME_MUREF])).to_numpy()
    diff[~np.isfinite(diff)] = 0
    cov = diff[:, None] @ diff[None, :]

    # Determine the gradient using simple linear regression
    reg = LinearRegression()
    weights = 1 / np.sqrt(0.003**2 + df1[VARNAME_MUERR].to_numpy()**2 \
                          + df2[VARNAME_MUERR].to_numpy()**2)
    mask = np.isfinite(weights)
    reg.fit(df1.loc[mask, [VARNAME_z]], diff[mask], 
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
        if len(np.argwhere(data['CIDstr'] == row['CID1'])) == 0:
            print(row['CID1'],'1 missing from output/cosmomc/data_wCID.txt')
            continue
        if len(np.argwhere(data['CIDstr'] == row['CID2'])) == 0:
            print(row['CID2'],'2 missing from output/cosmomc/data_wCID.txt')
            continue
        ww1 = np.argwhere(data['CIDstr'] == row['CID1'])[0][0]
        ww2 = np.argwhere(data['CIDstr'] == row['CID2'])[0][0]
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
    # Covopts will come in looking like "[cal] [+cal,=DEFAULT]"
    # We have to parse this. Eventually can make this structured and 
    # move away from legacy, but dont want to make too many people 
    # change how they are doing things in one go

    # split covopt into two terms so that extra pad spaces 
    # doesn't break findall command (RK May 14 2021)  
    tmp0 = covopt.split()[0]
    tmp1 = covopt.split()[1]
    label                       = re.findall(r"\[(.*)\]",      tmp0)[0]
    fitopt_filter, muopt_filter = re.findall(r"\[(.*),(.*)\]", tmp1)[0]

    fitopt_filter = fitopt_filter.strip()
    muopt_filter  = muopt_filter.strip()
    logging.debug(f"Computing COV({label}) with FITOPT filter " \
                  f"'{fitopt_filter}' and MUOPT filter '{muopt_filter}'")

    final_cov = None

    if calibrators:
        mask_calib = base.reset_index()["CID"].isin(calibrators)

    for key, cov in contributions.items():
        fitopt_label, muopt_label = key.split("|")

        if apply_filter(fitopt_label, fitopt_filter) and \
           apply_filter(muopt_label, muopt_filter):
            if final_cov is None:
                final_cov = cov.copy()*0
            if True:
                # If we have calibrators and this is VPEC term, filter out calib
                if calibrators and (apply_filter(fitopt_label, "+VPEC") or apply_filter(muopt_label, "+VPEC") or \
                                    apply_filter(fitopt_label, "+ZSHIFT") or apply_filter(muopt_label, "+ZSHIFT")):
                    print(f'FITOPT {fitopt_label} MUOPT {muopt_label} ignored for calibrators...')
                    cov2 = cov.copy()
                    cov2[mask_calib, :] = 0
                    cov2[:, mask_calib] = 0
                    final_cov += cov2
                else:
                    final_cov += cov

    assert final_cov is not None, f"No systematics matched COVOPT {label} with " \
        f"FITOPT filter '{fitopt_filter}' and MUOPT filter '{muopt_filter}'!"

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
            logging.info(f"WARNING: {label} Matrix is not Positive-Definite")
        
        # check that COV is well conditioned to deal with float precision
        epsilon = sys.float_info.epsilon
        cond = np.linalg.cond(effective_cov)
        assert cond < 1 / epsilon, "Cov matrix is ill-conditioned and cannot be inverted"
        logging.info(f"Covar condition for COVOPT {label} is {cond:.3f}")

        # May 2 2022 R.Kessler - mark delete here to save time re-inverting.
        # Finally, re-invert the precision matrix and ensure its within 
        # tolerance of the original covariance
        # xxx mark delete cov2 = np.linalg.inv(precision)
        # xxx mark delete assert np.all(np.isclose(effective_cov, cov2)), 
        # xxx mark delete "Double inversion does not give original covariance, matrix is unstable"

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

def write_standard_output(config, unbinned, covs, base):
    # Created 9.22.2021 by R.Kessler
    # Write standard cov matrices and HD for cosmology fitting programs;
    # e.g., wfit, CosmoSIS, firecrown ...
    # Note that CosmoMC uses a more specialized output created
    # by write_cosmomc_output().

    logging.info("")
    logging.info("   OUTPUT  ")

    outdir = Path(config["OUTDIR"])
    os.makedirs(outdir, exist_ok=True)

    data_file = outdir / HD_FILENAME

    # Apr 30 2022: get array of muerr_sys(ALL) for output
    muerr_sys_list = get_muerr_sys(covs)
            
    # - - - -
    if unbinned :
        write_HD_unbinned(data_file, base, muerr_sys_list)
    else:
        write_HD_binned(data_file, base, muerr_sys_list)

    # Create covariance matrices and datasets
    opt_cov = 1  # tag rows and diagonal elements
    for i, (label, cov) in enumerate(covs):
        base_file   = f"{PREFIX_COVSYS}_{i:03d}.txt" 
        base_file  += '.gz'  # force gzip file, Apr 22 2022
        covsys_file = outdir / base_file
        write_covariance(covsys_file, cov, opt_cov)

    return

    # end write_standard_output

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
            varlist = [VARNAME_CID, VARNAME_IDSURVEY, VARNAME_z, VARNAME_MU, VARNAME_MUERR]
        else:
            varlist = [VARNAME_z, VARNAME_MU, VARNAME_MUERR]

        base[varlist].to_csv(path, sep=" ", index=False, float_format="%.5f")
        return

    zs   = base[VARNAME_z].to_numpy()
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

    logging.info(f"Write binned HD to {path}")

    wrflag_nevt   = (VARNAME_NEVT_BIN in base)
    wrflag_syserr = (muerr_sys_list is not None)

    keyname_row = f"{VARNAME_ROW}:"
    varlist = f"{VARNAME_ROW} zCMB zHEL {VARNAME_MU} {VARNAME_MUERR}"

    name_list   = base[VARNAME_ROW].to_numpy()
    z_list      = base[VARNAME_z].to_numpy()
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
        write_HD_comments(f,wrflag_syserr)
        f.write(f"VARNAMES: {varlist}\n")
        for (name, z, mu, muerr, nevt, syserr) in \
            zip(name_list, z_list, mu_list, muerr_list, 
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

    logging.info(f"Write unbinned HD to {path}")
    
    #print(f"\n xxx base=\n{base} \n")

    varname_row = "CID"
    keyname_row = "SN:"

    varlist = f"{varname_row} {VARNAME_IDSURVEY} " \
              f"zCMB zHEL " \
              f"{VARNAME_MU} {VARNAME_MUERR}"

    name_list   = base[VARNAME_CID].to_numpy()
    idsurv_list = base[VARNAME_IDSURVEY].to_numpy()
    z_list      = base[VARNAME_z].to_numpy()
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
        write_HD_comments(f,found_muerr_sys)
        f.write(f"VARNAMES: {varlist}\n")
        for (name, idsurv, z, mu, muerr, muerr2, syserr) in \
            zip(name_list, idsurv_list, z_list, 
                mu_list, muerr_list, muerr2_list, syserr_list ):
            val_list = f"{name:<10} {idsurv:3d} " \
                       f"{z:6.5f} {z:6.5f} " \
                       f"{mu:8.5f} {muerr:8.5f}"
            if found_muerr_vpec: val_list += f" {muerr2:8.5f}"
            if found_muerr_sys:  val_list += f" {syserr:8.5f}"

            f.write(f"{keyname_row} {val_list}\n")
    return
    # end write_HD_unbinned

def write_HD_comments(f,wrflag_syserr):
    f.write(f"# MU        = distance modulus corrected for bias and " \
            "contamination\n")
    f.write(f"# MUERR     = stat-uncertainty on MU \n")

    if wrflag_syserr:
        f.write(f"# MUERR_SYS = sqrt(COVSYS_DIAG) for 'ALL' sys " \
                "(diagnostic)\n")

    f.write(f"#\n")
    return
    # end write_HD_comments

def write_covariance(path, cov, opt_cov):

    add_labels     = (opt_cov == 1) # label some elements for human readability
    file_base      = os.path.basename(path)
    covdet         = np.linalg.det(cov)
    nrow           = cov.shape[0]

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
        f.write(f"{c:12.8f}  {label}\n")

    f.close()

    # end write_covariance

def write_summary_output(config, covariances, base):
    out = Path(config["OUTDIR"])
    info = {}
    cov_info = {}
    for i, (label, cov) in enumerate(covariances):
        cov_info[i] = label
    info["COVOPTS"] = cov_info

    info[KEYNAME_ISDATA] = config[KEYNAME_ISDATA]

    logging.info("Write INFO.YML")
    with open(out / "INFO.YML", "w") as f:
        yaml.safe_dump(info, f)

    # end write_summary_output

def write_correlation(path, label, base_cov, diag, base):
    logging.debug(f"\tWrite out cov for COVOPT {label}")

    zs = base[VARNAME_z].round(decimals=5)
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
        sys_scale     = read_yaml(sys_scale_file)
    else:
        sys_scale = {0: (None,1.0) }

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
    covopts = ["[ALL] [,]"] + config.get("COVOPTS",[])  

    covariances = [get_cov_from_covopt(c, contributions, base, 
                                       config.get("CALIBRATORS")) for c in covopts]

    # write standard output (9.22.2021)
    write_standard_output(config, args.unbinned, covariances, base)

    # write specialized output for cosmoMC sampler
    if use_cosmomc :
        write_cosmomc_output(config, args, covariances, base)

    write_summary_output(config, covariances, base)

    # xxxxx mark delete Sep 2021 R.K
    #  write_debug_output(config, covariances, base, summary)
    # xxxxxxxxxx

    # end create_covariance

def prep_config(config,args):

    path_list = [ 'INPUT_DIR', 'OUTDIR', 'SYS_SCALE_FILE', 
                  'COSMOMC_TEMPLATES_PATH' ]
    
    # 9.22.2021 RK - check legacy keys
    key_legacy_list = [ 'COSMOMC_TEMPLATES',      'DATASET_FILE', 
                        'SYSFILE' ]
    key_update_list = [ 'COSMOMC_TEMPLATES_PATH', 'COSMOMC_DATASET_FILE',
                        'SYS_SCALE_FILE' ]

    for key_legacy,key_update in zip(key_legacy_list,key_update_list):
        if key_legacy in config:
            msg = f"Replace legacy key '{key_legacy}' with {key_update}"
            logging.info(msg)
            config[key_update] = config[key_legacy] 

    for path in path_list:
        if path in config:
            config[path] = os.path.expandvars(config[path]) ;   

    # check special/legacy features for cosmoMC/JLA
    config['use_cosmomc'] = False
    if 'COSMOMC_TEMPLATES_PATH' in config: 
        config['use_cosmomc'] = True

    # WARNING: later add option to read from input file
    #sys.exit(f" xxx nbin(x1,c) = {args.nbin_x1} {args.nbin_c} ")
    config['nbin_x1'] = args.nbin_x1
    config['nbin_c']  = args.nbin_c

    # check override args (RK, Feb 15 2021)
    if args.input_dir is not None :
        config["INPUT_DIR"] = args.input_dir
        logging.info(f"OPTION: override INPUT_DIR with {args.input_dir}")

    if args.outdir is not None :
        config["OUTDIR"] = args.outdir
        logging.info(f"OPTION: override OUTDIR with {args.outdir}")        

    if args.version is not None:
        config["VERSION"] = args.version
        logging.info(f"OPTION: override VERSION with {args.version}")

    if args.muopt >= 0 :
        global m_REF ; m_REF = args.muopt  # RK, Feb 2021
        logging.info(f"OPTION: use only MUOPT{m_REF:03d}")        
        
    # end prep_config

# ===================================================
if __name__ == "__main__":
    try:
        setup_logging()
        logging.info("# ========== BEGIN create_covariance ===============")
        args   = get_args()
        config = read_yaml(args.input_file)
        prep_config(config,args)  # expand vars, set defaults, etc ...
        create_covariance(config, args)
    except Exception as e:
        logging.exception(e)
        raise e
