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
# ===============================================

import argparse
import logging
import shutil
from functools import reduce
import numpy as np
import os
import pandas as pd
from pathlib import Path
import re
import yaml
import sys
from sklearn.linear_model import LinearRegression
import seaborn as sb
import matplotlib.pyplot as plt


SUFFIX_M0DIF  = "M0DIF"
SUFFIX_FITRES = "FITRES"

KEYNAME_COSMOMC_METHOD = 'COSMOMC_METHOD'
COSMOMC_METHOD_JLA = "JLA"
COSMOMC_METHOD_BBC = "BBC"

VARNAME_MU        = "MU"
VARNAME_M0DIF     = "M0DIF"
VARNAME_MUDIF     = "MUDIF"
VARNAME_MUDIFERR  = "MUDIFERR"
VARNAME_MUREF     = "MUREF"
VARNAME_MURES     = "MURES"
VARNAME_MUERR     = "MUERR"
VARNAME_iz     = "IZBIN"
VARNAME_z      = "z"  # note that zHD is internally renamed z
VARNAME_x1     = "x1"
VARNAME_c      = "c"

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
    parser.add_argument("input_file", help=msg, nargs="?", default=None)

    msg = "override INPUT_DIR" 
    parser.add_argument("-i", "--input_dir", help=msg, nargs='?',
                        type=str, default=None)
    
    msg = "override OUTDIR" 
    parser.add_argument("-o", "--outdir", help=msg, nargs='?', type=str, default=None)

    msg = "override COSMOMC_METHOD" 
    parser.add_argument("--method", help=msg, nargs='?', type=str, default=None)
    
    msg = "override data VERSION"
    parser.add_argument("-v", "--version", help=msg, nargs='?', type=str, default=None)
    
    msg = "use only this muopt number"
    parser.add_argument("-m", "--muopt", help=msg, nargs='?', type=int, default=-1 ) 
    
    msg = "Use each SN instead of BBC binning"
    parser.add_argument("-u", "--unbinned", help=msg, action="store_true")

    msg = "number of x1 bins (default=1)"
    parser.add_argument("--nbin_x1", help=msg, nargs='?', type=int, default=1 )

    msg = "number of c bins (default=1)"
    parser.add_argument("--nbin_c", help=msg, nargs='?', type=int, default=1 )

    #msg = "rebin args; 0 (unbinned) or e.g., c:5,x1:2 (5 c bins,2 x1 bins)"
    #parser.add_argument("--rebin", help=msg, nargs='+', type=str )

    msg = "Subtract MUERR(VPEC) from MUERR. Forces unbinned."
    parser.add_argument("-s", "--subtract_vpec", help=msg, action="store_true")
    args = parser.parse_args()
    if args.subtract_vpec: args.unbinned = True

    #if len(sys.argv) == 1:
    #    parser.print_help()
    #    sys.exit()

    if args.HELP : print_help_menu()

    return args


def print_help_menu():
    menu = """

   HELP MENU for create_covariance CONFIG
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

COSMOMC_TEMPLATES: <location of ini files for Cosmomc>

DATASET_FILE: <file name with cosmomc instructions>

COSMOMC_METHOD: JLA   # default: legacy option with old CosmoMC code
   or
COSMOMC_METHOD: BBC   # BBC option with new CosmoMC code

INPUT_DIR:  <BBC output dir from Pippin>
VERSION:    <subDir under INPUT_DIR>  

SYSFILE: <yaml file with systematic scales>

NAME:      <name> # ???

#  [CosmoMC-label]  [FITOPT-label,MUOPT-label]
COVOPTS:
- '[NOSYS]  [=DEFAULT,=DEFAULT]' # no syst (for JLA method only)
- '[PS1CAL] [+PS1CAL,=DEFAULT]'  # PS1CAL syst only MUOPT=0
- '[ALLCAL] [+CAL,=DEFAULT]'     # All syst with 'CAL' in name, MUOPT=0
- '[SALT2]  [+SALT2,=DEFAULT]'   # SALT2 syst only, MUOPT=0
- '[NOCAL]  [-CAL,=DEFAULT]'     # all syst except for calib, MUOPT=0
- '[SCAT]   [=DEFAULT,+SCAT]     # FITOPT=0, MUOPTs with SCAT in label

MUOPT_SCALES:  # replace scales in SYSFILE
  CLAS_SNIRF:  1.0
  SCATTER_C11: 1.0

OUTDIR:  <name of output directory>

    """
    print(f"{menu}")

    sys.exit()

def load_hubble_diagram(path, args, config):

    # read single M0DIF or FITRES file from BBC output,
    # and return contents.

    if not os.path.exists(path):
        raise ValueError(f"Cannot load data from {path} - it doesnt exist")

    df = pd.read_csv(path, delim_whitespace=True, comment="#")
    logging.debug(f"\tLoaded data with Nrow x Ncol {df.shape} from {path}")

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

    # Sort by CID for unbiined; sort by z for M0DIF
    # --> ensure direct subtraction comparison
    if "CID" in df.columns:
        df["CID"] = df["CID"].astype(str)
        df = df.sort_values(["zHD", "CID"])
        df = df.rename(columns={"zHD": VARNAME_z, "MUMODEL": VARNAME_MUREF})
        if args.subtract_vpec:
            msgerr = f"Cannot subtract VPEC because MUERR_VPEC " \
                     f"doesn't exist in {path}"
            assert "MUERR_VPEC" in df.columns, msgerr
            # xxx mark delete df[VARNAME_MUERR] = np.sqrt(df[VARNAME_MUERR] ** 2 - df['biasScale_muCOV']*df["MUERR_VPEC"] ** 2)
            df[VARNAME_MUERR] = np.sqrt(df[VARNAME_MUERR] ** 2 - df["MUERR_VPEC"] ** 2) # removed biasScale_muCOV
            logging.debug("Subtracted MUERR_VPEC from MUERR")
        #elif config.get("CALIBRATORS"):
        #    calib_mask = df["CID"].isin(config.get("CALIBRATORS"))
        #    df.loc[calib_mask, VARNAME_MUERR] = \
        #        np.sqrt(df.loc[calib_mask, VARNAME_MUERR] ** 2 - \
        #        df.loc[calib_mask, "MUERR_VPEC"] ** 2)
        df['CIDstr'] = df['CID'].astype(str)+"_"+df['IDSURVEY'].astype(str)
        df = df.set_index(["IDSURVEY", "CID"])
    elif VARNAME_z in df.columns:
        # should not be needed here for M0DIF, but what the heck
        df = df.sort_values(VARNAME_z)

    # - - - -  -
    # strip out only what is needed for CosmoMC; ignore the rest
    # .xyz BEWARE: do not strip for rebin
    # mark delete  df = df.loc[:, ["z", "MU", "MUERR", "MUREF"]]
    
    return df
    # end load_hubble_diagram


def get_hubble_diagrams(folder, args, config):

    # return table for each Hubble diagram 

    is_rebin    = config['nbin_x1'] > 1 or config['nbin_c'] > 1
    is_unbinned = args.unbinned 
    is_binned   = not (is_unbinned or is_rebin)

    folder_expand = Path(os.path.expandvars(folder))
    logging.debug(f"Loading all data files in {folder_expand}")
    HD_list = {}
    infile_list = []
    label_list  = []

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

            # grab contents of every M0DIF(binned) or FITRES(unbinned) file 
            HD_list[label] = load_hubble_diagram(folder_expand/infile,
                                                 args, config)


    #sys.exit(f"\n xxx\n result = {result} \n")

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
            HD_list[label]['iHD'] = config['col_iHD']
            HD_list[label] = rebin_hubble_diagram(config,HD_list[label])

    return HD_list


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

    logging.info(f"\t z (min,max) = {zmin} {zmax}   (nbin={nbin_z}) ")
    logging.info(f"\t x1(min,max) = {x1min} {x1max}  (nbin={nbin_x1}) ")
    logging.info(f"\t c (min,max) = {cmin} {cmax}  (nbin={nbin_c}) ") 
    
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
    config['col_ix1']        = col_ix1
    config['col_ic']         = col_ic
    config['col_iHD']        = col_iHD

    #sys.exit(f"\n xxx DEBUG DIE from get_rebins xxx ")

def rebin_hubble_diagram(config, HD_unbinned):

    # Dec 29 2020
    # rebin unbinned results based on nbin_xxx user inputs
    
    nbin_x1      = config['nbin_x1']
    nbin_c       = config['nbin_c']
    nbin_HD      = config['nbin_HD']
    col_iHD      = HD_unbinned['iHD']
    col_c        = HD_unbinned[VARNAME_c]
    col_z        = HD_unbinned[VARNAME_z]
    col_mures    = HD_unbinned[VARNAME_MURES]
    col_muerr    = HD_unbinned[VARNAME_MUERR]
    HD_rebin = {}

    wgt      = 1.0/(col_muerr*col_muerr)
    wgtmures = wgt * col_mures

    for i in range(0,nbin_HD):
        binmask  = (col_iHD == i)
        m0dif    = np.sum(wgtmures[binmask]) / np.sum(wgt[binmask])

        print(f"  xxx rebinned m0dif[{i:2d}] = {m0dif:7.4f}")

        if i == 1:
            print(f"\n xxx z   = \n{col_z[binmask]}\n")
            print(f"\n xxx wgt = \n{wgt[binmask]}\n")
            sys.exit(f"\n xxx DEBUG DIE xxx \n")
            
# .xyz

#    print(f"\n xxx HD_unbinned = \n{HD_unbinned}")

    sys.exit("\n xxx DEBUG DIE xxx \n")
    return HD_rebin

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
    """ Returns both the covariance contribution and summary stats (slope and mean abs diff) """
    assert df1[VARNAME_MU].size == df2[VARNAME_MU].size, "Oh no, looks like you have a different number of bins/supernova for your systematic and this is not good."
    diff = scale * ((df1[VARNAME_MU] - df1[VARNAME_MUREF]) - (df2[VARNAME_MU] - df2[VARNAME_MUREF])).to_numpy()
    diff[~np.isfinite(diff)] = 0
    cov = diff[:, None] @ diff[None, :]

    # Determine the gradient using simple linear regression
    reg = LinearRegression()
    weights = 1 / np.sqrt(0.003 ** 2 + df1[VARNAME_MUERR].to_numpy() ** 2 + df2[VARNAME_MUERR].to_numpy() ** 2)
    mask = np.isfinite(weights)
    reg.fit(df1.loc[mask, [VARNAME_z]], diff[mask], sample_weight=weights[mask])
    coef = reg.coef_[0]

    mean_abs_deviation = np.average(np.abs(diff), weights=weights)
    max_abs_deviation = np.max(np.abs(diff))
    return cov, (coef, mean_abs_deviation, max_abs_deviation)


def get_cov_from_covfile(data, covfile, scale):
    covindf = pd.read_csv(covfile,float_precision='high',low_memory=False)
    covindf['CID1'] = covindf['CID1'].astype(str)+"_"+covindf['IDSURVEY1'].astype(str)
    covindf['CID2'] = covindf['CID2'].astype(str)+"_"+covindf['IDSURVEY2'].astype(str)
    covindf = covindf.set_index(["CID1", "CID2"])
    covout = np.zeros((len(data),len(data)))
    for i,cid1,ra1,dec1 in zip(range(len(data)),data['CIDstr'],data['RA'],data['DEC']):
        nmissing = 0
        for j,cid2,ra2,dec2 in zip(range(len(data)),data['CIDstr'],data['RA'],data['DEC']):
            try:
                covout[i,j] = covindf.loc[(cid1,cid2), "MU_COV"]
            except KeyError:
                nmissing+=1
        if nmissing > 100:
            logging.info(f'{i} {cid1} RA {ra1} DEC {dec1} \t MISSING FROM COVFILE')
        
    return covout, (0, 0, 0)

def get_contributions(m0difs, fitopt_scales, muopt_labels, muopt_scales, extracovdict):
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

    # xxxx mark delete May 14 2021 (RK) xxxxxxxx
    #label, fitopt_filter, muopt_filter = \
    #    re.findall(r"\[(.*)\] \[(.*),(.*)\]", covopt)[0]
    # xxxxxxxx end mark xxxxxxxxx

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
                final_cov = cov.copy()
            else:
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
        # CosmoMC will add the diag terms, so lets do it here and make sure its all good
        effective_cov = final_cov + np.diag(base[VARNAME_MUERR] ** 2)

        # First just try and invert it to catch singular matrix errors
        precision = np.linalg.inv(effective_cov)

        # Then check that the matrix is well conditioned to deal with float precision
        epsilon = sys.float_info.epsilon
        cond = np.linalg.cond(effective_cov)
        assert cond < 1 / epsilon, "Cov matrix is ill-conditioned and cannot be inverted"
        logging.info(f"Covar condition for COVOPT {label} is {cond:.3f}")

        # Finally, re-invert the precision matrix and ensure its within 
        # tolerance of the original covariance
        cov2 = np.linalg.inv(precision)
        assert np.all(np.isclose(effective_cov, cov2)), "Double inversion does not give original covariance, matrix is unstable"

    except np.linalg.LinAlgError as ex:
        logging.exception(f"Unable to invert covariance matrix for COVOPT {label}")
        raise ex

    return label, final_cov


def write_dataset(path, data_file, cov_file, template_path):
    with open(template_path) as f_in:
        with open(path, "w") as f_out:
            f_out.write(f_in.read().format(data_file=data_file, cov_file=cov_file))


def write_data_jla(path, base, cosmomc_format=True):

    logging.info(f"Write Hubble diagram for cosmomc-jla method")
    if not cosmomc_format:
        base[[VARNAME_z, VARNAME_MU, VARNAME_MUERR]].to_csv(path, sep=" ", index=True, float_format="%.5f")
        return

    zs   = base[VARNAME_z].to_numpy()
    mu   = base[VARNAME_MU].to_numpy()
    mbs  = -19.36 + mu
    mbes = base[VARNAME_MUERR].to_numpy()

    # I am so sorry about this, but CosmoMC is very particular
    logging.info(f"Writing out data to {path}")
    with open(path, "w") as f:
        f.write("#name zcmb    zhel    dz mb        dmb     x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor\n")
        for i, (z, mb, mbe) in enumerate(zip(zs, mbs, mbes)):
            f.write(f"{i:5d} {z:6.5f} {z:6.5f} 0  {mb:8.5f} {mbe:8.5f} 0 0 0 0 0 0 0 0 0 0 0 0\n")



def write_data_bbc(path, base):
    # Dec 2020
    # Write format for new BBC method in cosmomc

    z_list      = base[VARNAME_z].to_numpy()
    mu_list     = base[VARNAME_MU].to_numpy()
    muerr_list  = base[VARNAME_MUERR].to_numpy()
    
    #sys.exit(f"\n xxx DEBUG DIE xxxx\n ")

    logging.info(f"Writing out data to {path}")
    with open(path, "w") as f:
        f.write("# name zcmb    zhel    dz  mu   muerr\n")
        for i, (z, mu, muerr) in enumerate(zip(z_list, mu_list, muerr_list)):
            f.write(f"{i:5d} {z:6.5f} {z:6.5f} 0  {mu:8.5f} {muerr:8.5f} \n")

def write_covariance(path, cov):

    cosmomc_method = config["COSMOMC_METHOD"]
    file_base      = os.path.basename(path)
    covdet         = np.linalg.det(cov)
    nrow           = cov.shape[0]

    logging.info(f"Write cov to {path}")

    # RK - write diagnostic to check if anything changes
    logging.info(f"    {file_base}: size={nrow}  |cov| = {covdet:.5e}")

    # - - - - -
    # Write out the matrix
    nwr = 0 ; rownum = 0; colnum=0

    with open(path, "w") as f:
        f.write(f"{nrow}\n")

        for c in cov.flatten():

            if cosmomc_method == COSMOMC_METHOD_JLA :
                f.write(f"{c:.14f}\n")
            else:
                # for bbc method write human-readable cov:
                # comment line for each new row, and off-diag elements
                # are indented with pad_space
                nwr += 1
                is_new_row = False
                pad_space  = "   "   # few spaces for off-diag
                if (nwr-1) % nrow == 0 : 
                    is_new_row = True ; rownum += 1 ; colnum = 0
                colnum += 1
                if rownum == colnum : 
                    pad_space = ""  # no pad space for diagonal

                if is_new_row:
                    f.write(f"# -------- Begin Row {rownum} of {nrow} " \
                            f"-----------\n")
                f.write(f"{pad_space}{c:11.8f}\n")


def write_cosmomc_output(config, covs, base):
    # Copy INI files. Create covariance matrices. Create .dataset. Modifying INI files to point to resources
    out = Path(config["OUTDIR"]) / "cosmomc"

    dataset_template = Path(config["COSMOMC_TEMPLATES"]) / config["DATASET_FILE"]    
    dataset_files = []
    cosmomc_method = config["COSMOMC_METHOD"]

    os.makedirs(out, exist_ok=True)

    # Create lcparam file
    if cosmomc_method == COSMOMC_METHOD_JLA :
        data_file      = out / f"data.txt"
        data_file_wCID = out / f"data_wCID.txt"
        prefix_covsys  = "sys"
        write_data_jla(data_file, base)
        write_data_jla(data_file_wCID, base, cosmomc_format=False)
    elif cosmomc_method == COSMOMC_METHOD_BBC :
        data_file      = out / f"hubble_diagram.txt"
        data_file_wCID = out / f"hubble_diagram_wCID.txt"
        prefix_covsys  = "covsyst"
        write_data_bbc(data_file, base)
    else:
        msg = f"Invalid COSMOMC method = {cosmomc_method}"
        raise ValueError(msg)

    # ?? Create supplementary file for people to merge covmat with fitres ??

    # Create covariance matrices and datasets
    for i, (label, cov) in enumerate(covs):
        dataset_file = out / f"dataset_{i}.txt"
        covsyst_file = out / f"{prefix_covsys}_{i}.txt" 

        write_covariance(covsyst_file, cov)
        write_dataset(dataset_file, data_file, covsyst_file, dataset_template)
        dataset_files.append(dataset_file)

    # Copy some INI files
    ini_files = [f for f in os.listdir(config["COSMOMC_TEMPLATES"]) if f.endswith(".ini") or f.endswith(".yml") or f.endswith(".md")]
    for ini in ini_files:
        op = Path(config["COSMOMC_TEMPLATES"]) / ini

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


def write_summary_output(config, covariances, base):
    out = Path(config["OUTDIR"])
    info = {}
    cov_info = {}
    for i, (label, cov) in enumerate(covariances):
        cov_info[i] = label
    info["COVOPTS"] = cov_info

    logging.info("Writing INFO.YML")
    with open(out / "INFO.YML", "w") as f:
        yaml.safe_dump(info, f)


def write_correlation(path, label, base_cov, diag, base):
    logging.debug(f"\tWriting out cov for COVOPT {label}")

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
    logging.info("Writing out summary.csv information")
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
    input_dir = Path(config["INPUT_DIR"])
    version   = config["VERSION"]
    data_dir  = input_dir / version
    sys_file  = Path(config["SYSFILE"])
    extra_covs = config.get("EXTRA_COVS",[])
    
    # Read in all the needed data
    submit_info = read_yaml(input_dir / "SUBMIT.INFO")
    sys_scale = read_yaml(sys_file)
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
        #muopt_labels = { args.muopt : tmp_label }
        muopt_labels = { args.muopt : "DEFAULT" }     
    
    muopt_scales = config["MUOPT_SCALES"]
    for extra in extra_covs:
        muopt_scales[extra.split()[0]] = extra.split()[1]

    # Load in all the hubble diagrams
    logging.info(f"Read Hubble diagrams for version = {version}")
    data = get_hubble_diagrams(data_dir, args, config)
    # Filter data to remove rows with infinite error
    data, base = remove_nans(data)
    # Now that we have the data, figure out how each much each FITOPT/MUOPT pair contributes to cov
    contributions, summary = get_contributions(data, fitopt_scales,
                                               muopt_labels, muopt_scales, extracovdict)
    # find contributions which match to construct covs for each COVOPT
    logging.info(f"Compute covariance for COVOPTS")
    covopts = ["[ALL] [,]"] + config.get("COVOPTS",[])  # Adds covopt to compute everything
    covariances = [get_cov_from_covopt(c, contributions, base, 
                                       config.get("CALIBRATORS")) for c in covopts]
    
    write_cosmomc_output(config, covariances, base)
    write_summary_output(config, covariances, base)
    write_debug_output(config, covariances, base, summary)

def prep_config(config,args):

    path_list = [ 'OUTDIR', 'SYSFILE', 'INPUT_DIR', 'COSMOMC_TEMPLATES' ]

    for path in path_list:
        config[path] = os.path.expandvars(config[path]) ;   

    # if cosmomc method is NOT specified, revert to legacy jla method
    # to preserve current/old inputs for pippin
    if KEYNAME_COSMOMC_METHOD not in config:
        config[KEYNAME_COSMOMC_METHOD] = COSMOMC_METHOD_JLA

    # WARNING: later add option to read from input file
    #sys.exit(f" xxx nbin(x1,c) = {args.nbin_x1} {args.nbin_c} ")
    config['nbin_x1'] = args.nbin_x1
    config['nbin_c']  = args.nbin_c

    # check override args (RK, Feb 15 2021)

    if args.method is not None :
        config["COSMOMC_METHOD"] = args.method
        logging.info(f"OPTION: override COSMOMC_METHOD with {args.method}")
        
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
