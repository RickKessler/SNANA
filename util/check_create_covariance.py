#!/usr/bin/env python
#
# Created Feb 2025 by A.Mitra and R.Kessler
# Compute one element of covariance matrix using input pair of {CID,IDSURVEY}.
# Then read element from output of create_covariance.py and compare.
#
# Sep 17 2025: refactor to work with ...
#   + reading text or npz covsys file.
#   + binned HD that has no IDSURVEY column
#   + cov_num > 0 (see get_bbc_fitres_list)
#
# ================================================================
import os, argparse, logging, shutil, time, datetime, subprocess
import re, yaml, sys, gzip, math, gc, glob
import numpy  as np
import pandas as pd

# ===================================

HD_METHOD_BINNED = 'BINNED'
HD_METHOD_REBIN  = 'REBINxxx'  # beware that check doesn't work for rebin
HD_METHOD_UNBIN  = 'UNBIN'


SUFFIX_FITOPT_DICT = {
    HD_METHOD_BINNED:  'M0DIF' ,
    HD_METHOD_UNBIN:   'FITRES'
}

VARNAME_IDROW_DICT = {
    HD_METHOD_BINNED:  'ROW' ,
    HD_METHOD_UNBIN:   'CID'
}

VARNAME_MU_DICT = {
    HD_METHOD_BINNED:  'MUDIF' ,
    HD_METHOD_UNBIN:   'MU'
}

VARNAME_IDROW    = None   # loaded in main after reading INFO.YML
VARNAME_IDSURVEY = 'IDSURVEY'
VARNAME_MU       = None   # loaded in main after reading INFO.YML
VARNAME_MUMODEL  = 'MUMODEL'  # for UNBIN only

OUTFILE_PREFIX      = "out_check_covnum"

# ===================================================

def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP menu for config options"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of the existing covariance directory"
    parser.add_argument("-d", "--cov_dir", help=msg,
                        nargs="?", type=str, default=None)

    msg = "covsys number (DEFAULT = 0 for all systematics) "
    parser.add_argument("-n", "--cov_num", help=msg,
                        nargs='?', type=int, default=0)

    msg = "CID,IDSURVEY pair (rg. 100,12 120,50 -> CID = 100 for LSST, CID = 120 for LOW-z)"
    parser.add_argument("-c", "--cidpair", help=msg,
                        nargs='+', type=str, default=None)

 
    args = parser.parse_args()

    return args  

def read_yaml(path):
    path_expand = os.path.expandvars(path)
    logging.debug(f"Reading YAML from {path_expand}")
    with open(path_expand) as f:
        return yaml.safe_load(f.read())

def read_covinfo(args):

    # read INFO.YAML from cov directory and store the following:
    # BBC_DIR, SIZE_HD

    
    info_yaml    = args.cov_dir + "/INFO.YML"
    contents     = read_yaml(info_yaml)
    contents['INFO_YAML_FILE'] = info_yaml

    hd_bin_method = contents['HD_BIN_METHOD'] 
    contents['IS_BINNED'] = (hd_bin_method == HD_METHOD_BINNED)
    contents['IS_UNBIN']  = (hd_bin_method == HD_METHOD_UNBIN)

    return contents
    # end read_covinfo

def get_sys_scale(args, covinfo):

    bbc_info_file  = covinfo['BBC_INFO_FILE']
    sys_scale_file = covinfo['SYS_SCALE_FILE']

    contents_sys_scale = read_yaml(sys_scale_file)
    contents_bbcinfo   = read_yaml(bbc_info_file)["FITOPT_OUT_LIST"]
 
    sys_scale_dict = {}  # at create_cov stage
    sys_label_dict = {}  # from bbc stage

    for line in contents_bbcinfo:
        fitopt = line[0] 
        label  = line[2]
        
        scale = contents_sys_scale.setdefault(label, 1.0)
        #print(' xxx Fitopt  label   scale = ', fitopt, label, scale)
        sys_scale_dict[fitopt] = scale
        sys_label_dict[fitopt] = label
    
    covinfo['sys_scale_dict'] = sys_scale_dict
    covinfo['sys_label_dict'] = sys_label_dict

    return 
    # end get_sys_scale

def parse_cidpair(args):
    cid_list      = []
    idsurvey_list = []
    
    for cidpair in args.cidpair:
        cid_split = cidpair.split(",")
        cid = str(cid_split[0])
        if len(cid_split) > 1:
            idsurvey = int(cidpair.split(",")[1])
        else:
            idsurvey = None

        cid_list.append(cid)
        idsurvey_list.append(idsurvey)
    args.cid_list = cid_list
    args.idsurvey_list = idsurvey_list
    
    return args

def get_uidrow(idrow,idsurvey):
    
    # return unique idrow combining idrow and idsurvey
    if idsurvey is None:
        uidrow = str(idrow)
    else:
        uidrow = str(idrow) + '+' + str(idsurvey)

    return uidrow

def compute_scaled_covpair(args, covinfo):
    """
    For each FITRES file (each FITOPT), this function:
      - First, loads the FITRES file with FITOPT identifier "FITOPT000" to obtain reference MU values 
        (one per CID/IDSURVEY pair). These reference values replace the MUMODEL values.
      - Then, for each FITRES file, it:
           * Finds the first matching row for each of two CID/IDSURVEY pairs.
           * Computes delta = current MU - (reference MU from FITOPT000).
           * Scales each delta by the systematic scale factor (sys_scale).
           * Computes the final delta as the product: (delta1 * sys_scale) * (delta2 * sys_scale).
           * Saves all the info in a dictionary with individual columns.
    
    Returns a DataFrame with columns:
      [ 'FITOPT', 'CID1', 'MU1', 'RefMU1', 'Delta1', 
        'CID2', 'MU2', 'RefMU2', 'Delta2', 'SysScale', 'FinalDelta' ]
    """

    sys_scale_dict = covinfo['sys_scale_dict']
    is_unbin       = covinfo['IS_UNBIN']
    is_binned      = covinfo['IS_BINNED']
    has_idsurvey = args.idsurvey_list[0] is not None

    # get list of FITRES (or M0DIF) files from BBC directory
    fitres_list = get_bbc_fitres_list(args,covinfo)

    # --- NEW CODE: Find the FITRES file with FITOPT identifier "FITOPT000" and load its reference MU values ---
    ref_file = None
    for f in fitres_list:
        fitopt_match = re.search(r'(FITOPT(\d+))', f)
        if fitopt_match and fitopt_match.group(2) == "000":
            ref_file = f
            break
    if not ref_file:
        sys.exit("\n ERROR: No FITOPT000 file found! Cannot obtain reference MU values.")

    # Load the reference file and extract MU for each CID/IDSURVEY pair (using the first matching row)
    df_ref = pd.read_csv(ref_file, comment='#', sep=r'\s+')
    df_ref[VARNAME_IDROW]  = df_ref[VARNAME_IDROW].astype(str)

    if has_idsurvey:    
        df_ref[VARNAME_IDSURVEY] = df_ref[VARNAME_IDSURVEY].astype(int)

    ref_mu      = {}  # dictionary mapping CID to its reference MU from FITOPT000
    ref_mumodel = {}

    for idrow, idsurvey in zip(args.cid_list, args.idsurvey_list):
        
        matching_rows = get_matching_rows(df_ref, idrow, idsurvey)

        if matching_rows.empty :
            msgerr = f"ERROR: no matching row for {VARNAME_IDROW}={idrow} IDSURVEY={idsurvey} ; " \
                     f"\n\t see {ref_file}"
            sys.exit(f"\n{msgerr}")
        else:
            uidrow              = get_uidrow(idrow,idsurvey)
            ref_mu[uidrow]      = matching_rows.iloc[0][VARNAME_MU]
            if is_unbin:
                ref_mumodel[uidrow] = matching_rows.iloc[0][VARNAME_MUMODEL]
            else:
                ref_mumodel[uidrow] = 0.0

    #sys.exit(f"\n xxx ref_mu = {ref_mu}")

    # ----------------------------------------------------------------------------------------------

    cov_info_list = []  # list to store each row's information as a dictionary

    cov_check = 0.0
    # Loop over each FITRES file
    for fitres in fitres_list:

        fitres_base = os.path.basename(fitres)

        # Extract the full FITOPT string (e.g. "FITOPT000") and the numeric part
        fitopt_match = re.search(r'(FITOPT(\d+))', fitres_base)

        fitopt_str = fitopt_match.group(1)         # e.g. FITOPT000
        fitopt_num = int(fitopt_match.group(2))    # For indexing sys_scale_list

        # xxx mark sys_scale = sys_scale_list[fitopt_num] if fitopt_num < len(sys_scale_list) else 1.0

        sys_scale = sys_scale_dict[fitopt_str]

        #print(f" xxx {fitres_base}  str={fitopt_str}  num={fitopt_num:2d}  scale={sys_scale}")  # .xyz

        # Load the current FITRES file into a DataFrame
        df = pd.read_csv(fitres, comment='#', sep=r'\s+')
        df[VARNAME_IDROW]    = df[VARNAME_IDROW].astype(str)
        if has_idsurvey:
            df[VARNAME_IDSURVEY] = df[VARNAME_IDSURVEY].astype(int)

        # Assume there are exactly two CID/IDSURVEY pairs in args.
        # We'll use the first matching row for each.
        pair_results = []  # list to store result for each CID        

        for idrow, idsurvey in zip(args.cid_list, args.idsurvey_list):

            matching_rows = get_matching_rows(df, idrow, idsurvey)

            # Use the first matching row
            row         = matching_rows.iloc[0]
            mu_val      = row[VARNAME_MU]
            if is_unbin:
                mumodel_val = row[VARNAME_MUMODEL]
            else:
                mumodel_val = 0.0

            # Use the reference MU from FITOPT000 as MU0.
            # Note that subtracting mumodel is needed for zshift sytematics;
            uidrow             = get_uidrow(idrow,idsurvey)
            reference_mu       = ref_mu[uidrow]
            reference_mumodel  = ref_mumodel[uidrow]
            delta_val    = (mu_val - mumodel_val) - \
                           (reference_mu - reference_mumodel)  # compute delta as current MU minus reference MU

            scaled_delta = delta_val * sys_scale
            pair_results.append((idrow, mu_val, reference_mu, delta_val, scaled_delta))

        # - - - - - - -
        # Unpack the two results (assumed order corresponds to the order in args)
        id1, mu1, ref_mu1, delta1, scaled_delta1 = pair_results[0]
        id2, mu2, ref_mu2, delta2, scaled_delta2 = pair_results[1]

        # Compute final delta as the product of the two scaled deltas
        final_delta    = scaled_delta1 * scaled_delta2
        cov_check     += final_delta

        # Build a dictionary with all the desired fields.
        row_dict = {
            "FITOPT"  : fitopt_str,
            "ID1"     : id1,
            "MU1"     : mu1,
            "RefMU1"  : ref_mu1,
            "Delta1"  : scaled_delta1,
            "ID2"     : id2,
            "MU2"     : mu2,
            "RefMU2"  : ref_mu2,
            "Delta2"  : scaled_delta2,
            "SysScale"   : sys_scale,
            "Covariance" : final_delta
        }
        cov_info_list.append(row_dict)

    # Create a DataFrame from the collected results
    df_cov_info = pd.DataFrame(cov_info_list)
    
    return df_cov_info, cov_check
    # end compute_scaled_covpair


def get_bbc_fitres_list(args,covinfo):
    # read list of FITOPTs from COV's INFO.YML file (covifo dict),
    # and select appropriate fitres files from BBC_DIR
    # If cov_num=0, then take them all.

    cov_num        = args.cov_num
    bbc_dir        = covinfo['BBC_DIR']
    hd_bin_method  = covinfo['HD_BIN_METHOD']
    sys_label_dict = covinfo['sys_label_dict']  # labels from BBC dir (passed from LCFIT)
    suffix         = SUFFIX_FITOPT_DICT[hd_bin_method]

    wildcard       = os.path.join(bbc_dir, f'FITOPT*.{suffix}.gz')
    fitres_list_all = sorted(glob.glob(wildcard))  # Load all FITRES files
    
    if cov_num == 0:
        fitres_list = fitres_list_all
    else:
        # go fishing for subset
        fitres_list = []
        covopt = covinfo['COVOPTS'][cov_num]
        cov_label = covopt.split()[0]  # label in create_cov input
        for ff in fitres_list_all:
            ff_base   = os.path.basename(ff)
            fitopt    = ff_base.split('_MUOPT')[0]  # e.g. FITOPT004
            bbc_label = sys_label_dict[fitopt]
            keep_ref = 'FITOPT000' in ff_base
            keep_sys = cov_label in bbc_label
            if keep_ref or keep_sys :
                fitres_list.append(ff)
                if keep_ref:
                    print(f"\t keep {fitopt} as reference.")
                else:
                    print(f"\t keep {fitopt}  {bbc_label} contains {cov_label} for cov_num={cov_num}")

    #sys.exit(f"\n xxx stop debug")

    return fitres_list
    # end get_bbc_fitres_list

def get_matching_rows(df, idrow, idsurvey):
    # Created Sep 16 2025 
    if idsurvey is None:
        matching_rows = df[ (df[VARNAME_IDROW] == idrow) ]
    else:
        matching_rows = df[ (df[VARNAME_IDROW]    == idrow) & \
                            (df[VARNAME_IDSURVEY] == idsurvey) ]

    return matching_rows
    
def get_cov_from_create_covariance(args, covinfo):

    cov_dir = args.cov_dir
    hd_file = os.path.join(cov_dir, "hubble_diagram.txt")

    is_unbin = covinfo['IS_UNBIN']

    # --- Step 1: Read the Hubble diagram file ---
    # The hubble diagram file is assumed to be space-delimited with commented header lines.

    print(f"\n Find hubble diagram rows in \n {hd_file}")
    sys.stdout.flush()

    df_hd = pd.read_csv(hd_file, comment='#', sep=r'\s+')    
    df_hd[VARNAME_IDROW]    = df_hd[VARNAME_IDROW].astype(str).str.strip()
    if is_unbin:
        df_hd[VARNAME_IDSURVEY] = df_hd[VARNAME_IDSURVEY].astype(int)

    # Find the row indices for each CID/IDSURVEY pair in your arguments.
    row_indices = []
    for idrow, idsurvey in zip(args.cid_list, args.idsurvey_list):

        # Note: Ensure the data types match (CID as string, IDSURVEY as integer)

        if idsurvey is None:
            idx = df_hd[(df_hd[VARNAME_IDROW] == (idrow)) ].index
        else:
            idx = df_hd[(df_hd[VARNAME_IDROW] == (idrow)) & (df_hd[VARNAME_IDSURVEY] == (idsurvey))].index

        if idx.empty:
            msgerr = f"ERROR: row not found for {VARNAME_IDROW}={idrow}, IDSURVEY={idsurvey} \n"\
                     f"\t in hubble diagram file {hd_file}"
            sys.exit(f"\n{msgerr}")
        else:
            # We take the first matching row
            row_indices.append(idx[0])
    
    print(" Found row indices in hubble diagram:")
    for idrow, idsurvey, row in zip(args.cid_list, args.idsurvey_list, row_indices):
        print(f"\t HD row = {row:5d} for {VARNAME_IDROW}={idrow:<10}  IDSURVEY={idsurvey}")
        sys.stdout.flush()
    
    # Check that you have at least two indices to compare.
    nrow = len(row_indices)
    if nrow < 2:
        msgerr = f"ERROR: found only {nrow} rows (expected 2) in " \
                 f"\n\t HD file {hd_file}"
        sys.exit(f"\n {msgerr}")
        return None

    # --- Step 2: Read the flattened covariance matrix file ---
    # The file is assumed to be gzipped, with the first line as the dimension (e.g. 1802)
    # and subsequent lines as the matrix elements (row-major order).

    cov_matrix = read_covsys_file(args)
    
    print("Covariance matrix shape:", cov_matrix.shape)
    sys.stdout.flush()

    # --- Step 3: Extract the covariance entry for the two selected rows ---
    # For example, if you have two indices (row_i, row_j):
    row_i = row_indices[0]
    row_j = row_indices[1]
    covariance_value = cov_matrix[row_i, row_j]
    print(f"         Covariance between row {row_i} and row {row_j}: {covariance_value}")
    print(f"Reversed Covariance between row {row_j} and row {row_i}: {covariance_value}")
    print(f"")
    sys.stdout.flush()

    return covariance_value
    # end get_cov_from_create_covariance 

def read_covsys_file(args):

    cov_dir = args.cov_dir
    cov_num = args.cov_num

    # check if covsys file is text of npz, and if it even exists
    covsys_file_txt  = os.path.join(cov_dir, f"covsys_{cov_num:03d}.txt.gz")
    covsys_file_npz  = os.path.join(cov_dir, f"covsys_{cov_num:03d}.npz")
    covsys_file_list = [ covsys_file_txt, covsys_file_npz ]
    covsys_file      = None
    for tmp_file in covsys_file_list:
        if os.path.exists(tmp_file): covsys_file = tmp_file

    if not covsys_file:
        base_txt = os.path.basename(covsys_file_txt)
        base_npz = os.path.basename(covsys_file_npz)
        sys.exit(f"\n ERROR: cannot find {base_txt} or {base_npz} in \n" \
                 f"\t {cov_dir}\n" \
                 f"\t Rerun create_covariance.py with option to write covsys in addition to covtot_inv")

    print(f"\n Read create_covariance element from \n\t {covsys_file}")
    sys.stdout.flush()

    if covsys_file == covsys_file_txt:
        covsys = read_covsys_txt(covsys_file)
    elif covsys_file == covsys_file_npz:
        covsys = read_covsys_npz(covsys_file)
    else:
        pass

    return covsys
    # end read_covsys_file

def read_covsys_npz(covsys_file):

    loaded_data     = np.load(covsys_file)
    dim             = loaded_data['nsn'][0]
    cov_upper       = loaded_data['cov']

    # restore full symmetric matrix from upper-triangle matrix

    cov_matrix    = np.empty((dim,dim)) 
    upper_indices = np.triu_indices(dim, k=0) # k=0 includes the diagonal, k=1 excludes it
    cov_matrix[upper_indices]   = cov_upper
    cov_matrix.T[upper_indices] = cov_upper

    return cov_matrix

def read_covsys_txt(covsys_file):

    # Created Sep 16 2025 [moved from covariance_value]
    # Read covsys from text file

    with gzip.open(covsys_file, 'rt') as f:
        # Read the first non-empty line which gives the dimension.
        first_line = f.readline().strip()
        if not first_line:
            sys.exit("\n ERROR: Covariance file is empty or does not contain dimension info.")
        try:
            dim = int(first_line)
            print(f"\t cov dimention to read: {dim} x {dim}")
            sys.stdout.flush()
        except ValueError:
            sys.exit("\n ERROR: Unable to parse cov_dim = '{first_line}' from first row of cov file.")

        # Now read the rest of the file and parse the covariance elements as floats.
        cov_elements = []
        for line in f:
            line = line.strip()
            if line:  # skip blank lines
                # The file may have one element per line or several;
                # we split by whitespace to be safe.
                parts = line.split()
                for part in parts:
                    cov_elements.append(float(part))
    
    # Check if we have the correct number of elements.
    if len(cov_elements) != dim * dim:
        sys.exit("\n ERROR: Expected", dim * dim, "elements, but got", len(cov_elements))

    cov_matrix = np.array(cov_elements).reshape(dim, dim)

    return cov_matrix

    # end read_covsys_txt

def print_cov_results(args, results_dict):

    cid_list   = args.cid_list
    idsrv_list = args.idsurvey_list

    cid0   = cid_list[0]
    cid1   = cid_list[1]
    idsrv0 = idsrv_list[0]
    idsrv1 = idsrv_list[1]
    if idsrv0 is None:
        id0  = cid0
        id1  = cid1
    else:
        id0  = f"{cid0}+{idsrv0}"
        id1  = f"{cid1}+{idsrv0}"

    dff_cov_info                  = results_dict['dff_cov_info']
    cov_check                     = results_dict['cov_check']
    cov_from_create_covariance    = results_dict['cov_from_create_covariance']

    # - - - - - - - - - - - 
    # write info to stdout and to file
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 0)  # Adjust to your console width if needed
    pd.set_option('display.expand_frame_repr', True)

    print(df_cov_info,'\n');         sys.stdout.flush()

    tmp_string = f"{id0}_{id1}"
    outfile    = f"{OUTFILE_PREFIX}{args.cov_num}_{tmp_string}.csv"

    df_cov_info.to_csv(outfile,  sep=' ', index=False)
    print(f"# Write cov_check terms to {outfile}") ; 
    sys.stdout.flush()

    # - - - - - - - - 
    comment_list  = [
        "covariance = (delta1 * sys scale) * (delta2 * sys scale)",
        "delta_i = MU_i - MU0 ",
        " ",
    ]


    cov_string = f"cov[ {id0} , {id1} ]"
    comment_list.append(f"Computed for crosscheck     : {cov_string} = {cov_check}")
    comment_list.append(f"Read from create_covariance : {cov_string} = {cov_from_create_covariance}")

    ratio = cov_check/cov_from_create_covariance
    comment_list.append(f"Ratio(computed/read) : {ratio:10.5f}  for {cov_string}")
    
    for comment in comment_list:
        print(f"# {comment}")
        sys.stdout.flush()

    return
    # end print_cov_results

# ============================================================    
# ============================================================    

if __name__ == "__main__":

    args = get_args()
    args = parse_cidpair(args)
    covinfo = read_covinfo(args)

    # set table column names based on HC binning method
    hd_bin_method = covinfo['HD_BIN_METHOD'] 
    VARNAME_IDROW = VARNAME_IDROW_DICT[hd_bin_method]  # CID or ROW
    VARNAME_MU    = VARNAME_MU_DICT[hd_bin_method]     # MU or MUDIF

    # read scale & label per systematic (store in covinfo
    get_sys_scale(args, covinfo)
    
    # compute cov(cid0,cid1); this is the cross check value
    df_cov_info, cov_check = compute_scaled_covpair(args, covinfo)


    # fetch cov(cid0,cid1) from output of create_covariance
    cov_from_create_covariance = get_cov_from_create_covariance(args, covinfo)

    #sys.exit(f"\n xxx covinfo = \n{covinfo}")

    # print results to stdout
    results_dict = {
        'dff_cov_info'  :  df_cov_info,
        'cov_check'     :  cov_check,
        'cov_from_create_covariance' : cov_from_create_covariance
    }
    print_cov_results(args, results_dict)

    # === END: ===
