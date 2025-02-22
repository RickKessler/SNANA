#!/usr/bin/env python
#
# Created Feb 2025 by A.Mitra and R.Kessler
# Compute one element of covariance matrix using input pair of {CID,IDSURVEY}.
# Then read element from output of create_covariance.py and compare.
#

# ================================================================
import os, argparse, logging, shutil, time, datetime, subprocess
import re, yaml, sys, gzip, math, gc, glob
import numpy  as np
import pandas as pd

# ===================================
VARNAME_CID      = 'CID'
VARNAME_IDSURVEY = 'IDSURVEY'
VARNAME_MU       = 'MU'
VARNAME_MUMODEL  = 'MUMODEL'

#xxx mark OUTFILE_COV_INFO    = "cov_info.csv"
OUTFILE_PREFIX      = "out_check_cov"

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

    info_yaml = args.cov_dir+"/INFO.YML"
    contents = read_yaml(info_yaml)
    args.bbc_dir = contents['BBC_DIR']
    args.sizehd  = contents['SIZE_HD']
    args.bbc_info_file  = contents['BBC_INFO_FILE']
    args.sys_scale_file = contents['SYS_SCALE_FILE']
    
    
    return args

def get_sys_scale(args):

    bbc_info_file  = args.bbc_info_file
    sys_scale_file = args.sys_scale_file

    contents_sys_scale = read_yaml(sys_scale_file)
    contents_bbcinfo   = read_yaml(bbc_info_file)["FITOPT_OUT_LIST"]
 
    syst_scale_list = []

    for line in contents_bbcinfo:
        fitopt = line[0] 
        label  = line[2]

        scale = contents_sys_scale.setdefault(label, 1.0)
        #print('XXX Fitopt label, scale = ', fitopt, label, scale)
        syst_scale_list.append(scale)
        
    return syst_scale_list

def parse_cidpair(args):
    cid_list      = []
    idsurvey_list = []
    
    for cidpair in args.cidpair:
        cid = str(cidpair.split(",")[0])
        idsurvey = int(cidpair.split(",")[1])
        cid_list.append(cid)
        idsurvey_list.append(idsurvey)
    args.cid_list = cid_list
    args.idsurvey_list = idsurvey_list
    
    return args

def get_idrow(cid,idsurvey):
    idrow = str(cid) + '+' + str(idsurvey)
    return idrow

def compute_scaled_covpair(args, sys_scale_list):  
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

    bbc_dir = args.bbc_dir
    wildcard = os.path.join(bbc_dir, 'FITOPT*FITRES.gz')
    fitres_list = sorted(glob.glob(wildcard))  # Load all FITRES files

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
    df_ref = pd.read_csv(ref_file, comment='#', sep='\s+')
    df_ref[VARNAME_CID]      = df_ref[VARNAME_CID].astype(str)
    df_ref[VARNAME_IDSURVEY] = df_ref[VARNAME_IDSURVEY].astype(int)
    ref_mu      = {}  # dictionary mapping CID to its reference MU from FITOPT000
    ref_mumodel = {}

    for cid, idsurvey in zip(args.cid_list, args.idsurvey_list):
        matching_rows = df_ref[(df_ref[VARNAME_CID] == cid) & \
                               (df_ref[VARNAME_IDSURVEY] == idsurvey)]
        if matching_rows.empty :
            msgerr = f"ERROR: no matching row for CID={cid} IDSURVEY={idsurvey} ; " \
                     f"\n\t see {ref_file}"
            sys.exit(f"\n{msgerr}")
        else:
            idrow              = get_idrow(cid,idsurvey)
            ref_mu[idrow]      = matching_rows.iloc[0][VARNAME_MU]
            ref_mumodel[idrow] = matching_rows.iloc[0][VARNAME_MUMODEL]

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

        # Get the systematic scale factor (default to 1.0 if out-of-range)
        sys_scale = sys_scale_list[fitopt_num] if fitopt_num < len(sys_scale_list) else 1.0

        #print(f" xxx {fitres_base}  str={fitopt_str}  num={fitopt_num:2d}  scale={sys_scale}")  # .xyz

        # Load the current FITRES file into a DataFrame
        df = pd.read_csv(fitres, comment='#', sep='\s+')
        df[VARNAME_CID]      = df[VARNAME_CID].astype(str)
        df[VARNAME_IDSURVEY] = df[VARNAME_IDSURVEY].astype(int)

        # Assume there are exactly two CID/IDSURVEY pairs in args.
        # We'll use the first matching row for each.
        pair_results = []  # list to store result for each CID        

        for cid, idsurvey in zip(args.cid_list, args.idsurvey_list):
            matching_rows = df[ (df[VARNAME_CID]      == cid) & \
                                (df[VARNAME_IDSURVEY] == idsurvey) ]

            # Use the first matching row
            row         = matching_rows.iloc[0]
            mu_val      = row[VARNAME_MU]
            mumodel_val = row[VARNAME_MUMODEL]

            # Use the reference MU from FITOPT000 as MU0.
            # Note that subtracting mumodel is needed for zshift sytematics;
            idrow             = get_idrow(cid,idsurvey)
            reference_mu      = ref_mu[idrow]
            reference_mumodel = ref_mumodel[idrow]
            delta_val    = (mu_val - mumodel_val) - \
                           (reference_mu - reference_mumodel)  # compute delta as current MU minus reference MU

            scaled_delta = delta_val * sys_scale
            pair_results.append((cid, mu_val, reference_mu, delta_val, scaled_delta))

        # - - - - - - -
        # Unpack the two results (assumed order corresponds to the order in args)
        cid1, mu1, ref_mu1, delta1, scaled_delta1 = pair_results[0]
        cid2, mu2, ref_mu2, delta2, scaled_delta2 = pair_results[1]

        # Compute final delta as the product of the two scaled deltas
        final_delta    = scaled_delta1 * scaled_delta2
        cov_check     += final_delta

        # Build a dictionary with all the desired fields.
        row_dict = {
            "FITOPT"  : fitopt_str,
            "CID1"    : cid1,
            "MU1"     : mu1,
            "RefMU1"  : ref_mu1,
            "Delta1"  : scaled_delta1,
            "CID2"    : cid2,
            "MU2"     : mu2,
            "RefMU2"  : ref_mu2,
            "Delta2"  : scaled_delta2,
            "SysScale"   : sys_scale,
            "Covariance" : final_delta
        }
        cov_info_list.append(row_dict)

    # Create a DataFrame from the collected results
    df_cov_info = pd.DataFrame(cov_info_list)
    #print('XXX df_result covar', df_results['Covariance'])
    
    
    return df_cov_info, cov_check
    # end compute_scaled_covpair


def get_cov_from_create_covariance(args):

    cov_dir = args.cov_dir
    cov_num = args.cov_num
    hd_file = os.path.join(cov_dir, "hubble_diagram.txt")
    covsys_file = os.path.join(cov_dir, f"covsys_{cov_num:03d}.txt.gz")

    # --- Step 1: Read the Hubble diagram file ---
    # The hubble diagram file is assumed to be space-delimited with commented header lines.

    print(f"\n Find hubble diagram rows in \n {hd_file}")
    sys.stdout.flush()

    df_hd = pd.read_csv(hd_file, comment='#', delim_whitespace=True)
    df_hd[VARNAME_CID]      = df_hd[VARNAME_CID].astype(str).str.strip()
    df_hd[VARNAME_IDSURVEY] = df_hd[VARNAME_IDSURVEY].astype(int)

    # Find the row indices for each CID/IDSURVEY pair in your arguments.
    row_indices = []
    for cid, idsurvey in zip(args.cid_list, args.idsurvey_list):
        # Note: Ensure the data types match (CID as string, IDSURVEY as integer)
        idx = df_hd[(df_hd['CID'] == (cid)) & (df_hd['IDSURVEY'] == (idsurvey))].index
        if idx.empty:
            msgerr = f"ERROR: No row found for CID={cid}, IDSURVEY={idsurvey} "\
                     f"in hubble diagram file {hd_file}"
            sys.exit(f"\n{msgerr}")
        else:
            # We take the first matching row
            row_indices.append(idx[0])
    
    print(" Found row indices in hubble diagram:")
    for cid, idsurvey, row in zip(args.cid_list, args.idsurvey_list, row_indices):
        print(f"\t HD row = {row:5d} for CID={cid:<10}  IDSURVEY={idsurvey}")
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

    print(f"\n Read create_covariance element from \n\t {covsys_file}")
    sys.stdout.flush()

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
    
    # Convert the list of elements to a NumPy array and reshape it.
    cov_matrix = np.array(cov_elements).reshape(dim, dim)
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


def print_cov_results(args, results_dict):

    cid_list   = args.cid_list
    idsrv_list = args.idsurvey_list
    id0        = f"{cid_list[0]}+{idsrv_list[0]}"
    id1        = f"{cid_list[1]}+{idsrv_list[1]}"

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
    outfile    = f"{OUTFILE_PREFIX}_{tmp_string}.csv"

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
    args = read_covinfo(args)

    # read scale per systematic
    sys_scale_list = get_sys_scale(args)
    
    # compute cov(cid0,cid1)
    df_cov_info, cov_check = compute_scaled_covpair(args, sys_scale_list)
    
    # fetch cov(cid0,cid1) from output of create_covariance
    cov_from_create_covariance = get_cov_from_create_covariance(args)

    # print results to stdout
    results_dict = {
        'dff_cov_info'  :  df_cov_info,
        'cov_check'     :  cov_check,
        'cov_from_create_covariance' : cov_from_create_covariance
    }
    print_cov_results(args, results_dict)

    # === END: ===
