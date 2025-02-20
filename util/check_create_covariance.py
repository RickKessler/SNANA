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


VARNAME_CID      = 'CID'
VARNAME_IDSURVEY = 'IDSURVEY'
VARNAME_MU       = 'MU'
VARNAME_MUMODEL  = 'MUMODEL'

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
    import pandas as pd  # ensure pandas is imported
    import re, os, glob

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
        print("No FITOPT000 file found! Cannot obtain reference MU values.")
        return pd.DataFrame()  # or raise an error

    # Load the reference file and extract MU for each CID/IDSURVEY pair (using the first matching row)
    df_ref = pd.read_csv(ref_file, comment='#', sep='\s+')
    df_ref[VARNAME_CID]      = df_ref[VARNAME_CID].astype(str)
    df_ref[VARNAME_IDSURVEY] = df_ref[VARNAME_IDSURVEY].astype(int)
    ref_mu      = {}  # dictionary mapping CID to its reference MU from FITOPT000
    ref_mumodel = {}

    for cid, idsurvey in zip(args.cid_list, args.idsurvey_list):
        matching_rows = df_ref[(df_ref[VARNAME_CID] == cid) & \
                               (df_ref[VARNAME_IDSURVEY] == idsurvey)]
        if matching_rows.empty:
            print(f"Reference: No matching row found for CID {cid} IDSURVEY {idsurvey} in FITOPT000")
        else:
            ref_mu[cid]      = matching_rows.iloc[0][VARNAME_MU]
            ref_mumodel[cid] = matching_rows.iloc[0][VARNAME_MUMODEL]

    # ----------------------------------------------------------------------------------------------

    results_list = []  # list to store each row's information as a dictionary

    sum_cov_terms = 0.
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
            #print('XXX matching rows ', matching_rows)
            #if matching_rows.empty:
            #    print(f"{fitopt_str}: No matching row found for CID {cid} IDSURVEY {idsurvey}")
            #    break  # Skip this FITRES file if one of the pairs is missing

            # Use the first matching row
            row         = matching_rows.iloc[0]
            mu_val      = row[VARNAME_MU]
            mumodel_val = row[VARNAME_MUMODEL]

            # Use the reference MU from FITOPT000 as MU0.
            # Note that subtracting mumodel is needed for zshift sytematics;
            reference_mu      = ref_mu.get(cid)
            reference_mumodel = ref_mumodel.get(cid)
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
        sum_cov_terms += final_delta

        # Build a dictionary with all the desired fields.
        row_dict = {
            "FITOPT": fitopt_str,
            "CID1"  : cid1,
            "MU1"   : mu1,
            "RefMU1": ref_mu1,
            "Delta1": scaled_delta1,
            "CID2"  : cid2,
            "MU2"   : mu2,
            "RefMU2": ref_mu2,
            "Delta2": scaled_delta2,
            "SysScale": sys_scale,
            "Covariance": final_delta
        }
        results_list.append(row_dict)

    # Create a DataFrame from the collected results
    df_results = pd.DataFrame(results_list)
    #print('XXX df_result covar', df_results['Covariance'])
    
    return df_results, sum_cov_terms


def compute_covpair(args, sys_scale_list):

    # xxxxxxx OBDSOLETE xxxxxxxx
    #
    # Compute thecovariance matrix elements 
    # Not the matrix, just a single element
    # Print this element to screen and return
    # element values

    cov_dir = args.cov_dir
    cov_num = args.cov_num

    bbc_dir = args.bbc_dir
    wildcard = f'{bbc_dir}/FITOPT*FITRES.gz'

    fitres_list = sorted(glob.glob(wildcard))  # Load all FITOPT files
    mu_dict = {}  # Dictionary to store MU values for each CID pair

    print("\n=== Processing FITRES Files ===")
    
    # xxxxxxx OBDSOLETE xxxxxxxx

    for fitres in fitres_list:
        print('\n ------')
        print(f"\nReading FITRES file: {fitres}")  
        

        # Load the FITRES data
        df = pd.read_csv(fitres, comment='#', sep='\s+')

        df['CID'] = df['CID'].astype(str)
        df['IDSURVEY'] = df['IDSURVEY'].astype(int)

        for cid, idsurvey in zip(args.cid_list, args.idsurvey_list):
            
            print('CID , IDSURVEY ', (cid,idsurvey))
            matching_rows = df[(df['CID'] == cid) & (df['IDSURVEY'] == idsurvey)]

            if not matching_rows.empty:
                print(f"\n CID {cid} found in {fitres}")
                print(matching_rows[['CID', 'MU', 'MUMODEL']])  # Debugging actual MU & MUMODEL values

                mu_values = matching_rows['MU'].tolist()
                muref_values = matching_rows['MUMODEL'].tolist()
                mudif_values = [(mu - muref) for mu, muref in zip(mu_values, muref_values)]

                print(f" Computed MUDIF values: {mudif_values}")

                # Ensure dictionary has correct structure
                if (cid, idsurvey) not in mu_dict:
                    mu_dict[(cid, idsurvey)] = {"all": [], "by_fitres": {}}

                # Store values by FITRES file and in a global list
                mu_dict[(cid, idsurvey)]["by_fitres"][fitres] = mudif_values
                mu_dict[(cid, idsurvey)]["all"].extend(mudif_values)

                # Debugging: Print dictionary updates
                print(f" Updated mu_dict[{(cid, idsurvey)}]['by_fitres'][{fitres}] = {mudif_values}")
    '''            
    # xxxxxxx OBDSOLETE xxxxxxxx
    print("\n=== Collected MUDIF values for CID pairs (All FITOPT Files) ===")
    for key, fitres_data in mu_dict.items():
        print(f"\n CID {key[0]}, IDSURVEY {key[1]}:")
        for fitres, mudif_vals in fitres_data["by_fitres"].items():
            print(f"   {fitres}: {mudif_vals}")
    '''
    
    # Compute covariance of 2 CID pairs
    cid_pairs = list(mu_dict.keys())
    if len(cid_pairs) >= 2:
        cid1, idsurvey1 = cid_pairs[0]
        cid2, idsurvey2 = cid_pairs[1]
        print('CID1',cid1,cid2)
        print('cid2', mu_dict[(cid1, idsurvey1)])
        mudif_values_1 = mu_dict[(cid1, idsurvey1)]["all"]
        mudif_values_2 = mu_dict[(cid2, idsurvey2)]["all"]

        # xxxxxxx OBDSOLETE xxxxxxxx
        if len(mudif_values_1) == len(mudif_values_2):
            covariance_matrix = np.cov(mudif_values_1, mudif_values_2)
            print('Covariance matrix \n', covariance_matrix)
            covariance = covariance_matrix[0, 1]
            print(f"\n Covariance between CID {cid1} and CID {cid2}: {covariance}")

            # Check matrix conditioning
            cond_number = np.linalg.cond(covariance_matrix)
            if cond_number > 1e10:
                print("Warning: Covariance matrix is poorly conditioned, inversion may be unstable.")
        else:
            print("\n Error: MUDIF lists have different lengths, covariance cannot be computed.")
    else:
        print("\n Error: Less than two CID pairs found, covariance cannot be computed.")

    # xxxxxxx OBDSOLETE xxxxxxxx
    return mu_dict


def get_cov_from_create_covariance(args):

    cov_dir = args.cov_dir
    cov_num = args.cov_num
    hd_file = os.path.join(cov_dir, "hubble_diagram.txt")
    covsys_file = os.path.join(cov_dir, f"covsys_{cov_num:03d}.txt.gz")

    print(f"\n Read create_covariance element from \n\t {covsys_file}")

    # --- Step 1: Read the Hubble diagram file ---
    # The hubble diagram file is assumed to be space-delimited with commented header lines.
    df_hd = pd.read_csv(hd_file, comment='#', delim_whitespace=True)
    #print("Columns:", df_hd.columns.tolist())
    #print("Data types:\n", df_hd.dtypes)
    #print(df_hd.head())

    df_hd[VARNAME_CID]      = df_hd[VARNAME_CID].astype(str).str.strip()
    df_hd[VARNAME_IDSURVEY] = df_hd[VARNAME_IDSURVEY].astype(int)

    # Find the row indices for each CID/IDSURVEY pair in your arguments.
    row_indices = []
    for cid, idsurvey in zip(args.cid_list, args.idsurvey_list):
        # Note: Ensure the data types match (CID as string, IDSURVEY as integer)
        idx = df_hd[(df_hd['CID'] == (cid)) & (df_hd['IDSURVEY'] == (idsurvey))].index
        #print('XXX idx =', idx)
        if idx.empty:
            print(f"No row found for CID {cid}, IDSURVEY {idsurvey} in the hubble diagram.")
        else:
            # We take the first matching row
            row_indices.append(idx[0])
    
    print("\t Found row indices in hubble diagram:", row_indices)
    
    # Check that you have at least two indices to compare.
    if len(row_indices) < 2:
        print("Not enough rows found to extract a covariance element.")
        return None

    # --- Step 2: Read the flattened covariance matrix file ---
    # The file is assumed to be gzipped, with the first line as the dimension (e.g. 1802)
    # and subsequent lines as the matrix elements (row-major order).
    with gzip.open(covsys_file, 'rt') as f:
        # Read the first non-empty line which gives the dimension.
        first_line = f.readline().strip()
        if not first_line:
            print("Covariance file is empty or does not contain dimension info.")
            return None
        try:
            dim = int(first_line)
        except ValueError:
            print("Unable to parse the dimension from the covariance file:", first_line)
            return None

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
        print("Warning: Expected", dim * dim, "elements, but got", len(cov_elements))
    
    # Convert the list of elements to a NumPy array and reshape it.
    cov_matrix = np.array(cov_elements).reshape(dim, dim)
    print("Covariance matrix shape:", cov_matrix.shape)
    
    # --- Step 3: Extract the covariance entry for the two selected rows ---
    # For example, if you have two indices (row_i, row_j):
    row_i = row_indices[0]
    row_j = row_indices[1]
    covariance_value = cov_matrix[row_i, row_j]
    print(f"         Covariance between row {row_i} and row {row_j}: {covariance_value}")
    print(f"Reversed Covariance between row {row_j} and row {row_i}: {covariance_value}")
    print(f"")

    return covariance_value




# ============================================================    
# ============================================================    

if __name__ == "__main__":

    args = get_args()
    args = parse_cidpair(args)
    args = read_covinfo(args)

    # read scale per systematic
    sys_scale_list = get_sys_scale(args)
    
    # compute cov(cid0,cid1)
    dff, sum_cov_terms = compute_scaled_covpair(args, sys_scale_list)
    
    # write info to file
    dff.to_csv("Covariance_deltas.csv",  sep=' ', index=False)

    # fetch cov(cid0,cid1) from output of create_covariance
    cov_from_create_covariance = get_cov_from_create_covariance(args)

    comment  = "# covariance = (delta1 * sys scale) * (delta2 * sys scale)\n"
    comment2 = "# delta_i = MU_i - MU0 " 
    print('\n', comment)
    print( "",comment2, '\n')
    print(dff,'\n')
    print('File saved as Covariance_deltas.csv')
    
    cov_string = f"cov({args.cid_list[0]}, {args.cid_list[1]})"
    print(f"Computed for crosscheck     : {cov_string} = {sum_cov_terms}")
    print(f"Read from create_covariance : {cov_string} = {cov_from_create_covariance}")
    
    ratio = sum_cov_terms/cov_from_create_covariance
    print(f"Ratio(computed/read) : {ratio}")

    # === END: ===
