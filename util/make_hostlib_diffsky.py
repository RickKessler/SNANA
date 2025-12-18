#!/usr/bin/env python
#
# Created Dec 2025 by R.Kessler
# Created HOSTLIB from diffksy catalog
#
# https://opencosmo.readthedocs.io/en/stable/main_api.html
#
# Ref code for Roman-DESC hostlib:
#   $ELASTICC_ROOT/roman+desc/inputs_hostlib/make_hostlib_from_GCRCat.py
#
# ==============================================

import os, argparse, logging, shutil, datetime, time, glob, random
import re, yaml, sys, gzip, math, subprocess
import pandas as pd

import opencosmo as oc
import numpy as np

# =======

KEY_CAT_DIR              = "CAT_DIR"
KEY_HOSTLIB_VARNAMES_MAP = "HOSTLIB_VARNAMES_MAP"
KEY_HOSTLIB_FILE         = "HOSTLIB_FILE"
KEY_CUTWIN               = "CUTWIN"
KEY_MAG_5SIG             = "MAG_5SIG"
KEY_HOSTLIB_GALID        = "GALID"

N_UPDATE_STDOUT = 10000  # stdout update after this many rows written

RANSEED = 234237905   # used to ensure same random subsets for smaller HOSTLIBs

REF_DICT = {
    'REF1': [ 'Heitmann+ 2021 (N-body sim)',
              'https://ui.adsabs.harvard.edu/abs/2021ApJS..252...19H' ]
    }


# hard-wire sersic indices for each galaxy component
SERSIC_INDEX_DICT = {  'n_disk': 1.0,    'n_bulge': 4.0 }


# ============================
def setup_logging():
    #logging.basicConfig(level=logging.DEBUG,
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)8s |%(lineno)4d]  %(message)s")
    # end setup_logging

def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP menu for config options"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of the yml config file to run"
    parser.add_argument("config_file", help=msg, nargs="?", default=None)

    msg = "print data columns containing substring (or all to list all columns), then quit"
    parser.add_argument("-p", "--print_column_substring", help=msg, type=str, default=None)


    # parse it
    args = parser.parse_args()

    if args.HELP : 
        print_help_menu()

    if args.config_file is None:
        sys.exit(f"\n ERROR: missing config_file input")
        
    return args
    # end get_args

def print_help_menu():

    help_menu = f"""

CAT_DIR: [catalog_dir]

# Define 1 more HOSTLIBs with different numbers of galaxies so that
# smaller HOSTLIBs can be used in the sim for testing with much
# shorter read/init time.
# The largest HOSTLIB defines the nominal sample, and the
# smaller HOSTLIBs are random subsets of nominal.

HOSTLIB_FILE:
- TEST_DIFFSKY_LSST_LARGE.HOSTLIB   100000
- TEST_DIFFSKY_LSST_MEDIUM.HOSTLIB   10000  # this smaller hostlib is a random subset
- TEST_DIFFSKY_LSST_SMALL.HOSTLIB     1000  # another evern smaller hostlib is a random subset


#   diffsky       HOSTLIB        HOSTLIB
#   varname       varname        format
HOSTLIB_VARNAMES_MAP:
-  galaxy_id       GALID           15d
-  redshift_true   ZTRUE_CMB       6.4f
-  ra              RA_GAL          11.6f
-  dec             DEC_GAL         11.6f
-  beta_disk       a0_Sersic       5.2f
-  alpha_disk      b0_Sersic       5.2f
#
-  lsst_u          u_obs           6.3f
-  lsst_g          g_obs           6.3f
-  lsst_r          r_obs           6.3f
-  lsst_i          i_obs           6.3f
-  lsst_z          z_obs           6.3f
-  lsst_y          y_obs           6.3f
#
-  lsst_u_err      u_obs_err         6.3f   # internally computed if MAG_5SIG are provided
-  lsst_g_err      g_obs_err         6.3f
-  lsst_r_err      r_obs_err         6.3f
-  lsst_i_err      i_obs_err         6.3f
-  lsst_z_err      z_obs_err         6.3f
-  lsst_y_err      y_obs_err         6.3f


MAG_5SIG:
- lsst_u  25  # use this 5 sigma depth to compute [band]_err columns
- lsst_g  27
- lsst_r  27
- lsst_i  27
- lsst_z  26
- lsst_y  25


# define selection cuts using diffsky-catalog column names
CUTWIN:
- logsm_obs        8   13  # logmass cut
- redshift_true    0  1.8
- lsst_r           0   28

    """
    print(f"{help_menu}\n")
    
    sys.exit()    
    # end print_help_menu
    
def read_yaml(path):
    path_expand = os.path.expandvars(path)
    logging.info(f"Reading YAML from {path_expand}")
    with open(path_expand) as f:
        return yaml.safe_load(f.read())
    # end read_yaml


def addcol_sersic_indices(cat_inp, config):

    n_row   = len(cat_inp)
    cat_out = cat_inp
    HOSTLIB_VARNAMES_MAP = config[KEY_HOSTLIB_VARNAMES_MAP]

    addcol_dict = {}
    string_sersic_list = list( SERSIC_INDEX_DICT.keys() )
    for row in HOSTLIB_VARNAMES_MAP:
        colname = row.split()[0]
        if colname in string_sersic_list:
            s_index = SERSIC_INDEX_DICT[colname] 
            addcol_dict[f"{colname}"] = np.array([ s_index ] * n_row)
            logging.info(f"\t append  {colname} = {s_index} ")
    # - - - -
    
    if len(addcol_dict) > 0 :
        cat_out = cat_out.with_new_columns(**addcol_dict)
        
    #sys.exit(f"\n xxx addcol_dict = \n{addcol_dict} ")

    return cat_out

def check_Sersic_columns(config):

    # sanity check to make sure that a,b,n are defined for each sersic component,
    # and that number of weights is n_profile - 1.
    # Check HOSTLIB names that well-defined.

    sersic_list_dict = get_Sersic_list(config)

    val_list_a = sersic_list_dict['a']
    parnames_a = [ f"a{x}_Sersic" for x in val_list_a ]
    for par_Sersic, val_list in sersic_list_dict.items() :
        if par_Sersic == 'w' : continue
        if val_list != val_list_a:
            parnames = [ f"{par_Sersic}{x}_Sersic" for x in val_list ]
            sys.exit(f"\n ERROR: {par_Sersic} has Sersic components {parnames} \n" \
                     f"\t but a has {parnames_a} \n" \
                     f"Each Sersic component must have a, b, n")

    # - - - - -
    num_w = len(sersic_list_dict['w'])
    num_a = len(sersic_list_dict['a'])

    if num_a-1 != num_w :
        sys.exit(f"\n ERROR: expected {num_a-1} Sersic weights (w#_Sersic), but found {num_w}\n")
        
    return

def get_Sersic_list(config):

    # return dictionary of lists for a, b, n, w
    # where a_list = [0,1] --> a0_Sersic and a1_Sersic are defined

    HOSTLIB_VARNAMES_MAP = config[KEY_HOSTLIB_VARNAMES_MAP]
    
    sersic_list_dict = { 'a': [],  'b': [],  'n': [],  'w':[]}

    for row in HOSTLIB_VARNAMES_MAP:
        varname_hostlib = row.split()[1]
        if '_Sersic' in varname_hostlib:
            par_Sersic = varname_hostlib[0:1] # a or b or n or w
            j      = varname_hostlib[1:2]  # 2nd char is 0, 1, etc 
            sersic_list_dict[par_Sersic].append(j)

    # sort each list in case user inputs strange ordering
    for par_Sersic, val_list in sersic_list_dict.items() :
        sersic_list_dict[par_Sersic] = sorted(val_list)
    
    return sersic_list_dict


def addcol_mag_errors(cat_inp, config):

    # if MAG_5SIG is defined, add columns with mag errors
    
    cat_out = cat_inp
    MAG_5SIG = config.setdefault(KEY_MAG_5SIG,None)

    if not MAG_5SIG:  return galaxy_cat_append

    logging.info('')
    logging.info('Append mag error columns:')

    t0 = time.time()
    
    mag_err_dict = {}
    for row in MAG_5SIG:
        str_band     = row.split()[0]
        str_band_err = str_band + '_err'
        m5sig        = float(row.split()[1])
        powm5        = 0.2*math.pow(10.0,-0.2*m5sig)

        logging.info(f"\t append {str_band_err} using m5sig = {m5sig}")

        mag_np = cat_out.select([str_band]).get_data().value
        mag_err_np     = powm5*np.power(10.0,+0.2*mag_np)
        mag_err_dict[f"{str_band_err}"] = mag_err_np

    # append the mag err columns
    cat_out = cat_out.with_new_columns(**mag_err_dict)    

    print_proc_time(t0, "ADD_MAGERR", None)
    
    return cat_out

def read_galaxy_cat(args, config):

    t0 = time.time()
    cat_dir  = config[KEY_CAT_DIR]
    
    if not os.path.exists(cat_dir):
        sys.exit(f"\n ERROR: cannot find galaxy catalog dir:\n\t{cat_dir}")

    # fragile-alert reading catalog files with specific prefix
    wildcard = f"{cat_dir}/lc_cores-*.hdf5"  # fragile alert; should read standard list file
    catalog_files = sorted( glob.glob(wildcard) )

    if args.print_column_substring:
        catalog_files = [ catalog_files[0] ]  # keep only first file for listing columns
        
    n_file = len(catalog_files)
    
    logging.info(f"Read {n_file} galaxy catalog files from {cat_dir}")

    ds = oc.open(*catalog_files)
    n_row = len(ds)
    logging.info(f"Done reading dataset with {n_row:,} rows")
    print_proc_time(t0,"READ_CATALOG", n_row)
        
    if args.print_column_substring:
        print_columns(args.print_column_substring, ds.columns)        

        
    # - - - - - - 
    return ds

    # end read_galaxy_cat

def apply_cuts(cat_inp, config):

    t0 = time.time()
    cat_out = cat_inp

    n_row_inp = len(cat_inp)
    CUTWIN  =  config.setdefault(KEY_CUTWIN,None)
    if CUTWIN is None: return cat_out

    logging.info(f"")
    logging.info(f"Apply cuts:")
    for row in CUTWIN:
        cutvar = row.split()[0]
        cutmin = float(row.split()[1])
        cutmax = float(row.split()[2])
        logging.info(f"\t Apply cut {cutmin} < {cutvar} < {cutmax} ")
        cat_out = cat_out.filter(oc.col(cutvar)<cutmax).filter(oc.col(cutvar)>cutmin)

    n_row_out = len(cat_out)
    logging.info(f" Cuts reduce {n_row_inp:,} rows to {n_row_out:,} rows")

    print_proc_time(t0, "APPLY_CUTS", None)
    
    return cat_out   # end apply_cuts

def print_columns(substring, column_list):

    # inputs:
    #   substring      : print columns containing this substring
    #   column_list    : list of all available columns
    
    if substring == 'all' :
        filtered_list = column_list
    else:
        filtered_list = [x for x in column_list if substring in x]
               
    print(f"\n DIFFSKY Columns: \n{filtered_list}")
    sys.exit("\nDone listing columns")
    
    return   # end print_columns

def parse_varname_map(config):

    logging.info('')
    logging.info(f"Parse {KEY_HOSTLIB_VARNAMES_MAP}" )

    # read map from user input and covert to dictionarys
    HOSTLIB_VARNAMES_MAP = config[KEY_HOSTLIB_VARNAMES_MAP]

    hostlib_varname_dict = {}
    hostlib_format_dict  = {}

    HOSTLIB_VARNAMES_STRING = ''
    FOUND_GALID = False
    
    for row in HOSTLIB_VARNAMES_MAP:
        logging.info(f"\t {row}")
        tmp_list = row.split()
        varname_diffsky = tmp_list[0]
        varname_hostlib = tmp_list[1]
        format_hostlib  = tmp_list[2]
        hostlib_varname_dict[varname_diffsky] = varname_hostlib
        hostlib_format_dict[varname_diffsky]  = format_hostlib
        HOSTLIB_VARNAMES_STRING += f"{varname_hostlib} "

        if varname_hostlib == KEY_HOSTLIB_GALID:  FOUND_GALID = True
            
    logging.info('')

    config['hostlib_varname_dict']    = hostlib_varname_dict
    config['hostlib_format_dict']     = hostlib_format_dict
    config['HOSTLIB_VARNAMES_STRING'] = HOSTLIB_VARNAMES_STRING
    config['FOUND_GALID']             = FOUND_GALID
    
    #sys.exit(f"\n xxx hostlib_varname_dict = \n{hostlib_varname_dict} \n\n xxx hostlib_format_dict =\n{hostlib_format_dict}")
    return config


def open_hostlib_files(config):

    HOSTLIB_VARNAMES_STRING = config['HOSTLIB_VARNAMES_STRING'] 
    hostlib_dict = {}

    ngal_list = []
    HOSTLIB_FILE_LIST = config[KEY_HOSTLIB_FILE]
    for row in HOSTLIB_FILE_LIST:
        hostlib_file = row.split()[0]
        ngal         = int(row.split()[1])
        fp = open(hostlib_file,"wt")
        write_hostlib_header(fp, ngal, config)
        hostlib_dict[hostlib_file] = [ ngal, fp ]
        ngal_list.append(ngal)

    nmax_gal  = max(ngal_list)
    n_hostlib = len(HOSTLIB_FILE_LIST)
    
    config['hostlib_dict'] = hostlib_dict
    config['nmax_gal']     = nmax_gal
    config['n_hostlib']    = n_hostlib

    if nmax_gal > 500000 :
        global N_UPDATE_STDOUT; N_UPDATE_STDOUT = 100000
        
    return config

def write_hostlib_header(fp, ngal, config):

    # write DOCUMENTATION block and VARNAMES

    HOSTLIB_VARNAMES_STRING =  config['HOSTLIB_VARNAMES_STRING'] 
    HOSTLIB_VARNAMES_MAP    =  config[KEY_HOSTLIB_VARNAMES_MAP]
    CUTWIN                  =  config.setdefault(KEY_CUTWIN,None)
    MAG_5SIG                =  config.setdefault(KEY_MAG_5SIG,None)
    
    cat_dir      = config[KEY_CAT_DIR]
    cat_dir_base = os.path.basename(cat_dir)
    
    USERNAME     = os.environ['USER']
    HOSTNAME     = os.uname()[1]
    
    tnow        = datetime.datetime.now()
    TSTAMP      = f"{tnow.year:04d}-{tnow.month:02d}-{tnow.day:02d}"
    CODE        = os.path.basename(sys.argv[0])
                           
    fp.write(f"DOCUMENTATION: \n")
    fp.write(f"  PURPOSE: host galaxy library for SNANA simulation \n")

    fp.write(f"  REF: \n")
    for ref, item_list in REF_DICT.items():
        AUTHOR = item_list[0]
        ADS    = item_list[1]
        fp.write(f"  - AUTHOR: {AUTHOR}\n")
        fp.write(f"    ADS:    {ADS}\n")
        
    fp.write(f"  USAGE_KEY:  HOSTLIB_FILE \n")
    fp.write(f"  USAGE_CODE: snlc_sim.exe \n")
    fp.write(f"  NOTES: \n")
    fp.write(f"  - extracted from diffsky catalog {cat_dir_base} \n")
    fp.write(f"  - ngal = {ngal} \n")

    fp.write(f"\n")
    fp.write(f"  HOSTLIB_VARNAMES_MAP: \n")
    for row in HOSTLIB_VARNAMES_MAP:
        fp.write(f"  - {row} \n")


    if MAG_5SIG:
        fp.write(f"\n")
        fp.write(f"  MAG_5SIG: \n")
        for row in MAG_5SIG:
            fp.write(f"  - {row} \n")


    if CUTWIN :
        fp.write(f"\n")
        fp.write(f"  CUTWIN: \n")
        for row in CUTWIN:
            fp.write(f"  - {row} \n")

    
    fp.write(f"\n")
    fp.write(f" PROVENANCE: \n")
    fp.write(f"  - creation code  {CODE} \n")
    fp.write(f"  - creation date  {TSTAMP} \n")
    fp.write(f"  - created by user={USERNAME} on node={HOSTNAME}  \n")
    fp.write(f"DOCUMENTATION_END: \n")
    fp.write(f"\n")


    fp.write(f"VARNAMES: {HOSTLIB_VARNAMES_STRING} \n")

    fp.flush()
    
    return  # end write_hostlib_header

def get_random_subset(args, config, galaxy_cat):

    t0 = time.time()
    nmax_gal             = config['nmax_gal']
    hostlib_varname_dict = config['hostlib_varname_dict']
    CUTWIN               = config[KEY_CUTWIN]
    
    cat_var_list = list(hostlib_varname_dict.keys())

    logging.info(f"")
    
    # fetch astropy table
    n_row     = 0
    ntake_gal = nmax_gal   # start with guess
    while n_row < nmax_gal:
        logging.info(f"Select random subset of {ntake_gal:,} rows ...")

        galaxy_subset = galaxy_cat \
            .take(ntake_gal, at="random") \
            .select(cat_var_list) \
            .get_data() 
            
        n_row = len(galaxy_subset)
        if n_row < nmax_gal:
            comment = f"too few; try again"
        else:
            comment = f">= {nmax_gal:,}; stop selecting rows"
            
        logging.info(f"\t --> {n_row:,} rows selected  ({comment})")
        ntake_gal = int(1.1 * ntake_gal * nmax_gal / n_row)
        
    # - - - -
        
    #.filter(oc.col("redshift_true")<.2) and oc.col("redshift_true")>.1)) \
    #.filter(oc.col("fof_halo_mass") > 1e13)
    
    # convert to data frame
    df = galaxy_subset.to_pandas()

    print_proc_time(t0, "RANDOM_SUBSET", None)
    
    return df

def get_hostlib_row_format(row_val_list , config):

    # for input row, return row for HOSTLIB with user-specified format per item

    hostlib_format_dict = config['hostlib_format_dict']
    var_list  = list(hostlib_format_dict.keys() )
    fmt_list  = list(hostlib_format_dict.values() )

    #sys.exit(f"\n xxx val_list = \n{row_val_list} \n xxx var_list = \n{var_list}  \n xxx fmt_list = \n {fmt_list}")
    row_hostlib = ''   
    for val, fmt in zip(row_val_list, fmt_list) :
        row_hostlib += f"{val:{fmt}} "

    #sys.exit(f"\n xxx row = \n{row_val_list} \n\t --> \n {row_hostlib} ")
    
    return row_hostlib  # end get_hostlib_row_format

def print_proc_time(t0, comment, N):
    t1 = time.time()
    dt = t1-t0

    key =  f"CPU({comment}):"

    msg = f"{key:<22} {dt:7.1f} sec"

    if N:
        rate = int(float(N)/dt)
        msg += f"  ({rate:,} per sec)"
        
    logging.info(f"\t {msg}")
    
    return

def write_hostlib_files(args, config, df_cat):

    t0 = time.time()
    
    hostlib_format_dict = config['hostlib_format_dict']
    hostlib_dict = config['hostlib_dict']
    FOUND_GALID  = config['FOUND_GALID']
    nmax_gal     = config['nmax_gal']
    n_hostlib    = config['n_hostlib']
    var_list  = list(hostlib_format_dict.keys() )
    row_list  = df_cat[var_list].values.tolist()
    nrow_tot  = len(df_cat)

    logging.info("")
    logging.info(f"Prepare writing {n_hostlib} HOSTLIB files with nmax_gal={nmax_gal:,}... ")
    
    # check for random subsets
    ngal_write_dict = {}
    random_row_mask = {}
    random.seed(RANSEED)   # ensure reproducible random rows
    for key, tmp in hostlib_dict.items():

        logging.info(f"\t prepare row mask for {key} ...")
        
        ngal = tmp[0]
        ngal_write_dict[key] = [ ngal, 0 ]  # n_write_expect   and n_write_actual

        frac_select = float(ngal / nrow_tot)
        random_row_mask[key] = np.random.choice([True, False], nrow_tot,
                                                p = [ frac_select, 1-frac_select] )
    # - - - -

    logging.info(f"  Start writing GAL rows ... ")
    nrow = -1
    nrow_wr = 0
    for item_list in row_list :  # each element of row_list is a list of items
        nrow += 1
        # convert comma-sep list into HOSTLIB-formatted row
        row_hostlib = get_hostlib_row_format(item_list,config)
        
        if not FOUND_GALID:
            GALID = 1000000 + nrow
            row_hostlib  = f"{GALID}  {row_hostlib}"

        WROTE_ROW = False
        for key, tmp in hostlib_dict.items():
            ngal = tmp[0]
            fp   = tmp[1]
            if random_row_mask[key][nrow] :
                fp.write(f"GAL:  {row_hostlib} \n")
                ngal_write_dict[key][1] += 1  # increment actual number of rows written
                WROTE_ROW = True
                
        if WROTE_ROW :
            nrow_wr += 1
            update = (nrow_wr % N_UPDATE_STDOUT == 0 or nrow_wr == nmax_gal) and nrow_wr > 0
            if update:
                logging.info(f"\t Finished writing {nrow_wr:12,} HOSTLIB rows")

    # summary expected and actual number of rows per output HOSTLIB
    logging.info('' )
    logging.info('Summary:' )
    for key in list ( hostlib_dict.keys() ):
        ngal_expect = ngal_write_dict[key][0]
        ngal_actual = ngal_write_dict[key][1]
        logging.info(f"  {key:32} : ngal[expect, actual] = {ngal_expect:10,}  {ngal_actual:10,}  ")

    # - - - -
    print_proc_time(t0, "WRITE_HOSTLIB", nmax_gal )
        
    return  # end write_hostlib_files


# ===================================================
if __name__ == "__main__":

    setup_logging()
    logging.info("# ========== BEGIN make_simsed_binaries ===============")

    t_start = time.time()
    args   = get_args() 
    config = read_yaml(args.config_file)

    # make sure Sersic column definitinos are consistent to avoid crash later
    check_Sersic_columns(config)
    
    # parse varname mapping between diffsky and hostlib
    config     = parse_varname_map(config)

    # open each hostlib file
    config       = open_hostlib_files(config)

    # fetch entire galaxy catalog
    galaxy_cat  = read_galaxy_cat(args, config)

    # apply cuts
    galaxy_cat  = apply_cuts(galaxy_cat, config)
        
    # check option to add and compute mag_error for each band
    galaxy_cat = addcol_mag_errors(galaxy_cat, config)

    # add sersic indices
    galaxy_cat = addcol_sersic_indices(galaxy_cat, config)
    
    # fetch data frame table for largest hostlib
    df_cat = get_random_subset(args, config, galaxy_cat)

    write_hostlib_files(args, config, df_cat)

    print_proc_time(t_start, "TOTAL", None)
    
    logging.info("Done.")
    # === END: ===
    
