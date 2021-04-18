#!/usr/bin/env python
#
# Created April 2021
# For fakes and/or SNANA sim, create table of efficiency vs. redshift.
# (ignores real data)
# User inputs name of output LCFIT directory created from submit_batch_jobs.
#
# Output eff table is vs. redshift, and includes
#   eff_found = efficiency of having data file -> found it
#   eff_lcfit = efficiency for passing SNANA cuts and LCFIT.
#
# ========================================

import os, sys, argparse, glob, yaml, math
import numpy as np
from   argparse import Namespace
import pandas as pd

JOBNAME_SNANA   = "snana.exe"
JOBNAME_COMBINE = "combine_fitres.exe"

# name of subdir created under lcfit_dir folder
SUBDIR_EFF = "EFFICIENCY"

# table names for effic numerator
TABLE_NAME_SNANA  = "SNANA"   # simple detection eff
TABLE_NAME_FITRES = "FITRES"  # effic of detection + decent LC fit
TABLE_NAME_LIST   = [ TABLE_NAME_SNANA, TABLE_NAME_FITRES ]

ITABLE_SNANA  = 0
ITABLE_FITRES = 1
NTABLE        = 2

# use these "combined" table variables to measure effic vs. z
zVARNAME_DENOM = "ZCMB"
zVARNAME_NUMER = [ "zCMB_2", "zCMB_3" ]  # evt lost if value less than 0

# subset of variables to keep when combining tables.
VARLIST_DENOM    = "ZCMB,TRESTMIN,TRESTMAX"
VARLIST_NUMER    = "zCMB"

# ========================================

def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP on input config file"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of directory with LCFIT output"
    parser.add_argument("lcfit_dir", help=msg, nargs="?", default=None)

    msg = "Trestmin denominator cut (w.r.t. peak, default=0)"
    parser.add_argument("--trestmin", help=msg, nargs="?", 
                        type=float, default=0.0)
    msg = "Trestmax denominator cut (w.r.t. peak, default=10 days)"
    parser.add_argument("--trestmax", help=msg, nargs="?", 
                        type=float, default=10.0)

    msg = "max redshift to print lost events in FITRES table"
    parser.add_argument("--zmax_lost", help=msg, nargs="?", 
                        type=float, default=0.15)

    msg = "do NOT remove temp tables"
    parser.add_argument("--noclean", help=msg, action="store_true")

    args = parser.parse_args()

    if len(sys.argv) == 1 :
        parser.print_help()
        sys.exit()

    return args

    # end get_args

def read_yaml(input_file):
    input_lines = []
    with open(input_file, 'r') as f :
        for line in f:
            input_lines.append(line)

    config_yaml = yaml.safe_load("\n".join(input_lines))
    return config_yaml
    # end read_yaml

def read_lcfit_info(lcfit_dir):
    
    lcfit_dir        = os.path.expandvars(lcfit_dir)
    merge_info_file  = f"{lcfit_dir}/MERGE.LOG"
    submit_info_file = f"{lcfit_dir}/SUBMIT.INFO"

    if os.path.exists(merge_info_file):
        print(f"  Read {merge_info_file} ")
        merge_info = read_yaml(merge_info_file)
    else:
        fatal_error(f"  Could not find {merge_info_file}" )

    if os.path.exists(submit_info_file):
        print(f"  Read {submit_info_file} ")
        submit_info = read_yaml(submit_info_file)
    else:
        sys.exit(f"\n ERROR: Could not find \n\t{submit_info_file}")

    sys.stdout.flush() 
    return merge_info, submit_info


def get_zbins(config):
    zmin = 0.0
    zmax = 2.0
    zbin = 0.10
    nzbin = int((zmax-zmin)/zbin)
    zbins = np.linspace(zmin, zmax, nzbin+1)
    return nzbin, zbins
    # end get_zbins

def make_eff_outdir(config):
    
    
    outdir = f"{config.args.lcfit_dir}/{SUBDIR_EFF}"
    if os.path.exists(outdir) :
        cmd_rm = f"rm -r {outdir}"
        os.system(cmd_rm)
    # - - - 

    cmd_mkdir = f"mkdir {outdir}"
    os.system(cmd_mkdir)

    print(f"  Created outdir: {outdir}")
    sys.stdout.flush() 

    return outdir
    # end make_eff_outdir

def eff_driver(config, row):

    # input row is from MERGE.LOG
    # table_name is either SNANA or FITRES
    # After combining dumpall and table from snana, 
    # Efficiency is NEVT(zCMB_2>0)/NEVT(ZCMB>0)

    submit_info = config.submit_info
    lcfit_dir   = config.args.lcfit_dir
    outdir      = config.outdir
    noclean     = config.args.noclean
    clean       = not noclean

    state   = row[0]  # should be DONE
    version = row[1]  # name of data version
    fitopt  = row[2]  # e..g, FITOPT000

    print('')
    print(f"       {version}_{fitopt} ")
    sys.stdout.flush() 

    # make sure job finished without error
    STATE_DONE = "DONE"
    if state != STATE_DONE :
        msgerr  = f"  state = {state} (should be '{STATE_DONE}') \n"
        msgerr += f"  check row = \n\t{row} \n   of MERGE.LOG file\n"
        fatal_error(msgerr)

    # - - - - - - - -
    # run quick snana.exe job to get version info: 
    # path for VERSION, and data type.
    version_info = get_version_info(config,version)
    version_path = version_info['SNDATA_PATH']
    datatype     = version_info['DATATYPE']
    if datatype == 'DATA' : return  # can't get eff for real data

    # locate DUMPALL file with info for every generated/overlaid event;
    # abort if not found
    dumpall_file = f"{version_path}/{version}.DUMP"
    exist = check_filegz_exists(dumpall_file,"Every generated event")
    print(f" Include DUMPALL  table from {version}")

    # start string list of table files for combine_fitres
    table_file_string = f"{dumpall_file} "
    n_table_file = 1

    # - - - - - - - 
    # check which tables exist      
    table_file_prefix = f"{lcfit_dir}/{version}/{fitopt}"
    table_exist_list = []
    for table_name in TABLE_NAME_LIST :
        table_file = f"{table_file_prefix}.{table_name}"
        exist      = check_filegz_exists(table_file,None)
        table_exist_list.append(exist)
        if exist :
            table_file_string += f"{table_file} "
            n_table_file += 1
            print(f" Include {table_name:<8} table from {version}")
            sys.stdout.flush() 

    config.table_exist_list = table_exist_list
    
    if n_table_file == 1:
        msgerr  = f"  Found no analysis tables; \n" \
                  f"  Expecting 1 or more from {TABLE_NAME_LIST}\n"
        fatal_error(msgerr)

    # - - - - - - - - -
    # get data frame for combined table
    df = read_combined_table(table_file_string,f"{version}_{fitopt}", clean)
    nrow_tot = len(df)

    # apply denominator cuts
    nrow_cuts, df = apply_denom_cuts(config,df)

    # table error returns zero rows -> bail
    if nrow_cuts == 0 : 
        print("  Table ERROR -> skip efficiency calc")
        return 0

    print(f"   Nrow = {nrow_tot} (total) -> {nrow_cuts} (cuts)")
    sys.stdout.flush() 

    # compute efficiencies vs. z
    eff_dict = eff_compute(config,df)

    # write eff vs. z to file
    eff_file = f"{outdir}/EFFz_{version}_{fitopt}.DAT"
    write_eff_file(config, eff_file, eff_dict)

    return 1

    # end eff_driver

def eff_compute(config,df):

    # Return eff_found (vs z) from SNANA table.
    # Return eff_lcfit (vs z) from FITRES table.
    # If a table is missing, return array of -9.
    # Logic is tricky to allow either table, or both;
    # beware that 2nd z variable is zCMB_2 regardless of whether
    # 2nd table is SNANA or FITRES.

    zbins = config.zbins
    nzbin = config.nzbin
    table_exist_list = config.table_exist_list # vs. SNANA, FITRES
    
    # init effiency arrays to -9 in case one of the tables doesn't exist.
    eff_found     = [ -9 ] * nzbin
    eff_lcfit     = [ -9 ] * nzbin
    eff_found_err = [ -9 ] * nzbin
    eff_lcfit_err = [ -9 ] * nzbin

    nall_unbinned = df[zVARNAME_DENOM]
    nall_digi     = np.digitize(nall_unbinned, zbins)
    nall_binned   = np.bincount(nall_digi)[1:]
    nzbin_eff     = len(nall_binned)
    zcen          = ((zbins[1:] + zbins[:-1])/2)[0:nzbin_eff]
    
    dump_flag = False

    if dump_flag :
        print(f" xxx ----------------------------------- ")
        print(f" xxx zbins = {zbins} ")
        print(f" xxx zcen  = {zcen}  ")
        print(f" xxx nall_binned = {nall_binned} ")

    itmp = 0
    for exist,table_name in  zip(table_exist_list,TABLE_NAME_LIST):
        if exist :
            zVARNAME = zVARNAME_NUMER[itmp] ; itmp+=1
            nana_unbinned  = df[zVARNAME]
            nana_digi      = np.digitize(nana_unbinned, zbins)
            nana_binned    = np.bincount(nana_digi,None,nzbin_eff)[1:]

            # compute efficiency
            effz           = np.divide(nana_binned,nall_binned)

            # compute binomial uncertainty
            ntmp      = abs(nall_binned-nana_binned)*nana_binned
            effz_err  = (ntmp/nall_binned**3)**0.5

            # get list of missing SN at low-z (FITRES table only)
            if table_name == TABLE_NAME_FITRES :
                zmax = config.args.zmax_lost
                df_lowz_lost = \
                    df[(df[zVARNAME_DENOM]<zmax) & (df[zVARNAME]<0.) ]

            if dump_flag :
                print(f" xxx table = {table_name}")
                print(f" xxx   nana_binned = {nana_binned} ")
                print(f" xxx   effz  = {effz} ")
                print(f" xxx   err   = {effz_err} ")

            if table_name == TABLE_NAME_SNANA :
                eff_found = effz ; eff_found_err = effz_err
            if table_name == TABLE_NAME_FITRES :
                eff_lcfit = effz ; eff_lcfit_err = effz_err

    # - - - - - - - - - - - - - - - 
    # load output dictionary
    eff_dict = {
        'nzbin'         : nzbin_eff ,
        'zcen'          : zcen ,
        'nall'          : nall_binned ,
        'eff_found'     : eff_found ,
        'eff_found_err' : eff_found_err ,
        'eff_lcfit'     : eff_lcfit ,
        'eff_lcfit_err' : eff_lcfit_err,
        'df_lowz_lost'  : df_lowz_lost
    }

    return eff_dict

    # end eff_compute

def write_eff_file(config, eff_file, eff_dict):

    args = config.args
    print(f"   write {eff_file}")
    sys.stdout.flush() 

    line_header = "VARNAMES: z  n_all  eff_found  " \
                  "eff_found_err  eff_lcfit  eff_lcfit_err"

    nzbin         = eff_dict['nzbin']
    zcen          = eff_dict['zcen']
    nall          = eff_dict['nall']
    eff_found     = eff_dict['eff_found']
    eff_found_err = eff_dict['eff_found_err']
    eff_lcfit     = eff_dict['eff_lcfit']
    eff_lcfit_err = eff_dict['eff_lcfit_err']
    df_lowz_lost  = eff_dict['df_lowz_lost']

    f = open(eff_file,"wt")

    f.write(f"# eff_found => found SN and created light curve. \n")
    f.write(f"# eff_lcfit => passes SNANA cuts and SALT2 LCFIT\n")
    f.write(f"# Trestmin < {args.trestmin} \n")
    f.write(f"# Trestmax > {args.trestmax} \n")
    f.write(f"\n")

    f.write(f"{line_header} \n")
    for z, n,eff0,err0, eff1,err1 in \
        zip(zcen,nall,eff_found,eff_found_err,eff_lcfit, eff_lcfit_err) :
        line = f"ROW: {z:.3f}   {n:5d}   " \
               f"{eff0:.4f} {err0:.4f}  {eff1:.4f} {err1:.4f}"
        f.write(f"{line} \n")
        
    # - - - - - - - -
    # write lowz-lost list from FITRES table
    nrow=0; nrow_max = 50
    zmax_lost = config.args.zmax_lost
    f.write(f"\n\n# Missing FITRES events with z < {zmax_lost}: \n")
    f.write(f"#           CID    zCMB   Trestmin  Trestmax \n")
    for index, row in df_lowz_lost.iterrows() :
        nrow += 1
        if nrow > nrow_max: 
            f.write(f"#    (stop writing after {nrow_max} events)\n")
            break
        cid = row['CID']
        z   = row['ZCMB']
        trestmin = row['TRESTMIN']
        trestmax = row['TRESTMAX']
        f.write(f"#    {cid:>12}  {z:.4f}  {trestmin:6.1f} {trestmax:6.1f}\n")


    f.close()
    # end write_eff_file

def read_combined_table(table_file_string, stringid, clean):

    # combine tables in "table_file_string", then read it.
    # Note that combine_fitres program automatically checks for .gz.
    # Function returns data frame.
    # Input stringid is part of temp file name so that in case of crash
    # it's a bit easier to know where the temp file is from.

    n_table_file = len(table_file_string.split())
    print(f"   Combine {n_table_file} table files ") 
    sys.stdout.flush() 

    prefix   = f"TEMP_COMBINED_{stringid}"
    log_file = f"{prefix}.LOG"
    combined_file = f"{prefix}.TEXT"
    cmd = f"{JOBNAME_COMBINE} {table_file_string} T " \
          f"--varnames {VARLIST_DENOM},{VARLIST_NUMER} " \
          f"--outprefix {prefix} > {log_file} "
    os.system(cmd)

    # read combined table into pandas data frame
    df = pd.read_csv(combined_file, comment="#", delim_whitespace=True)

    # remove TEMP table file
    if clean:
        cmd_rm = f"rm {log_file} {combined_file}"
        os.system(cmd_rm)

    return df
    # end read_combined_table

def apply_denom_cuts(args,df):

    trestmin = config.args.trestmin
    trestmax = config.args.trestmax

    # check that requested variables exist
    nerr = 0
    varlist_check = list(VARLIST_DENOM.split(',')) 
    for varname in varlist_check :
        if varname not in df:
            nerr += 1
            print(f"  ERROR: missing required column '{varname}' ")

    if nerr > 0 : return 0, df

    df = df[ df['TRESTMIN'] < trestmin ]
    df = df[ df['TRESTMAX'] > trestmax ]
    
    nrow_cuts = len(df)
    return nrow_cuts, df

    # end apply_denom_cuts

def check_filegz_exists(file_name, comment_abort ):
    # return true if file_name exists, or if file_name.gz exists
    file_namegz = f"{file_name}.gz"
    found = False
    if os.path.exists(file_name):   found = True
    if os.path.exists(file_namegz): found = True

    if not found and comment_abort is not None :
        msgerr = f" could not find required file :\n" \
                 f"  {file_name}\n" \
                 f" File comment: {comment_abort} \n"
        fatal_error(msgerr)
        
    return found
    # end check_filegz_exists

def get_version_info(config,version):

    # run snana job with 'GETINFO' arg to get path to folder
    # read and return YAML output in log file.
    # Beware that tail -10 is fragile.

    private_data_path = config.submit_info['PRIVATE_DATA_PATH']
    version_arg       = version  # default if SIM or in $SNDATA_ROOT/lcmerge

    # if private_data_path is set, check if version is there
    if private_data_path is not None :
        full_path   = os.path.expandvars(private_data_path)
        version_tmp = f"{private_data_path}/{version}" 
        if os.path.exists(f"{full_path}/{version}"): 
            version_arg = version_tmp   # include ENV for print and snana.exe

    version  = os.path.basename(version_arg)
    log_file = f"TEMP_GETINFO_{version}.log"
    cmd = f"{JOBNAME_SNANA} GETINFO {version_arg} | tail -10 > {log_file} "
    os.system(cmd)    
    version_info = read_yaml(log_file)

    cmd_rm = f"rm {log_file}"
    os.system(cmd_rm)
    return version_info

    # end get_version_info

def fatal_error(msgerr):
    print(f"\n FATAL ERROR:")
    sys,exit(msgerr)

# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    config      = Namespace()
    config.args = get_args()

    # get info from MERGE.LOG
    merge_info,submit_info = read_lcfit_info(config.args.lcfit_dir)
    config.merge_info  = merge_info
    config.submit_info = submit_info

    # get z bins
    config.nzbin, config.zbins = get_zbins(config)

    # create outdir under folder
    config.outdir = make_eff_outdir(config)

    print('')

    nrow_tot = 0;  nrow_eff=0
    for row in merge_info['MERGE']:
        nrow_tot += 1
        nrow_eff += eff_driver( config, row)

    print(f"\n Done computing eff table for " \
          f"{nrow_eff} of {nrow_tot} MERGE.LOG rows.\n")

    # END
