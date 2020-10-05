#!/usr/bin/env python
#
# Created Oct 2020 [replace old perl script get_fitresValue.pl]
#
# Extract table value(s) for cid(s) and print to screen.
# Intended for visual debugging.
#

import os, sys, argparse
import pandas as pd

# =====================================
def get_args():
    parser = argparse.ArgumentParser()

    
    msg = "name of fitres file"
    parser.add_argument("-f", "--file", help=msg, type=str, default="")

    msg = "comma separated list of CIDs"
    parser.add_argument("-c", "--cid", help=msg, type=str, default="")

    msg = "comma separated list of variable names"
    parser.add_argument("-v", "--varname", help=msg, type=str, default="")

    args = parser.parse_args()
    return args

    # end get_args

def parse_inputs(inputs_dict):

    # parse comma-sep list of cid(s) and varname(s)

    cid_string = inputs_dict['args'].cid
    var_string = inputs_dict['args'].varname

    if cid_string == '' :
        sys.exit(f"\n ERROR: missing -c argument.")

    if var_string == '' :
        sys.exit(f"\n ERROR: missing -v argument.")

    cid_list   = cid_string.split(",")
    var_list   = var_string.split(",")
    
    inputs_dict['cid_list'] = cid_list
    inputs_dict['var_list'] = var_list

    # end parse_inputs

def read_fitres_file(inputs_dict):

    ffile    = inputs_dict['args'].file
    var_list = inputs_dict['var_list']

    var_list_local = ['CID'] + var_list
    df  = pd.read_csv(ffile, comment="#", delim_whitespace=True, 
                      usecols=var_list_local)
    df["CID"] = df["CID"].astype(str)
    df        = df.set_index("CID", drop=False)
    inputs_dict['df'] = df 

    # end read_fitres_file

def print_info(inputs_dict):

    cid_list = inputs_dict['cid_list']
    var_list = inputs_dict['var_list']
    df       = inputs_dict['df']

    pd.set_option("display.max_columns", len(var_list) + 1, "display.width", 1000)
    print(df.loc[cid_list, var_list].__repr__())


    #for cid in cid_list:
    #    row = df.loc[cid, :]
    #    print("# - - - - - - - - - - - - - - -")
    #    for var in var_list:
    #        val = row[var]
    #        print(f" CID={cid}  {var:<12} = {val} ")

    # end print_info

# =============================================
if __name__ == "__main__":

    inputs_dict = {}

    # read command-line arguments
    args  = get_args()
    inputs_dict.update( { 'args' : args } )

    # parse comma-sep lists
    parse_inputs(inputs_dict)

    # read fitres file
    read_fitres_file(inputs_dict)

    # print requested info
    print_info(inputs_dict)

    # end main


