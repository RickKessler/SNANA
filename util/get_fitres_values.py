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
    parser.add_argument("-c", "--cid", help=msg, type=str, default=None)

    msg = "number of rows to fetch CIDs (no need to know CIDs)"
    parser.add_argument("--nrow", help=msg, type=int, default=0)

    msg = "comma separated list of variable names"
    parser.add_argument("-v", "--varname", help=msg, type=str, default=None)

    args = parser.parse_args()
    return args

    # end get_args

def parse_inputs(args):

    # parse comma-sep list of cid(s) and varname(s)

    exist_ff   = os.path.exists(args.file)
    msgerr_ff  = (f"fitres file {args.file} does not exist.")

    exist_cid_list = (args.cid is not None) or (args.nrow > 0 )
    msgerr_cid     = "Must cid list using --cid or --nrow arg"

    exist_var_list = args.varname is not None
    msgerr_var     = "Must specify varnames using -v arg"

    assert exist_ff,       msgerr_ff
    assert exist_cid_list, msgerr_cid
    assert exist_var_list, msgerr_var

    cid_list   = args.cid.split(",")
    nrow       = args.nrow
    var_list   = args.varname.split(",")

    return args.file, cid_list, nrow, var_list

    # end parse_inputs

def read_fitres_file(ff, var_list):

    # inputs:
    #   ff       = name of fitres file
    #   var_list = list of variables to parse

    var_list_local = ['CID'] + var_list
    df  = pd.read_csv(ff, comment="#", delim_whitespace=True, 
                      usecols=var_list_local)
    df["CID"] = df["CID"].astype(str)
    df        = df.set_index("CID", drop=False)
    return df

    # end read_fitres_file

def print_info(df, cid_list, nrow):

    # Inputs:
    #   df = data frame for input fitres file
    #   cid_list = comma-sep list of cids
    #   nrow     = number of rows to fetch CID

    # check option to use CIDs from first 'nrow' rows
    if nrow > 0 :
        cid_rows = df['CID'].iloc[:nrow].values
        print(f" xxx cid_rows = {cid_rows} ")
        #cid_list += cid_rows

    print(" xxx - - - - - - - - \n ")

    pd.set_option("display.max_columns", len(df.columns) + 1, 
                  "display.width", 1000)
    print(df.loc[cid_list, var_list].__repr__())
    # end print_info

# =============================================
if __name__ == "__main__":

    # read command-line arguments
    args  = get_args()

    # parse comma-sep lists
    ff, cid_list, nrow, var_list = parse_inputs(args)

    # read fitres file
    df = read_fitres_file(ff, var_list)

    # print requested info
    print_info(df,cid_list, nrow)

    # end main


