#!/usr/bin/env python
#
# Created Oct 2020 [replace old perl script get_fitresValue.pl]
#
# Extract table value(s) for cid(s) and print to screen.
# Intended for visual debugging.
#
# For "-f bla.fitres", script checks for both bla.fitres and bla.fitres.gz.
# for "-f bla.fitres.gz", script checks only for this file name.
#
# Feb 24 2021: print table with sorted(cid_list) to avoid random ordering.
# Sep 23 2021: add -g (--galid) option to work for HOSTLIB, and 
#              count DOCANA rows to skip
#
# Nov 22 2021: fix to work with gzip files
#
import os, sys, argparse, gzip
import pandas as pd

KEYLIST_DOCANA = [ 'DOCUMENTATION:', 'DOCUMENTATION_END:' ]

# =====================================
def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "name of fitres file (automatically checks for .gz extension)"
    parser.add_argument("-f", "--fitres_file", help=msg, type=str, default="")

    msg = "comma separated list of CIDs for FITRES file"
    parser.add_argument("-c", "--cid", help=msg, type=str, default=None)

    msg = "comma sep list of GALIDs for HOSTLIB"
    parser.add_argument("-g", "--galid", help=msg, type=str, default=None)

    msg = "number of rows to fetch CIDs (no need to know CIDs)"
    parser.add_argument("--nrow", help=msg, type=int, default=0)

    msg = "comma separated list of variable names"
    parser.add_argument("-v", "--varname", help=msg, type=str, default=None)

    if len(sys.argv) == 1:  parser.print_help(); sys.exit()

    args = parser.parse_args()

    return args

    # end get_args

def parse_inputs(args):

    # parse comma-sep list of cid(s) and varname(s)
    fitres_file = args.fitres_file

    # - - - - - -
    # check fitres_file and fitres_file.gz
    msgerr_ff  = (f"fitres file {args.fitres_file} does not exist.")
    ff_list    = [ fitres_file ]
    if ".gz" not in fitres_file:
        ff_list.append(f"{fitres_file}.gz")
    
    n_ff = 0
    exist_ff = False
    for ff in ff_list:
        if os.path.exists(ff):
            n_ff += 1
            exist_ff = True
            args.fitres_file = ff

    # - - - - - - -
    exist_cid_list = (args.cid   is not None) or \
                     (args.galid is not None) or \
                     (args.nrow > 0)

    msgerr_cid     = "Must cid list using --cid or --nrow arg"
    exist_var_list = args.varname is not None
    msgerr_var     = "Must specify varnames using -v arg"

    assert exist_ff,       msgerr_ff
    assert exist_cid_list, msgerr_cid
    assert exist_var_list, msgerr_var
    
    
    id_list = []
    keyname_id = None

    if args.cid is not None :    
        id_list = args.cid.split(",")
        keyname_id = 'CID'
    if args.galid is not None :    
        id_list = args.galid.split(",")
        keyname_id = 'GALID'

    nrow       = args.nrow
    var_list   = args.varname.split(",")

    info_fitres = {
        'fitres_file' : args.fitres_file,
        'id_list'     : id_list,
        'nrow'        : nrow,
        'var_list'    : var_list,
        'keyname_id'  : keyname_id
        
    }
    return info_fitres

    # end parse_inputs

def read_fitres_file(info_fitres):

    # strip inputs
    ff         = info_fitres['fitres_file']  # fitres file or hostlib
    var_list   = info_fitres['var_list']     # list of variables
    keyname_id = info_fitres['keyname_id']
 
    print(f"\n Read {ff}")

    nskip_row = 0
    isrow_docana = False

    if ".gz" in ff:
        f = gzip.open(ff,"rt", encoding='utf-8')
    else:
        f = open(ff,"rt")

    # - - - - - 
    # check keyname of id; e.g., CID, ROW, GALID, etc ...
    # read first VARNAMES element and number of rows to
    # skip for DOCANAN keys

    for line in f:       
        wdlist = line.split()
        if len(wdlist) < 1 : 
            nskip_row += 1
            continue

        if wdlist[0] == KEYLIST_DOCANA[0] : isrow_docana = True
        if wdlist[0] == KEYLIST_DOCANA[1] : isrow_docana = False
        if isrow_docana : 
            nskip_row += 1
            continue

        if wdlist[0] == 'VARNAMES:' : 
            keyname_id = wdlist[1]
            print(f" Found keyname_id = {keyname_id}")
            info_fitres['keyname_id'] = keyname_id
            break

    f.close()
    # - - - - - - 

    var_list_local =  [ keyname_id ] + var_list
    df  = pd.read_csv(ff, comment="#", delim_whitespace=True, 
                      skiprows=nskip_row,
                      usecols=var_list_local)

    df[keyname_id] = df[keyname_id].astype(str)
    df             = df.set_index(keyname_id, drop=False)

    # load data frame
    info_fitres['df'] = df
    return

    # end read_fitres_file

def print_info(info_fitres):

    # Inputs:
    #   df       = data frame for input fitres file
    #   id_list  = comma-sep list of cids or galids
    #   nrow     = number of rows to fetch CID

    df         = info_fitres['df']
    id_list    = info_fitres['id_list']
    nrow       = info_fitres['nrow']
    var_list   = info_fitres['var_list']
    keyname_id = info_fitres['keyname_id']

    # check option to use IDs from first 'nrow' rows
    if nrow > 0 :
        id_rows = df[keyname_id].head(nrow).tolist()
        #print(f" xxx cid_rows = {cid_rows}   ty={type(cid_rows)}")
        id_list += id_rows
        id_list = list(set(id_list))

    pd.set_option("display.max_columns", len(df.columns) + 1, 
                  "display.width", 1000)
    print(df.loc[sorted(id_list), var_list].__repr__())
    # end print_info

# =============================================
if __name__ == "__main__":

    # read command-line arguments
    args  = get_args()

    info_fitres = parse_inputs(args)
    # parse comma-sep lists
    ## xxx ff, id_list, nrow, var_list = parse_inputs(args)

    # read fitres file, tack on info_fitres['df']
    read_fitres_file(info_fitres)

    # print requested info
    print_info(info_fitres)

    # end main


