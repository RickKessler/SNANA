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
# Feb 07 2022: fix bug skipping comment lines before VARNAMES (see nrow_skip)
# Jul 16 2022: [Patrick Armstrong] add generic function (mean, min, max, etc...) and command line plotting

import os, sys, argparse, gzip
import plotext as plt
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

    msg = "comma seperated list of functions (mean, min, or max)"
    parser.add_argument("--func", help=msg, type=str, default=None)

    msg = "type of plot. Options include hist (histogram of each variable), and scatter (scatter plot, must specify 2 variables)"
    parser.add_argument("--plot", help=msg, type=str, default=None)

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

    exist_func = args.func is not None
    exist_plot = args.plot is not None
    msgerr_cid     = "If --func or --plot are not defined, must get cid list using --cid or --nrow arg"

    exist_var_list = args.varname is not None
    msgerr_var     = "Must specify varnames using -v arg"

    assert exist_ff,       msgerr_ff
    assert exist_cid_list or exist_func or exist_plot, msgerr_cid
    assert exist_var_list, msgerr_var

    if exist_func:
        # Check func is min, max, or mean
        correct_func = True in [x in ['min', 'max', 'mean'] for x in args.func.split(",")]
        msgerr_func = f"--func must be either min, mean, or max, not {args.func}"
        assert correct_func, msgerr_func 

    var_list   = args.varname.split(",")

    if exist_plot:
        # Check plot is hist or scatter 
        correct_plot = args.plot in ['hist', 'scatter']
        msgerr_plot = f"--plot must be either hist or scatter, not {args.plot}"
        assert correct_plot, msgerr_plot

        if args.plot == 'scatter':
            assert len(var_list) == 2, f"Only 2 variables can be defined when creating a scatter plot, not {len(var_list)}"

    
    id_list = []
    func_list = []
    keyname_id = None

    if args.cid is not None :    
        id_list = args.cid.split(",")
        keyname_id = 'CID'
    elif args.galid is not None :    
        id_list = args.galid.split(",")
        keyname_id = 'GALID'
    if args.func is not None :
        func_list = args.func.split(",")

    nrow       = args.nrow

    info_fitres = {
        'fitres_file' : args.fitres_file,
        'id_list'     : id_list,
        'nrow'        : nrow,
        'var_list'    : var_list,
        'keyname_id'  : keyname_id,
        'func_list'   : func_list,
        'plot'        : args.plot
        
    }
    return info_fitres

    # end parse_inputs

def read_fitres_file(info_fitres):

    # strip inputs
    ff         = info_fitres['fitres_file']  # fitres file or hostlib
    var_list   = info_fitres['var_list']     # list of variables
    keyname_id = info_fitres['keyname_id']
 
    print(f"\n Read {ff}")

    nrow_skip   = 0
    nrow_docana = 0
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
            nrow_skip += 1
            continue

        if line[0] == '#' :
            nrow_skip += 1
            continue

        if wdlist[0] == KEYLIST_DOCANA[0] : isrow_docana = True
        if wdlist[0] == KEYLIST_DOCANA[1] : isrow_docana = False
        if isrow_docana : 
            nrow_skip   += 1
            nrow_docana += 1
            continue

        if wdlist[0] == 'VARNAMES:' : 
            keyname_id = wdlist[1]
            print(f" Found keyname_id = {keyname_id}")
            if nrow_docana > 0:
                print(f"\t (skipped {nrow_docana} DOCANA rows)")

            info_fitres['keyname_id'] = keyname_id
            break

    f.close()
    # - - - - - - 

    var_list_local =  [ keyname_id ] + var_list
    if info_fitres['nrow'] > 0:
        df  = pd.read_csv(ff, comment="#", delim_whitespace=True, 
                          skiprows=nrow_skip,
                          usecols=var_list_local,
                          nrows = info_fitres['nrow'])
    else:
        df  = pd.read_csv(ff, comment="#", delim_whitespace=True, 
                          skiprows=nrow_skip,
                          usecols=var_list_local)


    df[keyname_id] = df[keyname_id].astype(str)
    df             = df.set_index(keyname_id, drop=False)

    # other columns to print as string to avoid int->float roundoff
    keyname_str_list = [ 'SIM_HOSTLIB_GALID', 'HOST_OBJID' ]
    for k in keyname_str_list:
        if k in df:
            df[k] = df[k].astype(str)

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
    func_list  = info_fitres['func_list']
    plot       = info_fitres['plot']

    pd.set_option("display.max_columns", len(df.columns) + 1, 
                  "display.width", 1000)

    # check option to use IDs from first 'nrow' rows
    if nrow > 0 :
        id_rows = df[keyname_id].head(nrow).tolist()
        #print(f" xxx cid_rows = {cid_rows}   ty={type(cid_rows)}")
        id_list += id_rows
        id_list = list(set(id_list))

    # Only print if id_list defined either by cid or nrow
    if len(id_list) > 0:
        df = df.loc[sorted(id_list), var_list]
        print(df.__repr__())
    else:
        df = df.loc[df[keyname_id].to_list(), var_list]

    for func in func_list: 
        if func == 'mean':
            print("\n",df.mean().to_frame('Mean'),"\n")
        if func == 'min':
            print("\n",df.min().to_frame('Min'),"\n")
        if func == 'max':
            print("\n",df.max().to_frame('Max'),"\n")

    if plot == "hist":
        for var in var_list:
            plt.hist(df[var].to_list(), label=var)
        plt.show()

    if plot == "scatter":
        x = df[var_list[0]].to_list()
        y = df[var_list[1]].to_list()
        plt.scatter(x, y)
        plt.xlabel(var_list[0])
        plt.ylabel(var_list[1])
        plt.show()

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


