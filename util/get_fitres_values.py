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
# Jul 16 2022: 
#   [Patrick Armstrong] add generic function (mean, min, max, etc...) 
#    and command line plotting
# Jul 19 2022 RK 
#     --filt replaced with --sel to avoid confusion with filters.
#     Put plotext import under try to avoid abort if not installed.
#
# Aug 21 2022 RK - fix bug counting lines before VARNAMES key
# Jul 31 2023 RK - new option --reformat
# Feb 02 2024 RK 
#   + improve stdout appearance
#   + new --outfile option
#
import os, sys, argparse, gzip, math
import numpy as np
import pandas as pd

try:  
    import plotext as plt
except ImportError:
    pass

KEYLIST_DOCANA = [ 'DOCUMENTATION:', 'DOCUMENTATION_END:' ]

STRING_FORMAT_EAZY  = "eazy"
ZP_nJy              = 31.4
ZP_SNANA            = 27.5
MAGERR_DEFAULT      = 0.05

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

    msg = "comma seperated list of selection functions: " \
          "combined via 'or'. Can use <, >, or =="
    parser.add_argument("--sel", help=msg, type=str, default=None)

    msg = "comma separated list of variable names to  print values"
    parser.add_argument("-v", "--varname", help=msg, type=str, default=None)

    msg = "comma seperated list of functions (mean, min, or max)"
    parser.add_argument("--func", help=msg, type=str, default=None)

    msg = "format to re-write fitres/hostlib file " \
          "(e.g, eazy*10 -> split into 10 files)"
    parser.add_argument("--reformat", help=msg, type=str, default=None)

    msg = "type of plot. Options include hist (histogram of each variable), and scatter (scatter plot, must specify 2 variables)"
    parser.add_argument("--plot", help=msg, type=str, default=None)

    msg = "name of output table in fitres format"
    parser.add_argument("-o", "--outfile", help=msg, type=str, default=None)

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
                     (args.sel   is not None) or \
                     (args.nrow > 0)

    exist_func = args.func is not None
    exist_plot = args.plot is not None
    msgerr_cid     = "If --func or --plot are not defined, must get cid list using --cid, --nrow, or --sel arg"

    exist_var_list = args.varname is not None or args.reformat is not None
    msgerr_var     = "Must specify varnames using -v arg"

    
    assert exist_ff,       msgerr_ff
    assert exist_var_list, msgerr_var

    exist_tmp = exist_cid_list or exist_func or exist_plot or args.reformat
    assert exist_tmp, msgerr_cid

    if exist_func:
        # Check func is min, max, or mean
        correct_func = True in [x in ['min', 'max', 'mean'] for x in args.func.split(",")]
        msgerr_func = f"--func must be either min, mean, or max, not {args.func}"
        assert correct_func, msgerr_func 

    if args.reformat is None:
        var_list   = args.varname.split(",")
    else:
        var_list = []

    if exist_plot:
        # Check plot is hist or scatter 
        correct_plot = args.plot in ['hist', 'scatter']
        msgerr_plot = f"--plot must be either hist or scatter, not {args.plot}"
        assert correct_plot, msgerr_plot

        if args.plot == 'scatter':
            assert len(var_list) == 2, f"Only 2 variables can be defined when creating a scatter plot, not {len(var_list)}"

    
    id_list    = []
    func_list  = []
    sel_list   = []
    keyname_id = None

    if args.cid is not None :    
        id_list = args.cid.split(",")
        keyname_id = 'CID'
    elif args.galid is not None :    
        id_list = args.galid.split(",")
        keyname_id = 'GALID'
    if args.func is not None :
        func_list = args.func.split(",")
    if args.sel is not None :
        for sel in args.sel.split(","):
            if "<" in sel:
                v, n = sel.split("<")
                sel_list.append([v, "<", float(n)])
            elif ">" in sel:
                v, n = sel.split(">")
                sel_list.append([v, ">", float(n)])
            elif "==" in sel:
                v, n = sel.split("==")
                sel_list.append([v, "==", float(n)])
            else:
                raise ValueError(f"Valid select functions are <, >, and ==. None of these were not found in {sel}")

    nrow       = args.nrow

    info_fitres = {
        'fitres_file' : args.fitres_file,   # input fitres file
        'outfile'     : args.outfile,       # optional output
        'id_list'     : id_list,
        'nrow'        : nrow,
        'var_list'    : var_list,
        'keyname_id'  : keyname_id,
        'func_list'   : func_list,
        'plot'        : args.plot,
        'sel_list'    : sel_list
    }
    return info_fitres

    # end parse_inputs

def read_fitres_file(info_fitres, reformat_option):

    # strip inputs
    ff         = info_fitres['fitres_file']  # fitres file or hostlib
    var_list   = info_fitres['var_list']     # list of variables
    keyname_id = info_fitres['keyname_id']
 
    print(f"\n# Read {ff}")

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
        if isrow_docana : 
            nrow_skip   += 1
            nrow_docana += 1

        if wdlist[0] == KEYLIST_DOCANA[1] : isrow_docana = False

        if wdlist[0] == 'VARNAMES:' : 
            keyname_id = wdlist[1]
            print(f"# Found keyname_id = {keyname_id} " \
                  f"after skipping {nrow_skip} rows")
            if nrow_docana > 0:
                print(f"#\t (skipped {nrow_docana} DOCANA rows)")

            info_fitres['keyname_id'] = keyname_id
            VARLIST_ALL = wdlist[2:]
            break

    f.close()
    # - - - - - - 

    # check reformat option(s) that require specific HOSTLIB variables
    if reformat_option:
        if STRING_FORMAT_EAZY in reformat_option:
            var_list = []
            for var in VARLIST_ALL:
                if '_obs' in var:
                    var_list.append(var)
            info_fitres['var_list'] = var_list

    # - - - -
    var_list_local =  [ keyname_id ] + var_list

    # - - - - 
    if info_fitres['nrow'] > 0:
        df  = pd.read_csv(ff, comment="#", delim_whitespace=True, 
                          skiprows=nrow_skip, dtype=str,
                          usecols=var_list_local,                          
                          nrows = info_fitres['nrow'])

    else:
        df  = pd.read_csv(ff, comment="#", delim_whitespace=True, 
                          skiprows=nrow_skip, dtype=str,
                          usecols=var_list_local)

    # xxx mark df[keyname_id] = df[keyname_id].astype(str)
    df             = df.set_index(keyname_id, drop=False)

    # xxx mark delete Feb 2024 xxxxxx
    # set a few other columns to print as string to avoid int->float roundoff
    #keyname_str_list = [ 'SIM_HOSTLIB_GALID', 'HOST_OBJID' ]
    #for k in keyname_str_list:
    #    if k in df:
    #        df[k] = df[k].astype(str)
    # xxxxxxxxx

    # load data frame
    info_fitres['df'] = df
    return

    # end read_fitres_file

def print_info(info_fitres):

    # Inputs:
    #   df       = data frame for input fitres file
    #   id_list  = comma-sep list of cids or galids
    #   nrow     = number of rows to fetch CID

    df           = info_fitres['df']
    id_list      = info_fitres['id_list']
    nrow         = info_fitres['nrow']
    var_list     = info_fitres['var_list']
    keyname_id   = info_fitres['keyname_id']
    func_list    = info_fitres['func_list']
    plot         = info_fitres['plot']
    sel_list     = info_fitres['sel_list']
    outfile      = info_fitres['outfile']

    pd.set_option("display.max_columns", len(df.columns) + 1, 
                  "display.width", 1000)

    # check option to use IDs from first 'nrow' rows
    if nrow > 0 :
        id_rows = df[keyname_id].head(nrow).tolist()
        #print(f" xxx cid_rows = {cid_rows}   ty={type(cid_rows)}")
        id_list += id_rows

    if len(sel_list) > 0:
        for sel in sel_list:
            if sel[1] == "<":
                id_rows = df[df[sel[0]] < sel[2]][keyname_id].tolist()
            elif sel[1] == ">":
                id_rows = df[df[sel[0]] > sel[2]][keyname_id].tolist()
            elif sel[1] == "==":
                id_rows = df[df[sel[0]] == sel[2]][keyname_id].tolist()
            id_list += id_rows

    id_list = list(set(id_list))
    # Only print if id_list defined either by cid or nrow
    if len(id_list) > 0:

        # xxx mark df = df.loc[sorted(id_list), var_list]
        # xxx mark delete print(df.__repr__())

        # add back keynam_id [e.g., CID or GALID] and then print
        # without index ... this avoids index name (CID) printed
        # to a separate row compared to other varnames
        df = df.loc[sorted(id_list), [keyname_id] + var_list] # RK Feb 2024
        print(df.to_string(index=False))                      # RK Feb 2024

        if outfile is not None:
            write_outfile_fitres(info_fitres)

    else:
        df = df.loc[df[keyname_id].to_list(), var_list]

    for func in func_list: 
        if func == 'mean':
            print("\n",df.mean().to_frame('Mean'),"\n")
        if func == 'min':
            print("\n",df.min().to_frame('Min'),"\n")
        if func == 'max':
            print("\n",df.max().to_frame('Max'),"\n")

    # warning: this plotter may be obsolete when plot_fitres.py is installed 
    if plot == "hist":
        for var in var_list:
            data = df[var]
            q1 = data.quantile(0.25)
            q3 = data.quantile(0.75)
            iqr = q3 - q1
            bin_width = (2 * iqr) / (len(data) ** (1 / 3))
            bin_count = int(np.ceil((data.max() - data.min()) / bin_width))
            plt.hist(df[var].to_list(), bins=bin_count, label=var)
        plt.show()

    if plot == "scatter":
        x = df[var_list[0]].to_list()
        y = df[var_list[1]].to_list()
        plt.scatter(x, y)
        plt.xlabel(var_list[0])
        plt.ylabel(var_list[1])
        plt.show()

    # end print_info

def write_outfile_fitres(info_fitres):

    df         = info_fitres['df']
    outfile    = info_fitres['outfile']
    var_list   = info_fitres['var_list']
    keyname_id = info_fitres['keyname_id']

    print(f"\n# Write selected fitres value to table in {outfile}")

    varnames_list   = [ keyname_id ] + var_list
    varnames_string = ' '.join(varnames_list)
    user_command    = ' '.join(sys.argv)
    nvar       = len(varnames_list)

    with open(outfile,"wt") as f:   
        f.write(f"# This table created by command \n# {user_command}\n#\n")
        f.write(f"VARNAMES: {varnames_string}\n")

        #f.write(df.to_string(index=False))
        for index, row in df.iterrows():
            row_values = ""
            for i in range(0,nvar):
                row_values += f"{row[i]} "  # gotta be a better way 

            f.write(f"SN: {row_values}\n")

    return
    # end write_outfile_fitres

def reformat(reformat, info_fitres):

    # reformat = eazy -> write entire hostlib in eazy format.
    # reformat = eazy*10 -> split output into 10 files

    ff         = info_fitres['fitres_file']  # fitres file or hostlib
    df         = info_fitres['df']
    id_list    = info_fitres['id_list']
    nrow       = info_fitres['nrow']
    var_list   = info_fitres['var_list']
    keyname_id = info_fitres['keyname_id']
    
    nrow_df = len(df)

    reformat_split = reformat.split('*')
    format_string  = reformat_split[0]
    if len(reformat_split) == 1 :
        nfile_split = 1
    else:
        nfile_split = int(reformat_split[1])
    
    print(f"\n Reformat {nrow_df} rows into " \
          f"{nfile_split} {format_string} files." )

    # - - - -

    info_fitres['nfile_split'] = nfile_split
    info_fitres['zp']          = ZP_nJy
    #print(f"\n xxx df = \n{df}\n")

    basename = os.path.basename(ff)

    # make sure basename has .gz extension
    if '.gz' not in basename:
        basename += '.gz'

    for i in range(0,nfile_split):
        outfile_reformat = f"{format_string}{i:02d}_{basename}"
        info_fitres['isplit']  = i
        info_fitres['outfile_reformat'] = outfile_reformat
        print(f"  Create {outfile} ")
        if format_string == STRING_FORMAT_EAZY :
            if i==0 :
                add_fluxcal_lists(info_fitres)
            reformat_eazy(info_fitres)
        
    return
    # end reformat

def reformat_eazy(info_fitres):

    # write gal id and fluxcal in format usable by eazy code.

    outfile_reformat       = info_fitres['outfile_reformat']
    isplit                 = info_fitres['isplit']
    nfile_split            = info_fitres['nfile_split']
    varname_fluxcal_list   = info_fitres['varname_fluxcal_list'] 
    varname_fluxcal_string = info_fitres['varname_fluxcal_string'] 
    
    nrow = len(info_fitres['GALID'])

    with gzip.open(outfile_reformat,"wt") as f:
        f.write(f"# id {varname_fluxcal_string}\n")
        for i in range(0,nrow):
            if (i % nfile_split) != isplit : continue
            GALID = info_fitres['GALID'][i]
            line  = str(GALID) + ' ' 
            for var in varname_fluxcal_list:
                val = info_fitres[var][i]
                line += f"{val:.3e} "
            f.write(f"{line}\n")

    return
    # end reformat_eazy

def add_fluxcal_lists(info_fitres):
    
    df       = info_fitres['df']
    zp       = info_fitres['zp']
    var_list = info_fitres['var_list']
    get_fluxcal_vectorized      = np.vectorize(get_fluxcal)
    get_fluxerr_frac_vectorized = np.vectorize(get_fluxerr_frac)

    info_fitres['GALID'] = df['GALID'].to_numpy()
    info_fitres['varname_fluxcal_list']   = []
    info_fitres['varname_fluxcal_string'] = ''

    for var_mag in var_list:
        if 'err' in var_mag: continue 

        # check if mag_err column exists
        var_magerr   = var_mag + '_err'
        exist_magerr = var_magerr in var_list

        band        = var_mag[0]
        var_flux    = 'fluxcal_'        + band # e.g., g_obs -> fluxcal_g
        var_fluxerr = 'fluxcalerr_'     + band
        info_fitres[var_flux]    = get_fluxcal_vectorized(zp,df[var_mag])

        if exist_magerr:
            info_fitres['fluxerr_frac'] = \
                get_fluxerr_frac_vectorized(df[var_magerr])
            info_fitres[var_fluxerr] = \
                    np.multiply(info_fitres[var_flux],
                                info_fitres['fluxerr_frac'] )
        else:
            # if there are no errors provided, use a default
            info_fitres[var_fluxerr] = info_fitres[var_flux] * MAGERR_DEFAULT

        info_fitres['varname_fluxcal_list'].append(var_flux)
        info_fitres['varname_fluxcal_list'].append(var_fluxerr)

        info_fitres['varname_fluxcal_string'] += f"{var_flux} "
        info_fitres['varname_fluxcal_string'] += f"{var_fluxerr} "

        print(f"\t compute {var_flux} list")

    return
    # end add_fluxcal_columns

def get_fluxcal(zp,mag):
    arg = -0.4*(mag-zp)
    fluxcal = math.pow(10.0,arg)
    return fluxcal
    # end

def get_fluxerr_frac(magerr):
    # return fluxerr/flux corresponding to input magerr

    if magerr > 10.0 : magerr = 10.0

    arg  = 0.4*magerr
    frac = math.pow(10.0,arg) - 1.0
    return frac
    # end

# =============================================
if __name__ == "__main__":

    # read command-line arguments
    args  = get_args()

    info_fitres = parse_inputs(args)

    # read fitres file, tack on info_fitres['df']
    read_fitres_file(info_fitres, args.reformat)

    # print requested info
    if args.reformat is not None:
        reformat(args.reformat,info_fitres)
    else:
        print_info(info_fitres)

    # end main


