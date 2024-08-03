#!/usr/bin/env python
#
# Created Aug 2 2024 by R.Kessler
#
# make plots from output of
#   snlc_sim.exe NOFILE SIMLIB_FILE <simlib>  SIMLIB_DUMP 1
#
# Beware that "convert" is used to combine png files into a single pdf.
#
# ===============================================

import os, sys, argparse, glob

PLOT_SCRIPT = "plot_table.py"

PREFIX_SIMLIB_DUMP = "SIMLIB_DUMP_SUMMARY_"
SUFFIX_SIMLIB_DUMP = ".TEXT"

# ==============
def get_args():

    parser = argparse.ArgumentParser()

    msg = "name of the SIMLIB DUMP file"
    parser.add_argument("input_file", help=msg,
                        nargs="?", default=None)

    args = parser.parse_args()
    args.prefix = args.input_file.split('.')[0]
    args.n_plot = 0
    
    return args

def next_filenames(args):
    args.n_plot  += 1
    prefix   = f"temp_plot{args.n_plot:02d}_{args.prefix}"
    pdf_file = f"{prefix}.png"
    log_file = f"{prefix}.log"
    return pdf_file, log_file

def get_survey_filters(args):

    inp_file = args.input_file

    j = len(PREFIX_SIMLIB_DUMP)
    survey_plus_filters = inp_file[j:].split('.')[0]
    survey   = survey_plus_filters.split('-')[0]
    filters  = survey_plus_filters.split('-')[1]    

    print(f" SURVEY:   {survey}")
    print(f" FILTERS:  {filters}")
    sys.stdout.flush()
    
    args.survey  = survey
    args.filters = filters

    return

def make_plots(args):

    tfile = args.input_file
    title = f"{args.survey} SIMLIB Properties"
    
    pdf_file, log_file = next_filenames(args)
    varname = "RA vs. DEC"
    print(f" Plot {varname:<12} to {pdf_file} ")
    sys.stdout.flush() 
    cmd = \
        f"{PLOT_SCRIPT} @@TFILE {tfile} @V RA:DEC @@ALPHA 0.2 " + \
        f"@@TITLE '{title}' " + \
        f"@@SAVE {pdf_file} " \
        f" >& {log_file}"

    os.system(cmd)

    # - - - - -

    # set up dictionary of varnames and bounds to plot
    var_dict = {
        'GAPMAX' : '0 30 1' ,
        'GAPAVG' : '0 30 1' ,
        'MJDMIN' : '55000 65000 100 ' ,
        'MJDMAX' : '55000 65000 100 ' ,        
    }
    filter_list = list(args.filters)
    for f in filter_list:
        var = f'ZPT_{f}'
        var_dict[var] = '15 35 0.5'

    for f in filter_list:
        var = f'M5SIG_{f}'
        var_dict[var] = '15 35 0.5'
        
    # - - - -
    for varname, bounds in var_dict.items():        
        pdf_file, log_file = next_filenames(args)
        print(f" Plot {varname:<12} to {pdf_file} ")
        sys.stdout.flush() 
        cmd = \
            f"{PLOT_SCRIPT} @@TFILE {tfile} " + \
            f"@V {varname}  @@BOUNDS {bounds} " + \
            f"@@TITLE '{title}'  @@LEGEND .   @@ALPHA 0.8  " + \
            f"@@OPT GRID NEVT MEAN STDDEV  " + \
            f"@@SAVE {pdf_file} " + \
            f" >& {log_file}"
        
        #print(f"\n xxx cmd = \n{cmd}\n")
        sys.stdout.flush() 
        os.system(cmd)
    
    return

# ===================================================
#   Add main, June 2024
# ========================================
if __name__ == "__main__":


    args = get_args()
    
    print(f" Make diagnostic plots for {args.input_file}")
    sys.stdout.flush()
    
    # get SURVEY and filters from input file name
    get_survey_filters(args)
    
    make_plots(args)

    # combine plots into one file
    pdf_file = f"plot_{args.prefix}.pdf"
    cmd = f"convert temp_plot* {pdf_file}"
    print(f"\n Combine plots into {pdf_file}\n\t {cmd}")
    sys.stdout.flush()
    os.system(cmd)
    
    # === END MAIN ===
