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
    pdf_file = f"temp_plot{args.n_plot}_{args.prefix}.png"
    log_file = f"temp_plot{args.n_plot}_{args.prefix}.log"
    return pdf_file, log_file

def make_plots(args):

    pdf_file, log_file = next_filenames(args)
    print(f" Plot RA vs. DEC to {pdf_file} ")
    cmd = \
        f"{PLOT_SCRIPT} @V DEC:RA @@ALPHA 0.2 " \
        f"@@SAVE {pdf_file} " \
        f" >& {log_file}"
    os.system(cmd)


    pdf_file, log_file = next_filenames(args)
    print(f" Plot MJD to {pdf_file} ")
    cmd = \
        f"{PLOT_SCRIPT} @V MJD @@ALPHA 0.8 " \
        f"@@SAVE {pdf_file} " \
        f" >& {log_file}"
    os.system(cmd)


    pdf_file, log_file = next_filenames(args)
    print(f" Plot NSEASON to {pdf_file} ")
    cmd = \
        f"{PLOT_SCRIPT} @V NSEASON @@ALPHA 0.8 " \
        f"@@SAVE {pdf_file} " \
        f" >& {log_file}"
    os.system(cmd)
    
    
    return

# ===================================================
#   Add main, June 2024
# ========================================
if __name__ == "__main__":


    args = get_args()
    
    print(f" Make diagnostic plots for {args.input_file}")

    make_plots(args)

    # combine plots into one file
    pdf_file = f"plot_{args.prefix}.pdf"
    cmd = "convert temp_plot* {pdf_file}"
    print(f"\n Combine plots into {pdf_file}")
    os.system(cmd)
    
    # === END MAIN ===
