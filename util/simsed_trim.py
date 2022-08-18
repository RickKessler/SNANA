#!/usr/bin/env python
#
# Created Auguest 2022 by R.Kessler
# Utility to trim SIMSED time-series model with cut on 
# DAY and/or LAM range.
# First usage is to remove DAY<0 for KN-K17 model where 
# SED fluxes are all zero.
#
# Example usage:
#  simsed_trim.py \
#     --path_model $SNDATA_ROOT/models/SIMSED/SIMSED.KNovae_BarnesKasen2013 \
#     --cutwin_day 4 50 \
#     --cutwin_lam 1250 25000
#
# will output SIMSED.KNovae_BarnesKasen2013_trim subdirectory
# and each SED time series will be truncated with DAY>4 and LAM>1250 A.
#
# ========================

import os, sys, argparse, glob, yaml, shutil, gzip

KEY_SED = "SED:"

# =========================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "path of SIMSED model "
    parser.add_argument("--path_model", help=msg, type=str, default=None)

    msg = "DAY cut-window to keep"
    parser.add_argument("--cutwin_day", help=msg, nargs=2, type=float)

    msg = "LAM cut-window to keep"
    parser.add_argument("--cutwin_lam", help=msg, nargs=2, type=float)


    args = parser.parse_args()

    args.path_model = os.path.expandvars(args.path_model)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return args

    # end get_args

def read_sed_list(model):
    # read SED.INFO file and return list of SED files

    sed_info_file = f"{model}/SED.INFO"

    if not os.path.exists(sed_info_file) :
        sys.exit(f"\n ERROR: Cannot find {sed_info_file}")

    sed_list = []
    with  open(sed_info_file,"rt") as f:
        for line in f:
            line  = line.rstrip() 
            if len(line) == 0 : continue
            wdlist = line.split()
            if wdlist[0] == KEY_SED :
                sed_list.append(wdlist[1])
         
    return sed_list
    # end read_sed_list

# ==========================================
def create_outdir(model):
    basename = os.path.basename(model)
    outdir   = f"{basename}_trim"

    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    print(f" Create {outdir}")
    os.mkdir(outdir)

    return outdir
    # end create_outdir

def trim_sed(args,outdir,sed):
    path_model = args.path_model

    sed_file_orig   = f"{path_model}/{sed}"
    
    if not os.path.exists(sed_file_orig):
        sed_file_orig += '.gz'

    sed_file_out = f"{outdir}/{sed}.gz"
    sout = gzip.open(sed_file_out,"wt")

    with gzip.open(sed_file_orig,"rt") as sin:
        for line_orig in sin:            
            line_new = line_orig.rstrip()
            if line_new[0] == '#' :
                keep = True
            else:
                wdlist   = line_new.split()
                day      = float(wdlist[0])
                lam      = float(wdlist[1])
                keep = apply_cutwin(args,day,lam)            
            if keep: 
                sout.write(f"{line_new}\n")

    sout.close()

    print(f"   Trim  {sed}")
    return
    # end trim_sed


def apply_cutwin(args,day,lam):

    keep = True

    if args.cutwin_day is not None:
        if day < args.cutwin_day[0] : keep = False
        if day > args.cutwin_day[1] : keep = False

    if args.cutwin_lam is not None:
        if lam < args.cutwin_lam[0] : keep = False
        if lam > args.cutwin_lam[1] : keep = False

    return keep
    # end apply_cutwin

# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    args    = get_args()
    sed_list = read_sed_list(args.path_model)
    outdir   = create_outdir(args.path_model)

    for sed in sed_list:
        trim_sed(args,outdir,sed)

    # end main
