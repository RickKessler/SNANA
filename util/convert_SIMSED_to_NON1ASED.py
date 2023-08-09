#!/usr/bin/env python
#
# Convert to SIMSED.XXX model into NON1ASED.XXX model because the latter can
# be used with TAKE_SPECTRUM and true SED options. SIMSED.XXX model works
# with pre-computed mags on a grid of redshift and filters, and thus the
# underlying SEDs are not available to produce spectra or true SEDs.
#
# Usage:
#  convert_SIMSED_to_NON1ASED.py --simsed_path <path> --model_name <short_name>
#
# ================

import os, sys, argparse

NON1A_LIST_FILE_BASENAME = "NON1A.LIST"
SED_INFO_FILE_BASENAME   = "SED.INFO"

PREFIX_NON1ASED = "NON1ASED"

# =====================================
def get_args():

    parser = argparse.ArgumentParser()
    
    msg = f"name of SIMSED path"
    parser.add_argument("--simsed_path", help=msg, type=str, default=None)

    msg = f"short name of model to write in {NON1A_LIST_FILE_BASENAME} file " \
          f"(e.g., Ia, SNII, etc...)"
    parser.add_argument("--model_name", help=msg, type=str, default=None)

    msg = f"clobber already existing {NON1A_LIST_FILE_BASENAME} file"
    parser.add_argument("--clobber", help=msg, action="store_true")

    args = parser.parse_args()

    args.simsed_path = os.path.expandvars(args.simsed_path)

    args.NON1A_LIST_FILE = f"{args.simsed_path}/{NON1A_LIST_FILE_BASENAME}"
    args.SED_INFO_FILE   = f"{args.simsed_path}/{SED_INFO_FILE_BASENAME}"

    return args

    # end get_args

def error_checks(args):
    if args.simsed_path is None:
        sys.exit(f"\n ERROR: must provide --simsed_path arg")

    if args.model_name is None:
        sys.exit(f"\n ERROR: must provide --model_name arg")

    if os.path.exists(args.NON1A_LIST_FILE) and not args.clobber:
        sys.exit(f"\n ERROR: already existing {NON1A_LIST_FILE_BASENAME} file in\n" \
                 f"  {args.NON1A_LIST_FILE}")

    # end error_checks

def write_NON1A_LIST_FILE(args):

    fp_sed   = open(args.SED_INFO_FILE,   "rt")
    fp_non1a = open(args.NON1A_LIST_FILE, "wt")
    
    print(f" Read  SED file names from {args.SED_INFO_FILE}")
    print(f" Write SED file names  to  {args.NON1A_LIST_FILE}")
    
    # write comments to out file
    argv_string = " ".join(sys.argv)
    fp_non1a.write(f"# Created by program\n")
    fp_non1a.write(f"# {argv_string}\n")
    fp_non1a.write(f"# \n")
    fp_non1a.write(f"# To use TAKE_SPECTRUM keys or write true SED, \n" \
                   f"# use this NON1ASED model instead of original SIMSED model\n")
    fp_non1a.write(f"# \n")

    sed_line_list = fp_sed.readlines()  
    nsed = 0 
    for line in sed_line_list:  
        line  = line.rstrip()  # remove trailing space and linefeed
        words = line.split()
        if len(words) == 0 : continue

        if words[0] == "SED:" :
            sed_file = words[1]
            nsed += 1
            line_out = f"NON1A:  {nsed:3d}  {args.model_name}  {sed_file}"
            fp_non1a.write(f"{line_out}\n")

    fp_sed.close()
    fp_non1a.close()
    
    print(f"\t Finished write {nsed} SED file names.\n")

    # end write_NON1A_LIST_FILE

def create_symlink(args):


    simsed_path = args.simsed_path
    base_simsed = os.path.basename(simsed_path)
    dir_simsed  = os.path.dirname(simsed_path)

    if dir_simsed == "" : dir_simsed = "."

    # construct name of NON1ASED model as follows.
    # If simsed_path is /path/SIMSED.XYZ then
    # base name for NON1ASED model is NON1ASED.XYZ

    suffix = base_simsed.split('.')[1]
    base_non1ased = f"{PREFIX_NON1ASED}.{suffix}"
    full_non1ased = f"{dir_simsed}/{base_non1ased}"

    #print(f" xxx full_non1ased = {full_non1ased}")

    if os.path.exists(full_non1ased):
        print(f" Symbolic link for {base_non1ased} already exists.")

    else:
        cmd = f"cd {simsed_path}/../ ; ln -s {base_simsed} {base_non1ased}"    
        print(f" Create symbolic link for {base_non1ased}")
        os.system(cmd)

    #print(f"\n xxx cmd(symlink) = \n {cmd}\n")

    # end create_symlink

# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":
    
    # read command-line arguments
    args  = get_args()

    error_checks(args)

    write_NON1A_LIST_FILE(args)

    create_symlink(args)

    # end main
