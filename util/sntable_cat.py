#!/usr/bin/env python
#
# Created May 2020 by R.Kessler
#
# Read input list of text-formatted SN table files (i.e., FITRES files)
# and catenate them after accounting for possible column differences
# among input files. See manual Sec 10.1 for explanations.
#
# This python script is just a wrapper that calls SALT2mu.exe
# with 'cat_only' argument that turns SALT2mu.exe into a utility 
# instead of a fitting program.
#
# See input arguments with
#    sntable_cat.py -h
#
# Jun 20 2025: minor refactor of inpfile_list to accept either comma-sep list as before,
#              or accept space-sep list or accept wildcard. Replace OptionParser with argparse.
#              Only changes are in parse_args().
# ------------------------------------------

import os, sys, argparse
#from optparse import OptionParser

# =============================
def parse_args():

    parser=argparse.ArgumentParser()

    msg = 'list of text-formatted SNANA tables (comma sep or space sep)'
    parser.add_argument('-i', help=msg, nargs="+", type=str,
                        dest='inpfile_list', default=None)

    msg = 'name of output/catenated sntable (text format)'
    parser.add_argument('-o',help=msg, type=str,
                        dest='outfile_cat',default=None)

    msg = 'list of column names to append if missing (wildcard allowed, e.g., PROB*)'
    parser.add_argument('-a',help=msg, type=str,
                        dest='append_varname_missing', default='PROB*,zPRIOR*')

    msg = 'optional integer prescale'
    parser.add_argument('-p',help=msg, type=int,
                        dest='prescale', default=1)

    msg = 'DEBUG ONLY: jobname'
    parser.add_argument('-j', help=msg, type=str,
                        dest='jobname', default='SALT2mu.exe')

    INPUTS  = parser.parse_args()

# abort on missing arguments
    if not INPUTS.inpfile_list:
        parser.error('Missing required arg: -i <inpfile_list>')

    # 6/20/2025: convert list to comma-sep string for C code input
    if len(INPUTS.inpfile_list) == 1:
        INPUTS.inpfile_list = INPUTS.inpfile_list[0] # convert list to string
    else:
        INPUTS.inpfile_list = ','.join(INPUTS.inpfile_list) 

    if not INPUTS.outfile_cat:  
        parser.error('Missing required arg: -o <outfile_cat>')

    return INPUTS


# ===================================
def exec_salt2mu(INPUTS):

    inpfile_list   = INPUTS.inpfile_list
    outfile_cat    = INPUTS.outfile_cat
    append_varname = INPUTS.append_varname_missing
    arg_cat        = "cat_only"
    jobname        = INPUTS.jobname

    if INPUTS.prescale > 1: 
        arg_cat += f"/{INPUTS.prescale}"

    arg_list = f"{arg_cat} " \
               f"datafile={inpfile_list} "   \
               f"catfile_out={outfile_cat} " \
               f"append_varname_missing={append_varname} "

#                % (arg_cat, inpfile_list, outfile_cat, append_varname) )

    command = f"{jobname} {arg_list}"
    print(f" {sys.argv[0]}: executes os.system command = \n  '{command}' \n")
    sys.stdout.flush()

    os.system(command)

    print(f"\n Done.")

# ===================================
# ============== MAIN ===============
# ===================================

if __name__ == "__main__":

    INPUTS = parse_args()

# execute SALT2mu.exe program to do the catenate
    exec_salt2mu(INPUTS)


