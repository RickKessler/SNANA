#!/usr/bin/env python
#
# Created March 2021 by R.Kessler
#   [re-write of perl script sntable_dump.pl]
#
# Utility to examine and dump analysis info from any SNTABLE
# in hbook (ntuple) or root (tree) format.
# The basic functions here are
#  - print list of all available variables in SNTABLE
#  - dump subset of variables for each event into text file.
#  - extract flux outliers w.r.t. sim truth or LC fit curve.
#
# This utility is aimed for those who prefer NOT to use 
# hbook or root, and prefer to extract analysis variables
# into text files.
#
# This python script is just a shell that calls underyling C program
# sntable_dump.exe .
#
# Apr 7 2021
#   + add flush statements
#   + exit if Cprogram aborts
#   + return error code for submit_batch_jobs
#
# Jul 14 2023
#   + call print_elapse_time for VERBOSE mode
#   + for combine_fitres.exe call in append mode, add 't' argument to 
#     write only text output and avoid writing useless ROOT file.
# ===================

import os, sys, argparse, yaml, datetime 
from   argparse import Namespace

VARNAME_CCID = "CCID"
FORMAT_ROOT  = "ROOT"
FORMAT_HBOOK = "HBOOK"

# define maps to correct user mistakes on table names because
# HBOOK takes table number while root takes tree name.
# Can specify SNANA or FITRES table name for any format,
# and the HBOOK map internally translates to integer ID
TABLE_NAME_MAP_HBOOK = { 
    'SNANA'   : '7100',    
    'FITRES'  : '7788',
    'OUTLIER' : '7800'
}
TABLE_NAME_MAP_ROOT  = { 
    '7100'  : 'SNANA',   
    '7788'  : 'FITRES',
    '7800'  : 'OUTLIER'
}

OUTFILE_PREFIX = "sntable_dump"
OUTFILE_SUFFIX = "fitres"

Cprogram_dump    = "sntable_dump.exe"
Cprogram_combine = "combine_fitres.exe"

# ===========
def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "name of table file to examine"
    parser.add_argument("input_table_file", help=msg, nargs="?", default=None)

    msg = "name of table (e.g., SNANA, FITRES, SPEC ...)"
    parser.add_argument("table_name", help=msg, nargs="?", default=None)

    msg = "List possible table names"
    parser.add_argument("-l", "--list_tables", help=msg, action="store_true")

    msg = "comma-sep list of variables to extract"
    parser.add_argument("-v", "--varlist", help=msg, nargs='?', 
                        type=str, default=None)

    msg = "Extract data-fit outliers with Nsig_min < |PULL| < Nsig_max"
    parser.add_argument("--outlier_fit", metavar=('Nsig_min', 'Nsig_max'),
                        help=msg, nargs=2, type=float, default=None)

    msg = "Extract data-sim outliers with Nsig_min < |PULL| < Nsig_max"
    parser.add_argument("--outlier_sim", metavar=('Nsig_min', 'Nsig_max'),
                        help=msg, nargs=2, type=float, default=None)

    msg = "append varlist to this already existing fitres text file"
    parser.add_argument("-a", "--append_file", help=msg, nargs='?', 
                        type=str, default=None)

    msg = "output file name for text formatted table"
    parser.add_argument("-o", "--outfile", help=msg, nargs='?', 
                        type=str, default=None)

    msg = "Print Number of events in table"
    parser.add_argument("-N", "--NEVT", help=msg, action="store_true")

    msg = "dump each observation, for FLUXERR map"
    parser.add_argument("-O", "--OBS", help=msg, action="store_true")

    msg = "Verbose mode"
    parser.add_argument("-V", "--VERBOSE", help=msg, action="store_true")

    args = parser.parse_args()

    return args
    # end get_args

def list_tables():

    print(" Possible table names (2nd arg for this script)")
    print(" from &SNLCINP input SNTABLE_LIST, and used by")
    print(" snana.exe, snlc_fit.exe, psnid.exe programs:")
    for table in TABLE_NAME_MAP_HBOOK:
        print(f"\t{table}")

    sys.exit("\nBye Bye.\n")
    # end list_tables

def insert_ccid_varlist(args):

    # append CCID to varlist, and remove commas

    # bail if CCID is already in varlist
    varlist = args.varlist
    if varlist is None : return None

    if VARNAME_CCID in varlist: return varlist

    varlist_new = f"{VARNAME_CCID},{varlist}"

    varlist_new = varlist_new.replace(',',' ')

    # replace commas with spaces for C program
    return varlist_new

    # end insert_ccid_varlist

def get_file_format(input_table_file):

    Format = None
    
    for suffix in [ 'ROOT', 'root' ]:
        if suffix in input_table_file: Format = FORMAT_ROOT

    for suffix in [ 'HBOOK', 'hbook' ]:
        if suffix in input_table_file: Format = FORMAT_HBOOK

    if Format is None :
        sys.exit(f"\n ERROR: unknown format for " \
                 f"input_table_file = {input_table_file}\n")

    return Format
    # end get_file_format

def check_table_name(table_name_orig,Format):

    table_name_new  = table_name_orig

    if Format == FORMAT_HBOOK:
        if table_name_orig in TABLE_NAME_MAP_HBOOK:
            table_name_new = TABLE_NAME_MAP_HBOOK[table_name_orig]

    if Format == FORMAT_ROOT:
        if table_name_orig in TABLE_NAME_MAP_ROOT:
            table_name_new = TABLE_NAME_MAP_ROOT[table_name_orig]

    return table_name_new

    # end check_table_name

def get_outfile_name(input_args,config):

    outfile    = input_args.outfile
    table_name = config.table_name 
   
    if outfile :  return outfile

    # default name is based on name of table
    outfile = f"{OUTFILE_PREFIX}_{table_name}.{OUTFILE_SUFFIX}"

    if input_args.OBS:
        outfile = f"{OUTFILE_PREFIX}_obs.{OUTFILE_SUFFIX}"

    if input_args.outlier_fit:
        outfile = f"{OUTFILE_PREFIX}_outlier_fit.{OUTFILE_SUFFIX}"

    if input_args.outlier_sim:
        outfile = f"{OUTFILE_PREFIX}_outlier_sim.{OUTFILE_SUFFIX}"

    return outfile

    # end get_outfile_name

def make_command_sntable_dump(input_args,config):

    # construct command for C program

    outfile_flag = False 
    varlist      = config.varlist
    table_file   = input_args.input_table_file
    outlier_fit  = input_args.outlier_fit
    outlier_sim  = input_args.outlier_sim
    table_name   = config.table_name
    outfile      = config.outfile

    command = f"{Cprogram_dump} {table_file} "
    
    if table_name is not None:
        command += f"{table_name} "
    else:
        return command

    if varlist is not None:
        command += f"-v {varlist} "
        outfile_flag = True

    if input_args.NEVT :
        command += f" NEVT "  

    if outlier_fit :
        outlier_list = " ".join(map(str,outlier_fit))
        command += f"--outlier_fit {outlier_list} "
        outfile_flag = True

    if outlier_sim :
        outlier_list = " ".join(map(str,outlier_sim))
        command += f"--outlier_sim {outlier_list} "
        outfile_flag = True

    if input_args.OBS:
        command += "obs "
        outfile_flag = True
        
    if outfile_flag :
        command += f"-o {outfile} "

    return command 
    # end make_command_sntable_dump


def append_fitres(input_args,config):

    # combine args.append file with config.outfile

    input_table_file = input_args.input_table_file
    append_file      = input_args.append_file
    varlist          = (input_args.varlist).replace(',',' ')
    outfile_dump     = config.outfile
    table_name       = config.table_name
    Format           = config.Format
    combine_arg      = Format[0]  # H or R to make HBOOK or ROOT file

    tmp = append_file.split('.')
    outfile_prefix = f"sntable_append_{tmp[0]}"
    outfile        = f"{outfile_prefix}.text"

    cmd = f"{Cprogram_combine} {append_file} {outfile_dump} {combine_arg} " \
          f"-outprefix {outfile_prefix} " \
          f"t"   # <== write only text output; no ROOT or hbook

    os.system(cmd)

    if input_args.VERBOSE:
        print(f" Finished {Cprogram_combine} with append mode.")
        t0 = config.t0
        print_elapse_time(t0)
        sys.stdout.flush()

    # append comments at top of appended fitres file
    comment_list = []
    comment_list.append(f"# Appended variables '{varlist}' ")
    comment_list.append(f"# from {input_table_file}  (table {table_name})")
    comment_list.append(f"# ")

    comment_file = "tmp_comments_append.text"
    with open(comment_file,"wt") as f:
        for comment in comment_list:
            f.write(f"{comment}\n")
    
    # catenate comments and append file
    cmd_cat = f"cat {outfile} >> {comment_file} ; mv {comment_file} {outfile}"
    os.system(cmd_cat)

    if input_args.VERBOSE:
        print(f" Finished adding comments to top of appended file.")
        t0 = config.t0
        print_elapse_time(t0)
        sys.stdout.flush()

    # end append_fitres

def print_elapse_time(t0):
    dt_sec = (datetime.datetime.now() - t0).total_seconds()
    print(f"\t elapse time: {dt_sec} sec")

# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":
    
    input_args   = get_args()
    config       = Namespace()

    if input_args.VERBOSE :
        print(f"\n# =============== BEGIN =============== ")
        command = ' '.join(sys.argv)
        print(f" Full python command: \n{command} \n")
        t0 = datetime.datetime.now()
        config.t0 = t0
        sys.stdout.flush()

    if input_args.list_tables: list_tables()

    # insert CCID as first variable
    config.varlist = insert_ccid_varlist(input_args)

    # check file format: root, hbook, ...
    config.Format = get_file_format(input_args.input_table_file)

    if input_args.VERBOSE :
        print(f"  File format is {config.Format}")
        print_elapse_time(t0)
        sys.stdout.flush()

    # make sure table name is appropriate for file format
    config.table_name = check_table_name(input_args.table_name,config.Format)

    config.outfile = get_outfile_name(input_args,config)

    # construct command for C code
    config.command = make_command_sntable_dump(input_args,config)

    if input_args.VERBOSE :
        print(f"\n C-program command = \n    {config.command}")
        print_elapse_time(t0) 
        print('\n')
        sys.stdout.flush()

    istat = os.system(config.command)

    if input_args.VERBOSE:
        print_elapse_time(t0) 

    if istat != 0 : sys.exit(istat)

    # check option to combine_fitres files
    if input_args.append_file : 
        append_fitres(input_args,config)

    sys.exit(0)

    # === END ===
