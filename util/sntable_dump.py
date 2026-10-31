#!/usr/bin/env python
#
# Created March 2021 by R.Kessler
#   [re-write of perl script sntable_dump.pl]
#
# Utility to examine and dump analysis info from any SNTABLE
# root (tree) format.
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
#
# May 30 2024
#   + Used program_exists to ensure exectuable exists.
#
# Sep 24 2026:  allow comma-sep or space sep varlist
# Oct 30 2026: 
#   +  add new arg --snana_dir to use private Cprog_dump and Cprog_combine
#   +  append header comments before VARNAMES instead of at top of file.
#
# ===================

import os, sys, argparse, yaml, datetime 
from   argparse import Namespace
import time

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

SNANA_DIR_DEFAULT        = os.getenv('SNANA_DIR')
Cprog_dump_base       = f"sntable_dump.exe"
Cprog_combine_base    = f"combine_fitres.exe"

Cprog_dump    = None
Cprog_combine = None

T_WAIT_MAX = 500.00 # Max time to wait for executable to exist

# ===========
def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "name of table file to examine"
    parser.add_argument("input_table_file", help=msg, nargs="?", default=None)

    msg = "name of table (e.g., SNANA, FITRES, SPEC ...)"
    parser.add_argument("table_name", help=msg, nargs="?", default=None)

    msg = "List possible table names"
    parser.add_argument("-l", "--list_tables", help=msg, action="store_true")

    msg = "list of variables to extract (comma-sep or space-sep)"
    parser.add_argument("-v", "--varlist", help=msg, nargs='+', 
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

    msg = f"SNANA_DIR for {Cprog_dump_base} and {Cprog_combine_base} " \
          f"(default = $SNANA_DIR = {SNANA_DIR_DEFAULT})"
    parser.add_argument("--snana_dir", help=msg, type=str, default=SNANA_DIR_DEFAULT)

    args = parser.parse_args()

    global Cprog_dump, Cprog_combine
    Cprog_dump    = f"{args.snana_dir}/bin/{Cprog_dump_base}"
    Cprog_combine = f"{args.snana_dir}/bin/{Cprog_combine_base}"

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
    # Sep 24 2026: allow comma-sep or space sep varlist

    # bail if CCID is already in varlist
    varlist = args.varlist
    if varlist is None : return None


    varlist_string = ' '.join(varlist)

    if VARNAME_CCID in varlist: return varlist

    varlist_new = f"{VARNAME_CCID},{varlist_string}"

    varlist_new = varlist_new.replace(',',' ')

    # replace commas with spaces for C program
    return varlist_new

    # end insert_ccid_varlist

def get_file_format(input_table_file):

    Format = None
    
    for suffix in [ 'ROOT', 'root' ]:
        if suffix in input_table_file: Format = FORMAT_ROOT

    if Format is None :
        sys.exit(f"\n ERROR: unknown format for " \
                 f"input_table_file = {input_table_file}\n")

    return Format
    # end get_file_format

def check_table_name(table_name_orig,Format):

    table_name_new  = table_name_orig


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

    command = f"{Cprog_dump} {table_file} "
    
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


def program_exists(program_name, t_wait_max, abort_not_exists):
    # Created May 23 2024 by R.Kessler  
    # Return True if profram_name exists within time t_wait_max (seconds).
    # If program_name does not exist after t_wait_time,
    # and abort_not_exists=True, the abort.

    t_wait_tot = 0
    t_wait     = 10  # wait this long before checking again
    
    while os.access(program_name, os.X_OK) is False:
        t_now   = datetime.datetime.now()
        tstr    = t_now.strftime("%Y-%m-%d %H:%M:%S")
        print(f"  waiting for {program_name} to exist ({tstr})")
        sys.stdout.flush()
        time.sleep(t_wait)
        t_wait_tot += t_wait
        if t_wait_tot > t_wait_max :
            print(f"WARNING: Could not find program {program_name}")
            print(f"\tafter waiting {t_wait_tot} seconds.")
            sys.stdout.flush()
            if abort_not_exists:
                log_assert(False, ["ABORT on missing program", program_name] )
            else:
                return False

    return True


    
def append_fitres(input_args,config):

    # combine args.append file with config.outfile
    # Oct 30 2026: remove combine_arg=R

    input_table_file = input_args.input_table_file
    append_file      = input_args.append_file

    # xxx mark varlist          = (input_args.varlist).replace(',',' ')
    varlist          = ' '.join(input_args.varlist).replace(',',' ')

    outfile_dump     = config.outfile
    table_name       = config.table_name
    Format           = config.Format
    # xxx mark combine_arg      = ''  # Format[0]  # R for ROOT file

    tmp = append_file.split('.')
    outfile_prefix = f"sntable_append_{tmp[0]}"
    outfile        = f"{outfile_prefix}.text"

    cmd = f"{Cprog_combine} {append_file} {outfile_dump} " \
          f"-outprefix {outfile_prefix} " \
          f"t"   # <== write only text output; no ROOT 

    # .check that Cprog_combine exists ...
    program_exists(Cprog_combine, T_WAIT_MAX, True)
    
    os.system(cmd)

    if input_args.VERBOSE:
        print(f" Finished {Cprog_combine} with append mode.")
        t0 = config.t0
        print_elapse_time(t0)

    # append comments at top of appended fitres file
    comment_list = []
    comment_list.append(f"# ")
    comment_list.append(f"# APPEND_VARLIST:   {varlist} ")
    comment_list.append(f"# APPEND_SOURCE:    {input_table_file} ")
    comment_list.append(f"# APPEND_PROGRAM:   {Cprog_combine}")
    comment_list.append(f"# ")

    sed_insert_text = ''
    for i, comment in enumerate(comment_list):
        if i == 0 :
            sed_insert_text = f"i {comment} "
        else:
            sed_insert_text += f"\\n{comment} "

    # Oct 30 2026: prepare unix sed command to insert comments before VARNAMES key
    sed_insert_command = f"sed -i '/VARNAMES/{sed_insert_text}' {outfile}"

    REFAC_INSERT_COMMENTS = True
    if REFAC_INSERT_COMMENTS:
        print(f"\n xxx sed_insert_command = \n{sed_insert_command}")
        os.system(sed_insert_command)
    else:
        # xxxxxxxx mark delete Oct 30 2026 xxxxxxxx
        # legacy comments that appear at the top of file
        comment_file = "tmp_comments_append.text"
        with open(comment_file,"wt") as f:
            for comment in comment_list:
                f.write(f"{comment}\n")
    
        # catenate comments and append file
        cmd_cat = f"cat {outfile} >> {comment_file} ; mv {comment_file} {outfile}"
        os.system(cmd_cat)
        # xxxxxxxx end mark xxxxxxxx

    # - - - - - - 
    if input_args.VERBOSE:
        print(f" Finished adding comments to top of appended file.")
        t0 = config.t0
        print_elapse_time(t0)

    # end append_fitres

def print_elapse_time(t0):
    dt_sec = (datetime.datetime.now() - t0).total_seconds()
    print(f"\t elapse time: {dt_sec} sec")
    sys.stdout.flush()
    
# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":
    
    input_args   = get_args()
    config       = Namespace()

    #sys.exit(f"\n xxx Cprog_dump = {Cprog_dump}")

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

    # check file format: root ...
    config.Format = get_file_format(input_args.input_table_file)

    if input_args.VERBOSE :
        print(f"  File format is {config.Format}")
        print_elapse_time(t0)

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

    # .xyz check that Cprogram exists ...
    program_exists(Cprog_dump, T_WAIT_MAX, True)
    
    istat = os.system(config.command)

    if input_args.VERBOSE:
        print_elapse_time(t0) 

    if istat != 0 : sys.exit(istat)

    # check option to combine_fitres files
    if input_args.append_file : 
        append_fitres(input_args,config)

    sys.exit(0)

    # === END ===
