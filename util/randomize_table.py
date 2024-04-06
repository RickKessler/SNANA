#!/usr/bin/env python
#
# Created Feb 2024
# read SNANA-formatted table file and re-write table rows in random order.
# Works on HOSTLIBs and FITRES files.
# Motivations:
#   + HOSTLIBs are sometimes created ordered by redshift or coordinates,
#     and thus processing the first "HOSTLIB_MAXREAD" events results in
#     heavily biased selection.
#   + FITRES files from BBC are ordered by survey, so again processing
#     the first N entries can result in excluding some of the surveys.
#
# usage:
#   randomize_table.py <tableFile>
# ->
#   outputs <tableFile>_randomized.
#
# Apr 6 2024: add prescale option to create smaller table
#             (e.g., reduced HOSTLIB)
#
# ===================================================

import os, sys, argparse, datetime, gzip
import random

tnow       = datetime.datetime.now()
DATE_STAMP = ('%4.4d-%2.2d-%2.2d' % (tnow.year,tnow.month,tnow.day) )
USERNAME         = os.environ['USER']

NROW_STDOUT_UPDATE = 200000

# ============
def get_args():
    parser = argparse.ArgumentParser()

    msg = "name of table file (FITRES or HOSTLIB) to randomize"
    parser.add_argument("input_table_file", help=msg, type=str, default=None)

    msg = "prescale (default=1)"
    parser.add_argument("--prescale", help=msg, type=int, default=1)

    args = parser.parse_args()

    return args
    # end get_args

def read_input(input_file):
    print(f" Process {input_file}")

    if '.gz' in input_file:
        f = gzip.open(input_file,"rt")
    else:
        f = open(input_file,"rt")
        
    input_lines = f.readlines()
    f.close()
    
    return input_lines

def write_output(args, input_lines, output_file):

    print(f" Write randomized table to {output_file}")
    nrow_tot = 0
    nrow_wr  = 0
    prescale = args.prescale
    base_code = os.path.basename(sys.argv[0])
    
    f = open(output_file,"wt") 

    f.write(f"# {base_code} with prescale={prescale} run by {USERNAME} on {DATE_STAMP}\n\n")

    for line in input_lines:
        line = line.rstrip()
        if is_rowkey(line): 
            f.write(f"\n")
            break
        f.write(f"{line}\n")

    # now randomize everything, including header stuff,
    # and then write only lines with valid row key
    
    print(f"\t randomize using random.shuffle ... ")
    sys.stdout.flush()
    random.shuffle(input_lines)

    print(f"\t write output ... ")
    sys.stdout.flush()

    for line in input_lines :
        line = line.rstrip()
        if is_rowkey(line):
            nrow_tot += 1
            if nrow_tot % args.prescale != 0 : continue
            
            f.write(f"{line}\n")
            nrow_wr += 1
            if (nrow_wr % NROW_STDOUT_UPDATE) == 0 :
                print(f"\t Write table row {nrow_wr}  (read {nrow_tot})")
                sys.stdout.flush()
    # - -  -
    f.close()
    return nrow_wr


def is_rowkey(line):
    # return True if this line is a table row that begins
    # with a standard row key.

    if len(line) == 0 : return False
    if line[0] == '#' : return False
    rowkey_list = [ 'SN:', 'GAL:', 'ROW:' ]
    for key in rowkey_list:
        if key in line : return True

    return False
    # end is_rowkey

# =================================================
# =================================================

if __name__ == "__main__":
    args   = get_args() 

    input_file  = args.input_table_file

    base_table_file = os.path.basename(args.input_table_file)
    output_file     = base_table_file.split('.gz')[0] + '_randomized' 

    input_lines = read_input(input_file)

    nrow = write_output(args, input_lines, output_file)
    
    print(f'\n Done writing {nrow} randomizing table rows in {output_file}')

    # === end main ===
