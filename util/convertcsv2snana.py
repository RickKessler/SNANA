#!/usr/bin/env python
#
# Translate csv format into snana-key format.
# If fisrt column is CID,GAL,STARTID,ROW ...,  then translate
# without modification; otherwise add ROW column for 
# SNANA-KEY format
#
# Apr 27 2021: replace '' with VOID (works only with first '')
#

import os, sys, argparse, csv

VALID_ROWID_LIST = [ "CID", "SNID", "ROW", "GALID", "STARID" ]
ROWID_DEFAULT    = "ROW" 

# =====================================
def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "name of input csv file"
    parser.add_argument("-i", "--csv_file", help=msg, type=str, default=None)

    msg = "optional name of output text file"
    parser.add_argument("-o", "--out_file", help=msg, type=str, default=None)

    msg = "write output VARNAMES with all CAPS"
    parser.add_argument("-C", "--CAPS", help=msg, action="store_true")

    args = parser.parse_args()

    if len(sys.argv) == 1:  parser.print_help(); sys.exit()

    if args.csv_file is None :
        sys.exit("ERROR: must specifiy input csv file")

    return args

# ====================================
def get_out_file_name(args):
    if args.out_file is not None :
        out_file = args.out_file
    else:
        infile   = (args.csv_file).split(".")[0]
        out_file = infile + ".text"
        
    print(f" Write SNANA-key format to {out_file}")
    return out_file

# =================================
def read_csv_file(csv_file):

    print(f" Read csv contents from: {csv_file}")

    nrow = 0; row_list = []
    with open(csv_file, 'r') as f:  
        contents = csv.reader(f)
        for row in contents :
            if len(row) == 0 : continue
            if row[0][0] == '#' : continue
            if nrow == 0 :
                varname_list = row
            else:
                row_list.append(row)
            nrow += 1

            if '' in row:
                j = row.index('')
                row[j] = 'VOID'

    #print(f" xxx varname_list = {varname_list}")
    #print(f" xxx row_list = {row_list}")

    return varname_list, row_list

# ==============================
def  write_out_file(out_file,varname_list,row_list):

    # check if first varname is a valid IDentifier

    add_rownum = False 
    varname0 = varname_list[0]

    #print(f" xxx varname_list = {varname_list}")
    #print(f" xxx varname0 = '{varname0}'")
    #print(f" xxx VALID_ROWID_LIST = {VALID_ROWID_LIST}")

    if varname0 in VALID_ROWID_LIST:
        header = varname_list
    else:
        header = [ ROWID_DEFAULT ] + varname_list
        add_rownum = True


    f = open(out_file,"wt")

    # write command in comment fields
    command = " ".join(sys.argv)
    f.write(f"# Created with command\n")
    f.write(f"#   {command} \n#\n")

    # write list of VARNAMES
    line = "VARNAMES: " + " ".join(header)
    f.write(f"{line}\n")

    # write each row
    rownum = 0 
    for row in row_list:
        rownum += 1
        line  = "ROW: "
        if add_rownum :  line += (f"{rownum:3d}  ")
        line += " ".join(row)
        f.write(f"{line}\n")

    f.close()


# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":

    # read command-line arguments
    args  = get_args()

    # read contents of csv file
    varname_list, row_list = read_csv_file(args.csv_file)

    out_file = get_out_file_name(args)

    # check all caps option
    if args.CAPS :
        VARNAME_LIST = [v.upper() for v in varname_list ]
    else:
        VARNAME_LIST = varname_list

    write_out_file(out_file,VARNAME_LIST,row_list)
    print(" Done.\n")

# END
