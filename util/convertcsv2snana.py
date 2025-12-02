#!/usr/bin/env python
#
# Translate csv format into snana-key format.
# If fisrt column is CID,GAL,STARTID,ROW ...,  then translate
# without modification; otherwise add ROW column for 
# SNANA-KEY format
#
# Apr 27 2021: replace '' with VOID (works only with first '')
# Mar 05 2024: 
#  + write "GAL:" row key if GALID is in header
#  + fix to work with either comma-sep or space-sep var list in header.
#
# Dec 01 2025: 
#   + replace blanks and dash (-) with -9
#   + if valid row key (e.g. SNID) is not in column 0, them automatically move it
#     to row 0
#
import os, sys, argparse, csv

VALID_ROWID_LIST = [ "CID", "SNID", "snid", "ROW", "GALID", "STARID" ]
ROWID_DEFAULT    = "ROW" 

ROWKEY = "ROW:"  # default row key

BLANK_LIST = [ '', '-' ]  # replace this nulls with VALUE_BLANK
VALUE_BLANK = '-9'   # replace blanks with this value

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

    # - - - - - -                                                               
    # if varname_list has only one item,                                        
    #     e.g., [ 'GALID ZPHOT ZPHOTERR' ], then split it into 
    # list of variables. This will enable reading both comma-sep 
    # and space-sep list of variables in header.                   
    if len(varname_list) == 1 :
        varname_list = varname_list[0].split()

    return varname_list, row_list
    # end read_csv_file

def get_colnum_rowid(varname_list):

    colnum = -9
    for valid_rowid in VALID_ROWID_LIST:
        if valid_rowid in varname_list:
            colnum = varname_list.index(valid_rowid)

    return colnum

# ==============================
def  write_out_file(out_file,varname_list,row_list):

    # check if first varname is a valid IDentifier

    add_rownum = False 
    # xxx varname0 = varname_list[0]

    colnum_rowid = get_colnum_rowid(varname_list)

    if colnum_rowid >= 0:
        
        varname_rowid = varname_list[colnum_rowid]
        if varname_rowid == "GALID": 
            global ROWKEY ; ROWKEY = "GAL:"

        if colnum_rowid == 0:
            header = varname_list
        else:
            print(f" Move {varname_rowid} column from colnum={colnum_rowid} to colnum=0")
            rowid =  varname_list.pop(colnum_rowid)
            del      varname_list[colnum_rowid] 
            header = [ rowid ] + varname_list
            #sys.exit(f"\n xxx header = \n{header}")
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
    for row_orig in row_list:

        row = row_orig
        rownum += 1
        line  = f"{ROWKEY} "
        if add_rownum : 
            line += (f"{rownum:3d}  ")

        if colnum_rowid > 0:  # move rowid value to 1st column
            rowid =  row.pop(colnum_rowid)
            del      row[colnum_rowid] 
            row    = [ rowid ] + row

        for blank in BLANK_LIST:
            if blank in row: row = [ VALUE_BLANK  if x == blank else x for x in row]

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

    write_out_file(out_file, VARNAME_LIST, row_list)
    print(" Done.\n")

# END
