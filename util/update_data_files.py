#!/usr/bin/env python 
#
# Created Oct 2020 [replace old perl script update_data_files.pl
#
# Update text-formatted data files with following options:
#   + change sing of vpec correction 
#   + replace single char filter names with [fullName]-[char]
#
# Sep 29 2022 RK - add new optoins
#   --kcor_file option to replace filtchar with full filt name.
#   --rename option for substrings in file names
# Oct 18 2022 RK  add --filters_remove option

import os, sys, argparse, shutil, gzip, re

SNDATA_ROOT        = os.environ['SNDATA_ROOT']
DATAKEY_VPEC       = "VPEC:"
DATAKEY_FILTERS    = "FILTERS:"
DATAKEY_VARLIST    = "VARLIST:"
DATAKEY_NOBS       = "NOBS:"
DATAKEY_OBS        = "OBS:"

# =====================================
def get_args():

    parser = argparse.ArgumentParser()
    
    msg = "Data version to modify (include full path)"
    parser.add_argument("-v", "--version", help=msg, type=str, default="")

    msg = "option to change sign of VPEC"
    parser.add_argument("--MINUS_VPEC", help=msg, action="store_true")

    msg = "kcor input file to get full filter names"
    parser.add_argument("-k", "--kcor_file", help=msg, type=str, default=None)

    msg = "file name substring to rename; e.g, SWFIT,SOUSA"
    parser.add_argument("-r", "--rename", help=msg, type=str, default=None)

    msg = "list of filters to remove; e.g., XNW"
    parser.add_argument("--filters_remove", help=msg, type=str, default=None)

    if len(sys.argv) == 1:  
        parser.print_help(); sys.exit()

    args = parser.parse_args()

    args.rename_list = []
    if args.rename:
        args.rename_list = args.rename.split(',')

    return args

    # end get_args

# ---------------------------------
def get_data_dir(args):
    version_user = args.version
    if '/' in version_user :
        path_data = version_user
        version   = os.path.basename(version_user)
    else:
        path_data = SNDATA_ROOT + '/lcmerge/' + version_user
        version   = version

    # append args as if user entered version and path separately
    args.version    = version
    args.path_data  = path_data

    print(f" Data path: \n   {path_data}")

    return path_data, version
    # end get_data_dir

# ------------
def get_list_data(path_data):
    version   = os.path.basename(path_data)
    LIST_FILE = (f"{path_data}/{version}.LIST")
    #print(f" xxx list_file = {list_file}")

    list_data_files = []
    with open(LIST_FILE,"rt") as f:
        list_data_files = f.read().splitlines()

    #print(f" xxx {list_data_files} ")
    return list_data_files
    # end get_list_data
    
def update_list_file(version, rename_list, list_data_files_remove):

    # use linux sed to modify file names in [VERSION].LIST file
    # inputs:
    #  + version      data version to update [VERSION].LIST file
    #  + rename_list  replace rename_list[0] with rename_list[1]
    #                  in each data file name
    #  + list_data_files_remove : 
    #       list original file names to remove from [VERSION].LIST

    LIST_FILE = (f"{version}/{version}.LIST")

    substr_orig = rename_list[0]
    subtrs_new  = rename_list[1]

    if len(list_data_files_remove) > 0:
        for data_file in list_data_files_remove:
            cmd = f"sed -i '/{data_file}/d' {LIST_FILE}"
            os.system(cmd)

    cmd = f"sed -i 's/{substr_orig}/{subtrs_new}/g' {LIST_FILE}"
    os.system(cmd)

    return
    # end update_list_file

def make_outdir(path_data):
    # strip version from path_data
    version   = os.path.basename(path_data)

    # remove output version dir if it already exists
    if os.path.exists(version) :
        shutil.rmtree(version)

    print(f" Create output dir: {version}")
    os.mkdir(version)
    
    # copy auxillary files
    cmd_copy = (f"rsync {path_data}/{version}.* {version}/")

    print(f" Copy auxillary files ... ")
    os.system(cmd_copy)

    return
    # end make_outdir

def get_filter_dict(kcor_file):

    # if kcor file contains row with
    #    FILTER: Bessell-U  Bessell_u.DAT
    # then store
    #    filter_dict['U'] = 'Bessell-U'
    #

    filter_dict = {}

    kcor_file = os.path.expandvars(kcor_file)
    if not os.path.exists(kcor_file):
        sys.exit(f"\n ERROR: cannot find kcor_file = \n\t {kcor_file}")

    print(f" Read filter names from \n   {kcor_file}")

    with open(kcor_file,"rt") as f:
        contents     = f.readlines()

    maxlen = 0
    for line in contents:
        wdlist = line.split()
        #print(f" xxx line = {wdlist}")
        if len(wdlist) < 2: continue
        if wdlist[0] == 'FILTER:' :
            filtname_full = wdlist[1]
            filtchar      = filtname_full[-1]
            filter_dict[filtchar] = filtname_full
            filter_dict[f"nobs_{filtchar}"] = 0  # init global counter

            flen = len(filtname_full)
            if flen > maxlen: maxlen = flen
            print(f"\t Store filter map: {filtchar} -> {filtname_full}")

    filter_dict['maxlen'] = maxlen
    return filter_dict
    # end get_filter_dict

# -----------------------------------------
def update_data_file(data_file, args):

    # Update data_file in outdir
    outdir          = args.version
    path_data       = args.path_data  # path to original data
    filter_dict     = args.filter_dict
    filters_remove  = args.filters_remove
    DO_MINUS_VPEC   = args.MINUS_VPEC
    DO_UPD_FILT     = filter_dict is not None
    
    DATA_FILE_ORIG = (f"{path_data}/{data_file}")
    if os.path.exists(DATA_FILE_ORIG) :
        f_orig = open(DATA_FILE_ORIG,"r")
    else:
        DATA_FILE_ORIG = DATA_FILE_ORIG.strip() + '.gz'
        f_orig = gzip.open(DATA_FILE_ORIG,"r")

    DATA_FILE_OUT  = (f"{version}/{data_file}.gz")

    nline_new  = 0
    iline_NOBS = -9
    NOBS_ORIG = 0; NOBS_NEW = 0

    contents     = f_orig.readlines()
    f_orig.close()

    contents_new = []
    
    if DO_UPD_FILT:
        iwd_band     = -1

    for line_orig in contents :
        line   = line_orig.decode('utf-8')
        wdlist = line.split()

        if len(wdlist) < 2:
            contents_new.append(line)
            nline_new += 1
            continue

        ISKEY_FILTERS = (wdlist[0] == DATAKEY_FILTERS)
        ISKEY_VARLIST = (wdlist[0] == DATAKEY_VARLIST)
        ISKEY_OBS     = (wdlist[0] == DATAKEY_OBS )
        ISKEY_NOBS    = (wdlist[0] == DATAKEY_NOBS )

        if ISKEY_FILTERS and filters_remove :
            f_remove_list = list(filters_remove)            
            filters_orig  = wdlist[1]
            filters_new   = re.sub("|".join(f_remove_list), "", filters_orig)
            line = line.replace(filters_orig,filters_new)

        if ISKEY_VARLIST:
            iwd_band = wdlist.index('BAND')

        if ISKEY_OBS:
            band_orig = wdlist[iwd_band]
            NOBS_ORIG += 1
        else:
            band_orig = "?"

        if filters_remove:
            if band_orig[-1] in filters_remove:
                continue

        if ISKEY_OBS: 
            NOBS_NEW += 1

        if DO_MINUS_VPEC and DATAKEY_VPEC in wdlist :
            j           = wdlist.index(DATAKEY_VPEC)
            vpec_orig   = float(wdlist[j+1])
            vpec_new    = -1.0 * vpec_orig
            wdlist[j+1] = str(vpec_new)
            line        = "  ".join(wdlist) + '\n'

        if DO_UPD_FILT and ISKEY_OBS and len(band_orig) == 1:
            filter_dict[f"nobs_{band_orig}"] += 1
            maxlen    = filter_dict['maxlen'] 
            fullname  = filter_dict[band_orig]
            # each filter name uses max string length for column align
            wdlist[iwd_band] = f"{fullname:<{maxlen}}"
            line   = "  ".join(wdlist) + '\n'

        contents_new.append(line)
        if ISKEY_NOBS:   iline_NOBS = nline_new
        nline_new += 1

    # - - - - - - - - - - - 
    # if NOBS has changed (e.g., remove filter bands), 
    # update line with NOBS key
    #print(f"xxx  NOBS = {NOBS_ORIG} -> {NOBS_NEW}  iline_OBS={iline_NOBS}")
    if NOBS_NEW < NOBS_ORIG:
        line   = contents_new[iline_NOBS]
        wdlist = line.split()
        key    = wdlist[0]
        if key != DATAKEY_NOBS:
            sys.exit(f"\n ERROR: key={key} but expect {DATAKEY_NOBS}\n" \
                     f"line = {line}")
        line = line.replace(str(NOBS_ORIG),str(NOBS_NEW))
        contents_new[iline_NOBS] = line

    # - - - - - - 
    # do NOT update data file if no bands pass filter cut
    if NOBS_NEW == 0:
        print(f"  WARNING: no output obs in {data_file} " \
              f"-> do not write data file")
        return False
        

    if args.rename is None:
        out_file = DATA_FILE_OUT
    else:
        # make string substituion in output file name
        substr_orig = args.rename_list[0]
        subtrs_new  = args.rename_list[1]
        out_file = DATA_FILE_OUT.replace(substr_orig,subtrs_new)

    with gzip.open(out_file, "wt") as o :
        o.write("".join(contents_new))

    return True

    # end update_data_file

def print_summary(args):
    
    if args.MINUS_VPEC :
        print(f"   --> VPEC signs are flipped")

    if args.filter_dict is not None:
        print(f"   --> each filter char replaced with full filter name")
        for key, val in args.filter_dict.items():
            if 'nobs' in key:
                txt_warn = ""
                if val == 0 : txt_warn = "# <=== WARNING"
                print(f"\t {key} = {val:5d}  {txt_warn}")

    return
    # end print_summary

# =============================================
#       MAIN
# =============================================
if __name__ == "__main__":

    # read command-line arguments
    args  = get_args()

    path_data, version = get_data_dir(args)

    list_data_files = get_list_data(path_data)

    make_outdir(path_data)

    if args.kcor_file is None:
        args.filter_dict = None
    else:
        args.filter_dict = get_filter_dict(args.kcor_file)

    nupd = 0
    n_data_file = len(list_data_files)
    list_data_files_remove = []
    for data_file in list_data_files :
        is_write = update_data_file(data_file, args)
        if is_write:
            nupd += 1
        else:
            list_data_files_remove.append(data_file)

    # if file names are modified, modify names in [VERSION].LIST file
    if args.rename:
        update_list_file(version,args.rename_list,list_data_files_remove)

    print(f"\n Done modifying {nupd} of {n_data_file} " \
          f"data files under {version}\n")

    # - - - - - - 
    print_summary(args)

    # END main


