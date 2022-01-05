#!/usr/bin/env python
#
# Created Oct 2020 [replace old perl script update_data_files.pl
#
# Update text-formatted data files with following options:
# 
#   update_data_files.pl -V <VERSION>  -u <updFile>  -v <varList>
#     -> update varList variables in updFile
#
#    or
#
#   update_data_files.pl -V <VERSION> --MINUS_VPEC
#     -> change sign of VPEC
#
# BEWARE: only --MINUS_VPEC works now ... 
#         later need to implement updFile and var_list.
#
# BEWARE: works only on gzipped input files
#
import os, sys, argparse, shutil, gzip
import pandas as pd

SNDATA_ROOT = os.environ['SNDATA_ROOT']

# =====================================
def get_args():

    parser = argparse.ArgumentParser()
    
    msg = "Data version to modify"
    parser.add_argument("-V", "--VERSION", help=msg, type=str, default="")

    msg = "fitres file with update values"
    parser.add_argument("-u", "--upd_file", help=msg, type=str, default=None)

    msg = "comma-sep list of variables to update"
    parser.add_argument("-v", "--var_list", help=msg, type=str, default=None)

    msg = "option to change sign of VPEC"
    parser.add_argument("--MINUS_VPEC", help=msg, action="store_true")

    if len(sys.argv) == 1:  parser.print_help(); sys.exit()

    args = parser.parse_args()

    return args

    # end get_args

# ---------------------------------
def get_data_dir(args):
    version_user = args.VERSION
    if '/' in version_user :
        path_data = version_user
        version   = os.path.basename(version_user)
        pass
    else:
        path_data = SNDATA_ROOT + '/lcmerge/' + version_user
        version   = version_user

    return path_data, version

# ---------------------------------
def read_upd_file(args):

    upd_file = args.upd_file
    if upd_file is None: return None

    var_list = args.var_list.split(",")

    if upd_file is None:
        return None

    var_list_local = ['CID'] + var_list
    df  = pd.read_csv(upd_file, comment="#", delim_whitespace=True, 
                      usecols=var_list_local)
    df["CID"] = df["CID"].astype(str)
    df        = df.set_index("CID", drop=False)

    #print(f" xxx upd_df = {df}")
    return df
    
# ------------
def get_list(path_data):
    version   = os.path.basename(path_data)
    LIST_FILE = (f"{path_data}/{version}.LIST")
    #print(f" xxx list_file = {list_file}")

    list_data_files = []
    with open(LIST_FILE,"rt") as f:
        list_data_files = f.read().splitlines()

    #print(f" xxx {list_data_files} ")
    return list_data_files

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

# --------------
def update_data_file(upd_dict):

    # Update data_file in outdir
    # upd_df is the list of variables to update from upd_file.

    outdir     = upd_dict['version']
    path_data  = upd_dict['path_data']
    data_file  = upd_dict['data_file']
    upd_df     = upd_dict['upd_df']
    var_list   = upd_dict['var_list']
    MINUS_VPEC = upd_dict['MINUS_VPEC']

    if MINUS_VPEC :
        print(f"   Change VPEC sign for {data_file}")
    else:
        print(f"   Update {var_list} for {data_file}")

    DATA_FILE_ORIG = (f"{path_data}/{data_file}")
    if not os.path.exists(DATA_FILE_ORIG) :
        DATA_FILE_ORIG = DATA_FILE_ORIG.strip() + '.gz'

    DATA_FILE_OUT  = (f"{version}/{data_file}.gz")
    KEY_SNID       = "SNID:"
    KEY_VPEC       = "VPEC:"

    if upd_df is not None :
        cid_upd_list   = upd_df['CID'].head().tolist()

    nline = 0
    with gzip.open(DATA_FILE_ORIG,"r") as d :
        contents     = d.readlines()
        contents_new = []
        for line_orig in contents :
            nline += 1
            line   = line_orig.decode('utf-8')
            wdlist = line.split()
            if KEY_SNID in wdlist :
                j    = wdlist.index(KEY_SNID)
                SNID = wdlist[j+1]
            
            # TO DO: replace value here from upd_file ???

            if MINUS_VPEC and KEY_VPEC in wdlist :
                j           = wdlist.index(KEY_VPEC)
                vpec_new    = -1.0 * float(wdlist[j+1])
                wdlist[j+1] = str(vpec_new)
                line        = " ".join(wdlist) + '\n'

            contents_new.append(line)


    out_file = DATA_FILE_OUT
    with gzip.open(out_file, "wt") as o :
        o.write("".join(contents_new))

    #sys.exit("\n xxx DEBUG DIE xxx \n")
    
# =============================================
#       MAIN
# =============================================
if __name__ == "__main__":

    # read command-line arguments
    args  = get_args()

    path_data,version = get_data_dir(args)

    upd_df = read_upd_file(args)

    #print(f"\n xxx upd_df = {upd_df}\n")
    list_data_files = get_list(path_data)

    make_outdir(path_data)

    upd_dict = {
        "version"    : version,
        "upd_df"     : upd_df,
        "var_list"   : args.var_list,
        "MINUS_VPEC" : args.MINUS_VPEC,
        "path_data"  : path_data,
        "data_file"  : "bla"
    }

    nupd = 0
    for data_file in list_data_files :
        upd_dict['data_file'] = data_file
        update_data_file(upd_dict)
        nupd += 1

    print(f"\n Done modifying {nupd} data files under {version}\n")

    # END

