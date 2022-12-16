#!/usr/bin/env python

# Created Oct 14 2022 by R.Kessler
# For $USER, find snana-related folders and tar files under
#    * all sim data paths in $SNDATA_ROOT/SIM/PATH_SNDATA_SIM.LIST
#    * $PIPPIN_OUTPUT
#    * ENVs of the form $SNANA_XXX_ROOT
#    * ENVs of the form $SNANA_XXX_USERS
#
# Use linux 'find' command to locate all folders owned by $USER,
# and to sum total disk usage. An output summary file is created
# with a summary by path (path n_folder size) and then a date-sorted
# list of all folders.
#
#  syntax:
#    find_my_snana_folders.py        [locates folders owned by you]
#    find_my_snana_folders.py --user <somebodyElse>
#
#
# Dec 14 2022: search tar files in addition to folders
#
# ==========================

import os, sys, glob, argparse, shutil, subprocess

USER_DEFAULT  = os.environ['USER']
SNDATA_ROOT   = os.environ['SNDATA_ROOT']
PATH_SIMDATA_DEFAULT = '$SNDATA_ROOT/SIM'  # patch expanded later

PATH_PIPPIN_OUTPUT   = '$PIPPIN_OUTPUT'

# define file containing list of paths to search for sim data
PATH_SNDATA_LIST_FILE =  f"{SNDATA_ROOT}/SIM/PATH_SNDATA_SIM.LIST"

# ==================================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "user name"
    parser.add_argument("--user", help=msg,
                        nargs="?", default=None)

    # parse it
    args = parser.parse_args()
    
    if args.user == None:
        args.user = USER_DEFAULT

    return args
    # end get_args


def get_snana_path_list():

    # return list of all SNANA-related paths to check

    path_list = []

    # - - - -
    # start with pippin output
    path_list.append(PATH_PIPPIN_OUTPUT)

    # - - - -
    # next append all possible simdata directories
    path_list.append(PATH_SIMDATA_DEFAULT)
    
    with open(PATH_SNDATA_LIST_FILE,"rt") as f:
        contents = f.readlines()

    for line in contents:
        path_name  = line.rstrip() 
        if len(path_name) > 0 :
            path_list.append(path_name)

    # - - - - 
    # ENVs of the form SNANA_XXX_ROOT
    for env in os.environ:
        if '_SIM' in env: continue
        if 'SNANA_' in env:
            path_list.append(f"${env}")        

    return path_list

    # end get_snana_path_list


def find_my_folders(path,user,s):

    # Inputs:
    #  path:  sim data path to examine
    #  user:  user name to search folders
    #  s   :  pointer to summary file

    print(f"# ------------------------------------------- ")
    print(f" Search folders under  {path}")

    path_exclude = f"{path}/Archive"

    unix_find_dir_base = f"find {path}//* " \
                         f"-maxdepth 0 " \
                         f"-user {user} -type d "

    unix_find_dir_size = f"{unix_find_dir_base} " \
                         f"-exec du -smc {{}} + "

    unix_find_dir_bydate = f"{unix_find_dir_base} " \
                           f"-exec ls -dlt {{}} +"

    # - - - - - - - - - - - -
    print(f"\t fetch number of folders and total size ...")
    ret_size = subprocess.run( [ unix_find_dir_size ], shell=True, 
                               capture_output=True, text=True )
    folder_size_list = (ret_size.stdout).split('\n')
    nerr_size        = len(ret_size.stderr)

    #debug_flag = True
    debug_flag = False
    if debug_flag :
        print(f" xxx ret_size = {ret_size}")
        print(f" xxx folder_size_list = {folder_size_list}")
        print(f" xxx nerr_size = {nerr_size}")

    # - - - - -
    if len(folder_size_list) > 0:
        print(f"\t fetch list sorted by date")
        ret_bydate = subprocess.run( [ unix_find_dir_bydate ], shell=True, 
                                     capture_output=True, text=True )
        folder_bydate_list = (ret_bydate.stdout).split('\n')
        nerr_bydate   = len(ret_bydate.stderr)
    else:
        folder_bydate_list = []
        nerr_bydate = nerr_size

    if folder_bydate_list[-1] == '' :
        del folder_bydate_list[-1]

    # - - - - -
    size_MB = 0.0
    for folder in folder_size_list :
        if 'total' in folder:
            size_MB = float(folder.split()[0])

    size_GB = size_MB/1000.0
    n_folder = len(folder_bydate_list) 
    # xxx if n_folder > 0: n_folder -= 1 # do not count blank line

    # check if string 'scratch' is in path name
    if n_folder > 0:
        path_expand = os.path.expandvars(path)
        str_scratch = ''
        if 'scratch' in path_expand:   str_scratch = 'scratch'
        
        str_line = (f"  - {path:<30} {n_folder:5d} folders   " \
                    f"{size_GB:8.2f}   {str_scratch}")
        s.write(f"{str_line}\n")
        s.flush()

        print(f"\t n_folder={n_folder}  size={size_GB:.2f} GB  " \
              f"  {str_scratch}")


    return folder_bydate_list

    # end find_my_folders

def find_my_tarfiles(path,user,s):

    # Inputs:
    #  path:  sim data path to examine
    #  user:  user name to search tar files
    #  s   :  pointer to summary file

    print(f"# ------------------------------------------- ")
    print(f" Search tar files under {path}")

    path_exclude = f"{path}/Archive"

    # - - -
    unix_find_tar_base = f"find {path}//*.tar* " \
                         f"-user {user} "

    unix_find_tar_size = f"{unix_find_tar_base} " \
                         f"-exec du -smc {{}} + "

    unix_find_tar_bydate = f"{unix_find_tar_base} " \
                           f"-exec ls -dlt {{}} +"

    # - - - - - - - - - - - -
    print(f"\t fetch number of tar files and total size ...")
    ret_size = subprocess.run( [ unix_find_tar_size ], shell=True, 
                               capture_output=True, text=True )
    tarfile_size_list = (ret_size.stdout).split('\n')
    nerr_size         = len(ret_size.stderr)

    #debug_flag = True
    debug_flag = False
    if debug_flag :
        print(f" xxx ret_size = {ret_size}")
        print(f" xxx tarfile_size_list = {tarfile_size_list}")
        print(f" xxx nerr_size = {nerr_size}")
        sys.exit(f"\n xxx DEBUG DIE xxx \n")
    # - - - - -
    if len(tarfile_size_list) > 0:
        print(f"\t fetch list sorted by date")
        ret_bydate = subprocess.run( [ unix_find_tar_bydate ], shell=True, 
                                     capture_output=True, text=True )
        tarfile_bydate_list = (ret_bydate.stdout).split('\n')
        nerr_bydate         = len(ret_bydate.stderr)
    else:
        tarfile_bydate_list = []
        nerr_bydate         = nerr_size

    if tarfile_bydate_list[-1] == '' :
        del tarfile_bydate_list[-1]

    # - - - - -
    size_MB = 0.0
    for tarfile in tarfile_size_list :
        if 'total' in tarfile:
            size_MB = float(tarfile.split()[0])

    size_GB = size_MB/1000.0
    n_tarfile = len(tarfile_bydate_list) 
    # xxx mark if n_tarfile > 0: n_tarfile -= 1 # do not count blank line

    # check if string 'scratch' is in path name
    if n_tarfile > 0:
        path_expand = os.path.expandvars(path)
        str_scratch = ''
        if 'scratch' in path_expand:   str_scratch = 'scratch'
        
        str_line = (f"  - {path:<30} {n_tarfile:5d} tarFiles  " \
                    f"{size_GB:8.2f}   {str_scratch}")
        s.write(f"{str_line}\n")
        s.flush()

        print(f"\t n_tarfile={n_tarfile}  size={size_GB:.2f} GB  " \
              f"  {str_scratch}")

    return tarfile_bydate_list

    # end find_my_tarfiles


def open_summary_file(user):
    summary_file = f"snana_folders_{args.user}.log"
    print(f" Open summary file: {summary_file}")
    s = open(summary_file,"wt")
    s.write("SUMMARY:   # path                     n_item         size(GB)\n")
    s.flush()
    return s, summary_file

def write_item_bydate(item_type, path, dict_item_bydate_list, s):

    if len(dict_item_bydate_list) == 0 : return

    s.write(f"\n# {path} {item_type}\n")
    for item_line in dict_item_bydate_list:
        wdlist = item_line.split()
        if len(wdlist) < 4 : 
            continue
        date             = " ".join(wdlist[5:8])
        item           = wdlist[8]
        item_base      = os.path.basename(item)
        #sys.exit(f"\n xxx wdlist = {wdlist}\n xxx folder = {folder}\n")
        s.write(f"  {date}   {item_base}\n")
    s.flush()
    return
    # end write_item_bydate

# ===================================================
if __name__ == "__main__":

    args   = get_args() 
    
    print(f" Find snana folders for user = {args.user}")

    s, summary_file = open_summary_file(args.user)

    snana_path_list = get_snana_path_list()

    #sys.exit(f"\n xxx snana_path_list = \n{snana_path_list}")

    dict_folder_lists  = {}
    dict_tarfile_lists = {}
    for path in snana_path_list:
        folder_bydate_list  = find_my_folders(path,args.user,s)
        tarfile_bydate_list = find_my_tarfiles(path,args.user,s)
        dict_folder_lists[path]  = folder_bydate_list
        dict_tarfile_lists[path] = tarfile_bydate_list

    # dump out date-sorted folders
    s.write("\n\nFOLDER_PLUS_TAR_LIST:   # sorted by path and date\n")
    s.flush()
    for path in snana_path_list:
        write_item_bydate("folders",  path, dict_folder_lists[path] , s)
        write_item_bydate("tarFiles", path, dict_tarfile_lists[path], s)

    s.close()

    print(f"\n See summary in {summary_file}")

    # == END ===


