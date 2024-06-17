#!/usr/bin/env python
#
# Created Jun 2024 by R.Kessler
#
# Install SNDATA_ROOT from a backup tar file (e.g., from zenodo).
# Prior to using this script, copy (download) SNDATA_ROOT_[data].tar.gz
# to the same directory containing the old SNDATA_ROOT.
#
# Running this script (no args) will do the following:
#  move old SNDATA_ROOT to SNDATA_ROOT_BACKUP_[date]
#  mkdir SNDATA_ROOT
#  unpack most recent SNDATA_ROOT_*.tar.gz file that is in
#    same dir as SNDATA_ROOT (hence no need for command line args)
#  rsync original SIM/ contents back to new $SNDATA_ROOT/SIM
#
# Some steps prompt user before executing
#
# ===========

import os, sys, glob, datetime
from argparse import Namespace

STRING_SNDATA_ROOT  = "SNDATA_ROOT"

SNDATA_ROOT          = os.environ[STRING_SNDATA_ROOT]
SNDATA_ROOT_TOPDIR   = os.path.dirname(SNDATA_ROOT)

# define date stamp to mark backup of current SNDATA_ROOT
tnow       = datetime.datetime.now()
DATE_STAMP = ('%4.4d-%2.2d-%2.2d' % (tnow.year,tnow.month,tnow.day) )

response_yes = 'y'

# ==================================================

def query_install_continue(question, tar_file):
    response = input(f"\n {question}\n y/n ? ")
    if response != response_yes:
        sys.exit(f"\n\t *** ABORT installation of {tar_file} ***")
    
def get_install_tar_file(topdir_search):
    # return name of most recent tar file in same base dir as SNDATA_ROOT.

    wildcard = f"{STRING_SNDATA_ROOT}_*.tar.gz"
    tar_list = sorted(glob.glob1(topdir_search, wildcard))
    n_tar = len(tar_list)
    
    if n_tar == 0 :
        msgerr = f"\nERROR: Found no {wildcard} backups to install under\n\t{topdir_search}"
        sys.exit(msgerr)
        
    tar_file = tar_list[n_tar-1]  # pick more recent -> last in sorted list

    query_install_continue(f"Install {tar_file} in {topdir_search}/{STRING_SNDATA_ROOT}", tar_file)
            
    return tar_file

def move_old_SNDATA_ROOT(info):
    tar_file = info.install_tar_file
    topdir   = SNDATA_ROOT_TOPDIR
    
    dirname_backup = f"{STRING_SNDATA_ROOT}_backup_" + f"{DATE_STAMP}"
    info.dirname_backup = dirname_backup
    
    full_path = f"{topdir}/{dirname_backup}"
    if os.path.exists(full_path):
        print(f"\n backup already exists in {full_path} \n --> Skip backup.")
    else:
        cmd       = f"cd {topdir} ; mv {STRING_SNDATA_ROOT} {dirname_backup}"
        cmd_query = f"cd {topdir} \n mv {STRING_SNDATA_ROOT} {dirname_backup}"
        query_install_continue(cmd_query, tar_file)
        os.system(cmd)

    return

def install(info):

    tar_file = info.install_tar_file
    topdir   = SNDATA_ROOT_TOPDIR
    
    # create SNDATA_ROOT
    if os.path.exists(SNDATA_ROOT) :
        sys.exit(f"\n ERROR: {STRING_SNDATA_ROOT} still exists under {topdir}")

    os.mkdir(SNDATA_ROOT)

    cmd_list=  []
    cmd_cd = f"cd {SNDATA_ROOT}"
    cmd_list.append(f"mv ../{tar_file} .")
    cmd_list.append(f"tar -xzf {tar_file}")
    cmd_list.append(f"mv {tar_file} ../")

    print()
    for cmd in cmd_list:
        print(f"   {cmd}")
        os.system(f"{cmd_cd} ; {cmd}")
    
    return

def restore_SIM(info):

    # since SNDATA_ROOT backup has nothing in SIM/ subdir,
    # restore local SIM/ contents.
    
    dirname_backup       = info.dirname_backup
    topdir               = SNDATA_ROOT_TOPDIR
    subdir_sim           = "SIM"
        
    full_path  = f"{topdir}/{dirname_backup}/SIM"
    cmd = f"rsync -rt {full_path}/* {SNDATA_ROOT}/SIM/"
    print(f"\n Restore old contents of {STRING_SNDATA_ROOT}/{subdir_sim}/")
    os.system(cmd)
    
    return

# ==================================================
if __name__ == "__main__":

    info = Namespace()
    info.install_tar_file = get_install_tar_file(SNDATA_ROOT_TOPDIR)
    
    move_old_SNDATA_ROOT(info)

    install(info)

    restore_SIM(info)

    print(f"\n Done with {STRING_SNDATA_ROOT} intall.")
    
    # == END: ===
    
 

