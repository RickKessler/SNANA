#!/usr/bin/env python
#
# Interface to numerous snana.exe commands that are hard to find/remember.
# Use "quick_commands.py -H" to get explicit examples on how to combine 
# arguments to perform specific tasks.
#
# =========================

import os, sys, argparse
import pandas as pd

# ----------------
snana_program = "snana.exe"

LOG_FILE = "quick_command.log"

HELP_COMMANDS = f"""
# translate TEXT format (from $SNDATA_ROOT/lcmerge) to FITS format
 quick_commands.py -v CFA3_KEPLERCAM --version_reformat_fits CFA3_FITS

# translate TEXT format from local dir to FITS format:
  quick_commands.py -v DES-SN3YR_LOWZ  -p ./  \
  --version_reformat_fits LOWZ_FITS

# extract a few events from FITS format back to TEXT format:
 quick_commands.py -v SDSS_allCandidates+BOSS --cidlist_text 1032,5944
 quick_commands.py -v SDSS_allCandidates+BOSS --cidlist_text TEN

# create SIMLIB file from LOWZ data
 quick_commands.py -v DES-SN3YR_LOWZ \
   -k $SNDATA_ROOT/kcor/DES/DES-SN3YR/kcor_LOWZ.fits \
   --simlib_outfile DES-SN3YR_LOWZ.SIMLIB

# make human-readable table of redshifts and vpec for list of CIDs
 quick_commands.py -v DES-SN3YR_LOWZ --cidlist_ztable 2002de,2007cq,2010dt
 quick_commands.py -v DES-SN3YR_LOWZ --cidlist_ztable TEN

# make SNANA-formatted table(s)
 quick_commands.py -v DES-SN3YR_LOWZ --sntable_list SNANA
 quick_commands.py -v DES-SN3YR_LOWZ --sntable_list 'SNANA(text:host)'
 quick_commands.py -v DESY1_forcePhoto_fake_snana_fits -p $DES_ROOT/lcmerge \
    --sntable_list 'SNANA OUTLIER(nsig:10)'

# extract info about a photometry version
  quick_commands.py -v SDSS_allCandidates+BOSS --get_info_phot

# extract info about SNANA code
  quick_commands.py --get_info_code

"""
# =============================


# =====================================
def get_args():
    parser = argparse.ArgumentParser()
    
    msg = "HELP on input config file"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of data version"
    parser.add_argument("-v", "--version", help=msg, type=str, default="")

    msg = "name of PRIVATE_DATA_PATH (if not a default path)"
    parser.add_argument("-p", "--path", help=msg, type=str, default=None)

    msg = "output FITS-format version"
    parser.add_argument("--version_reformat_fits",help=msg,type=str,default="")

    msg = "name of kcor/calib file (e.g., for SIMLIB_OUT)"
    parser.add_argument("-k", "--kcor_file", help=msg, type=str, default="")
    msg = "output SIMLIB file created from data"
    parser.add_argument("--simlib_outfile",help=msg,type=str,default="")

    msg = "Table(s) to create for all events"
    parser.add_argument("--sntable_list",help=msg,type=str,default="")

    msg = "CID list to extract into text format"
    parser.add_argument("--cidlist_text",help=msg,type=str,default="")

    msg = "CID list to make table of redshifts and vpec"
    parser.add_argument("--cidlist_ztable", help=msg, type=str, default="")

    msg = "get info for photometry version"
    parser.add_argument("--get_info_phot", help=msg, action="store_true")

    msg = "get info for code version"
    parser.add_argument("--get_info_code", help=msg, action="store_true")

    msg = "convert this simgen_dump file to fitres file for SALT2mu"
    parser.add_argument("--simgen_dump_file", help=msg, type=str, default="")

# SIMLIB_OUT ...

    #msg = "comma-sep list of CIDs (for list args)"
    #parser.add_argument("--c", help=msg, type=str, default="")

    #msg = "Use each SN instead of BBC binning"
    #parser.add_argument("-u", "--unbinned", help=msg, action="store_true")


    if len(sys.argv) == 1:  parser.print_help(); sys.exit()

    args = parser.parse_args()

    return args

    # end get_args

def reformat_fits(args):
    # Use SNANA program to translate text format to FITS format
    vout = args.version_reformat_fits
    rmdir_check(vout)

    command = snana_command_plus_version(args)
    command += f"VERSION_REFORMAT_FITS {vout} "
    exec_command(command,args,0)
    # end reformat_fits

def make_ztable(args):
    cidlist     = args.cidlist_ztable
    prefix_temp = "OUT_TEMP_SNANA"
    command = snana_command_plus_version(args)
    command += f"SNTABLE_LIST SNANA "
    command += f"TEXTFILE_PREFIX {prefix_temp} "
    command += arg_cidlist(cidlist)

    # extract info from SNANA table and present in human-readable form
    table_script = "get_fitres_values.py"
    table_file   = f"{prefix_temp}.SNANA.TEXT"
    zvar_list    = "zHEL,zHELERR,zCMB,zCMBERR,zHD,zHDERR,VPEC,VPECERR"
    cmd_table    = f"{table_script} -f {table_file} --nrow 10 -v {zvar_list}"

    command_list = [ command, cmd_table ]

    istat0, istat1 = exec_command(command_list,args,0)

    # end make_ztable

def extract_text_format(args):

    vin     = args.version
    vout    = f"{vin}_TEXT"
    cidlist = args.cidlist_text

    rmdir_check(vout)
        
    print(f"\n Create new data folder: {vout}")

    command = snana_command_plus_version(args)
    command += f"VERSION_REFORMAT_TEXT {vout} "
    command += arg_cidlist(cidlist)
    exec_command(command,args,0)
    # end extract_text_format

def make_simlib(args):

    kcor_file = args.kcor_file
    outfile   = args.simlib_outfile

    command = snana_command_plus_version(args)
    command += f"KCOR_FILE {kcor_file} "
    command += f"SIMLIB_OUTFILE {outfile} "
    exec_command(command, args, 0)

    # end make_simlib

def make_sntable(args):

    sntable_list = args.sntable_list
    if 'FITRES' in sntable_list:
        sys.exit("\n FITRES option not valid here. " \
                 "\n Must run snlc_fit.exe for FITRES table")

    ntail = 0
    if 'OUTLIER' in sntable_list: ntail = 30

    prefix = "SURVEY" # key for snana.exe to substitute with survey name

    command = snana_command_plus_version(args)
    command += f"SNTABLE_LIST '{sntable_list}' "
    command += f"TEXTFILE_PREFIX {prefix} "
    exec_command(command, args, ntail)

    # grep out name of table file from log
    cmd_grep = f"grep Close {LOG_FILE}"
    os.system(cmd_grep)

    # end make_sntable

def arg_cidlist(cidlist):
    # return snana.exe args for cidlist
    if cidlist == "TEN" :
        arg_list = f"MXEVT_PROCESS 10 "
    else:
        arg_list = f"SNCCID_LIST {cidlist} "
        
    return arg_list
    # end arg_cidlist

def get_info_photometry(args):

    command  = f"{snana_program} GETINFO {args.version}"
    exec_command(command,args,9)

    # end get_info_photometry

def get_info_code(args):
    command = f"snana.exe --snana_version"
    exec_command(command,args,4)

def rmdir_check(dir_name):
    if os.path.exists(dir_name) :
        cmd_rm = f"rm -r {dir_name}"
        os.system(cmd_rm)
    

def snana_command_plus_version(args):
    # create part of snana.exe command that includes photometry
    # version and option private data path
    command = f"{snana_program} "
    command += f"NOFILE "
    command += f"VERSION_PHOTOMETRY {args.version} "
    if args.path:
        command += f"PRIVATE_DATA_PATH {args.path} "
    return command
    # end snana_command_plus_version

def exec_command(command,args,ntail):

    if isinstance(command, list):
        cmd_snana    = command[0]
        cmd_postproc = command[1]
    else:
        cmd_snana = command
        cmd_postproc = ""

    istat0 = 0 ; istat1 = 0

    cmd_plus_log = f"{cmd_snana} > {LOG_FILE}"
    print(f"\n Run command: \n  {cmd_snana}")
    istat0 = os.system(cmd_plus_log)

    if cmd_postproc:
        print(f"\n Run command: \n  {cmd_postproc}")
        istat1 = os.system(cmd_postproc)

    if istat0 == 0 :
        print(f"\n Quick command SUCCESS; check output.")
    else:
        print(f"\n Quick command FAILED; tail -50 {LOG_FILE}")

    if ntail > 0:
        cmd_tail = f"tail -{ntail} {LOG_FILE}"
        os.system(cmd_tail)

    return istat0, istat1
    # end exec_command

def translate_simgen_dump_file(args):
    # Aug 27 2021
    # translate simgen dump file (form SNANA sim) into an ideal 
    # fitres file where fitted parameters are true values.
    # Enables using SALT2mu on true values.

    simgen_dump_file = args.simgen_dump_file
    out_table_file   = "simgen_dump.fitres"

    # define artificially small errors
    zHDERR = 0.0001
    mBERR  = 0.0001
    cERR   = 0.0001/3.0
    x1ERR  = 0.0001/0.14
    IDSURVEY = 1  # anything in SURVEY.DEF file

    df  = pd.read_csv(simgen_dump_file, comment="#", delim_whitespace=True)
    df["CID"] = df["CID"].astype(str)

    FOUND_SIM_mB  = ("SIM_mB" in df)
    FOUND_SIM_x0  = ("SIM_x0" in df)
    FOUND_MAGSMEAR_COH = ("MAGSMEAR_COH" in df)
    FOUND_gammaDM      = ("SALT2gammaDM" in df)

    VARNAMES_STRING = f"CID IDSURVEY zHD zHDERR mB mBERR " \
                      f"x0 x0ERR x1 x1ERR c cERR  COVx0x COVx0c COVx1c"


    nrow = 0 
    with open(out_table_file,"wt") as o:
        o.write(f"VARNAMES: {VARNAMES_STRING} \n")
        for index, row in df.iterrows():
            nrow += 1
            line = "SN: "
            line += f"{row['CID']:14s} "

            line += f"{IDSURVEY} "

            line += f"{row['ZCMB']:.4f} "
            line += f"{zHDERR:0.4f} "

            magsmear_coh=0.0; gammaDM=0.0
            if FOUND_MAGSMEAR_COH:               
                magsmear_coh = row['MAGSMEAR_COH']
                if magsmear_coh == 0.0 : continue # PEAKMJD far from cadence

            if FOUND_gammaDM:
                gammaDM = row['SALT2gammaDM']

            SIM_mB = row['SIM_mB'] 
            if FOUND_SIM_x0 :
                SIM_x0 = row['SIM_x0'] 
            else:
                # mB = -2.5*log10(x0) + 10.63
                SIM_x0 = 10**(-0.4*(SIM_mB-10.63))

            extra_mB = magsmear_coh + gammaDM
            SIM_mB += extra_mB
            SIM_x0 *= 10**(-0.4*extra_mB)
            x0ERR   = SIM_x0 * mBERR

            line += f"{SIM_mB:.5f} "
            line += f"{mBERR:0.5f} "

            line += f"{SIM_x0:11.4e} "
            line += f"{x0ERR:11.4e} "

            line += f"{row['SIM_x1']:.5f} "
            line += f"{x1ERR:0.5f} "

            line += f"{row['SIM_c']:.5f} "
            line += f"{cERR:0.5f} "

            line += f"0.0 0.0 0.0" # 3 covariances

            o.write(f"{line}\n")

    # end translate_simgen_dump_file

def print_HELP():
    see_me = (f" !!! ************************************************ !!!")
    print(f"\n{see_me}\n{see_me}\n{see_me}")
    print(f"{HELP_COMMANDS}")
    sys.exit(' Scroll up to see full HELP menu.\n Done: exiting Main.')
    
# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    args = get_args()
    
    # option for long HELP menus
    if args.HELP : 
        print_HELP()

    if args.version_reformat_fits :
        reformat_fits(args)

    if args.cidlist_text :
        extract_text_format(args)

    if args.simlib_outfile :
        make_simlib(args)

    if args.sntable_list :
        make_sntable(args)

    if args.cidlist_ztable :
        make_ztable(args)

    if args.get_info_phot :
        get_info_photometry(args)

    if args.get_info_code :
        get_info_code(args)

    if args.simgen_dump_file :
        translate_simgen_dump_file(args)

    # END
