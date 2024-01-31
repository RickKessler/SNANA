#!/usr/bin/env python
#
# Interface to numerous snana.exe commands that are hard to find/remember.
# Use "quick_commands.py -H" to get explicit examples on how to combine 
# arguments to perform specific tasks.
#
# Jan 22 2022: add --diff_fitres option 
# Apr 22 2022: for -d option, include 'MU' if it exists
# Sep 12 2022: fix --extract_sim_input for sims run in batch mode
# Oct 13 2022: add : --cov_file option. Rewrites diagonal and
#                    begining of row positions and indices in comments
#                       A. Mitra
#
# Mar 02 2023: if output version name exceeds MXCHAR_VERSION=72, 
#              truncate output version name to fit in max allowed string len
#              for snana.exe.
#
# Jun 26 2022: add --diff_data option 
#
# Jul 19 2023: --extract_spectra_format
# Jan 31 2024: fix bug finding DOCANA keys for --extract_sim_input  option
#
# =========================

import os, sys, argparse, subprocess, yaml, tarfile, fnmatch, glob
import pandas as pd
import numpy as np
import gzip

# ----------------
snana_program          = "snana.exe"
combine_fitres_program = "combine_fitres.exe"
combine_fitres_file    = "combine_fitres.text"

LOG_FILE = "quick_command.log"
MXCHAR_VERSION = 72 # should match same parameter in snana.car

ISTYPE_DIFF_DICT = {
    "FITRES" : "file",
    "DATA"  : "folder"
}

SPECTRA_FORMAT_ASTRODASH = "astrodash"
SPECTRA_FORMAT_SNID      = "snid"

SNDATA_ROOT = os.getenv('SNDATA_ROOT')


HELP_COMMANDS = f"""
# translate TEXT format (from $SNDATA_ROOT/lcmerge) to FITS format
 quick_commands.py -v CFA3_KEPLERCAM --version_reformat_fits CFA3_FITS

# translate TEXT format from local dir to FITS format:
  quick_commands.py -v DES-SN3YR_LOWZ  -p ./  \\
  --version_reformat_fits LOWZ_FITS

# extract a few events from FITS format back to TEXT format:
 quick_commands.py -v SDSS_allCandidates+BOSS --cidlist_text 1032,5944
 quick_commands.py -v SDSS_allCandidates+BOSS --cidlist_text TEN
 quick_commands.py -v SDSS_allCandidates+BOSS --cidlist_text TEN --nospectra

# create SIMLIB file from LOWZ data
 quick_commands.py -v DES-SN3YR_LOWZ \\
   -k $SNDATA_ROOT/kcor/DES/DES-SN3YR/kcor_LOWZ.fits \\
   --simlib_outfile DES-SN3YR_LOWZ.SIMLIB

# make human-readable table of redshifts and vpec for list of CIDs
 quick_commands.py -v DES-SN3YR_LOWZ --cidlist_ztable 2002de,2007cq,2010dt
 quick_commands.py -v DES-SN3YR_LOWZ --cidlist_ztable TEN

# make SNANA-formatted table(s)
 quick_commands.py -v DES-SN3YR_LOWZ --sntable_list SNANA
 quick_commands.py -v DES-SN3YR_LOWZ --sntable_list 'SNANA(text:host)'
 quick_commands.py -v DESY1_forcePhoto_fake_snana_fits -p $DES_ROOT/lcmerge \\
    --sntable_list 'SNANA OUTLIER(nsig:10)'

# extract info about a photometry version
  quick_commands.py -v SDSS_allCandidates+BOSS --get_info_phot

# extract info about SNANA code
  quick_commands.py --get_info_code

# extract spectra from SALT3-K21-Frag (astrodash format with lam Flam Flamerr):
  quick_commands.py --extract_spectra_format astrodash \\
                    -v $SNDATA_ROOT/lcmerge/SALT3TRAIN_K21-Frag/SAL\*

# analyze stat differences (z,mB,x1,c) between two fitres files
# run on same events (e.g., to validate updated data set);
# Computes mean diff, RMS(diff), max outliers ...
  quick_commands.py --diff_fitres  lcfit_ref.fitres  lcfit_test.fitres

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

    msg = "suppress reading/reformattig spectra"
    parser.add_argument("--nospectra", help=msg, action="store_true")

    msg = "TYPE list to extract into text format (e.g., 118,119,120)"
    parser.add_argument("--typelist_text",help=msg,type=str,default="")

    msg = "extra snana.exe args for reformatting (e.g., OPT_MWEBV, ...)"
    parser.add_argument("--reformat_args", help=msg, type=str,default="")

    msg = "CID list to make table of redshifts and vpec"
    parser.add_argument("--cidlist_ztable", help=msg, type=str, default="")

    msg = "Name of covariance matrix to display for debug"
    parser.add_argument("--cov_file", help=msg, type=str, default="")
    
    msg = "get info for photometry version"
    parser.add_argument("--get_info_phot", help=msg, action="store_true")

    msg = "get info for code version"
    parser.add_argument("--get_info_code", help=msg, action="store_true")

    msg = "convert this simgen_dump file to fitres file for SALT2mu"
    parser.add_argument("--simgen_dump_file", help=msg, type=str, default="")

    msg = "extract sim-input file from sim VERSION.README (-v arg)"
    parser.add_argument("--extract_sim_input", help=msg, action="store_true")

    msg = "format/code to extract spectra from data VERSION(s); astrodash, ..."
    parser.add_argument("--extract_spectra_format", help=msg, type=str, default=None)

    msg = "two fitres files to analyse stat difference for SALT2 fit params"
    parser.add_argument("-d", "--diff_fitres", nargs='+', 
                        help=msg, type=str,default=None)

    msg = "two data folders (full paths) to analyse stat difference of contents"
    parser.add_argument("-D", "--diff_data", nargs='+', 
                        help=msg, type=str,default=None)

# SIMLIB_OUT ...

    if len(sys.argv) == 1:  
        parser.print_help(); sys.exit()

    args = parser.parse_args()

    # - - - - - - - - - - - - - - - 
    # if version includes full path, separate here into version
    # and private path args
    if args.version:
        v = args.version
        if '/' in v:
            args.version = os.path.basename(v)
            args.path    = os.path.dirname(v)

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
    typelist    = ""
    prefix_temp = "OUT_TEMP_SNANA"
    command = snana_command_plus_version(args)
    command += f"SNTABLE_LIST SNANA "
    command += f"TEXTFILE_PREFIX {prefix_temp} "
    command += arg_cidlist(cidlist,typelist)

    # extract info from SNANA table and present in human-readable form
    table_script = "get_fitres_values.py"
    table_file   = f"{prefix_temp}.SNANA.TEXT"
    zvar_list    = "zHEL,zHELERR,zCMB,zCMBERR,zHD,zHDERR,VPEC,VPECERR"
    cmd_table    = f"{table_script} -f {table_file} --nrow 10 -v {zvar_list}"

    command_list = [ command, cmd_table ]

    istat0, istat1 = exec_command(command_list,args,0)

    # end make_ztable

def extract_text_format(args):

    vin      = args.version
    vout     = create_vout_string(vin,"TEXT")
    cidlist  = args.cidlist_text
    typelist = args.typelist_text
    reformat_args = args.reformat_args  # Nov 2023
    rmdir_check(vout)
        
    print(f"\n Create new data folder: {vout}")

    command = snana_command_plus_version(args)
    command += f"VERSION_REFORMAT_TEXT {vout} "
    command += arg_cidlist(cidlist,typelist)
    command += f"{reformat_args} "

    if args.nospectra :
        command += ' DEBUG_FLAG -333'

    exec_command(command,args,0)
    # end extract_text_format

def create_vout_string(vin,suffix):
    # Created Mar 2 2023
    # Default is that vout = [vin]_[suffix] ;
    # however, if vout exceeds MXCHAR_VERSION length, then
    # truncate chars from vin so that vout strlen < MXCHAR_VERSION

    vout = f"{vin}_{suffix}"
    len_vout = len(vout)
    if len_vout > MXCHAR_VERSION:
        nchar_remove = len_vout - MXCHAR_VERSION + 1
        vout_orig    = vout
        vin_truncate = vin[:-nchar_remove]
        vout         = f"{vin_truncate}_{suffix}" 
        print(f"\n WARNING:")
        print(f" len({vout_orig}) = {len_vout} \n" \
              f"\t exceeds MXCHAR_VERSION={MXCHAR_VERSION}")
        print(f" Truncate vout -> \n\t {vout}\n")

    return vout
    # end create_vout_string

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

def arg_cidlist(cidlist,typelist):
    # return snana.exe args for cidlist
    arg_string = ''

    if cidlist == "TEN" :
        arg_string = f"MXEVT_PROCESS 10 "
    elif len(cidlist) > 0 :
        arg_string = f"SNCCID_LIST {cidlist} "
        
    if len(typelist) > 0:
        arg_string += f"SNTYPE_LIST {typelist} "

    return arg_string

    # end arg_cidlist

def get_info_photometry(args):

    # if private_data_path is set, then glue it back to version
    if args.path :
        version = f"{args.path}/{args.version}"
    else:
        version = args.version 

    command  = f"{snana_program} GETINFO {version} "
    exec_command(command,args,9)

    # end get_info_photometry

def get_info_code(args):
    command = f"{snana_program} --snana_version"
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

    FOUND_SIM_mB       = ("SIM_mB" in df)
    FOUND_SIM_x0       = ("SIM_x0" in df)
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

def get_README_contents(version):
    # find README file(s)
    # if sim was run interactivly, this(or these) should just be VERSION.README
    # if sim was run with submit_batch, this(or these) should be in misc.tar.gz
    
    # find directory with sim data using snana uti;
    command  = f"{snana_program} GETINFO {version}"
    ret = subprocess.run( [ command ], cwd=os.getcwd(),
                          shell=True, capture_output=True, text=True )

    ret_list = (ret.stdout).split()
    j        = ret_list.index("SNDATA_PATH:")
    path_simdata = ret_list[j+1]
    print(f"\n Found sim data in : {path_simdata} ")

    dict_yaml={}
    readme_file = f"{path_simdata}/{version}.README"
    misc_file = f"{path_simdata}/misc.tar.gz"
    if os.path.exists(misc_file):
        members = tarfile.open(misc_file).getmembers()
        members_name = [m.name for m in members]
#        print (f"len members {len(members)}, len members_name {len(members_name)}")
        wildcard = f"*.README"
        match_list = sorted(fnmatch.filter(members_name, wildcard))
        for imodel in range(0,10):
            key_model = f"MODEL{imodel}"
            model_match = [match for match in match_list if key_model in match]
            if len(model_match)==0: break
            member = [m for m in members if m.name==model_match[0]][0]
#            print (f"member {member.name}, model_match[0] {model_match[0]}")
            with tarfile.open(misc_file).extractfile(member) as r:
                docana_yaml = yaml.safe_load(r)
                dict_yaml[f"{version}_{key_model}"] = docana_yaml
    else:
        # read README created by simulation
        with open(readme_file, 'rt') as r:
            docana_yaml = yaml.safe_load(r)
            dict_yaml[version] = docana_yaml

    return dict_yaml

def write_sim_input_file(sim_input_file, sim_readme_yaml, version_repeat):
    nkey_write = 0

    key_list = [ 'INPUT_KEYS', 'INPUT_KEYS_SNIaMODEL0', 'INPUT_KEYS_NONIaMODEL0' ]
    # xxx makr key = 'INPUT_KEYS'  # INPUT_KEYS_SNIaMODEL0  INPUT_KEYS_NONIaMODEL0

    DOCANA = sim_readme_yaml['DOCUMENTATION']
    INPUT_KEYS = None

    for key in key_list:
        if key in DOCANA:
            INPUT_KEYS = DOCANA[key]

    if INPUT_KEYS is None:
        msgerr = f"\n ERROR: could not find {key_list} key in README "
        assert False, msgerr 

    # make list of special keys that will cause sim-abort or other problems.
    KEYLIST_SUPPRESS = [ 'SIMGEN_DUMPADD' ]

    with open(sim_input_file,"wt") as f :
        for key, val_orig in INPUT_KEYS.items() :

            SKIP_KEY = False
            for key_tmp in KEYLIST_SUPPRESS:
                if key == key_tmp: SKIP_KEY = True
            if SKIP_KEY: continue

            key_plus_colon = key + ':'
            val_out = val_orig
            if key == "GENVERSION" : val_out = version_repeat

            # check for key list with multiple entries 
            if isinstance(val_out,list):
                val_list = val_out  # val_out is already a list
            else:
                val_list = [ val_out ] # convert item into list

            for val in val_list:
                f.write(f"{key_plus_colon:<28}  {val}\n")
                nkey_write += 1
    print(f"\t Done writing {nkey_write} keys for {sim_input_file}")
    return None

def extract_sim_input_file(args):

    # read INPUT_KEYS from DOCUMENTATION in VERSION.README,
    # and create a sim-input file. Modify the GENVERSION
    # to be {version}_REPEAT to avoid clobbering orginal output.
    # Command to implement this feature:
    #   quick_commands.py -v [version] --extract_sim_input                             

    version_orig   = args.version
    version_repeat = f"{version_orig}_REPEAT"

    sim_input_file = f"sim_input_{version_orig}.input"

    dict_yaml = get_README_contents(args.version)
    for model_name,sim_readme_yaml in dict_yaml.items():
        sim_input_file = f"sim_input_{model_name}.input"
        print(f" Create sim-input file: {sim_input_file}")       
        write_sim_input_file(sim_input_file, sim_readme_yaml, version_repeat)

    return

    # end extract_sim_input_file

def extract_spectra(args):

    fmt = args.extract_spectra_format


    if args.path :
        input_data_path = f"{args.path}/{args.version}"
    else:
        input_data_path = f"{SNDATA_ROOT}/lcmerage/{args.version}"


    data_path_list = glob.glob(input_data_path)

    if fmt == SPECTRA_FORMAT_ASTRODASH:
#        cmd_fmt = f"OPT_REFORMAT_SPECTRA 1  " # lam flam
        cmd_fmt = f"OPT_REFORMAT_SPECTRA 2  "  # lam flux
    else:
        msg = f" '{fmt}' is invalid format for --extract_spectra_format argument"
        assert False, msg

    for data_path in data_path_list:
        version           = os.path.basename(data_path)
        private_data_path = os.path.dirname(data_path)
        log_file          = f"EXTRACT_SPECTRA_{version}.log"
        cmd_snana         = f"{snana_program} NOFILE VERSION_PHOTOMETRY {version}  "
        if len(private_data_path) > 2:
            cmd_snana += f"PRIVATE_DATA_PATH {private_data_path}  "

        cmd_snana += cmd_fmt
        cmd_snana += f"> {log_file}"

        print(f" Extract spectra from data version {version} ... ")
        os.system(cmd_snana)



    return
    # end extract_spectra

def util_analyze_diff_INIT(diff_list, WHAT):
    
    # Generic init util for --diff_fitres or --diff_data
    # inputs:
    #   diff_list = list of fitres files, or list of data folders
    #   WHAT      = "FITRES" or "DATA"

    # suppress strange pandas warnings
    pd.options.mode.chained_assignment = None 

    # local variables with names of fitres files
    f_ref  = os.path.expandvars(diff_list[0])
    f_test = os.path.expandvars(diff_list[1])

    istype = ISTYPE_DIFF_DICT[WHAT]  # "file" or "folder"
    what_type = f"{WHAT} {istype}"

    if not os.path.exists(f_ref):
        msgerr = f"Could not find REF-input {what_type}: \n\t {f_ref}"
        assert False, msgerr

    if not os.path.exists(f_test):
        msgerr = f"Could not find TEST-input {what_type}: \n\t {f_test}"
        assert False, msgerr

    print(f"\n Analyze statistical differences between")
    print(f"   REF  {what_type}: {f_ref}")
    print(f"   TEST {what_type}: {f_test}")
    print(f"\t Definition: dif_X = X(TEST) - X(REF)")
    print(f"")

    diff_list_expand = [ f_ref, f_test]
    return diff_list_expand

    sys.stdout.flush()
    
    # end util_analyze_diff_INIT

def util_analyze_diff_EXEC(diff_list, var_list_require, var_list_optional):

    # Inputs
    #   diff_list : list of two FITRES-formatted table files (REF and TEST)
    #   var_list_require : required list of variables to check
    #   var_list_optional : optional list of variables to check (if they exist)

    ff_ref  = diff_list[0]
    ff_test = diff_list[1]

    print(f"\n  Execute analyze_diff on:")
    print(f"   Required var_list = {var_list_require}")
    print(f"   Optional var_list = {var_list_optional}")
    sys.stdout.flush()

    cmd = f"{combine_fitres_program} "
    cmd += f"{ff_ref} {ff_test} "
    cmd += f"t "    # only text output; no HBOOK or ROOT
    

    ret = subprocess.run( [ cmd ], cwd=os.getcwd(),
                          shell=True, capture_output=True, text=True )

    if not os.path.exists(combine_fitres_file):
        msgerr = f"Could not find combined fitres file: {combine_fitres_file}"
        assert False, msgerr

    df  = pd.read_csv(combine_fitres_file, comment="#", delim_whitespace=True)
    df["CID"] = df["CID"].astype(str)

    # define ref variables to check; test var name is {var}_2
    var_check_list = var_list_require

    var_list_missing = []
    for var_name in var_list_require:
        if var_name not in df:
            print(f"ERROR: required varname = {var_name} is not in table.")
            var_list_missing.append(var_name)
        
    if len(var_list_missing) > 0:
        msgerr = f"Missing varnames in table {ff_ref}: \n  {var_list_missing}"
        assert False, msgerr

    # tack on optional varialbes
    for var_name in var_list_optional:
        if var_name in df:
            var_check_list.append(var_name)

    # define dfsel = table rows where both ref and test are defined
    ISTABLE_HEAD = False; ISTABLE_PHOT = False
    if 'zHD_2' in df:
        # regular FITRES file with one row per SN (header info)
        dfsel        = df.loc[df['zHD_2']>-8.0]
        dfcut        = df.loc[df['zHD_2']<-8.0]
        ISTABLE_HEAD = True
    elif 'MJD' in df:
        # LCPLOT file with one row per observation (phot info)
        ISTABLE_PHOT = True
        dfsel        = df.loc[df['MJD_2']>-8.0]
        dfcut        = df.loc[df['MJD_2']<-8.0]        


    len_tot = len(df)
    len_sel = len(dfsel)

    if ISTABLE_HEAD:

        CID_lost_list = dfcut['CID'].to_numpy()

        print(f" TEST table contains {len_sel} of {len_tot} REF events ")
        print(f" CIDs missing in TEST: {CID_lost_list[0:10]}")

    # - - - - - 
    print("")
    print("   quantity          avg       median      std          " \
          f"min/max      CIDmin/CIDmax")
    print("# --------------------------------------------------" \
          "--------------------------- ")
    for var in var_check_list:
        var_2   = f"{var}_2"
        var_dif = f"dif_{var}"
        dfsel[var_dif] = dfsel[var_2] - dfsel[var]
        mean  = dfsel[var_dif].mean()
        med   = dfsel[var_dif].median()
        std = 0.0
        if len_sel > 1: std   = dfsel[var_dif].std()
        mn    = dfsel[var_dif].min()
        mx    = dfsel[var_dif].max()

        CIDmin = None ; CIDmax=None
        if mn < 0.0 :
            CIDmin = dfsel.loc[dfsel[var_dif].idxmin()]['CID']
            if ISTABLE_PHOT:
                MJDmin = dfsel.loc[dfsel[var_dif].idxmin()]['MJD']   
                CIDmin += f"({MJDmin})"
        if mx > 0.0 :
            CIDmax = dfsel.loc[dfsel[var_dif].idxmax()]['CID']
            if ISTABLE_PHOT:
                MJDmax = dfsel.loc[dfsel[var_dif].idxmax()]['MJD']   
                CIDmax += f"({MJDmax})"

        print(f"  {var_dif:16} {mean:8.5f}  {med:8.5f}   {std:8.5f}  " \
              f" {mn:8.5f}/{mx:8.5f}  {CIDmin}/{CIDmax}")

        sys.stdout.flush()        
    
    # end util_analyze_diff_EXEC
    
def analyze_diff_data(args):

    folder_expand_list = util_analyze_diff_INIT(args.diff_data,"DATA")
    
    # construct snana command to produce SNANA+LCPLOT table for
    # each data foloder

    outname_list = [ 'REF', 'TEST' ]
    snana_table_list = []
    phot_table_list  = []
    snana_log_list   = []
    out_prefix_list  = []

    # get list of filters to show photometry results by band
    cmd = f"{snana_program} GETINFO {folder_expand_list[0]}"
    ret = subprocess.run( [ cmd ], cwd=os.getcwd(),
                          shell=True, capture_output=True, text=True )
    ret_list = (ret.stdout).split()
    j = ret_list.index('FILTERS:')
    filter_string = ret_list[j+1]
    filter_list   = [ *filter_string ]

    for data_folder, outname in zip(folder_expand_list,outname_list):

        data_dir         = os.path.dirname(data_folder)
        version          = os.path.basename(data_folder)
        out_prefix       = f"OUT_TABLE_{version}_{outname}"
        out_prefix_list.append(out_prefix)

        snana_log_file   = f"{out_prefix}.LOG"
        snana_table_file = f"{out_prefix}.SNANA.TEXT"
        phot_table_file  = f"{out_prefix}.LCPLOT.TEXT"
        snana_table_list.append(snana_table_file)
        phot_table_list.append(phot_table_file)
        snana_log_list.append(snana_log_file)

        # remove existing table file to avoid process stale file if
        # snana job fails.
        if os.path.exists(snana_table_file):
            os.remove(snana_table_file)
        if os.path.exists(phot_table_file):
            os.remove(phot_table_file)

        print(f"  Extract SNANA tables for {version} ... ")
        sys.stdout.flush()

        cmd  = f"{snana_program} NOFILE " \
               f"PRIVATE_DATA_PATH {data_dir} " \
               f"VERSION_PHOTOMETRY {version} " \
               f"SNTABLE_LIST 'SNANA(text:key,text:host) LCPLOT(text:key)' "\
               f"TEXTFILE_PREFIX {out_prefix} "
        
        arg_select =  arg_cidlist(args.cidlist_text, args.typelist_text)
        if len(arg_select) > 0 :
            cmd += f"{arg_select} "
        #sys.exit(f"\n arg_select = {arg_select}  [cidlist_text={args.cidlist_text}")

        cmd += f" > {snana_log_file} "

        ret = subprocess.run( [ cmd ], cwd=os.getcwd(),
                              shell=True, capture_output=True, text=True )
        
        if not os.path.exists(snana_table_file):
            msgerr = f"Could not find expected output table {snana_table_file}\n" \
                     f"Check {snana_log_file}" 
            assert False, msgerr

    # - - - - - - - - - - - - - - - 
    # analyze snana table
    var_list_require  = [ 'zHD', 'zHDERR', 'VPEC', 'VPECERR', 'MWEBV',
                          'SNRMAX1', 'SNRMAX3',
                          'HOST_ZPHOT', 'HOST_DDLR' 
    ]  
    var_list_optional = [ ]
    util_analyze_diff_EXEC(snana_table_list, var_list_require, var_list_optional)
    
    # analyze LCPLOT(PHOT) table for each band
    var_list_require = [ 'MJD', 'FLUXCAL', 'FLUXCAL_ERR' ]
    var_list_optional = [ ]
    util_analyze_diff_EXEC(phot_table_list, var_list_require, var_list_optional)


    # - - - - - -  - -
    # clean up mess
    CLEAN = False
    if CLEAN:
        print(f" Remove temporary snana output files.")
        os.remove(combine_fitres_file)
        for out_prefix in out_prefix_list:
            os.system(f"rm {out_prefix}*")

    # end analyze_diff_data

def analyze_diff_fitres(args):

    diff_fitres_expand = util_analyze_diff_INIT(args.diff_fitres,"FITRES")

    # remove pre-existing combine-fitres file to avoid confusion
    # if new file fails to be created.
    if os.path.exists(combine_fitres_file):
        os.remove(combine_fitres_file)

    var_list_require  = [ 'zHD', 'PKMJD', 'mB', 'x1', 'c' ]  
    var_list_optional = [ 'MU', 'MUERR' ]

    util_analyze_diff_EXEC(diff_fitres_expand, var_list_require, var_list_optional)

    # end analyze_diff_fitres

def analyze_diff_fitres_legacy(args):

    # suppress strange pandas warnings
    pd.options.mode.chained_assignment = None 

    # local variables with names of fitres files
    ff_ref  = os.path.expandvars(args.diff_fitres[0])
    ff_test = os.path.expandvars(args.diff_fitres[1])

    if not os.path.exists(ff_ref):
        msgerr = f"Could not find REF-input fitres file: \n\t {ff_ref}"
        assert False, msgerr

    if not os.path.exists(ff_test):
        msgerr = f"Could not find TEST-input fitres file: \n\t {ff_test}"
        assert False, msgerr

    print(f"\n Analyze statistical differences between")
    print(f"\t REF  fitres file: {ff_ref}")
    print(f"\t TEST fitres file: {ff_test}")
    print(f"\t Definition: dif_X = X(TEST) - X(REF)")
    sys.stdout.flush()

    # remove pre-existing combine-fitres file to avoid confusion
    # if new file fails to be created.
    if os.path.exists(combine_fitres_file):
        os.remove(combine_fitres_file)

    cmd = f"{combine_fitres_program} "
    cmd += f"{ff_ref} {ff_test} "
    cmd += f"t "    # only text output; no HBOOK or ROOT

    ret = subprocess.run( [ cmd ], cwd=os.getcwd(),
                          shell=True, capture_output=True, text=True )

    if not os.path.exists(combine_fitres_file):
        msgerr = f"Could not find combined fitres file: {combine_fitres_file}"
        assert False, msgerr

    df  = pd.read_csv(combine_fitres_file, comment="#", delim_whitespace=True)
    df["CID"] = df["CID"].astype(str)

    # define ref variables to check; test var name is {var}_2
    var_check_list = [ 'zHD', 'PKMJD', 'mB', 'x1', 'c' ]  

    # Apr 22 2022: add MU if it's there (e.g., output of BBC)
    if 'MU' in df:
        var_check_list.append('MU')
        var_check_list.append('MUERR')

    # define dfsel = table rows where both ref and test are defined
    #dfsel        = df.loc[df['c_2']>-8.0]
    dfsel        = df.loc[df['c_2']>-8.0]
    dfcut        = df.loc[df['c_2']<-8.0]

    len_tot = len(df)
    len_sel = len(dfsel)
    CID_lost_list = dfcut['CID'].to_numpy()

    print(f" TEST table contains {len_sel} of {len_tot} REF events ")
    print(f" CIDs missing in TEST: {CID_lost_list[0:10]}")

    print("")
    print("   quantity          avg       median      std          " \
          f"min/max      CIDmin/CIDmax")
    print("# --------------------------------------------------" \
          "--------------------------- ")
    for var in var_check_list:
        var_2   = f"{var}_2"
        var_dif = f"dif_{var}"
        dfsel[var_dif] = dfsel[var_2] - dfsel[var]
        mean  = dfsel[var_dif].mean()
        med   = dfsel[var_dif].median()
        std   = dfsel[var_dif].std()
        mn    = dfsel[var_dif].min()
        mx    = dfsel[var_dif].max()

        CIDmin = None ; CIDmax=None
        if mn < 0.0 :
            CIDmin = dfsel.loc[dfsel[var_dif].idxmin()]['CID']
        if mx > 0.0 :
            CIDmax = dfsel.loc[dfsel[var_dif].idxmax()]['CID']

        print(f"  {var_dif:16} {mean:8.5f}  {med:8.5f}   {std:8.5f}  " \
              f" {mn:8.5f}/{mx:8.5f}  {CIDmin}/{CIDmax}")
        sys.stdout.flush()

    return
    # end analyze_diff_fitres_legacy

def rewrite_cov_file(args):

    # Created 13 Oct 2022 by A.Mitra
    # 1. Rewrite cov with row,column labels
    # 2. Add "Start row" for readibility 0.20043.  # (0,2)  START_ROW
    # 3. Add Diagonal" for readibility : 0.23243.  # (2,2)  DIAGONAL

    cov_file = os.path.expandvars(args.cov_file)
    data     = args.cov_file

    X=[];comment_1 = []; comment_2=[];
    com_row = 'START ROW';  com_d = ' DIAGONAL'; com_null=' '
    cc = 0;index_elements = [];

    cov_basename  = os.path.basename(cov_file) 
    out_cov_file  = f"DISPLAY_{cov_basename}"
    print(f"Input cov matrix file: {cov_file}")
    print(f"rewrite cov matrix to: {out_cov_file}")

    c     = pd.read_csv(data,compression='gzip',sep='\s+',comment="#")
    c0    = np.array(c); 
    shape = int(np.sqrt(np.shape(c)[0]))
    c1    = np.reshape(c0,(-1,shape));
    D     = np.diag(c1);

    for i in np.ndindex(c1.shape):
        tmp = "#" + str(i)
        X.append(tmp)
        index_elements.append((c1[i],i))

        if(cc%shape == 0):
            comment_1.append(com_row)
        else :
            comment_1.append(com_null)

        if (D.__contains__(c1[i])==True):     
            comment_2.append(com_d)
        else :
            comment_2.append(com_null)
            cc += 1

    X = pd.DataFrame(X) 
    X.columns = (["Index"])

    comment_1= pd.DataFrame(comment_1) 
    comment_1.columns = (["Comments"])

    comment_2= pd.DataFrame(comment_2) 
    comment_2.columns = (["Comments"])

    comments = comment_1 + comment_2

    cov_m = pd.concat([pd.DataFrame(c),X], axis=1)
    cov_m = pd.concat([cov_m,comments],    axis=1)
    #print(cov_m)

    cov_m.to_csv(out_cov_file, sep='\t', encoding='utf-8',
                 header=True, index=False, compression='gzip')

    return 
    # end rewrite_cov_file
   


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

    if args.cov_file :
        rewrite_cov_file(args)

    if args.get_info_phot :
        get_info_photometry(args)    
        
    if args.get_info_code :
        get_info_code(args)

    if args.simgen_dump_file :
        translate_simgen_dump_file(args)

    if args.extract_sim_input :
        extract_sim_input_file(args)

    if args.extract_spectra_format :
        extract_spectra(args)

    if args.diff_fitres :
        #analyze_diff_fitres_legacy(args)
        analyze_diff_fitres(args)

    if args.diff_data :
        analyze_diff_data(args)

    # END main

