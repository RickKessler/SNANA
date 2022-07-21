# Created Oct 2021 [move utils out of base]
# utilities to write data in snana format
#
# Jun 24 2022 RK - avoid writing HOSTGAL2_MAG[ERR] of there is no 2nd host.
#
import datetime
import glob
import logging
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import yaml

import makeDataFiles_params as gpar
import makeDataFiles_util as util

#from makeDataFiles_params import *

# ====================================================

def write_event_text_snana(args, config_data, data_event_dict):

    # Inputs:
    #   args : user command line inputs
    #   config_data       : info about data units and phot varnames
    #   data_event_dict   : current event: header, phot, spec

    # create one text file and write one event described by
    # dictionary data_event_dict.
    # Input data_unit_name is used to determine name of folder.
    ISTEXT = True
    outdir                = args.outdir_snana
    data_unit_name_list   = config_data['data_unit_name_list']
    data_unit_nevent_list = config_data['data_unit_nevent_list']
    prefix                = config_data['data_folder_prefix']

    data_unit_name    = data_event_dict['data_unit_name']
    index_unit        = data_event_dict['index_unit']

    nevent         = data_unit_nevent_list[index_unit]
    folder         = output_data_folder_name(config_data,
                                             data_unit_name, ISTEXT)
    data_dir       = f"{outdir}/{folder}"
    tar_file       = f"{outdir}/{folder}.tar.gz"
    exist_folder   = os.path.exists(data_dir)
    exist_tar_file = os.path.exists(tar_file)

    if nevent == 0 :
        util.create_output_folder(data_dir)

        # xxxx mark delete Jul 20 2022 xxxxxxxxx
        # remove folder if it already exists
        #if exist_folder :
        #    cmd_rm = f"rm -r {data_dir}"
        #    os.system(cmd_rm)
        #if exist_tar_file :
        #    cmd_rm = f"rm -r {tar_file}"
        #    os.system(cmd_rm)
        # create folder
        #logging.info(f"\t Create folder {folder}")
        #sys.stdout.flush()
        #os.mkdir(data_dir)

    # - - - - -
    head_raw  = data_event_dict['head_raw']
    head_calc = data_event_dict['head_calc']

    # check optional pieces of header
    head_private = None
    if 'head_private' in data_event_dict :
        head_private = data_event_dict['head_private']

    head_sim     = None
    if 'head_sim' in data_event_dict :
        head_sim  = data_event_dict['head_sim']

    phot_raw  = data_event_dict['phot_raw']
    spec_raw  = data_event_dict['spec_raw']

    SNID      = head_raw[gpar.DATAKEY_SNID]

    str_SNID     = SNID
    if SNID.isdigit(): str_SNID = f"{int(SNID):010d}"

    data_file     = f"{data_dir}/{prefix}_{str_SNID}.DAT"

    with open(data_file, "wt") as f :

        # write header info
        write_header_snana(f,head_raw)

        f.write("\n# computed quantities \n")

        write_header_snana(f,head_calc)

        if head_private :
            f.write("\n# PRIVATE (non-standard) variables \n")
            write_header_snana(f,head_private)

        if head_sim :
            f.write("\n# sim/truth quantities \n")
            write_header_snana(f,head_sim)

        # write epoch/phot info
        write_phot_snana(f, head_raw, phot_raw, config_data)

        # write optional spectra
        write_spec_snana(f, head_raw, spec_raw)

    # end write_event_text_snana

def write_header_snana(f, data_head):

    # write list of header key,val pairs in data_head.
    # If XXX and XXX_ERR both exist, write as
    #   KEYXXX:  XXX +_ XXX_ERR
    # Be careful for HOSTGAL_XXX keys that are filter dependent;
    # filter-dependent values are written on same line as key,
    # not separate line per filter.
    #
    # Some tricks are used for human-readability;
    # e.g,, -1.00000 is written as -1, etc ...

    # init filter-dependent list of hostgal values
    hostgal_string_vals = {}
    for prefix  in gpar.HOSTKEY_PREFIX_LIST:
        hostgal_string_vals[prefix] = ""

    n_hostkey = 0

    for key in data_head:
        #print(f" xxx header key = {key}")
        if '_ERR' in key: continue
        key_plus_err   = f"{key}_ERR"
        key_plus_colon = f"{key}:"
        val            = data_head[key]

        string_val     = f"{val}"
        if isinstance(val,float):  
            ival = int(val)
            if val == ival : 
                string_val = f"{ival}"
            else:
                string_val = f"{val:.6f}" # Jun 24 2022 RK doesn't work??

        is_redshift    = "REDSHIFT" in key
        is_private     = "PRIVATE"  in key
        is_hostgal     = "HOSTGAL"  in key

        # skip keys with NULL value, except for private keys and 
        # keys with REDSHIFT ... must always write REDSHIFT.
        skip = False
        force_write = is_redshift or is_private
        if not force_write:
            skip = val in gpar.VAL_NULL_LIST

        if skip: continue

        # write blank line before host info (for visual aid)
        if is_hostgal:
            n_hostkey += 1
            if n_hostkey == 1: f.write(f"\n")

        # increment filter-dependent host strings
        skip_hostkey = False
        for prefix  in gpar.HOSTKEY_PREFIX_LIST:
            prefix_ = prefix + "_"
            if prefix_ in key:
                skip_hostkey = True
                hostgal_string_vals[prefix] += f"{val:.3f} "
        if skip_hostkey : continue

        if key_plus_err in data_head:
            err  = data_head[key_plus_err]
            if err > 0.0:
                string_val = f"{val:9.6f} +- {err:9.6f}"
            else:
                string_val = f"{int(val)} +- {int(err)}"  # -9 +_ -9

        f.write(f"{key_plus_colon:<20s}  {string_val} \n")

    # - - - - - -
    # write mag-dependent host info where n_filter values
    # are written on one line.
    for prefix in gpar.HOSTKEY_PREFIX_LIST:
        string_vals = hostgal_string_vals[prefix]
        if len(string_vals) > 0:
            f.write(f"{prefix}: {string_vals}\n")

    f.flush()
    return

    # end write_header_snana

def write_phot_snana(f, head_raw, phot_raw, config_data):

    # write photometry (phot_raw) in SNANA format to text file
    # poitner f.
    nvar_obs      = config_data['nvar_obs']
    varlist_obs   = config_data['varlist_obs']
    varlist_fmt   = config_data['varlist_fmt']
    vallist_undef = config_data['vallist_undef']
    varstring_obs = ' '.join(varlist_obs)
    msgerr   = []
    SNID     = head_raw[gpar.DATAKEY_SNID]
    FILTERS  = head_raw[gpar.DATAKEY_FILTERS]
    NOBS     = phot_raw[gpar.DATAKEY_NOBS]

    f.write(f"\n# -------------------------------------- \n" \
            f"# obs info\n")
    f.write(f"NOBS: {NOBS}\nNVAR: {nvar_obs} \n"
            f"VARLIST: {varstring_obs}\n")

    for obs in range(0,NOBS):
        LINE = "OBS:"
        for varname,fmt,val_undef in \
            zip(varlist_obs,varlist_fmt,vallist_undef):
            val = phot_raw[varname][obs]
            if val == None:  val = val_undef

            if val == 12345.333 :  # gpar.VAL_ABORT :  # problem for DES??
                msgerr.append(f"Missing required PHOT column {varname}")
                msgerr.append(f"Check SNID = {SNID}")
                util.log_assert(False,msgerr)

            if varname == 'BAND' :
                band = val[-1]
                if band not in FILTERS:
                    msgerr.append(f"Unknown band {band} is not in "\
                                  f"{FILTERS} for SNID={SNID}")
                    msgerr.append(f"Check SURVEY_INFO[FILTERS] ")
                    util.log_assert(False,msgerr)

            LINE += f" {val:{fmt}}"
        f.write(f"{LINE}\n")

    # - - - - -
    f.write(f"END:\n")
    return

    # end write_phot_snana


def write_spec_snana(f, head_raw, spec_raw):
    SNID     = head_raw[gpar.DATAKEY_SNID]

    NSPEC = len(spec_raw)
    if NSPEC == 0 : return

    VARLIST = "LAMMIN LAMMAX FLAM FLAMERR"
    NVAR    = len(VARLIST.split())
    f.write(f"\nNSPECTRA: {NSPEC} \n\n");

    f.write(f"NVAR_SPEC:      {NVAR} \n")
    f.write(f"VARNAMES_SPEC:  {VARLIST} \n")

    ID = 0
    for mjd in spec_raw.keys():
        spec_data      = spec_raw[mjd]
        wave_list      = spec_data['wave']
        flux_list      = spec_data['flux']
        fluxerr_list   = spec_data['fluxerr']
        ID    += 1
        nblam  = len(wave_list)

        wave_diff_list  = np.diff(wave_list)
        diff_last       = wave_diff_list[nblam-2]
        wave_diff_list  = np.append(wave_diff_list,diff_last)

        wave_min_list   = wave_list - wave_diff_list/2.0
        wave_max_list   = wave_list + wave_diff_list/2.0

        f.write(f"SPECTRUM_ID:    {ID} \n")
        f.write(f"SPECTRUM_MJD:   {mjd:.3f} \n")
        f.write(f"SPECTRUM_NLAM:  {nblam} \n")

        for wave_min, wave_max, flux,fluxerr in \
            zip(wave_min_list, wave_max_list, flux_list,fluxerr_list):
            f.write(f"SPEC: {wave_min:.2f}  {wave_max:.2f} " \
                    f"{flux:12.3e} {fluxerr:12.3e}\n")

        f.write(f"\n")

    return

    # end write_spec_snana

def output_data_folder_name(config_data, data_unit_name, ISTEXT):

    prefix              = config_data['data_folder_prefix']
    data_unit_name_list = config_data['data_unit_name_list']

    if data_unit_name not in data_unit_name_list:
        msgerr = []
        msgerr.append(f" Invalid data unit '{data_unit_name}")
        msgerr.append(f" Valid data units are : ")
        msgerr.append(f"   {data_unit_name_list}")  # <<< I assume this is correct
        #msgerr.append(f"   {data_unit_list}")      # <<< instead of this
        util.log_assert(False,msgerr)

    folder = f"{data_unit_name}"
    if ISTEXT:
        folder = f"{gpar.FORMAT_TEXT}_{folder}"

    return folder

    # end output_data_folder_name

def write_aux_files_snana(name, args, config_data):

    # write auxilary files (.LIST and .README) for data unit "name"
    # with TEXT formatted files.

    outdir        = args.outdir_snana
    nevent_list   = config_data['data_unit_nevent_list']
    name_list     = config_data['data_unit_name_list']
    prefix        = config_data['data_folder_prefix']
    readme_stats_list = config_data['readme_stats_list']

    folder_out  = output_data_folder_name(config_data, name, True)
    index_unit  = name_list.index(name)

    msg = f" Create aux files for {folder_out} and gzip " \
          f"{gpar.TEXTFILE_SUFFIX} files."
    logging.info(msg)

    data_dir      = f"{outdir}/{folder_out}"
    search_string = f"{prefix}*{gpar.TEXTFILE_SUFFIX}"

    data_file_list = glob.glob1(data_dir, f"{search_string}" )
    list_file      = f"{data_dir}/{folder_out}.LIST"
    readme_file    = f"{data_dir}/{folder_out}.README"

    with open(list_file,"wt") as f:
        for data_file in data_file_list:
            f.write(f"{data_file}\n")

    # prepare standard readme dictionary for global utility write_readme
    readme_dict = {
        'readme_file'  : readme_file,
        'readme_stats' : readme_stats_list[index_unit],
        'data_format'  : gpar.FORMAT_TEXT,
        'docana_flag'  : True
    }
    util.write_readme(args, readme_dict)

    # gzip data files
    cmd = f"cd {data_dir} ; gzip {search_string}"
    os.system(cmd)

    # end write_aux_files_snana

def convert2fits_snana(args, config_data):

    # loop over newly created TEXT file versions and convert
    # to fits format ... then tar up TEXT folder.

    outdir        = args.outdir_snana
    text          = args.text
    nevent_list   = config_data['data_unit_nevent_list']
    name_list     = config_data['data_unit_name_list']
    prefix        = config_data['data_folder_prefix']
    readme_stats_list = config_data['readme_stats_list']

    NEVT_SPECTRA  = config_data['NEVT_SPECTRA']
    write_spectra = False

    opt_snana = gpar.OPTIONS_TEXT2FITS_SNANA
    if NEVT_SPECTRA > 0 :  # global counter over all data units
        opt_snana += f"  {gpar.OPTION_TEXT2FITS_SPECTRA_SNANA}"
        write_spectra = True

    print(f"")
    sys.stdout.flush()

    for nevent, name in zip(nevent_list, name_list):
        if nevent == 0 : continue

        folder_text    = output_data_folder_name(config_data, name, True)
        folder_fits    = output_data_folder_name(config_data, name, False)
        index_unit     = name_list.index(name)
        log_file       = f"{folder_text}/convert2fits_{folder_fits}.log"
        yaml_file      = f"{outdir}/{folder_text}.YAML"  # expected output

        msg = f"  Convert TEXT -> FITS for {folder_fits}" \
              f" NEVT={nevent}  (write spectra: {write_spectra})"
        logging.info(msg)
        sys.stdout.flush()

        time_0 = datetime.datetime.now()
        outdir_text  = f"{outdir}/{folder_text}"
        outdir_fits  = f"{outdir}/{folder_fits}"

        # rm fits folder if still there from previous job
        if os.path.exists(outdir_fits):
            cmd_rm = f"cd {outdir} ; rm -r {folder_fits}"
            os.system(cmd_rm)

        cmd_snana   = f"{gpar.PROGRAM_SNANA} NOFILE " \
                      f"PRIVATE_DATA_PATH ./ " \
                      f"VERSION_PHOTOMETRY    {folder_text} " \
                      f"VERSION_REFORMAT_FITS {folder_fits} " \
                      f"{opt_snana} "
        cmd = f"cd {outdir}; {cmd_snana} > {log_file}"
        os.system(cmd)

        # - - - -
        # if YAML file doesn't exist, abort with message that
        # convert job probably aborted or crashed.
        if not os.path.exists(yaml_file):
            msgerr = []
            msgerr.append(f"Cannot find expected yaml file:")
            msgerr.append(f"    {yaml_file}")
            msgerr.append(f"TEXT->FITS convert job probably " \
                          f"aborted or crashed;")
            msgerr.append(f"See convert-log file:")
            msgerr.append(f"    {outdir}/{log_file} ")
            util.log_assert(False,msgerr)

        # - - - - -
        # clean up

        # xxxxxxxxx Jan 6 mark delete since translate code does gzip
        # gzip FITS files and make compressed tar file from TEXT dir
        ### cmd_gzip_fits = f"cd {outdir_fits} ; gzip *.FITS"
        ### os.system(cmd_gzip_fits)
        # xxxxxxxxxxxxx end mark xxxxxxxxxx

        tar_file = f"{folder_text}.tar"
        cmd_tar_text  = f"cd {outdir} ; " \
                        f"tar -cf {tar_file} {folder_text} ; " \
                        f"gzip {tar_file} ; " \
                        f"rm -r {folder_text} "

        if not text:
            os.system(cmd_tar_text)

        # remove YAML file
        cmd_rm = f"rm {yaml_file}"
        os.system(cmd_rm)

        # re-write readme in FITS data folder
        readme_file    = f"{outdir_fits}/{folder_fits}.README"
        readme_dict = {
            'readme_file'  : readme_file,
            'readme_stats' : readme_stats_list[index_unit],
            'data_format'  : gpar.FORMAT_FITS,
            'docana_flag'  : True
        }
        util.write_readme(args, readme_dict)

        time_1   = datetime.datetime.now()
        time_dif = (time_1 - time_0).total_seconds()
        rate     = int(float(nevent)/float(time_dif))
        logging.info(f"\t Rate(convert+cleanup): {rate}/sec ")
        sys.stdout.flush()

    # - - - - -

    return

    # end convert2fits_snana


# ================================================
#    MERGE PROCESS
# ================================================

def merge_snana_driver(args):

    # called by main after processing to perform optional
    # merge process; e.g., convert text to other format,
    # re-organize files, etc...

    outdir = args.outdir_snana
    survey = args.survey

    if outdir is None:
        sys.exit(f"\n ERROR: must specify --outdir_snana\n")

    print(f"\n Merge {outdir}")

    # merge SPLIT folders; get all prefixes by scooping up all SPLIT001 job

    search_string = f"{survey}*{gpar.PREFIX_SPLIT}001"
    split_dir_list = sorted(glob.glob1(outdir, search_string ))
    for split_dir in split_dir_list:
        # if split_dir = LSST_WFDY01_SPLIT001, base_name=LSST_WFDY01
        merge_folder = split_dir.split(f"_{gpar.PREFIX_SPLIT}")[0]
        search_string = f"{merge_folder}_{gpar.PREFIX_SPLIT}*"
        merge_snana_folders(gpar.MODE_MERGE_MOVE,
                            outdir, search_string, merge_folder)


    # merge Y## folders into folder with all seasons
    search_string = f"{survey}_*{gpar.PREFIX_SEASON}*"
    year_dir_list = sorted(glob.glob1(outdir, search_string ))
    merge_folder  = year_dir_list[0].split(f"_{gpar.PREFIX_SEASON}")[0]
    merge_snana_folders(gpar.MODE_MERGE_LINK,
                        outdir, search_string, merge_folder)

    # archive TEXT versions
    TEXT_archive_dir = "TEXT_archive"
    TEXT_list = glob.glob1(outdir, f"TEXT*" )
    if len(TEXT_list) > 0 :
        os.mkdir(f"{outdir}/{TEXT_archive_dir}")
        cmd_mv = f"cd {outdir} ; mv TEXT*.tar.gz {TEXT_archive_dir}"
        os.system(cmd_mv)

    return
    # end merge_snana_driver

def merge_snana_folders(MODE, outdir, folder_list_string, merge_folder):

    # e.g., folder_list_string = LSST_WFD01_SPLIT*
    #       merge_folder       = LSST_WFD01
    #    ->
    #      combine data from all LSST_WFD01_SPLIT* folders into
    #      single folder LSST_WFD01

    msgerr = []
    logging.info(f"\t Create merge-folder {merge_folder}")
    merge_folder_full = f"{outdir}/{merge_folder}"
    if os.path.exists(merge_folder_full):
        msgerr.append(f"{merge_folder_full} already exists ?!?!?!")
        msgerr.append(f"Something is wacky.")
        sys.exit(msgerr)

    os.mkdir(merge_folder_full)

    n_move = 0
    statsum_dict = {}
    for key in gpar.KEYLIST_README_STATS:   statsum_dict[key] = 0

    folder_list = glob.glob1(outdir, f"{folder_list_string}" )
    for folder in folder_list :
        n_move += 1
        FOLDER = f"{outdir}/{folder}"

        if MODE == gpar.MODE_MERGE_MOVE :
            mv_string = f"*.FITS.gz"
            cmd_mv = f"mv {mv_string} ../{merge_folder}"
            cmd    = f"cd {FOLDER}; {cmd_mv}"
        else:
            cmd = f"cd {outdir}/{merge_folder} ; "
            fits_file_list = glob.glob1(FOLDER, f"*.FITS.gz" )
            for fits_file in fits_file_list :
                cmd += f"ln -s ../{folder}/{fits_file} {fits_file} ; "
        os.system(cmd)

        if n_move == 1 :
            cmd_mv = f"mv *.IGNORE ../{merge_folder}/{merge_folder}.IGNORE"
            cmd    = f"cd {FOLDER}; {cmd_mv}"
            os.system(cmd)

        # increment sum stats from readme
        README_file  = f"{FOLDER}/{folder}.README"
        README_yaml  = util.read_yaml(README_file)
        for key in gpar.KEYLIST_README_STATS:
            NEVT  = README_yaml[gpar.DOCANA_KEY][key]
            statsum_dict[key] += NEVT

    # - - - - - - - -
    # update sum stats and re-write readme
    README_file = f"{merge_folder_full}/{merge_folder}.README"
    for key in gpar.KEYLIST_README_STATS:
        NEVT = statsum_dict[key]
        README_yaml[gpar.DOCANA_KEY][key] = statsum_dict[key]
        util.write_yaml(README_file,README_yaml)

    # - - - -
    # remove original folders
    if MODE == gpar.MODE_MERGE_MOVE:
        cmd_rm = f"rm -r {folder_list_string}"
        cmd    = f"cd {outdir}; {cmd_rm}"
        os.system(cmd)

    # create merged LIST file
    HEAD_list = glob.glob1(merge_folder_full, "*HEAD*.FITS.gz" )
    LIST_file = f"{merge_folder_full}/{merge_folder}.LIST"
    with open(LIST_file,"wt") as l :
        for item in HEAD_list:
            HEAD_file = item
            # strip off .gz extension
            if '.gz' in item: HEAD_file = item.split('.gz')[0]
            l.write(f"{HEAD_file}\n")

    #print(f" xxx ---------------------------- ")
    #print(f" xxx base_name = {base_name} ")
    #print(f" xxx split_dir_list = {split_dir_list} ")

    return
    # end merge_snana_folders

# end:
