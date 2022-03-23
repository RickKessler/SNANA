#!/usr/bin/env python
#
# Created Jan 22 2022 by R.Kessler
# Read data version created by combine_dataVersions.py,
# and break it up into original data versions. The following
# are modified in each output data file:
#
#  + SURVEY: VERSION_COMBINED(VERSION) -> SURVEY: VERSION
#  + restore original filters
#  + output folder name is  UNCOMBINED_[SURVEY]
#
# Syntax:
#    uncombine_dataVersions.py \
#        -v [version] 
#        -p [private_data_path]  # optional
#
# Notes:
#   + This script does NOT implement HEADER overrides nor use snana.exe.
#   + works only on TEXT format and real data
#   + original motivaton is to organize Pantheon+. Need uncombined data files
#      to combine with newer samples.
#
# =========================
import os, sys, shutil, glob, logging, argparse

SNDATA_ROOT        = os.getenv('SNDATA_ROOT')
PATH_DATA_DEFAULT = f"{SNDATA_ROOT}/lcmerge"

PREFIX_VERSION_UNCOMBINED = "UNCOMBINED"

KEY_SURVEY  = "SURVEY:"
KEY_FILTERS = "FILTERS:"
KEY_VARLIST = "VARLIST:"
KEY_NOBS    = "NOBS:"
KEY_OBS     = "OBS:"

COLNAME_BAND_LIST = [ "FLT", "BAND" ]

# ==========================
# ============================
def setup_logging():

    #logging.basicConfig(level=logging.DEBUG,
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)

def get_args():
    parser = argparse.ArgumentParser()

    msg = "name of combined data version to uncombine"
    parser.add_argument("-v", "--version", 
                        help=msg, nargs="?", default=None)

    msg = "optional private data path"
    parser.add_argument("-p", "--path_data",
                        help=msg, nargs="?", default=PATH_DATA_DEFAULT)
    
    args = parser.parse_args()

    if '/' in args.version:
        v = args.version
        args.version      = os.path.basename(v)
        args.path_data    = os.path.dirname(v)
    
    args.folder_version  = f"{args.path_data}/{args.version}"

    return args
    # end get_args

# ==============================
def get_file_list(args):
    v = args.version
    folder = args.folder_version
    list_file = f"{folder}/{v}.LIST"
    
    if not os.path.exists(list_file):
        msg = f"ERROR: cannot find list file:\n{list_file}"
        assert False,  msg

    data_file_list = []
    with  open(list_file,"rt") as l:
        for line in l:
            if line[0] == '#' :  continue
            fnam = line.rstrip('\n')
            fnam = fnam.rstrip(' ')
            data_file_list.append(fnam)

    return data_file_list
    # end get_file_list

def extract_data_file(data_file, args, uncombined_dict ):

    # copy data file to uncompressed version
    # return uncombined_dict = updated list of uncombined versions

    v      = args.version
    folder = args.folder_version
    data_file_path = f"{folder}/{data_file}"

    with  open(data_file_path,"rt")  as f:
        contents = f.readlines()

    # - - - - - - - - - -  -
    survey_name_combined   = None
    survey_name_uncombined = None
    icol_band = -9
    # search contents to extract SURVEY and FILTERS and column with FLT/BAND
    for line in contents:
        words = line.split()
        if len(words) == 0 : continue 

        if words[0] == KEY_SURVEY :
            survey_name_combined = words[1]

        # check for FILTERS: [all] # uncombined -> combined
        if words[0] == KEY_FILTERS :
            filt_list_all            = list(words[1])
            filt_list_uncombined     = list(words[3])
            filt_list_combined       = list(words[5])
            
        # find column with FLT or BAND
        if words[0] == KEY_VARLIST :
            for colname in COLNAME_BAND_LIST :
                if colname in words :
                    icol_band = words.index(colname)
                    
        if words[0] == KEY_OBS :
            break

    # - - - - - - - - - - 
    if icol_band < 0 :
        assert False, f"ERROR: Could not find {COLNAME_BAND_LIST}"

    # create dictionary to map filter char
    filter_map = {}
    for f_comb, f_uncomb in zip(filt_list_combined, filt_list_uncombined):
        filter_map[f_comb] = f_uncomb

    # strip uncombined survey name out of SURVEY_COMBINED(SURVEY_UNCOMBINED)
    survey_name_uncombined = survey_name_combined.split('(')[1][:-1]

    # construct name of output uncombined data version
    version_uncombined = f"{PREFIX_VERSION_UNCOMBINED}_{survey_name_uncombined}"
    
    # if there is no directory, create it
    FIRST  = version_uncombined not in uncombined_dict['version_list']
    EXISTS = os.path.exists(version_uncombined)

    # create uncombined folder on first data file
    if FIRST :
        if EXISTS:
            cmd_rm = f"rm -r {version_uncombined}"
            os.system(cmd_rm)
        os.mkdir(version_uncombined)
        uncombined_dict['version_list'].append(version_uncombined)
        uncombined_dict['survey_list'].append(survey_name_uncombined)
        uncombined_dict['filter_list'].append("".join(filt_list_uncombined))

    # loop again over contents and make modifications for uncombined version
    contents_uncombined = []
    for line in contents:
        line_uncombined = line
        words = line.split()
        if len(words) == 0 : 
            contents_uncombined.append(line_uncombined)
            continue         
        if words[0] == KEY_SURVEY :
            line_uncombined = line.replace(survey_name_combined,survey_name_uncombined)
        if words[0] == KEY_FILTERS :
            filt_str        = "".join(filt_list_uncombined)
            line_uncombined = f"{KEY_FILTERS}  {filt_str}  \n"

        if words[0] == KEY_OBS :
            filt_str_comb   = words[icol_band]
            f_comb          = filt_str_comb[-1]      # strip last char
            f_uncomb        = filter_map[f_comb]
            filt_str_uncomb = filt_str_comb[0:-1] + f_uncomb

            # use pad spacing for filter-string replace to avoid modifying FIELD
            line_uncombined = line.replace(f" {filt_str_comb} ",
                                           f" {filt_str_uncomb} " )

        contents_uncombined.append(line_uncombined)

    # - - - - - - - 
    dump = False
    if dump :
        print(f" xxx survey_combined      = {survey_name_combined}")
        print(f" xxx survey_uncombined    = {survey_name_uncombined}")
        print(f" xxx filt_list_uncombined = {filt_list_uncombined}")
        print(f" xxx filt_list_combined   = {filt_list_combined}")
        print(f" xxx filter_map  = {filter_map}")
        print(f" xxx icol_band   = {icol_band}")
        sys.exit(" xxx BYE BYE")

    # write output contents to output data file
    data_file_uncombined = f"{version_uncombined}/{data_file}"
    with open(data_file_uncombined,"wt") as u:
        print(f"\t Create {data_file_uncombined}")
        for line in contents_uncombined:
            u.write(f"{line}")

    #sys.exit(" xxx BYE BYE")
    # end extract_data_file
    
def create_aux_files(version, survey, filters, args):
    
    # write LIST and README file for final uncombined data version

    logging.info(f" Create AUXILARY files for {version}")

    # start with LIST file
    list_file      = f"{version}/{version}.LIST"
    data_file_list = glob.glob1(version, "*.*")
    n_file = len(data_file_list)

    with open(list_file,"wt") as l:
        for data_file in data_file_list:
            l.write(f"{data_file}\n")

    # create README with bare-bones and request to add more info
    readme_file = f"{version}/{version}.README"
    cwd = os.getcwd() 

    with open(readme_file,"wt") as r:
        r.write(f"DOCUMENTATION:\n")
        r.write(f"  PURPOSE: uncombined lightcurve data files for {survey} \n")
        r.write(f"  FILTERS: {filters}\n")
        r.write(f"  CWD:     {cwd} \n")
        r.write(f"  NOTES: \n")
        r.write(f"  - {n_file} data files \n")
        r.write(f"  - extracted from {args.folder_version}\n")
        r.write(f"DOCUMENTATION_END:\n\n")
    # end create_aux_files

# ==============================
if __name__ == "__main__":

    #print(f"HEADER_OVERRIDE_FILE = \n{HEADER_OVERRIDE_FILE}\n")

    setup_logging()
    logging.info("# ========== BEGIN uncombine_dataVersions  ===============")
    args   = get_args() 
 
    data_file_list = get_file_list(args)
    
    uncombined_dict = { 
        'version_list' : [],  
        'survey_list'  : [],
        'filter_list'  : [] 
    }

    for data_file in data_file_list:
        extract_data_file(data_file, args, uncombined_dict)

    print('\n')

    # create aux files for each new version
    version_uncombined_list = uncombined_dict['version_list']
    survey_uncombined_list  = uncombined_dict['survey_list']
    filter_uncombined_list  = uncombined_dict['filter_list']
    for version, survey, filters in \
        zip(version_uncombined_list, survey_uncombined_list, filter_uncombined_list) :
        create_aux_files(version, survey, filters, args)

    logging.info(' Done')
    logging.info(' Remember to update README files')

# === end ===
       
