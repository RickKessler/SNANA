#!/usr/bin/env python
#
# Created July 2021 by R.Kessler
# Program to
#   + read data from user-specified source
#       (e.g., LSST-AP, LSST-DRP, SIRAH, ZTF, or test with data folder)
#   + write data files in human-readable TEXT format (1 file per event)
#   + convert TEXT -> FITS format
#   + tar up TEXT files
#
# For batch distribution, a separate submit-script should be able to
# distribute jobs that are split by
#   + field (e.g., DDF or WFD)
#   + season (1-10)
#   + random sub-sample (see --nsplit and --isplit  inputs)
#
# A "data_unit" corresponds to a single FIELD-SEASON-SUBSAMPLE.
# A merge process combines data_units into more useful samples
# such as entire season and ALL data.
#
# To-DO list:                    
#    - write PHOTFLAG_REJECT_BITS to IGNORE files      
#    - add makeDataFile class to submit_batch framework 
#    - update snlc_sim to create season map from SIMLIB;  
#         read season map here to determin integer YY index   
#    - modify framework to write output directly without  
#       intermediate text files.
#
# ============================================

import sys, yaml, argparse, subprocess, logging, glob
import makeDataFiles_util  as util

from   makeDataFiles_params    import *
from   read_data_lsst_ap       import data_lsst_ap
from   read_data_lsst_drp      import data_lsst_drp
from   read_data_sirah_folder  import data_sirah_folder
from   read_data_snana_folder  import data_snana_folder

# =====================================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "Data source: LSST AP"
    parser.add_argument("--lsst_ap", help=msg, action="store_true")

    msg = "Data source: LSST_DRP"
    parser.add_argument("--lsst_drp", help=msg, action="store_true")    

    msg = "Data source: SIRAH pkl folder"
    parser.add_argument("--sirah_folder", help=msg, type=str, default=None )

    msg = "Data source: SNANA sim-data folder (for testing)"
    parser.add_argument("--snana_folder", help=msg, type=str, default=None )

    msg = "field name (e.g., SHALLOW, DEEP, etc ..)"
    parser.add_argument("--field", help=msg, type=str, default=FIELD_VOID )

    msg = "SNANA format: output top-directory for data"
    parser.add_argument("--outdir_snana", help=msg, type=str, default=None )

    msg = "XYZ format: output top-directory for data"
    parser.add_argument("--outdir_XYZ", help=msg, type=str, default=None )    
    
    msg = "number of random sub-samples (default=1)"
    parser.add_argument("--nsplitran", help=msg, type=int, default=1 )

    msg = "select isplitran (1-nsplitran): default=-1 -> all"
    parser.add_argument("--isplitran", help=msg, type=int, default=-1 )

    msg = "select year index (1-Nyear); default = -1 -> all"
    parser.add_argument("-y", "--year", help=msg, type=int, default=-1 )

    msg = "number of events to process (default is all)"
    parser.add_argument("--nevt", help=msg, type=int, default=99999999 )    

    msg = "increase output verbosity (default=True)"
    parser.add_argument("-v", "--verbose", help=msg, action="store_true")

    msg = "leave text files uncompressed"
    parser.add_argument("--text", help=msg, action="store_true")

    msg = "merge SNANA files (after split jobs finish)"
    parser.add_argument("--merge_snana", help=msg, action="store_true")

    msg = "process fakes (default is real data)"
    parser.add_argument("--fake", help=msg, action="store_true")    

    # - - - -
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()

    # end get_args

def restore_args_from_readme(args,readme_yaml):

    # restore user args from readme_yaml that was read from README file.
    
    args.lsst_ap      = False
    args.lsst_drp     = False
    args.sirah_folder = None
    args.snana_folder = None
            
    key = 'SOURCE_LSST_AP'
    if key in readme_yaml:
        args.lsst_ap   = readme_yaml[key]

    key = 'SOURCE_LSST_DRP'
    if key in readme_yaml:
        args.lsst_drp  = readme_yaml[key]

    key = 'SOURCE_SIRAH_FOLDER'
    if key in readme_yaml:
        args.sirah_folder  = readme_yaml[key]

    key = 'SOURCE_SNANA_FOLDER'
    if key in readme_yaml:
        args.snana_folder  = readme_yaml[key]

    args.field        = readme_yaml['FIELD']

    
    # end restore_args_from_readme
    
def which_program_class(args):
    program_class = None 

    # if merge process, read any one of the README files to 
    # recover the user-input logicals
    if args.merge_snana:
        # restore args for merge process.
        outdir      = args.outdir_snana
        folder      = glob.glob1(outdir, f"[!_TEXT]*" )[0]        
        readme_file = f"{outdir}/{folder}/{folder}.README"
        readme_yaml = util.read_yaml(readme_file)
        restore_args_from_readme(args,readme_yaml[DOCANA_KEY])
        
    if args.lsst_ap:
        program_class = data_lsst_ap
        args.survey = "LSST"
    elif args.lsst_drp :
        program_class = data_lsst_drp
        args.survey = "LSST"
    elif args.sirah_folder is not None :        
        program_class = data_sirah_folder
        args.survey = "SIRAH"
    elif args.snana_folder is not None :
        program_class     = data_snana_folder
        snana_folder_base = os.path.basename(args.snana_folder)
        args.survey   = util.get_survey_snana(snana_folder_base)
    else:
        sys.exit("\nERROR: Could not determine program_class")

    
    return program_class

    # end which_program_class    

# =============================================
if __name__ == "__main__":

    args  = get_args()
    store = util.setup_logging(args)

    # determine which program class (AP, DRP, test data)
    program_class  = which_program_class(args)

    # store inputs; for now just include command-line args, but later
    # might include yaml contents from input file.
    config_inputs =  { 'args' : args }

    # init the data-source class
    program = program_class(config_inputs)  # calls __init only

    if args.merge_snana:
        program.merge_snana_driver()
        sys.exit(' Done with merge: exiting Main.')

    # read data and write each event to text-format data files;
    # these intermediate text files are useful for visual debugging,
    # and they will be translated to binary below.
    program.read_data_driver()

    # translate TEXT -> BINARY; allow multiple output formats
    if args.outdir_snana is not None:
        program.convert2fits_snana()

    #if args.outdir_XYZ is not None:
    #    program.convert2XYZ()        

    # final summary
    program.final_summary()
    
    # === END ===
