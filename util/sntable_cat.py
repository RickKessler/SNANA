#!/usr/bin/env python
#
# Created May 2020 by R.Kessler
#
# Read input list of text-formatted SN table files (i.e., FITRES files)
# and catenate them after accounting for possible column differences
# among input files. See manual Sec 10.1 for explanations.
#
# This python script is just a wrapper that calls SALT2mu.exe
# with 'cat_only' argument that turns SALT2mu.exe into a utility 
# instead of a fitting program.
#
# See input arguments with
#    sntable_cat.py -h
# ------------------------------------------

import os,sys
from optparse import OptionParser

# =============================
def parse_args():
    parser = OptionParser()

    msg = 'comma-sep list of text-formatted SNANA tables'
    parser.add_option('-i',help=msg,action='store',type='string',
                      dest='inpfile_list',default=None)

    msg = 'name of output/catenated sntable (text format)'
    parser.add_option('-o',help=msg,action='store',type='string',
                      dest='outfile_cat',default=None)

    msg = 'list of column names to append if missing'
    parser.add_option('-a',help=msg,action='store',type='string',
                      dest='append_varname_missing',default='PROB*')

    (INPUTS,ARGS)=parser.parse_args()

# abort on missing arguments
    if (not INPUTS.inpfile_list ):
        parser.error('Missing required arg: -i <inpfile_list>')

    if (not INPUTS.outfile_cat ):  
        parser.error('Missing required arg: -o <outfile_cat>')

    return INPUTS


# ===================================
def exec_salt2mu(INPUTS):

    inpfile_list   = INPUTS.inpfile_list
    outfile_cat    = INPUTS.outfile_cat
    append_varname = INPUTS.append_varname_missing
 
    arg_list = ("datafile=%s catfile_out=%s append_varname_missing=%s"
                % (inpfile_list,outfile_cat,append_varname) )

    command = ("SALT2mu.exe cat_only %s" % arg_list )
    print(" Execute command = '%s' " % command)
    os.system(command)

# ===================================
# ============== MAIN ===============
# ===================================

if __name__ == "__main__":

    INPUTS = parse_args()

# execute SALT2mu.exe program to do the catenate
    exec_salt2mu(INPUTS)



