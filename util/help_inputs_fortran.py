#!/usr/bin/env python
#
# Created Auguest 2022 by R.Kessler
#
# ========================

import os, sys, argparse, glob, yaml
from   argparse import Namespace

LANG_FORTRAN = "fortran"
LANG_C       = "C"
LANG_PYTHON  = "python"

COMMENT_CHAR_C       = '//'
COMMENT_CHAR_FORTRAN = '!'
COMMENT_CHAR_PYTHON  = '#'

KEY_INPUT = "I:"

SNANA_DIR = os.environ['SNANA_DIR']

# =========================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "name of source code (snana.car snlc_fit.car  psnid.car snlc_sim.c SALT2mu.c"
    parser.add_argument("code_file", help=msg, nargs="?", default=None)

    #msg = "clobber everything and start over"
    #parser.add_argument("--clobber", help=msg, action="store_true")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    # if code_file does not include path, automatically append $SNANA_DIR/src
    if '/' not in args.code_file:
        code_basename = args.code_file
        args.code_file = f"{SNANA_DIR}/src/{code_basename}"
        pass

    return args

    # end get_args

def read_code_contents(code_file):

    if '.car' in code_file or '.f' in code_file:
        lang = LANG_FORTRAN
        comment_char = COMMENT_CHAR_FORTRAN
    elif '.c' in code_file or '.h' in code_file:
        lang = LANG_C
        comment_char = COMMENT_CHAR_C
    elif '.py' in code_file:
        lang = LANG_PYTHON
        comment_char = COMMENT_CHAR_PYTHON
    else:
        msgerr = "Unknown code language for {code_file}"
        sys.exit(msgerr)

    line_list = []
    nkeep = 0
    with  open(code_file,"rt") as f:
        for line in f:
            line  = line.rstrip() 
            if len(line) == 0 : continue
            if 'SUBROUTINE' in line: break
            if '+KEEP' in line: nkeep += 1
            line_list.append(line)
    
    nline = len(line_list)
    print(f" Read {nline} code lines from {code_file}")
    print(f" Code language is {lang}")
    print(f" Comment char is {comment_char}")
    if nkeep > 0 :
        print(f" Found {nkeep} common block sections")

    return line_list, lang, comment_char
    # end read_code_contents

def print_inputs(args,config):

    # after each 
    # +KEEP,XYZ
    # print every variable that has comment beginning with I:
    # Beware that many +KEEP blocks have no user input, so print
    # info only if some inputs are found.

    lang         = config['lang'] 
    contents     = config['code_contents']
    comment_char = config['comment_char']
    n_line = 0 
    n_input = 0
    for line in contents:
        n_line += 1
        if '+DECK' in line:
            break
        if '+KEEP' in line:
            block = line.split(',')[1].replace('.','')
            #print(f" xxx found KEEP {block}")
        if KEY_INPUT in line:
            tmp = line.split(KEY_INPUT)
            varname = tmp[0]
            varname = varname.replace('&','')
            varname = varname.replace(',','')
            comment = tmp[1]
            print(f" {varname} {comment}")
            n_input += 1

    print(f"\n Found {n_input} input variables from {n_line} code lines")
    return
    # end print_inputs

# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    args    = get_args()
    config      = {}

    code_contents, lang, comment_char = \
                read_code_contents(args.code_file)
    
    config['code_contents'] = code_contents
    config['lang']          = lang
    config['comment_char']  = comment_char

    print_inputs(args,config)

    # end main
