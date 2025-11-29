#!/usr/bin/env python   
#
# Custom script to convert
#   snana.car    -> snana.F90    [no longer need snana.cra]
#   snlc_fit.car -> snlc_fit.F90 [no longer need snlc_fit.cra]
#   psnid.car    -> psnid.F90 [no longer need psnid.cra]
#
# with usage
#   convert_f77_to_f90.py $SNANA_DIR/src/snana.car
#   convert_f77_to_f90.py $SNANA_DIR/src/snlc_fit.car
#   convert_f77_to_f90.py $SNANA_DIR/src/psnid.car
#
# This script has lots of custom accomodations, so it is not
# a general script to convert f77 code to f90.
#
# After running this script, the same compilation and link commands
# are used, but there is no need to run ypatchy.pl.
#
import os, sys, argparse, datetime
import subprocess

CAST_LIST = [ 'INTEGER', 'REAL', 'DOUBLE', 'CHARACTER', 'LOGICAL' ]
INDENT     = '  '
INDENT2    = '    '

PREFIX_MODULE   = ''  
MODULE_SNPAR    =  PREFIX_MODULE + 'SNPAR'
MODULE_SNCUTS   =  PREFIX_MODULE + 'SNCUTS'
MODULE_SNLCINP  =  PREFIX_MODULE + 'SNLCINP_NML'
MODULE_SNDATCOM =  PREFIX_MODULE + 'SNDATCOM'
MODULE_SNFITCOM =  PREFIX_MODULE + 'SNFITCOM'
MODULE_PSNIDPAR =  PREFIX_MODULE + 'PSNIDPAR'
MODULE_PSNIDCOM =  PREFIX_MODULE + 'PSNIDCOM'

# define dictionary of modules to rename to avoid conflicts
RENAME_MODULE_DICT = {
    'SNLCINP'   :  'SNLCINP_NML' ,
    'SNFITINP'  :  'FITINP_NML' ,
    'PSNIDINP'  :  'PSNIDINP_NML'
}

SNANA_DIR = os.getenv('SNANA_DIR')

BASE_F90_CODE  = "snana.F90"
BASE_CAR_CODE  = f"{SNANA_DIR}/src/snana.car"

tnow       = datetime.datetime.now()
DATE_STAMP = ('%4.4d-%2.2d-%2.2d' % (tnow.year,tnow.month,tnow.day) )


# ========================================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "name of f77 .car file"
    parser.add_argument("input_code_file", help=msg,
                        nargs="?", default=None)

    args = parser.parse_args()

    args.input_code_base = os.path.basename(args.input_code_file)

    args.is_snana =  'snana'    in args.input_code_file
    args.is_snfit =  'snlc_fit' in args.input_code_file
    args.is_psnid =  'psnid'    in args.input_code_file

    return args

def read_orignal_code(args):
    with open(args.input_code_file,"r") as c:
        code_lines = c.readlines()
        
    nline = len(code_lines)
    print(f"Read {nline} lines of original code from {args.input_code_file}")
    return code_lines

def open_output_code_file(args):
    base = os.path.basename(args.input_code_file)
    output_code_file = base.split('.')[0] + '.F90'

    print(f"Open output code file: {output_code_file}")
    o  = open(output_code_file,"wt")

    # write header comments at top
    full_command = ' '.join(sys.argv)
    o.write(f"! F77 -> F90 translation created {DATE_STAMP} with command:\n")
    o.write(f"!   {full_command}\n")
    o.write(f"\n")

    if not args.is_snana:
        o.write(f'! Include base snana code to read/write data and apply cuts\n')
        o.write(f'#include "{BASE_F90_CODE}" \n')
    o.write(f"\n")

    return o


def get_patchy_name(line):
    # if line = '+CDE,SNPAR' return SNPAR
    # if line = '+KEEP,SNPAR' return SNPAR
    tmp =  line.split('.')[0]
    name = tmp.split(',')[1]
    return name

def fix_line_CDE(line, keep_name_list):

    # if line = "+CDE,SNPAR" then return   "    USE SNPAR"

    name = get_patchy_name(line)

    if name not in keep_name_list:
        sys.exit(f"\n ERROR: found name = {name} from \n\t line = {line}\n\t but expected name from {keep_name_list}")

    if name in RENAME_MODULE_DICT:
        name = RENAME_MODULE_DICT[name]

    new_line = f"{INDENT2}USE {PREFIX_MODULE}{name}"
    return new_line

def fix_line_amper(line):
    # add & at end of line, but before ! if there is a comment

    if len(line) == 0 : return line
    if line[0] == '!' : return line

    if '!' in line:
        j    = line.index('!')
        new_line = line[0:j-1] + '  &  ' + line[j:]
    else:
        new_line = line + '  & '

    return new_line

def fix_line_comment(line):

    nch = len(line)
    new_line = line
    if nch > 0:
        if line[0] == 'c' or line[0] == 'C' : 
            if nch == 1 :
                new_line = '! '  
            else : 
                new_line='! ' + line[2:]

    return new_line

def fix_line_SELF(line, last_line, preproc_keylist ):

    #print(f" xxx process SELF for line = {line}")
    line_nodot    = line.split('.')[0]  # remove dot and anything after

    if 'IF=' in line:
        var_string    = line_nodot.split('IF=')[1]
        var_list      = var_string.split(',')

        macro_list = []
        for var in var_list:   
            macro_list.append(f"defined({var})")
            if var not in preproc_keylist: preproc_keylist.append(var)

        if 'endif' in last_line:
            new_line = '#if ' + ' || '.join(macro_list)
        else:
            new_line = '#elif ' + ' || '.join(macro_list)

        #print(f"xxx \t {line}   -->  {new_line}")
    else:
        new_line = '#endif'

    return new_line

def fix_continuations(line_list):
    # fix continuations by moving & to previous line; e.g.
    #   INTEGER
    #  &  J   ! comment about J
    #  &  K   ! comment about K
    #   is replaced with
    #   INTEGER &
    #    J      &  ! comment about J
    #    K         ! comment about K

    #new_line_list = []
    n_line = len(line_list)

    for i in range(0,n_line-1):
        line      = line_list[i]

        if len(line) == 0   : continue
        if line[0] == '!'   : continue

        i2 = i
        while ( i2 < n_line-1 ) :
            i2 += 1
            line_tmp = line_list[i2] + ' '
            is_comment = line_tmp[0]  == '!'  or (line_tmp.split('!')[0]).isspace() 
            if not is_comment: break

        line_next = line_list[i2]
        if len(line_next) == 0 : continue

        #print(f" xxx check line = {line} ")
        #sys.stdout.flush()

        amper      = line.lstrip()[0]      == '&'
        amper_next = line_next.lstrip()[0] == '&'
        if amper_next:
            line      = fix_line_amper(line)
            line_next = line_next.replace('&',' ')
            line_list[i]   = line
            line_list[i2]  = line_next

    return line_list

def get_name_subroutine(line):
    tmp = line.split('SUBROUTINE ')[1]
    if '(' in tmp:
        name = tmp.split('(')[0]
    else:
        name = tmp
    return name

def get_name_function(line):
    tmp = line.split('FUNCTION ')[1]
    if '(' in tmp:
        name = tmp.split('(')[0]
    else:
        name = tmp
    return name


def get_lines_append_module(args, name_module):

    # for newly determined module, determine line to append using a previous module.
    # This is all hard-wired.

    lines = []
    
    if name_module != MODULE_SNPAR and name_module != MODULE_SNDATCOM and name_module != MODULE_SNFITCOM :  
        lines.append(f"{INDENT2}USE {MODULE_SNPAR}")
                
    if name_module == MODULE_SNLCINP:
        lines.append(f"{INDENT2}USE {MODULE_SNCUTS}")

    if name_module == f"{PREFIX_MODULE}FILTUPDCM" :
        lines.append(f"{INDENT2}USE {PREFIX_MODULE}FILTCOM")

    # for SNFIT
    if name_module == f"{PREFIX_MODULE}SNFITVAR" :
        lines.append(f"{INDENT2}USE {PREFIX_MODULE}SNFITPAR")
             
    if name_module == f'{PREFIX_MODULE}MAGDIFCOM'  or name_module == f'{PREFIX_MODULE}FITRESTCOM' :
        lines.append(f"{INDENT2}USE {PREFIX_MODULE}ALLFILTCOM")

    # for PSNID
    if args.is_psnid :
        if name_module != MODULE_PSNIDPAR and name_module != MODULE_PSNIDCOM :
            lines.append(f"{INDENT2}USE {MODULE_PSNIDPAR}")
        
    # - - - 
    is_group = (name_module == MODULE_SNDATCOM or name_module == MODULE_SNFITCOM or \
                name_module == MODULE_PSNIDCOM)
    if not is_group:
        lines.append(f"{INDENT2}IMPLICIT NONE")

    debug_print = False
    if debug_print:
        print(f" module {name_module} --> append {lines}")

    return lines
    
def translate_to_F90(KEEP_LIST, code_lines,o):

    # convert all common blocks under KEEP to MODULE
    
    STR_KEEP = 'KEEP'
    STR_DECK = 'DECK'

    n_line_tot = {  STR_KEEP : 0,        STR_DECK  : 0    }
    n_line_wr  = {  STR_KEEP : 0,        STR_DECK  : 0    }

    n_KEEP = 0
    n_DECK = 0
    n_SUB  = 0
    n_FUN  = 0
    n_MAIN = 0

    F90_line_list  = []  # store all new lines in a list before writing them out to F90 file
    skip_COMMON = False  # reset skip logical
    sub_name     = None
    fun_name     = None
    found_CDE_last  = False
    PREPROC_KEYLIST = []
    CODE_SECTION = STR_KEEP

    debug_flag = False

    line_equal = '! ====================================================================='
    last_line_self  = 'endif'

    for line in code_lines:

        n_line_tot[CODE_SECTION] += 1

        line = line.rstrip()  # remove trailing space and linefeed
        line = fix_line_comment(line)

        keep_line  = True
        blank_line = len(line)==0 or line.isspace()

        line_nocomment = line.split('!')[0]
        found_KEEP     = line[0:5]  == "+KEEP"
        found_END_KEEP = 'END_KEEP' in line
        found_PATCH   = line[0:6]  == "+PATCH" 
        found_CDE     = line[0:4]  == "+CDE"
        found_COMMON  = ('COMMON' in line or 'common' in line) and '/' in line
        found_DECK    = line[0:5]  == "+DECK"
        found_USE     = line[0:4]  == "+USE"
        found_SUB     = line[6:16] == 'SUBROUTINE'
        found_MAIN    = 'MAIN' in line_nocomment and 'PROGRAM' in line_nocomment
        found_CAST    = any(item in line for item in CAST_LIST)
        found_FUN     = 'FUNCTION' in line and found_CAST
        found_END     = line[6:9]  == 'END'   \
                        and 'DO' not in line_nocomment   \
                        and 'IF' not in line_nocomment
        found_USR     = 'USRINI' in line or 'USRANA' in line or 'USREND' in line
        
        found_SELF    = line[0:5]  == '+SELF'
        found_IMPNONE = 'IMPLICIT' in line and 'NONE' in line
        found_DELETE  = 'DELETE_F90' in line

        write_end_module = n_KEEP > 0 and (found_KEEP or found_END_KEEP)

        if write_end_module: 
            F90_line_list.append(f"{INDENT}END MODULE {name_F90}")
            
        if found_KEEP:
            n_KEEP += 1
            name_orig  = get_patchy_name(line)

            if name_orig in RENAME_MODULE_DICT :
                tmp = name_orig
                name_orig = RENAME_MODULE_DICT[name_orig]
                print(f" rename {tmp} --> {name_orig}")

            name_F90   = PREFIX_MODULE + name_orig
            F90_line_list.append(f"")
            F90_line_list.append(f"{line_equal}")
            F90_line_list.append(f"{INDENT}MODULE {name_F90}")

            lines_append = get_lines_append_module(args, name_F90)
            if len(lines_append) > 0 :
                F90_line_list += lines_append

            keep_line   = False  # do not write this line to F90 file

        if found_COMMON:
            skip_COMMON = True
        if skip_COMMON and blank_line:
            skip_COMMON = False            


        # remove  stuff that is not needed in F90 file
        if found_IMPNONE or found_PATCH or found_END_KEEP or skip_COMMON or \
           found_USE or found_USR or found_DELETE :
            keep_line = False

        if found_CDE: 
            line = fix_line_CDE(line, KEEP_LIST) # replace +CDE,XXX with USE XXX

        # IMPLICIT NONE statements in input code are before +CDE (USE) statements,
        # but they need to come after for F90
        if blank_line and found_CDE_last and n_DECK > 0:
            F90_line_list.append(f"")
            F90_line_list.append(f"{INDENT2}IMPLICIT NONE")
        found_CDE_last = found_CDE


        if found_DECK:
            CODE_SECTION = STR_DECK
            keep_line = False
            n_DECK += 1
            line_DECK = line
            if n_DECK == 1:
                F90_line_list.append(f"{line_equal}")
                F90_line_list.append(f"{line_equal}")
                F90_line_list.append(f"{line_equal}")            
            if debug_flag: print(f"  {line_DECK} -> ")

        if found_MAIN: 
            n_MAIN += 1
            if debug_flag : print(f"\t MAIN")

        if found_SUB:
            sub_name = get_name_subroutine(line)
            n_SUB += 1
            if debug_flag : print(f"\t SUBROUTINE {sub_name}")

        if found_FUN:
            fun_name = get_name_function(line)
            n_FUN += 1
            if debug_flag : print(f"\t FUNCTION {fun_name}")

        if found_SELF:
            line = fix_line_SELF(line, last_line_self, PREPROC_KEYLIST)
            last_line_self = line

        if found_END :
            if sub_name is not None:
                line = f"{INDENT}END SUBROUTINE {sub_name}"
                sub_name = None
            if fun_name is not None:
                line = f"{INDENT}END FUNCTION {fun_name}"
                fun_name = None

        # - - - - - - -  -
        if keep_line:
            n_line_wr[CODE_SECTION] += 1
            if line[0:6] == '      ' :   line = line[2:]
            F90_line_list.append(line)
            

    # -------------------------------------------------------------
    F90_line_list = fix_continuations(F90_line_list)

    for line in F90_line_list:
        o.write(f"{line}\n")

    n_TOT = n_MAIN + n_SUB + n_FUN
    print(f"Done converting ")
    print(f"\t {n_KEEP:3d} KEEPs to MODULEs  (wrote {n_line_wr[STR_KEEP]} of {n_line_tot[STR_KEEP]} lines)")
    print(f"\t {n_DECK:3d} DECKs to MODULEs  (wrote {n_line_wr[STR_DECK]} of {n_line_tot[STR_DECK]} lines)")
    print(f"  PREPROC directives: {PREPROC_KEYLIST}")
    print(f"  Number of MAIN + SUBROUTINE + FUNCTION : {n_MAIN} + {n_SUB} + {n_FUN} = {n_TOT}")
    if n_TOT != n_DECK:
        print(f"\t WARNING: n_DECK != Sum of MAIN+SUBR+FUNC")
    return



def grep_KEEP_list(code_file):
    
    # grep out +KEEP and return list of KEEP names
    KEEP_LIST = []

    command = ["grep", "+KEEP", code_file ]
    
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    temp_list = (result.stdout).split()

    for temp in temp_list:
        name = get_patchy_name(temp)
        KEEP_LIST.append(name)


    return KEEP_LIST

# ===================================================
if __name__ == "__main__":

    args = get_args()

    # use grep to fish out expected KEEP names
    KEEP_LIST = grep_KEEP_list(BASE_CAR_CODE)
    if args.input_code_base not in BASE_CAR_CODE:
        KEEP_LIST += grep_KEEP_list(args.input_code_file)

    #sys.exit(f"\n xxx KEEP_LIST = \n{KEEP_LIST}")

    # - - - - -
    code_lines = read_orignal_code(args)
    o          = open_output_code_file(args)
    
    translate_to_F90(KEEP_LIST, code_lines, o)

    # === END:


    
