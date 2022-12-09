#!/usr/bin/env python
#
# April 2019, R. Kessler
# Called by makeFile to set C preprocessor flags based
# on the existance or absence of environment variables.
# This script replaces the older defineOutputFlags.pl
#
# Note that all input names are hard-coded here as globals.
#
# Apr 29 2019: replace USE_BYOSED ENV with SNANA_PYTHON_DIR
# May 16 2019: use relative path instead of $SNANA_DIR to allow
#              private/debug compile without overwriting
#              official SNANA install.
#
# Sep 10 2020: minor refactor with BYOSED -> PySEDMODEL
# Jun 10 2022: check VERSION_LIBPYTHON instead of SNANA_PYTHON_DIR
#
import os
import sys

# GLOBALS
SNANA_DIR          = os.environ['SNANA_DIR']
INCFILE_OUTPUT     = '../src/sntools_output.h'
INCFILE_PySEDMODEL = '../src/genmag_PySEDMODEL.h'
INCFILE_BAYESN     = '../src/genmag_BAYESN.h'

LIST_CFLAG     = [ 'USE_HBOOK' , 'USE_ROOT' , 'USE_PYTHON', 'USE_BAYESN' ]
# xxx mark LIST_ENV       = [ 'CERN_DIR'  , 'ROOT_DIR' , 'SNANA_PYTHON_DIR' ]
LIST_ENV       = [ 'CERN_DIR'  , 'ROOT_DIR' , 'VERSION_LIBPYTHON' , 'YAML_DIR']
LIST_INCFILE   = [ INCFILE_OUTPUT, INCFILE_OUTPUT, INCFILE_PySEDMODEL, INCFILE_BAYESN ]
NCFLAG         = len(LIST_CFLAG)


def check_Cflag(iflag):

    import subprocess
    import shutil
    ENV     =   LIST_ENV[iflag]
    ENV_val =   os.environ.get(ENV)
    CFLAG   =   LIST_CFLAG[iflag]
    INCFILE =   LIST_INCFILE[iflag]

    if ( ENV_val == None ):
        SET_CFLAG = False
    else:
        SET_CFLAG = True

    print('# --------------------------------------------- ')
    print(' Check C preproc flag: %s'%CFLAG)
    print(' File: %s'%INCFILE)
    print(' ENV(%s) = %s ' % (ENV, ENV_val))
    print(' SET_CFLAG = %s'%SET_CFLAG)

    # grep out CFLAG
    cmd     = 'grep define %s | grep %s ' % (INCFILE,CFLAG)
    if sys.version_info > (3,0):
        oldLine = subprocess.check_output(cmd, shell=True).rstrip().decode('utf-8')
    else:
        oldLine = subprocess.check_output(cmd, shell=True).rstrip()
    words   = oldLine.split()

    if ( CFLAG == words[1] ) :
        MATCH = True
    else:
        MATCH = False

    # return if flag is already correcly set, or correcly not set
    if ( SET_CFLAG and MATCH ):
        print('\t %s already set -> return'%CFLAG)
        return

    if ( not SET_CFLAG and not MATCH):
        print('\t %s not set -> return'%CFLAG)
        return

    # if we get here, change the flag in the file
    if ( SET_CFLAG and not MATCH ) :
        print(' Set %s pre-processor flag'%CFLAG)
        newLine = '#define %s ' % CFLAG
    else:
        print(' UNset %s pre-processor flag'%CFLAG)
        newLine = '#define %sxxx' % CFLAG


    # use sed and system call to make the substution
    sedCmd   = "sed -e 's/%s/%s/1'"%(oldLine,newLine)

    TEMPFILE = 'TEMP_INCLUDE_%s.h' % CFLAG
    os_sed   = '%s %s > %s ' % (sedCmd,INCFILE,TEMPFILE)
    os_mv    = 'mv %s %s' % (TEMPFILE,INCFILE)
    os_cmd   = os_sed + ';' + os_mv


#    print('\n xxx os_cmd = %s'%os_cmd )
    os.system(os_cmd)


# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":

    print('LIST_CFLAG = %s'%LIST_CFLAG )
    print('Number of C preproc flags: %s'%NCFLAG)

    for i in range(NCFLAG):
        check_Cflag(i)

    print('\n Done.')

# ========= END MAIN ================

