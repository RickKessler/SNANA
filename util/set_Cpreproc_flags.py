#!/usr/bin/env python
#
# April 2019, R. Kessler
# Called by makeFile to set C preprocessor flags based
# on the existance or absence of environment variables.
# This script replaces the older defineOutputFlags.pl
#
# Note that all input names are hard-coded here as globals.
#

import os
import sys 

# GLOBALS
SNANA_DIR      = os.environ['SNANA_DIR']
#SNANA_DIR      = '/home/s1/rkessler/snana_debug/snana'
INCFILE_OUTPUT = SNANA_DIR + '/src/sntools_output.h'
INCFILE_BYOSED = SNANA_DIR + '/src/genmag_BYOSED.h'

LIST_CFLAG     = [ 'USE_HBOOK' , 'USE_ROOT' , 'USE_PYTHON' ]
LIST_ENV       = [ 'CERN_DIR'  , 'ROOT_DIR' , 'USE_BYOSED' ] 
LIST_INCFILE   = [ INCFILE_OUTPUT, INCFILE_OUTPUT, INCFILE_BYOSED ]
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


    print('\n xxx os_cmd = %s'%os_cmd )
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

