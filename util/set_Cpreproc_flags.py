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
# Jun 20 2023: remove USE_BAYESN/YAML_DIR/INCFILE_BAYESN from lists
# May 15 2025: remove obsolete hbook dependence (preparing for automake)
#
import os
import sys
import subprocess
import shutil

# GLOBALS
SNANA_DIR          = os.environ['SNANA_DIR']
INCFILE_OUTPUT     = '../src/sntools_output.h'
INCFILE_PySEDMODEL = '../src/genmag_PySEDMODEL.h'
INCFILE_BAYESN     = '../src/genmag_BAYESN.h'
INCFILE_ZPDF_SPLINE= '../src/sntools_zPDF_spline.h'


# xxxxxxxxxxxxxx mark delete May 15 2025 (get rid of hbook) 
#xxxLIST_CFLAG     = [ 'USE_HBOOK' , 'USE_ROOT' , 'USE_PYTHON', 'GSL_INTERP_STEFFEN' ]
#xxxLIST_ENV       = [ 'CERN_DIR'  , 'ROOT_DIR' , 'VERSION_LIBPYTHON', 'GSL_DIR 2.0' ]
#xxxLIST_INCFILE   = [ INCFILE_OUTPUT, INCFILE_OUTPUT, INCFILE_PySEDMODEL, INCFILE_ZPDF_SPLINE ]
#xxxxxxxxxxx end mark 

LIST_CFLAG     = [ 'USE_ROOT' , 'USE_PYTHON', 'GSL_INTERP_STEFFEN' ]
LIST_ENV       = [ 'ROOT_DIR' , 'VERSION_LIBPYTHON', 'GSL_DIR 2.0' ]
LIST_INCFILE   = [ INCFILE_OUTPUT, INCFILE_PySEDMODEL, INCFILE_ZPDF_SPLINE ]
NCFLAG         = len(LIST_CFLAG)



def fetch_version(ENV):
    '''
    May 28/2024
    '''
    if 'GSL' in ENV :
        command = 'gsl-config --version'
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output  = process.communicate()[0]
        version =  output.decode()
    elif 'ROOT' in ENV :
        version = None  # Maybe Later
    else :
        version = None

    return version    
        

def check_Cflag(iflag):

    ENV_ITEMS      =   LIST_ENV[iflag].split()
    NITEMS         =   len(ENV_ITEMS)
    ENV            =   ENV_ITEMS[0]
    ENV_val        =   os.environ.get(ENV)
    CFLAG          =   LIST_CFLAG[iflag]
    INCFILE        =   LIST_INCFILE[iflag]
    SET_CFLAG      = False
    
    if (NITEMS == 2 ):
        VERSION_WANT = ENV_ITEMS[1]
        VERSION_ACTUAL = fetch_version(ENV)
    else :
        VERSION_ACTUAL = 'unknown'

    if ( ENV_val is not  None ):
        #SET_CFLAG = True
        if (NITEMS == 2 ):
            # Require version to SET_CFLAG
            VERSION_WANT   = ENV_ITEMS[1]
            VERSION_ACTUAL = fetch_version(ENV)
            VERSION_LIST   = [VERSION_WANT, VERSION_ACTUAL]
            VERSION_SORT   = sorted(VERSION_LIST)
            if (VERSION_SORT[1] == VERSION_ACTUAL ):
                SET_CFLAG = True
        else :
            SET_CFLAG = True # ENV exists so SET_CFLAG 

    print('# --------------------------------------------- ')
    print(' Check C preproc flag: %s'%CFLAG)
    print(' File: %s'%INCFILE)
    print(' ENV(%s) = %s ' % (ENV, ENV_val))
    if (NITEMS == 2):
        print('%s Version = %s'%(ENV, VERSION_ACTUAL))
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

