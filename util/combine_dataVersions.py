#!/usr/bin/env python
#
# Created Feb 2017 by R.Kessler
#
# @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
#
#  WARNING: Bug in SUBSURVEY_LIST key in output SIMLIB  file 
#
# @!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
#
#
# Combine multiple data versions into a single data version.
# Also produced combined kcor file and simLib file, and convert
# text data into FITS format.
#
# Since the same filter char can appear in different data versions,
# replace filters with sequential 'abcdef ... ABCDEF ... 0123456789'
# Intended use is for low-z, but will work with any data samples
# in TEXT format.
#
# Usage:
#   combine_dataVersions.py <inFile>
#
# where <inFile> contains
#   PRIVATE_DATA_PATH:  <path>   # this is optional
#   VERSION:  <ver1>  <kcor_inFile1>
#   VERSION:  <ver2>  <kcor_inFile2>
#     etc ...
#   SURVEY_OUT:  <name of combined survey>
#   VPEC_FILE:   <name file file with VPEC & VPEC_ERR>
#
# Outputs:
#    <SURVEY_OUT>_TEXT/       ! combined data directory, TEXT format
#    <SURVEY_OUT>_FITS/       ! combined data directory, FITS format
#    kcor_<SURVEY_OUT>.fits   ! combined kcor file
#    <SURVEY_OUT>.SIMLIB      ! combined SIMLIB file
#
#      History
#  Apr 12 2017: D.Scolnic added key SIMLIB_ZPERR_LIST
#
#  Oct 25 2017: RK comment out all SIMLIB_ZPERR_LIST code since it
#               causes code to crash.
#
#  Dec 8 2017: S.Hinton - merge duplicates
#  Jan 8 2018: RK - add VPEC_FILE option
#
#  Nov 16 2019 RK - bug fix writing AB_SED and BD17_SED
#
# ====================================

import os
import sys 
import numpy as np
import getpass
import time
import shutil
import string
# globals

SNDATA_ROOT = os.environ['SNDATA_ROOT']
HOSTNAME    = os.environ['HOSTNAME']
NOW         = time.strftime("%c")
CWD         = os.getcwd() 
MXFILTERS   = 62
TOPDIR_DATA = SNDATA_ROOT + '/lcmerge'
TOPDIR_KCOR = SNDATA_ROOT + '/kcor'

FILTER_CHARLIST = 'abcdefghijklmnopqrstuvwxyz' + 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' + '0123456789'


KEYNAME_VERSION        = "VERSION:"
KEYNAME_SURVEY_OUT     = "SURVEY_OUT:"
KEYNAME_VPEC_FILE      = "VPEC_FILE:"
KEYNAME_SIMLIB_ZPERR_LIST = "SIMLIB_ZPERR_LIST:"

# ========== BEGIN ===========

def parseLines(Lines,key,narg,vbose):
    # Lines is input array of lines in file
    # key is key to search
    # narg is number of args to return after key

    arg     = []
    rowList = Lines[np.char.startswith(Lines,key)]
    nrow    = len(rowList)

    if ( nrow == 1 ):
        if ( narg==1 ):
            arg = rowList[0].split()[1]
        else:
            arg = rowList[0].split()[1:narg+1]
    elif ( nrow > 1 ):
        for row in rowList:
            arg.append(row.split()[1:narg+1])

    if ( vbose > 0 ):
        print '\t ', key, arg
        
    return(arg)


def change_filterChar(versionInfo,kcorInfo):
    nver=0    # version index = kcor index
    for kcor in kcorInfo :
        # check of this filter/system has changed since last one
        new = 1
        if ( nver > 0 ):
            new=0
            if kcor.MAGSYSTEM != kcor_last.MAGSYSTEM:
                new=1
            if kcor.FILTSYSTEM != kcor_last.FILTSYSTEM:
                new=1
            if kcor.FILTPATH != kcor_last.FILTPATH:
                new=1

        kcorInfo[nver].NEW = new   
        kcor_last = kcor
        
        for line in kcor.FILTER :
            filter_name  = line[0]  # full name of filter
            N = versionInfo.NFILTER_TOT
            filter_oldChar = filter_name[-1]     # current filter char
            filter_newChar = FILTER_CHARLIST[N]  # new new character for filter
            versionInfo.FILTER_CHARLIST_OLD += filter_oldChar
            versionInfo.FILTER_CHARLIST_NEW += filter_newChar
            versionInfo.NFILTER_TOT    += 1

            # append new char to end of filter name to ensure
            # unique filter-char for each band
            filter_name_new = filter_name + '/' + filter_newChar
            kcorInfo[nver].FILTER_NEWNAME.append(filter_name_new) # full name
            kcorInfo[nver].FILTER_CHARLIST_OLD += filter_oldChar # 1-char
            kcorInfo[nver].FILTER_CHARLIST_NEW += filter_newChar


        # print band replacement for each version
        OLD = kcorInfo[nver].FILTER_CHARLIST_OLD
        NEW = kcorInfo[nver].FILTER_CHARLIST_NEW
        tmp = OLD + ' --> ' + NEW
        print '\t', versionInfo.NAME[nver], ' bands: ', tmp
        
        nver += 1

        


def write_kcor_inputFile(versionInfo,kcorInfo):

    fname = versionInfo.kcor_inFile
    print '\n Create ', fname
    f = open(fname,"wt")

    f.write("# Combined kcor file create by\n#   %s\n" % sys.argv )
    f.write("# User = %s\n" % ( getpass.getuser())  )
    f.write("# Host = %s\n" % ( HOSTNAME) )
    f.write("# %s\n" % NOW )
    f.write("\n")
        
    f.write( "SN_SED:       %s\n" % ( kcorInfo[0].SN_SED ) )    

    # check all kcor files for primary(s)
    nkcor=0
    USE_BD17 = 0 ; USE_AB = 0
    for kcor in kcorInfo :
        if ( len(kcorInfo[nkcor].BD17_SED) > 0 and USE_BD17==0 ) :
            f.write("BD17_SED:     %s\n" % (kcorInfo[nkcor].BD17_SED ) )
            USE_BD17 = 1

        if ( len(kcorInfo[nkcor].AB_SED) > 0  and  USE_AB==0) :
            f.write("AB_SED:       %s\n" % (kcorInfo[nkcor].AB_SED ) )
            USE_AB = 1
        nkcor += 1


    f.write( "LAMBDA_RANGE: %s %s \n" %
             ( kcorInfo[0].LAMRANGE[0], kcorInfo[0].LAMRANGE[1] ) )

    f.write("OUTFILE:      %s\n" % ( versionInfo.kcor_outFile ) )
    
    # loop over filter sets
    nkcor=0
    for kcor in kcorInfo :

        V = versionInfo.NAME[nkcor]
        nkcor += 1
                
        f.write("\n" )
        f.write("# Start filters for VERSION = %s\n" % V)
        if kcor.NEW == 1:
            f.write("MAGSYSTEM:   %s\n" %  ( kcor.MAGSYSTEM  ) )
            f.write("FILTSYSTEM:  %s\n" %  ( kcor.FILTSYSTEM ) )
            f.write("FILTPATH:    %s\n" %  ( kcor.FILTPATH   ) )
        kcor_last = kcor

        nfilt=0
        for line in kcor.FILTER :
            filter_name  = line[0]
            filter_file  = line[1]
            filter_zpoff = line[2]
            filter_name_new = kcor.FILTER_NEWNAME[nfilt]
            f.write("FILTER:  %-20.20s %s  %s \n"
                    % (filter_name_new,filter_file,filter_zpoff) )
            nfilt += 1
    f.close

def run_kcor(versionInfo):
    inFile  = versionInfo.kcor_inFile
    logFile = versionInfo.kcor_logFile
    cmd = 'kcor.exe ' + inFile + ' > ' + logFile
    print ' Run kcor program ... '
    os.system(cmd)

    # check for fatal error in log file
    f      = open(logFile,"rt")
    Lines  = f.readlines()
    f.close
    if any("ABORT" in s for s in Lines):
        msg = "\nFATAL error running kcor program:\n   Check %s\n" % inFile
        sys.exit(msg)
    
class VERSION_INFO:
    def __init__(self,filename):
        f = open(filename,"rt")
        Lines = f.readlines()

        # read name output version
        reader = np.array([x.split() for x in Lines
                           if x.startswith(KEYNAME_SURVEY_OUT)])
                
        SOUT   = reader[0:,1][0]
        self.SURVEY_OUT         = SOUT
        self.VERSION_OUT_TEXT   = SOUT + '_TEXT'
        self.VERSION_OUT_FITS   = SOUT + '_FITS'
        self.AUXFILE_README     = SOUT + '_TEXT.README'
        self.AUXFILE_IGNORE     = SOUT + '_TEXT.IGNORE'
        self.AUXFILE_LIST       = SOUT + '_TEXT.LIST'

        
        self.kcor_inFile  = 'kcor_' + SOUT + '.input'
        self.kcor_outFile = 'kcor_' + SOUT + '.fits'
        self.kcor_logFile = 'kcor_' + SOUT + '.log'
        self.simlibFile   = SOUT + '.SIMLIB'
        
        reader = np.array([x.split() for x in Lines
                           if x.startswith('PRIVATE_DATA_PATH:')])
        if ( len(reader) > 0 ):
            tmpDir = reader[0:,1][0]
            self.PRIVATE_DATA_PATH = os.path.expandvars(tmpDir)
        else:
            self.PRIVATE_DATA_PATH = ""


        # read optional name VPEC_FILE
        reader = np.array([x.split() for x in Lines
                           if x.startswith(KEYNAME_VPEC_FILE)])
        if ( len(reader) > 0 ):
            tmpFile = reader[0:,1][0]
            self.VPEC_FILE = tmpFile
        else:
            self.VPEC_FILE = ""

        
        # define strings of old and new filter char;
        # to be filled later.
        self.FILTER_CHARLIST_OLD = ""
        self.FILTER_CHARLIST_NEW = ""
        self.NFILTER_TOT    = 0
        self.NFILE          = []   # to be filled later
        
        # read input versions
        reader = np.array([x.split() for x in Lines
                           if x.startswith(KEYNAME_VERSION)])
        self.NAME        = reader[:,1]
        self.INFILE_KCOR = reader[:,2]
        
        f.close()

        
class KCOR_INFO:
        def __init__(self,filename):
            print '   Parse kcor input file: ', filename

            filename_expandvars = os.path.expandvars(filename)
            # check local dir first; then check $SNDATA_ROOT/kcor
            if os.path.isfile(filename_expandvars):
                fname_local = filename_expandvars
            else:
                fname_local = TOPDIR_KCOR + '/' + filename_expandvars

            # open file file and read all lines into Lines
            f = open(fname_local,"rt")
            Lines  = np.array(f.readlines())

            # parse SED stuff
            self.SN_SED   = parseLines(Lines,'SN_SED:',  1, 1)  
            self.BD17_SED = parseLines(Lines,'BD17_SED:',1, 1)
            self.AB_SED   = parseLines(Lines,'AB_SED:',  1, 1)
            
            # wavelength range
            self.LAMRANGE = parseLines(Lines,'LAMBDA_RANGE:',2, 1)
            
            # filter stuff
            self.MAGSYSTEM  = parseLines(Lines,'MAGSYSTEM:', 1, 1)
            self.FILTSYSTEM = parseLines(Lines,'FILTSYSTEM:',1, 1)
            self.FILTPATH   = parseLines(Lines,'FILTPATH:',  1, 1)
            self.FILTER     = parseLines(Lines,'FILTER:',    3, 1)
            self.NFILTER    = len(self.FILTER)

            if ( self.NFILTER > MXFILTERS ):
                errMsg = "%d filters exceeds bound of %d" % ( self.NFILTERS,MXFILTERS)
                sys.exit(errMsg)
                
            # define things to change later
            self.FILTER_NEWNAME = []     # to be changed later
            self.NEW            = 1      # assume new filter/mag system
            self.FILTER_CHARLIST_OLD   = ""
            self.FILTER_CHARLIST_NEW   = ""
            f.close()

def create_newVersion(versionInfo):

    VOUT = versionInfo.VERSION_OUT_TEXT  # name of output version

    # create new version-subDir
    if (len(VOUT) == 0 ):
        sys.exit("Output directory not defined\n Check SURVEY_OUT key")

    if ( os.path.exists(VOUT) ):
        shutil.rmtree(VOUT)
        
    os.mkdir(VOUT)

    # create auxillary files
    cdV = 'cd ' + VOUT + ' ; '

    cmd = cdV + ' touch ' + versionInfo.AUXFILE_README
    os.system(cmd)
    cmd = cdV + ' touch ' + versionInfo.AUXFILE_LIST
    os.system(cmd)
    cmd = cdV + ' touch ' + versionInfo.AUXFILE_IGNORE
    os.system(cmd)

    # write stuff in README
    README_OUTFILE = VOUT + '/' + versionInfo.AUXFILE_README
    f = open(README_OUTFILE,"wt")
    f.write("# Combined data files with command:\n")
    f.write("#    %s %s \n" % (sys.argv[0], sys.argv[1]) )
    f.write("# Host      = %s\n" % HOSTNAME )
    f.write("# User      = %s\n" % getpass.getuser()  )
    f.write("# Directory = %s\n" % CWD   )
    f.write("# Time      = %s\n" % NOW   )
    f.write("\n")
    f.close

    
def add_newVersion(VIN,versoinInfo,kcorInfo):

    # add new version to combined version;
    # copy all data files and use 'sed' utility
    # to make changes to filters and SURVEY name.
    
    dataDir   = TOPDIR_DATA + '/' + VIN
    SOUT      = versionInfo.SURVEY_OUT
    VOUT_TEXT = versionInfo.VERSION_OUT_TEXT
    
    # read list file
    LISTFILE_IN  = dataDir + '/' + VIN + '.LIST'
    LISTFILE_OUT = VOUT_TEXT + '/' + versionInfo.AUXFILE_LIST
    PTR_L        = open(LISTFILE_IN,"rt")
    fileList     = PTR_L.readlines()
    PTR_L.close

    # read contents of first file
    first_fileName = dataDir + '/' + fileList[0]
    first_fileName = first_fileName.replace("\n", "")
    f0    = open(first_fileName,"rt")
    fileContents = np.array(f0.readlines())
    f0.close
    
    # read FILTER string from first data file
    FILTERSTRING_OLD = parseLines(fileContents, 'FILTERS:', 1, 0)

    # get full filter lists from kcor file
    FILTERLIST_OLD = kcorInfo.FILTER_CHARLIST_OLD     # only this version
    FILTERLIST_NEW = kcorInfo.FILTER_CHARLIST_NEW     # only this version
    FILTERLIST_ALL = versionInfo.FILTER_CHARLIST_NEW  # all filters
    NFILTER        = kcorInfo.NFILTER
    print '\t ', FILTERLIST_OLD, ' -> ', FILTERLIST_NEW

    # read name of survey from first file
    SURVEY = parseLines(fileContents, 'SURVEY:', 1, 0)
    #print '\t  SURVEY_NAME = ', SURVEY
    
    # open new list file in append mode
    PTR_NEWLIST = open(LISTFILE_OUT,"at")

    # - - - - - start constructino of 'sed' command - - - - - - -
    sedcmd = "sed "
    
    # replace SURVEY
    sedAdd  = "-e 's/%s/ %s(%s)/g' " % (SURVEY,SOUT,SURVEY) 
    sedcmd += sedAdd
    
    # Replace global filter string.
    OLD  = FILTERSTRING_OLD  # from data file
    OLD2 = FILTERLIST_OLD    # from kcor file
    NEW  = FILTERLIST_NEW    # new kcor list
    ALL  = FILTERLIST_ALL    # all filters from all files
    sedAdd  = "-e 's/%s/%s  # %s -> %s/g' " % ( OLD, ALL, OLD2, NEW ) 
    sedcmd += sedAdd

    # Replace each single-char band.
    # Add '??' before each band to avoid removing new band
    # that matches an old band. Then remove ?? separately.

    for i in range(NFILTER):
        old = FILTERLIST_OLD[i]
        new = FILTERLIST_NEW[i]
        sedAdd = "-e 's/ %c / ??%c /g' " % ( old, new )
        sedcmd += sedAdd
        
    # remove the temporary '?s?'
    sedcmd += "-e 's/??//g' "
    
    # loop over all files and run 'sedcmd' 
    nfile=0
    for fin in fileList:
        fin    = fin.replace("\n", " ")
        FIN    = dataDir   + '/' + fin
        FOUT   = VOUT_TEXT + '/' + fin
        SEDCMD = sedcmd + FIN + ' > ' + FOUT        
        os.system(SEDCMD)
        PTR_NEWLIST.write("%s\n" % (fin) )
        nfile += 1

    PTR_NEWLIST.close

    # update README file
    README_OUTFILE = VOUT_TEXT + '/' + versionInfo.AUXFILE_README
    PTR_README = open(README_OUTFILE,"at")
    txt1 = "%3d data files from %-28.28s" % (nfile,VIN)
    txt2 = "%s -> %s" % (FILTERLIST_OLD,FILTERLIST_NEW)
    PTR_README.write( " %s %s \n" % (txt1,txt2) )
    PTR_README.close


def merge_duplicates(versionInfo):

    VOUT_TEXT  = versionInfo.VERSION_OUT_TEXT 

    # get list of SNID using grep :  file  SNID
    cmd = 'cd %s ; grep SNID: *.*' % VOUT_TEXT
    import subprocess
    import shutil
    output = subprocess.check_output(cmd, shell=True)
    import logging
    logging.basicConfig(level=logging.DEBUG, format="[%(funcName)20s()] %(message)s")
    output = output.split('\n')[:-1]
    d = {}
    for line in output:
        s = line.split()
        key = s[-1]
        val = s[0][:s[0].index("SNID:")-1]
        if key not in d:
            d[key] = []
        d[key].append(val)
    
    duplicates = [k for k in d.keys() if len(d[k]) > 1]
    logging.warn("Have %d duplicates for %s" % (len(duplicates), duplicates))
    for sn in duplicates:
        logging.warn("    %s %s" % (sn, d[sn])) 

    dup_dir = VOUT_TEXT + "/DUPLICATES"
    logging.info("Merging into %s" % dup_dir)
    if not os.path.exists(dup_dir):
        logging.debug("Creating directory %s" % dup_dir)
        os.makedirs(dup_dir)
    
    to_add = []
    for sn in duplicates:
        files = ["%s/%s" % (VOUT_TEXT, f) for f in d[sn]]
        
        output_file = "%s/merged_%s.dat" % (VOUT_TEXT, sn)
        to_add.append(os.path.basename(output_file) +"\n")
        buf = [] 
        with open(files[0]) as scaffold:
            for line in scaffold:
                if line.startswith("END:"):
                    continue
                buf.append(line)
        for f in files[1:]:
            with open(f) as obs:
                for line in obs:
                    if not line.startswith("OBS:"):
                        continue
                    buf.append(line)
        buf.append("END:") 
        n_obs = len([None for l in buf if l.startswith("OBS")])
        logging.debug("Merging to %s with %d obs" % (output_file, n_obs))
        for i, l in enumerate(buf):
            if l.startswith("NOBS:"):
                buf[i] = "NOBS: %d\n" % n_obs
        with open(output_file, "w") as output:
            for line in buf:
                output.write(line)
        for f in files:
            shutil.move(f, "%s/%s" % (dup_dir, os.path.basename(f)))

    # Update list
    list_file = "%s/%s.LIST" % (VOUT_TEXT, VOUT_TEXT)
    with open(list_file) as f:
        list_sn = [l for l in f]
    to_remove = [os.path.basename(filename) for sn in duplicates for filename in d[sn]]
    final_list = to_add + [f for f in list_sn if "".join(f.split()) not in to_remove]
    with open(list_file, "w") as f:
        f.writelines(final_list)
    logging.info("Updated list file %s" % list_file) 


    # For each pair of duplicates:
    # Get list of duplicates : filename & SNID
    # Merge "OBS:" lines into one file, keeping header of first file.
    # Update NOBS key-value and make sure that "END:" is at the end
    # Move duplicates into /DUPLICATES subDir
    

def add_vpec(versionInfo):

    # Created Jan 2018 by RK
    # call external script to read peculiar velocity correction 
    # (VPEC) and its error (VPEC_ERR) from a fitres-formatted 
    # text file, and  append each data file header. 
    # Use CLOBBER mode so that we don't have to move or copy 
    # anything when it's done.
    # First argument is ./VERSION (instead of just VERSION)
    # so that the script will see a slash and use this 
    # directory instead of searching under $SNDATA_ROOT/lcmerge.

    VPEC_FILE  = versionInfo.VPEC_FILE

    if ( len(VPEC_FILE) > 1 ):
        print "\n Add VPEC and VPEC_ERR from %s\n" % VPEC_FILE
        VOUT_TEXT  = versionInfo.VERSION_OUT_TEXT
        LOGFILE    = "ADD_VPEC.LOG"
        cmd  = "update_data_files.pl ./%s %s 'VPEC VPEC_ERR' CLOBBER > %s " % (VOUT_TEXT,VPEC_FILE,LOGFILE)
        os.system(cmd)




def make_simlib(versionInfo):

    print "\n Create SIMLIB "
    
    SOUT       = versionInfo.SURVEY_OUT
    VOUT_TEXT  = versionInfo.VERSION_OUT_TEXT
#    ZP_TEXT = versionInfo.ZPERR 
    nmlFile    = 'make_simlib_' + SOUT + '.nml'
    logFile    = 'make_simlib_' + SOUT + '.log'
    simlibFile = versionInfo.simlibFile

    # open with 1 line buffer so that it is ready to run snana.exe
    PTR_NML = open(nmlFile,"wt",1)

    PTR_NML.write(" &SNLCINP\n")
    PTR_NML.write("   VERSION_PHOTOMETRY = '%s' \n" % VOUT_TEXT )
    PTR_NML.write("   PRIVATE_DATA_PATH  = './' \n")
    PTR_NML.write("   KCOR_FILE    = '%s' \n" % versionInfo.kcor_outFile)

#    PTR_NML.write("   SIMLIB_ZPERR_LIST    = '%s' \n" % ZP_TEXT)
        
    
    PTR_NML.write("   SNTABLE_LIST = '' \n")
    PTR_NML.write("   SIMLIB_OUT   = '%s' \n" % simlibFile)
    
    PTR_NML.write(" \n")
    PTR_NML.write(" &END\n\n")
    PTR_NML.flush
    PTR_NML.close

    cmd  = "snana.exe %s > %s" % (nmlFile,logFile)
    os.system(cmd)

    
def convert2FITS(versionInfo):

    # convert TEXT format into FITS format
    S      = versionInfo.SURVEY_OUT
    V_TEXT = versionInfo.VERSION_OUT_TEXT
    V_FITS = versionInfo.VERSION_OUT_FITS

    print "\n Convert TEXT formatted data files into FITS: "
    print "\t %s -> %s\n" % (V_TEXT,V_FITS)
    
    # create new version directory
    if ( os.path.exists(V_FITS) ):
        shutil.rmtree(V_FITS)
        
    os.mkdir(V_FITS)

    # construct snana nmlFile
    nmlFile   = "convert2FITS_%s.nml" % S
    NMLFILE   = "%s/%s" % (V_FITS,nmlFile)
    logFile   = "convert2FITS_%s.log" % S

    PTR_NML   = open(NMLFILE,"wt",1)
    PTR_NML.write(" &SNLCINP \n")
    PTR_NML.write("   VERSION_PHOTOMETRY = '%s' \n" % V_TEXT )
    PTR_NML.write("   PRIVATE_DATA_PATH  = '%s' \n"  % CWD   )
    PTR_NML.write("   VERSION_REFORMAT_FITS = '%s' \n" % V_FITS )
    PTR_NML.write(" &END \n\n")
    PTR_NML.flush
    PTR_NML.close
    cmd = "cd %s ; snana.exe %s > %s " % (V_FITS,nmlFile, logFile)
    os.system(cmd)

    
def printSummary(versionInfo):
    print '\n# =============================================== '
    print ' Summary of outputs: '

    V = versionInfo.VERSION_OUT_TEXT
    print "\t %s/    (combined Data version)" % V

    V = versionInfo.VERSION_OUT_FITS
    print "\t %s/    (combined Data version)" % V

    file = versionInfo.kcor_outFile
    print "\t %s    (combined kcor file)" % file

    file = versionInfo.simlibFile
    print "\t %s    (combined SIMLIB file)" % file
    
    
# =========================
# ======= MAIN ============
# =========================


if __name__ == "__main__":

# parse input argument(s)
    if ( len(sys.argv) < 2 ):
        sys.exit("Must give INFILE arguent\n-->ABORT")
    else:
        INFILE = sys.argv[1]
        print 'Input file: ', INFILE

    print 'SNDATA_ROOT = ', SNDATA_ROOT 

    versionInfo = VERSION_INFO(INFILE)

    # ----------------------------------------
    # prepare kcor file with all filters
    # ----------------------------------------
    kcorInfo = []
    nkcor=0
    for kcorFile in versionInfo.INFILE_KCOR:
        print ' '
        nkcor += 1
        kcorInfo.append(KCOR_INFO(kcorFile))

    print "\n Done parsing ", nkcor, " kcor-input files "

    # change filter char
    change_filterChar(versionInfo,kcorInfo)
            
    # write new combined kcor-input file
    write_kcor_inputFile(versionInfo,kcorInfo) 
    run_kcor(versionInfo)
    
    # ---------------------------------------------------------
    # now re-write data files with different filter strings
    # ---------------------------------------------------------

    print "\n# - - - - - - - - - - - - - - - - - - - - - - "
    
    # check for private data path
    if ( len(versionInfo.PRIVATE_DATA_PATH) > 0 ):
        global TOPDIR_DATA
        TOPDIR_DATA = versionInfo.PRIVATE_DATA_PATH
        print " Use PRIVATE_DATA_PATH = ", TOPDIR_DATA

    create_newVersion(versionInfo)  
    nver=0
    for vname in versionInfo.NAME :
        print ' Swap filter strings in DATA-VERSION: ', vname
        add_newVersion(vname,versionInfo,kcorInfo[nver])
        nver += 1

    # -----------------
    # merge duplicate light curves from different instruments (Dec 2017)
    merge_duplicates(versionInfo)

    # --------------------
    # check option to add VPEC & VPEC_ERR (Jan 2018)
    add_vpec(versionInfo)
    
    # -----------------
    make_simlib(versionInfo)

    # -----------------
    convert2FITS(versionInfo)
    
    # --------------------
    # print summary of outputs
    printSummary(versionInfo)
    
    
# ========= END MAIN ================

