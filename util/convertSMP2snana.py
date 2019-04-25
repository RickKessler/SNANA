#!/usr/bin/env python
#
# Created Mar 2017 by R.Kessler & D.Brout
#
# Read output from Scene-Model-Photometry (SMP) and
# convert into SNANA-ascii-formatted light curves.
# Also run converter to make FITS-formatted SNANA files.
#
# Usage:
#   convertSMP2snana.py inpDir outPrefix
#
#  If inpDir = /project/data/bla/mylc, then the mylc.LIST
#  must exist inside the inpDir, and mylc.LIST contains a
#  list of SMP data files.
#  An optional file, 'mylc.HEADER' , contains additional
#  HEADER variables to include in the output SNANA data files.
#  outPrefix is the output directory prefix, and is also the
#  prefix of the SNANA version. The outDir
#  names are
#    [outPrefix]_snana_text
#    [outPrefix]_snana_fits
#
# DEBUG command:
# convertSMP2snana.py SMP_RAW_SPECTYPE_v1_1 DESALL_specType_SMP_real
#
# ============================================================

import os
import sys 
import numpy as np
import getpass
import time
import shutil
from copy   import copy
from shutil import copyfile

# globals

SNDATA_ROOT = os.environ['SNDATA_ROOT']
HOSTNAME    = os.environ['HOSTNAME']
NOW         = time.strftime("%c")
CWD         = os.getcwd() 
ZP_SNANA    = 27.5
SMP_FLAG_REJECT = 32768 + 65536  # reject exposure if this bit is set

SUFFIX_TEXT = "_snana_text";
SUFFIX_FITS = "_snana_fits";

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
            
            
def PREPARE_DATAFILE(fileName,SMP_INFO,HEADER_INFO,SNANA_INFO):

    # construct name of input SMP file, and output SNANA file.
    INFILE  = "%s/%s" % (SMP_INFO.INPDIR,fileName)    
    OUTFILE = "%s/%s" % (SNANA_INFO.VERSION_TEXT,fileName)
    outFile = fileName
    
    # read entire SMP file into memory
    fin = open(INFILE,"rt")
    Lines  = np.array(fin.readlines())
    fin.close();
        
    # get SNID & number of filters
    VALUE = parseLines(Lines, "FILTERS:", 1, 0 )
    NFILTER = len(VALUE)
    
    SNID   = parseLines(Lines, "SNID:",   1, 0 )
    SURVEY = parseLines(Lines, "SURVEY:", 1, 0 )
    
    # get header row index for this SNID
    try:
        iHead = HEADER_INFO.CID.tolist().index(SNID)
        SNANA_INFO.NFILE_CREATE += 1
        SNANA_INFO.OUTFILENAME_LIST += ("%s\n" % outFile)
    except:
        iHead = -1   # not in header file
        SNANA_INFO.NFILE_SKIP += 1
        SNANA_INFO.SNID_MISSING += ("%s " % SNID )
        return
    
    # open output file 
    fout  = open(OUTFILE,"wt") 
    found = copyHeaderKey(Lines, 1,"SURVEY:",   HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"SNID:",     HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"FILTERS:",  HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"FAKE:",     HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"RA:",       HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"DECL:",     HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"TYPE:",     HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"PEAKMJD:",  HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"PSF_UNIT:", HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 3,"REDSHIFT_HELIO:", HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 3,"REDSHIFT_FINAL:", HEADER_INFO,iHead,fout)
            
    fout.write("MWEBV:       0.0   # WARNING: see OPT_MWEBV in manual\n")
    fout.write("MWEBV_ERR:   0.0 \n")
    fout.write("\n")
    
    found = copyHeaderKey(Lines, 1,"HOSTGAL_OBJID:",    HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 3,"HOSTGAL_PHOTOZ:",   HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 3,"HOSTGAL_SPECZ:",    HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 1,"HOSTGAL_SNSEP:",    HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, 3,"HOSTGAL_LOGMASS:",  HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, NFILTER,"HOSTGAL_MAG:", HEADER_INFO,iHead,fout)
    found = copyHeaderKey(Lines, NFILTER,"HOSTGAL_SB_FLUXCAL:", HEADER_INFO,iHead,fout)

    fout.write("\n\n# Appended from HEADER-override file. \n")
    # write out unused HEADER variables
    HEADER_INFO.CID[iHead] = 'USED'
    for keyName in HEADER_INFO.__dict__.keys():
        value = HEADER_INFO.__dict__[keyName][iHead]
        SKIP  = 0
        if ( value == 'USED' or keyName.find('_ERR') > 0 ):
            SKIP = 1
        if ( SKIP == 0 ):
            KEY = ("PRIVATE(%s):" % keyName)
            fout.write("%-30.30s   %s \n" % (KEY,value) )
    
    # - - - - - - - - - - - - - - - - - -  -
    fout.write("\n# ---------------------- \n")

    # - - - - - - - -
    # read and write epochs
    write_SNANA_OBS(SURVEY,SNID,INFILE,fout)
    
    
    # close SNANA output file
    fout.close()
    return 
    
def copyHeaderKey(Lines,narg,keyName,HEADER_INFO, iHead, fout):
    
    # always check optional header file to override
    # whatever is in the SMP file
    if ( iHead >= 0 ):
        found = checkHeaderFile(keyName, narg, HEADER_INFO, iHead, fout)
        if ( found == 1 ):
            return(1)
        
    # read narg values after keyname
    VALUE = parseLines(Lines, keyName, narg, 0 )
    
    if ( len(VALUE) == 0 ):
        return(0)

    outLine = "%-20.20s" % (keyName)
    if ( narg == 1 ):
        outLine = ("%s %s" % (outLine,VALUE) )
    else:
        for i in range(0,narg):
            val     = VALUE[i]
            outLine = ("%s %s" % (outLine,val) )
        
    fout.write("%s \n" % (outLine ) )
    return(1)


def checkHeaderFile(keyName, narg, HEADER_INFO, iHead, fout):

    # if keyName is specified in HEADER_INFO, then
    # write it in SNANA data (fout) file and ignore the
    # value in the input SMP data file.
    # If narg==3, then check for keyName_ERR
    # Input iHead is the row in header file.
    #
    # Returns 1 if keyName is found and written; return 0 otherwise
    
    key    = keyName.replace(":","")  # remove colon
    keyErr = key + "_ERR"             # in case narg==3
    
    if key in HEADER_INFO.__dict__.keys():
        value = HEADER_INFO.__dict__[key][iHead]
        VALUE = value
        HEADER_INFO.__dict__[key][iHead] = 'USED'
    else:
        return(0)


    if ( narg == 3 and keyErr in HEADER_INFO.__dict__.keys() ):
        err = HEADER_INFO.__dict__[keyErr][iHead]
        VALUE = "%s +- %s" % ( value,err )
        HEADER_INFO.__dict__[keyErr][iHead] = 'USED'
    
    fout.write("%-20.20s  %s    #  OVERRIDE \n" % (keyName,VALUE) )
    
    return(1)

def write_SNANA_OBS(SURVEY,SNID,INFILE,fout):

    VARLIST = [ "MJD", "BAND", "FIELD", "FLUXCAL", "FLUXCALERR",
                "ZPFLUX", "PSF", "SKYSIG", "PHOTFLAG", "PHOTPROB" ]
    
    NVAR    = len(VARLIST)

    
    OBSTABLE = READ_TABLE( "OBS:", INFILE, 0)

    # convert string arrays into int and float
    OBSTABLE.ID_OBS     = np.array(OBSTABLE.ID_OBS,dtype='int')
    OBSTABLE.ID_COADD   = np.array(OBSTABLE.ID_COADD,dtype='int')
    OBSTABLE.FLUXERR    = np.array(OBSTABLE.FLUXERR,dtype='float')
    OBSTABLE.MJD        = np.array(OBSTABLE.MJD,dtype='float')
    OBSTABLE.ZP         = np.array(OBSTABLE.ZP,dtype='float')
    OBSTABLE.PSF        = np.array(OBSTABLE.PSF,dtype='float')
    OBSTABLE.SKYSIG     = np.array(OBSTABLE.SKYSIG,dtype='float')
    
    OBSTABLE.SMP_FLAG    = np.array(OBSTABLE.SMP_FLAG,dtype='int')
    OBSTABLE.SMP_CHISQ   = np.array(OBSTABLE.SMP_CHISQ,dtype='float')
    OBSTABLE.SMP_FLUX    = np.array(OBSTABLE.SMP_FLUX,dtype='float')
    OBSTABLE.SMP_FLUXERR = np.array(OBSTABLE.SMP_FLUXERR,dtype='float')
    OBSTABLE.SMP_FLUX_ZPT = np.array(OBSTABLE.SMP_FLUX_ZPT,dtype='float')

    # avoid float overflow when ZPT value is crazy (or any
    SMP_FLUX_ZPT = copy(OBSTABLE.SMP_FLUX_ZPT)
    SMP_FLUX_ZPT[(SMP_FLUX_ZPT<0.0) | (SMP_FLUX_ZPT>80.)] = ZP_SNANA
    
    # get weight for wgted avg. If SMP_FLAG<0, set WGT=0
    # so that bad epochs are excluded
    WGT  = 1.0/OBSTABLE.FLUXERR**2
    WGT[(OBSTABLE.SMP_FLAG<0.0)] = 0.0
    WGT[(OBSTABLE.SMP_FLAG & SMP_FLAG_REJECT)>0] = 0.0

    # compute SNANA FLUXCAL[ERR] using SNANA ZP_SNANA
    ZPDIF      = SMP_FLUX_ZPT - ZP_SNANA
    arg        = 10.0**(-0.4*ZPDIF)
    FLUXCAL    = OBSTABLE.SMP_FLUX * arg
    FLUXCALERR = OBSTABLE.SMP_FLUXERR * arg

    NOBS = 0
    OBSLINE_ALL = ""
    for ID_COADD in np.unique(OBSTABLE.ID_COADD):
        ww     = (OBSTABLE.ID_COADD == ID_COADD)  # ID_OBS mask
        WGTSUM = np.sum(WGT[ww])
        if ( WGTSUM == 0.0 ):  continue
                
        MJD     = np.sum(WGT[ww]*OBSTABLE.MJD[ww])/WGTSUM
        BAND    = OBSTABLE.BAND[ww][0]

        try:
            FIELD   = OBSTABLE.FIELD[ww][0]
        except:
            IMAGE_NAME = OBSTABLE.IMAGE_NAME_SEARCH[ww][0]
            FIELD      = getField(SURVEY,IMAGE_NAME)
                
        F       = np.sum(WGT[ww]*FLUXCAL[ww])/WGTSUM
        FERR    = np.sum(WGT[ww]*FLUXCALERR[ww])/WGTSUM

        if ( FERR   == 0.0 ):  continue
                
        ZP      = np.sum(WGT[ww]*OBSTABLE.ZP[ww])/WGTSUM
        PSF     = np.sum(WGT[ww]*OBSTABLE.PSF[ww])/WGTSUM
        SKYSIG  = np.sum(WGT[ww]*OBSTABLE.SKYSIG[ww])/WGTSUM

        PHOTFLAG=0
        for flag in OBSTABLE.SMP_FLAG[ww]:  PHOTFLAG |= flag
            
        PHOTPROB   = np.sum(WGT[ww]*OBSTABLE.SMP_CHISQ[ww])/WGTSUM
        
        STR_FLUX   = "%8.3f %8.3f " % (F,FERR)
        STR_OBS    = "%6.3f %5.3f %6.2f " % (ZP,PSF,SKYSIG)
        STR_PHOT   = "%5d %.2f " % (PHOTFLAG,PHOTPROB)
        STR_ALL    = STR_FLUX + STR_OBS + STR_PHOT
        OBSLINE_ALL += "OBS: %.3f %s %s %s\n" % (MJD,BAND,FIELD,STR_ALL)
        NOBS += 1

    # --------------
    # write NOBS key and then all OBS lines
    fout.write("NOBS: %d\n" % (NOBS) )
    fout.write("NVAR: %d\n" % NVAR)
    line = "VARLIST: "
    for  var in VARLIST:  line = ( "%s %s" % (line,var) )
    fout.write("%s\n" % line);

    fout.write("%s" % (OBSLINE_ALL) )
    return

#def sdfs(ID_COADD
    

         
def write_auxFiles(SMP_INFO):

    V = SNANA_INFO.VERSION_TEXT
    
    # copy list file xyz
    TMPFILE = "%s/%s.LIST" % (V,V)
    f = open(TMPFILE,"wt")
    f.write(SNANA_INFO.OUTFILENAME_LIST)
    f.close()

    # create blank IGNORE file
    TMPFILE = "%s/%s.IGNORE" % (V,V)
    f = open(TMPFILE,"w")
    f.close()

    # create README and add info at top of it
    f = open(SMP_INFO.DOCFILE,'rt')
    temp = f.read()
    f.close()

    TMPFILE = "%s/%s.README" % (V,V)
    f = open(TMPFILE, 'w')
    cmd = "%s %s %s" % (sys.argv[0],sys.argv[1],sys.argv[2])
    f.write("Created with command \n  %s\n" % cmd )
    f.write("HOST: %s\n" % HOSTNAME)
    f.write("TIME: %s\n" % NOW )
    f.write("\n ------------- SMP README ------------- \n")
    f.write(temp)
    f.close()
    
    return


def getField(SURVEY,IMAGE_NAME):

    # kludge function to extract FIELD from IMAGE_NAME;
    # called only if FIELD is not in SMP file.
    
    if SURVEY == 'DES':
        # grab 2 characters after 'SN-'
        FIELD = IMAGE_NAME.split('SN-')[1][:2]
    else:
        FIELD = 'UNKNOWN'
        
    return(FIELD)

def convert_text_to_fits(SNANA_INFO):

    # convert SNANA TEXT files into FITS format.
    
    V_TEXT = SNANA_INFO.VERSION_TEXT  # input TEXT version
    V_FITS = SNANA_INFO.VERSION_FITS  # output FITS version
    
    # construct argList for snana.exe
    argVin  = "VERSION_PHOTOMETRY "    + V_TEXT
    argVout = "VERSION_REFORMAT_FITS " + V_FITS
    argPath = "PRIVATE_DATA_PATH ../ "
    argList = "NOFILE %s %s %s " % (argVin,argVout,argPath)
    cmd     = "cd %s ; snana.exe %s >& snana.log" % (V_FITS,argList)

    print "\n Creating FITS-formatted version ", V_FITS
    os.system(cmd)
    
class SMP_INPUT_INFO:
    # return INPUT_PREFIX, names of auxilary files, and name
    # of each data file
    def __init__(self,INPDIR):

        PREFIX = INPDIR

        LISTFILE = INPDIR + '/' + PREFIX + '.LIST'
        HEADFILE = INPDIR + '/' + PREFIX + '.HEADER'
        DOCFILE  = INPDIR + '/' + PREFIX + '.README'
        
        f = open(LISTFILE,"rt")
        Lines = np.array(f.readlines())
        f.close()

        NFILE = len(Lines)

        print " Found ", NFILE," SMP files in ", INPDIR
#        print Lines
        
        # check optional header file

        self.INPDIR   = INPDIR 
        self.PREFIX   = PREFIX
        self.LISTFILE = LISTFILE
        self.HEADFILE = HEADFILE
        self.DOCFILE  = DOCFILE
        self.NFILE    = NFILE
        self.FILE_LIST = Lines

        
class READ_TABLE:
    def __init__(self,LINEKEY,inFile,VBOSE):

        if ( os.path.exists(inFile) == 0 ):
            return

        if ( VBOSE > 0 ):
            print " Read tables rows starting with ", LINEKEY, "\n    from file: ", inFile

        # read SNANA-formatted file using code taken from
        # D.Jones' ovdatamc.py

        fin   = open(inFile,'rt')
        lines = fin.readlines()
        for l in lines:
            if l.startswith('VARNAMES:'):
                l = l.replace('\n','')
                coldefs = l.split()
                break
            
        with open(inFile) as f:
            reader = [x.split() for x in f if x.startswith(LINEKEY)]
                
        i = 0
        for column in zip(*reader):
            if ( coldefs[i] != 'VARNAMES:' ):  # skip VARNAMES
                self.__dict__[coldefs[i]] = np.array(column[:],dtype='object')
            i += 1
            
        fin.close()
        
class SNANA_OUTPUT_INFO:
    # return INPUT_PREFIX, names of auxilary files, and name
    # of each data file
    def __init__(self,OUTPUT_PREFIX):

        V_TEXT = OUTPUT_PREFIX + SUFFIX_TEXT
        V_FITS = OUTPUT_PREFIX + SUFFIX_FITS
        self.VERSION_TEXT = V_TEXT
        self.VERSION_FITS = V_FITS
        self.NFILE_CREATE = 0
        self.NFILE_SKIP   = 0
        self.SNID_MISSING = ""
        self.OUTFILENAME_LIST = ""
        
        print ""
        print "Create output SNANA directories: "
        print "   ", V_TEXT
        print "   ", V_FITS

        if ( os.path.exists(V_TEXT) ):
            shutil.rmtree(V_TEXT)

        if ( os.path.exists(V_FITS) ):
            shutil.rmtree(V_FITS)
                        
        os.mkdir(V_TEXT)
        os.mkdir(V_FITS)

        
# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":

# parse input argument(s)
    if ( len(sys.argv) < 3 ):
        sys.exit("Must give inpDir and outDir arguents \n-->ABORT")
    else:
        INPDIR         = sys.argv[1]
        OUTPUT_PREFIX  = sys.argv[2]
        print 'Input  SMP   Directory: ', INPDIR
        print 'output SNANA Prefix   : ', OUTPUT_PREFIX

    SMP_INFO    = SMP_INPUT_INFO(INPDIR)
    HEADER_INFO = READ_TABLE("SN:",SMP_INFO.HEADFILE,1) 
    SNANA_INFO  = SNANA_OUTPUT_INFO(OUTPUT_PREFIX)

    NFILE=0;
    for file in SMP_INFO.FILE_LIST:
        ff = file.replace('\n','')    # remove <CR>
        print "\t Process SMP file ", ff
        PREPARE_DATAFILE(ff,SMP_INFO,HEADER_INFO,SNANA_INFO)
        NFILE += 1 ;
#        sys.exit("\n xxx debug die xxx")

    write_auxFiles(SMP_INFO)
    print " SNANA TEXT files under ", SNANA_INFO.VERSION_TEXT
    
    # convert text file format into FITS format
    convert_text_to_fits(SNANA_INFO)
    
    print "%5d events created. " % SNANA_INFO.NFILE_CREATE
    print "%5d events skipped. " % SNANA_INFO.NFILE_SKIP

    if SNANA_INFO.NFILE_SKIP>0 :
        print "\n Missing SNID in ", SMP_INFO.HEADFILE
        print "\t ", SNANA_INFO.SNID_MISSING
        
# ========= END MAIN ===============

