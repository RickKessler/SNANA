# =========================================================
# Created April 2020 by R.Kessler (and help from D.Jones)
# Translate CSPDR3 data files into SNANA format.
#
# Usage:
#   python translate_CSPDR3.py
#
# CSP data directory must exist as
#   DR3/
# that is downloaded and untarred from
#   https://csp.obs.carnegiescience.edu/data
#     [click  CSP_Photomery_DR3.tar.gz]
#
# CSP DR3 does not include LOGMASS values, but user-computed
# LOGMASS vales can be read [optionally] from CSPDR3_LOGMASS.DAT
# if this file exists with following format:
#    SN: <SNID> <LOGMASS> <LOGMASS_ERR>
#
# Output SNANA directory is CSPDR3/
#
# =========================================================

import os
import sys 
import time
import glob
import gzip
import shutil
import datetime
import numpy as np

# get SFD98 map
from astropy.coordinates import SkyCoord
import astropy.units as u
from dustmaps.sfd import SFDQuery
sfd = SFDQuery()

INPDIR             = 'DR3'
VERSION_SNANA      = 'CSPDR3'  # output dir name
LOGMASS_FILE       = 'CSPDR3_LOGMASS.DAT'
zERR               = 0.0001    # no z-errors in DR, so this is a fudge
MJDOFF             = 53000.0   # add this back to CSP dates

# define auxilary SNANA files
VV                = ( "%s/%s" % (VERSION_SNANA,VERSION_SNANA) )
SNANA_IGNORE_FILE = ( "%s.IGNORE" %  VV )
SNANA_README_FILE = ( "%s.README" %  VV )
SNANA_LIST_FILE   = ( "%s.LIST"   %  VV )

ZP_SNANA      = 27.5        # FLUXCAL = 10^[-0.4*MAG-ZP_SNANA]

# misc.
HOSTNAME      = os.environ['HOSTNAME']
CWD           = os.getcwd()
SCRIPTNAME    = sys.argv[0]
tnow          = datetime.datetime.now()
TSTAMP        = ('%4.4d-%2.2d-%2.2d' %  (tnow.year,tnow.month,tnow.day) )


# ========================================================
# ========================================================
def makeFilterMap():

    # define map between CSP defintions (left) and SNANA charName (right)
    FILTERMAP_SNANA = {
        "u"    : "u" ,
        "g"    : "g" ,
        "r"    : "r" ,
        "i"    : "i" ,
        "B"    : "B" ,
        "V"    : "n" ,
        "V0"   : "m" ,
        "V1"   : "o" ,
        "J"    : "J" ,
        "Jrc2" : "j" ,
        "Y"    : "Y" ,
        "Ydw"  : "y" ,
        "H"    : "H"
    }

    FLIST = ''
    for f in FILTERMAP_SNANA:
        FLIST += FILTERMAP_SNANA[f]

    print(' Init SNANA Filter list: %s ' % FLIST )
    FILTERMAP_SNANA['ALL'] = FLIST
    return FILTERMAP_SNANA

def znew(ra, dec, z):
    c_icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame = 'icrs')
    c_icrs = c_icrs.galactic
    b = c_icrs.b.degree
    l = c_icrs.l.degree
    b = np.radians(b)
    l = np.radians(l)
    l_0 = np.radians(264.14)
    b_0 = np.radians(48.26)
    v = (float(z)*3*10**5 + 371 * (np.sin(b) * np.sin(b_0) + np.cos(b) * np.cos(b_0) * np.cos(l-l_0)))/(3*10**5)
    return v

# ====================================
def get_SNFILE_LIST():

     LIST = []
#     os.chdir(INPDIR)
     for file in glob.glob("DR3/SN*snpy.txt"):
     #     print(" file = '%s' " % file )
          LIST.append(file)

     return LIST

# ==================================
def create_OUTDIR():
     # create out dir for SNANA-formatted files.

     if ( os.path.exists(VERSION_SNANA) ):
          shutil.rmtree(VERSION_SNANA)

     print(" Create %s" % VERSION_SNANA)
     os.mkdir(VERSION_SNANA)

# ===================================
def translate_file(INPFILE, FILTERMAP_SNANA, LOGMASS_LIST):

# read input file from CSP and create output SNANA-formatted file.
# Note that INPFILE includes DR3/ path.

# get SN name between SN and _snpy
     jstart  = INPFILE.find("SN") + 2
     jend    = INPFILE.find("_snpy")
     SNID    = INPFILE[jstart:jend]

# construct file names: local, and OUTFILE that includes path
     outFile   = ( "%s_%s.DAT" % (VERSION_SNANA,SNID) )
     OUTFILE   = ( "%s/%s"     % (VERSION_SNANA,outFile) )
     OUTFILEgz = ( "%s.gz"     % OUTFILE )

# open output SNANA-formatted file
#     fout = open(OUTFILE,"wt")
     fout = gzip.open(OUTFILEgz,'wt')

     print(" translate file = %s -> %s  " % (INPFILE,OUTFILEgz) )

     with open(INPFILE, 'rt') as finp:
          lineList = finp.readlines()    

     NLINE=0;  MJD_LIST = [];  MAG_LIST=[]; MAGERR_LIST=[]; BAND_LIST=[]
     PEAKMJD = -9.0; PEAKMAG = 99.0

     for line in lineList:          
          NLINE+=1 
          line  = line.rstrip()  # remove trailing space and linefeed
          words = line.split()
          if ( NLINE == 1 ):
               zHELIO = float(words[1])
               RA     = float(words[2])
               DEC    = float(words[3])

               sc = SkyCoord(RA,DEC,frame="fk5",unit=u.deg) # for decimal deg
               SFD98     = sfd(sc)
               MWEBV     = SFD98*0.86  # include SF11 corr
               MWEBVERR  = MWEBV/6.0
               zCMB      = znew(RA,DEC,zHELIO)

          elif ( words[0] == 'filter' ) :
               band_CSP   = words[1]
               if ( band_CSP not in FILTERMAP_SNANA ):
                   sys.exit('\n ERROR: unknown CSP band = %s \n' % band_CSP)
               band_SNANA = FILTERMAP_SNANA[band_CSP] 
          else :
               mjd = float(words[0])
               mag = float(words[1])
               magerr = float(words[2])
               MJD_LIST.append(mjd)
               MAG_LIST.append(mag)
               MAGERR_LIST.append(magerr)
               BAND_LIST.append(band_SNANA)
               if ( mag < PEAKMAG ) :
                    PEAKMAG = mag
                    PEAKMJD = mjd

# write to output file
     fout.write('SURVEY:   CSP \n')
     fout.write('SNID:     %s \n' % SNID)
     fout.write('FILTERS:  %s \n' % FILTERMAP_SNANA["ALL"] )
     fout.write('RA:       %10.6f  # deg \n' % RA)
     fout.write('DEC:      %10.6f  # deg \n' % DEC)
     fout.write('REDSHIFT_HELIO: %6.4f +- %6.4f  # Helio \n'
                % (zHELIO,zERR) )
     fout.write('REDSHIFT_FINAL: %6.4f +- %6.4f  # CMB \n'
                % (zCMB,zERR) )
     fout.write('MWEBV:    %6.3f +- %6.3f        # SFD98 * 0.86 \n' 
                % (MWEBV,MWEBVERR) )

     PEAKMJD += MJDOFF
     fout.write('PEAKMJD:  %8.2f                # at brightest obs \n' % PEAKMJD )

     # check for logmass value
     if SNID in LOGMASS_LIST :
         LOGMASS_STRING = LOGMASS_LIST[SNID]
         LOGMASS_LIST["NUSE"]+=1
     else:
         LOGMASS_STRING = "10.00 +- 3.00"
     fout.write('HOSTGAL_LOGMASS:  %s  #  \n' % LOGMASS_STRING )


     NOBS    = len(MJD_LIST)
     NVAR    = 7
     VARLIST = 'MJD  FLT FIELD   FLUXCAL   FLUXCALERR    MAG   MAGERR'

     fout.write('\n# ----------------------------------- \n')
     fout.write('NOBS: %d \n' % NOBS)
     fout.write('NVAR: %d \n' % NVAR)
     fout.write('VARLIST: %s \n' % VARLIST )
     
     for i in range(0, NOBS):
          mjd     = MJD_LIST[i] + MJDOFF
          mag     = MAG_LIST[i]
          magerr  = MAGERR_LIST[i]
          band    = BAND_LIST[i]
          arg0    = -0.4*(mag-ZP_SNANA)
          fluxcal = pow(10.0,arg0)
          arg1    = 0.4*magerr
          fluxcalerr = fluxcal * (pow(10.0,arg1)-1.0)
          fout.write('OBS: %9.3f  %s NULL  %10.5e %10.5e   %6.3f %6.3f \n' % 
                     (mjd,band,fluxcal, fluxcalerr, mag,magerr) )
          
     fout.write('END: \n')
     fout.close
     finp.close

     return outFile

# =========================================
def makeAuxFile_LIST(OUTFILE_LIST):

     f = open(SNANA_LIST_FILE,"wt")

     for file in OUTFILE_LIST :
          f.write('%s\n' % file)
     f.close

def makeAuxFile_README(OUTFILE_LIST,LOGMASS_LIST,FILTERMAP_SNANA):
     f = open(SNANA_README_FILE,"wt")

     NFILE = len(OUTFILE_LIST)
     f.write(' %s \n' % TSTAMP)
     f.write(' Translated %d CSP data files from \n' % NFILE)
     f.write('   https://csp.obs.carnegiescience.edu/data \n')
     f.write(' to SNANA format using %s\n' % SCRIPTNAME )
     f.write('\n')
     f.write(' If using these data, please reference\n')
     f.write('   https://ui.adsabs.harvard.edu/abs/2017AJ....154..211K\n')
     f.write('\n')
     f.write(' Extra quantities in SNANA data files that are not in DR3:\n')
     f.write('   REDSHIFT_FINAL computed from CMB dipole.\n')
     f.write('   zERR    = %.4f \n' % zERR)
     f.write('   PEAKMJD is from brightest observation.\n')
     f.write('   MWEBV   is from SFD98 x 0.86 correction from SF2011. \n' )

     if ( "COMMENT" in LOGMASS_LIST ) :
         f.write('   LOGMASS is from %s\n' % LOGMASS_LIST["COMMENT"] )

     f.write('\n')
     f.write('FILTER MAP: \n')
     for f_CSP in FILTERMAP_SNANA :
         f_SNANA = FILTERMAP_SNANA[f_CSP]
         if ( f_CSP != "ALL" ) :
             f.write('   SNANA %s -> CSP %s \n' % (f_SNANA,f_CSP) )

     f.write('\nEND:\n')
     f.close

def makeAuxFile_IGNORE():
     f = open(SNANA_IGNORE_FILE,"wt")
     f.close

# ========================================
def read_LOGMASS():

    LOGMASS_LIST = {}
    if ( not os.path.isfile(LOGMASS_FILE) ) :
        print(' Cannot find %s -> skip logmass \n' % LOGMASS_FILE )
        return LOGMASS_LIST

    with open(LOGMASS_FILE, 'rt') as f:
        lineList = f.readlines()    
    
    NRD = 0
    for line in lineList:
        line  = line.rstrip()  # remove trailing space and linefeed
        words = line.split()
        if ( len(words) == 0 ) :
            continue
        KEY   = words[0]
        if ( KEY == "COMMENT:" ):
            STRING = words[1:]
            LOGMASS_LIST["COMMENT"] = line[8:]
        if ( KEY == 'SN:' ):
            SNID        = words[1]
            LOGMASS     = float(words[2])
            LOGMASS_ERR = float(words[3])
            STRING      = ( "%5.2f +- %5.2f" % (LOGMASS,LOGMASS_ERR) )
            LOGMASS_LIST[SNID] = STRING
            NRD+=1

    f.close

    LOGMASS_LIST["NREAD"] = NRD;
    LOGMASS_LIST["NUSE"]  = 0 ;
    print(' Read %d LOGMASS values from  %s \n' % (NRD,LOGMASS_FILE) )    
    return LOGMASS_LIST

# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":

      # make filter list for data file header
      FILTERMAP_SNANA = makeFilterMap()

      SNFILE_LIST_INP = get_SNFILE_LIST()
      NFILE = len(SNFILE_LIST_INP)

      create_OUTDIR()


      LOGMASS_LIST = read_LOGMASS()

      NPROC = 0
      OUTFILE_LIST = []
      for file in SNFILE_LIST_INP :
           outFile = translate_file(file,FILTERMAP_SNANA,LOGMASS_LIST)
           OUTFILE_LIST.append(outFile)
           NPROC+=1;

      makeAuxFile_LIST(OUTFILE_LIST)
      makeAuxFile_README(OUTFILE_LIST,LOGMASS_LIST,FILTERMAP_SNANA)
      makeAuxFile_IGNORE()

      print("\n Done translating %d data files. " % NFILE)
      print(' SNANA Filter list: %s ' % FILTERMAP_SNANA["ALL"] )
      print(' %d LOGMASS values included.' % LOGMASS_LIST["NUSE"])
      print(' See %s for more info.' % SNANA_README_FILE)

# ==== END ====

