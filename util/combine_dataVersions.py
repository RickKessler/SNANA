#!/usr/bin/env python
#
# Created Feb 2017 by R.Kessler
# Nov 10 2021: Major update/refactor (see history below)
#  
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
#   PRIVATE_DATA_PATH:  <path>   # optional path to data
#   KCOR_PATH:          <path>   # optional path to kcor files
#   SURVEY_OUT:  <name of combined survey>  # required
#   VPEC_FILE:   <name file file with VPEC & VPEC_ERR>
#   MAGSYSTEM_FINAL:    <system>  # needed if there are multiple systems
#   CHANGE_FILTER_CHAR: False     # optional preserve filter char names
#
#   VERSION:  <ver1>  <kcor_inFile1>
#   VERSION:  <ver2>  <kcor_inFile2>
#     etc ...
#
# Outputs:
#    <SURVEY_OUT>_TEXT/       ! combined data directory, TEXT format
#    <SURVEY_OUT>_FITS/       ! combined data directory, FITS format
#    kcor_<SURVEY_OUT>.fits   ! combined kcor file
#    <SURVEY_OUT>.SIMLIB      ! combined SIMLIB file
#
#      History
#
#  Dec 8 2017: S.Hinton - merge duplicates
#  Jan 8 2018: RK - add VPEC_FILE option
#
#  Nov 16 2019 RK - bug fix writing AB_SED and BD17_SED
#
#  Dec 08 2020 RK - fix to properly handle multiple surveys using
#                    same filters
#
#  Dec 17 2020 RK - write full filter names to text files (e.g., CFA3K/l)
#      for visual convenience. snana.car was modified to strip last
#      char of filter column. FITS format still single char.
#
# Oct 13 2021 RK - fix snana.exe arg bugs in get_survey() function
#
# Nov 11 2021 RK
#  + add optional "CHANGE_FILTER_CHAR: False" key to not add extra char
#  + add optional KEYLIST_REMOVE key to remove obsolete/garbage keys
#  + in combined kcor input, don't duplicate filter paths
#  + parse mulitple filter paths in each input kcor
#     (don't assume just 1 filter set like PS1 or SDSS)
#  * NOT DONE: when merging SN with same name, check spectra
#         (this might work, but didn't check it)
#
# Jan 18 2022 RK - few fixes for Pantheon+
#    + addd new KCOR_PATH arg
#    + all comment lines in [VERSION].LIST
#
# ====================================

import os, sys
import numpy as np
import time, string, getpass, yaml, gzip
import subprocess, shutil, logging, datetime


# globals

SNDATA_ROOT = os.environ['SNDATA_ROOT']
HOSTNAME    = os.environ['HOSTNAME']
CWD         = os.getcwd() 
USERNAME    = getpass.getuser()
tnow        = datetime.datetime.now()
TSTAMP      = f"{tnow.year:04d}-{tnow.month:02d}-{tnow.day:02d}"
MXFILTERS   = 62
TOPDIR_DATA = SNDATA_ROOT + '/lcmerge'
SUBDIR_DUPLICATES = "DUPLICATES"

# define list of new characters for combined data
FILTER_CHARLIST = 'abcdefghijklmnopqrstuvwxyz' + \
                  'ABCDEFGHIJKLMNOPQRSTUVWXYZ' + \
                  '0123456789'

# keys read from input file
KEYNAME_VERSION        = "VERSION_LIST"     # for input data versions
KEYNAME_SURVEY_OUT     = "SURVEY_OUT"       # name of combined data version
KEYNAME_VPEC_FILE      = "VPEC_FILE"        # table file with VPEC and VPEC_ERR
KEYNAME_PRIVATE        = "PRIVATE_DATA_PATH"   # location of data
KEYNAME_KCOR_PATH      = "KCOR_PATH"
KEYNAME_CHANGE_CHAR    = "CHANGE_FILTER_CHAR"  # True or False
KEYNAME_KEYLIST_REMOVE = "KEYLIST_REMOVE"      # e.g., PRIVATE VERSION_PHOTOMETRY
KEYNAME_MAGSYS_FINAL   = "MAGSYSTEM_FINAL"     # e.g., pick AB or BD17 if both exist


# if update_data_files.py is in cwd, use it for debugging;
# else use public version
PROGRAM_UPDATE_DATA_FILES = "update_data_files.py"
if os.path.exists(PROGRAM_UPDATE_DATA_FILES):
    PROGRAM_UPDATE_DATA_FILES = './' + PROGRAM_UPDATE_DATA_FILES

OUTDIR_TEMP = "TEMP_DATA_FILES" # write updated data files with VPEC here

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

#    if ( vbose > 0 ):   print '\t ', key, arg
        
    return(arg)


def change_filterChar(versionInfo, kcorInfo_list, k ):

    kcorInfo           = kcorInfo_list[k]
    input_config       = versionInfo.input_config
    CHANGE_FILTER_CHAR = input_config[KEYNAME_CHANGE_CHAR]

    if CHANGE_FILTER_CHAR:
        print(f"  Change filter characters ")
    else:
        print(f"  Store filter characters ")

    sys.stdout.flush() 

    # - - - - - 
    s = -1 
    for survey, filter_list, repeat_filt, repeat_survey in \
        zip(kcorInfo.SURVEY_LIST,
            kcorInfo.FILTER_LIST, 
            kcorInfo.REPEAT_FILTPATH_LIST,
            kcorInfo.REPEAT_SURVEY_LIST  ):

        s += 1
        repeat = (repeat_filt and repeat_survey)
        new    = True
        # xxx if repeat : continue

        ifilt = 0
        for line in filter_list :
            wdlist = line.split()
            filter_name  = wdlist[0]  # full name of filter

            N = versionInfo.NFILTER_TOT
            filter_oldChar = filter_name[-1]     # current filter char        

            if repeat : 
                kcorInfo.FILTER_CHARLIST_OLD += filter_oldChar # 1-char
                continue

            if new :
                if CHANGE_FILTER_CHAR:
                    filter_newChar = FILTER_CHARLIST[N]  # new char for filter
                else:
                    filter_newChar = filter_oldChar

                versionInfo.FILTER_CHARLIST_OLD += filter_oldChar
                versionInfo.FILTER_CHARLIST_NEW += filter_newChar
                versionInfo.NFILTER_TOT    += 1
            else:
                filter_newChar = kcorInfo.FILTER_CHARLIST_NEW[ifilt]

            # append new char to end of filter name to ensure
            # unique filter-char for each band
            if CHANGE_FILTER_CHAR :
                filter_name_new = filter_name + '/' + filter_newChar
            else:
                filter_name_new = filter_name 

            #print(f" xxx load '{filter_name_new} for s={s} and ifilt={ifilt}")
            kcorInfo.FILTER_NEWNAME.append(filter_name_new)     # full name
            kcorInfo.FILTER_NEWNAME_LIST[s][ifilt] = filter_name_new

            kcorInfo.FILTER_CHARLIST_OLD += filter_oldChar # 1-char
            kcorInfo.FILTER_CHARLIST_NEW += filter_newChar
            ifilt += 1

        # - - - - - - - - 
        # if NEW filter list isn't set, search previous kcorInfo lists
        if len(kcorInfo.FILTER_CHARLIST_NEW) == 0:
            OLD = kcorInfo.FILTER_CHARLIST_OLD  # OLD means original name            
            for kk in range(0,k):
                old_kk = kcorInfo_list[kk].FILTER_CHARLIST_OLD
                if old_kk == OLD :
                    kcorInfo.FILTER_CHARLIST_NEW = \
                        kcorInfo_list[kk].FILTER_CHARLIST_NEW
                    kcorInfo.FILTER_NEWNAME      = \
                        kcorInfo_list[kk].FILTER_NEWNAME

        # print band replacement for each version
        OLD = kcorInfo.FILTER_CHARLIST_OLD
        NEW = kcorInfo.FILTER_CHARLIST_NEW
        tmp = OLD + ' --> ' + NEW
        
        print(f"    {survey:<28.28} bands: {tmp}   (new={new})")        
        sys.stdout.flush() 

    #sys.exit("\n xxx DEBUG EXIT xxx \n")
    return

    # end change_filterChar

def write_kcor_inputFile(versionInfo, kcorInfo_list):

    # write kcor input file for combined data versions.
    fname             = versionInfo.kcor_inFile
    input_config      = versionInfo.input_config
    PRIVATE_DATA_PATH = input_config[KEYNAME_PRIVATE]
    MAGSYS_FINAL      = input_config[KEYNAME_MAGSYS_FINAL]

    print(f"\n Create {fname}")

    if MAGSYS_FINAL is not None:
        print(f"\t Transform all MAGSYSTEMs to {MAGSYS_FINAL}")
    sys.stdout.flush() 

    f = open(fname,"wt")

    f.write(f"# Combined kcor file create by\n#   {sys.argv}\n" )
    f.write(f"# User = {getpass.getuser()} \n")
    f.write(f"# Host = {HOSTNAME}\n"  )
    f.write(f"# Date = {TSTAMP}\n"  )
    f.write(f"\n")
        
    f.write(f"SN_SED:       {kcorInfo_list[0].SN_SED} \n")

    # check all kcor files for primary(s)
    USE_BD17 = False ; USE_AB = False
    n_prim = 0
    for kcorInfo in kcorInfo_list :
        if ( len(kcorInfo.BD17_SED) > 0 and  not USE_BD17 ) :
            f.write(f"BD17_SED:     {kcorInfo.BD17_SED}\n")
            USE_BD17 = True
            n_prim += 1

        if ( len(kcorInfo.AB_SED) > 0  and  not USE_AB) :
            f.write(f"AB_SED:       {kcorInfo.AB_SED} \n")
            USE_AB = True
            n_prim += 1

    # - - - - -
    L0 = kcorInfo_list[0].LAMRANGE[0]
    L1 = kcorInfo_list[0].LAMRANGE[1] 
    f.write(f"LAMBDA_RANGE: {L0} {L1} \n" )

    
    # - - - - -
    # check for spectrograph (Jan 2021)
    # (warning: need to abort of spectcro_file names are different)
    n_spectro = 0
    for kcorInfo in kcorInfo_list :
        spectro_file = kcorInfo.SPECTROGRAPH
        if spectro_file is not None and n_spectro == 0 :
            f.write(f"SPECTROGRAPH: {spectro_file} \n")
            n_spectro += 1

    # - - - - 
    f.write(f"OUTFILE:      {versionInfo.kcor_outFile} \n")
    
    # - - - - - - - - - - - - - - -
    # loop over filter sets
    k = 0
    for kcorInfo in kcorInfo_list :

        V           = versionInfo.NAME[k]
        survey_name = versionInfo.SURVEY_INP[k]
        k += 1

        nfsys = len(kcorInfo.MAGSYSTEM_LIST)
        for s in range(0,nfsys):        
            magsys          = kcorInfo.MAGSYSTEM_LIST[s]
            filtsys         = kcorInfo.FILTSYSTEM_LIST[s]
            filtpath        = kcorInfo.FILTPATH_LIST[s]
            survey          = kcorInfo.SURVEY_LIST[s]  # redundant with survey_name?
            filter_orig     = kcorInfo.FILTER_LIST[s]
            filter_newname  = kcorInfo.FILTER_NEWNAME_LIST[s]
            repeat_filtpath = kcorInfo.REPEAT_FILTPATH_LIST[s]
            repeat_survey   = kcorInfo.REPEAT_SURVEY_LIST[s]
            nfilter         = len(filter_orig)

            if MAGSYS_FINAL is not None:
                if magsys != MAGSYS_FINAL :
                    magsys += f"->{MAGSYS_FINAL}"

            if len(filter_orig) != len(filter_newname) : continue

            #if repeat_survey : continue  #.xyz
            #if repeat_filtpath : continue
            
            if filter_newname[0] == '' : continue
            if "IGNORE" in survey_name : continue 
                
            f.write(f"\n" )
            f.write(f"# Start filters for VERSION = {V}\n" )
            f.write(f"MAGSYSTEM:   {magsys} \n")
            f.write(f"FILTSYSTEM:  {filtsys} \n")
            f.write(f"FILTPATH:    {filtpath} \n")
            f.write(f"SURVEY:      {survey} \n")
        
            nfilt=0
            for filter_line, filter_name_new in \
                zip(filter_orig, filter_newname ):
                filter_name  = filter_line.split()[0]  # old filter name
                filter_file  = filter_line.split()[1]
                filter_zpoff = filter_line.split()[2]
                f.write(f"FILTER:  {filter_name_new:<20} " \
                        f"{filter_file}  {filter_zpoff} \n")

    f.close()
    # end write_kcor_inputFile

def get_survey(PRIVATE_DATA_PATH,VERSION):

    # run snana.exe on version, then read yaml file to get SURVEY name.

    survey = None
    prefix = (f"OUT_{VERSION}")
    cmd = "snana.exe NOFILE "
    cmd += f"VERSION_PHOTOMETRY={VERSION} "
    cmd += f"SNTABLE_LIST=NONE "
    cmd += f"OPT_YAML=1 " 
    # xxx ??? no such key in snana.exe cmd += (f"READ_SPECTRA=F " )
    cmd += f"TEXTFILE_PREFIX={prefix} " 
    if len(PRIVATE_DATA_PATH) > 2 :
        cmd += f"PRIVATE_DATA_PATH={PRIVATE_DATA_PATH} "
    cmd += f" > {prefix}.LOG"
    #print(f" xxx cmd = {cmd}\n")
    os.system(cmd)

    # read from YAML file
    yaml_file = (f"{prefix}.YAML")
    with open(yaml_file,"rt") as y :
        for line in y:
            word_list = line.split()
            if word_list[0] == "SURVEY:" : survey = word_list[1]

    # remove junk files, and be careful not to   rm .* by accident
    if len(prefix) > 2 :
        cmd_rm = (f"rm {prefix}.*")
        os.system(cmd_rm)

    #sys.exit(f"\n xxx {cmd_rm}\n")

    return survey 

    # end get_survey

def run_kcor(versionInfo):
    inFile  = versionInfo.kcor_inFile
    logFile = versionInfo.kcor_logFile
    cmd = 'kcor.exe ' + inFile + ' > ' + logFile
    print(f" Run kcor program ... ")
    sys.stdout.flush() 
    os.system(cmd)

    # check for fatal error in log file
    f      = open(logFile,"rt")
    Lines  = f.readlines()
    f.close
    if any(" ABORT " in s for s in Lines):
        msg = "\nFATAL error running kcor program:\n   Check %s\n" % inFile
        sys.exit(msg)
    
class VERSION_INFO:
    def __init__(self,filename):

        line_list = []
        with  open(filename,"rt") as f:
            for line in f:
                if len(line) == 0 : continue
                line_list.append(line)

        input_config = yaml.safe_load("\n".join(line_list))

        # - - - -
        # load a few defaults if not specified in input config file
        if KEYNAME_PRIVATE not in input_config:
            input_config[KEYNAME_PRIVATE] = ""

        if KEYNAME_KCOR_PATH not in input_config:
            input_config[KEYNAME_KCOR_PATH] = None

        if KEYNAME_CHANGE_CHAR not in input_config:
            input_config[KEYNAME_CHANGE_CHAR] = True

        if KEYNAME_VPEC_FILE not in input_config:
            input_config[KEYNAME_VPEC_FILE] = ""
        
        if KEYNAME_KEYLIST_REMOVE not in input_config:
            input_config[KEYNAME_KEYLIST_REMOVE] = ""

        if KEYNAME_MAGSYS_FINAL not in input_config:
            input_config[KEYNAME_MAGSYS_FINAL] = None

        # - - - - - - 
        SOUT   = input_config[KEYNAME_SURVEY_OUT]
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
        
        # define strings of old and new filter char;
        # to be filled later.
        self.FILTER_CHARLIST_OLD = ""
        self.FILTER_CHARLIST_NEW = ""
        self.NFILTER_TOT    = 0
        self.NDATA_FILE     = []   # to be filled later
        
        self.NAME = []
        self.INFILE_KCOR = []
        KCOR_PATH  = input_config[KEYNAME_KCOR_PATH]
        for tmpLine in input_config[KEYNAME_VERSION]:
            name        = tmpLine.split()[0]
            infile_kcor = tmpLine.split()[1]
            if KCOR_PATH is not None:
                infile_kcor = f"{KCOR_PATH}/{infile_kcor}"
            self.NAME.append(name)
            self.INFILE_KCOR.append(infile_kcor)
                
        # get survey name with snana job
        print(f"")
        self.SURVEY_INP = []
        PRIVATE_DATA_PATH = input_config[KEYNAME_PRIVATE]
        for V in self.NAME :
            survey_name = get_survey(PRIVATE_DATA_PATH,V)
            print(f"    VERSION {V:<28.28} -> {survey_name} ")
            sys.stdout.flush() 
            self.SURVEY_INP.append(survey_name)

        self.input_config = input_config
        print(f"")
        sys.stdout.flush() 

        return

        
class KCOR_INFO:
    def __init__(self, filename, FILTPATH_LIST_GLOBAL):

        print(f"   Parse kcor input file: {filename}")
        sys.stdout.flush() 
        filename_expandvars = os.path.expandvars(filename)
        fname_local         = filename_expandvars

        # open file file and read all lines into Lines
        with open(fname_local,"rt") as f:
            line_list = f.read().splitlines() 
            
        NFSYS = 0
        self.SPECTROGRAPH      = None
        self.BD17_SED          = ""
        self.AB_SED            = ""

        self.MAGSYSTEM_LIST    = [] 
        self.FILTSYSTEM_LIST   = [] 
        self.FILTPATH_LIST     = [] 
        self.SURVEY_LIST       = [] 
        self.FILTER_LIST          = []   # list of lists
        self.FILTER_NEWNAME_LIST  = []
        self.NFILTER_LIST      = []

        self.REPEAT_FILTPATH_LIST = [] # keep track of repeated filtpath
        self.REPEAT_SURVEY_LIST   = [] # keep track of repeated SURVEY

        for line in line_list:
            wdlist = line.split()
            if len(wdlist) == 0   : continue
            if wdlist[0][0] == '#' : continue
            key = wdlist[0]            

            if key == 'SN_SED:' :
                self.SN_SED = wdlist[1]
            elif key == 'BD17_SED:' :
                self.BD17_SED = wdlist[1]
            elif key == 'AB_SED:' :
                self.AB_SED = wdlist[1]
            elif key == 'LAMBDA_RANGE:' :
                self.LAMRANGE = wdlist[1:3]
            elif key == 'SPECTROGRAPH:' :
                self.SPECTROGRAPH = wdlist[1]
                if len(self.SPECTROGRAPH) == 0 :  self.SPECTROGRAPH = None

            elif key == 'MAGSYSTEM:' :
                self.MAGSYSTEM_LIST.append(wdlist[1])
                NFSYS += 1
                self.FILTER_LIST.append([])
                self.FILTER_NEWNAME_LIST.append([])
                self.NFILTER_LIST.append(0)
            elif key == 'FILTSYSTEM:' :
                self.FILTSYSTEM_LIST.append(wdlist[1])
            elif key == 'SURVEY:' :
                survey = wdlist[1]
                self.SURVEY_LIST.append(survey)

                REPEAT = survey in GLOBAL_DICT['SURVEY_LIST']
                if not REPEAT: GLOBAL_DICT['SURVEY_LIST'].append(survey)
                self.REPEAT_SURVEY_LIST.append(REPEAT)

            elif key == 'FILTPATH:' :
                filtpath = wdlist[1]
                self.FILTPATH_LIST.append(filtpath)

                REPEAT = filtpath in GLOBAL_DICT['FILTPATH_LIST']
                if not REPEAT: GLOBAL_DICT['FILTPATH_LIST'].append(filtpath)
                self.REPEAT_FILTPATH_LIST.append(REPEAT)
            elif key == 'FILTER:' :
                str_tmp = '  '.join(wdlist[1:])
                self.FILTER_LIST[NFSYS-1].append(str_tmp)
                self.FILTER_NEWNAME_LIST[NFSYS-1].append("")

                NF = len(self.FILTER_LIST[NFSYS-1])
                self.NFILTER_LIST[NFSYS-1] = NF

        # - - - - - - - - - -
        self.N_FILTERSYSTEM = len(self.FILTSYSTEM_LIST)
        NFSYS = self.N_FILTERSYSTEM

        # - - - - - - -
        dmp_flag =  False # True 
        if dmp_flag :
            print(f" xxx -------------------------------------- ")
            print(f" xxx MAGSYSTEM_LIST  = {self.MAGSYSTEM_LIST}")
            print(f" xxx FILTSYSTEM_LIST = {self.FILTSYSTEM_LIST}")
            print(f" xxx FILTPATH_LIST   = {self.FILTPATH_LIST}")
            for f_list, fnew_list in zip(self.FILTER_LIST, self.FILTER_NEWNAME_LIST):
                print(f" xxx FILTER_LIST     = {f_list}")
                print(f" xxx FILTER_NEWNAME_LIST = {fnew_list}")
            print(f" xxx NFILTER_LIST         = {self.NFILTER_LIST} ")
            print(f" xxx REPEAT_FILTPATH_LIST = {self.REPEAT_FILTPATH_LIST}")
            print(f" xxx REPEAT_SURVEY_LIST   = {self.REPEAT_SURVEY_LIST}")
            sys.stdout.flush() 
        # define things to change later
        self.FILTER_NEWNAME        = []  # one big list
        # xxx self.new            = True     # assume new filter/mag system
        self.FILTER_CHARLIST_OLD   = ""
        self.FILTER_CHARLIST_NEW   = ""

        return

        # end KCOR_INFO

def create_newVersion(versionInfo):

    VOUT = versionInfo.VERSION_OUT_TEXT  # name of output version

    # create new version-subDir
    if (len(VOUT) == 0 ):
        sys.exit("Output directory not defined\n Check SURVEY_OUT key")

    if ( os.path.exists(VOUT) ):
        print(f" Remove pre-existing {VOUT}")
        sys.stdout.flush() 
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

    return
    # end create_newVersion
    
def add_newVersion(VIN,versoinInfo,kcorInfo):

    # add new data version to combined version;
    # copy all data files and use 'sed' utility
    # to make changes to filters and SURVEY name.
    
    dataDir   = os.path.expandvars(TOPDIR_DATA) + '/' + VIN
    SOUT      = versionInfo.SURVEY_OUT
    VOUT_TEXT = versionInfo.VERSION_OUT_TEXT
    
    # read list file
    LISTFILE_IN  = dataDir + '/' + VIN + '.LIST'
    LISTFILE_OUT = VOUT_TEXT + '/' + versionInfo.AUXFILE_LIST    
    PTR_L        = open(LISTFILE_IN,"rt")
    fileList     = PTR_L.readlines()
    PTR_L.close

    # read contents of first file
    fileList[0] = fileList[0].strip()
    first_fileName = dataDir + '/' + fileList[0]
    first_fileName = first_fileName.replace("\n", "")
    first_fileName_gz = first_fileName + ".gz"

    try:
        f0    = open(first_fileName,"rt")
        gz_flag = False
    except:
        f0    = gzip.open(first_fileName_gz,"rt")
        gz_flag = True

    fileContents = np.array(f0.readlines())
    f0.close
    
    # read FILTER string from first data file
    FILTERSTRING_OLD = parseLines(fileContents, 'FILTERS:', 1, 0)
    SURVEY           = parseLines(fileContents, 'SURVEY:',  1, 0)
    VARLIST          = parseLines(fileContents, 'VARLIST:', 1, 0) # ??

    # get full filter lists from kcor file
    FILTERLIST_OLD = kcorInfo.FILTER_CHARLIST_OLD     # only this version
    FILTERLIST_NEW = kcorInfo.FILTER_CHARLIST_NEW     # only this version
    FILTERLIST_ALL = versionInfo.FILTER_CHARLIST_NEW  # all filters
    FILTER_NEWNAME = kcorInfo.FILTER_NEWNAME

    print(f"\t {FILTERLIST_OLD} -> {FILTERLIST_NEW} " )
    sys.stdout.flush() 

    # open new list file in append mode
    PTR_NEWLIST = open(LISTFILE_OUT,"at")

    # - - - - - start construction of 'sed' command - - - - - - -
    sedcmd = "sed "
    
    # replace SURVEY ... only first occurance !
    sedAdd  = f"-e '0,/SURVEY:/s/{SURVEY}/ {SOUT}({SURVEY})/' "
    sedcmd += sedAdd
    
    # Replace global filter string.
    OLD  = FILTERSTRING_OLD  # from data file
    OLD2 = FILTERLIST_OLD    # from kcor file
    NEW  = FILTERLIST_NEW    # new kcor list
    ALL  = FILTERLIST_ALL    # all filters from all files
    # xxx mark sedAdd = "-e 's/%s/%s  # %s -> %s/g' " % (OLD,ALL,OLD2,NEW) 
    sedAdd  = f"-e 's/{OLD}/{ALL}  # {OLD2} -> {NEW}/g' " 
    sedcmd += sedAdd

    # Replace each single-char band with new full filter name
    # Add '??' before each band to avoid removing new band
    # that matches an old band. Then remove ?? separately.

    for old, new, newname in \
        zip(FILTERLIST_OLD, FILTERLIST_NEW, FILTER_NEWNAME):
        filter_newname = newname.replace("/","\/") # new full name
        sedAdd = f"-e 's/ {old} /  ??{filter_newname}  /g' "
        sedcmd += sedAdd
        
    #sys.exit("\n xxx DEBUG STOP xxx \n")

    # remove the temporary '?s?'
    sedcmd += "-e 's/??//g' "
    
    # check for keys to remove
    input_config = versionInfo.input_config
    KEYLIST_REMOVE = input_config[KEYNAME_KEYLIST_REMOVE].split()
    for key in KEYLIST_REMOVE :
        sedcmd += f"-e '/{key}/d' "

    # loop over all files and run 'sedcmd' 
    nfile = 0
    for infile in fileList:
        if infile[0] == '#' : continue 
        infile      = infile.strip()
        infile_base = infile.replace("\n", "")
        infile_copy = infile_base 
        if gz_flag: 
            infile_copy += '.gz'
        FIN    = dataDir   + '/' + infile_copy
        FOUT   = VOUT_TEXT + '/' + infile_base

        if gz_flag :
            # run sed on gzipped files, output regular files.
            # The output here is gzipped later after merge_duplicate
            SEDCMD = f"zcat {FIN} | {sedcmd} > {FOUT}"
        else:
            SEDCMD = f"{sedcmd} {FIN} > {FOUT}"
                    
        os.system(SEDCMD)
        PTR_NEWLIST.write(f"{infile_base}\n")
        nfile += 1

    versionInfo.NDATA_FILE.append(nfile)
    PTR_NEWLIST.close()

    return
    # end add_newVersion

def write_readme(versionInfo, kcorInfo_list):

    # write DOCUMENTATION block in README

    input_config     = versionInfo.input_config
    VOUT             = versionInfo.VERSION_OUT_TEXT  # name of output version
    AUXFILE_README   = versionInfo.AUXFILE_README
    CHANGE_FILTER_CHAR = input_config[KEYNAME_CHANGE_CHAR]
    VPEC_FILE          = input_config[KEYNAME_VPEC_FILE] 

    README_OUTFILE   = VOUT + '/' + AUXFILE_README
    f = open(README_OUTFILE,"wt")

    f.write(f"DOCUMENTATION:\n")
    f.write(f"  PURPOSE:  combined data set for analysis\n")    
    f.write(f"  USAGE_KEY: VERSION_PHOTOMETRY\n")    
    f.write(f"  USAGE_CODE: snlc_fit.exe\n")    
    f.write(f"  NOTES:\n")    

    cmd = ' '.join(sys.argv)
    f.write(f"  - created with command {cmd}\n") 
    if len(VPEC_FILE) > 1:
        f.write(f"  - VPEC included from {VPEC_FILE}\n")

    k = 0
    for vname, ndata_file in   zip(versionInfo.NAME, versionInfo.NDATA_FILE ):
        kcorInfo  = kcorInfo_list[k]
        OLD       = kcorInfo.FILTER_CHARLIST_OLD 
        NEW       = kcorInfo.FILTER_CHARLIST_NEW 
        k += 1

        txt1  = f"{ndata_file:3d} data files from {vname:>28}"        
        if CHANGE_FILTER_CHAR:
            txt2  = f"{OLD} -> {NEW}"
        else:
            txt2 = ""
        f.write(f"  - {txt1} {txt2}\n")

    nfile_orig  = versionInfo.NFILE_ORIG
    nfile_final = versionInfo.NFILE_FINAL
    nfile_dupl  = versionInfo.NFILE_DUPL
    f.write(f"  - {nfile_orig} original data files were processed.\n")
    f.write(f"  - {nfile_final} combined data files includes " \
            f"{nfile_dupl} merged duplicates.\n")

    f.write(f"\n")    
    f.write(f"  VERSION:\n")    
    f.write(f"  - DATE:      {TSTAMP}\n")
    f.write(f"    USERNAME:  {USERNAME} \n")
    f.write(f"    HOSTNAME:  {HOSTNAME}\n")
    f.write(f"    DIRECTORY: {CWD}\n")
    f.write(f"DOCUMENTATION_END: \n\n")
    f.close()

    return
    # end write_readme

def merge_duplicates(versionInfo):

    # For each pair of duplicates:
    # Get list of duplicates : filename & SNID
    # Merge "OBS:" lines into one file, keeping header of first file.
    # Update NOBS key-value and make sure that "END:" is at the end
    # Move duplicates into /DUPLICATES subDir
    #
    # Jan 21 2022 RK - skip END_PHOTOMETRY keys (used in PPLUS)

    VOUT_TEXT  = versionInfo.VERSION_OUT_TEXT 

    # get list of SNID using grep :  file  SNID
    cmd    = (f"cd {VOUT_TEXT} ; grep SNID: *.* ")
    output = subprocess.check_output(cmd, shell=True)
    logging.basicConfig(level=logging.DEBUG, 
                        format="[%(funcName)20s()] %(message)s")

    output = output.decode('utf-8')
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
    ndupl      = len(duplicates)
    print(f"\n Found {ndupl} duplicates:" )
    for sn in duplicates :
        print(f"    {sn} {d[sn]} ") 

    
    dup_dir = VOUT_TEXT + f"/{SUBDIR_DUPLICATES}"
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
                if line.startswith("END:"):             continue
                if line.startswith("END_PHOTOMETRY:"):  continue
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

    # compress  duplicate dir (Jan 25 2022)
    cmd_tar = f"cd {VOUT_TEXT} ; " \
              f"tar -czf {SUBDIR_DUPLICATES}.tar.gz {SUBDIR_DUPLICATES}; " \
              f"rm -r {SUBDIR_DUPLICATES}"
    os.system(cmd_tar)

    # Update list
    list_file = f"{VOUT_TEXT}/{VOUT_TEXT}.LIST"
    with open(list_file) as f:
        list_sn = [l for l in f]

    to_remove = [os.path.basename(filename) \
                 for sn in duplicates for filename in d[sn]]
    final_list = to_add + [f for f in list_sn \
                           if "".join(f.split()) not in to_remove]

    nfile_orig  = len(list_sn)
    nfile_final = len(final_list)

    with open(list_file, "w") as f:
        f.writelines(final_list)

    logging.info(f"Updated list file {list_file}: " \
                 f" {nfile_orig} -> {nfile_final} data files." ) 
    print('')
    versionInfo.NFILE_FINAL = nfile_final
    versionInfo.NFILE_ORIG  = nfile_orig
    versionInfo.NFILE_DUPL  = ndupl

    return 

    # end merge_duplicates

def gzip_newVersion(versionInfo):
    VOUT_TEXT  = versionInfo.VERSION_OUT_TEXT 
    
    file_list = glob.glob1(VOUT_TEXT,"*.gz")
    n_file    = len(file_list)

    print(f" gzip {n_file} files in {VOUT_TEXT}")
    sys.stdout.flush() 

    cmd = f"cd {VOUT_TEXT}; gzip *dat *.DAT 2>/dev/null"
    os.system(cmd)
    return
    # end gzip_newVersion

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

    VOUT_TEXT    = versionInfo.VERSION_OUT_TEXT
    input_config = versionInfo.input_config
    VPEC_FILE    = input_config[KEYNAME_VPEC_FILE]

    if ( len(VPEC_FILE) < 2 ): return

    print(f"\n Add VPEC and VPEC_ERR from {VPEC_FILE}")
    sys.stdout.flush() 
    LOGFILE    = "ADD_VPEC.LOG"

    cmd  = f"{PROGRAM_UPDATE_DATA_FILES} " \
           f"-V ./{VOUT_TEXT} " \
           f"-u {VPEC_FILE}  " \
           f"-o {OUTDIR_TEMP} " \
           f"-v VPEC,VPEC_ERR " \
           f" > {LOGFILE} " 
    #sys.exit(f"\n xxx cmd = \n{cmd}")
    os.system(cmd)

    # - - - - - - 
    # move updated files out of OUTDIR_TEMP 
    cmd = f"mv {OUTDIR_TEMP}/{VOUT_TEXT}/*.* {VOUT_TEXT}/ ; " \
          f"rm -r {OUTDIR_TEMP}"
    os.system(cmd)

    return

    # end add_vpec


def make_simlib(versionInfo):

    print(f"\n Create SIMLIB ")
    sys.stdout.flush() 
    SOUT       = versionInfo.SURVEY_OUT
    VOUT_TEXT  = versionInfo.VERSION_OUT_TEXT
#    ZP_TEXT = versionInfo.ZPERR 
    nmlFile    = 'make_simlib_' + SOUT + '.nml'
    logFile    = 'make_simlib_' + SOUT + '.log'
    simlibFile = versionInfo.simlibFile

    # open with 1 line buffer so that it is ready to run snana.exe
    PTR_NML = open(nmlFile,"wt",1)

    PTR_NML.write(f" &SNLCINP\n")
    PTR_NML.write(f"   VERSION_PHOTOMETRY = '{VOUT_TEXT}' \n" )
    PTR_NML.write(f"   PRIVATE_DATA_PATH  = './' \n")
    PTR_NML.write(f"   KCOR_FILE    = '{versionInfo.kcor_outFile}' \n")

#    PTR_NML.write("   SIMLIB_ZPERR_LIST    = '%s' \n" % ZP_TEXT)
        
    
    PTR_NML.write(f"   SNTABLE_LIST = '' \n")
    PTR_NML.write(f"   SIMLIB_OUT   = '{simlibFile}' \n" )
    
    PTR_NML.write(f" \n")
    PTR_NML.write(f" &END\n\n")
    PTR_NML.flush
    PTR_NML.close

    cmd  = "snana.exe %s > %s" % (nmlFile,logFile)
    os.system(cmd)

    
def convert2FITS(versionInfo):

    # convert TEXT format into FITS format
    S      = versionInfo.SURVEY_OUT
    V_TEXT = versionInfo.VERSION_OUT_TEXT
    V_FITS = versionInfo.VERSION_OUT_FITS

    print(f"\n Convert TEXT formatted data files into FITS: ")
    print(f"\t {V_TEXT} -> {V_FITS}\n")
    sys.stdout.flush() 

    # remove output FITS dir if it exits
    if ( os.path.exists(V_FITS) ):    shutil.rmtree(V_FITS)

    # construct snana nmlFile
    prefix     = "convert2FITS"
    nml_file   = f"{prefix}_{S}.nml"
    log_file   = f"{prefix}_{S}.log"

    PTR_NML   = open(nml_file,"wt",1)
    PTR_NML.write(f" &SNLCINP \n")
    PTR_NML.write(f"   VERSION_PHOTOMETRY = '{V_TEXT}' \n")
    PTR_NML.write(f"   PRIVATE_DATA_PATH  = '{CWD}' \n"  )
    PTR_NML.write(f"   VERSION_REFORMAT_FITS = '{V_FITS}' \n")
    PTR_NML.write(f" &END \n\n")
    PTR_NML.flush()
    PTR_NML.close()

    cmd = f"snana.exe {nml_file} > {log_file} ; mv {prefix}* {V_FITS}"

    os.system(cmd)

    # end convert2FITS
    
def printSummary(versionInfo):
    print(f"\n# =============================================== ")
    print(f"  Summary of outputs: ")

    V = versionInfo.VERSION_OUT_TEXT
    print(f"\t {V}/    (combined Data version)" )

    V = versionInfo.VERSION_OUT_FITS
    print(f"\t {V}/    (combined Data version)" )

    out_file = versionInfo.kcor_outFile
    print(f"\t {out_file}    (combined kcor file)" )

    simlib_file = versionInfo.simlibFile
    print(f"\t {simlib_file}    (combined SIMLIB file)" )
    sys.stdout.flush() 

    # end printSummary

# =========================
# ======= MAIN ============
# =========================

if __name__ == "__main__":

# parse input argument(s)
    if ( len(sys.argv) < 2 ):
        sys.exit("Must give INFILE arguent\n-->ABORT")
    else:
        INFILE = sys.argv[1]
        print(f"Input file: {INFILE}")

    print(f" SNDATA_ROOT = {SNDATA_ROOT}")
    sys.stdout.flush() 

    versionInfo = VERSION_INFO(INFILE)
    input_config = versionInfo.input_config

    # ----------------------------------------
    # prepare kcor file with all filters
    # ----------------------------------------
    kcorInfo_list        = []
    GLOBAL_DICT = { 'SURVEY_LIST': [],    'FILTPATH_LIST' : []  }

    for kcorFile in versionInfo.INFILE_KCOR:
        kcorInfo = KCOR_INFO(kcorFile, GLOBAL_DICT)
        kcorInfo_list.append(kcorInfo)
        
    nkcor = len(kcorInfo_list)
    print(f"\n Done parsing {nkcor} kcor-input files ")
    sys.stdout.flush() 

    # combine multuple surveys with same FILTPATH (Dec 8 2020)
    # xxx mark delete ? combine_duplicate_filtpath(versionInfo, kcorInfo)

    # change filter char
    print(' ')
    for k in range(0,nkcor):
        change_filterChar(versionInfo, kcorInfo_list, k)

    # write new combined kcor-input file
    write_kcor_inputFile(versionInfo, kcorInfo_list) 
    run_kcor(versionInfo)
    
    # ---------------------------------------------------------
    # now re-write data files with different filter strings
    # ---------------------------------------------------------

    print(f"\n# - - - - - - - - - - - - - - - - - - - - - - ")
    
    # check for private data path
    PRIVATE_DATA_PATH = input_config[KEYNAME_PRIVATE]
    if ( len(PRIVATE_DATA_PATH) > 0 ):
        TOPDIR_DATA = PRIVATE_DATA_PATH
        print(f" Use PRIVATE_DATA_PATH = {TOPDIR_DATA}")

    create_newVersion(versionInfo)  
    for vname, kcorInfo in zip(versionInfo.NAME, kcorInfo_list):        
        print(f" Swap filter strings in DATA-VERSION:  {vname}")
        add_newVersion(vname, versionInfo, kcorInfo)

    sys.stdout.flush() 

    # -----------------
    # merge duplicate light curves from different instruments (Dec 2017)
    merge_duplicates(versionInfo)

    # gzip files in new version (after merging duplicates
    gzip_newVersion(versionInfo)

    # write readme file using DOCANA style
    write_readme(versionInfo, kcorInfo_list)

    # check option to add VPEC & VPEC_ERR from a supplemental table file
    add_vpec(versionInfo)
    
    # - - - - - -
    convert2FITS(versionInfo)
    
    # ------------------------------------------
    make_simlib(versionInfo)

    # --------------------
    # print summary of outputs
    printSummary(versionInfo)
    
    print('\n Done.\n')

# ========= END MAIN ================

