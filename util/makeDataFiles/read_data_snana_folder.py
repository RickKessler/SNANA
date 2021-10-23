# Read data from already existing folder in SNANA-FITS format.
# Intended only for testing makeDataFiles; not for production.

import os,sys,glob,yaml,shutil
import logging # , coloredlogs

import makeDataFiles_util  as    util

from   makeDataFiles_params  import *
from   makeDataFiles_base    import Program
from   astropy.io import fits


# - - - - - - - - - - - - - - - - - - -     -
class data_snana_folder(Program):
    def __init__(self, config_inputs) :
        config_data = {}  
        print(" Init data_snana_folder class.")
        super().__init__(config_inputs,config_data)

    def init_read_data(self):
        args = self.config_inputs['args']  # command line args

        data_folder  = os.path.expandvars(args.snana_folder)
        version      = os.path.basename(data_folder)
        list_file    = f"{data_folder}/{version}.LIST"
        
        logging.info(f" Read data version = {version}")
        logging.info(f" from data_foler = {data_folder}")

        # scoop up contents of LIST file
        with open(list_file, 'r') as f:
            HEAD_file_list = f.read().replace('\n', ' ').split()

        self.config_data['HEAD_file_list'] = HEAD_file_list
        self.config_data['n_HEAD_file']    = len(HEAD_file_list)
        self.config_data['version']        = version
        self.config_data['data_folder']    = data_folder
        
        # end init_read_data

    def prep_read_data_subgroup(self, i_subgroup):

        data_folder      = self.config_data['data_folder']
        HEAD_file_base   = self.config_data['HEAD_file_list'][i_subgroup]
        n_HEAD_file      = self.config_data['n_HEAD_file']

        if i_subgroup == n_HEAD_file - 1 :
            return 0 # done reading
        
        HEAD_file       = f"{data_folder}/{HEAD_file_base}"
        nevt, hdu_head  = self.open_snana_fits(HEAD_file)
        
        PHOT_file_base  = hdu_head[0].header['PHOTFILE']
        PHOT_file       = f"{data_folder}/{PHOT_file_base}"
        NROW, hdu_phot  = self.open_snana_fits(PHOT_file)

        logging.info(f"\t Read {nevt} events from {HEAD_file_base}")
        table_head = hdu_head[1].data
        table_phot = hdu_phot[1].data
            
        table_dict = {
            'head_file'  : HEAD_file_base,
            'table_head' : table_head,
            'table_phot' : table_phot
        }

        self.config_data['nevt_subgroup'] = nevt
        self.config_data['table_dict'] = table_dict
        self.config_data['hdu_head']   = hdu_head
        self.config_data['hdu_phot']   = hdu_phot

        return nevt
    
        #end prep_read_data_subgroup

    def end_read_data_subgroup(self):
        self.config_data['hdu_head'].close()
        self.config_data['hdu_phot'].close()

    def end_read_data(self):
        # global end for reading data
        pass
        
    def open_snana_fits(self,file_name):
       # check file_name and file_name.gz, and open the file that exists.
       # Function returns hdu pointer and number of rows in table.

        msgerr = []
        file_namegz = f"{file_name}.gz"
        if os.path.exists(file_namegz) :
            hdul = fits.open(file_namegz)
        elif os.path.exists(file_name):
            hdul = fits.open(file_name)
        else:
            msgerr.append(f"Cannot find fits file")
            msgerr.append(f" {file_name}   not")
            msgerr.append(f" {file_namegz} ")
            util.log_assert(False,msgerr)
            
        NROW = hdul[1].header['NAXIS2']
        return NROW, hdul
    
    # end open_snana_fits

    def read_event(self, evt ):

        msgerr     = []
        table_dict = self.config_data['table_dict']
        args       = self.config_inputs['args']  # command line args

        # read and store one event for row "evt" and return data_dict.
        
        varlist_obs = self.config_data['varlist_obs']

        # define local pointers to head and phot tables from FITS file
        table_head = table_dict['table_head']
        table_phot = table_dict['table_phot']
        
        # init output dictionaries
        head_raw, head_calc = self.reset_data_event_dict()

        try:
            SNID = table_head.SNID[evt].decode('utf-8').replace(' ','')
        except:
            SNID = table_head.SNID[evt]
        head_raw[DATAKEY_SNID]  = SNID
        
        head_raw[DATAKEY_RA]    = table_head.RA[evt]

        # check 'DEC' and legacy column name 'DECL'
        head_raw[DATAKEY_DEC] = \
            self.get_table_value(['DEC','DECL'],evt,table_head)

        head_raw[DATAKEY_zHEL]       = table_head.REDSHIFT_HELIO[evt]
        head_raw[DATAKEY_zHEL_ERR]   = table_head.REDSHIFT_HELIO_ERR[evt]

        # strip off calculated values
        head_calc[DATAKEY_zCMB]          = table_head.REDSHIFT_FINAL[evt]
        head_calc[DATAKEY_zCMB_ERR]      = table_head.REDSHIFT_FINAL_ERR[evt]
        head_calc[DATAKEY_MWEBV]         = table_head.MWEBV[evt]
        head_calc[DATAKEY_MWEBV_ERR]     = table_head.MWEBV_ERR[evt]

        # lightcurve-MJD info. Note that MJD_DETECT_FIRST is optional
        head_calc[DATAKEY_PEAKMJD]       = int(table_head.PEAKMJD[evt])

        if DATAKEY_MJD_DETECT in vars(table_head):  # could be very slow??
            head_calc[DATAKEY_MJD_DETECT] = table_head.MJD_DETECT_FIRST
        else:
            if args.mjd_detect_range is not None:
                msgerr.append(f"Cannot implement args.mjd_detect_range = " \
                              f"{args.mjd_detect_range}")
                msgerr.append(f"Because {DATAKEY_MJD_DETECT} is not in " \
                              f"data header")
                util.log_assert(False,msgerr)


        head_raw[HOSTKEY_OBJID]         = table_head.HOSTGAL_OBJID[evt]
        head_raw[HOSTKEY_SPECZ]         = table_head.HOSTGAL_SPECZ[evt]
        head_raw[HOSTKEY_SPECZ_ERR]     = table_head.HOSTGAL_SPECZ_ERR[evt] 
        head_raw[HOSTKEY_SNSEP]         = table_head.HOSTGAL_SNSEP[evt]

        head_calc[HOSTKEY_PHOTOZ]       = table_head.HOSTGAL_PHOTOZ[evt]
        head_calc[HOSTKEY_PHOTOZ_ERR]   = table_head.HOSTGAL_PHOTOZ_ERR[evt]
        
        # ugly/embarassing hack to get field (DDF or WFD) from filename
        # because original plasticc data didn't store field.
        head_file  = table_dict['head_file']
        if FIELD_DDF in head_file:
            head_raw[DATAKEY_FIELD] = FIELD_DDF
        elif FIELD_WFD in head_file:
            head_raw[DATAKEY_FIELD] = FIELD_WFD
        else:
            msgerr.append(f"Unable to determine FIELD for")
            msgerr.append(f"{head_file}")
            util.log_assert(False,msgerr)
            
        # - - - - - - -
        # get pointers to PHOT table
        ROWMIN = table_head.PTROBS_MIN[evt]
        ROWMAX = table_head.PTROBS_MAX[evt]        
        
        # note that there are NOBS+1 rows, but last row is pad word -777
        NOBS = ROWMAX - ROWMIN
        phot_raw = self.init_phot_dict(NOBS)

        table_column_names = table_phot.columns.names

        for varname in varlist_obs:
            varname_table = varname
            if varname == 'BAND' : varname_table = 'FLT'
            if varname_table not in table_column_names : continue
            phot_raw[varname] = \
                table_phot[varname_table][ROWMIN:ROWMAX].copy()
            
        # - - - -
        spec_raw = {}

        # - - - - -
        # load output dictionary
        data_dict = {
            'head_raw'  : head_raw,
            'head_calc' : head_calc,
            'phot_raw'  : phot_raw,
            'spec_raw'  : spec_raw
        }
        
        return data_dict
    
        # end read_event

    def get_table_value(self, varlist, irow, table):

        # return "irow" table value for varlist,
        # where varlist = ['NAME1', 'NAME2', etccc]
        # will sequentially check NAME1, NAME2, etc ...

        value = None
        for varname in varlist:
            try:
                value = table[varname][irow]
                return value
            except:
                pass  # just try next varname
            
        return value
        # end get_table_value
        

    def set_dump_flag(self, isn, data_event_dict):
        d_raw = data_event_dict['head_raw']
        zhel  = d_raw['REDSHIFT_HELIO']
        dump_flag = False # isn<200 and zhel < 0.25
        return dump_flag
        # set_dump_flag

   
