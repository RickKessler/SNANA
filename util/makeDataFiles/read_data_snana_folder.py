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
        sys.stdout.flush()

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
        
        if i_subgroup == n_HEAD_file -1 :
            return 0 # done reading

        HEAD_file       = f"{data_folder}/{HEAD_file_base}"
        nevt, hdu_head  = self.open_snana_fits(HEAD_file)

        PHOT_file_base  = hdu_head[0].header['PHOTFILE']
        PHOT_file       = f"{data_folder}/{PHOT_file_base}"
        NROW, hdu_phot  = self.open_snana_fits(PHOT_file)

        logging.info(f"   Read {nevt} events from {HEAD_file_base}")
        sys.stdout.flush()
        
        table_head = hdu_head[1].data
        table_phot = hdu_phot[1].data

        head_names = table_head.columns.names
        phot_names = table_phot.columns.names

        # on first subgroup, check for true mag in PHOT table
        if i_subgroup == 0 and  VARNAME_TRUEMAG in phot_names:
            self.append_truemag_obs()
        
        table_dict = {
            'head_file'  : HEAD_file_base,
            'table_head' : table_head,
            'table_phot' : table_phot,
            'head_names' : head_names,
            'phot_names' : phot_names
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
        head_names = table_dict['head_names']
        phot_names = table_dict['phot_names']

        # init output dictionaries
        head_raw, head_calc, head_sim = self.reset_data_event_dict()

        try:
            SNID = table_head.SNID[evt].decode('utf-8').replace(' ','')
        except:
            SNID = table_head.SNID[evt]
        head_raw[DATAKEY_SNID]  = SNID

        head_raw[DATAKEY_RA]    = table_head.RA[evt]

        # check 'DEC' and legacy column name 'DECL'
        head_raw[DATAKEY_DEC] = \
            self.get_table_value(['DEC','DECL'],evt,table_head)

        # lightcurve-MJD info. Note that MJD_DETECT_FIRST is optional
        head_calc[DATAKEY_PEAKMJD]   = int(table_head.PEAKMJD[evt])

        if DATAKEY_MJD_DETECT_FIRST in head_names:
            head_calc[DATAKEY_MJD_DETECT_FIRST] = \
                table_head.MJD_DETECT_FIRST[evt]
            head_calc[DATAKEY_MJD_DETECT_LAST] = \
                table_head.MJD_DETECT_LAST[evt]
        else:
            if args.nite_detect_range is not None:
                msgerr.append(f"Cannot implement args.nite_detect_range = " \
                              f"{args.nite_detect_range}")
                msgerr.append(f"Because {DATAKEY_MJD_DETECT_FIRST} is not in "\
                              f"data header")
                util.log_assert(False,msgerr)

        # - - - - - - -
        # check user sub-sample selection here to avoid reading
        # remainder of header and photometry for rejected events.
        apply_select = True
        if apply_select :
            var_dict = {
                DATAKEY_SNID       : int(SNID),
                DATAKEY_RA         : head_raw[DATAKEY_RA],
                DATAKEY_DEC        : head_raw[DATAKEY_DEC],
                DATAKEY_PEAKMJD    : head_calc[DATAKEY_PEAKMJD],
                DATAKEY_MJD_DETECT_FIRST : head_calc[DATAKEY_MJD_DETECT_FIRST]
            }
            sel = util.select_subsample(args,var_dict)
            if sel is False :
                data_dict = {
                    'head_raw'  : head_raw,
                    'head_calc' : head_calc,
                    'select'    : False
                }
                return data_dict

        # - - - - - - -
        head_raw[DATAKEY_zHEL]       = table_head.REDSHIFT_HELIO[evt]
        head_raw[DATAKEY_zHEL_ERR]   = table_head.REDSHIFT_HELIO_ERR[evt]

        # strip off calculated values
        head_calc[DATAKEY_zCMB]          = table_head.REDSHIFT_FINAL[evt]
        head_calc[DATAKEY_zCMB_ERR]      = table_head.REDSHIFT_FINAL_ERR[evt]
        head_calc[DATAKEY_MWEBV]         = table_head.MWEBV[evt]
        head_calc[DATAKEY_MWEBV_ERR]     = table_head.MWEBV_ERR[evt]

        # - - - - - -
        # store HOSTGAL and HOSTGAL2 keys in head_raw[calc]
        self.store_hostgal(DATAKEY_LIST_RAW,  evt, head_raw ) # return head_raw
        self.store_hostgal(DATAKEY_LIST_CALC, evt, head_calc)

        # check for true sim type (sim or fakes), Nov 14 2021
        if SIMKEY_TYPE_INDEX in head_names:
            head_sim[SIMKEY_TYPE_INDEX] = table_head[SIMKEY_TYPE_INDEX][evt]
 
        # - - - - - - - - - - -
        # get pointers to PHOT table
        ROWMIN = table_head.PTROBS_MIN[evt]
        ROWMAX = table_head.PTROBS_MAX[evt]

        # note that there are NOBS+1 rows, but last row is pad word -777
        NOBS = ROWMAX - ROWMIN
        phot_raw = self.init_phot_dict(NOBS)

        table_column_names = table_phot.columns.names

        if 'FLT' in table_column_names:
            LEGACY_FLT = True  # legacy column name is FLT for band
        else:
            LEGACY_FLT = False

        for varname in varlist_obs:
            varname_table = varname
            if LEGACY_FLT:
                if varname == 'BAND' : varname_table = 'FLT'

            if varname_table in table_column_names :
                phot_raw[varname] = \
                    table_phot[varname_table][ROWMIN:ROWMAX].copy()

        # - - - - -
        # get field from from first observation,
        # Beware that event can overlap multiple fields.
        field = phot_raw[DATAKEY_FIELD][0]
        missing_field = (field == FIELD_NULL or field == FIELD_VOID )
        if missing_field  and args.survey == 'LSST' :
            field = self.field_plasticc_hack(table_dict['head_file'])

        head_raw[DATAKEY_FIELD] = field

        # - - - -
        spec_raw = {}

        # - - - - -
        # load output dictionary
        data_dict = {
            'head_raw'  : head_raw,
            'head_calc' : head_calc,
            'phot_raw'  : phot_raw,
            'spec_raw'  : spec_raw,
        }
        if len(head_sim) > 0:
            data_dict['head_sim'] =  head_sim
            
        if apply_select :
            data_dict['select'] = True

        return data_dict

        # end read_event

    def store_hostgal(self, datakey_list, evt, head_store):

        # store hostgal info in head_store dictionary.
        # Note that input head_store is modified here.

        table_dict = self.config_data['table_dict']
        table_head = table_dict['table_head']
        head_names = table_dict['head_names']

        len_base = len(HOSTKEY_BASE)
        for key in datakey_list:
            if HOSTKEY_BASE not in key: continue
            key2 = HOSTKEY_BASE + '2' + key[len_base:] # neighbor host
            key3 = HOSTKEY_BASE + '3' + key[len_base:]
            key_list = [ key, key2, key3]
            for k in key_list:
                if k in head_names :
                    head_store[k] = table_head[k][evt]

        # end store_hostgal


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

    def field_plasticc_hack(self,head_file_name):

        # ugly/embarassing hack to get field (DDF or WFD) from filename
        # because original plasticc data didn't store field.

        if FIELD_DDF in head_file_name:
            field = FIELD_DDF
        elif FIELD_WFD in head_file_name:
            field = FIELD_WFD
        else:
            msgerr.append(f"Unable to determine FIELD for")
            msgerr.append(f"{head_file_name}")
            util.log_assert(False,msgerr)

        return field
        # end field_plasticc_hack

    def set_dump_flag(self, isn, data_event_dict):
        d_raw = data_event_dict['head_raw']
        zhel  = d_raw['REDSHIFT_HELIO']
        dump_flag = False # isn<200 and zhel < 0.25
        return dump_flag
        # set_dump_flag


