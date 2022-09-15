#
# Created July 2022 by R.Kessler
# Write data files in csv format similar to data released
# on Kaggle platform for plasticc challenge.
# Original intent is to release ELASTICC data to public, but this
# csv format may be useful for future public DRs.
#
# TO DO:
#  + write detect flag (auto read from FITS header?)
#  - provide separate table of flux correction vs. MWEBV and band
#  + flush csv files
#  - sort each file by SNID

import datetime, glob, yaml
import logging
import os, sys, math
import shutil, subprocess
import numpy as np
import makeDataFiles_params as gpar
import makeDataFiles_util as util

# ==============================================

class csvWriter:
    METADATA_FILENAME      = "metadata.csv"
    LIGHTCURVE_FILE_PREFIX = "lightcurve"

    def __init__( self, args, config_data):

        # store info passed from base
        self.args = args           
        self.config_data = config_data

        if args.merge :  return

        # - - - - - - - - - - - - - - - -
        # create output folder
        outdir = args.outdir_csv
        util.create_output_folder(outdir)

        # store file pointer for each possible light curve csv file,        
        n_unit = len(config_data['data_unit_name_list'])
        fp_lc_list = [ None ] * n_unit
        self.fp_lc_list = fp_lc_list

        # store map between SNANA varname and output csv name
        self.mapvar_meta = {
            gpar.DATAKEY_SNID            : ['snid',  '9', 
                                            'SN identifier'] ,
            gpar.DATAKEY_NOBS            : ['nobs',  '3',
                                            'Nobs in lightcurve file' ] ,
            gpar.DATAKEY_SNTYPE          : ['type',  '2d',
                                            'spectroscopic type (integer)'] ,
            gpar.DATAKEY_RA              : ['ra',    '.6f',
                                            'Right Ascension (degrees)' ] ,
            gpar.DATAKEY_DEC             : ['dec',   '.6f',
                                            'Declination (degrees)' ] ,
            gpar.DATAKEY_MWEBV           : ['mwebv', '.3f',
                                            'Galactic E(B-V)' ] ,
            gpar.HOSTKEY_SPECZ           : ['hostgal_zspec',      '.3f',
                                            'host spec-redshift' ] ,
            gpar.HOSTKEY_SPECZ_ERR       : ['hostgal_zspec_err',  '.3f',
                                            'error on host spec-z' ] ,
            gpar.HOSTKEY_PHOTOZ          : ['hostgal_zphot',      '.3f',
                                            'host photo-z' ] ,
            gpar.HOSTKEY_PHOTOZ_ERR      : ['hostgal_zphot_err',  '.3f',
                                            'error on host photo-z' ] ,
            gpar.HOSTKEY_SNSEP           : ['hostgal_snsep',      '.3f',
                                            'host-SN sep, arcsec' ] 
        }

        self.mapvar_lc = {
            gpar.DATAKEY_SNID         : [ 'snid',         '9',
                                          'SN identifier' ],
            gpar.DATAKEY_MJD          : [ 'mjd'  ,        '.4f',
                                          'Modified Julien Date (days)' ],
            gpar.DATAKEY_BAND         : [ 'band' ,        '1',
                                          'single-char passband; e.g., r' ],
            gpar.DATAKEY_PHOTFLAG     : [ 'detect',       '1',
                                          '0=not detcted, 1=detected' ],
            gpar.DATAKEY_FLUXCAL      : [ 'fluxcal' ,     '11.4e',
                                          'calibrated flux with ZP=27.5' ],
            gpar.DATAKEY_FLUXCALERR   : [ 'fluxcal_err' , '10.4e',
                                          'uncertainty on fluxcal' ]
        }

        # there is only one metadata file, so create it now during init
        # and save pointer so we don't have to close and reopen this file
        meta_file = f"{outdir}/{self.METADATA_FILENAME}"
        logging.info(f" Open meta file: {meta_file}")
        self.fp_meta = open(meta_file,"wt")

        # write header
        self.write_csv_header( self.fp_meta, self.mapvar_meta)

        # init a few counters
        self.nevt_all_write   = 0
        self.nevt_type_write  = 0
        self.nobs_write       = 0
        self.readme_file = f"{outdir}/DATA.README"

        return

    # end __init__

    def append_readme_vardef_csv(self):

        # append variable definitions to readme file
        readme_file = self.readme_file

        logging.info(f" Define variables in {readme_file}")

        # write with yaml format
        with open(readme_file,"at") as r:
            
            r.write("\nVARDEF_METADATA:\n")
            for var_snana, var_csv in self.mapvar_meta.items() :
                var     = var_csv[0] + ':'
                comment = var_csv[2]
                r.write(f"  {var:<20} {comment} \n")

            r.write("\nVARDEF_LIGHTCURVE:\n")
            for var_snana, var_csv in self.mapvar_lc.items() :
                var     = var_csv[0] + ':'
                comment = var_csv[2]
                r.write(f"  {var:<20} {comment} \n")

            r.write('\n')
        return
        # end append_readme_vardef_csv


# =============================================
    def write_csv_header(self, fp, mapvar):
        # write csv header to file.
        # Inputs: 
        #   fp= file pointer
        #   mapvar = map between SNANA varname and csv varname

        line_header = ""
        for var_snana, var_csv in mapvar.items() :
            var = var_csv[0]
            line_header += f"{var}, "

        line_header = line_header[:-2]  # remove last comma
        fp.write(f"{line_header}\n")
        fp.flush()
        return

        # end write_csv_header

# ========================================================
    def write_event_csv(self, data_event_dict):

        # Inputs:
        #  data_event_dict: event dictionary

        outdir          = self.args.outdir_csv

        data_unit_name_list   = self.config_data['data_unit_name_list']
        data_unit_nevent_list = self.config_data['data_unit_nevent_list']
        index_unit  = data_event_dict['index_unit']
        # xx ?? nevent      = data_unit_nevent_list[index_unit]
        
        if self.fp_lc_list[index_unit] is None:
            prefix = self.LIGHTCURVE_FILE_PREFIX
            lc_file = f"{outdir}/{prefix}{index_unit:02d}.csv" 
            logging.info(f" Open lc_file: {lc_file}") 
            fp = open(lc_file,"wt")
            self.fp_lc_list[index_unit] = fp
            self.write_csv_header( fp, self.mapvar_lc) # write header

        # update metadata 
        fp = self.fp_meta
        self.append_metadata_csv(fp, data_event_dict)

        # update lightcurve data
        fp = self.fp_lc_list[index_unit]
        self.append_lightcurve_csv(fp, data_event_dict)
        
        return

    # end write_event_csv

# ====================================================
    def append_metadata_csv(self, fp, data_event_dict):

        head_raw  = data_event_dict['head_raw']
        head_calc = data_event_dict['head_calc']
        phot_raw  = data_event_dict['phot_raw']  # for NOBS
        msgerr = []
        line = ""

        for var_snana, var_csv  in self.mapvar_meta.items() :
            if var_snana in head_raw :
                value = head_raw[var_snana]
            elif var_snana in head_calc :
                value = head_calc[var_snana]                
            elif var_snana in phot_raw :
                value = phot_raw[var_snana]
            else:
                msgerr.append(f" Cannot find SNANA header var {var_snana}")
                util.log_assert(False,msgerr)

            fmt = var_csv[1]
            if value == -9.0 : 
                value = -9 ;  fmt='d'

            line += f"{value:{fmt}}, "

        line = line[:-2] # remove last comma
        fp.write(f"{line}\n")
        fp.flush()
        self.nevt_all_write += 1
        
        return
        # end append_metadata_csv

# ====================================================
    def append_lightcurve_csv(self, fp, data_event_dict):

        head_raw = data_event_dict['head_raw']
        phot_raw = data_event_dict['phot_raw']
        nobs     = phot_raw[gpar.DATAKEY_NOBS]
        snid     = head_raw[gpar.DATAKEY_SNID]
        fmt_snid        = self.mapvar_lc[gpar.DATAKEY_SNID][1]
        photflag_detect = self.args.photflag_detect

        lines = ''
        for o in range(0,nobs):
            line_tmp = f"{snid:{fmt_snid}}, "
            for var_snana, var_csv  in self.mapvar_lc.items() :
                if var_snana == gpar.DATAKEY_SNID: continue

                value = phot_raw[var_snana][o]
                fmt   = var_csv[1]

                # replace photflag with detect = 0 or 1 flag
                if var_snana == gpar.DATAKEY_PHOTFLAG:
                    detect = (value & photflag_detect)/photflag_detect
                    value  = int(detect) ;  fmt   = '1'

                line_tmp += f"{value:{fmt}}, "
            line_tmp = line_tmp[:-2]  # remove last comma
            lines += f"{line_tmp}\n"
            self.nobs_write += 1

        # write all photometry lines with one call
        fp.write(f"{lines}")
        fp.flush()

        return
        # end append_lightcurve_csv

    def end_write(self, index_unit):

        # close file pointer for this data unit
        self.fp_lc_list[index_unit].close()
        
        # continue only for 1st data unit
        if index_unit != 1 :
            return

        # - - - - - - 
        readme_stats_sum = self.config_data['readme_stats_sum']
        readme_dict = {
            'readme_file'  : self.readme_file,
            'readme_stats' : readme_stats_sum,
            'data_format'  : gpar.FORMAT_CSV,
            'docana_flag'  : True
        }
        util.write_readme(self.args, readme_dict)

        # append csv variable definitions
        self.append_readme_vardef_csv()

        nevt = self.nevt_all_write
        nobs = self.nobs_write        
        logging.info(f"\n csv-write summary:")
        logging.info(f"\t Finished writing {nevt} events to csv files.")
        logging.info(f"\t Total number of observations: {nobs}")
        self.fp_meta.close()

        return
        # end of end_read_data

# =======================================================
    def merge_csv_driver(self):
        outdir = self.args.outdir_csv
            
        print(f" xxx hello from merge. outdir = {outdir}")
        return
    # end merge_csv_driver
