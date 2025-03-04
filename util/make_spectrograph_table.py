#!/usr/bin/env python
#
# Created Feb 25 2025                
# Read 'wave flux fluxerr' from multiple files, presumably created by
# an ETC for a spectroscopic instrument,
# and create an SNANA-formatted SPECTROGRAPH table for the SNANA kcor.exe
# program and simulation.


import os, sys, yaml, argparse, logging, datetime
import pandas as pd
import numpy  as np

SNR_MIN = 0.1  # to write specbin, require SNR > SNR_MIN for all Texpose

tnow       = datetime.datetime.now()
DATE_STAMP = ('%4.4d-%2.2d-%2.2d' % (tnow.year,tnow.month,tnow.day) )
USERNAME   = os.environ['USER']


UNIT_CONVERT_DICT = {
    'nm'  : 10.0,
    'um'  : 10000.0 ,
    'A'   : 1.0
}

# define config key names
KEYNAME_INSTRUMENT      = 'INSTRUMENT'
KEYNAME_WAVE_R          = 'WAVE_R'
KEYNAME_WAVE_UNIT       = 'WAVE_UNIT'
KEYNAME_WAVE_BIN_SIZE   = 'WAVE_BIN_SIZE'
KEYNAME_OUTFILE_SNANA   = 'OUTFILE_SNANA'
KEYNAME_SEDFLUX_TABLES  = 'SEDFLUX_TABLES'


KEYLIST_DICT_REQUIRE = { 
    KEYNAME_INSTRUMENT     :   'name of instrument', 
    KEYNAME_WAVE_R         :   'resolving power, wave/FWHM_wave',  
    KEYNAME_WAVE_UNIT      :   'unit of input wavelength (nm, A, um)', 
    KEYNAME_WAVE_BIN_SIZE  :   'defines SPECTROGRAPH wavelength binning, A' ,
    KEYNAME_OUTFILE_SNANA  :   'name of output SPECTROGRAPH table file for SNANA sim', 
    KEYNAME_SEDFLUX_TABLES :   'list of flux tables'
}

# define comments for --HELP option
#COMMENT_KEY_LIST = [ 'name of instrument

KEYNAME_SPECBIN         = 'SPECBIN'  # for output spectrographi file

# =======================================                                                              
def get_args():
    parser = argparse.ArgumentParser()

    msg = "HELP on input config file"
    parser.add_argument("-H", "--HELP", help=msg, action="store_true")

    msg = "name of the yml config file to run"
    parser.add_argument("input_config", help=msg, nargs="?", default=None)

    args = parser.parse_args()

    if args.HELP:        print_help()
    args.input_config = os.path.expandvars(args.input_config)

    return args

def read_config(args):
    with open(args.input_config) as f:
        config = yaml.safe_load(f.read())

    for key in KEYLIST_DICT_REQUIRE.keys():
        if key not in config:
            sys.exit(f"\n ERROR: missing required config key {key}")
    return config
    # end read_config

def print_help():

    print('')
    print(f"# Description of input config file used for command")
    print(f"#    make_spectrograph_table.py <config_file>")
    print('')

    for key, comment in KEYLIST_DICT_REQUIRE.items():
        if 'SEDFLUX' not in key:
            key_plus_colon = key + ':'
            print(f"{key_plus_colon:<16}  <{comment}> ")

    # - - - - -
    print('')
    print('# The example below has separate blue and red spectra that are')
    print('# combined into a single SNR map in the output SPECTROGRAPH file,')
    print('# which results in a single blue+red spectrum from the simulation. ')
    print('# Beware that the Texpose range below must cover Texpose requests')
    print('# in the TEXPOSE argument of sim-input TAKE_SPECTRUM keys.')
    print('')
    print(f"{KEYNAME_SEDFLUX_TABLES}:   ")
    print(f"#   magref Texpose   spec_table_file")
    print(f"  -   19     500    spec_mag19_blue_texpose500.dat     # table contents: wave flux fluxerr")
    print(f"  -   19     500    spec_mag19_red_texpose500.dat")
    print(f"  -   22     500    spec_mag22_blue_texpose500.dat")
    print(f"  -   22     500    spec_mag22_red_texpose500.dat")
    print(f"#")
    print(f"  -   19    1200    spec_mag19_blue_texpose1200.dat")
    print(f"  -   19    1200    spec_mag19_red_texpose1200.dat")
    print(f"  -   22    1200    spec_mag22_blue_texpose1200.dat")
    print(f"  -   22    1200    spec_mag22_red_texpose1200.dat")
    print(f"#")
    print(f"  -   19    3000    spec_mag19_blue_texpose3000.dat")
    print(f"  -   19    3000    spec_mag19_red_texpose3000.dat")
    print(f"  -   22    3000    spec_mag22_blue_texpose3000.dat")
    print(f"  -   22    3000    spec_mag22_red_texpose3000.dat")
    print(f"#")
    print(f"#")
    print(f"# If you have SNR information for only a single exposure time,")
    print(f"# and want to scale SNR to other exposure times, re-use each flux-table")
    print(f"# with different snr scales:")
    print(f"#")
    print(f"SEDFLUX_TABLES:   ")
    print(f"#   magref Texpose   spec_table_file")
    print(f"  -   19      600    spec_snr_19.txt*0.577    # divide snr by sqrt(3) ")
    print(f"  -   21      600    spec_snr_21.txt*0.577")
    print(f"  -   19     1800    spec_snr_19.txt          # info only for this Texpose")
    print(f"  -   21     1800    spec_snr_21.txt          # info only for this Texpose" )
    print(f"  -   19     3600    spec_snr_19.txt*1.41     # multiply SNR by sqrt(2)")
    print(f"  -   21     3600    spec_snr_21.txt*1.41" )
    

    sys.exit(0)

def read_single_flux_table(flux_table, wave_scale ):

    # Read 3-column file of "wave flux fluxerr",
    # Next, scale wave.
    
    colname_wave    = 'wave'
    colname_flux    = 'flux'
    colname_fluxerr = 'fluxerr'
    colname_snr     = 'snr'

    df = pd.read_csv(flux_table, comment="#", delim_whitespace=True)
    key_list = df.keys().to_list()
    
    found_wave    = colname_wave in key_list[0].lower()
    found_flux    = colname_flux in key_list[1].lower()
    found_snr     = colname_snr  in key_list[1].lower()
    if len(key_list) > 2:
        found_fluxerr = colname_fluxerr in key_list[2].lower()

    if not found_wave:
        sys.exit(f"\n ERROR: did not find 'wave' in first column of {flux_table}")

    
    logging.info(f"\t Read {flux_table}")
    logging.info(f"\t Found columns {key_list}")
    
    wave_list     = df.iloc[:, 0].to_list()
    nbin          = len(wave_list)
    wave_list     = [wave * wave_scale for wave in wave_list]
    wave_min      = wave_list[0]
    wave_max      = wave_list[-1]

    if found_flux :
        flux_list     = df.iloc[:, 1].to_list()
        fluxerr_list  = df.iloc[:, 2].to_list()    
    elif found_snr:
        snr_list      = df.iloc[:, 1].to_list()
        flux_list     = [100.0] * nbin   # flux unit doesn't matter here
        fluxerr_list  = [x/(y+1.0e-9) for x,y in zip(flux_list,snr_list) ]
    else:
        sys.exit(f"\n ERROR: did not find {colname_flux} or {colname_snr} column in {flux_table}")

    
    flux_dict = {
        'nbin'           : nbin,
        'wave_list'      : wave_list,
        'flux_list'      : flux_list,
        'fluxerr_list'   : fluxerr_list,
        'wave_min'       : wave_min,
        'wave_max'       : wave_max
    }

    return flux_dict
    # end read_single_flux_table

def sanity_check(magref, texpose, flux_table):

    # abort on crazy value or missing file

    
    msgerr = None

    if magref < 10 or magref > 30 :
        msgerr = f"Insane value for magref={magref} " \
                 f"\n\t (flux_table={flux_table})"

    if texpose < 10 or texpose > 100000 :
        msgerr = f"Insane value for texpose={texpose}" \
                 f"\n\t (flux_table={flux_table})"

    if not os.path.exists(flux_table):
        msgerr  = f"cannot find flux-table file: \n\t {flux_table}"

    # - - - - - - 
    if msgerr:
        sys.exit(f"\nERROR: {msgerr}")

    return


def read_sedflux_tables(args, config, spectro_data):

    # read sed_flux tables and store rebinned contents in dictionary 
    # attached to config
    logging.info('')
    logging.info(f"# =====================================================")
    logging.info(f" read_sedflux_tables: read flux, flxuerr from tables")
    
    magref_dict = {}
    SEDFLUX_TABLES     = config[KEYNAME_SEDFLUX_TABLES]
    wave_bin_size      = config[KEYNAME_WAVE_BIN_SIZE]
    wave_unit          = config[KEYNAME_WAVE_UNIT]
    wave_scale         = UNIT_CONVERT_DICT[wave_unit]

    flux_dict_list = []
    wave_min_global = 9.0E9
    wave_max_global = 0.0
    
    texpose_dict = {}

    for row in SEDFLUX_TABLES:
        magref     = row.split()[0]
        texpose    = row.split()[1]
        flux_table = row.split()[2]

        # check for <flux_table>*<snr_scale>
        if '*' in flux_table:
            tmp        = flux_table.split('*')  # [ flux_table, snr_scale ]
            flux_table = tmp[0]
            snr_scale  = float(tmp[1])
        else:
            snr_scale  = 1.0
        
        flux_table = os.path.expandvars(flux_table)
        
        sanity_check( float(magref), float(texpose), flux_table)

        if magref not in texpose_dict :  texpose_dict[magref] = []
        texpose_dict[magref].append(texpose)

        logging.info(f"Process flux table for magref={magref} | Texpose={texpose} | snr_scale={snr_scale}: ")
        
        flux_dict = read_single_flux_table(flux_table, wave_scale)
        flux_dict['magref']  = magref
        flux_dict['texpose'] = texpose
        flux_dict['flux_table_file'] = os.path.basename(flux_table)
        flux_dict['snr_scale']       = snr_scale
        
        flux_dict_list.append(flux_dict)
        
        nbin     = flux_dict['nbin']
        wave_min = flux_dict['wave_min']
        wave_max = flux_dict['wave_max']
        
        logging.info(f"\t nbin_raw={nbin}  " \
                     f"wave[min,max] = {wave_min:.1f} to {wave_max:.1f} A")


        if wave_min < wave_min_global:  wave_min_global= wave_min
        if wave_max > wave_max_global:  wave_max_global= wave_max
        logging.info(f"")
        
    logging.info(f" Global wave[min,max] = " \
                 f"{wave_min_global:.1f} to {wave_max_global:.1f}")
    

    check_texpose(texpose_dict)

    spectro_data['flux_dict_list']  = flux_dict_list
    spectro_data['wave_min_global'] = wave_min_global
    spectro_data['wave_max_global'] = wave_max_global
    
    return
    # end read_sedflux_tables

def check_texpose(texpose_dict):

    # check that the same Texpose values are defined for each magref;
    # if not, then abort
 
    n_magref = len(texpose_dict)
    magref_list = list(texpose_dict.keys())
    if n_magref != 2:
        msgerr = f"\n ERROR: Found {n_magref} magref values {magref_list}; " \
                 f"must be exactly 2 values"
        sys.exit(msgerr)

    unique_texpose_dict = {}
    for magref, texpose_list in texpose_dict.items():
        unique_texpose_list = sorted(list(set(texpose_list)))
        unique_texpose_dict[magref] = unique_texpose_list

    magref0 = magref_list[0]
    magref1 = magref_list[1]
    texpose_list0 = unique_texpose_dict[magref0]
    texpose_list1 = unique_texpose_dict[magref1]

    if texpose_list0 != texpose_list1:
        print(f"\n\n")
        print(f"ERROR: Texpose lists are different for each magref:")
        print(f"\t magref={magref0} -> Texpose_list = {texpose_list0}")
        print(f"\t magref={magref1} -> Texpose_list = {texpose_list1}")
        sys.exit("\n ABORT")

    return

def rebin_sedflux_tables(args, config, spectro_data):

    # rebin on user-define wave grid, and sum flux,fluxvar
    # for files with same magref and texpose; e.g. sum blue and red
    # regions of a spectrograph
    
    wave_bin_size        = config[KEYNAME_WAVE_BIN_SIZE]

    logging.info('')
    logging.info(f"# =====================================================")
    logging.info(f" rebin_sedflux_tables: " \
                 f"rebin flux,fluxerr on {wave_bin_size:.0f} A grid")

    flux_dict_list  = spectro_data['flux_dict_list'] 
    wave_min_global = spectro_data['wave_min_global']
    wave_max_global = spectro_data['wave_max_global'] 

    ibin_min      = int(wave_min_global/wave_bin_size)
    ibin_max      = int(wave_max_global/wave_bin_size) + 1
    nbin_grid     = ibin_max - ibin_min
    wave_min_grid = ibin_min * wave_bin_size
    wave_max_grid = ibin_max * wave_bin_size

    logging.info(f" output wave grid: {wave_min_grid} to {wave_max_grid}  "\
                 f"nbin_grid={nbin_grid}")
    
    # define lower edge of each wave bin
    wave_grid_array = np.linspace(wave_min_grid, wave_max_grid, nbin_grid,
                                  endpoint=False)
    
    #print(f"\n xxx wave_bin_array = \n{wave_bin_array}")
    flux_rebin_dict_list = []
    key_unique_list      = []
    
    for flux_dict in flux_dict_list:
        magref          = flux_dict['magref']
        texpose         = flux_dict['texpose']
        flux_table_file = flux_dict['flux_table_file']
        snr_scale       = flux_dict['snr_scale']
        
        wave_list    = flux_dict['wave_list']
        flux_list    = flux_dict['flux_list']
        fluxerr_list = flux_dict['fluxerr_list']
        fluxvar_list = [x*x for x in fluxerr_list]

        key_unique         = str(magref) + '_' + str(texpose)
        flux_grid_array    = [ 0.0 ] * nbin_grid
        fluxvar_grid_array = [ 0.0 ] * nbin_grid 
        
        for wave, flux, fluxvar in zip(wave_list, flux_list, fluxvar_list):
            ibin = int((wave - wave_min_grid)/wave_bin_size)
            flux_grid_array[ibin]    += flux
            fluxvar_grid_array[ibin] += fluxvar
            
        flux_rebin_dict = {
            'nbin'           : nbin_grid,
            'wave_list'      : wave_grid_array,
            'flux_list'      : flux_grid_array,
            'fluxvar_list'   : fluxvar_grid_array
        }


        if key_unique in key_unique_list:
            # add flux and fluxvar to existing flux_rebin_dict
            j = key_unique_list.index(key_unique)
            tmp0 = flux_rebin_dict_list[j]['flux_list']
            tmp1 = flux_rebin_dict['flux_list']
            flux_rebin_dict_list[j]['flux_list']  =  \
                [ x + y for x, y in zip(tmp0, tmp1)]

            tmp0 = flux_rebin_dict_list[j]['fluxvar_list']
            tmp1 = flux_rebin_dict['fluxvar_list']
            flux_rebin_dict_list[j]['fluxvar_list']  =  \
                [ x + y for x, y in zip(tmp0, tmp1)]

            flux_rebin_dict_list[j]['flux_table_file'] += f"+{flux_table_file}"
            logging.info(f"\t append {key_unique} for {flux_table_file} ")            
        else:
            # add new dictionary
            flux_rebin_dict['key_unique'] = key_unique
            flux_rebin_dict['magref']     = magref
            flux_rebin_dict['texpose']    = texpose
            flux_rebin_dict['snr_scale']  = snr_scale            
            flux_rebin_dict['flux_table_file'] = flux_table_file
            logging.info(f"\t load   {key_unique} for {flux_table_file} ")
            key_unique_list.append(key_unique)
            flux_rebin_dict_list.append(flux_rebin_dict)

    # - - - - - 
    # compute snr_list for each dictionary    
    for j, flux_rebin_dict in enumerate(flux_rebin_dict_list):
        flux_list    = flux_rebin_dict['flux_list']
        fluxvar_list = flux_rebin_dict['fluxvar_list']
        snr_scale    = flux_rebin_dict['snr_scale'] 
        snr_list     = [ snr_scale*x/np.sqrt(y+1.0E-9) for x, y in zip(flux_list,fluxvar_list) ]
        flux_rebin_dict_list[j]['snr_list'] = snr_list

    
    spectro_data['flux_rebin_dict_list'] = flux_rebin_dict_list
    return
    # end rebin_sedflux_tables

def open_outfile(args, config):

    # open spectrograph output file, write DOCANA, and write some header info
    
    outfile = config['OUTFILE_SNANA']
    instr   = config[KEYNAME_INSTRUMENT]
    CWD     = os.getcwd()
    command = ' '.join(sys.argv)
    
    logging.info('')
    logging.info('# ===========================================================')    
    logging.info(f" Open output spectrograph file for SNANA sim: {outfile}")

    o = open(outfile,"wt")

    o.write(f'DOCUMENTATION:\n')
    o.write(f'  PURPOSE:   SNR-vs-wavelength and Texpose; to simulate {instr} spectra:\n')
    o.write(f'  USAGE_KEY:      SPECTROGRAPH \n')
    o.write(f'  USAGE_CODE:     kcor.exe \n')
    o.write(f'  CREATE_DATE:    {DATE_STAMP}\n')
    o.write(f'  CREATE_USER:    {USERNAME} \n')
    o.write(f'  CREATE_COMMAND: {command} \n')    
    o.write(f'  CREATE_DIRNAME: {CWD} \n')
    o.write(f'DOCUMENTATION_END:\n')

    o.write('\n\n')


    o.write(f'INSTRUMENT:  {instr}\n')

    return o
    
def write_specbins(o, args, config, specto_data):

    flux_rebin_dict_list = specto_data['flux_rebin_dict_list']

    magref_unique = []
    texpose_unique = []
    wave_list      = flux_rebin_dict_list[0]['wave_list']
    
    # make list of unique magref and texpose    
    for flux_rebin_dict in flux_rebin_dict_list:
        magref  = flux_rebin_dict['magref']
        texpose = flux_rebin_dict['texpose']
        flux_table_file = flux_rebin_dict['flux_table_file']

        if magref not in magref_unique:
            magref_unique.append(str(magref))
            
        if texpose not in texpose_unique:
            texpose_unique.append(str(texpose))

    string_magref  = "MAGREF_LIST:  "  + ' '.join(magref_unique) + \
        '      # defines SNR0 & SNR1'
    string_texpose = "TEXPOSE_LIST: "  + ' '.join(texpose_unique) 
    o.write(f"{string_magref}\n")
    o.write(f"{string_texpose}\n")
    o.write(f"SNR_POISSON_RATIO_ABORT_vsTEXPOSE: 2.0  " \
            f"# allow SNR ~ sqrt[Texpose] within x2 factor\n")
    o.write(f"\n")    

    # - - - - - -
    # write spec bin info
    wave_bin_size = config[KEYNAME_WAVE_BIN_SIZE]
    wave_r        = config[KEYNAME_WAVE_R]
    str_magref0   = str(magref_unique[0])
    str_magref1   = str(magref_unique[1])

    comment_snr = '       '.join(['SNR0   SNR1' ] * len(texpose_unique))
    o.write(f"#          wave_min   wave_max wave_res   {comment_snr}\n")
    
    for iwave, wave_lo in enumerate(wave_list) :
        wave_hi = wave_lo + wave_bin_size

        wave_res = (wave_lo/wave_r)/2.235  # convert FWHM to sigma
        line = f"{KEYNAME_SPECBIN}:  {wave_lo:9.3f}  {wave_hi:9.3f}  {wave_res:6.2f}  "

        n_snr_cut = 0
        for texpose in texpose_unique:
            string_snr_dict = {str_magref0 : '' , str_magref1 : '' }
            for flux_rebin_dict in flux_rebin_dict_list :
                magref_tmp  = flux_rebin_dict['magref']
                texpose_tmp = flux_rebin_dict['texpose']
                if texpose_tmp == texpose :
                    snr     = flux_rebin_dict['snr_list'][iwave]
                    if snr <= SNR_MIN : n_snr_cut += 1 
                    string_snr_dict[ str(magref_tmp)]  += f"{snr:6.3f} "
                    
            line += string_snr_dict[str_magref0]
            line += string_snr_dict[str_magref1]
            line += '    '
        if n_snr_cut == 0:
            o.write(f"{line}\n")
        
    #sys.exit(f"\n xxx flux_rebin_dict_list = \n{flux_rebin_dict_list}")
    
    return
    # end write_specbins
    
# =============================================
#       MAIN                                                                     
# =============================================

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO,
                        format="[%(levelname)4s| %(message)s")

    args   = get_args()

    config = read_config(args)

    spectro_data = {}
    
    read_sedflux_tables(args, config, spectro_data)

    rebin_sedflux_tables(args, config, spectro_data)

    o = open_outfile(args, config)    
    write_specbins(o, args, config, spectro_data)

    o.close()

    # === END: ===



    
