#!/usr/bin/env python
#
# Created June 29 2023
#
# Author: Amanda Wasserman
#
# Purpose: Resimulate a given SN with its spectrum (with a user specified spectrograph) 
#          for a specific MJD utilizing the original dump file and input file.
#

import argparse, yaml, sys, os, re, glob, logging, pandas as pd

def parse_yaml(args):
    # parse yaml file, return dictionary of values

    path = args.file_name

    with open(path, "r") as stream:
        try:
            config_dic = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            logging.error("%s",exc)
    return(config_dic)

def get_args():
    # parse command line arguments, return arguments

    parser = argparse.ArgumentParser()
    
    msg = "yaml file name"
    msg1 = "flag to not resimulate"
    msg2 = "increase output verbosity"
    msg3 = "turn off noise in spectra"
    msg4 = "print help menu for config parameters"
    msg5 = "extra sim args to append to snlc_sim.exe"
    parser.add_argument("file_name", help=msg, type=str, default='')
    parser.add_argument("-n", "--no_resim", help=msg1, action='store_true', default=False)
    parser.add_argument("-v", "--verbose", help=msg2, action="store_true", default=False)
    parser.add_argument("-p", "--perfect", help=msg3, action="store_true", default=False)
    parser.add_argument("-H", "--help_config", help=msg4, action="store_true", default=False)
    parser.add_argument("-e", "--extra_sim_args",help=msg5, type=str, default='')
    
    args = parser.parse_args()
    Help = args.help_config
    if Help==True:
        print_help()
        sys.exit()

    file = args.file_name
    if file == '':
        logging.error('No yaml file name given.')
        sys.exit()

    no_resim = args.no_resim
    verbose = args.verbose
    perfect = args.perfect
    extra_sim_args = args.extra_sim_args
    config_dic = parse_yaml(args)

    if len(sys.argv) == 1:
        parser.print_help(); sys.exit()
        
    if verbose==True:
        logging.info('Retrieved arguments.')

    return no_resim, verbose, perfect, config_dic, Help, extra_sim_args



def find_dump_file(cid, simfolder):
    #find and return dump file for a given cid
    exp_var = os.path.expandvars(simfolder)
    all_dump_files = exp_var+'/*DUMP'
    list_dumps = glob.glob(all_dump_files)
    if len(list_dumps)==0:
        logging.error('No dump files found in %s', exp_var)
        sys.exit()
    counter = 0
    for dump in list_dumps:
        if "REPEAT" not in dump:
            with open(dump) as file:
                data = file.read()
            if re.search(str(cid) + " ",data):
                dump_path = dump
                counter += 1
    if counter !=1:
        logging.error('Number of dump files found with the cid is expected 1, but is: %s', counter)
        sys.exit()

    return(dump_path)

def retrieve_values_from_dump(dump_row, dump_key):
    #retrieve and return values from dump file for a given cid and key
    #if verbose == True:
    #    logging.info('Accessing dump file: %s', dump_file_path)
    #command = f"get_fitres_values.py -f {dump_file_path} -c {cid} -v "+dump_key+" -o out_simgen_resim.dump"
    #if verbose == True:
    #    logging.info('Getting fitres values using command: %s ', command)
    #os.system(command) 
    #df = pd.read_csv('out_simgen_resim.dump', comment='#', sep=r'\s+')
    value = str(dump_row[dump_key].values[0])
    return(value)

def is_it_Ia(lines):
    #check if the input file is for a Ia. Returns true if it is, false if it isn't
    word_to_check = 'SALT'

    result = any(word_to_check in line for line in lines)

    return(result)

def get_genversion(dump_file_path):
    # extract genversion string from simgen-dump file path (RK)
    genversion = (((dump_file_path.split("/"))[-1]).split("."))[0]
    return genversion
    # end get_genversion

def edit_input(cid, mjd, zcmb, peak_mjd, file_name, dump_file_path, simgen_eazy, hostlib_file, kcor, outdir, dump_cid_row):
    #edit input file for a given cid and mjd
    #returns path of new input file

    version = get_genversion(dump_file_path)
    # xxx mark delete RK version = file_name[file_name.find('WFD_')+4:file_name.find('_MODEL0')]  #change to be more general for SNANA, WFD_ is for wide fast deep

    ACTION_CHANGE = 'CHANGE'
    ACTION_REMOVE = 'REMOVE'
    ACTION_ADD = 'ADD'
    KCOR = kcor
    TAKE_SPECTRUM_VALUE = 'MJD('+str(mjd)+':'+str(mjd+2)+') TEXPOSE_ZPOLY(1200)'
    TAKE_HOST = 'HOST                     TEXPOSE_ZPOLY(1200)'
    GENVERSION_VALUE = 'ADD_SPEC_'+version+'_cid'+str(cid)+'_mjd'+str(mjd)    
    CIDOFF_VALUE = str(cid-1) 
    GENRANGE_REDSHIFT_VALUE = str(zcmb)+' '+str(zcmb)
    GENRANGE_PEAKMJD_VALUE = str(peak_mjd)+' '+str(peak_mjd)
    SIMGEN_EAZY = simgen_eazy
    HOSTLIB_FILE = hostlib_file   

    KEY_CHANGE_DICT = {
            # sim-input key                  simgen dump key        action           value                          only for Ia?
            'CIDRAN_MIN:'                  :    [  None,              ACTION_REMOVE,   None,                           False],
            'CIDRAN_MAX:'                  :    [  None,              ACTION_REMOVE,   None,                           False],
            'GENMAG_SMEAR_MODELNAME:'      :    [  None,              ACTION_REMOVE,   None,                           False],
            'LENSING_PROBMAP_FILE:'        :    [  None,              ACTION_REMOVE,   None,                           False],
            'SIMLIB_NREPEAT:'              :    [  None,              ACTION_REMOVE,   None,                           False],
            'GENSIGMA_VPEC:'               :    [  None,              ACTION_REMOVE,   None,                           False],
            'VPEC_ERR:'                    :    [  None,              ACTION_REMOVE,   None,                           False],
            'HOSTLIB_WGTMAP_FILE:'         :    [  None,              ACTION_REMOVE,   None,                           False],
            'SEARCHEFF_PIPELINE_EFF_FILE:' :    [  None,              ACTION_REMOVE,   None,                           False],
            'NON1A_KEYS:'                  :    [  None,              ACTION_REMOVE,   None,                           False],
            'NON1A:'                       :    [  None,              ACTION_REMOVE,   None,                           False],
            'SIMLIB_IDLOCK:'               :    [  'LIBID',           ACTION_ADD,      None,                           False],    
            'TAKE_SPECTRUM:'               :    [  None,              ACTION_ADD,      TAKE_HOST,                      False],
            'TAKE_SPECTRUM: '              :    [  None,              ACTION_ADD,      TAKE_SPECTRUM_VALUE,            False],
            'HOSTLIB_GALID_FORCE:'         :    [ 'GALID',            ACTION_ADD,      None,                           False],
            'GENMAG_OFF_GLOBAL:'           :    [ 'MAGSMEAR_COH',     ACTION_ADD,      None,                           False],
            'HOSTLIB_SPECBASIS_FILE:'      :    [  None,              ACTION_ADD,      SIMGEN_EAZY,                    False],
            'GENPREFIX:'                   :    [  None,              ACTION_ADD,      'ADD_SPEC',                     False], 
            'GENVERSION:'                  :    [  None,              ACTION_CHANGE,   GENVERSION_VALUE,               False],
            'NGENTOT_LC:'                  :    [  None,              ACTION_CHANGE,   '1',                            False],
            'FORMAT_MASK:'                 :    [  None,              ACTION_CHANGE,   '2',                            False],
            'CIDOFF:'                      :    [  None,              ACTION_CHANGE,   CIDOFF_VALUE,                   False],
            'KCOR_FILE:'                   :    [  None,              ACTION_CHANGE,   KCOR,                           False],
            'GENRANGE_REDSHIFT:'           :    [  None,              ACTION_CHANGE,   GENRANGE_REDSHIFT_VALUE,        False],
            'GENRANGE_PEAKMJD:'            :    [  None,              ACTION_CHANGE,   GENRANGE_PEAKMJD_VALUE,         False],
            'GENPEAK_SALT2x1:'             :    [  'SALT2x1',         ACTION_CHANGE,   None,                           True ],
            'GENSIGMA_SALT2x1:'            :    [  None,              ACTION_CHANGE,   '0 0',                          True ],
            'GENPEAK_SALT2c:'              :    [  'SALT2c',          ACTION_CHANGE,   None,                           True ],  
            'GENSIGMA_SALT2c:'             :    [  None,              ACTION_CHANGE,   '0 0',                          True ],
            'HOSTLIB_FILE:'                :    [  None,              ACTION_CHANGE,   HOSTLIB_FILE,                   False],
            }

    if perfect==True:
        KEY_CHANGE_DICT['SPECTROGRAPH_OPTMASK:'] = [  None,              ACTION_ADD,      '32768',                        False]
        logging.info('Perfect flag was raised, so turn off noise in spectrum.')

    lines = []
    with open(file_name, 'r') as f:
        for line in f:
            lines.append(line.strip())
    
    is_Ia = is_it_Ia(lines)
    
    keys = list(KEY_CHANGE_DICT.keys())
    final_lines = []

    #go through each line
    for line in lines:
        #if nothing changes the new string is the same as the original
        new_string = line
        for key in keys:
            #check if one of the keys we need to change/remove is in the line
            if key in line:
                new_string = line.split(':')[0]+":  "
                #if the action is remove, remove the line from the input file
                if KEY_CHANGE_DICT[key][1] == ACTION_REMOVE:
                    new_string = "#"+ key+" REMOVED BY SCRIPT"
                    if verbose==True:
                        logging.info('REMOVED BY SCRIPT: %s', key[0:-1])
                #if the action is change
                elif KEY_CHANGE_DICT[key][1] == ACTION_CHANGE:
                    #if the key is associated with a value and the SN is Ia, input the SALT values into the line
                    if KEY_CHANGE_DICT[key][0] == None and KEY_CHANGE_DICT[key][3]==is_Ia:
                        new_string = new_string + KEY_CHANGE_DICT[key][2] + "   # CHANGED BY SCRIPT"
                        if verbose==True:
                            logging.info('CHANGED BY SCRIPT: %s', key[0:-1])
                    #if the key is associated with a value, input that into the line
                    elif KEY_CHANGE_DICT[key][0] == None and KEY_CHANGE_DICT[key][3]==False:
                        new_string = new_string + KEY_CHANGE_DICT[key][2] + "   # CHANGED BY SCRIPT"
                        if verbose==True:
                            logging.info('CHANGED BY SCRIPT: %s', key[0:-1])      
                    #otherwise the key is associated with a dump file key, so retrieve that value from the dump file
                    #SALT values if Ia and then the rest
                    elif KEY_CHANGE_DICT[key][3]==is_Ia:
                        from_dump = retrieve_values_from_dump(dump_cid_row,KEY_CHANGE_DICT[key][0])
                        new_string = new_string + from_dump + "   # CHANGED BY SCRIPT"
                        if verbose==True:
                            logging.info('CHANGED BY SCRIPT: %s', key[0:-1]) 
                    else:
                        from_dump = retrieve_values_from_dump(dump_cid_row,KEY_CHANGE_DICT[key][0])
                        new_string = new_string + from_dump + "   # CHANGED BY SCRIPT"
                        if verbose==True:
                            logging.info('CHANGED BY SCRIPT: %s', key[0:-1])   
        final_lines.append(new_string)
    
    #add the action change keys
    for key in keys:
        value = KEY_CHANGE_DICT[key]
        if value[1] == ACTION_ADD:
            if value[0] == None:
                new_string = key + '  ' + value[2] + "   # ADDED BY SCRIPT"
                if verbose==True:
                    logging.info('ADDED BY SCRIPT: %s', key)
                final_lines.append(new_string)
            else:
                from_dump = retrieve_values_from_dump(dump_cid_row,KEY_CHANGE_DICT[key][0])
                new_string = key + ' ' + from_dump + "   # ADDED BY SCRIPT"
                if verbose==True:
                    logging.info('ADDED BY SCRIPT: %s', key)
                final_lines.append(new_string)
    
    #add the nonIa keys
    if is_Ia==False:
        new_string = 'NON1A_KEYS:     5 INDEX WGT MAGOFF MAGSMEAR SNTYPE   # ADDED BY SCRIPT'
        final_lines.append(new_string)


        nonIa_index = retrieve_values_from_dump(dump_cid_row, 'NON1A_INDEX')
        gentype = retrieve_values_from_dump(dump_cid_row, 'GENTYPE')

        new_string = 'NON1A:    '+ nonIa_index + ' 1 0 0  ' + gentype + '  # ADDED BY SCRIPT'
        final_lines.append(new_string)
        if verbose==True:
            logging.info('ADDED BY SCRIPT: NON1A_KEYS, NON1A')

    if verbose==False:
        logging.info('Keys have been changed, removed, and added as needed.')

    filename = version+"_cid"+str(cid)+'_mjd'+str(mjd)+"_new.input"   # name of the file to create 
    logging.info('Final input file name stored in %s/%s',outdir,filename)
    #join the sub-directory and filename using os.path.join() 

    path = os.path.join(outdir, filename) 
    # create the file using the built-in open() function 
    
    with open(path, 'w') as f: 
    #with open(filename, 'w') as f: 
        for i in final_lines:
            f.write(i+'\n')
    return(path)
    #return(filename)

def process_cid(cid,mjds,config_dic):
    #cid (string):
        #cid of SN to resimulate
    #mjds (list of strings):
        #list of mjds or 'peak' to simulate
    #config_dic (dictionary):
        #dictionary of config parameters

    #creates input file and runs simulation (unless no_resim flag is raised)

    simfolder = config_dic['SIMFOLDER']
    kcor = config_dic['KCOR_PATH']
    simgen_eazy = config_dic['SIMGEN_EAZY']
    hostlib_file = config_dic['HOSTLIB_FILE']
    if 'OUTDIR_SIM_INPUT' in config_dic:
        outdir = config_dic['OUTDIR_SIM_INPUT']
    else:
        outdir = 'store_resim_inputs'

    if outdir == 'store_resim_inputs' and not os.path.exists(outdir):
            os.makedirs(outdir)

    if verbose==True:
        logging.info('Simfolder:  %s', simfolder)
        logging.info('kcor file: %s', kcor)
        logging.info('simgen_eazy file: %s', simgen_eazy)
        logging.info('hostlib file: %s', hostlib_file)

    dump_file_path = find_dump_file(cid, simfolder)
    dump_df = pd.read_csv(dump_file_path, comment='#', sep=r'\s+')
    dump_cid_row = dump_df.loc[dump_df['CID'] == int(cid)]
    #store z and peak_mjd
    zcmb = retrieve_values_from_dump(dump_cid_row,'ZCMB')
    if verbose==True:
        logging.info('ZCMB is: %s', zcmb)
    peak_mjd = retrieve_values_from_dump(dump_cid_row,'PEAKMJD')
    if verbose==True:
        logging.info('Peak MJD is: %s', peak_mjd)

    #convert string 'peak' to peak_mjd
    if 'peak' in mjds:
        mjds.remove('peak')
        mjds.append(peak_mjd)

    logging.info('Processing for CID: %s', cid)
    mjds = [float(i) for i in mjds]
    cid = int(cid)
    logging.info('Processing for the following MJDs: %s', mjds)

    for date in mjds:
        #create file for each date
            #name structure for inpute file: version+"_cid"+str(cid)+'_mjd'+str(mjd)+"_new.input

        #create input file
        cwd = os.getcwd()

        # extract version string from simgen-dump file
        genversion = get_genversion(dump_file_path)

        command = f"quick_commands.py -v {genversion} --extract_sim_input"

        if outdir == 'store_resim_inputs' :
            os.chdir(cwd+ "/store_resim_inputs")
            os.system(command)
            os.chdir(cwd)
            file_name = cwd+"/store_resim_inputs/sim_input_"+genversion+"_MODEL0.input"
        else:
            os.chdir(outdir)
            os.system(command)
            os.chdir(cwd)
            file_name = outdir+ "/sim_input_" + genversion + "_MODEL0.input"
        
        if verbose==True:
            logging.info('Run the following command in %s: %s', outdir, command)

        
        logging.info('Input file to edit: %s', file_name)
        #edit input file
        if verbose==True:
            logging.info('Create edited input file.')

        path = edit_input(cid, date, zcmb, peak_mjd, file_name, dump_file_path, simgen_eazy, hostlib_file, kcor, outdir,dump_cid_row) 

        if no_resim==False:
            #run sim
            logging.info('Running simulation.')
            command = f"snlc_sim.exe "+path.split('/')[1]+extra_sim_args
            logging.info('PATH: %s', path)
            logging.info('Run the following command: %s', command)
            if outdir == 'store_resim_inputs':
                os.chdir(cwd+ "/store_resim_inputs")
                os.system(command)
                os.chdir(cwd)
            else:
                os.chdir(outdir)
                os.system(command)
                os.chdir(cwd)
        else:
            logging.info("No resim flag was raised, so not resimulating.")

def print_help():
    #print help menu for config parameters
    print("\nConfig File Set Up: \n ")
    sys.stdout.flush()
    print("CID_MJD_LIST: ")
    sys.stdout.flush()
    print("  - CID1 MJD1 peak ... MJDn   #can also include 'peak' to substitute true peak mjd")
    sys.stdout.flush()
    print("  - CID2 MJD1 MJD2 ... MJDn \n    . \n    . \n    .")
    sys.stdout.flush()
    print("  - CIDn MJD1 MJD2 ... MJDn")
    sys.stdout.flush()
    print("SIMFOLDER: /path/to/simfolder   #add a wildcard to use multiple folders with the same prefix") 
    sys.stdout.flush()
    print("KCOR_PATH: /path/to/kcor")
    sys.stdout.flush()
    print("SIMGEN_EAZY: simgen_eazy_templates_file.INPUT")
    sys.stdout.flush()
    print("HOSTLIB_FILE: /path/to/hostlib")
    sys.stdout.flush()
    print("OUTDIR_SIM_INPUT: /path/to/store_resim_inputs #optional OUTDIR for new inputs, otherwise create 'store inputs' folder in current dir" )
    sys.stdout.flush()
    print("")
    sys.stdout.flush()

def print_verbose_main(config_dic, cids_mjd_rows):
    #print all verbose statements in main method
    logging.info('YAML contained: %s', config_dic)
    logging.info('CIDS + mjds are: %s', cids_mjd_rows)

# =============================================
#       MAIN
# =============================================
if __name__ == "__main__":
    logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    no_resim, verbose, perfect, config_dic, Help, extra_sim_args = get_args()
    cid_mjd_rows = config_dic['CID_MJD_LIST'] 

    if verbose==True:
        print_verbose_main(config_dic, cid_mjd_rows)
    
    for row in cid_mjd_rows:
        cid = row.split()[0]
        mjds = row.split()[1:] 
        process_cid(cid, mjds, config_dic) 
