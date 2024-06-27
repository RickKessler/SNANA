#!/usr/bin/env python
#
# Created by B.Popovic during his graduate career at Duke University.
# Installed into SNANA Jun 24 2024 by R.Kessler
# Refactor to have main, and start translate_VARIABLE and tranlate_CUT
# methods to automatically append data frame df.
# 
import os, sys, gzip, copy, logging, math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import argparse
from argparse import RawTextHelpFormatter
from argparse import Namespace

parser=argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, prefix_chars='@')
from scipy.stats import binned_statistic
import distutils.util 

# ====================
# globals

BOUNDS_AUTO    = "AUTO"
DELIMITER_VAR_LIST  = [ '+', '-', '/', '*', ':' ]
DELIMITER_CUT_LIST  = [ '&', '|', '>', '<', '=' ]

STR_DF         = 'df.'
STR_DF_LOC     = 'df.loc'

"""
This plot unility works on 
  * FITRES table files create by light curve fitting code (snlc_fit.exe) and BBC (SALT2mu.exe)
  * M0DIF files from BBC
  * HOSTLIB files used in simulation
  * any file with same format that has VARNAMES key

BEWARE that conventional dashes in input keys are replaced with @@
E.g., --VARIABLE in any other python code is @@VARIABLE here ...
this is because dashes are confused with minus signs when plotting differences.

Next, specify what variables to plot ussing @@VARIABLE for  
  * 1D histogram or
  * 2D scatter plot or 
  * 1D or 2D function of variables. 

Input table files are loaded with pandas dataframes, and therefore you need to be familiar 
with that syntax. For instance, to plot redshift distribution require user input
    @@VARIALBE df.zHD
In short, need to use df pandas notation. For future, we are working on interpreter
to automatically append df internally, but unclear how long this will take.

2D plots are parsed by a colon (:) to separate the x and y variables. The syntax is generally x:y
For instance, @@VARIABLE df.zHD:df.x1  displays x1 vs z. 

You can also do algebra, e.g,
    @@VARIABLE df.zHD:df.mB - 3.1*df.c + 0.16*df.x1
will plot the distance for Tripp estimator. Algebra works in 1D and 2D.

Custom axis boundaries are input with 
   @@BOUNDS minimum maximum binsize          # 1D
   @@BOUNDS xmin xmax xbin: ymin ymax ybin   # 2D
Mean, Median, stdev only include entries in the plot; overflows are ignored.

Use  @@SAVE to save figure as pdf or png or jpeg (based on user-suppled extension)

To compare and contrast values for the same CID, enable the @@DIFF option! 
There are two valid DIFF options. 'ALL' will compare the difference in median values. 
'CID' will compare the difference for the same CIDs. 

@@CUT applies selection cuts on the sample. This is usually a 'df.loc[]' option of some capacity. 
This option MUST have quotation marks around it (due to bash issues)

Finally, @@ALPHA adjusts the matplot alpha values to adjust transparency (0=transparent, 1=solid)

Examples:

plot_table.py @@FITRES File1.FITRES File2.FITRES \
   @@VARIABLE df.zHD:df.mB - 3.1*df.c + 0.16*df.x1 - df.MU \
   @@CUT "df.loc[df.IDSURVEY < 15]"

plot_table.py @@FITRES scone_predict_diff.text \
   @@VARIABLE df.PROB_SCONE_2:df.PROB_SCONE_1 \
   @@CUT "df.loc[df['PROB_SCONE_2']>0]"

"""

def setup_logging():

    #logging.basicConfig(level=logging.DEBUG,
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)8s | %(message)s")
    # xxx mark format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)

    return


def get_args():
    parser.add_argument('@@FITRES', help='The location of your FITRES file that you need plotted. You can give a space delineated list of different FITRES files. Please make sure they are named differently, as their names (not full filepath) will be used for labeling.', nargs="+")
    
    parser.add_argument('@@VARIABLE', help="""The variable you want to plot. Needs to be a valid FITRES parameter. \n
If you give only one value, this will generate a histogram. If you give two colon delimited values, such as df.zHD:df.mB it will plot zHD (x) vs mB (y). Please note that for the histogram, counts will be normalised to the first file given.""", nargs="+")
    
    parser.add_argument('@@BOUNDS', default=BOUNDS_AUTO, help=
"""AUTO (default): bounds are maximum and minimum values of specified parameter. \n 
Custom: Give a set of numbers, colon delimited. An example is min:max:binsize \n
Please note, if you are doing a two dimensional plot, you need to specify both x and y sets in order. The binsize for y will not be used. """, nargs = '+')
    
    parser.add_argument('@@SAVE', default='None', help=
"""Filename to save image under. Can give a custom filepath, otherwise saves in the working directory. Default does not save images.""")
    
    parser.add_argument('@@DIFF', default=None, type=str, help=
"""Plot the difference in the y-axis between files. Valid options are None, ALL, and CID. \n
ALL will plot the median difference between the first file and subsequent ones. \n
CID will plot the per-CID difference.
""")
    
    parser.add_argument('@@ALPHA', default=0.3, type=float, help='Alpha value for plotting. Set to 0 if you just want to see averages. If you set ALPHA = 0 and DIFF = True, you can compare the average difference between the two files even if there are no overlapping CIDS.')

    parser.add_argument("@@CUT", help="NEEDS TO BE GIVEN IN QUOTATION MARKS!!! SUPER IMPORTANT!!! This takes the form of a df.loc[] option, typically. Any sort of cuts you want to make.", nargs="+")

    parser.add_argument("@@NROWS", help="choose number of rows to read in for larger files", type=int, default=0)

    args    = parser.parse_args()
        
    # - - - - - - -
    # fix inputs under the hood; mianly remove pad spacing
    if args.DIFF:
        args.DIFF = args.DIFF.strip()  # remove pad spacing

    if args.VARIABLE:
        args.VARIABLE = ''.join([str(elem) for elem in args.VARIABLE])
        args.VARIABLE_ORIG = args.VARIABLE
    
    if args.CUT:
        args.CUT = ''.join([str(elem) for elem in args.CUT])
        args.CUT_ORIG      = args.CUT
        
    if args.BOUNDS != BOUNDS_AUTO:
        args.BOUNDS = ' '.join([str(elem) for elem in args.BOUNDS])

    # tack on new name space elements that are trivially dependent on user input

    table_list      = []
    table_base_list = []  # base names only; for plot legend
    for table_file in args.FITRES:
        table_list.append( os.path.expandvars(table_file) )
        table_base_list.append( os.path.basename(table_file) )
        
    # xxx fitres_list = [l.split("/")[-1] for l in FILENAME]
    if (any(table_list.count(x) > 1 for x in table_list)):
        sys.exit(f"\n ERROR: found duplicate table file; check {table_list}")

    args.table_list      = table_list  # ENVs are expanded
    args.table_base_list = table_base_list
    
    return args

    # end get_args()

def get_var_list(VARIABLE, DELIMITER_LIST):
    # if VARIABLE = 'zHD-zHD_2:SNRMAX' -> return var_list = ['zHD', 'zHD_2', 'SNRMAX']
    # which is a list of variables wihtout symbols

    # first replace any valid delimiter with '!' so that we can split on
    # single ! char
    VAR_TMP = copy.copy(VARIABLE)
    VAR_TMP = VAR_TMP.replace(' ','')

    for delim in DELIMITER_LIST:
        VAR_TMP = VAR_TMP.replace(delim,' ')

    var_list = VAR_TMP.split(' ')
    return  sorted(var_list)
    # end get_var_list
    
def translate_VARIABLE(args):
    # add df. as needed to args.VARIABLE
    # assume that all df. have been provided by user, or none;
    # will not handle in-between cases.
    
    VARIABLE = args.VARIABLE
    
    if STR_DF in VARIABLE: return

    # break down VARIABLE into a list of variables without any symbols
    var_list  = get_var_list(VARIABLE, DELIMITER_VAR_LIST)

    # next, put df. in front of each variable ... but be careful about
    # variable names that are subsets of others. E.g., naively prepending
    # df in from of z and z_2 results in df.z and df.df.z_2.
    for var in var_list:
        df_var = STR_DF + var
        if df_var not in VARIABLE:  # modify only if not already modified
            VARIABLE = VARIABLE.replace(var,df_var)

    logging.info(f"Translate VARIABLE {args.VARIABLE_ORIG}  ->  {VARIABLE}")
    
    args.VARIABLE = VARIABLE
    return

def translate_CUT(args):
    # add df. as needed to args.CUT

    CUT = args.CUT
    if not CUT: return
    if STR_DF in CUT: return

    var_list  = get_var_list(CUT, DELIMITER_CUT_LIST)
    
    for var in var_list:
        var = var.replace('(','')
        var = var.replace(')','')                
        if is_number(var) : continue
        if "'"  in var    : continue  # e.g beware of FIELD='DEEP'
        df_var = STR_DF + var        
        if df_var not in CUT:  # modify only if not already modified
            CUT = CUT.replace(var,df_var)

    CUT = STR_DF_LOC + '[' + CUT + ']'
    
    args.CUT = CUT

    logging.info(f"Translate CUT {args.CUT_ORIG}  ->  {CUT}")    
    #print(f"\n xxx var_list = {var_list}")
    
    return

def is_number(string): 
    try: 
        float(string) 
        return True
    except ValueError: 
        return False
    return
        
def set_var_dict(args, plot_info):

    # strip off user args
    VARIABLE = args.VARIABLE
    BOUNDS   = args.BOUNDS

    # store x and [optional] y variable names
    plotdic = {}
    for n,VAR in enumerate(VARIABLE.split(":")):
        if n == 0:
            plotdic['x'] = str(VAR)
        else:
            plotdic['y'] = str(VAR)

    # check for user (custom) bounds in plot
    custom_bounds = False
    boundsdic     = {}
    if BOUNDS != BOUNDS_AUTO:
        for n,BND in enumerate(BOUNDS.split(':')):
            custom_bounds = True
            if n == 0:
                boundsdic['x'] = [float(i) for i in BND.split()]
            else:
                boundsdic['y'] = [float(i) for i in BND.split()]

    #print(f" xxx plotdic = {plotdic}")
    
    # load output namespace
    plot_info.plotdic       = plotdic
    plot_info.boundsdic     = boundsdic
    plot_info.custom_bounds = custom_bounds
    
    return
    # end set_var_dict

def read_tables(args, plot_info):

    table_list      = args.table_list
    table_base_list = args.table_base_list
    CUT             = args.CUT
    NROWS           = args.NROWS

    plotdic    = plot_info.plotdic
    boundsdic  = plot_info.boundsdic
    
    MASTER_VAR_DICT = {}  # dictionary of variables to plot (was MASTERLIST)
    
    for l in table_list:
        l_base = os.path.basename(l)
        logging.info(f"Loading {l}")
        if not os.path.exists(l):
            sys.exit(f"\n ERROR: cannot find {l}")
            
        df  = pd.read_csv(l, comment="#", sep=r"\s+")
        if NROWS > 0 :
            # read NROWS subset
            df  = pd.read_csv(l, comment="#", sep=r"\s+", nrows=NROWS)
        else:
            # read all
            df  = pd.read_csv(l, comment="#", sep=r"\s+")

        try:
            df['CID'] = df['CID'].astype(str)
        except KeyError:
            logging.warn("No CIDs present in this file. OK for some file types.")

        # apply user cuts
        if CUT: df = eval(CUT)

        # increment MASTER_VAR_DICT dictionary; note that filename is dict key.

        MASTER_VAR_DICT[l] = df
        name_legend = get_name_legend(l, table_list)
        logging.info(f"\t --> name_legend = {name_legend}")

        MASTER_VAR_DICT[l]['name_legend'] = name_legend 

        try:
            MASTER_VAR_DICT[l]['x_plot_val'] = eval(plotdic['x'])
            boundsdic[l + "_min"] = np.amin(MASTER_VAR_DICT[l]['x_plot_val'])
            boundsdic[l + "_max"] = np.amax(MASTER_VAR_DICT[l]['x_plot_val'])
            if len(plotdic) == 2:
                MASTER_VAR_DICT[l]['y_plot_val'] = eval(plotdic['y'])
        except AttributeError:
            sys.exit(f"\n ERROR: Couldn't set bounds for {plotdic} and {l}")

            
    # - - - - - - - -
    logging.info("Done loading all table files.")

    # load output namespace
    plot_info.MASTER_VAR_DICT = MASTER_VAR_DICT
    
    return
    # end read_tables

    
def get_name_legend(table_file, table_list):

    # for input table_file, return name to put in plot legend.
    # If there are duplicate base names in table_file_list, then
    # modify legend name accordinngly.

    base        =  os.path.basename(table_file)
    name_legend =  base  # default 

    if table_file == table_list[0] : return name_legend

    # - - - - -
    # alter name_legend for duplicate base names
    n_base_match = 0
    for t in table_list :
        b    = os.path.basename(t)
        if b == base:
            n_base_match += 1
            if t != table_file:
                name_legend += f"-{n_base_match}"
    
    return name_legend
    # end get_name_legend
    

def poisson_interval(k, alpha=0.32):
    """  
    uses chisquared info to get the poisson interval. Uses scipy.stats                   
    (imports in function).                                                               
    (http://stackoverflow.com/questions/14813530/poisson-confidence-interval-with-numpy) 
    """
    from scipy.stats import chi2
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    low[k == 0] = 0.0
    #if k == 0:                                                                          
    #    low = 0.0                                                                       
    return low, high

def plotter_func(args, plot_info):

    # utility to create the plot (but doesn't show it)

    # strip off local args from input name spaces
    DIFF      = args.DIFF
    CUT       = args.CUT
    ALPHA     = args.ALPHA

    MASTER_VAR_DICT  = plot_info.MASTER_VAR_DICT
    plotdic          = plot_info.plotdic
    boundsdic        = plot_info.boundsdic
    custom_bounds    = plot_info.custom_bounds
    
    if custom_bounds:                       
        bins = np.arange(boundsdic['x'][0],boundsdic['x'][1],boundsdic['x'][2]) 
        plt.xlim([boundsdic['x'][0], boundsdic['x'][1]])  
    else:                                   
        bins = np.linspace(boundsdic[min(boundsdic, key=boundsdic.get)],
                           boundsdic[max(boundsdic, key=boundsdic.get)], 30)

        
    if len(plotdic) == 1:                           
        if (DIFF == 'ALL') or (DIFF == 'CID'):
            sys.exit("The DIFF feature does not work for histograms. ABORT to avoid confusion.")

        for n, key_name in enumerate(MASTER_VAR_DICT): 
            k = MASTER_VAR_DICT[key_name]   # recall that key_name is file name
            msg = f"The upper and lower bounds are: " \
                f"{np.around(bins[0],4)}  {np.around(bins[-1],4)} respectively"
            logging.info(msg)
            sb = binned_statistic(k.x_plot_val, k.x_plot_val, bins=bins,
                                  statistic='count')[0] #Get counts            
            errl,erru = poisson_interval(sb) # And error for those counts
            nevt = np.sum(sb)  # number of events before normalization
            
            if n==0 :
                sb0 = copy.deepcopy(sb)  # preserve 1st file contents to normalize other files
            else:
                sb *= np.sum(sb0) / np.sum(sb) # normalize integral to match file 0

            name_legend = k['name_legend'][0]
            plt.errorbar((bins[1:] + bins[:-1])/2., sb, label=name_legend,
                         yerr=[sb-errl, erru-sb], fmt='o')                 

            stat_dict = {
                'nevt'   : nevt,
                'mean'   : np.mean(k.x_plot_val),
                'median' : np.median(k.x_plot_val),
                'stdev'  : np.std(k.x_plot_val)
                # overflow/underflow ??
            }
            for str_stat, val_stat in stat_dict.items():
                logging.info(f" {str_stat:8} value for {name_legend}:  {val_stat:.3f}")
                
        plt.xlabel(plotdic['x'])  #In this case, VAR = [string], so we're going to strip the list.                     
        plt.legend()                     
        if CUT:
            plt.title(CUT)
            
    elif (DIFF == 'CID') or (DIFF == 'ALL'):
        # 2D plot of difference between two files
        logging.info("plotting DIFF between two files.")
        try:         
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]])            
        except KeyError:     
            pass  # auto scale y axis
            
        keylist    = list(MASTER_VAR_DICT.keys()) # really, it's a file list
        df_ref = MASTER_VAR_DICT[keylist[0]]  # reference df  for difference
        # .xyz df_ref = differator
        for k in keylist[1:]:
            k = MASTER_VAR_DICT[k]
            if DIFF == 'CID':
                #need to do an inner join with each entry in dic, then plot the diff
                # (join logic thanks to Charlie Prior)
                join = df_ref.join(k.set_index('CID'), on='CID', how='inner', lsuffix='_1', rsuffix='_2')
                plt.scatter(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values,
                            alpha=ALPHA, label='Diff')
                avgdiff = binned_statistic(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values, bins=bins, statistic='median')[0]                                     
                plt.scatter((bins[1:] + bins[:-1])/2, avgdiff, label="Mean Difference", color='k')
            elif (DIFF == 'ALL'):
                try:
                    plt.scatter(df_ref.x_plot_val, df_ref.y_plot_val - k.y_plot_val, label=df_ref.name.values[0]+" - "+k.name.values[0], alpha=ALPHA)
                except ValueError:
                    pass
                median_ref = binned_statistic(df_ref.x_plot_val, df_ref.y_plot_val, bins=bins, statistic='median')[0]
                median = binned_statistic(k.x_plot_val, k.y_plot_val, bins=bins, statistic='median')[0]
                plt.scatter((bins[1:] + bins[:-1])/2., median_ref - median, 
                            label=df_ref.name.values[0]+" - "+k.name.values[0]+" median", marker="^", zorder=10)
            else:  
                sys.exit(f"\n ERROR: {str(DIFF)} is not a valid DIFF option -> ABORT.")         

            plt.xlabel(plotdic['x'])                    
            plt.ylabel(plotdic['y'] + " diff")                    
            plt.legend()                                
            if CUT: plt.title(CUT)       
    else:
        # 2D for each file
        try:
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]]) 
        except KeyError:
            pass  # auto scale axis
        
        for key_name, k in MASTER_VAR_DICT.items():
            name_legend = k['name_legend'][0]
            plt.scatter(k.x_plot_val, k.y_plot_val, alpha=ALPHA, label=name_legend, zorder=0)

            # overlay information on plot
            median = binned_statistic(k.x_plot_val, k.y_plot_val, bins=bins, statistic='median')[0] 
            plt.scatter((bins[1:] + bins[:-1])/2., median, label= name_legend+" median",
                        marker="^", zorder=10)
            
        plt.xlabel(plotdic['x'])                    
        plt.ylabel(plotdic['y'])                    
        plt.legend()                                
        if CUT:
            plt.title(CUT) 
    return 
    # end plotter_func


# ===================================================
#   Add main, June 2024
# ========================================
if __name__ == "__main__":

    setup_logging()
    logging.info("# ========== BEGIN plot_table.py  ===============")

    args = get_args()

    translate_VARIABLE(args) # add df. as needed to args.VARIABLE
    translate_CUT(args)      # add df. and df.loc as needed to args.CUT

    #sys.exit(f"\n xxx bye.")
    
    plot_info = Namespace()  # someplace to store internally computed info
    
    set_var_dict(args, plot_info) # set plot bounds and axis info

    read_tables(args, plot_info)  # read each input file and store data frames

    plt.figure()  # initialize matplotlib figure
    
    plotter_func(args, plot_info)  # prepare plot in matplotlib

    # check for user-define output (e.g. myplot.pdf, myplot.png, etc ...)
    # Note that we either save to file, or show in pop-up window ... but not both.
    # This enables creating plots as part of pipelines.
    SAVE = args.SAVE
    if SAVE !="None":
        logging.info(f"Save plot to {SAVE}")
        plt.savefig(SAVE, bbox_inches="tight", format=SAVE.split(".")[-1])
    else:
        plt.show()  # show plot in pop-up window

    logging.info('Done.')
    
    # === END: ====
    

    
