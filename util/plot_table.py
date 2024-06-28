#!/usr/bin/env python
#
# Created by B.Popovic during his graduate career at Duke University.
# Installed into SNANA Jun 24 2024 by R.Kessler
# Refactor to have __main__, and add translate_VARIABLE and tranlate_CUT
# methods to automatically append data commands to simplified user input.
# 
# TO DO:
#   - translate_VARIABLE needs to catch and fix np.math(var+var2...)
#      and has to work with combinaton of var and number
#
#   + @@HELP to call print_help(); shorten strings in python help
#   + @@TITLE
#   - write underflow/overflow stats
#   - mask of stat info to display on plot
#
# ==============================================
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

BOUNDS_AUTO    = "AUTO"  # indicates default axis bounds for plot


DELIMITER_VAR_LIST  = [ '+', '-', '/', '*', ':' ]  # for @@VARIABLE 
DELIMITER_CUT_LIST  = [ '&', '|', '>', '<', '=' ]  # for @@CUT

# define pandas strings to add into VARIABLE and CUT strings.
STR_DF         = 'df.'
STR_DF_LOC     = 'df.loc'

# internal flag to exit after translating VARIALBE and CUT
#DEBUG_TRANSLATE = True   
DEBUG_TRANSLATE = False


# ================================
# HELP

def print_help():

    help_string = \
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

Input table files are loaded with pandas dataframes, and therefore input variables
are internally converted to pandas notation. For instance, to plot redshift 
distribution require user input
    @@VARIALBE zHD
      or
    @@VARIALBE df.zHD

2D plots are parsed by a colon (:) to separate the x and y variables. 
The syntax is generally x:y For instance, 
   @@VARIABLE zHD:x1      # displays x1 vs z. 
      or
   @@VARIABLE df.zHD:df.x1      # displays x1 vs z. 

In short, this code internally prepends pandas notation (df.) to simplify
command-line inputs, or you can input the pandas notation with the df symobls. 
The @@VARIALBE must have all df or none of them; cannot mix variables with and
without df. This script writes out the pandas-translated @@VARIABLE string 
that can be scooped up for other plotting purposes.

If a command without df fails, please post github issue.

Algebraic and numpy functions are also allowed, e.g., 
    @@VARIABLE zHD:mB - 3.1*c + 0.16*x1
        or
    @@VARIABLE df.zHD:df.mB - 3.1*df.c + 0.16*df.x1

    @@VARIABLE np.sqrt(MUERR**2 + .05**2):zHD
         or
    @@VARIABLE np.sqrt(df.MUERR**2 + .05**2):df.zHD


Custom axis boundaries are input with 
   @@BOUNDS xmin xmax xbin                   # 1D
   @@BOUNDS xmin xmax xbin: ymin ymax ybin   # 2D
Mean, Median, stdev only include entries within the plot bounds; 
overflows are ignored.


Use  @@SAVE to save figure as pdf or png or jpeg (based on user-suppled extension); 
e.g.
    @@SAVE my_first_table_plot.pdf
       or
    @@SAVE my_first_table_plot.png

To compare and contrast values for the same CID, enable the @@DIFF option! 
There are two valid DIFF options:
    @@DIFF ALL    #  compare the difference in median values. 
      or
    @@DIFF CID    #  compare the difference for the same CIDs. 

@@CUT applies selection cuts on the sample. As with @@VARIABLE, CUT is internally 
translated to append df and df.loc as needed, or the user can imput the pandas 
notation; e.g,
   @@CUT "zHD > 0.2 &  zHDERR<.1 & SNRMAX1 < 100 & FIELD='C3'"
      or
   @@CUT "df.loc[(df.zHD > 0.2) &  (df.zHDERR<.1) & (df.SNRMAX1 < 100) & (df.FIELD=='C3')]"

Finally, @@ALPHA adjusts the matplot alpha values to adjust transparency 
(0=transparent, 1=solid)


Examples:

plot_table.py @@TFILE File1.FITRES File2.FITRES \\
   @@VARIABLE zHD:mB - 3.1*c + 0.16*x1 - MU \\
   @@CUT "IDSURVEY < 15 & zHD>0.1" 

plot_table.py @@TFILE scone_predict_diff.text \\
   @@VARIABLE PROB_SCONE_2:PROB_SCONE_1 \\
   @@CUT "PROB_SCONE_2 > 0"


"""

    sys.exit(f"\n{help_string}\n")
    return

def setup_logging():

    #logging.basicConfig(level=logging.DEBUG,
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)8s | %(message)s")
    # xxx mark format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)

    return


def get_args():
    parser.add_argument('@@TFILE', help="Name of TABLE file or space delineated list of different TABLE files. If >1 TFILE, plots are overlaid", nargs="+")
    
    parser.add_argument('@@VARIABLE', help="""Variable(s) to plot from table file, or function of variables. \n
One variable results in 1D histogram. Two colon delimited variables, such as zHD:mB, plots mB (y) vs zHD (x). For the histogram, counts are normalised to the first table file.""", nargs="+")
    
    parser.add_argument('@@BOUNDS', default=BOUNDS_AUTO, help=
"""AUTO (default): optional bounds are maximum and minimum values of specified parameter. \n 
An example is min:max:binsize . \n
For 2D plot, must specify both x and y bounds. The binsize for y is ignored. """, nargs = '+')

    parser.add_argument('@@TITLE', default=None, help="Override default plot title = arg of @@CUTS")
    
    parser.add_argument('@@SAVE', default='None', help=
"""Filename to save image under. Can give a custom filepath, otherwise saves in the working directory. Default does not save images.""")
    
    parser.add_argument('@@DIFF', default=None, type=str, help=
"""Plot the difference in the y-axis between files. Valid options are None, ALL, and CID. \n
ALL will plot the median difference between the first file and subsequent ones. \n
CID will plot the per-CID difference.
""")
    
    parser.add_argument('@@ALPHA', default=0.3, type=float, help='Alpha value for plotting. Set to 0 if you just want to see averages. If you set ALPHA = 0 and DIFF = True, you can compare the average difference between the two files even if there are no overlapping CIDS.')

    parser.add_argument("@@CUT", help="cuts with boolean and algegraic operations; see @@HELP:", nargs="+")

    parser.add_argument("@@NROWS", help="number of rows to read (for large files)", type=int, default=0)

    parser.add_argument("@@HELP", help="Full help menu printed to stdout", action="store_true")

    args    = parser.parse_args()
    
    if args.HELP:
        print_help()
    
    # - - - - - - -
    # fix inputs under the hood; mianly remove pad spacing
    if args.DIFF:
        args.DIFF = args.DIFF.strip()  # remove pad spacing

    args.VARIABLE_ORIG = args.VARIABLE    
    if args.VARIABLE:
        args.VARIABLE = ''.join([str(elem) for elem in args.VARIABLE])
        
    args.CUT_ORIG  = args.CUT
    if args.CUT:
        args.CUT = ''.join([str(elem) for elem in args.CUT])
        args.CUT_ORIG = args.CUT_ORIG[0]

    if args.BOUNDS != BOUNDS_AUTO:
        args.BOUNDS = ' '.join([str(elem) for elem in args.BOUNDS])

    # tack on new name space elements that are trivially dependent on user input

    table_list      = []
    table_base_list = []  # base names only; for plot legend
    for table_file in args.TFILE:
        table_list.append( os.path.expandvars(table_file) )
        table_base_list.append( os.path.basename(table_file) )
        
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

    var_list = sorted(VAR_TMP.split(' '))

    isnum_list = []
    for var in var_list:
        isnum_list.append(is_number(var))

    return  var_list, isnum_list

    # end get_var_list
    
def translate_VARIABLE(args):
    # add df. as needed to args.VARIABLE
    # assume that all df. have been provided by user, or none;
    # will not handle in-between cases.
    
    VARIABLE = args.VARIABLE
    
    if STR_DF in VARIABLE: return

    # break down VARIABLE into a list of variables without any symbols
    var_list, isnum_list  = get_var_list(VARIABLE, DELIMITER_VAR_LIST)

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

    # add df.loc, df. and () as needed to args.CUT
    # Add "(df." in front of each var_list element that is NOT a number.
    # Add ")" after each var_list eleement that is a number.

    CUT = args.CUT
    if not CUT: return
    if STR_DF in CUT: return

    # '=' is the only delimeter where user might use '==' instead,
    # and 2-char delimiter totally breaks the logic below. Rather 
    # than abort, just fix it here so that FIELD='C3' or FIELD=='C3' 
    # will both work.
    if '==' in CUT:
        CUT = CUT.replace('==', '=')

    # examine CUT string and return list of var names and numbers that
    # represent cut values.
    var_list, isnum_list  = get_var_list(CUT, DELIMITER_CUT_LIST)
    
    if DEBUG_TRANSLATE:
        print(f" xxx CUT var_list   = {var_list}")
        print(f" xxx CUT isnum_list = {isnum_list}")

    for var, isnum in zip(var_list, isnum_list):
        has_quotes = "'" in var  # e.g., FIELD='C3'

        if isnum or has_quotes:
            var_parenth = var + ')'
            if var_parenth not in CUT:
                CUT = CUT.replace(var,var_parenth)

        else:
            df_var = STR_DF + var        
            if df_var not in CUT:  # modify only if not already modified
                CUT = CUT.replace(var, '(' + df_var)

    # replace input for '=' with '=='
    if '=' in CUT:
        CUT = CUT.replace('=', '==')

    # finally, wrap entire cut in df.loc[ CUT ]
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
    plotdic            = {}  # has df for making plots
    plotdic_axis_label = {}  # cleaned up for axis labels

    if args.TITLE:
        plot_title    = args.TITLE
    elif args.CUT_ORIG :
        plot_title    = args.CUT_ORIG
    else:
        plot_title    = args.VARIABLE_ORIG

    for n, VAR in enumerate(VARIABLE.split(":")):
        STR_VAR   = str(VAR)
        STR_LABEL = STR_VAR.replace('df.','') 
        if n == 0:
            plotdic['x'] = STR_VAR
            plotdic_axis_label['x'] = STR_LABEL
        else:
            plotdic['y'] = STR_VAR
            plotdic_axis_label['y'] = STR_LABEL

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
    plot_info.plotdic             = plotdic
    plot_info.plotdic_axis_label  = plotdic_axis_label
    plot_info.plot_title          = plot_title
    plot_info.boundsdic           = boundsdic
    plot_info.custom_bounds       = custom_bounds
    
    return
    # end set_var_dict

def read_tables(args, plot_info):

    table_list      = args.table_list
    table_base_list = args.table_base_list
    CUT             = args.CUT
    NROWS           = args.NROWS

    plotdic    = plot_info.plotdic
    boundsdic  = plot_info.boundsdic
    
    MASTER_DF_DICT = {}  # dictionary of variables to plot (was MASTERLIST)
    
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

        # increment MASTER_DF_DICT dictionary; note that filename is dict key.

        MASTER_DF_DICT[l] = df
        nrow = len(df)
        name_legend = get_name_legend(l, table_list)
        logging.info(f"\t --> name_legend = {name_legend}")

        MASTER_DF_DICT[l] = {
            'df'           : df,
            'name_legend'  : name_legend
        }
        # xxx mark MASTER_DF_DICT[l]['name_legend'] = eval(name_legend)


        try:
            MASTER_DF_DICT[l]['df']['x_plot_val'] = eval(plotdic['x'])
            boundsdic[l + "_min"] = np.amin(MASTER_DF_DICT[l]['df']['x_plot_val'])
            boundsdic[l + "_max"] = np.amax(MASTER_DF_DICT[l]['df']['x_plot_val'])
            if len(plotdic) == 2:
                MASTER_DF_DICT[l]['df']['y_plot_val'] = eval(plotdic['y'])
        except AttributeError:
            sys.exit(f"\n ERROR: Couldn't set bounds for {plotdic} and {l}")

            
    # - - - - - - - -
    logging.info("Done loading all table files.")

    # load output namespace
    plot_info.MASTER_DF_DICT = MASTER_DF_DICT
    
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
    CUT_ORIG  = args.CUT_ORIG
    ALPHA     = args.ALPHA

    MASTER_DF_DICT       = plot_info.MASTER_DF_DICT
    plotdic              = plot_info.plotdic
    plotdic_axis_label   = plot_info.plotdic_axis_label
    boundsdic            = plot_info.boundsdic
    custom_bounds        = plot_info.custom_bounds
    plot_title           = plot_info.plot_title

    if custom_bounds:                       
        bins = np.arange(boundsdic['x'][0],boundsdic['x'][1],boundsdic['x'][2]) 
        plt.xlim([boundsdic['x'][0], boundsdic['x'][1]])  
    else:                                   
        bins = np.linspace(boundsdic[min(boundsdic, key=boundsdic.get)],
                           boundsdic[max(boundsdic, key=boundsdic.get)], 30)

        
    if len(plotdic) == 1:                           
        if (DIFF == 'ALL') or (DIFF == 'CID'):
            sys.exit("\nERROR: DIFF does not work for histograms. ABORT to avoid confusion.")

        for n, key_name in enumerate(MASTER_DF_DICT): 
            df = MASTER_DF_DICT[key_name]['df']   # recall that key_name is file name
            msg = f"The upper and lower bounds are: " \
                f"{np.around(bins[0],4)}  {np.around(bins[-1],4)} respectively"
            logging.info(msg)
            sb = binned_statistic(df.x_plot_val, df.x_plot_val, bins=bins, statistic='count')[0] #Get counts            
            errl,erru = poisson_interval(sb) # And error for those counts
            nevt = np.sum(sb)                # number of events before normalization
            
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
                
        plt.xlabel(plotdic_axis_label['x'])
        plt.legend()                     
        plt.title(plot_title)
            
    elif (DIFF == 'CID') or (DIFF == 'ALL'):
        # 2D plot of difference between two files
        logging.info("plotting DIFF between two files.")
        try:         
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]])            
        except KeyError:     
            pass  # auto scale y axis
            
        keylist    = list(MASTER_DF_DICT.keys()) # really, it's a file list
        df_ref = MASTER_DF_DICT[keylist[0]]['df']  # reference df  for difference
        for k in keylist[1:]:
            df = MASTER_DF_DICT[k]['df']
            if DIFF == 'CID':
                #need to do an inner join with each entry in dic, then plot the diff
                # (join logic thanks to Charlie Prior)
                join = df_ref.join(df.set_index('CID'), on='CID', how='inner', lsuffix='_1', rsuffix='_2')
                plt.scatter(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values,
                            alpha=ALPHA, label='Diff')
                avgdiff = binned_statistic(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values, bins=bins, statistic='median')[0]                                     
                plt.scatter((bins[1:] + bins[:-1])/2, avgdiff, label="Mean Difference", color='k')
            elif (DIFF == 'ALL'):
                try:
                    plt.scatter(df_ref.x_plot_val, df_ref.y_plot_val - df.y_plot_val, 
                                label=df_ref.name.values[0]+" - "+df.name.values[0], alpha=ALPHA)
                except ValueError:
                    pass
                median_ref = binned_statistic(df_ref.x_plot_val, df_ref.y_plot_val, bins=bins, statistic='median')[0]
                median = binned_statistic(df.x_plot_val, df.y_plot_val, bins=bins, statistic='median')[0]
                plt.scatter((bins[1:] + bins[:-1])/2., median_ref - median, 
                            label=df_ref.name.values[0]+" - "+k.name.values[0]+" median", marker="^", zorder=10)
            else:  
                sys.exit(f"\n ERROR: {str(DIFF)} is not a valid DIFF option -> ABORT.")         

            plt.xlabel(plotdic_axis_label['x'])
            plt.ylabel(plotdic_axis_label['y'] + " diff")                    
            plt.legend()                                
            plt.title(plot_title)       
    else:
        # 2D for each file
        try:
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]]) 
        except KeyError:
            pass  # auto scale axis
        
        for key_name, df_dict in MASTER_DF_DICT.items():
            df          = df_dict['df']
            name_legend = df_dict['name_legend']
            plt.scatter(df.x_plot_val, df.y_plot_val, alpha=ALPHA, label=name_legend, zorder=0)

            # overlay information on plot
            median = binned_statistic(df.x_plot_val, df.y_plot_val, bins=bins, statistic='median')[0] 
            plt.scatter((bins[1:] + bins[:-1])/2., median, label= name_legend+" median",
                        marker="^", zorder=10)
            
        plt.xlabel(plotdic_axis_label['x'])                    
        plt.ylabel(plotdic_axis_label['y'])                    
        plt.legend()                                
        plt.title(plot_title) 
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

    if  DEBUG_TRANSLATE :
        sys.exit(f"\n xxx bye.")
    
    plot_info = Namespace()  # someplace to store internally computed info
    
    set_var_dict(args, plot_info) # set plot bounds and axisinfo

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
    

    
