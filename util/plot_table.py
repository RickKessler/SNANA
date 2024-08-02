#!/usr/bin/env python
#
# Created by B.Popovic during his graduate career at Duke University.
# Installed into SNANA Jun 24 2024 by R.Kessler
# Refactor to have __main__, and add translate_VARIABLE and tranlate_CUT
# methods to automatically append data commands to simplified user input.
# 
#
# ==============================================
import os, sys, gzip, copy, logging, math, re
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

DELIMITER_VAR_LIST  = [ '+', '-', '/', '*', ':', '(', ')' ]  # for @@VARIABLE 
DELIMITER_CUT_LIST  = [ '&', '|', '>', '<', '=', '*', '+', '-', '/' ]  # for @@CUT

COLON = ':'

# define pandas strings to add into VARIABLE and CUT strings.
STR_df         = 'df.'
STR_df_loc     = 'df.loc'
STR_np         = 'np.'

# define strings for @@OPT
OPT_NEVT      = "NEVT"       # append N={nevt} to legend
OPT_CHI2      = "CHI2"       # show tfile1/tfile2 chi2/dof and scale tfile2 to match tfile1
OPT_MEDIAN    = "MEDIAN"
OPT_DIAG_LINE = "DIAG_LINE"  # draw diagonal line on plot
OPT_LOGY      = "LOGY"       # log scale along Y axis
OPT_GRID      = "GRID"
OPT_LIST_CID  = "LIST_CID"   # list CIDs passing cuts

VALID_OPT_LIST = [ OPT_NEVT, OPT_CHI2, OPT_MEDIAN, OPT_DIAG_LINE, OPT_LOGY,
                   OPT_GRID, OPT_LIST_CID ]

NMAX_CID_LIST = 20  # max number of CIDs to print for @@OPT CID_LIST

# internal flag to exit after translating VARIALBE and CUT
#DEBUG_TRANSLATE = True   
DEBUG_TRANSLATE = False


# ================================

def print_help():

    help_string = \
f"""
This plot unility works on 
  * FITRES table files create by light curve fitting code (snlc_fit.exe) and BBC (SALT2mu.exe)
  * M0DIF files from BBC
  * HOSTLIB files used in simulation
  * any file with same format that has VARNAMES key

BEWARE that conventional dashes in input keys are replaced with @@
E.g., --VARIABLE in any other python code is @@VARIABLE here ...
this is because dashes are confused with minus signs when plotting differences.

       IF A COMMAND FAILS, PLEAE POST SNANA-GITHUB ISSUE !!! 

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

Algebraic and numpy functions are also allowed, e.g., 
    @@VARIABLE zHD:mB - 3.1*c + 0.16*x1
        or
    @@VARIABLE df.zHD:df.mB - 3.1*df.c + 0.16*df.x1

    @@VARIABLE 'np.sqrt(MUERR**2 + .05**2):zHD'
         or
    @@VARIABLE 'np.sqrt(df.MUERR**2 + .05**2):df.zHD'

Single quotes are numpy expressions are needed to avoid linux problems parsing ().

Custom axis boundaries are input with 
   @@BOUNDS xmin xmax xbin                   # 1D
   @@BOUNDS xmin xmax xbin: ymin ymax ybin   # 2D
Mean, Median, stdev only include entries within the plot bounds; 
overflows are ignored.

Units are defined with
   @@UNITS days       # units for 1D plot or x-axis of 2D plot
   @@UNITS :deg       # units for y-axis of 2D plot
   @@UNITS days:deg   # units for both axes of 2D plot
     or
   @U <arg>

Use @@SAVE to save figure as pdf or png; e.g.
    @@SAVE my_first_table_plot.pdf
       or
    @@SAVE my_first_table_plot.png
Plot is either displayed in pop-up window or stored in @@SAVE file ...
but not both. This enables pipelines to make post-processing plots
without worrying about pop-windows in slurm jobs.

To compare and contrast values for the same CID, enable the @@DIFF option! 
There are two valid DIFF options:
    @@DIFF ALL    #  compare the difference in median values. 
      or
    @@DIFF CID    #  compare the difference for the same CIDs. 

@@CUT applies selection cuts on the sample. As with @@VARIABLE, CUT is internally 
translated to append df and df.loc as needed. Beware that quotes (single or double)
are required around all cuts: e.g,

   @@CUT "zHD > 0.2 &  zHDERR<.1 & SNRMAX1 < 100 & FIELD='C3'"

If a cut includes a string matches (e.g., FIELD='C3'), single quotes must be 
used around the string to match, and double quotes around the entire CUT arg.
To overlay the same variable with different cuts, provide a list of cuts,
   @@CUT  'SNRMAX1>10'  'SNRMAX1>20'  'SNRMAX1>40'
results in 3 overlaid plots, and default legend shows each cut.

@@ALPHA adjusts the matplot alpha values to adjust transparency 
(0=transparent, 1=solid).
For 2D overlays of multuple files or cuts, can specify multiple 
alpha values, e.g.,
   @@ALPHA 0.9 0.1
so that primary plot is dark and overlay plot is nearly transparent.


  @@OPT  {' '.join(VALID_OPT_LIST)}
    where
      {OPT_NEVT:<12} ==> append N=Nevt on each legend
      {OPT_CHI2:<12} ==> display chi2/dof on plot for two table files
      {OPT_MEDIAN:<12} ==> plot median (vertical axis) for 2D plot
      {OPT_DIAG_LINE:<12} ==> draw line with slope=1 for 2D plot
      {OPT_LOGY:<12} ==> log scale for vertical axis
      {OPT_GRID:<12} ==> display grid on plot 
      {OPT_LIST_CID:<12} ==> print up to 100 CIDs passing cuts (1D plot only)

Examples:

plot_table.py @@TFILE File1.FITRES File2.FITRES \\
   @@VARIABLE zHD:mB - 3.1*c + 0.16*x1 - MU \\
   @@CUT "IDSURVEY < 15 & zHD>0.1" 

plot_table.py @@TFILE scone_predict_diff.text \\
   @@VARIABLE PROB_SCONE_2:PROB_SCONE_1 \\
   @@CUT "PROB_SCONE_2 > 0"


       IF A COMMAND FAILS, PLEAE POST SNANA-GITHUB ISSUE !!! 
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

    msg = "required: Name of TABLE file or space delineated list of " \
          "different TABLE files. If >1 TFILE, plots are overlaid"
    parser.add_argument('@@TFILE', help=msg, nargs="+")
    
    msg = "required: Variable(s) to plot from table file, or function of variables.\n" \
          "One variable results in 1D histogram. Two colon delimited variables, " \
          "such as zHD:mB, plots mB (y) vs zHD (x). " \
          "For the histogram, counts are normalised to the first table file."
    parser.add_argument('@V', '@@VARIABLE', help=msg, nargs="+")
    
    msg = "AUTO (default): optional bounds min, max and binsize of plot parameters.\n" \
          "For 2D plot, must specify both x and y bounds. y-binsize is ignored."
    parser.add_argument('@@BOUNDS', default=BOUNDS_AUTO, help=msg, nargs = '+')

    msg = "Units to show for each axis; see @@HELP"
    parser.add_argument('@U', '@@UNITS', default=None, help=msg, nargs="+")
    
    msg = "Override default plot title = arg of @@CUTS"
    parser.add_argument('@@TITLE',  default=None, help=msg )

    msg = "Override default legend on plot (space sep list per TFILE)"
    parser.add_argument('@@LEGEND', default=None, help=msg, nargs="+")

    msg = "Alpha value for plot. Set to 0 to see only averages." \
          "ALPHA=0 and DIFF=True compares average difference between two files, " \
          "even if there are no overlapping CIDS."
    parser.add_argument('@@ALPHA', default=[0.3], type=float, help=msg, nargs="+" )
    
    msg = "Filename to save plot.  Default: does not save plots."
    parser.add_argument('@@SAVE', default='None', help=msg )

    msg = "Plot difference in y-axis between files. Valid options are None, ALL, and CID.\n" \
          "ALL plots median difference between first file and subsequent ones.\n" \
          "CID plots per-CID difference."
    parser.add_argument('@@DIFF', default=None, type=str, help=msg )
    
    msg = "cuts with boolean and algegraic operations; see @@HELP"
    parser.add_argument("@@CUT", help=msg, nargs="+", default =[None])

    msg = "number of rows to read (for large files)"
    parser.add_argument("@@NROWS", help=msg, type=int, default=0)

    msg = "options; see @@HELP"
    parser.add_argument("@@OPT", help=msg, nargs="+", default = [])

    msg = "Full help menu printed to stdout"
    parser.add_argument("@H", "@@HELP", help=msg, action="store_true")

    args    = parser.parse_args()
    
    if args.HELP:
        print_help()    

    return args

    # end get_args()

def process_args(args):

    # process arguments such as remove pad spacing and forcing
    # number or @@TFILE args to match number of @@CUT args
    
    # - - - - - - -
        
    args.VARIABLE_ORIG = args.VARIABLE    
    if args.VARIABLE:
        args.VARIABLE      = ''.join([str(elem) for elem in args.VARIABLE])        
        args.VARIABLE_ORIG = args.VARIABLE_ORIG[0]

    # CUT is tricky. Make sure that length of cut list matchs length
    # of table-file list ... or vice versa ... make sure that length of
    # tfile list matches length of CUT list.
    if DEBUG_TRANSLATE:
        print(f" xxx proc_args: args.CUT = {args.CUT}  (before modifications)")

    n_tfile_orig = len(args.TFILE)
    n_cut_orig   = len(args.CUT)
    tfile_list   = copy.copy(args.TFILE)
    cut_list     = copy.copy(args.CUT)

    if n_cut_orig == 1 and n_tfile_orig > 1:
        cut_list = [ cut_list[0] ] * n_tfile_orig # same cut per tfile
        
    if n_cut_orig > 1 and n_tfile_orig == 1:
        tfile_list = [ tfile_list[0] ] * n_cut_orig # same tfile per cut

    if DEBUG_TRANSLATE:
        print(f" xxx proc_args: args.CUT   -> {args.CUT}")
        print(f" xxx proc_args: args.TFILE -> {args.TFILE}")
    
    table_list      = []
    table_base_list = []  # base names only; for plot legend
    for tfile in tfile_list :
        table_list.append( os.path.expandvars(tfile) )
        table_base_list.append( os.path.basename(tfile) )

    # xxx mark delete xxxx
    #if (any(table_list.count(x) > 1 for x in table_list)):
    #    sys.exit(f"\n ERROR: found duplicate table file; check {table_list}")
    # xxxxxxxxx
    
    args.tfile_list      = table_list           # ENVs are expanded
    args.tfile_base_list = table_base_list
    args.cut_list        = cut_list
    
    # - - - - -

    # make sure there is a colon in UNITS atg
    # e.g  "degrees" -> "degrees:"
    n_var = len(args.VARIABLE.split(COLON))
    if args.UNITS:
        args.UNITS  = ''.join([str(elem) for elem in args.UNITS])
        if COLON not in args.UNITS:
            args.UNITS += COLON
    else:
        args.UNITS = COLON    

        
    if args.BOUNDS != BOUNDS_AUTO:
        args.BOUNDS = ' '.join([str(elem) for elem in args.BOUNDS])

    # tack on new name space elements that are trivially dependent on user input

    
    if not args.LEGEND:
        args.LEGEND = [ None ] * len(args.TFILE)
        if len(args.CUT) > 1 :
            args.LEGEND = args.CUT
        else :
            args.LEGEND = args.tfile_base_list
        
    narg_legend       = len(args.LEGEND)
    args.legend_list  = args.LEGEND
    
    narg_tfile  = len(args.tfile_list)
    if narg_tfile != narg_legend:
        sys.exit(f"ERROR: narg_tfile={narg_tfile} but narg_legend={narg_legend}; " \
                 f"narg_legend must match number of table files.")

    # if only 1 alpha, make sure there is alpha for each file/cut .xyz
    narg_alpha = len(args.ALPHA)
    if narg_alpha < narg_tfile:
        a          = args.ALPHA[0]
        alpha_list = [a]*narg_tfile
    else:
        alpha_list = copy.copy(args.ALPHA)
    args.alpha_list = alpha_list
    
    if args.DIFF:
        args.DIFF = args.DIFF.replace(' ','')  # remove pad spacing

    # abort on invalid @@OPT
    invalid_OPT_list = []
    if args.OPT:
        for opt in args.OPT:
            if opt not in  VALID_OPT_LIST:
                invalid_OPT_list.append(opt)
        if len(invalid_OPT_list) > 0:
            sys.exit(f"\n ERROR: Invalid @@OPT {invalid_OPT_list}")
                
    return

def get_var_list(VARIABLE, DELIMITER_LIST):
    # if VARIABLE = 'zHD-zHD_2:SNRMAX' -> return var_list = ['zHD', 'zHD_2', 'SNRMAX']
    # which is a list of variables wihtout symbols

    # first replace any valid delimiter with '!' so that we can split on
    # single ! char
    VAR_TMP = copy.copy(VARIABLE)

    VAR_TMP = VAR_TMP.strip()  # remove pad spaces

    # replace algabraic delimiters with pad space for easier split
    for delim in DELIMITER_LIST:
        VAR_TMP = VAR_TMP.replace(delim,' ')

    var_list_all = sorted(VAR_TMP.split())

    # keep elements that do NOT have numpy "np."
    var_list = []
    for var in var_list_all:
        if STR_np not in var:
            var_list.append(var)

    # create supplemental list of logicals indicating which
    # var_list elements are numbers. A string in quotes is
    # treated like a number.
    isnum_list = []
    for var in var_list:
        isnum = is_number(var) or "'" in var
        isnum_list.append(isnum)

    return  var_list, isnum_list

    # end get_var_list
    
def translate_VARIABLE(VARIABLE):
    # add df. as needed to VARIABLE
    # assume that all df. have been provided by user, or none;
    # will not handle in-between cases.

    VARIABLE_ORIG = VARIABLE
    
    if STR_df in VARIABLE:
        return VARIABLE

    # break down VARIABLE into a list of variables without any symbols
    var_list, isnum_list  = get_var_list(VARIABLE, DELIMITER_VAR_LIST)

    if DEBUG_TRANSLATE:
        print(f" xxx VARIABLE var_list   = {var_list}")
        print(f" xxx VARIABLE isnum_list = {isnum_list}")

    # next, put df. in front of each variable ... but be careful about
    # variable names that are subsets of others. E.g., naively prepending
    # df in from of z and z_2 results in df.z and df.df.z_2.
    for var, isnum in zip(var_list,isnum_list):
        if isnum : continue
        df_var = STR_df + var
        if df_var not in VARIABLE:  # modify only if not already modified

            # .xyz PROBLEM c-SIM_c -> df.c-SIM_df.c unless makig only 1 replace;
            #     but then SIM_c-c  -> SIM_df.c-c fails ???
            #     need better logic when 1 variable name is substring of another.
            VARIABLE = VARIABLE.replace(var,df_var,1)

    logging.info(f"Translate VARIABLE {VARIABLE_ORIG}  ->  {VARIABLE}")
    
    return VARIABLE

def translate_CUT(CUT):

    # add df.loc, df. and () as needed to input cut
    # Add "(df." in front of each var_list element that is NOT a number.
    # Add ")" after each var_list eleement that is a number.
    #
    # July 24 2024 rewrite logic to split by bool and (),
    #   then wrap each item as '(df.' + item + ')'

    CUT_ORIG = CUT
    if not CUT:       return CUT
    if STR_df in CUT: return CUT

    CUT = CUT.replace(' ','')
    
    # '=' is the only delimeter where user might use '==' instead,
    # and 2-char delimiter totally breaks the logic below. Rather 
    # than abort, just fix it here so that FIELD='C3' or FIELD=='C3' 
    # will both work.
    if '==' in CUT:
        CUT = CUT.replace('==', '=')


    # split into sections separated by boolean &, |, or ()
    cut_list = re.split(r"\(|\)|\&|\|", CUT.replace(' ',''))
    cut_list = list(filter(None, cut_list))

    if DEBUG_TRANSLATE:
        print(f" xxx CUT      = {CUT}" )
        print(f" xxx cut_list = {cut_list}")

    # for each cut item, append parentheses and df.
    cut_list_df = []

    REFAC = True
    if REFAC:
        for cut in cut_list:
            # append df. in from of each string to allow cutting on
            # math functions of variables.
            var_list, isnum_list  = get_var_list(cut, DELIMITER_CUT_LIST)
            cut_df = cut
            for var, isnum in zip(var_list, isnum_list):
                if not isnum:
                    cut_df = cut_df.replace(var, STR_df + var)
            cut_df = '(' + cut_df + ')'
            cut_list_df.append(cut_df)
    else:
        # legac/obsolete
        for cut in cut_list:
            cut_df = '(' + STR_df + cut + ')'
            cut_list_df.append(cut_df)

    if DEBUG_TRANSLATE:
        print(f" xxx cut_list_df = {cut_list_df}")
        
    # replace each original cut with cut_df in CUT
    CUT_df = CUT
    for cut, cut_df in zip(cut_list, cut_list_df):
        CUT_df = CUT_df.replace(cut,cut_df)
        if DEBUG_TRANSLATE:
            print(f" xxx \t replace user cut {cut} --> {cut_df} ")

    # replace input for '=' with '=='
    if '=' in CUT_df:
        CUT_df = CUT_df.replace('=', '==')
            

    # finally, wrap entire cut in df.loc[ CUT ]
    CUT_df = STR_df_loc + '[' + CUT_df + ']'    

    logging.info(f"Translate CUT {CUT}  ->  {CUT_df}")
    #print(f"\n xxx var_list = {var_list}")
    
    return CUT_df


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
    UNITS    = args.UNITS
    
    # store x and [optional] y variable names
    plotdic            = {}  # has df for making plots
    plotdic_axis_label = {}  # cleaned up for axis labels

    if args.TITLE:
        plot_title    = str(args.TITLE)
    elif args.CUT :
        plot_title    = str(args.CUT)
    else:
        plot_title    = str(args.VARIABLE_ORIG)
    
    if STR_np in plot_title:
        plot_title = plot_title.replace(STR_np,'')

    # - - - - 
    VAR_LIST   = VARIABLE.split(COLON)
    UNITS_LIST = UNITS.split(COLON)
    
    for n, VAR in enumerate(VAR_LIST):
        STR_VAR   = str(VAR)
        STR_LABEL = STR_VAR.replace(STR_df,'') 
        if STR_np in STR_LABEL:  STR_LABEL = STR_LABEL.replace(STR_np,'')

        U = UNITS_LIST[n]
        if len(U) > 0:  STR_LABEL += '  (' + U + ')'
            
        if n == 0:
            plotdic['x']            = STR_VAR
            plotdic_axis_label['x'] = STR_LABEL
            plotdic_axis_label['y'] = None

        else:
            plotdic['y']            = STR_VAR
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

    tfile_list      = args.tfile_list
    tfile_base_list = args.tfile_base_list
    cut_list        = args.cut_list
    legend_list     = args.legend_list
    alpha_list      = args.alpha_list
    NROWS           = args.NROWS

    plotdic    = plot_info.plotdic
    boundsdic  = plot_info.boundsdic
    
    MASTER_DF_DICT = {}  # dictionary of variables to plot (was MASTERLIST)
    nf = 0
    
    for tfile, cut, legend, alpha in zip(tfile_list, cut_list, legend_list, alpha_list):
        tfile_base = os.path.basename(tfile)
        logging.info(f"Loading {tfile_base}")
        if not os.path.exists(tfile):
            sys.exit(f"\n ERROR: cannot find {tfile}")
            
        df  = pd.read_csv(tfile, comment="#", sep=r"\s+")

        if NROWS > 0 :
            # read NROWS subset
            df  = pd.read_csv(tfile, comment="#", sep=r"\s+", nrows=NROWS)
        else:
            # read all
            df  = pd.read_csv(tfile, comment="#", sep=r"\s+")

        try:
            df['CID'] = df['CID'].astype(str)
        except KeyError:
            logging.warn("No CIDs present in this file. OK for some file types.")

        # apply user cuts
        if cut:
            df = eval(cut)

        # increment MASTER_DF_DICT dictionary; note that filename index is dict key.
        key = f"tf{nf}"
        nf += 1
        
        MASTER_DF_DICT[key] = df
        nrow = len(df)
        name_legend = legend
            
        logging.info(f"\t --> nrow={nrow}   name_legend = {name_legend}")

        if nrow == 0:
            sys.exit(f"\n ERROR: zero rows read for {name_legend}")

        MASTER_DF_DICT[key] = {
            'df'           : df,
            'name_legend'  : name_legend,
            'alpha'        : alpha
        }

        try:
            MASTER_DF_DICT[key]['df']['x_plot_val'] = eval(plotdic['x'])
            boundsdic[key + "_min"] = np.amin(MASTER_DF_DICT[key]['df']['x_plot_val'])
            boundsdic[key + "_max"] = np.amax(MASTER_DF_DICT[key]['df']['x_plot_val'])
            if len(plotdic) == 2:
                MASTER_DF_DICT[key]['df']['y_plot_val'] = eval(plotdic['y'])
        except AttributeError:
            sys.exit(f"\n ERROR: Couldn't set bounds for {plotdic} and {tfile}")

            
    # - - - - - - - -
    logging.info("Done loading all table files.")
    logging.info("# - - - - - - - - - - ")
    
    # load output namespace
    plot_info.MASTER_DF_DICT = MASTER_DF_DICT
    
    return
    # end read_tables

    
def get_name_legend_default(table_file, table_list):

    # xxxxxx OBSOLETE xxxxxxx
    
    # for input table_file, return name to put in plot legend.
    # If there are duplicate base names in table_file_list, then
    # modify legend name accordinngly.

    # xxxxxx OBSOLETE xxxxxxx 
    base        =  os.path.basename(table_file)
    name_legend =  base  # default 

    if table_file == table_list[0] :
        return name_legend

    # - - - - -
    # alter name_legend for duplicate base names
    n_base_match = 0
    for t in table_list :
        b    = os.path.basename(t)
        if b == base:
            n_base_match += 1
            if t != table_file:
                name_legend += f"-{n_base_match}"
    # xxxxxx OBSOLETE xxxxxxx
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

    low[k == 0]  = 0.0
    high[k == 0] = 0.0    
    
    return low, high

def plotter_func(args, plot_info):

    # utility to create the plot (but doesn't show it)

    # strip off local args from input name spaces
    DIFF      = args.DIFF
    # xxx mark ALPHA     = args.ALPHA
    OPT       = args.OPT

    do_chi2      = OPT_CHI2      in OPT
    do_median    = OPT_MEDIAN    in OPT
    do_diag_line = OPT_DIAG_LINE in OPT
    do_list_cid  = OPT_LIST_CID  in OPT
    do_nevt      = OPT_NEVT      in OPT
    
    MASTER_DF_DICT       = plot_info.MASTER_DF_DICT
    plotdic              = plot_info.plotdic
    plotdic_axis_label   = plot_info.plotdic_axis_label
    xlabel               = plotdic_axis_label['x']
    ylabel               = plotdic_axis_label['y']
    boundsdic            = plot_info.boundsdic
    plotdic              = plot_info.plotdic
    custom_bounds        = plot_info.custom_bounds
    plot_title           = plot_info.plot_title

    
    if custom_bounds:                       
        xmin = boundsdic['x'][0]; xmax = boundsdic['x'][1]
        xbin = boundsdic['x'][2]
        bins = np.arange(xmin, xmax, xbin)
        plt.xlim(xmin, xmax)
        if len(plotdic) == 2:
            ymin = boundsdic['y'][0]; ymax = boundsdic['y'][1]
            ybin = boundsdic['y'][2]
            ybins = np.arange(ymin, ymax, ybin) # ignored 
            plt.ylim(ymin, ymax)
    else:                                   
        bins = np.linspace(boundsdic[min(boundsdic, key=boundsdic.get)],
                           boundsdic[max(boundsdic, key=boundsdic.get)], 30)
        xmin = bins[0];  xmax = bins[-1]


    msg = f"The min & max x-axis bounds are: " \
          f"{np.around(bins[0],4)}  {np.around(bins[-1],4)} respectively"
    logging.info(msg)

    
    if len(plotdic) == 1:                           
        if (DIFF == 'ALL') or (DIFF == 'CID'):
            sys.exit("\nERROR: DIFF does not work for histograms. ABORT to avoid confusion.")

        for n, key_name in enumerate(MASTER_DF_DICT): 
            df_dict     = MASTER_DF_DICT[key_name]
            df          = df_dict['df']
            name_legend = df_dict['name_legend']        
            plt_legend  = name_legend
            
            # get counts sb
            sb = binned_statistic(df.x_plot_val, df.x_plot_val, 
                                  bins=bins,
                                  statistic='count')[0]
            errl,erru = poisson_interval(sb) # And error for those counts
            nevt = np.sum(sb)                # nevt before normalization
            if do_nevt: plt_legend += f'  N={int(nevt)}'
            
            if n == 0 :
                sb0 = copy.deepcopy(sb)  # preserve 1st file contents to normalize other files
                logging.info(f"Plot {name_legend}")
                name0_legend = name_legend
            else:
                if do_chi2:
                    # re-scale overlay plot only if chi2 option is requested
                    scale = np.sum(sb0) / np.sum(sb)
                else:
                    scale = 1.0  # overlay plot is not scaled
                sb   *= scale # normalize integral to match file 0 
                errl *= scale
                erru *= scale
                logging.info(f"Overlay {name_legend} scaled by {scale:.3e}")

            # determine plot style: 
            # default is solid-filled circles with error bars
            do_errorbar = True 
            do_ovsim    = False  # overlay sim
            chi2red     = 0.0

            if n > 0 and do_chi2:
                # prepare for sim overlay with histogram
                do_errorbar = False; do_ovsim = True 
                
            if do_errorbar :
                plt.errorbar((bins[1:] + bins[:-1])/2., sb, label=plt_legend,
                             yerr=[sb-errl, erru-sb], fmt='o')  
            elif do_ovsim :
                x_val = df.x_plot_val
                wgt   = [ scale ] * len(x_val)
                plt.hist(x_val, bins, alpha=0.25, weights = wgt, 
                         label=plt_legend)
            else:
                sys.exit(f"\n ERROR: cannot determine which plot type: " \
                         f"errorbar or hist")

            stat_dict = {
                'nevt'    : nevt,
                'mean'    : np.mean(df.x_plot_val),
                'median'  : np.median(df.x_plot_val),
                'stdev'   : np.std(df.x_plot_val),
                # overflow/underflow ??
            }
            for str_stat, val_stat in stat_dict.items():
                logging.info(f"\t {str_stat:8} value for {name_legend}:  {val_stat:.3f}")
            if do_list_cid:
                print_cid_list(df, name_legend)
                
        setup_plot(args, plot_title, xlabel, ylabel) 

        # check option to compute and print chi2/dof info on plot
        # Froce min error =1 in chi2 calc so that it's ok to plot
        # error=0 for bins with zero events.
        if do_ovsim :
            sqdif = (sb0-sb)**2
            sqerr = np.maximum((sb0+sb*scale*scale), 1.0)
            chi2  = np.sum( sqdif / sqerr )
            ndof  = len(bins) - 1 
            text_chi2 = f"chi2/dof = {chi2:.1f}/{ndof}"
            logging.info(f"{name0_legend}/{name_legend} {text_chi2}") 
            x_text = xmin + 0.7*(xmax-xmin)
            y_text = 1.0 * np.max(sb0)  # warning; fragile coord calc
            plt.text(x_text, y_text, text_chi2 )

    elif (DIFF == 'CID') or (DIFF == 'ALL'):
        # 2D plot of difference between two files
        logging.info("plotting DIFF between two files.")
        try:         
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]])            
        except KeyError:     
            pass  # auto scale y axis
            
        keylist    = list(MASTER_DF_DICT.keys())   # tf0, tf1 ...
        df_ref = MASTER_DF_DICT[keylist[0]]['df']  # reference df  for difference
        for k in keylist[1:]:
            df     = MASTER_DF_DICT[k]['df']
            plt_alpha  = df_dict['alpha']
            if DIFF == 'CID':
                #need to do an inner join with each entry in dic, then plot the diff
                # (join logic thanks to Charlie Prior)
                join = df_ref.join(df.set_index('CID'), on='CID', how='inner', lsuffix='_1', rsuffix='_2')
                plt.scatter(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values,
                            alpha=plt_alpha, label='Diff')
                avgdiff = binned_statistic(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values, bins=bins, statistic='median')[0]                                     
                plt.scatter((bins[1:] + bins[:-1])/2, avgdiff, label="Mean Difference", color='k')
            elif (DIFF == 'ALL'):
                text_label = df_ref.name.values[0]+ " - " + df.name.values[0]
                try:
                    plt.scatter(df_ref.x_plot_val, df_ref.y_plot_val - df.y_plot_val, 
                                label=text_label, alpha=plt_alpha)
                except ValueError:
                    pass

                if do_median:
                    median_ref = binned_statistic(df_ref.x_plot_val, df_ref.y_plot_val,
                                                  bins=bins, statistic='median')[0]
                    median = binned_statistic(df.x_plot_val, df.y_plot_val,
                                              bins=bins, statistic='median')[0]
                    plt.scatter((bins[1:] + bins[:-1])/2., median_ref - median, 
                                label=text_label+" median", marker="^", zorder=10)
                    
            else:  
                sys.exit(f"\n ERROR: {str(DIFF)} is not a valid DIFF option -> ABORT.")         

            setup_plot(args, plot_title, xlabel, ylabel+" diff")  
    else:
        # 2D for each file
        try:
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]]) 
        except KeyError:
            pass  # auto scale axis
        
        for key_name, df_dict in MASTER_DF_DICT.items():

            df          = df_dict['df']
            name_legend = df_dict['name_legend']
            plt_alpha   = df_dict['alpha']
            nevt        = len(df) 
            size        = 20 / math.log10(nevt)  # dot size gets smaller with nevt ??

            plt_legend = name_legend
            if do_nevt: plt_legend += f'  N={nevt}'
            
            #print(f"\n xxx nevt= {nevt}  size={size}\n")
            plt.scatter(df.x_plot_val, df.y_plot_val, alpha=plt_alpha, label=plt_legend,
                        zorder=0, s=size)

            # overlay information on plot
            if do_median:
                median = binned_statistic(df.x_plot_val, df.y_plot_val,
                                          bins=bins, statistic='median')[0] 
                plt.scatter((bins[1:] + bins[:-1])/2., median,
                            label= name_legend+" median", marker="^", zorder=10)

            if do_diag_line:
                x = np.linspace(xmin,xmax,100);  y = x
                plt.plot(x,y)

            if do_list_cid:
                print_cid_list(df, name_legend)
            
        setup_plot(args, plot_title, xlabel, ylabel)    

        
    return 
    # end plotter_func

def setup_plot(args, plot_title, xlabel, ylabel):

    # Created July 21 2024 by R.Kessler
    # wrapper for lots of plt. calls that are repeated several times.
    
    OPT          = args.OPT    
    do_logy      = OPT_LOGY    in OPT
    do_grid      = OPT_GRID    in OPT
    
    if do_logy:    plt.yscale("log")
    if do_grid:    plt.grid()

    plt.xlabel(xlabel)  
    if ylabel is not None: plt.ylabel(ylabel)
    
    plt.legend()                                
    plt.title(plot_title) 

def print_cid_list(df, name_legend) :
    # print list of cids to stdout    
    cid_list = sorted(df['CID'].to_numpy())[0:NMAX_CID_LIST]
    print(f"\n CIDs passing cuts for '{name_legend}' : \n{cid_list}" )
    sys.stdout.flush()
    return


# ===================================================
#   Add main, June 2024
# ========================================
if __name__ == "__main__":

    setup_logging()
    logging.info("# ========== BEGIN plot_table.py  ===============")

    args = get_args()

    # process/modify input args
    process_args(args)

    V = args.VARIABLE
    args.VARIABLE = translate_VARIABLE(V)   # add df. as needed to args.VARIABLE
    
    for icut, cut in enumerate(args.cut_list):
        args.cut_list[icut] = translate_CUT(cut)  # add df. and df.loc as needed
        
    if  DEBUG_TRANSLATE :
        sys.exit(f"\n xxx bye .")
    
    plot_info = Namespace()  # someplace to store internally computed info
    
    set_var_dict(args, plot_info) # set plot bounds and axisinfo

    read_tables(args, plot_info)  # read each input file and store data frames
    #sys.exit(f"\n xxx DEBUG STOP xxx")
    
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
    

    
