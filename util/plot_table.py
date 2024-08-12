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
OPT_NEVT      = "NEVT"      # append N={nevt} to legend
OPT_AVG       = "AVG"       # 1D->append avg to legend; 2D->overlay avg in x-bins
OPT_MEAN      = "MEAN"      # same as AVG
OPT_STDDEV    = "STDDEV"    # append stddev to legend
OPT_CHI2      = "CHI2"      # show tfile1/tfile2 chi2/dof and scale tfile2 to match tfile1
OPT_MEDIAN    = "MEDIAN"
OPT_DIAG_LINE = "DIAG_LINE"  # draw diagonal line on plot
OPT_LOGY      = "LOGY"       # log scale along Y axis
OPT_GRID      = "GRID"
OPT_LIST_CID  = "LIST_CID"   # list CIDs passing cuts

VALID_OPT_LIST = [ OPT_NEVT, OPT_AVG, OPT_MEAN, OPT_STDDEV, OPT_CHI2, OPT_CHI2,
                   OPT_MEDIAN, OPT_DIAG_LINE,
                   OPT_LOGY, OPT_GRID, OPT_LIST_CID ]

NMAX_CID_LIST = 20  # max number of CIDs to print for @@OPT CID_LIST

ARG_DIFF_CID = "CID"
ARG_DIFF_ALL = "ALL"
VALID_ARG_DIFF_LIST = [ ARG_DIFF_CID, ARG_DIFF_ALL ]

# define dictuionay of user-defined functions (key)
# and the np.xxx replacement for pandas. These functions
# can be used in variables (@V) and weigts (@@WEIGHT).
NUMPY_FUNC_DICT = {
    'exp'         :  'np.exp'  ,
    'log'         :  'np.log'  ,   # works for log and log10
    'sqrt'        :  'np.sqrt' ,
    'abs'         :  'np.abs'  ,
    'heaviside'   :  'np.heaviside'
}

        
# internal DEBUG flags
DEBUG_FLAG_REFAC          =  2
DEBUG_FLAG_LEGACY         = -2
DEBUG_FLAG_DUMP_TRANSLATE =  3


# ================================

def print_help():

    help_string = \
f"""
This plot unility works on 
  * FITRES table files create by LC-fitting code (snlc_fit.exe) and BBC (SALT2mu.exe)
  * M0DIF files from BBC
  * HOSTLIB files used in simulation
  * any file with same format that has VARNAMES key

BEWARE that conventional dashes for command-line input keys are replaced 
with @@ ; e.g., --VARIABLE in any other python code is @@VARIABLE here ...
this change avoids confusing dashes and minus signs.

       IF A COMMAND FAILS, PLEAE POST SNANA-GITHUB ISSUE !!! 

There are two general types of comamnd-line input:
  1. plot content
  2. presentation style

and two types of command-line input delimeters
  1. colon separates args for x and y axes; 
      (e.g.,  @V zHD:c   or    @@UNITS day:degrees)
  2. space separates args for mutliple plots 
      (e.g., @@TFILE A.TXT B.TXT  or  @@CUT 'zHD<0.2' 'zHD<0.4' )


      INPUTS FOR PLOT CONTENT 
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~

@@VARIABLE or @V:
  Input table files are loaded with pandas dataframes, and therefore input variables
  are internally converted to pandas notation. For instance, to plot redshift 
  distribution requires user input
     @@VARIALBE zHD

  2D plots are parsed by a colon (:) to separate the x and y variables. 
  The syntax is generally x:y For instance, 
     @@VARIABLE zHD:x1      # displays x1 vs z. 

  In short, this code internally prepends pandas notation (df.) to simplify
  command-line inputs, This script prints the pandas-translated @@VARIABLE 
  string (with df.) to stdout.

  Algebraic and numpy functions are also allowed, e.g., 
      @@VARIABLE zHD:mB - 3.1*c + 0.16*x1
      @@VARIABLE 'sqrt(MUERR**2 + .05**2):zHD'
      @@VARIABLE 'np.sqrt(MUERR**2 + .05**2):zHD'  # can explicitly define with np.

  Single quotes around expressions with parenetheses are needed to avoid 
  linux problems parsing (). Math functions (sqrt, abs, exp, log, log10)
  are internally updated with 'np.' prefix. There is a hard-wired list of 
  functions to check for missing np, so if using an undefined function
  you can explicitly prepend np (and please post github issue about 
  missing function).

@@ERROR @E
  for 2D plot, specify variable(s) to plot as error bar; 
  e.g, plot hubble diagram as
    @@VARIABLE zHD:MU  @@ERROR zHDERR:MUERR
      or
    @@VARIABLE zHD:MU  @@ERROR  :MUERR  # include MU error bar but no zHD error bar
    
  When @@ERROR defines y-axis error, '@@OPT MEAN' (see below) replaces the default 
  arithmetic mean with 1/ERR^2-weighted mean. 

@@error @e
  For plots with many data points, showing many error bars can result in a plot that 
  looks like a painted wall and no detail can be seen among the individual data points. 
  To use errors for weight-average computation and plot data points WITHOUT error bars,
  replace @@ERROR @E --> @@error @e.

  
@@CUT 
  Apply selection cuts on the sample. As with @@VARIABLE, CUT is internally 
  translated to append df and df.loc as needed. Beware that quotes (single or 
  double) are required around all cuts: e.g,

     @@CUT "zHD > 0.2 &  zHDERR<.1 & SNRMAX1 < 100 & FIELD='C3'"

  If a cut includes a string matches (e.g., FIELD='C3'), single quotes must be 
  used around the string to match, and double quotes around the entire CUT arg.
  To overlay the same variable with different cuts, provide a list of cuts,
     @@CUT  'SNRMAX1>10'  'SNRMAX1>20'  'SNRMAX1>40'
  results in 3 overlaid plots, and default legend shows each cut.
  While math functions (sqrt, exp, log ...) are allowed for variables,
  they cannot be used for cuts.

@@WEIGHT
  Reweight 1D bin contents with arbitrary math function of x-axis using
  syntax with 'x' as variable:
     @@WEIGHT '1-x+x**2'
     @@WEIGHT 'exp(-0.4*((x+.3)/.2)**2)'
     @@WEIGHT 'exp(-0.4*((x+.3)/.2)**2)*heaviside(x-2,0.5)'

  A weighted plot can be overlaid on original plot (single @@TFILE arg)
  with list of two weight functions, 
     @@TFILE A.TEXT  @@WEIGHT  1   '1-x+x**2'
  where the "1" arg means that first plot is not modified. Finally,
  when overlaying plots from two different files, providing two 
  weights applies a separate weight to each file, e.g.
     @@TFILE A.TEXT B.TEXT  @@WEIGHT 1  'exp(-x/2)'
  aplplies unit weight to A.TEXT and exp(-x/2) weight to B.TEXT.

  Vertical axis label shows Counts * WEIGHT.  Quotes are required 
  around @@WEIGHT arg to avoid conflicts with unix commands.
  Internally, the plot script applies "np." where needed; e.g. 
  exp is replaced with np.exp, and heaviside -> np.heaviside.
  

@@BOUNDS
  Custom axis boundaries are input with 
     @@BOUNDS xmin xmax xbin                   # 1D
     @@BOUNDS xmin xmax xbin: ymin ymax ybin   # 2D; ybin ignored
  Mean, Median, stdev only include entries within the plot bounds; 
  overflows are ignored.

@@DIFF
  For 2D plots with 2 table files (or 2 sets of cuts), compare 
  values using;
    @@DIFF ALL    #  plot y-axis difference in median values
                  #  (can compare correlated or independent samples)
      or
    @@DIFF CID    #  plot y-axis difference for each CID.
                  #   (compare correlated samples only; must have CID overlap)


      INPUTS FOR PLOT STYLE
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~

@UNITS
  Default x- and y-axis labels show variable name only.
  To display units in (),
     @@UNITS days       # units for 1D plot or x-axis of 2D plot
     @@UNITS :deg       # units for y-axis of 2D plot (no unit for x-axis)
     @@UNITS days:deg   # units for both axes of 2D plot
       or
     @U <arg>

@@ALPHA
  Aadjust the matplotlib alpha values to adjust transparency 
  (0=transparent, 1=solid). For 2D overlays of multuple files or cuts, 
  can specify multiple alpha values, e.g.,
     @@ALPHA 0.9 0.1
  so that primary plot is dark and overlay plot is nearly transparent.

@@MARKER
  Override default solid circle markers with
     @@MARKER s ^  # square and triangle for 1st and 2nd file/cut
     @@MARKER s x  # square and X        for 1st and 2nd file/cut
     etc ... (see online doc on matplotlib markers)

@@SAVE
  Save figure as pdf or png; e.g.
    @@SAVE my_first_table_plot.pdf
       or
    @@SAVE my_first_table_plot.png
  Plot is either displayed in pop-up window or stored in @@SAVE file ...
  but not both. This enables pipelines to make post-processing plots
  without worrying about pop-windows in slurm jobs.


@@OPT   {' '.join(VALID_OPT_LIST)}
    where
      {OPT_NEVT:<12} ==> append N=Nevt on each legend (1D and 2D).
      {OPT_MEAN:<12} ==> 1D->append mean on legend; 2D->overlay mean in x-bins.
                         If @@ERROR defines y-axis error, replace arithmetic
                         mean with 1/ERR^2-weighted mean.
      {OPT_AVG:<12} ==> same as {OPT_MEAN}.
      {OPT_MEDIAN:<12} ==> overlay median in x-bins (2D only).
      {OPT_STDDEV:<12} ==> append stddev on each legend (1D only).
      {OPT_CHI2:<12} ==> display chi2/dof on plot for two table files (1D only).
      {OPT_DIAG_LINE:<12} ==> draw line with slope=1 for 2D plot.
      {OPT_LOGY:<12} ==> log scale for vertical axis.
      {OPT_GRID:<12} ==> display grid on plot.
      {OPT_LIST_CID:<12} ==> print up to 100 CIDs passing cuts.

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
    parser.add_argument('@@TFILE', '@@tfile', help=msg, nargs="+")
    
    msg = "required: Variable(s) to plot from table file, or function of variables." \
          "For the histogram, counts are normalised to the first table file."
    parser.add_argument('@V', '@@VARIABLE', '@v', help=msg, nargs="+")

    msg = 'variable to use for error bar in 2D plot, and for weighted avg per bin'
    parser.add_argument('@E', '@@ERROR', help=msg, nargs="+", default=None)

    msg = 'variable to only for weight avg per bin; do not show error bard in 2D plot'
    parser.add_argument('@e', '@@error', help=msg, nargs="+", default=None)        

    msg = "cuts with boolean and algegraic operations. If >1 CUT, plots are overlaid."
    parser.add_argument("@@CUT", "@@cut", help=msg, nargs="+", default =[None])
    
    msg = "Plot difference in y-axis between files.\n" \
          "ALL plots median difference between first file and subsequent ones.\n" \
          "CID plots per-CID difference."
    parser.add_argument('@@DIFF', '@@diff', default=None, type=str, help=msg )

    
    msg = "WEIGHT function(s) for 1D hist only; e.g., 1-0.5*x+3*x**2"
    parser.add_argument("@@WEIGHT", "@@weight", help=msg, nargs="+", default=[None] ) 

    msg = "number of rows to read (for large files)"
    parser.add_argument("@@NROWS", "@@nrows", help=msg, type=int, default=0)
    
    # -------- decorative options below ---------
    
    msg = "AUTO (default): optional bounds min, max and binsize of plot parameters.\n" \
          "For 2D plot, must specify both x and y bounds. y-binsize is ignored."
    parser.add_argument('@@BOUNDS', '@@bounds', default=BOUNDS_AUTO, help=msg, nargs = '+')

    msg = "Units to show for each axis; see @@HELP"
    parser.add_argument('@U', '@@UNITS', '@@units', default=None, help=msg, nargs="+")
    
    msg = "Override default plot title"
    parser.add_argument("@@TITLE", "@@title", default=None, help=msg )

    msg = "Override default legend on plot (space sep list per TFILE)"
    parser.add_argument('@@LEGEND', '@@legend', default=None, help=msg, nargs="+")

    msg = "Override default marker='o'"
    parser.add_argument('@@MARKER', '@@marker', default=['o'], help=msg, nargs="+")    

    msg = "Alpha value for plot. Set to 0 to see only averages." \
          "ALPHA=0 and DIFF=True compares average difference between two files, " \
          "even if there are no overlapping CIDS."
    parser.add_argument('@@ALPHA', '@@alpha',  default=[0.8], type=float, help=msg, nargs="+" )
    
    msg = "Filename to save plot.  Default: does not save plots."
    parser.add_argument('@@SAVE', '@@save', default=None, help=msg )

    msg = "options; see @@HELP"
    parser.add_argument("@@OPT", "@@opt", help=msg, nargs="+", default = [])

    msg = "debug options (for development only)"
    parser.add_argument("@@DEBUG_FLAG", "@@debug_flag", help=msg, type=int, default=0)    

    msg = "Full help menu printed to stdout"
    parser.add_argument("@H", "@@HELP", help=msg, action="store_true")

    args    = parser.parse_args()
    
    if args.HELP:
        print_help()    

    return args

    # end get_args()

def arg_prep_driver(args):

    # process arguments such as remove pad spacing and forcing
    # number or @@TFILE args to match number of @@CUT args
    # (and vice versa)

    # ---------------------------------------
    # check debug args (for develepers only)
    arg_prep_DEBUG_FLAG(args)
    
    # - - - - - - -
    # for variable(s), remove pad spacing and add np. if needed
    # for functions
    args.VARIABLE_ORIG = args.VARIABLE 
    if args.VARIABLE:
        args.VARIABLE      = ''.join([str(elem) for elem in args.VARIABLE])        
        args.VARIABLE_ORIG = args.VARIABLE_ORIG[0]
        args.VARIABLE = numpy_fun_replace(args.VARIABLE)
    else:
        sys.exit(f"\n ERROR: must define variable(s) to plot with @@VARIABLE or @@V")

    # CUT is tricky. Make sure that length of cut list matchs length
    # of table-file list ... or vice versa ... make sure that length of
    # tfile list matches length of CUT list.
    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx proc_args: args.CUT = {args.CUT}  (before modifications)")

    n_tfile_orig    = len(args.TFILE)
    n_cut_orig      = len(args.CUT)
    n_weight_orig   = len(args.WEIGHT)
    tfile_list      = copy.copy(args.TFILE)
    cut_list        = copy.copy(args.CUT)
    weight_list     = copy.copy(args.WEIGHT)

    name_arg_list  = [ '@@TFILE', '@@CUT', '@@WEIGHT' ]
    n_orig_list    = [ n_tfile_orig, n_cut_orig, n_weight_orig ]
    arg_list_list  = [ tfile_list,   cut_list,   weight_list   ]

    # abort if more than 1 table file and more than 1 cut are requested.
    # However, allow 2 table files and 2 WEIGHTs 
    n_multiple = sum(n > 1 for n in n_orig_list)
    if n_tfile_orig> 1 and n_cut_orig > 1:
        sys.exit("\n ERROR: cannot specify >1 tfile and >1 cut: \n" \
                 f"\t n_plots = {n_orig_list} for {name_arg_list}")

    if n_multiple >= 1:
        n_plot = max(n_orig_list)
        for j, arg_list in enumerate(arg_list_list):
            if len(arg_list) == 1:
                arg_list_list[j] = [ arg_list_list[j][0] ] * n_plot

    # update lists that will be used later
    tfile_list   = arg_list_list[0]
    cut_list     = arg_list_list[1]
    weight_list  = arg_list_list[2]    
        
    # - - - -

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx proc_args: args.TFILE  -> {args.TFILE}")
        print(f" xxx proc_args: args.CUT    -> {args.CUT}")
    
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
    narg_tfile           = len(args.tfile_list)
    
    args.tfile_base_list = table_base_list
    args.cut_list        = cut_list
    args.weight_list     = weight_list
    # - - - - -

    # make sure there is a colon in ERROR and UNITS args
    # e.g  "degrees" -> "degrees:" so that user can specify only x-axis
    # without colon and code here adds required colon. For y-axis,
    # user must provide colon; e.g.,  :degress applies y-axis label.
    args.use_err_list = True  # default
    if args.error:
        args.ERROR = args.error
        args.use_err_list = False  # disable error bar on 2D plot
    
    args.UNITS = arg_prep_axis(args.UNITS)
    args.ERROR = arg_prep_axis(args.ERROR)
    
    if args.BOUNDS != BOUNDS_AUTO:
        args.BOUNDS = ' '.join([str(elem) for elem in args.BOUNDS])


    args.legend_list = arg_prep_legend(args)
    
    # if only 1 alpha, make sure there is alpha for each file/cut
    args.alpha_list  = arg_prep_extend_list(narg_tfile, args.ALPHA)
    args.marker_list = arg_prep_extend_list(narg_tfile, args.MARKER) 

    # if user has not specified any kind of statistical average
    # for DIFF, then set default MEDIAN
    if args.DIFF:
        OPT = args.OPT
        do_tmp = OPT_MEDIAN in OPT  or OPT_MEAN in OPT or OPT_AVG in OPT
        if not do_tmp:
            args.OPT.append(OPT_MEDIAN) 

    # - - - - - - -
    args.OPT = arg_prep_OPT(args)            
    
    return  # end process_args

def arg_prep_OPT(args):

    # abort on invalid @@OPT, and change all OPT to upper case    
    invalid_OPT_list = []
    args_upper_list = []
    OPT_orig = args.OPT
    OPT_out  = OPT_orig  # default
    
    if OPT_orig:
        for opt in OPT_orig:
            opt = opt.upper()
            args_upper_list.append(opt)
            if opt not in  VALID_OPT_LIST:
                invalid_OPT_list.append(opt)
        if len(invalid_OPT_list) > 0:
            sys.exit(f"\n ERROR: Invalid @@OPT {invalid_OPT_list} ; " \
                     f"Valid @@OPT args are\n\t {' '.join(VALID_OPT_LIST)}")

        OPT_out = args_upper_list

    return OPT_out

def arg_prep_DEBUG_FLAG(args):

    args.DEBUG_FLAG_REFAC          = True
    args.DEBUG_FLAG_LEGACY         = False
    args.DEBUG_FLAG_DUMP_TRANSLATE = False
    if args.DEBUG_FLAG != 0:
        args.DEBUG_FLAG_LEGACY         = args.DEBUG_FLAG == DEBUG_FLAG_LEGACY
        args.DEBUG_FLAG_REFAC          = not args.DEBUG_FLAG_LEGACY        
        args.DEBUG_FLAG_DUMP_TRANSLATE = args.DEBUG_FLAG == DEBUG_FLAG_DUMP_TRANSLATE
        logging.info(f"# \t DEBUG_FLAG={args.DEBUG_FLAG} " )
        logging.info(f"# \t DEBUG_FLAG_REFAC          = {args.DEBUG_FLAG_REFAC}")
        logging.info(f"# \t DEBUG_FLAG_DUMP_TRANSLATE = {args.DEBUG_FLAG_DUMP_TRANSLATE}")
        print(f"")        
    
    return

def arg_prep_legend(args):

    LEGEND_orig = args.LEGEND
    LEGEND_out  = LEGEND_orig   # default is user input
    
    if LEGEND_orig is None:
        # no user supplied legend, so make up reasonable legend
        LEGEND_out = [ None ] * len(args.TFILE)
        if len(args.CUT) > 1 :
            LEGEND_out = args.CUT
        else :
            LEGEND_out = []
            for t in args.tfile_base_list:
                legend = t.split('.')[0]  # base file name without extension after dot
                LEGEND_out.append(legend)
    # - - - - - 
    narg_legend  = len(LEGEND_out)
    narg_tfile   = len(args.tfile_list)    
    match_narg   = narg_tfile == narg_legend
    if  not match_narg:
        sys.exit(f"ERROR: narg_tfile={narg_tfile} but narg_legend={narg_legend}; " \
                 f"narg_legend must match number of table files.")

    return LEGEND_out

def arg_prep_extend_list(narg_tfile, arg_list_orig):

    # if arg_list_orig has fewer elements than number of table files,
    # extend list to be same length as number or table files;
    # else return original list.  This allows users to specify one
    # value (e.g. for alpha or maker) that is used on all plots.
    
    narg_orig = len(arg_list_orig)
    if narg_orig < narg_tfile:
        val          = arg_list_orig[0]
        arg_list_out = [val]*narg_tfile
    else:
        arg_list_out = copy.copy(arg_list_orig)

    return arg_list_out

    
def arg_prep_axis(arg_axis):
    # if input arg_axis = "degress" then return "degrees:"
    # to ensure that splitting by colon works later.
    # If arg_axis is None, return ":" to indicate nothing for either axis.
    if arg_axis:
        arg_axis  = ''.join([str(elem) for elem in arg_axis])
        if COLON not in arg_axis:
            arg_axis += COLON
    else:
        arg_axis = COLON    
    return arg_axis


def numpy_fun_replace(var)  :
    # replace user functions (e.g., exp or sqrt) with numpy
    # functions (e.g. np.exp or np.sqrt). If np.xxx is already
    # defined, then avoid replacement.

    for fun, np_fun in NUMPY_FUNC_DICT.items():            
        if fun in var and np_fun not in var:
            var = var.replace(fun,np_fun)
    return var

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

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
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

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx CUT      = {CUT}" )
        print(f" xxx cut_list = {cut_list}")

    # for each cut item, append parentheses and df.
    cut_list_df = []

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

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx cut_list_df = {cut_list_df}")
        
    # replace each original cut with cut_df in CUT
    CUT_df = CUT
    for cut, cut_df in zip(cut_list, cut_list_df):
        CUT_df = CUT_df.replace(cut,cut_df)
        if args.DEBUG_FLAG_DUMP_TRANSLATE:
            print(f" xxx \t replace user cut {cut} --> {cut_df} ")

    # replace input for '=' with '==', but be careful not change '!='
    # This logic is a bit goofy; is there a  simpler logic?
    if '=' in CUT_df:
        CUT_df = CUT_df.replace('=', '==')
    if '!==' in CUT_df:
        CUT_df = CUT_df.replace('!==', '!=')        
            

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
    ERROR    = args.ERROR
    BOUNDS   = args.BOUNDS
    UNITS    = args.UNITS
    
    # store x and [optional] y variable names
    plotdic            = {}  # x & y varnames, err varnames, axis labels

    if args.TITLE:
        plot_title    = str(args.TITLE)
    elif args.CUT :
        plot_title    = str(args.CUT)
    else:
        plot_title    = str(args.VARIABLE_ORIG)
    
    if STR_np in plot_title:
        plot_title = plot_title.replace(STR_np,'')

    # - - - - 
    VAR_LIST     = VARIABLE.split(COLON)
    VARERR_LIST  = ERROR.split(COLON)
    UNITS_LIST   = UNITS.split(COLON)

    for n, VAR in enumerate(VAR_LIST):

        STR_VAR      = str(VAR)
        STR_LABEL    = STR_VAR.replace(STR_df,'') 
        if STR_np in STR_LABEL:
            STR_LABEL = STR_LABEL.replace(STR_np,'')

        VARERR  = VARERR_LIST[n]
        if len(VARERR) > 0:
            STR_VARERR   = STR_df + str(VARERR)
        else:
            STR_VARERR   = None            
        
        U = UNITS_LIST[n]
        if len(U) > 0:  STR_LABEL += '  (' + U + ')'
            
        if n == 0:
            xlabel = STR_LABEL
            ylabel = 'Counts'
            if args.WEIGHT[0]: ylabel = f'Counts * {args.WEIGHT}'            
            plotdic['x']            = STR_VAR
            plotdic['xerr']         = STR_VARERR
            plotdic['xaxis_label']  = xlabel
            plotdic['yaxis_label']  = ylabel
            plotdic['ndim']         = 1
            
        else:
            plotdic['y']            = STR_VAR
            plotdic['yerr']         = STR_VARERR
            plotdic['yaxis_label']  = STR_LABEL
            plotdic['ndim']         = 2
            
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

    # - - - - - -  - - 
    # error checks 
    if args.DIFF:
        if plotdic['ndim'] == 1 :
            sys.exit("\nERROR: DIFF does not work for 1D histograms." \
                     " ABORT to avoid confusion.")
        if args.DIFF not in VALID_ARG_DIFF_LIST :
            sys.exit(f"\n ERROR: @@DIFF {args.DIFF} is not a valid option; " \
                     f" valid options are {VALID_ARG_DIFF_LIST}")
        if len(args.TFILE) == 1:
            sys.exit(f"\n ERROR '@@DIFF {args.DIFF}' does not work with 1 table file;\n"\
                     f"\t need 2 or more table files specified after @@TFILE key.")
            
    # load output namespace
    plot_info.plotdic             = plotdic
    plot_info.plot_title          = plot_title
    plot_info.boundsdic           = boundsdic
    plot_info.custom_bounds       = custom_bounds
    
    return
    # end set_var_dict

def read_tables(args, plot_info):

    tfile_list      = args.tfile_list
    tfile_base_list = args.tfile_base_list
    cut_list        = args.cut_list
    weight_list     = args.weight_list
    legend_list     = args.legend_list
    alpha_list      = args.alpha_list
    marker_list     = args.marker_list
    NROWS           = args.NROWS

    plotdic    = plot_info.plotdic
    boundsdic  = plot_info.boundsdic
    
    MASTER_DF_DICT = {}  # dictionary of variables to plot (was MASTERLIST)
    nf = 0
    
    for tfile, cut, weight, legend, alpha, marker in \
        zip(tfile_list, cut_list, weight_list, legend_list, alpha_list, marker_list):
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
            'weight'       : weight,
            'name_legend'  : name_legend,
            'alpha'        : alpha,
            'marker'       : marker
        }

        try:
            MASTER_DF_DICT[key]['df']['x_plot_val'] = eval(plotdic['x'])
            if plotdic['xerr'] is not None:
                MASTER_DF_DICT[key]['df']['x_plot_err'] = eval(plotdic['xerr'])       
            boundsdic[key + "_min"] = np.amin(MASTER_DF_DICT[key]['df']['x_plot_val'])
            boundsdic[key + "_max"] = np.amax(MASTER_DF_DICT[key]['df']['x_plot_val'])
            if plotdic['ndim'] == 2:
                MASTER_DF_DICT[key]['df']['y_plot_val'] = eval(plotdic['y'])
                if plotdic['yerr'] is not None:
                    MASTER_DF_DICT[key]['df']['y_plot_err'] = eval(plotdic['yerr'])
                
        except AttributeError:
            sys.exit(f"\n ERROR: Couldn't set bounds for plotdic={plotdic} and {tfile}")

            
    # - - - - - - - -
    logging.info("Done loading all table files.")
    logging.info("# - - - - - - - - - - ")
    
    # load output namespace
    plot_info.MASTER_DF_DICT = MASTER_DF_DICT
    
    return
    # end read_tables

    

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

def plotter_func_driver(args, plot_info):

    # REFAC Aug 10 2024
    # utility to create the plot (but doesn't show it)

    # strip off local args from input name spaces
    DIFF      = args.DIFF
    OPT       = args.OPT

    do_chi2      = OPT_CHI2      in OPT
    do_median    = OPT_MEDIAN    in OPT
    do_diag_line = OPT_DIAG_LINE in OPT
    do_list_cid  = OPT_LIST_CID  in OPT
    do_nevt      = OPT_NEVT      in OPT
    do_avg       = OPT_AVG       in OPT or OPT_MEAN in OPT
    do_stddev    = OPT_STDDEV    in OPT    
    
    MASTER_DF_DICT       = plot_info.MASTER_DF_DICT
    plotdic              = plot_info.plotdic
    boundsdic            = plot_info.boundsdic
    custom_bounds        = plot_info.custom_bounds
    plot_title           = plot_info.plot_title
    
    xlabel               = plotdic['xaxis_label']
    ylabel               = plotdic['yaxis_label']
    NDIM_PLOT            = plotdic['ndim']  # 1 or 2

    xbins = None
    ybins = None    
    if custom_bounds:                       
        xmin  = boundsdic['x'][0]; xmax = boundsdic['x'][1]
        xbin  = boundsdic['x'][2]
        xbins = np.arange(xmin, xmax, xbin)
        if NDIM_PLOT == 2:
            ymin = boundsdic['y'][0]; ymax = boundsdic['y'][1]
            ybin = boundsdic['y'][2]
            ybins = np.arange(ymin, ymax, ybin) # ignored 
    else:                                   
        xbins = np.linspace(boundsdic[min(boundsdic, key=boundsdic.get)],
                            boundsdic[max(boundsdic, key=boundsdic.get)], 30)
        xmin = xbins[0];  xmax = xbins[-1]

    xbins_cen   = ( xbins[1:] + xbins[:-1] ) / 2.  # central value for each xbin 
    # - - - - - - - 
    msg = f"The min & max x-axis bounds are: " \
          f"{np.around(xbins[0],4)}  {np.around(xbins[-1],4)} respectively"
    logging.info(msg)

    numplot = 0 
    info_plot_dict = { }
    
    for key_name, df_dict in MASTER_DF_DICT.items():

        df          = df_dict['df'] 
        name_legend = df_dict['name_legend']        
        logging.info(f"Plot {name_legend}")

        info_plot_dict['do_plot_errorbar']    = False
        info_plot_dict['do_ov1d_hist']        = False
        info_plot_dict['do_ov2d_binned_stat'] = False
        info_plot_dict['numplot']             =  numplot
        info_plot_dict['xbins']               =  xbins
        info_plot_dict['ybins']               =  ybins
        info_plot_dict['xbins_cen']           =  xbins_cen        
        info_plot_dict['df_dict']             =  df_dict        
                
        if NDIM_PLOT == 1:
            # 1D
            get_info_plot1d(args, info_plot_dict)        
        else:
            # 2D
            get_info_plot2d(args, info_plot_dict)      

        # strip off arguments to pass to matplotlib ...
        do_plot_errorbar    = info_plot_dict['do_plot_errorbar']
        do_ov1d_hist        = info_plot_dict['do_ov1d_hist']
        do_ov2d_binned_stat = info_plot_dict['do_ov2d_binned_stat']
        
        xval_list   = info_plot_dict['xval_list'] 
        yval_list   = info_plot_dict['yval_list']

        if args.use_err_list:        
            xerr_list   = info_plot_dict['xerr_list'] 
            yerr_list   = info_plot_dict['yerr_list']
        else:
            xerr_list = None
            yerr_list = None
                
        plt_size    = info_plot_dict['plt_size']   # depends on nevt for 2D
        plt_legend  = info_plot_dict['plt_legend'] # can be appended with more info
        plt_text_dict = info_plot_dict['plt_text_dict']
        
        plt_alpha   = df_dict['alpha']  # fixed by user
        plt_marker  = df_dict['marker'] # fixed by user
        # - - - - -
        
        if custom_bounds:
            plt.xlim(xmin, xmax)
            if NDIM_PLOT==2:  plt.ylim(ymin, ymax)
        
        if do_plot_errorbar :
            # nominal
            plt.errorbar(xval_list, yval_list, xerr=xerr_list, yerr=yerr_list, 
                         fmt=plt_marker, label=plt_legend,
                         markersize=plt_size, alpha=plt_alpha )
            
            if NDIM_PLOT==2 and do_ov2d_binned_stat :
                overlay2d_binned_stat(args, info_plot_dict)
            
        elif NDIM_PLOT==1 and do_ov1d_hist:
            # 1D, typically sim overlaid on data
            wgt_ov = info_plot_dict['wgt_ov']
            plt.hist(df.x_plot_val, xbins, alpha=0.25, weights = wgt_ov,
                     label = plt_legend)
        else:
            # nothing to plot; e.g, 1st file for DIFF option
            numplot += 1
            continue

        # - - - - -
            
        # - - - -
        # check for misc options
        if do_list_cid:
            print_cid_list(df, name_legend)

        if do_diag_line:
            x = np.linspace(xmin,xmax,100);  y = x
            plt.plot(x,y)

        apply_plt_misc(args, plot_title, xlabel, ylabel, plt_text_dict)
        numplot += 1
        
        
    return   # end plotter_func_driver


def get_info_plot1d(args, info_plot_dict):

    do_chi2   = OPT_CHI2 in args.OPT
    do_nevt   = OPT_NEVT in args.OPT
    do_avg    = OPT_AVG  in args.OPT or OPT_MEAN in args.OPT
    do_stddev = OPT_STDDEV in args.OPT
    
    numplot   = info_plot_dict['numplot']
    xbins     = info_plot_dict['xbins']
    xbins_cen = info_plot_dict['xbins_cen']    
    df_dict   = info_plot_dict['df_dict']
    
    df             = df_dict['df']
    weight         = df_dict['weight']
    name_legend    = df_dict['name_legend']
    plt_legend     = name_legend
    plt_text_dict  = None
    
    xval_list = xbins_cen
    xerr_list = None
    yval_list = binned_statistic(df.x_plot_val, df.x_plot_val, 
                                 bins=xbins, statistic='count')[0]                

    do_plot_errorbar = True
    do_ov1d_hist     = False

    # apply option user-weight function (see @@WEIGHT arg)
    if weight:
        wgt_user   = get_weights_user(xbins_cen, weight)
        yval_list *= wgt_user

    nevt = np.sum(yval_list)  # sum before re-normalizing (for printing only)
    errl,erru = poisson_interval(yval_list) # compute poisson errors
        
    if numplot == 0 :
        info_plot_dict['yval0_list']   = yval_list # save for normalizing overlay
        info_plot_dict['name0_legend'] = name_legend
    else:
        if do_chi2:
            # re-scale overlay plot only if chi2 option is requested
            yval0_list   = info_plot_dict['yval0_list']
            name0_legend = info_plot_dict['name0_legend'] 
            scale = np.sum(yval0_list) / np.sum(yval_list)
            do_plot_errorbar = False
            do_ov1d_hist     = True
        else:
            scale = 1.0  # overlay plot is not scaled

        yval_list  *= scale # normalize integral to match file 0 
        errl       *= scale
        erru       *= scale
        logging.info(f"Overlay {name_legend} scaled by {scale:.3e}")

    # ---------------

    # check for overlay weights with chi2
    if numplot > 0 and do_chi2:
        wgt_ov = [ scale ] * len(df.x_plot_val)
        if weight :
            wgt_user = get_weights_user(df.x_plot_val,weight) 
            wgt_ov   = np.multiply(wgt_ov,wgt_user)

            len_wgt_ov = len(wgt_ov)
            #sys.exit(f"\n xxx wgt_ov({len_wgt_ov}) = \n{wgt_ov}\n xxx x=\n{df.x_plot_val}\n")
        # check option to compute and print chi2/dof info on plot
        # Froce min error =1 in chi2 calc so that it's ok to plot
        # error=0 for bins with zero events.
        sqdif = (yval0_list-yval_list)**2
        sqerr = np.maximum((yval0_list+yval_list*scale*scale), 1.0)
        chi2  = np.sum( sqdif / sqerr )
        ndof  = len(xbins) - 1 
        text_chi2 = f"chi2/dof = {chi2:.1f}/{ndof}"
        logging.info(f"{name0_legend}/{name_legend} {text_chi2}")
        xmin = xbins[0];  xmax = xbins[-1]
        x_text = xmin + 0.7*(xmax-xmin)
        y_text = 1.0 * np.max(yval0_list)  # warning; fragile coord calc
        plt_text_dict = { 'x_text': x_text, 'y_text': y_text, 'text': text_chi2 }
        #plt.text(x_text, y_text, text_chi2 )
            
    else:
        wgt_ov = None    

    yerr_list = [yval_list-errl, erru-yval_list]
    nevt    = nevt
    avg     = np.mean(df.x_plot_val)
    median  = np.median(df.x_plot_val)
    stddev  = np.std(df.x_plot_val)

    stat_dict = {
        'N'       : [ int(nevt),  do_nevt,   'd'     ],
        'avg'     : [ avg,        do_avg,   '.2f'   ],
        'median'  : [ median,     False,     '.2f'   ],
        'stddev'  : [ stddev,     do_stddev, '.2f'   ]
        # overflow/underflow ??
    }
    for str_stat, tmp_list in stat_dict.items():
        val           = tmp_list[0]
        add_to_legend = tmp_list[1]
        fmt_legend    = tmp_list[2]
        logging.info(f"\t {str_stat:8} value for {name_legend}:  {val:.3f}")
        if add_to_legend:
            plt_legend += f'  {str_stat}={val:{fmt_legend}}'            

    # - - - - - -
    info_plot_dict['do_plot_errorbar']  = do_plot_errorbar
    info_plot_dict['do_ov1d_hist']      = do_ov1d_hist
    info_plot_dict['xval_list']    = xval_list
    info_plot_dict['yval_list']    = yval_list
    info_plot_dict['xerr_list']    = None
    info_plot_dict['yerr_list']    = yerr_list
    info_plot_dict['plt_size']     = None
    info_plot_dict['plt_legend']   = plt_legend
    info_plot_dict['plt_text_dict'] = plt_text_dict
    info_plot_dict['wgt_ov']       = wgt_ov
     
    return  # end get_info_plot1d



def get_info_plot2d(args, info_plot_dict):

    # prepare arguments for matplotlib's plt.errobar.
    
    do_nevt      = OPT_NEVT in args.OPT

    numplot      = info_plot_dict['numplot']
    df_dict      = info_plot_dict['df_dict']
    df           = df_dict['df']
    name_legend  = df_dict['name_legend']
    plt_legend   = name_legend

    nevt         = len(df) 
    plt_size     = 5 / math.log10(nevt)  # dot size gets smaller with nevt 
    if do_nevt: plt_legend += f'  N={nevt}'
    
    xval_list    = df.x_plot_val
    yval_list    = df.y_plot_val
    xerr_list    = None
    yerr_list    = None

    if 'x_plot_err' in df:  xerr_list = df.x_plot_err
    if 'y_plot_err' in df:  yerr_list = df.y_plot_err    

    info_plot_dict['do_plot_errorbar']    = True
    info_plot_dict['do_ov1d_hist']        = False
    info_plot_dict['do_ov2d_binned_stat'] = True
    info_plot_dict['plt_size']      = plt_size
    info_plot_dict['plt_legend']    = plt_legend
    info_plot_dict['plt_text_dict'] = None

    info_plot_dict['xval_list']     = xval_list
    info_plot_dict['yval_list']     = yval_list
    info_plot_dict['xerr_list']     = xerr_list
    info_plot_dict['yerr_list']     = yerr_list

    if args.DIFF:
            
        if numplot == 0:
            # nothing to plot on first file; store references
            info_plot_dict['do_plot_errorbar']   = False  # disable making plot
            info_plot_dict['df_ref']             = df     # store ref table for next plot
            info_plot_dict['name_legend_ref']    = name_legend
            return

        # strip off reference values from first plot
        df_ref          = info_plot_dict['df_ref']
        name_legend_ref = info_plot_dict['name_legend_ref']
            
        # if there is no user-supplied legend, construct legend using
        # auto-generated legend from each plot
        if args.LEGEND is None:
            info_plot_dict['plt_legend']  = name_legend_ref + ' - ' + name_legend
        
        if args.DIFF == ARG_DIFF_CID :
            #need to do an inner join with each entry in dic, then plot the diff
            # (join logic thanks to Charlie Prior)
            join   = df_ref.join(df.set_index('CID'), on='CID', how='inner',
                                 lsuffix='_0', rsuffix='_1')
            xval_list = join.x_plot_val_0.values
            yval_list = join.y_plot_val_0.values - join.y_plot_val_1.values
            info_plot_dict['xval_list']     = xval_list
            info_plot_dict['yval_list']     = yval_list
            if yerr_list is not None:
                yerr_list_0 = join.y_plot_err_0.values
                info_plot_dict['yerr_list']  = yerr_list_0
                
        elif args.DIFF == ARG_DIFF_ALL :
            sys.exit(f"\n ERROR: refactored @@DIFF ALL is not ready ")
            
            # TO DO : ALL option: don't make standard 2D plotl
            #   only plot median/mean difference
        
    return  # end  get_info_plot2d

def overlay2d_binned_stat(args, info_plot_dict):

    # prepare 2D plot overlay of median, avg or wgt-avg in each x-bin.
    
    xbins        = info_plot_dict['xbins']
    xbins_cen    = info_plot_dict['xbins_cen']     
    df_dict      = info_plot_dict['df_dict']
    xval_list    = info_plot_dict['xval_list']
    yval_list    = info_plot_dict['yval_list']
    xerr_list    = info_plot_dict['xerr_list']
    yerr_list    = info_plot_dict['yerr_list']    
    plt_legend   = info_plot_dict['plt_legend'] 
    df           = df_dict['df']

    y_err_stat = None
    
    OPT = args.OPT
    do_median    = OPT_MEDIAN  in OPT
    do_avg       = OPT_MEAN    in OPT or OPT_AVG in OPT
    do_wgtavg    = yerr_list is not None

    DO_DUMP    = False
    
    OPT        = args.OPT
    which_stat = None
    if do_median:
        which_stat = 'median'
    elif do_avg :
        which_stat = 'mean'
    
    if which_stat :
        stat_legend = which_stat  # default for legend
        y_stat  = binned_statistic(xval_list, yval_list,
                                  bins=xbins, statistic=which_stat)[0]

        if do_wgtavg:
            stat_legend = 'wgtavg'
            y_avg   = copy.deepcopy(y_stat)  # for diagnostic print

            # compute error on the mean per bin: stddev/sqrt(N)
            y_std   = binned_statistic(xval_list, yval_list,
                                       bins=xbins, statistic='std')[0]
            y_count = binned_statistic(xval_list, yval_list,
                                       bins=xbins, statistic='count')[0]
            y_err_stat = y_std/np.sqrt(y_count)  # error on mean per bin

            # compute wgted avg using individual errors as weight
            y_wgt_list  = yval_list/yerr_list**2
            wgt_list    = 1.0/yerr_list**2
            sum_y_wgt   = binned_statistic(xval_list, y_wgt_list,
                                           bins=xbins, statistic='sum')[0]
            sum_wgt     = binned_statistic(xval_list, wgt_list,
                                           bins=xbins, statistic='sum')[0]
            y_stat  = sum_y_wgt/sum_wgt

            if DO_DUMP:
                n_val = len(yval_list)
                print(f"")
                print(f" xxx n_val={n_val}")
                print(f" xxx yval_list  = {yval_list[1:10]}")
                print(f" xxx yerr_list  = {yerr_list[1:10]}")   
                print(f" xxx y_wgt_list = {y_wgt_list[1:10]}")
                print(f" xxx wgt_list   = {wgt_list[1:10]}")    
                print(f" xxx binned y_count   = {y_count}")
                print(f" xxx binned sum_y_wgt = {sum_y_wgt}")
                print(f" xxx binned sumwgt    = {sum_wgt}")
                print(f" xxx binned y_avg     = {y_avg}")
                print(f" xxx binned y_wgtavg  = {y_stat}")
                print(f"")            
        # TO DO: error on mean/median ??
        legend  = plt_legend + ' ' + stat_legend
        plt.errorbar(xbins_cen, y_stat, yerr=y_err_stat, fmt='^', label=legend,
                     zorder=5 ) 
            
    return  # end of overlay2d_binned_stat

def plotter_func_legacy(args, plot_info):

    # utility to create the plot (but doesn't show it)

    # strip off local args from input name spaces
    DIFF      = args.DIFF
    OPT       = args.OPT

    do_chi2      = OPT_CHI2      in OPT
    do_median    = OPT_MEDIAN    in OPT
    do_diag_line = OPT_DIAG_LINE in OPT
    do_list_cid  = OPT_LIST_CID  in OPT
    do_nevt      = OPT_NEVT      in OPT
    do_avg       = OPT_AVG       in OPT or OPT_MEAN in OPT
    do_stddev    = OPT_STDDEV    in OPT    

    # XXXXXXXXX LEGACY XXXXXXXXXXXX
    MASTER_DF_DICT       = plot_info.MASTER_DF_DICT
    plotdic              = plot_info.plotdic
    xlabel               = plotdic['xaxis_label']
    ylabel               = plotdic['yaxis_label']
    boundsdic            = plot_info.boundsdic

    custom_bounds        = plot_info.custom_bounds
    plot_title           = plot_info.plot_title

    
    if custom_bounds:                       
        xmin = boundsdic['x'][0]; xmax = boundsdic['x'][1]
        xbin = boundsdic['x'][2]
        bins = np.arange(xmin, xmax, xbin)
        plt.xlim(xmin, xmax)
        if plotdic['ndim'] == 2:
            ymin = boundsdic['y'][0]; ymax = boundsdic['y'][1]
            ybin = boundsdic['y'][2]
            ybins = np.arange(ymin, ymax, ybin) # ignored 
            plt.ylim(ymin, ymax)
    else:                                   
        bins = np.linspace(boundsdic[min(boundsdic, key=boundsdic.get)],
                           boundsdic[max(boundsdic, key=boundsdic.get)], 30)
        xmin = bins[0];  xmax = bins[-1]

    # XXXXXXXXX LEGACY XXXXXXXXXXXX
    
    msg = f"The min & max x-axis bounds are: " \
          f"{np.around(bins[0],4)}  {np.around(bins[-1],4)} respectively"
    logging.info(msg)
    xcen     = ( bins[1:] + bins[:-1] ) / 2.
    
    if plotdic['ndim'] == 1:
        # 1D plot(s)

        for n, key_name in enumerate(MASTER_DF_DICT): 
            df_dict     = MASTER_DF_DICT[key_name]
            df          = df_dict['df']
            weight      = df_dict['weight']  # optoinal weight function
            plt_marker  = df_dict['marker']
            name_legend = df_dict['name_legend']            
            plt_legend  = name_legend

            # get counts sb
            sb = binned_statistic(df.x_plot_val, df.x_plot_val, 
                                  bins=bins,
                                  statistic='count')[0]
                
            if weight:
                wgt_user   = get_weights_user(xcen,weight)
                sb        *= wgt_user

            # XXXXXXXXX LEGACY XXXXXXXXXXXX                
            errl,erru = poisson_interval(sb) # And error for those counts
            nevt   = np.sum(sb)              # nevt before normalization
            avg    = np.mean(df.x_plot_val)
            median = np.median(df.x_plot_val)
            stdev  = np.std(df.x_plot_val)
            
            stat_dict = {
                'nevt'    : nevt,
                'avg'     : avg,
                'median'  : median,
                'stdev'   : stdev
                # overflow/underflow ??
            }

            if do_nevt:   plt_legend += f'  N={int(nevt)}'
            if do_avg:    plt_legend += f'  avg={avg:.2f}'
            if do_stddev: plt_legend += f'  stdev={stdev:.2f}'  
            
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

            # XXXXXXXXX LEGACY XXXXXXXXXXXX                
            # determine plot style: 
            # default is solid-filled circles with error bars
            do_errorbar = True 
            do_ovsim    = False  # overlay sim
            chi2red     = 0.0

            if n > 0 and do_chi2:
                # prepare for sim overlay with histogram
                do_errorbar = False; do_ovsim = True 

            if do_errorbar :
                yval     = sb 
                yval_err = [sb-errl, erru-sb] 
                plt.errorbar(xcen, yval, yerr=yval_err, 
                             fmt=plt_marker, label=plt_legend )
            elif do_ovsim :
                x_val  = df.x_plot_val
                wgt_ov = [ scale ] * len(x_val)
                if weight:
                    wgt_user = get_weights_user(x_val,weight) 
                    wgt_ov   = np.multiply(wgt_ov,wgt_user)
                    
                plt.hist(x_val, bins, alpha=0.25, weights = wgt_ov,
                         label=plt_legend)
            else:
                sys.exit(f"\n ERROR: cannot determine which plot type: " \
                         f"errorbar or hist")

            for str_stat, val_stat in stat_dict.items():
                logging.info(f"\t {str_stat:8} value for {name_legend}:  {val_stat:.3f}")
            if do_list_cid:
                print_cid_list(df, name_legend)

        apply_plt_misc(args, plot_title, xlabel, ylabel, None)      

        # XXXXXXXXX LEGACY XXXXXXXXXXXX        
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

    elif DIFF :
        # 2D plot of difference between two files
        logging.info("plotting DIFF between two files.")
        try:         
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]])            
        except KeyError:     
            pass  # auto scale y axis
            
        keylist    = list(MASTER_DF_DICT.keys())   # tf0, tf1 ...
        df_ref_dict = MASTER_DF_DICT[keylist[0]]
        df_ref      = MASTER_DF_DICT[keylist[0]]['df']  # reference df  for difference
        for k in keylist[1:]:
            df_dict    = MASTER_DF_DICT[k] 
            df         = df_dict['df']            
            plt_alpha  = df_dict['alpha']
            plt_marker = df_dict['marker']
            
            if DIFF == ARG_DIFF_CID :
                #need to do an inner join with each entry in dic, then plot the diff
                # (join logic thanks to Charlie Prior)
                join = df_ref.join(df.set_index('CID'), on='CID', how='inner',
                                   lsuffix='_1', rsuffix='_2')
                plt.scatter(join.x_plot_val_1.values,
                            join.y_plot_val_1.values - join.y_plot_val_2.values,
                            alpha=plt_alpha, label='Diff')
                avgdiff = binned_statistic(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values, bins=bins, statistic='median')[0]                                     
                plt.scatter((bins[1:] + bins[:-1])/2, avgdiff, label="Mean Difference",
                            color='k', marker=plt_marker)
                
            elif DIFF == ARG_DIFF_ALL :
                text_label =  df_ref_dict['name_legend'] + " - " + df_dict['name_legend']
                try:
                    plt.scatter(df_ref.x_plot_val, df_ref.y_plot_val - df.y_plot_val, 
                                label=text_label, alpha=plt_alpha, marker=plt_marker)
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
                sys.exit(f"\n ERROR: Invalid DIFF option:  {DIFF} \n")

            # XXXXXXXXX LEGACY XXXXXXXXXXXX
            apply_plt_misc(args, plot_title, xlabel, ylabel+" diff", None)  
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
            plt_marker  = df_dict['marker']
            nevt        = len(df) 
            size        = 20 / math.log10(nevt)  # dot size gets smaller with nevt ??

            plt_legend = name_legend
            if do_nevt: plt_legend += f'  N={nevt}'

            plt.scatter(df.x_plot_val, df.y_plot_val, alpha=plt_alpha,
                        label=plt_legend, zorder=0, s=size, marker=plt_marker)

            # XXXXXXXXX LEGACY XXXXXXXXXXXX                
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

        apply_plt_misc(args, plot_title, xlabel, ylabel, None)  


    # XXXXXXXXX LEGACY XXXXXXXXXXXX        
    return 
    # end plotter_func_legacy

    
def apply_plt_misc(args, plot_title, xlabel, ylabel, plt_text_dict):

    # Created July 21 2024 by R.Kessler
    # wrapper for lots of plt. calls that are repeated several times.
    
    OPT          = args.OPT    
    do_logy      = OPT_LOGY    in OPT
    do_grid      = OPT_GRID    in OPT

    if do_logy:    plt.yscale("log")
    if do_grid:    plt.grid(zorder=3)

    plt.xlabel(xlabel)  
    if ylabel is not None: plt.ylabel(ylabel)
    
    plt.legend()                                
    plt.title(plot_title)

    if plt_text_dict is not None:
        x_text = plt_text_dict['x_text']
        y_text = plt_text_dict['y_text']
        text   = plt_text_dict['text']
        plt.text(x_text, y_text, text )
    
    return

def get_weights_user(xcen,weight):

    n_val     = len(xcen)

    if weight == '1' :
        wgt_vals  = [1.0] * n_val   # default weights are 1
    elif weight:

        # create numpy array x (from xcen) and x varname must be used to
        # match variable in user-define WEIGHT (e.g., 1-x).
        x   = np.array(xcen)
        weight_plus_np = weight        
        for str_orig, str_final in NUMPY_FUNC_DICT.items(): 
            weight_plus_np = weight_plus_np.replace(str_orig,str_final)

        wgt_vals = eval(weight_plus_np)
        
    return wgt_vals

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
                              
    # prepare input args
    arg_prep_driver(args)

    V = args.VARIABLE
    args.VARIABLE = translate_VARIABLE(V)   # add df. as needed to args.VARIABLE
    
    for icut, cut in enumerate(args.cut_list):
        args.cut_list[icut] = translate_CUT(cut)  # add df. and df.loc as needed
        
    if args.DEBUG_FLAG_DUMP_TRANSLATE :
        sys.exit(f"\n xxx bye .")
    
    plot_info = Namespace()  # someplace to store internally computed info
    
    set_var_dict(args, plot_info) # set plot bounds and axisinfo

    read_tables(args, plot_info)  # read each input file and store data frames
    #sys.exit(f"\n xxx DEBUG STOP xxx")
    
    plt.figure()  # initialize matplotlib figure

    if args.DEBUG_FLAG_LEGACY:
        plotter_func_legacy(args, plot_info)  # prepare plot in matplotlib
    else:
        plotter_func_driver(args, plot_info)  # prepare plot in matplotlib
        
    # check for user-define output (e.g. myplot.pdf, myplot.png, etc ...)
    # Note that we either save to file, or show in pop-up window ... but not both.
    # This enables creating plots as part of pipelines.
    SAVE = args.SAVE
    if SAVE is None :
        plt.show()  # show plot in pop-up window
    else:
        # save to file
        logging.info(f"Save plot to {SAVE}")
        fmt = SAVE.split(".")[-1]  # test.png -> fmt = png
        plt.savefig(SAVE, bbox_inches="tight", format = fmt)

    logging.info('Done.')
    
    # === END: ====
    

    
