#!/usr/bin/env python
#
# Created by B.Popovic during his graduate career at Duke University.
# Installed into SNANA Jun 24 2024 by R.Kessler
# Refactor to have __main__, and add translate_VARIABLE and tranlate_CUT
# methods to automatically append data commands to simplified user input.
# 
#
# ==============================================
import os, sys, gzip, copy, logging, math, re, gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import argparse
from argparse import RawTextHelpFormatter
from argparse import Namespace

parser=argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, prefix_chars='@')
from scipy.stats import binned_statistic
from collections import Counter
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

OPT_DIFF_CID  = "DIFF_CID"   # 2 files and 2D: plot y-axis diff for each CID
OPT_DIFF_ALL  = "DIFF_ALL"   # 2 files and 2D: plot y-axis diff between means
OPT_RATIO     = "RATIO"      # 1D: plot ratio between 2 files or 2 cuts

VALID_OPT_LIST = [ OPT_NEVT, OPT_AVG, OPT_MEAN, OPT_STDDEV, OPT_CHI2, OPT_CHI2,
                   OPT_MEDIAN, OPT_DIAG_LINE,
                   OPT_LOGY, OPT_GRID, OPT_LIST_CID,
                   OPT_DIFF_CID, OPT_DIFF_ALL, OPT_RATIO]

NMAX_CID_LIST = 20  # max number of CIDs to print for @@OPT CID_LIST

ARG_DIFF_CID = "CID"
ARG_DIFF_ALL = "ALL"
VALID_ARG_DIFF_LIST = [ ARG_DIFF_CID, ARG_DIFF_ALL ]

# define dictuionay of user-defined functions (key)
# and the np.xxx replacement for pandas. These functions
# can be used in variables (@V) and weigts (@@WEIGHT).
NUMPY_FUNC_DICT = {
    'exp'         :  'np.exp'  ,
    'log10'       :  'np.log'  ,   # works for log and log10        
    'log'         :  'np.log'  ,   # works for log and log10
    'sqrt'        :  'np.sqrt' ,
    'abs'         :  'np.abs'  ,
    'heaviside'   :  'np.heaviside'
}
NUMPY_FUNC_LIST = list(NUMPY_FUNC_DICT.keys())

# list possible VARNAME to identify row
VARNAME_IDROW_LIST = [ 'CID', 'GALID', 'ROW' ]

# internal strings to identify type of string in @V or @@CUT
STRTYPE_VAR   = "VARIABLE"
STRTYPE_NUM   = "NUMBER"
STRTYPE_DELIM = "DELIMETER"   # e.g., + - : / *
STRTYPE_FUNC  = "FUNCTION"    # e.g., log10, sqrt 
STRTYPE_LIST = [ STRTYPE_VAR, STRTYPE_NUM, STRTYPE_DELIM, STRTYPE_FUNC]

FIGSIZE = [ 6.4, 4.8 ]  # default figure size, inches

# internal DEBUG flags
DEBUG_FLAG_REFAC           =  2
DEBUG_FLAG_LEGACY          = -2
DEBUG_FLAG_DUMP_TRANSLATE  =  3
DEBUG_FLAG_DUMP_TRANSLATE2 =  33  # info dump for each char of @V and @@CUT string


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

  For 2D plots, code aborts if arg is missing colon;
  e.g., "@@ERROR MUERR" results in abort for 2D plot.

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

@@XLABEL 
   Override default x-axis label = variable name
@@YLABEL  
   Override default y-axis label = variable name

@@LEGEND
   Overwrite default legend; must give one arg per plot.
   For 2 files and @@DIF CID option, only need to give one @@LEGEND arg.

@@UNITS @U
  Default x- and y-axis labels show variable name only without units.
  To display units in (),
     @@UNITS days       # units for 1D plot; e.g., "PEAKMJD (days)"
                        #  (aborts for 2D plot because of missing colon)
     @@UNITS days:      # x-axis units for 2D plot (no y-axis units)
     @@UNITS :deg       # units for y-axis of 2D plot (no units for x-axis)
     @@UNITS days:deg   # units for both axes of 2D plot
       or
     @U <arg>

   For 2D plots, code aborts if @@UNITS arg is missing colon.

@@ALPHA
  Aadjust the matplotlib alpha values to adjust transparency 
  (0=transparent, 1=solid). For 2D overlays of multuple files or cuts, 
  can specify multiple alpha values, e.g.,
     @@ALPHA 0.9 0.1
  so that primary plot is dark and overlay plot is nearly transparent.

  For 2D plot with @@OPT MEDIAN or MEAN, setting ALPHA=0 will suppress 
  individual points and show only MEDIAN/MEAN.  

@@MARKER
  Override default solid circle markers with
     @@MARKER s ^  # square and triangle for 1st and 2nd file/cut
     @@MARKER s x  # square and X        for 1st and 2nd file/cut
     etc ... (see online doc on matplotlib markers)

@@XSIZE_SCALE
@@YSIZE_SCALE
  scale horizontal or vertical size of plot; e.g, 
    @@XSIZE_SCALE 0.5 -> make plot horizontallly narrower
    @@YSIZE_SCALE 0.5 -> make plot vertically shorter
    
 
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
      {OPT_STDDEV:<12} ==> 1D -> append stddev on each legend;
                       2D -> mean error bar is STDDEV instead of STDDEV/sqrt(N)
      {OPT_CHI2:<12} ==> display chi2/dof on plot for two table files (1D only).
      {OPT_DIAG_LINE:<12} ==> draw line with slope=1 for 2D plot.
      {OPT_LOGY:<12} ==> log scale for vertical axis.
      {OPT_GRID:<12} ==> display grid on plot.
      {OPT_LIST_CID:<12} ==> print up to 100 CIDs passing cuts.      
      {OPT_DIFF_ALL:<12} ==> 2D: plot y-axis difference in median values between 2 plots
                       (can compare correlated or independent samples)
      {OPT_DIFF_CID:<12} ==> 2D: plot y-axis difference for each CID.
                       (compare correlated samples only; must have CID overlap)
      {OPT_RATIO:<12} ==> 1D: plot historgram ratio file0/file1 or cut0/cut1
                       (binomial error if 1st file(cut) is subset of 2nd file)


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
    parser.add_argument("@@NROWS", "@@nrows", help=msg, type=int, default=None)
    
    # -------- decorative options below ---------
    
    msg = "AUTO (default): optional bounds min, max and binsize of plot parameters.\n" \
          "For 2D plot, must specify both x and y bounds. y-binsize is ignored."
    parser.add_argument('@@BOUNDS', '@@bounds', default=BOUNDS_AUTO, help=msg, nargs = '+')

    msg = "Units to show for each axis; see @@HELP"
    parser.add_argument('@U', '@@UNITS', '@@units', default=None, help=msg, nargs="+")

    msg = "x-axis label (override default)"
    parser.add_argument('@@XLABEL', '@@xlabel', default=None, help=msg, type=str)
    msg = "y-axis label (override default)"
    parser.add_argument('@@YLABEL', '@@ylabel', default=None, help=msg, type=str)        
    
    msg = "Override default plot title"
    parser.add_argument("@@TITLE", "@@title", default=None, help=msg )

    msg = "Override default legend on plot (space sep list per TFILE)"
    parser.add_argument('@@LEGEND', '@@legend', default=None, help=msg, nargs="+")

    msg = "Override default marker='o'"
    parser.add_argument('@@MARKER', '@@marker', default=['o'], help=msg, nargs="+")    

    msg = "scale horizontal size of plot"
    parser.add_argument('@@XSIZE_SCALE', '@@xsize_scale',
                        default=1.0, help=msg, type=float)
    msg = "scale vertical size of plot"
    parser.add_argument('@@YSIZE_SCALE', '@@ysize_scale',
                        default=1.0, help=msg, type=float)    
    
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

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx args.VARIABLE = {args.VARIABLE}  (before modifications)")  
        print(f" xxx args.CUT      = {args.CUT}  (before modifications)")
        
    # - - - - - - -
    # for variable(s), remove pad spacin
    args.VARIABLE_ORIG = args.VARIABLE 
    if args.VARIABLE:
        args.VARIABLE      = ''.join([str(elem) for elem in args.VARIABLE])        
        args.VARIABLE_ORIG = args.VARIABLE_ORIG[0]
    else:
        sys.exit(f"\n ERROR: must define variable(s) to plot with @@VARIABLE or @@V")

    # CUT is tricky. Make sure that length of cut list matchs length
    # of table-file list ... or vice versa ... make sure that length of
    # tfile list matches length of CUT list.

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

    args.use_err_list = True  # default
    if args.error:
        args.ERROR = args.error
        args.use_err_list = False  # disable error bar on 2D plot

    ndim       = len(args.VARIABLE_ORIG.split(COLON))
    args.UNITS = arg_prep_axis(ndim, 'UNITS', args.UNITS)
    args.ERROR = arg_prep_axis(ndim, 'ERROR', args.ERROR)

    
    if args.BOUNDS != BOUNDS_AUTO:
        args.BOUNDS = ' '.join([str(elem) for elem in args.BOUNDS])
    
    # if only 1 alpha, make sure there is alpha for each file/cut
    args.alpha_list  = arg_prep_extend_list(narg_tfile, args.ALPHA)
    args.marker_list = arg_prep_extend_list(narg_tfile, args.MARKER) 

    args = args_prep_DIFF(args)    

    args.legend_list = arg_prep_legend(args) # must be after prep_DIFF
    
    # - - - - - - -
    args.OPT = arg_prep_OPT(args)            
    
    return  # end arg_prep_driver


def get_list_match(list0, list1):    
    # if any element of opt_subset_list is in OPT_LIST,
    # return that element; else return None

    if list0 is None: return None
    if list1 is None: return None
    
    for opt in list0:
        if opt in list1:
            return opt

    # if we get here, return None
    return None

def args_prep_DIFF(args):

    # set args.DIFF internally, and append args.OPT with MEDIAN
    # if no stat option is given.
    
    if args.DIFF:
        sys.exit(f"\n ERROR: @@DIFF {args.DIFF} is obsolete and no longer valid;\n" \
                 f"\t instead, use @@OPT DIFF_{args.DIFF}" )
    
    # if user has not specified any kind of statistical average
    # for DIFF, then set default MEDIAN
    OPT        = args.OPT
    opt_diff   = get_list_match( [OPT_DIFF_CID, OPT_DIFF_ALL], args.OPT)
    opt_stat   = get_list_match( [OPT_MEDIAN,OPT_MEAN,OPT_AVG], args.OPT)
    args.DIFF  = opt_diff  # e.g., 'CID' or 'ALL' or None
        
    if opt_diff and  opt_stat is None :
        args.OPT.append(OPT_MEDIAN) 

    return args

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

    args.DEBUG_FLAG_REFAC           = True
    args.DEBUG_FLAG_LEGACY          = False
    args.DEBUG_FLAG_DUMP_TRANSLATE  = False
    args.DEBUG_FLAG_DUMP_TRANSLATE2 = False    
    if args.DEBUG_FLAG != 0:
        args.DEBUG_FLAG_LEGACY          = args.DEBUG_FLAG == DEBUG_FLAG_LEGACY
        args.DEBUG_FLAG_REFAC           = not args.DEBUG_FLAG_LEGACY
        
        args.DEBUG_FLAG_DUMP_TRANSLATE2 = args.DEBUG_FLAG == DEBUG_FLAG_DUMP_TRANSLATE2        
        args.DEBUG_FLAG_DUMP_TRANSLATE  = \
            args.DEBUG_FLAG == DEBUG_FLAG_DUMP_TRANSLATE   or \
            args.DEBUG_FLAG == DEBUG_FLAG_DUMP_TRANSLATE2
        
        logging.info(f"# \t DEBUG_FLAG={args.DEBUG_FLAG} " )
        logging.info(f"# \t DEBUG_FLAG_REFAC          = {args.DEBUG_FLAG_REFAC}")
        logging.info(f"# \t DEBUG_FLAG_DUMP_TRANSLATE = {args.DEBUG_FLAG_DUMP_TRANSLATE}")
        print(f"")        
    
    return

def arg_prep_legend(args):

    LEGEND_orig = args.LEGEND
    LEGEND_out  = LEGEND_orig   # default is user input

    do_diff  = args.DIFF
    do_ratio = OPT_RATIO in args.OPT
    narg_cut = len(args.CUT)
    
    if LEGEND_orig is None:
        # no user supplied legend, so make up reasonable legend
        LEGEND_out = [ None ] * len(args.TFILE)
        if narg_cut > 1 :
            LEGEND_out = args.CUT
            #if do_ratio:
            #    legend_ratio = args.CUT[0] + ' / ' + args.CUT[1]
            #    LEGEND_out   = [legend_ratio] * narg_cut
        else :
            LEGEND_out = []
            for t in args.tfile_base_list:
                legend = t.split('.')[0]  # base file name without extension after dot
                LEGEND_out.append(legend)


    # - - - - - 
    narg_legend  = len(LEGEND_out)
    narg_tfile   = len(args.tfile_list)

    allow_single_legend = args.DIFF or OPT_RATIO in args.OPT
    if allow_single_legend and narg_legend == 1:
        LEGEND_out.append(LEGEND_out[0])
        narg_legend  = len(LEGEND_out)

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

    
def arg_prep_axis(ndim,name_arg, arg_axis):
    
    # If arg_axis is None, return ':'
    # If ndim==2 and there is no colon, abort.
    if arg_axis:
        arg_axis  = ''.join([str(elem) for elem in arg_axis])
        if ndim == 2 and COLON not in arg_axis:
            sys.exit(f"\n ERROR: ndim={ndim} but {name_arg} = {arg_axis} has no colon.")
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

def split_var_string(STRING, DELIM_LIST, FUNC_LIST):

    # Created Aug 21 2024
    # This is a custom split that splits input string into list of
    #  * variables names
    #  * numbers
    #  * delimeters (from input DELIM_LIST)
    #  * functions  (from input FUNC_LIST)
    #
    # Example:
    #  STRING = 'log10(zcmb):SNRMAX'
    #   returns ['log10', '(',     'zcmb',    ')',    ':',   'SNRMAX']
    #   and     ['FUNC' , 'DELIM', 'VAR',  'DELIM', 'DELIM', 'VAR' ]
    #
    #  STRING = 'zcmb:z-zcmb + .003'
    #    return ['zcmb', ':', 'z', '-', 'zcmb', '+', '.003']
    #
    split_list   = []
    strtype_list = [] 

    STRING_local = STRING
    
    # remove 'np.' that isn't needed because later code adds it.
    if STR_np in STRING:        STRING_local = STRING.replace(STR_np,'')
    len_str = len(STRING_local)
    
    if args.DEBUG_FLAG_DUMP_TRANSLATE2:
        print(f"\n split_var_string: prepare dump for each char in \n" \
              f"\t STRING_local({len_str}) = {STRING_local}")
        
            
    j_last = 0
    for j in range(0,len_str):
        ch         = STRING_local[j:j+1]  # current char
        is_last    = (j == len_str-1)
        is_delim   = (ch in DELIM_LIST)
        str_last = None
        j_last_dump = j_last
        if is_delim:
            str_last = STRING_local[j_last:j].replace(' ','')
            j_last = j+1
        elif is_last:
            str_last = STRING_local[j_last:j+1].replace(' ','')
            
        valid_str_last = (str_last is not None and len(str_last)>0)

        if valid_str_last:
            split_list.append(str_last)
            if is_number(str_last) or "'" in str_last:
                strtype_list.append(STRTYPE_NUM)
            elif str_last in FUNC_LIST:
                strtype_list.append(STRTYPE_FUNC)
            else:
                strtype_list.append(STRTYPE_VAR)

                
        if is_delim:
            split_list.append(ch)
            strtype_list.append(STRTYPE_DELIM)

        if args.DEBUG_FLAG_DUMP_TRANSLATE2 :
            print(f"\t xxx - - - - - - - - ")
            print(f"\t xxx j={j}  j_last={j_last_dump}  ch={ch}  " \
                  f"is_[last,delim]={is_last},{is_delim}  " \
                  f"str_last={str_last}  valid_str_last={valid_str_last}")
            print(f"\t xxx split_list -> {split_list}")
        
    # - - - -
    return split_list, strtype_list
    # end split_var_string
    

def get_var_list_legacy(VARIABLE, DELIMITER_LIST):
    # if VARIABLE = 'zHD-zHD_2:SNRMAX' -> return var_list = ['zHD', 'zHD_2', 'SNRMAX']
    # which is a list of variables wihtout symbols

    # first replace any valid delimiter with '!' so that we can split on
    # single ! char
    VAR_TMP = copy.copy(VARIABLE)

    VAR_TMP = VAR_TMP.strip()  # remove pad spaces

    # XXXXXXXXX LEGACY/OBSOLETE XXXXXXXXXXXXXX
    
    # replace algabraic delimiters with pad space for easier split
    for delim in DELIMITER_LIST:
        VAR_TMP = VAR_TMP.replace(delim,' ')

    var_list_all = sorted(VAR_TMP.split())

    # keep elements that do NOT have numpy "np."
    var_list = []
    for var in var_list_all:
        if STR_np not in var:
            var_list.append(var)

    # XXXXXXXXX LEGACY/OBSOLETE XXXXXXXXXXXXXX
    
    # create supplemental list of logicals indicating which
    # var_list elements are numbers. A string in quotes is
    # treated like a number.
    isnum_list = []
    for var in var_list:
        isnum = is_number(var) or "'" in var
        isnum_list.append(isnum)

    # XXXXXXXXX LEGACY/OBSOLETE XXXXXXXXXXXXXX        
    return  var_list, isnum_list
    # end get_var_list_legacy
    
def translate_VARIABLE(VARIABLE):
    # add df. as needed to VARIABLE
    # assume that all df. have been provided by user, or none;
    # will not handle in-between cases.

    VARIABLE_ORIG = VARIABLE
    
    if STR_df in VARIABLE:
        return VARIABLE

    # break down VARIABLE into a list of variables without any symbols

    split_list, strtype_list = \
        split_var_string(VARIABLE, DELIMITER_VAR_LIST, NUMPY_FUNC_LIST)

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx User input VARIABLE = {VARIABLE}")
        for tmp_str, tmp_type in zip(split_list,strtype_list):
            tmp = f"'{tmp_str}'"
            print(f"\t xxx var sub-string = {tmp:<14} is {tmp_type}")

        
    # next, put df. in front of each variable, and np. in front of functions
    VARIABLE_df = ''
    for tmp_str, tmp_type in zip(split_list, strtype_list):
        tmp_append = tmp_str  # default string to append
        if tmp_type == STRTYPE_VAR:
            tmp_append = STR_df + tmp_str
        elif tmp_type == STRTYPE_FUNC :
            tmp_append = STR_np + tmp_str
        else:
            pass
        
        VARIABLE_df += tmp_append    
        
    logging.info(f"Translate VARIABLE {VARIABLE_ORIG}")
    logging.info(f"             ->    {VARIABLE_df}")
    
    return VARIABLE_df
    # end translate_VARIABLE
    
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
        split_list, strtype_list = \
            split_var_string(cut, DELIMITER_CUT_LIST, NUMPY_FUNC_LIST)

        cut_df = ''
        for tmp_str, tmp_type in zip(split_list, strtype_list ):
            if args.DEBUG_FLAG_DUMP_TRANSLATE:
                tmp = f"'{tmp_str}'"
                print(f"\t xxx cut sub-string = {tmp:<14} is {tmp_type}")
                
            tmp_append = tmp_str
            if tmp_type == STRTYPE_VAR:
                tmp_append = STR_df + tmp_str
            elif tmp_type == STRTYPE_FUNC:
                tmp_append = STR_np + tmp_str
                sys.exit(f"\n ERROR: cannot use numpy functons in @@CUT; " \
                         f"remove '{tmp_str}' from @@CUT arg")
            else:
                pass
            cut_df += tmp_append
            

        if STR_np not in cut_df:
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

    logging.info(f"Translate CUT {CUT}  ")
    logging.info(f"          ->  {CUT_df}")
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

    all_var_list = []  # list of all vars and varerr used to check existence
    
    for n, VAR in enumerate(VAR_LIST):

        all_var_list.append(VAR)
        
        STR_VAR      = str(VAR)
        STR_LABEL    = STR_VAR.replace(STR_df,'') 
        if STR_np in STR_LABEL:
            STR_LABEL = STR_LABEL.replace(STR_np,'')

        VARERR  = VARERR_LIST[n]
        if len(VARERR) > 0:
            all_var_list.append(VARERR)
            STR_VARERR   = STR_df + str(VARERR)
        else:
            STR_VARERR   = None            
        
        U = UNITS_LIST[n]
        if len(U) > 0:  STR_LABEL += '  (' + U + ')'
            
        if n == 0:
            xlabel = STR_LABEL
            ylabel = 'Counts'
            if args.WEIGHT[0]        : ylabel = f'Counts * {args.WEIGHT}'
            if OPT_RATIO in args.OPT : ylabel = 'Ratio'
            if args.XLABEL   : xlabel = args.XLABEL # user override
            if args.YLABEL   : ylabel = args.YLABEL # user override
            plotdic['x']            = STR_VAR
            plotdic['xerr']         = STR_VARERR
            plotdic['xaxis_label']  = xlabel
            plotdic['yaxis_label']  = ylabel
            plotdic['ndim']         = 1
            
        else:
            ylabel = STR_LABEL
            if args.DIFF   : ylabel += ' Diff'            
            if args.YLABEL : ylabel = args.YLABEL
            plotdic['y']            = STR_VAR
            plotdic['yerr']         = STR_VARERR
            plotdic['yaxis_label']  = ylabel
            plotdic['ndim']         = 2

    plotdic['set_ylim'] = plotdic['ndim'] == 2 or OPT_RATIO in args.OPT
    
    # check for user (custom) bounds in plot
    custom_bounds = False
    boundsdic     = {}
    if BOUNDS != BOUNDS_AUTO:
        for n,BND in enumerate(BOUNDS.split(':')):
            custom_bounds = True
            load_xbound = (n==0)
            load_ybound = (n >0)
            if load_xbound:
                boundsdic['x'] = [float(i) for i in BND.split()]
            if load_ybound:
                boundsdic['y'] = [float(i) for i in BND.split()]

    # - - - - - -  - - 
    # error checks 
    if args.DIFF:
        if plotdic['ndim'] == 1 :
            sys.exit("\nERROR: DIFF does not work for 1D histograms." \
                     " ABORT to avoid confusion.")

        if len(args.TFILE) == 1:
            sys.exit(f"\n ERROR '@@DIFF {args.DIFF}' does not work with 1 table file;\n"\
                     f"\t need 2 or more table files specified after @@TFILE key.")
            
    # load output namespace
    plot_info.plotdic             = plotdic
    plot_info.plot_title          = plot_title
    plot_info.boundsdic           = boundsdic
    plot_info.custom_bounds       = custom_bounds
    plot_info.all_var_list        = all_var_list
    
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

        # check for DOCANA block and count number of lines to skip before VARNAMES;
        # e.g., skip DOCANA for HOSTLIB 
        nrow_skip = count_rows_to_skip(tfile)
        
        df  = pd.read_csv(tfile, comment="#", sep=r"\s+",
                          skiprows=nrow_skip, nrows=NROWS)

        # xxxxxxx mark delete Aug 23 2024 xxxxxxx
        #if NROWS > 0 :
        #    # read NROWS subset
        #    df  = pd.read_csv(tfile, comment="#", sep=r"\s+", nrows=NROWS)
        #else:
        #    # read all
        #    df  = pd.read_csv(tfile, comment="#", sep=r"\s+")
        # xxxxxxxxxxxxxxxx
        
        # figure out KEYNAME to id events
        plot_info.varname_idrow = None
        for varname_idrow in VARNAME_IDROW_LIST:
            if varname_idrow in df:
                plot_info.varname_idrow            = varname_idrow
                
        varname_idrow = plot_info.varname_idrow
        if varname_idrow is None:
            sys.exit(f"\n ERROR: could not find valid VARNAME_IDROW " \
                     f"among {VARNAME_IDROW_LIST}")

        # make sure that all requested variables (and error-vars) exist in df
        # xxx ?? check_vars_exist(args, df, plot_info)
        
        try:
            df[varname_idrow] = df[varname_idrow].astype(str)
        except KeyError:
            logging.warn(f"No {varname_idrow} present in this file. OK for some file types.")

        # apply user cuts
        if cut:
            df = eval(cut)

        # drop duplicates for DIFF CID matching
        if args.DIFF == ARG_DIFF_CID:
            df.drop_duplicates(varname_idrow, inplace=True)

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

def count_rows_to_skip(tfile):

    nrow_skip = 0
    if '.gz' in tfile:
        t = gzip.open(tfile,'rt')
    else:
        t = open(tfile,'rt')

    nrow_skip = 0
    
    for line in t:
        wdlist = line.split()
        if len(wdlist) == 0 :
            nrow_skip += 1
            continue

        if wdlist[0] == 'VARNAMES:' :
            break
        else:
            nrow_skip += 1            
            continue

    t.close()
    logging.info(f"\t found {nrow_skip} rows to skip before 'VARNAMES:' key")

    return nrow_skip
    # end count_rows_to_skip

def check_vars_exist(args, df, plot_info):

    # Tricky because all_var_list can include math symbols.
    varlist_missing = []
    for tmp_var_df in plot_info.all_var_list:
        tmp_var = tmp_var_df.split(STR_df)[1]  # remove df.
        if tmp_var not in df:
            varlist_missing.append(tmp_var)
            logging.warning(f"Could not find requested {tmp_var} for {tfile}")
    assert (len(varlist_missing) == 0),  \
        f"\n ERROR: Missing {varlist_missing} in df; " + \
        f"check @@VARIABLE and  @@ERROR args" 
    return

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
    varname_idrow        = plot_info.varname_idrow
    
    xlabel               = plotdic['xaxis_label']
    ylabel               = plotdic['yaxis_label']
    NDIM_PLOT            = plotdic['ndim']  # 1 or 2
    set_ylim             = plotdic['set_ylim']
    
    xbins = None
    if custom_bounds:                       
        xmin  = boundsdic['x'][0]; xmax = boundsdic['x'][1]
        xbin  = boundsdic['x'][2]
        xbins = np.arange(xmin, xmax, xbin)
        if set_ylim:
            ymin = boundsdic['y'][0]; ymax = boundsdic['y'][1]
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
        info_plot_dict['xbins_cen']           =  xbins_cen        
        info_plot_dict['df_dict']             =  df_dict        
        info_plot_dict['varname_idrow']       =  varname_idrow
        
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
            if set_ylim: plt.ylim(ymin, ymax)
        
        if do_plot_errorbar :
            # nominal
            if plt_alpha > 0:
                plt.errorbar(xval_list, yval_list,
                             xerr=xerr_list, yerr=yerr_list, 
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
            plt.plot(x,y, zorder=10)

        apply_plt_misc(args, plot_title, xlabel, ylabel, plt_text_dict)
        numplot += 1
        
        
    return   # end plotter_func_driver


def get_info_plot1d(args, info_plot_dict):

    do_chi2   = OPT_CHI2   in args.OPT
    do_nevt   = OPT_NEVT   in args.OPT
    do_avg    = OPT_AVG    in args.OPT  or  OPT_MEAN in args.OPT
    do_ratio  = OPT_RATIO  in args.OPT
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
    lo_list, up_list = poisson_interval(yval_list) # compute poisson interval
    errl_list        = yval_list - lo_list         # low-side error bar
    erru_list        = up_list   - yval_list       # high-side error bar
    
    wgt_ov = None
    ov_scale = 1.0  # default overlay plot is not scaled
    
    if numplot == 0 :
        info_plot_dict['df0']          = df
        info_plot_dict['yval0_list']   = yval_list # save for norm overlay or ratio
        info_plot_dict['errl0_list']   = errl_list
        info_plot_dict['erru0_list']   = erru_list
        info_plot_dict['name0_legend'] = name_legend        
        if do_ratio:
            do_plot_errorbar = False
    else:
        df0          = info_plot_dict['df0'] 
        yval0_list   = info_plot_dict['yval0_list']
        errl0_list   = info_plot_dict['errl0_list']
        erru0_list   = info_plot_dict['erru0_list']        
        name0_legend = info_plot_dict['name0_legend']         
        if do_chi2:
            # re-scale overlay plot only if chi2 option is requested
            ov_scale = np.sum(yval0_list) / np.sum(yval_list)
            do_plot_errorbar = False
            do_ov1d_hist     = True
        elif do_ratio:
            # check of df0 is a subset of df --> binomial errors
            is_subset = is_cid_subset(info_plot_dict, df0, df)
            ratio, errl_ratio, erru_ratio = \
                compute_ratio(yval0_list, errl0_list, erru0_list,
                              yval_list,  errl_list,  erru_list, is_subset  )        
            # load arrays that will be store in dictionary below              
            yval_list  = ratio
            errl_list  = errl_ratio
            erru_list  = erru_ratio
            ratio = len(df0) / len(df)
            plt_legend = f"<ratio> = {ratio:.2f}"
        else:
            pass

        yval_list  *= ov_scale # normalize integral to match file 0 
        errl_list  *= ov_scale
        erru_list  *= ov_scale
        logging.info(f"Overlay {name_legend} scaled by {ov_scale:.3e}")

    # ---------------

    # check for overlay weights with chi2
    if numplot > 0 and do_chi2:
        wgt_ov = [ ov_scale ] * len(df.x_plot_val)
        if weight :
            wgt_user = get_weights_user(df.x_plot_val,weight) 
            wgt_ov   = np.multiply(wgt_ov,wgt_user)
            len_wgt_ov = len(wgt_ov)

        # check option to compute and print chi2/dof info on plot
        # Froce min error =1 in chi2 calc so that it's ok to plot
        # error=0 for bins with zero events.
        sqdif = (yval0_list-yval_list)**2
        sqerr = np.maximum((yval0_list + yval_list*ov_scale*ov_scale), 1.0)
        chi2  = np.sum( sqdif / sqerr )
        ndof  = len(xbins) - 1 
        text_chi2 = f"chi2/dof = {chi2:.1f}/{ndof}"
        logging.info(f"{name0_legend}/{name_legend} {text_chi2}")
        xmin = xbins[0];  xmax = xbins[-1]
        x_text = xmin + 0.7*(xmax-xmin)
        y_text = 1.0 * np.max(yval0_list)  # warning; fragile coord calc
        plt_text_dict = { 'x_text': x_text, 'y_text': y_text, 'text': text_chi2 }


    # - - - - - - -
    yerr_list = [errl_list, erru_list]
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
    info_plot_dict['xval_list']         = xval_list
    info_plot_dict['yval_list']         = yval_list
    info_plot_dict['xerr_list']         = None
    info_plot_dict['yerr_list']         = yerr_list
    info_plot_dict['plt_size']          = None
    info_plot_dict['plt_legend']        = plt_legend
    info_plot_dict['plt_text_dict']     = plt_text_dict
    info_plot_dict['wgt_ov']            = wgt_ov
     
    return  # end get_info_plot1d


def is_cid_subset(plot_dict, df0, df1):
    
    # return True if cids (or galid or rownum) in df0 are a subset of df1;
    # else return False. Rather than require a strict overlap, 2% of the
    # df0 elements are allowed to be outside of df1, and still return True.
    # This handles cases where a naively expected subset has a few outliers;
    # e.g, light curve fitting with true redshift vs. photo-z, we naively
    # expect the zphot sample to be a subset of ztrue set ... but in reality
    # there can be a very small fraction of zphot fits that fail ztrue fit,
    # and thus zphot is not a literal subset of ztrue. For plotting zphot/ztrue
    # efficiency, however, binomial errors are appropriate.
    #
    # This result is used to determine if errors on a ratio are independent
    # or correlated-binomial.

    MIN_OVERLAP_SUBSET = 0.98
    is_subset     = False
    varname_idrow = plot_dict['varname_idrow']

    n0 = len(df0)
    n1 = len(df1)
    temp = Counter(df0[varname_idrow]) - Counter(df1[varname_idrow])
    ndif = len(temp)
    overlap = 1.0 - ndif/n0

    is_subset = ( overlap > MIN_OVERLAP_SUBSET)
    
    # xxx mark delete
    #id_set0       = set(df0[varname_idrow])
    #id_set1       = set(df1[varname_idrow])
    #is_subset     = id_set0.issubset(id_set1)
    # xxx end
    
    return is_subset

def compute_ratio(yval0_list, errl0_list, erru0_list,
                  yval1_list, errl1_list, erru1_list, is_subset ):

    # for two sets of input y-axis values and errors (lo and upper),
    # return ratio and error (lo,up) on ratio.
    # If is_subset=True, compute binomial error instead of default
    # independent error.

    ratio  = yval0_list / yval1_list
        
    if is_subset:
        # binomial error with set0 being subset of set1,
        # and input errors are ignored;
        # COV = p ( 1 − p ) * NTOT (per bin)
        logging.info("\t Compute binomial error for ratio.")
        p          = yval0_list/yval1_list
        cov_ratio  = p*(1-p) / yval1_list
        errl_ratio = np.sqrt(cov_ratio)
        erru_ratio = errl_ratio
        
    else:
        # independent error
        logging.info("\t Compute independent error for ratio.")
        sqrelerrl  = (errl0_list/yval0_list)**2 + (errl1_list/yval1_list)**2
        sqrelerru  = (erru0_list/yval0_list)**2 + (erru1_list/yval1_list)**2
        errl_ratio = ratio * np.sqrt(sqrelerrl)
        erru_ratio = ratio * np.sqrt(sqrelerru)
    
    return ratio, errl_ratio, erru_ratio

    # end compute_ratio
    
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
            varname_idrow = plot_info.varname_idrow # e.g., CID or GALID or ROW
            join   = df_ref.join(df.set_index(varname_idrow), on=varname_idrow,
                                 how='inner', lsuffix='_0', rsuffix='_1')
            xval_list = join.x_plot_val_0.values
            yval_list = join.y_plot_val_0.values - join.y_plot_val_1.values
            info_plot_dict['xval_list']     = xval_list
            info_plot_dict['yval_list']     = yval_list
            if yerr_list is not None:
                yerr_list_0 = join.y_plot_err_0.values
                info_plot_dict['yerr_list']  = yerr_list_0
                
        elif args.DIFF == ARG_DIFF_ALL :

            if OPT_MEDIAN in args.OPT:
                which_stat = 'median'
            elif OPT_MEAN in args.OPT:
                which_stat = 'mean'                

            xbins  = info_plot_dict['xbins']
            y_ref  = binned_statistic(df_ref.x_plot_val, df_ref.y_plot_val,
                                      bins=xbins, statistic=which_stat)[0]
            y      = binned_statistic(df.x_plot_val, df.y_plot_val,
                                      bins=xbins, statistic=which_stat)[0]            

            # convert NaN (empty bins) to zero to avoid crash when plotting
            y_ref[np.isnan(y_ref)] = 0
            y[np.isnan(y)] = 0                        
            #sys.exit(f"\n xxx y_ref = \n{y_ref} \n xxx y = \n{y}")
            info_plot_dict['xval_list'] = info_plot_dict['xbins_cen']
            info_plot_dict['yval_list'] = y_ref - y
            info_plot_dict['yerr_list'] = None            
        
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
    plt_marker   = df_dict['marker']
    
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

        # compute error on the mean per bin: stddev/sqrt(N), or just STDDEV
        y_std   = binned_statistic(xval_list, yval_list,
                                   bins=xbins, statistic='std')[0]
        y_count = binned_statistic(xval_list, yval_list,
                                   bins=xbins, statistic='count')[0]
        
        if OPT_STDDEV in OPT:
            y_err_stat = y_std  # stddev per bin (user input @@OPT STTDEV)
        else:
            y_err_stat = y_std/np.sqrt(y_count)  # error on mean per bin (default)
        
        if do_wgtavg:
            stat_legend = 'wgtavg'
            y_avg   = copy.deepcopy(y_stat)  # for diagnostic print
               
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

        # - - - - - - -
        plt_legend  += ' ' + stat_legend
        plt.errorbar(xbins_cen, y_stat, yerr=y_err_stat,
                     fmt=plt_marker, label=plt_legend, zorder=5 ) 
            
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

    fsize = 12
    plt.xlabel(xlabel, fontsize=fsize)  
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=fsize)
    
    plt.legend()                                
    plt.title(plot_title, fontsize=14)

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
    #varname_idrow = plot_info.varname_idrow # e.g., CID or GALID or ROW 

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
    logging.info('')
    
    for icut, cut in enumerate(args.cut_list):
        args.cut_list[icut] = translate_CUT(cut)  # add df. and df.loc as needed
    logging.info('')
    
    if args.DEBUG_FLAG_DUMP_TRANSLATE :
        sys.exit(f"\n xxx bye .")
    
    plot_info = Namespace()  # someplace to store internally computed info
    
    set_var_dict(args, plot_info) # set plot bounds and axisinfo

    read_tables(args, plot_info)  # read each input file and store data frames
    #sys.exit(f"\n xxx DEBUG STOP xxx")

    # initialize matplotlib figure
    plt.figure( figsize=(FIGSIZE[0]*args.XSIZE_SCALE,
                         FIGSIZE[1]*args.YSIZE_SCALE) ) 

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
    

    
