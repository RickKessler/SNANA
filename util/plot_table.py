#!/usr/bin/env python
#
# Created by B.Popovic during his graduate career at Duke University.
# Installed into SNANA Jun 24 2024 by R.Kessler
# Refactor to have __main__, and add translate_VARIABLE and tranlate_CUT
# methods to automatically append df, df.loc and np to input variables
# to simplified user input.
#
# Jan 7 2025: input "@@LEGEND NONE" suppresses legend.
#
# ==============================================
import os, sys, gzip, copy, logging, math, re, gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import colors

import argparse
from argparse import RawTextHelpFormatter
from argparse import Namespace

parser=argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, prefix_chars='@')
from scipy.stats import binned_statistic
from collections import Counter
#import distutils.util 

# ====================
# globals

NXBIN_AUTO  = 30        # number of x-bins of user does not provide @@BOUNDS arg
NXYBIN_AUTO = 20        # auto nbin for hist2d (if @@BOUNDS is not given)

DELIMITER_VAR_LIST  = [ '+', '-', '/', '*', ':', '(', ')' ]  # for @@VARIABLE 
DELIMITER_CUT_LIST  = [ '&', '|', '>=', '<=', '>', '<',
                        '==', '!=', '=', '*', '+', '-', '/' ]  # for @@CUT

COLON = ':'
BAR   = '|'

DEFAULT  = 'DEFAULT'
SUPPRESS = 'SUPPRESS'

# define pandas strings to add into VARIABLE and CUT strings.
STR_df         = 'df.'
STR_df_loc     = 'df.loc'
STR_df_iloc    = 'df.iloc'
STR_np         = 'np.'

# define strings for @@OPT
OPT_HIST      = "HIST"      # 1D or 2D histogram instead of points with error bars
OPT_HISTFILL  = "HISTFILL"  # same, but stepfilled 1Dhist
OPT_NEVT      = "NEVT"      # append N={nevt} to legend
OPT_AVG       = "AVG"       # 1D->append avg to legend; 2D->overlay avg in x-bins
OPT_MEAN      = "MEAN"      # same as AVG
OPT_STDDEV    = "STDDEV"    # append stddev to legend
OPT_OV        = "OV"        # show tfile1/tfile2 chi2/dof and scale tfile2 to match tfile1;
OPT_OVCHI2    = "OVCHI2"    # exec OPT_OV and also print CHI2/dof on plot
OPT_CHI2      = "CHI2"      # legacy option; same as OVCHI2
OPT_MEDIAN    = "MEDIAN"
OPT_DIAG_LINE = "DIAG_LINE"  # draw diagonal line on plot
OPT_LOGY      = "LOGY"       # log scale along Y axis (1D or 2D)
OPT_LOGZ      = "LOGZ"       # log scale along Z axis (2D only)
OPT_GRID      = "GRID"       # draw grid on plot
OPT_LIST_CID  = "LIST_CID"   # list CIDs passing cuts
OPT_LIST_ROW  = "LIST_ROW"   # list ROWS passing cuts (same as LIST_CID)

OPT_DIFF_CID  = "DIFF_CID"   # 2 files and 2D: plot y-axis diff for each CID
OPT_DIFF_ALL  = "DIFF_ALL"   # 2 files and 2D: plot y-axis diff between means
OPT_RATIO     = "RATIO"      # 1D: plot ratio between 2 files or 2 cuts

VALID_OPT_LIST = [ OPT_HIST, OPT_HISTFILL,
                   OPT_NEVT, OPT_AVG, OPT_MEAN, OPT_STDDEV, OPT_OV, OPT_OVCHI2, OPT_CHI2,
                   OPT_MEDIAN, OPT_DIAG_LINE,
                   OPT_LOGY, OPT_LOGZ, OPT_GRID, OPT_LIST_CID, OPT_LIST_ROW,
                   OPT_DIFF_CID, OPT_DIFF_ALL, OPT_RATIO]

NMAX_CID_LIST = 20  # max number of CIDs to print for @@OPT CID_LIST


# define dictuionay of user-defined functions (key)
# and the np.xxx replacement for pandas. These functions
# can be used in variables (@V) and weigts (@@WGTFUN).
NUMPY_FUNC_DICT = {
    'exp'         :  'np.exp'  ,
    'log10'       :  'np.log'  ,   # works for log and log10        
    'log'         :  'np.log'  ,   # works for log and log10
    'sqrt'        :  'np.sqrt' ,
    'abs'         :  'np.abs'  ,
    'cos'         :  'np.cos'  ,
    'sin'         :  'np.sin'  ,
    'tan'         :  'np.tan'  ,        
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

# Define default sequence of linewidth and linestyle.
# First 3 histograms are solid with different line thickness;
# next 3 are dashed; next 3 are dot-dashed.
HIST_LINE_ARGS = [
    (2.0,'-' ), (1.6, '-' ), (1.2, '-' ),   # solid
    (2.0,'--'), (1.6, '--'), (1.2, '--'),   # dashed
    (2.0,':' ), (1.6, ':' ), (1.2, ':' ),   # dot-dashed
    (None,None)  # dummy that has no comma
]


STAT_NAME_N       = 'N'
STAT_NAME_AVG     = 'avg'
STAT_NAME_MEDIAN  = 'median'
STAT_NAME_STDDEV  = 'stddev'

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
      @@VARIABLE 'zHD:mB - 3.1*c + 0.16*x1'
      @@VARIABLE 'sqrt(MUERR**2 + .05**2):zHD'
      @@VARIABLE 'np.sqrt(MUERR**2 + .05**2):zHD'  # can explicitly define with np.

  Single quotes around expressions with parenetheses are needed to avoid 
  linux problems parsing (). Math functions (sqrt, abs, exp, log, log10)
  are internally updated with 'np.' prefix. There is a hard-wired list of 
  functions to check for missing np, so if using an undefined function
  you can explicitly prepend np (and please post github issue about 
  missing function).

  Multiple variables are overlaid on same plot; e.g.
     @V SNRMAX1  SNRMAX2  SNRMAX3
  If using math symbols or pad spaces, use single quotes, e.g.
     @V  'PEAKMAG_r - PEAKMAG_i'  'PEAKMAG_i - PEAKMAG_z'  'PEAKMAG_Y - PEAKMAG_z'


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

@@PRESCALE
  Select prescaled subset of rows; useful to make plots more quickly
  from very large tables, or to see more detail in a 2D plot.
      @@PRESCALE 7    # select every 7th row for all tables
      @@PRESCALE 1 7  # select every 7th row for 2nd table (e.g., large sim)

@@WEIGHT
  Reweight 1D bin contents with arbitrary math function of x-axis using
  syntax with 'x' as variable:
     @@WEIGHT '1-x+x**2'
     @@WEIGHT 'exp(-0.4*((x+.3)/.2)**2)'
     @@WEIGHT 'exp(-0.4*((x+.3)/.2)**2)*heaviside(x-2,0.5)'

  A weighted plot can be overlaid on original plot (single @@TFILE arg)
  with list of two weight functions, 
     @@TFILE A.TXT  @@WEIGHT  1   '1-x+x**2'
  where the "1" arg means that first plot is not modified. Finally,
  when overlaying plots from two different files, providing two 
  weights applies a separate weight to each file, e.g.
     @@TFILE A.TXT B.TXT  @@WEIGHT 1  'exp(-x/2)'
  aplplies unit weight to A.TXT and exp(-x/2) weight to B.TXT.

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

@@NBIN_AUTO_SCALE
  scale number of bins with automatic bounds; must be integer

@@OPT   {' '.join(VALID_OPT_LIST)}

   {OPT_HIST:<12} ==> 1D or 2D histogram (instead of points with error bars)
                      Default histtype='step' -> unshaded 
   {OPT_HISTFILL:<12} ==> histtype='stepfill' --> shaded histogram
                    For overlaying plots, this option is recommended only if each
                    sequential plot is a subset of previous plot such as, e.g.
                    @@CUT 'x1<3'  @@CUT 'x1<1'  @@CUT 'x1<-1'
   {OPT_NEVT:<12} ==> append N=Nevt on each legend (1D and 2D).
   {OPT_MEAN:<12} ==> 1D->append mean on legend; 2D->overlay y-axis mean in x-bins.
                    @@ERROR or @E (y-axis) -> replace arithmetic mean with 1/ERR^2-weighted mean.
                    @@error or @e -> same as @@ERROR, but suppress error bars on plot.
   {OPT_AVG:<12} ==> same as {OPT_MEAN}.
   {OPT_MEDIAN:<12} ==> overlay y-axis median in x-bins (2D only).
   {OPT_STDDEV:<12} ==> 1D -> append stddev on each legend;
                    2D -> mean error bar is STDDEV instead of STDDEV/sqrt(N)
   {OPT_OV:<12} ==> overlay two 1D plots (2 files or 2 cuts); do not print chi2/dof on plot
   {OPT_OVCHI2:<12} ==> execut OPT_OV and display chi2/dof on plot
   {OPT_DIAG_LINE:<12} ==> draw line with slope=1 for 2D plot.

   {OPT_LOGY:<12} ==> log scale for vertical axis (1D or 2D)
   {OPT_LOGZ:<12} ==> log scale for Z axis (2D HIST only)
   {OPT_GRID:<12} ==> display grid on plot.
   {OPT_LIST_CID:<12} ==> print up to 100 CIDs passing cuts.      
   {OPT_LIST_ROW:<12} ==> print up to 100 ROWs passing cuts.      
   {OPT_DIFF_ALL:<12} ==> 2D: plot y-axis difference in median (or mean) values between plots.
                    Overlay file1-file0, file2-file0, file3-file0, etc, or
                    overlay cut1-cut0, cut2-cut0, cut3-cut0, etc;
                    e.g., @@CUT ''  'c<.1'  'c<0'  'c<-0.1'  
                    where blank cut '' means reference has no cuts.
                    Only median (or mean) is shown; individual data points are suppressed.
   {OPT_DIFF_CID:<12} ==> 2D: plot y-axis difference for each CID between samples.
                    (compare correlated samples only; must have CID overlap).
                    As with {OPT_DIFF_ALL}, can compare multiple files or cuts.
                    Median (or mean) is overlaid over individual data points;
                    to suppress points and only show median(mean), set "@@ALPHA 0".
   {OPT_RATIO:<12} ==> 1D: plot histogram ratio file1/file0 or cut1/cut0.
                    Computes binomial error if numerator is subset of denominator;
                    else adds denom+numerator Poisson uncertainties in quadrature.
                    Can overlay multiple ratios with multiple files or cuts; e.g.
                      @@CUT ''  'c<.1'  'c<0'  'c<-0.1'  
                    where the blank cut '' selects entire sample for denominator.



      INPUTS FOR PLOT STYLE & LABELS
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@@XLABEL 
   Override default x-axis label = variable name
@@YLABEL  
   Override default y-axis label = variable name

@@LEGEND
   Overwrite default legend; must give one arg (in quotes) per plot.
   E.g., For 3 sets of cuts, 
      @@LEGEND 'No cuts'  'SALT2c < 0'  'SALT2c>0'
   For 2 files and "@@OPT DIFF_CID", only need to give one @@LEGEND arg.
   To suppress legend,  @@LEGEND NONE / @@legend none

@@LEGEND_SIDE
  For a busy plot with insufficient space for legend, this option is the
  same as @@LEGEND, except that the legend appears on the right side of 
  the plot (outside the plot box).

@@TITLE
  Text of title to dispaly above plot. For text length > 50 chars,
  the font size is scaled by 50/len(text) so that it doesn't run off
  the plot page.

@@TEXT
  Provide relative coordinates (0 to 1) per axis and text to add on plot.
  Relative coord 0 puts text at left (x-axis) or bottom (y-axis); 
  relative coord of 0.5 puts text in the middle, etc ...
  The arg syntax is groups of   x_rel y_rel 'text' ; e.g, 
  e.g.
    @@TEXT 0.05 0.90 '(a) comment bla bla $alpha = \beta$' 
        or
    @@TEXT 0.05 0.90 '(a) comment 1'   0.10 0.80 'comment 2'  0.10 0.70 'comment 3'

  There is no explicit delimiter and thus each text strimg must be enclosed in quotes
  so that it counts as a single item when parsed. Script aborts if number of TEXT 
  args is not a multiple of 3.

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
    
@@FONTSIZE_SCALE
  scale font sizes for title and axis labels

@@SAVE
  Save figure as pdf or png; e.g.
    @@SAVE my_first_table_plot.pdf
       or
    @@SAVE my_first_table_plot.png
  Plot is either displayed in pop-up window or stored in @@SAVE file ...
  but not both. This enables pipelines to make post-processing plots
  without worrying about pop-windows in slurm jobs.

@@DUMP
  Dump histogram contentes to text file; e.g.
     @@DUMP plot_hist1d.text
  Works only with @@OPT HIST option (1D or 2D)

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
    parser.add_argument('@E', '@@ERROR', default=None, help=msg, nargs="+")

    msg = 'variable to only for weight avg per bin; do not show error bard in 2D plot'
    parser.add_argument('@e', '@@error', default=None, help=msg, nargs="+" )

    msg = 'variable to reweight 1D plot'
    parser.add_argument('@W', '@@WGTVAR', '@w', default=None, help=msg, nargs="+")
    
    msg = "cuts with boolean and algegraic operations. If >1 CUT, plots are overlaid."
    parser.add_argument('@@CUT', '@@cut', default =[None], help=msg, nargs="+")

    msg = "integer prescale for selecting table rows"
    parser.add_argument('@@PRESCALE', '@@prescale', default =[1],
                        type=int, help=msg, nargs="+" )
        
    msg = "WEIGHT function(s) for 1D hist only; e.g., 1-0.5*x+3*x**2"
    parser.add_argument('@@WGTFUN', '@@wgtfun', default=[None], help=msg, nargs="+" )

    msg = "number of rows to read (for large files)"
    parser.add_argument('@@NROWS', '@@nrows', default=None, help=msg, type=int )
    
    # -------- decorative options below ---------
    
    msg = "AUTO (default): optional bounds min, max and binsize of plot parameters.\n" \
          "For 2D plot, must specify both x and y bounds. y-binsize is ignored."
    parser.add_argument('@@BOUNDS', '@@bounds', default=None, help=msg, nargs = '+' )

    msg = "scale auto NBIN  (no @@BOUNDS input)"
    parser.add_argument('@@NBIN_AUTO_SCALE', '@@nbin_auto_scale',
                        default=1, help=msg, type=int)    
    
    msg = "Units to show for each axis; see @@HELP"
    parser.add_argument('@U', '@@UNITS', '@@units', default=None, help=msg, nargs="+")

    msg = "x-axis label (override default)"
    parser.add_argument('@@XLABEL', '@@xlabel', default=None, help=msg, type=str)
    msg = "y-axis label (override default)"
    parser.add_argument('@@YLABEL', '@@ylabel', default=None, help=msg, type=str)        
    
    msg = "Override default plot title"
    parser.add_argument('@@TITLE', '@@title', default=None, help=msg, type=str )

    msg = "Extra text on plot"
    parser.add_argument('@@TEXT', '@@text', default=None, help=msg, nargs="+")    

    msg = "Override default legend on plot (space sep list per TFILE)"
    parser.add_argument('@@LEGEND', '@@legend', default=None, help=msg, nargs="+")

    msg = "Same as @@LEGEND, but place outside plot on right side"
    parser.add_argument('@@LEGEND_SIDE', '@@legend_side', default=None, help=msg, nargs="+")

    msg = "Override default marker='o'"
    parser.add_argument('@@MARKER', '@@marker', default=['o'], help=msg, nargs="+")    

    msg = "scale horizontal size of plot"
    parser.add_argument('@@XSIZE_SCALE', '@@xsize_scale',
                        default=1.0, help=msg, type=float)
    msg = "scale vertical size of plot"
    parser.add_argument('@@YSIZE_SCALE', '@@ysize_scale',
                        default=1.0, help=msg, type=float)    

    msg = "scale font size of axis labels and title"
    parser.add_argument('@@FONTSIZE_SCALE', '@@fontsize_scale',
                        default=1.0, help=msg, type=float)
    
    msg = "Alpha value for plot. Set to 0 to see only averages." \
          "ALPHA=0 and DIFF=True compares average difference between two files, " \
          "even if there are no overlapping CIDS."
    parser.add_argument('@@ALPHA', '@@alpha',  default=[0.8], type=float, help=msg, nargs="+" )
    
    msg = "Filename to save plot.  Default: does not save plots."
    parser.add_argument('@@SAVE', '@@save', default=None, help=msg )

    msg = "Filename to dump histogram contents (works only with @@OPT HIST)"
    parser.add_argument('@@DUMP', '@@dump', default=None, help=msg )
    
    msg = "options; see @@HELP"
    parser.add_argument('@@OPT', '@@opt', help=msg, nargs="+", default = [])

    msg = "debug options (for development only)"
    parser.add_argument('@@DEBUG_FLAG', "@@debug_flag", help=msg, type=int, default=0)    

    msg = "Full help menu printed to stdout"
    parser.add_argument('@H', '@@HELP', help=msg, action="store_true")

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
    # for variable(s)
    if args.VARIABLE:
        args.VARIABLE_ORIG = args.VARIABLE[0]
    else:
        sys.exit(f"\n ERROR: must define variable(s) to plot with @@VARIABLE or @@V")

    # store plot dimension as if it were passed on command line
    args.NDIM = len(args.VARIABLE[0].split(COLON))

    if args.NDIM == 2 and OPT_RATIO in args.OPT:
        sys.exit(f"\n ERROR: '@@OPT RATIO' is not valid for 2D plot")
        
    # CUT is tricky. Make sure that length of cut list matchs length
    # of table-file list ... or vice versa ... make sure that length of
    # tfile list matches length of CUT list.

    n_tfile_orig    = len(args.TFILE)
    n_var_orig      = len(args.VARIABLE)
    n_cut_orig      = len(args.CUT)
    n_wgtfun_orig   = len(args.WGTFUN)
    tfile_list      = copy.copy(args.TFILE)
    var_list        = copy.copy(args.VARIABLE)
    cut_list        = copy.copy(args.CUT)
    wgtfun_list     = copy.copy(args.WGTFUN)

    name_arg_list  = [ '@@TFILE', '@@VARIABLE', '@@CUT', '@@WGTFUN' ]
    n_orig_list    = [ n_tfile_orig, n_var_orig, n_cut_orig, n_wgtfun_orig ]
    arg_list_list  = [ tfile_list,   var_list,   cut_list,   wgtfun_list   ]

    # abort if more than 1 table file and more than 1 cut are requested.
    # However, allow 2 table files and 2 WGTFUNs
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
    var_list     = arg_list_list[1]
    cut_list     = arg_list_list[2]
    wgtfun_list  = arg_list_list[3]    
        
    # - - - -

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx proc_args: args.TFILE     -> {args.TFILE}")
        print(f" xxx proc_args: args.VARIABLE  -> {args.VARIABLE}")        
        print(f" xxx proc_args: args.CUT       -> {args.CUT}")
    
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
    narg_tfile           = len(args.tfile_list)    

    args.var_list        = var_list
    args.cut_list        = cut_list
    args.wgtfun_list     = wgtfun_list
    
    # - - - - -

    args.use_err_list = True  # default
    if args.error:
        args.ERROR = args.error
        args.use_err_list = False  # disable error bar on 2D plot

    ndim       = args.NDIM
    args.UNITS = arg_prep_axis(ndim, 'UNITS', args.UNITS)
    args.ERROR = arg_prep_axis(ndim, 'ERROR', args.ERROR)

    
    if args.BOUNDS:
        args.BOUNDS = ' '.join([str(elem) for elem in args.BOUNDS])

    args = args_prep_DIFF(args)
    
    # if only 1 alpha, make sure there is alpha for each file/cut
    args.alpha_list    = arg_prep_extend_list(args, narg_tfile, args.ALPHA,  1.0   )
    args.marker_list   = arg_prep_extend_list(args, narg_tfile, args.MARKER, 'XXX' )
    args.prescale_list = arg_prep_extend_list(args, narg_tfile, args.PRESCALE, 1   )
    
    args.TITLE       = arg_prep_TITLE(args)
    args.legend_list = arg_prep_legend(args) # must be after prep_DIFF

    args.text_list   = arg_prep_TEXT(args)

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

    narg_tfile = len(args.tfile_list)
    
    # if user has not specified any kind of statistical average
    # for DIFF, then set default MEDIAN
    OPT        = args.OPT
    opt_diff   = get_list_match( [OPT_DIFF_CID, OPT_DIFF_ALL], args.OPT)
    opt_stat   = get_list_match( [OPT_MEDIAN,OPT_MEAN,OPT_AVG], args.OPT)
    args.DIFF  = opt_diff  # e.g., 'DIFF_CID' or 'DIFF_ALL' or None
    narg_tfile = len(args.tfile_list)
    
    # if user does not specify statistic for DIFF, set MEDIAN as default    
    if opt_diff and  opt_stat is None :
        args.OPT.append(OPT_MEDIAN) 


    # for DIFF_ALL option, set ALPHA=0 which results in showing
    # only the MEDAN (or MEAN) 
    if opt_diff == OPT_DIFF_ALL:
        args.ALPHA = [ 0.0 ] * narg_tfile

    
    if args.DIFF and narg_tfile == 1:
            sys.exit(f"\n ERROR '@@OPT {args.DIFF}' does not work with 1 table file;\n" \
                     f"\t need 2 or more table files specified after @@TFILE key.")
        
    return args
    # end args_prep_DIFF

def arg_prep_TITLE(args):

    if  args.TITLE:
        plot_title = args.TITLE
    else:
        if args.CUT[0] :
            plot_title = str(args.CUT)
        else:
            plot_title = more_human_readable(str('  '.join(args.VARIABLE)))
    
        if STR_np in plot_title:
            plot_title = plot_title.replace(STR_np,'')    

    return plot_title
    # end arg_prep_TITLE

def arg_prep_TEXT(args):

    # if user input TEXT is
    #   .05 .90 '(a) text1'   0.06 0.80 'text2'   0.06 0.70 'text3'
    # then return
    #   text_list = [
    #                 [ 0.06, 0.90, '(a) text1'] ,
    #                 [ 0.06, 0.80, 'text2' ] ,
    #                 [ 0.06, 0.70, 'text3' ]
    #                ]
    #
    # The user text must be inside quotes since there is no explicit
    # delimiter between text times. The number of args must  be a multiple of 3.

    TEXT = args.TEXT
    if TEXT is None: return None

    text_list = []
    n_arg = len(TEXT)
    if n_arg % 3 != 0:
        sys.exit(f"\nERROR: {n_arg} TEXT args is not a multiple of 3 " \
                 f"-> invalid\n\t Check for missing quotes.\n\t TEXT={TEXT}")

    for i in range(0,n_arg,3):
        x_rel = float(TEXT[i+0])
        y_rel = float(TEXT[i+1])
        text  = TEXT[i+2]
        text_list.append( [ x_rel, y_rel, text] )
            
    return text_list  # end arg_prep_TEXT

def arg_prep_OPT(args):

    # abort on invalid @@OPT, and change all OPT to upper case    
    invalid_OPT_list = []
    args_upper_list = []
    OPT_orig = args.OPT
    OPT_out  = OPT_orig  # default

    # fix lower case OPTs to be upper case; e.g, median -> MEDIAN
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

    # LIST_ROW and LIST_CID are the same, but internally we only use
    # the LIST_CID flag
    if OPT_LIST_ROW in OPT_out and OPT_LIST_CID not in OPT_out:
        OPT_out.append(OPT_LIST_CID)

    # check OV[CHI2] options
    if OPT_CHI2 in OPT_out: # check for legacy option that means OVCHI2
        OPT_out.append(OPT_OVCHI2)

    if OPT_OVCHI2 in OPT_out and OPT_OV not in OPT_out:
        OPT_out.append(OPT_OV)

    if OPT_HISTFILL in OPT_out :
        OPT_out.append(OPT_HIST)

    if args.DUMP and OPT_HIST not in OPT_out:
        sys.exit(f"\n ERROR: invalid @@DUMP without @@OPT HIST ")
        
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

    if args.LEGEND_SIDE:
        args.LEGEND = args.LEGEND_SIDE
        
    LEGEND_orig = args.LEGEND
    LEGEND_out  = LEGEND_orig   # default is user input

    NDIM     = args.NDIM
    do_diff  = args.DIFF
    do_ratio = OPT_RATIO in args.OPT
    narg_var = len(args.VARIABLE)
    narg_cut = len(args.CUT)
    
    if LEGEND_orig is None:
        # no user supplied legend, so make up reasonable default legend
        LEGEND_out = [ None ] * len(args.TFILE)

        if narg_var > 1:
            LEGEND_out = args.VARIABLE
        elif narg_cut > 1 :
            LEGEND_out = args.CUT
        else :
            LEGEND_out = []
            for t in args.tfile_base_list:
                legend = t.split('.')[0]  # base file name without extension after dot
                LEGEND_out.append(legend)

        LEGEND_human_readable = []
        for legend in LEGEND_out:
            LEGEND_human_readable.append(more_human_readable(legend))
        LEGEND_out = LEGEND_human_readable
        
    elif LEGEND_orig[0].upper() == "NONE" :
        args.LEGEND = SUPPRESS

    # - - - - - 
    narg_legend  = len(LEGEND_out)
    narg_tfile   = len(args.tfile_list)

    # check option to allow NTFILE-1 legends 
    allow_missing_legend = do_diff or do_ratio
    if allow_missing_legend and narg_legend == narg_tfile-1:
        LEGEND_out = ['dummy'] + LEGEND_out
        narg_legend  = len(LEGEND_out)

    match_narg   = narg_tfile == narg_legend
    if  not match_narg:
        sys.exit(f"ERROR: narg_tfile={narg_tfile} but narg_legend={narg_legend}; " \
                 f"narg_legend must match number of table files.")

    return LEGEND_out
    # end arg_prep_legend

    
def arg_prep_extend_list(args, narg_tfile, arg_list_orig, arg_missing):

    # if arg_list_orig has fewer elements than number of table files,
    # extend list to be same length as number or table files;
    # else return original list.  This allows users to specify one
    # value (e.g. for alpha or maker) that is used on all plots.
    # Input arg_missing is used to fill missing arg.
    
    narg_orig = len(arg_list_orig)
    
    # check option to allow one missing 
    allow_missing = args.DIFF or  OPT_RATIO in args.OPT
    if allow_missing and narg_orig == narg_tfile-1:
        arg_list_out = [ arg_missing ] + arg_list_orig    
        return arg_list_out
    
    if narg_orig == 1 and narg_tfile > 1:
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
        
            
    j_last = 0 ; j=0
    while j < len_str :
        ch         = STRING_local[j:j+1]  # current char
        ch2        = STRING_local[j:j+2]  # allow 2-char delim; e.g., >=
        is_last    = (j == len_str-1)

        delim = None
        if ch2 in DELIM_LIST:            
            delim = ch2  # 2-char delims take priority over single-char delim
        elif ch in DELIM_LIST:
            delim = ch
            
        is_delim  = delim is not None        
        
        str_last = None
        j_last_dump = j_last
        if is_delim:
            len_delim = len(delim)            
            str_last = STRING_local[j_last:j].replace(' ','')
            # xxx mark delete j_last = j+1
            j_last = j + len_delim 
        elif is_last:
            str_last = STRING_local[j_last:j+1].replace(' ','')
            
        valid_str_last = (str_last is not None and len(str_last)>0)

        if valid_str_last:
            split_list.append(str_last)
            if is_number(str_last) or "'" in str_last or '"' in str_last :
                strtype_list.append(STRTYPE_NUM)
            elif str_last in FUNC_LIST:
                strtype_list.append(STRTYPE_FUNC)
            else:
                strtype_list.append(STRTYPE_VAR)

                
        if is_delim:
            if delim == '=' : delim = '=='  # allow user to specify single '='
            split_list.append(delim)
            strtype_list.append(STRTYPE_DELIM)

        if args.DEBUG_FLAG_DUMP_TRANSLATE2 :
            print(f"\t xxx - - - - - - - - ")
            print(f"\t xxx j={j}  j_last={j_last_dump}  ch='{ch}' or '{ch2}'  " \
                  f"is_[last,delim]={is_last},{is_delim}  " \
                  f"str_last={str_last}  valid_str_last={valid_str_last}")
            print(f"\t xxx split_list -> {split_list}")

        # increment char counter by 1, unless delimiter has 2 chars
        if is_delim:
            j += len_delim  # advance 1 or 2 chars for delimeter
        else:
            j += 1          # advance 1 char for everything else
        
    # - - - -
    return split_list, strtype_list
    # end split_var_string
    

def translate_driver(args):

    # add df. as needed to args.VARIABLE
    args.raw_var_list = []
    for ivar, var in enumerate(args.var_list):
        df_var, raw_plotvar_list = translate_VARIABLE(var)
        args.var_list[ivar]      = df_var
        args.raw_var_list       += raw_plotvar_list  # used to check existence of vars
    logging.info('')
    
    for icut, cut in enumerate(args.cut_list):
        cut_list, raw_cutvar_list = translate_CUT(cut)  # add df. and df.loc as needed
        args.cut_list[icut]       = cut_list
        args.raw_var_list        += raw_cutvar_list
    logging.info('')
    
    # tack on error vars
    if args.ERROR :
        var_err_list = args.ERROR.split(COLON)
        for var in var_err_list:
            if len(var) > 0: args.raw_var_list += [var]
    
    # remove duplicates
    args.raw_var_list = list(set(args.raw_var_list))
    
    if args.DEBUG_FLAG_DUMP_TRANSLATE :
        print(f" xxx raw_var_list = {args.raw_var_list}")
        sys.exit(f"\n xxx bye .")
    
    return  # end translate_driver

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
    table_var_list = []  # list of all collumn var names used in VARIABLE.
    
    for tmp_str, tmp_type in zip(split_list, strtype_list):
        tmp_append = tmp_str  # default string to append
        if tmp_type == STRTYPE_VAR:
            tmp_append = STR_df + tmp_str
            table_var_list.append(tmp_str)
        elif tmp_type == STRTYPE_FUNC :
            tmp_append = STR_np + tmp_str
        else:
            pass
        
        VARIABLE_df += tmp_append    
        
    logging.info(f"Translate VARIABLE {VARIABLE_ORIG}")
    logging.info(f"             ->    {VARIABLE_df}")
    
    return VARIABLE_df, table_var_list

    # end translate_VARIABLE
    
def translate_CUT(CUT):

    # add df.loc, df. and () as needed to input cut
    # Add "(df." in front of each var_list element that is NOT a number.
    # Add ")" after each var_list element that is a number.
    #
    # July 24 2024 rewrite logic to split by bool and (),
    #   then wrap each item as '(df.' + item + ')'

    CUT_ORIG = CUT
    if not CUT:       return CUT, []
    if STR_df in CUT: return CUT, []

    CUT = CUT.replace(' ','')

    
    # split into sections separated by boolean &, |, or ()
    cut_list = re.split(r"\(|\)|\&|\|", CUT.replace(' ',''))
    cut_list = list(filter(None, cut_list))

    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx CUT      = {CUT}" )
        print(f" xxx cut_list = {cut_list}")

    # for each cut item, append parentheses and df.
    cut_list_df = []
    raw_var_list = []
    is_func_list   = []
    icut = 0
    icut_func = -9
    
    for cut in cut_list:

        if cut in NUMPY_FUNC_DICT:
            is_func = True
            icut_func = icut
        else:
            is_func = False
        is_func_list.append(is_func)  
        
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
                if tmp_str not in raw_var_list: raw_var_list.append(tmp_str)
            elif tmp_type == STRTYPE_FUNC:
                tmp_append = tmp_str
                #sys.exit(f"\n ERROR: cannot use numpy functons in @@CUT; " \
                #         f"remove '{tmp_str}' from @@CUT arg")
            else:
                pass
            cut_df += tmp_append
            

        #if STR_np not in cut_df:
        add_parenth = True

        # Dec 26 2024: special logic for functions; add np. instead of ()
        if is_func:
            add_parenth = False
            cut_df = STR_np + cut_df
        if icut - icut_func <= 2:
            add_parenth = False
        
        if add_parenth:
            cut_df = '(' + cut_df + ')'

        cut_list_df.append(cut_df)
        icut += 1
        
    # - - - - -
    if args.DEBUG_FLAG_DUMP_TRANSLATE:
        print(f" xxx cut_list_df = {cut_list_df}")
        
    # replace each original cut with cut_df in CUT
    CUT_df = CUT

    # xxxxxxx mark delete Dec 26 2024 xxxxxxx
    #n_cut = len(cut_list)
    #for icut in range(0,n_cut):
    #    cut     = cut_list[icut]
    #    cut_df  = cut_list_df[icut]
    #    is_func = is_func_list[icut]        
    #    
    #    CUT_df  = CUT_df.replace(cut,cut_df) 
    #    if args.DEBUG_FLAG_DUMP_TRANSLATE:
    #        print(f" xxx \t {icut}: replace user cut {cut} --> {cut_df} " \
    #              f" (is_func={is_func})")
    #    icut += 1
    # xxxxxxxxxxx end mark xxxxxxxxx
    

    for cut, cut_df, is_func in zip(cut_list, cut_list_df, is_func_list):
        CUT_df = CUT_df.replace(cut,cut_df)
        if args.DEBUG_FLAG_DUMP_TRANSLATE:
            print(f" xxx \t replace user cut {cut} --> {cut_df}  (is_func={is_func})")

    # finally, wrap entire cut in df.loc[ CUT ]
    CUT_df = STR_df_loc + '[' + CUT_df + ']'    

    logging.info(f"Translate CUT {CUT}  ")
    logging.info(f"          ->  {CUT_df}")
    
    return CUT_df, raw_var_list


def is_number(string): 
    try: 
        float(string) 
        return True
    except ValueError: 
        return False
    return

def set_axis_dict(args, plot_info, var):

    # copute and store plot_info.axis_dict
    
    axis_dict       = {}  # x & y varnames, err varnames
    
    # - - - -
    NDIM         = args.NDIM
    VAR_LIST     = var.split(COLON)  
    VARERR_LIST  = args.ERROR.split(COLON)
    
    for n, VAR in enumerate(VAR_LIST):

        
        STR_VAR      = str(VAR)
        VARERR  = VARERR_LIST[n]
        if len(VARERR) > 0:
            STR_VARERR   = STR_df + str(VARERR)
        else:
            STR_VARERR   = None            
            
        if n == 0:
            axis_dict['x']       = STR_VAR
            axis_dict['xerr']    = STR_VARERR            
        else:
            axis_dict['y']       = STR_VAR
            axis_dict['yerr']    = STR_VARERR

    axis_dict['set_ylim'] = (NDIM == 2 or OPT_RATIO in args.OPT)
    
    plot_info.axis_dict_list.append(axis_dict)
    
    #sys.exit(f"\n xxx set_axis_dict: axis_dict = {axis_dict}")
    
    return  # end set_axis_dict

def set_axis_labels(args, plot_info):

    # set x and y axis labels in plot_info.axis_dict.
    # If user does not define @@XAXIS (@@YAXIS), define reasonable defaults
    # based on user options.
    
    NDIM          = args.NDIM
    OPT           = args.OPT
    UNITS_LIST2D  = args.UNITS.split(COLON)

    VAR_LIST2D    = args.VARIABLE[0].split(COLON)  # x or x & y varnames of first set of vars
    
    WGTFUN       = args.WGTFUN[0]  # remove [] for axis title 
    if len(args.WGTFUN) > 1:
        WGTFUN = '[' + ' , '.join(args.WGTFUN) +']'    # show list of all weights for overlay
            
    axis_dict = plot_info.axis_dict_list[0]  # ?? fragile ??
    xbin      = plot_info.bounds_dict['xbin']

    # - - - - 
    # avoid crazy number of digits for ylabel 'counts per {xbin}'
    str_xbin = str(xbin)
    if '.' in str_xbin:
        n_digit = len(str_xbin.split('.')[1])  # number of digits after decimal
        if n_digit > 4 :
            str_xbin = f'{xbin:.4f}'

    xlabel = None
    ylabel = None
    
    for n, VAR in enumerate(VAR_LIST2D):
            
        LABEL    = VAR.replace(STR_df,'') 
        LABEL    = LABEL.replace(STR_np,'')

        U = UNITS_LIST2D[n]
        if len(U) > 0:  LABEL += '  (' + U + ')'
    
        if n == 0:
            # set default x-axis label
            xlabel = more_human_readable(LABEL)

            # for 1D plots, set default y-axis label based on user optons
            if NDIM == 1:
                ylabel  = f'Counts per {str_xbin}' 
                if WGTFUN           : ylabel = f'Counts * {WGTFUN} per {str_xbin}'
                if OPT_RATIO in OPT : ylabel = 'Ratio' 

            # check user overrides
            if args.XLABEL   : xlabel = args.XLABEL 
            if args.YLABEL   : ylabel = args.YLABEL 
        else:
            # for 2D plots, default axis label is varname
            ylabel = more_human_readable(LABEL)
            if args.DIFF      : ylabel += ' Diff'            
            if args.YLABEL    : ylabel = args.YLABEL  # user overrride

    # - - - - - - 
    axis_dict['xaxis_label']  = xlabel
    axis_dict['yaxis_label']  = ylabel            
    plot_info.axis_dict_list[0]    = axis_dict
    
    return  # end set_axis_labels

def more_human_readable(label_orig):

    # For input string label_orig, replace math symbols with
    # math mode (e.g, replace dash with real minus sign) and
    # add pad spacing around some math symbols for better readability.

    label_out = label_orig
    math_symbol_list = [ '+', '-', '/', ':' ]

    for sym in math_symbol_list:
        label_out = label_out.replace(sym,f' ${sym}$ ')
    
    return label_out
    # end more_human_readable
    
def set_hist2d_args(args, plot_info):

    
    do_hist2d   = args.NDIM == 2 and (OPT_HIST in args.OPT)
    if not do_hist2d:
        plot_info.hist2d_args = None
        return None

    bounds_dict = plot_info.bounds_dict

    hist2d_args         = Namespace()
    if bounds_dict['custom'] :        
        x_range = [ bounds_dict['xmin'], bounds_dict['xmax'] ]
        y_range = [ bounds_dict['ymin'], bounds_dict['ymax'] ]
        hist2d_args.range =  [ x_range, y_range ]
        hist2d_args.bins  =  ( bounds_dict['nxbin'] , bounds_dict['nybin'] )
                    
    else:
        # load None values to pass to plt.hist2d
        NBIN_SCALE = args.NBIN_AUTO_SCALE
        hist2d_args.bins    = (NXYBIN_AUTO*NBIN_SCALE, NXYBIN_AUTO*NBIN_SCALE)
        hist2d_args.range   = None
        
    # logz option applies for custom or auto bounds
    if OPT_LOGZ in args.OPT:
        hist2d_args.norm = colors.LogNorm()
    else:
        hist2d_args.norm = None
        
    logging.info(f"hist2d_args = {hist2d_args}")
    plot_info.hist2d_args = hist2d_args
                    
    return   # end set_hist2d_args
    
def set_custom_bounds_dict(args,plot_info):

    # set plot_info.bounds_dict for @@BOUNDS input.
    # Also set hist2d_args if this option is set.
    
    bounds_dict     = {}  # user/custom bounds
    bounds_dict['custom'] = False
    
    BOUNDS   = args.BOUNDS

    if BOUNDS is None:
        plot_info.bounds_dict  = bounds_dict
        return
    
    # - - - - - 
    for n, axis_bounds in enumerate(BOUNDS.split(':')):
        bounds_dict['custom'] = True
        axis_bounds_values  = [float(i) for i in axis_bounds.split() ]

        amin = float(axis_bounds_values[0])   # axis min val
        amax = float(axis_bounds_values[1])   # axis max val
        abin = float(axis_bounds_values[2])   # axis bin size
        if abin > 0.0:
            nbin = (amax-amin+1.0e-8)/abin        # axis nbin
        else:
            nbin = 0  # allowed for y-axis if not used
                
        load_xbound = (n == 0)
        load_ybound = (n  > 0)
            
        # check that bin size is compatible with range for x-axis
        # or for hist2d y-axis
        do_check_nbin = load_xbound or (load_ybound and OPT_HIST in args.OPT)
        if do_check_nbin:
            if abs(nbin-int(nbin)) > 1.0e-5 :
                sys.exit(f"\n ERROR: invalid bin size for  " \
                         f"min max bin = '{axis_bounds}' in @@BOUNDS {BOUNDS}")
            
        if load_xbound:
            bounds_dict['x']     = axis_bounds_values
            bounds_dict['xmin']  = amin
            bounds_dict['xmax']  = amax
            bounds_dict['xbin']  = abin
            bounds_dict['nxbin'] = int(nbin)

        if load_ybound:
            bounds_dict['y']     = axis_bounds_values
            bounds_dict['ymin']  = amin
            bounds_dict['ymax']  = amax
            bounds_dict['ybin']  = abin
            bounds_dict['nybin'] = int(nbin)

    # - - - - - -  - - 
    # error checks 
    if args.DIFF and  args.NDIM == 1 :
        sys.exit("\nERROR: @@OPT {args.DIFF}  does not work for 1D histograms." \
                 " ABORT to avoid confusion.")

    # load output namespace
    plot_info.bounds_dict         = bounds_dict
    
    return # end set_custom_bounds_dict

def read_tables(args, plot_info):

    tfile_list      = args.tfile_list
    tfile_base_list = args.tfile_base_list
    var_list        = args.var_list
    cut_list        = args.cut_list
    ps_list         = args.prescale_list
    wgtfun_list     = args.wgtfun_list
    legend_list     = args.legend_list
    alpha_list      = args.alpha_list
    marker_list     = args.marker_list
    NROWS           = args.NROWS

    axis_dict_list = plot_info.axis_dict_list

            
    MASTER_DF_DICT = {}  # dictionary of variables to plot (was MASTERLIST)
    nf = 0

    same_tfiles = (len(set(tfile_list)) == 1)
    nrow_tot = 0
    
    for tfile, var, axis_dict, cut, prescale, wgtfun, legend, alpha, marker in \
        zip(tfile_list, var_list, axis_dict_list,
            cut_list, ps_list, wgtfun_list, legend_list, alpha_list, marker_list):

        varname_x      = axis_dict['x']
        varname_nodf_x = varname_x.replace(STR_df,'')  # for diagnostic print 
        varname_xerr   = axis_dict['xerr']
        if args.NDIM == 2:
            varname_y      = axis_dict['y']
            varname_nodf_y = varname_y.replace(STR_df,'')  # for diagnostic print         
            varname_yerr   = axis_dict['yerr']    
        
        tfile_base = os.path.basename(tfile)
        if not os.path.exists(tfile):
            sys.exit(f"\n ERROR: cannot find {tfile}")

        store_df_ref = same_tfiles and nf == 0
        copy_df_ref  = same_tfiles and nf >  0

        if copy_df_ref:
            logging.info(f"Copy duplicate df for {tfile_base}")
            df = copy.deepcopy(df_ref)
        else:
            logging.info(f"Read {tfile_base}")
            # Get varname info including number of rows to skip (e.g. skip DOCANA),
            # name of first column identifier, and check that requested
            # variables exist.
            # count number of rows to skip before VARNAMES; e.g., skip DOCANA for HOSTLIB 
            varname_idrow, nrow_skip = check_table_varnames(tfile,args.raw_var_list)
            plot_info.varname_idrow  = varname_idrow            
            usecol_list = [ varname_idrow ] + args.raw_var_list

            # read table and store in data frame.
            # only read needed columms to reduce memory consumption.
            df  = pd.read_csv(tfile, comment="#", sep=r"\s+",
                              usecols  = usecol_list,
                              skiprows = nrow_skip,
                              nrows    = NROWS )

        if store_df_ref:    df_ref = copy.deepcopy(df)
                
        try:
            df[varname_idrow] = df[varname_idrow].astype(str)
        except KeyError:
            logging.warn(f"No {varname_idrow} present in this file. OK for some file types.")

        # apply user cuts

        if cut:
            df = eval(cut)

        if prescale > 1:
            df = df.iloc[::prescale]
            
        # drop duplicates for DIFF CID matching
        if args.DIFF == OPT_DIFF_CID:
            df.drop_duplicates(varname_idrow, inplace=True)

        # increment MASTER_DF_DICT dictionary; note that filename index is dict key.
        key = f"tf{nf}"
        nf += 1
        
        MASTER_DF_DICT[key] = df
        nrow = len(df)
        nrow_tot   += nrow
        name_legend = legend
            
        logging.info(f"\t Read nrow={nrow}  for {name_legend}")

        if nrow == 0: 
            logging.warning(f"zero rows read for {name_legend}")

        MASTER_DF_DICT[key] = {
            'df'           : df,
            'wgtfun'       : wgtfun,
            'name_legend'  : name_legend,
            'alpha'        : alpha,
            'marker'       : marker
        }

        try:
            # xxx MASTER_DF_DICT[key]['df']['x_plot_val'] = eval(varname_x)
            MASTER_DF_DICT[key]['df'].loc[:,'x_plot_val'] = eval(varname_x)
            xmin = np.amin(df['x_plot_val'])
            xmax = np.amax(df['x_plot_val'])
            logging.info(f"\t x-axis({varname_nodf_x}) range : {xmin} to {xmax}") 
            if varname_xerr is not None:
                MASTER_DF_DICT[key]['df'].loc[:,'x_plot_err'] = eval(varname_xerr)
            
            if args.NDIM == 2:
                MASTER_DF_DICT[key]['df'].loc[:,'y_plot_val'] = eval(varname_y)
                ymin = np.amin(df['y_plot_val'])
                ymax = np.amax(df['y_plot_val'])
                logging.info(f"\t y-axis({varname_nodf_y}) range : {ymin} to {ymax}")  
                if varname_yerr is not None:
                    # xxx MASTER_DF_DICT[key]['df']['y_plot_err'] = eval(varname_yerr)
                    MASTER_DF_DICT[key]['df'].loc[:,'y_plot_err'] = eval(varname_yerr)
            else:
                ymin = None
                ymax = None
            
        except AttributeError:
            sys.exit(f"\n ERROR: Couldn't set bounds for axis_dict={axis_dict} and {tfile}")


        # store min and max for each variable to plot
        table_bounds_key_list = [ 'xmin', 'xmax', 'ymin', 'ymax' ]
        table_bounds_val_list = [  xmin,   xmax,   ymin,   ymax  ]
        for tmp_key, tmp_val in zip(table_bounds_key_list, table_bounds_val_list):
            MASTER_DF_DICT[key][tmp_key] = tmp_val
    
    # - - - - - - - -
    logging.info("Finished loading all table files.")
    logging.info("# - - - - - - - - - - ")

    if nrow_tot == 0:
        logging.fatal(f"zero rows read among all tables ")
        sys.exit('\n\t ABORT')
            
    # load output namespace
    plot_info.MASTER_DF_DICT = MASTER_DF_DICT
    
    return
    # end read_tables

def set_xbins(args, plot_info):

    # set xbins and xbins_cen to pass later as arguments to matplot routines
    
    bounds_dict    = plot_info.bounds_dict
    MASTER_DF_DICT = plot_info.MASTER_DF_DICT

    custom      = bounds_dict['custom']

    if custom:
        # most info already set in set_custom_bounds_dict
        xmin  = bounds_dict['xmin']
        xmax  = bounds_dict['xmax']
        nxbin = bounds_dict['nxbin']

    else:
        # fetch xmin/maxfrom the data that was read in read_tables()
        # Make sure get min/max among all variables plotted
        xmin =  1.0E20
        xmax = -1.0e20
        for key_tf  in MASTER_DF_DICT:
            xmin = min(xmin,MASTER_DF_DICT[key_tf]['xmin'])
            xmax = max(xmax,MASTER_DF_DICT[key_tf]['xmax'])
                
        if xmin == xmax :
            xmin -= 1.0  # hack to avoid crash
            xmax += 1.0
            
        nxbin = NXBIN_AUTO * args.NBIN_AUTO_SCALE
        
    #  - - - -

    xbins       = np.linspace(xmin, xmax, nxbin+1)
    xbins_cen   = ( xbins[1:] + xbins[:-1] ) / 2.  # central value for each xbin
        
    bounds_dict['xbins']     = xbins
    bounds_dict['xbins_cen'] = xbins_cen

    if not custom :
        bounds_dict['xmin']  = xbins[0]                                 
        bounds_dict['xmax']  = xbins[-1]
        bounds_dict['xbin']  = xbins[1] - xbins[0]
        bounds_dict['nxbin'] = nxbin
    
    # - - - - - 

    #sys.exit(f"\n xxx xbins(arange) = \n{xbins}\n xbins_cen=\n{xbins_cen} \n")
        
    plot_info.bounds_dict = bounds_dict

    msg = f"x-axis has {nxbin} bins from {xmin} to {xmax} "
    logging.info(msg)
    
    return   # end set_xbins

def check_table_varnames(tfile, var_list):

    # count number of rows to skip before VARNAMES key,
    # and read name of first column after VARNAMES key
    # Also check that all variables in var_list are specified in VARNAMES.
    
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
            varname_idrow = wdlist[1]
            n_missing = 0

            for var in var_list:
                if var not in wdlist:
                    n_missing += 1
                    logging.error(f"Missing {var} in {tfile}")
            if n_missing > 0:
                logging.fatal(f"Missing {n_missing} variables; check @V and @@ERROR args.")
                sys.exit(f"\n\t ABORT")
            break
        else:
            nrow_skip += 1            
            continue

    t.close()
    logging.info(f"\t Found {nrow_skip} rows to skip before 'VARNAMES:' key")

    return varname_idrow, nrow_skip

    # end count_rows_to_skip

def check_vars_exist(args, df, tfile):

    # xxxxxx OBSOLETE xxxxxxxx
    
    # strip off list of variables used for plotting (@V) or cut (@@CUT).
    # This list has no df. so it's easy to check existence.
    
    raw_var_list = args.raw_var_list

    varlist_missing = []
    for tmp_var in raw_var_list:
        if tmp_var not in df:
            varlist_missing.append(tmp_var)
            logging.error(f"Could not find requested {tmp_var} for {tfile}")

    # xxxxxx OBSOLETE xxxxxxxx            
    n_var_missing = len(varlist_missing)
    if  n_var_missing > 0 :
        logging.fatal(f"Missing {n_var_missing} variables in table file(s); " + \
                      f"check @@VARIABLE and  @@ERROR args")
        sys.exit(f"\n\t ABORT")
    # xxxxxx OBSOLETE xxxxxxxx
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
    NDIM          = args.NDIM    
    OPT           = args.OPT
    do_list_cid   = OPT_LIST_CID  in OPT 
    
    MASTER_DF_DICT       = plot_info.MASTER_DF_DICT
    bounds_dict          = plot_info.bounds_dict
    varname_idrow        = plot_info.varname_idrow  # e.g., 'CID' or 'GALID' 
    
    xbins          = bounds_dict['xbins']
    xbins_cen      = bounds_dict['xbins_cen']    
    # - - - - - - - 

    numplot = 0
    numplot_tot = len(MASTER_DF_DICT)
    info_plot_dict = { }
    
    for key_name, df_dict in MASTER_DF_DICT.items():

        df          = df_dict['df'] 
        name_legend = df_dict['name_legend']        
        logging.info(f"Plot {name_legend}")

        info_plot_dict['do_plot_errorbar']    =  False
        info_plot_dict['do_plot_hist']        =  False
        info_plot_dict['do_ov1d_hist']        =  False
        info_plot_dict['do_ov2d_binned_stat'] =  False
        info_plot_dict['numplot']             =  numplot
        info_plot_dict['xbins']               =  xbins
        info_plot_dict['xbins_cen']           =  xbins_cen        
        info_plot_dict['df_dict']             =  df_dict        
        info_plot_dict['varname_idrow']       =  varname_idrow
        
        if NDIM == 1:
            # 1D
            get_info_plot1d(args, info_plot_dict)        
        else:
            # 2D
            get_info_plot2d(args, info_plot_dict)      

        # strip off arguments to pass to matplotlib ...
        do_plot_errorbar    = info_plot_dict['do_plot_errorbar']
        do_plot_hist        = info_plot_dict['do_plot_hist']
        do_plot_hist_ov1d   = info_plot_dict['do_plot_hist_ov1d']
        do_ov2d_binned_stat = info_plot_dict['do_ov2d_binned_stat']
            
        xval_list   = info_plot_dict['xval_list'] 
        yval_list   = info_plot_dict['yval_list']
        wgt_ov      = info_plot_dict['wgt_ov']
        
        if args.use_err_list:        
            xerr_list   = info_plot_dict['xerr_list'] 
            yerr_list   = info_plot_dict['yerr_list']
        else:
            xerr_list = None
            yerr_list = None

        plt_histtype  = info_plot_dict['plt_histtype']                 
        plt_size      = info_plot_dict['plt_size']   # depends on nevt for 2D
        plt_legend    = info_plot_dict['plt_legend'] # can be appended with more info
        plt_text_dict = info_plot_dict['plt_text_dict']
        
        plt_alpha   = df_dict['alpha']  # fixed by user
        plt_marker  = df_dict['marker'] # fixed by user

        lwid = HIST_LINE_ARGS[numplot][0]  # for 1D hist only
        lsty = HIST_LINE_ARGS[numplot][1]          
        # - - - - -        
        
        if do_plot_errorbar :
            # default
            if plt_alpha > 0 :
                plt.errorbar(xval_list, yval_list,
                             xerr=xerr_list, yerr=yerr_list, 
                             fmt=plt_marker, label=plt_legend,
                             markersize=plt_size, alpha=plt_alpha )
            
            if do_ov2d_binned_stat and NDIM == 2  :
                overlay2d_binned_stat(args, info_plot_dict, None)

        elif do_plot_hist and NDIM == 1 :
            (contents_1d, xedges, patches) = \
                plt.hist(df.x_plot_val, xbins, alpha=plt_alpha, histtype=plt_histtype,
                     label = plt_legend, linewidth=lwid, linestyle=lsty)

            dump_hist1d_contents(args, xedges, contents_1d)
            
        elif do_plot_hist and NDIM == 2 :
            hist2d_args = plot_info.hist2d_args
            contents_2d, xedges, yedges, im = \
                plt.hist2d(df.x_plot_val, df.y_plot_val, label=plt_legend,
                           cmin=.1, alpha=plt_alpha, cmap='rainbow_r', 
                           bins  = hist2d_args.bins,
                           range = hist2d_args.range,
                           norm  = hist2d_args.norm  )
            plt.colorbar(im)
            dump_hist2d_contents(args, xedges, yedges, contents_2d) 
            
            if do_ov2d_binned_stat and NDIM == 2 :
                overlay2d_binned_stat(args, info_plot_dict, 'ORANGE')

            
        elif do_plot_hist_ov1d and NDIM == 1:
            # 1D, typically sim overlaid on data
            plt.hist(df.x_plot_val, xbins, alpha=plt_alpha, histtype='step',
                     weights = wgt_ov, label = plt_legend,
                     linewidth=lwid, linestyle=lsty)
            
        else:
            # nothing to plot; e.g, 1st file for DIFF or RATIO option
            numplot += 1
            continue

        # - - - - -
        # check option to dump CIDs (or ROWs) for events passing cuts
        if do_list_cid:
            print_cid_list(df, name_legend)

        # check for misc plt options (mostly decoration)
        if numplot == numplot_tot-1:
            apply_plt_misc(args, plot_info, plt_text_dict)            

        numplot += 1
        
        
    return   # end plotter_func_driver


def dump_hist1d_contents(args, xedges, contents_1d):
    # dump 1d-histogram contents to text file.

    dump_file = args.DUMP
    if dump_file is None: return

    logging.info(f"Dump HIST1D bin contents to {dump_file}")
    f = open(dump_file, "wt")
    
    xbins_cen   = ( xedges[1:] + xedges[:-1] ) / 2.  # central value for each xbin
    contents_1d[np.isnan(contents_1d)] = 0
    sum_contents = sum(contents_1d)
    nbin         = len(contents_1d)
    
    var_orig = args.VARIABLE[0]
    varname = get_varname_dump(var_orig)
    f.write(f"# Sum of contents: {sum_contents}\n")
    f.write(f"# Number of {var_orig} bins:  {nbin}\n")    
    f.write(f"# x-axis: {var_orig}\n")
    f.write(f"\n")
    f.write(f"VARNAMES:  ROW  {varname}  CONTENTS\n")
    
    rownum = 0
    for x, x_content in zip(xbins_cen, contents_1d) :
        rownum += 1
        f.write(f"ROW:  {rownum:3d}   {x:.5f}  {x_content}\n")
        
    f.close()
    
    return  # end dump_hist1d_contents

def dump_hist2d_contents(args, xedges, yedges, contents_2d) :

    # dump 2D histogram contents to text file
    dump_file = args.DUMP
    if dump_file is None: return

    logging.info(f"Dump HIST2D bin contents to {dump_file}")

    var_list  = args.VARIABLE[0].split(COLON)
    varname_x = get_varname_dump(var_list[0])
    varname_y = get_varname_dump(var_list[1])
    
    xbins_cen   = ( xedges[1:] + xedges[:-1] ) / 2.  # central value for each xbin
    ybins_cen   = ( yedges[1:] + yedges[:-1] ) / 2.  # central value for each ybin
    contents_2d[np.isnan(contents_2d)] = 0  # set NaN contents to zero
    sum_contents = sum(sum(contents_2d))
    nxbin        = len(contents_2d)
    nybin        = len(contents_2d[0])
    nbin         = nxbin * nybin
    rownum = 0

    f = open(dump_file, "wt")
    f.write(f"# Sum of contents: {sum_contents}\n")
    f.write(f"# Number of {var_list[0]:<12} bins:  {nxbin}\n")
    f.write(f"# Number of {var_list[1]:<12} bins:  {nybin}\n")    
    f.write(f"# Number of 2D bins:  {nbin}\n")    
    f.write(f"# x-axis: {var_list[0]}\n")
    f.write(f"# y-axis: {var_list[1]}\n")
    f.write(f"\n")
    f.write(f"VARNAMES:  ROW  {varname_x}  {varname_y}  CONTENTS\n")
    
    for x, y_slice in zip(xbins_cen, contents_2d) :
        for y, xy_content in zip(ybins_cen, y_slice):
            rownum += 1
            f.write(f"ROW:  {rownum:3d}   {x:10.5f}  {y:10.5f}    {xy_content}\n")
        
    f.close()
    
    return   # end dump_hist2d_contents 

def get_varname_dump(varname_orig):
    varname_dump = f"BINCENTER_{varname_orig}"

    # replace or remove special symbols that should not be used for
    # variables name
    # e.g.,
    #   BINCENTER_sqrt(q+3)  -> BINCENTER_sqrt_q+3
    #   BINCENTER_q*5        -> BINCENTER_sqrt_q_x_5

    replace_dict = {
        ' ' : ''  ,
        '(' : '_' ,
        ')' : ''  ,
        '*' : '_x_'
    }

    for ch0, ch1 in replace_dict.items():   
        varname_dump = varname_dump.replace(ch0,ch1)    
    
    return varname_dump

def get_info_plot1d(args, info_plot_dict):

    do_hist   = OPT_HIST   in args.OPT
    do_ov     = OPT_OV     in args.OPT
    do_chi2   = OPT_OVCHI2 in args.OPT
    do_nevt   = OPT_NEVT   in args.OPT
    do_avg    = OPT_AVG    in args.OPT  or  OPT_MEAN in args.OPT
    do_ratio  = OPT_RATIO  in args.OPT
    do_stddev = OPT_STDDEV in args.OPT
    
    numplot   = info_plot_dict['numplot']
    xbins     = info_plot_dict['xbins']
    xbins_cen = info_plot_dict['xbins_cen']    
    df_dict   = info_plot_dict['df_dict']
    
    df             = df_dict['df']
    wgtfun         = df_dict['wgtfun']
    name_legend    = df_dict['name_legend']
    plt_legend     = name_legend
    plt_text_dict  = None
    nxbins         = len(xbins)
    
    xval_list = xbins_cen
    xerr_list = None

    is_empty = len(df.x_plot_val) == 0
    
    if is_empty:
        yval_list = [ 0.0 ]  * nxbins  # allow for empty plot        
    else:
        yval_list = binned_statistic(df.x_plot_val, df.x_plot_val, 
                                     bins=xbins, statistic='count')[0]
    
    plt_histtype = 'step'  # default for HIST
    
    if do_hist :
        do_plot_errorbar = False
        do_plot_hist     = True
        if OPT_HISTFILL in args.OPT:
            plt_histtype = 'stepfilled'
    else:
        # default
        do_plot_errorbar = True
        do_plot_hist     = False
        
    do_plot_hist_ov1d     = False

    # apply option user-weight function (see @@WGTFUN arg)
    if wgtfun:
        wgtfun_user   = get_wgtfun_user(xbins_cen, wgtfun)
        yval_list    *= wgtfun_user

    nevt = np.sum(yval_list)  # sum before re-normalizing (for printing only)
    if is_empty:
        errl_list = [1] * nxbins
        erru_list = [1] * nxbins        
    else:
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
            do_plot_hist     = False            
    else:
        df0          = info_plot_dict['df0'] 
        yval0_list   = info_plot_dict['yval0_list']
        errl0_list   = info_plot_dict['errl0_list']
        erru0_list   = info_plot_dict['erru0_list']        
        name0_legend = info_plot_dict['name0_legend']         
        if do_ov:
            # re-scale overlay plot only if chi2 option is requested
            ov_scale = np.sum(yval0_list) / np.sum(yval_list)
            do_plot_errorbar   = False
            do_plot_hist       = False
            do_plot_hist_ov1d  = True
        elif do_ratio:
            # plot ratio of 2nd plot over 1st plot, or 3rd plot over 1st plot; etc..
            
            # check of df is a subset of df0 --> binomial errors
            is_subset = is_cid_subset(info_plot_dict, df, df0)
            ratio, errl_ratio, erru_ratio = \
                compute_ratio(yval_list,  errl_list,  erru_list,
                              yval0_list, errl0_list, erru0_list, is_subset ) 
            # load arrays that will be store in dictionary below              
            yval_list  = ratio
            errl_list  = errl_ratio
            erru_list  = erru_ratio
            ratio = len(df) / len(df0)
            plt_legend = f"<ratio> = {ratio:.2f}"
            if args.legend_list:
                legend = args.legend_list[numplot]
                plt_legend += f' for {legend} '                
            elif args.CUT:  # user cut without df.
                cut = args.CUT[numplot]
                plt_legend += f' for {cut} '
                
        else:
            pass

        if not is_empty:
            yval_list  *= ov_scale # normalize integral to match file 0 
            errl_list  *= ov_scale
            erru_list  *= ov_scale
        logging.info(f"Overlay {name_legend} scaled by {ov_scale:.3e}")

    # ---------------

    # check for overlay weights with chi2
    if numplot > 0 and do_ov:
        wgt_ov = [ ov_scale ] * len(df.x_plot_val)
        if wgtfun :
            wgtfun_user = get_wgtfun_user(df.x_plot_val,wgtfun) 
            wgt_ov      = np.multiply(wgt_ov,wgtfun_user)
            len_wgt_ov  = len(wgt_ov)

        # check option to compute and print chi2/dof info on plot
        # Froce min error =1 in chi2 calc so that it's ok to plot
        # error=0 for bins with zero events.
        sqdif = (yval0_list-yval_list)**2
        sqerr = np.maximum((yval0_list + yval_list*ov_scale*ov_scale), 1.0)
        chi2  = np.sum( sqdif / sqerr )
        ndof  = len(xbins) - 1 
        text_chi2 = f"chi2/dof = {chi2:.1f}/{ndof}"
        logging.info(f"{name0_legend}/{name_legend} {text_chi2}")
        if do_chi2:
            xmin = xbins[0];  xmax = xbins[-1]
            x_text = xmin + 0.6*(xmax-xmin)
            y_text = 0.9 * np.max(yval0_list)  # warning; fragile coord calc
            plt_text_dict = { 'x_text': x_text, 'y_text': y_text,
                              'text': text_chi2 }


    # - - - - - - -
    # remove bins with zero counts to avoid strange looking x-axis in plot
    suppress_zero_count_bins = True
    if suppress_zero_count_bins:
        plot_lists = [ [], [], [], [] ]
        for x, y, errl, erru in zip(xval_list, yval_list, errl_list, erru_list):
            for n, val in enumerate([x,y, errl, erru]):
                if y != 0:
                    plot_lists[n].append(val)
        xval_list = plot_lists[0]
        yval_list = plot_lists[1]
        errl_list = plot_lists[2]
        erru_list = plot_lists[3]  
    
        
    # - - - -
    yerr_list = [errl_list, erru_list]
    nevt    = nevt
    avg     = np.mean(df.x_plot_val)
    median  = np.median(df.x_plot_val)
    stddev  = np.std(df.x_plot_val)

    stat_dict = {  # value        add legend    format
        STAT_NAME_N       : [ int(nevt),  do_nevt,
                              stat_format(STAT_NAME_N,nevt)         ],
        STAT_NAME_AVG     : [ avg,        do_avg,
                              stat_format(STAT_NAME_AVG,avg)        ],
        STAT_NAME_MEDIAN  : [ median,     False,
                              stat_format(STAT_NAME_MEDIAN,median)  ],
        STAT_NAME_STDDEV  : [ stddev,     do_stddev,
                              stat_format(STAT_NAME_STDDEV,stddev)  ]
        # overflow/underflow ??
    }
    
    for str_stat, tmp_list in stat_dict.items():
        val           = tmp_list[0]
        add_to_legend = tmp_list[1]
        logging.info(f"\t {str_stat:8} value for {name_legend}:  {val:.3f}")
        if add_to_legend:
            fmt_legend  = tmp_list[2]                
            plt_legend += f'  {str_stat}={val:{fmt_legend}}'            

    # - - - - - -
    info_plot_dict['do_plot_errorbar']  = do_plot_errorbar
    info_plot_dict['do_plot_hist']      = do_plot_hist
    info_plot_dict['do_plot_hist_ov1d'] = do_plot_hist_ov1d 
    info_plot_dict['xval_list']         = xval_list
    info_plot_dict['yval_list']         = yval_list
    info_plot_dict['xerr_list']         = None
    info_plot_dict['yerr_list']         = yerr_list
    info_plot_dict['plt_histtype']      = plt_histtype
    info_plot_dict['plt_size']          = None
    info_plot_dict['plt_legend']        = plt_legend
    info_plot_dict['plt_text_dict']     = plt_text_dict
    info_plot_dict['wgt_ov']            = wgt_ov
     
    return  # end get_info_plot1d

def stat_format(stat_name, stat_value):

    # Created Dec 27 2024
    # return format string to use on plot.
    # Format can depend on name of variables (stat_name)
    # and value (stat_value).

    if stat_name == STAT_NAME_N :
        return 'd'

    if abs(stat_value) > 1000:  # e.g., peakmjd
        fmt = '.1f'        
    elif abs(stat_value) < 5:   # e.g., mean color, stretch, or resid
        fmt = '.3f'
    else:
        fmt = '.2f'
        
    return fmt


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

    do_hist      = OPT_HIST in args.OPT 
    do_nevt      = OPT_NEVT in args.OPT

    numplot      = info_plot_dict['numplot']
    df_dict      = info_plot_dict['df_dict']
    df           = df_dict['df']
    name_legend  = df_dict['name_legend']
    plt_legend   = name_legend

    nevt         = len(df)
    plt_size     = 5 / math.log10(max(10,nevt))  # dot size gets smaller with nevt 
    if do_nevt: plt_legend += f'  N={nevt}'
    
    xval_list    = df.x_plot_val
    yval_list    = df.y_plot_val
    xerr_list    = None
    yerr_list    = None

    if 'x_plot_err' in df:  xerr_list = df.x_plot_err
    if 'y_plot_err' in df:  yerr_list = df.y_plot_err    

    if do_hist:
        info_plot_dict['do_plot_errorbar']    = False
        info_plot_dict['do_plot_hist']        = True
    else:
        # default
        info_plot_dict['do_plot_errorbar']    = True
        info_plot_dict['do_plot_hist']        = False
    
    info_plot_dict['do_plot_hist_ov1d']   = False
    info_plot_dict['do_ov2d_binned_stat'] = True
    info_plot_dict['wgt_ov']              = None
    info_plot_dict['plt_size']      = plt_size
    info_plot_dict['plt_legend']    = plt_legend
    info_plot_dict['plt_text_dict'] = None
    info_plot_dict['plt_histtype']  = None    

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
            info_plot_dict['plt_legend'] = name_legend_ref + ' - ' + name_legend 
            # xxx mark delete info_plot_dict['plt_legend'] = f"{name_legend}-{name_legend_ref}"
        
        if args.DIFF == OPT_DIFF_CID :
            #need to do an inner join with each entry in dic, then plot the diff
            # (join logic thanks to Charlie Prior)
            varname_idrow = plot_info.varname_idrow # e.g., CID or GALID or ROW
            join   = df_ref.join(df.set_index(varname_idrow), on=varname_idrow,
                                 how='inner', lsuffix='_0', rsuffix='_1')
            xval_list = join.x_plot_val_0.values
            yval_list = join.y_plot_val_0.values - join.y_plot_val_1.values
            # xxx mark delete yval_list = join.y_plot_val_1.values-join.y_plot_val_0.values 

            info_plot_dict['xval_list']     = xval_list
            info_plot_dict['yval_list']     = yval_list
            if yerr_list is not None:
                yerr_list_0 = join.y_plot_err_0.values
                info_plot_dict['yerr_list']  = yerr_list_0
                
        elif args.DIFF == OPT_DIFF_ALL :

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
            y[np.isnan(y)]         = 0                        
            info_plot_dict['xval_list'] = info_plot_dict['xbins_cen']
            info_plot_dict['yval_list'] = y_ref - y
            info_plot_dict['yerr_list'] = None            
        
    return  # end  get_info_plot2d

def overlay2d_binned_stat(args, info_plot_dict, ovcolor ):

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
            y_count    = replace_zeros(y_count)  # avoid divide by zero
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

            sum_wgt = replace_zeros(sum_wgt) # avoid divide by zero
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
                     fmt=plt_marker, label=plt_legend, color=ovcolor, zorder=5 ) 
            
    return  # end of overlay2d_binned_stat

def replace_zeros(float_list):

    float_list_nozeros = []
    val_replace = 1.0e-20
    
    for val in float_list:
        if val == 0.0 : val = val_replace
        float_list_nozeros.append(val)
    
    return float_list_nozeros

def apply_plt_misc(args, plot_info, plt_text_dict):

    # Created July 21 2024 by R.Kessler
    # wrapper for lots of plt. calls that are repeated several times.
    
    OPT          = args.OPT    
    do_logy      = OPT_LOGY      in OPT
    do_grid      = OPT_GRID      in OPT
    do_diag_line = OPT_DIAG_LINE in OPT
    
    bounds_dict   = plot_info.bounds_dict
    axis_dict     = plot_info.axis_dict_list[0] # ?? fragile?
    
    custom_bounds = bounds_dict['custom']
    xmin          = bounds_dict['xmin']
    xmax          = bounds_dict['xmax']
    set_ylim      = axis_dict['set_ylim']
    
    if do_logy:
        plt.yscale("log")
        
    if do_grid:
        plt.grid(zorder=10)

    if custom_bounds:
        plt.xlim(xmin, xmax)
        if set_ylim:
            ymin  = bounds_dict['ymin']
            ymax  = bounds_dict['ymax']            
            plt.ylim(ymin, ymax)

    # - - - -
    ymin, ymax    = plt.ylim() # valid for auto or custom bounds

    if do_diag_line:
        x = np.linspace(xmin,xmax,100);  y = x
        plt.plot(x,y, zorder=10)
    
    fsize_label = 12 * args.FONTSIZE_SCALE
    xlabel  = axis_dict['xaxis_label']
    ylabel  = axis_dict['yaxis_label']    
    plt.xlabel(xlabel, fontsize=fsize_label)  
    if ylabel is not None:
        plt.ylabel(ylabel, fontsize=fsize_label)

    # adjust tick label sizes (Dec 15 2024)
    fsize_ticklabel = 10*args.FONTSIZE_SCALE
    plt.xticks(fontsize=fsize_ticklabel)  # numbers under x-axis tick marks
    plt.yticks(fontsize=fsize_ticklabel)  # numbers left of y-axis tick marks
    
    if args.LEGEND_SIDE:
        # push legend outside of box on right side
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5) )
    elif args.LEGEND == SUPPRESS:
        pass  # no legend
    else:
        # default; let matplotlib find best place for legend
        fsize_legend = 10 * args.FONTSIZE_SCALE
        plt.legend(loc=None, bbox_to_anchor=None, fontsize = fsize_legend)

    len_title = len(args.TITLE)
    fsize_title   = 14 * args.FONTSIZE_SCALE
    if len_title > 50:
        fsize_title *= (50.0/len_title)**0.9 # reduce font size to fit on page
        
    plt.title(args.TITLE, fontsize=fsize_title)

    # auto-generated text; eg., chi2/dof for chi2 option
    fsize_text = 14 * args.FONTSIZE_SCALE
    if plt_text_dict is not None:
        x_text = plt_text_dict['x_text']
        y_text = plt_text_dict['y_text']
        text   = plt_text_dict['text']
        plt.text(x_text, y_text, text, fontsize=fsize_text )

    # optional display text from user input @@TEXT
    if args.text_list:
        for [ x_rel, y_rel, text ] in args.text_list:
            x     = xmin + x_rel*(xmax-xmin)
            y     = ymin + y_rel*(ymax-ymin)
            plt.text(x, y, text, fontsize=fsize_text)
        
        
    return
    # end apply_plt_misc
    

def get_wgtfun_user(xcen,wgtfun):

    n_val     = len(xcen)

    if wgtfun == '1' :
        wgt_vals  = [1.0] * n_val   # default weights are 1
    elif wgtfun:

        # create numpy array x (from xcen) and x varname must be used to
        # match variable in user-define WGTFUN (e.g., 1-x).
        x   = np.array(xcen)
        wgt_plus_np = wgtfun
        for str_orig, str_final in NUMPY_FUNC_DICT.items(): 
            wgt_plus_np = wgt_plus_np.replace(str_orig,str_final)

        wgt_vals = eval(wgt_plus_np)
        
    return wgt_vals

def print_cid_list(df, name_legend) :
    
    # print list of cids to stdout
    varname_idrow = plot_info.varname_idrow # e.g., CID or GALID or ROW 
    id_list = sorted(df[varname_idrow].to_numpy())[0:NMAX_CID_LIST]
    print(f"\n {varname_idrow}s passing cuts for '{name_legend}' : \n{id_list}" )
    sys.stdout.flush()
    return


# ===================================================
#   Add main, June 2024
# ========================================
if __name__ == "__main__":

    setup_logging()
    logging.info("# ========== BEGIN plot_table.py  ===============")
    logging.info(f"Full command: {' '.join(sys.argv)}")
    
    args = get_args()
                              
    # prepare input args
    arg_prep_driver(args)

    # translate @@VARIABLE, @@CUT, @@ERROR to dataframe syntaix
    translate_driver(args)
    
    
    plot_info = Namespace()  # someplace to store internally computed info
    
    # set axis info (var names, labels, ...)
    plot_info.axis_dict_list = []
    for var in  args.var_list: 
        set_axis_dict(args, plot_info, var)

    # set custom bounds if user passes @@BOUNDS
    set_custom_bounds_dict(args,plot_info)

    # set args for hist2d if 2d histogram is requested
    set_hist2d_args(args, plot_info)

    read_tables(args, plot_info)  # read each input file and store data frames

    # set x bins based on custom bounds or auto bounds
    set_xbins(args, plot_info)

    # set axis labels after reading tables and setting xbins
    set_axis_labels(args, plot_info)
    
    # initialize matplotlib figure
    plt.figure( figsize=(FIGSIZE[0]*args.XSIZE_SCALE,
                         FIGSIZE[1]*args.YSIZE_SCALE) ) 


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

    logging.info('Done making plots')
    
    # === END: ====
    

    
