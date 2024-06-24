#!/usr/bin/env python
#
# Created by B.Popovic during his graduate career at Duke University.
# Installed into SNANA Jun 24 2024 by R.Kessler
#
# 
import os
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

import argparse
from argparse import RawTextHelpFormatter
parser=argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, prefix_chars='@')
from scipy.stats import binned_statistic
import distutils.util 

"""
How do I use this plotter? Well, first you're going to need to present your M0DIF/FITRES/HOSTLIB files!

Then you'll need to specify what you want to plot, eg, @@VARIABLE. For this, you can do a 1D histogram or a 2D function. This is formatted with pandas dataframes, so you'll need to be familiar with that syntax. For instance, if you want to see the redshift distribution, df.zHD or df.zHD.values are both acceptable! Need to use df, this is very important.

You can also do 2D functions. These are parsed by using a colon (:) to separate the x and y variables. The syntax is generally x:y
   For instance, 'df.zHD.values:df.x1.values' will show z vs x1. You can also do algebra! 'df.zHD.values:df.mB.values - 3.1*df.c.values + 0.16*df.x1.values' will give you the Tripp estimator! 


Then you can also set custom boundaries and bins! This is done by giving '@@BOUNDS minimum maximum binsize'. For 2D plots, same colon separation deal. 

You can save your figure using @@SAVE.

If you want to compare and contrast values for the same CID, you enable the @@DIFF option! There are two valid options. 'ALL' will compare the difference in median values. 'CID' will compare the difference for the same CIDs. 

CUT is for when you want to place some sort of cut on the sample. This is usually a 'df.loc[]' option of some capacity. Please note that this NEEDS to have quotation marks around it. This is a problem with bash and there's nothing I can do about that.

Finally, ALPHA exists to tinker with the alpha values. Pretty simple! 

Anyways, an example of the usage is like so:


plotter-class.py @@FITRES File1.FITRES File2.FITRES @@VARIABLE df.zHD.values:df.mB.values - 3.1*df.c.values + 0.16*df.x1.values - df.MU.values @@CUT "df.loc[df.IDSURVEY < 15]"

That'll plot z vs Tripp Estimator for your two files!

"""

parser.add_argument('@@FITRES', help='The location of your FITRES file that you need plotted. You can give a space delineated list of different FITRES files. Please make sure they are named differently, as their names (not full filepath) will be used for labeling.', nargs="+")
parser.add_argument('@@VARIABLE', help="""The variable you want to plot. Needs to be a valid FITRES parameter. \n
If you give only one value, this will generate a histogram. If you give two colon delimited values, such as df.zHD:df.mB it will plot zHD (x) vs mB (y). Please note that for the histogram, counts will be normalised to the first file given.""", nargs="+")
parser.add_argument('@@BOUNDS', default='loose', help=
"""loose (default): bounds are maximum and minimum values of specified parameter. \n 
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

args = parser.parse_args()
VARIABLE = args.VARIABLE
FILENAME = args.FITRES
BOUNDS = args.BOUNDS
FORMAT = args.SAVE
DIFF = args.DIFF
ALPHA = args.ALPHA
CUT = args.CUT
NROWS = args.NROWS

if DIFF: DIFF = DIFF.strip()

VARIABLE = ''.join([str(elem) for elem in VARIABLE])
if CUT: CUT = ''.join([str(elem) for elem in CUT])


def NAndR(filename):
    with open(filename) as fp:   
        for i, line in enumerate(fp):     
            if line.startswith('VARNAMES:'):         
                line = line.replace(',',' ')   
                line = line.replace('\n','')    
                Names = line.split()     
            elif (line.startswith('SN')) or (line.startswith('ROW')) or (line.startswith('GAL:')):     
                Startrow = i 
                break   
    return Names, Startrow  


def NAndRzip(filename):
    import gzip
    with gzip.open(filename) as fp:
        for i, bline in enumerate(fp):
            line = bline.decode("utf-8")
            if line.startswith('VARNAMES:'):       
                line = line.replace(',',' ')            
                line = line.replace('\n','')         
                Names = line.split()           
            elif (line.startswith('SN')) or (line.startswith('ROW')) or (line.startswith('GAL:')):                
                Startrow = i            
                break       
    return Names, Startrow  

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

filenames = [l.split("/")[-1] for l in FILENAME]
if (any(filenames.count(x) > 1 for x in filenames)):
    "flag this"

MASTERLIST = {}

def plotter_func(dictionary, DIFF, plotdic, boundsdic, CUT, ALPHA):
    if custom_bounds:                       
        bins = np.arange(boundsdic['x'][0],boundsdic['x'][1],boundsdic['x'][2])                                                           
        plt.xlim([boundsdic['x'][0], boundsdic['x'][1]])  
    else:                                   
        bins = np.linspace(boundsdic[min(boundsdic, key=boundsdic.get)], boundsdic[max(boundsdic, key=boundsdic.get)], 30)            
    if len(plotdic) == 1:                           
        if (DIFF == 'ALL') or (DIFF == 'CID'):
            print("The DIFF feature does not work for histograms. Quitting to avoid confusion.")
            quit()

        for n,k in enumerate(dictionary):           
            k = dictionary[k]                       
            print('The upper and lower bounds are:', np.around(bins[0],4), 'and', np.around(bins[-1],4), 'respectively')    
            if n == 1:                              
                sb = binned_statistic(k.x_plot_val, k.x_plot_val, bins=bins, statistic='count')[0] #Get counts                                
                errl,erru = poisson_interval(np.sum(db)*sb/np.sum(sb)) 
                sb = np.sum(db)*sb/np.sum(sb)
                plt.errorbar((bins[1:] + bins[:-1])/2., sb, yerr=[sb-errl, erru-sb], label=k.name.values[0], fmt='o')                    
            else:                                   
                db = binned_statistic(k.x_plot_val, k.x_plot_val, bins=bins, statistic='count')[0] #Get counts                                
                errl,erru = poisson_interval(db) #And error for those counts                                                
                plt.errorbar((bins[1:] + bins[:-1])/2., db, label=k.name.values[0], yerr=[db-errl, erru-db], fmt='o')                 
            print("The Mean value for ", k.name.values[0], " is:", np.mean(k.x_plot_val))                                                      
            print("The Median value for ", k.name.values[0], " is:", np.median(k.x_plot_val))                                                  
            print("The standard deviation value for ", k.name.values[0], " is:", np.std(k.x_plot_val))                                         
        plt.xlabel(plotdic['x'])  #In this case, VAR = [string], so we're going to strip the list.                     
        plt.legend()                                
        if CUT: plt.title(CUT)
    elif (DIFF == 'CID') or (DIFF == 'ALL'):
        print("plotting DIFF option now")
        if len(plotdic) == 2: 
            try:         
                plt.ylim([boundsdic['y'][0], boundsdic['y'][1]])            
            except KeyError:     
                pass  
        keylist = list(dictionary.keys())
        differator = dictionary[keylist[0]]
        for k in keylist[1:]:
            k = dictionary[k]
            if DIFF == 'CID':
                #need to do an inner join with each entry in dic, then plot the diff
                join = differator.join(k.set_index('CID'), on='CID', how='inner', lsuffix='_1', rsuffix='_2') #Thank you to Charlie for this
                plt.scatter(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values, alpha=ALPHA, label='Diff')
                avgdiff = binned_statistic(join.x_plot_val_1.values, join.y_plot_val_1.values - join.y_plot_val_2.values, bins=bins, statistic='median')[0]                                     
                plt.scatter((bins[1:] + bins[:-1])/2, avgdiff, label="Mean Difference", color='k')
            elif (DIFF == 'ALL'):
                try:
                    plt.scatter(differator.x_plot_val, differator.y_plot_val - k.y_plot_val, label=differator.name.values[0]+" - "+k.name.values[0], alpha=ALPHA)
                except ValueError:
                    pass
                median_differator = binned_statistic(differator.x_plot_val, differator.y_plot_val, bins=bins, statistic='median')[0]
                median = binned_statistic(k.x_plot_val, k.y_plot_val, bins=bins, statistic='median')[0]
                plt.scatter((bins[1:] + bins[:-1])/2., median_differator - median, 
                            label=differator.name.values[0]+" - "+k.name.values[0]+" median", marker="^", zorder=10)
            else:  
                print("You gave "+str(DIFF))                               
                print("That is not a valid DIFF option. Quitting.")               
                quit()  
            plt.xlabel(plotdic['x'])                    
            plt.ylabel(plotdic['y']+" diff")                    
            plt.legend()                                
            if CUT: plt.title(CUT)       
    else:    
        try:
            plt.ylim([boundsdic['y'][0], boundsdic['y'][1]]) 
        except KeyError:
            pass
        for k in dictionary:                        
            k = dictionary[k]                       
            plt.scatter(k.x_plot_val, k.y_plot_val, alpha=ALPHA, label=k.name.values[0], zorder=0)            
            median = binned_statistic(k.x_plot_val, k.y_plot_val, bins=bins, statistic='median')[0] #Get counts                
            plt.scatter((bins[1:] + bins[:-1])/2., median, label=k.name.values[0]+" median", marker="^", zorder=10)           
        plt.xlabel(plotdic['x'])                    
        plt.ylabel(plotdic['y'])                    
        plt.legend()                                
        if CUT: plt.title(CUT) 
    return 

plotdic = {}
boundsdic = {} 


for n,VAR in enumerate(VARIABLE.split(":")):
    if n == 0:
        plotdic['x'] = str(VAR) 
    else:
        plotdic['y'] = str(VAR) 


for l in (FILENAME): 
    print("Loading ",l.split("/")[-1], "...") #Inform that we're loading the first file.
    try:
        Names1, StartRow1 = NAndR(l) #Get info on where to start reading the file
    except (FileNotFoundError, NameError, UnicodeDecodeError):
        if UnicodeDecodeError:
            Names1, StartRow1 = NAndRzip(l)
        else:
            print('Could not find the FITRES you specified!')
            print("You were pointing to: ", FILENAME)
            sys.stdout.flush() # "dad! dad! look what got caught in the snare!" "good work, timmy, its AttributeError for dinner tonight" -Ross
            quit() #Quits if one or more files is missing
    if NROWS != 0: df = pd.read_csv(l, header=None, skiprows=StartRow1,names=Names1, sep=r"\s+", skip_blank_lines=True, comment='#', nrows=NROWS)
    else: df = pd.read_csv(l, header=None, skiprows=StartRow1,names=Names1, sep=r"\s+", skip_blank_lines=True,  comment='#')
    try:
        df['CID'] = df['CID'].astype(str)
    except KeyError:
        print("No CIDs present in this file. Making note of that here.")

    if CUT: df = eval(CUT)
    #Will need to do all the x and y assignments here if we use this setup 
    keycount = -1
    keyname = l.split("/")[keycount]
    while keyname in list(MASTERLIST.keys()):
        keycount -= 1
        keyname = l.split("/")[keycount]

    MASTERLIST[keyname] = df
    MASTERLIST[keyname]['name'] = keyname
    try:
        MASTERLIST[keyname]['x_plot_val'] = eval(plotdic['x'])
        boundsdic[keyname+"_min"] = np.amin(MASTERLIST[keyname]['x_plot_val'])
        boundsdic[keyname+"_max"] = np.amax(MASTERLIST[keyname]['x_plot_val'])
        if len(plotdic) == 2:
            MASTERLIST[keyname]['y_plot_val'] = eval(plotdic['y'])
    except AttributeError:
        print("Couldn't process this command! One of the things you are trying to plot is not present in one or more of the files!")
        quit() # "oh the sweet turgid flesh of access discrepancy" - Ross 
    #need to comb through existing keys in dictionary and make sure there are no overlaps
    print("Done loading that file!")

print("Done loading all files!")

#Now I've got all the things loaded in as a dataframe rather than a class.

if (any(list(MASTERLIST.keys()).count(x) > 1 for x in list(MASTERLIST.keys()))):
    print("Your filenames are identical and the directory above them also has the same name. Quitting...")
    quit()

custom_bounds = False                                        
if (len(BOUNDS) != 5):
    BOUNDS = ' '.join([str(elem) for elem in BOUNDS])         
    custom_bounds = True                                                 
    for n,BND in enumerate(BOUNDS.split(':')):           
        if n == 0:     
            boundsdic['x'] = [float(i) for i in BND.split()]            
        else:                                                   
            boundsdic['y'] = [float(i) for i in BND.split()]


# if DIFF == True: #Thank you to Charlie for this part to ensure that join has only shared CIDs in it! 
#     if FILENAME[0].endswith('M0DIF'): 
#         pass
#     elif ALPHA == 0:
#         pass
#     else:
#         df1['CID'] = df1['CID'].astype(str)
#         df2['CID'] = df2['CID'].astype(str)
#         try:
#             join = df1.join(df2.set_index('CID'), on='CID', how='inner', lsuffix='_df1', rsuffix='_df2') #creates a single shared fitres file for use later
#         except NameError:
#             print("For DIFF to work, you need to give me two files! You only gave one. Quitting...")
#             quit()


#now to rewrite the BOUNDS bits.

print("Plotting now!")


plt.figure()
plotter_func(MASTERLIST, DIFF, plotdic, boundsdic, CUT, ALPHA)
if FORMAT !="None":                  
    plt.savefig(FORMAT, bbox_inches="tight", format=FORMAT.split(".")[-1])                                             
plt.show()         
