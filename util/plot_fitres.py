#!/usr/bin/env python
"""
Script for Comparative Visualization of FITRES Data Files

This script is for visualization of data contained within (any two) FITRES files. It supports 
plotting various data metrics, including scatter plots and histograms, with options for binning, scaling and styling. The script allows for the comparison between two datasets, with features to scale the second 
dataset automatically or according to a user-specified factor.


Features:
- Supports reading FITRES files while ignoring comment lines.
- Identifies numeric columns for potential plotting.
- Customizable plot aesthetics including figure size, tick style, and log scaling.
- Comparative plotting with automatic or manual scaling of the second dataset.
- Outputs plots to a single PDF file for easy sharing and viewing.

Usage:
The script is executed from the command line, with various options to specify the input files, output file, 
and plotting parameters. Example usage:

python script_name.py <path_to_fitres_file1> -I <path_to_fitres_file2> -O <output_path.pdf> -l1 "Dataset 1" -l2 "Dataset 2" -c1 <column_name> -c2 <column_name> -p scatter -S auto

Example:
python compare_fitres.py /global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SURVEYS/LSST/USERS/kessler/debug/plot_fitres/FITOPT000_BCOR_PRESCALE10.FITRES.gz --input2 /global/cfs/cdirs/lsst/groups/TD/SN/SNANA/SURVEYS/LSST/USERS/kessler/debug/plot_fitres/FITOPT000_BCOR.FITRES.gz -O ~/compare_fitres_plot_test.pdf  -l1 "DATA" -l2 "SIM" -c1 "c"  -p "scatter" -S

Arguments:
- `file_path`: Path to the first FITRES file (mandatory).
- `-I` / `--input2`: Path to the second FITRES file (optional).
- `-O` / `--output`: Path for the output PDF file (optional, default: current directory).
- `-l1` / `--label1`: Label for the first dataset (optional).
- `-l2` / `--label2`: Label for the second dataset (optional).
- `-c1` / `--custom_col1`: Custom column name for the x-axis (optional).
- `-c2` / `--custom_col2`: Custom column name for the y-axis (optional).
- `-p` / `--plot_type`: Type of plot: "scatter" or "histogram" (optional).
- `-S` / `--scale2`: Apply scaling to the second dataset (optional).

Note:
- This script suppresses warnings to improve output readability. This behavior can be modified as needed.


Ayan Mitra
2024
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import argparse
from matplotlib.ticker import ScalarFormatter
import sys
import warnings
warnings.filterwarnings("ignore")

column_pairs = [
    {'columns': ('zHD', 'MU'), 'xlim': (0, 1.5), 'ylim': (34, 46), 'color': 'blue', 'logscale': None},
    {'columns': ('zHD', 'MUERR'), 'xlim': (0, 1.5), 'ylim': (0.04, 30), 'color': 'green', 'logscale': None},
    {'columns': ('zHD', 'biasCor_mu'), 'xlim': (0, 1.5), 'ylim': (10, 20), 'color': 'red', 'logscale': None}, # 'x'
    {'columns': ('zHD', 'SIM_DLMAG'),  'color': 'red'} #'xy'
]


# For histograms
individual_columns = ['zHD', 'SIM_TEMPLATE_INDEX']

def read_fitres_file(file_path):
    data = pd.read_csv(file_path, comment='#', delim_whitespace=True)
    return data




def find_suitable_plot_columns(data):
    numeric_cols = data.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) < 2:
        print("Error: Not enough numeric columns for plotting.")
        return None, None
    return numeric_cols[0], numeric_cols[1]  # Adjust as needed



def set_plot_aesthetics(ax):
    plt.rcParams['figure.figsize'] = [8, 4]
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in', labelsize=26)
    ax.tick_params(axis='both', which='minor', length=4)  # Minor ticks
    ax.tick_params(axis='both', which='major', length=6)  # Major ticks
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

def plot_data(data, x_col, y_col, label=None, ax=None, plot_aesthetics=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=[8, 4], dpi=150)  # Create a new figure if ax is not provided

    if plot_aesthetics is not None:
        ax.set_xlim(plot_aesthetics.get('xlim', (data[x_col].min(), data[x_col].max())))
        ax.set_ylim(plot_aesthetics.get('ylim', (data[y_col].min(), data[y_col].max())))
        color = plot_aesthetics.get('color', 'b')  # Default to blue if not specified
        
        # Check if logscale is not None before checking 'x' or 'y' in logscale
        logscale = plot_aesthetics.get('logscale', '')
        if logscale:  # This will be False if logscale is None or an empty string
            if 'x' in logscale:
                ax.set_xscale('log')
            if 'y' in logscale:
                ax.set_yscale('log')

    try:
        data.plot(x=x_col, y=y_col, kind='scatter', color=color, ax=ax)
        set_plot_aesthetics(ax)

        if label:
            ax.set_title(label)
        ax.set_xlabel(x_col, fontsize=18)
        ax.set_ylabel(y_col, fontsize=18)

    except Exception as e:
        print(f"Error in plotting data: {e}")


def plot_binned_data(data, x_col, y_col, label=None, ax=None, x_bins=10, y_bins=10, plot_aesthetics=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=[8, 4], dpi=150)
        own_figure = True
    else:
        own_figure = False

    try:
        # Clean data to remove NaNs or infinite values
        data = data.dropna(subset=[x_col, y_col])
        data = data[np.isfinite(data[x_col]) & np.isfinite(data[y_col])]

        # Define bin edges
        x_bin_edges = np.linspace(data[x_col].min(), data[x_col].max(), x_bins + 1)
        y_bin_edges = np.linspace(data[y_col].min(), data[y_col].max(), y_bins + 1)

        # Assign data points to bins
        data['x_bin'] = pd.cut(data[x_col], bins=x_bin_edges, labels=False, include_lowest=True)
        data['y_bin'] = pd.cut(data[y_col], bins=y_bin_edges, labels=False, include_lowest=True)

        # Group by bins and calculate mean
        binned_data = data.groupby(['x_bin', 'y_bin']).agg({x_col: 'mean', y_col: 'mean'})

        # Calculate standard error separately
        std_err = data.groupby(['x_bin', 'y_bin']).agg({x_col: lambda x: np.std(x) / np.sqrt(len(x)),
                                                        y_col: lambda y: np.std(y) / np.sqrt(len(y))})
        binned_data[x_col + '_std'] = std_err[x_col]
        binned_data[y_col + '_std'] = std_err[y_col]

        # Remove rows with NaNs after aggregation
        binned_data = binned_data.dropna()

        if binned_data.empty:
            print("No data to plot after binning.")
            return

        # Apply custom aesthetics
        if plot_aesthetics is not None:
            if 'xlim' in plot_aesthetics:
                ax.set_xlim(plot_aesthetics['xlim'])
            if 'ylim' in plot_aesthetics:
                ax.set_ylim(plot_aesthetics['ylim'])
            color = plot_aesthetics.get('color', 'r')  # Default color to red if not specified
            logscale = plot_aesthetics.get('logscale', '')
            if 'x' in logscale:
                ax.set_xscale('log')
            else:
                ax.set_xscale('linear')
            if 'y' in logscale:
                ax.set_yscale('log')
            else:
                ax.set_yscale('linear')
        else:
            color = 'r'  # Default color to red if plot_aesthetics is not provided

        # Plotting
        ax.errorbar(binned_data[x_col], binned_data[y_col],
                    xerr=binned_data[x_col + '_std'], yerr=binned_data[y_col + '_std'],
                    fmt='o', ecolor='gray', color=color, capsize=3)

        set_plot_aesthetics(ax)

        if label:
            ax.set_title(label)
        ax.set_xlabel(x_col, fontsize=18)
        ax.set_ylabel(y_col, fontsize=18)

        # Show the plot if we created our own figure
        if own_figure:
            plt.show()

    except Exception as e:
        print(f"Error in plotting data: {e}")
        
def plot_data_for_each_dataset(data1, data2, args, axes):
    # Plot for data1
    ax1 = axes[0]
    plot_custom_data(data1, args, ax1)

    # Plot for data2
    ax2 = axes[1]
    args.label1 = args.label2  # Adjust label for data2
    plot_custom_data(data2, args, ax2)
    #plt.legend()
    set_plot_aesthetics(ax1)
    set_plot_aesthetics(ax2)


def plot_data_for_single_dataset(data1, args, axes):
    # Plot for data1
    ax1 = axes[0]
    plot_custom_data(data1, args, ax1)

    set_plot_aesthetics(ax1)

    
def plot_custom_data(data, args, ax):
    if args.custom_col2 and args.custom_col2 in data.columns:
        # If two input parameters are passed, make a scatter plot
        plot_data(data, args.custom_col1, args.custom_col2, label=f"{args.label1}, ' '", ax=ax)
        ax.legend()
    # If only 1 parameter is passed, there is a choice b/w
    # histogram or scatter plot for each
    elif args.plot_type == 'histogram':
        data[args.custom_col1].plot(kind='hist', ax=ax,color='red',bins = 15, histtype='step',lw=3)
        ax.set_title(f"{args.label1}")
        ax.set_xlabel(args.custom_col1, fontsize=16)
        ax.set_ylabel('Frequency', fontsize=16)
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    else:
        num_bins = 10
        # Use np.histogram to get the bin counts and edges
        counts, bin_edges = np.histogram(data[args.custom_col1], bins=num_bins)
        # Calculate bin centers from edges for plotting purposes
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        
        # Step 2: Plot the binned data
        # Now, ax.plot() plots bin_centers vs counts
        ax.plot(bin_centers, counts, linestyle='-', marker='o', color='red')  # Line plot with markers
        ax.set_ylabel('Counts', fontsize=16)  # Set the x-axis label
        ax.set_xlabel(f"Binned {args.custom_col1}", fontsize=16)  # Set the y-axis label to the name of custom_col1
        ax.set_title(f"{args.label1}")  # Set the title of the plot
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #ax.legend()
        
def main():
    parser = argparse.ArgumentParser(description='Plot data from FITRES files and save to a PDF.')
    parser.add_argument('file_path', type=str, help='Path to the first FITRES file.')
    parser.add_argument('-I', '--input2', type=str, help='Path to the second FITRES file (optional).')
    parser.add_argument('-O', '--output', type=str, default=os.path.join(os.getcwd(), 'compare_fitres_plot_default.pdf'),
                        help='Path for the output PDF file (optional).')
    parser.add_argument('-l1', '--label1', type=str, default="Original", help='Label for the first file (optional).')
    parser.add_argument('-l2', '--label2', type=str, default="Reference", help='Label for the second file (optional).')
    parser.add_argument('-c1', '--custom_col1', type=str, help='Custom column name for the x-axis (optional).')
    parser.add_argument('-c2', '--custom_col2', type=str, help='Custom column name for the y-axis (optional).')
    parser.add_argument('-p', '--plot_type', type=str, choices=['scatter', 'histogram'], default='scatter', help='Type of plot: "scatter" or "histogram" (optional).')
    parser.add_argument('-S', '--scale2', nargs='?', const='auto', default=1, help='Apply scaling. If no value is provided, scale factor is calculated as N_row Data1/N_row Data2. A specific value can also be provided (optional), which will override all other values.')


    args = parser.parse_args()
    try:
        data1 = read_fitres_file(args.file_path.strip())
        data2 = None
        factor = 1  
        if args.input2:
            data2 = read_fitres_file(args.input2.strip())
            print('Shape of data1 ', np.shape(data1))
            print('Shape of data2 ', np.shape(data2))
            print('** 1 row',data1.iloc[3],'\n',data2.iloc[3])
            
            if args.scale2 == 'auto':
                # Compute factor only if scaling is enabled
                print('Using Auto scaling')
                nume, deno = sorted((df.shape[0] for df in (data1, data2)))
                factor = nume / deno
                print("** Scale Factor %s/%s = %s"%(nume, deno,factor))
            else:
                # Use the user-provided factor
                print('Using User input scaling')
                factor = 1/float(args.scale2)
            print("Factor:", factor)
            
            cols1 = data1.columns
            cols2 = data2.columns
            common_cols = set(cols1).intersection(set(cols2))
            unique_cols1 = set(cols1).difference(set(cols2))                          
            unique_cols2 = set(cols2).difference(set(cols1))                           
            total_unique = len(unique_cols1) + len(unique_cols2)                          
            total_cols = len(set(cols1).union(set(cols2)))                             
            print(f"Common columns: {len(common_cols)}")                                                        
            print(f"Unique columns in first file: {len(unique_cols1)} {unique_cols1}")                          
            print(f"Unique columns in second file: {len(unique_cols2)} {unique_cols2}")
            print(f"Common columns \n\n {common_cols} \n\n")
            if total_unique > total_cols * 0.5:
                print("More than 50% columns are different. Aborting.")
                sys.exit()


    except FileNotFoundError:
        print("File not found")
        return
    except  Exception as e:
        print(f"An error occured: {e}")
        return

    with PdfPages(args.output) as pdf:
        # Plot default column pairs

        for pair in column_pairs:
            x_col, y_col = pair['columns']  # Extract the column names from the dictionary
            plot_aesthetics = {k: v for k, v in pair.items() if k != 'columns'}  # Extract additional plotting aesthetics
            
            if x_col in data1.columns and y_col in data1.columns:
                # Create a figure for the pair
                fig, axes = plt.subplots(nrows=1, ncols=2, figsize=[8, 4], dpi=150)
                if data2 is not None and x_col in data2.columns and y_col in data2.columns:
                    # Plot for data1 on the first axis
                    plot_binned_data(data1, x_col, y_col, label=f"{args.label1} ", ax=axes[0], plot_aesthetics=plot_aesthetics)
                    # Plot for data2 on the second axis
                    plot_binned_data(data2, x_col, y_col, label=f"{args.label2}", ax=axes[1], plot_aesthetics=plot_aesthetics)
                else:
                    # If data2 is not provided or does not have the columns, plot only data1
                    plot_binned_data(data1, x_col, y_col, label=f"{args.label1} ", ax=axes[0], plot_aesthetics=plot_aesthetics)
                    # Hide the second subplot if not used
                    axes[1].set_visible(False)

                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                

        for column in individual_columns:
            fig, axes = None, None
            if data2 is not None and column in data2.columns:
                # If data2 is provided, plot side by side
                fig, ax = plt.subplots(dpi=150)
                if column in data1.columns:
                    #ax1 = axes[0]
                    counts, bin_edges = np.histogram(data1[column], bins=15)
                    ax.bar(bin_edges[:-1], counts, width=np.diff(bin_edges), 
                            edgecolor='red',facecolor='none', align='edge', lw=1,label=f"{args.label1}")
                    ax.set_ylabel('Frequency', fontsize=16)
                    ax.set_xlabel(column, fontsize=16)
                    set_plot_aesthetics(ax)
                if column in data2.columns:
                    #ax2 = axes[1]
                    # Calculate histogram for data2 without scaling data values
                    counts, bin_edges = np.histogram(data2[column], bins=15)
                    # Apply factor to the frequencies
                    scaled_counts = counts * factor
                    # Plot the scaled histogram
                    ax.bar(bin_edges[:-1], scaled_counts, width=np.diff(bin_edges), 
                            edgecolor='black',alpha=0.5, facecolor='none',
                           align='edge', lw=1, label=f'{args.label2}, SCALE={factor}')
                    ax.legend(fontsize=10)
            else:
                # If only data1 is provided, plot single column
                fig, ax = plt.subplots(dpi=150)
                counts, bin_edges = np.histogram(data1[column], bins=15)
                ax.bar(bin_edges[:-1], counts, width=np.diff(bin_edges),
                       edgecolor='black',facecolor='none', align='edge', lw=3)
                ax.set_title(f"{args.label1}, Index vs {column}")
                ax.set_xlabel(column, fontsize=18)
                ax.set_ylabel('Frequency', fontsize=18)
                set_plot_aesthetics(ax)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)


            

        # Custom column plotting
        if args.custom_col1:
            fig, axes = None, None
            if data2 is not None and args.custom_col1 in data2.columns:
                # If data2 is provided, plot side by side
                fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), dpi=150)
                plot_data_for_each_dataset(data1, data2, args, axes)
            else:
                # If only data1 is provided, plot single column
                fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150)
                axes = [ax]
                plot_data_for_single_dataset(data1,  args, axes)

            #plot_data_for_each_dataset(data1, data2, args, axes)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    print(f"Plots saved to {args.output}")

    
if __name__ == "__main__":
    main()
