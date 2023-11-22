#!/usr/bin/env python
"""A diagnostic tool that takes SNANA .DUMP files and plots PEAKMAG* vs ZCMB.
It will iterate over all PEAKMAG* columns in the .DUMP file and put them in 
a single multiplot output.

Maintained by Benjamin Rose <Ben_Rose@baylor.edu>
"""

history = """# History

* Nov 21, 2023 - v1.1.1 - fix figure labels so they can't start with "_"
* Oct 31, 2023 - v1.1 - `files` & `--reference` can now take paths that contain one DUMP file, BMR
* Oct 27, 2023 - no longer require --reference, BMR
* Oct 27, 2023 - initial verison, Benjamin Rose
"""

import argparse
from pathlib import Path
import re
from sys import exit

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

__version__ = "1.1.1"
SNANA_NA_VALUES = [-9, -99, 99]

parser = argparse.ArgumentParser(prog="mag-redshift.py", description=__doc__, 
								 formatter_class=argparse.MetavarTypeHelpFormatter)
parser.add_argument("files", type=Path, nargs="+",
					help="SNANA DUMP files to plot. Skips an argument if it is also passed as --reference.\
						  Can be paths, if DUMP files are one per path.")
parser.add_argument("--reference", type=Path, 
					help="DUMP file to use as a plotting reference. Can be a path, if it contains only one DUMP files.")
parser.add_argument("--id", type=str, 
					help="reference string to be used in the figure title and filename (default: DUMP file's parent folder name)")
parser.add_argument("--output_format", default="png", type=str, 
					help="file format for resulting plots (default: %(default)s)")
parser.add_argument("--zp", default=40, type=float, 
					help="plotting value for non-detections (default: %(default)s)")
parser.add_argument("--n_columns", default=4, type=int, 
					help="number of subplot columns (default: %(default)s)")
parser.add_argument("--history", action="store_true", help="print out development hisotry and exit")
parser.add_argument("--version", action="version", version='%(prog)s ' + __version__)
args = parser.parse_args()

if args.history:
	print(history)
	exit()

files = args.files
output_format = args.output_format
run_id = args.id
no_detection_mag = args.zp
columns = args.n_columns
	
if args.reference is not None:
	if args.reference in files:
		files.remove(args.reference)
	
	if args.reference.is_file():
		reference_file = args.reference
	else:
		reference_file = next(args.reference.glob("*.DUMP"))  #case_sensitive=False in python 3.12

if run_id is None:
	NO_ID = True
else:
	NO_ID = False


# import SNIa data
if args.reference is not None:
	with open(reference_file, "r") as f:
		sn_data = pd.read_csv(f, sep="\s+", comment="#", na_values=[-9, -99])
	sn_data.drop(columns="VARNAMES:", inplace=True)
	sn_data.fillna(no_detection_mag, inplace=True)


# Iterate over model, then filter and make lots of plots
for file in files:
	if NO_ID:
		run_id = file.resolve().parts[-2]

	if file.is_dir():
		file = next(file.glob("*.DUMP"))

	with open(file, "r") as f:
		data = pd.read_csv(f, sep="\s+", comment="#", na_values=[-9, -99, 99])
	data.drop(columns="VARNAMES:", inplace=True)
	data.fillna(no_detection_mag, inplace=True)

	y_vars = [column for column in data.columns if "PEAKMAG" in column]
	rows = int(np.ceil(len(y_vars)/columns))

	if len(data) < 100:
		ms = None
		markerscale = 1
	else:
		ms = 0.5
		markerscale = 30

	fig, axs = plt.subplots(rows, columns, tight_layout=True, figsize=(3*columns,3*rows))
	
	for y_var, ax in zip(y_vars, axs.flat):
		if args.reference is not None:
			ax.plot(sn_data["ZCMB"], sn_data[y_var], ".", ms=ms, label="SNIa", alpha=0.95)
		label = file.stem[17:]
		label = re.sub("^_", "", label)    # labels can't start with _
		ax.plot(data["ZCMB"], data[y_var], ".", ms=ms, label=label, alpha=0.95)
		
		ax.set_xlabel("redshift")
		ax.set_ylim([17, 42])
		ax.set_ylabel("mag")
		ax.invert_yaxis()
		ax.set_title(y_var)
	
	
	axs.flat[columns - 1].legend(loc="center left", bbox_to_anchor=(1, 0.5), markerscale=markerscale, fontsize="large")
	
	for ax in axs.flat:
		if not bool(ax.has_data()):
			fig.delaxes(ax)
	
	fig.suptitle(run_id)
	plt.savefig(
		run_id + "_" + file.stem[17:] + "." + output_format, bbox_inches="tight"
	)
	plt.close()
