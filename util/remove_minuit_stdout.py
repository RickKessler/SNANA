#!/usr/bin/env python

# Created Apr 20220 by R.Kessler
# Remove garbage stdout produced by MINUIT fitting program.
#
# Jul 15 2023: replace os.system with subprocess.run to hopefully
#       get better performance on batch core.
#
# Dec 5 2023: re-write script to brute-force write stdout file
#      instead of using subprocess ... hopefully less hangup.
#
# May 30 3024: remove MINOs instead of MINOS so that snana's MINOS comments remain.
#       Previous submit_batch_jobs output with USE_MINOS=T has no record of MINOS
#       in the log files. Most MINUIT output has been suppressed anyway in minuit.F
#

import os, sys, subprocess, datetime

STRING_REMOVE_LIST = [ 'MINUIT', 'MINOs', '=======' ,
                       'EIGENVALUES', 'MINIMIZ', 'VV' ]

# ===================================================
if __name__ == "__main__":

    t0 = datetime.datetime.now()

    stdout_file = sys.argv[1]

    with open(stdout_file,"rt") as f:
        lineList = f.readlines() 

    # re-write stdout file excluding some lines
    f =  open(stdout_file,"wt") 
    for line in lineList:
        keep = True
        for str in STRING_REMOVE_LIST:
            if str in line: 
                keep = False
        if keep:
            f.write(f"{line}")

    f.close()

    t1 = datetime.datetime.now()
    dt = (t1-t0).total_seconds()
    print(f" Remove minuit stdout from {stdout_file}  ({dt} seconds)")

# END
