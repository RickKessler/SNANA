#!/usr/bin/env python

# Created Apr 20220 by R.Kessler
# Remove garbage stdout produced by MINUIT fitting program.
#
# Jul 15 2023: replace os.system with subprocess.run to hopefully
#       get better performance on batch core.
#

import os, sys, subprocess, datetime

STRING_REMOVE_LIST = [ 'MINUIT', 'MINOS', '=======' ,
                       'EIGENVALUES', 'MINIMIZ', 'VV' ]

# ===================================================
if __name__ == "__main__":

    t0 = datetime.datetime.now()

    stdout_file = sys.argv[1]
    tmp_file    = f"TMP_{stdout_file}"

    cmd_sed = "sed "
    for string in STRING_REMOVE_LIST:
        cmd_sed += f"'-e /{string}/d' "

    cmd_mv   = f"mv {tmp_file} {stdout_file}"
    cmd_sed += f"{stdout_file} > {tmp_file} ; {cmd_mv}"

    # xxx mark delete 7/15/2023 os.system(cmd_sed)

    ret = subprocess.run( [ cmd_sed ], shell=True,
                          capture_output=False, text=True )

    t1 = datetime.datetime.now()
    dt = (t1-t0).total_seconds()
    print(f" Remove minuit stdout from {stdout_file}  ({dt} seconds)")

# END
