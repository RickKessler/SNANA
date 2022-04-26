#!/usr/bin/env python

# Created Apr 20220 by R.Kessler
# Remove garbage stdout produced by MINUIT fitting program.
#

import os, sys

STRING_REMOVE_LIST = [ 'MINUIT', 'MINOS', '=======' ,
                       'EIGENVALUES', 'MINIMIZ', 'VV' ]

# ===================================================
if __name__ == "__main__":
 
    stdout_file = sys.argv[1]
    tmp_file    = f"TMP_{stdout_file}"
    print(f" Remove minuit stdout from {stdout_file}")

    cmd_sed = "sed "
    for string in STRING_REMOVE_LIST:
        cmd_sed += f"'-e /{string}/d' "

    cmd_mv   = f"mv {tmp_file} {stdout_file}"
    cmd_sed += f"{stdout_file} > {tmp_file} ; {cmd_mv}"
    os.system(cmd_sed)

# END
