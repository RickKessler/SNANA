#!/usr/bin/env python

# Created Oct 2022 by R.Kessler
#  [translate remove_locf_messages.pl to python]
#
# Remove annoying locf messages from [log] file.
# These messages come from 64-bit CERNLIB calls .
#
# Usage:
#   remove_locf_messages.py  <file>
#      or
#   remove_locf_messages.py  <file>  QUIET
#
# Note that <file> is modified, so if you want to
# keep the original then copy it first .
#
# The sed command is
#  sed -e '/locf/d' -e '/crash/d' -e '/!!!!!!/d' $INFILE 
#
# Jul 15 2023: replace os.system with subprocess.run

import os, sys, argparse, subprocess, datetime

# =================================================== 
def get_args():
    parser = argparse.ArgumentParser()

    msg = "name of file to remove locf messages"
    parser.add_argument("input_file", help=msg,
                        nargs="?", default=None)

    msg = "Quiet mode"
    parser.add_argument("--quiet", "-q", help=msg,
                        action="store_true")
     
    args = parser.parse_args()
    return args
    # end get_args


# =======================================
if __name__ == "__main__":

    args = get_args()
    input_file = args.input_file

    t0 = datetime.datetime.now()

    tmp_file = f"{input_file}.TEMP" ;

    cmd = f"sed " \
          f"-e '/locf/d' " \
          f"-e '/crash/d' " \
          f"-e '/!!!!!!/d' " \
          f"{input_file} > " \
          f"{tmp_file} ; " \
          f"mv {tmp_file}  {input_file} "

    # xxx mark delete 7/15/2023 os.system(cmd)
    
    ret = subprocess.run( [ cmd ], shell=True,
                          capture_output=False, text=True )
    if not args.quiet :
        t1 = datetime.datetime.now()
        dt = (t1-t0).total_seconds()
        print(f"  locf messages removed from {input_file} ({dt} seconds)")

# END
