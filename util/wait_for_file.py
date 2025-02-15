#!/usr/bin/env python

# Created Feb 14 2025
# wait for input file to exist;
# Check for optional string to search inside file.
#
# Examples:
#   wait_for_file.py ALL.DONE
#    -> Script returns 0 when ALL.DONE exists, regardless of contents
#    -> Script holds while waiting for ALL.DONE to exist
#
#   wait_for_file.py ALL.DONE --string SUCCESS
#    -> Script returns 0  when ALL.DONE exists and SUCCESS string is found
#    -> Script returns 6  when ALL.DONE exists and SUCCESS string is not found
#
# ========================================================

import os, sys, datetime, time, argparse, logging

WAIT_TIME_DEFAULT = 10  # default wait time is 10 sec


# ============================
def setup_logging():

    #logging.basicConfig(level=logging.DEBUG, 
    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)6s |%(filename)10s] %(message)s")                        
    # xxx mark format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

# =====================================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "name of file to wait for"
    parser.add_argument("input_file", help=msg, nargs="?", default=None)

    # misc user args

    msg = "optional string to require inside file"
    parser.add_argument("-s", "--string", help=msg, type=str, default=None)

    msg = "wait time between each file check (seconds)"
    parser.add_argument("-w", "--wait_time", help=msg, type=int, default=WAIT_TIME_DEFAULT)

    msg = "verbose mode; print update for each file check"
    parser.add_argument("-v", "--verbose", help=msg, action="store_true")

    args = parser.parse_args()
    return args
    # end get_args

def check_string_exist(args):
    if args.string is None:
        ierr = 0
    else:
        ierr = 6
        with open(input_file_exand,"rt") as f:
            line_list = f.readlines()
        for line in line_list:
            if args.string in line:
                ierr = 0
                        
        if ierr == 0:
            logging.info(f"--> Found string '{args.string}' in {args.input_file}")
        else:
            logging.info(f"Did not find string '{args.string}' in {args.input_file}")


    return ierr

# =============================================
if __name__ == "__main__":
    
    args  = get_args()
    setup_logging()

    input_file_exand = os.path.expandvars(args.input_file)
    t0 = datetime.datetime.now()

    logging.info(f"Begin search every {args.wait_time} sec for {input_file_exand}")
    if not args.verbose:
        logging.info("There will be no update until file is found.")

    while not os.path.exists(input_file_exand):
        if args.verbose:
            t_now   = datetime.datetime.now()
            tstr    = t_now.strftime("%Y-%m-%d %H:%M:%S")
            logging.info(f"Wait for {args.input_file}  ({tstr})")

        time.sleep(args.wait_time)

    # - - - - - 
    # check for optional string inside file

    t_now   = datetime.datetime.now()
    t_search = (t_now - t0).seconds /60.0

    tstr    = t_now.strftime("%Y-%m-%d %H:%M:%S")
    logging.info(f"--> Found {args.input_file} at {tstr}")
    logging.info(f"Search time was {t_search:.2f} minutes")

    ierr = check_string_exist(args)
    logging.info(f"Return ierr={ierr}")


    sys.exit(ierr)
    # ==== END: =====

