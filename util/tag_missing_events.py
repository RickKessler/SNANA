#!/usr/bin/env python
#
# Created May 2025 by R.Kessler
#
# Read in list of tables and count number of tables missing each each event;
# Append resulting N_MISSING colunm to original reference table.
# Original use is submit_batch_jobs with BBC task to flag the following
# for each event:
#   NFITOPT_REJECT_LCFIT
#   NFITOPT_REJECT_BIASCOR
#
# and enable monitoring loss from BBC's CUTWIN separately from the biasCor loss.
#
#

import os, sys, logging, argparse, glob
import pandas as pd
import numpy  as np

COLNAME_UCID        = 'ucid'  # temporray col name for unique cid__idsurvey__field
DEFAULT_COLNAME_ADD = "N_MISSING"
DEFAULT_OUTFILE     = "tag_missing_events.fitres"

# ============================
def setup_logging():

    logging.basicConfig(level=logging.INFO,
        format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")

def get_args():

    parser = argparse.ArgumentParser()

    msg = "name of reference table"
    parser.add_argument("--tfile_ref", help=msg, nargs="?", default=None)

    
    msg = "list of tables to check for missing events (can use wildcard)"
    parser.add_argument("--tfile_list", help=msg, nargs="+", default=None)

    msg = "name of column to add to ref table (default={DEFAULT_COLNAME_ADD})"
    parser.add_argument("--ref_colname_add", help=msg, nargs="?", default=DEFAULT_COLNAME_ADD)

    msg = "name of output table file with added column"
    parser.add_argument("--outfile", help=msg, nargs="?", default=DEFAULT_OUTFILE )

    msg = "List of comment string to write at top of output file"
    parser.add_argument("--comments", help=msg, nargs="+", default=None )


    # parse it
    args = parser.parse_args()

    return args
    # end get_args

def read_input_tables(args):

    df_list   = []
    ucid_list = []
    tfile_local_list = [ args.tfile_ref ] + args.tfile_list

    ntf = len(args.tfile_list)
    logging.info(f"Read ref table + {ntf} tables to check for missing events.")

    for tf in tfile_local_list: 
        if '.gz' in tf:
            df = pd.read_csv(tf, compression='gzip',sep='\s+',comment="#")
        else:
            df = pd.read_csv(tf, comment='#', delim_whitespace=True)

        # define unique cid (ucid) base on CID, IDSURVYE, FIELD
        df[COLNAME_UCID] = df.CID.astype(str) + "__" + df.IDSURVEY.astype(str) + "__" + \
                           df.FIELD.astype(str)
        
        # keep subset of rows in df_ref
        if len(df_list) == 0 :
            df_ref = df
        else:
            df = df[df[COLNAME_UCID].isin(df_ref[COLNAME_UCID])]

        df_list.append(df)    

    return df_list

def write_output(args, df):

    logging.info(f"Write {args.outfile} ")

    with  open(args.outfile, "wt") as f:
        if args.comments is not None:
            for comment in args.comments:
                f.write(f"# {comment} \n")
    
    df_ref.to_csv(args.outfile, mode='a', sep=' ', index=False)

    return

# ===================================================
if __name__ == "__main__":

    setup_logging()
    logging.info("# ========== BEGIN tag_missing_events.py ===============")
    logging.info("# full command: {sys.argv} ")

    args   = get_args()

    # store all tables in list of data frames
    df_list = read_input_tables(args) 
    n_tf    = len(df_list)
    nevt    = len(df_list[0])

    # combine all UCID lists together
    ucid_list   = []
    for df in df_list:
        df_id    = df[COLNAME_UCID]
        ucid_list = np.concatenate( (ucid_list, df_id) )  # list over all files

    # count number of times each UCID appears
    ucid_unique, ucid_index, n_count = np.unique(ucid_list, return_index=True, return_counts=True)
    n_reject_sorted                  = n_tf - n_count
    n_reject = np.array( [0] * nevt )
    for j, ind in enumerate(ucid_index):
        n_reject[ind] = n_reject_sorted[j]  # aligns with orginal data frame

    # - - - -

    # strip off reference data frame
    df_ref = df_list[0]

    # add n_reject colunm 
    logging.info(f"Append {args.ref_colname_add} column")
    df_ref[args.ref_colname_add] = n_reject

    # remove UCID column before writing it to csv file
    df_ref = df_ref.drop(COLNAME_UCID, axis=1)

    write_output(args,df_ref)

    # - - - - - -
    debug_dump = False
    if debug_dump :
        len_count = len(n_count)
        len_rej   = len(n_reject)
        len_ucid  = len(ucid_unique) 
        print(f"\n xxx n_count({len_count})  = \n{n_count}")
        print(f"\n xxx n_reject({len_rej}) = \n{n_reject}")
        print(f"\n xxx len_ucid={len_ucid}")
        print(f"\n xxx n_tf = {n_tf}")

        ucid_list_orig = df_list[0][COLNAME_UCID].to_numpy()
        print(f"\n xxx ucid_list_orig = \n{ucid_list_orig[:20]}")
        print(f"\n xxx ucid_list {len(ucid_list)}      = \n{ucid_list[:20]}")
        print(f"\n xxx ucid_uniq = \n{ucid_unique[:20]}")
        print(f"\n xxx ucid_index ({len(ucid_index)}) = \n{ucid_index[:20]}" )
        print(f"\n xxx n_rej(sorted/unsorted) = \n{n_reject_sorted}\n\n{n_reject}")
        
    # === END:
