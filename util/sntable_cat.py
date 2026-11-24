#!/usr/bin/env python
#
# Created May 2020 by R.Kessler
#
# Read input list of text-formatted SN table files (i.e., FITRES files)
# and catenate them after accounting for possible column differences
# among input files. See manual Sec 10.1 for explanations.
#
# This python script is just a wrapper that calls SALT2mu.exe
# with 'cat_only' argument that turns SALT2mu.exe into a utility 
# instead of a fitting program.
#
# See input arguments with
#    sntable_cat.py -h
#
# Jun 20 2025: minor refactor of inpfile_list to accept either comma-sep list as before,
#              or accept space-sep list or accept wildcard. Replace OptionParser with argparse.
#              Only changes are in parse_args().
# Nov 24 2026: default is now pure python+pandas (no SALT2mu.exe).
#              Inputs may be plain text, text.gz (.gz) or NPZ (.npz/.NPZ).
#              Output format is set by the -o extension (.npz/.NPZ, .gz, or text).
#              Use --legacy to run the old SALT2mu.exe catenation.
# ------------------------------------------

import os, sys, argparse, gzip, glob, fnmatch, json, zipfile
import numpy as np
import pandas as pd
from io import StringIO
#from optparse import OptionParser

DEFAULT_APPEND_VARNAME_MISSING = 'PROB*,zPRIOR*,COV*'

# =============================
def parse_args():

    parser=argparse.ArgumentParser()

    msg = 'list of text-formatted SNANA tables (comma sep or space sep)'
    parser.add_argument('-i', help=msg, nargs="+", type=str,
                        dest='inpfile_list', default=None)

    msg = 'name of output/catenated sntable (text format)'
    parser.add_argument('-o',help=msg, type=str,
                        dest='outfile_cat',default=None)

    msg = f'list of column names to append if missing (wildcard allowed, ' \
          f'default = {DEFAULT_APPEND_VARNAME_MISSING}'
    parser.add_argument('-a',help=msg, type=str,
                        dest='append_varname_missing', default=DEFAULT_APPEND_VARNAME_MISSING)

    msg = 'optional integer prescale'
    parser.add_argument('-p',help=msg, type=int,
                        dest='prescale', default=1)

    # - - - - - 
    msg = 'legacy option using SALT2mu.exe (C code)'
    parser.add_argument('--legacy', help=msg, action="store_true")

    msg = 'LEGACY ONLY: SALT2mu jobname'
    parser.add_argument('-j', help=msg, type=str,
                        dest='jobname', default='SALT2mu.exe')

    INPUTS  = parser.parse_args()

# abort on missing arguments
    if not INPUTS.inpfile_list:
        parser.error('Missing required arg: -i <inpfile_list>')

    # 6/20/2025: convert list to comma-sep string for C code input
    if len(INPUTS.inpfile_list) == 1:
        INPUTS.inpfile_list = INPUTS.inpfile_list[0] # convert list to string
    else:
        INPUTS.inpfile_list = ','.join(INPUTS.inpfile_list) 

    if not INPUTS.outfile_cat:  
        parser.error('Missing required arg: -o <outfile_cat>')

    return INPUTS


# ===================================
def exec_salt2mu(INPUTS):

    inpfile_list   = INPUTS.inpfile_list
    outfile_cat    = INPUTS.outfile_cat
    append_varname = INPUTS.append_varname_missing
    arg_cat        = "cat_only"
    jobname        = INPUTS.jobname

    if INPUTS.prescale > 1: 
        arg_cat += f"/{INPUTS.prescale}"

    arg_list = f"{arg_cat} " \
               f"datafile={inpfile_list} "   \
               f"catfile_out={outfile_cat} " \
               f"append_varname_missing={append_varname} "

#                % (arg_cat, inpfile_list, outfile_cat, append_varname) )

    command = f"{jobname} {arg_list}"
    print(f" {sys.argv[0]}: executes os.system command = \n  '{command}' \n")
    sys.stdout.flush()

    os.system(command)

    print(f"\n Done.")

def is_npz(fname):
    return fname.lower().endswith('.npz')

def is_gz(fname):
    return fname.lower().endswith('.gz')

# keys written by this program (re-created on output, so never inherited);
# the legacy names are the keys written by SALT2mu.exe cat_only.
KEYS_NEW = ['PROGRAM_TABLE_CAT', 'CATENATED_TABLE_FILES', 'APPENDED_COLUMNS', 'DROPPED_COLUMNS', 'PRESCALE']
KEYS_SKIP = KEYS_NEW + ['CATENATED DATA FILES', 'APPENDED COLUMNS',
                        'DROPPED COLUMNS', 'DOCUMENTATION']
KEY_VERSION_PHOTOMETRY = 'VERSION_PHOTOMETRY'

# ===================================
def parse_comments(comment_lines):
    """'# KEY: value' lines -> dict. Continuation lines such as '#   - file'
    or '#   + file' and free-text comments without a key are ignored."""
    header = {}
    for line in comment_lines:
        body = line.lstrip('#').strip()
        if line.startswith('#  ') and body[:1] in '+-':  # list item
            continue
        if ':' not in body:
            continue
        key, val = body.split(':', 1)
        key = key.strip()
        if not key or key.upper() in KEYS_SKIP or ' ' in key:
            continue
        header[key] = val.strip()
    return header

# ===================================
def read_table(fname):
    """Return (DataFrame, header dict) for text, text.gz, or npz.
    Header dict = {KEY: value} taken from '# KEY: value' comment lines (text)
    or from the zip-archive comment (JSON string:string dict) for npz."""
    if is_npz(fname):
        z = np.load(fname, allow_pickle=True)
        n = None
        cols = {}
        for key in z.files:
            arr = z[key]
            if arr.ndim != 1:
                continue
            if n is None:
                n = len(arr)
            if len(arr) != n:
                continue
            if arr.dtype.kind == 'S':
                arr = np.char.strip(np.char.decode(arr, 'utf-8'))
            cols[key] = arr
        header = {}
        try:
            raw = json.loads(zipfile.ZipFile(fname).comment.decode('utf-8') or '{}')
            header = {k.strip().rstrip(':'): str(v).strip() for k, v in raw.items()
                      if k.strip().rstrip(':').upper() not in KEYS_SKIP}
        except Exception:
            pass
        return pd.DataFrame(cols), header

    # text or gzipped text
    opener = gzip.open if is_gz(fname) else open
    with opener(fname, 'rt') as f:
        comments = [l.rstrip('\n') for l in f if l.startswith('#')]
    header = parse_comments(comments)
    compression = 'gzip' if is_gz(fname) else None
    df = pd.read_csv(fname, compression=compression, sep=r'\s+',
                     comment='#', dtype={'CID': str})
    # first column is the 'VARNAMES:' key (rows start with 'SN:'): remove it
    key = df.columns[0]
    if key.upper().startswith('VARNAMES'):
        df = df[df[key].astype(str).str.upper() == 'SN:'].drop(columns=key)
        df = df.reset_index(drop=True)
    return df, header

# ===================================
def write_table(df, fname, header, new_keys):
    """header: {KEY: str} inherited keys. new_keys: {KEY: list or str}
    written after the inherited keys.
    Text: '# KEY: value', lists as one '#   - item' line per item.
    NPZ : flat JSON {"KEY": "value"} in the zip-archive comment (same as
          snana.exe NPZ tables, so C code and python can both read it);
          lists are written as the string '[ item1, item2, ... ]'."""
    if is_npz(fname):
        arrays = {}
        for c in df.columns:
            if df[c].dtype == object:
                arrays[c] = df[c].astype(str).to_numpy().astype('S')
            else:
                arrays[c] = df[c].to_numpy()
        # use file object so numpy does not append '.npz' to a '.NPZ' name
        with open(fname, 'wb') as f:
            np.savez_compressed(f, **arrays)

        hdr = dict(header)
        for k, v in new_keys.items():
            hdr[k] = '[ ' + ', '.join(v) + ' ]' if isinstance(v, list) else str(v)
        comment = json.dumps(hdr)
        if len(comment) > 65535:   # zip comment size limit
            hdr["CATENATED_TABLE_FILES"] = f"[ {len(new_keys['CATENATED_TABLE_FILES'])} files ]"
            comment = json.dumps(hdr)
        with zipfile.ZipFile(fname, 'a') as zf:
            zf.comment = comment.encode('utf-8')
        return

    hdr_lines = [f"# {k}: {v}" for k, v in header.items()]
    for k, v in new_keys.items():
        if isinstance(v, list):
            if v:
                hdr_lines += [f"# {k}:"] + [f"#   - {item}" for item in v]
        else:
            hdr_lines.append(f"# {k}: {v}")
    body = df.to_csv(sep=' ', header=False, index=False, float_format='%.8g')
    body = ''.join('SN: ' + l + '\n' for l in body.splitlines())
    text = '\n'.join(hdr_lines) + '\n#\n' \
         + 'VARNAMES: ' + ' '.join(df.columns) + '\n\n' + body
    opener = gzip.open if is_gz(fname) else open
    with opener(fname, 'wt') as f:
        f.write(text)

# ===================================
def exec_table_cat(INPUTS):
    """Catenate tables (text, .gz, .npz) using only python+pandas.
    Columns in every file are kept. Columns matching the -a patterns that are
    missing in some files are appended with value -9. Other columns that are
    not in all files are dropped. Prescale (-p) is applied after catenation."""

    # expand comma-sep list and wildcards
    flist = []
    for item in INPUTS.inpfile_list.split(','):
        item = item.strip()
        if not item:
            continue
        flist += sorted(glob.glob(item)) if any(c in item for c in '*?[') else [item]
    if not flist:
        sys.exit(f"\n ERROR: no input files found in '{INPUTS.inpfile_list}'")

    patterns = [a.strip() for a in INPUTS.append_varname_missing.split(',') if a.strip()]
    def match_append(col):
        return any(fnmatch.fnmatch(col, pat) for pat in patterns)

    print(f" {sys.argv[0]}: catenate {len(flist)} tables (pure python)")
    dfs = []
    headers = []
    for f in flist:
        if not os.path.exists(f):
            sys.exit(f"\n ERROR: cannot find input file {f}")
        df, header = read_table(f)
        print(f"   read {len(df):6d} rows x {len(df.columns):3d} cols from {f}")
        dfs.append(df)
        headers.append(header)

    # column bookkeeping (keep order of first appearance)
    all_cols = list(dict.fromkeys(c for df in dfs for c in df.columns))
    common   = [c for c in all_cols if all(c in df.columns for df in dfs)]
    partial  = [c for c in all_cols if c not in common]
    append   = [c for c in partial if match_append(c)]
    dropped  = [c for c in partial if c not in append]

    for i, df in enumerate(dfs):
        for c in append:
            if c not in df.columns:
                df[c] = -9
        dfs[i] = df[common + append]

    df_cat = pd.concat(dfs, ignore_index=True)

    if INPUTS.prescale > 1:
        df_cat = df_cat.iloc[::INPUTS.prescale].reset_index(drop=True)

    # header: inherited keys (first file wins), except VERSION_PHOTOMETRY
    # which becomes a comma-separated list over all input files
    header = {}
    for h in headers:
        for k, v in h.items():
            header.setdefault(k, v)
    versions = []
    for h in headers:
        versions += [x.strip() for x in h.get(KEY_VERSION_PHOTOMETRY, '').split(',')]
    versions = list(dict.fromkeys(x for x in versions if x))
    if versions:
        header[KEY_VERSION_PHOTOMETRY] = ','.join(versions)

    new_keys = {'PROGRAM_TABLE_CAT': sys.argv[0],
                'CATENATED_TABLE_FILES': flist,
                'APPENDED_COLUMNS': patterns,
                'DROPPED_COLUMNS': dropped}
    if INPUTS.prescale > 1:
        new_keys['PRESCALE'] = INPUTS.prescale

    write_table(df_cat, INPUTS.outfile_cat, header, new_keys)

    print(f"   appended (with -9 where missing): {append}")
    print(f"   dropped (not in all files)      : {dropped}")
    print(f"\n Wrote {len(df_cat)} rows x {len(df_cat.columns)} cols to {INPUTS.outfile_cat}")
    print(f"\n Done.")

# ===================================
# ============== MAIN ===============
# ===================================

if __name__ == "__main__":

    INPUTS = parse_args()

    if INPUTS.legacy:
        # execute SALT2mu.exe program to do the catenate
        exec_salt2mu(INPUTS)
    else:
        # Nov 24 2026: switch it all over to python, and work with gz or npz
        exec_table_cat(INPUTS)

