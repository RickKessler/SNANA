#!/usr/bin/env python
#
# Created August 2022 by R.Kessler
# Utility to modif SED time-series model (either SIMSED or NON1ASED).
#
# CUT-Modification:
#   Cut on DAY and/or LAM range.
#   First usage is to remove DAY<0 for KN-K17 model where 
#   SED fluxes are all zero.
#
#   Example usage:
#    simsed_trim.py \
#     --path_model $SNDATA_ROOT/models/SIMSED/SIMSED.KNovae_BarnesKasen2013 \
#     --cutwin_day 4 50 \
#     --cutwin_lam 1250 25000
#
# will output SIMSED.KNovae_BarnesKasen2013_trim subdirectory
# and each SED time series will be truncated with DAY>4 and LAM>1250 A.
#
# EXTRAP-Modification
#  --early_extrap <Fearly/Fpeak> extrapolates rise time until
#  Fearly/Fpeak reaches desired level.
#
#
# ========================

import os, sys, argparse, glob, yaml, shutil, gzip, math

MODEL_CLASS_SIMSED   = "SIMSED"
MODEL_CLASS_NON1ASED = "NON1ASED"

FLAG_MODIFY = +1
FLAG_COPY   =  0
FLAG_REMOVE = -1

# define info file key that has sed info
SEDKEY_DICT = { 
    MODEL_CLASS_SIMSED:   [ "SED:" ],
    MODEL_CLASS_NON1ASED: [ "NON1A:", "NONIA:" ] 
}

IWD_DICT = {
    MODEL_CLASS_SIMSED:   1,  # sed file name is 1st word after key
    MODEL_CLASS_NON1ASED: 3   # sed file name is 3rd word after key
}

INFO_FILE_DICT = {
    MODEL_CLASS_SIMSED:   "SED.INFO",
    MODEL_CLASS_NON1ASED: "NON1A.LIST"
}

# =========================
def get_args():
    parser = argparse.ArgumentParser()

    msg = "path of SIMSED model "
    parser.add_argument("--path_model", help=msg, type=str, default=None)

    msg = "DAY cut-window to keep"
    parser.add_argument("--cutwin_day", help=msg, nargs=2, type=float)

    msg = "LAM cut-window to keep"
    parser.add_argument("--cutwin_lam", help=msg, nargs=2, type=float)

    # - - -
    msg = "early extrap to get Flux(early)/Flux(peak) < this value"
    parser.add_argument("--early_extrap", help=msg, type=float, default=None)

    msg = "max number of seds to process"
    parser.add_argument("--mxsed", help=msg, type=int, default=9999)

    args = parser.parse_args()

    args.path_model = os.path.expandvars(args.path_model)

    args.trim = False
    if args.cutwin_day or args.cutwin_lam:
        args.trim = True

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return args

    # end get_args

def read_sed_list(model):

    # read SED.INFO file and return list of SED files

    if MODEL_CLASS_SIMSED in model:
        model_class = MODEL_CLASS_SIMSED
    elif MODEL_CLASS_NON1ASED in model:
        model_class = MODEL_CLASS_NON1ASED
    else:
        sys.exit(f"\n ERROR: Cannot find info file for model = {model}")
    
    info_base = INFO_FILE_DICT[model_class]
    info_file = f"{model}/{info_base}"
    iwd_sed   = IWD_DICT[model_class]
    KEY_LIST  = SEDKEY_DICT[model_class]
    sed_list = []

    #print(f" xxx model_class={model_class}   info_file={info_file}")

    with  open(info_file,"rt") as f:
        for line in f:
            line  = line.rstrip() 
            if len(line) == 0 : continue
            wdlist   = line.split()
            key      = wdlist[0]
            if key in KEY_LIST:
                sed_file = os.path.basename(wdlist[iwd_sed])
                sed_list.append(sed_file)
         
    return model_class, sed_list

    # end read_sed_list

# ==========================================
def create_outdir(model_class, path_model):

    basename = os.path.basename(path_model)
    outdir   = f"{basename}_modified"

    if os.path.exists(outdir):
        shutil.rmtree(outdir)

    print(f" Create {outdir}")
    os.mkdir(outdir)

    # copy original INFO file and some others
    info_base = INFO_FILE_DICT[model_class]
    if model_class == MODEL_CLASS_SIMSED:
        copy_list = [ info_base ]
    else:
        copy_list = [ info_base, "SIMGEN*", "SUGGEST*" ]

    for item in copy_list:
        copy_file = f"{path_model}/{item}"
        print(f"  Copy {item} to output dir")
        if item == info_base:
            # re-write info file with path removed from sed file name
            write_info_file(model_class,copy_file,outdir)
        else:
            # simple linux copy
            cmd       = f"cp {copy_file} {outdir}"
            os.system(cmd)

    return outdir
    # end create_outdir

def write_info_file(model_class, info_file_orig, outdir):

    info_base = INFO_FILE_DICT[model_class]
    iwd_sed   = IWD_DICT[model_class]
    KEY_LIST  = SEDKEY_DICT[model_class]

    line_list_out = []

    with open(info_file_orig,"rt") as f:
        for line in f:
            line  = line.rstrip() 
            line_out = line
            if len(line) > 0 : 
                wdlist   = line.split()
                key      = wdlist[0]
                if key in KEY_LIST:
                    sed_file_orig = wdlist[iwd_sed]
                    sed_file_base = os.path.basename(sed_file_orig)
                    line_out = line_out.replace(sed_file_orig,sed_file_base)
            line_list_out.append( line_out )

    # - - - -
    info_file_out = f"{outdir}/{info_base}"
    with open(info_file_out,"wt") as f:
        for line in line_list_out:
            f.write(f"{line}\n")
         
    return
    # end write_info_file

def open_sed_outfile(outdir,sed):
    sed_file_out = f"{outdir}/{sed}.gz"
    sout = gzip.open(sed_file_out, "wt", compresslevel=6)
    return sout

def open_log_file(outdir,args):
    
    log_file = f"{outdir}/simsed_modify.log"
    f_log    = open(log_file,"wt")

    print(f" Open {log_file} to record modifications")
    
    if args.early_extrap:
        f_log.write(f"Extrapolate early time fluxes until " \
                    f"bolometric FluxEarly/FluxMax < " \
                    f"{args.early_extrap} \n\n")

    if args.trim:
        pass  # write something ??

    f_log.flush()

    return f_log

    #end open_log_file

def read_sed_file(args,sed):

    # read sed file and return contents (list of strings)

    path_model = args.path_model

    sed_file_orig   = f"{path_model}/{sed}"
    
    if not os.path.exists(sed_file_orig):
        sed_file_orig += '.gz'

    sin  = gzip.open(sed_file_orig,"rt")
    sed_contents = sin.readlines()

    sin.close()

    return sed_contents
    # end read_sed_file

def early_extrap(sed_file, sed_contents, args, sout, f_log):

    # Extrapolate earch sed and return
    #  +1 --> modified sed
    #   0 --> original sed is ok, so copy without modification
    #  -1 --> unable to extrapolate sed, hence remove it
 
    text = f"\n   Check early-epoch extrapolation for  {sed_file}"
    print_stdout_and_log(f_log,text)

    line_dict = {} # each day contains list of lines
    line_comment_list = []
    day_list          = []  # list of SED days
    lam_list          = []  # list of SED wavelengths
    flux_lam0_list    = []  # flux vs. lam in day bin 0
    flux_lam1_list    = []  # flux vs. lam in day bin 1
    flux_lam2_list    = []  # flux vs. lam in day bin 2
    fluxInteg         = {}
    nbin_day  = 0


    #  - - - - -
    for line_orig in sed_contents:
        line = line_orig.rstrip()
        if line[0] == '#' : 
            line_comment_list.append(line)
            continue

        day, lam, flux = parse_sed_line(line)

        if day not in day_list:
            day_list.append(day)
            line_dict[day] = []
            fluxInteg[day] = 0.0
            nbin_day += 1

        line_dict[day].append(line)
        fluxInteg[day] += flux 

        if nbin_day == 1:
            flux_lam0_list.append(flux)
            lam_list.append(lam)
        if nbin_day == 2:
            flux_lam1_list.append(flux)
        if nbin_day == 3:
            flux_lam2_list.append(flux)

    # - - - - - - - - - 
    # check bolometric flux-ratio at Trestmin compared to peak
    daymin_orig = day_list[0]
    day_at_peak = max(fluxInteg, key=fluxInteg.get)    
    fratio      = fluxInteg[daymin_orig] / fluxInteg[day_at_peak]
    text = f"\t FluxEarly[{daymin_orig}d] / FluxMax[{day_at_peak}d] " \
           f"= {fratio:.3f}"
    print_stdout_and_log(f_log,text)

    if fratio < args.early_extrap :
        write_sed(sout, line_comment_list, day_list, line_dict)
        return FLAG_COPY

    # prepend early extrapolation

    mag0_list  = [ -2.5*math.log10(x) for x in flux_lam0_list]
    mag1_list  = [ -2.5*math.log10(x) for x in flux_lam1_list]
    mag2_list  = [ -2.5*math.log10(x) for x in flux_lam2_list]

    # if any mag0-mag1 has wrong slope, use slope in next bin
    magdif_list = [ x - y for x,y in zip(mag0_list,mag1_list) ]
    if any(t < 0 for t in magdif_list): 
        magdif_list = [ x - y for x,y in zip(mag1_list,mag2_list) ]
        if any(t < 0 for t in magdif_list): 
            text = f" wrong mag-slope in first two phase bins --> " \
                   f" cannot extrapolate rise time for {sed_file}"
            print_stdout_and_log(f_log,text)
            return FLAG_REMOVE

        mag0_list = mag1_list.copy()
        mag1_list = mag2_list.copy()
        day_list.remove(day_list[0])
        text = f"\t wrong mag-slope in first phase bin of {sed_file}; " \
               f"use next bin"
        print_stdout_and_log(f_log,text)

    day0      = day_list[0]
    day_bin   = day_list[1] - day_list[0]
    day_early = day0
    fac       = 1.0

    while fratio > args.early_extrap :
        day_early -= day_bin
        day_ratio = (day0-day_early)/day_bin

        day_ratio *= fac
        fac *= 1.1

        # extrapolate with list vs. wavelength
        mag_shift_list = [ (x - y)*day_ratio for x, y in \
                           zip(mag0_list,mag1_list)]  
        mag_early_list = [ x+y for x,y in zip(mag0_list,mag_shift_list) ]
        flux_extrap_list = [ pow(10.0,-0.4*x) for x in mag_early_list ]

        fluxInteg[day_early] = sum(flux_extrap_list) 
        fratio = fluxInteg[day_early] / fluxInteg[day_at_peak]

        text = f"\t prepend extrap flux at day = {day_early:7.3f} " \
              f" F/Fpeak={fratio:.3f}"
        print_stdout_and_log(f_log,text)

        day_list = [day_early] + day_list
        line_dict[day_early] = []

        # loop over lam to add lines to output sed file
        for lam, flux_extrap in zip(lam_list,flux_extrap_list):
            line_extrap = f"   {day_early} {lam} {flux_extrap:.4e} "
            line_dict[day_early].append(f"{line_extrap}")

    # - - - - - - 
    # write output lines ordered by day
    write_sed(sout, line_comment_list, day_list, line_dict)

    return FLAG_MODIFY
    # end early_extrap

def write_sed(sout, line_comment_list, day_list, line_dict):

    for line in line_comment_list:
        sout.write(f"{line}\n")

    for day in day_list:
        for line in line_dict[day]:
            sout.write(f"{line}\n")
    return
    # end write_sed

def early_extrap_obsolete(sed_name, sed_contents, args, sout):

    fluxInteg = {}
    line_comment_list = []
    daymin_orig = 99999.0
    day_list = []
    day_last = 999999.
    nbin_day = 0
    flux_lam0 = []
    flux_lam1 = []

    line_dict = {} # each day contains list of lines

    for line_orig in sed_contents:
        line = line_orig.rstrip()
        if line[0] == '#' : 
            line_comment_list.append(line)
            continue

        day, lam, flux = parse_sed_line(line)

        if day not in day_list:
            day_list.append(day)
            line_dict[day] = []
            nbin_day += 1

        line_dict[day].append(line)

        if nbin_day == 1:
            flux_lam0.append(flux)
        if nbin_day == 2:
            flux_lam1.append(flux)

        # integrate flux over wavelength in each day bin
        if day in fluxInteg:
            fluxInteg[day] += flux            
        else:
            fluxInteg[day] = flux

    # - - - - - - - -
    # find day at peak flux
    daymin_orig = day_list[0]
    flux_max    = 0.0
    day_at_peak = max(fluxInteg, key=fluxInteg.get)    
    fratio      = fluxInteg[daymin_orig] / fluxInteg[day_at_peak]

    day_dif    = day_list[1]-day_list[0]
    day_new    = day_list[0] - day_dif
    flux_slope = (flux_lam1 - flux_lam0) / day_dif
    flux_new   = flux_lam0  - flux_slope * day_dif

    # estimate mag/day
    day0   =  day_list[0]
    day1   =  day_list[1]
    day2   =  day_list[2]
    f0     =  fluxInteg[day0]
    f1     =  fluxInteg[day1]
    f2     =  fluxInteg[day2]
    magdif1 = -2.5*math.log10(f1/f0)  
    magdif2 = -2.5*math.log10(f2/f1)  
    daydif1 = day1 - day0
    daydif2 = day2 - day1
    slope1  = magdif1 / daydif1
    slope2  = magdif2 / daydif2
    do_dump = True
    if do_dump:
        print(f" {sed_name}: day_at_peak = " \
              f"{day_at_peak:5.1f}  Frat={fratio:6.3f}  " \
              f"rise= {slope1:.2f},{slope2:.2f} mag/day")

    # NOT FINISHED: not clear how to extrapolate in UV when there is no rise 
    
    return
    # end early_extrap_obsolete

def print_stdout_and_log(f_log,text):
    print(f"{text}")
    f_log.write(f"{text}\n")
    f_log.flush()

    return

# =======================================================
def trim_sed(sed_file, sed_contents, args, sout, f_log):

    print_stdout_and_log(f_log, f"   Trim  {sed_file}")
    flag = FLAG_COPY

    for line_orig in sed_contents:            
        line_new = line_orig.rstrip()
        if line_new[0] == '#' :
            keep = True
        else:
            day, lam, flux = parse_sed_line(line_new)
            keep = apply_cutwin(args,day,lam)            
        if keep: 
            sout.write(f"{line_new}\n")
        else:
            flag = FLAG_MODIFY

    return flag

    # end trim_sed

def parse_sed_line(line):
    wdlist   = line.split()
    day      = float(wdlist[0])
    lam      = float(wdlist[1])
    flux     = float(wdlist[2])
    if flux < 0: flux = -flux

    return day, lam, flux

def apply_cutwin(args,day,lam):

    keep = True

    if args.cutwin_day is not None:
        if day < args.cutwin_day[0] : keep = False
        if day > args.cutwin_day[1] : keep = False

    if args.cutwin_lam is not None:
        if lam < args.cutwin_lam[0] : keep = False
        if lam > args.cutwin_lam[1] : keep = False

    return keep
    # end apply_cutwin

def remove_sed(info_file,sed_file):

    # remove line containing {sed} from info_file.
    # This function modifies an info file; no file is actually removed
    # Note this function uses linux sed utility (streamed editori),
    # which has the same acronym as SED for spectral energy distribution.

    cmd_remove = "sed -i '/{sed_file}/d' {info_file}" 
    os.system(cmd_remove)

    return

# =====================================
#
#      MAIN
#
# =====================================

if __name__ == "__main__":

    args     = get_args()
    model_class, sed_list = read_sed_list(args.path_model)
    outdir   = create_outdir(model_class, args.path_model)

    f_log = open_log_file(outdir,args)

    info_file_out = f"{outdir}/{INFO_FILE_DICT[model_class]}"

    n_sed = len(sed_list)
    n_modify = 0
    n_remove = 0
    n_copy   = 0 
    n_sed_proc = 0

    # sed_list = [ 'pycoco_SN2008D.SED' ]  # debug

    for sed_file in sed_list:  # sed is base name
        sout         = open_sed_outfile(outdir,sed_file)
        sed_contents = read_sed_file(args,sed_file)


        if args.trim:
            flag = trim_sed(sed_file, sed_contents, args, sout, f_log)
        if args.early_extrap :
            flag = early_extrap(sed_file, sed_contents, args, sout, f_log)

        sout.close()

        # - - - update log file and stats based on flag - - - -
        if flag == FLAG_REMOVE:
            n_remove += 1
            remove_sed(info_file_out,sed_file)
            print_stdout_and_log(f_log,f"\t ACTION: Removed {sed_file} ")
        elif flag == FLAG_MODIFY:
            n_modify += 1
            print_stdout_and_log(f_log,f"\t ACTION: Modified {sed_file}")
        elif flag == FLAG_COPY:
            n_copy += 1 
            print_stdout_and_log(f_log,f"\t ACTION: Copied {sed_file} (no change)")

        else:
            sys.exit(f"\n ERROR: Invalid flag = {flag}")

        n_sed_proc += 1
        if n_sed_proc >= args.mxsed : break

    # - - - - -

    text = f"\n SUMMARY:"
    print_stdout_and_log(f_log,text)

    n_list   = [ n_modify, n_copy, n_remove ]
    str_list = [ 'Modified' , 'Copied', 'Removed' ]

    for n, str in zip(n_list, str_list):
        text = f"   {str:<12} {n} of {n_sed} SEDs"
        print_stdout_and_log(f_log,text)

    f_log.close()

    # end main
