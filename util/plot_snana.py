#!/usr/bin/env python
#
#  June 2019 J. Pierel
#  Aug 2020 B. Sanchez
#  Plotter tool for SNANA LCs and Spectra

from __future__ import print_function

import matplotlib as mpl
import numpy as np

import glob
import math
import os
import sys
import textwrap
from optparse import SUPPRESS_HELP, OptionParser

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d

mpl.use('Agg')


BANDS = ["u", "b", "g", "v", "r", "i", "z", "y", "j", "h", "k"]
__band_order__ = np.append(BANDS, [x.upper() for x in BANDS])


def read_spec(cid, base_name):
    """Function that reads the spectra of a single object

    Parameters
    ----------
    cid:

    base_name:

    Returns
    -------
    sn: a dictionary, with following keys
    - 'wave': wavelength
    - 'flux': fluxes in 
    - 'fluxerr': flux error
    - 'tobs': time of observation
    - 'mjd': Median Julian Date

    """
    mjds = []
    names = ["wave", "flux", "fluxerr", "tobs", "mjd"]
    sn = {k: [] for k in names}
    id_to_obs = dict([])

    with open(base_name + ".SPECLIST.TEXT", "rb") as f:
        dat = f.readlines()
    for line in dat:
        temp = line.split()
        if len(temp) > 0 and b"VARNAMES:" in temp:
            varnames = [str(x.decode("utf-8")) for x in temp]
        else:
            the_id = int(temp[varnames.index("ID")])
            id_to_obs[the_id] = float(temp[varnames.index("TOBS")])
            mjds.append(float(temp[varnames.index("MJD")]))

    with open(base_name + ".SPECPLOT.TEXT", "rb") as f:
        dat = f.readlines()
    temp_id = None
    mjd_ind = 0

    for line in dat:
        temp = line.split()
        if len(temp) <= 0:
            continue
        if b"VARNAMES:" in temp:
            varnames = [str(x.decode("utf-8")) for x in temp]
            i_cid = varnames.index("CID")

        elif (b"OBS:" in temp and str(temp[i_cid].decode("utf-8")) in cid):

            if temp_id is None:
                temp_id = int(temp[varnames.index("ID")])

            if temp_id != int(temp[varnames.index("ID")]):
                mjd_ind += 1

            temp_id = int(temp[varnames.index("ID")])

            sn["flux"].append(float(temp[varnames.index("FLAM")]))
            sn["fluxerr"].append(float(temp[varnames.index("FLAMERR")]))
            sn["tobs"].append(id_to_obs[int(temp[varnames.index("ID")])])
            sn["mjd"].append(mjds[mjd_ind])
            sn["wave"].append(
                (
                    float(temp[varnames.index("LAMMAX")])
                    + float(temp[varnames.index("LAMMIN")])
                )
                / 2.0
            )

    sn = {k: np.array(sn[k]) for k in sn.keys()}

    return sn


def read_lc(cid, base_name, plotter_choice, tmin, tmax, filter_list):
    """Function to read a light curve.

    Parameters
    ----------


    Returns
    -------

    """
    peak = None
    fitted = False
    mintime = np.inf
    maxtime = -np.inf

    names = ["time", "flux", "fluxerr", "filter", "chi2"]
    sn = {k: [] for k in names}
    fit = {k: [] for k in ["time", "flux", "filter"]}

    with open(base_name + ".LCPLOT.TEXT", "rb") as f:
        dat = f.readlines()

    for line in dat:
        temp = line.split()
        if len(temp) <= 0:
            continue

        if b"VARNAMES:" in temp:
            varnames = [str(x.decode("utf-8")) for x in temp]
            i_cid = varnames.index("CID")

        elif (b"OBS:" in temp and str(temp[i_cid].decode("utf-8")) in cid):
            tCid = str(temp[i_cid].decode("utf-8"))
            tObs = float(temp[varnames.index("Tobs")])
            band = str(temp[varnames.index("BAND")].decode("utf-8"))
            dataflag = int(temp[varnames.index("DATAFLAG")])

            if not (tmin < tObs < tmax):
                continue
            if (filter_list is not None and band not in filter_list):
                continue

            if dataflag == 1:
                if peak is None:
                    peak = float(temp[varnames.index("MJD")]) - tObs
                if tObs > maxtime:
                    maxtime = tObs
                if tObs < mintime:
                    mintime = tObs
                sn["time"].append(tObs)
                sn["flux"].append(float(temp[varnames.index("FLUXCAL")]))
                sn["fluxerr"].append(float(temp[varnames.index("FLUXCAL_ERR")]))
                sn["filter"].append(band)
                sn["chi2"].append(float(temp[varnames.index("CHI2")]))

            elif dataflag == 0:
                fitted = True
                fit["time"].append(tObs)
                fit["flux"].append(float(temp[varnames.index("FLUXCAL")]))
                fit["filter"].append(band)

    if fitted and plotter_choice == 'salt2':
        with open(base_name+".FITRES.TEXT",'rb') as f:
            dat = f.readlines()
        for line in dat:
            temp = line.split()
            if len(temp) <= 0:
                continue
            if b'VARNAMES:' in temp:
                varnames = [str(x.decode('utf-8')) for x in temp]
                tCid = str(temp[varnames.index('CID')].decode('utf-8'))

            elif b'SN:' in temp and str(temp[i_cid].decode("utf-8")) in cid:
                datadict = {}
                datadict['NDOF'] = float(temp[varnames.index('NDOF')])
                datadict['FITCHI2'] = float(temp[varnames.index('FITCHI2')])
                datadict['x0'] = (
                        float(temp[varnames.index('x0')]),
                        float(temp[varnames.index('x0ERR')]),
                ) 
                datadict['x1'] = (
                        float(temp[varnames.index('x1')]),
                        float(temp[varnames.index('x1ERR')]),
                ) 
                datadict['c'] = (
                        float(temp[varnames.index('c')]),
                        float(temp[varnames.index('cERR')]),
                )
                fit['params'] = datadict
                break

    sn = {k: np.array(sn[k]) for k in sn.keys()}
    fit_arrays = {}
    for k, v in fit.items():
        if k == 'params':
            fit_arrays[k] = v
        else:
            fit_arrays[k] = np.array(v)
    fit = fit_arrays

    if len(fit['filter']) > 0:
        fits = {}
        tranges = {}
        for afilter in np.unique(fit['filter']):
            fits[afilter] = interp1d(
                fit['time'][fit['filter'] == afilter],
                fit['flux'][fit['filter'] == afilter]
            )
            tranges[afilter] = [
                np.min(fit['time'][fit['filter'] == afilter]),
                np.max(fit['time'][fit['filter'] == afilter])
            ]
        fits['trange'] = tranges

        if 'params' in fit.keys():
            fits['params'] = fit['params']
        else:
            fits['params'] = {}
    else:
        fits = []
    return(sn, fits, peak, (mintime, maxtime))


def read_fitres(fitres_filename, param):
    fit = {}
    with open(fitres_filename, "rb") as f:
        dat = f.readlines()
    for line in dat:
        temp = line.split()
        if len(temp) > 0 and b"VARNAMES:" in temp:
            varnames = [str(x.decode("utf-8")) for x in temp]
        elif len(temp) > 0 and b"SN:" in temp:
            datadict = {}

            datadict['NDOF'] = float(temp[varnames.index('NDOF')])
            datadict['FITCHI2'] = float(temp[varnames.index('FITCHI2')])
            datadict['x0'] = (
                    float(temp[varnames.index('x0')]),
                    float(temp[varnames.index('x0ERR')]),
            ) 
            datadict['x1'] = (
                    float(temp[varnames.index('x1')]),
                    float(temp[varnames.index('x1ERR')]),
            ) 
            datadict['c'] = (
                    float(temp[varnames.index('c')]),
                    float(temp[varnames.index('cERR')]),
            ) 

            cid = str(temp[varnames.index("CID")].decode("utf-8"))
            fit[cid] = datadict
            
            if param is not None:
                if param not in varnames:
                    raise RuntimeError(
                        "Parameter %s given for joint distribution but not found in FITRES file %s"
                        % (param, fitres_filename)
                    )
                fit[cid][param] = float(temp[varnames.index(param)])
    return fit


def read_snana(snana_filename, cid, param):
    if isinstance(cid, (list, tuple, np.ndarray)) and len(cid) == 1:
        cid = cid[0]

    with open(snana_filename + ".SNANA.TEXT", "r") as f:
        dat = f.readlines()
    for line in dat:
        temp = line.split()
        if len(temp) > 0 and "VARNAMES:" in temp:
            varnames = temp
        if len(temp) > 0 and "SN:" in temp and str(cid) in temp:
            try:
                ind = varnames.index(param)
            except:
                print(
                    "Either format of SNANA.TEXT file is wrong and no varnames, or %s not in file"
                    % param
                )
            return (temp[ind], temp[ind + 1])

    print("%s not found in SNANA.TEXT file." % (str(cid)))
    return


def plot_spec(cid, bin_size, base_name, noGrid, zname):
    sn = read_spec(cid, base_name)

    if len(sn["tobs"]) == 0:
        return []
    uniq_tobs = np.unique(sn["tobs"])

    if len(uniq_tobs) > 1:
        figs = []
        m = 0
        for nfig in range(math.ceil(len(uniq_tobs) / 4.0)):
            fig, ax = plt.subplots(
                nrows=min(len(uniq_tobs), 4),
                ncols=1,
                figsize=(8, 8),
                sharex=True,
            )
            try:
                ax[0].set_title('SNID=%s' % cid[0], fontsize=16)
            except:
                ax = [ax]
                ax[0].set_title('SNID=%s' % cid[0], fontsize=16)
            for j in range(min(len(uniq_tobs[m:]), 4)):
                temp_sn = np.where(sn["tobs"] == uniq_tobs[m])[0]
                if bin_size != 0:
                    #  bins = np.digitize(np.array(sn['wave']),
                    #   np.arange(sn['wave'][0], sn['wave'][-1], bin_size))
                    binned_wave = []
                    binned_flux = []
                    binned_fluxerr = []
                    bins = np.trunc(sn["wave"][temp_sn] / bin_size)
                    for i in np.unique(bins):
                        binned_wave = np.append(
                            binned_wave, 
                            np.mean(sn["wave"][temp_sn][bins == i])
                        )
                        binned_flux = np.append(
                            binned_flux, 
                            np.mean(sn["flux"][temp_sn][bins == i])
                        )
                        binned_fluxerr = np.append(
                            binned_fluxerr, 
                            np.mean(sn["fluxerr"][temp_sn][bins == i])
                        )
                else:
                    binned_wave = sn["wave"][temp_sn]
                    binned_flux = sn["flux"][temp_sn]
                    binned_fluxerr = sn["fluxerr"][temp_sn]
                    # sn = (
                    #   sn.group_by(np.trunc(sn['wave']/bin_size))
                    # ).groups.aggregate(np.mean)

                if np.unique(sn["mjd"])[j] < 0:
                    spec_label = "HOST"
                else:
                    spec_label = "SN:%.2f" % np.unique(sn["tobs"])[j]
                ax[j].plot(binned_wave, binned_flux, color="k", label=spec_label)
                ylim = ax[j].get_ylim()
                ax[j].fill_between(
                    binned_wave,
                    binned_flux - binned_fluxerr,
                    binned_flux + binned_fluxerr,
                    color="r",
                    alpha=0.3,
                    label=r"$1\sigma$ Error",
                )
                ax[j].plot([binned_wave[0], binned_wave[-1]],
                           [0, 0], "k--", alpha=0.5)
                ax[j].set_ylim(ylim)
                ax[j].legend(fontsize=12)

                ax[j].set_ylabel(r"$F_\lambda$", fontsize=16)
                if not noGrid:
                    ax[j].grid()
                m += 1

            ax[j].set_xlabel("Observer Frame Wavelength ($\AA$)", fontsize=16)

            figs.append(fig)
            plt.close()
    else:
        fig = plt.figure(figsize=(10, 8))

        if bin_size != 0:
            #  bins = np.digitize(np.array(sn['wave']),
            #   np.arange(sn['wave'][0], sn['wave'][-1], bin_size))     
            binned_wave = []
            binned_flux = []
            binned_fluxerr = []
            bins = np.trunc(sn['wave'] / bin_size)
            for i in np.unique(bins):
                binned_wave = np.append(
                    binned_wave, np.mean(sn['wave'][bins == i])
                    )
                binned_flux = np.append(
                    binned_flux, np.mean(sn['flux'][bins == i])
                    )
                binned_fluxerr = np.append(
                    binned_fluxerr, np.mean(sn['fluxerr'][bins == i])
                    )
        else:
            binned_wave = sn['wave']
            binned_flux = sn['flux']
            binned_fluxerr = sn['fluxerr']
            # sn = (sn.group_by(np.trunc(sn['wave']/bin_size))
            #       ).groups.aggregate(np.mean)
        ind = np.argsort(binned_wave)
        plt.plot(
            binned_wave[ind],
            binned_flux[ind],
            color="k",
            label="TOBS:%.2f" % np.unique(sn["tobs"])[0],
        )
        ylim = plt.ylim()
        plt.fill_between(
            binned_wave[ind],
            binned_flux[ind] - binned_fluxerr[ind],
            binned_flux[ind] + binned_fluxerr[ind],
            color="r",
            alpha=0.3,
            label=r"$1\sigma$ Error",
        )
        plt.plot([binned_wave[0], binned_wave[-1]], [0, 0], "k--", alpha=0.5)
        plt.ylim(ylim)
        plt.legend(fontsize=12)
        plt.xlabel("Observer Frame Wavelength ($\AA$)", fontsize=16)
        plt.ylabel("Flux", fontsize=16)
        plt.title("SN%s" % cid[0], fontsize=16)
        if not noGrid:
            plt.grid()
        figs = [fig]
        plt.close()
    # plt.savefig('SNANA_SPEC_%s.pdf'%'_'.join(cid),format='pdf',overwrite=True)
    return(figs)


def plot_lc(cid, base_name, noGrid, plotter_choice,
            tmin, tmax, filter_list, plot_all, zname):
    if tmin is None:
        tmin = -np.inf
    if tmax is None:
        tmax = np.inf

    sn, fits, peak, minmaxtime = read_lc(
        cid, base_name, plotter_choice, tmin, tmax, filter_list
    )
    z = read_snana(base_name, cid, zname)
    if len(sn["time"]) == 0:
        return [[], []]
    rows = int(math.ceil(len(np.unique(sn["filter"]))))
    figs = []
    all_bands = np.append(
        [x for x in __band_order__ if x in np.unique(sn["filter"])],
        [x for x in np.unique(sn["filter"]) if x not in __band_order__],
    )

    j = 0
    minx = np.min(sn["time"])
    maxx = np.max(sn["time"])
    if minx < 0:
        minx = min(minx * 1.1, minx - 5)
    else:
        minx = min(minx * 0.9, minx - 5)
    if maxx < 0:
        maxx = max(maxx * 0.9, maxx + 5)
    else:
        maxx = max(maxx * 1.1, maxx + 5)

    xlims = (minx, maxx)
    sharedx = True
    for nfig in range(int(math.ceil(rows / 4.0))):
        fig, ax = plt.subplots(
            nrows=min(len(all_bands), 4), ncols=1,
            figsize=(8, 8), sharex=sharedx
        )
        try:
            ax[0].set_title("SNID=%s" % cid[0], fontsize=16)
        except:
            ax = [ax]
            ax[0].set_title("SNID=%s" % cid[0], fontsize=16)
        fit_print = False

        for i in range(min(len(all_bands[j:]), 4)):
            idx = np.where(sn["filter"] == all_bands[j])[0]
            temp_sn = {k: sn[k][idx] for k in sn.keys()}
            chi2 = np.mean(temp_sn["chi2"])
            if chi2 > 0:
                lab = r"%s: $\chi^2_{red}$=%.1f" % (
                    all_bands[j],
                    np.mean(temp_sn["chi2"]),
                )
                leg_size = 10
            else:
                lab = all_bands[j]
                leg_size = 12

            ax[i].errorbar(
                temp_sn["time"],
                temp_sn["flux"],
                yerr=temp_sn["fluxerr"],
                fmt=".",
                markersize=8,
                color="k",
                label=lab,
            )
            if len(fits) > 0:
                if not plot_all:
                    fit_time = np.arange(
                        np.min(temp_sn["time"]) - 5,
                        np.max(temp_sn["time"]) + 5, 1)
                else:
                    fit_time = np.arange(
                        fits["trange"][all_bands[j]][0],
                        fits["trange"][all_bands[j]][1], 1)

                fit_time = fit_time[
                    np.where(
                        np.logical_and(np.logical_and(fit_time >= minx, fit_time <= maxx),
                                       np.logical_and(fit_time >= fits["trange"][all_bands[j]][0],
                                                      fit_time <= fits["trange"][all_bands[j]][1]))
                        )[0]
                ]
                ax[i].plot(
                    fit_time,
                    fits[all_bands[j]](fit_time),
                    color="r",
                    label="Best Fit",
                    linewidth=3,
                )
                if not fit_print:
                    to_print = []
                    for fit_key in fits["params"].keys():
                        if fit_key == "x0":
                            val, err = fits["params"][fit_key]
                            to_print.append(
                                ["$%s: %.2e" % (fit_key, val), "%.2e$\n" % err]
                            )
                        elif fit_key in ['x1', 'c']:
                            val, err = fits["params"][fit_key]
                            to_print.append(
                                ['$%s: %.2f' % (fit_key, val), '%.2f$\n' % err]
                            )
                        elif fit_key == 'NDOF':
                            chi2 = fits['params']['FITCHI2']
                            ndof = fits['params'][fit_key]
                            to_print.append(
                                'CHI2/NDOF: %.2f/%.2f\n' % (chi2, ndof)
                            )
                        else:
                            pass
                    if z is not None:
                        fz, fzerr = float(z[0]), float(z[1])
                        sz, szerr = str(z[0]), str(z[1])
                        if np.any([
                            d!='0' for d in szerr[
                                szerr.find('.')+1:szerr.find('.')+4
                                ]
                            ]):
                            print(szerr)
                            to_print.append('z: %.2f'%fz+r'$\pm$'+'%.3f'%fzerr)
                        else:
                            to_print.append('z: %.2f'%fz)

                    bits = []
                    for aval in to_print:
                        if isinstance(aval, list):
                            bits.append(aval[0]+r'\pm'+aval[1])
                        else:
                            bits.append(aval)
                    annotation = "".join(bits)
                    ax[i].annotate(annotation, xy=(0.015, 0.72),
                                   xycoords='axes fraction', fontsize=6)
                fit_print = True

            ax[i].legend(fontsize=leg_size)
            ax[i].set_ylabel('Flux', fontsize=16)

            if len(fits) > 0:
                try:
                    maxFlux = max(
                        np.max(temp_sn["flux"]),
                        np.max(fits[all_bands[j]](fit_time))
                    )
                except:
                    maxFlux = np.max(temp_sn["flux"])
            else:
                maxFlux = np.max(temp_sn["flux"])

            ax[i].set_ylim((-0.1 * np.max(temp_sn["flux"]), 1.2 * maxFlux))
            if not noGrid:
                ax[i].grid()
            j += 1
            # i+=1
        for k in range(i + 1, min(len(all_bands), 4)):
            fig.delaxes(ax[k])
        ax[i].tick_params(axis="x", labelbottom=True, bottom=True)
        ax[i].set_xlabel("MJD-%.2f" % peak, fontsize=16)
        ax[i].set_xlim(xlims)
        figs.append(fig)
        plt.close()

    # fig.text(0.5, 0.02, 'Time (Rest Frame Days)', ha='center',fontsize=16)
    # fig.text(0.04, .5, 'Flux', va='center', rotation='vertical',fontsize=16)
    # plt.savefig('SNANA_LC_%s.pdf'%'_'.join(cid),format='pdf',overwrite=True)

    return (figs, fits)


def plot_cmd(genversion, cid_list, nml, isdist, private):
    plotter = "normal"
    if nml is not None:
        if os.path.splitext(nml)[1].upper() != ".NML":
            nml = os.path.splitext(nml)[0] + ".NML"
        with open(nml, "r") as f:
            p = f.readlines()

        for line in p:
            if "FITMODEL_NAME" in line:
                if "SALT2" in line:
                    plotter = "salt2"

    rand = str(np.random.randint(10000, 100000))
    if private is not None:
        private_path = " PRIVATE_DATA_PATH %s" % private
    else:
        private_path = ""
    genversion += private_path
    
    if nml is not None:
        if cid_list is not None:
            cmd = (
                "snlc_fit.exe "
                + nml
                + " VERSION_PHOTOMETRY "
                + genversion
                + " SNCCID_LIST "
                + cid_list
                + " CUTWIN_CID 0 0 SNTABLE_LIST 'FITRES(text:key) SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX 'OUT_TEMP_"
                + rand
                + "' > OUT_TEMP_"
                + rand
                + ".LOG"
            )
        elif isdist:
            cmd = (
                "snlc_fit.exe "
                + nml
                + " VERSION_PHOTOMETRY "
                + genversion
                + " SNTABLE_LIST "
                + "'FITRES(text:key) SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX OUT_TEMP_"
                + rand
                + " > OUT_TEMP_"
                + rand
                + ".LOG"
            )
        else:
            cmd = (
                "snlc_fit.exe "
                + nml
                + " VERSION_PHOTOMETRY "
                + genversion
                + " MXEVT_PROCESS 5 SNTABLE_LIST "
                + "'FITRES(text:key) SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX OUT_TEMP_"
                + rand
                + " > OUT_TEMP_"
                + rand
                + ".LOG"
            )
    elif cid_list is None:
        cmd = (
            "snana.exe NOFILE VERSION_PHOTOMETRY "
            + genversion
            + " MXEVT_PROCESS 5 SNTABLE_LIST 'SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)'"
            + " TEXTFILE_PREFIX 'OUT_TEMP_"
            + rand
            + "' > OUT_TEMP_"
            + rand
            + ".LOG"
        )
    else:
        cmd = (
            "snana.exe NOFILE VERSION_PHOTOMETRY "
            + genversion
            + " SNCCID_LIST "
            + cid_list
            + " CUTWIN_CID 0 0 SNTABLE_LIST 'SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PRE\
FIX 'OUT_TEMP_"
            + rand
                + "' > OUT_TEMP_"
                + rand
            + ".LOG"
        )
    os.system(cmd)
    with open("OUT_TEMP_" + rand + ".LOG", "rb+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(
            b"SNANA COMMAND:\n\n"
            + bytearray(textwrap.fill(cmd, 80), encoding="utf-8")
            + b"\n"
            + content
        )
    if len(glob.glob("OUT_TEMP_" + rand + "*.TEXT")) == 0:
        print("There was an error in retrieving your SN")
        sys.exit()

    if cid_list is None:
        with open("OUT_TEMP_" + rand + ".LCLIST.TEXT", "rb") as f:
            all_dat = f.readlines()
        all_cids = []
        for line in all_dat:
            temp = line.split()
            if len(temp) > 0 and b"VARNAMES:" in temp:
                varnames = [str(x.decode("utf-8")) for x in temp]
            elif len(temp) > 0 and b"SN:" in temp:
                all_cids.append(str(temp[varnames.index("CID")].decode("utf-8")))
        all_cids = ",".join(all_cids)
    else:
        all_cids = cid_list
    return (plotter, "OUT_TEMP_" + rand, all_cids)


def read_existing(nml):
    files = glob.glob("OUT_TEMP_*")
    unq_ids = np.unique([f[f.rfind("_") + 1 : f.find(".")] for f in files])

    if len(unq_ids) > 1:
        print(
            "For existing files, must have only one set of OUT_TEMP_ files in directory..."
        )
        sys.exit(1)
    plotter = "normal"
    if nml is not None:
        if os.path.splitext(nml)[1].upper() != ".NML":
            nml = os.path.splitext(nml)[0] + ".NML"
        with open(nml, "r") as f:
            p = f.readlines()

        for line in p:
            if "FITMODEL_NAME" in line:
                if "SALT2" in line:
                    plotter = "salt2"
    base_name = "OUT_TEMP_" + unq_ids[0]
    with open(base_name + ".LCLIST.TEXT", "r") as f:
        dat = f.readlines()

    cids = []
    for line in dat:
        temp = line.split()
        if len(temp) > 0 and temp[0] == "SN:":
            cids.append(temp[1])

    with open(base_name + ".LOG", "r") as f:
        dat = f.readlines()
    for line in dat:
        temp = line.split()
        if len(temp) > 0 and "snlc" in line:
            split_line = line.split()
            genversion = split_line[split_line.index("VERSION_PHOTOMETRY") + 1]

    return (plotter, base_name, ",".join(cids), genversion)


def output_fit_res(fitres, filename):
    with open(os.path.splitext(filename)[0] + ".fitres", "w") as f:
        f.write("VARNAMES: CID x0 x0err x1 x1err c cerr\n")
        for cid in fitres.keys():
            f.write(
                "SN: %s %f %f %f %f %f %f\n"
                % (
                    cid,
                    fitres[cid]["x0"][0],
                    fitres[cid]["x0"][1],
                    fitres[cid]["x1"][0],
                    fitres[cid]["x1"][1],
                    fitres[cid]["c"][0],
                    fitres[cid]["c"][1],
                )
            )


def create_dists(fitres, param, joint_type):
    try:
        import seaborn as sns
    except:
        print("For distribution plotting, seaborn is needed.")
        sys.exit(1)  # really necessary to kill the whole VMachine ??

    res = {p: [] for p in ["x0", "x1", "c"]}
    reserr = {p: [] for p in ["x0", "x1", "c"]}
    if param is not None:
        res[param] = []

    for cid in fitres.keys():
        for p in ["x0", "x1", "c"]:
            try:
                res[p].append(fitres[cid][p][0])
                reserr[p].append(fitres[cid][p][1])
            except RuntimeError:
                print(f"Skipping {cid} for distributions...")
        if param is not None:
            res[param].append(fitres[cid][param])

    figs = []
    for p in ["x0", "x1", "c"]:
        if param is not None:
            mean_valx = np.mean(res[p])
            mean_valy = np.mean(res[param])
            std_valx = np.std(res[p])
            std_valy = np.std(res[param])
            ax = sns.jointplot(x=res[p], y=res[param], kind=joint_type)
            fig = plt.gcf()

            if joint_type in ["reg", "scatter"]:
                plt.errorbar(res[p], res[param], 
                             xerr=reserr[p], fmt=".", markersize=5)
            fig.set_size_inches(10, 8)

            if joint_type == "kde":
                ax.ax_marg_x.set_xlim(
                    mean_valx - 3 * std_valx, mean_valx + 3 * std_valx
                )
                ax.ax_marg_y.set_ylim(
                    mean_valy - 3 * std_valy, mean_valy + 3 * std_valy
                )
            ax.set_axis_labels(
                "%s Parameter" % p,
                "Simulated %s"
                % (" ".join([x[0] + x[1:].lower() for x in param.split("_")])),
                fontsize=16,
            )

        else:
            fig = plt.figure(figsize=(10, 8))
            plt.hist(res[p])
            plt.xlabel("%s Parameter" % p, fontsize=16)
            plt.ylabel("N SN", fontsize=16)
        figs.append(fig)

    return figs


def find_files(version, cid_list=[]):
    for dirpath, dirnames, filenames in os.walk(PATH):
        dirpath = Path(dirpath)
        for sub_dirname in dirnames:
            if sub_dirname == version:
                listpath = dirpath / sub_dirname / f"{version}.LIST"

                list_files = {}
                for x in np.loadtxt(str(listpath), dtype=str):
                    head = os.path.splitext(x)[0]
                    idx = head.rfind("_SN") + 3
                    val = head[idx:].lstrip("0")
                    if val in cid_list or len(cid_list) == 0:
                        list_files[val] = dirpath / sub_dirname / x
                return list_files


def main():
    import argparse

    DESC = "Plot utility for SNANA lightcurves and parameter distributions"
    EPIL = "This produces graphic plots in various formats"
    parser = argparse.ArgumentParser(description=DESC, epilog=EPIL)

    #parser = OptionParser()

    parser.add_argument(
        "-B",
        help="Spectra: Bin size for spectral plotting",
        action="store",
        type=float,
        dest="bin_size",
        default=0,
    )

    parser.add_argument(
        "--spec",
        help="LC and Spectra: Plot only spectra",
        action="store_true",
        dest="spec",
        default=False,
    )
    parser.add_argument(
        "--lc",
        help="LC and Spectra: Plot only LC",
        action="store_true",
        dest="lc",
        default=False,
    )
    parser.add_argument(
        "--nogrid",
        help="LC and Spectra: Do not add a grid to the plots.",
        action="store_true",
        dest="noGrid",
        default=False,
    )

    parser.add_argument(
        "-a",
        help="LC: tmin for lc plotting (Phase).",
        action="store",
        dest="tmin",
        type=float,
        default=None,
    )
    parser.add_argument(
        "-b",
        help="LC: tmax for lc plotting (Phase).",
        action="store",
        dest="tmax",
        type=float,
        default=None,
    )
    parser.add_argument(
        "-t",
        help="LC: comma separated list of filters to plot (default all)",
        action="store",
        dest="filt_list",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--plotDataRange",
        help="LC: plots data range only (default full fit range when present)",
        action="store_true",
        dest="plot_all2",
        default=False,
    )
    parser.add_argument(
        "--plotAll",
        help=SUPPRESS_HELP,
        action="store_true",
        dest="plot_all",
        default=False,
    )
    parser.add_argument(
        "-n",
        help="LC: comma separated list of variables to print to plot.",
        action="store",
        dest="plot_text",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-f",
        help="LC and Distributions: .NML filename",
        action="store",
        type=str,
        dest="nml_filename",
        default=None,
    )
    parser.add_argument(
        "--fitres",
        help="LC and Distributions: Output a file containing fit results.",
        action="store_true",
        dest="res_out",
        default=False,
    )

    parser.add_argument(
        "-F",
        help="Distributions: fitres filename, used to create distributions of SALT2 fitting parameters.",
        action="store",
        type=str,
        dest="fitres_filename",
        default=None,
    )
    parser.add_argument(
        "-s",
        help="Distributions: Name of parameter to view in joint distribution with SALT2 fitting parameters",
        action="store",
        type=str,
        dest="joint_param",
        default=None,
    )
    parser.add_argument(
        "-k",
        help="Distributions: Joint plot type (scatter,reg,resid,kde,hex)",
        action="store",
        type=str,
        dest="joint_type",
        default="kde",
    )
    parser.add_argument(
        "--dist",
        help="Distributions: Fit and then plot the distributions of fitting parameters.",
        action="store_true",
        dest="dist",
        default=False,
    )

    parser.add_argument(
        "-i",
        help="All: CID(s) as comma separated list or range (1-10)",
        action="store",
        type=str,
        dest="CID",
        default="None",
    )
    parser.add_argument(
        "-p",
        help="All: Private Data Path",
        action="store",
        type=str,
        dest="private_path",
        default=None,
    )
    parser.add_argument(
        "-v",
        help="All: Version",
        action="store",
        type=str,
        dest="version",
        default=None,
    )
    parser.add_argument(
        "-z",
        help="All: Name of redshift parameter to read",
        action="store",
        type=str,
        dest="zname",
        default="zCMB",
    )
    parser.add_argument(
        "--existing",
        help="All: Use existing output files for plotting (Must have only 1 set in directory)",
        action="store_true",
        dest="existing",
        default=False,
    )
    parser.add_argument(
        "--noclean",
        help="All: Leave intermediate files for debugging",
        action="store_true",
        dest="noclean",
        default=False,
    )
    parser.add_argument(
        "--silent",
        help="All: Do not print anything",
        action="store_true",
        dest="silent",
        default=False,
    )

    options = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    options.plot_all = options.plot_all or not options.plot_all2
    
    if not options.existing:
        if options.version is None:
            if options.fitres_filename is None:
                raise RuntimeError("Need to define genversion")
        if options.CID == "None":
            if options.dist or options.fitres_filename is not None:
                print(
                    """No CID given, assuming all for distributions, 
                    then first 5 for LC/SPEC plotting..."""
                )
            else:
                print("No CID given, assuming first 5...")
            options.CID = None
            all_cid = True
        elif "-" in options.CID:
            options.CID = ",".join(
                [
                    str(i)
                    for i in range(
                        int(options.CID[: options.CID.find("-")]),
                        int(options.CID[options.CID.find("-") + 1 :]) + 1,
                    )
                ]
            )
            all_cid = True
        else:
            all_cid = False
    else:
        all_cid = True

    if options.dist and options.nml_filename is None:
        raise RuntimeError(
            """If you use the 'dist' option, you must provide an 
            NML filename with the -f flag."""
        )

    if options.joint_type not in ["scatter", "reg", "resid", "kde", "hex"]:
        print("Joint plot type not recognized (see help), setting to kde")
        options.joint_type = "kde"
    if options.fitres_filename is None:
        if not options.existing:
            plotter_choice, options.base_name, options.CID = plot_cmd(
                options.version,
                options.CID,
                options.nml_filename,
                options.dist,
                options.private_path
            )
        else:
            (
                plotter_choice,
                options.base_name,
                options.CID,
                options.version,
            ) = read_existing(options.nml_filename)
        options.CID = options.CID.split(",")
        if options.filt_list is not None:
            options.filt_list = options.filt_list.split(",")
        filename = 'lcplot_' + options.version + ".pdf"
        num = 0
        if os.path.exists(filename):
            filename = os.path.splitext(filename)[0] + "_" + str(num) + ".pdf"
        while os.path.exists(filename):
            num += 1
            filename = os.path.splitext(filename)[0][:-1] + str(num) + ".pdf"
        if options.dist:
            fitres = read_fitres(
                options.base_name + ".FITRES.TEXT", options.joint_param
            )
            figs = create_dists(fitres, options.joint_param, options.joint_type)
        else:
            figs = []
        if all_cid:
            options.CID = options.CID[:5]
        with PdfPages(filename) as pdf:
            for f in figs:
                pdf.savefig(f)
            for cid in options.CID:
                if not options.silent:
                    print("Plotting SN %s" % cid)
                if options.spec:
                    figs = plot_spec(
                        [cid],
                        options.bin_size,
                        options.base_name,
                        options.noGrid,
                        options.zname,
                    )
                    for f in figs:
                        pdf.savefig(f)
                elif options.lc:
                    figs, fits = plot_lc(
                        [cid],
                        options.base_name,
                        options.noGrid,
                        plotter_choice,
                        options.tmin,
                        options.tmax,
                        options.filt_list,
                        options.plot_all,
                        options.zname,
                    )
                    for f in figs:
                        pdf.savefig(f)
                else:
                    figs = plot_spec(
                        [cid],
                        options.bin_size,
                        options.base_name,
                        options.noGrid,
                        options.zname,
                    )
                    for f in figs:
                        pdf.savefig(f)
                    figs, fits = plot_lc(
                        [cid],
                        options.base_name,
                        options.noGrid,
                        plotter_choice,
                        options.tmin,
                        options.tmax,
                        options.filt_list,
                        options.plot_all,
                        options.zname,
                    )
                    for f in figs:
                        pdf.savefig(f)

        if options.res_out:
            output_fit_res(fitres, filename)
        if not options.noclean:
            for x in glob.glob(options.base_name + "*"):
                os.remove(x)
    else:
        print("Creating distributions from FITRES file...")
        filename = (
            'lcplot_'+options.version + ".pdf"
            if options.version is not None
            else os.path.splitext(options.fitres_filename)[0] + ".pdf"
        )
        num = 0
        if os.path.exists(filename):
            filename = os.path.splitext(filename)[0] + "_" + str(num) + ".pdf"
        while os.path.exists(filename):
            num += 1
            filename = os.path.splitext(filename)[0][:-1] + str(num) + ".pdf"
        fitres = read_fitres(options.fitres_filename, options.joint_param)
        figs = create_dists(fitres, options.joint_param, options.joint_type)
        with PdfPages(filename) as pdf:
            for f in figs:
                pdf.savefig(f)
    if not options.silent:
        print("Done.")


if __name__ == "__main__":
    main()
