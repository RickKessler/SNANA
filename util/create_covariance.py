import argparse
import logging
import shutil
from functools import reduce
import numpy as np
import os
import pandas as pd
from pathlib import Path
import re
import yaml
import sys
from sklearn.linear_model import LinearRegression
import seaborn as sb
import matplotlib.pyplot as plt


def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format="[%(levelname)8s |%(filename)21s:%(lineno)3d]   %(message)s")
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("seaborn").setLevel(logging.ERROR)


def read_yaml(path):
    logging.debug(f"Reading YAML from {path}")
    with open(path) as f:
        return yaml.safe_load(f.read())


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="the name of the yml config file to run.")
    parser.add_argument("-u", "--unbinned", help="Utilise individual SN instead of binning", action="store_true")
    parser.add_argument("-s", "--subtract_vpec", help="Subtract VPEC from MUERR. Forces unbinned", action="store_true")
    args = parser.parse_args()
    if args.subtract_vpec:
        args.unbinned = True
    return args


def load_data(path, args, config):
    if not os.path.exists(path):
        raise ValueError(f"Cannot load data from {path} - it doesnt exist")
    df = pd.read_csv(path, delim_whitespace=True, comment="#")
    logging.debug(f"\tLoaded data with shape {df.shape} from {path}")

    # Do a bit of data cleaning
    df = df.replace(999.0, np.nan)
    # M0DIF doesnt have MU column, so add it back in
    if "MU" not in df.columns:
        df["MU"] = df["MUREF"] + df["MUDIF"]
    if "MUERR" not in df.columns:
        df["MUERR"] = df["MUDIFERR"]

    # Sort to ensure direct subtraction comparison
    if "CID" in df.columns:
        df["CID"] = df["CID"].astype(str)
        df = df.sort_values(["zHD", "CID"])
        df = df.rename(columns={"zHD": "z", "MUMODEL": "MUREF"})
        if args.subtract_vpec:
            assert "MUERR_VPEC" in df.columns, f"Cannot subtract VPEC contribution as column 'MUERR_VPEC' doesn't exist in path {path}"
            df["MUERR"] = np.sqrt(df["MUERR"] ** 2 - df["MUERR_VPEC"] ** 2)
            logging.debug("Subtracted MUERR_VPEC from MUERR")
        elif config.get("CALIBRATORS"):
            calib_mask = df["CID"].isin(config.get("CALIBRATORS"))
            df.loc[calib_mask, "MUERR"] = np.sqrt(df.loc[calib_mask, "MUERR"] ** 2 - df.loc[calib_mask, "MUERR_VPEC"] ** 2)
        df = df.set_index(["IDSURVEY", "CID"])
    elif "z" in df.columns:
        df = df.sort_values("z")

    df = df.loc[:, ["z", "MU", "MUERR", "MUREF"]]
    return df


def get_data_files(folder, args, config):
    logging.debug(f"Loading all data files in {folder}")
    result = {}
    for file in sorted(os.listdir(folder)):
        if (not args.unbinned and ".M0DIF" in file) or (args.unbinned and ".FITRES" in file and "MUOPT" in file):
            label = file.replace(".gz", "").replace(".M0DIF", "").replace(".FITRES", "")
            result[label] = load_data(folder / file, args, config)

    if args.unbinned:
        result = get_common_set_of_sne(result)

    return result


def get_common_set_of_sne(datadict):

    combined = None
    for label, df in datadict.items():
        if combined is None:
            combined = df.index
        else:
            combined = combined.intersection(df.index)
        logging.debug(f"Common set from {label} has {combined.shape[0]} elements")
        assert combined.shape[0], "You have no common supernova. Please check this"
    logging.info(f"Common set of SN have {combined.shape[0]} events")
    for label, df in datadict.items():
        datadict[label] = df.loc[combined, :]

    return datadict


def get_fitopt_muopt_from_name(name):
    f = int(name.split("FITOPT")[1][:3])
    m = int(name.split("MUOPT")[1][:3])
    return f, m


def get_name_from_fitopt_muopt(f, m):
    return f"FITOPT{f:03d}_MUOPT{m:03d}"


def get_fitopt_scales(lcfit_info, sys_scales):
    """ Returns a dict mapping FITOPT numbers to (label, scale) """
    fitopt_list = lcfit_info["FITOPT_OUT_LIST"]
    result = {}
    for number, _, label, _ in fitopt_list:
        if label != "DEFAULT" and label is not None:
            if label not in sys_scales:
                logging.warning(f"No FITOPT scale found for {label}")
        scale = sys_scales.get(label, 1.0)
        d = int(number.replace("FITOPT", ""))
        result[d] = (label, scale)
    return result


def get_cov_from_diff(df1, df2, scale):
    """ Returns both the covariance contribution and summary stats (slope and mean abs diff) """
    assert df1["MU"].size == df2["MU"].size, "Oh no, looks like you have a different number of bins/supernova for your systematic and this is not good."
    diff = scale * ((df1["MU"] - df1["MUREF"]) - (df2["MU"] - df2["MUREF"])).to_numpy()
    diff[~np.isfinite(diff)] = 0
    cov = diff[:, None] @ diff[None, :]

    # Determine the gradient using simple linear regression
    reg = LinearRegression()
    weights = 1 / np.sqrt(0.003 ** 2 + df1["MUERR"].to_numpy() ** 2 + df2["MUERR"].to_numpy() ** 2)
    mask = np.isfinite(weights)
    reg.fit(df1.loc[mask, ["z"]], diff[mask], sample_weight=weights[mask])
    coef = reg.coef_[0]

    mean_abs_deviation = np.average(np.abs(diff), weights=weights)
    max_abs_deviation = np.max(np.abs(diff))
    return cov, (coef, mean_abs_deviation, max_abs_deviation)


def get_contributions(m0difs, fitopt_scales, muopt_labels, muopt_scales):
    """ Gets a dict mapping 'FITOPT_LABEL|MUOPT_LABEL' to covariance)"""
    result, slopes = {}, []

    for name, df in m0difs.items():
        f, m = get_fitopt_muopt_from_name(name)
        logging.debug(f"Determining contribution for FITOPT {f}, MUOPT {m}")

        # Get label and scale for FITOPTS and MUOPTS. Note 0 index is DEFAULT
        if f == 0:
            fitopt_label = "DEFAULT"
            fitopt_scale = 1.0
        else:
            fitopt_label, fitopt_scale = fitopt_scales[f]
        muopt_label = muopt_labels[m] if m else "DEFAULT"
        muopt_scale = muopt_scales.get(muopt_label, 1.0)
        logging.debug(f"FITOPT {f} has scale {fitopt_scale}, MUOPT {m} has scale {muopt_scale}")

        # Depending on f and m, compute the contribution to the covariance matrix
        if f == 0 and m == 0:
            # This is the base file, so don't return anything. CosmoMC will add the diag terms itself.
            cov = np.zeros((df["MU"].size, df["MU"].size))
            summary = 0, 0, 0
        elif m > 0:
            # This is a muopt, to compare it against the MUOPT000 for the same FITOPT
            df_compare = m0difs[get_name_from_fitopt_muopt(f, 0)]
            cov, summary = get_cov_from_diff(df, df_compare, fitopt_scale * muopt_scale)
        else:
            # This is a fitopt with MUOPT000, compare to base file
            df_compare = m0difs[get_name_from_fitopt_muopt(0, 0)]
            cov, summary = get_cov_from_diff(df, df_compare, fitopt_scale)

        result[f"{fitopt_label}|{muopt_label}"] = cov
        slopes.append([name, fitopt_label, muopt_label, *summary])

    summary_df = pd.DataFrame(slopes, columns=["name", "fitopt_label", "muopt_label", "slope", "mean_abs_deviation", "max_abs_deviation"])
    summary_df = summary_df.sort_values(["slope", "mean_abs_deviation", "max_abs_deviation"], ascending=False)
    return result, summary_df


def apply_filter(string, pattern):
    """ Used for matching COVOPTs to FITOPTs and MUOPTs"""
    string, pattern = string.upper(), pattern.upper()
    if pattern.startswith("="):
        return string == pattern[1:]
    elif pattern.startswith("+"):
        return pattern[1:] in string
    elif pattern.startswith("-"):
        return pattern[1:] not in string
    elif pattern == "":
        return True
    else:
        raise ValueError(f"Unable to parse COVOPT matching pattern {pattern}")


def get_cov_from_covopt(covopt, contributions, base, calibrators):
    # Covopts will come in looking like "[cal] [+cal,=DEFAULT]"
    # We have to parse this. Eventually can make this structured and move away from
    # legacy, but dont want to make too many people change how they are doing things
    # in one go
    label, fitopt_filter, muopt_filter = re.findall(r"\[(.*)\] \[(.*),(.*)\]", covopt)[0]
    fitopt_filter = fitopt_filter.strip()
    muopt_filter = muopt_filter.strip()
    logging.debug(f"Computing cov for COVOPT {label} with FITOPT filter '{fitopt_filter}' and MUOPT filter '{muopt_filter}'")

    final_cov = None

    if calibrators:
        mask_calib = base.reset_index()["CID"].isin(calibrators)

    for key, cov in contributions.items():
        fitopt_label, muopt_label = key.split("|")
        if apply_filter(fitopt_label, fitopt_filter) and apply_filter(muopt_label, muopt_filter):
            if final_cov is None:
                final_cov = cov.copy()
            else:
                # If we have calibrators and this is a VPEC term, filter out calib
                if calibrators and (apply_filter(fitopt_label, "+VPEC") or apply_filter(muopt_label, "+VPEC")):
                    cov2 = cov.copy()
                    cov2[mask_calib, :] = 0
                    cov2[:, mask_calib] = 0
                    final_cov += cov2
                else:
                    final_cov += cov

    assert final_cov is not None, f"No systematics matched COVOPT {label} with FITOPT filter '{fitopt_filter}' and MUOPT filter '{muopt_filter}'!"

    # Validate that the final_cov is invertible
    try:
        # CosmoMC will add the diag terms, so lets do it here and make sure its all good
        effective_cov = final_cov + np.diag(base["MUERR"] ** 2)

        # First just try and invert it to catch singular matrix errors
        precision = np.linalg.inv(effective_cov)

        # Then check that the matrix is well conditioned to deal with float precision
        episolon = sys.float_info.epsilon
        cond = np.linalg.cond(effective_cov)
        assert cond < 1 / episolon, "Covariance matrix is ill-conditioned and cannot be stably inverted"
        logging.info(f"Condition for covariance COVOPT {label} is {cond}")

        # Finally, re-invert the precision matrix and ensure its within tolerance of the original covariance
        cov2 = np.linalg.inv(precision)
        assert np.all(np.isclose(effective_cov, cov2)), "Double inversion does not give original covariance, matrix is unstable"

    except np.linalg.LinAlgError as ex:
        logging.exception(f"Unable to invert covariance matrix for COVOPT {label}")
        raise ex

    return label, final_cov


def write_dataset(path, data_file, cov_file, template_path):
    with open(template_path) as f_in:
        with open(path, "w") as f_out:
            f_out.write(f_in.read().format(data_file=data_file, cov_file=cov_file))


def write_data(path, base, cosmomc_format=True):
    if not cosmomc_format:
        base[["z", "MU", "MUERR"]].to_csv(path, sep=" ", index=True, float_format="%.5f")
        return

    zs = base["z"].to_numpy()
    mu = base["MU"].to_numpy()
    mbs = -19.36 + mu
    mbes = base["MUERR"].to_numpy()

    # I am so sorry about this, but CosmoMC is very particular
    logging.info(f"Writing out data to {path}")
    with open(path, "w") as f:
        f.write("#name zcmb    zhel    dz mb        dmb     x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor\n")
        for i, (z, mb, mbe) in enumerate(zip(zs, mbs, mbes)):
            f.write(f"{i:5d} {z:6.5f} {z:6.5f} 0  {mb:8.5f} {mbe:8.5f} 0 0 0 0 0 0 0 0 0 0 0 0\n")


def write_covariance(path, cov):
    logging.info(f"Writing covariance to {path}")

    # Write out the slopes
    with open(path, "w") as f:
        f.write(f"{cov.shape[0]}\n")
        for c in cov.flatten():
            f.write(f"{c:0.8f}\n")


def write_cosmomc_output(config, covs, base):
    # Copy INI files. Create covariance matrices. Create .dataset. Modifying INI files to point to resources
    out = Path(config["OUTDIR"]) / "cosmomc"
    data_file = out / f"data.txt"
    data_file_wCID = out / f"data_wCID.txt"
    dataset_template = Path(config["COSMOMC_TEMPLATES"]) / config["DATASET_FILE"]
    dataset_files = []
    os.makedirs(out, exist_ok=True)

    # Create lcparam file
    write_data(data_file, base)

    # Create supplementary file for people to merge covmat with fitres
    write_data(data_file_wCID, base, cosmomc_format=False)

    # Create covariance matrices and datasets
    for i, (label, cov) in enumerate(covs):
        cov_file = out / f"sys_{i}.txt"
        dataset_file = out / f"dataset_{i}.txt"

        write_covariance(cov_file, cov)
        write_dataset(dataset_file, data_file, cov_file, dataset_template)
        dataset_files.append(dataset_file)

    # Copy some INI files
    ini_files = [f for f in os.listdir(config["COSMOMC_TEMPLATES"]) if f.endswith(".ini") or f.endswith(".yml") or f.endswith(".md")]
    for ini in ini_files:
        op = Path(config["COSMOMC_TEMPLATES"]) / ini

        if ini in ["base.ini"]:
            # If its the base.ini, just copy it
            npath = out / ini
            shutil.copy(op, npath)
        else:
            # Else we need one of each ini per covopt
            for i, _ in enumerate(covs):
                # Copy with new index
                npath = out / ini.replace(".ini", f"_{i}.ini")
                shutil.copy(op, npath)

                basename = os.path.basename(npath).replace(".ini", "")

                # Append the dataset info
                with open(npath, "a+") as f:
                    f.write(f"\nfile_root={basename}\n")
                    f.write(f"jla_dataset={dataset_files[i]}\n")


def write_summary_output(config, covariances, base):
    out = Path(config["OUTDIR"])
    info = {}
    cov_info = {}
    for i, (label, cov) in enumerate(covariances):
        cov_info[i] = label
    info["COVOPTS"] = cov_info

    logging.info("Writing INFO.YML")
    with open(out / "INFO.YML", "w") as f:
        yaml.safe_dump(info, f)


def write_correlation(path, label, base_cov, diag, base):
    logging.debug(f"\tWriting out cov for COVOPT {label}")

    zs = base["z"].round(decimals=5)
    cov = diag + base_cov

    diag = np.sqrt(np.diag(cov))
    corr = cov / (diag[:, None] @ diag[None, :])

    np.fill_diagonal(corr, 1.0)
    np.savetxt(path, corr, fmt="%5.2f")
    np.savetxt(str(path).replace("corr", "cov"), cov, fmt="%9.6f")
    covdf = pd.DataFrame(cov, columns=zs, index=zs)
    base_covdf = pd.DataFrame(base_cov, columns=zs, index=zs)
    corr = pd.DataFrame((corr * 100).astype(int), columns=zs, index=zs)
    precision = pd.DataFrame(np.arcsinh(np.linalg.inv(cov)), columns=zs, index=zs)

    if corr.shape[0] < 1000:
        logging.debug("\tCreating precision and correlation matrix plots. Sit tight.")
        fig, axes = plt.subplots(figsize=(16, 14), ncols=2, nrows=2)
        axes = axes.flatten()
        annot = corr.shape[0] < 30

        args = {"shrink": 0.9}
        cm = np.nanmax(np.abs(cov))
        pm = np.nanmax(np.abs(precision.to_numpy()))
        bcm = np.nanmax(np.abs(base_covdf.to_numpy()))
        sb.heatmap(covdf.replace([np.inf, -np.inf], np.nan), annot=False, ax=axes[0], vmin=-cm, vmax=cm, cmap="PuOr", square=True, cbar_kws=args)
        sb.heatmap(precision, annot=False, ax=axes[1], cmap="PuOr", square=True, vmin=-pm, vmax=pm, cbar_kws=args)
        sb.heatmap(corr, annot=annot, fmt="d", ax=axes[2], cmap="RdBu", vmin=-100, vmax=100, square=True, cbar_kws=args)
        sb.heatmap(base_covdf.replace([np.inf, -np.inf], np.nan), annot=False, ax=axes[3], cmap="PuOr", vmin=-bcm, vmax=bcm, square=True, cbar_kws=args)
        axes[0].set_title("Covariance matrix")
        axes[1].set_title(f"Arcsinh(precision) ")
        axes[2].set_title("Correlation matrix (percent)")
        axes[3].set_title("Covariance matrix without diagonal contributions")
        plt.tight_layout()
        fig.savefig(path.with_suffix(".png"), bbox_inches="tight", dpi=300)
    else:
        logging.info("\tMatrix is large, skipping plotting")


def write_debug_output(config, covariances, base, summary):
    # Plot correlation matrix
    out = Path(config["OUTDIR"])

    # The slopes can be used to figure out what systematics have largest impact on cosmology
    logging.info("Writing out summary.csv information")
    with open(out / "summary.csv", "w") as f:
        with pd.option_context("display.max_rows", 100000, "display.max_columns", 100, "display.width", 1000):
            f.write(summary.__repr__())

    logging.info("Showing correlation matrices:")
    diag = np.diag(base["MUERR"] ** 2)
    for i, (label, cov) in enumerate(covariances):
        write_correlation(out / f"corr_{i}_{label}.txt", label, cov, diag, base)


def get_lcfit_info(submit_info):
    path = Path(submit_info["INPDIR_LIST"][0]) / "SUBMIT.INFO"
    logging.info(f"Loading LCFIT SUBMIT.INFO from {path}")
    return read_yaml(path)


def filter_nans(data):
    base_name = get_name_from_fitopt_muopt(0, 0)
    base = data[base_name]
    mask = ~(base.isnull().any(axis=1))
    for name, df in data.items():
        data[name] = df.loc[mask, :]
    return data, base.loc[mask, :]


def create_covariance(config, args):
    # Define all our pathing to be super obvious about where it all is
    input_dir = Path(config["INPUT_DIR"])
    data_dir = input_dir / config["VERSION"]
    sys_file = Path(config["SYSFILE"])

    # Read in all the needed data
    submit_info = read_yaml(input_dir / "SUBMIT.INFO")
    sys_scale = read_yaml(sys_file)
    fitopt_scales = get_fitopt_scales(submit_info, sys_scale)
    # Also need to get the FITOPT labels from the original LCFIT directory
    muopt_labels = {int(x.replace("MUOPT", "")): l for x, l, _ in submit_info.get("MUOPT_OUT_LIST", [])}
    muopt_scales = config["MUOPT_SCALES"]

    # Load in all the data
    data = get_data_files(data_dir, args, config)

    # Filter data to remove rows with infinite error
    data, base = filter_nans(data)

    # Now that we have the data, figure out how each much each FITOPT/MUOPT pair contributes to cov
    contributions, summary = get_contributions(data, fitopt_scales, muopt_labels, muopt_scales)

    # For each COVOPT, we want to find the contributions which match to construct covs for each COVOPT
    logging.info("Computing covariance for COVOPTS")
    covopts = ["[ALL] [,]"] + config["COVOPTS"]  # Adds the covopt to compute everything
    covariances = [get_cov_from_covopt(c, contributions, base, config.get("CALIBRATORS")) for c in covopts]

    write_cosmomc_output(config, covariances, base)
    write_summary_output(config, covariances, base)
    write_debug_output(config, covariances, base, summary)


if __name__ == "__main__":
    try:
        setup_logging()
        args = get_args()
        config = read_yaml(args.input)
        create_covariance(config, args)
    except Exception as e:
        logging.exception(e)
        raise e
