import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import median_abs_deviation as mad


mad_metrics = [
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
    "pct_counts_in_top_20_genes",
]

qc_metrics = mad_metrics + ["pct_counts_mt"]


def assign_adatas(data_dir):
    adatas = []
    for file in sorted(Path(data_dir).glob("*.h5")):
        sample_id, sex, genotype = file.name.split("_")[:3]

        adata = sc.read_10x_h5(file)
        adata.obs["sample_id"] = sample_id
        adata.obs["sex"] = sex
        adata.obs["genotype"] = genotype
        adata.obs["mat_type"] = "filtered"
        adata.obs.index = adata.obs.index + "_" + sample_id
        adata.var_names_make_unique()

        adatas.append(adata)

    return adatas

def calculate_qc_metrics(adata):
    adata.var["mt"] = adata.var_names.str.startswith("mt-")  
    adata.var["hb"] = adata.var_names.str.contains(r"^Hb[ab]") 
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "hb"],
        percent_top=[20], log1p=True, inplace=True,
    )
    return adata


def is_mad_outlier(series, nmads):
    median = np.median(series)
    deviation = nmads * mad(series)

    return (series < (median - deviation)) | (series > (median + deviation))


def flag_outliers(adata, metrics=mad_metrics, nmads=5, pct_mt_max=20):
    agg = np.zeros(adata.n_obs, dtype=bool)
    for metric in metrics:
        mask = is_mad_outlier(adata.obs[metric], nmads)
        adata.obs[f"{metric}_outlier"] = mask
        agg |= np.asarray(mask)

    if pct_mt_max is not None:
        mask = adata.obs["pct_counts_mt"] > pct_mt_max
        adata.obs["pct_counts_mt_outlier"] = mask
        agg |= np.asarray(mask)

    adata.obs["outlier"] = agg
    return adata


def plot_qc_metrics(adatas, metrics=qc_metrics, groupby="sample_id", hue=None, palette=None):
    obs = pd.concat([ad.obs for ad in adatas]) if isinstance(adatas, (list, tuple)) else adatas.obs

    n_samples = obs[groupby].nunique()
    fig, axes = plt.subplots(
        len(metrics), 1, sharex=True,
        figsize=(0.45 * n_samples + 2, 2.6 * len(metrics)),
    )
    for ax, metric in zip(np.atleast_1d(axes), metrics):
        sns.violinplot(
            data=obs, x=groupby, y=metric,
            hue=hue or groupby,
            split=hue is not None,
            palette=palette or (sns.color_palette("husl" if hue is None else None)),
            inner="quartile",
            density_norm="width", cut=0, ax=ax,
        )
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    return fig


def filter_outliers(adata):
    if "outlier" not in adata.obs:
        raise ValueError("run flag_outliers(adata) first")
    discard = adata.obs["outlier"].to_numpy().copy()
    if "predicted_doublet" in adata.obs:
        discard |= adata.obs["predicted_doublet"].to_numpy()
    n_before = adata.n_obs
    filtered = adata[~discard].copy()
    filtered.uns["n_removed_cells"] = n_before - filtered.n_obs
    return filtered
