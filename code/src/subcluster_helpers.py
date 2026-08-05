import numpy as np
import scanpy as sc
import seaborn as sns
import harmonypy as hm
from matplotlib import pyplot as plt

from src.dgs_helper import select_deviant_genes


def deviant_hvg(adata, n_top=2000, layer="counts"):
    adata = select_deviant_genes(adata, n_top=n_top, layer=layer)
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    return adata


def plot_deviance(adata, figsize=(15, 4)):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharey=False)
    dev = np.sort(adata.var["binomial_deviance"])[::-1]
    ax1.plot(dev)
    ax1.axvline(adata.var["highly_deviant"].sum(), color="red", linestyle="--")
    ax1.set(xlabel="Gene rank", ylabel="Binomial deviance", title="Deviance of genes")

    adata.var["mean_counts"] = np.asarray(adata.layers["soupX_counts"].mean(axis=0)).ravel()
    sns.scatterplot(
        data=adata.var,
        x="mean_counts", y="binomial_deviance",
        hue="highly_deviant", s=6, linewidth=0,
        palette={True: "crimson", False: "gray"}, ax=ax2,
    )
    ax1.set_yscale("log")
    ax2.set_xscale("log")


def harmony_embed(adata, batch_key="sample_id", n_comps=50):
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=True, svd_solver="arpack")
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, vars_use=[batch_key])
    adata.obsm["X_pca_harmony"] = ho.Z_corr


def cluster_subset(adata, n_pcs, resolution, min_dist=0.5, random_state=1):
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=n_pcs, random_state=random_state)
    sc.tl.umap(adata, maxiter=500, min_dist=min_dist, random_state=random_state)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph",
                 n_iterations=2, random_state=random_state)


def apply_celltypes(adata, celltypes, key="leiden", out="celltype"):
    cluster_to_label = {cl: label for label, clusters in celltypes.items() for cl in clusters}
    adata.obs[out] = adata.obs[key].map(cluster_to_label).astype("category")
    assert not adata.obs[out].isna().any()
