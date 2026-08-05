import scanpy as sc
import rpy2.robjects as ro
from pathlib import Path
from rpy2.robjects import numpy2ri
from anndata2ri import scipy2ri
from rpy2.robjects.conversion import localconverter


def get_soupx_groups(adata):
    pp = adata.copy()
    sc.pp.normalize_total(pp, target_sum=1e4)
    sc.pp.log1p(pp)
    sc.tl.pca(pp, svd_solver="arpack")
    sc.pp.neighbors(pp)
    sc.tl.leiden(pp, key_added="soupx_groups", flavor="igraph", n_iterations=2)
    return pp.obs["soupx_groups"]


def prepare_broth(adata, raw_dir):
    sample_id = adata.obs["sample_id"].iloc[0]
    sex = adata.obs["sex"].iloc[0]
    genotype = adata.obs["genotype"].iloc[0]

    raw_file = Path(raw_dir) / f"{sample_id}_{sex}_{genotype}_sample_raw_feature_bc_matrix.h5"
    raw_adata = sc.read_10x_h5(raw_file)
    raw_adata.var_names_make_unique()

    common = adata.var_names[adata.var_names.isin(raw_adata.var_names)]
    adata = adata[:, common].copy()
    raw_adata = raw_adata[:, common].copy()

    data = adata.X.T          
    raw = raw_adata.X.T        
    genes = adata.var_names
    cells = adata.obs_names
    soupx_groups = get_soupx_groups(adata)

    return adata, data, raw, genes, cells, soupx_groups


def cook_soup(adata, raw_dir):
    if "make_soup" not in ro.globalenv:
        raise RuntimeError(
            "make_soup() is not defined in R. Run the %%R cell that defines it first."
        )

    adata, data, raw, genes, cells, soupx_groups = prepare_broth(adata, raw_dir)
    make_soup = ro.globalenv["make_soup"]

    cv = ro.default_converter + numpy2ri.converter + scipy2ri.converter
    with localconverter(cv):
        out = make_soup(
            data, raw,
            ro.StrVector(genes.tolist()),
            ro.StrVector(cells.tolist()),
            ro.StrVector(soupx_groups.tolist()),
        )
        out = ro.conversion.get_conversion().rpy2py(out)   

    adata.layers["raw_counts"] = adata.X.copy()
    adata.layers["soupX_counts"] = out.T                  
    adata.X = adata.layers["soupX_counts"].copy()
    return adata
