import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from anndata2ri import scipy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri

def flag_doublets(adata):
    if "find_doublets" not in ro.globalenv:
        raise RuntimeError("Run %%R cell that defines find_doublets()")
    counts = adata.X.T
    cv = (ro.default_converter + numpy2ri.converter
          + scipy2ri.converter + pandas2ri.converter)
    with localconverter(cv):
        df = ro.globalenv["find_doublets"](
            counts,
            ro.StrVector(adata.var_names.tolist()),
            ro.StrVector(adata.obs_names.tolist()),
        )
    adata.obs["scDblFinder_score"] = df["score"].values
    adata.obs["predicted_doublet"] = (df["class"] == "doublet").values
    return adata