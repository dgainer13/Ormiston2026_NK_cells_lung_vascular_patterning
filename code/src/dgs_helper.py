import numpy as np
import rpy2.robjects as ro
import anndata2ri
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

def select_deviant_genes(adata, n_top=4000, layer="X"):
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv["sce"] = adata
    
    ro.globalenv["assay_name"] = layer
    ro.r('''
         suppressMessages(library(scry))
         sce <- devianceFeatureSelection(sce, assay = assay_name)
         deviance <- rowData(sce)$binomial_deviance
    ''')

    with localconverter(ro.default_converter + numpy2ri.converter):
        deviance = np.asarray(ro.globalenv["deviance"], dtype=float)

    deviance = np.nan_to_num(deviance, nan=0.0)
    adata.var["binomial_deviance"] = deviance

    idx = np.argsort(deviance)[-n_top:]
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask
    return adata