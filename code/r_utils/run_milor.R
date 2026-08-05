suppressPackageStartupMessages({
  library(miloR)
  library(SingleCellExperiment)
  library(tidyverse)
})

run_milo <- function(sce, k, d, prop, reduced_dim, sample, condition, ref) {
  milo <- Milo(sce)

  milo <- buildGraph(milo,
    k = k,
    d = d,
    reduced.dim = reduced_dim
  )

  milo <- makeNhoods(milo,
    prop = prop,
    k = k,
    d = d,
    reduced_dims = reduced_dim,
    refined = TRUE
  )

  meta <- as.data.frame(colData(milo))
  milo <- countCells(milo,
    meta.data = meta,
    samples = sample
  )

  milo_design <- data.frame(colData(milo))[, c(sample, condition)] |>
    distinct() |>
    remove_rownames() |>
    column_to_rownames(sample)
  
  milo_design[[condition]] <- relevel(factor(milo_design[[condition]]), ref = ref)

  da_res <- testNhoods(milo,
    design = reformulate(condition),
    design.df = milo_design,
    reduced.dim = reduced_dim,
    fdr.weighting = "graph-overlap"
  )

  milo <- buildNhoodGraph(milo)

  return(list(
    milo = milo,
    design = milo_design,
    da_res = da_res
  ))
}