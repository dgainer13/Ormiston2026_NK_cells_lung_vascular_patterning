suppressPackageStartupMessages({
  library(gprofiler2)
  library(tidyverse)
})

symbol_to_ensembl <- function(object,
                             query,
                             organism = "mmusculus",
                             target = "ENSG",
                             filter_na = TRUE) {
  mapping <- gconvert(query,
    organism = organism,
    target = target,
    filter_na = filter_na
  )

  mapping <- mapping[!duplicated(mapping$input), ]

  ens <- mapping$target[match(rownames(object), mapping$input)]
  keep <- !is.na(ens) & !duplicated(ens)

  rowData(object)$gene <- rownames(object)
  object <- object[keep, ]
  rownames(object) <- ens[keep]
  

  return(object)
}
