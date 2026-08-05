h5ad_to_seurat <- function(h5ad_path, counts_layer = "counts") {
  sce <- zellkonverter::readH5AD(h5ad_path, reader = "R")

  rd <- SingleCellExperiment::reducedDimNames(sce)
  if ("X_umap" %in% rd) {
    rd[rd == "X_umap"]  <- "UMAP"
    SingleCellExperiment::reducedDimNames(sce) <- rd
  }

  seurat <- Seurat::as.Seurat(sce, counts = counts_layer, data = "X")

  if ("UMAP" %in% Seurat::Reductions(seurat)) {
    Seurat::Key(seurat[["UMAP"]]) <- "UMAP_"
  }

  seurat
}