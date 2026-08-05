suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(edgeR)
  library(GSVA)
  library(limma)
  library(dplyr)
})

run_gsva <- function(sce, gene_sets,
                     celltype_col = "celltype",
                     sample_col = "sample_id",
                     genotype_col = "genotype",
                     assay_type = "counts") {

  gsva_results <- list()
  de_results   <- list()
  fc_list      <- list()

  cell_types <- unique(as.character(colData(sce)[[celltype_col]]))

  for (cell_type in cell_types) {

    sce_ct <- sce[, as.character(colData(sce)[[celltype_col]]) == cell_type]

    agg <- aggregateAcrossCells(
      sce_ct,
      ids = colData(sce_ct)[[sample_col]],
      use.assay.type = assay_type
    )

    counts  <- assay(agg, assay_type)
    samples <- data.frame(
      sample   = as.character(colData(agg)[[sample_col]]),
      genotype = as.character(colData(agg)[[genotype_col]])
    )
    colnames(counts)  <- samples$sample
    rownames(samples) <- samples$sample
    stopifnot(identical(colnames(counts), rownames(samples)))

    dge    <- calcNormFactors(DGEList(counts = counts))
    logcpm <- cpm(dge, log = TRUE, prior.count = 1)

    gsvapar     <- gsvaParam(logcpm, geneSets = gene_sets)
    gsva_scores <- gsva(gsvapar)
    gsva_results[[cell_type]] <- list(gsva_scores = gsva_scores,
                                      samples = samples)

    genotype <- factor(samples$genotype)

    if (!all(c("KO", "WT") %in% levels(genotype))) {
      message(sprintf(
        "Skipping '%s': need both KO and WT (found: %s).",
        cell_type, paste(levels(genotype), collapse = ", ")
      ))
      next
    }

    design <- model.matrix(~ 0 + genotype)
    colnames(design) <- levels(genotype)

    fit  <- lmFit(gsva_scores, design)
    fit2 <- contrasts.fit(fit, makeContrasts(KOvsWT = KO - WT, levels = design))
    fit2 <- eBayes(fit2)

    de_result <- topTable(fit2, number = Inf, adjust.method = "bonferroni")

    de_results[[cell_type]] <- de_result
    fc_list[[cell_type]] <- data.frame(
      CellType  = cell_type,
      GeneSet   = rownames(de_result),
      LogFC     = de_result$logFC,
      adj.P.Val = de_result$adj.P.Val
    )
  }

  fc_df_all <- bind_rows(fc_list)

  list(
    gsva_results = gsva_results,
    de_results   = de_results,
    fc_df        = fc_df_all
  )
}
