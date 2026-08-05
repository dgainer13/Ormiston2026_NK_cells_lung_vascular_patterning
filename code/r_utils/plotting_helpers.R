suppressPackageStartupMessages({
  library(EnhancedVolcano)
  library(ggplot2)
  library(data.table)
  library(tidyverse)
})

mm_lineage_colors <- c(
  "immune" = "#56c5cc",
  "endothelial" = "#f14f7b",
  "epithelial" = "#f7aa1c",
  "mesenchymal" = "#33658A"
)

mm_immune_colors <- c(
  "alv_macs" = "#56c5cc",
  "b cells" = "#f14f7b",
  "basophils" = "#f7aa1c",
  "cdcs" = "#33658A",
  "classical monocytes" = "#6B4C9A",
  "gd_t_cells" = "#C03A83",
  "ilc2s" = "#FFB347",
  "ilc3s" = "#FF8C8C",
  "inter_macs" = "#FDB863",
  "mast cells/megakaryocytes" = "#779ECB",
  "migratory dcs" = "#B19CD9",
  "neutrophils" = "#80ED99",
  "nks" = "#6EE7B7",
  "nonclassical monocytes" = "#F7A072",
  "pdcs" = "#C08497",
  "prolif alv_macs" = "#A1983A",
  "prolif b cells" = "#9E5E9F",
  "prolif inter_macs" = "#7FB3D5",
  "t cells" = "#9184FF"
)

mm_nk_ilc_colors <- c(
  "nks" = "#56c5cc",
  "ilc2s" = "#f14f7b",
  "ilc3s" = "#f7aa1c"
)

mm_endo_colors <- c(
  "gen_cap" = "#56c5cc",
  "aerocytes" = "#f14f7b",
  "prolif gen_cap" = "#f7aa1c",
  "arterial" = "#33658A",
  "venous" = "#6B4C9A",
  "lymphatic" = "#C03A83"
)

mm_epithelial_colors <- c(
  "at1" = "#56c5cc",
  "at2" = "#f14f7b",
  "ciliated" = "#f7aa1c",
  "club" = "#33658A",
  "prolif at2" = "#6B4C9A",
  "transitional" = "#C03A83"
)

hs_immune_colors <- c(
  "T Cell" = "#56c5cc",
  "NK Cell" = "#f14f7b",
  "NKT Cell" = "#f7aa1c",
  "Alveolar macrophage" = "#33658A",
  "Monocyte" = "#6B4C9A",
  "pDC" = "#C03A83",
  "cDC" = "#80ED99",
  "B Cell" = "#FF8C8C",
  "Mast cell" = "#FDB863",
  "Plasma cell" = "#779ECB",
  "Basophil" = "#B19CD9",
  "Neutrophil" = "#9E5E9F"
)

genotype_colors <- c(
  "WT" = "#00AFBB",
  "KO" = "#D32F2F"
)


plot_save_dimplot <- function(seurat_obj,
                              group_by,
                              filename,
                              output_dir,
                              colors = NULL,
                              pt_size = 0.5,
                              split_by = NULL,
                              width = 6,
                              height = 6,
                              dpi = 600) {
  p <- scplotter::CellDimPlot(seurat_obj,
    group_by = group_by,
    palcolor = colors,
    pt_size = pt_size,
    facet_by = split_by,
    legend.position = "none",
    show_stat = FALSE,
    bg_color = NA,
    theme = "theme_void"
  ) &
    theme(
      plot.title = element_blank(),
      strip.text = if (!is.null(split_by)) element_blank() else element_text()
    )
  ggsave(file.path(output_dir, filename), p, width = width, height = height, dpi = dpi)
}

plot_save_featureplot <- function(seurat_obj,
                                  features,
                                  filename,
                                  output_dir,
                                  colors,
                                  pt_size = 0.5,
                                  split_by = NULL,
                                  width = 5.5,
                                  height = 5.5,
                                  dpi = 600,
                                  order = TRUE,
                                  max_cutoff = 0.99,
                                  highlight = NULL,
                                  highlight_color = "black",
                                  highlight_size = 1,
                                  highlight_stroke = 0.8,
                                  highlight_alpha = 1) {
  p <- scplotter::FeatureStatPlot(
    seurat_obj,
    features = features,
    plot_type = "dim",
    pt_size = pt_size,
    facet_by = split_by,
    palcolor = colorRampPalette(colors)(100),
    order = if (isTRUE(order)) "high-top" else if (isFALSE(order)) "as-is" else order,
    upper_quantile = max_cutoff,
    bg_cutoff = NULL,
    bg_color = NA,
    theme = "theme_void",
    highlight = highlight,
    highlight_color = highlight_color,
    highlight_size = highlight_size,
    highlight_stroke = highlight_stroke,
    highlight_alpha = highlight_alpha
  ) &
    theme(
      plot.title = element_blank(),
      strip.text = if (!is.null(split_by)) element_blank() else element_text()
    )
  ggsave(file.path(output_dir, filename), p, width = width, height = height, dpi = dpi)
}

scale_marker_matrix <- function(seurat_obj, features, group_by, layer = "data") {
  expr <- Seurat::GetAssayData(seurat_obj, layer = layer)
  cells <- split(colnames(seurat_obj), droplevels(factor(seurat_obj@meta.data[[group_by]])))
  mat <- sapply(cells, function(cl) Matrix::rowMeans(expr[features, cl, drop = FALSE]))
  mat <- matrix(mat, nrow = length(features), dimnames = list(features, names(cells)))
  mat <- t(apply(mat, 1, function(x) { r <- max(x) - min(x); if (r == 0) x * 0 else (x - min(x)) / r }))
  df  <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  setNames(df, c("Feature", group_by, "value"))
}

assign_cols <- function(de_res, logfc_thresh, p_val_adj_thresh) {
  keyvals <- ifelse(de_res$p_val_adj < p_val_adj_thresh &
                    de_res$avg_log2FC > logfc_thresh, "#D32F2F",
                        ifelse(de_res$p_val_adj < p_val_adj_thresh &
                        de_res$avg_log2FC < -logfc_thresh, "#00AFBB", "darkgrey"))

  names(keyvals)[keyvals == "#D32F2F"] <- "Upregulated"
  names(keyvals)[keyvals == "#00AFBB"] <- "Downregulated"
  names(keyvals)[keyvals == "darkgrey"] <- "Not significant"

  return(keyvals)
}

plot_save_volcano <- function(de_res,
                              filename,
                              output_dir,
                              xlim,
                              col_custom,
                              select_lab,
                              ylim = c(0, NA),
                              p_cutoff = 0.05,
                              fc_cutoff = 0.25,
                              flip = TRUE,
                              pseudo_log = TRUE,
                              width = 10,
                              height = 8,
                              dpi = 600) {
  p <- EnhancedVolcano(de_res,
                       lab = de_res$primerid,
                       x = "avg_log2FC",
                       y = "p_val_adj",
                       pCutoff = p_cutoff,
                       FCcutoff = fc_cutoff,
                       xlim = xlim,
                       ylim = ylim,
                       colAlpha = 1,
                       pointSize = 3.0,
                       colCustom = col_custom,
                       selectLab = select_lab,
                       boxedLabels = TRUE,
                       drawConnectors = TRUE,
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE)

  if (pseudo_log) {
    p <- p + scale_y_continuous(transform = scales::pseudo_log_trans(sigma = 5, base = 10))
  }
  if (flip) {
    p <- p + coord_flip()
  }

  ggsave(file.path(output_dir, filename), p, width = width, height = height, dpi = dpi)
  return(p)
}

plot_gsea_res <- function(gsea_res, top_n = 15) {
  gsea_res$result |>
  as_tibble() |>
  slice_min(p.adjust, n = top_n, with_ties = FALSE) |>
  ggplot(aes(x = fct_reorder(Description, NES),
             y = NES,
             fill = p.adjust)) +
      geom_col() +
      scale_fill_gradient(low = "red", high = "blue") +
      coord_flip() +
      theme_light() +
      theme(panel.border = element_rect(color = "black",
                            fill = NA,
                            linewidth = 0.5),
                            axis.title.y = element_blank(),
                            panel.grid = element_blank())
}

plot_geneexpr_pt <- function(
  cds,
  gene,
  span = 1
) {
  if (!gene %in% rownames(cds))
    stop("'", gene, "' not in rownames(cds).")

  pt <- pseudotime(cds)
  keep <- is.finite(pt)
  if (any(!keep))
    message("Dropping ", sum(!keep), " cells with non-finite pseudotime.")

  expr <- normalized_counts(cds, norm_method = "size_only")[gene, keep]

  dt <- data.table(
    pseudotime = pt[keep],
    expression = log1p(expr),
    genotype = droplevels(factor(cds$genotype[keep]))
  )

  message("Plotting target gene log1p(expression) over pseudotime")
  ggplot(dt, aes(x = pseudotime, y = expression, color = genotype)) +
    geom_point(alpha = 0.2) +
    geom_smooth(
      method = "loess",
      span = span,
      fill = "grey70",
      key_glyph = "smooth"
    ) +
    scale_color_manual(values = genotype_colors) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    xlab("Pseudotime") +
    ylab("log1p(normalized counts)") +
    ggtitle(gene)
}