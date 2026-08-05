suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(ReactomePA)
})

run_pathway_analysis <- function(de_df,
                                 pvalueCutoff = 0.05,
                                 min_gs_size = 10,
                                 max_gs_size = 500,
                                 collection = c("GO", "Reactome"),
                                 analysis = c("GSEA", "ORA"),
                                 org = c("mouse", "human"),
                                 score_col = NULL,
                                 do_simplify = TRUE) {
  collection <- match.arg(collection)
  analysis <- match.arg(analysis)
  org <- match.arg(org)
  orgdb <- if (org == "mouse") org.Mm.eg.db else org.Hs.eg.db
  reactome_org <- if (org == "mouse") "mouse" else "human"

  ranked <- de_df |>
    filter(!is.na(primerid)) |>
    mutate(
      signed_rank_metric = if(!is.null(score_col)) {
        .data[[score_col]]
      } else {
        ifelse(p_val == 0, sign(avg_log2FC) * 300,
               sign(avg_log2FC) * -log10(p_val))
      }) |>
    filter(!is.na(signed_rank_metric))

    mapping <- bitr(ranked$primerid,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = orgdb,
                    drop = TRUE)

    annotated <- ranked |>
      inner_join(mapping, by = c("primerid" = "SYMBOL")) |>
      group_by(ENTREZID) |>
      slice_max(abs(signed_rank_metric), n = 1, with_ties = FALSE) |>
      ungroup()

    gene_list <- annotated |>
      arrange(desc(signed_rank_metric)) |>
      dplyr::select(ENTREZID, signed_rank_metric) |>
      deframe()

    if (analysis == "ORA") {
      genes <- annotated |>
        filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) |>
        pull(ENTREZID)
      universe <- annotated$ENTREZID
      if (collection == "GO") {
        res  <- enrichGO(gene = genes,
                         universe = universe,
                         OrgDb = orgdb,
                         ont = "ALL",
                         pvalueCutoff = pvalueCutoff,
                         minGSSize = min_gs_size,
                         maxGSSize = max_gs_size)
        if (isTRUE(do_simplify) && !is.null(res) && nrow(res) > 0) {
          res <- clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust", select_fun = min)
        }
      } else if (collection == "Reactome") {
        res <- enrichPathway(gene = genes,
                             universe = universe,
                             organism = reactome_org,
                             pvalueCutoff = pvalueCutoff,
                             minGSSize = min_gs_size,
                             maxGSSize = max_gs_size)
      }
    }
    else if (analysis == "GSEA") {
      if (collection == "GO") {
        res <- gseGO(geneList = gene_list,
                     ont = "BP",
                     OrgDb = orgdb,
                     pvalueCutoff = pvalueCutoff,
                     minGSSize = min_gs_size,
                     maxGSSize = max_gs_size,
                     eps = 0
        )

        if (isTRUE(do_simplify) && !is.null(res) && nrow(res) > 0) {
          res <- clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust", select_fun = min)
        }
        
      } else if (collection   == "Reactome") {
        res <- gsePathway(geneList = gene_list,
                          organism = reactome_org,
                          pvalueCutoff = pvalueCutoff,
                          minGSSize = min_gs_size,
                          maxGSSize = max_gs_size,
                          eps = 0
        )
      }
    }
    if (is.null(res) || nrow(res) == 0) {
      warning("No significant pathways found.")
      return(list(
        result = res,
        table = tibble(),
        ranking = gene_list
      ))
    }
    res <- setReadable(res, OrgDb = orgdb, keyType = "ENTREZID")

    list(
      result = res,
      table = as_tibble(res@result),
      ranking = gene_list
    )
}