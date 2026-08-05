suppressPackageStartupMessages({
  library(monocle3)
  library(tidyverse)
})

run_pseudotime_diff_exp <- function(obj,
                                    min_pct = 0.05,
                                    cores = 2,
                                    model_formula_str =
                                      "~ genotype * pseudotime",
                                    contrast = "genotypeKO:pseudotime",
                                    q_value_cutoff = 0.05,
                                    p_adjust_method = "BH",
                                    expression_family = c(
                                      "quasipoisson",
                                      "negbinomial",
                                      "poisson",
                                      "binomial",
                                      "gaussian",
                                      "zipoisson",
                                      "zinegbinomial",
                                      "mixed-negbinomial"
                                    )) {
  expression_family <- match.arg(expression_family)

  if ("genotype" %in% colnames(colData(obj)) &&
        !is.factor(colData(obj)$genotype))
    warning("'genotype' is not a factor; its reference level will be chosen ",
            "alphabetically. Set e.g. relevel(factor(genotype), ref = 'WT').")

  finite_pt <- is.finite(pseudotime(obj))
  if (any(!finite_pt)) {
    message("Dropping ", sum(!finite_pt), " cells with non-finite pseudotime.")
    obj <- obj[, finite_pt]
  }

  message("Keeping genes expressed in > ", min_pct * 100, "% of cells")
  keep <- Matrix::rowSums(counts(obj) > 0) / ncol(obj) > min_pct
  obj_filt <- obj[keep, ]
  obj_genes <- rownames(obj)[keep]
  message(sum(keep), "/", length(keep), " genes retained.")

  message("Fitting models (", expression_family, "): ", model_formula_str)
  de_res <- fit_models(obj_filt, model_formula_str = model_formula_str,
                       expression_family = expression_family,
                       verbose = TRUE, cores = cores)

  de_res_table <- coefficient_table(de_res) |>
    dplyr::filter(status == "OK") |>
    dplyr::group_by(term) |>
    dplyr::mutate(q_value = p.adjust(p_value, method = p_adjust_method)) |>
    dplyr::ungroup()

  avail <- sort(unique(de_res_table$term))
  if (!contrast %in% avail)
    stop("contrast '", contrast, "' not in fitted terms: ",
         paste(avail, collapse = ", "))

  de_res_table_sig <- de_res_table |>
    dplyr::filter(term == contrast) |>
    dplyr::filter(q_value < q_value_cutoff) |>
    dplyr::select(term, gene_id, estimate, std_err, p_value, q_value) |>
    dplyr::arrange(q_value)

  list(filtered_cds = obj_filt,
       filtered_genes = obj_genes,
       de_res_allgenes = de_res_table,
       de_res_sig = de_res_table_sig)
}
