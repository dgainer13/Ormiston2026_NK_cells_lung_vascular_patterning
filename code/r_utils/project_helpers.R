make_paths <- function(script_id, data = list()) {
  out <- list(
    figures = here::here("results", "figures", script_id),
    tables  = here::here("results", "tables",  script_id),
    outputs = here::here("results", "outputs", script_id)
  )
  purrr::walk(out, dir.create, recursive = TRUE, showWarnings = FALSE)
  list(data = data, out = out)
}

kegg_symbols <- function(pathway_id) {
  genes <- KEGGREST::keggGet(pathway_id)[[1]]$GENE
  genes <- genes[seq(2, length(genes), by = 2)]
  gsub("\\;.*", "", genes)
}

go_symbols <- function(goid) {
  unique(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = goid,
                               keytype = "GOALL", columns = "SYMBOL")$SYMBOL)
}