build_cellchat <- function(obj_list, db, min_cells, group_by = "celltype", assay = NULL) {
  lapply(obj_list, function(x) {
    cc <- createCellChat(object = x, group.by = group_by, assay = assay)
    cc@DB <- db

    cc |>
      subsetData() |>
      identifyOverExpressedGenes() |>
      identifyOverExpressedInteractions() |>
      computeCommunProb(
        type = "triMean",
        trim = NULL,
        raw.use = TRUE,
        population.size = TRUE
      ) |>
      filterCommunication(min.cells = min_cells) |>
      computeCommunProbPathway() |>
      aggregateNet() |>
      netAnalysis_computeCentrality()
  })
}

cellchat_diff_genes <- function(cc, pos_dataset, neg_dataset, thresh.fc, ligand.logFC = 0.1) {
  cc <- identifyOverExpressedGenes(cc,
    group.dataset = "datasets",
    pos.dataset = pos_dataset,
    features.name = pos_dataset,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = thresh.fc
  )

  net <- netMappingDEG(cc, features.name = pos_dataset)
  net_up <- subsetCommunication(cc,
    net = net,
    datasets = pos_dataset,
    ligand.logFC = ligand.logFC,
    receptor.logFC = NULL
  )
  net_down <- subsetCommunication(cc,
    net = net,
    datasets = neg_dataset,
    ligand.logFC = -ligand.logFC,
    receptor.logFC = NULL
  )

  list(cc = cc, net_up = net_up, net_down = net_down)
}

cellchat_chord <- function(cc, sources.use, targets.use, net, color.use, show.legend) {
netVisual_chord_gene(cc,
    sources.use = sources.use,
    targets.use = targets.use,
    slot.name = "net",
    net = net,
    lab.cex = 0.8,
    small.gap = 3.5,
    color.use = color.use,
    show.legend = show.legend
  )
}
