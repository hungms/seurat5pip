#' Calculate cell cycle phase
#'
#' This function calculates the cell cycle phase.
#' @param obj Seurat object
#' @param org Organism, either "human" or "mouse"
#' @param assay assay name, defaults to "RNA"
#' @return Seurat object
#' @export
cc_by_seurat <- function(obj, split.by = NULL, assay = "RNA"){

    # validate split.by
    split.by <- validate_split.by(split.by, obj)

    # split object
    obj.list <- split_object(obj, split.by)

    # for each object
    for(o in 1:length(obj.list)){

        # get organism
        org <- get_org(obj.list[[o]], assay)

        # get cell cycle genes
        cc_list <- get_cc_genes(org = org)

        # calculate cell cycle phase
        obj.list[[o]] <- CellCycleScoring(obj.list[[o]], s.features = cc_list$Seurat_S, g2m.features = cc_list$Seurat_G2M, set.ident = F, assay = assay)
    }

    # merge objects
    obj <- merge_objects(obj.list)

    obj$Phase <- factor(obj$Phase, levels = c("G1", "S", "G2M"))

    # log
    log_function()

    # return object
    return(obj)}



#'
#' 
cc_by_tricycle <- function(obj, assay = "RNA"){

    # get organism
    org <- get_org(obj, assay)

    # get cell cycle genes
    cc_list <- get_cc_genes(org = org)

    # create SingleCellExperiment object
    

    sce <- tricycle::project_cycle_space(sce, species = 'human', gname.type = gname.type)
  print(scater::plotReducedDim(sce, dimred = "tricycleEmbedding") +
          labs(x = "Projected PC1", y = "Projected PC2") +
          ggtitle(sprintf("Projected cell cycle space (n=%d)", ncol(sce))) +
          theme_bw(base_size = 14)
  )

  sce <- tricycle::estimate_cycle_position(sce)
  print(scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "tricyclePosition") +
          labs(x = "Projected PC1", y = "Projected PC2") +
          ggtitle(sprintf("Projected cell cycle space (n=%d), tricyclePosition", ncol(sce))) +
          theme_bw(base_size = 14)
  )

  sce <- tricycle::estimate_Schwabe_stage(sce, gname.type = gname.type, species = species)
  print(scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "CCStage") +
          labs(x = "Projected PC1", y = "Projected PC2",
          title = paste0("Projected cell cycle space (n=", ncol(sce), "), CCStage")) +
          theme_bw(base_size = 14)
  )

  print(tricycle::plot_emb_circle_scale(sce, dimred = 1,
      point.size = 3.5, point.alpha = 0.9) +
      theme_bw(base_size = 14)
  )

  toAdd <- data.frame(row.names = rownames(sce@colData), tricyclePosition = sce@colData$tricyclePosition, CCStage = sce@colData$CCStage)
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = toAdd)

  print(Seurat::FeaturePlot(seuratObj, features = 'tricyclePosition'))
  print(Seurat::DimPlot(seuratObj, group.by = 'CCStage'))

  return(seuratObj)
}