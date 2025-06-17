#' Calculate cell cycle phase with Seurat
#'
#' This function calculates the cell cycle phase.
#' @param obj Seurat object
#' @param org Organism, either "human" or "mouse"
#' @param assay assay name, defaults to "RNA"
#' @return Seurat object
#' @export
cc_by_seurat <- function(obj, split.by = NULL, assay = "RNA", output_dir = NULL){

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

    # plot cell cycle
    plot_cc(obj, var = "Phase", split.by = split.by, output_dir = output_dir)

    # log
    log_function()

    # return object
    return(obj)}



#' Calculate cell cycle phase with tricycle
#' 
cc_by_tricycle <- function(obj, assay = "RNA", output_dir = NULL){

    # get organism
    org <- get_org(obj, assay)

    # scale data
    if(!"data" %in% Layers(obj, assay = assay)){
        obj <- NormalizeData(obj, assay = assay)
    }
    if(!"scale.data" %in% Layers(obj, assay = assay)){
        obj <- ScaleData(obj, assay = assay)
    }

    # convert to SingleCellExperiment object
    sce <- as.SingleCellExperiment(obj, assay = assay)
    print(sce)

    # create SingleCellExperiment object
    sce <- tricycle::project_cycle_space(sce, species = org, gname.type = "SYMBOL")
    sce <- tricycle::estimate_cycle_position(sce)
    p1 <- scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "tricyclePosition") +
        labs(x = "PC1", y = "PC2", title = "tricyclePosition") +
        ggprism::theme_prism(border = T) +
        coord_fixed()

    #top2a.idx <- which(rowData(sce)$Gene %in% c('Top2a', 'TOP2A'))
    #fit.l <- tricycle::fit_periodic_loess(
    #    sce$tricyclePosition,
    #    assay(sce, 'logcounts')[top2a.idx,],
    #    plot = TRUE,
   #     x_lab = "Cell cycle position", y_lab = "log2(Top2a)")
   # p2 <- fit.l$fig + ggprism::theme_prism(border = T)

    sce <- tricycle::estimate_Schwabe_stage(sce, gname.type = "SYMBOL", species = org)
    p3 <- scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "CCStage") +
        labs(x = "PC1", y = "PC2", title = "CCStage") +
        ggprism::theme_prism(border = T) +
        coord_fixed()

    p4 <- tricycle::plot_emb_circle_scale(sce, dimred = 1,
        point.size = 3.5, point.alpha = 0.9) +
        ggprism::theme_prism(border = T) + 
        coord_fixed()
    legend <- tricycle::circle_scale_legend(text.size = 5, alpha = 0.9)
    p4 <- plot_grid(p4, legend, ncol = 2, rel_widths = c(1, 0.4))

    # save plots
    write_png(p1, output_dir = output_dir, filename = "tricycle_position.png", width = 1500, height = 1200)
    #write_png(p2, output_dir = output_dir, filename = "tricycle_TOP2A.png", width = 1200, height = 1000)
    write_png(p3, output_dir = output_dir, filename = "tricycle_CCstage.png", width = 1500, height = 1200)
    write_png(p4, output_dir = output_dir, filename = "tricycle_circle.png", width = 1800, height = 1500)


    toAdd <- data.frame(row.names = rownames(sce@colData), tricyclePosition = sce@colData$tricyclePosition, CCStage = sce@colData$CCStage)
    write.csv(toAdd, file.path(output_dir, "seuratobj_metadata_tricycle.csv"))

    # add metadata
    if(any(colnames(obj@meta.data) %in% colnames(toAdd))){
        obj@meta.data[which(colnames(obj@meta.data) %in% colnames(toAdd))] <- NULL}
    obj <- Seurat::AddMetaData(obj, metadata = toAdd)

    log_function()

    return(obj)
}