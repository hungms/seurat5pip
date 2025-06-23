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
#' This function calculates the cell cycle phase with tricycle.
#' @param obj Seurat object
#' @param org Organism, either "human" or "mouse"
#' @param assay assay name, defaults to "RNA"
#' @param output_dir output directory
#' @return Seurat object
#' @export
cc_by_tricycle <- function(obj, split.by = NULL, assay = "RNA", output_dir = NULL){

    # get organism
    org <- get_org(obj, assay)
    split.by <- validate_split.by(split.by, obj)

    obj.list <- SplitObject(obj, split.by = split.by)

    toAddlist <- list()

    for(i in seq_along(obj.list)){
        tmp <- obj.list[[i]]

        # scale data
        stopifnot(all(c("data", "scale.data") %in% Layers(tmp, assay = assay)))

        # convert to SingleCellExperiment object
        sce <- as.SingleCellExperiment(tmp, assay = assay)
        print(sce)

        # run tricycle
        sce <- tricycle::project_cycle_space(sce, species = org, gname.type = "SYMBOL")
        sce <- tricycle::estimate_cycle_position(sce)
        sce <- tricycle::estimate_Schwabe_stage(sce, species = org, gname.type = 'SYMBOL')

        p1 <- scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "tricyclePosition") +
            labs(x = "PC1", y = "PC2", title = "tricyclePosition") +
            ggprism::theme_prism(border = T) +
            coord_fixed()

        p2 <- scater::plotReducedDim(sce, dimred = "tricycleEmbedding", colour_by = "CCStage") +
            labs(x = "PC1", y = "PC2", title = "CCStage") +
            ggprism::theme_prism(border = T) +
            coord_fixed()

        p3 <- tricycle::plot_emb_circle_scale(sce, dimred = 1,
            point.size = 3.5, point.alpha = 0.9) +
            ggprism::theme_prism(border = T) + 
            coord_fixed()
        legend <- tricycle::circle_scale_legend(text.size = 5, alpha = 0.9)
        p3 <- plot_grid(p3, legend, ncol = 2, rel_widths = c(1, 0.4))


        # save plots
        dir.create(paste0(output_dir, "/", names(obj.list)[i], "/"), showWarnings = F, recursive = T)
        write_png(p1, output_dir = paste0(output_dir, "/", names(obj.list)[i], "/"), filename = "tricycle_position.png", width = 1500, height = 1200)
        write_png(p2, output_dir = paste0(output_dir, "/", names(obj.list)[i], "/"), filename = "tricycle_CCstage.png", width = 1500, height = 1200)
        write_png(p3, output_dir = paste0(output_dir, "/", names(obj.list)[i], "/"), filename = "tricycle_circle.png", width = 1800, height = 1500)

        toAddlist[[i]] <- data.frame(row.names = rownames(sce@colData), tricyclePosition = sce@colData$tricyclePosition, CCStage = sce@colData$CCStage)}

    toAdd <- bind_rows(toAddlist)
    write.csv(toAdd, file.path(output_dir, "seuratobj_metadata_tricycle.csv"))

    # add metadata
    if(any(colnames(obj@meta.data) %in% colnames(toAdd))){
        obj@meta.data[which(colnames(obj@meta.data) %in% colnames(toAdd))] <- NULL}
    obj <- Seurat::AddMetaData(obj, metadata = toAdd)

    log_function()

    return(obj)
}