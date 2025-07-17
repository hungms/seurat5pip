#' hvf_by_seurat_vst
#' 
#' @param obj Seurat object
#' @param normalization.method Normalization method
#' @param split.by Split by
#' @param nfeatures Number of features
#' @param output_dir Output directory
#' @return Seurat object
#' @export
hvf_by_seurat_vst <- function(obj, normalization.method = "SCT", split.by = NULL, nfeatures = 2000, output_dir = NULL){

    stopifnot(normalization.method %in% c("LogNormalize", "SCT"))

    split.by <- validate_split.by(split.by, obj)
    obj.list <- SplitObject(obj, split.by = split.by)

    for(i in seq_along(obj.list)){
        if(normalization.method == "SCT"){
            obj.list[[i]] <- SCTransform(obj.list[[i]], verbose = FALSE)}
        obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], nfeatures = nfeatures)}

    features <- SelectIntegrationFeatures(obj.list, nfeatures = nfeatures)

    if(normalization.method == "SCT"){
        obj <- SCTransform(obj, return.only.var.genes = F, verbose = FALSE)}
    else {
	obj <- FindVariableFeatures(obj, nfeatures = nfeatures)}

    VariableFeatures(obj) <- features

    p1 <- VariableFeaturePlot(obj) +
        ggprism::theme_prism(border = TRUE) +
        geom_text_repel(
            data = VariableFeaturePlot(obj)$data %>%
                rownames_to_column("gene") %>%
                filter(gene %in% features) %>%
                as.data.frame(),
                aes(label = gene),
                size = 3)

    dir.create(paste0(output_dir, "/seurat-vst/"))
    write.table(features, paste0(output_dir, "/seurat-vst/", normalization.method, "_", nfeatures, ".txt"), row.names = F)
    write_png(p1, output_dir = paste0(output_dir, "/seurat-vst/"), filename = paste0(normalization.method, "_", nfeatures, ".png"), width = 2000, height = 1600)

    return(obj)
}

#' find_npcs
#' 
#' @param obj Seurat object
#' @param reduction Reduction name
#' @return Number of PCs
#' @export
find_npcs <- function(obj, reduction = "pca"){
    pca_variance <- cumsum(obj[[reduction]]@stdev ^ 2) / sum(obj[[reduction]]@stdev ^ 2)
    npcs <- min(which(pca_variance > 0.95))
    return(npcs)}


#' run_pca
#' 
#' @param obj Seurat object
#' @param split.by Split by
#' @param reduction Reduction name
#' @param output_dir Output directory
#' @return Seurat object
#' @export
run_pca <- function(obj, split.by = NULL, pca.name = "pca", output_dir = NULL){

    stopifnot(length(VariableFeatures(obj)) > 0)

    split.by <- validate_split.by(split.by, obj)
    
    obj <- RunPCA(obj, features = VariableFeatures(obj), ndims = 50, reduction.name = reduction)

    # save embeddings   
    pca_embeddings <- obj[[pca.name]]@cell.embeddings
    dir.create(paste0(output_dir, "/", pca.name, "/"), showWarnings = FALSE)
    write.csv(pca_embeddings, file = paste0(output_dir, "/", pca.name, "/embeddings.csv"))

    # save loadings
    pca_loadings <- obj[[pca.name]]@feature.loadings
    write.csv(pca_loadings, file = paste0(output_dir, "/", pca.name, "/loadings.csv"))

    # save variance
    pca_variance <- data.frame(row.names = paste0("PC", 1:50), obj[[pca.name]]@stdev^2 / sum(obj[[pca.name]]@stdev^2))
    write.csv(pca_variance, file = paste0(output_dir, "/", pca.name, "/variance.csv"))

    # calculate cumulative variance
    npcs <- find_npcs(obj, pca.name)

    # plot elbow, feature loadings, and PCA
    p1 <- ElbowPlot(obj, ndims = 50, reduction = pca.name) + geom_vline(xintercept = npcs, linetype = "dashed", color = "red")
    p2 <- plot_featureloadings(obj, reduction = pca.name) + geom_hline(yintercept = npcs, linetype = "dashed", color = "grey30")
    p3 <- plot_reduction(obj, group.by = split.by, reduction = pca.name, reduction.type = "PCA", label = F, mask = F, count.groups = F, shuffle = T)

    # save plots
    write_png(p1, output_dir = paste0(output_dir, "/", pca.name, "/"), filename = "plot_elbow.png", width = 1200, height = 1000)
    write_png(p2, output_dir = paste0(output_dir, "/", pca.name, "/"), filename = "plot_feature_loadings.png", width = 4000, height = 6000)
    write_png(p3, output_dir = paste0(output_dir, "/", pca.name, "/"), filename = "plot_pca.png", width = 2000, height = 1600)

    return(obj)
}


#' run_umap
#' 
#' @param obj Seurat object
#' @param split.by Split by
#' @param reduction Reduction name
#' @param output_dir Output directory
#' @return Seurat object
#' @export
run_umap <- function(obj, split.by = NULL, reduction = "pca", umap.name = "umap", npcs = NULL, output_dir = NULL){

    split.by <- validate_split.by(split.by, obj)

    if (is.null(npcs)) {
        npcs <- find_npcs(obj, reduction)
        npcs_label <- paste0("95% variation : first ", npcs, " PCs")
    } else {
        npcs_label <- paste0("User-defined : first ", npcs, " PCs")
    }
    obj <- RunUMAP(obj, reduction = reduction, dims = 1:npcs, reduction.name = umap.name)
    write.csv(obj[[umap.name]]@cell.embeddings[1:2], paste0(output_dir, "/", umap.name, "/embeddings.csv"), row.names = T)
    plot <- plot_reduction(obj, group.by = split.by, reduction = umap.name, label = F, mask = F, count.groups = F, shuffle = T) +
        labs(title = npcs_label)
    if (!is.null(output_dir)) {
        write_png(plot, output_dir = paste0(output_dir, "/", umap.name, "/"), filename = paste0("plot_", umap.name, "_", npcs, ".png"), width = 2000, height = 1600)
    }
    return(obj)
}

