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
    write_png(p1, output_dir = paste0(output_dir, "/seurat-vst/"), filename = paste0(normalization.method, "_", nfeatures, ".png"), width = 4000, height = 3500)

    return(obj)
}

#' run_pca
#' 
#' @param obj Seurat object
#' @param reduction Reduction name
#' @param output_dir Output directory
#' @return Seurat object
#' @export
run_pca <- function(obj, reduction = "pca", output_dir = NULL){

    stopifnot(length(VariableFeatures(obj)) > 0)
    
    obj <- RunPCA(obj, features = VariableFeatures(obj), ndims = 50, reduction.name = reduction)

    # save embeddings
    pca_embeddings <- obj[[reduction]]@cell.embeddings
    dir.create(paste0(output_dir, "/", reduction, "/"), showWarnings = FALSE)
    write.csv(pca_embeddings, file = paste0(output_dir, "/", reduction, "/embeddings.csv"))

    # save loadings
    pca_loadings <- obj[[reduction]]@feature.loadings
    write.csv(pca_loadings, file = paste0(output_dir, "/", reduction, "/loadings.csv"))

    # save variance
    pca_variance <- data.frame(row.names = paste0("PC", 1:50), obj[[reduction]]@stdev^2 / sum(obj[[reduction]]@stdev^2))
    write.csv(pca_variance, file = paste0(output_dir, "/", reduction, "/variance.csv"))

    # plot elbow, feature loadings, and PCA
    p1 <- ElbowPlot(obj, ndims = 50, reduction = reduction)
    p2 <- plot_featureloadings(obj, reduction = reduction, output_dir = output_dir)
    p3 <- plot_reduction(obj, reduction = reduction, reduction.type = "PCA", output_dir = output_dir)

    # save plots
    write_png(p1, output_dir = output_dir, filename = "elbow.png", width = 1000, height = 800)
    write_png(p2, output_dir = output_dir, filename = "feature_loadings.png", width = 1000, height = 800)
    write_png(p3, output_dir = output_dir, filename = "pca.png", width = 1000, height = 800)

    return(obj)
}