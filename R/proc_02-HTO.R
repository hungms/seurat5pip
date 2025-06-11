#' process_hto
#'
#' Run clustering algorithms on HTO counts
#' @param obj Seurat object
#' @param assay name of HTO assay, defaults to "HTO"
#' @export
process_hto <- function(obj, assay = "HTO"){

    assayname <- tolower(assay)
    if(!is.list(obj)){
        obj <- list(obj)}

    for(i in seq_along(obj)){
        nfeatures <- nrow(obj[[i]][[assay]])
        if(nfeatures < 5){
            next}
        features <- rownames(obj[[i]][[assay]])
        DefaultAssay(obj[[i]]) <- assay
        obj[[i]] <- ScaleData(obj[[i]], features = features, verbose = F)
        obj[[i]] <- RunPCA(obj[[i]], features = features, approx = FALSE, reduction.name = paste0("pca_", assayname), verbose = F)
        obj[[i]] <- RunUMAP(obj[[i]], reduction = paste0("pca_", assayname), dims = 1:nfeatures, reduction.name = paste0("umap_", assayname), verbose = F)}

    if(length(obj) == 1){
        obj <- obj[[1]]}

    return(obj)
}



#' run_htodemux
#'
#' HTODemux algorithm to demultiplex cell hashtags
#' @param obj Seurat object
#' @param process if TRUE, run cluster algorithms on HTO counts
#' @export
run_htodemux <- function(obj, assay = "HTO", save_dir){
    
    # Validate the assay
    validate_assay(obj, assay)
    DefaultAssay(obj) <- assay
    
    # Add pseudocount of 1 to each hashtag oligo in each cell
    obj@assays$HTO$counts <- as.matrix(obj@assays$HTO$counts) + 1
    obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2)
    obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)

    # Process HTO counts
    obj <- process_hto(obj, min.features = 5)

    # Subtract the pseudocount after demultiplexing to maintain original values
    obj@assays$HTO$counts <- as.matrix(obj@assays$HTO$counts) - 1

    # Save the demultiplexed object
    save_matrix(obj@assays$HTO$counts, paste0(save_dir, "/counts/"))
    save_matrix(obj@assays$HTO$data, paste0(save_dir, "/data/"))
    return(obj)
}

#' run_multiseqdemux
#'
#' MULTIseqDemux algorithm to demultiplex cell hashtags
#' @param obj Seurat object
#' @param process if TRUE, run cluster algorithms on HTO counts
#' @export
run_multiseqdemux <- function(obj, process = F){
    if("HTO" %in% names(obj@assays)){
        DefaultAssay(obj) <- "HTO"
        obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR", margin = 2)
        obj <- MULTIseqDemux(obj, assay = "HTO")
        obj@meta.data$MULTI.global <- ifelse(obj@meta.data$MULTI_ID == "Doublet", "Doublet", ifelse(obj@meta.data$MULTI_ID == "Negative", "Negative", "Singlet"))
        if(process){
            obj <- process_hto(obj)}
        DefaultAssay(obj) <- "RNA"}
    return(obj)}



    