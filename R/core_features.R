#' Get all features from a Seurat object
#'
#' This function returns a vector of all features from a Seurat object, including
#' both gene names and cell names.
#'
#' @param obj A Seurat object
#' @return A vector of all features from the Seurat object
get_all_features <- function(obj){
    all_features <- c()
    assay_names <- names(obj@assays)
    for(i in seq_along(assay_names)){
        assay_key <- obj[[assay_names[i]]]@key
        assay_features <- rownames(obj[[assay_names[i]]])
        names(assay_features) <- paste0(assay_key, assay_features)
        all_features <- c(all_features, assay_features)}
    return(all_features)}

#' remove_features
#'
#' Remove features from current assay to new assay
#' @param obj Seurat object
#' @param features a vector features to remove
#' @param from.assay current assay name, defaults to "RNA"
#' @param to.assay new assay name to store removed features
#' @export
remove_features <- function(obj, features = NULL, from.assay = "RNA", to.assay){

    if(from.assay == "RNA"){
        stopifnot("data" %in% Layers(obj, assay = "RNA"))}

    DefaultAssay(obj) <- from.assay
    keep <- which(rownames(obj) %in% c(features))

    if(length(keep) <= 1){
        message("Detected features <= 1, stop subsetting features")
        return(obj)}
    else{
        message(paste0("Detected features = ", length(keep)))}

    feature.counts <- obj[[from.assay]]$counts[keep,]
    new.counts <- obj[[from.assay]]$counts[-c(keep),]

    if("data" %in% Layers(obj, assay = from.assay)){
        feature.data <- obj[[from.assay]]$data[keep,]
        new.data <- obj[[from.assay]]$data[-c(keep),]
        obj[[to.assay]] <- CreateAssay5Object(counts = as(feature.counts, "sparseMatrix"), data = as(feature.data, "sparseMatrix"), min.features = 0, min.cells = 0)
        obj[[from.assay]] <- CreateAssay5Object(counts = as(new.counts, "sparseMatrix"), data = as(new.data, "sparseMatrix"), min.features = 0, min.cells = 0)}
    else{
        obj[[to.assay]] <- CreateAssay5Object(counts = as(feature.counts, "sparseMatrix"), min.features = 0, min.cells = 0)
        obj[[from.assay]] <- CreateAssay5Object(counts = as(new.counts, "sparseMatrix"), min.features = 0, min.cells = 0)}

    return(obj)}

#' return_features
#'
#' Return features from an assay to a specific assay
#' @param obj Seurat object
#' @param from.assay assay name from which features are taken
#' @param to.assay assay name to which features are stored, defaults to "RNA"
#' @export
return_features <- function(obj, from.assay, to.assay = "RNA"){
    counts <- rbind(obj[[from.assay]]$counts, obj[[to.assay]]$counts)
    data <- rbind(obj[[from.assay]]$data, obj[[to.assay]]$data)
    obj[[to.assay]] <- CreateAssay5Object(counts = counts, data = data, min.cells = 0, min.features = 0)
    obj[[from.assay]] <- NULL
    return(obj)}

#' remove_vdj_features
#'
#' Remove VDJ features from an assay
#' @param obj Seurat object
#' @param bcr if TRUE, remove BCR-VDJ features. Defaults to FALSE
#' @param tcr if TRUE, remove TCR-VDJ features. Defaults to FALSE
#' @param from.assay current assay name, defaults to "RNA"
#' @export
remove_vdj_features <- function(obj, bcr = T, tcr = T, from.assay = "RNA"){
    if(bcr){
        bcr.features <- rownames(obj[[from.assay]])[which(str_detect(rownames(obj[[from.assay]]), bcr.string))]
        obj <- remove_features(obj, features = bcr.features, from.assay = from.assay, to.assay = "BCR")}
    if(tcr){
        tcr.features <- rownames(obj[[from.assay]])[which(str_detect(rownames(obj[[from.assay]]), tcr.string))]
        obj <- remove_features(obj, features = tcr.features, from.assay = from.assay, to.assay = "TCR")}
    return(obj)}

#' intersect_features
#'
#' intersect common features between a list of seurat objects
#' @param obj.list a list of Seurat object
#' @param assay assay name
#' @param min.cells minimum number of cells expressing the gene
#' @export
intersect_features <- function(obj.list, assay = "RNA", min.cells = 3){
    intersect.list <- list()
    for(i in seq_along(obj.list)){
        keep <- rowSums(obj.list[[i]][[assay]]$counts >= 1) >= min.cells
        intersect.list[[i]] <- rownames(obj.list[[i]][[assay]]$counts)[keep]}
    features <- Reduce(intersect, intersect.list)
    return(features)}