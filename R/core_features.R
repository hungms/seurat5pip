#' Get all features from a Seurat object
#'
#' This function returns a vector of all features from a Seurat object, including
#' both gene names and cell names.
#'
#' @param obj A Seurat object
#' @return A vector of all features from the Seurat object, with assay keys included in vector names
get_all_features <- function(obj){
    all_features <- c()
    assay_names <- names(obj@assays)
    for(i in seq_along(assay_names)){
        assay_key <- obj[[assay_names[i]]]@key
        assay_features <- rownames(obj[[assay_names[i]]])
        names(assay_features) <- paste0(assay_key, assay_features)
        all_features <- c(all_features, assay_features)}
    return(all_features)}

#' subset_features
#'
#' Subset features from current assay to new assay
#' @param obj Seurat object
#' @param features a vector features to remove
#' @param assay.from current assay name
#' @param assay.to new assay name to store removed features
#' @export
subset_features <- function(obj, features = NULL, assay.from, assay.to){

    # validate data layer
    if(class(obj[[assay.from]]) == "Assay5" & !"data" %in% Layers(obj, assay = assay.from)){
        stop("Please normalize your data before removing features")}

    # set default assay
    DefaultAssay(obj) <- assay.from

    # find index of features to keep
    index.to.keep <- which(rownames(obj) %in% c(features))
    if(length(index.to.keep) <= 1){
        warning("Number of features detected <= 1, stop subsetting")
        return(obj)}

    # subset counts
    assay.to.counts <- obj[[assay.from]]$counts[index.to.keep,]
    assay.to.data <- obj[[assay.from]]$data[index.to.keep,]
    assay.from.counts <- obj[[assay.from]]$counts[-c(index.to.keep),]
    assay.from.data <- obj[[assay.from]]$data[-c(index.to.keep),]

    # store assays
    obj[[assay.to]] <- CreateAssay5Object(counts = as(assay.to.counts, "sparseMatrix"), data = as(assay.to.data, "sparseMatrix"), min.features = 0, min.cells = 0)
    obj[[assay.from]] <- CreateAssay5Object(counts = as(assay.from.counts, "sparseMatrix"), data = as(assay.from.data, "sparseMatrix"), min.features = 0, min.cells = 0)
    
    # log
    m1 <- paste0("number of features to subset = ", length(index.to.keep))
    log_function(m1)

    return(obj)}

#' subset_vdj_features
#'
#' subset VDJ features from an assay
#' @param obj Seurat object
#' @param bcr if TRUE, remove BCR-VDJ features. Defaults to FALSE
#' @param tcr if TRUE, remove TCR-VDJ features. Defaults to FALSE
#' @param from.assay current assay name, defaults to "RNA"
#' @export
subset_vdj_features <- function(obj, assay, bcr = T, tcr = T){

    # if subset BCR features
    if(bcr){
        message("LOG : subset_vdj_features for BCR")
        bcr.features <- rownames(obj[[assay]])[which(str_detect(rownames(obj[[assay]]), bcr.string))]
        obj <- subset_features(obj, features = bcr.features, assay.from = assay, assay.to = "BCR")}

    # if subset TCR features
    if(tcr){
        message("LOG : subset_vdj_features for TCR")
        tcr.features <- rownames(obj[[assay]])[which(str_detect(rownames(obj[[assay]]), tcr.string))]
        obj <- subset_features(obj, features = tcr.features, assay.from = assay, assay.to = "TCR")}
    
    return(obj)}


#' remove_features
#'
#' remove features from a Seurat object
#' @param obj Seurat object
#' @param features a vector of features to remove
#' @export
remove_features <- function(obj, assay, features){

    # 
    stopifnot(class(obj[[assay]]) == "Assay5")

    # remove features
    obj <- obj[-index.to.remove]

    return(obj)}


#' intersect_features
#'
#' intersect common features between a list of seurat objects
#' @param obj.list a list of Seurat object
#' @param assay assay name
#' @param min.cells minimum number of cells expressing the gene
#' @export
intersect_features <- function(obj, split.by, assay){

    # split object
    obj.list <- split_object(obj, split.by)

    # initialize list
    intersect.list <- list()

    # get features from each object
    for(i in seq_along(obj.list)){
        intersect.list[[i]] <- rownames(obj.list[[i]][[assay]]$counts)}

    # intersect features across all objects
    features <- Reduce(intersect, intersect.list)
    return(features)}