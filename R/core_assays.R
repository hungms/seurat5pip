#' set_assay_keys
#'
#' Bugfix for setting assay keys when creating Seurat object
#' @param obj Seurat object
#' @export
set_assay_keys <- function(obj){
    for(assay in names(obj@assays)){
        obj[[assay]]@key <- paste0(tolower(assay), "_")}
    return(obj)}

#' get_assay_keys
#'
#' Get assay keys from a Seurat object
#' @param obj Seurat object
#' @return A vector of assay keys
#' @export
get_assay_keys <- function(obj){
    return(sapply(names(obj@assays), function(assay) obj[[assay]]@key))}

#' join_layers
#'
#' Join layers for all Assay5 objects in a Seurat object
#' @param obj Seurat object
#' @return Seurat object
#' @export
join_layers <- function(obj){
    
    # which assays are Assay5 objects?
    classes <- unlist(lapply(obj@assays, class))
    assays.to.join <- names(classes)[which(classes == "Assay5")]
    
    for(i in assays.to.join){
        obj[[i]] <- JoinLayers(obj[[i]])}
    return(obj)}



