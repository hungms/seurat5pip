#' Validate assay
#'
#' This function validates that the assay is present in the Seurat object.
#' @param obj Seurat object
#' @param assay Assay to validate
#' @return Assay
#' @export
validate_assays <- function(obj){

    # get class of assays
    assay.class <- unlist(lapply(obj@assays, class))

    # validate assays
    if(!all(assay.class %in% c("Assay", "Assay5", "SCTAssay"))){
        stop("Assays must be of class Assay, Assay5, or SCTAssay")}

    # for each assay
    assay.name <- names(assay.class)[which(assay.class == "Assay")]

    # for all Assay3
    for(assay in assay.name){
        # convert to Assay5
        obj[[assay]] <- as(obj[[assay]], "Assay5")}

    # log
    paste0("Converted ", paste0(names(assay.class), collapse = ", "), " from Assay3 to Assay5")
    return(obj)}

#' Validate object
#'
#' This function validates that the object is a Seurat object.
#' @param obj Seurat object
#' @return Seurat object
#' @export
validate_object <- function(obj){

    # validate assays
    obj <- validate_assays(obj)

    # project must be in metadata
    if(!"project" %in% names(obj@meta.data)){
        stop("Please provide a metadata column named 'project' specifying the sequencing runs")}
    else{
        # make sure project is a character vector
        obj$project <- as.character(obj$project)}

    return(obj)}



#' Validate split.by
#'
#' This function validates that the split.by is a character vector.
#' @param split.by Split.by to validate
#' @return Split.by
#' @export
validate_split.by <- function(split.by, obj){
    if(is.null(split.by)){
        split.by <- "project"
        validate_object(obj)
        return(split.by)}
    else{
        if(!split.by %in% colnames(obj@meta.data)){
            stop("split.by must be a column in the metadata")}
        return(split.by)}
}