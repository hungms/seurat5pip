#' Validate assay
#'
#' This function validates that the assay is present in the Seurat object.
#' @param obj Seurat object
#' @param assay Assay to validate
#' @return Assay
#' @export
validate_assays <- function(obj){

    if(!all(class(obj@assays) %in% c("Assay", "Assay5", "SCTAssay"))){
        stop("Assays must be of class Assay, Assay5, or SCTAssay")}

    # for each assay
    assay.name <- names(obj@assays)[which(class(obj@assays) == "Assay")]

    # for all Assay3
    for(assay in which.assay){
        # convert to Assay5
        obj[[assay]] <- as(obj[[assay]], "Assay5")}

    # log
    m1 <- paste0("Converted ", paste0(assay.name, collapse = ", "), " from Assay3 to Assay5")
    log_function(m1)

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