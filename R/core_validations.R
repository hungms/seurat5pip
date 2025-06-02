#' Validate assay
#'
#' This function validates that the assay is present in the Seurat object.
#' @param obj Seurat object
#' @param assay Assay to validate
#' @return Assay
#' @export
validate_assay <- function(obj, assay){
    if(!assay %in% names(obj@assays)){
        stop("Assay not found in Seurat object")}
    else{
        if(class(obj[[assay]]) == "Assay"){
            obj[[assay]] <- as(obj[[assay]], "Assay5")}
        stopifnot(class(obj[[assay]]) %in% c("Assay5", "SCTAssay"))}
    return(obj)}


#' Validate Seurat object
#'
#' This function validates that the Seurat object is valid.
#' @param obj Seurat object
#' @return Seurat object
#' @export
validate_obj <- function(obj, assay = "RNA", min.features = 200, min.pct = 0.01, min.cells = NULL){

    # validate min.pct and min.cells
    stopifnot(is.null(min.cells) != is.null(min.pct))
    stopifnot(min.pct > 0)

    # validate assay
    obj <- validate_assay(obj, assay)

    # get old dimensions
    old_dim <- dim(obj[[assay]]$counts)

    # if min.pct
    if(!is.null(min.pct)){
        
        # get features with min.pct expression
        features.to.keep <- rownames(obj[[assay]]$counts)[rowSums(obj[[assay]]$counts > 0) / ncol(obj[[assay]]$counts) >= min.pct]

        # create new assay
        obj[[assay]] <- CreateAssay5Object(
            counts = obj[[assay]]$counts[features.to.keep, ], 
            data = obj[[assay]]$data[features.to.keep, ], 
            min.features = min.features, 
            min.cells = 0)
        
        # min.cells.message
        min.cells.message <- paste0(min.pct, "%")
        }

    # if min.cells
    else if(!is.null(min.cells)){

        # create new assay
        obj[[assay]] <- CreateAssay5Object(
            counts = obj[[assay]]$counts, 
            data = obj[[assay]]$data, 
            min.features = min.features, 
            min.cells = min.cells)
        
        # min.cells.message
        min.cells.message <- min.cells
        } 

    # get new dimensions
    new_dim <- dim(obj[[assay]]$counts)

    # print message
    m1 <- paste0(old_dim[1] - new_dim[1], " features removed from ",  assay, " assay")
    m2 <- paste0(old_dim[2] - new_dim[2], " cells removed from ",  assay, " assay")
    log_function(m1, m2)

    # return object
    return(obj)}