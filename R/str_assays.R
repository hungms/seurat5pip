#' join_layers
#'
#' Joins all layers for each Assay5 object in a Seurat object.
#' This function identifies all assays of class Assay5 and applies JoinLayers to each one.
#' 
#' @param obj A Seurat object containing one or more Assay5 objects
#' @return The Seurat object with all Assay5 layers joined
#' @export
join_layers <- function(obj){
    
    # find all Assay5 class
    classes <- unlist(lapply(obj@assays, class))
    assays.to.join <- names(classes)[which(classes == "Assay5")]
    
    # join every Assay5 class
    for(i in assays.to.join){
        obj[[i]] <- JoinLayers(obj[[i]])}
    return(obj)}

#' set_assay_keys
#'
#' Sets standardized keys for all assays in a Seurat object.
#' This function addresses inconsistencies in assay key naming by
#' setting each assay key to the lowercase assay name followed by an underscore.
#'
#' @param obj A Seurat object
#' @return The Seurat object with standardized assay keys
#' @keywords internal
set_assay_keys <- function(obj){
    for(assay in names(obj@assays)){
        obj[[assay]]@key <- paste0(tolower(assay), "_")}
    return(obj)}

#' get_assay_keys
#'
#' Retrieves all assay keys from a Seurat object.
#' This function extracts the key attribute from each assay in the Seurat object.
#'
#' @param obj A Seurat object
#' @return A named vector of assay keys
#' @keywords internal
get_assay_keys <- function(obj){
    return(sapply(names(obj@assays), function(assay) obj[[assay]]@key))}

#' merge_assays
#'
#' Merge two Assay5 objects in a Seurat object, combining their features
#' @param obj Seurat object containing both assays
#' @param assay1 Name of the assay that will contain the merged data
#' @param assay2 Name of the assay that will be merged into assay1 and then removed
#' @return Seurat object with merged assays
#' @export
merge_assays <- function(obj, assay1, assay2){

    # validate assay5
    if(all(c(class(obj[[assay1]]), class(obj[[assay2]])) != "Assay5")){
        stop("assay1 and assay2 must be Assay5 objects")}

    # validate data layer
    if(all(Layers(obj, assay = assay2) %in% Layers(obj, assay = assay1))){
        stop("assay2 must contain the same layers as assay1")}

    # get counts and data
    counts <- rbind(obj[[assay1]]$counts, obj[[assay2]]$counts)

    # if data layer is not present in assay1 and assay2
    if(!"data" %in% Layers(obj, assay = assay1) | !"data" %in% Layers(obj, assay = assay2)){

        # store counts only
        obj[[assay1]] <- CreateAssay5Object(counts = counts, min.cells = 0, min.features = 0)}

    # if data layer is present in both assay1 and assay2
    else{

        # store counts and data
        data <- rbind(obj[[assay1]]$data, obj[[assay2]]$data)
        obj[[assay1]] <- CreateAssay5Object(counts = counts, data = data, min.cells = 0, min.features = 0)}

    # log message
    m1 <- paste("layers returned = ", paste(Layers(obj, assay = assay2), collapse = ", "))
    m2 <- paste("number of features returned = ", nrow(obj[[assay2]]$counts))
    log_function(m1, m2)

    # remove assay2
    obj[[assay2]] <- NULL

    return(obj)}


#' get_org
#'
#' Retrieves the organism from a Seurat object.
#'
#' @param obj A Seurat object
#' @return The organism
#' @export
get_org <- function(obj, assay = "RNA"){
    # validate
    org <- ifelse(all(str_detect(rownames(obj[[assay]]), "[a-z]")), "mouse", "human")
    return(org)
   }