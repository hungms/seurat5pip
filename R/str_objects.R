#' create_seurat_object
#'
#' Creates a list of Seurat objects from 10x Genomics cellranger/h5 outputs.
#' This function processes RNA, ADT (Antibody-Derived Tags), and HTO (Hashtag Oligonucleotides)
#' modalities if present in the cellranger/h5 outputs.
#' 
#' @param dir A vector of directories containing cellranger/h5 outputs
#' @param project A vector of run IDs for each directory
#' @param file_format The format of the data to be loaded. Either "cellranger" or "h5". Defaults to "cellranger"
#' @param hto_prefix Prefix to identify HTO tag names if HTO library is present. Defaults to NULL
#' @return A list of Seurat objects, with each entry corresponding to a sample
#' @export
create_seurat_object <- function(dir, project, file_format = "cellranger", hto_prefix = NULL, merge = TRUE){    # , denoise_adt = TRUE

    # validations
    stopifnot(file_format %in% c("cellranger", "h5")) 
    stopifnot(length(project) == length(dir))

    # initialize list
    obj.list <- list()
    
    # setup for loop
    for(i in seq_along(project)){
        
        # log progress bar
        message(paste0("Processing ", i, " out of ", length(project), " runs..."))
        log_progress_bar(i, length(project), pb = NULL)
        
        # read cellranger output
        if(file_format == "cellranger"){

            # validate directory and cellranger output
            stopifnot(dir.exists(dir[i]))
            stopifnot(c("barcodes.tsv", "features.tsv", "matrix.mtx") %in% list.files(dir[i]))

            # read cellranger output
            matrix_data <- Read10X(dir[i])}

        # read h5 output
        else if(file_format == "h5"){

            # validate directory and h5 file
            stopifnot(file.exists(dir[i]))
            stopifnot(grepl(".h5$", dir[i]))

            # read h5 file
            matrix_data <- Read10X_h5(dir[i])

            # replace underscores with dashes in feature names
            for(k in seq_along(matrix_data)){
                rownames(matrix_data[[k]]) <- gsub("_", "-", rownames(matrix_data[[k]]))}} 


        # check if ADT assay is present
        has_adt <- "Antibody Capture" %in% names(matrix_data)
        if(has_adt){
            rna_counts <- matrix_data$`Gene Expression`
        } else {
            rna_counts <- matrix_data
        }

        # create RNA assay
        obj <- suppressWarnings(CreateSeuratObject(
            as(rna_counts, "sparseMatrix"), 
            assay = "RNA", 
            min.cells = 5, 
            min.features = 200
        ))

        # check if HTO assay is present
        has_hto <- FALSE
        if(has_adt){
            has_hto <- any(str_detect(rownames(matrix_data$`Antibody Capture`), hto_prefix))
            has_adt <- any(!(str_detect(rownames(matrix_data$`Antibody Capture`), hto_prefix)))}

       	# if HTO is present
        if(has_hto){

            # extract HTO from ADT assay
            hto_data <- matrix_data$`Antibody Capture`[which(str_detect(rownames(matrix_data$`Antibody Capture`), hto_prefix)),]

            # create HTO assay
            obj[["HTO"]] <- CreateAssay5Object(
                as(hto_data[, colnames(obj[["RNA"]]$counts)], "sparseMatrix"), 
                min.cells = 0, 
                min.features = 0)}

       	# if ADT is present
        if(has_adt){
            
            # remove HTO from ADT assay
            adt_data <- matrix_data$`Antibody Capture`
            if(has_hto){
                adt_data <- adt_data[-which(str_detect(rownames(adt_data), hto_prefix)), ]}
            
            # create ADT assay
            obj[["ADT"]] <- CreateAssay5Object(
                as(adt_data[, colnames(obj[["RNA"]]$counts)], "sparseMatrix"), 
                min.cells = 0, 
                min.features = 0)
            
            # denoise ADT counts
            #if(denoise_adt){
            #    obj <- run_dsb(
            #        obj, 
            #        dir = dir[i], 
            #        denoise.counts = TRUE, 
            #        use.isotype.control = FALSE, 
            #        isotype.control.name.vec = NULL, 
            #        min.count = NULL)
            # }
        }

       	# add run names to metadata
        obj <- set_assay_keys(obj)
        obj@meta.data$project <- rep(project[i], ncol(obj))
        obj <- RenameCells(
            obj, 
            new.names = paste0(project[i], "_", colnames(obj)))

        # validate object
        obj <- validate_object(obj)

        # add to list
        obj.list[[i]] <- obj
    }

    # add run names to list
    names(obj.list) <- paste0(project)

    # logging
    m1 <- paste("project = ", paste0(project, collapse = ", "))
    m2 <- paste("Number of cells = ", sum(sapply(obj.list, function(x) ncol(x))))
    m3 <- paste("Number of RNA features = ", sum(sapply(obj.list, function(x) if(!"RNA" %in% names(x@assays)){0} else {nrow(x[["RNA"]])})))
    m4 <- paste("Number of ADT features = ", sum(sapply(obj.list, function(x) if(!"ADT" %in% names(x@assays)){0} else {nrow(x[["ADT"]])})))
    m5 <- paste("Number of HTO features = ", sum(sapply(obj.list, function(x) if(!"HTO" %in% names(x@assays)){0} else {nrow(x[["HTO"]])})))
    #m6 <- paste("Number of common features = ", length(intersect_features(obj, split.by = "project")))
    log_function(m1, m2, m3, m4, m5)
    
    # return object
    print(obj.list)

    # merge objects
    if(merge){
        obj <- merge_objects(obj.list)}
    else{
        obj <- obj.list}

    return(obj)
}

#' split_object
#'
#' Splits a Seurat object into a list of Seurat objects based on a metadata column.
#'
#' @param obj A Seurat object.
#' @param split.by A metadata column name to split the object by.
#' @return A list of Seurat objects.
#' @export
split_object <- function(obj, split.by){

    # validate object
    split.by <- validate_split.by(split.by = split.by, obj = obj)

    # split object
    obj.list <- SplitObject(obj, split.by = split.by)

    # return list
    return(obj.list)
}

#' merge_objects
#'
#' Merges a list of Seurat objects into a single Seurat object.
#'
#' @param obj.list A list of Seurat objects.
#' @return A single Seurat object.
#' @export
merge_objects <- function(obj.list){

    # merge objects
    if(length(obj.list) > 1){
        obj <- merge(obj.list[[1]], obj.list[c(2:length(obj.list))], merge.dr = TRUE)}
    else{
        obj <- obj.list[[1]]}

    # join layers
    obj <- join_layers(obj)

    # return object
    return(obj)
}