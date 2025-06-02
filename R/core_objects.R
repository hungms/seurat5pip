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
            has_hto <- any(str_detect(rownames(matrix_data$`Antibody Capture`), hto_pattern))
            has_adt <- any(!(str_detect(rownames(matrix_data$`Antibody Capture`), hto_pattern)))}

       	# if HTO is present
        if(has_hto){

            # extract HTO from ADT assay
            hto_data <- matrix_data$`Antibody Capture`[which(str_detect(rownames(matrix_data$`Antibody Capture`), hto_pattern)),]

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
                adt_data <- adt_data[-which(str_detect(rownames(adt_data), hto_pattern)), ]}
            
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
        obj.list[[i]] <- obj
    }

    # add run names to list
    names(obj.list) <- paste0(project)

    # logging
    m1 <- paste("project = ", paste0(project, collapse = ", "))
    m2 <- paste("Number of cells = ", sum(sapply(obj.list, function(x) ncol(x))))
    m3 <- paste("Number of RNA features = ", sum(sapply(obj.list, function(x) nrow(x[["RNA"]]))))
    m4 <- paste("Number of ADT features = ", sum(sapply(obj.list, function(x) nrow(x[["ADT"]]))))
    m5 <- paste("Number of HTO features = ", sum(sapply(obj.list, function(x) nrow(x[["HTO"]]))))
    m6 <- paste("Number of common features = ", length(intersect_features(obj.list)))
    log_function(m1, m2, m3, m4, m5, m6)
    
    # return object
    print(obj.list)

    # merge objects
    if(merge){
        obj <- merge(obj.list[[1]], obj.list[c(2:length(obj.list))])}

    return(obj)
}

#' split_obj
#'
#' Splits a Seurat object into a list of Seurat objects based on a metadata column.
#'
#' @param obj A Seurat object.
#' @param split.by A metadata column name to split the object by.
#' @return A list of Seurat objects.
split_obj <- function(obj, split.by){
    if(!is.null(split.by)){
        stopifnot(split.by %in% colnames(obj@meta.data))
        obj.list <- SplitObject(obj, split.by = split.by)}
    else{
        obj.list <- list(obj)
        if(project %in% colnames(obj@meta.data)){
            if(length(unique(obj$project)) == 1){
                names(obj.list) <- unique(obj$project)}}
    }
    return(obj.list)
}