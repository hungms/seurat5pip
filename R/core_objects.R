#' create_seurat_object
#'
#' Creates a list of Seurat objects from 10x Genomics cellranger-multi outputs.
#' This function processes RNA, ADT (Antibody-Derived Tags), and HTO (Hashtag Oligonucleotides)
#' modalities if present in the cellranger outputs.
#' 
#' @param data_dirs A vector of directories containing cellranger outputs
#' @param samples A vector of sample names for each directory
#' @param format The format of the data to be loaded. Either "cellranger" or "h5". Defaults to "cellranger"
#' @param hto_pattern Prefix to identify HTO tag names if HTO library is present. Defaults to NULL
#' @param denoise_adt If TRUE, denoise ADT counts using DSB normalization. Defaults to TRUE
#' @return A list of Seurat objects, with each entry corresponding to a sample
#' @export
create_seurat_object <- function(data_dirs, samples, format = "cellranger", hto_pattern = NULL, denoise_adt = TRUE){    

    stopifnot(format %in% c("cellranger", "h5")) 
    stopifnot(length(samples) == length(data_dirs))

    obj.list <- list()

    for(i in seq_along(samples)){
        message(paste0("\033[1mLoading sample ", i, " for ", samples[i], "\033[0m"))

        # read 10x output
        message("Extracting RNA expression...")
        if(format == "cellranger"){
            matrix_data <- Read10X(data_dirs[i])}
        else if(format == "h5"){
            matrix_data <- Read10X_h5(data_dirs[i])
            for(k in seq_along(matrix_data)){
                rownames(matrix_data[[k]]) <- gsub("_", "-", rownames(matrix_data[[k]]))}} 

        has_adt <- "Antibody Capture" %in% names(matrix_data) # check for ADT
        if(has_adt){
            rna_counts <- matrix_data$`Gene Expression`
        } else {
            rna_counts <- matrix_data
        }
        
        # Create Seurat object with RNA data 
        obj <- suppressWarnings(CreateSeuratObject(
            as(rna_counts, "sparseMatrix"), 
            assay = "RNA", 
            min.cells = 5, 
            min.features = 200
        ))

        # check if other modalities are present
        message("Check if HTO and ADT modalities are present...")
        has_hto <- FALSE
        if(has_adt){
            has_hto <- any(str_detect(rownames(matrix_data$`Antibody Capture`), hto_pattern))
            has_adt <- any(!(str_detect(rownames(matrix_data$`Antibody Capture`), hto_pattern)))}

       	# if HTO is present, add it to the Seurat object
        if(has_hto){
            # extract HTO data
            message("Extracting HTO expression...")
            hto_data <- matrix_data$`Antibody Capture`[which(str_detect(rownames(matrix_data$`Antibody Capture`), hto_pattern)),]
            obj[["HTO"]] <- CreateAssay5Object(
                as(hto_data[, colnames(obj[["RNA"]]$counts)], "sparseMatrix"), 
                min.cells = 0, 
                min.features = 0)}

       	# if ADT is present, add it to the Seurat object
        if(has_adt){
            message("Extracting ADT expression...")
            adt_data <- matrix_data$`Antibody Capture`
            if(has_hto){
                adt_data <- adt_data[-which(str_detect(rownames(adt_data), hto_pattern)), ]}
            
            obj[["ADT"]] <- CreateAssay5Object(
                as(adt_data[, colnames(obj[["RNA"]]$counts)], "sparseMatrix"), 
                min.cells = 0, 
                min.features = 0)
            
            if(denoise_adt){
                obj <- run_dsb(
                    obj, 
                    dir = data_dirs[i], 
                    denoise.counts = TRUE, 
                    use.isotype.control = FALSE, 
                    isotype.control.name.vec = NULL, 
                    min.count = NULL)
            }
        }

       	# add sample id in metadata
        message("Modifying Seurat object metadata...")
        obj <- set_assay_keys(obj)
        obj@meta.data$samples <- rep(samples[i], ncol(obj))
        obj <- RenameCells(
            obj, 
            new.names = paste0(samples[i], "_", colnames(obj)))
        obj.list[[i]] <- obj
    }

    names(obj.list) <- paste0(samples)
    return(obj.list)
}