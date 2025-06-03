#' qc standard
#'
#' This function performs standard quality control on a Seurat object.
#' @param obj Seurat object
#' @param split.by Metadata column to split the object by
#' @param assay Assay to validate
#' @param min.features Minimum number of features to keep
#' @param min.pct Minimum percentage of cells expressing a feature
#' @param min.cells Minimum number of cells expressing a feature
#' @return Seurat object
#' @export
qc_standard <- function(obj, split.by = NULL, assay = "RNA", min.features = 200, min.pct = 0.01, min.cells = NULL){

    # validate min.pct and min.cells
    stopifnot(is.null(min.cells) != is.null(min.pct))
    stopifnot(min.pct > 0)

    # validate assays
    obj <- validate_assays(obj)

    # get old dimensions
    old_dim <- dim(obj[[assay]]$counts)

    # split object
    obj.list <- split_obj(obj = obj, split.by = split.by)

    # for each object
    for(o in seq_along(obj.list)){

        # if min.pct
        if(!is.null(min.pct)){
            
            # get features with min.pct expression
            ncells_with_feature <- rowSums(obj.list[[o]][[assay]]$counts > 0)
            ncells <- ncol(obj.list[[o]])
            features.to.keep <- rownames(obj.list[[o]][[assay]]$counts)[ncells_with_feature / ncells >= min.pct]

            # create new assay
            obj.list[[o]][[assay]] <- CreateAssay5Object(
                counts = obj.list[[o]][[assay]]$counts[features.to.keep, ], 
                data = obj.list[[o]][[assay]]$data[features.to.keep, ], 
                min.features = min.features, 
                min.cells = 0)
            
            # min.cells.message
            min.cells.message <- paste0(min.pct, "%")
            }

        # if min.cells
        else if(!is.null(min.cells)){

            # create new assay
            obj.list[[o]][[assay]] <- CreateAssay5Object(
                counts = obj.list[[o]][[assay]]$counts, 
                data = obj.list[[o]][[assay]]$data, 
                min.features = min.features, 
                min.cells = min.cells)
            
            # min.cells.message
            min.cells.message <- min.cells
            }
        }
    
    # merge objects
    if(length(obj.list) > 1){
        obj <- merge(obj.list[[1]], obj.list[c(2:length(obj.list))])}
    else{
        obj <- obj.list[[1]]}

    # get new dimensions
    new_dim <- dim(obj[[assay]]$counts)

    # print message
    m1 <- paste0(old_dim[1] - new_dim[1], " features removed from ",  assay, " assay")
    m2 <- paste0(old_dim[2] - new_dim[2], " cells removed from ",  assay, " assay")
    log_function(m1, m2)

    # return object
    return(obj)}


#' feature percentages
#'
#' This function calculates the percentage of features in the Seurat object.
#' @param obj Seurat object
#' @return Seurat object
#' @export
feature_percentages <- function(obj, assay = "RNA"){

    # get strings to detect MT-, RB-, HB- genes etc...
    str <- grep("\\.string$", ls("package:strpip"), value = TRUE)
    str_dict <- mget(str, envir = asNamespace("strpip"))
    
    # calculate percentage of each category
    for(i in seq_along(str_dict)){
        colname <- paste0("percent.", gsub("\\.string$", "", str[i]))
        obj <- PercentageFeatureSet(obj, pattern = str_dict[[i]], col.name = colname)}
    
    # validate
    stopifnot(any(grepl("^percent\\.", colnames(obj@meta.data))))

    # return object
    return(obj)
    }

#' qc by MAD
#'
#' This function calculates the MAD of features in the Seurat object.
#' @param obj Seurat object
#' @param var Variable columns from meta.data to calculate MAD for
#' @param dev Deviation from median
#' @param output_dir Directory to save MAD threshold table and outcomes
#' @return Seurat object
#' @export
qc_by_mad <- function(obj, split.by = NULL, var = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"), dev = 5, output_dir = NULL){

    # validate variables as numeric
    stopifnot(all(var %in% colnames(obj@meta.data)))
    for(v in var){
        obj@meta.data[[v]] <- as.numeric(obj@meta.data[[v]])}

    # split object
    obj.list <- split_obj(obj = obj, split.by = split.by)

    # for each object
    for(o in seq_along(obj.list)){

        # progress bar
        message("Processing ", o, "/", length(obj.list), " objects")
        log_progress_bar(o, length(obj.list), pb = NULL)

        # get meta.data
        meta <- obj.list[[o]]@meta.data

        # initialize lists
        threshold_list <- list()
        outcome_list <- list()
        
        # calculate MAD for each variable
        for(i in seq_along(var)){

            # get variable
            selected <- meta[[var[i]]]

            # calculate MAD
            median <- median(selected)
            mad <- median(abs(selected - median))
            upper <- median + dev*mad
            lower <- median - dev*mad

            # make dataframe of threshold values
            threshold_list[[i]] <- c(obj = names(obj.list)[o], variable = var[i], median = median, mad = mad, dev = dev, upper = upper, lower = lower)
            outcome_list[[i]] <- ifelse(selected >= (median - dev*mad) & selected <= (median + dev*mad), "Pass", "Fail")
            names(outcome_list)[i] <- paste0("MAD_", var[i])
            }
    
        # combine threshold dataframes
        threshold_df <- bind_rows(threshold_list)

        # combine outcome dataframes
        outcome_list <- as.data.frame(outcome_list)
        outcome_list$MAD_outcome <- apply(outcome_list, 1, function(row) ifelse(all(row == "Pass"), "Pass", "Fail"))
        rownames(outcome_list) <- rownames(meta)
        outcome_df <- bind_rows(outcome_list)[colnames(obj.list[[o]]@meta.data),]

        # save threshold table
        if(!is.null(output_dir)){
            write.csv(threshold_df, file = paste0(output_dir, "/", names(obj.list)[o], "_MAD_filter_thresholds.csv"), row.names = FALSE)
            write.csv(outcome_df, file = paste0(output_dir, "/", names(obj.list)[o], "_MAD_filter_outcomes.csv"), row.names = FALSE)}

        # update meta.data
        obj.list[[o]]@meta.data <- cbind(obj.list[[o]]@meta.data, outcome_df)
        }

    # merge objects
    if(length(obj.list) > 1){
        obj <- merge(obj.list[[1]], obj.list[c(2:length(obj.list))])}
    else{
        obj <- obj.list[[1]]}

    # log
    m1 <- paste("Total cells = ", ncol(obj))
    m2 <- paste("Passed cells = ", sum(obj@meta.data$MAD_outcome == "Pass"))
    m3 <- paste("Failed cells = ", sum(obj@meta.data$MAD_outcome == "Fail"))
    log_function(m1, m2, m3)

    # return object
    return(obj)
}

