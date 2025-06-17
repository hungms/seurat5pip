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
qc_by_standard <- function(obj, split.by = NULL, assay = "RNA", min.features = 200, min.pct = NULL, min.cells = 5, output_dir = NULL){

    # validate min.pct and min.cells
    stopifnot(is.null(min.cells) != is.null(min.pct))
    stopifnot(min.pct > 0)

    # get old dimensions
    old_dim <- dim(obj[[assay]]$counts)
    obj@meta.data[[paste0("nCount_", assay)]] <- NULL
    obj@meta.data[[paste0("nFeature_", assay)]] <- NULL

    # split object
    obj.list <- split_object(obj = obj, split.by = split.by)

    # for each object
    for(o in seq_along(obj.list)){

        # if min.pct
        if(!is.null(min.pct)){
            
            # get features with min.pct expression
            ncells_with_feature <- rowSums(as.matrix(obj.list[[o]][[assay]]$counts) > 0)
            ncells <- ncol(obj.list[[o]])
            features.to.keep <- rownames(obj.list[[o]][[assay]]$counts)[(ncells_with_feature / ncells) >= min.pct]

            # create new assay
            obj.list[[o]][[assay]] <- CreateAssay5Object(
                counts = obj.list[[o]][[assay]]$counts[features.to.keep, ], 
                #data = obj.list[[o]][[assay]]$data[features.to.keep, ], 
                min.cells = 0)
            
            # min.cells.message
            min.cells.message <- paste0(min.pct, "%")
            }

        # if min.cells
        else if(!is.null(min.cells)){

            # create new assay
            obj.list[[o]][[assay]] <- CreateAssay5Object(
                counts = obj.list[[o]][[assay]]$counts, 
                #data = obj.list[[o]][[assay]]$data, 
                min.cells = min.cells)
            
            # min.cells.message
            min.cells.message <- min.cells
            }
        }
    
    # merge objects
    obj <- merge_objects(obj.list)

    # subset object
    obj <- subset(obj, subset = nFeature_RNA >= min.features)

    # get new dimensions
    new_dim <- dim(obj[[assay]]$counts)

    # generate relevant plots
    #plot_qc_by_standard(obj = obj, split.by = split.by, assay = assay, output_dir = output_dir)
    
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
#' @param assay Assay to validate
#' @return Seurat object
#' @export
qc_find_percentages <- function(obj, assay = "RNA"){

    # get strings to detect MT-, RB-, HB- genes etc...
    str <- grep("\\.string$", ls("package:strpip"), value = TRUE)
    str_dict <- mget(str, envir = asNamespace("strpip"))
    
    # calculate percentage of each category
    for(i in seq_along(str_dict)){
        colname <- paste0("percent.", gsub("\\.string$", "", str[i]))
        obj <- PercentageFeatureSet(obj, pattern = str_dict[[i]], col.name = colname)}
    
    # validate
    stopifnot(any(grepl("^percent\\.", colnames(obj@meta.data))))

    # log
    log_function()

    # return object
    return(obj)
    }

#' qc by MAD
#'
#' This function calculates the MAD of features in the Seurat object.
#' @param obj Seurat object
#' @param assay Assay to validate
#' @param split.by Metadata column to split the object by
#' @param var Variable columns from meta.data to calculate MAD for
#' @param dev Deviation from median
#' @param output_dir Directory to save MAD threshold table and outcomes
#' @return Seurat object
#' @export
qc_by_mad <- function(obj, assay = "RNA", split.by = NULL, var = c("percent.mt", "percent.hb"), dev = 5, output_dir = NULL){

    # validate variables as numeric
    split.by <- validate_split.by(split.by = split.by, obj = obj)
    stopifnot(all(var %in% c("percent.mt", "percent.hb", "percent.rb", "percent.tcr", "percent.bcr")))
    var <- c(paste0(c("nFeature_", "nCount_"), assay), var)
    stopifnot(all(var %in% colnames(obj@meta.data)))
    for(v in var){
        obj@meta.data[[v]] <- as.numeric(obj@meta.data[[v]])}

    # split object
    obj.list <- split_object(obj = obj, split.by = split.by)

    # initialize list
    agg_threshold <- list()

    # for each object
    for(o in seq_along(obj.list)){

        # get meta.data
        meta <- obj.list[[o]]@meta.data

        # initialize lists
        list_threshold <- list()
        list_outcome <- list()
        
        # calculate MAD for each variable
        for(i in seq_along(var)){   

            # get variable
            selected <- meta[[var[i]]]

            # calculate MAD
            median <- median(selected, na.rm = T)
            mad <- median(abs(selected - median))
            upper <- median + dev*mad
            lower <- median - dev*mad

            # make dataframe of threshold values
            list_threshold[[i]] <- data.frame(obj = names(obj.list)[o], variable = var[i], median = median, mad = mad, dev = dev, upper = upper, lower = lower)
            colnames(list_threshold[[i]])[1] <- split.by

            list_outcome[[i]] <- rep("Fail", length(selected))
            list_outcome[[i]] <- ifelse(selected >= (lower) & selected <= (upper), "Pass", list_outcome[[i]])
            stopifnot(all(list_outcome[[i]] %in% c("Pass", "Fail")))
            names(list_outcome)[i] <- paste0(var[i], "_by_mad")
            }
    
        # combine threshold dataframes
        agg_threshold[[o]] <- bind_rows(list_threshold)

        # combine outcome dataframes
        df_outcome <- as.data.frame(list_outcome)
        df_outcome$qc_by_mad <- apply(df_outcome, 1, function(row) ifelse(all(row == "Pass"), "Pass", "Fail"))
        rownames(df_outcome) <- rownames(meta)

        # update meta.data
        obj.list[[o]]@meta.data <- cbind(obj.list[[o]]@meta.data, df_outcome)
        }

    # merge objects
    obj <- merge_objects(obj.list)

    # get stats
    var <- c(paste0(var, "_by_mad"), "qc_by_mad")
    list_stats <- list()  # Initialize the list_stats variable
    
    for(i in seq_along(var)){
        stopifnot(all(obj@meta.data[[var[i]]] %in% c("Pass", "Fail")))
        list_stats[[i]] <- obj@meta.data %>%
            group_by(!!sym(split.by), !!sym(var[i])) %>%
            summarize(count = n()) %>%
            group_by(!!sym(split.by)) %>%
            mutate(
                !!sym(var[i]) := factor(!!sym(var[i]), levels = c("Pass", "Fail")),
                Variable = gsub("_by_mad", "", var[i]),
                Total = sum(count)) %>%
            ungroup() %>%
            complete(!!sym(var[i]), nesting(!!sym(split.by)), fill = list(count = 0)) %>%
            ungroup() %>%
            pivot_wider(names_from = !!sym(var[i]), values_from = count) %>%
            rowwise() %>%
            mutate(
                `Pass (%)` = round(Pass / Total * 100, 2),
                `Fail (%)` = round(Fail / Total * 100, 2),
                Flag = ifelse(Pass / Total < 0.8, T, F))
    }

    stats <- bind_rows(list_stats)
    
    # save threshold table
    if(!is.null(output_dir)){
        write.csv(bind_rows(agg_threshold), file = paste0(output_dir, "/qc_by_mad_threshold.csv"), row.names = FALSE)
        write.csv(stats, file = paste0(output_dir, "/qc_by_mad_stats.csv"), row.names = FALSE)
        mad_cols <- colnames(obj@meta.data)[grepl("_by_mad", colnames(obj@meta.data))]
        write.csv(obj@meta.data %>% dplyr::select(all_of(mad_cols)), file = paste0(output_dir, "/qc_by_mad_metadata.csv"), row.names = T)}

    # generate relevant plots
    #plot_qc_by_mad(obj = obj, split.by = split.by, assay = assay, output_dir = output_dir)

    # log
    m1 <- paste("Total cells = ", ncol(obj))
    m2 <- paste("Passed cells = ", sum(obj@meta.data$qc_by_mad == "Pass"))
    m3 <- paste("Failed cells = ", sum(obj@meta.data$qc_by_mad == "Fail"))
    log_function(m1, m2, m3)

    # return object
    return(obj)
}


