#' Format gene list for dotplot
#'
#' Processes a list or vector of gene names into a formatted dataframe for dot plots
#'
#' @param input Features as character vector or named list
#' @param unique Logical; whether to keep only unique features across gene sets
#'
#' @return Dataframe with formatted features
#' @keywords internal
dotplot_format_features <- function(input, unique = FALSE, coord_flip = FALSE) {
    
    stopifnot(is.character(input) | is.list(input))

    # Convert character vector to list if needed
    if(is.character(input)){
        input <- list(` ` = input)}

    # Process list of features
    input_df <- data.frame(
            FeatureSet = rep(names(input), sapply(input, length)),
            Feature = unlist(input, use.names = FALSE)) %>%
        mutate(FeatureSet = factor(FeatureSet, unique(names(input)))) %>%
        group_by(FeatureSet) %>%
        mutate(Feature = factor(Feature, if(coord_flip) rev(unique(Feature)) else unique(Feature))) %>%
        arrange(Feature) %>%
        ungroup() %>%
        mutate(order = factor(1:nrow(.)))
        
    # Keep only unique features if requested
    if(unique){
        input_df <- input_df %>% distinct(FeatureSet, Feature, .keep_all = TRUE)}

    return(input_df)}

#' Validate features
#'
#' Validate features for dotplot
#'
#' @param obj Seurat object
#' @param input_df Dataframe with formatted features
#'
#' @return Dataframe with validated features
#' @keywords internal
dotplot_validate_features <- function(obj, input_df){

    assay_keys <- paste0(get_assay_keys(obj), collapse = "|")
    
    # which features do not have an assay key?
    features.to.validate <- input_df %>% 
        mutate(
            in_assay = ifelse(str_detect(Feature, assay_keys), T, F),
            in_meta = ifelse(Feature %in% colnames(obj@meta.data), T, F)) %>%
        filter(in_meta == F) %>%
        .$Feature

    # get a list of all assay and metadata features
    all_features <- get_all_features(obj)
    meta_features <- colnames(obj@meta.data)
    names(meta_features) <- meta_features
    all_features <- c(all_features, meta_features)
    
    # remove features not present in assay or metadata
    remove_features <- features.to.validate[!features.to.validate %in% all_features]
    message(paste0("Removing ", length(remove_features), " features not present in assay or metadata: ", paste(remove_features, collapse = ", ")))

    # update features with correct assay key
    input_df <- input_df %>% 
        filter(!Feature %in% remove_features)
    
    # Only apply case_when if there are features to edit and the feature list is properly named
    if(length(features.to.validate) > 0 && !is.null(names(all_features))) {
        input_df <- input_df %>%
            mutate(
                Feature = case_when(
                    Feature %in% all_features ~ names(all_features)[match(Feature, all_features)],
                    .default = as.character(Feature)
                )
            )
    }

    return(input_df)
}

#' Get gene expression data
#'
#' Extracts and summarizes gene expression data from a Seurat object
#'
#' @param obj Seurat object
#' @param features Vector of features to include
#' @param group.by Column to group by
#' @param split.by Column to split by (optional)
#' @param scale Logical; whether to scale data
#' @param mode 1 or 2; 1 is to replace them with NA, 2 is to replace them with 0
#'
#' @return Dataframe with summarized expression values
#' @keywords internal
dotplot_get_expression <- function(obj, features, group.by, split.by = NULL, scale = TRUE, mode = 1) {

    # Extract expression data
    assay_keys <- paste0(get_assay_keys(obj), "|")
    exprs <- FetchData(obj, vars = c(features), slot = "data") %>% as.data.frame(.)

    # As probir recommended, remove 0 counts
    exprs[exprs == 0] <- NA
    
    # Add grouping variables
    exprs[["Group"]] <- obj@meta.data[[group.by]]
    exprs[["Split"]] <- if (!is.null(split.by)) obj@meta.data[[split.by]] else ""
    
    # Reshape and summarize
    exprs <- exprs %>%
        as.data.frame(.) %>%
        pivot_longer(
            cols = all_of(features), 
            names_to = "Feature", 
            values_to = "Expression")
    
    # Calculate average expression and percentage of cells expressing each feature
    exprs <- exprs %>%
        group_by(Group, Feature, Split) %>%
        summarise(
            Avg = sum(Expression, na.rm = TRUE) / if(mode == 1) sum(Expression > 0, na.rm = TRUE) else length(Expression), # sum of (scaled) expression divided by number of non-zero expression values
            Pct = sum(Expression > 0, na.rm = TRUE) / length(Expression) * 100,
            .groups = "drop") %>%
        arrange(Feature)
    
    # Replace NAs with zeros
    exprs[is.na(exprs)] <- 0

    # Scale data if requested
    if(scale){
        exprs <- exprs %>%
            group_by(Feature) %>%
            mutate(Avg = scale(Avg)) %>%
            ungroup()}
    
    return(exprs)
}

#' Add differential expression significance to expression data
#'
#' Adds statistical significance information from differential expression results
#' 
#' @param obj Seurat object
#' @param exprs Expression dataframe
#' @param diffexp Differential expression results
#' @param split.by Column to split by (optional)
#'
#' @return Annotated expression dataframe
#' @keywords internal
dotplot_add_statistics <- function(obj, exprs, diffexp, fc.thresh = 0, padj.thresh = 0.05, min.pct = 0.1, split.by = NULL) {

    # Load differential expression data if it's a file path
    if (!is.data.frame(diffexp)) {
        stopifnot(file.exists(diffexp))
        diffexp <- read.csv(diffexp, row.names = 1)}
    
    # Process differential expression data
    diffexp <- diffexp %>%
        mutate(
            signif = case_when(
                avg_log2FC > fc.thresh & p_val_adj < padj.thresh & pct.1 > min.pct ~ "< 0.05",
                .default = "ns"),
            signif = factor(signif, c("< 0.05", "ns")),
            Feature = gene,
            Group = cluster
        )
    
    # Add split information if provided
    if (!is.null(split.by)) {
        group_cols <- c("Group", "Feature", "Split")
        diffexp$Split <- diffexp[[split.by]]
    } else {
        group_cols <- c("Group", "Feature")
    }
    
    # which features do not have an assay key?
    assay_keys <- paste0(get_assay_keys(obj), collapse = "|")
    
    features.to.validate <- diffexp %>% 
        mutate(
            in_assay = ifelse(str_detect(Feature, assay_keys), T, F)) %>%
        .$Feature

    all_features <- get_all_features(obj)

    if(length(features.to.validate) > 0 && !is.null(names(all_features))) {
        diffexp <- diffexp %>%
                mutate(
                    Feature = case_when(
                        Feature %in% all_features ~ names(all_features)[match(Feature, all_features)],
                        .default = as.character(Feature)))}

    # Merge with expression data
    exprs <- exprs %>%
        merge(., diffexp, by = group_cols, all.x = TRUE) %>%
        dplyr::select(!c(gene, cluster))
    
    # Replace NAs with "ns"
    exprs[is.na(exprs)] <- "ns"
    
    return(exprs)
}

#' Create a dot plot with customizable settings
#'
#' @param exprs Dataframe with expression values
#' @param x X-axis column name
#' @param y Y-axis column name
#' @param coord_flip Logical; whether the plot coordinates are flipped
#' @param palette Color palette for expression values
#' @param scale Logical; whether data was scaled
#'
#' @return A ggplot object
#' @keywords internal
gg_dotplot <- function(exprs, x, y, coord_flip = TRUE, palette, scale = TRUE) {
    # Check if significance column exists and has meaningful values
    has_significance <- "signif" %in% colnames(exprs) && !all(exprs$signif == "ns")
    
    # Set common parameters
    direction <- "horizontal"
    colorbar_title <- if(scale) "Scaled Expression" else "Average Expression"
    legend_title_list <- c(colorbar_title, "Percent Expressed", "FDR of\nupregulated genes")
    colorbar_height <- 0.6
    
    # Apply orientation-specific settings
    if(coord_flip) {
        direction <- "vertical"
        colorbar_title <- if(scale) "Scaled\nExpression" else "Average\nExpression"
        legend_title_list[1] <- colorbar_title
        legend_title_list[2] <- "Percent\nExpressed"
        colorbar_height <- 2
    } else {
        exprs <- exprs %>% 
            mutate(
                Group = factor(Group, levels = rev(levels(Group))))}
    
    # Create plot base
    p <- ggplot(exprs, aes_string(x = x, y = y))
    
    # Add appropriate point aesthetics
    if (has_significance) {
        p <- p + 
            geom_point(aes_string(size = "Pct", fill = "Avg", stroke = "signif"), shape = 21) +
            scale_discrete_manual(aesthetics = "stroke", values = c(1.3, 0)) +
            guides(stroke = guide_legend(
                direction = direction,
                order = 3,
                keyheight = unit(0.3, "cm"),
                keywidth = unit(0.3, "cm"),
                override.aes = list(fill = palette[6], size = 6),
                title = legend_title_list[3],
                title.position = "top"
            ))
    } else {
        p <- p + geom_point(aes_string(size = "Pct", fill = "Avg"), color = "white", shape = 21)
    }
    
    # Add common styling (orientation-independent)
    p <- p +
        scale_color_manual(values = c("black", "white")) +
        scale_fill_gradientn(colors = palette) +
        scale_size(range = c(0, 10), breaks = c(25, 50, 75), 
                  labels = c("25%", "50%", "75%")) +
        labs(x = "", y = "") +
        theme_border() +
        theme_text(axis_text_size = 14)
    
    # Add guides (using orientation-dependent parameters)
    p <- p + guides(
        fill = guide_colorbar(
            title = legend_title_list[1],
            title.position = "top",
            direction = direction,
            frame.colour = "black",
            ticks.colour = "black",
            barheight = unit(colorbar_height, "cm"),
            barwidth = if(direction == "vertical") unit(0.5, "cm") else unit(4, "cm"),
            order = 1
        ),
        size = guide_legend(
            direction = direction,
            order = 2,
            override.aes = list(fill = "black"),
            title = legend_title_list[2],
            title.position = "top"
        )
    )

    if(coord_flip){
        p <- p + coord_flip()}
    else{
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
        
    return(p)
}


#' Create dot plot visualization of gene expression
#'
#' Creates a customizable dot plot visualization from a Seurat object, where dot size 
#' represents percentage of cells expressing a gene and color intensity represents 
#' average expression level.
#'
#' @param obj Seurat object
#' @param group.by Metadata column to group by 
#' @param split.by Metadata column to split display (optional)
#' @param features Features as character vector or named list
#' @param scale Logical; if TRUE, scale gene expression across cells
#' @param mode 1 or 2; 1 is to replace them with NA, 2 is to replace them with 0
#' @param pal A vector of hexcode colors for the expression gradient
#' @param diffexp Dataframe or file path for differential expression results
#' @param coord_flip Logical; if TRUE, plot features on y-axis
#' @param unique Logical; if TRUE, do not repeat feature names
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_dotplot(seurat_obj, group.by = "seurat_clusters", features = c("CD3E", "CD4", "CD8A"))
#'
#' # Using named gene lists and showing differential expression
#' features_list <- list(
#'   "T cells" = c("CD3E", "CD4", "CD8A"),
#'   "B cells" = c("CD19", "MS4A1")
#' )
#' plot_dotplot(seurat_obj, group.by = "cell_type", features = features_list, 
#'              diffexp = "de_results.csv")
#' }
plot_dotplot <- function(
    obj, group.by, split.by = NULL, 
    features, diffexp = NULL, mode = 1,
    scale = TRUE, pal = NULL, coord_flip = TRUE, 
    unique = FALSE) {

    # Process data pipeline - orientation-independent
    input_df <- dotplot_format_features(features, unique, coord_flip)
    input_df <- dotplot_validate_features(obj, input_df)
    feature_names <- input_df$Feature
    exprs <- dotplot_get_expression(obj, feature_names, group.by, split.by, scale, mode = mode)
    
    # Add differential expression if provided
    if (!is.null(diffexp)) {
        exprs <- dotplot_add_statistics(obj, exprs, diffexp, split.by=split.by)
    }
    
    # Merge feature information with expression data
    exprs <- input_df %>%
        merge(., exprs, by = "Feature", all.x = TRUE) %>%
        arrange(order)

    # Format feature labels - orientation-independent
    assay_keys <- paste0(get_assay_keys(obj), collapse ="|")
    breaks <- exprs %>% distinct(order, Feature) %>% .$order
    labels <- exprs %>% distinct(order, Feature) %>% .$Feature %>% gsub(assay_keys, "", .)
    
    # Set color palette
    if(is.null(pal)){
        pal <- plotr::get_palette("Reds", 9)
    }
    
    # Default parameters (non-flipped orientation)
    x_var <- "order"
    y_var <- "Group"
    row_var <- "FeatureSet"
    col_var <- "Split"
    
    x_scale <- scale_x_discrete(position = ifelse(coord_flip, "top", "bottom"), breaks = breaks, labels = labels)
    y_scale <- scale_y_discrete(position = if (is.vector(features)) "left" else "right")
    facet_strip_theme <- theme(
        strip.text.x = element_text(face = "bold", size = 14),
        strip.text.y = element_text(face = "bold", size = 16),
        strip.background.x = if (is.list(features)) element_rect() else element_blank(),
        strip.background.y = element_blank(),
        strip.placement = "outside"
    )
    
    # Handle orientation-specific settings in a single block
    if(coord_flip){
        # Modify parameters for flipped orientation
        row_var <- "Split"
        col_var <- "FeatureSet"
        facet_strip_theme <- theme(
            strip.text.x = element_text(face = "bold", size = 16),
            strip.text.y = element_text(face = "bold", size = 14),
            strip.background.x = element_blank(),
            strip.background.y = if (is.list(features)) element_rect() else element_blank(),
            strip.placement = "outside"
        )}

    # Create the plot
    plot <- gg_dotplot(
        exprs, 
        x = x_var, 
        y = y_var, 
        palette = pal, 
        scale = scale, 
        coord_flip = coord_flip)
    
    # Apply faceting
    plot <- plot + ggh4x::facet_grid2(
        vars(!!sym(col_var)), vars(!!sym(row_var)),
        scales = "free",
        independent = "none",
        space = "free",
        switch = "y"
    )
    
    # Apply scales
    plot <- plot + x_scale
    if(!is.null(y_scale)) {
        plot <- plot + y_scale}
    
    # Apply facet strip theme and legend formatting
    plot <- plot + 
        facet_strip_theme +
        theme(
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 10, face = "bold"),
            legend.spacing.x = unit(0, "cm"),
            legend.spacing.y = unit(0, "cm")
        )
    
    return(plot)
}
