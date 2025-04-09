#' Create a dot plot with customizable settings
#'
#' @param df Dataframe with expression values
#' @param x X-axis column name
#' @param y Y-axis column name
#' @param coord_flip Logical; whether the plot coordinates are flipped
#' @param stroke_legend Whether to include significance stroke legend
#' @param palette Color palette for expression values
#' @param scale Logical; whether data was scaled
#'
#' @return A ggplot object
#' @keywords internal
gg_dotplot <- function(df, x, y, coord_flip = TRUE, stroke_legend = FALSE, palette, scale = TRUE) {
    # Set legend direction and titles based on coordinate orientation
    if(coord_flip) {
        direction <- "vertical"
        colorbar_title <- if(scale) "Scaled\nExpression" else "Average\nExpression"
        legend_title_list <- c(colorbar_title, "Percent\nExpressed", "FDR of\nupregulated\ngenes")
        colorbar_height <- 2 # Reduced height for vertical colorbar when coord_flip=TRUE
    } else {
        direction <- "horizontal"
        colorbar_title <- if(scale) "Scaled Expression" else "Average Expression"
        legend_title_list <- c(colorbar_title, "Percent Expressed", "FDR of\nupregulated genes")
        colorbar_height <- 0.5 # Default height for horizontal colorbar
    }
    
    # Set up significance stroke if available
    if ("signif" %in% colnames(df)) {
        stroke <- "signif"
        stroke_guide <- guide_legend(
            direction = direction,
            order = 3,
            keyheight = unit(0.3, "cm"),
            keywidth = unit(0.3, "cm"),
            override.aes = list(fill = palette[6], size = 6),
            title = legend_title_list[3],
            title.position = "top"
        )
    } else {
        stroke <- NULL
        stroke_guide <- guide_none()
    }

    # Create value range for fill scale
    breaks <- c(min(df[["Avg"]]), 0, min(df[["Avg"]])/2, max(df[["Avg"]]))
    
    # Create the dot plot
    ggplot(df, aes_string(x = x, y = y)) +
        geom_point(aes_string(size = "Pct", fill = "Avg", stroke = stroke), shape = 21) +
        scale_discrete_manual(aesthetics = "stroke", values = c(1.3, 0)) +
        scale_color_manual(values = c("black", "white")) +
        scale_fill_gradientn(colors = palette) +
        scale_size(range = c(0, 10), breaks = c(25, 50, 75), 
                  labels = c("25%", "50%", "75%")) +
        labs(x = "", y = "") +
        guides(
            fill = guide_colorbar(
                title = legend_title_list[1],
                title.position = "top",
                direction = direction,
                frame.colour = "black",
                ticks.colour = "black",
                barheight = unit(colorbar_height, "cm"), # Modified colorbar height
                barwidth = if(direction == "vertical") unit(0.5, "cm") else unit(4, "cm"),
                order = 1
            ),
            size = guide_legend(
                direction = direction,
                order = 2,
                override.aes = list(fill = "black"),
                title = legend_title_list[2],
                title.position = "top"
            ),
            stroke = stroke_guide
        ) +
        theme_border() +
        theme(
            axis.text.x = element_text(size = 16, color = "black", 
                                      angle = 90, 
                                      hjust = 1,
                                      vjust = 0.5),
            axis.text.x.top = element_text(size = 16, color = "black",
                                          angle = 90,
                                          hjust = 0,
                                          vjust = 0.5),
            axis.text.y = element_text(size = 16, color = "black"),
            axis.ticks.length = unit(0.15, "cm"))
}

#' Format feature list for dotplot
#'
#' Processes a list or vector of gene names into a formatted dataframe for dot plots
#'
#' @param input Genes as character vector or named list
#' @param unique Logical; whether to keep only unique genes across gene sets
#'
#' @return Dataframe with formatted features
#' @keywords internal
dotplot_convert_genes <- function(input, unique = FALSE) {
    # Convert character vector to list if needed
    if (is.character(input)) {
        input <- list(` ` = input)
    }
    
    # Process list of features
    if (is.list(input)) {
        input_df <- data.frame(
            GeneSet = rep(names(input), sapply(input, length)),
            Gene = unlist(input, use.names = FALSE)
        ) %>%
        mutate(GeneSet = factor(GeneSet, unique(names(input)))) %>%
        group_by(GeneSet) %>%
        mutate(Gene = factor(Gene, rev(unique(Gene)))) %>%
        arrange(Gene) %>%
        ungroup() %>%
        mutate(order = factor(1:nrow(.)))
        
        # Keep only unique genes if requested
        if (unique) {
            input_df <- input_df %>% distinct(GeneSet, Gene, .keep_all = TRUE)
        }
    }
    
    return(input_df)
}

#' Get gene expression data
#'
#' Extracts and summarizes gene expression data from a Seurat object
#'
#' @param obj Seurat object
#' @param genes Vector of genes to include
#' @param group.by Column to group by
#' @param split.by Column to split by (optional)
#' @param assay Assay to use
#' @param slot Data slot to use
#' @param scale Logical; whether to scale data
#'
#' @return Dataframe with summarized expression values
#' @keywords internal
dotplot_get_expression <- function(obj, genes, group.by, split.by = NULL, 
                         assay = "RNA", slot = "data", scale = TRUE) {
    # Extract expression data
    DefaultAssay(obj) <- assay
    df <- FetchData(obj, vars = c(genes), assay = assay, slot = slot)
    genes <- genes[which(genes %in% colnames(df))]
    
    # Scale data if requested
    if (scale) {
        df <- as.data.frame(scale(as.matrix(df)))
    }
    
    # Add grouping variables
    df[["Group"]] <- obj@meta.data[[group.by]]
    df[["Split"]] <- if (!is.null(split.by)) obj@meta.data[[split.by]] else ""
    
    # Reshape and summarize
    exprs_df <- df %>%
        pivot_longer(
            cols = all_of(genes), 
            names_to = "Gene", 
            values_to = "Expression"
        ) %>%
        group_by(Group, Gene, Split) %>%
        summarise(
            Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100,
            .groups = "drop"
        ) %>%
        arrange(Gene)
    
    # Replace NAs with zeros
    exprs_df[is.na(exprs_df)] <- 0
    
    return(exprs_df)
}

#' Add differential expression significance to expression data
#'
#' Adds statistical significance information from differential expression results
#'
#' @param exprs_df Expression dataframe
#' @param diffexp Differential expression results
#' @param split.by Column to split by (optional)
#'
#' @return Annotated expression dataframe
#' @keywords internal
dotplot_add_statistics <- function(exprs_df, diffexp, split.by = NULL) {
    # Load differential expression data if it's a file path
    if (!is.data.frame(diffexp)) {
        stopifnot(file.exists(diffexp))
        diffexp <- read.csv(diffexp, row.names = 1)
    }
    
    # Process differential expression data
    diffexp <- diffexp %>%
        mutate(
            signif = case_when(
                avg_log2FC > 0 & p_val_adj < 0.05 & pct.1 > 0.1 ~ "< 0.05",
                .default = "ns"
            ),
            signif = factor(signif, c("< 0.05", "ns")),
            Gene = gene,
            Group = cluster
        )
    
    # Add split information if provided
    if (!is.null(split.by)) {
        group_cols <- c("Group", "Gene", "Split")
        diffexp$Split <- diffexp[[split.by]]
    } else {
        group_cols <- c("Group", "Gene")
    }
    
    # Merge with expression data
    exprs_df <- exprs_df %>%
        merge(., diffexp, by = group_cols, all.x = TRUE) %>%
        dplyr::select(!c(gene, cluster))
    
    # Replace NAs with "ns"
    exprs_df[is.na(exprs_df)] <- "ns"
    
    return(exprs_df)
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
#' @param genes Gene names as vector or named list
#' @param assay Assay name, defaults to "RNA"
#' @param slot Slot name, defaults to "data"
#' @param scale Logical; if TRUE, scale gene expression across cells
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
#' plot_dotplot(seurat_obj, group.by = "seurat_clusters", genes = c("CD3E", "CD4", "CD8A"))
#'
#' # Using named gene lists and showing differential expression
#' genes_list <- list(
#'   "T cells" = c("CD3E", "CD4", "CD8A"),
#'   "B cells" = c("CD19", "MS4A1")
#' )
#' plot_dotplot(seurat_obj, group.by = "cell_type", genes = genes_list, 
#'              diffexp = "de_results.csv")
#' }
plot_dotplot <- function(
    obj, group.by, split.by = NULL, genes,
    assay = "RNA", slot = "data", scale = TRUE, 
    pal = NULL, diffexp = NULL,
    coord_flip = TRUE, unique = FALSE) {

    # Process genes and get expression data
    genes_df <- dotplot_convert_genes(genes, unique)
    var <- genes_df$Gene
    
    # Get expression summaries
    exprs_df <- dotplot_get_expression(
        obj, var, group.by, split.by, assay, slot, scale)
    
    # Add differential expression information if provided
    if (!is.null(diffexp)) {
        exprs_df <- dotplot_add_statistics(exprs_df, diffexp, split.by)
        stroke_legend <- TRUE
    } else {
        stroke_legend <- FALSE
    }
    
    # Combine feature information with expression data
    exprs_df <- genes_df %>%
        merge(., exprs_df, by = "Gene", all.x = TRUE) %>%
        filter(Gene %in% rownames(obj[[assay]])) %>%
        arrange(order)
    
    # Get color palette
    if(is.null(pal)){
        pal <- plotr::get_palette("Reds", 9)
    }
    
    # Create the plot with appropriate orientation
    if (coord_flip) {
        plot <- gg_dotplot(
            exprs_df, 
            x = "Group", 
            y = "order", 
            stroke_legend = stroke_legend, 
            palette = pal, 
            scale = scale, 
            coord_flip = coord_flip
        ) +
        ggh4x::facet_grid2(
            vars(GeneSet), vars(Split), 
            scales = "free", 
            independent = "none", 
            space = "free", 
            switch = "y"
        ) +
        scale_x_discrete(position = "top") +
        scale_y_discrete(
            position = "right",
            breaks = exprs_df %>% distinct(order, Gene) %>% .$order,
            labels = exprs_df %>% distinct(order, Gene) %>% .$Gene
        ) +
        theme(
            strip.text.y = element_text(face = "bold", size = 16),
            strip.background.x = element_blank(),
            strip.placement = "outside"
        )
    } else {
        exprs_df <- exprs_df %>%
            mutate(Group = factor(Group, levels = rev(levels(Group))))
        plot <- gg_dotplot(
            exprs_df, 
            x = "order", 
            y = "Group", 
            stroke_legend = stroke_legend, 
            palette = pal, 
            scale = scale, 
            coord_flip = coord_flip
        ) +
        ggh4x::facet_grid2(
            vars(Split), vars(GeneSet), 
            scales = "free", 
            independent = "none", 
            space = "free", 
            switch = "y"
        ) +
        scale_x_discrete(
            breaks = exprs_df %>% distinct(order, Gene) %>% .$order,
            labels = exprs_df %>% distinct(order, Gene) %>% .$Gene
        ) +
        theme(
            strip.text.x = element_text(face = "bold", size = 16),
            strip.background.y = element_blank(),
            strip.placement = "outside"
        )
    }
    
    # Add compact legend theme
    plot <- plot + theme(
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10, face = "bold"),
        legend.spacing.x = unit(0, "cm"),
        legend.spacing.y = unit(0, "cm")
    )
    
    return(plot)
}
