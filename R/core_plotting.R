#' Plot percent
#'
#' @param df Data frame
#' @param group.by Metadata column to group by
#' @param var Metadata column to plot
#' @param split.by Metadata column to split the object by
#' @param pal Palette to use
#' @param legend.ncol Number of columns in the legend
#' @param count Whether to show the count of each group
#' @param filename Name of the file to save
#' @param output_dir Directory to save the output
#' @return ggplot object
#' @export
plot_percent <- function(df, group.by, var, split.by = NULL, pal = NULL, legend.ncol = 1, count = T, filename = NULL, output_dir = NULL){
    
    # validate arguments
    stopifnot(c(group.by, var, split.by) %in% colnames(obj@meta.data))
    group <- c(group.by, split.by)

    # get metadata
    metadata <- df
    group.levels <- unique(metadata[[group.by]])

    print(pal)
    
    # plot
    plot <- metadata %>%
        filter(!is.na(!!sym(var))) %>%
        group_by_at(c(group, var)) %>%
        summarize(count = n()) %>%
        group_by_at(group) %>%
        mutate(sum = paste0("n=", sum(count)), pct = count*100/sum(count), pct.round = paste0(round(pct,0), "%"), y.label = 102) %>%
        ggplot(aes_string(y = "pct", x = group.by)) +
        geom_col(aes_string(fill = var), width = 0.85, position = "stack", linewidth = 0.9) +
        guides(fill = guide_legend(title = "", ncol = legend.ncol)) +
        labs(y = "Frequency (%)", x = "") +
        ggprism::theme_prism(border = F) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_y_continuous(
            guide = "prism_minor",
            expand = expansion(mult = c(0, 0.1)),
            breaks = seq(0, 100, 25),
            ) +
        guides(y = "prism_offset_minor")
    
    if(count){
        plot <- plot +
            geom_text(aes_string(y = "y.label", fill = var, label = "sum"), size = 2) + 
            geom_text(aes_string(fill = var, label = "pct.round"), color = "white", position = position_stack(vjust = 0.5), size = 3)
            }

    if(!is.null(pal)){
        plot <- plot +
            scale_fill_manual(values = pal)}

    if(!is.null(split.by)){
        plot <- plot +
            facet_wrap(as.formula(paste0("~ ", split.by)), ncol = 1, axes = "all")}

    # save plot
    if(!is.null(output_dir) & !is.null(filename)){
        w.base = 800
        h.base = 800
        w.increment = 150
        h.increment = 50
        h.scale <- max(nchar(as.character(group.levels)))
        h.scale <- ifelse(h.scale < 10, 10, h.scale)
        w = w.base + (w.increment * length(group.levels))
        h = h.base + (h.increment * h.scale)
        
        
        if(!is.null(split.by)){
            facet.levels = unique(metadata[[split.by]])
            h = h * length(facet.levels)
            }
    
        write_png(plot, output_dir = output_dir, filename = filename, width = w, height = h)}


    return(plot)
}

#' Plot boxplot
#'
#' @param df Data frame
#' @param var Metadata column to plot
#' @param group.by Metadata column to group by
#' @param split.by Metadata column to split the object by
#' @param filename Name of the file to save
#' @param output_dir Directory to save the output
#' @return ggplot object
#' @export
plot_boxplot <- function(df, var, group.by, split.by = NULL, filename = NULL, output_dir = NULL){

    # validate arguments
    stopifnot(c(group.by, var, split.by) %in% colnames(df))
    stopifnot(dir.exists(output_dir))
    group.levels <- unique(df[[group.by]])
    if(!is.null(split.by)){
        facet.levels <- unique(df[[split.by]])}

    # plot
    plot <- ggplot(df, aes(x = !!sym(group.by), y = !!sym(var))) +
        geom_jitter(color = "grey60", size = 0.1, width = 0.25) +
        geom_boxplot(aes(fill = !!sym(group.by)), width = 0.5, linewidth = 0.9, outliers = F) +
        ggprism::theme_prism(border = T) +
        ylim(0, 1) +
        labs(x = group.by, y = var)

    if(!is.null(split.by)){
        plot <- plot +
            facet_wrap(as.formula(paste0("~ ", split.by)), ncol = 1)}


    # save plot
    if(!is.null(output_dir) & !is.null(filename)){
        w.base = 800
        h.base = 800
        w.increment = 150
        h.increment = 50
        h.scale <- max(nchar(as.character(group.levels)))
        h.scale <- ifelse(h.scale < 10, 10, h.scale)
        w = w.base + (w.increment * length(group.levels))
        h = h.base + (h.increment * h.scale)

        
        if(!is.null(split.by)){
            facet.levels = unique(metadata[[split.by]])
            h = h * length(facet.levels)
            }
    
        write_png(plot, output_dir = output_dir, filename = filename, width = w, height = h)}

    return(plot)
}
