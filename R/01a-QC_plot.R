#' Plot QC
#'
#' @param obj Seurat object
#' @param qc_mode QC mode to plot
#' @param split.by Metadata column to split the object by
#' @param output_dir Directory to save the output
#' @return List of plots
#' @export
plot_qc <- function(obj, qc_mode = NULL, split.by = NULL, output_dir = NULL){

    # get columns to check
    split.by <- validate_split.by(split.by, obj)
    var = c("nFeature_RNA", "nCount_RNA", "percent.hb", "percent.mt", "percent.rb", "percent.bcr", "percent.tcr")
    cols_to_check <- c(split.by, var)
    
    stopifnot(all(cols_to_check %in% colnames(obj@meta.data)))
    metadata <- obj@meta.data

    # if qc_mode is not NA, add it to the columns to check
    if(!is.null(qc_mode)){
        qc_mode <- paste0("qc_by_", qc_mode)
        stopifnot(qc_mode %in% colnames(obj@meta.data))
        metadata <- obj@meta.data %>%
            filter(!!sym(qc_mode) == "Pass")}

    # 1. All variables
    plots <- list()
    boxplot <- list()
    scatter <- list()
    levels <- unique(metadata[[split.by]])

    for(i in seq_along(levels)){
        
        p1 <- list()
        # plot each variable
        p1[[i]] <- obj@meta.data %>%
            filter(!!sym(split.by) == levels[i]) %>%
            pivot_longer(cols = all_of(var[3:length(var)]), names_to = "variable", values_to = "value") %>%
            mutate(variable = factor(variable, levels = var[3:length(var)])) %>%
            ggplot(aes(x = variable, y = value)) +
            labs(title = levels[i], x = NULL, y = NULL)

        if(!is.null(qc_mode)){
            p1[[i]] <- p1[[i]] + 
                geom_point(aes(fill = !!sym(qc_mode)), color = "grey60", size = 0.1, position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.85)) + 
                geom_boxplot(aes(fill = !!sym(qc_mode)), outliers = F, width = 0.5, position = position_dodge(width = 0.85), linewidth = 0.9) +
                ggprism::theme_prism(border = T)
                }
        else{
            p1[[i]] <- p1[[i]] + 
                geom_jitter(width = 0.35, size = 0.1, color = "grey60") + 
                geom_boxplot(outliers = F, width = 0.5, linewidth = 0.9) +
                ggprism::theme_prism(border = T)
                }
        boxplot[[i]] <- p1[[i]]

        # 2. nCount vs nFeature
        p2 <- ggplot(metadata, aes(x = log10(nCount_RNA), y = nFeature_RNA)) +
            geom_point(size = 0.1) +
            geom_smooth(method = "lm", color = "red", linewidth = 0.9) +
            ggprism::theme_prism(border = T) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
            stat_cor(method = "pearson", color = "red", label.y = max(metadata$nFeature_RNA) * 1.05) +
            labs(subtitle = if(is.null(qc_mode)) NULL else paste0(qc_mode, " = Pass"))

        # 3. nCount vs percent.mt
        p3 <- ggplot(metadata, aes(x = log10(nCount_RNA), y = percent.mt)) +
            geom_point(size = 0.1) +
            geom_smooth(method = "lm", color = "red", linewidth = 0.9) +
            ggprism::theme_prism(border = T) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
            stat_cor(method = "pearson", color = "red", label.y = max(metadata$percent.mt) * 1.05) +
            labs(subtitle = if(is.null(qc_mode)) NULL else paste0(qc_mode, " = Pass"))

        # 4. percent.rb vs percent.mt
        p4 <- ggplot(metadata, aes(x = percent.rb, y = percent.mt)) +
            geom_point(size = 0.1) +
            geom_smooth(method = "lm", color = "red", linewidth = 0.9) +
            ggprism::theme_prism(border = T) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
            stat_cor(method = "pearson", color = "red", label.y = max(metadata$percent.mt) * 1.05) +
            labs(subtitle = if(is.null(qc_mode)) NULL else paste0(qc_mode, " = Pass"))

        # 5. percent.tcr vs percent.bcr 
        p5 <- ggplot(metadata, aes(x = percent.tcr, y = percent.bcr)) +
            geom_point(size = 0.1) +
            ggprism::theme_prism(border = T) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.075))) +
            stat_cor(method = "pearson", color = "red", label.y = max(metadata$percent.bcr) * 1.05) +
            labs(subtitle = if(is.null(qc_mode)) NULL else paste0(qc_mode, " = Pass"))

        # add title
        title <- ggdraw() + 
            draw_label(levels[i], fontface = 'bold', size = 16, fontfamily = "sans") +
            theme(plot.margin = margin(0, 0, 0, 0))

        plot_rows <- plot_grid(p2, p3, p4, p5, ncol = 4, align = "vh")
        scatter[[i]] <- plot_grid(title, plot_rows, ncol = 1, rel_heights = c(0.1, 1), align = "v")

        boxplot_dir <- paste0(output_dir, "/plot_qc_boxplot", gsub("qc", "", qc_mode))
        scatter_dir <- paste0(output_dir, "/plot_qc_scatter", gsub("qc", "", qc_mode))
        dir.create(boxplot_dir, recursive = T)
        dir.create(scatter_dir, recursive = T)

        write_png(boxplot[[i]], output_dir = boxplot_dir, filename = paste0(levels[i], ".png"), width = 3000, height = 1200)
        write_png(scatter[[i]], output_dir = scatter_dir, filename = paste0(levels[i], ".png"), width = 6000, height = 1200)

        plots[[i]] <- list(p1, p2, p3, p4, p5)
    }

    names(plots) <- levels
    return(plots)

}

#' Plot cell count
#'
#' @param obj Seurat object
#' @param qc_mode QC mode to plot
#' @param split.by Metadata column to split the object by
#' @param output_dir Directory to save the output
#' @return ggplot object
#' @export
plot_qc_cell_count <- function(obj, qc_mode = NULL, split.by = NULL, output_dir = NULL){

    split.by <- validate_split.by(split.by, obj)

    if(!is.null(qc_mode)){
        stopifnot(qc_mode %in% c("manual", "mad"))
        group_cols <- c(split.by, paste0("qc_by_", qc_mode))}
    else{
        group_cols <- c(split.by)}
    
    plot <- obj@meta.data %>%
        group_by_at(group_cols) %>%
        summarize(n = n()) %>%
        ggplot(aes(x = project, y = n)) +
        labs(x = "Project", y = "No. of cells") +
        ggprism::theme_prism(border = F) +
        scale_y_continuous(
            guide = "prism_offset",
            expand = expansion(mult = c(0, 0.075))) +
        guides(y = "prism_offset_minor") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        

    if(!is.null(qc_mode)){
        plot <- plot +
            geom_col(aes(fill = !!sym(paste0("qc_by_", qc_mode))), width = 0.85, stat = "identity", position = "stack", linewidth = 0.9) +
            geom_text(aes(fill = !!sym(paste0("qc_by_", qc_mode)), label = n), position = position_stack(vjust = 0.5))
        plotname = paste0("plot_qc_cell_count_", qc_mode, ".png")}
    else{
        plot <- plot +
            geom_col(width = 0.85, position = "stack", linewidth = 0.9) +
            geom_text(aes(label = n), vjust = -0.5)
        plotname = paste0("plot_qc_cell_count.png")}
    

    if(!is.null(output_dir)){
        group.levels = unique(obj@meta.data[[split.by]])
        w.base = 800
        h.base = 800
        w.increment = 150
        h.increment = 50
        h.scale <- max(nchar(as.character(group.levels)))
        h.scale <- ifelse(h.scale < 10, 10, h.scale)
        w = w.base + (w.increment * length(group.levels))
        h = h.base + (h.increment * h.scale)

        write_png(plot, output_dir = output_dir, plotname, width = w, height = h)}
}


#' Plot feature percentages
#'
#' @param obj Seurat object
#' @param assay Assay name
#' @param threshold Threshold for each feature
#' @param split.by Metadata column to split the object by
#' @param output_dir Directory to save the output
#' @return List of plots
#' @export
plot_feature_percentages <- function(obj, assay = "RNA", threshold = NULL, split.by = NULL, output_dir = NULL){

    assay_var <- paste0(c("nCount_", "nFeature_"), assay)
    var<- c(assay_var, "percent.mt", "percent.hb", "percent.rb", "percent.bcr", "percent.tcr", "percent.mhc")

    split.by <- validate_split.by(split.by,obj)
    stopifnot(all(var %in% colnames(obj@meta.data)))

    if(!is.null(threshold)){
        stopifnot(all(names(threshold) %in% var))
        stopifnot(all(is.numeric(threshold)))}

    plist <- list()


    for(i in seq_along(var)){

        p1 <- obj@meta.data %>% 
            ggplot(aes(x = !!sym(split.by), y = !!sym(var[i]))) +
            geom_violin(aes(fill = !!sym(split.by)), linewidth = 0.9) +
            ggprism::theme_prism(border = T) +
            guides(fill = guide_none()) +
            labs(x = split.by, y = NULL, title = var[i]) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

        if(!is.null(threshold)){
            p1 <- p1 +
                geom_hline(yintercept = threshold[var[i]], color = "black", linetype = "dashed")}


        p2 <- obj@meta.data %>% 
            ggplot(aes(x = !!sym(var[i]), y = !!sym(split.by))) +
            ggridges::geom_density_ridges(aes(fill = !!sym(split.by))) +
            geom_vline(xintercept = threshold[var[i]], color = "red", linetype = "dashed") +
            ggprism::theme_prism(border = T) +
            guides(fill = guide_none()) +
            labs(title = var[i])

        dir.create(paste0(output_dir, "/plot_feature_percentages/"), showWarnings = F, recursive = T)

        if(!is.null(output_dir)){
            group.levels = unique(obj@meta.data[[split.by]])
            w.base = 800
            h.base = 800
            w.increment = 150
            h.increment = 50
            h.scale <- max(nchar(as.character(group.levels)))
            h.scale <- ifelse(h.scale < 10, 10, h.scale)
            w = w.base + (w.increment * length(group.levels))
            h = h.base + (h.increment * h.scale)
            write_png(p1, output_dir = paste0(output_dir, "/plot_feature_percentages/"), filename = paste0(var[i], ".png"), width = w, height = h)

            group.levels = unique(obj@meta.data[[split.by]])
            w.base = 1000
            h.base = 800
            h.increment = 150
            w.increment = 50
            w.scale <- max(nchar(as.character(group.levels)))
            w.scale <- ifelse(w.scale < 10, 10, w.scale)
            w = w.base + (w.increment * w.scale)
            h = h.base + (h.increment * length(group.levels))
            write_png(p2, output_dir = paste0(output_dir, "/plot_feature_percentages/"), filename = paste0(var[i], "_ridge_density.png"), width = w, height = h)
            }

        plist[[i]] <- list(p1, p2)
        names(plist)[i] <- var[i]
    }

    return(plist)
}