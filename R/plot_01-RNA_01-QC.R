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
        plots[[i]] <- list(p1, p2, p3, p4, p5)
    }


    boxplot <- plot_grid(plotlist = boxplot, ncol = 1, align = "vh")
    scatter <- plot_grid(plotlist = scatter, ncol = 1, align = "vh")

    if(!is.null(output_dir) & dir.exists(output_dir)){

        boxplot_name <- paste0(c("/plot_qc_boxplot", qc_mode), collapse = "_")
        scatter_name <- paste0(c("/plot_qc_scatter", qc_mode), collapse = "_")

        scale = length(levels)
        Cairo::CairoPDF(paste0(output_dir, paste0(boxplot_name, ".pdf")), width = 30, height = 20 * scale, res = 300)
        print(boxplot)
        dev.off()

        Cairo::CairoPDF(paste0(output_dir, paste0(scatter_name, ".pdf")), width = 80, height = 20 * scale, res = 300)
        print(scatter)
        dev.off()
    }

    names(plots) <- levels
    return(plots)



}




