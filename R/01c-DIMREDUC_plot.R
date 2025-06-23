#' plot_featureloadings
#' 
#' 
#' @param obj seurat object
#' @param reduction reduction name
#' @param output_dir output directory
#' @return ggplot
#' @export

plot_featureloadings <- function(obj, reduction = "pca", output_dir = NULL){
    
    loadings_df <- obj[[reduction]]@feature.loadings %>%
        as.data.frame() %>%
        rownames_to_column(var = "gene") %>%
        pivot_longer(cols = !gene, names_to = "PC", values_to = "loading") %>%
        mutate(
            PC = as.numeric(gsub("PC_", "", PC)),
            PC = factor(PC, levels = rev(1:50)),
            dir = ifelse(loading > 0, "pos", "neg")) %>%
        group_by(PC)
    
    loadings_df_labels <- loadings_df %>%
        group_by(PC, dir) %>%
        slice_max(order_by = abs(loading), n = 10)

    p1 <- loadings_df %>%
        ggplot(aes(x = loading, y = PC)) +
        geom_point(size = 2, aes(color = dir)) +
        scale_color_manual(values = c("pos" = "red", "neg" = "blue")) +
        geom_text_repel(
            data = loadings_df_labels,aes(label = gene), size = 2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
        ggprism::theme_prism(border = T) +
        guides(color = "none") +
        labs(x = "Loading", y = "PC")

    if (!is.null(output_dir)) {
        write_png(p1, output_dir = output_dir, filename = "plot_feature_loadings.png", width = 4000, height = 6000)}

    return(p1)
    
}