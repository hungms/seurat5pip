#' get_embeddings
#'
#' Get embeddings from Seurat object
#' @param obj Seurat object
#' @param reduction Reduction name
#' @return Data frame with embeddings
#' @keywords internal
get_embeddings <- function(obj, reduction){
    # requirements
    stopifnot(reduction %in% names(obj@reductions))

    # get embeddings
    emb <- as.data.frame(Embeddings(obj@reductions[[reduction]]))
    colnames(emb) <- paste0(toupper(reduction), 1:ncol(emb))
    emb <- merge(emb, obj@meta.data, by = 0, all.x = T) %>% column_to_rownames("Row.names")
    return(emb)}


#' gg_reduction
#'
#' Get reduction from Seurat object
#' @param obj Seurat object
#' @param reduction Reduction name
#' @return Data frame with reduction
#' @keywords internal
gg_reduction <- function(df, group.by = NULL, split.by = NULL, reduction.type = "UMAP", pt.size = 0.5, pt.alpha = 0.5, facet.ncol = 2, legend.ncol = 1, shuffle = T){
        # shuffle points
    if(shuffle){
        df <- df %>% sample_frac(1)}

    # main ggplot
    plot <- ggplot(df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
        scale_x_continuous(expand = c(0.05, 0.05)) +
        scale_y_continuous(expand = c(0.05, 0.05)) +
        theme_border() +
        theme_text() +
        umap_aes(type = reduction.type) +
        theme(legend.key = element_blank())

    # add geom_point
    if(length(group.by) > 0){
        plot <- plot +
            geom_point(aes_string(color = group.by), size = pt.size, alpha = pt.alpha) +
            guides(color = guide_legend(title = "", override.aes = list(size = 4, alpha = 1, borders = F), ncol = legend.ncol))}
    else{
        plot <- plot +
            geom_point(color = "grey80", size = pt.size, alpha = pt.alpha)}

    # add split
    if(length(split.by) > 0){
        plot <- plot + facet_wrap(as.formula(paste("~", split.by)), ncol = facet.ncol) + facet_aes()}

    return(plot)}

#' gg_count
#'
#' Add total cell count to ggplot
#' @param plot ggplot object
#' @param df Data frame
#' @param split.by Metadata column to split by
#' @param count.text.size Text size for cell count
#' @return ggplot object
#' @keywords internal
gg_count <- function(plot, df, split.by = NULL, count.text.size = 5){

    # count cell numbers w/o split
    if(length(split.by) > 0){
        n.total <- df %>% group_by(!!sym(split.by)) %>% summarize(n = paste0("n = ", scales::comma(n()))) %>% ungroup()}
    else{
        n.total <- df %>% summarize(n = paste0("n = ", scales::comma(n())))}
    
    # add cell number to the plot
    n.total[[colnames(df)[1]]] <-  max(df[[colnames(df)[1]]])
    n.total[[colnames(df)[2]]] <-  min(df[[colnames(df)[2]]])
    plot <- plot + 
        geom_text(data = n.total, aes(label = n), size = count.text.size, color = "black", vjust="inward", hjust="inward")
        
    return(plot)}



#' gg_colorlabel
#' 
#' Add labels to ggplot
#' @param plot ggplot object
#' @param df Data frame
#' @param group.by Metadata column to group by
#' @param label Metadata column to label
#' @param label.size Label size
#' @param label.color Label color
#' @return ggplot object
#' @keywords internal
gg_colorlabel <- function(plot, df, group.by, label.size = 4){
    centroids <- df %>%
            group_by(!!sym(group.by)) %>%
            do({
                kde <- kde2d(.[[1]], .[[2]], n = 100) # Kernel density estimation
                max_density <- which(kde$z == max(kde$z), arr.ind = TRUE) # Find peak density
                data.frame(
                    density_x = kde$x[max_density[1]],
                    density_y = kde$y[max_density[2]])}) %>%
            ungroup()

        plot <- plot +
            geom_text(data = centroids, aes(x = density_x, y = density_y, label = !!sym(group.by)), size = label.size, color = "black", fontface = "bold")
    
    return(plot)}


#' plot_reduction
#'
#' Plot reduction from Seurat object
#' @param obj Seurat object
#' @param group.by Metadata column to group by
#' @param reduction Reduction name
#' @param reduction.type Reduction type
#' @param split.by Metadata column to split by
#' @param pal Colors to use for groups
#' @param pt.size Point size
#' @param pt.alpha Point alpha
#' @param count Whether to count cells
#' @param count.text.size Text size for cell count
#' @param count.groups Whether to count cells per group
#' @param label Whether to label groups
#' @param label.size Label size
#' @param facet.ncol Number of columns for facets
#' @param legend.ncol Number of columns for legend
#' @param shuffle Whether to shuffle points
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_reduction(obj, group.by = "seurat_clusters")
#' }
plot_reduction <- function(obj, group.by, reduction = "umap", reduction.type = "UMAP", split.by = NULL, pal = NULL, pt.size = 0.5, pt.alpha = 0.5, count = T, count.text.size = 5, count.groups = F, label = T, label.size = 4, facet.ncol = 2, legend.ncol = 1, shuffle = T){

    # requirements
    stopifnot(all(is.numeric(c(pt.size, pt.alpha, legend.ncol))))
    stopifnot(all(is.logical(c(count, count.groups, shuffle))))


    # get embeddings
    df <- get_embeddings(obj, reduction)

    # set group.by as factor
    if(!is.factor(df[[group.by]])){
        df[[group.by]] <- factor(df[[group.by]], sort(unique(df[[group.by]])))}

    # main ggplot
    plot <- gg_reduction(df, group.by = group.by, split.by = split.by, reduction.type = reduction.type, pt.size = pt.size, pt.alpha = pt.alpha, facet.ncol = facet.ncol, legend.ncol = legend.ncol, shuffle = shuffle)

    # add labels
    if(label){
        plot <- gg_colorlabel(plot, df, group.by = group.by, label.size = label.size)}

    # count total cell numbers
    if(count){
        plot <- gg_count(plot, df, split.by = split.by, count.text.size = count.text.size)}

    # count cell numbers per group
    if(count.groups){
        group_labels <- df %>%
            group_by(!!sym(group.by)) %>%
            mutate(n.group = factor(paste0(!!sym(group.by), " (", n(), ") "))) %>%
            ungroup() %>%
            arrange(!!sym(group.by)) %>%
            .$n.group %>%
            unique(.)}
    else{
        group_labels <- df %>% 
            pull(!!sym(group.by)) %>% 
            sort(.) %>%
            unique(.)}

    # add colors
    if(length(pal) > 0){
        # add names to colors if not provided
        if(length(names(pal)) == 0){
            names(pal) <- levels(df[[group.by]])}
        # add colors and labels
        plot <- plot + scale_color_manual(values = pal, labels = group_labels)}

    return(plot)}


#' plot_reduction_mask
#'
#' Plot reduction from Seurat object with masks
#' @param obj Seurat object
#' @param group.by Metadata column to group by
#' @param reduction Reduction name
#' @param mask.size Mask size
#' @param ... Additional arguments to plot_reduction
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_reduction_mask(obj, group.by = "seurat_clusters")
#' }
plot_reduction_mask <- function(obj, group.by, reduction = "umap", mask.size = 1, ...){

    # get embeddings
    df <- get_embeddings(obj, reduction)

    # plot reduction base
    plot <- plot_reduction(obj = obj, group.by = group.by, reduction = reduction, ...)

    # calculate group masks
    df.mask <- df %>%
        filter(
            abs(scale(!!sym(colnames(df)[1]))) <= 4,
            abs(scale(!!sym(colnames(df)[2]))) <= 4)
    maskTable <- generateMask(dims=df.mask[,c(1,2)], clusters=df.mask[[group.by]])

    # plot group masks
    plot <- plot + geom_path(data=maskTable, aes(group=group), linetype = "dashed", size = mask.size)
    return(plot)
}


#' plot_reduction_density
#'
#' Plot reduction from Seurat object with density
#' @param obj Seurat object
#' @param split.by Metadata column to split by
#' @param reduction Reduction name
#' @param reduction.type Reduction type
#' @param adjust Adjustment for density
#' @param pal Color palette
#' @param ... Additional arguments to plot_reduction
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
plot_reduction_density <- function(obj, split.by = NULL, reduction = "umap", reduction.type = "UMAP", adjust = 5, pal = "Greys",...){
    
    # get embeddings
    df <- get_embeddings(obj, reduction)

    # plot reduction base
    plot <- gg_reduction(df, group.by = NULL, split.by = split.by, ...)

    # add contour 
    plot <- plot + 
        geom_density_2d(aes(color = after_stat(level)), bins = adjust) +
        scale_color_distiller(palette = pal, direction = 1, name = "Density")
    
    return(plot)
}
