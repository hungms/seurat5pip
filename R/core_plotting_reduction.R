#' get_embeddings
#'
#' Get embeddings from Seurat object
#' @param obj Seurat object
#' @param reduction Reduction name
#' @return Data frame with embeddings
#' @export
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
#' @param df Data frame
#' @param reduction.type Reduction type
#' @param group.by Metadata column to group by
#' @param split.by Metadata column to split by
#' @param facet.ncol Number of columns for facets
#' @param shuffle Whether to shuffle points
#' @return Data frame with reduction
#' @export
gg_reduction <- function(df, reduction.type = NULL, group.by = NULL, split.by = NULL, facet.ncol = 2, shuffle = T){
    
    # shuffle points
    if(shuffle){
        df <- df %>% sample_frac(1)}

    # add reduction type
    if(is.null(reduction.type)){
        reduction.type <- gsub("(.*)\\s\\d+$", "", colnames(df)[1])}

    # main ggplot
    plot <- ggplot(df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
        scale_x_continuous(expand = c(0.05, 0.05)) +
        scale_y_continuous(expand = c(0.05, 0.05)) +
        theme_border() +
        theme_text() +
        umap_aes(type = reduction.type) +
        theme(legend.key = element_blank())

    # add split
    if(length(split.by) > 0){
        plot <- plot + facet_wrap(as.formula(paste("~", split.by)), ncol = facet.ncol) + facet_aes()}

    return(plot)}

#' gg_point
#' 
#' Add points to ggplot
#' @param plot ggplot object
#' @param df Data frame
#' @param group.by Metadata column to group by  
#' @param pt.size Point size
#' @param pt.alpha Point alpha
#' @param legend.ncol Number of columns for legend
#' @return ggplot object
#' @export
gg_point <- function(plot, df, group.by = NULL, pt.size = 0.8, pt.alpha = 0.5, legend.ncol = 1){
    if(!is.null(group.by)){
        plot <- plot +
            geom_point(data = df, aes_string(color = group.by), size = pt.size, alpha = pt.alpha) +
            guides(color = guide_legend(title = "", override.aes = list(size = 4, alpha = 1, borders = F), ncol = legend.ncol))}
    else{
        plot <- plot +
            geom_point(data = df, color = "grey85", size = pt.size, alpha = pt.alpha)}
    return(plot)}

#' gg_count
#'
#' Add total cell count to ggplot
#' @param plot ggplot object
#' @param df Data frame
#' @param split.by Metadata column to split by
#' @param count.text.size Text size for cell count
#' @return ggplot object
#' @export
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
#' @export
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

#' gg_mask
#' 
#' Add mask to ggplot
#' @param plot ggplot object
#' @param df Data frame
#' @param group.by Metadata column to group by
#' @param mask.size Mask size
#' @return ggplot object
#' @export
gg_mask <- function(plot, df, group.by, mask.size = 1){
    
    # calculate group masks
    df.mask <- df %>%
        filter(
            abs(scale(!!sym(colnames(df)[1]))) <= 4,
            abs(scale(!!sym(colnames(df)[2]))) <= 4)
    maskTable <- generateMask(dims=df.mask[,c(1,2)], clusters=df.mask[[group.by]])

    # plot group masks
    plot <- plot + geom_path(data=maskTable, aes(group=group), linetype = "dashed", size = mask.size)

    return(plot)}

#' gg_features
#' 
#' Add features to ggplot
#' @param plot ggplot object
#' @param df Data frame
#' @param features Features to add
#' @param feature.size Feature size
#' @param pal Colors to use for features
#' @param pt.size Point size
#' @param pt.alpha Point alpha
#' @param legend.ncol Number of columns for legend
#' @param label Whether to label features
#' @param count Whether to count features
#' @param count.groups Whether to count features per group
#' @param mask Whether to mask features
#' @param count.text.size Text size for cell count
#' @param label.size Label size
#' @param mask.size Mask size
#' @return ggplot object
#' @export
#'
#' 
gg_features <- function(
    plot, 
    df, 
    group.by = NULL, 
    split.by = NULL, 
    pal = NULL, 
    pt.size = 0.5, 
    pt.alpha = 0.5, 
    legend.ncol = 1,
    label = T, 
    count = T, 
    count.groups = T, 
    mask = T, 
    count.text.size = 5, 
    label.size = 4, 
    mask.size = 1){

    # add points
    plot <- gg_point(plot, df, group.by = group.by, pt.size = pt.size, pt.alpha = pt.alpha, legend.ncol = legend.ncol)

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
        if(length(names(pal)) == 0){
            names(pal) <- levels(df[[group.by]])}
        plot <- plot + scale_color_manual(values = pal, labels = group_labels)}

    # add mask
    if(mask){
        plot <- gg_mask(plot, df, group.by = group.by, mask.size = mask.size)}

    return(plot)
}

#' plot_reduction
#'
#' Plot reduction from Seurat object
#' @param obj Seurat object
#' @param reduction Reduction name
#' @param reduction.type Reduction type
#' @param group.by Metadata column to group by
#' @param split.by Metadata column to split by
#' @param facet.ncol Number of columns for facets
#' @param shuffle Whether to shuffle points
#' @param ... Additional arguments to gg_features
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_reduction(obj, group.by = "seurat_clusters")
#' }
plot_reduction <- function(obj, reduction = "umap", reduction.type = "UMAP", group.by, split.by = NULL, facet.ncol = 2, shuffle = T, ...){

    # get embeddings
    df <- get_embeddings(obj, reduction)

    # set group.by as factor
    if(!is.factor(df[[group.by]])){
        df[[group.by]] <- factor(df[[group.by]], sort(unique(df[[group.by]])))}

    # main ggplot
    plot <- gg_reduction(df, group.by = group.by, split.by = split.by, reduction.type = reduction.type, facet.ncol = facet.ncol, shuffle = shuffle)

    # add features
    plot <- gg_features(plot, df, group.by = group.by, split.by = split.by, ...)

    # add pca variance
    if(reduction.type == "PCA"){
        pca_variance <- 100*obj[[reduction]]@stdev^2 / sum(obj[[reduction]]@stdev^2)
        plot <- plot +
            labs(x = paste0("PC1 (", round(pca_variance[1], 0), "% variation)"), y = paste0("PC2 (", round(pca_variance[2], 0), "% variation)"))}

    return(plot)}

#' plot_reduction_mapquery()
#' 
#' Plot reduction from Seurat object with mapquery
#' @param obj Seurat object
#' @param group.by Metadata column to group by
#' @param reduction Reduction name
#' @param facet.ncol Number of columns for facets
#' @param shuffle Whether to shuffle points
#' @param ... Additional arguments to plot_reduction
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_reduction_mapquery(obj, group.by = "seurat_clusters")
#' }
plot_reduction_mapquery <- function(obj, ref.obj, group.by = NULL, split.by = NULL, facet.ncol = 2, reduction = "umap", reduction.type = "UMAP", shuffle = T, ...){

    # get ref embeddings
    ref.df <- get_embeddings(ref.obj, reduction)
    
    # split obj embeddings
    if(!is.null(split.by)){
        split.levels <- obj@meta.data[[split.by]]
        if(is.factor(split.levels)){
            split.levels <- levels(split.levels)}
        else{
            split.levels <- sort(unique(split.levels))}

        ref.df.list <- list()
        for(i in seq_along(split.levels)){
            ref.df.list[[i]] <- ref.df
            ref.df.list[[i]][[split.by]] <- split.levels[i]}
        ref.df <- bind_rows(ref.df.list)
        ref.df[[split.by]] <- factor(ref.df[[split.by]], levels = split.levels)
    }

    # main ggplot
    plot <- gg_reduction(ref.df, group.by = NULL, split.by = split.by, reduction.type = reduction.type, facet.ncol = facet.ncol)
    plot <- gg_point(plot, ref.df, group.by = NULL, pt.size = 0.3, pt.alpha = 0.3)

    # get obj embeddings
    obj.df <- get_embeddings(obj, reduction)

    # add features
    plot <- gg_features(plot, obj.df, group.by = group.by, split.by = split.by, ...)
    return(plot)
}

#' plot_reduction_density
#'
#' Plot reduction from Seurat object with density
#' @param obj Seurat object
#' @param split.by Metadata column to split by
#' @param facet.ncol Number of columns for facets
#' @param reduction Reduction name
#' @param reduction.type Reduction type
#' @param adjust Adjustment for density
#' @param pal Color palette
#' @param ... Additional arguments to plot_reduction
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_reduction_density(obj, split.by = "seurat_clusters")
#' }
plot_reduction_density <- function(obj, split.by = NULL, facet.ncol = 2, reduction = "umap", reduction.type = "UMAP", adjust = 5, pal = "Greys",...){
    
    # get embeddings
    df <- get_embeddings(obj, reduction)

    # plot reduction base
    plot <- gg_reduction(df, group.by = NULL, split.by = split.by, reduction.type = reduction.type, facet.ncol = facet.ncol)

    # add points
    plot <- gg_point(plot, df, group.by = NULL, ...)

    # add contour 
    plot <- plot + 
        geom_density_2d(aes(color = after_stat(level)), bins = adjust) +
        scale_color_distiller(palette = pal, direction = 1, name = "Density")
    
    return(plot)
}
