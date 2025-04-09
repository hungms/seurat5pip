#' Select significant genes from a dataframe
#' 
#' This function selects significant genes from a dataframe based on the padj and log2FoldChange columns.
#' 
#' @param df A dataframe with padj and log2FoldChange columns
#' @param thresh A numeric threshold for the padj column
#' @param gene_column A character string specifying the column name of the gene IDs
#' @param n An optional numeric value specifying the number of genes to return
#' 
#' @return A vector of gene IDs
#' @export
#' @examples
#' \dontrun{
#' select_significant_genes(df, thresh = 0.5, gene_column = "gene", n = NULL)
#' }
#' 
select_significant_genes <- function(df, thresh = 0.5, gene_column = "gene", n = NULL){
    df <- df %>%
        filter(padj < 0.05 & log2FoldChange^2 > thresh^2) %>%
        group_by(!!sym(gene_column)) %>%
        mutate(
            direction = ifelse(log2FoldChange > 0, "Up", "Down"),
            count = n()) %>%
        filter(count >= 2 & length(unique(direction)) == 1) %>%
        summarise(
            direction = unique(direction),
            rank = mean(log2FoldChange))

    if(length(n) > 0){
        genes <- df %>%
            group_by(direction) %>%
            slice_max(order_by = rank^2, n = n) %>%
            pull(!!sym(gene_column))}
    else{
        genes <- df %>%
            pull(!!sym(gene_column))
    }
    
    return(genes)}

#' Plot differential expression comparison
#' 
#' This function plots differential expression results between two conditions
#' and highlights significant genes in a scatter plot.
#' 
#' @param df A dataframe with differential expression results
#' @param x The column name of the first condition
#' @param y The column name of the second condition
#' @param thresh A numeric threshold for the padj column
#' @param highlight.genes A vector of gene IDs to highlight in the plot
#' @param tr The label for the top right quadrant
#' @param bl The label for the bottom left quadrant
#' @param tl The label for the top left quadrant
#' @param br The label for the bottom right quadrant
#' 
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_diffexp_comparison(df, x = "log2FoldChange_condition1", y = "log2FoldChange_condition2", thresh = 0.5, highlight.genes = c("gene1", "gene2"), tr = "Both Up", bl = "Both Down", tl = "Opposite", br = "Opposite")
#' }
#' 

plot_diffexp_comparison <- function(df, x, y, thresh = 0.5, highlight.genes = c(), tr = "Both Up", bl = "Both Down", tl = "Opposite", br = "Opposite"){
    df <- df %>%
        mutate(
            quadrant = case_when(
                !!sym(x) > thresh & !!sym(y) > thresh ~ tr,
                !!sym(x) < -thresh & !!sym(y) < -thresh ~ bl,
                !!sym(x) > thresh & !!sym(y) < -thresh ~ tl,
                !!sym(x) < -thresh & !!sym(y) > thresh ~ br,
                .default = "Non-DE"),
            quadrant = factor(quadrant, unique(c(tr, bl, tl, br, "Non-DE"))))

    if(length(highlight.genes) > 0){
        df <- df %>% 
            mutate(
                label = ifelse(gene %in% highlight.genes, gene, ""),
                size = ifelse(label == "", 0.25, 1.5)
                )}
    else{
        df <- df %>% 
            mutate(
                label = ifelse(quadrant != "Non-DE", gene, ""),
                size = ifelse(label == "", 0.25, 1.5))}
        
    df <- df %>%
        group_by(quadrant) %>%
        mutate(pct = ifelse(label != "", paste0(round(n()*100/nrow(df), 1), "%"), "")) %>%
        mutate(
            quad.x = case_when(
                !!sym(x) > thresh & !!sym(y) > thresh ~ max(!!sym(x))*1.2,
                !!sym(x) > thresh & !!sym(y) < -thresh ~ max(!!sym(x))*1.2,
                !!sym(x) < -thresh & !!sym(y) > thresh ~ min(!!sym(x))*1.2,
                !!sym(x) < -thresh & !!sym(y) < -thresh ~ min(!!sym(x))*1.2),
            quad.y = case_when(
                !!sym(x) > thresh & !!sym(y) > thresh ~ max(!!sym(y))*1.2,
                !!sym(x) < -thresh & !!sym(y) > thresh ~ max(!!sym(y))*1.2,
                !!sym(x) > thresh & !!sym(y) < -thresh ~ min(!!sym(y))*1.2,
                !!sym(x) < -thresh & !!sym(y) < -thresh ~ min(!!sym(y))*1.2)
            ) %>%
        ungroup()
    p <- df %>%
        ggplot(aes_string(x = x, y = y)) +
        geom_point(aes(size = size, color = quadrant)) +
        geom_text(aes(label = pct, x = quad.x, y= quad.y), color = "black", size = 5, fontface = "bold") +
        geom_text_repel(aes(label = label), size = 3, force = 20, force_pull = 20) +
        theme_scatter() +
        ylab(y) +
        xlab(x) +
        ggtitle(paste0(x, " vs ", y)) +
        scale_color_manual(values = c("red", "blue", "grey30", "grey60")) +
        scale_size(range=c(0.25, 1.5),breaks = c(0.25, 1.5), labels = c(" ", "Top 25")) +
        scale_x_continuous(expand = c(0.075, 0.075)) +
        scale_y_continuous(expand = c(0.025, 0.025)) +
        guides(size = guide_none(), color = guide_legend(title = "", override.aes = list(size = 5), order = 1))
    return(p)
}