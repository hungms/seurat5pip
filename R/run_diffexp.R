#' Run differential expression analysis on all cells in a Seurat object
#' 
#' This function performs differential expression analysis on all cells in a Seurat object
#' using the specified test.use method.
#' 
#' @param obj A Seurat object
#' @param group.by The column name in the metadata to use for grouping cells
#' @param test.use The test to use for differential expression analysis
#' @param batch The batch variable to use for batch correction
#' @param min.pct The minimum percentage of cells expressing the gene in either group
#' @param logfc.threshold The log fold change threshold for differential expression
#' @param fc.slot The slot to use for fold change calculation
#' @param assay The assay to use for differential expression analysis
#' @param save_data Whether to save the differential expression results
#' @param save_dir The directory to save the differential expression results
#' @param ... Additional arguments to pass to FindAllMarkers
#' 
#' @return A dataframe with differential expression results
#' @export
#' @examples
#' \dontrun{
#' run_diffexp_onetoall(obj, group.by = "seurat_clusters")
#' }
run_diffexp_onetoall <- function(obj, group.by, test.use = "DESeq2", batch = NULL, min.pct = 0.1, logfc.threshold = -Inf, fc.slot = "counts", assay = "RNA", save_data = TRUE, save_dir = getwd(), ...){

    Idents(obj) <- group.by

    markers <- FindAllMarkers(
        obj,
        only.pos = FALSE,
        test.use = test.use,
        min.pct = min.pct,
        fc.slot = fc.slot,
        logfc.threshold = logfc.threshold,
        latent.vars = batch,
        ...)

    markers.list <- markers %>% 
        convert_seurat5_to_deseq2(.) %>%
        split(., .$comparison)

    if(length(batch) > 0){
        batch_label <- paste0("_CORRECTED-FOR-", paste0(batch, collapse = "-"))}
    else {
        batch_label <- ""}

    if(save_data){

        for(i in seq_along(markers.list)){
            save_tsv(markers.list[[i]], tsv_name = paste0("diffexp_", test.use, ".tsv"), save_dir = paste0(save_dir, "/diffexp/", paste0(assay, "_", test.use, "_", "GROUP-BY-", group.by, batch_label, "/", names(markers.list)[i])), row.names = FALSE)}
    
    return(markers)
}

#' Convert Seurat5 differential expression results to DESeq2 format
#' 
#' This function converts differential expression results from Seurat5 to DESeq2 format.
#' 
#' @param diffexp A dataframe with differential expression results
#' @return A dataframe with differential expression results in DESeq2 format
#' @export
#' @examples
#' \dontrun{
#' convert_seurat5_to_deseq2(diffexp)
#' }
convert_seurat5_to_deseq2 <- function(diffexp){
    columns <- c("pvalue", "padj", "log2FoldChange", "gene", "cluster")
    diffexp <- diffexp %>%
        mutate(
            pvalue = p_val, 
            padj = p_val_adj, 
            log2FoldChange = avg_log2FC,
            comparison = ifelse(str_detect(cluster, "_vs_"), cluster, paste0(cluster, "_vs_All"))) %>%
        dplyr::select(-c("p_val", "p_val_adj", "avg_log2FC", "cluster"))
    return(diffexp)}