#' Process RNA QC
#'
#' @param obj Seurat object
#' @param assay Assay to validate
#' @param split.by Metadata column to split the object by
#' @param output_dir Directory to save the output
#' @return Seurat object
#' @export
process_RNA_QC <- function(obj, assay = "RNA", split.by = NULL, output_dir = NULL){

    # remove empty droplets & lowly expressed genes
    obj <- qc_by_standard(obj, split.by = split.by, min.features = 200, min.pct = 0.01, output_dir = output_dir)

    # calculate feature percentages
    obj <- qc_find_percentages(obj, assay = assay)

    # identify outlier cells
    obj <- qc_by_mad(obj, assay = assay, split.by = split.by, var = c("percent.mt", "percent.hb"), dev = 5, output_dir = output_dir)

    # plot qc plots per object
    plot_qc(obj, split.by = split.by, output_dir = output_dir)
    plot_qc(obj, qc_mode = "mad", split.by = split.by, output_dir = output_dir)

    return(obj)
}


#process_RNA_CELLCYCLE <- function(obj, assay = "RNA", split.by = NULL, output_dir = NULL){
