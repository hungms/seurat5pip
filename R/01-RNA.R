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
    p1 <- plot_qc_cell_count(obj, qc_mode = NULL, split.by = split.by, output_dir = paste0(output_dir, "/standard/"))

    # calculate feature percentages
    obj <- qc_find_percentages(obj, assay = assay)
    p2 <- plot_feature_percentages(obj, assay = assay, split.by = split.by, output_dir = paste0(output_dir, "/standard/"))

    # identify outlier cells
    obj <- qc_by_mad(obj, assay = assay, split.by = split.by, var = c("percent.mt", "percent.hb"), dev = 5, output_dir = output_dir)
    p3 <- plot_qc_cell_count(obj, qc_mode = "mad", split.by = split.by,  output_dir = paste0(output_dir, "/mad/"))
    p4 <- plot_percent(obj@meta.data, group.by = split.by, var = "qc_by_mad", filename = "plot_qc_by_mad_percent.png", output_dir = paste0(output_dir, "/mad/"))
    p5 <- plot_qc(obj, qc_mode = "mad", split.by = split.by, output_dir = paste0(output_dir, "/mad/"))

    # plot qc plots per object
    plot_qc(obj, split.by = split.by, output_dir = output_dir)
    plot_qc(obj, qc_mode = "mad", split.by = split.by, output_dir = output_dir)

    return(obj)
}

#' Process cell cycle scoring
#'
#' @param obj Seurat object
#' @param assay Assay to validate
#' @param split.by Metadata column to split the object by
#' @param output_dir Directory to save the output
#' @return Seurat object
process_RNA_CELLCYCLE <- function(obj, assay = "RNA", split.by = NULL, output_dir = NULL){

    # validate arguments
    split.by <- validate_split.by(split.by, obj)

    # cell cycle scoring
    obj <- cc_by_seurat(obj, assay = assay, split.by = split.by)

    # plot cell cycle
    plot_cc(obj, var = "Phase", split.by = split.by, output_dir = output_dir)

    return(obj)
}
