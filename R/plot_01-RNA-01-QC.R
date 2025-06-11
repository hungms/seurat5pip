plot_qc_by_standard <- function(obj, split.by = NULL, assay = "RNA", output_dir = NULL){

    # average expression
    average.exp <- AverageExpression(obj, assays = assay, features = features, return.seurat = FALSE)

}
