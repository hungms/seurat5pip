#' Write a matrix to a directory
#'
#' This function writes a matrix to a directory, including the matrix, genes, and barcodes.
#'
#' @param matrix The matrix to write
#' @param dir The directory to write the matrix to
#' @return None. Creates a matrix.mtx, genes.tsv, and barcodes.tsv file in the specified directory
write_matrix <- function(obj, assay, dir){
    stopifnot(dir.exists(dir))
    layers <- Layers(obj[[assay]])

    for(l in layers){
        dir.create(paste0(dir, "/", l), recursive = T)
        mat <- obj[[assay]][l]
        Matrix::writeMM(obj = mat, file=paste0(dir, "/", l, "/matrix.mtx"))
        write(x = rownames(mat), file = paste0(dir, "/", l, "/genes.tsv"))
        write(x = colnames(mat), file = paste0(dir, "/", l, "/barcodes.tsv")) 
    }   
}