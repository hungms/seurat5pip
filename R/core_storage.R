#' Write a matrix to a directory
#'
#' This function writes a matrix to a directory, including the matrix, genes, and barcodes.
#'
#' @param matrix The matrix to write
#' @param dir The directory to write the matrix to
#' @return None. Creates a matrix.mtx, genes.tsv, and barcodes.tsv file in the specified directory
write_matrix <- function(matrix, dir){

    # Validate the directory
    validate_dir(dir)

    # Write the matrix
    writeMM(obj = matrix, file=paste0(dir, "/matrix.mtx"))

    # Write the genes
    write(x = rownames(matrix), file = paste0(dir, "/genes.tsv"))

    # Write the barcodes
    write(x = colnames(matrix), file = paste0(dir, "/barcodes.tsv"))

}