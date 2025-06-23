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

#' Write a png file
#'
#' This function writes a png file to a directory.
#'
#' @param plot The plot to write
#' @param output_dir The directory to write the plot to
#' @param filename The name of the plot
#' @param ... Additional arguments to pass to the png function
#' @return None. Creates a png file in the specified directory
#' @export
write_png <- function(plot, output_dir, filename, ...){

    stopifnot(all(is.character(c(output_dir, filename))))
    stopifnot(str_detect(filename, "\\.png$"))
    stopifnot(dir.exists(output_dir))

    if(!is.null(output_dir)){
        dir.create(output_dir, recursive = T)}

    Cairo::CairoPNG(file.path(output_dir, filename), res = 300, ...)
    print(plot)
    dev.off()
}




