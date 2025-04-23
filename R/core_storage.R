#' Save Data Frame as TSV File
#'
#' This function saves a data frame as a tab-separated values (TSV) file.
#'
#' @param input Data frame to be saved
#' @param experiment Name of the experiment. Default is NULL
#' @param tsv_name Name of the output TSV file
#' @param save_dir Directory where the TSV file will be saved
#' @param row.names Logical. If TRUE, row names will be included in the output. Default is FALSE
#' @return None. Creates a TSV file in the specified directory
#' @examples
#' \dontrun{
#' # Save data frame without row names
#' save_tsv(my_data, "output.tsv", "results")
#' 
#' # Save data frame with row names
#' save_tsv(my_data, "output.tsv", "results", row.names = TRUE)
#' }
#' @export
save_tsv <- function(input, tsv_name, save_dir, row.names = F){
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)}
    write.table(input, paste0(save_dir, "/", tsv_name), sep = "\t", row.names = row.names, col.names = T, quote = F)}

#' Save ggplot Object as PDF
#'
#' This function saves a ggplot object as a PDF file with specified dimensions.
#' It tries to use Cairo PDF device if available, or falls back to standard PDF if not.
#'
#' @param input ggplot object to be saved
#' @param experiment Name of the experiment. Default is NULL
#' @param plot_name Name of the output PDF file
#' @param save_dir Directory where the PDF file will be saved
#' @param w Width of the PDF in inches
#' @param h Height of the PDF in inches
#' @return The file path (invisibly)
#' @examples
#' \dontrun{
#' # Save plot with default dimensions
#' save_plot(my_plot, "plot.pdf", "figures", w = 8, h = 6)
#' 
#' # Save plot with custom dimensions
#' save_plot(my_plot, "plot.pdf", "figures", w = 10, h = 8)
#' }
#' @export
save_plot <- function(input, plot_name, save_dir, w, h){
    # Create directory if it doesn't exist
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }
    
    # Construct file path
    file_path <- file.path(save_dir, paste0(plot_name))
    
    # Try using Cairo if available, otherwise use standard PDF
    tryCatch({
        if (requireNamespace("Cairo", quietly = TRUE)) {
            Cairo::CairoPDF(file = file_path, width = w, height = h)
        } else {
            pdf(file = file_path, width = w, height = h)
        }
        print(input)
        dev.off()
    }, error = function(e) {
        # If Cairo fails, try standard PDF
        message("Cairo PDF failed, using standard PDF device instead: ", e$message)
        pdf(file = file_path, width = w, height = h)
        print(input)
        dev.off()
    })
    
    # Return the file path invisibly
    invisible(file_path)
}


#' Get the name of an object
#'
#' This function searches the global environment for an object and returns its name.
#'
#' @param obj The object to search for
#' @return The name of the object, or NULL if not found
#' @examples
#' \dontrun{
#' get_name(my_plot)
#' }
#' @export
get_name <- function(obj) {
  var_names <- ls(envir = .GlobalEnv) # List all variable names in the global environment
  for (var in var_names) {
    if (identical(get(var, envir = .GlobalEnv), obj)) {
      return(var)}}
  return(NULL)}


#' Save a plot in Jupyter notebook
#'
#' This function saves a plot as a PDF file in the specified directory.
#' It uses the Cairo PDF device if available, or falls back to standard PDF if not.
#'
#' @param plot The plot to save
#' @param save_dir The directory to save the plot
#' @param width The width of the plot
#' @param height The height of the plot
#' @return None. Creates a PDF file in the specified directory
#' @examples
#' \dontrun{
#' save_jupyter_plot(my_plot, "figures")
#' }
#' @export
save_jupyter_plot <- function(plot, save_dir, width = getOption("repr.plot.width"), height = getOption("repr.plot.height")){
    # Create directory if it doesn't exist
    if(!dir.exists(save_dir)){
        dir.create(save_dir, recursive = TRUE)
    }
    
    # Construct file path
    file_path <- paste0(save_dir, "/", get_name(plot), ".pdf")
    print(plot)

    # Try using Cairo if available, otherwise use standard PDF
    tryCatch({
        if (requireNamespace("Cairo", quietly = TRUE)) {
            Cairo::CairoPDF(file = file_path, width = width, height = height)
        } else {
            pdf(file = file_path, width = width, height = height)
        }
        print(plot)
        dev.off()
    }, error = function(e) {
        # If Cairo fails, try standard PDF
        message("Cairo PDF failed, using standard PDF device instead: ", e$message)
        pdf(file = file_path, width = width, height = height)
        print(plot)
        dev.off()
    })
    
    # Return the file path invisibly
    invisible(file_path)
}

