#' Log progress bar
#'
#' @param i The current iteration
#' @param n The total number of iterations
#' @param pb The progress bar
#' @return The progress bar
#' @export
log_progress_bar <- function(i, n, pb = NULL) {
  if (is.null(pb)) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    setTxtProgressBar(pb, i)
    return(pb)
  } else {
    setTxtProgressBar(pb, i)
    if (i == n) close(pb)
    return(invisible(NULL))
  }
}


#' Log message
#'
#' @return The message
#' @export
log_function <- function(...) {
  
  # Capture the calling function's name and arguments
  call <- match.call()

  # Get the function name
  func_name <- as.character(call[[1]])

  # Get the arguments
  args_list <- as.list(call)[-1]
  
  # Get other messages
  other_messages <- list(...)
  
  # Print the log messages
  message("LOG : ", func_name)
  message("#============================================")
  
  # Automatically print each argument's name and value
  for (arg_name in names(args_list)) {
    message(paste(arg_name, "=", args_list[[arg_name]], sep = " "))}

  # Print other messages
  for (m in other_messages) {
    message(as.character(m))}

  message("#============================================")
}


# to capture output
# sink(file_name, append = append)
# sink()