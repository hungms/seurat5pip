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
  
  # Capture the calling function's name and arguments (from parent frame)
  call <- sys.call(-1)  # Get call from parent frame instead of current frame

  # Get the function name from the calling function
  func_name <- as.character(call[[1]])

  # Get the arguments from the calling function
  args_list <- as.list(call)[-1]
  
  
  # Get other messages passed to log_function
  other_messages <- list(...)
  
  # Print the log messages
  message("\n#LOG : ", func_name)
  message("#============================================")
  
  # Automatically print each argument's name and value as character strings
  for (arg_name in names(args_list)) {
    arg_value_str <- deparse(args_list[[arg_name]], width.cutoff = 500)
    if(!is.null(arg_value_str)){
      message(paste0("#", arg_name, " = ", arg_value_str))}}
  
  # new line
  message(paste0("\n"))

  # Print other messages
  for (m in other_messages) {
    message(paste0(as.character(m)))}
}


# to capture output
# sink(file_name, append = append)
# sink()