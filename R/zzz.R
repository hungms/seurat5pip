.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
        "MASS", "scales", "stringr", "stats", "magrittr", "rlang", "plotr", "strpip", "mascarade")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
  errorMessage <- NULL
}