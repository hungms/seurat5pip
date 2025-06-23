.onLoad <- function(libname, pkgname) {
  ### load default packages
  packages <- c(
        "MASS", "scales", "stringr", "stats", "ggpubr", "magrittr", "rlang", "plotr", "strpip", "mascarade", "ggprism", "deseq2pip", "qs", "Seurat", "ggprism")
  invisible(lapply(packages, library, character.only = TRUE))

  ### start up settings
  options(dplyr.summarise.inform = FALSE)
  set.seed(123)
  errorMessage <- NULL
}