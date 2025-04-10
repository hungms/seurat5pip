pkgs <- c(
    "Seurat", # Core single-cell analysis package (v5+)
    "MASS", # Statistical functions
    "scales", # Scale functions for visualization
    "stringr", # String manipulation
    "rlang", # String manipulation
    "stats", # Statistical functions
    "qs", # Quick serialization
    "magrittr", # Pipe operator
    "plotr", # Custom plotting package
    "mascarade" # Cell type annotation package
    )

for(x in pkgs){
    usethis::use_package(x, type = "depends")} #, type = "depends"