pkgs <- c(
    "Seurat", # Core single-cell analysis package (v5+)
    "MASS", # Statistical functions
    "rlang", # String manipulation
    "stats", # Statistical functions
    "ggpubr", # Publication-ready ggplot2 graphics
    "qs", # Quick serialization
    "magrittr", # Pipe operator
    "plotr", # Custom plotting package
    "deseq2pip", # DESeq2 pipeline
    "strpip", # String manipulation
    "mascarade", # Cell type annotation package
    "ggprism" # Publication-ready ggplot2 graphics
    )

for(x in pkgs){
    usethis::use_package(x, type = "depends")} #, type = "depends"