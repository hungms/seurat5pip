pkgs <- c(
    "tidyverse", "ggrepel", "patchwork", "cowplot", "scales", 
    "grid", "gridExtra", "viridis", "wesanderson", "RColorBrewer", "qs", 
    "readxl", "writexl", "rstatix", "ggpubr", "Seurat", "harmony", "ggh4x", 
    "biomaRt", "BiocParallel", "UCell", "scDblFinder", "dsb",
    "DoubletFinder", "ComplexHeatmap")

for(x in pkgs){
    usethis::use_package(x, type = "depends")} #, type = "depends"