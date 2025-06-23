# 
scripts <- list.files("/Users/hungm/The Francis Crick Dropbox/CaladoD/Matthew/github-repos/r-pipelines/seurat5pip/R/")
scripts <- gsub(".R$", " = ", scripts)
write.table(scripts, file = "/Users/hungm/The Francis Crick Dropbox/CaladoD/Matthew/github-repos/r-pipelines/seurat5pip/checklist/R/all_functions.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


