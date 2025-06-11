# seurat5pip


```
run_seurat5pip(
    csv = csv,                               # csv file containing file paths, run ID
    rna_meta = rna_meta,                     # RNA metadata file
    hto_meta = hto_meta,                     # HTO metadata file
    batch_vars = c("gender"),                # batch effects to correct
    celltype = "B_CELLS"                     # specify celltypes
    )
```



```
seuratpip
├── 00-LOGS                    
├── 01-RNA
|   ├── 00-OUTPUT              # normalized expression, metadata, evaluations
|   ├── 01-QC                  # basic cell/feature filtering
|   ├── 02-CELLCYCLE           # cell cycle phase inference
|   ├── 03-DIMREDUC            # feature select, pca, umap
|   ├── 04-INTEGRATION         # integration methods
|   └── 05-CLUSTERS            # optimal leiden clustering
├── 02-HTO
|   ├── 00-OUTPUT              # normalized expression, metadata, evaluations
|   ├── 01-DEMULTIPLEX         # demultiplexing methods
|   └── 02-DIMREDUC            # pca, umap
├── 03-CLEAN
|   ├── 00-OUTPUT              # clean expression, metadata, evaluations
|   ├── 01-DECONTAMINATION     # RNA decontamination methods
|   └── 02-DOUBLET             # doublet detection methods
├── 04-FINAL
|   ├── 00-OUTPUT              # normalized expression, metadata, evaluations
|   ├── 01-QC                  # basic cell/feature filtering
|   ├── 02-CELLCYCLE           # cell cycle phase inference
|   ├── 03-DIMREDUC            # feature select, pca, umap
|   ├── 04-INTEGRATION         # integration methods
|   ├── 05-CLUSTERS            # optimal leiden clustering
|   └── 06-ANNOTATION          # label transfer methods
├── 05-ADT
|   ├── 00-OUTPUT              # normalized expression, metadata, evaluations
|   ├── 01-QC                  # basic cell/feature filtering
|   ├── 02-CELLCYCLE           # cell cycle phase inference
|   ├── 03-DIMREDUC            # feature select, pca, umap
|   ├── 04-INTEGRATION         # integration methods
|   ├── 05-CLUSTERS            # optimal leiden clustering
|   └── 06-ANNOTATION          # label transfer methods
├── 06-VDJ
|   ├── 00-OUTPUT
|   ├── 01-DANDELION           # dandelion processing & clonotype calling
|   └── 02-IMMCANTATION        # dandelion processing & clonotype calling
                  
```