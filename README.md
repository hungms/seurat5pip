# seurat5pip


```
run_seurat5pip(
    csv = csv,                               # csv file containing file paths, run ID
    rna_meta = rna_meta,                     # RNA metadata file
    hto_meta = hto_meta,                     # HTO metadata file
    batch_vars = c("gender"),                # batch effects to correct
    celltype = "B_CELLS"                     # specify celltypes
    qc_mode = "auto",                        # specify qc mode c("auto", "manual")
    qc_manual_thresh = NA                    # specify qc manual threshold
)
```



```
seuratpip
├── 00-LOG                      #
├── 01-RNA
|   ├── checkpoint              # qs, metadata
|   ├── eval                    # evaluations
|   ├── 01a-QC                  # basic cell/feature filtering
|   ├── 01b-CELLCYCLE           # cell cycle phase inference
|   ├── 01c-DIMREDUC            # feature select, pca, umap
|   ├── 01d-INTEGRATION         # integration methods
|   └── 01e-CLUSTERS            # optimal leiden clustering
├── 02-HTO
|   ├── checkpoint              # qs, metadata
|   ├── eval                    # evaluations
|   ├── 02a-DEMULTIPLEX         # demultiplexing methods
|   └── 02b-DIMREDUC            # pca, umap
├── 03-CLEAN
|   ├── checkpoint              # qs, metadata
|   ├── eval                    # evaluations
|   ├── 03a-DECONTAMINATION     # RNA decontamination methods
|   └── 03b-DOUBLET             # doublet detection methods
├── 04-FINAL
|   ├── checkpoint              # qs, metadata
|   ├── eval                    # evaluations
|   ├── 04a-QC                  # basic cell/feature filtering
|   ├── 04b-CELLCYCLE           # cell cycle phase inference
|   ├── 04c-DIMREDUC            # feature select, pca, umap
|   ├── 04d-INTEGRATION         # integration methods - use previous integrated values where possible
|   └── 04e-CLUSTERS            # optimal leiden clustering
├── 05-ADT
|   ├── checkpoint              # qs, metadata
|   ├── eval                    # evaluations
|   ├── 05a-QC                  # basic cell/feature filtering
|   ├── 05b-CELLCYCLE           # cell cycle phase inference
|   ├── 05c-DIMREDUC            # feature select, pca, umap
|   ├── 05d-INTEGRATION         # integration methods
|   └── 05e-CLUSTERS            # optimal leiden clustering  
```

```
01a-QC
|   ├── STANDARD           # qc scatter plots (per sequencing run)
|   └── MAD                # qc boxplot, metdata, stats, threshold (per sequencing run)

02a-QC
|   ├── STANDARD           # qc scatter plots (per sequencing run)
|   └── MAD                # qc boxplot, metdata, stats, threshold (per sequencing run)


```

|   └── 06-ANNOTATION          # label transfer methods
├── 06-VDJ
|   ├── 00-OUTPUT
|   ├── 01-DANDELION           # dandelion processing & clonotype calling
|   └── 02-IMMCANTATION        # dandelion processing & clonotype calling