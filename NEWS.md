# seurat5pip News

## seurat5pip 0.1.0 (2023-11-01)

### Features

* Initial release of seurat5pip, a package for Seurat v5 pipelines
* Comprehensive tools for single-cell multi-modal data analysis and visualization
* Core visualization functions:
  * `plot_reduction()` - Visualize dimensionality reduction (UMAP, tSNE) with customizable aesthetics
  * `plot_reduction_mask()` - Add cluster outline masks to dimensionality reduction plots
  * `plot_reduction_density()` - Add density contours to dimensionality reduction plots

* Differential expression analysis:
  * `run_diffexp_onetoall()` - Perform differential expression analysis on all cell clusters
  * `select_significant_genes()` - Filter and select significant differentially expressed genes
  * `convert_seurat5_to_deseq2()` - Convert Seurat5 differential expression results to DESeq2 format

* Data processing:
  * `process_seurat()` - Streamlined processing of Seurat objects with sensible defaults
  * `integrate_v4()` - Integration of multiple samples using Seurat v4 methods

* Utility functions:
  * `save_tsv()` - Save data frames as tab-separated files
  * `save_plot()` - Save ggplot objects as PDF files with Cairo support
  * `save_jupyter_plot()` - Specialized plotting for Jupyter notebook environments

### Dependencies

* Compatible with Seurat v5.0.0 and above
* Utilizes MASS for statistical functions
* Optimized for R 4.0.0 and above

### Notes

* This is the first release of the package, intended for use with Seurat v5
* Custom plotting themes and aesthetics for publication-ready visualizations
* Support for various visual customizations through parameters 