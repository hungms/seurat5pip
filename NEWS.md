# seurat5pip News

## v1.0.0 
### Official Release of seurat5pip 


# Development Stage
## v0.1.1
### Overview
* Rebuild pipeline structure for all plot functions for readability
* Implement core functionalities for Seurat objects
* Rename files from `R/` directories for pipeline development

### New Features
* New core functionalities to modify Seurat Object: 
  * assays e.g. `set_assay_keys()`, `get_assay_keys`, `join_layers()`
  * features e.g. `remove_features()`, `return_features()`, `intersect_features()`
* `plot_reduction_mapquery()` to visualize query data on original reference map dataset

### Enhancements
* Rebuild structure for `plot_reduction()` and `plot_dotplot()`
* `plot_dotplot()` now remove cells with 0 expression of selected genes to calculate average expression and improve color visualization

## v0.1.0
### Overview
* This is the first development stage of the package, intended for use with Seurat v5
* Custom plotting themes and aesthetics for publication-ready visualizations
* Support for various visual customizations through parameters 

### New Features

* Initial release of seurat5pip, a package for Seurat v5 pipelines
* Comprehensive tools for single-cell multi-modal data analysis and visualization
* Core visualization functions:
  * `plot_reduction()` - Visualize dimensionality reduction (UMAP, tSNE) with customizable aesthetics
  * `plot_reduction_mask()` - Add cluster outline masks to dimensionality reduction plots
  * `plot_reduction_density()` - Add density contours to dimensionality reduction plots


* Data processing:
  * `process_seurat()` - Streamlined processing of Seurat objects with sensible defaults
  * `integrate_v4()` - Integration of multiple samples using Seurat v4 methods

* Utility functions:
  * `save_tsv()` - Save data frames as tab-separated files
  * `save_plot()` - Save ggplot objects as PDF files with Cairo support
  * `save_jupyter_plot()` - Specialized plotting for Jupyter notebook environments
