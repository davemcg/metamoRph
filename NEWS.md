# metamoRph 0.2.0

  - two new functions, `select_HVG` and `normalize_data` have been pulled
  out of partially duplicated code in `run_pca` and `metamoRph` to ease future 
  normalization additions and to modularize the code-base
  - irlba added as a prcomp alternative to speed up large single cell dataset 
  dimensionality reduction
  - add a new sample normalization method: Seurat's "LogNormalize"
  - two new scRNA-seq vignettes added

# metamoRph 0.1.3

  - `run_pca` and `metamoRph` now support sparseMatrix
  
# metamoRph 0.1.2

  - add tests with `testthat`
  - tweak `Metadata Label Transfer on Real Data` vignette plots to mark
  tissue grouping with `ggforce::geom_mark_ellipse`
  
# metamoRph 0.1.1

  - improve the vignettes
  - add `extract_prcomp_scaling()`function to yank center/scale information
  from `prcomp` object
  
# metamoRph 0.1.0

  - release