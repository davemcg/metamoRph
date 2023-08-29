#' normalize_data
#'
#' Applies *sample* level normalization and high variable feature (gene)
#' selection. If `cpm` scaling selected, then `scuttle::calculateCPM` is 
#' applied before log1p (if `TRUE`).
#' `zscale` will apply `log1p` first (if `TRUE`).
#' 
#' This function is used by both metamoRph::metamoRph and metamoRph::run_pca
#'
#' @param feature_by_sample features (row) by sample (column) matrix 
#' @param sample_scale Options include `cpm`, `seurat`, `zscale`,  and `none`
#' @param log1p Default is TRUE
#' @return  a list object containing the prcomp "center" and "scale"
#' values for use in the metamoRph function
#' @importFrom scuttle calculateCPM
#' @export
normalize_data <- function(feature_by_sample, sample_scale = 'cpm', log1p = TRUE){
  new_scale <- feature_by_sample # in case no applicable sample scaling is selected
  
  if (sample_scale == 'cpm'){
    message("Sample CPM scaling")
    new_scale <- scuttle::calculateCPM(feature_by_sample)
  } else if (sample_scale == 'seurat'){
    message("Seurat LogNormalize scaling")
    new_scale <- ((feature_by_sample) / 
                    colSums(feature_by_sample)) * 10000
  } else {
    message("no sample scaling")
  }
  
  if (log1p){
    message("log1p scaling")
    new_scale <- new_scale |> log1p()
  } else {
    message("no log1p scaling")
  }
  
  # the standard for zscaling is log1p, then zscale 
  if (sample_scale == 'zscale'){
    message("zscale sample scaling")
    new_scale <- scale(log1p(feature_by_sample))
  } 
  
  return(new_scale)
}

#' select_HVG
#'
#' Takes as input a normalized matrix (see metamoRph::normalize_data).
#' Returns the `ntop` high variable genes either with the "classic" method
#' which orders the features by variance or "scran" which uses the scran package's
#' strategy of scaling variance by value (as larger features/genes) will
#' also have higher variance and thus may be less useful for sample distinction.
#' 
#' This function is used by metamoRph::run_pca
#' 
#' @param feature_by_sample Raw feature (gene) count matrix (where genes/features are rows and samples are
#' columns).
#' @param ntop Number of highly variable genes/features to use in the prcomp PCA. Defaults to 1000.
#' @param hvg_selection Either "classic" or "scran" to select the "ntop" features. "classic"
#' will simply use the top n features by variance, and "scran" will use the scran package's
#' strategy of scaling variance by expression (as highly expressed features/genes) will
#' also have higher variance and thus may be less useful for sample distinction.
#' @param remove_regex Default regex pattern is '^MT|^RPS|^RPL'. Set to '' to skip.
#' @param hvg_force Optional vector of features / genes that must be in the stats::promp
#' @return character vector of the ntop features
#' @importFrom sparseMatrixStats rowVars
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
select_HVG <- function(feature_by_sample,
                       ntop = 1000, 
                       hvg_selection = 'scran',
                       remove_regex = '^MT|^RPS|^RPL',
                       hvg_force = NULL){
  bio <- Gene <- NULL
  # Remove features/genes with a regex
  # if given as empty, then skip
  if (remove_regex != ''){
    feature_by_sample <- feature_by_sample[!grepl(remove_regex, row.names(feature_by_sample)),]
  }
  
  # Remove zero count features/genes
  row_sums <- Matrix::rowSums(feature_by_sample)
  feature_by_sample <- feature_by_sample[row_sums > 0,]
  
  if (hvg_selection == 'classic') {
    if ('dgCMatrix' %in% class(feature_by_sample)) {
      Pvars <- sparseMatrixStats::rowVars(feature_by_sample)
    } else {
      Pvars <- rowVars(feature_by_sample)
    }
    names(Pvars) <- row.names(feature_by_sample)
    select <- sort(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))] |> names()
  } else if (hvg_selection == 'scran') {
    sce_internal <- SingleCellExperiment(list(logcounts = (feature_by_sample)))
    feature_var <- modelGeneVar(sce_internal)
    select <- feature_var |> as_tibble(rownames = 'Gene') |>
      arrange(-bio) |> head(ntop) |> pull(Gene)
  } else {
    stop("Select either classic or scran for hvg_selection please")
  }
  
  if (length(hvg_force) > 0){
    # add in user required features
    select <- c(select, hvg_force)
    # if (!hvg_force %in% row.names(feature_by_sample)){
    #   stop("Requested feature / gene not in your input matrix!")
    # }
  }
  
  return(select)
}



#' metamoRph
#'
#' Takes in a count matrix (where genes (features) are rows and samples are
#' columns) as well as a named vector with the eigenvalues (see [metamoRph::run_pca()])
#' and pulls the gene (feature) information from the rotation vector and cuts down
#' the new_counts matrix to match the rotation vector gene (feature) names. Any
#' genes (features) missing from the input new_counts matrix will be replaced with zeros.
#'
#' The function will scale the new_counts matrix in the same manner as [metamoRph::run_pca()]
#' and matrix multiply by the rotation vector. The output is equivalent
#' to the prcomp "$x" matrix.
#'
#' @param new_counts raw gene count matrix (where genes are rows and samples are
#' columns)
#' @param rotation matrix where the row names are genes and the col names are the
#' principal components. If you used metamoRph::run_pca() then this would be
#' in the output$PCA$rotation slot.
#' @param center_scale list object where the $center slot has the center values
#' and the $scale slot has the scale value for the "scale" function. If you do not
#' give a value here, then feature/gene scaling WILL NOT HAPPEN.
#' @param sample_scale Options include `cpm`, `seurat`, `zscale`,  and `none`. This MUST match
#' the normalization used in the data for `run_pca` (or where your rotation came from). 
#' @param log1p Default is TRUE. Again, must match the normalization used in the data for 
#' `run_pca` (or where your rotation came from).
#' @param normalization Default is TRUE, if set to FALSE will override `sample_scale` and `log1p` 
#' and not do any sample scaling
#' @param feature_scale Default is FALSE, if TRUE will apply feature (gene) scaling to the input (query)
#' data (if center_scale is left empty). VERY DANGEROUS as this will likely put the query data in a different scale than your
#' reference data. Use with deliberate intent.
#' @return A matrix with the transformated eigenvalue matrix which should be equivalent
#' to the original rotation matrix's eigenvalue/pattern matrix (The $x slot from
#' the output of prcomp)
#' @import dplyr
#' @import tidyr
#' @importFrom Matrix t colSums
#' @importFrom tibble enframe
#' @importFrom stats prcomp
#' @export
metamoRph <- function(new_counts,
                      rotation,
                      center_scale = NULL,
                      normalization = TRUE,
                      sample_scale = 'cpm',
                      feature_scale = FALSE,
                      log1p = TRUE){
  value <- ensgene <- name <- NULL # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  
  if (is.null(center_scale)){
    warning(
"No center and scale values provided. Projection may be erroneous 
(unless you already ensured your input data
is scaled the same as the rotation/PCA data input data")
  }
  if (normalization){
    new_scale <- normalize_data(new_counts, 
                                sample_scale = sample_scale,
                                log1p = log1p)
  } else {
    message("Normalization skipped")
    new_scale <- new_counts 
  }
  
  # extract gene names from new data and make upper case
  row_genes <- row.names(new_scale) |> toupper()
  row.names(new_scale) <- row_genes
  # this bit is for the mcgaughey (sc)EiaD data which uses a gene naming
  # scheme that pastes together the common gene name with the ENSEMBL id
  suppressWarnings(feature_id_table <- row.names(rotation) |>
                     toupper() |>
                     enframe() |>
                     separate(value, c("feature_id", "ensgene"), sep = " \\(") |>
                     mutate(ensgene = gsub(")", "", ensgene)) |> dplyr::select(-name))
  overlap_with_ID <- row_genes[row_genes %in% feature_id_table$feature_id]
  overlap_with_ens <- row_genes[row_genes %in% feature_id_table$ensgene]
  
  if ((length(overlap_with_ID) < 2) &
      (length(overlap_with_ens) < 2)){
    stop("Failure to align gene (feature) names, please check your
    input matrix rownames against the names of your rotation vector")
  } else {
    # select column ID type to use
    if (length(overlap_with_ID) >
        length(overlap_with_ens)){
      column_val <- 1
    } else {
      column_val <- 2
    }
  }
  
  # Only genes that match the rotation/eigenvalues genes used
  feature_universe <- feature_id_table |> pull(column_val) |> make.unique()
  feature_universe[feature_universe == ''] = 'X'
  cutdown <- new_scale[feature_universe[feature_universe %in% row.names(new_scale)], , drop = FALSE]
  message(paste0("Aligned ", max(nrow(cutdown)), " features/genes."))
  # add in missing genes
  # by building a matrix composed of the missing genes and the samples
  # then rbind it with the existing data
  not_included <- feature_universe[!feature_universe %in% intersect(row.names(cutdown), feature_universe)]
  data <- matrix(nrow=length(not_included), ncol=ncol(cutdown))
  row.names(data) <- not_included
  colnames(data) <- colnames(cutdown)
  cutdown <- rbind(cutdown, data)
  
  # Match order of feature_ids from reference PCA
  ## ensure the input data is in the
  ## same order as the rotation matrix
  cutdown <- cutdown[feature_universe, , drop = FALSE]
  
  if (!missing(center_scale)){
    message("Applying feature scaling from reference data")
    scaled_cutdown <- scale(Matrix::t(cutdown), center_scale$center,
                            center_scale$scale)
  } else if (feature_scale) {
    message("Running feature scaling on query data")
    scaled_cutdown <- Matrix::t(scale(cutdown))
  } else {
    scaled_cutdown <- Matrix::t(cutdown)
  }
  # Replace NAs with 0
  scaled_cutdown[is.na(scaled_cutdown)] <- 0
  
  # project new data onto the PCA space
  ## matrix multiply the expression of the scaled new data against
  ## the supplied eigenvector
  projected_PC_space <- scaled_cutdown %*% rotation |> as.data.frame()
  
  return((projected_PC_space))
}

#' run_pca
#'
#' This function takes in a count matrix (where genes (features) are rows and samples are
#' columns) and sample level metadata and returns a list object with an R::prcomp calculated
#' object, the metadata, the percent variance explained for each principal component, and
#' the genes (features) chosen for the PCA
#'
#' @param feature_by_sample Raw feature (gene) count matrix (where genes/features are rows and samples are
#' columns).
#' @param meta Metadata for the samples. The rows must match the columns
#' for `feature_by_sample`.
#' @param method Defaults to prcomp, use irlba for large matrices for speed improvement
#' @param irlba_n Default 50. Only used for irlba - will return this many principal components.
#' @param ntop Number of highly variable genes/features to use in the prcomp PCA. Defaults to 1000.
#' @param hvg_selection Either "classic" or "scran" to select the "ntop" features. "classic"
#' will simply use the top n features by variance, and "scran" will use the scran package's
#' strategy of scaling variance by expression (as highly expressed features/genes) will
#' also have higher variance and thus may be less useful for sample distinction.
#' @param hvg_force Optional vector of features / genes that must be in the stats::promp
#' input
#' @param feature_scale Default is TRUE, which means features (genes) are scaled
#' with the R::scale function.
#' @param feature_center Default is TRUE, which means features (genes) are centered
#' @param sample_scale Default is cpm; performs cpm scaling on the samples with the scuttle::calculateCPM() function.
#' @param log1p Default is TRUE; applies log1p scaling to the input count matrix.
#' @param normalization Default is TRUE, if set to FALSE will override `sample_scale` and `log1p` 
#' and not do any sample scaling
#' @param remove_regex Default regex pattern is '^MT|^RPS|^RPL'. Set to '' to skip.
#' @return A named list object with the prcomp output returned under the $PCA slot, the given
#' metadata under the $meta slot, the percent variance of each PC as the
#' $percentVar slot, a list object containing the scaled data's "center" and "scale"
#' values for use in the metamoRph function, and the used parameters under the $params slot.\
#' @importFrom Matrix t
#' @importFrom irlba prcomp_irlba
#' @importFrom utils head
#' @importFrom scran modelGeneVar
#' @importFrom sparseMatrixStats rowVars
#' @importFrom stats prcomp
#' @export
run_pca <- function(feature_by_sample,
                    meta,
                    method = 'prcomp',
                    ntop = 1000,
                    hvg_selection = 'scran',
                    hvg_force = NULL,
                    feature_scale = TRUE,
                    feature_center = TRUE,
                    normalization = TRUE,
                    sample_scale = 'cpm',
                    log1p = TRUE,
                    remove_regex = '^MT|^RPS|^RPL',
                    irlba_n = 50) {
  bio <- Gene <- NULL # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  
  if (normalization){
    normalized_feature_by_sample <- normalize_data(feature_by_sample,
                                                   sample_scale = sample_scale,
                                                   log1p = log1p)
  } else {
    message("Normalization skipped")
    normalized_feature_by_sample <- feature_by_sample
  }
  
  select <- select_HVG(normalized_feature_by_sample,
                       ntop = ntop, 
                       hvg_selection = hvg_selection,
                       remove_regex = remove_regex,
                       hvg_force = hvg_force)
  
  # Rotate
  sample_by_feature <- Matrix::t(normalized_feature_by_sample[select,])
  
  # finally doing the PCA!
  if (method == 'prcomp'){
    message('prcomp start')
    PCA <- prcomp(sample_by_feature, scale = feature_scale, center = feature_center)
    
  } else if (method == 'irlba'){
    
    irlba_n <- min(nrow(sample_by_feature),
                   ncol(sample_by_feature),
                   irlba_n)
    message('irlba start')
    PCA <- irlba::prcomp_irlba(sample_by_feature, n = irlba_n,
                               scale = feature_scale, center = feature_center)
    # add in the feature and sample names
    row.names(PCA$x) <- row.names(sample_by_feature)
    row.names(PCA$rotation) <- colnames(sample_by_feature)
  }
  
  # yank out the center and scale values to place
  # at the top level of the output list object
  center_scale <- extract_prcomp_scaling(PCA)
  # Calculate percent variance explained by each PC
  percentVar <- round(100 * PCA$sdev^2 / sum(PCA$sdev^2), 1)
  
  out <- list(PCA = PCA,
              meta = meta,
              percentVar = percentVar,
              center_scale = center_scale,
              params = list('ntop' = ntop,
                            'hvg_selection' = hvg_selection,
                            'feature_scale' = feature_scale,
                            'sample_scale' = sample_scale,
                            'log1p' = log1p,
                            'remove_regex' = remove_regex))
  out
}

#' extract_prcomp_scaling
#'
#' This function takes in a prcomp object and returns the center and scale
#' vectors a list that can be used in metamoRph::metamoRph.
#'
#' This function is used by metamoRph::run_pca
#' 
#' @param prcomp_object A precomputed prcomp object
#' @return  a list object containing the prcomp "center" and "scale"
#' values for use in the metamoRph function

#' @export
extract_prcomp_scaling <- function(prcomp_object){
  out <- list(center = prcomp_object$center,
              scale = prcomp_object$scale)
  return(out)
}
