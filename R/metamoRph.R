#' metamoRph
#'
#' This function takes in a count matrix (where genes (features) are rows and samples are
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
#' @param sample_cpm_scale Default is TRUE; This should match the value given
#' to metamoRph::run_pca(). If you are using your own rotations and don't know whether you have
#' done cpm scaling then...you probably should set this to FALSE.
#' @param log1p Default is TRUE; log1p scales the input count matrix
#' @return A matrix with the transformated eigenvalue matrix which should be equivalent
#' to the original rotation matrix's eigenvalue/pattern matrix (The $x slot from
#' the output of prcomp)
#' @import dplyr
#' @import tidyr
#' @import Matrix
#' @importFrom tibble enframe
#' @importFrom stats prcomp
#' @export
metamoRph <- function(new_counts,
                      rotation,
                      center_scale = NULL,
                      sample_cpm_scale = TRUE,
                      log1p = TRUE){
  value <- ensgene <- name <- NULL # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  
  new_cpm <- new_counts # in case no sample scaling is selected
  
  if (sample_cpm_scale){
    message("Sample CPM scaling")
    new_cpm <- scuttle::calculateCPM(new_counts)
  }
  
  if (log1p){
    message("Sample log1p scaling")
    new_cpm <- new_cpm |> log1p()
  }
  
  # extract gene names from new data and make upper case
  row_genes <- row.names(new_cpm) |> toupper()
  row.names(new_cpm) <- row_genes
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
  cutdown <- new_cpm[feature_universe[feature_universe %in% row.names(new_cpm)], , drop = FALSE]
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
    message("Applying custom scaling")
    scaled_cutdown <- scale(Matrix::t(cutdown), center_scale$center,
                            center_scale$scale)
  } else{
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
#' @param sample_cpm_scale Default is TRUE; performs cpm scaling on the samples with the scuttle::scuttle::calculateCPM() function.
#' @param log1p Default is TRUE; applies log1p scaling to the input count matrix.
#' @param remove_regex Default regex pattern is '^MT|^RPS|^RPL'. Set to '' to skip.
#' @return A named list object with the prcomp output returned under the $PCA slot, the given
#' metadata under the $meta slot, the percent variance of each PC as the
#' $percentVar slot, a list object containing the scaled data's "center" and "scale"
#' values for use in the metamoRph function, and the used parameters under the $params slot.\
#' @import Matrix 
#' @importFrom utils head
#' @importFrom matrixStats rowVars
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scran modelGeneVar
#' @importFrom sparseMatrixStats rowVars
#' @importFrom stats prcomp
#' @export
run_pca <- function(feature_by_sample,
                    meta,
                    ntop = 1000,
                    hvg_selection = 'scran',
                    hvg_force = NULL,
                    feature_scale = TRUE,
                    feature_center = TRUE,
                    sample_cpm_scale = TRUE,
                    log1p = TRUE,
                    remove_regex = '^MT|^RPS|^RPL') {
  bio <- Gene <- NULL # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  
  # Remove genes with a regex
  # if given as empty, then skip
  if (remove_regex != ''){
    feature_by_sample <- feature_by_sample[!grepl(remove_regex, row.names(feature_by_sample)),]
  }
  
  # Remove zero count features/genes
  row_sums <- Matrix::rowSums(feature_by_sample)
  feature_by_sample <- feature_by_sample[row_sums > 0,]
  
  if (sample_cpm_scale) {
    feature_by_sample <- scuttle::calculateCPM(feature_by_sample)
  }
  
  if (log1p) {
    feature_by_sample <- feature_by_sample |> log1p()
  }
  
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
  
  # Rotate
  sample_by_feature <- Matrix::t(feature_by_sample[select,])
  
  # finally doing the PCA!
  PCA <- prcomp(sample_by_feature, scale = feature_scale, center = feature_center)
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
                            'sample_cpm_scale' = sample_cpm_scale,
                            'log1p' = log1p,
                            'remove_regex' = remove_regex))
  out
}

#' extract_prcomp_scaling
#'
#' This function takes in a prcomp object and returns the center and scale
#' vectors a list for direct use in metamoRph::metamoRph.
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
