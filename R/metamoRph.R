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
#' @param cpm_scale Default is TRUE; This should match the value given
#' to metamoRph::run_pca(). If you are using your own rotations and don't know whether you have
#' done cpm scaling then...you probably should set this to FALSE.
#' @param log1p Default is TRUE; log1p scales the input count matrix
#' @return A matrix with the transformated eigenvalue matrix which should be equivalent
#' to the original rotation matrix's eigenvalue/pattern matrix (The $x slot from
#' the output of prcomp)
#' @import dplyr
#' @import tidyr
#' @importFrom tibble enframe
#' @importFrom stats prcomp
#' @export
metamoRph <- function(new_counts,
                      rotation,
                      center_scale = NULL,
                      cpm_scale = TRUE,
                      log1p = TRUE){
  value <- ensgene <- name <- NULL # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  if (cpm_scale){
    new_cpm <- edgeR::cpm(new_counts, log = FALSE)
  }

  if (log1p){
    new_cpm <- new_cpm |> log1p()
  }

  # extract gene names from new data
  ## remove .digit endings (if they are on ensgene then trouble if the
  ## new ensgene are a different version than the old one)
  row_genes <- gsub('\\.\\d+','',row.names(new_cpm)) |> toupper()
  row.names(new_cpm) <- row_genes
  # this bit is for the mcgaughey (sc)EiaD data which uses a gene naming
  # scheme that pastes together the common gene name with the ENSEMBL id
  suppressWarnings(feature_id_table <- row.names(rotation) |> enframe() |>
                     separate(value, c("feature_id", "ensgene"), sep = " \\(") |>
                     mutate(ensgene = gsub(")", "", ensgene)) |> dplyr::select(-name))
  overlap_with_ID <- row_genes[row_genes %in% feature_id_table$feature_id]
  overlap_with_ens <- row_genes[row_genes %in% feature_id_table$ensgene]
  print(paste0("Aligned ", max(length(overlap_with_ID), length(overlap_with_ens)), " features/genes."))
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

  scaled_cutdown <- t(cutdown)
  if (!missing(center_scale)){
    # custom feature scaling
    print("Scaling on given values")
    scaled_cutdown <- scale(t(cutdown), center_scale$center,
                            center_scale$scale)
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
#' columns) and sample level metadata and returns a list object with a R::prcomp calculated
#' object, the metadata, the percent variance explained for each principal component, and
#' the genes (features) chosen for the PCA
#'
#' @param feature_by_sample raw feature (gene) count matrix (where genes/features are rows and samples are
#' columns)
#' @param meta metadata for the samples. The rows must match the columns
#' for `feature_by_sample`
#' @param ntop Number of genes/features to use in the prcomp PCA. Defaults to 1000.
#' @param HVG_selection either "classic" or "scran" to select the "ntop" features. "classic"
#' will simply the top n features by variance and "scran" will use the scran package's
#' strategy of scaling variance by expression (as highly expressed features/genes) will
#' also have higher variance and thus may be less useful for sample distinction
#' @param feature_scale Default is TRUE, which means features (genes) are center and scaled
#' with the R::scale function.
#' @param sample_cpm_scale Default is TRUE; does cpm scaling on the samples with the edgeR::cpm function
#' @param log1p Default is TRUE; log1p scales the input count matrix
#' @param remove_mito_rpl Default is TRUE, which means genes with matching the regex
#' '^MT|^RPS|^RPL' are removed.
#' @return named list object with the prcomp output returned under the $PCA slot, the given
#' metadata under the $meta slot, the percent variance of each PC as the
#' $percentVar slot, a list object containing the scaled data's "center" and "scale"
#' values for use in the metamoRph function, and the used parameters under the $param slot
#' @export
run_pca <- function(feature_by_sample,
                    meta,
                    ntop = 1000,
                    HVG_selection = 'scran',
                    feature_scale = TRUE,
                    sample_cpm_scale = TRUE,
                    log1p = TRUE,
                    remove_mito_rpl = TRUE){
  if (remove_mito_rpl){
    # remove mito and rpl/rps genes
    feature_by_sample <- feature_by_sample[!grepl('^MT|^RPS|^RPL', row.names(feature_by_sample)),]
  }
  # remove zero count features / genes
  rsums <- rowSums(feature_by_sample)
  feature_by_sample <- feature_by_sample[rsums>0,]

  if (sample_cpm_scale){
    feature_by_sample <- edgeR::cpm(feature_by_sample, log = FALSE)
  }

  if (log1p){
    feature_by_sample <- feature_by_sample |> log1p()
  }

  if (HVG_selection == 'classic'){
    Pvars <- matrixStats::rowVars(feature_by_sample)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                          length(Pvars)))]
  } else if (HVG_selection == 'scran'){
    sce_internal <- SingleCellExperiment::SingleCellExperiment(list(logcounts = (feature_by_sample)))
    feature_var <- scran::modelGeneVar(sce_internal)
    select <- feature_var |> as_tibble(rownames = 'Gene') |>
      arrange(-bio) |> head(ntop) |> pull(Gene)
  }
  # rotate
  sample_by_feature <- t(feature_by_sample[select,])
  if (feature_scale){
    print("Feature scaling")
    # run PCA
    PCA <- prcomp(sample_by_feature, scale = TRUE, center = TRUE)
    center_scale <- list(center = PCA$center,
                         scale = PCA$scale)
  } else{
    PCA <- prcomp(sample_by_feature, scale = FALSE, center = FALSE)
    center_scale = NULL
  }

  # calculate percent variance explained by each PC
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

  out <- list(PCA = PCA,
              meta = meta,
              percentVar = percentVar,
              center_scale = center_scale,
              params = list('ntop' = ntop,
                            'HVG_selection' = HVG_selection,
                            'feature_scale' = feature_scale,
                            'sample_cpm_scale' = sample_cpm_scale,
                            'log1p' = log1p,
                            'remove_mito_rpl' = remove_mito_rpl
              ))
  out
}








