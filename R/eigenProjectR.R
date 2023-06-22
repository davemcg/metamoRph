#' eigenProjectR
#'
#' This function takes in a count matrix (where genes (features) are rows and samples are
#' columns) as well as a named vector with the eigenvalues (see [eigenProjectR::run_pca()])
#' and pulls the gene (feature) information from the rotation vector and cuts down
#' the new_counts matrix to match the rotation vector gene (feature) names. Any
#' genes (features) missing from the input new_counts matrix will be filled with
#' zeros.
#'
#' The function will scale the new_counts matrix in the same manner as [eigenProjectR::run_pca()]
#' and matrix multiply by the rotation vector. The matrix output is equivalent
#' to prcomp's "$x" matrix.
#'
#' @param new_counts raw gene count matrix (where genes are rows and samples are
#' columns)
#' @param rotation matrix where the row names are genes and the col names are the
#' principal components. If you used eigenProjectR::run_pca() then this would be
#' in the output$PCA$rotation slot.
#' @param rotation_mean_sd The gene (feature) mean and sd
#' stats requried to scale the new data. If you used eigenProjectR::run_pca()
#' then this would be in the output$feature_mean_sd slot.
#' @return The new_counts matrix with the eigen transformation
#' @export
eigenProjectR <- function(new_counts, rotation, scaled_data){

  # cpm norm, then log1p
  new_cpm <- edgeR::cpm(new_counts, log = FALSE) %>% log1p()
  # extract gene names from new data
  ## remove .digit endings (if they are on ensgene then trouble if the
  ## new ensgene are a different version than the old one)
  row_genes <- row.names(new_cpm) %>% gsub('\\.\\d+','',.) %>% toupper()
  row.names(new_cpm) <- row_genes
  # # special handling for the McGaughey EiaD resources which use
  # # the gene naming scheme is: GENE (ENSG)
  # if (grepl(' \\(ENS', row.names(new_counts)[1])){
  suppressWarnings(feature_id_table <- row.names(rotation) %>% enframe() %>%
                     separate(value, c("feature_id", "ensgene"), sep = " \\(") %>%
                     mutate(ensgene = gsub(")", "", ensgene)) %>% dplyr::select(-name))

  # the internal scEiaD precalculated PCA rotations have a naming scheme as follows:
  ## GENE ID (ENSEMBL ID)
  ## for example: RHO (ENSG00000163914)
  overlap_with_ID <- row_genes[row_genes %in% feature_id_table$feature_id]
  overlap_with_ens <- row_genes[row_genes %in% feature_id_table$ensgene]

  if ((length(overlap_with_ID) < 200) &
      (length(overlap_with_ens) < 200)){
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


  # Only genes that match the rotation/eignevalues genes used
  feature_universe <- feature_id_table %>% pull(column_val) %>% make.unique()
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
  cutdown <- cutdown[feature_universe, , drop = FALSE] %>% log1p()

  # custom feature scaling
  # scaled_cutdown <- scale(t(cutdown), attr(scaled_data, "scaled:center"),
  #                        attr(scaled_data, "scaled:scale"))
  scaled_cutdown <- scale(t(cutdown))
  # Replace NAs with 0
  scaled_cutdown[is.na(scaled_cutdown)] <- 0

  # project new data onto the PCA space
  ## matrix multiply the expression of the scaled new data against
  ## the supplied eigenvector
  projected_PC_space <- scaled_cutdown %*% rotation %>% as.data.frame()

  return(projected_PC_space)
}


#' run_pca
#'
#' This function takes in a count matrix (where genes (features) are rows and samples are
#' columns) and sample level metadata and returns a list object with a R::prcomp calculated
#' object, the metadata, the percent variance explained for each principal component, and
#' the genes (features) chosen for the PCA
#'
#' @param feature_by_sample raw gene count matrix (where genes/features are rows and samples are
#' columns)
#' @param meta metadata for the samples. The rows must match the columns
#' for `feature_by_sample`
#' @param ntop Number of genes/features to use in the prcomp PCA. Defaults to 1000.
#' @param HVG_selection either "classic" or "scran" to select the "ntop" features. "classic"
#' will simply the top n features by variance and "scran" will use the scran package's
#' strategy of scaling variance by expression (as highly expressed features/genes) will
#' also have higher variance and thus may be less useful for sample distinction
#' @return named list object
#' @export
run_pca <- function(feature_by_sample,
                    meta,
                    ntop = 1000,
                    HVG_selection = 'scran'){
  feature_by_sample_original <- feature_by_sample

  # remove mito and rpl/rps genes
  feature_by_sample <- feature_by_sample[!grepl('^MT|^RPS|^RPL', row.names(feature_by_sample)),]

  # remove zero count genes
  rsums <- rowSums(feature_by_sample_original)
  feature_by_sample <- feature_by_sample_original[rsums>0,]


  # cpm scale by sample (all counts for a sample sum to 1e6)
  # then log1p transform
  feature_by_sample <- edgeR::cpm(feature_by_sample, log = FALSE) %>% log1p()

  if (HVG_selection == 'classic'){
    Pvars <- rowVars(feature_by_sample)
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                          length(Pvars)))]
  } else if (HVG_selection == 'scran'){
    sce_internal <- SingleCellExperiment::SingleCellExperiment(list(logcounts = (feature_by_sample)))
    feature_var <- scran::modelGeneVar(sce_internal)
    #select <- getTopHVGs(feature_var, n = ntop)
    select <- feature_var %>% as_tibble(rownames = 'Gene') %>%
      arrange(-bio) %>% head(ntop) %>% pull(Gene)
  }
  # log and scale (by gene)
  sample_by_feature <- scale(t((feature_by_sample[select,])))
  # run PCA
  PCA <- prcomp(sample_by_feature, scale = FALSE, center = FALSE)
  # grab HVG names
  HVG <- colnames(sample_by_feature)
  # calculate percent variance explained by each PC
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

  out <- list(PCA = PCA,
              meta = meta,
              percentVar = percentVar,
              HVG = HVG,
              scaled_select = sample_by_feature)
  out
}
############








