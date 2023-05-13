#' model_build
#'
#' This function uses the training_data and the training_labels to build a lm
#' based model for each label type which can be used in [eigenProjectR::model_apply()].
#' The training_data is intended to be the sample x PC (principal component)
#' row x column matrix. Which is the `$x` output from base R prcomp. We provide
#' precomputed prcomp PCA outputs from the \url{plae.nei.nih.gov} resource
#' for adult human eye, adult mouse eye, fetal human eye, and fetal mouse eye (
#' see \code{vignette("pca_download", package = "eigenProjectR")})
#'
#' @param training_data sample (row) by principal component (column) matrix
#' @param training_labels vector which has the row-matched labels (e.g. cell types)
#' for each sample.
#' @param num_PCs number of principal components to use from the training_data.
#' Defaults to the first (top) 200.
#' @param BPPARAM The BiocParallel class
#' @param model lm or glm. We strongly recommend using lm.
#' @param verbose Print training status for each label type
#' @return A list of models for each individual label type
#' @export
model_build <- function(training_data,
                        training_labels,
                        num_PCs = 200,
                        BPPARAM = BiocParallel::SerialParam(),
                        model = 'lm',
                        verbose = TRUE){
  # cut down to requested PC
  ## issue warning if num_PCs > number of actual PC
  if (ncol(training_data) < num_PCs){
    warning(paste0("Dropping num PCs requested from ",
                   num_PCs, " to number of columns in training data"))
  }
  num_PCs <- min(ncol(training_data), num_PCs)
  training_data <- training_data[,1:num_PCs]
  require(BiocParallel)
  models <- lapply(unique(sort(training_labels)), function(target) {
    if (verbose){
      print(paste0("Model training for ", target))
    }
    binary_training_labels <- training_labels
    binary_training_labels[binary_training_labels != target] <- 'Other'

    labels <- binary_training_labels %>%
      factor(levels = c(target, "Other")) %>%
      as.integer() - 1
    # if (model == 'xgboost'){
    # nrounds = 200
    # early_stopping_rounds = 10
    # subsample = 0.5
    # xg_verbose = FALSE
    # eta = 0.3
    #   model <- xgboost(data = training_data,
    #                    label = labels,
    #                    nrounds = nrounds,
    #                    early_stopping_rounds = early_stopping_rounds,
    #                    subsample = subsample,
    #                    eta = eta,
    #                    objective = 'binary:logistic',
    #                    verbose = xg_verbose)
    # }
    if (model == 'lm'){
      model <- lm(labels ~ ., data.frame(cbind(labels,training_data)))
    } else if (model == 'glm'){
      mdata <- data.frame(cbind(labels,training_data))
      mdata$labels <- as.factor(mdata$labels)
      model <- glm(labels ~ ., data = mdata, family = binomial(link='log'))
    }
    return(list(target = target, model = model))
  })
  model_out <- list()
  for (i in 1:length(models)){
    model_out[[models[[i]]$target]] <- models[[i]]$model
  }
  model_out
}


#' model_apply
#'
#' This function uses the output from [eigenProjectR::model_build()]
#'
#' @param list_of_models sample (row) by principal component (column) matrix
#' @param experiment_data Projected data [eigenProjectR::eigenProjectR()]
#' @param experiment_labels Optional labels for the users experiment_data
#' @param return_predictions By default the predicted labels are returned.
#' If set to TRUE, then the entire matrix of probabilites for each sample (against
#' each label type) is returned.
#' @return By default, a table of predicted labels for your data
#' @export
model_apply <- function(list_of_models,
                        experiment_data,
                        experiment_labels = '',
                        return_predictions = FALSE){

  predictions <- list()
  for (i in names(list_of_models)){
    predictions[[i]] <- (predict(list_of_models[[i]], (experiment_data)))
  }
  # turn predictions into a tibble
  predictions <- predictions %>% bind_cols()
  colnames(predictions) <- names(list_of_models)
  # identify the most likely celltype call (lowest value) for each sample
  calls <- cbind(row.names(experiment_data),
                 colnames(predictions)[apply(predictions, 1, which.min)],
                 apply(predictions, 1, min)) %>%
    data.frame()
  colnames(calls) <- c('sample_id', 'predict', 'min_score')
  # optionally put in true labels for input data (if given)
  if (experiment_labels[1] != ''){
    calls$sample_label <- experiment_labels
    colnames(calls) <- c('sample_id', 'predict', 'min_score', 'sample_label')
    calls <- calls %>% relocate(sample_id, sample_label)
  }

  calls <- as_tibble(calls) %>%
    mutate(min_score = as.numeric(min_score),
           predict_stringent = case_when(min_score >0.5 ~ 'Unknown',
                                         TRUE ~ predict)) %>%
    relocate(min_score, .after = predict_stringent)
  if (return_predictions){
    predictions
  } else{
    calls
  }
}


