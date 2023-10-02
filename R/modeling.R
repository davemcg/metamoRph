#' model_build
#'
#' This function uses the training_data and the training_labels to build a lm
#' based model for each label type which can be used in [metamoRph::model_apply()].
#' The training_data is intended to be the sample x PC (principal component)
#' row x column matrix. Which is the `$x` output from base R prcomp. We provide
#' precomputed prcomp PCA outputs from the \url{plae.nei.nih.gov} resource
#' for adult human eye, adult mouse eye, fetal human eye, and fetal mouse eye (
#' see \code{vignette("pca_download", package = "metamoRph")})
#'
#' @param training_data sample (row) by principal component (column) matrix
#' @param training_labels vector which has the row-matched labels (e.g. cell types)
#' for each sample.
#' @param num_PCs number of principal components to use from the training_data.
#' Defaults to the first (top) 50.
#' @param BPPARAM The BiocParallel class
#' @param model Default is lm. We also support xgboost, glm, rf, and svm.
#' In our experience we find lm and svm to be the best performers.
#' @param verbose Print training status for each label type
#' @return A list of models for each individual label type
#' @import BiocParallel
#' @importFrom stats lm glm binomial
#' @export
model_build <- function(training_data,
                        training_labels,
                        num_PCs = 50,
                        BPPARAM = BiocParallel::SerialParam(),
                        model = 'lm',
                        verbose = TRUE){
  # ensure the training labels do not start with a digit, as the 
  # model_apply relies on column names
  # to align the data ...which cannot start with a digit
  training_labels = paste0('metamoRphPrefix_', training_labels)
  # cut down to requested PC
  ## issue warning if num_PCs > number of actual PC
  if (ncol(training_data) < num_PCs){
    warning(paste0("Dropping num PCs requested from ",
                   num_PCs, " to number of columns in training data (", ncol(training_data), ")")
    )
    num_PCs <- ncol(training_data)
  }
  training_data <- training_data[,1:num_PCs]

  models <- bplapply(unique(sort(training_labels)), function(target) {
    if (verbose){
      message(paste0("Model training for ", gsub('metamoRphPrefix_','',target)))
    }
    binary_training_labels <- training_labels
    binary_training_labels[binary_training_labels != target] <- 'Other'

    labels <- binary_training_labels %>%
      factor(levels = c("Other", target)) %>%
      as.integer() - 1

    if (model == 'lm'){
      model <- lm(labels ~ ., data.frame(cbind(labels,training_data)))
    } else if (model == 'glm'){
      model <- glm(labels ~ ., data = data.frame(cbind(labels, training_data)),family = binomial("logit"))
    } else if (model == 'xgboost'){
      if(!requireNamespace("parsnip")){
        message("This requires the 'parsnip' package.")
        return(invisible())
      }
      model <- parsnip::boost_tree(mode = "regression", trees = 200) %>%
        parsnip::fit(labels ~ .,
            data = data.frame(cbind(labels,training_data)))
      # } else if (model == 'nnet'){
      #   model <- mlp() %>%
      #     set_engine("nnet") %>%
      #     set_mode("regression") %>%
      #     fit(labels ~ .,
      #         data = data.frame(cbind(labels,training_data)))
    } else if (model == 'rf'){
      if(!requireNamespace("parsnip")){
        message("This requires the 'parsnip' package.")
        return(invisible())
      }
      model <-  parsnip::rand_forest(trees = 500, min_n = 5) %>%
        parsnip::set_mode("regression") %>%
        parsnip::set_engine("ranger", importance = "impurity") %>%
        parsnip::fit(labels ~ .,
            data = data.frame(cbind(labels,training_data)))
    } else if (model == 'svm'){
      if(!requireNamespace("parsnip")){
        message("This requires the 'parsnip' package.")
        return(invisible())
      }
      model <- parsnip::svm_linear(cost = 1, margin = 0.1) %>%
        parsnip::set_mode("regression") %>%
        parsnip::set_engine("LiblineaR") %>%
        parsnip::fit(labels ~ .,
            data = data.frame(cbind(labels,training_data)))
    }
    return(list(target = target, model = model))
  })
  model_out <- list()
  for (i in 1:length(models)){
    model_out[[models[[i]]$target]] <- models[[i]]$model
  }
  
  # remove prefix
  names(model_out) <- gsub("metamoRphPrefix_","",names(model_out))
  return(model_out)
}


#' model_apply
#'
#' This function uses the output from [metamoRph::model_build()]
#'
#' @param list_of_models list object containing one model per sample type (e.g. photoreceptors vs not-photoreceptors)
#' @param experiment_data Projected data [metamoRph::metamoRph()]
#' @param experiment_labels Optional labels for the users experiment_data
#' @param return_predictions By default the predicted labels are returned.
#' If set to TRUE, then the entire matrix of probabilities for each sample (against
#' each label type) is returned.
#' @return By default, a table of predicted labels for your data, if return_predictions set to TRUE then
#' a matrix of probabilities for each sample type is returned (where rows are samples and columns
#' are each sample type probability)
#' @importFrom stats predict
#' @importFrom dplyr bind_cols as_tibble mutate relocate case_when
#' @export
model_apply <- function(list_of_models,
                        experiment_data,
                        experiment_labels = '',
                        return_predictions = FALSE){
  
  sample_id <- sample_label <- max_score <- predict_stringent <- NULL
  predictions <- list()
  for (i in names(list_of_models)){
    predictions[[i]] <- (predict(list_of_models[[i]], (experiment_data)))
  }
  # turn predictions into a tibble
  predictions <- predictions %>% bind_cols()
  colnames(predictions) <- names(list_of_models)

  # identify the most likely celltype call (highest value) for each sample
  # and the 2nd
  calls <- cbind(row.names(experiment_data),
                 colnames(predictions)[apply(predictions, 1, which.max)],
                 apply(predictions, 1, max)) %>%
    data.frame()
  colnames(calls) <- c('sample_id', 'predict', 'max_score')
  # optionally put in true labels for input data (if given)
  if (experiment_labels[1] != ''){
    calls$sample_label <- experiment_labels
    colnames(calls) <- c('sample_id', 'predict', 'max_score', 'sample_label')
    calls <- calls %>% relocate(sample_id, sample_label)
  }

  calls <- as_tibble(calls) %>%
    mutate(max_score = as.numeric(max_score),
           predict_stringent = case_when(max_score < 0.5 ~ 'Unknown',
                                         TRUE ~ predict)) %>%
    relocate(max_score, .after = predict_stringent)


  if (return_predictions){
    # remove prefix before outputting results
    colnames(predictions) <- gsub('metamoRphPrefix_','',colnames(predictions))
    predictions
  } else{

    # remove prefix before outputting results
    calls <- calls %>% mutate(predict = gsub('metamoRphPrefix_','',predict), 
                                          predict_stringent = gsub('metamoRphPrefix_','',predict_stringent))
    calls
  }
}


