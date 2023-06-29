
library(eigenProjectR)
library(tidyverse)
library(data.table)
library(matrixStats)


load('inst/data/scEiaD_seacell_pca_objs_2023_06_24.Rdata')

sample_frac <- 0.1
sample_min <- 5
rounds <- 20
label_guesses_list <- list()
label_val_list <- list()

training_mat <- gene_counts
training_metadata <- meta

projection_mat <- gene_counts
projection_metadata <- meta

# e.g. "CellType" or "Tissue"
sample_label <- 'Tissue'

training_sample_id <- 'sample_accession'
# column that contains the UNIQUE sample identifier in
# like SRS12341234 or "BrainSample1"
projection_sample_id <- 'sample_accession'


na_if_sample_in_training <- FALSE


for (bin in seq(1:rounds)) {
  set.seed(sample(10))
  print(paste(bin, ' of ', rounds))
  training_metadata$sample_label <- training_metadata %>% pull(sample_label)
  meta_for_training <- training_metadata %>%
    group_by(sample_label) %>%
    slice_sample(n = sample_min) %>%
    unique()

  # line up gene counts to metadata
  gene_counts_train <- training_mat[,meta_for_training %>% pull(training_sample_id)]

  # run pca
  pca_output <- run_pca((gene_counts_train %>% as.matrix()), meta_for_training)
  # train the model
  trained_model <- model_build(pca_output$PCA$x,
                               pca_output$meta %>% pull(sample_label),
                               model = 'lm', verbose = FALSE)

  # project data from the pca
  rownames(projection_mat) <- gsub(' \\(.*','',rownames(projection_mat))
  projected_data <- eigenProjectR((projection_mat[,projection_metadata %>% pull(projection_sample_id)]),
                                  pca_output$PCA$rotation)
  # apply model
  label_guesses <- model_apply(trained_model,
                               projected_data,
                               projection_metadata %>% pull(sample_label)

  )

  #label_guesses %>% group_by(sample_label, predict) %>% summarise(Count = n()) %>% data.frame
  #label_guesses %>% group_by(sample_label, predict) %>% summarise(Count = n()) %>% data.frame %>% filter(sample_label == predict) %>% pull(Count) %>% sum()
  label_mat <- model_apply(trained_model,
                           projected_data,
                           projection_metadata %>% pull(sample_label),
                           return_predictions = TRUE)
  if (na_if_sample_in_training){
    label_mat[(projection_metadata[training_sample_id] %>% pull(1)) %in% (meta_for_training[training_sample_id] %>% pull(1)),] <- NA
  }

  #label_guesses_list[[paste0(age_group, "_", study)]] <- label_guesses
  label_guesses_list[[bin]] <- label_guesses
  label_val_list[[bin]] <- label_mat
  #cat("***", "Study", study, "has been processed for the", age_group, "age group", "***", "\n")
}

# combine label guesses into a single data frame
#predictions <- bind_rows(label_guesses_list, .id = "study_accession")
#pred_values <- bind_rows(label_val_list, .id = "study_accession")

model_scores <- (Reduce("+", label_val_list) / length(label_val_list))

# create 3D array to stack all of the model scores
model_score_3D <- array(unlist(label_val_list),
                        c(nrow(label_val_list[[1]]),
                          ncol(label_val_list[[1]]),
                          length(label_val_list)))
model_scores <- rowMeans(model_score_3D, dims = 2, na.rm = TRUE)
colnames(model_scores) <- colnames(label_val_list[[1]])
calls <- cbind(label_guesses_list[[1]] %>% pull(sample_id),
               label_guesses_list[[1]] %>% pull(sample_label),
               colnames(model_scores)[apply(model_scores, 1, which.min)],
               apply(model_scores, 1, min)) %>%
  data.frame()
calls %>% filter(X2 == X3) %>% dim()
calls %>% filter(X2 != X3) %>% dim()
