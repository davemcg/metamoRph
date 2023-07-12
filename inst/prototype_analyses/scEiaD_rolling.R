
library(metamoRph)
library(tidyverse)
library(data.table)
library(matrixStats)

load('inst/data/scEiaD_seacell_pca_objs_2023_07_07.Rdata')

sample_frac <- 0.1
sample_min <- 10
rounds <- 5
label_guesses_list <- list()
label_val_list <- list()

training_mat <- t(seacell_mat)
training_metadata <- meta_mat_human_adult %>% filter(CellType != 'NA')

projection_mat <- t(seacell_mat)
projection_metadata <- meta_mat_human_adult
# SCP1755 example (see SCP1755_broad_attempt for loading data)
#projection_mat <- scp1755_sum_counts
#projection_metadata <- colnames(scp1755_sum_counts) %>% enframe()

# e.g. "CellType" or "Tissue"
sample_label <- 'CellType'
projection_id <- 'CellType'
# SCP1755 example (see SCP1755_broad_attempt for loading data)
#projection_id <- 'value'

# column that contains the UNIQUE sample identifier in
# like SRS12341234 or "BrainSample1"
training_sample_id <- 'seacell_id'
#projection_sample_id <- 'seacell_id'
# SCP1755 example (see SCP1755_broad_attempt for loading data)
#projection_sample_id <- 'value'


na_if_sample_in_training <- FALSE
for (bin in seq(1:rounds)) {
  #set.seed(sample(10))
  print(paste(bin, ' of ', rounds))
  training_metadata$sample_label <- training_metadata %>% pull(sample_label)
  # run twice, first by sample_frac
  meta_for_training_by_frac <- training_metadata %>%
    group_by(sample_label) %>%
    slice_sample(prop = sample_frac) %>%
    mutate(SCount = n()) %>%
    unique()
  # second time by sample_n
  meta_for_training_by_n <- training_metadata %>%
    group_by(sample_label) %>%
    slice_sample(n = sample_min) %>%
    mutate(SCount = n()) %>%
    unique()
  # combine, remove any meta_for_training_by_frac samples with
  # a Scount < sample_min and replace with the meta_for_training_by_n
  meta_for_training <- meta_for_training_by_frac %>% filter(SCount > sample_min)
  meta_for_training <- bind_rows(meta_for_training,
                                 meta_for_training_by_n %>%
                                   filter(!sample_label %in% (meta_for_training %>% pull(sample_label))))
  # line up gene counts to metadata
  gene_counts_train <- training_mat[,meta_for_training %>% pull(training_sample_id)]

  # run pca
  pca_output <- run_pca(gene_counts_train, meta_for_training)
  # train the model
  trained_model <- model_build(pca_output$PCA$x,
                               pca_output$meta %>% pull(sample_label),
                               model = 'lm', verbose = FALSE, 150)

  # project data from the pca
  rownames(projection_mat) <- gsub(' \\(.*','',rownames(projection_mat))
  projected_data <- metamoRph((projection_mat[,projection_metadata %>% pull(projection_sample_id)]),
                                  pca_output$PCA$rotation,
                                  pca_output$center_scale)
  # apply model
  label_guesses <- model_apply(trained_model,
                               projected_data,
                               projection_metadata %>% pull(projection_id)

  )
  label_mat <- model_apply(trained_model,
                           projected_data,
                           projection_metadata %>% pull(projection_id),
                           return_predictions = TRUE)
  if (na_if_sample_in_training){
    label_mat[(projection_metadata[training_sample_id] %>% pull(1)) %in%
                (meta_for_training[training_sample_id] %>% pull(1)),] <- NA
  }

  label_guesses_list[[bin]] <- label_guesses
  label_val_list[[bin]] <- label_mat
}

# combine label guesses into a single data frame
predictions <- label_guesses_list %>% bind_rows() %>%
  group_by(sample_id) %>% summarise(miscall_count = sum(sample_label != predict),
                                    sample_label = unique(sample_label),
                                    predict = paste(unique(predict), collapse = ', '),
                                    mean_min_score = mean(min_score),
                                    min = min(min_score),
                                    max = max(min_score))
#pred_values <- bind_rows(label_val_list, .id = "study_accession")

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
               apply(model_scores, 1, min),
               model_scores) %>%
  data.frame()
colnames(calls)[1:4] <- c("sample_id","sample_label","predict","mean_min_score")
calls$mean_min_score <- as.numeric(calls$mean_min_score)
calls %>% filter(sample_label == predict) %>% dim()
calls %>% filter(sample_label != predict) %>% dim()

predictions <- predictions %>% dplyr::select(-predict) %>%
  left_join(calls %>% dplyr::select(sample_id, predict))

predictions %>% data.frame() %>% dplyr::select(sample_id, predict, mean_min_score)
