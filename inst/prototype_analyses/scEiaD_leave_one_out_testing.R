#remotes::install_github("davemcg/eigenProjectR")
library(metamoRph)
library(tidyverse)
library(data.table)
library(matrixStats)

load('inst/data/scEiaD_seacell_pca_objs_2023_07_10.Rdata')
load('inst/data/scEiaD_seacell_mat_2023_07_10.Rdata')
# create list to store label guesses for each age group
label_guesses_list <- list()
label_val_list <- list()

# loop through excluded studies
for (study in unique(meta_mat_human_adult$study_accession)) {
  print(study)
  # filter to exclude one study
  meta_filtered_excluded <- meta_mat_human_adult %>% filter(study_accession != study)

  # line up gene counts to metadata
  gene_counts_filtered <- seacell_mat[meta_filtered_excluded$seacell_id,]

  # run pca
  pca_output <- run_pca(t(gene_counts_filtered %>% as.matrix()),
                        meta_filtered_excluded)

  # new counts to be labeled
  new_counts <- seacell_mat[meta_mat_human_adult %>% filter(study_accession == study) %>%
                              pull(seacell_id), , drop = FALSE] %>% as.matrix()

  # apply to new data (excluded study)
  colnames(new_counts) <- gsub(' \\(.*','',colnames(new_counts))
  projected_new_data <- metamoRph(new_counts = t(new_counts),
                                  rotation = pca_output$PCA$rotation,
                                  center_scale = pca_output$center_scale)

  trained_model <- model_build(pca_output$PCA$x[,1:200],
                               pca_output$meta$CellType,
                               model = 'lm',
                               BPPARAM = MulticoreParam(10))

  label_guesses <- model_apply(trained_model,
                               projected_new_data,
                               meta_mat_human_adult %>% filter(study_accession == study) %>%
                                 pull(CellType)

  )

  #label_guesses %>% group_by(sample_label, predict) %>% summarise(Count = n()) %>% data.frame
  #label_guesses %>% group_by(sample_label, predict) %>% summarise(Count = n()) %>% data.frame %>% filter(sample_label == predict) %>% pull(Count) %>% sum()
  label_mat <- model_apply(trained_model, projected_new_data, meta_mat_human_adult %>% filter(study_accession == study) %>% pull(CellType),
                           return_predictions = TRUE)

  #label_guesses_list[[paste0(age_group, "_", study)]] <- label_guesses
  label_guesses_list[[study]] <- label_guesses
  label_val_list[[study]] <- label_mat
  #cat("***", "Study", study, "has been processed for the", age_group, "age group", "***", "\n")
}

# combine label guesses into a single data frame
predictions <- bind_rows(label_guesses_list, .id = "study_accession")
pred_values <- bind_rows(label_val_list, .id = "study_accession")

write_tsv(predictions, 'data/2023_07_12_ML_tissue_predictions.tsv.gz')
write_tsv(pred_values, 'data/2023_07_12_ML_tissue_values.tsv.gz')
