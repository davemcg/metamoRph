# build test metamoRph data from EiaD
library(dplyr)
library(readr)
mat_all <- data.table::fread('~/git/EiaD_build/data/EiaD_pca_analysis_mat_all.csv.gz')
genes <- mat_all$Gene
mat_all <- mat_all[,-1] %>% as.matrix()
row.names(mat_all) <- genes
meta_mat_all <- data.table::fread('~/git/EiaD_build/data/EiaD_pca_analysis_meta_mat_all.csv.gz')

set.seed(984)
meta_ref <- meta_mat_all %>% filter(!grepl("AMD", Perturbation),
                                    study_accession != 'SRP310948',
                                    Age == 'Adult',
                                    study_accession != 'SRP097696', # unusual endothelial retina set
                                    Tissue %in% c("Conjunctiva", "Cornea","Lens", "RPE","Retina", "Trabecular Meshwork")) %>%
  group_by(Tissue) %>%
  sample_n(20,replace = TRUE) %>% unique()

meta_project <- meta_mat_all %>% filter(study_accession == 'SRP310948')

out_counts <- bind_cols(mat_all[,meta_ref %>% pull(sample_accession)]%>%
                          as_tibble(rownames = 'Gene'),
                        mat_all[, meta_project %>% pull(sample_accession)])

Pvars <- rowVars(as.matrix(out_counts[,2:ncol(out_counts)] ))
names(Pvars) <- out_counts$Gene
select <- sort(Pvars, decreasing = TRUE)[seq_len(min(2000, length(Pvars)))] |> names()


write_csv(out_counts |> filter(Gene %in% select),
          '~/git/metamoRph/inst/test_data//EiaD__eye_samples_counts.csv.gz')

write_csv(meta_project,
          '~/git/metamoRph/inst/test_data/EiaD__eye_samples_metaProject.csv.gz')

write_csv(meta_ref,
          '~/git/metamoRph/inst/test_data/EiaD__eye_samples_metaRef.csv.gz')
