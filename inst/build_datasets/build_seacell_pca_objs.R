library(matrixStats)
library(tidyverse)
#library(umap)

# URL is hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/seacells/seacells_2023_06_15.tar
# you have to manually download and extract this data from [URL]
# on the command line (bash)
# e.g. wget [URL]
# then tar -xzvf seacells_2023_06_15.tar
meta_files <- list.files('~/data/eiad_seacells/', full.names = TRUE) %>% grep('obs', ., value=TRUE)
seacell_files <- list.files('~/data/eiad_seacells/', full.names = TRUE) %>% grep('aggr', ., value=TRUE)

meta_l <- list()
for (i in meta_files){
  meta_l[[i]] <- data.table::fread(i) %>% mutate_all(as.character)
}
meta <- bind_rows(meta_l)

seacell_l <- list()
for (i in seacell_files){
  seacell_l[[i]] <- data.table::fread(i)
}
seacell <- bind_rows(seacell_l)

seacell_summarise <- meta %>%
  mutate(CellType_predict = case_when(!is.na(TabulaMurisCellType_predict) ~ TabulaMurisCellType_predict,
                                      TRUE ~ CellType_predict)) %>%
  mutate(seacell_id = id, UMAP_1 = as.numeric(UMAP_1), UMAP_2 = as.numeric(UMAP_2)) %>%
  group_by(seacell_id, Compartment, Tissue, Source, study_accession, organism, integration_group, CellType_predict) %>%
  summarise(Count = n(), UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2)) %>%
  mutate(Ratio = Count / sum(Count)) %>%
  filter(Ratio > 0.15) %>%
  arrange(-Ratio) %>%
  mutate(CT = paste0(CellType_predict, ' (', round(Ratio, 2), ')')) %>%
  summarise(CTP = paste0(CT, collapse =', '),  UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))


seacell_summarise <- seacell_summarise %>%
  left_join(meta %>%
              mutate(seacell_id =id) %>%
              filter(!is.na(Paper)) %>%
              group_by(seacell_id) %>%
              summarise(Paper = paste(unique(Paper), collapse = ', '))) %>%
  mutate(CellType = gsub(' \\(.*','',CTP))
meta %>%
  mutate(seacell_id = id) %>%
  dplyr::select(seacell_id, organism, integration_group) %>% unique()

seacell_summarise %>% ungroup() %>% sample_n(10)


seacell_mat <- seacell[,2:ncol(seacell)] %>% as.matrix()
row.names(seacell_mat) <- seacell[,1] %>% pull(1)


#####################################
# add in human readable gene names
# gene_annotations <-
#   rtracklayer::import("data/human.gene_sums.G029.gtf.gz") %>%
#   as.data.frame()
library('biomaRt')
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
gene_table <- getBM(
  mart = mart,
  attributes = c( 'ensembl_gene_id', 'hgnc_symbol'),
  filter = 'ensembl_gene_id',
  values = colnames(seacell_mat),
  uniqueRows = TRUE)

# new gene names
colnames(seacell_mat) <- colnames(seacell_mat) %>%
  enframe() %>%
  left_join(gene_table %>%
              group_by(ensembl_gene_id) %>%
              # if there are multiple ensembl <-> hgnc relations, then
              # just pick the first one. YOLO.
              summarise(hgnc_symbol = head(hgnc_symbol,1)) %>%
              mutate(value = ensembl_gene_id) ) %>%
  mutate(nname = paste0(hgnc_symbol, ' (', value, ')')) %>%
  pull(nname)

# add counts info to meta
seacell_summarise <- seacell_summarise %>% left_join(seacell_mat %>% rowSums() %>% enframe(name = 'seacell_id', value = 'sum_counts'))


library(metamoRph)

##### ensure metadata is in same order as count matrix
meta_mat <- row.names(seacell_mat) %>% enframe(value = 'seacell_id') %>% left_join(seacell_summarise)
meta_mat <- meta_mat %>% ungroup() %>%
  mutate(Compartment = case_when(Tissue == 'Sclera' ~ 'Front Eye',
                                 TRUE ~ Compartment))

# human adult pca
set.seed(20230710)
meta_mat_human_adult <- meta_mat %>%
  filter(organism == 'Homo sapiens', Compartment != 'Body',
         study_accession != 'SRP238587',
         integration_group %in% c("Mature","Maturing"),
         Source == 'Tissue', Compartment == 'Back Eye',
         !grepl("RPC", CellType)) %>%
  group_by(CellType) %>% sample_n(50, replace = TRUE) %>%
  unique()

seacell_pca_human_adult <- run_pca(t(seacell_mat[meta_mat_human_adult$seacell_id,]),
                                   meta_mat_human_adult,
                                   hvg_selection = 'scran')

# human fetal pca
meta_mat_human_fetal <- meta_mat %>% filter(organism == 'Homo sapiens',
                                            !integration_group %in% c("Mature","Maturing"),
                                            Source == 'Tissue') %>%
  group_by(CellType) %>% sample_n(50, replace = TRUE) %>%
  unique()

seacell_pca_human_fetal <- run_pca(t(seacell_mat[meta_mat_human_fetal$seacell_id,]) ,
                                   meta_mat_human_fetal,
                                   hvg_selection = 'scran')

# mouse adult pca
meta_mat_mouse_adult <- meta_mat %>%
  filter(organism == 'Mus musculus',
         study_accession != 'SRP238587',
         integration_group %in% c("Mature","Maturing"), Source == 'Tissue',
         !grepl("RPC|Neuro|Precu|NA", CellType)) %>%
  group_by(CellType) %>% sample_n(50, replace = TRUE) %>%
  unique()
seacell_pca_mouse_adult <- run_pca(t(seacell_mat[meta_mat_mouse_adult$seacell_id,]),
                                   meta_mat_mouse_adult,
                                   hvg_selection = 'scran')

# mouse fetal pca
meta_mat_mouse_fetal <- bind_rows(meta_mat %>% filter(organism == 'Mus musculus',  !integration_group %in% c("Mature","Maturing"), Source == 'Tissue'),
                                  meta_mat %>% filter(organism == 'Mus musculus',  grepl("RPC|Neuro|Precu", CellType, Source == 'Tissue'))) %>%
  unique() %>%
  group_by(CellType) %>% sample_n(50, replace = TRUE) %>%
  unique()
seacell_pca_mouse_fetal <- run_pca(t(seacell_mat[meta_mat_mouse_fetal$seacell_id,]),
                                   meta_mat_mouse_fetal,
                                   hvg_selection = 'scran')

# # all
# meta_mat_all_downsampled <- meta_mat %>% group_by(organism, CellType) %>% sample_n(50, replace = TRUE) %>% unique()
# seacell_pca <- run_pca(t(seacell_mat[meta_mat_all_downsampled$seacell_id,]),
#                        meta_mat_all_downsampled,
#                        hvg_selection = 'scran',
#                        ntop = 1000)


save(
     meta_mat,
     meta_mat_human_adult,
     meta_mat_human_fetal,
     meta_mat_mouse_adult,
     meta_mat_mouse_fetal,
    # meta_mat_all_downsampled,
     seacell_pca_human_adult,
     seacell_pca_human_fetal,
     seacell_pca_mouse_adult,
     seacell_pca_mouse_fetal,
     # seacell_pca,
    file = 'inst/data/scEiaD_seacell_pca_objs_2023_07_10.Rdata')
save(meta_mat, seacell_mat, file = 'inst/data/scEiaD_seacell_mat_2023_07_10.Rdata')

