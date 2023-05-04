library(matrixStats)
library(tidyverse)
library(umap)

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
  mutate(seacell_id = paste0(sample_accession, '__', str_extract(SEACell, '\\d+'))) %>%
  group_by(seacell_id, Compartment, Tissue, Source, study_accession, organism, integration_group, CellType_predict) %>%
  summarise(Count = n()) %>%
  mutate(Ratio = Count / sum(Count)) %>%
  filter(Ratio > 0.15) %>%
  arrange(-Ratio) %>%
  mutate(CT = paste0(CellType_predict, ' (', round(Ratio, 2), ')')) %>%
  summarise(CTP = paste0(CT, collapse =', '))

seacell_summarise <- seacell_summarise %>%
  left_join(meta %>%
              mutate(seacell_id = paste0(sample_accession, '__', str_extract(SEACell, '\\d+'))) %>%
              filter(!is.na(Paper)) %>%
              group_by(seacell_id) %>%
              summarise(Paper = paste(unique(Paper), collapse = ', '))) %>%
  mutate(CellType = gsub(' \\(.*','',CTP))
meta %>%
  mutate(seacell_id = paste0(sample_accession, '__', str_extract(SEACell, '\\d+'))) %>%
  select(seacell_id, organism, integration_group) %>% unique()
seacell_summarise %>% head()



#####################################
# add in human readable gene names
gene_annotations <-
  rtracklayer::import("data/human.gene_sums.G029.gtf.gz") %>%
  as.data.frame()

seacell_mat <- seacell[,2:ncol(seacell)] %>% as.matrix()
row.names(seacell_mat) <- seacell[,1] %>% pull(1)

# new gene names
colnames(seacell_mat) <- colnames(seacell_mat) %>%
  enframe() %>%
  left_join(gene_annotations %>% mutate(value = gsub('\\.\\d+','',gene_id)) ) %>%
  mutate(nname = paste0(gene_name, ' (', value, ')')) %>%
  pull(nname)

# add counts info to meta
seacell_summarise <- seacell_summarise %>% left_join(seacell_mat %>% rowSums() %>% enframe(name = 'seacell_id', value = 'sum_counts'))

# write
#write_tsv(meta, file = '../data/seacell_meta.tsv.gz')
#write_tsv(seacell_summarise, file = '../data/seacell_summarise.tsv.gz')
#write_tsv(seacell_mat %>% as_tibble(rownames = 'seacell_id'), file = '../data/seacell_mat.tsv.gz')

source('src/run_pca.R')

##### ensure metadata is in same order as count matrix
meta_mat <- row.names(seacell_mat) %>% enframe(value = 'seacell_id') %>% left_join(seacell_summarise)


# human adult pca
meta_mat_human_adult <- meta_mat %>% filter(organism == 'Homo sapiens', study_accession != 'SRP238587', integration_group %in% c("Mature","Maturing"), Source == 'Tissue', !grepl("RPC|Neuro|Precu", CellType))
seacell_pca_human_adult <- run_pca(t(seacell_mat[meta_mat_human_adult$seacell_id,]),
                                  meta_mat_human_adult,
                                  HVG_selection = 'scran',
                                  ntop = 1000)

# human fetal pca
meta_mat_human_fetal <- meta_mat %>% filter(organism == 'Homo sapiens', !integration_group %in% c("Mature","Maturing"), Source == 'Tissue')
seacell_pca_human_fetal <- run_pca(t(seacell_mat[meta_mat_human_fetal$seacell_id,]) ,
                                   meta_mat_human_fetal,
                                   HVG_selection = 'scran',
                                   ntop = 1000)

# human pca
# meta_mat_human <- meta_mat %>% filter(organism == 'Homo sapiens', Source == 'Tissue')
# seacell_pca_human <- run_pca(t(seacell_mat[meta_mat_human$seacell_id,]),
#                              meta_mat_human,
#                              HVG_selection = 'scran',
#                              ntop = 1000)

# mouse adult pca
meta_mat_mouse_adult <- meta_mat %>% filter(organism == 'Mus musculus', study_accession != 'SRP238587', integration_group %in% c("Mature","Maturing"), Source == 'Tissue', !grepl("RPC|Neuro|Precu", CellType))
seacell_pca_mouse_adult <- run_pca(t(seacell_mat[meta_mat_mouse_adult$seacell_id,]),
                                   meta_mat_mouse_adult,
                                   HVG_selection = 'scran',
                                   ntop = 1000)

# mouse fetal pca
meta_mat_mouse_fetal <- bind_rows(meta_mat %>% filter(organism == 'Mus musculus',  !integration_group %in% c("Mature","Maturing"), Source == 'Tissue'),
                                  meta_mat %>% filter(organism == 'Mus musculus',  grepl("RPC|Neuro|Precu", CellType, Source == 'Tissue'))) %>%
  unique()
seacell_pca_mouse_fetal <- run_pca(t(seacell_mat[meta_mat_mouse_fetal$seacell_id,]),
                                   meta_mat_mouse_fetal,
                                   HVG_selection = 'scran',
                                   ntop = 1000)
# mouse pca
# meta_mat_mouse <- meta_mat %>% filter(organism == 'Mus musculus', Source == 'Tissue')
# seacell_pca_mouse <- run_pca(t(seacell_mat[meta_mat_mouse$seacell_id,]),
#                              meta_mat_mouse,
#                              HVG_selection = 'scran',
#                              ntop = 1000)


# all
meta <- meta_mat
seacell_pca <- run_pca(t(seacell_mat[meta$seacell_id,]),
                       meta,
                       HVG_selection = 'scran',
                       ntop = 2500)


save(seacell_pca_human_adult, seacell_pca_human_fetal, seacell_pca_mouse_adult, seacell_pca_mouse_fetal,  file = 'data/scEiaD_seacell_pca_objs.Rdata')
save(seacell_pca, file = 'data/scEiaD_seacell_pca_all.Rdata')
