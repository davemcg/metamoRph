library(TENxIO)
library(tidyverse)
library(scran)
source('src/pca_projector.R')

# read in scp1755
scp1755 <- import(TENxFile('~/Downloads/SCP1755.tar.gz'))
scp1755_meta <- read_tsv('~/Downloads/SCP1755_scPortal_metadata_file.txt')

scp1755_meta <- scp1755_meta[2:nrow(scp1755_meta),]

colData(scp1755)$NAME <- colnames(assay(scp1755))

colData(scp1755)$celltype <- colData(scp1755) %>% as_tibble() %>% left_join(scp1755_meta) %>% pull(cell_type__custom)

scp1755_sum_counts <- aggregateAcrossCells(scp1755, colData(scp1755)$celltype) %>% assay()
row.names(scp1755_sum_counts) <- rowData(scp1755)$Symbol

# project onto our adult human pca scEiaD data
scp1755_project  <- eyeProjector(rotation = 'adult_human', scp1755_sum_counts)


# first no harmony
custom.settings = umap.defaults
custom.settings$n_neighbors = 30
custom.settings$metric <- 'euclidean'
custom.settings$min_dist <- 0.5
seacell_umap <- umap::umap(seacell_pca_human_adult$PCA$x)
p <- seacell_umap$layout %>% as_tibble(rownames = 'info') %>%
  left_join(seacell_pca_human_adult$meta, by = c('info' = 'seacell_id')) %>%
  mutate(label = paste0(study_accession, '\n',CTP)) %>%
  mutate(CellType = case_when(is.na(CellType) ~ tolower(info), TRUE ~ CellType)) %>%
  ggplot(aes(x=V1,y=V2,color=CellType, label = label)) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2(), pals::alphabet2()) %>% unname())
plotly::ggplotly(p)

p <- harmony_test %>% as_tibble(rownames = 'info') %>%
  left_join(seacell_summarise, by = c('info' = 'seacell_id')) %>%
  mutate(label = paste0(study_accession, '\n',CTP)) %>%
  ggplot(aes(x=PC1,y=PC2,color=CellType, label = label)) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2(), pals::alphabet2()) %>% unname())
plotly::ggplotly(p)



#######
set.seed(1234)
harmony_test <- harmony::HarmonyMatrix(rbind(scp1755_project[[1]][,1:200],
                                             seacell_pca_human_adult[[1]]$x[,1:200]),
                                       do_pca = FALSE,
                                       meta_data = c(rep('scp1755', 13), seacell_pca_human_adult[[1]]$x %>% row.names() %>% gsub('__\\d+','',.)),
                                       max.iter.harmony = 20)

custom.settings = umap.defaults
custom.settings$n_neighbors = 25
custom.settings$metric <- 'euclidean'
custom.settings$min_dist <- 0.5
seacell_umap <- umap::umap(harmony_test, config = custom.settings)
p <- seacell_umap$layout %>% as_tibble(rownames = 'info') %>%
  left_join(seacell_summarise, by = c('info' = 'seacell_id')) %>%
  mutate(label = paste0(study_accession, '\n',CTP)) %>%
  mutate(CellType = case_when(is.na(CellType) ~ tolower(info), TRUE ~ CellType)) %>%
  ggplot(aes(x=V1,y=V2,color=CellType, label = label)) +
  geom_point(size = 0.5) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2(), pals::alphabet2()) %>% unname())
plotly::ggplotly(p)



#big_pc <- rbind(scp1755_project[[1]], scp1755_project[[2]][[1]]$x)
dist_big <- dist(harmony_cor)

(dist_big %>% as.matrix()) %>%
  as_tibble(rownames = 'A') %>%
  pivot_longer(-A, names_to = 'B') %>%
  filter(A == 'Muller glia') %>%
  data.frame() %>% arrange(value) %>%
  left_join(seacell_summarise %>% ungroup() %>%
              select(seacell_id,CTP), by = c('B' = 'seacell_id')) %>%
  filter(!is.na(CTP)) %>%
  head(9)

