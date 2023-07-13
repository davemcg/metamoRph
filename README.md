# metamoRph

A framework (or "guardrails") for projecting new data onto an reference PCA space. This can be a two step process where the user runs our `run_pca` function (which wraps `prcomp` and provides some sensible defaults and enhanced outputs. After `run_pca` the user can then use `metamoRph` to project new data onto the existing PCA. 

# Quick Start
```
remotes::install_github("davemcg/metamoRph")
library(metamoRph)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)

genes <- c('RHO','RPE65','ABCA4','PMEL', 'KRT5')
samples <- c('Retina1','Retina2','RPE1','RPE2','Cornea1')

faux_mat <- cbind(c(560,650,5,6,0), # rho
                  c(42,32, 1103,1201,2), #rpe65
                  c(810,903,202,205,45),  #abca4
                  c(100,105,2004,1980,101),# pmel
                  c(3,32,101,202,1567)) |> data.frame()# krt5
colnames(faux_mat) <- genes
row.names(faux_mat) <- samples

new_data <- cbind(c(5,6,3), # rho
                  c(32,23,54), #rpe65
                  c(65,22,10),  #abca4
                  c(122,101,57),# pmel
                  c(2567,1755,2218)) * 100 # krt5
colnames(new_data) <- genes
row.names(new_data) <- c("Cornea2","Cornea3","Cornea4")


mm_pca <- run_pca(t(faux_mat), meta = samples |> data.frame())
projected_pca <- metamoRph(t(new_data), 
                               mm_pca$PCA$rotation, 
                               center_scale = mm_pca$center_scale)

# bind_rows(as_tibble(mm_pca$PCA$x, rownames = 'samples'),
#       data.frame((projected_pca)) |>mutate(samples = 'Cornea')) |>
#   mutate(samples = gsub('\\d','',samples)) |>
#   ggplot(aes(x=PC1,y=PC2,color=samples)) + 
#   geom_point() +
#   cowplot::theme_cowplot()
```

## Label projection
```
# continue from code chunk above
## WARNING: Use many more num_PCs (20+) for "real" genomic data
trained_model <- model_build(mm_pca$PCA$x,
                             mm_pca$meta$samples,
                             model = 'lm', num_PCs = 2)

label_guesses <- model_apply(trained_model,
                             projected_pca,
                             c("Cornea","Cornea","Cornea")
)

# label_guesses
```

It is also possible to use a pre-existing prcomp object directly with `metamoRph`.

## Existing prcomp object projection example (pseudocode-ish)
```
library(metamoRph)
projected_pca <- metamoRph(t(new_data), 
                               your_prcomp_object$PCA$rotation, 
                               center_scale = extract_prcomp_scaling(your_prcomp_object))
```



