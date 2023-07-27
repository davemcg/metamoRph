# metamoRph

A framework for projecting query (new) data onto an reference PCA space. This is easiest as a two step process where the user runs our `run_pca` function (which wraps `prcomp` and provides some sensible defaults and enhanced outputs. After `run_pca` the user can then use `metamoRph` to normalize the query (new) data in the same manner as the reference data and then project it onto the existing reference PCA. 

We also provides two functions, `model_build` and `model_apply` to very quickly transfer metadata onto new data using the shared PCA space between the reference and the query.

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

# example plot
# bind_rows(as_tibble(mm_pca$PCA$x, rownames = 'samples'),
#       as_tibble(projected_pca, rownames = 'Cornea')) |>
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
                             gsub('\\d+','', mm_pca$meta$samples), #remove trailing digit 
                             model = 'lm', num_PCs = 2)

label_guesses <- model_apply(trained_model,
                             projected_pca,
                             c("Cornea","Cornea","Cornea")
)

# label_guesses
```



## Existing prcomp object projection example (pseudocode-ish)

It is also possible to use a pre-existing prcomp object with `metamoRph`, skipping
the `run_pca` step. You will have to ensure that you sample scale your `new_data`
in the same manner that you used for your `prcomp` run (here we use `log2` - but replace
with whatever scaling your have done. Maybe none?). Note how we have turned
off the cpm and log1p scaling in the `metamoRph` function with `normaliztion = FALSE`.

```
library(metamoRph)
projected_pca <- metamoRph(t(log2(new_data)), 
                               your_prcomp_object$PCA$rotation, 
                               center_scale = extract_prcomp_scaling(your_prcomp_object),
                               normalization = FALSE)
```



# Some Related Tools

  - [ProjecTILs](https://github.com/carmonalab/ProjecTILs)
    - scRNA/Seurat based tool that uses a similar approach to project data (essentialy multiplying the new against the rotation matrix)
  - [Seurat](https://satijalab.org/seurat/articles/integration_mapping.html)
    - The big kahuna itself. scRNA analysis playground.
  - [projectR](https://www.bioconductor.org/packages/release/bioc/html/projectR.html)
    - Can transfer data onto PCA *and* NMF 
    - At least for PCA does not provide any guardrails for properly scaling the query data
