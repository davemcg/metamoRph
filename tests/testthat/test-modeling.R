library(testthat)
library(metamoRph)



test_that("predict tissues", {
  
  genes <- c('RHO','RPE65','ABCA4','PMEL', 'KRT5')
  samples <- c('Retina1','Retina2','RPE1','RPE2','Cornea1')
  
  faux_mat <- cbind(c(560,650,5,6,0), # rho
                    c(42,32, 1103,1201,2), #rpe65
                    c(810,903,202,205,45),  #abca4
                    c(100,105,2004,1980,101),# pmel
                    c(3,32,101,202,1567)) |> data.frame()# krt5
  colnames(faux_mat) <- genes
  row.names(faux_mat) <- samples
  
  new_data <- cbind(c(0,6,3), # rho
                    c(2,23,54), #rpe65
                    c(45,22,10),  #abca4
                    c(101,101,57),# pmel
                    c(1567,1755,2218)) # krt5
  colnames(new_data) <- genes
  row.names(new_data) <- c("Cornea2","Cornea3","Cornea4")
  
  
  pca_result <- run_pca(t(faux_mat), meta = samples |> data.frame(),
                        sample_scale = 'cpm',
                        log1p = TRUE)
  
  morph_result <- metamoRph(t(new_data) * 100, 
                            pca_result$PCA$rotation, 
                            center_scale = pca_result$center_scale,
                            sample_scale = 'cpm',
                            log1p = TRUE)
  
  trained_model <- model_build(pca_result$PCA$x,
                               gsub("\\d","",pca_result$meta$samples),
                               model = 'lm', num_PCs = 2, verbose = FALSE)
  
  label_guesses <- model_apply(trained_model,
                               morph_result,
                               gsub("\\d","",row.names(morph_result))
  )
  
  label_guesses
  
  expect_equal(sum(label_guesses$sample_label == label_guesses$predict),
               3)
})
