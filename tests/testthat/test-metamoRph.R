library(testthat)
library(metamoRph)

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


test_that("cpm scaling only, metamoRph projects new cornea approx the same as original cornea", {
  # Define input data for the test
  new_counts <- faux_mat
  # Call the function
  pca_result <- run_pca(t(faux_mat), meta = samples |> data.frame(),
                    sample_cpm_scale = TRUE,
                    log1p = FALSE,
                    feature_scale = FALSE,
                    feature_center = FALSE)
  
  morph_result <- metamoRph(t(new_data) * 100, 
                            pca_result$PCA$rotation, 
                            center_scale = pca_result$center_scale,
                            sample_cpm_scale = TRUE,
                            log1p = FALSE)
                            
  expect_equal(pca_result$PCA$x[5,1],  morph_result[1,1])
  # Add more assertions as needed
})

test_that("cpm and log1p scaling only, metamoRph projects new cornea approx the same as original cornea", {
  # Define input data for the test
  new_counts <- faux_mat
  # Call the function
  pca_result <- run_pca(t(faux_mat) , meta = samples |> data.frame(),
                        sample_cpm_scale = TRUE,
                        log1p = TRUE,
                        feature_scale = FALSE,
                        feature_center = FALSE)
  
  morph_result <- metamoRph(t(new_data) * 100, 
                            pca_result$PCA$rotation, 
                            center_scale = pca_result$center_scale,
                            sample_cpm_scale = TRUE,
                            log1p = TRUE)
  
  # Make assertions to check if the output is as expected
  expect_equal(pca_result$PCA$x[5,1],  morph_result[1,1])
})

test_that("log1p scaling only, metamoRph projects new cornea as original cornea", {
  # Define input data for the test
  new_counts <- faux_mat
  # Call the function
  pca_result <- run_pca(t(faux_mat), meta = samples |> data.frame(),
                        sample_cpm_scale = FALSE,
                        log1p = TRUE,
                        feature_scale = FALSE,
                        feature_center = FALSE)
  
  morph_result <- metamoRph(t(new_data), 
                            pca_result$PCA$rotation, 
                            center_scale = pca_result$center_scale,
                            sample_cpm_scale = FALSE,
                            log1p = TRUE)
  
  # Make assertions to check if the output is as expected
  expect_equal(pca_result$PCA$x[5,1],  morph_result[1,1])
})

test_that("cpm,log1p, and center-scale with features metamoRph projects new cornea approx the same as original cornea", {
  # Define input data for the test
  new_counts <- faux_mat
  # Call the function
  pca_result <- run_pca(t(faux_mat), meta = samples |> data.frame(),
                        sample_cpm_scale = TRUE,
                        log1p = TRUE,
                        feature_scale = TRUE,
                        feature_center = TRUE)
  
  morph_result <- metamoRph(t(new_data) * 100, 
                            pca_result$PCA$rotation, 
                            center_scale = pca_result$center_scale,
                            sample_cpm_scale = TRUE,
                            log1p = TRUE)
  
  # Make assertions to check if the output is as expected
  expect_equal(pca_result$PCA$x[5,1],  morph_result[1,1])
})
