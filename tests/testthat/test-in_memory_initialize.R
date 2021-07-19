if (FALSE) {

# Initialize metadata odms from in-memory R objects in different directories
cov_odms_from_memory <- lapply(r_mats_plus_metadata, function(l) {
  file_dir <- create_new_directory()
  metadata_odm <- create_ondisc_matrix_from_R_matrix(r_matrix = as.matrix(l$r_mat),
                                                     barcodes = l$barcodes,
                                                     features_df = l$features_df,
                                                     return_metadata_ondisc_matrix = TRUE,
                                                     on_disk_dir = file_dir)
})


test_that("Compare ground truth matrix to on_disc_matrix initialized in memory", {
  for (i in seq(1, n_datasets)) {
    m1 <- r_mats[[i]]
    on_disc_mat <- cov_odms_from_memory[[i]]@ondisc_matrix
    m2 <- on_disc_mat[[,1:ncol(on_disc_mat)]]
    m3 <- on_disc_mat[[1:nrow(on_disc_mat),]]
    m4 <- on_disc_mat[[1:nrow(on_disc_mat),1:ncol(on_disc_mat)]]
    expect_true(all(m1 == m2))
    expect_true(all(m2 == m3))
    expect_true(all(m3 == m4))
  }
})


test_that("Check dimension",{
  for (i in seq(1, n_datasets)) {
    original_dim <- dim(r_mats[[i]])
    test_dim <- dim(cov_odms_from_memory[[i]]@ondisc_matrix)
    expect_true(all(original_dim == test_dim))
  }
})

# add test to check covariates here


}
