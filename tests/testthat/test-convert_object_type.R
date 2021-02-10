test_that("Convert to Seurat", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    s_obj <- convert_to_seurat_object(on_disc_mat)
    expect_true(all(Mat == s_obj@assays$RNA@counts))
  }
})
