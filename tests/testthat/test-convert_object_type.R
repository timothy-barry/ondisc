test_that("Convert to Seurat", {
  for (i in 1:n_datasets) {
    if (n_datasets > 1) cat(paste0("Running test ", i, ".\n"))
    test_obj <- load_on_disc_and_mat(data_dir = temp_test_dir, idx = i)
    on_disc_mat <- test_obj$on_disc_matrix
    m <- test_obj$r_Matrix
    s_obj <- convert_to_seurat_object(on_disc_mat)
    expect_true(all(m == s_obj@assays$RNA@counts))
  }
})
