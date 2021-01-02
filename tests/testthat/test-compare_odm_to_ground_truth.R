test_that("Compare ground truth matrix to on_disc_matrix", {
  for (i in 1:n_datasets) {
    test_obj <- load_on_disc_and_mat(data_dir = temp_test_dir, idx = i)
    m1 <- test_obj$r_Matrix
    on_disc_obj <- test_obj$on_disc_matrix
    m2 <- on_disc_obj[[,1:ncol(on_disc_obj)]]
    m3 <- on_disc_obj[[1:nrow(on_disc_obj),]]
    expect_true(all(m1 == m2))
    expect_true(all(m2 == m3))
  }
})
