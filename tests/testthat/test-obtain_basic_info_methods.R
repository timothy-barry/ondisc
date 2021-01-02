test_that("show and head", {
  for (i in 1:n_datasets) {
    test_obj <- load_on_disc_and_mat(data_dir = temp_test_dir, idx = i)
    on_disc_mat <- test_obj$on_disc_matrix
    # show and head methods -- verify they works without error
    show(on_disc_mat)
    head(on_disc_mat)
  }
  expect_true(TRUE)
})

