test_that("show and head", {
  for (i in 1:n_datasets) {
    Mat <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    # show and head methods -- verify they run without error
    show(on_disc_mat)
    head(on_disc_mat)
  }
  expect_true(TRUE)
})
