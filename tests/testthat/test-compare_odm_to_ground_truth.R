test_that("Compare ground truth matrix to on_disc_matrix", {
  for (i in seq(1, n_datasets)) {
    m1 <- r_mats[[i]]
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    m2 <- on_disc_mat[[,1:ncol(on_disc_mat)]]
    m3 <- on_disc_mat[[1:nrow(on_disc_mat),]]
    m4 <- on_disc_mat[[1:nrow(on_disc_mat),1:ncol(on_disc_mat)]]
    expect_true(all(m1 == m2))
    expect_true(all(m2 == m3))
    expect_true(all(m3 == m4))
  }
})
