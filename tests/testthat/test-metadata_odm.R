test_that("metadata_odm show", {
  for (i in seq(1, n_datasets)) {
    metadata_odm <- cov_odms[[i]]
    show(metadata_odm)
  }
  expect_true(TRUE)
})


test_that("metadata_odm get feature ids, names, and cell barcodes", {
  for (i in seq(1, n_datasets)) {
    metadata_odm <- cov_odms[[i]]
    expect_equal(get_feature_ids(metadata_odm),
                 get_feature_ids(metadata_odm@ondisc_matrix))
    expect_equal(get_feature_names(metadata_odm),
                 get_feature_names(metadata_odm@ondisc_matrix))
    expect_equal(get_cell_barcodes(metadata_odm),
                 get_cell_barcodes(metadata_odm@ondisc_matrix))
  }
})


test_that("subset metadata_odm", {
  for (i in seq(1, n_datasets)) {
    metadata_odm <- cov_odms[[i]]
    orig_dim <- dim(metadata_odm)
    row_idxs <- get_random_subset(orig_dim[1])
    col_idxs <- get_random_subset(orig_dim[2])

    # subset by feature
    x <- metadata_odm[row_idxs,]
    expect_equal(nrow(x@ondisc_matrix), length(row_idxs))
    expect_equal(x@cell_covariates, metadata_odm@cell_covariates)
    expect_equal(x@feature_covariates, metadata_odm@feature_covariates[row_idxs,,drop=F])

    # subset by column
    y <- metadata_odm[,col_idxs]
    expect_equal(ncol(y), length(col_idxs))
    expect_equal(y@cell_covariates, metadata_odm@cell_covariates[col_idxs,,drop=FALSE])
    expect_equal(y@feature_covariates, metadata_odm@feature_covariates)

    # subset by both
    z <- metadata_odm[row_idxs, col_idxs]
    expect_equal(nrow(z), length(row_idxs))
    expect_equal(ncol(z), length(col_idxs))

    # subset by neither
    expect_equal(metadata_odm[], metadata_odm)
  }
})


test_that("extract metadata_odm", {
  for (i in 1:n_datasets) {
    metadata_odm <- cov_odms[[i]]
    d <- dim(metadata_odm)
    expect_equal(metadata_odm[[d[1],]], metadata_odm@ondisc_matrix[[d[1],]])
    expect_equal(metadata_odm[[,d[2]]], metadata_odm@ondisc_matrix[[,d[2]]])
    expect_equal(metadata_odm[[d[1],d[2]]], metadata_odm@ondisc_matrix[[d[1],d[2]]])
  }
})
