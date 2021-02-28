test_that("subset metadata_odm", {
  for (i in 1:n_datasets) {
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


test_that("subset multimodal_odm", {
  # subsetting by feature throws error
  expect_error(multimodal_mat[1,])
  expect_error(multimodal_mat[1,1])
  # subset by neither
  expect_equal(multimodal_mat[], multimodal_mat)
  # subset by cell
  col_idxs <- get_random_subset(ncol(multimodal_mat))
  x <- multimodal_mat[,col_idxs]
  expect_equal(x@global_cell_covariates, multimodal_mat@global_cell_covariates[col_idxs,])
  for (i in seq(1, length(multimodal_mat@modalities))) {
    modality_sub <- x@modalities[[i]]
    modality <-multimodal_mat@modalities[[i]]
    expect_equal(ncol(modality_sub), length(col_idxs))
    expect_equal(modality_sub@feature_covariates, modality@feature_covariates)
    expect_equal(ncol(modality_sub@ondisc_matrix), length(col_idxs))
    expect_equal(modality_sub@cell_covariates, modality@cell_covariates[col_idxs,,drop=FALSE])
  }
})


test_that("extract metadata_odm, multimodal_odm", {
  for (i in 1:n_datasets) {
    metadata_odm <- cov_odms[[i]]
    expect_error(metadata_odm[[1,]])
    expect_error(metadata_odm[[,1]])
    expect_error(metadata_odm[[1, 1]])
  }
  expect_error(multimodal_mat[[1,]])
  expect_error(multimodal_mat[[,1]])
  expect_error(multimodal_mat[[1,1]])
})
