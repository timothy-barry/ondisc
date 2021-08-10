# set parameters
multimodal_ncol <- sample(x = seq(500, 1500), size = 1, replace = TRUE)
mulitmodal_nrow <- sample(x = seq(500, 1500), size = 2)
multimodal_logical <- sample(c(TRUE, FALSE), size = 2, replace = TRUE)

# create metadata odm objects
metadata_odms <- lapply(c(1, 2), function(i) {
  fp <- create_synthetic_data(n_row = mulitmodal_nrow[i], n_col = multimodal_ncol,
                              logical_mat = multimodal_logical[i], write_as_mtx_to_disk = TRUE)
  file_dir <- create_new_directory()
  metadata_odm <- create_ondisc_matrix_from_mtx(mtx_fp = fp$matrix_fp, barcodes_fp = fp$barcodes_fp,
                                                features_fp = fp$features_fp, odm_fp = file_dir, progress = FALSE)
  return(metadata_odm)
})
names(metadata_odms) <- c("modality_1", "modality_2")
multimodal_mat <- multimodal_ondisc_matrix(covariate_ondisc_matrix_list = metadata_odms)

# run tests
test_that("multimodal_odm show", {
  show(multimodal_mat)
  expect_true(TRUE)
})


test_that("multimodal_odm get feature ids, names, and cell barcodes", {
  expect_equal(get_feature_ids(multimodal_mat),
               lapply(multimodal_mat@modalities, get_feature_ids))
  expect_equal(get_feature_ids(multimodal_mat),
               lapply(multimodal_mat@modalities, get_feature_ids))
  expect_equal(get_cell_barcodes(multimodal_mat),
               get_cell_barcodes(multimodal_mat@modalities[[1]]))
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


test_that("extract multimodal_odm", {
    expect_error(multimodal_mat[[1,]])
    expect_error(multimodal_mat[[,1]])
    expect_error(multimodal_mat[[1,1]])
})


test_that("covariate matrices", {
  expect_error(get_feature_covariates(multimodal_mat))
  expect_error(mutate_feature_covariates(multimodal_mat))
  expect_equal(multimodal_mat@global_cell_covariates,
                        get_cell_covariates(multimodal_mat))
  expect_true("test_col" %in% (mutate_cell_covariates(multimodal_mat, test_col = 1) %>%
                                 get_cell_covariates() %>% colnames()))
})


test_that("get modality", {
  expect_identical(get_modality(multimodal_mat, "modality_1"), multimodal_mat@modalities$modality_1)
})
