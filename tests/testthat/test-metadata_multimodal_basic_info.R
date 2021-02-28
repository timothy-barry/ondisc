test_that("metadata_odm, multimodal_odm show", {
  for (i in 1:n_datasets) {
    metadata_odm <- cov_odms[[i]]
    show(metadata_odm)
  }
  show(multimodal_mat)
  expect_true(TRUE)
})


test_that("metadata_odm, multimodal_odm get feature ids, names, and cell barcodes", {
  for (i in 1:n_datasets) {
    metadata_odm <- cov_odms[[i]]
    expect_equal(get_feature_ids(metadata_odm),
                 get_feature_ids(metadata_odm@ondisc_matrix))
    expect_equal(get_feature_names(metadata_odm),
                 get_feature_names(metadata_odm@ondisc_matrix))
    expect_equal(get_cell_barcodes(metadata_odm),
                 get_cell_barcodes(metadata_odm@ondisc_matrix))
  }
  expect_equal(get_feature_ids(multimodal_mat),
               lapply(multimodal_mat@modalities, get_feature_ids))
  expect_equal(get_feature_ids(multimodal_mat),
               lapply(multimodal_mat@modalities, get_feature_ids))
  expect_equal(get_cell_barcodes(multimodal_mat),
               get_cell_barcodes(multimodal_mat@modalities[[1]]))
})
