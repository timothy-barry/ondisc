test_that("save and read odm", {
  cov_odm_1 <- cov_odms[[1]]
  cov_odm_sub <- cov_odm_1[get_random_subset(nrow(cov_odm_1)),
                           get_random_subset(ncol(cov_odm_1))]
  metadata_fp <- paste0(create_new_directory(), "/metadata.rds")
  save_odm(odm = cov_odm_sub, metadata_fp = metadata_fp)
  loaded_cov_odm_sub <- read_odm(odm_fp = cov_odm_sub@ondisc_matrix@h5_file,
                                 metadata_fp = metadata_fp)
  expect_true(identical(loaded_cov_odm_sub, cov_odm_sub))
  expect_equal(loaded_cov_odm_sub[[1:nrow(loaded_cov_odm_sub)]],
              cov_odm_sub[[1:nrow(cov_odm_sub)]])

  odm_1 <- cov_odm_1@ondisc_matrix
  save_odm(odm = odm_1, metadata_fp = metadata_fp)
  odm_loaded <- read_odm(odm_fp = odm_1@h5_file)
  expect_equal(odm_1[[1:nrow(odm_1)]],
              odm_loaded[[1:nrow(odm_loaded),]])
})


test_that("wrong metadata file", {
  cov_odm_1 <- cov_odms[[1]]
  cov_odm_2 <- cov_odms[[2]]
  metadata_fp_1 <- paste0(create_new_directory(), "/metadata_1.rds")
  save_odm(cov_odm_1, metadata_fp_1)
  expect_error(read_odm(odm_fp = cov_odm_2@ondisc_matrix@h5_file,
           metadata_fp = metadata_fp_1))
})
