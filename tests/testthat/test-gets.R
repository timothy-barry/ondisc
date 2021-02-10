test_that("gets, no subsets", {
  for (i in 1:n_datasets) {
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
    df <- readr::read_tsv(file = paste0(temp_test_dir, "/features_", i,".tsv"), col_names = c("id", "name"), col_types = c("cc"))
    all(df$name == get_feature_names(on_disc_mat)) %>% expect_true()
    all(paste0("ENSG000", 1:nrow(on_disc_mat)) == get_feature_ids(on_disc_mat)) %>% expect_true()
    all(paste0("cell_", 1:ncol(on_disc_mat)) == get_cell_barcodes(on_disc_mat)) %>% expect_true()
    }
})


test_that("gets after subset", {
  for (i in 1:n_datasets) {
    on_disc_mat <- cov_odms[[i]]@ondisc_matrix
      for (j in 1:n_reps) {
        subset_size_col <- sample(1:(ceiling(ncol(on_disc_mat)/10)), 1)
        subset_size_row <- sample(1:(ceiling(nrow(on_disc_mat)/10)), 1)
        col_names <- get_cell_barcodes(on_disc_mat)
        row_names <- get_feature_ids(on_disc_mat)
        # subset a first time
        t1 <- on_disc_mat[,col_names]
        expect_true(all(get_cell_barcodes(t1) == col_names))
        t2 <- on_disc_mat[row_names,]
        expect_true(all(get_feature_ids(t2) == row_names))
        # subset a second time
        subset_size_col_2 <- sample(1:length(col_names), 1)
        subset_size_row_2 <- sample(1:length(row_names), 1)
        col_names <- sample(col_names, subset_size_col_2)
        row_names <- sample(row_names, subset_size_row_2)
        t1 <- t1[,col_names]
        expect_true(all(get_cell_barcodes(t1) == col_names))
        t2 <- t2[row_names,]
        expect_true(all(get_feature_ids(t2) == row_names))
      }
    }
})
