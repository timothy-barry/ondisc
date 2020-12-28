test_params <- get_test_parameters(get_test_type())
set.seed(test_params$seed)

test_that("gets, no subsets", {
  for (i in 1:test_params$n_datasets) {
    if (test_params$n_datasets > 1) cat(paste0("Running test ", i, ".\n"))
    test_obj <- load_on_disc_and_mat(data_dir = test_params$synthetic_data_dir, idx = i)
    on_disc_mat <- test_obj$on_disc_matrix
    all(paste0("gene_", 1:nrow(on_disc_mat)) == get_gene_names(on_disc_mat)) %>% expect_true()
    all(paste0("ENSG000", 1:nrow(on_disc_mat)) == get_gene_ids(on_disc_mat)) %>% expect_true()
    all(paste0("cell_", 1:ncol(on_disc_mat)) == get_cell_barcodes(on_disc_mat)) %>% expect_true()
    }
})

test_that("gets after subset", {
  for (i in 1:test_params$n_datasets) {
    if (test_params$n_datasets > 1) cat(paste0("Running test ", i, ".\n"))
    test_obj <- load_on_disc_and_mat(data_dir = test_params$synthetic_data_dir, idx = i)
    on_disc_mat <- test_obj$on_disc_matrix
      for (j in 1:test_params$n_reps_per_dataset) {
        cat(paste0("\tRunning sub-test ", j, ".\n"))
        subset_size_col <- sample(1:(ceiling(ncol(on_disc_mat)/10)), 1)
        subset_size_row <- sample(1:(ceiling(nrow(on_disc_mat)/10)), 1)
        col_names <- paste0("cell_",sample(x = 1:ncol(on_disc_mat), size = subset_size_col))
        row_names <- paste0("ENSG000", sample(x = 1:nrow(on_disc_mat), size = subset_size_row))
        # subset a first time
        t1 <- on_disc_mat[,col_names]
        expect_true(all(get_cell_barcodes(t1) == col_names))
        t2 <- on_disc_mat[row_names,]
        expect_true(all(get_gene_ids(t2) == row_names))
        expect_true(all(get_gene_names(t2) == paste0("gene_", t2@gene_subset)))
        # subset a second time
        subset_size_col_2 <- sample(1:length(col_names), 1)
        subset_size_row_2 <- sample(1:length(row_names), 1)
        col_names <- sample(col_names, subset_size_col_2)
        row_names <- sample(row_names, subset_size_row_2)
        t1 <- t1[,col_names]
        expect_true(all(get_cell_barcodes(t1) == col_names))
        t2 <- t2[row_names,]
        expect_true(all(get_gene_ids(t2) == row_names))
        expect_true(all(get_gene_names(t2) == paste0("gene_", t2@gene_subset)))
      }
    }
})
