test_that("summarize expression matrix", {
  for (i in 1:n_datasets) {
    test_obj <- load_on_disc_and_mat(data_dir = temp_test_dir, idx = i)
    m <- test_obj$r_Matrix
    on_disc_obj <- test_obj$on_disc_matrix
    cov_mats <- summarize_expression_matrix(x = on_disc_obj, chunk_size = floor(ncol(on_disc_obj)/3))

    cell_lib_size <- Matrix::colSums(m)
    n_genes_expressed <- Matrix::colSums(m >= 1)
    gene_lib_size <- Matrix::rowSums(m)
    n_cells_expressed <- Matrix::rowSums(m >= 1)

    expect_equal(cell_lib_size, cov_mats$cell_covariate_matrix[,"total_umis"])
  }
})


test_that("apply", {
  for (i in 1:n_datasets) {
    test_obj <- load_on_disc_and_mat(data_dir = temp_test_dir, idx = i)
    m <- test_obj$r_Matrix
    on_disc_obj <- test_obj$on_disc_matrix

    gene_sds_1 <- apply(X = on_disc_obj, MARGIN = 1, FUN = sd, chunk_size = floor(nrow(on_disc_obj)/3))
    gene_sds_2 <- apply(X = m, MARGIN = 1, FUN = sd)
    expect_equal(gene_sds_1, gene_sds_2)

    cell_sds_1 <- apply(X = on_disc_obj, MARGIN = 2, FUN = sd, chunk_size = floor(ncol(on_disc_obj)/3))
    cell_sds_2 <- apply(X = m, MARGIN = 2, FUN = sd)
    expect_equal(cell_sds_1, cell_sds_2)
  }
})
