test_that("summarize expression matrix", {
  for (i in 1:n_datasets) {
    test_obj <- load_on_disc_and_mat(data_dir = temp_test_dir, idx = i)
    m <- test_obj$r_Matrix
    on_disc_obj <- test_obj$on_disc_matrix

    # if i is even, exclude mitochondrial genes by subsetting
    gene_ids <- get_gene_names(on_disc_obj)
    mt_genes <- grep(pattern = "^MT-", x = gene_ids)
    if (i %% 2 == 0) {
      on_disc_obj <- on_disc_obj[-mt_genes,]
      m <- m[-mt_genes,,drop=FALSE]
      p_mito <- rep(0, nrow(m))
    }
    cov_mats <- summarize_expression_matrix(x = on_disc_obj, n_cells_to_process_at_time = floor(ncol(on_disc_obj)/3))
    cell_lib_size <- Matrix::colSums(m)
    n_genes_expressed <- Matrix::colSums(m >= 1)
    if (i %% 2 == 0) {
      p_mito <- rep(0, ncol(m))
    } else {
      cell_mt_count <- Matrix::colSums(m[mt_genes,,drop=FALSE])
      p_mito <- cell_mt_count/cell_lib_size * 100
    }
    gene_lib_size <- Matrix::rowSums(m)
    n_cells_expressed <- Matrix::rowSums(m >= 1)

    expect_true(all(cbind(cell_lib_size, n_genes_expressed, p_mito) == cov_mats$cell_covariate_matrix))
    expect_true(all(cbind(gene_lib_size, n_cells_expressed) == cov_mats$gene_covariate_matrix))
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

    cell_stats <- apply(X = on_disc_obj, MARGIN = 2, FUN = function(col) c(m = mean(col), sd = sd(col)), chunk_size = floor(ncol(on_disc_obj)/3))
    cell_means <- apply(X = m, MARGIN = 2, FUN = mean)

    expect_true(all(rbind(cell_means, cell_sds_2) == cell_stats))
  }
})
