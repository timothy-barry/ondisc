test_that("pairwise QC stops at the end of pair vectors", {
  matrix_in <- Matrix::Matrix(
    matrix(c(1, 0, 0, 0, 1, 0), nrow = 2L, byrow = TRUE),
    sparse = TRUE
  )
  rownames(matrix_in) <- c("gene_1", "gene_2")
  odm_file <- tempfile(fileext = ".odm")
  on.exit(unlink(odm_file))
  odm <- create_odm_from_r_matrix(
    mat = matrix_in,
    file_to_write = odm_file,
    chunk_size = 1L
  )
  row_pointer <- ondisc:::read_row_ptr(odm@h5_file, nrow(odm))

  expect_silent(
    full_result <- ondisc:::compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc(
      file_name_in = odm@h5_file,
      f_row_ptr = row_pointer,
      n_genes = 2L,
      n_cells_orig = 3L,
      n_cells_sub = 3L,
      grna_group_idxs = list(1L),
      indiv_nt_grna_idxs = list(1L),
      all_nt_idxs = 2L,
      to_analyze_response_idxs = 1L,
      to_analyze_grna_idxs = 1L,
      control_group_complement = TRUE,
      cells_in_use = 1:3
    )
  )
  expect_identical(full_result$n_nonzero_trt, 1L)
  expect_identical(full_result$n_nonzero_cntrl, 0L)

  expect_silent(
    trans_result <- ondisc:::compute_n_ok_pairs_ondisc(
      file_name_in = odm@h5_file,
      f_row_ptr = row_pointer,
      n_genes = 2L,
      n_cells_orig = 3L,
      n_cells_sub = 3L,
      grna_group_idxs = list(1L),
      all_nt_idxs = 2L,
      to_analyze_response_idxs = 1L,
      to_analyze_grna_idxs = 1L,
      control_group_complement = TRUE,
      cells_in_use = 1:3,
      unique_response_idxs = 1L
    )
  )
  expect_identical(trans_result$n_nonzero_trt, 1L)
  expect_identical(trans_result$n_nonzero_cntrl, 0L)
})
