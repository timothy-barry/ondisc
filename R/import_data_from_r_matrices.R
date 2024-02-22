create_odm_from_r_matrix <- function(mat, file_name_in, integer_id = 0) {
  if (!is(mat, "dgRMatrix")) stop("`mat` must be an object of class `dgRMatrix`.")
  # create odm
  create_odm(file_name_in = file_name_in,
             n_nonzero_features = n_nonzero_features,
             feature_ids = rownames(mat),
             n_cells = ncol(mat),
             integer_id = integer_id,
             chunk_size = 1000L,
             compression_level = 3L)
  # write data to odm
  # mat@p, mat@j, mat@Dim, mat@x
}
