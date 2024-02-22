create_odm_from_r_matrix <- function(mat, file_to_write, chunk_size = 1000L, compression_level = 3L, integer_id = 0) {
  if (!is(mat, "dgRMatrix")) stop("`mat` must be an object of class `dgRMatrix`.")
  # create odm
  create_odm_r_matrix_cpp(file_name_in = file_to_write,
                          feature_ids = rownames(mat),
                          n_features = nrow(mat),
                          n_cells = ncol(mat),
                          integer_id = integer_id,
                          chunk_size = chunk_size,
                          compression_level = compression_level,
                          j = mat@j, x = as.integer(mat@x), p = mat@p)
}
