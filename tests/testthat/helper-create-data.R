# 1. create a new directory
create_new_directory <- function() {
  f <- tempfile()
  while (dir.exists(f)) f <- tempfile()
  dir.create(f)
  return(f)
}

# 2. create a random matrix
create_random_matrix <- function(n_row, n_col, p_zero = 0.9, p_set_col_zero = 0.05, p_set_row_zero = 0.05, matrix_values = seq(1L, 10L)) {
  # sample binary matrix
  r <- matrix(data = stats::rbinom(n = n_row * n_col, size = 1, prob = 1 - p_zero), nrow = n_row, ncol = n_col)
  # randomly set some rows or columns to all zero
  zero_rows <- sample(seq(1, n_row), size = floor(p_set_row_zero * n_row), replace = FALSE)
  zero_cols <- sample(seq(1, n_col), size = floor(p_set_row_zero * n_col), replace = FALSE)
  r[zero_rows,] <- 0; r[,zero_cols] <- 0
  # initialize integer or logical matrix
  m <- sample(x = matrix_values, size = sum(r), replace = TRUE)
  r[r == 1] <- m
  # convert the matrix to sparse format and return
  out <- Matrix::Matrix(r, sparse = TRUE)
  return(out)
}

# 3. add row names
add_row_names <- function(mat, gene_modality = TRUE) {
  if (gene_modality) {
    rownames(mat) <- paste0("ENSG00000", seq(1, nrow(mat)))
  } else {
    rownames(mat) <- paste0("grna_", seq(1, nrow(mat)))
  }
  return(mat)
}

# 4. generate gene names
generate_gene_names <- function(mat, frac_mito = 0.1) {
  n_row <- nrow(mat)
  mt_idxs <- sample(x = seq(1L, n_row), size = floor(frac_mito * n_row), replace = FALSE)
  paste0(ifelse(seq(1L, n_row) %in% mt_idxs, "MT-gene_", "gene_"), seq(1, n_row))
}

# 5. generate batch
generate_batch <- function(mat, n_batches) {
  if (n_batches == 1) {
    rep("batch_1", ncol(mat)) |> factor()
  } else {
    cut(seq(1L, ncol(mat)), breaks = n_batches, labels = paste0("batch_", seq(1L, n_batches)))
  }
}

# 6. write to cellranger format
write_sceptre_object_to_cellranger_format <- function(gene_matrix, gene_names, grna_matrix, directory, batch) {
  # 0. create directory
  if (dir.exists(directory)) unlink(directory, recursive = TRUE)
  dir.create(directory)

  # 1. combine matrices across mod alities
  combined_mat <- rbind(gene_matrix, grna_matrix)

  # 2. construct the features df
  gene_ids <- rownames(gene_matrix)
  grna_ids <- rownames(grna_matrix)
  feature_df <- data.frame(gene_id = c(gene_ids, grna_ids),
                           gene_name = c(gene_names, grna_ids),
                           modality = c(rep("Gene Expression", length(gene_ids)),
                                        rep("CRISPR Guide Capture", length(grna_ids))))

  # 3. split the matrices according to batch; loop over the batches and save the matrix and features file
  batch_levels_v <- as.character(unique(batch)) |> sort()

  for (i in seq_along(batch_levels_v)) {
    batch_level <- batch_levels_v[i]
    mat_sub <- combined_mat[, batch_level == batch]
    dir_name <- paste0(directory, "/", batch_level)
    dir.create(dir_name)
    Matrix::writeMM(obj = mat_sub, file = paste0(dir_name, "/matrix.mtx"))
    readr::write_tsv(x = feature_df, file = paste0(dir_name, "/features.tsv"), col_names = FALSE)
    curr_files <- list.files(dir_name, full.names = TRUE)
    for (curr_file in curr_files) {
      R.utils::gzip(filename = curr_file, destname = paste0(curr_file, ".gz"))
    }
  }
  return(NULL)
}
