# 1. create a new directory
create_new_directory <- function() {
  f <- tempfile()
  while (dir.exists(f)) f <- tempfile()
  dir.create(f)
  return(f)
}

# 2. create a random matrix
create_random_matrix <- function(n_row, n_col, p_zero = 0.9, p_set_col_zero = 0.05, p_set_row_zero = 0.05, matrix_values = seq(1L, 10L), column_access = TRUE) {
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
  out <- if (column_access) Matrix::Matrix(r, sparse = TRUE) else as(r, "RsparseMatrix")
  return(out)
}

# 3. add row names
add_row_names <- function(mat, modality) {
  if (modality == "gene") {
    rownames(mat) <- paste0("ENSG00000", seq(1, nrow(mat)))
  } else if (modality == "grna") {
    rownames(mat) <- paste0("grna_", seq(1, nrow(mat)))
  } else if (modality == "protein") {
    rownames(mat) <- paste0("protein_", seq(1, nrow(mat)))
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
generate_batch <- function(n_cells, n_batches) {
  if (n_batches == 1) {
    rep("batch_1", n_cells) |> factor()
  } else {
    cut(seq(1L, n_cells), breaks = n_batches, labels = paste0("batch_", seq(1L, n_batches)))
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


# 6. write to cellranger format
write_sceptre_object_to_cellranger_format_v2 <- function(mats, gene_names, directory, batch) {
  # 0. create directory
  if (dir.exists(directory)) unlink(directory, recursive = TRUE)
  dir.create(directory)

  # 1. combine the matrices
  combined_mat <- Reduce(function(x, y) rbind(x, y), mats)

  # 2. get the ids
  ids_df <- rownames(combined_mat)

  # 3. construct the modality names
  modality_names <- names(mats)
  modality_names_df <- sapply(X = seq_along(modality_names), FUN = function(i) {
    modality_name <- modality_names[i]
    n_feats <- nrow(mats[[i]])
    switch(EXPR = modality_name,
           "gene" = "Gene Expression",
           "protein" = "Antibody Capture",
           "grna" = "CRISPR Guide Capture") |> rep(n_feats)
  }) |> unlist()


  # 3. construct the names
  names_df <- ids_df
  if ("gene" %in% modality_names) {
    gene_idxs <- which(modality_names_df == "Gene Expression")
    names_df[gene_idxs] <- gene_names
  }

  # 4. construct the features df
  feature_df <- data.frame(id = ids_df, name = names_df,
                           modality = modality_names_df)

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


#' Write example cellranger dataset
#'
#' @param n_feature_vect number of features per matrix
#' @param n_cells number of cells
#' @param n_batch number of batches
#' @param modalities string indicating the modalities to create
#' @param dir_to_write the directory in which to write the files
#'
#' @return NULL
#' @export
#'
#' @examples
#' n_features <- c(1000, 40, 400)
#' modalities <- c("gene", "protein", "grna")
#' n_cells <- 10000
#' n_batch <- 2
#' dir_to_write <- tempdir()
#' write_example_cellranger_dataset(n_features, n_cells, n_batch, modalities, dir_to_write)
write_example_cellranger_dataset <- function(n_features, n_cells, n_batch, modalities, dir_to_write,
                                             p_zero = NULL, p_set_col_zero = NULL, p_set_row_zero = NULL) {
  if (!all(modalities %in% c("gene", "grna", "protein"))) {
    stop("`modalities` must be a vector containing the strings 'gene', 'grna', or 'protein'.")
  }
  if (length(unique(modalities)) != length(modalities)) {
    stop("`modalities` must not contain duplicate elements.")
  }
  if (length(modalities) != length(n_features)) {
    stop("`The lengths of `modalities` and `n_features` must coincide.")
  }

  # set the p_zero, p_set_col_zero, and p_set_row_col_zero parameters
  set.seed(4)
  if (is.null(p_zero)) p_zero <- min(runif(1), 0.9)
  if (is.null(p_set_col_zero)) p_set_col_zero <- min(runif(1), 0.9)
  if (is.null(p_set_row_zero)) p_set_row_zero <- min(runif(1), 0.9)
  # generate the example matrices
  mats <- lapply(seq_along(modalities), function(i) {
    modality <- modalities[i]
    mat <- create_random_matrix(n_row = n_features[i], n_col = n_cells, p_zero = p_zero,
                                p_set_col_zero = p_set_col_zero, p_set_row_zero = p_set_row_zero) |>
      add_row_names(modality)
  }) |> stats::setNames(modalities)
  # generate the gene names (if gene modality supplied)
  if ("gene" %in% modalities) {
    gene_names <- generate_gene_names(mat = mats[["gene"]],
                                      frac_mito = runif(1, min = 0, max = 0.5))
  }
  # generate the batch
  batch <- generate_batch(n_cells, sample(x = seq(1, n_batch), size = 1))

  # write the data
  write_sceptre_object_to_cellranger_format_v2(mats = mats,
                                               gene_names = gene_names,
                                               directory = dir_to_write,
                                               batch = batch)
}
