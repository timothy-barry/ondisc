# 1. create a new directory
create_new_directory <- function() {
  f <- tempfile()
  while (dir.exists(f)) f <- tempfile()
  dir.create(f)
  return(f)
}

# 2. create a random matrix
create_random_matrix <- function(n_row, n_col, p_zero = 0.9, p_set_col_zero = 0.05, p_set_row_zero = 0.05, matrix_values = seq(1L, 10L), matrix_class = "CsparseMatrix") {
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
  if (matrix_class == "matrix") {
    out <- r
  } else if (matrix_class == "CsparseMatrix") {
    out <- Matrix::Matrix(r, sparse = TRUE)
  } else if (matrix_class == "RsparseMatrix") {
    out <- methods::as(r, "RsparseMatrix")
  } else {
    stop("`matrix_class` not recognized.")
  }
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
  modality_names_df <- lapply(X = seq_along(modality_names), FUN = function(i) {
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

  # 5. construct the barcoes df
  barcode_df <- data.frame(cell_id = paste0("cell_", seq_len(ncol(combined_mat))))

  # 6. split the matrices according to batch; loop over the batches and save the matrix and features file
  batch_levels_v <- as.character(unique(batch)) |> sort()

  for (i in seq_along(batch_levels_v)) {
    batch_level <- batch_levels_v[i]
    mat_sub <- combined_mat[, batch_level == batch]
    barcode_df_sub <- barcode_df[batch_level == batch,,drop = FALSE]
    dir_name <- paste0(directory, "/", batch_level)
    dir.create(dir_name)
    Matrix::writeMM(obj = mat_sub, file = paste0(dir_name, "/matrix.mtx"))
    readr::write_tsv(x = feature_df, file = paste0(dir_name, "/features.tsv"), col_names = FALSE)
    readr::write_tsv(x = barcode_df_sub, file = paste0(dir_name, "/barcodes.tsv"), col_names = FALSE)
    curr_files <- list.files(dir_name, full.names = TRUE)
    for (curr_file in curr_files) {
      R.utils::gzip(filename = curr_file, destname = paste0(curr_file, ".gz"))
    }
  }
  return(NULL)
}


#' Write example Cell Ranger dataset
#'
#' `write_example_cellranger_dataset()` creates an example dataset and writes the dataset to disk in Cell Ranger format. The dataset can be unimodal or multimodal, containing some subset of the gene expression, gRNA expression, and protein expression modalities.
#'
#' @param n_features an integer vector specifying the number of features to simulate for each modality.
#' @param n_cells an integer specifying the number of cells to simulate.
#' @param n_batch an integer specifying the number of batches to simulate. Cells from different batches are written to different directories.
#' @param modalities a character vector indicating the modalities to simulate. The vector should contain some subset of the strings "gene", "grna", and "protein".
#' @param directory_to_write the directory in which to write the Cell Ranger feature barcode files.
#' @param p_zero (optional; default random) the fraction of entries to set randomly to zero.
#' @param p_set_col_zero (optional; default random) the fraction of columns to set randomly to zero.
#' @param p_set_row_zero (optional; default random) the fraction of rows to set randomly to zero.
#'
#' @return a list containing (i) a list of the simulated expression matrices, (ii) a character vector containing the names of the genes (or `NULL` if the gene modality was not simulated), and (iii) a character vector specifying the batch of each cell.
#' @export
#'
#' @examples
#' set.seed(4)
#' n_features <- c(1000, 40, 400)
#' modalities <- c("gene", "protein", "grna")
#' n_cells <- 10000
#' n_batch <- 2
#' directory_to_write <- tempdir()
#' p_set_col_zero <- 0
#' out <- write_example_cellranger_dataset(
#'   n_features = n_features,
#'   n_cells = n_cells,
#'   n_batch = n_batch,
#'   modalities = modalities,
#'   directory_to_write = directory_to_write,
#'   p_set_col_zero = p_set_col_zero
#' )
#'
#' # directories written to directory_to_write
#' fs <- list.files(directory_to_write, pattern = "batch*", full.names = TRUE)
#' # files contained within the directories
#' list.files(fs[1])
#' list.files(fs[2])
write_example_cellranger_dataset <- function(n_features, n_cells, n_batch, modalities, directory_to_write,
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
  if (is.null(p_zero)) p_zero <- min(stats::runif(1), 0.9)
  if (is.null(p_set_col_zero)) p_set_col_zero <- min(stats::runif(1), 0.9)
  if (is.null(p_set_row_zero)) p_set_row_zero <- min(stats::runif(1), 0.9)
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
                                      frac_mito = stats::runif(1, min = 0, max = 0.5))
  } else {
    gene_names <- NULL
  }

  # generate the batch
  batch <- generate_batch(n_cells, n_batch)

  # write the data
  write_sceptre_object_to_cellranger_format_v2(mats = mats, gene_names = gene_names,
                                               directory = directory_to_write, batch = batch)

  # return the data
  return(list(matrix_list = mats, gene_names = gene_names, batch = batch))
}
