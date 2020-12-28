#' Create on_disc_matrix from 10x .mtx file
#' Creates an on_disc_matrix from a 10x .mtx file.
#' @param mtx_fp file path to the .mtx file
#' @param barcode_fp file path to the barcode.tsv file.
#' @param features_fp file path to the features.tsv file.
#' @param on_disc_dir directory in which to store the initialized HDF5 object
#' @export
#' @return an on_disc_matrix object
#' @examples
#' mtx_fp <- system.file("extdata", "matrix.mtx", package = "ondisc")
#' barcode_fp <- system.file("extdata", "barcodes.tsv", package = "ondisc")
#' features_fp <- system.file("extdata", "features.tsv", package = "ondisc")
#' on_disc_dir <- system.file("extdata", package = "ondisc")
#' if (system.file("extdata", "on_disc_matrix_1.h5", package = "ondisc") == "") {
#' # if on_disc_matrix does not exist, create it.
#' exp_mat <- create_on_disc_matrix_from_10x_mtx(mtx_fp, barcode_fp, features_fp, on_disc_dir)
#' }
create_on_disc_matrix_from_10x_mtx <- function(mtx_fp, barcode_fp, features_fp, on_disc_dir) {
  # First, work with the .mtx file. Determine the number of rows containing comments.
  n_rows_with_comments <- 0
  repeat {
    curr_row <- utils::read.table(mtx_fp, nrows = 1, skip = n_rows_with_comments, header = FALSE, sep = "\n") %>% dplyr::pull()
    is_comment <- substr(curr_row, start = 1, stop = 1) == "%"
    if (!is_comment) {
      break()
    } else {
      n_rows_with_comments <- n_rows_with_comments + 1
    }
  }

  # Extract the number of rows, columns, and total data points in the expression matrix.
  metadata <- utils::read.table(mtx_fp, nrows = 1, skip = n_rows_with_comments, header = FALSE)
  n_genes <- metadata %>% dplyr::pull(1)
  n_cells <- metadata %>% dplyr::pull(2)
  n_data_points <- metadata %>% dplyr::pull(3)

  # Grab the cell barcodes, gene_ids, and gene_names.
  cell_barcodes <- invisible(readr::read_tsv(file = barcode_fp, col_types = "c", col_names = FALSE) %>% dplyr::pull())
  all_genes <- invisible(readr::read_tsv(file = features_fp, col_types = "ccc", col_names = c("gene_id", "gene_name", "feature_type")))
  gene_ids <- all_genes %>% dplyr::pull("gene_id")
  gene_names <- all_genes %>% dplyr::pull("gene_name")

  # Initialize the h5 file on disk
  h5_loc <- create_h5_file_on_disk(on_disc_dir, n_genes, n_cells, n_data_points, cell_barcodes, gene_ids, gene_names)
  rm(cell_barcodes, all_genes, gene_ids, gene_names); invisible(gc())

  # Load expression data into matrix in compressed sparse column format; output from this function the number of nonzero entries in each row.
  convert_mtx_to_csc(mtx_fp, h5_loc, n_rows_with_comments, chunk_size = 1e7)

  # Finally, transpose the CSC matrix to create a CSR matrix
  transpose_on_disc_csc_matrix(h5_loc, cell_chunk_size = 5000)

  # Return the newly created on_disc_matrix object
  ret <- on_disc_matrix(h5_file = h5_loc)
  return(ret)
}
