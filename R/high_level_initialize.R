#' Create an ondisc_matrix from a .mtx file.
#'
#' Initializes an ondisc_matrix object from a .mtx file, a features.tsv file, and a barcodes.tsv file. Returns an ondisc_matrix R object along with cell-specific and feature-specific covariate matrices.
#'
#' @param mtx_fp file path to a .mtx file storing the expression data. The .mtx file can represent either an integer matrix or a logical (i.e., binary) matrix. If the .mtx file contains only two columns (after the initial three-column row of metadata), then the .mtx file is assumed to represent a logical matrix.
#' @param barcodes_fp file path to the .tsv file containing the cell barcodes.
#' @param features_fp file path to the features.tsv file. The first column (required) should contain the feature IDs (e.g., ENSG00000186092), and the second column (optional) should contain the human-readable feature names (e.g., OR4F5). Any subsequent columns are discarded.
#' @param n_gb_per_chunk (optional) amount of data (in GB) to process per chunk. Defaults to 4 GB.
#' @param on_disc_dir (optional) directory in which to store the on-disk portion of the ondisc_matrix. Defaults to the directory in which the .mtx file is located.
#' @param file_name (optional) name of the file in which to store the .h5 data on-disk. Defaults to ondisc_matrix_x.h5, where x is a unique integer starting at 1.
#'
#' @return
#' @export
create_ondisc_matrix_from_mtx <- function(mtx_fp, barcodes_fp, features_fp, n_gb_per_chunk = 4, on_disc_dir = NULL, file_name = NULL) {
  # Define "bag_of_variables" environment for storing args
  bag_of_variables <- new.env()

  # extract .mtx metadata
  n_rows_with_comments <- get_n_rows_with_comments_mtx(mtx_fp)
  mtx_metadata <- get_mtx_metadata(mtx_fp, n_rows_with_comments)
  bag_of_variables[[arguments_enum()$n_cells]] <- mtx_metadata$n_cells

  # extract features.tsv metadata; as a side-effect, if there are MT genes, put the locations of those genes into the bag_of_vars.
  features_metadata <- get_features_metadata(features_fp, bag_of_variables)

  # set the on_disc_dir, if necessary
  if (is.null(on_disc_dir)) on_disc_dir <- gsub(pattern = '/[^/]*$', replacement = "", x = mtx_fp)

  # Generate a name for the ondisc_matrix .h5 file, if necessary
  if (is.null(file_name)) h5_fp <- generate_on_disc_matrix_name(on_disc_dir)

  # Initialize the .h5 file on-disk (side-effect)
  initialize_h5_file_on_disk(h5_fp, mtx_metadata, features_metadata, barcodes_fp, features_fp)

  # Determine which covariataes to compute
  covariates <- map_inputs_to_covariates(mtx_metadata, features_metadata)

  # Obtain grammar
  grammar <- initialize_grammar()

  # Determine which terminal symbols to compute
  terminal_symbols <- lapply(unlist(covariates),
                             get_terminals_for_covariate, grammar = grammar) %>% unlist() %>% unique()


}
