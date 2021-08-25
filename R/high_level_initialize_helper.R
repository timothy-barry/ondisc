#' Get metadata for features.tsv file
#'
#' Gets metadata from a features.tsv file. As a side-effect, if MT genes are present, puts into the bag_of_variables a logical vector indicating the positions of those genes.
#'
#' @param features_fp file path to features.tsv file
#' @param bag_of_variables the bag of variables to which to add the mt_genes logical vector (if applicable)
#'
#' @return a list containing elements feature_names (logical), n_cols (integer), and whether MT genes are present (logical)
get_features_metadata <- function(features_fp, bag_of_variables) {
  first_row <- readr::read_tsv(file = features_fp, n_max = 1, col_names = FALSE, col_types = readr::cols())
  n_cols <- ncol(first_row)
  feature_names <- ncol(first_row) >= 2
  mt_genes_present <- FALSE
  if (feature_names) {
    gene_names <- read_given_column_of_tsv(col_idx = 2, n_cols = n_cols, tsv_file = features_fp)
    mt_genes <- grepl(pattern = "^MT-", x = gene_names)
    if (any(mt_genes)) {
      mt_genes_present <- TRUE
      bag_of_variables[[arguments_enum()$mt_gene_bool]] <- mt_genes
    }
  }
  return(list(feature_names = feature_names, n_cols = n_cols, mt_genes_present = mt_genes_present))
}


#' Read given column of tsv
#'
#' @param col_idx index of column to read
#' @param n_cols number of columns in file
#' @param tsv_file file path to .tsv file
#'
#' @return contents of the specified column in vector form
read_given_column_of_tsv <- function(col_idx, n_cols, tsv_file) {
  type_pattern <- c(rep("_", col_idx - 1), "c", rep("_", n_cols - col_idx)) %>% paste0(collapse = "")
  dplyr::pull(readr::read_tsv(file = tsv_file, col_names = FALSE, col_types = type_pattern))
}


#' Get mtx metadata
#'
#' @param mtx_fp path to the mtx file
#'
#' If total number of data points exceeds maximum integer value in R (.Machine$integer.max), throw error.
#'
#' @return a list containing (i) n_genes, (ii) n_cells, (iii) the number of
#'     data points (i.e., fraction of entries that are zero),
#'     (iv) (TRUE/FALSE) matrix is logical,
#'     (v) number of rows to skip before reading the data
#'     (vi) number of data points in each file,
get_mtx_metadata <- function(mtx_fp) {
  if (length(mtx_fp) == 1) {
    out <- .Call(`_ondisc_get_mtx_metadata`, mtx_fp)
  } else {
    # Handle overflow: R will type cast integer to double if there's an integer overflow
    cumulative_n_cells <- 0L
    cumulative_n_data_points <- 0L
    mtx_metadata_list <- vector(mode = "list", length = length(mtx_fp))
    n_rows_to_skip_list <- vector(mode = "integer", length = length(mtx_fp))
    for (i in 1:length(mtx_fp)) {
      mtx_metadata_list[[i]] <- .Call(`_ondisc_get_mtx_metadata`, mtx_fp[i])
      cumulative_n_cells <- cumulative_n_cells + mtx_metadata_list[[i]]$n_cells
      cumulative_n_data_points <- cumulative_n_data_points + mtx_metadata_list[[i]]$n_data_points
      n_rows_to_skip_list[[i]] <- mtx_metadata_list[[i]]$n_rows_to_skip
    }
    out <- mtx_metadata_list[[1]]
    out$n_cells <- cumulative_n_cells
    out$n_data_points <- cumulative_n_data_points
    out$n_rows_to_skip <- n_rows_to_skip_list
    out$n_cells_in_files <- sapply(mtx_metadata_list, function(l) l$n_cells)
  }
  return(out)
}


#' Get h5 full name by the dataset name
#'
#' @param h5_info a dataframe that contains the information of the  h5 file
#' @param name the dataset name of
#'
#' @return full name of the dataset
#'
#' @noRd
get_h5_full_name <- function(h5_info, name) {
  idx <- which(h5_info$name == name)
  return(paste(h5_info$group[idx], name, sep = "/"))
}


#' Get h5 barcodes
#'
#' @param h5_list a vector of paths to the h5 file
#' @param barcode_suffixes a vector of suffix that appended to the barcodes to for each input barcode_fp. If NULL, append file index.
#'
#' @return a list of barcodes for each h5 file
get_h5_barcodes <- function(h5_list, barcode_suffixes) {
  if (is.null(barcode_suffixes) && length(h5_list) > 1L) {
    barcode_suffixes <- seq(1,length(h5_list))
  }
  barcodes_list <- vector(mode = "list", length = length(h5_list))
  for (i in seq(1,length(h5_list))) {
    h5_info <- rhdf5::h5ls(h5_list[i])
    barcodes_name <- get_h5_full_name(h5_info, "barcodes")
    barcodes_list[[i]] <- rhdf5::h5read(h5_list[i], barcodes_name)
    if (!is.null(barcode_suffixes)) {
      barcodes_list[[i]] <- paste(barcodes_list[[i]], barcode_suffixes[i], sep = "_")
    }
  }
  return(barcodes_list)
}


#' Get h5 features
#'
#' @param h5_list a vector of paths to the h5 file
#'
#' @return features_df data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded. Gene names starting with "MT-" are assumed to be mitochondrial genes and will be used to compute the p_mito covariate.
get_h5_features <- function(h5_list) {
  h5_info <- rhdf5::h5ls(h5_list[1])
  genes_name <- get_h5_full_name(h5_info, "genes")
  id <- rhdf5::h5read(h5_list[1], genes_name)
  gene_names_name <- get_h5_full_name(h5_info, "gene_names")
  name <- rhdf5::h5read(h5_list[1], gene_names_name)
  features_df <- data.frame(id, name)
  return(features_df)
}


#' Get h5 cells metadata
#'
#' @param h5_list a vector of paths to the h5 file
#'
#' @return a list containing (i)n_cells, (ii) the number of
#'     data points (i.e., fraction of entries that are zero),
#'     (iii) number of data points in each file,
#'     (iv) (always FALSE) matrix is logical
#' @noRd
get_h5_cells_metadata <- function(h5_list) {
  cumulative_n_cells <- 0L
  cumulative_n_data_points <- 0L
  n_cells_in_files <- vector(mode = "integer", length = length(h5_list))
  for (i in seq(1,length(h5_list))) {
    h5_info <- rhdf5::h5ls(h5_list[i])
    shape_name <- get_h5_full_name(h5_info, "shape")
    shape <- rhdf5::h5read(h5_list[i], shape_name)
    cumulative_n_cells <- cumulative_n_cells + shape[2]
    n_cells_in_files[i] <- shape[2]
    idx <- which(h5_info$name == "data")
    n_data_points <- h5_info$dim[idx]
    cumulative_n_data_points <- cumulative_n_data_points + as.integer(n_data_points)
  }
  out <- list()
  out$n_cells <- cumulative_n_cells
  out$n_data_points <- cumulative_n_data_points
  out$n_cells_in_files <- n_cells_in_files
  out$is_logical <- FALSE
  return(out)
}


#' Get string
#'
#' @param barcodes_fp path to barcodes file(s)
#' @param features_fp path the features file
#' @param features_metadata the features_metadata list
#' @param barcode_suffixes a vector of suffix that appended to the barcodes to for each input barcode_fp. If NULL, append file index.
#'
#' @return list containing (i) cell barcodes, (ii) feature IDs, and (iii) feature names (NA_character_ if absent)
get_string_arrays <- function(barcodes_fp, features_fp, features_metadata, barcode_suffixes) {
  if (is.null(barcode_suffixes) && length(barcodes_fp) > 1L) {
    barcode_suffixes <- seq(1,length(barcodes_fp))
  }
  # barcodes first
  cell_barcodes <- vector(mode = "character")
  for (i in seq(1L, length(barcodes_fp))) {
    single_cell_barcodes <- dplyr::pull(readr::read_tsv(file = barcodes_fp[i], col_names = FALSE, col_types = "c"))
    if (!is.null(barcode_suffixes)) {
      single_cell_barcodes <- paste(single_cell_barcodes, barcode_suffixes[i], sep = "_")
    }
    # concatenate cell_barcodes in each file to a single barcodes vector
    cell_barcodes <- c(cell_barcodes, single_cell_barcodes)
  }
  # feature IDs and names next
  feature_ids <- read_given_column_of_tsv(col_idx = 1, n_cols = features_metadata$n_cols, tsv_file = features_fp)
  if (features_metadata$feature_names) {
    feature_names <- read_given_column_of_tsv(col_idx = 2, n_cols = features_metadata$n_cols, tsv_file = features_fp)
  } else {
    feature_names <- NA_character_
  }
  return(list(cell_barcodes = cell_barcodes, feature_ids = feature_ids, feature_names = feature_names))
}
