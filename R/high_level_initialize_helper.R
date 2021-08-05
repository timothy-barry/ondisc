#' Get metadata for features.tsv file
#'
#' Gets metadata from a features.tsv file. As a side-effect, if MT genes are present, puts into the bag_of_variables a logical vector indicating the positions of those genes.
#'
#' @param features_fp file path to features.tsv file
#' @param bag_of_variables the bag of variables to which to add the mt_genes logical vector (if applicable)
#'
#' @return a list containing elements feature_names (logical), n_cols (integer), and whether MT genes are present (logical)
#' @noRd
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


#' Generate on disc_matrix_name
#'
#' Generates the name of an on_disc_matrix object given a directory.
#' This function searches for files named on_disc_matrix_id.h5 in the specified directory.
#' If none exists, it returns on_disc_matrix_1.h5. Else, it returns n_disc_matrix_id.h5
#' with a unique integer in place of id.
#'
#' @param on_disc_dir directory in which to store the on_disc_matrix.
#' @return a new name for an on_disc_matrix.
#' @noRd
generate_on_disc_matrix_name <- function(on_disc_dir) {
  # Only list ondisc_matrix_<id>.h5 files
  base_name <- "ondisc_matrix_"
  fs <- list.files(on_disc_dir, pattern = paste0(base_name, "[0-9]+\\.h5"))
  if (length(fs) == 0) {
    name <- paste0(base_name, "1.h5")
  } else {
    ints_in_use <- gsub(pattern = paste0(base_name, "(\\d+)\\.h5"), replacement = "\\1", x = fs) %>% as.integer()
    new_int <- max(ints_in_use) + 1
    name <- paste0(base_name, new_int, ".h5")
  }
  return(name)
}


#' Read given column of tsv
#'
#' @param col_idx index of column to read
#' @param n_cols number of columns in file
#' @param tsv_file file path to .tsv file
#'
#' @return contents of the specified column in vector form
#' @noRd
read_given_column_of_tsv <- function(col_idx, n_cols, tsv_file, progress = FALSE) {
  type_pattern <- c(rep("_", col_idx - 1), "c", rep("_", n_cols - col_idx)) %>% paste0(collapse = "")
  dplyr::pull(readr::read_tsv(file = tsv_file, col_names = FALSE, col_types = type_pattern, progress = progress))
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


#' Verify fp
#'
#' Verifies that a file path to be used as the storage location for an ondisc matrix is acceptable.
#'
#' There are three cases:
#' - case 1. The directory does not exist; return TRUE, and as a side-effect, create the empty directory.
#' - case 2. The directory does exist and is empty; return TRUE.
#' - case 3. The directory does exist and is nonempty; return FALSE.
#'
#' @param fp a file path
#'
#' @return boolean; if TRUE, no problem; if FALSE, problem.
verify_fp <- function(fp) {
  if (!dir.exists(fp)) {
    # side effect: create directory
    dir.create(path = fp, recursive = TRUE)
    OK <- TRUE
  } else {
    empty_dir <- length(list.files(fp)) == 0
    if (empty_dir) {
      OK <- TRUE
    } else {
      OK <- FALSE
    }
  }
  return(OK)
}


#' Get metadata for features_df,a data frame giving the names of the features.
#'
#' Gets metadata from a features data frame features_df.
#'
#' @param features_df a data frame giving the names of the features.
#' @param bag_of_variables the bag of variables to which to add the mt_genes logical vector (if applicable)
#'
#' @return a list containing elements feature_names (logical), n_cols (integer), and whether MT genes are present (logical)
#' @noRd
get_features_metadata_from_table <- function(features_df, bag_of_variables = NULL) {
  n_cols <- ncol(features_df)
  feature_names <- n_cols >= 2
  mt_genes_present <- FALSE
  if (feature_names) {
    # Assume the second column is always feature_name. Or we can extract by the col name
    gene_names <- dplyr::pull(features_df, 2)
    mt_genes <- grepl(pattern = "^MT-", x = gene_names)
    if (any(mt_genes)) {
      mt_genes_present <- TRUE
      if (!is.null(bag_of_variables)) {
        bag_of_variables[[arguments_enum()$mt_gene_bool]] <- mt_genes
      }
    }
  }
  return(list(feature_names = feature_names, n_cols = n_cols, mt_genes_present = mt_genes_present))
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
  return(paste(h5_info$group[idx], name,sep = "/"))
}


#' Get h5 barcodes
#'
#' @param h5_list a vector of paths to the h5 file
#'
#' @return a list of barcodes for each h5 file
#' @noRd
get_h5_barcodes <- function(h5_list) {
  barcodes_list <- vector(mode = "list", length = length(h5_list))
  for (i in seq(1,length(h5_list))) {
    h5_info <- rhdf5::h5ls(h5_list[i])
    barcodes_name <- get_h5_full_name(h5_info, "barcodes")
    barcodes_list[[i]] <- rhdf5::h5read(h5_list[i], barcodes_name)
  }
  return(barcodes_list)
}


#' Get h5 features
#'
#' @param h5_list a vector of paths to the h5 file
#'
#' @return features_df data frame giving the names of the features. The first column (required) contains the feature IDs (e.g., ENSG00000186092), and the second column (optional) contains the human-readable feature names (e.g., OR4F5). Subsequent columns are discarded. Gene names starting with "MT-" are assumed to be mitochondrial genes and will be used to compute the p_mito covariate.
#' @noRd
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
#' @return a list containing (i) n_genes, (ii) n_cells, (iii) the number of
#'     data points (i.e., fraction of entries that are zero)
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
  return(out)
}
