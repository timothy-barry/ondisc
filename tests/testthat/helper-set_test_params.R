cat("Running test setup script.\n")

# 1. Set the hyperparameters that define the test.
n_datasets <- 3
n_reps <- 2
n_row <- sample(x = seq(500, 1500), size = n_datasets, replace = TRUE)
n_col <- sample(x = seq(500, 1500), size = n_datasets, replace = TRUE)
logical_mat <- sample(c(TRUE, FALSE), size = n_datasets, replace = TRUE)


# 2. Create n_datasets synthetic matrices; write the mtx file, barcodes.tsv file, and features.tsv file to disk.
r_mats_plus_metadata <- lapply(seq(1, n_datasets), function(i) {
  create_synthetic_data(n_row = n_row[i], n_col = n_col[i],
                        logical_mat = logical_mat[i], write_as_mtx_to_disk = TRUE)
})
r_mats <- lapply(r_mats_plus_metadata, function(l) l$r_mat)


# 3. Initialize metadata_ondisc_matrix objects for each dataset (in a separate directory). Use variable chunk size.
cov_odms <- lapply(seq(1, n_datasets), function(i) {
  n_nonzero <- r_mats_plus_metadata[[i]]$n_nonzero
  bern <- as.logical(rbinom(1, 1, 0.5))
  chunk_size <- if (bern) n_nonzero + 1 else sample(x = seq(2, n_nonzero), 1)
  file_dir <- create_new_directory()
  metadata_odm <- create_ondisc_matrix_from_mtx(mtx_fp = r_mats_plus_metadata[[i]]$matrix_fp,
                                                barcodes_fp =  r_mats_plus_metadata[[i]]$barcodes_fp,
                                                features_fp = r_mats_plus_metadata[[i]]$features_fp,
                                                on_disk_dir = file_dir,
                                                n_lines_per_chunk = chunk_size,
                                                return_metadata_ondisc_matrix = TRUE,
                                                progress = FALSE)
})
