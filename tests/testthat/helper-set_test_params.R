cat("Running test setup script.\n")

# First, define the test type "small" or "large." Also set the temp_test_dir.
test_type <- "small"
temp_test_dir <- tempdir()

# Initialize the parameters.
if (test_type == "small") {
  n_datasets <- 3
  n_reps <- 2
  n_row <- NULL
  n_col <- NULL
} else {
  n_datasets <- 15
  n_reps <- 10
  n_row <- NULL
  n_col <- NULL
}

# Create synthetic data; save the corresponding list of r matrices in r_mats.
r_mats_plus_data <- create_synthetic_data(n_datasets = n_datasets, simulated_data_dir = temp_test_dir, n_row = n_row, n_col = n_col)
r_mats <- lapply(r_mats_plus_data, function(l) l$r_matrix )

# Initialize the ondisc_covariate_matrix objects
cov_odms <- get_metadata_odm_list(r_mats, 1L, temp_test_dir)

# initialize the multimodal_odm corresponding to matrix 1
n_cells_m1 <- ncol(r_mats[[1]])
new_modality <- create_synthetic_data(n_datasets = 1, simulated_data_dir = temp_test_dir, n_row = n_row, n_col = n_cells_m1, idx_start = length(r_mats) + 1L)
new_modality_metadata_odm <- get_metadata_odm_list(list(new_modality[[1]]$r_matrix), length(r_mats) + 1L, temp_test_dir)
multimodal_mat <- multimodal_ondisc_matrix(list(modality_1 = cov_odms[[1]],
                                                modality_2 = new_modality_metadata_odm[[1]]))

# write r_mats to .h5 files using create_ondisc_matrix_from_R_matrix
cov_odms_from_memory <- list()
for (i in seq(1, n_datasets)) {
  r_matrix <- r_mats_plus_data[[i]]$r_matrix
  barcodes <- r_mats_plus_data[[i]]$cell_barcodes
  features_df <- r_mats_plus_data[[i]]$features_df
  odm <- create_ondisc_matrix_from_R_matrix(r_matrix = r_matrix,
                                     barcodes = barcodes,
                                     features_df = features_df,
                                     on_disk_dir = temp_test_dir)
  cov_odms_from_memory[[i]] <- odm
}
