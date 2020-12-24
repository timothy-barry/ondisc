# This is a setup script. Its purpose is to ensure the synthetic data are in in place.

# Obtain test parameters given test type.
test_params <- get_test_parameters(get_test_type())

# Check if the synthetic data are already available; if so, no need to do anything
fs <- list.files(test_params$synthetic_data_dir)
data_types <- c("barcodes", "features", "matrix", "r_matrix")
file_extensions <- c(".tsv", ".tsv", ".mtx", ".rds")
synthetic_data_available <- sapply(1:length(data_types), function(i) {
  curr_fs <- paste0(data_types[i], "_", 1:test_params$n_datasets, file_extensions[i])
  all(curr_fs %in% fs)
  }) %>% all()

# If the synthetic data are not available, then create them
if (!synthetic_data_available) {
 create_synthetic_data(n_datasets = test_params$n_datasets, simulated_data_dir = test_params$synthetic_data_dir, seed = test_params$seed, n_row = test_params$n_row, n_col = test_params$n_col)
}
