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
r_mats <- create_synthetic_data(n_datasets = n_datasets, simulated_data_dir = temp_test_dir, n_row = n_row, n_col = n_col)

# Initialize the ondisc_covariate_matrix objects
cov_odms <- vector(mode = "list", length = n_datasets)
for (i in seq(1, n_datasets)) {
  fps <- get_simulation_data_fps(data_dir = temp_test_dir, idx = i)
  # check if on_disc_matrix already has been created; if so, delete that as well as the .h5 file
  if (file.exists(fps[["on_disc_matrix_h5"]])) file.remove(fps[["on_disc_matrix_h5"]]) %>% invisible() # h5 file
  m <- r_mats[[i]]
  n_data_points <- length(m@x)
  if (stats::rbinom(1, 1, 0.5)) {
    chunk_size <- sample(x = seq(2, n_data_points - 1), size = 1) # choose chunk size less than n_data_points
  } else {
    chunk_size <- sample(x = seq(n_data_points + 1, 2 * n_data_points), size = 1)
  }
  cov_odm_obj <- create_ondisc_matrix_from_mtx(mtx_fp = fps[["mtx"]],
                                               barcodes_fp = fps[["barcodes"]],
                                               features_fp = fps[["features"]],
                                               n_lines_per_chunk = chunk_size,
                                               on_disc_dir = temp_test_dir,
                                               return_covariate_ondisc_matrix = TRUE)
  cov_odms[[i]] <- cov_odm_obj
}
