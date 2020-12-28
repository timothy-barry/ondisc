# This script converts the .mtx file into an on_disc_matrix for use in subsequent tests.

cat("\nCreate on_disc_matrix from .mtx file.\n")

test_params <- get_test_parameters(get_test_type())
set.seed(test_params$seed)
for (i in 1:(test_params$n_datasets)) {
  if (test_params$n_datasets > 1) cat(paste0("Working on dataset ", i, ".\n"))
  fps <- get_simulation_data_fps(data_dir = test_params$synthetic_data_dir, idx = i)
  # check if on_disc_matrix already has been created; if so, delete that as well as the .h5 file
  if (file.exists(fps[["on_disc_matrix"]])) file.remove(fps[["on_disc_matrix"]]) %>% invisible()
  if (file.exists(fps[["on_disc_matrix_h5"]])) file.remove(fps[["on_disc_matrix_h5"]]) %>% invisible()
  m <- readRDS(fps[["r_matrix"]])
  n_col <- ncol(m)
  chunk_size <- sample(x = 2:(2 * n_col), size = 1)
  on_disc_obj <- create_on_disc_matrix_from_10x_mtx(mtx_fp = fps[["mtx"]], barcode_fp = fps[["barcodes"]], features_fp = fps[["features"]], on_disc_dir = test_params$synthetic_data_dir)
  saveRDS(object = on_disc_obj, file = fps[["on_disc_matrix"]])
}
