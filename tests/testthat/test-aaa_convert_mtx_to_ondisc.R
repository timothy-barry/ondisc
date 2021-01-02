# This script converts the .mtx file into an on_disc_matrix for use in subsequent tests. No testthat statement actually is evaluated here; the purpose simply is to verify that create_on_disc_matrix_from_10x_mtx executes without error.

cat("\nCreate on_disc_matrix from .mtx file.\n")

test_that("Create on_disc_matrices", {
  for (i in 1:(n_datasets)) {
    if (n_datasets > 1) cat(paste0("Working on dataset ", i, ".\n"))
    fps <- get_simulation_data_fps(data_dir = temp_test_dir, idx = i)
    # check if on_disc_matrix already has been created; if so, delete that as well as the .h5 file
    if (file.exists(fps[["on_disc_matrix"]])) file.remove(fps[["on_disc_matrix"]]) %>% invisible()
    if (file.exists(fps[["on_disc_matrix_h5"]])) file.remove(fps[["on_disc_matrix_h5"]]) %>% invisible()
    m <- readRDS(fps[["r_matrix"]])
    n_col <- ncol(m)
    if (i %% 2 == 0) { # if i is even
      chunk_size <- sample(x = 2:((n_col) - 1), size = 1) # pick a chunk size less than n_col
    } else { # if i is odd
      chunk_size <- sample(x = (n_col + 1):(2 * n_col), size = 1) # pick a chunk size greater than n_col
    }
    on_disc_obj <- create_on_disc_matrix_from_10x_mtx(mtx_fp = fps[["mtx"]], barcode_fp = fps[["barcodes"]], features_fp = fps[["features"]], on_disc_dir = temp_test_dir, n_cells_to_process_at_time = chunk_size)
    saveRDS(object = on_disc_obj, file = fps[["on_disc_matrix"]])
  }
  expect_true(TRUE)
})
