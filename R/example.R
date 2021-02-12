if (FALSE) {
  f_name <- "example.h5"
  full_f <- paste0(tempdir(), "/", f_name)
  rhdf5::h5createFile(file = full_f)
  rhdf5::h5createDataset(file = full_f,
                  dataset = "example_vector",
                  dims = 10L,
                  storage.mode = "integer",
                  level = 0,
                  chunk = 10)
  rhdf5::h5read(file = full_f, name = "example_vector")
  read_data_h5(file_name_in = full_f, dataset_name_in = "example_vector")
  to_write <- as.integer(c(5,2,8))
  write_data_h5(file_name_in = full_f, dataset_name_in = "example_vector", buffer = to_write, start_pos = 4)
}
