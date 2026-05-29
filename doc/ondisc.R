# # install.packages("BiocManager"); install.packages("devtools")
# BiocManager::install("Rhdf5lib", type = "source") # Rhdf5lib
# devtools::install_github("timothy-barry/ondisc") # ondisc

library(ondisc)

directories_to_load <- paste0(
  system.file("extdata", "highmoi_example", package = "ondisc"), 
  "/gem_group_", 1:2
)
directories_to_load # file paths to the example data on your computer

list.files(directories_to_load[1])

list.files(directories_to_load[2])

temp_dir <- tempdir()

out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = temp_dir 
)

list.files(temp_dir, pattern = "*.odm")

gene_odm <- out_list[["gene"]]

gene_odm

n_features <- nrow(gene_odm)
n_features

n_cells <- ncol(gene_odm)
n_cells

feature_ids <- rownames(gene_odm)
head(feature_ids)

expression_vector <- gene_odm[2,]
head(expression_vector)

expression_vector <- gene_odm[rownames(gene_odm)[2],]
head(expression_vector)

object.size(gene_odm) |> format(units = "Kb")

set.seed(4)
example_data <- write_example_cellranger_dataset(
  n_features = c(100, 20, 10),
  n_cells = 500,
  n_batch = 3,
  modalities = c("gene", "grna", "protein"),
  directory_to_write = temp_dir ,
  p_set_col_zero = 0
)

directories_to_load <- list.files(
  temp_dir,
  pattern = "batch_",
  full.names = TRUE
)
directories_to_load

list.files(directories_to_load[1])

out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = temp_dir
)

names(out_list)

list.files(temp_dir, pattern = "*.odm")

cellwise_covariates <- out_list[["cellwise_covariates"]]
head(cellwise_covariates)

rm(list = ls()) # delete all variables
temp_dir <- tempdir()
gene_odm <- initialize_odm_from_backing_file(
  paste0(temp_dir, "/gene.odm")
)
gene_odm

gene_mat <- matrix(
  c(1L, 0L, 3L, 0L, 2L, 0L, 4L, 5L, 0L, 6L, 0L, 7L),
  nrow = 3L,
  dimnames = list(paste0("gene_", 1:3), paste0("cell_", 1:4))
)

file_to_write <- paste0(temp_dir, "/gene.odm")
gene_odm <- create_odm_from_r_matrix(
  mat = gene_mat,
  file_to_write = file_to_write,
  chunk_size = 5L
)

gene_odm

library(sessioninfo); session_info()
