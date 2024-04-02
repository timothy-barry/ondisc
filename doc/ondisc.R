## ----eval=FALSE---------------------------------------------------------------
#  # install.packages("BiocManager"); install.packages("devtools")
#  BiocManager::install("Rhdf5lib", type = "source") # Rhdf5lib
#  devtools::install_github("timothy-barry/ondisc") # ondisc
#  devtools::install_github("katsevich-lab/sceptredata") # sceptredata

## -----------------------------------------------------------------------------
library(ondisc)
library(sceptredata)

## -----------------------------------------------------------------------------
directories_to_load <- paste0(
  system.file("extdata", package = "sceptredata"), 
  "/highmoi_example/gem_group_", 1:2
)
directories_to_load # file paths to the example data on your computer

## -----------------------------------------------------------------------------
list.files(directories_to_load[1])

## -----------------------------------------------------------------------------
list.files(directories_to_load[2])

## ----results="hide"-----------------------------------------------------------
temp_dir <- tempdir()

out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = temp_dir 
)

## -----------------------------------------------------------------------------
list.files(temp_dir, pattern = "*.odm")

## -----------------------------------------------------------------------------
gene_odm <- out_list[["gene"]]

## -----------------------------------------------------------------------------
gene_odm

## -----------------------------------------------------------------------------
n_features <- nrow(gene_odm)
n_features

n_cells <- ncol(gene_odm)
n_cells

## -----------------------------------------------------------------------------
feature_ids <- rownames(gene_odm)
head(feature_ids)

## -----------------------------------------------------------------------------
expression_vector <- gene_odm[2,]
head(expression_vector)

expression_vector <- gene_odm["ENSG00000117222",]
head(expression_vector)

## -----------------------------------------------------------------------------
object.size(gene_odm) |> format(units = "Kb")

## -----------------------------------------------------------------------------
set.seed(4)
example_data <- write_example_cellranger_dataset(
  n_features = c(500, 50, 20),
  n_cells = 10000,
  n_batch = 3,
  modalities = c("gene", "grna", "protein"),
  directory_to_write = temp_dir ,
  p_set_col_zero = 0
)

## -----------------------------------------------------------------------------
directories_to_load <- list.files(
  temp_dir,
  pattern = "batch_",
  full.names = TRUE
)
directories_to_load

## -----------------------------------------------------------------------------
list.files(directories_to_load[1])

## ----results="hide"-----------------------------------------------------------
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = temp_dir
)

## -----------------------------------------------------------------------------
names(out_list)

## -----------------------------------------------------------------------------
list.files(temp_dir, pattern = "*.odm")

## -----------------------------------------------------------------------------
cellwise_covariates <- out_list[["cellwise_covariates"]]
head(cellwise_covariates)

## -----------------------------------------------------------------------------
rm(list = ls()) # delete all variables
temp_dir <- tempdir()
gene_odm <- initialize_odm_from_backing_file(
  paste0(temp_dir, "/gene.odm")
)
gene_odm

## -----------------------------------------------------------------------------
data(lowmoi_example_data)
gene_mat <- lowmoi_example_data$response_matrix

## -----------------------------------------------------------------------------
file_to_write <- paste0(temp_dir, "/gene.odm")
gene_odm <- create_odm_from_r_matrix(
  mat = gene_mat,
  file_to_write = file_to_write
)

## -----------------------------------------------------------------------------
gene_odm

## -----------------------------------------------------------------------------
library(sessioninfo); session_info()

