## Example, small, synthetic expression matrix in .mtx format
## Note: it is assumed that the PBMC data have been downloaded from 10X; see the "Using ondisc" vignette for code to do this.

library(tidyverse)
data_dir <- "/Users/timbarry/Box/onDisc_all/onDisc_offsite/raw_data" # change this filepath to reproduce!
package_dir <- "/Users/timbarry/Box/onDisc_all/ondisc/" # change this filepath as well!

# Load package via load_all (to access non-exported functions)
load_all(path = package_dir, helpers = FALSE)
ext_data_dir <- paste0(package_dir, "inst/extdata")

# Load barcodes and gene ids, names
sub_data_dir <- paste0(data_dir, "/filtered_feature_bc_matrix")
all_barcodes <- read_tsv(file = paste0(sub_data_dir, "/barcodes.tsv"), col_types = "c", col_names = FALSE) %>% pull()
all_gene_ids_and_names <- read_tsv(file = paste0(sub_data_dir, "/features.tsv"), col_types = "ccc", col_names = c("gene_id", "gene_name", "feature_type"))

# generate data
n_row <- 300
n_col <- 900
set.seed(4)
my_cell_barcodes <- sample(x = all_barcodes, size = n_col, replace = FALSE)
samped_gene_ids_and_names <- all_gene_ids_and_names %>% sample_n(n_row)
my_gene_ids <- samped_gene_ids_and_names %>% pull(gene_id)
my_gene_names <- samped_gene_ids_and_names %>% pull(gene_name)

m <- create_random_matrix(n_row = 300, n_col = 900)
save_random_matrix_as_10x(m = m, data_dir = ext_data_dir, idx = NULL, cell_barcodes = my_cell_barcodes, gene_names = my_gene_names, gene_ids = my_gene_ids, save_r_matrix = FALSE)
