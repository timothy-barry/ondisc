## Example, small, synthetic expression matrix in .mtx format

load_all()
package_dir <- system.file("extdata", package = "ondisc")
library(tidyverse)
library(Matrix)

set.seed(4)
m <- create_random_matrix(n_row = 300, n_col = 900)
save_random_matrix_as_10x(m = m, data_dir = package_dir, idx = 1)
