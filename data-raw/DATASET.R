## Example, small, synthetic expression matrix in .mtx format

package_dir <- "/Users/timbarry/Box/onDisc_all/ondisc/"
library(tidyverse)
library(Matrix)
load_all()

set.seed(4)
m <- create_random_matrix(n_row = 300, n_col = 900)
save_random_matrix_as_10x(m = m, data_dir = paste0(package_dir, "inst/extdata"), idx = 1)
