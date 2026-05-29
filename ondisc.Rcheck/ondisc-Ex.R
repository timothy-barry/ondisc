pkgname <- "ondisc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ondisc-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ondisc')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("create_odm_from_cellranger")
### * create_odm_from_cellranger

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_odm_from_cellranger
### Title: Create 'odm' object from Cell Ranger
### Aliases: create_odm_from_cellranger

### ** Examples

library(sceptredata)
directories_to_load <- paste0(
 system.file("extdata", package = "sceptredata"),
 "/highmoi_example/gem_group_", c(1, 2)
)
directory_to_write <- tempdir()
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = directory_to_write,
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_odm_from_cellranger", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("create_odm_from_r_matrix")
### * create_odm_from_r_matrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create_odm_from_r_matrix
### Title: Create 'odm' object from R matrix
### Aliases: create_odm_from_r_matrix

### ** Examples

library(sceptredata)
data(lowmoi_example_data)
gene_matrix <- lowmoi_example_data$response_matrix
file_to_write <- paste0(tempdir(), "/gene.odm")
odm_object <- create_odm_from_r_matrix(
  mat = gene_matrix,
  file_to_write = file_to_write
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create_odm_from_r_matrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dim-odm-method")
### * dim-odm-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dim,odm-method
### Title: Return the number of columns and rows of an 'odm' object
### Aliases: dim,odm-method

### ** Examples

library(sceptredata)
directories_to_load <- paste0(
 system.file("extdata", package = "sceptredata"),
 "/highmoi_example/gem_group_", c(1, 2)
)
directory_to_write <- tempdir()
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = directory_to_write,
)
gene_odm <- out_list$gene
# return the dimension, number of rows, and number of columns
dim(gene_odm)
nrow(gene_odm)
ncol(gene_odm)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dim-odm-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("initialize_odm_from_backing_file")
### * initialize_odm_from_backing_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: initialize_odm_from_backing_file
### Title: Initialize an 'odm' object
### Aliases: initialize_odm_from_backing_file

### ** Examples

library(sceptredata)
directories_to_load <- paste0(
 system.file("extdata", package = "sceptredata"),
 "/highmoi_example/gem_group_", c(1, 2)
)
directory_to_write <- tempdir()
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = directory_to_write,
)
gene_odm <- out_list$gene
gene_odm

# delete the gene_odm object
rm(gene_odm)

# reinitialize the gene_odm object
gene_odm <- initialize_odm_from_backing_file(
  paste0(tempdir(), "/gene.odm")
)
gene_odm



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("initialize_odm_from_backing_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ondisc-package")
### * ondisc-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ondisc-package
### Title: ondisc: Algorithms and data structures for large single-cell
###   expression matrices
### Aliases: ondisc ondisc-package

### ** Examples

# initialize odm objects from Cell Ranger output; also, compute the cellwise covariates
library(sceptredata)
directories_to_load <- paste0(
 system.file("extdata", package = "sceptredata"),
 "/highmoi_example/gem_group_", c(1, 2)
)
directory_to_write <- tempdir()
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = directory_to_write,
)

# extract the odm corresponding to the gene modality
gene_odm <- out_list$gene
gene_odm

# obtain dimension information
dim(gene_odm)
nrow(gene_odm)
ncol(gene_odm)

# obtain rownames (i.e., the feature IDs)
rownames(gene_odm) |> head()

# extract row into memory, first by integer and then by string
expression_vector_1 <- gene_odm[10,]
expression_vector_2 <- gene_odm["ENSG00000135046",]

# delete the gene_odm object
rm(gene_odm)

# reinitialize the gene_odm object
gene_odm <- initialize_odm_from_backing_file(
  paste0(tempdir(), "/gene.odm")
)
gene_odm



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ondisc-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rownames-odm-method")
### * rownames-odm-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dimnames,odm-method
### Title: Return the rownames of an 'odm' object
### Aliases: dimnames,odm-method rownames

### ** Examples

library(sceptredata)
directories_to_load <- paste0(
 system.file("extdata", package = "sceptredata"),
 "/highmoi_example/gem_group_", c(1, 2)
)
directory_to_write <- tempdir()
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = directory_to_write,
)
gene_odm <- out_list$gene
# return the rownames
rownames(gene_odm) |> head()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rownames-odm-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sub-odm-ANY-missing-missing-method")
### * sub-odm-ANY-missing-missing-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: [,odm,ANY,missing,missing-method
### Title: Load a row of an 'odm' object into memory
### Aliases: [,odm,ANY,missing,missing-method

### ** Examples

library(sceptredata)
directories_to_load <- paste0(
 system.file("extdata", package = "sceptredata"),
 "/highmoi_example/gem_group_", c(1, 2)
)
directory_to_write <- tempdir()
out_list <- create_odm_from_cellranger(
  directories_to_load = directories_to_load,
  directory_to_write = directory_to_write,
)
gene_odm <- out_list$gene
# extract rows into memory by index and ID
v1 <- gene_odm[10L,]
v2 <- gene_odm["ENSG00000173825",]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sub-odm-ANY-missing-missing-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("write_example_cellranger_dataset")
### * write_example_cellranger_dataset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: write_example_cellranger_dataset
### Title: Write example Cell Ranger dataset
### Aliases: write_example_cellranger_dataset

### ** Examples

set.seed(4)
n_features <- c(1000, 40, 400)
modalities <- c("gene", "protein", "grna")
n_cells <- 10000
n_batch <- 2
directory_to_write <- tempdir()
p_set_col_zero <- 0
out <- write_example_cellranger_dataset(
  n_features = n_features,
  n_cells = n_cells,
  n_batch = n_batch,
  modalities = modalities,
  directory_to_write = directory_to_write,
  p_set_col_zero = p_set_col_zero
)

# directories written to directory_to_write
fs <- list.files(directory_to_write, pattern = "batch*", full.names = TRUE)
# files contained within the directories
list.files(fs[1])
list.files(fs[2])



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("write_example_cellranger_dataset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
