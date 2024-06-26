% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ondisc.R
\docType{package}
\name{ondisc-package}
\alias{ondisc}
\alias{ondisc-package}
\title{ondisc: Algorithms and data structures for large single-cell expression matrices}
\description{
Single-cell datasets are growing in size, posing challenges as well as opportunities for genomics researchers. `ondisc` is an R package that facilitates analysis of large-scale single-cell data out-of-core on a laptop or distributed across tens to hundreds processors on a cluster or cloud. In both of these settings, `ondisc` requires only a few gigabytes of memory, even if the input data are tens of gigabytes in size. `ondisc` mainly is oriented toward single-cell CRISPR screen analysis, but ondisc also can be used for single-cell differential expression and single-cell co-expression analyses. ondisc is powered by several new, efficient algorithms for manipulating and querying large, sparse expression matrices.
}
\examples{
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
}
\seealso{
Useful links:
\itemize{
  \item \url{https://timothy-barry.github.io/ondisc/}
  \item \url{https://timothy-barry.github.io/sceptre-book/}
  \item Report bugs at \url{https://github.com/timothy-barry/ondisc/issues}
}

}
\author{
\strong{Maintainer}: Timothy Barry \email{tbarry@hsph.harvard.edu} (\href{https://orcid.org/0000-0002-4356-627X}{ORCID})

Authors:
\itemize{
  \item Eugene Katsevich \email{ekatsevi@wharton.upenn.edu} [thesis advisor]
}

Other contributors:
\itemize{
  \item Songcheng Dai [contributor]
  \item Yixuan Qiu [contributor]
}

}
