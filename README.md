
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ondisc` <img src="man/figures/hex.png" align="right" alt="" width="180" />

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Codecov test
coverage](https://codecov.io/gh/Timothy-Barry/ondisc/branch/main/graph/badge.svg)](https://codecov.io/gh/Timothy-Barry/ondisc?branch=main)
\[![R build
status](https://travis-ci.com/timothy-barry/ondisc.svg?branch=main)
<!-- badges: end -->

Single-cell datasets can be large and are growing in size as sequencing
costs drop. `ondisc` (short for **on-di**sk **s**ingle-**c**ell ) is an
R package that powers large-scale analyses of single-cell data. `ondisc`
has two main use cases: providing access to expression matrices too
large to fit in memory, and powering analyses on computer clusters (or
supercomputers) in which different processors operate on different
subsets of data. `ondisc` is functional (i.e., objects are persistent)
and efficient (i.e., algorithms are theoretically optimal), making the
package user-friendly. The vignette [Using
ondisc](https://timothy-barry.github.io/ondisc/articles/using_ondisc.html)
(10 min. read) gives a quick overview of the package.

## Installation

You can install development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("Timothy-Barry/ondisc")
```
