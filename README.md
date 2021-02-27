
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ondisc` <img src="man/figures/hex.png" align="right" alt="" width="180" />

<!-- badges: start -->

[![R build
status](https://travis-ci.com/timothy-barry/ondisc.svg?branch=main)](https://travis-ci.com/timothy-barry/ondisc)
[![Codecov test
coverage](https://codecov.io/gh/Timothy-Barry/ondisc/branch/main/graph/badge.svg)](https://codecov.io/gh/Timothy-Barry/ondisc?branch=main)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

Single-cell datasets are growing in size rapidly, posing a challenge for
biology researchers. `ondisc` (short for “on-disk single-cell”) makes
computing on large-scale single-cell data **FUN**:

  - **Fast**: `ondisc` is supported by several novel, highly efficient
    algorithms and data structures. All low-level code is written in C++
    or C for maximum performance.
  - **Universal**: `ondisc` runs on all platforms, from laptops to
    supercomputers. `ondisc` works seamlessly when the size of the data
    exceeds the amount of available memory.
  - **Ntuitive**: `ondisc` is easy to use. `ondisc` leverages ideas from
    functional programming, making it similar to the popular
    [tidyverse](https://www.tidyverse.org) suite of packages.

Take a look at the
[tutorials](https://timothy-barry.github.io/ondisc/articles/tutorial_odm_class.html)
on the [package
website](https://timothy-barry.github.io/ondisc/index.html) to get
up-and-running quickly.

## Installation

You can install the development version from GitHub with:

``` r
install.packages("devtools")
devtools::install_github("timothy-barry/ondisc")
```
