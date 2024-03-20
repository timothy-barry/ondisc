library(testthat)
library(ondisc)

my_seed <- as.integer(Sys.time() |> format("%H%M%S"))
test_check("ondisc")
