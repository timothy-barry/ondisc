library(testthat)
library(ondisc)

# Note: the change test type from "small" to "big", edit the function get_test_type. In the future, it may be possible to somehow change the test type using a simlar approach, like passing a command-line argument or sourcing from another script, but for now, this approach is safest.
test_check("ondisc")
