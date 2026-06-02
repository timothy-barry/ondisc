This is a resubmission of a previously archived package, ondisc. The package was archived because CRAN issues were not fixed in time.

The 1.3.x updates remove the GitHub packages `sceptre` and `sceptredata` as dependencies and reduce the runtime of the vignette build, examples, and tests.

Version 1.3.1 fixed an AddressSanitizer heap-buffer-overflow in `write_to_csr()` triggered when importing Cell Ranger data with collapsed gRNA vector counts.

Version 1.3.2 further reduces test runtime by shrinking randomized test matrices in response to CRAN feedback.
