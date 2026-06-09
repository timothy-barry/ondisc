This is a resubmission of a previously archived package, ondisc. The package was archived because CRAN issues were not fixed in time.

The 1.3.x updates remove the GitHub packages `sceptre` and `sceptredata` as dependencies and reduce the runtime of the vignette build, examples, and tests.

Version 1.3.1 fixed an AddressSanitizer heap-buffer-overflow in `write_to_csr()` triggered when importing Cell Ranger data with collapsed gRNA vector counts.

Version 1.3.2 further reduces test runtime by shrinking randomized test matrices in response to CRAN feedback.

Version 1.3.3 ensures that Cell Ranger import examples run with a single data.table thread, addressing a Debian NOTE about CPU time exceeding elapsed time in examples.

Version 1.3.4 addresses CRAN manual review comments by using CRAN-compatible quotation marks and DOI formatting in DESCRIPTION, replacing unsuppressible progress output with message(), and removing a global-environment cleanup line from the vignette.
