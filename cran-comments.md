## Test environments

* local Windows install, R 4.4.2

## R CMD check results

This submission updates `scMetaTraj` to version 0.1.2.

Key changes in this release:

* restored package-level GitHub/CRAN documentation (`README`, `NEWS`, vignette)
* exported the core workflow functions used for scoring, embedding, cluster
  profiling, and metabolic pseudotime inference
* clarified package metadata, URLs, and build exclusions

Local `R CMD check --as-cran` results:

* 0 ERROR
* 0 WARNING
* 1 NOTE

NOTE details:

* `unable to verify current time`
  This is a standard environment/system-time note on the local Windows check
  environment and does not reflect a package defect.

## Downstream dependencies

There are currently no known downstream dependencies for this package.
