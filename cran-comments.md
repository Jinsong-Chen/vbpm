# cran-comments

New submission. No earlier version of vbpm has been on CRAN. This is a
resubmission of 0.9.0, which failed the Debian incoming check.

## The 0.9.0 failure and the fix

One test compared `vbfa()` output to a stored reference with
`expect_identical()`. The estimator runs through compiled BLAS/LAPACK, so the
last one or two ulps are not portable; Debian disagreed in the sixteenth
significant digit on four values (797 tests passed). The test now asserts
agreement at a strict relative `tolerance = 1e-8`, which remains sensitive to
material regression. Bitwise identity is no longer asserted anywhere, being a
property of a single build. No estimator or user-visible behaviour changed.

## Test environments

* Windows 11 x64, R 4.5.1, `R CMD check --as-cran --run-donttest`
* Linux, x86_64-conda-linux-gnu, R 4.5.1

## R CMD check results

0 ERRORs, 0 WARNINGs, 1 NOTE (new submission).

## Notes for the reviewer

* Possibly mis-spelled words in DESCRIPTION. `Guo`, `Jin`, `Yan` and `Zhang`
  are author surnames in the cited references. `vbfa` and `vbmimic` are the
  package's estimator functions. `QUIC` is the sparse inverse-covariance
  algorithm behind the optional local-dependence branch, `MIMIC` is the
  standard multiple-indicators multiple-causes acronym, and `Variational` and
  `bifactor` are the usual spellings of those terms in this literature.

* The included `nlsy27` data are derived from the public National
  Longitudinal Survey of Youth 1997 (U.S. Bureau of Labor Statistics) and were
  previously distributed by the same maintainer in the CRAN package LAWBL
  (GPL-3). Provenance is documented in `?nlsy27`.
