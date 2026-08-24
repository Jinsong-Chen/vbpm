# cran-comments

New submission. No earlier version of vbpm has been on CRAN.

## Test environments

* Windows 11 x64, R 4.5.1, `R CMD check --as-cran --run-donttest`
* Linux, x86_64-conda-linux-gnu, R 4.5.1

## R CMD check results

0 ERRORs, 0 WARNINGs.

The only NOTEs were the expected "New submission" and a local one reporting
that pandoc was not visible to the checking subprocess, so `README.md` and
`NEWS.md` were skipped; both render cleanly under pandoc 3.6.3 and the note
should not arise on your machines.

## Notes for the reviewer

* Possibly mis-spelled words in DESCRIPTION: `vbfa` and `vbmimic` are the
  package's estimator functions, `QUIC` is the sparse inverse-covariance
  algorithm behind the optional local-dependence branch, and `MIMIC` is the
  standard multiple-indicators multiple-causes acronym.

* The included `nlsy27` data are derived from the public National
  Longitudinal Survey of Youth 1997 (U.S. Bureau of Labor Statistics) and were
  previously distributed by the same maintainer in the CRAN package LAWBL
  (GPL-3). Provenance is documented in `?nlsy27`.
