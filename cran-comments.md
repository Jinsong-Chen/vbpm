# cran-comments

## Submission summary

This is a new submission of **vbpm** (Variational Bayes Psychometric Models).

The package provides variational Bayes estimation for partially confirmatory
psychometric measurement models: `vbfa()` (factor analysis with spike-and-slab
selection of unspecified loadings, an orthogonal bifactor option, and optional
sparse residual/local-dependence estimation) and `vbmimic()` (the MIMIC
extension, regularizing both the measurement and the structural part), together
with fit statistics, factor-number selection, and two data simulators.

## Test environments

* Windows 11 x64, R 4.5.1 (local),
  `R CMD check --as-cran --run-donttest` (PDF manual included)

## R CMD check results

Status: 2 NOTEs, 0 ERRORs, 0 WARNINGs.

* The sandbox could not verify the current time.
* The optional package `V8` was unavailable, so HTML-manual math rendering was
  not checked.

Compiled installation, examples, all tests, Rd validation, vignette rebuilding,
and indexed PDF-manual generation passed.

## Notes for the reviewer

* **Compiled code.** `src/` contains a QUIC (quadratic approximation for sparse
  inverse covariance) solver used only by the optional local-dependence branch,
  built with Rcpp/RcppArmadillo. `src/Makevars` links LAPACK/BLAS explicitly.

* **Included data.** `nlsy27` (3,458 x 27 responses) is derived from the public
  National Longitudinal Survey of Youth 1997 (Bureau of Labor Statistics, U.S.
  Department of Labor) and was previously distributed by the same maintainer in
  the CRAN package **LAWBL** (GPL-3). Provenance is documented in
  `?nlsy27`.

* **Relationship to LAWBL.** This package supersedes the variational-estimation
  side of the maintainer's earlier LAWBL package. `sim_lvm()` is carried over
  from LAWBL with corrections (notably a missingness mask that previously
  recycled across the covariate block) and extensions.

* **Deliberate `NA` values.** For local-dependence fits, `ELBO` is returned as
  `NA_real_` by design: the residual precision and mixing proportion are point
  updates rather than variational factors, so no single joint mean-field bound
  is defined. The quantities that are defined are returned under their own
  names (`objective`, `objective_type`, `ELBO_conditional`). Likewise
  `fit_stats()` returns typed `NA` for covariance indices that are not yet
  derived for `vbmimic` fits, rather than substituting a plausible but wrong
  number. This is documented in `?vbfa`, `?fit_stats` and the vignettes.

* **Reproducibility.** `vbfa()` and `vbmimic()` consume no random numbers and
  take no seed argument; repeated runs are bit-identical. A stored-reference
  regression test pins `vbfa(v0 = 0.001)` to the published fixed-spike
  estimator.
