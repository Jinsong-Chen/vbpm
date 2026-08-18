# cran-comments

## Submission summary

This is a new submission of **vbpm** 0.9.0 (Variational Bayes Psychometric
Models). No earlier version was released on CRAN.

The package provides variational Bayes estimation for partially confirmatory
psychometric measurement models:

* `vbfa()` — factor analysis with spike-and-slab selection of unspecified
  loadings, a warm-started regularization path, an orthogonal bifactor
  parameterization, and optional sparse residual (local-dependence) estimation;
* `vbmimic()` — the MIMIC extension, regularizing both the measurement and the
  structural design;
* `fit_stats()` and `special_effects()` — SEM-like fit statistics with nominal,
  hard-selected, and soft/effective parameter counts, and post-processing of a
  fitted bifactor model;
* `pefa()` — a fixed-window factor-count sweep that returns an **evidence
  object** and makes no decision. It selects no factor count, applies no
  threshold, and vetoes no comparison; it returns candidate fit evidence
  (`$sweep`), adjacent criterion gains (`$transitions`), and all-pairs
  structural-persistence evidence (`$persistence`), together with the candidate
  loading and PIP matrices and an immutable SHA-256 `$provenance` record;
* `persistence()`, `ssl()`, and `verify_pefa()` — accessors that pivot one
  structural metric into a horizon triangle, report per-candidate sums of
  squared loadings, and re-verify a stored evidence object against the
  installed package;
* `select_K_elbow()` — an optional criterion-generic count calculator that
  applies a user-supplied relative-gain rule to any larger-is-better path.
  `pefa()` never invokes it;
* `sim_fa()` and `sim_lvm()` — two data generators.

## Test environments

<!-- PLACEHOLDER: replace with the actual submission-environment run before
     tagging. Do not carry a candidate-run result into the final file. -->

* PLACEHOLDER — submission environment (OS, R version, check invocation)

## R CMD check results

<!-- PLACEHOLDER: fill from the exact-commit `R CMD check --as-cran
     --run-donttest` run against the retained tarball. Every count and NOTE
     text below must be copied from that log, not restated from memory. -->

Status: PLACEHOLDER — NOTEs / ERRORs / WARNINGs from the exact-commit run.

* PLACEHOLDER — verbatim text of each NOTE, with a one-line explanation of
  whether it is package-attributable or environmental.

Compiled installation, examples, the full test suite, Rd validation, all four
vignettes, and indexed PDF-manual generation must all pass on that run before
the tag is created.

## Notes for the reviewer

* **New submission.** The incoming-feasibility "New submission" NOTE is
  expected. There is no archived CRAN version of this package to reference.

* **Possibly mis-spelled words in DESCRIPTION.** `vbfa` and `vbmimic` are the
  package's estimator function names; `QUIC` is the sparse inverse-covariance
  algorithm used by the local-dependence branch; `MIMIC` is the standard
  multiple-indicators multiple-causes acronym. All method references in the
  Description carry DOIs, including the preprints
  <doi:10.31234/osf.io/dehtv_v2> for the warm-started regularization path and
  the local-dependence branch, and <doi:10.48550/arXiv.2607.07159> for the fit
  and factor-count evidence.

* **Compiled code.** `src/` contains a QUIC (quadratic approximation for sparse
  inverse covariance) solver used only by the optional local-dependence branch,
  built with Rcpp/RcppArmadillo. `src/Makevars` links LAPACK/BLAS explicitly.
  There is no other compiled code and no external system dependency.

* **The `digest` dependency.** `digest` is a hard Import used solely to compute
  the SHA-256 fields in `pefa()`'s `$provenance` record and to re-verify them in
  `verify_pefa()`. Serialization options are fixed in one internal helper
  (canonical XDR, serialization version 3, with the writer-dependent stream
  header excluded), so an unchanged scientific payload hashes identically across
  R versions and locales. `tools::sha256sum(bytes =)` would remove the
  dependency but was added in R 4.5.0, which would raise this package's floor
  from `R (>= 4.1)`; `digest` itself requires only `R (>= 3.3.0)` and imports
  `utils` alone.

* **Included data.** `nlsy27` (3,458 x 27 responses) is derived from the public
  National Longitudinal Survey of Youth 1997 (Bureau of Labor Statistics, U.S.
  Department of Labor) and was previously distributed by the same maintainer in
  the CRAN package **LAWBL** (GPL-3). Provenance is documented in `?nlsy27`.

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
  derived for `vbmimic` fits, and `pefa()` returns a typed `NA` for a
  structural metric whose mathematical inputs are undefined, rather than
  substituting a plausible but wrong number. Each case is documented in
  `?vbfa`, `?fit_stats`, `?pefa`, and the vignettes.

* **Reproducibility.** `vbfa()` and `vbmimic()` consume no random numbers and
  take no seed argument; repeated runs are bit-identical. A stored-reference
  regression test pins `vbfa(v0 = 0.001)` to the published fixed-spike
  estimator. The test suite also re-verifies stored `pefa()` evidence objects
  through `verify_pefa()`, so a change in matcher or hashing semantics fails
  the tests rather than passing silently.
