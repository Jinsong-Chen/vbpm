# vbpm 0.2.0

* `vbmimic()` — regularized variational Bayes for extended multiple-indicators
  multiple-causes (MIMIC) models, with spike-and-slab priors on both the
  measurement (`Q_A`) and the structural (`Q_B`) part. Supports a warm-started
  `v0` path. The defaults (`v0 = 0.001`, `standardize = FALSE`) reproduce the
  estimator of Jin & Chen (2025, *Multivariate Behavioral Research*)
  bit-for-bit. The ELBO is not yet available, so `vb_fit()` and `pefa()` do
  not accept `vbmimic` fits.
* `vbfa()` gains an explicit `ld` switch for local-dependence estimation.
  `ld = FALSE` (default) fits diagonal residuals and ignores `Qe` (with a
  warning if one is supplied); `ld = TRUE` with `Qe = NULL` means fully
  exploratory local dependence (`Qe = matrix(-1, J, J)`); an explicit `Qe` is
  validated as before. Previously the mode was toggled implicitly by whether
  `Qe` was supplied.
* Fixed: `vb_fit()` defaulted `orthogonal = FALSE` regardless of how the model
  was fit, so calling it on a bifactor (`orthogonal = TRUE`) fit without
  repeating the argument silently counted `K(K-1)/2` nonexistent factor
  correlations, biasing `t`, AIC, BIC, RMSEA, CFI, and TLI. `vbfa()` now
  records `orthogonal` (and `ld`) in the fit object and `vb_fit()` reads it
  from there; a contradictory explicit value is an error.
* `nlsy27` — the NLSY 1997 empirical illustration dataset (3,458 x 27,
  mixed-type), carried over from LAWBL.
* Examples added or extended for `vbfa()` (bifactor, local dependence,
  empirical data), `vb_fit()`, `pefa()`, and `select_K_elbow()`.
* `sim_lvm()` — latent-variable-model data generator, absorbed from the LAWBL
  package, which `vbpm` supersedes. Simulates observed predictors (`P`, `b`)
  and latent predictors (`K1`, `ph1`, `b1`) in addition to the factor-analytic
  features. It complements `sim_fa()`: `sim_lvm()` is driven by a loading
  matrix and covers predictor designs and mixed response formats, while
  `sim_fa()` is driven by the loading *pattern* (items per factor, alternating
  cross-loadings, minor factors).
  * New `phd` argument: correlated structural disturbances, so designs with an
    unrestricted factor residual covariance can be simulated.
  * New `phx` argument: the covariate correlation, previously tied to `phi`.
  * Fixed: with covariates present, the missingness mask was recycled across
    the covariate block, so `NA`s appeared in `X`. Missingness is now confined
    to the `J` items.
  * Added validation for the predictor coefficient matrices and for designs
    that imply a non-positive disturbance variance.

# vbpm 0.1.0

First release.

* `vbfa()` — regularized variational Bayes factor analysis across the
  confirmatory–exploratory continuum. Spike-and-slab priors on unspecified
  loadings; optional warm-started dynamic regularization path (`v0` a vector);
  orthogonal bifactor parameterization (`orthogonal = TRUE`); optional sparse
  residual / local-dependence estimation (`Qe`) via a graphical spike-and-slab
  prior solved by QUIC. A scalar `v0 = 0.001` reproduces the fixed-spike
  estimator exactly.
* `vb_fit()` — SEM-like fit statistics (RMSEA, SRMR, CFI, TLI, AIC, BIC,
  ELBO) with hard selection as the default.
* `pefa()` and `select_K_elbow()` — window sweep over the number of factors
  with a scale-free gain rule for factor-number selection.
* `sim_fa()` — simulate factor-analytic and IRT data with cross-loadings,
  local dependence, minor factors, categorical responses, and missingness.

## Not yet implemented (planned for 0.1.1 / 0.2)

* Missing-data support in `vbfa()` (currently errors on `NA`).
* ELBO under the local-dependence model (currently `NA`).
* LD-aware `vb_fit()`.
