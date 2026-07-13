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
