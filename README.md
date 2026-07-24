# vbpm

**Variational Bayes Psychometric Models**

`vbpm` is a home for variational Bayes estimation of psychometric measurement
models. Its models are *partially confirmatory*: they span the
confirmatory-to-exploratory continuum, with parameters on *specified* paths
anchored by a design matrix and parameters on *unspecified* paths selected from
the data by continuous spike-and-slab priors.

The package provides **variational Bayes factor analysis** (`vbfa`), where a
design matrix `Q` governs the item loadings, and **regularized MIMIC**
(`vbmimic`), which adds covariates predicting the factors and regularizes both
the measurement part (`Q_A`) and the structural part (`Q_B`).

## Features

- **`vbfa()`** — the factor-analysis estimator.
  - Spike-and-slab regularization of unspecified loadings.
  - **Dynamic exploration**: pass `v0` a decreasing vector for a warm-started
    regularization path (Ročková & George, 2018); a scalar `v0` is a single
    fixed spike. `v0 = 0.001` reproduces the fixed-spike estimator exactly.
  - **Bifactor**: `bifactor = TRUE` fits one general factor plus `K`
    orthogonal group factors, with `Q` supplying only the group design --
    so `K` counts group factors, comparable to an oblique `K`-factor
    model. The `general` argument designs the general column (including
    bifactor-(S-1)-style zeros). `orthogonal = TRUE` remains available
    for plain orthogonal models (or hand-built bifactor designs).
  - **Local dependence**: set `ld = TRUE` to estimate a sparse residual
    precision by a graphical spike-and-slab prior (QUIC) — fully exploratory
    by default, or restricted by a residual design matrix `Qe`.
- **`vbmimic()`** — the MIMIC estimator.
  - Spike-and-slab regularization of **both** the measurement design (`Q_A`,
    items on factors) and the structural design (`Q_B`, factors on covariates).
  - Supports the same warm-started `v0` path. The defaults reproduce the
    published estimator exactly.
  - The ELBO is not yet available for this model, so `fit_stats()` returns a
    clearly marked limited result and `pefa()` does not accept `vbmimic` fits.
- **`fit_stats()`** — SEM-like fit statistics (RMSEA, SRMR, CFI, TLI, AIC, BIC)
  plus the model's objective. An S3 generic over fitted models: it builds the
  implied covariance from the residual covariance, so it is correct for
  diagonal *and* local-dependence `vbfa` fits, and it reads `orthogonal`/`ld`
  from the fit rather than asking you to repeat them.
- **`pefa()`** / **`select_K_elbow()`** — sweep the number of factors and
  select it with a scale-free gain rule. `pefa()` returns a sweep object with
  `summary()`, `plot()`, and `selected_fit()`.
- **`sim_fa()`** and **`sim_lvm()`** — two complementary data generators.
  `sim_fa()` is driven by the loading *pattern* (items per factor, alternating
  cross-loadings, minor factors) and can generate **higher-order/testlet**
  data directly (`gamma =` second-order loadings, via the Schmid-Leiman
  equivalence); `sim_lvm()` is driven by a loading *matrix* and covers
  predictor designs (observed and latent) and mixed response formats. Both
  simulate local dependence, categorical items, and missingness.

## Installation

A C++ toolchain is required to compile the QUIC solver (Rtools on Windows,
Xcode command line tools on macOS).

```r
# install.packages("remotes")
remotes::install_github("Jinsong-Chen/vbpm", build_vignettes = TRUE)
```

Or clone and install from source:

```sh
git clone https://github.com/Jinsong-Chen/vbpm.git
R CMD INSTALL vbpm
```

Then:

```r
library(vbpm)
vignette("vbfa")         # start here
vignette("pefa")         # choosing the number of factors
vignette("bifactor")     # bifactor, higher-order, and testlet structures
vignette("vbmimic")      # covariates predicting the factors
citation("vbpm")
```

Building the vignettes needs `knitr`, `rmarkdown`, and pandoc (bundled with
RStudio; from a plain shell, point `RSTUDIO_PANDOC` at it).

## Quick start

```r
library(vbpm)

sim <- sim_fa(N = 500, K = 3, ipf = 6, lam = .7, lac = .3)
Y   <- sim$dat

## partially confirmatory Q: 1 = anchor, 0 = fixed zero, -1 = estimate
Q <- matrix(-1, ncol(Y), 3)
Q[1:2, ] <- 0; Q[1:2, 1] <- 1        # two anchors on factor 1, etc.

fit <- vbfa(Y, Q)                    # dynamic path on by default; quiet
fit                                  # compact summary (S3 print method)
idx <- fit_stats(fit)                # reads data/design/settings from the fit
```

Fits are ordinary named lists with a shared `vbpm_fit` class attached — fields
such as `fit$Lam`, `fit$pi`, and `fit$Phi` remain public API. Focused vignettes
cover `vbfa`, factor-number selection (`pefa`), bifactor/higher-order/testlet
structures, and `vbmimic`.

MIMIC, with covariates predicting the factors:

```r
B <- matrix(0, 3, 9)                 # sparse structural design: 3 covariates
for (k in 1:3) B[k, (k*3-2):(k*3)] <- .3   # per factor

sim <- sim_lvm(N = 500, K = 3, J = 18, P = 9, b = B, phx = 0)
Y   <- sim$dat[, 1:18]               # items first ...
X   <- sim$dat[, 19:27]              # ... covariates last

## anchor two items per factor; leave the structural part exploratory.
## (With BOTH parts fully exploratory the solution can be rotationally
## ambiguous -- anchor at least one part.)
Q_A <- matrix(-1L, 18, 3)
for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
                 Q_A[a, ] <- 0L; Q_A[a, k] <- 1L }
Q_B <- matrix(-1L, 3, 9)

fit <- vbmimic(Y, X, Q_A, Q_B)
round(fit$B, 2)                      # which covariates predict which factors
```

## References

- Partially confirmatory factor analysis with the Bayesian Lasso: Chen, J.,
  Guo, Z., Zhang, L., & Pan, J. (2021). A partially confirmatory approach to
  scale development with the Bayesian Lasso. *Psychological Methods*, 26(2),
  210–235. <https://doi.org/10.1037/met0000293>
- Generalized PCFA (mixed response formats; missing data under MAR): Chen, J.
  (2021). A generalized partially confirmatory factor analysis framework with
  mixed Bayesian Lasso methods. *Multivariate Behavioral Research*, 57(6),
  879–894. <https://doi.org/10.1080/00273171.2021.1925520>
- Bifactor, higher-order, and testlet effects within GPCFA (including effect
  sizes for testlet effects): Zhang, Y., & Chen, J. (2024). Accommodating and
  extending various models for special effects within the generalized
  partially confirmatory factor analysis framework. *Applied Psychological
  Measurement*, 48(4–5), 208–229. <https://doi.org/10.1177/01466216241261704>
- Fully and partially exploratory factor analysis: Chen, J. (2023). Fully and
  partially exploratory factor analysis with bi-level Bayesian regularization.
  *Behavior Research Methods*, 55(4), 2125–2142.
  <https://doi.org/10.3758/s13428-022-01884-7>
- Variational estimator: Jin, Y., & Chen, J. (2025). Regularized variational
  approximation for partially confirmatory factor analysis. *Structural
  Equation Modeling*, 32(3), 437–449.
  <https://doi.org/10.1080/10705511.2024.2432612>
- MIMIC: Jin, Y., & Chen, J. (2025). Regularized variational Bayesian
  approximations for variable selection in extended multiple-indicators
  multiple-causes models. *Multivariate Behavioral Research*.
  <https://doi.org/10.1080/00273171.2025.2483253>
- Fit assessment and factor-number selection: Chen, J., & Jin, Y. (2026).
  Recovering latent structures after variational Bayesian variable selection:
  Fit assessment and factor-number selection in partially exploratory factor
  analysis. *arXiv preprint* arXiv:2607.07159.
  <https://arxiv.org/abs/2607.07159>

## Known limitations

As of 0.5.0:

- **Response types.** The estimators model **continuous Gaussian responses**,
  and in-loop missing-data support covers that case in both `vbfa()` and
  `vbmimic()`. Categorical and mixed responses are **out of scope for the
  current line of releases** — supporting them means adding a threshold model
  and augmenting latent responses, which is a modelling extension rather than
  a feature. Simulating them is supported (`cati`/`noc` in `sim_fa()`,
  `ilvl` in `sim_lvm()`); estimating from them is not.
- **Complete covariates.** `vbmimic()` requires complete `X`: covariates are
  conditioned on rather than modelled, so the model supplies no distribution
  to impute them from.
- **LD objective.** Local-dependence `vbfa()` fits return the
  terminal VECM objective and LD-aware fit statistics. `ELBO` is intentionally
  `NA`: calling the point-updated precision objective a joint variational lower
  bound would be misleading.
- **ELBO for MIMIC.** `vbmimic()` converges on a residual criterion rather than
  the bound, so `ELBO` is `NA`. `fit_stats()` describes the available parameter
  counts but covariance fit indices await a fully specified joint covariance.
- **Fully exploratory MIMIC.** With both `Q_A` and `Q_B` left entirely
  unspecified the model converges cleanly but the solution can be rotationally
  ambiguous. Anchor at least one part.

Both estimators reproduce their published references bit-for-bit, and the
`vbfa` guarantee is enforced by a stored-reference regression test that runs on
every `R CMD check`.

See `NEWS.md` for the per-version changelog.

## License

GPL-3.

## Authors

- Jinsong Chen (maintainer) — [\@Jinsong-Chen](https://github.com/Jinsong-Chen)
- Yi Jin — [\@jinyi1998](https://github.com/jinyi1998), first author of the
  method papers behind `vbfa()` and `vbmimic()`
