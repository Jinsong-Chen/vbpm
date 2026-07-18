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
  - **Bifactor**: `orthogonal = TRUE` fixes the factor correlation to the
    identity (one general + orthogonal group factors).
  - **Local dependence**: set `ld = TRUE` to estimate a sparse residual
    precision by a graphical spike-and-slab prior (QUIC) — fully exploratory
    by default, or restricted by a residual design matrix `Qe`.
- **`vbmimic()`** — the MIMIC estimator.
  - Spike-and-slab regularization of **both** the measurement design (`Q_A`,
    items on factors) and the structural design (`Q_B`, factors on covariates).
  - Supports the same warm-started `v0` path. The defaults reproduce the
    published estimator exactly.
  - The ELBO is not yet available for this model, so `vb_fit()` and `pefa()`
    do not accept `vbmimic` fits.
- **`vb_fit()`** — SEM-like fit statistics (RMSEA, SRMR, CFI, TLI, AIC, BIC,
  ELBO).
- **`pefa()`** / **`select_K_elbow()`** — sweep the number of factors and
  select it with a scale-free gain rule.
- **`sim_fa()`** and **`sim_lvm()`** — two complementary data generators.
  `sim_fa()` is driven by the loading *pattern* (items per factor, alternating
  cross-loadings, minor factors); `sim_lvm()` is driven by a loading *matrix*
  and covers predictor designs (observed and latent) and mixed response
  formats. Both simulate local dependence, categorical items, and missingness.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Jinsong-Chen/vbpm")
```

A C++ toolchain is required to compile the QUIC solver (Rtools on Windows).

## Quick start

```r
library(vbpm)

sim <- sim_fa(N = 500, K = 3, ipf = 6, lam = .7, lac = .3)
Y   <- sim$dat

## partially confirmatory Q: 1 = anchor, 0 = fixed zero, -1 = estimate
Q <- matrix(-1, ncol(Y), 3)
Q[1:2, ] <- 0; Q[1:2, 1] <- 1        # two anchors on factor 1, etc.

fit <- vbfa(Y, Q)                    # dynamic path on by default
idx <- vb_fit(fit, Y, Q)
```

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

## Methods

- Partially confirmatory factor analysis with the Bayesian Lasso: Chen, J.,
  Guo, Z., Zhang, L., & Pan, J. (2021). A partially confirmatory approach to
  scale development with the Bayesian Lasso. *Psychological Methods*, 26(2),
  210–235. <https://doi.org/10.1037/met0000293>
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

See `NEWS.md` for the current feature set and roadmap.

## License

GPL-3. Authors: Jinsong Chen (maintainer) and Yi Jin.
