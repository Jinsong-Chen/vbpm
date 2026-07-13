# vbpm

**Variational Bayes Partially Confirmatory Measurement Models**

`vbpm` provides regularized variational Bayes estimation of partially
confirmatory factor-analytic measurement models — a single framework spanning
the confirmatory-to-exploratory continuum. Item loadings on *specified*
factors are anchored by a design matrix `Q`; loadings on *unspecified* factors
are selected from the data by continuous spike-and-slab priors.

## Features

- **`vbfa()`** — the estimator.
  - Spike-and-slab regularization of unspecified loadings.
  - **Dynamic exploration**: pass `v0` a decreasing vector for a warm-started
    regularization path (Ročková & George, 2018); a scalar `v0` is a single
    fixed spike. `v0 = 0.001` reproduces the fixed-spike estimator exactly.
  - **Bifactor**: `orthogonal = TRUE` fixes the factor correlation to the
    identity (one general + orthogonal group factors).
  - **Local dependence**: supply a residual design matrix `Qe` to estimate a
    sparse residual precision by a graphical spike-and-slab prior (QUIC).
- **`vbfa_fit()`** — SEM-like fit statistics (RMSEA, SRMR, CFI, TLI, AIC, BIC,
  ELBO).
- **`pefa()`** / **`select_K_elbow()`** — sweep the number of factors and
  select it with a scale-free gain rule.
- **`sim_fa()`** — simulate FA/IRT data with cross-loadings, local dependence,
  minor factors, categorical items, and missingness.

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
idx <- vbfa_fit(fit, Y, Q)
```

## Methods

- Partially confirmatory factor analysis with the Bayesian Lasso: Chen, J.,
  Guo, Z., Zhang, L., & Pan, J. (2021). A partially confirmatory approach to
  scale development with the Bayesian Lasso. *Psychological Methods*, 26(2),
  210–235. <https://doi.org/10.1037/met0000293>
- Variational estimator: Jin, Y., & Chen, J. (2025). Regularized variational
  approximation for partially confirmatory factor analysis. *Structural
  Equation Modeling*, 32(3), 437–449.
  <https://doi.org/10.1080/10705511.2024.2432612>
- Fit assessment and factor-number selection: Chen, J., & Jin, Y. (2026).
  Recovering latent structures after variational Bayesian variable selection.
  *arXiv:2607.07159*.

See `NEWS.md` for the current feature set and roadmap.

## License

GPL-3. Authors: Yi Jin and Jinsong Chen.
