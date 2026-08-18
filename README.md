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
    regularization path (Ročková & George, 2018; Jin, Chen, Yan, & Zhang,
    2026); a scalar `v0` is a single fixed spike. `v0 = 0.001` reproduces the
    fixed-spike estimator exactly.
  - **Bifactor**: `bifactor = TRUE` fits one general factor plus `K`
    orthogonal group factors, with `Q` supplying only the group design --
    so `K` counts group factors, comparable to an oblique `K`-factor
    model. The `general` argument designs the general column (including
    bifactor-(S-1)-style zeros). `orthogonal = TRUE` remains available
    for plain orthogonal models (or hand-built bifactor designs).
  - **Local dependence**: set `ld = TRUE` to estimate a sparse residual
    precision by a graphical spike-and-slab prior (QUIC) — fully exploratory
    by default, or restricted by a residual design matrix `Qe`
    (Jin, Chen, Yan, & Zhang, 2026).
- **`vbmimic()`** — the MIMIC estimator.
  - Spike-and-slab regularization of **both** the measurement design (`Q_A`,
    items on factors) and the structural design (`Q_B`, factors on covariates).
  - Supports the same warm-started `v0` path. The defaults reproduce the
    published estimator exactly.
  - The ELBO is not yet available for this model, so `fit_stats()` returns a
    clearly marked limited result and `pefa()` sweeps `vbfa()` candidates only.
- **`fit_stats()`** — SEM-like fit statistics (RMSEA, SRMR, CFI, TLI, AIC, BIC)
  plus the model's objective. An S3 generic over fitted models: it builds the
  implied covariance from the residual covariance, so it is correct for
  diagonal *and* local-dependence `vbfa` fits, and it reads `orthogonal`/`ld`
  from the fit rather than asking you to repeat them. Information criteria use
  the nominal hard-selected parameter count by default; the numerical
  Jacobian-rank count is an explicit `rank_adjust = TRUE` sensitivity.
- **`special_effects()`** — post-process a bifactor fit: Schmid-Leiman
  proportionality CVs, approximate higher-order parameters (second-order
  loadings per group factor), and testlet/special effect sizes
  (Zhang & Chen, 2024, Eq. 16).
- **`pefa()`** — fit one partially exploratory model for every `K` in a fixed
  consecutive window and report what was measured. It selects no factor count
  and carries no cut, threshold, persistence horizon, collision policy, or
  stopping rule; composing the numbers into a conclusion is yours. Three
  tables:
  - `$sweep` — one row per candidate: `K`, `ELBO`, `AIC`, `BIC`,
    `RMSEA`, `SRMR`, `CFI`, `TLI`, `t`, `iter`, `converged`.
  - `$transitions` — one row per *adjacent* pair: percent-of-largest-positive
    ELBO and BIC gains, plus the seven structural facts `phi_min`, `rmsd`,
    `rmsd_max`, `ari`, `pip_rmsd`, `unmatched_ssl`, `collision`.
  - `$persistence` — three upper-triangular matrices, `phi` (holding
    `phi_min`), `rmsd` (holding `rmsd_max`), and a logical `collision` mask,
    covering *every* ordered pair. Read them directly as
    `x$persistence$phi`; there is no accessor, because selecting a stored
    triangle would only be a synonym for `$`.

  Columns correspond by **minimum sign-aligned squared distance**: each source
  column takes the target minimizing the squared length of `a - s b`, with
  `s = 1` when `a'b >= 0` and `-1` otherwise — Equation 17 of Chen (2023) made
  invariant to whole-column reflection. Backbone columns pair by position, each
  exploratory column then chooses independently, and congruence only breaks
  exact ties. The assignment is deliberately not one-to-one: when two columns
  choose the same target that reuse is a **collision**, reported as evidence
  and never repaired into a forced matching. Every ordered pair is compared
  directly from its two endpoints, never chained through the candidates
  between them. Candidate matrices live in K-named `$loadings` and `$pips`
  lists (`result$loadings[["4"]]`). Convergence is *recorded*, never used to
  suppress a number that is defined: a candidate that stopped at `max_it`
  keeps its matrices and all of its comparisons.
- **`ssl()`** — per-candidate sums of squared loadings, `colSums(Lam^2)`, on
  the same scale as `unmatched_ssl`. It screens nothing and enters no
  correspondence rule; it is what distinguishes a low `phi_min` caused by a
  near-empty column from one caused by structure actually moving. Compact
  `print()`, `summary()`, and `plot()` methods complete the set; none of them
  announces a selected count, because none is chosen.
- **`sim_fa()`** and **`sim_lvm()`** — two complementary data generators.
  `sim_fa()` is driven by the loading *pattern* (items per factor, alternating
  cross-loadings, minor factors) and can generate **higher-order/testlet**
  data directly (`gamma =` second-order loadings, via the Schmid-Leiman
  equivalence); `sim_lvm()` is driven by a loading *matrix* and covers
  predictor designs (observed and latent) and mixed response formats. Both
  simulate local dependence, categorical items, and missingness.

## Installation

### 1. A C++ toolchain (one-time setup)

`vbpm` ships C++ code — the QUIC solver behind the local-dependence branch —
so installing it compiles that code and a C++ toolchain must be present.
Many R users already have one (any package installed *from source* needs it);
check before installing anything:

```r
# install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)   # TRUE = you are ready, skip to step 2
```

If that returns `FALSE`:

| system | what to install | where |
|---|---|---|
| **Windows** | **Rtools**, matched to your R version — Rtools45 for R 4.5.x, Rtools44 for R 4.4.x (check with `getRversion()`). Run the installer with its defaults; ≈500 MB installed. | <https://cran.r-project.org/bin/windows/Rtools/> |
| **macOS** | Xcode **command line tools** (not the full Xcode app): run `xcode-select --install` in Terminal and accept the prompt; ≈1 GB. | <https://mac.r-project.org/tools/> |
| **Linux** | Your distribution's R development packages: `sudo apt install r-base-dev` (Debian/Ubuntu) or `sudo dnf install R-devel gcc-c++` (Fedora/RHEL); a few tens of MB. | — |

Restart R afterwards, then re-run the `has_build_tools()` check.

### 2. The package

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
vignette("pefa")         # the factor-count sweep and its three tables
vignette("bifactor")     # bifactor structures, AO/AZ backbones, and sweeps
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
idx <- fit_stats(fit)                # nominal hard-selected count by default
```

Fits are ordinary named lists with a shared `vbpm_fit` class attached — fields
such as `fit$Lam`, `fit$pi`, and `fit$Phi` remain public API. Focused vignettes
cover `vbfa`, one `pefa()` sweep walked table by table under overextraction,
bifactor/higher-order/testlet structures with AO and AZ backbones, and
`vbmimic`.

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
- Warm-started regularization path and sparse residual/local-dependence
  estimation: Jin, Y., Chen, J., Yan, Z., & Zhang, Y. (2026). Sparse residual
  estimation in partially confirmatory factor analysis. *PsyArXiv preprint*.
  <https://doi.org/10.31234/osf.io/dehtv_v2>
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

As of 0.9.0:

- **Response types.** The estimators model **continuous Gaussian responses**,
  and in-loop missing-data support covers that case in both `vbfa()` and
  `vbmimic()`. Categorical and mixed responses are **out of scope for the
  current line of releases** — supporting them means adding a threshold model
  and augmenting latent responses, which is a modelling extension rather than
  a feature. Simulating them is supported (`cati`/`noc` in `sim_fa()`,
  `ilvl` in `sim_lvm()`); estimating from them is not.
- **Bifactor `pefa()` sweeps are a research extension.** Their `K` counts group
  factors, and every comparison removes the labelled general column before
  applying exactly the loading/PIP correspondence an ordinary oblique sweep
  uses. General and group columns can nevertheless redistribute common variance
  as `K` changes, so a study that needs both modes should run them on the same
  data, backbone, and window rather than assume that either route must agree.
  The `bifactor` vignette shows one such pair of sweeps, where the bifactor
  criterion path and persistence triangle separate the candidates far less than
  the ordinary ones do. That is a property of that example, not a general
  ranking of the two modes.
- **No calibrated thresholds.** Nothing in the package tells you how large a
  congruence, how small an RMSD, or how steep a gain has to be. The vignettes
  show readings composed in local code and label them as such; none of those
  cuts has been validated against an error rate here.
- **No DIF paths in `vbmimic()`.** Covariates predict the *factors*
  (`Q_B`, i.e. impact); there are no direct covariate-to-item paths, so the
  model cannot currently express or detect differential item functioning.
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
- Yi Jin — [\@jinyi1998](https://github.com/jinyi1998)
