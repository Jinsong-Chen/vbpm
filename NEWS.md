# vbpm 0.6.1

* **New `special_effects()`** post-processes an orthogonal bifactor fit
  into the three summaries of the `bifactor` vignette in one call: the
  Schmid-Leiman proportionality check (per-group and item-weighted mean CV
  of the group/general loading ratio), the approximate higher-order
  parameterization (`gamma` per group factor, `lambda` per item), and the
  testlet/special effect sizes of Zhang & Chen (2024, Eq. 16) with
  posterior-inclusion gating of unspecified loadings. Group membership is
  inferred from anchors and inclusion probabilities; unassigned,
  cross-loaded, and no-general-loading items are flagged rather than
  silently pooled. Deliberately descriptive: no cutoff argument and no
  verdict -- the conventional `.1` reading of the mean CV is mentioned as
  suggestive, and all quantities are reported regardless. Accepts
  `bifactor = TRUE` fits natively and legacy hand-built designs (one
  all-specified column under `orthogonal = TRUE`). Tests pin the function
  to an independent reimplementation of the vignette's plain-R code.

# vbpm 0.6.0

* **New `bifactor` mode in `vbfa()` and `pefa()`.** Previously a bifactor
  was expressed by hand-building a general column into `Q` under
  `orthogonal = TRUE`, which made the user-facing factor count include the
  general factor -- off by one against the comparable oblique model, a
  recurring source of confusion in sweeps ("`Kmin = 3` means 2 group
  factors"). With `bifactor = TRUE`, `Q` (and `pefa()`'s `Q0`) supply only
  the **group** design; the general column is added internally, and `K`
  everywhere means group factors, directly comparable to an oblique
  `K`-factor model. Details:
  - `bifactor = TRUE` implies orthogonal factors and **overrides
    `orthogonal = FALSE`** (with a message when the latter was passed
    explicitly).
  - A new `general` argument designs the general column (`1`/`0`/`-1` per
    item, default all `1`), so bifactor-(S-1)-style designs (some items
    with no general loading) and partially exploratory general columns are
    now expressible.
  - Passing a `Q` that already contains an all-specified column together
    with `bifactor = TRUE` is an error rather than a silent double count.
  - Fits carry `bifactor`, `n_general`, and `general`; `print()` shows
    "1 general + K group factors". The pefa sweep table gains a `K_total`
    column and the print methods label the window as group factors.
  - The long form (`cbind(1L, Q)` with `orthogonal = TRUE`) still works
    and is bit-identical -- the new mode is bookkeeping, not a new
    estimator -- and a regression test enforces that equality.
  - The relabeling changes none of the identification facts: bifactor
    sweeps remain weakly identified and the `?pefa` warnings keep full
    force.
  - Non-bifactor `orthogonal = TRUE` fits now print as "orthogonal"
    rather than "orthogonal/bifactor".

# vbpm 0.5.0

Documentation-and-parallelism release preparing the CRAN submission, plus one
estimator bug fix.

* **Removed: `vb_fit()`.** The pre-0.4.0 name for `fit_stats()` survived
  0.4.x as a thin wrapper; it is gone as of this release (before any CRAN
  submission, so no released API is broken). Use `fit_stats()`.
* `sim_fa()` gains `gamma` for **direct higher-order (testlet) generation**:
  a scalar or per-factor vector of second-order loadings in (0, 1) makes each
  group factor load on one general factor, and the data are generated from
  the equivalent Schmid-Leiman orthogonal bifactor matrix (returned as `MLA`,
  with the first-order matrix and `gamma` as `MLA1`/`gamma`). This replaces
  the hand-built-matrix recipe previously shown in the bifactor vignette.
* `plot.pefa()` redesigned. The objective panel gives the objective and BIC
  their own axes (one shared axis flattened both curves into unreadable
  lines); the fit panel splits absolute misfit (RMSEA/SRMR) and incremental
  fit (CFI/TLI) into side-by-side panels with conventional cutoffs; all
  panels use real plotting symbols that match their legends (previously
  `matplot()` drew digit characters while the legends showed symbols) and
  mark the selected `K` with a dotted vertical line.
* `print.summary.pefa()` now prints one compact comparison table with the raw
  objective and BIC, their marginal **gains as % of the largest gain** (the
  quantity the selection rule actually thresholds), the fit indices, and a
  marker on the selected row. Raw gains and timing remain in `$comparison`.
* **Orthogonal (bifactor) sweeps are now flagged as weakly identified.**
  Simulations from true bifactor and true higher-order populations show the
  ordinary bifactor sweep (general column + anchored groups, rest `-1`)
  systematically selects the lower window boundary regardless of the true
  `K`: an omitted group factor's rank-one covariance block is absorbed by
  the unspecified entries of the remaining columns, making the candidates
  covariance-equivalent -- and the truth-sized candidate inside such a sweep
  is itself not interpretable. `?pefa` gains a dedicated section; the
  `print()` methods replace the generic "extend the window" advice with an
  explicit caution when an orthogonal sweep lands on the lower boundary.
* **The bifactor vignette is rebuilt as a small tutorial paper.** It opens
  with the covariance geometry that drives everything (cross-cluster blocks
  are rank-one under both bifactor and oblique models, so a general factor
  is identified only by within-cluster rank-two structure), demonstrates
  the sweep failure on both true-bifactor and higher-order data, and closes
  with a validated **two-step decision route**: (1) an ordinary oblique
  `pefa()` sweep for the first-order factor count -- interior selections
  are trustworthy there even when the data are bifactor or higher-order;
  (2) anchored oblique-K versus anchored bifactor-(K+1), compared by BIC
  (not ELBO, whose prior mass differs across the two designs); (3) within
  the winner, the Schmid-Leiman mean-CV check (bifactor side) or a triad
  factorization of the fitted factor correlations (oblique side) settles
  bifactor versus higher-order versus plain correlated factors. Every
  population in the vignette routes correctly, and the route's power limits
  (a near-homogeneous true bifactor needs a large sample at step 2) are
  reported from a 3 x 3 x 3 simulation grid. `?pefa` summarizes the route.
* **Fixed: inflated residual variances for anchored items under local
  dependence.** The posterior variance of fixed-zero (`Q == 0`) loadings was
  left at its initialization value 1 instead of 0, and the LD branch's
  expected residual covariance summed loading variances over all columns --
  adding roughly +1 per fixed-zero column to the residual variance of every
  item that had any. On the standard test design, anchored items' residual
  variances came out ~5x too large (loadings and residual-edge selection were
  unaffected). Diagonal (`ld = FALSE`) results are unchanged -- the r1.1
  bit-exactness fixture still passes -- and the reported `Lam_var` is now 0
  for fixed-zero entries in both modes. The defect is inherited from the
  reference implementation `ld-vb_r.R`, which computes the same term the same
  way.
* `sim_fa()` gains `mla`: supplying a population loading matrix switches it to
  matrix-driven generation (pattern arguments are then ignored), with the same
  convention as `sim_lvm()`. Its factor-1-2 correlation argument is renamed
  `ph1` -> `ph12` to match `sim_lvm()`, so `ph1` now exists only in
  `sim_lvm()`, where it is the correlation among latent predictors. Both help
  pages gained matrix-driven examples and now read as a matched pair.
* New vignette **bifactor**: orthogonal bifactor estimation on true bifactor
  data, a `pefa()` sweep over bifactor candidates -- including an explicit
  note that `Kmin`/`Kmax` count the general factor, and a demonstration of
  why factor-number discovery is intrinsically hard for orthogonal bifactor
  structure (an omitted group factor is absorbed by the unspecified entries
  of the remaining columns) -- and the higher-order/testlet connection:
  simulating a higher-order CFA directly via `sim_fa(gamma = )`, checking the
  Schmid-Leiman constraint with a mean-CV rule (approximately higher-order
  when mean CV < .1), transforming back to higher-order parameters via the
  per-cluster mean ratio, and computing special-effect sizes by Eq. 16 of
  Zhang & Chen (2024) with posterior-inclusion-based selection of unspecified
  loadings. The testlet and higher-order models are noted as mathematically
  equivalent.
* `vbfa` vignette reorganized: restricted-`Qe` local dependence (testlet-block
  design matrices), a diagonal-vs-LD model comparison via `fit_stats()`,
  executed `v0`-path comparisons, an `ld_control` demonstration, and a
  "Missing data" section stating the missing-at-random assumption
  (Chen, 2021) with simulated and empirical demonstrations. The factor-number
  and bifactor material moved to the `pefa` and `bifactor` vignettes.
* `pefa` vignette expanded: the lead sweep now runs `vbfa()`'s default
  warm-started `v0` path on a true K = 4 simulation (N = 500, window 2:6 on a
  K0 = 2 backbone) and recovers it in the interior; checkpoint/resume with
  the provenance manifest (including the refusal on changed settings),
  executed `plot()` panels, a boundary demonstration on the same data, and a
  sweep on data with minor factors showing the gain rule ignoring dust
  dimensions.
* `vbmimic` vignette: missing responses in `Y` demonstrated.
* Housekeeping: stale `tmp/` artifacts and the sibling `vbpm.Rcheck/` removed;
  references across README and vignettes made consistent, adding Chen (2021,
  Multivariate Behavioral Research) and Zhang & Chen (2024, Applied
  Psychological Measurement).

# vbpm 0.4.0

* `vbmimic()` now also accepts missing responses in `Y`. Its residual
  covariance is diagonal, so the joint Gaussian conditional collapses to the
  factor-model mean `eta_i A'`, with the conditional variance `1 / V_j` carried
  into the residual sum of squares; the covariate information is already
  carried by `eta_i`. Missing values in `X` remain an error: covariates are
  conditioned on rather than modelled, so imputing them would require a
  distributional assumption the MIMIC model does not make.
* `vbfa()` now handles incomplete continuous-response data by deterministic,
  joint Gaussian conditional-moment updates. The model and conditional moments
  match the data-augmentation formulation used by LAWBL and the related PCFA
  MCMC work; rows or items with no observed responses are rejected.
* Local-dependence fits now expose the terminal VECM objective derived in the
  LD-PCFA-VB V4 manuscript. It is stored as `objective` with
  `objective_type = "vecm"`; `ELBO` remains `NA` because the residual precision
  and mixing probability are point updates rather than variational factors.
* New `fit_stats()` S3 generic works for all `vbpm_fit` objects. It computes
  covariance fit statistics for diagonal and LD `vbfa` fits and returns a
  clearly marked limited result for `vbmimic`; `vb_fit()` remains as a
  compatibility wrapper.
* `pefa()` now returns a `c("pefa", "vbpm_sweep")` object with candidate fits,
  a sweep table, `summary()`, `plot()`, and `selected_fit()`. Summaries report
  objective, objective-gain, ELBO, ELBO-gain, BIC, and BIC-gain selections.
* The package now has focused `vbfa`, `vbmimic`, and `pefa` vignettes.
* Verification for the above: the Gaussian conditional-moment kernel is now an
  internal function (`.cond_moments()`) tested directly against textbook
  multivariate-normal conditional means and covariances, for diagonal and full
  residual covariance; missing-data fits are checked for mask integrity and
  row-permutation invariance; LD fit statistics are checked to be built from
  the residual covariance `W` rather than the precision; and missing data is
  checked to compose with local dependence.

# vbpm 0.3.0

Fits are now classed objects, output names are plain-language, the default is
quiet, and the package has vignettes. One release, one breaking-change window,
before any public release — no public API ever carried the old behaviour.

* **S3 class `vbpm_fit`.** `vbfa()` returns `c("vbfa", "vbpm_fit")` and
  `vbmimic()` returns `c("vbmimic", "vbpm_fit")`. A fit is still a plain named
  list — `fit$Lam` etc. are unchanged and remain public API — but `print(fit)`
  now shows a compact summary instead of dumping every matrix, and `coef()`
  works. See `?vbpm_fit` for exactly what the class does and does not change.
* **Fits carry their own settings and designs.** `Q` (and `Qe`, `Q_A`, `Q_B`),
  `orthogonal`, and `ld` are stored in the fit, which is how `vb_fit()` and
  `print()` know what they are looking at without being told.
* **Output renames** (the variational `.q.` notation stays *inside* the
  estimators, where it mirrors the papers' derivations; the returned names are
  now plain language):
  * `sigsq.q.Lam` -> `Lam_var`, `PHI.q.eta` -> `eta_cov`,
    `M.q.PHIiver` -> `Phi_inv_mean` (in `vbfa()`);
  * `sigsq.q.A`/`sigsq.q.B` -> `A_var`/`B_var` (in `vbmimic()`);
  * `plotDat` dropped (always-`NULL` placeholder).
* **Quiet by default.** `vbfa(convChk = FALSE)` is now the default, and
  progress goes through `message()` (suppressible) instead of `cat()`.
  `vbmimic()` already behaved this way.
* **Vignettes**: "Getting started with vbpm" (Q design, fitting, fit
  statistics, factor-number selection, bifactor, local dependence, the NLSY
  empirical example, and turning graded prior knowledge into a Q matrix) and
  "MIMIC models with vbmimic".
* **Stored-reference regression test**: `vbfa(v0 = 0.001)` is now proven
  bit-identical to the frozen fixed-spike estimator behind the 2025 SEM paper
  (fixture generated from the original source; 79 iterations, all deltas
  exactly zero). This is the guarantee that lets analysis code migrate from
  the historical scripts to the package.

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

## Not implemented in this first release

All three items below were addressed in 0.4.0; they are recorded here as the
state of 0.1.0. See "Known limitations" in `README.md` for current status.

* Missing-data support in `vbfa()` (errored on `NA`). *Added in 0.4.0 for
  continuous responses, by in-loop Gaussian conditional moments.*
* ELBO under the local-dependence model (`NA`). *Resolved in 0.4.0 by
  supplying the quantities that are actually defined — the terminal VECM
  `objective` and the conditional bound `ELBO_conditional`. `ELBO` remains
  `NA` by design, because the precision and mixing proportion are point
  updates rather than variational factors, so no single joint mean-field
  bound exists to report.*
* LD-aware `vb_fit()`. *Added in 0.4.0 via `fit_stats()`, which builds the
  implied covariance from the residual covariance `W` and counts selected
  residual edges.*
