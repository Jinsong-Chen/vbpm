# vbpm 0.8.3

This release makes the compact PEFA object easier to inspect and removes two
uninformative PIP summaries. It is a breaking refinement of the 0.8.2 schema;
the 0.8.2 release history below is unchanged.

* `$sweep` is now an ordinary scalar-only data frame. Candidate loading and
  PIP matrices move to top-level, K-named `$loadings` and `$pips` lists, so a
  candidate is read directly as `x$loadings[["4"]]` or `x$pips[["4"]]`.
  Their names always equal `as.character(x$sweep$K)`; full candidate fits are
  still not retained.
* `pip_sum` and `pip_product` are removed from `$transitions`. Both thresholded
  PIPs at `.5` and were diluted by fixed design cells; the product also returned
  zero whenever all changes ran in one direction. The replacement `pip_rmsd`
  is continuous: after loading-column matching, it is the RMSD between PIPs
  over matched cells that are selectable at both endpoints. Fixed backbone
  cells and unmatched columns are excluded. Loading-based metrics remain
  available when only a PIP endpoint is malformed.
* Loading alignment now follows Equation 17 directly: each retained
  non-backbone column in the smaller fit independently selects its nearest
  retained column in the larger fit. Reuse is allowed and the logical
  `$transitions$collision` flag reports when a larger-fit column is selected
  more than once; it is not silently resolved by a one-to-one assignment.
* Bifactor transitions no longer suppress all stability fields. The labelled
  general column is removed and the same alignment, screening, loading, PIP,
  collision, and unmatched-column definitions are applied to the group block,
  making same-data oblique and bifactor sweeps directly comparable. The full
  `K + 1` matrices remain available for candidate-level inspection.
  This supersedes the earlier blanket implementation warning: relative sweep
  reliability is now assessed empirically rather than decided from model class.
* The remaining `unmatched_max` stays a separate continuous,
  descriptive size diagnostic. It is not a loading-salience cutoff or an
  automatic selection rule.
* The bifactor vignette now runs ordinary-oblique and bifactor LPC sweeps on
  the same heterogeneous bifactor population and `K = 2:5` window. In this
  reproducible example the ordinary sweep selects the true three groups at
  both cuts, whereas the bifactor sweep selects its lower boundary. A derived
  SS-loading table (`colSums(Lam^2)`) exposes weak columns without adding a
  result field or selection rule; replicated R6 work will test whether this
  example generalizes.
* `summary()` derives a per-candidate `ssl` component -- the sums of squared
  loadings `colSums(Lam^2)` -- from `$loadings` at summary time, and its print
  method displays one row per candidate. In bifactor mode the values cover
  all `K + 1` columns with the general factor first, unlike the group-block
  transition diagnostics. SS loadings are descriptive column sizes, not Gram
  eigenvalues, and enter no selection rule; the sweep object itself stores no
  derived state.
* The top-level object is reduced to `$sweep`, `$transitions`, `$loadings`,
  `$pips`, `$selected_K`, `$boundary`, `$Q0`, and `$settings`. The lean
  `$settings` list contains only `cuts`, `sustain`, `bifactor`, `general`,
  `tau`, `stability_eps`, and `rank_adjust`. The requested window is derived
  from `$sweep$K`, and the backbone size from `ncol(x$Q0)`.

# vbpm 0.8.2

This release deliberately simplifies PEFA before 1.0 and changes both its
object schema and its adjacent-solution diagnostics.

* `fit_stats.vbfa()` now uses the deterministic nominal hard-selected
  parameter count by default (`rank_adjust = FALSE`). The numerical
  Jacobian-rank count remains available only by explicit request and requires
  the suggested package `numDeriv`; `rank_max_J` is an explicit computational
  guard. `pefa()` exposes and records the same controls, computing the chosen
  statistics while each candidate fit still exists.
* A `pefa` result now has exactly two data frames. `$sweep` retains scalar
  candidate statistics and adds `loading` and `pip` matrix list-columns;
  `$transitions` retains adjacent ELBO/BIC gains and the revised diagnostics.
  Named `$selected_K` and parallel `$boundary` vectors report the ELBO-gain
  result for every named cut. Raw ELBO/BIC and BIC-gain selections are no
  longer returned automatically, but remain directly reproducible from the
  retained scalar paths.
* Full candidate fits, `$selection`, and `selected_fit()` are removed. Extract
  loading/PIP matrices by matching a desired `K` to `$sweep$K`; explicitly
  refit that design when downstream work needs a complete `vbfa` object.
* Adjacent oblique candidates are compared with a deterministic,
  backbone-aware loading/PIP matcher. It holds the `Q0` columns fixed by
  position, screens free columns by squared L2 norm, and matches the retained
  free columns by minimum SSE. The transition fields are now `rmsd`,
  `rmsd_max`, signed `phi_min`, `ari`, `pip_sum`, `pip_product`, and the
  separate size diagnostic `max_unmatched_loading`; the previous `phi_mean`
  and ambiguous `max_loading` fields are removed.
* All transition diagnostics remain descriptive and do not enter automatic
  selection. Signed congruence can be unstable for near-zero-mean contrast
  columns; PIP summaries can be diluted by fixed design cells; and
  `max_unmatched_loading` is not a `.3` salient-loading rule. Bifactor edges,
  failed endpoints, and inapplicable or incomplete comparisons return typed
  `NA` diagnostics.
* General sweep mechanics and stability interpretation now live in the PEFA
  vignette. The bifactor vignette stays focused on the compact two-step
  first-order-count and matched-model comparison.
* Version 0.8.2 still has no persistence, checkpoint, resume, extension, or
  interrupt subsystem. Its internal sweep seam is groundwork only: this
  release does not export `extend.pefa()`, `vblvm()`, MIMIC sweeping, ESEM
  support, or a generic sweep engine.

# vbpm 0.8.1

`pefa()` is reduced to a compact, in-memory factor-number sweep. This is a
breaking cleanup of the 0.8.0 interface; the 0.8.0 release history below is
left unchanged.

* The paper-defined workflow now fits only ordinary oblique candidates with
  diagonal residuals. The research-extension `bifactor = TRUE` mode remains,
  with `K` consistently denoting the number of group factors. Standalone
  orthogonal and local-dependence models remain available through `vbfa()` but
  are no longer PEFA sweep modes.
* ELBO and BIC each receive raw and gain selections. AIC remains in `$sweep`
  as a descriptive fit statistic but is not automated as a gain rule.
  `select_K_elbow()` and `pefa()` now use the same full sustained-drop scan: a
  stop requires all `sustain` following gains, including near the right edge.
* RMSEA, SRMR, CFI, and TLI are descriptive trajectories, not eligibility
  gates. A selection requires one complete, converged, finite criterion path;
  failed paths are not bridged or split into later components.
* The result has three compact tables: candidate-level `$sweep`, edge-level
  `$transitions`, and rule-level `$selection`. Raw adjacent ELBO and BIC gains
  are stored with their transitions, while each rule records its own gain
  normalizer and threshold.
* Every candidate fit is retained in `$fits`, and `selected_fit()` retrieves
  the primary selected fit. Removed are `keep_fits`, duplicate selected-fit
  fields, compatibility aliases and schema fields, the secondary
  `vbpm_sweep` class, and `transition_detail()`.
* Removed without deprecated aliases are the 0.8.0 `delta`, `fit_cut`,
  `cutoffs`, `orthogonal`, and `keep_fits` arguments. Named `cuts` remains the
  only multi-rule interface; PEFA rejects local-dependence and orthogonal mode
  controls passed through `...`.
* Adjacent-count stability is deliberately small and descriptive. Consecutive
  oblique solutions are aligned by a simple best-unused-column convention and
  report `phi_min`, `phi_mean`, `ari`, `rmsd`, and the unmatched column's
  `max_loading`. There is no assignment search, tie policy, fixed stability
  cutoff, status vocabulary, or public matching detail.
* `pefa()` has no persistence, checkpoint, resume, interrupt, or overlap-reuse
  subsystem. A boundary selection asks for a complete rerun with a better
  bracketed window.
* PEFA fit statistics now use the nominal hard-selected parameter count
  consistently. Hard-selected CFI returns one when its comparison denominator
  is zero, matching the paper's stated convention. **This changes reported
  numbers, not just interfaces**: on identical data, `$sweep` AIC/BIC (and
  potentially a BIC selection) can differ from a 0.8.0 run that silently used
  the optional Jacobian rank adjustment, and the new column alignment changes
  every stability column relative to 0.8.0's matcher.
* The PEFA vignette now follows the paper's compact workflow and labels
  adjacent-count stability as a package extension. The bifactor vignette is
  reduced to definitions, direct anchored fitting, `special_effects()`, and a
  two-step first-order-count then matched BIC comparison.

# vbpm 0.8.0

`pefa()` is reorganized around what each quantity is indexed by, and gains
adjacent-count stability diagnostics. **This release changes results as well
as structure** -- see the two breaking notes below.

* **Three tables instead of one.** Everything a sweep produces is indexed by
  a candidate, a transition between adjacent candidates, or a selection rule:
  - `$sweep` keeps one row per candidate (ELBO, Objective, AIC, BIC, fit
    indices, convergence, timing) and **loses all six gain columns**;
  - **`$transitions`** is new -- one row per adjacent pair, holding the gains
    (which are properties of a *step*, not of a candidate) together with the
    stability indexes below;
  - **`$selection`** becomes a long-form table with one row per criterion x
    form x cut x eligibility scope, recording the normalizer `gain_max`, the
    absolute `threshold`, `selected_K`, and both boundary readouts.
  No percentage-of-largest-gain column is stored anywhere: the normalizer
  depends on the eligibility scope, so percentages are computed for display.
* **Adjacent-count stability.** Each transition reports how much of the
  `K`-factor solution survives inside the `K+1`-factor solution the sweep
  already fitted: `phi_min` (the binding index), `phi_mean`, `n_phi90`, `ari`,
  `rmsd`, `surplus_max`, per-endpoint assignment-coverage counts, and
  `stability_status`. Columns are matched by **bottleneck (max-min) Tucker
  congruence**, which -- unlike the greedy stored-order matching used in
  hand-rolled versions -- is permutation invariant and actually maximizes the
  binding statistic. `transition_detail()` returns the per-column material
  behind one row.
* **The package computes, it does not render verdicts.** No `stable` column,
  no comparison against any threshold. The manuscript's `.85` operating point
  is documented as provisional and specific to its own workflow.
* **Stability is gated where loadings are not interpretable.** Under
  `orthogonal`/`bifactor` sweeps the padded columns are rotationally
  indeterminate, and a non-converged endpoint is not a solution; in both cases
  the loading-based columns are `NA` with the reason in `stability_status`.
  Gains remain valid in every mode.
* **New `cuts` argument** replaces the overloaded scalar `delta`: a named
  vector of gain thresholds with exactly one `primary` entry, so a primary
  rule and any number of sensitivity rules are evaluated in one call. The
  default stays `c(primary = 10)`; `delta` continues to work as a deprecated
  alias, and supplying both is an error.
* **New `keep_fits`** (`"all"`, `"selected"`, `"none"`) for callers who want
  the tables without the fit payload. `selected_fit` is always returned.
  `summary.pefa()` now carries the three tables and **no** fit payload.
* **Breaking, structural:** code reading `ELBO_gain`, `BIC_gain`, or any
  `*_gain_pct` column from `$sweep` must read `$transitions` instead;
  `summary()$comparison` is now `summary()$sweep`; `$selection` is a data
  frame rather than a named vector. `$boundary` and `$delta` survive as
  deprecated aliases of `$window_boundary` and `$cuts[["primary"]]`.
* **Breaking, in results:** selection no longer bridges a gap. When an
  interior candidate is ineligible (non-converged, or failing the absolute-fit
  gate), the old code differenced *across* it and reported a two-factor jump
  as a one-factor marginal gain. Now an unusable edge is skipped and resets the
  sustained-run counter, so `selected_K` can differ from 0.7.1 on any sweep
  with an interior gap. Sweeps without a gap are unaffected.
* **Removed: the `save_path` checkpoint/resume subsystem.** `pefa()` performs
  one in-memory sweep and returns one complete object; it writes nothing and
  cannot resume an interrupted call. Saving a result is ordinary
  `saveRDS()` workflow. The manifest that made resuming trustworthy is gone
  with it and is not replaced, because a single in-memory sweep cannot pool
  candidates fit under different settings in the first place.
* Fixed: `print.pefa()` inferred the fitting mode from `$fits[[1]]`, so its
  orthogonal weak-identification caution silently vanished whenever the first
  candidate's fit was absent. Both print methods now read explicit top-level
  mode settings. Reported also in `?pefa`.
* `pefa` vignette: the checkpoint walkthrough is replaced by a section on
  saving a sweep and widening its window, which shows that re-running a wider
  window *reproduces* the shared candidates bit-for-bit rather than wasting
  them, and that `Kmin` cannot go below `K0`.
* `bifactor` vignette, Section 6: three corrections to the two-step process.
  Cut disagreement is no longer described as the workflow's protection
  against absorbed local dependence -- it is not a dependence detector, and
  the hard case leaves every count-level indicator clean; residual inspection
  is now an unconditional step, labelled as the exploratory sensitivity it is.
  The doublet repair no longer adds a two-item factor, which perturbs the two
  competitors *differently* (the oblique arm estimates the new column's
  correlations while `bifactor = TRUE` forces it orthogonal); it now frees the
  same residual edge in both arms with `ld = TRUE` and a shared `Qe`, which
  reverses the contaminated verdict at both anchoring depths. The thin `-9`
  LPC margin is read against the vignette's own small-margin warning, and the
  accompanying manuscript is now cited by name rather than as "the
  accompanying manuscript."

# vbpm 0.7.1

* Documentation currency pass. The README's "Known limitations" header was
  still dated 0.5.0, and the list was missing two limitations established
  since: **bifactor factor-number sweeps are weakly identified** (an omitted
  group factor is absorbed by the unspecified entries of the remaining
  columns, so the selection is uninformative -- with a pointer to the
  two-step route that does work), and **`vbmimic()` has no direct
  covariate-to-item paths**, so it models impact but cannot express or detect
  differential item functioning. A quick-start comment that named the
  `$design`/`$settings` fields removed in 0.7.0 was reworded.
* Dropped three unused namespace imports (`stats::rnorm`, `stats::var`,
  `stats::quantile`).

# vbpm 0.7.0

Slimming release: fitted objects now store each quantity exactly once. No
estimator changes -- every number a model produces is unchanged.

* **Removed the mirror lists `$coefficients`, `$design`, `$posterior`, and
  `$settings`** from `vbfa()` and `vbmimic()` fits. Each of their elements
  duplicated a top-level component (`Lam`, `Phi`, `pi`, `Lam_var`,
  `eta_cov`, `Q`, `Qe`, `orthogonal`, `ld`, `bifactor`, `v0`); they had no
  consumers anywhere in the package, its tests, or its vignettes, and
  `coef()` never read them. `vbmimic()`'s `standardize` setting, the only
  element that lived solely in a mirror, is now a top-level component.
* **`preprocess$missing_mask` is stored only when the data actually contain
  missing values** (`NULL` otherwise). For complete data it was an all-`FALSE`
  `N x J` matrix and the single largest part of a fit -- 65% of the object at
  `N = 3000`, `J = 27`. `preprocess$n_missing` still reports the count in all
  cases.
* **`pefa()` no longer returns the `$fit` alias**, which was a third copy of
  the selected model (already available as `$selected_fit` and inside
  `$fits`) with no consumers.
* Net effect: a saved `vbfa` fit is about **half** its previous size
  (144 KB -> 74 KB for `N = 3000`, `J = 27`), and a saved sweep about 20%
  smaller. A new test pins these invariants so the duplication cannot
  silently return.
* Internal: the QUIC wrapper no longer copies the solver's output array into
  an identical new array.
* **README: proper C++ toolchain instructions.** Because the package ships
  compiled code, installing it requires a toolchain -- the one real setup
  barrier. The installation section now opens with a one-line check
  (`pkgbuild::has_build_tools(debug = TRUE)`, a genuine compile test on every
  platform) and, if it fails, per-platform instructions and links: Rtools
  matched to the R version on Windows, `xcode-select --install` on macOS, and
  the distribution's R development packages on Linux.

Nothing was on CRAN when these fields were removed, so no released API is
affected; scripts that read the documented top-level components
(`fit$Lam`, `fit$pi`, `fit$Q`, ...) are unaffected.

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
