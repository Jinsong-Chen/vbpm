## Partially exploratory factor analysis: a window sweep over the number of
## factors with a scale-free gain rule for factor-number selection. See the
## roxygen documentation of pefa() and select_K_elbow() below.

#' Scale-free gain rule for choosing the number of factors
#'
#' Selects a number of factors from a criterion curve by thresholding the
#' marginal gain of each added factor against the largest gain in the candidate
#' range. Because the threshold is a fraction of the maximum gain, the rule is
#' dimensionless and invariant to sample size, item count, and window.
#'
#' A sustained-drop guard prevents a transient one-step dip from firing the
#' elbow prematurely: a step counts as the elbow only if the next
#' `sustain - 1` steps are also below threshold (or it is the last step).
#' `sustain = 1` reproduces a first-crossing rule.
#'
#' @param K Candidate numbers of factors (any order).
#' @param score A criterion for which larger is better. For an information
#'   criterion (smaller is better), pass `score = -IC`.
#' @param delta Threshold as a percentage (0-100) of the largest marginal gain
#'   (default `10`).
#' @param sustain Consecutive sub-threshold steps required to call the elbow
#'   (default `2`).
#'
#' @return The selected number of factors (a scalar from `K`).
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' Chen, J. (2023). Fully and partially exploratory factor analysis with
#' bi-level Bayesian regularization. *Behavior Research Methods*, 55(4),
#' 2125-2142. \doi{10.3758/s13428-022-01884-7}
#'
#' @examples
#' ## ELBOs from a window sweep over K = 2:6 -- large gains up to K = 4,
#' ## then only dust. The rule fires at the elbow:
#' elbo <- c(-5400, -5100, -4900, -4895, -4892)
#' select_K_elbow(K = 2:6, score = elbo)          # 4
#'
#' ## for an information criterion (smaller is better), pass its negative
#' bic <- c(10900, 10500, 10310, 10395, 10480)
#' select_K_elbow(K = 2:6, score = -bic)          # 4
#'
#' ## sustain = 1 is a plain first-crossing rule; the default sustain = 2
#' ## refuses to fire on a single transient dip
#' jag <- c(-5400, -5100, -5095, -4900, -4897, -4895)
#' select_K_elbow(K = 2:7, score = jag, sustain = 1)  # 3 (fooled by the dip)
#' select_K_elbow(K = 2:7, score = jag, sustain = 2)  # 5 (rides through it)
#'
#' @seealso [pefa()]
#' @export
select_K_elbow <- function(K, score, delta = 10, sustain = 2) {
  ## ---- input contract -------------------------------------------------
  if (length(K) != length(score))
    stop("K and score must have the same length (got ", length(K), " and ",
         length(score), ").", call. = FALSE)
  if (!length(K)) stop("K must not be empty.", call. = FALSE)
  if (!is.numeric(K) || anyNA(K) || any(!is.finite(K)))
    stop("K must be finite and numeric.", call. = FALSE)
  if (any(K != round(K)))
    stop("K must contain whole numbers.", call. = FALSE)
  if (anyDuplicated(K))
    stop("K must not contain duplicates.", call. = FALSE)
  if (!is.numeric(score) || anyNA(score) || any(!is.finite(score)))
    stop("score must be finite and numeric; drop or impute non-finite ",
         "candidates before selecting.", call. = FALSE)
  if (length(delta) != 1L || !is.finite(delta) || delta < 0 || delta > 100)
    stop("delta must be a single number in [0, 100].", call. = FALSE)
  if (length(sustain) != 1L || !is.finite(sustain) || sustain < 1 ||
      sustain != round(sustain))
    stop("sustain must be a single positive whole number.", call. = FALSE)

  o <- order(K); K <- K[o]; score <- score[o]
  if (length(K) == 1L) return(K)
  g <- diff(score)                       # marginal gain of adding each factor
  gmax <- max(g)
  if (gmax <= 0) return(K[1])            # nothing worth adding
  thr <- (delta / 100) * gmax
  below <- g < thr; n <- length(g)
  for (i in seq_len(n)) {                # first SUSTAINED sub-threshold run
    win <- i:min(i + sustain - 1L, n)
    if (all(below[win])) return(K[i])
  }
  K[length(K)]                           # never sustained -> upper bound
}

#' Partially exploratory factor analysis with gain-rule factor selection
#'
#' Fits [vbfa()] over a window of factor numbers, holding a fixed partially
#' confirmatory backbone `Q0` and padding each model with fully exploratory
#' columns, then selects the number of factors with the scale-free gain rule
#' ([select_K_elbow()]). Raw ELBO and BIC optima are reported for reference.
#'
#' @param Q0 The `J x K0` backbone design matrix (`1`/`0`/`-1` as in [vbfa()]);
#'   `K0 = ncol(Q0)` specified factors.
#' @param Y The `N x J` data matrix.
#' @param Kmin,Kmax Inclusive factor-number window; `Kmin >= K0`.
#' @param delta **Deprecated**; use `cuts`. A scalar maps to
#'   `cuts = c(primary = delta)`. Supplying both is an error.
#' @param sustain Consecutive sub-threshold steps required to call the elbow.
#' @param cuts Named numeric vector of gain thresholds, as percentages of the
#'   largest usable marginal gain. Exactly one entry must be named `primary`
#'   (it drives the returned `selected_K`); any further named entries are
#'   evaluated as sensitivity rules and appear in `$selection`. Default
#'   `c(primary = 10)`. The manuscript workflow uses
#'   `cuts = c(primary = 20, sensitivity = 10)`; that 20% figure is calibrated
#'   within its own populations and is therefore *not* the package default.
#' @param keep_fits Which candidate fits to retain in `$fits`: `"all"`
#'   (default), `"selected"`, or `"none"`. `selected_fit` is returned in every
#'   case; the candidate slots exist for *inspection*, which summary columns
#'   cannot replace.
#' @param fit_cut Logical; if `TRUE`, gate selection on absolute fit passing
#'   `cutoffs`. Default `FALSE`.
#' @param cutoffs Named absolute-fit cutoffs used when `fit_cut = TRUE`.
#' @param max_it,tau Passed to [vbfa()] / [fit_stats()] (`max_it` is per stage).
#' @param orthogonal Passed to [vbfa()] and [fit_stats()].
#' @param bifactor Logical. If `TRUE`, sweep **bifactor** candidates with
#'   `Q0` supplying only the group backbone: every candidate gets a general
#'   column internally (see `general`), and `Kmin`/`Kmax` and the reported
#'   `K` count **group factors** -- directly comparable to an oblique sweep
#'   over the same window. The sweep table gains a `K_total` column
#'   (`K + 1`). Implies orthogonal factors; `bifactor = TRUE` overrides
#'   `orthogonal = FALSE` (with a message). The identification warning below
#'   applies with full force.
#' @param general Design of the general column when `bifactor = TRUE`;
#'   passed to [vbfa()] (scalar or length-`J` vector of `1`/`0`/`-1`,
#'   default `1`).
#' @param verbose Logical; print per-model progress.
#' @param ... Further arguments passed to [vbfa()] (e.g. `v0`, `Qe`,
#'   `ld_control`, `tolVal`).
#'
#' @return An object of class `c("pefa", "vbpm_sweep")`. Everything a sweep
#'   produces is indexed by one of three things, and each has its own table:
#'   \describe{
#'     \item{`$sweep`}{One row per candidate: `ELBO`, `Objective`, `AIC`,
#'       `BIC`, the fit indices, `t`, `iter`, `flag`, `secs`, `converged`,
#'       `pass_fit`, `eligible`. Properties of a single fitted model.}
#'     \item{`$transitions`}{One row per **adjacent pair** (`W - 1` rows):
#'       the marginal gains (`Objective_gain`, `ELBO_gain`, `BIC_gain`,
#'       oriented so larger is better) together with the adjacent-count
#'       stability indexes `phi_min`, `phi_mean`, `n_phi90`, `ari`, `rmsd`,
#'       `surplus_max`, the four assignment-coverage counts, and
#'       `stability_status`. Gains live here, not in `$sweep`: they are
#'       properties of a step, not of a candidate.}
#'     \item{`$selection`}{One row per rule -- criterion x form x cut x
#'       eligibility scope -- with the normalizer `gain_max`, the absolute
#'       `threshold`, `selected_K`, and both boundary readouts.}
#'   }
#'   Plus `selected_K`, `selected_fit`, `fit_stats`, `Q_selected`, `$fits`
#'   (per `keep_fits`), `$cuts`, `window_boundary`, `support_boundary`, and
#'   the normalized `orthogonal`/`bifactor`/`ld` settings. `$boundary` and
#'   `$delta` remain as deprecated aliases for the 0.8.x line.
#'
#'   No percentage-of-largest-gain column is stored: the normalizer depends on
#'   the eligibility scope, so percentages are computed for display from
#'   `$transitions` and the matching `$selection` row.
#'
#'   Use [summary.pefa()] for the aggregate view, [selected_fit()] for the
#'   selected model, and [transition_detail()] for the per-column material
#'   behind one transition row.
#'
#' @section Orthogonal (bifactor) sweeps are weakly identified:
#' With `bifactor = TRUE` (or, equivalently, `orthogonal = TRUE` and a
#' hand-built general column in the backbone), do **not**
#' expect a sweep to discover the number of group factors. An omitted
#' orthogonal group factor contributes a rank-one within-cluster covariance
#' block that the unspecified (`-1`) entries of the *remaining* columns can
#' reproduce almost exactly, so the smaller candidate is covariance-equivalent
#' to the truth: fit indices sit at their ceilings for every candidate, the
#' objective typically *falls* with `K`, and the selection lands on the lower
#' window boundary regardless of the true `K`. In simulations from true
#' bifactor and true higher-order populations this failure is systematic (it
#' is an identification property, not sampling noise), and the truth-sized
#' candidate inside such a sweep is itself rotationally scrambled -- its
#' loadings must not be interpreted. Encoding every loading that is *known*
#' to be absent as a fixed `0` rather than `-1` restores an interior
#' selection for the factor count, but the sweep's fits remain
#' uninterpretable; all interpretation belongs to an anchored refit.
#'
#' The decision route that works runs through the *ordinary* sweep instead:
#' (1) sweep **without** a general column (`orthogonal = FALSE`) for the
#' first-order factor count `K` -- for that question an omitted factor
#' produces honest misfit, and interior selections are trustworthy even when
#' the data are truly bifactor or higher-order (the general factor is then
#' absorbed into the factor correlations, which is harmless at this step);
#' (2) at the selected `K`, fit the anchored oblique `K`-factor model and
#' the anchored orthogonal bifactor with one added, fully specified general
#' column, and compare their `BIC` from [fit_stats()] -- compare `BIC`, not
#' `ELBO`, because the two designs carry spike-and-slab prior mass on
#' different numbers of entries; (3) if the bifactor wins, check
#' Schmid-Leiman proportionality to decide between a genuine bifactor and a
#' higher-order collapse; if the oblique model wins, a higher-order layer
#' may still be read from the factor correlations. The `bifactor` vignette
#' works through the failure, the route, and its power limits (a true
#' bifactor with near-homogeneous loadings needs a large sample at step 2).
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' Chen, J. (2023). Fully and partially exploratory factor analysis with
#' bi-level Bayesian regularization. *Behavior Research Methods*, 55(4),
#' 2125-2142. \doi{10.3758/s13428-022-01884-7}
#'
#' @examples
#' \donttest{
#' ## three-factor data; the backbone Q0 anchors only the first two factors,
#' ## and the sweep asks how many factors the data actually support
#' sim <- sim_fa(N = 400, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' J   <- ncol(sim$dat)
#' Q0  <- matrix(-1L, J, 2)
#' for (k in 1:2) { a <- which(rep(1:3, each = 6) == k)[1:2]
#'                  Q0[a, ] <- 0L; Q0[a, k] <- 1L }
#'
#' r <- pefa(Q0, sim$dat, Kmin = 2, Kmax = 4, verbose = FALSE)
#' summary(r)            # objective, objective-gain, BIC and BIC-gain selections
#' r$sweep               # the full sweep table (ELBO, gains, fit, timing)
#' selected_fit(r)$flag  # the refitted selected model is returned
#'
#' ## bifactor sweep: Q0 is the GROUP backbone and the window counts GROUP
#' ## factors (the general column is added internally to every candidate).
#' ## Read the identification warning above before interpreting one.
#' rb <- pefa(Q0, sim$dat, Kmin = 2, Kmax = 4, bifactor = TRUE, v0 = .001,
#'            verbose = FALSE)
#' rb$sweep[, c("K", "K_total", "BIC")]
#' }
#'
#' @seealso [vbfa()], [select_K_elbow()]
#' @export
pefa <- function(Q0, Y, Kmin, Kmax, delta = NULL, sustain = 2, fit_cut = FALSE,
                    cutoffs = c(RMSEA = .06, SRMR = .10, CFI = .90, TLI = .90),
                    max_it = 10000, tau = 0.50,
                    orthogonal = FALSE, bifactor = FALSE, general = 1,
                    verbose = TRUE, ..., cuts = c(primary = 10),
                    keep_fits = c("all", "selected", "none")) {
  call <- match.call()
  cuts <- .pefa_cuts(cuts, delta, !missing(cuts), !missing(delta))
  keep_fits <- match.arg(keep_fits)
  Y  <- as.matrix(Y); Q0 <- as.matrix(Q0)
  J  <- ncol(Y); K0 <- ncol(Q0)
  if (nrow(Q0) != J) stop("nrow(Q0) must equal ncol(Y).", call. = FALSE)
  ## bifactor sweep: Q0 is the GROUP backbone, the window counts GROUP
  ## factors, and every candidate gets the general column from vbfa()
  ## internally. bifactor = TRUE overrides orthogonal = FALSE.
  bifactor <- isTRUE(bifactor)
  if (bifactor) {
    if (any(apply(Q0 == 1, 2, all)))
      stop("Q0 contains an all-specified column. With bifactor = TRUE, ",
           "supply only the GROUP backbone in Q0; the general column is ",
           "added internally (see the `general` argument).", call. = FALSE)
    if (!missing(orthogonal) && !isTRUE(orthogonal))
      message("bifactor = TRUE implies orthogonal factors; ",
              "orthogonal = FALSE is overridden.")
    orthogonal <- TRUE
  }
  ## Validate the window before fitting: Kmin:Kmax silently reverses when
  ## Kmax < Kmin, which would sweep a different set than the user asked for.
  for (nm in c("Kmin", "Kmax")) {
    v <- get(nm)
    if (length(v) != 1L || !is.finite(v) || v != round(v) || v < 1)
      stop(nm, " must be a single positive whole number.", call. = FALSE)
  }
  if (Kmax < Kmin)
    stop(sprintf("Kmax (%d) < Kmin (%d): the window must be non-decreasing.",
                 Kmax, Kmin), call. = FALSE)
  if (Kmin < K0)
    stop(sprintf("Kmin (%d) < K0 (%d): the window must be at least the number of specified factors.",
                 Kmin, K0), call. = FALSE)
  if (length(tau) != 1L || !is.finite(tau) || tau < 0 || tau > 1)
    stop("tau must be a single number in [0, 1].", call. = FALSE)
  Ks  <- Kmin:Kmax
  pad <- function(K) if (K == K0) Q0 else cbind(Q0, matrix(-1L, J, K - K0))

  fits  <- setNames(vector("list", length(Ks)), as.character(Ks))
  rows  <- vector("list", length(Ks))

  for (i in seq_along(Ks)) {
    K <- Ks[i]
    t0  <- proc.time()[["elapsed"]]
    fit <- if (bifactor)
      vbfa(Y, pad(K), max_it = max_it, convChk = FALSE,
           bifactor = TRUE, general = general, ...)
    else
      vbfa(Y, pad(K), max_it = max_it, convChk = FALSE,
           orthogonal = orthogonal, ...)
    fs  <- fit_stats(fit, Y, fit$Q, tau = tau, orthogonal = orthogonal)
    fits[[i]] <- fit
    rows[[i]] <- data.frame(
      K = K, K_total = K + as.integer(bifactor),
      ELBO = as.numeric(fit$ELBO),
      Objective = as.numeric(fit$objective),
      objective_type = fit$objective_type,
      AIC = fs["AIC"], BIC = fs["BIC"],
      RMSEA = fs["RMSEA"], SRMR = fs["SRMR"], CFI = fs["CFI"], TLI = fs["TLI"],
      t = fs["t"], iter = fit$iter, flag = fit$flag,
      secs = max(0, proc.time()[["elapsed"]] - t0), row.names = NULL)
    if (verbose) cat(sprintf("K=%2d  ELBO=%.1f  BIC=%.0f  CFI=%.3f  conv=%d  %.0fs\n",
                             K, rows[[i]]$ELBO, rows[[i]]$BIC, rows[[i]]$CFI,
                             as.integer(fit$flag == 1), rows[[i]]$secs))
  }
  sweep <- do.call(rbind, rows); sweep <- sweep[order(sweep$K), ]

  types <- unique(sweep$objective_type)
  if (length(types) != 1L)
    stop("PEFA candidates have mixed objective types: ",
         paste(types, collapse = ", "), call. = FALSE)

  ## ---- eligibility: converged within max_it AND passes absolute fit -------
  sweep$converged <- sweep$flag == 1L
  sweep$pass_fit  <- if (!isTRUE(fit_cut)) rep(TRUE, nrow(sweep)) else with(sweep,
      !is.na(RMSEA) & RMSEA <= cutoffs[["RMSEA"]] &
      !is.na(SRMR)  & SRMR  <= cutoffs[["SRMR"]] &
      !is.na(CFI)   & CFI   >= cutoffs[["CFI"]] &
      !is.na(TLI)   & TLI   >= cutoffs[["TLI"]])
  sweep$eligible  <- sweep$converged & sweep$pass_fit

  ## ---- edge table: gains plus adjacent-count stability --------------------
  ## Loading-based stability is meaningful only for oblique sweeps; under
  ## orthogonal/bifactor rotation the padded columns are indeterminate.
  transitions <- .build_transitions(sweep, fits, unsupported = orthogonal)

  ## ---- rule table ---------------------------------------------------------
  selection <- .pefa_selection(sweep, transitions, cuts, sustain, types)
  prim  <- .pefa_primary_row(selection)
  K_sel <- selection$selected_K[prim]

  sel_fit <- sel_stats <- NULL
  if (!is.na(K_sel)) {
    sel_fit <- fits[[as.character(K_sel)]]
    sel_stats <- fit_stats(sel_fit, Y, sel_fit$Q, tau = tau, orthogonal = orthogonal)
  }
  ## keep_fits governs the candidate slots only; selected_fit is always returned
  if (keep_fits != "all") {
    keep <- if (keep_fits == "selected" && !is.na(K_sel)) as.character(K_sel) else character(0)
    for (nm in setdiff(names(fits), keep)) fits[[nm]] <- list(NULL)[[1]]
  }

  if (verbose) {
    n_fail <- sum(!sweep$converged); n_unfit <- sum(sweep$converged & !sweep$pass_fit)
    gate <- if (isTRUE(fit_cut)) "ON" else "OFF"
    cat(sprintf("\n--- selection (primary cut = %g%%; convergence gate ON, absolute-fit gate %s) ---\n",
                cuts[["primary"]], gate))
    cat(sprintf("  excluded: %d non-converged, %d converged-but-failing-fit; %d eligible\n",
                n_fail, n_unfit, sum(sweep$eligible)))
    if (is.na(K_sel)) cat("  NO eligible model -- selection is NA (relax fit_cut or extend/converge).\n")
    else cat(sprintf("  %s gain -> K = %d   (selected; model returned)\n", types, K_sel))
    un <- attr(selection, "unavailable")
    if (length(un))
      cat("  unavailable rules: ", paste(un, collapse = ", "), "\n", sep = "")
  }
  out <- list(call = call, object_schema_version = 2L,
       selected_K = K_sel, selection = selection,
       selected_fit = sel_fit, fits = fits,
       fit_stats = sel_stats,
       Q_selected = if (is.na(K_sel)) NULL else
                    if (bifactor) sel_fit$Q else pad(K_sel),
       sweep = sweep, transitions = transitions,
       window = list(K = Ks, Kmin = min(Ks), Kmax = max(Ks), K0 = K0),
       K0 = K0, cuts = cuts, sustain = sustain,
       orthogonal = orthogonal, bifactor = bifactor, ld = isTRUE(list(...)$ld),
       n_general = if (bifactor) 1L else 0L, keep_fits = keep_fits,
       window_boundary  = selection$window_boundary[prim],
       support_boundary = selection$support_boundary[prim],
       objective_type = types,
       fit_cut = isTRUE(fit_cut), cutoffs = cutoffs)
  ## deprecated 0.8.x aliases; methods read $cuts and $window_boundary
  out$delta    <- unname(cuts[["primary"]])
  out$boundary <- out$window_boundary
  class(out) <- c("pefa", "vbpm_sweep")
  out
}


#' Methods for a PEFA factor-number sweep
#'
#' A [pefa()] result is a sweep over a candidate factor-number window, not a
#' separate fitted probability model. `summary()` aggregates the candidate
#' fits and reports the raw objective, objective-gain, raw BIC, and BIC-gain
#' selections. The objective is the ELBO for diagonal VBFA sweeps and the
#' terminal VECM objective for LD sweeps; ELBO-specific columns are reported
#' only when every candidate has `objective_type = "elbo"`.
#' [selected_fit()] extracts the actual fitted [vbfa()] object selected by the
#' primary objective-gain rule.
#'
#' @param object,x A result returned by [pefa()].
#' @param ... Further arguments. For `plot()`, graphical arguments passed to
#'   base plotting functions.
#' @param type Plot type. `"objective"` draws the objective (left axis) and
#'   BIC (right axis) on their own scales; `"gain"` draws the marginal gains
#'   as a percentage of the largest gain, with the `delta` threshold
#'   (gains are undefined at `Kmin`, so those curves start one step into the
#'   window); `"fit"` draws absolute misfit (RMSEA/SRMR) and incremental fit
#'   (CFI/TLI) in separate panels with conventional cutoffs. All types mark
#'   the selected `K` with a dotted vertical line.
#' @param digits Number of decimals for the fit indices in the printed
#'   comparison table (gain percentages always print with one decimal).
#' @return `summary.pefa()` returns an object of class `summary.pefa` carrying
#'   the sweep's three tables under their own names (`sweep`, `transitions`,
#'   `selection`) plus window and tuning metadata, convergence counts, and both
#'   boundary readouts. It deliberately carries **no fit payload**: use
#'   [selected_fit()] on the `pefa` object for the fitted model, so summaries
#'   stay small and printable. Print methods return their input invisibly;
#'   `plot.pefa()` returns `x` invisibly.
#'
#'   Displayed gain percentages are computed on demand, normalized by the
#'   `converged+fit` scope's `gain_max`, and are `NA` for an edge that scope
#'   cannot use; no percentage column is stored.
#' @name pefa-methods
NULL

#' @rdname pefa-methods
#' @export
selected_fit <- function(object, ...) UseMethod("selected_fit")

#' @rdname pefa-methods
#' @export
selected_fit.pefa <- function(object, ...) object$selected_fit

#' @rdname pefa-methods
#' @export
print.pefa <- function(x, ...) {
  nconv <- sum(x$sweep$converged)
  klab  <- if (isTRUE(x$bifactor)) " group factors (+ 1 general)" else ""
  cat(sprintf("PEFA sweep: K = %d:%d%s (%d of %d converged)
",
              x$window$Kmin, x$window$Kmax, klab, nconv, nrow(x$sweep)))
  cat(sprintf("  %s gain (primary cut = %g%%) -> K = %s
",
              toupper(x$objective_type), x$cuts[["primary"]],
              if (is.na(x$selected_K)) "NA" else x$selected_K))
  cat("  objective:", x$objective_type, "
")
  .pefa_boundary_advice(x$window_boundary, x$support_boundary,
                        isTRUE(x$orthogonal), x$sweep, x$window)
  cat("Use summary() for all selections and the per-K comparison.
")
  invisible(x)
}

## Shared boundary advice. The requested window and the usable support are
## different facts: extending helps only the first.
.pefa_boundary_advice <- function(wb, sb, orth, sweep, window) {
  if (identical(wb, "lower") || identical(wb, "upper")) {
    if (orth && identical(wb, "lower")) {
      cat("  boundary selection (orthogonal sweep): WEAK EVIDENCE about K --
",
          "  an omitted group factor is absorbed by unspecified entries of the
",
          "  remaining columns (covariance-equivalent candidates). Encode known
",
          "  zeros as 0, not -1; see ?pefa and vignette('bifactor').
", sep = "")
    } else {
      direction <- if (identical(wb, "lower")) "downward" else "upward"
      cat(sprintf("  window boundary: extend the window %s
", direction))
    }
  }
  if (sb %in% c("lower", "upper", "single") && !identical(wb, sb)) {
    excl <- sweep$K[!sweep$eligible]
    cat(sprintf("  support boundary (%s): the selection sits at the edge of the
", sb),
        sprintf("  usable candidates%s. Address convergence or the absolute-fit
",
                if (length(excl)) paste0(" (excluded K: ", paste(excl, collapse = ", "), ")") else ""),
        "  gate rather than extending the window.
", sep = "")
  }
  invisible(NULL)
}

#' @rdname pefa-methods
#' @export
summary.pefa <- function(object, ...) {
  z <- list(
    call = object$call,
    object_schema_version = object$object_schema_version,
    orthogonal = isTRUE(object$orthogonal),
    bifactor = isTRUE(object$bifactor),
    ld = isTRUE(object$ld),
    window = object$window,
    cuts = object$cuts,
    sustain = object$sustain,
    sweep = object$sweep,
    transitions = object$transitions,
    selection = object$selection,
    selected_K = object$selected_K,
    selected_rule = "objective gain at the primary cut (converged+fit)",
    window_boundary = object$window_boundary,
    support_boundary = object$support_boundary,
    objective_type = object$objective_type,
    convergence = c(converged = sum(object$sweep$converged),
                    total = nrow(object$sweep)))
  z$delta <- object$delta          # deprecated alias
  z$boundary <- object$window_boundary
  class(z) <- "summary.pefa"
  z
}

## Percentage display rule, shared by print and plot: normalize by the
## converged+fit gain_max for that criterion, and show NA where that
## denominator is absent, non-positive, or did not apply to the edge.
.pefa_gain_pct <- function(x, criterion) {
  tr <- x$transitions
  sel <- x$selection
  i <- which(sel$criterion == criterion & sel$form == "gain" &
             sel$eligibility_scope == "converged+fit")[1L]
  gmax <- if (length(i) && !is.na(i)) sel$gain_max[i] else NA_real_
  g <- tr[[paste0(if (criterion == "bic") "BIC" else
                  if (criterion == "elbo") "ELBO" else "Objective", "_gain")]]
  if (!length(g) || !is.finite(gmax) || gmax <= 0) return(rep(NA_real_, length(g)))
  ok <- x$sweep$eligible[match(tr$K_from, x$sweep$K)] &
        x$sweep$eligible[match(tr$K_to,   x$sweep$K)] & is.finite(g)
  ifelse(ok, 100 * g / gmax, NA_real_)
}

#' @rdname pefa-methods
#' @export
print.summary.pefa <- function(x, digits = 3, ...) {
  klab <- if (isTRUE(x$bifactor)) " group factors (+ 1 general)" else ""
  cat(sprintf("PEFA sweep summary: K = %d:%d%s; primary cut = %g%%; sustain = %d
",
              x$window$Kmin, x$window$Kmax, klab, x$cuts[["primary"]], x$sustain))
  sel <- x$selection
  show <- sel[sel$eligibility_scope == "converged+fit", , drop = FALSE]
  for (i in seq_len(nrow(show))) {
    lab <- if (show$form[i] == "raw") sprintf("%s raw", toupper(show$criterion[i]))
           else sprintf("%s gain [%s = %g%%]", toupper(show$criterion[i]),
                        show$cut[i], show$cut_value[i])
    cat(sprintf("  %-28s -> K = %s
", lab,
                if (is.na(show$selected_K[i])) "NA" else show$selected_K[i]))
  }
  cat(sprintf("  converged: %d of %d; window boundary: %s; support boundary: %s
",
              x$convergence["converged"], x$convergence["total"],
              x$window_boundary, x$support_boundary))
  if (isTRUE(x$orthogonal) && identical(x$window_boundary, "lower"))
    cat("  (orthogonal sweep at the lower boundary: weak evidence about K --
",
        "   see ?pefa, section on orthogonal sweeps, and vignette('bifactor'))
",
        sep = "")

  cat("
Candidates:
")
  d <- x$sweep; lab <- toupper(x$objective_type)
  tab <- data.frame(K = d$K)
  if (isTRUE(x$bifactor) && "K_total" %in% names(d)) tab$K_total <- d$K_total
  tab[[lab]] <- round(d$Objective, 1)
  tab$BIC <- round(d$BIC, 1)
  for (nm in c("RMSEA", "SRMR", "CFI", "TLI"))
    if (nm %in% names(d)) tab[[nm]] <- round(d[[nm]], digits)
  tab$converged <- d$converged
  tab$` ` <- ifelse(!is.na(x$selected_K) & d$K == x$selected_K, "<- selected", "")
  print(tab, row.names = FALSE)

  tr <- x$transitions
  if (nrow(tr)) {
    cat("
Transitions (gain % of the largest usable gain, converged+fit scope):
")
    tt <- data.frame(edge = paste(tr$K_from, "->", tr$K_to))
    tt[[paste0(lab, "_gain")]] <- round(tr$Objective_gain, 1)
    tt[[paste0(lab, "_gain%")]] <- round(.pefa_gain_pct(x, "objective"), 1)
    tt$`BIC_gain%` <- round(.pefa_gain_pct(x, "bic"), 1)
    tt$phi_min <- round(tr$phi_min, 3)
    tt$ari <- round(tr$ari, 3)
    tt$surplus_max <- round(tr$surplus_max, 3)
    tt$status <- tr$stability_status
    print(tt, row.names = FALSE)
  }
  cat("(Raw gains and all stability columns are in $transitions;",
      "every rule is in $selection.)
")
  invisible(x)
}

#' @rdname pefa-methods
#' @importFrom graphics abline axis legend lines mtext par
#' @export
plot.pefa <- function(x, type = c("objective", "gain", "fit"), ...) {
  type <- match.arg(type)
  d    <- x$sweep
  K    <- d$K
  Ksel <- x$selected_K
  lab  <- toupper(x$objective_type)
  col1 <- "#1f77b4"; col2 <- "#d62728"   # objective/absolute vs BIC/incremental
  ## headroom at the top of a y-range so the legend never sits on a curve
  pad <- function(v, top = 0.30) {
    r <- range(v, na.rm = TRUE); s <- diff(r)
    if (s == 0) s <- abs(r[1]) * 0.01 + 1
    c(r[1] - 0.05 * s, r[2] + top * s)
  }
  mark_sel <- function() if (!is.na(Ksel)) abline(v = Ksel, lty = 3, col = "grey40")

  if (type == "objective") {
    ## The objective (higher is better) and BIC (lower is better) live on
    ## different scales; forcing them onto one axis flattens both curves,
    ## so each gets its own axis.
    op <- par(mar = c(5, 4.2, 4, 4.2) + 0.1); on.exit(par(op))
    plot(K, d$Objective, type = "b", pch = 16, lty = 1, col = col1,
         ylim = pad(d$Objective), xaxt = "n", xlab = "Number of factors (K)",
         ylab = paste0(lab, " (higher is better)"),
         main = "PEFA objective sweep", ...)
    axis(1, at = K)
    mark_sel()
    par(new = TRUE)
    plot(K, d$BIC, type = "b", pch = 17, lty = 2, col = col2,
         ylim = pad(d$BIC), axes = FALSE, xlab = "", ylab = "")
    axis(4, col = col2, col.axis = col2)
    mtext("BIC (lower is better)", side = 4, line = 2.8, col = col2)
    legend("top", c(paste(lab, "(left axis)"), "BIC (right axis)",
                    if (!is.na(Ksel)) sprintf("selected K = %d", Ksel)),
           col = c(col1, col2, if (!is.na(Ksel)) "grey40"),
           lty = c(1, 2, if (!is.na(Ksel)) 3),
           pch = c(16, 17, if (!is.na(Ksel)) NA),
           ncol = 2, bty = "n", cex = 0.9)
  } else if (type == "gain") {
    ## Edge quantities: plotted at K_to so each point sits on the candidate the
    ## step arrives at. Percentages use the converged+fit normalizer (see
    ## ?pefa-methods); unusable edges are gaps, not zeros.
    tr <- x$transitions
    og <- .pefa_gain_pct(x, "objective"); bg <- .pefa_gain_pct(x, "bic")
    cutv <- x$cuts
    ylim <- pad(c(og, bg, unname(cutv), 0))
    plot(tr$K_to, og, type = "b", pch = 16, lty = 1, col = col1,
         ylim = ylim, xlim = range(K), xaxt = "n",
         xlab = "Number of factors (K)",
         ylab = "Marginal gain (% of largest usable gain)",
         main = "PEFA gain sweep", ...)
    axis(1, at = K)
    lines(tr$K_to, bg, type = "b", pch = 17, lty = 2, col = col2)
    abline(h = 0, col = "grey85")
    for (i in seq_along(cutv)) abline(h = cutv[[i]], lty = 3)
    mark_sel()
    legend("top", c(paste(lab, "gain"), "BIC gain",
                    paste0(names(cutv), " = ", unname(cutv), "%"),
                    if (!is.na(Ksel)) sprintf("selected K = %d", Ksel)),
           col = c(col1, col2, rep("black", length(cutv)),
                   if (!is.na(Ksel)) "grey40"),
           lty = c(1, 2, rep(3, length(cutv)), if (!is.na(Ksel)) 3),
           pch = c(16, 17, rep(NA, length(cutv)), if (!is.na(Ksel)) NA),
           ncol = 2, bty = "n", cex = 0.9)
  } else {
    ## Absolute misfit (RMSEA/SRMR, near 0) and incremental fit (CFI/TLI,
    ## near 1) are an order of magnitude apart; one panel each, with the
    ## conventional cutoffs drawn for orientation.
    op <- par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1) + 0.1,
              oma = c(0, 0, 2, 0)); on.exit(par(op))
    plot(K, d$RMSEA, type = "b", pch = 16, lty = 1, col = col1,
         ylim = pad(c(0, d$RMSEA, d$SRMR, 0.1), top = 0.35), xaxt = "n",
         xlab = "Number of factors (K)", ylab = "Index (lower is better)",
         main = "Absolute misfit", ...)
    axis(1, at = K)
    lines(K, d$SRMR, type = "b", pch = 17, lty = 2, col = col2)
    abline(h = c(0.06, 0.08), lty = 3, col = "grey60")
    mark_sel()
    legend("top", c("RMSEA (cutoff .06)", "SRMR (cutoff .08)"),
           col = c(col1, col2), lty = 1:2, pch = c(16, 17), bty = "n", cex = 0.9)
    plot(K, d$CFI, type = "b", pch = 16, lty = 1, col = col1,
         ylim = c(min(c(d$CFI, d$TLI, 0.9), na.rm = TRUE) - 0.02, 1.02),
         xaxt = "n", xlab = "Number of factors (K)",
         ylab = "Index (higher is better)", main = "Incremental fit", ...)
    axis(1, at = K)
    lines(K, d$TLI, type = "b", pch = 17, lty = 2, col = col2)
    abline(h = 0.95, lty = 3, col = "grey60")
    mark_sel()
    legend("bottom", c("CFI (cutoff .95)", "TLI (cutoff .95)"),
           col = c(col1, col2), lty = 1:2, pch = c(16, 17), bty = "n", cex = 0.9)
    mtext("PEFA fit-index sweep", outer = TRUE, font = 2)
  }
  invisible(x)
}

#' Per-column detail behind one transition
#'
#' Returns the ragged material summarized by one row of a sweep's
#' `$transitions` table: the per-column congruences, which larger-solution
#' column each smaller-solution column was matched to, the sign multipliers,
#' and the unmatched (surplus) column with its full loading vector.
#'
#' The accessor never bypasses the gates of [pefa()]: a transition inside an
#' orthogonal/bifactor sweep, or one with a non-converged endpoint, returns
#' `NULL` with the reason, because those loadings are not interpretable even
#' when the fits are retained. It also returns `NULL` when `keep_fits`
#' discarded either endpoint.
#'
#' @param object A [pefa()] result.
#' @param K_from The lower factor number of the transition (`K_from -> K_from + 1`).
#' @param ... Unused.
#'
#' @return A list with `K_from`, `K_to`, `phi` (per matched column),
#'   `matched` (larger-solution column indices), `sign`, `surplus` (index) and
#'   `surplus_loadings`, plus the tolerances used; or `NULL` when detail is
#'   unavailable.
#' @seealso [pefa()]
#' @export
transition_detail <- function(object, K_from, ...) UseMethod("transition_detail")

#' @rdname transition_detail
#' @export
transition_detail.pefa <- function(object, K_from, ...) {
  tr <- object$transitions
  i  <- which(tr$K_from == K_from)
  if (!length(i))
    stop("K_from = ", K_from, " is not a transition in this sweep (K_from in ",
         paste(tr$K_from, collapse = ", "), ").", call. = FALSE)
  if (!identical(tr$stability_status[i], "ok") &&
      grepl("unsupported|not converged|not retained", tr$stability_status[i])) {
    message("No transition detail: ", tr$stability_status[i], ".")
    return(invisible(NULL))
  }
  det <- attr(tr, "detail")[[i]]
  if (is.null(det)) { message("No transition detail: fit not retained."); return(invisible(NULL)) }
  B <- object$fits[[as.character(tr$K_to[i])]]$Lam
  list(K_from = tr$K_from[i], K_to = tr$K_to[i],
       phi = det$phi, matched = det$idx, sign = det$sign,
       surplus = det$surplus,
       surplus_loadings = if (length(det$surplus)) B[, det$surplus] else NULL,
       tolerances = .pefa_tol)
}
