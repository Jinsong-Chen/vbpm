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
#' @param delta Threshold as a percentage (0–100) of the largest marginal gain
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
#' @seealso [pefa()]
#' @export
select_K_elbow <- function(K, score, delta = 10, sustain = 2) {
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
#' @param delta,sustain Passed to [select_K_elbow()].
#' @param fit_cut Logical; if `TRUE`, gate selection on absolute fit passing
#'   `cutoffs`. Default `FALSE`.
#' @param cutoffs Named absolute-fit cutoffs used when `fit_cut = TRUE`.
#' @param max_it,rseed,tau Passed to [vbfa()] / [vbfa_fit()] (`max_it` is per
#'   stage).
#' @param orthogonal Passed to [vbfa()] and [vbfa_fit()].
#' @param save_path Optional CSV checkpoint path (resume-safe).
#' @param verbose Logical; print per-model progress.
#' @param ... Further arguments passed to [vbfa()] (e.g. `v0`, `Qe`,
#'   `ld_control`, `tolVal`).
#'
#' @return A list with `selected_K`, the `selection` vector, the fitted model
#'   and `fit_stats` at the selected `K`, `Q_selected`, and the full `sweep`
#'   data frame (carrying `converged`, `pass_fit`, `eligible` flags).
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
#' @seealso [vbfa()], [select_K_elbow()]
#' @export
pefa <- function(Q0, Y, Kmin, Kmax, delta = 10, sustain = 2, fit_cut = FALSE,
                    cutoffs = c(RMSEA = .06, SRMR = .10, CFI = .90, TLI = .90),
                    max_it = 10000, rseed = 12345, tau = 0.50,
                    orthogonal = FALSE, save_path = NULL, verbose = TRUE, ...) {
  Y  <- as.matrix(Y); Q0 <- as.matrix(Q0)
  J  <- ncol(Y); K0 <- ncol(Q0)
  if (nrow(Q0) != J) stop("nrow(Q0) must equal ncol(Y).", call. = FALSE)
  if (Kmin < K0)
    stop(sprintf("Kmin (%d) < K0 (%d): the window must be at least the number of specified factors.",
                 Kmin, K0), call. = FALSE)
  Ks  <- Kmin:Kmax
  pad <- function(K) if (K == K0) Q0 else cbind(Q0, matrix(-1L, J, K - K0))

  done  <- if (!is.null(save_path) && file.exists(save_path)) read.csv(save_path) else NULL
  fits  <- setNames(vector("list", length(Ks)), as.character(Ks))
  rows  <- vector("list", length(Ks))

  for (i in seq_along(Ks)) {
    K <- Ks[i]
    if (!is.null(done) && K %in% done$K) {
      rows[[i]] <- done[done$K == K, , drop = FALSE]
      if (verbose) cat(sprintf("K=%2d  (from checkpoint)\n", K)); next
    }
    t0  <- Sys.time()
    fit <- vbfa(Y, pad(K), max_it = max_it, convChk = FALSE, rseed = rseed,
                orthogonal = orthogonal, ...)
    fs  <- vbfa_fit(fit, Y, pad(K), tau = tau, orthogonal = orthogonal)
    fits[[i]] <- fit
    rows[[i]] <- data.frame(
      K = K, ELBO = as.numeric(fit$ELBO), AIC = fs["AIC"], BIC = fs["BIC"],
      RMSEA = fs["RMSEA"], SRMR = fs["SRMR"], CFI = fs["CFI"], TLI = fs["TLI"],
      t = fs["t"], iter = fit$iter, flag = fit$flag,
      secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), row.names = NULL)
    if (!is.null(save_path))
      write.table(rows[[i]], save_path, sep = ",", row.names = FALSE,
                  col.names = !file.exists(save_path), append = file.exists(save_path))
    if (verbose) cat(sprintf("K=%2d  ELBO=%.1f  BIC=%.0f  CFI=%.3f  conv=%d  %.0fs\n",
                             K, rows[[i]]$ELBO, rows[[i]]$BIC, rows[[i]]$CFI,
                             as.integer(fit$flag == 1), rows[[i]]$secs))
  }
  sweep <- do.call(rbind, rows); sweep <- sweep[order(sweep$K), ]

  ## ---- eligibility: converged within max_it AND passes absolute fit -------
  sweep$converged <- sweep$flag == 1L
  sweep$pass_fit  <- if (!isTRUE(fit_cut)) rep(TRUE, nrow(sweep)) else with(sweep,
      !is.na(RMSEA) & RMSEA <= cutoffs[["RMSEA"]] &
      !is.na(SRMR)  & SRMR  <= cutoffs[["SRMR"]] &
      !is.na(CFI)   & CFI   >= cutoffs[["CFI"]] &
      !is.na(TLI)   & TLI   >= cutoffs[["TLI"]])
  sweep$eligible  <- sweep$converged & sweep$pass_fit

  elig <- sweep[sweep$eligible,  , drop = FALSE]
  conv <- sweep[sweep$converged, , drop = FALSE]
  pick <- function(d, score) if (nrow(d) == 0) NA_integer_ else select_K_elbow(d$K, score, delta, sustain)

  selection <- c(
    elbo_gain         = pick(elig, elig$ELBO),                       # primary
    elbo_raw          = if (nrow(elig)) elig$K[which.max(elig$ELBO)] else NA_integer_,
    bic_gain          = pick(elig, -elig$BIC),                       # reference
    bic_raw           = if (nrow(elig)) elig$K[which.min(elig$BIC)]  else NA_integer_,
    elbo_gain_allconv = pick(conv, conv$ELBO))                       # no fit gate
  K_sel <- unname(selection["elbo_gain"])

  ## ---- attach model + fit for the SELECTED K (refit if from checkpoint) ---
  sel_fit <- sel_stats <- NULL
  if (!is.na(K_sel)) {
    sel_fit <- fits[[as.character(K_sel)]]
    if (is.null(sel_fit)) sel_fit <- vbfa(Y, pad(K_sel), max_it = max_it, convChk = FALSE,
                                          rseed = rseed, orthogonal = orthogonal, ...)
    sel_stats <- vbfa_fit(sel_fit, Y, pad(K_sel), tau = tau, orthogonal = orthogonal)
  }

  if (verbose) {
    n_fail <- sum(!sweep$converged); n_unfit <- sum(sweep$converged & !sweep$pass_fit)
    gate <- if (isTRUE(fit_cut)) "ON" else "OFF"
    cat(sprintf("\n--- selection (delta = %g%%; convergence gate ON, absolute-fit gate %s) ---\n",
                delta, gate))
    cat(sprintf("  excluded: %d non-converged, %d converged-but-failing-fit; %d eligible\n",
                n_fail, n_unfit, nrow(elig)))
    if (is.na(K_sel)) cat("  NO eligible model -- selection is NA (relax fit_cut or extend/converge).\n")
    else cat(sprintf("  ELBO gain  -> K = %d   (selected; model returned)\n", K_sel))
    cat(sprintf("  ELBO raw   -> K = %s\n",  selection["elbo_raw"]))
    cat(sprintf("  BIC  gain  -> K = %s\n",  selection["bic_gain"]))
    cat(sprintf("  BIC  raw   -> K = %s\n",  selection["bic_raw"]))
    cat(sprintf("  ELBO gain (all converged, no fit gate) -> K = %s\n", selection["elbo_gain_allconv"]))
  }
  list(selected_K = K_sel, selection = selection, fit = sel_fit,
       fit_stats = sel_stats, Q_selected = if (!is.na(K_sel)) pad(K_sel) else NULL,
       sweep = sweep, K0 = K0, delta = delta,
       fit_cut = isTRUE(fit_cut), cutoffs = cutoffs)
}
