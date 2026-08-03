## Partially exploratory factor analysis: an in-memory window sweep with
## ELBO/BIC factor-number selection.

#' Gain rule for choosing the number of factors
#'
#' Selects a factor count by comparing each marginal gain with a percentage of
#' the largest gain on the same fixed criterion path. The rule is invariant to
#' positive affine transformations of that path; candidates must otherwise use
#' the same data, priors, tuning, and optimization setup.
#'
#' A candidate is selected only when all `sustain` gains following it are
#' strictly below the threshold. The complete look-ahead is required: a short
#' suffix at the right edge does not qualify.
#'
#' @param K Consecutive candidate factor counts (any order).
#' @param score A finite criterion for which larger is better. Pass the negative
#'   of an information criterion for which smaller is better.
#' @param delta Threshold as a percentage of the largest marginal gain.
#' @param sustain Number of consecutive sub-threshold gains required.
#'
#' @return One factor count from `K`.
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' @examples
#' elbo <- c(-5400, -5100, -4900, -4895, -4892)
#' select_K_elbow(2:6, elbo)
#'
#' bic <- c(10900, 10500, 10310, 10395, 10480)
#' select_K_elbow(2:6, -bic)
#'
#' jagged <- c(-5400, -5100, -5095, -4900, -4897, -4895)
#' select_K_elbow(2:7, jagged, sustain = 1)
#' select_K_elbow(2:7, jagged, sustain = 2)
#'
#' @seealso [pefa()]
#' @export
select_K_elbow <- function(K, score, delta = 10, sustain = 2) {
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
    stop("score must be finite and numeric.", call. = FALSE)
  if (length(delta) != 1L || !is.numeric(delta) || !is.finite(delta) ||
      delta < 0 || delta > 100)
    stop("delta must be a single number in [0, 100].", call. = FALSE)
  sustain <- .pefa_sustain(sustain)

  o <- order(K)
  K <- as.integer(K[o])
  score <- score[o]
  if (length(K) > 1L && any(diff(K) != 1L))
    stop("K must contain consecutive whole numbers.", call. = FALSE)
  .scan_gain(K, diff(score), delta, sustain)$selected_K
}

#' Partially exploratory factor analysis with gain-rule selection
#'
#' Fits [vbfa()] over a consecutive factor-count window while holding the
#' partially specified backbone `Q0` fixed and adding exploratory columns. It
#' reports raw and gain-rule selections for ELBO and BIC. AIC and absolute fit
#' indices remain descriptive trajectories and do not gate selection.
#'
#' Ordinary `pefa()` candidates are oblique and have diagonal residuals. The
#' optional bifactor sweep is retained as a research extension: its reported
#' `K` is the number of group factors and [vbfa()] adds one general factor.
#' Bifactor sweeps are not presented as a reliable way to discover the number
#' of group factors.
#'
#' @param Q0 A `J` by `K0` backbone matrix with entries `1`, `0`, or `-1`.
#' @param Y An `N` by `J` data matrix.
#' @param Kmin,Kmax Inclusive factor-count window, with `Kmin >= ncol(Q0)`.
#' @param cuts A named numeric vector of thresholds in `[0, 100]`. Names must
#'   be unique after trimming and exactly one must be `primary`.
#' @param sustain Number of consecutive sub-threshold gains required.
#' @param bifactor Logical; fit the labelled bifactor extension. In this mode
#'   `Kmin` and `Kmax` count group factors.
#' @param general General-column design passed to [vbfa()] in bifactor mode.
#' @param max_it Maximum iterations per [vbfa()] path stage.
#' @param tau Hard-selection threshold passed to [fit_stats()].
#' @param verbose Logical; print per-candidate progress and the primary result.
#' @param ... Named estimation controls passed to [vbfa()]. PEFA controls its
#'   own model class: `ld`, `Qe`, `ld_control`, and `orthogonal` are rejected.
#'
#' @return An object of class `pefa` with these three exact data tables:
#'   \describe{
#'     \item{`$sweep`}{Candidate properties: `K`, `ELBO`, `AIC`, `BIC`,
#'       `RMSEA`, `SRMR`, `CFI`, `TLI`, `t`, `iter`, `secs`, `converged`.}
#'     \item{`$transitions`}{Adjacent-edge gains and descriptive stability:
#'       `K_from`, `K_to`, `ELBO_gain`, `BIC_gain`, `phi_min`, `phi_mean`,
#'       `ari`, `rmsd`, `max_loading`. Stability is unavailable for bifactor
#'       endpoints, failed endpoints, and edges where either solution
#'       contains a zero-norm loading column.}
#'     \item{`$selection`}{Rule results: `criterion`, `form`, `cut_name`,
#'       `cut`, `gain_max`, `threshold`, `selected_K`, `boundary`.}
#'   }
#'   Candidate fits are always retained in `$fits`. The remaining fields are
#'   `call`, `selected_K`, `window`, `cuts`, `sustain`, `bifactor`, `general`,
#'   and `tau`.
#'
#' @section Complete paths and boundaries:
#' Selection requires every candidate on a criterion path to have converged
#' and to have a finite value. A failed path produces `NA` raw and gain
#' selections; candidates on opposite sides of a gap are never joined. A
#' lower or upper boundary result calls for a complete rerun over a wider
#' window. If the lower boundary equals the backbone size, reduce or revise
#' the backbone instead of requesting fewer factors.
#'
#' @section Comparability:
#' Positive affine transformations of one fixed criterion path leave the gain
#' rule unchanged. This does not make sweeps with different preprocessing,
#' priors, hyperparameters, `tau`, or optimization settings comparable.
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' Chen, J. (2023). Fully and partially exploratory factor analysis with
#' bi-level Bayesian regularization. *Behavior Research Methods*, 55(4),
#' 2125--2142. \doi{10.3758/s13428-022-01884-7}
#'
#' @examples
#' \donttest{
#' sim <- sim_fa(N = 400, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' Q0 <- matrix(-1L, ncol(sim$dat), 2)
#' groups <- rep(1:3, each = 6)
#' for (k in 1:2) {
#'   anchors <- which(groups == k)[1:2]
#'   Q0[anchors, ] <- 0L
#'   Q0[anchors, k] <- 1L
#' }
#' result <- pefa(Q0, sim$dat, 2, 4, verbose = FALSE)
#' result$sweep
#' result$selection
#' selected_fit(result)
#' }
#'
#' @seealso [vbfa()], [fit_stats()], [select_K_elbow()]
#' @export
pefa <- function(Q0, Y, Kmin, Kmax,
                 cuts = c(primary = 10), sustain = 2,
                 bifactor = FALSE, general = 1,
                 max_it = 10000, tau = 0.50, verbose = TRUE, ...) {
  call <- match.call()
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names) || anyNA(dot_names) ||
        any(!nzchar(trimws(dot_names))))
      stop("Every argument in ... must be named.", call. = FALSE)
    if (anyDuplicated(dot_names))
      stop("Arguments in ... must have unique names.", call. = FALSE)

    forbidden <- c("ld", "Qe", "ld_control", "orthogonal")
    changes_mode <- vapply(dot_names, function(nm)
      any(startsWith(forbidden, nm)), logical(1))
    if (any(changes_mode))
      stop("pefa() does not allow mode controls in ...: ",
           paste(dot_names[changes_mode], collapse = ", "),
           ". Use vbfa() directly for orthogonal or local-dependence fits.",
           call. = FALSE)

    controlled <- c("Y", "Q", "bifactor", "general", "max_it")
    collisions <- vapply(dot_names, function(nm)
      any(startsWith(controlled, nm)), logical(1))
    if (any(collisions))
      stop("pefa() controls these vbfa() arguments: ",
           paste(dot_names[collisions], collapse = ", "), ".", call. = FALSE)
  }

  cuts <- .pefa_cuts(cuts)
  sustain <- .pefa_sustain(sustain)
  Y <- as.matrix(Y)
  Q0 <- as.matrix(Q0)
  J <- ncol(Y)
  K0 <- ncol(Q0)

  if (nrow(Q0) != J)
    stop("nrow(Q0) must equal ncol(Y).", call. = FALSE)
  if (!all(Q0 %in% c(-1, 0, 1)))
    stop("Q0 entries must be in {-1, 0, 1}.", call. = FALSE)
  for (nm in c("Kmin", "Kmax")) {
    value <- get(nm)
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
        value < 1 || value != round(value))
      stop(nm, " must be a single positive whole number.", call. = FALSE)
  }
  Kmin <- as.integer(Kmin)
  Kmax <- as.integer(Kmax)
  if (Kmax < Kmin)
    stop(sprintf("Kmax (%d) < Kmin (%d): the window must be non-decreasing.",
                 Kmax, Kmin), call. = FALSE)
  if (Kmin < K0)
    stop(sprintf("Kmin (%d) < K0 (%d): the window must be at least the backbone size.",
                 Kmin, K0), call. = FALSE)
  if (length(max_it) != 1L || !is.numeric(max_it) || !is.finite(max_it) ||
      max_it < 1 || max_it != round(max_it))
    stop("max_it must be a single positive whole number.", call. = FALSE)
  if (length(tau) != 1L || !is.numeric(tau) || !is.finite(tau) ||
      tau < 0 || tau > 1)
    stop("tau must be a single number in [0, 1].", call. = FALSE)
  if (length(bifactor) != 1L || !is.logical(bifactor) || is.na(bifactor))
    stop("bifactor must be TRUE or FALSE.", call. = FALSE)
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)

  if (bifactor && K0 > 0L && any(apply(Q0 == 1, 2L, all)))
    stop("With bifactor = TRUE, Q0 must contain only the group backbone; ",
         "the general column is added internally.", call. = FALSE)

  Ks <- seq.int(Kmin, Kmax)
  pad <- function(K) {
    if (K == K0) Q0 else cbind(Q0, matrix(-1L, J, K - K0))
  }
  fits <- setNames(vector("list", length(Ks)), as.character(Ks))
  rows <- vector("list", length(Ks))

  for (i in seq_along(Ks)) {
    K <- Ks[i]
    started <- proc.time()[["elapsed"]]
    args <- list(Y = Y, Q = pad(K), max_it = as.integer(max_it))
    if (bifactor)
      args <- c(args, list(bifactor = TRUE, general = general))
    ## A failed candidate is reported once for the whole criterion path below;
    ## suppress only vbfa()'s per-candidate iteration-limit warning.
    fit <- withCallingHandlers(
      do.call(vbfa, c(args, dots)),
      warning = function(w) {
        if (startsWith(conditionMessage(w),
                       "Maximum number of MFVB iterations exceeded"))
          invokeRestart("muffleWarning")
      }
    )
    fs <- fit_stats(fit, Y = Y, Q = fit$Q, tau = tau,
                    rank_adjust = FALSE)
    elapsed <- max(0, proc.time()[["elapsed"]] - started)
    converged <- isTRUE(fit$converged)

    fits[[i]] <- fit
    rows[[i]] <- data.frame(
      K = as.integer(K),
      ELBO = as.numeric(fit$ELBO),
      AIC = unname(fs[["AIC"]]),
      BIC = unname(fs[["BIC"]]),
      RMSEA = unname(fs[["RMSEA"]]),
      SRMR = unname(fs[["SRMR"]]),
      CFI = unname(fs[["CFI"]]),
      TLI = unname(fs[["TLI"]]),
      t = unname(fs[["t"]]),
      iter = as.integer(fit$iter),
      secs = as.numeric(elapsed),
      converged = converged,
      stringsAsFactors = FALSE
    )
    if (verbose)
      cat(sprintf("K=%2d  ELBO=%.1f  BIC=%.0f  CFI=%.3f  conv=%d  %.1fs\n",
                  K, rows[[i]]$ELBO, rows[[i]]$BIC, rows[[i]]$CFI,
                  as.integer(converged), elapsed))
  }

  sweep <- do.call(rbind, rows)
  rownames(sweep) <- NULL
  transitions <- .build_transitions(sweep, fits, bifactor = bifactor)
  selection <- .pefa_selection(sweep, transitions, cuts, sustain)
  primary <- which(selection$criterion == "elbo" &
                   selection$form == "gain" &
                   selection$cut_name == "primary")[1L]
  selected_K <- selection$selected_K[primary]

  if (verbose) {
    shown <- if (is.na(selected_K)) "NA" else as.character(selected_K)
    cat(sprintf("PEFA primary ELBO gain selection: K = %s\n", shown))
  }

  out <- list(
    call = call,
    selected_K = selected_K,
    fits = fits,
    sweep = sweep,
    transitions = transitions,
    selection = selection,
    window = list(K = Ks, Kmin = Kmin, Kmax = Kmax, K0 = K0),
    cuts = cuts,
    sustain = sustain,
    bifactor = bifactor,
    general = general,
    tau = tau
  )
  class(out) <- "pefa"
  out
}

#' Methods for a PEFA factor-number sweep
#'
#' `selected_fit()` indexes the always-retained candidate fits using the
#' primary ELBO gain selection. `summary()` and the print methods read the
#' sweep, transition, and selection tables directly.
#'
#' @param object,x A result returned by [pefa()].
#' @param ... Further arguments. For `plot()`, graphical arguments passed to
#'   the first base plot in the selected display.
#' @param type Plot type: criterion trajectories (`"objective"`), one gain
#'   trajectory (`"gain"`), or descriptive absolute-fit trajectories
#'   (`"fit"`).
#' @param criterion Criterion used when `type = "gain"`: `"elbo"` (default)
#'   or `"bic"`.
#' @param digits Number of decimals used when printing fit and stability
#'   summaries.
#'
#' @return `selected_fit()` returns the selected [vbfa()] fit, or `NULL` when
#'   selection is unavailable. `summary()` returns a compact `summary.pefa`
#'   object without fit payloads. Print and plot methods return their input
#'   invisibly.
#'
#' @name pefa-methods
NULL

#' @rdname pefa-methods
#' @export
selected_fit <- function(object, ...) UseMethod("selected_fit")

#' @rdname pefa-methods
#' @export
selected_fit.pefa <- function(object, ...) {
  if (length(object$selected_K) != 1L || is.na(object$selected_K)) return(NULL)
  object$fits[[as.character(object$selected_K)]]
}

.print_pefa_rules <- function(selection, indent = "  ") {
  for (i in seq_len(nrow(selection))) {
    row <- selection[i, , drop = FALSE]
    label <- if (row$form == "raw") {
      paste(toupper(row$criterion), "raw")
    } else {
      sprintf("%s gain [%s = %g%%]", toupper(row$criterion),
              row$cut_name, row$cut)
    }
    selected <- if (is.na(row$selected_K)) "NA" else as.character(row$selected_K)
    cat(sprintf("%s%-29s -> K = %-3s (%s)\n", indent, label, selected,
                row$boundary))
  }
  invisible(NULL)
}

.pefa_boundary_advice <- function(x) {
  i <- which(x$selection$criterion == "elbo" &
             x$selection$form == "gain" &
             x$selection$cut_name == "primary")[1L]
  boundary <- x$selection$boundary[i]
  if (identical(boundary, "upper"))
    cat("  Upper-window selection: rerun the complete sweep with a larger Kmax.\n")
  if (identical(boundary, "lower")) {
    if (x$window$Kmin > x$window$K0)
      cat("  Lower-window selection: rerun the complete sweep with a smaller Kmin.\n")
    else
      cat("  Lower selection equals the backbone size; reduce or revise Q0 to search lower.\n")
  }
  invisible(NULL)
}

#' @rdname pefa-methods
#' @export
print.pefa <- function(x, ...) {
  mode <- if (isTRUE(x$bifactor)) " group factors (+ 1 general)" else ""
  cat(sprintf("PEFA sweep: K = %d:%d%s; %d of %d converged\n",
              x$window$Kmin, x$window$Kmax, mode,
              sum(x$sweep$converged), nrow(x$sweep)))
  .print_pefa_rules(x$selection)
  .pefa_boundary_advice(x)
  invisible(x)
}

#' @rdname pefa-methods
#' @export
summary.pefa <- function(object, ...) {
  out <- list(
    call = object$call,
    bifactor = object$bifactor,
    general = object$general,
    window = object$window,
    cuts = object$cuts,
    sustain = object$sustain,
    tau = object$tau,
    selected_K = object$selected_K,
    sweep = object$sweep,
    transitions = object$transitions,
    selection = object$selection,
    convergence = c(converged = sum(object$sweep$converged),
                    total = nrow(object$sweep))
  )
  class(out) <- "summary.pefa"
  out
}

#' @rdname pefa-methods
#' @export
print.summary.pefa <- function(x, digits = 3, ...) {
  mode <- if (isTRUE(x$bifactor)) " group factors (+ 1 general)" else ""
  cat(sprintf("PEFA sweep summary: K = %d:%d%s; sustain = %d\n",
              x$window$Kmin, x$window$Kmax, mode, x$sustain))
  .print_pefa_rules(x$selection)

  cat("\nCandidates:\n")
  candidates <- x$sweep
  if (isTRUE(x$bifactor))
    candidates <- cbind(K = candidates$K,
                        K_total = candidates$K + 1L,
                        candidates[, setdiff(names(candidates), "K"), drop = FALSE])
  for (nm in c("ELBO", "AIC", "BIC"))
    candidates[[nm]] <- round(candidates[[nm]], 1L)
  for (nm in c("RMSEA", "SRMR", "CFI", "TLI"))
    candidates[[nm]] <- round(candidates[[nm]], digits)
  candidates$selected <- ifelse(!is.na(x$selected_K) &
                                  candidates$K == x$selected_K, "*", "")
  print(candidates, row.names = FALSE)

  if (nrow(x$transitions)) {
    cat("\nTransitions:\n")
    transitions <- x$transitions
    numeric_cols <- vapply(transitions, is.numeric, logical(1))
    transitions[numeric_cols] <- lapply(transitions[numeric_cols], round, digits)
    print(transitions, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname pefa-methods
#' @importFrom graphics abline axis legend lines mtext par plot
#' @export
plot.pefa <- function(x, type = c("objective", "gain", "fit"),
                      criterion = c("elbo", "bic"), ...) {
  type <- match.arg(type)
  criterion <- match.arg(criterion)
  d <- x$sweep
  K <- d$K
  selected <- x$selected_K
  col1 <- "#1f77b4"
  col2 <- "#d62728"

  padded_range <- function(values, top = 0.20) {
    values <- values[is.finite(values)]
    if (!length(values)) return(c(0, 1))
    r <- range(values)
    span <- diff(r)
    if (span == 0) span <- abs(r[1L]) * 0.01 + 1
    c(r[1L] - 0.05 * span, r[2L] + top * span)
  }
  mark_selected <- function(value = selected) {
    if (length(value) == 1L && !is.na(value))
      abline(v = value, lty = 3, col = "grey40")
  }

  if (type == "objective") {
    op <- par(mar = c(5, 4.2, 4, 4.2) + 0.1)
    on.exit(par(op))
    plot(K, d$ELBO, type = "b", pch = 16, col = col1,
         ylim = padded_range(d$ELBO), xaxt = "n",
         xlab = "Number of factors (K)",
         ylab = "ELBO (higher is better)", main = "PEFA criteria", ...)
    axis(1, at = K)
    mark_selected()
    par(new = TRUE)
    plot(K, d$BIC, type = "b", pch = 17, lty = 2, col = col2,
         ylim = padded_range(d$BIC), axes = FALSE, xlab = "", ylab = "")
    axis(4, col = col2, col.axis = col2)
    mtext("BIC (lower is better)", side = 4, line = 2.8, col = col2)
    legend("top", c("ELBO (left axis)", "BIC (right axis)"),
           col = c(col1, col2), lty = c(1, 2), pch = c(16, 17),
           bty = "n", cex = 0.9)
  } else if (type == "gain") {
    transition_col <- if (criterion == "elbo") "ELBO_gain" else "BIC_gain"
    gains <- x$transitions[[transition_col]]
    rules <- x$selection[x$selection$criterion == criterion &
                           x$selection$form == "gain", , drop = FALSE]
    primary_selected <- rules$selected_K[rules$cut_name == "primary"][1L]
    finite_threshold <- is.finite(rules$threshold)
    y_values <- c(gains, rules$threshold[finite_threshold], 0)
    if (nrow(x$transitions)) {
      plot(x$transitions$K_to, gains, type = "b", pch = 16, col = col1,
           xlim = range(K), ylim = padded_range(y_values), xaxt = "n",
           xlab = "Number of factors (K)", ylab = "Marginal gain",
           main = paste(toupper(criterion), "gain sweep"), ...)
    } else {
      plot(K, rep(NA_real_, length(K)), type = "n",
           xlim = range(K), ylim = padded_range(y_values), xaxt = "n",
           xlab = "Number of factors (K)", ylab = "Marginal gain",
           main = paste(toupper(criterion), "gain sweep"), ...)
    }
    axis(1, at = K)
    abline(h = 0, col = "grey85")
    if (any(finite_threshold))
      for (value in rules$threshold[finite_threshold])
        abline(h = value, lty = 3, col = "grey40")
    mark_selected(primary_selected)
    labels <- toupper(paste(criterion, "gain"))
    lty <- 1
    pch <- 16
    cols <- col1
    if (any(finite_threshold)) {
      labels <- c(labels,
                  sprintf("%s = %g", rules$cut_name[finite_threshold],
                          rules$threshold[finite_threshold]))
      lty <- c(lty, rep(3, sum(finite_threshold)))
      pch <- c(pch, rep(NA, sum(finite_threshold)))
      cols <- c(cols, rep("grey40", sum(finite_threshold)))
    }
    legend("top", labels, col = cols, lty = lty, pch = pch,
           bty = "n", cex = 0.9)
  } else {
    op <- par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1) + 0.1,
              oma = c(0, 0, 2, 0))
    on.exit(par(op))
    plot(K, d$RMSEA, type = "b", pch = 16, col = col1,
         ylim = padded_range(c(d$RMSEA, d$SRMR)), xaxt = "n",
         xlab = "Number of factors (K)", ylab = "Index (lower is better)",
         main = "Absolute misfit", ...)
    axis(1, at = K)
    lines(K, d$SRMR, type = "b", pch = 17, lty = 2, col = col2)
    mark_selected()
    legend("top", c("RMSEA", "SRMR"), col = c(col1, col2),
           lty = c(1, 2), pch = c(16, 17), bty = "n", cex = 0.9)

    plot(K, d$CFI, type = "b", pch = 16, col = col1,
         ylim = padded_range(c(d$CFI, d$TLI)), xaxt = "n",
         xlab = "Number of factors (K)", ylab = "Index (higher is better)",
         main = "Incremental fit")
    axis(1, at = K)
    lines(K, d$TLI, type = "b", pch = 17, lty = 2, col = col2)
    mark_selected()
    legend("bottom", c("CFI", "TLI"), col = c(col1, col2),
           lty = c(1, 2), pch = c(16, 17), bty = "n", cex = 0.9)
    mtext("PEFA descriptive fit trajectories", outer = TRUE, font = 2)
  }
  invisible(x)
}
