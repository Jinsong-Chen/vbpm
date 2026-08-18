## Partially exploratory factor analysis: a fixed-window candidate sweep.
##
## pefa() fits every K in a consecutive window and reports candidate and
## between-candidate measurements.  It chooses no factor count, threshold,
## horizon, collision policy, or stopping rule; those belong to the analysis
## that consumes the object.
##
## The failure contract is deliberately simple.  A candidate that cannot be
## fitted, or whose required quantities are malformed, aborts the whole call
## with an error naming the K: there is no status/reason taxonomy and no
## partially populated object.  A candidate that merely failed to converge is
## retained, because its returned quantities are still numbers.


# ---- canonical naming -------------------------------------------------

## Item names come from the data when usable, otherwise from a fixed pattern.
## A supplied rownames(Q0) must agree rather than silently winning.
.pefa_item_names <- function(Y, Q0, J) {
  nm <- colnames(Y)
  nm <- if (!is.null(nm) && length(nm) == J && !anyNA(nm) &&
            all(nzchar(nm)) && !anyDuplicated(nm)) {
    as.character(nm)
  } else {
    paste0("item_", seq_len(J))
  }
  qn <- rownames(Q0)
  if (!is.null(qn) && !identical(as.character(qn), nm)) {
    stop("rownames(Q0) must match the item names taken from colnames(Y).",
         call. = FALSE)
  }
  nm
}

.pefa_backbone_names <- function(Q0, K0) {
  nm <- colnames(Q0)
  if (!is.null(nm) && length(nm) == K0 && !anyNA(nm) &&
      all(nzchar(nm)) && !anyDuplicated(nm)) {
    return(as.character(nm))
  }
  if (K0 == 0L) character(0) else paste0("F", seq_len(K0))
}

## Group-block factor names for candidate K: backbone first, then appended
## exploratory columns continuing the F-numbering.  Bifactor prepends "G".
.pefa_dimnames <- function(item_names, backbone_names, K0, K, bifactor) {
  extra <- if (K > K0) paste0("F", seq.int(K0 + 1L, K)) else character(0)
  cols <- c(backbone_names, extra)
  if (isTRUE(bifactor)) cols <- c("G", cols)
  list(item_names, cols)
}


# ---- one candidate ----------------------------------------------------

.pefa_stop_K <- function(K, what) {
  stop(sprintf("PEFA candidate K = %d: %s.", K, what), call. = FALSE)
}

## "Required" means the loading and PIP matrices plus the criterion scalars a
## downstream gain or count consumes.  RMSEA/SRMR/CFI/TLI are descriptive and
## may legitimately be NA -- RMSEA and TLI are undefined at df <= 0 -- so they
## are normalized rather than treated as candidate failures.
.pefa_required_stats <- c("AIC", "BIC", "t")
.pefa_descriptive_stats <- c("RMSEA", "SRMR", "CFI", "TLI")

.pefa_matrix_or_stop <- function(x, dn, kind, K) {
  if (is.null(x)) .pefa_stop_K(K, paste("the", kind, "matrix is missing"))
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x)) {
    .pefa_stop_K(K, paste("the", kind, "matrix is not a real numeric matrix"))
  }
  if (!identical(dim(x), c(length(dn[[1L]]), length(dn[[2L]])))) {
    .pefa_stop_K(K, sprintf("the %s matrix is %d by %d, expected %d by %d",
                            kind, nrow(x), ncol(x),
                            length(dn[[1L]]), length(dn[[2L]])))
  }
  value <- matrix(as.double(x), nrow(x), ncol(x), dimnames = dn)
  if (!all(is.finite(value))) {
    .pefa_stop_K(K, paste("the", kind, "matrix has nonfinite entries"))
  }
  if (identical(kind, "PIP") && (any(value < 0) || any(value > 1))) {
    .pefa_stop_K(K, "the PIP matrix has values outside [0, 1]")
  }
  value
}

.pefa_candidate <- function(fit, fs, K, dn) {
  lam <- .pefa_matrix_or_stop(fit[["Lam"]], dn, "loading", K)
  pip <- .pefa_matrix_or_stop(fit[["pi"]], dn, "PIP", K)

  elbo <- fit[["ELBO"]]
  if (is.null(elbo) || length(elbo) != 1L || !is.numeric(elbo) ||
      is.complex(elbo) || !is.finite(elbo)) {
    .pefa_stop_K(K, "ELBO is not one finite real number")
  }

  need <- c(.pefa_required_stats, .pefa_descriptive_stats)
  if (!is.numeric(fs) || is.complex(fs) || is.null(names(fs)) ||
      !all(need %in% names(fs))) {
    .pefa_stop_K(K, "fit_stats() did not return the required named scalars")
  }
  value <- vapply(need, function(f) as.double(unname(fs[[f]])), numeric(1))
  if (!all(is.finite(value[.pefa_required_stats]))) {
    .pefa_stop_K(K, "AIC, BIC, or t is nonfinite")
  }
  ## A nonfinite descriptive index becomes NA_real_.  That is a reported
  ## value, not a defect; see the Identification section of ?pefa.
  descriptive <- value[.pefa_descriptive_stats]
  descriptive[!is.finite(descriptive)] <- NA_real_

  flag <- fit[["flag"]]
  converged <- length(flag) == 1L && is.numeric(flag) && !is.na(flag) &&
    is.finite(flag) && flag == 1
  iter <- fit[["iter"]]
  iter <- if (length(iter) == 1L && is.numeric(iter) && !is.na(iter) &&
              is.finite(iter)) as.integer(iter) else NA_integer_

  row <- data.frame(
    K = as.integer(K),
    ELBO = as.double(elbo),
    AIC = value[["AIC"]], BIC = value[["BIC"]],
    RMSEA = descriptive[["RMSEA"]], SRMR = descriptive[["SRMR"]],
    CFI = descriptive[["CFI"]], TLI = descriptive[["TLI"]],
    t = value[["t"]],
    iter = iter,
    converged = converged,
    stringsAsFactors = FALSE
  )
  list(row = row, loading = lam, pip = pip)
}


# ---- derived views ----------------------------------------------------

## Percent of the largest positive gain on that path.  Dividing before scaling
## keeps gains near the largest finite double representable.
.pefa_gain_pct <- function(gain) {
  if (!length(gain)) return(numeric(0))
  gmax <- max(gain)
  if (!(gmax > 0)) return(rep(NA_real_, length(gain)))
  100 * (gain / gmax)
}

.pefa_fact_names <- c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                      "unmatched_ssl", "collision")

.pefa_transitions <- function(sweep, pair_facts) {
  n <- nrow(sweep) - 1L
  if (n <= 0L) {
    out <- data.frame(K_from = integer(0), K_to = integer(0),
                      ELBO_gain_pct = numeric(0), BIC_gain_pct = numeric(0),
                      stringsAsFactors = FALSE)
    for (f in .pefa_fact_names) {
      out[[f]] <- if (f == "collision") logical(0) else numeric(0)
    }
    return(out)
  }
  from <- seq_len(n)
  to <- from + 1L

  ## Positive favors K_to on both paths.
  g_elbo <- as.double(sweep$ELBO[to] - sweep$ELBO[from])
  g_bic <- as.double(sweep$BIC[from] - sweep$BIC[to])
  if (!all(is.finite(g_elbo)) || !all(is.finite(g_bic))) {
    stop("PEFA adjacent criterion differences must be finite; a ",
         "nonrepresentable difference is a construction error rather than a ",
         "partially populated path.", call. = FALSE)
  }

  out <- data.frame(
    K_from = as.integer(sweep$K[from]), K_to = as.integer(sweep$K[to]),
    ELBO_gain_pct = .pefa_gain_pct(g_elbo),
    BIC_gain_pct = .pefa_gain_pct(g_bic),
    stringsAsFactors = FALSE
  )
  ## The adjacent rows reuse the same in-memory pair results the all-pairs
  ## view projects; matching never runs a second way for the two outputs.
  key <- paste(out$K_from, out$K_to, sep = "-")
  for (f in .pefa_fact_names) {
    out[[f]] <- vapply(key, function(k) pair_facts[[k]][[f]],
                       if (f == "collision") logical(1) else numeric(1),
                       USE.NAMES = FALSE)
  }
  out
}

## $rmsd stores rmsd_max, the weakest-column value, while $transitions$rmsd is
## the pooled one.  The two names differ on purpose and must not be crossed.
.pefa_persistence <- function(Ks, pair_facts) {
  W <- length(Ks)
  nm <- as.character(Ks)
  blank <- function(fill) matrix(fill, W, W, dimnames = list(nm, nm))
  out <- list(phi = blank(NA_real_), rmsd = blank(NA_real_),
              collision = blank(NA))
  src <- c(phi = "phi_min", rmsd = "rmsd_max", collision = "collision")
  for (i in seq_len(W)) {
    for (j in seq_len(W)) {
      if (i >= j) next
      f <- pair_facts[[paste(Ks[i], Ks[j], sep = "-")]]
      for (m in names(src)) out[[m]][i, j] <- f[[src[[m]]]]
    }
  }
  out
}


# ---- argument validation ----------------------------------------------

.pefa_whole <- function(value, nm, min_value) {
  if (length(value) != 1L || !is.numeric(value) || is.complex(value) ||
      is.na(value) || !is.finite(value) || value < min_value ||
      value != round(value) || value > .Machine$integer.max) {
    stop(nm, " must be a single whole number in [", min_value, ", ",
         .Machine$integer.max, "].", call. = FALSE)
  }
  as.integer(value)
}

.pefa_flag <- function(value, nm) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(nm, " must be TRUE or FALSE.", call. = FALSE)
  }
  as.logical(value)
}

.pefa_bounded <- function(value, nm, lo, hi) {
  if (length(value) != 1L || !is.numeric(value) || is.complex(value) ||
      is.na(value) || !is.finite(value) || value < lo || value > hi) {
    stop(nm, " must be a single finite number in [", lo, ", ", hi, "].",
         call. = FALSE)
  }
  as.double(value)
}


# ---- constructor ------------------------------------------------------

#' Partially exploratory factor analysis over a factor-count window
#'
#' Fits [vbfa()] for every integer `K` in a fixed consecutive window, holding
#' the partially specified backbone `Q0` fixed and appending regularized
#' exploratory columns. It reports candidate summaries and between-candidate
#' loading measurements.
#'
#' `pefa()` **chooses nothing**: no factor count, threshold, persistence
#' horizon, collision policy, or stopping rule. Those belong to the analysis
#' that consumes the object.
#'
#' Ordinary candidates are oblique with diagonal residuals. In bifactor mode
#' `K` counts group factors and [vbfa()] adds one general column, so a stored
#' matrix has `K + 1` columns; every between-candidate comparison removes that
#' column first and works on the `K` group columns.
#'
#' @section The returned object:
#' A list of class `"pefa"` with exactly these components, in order.
#' \describe{
#'   \item{`sweep`}{One row per candidate: `K`, `ELBO`, `AIC`, `BIC`, `RMSEA`,
#'     `SRMR`, `CFI`, `TLI`, `t`, `iter`, `converged`. `AIC` and `BIC` are
#'     hard-selection values at the resolved `tau` and `t` is their effective
#'     parameter count, rank-adjusted when `rank_adjust = TRUE`.}
#'   \item{`transitions`}{One row per adjacent pair: `K_from`, `K_to`,
#'     `ELBO_gain_pct`, `BIC_gain_pct`, then the seven structural facts
#'     `phi_min`, `rmsd`, `rmsd_max`, `ari`, `pip_rmsd`, `unmatched_ssl`,
#'     `collision`. Positive gains favor `K_to`. Each percentage divides by the
#'     largest positive gain on its own path, or is `NA_real_` when no positive
#'     gain exists. Raw criteria stay in `sweep`.}
#'   \item{`persistence`}{A list of three upper-triangular `W` by `W` matrices
#'     over the fitted `K` values: `phi` holds `phi_min`, `rmsd` holds
#'     **`rmsd_max`**, and `collision` is the logical reuse mask. Diagonal and
#'     lower triangle are `NA`. Every cell is a direct endpoint comparison,
#'     never a chain through intervening candidates.}
#'   \item{`loadings`, `pips`}{Named `K`-indexed lists of candidate matrices,
#'     retaining the general column in bifactor mode.}
#'   \item{`Q0`, `settings`}{The canonical backbone, and the resolved fitting
#'     and fit-statistic controls. `settings` holds no analysis rule.}
#' }
#'
#' @section Loading correspondence:
#' For source column `a` and target column `b`, let `s = 1` when `a'b >= 0` and
#' `-1` otherwise, and let `d(a, b) = ||a - s b||^2`. Backbone columns pair by
#' position. Each exploratory source independently selects the exploratory
#' target minimizing `d`, which is Equation 17 of Chen (2023) made invariant to
#' whole-column reflection. Ties break by larger finite absolute Tucker
#' congruence, then by smaller target index.
#'
#' The assignment is independent, not one-to-one, so two sources may select the
#' same target. That reuse is a **collision**, recorded rather than repaired.
#' There is no salience screen: every column participates, and a comparison
#' always contains exactly `K_from` scored pairs.
#'
#' @section Identification:
#' `pefa()` imposes no upper bound on `K` relative to `J`. The Ledermann bound
#' `K <= (2J + 1 - sqrt(8J + 1)) / 2` orients the choice of `Kmax`, but it is
#' derived for an unrestricted factor model and does not constrain a
#' regularized partially confirmatory candidate, so applying it as a package
#' check would be wrong.
#'
#' Beyond that bound expect candidates that do not converge and fit indices
#' that become `NA`: `RMSEA` and `TLI` are undefined at non-positive degrees of
#' freedom. Both are reported evidence about the window rather than failures.
#' `rank_adjust = TRUE` is the sharper diagnostic, since it counts the rank of
#' the Jacobian of the unique covariance elements and therefore tests local
#' identification directly; that count cannot exceed `J(J + 1)/2`, so its
#' degrees of freedom saturate at zero instead of going negative.
#'
#' Because `phi_min` is a minimum over *every* column, an overextracted
#' candidate carrying a near-empty column will report a low `phi_min` and often
#' `collision = TRUE`. Read that with [ssl()]: it usually means the candidate
#' contains a column with nothing in it, not that the retained structure moved.
#'
#' @param Q0 A `J` by `K0` backbone with entries `1`, `0`, or `-1`. Zero
#'   columns (`K0 = 0`, fully exploratory) are allowed.
#' @param Y A numeric `N` by `J` matrix, optionally with `NA`.
#' @param Kmin,Kmax Inclusive window, whole numbers with
#'   `1 <= Kmin <= Kmax` and `Kmin >= ncol(Q0)`.
#' @param bifactor Fit an orthogonal bifactor parameterization.
#' @param general Scalar or length-`J` general design in `{-1, 0, 1}`, used
#'   only when `bifactor = TRUE`.
#' @param v0,max_it,convChk,tolVal Optimizer controls passed to [vbfa()].
#' @param tau Loading-inclusion threshold used by [fit_stats()].
#' @param rank_adjust,rank_max_J Rank-adjustment controls for [fit_stats()].
#' @param verbose Print one progress line per candidate.
#'
#' @return An object of class `"pefa"`; see Details.
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
#' for (k in 1:2) Q0[which(groups == k)[1:2], k] <- 1L
#'
#' fit <- pefa(Q0, sim$dat, 2, 5, verbose = FALSE)
#' fit$sweep
#' fit$transitions[, c("K_from", "K_to", "ELBO_gain_pct", "phi_min")]
#' fit$persistence$phi
#' ssl(fit)
#'
#' ## A count rule is the analysis's, never the package's. This one reads a
#' ## 20% cut on the ELBO path with one candidate of look-ahead.
#' g <- fit$transitions$ELBO_gain_pct
#' fit$transitions$K_from[which(g < 20)[1]]
#' }
#'
#' @seealso [vbfa()], [fit_stats()], [ssl()]
#' @export
pefa <- function(Q0, Y, Kmin, Kmax,
                 bifactor = FALSE, general = 1,
                 v0 = c(.01, .005, .002, .001), max_it = 10000,
                 convChk = FALSE, tolVal = 1e-4, tau = .50,
                 rank_adjust = FALSE, rank_max_J = 100, verbose = TRUE) {
  Y <- as.matrix(Y)
  Q0 <- as.matrix(Q0)
  J <- ncol(Y)
  K0 <- ncol(Q0)
  if (!is.numeric(Y) || is.complex(Y) || nrow(Y) < 1L || J < 1L) {
    stop("Y must be a nonempty numeric N by J matrix.", call. = FALSE)
  }
  if (any(is.nan(Y)) || any(is.infinite(Y))) {
    stop("Y may contain literal NA but not NaN or infinite values.",
         call. = FALSE)
  }
  if (nrow(Q0) != J) stop("nrow(Q0) must equal ncol(Y).", call. = FALSE)
  if (K0 > 0L && (anyNA(Q0) || !all(Q0 %in% c(-1, 0, 1)))) {
    stop("Q0 entries must be -1, 0, or 1 with no missing values.",
         call. = FALSE)
  }

  Kmin <- .pefa_whole(Kmin, "Kmin", 1)
  Kmax <- .pefa_whole(Kmax, "Kmax", 1)
  max_it <- .pefa_whole(max_it, "max_it", 1)
  rank_max_J <- .pefa_whole(rank_max_J, "rank_max_J", 1)
  if (Kmax < Kmin) {
    stop(sprintf("Kmax (%d) < Kmin (%d): the window must be non-decreasing.",
                 Kmax, Kmin), call. = FALSE)
  }
  if (Kmin < K0) {
    stop(sprintf("Kmin (%d) < K0 (%d): the window must reach the backbone.",
                 Kmin, K0), call. = FALSE)
  }
  tau <- .pefa_bounded(tau, "tau", 0, 1)
  tolVal <- .pefa_bounded(tolVal, "tolVal", 0, Inf)
  bifactor <- .pefa_flag(bifactor, "bifactor")
  convChk <- .pefa_flag(convChk, "convChk")
  rank_adjust <- .pefa_flag(rank_adjust, "rank_adjust")
  verbose <- .pefa_flag(verbose, "verbose")
  if (!is.numeric(v0) || is.complex(v0) || !length(v0) || anyNA(v0) ||
      any(!is.finite(v0)) || any(v0 <= 0)) {
    stop("v0 must be a nonempty vector of finite positive numbers.",
         call. = FALSE)
  }
  v0 <- as.double(unname(v0))
  .validate_rank_controls(rank_adjust, rank_max_J)
  if (isTRUE(rank_adjust)) {
    if (J > rank_max_J) {
      stop("rank_adjust = TRUE requires J <= rank_max_J (", J, " > ",
           rank_max_J, ").", call. = FALSE)
    }
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("rank_adjust = TRUE requires the optional numDeriv package.",
           call. = FALSE)
    }
  }

  general_design <- NULL
  if (isTRUE(bifactor)) {
    if (!is.numeric(general) || is.complex(general) || anyNA(general) ||
        !length(general) || !all(general %in% c(-1, 0, 1)) ||
        !(length(general) %in% c(1L, J))) {
      stop("general must be a scalar or length-J vector of -1, 0, or 1.",
           call. = FALSE)
    }
    general_design <- rep_len(as.integer(general), J)
    if (!any(general_design == 1L)) {
      stop("general must specify at least one general loading (1).",
           call. = FALSE)
    }
    if (K0 > 0L && any(apply(Q0 == 1, 2L, all))) {
      stop("With bifactor = TRUE, Q0 supplies only the group backbone.",
           call. = FALSE)
    }
  }

  item_names <- .pefa_item_names(Y, Q0, J)
  backbone_names <- .pefa_backbone_names(Q0, K0)
  Q0 <- matrix(as.integer(Q0), J, K0,
               dimnames = list(item_names, backbone_names))

  Ks <- seq.int(Kmin, Kmax)
  settings <- list(bifactor = bifactor, general = general_design, v0 = v0,
                   max_it = max_it, convChk = convChk, tolVal = tolVal,
                   tau = tau, rank_adjust = rank_adjust,
                   rank_max_J = rank_max_J)

  rows <- vector("list", length(Ks))
  loadings <- vector("list", length(Ks))
  pips <- vector("list", length(Ks))

  for (i in seq_along(Ks)) {
    K <- Ks[i]
    dn <- .pefa_dimnames(item_names, backbone_names, K0, K, bifactor)
    Q <- if (K > K0) cbind(Q0, matrix(-1L, J, K - K0)) else Q0
    args <- list(Y = Y, Q = Q, v0 = v0, max_it = max_it,
                 convChk = convChk, tolVal = tolVal)
    if (isTRUE(bifactor)) {
      args <- c(args, list(bifactor = TRUE, general = general_design))
    }

    ## The per-candidate iteration-limit warning is muffled: nonconvergence is
    ## a recorded fact reported once at the end, not a repeated console event.
    fit <- tryCatch(
      withCallingHandlers(
        do.call(vbfa, args),
        warning = function(w) {
          if (startsWith(conditionMessage(w),
                         "Maximum number of MFVB iterations exceeded")) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) {
        .pefa_stop_K(K, paste("vbfa() failed:", conditionMessage(e)))
      }
    )
    fs <- tryCatch(
      fit_stats(fit, Y = Y, Q = fit$Q, tau = tau,
                rank_adjust = rank_adjust, rank_max_J = rank_max_J),
      error = function(e) {
        .pefa_stop_K(K, paste("fit_stats() failed:", conditionMessage(e)))
      }
    )

    candidate <- .pefa_candidate(fit, fs, K, dn)
    rows[[i]] <- candidate$row
    loadings[[i]] <- candidate$loading
    pips[[i]] <- candidate$pip
    if (isTRUE(verbose)) {
      cat(sprintf("K=%2d  ELBO=%-12s BIC=%-10s CFI=%-7s conv=%s\n", K,
                  format(candidate$row$ELBO, digits = 7),
                  format(candidate$row$BIC, digits = 7),
                  format(candidate$row$CFI, digits = 3),
                  candidate$row$converged))
    }
  }

  sweep <- do.call(rbind, rows)
  rownames(sweep) <- NULL
  names(loadings) <- as.character(Ks)
  names(pips) <- as.character(Ks)

  ## One comparison per ordered pair, computed directly from its endpoints.
  ## Both derived views read these same results.
  offset <- as.integer(isTRUE(bifactor))
  group <- function(m) if (offset == 0L) m else m[, -1L, drop = FALSE]
  reg <- function(K) {
    Q <- if (K > K0) cbind(Q0, matrix(-1L, J, K - K0)) else Q0
    matrix(Q == -1L, J, K)
  }
  pair_facts <- list()
  for (i in seq_along(Ks)) {
    for (j in seq_along(Ks)) {
      if (i >= j) next
      a <- as.character(Ks[i])
      b <- as.character(Ks[j])
      pair_facts[[paste(Ks[i], Ks[j], sep = "-")]] <- .pefa_compare_pair(
        lam_from = group(loadings[[a]]), lam_to = group(loadings[[b]]),
        pip_from = group(pips[[a]]), pip_to = group(pips[[b]]),
        reg_from = reg(Ks[i]), reg_to = reg(Ks[j]), K0 = K0
      )
    }
  }

  out <- list(sweep = sweep,
              transitions = .pefa_transitions(sweep, pair_facts),
              persistence = .pefa_persistence(Ks, pair_facts),
              loadings = loadings, pips = pips, Q0 = Q0,
              settings = settings)
  class(out) <- "pefa"

  bad <- sweep$K[!sweep$converged]
  if (length(bad)) {
    warning("PEFA candidates did not converge at K = ",
            paste(bad, collapse = ", "),
            "; their returned quantities are retained.", call. = FALSE)
  }
  out
}
