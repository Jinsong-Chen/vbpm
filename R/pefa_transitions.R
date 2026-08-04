## Adjacent-count summaries for pefa() sweeps.
##
## Loading/PIP matching is deliberately separate from transition-table
## construction.  The matcher is deterministic and uses base R only; the
## builder adds PEFA-specific convergence gates and bifactor group extraction.

.transition_metric_names <- c(
  "rmsd", "rmsd_max", "phi_min", "ari", "pip_rmsd",
  "collision", "unmatched_max"
)

.transition_numeric_metric_names <- setdiff(
  .transition_metric_names, "collision"
)

.na_transition_metrics <- function() {
  stats::setNames(rep(NA_real_, length(.transition_metric_names)),
                  .transition_metric_names)
}

.match_sweep_loadings <- function(lam1, lam2, K0, eps = 0.1,
                                  pip1 = NULL, pip2 = NULL,
                                  regularized1 = NULL,
                                  regularized2 = NULL) {
  lam1 <- as.matrix(lam1)
  lam2 <- as.matrix(lam2)

  valid_loading <- function(x) {
    is.numeric(x) && !is.complex(x) && all(is.finite(x))
  }
  scalar_whole <- function(x) {
    length(x) == 1L && is.numeric(x) && !is.complex(x) &&
      !is.na(x) && is.finite(x) && x == floor(x)
  }
  orient <- function(x) {
    direction <- ifelse(colMeans(x) < 0, -1, 1)
    sweep(x, 2L, direction, "*")
  }
  free_columns <- function(K, K0) {
    if (K0 < K) seq.int(K0 + 1L, K) else integer(0)
  }
  log_norm2 <- function(x) {
    scale <- max(abs(x))
    if (scale == 0) return(-Inf)
    z <- x / scale
    2 * log(scale) + log(sum(z^2))
  }
  norm2_at_least <- function(x, log_threshold) {
    value <- log_norm2(x)
    if (!is.finite(value)) return(FALSE)
    ## Admit representational equality at the inclusive cutoff without giving
    ## up the log scale's overflow/underflow protection.
    tolerance <- 4 * .Machine$double.eps *
      max(1, abs(value), abs(log_threshold))
    value >= log_threshold - tolerance
  }
  distance_key <- function(a, b) {
    scale <- max(max(abs(a)), max(abs(b)))
    if (scale == 0) return(-Inf)
    z <- a / scale - b / scale
    z_scale <- max(abs(z))
    if (z_scale == 0) return(-Inf)
    z <- z / z_scale
    2 * log(scale) + 2 * log(z_scale) + log(sum(z^2))
  }
  stable_rmsd <- function(a, b) {
    scale <- max(max(abs(a)), max(abs(b)))
    if (scale == 0) return(0)
    z <- a / scale - b / scale
    z_scale <- max(abs(z))
    if (z_scale == 0) return(0)
    z <- z / z_scale
    exp(log(scale) + log(z_scale) + 0.5 * log(mean(z^2)))
  }
  tucker <- function(a, b) {
    scale_a <- max(abs(a))
    scale_b <- max(abs(b))
    if (scale_a == 0 || scale_b == 0) return(NA_real_)
    unit_a <- a / scale_a
    unit_b <- b / scale_b
    unit_a <- unit_a / sqrt(sum(unit_a^2))
    unit_b <- unit_b / sqrt(sum(unit_b^2))
    value <- sum(unit_a * unit_b)
    max(-1, min(1, value))
  }
  adjusted_rand <- function(p, q) {
    tab <- table(p, q)
    n <- sum(tab)
    if (n < 2L) return(NA_real_)

    pairs <- choose(n, 2)
    same_both <- sum(choose(tab, 2))
    same_p <- sum(choose(rowSums(tab), 2))
    same_q <- sum(choose(colSums(tab), 2))
    expected <- same_p * same_q / pairs
    den <- (same_p + same_q) / 2 - expected
    if (!is.finite(den) || abs(den) <= .Machine$double.eps^0.5) {
      return(NA_real_)
    }
    value <- (same_both - expected) / den
    max(-1, min(1, value))
  }

  if (!valid_loading(lam1) || !valid_loading(lam2)) {
    stop("lam1 and lam2 must be real, finite numeric matrices.",
         call. = FALSE)
  }
  if (nrow(lam1) == 0L || nrow(lam1) != nrow(lam2)) {
    stop("lam1 and lam2 must have the same positive number of rows.",
         call. = FALSE)
  }
  if (ncol(lam1) == 0L) {
    stop("lam1 and lam2 must each contain at least one factor column.",
         call. = FALSE)
  }
  if (ncol(lam1) > ncol(lam2)) {
    stop("ncol(lam1) must be no greater than ncol(lam2).", call. = FALSE)
  }
  if (!is.null(rownames(lam1)) && !is.null(rownames(lam2)) &&
      !identical(rownames(lam1), rownames(lam2))) {
    stop("lam1 and lam2 row names must have the same order.", call. = FALSE)
  }
  if (!scalar_whole(K0) || K0 < 0 || K0 > ncol(lam1)) {
    stop("K0 must be one whole number from 0 through ncol(lam1).",
         call. = FALSE)
  }
  K0 <- as.integer(K0)
  if (length(eps) != 1L || !is.numeric(eps) || is.complex(eps) ||
      is.na(eps) || !is.finite(eps) || eps <= 0) {
    stop("eps must be one finite, real, positive number.", call. = FALSE)
  }

  have_pip1 <- !is.null(pip1)
  have_pip2 <- !is.null(pip2)
  if (xor(have_pip1, have_pip2)) {
    stop("Supply both pip1 and pip2, or neither.", call. = FALSE)
  }
  have_mask1 <- !is.null(regularized1)
  have_mask2 <- !is.null(regularized2)
  if (xor(have_mask1, have_mask2)) {
    stop("Supply both regularized1 and regularized2, or neither.",
         call. = FALSE)
  }
  if (have_pip1 && !have_mask1) {
    stop("Supply regularized1 and regularized2 with the PIP matrices.",
         call. = FALSE)
  }
  if (!have_pip1 && have_mask1) {
    stop("Regularized masks are used only with PIP matrices.",
         call. = FALSE)
  }
  if (have_pip1) {
    pip1 <- as.matrix(pip1)
    pip2 <- as.matrix(pip2)
    regularized1 <- as.matrix(regularized1)
    regularized2 <- as.matrix(regularized2)
    compatible_names <- function(x, reference) {
      (is.null(rownames(x)) || is.null(rownames(reference)) ||
         identical(rownames(x), rownames(reference))) &&
        (is.null(colnames(x)) || is.null(colnames(reference)) ||
           identical(colnames(x), colnames(reference)))
    }
    valid_pip <- function(x, reference) {
      is.numeric(x) && !is.complex(x) && all(is.finite(x)) &&
        identical(dim(x), dim(reference)) && all(x >= 0 & x <= 1) &&
        compatible_names(x, reference)
    }
    if (!valid_pip(pip1, lam1) || !valid_pip(pip2, lam2)) {
      stop(paste0(
        "PIP matrices must match their loading dimensions and supplied ",
        "dimnames, and contain only finite values in [0, 1]."
      ), call. = FALSE)
    }
    valid_mask <- function(x, reference) {
      is.logical(x) && !anyNA(x) && identical(dim(x), dim(reference)) &&
        compatible_names(x, reference)
    }
    if (!valid_mask(regularized1, lam1) ||
        !valid_mask(regularized2, lam2)) {
      stop(paste0(
        "Regularized masks must be logical matrices matching their loading ",
        "dimensions and supplied dimnames, without missing values."
      ), call. = FALSE)
    }
  }

  lam1 <- orient(lam1)
  lam2 <- orient(lam2)

  free1 <- free_columns(ncol(lam1), K0)
  free2 <- free_columns(ncol(lam2), K0)
  log_eps <- log(eps)
  keep1 <- free1[vapply(
    free1, function(k) norm2_at_least(lam1[, k], log_eps), logical(1)
  )]
  keep2 <- free2[vapply(
    free2, function(k) norm2_at_least(lam2[, k], log_eps), logical(1)
  )]
  backbone <- seq_len(K0)

  retained1 <- c(backbone, keep1)
  retained2 <- c(backbone, keep2)
  ari <- if (length(retained1) && length(retained2)) {
    part1 <- max.col(abs(lam1[, retained1, drop = FALSE]),
                     ties.method = "first")
    part2 <- max.col(abs(lam2[, retained2, drop = FALSE]),
                     ties.method = "first")
    adjusted_rand(part1, part2)
  } else {
    NA_real_
  }

  out <- c(
    rmsd = NA_real_, rmsd_max = NA_real_, phi_min = NA_real_, ari = ari,
    pip_rmsd = NA_real_,
    collision = NA_real_,
    unmatched_max = NA_real_
  )

  ## ARI needs no factor match, so it remains available if a retained smaller
  ## factor lacks any eligible candidate in the larger fit.
  if (length(keep1) && !length(keep2)) return(out)

  ## Equation 17 matches each retained non-backbone column independently.
  ## Consequently, more than one smaller-fit column may select the same
  ## larger-fit column; the collision flag makes that reuse explicit.
  matched2 <- vapply(keep1, function(k) {
    distance <- vapply(
      keep2,
      function(j) distance_key(lam1[, k], lam2[, j]),
      numeric(1)
    )
    keep2[which.min(distance)]
  }, integer(1))
  out["collision"] <- as.numeric(anyDuplicated(matched2) > 0L)

  idx1 <- c(backbone, keep1)
  idx2 <- c(backbone, matched2)

  ## Screened columns remain eligible for this raw size diagnostic.  It does
  ## not identify a unique "added factor" or impose a salience cutoff.
  unmatched2 <- setdiff(free2, idx2)
  out["unmatched_max"] <- if (length(unmatched2)) {
    max(abs(lam2[, unmatched2, drop = FALSE]))
  } else {
    0
  }

  if (!length(idx1)) return(out)

  common1 <- lam1[, idx1, drop = FALSE]
  common2 <- lam2[, idx2, drop = FALSE]
  phi <- vapply(
    seq_along(idx1),
    function(i) tucker(common1[, i], common2[, i]),
    numeric(1)
  )

  out["rmsd"] <- stable_rmsd(common1, common2)
  column_rmsd <- vapply(
    seq_along(idx1),
    function(i) stable_rmsd(common1[, i], common2[, i]),
    numeric(1)
  )
  out["rmsd_max"] <- max(column_rmsd)
  if (!anyNA(phi)) out["phi_min"] <- min(phi)

  if (have_pip1) {
    eligible <- regularized1[, idx1, drop = FALSE] &
      regularized2[, idx2, drop = FALSE]
    if (any(eligible)) {
      common_pip1 <- pip1[, idx1, drop = FALSE]
      common_pip2 <- pip2[, idx2, drop = FALSE]
      out["pip_rmsd"] <- stable_rmsd(
        common_pip1[eligible], common_pip2[eligible]
      )
    }
  }

  out
}

.transition_table <- function(n = 0L) {
  data.frame(
    K_from = integer(n),
    K_to = integer(n),
    ELBO_gain = numeric(n),
    BIC_gain = numeric(n),
    rmsd = rep(NA_real_, n),
    rmsd_max = rep(NA_real_, n),
    phi_min = rep(NA_real_, n),
    ari = rep(NA_real_, n),
    pip_rmsd = rep(NA_real_, n),
    collision = rep(NA, n),
    unmatched_max = rep(NA_real_, n),
    stringsAsFactors = FALSE
  )
}

.build_transitions <- function(sweep, loadings, pips, Q0,
                               stability_eps = 0.1,
                               bifactor = FALSE) {
  n <- nrow(sweep) - 1L
  if (n <= 0L) return(.transition_table())

  out <- .transition_table(n)
  from <- seq_len(n)
  to <- from + 1L
  out$K_from <- as.integer(sweep$K[from])
  out$K_to <- as.integer(sweep$K[to])

  out$ELBO_gain <- as.numeric(sweep$ELBO[to] - sweep$ELBO[from])
  out$BIC_gain <- as.numeric(sweep$BIC[from] - sweep$BIC[to])
  out$ELBO_gain[!is.finite(out$ELBO_gain)] <- NA_real_
  out$BIC_gain[!is.finite(out$BIC_gain)] <- NA_real_

  edge_converged <- vapply(from, function(i) {
    isTRUE(sweep$converged[i]) && isTRUE(sweep$converged[i + 1L])
  }, logical(1))
  malformed_loading <- logical(n)
  malformed_pip <- logical(n)
  K0 <- ncol(Q0)

  endpoint <- function(x, K) {
    value <- x[[as.character(K)]]
    if (is.null(value)) {
      stop("missing K-named sweep endpoint", call. = FALSE)
    }
    if (isTRUE(bifactor)) {
      expected <- as.integer(K) + 1L
      if (!is.matrix(value) || ncol(value) != expected) {
        stop(
          "bifactor sweep endpoints must be matrices with K + 1 columns",
          call. = FALSE
        )
      }
      ## The first column is the labelled general factor.  Stability is
      ## defined on the group block so the ordinary backbone-aware matcher and
      ## every diagnostic retain exactly the same meaning in both modes.
      value <- value[, -1L, drop = FALSE]
    }
    value
  }
  regularized_mask <- function(lam) {
    if (!is.matrix(Q0) || nrow(Q0) != nrow(lam) || K0 > ncol(lam)) {
      stop("Q0 is incompatible with a sweep endpoint", call. = FALSE)
    }
    added <- ncol(lam) - K0
    mask <- cbind(
      Q0 == -1,
      matrix(TRUE, nrow(lam), added)
    )
    dimnames(mask) <- dimnames(lam)
    mask
  }

  for (i in from) {
    if (!edge_converged[i]) next

    metrics <- tryCatch({
      lam1 <- endpoint(loadings, out$K_from[i])
      lam2 <- endpoint(loadings, out$K_to[i])
      .match_sweep_loadings(
        lam1, lam2, K0 = K0, eps = stability_eps
      )
    }, error = function(e) {
      malformed_loading[i] <<- TRUE
      .na_transition_metrics()
    })

    ## A malformed PIP endpoint should not erase valid loading diagnostics.
    if (!malformed_loading[i]) {
      pip_metric <- tryCatch({
        pip1 <- endpoint(pips, out$K_from[i])
        pip2 <- endpoint(pips, out$K_to[i])
        if (!is.matrix(pip1) || !is.matrix(pip2)) {
          stop("PIP endpoints must be matrices", call. = FALSE)
        }
        with_pip <- .match_sweep_loadings(
          lam1, lam2, K0 = K0, eps = stability_eps,
          pip1 = pip1, pip2 = pip2,
          regularized1 = regularized_mask(lam1),
          regularized2 = regularized_mask(lam2)
        )
        with_pip["pip_rmsd"]
      }, error = function(e) {
        malformed_pip[i] <<- TRUE
        NA_real_
      })
      metrics["pip_rmsd"] <- pip_metric
    }
    out[i, .transition_numeric_metric_names] <-
      as.list(metrics[.transition_numeric_metric_names])
    out$collision[i] <- as.logical(metrics["collision"])
  }

  malformed <- malformed_loading | malformed_pip
  if (any(malformed)) {
    edges <- paste0(out$K_from[malformed], "->", out$K_to[malformed])
    warning(
      paste0(
        "Some stability diagnostics are unavailable for ", sum(malformed),
        " PEFA transition edge", if (sum(malformed) == 1L) "" else "s",
        " (K ", paste(edges, collapse = ", "),
        ") because a converged endpoint has malformed or nonfinite ",
        "loading/PIP data; finite gains and any valid loading diagnostics ",
        "were retained."
      ),
      call. = FALSE
    )
  }

  out
}
