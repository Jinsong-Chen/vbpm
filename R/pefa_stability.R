## Adjacent-count summaries for pefa() sweeps.
##
## This is deliberately a descriptive extension to PEFA selection.  Columns
## from the smaller solution are considered in fitted order and matched to the
## most congruent unused column in the larger solution.

.tucker <- function(a, b) {
  den <- sqrt(sum(a^2) * sum(b^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(a * b) / den
}

.match_columns <- function(A, B) {
  used <- integer(0)
  idx <- integer(ncol(A))
  phi <- numeric(ncol(A))

  for (i in seq_len(ncol(A))) {
    available <- setdiff(seq_len(ncol(B)), used)
    congruence <- vapply(
      available,
      function(j) abs(.tucker(A[, i], B[, j])),
      numeric(1)
    )
    best <- which.max(congruence)
    pick <- available[best]
    idx[i] <- pick
    phi[i] <- congruence[best]
    used <- c(used, pick)
  }

  inner <- vapply(
    seq_len(ncol(A)),
    function(i) sum(A[, i] * B[, idx[i]]),
    numeric(1)
  )
  list(
    idx = idx,
    phi = phi,
    sign = ifelse(inner < 0, -1, 1),
    unmatched = setdiff(seq_len(ncol(B)), idx)
  )
}

.partition_dominant <- function(L) {
  max.col(abs(L), ties.method = "first")
}

## Adjusted Rand index; return a typed NA when its denominator degenerates.
.ari <- function(p, q) {
  tab <- table(p, q)
  n <- sum(tab)
  if (n < 2L) return(NA_real_)

  pairs <- choose(n, 2)
  same_both <- sum(choose(tab, 2))
  same_p <- sum(choose(rowSums(tab), 2))
  same_q <- sum(choose(colSums(tab), 2))
  expected <- same_p * same_q / pairs
  den <- (same_p + same_q) / 2 - expected
  if (!is.finite(den) || abs(den) <= .Machine$double.eps^0.5)
    return(NA_real_)
  (same_both - expected) / den
}

.na_stability <- function() {
  c(
    phi_min = NA_real_, phi_mean = NA_real_, ari = NA_real_,
    rmsd = NA_real_, max_loading = NA_real_
  )
}

.transition_stability <- function(A, B) {
  valid <- is.matrix(A) && is.matrix(B) && is.numeric(A) && is.numeric(B) &&
    nrow(A) == nrow(B) && ncol(B) == ncol(A) + 1L &&
    all(is.finite(A)) && all(is.finite(B))
  if (!valid) return(.na_stability())

  norms <- c(sqrt(colSums(A^2)), sqrt(colSums(B^2)))
  if (any(norms == 0)) return(.na_stability())

  matched <- .match_columns(A, B)
  B_matched <- B[, matched$idx, drop = FALSE]
  B_matched <- sweep(B_matched, 2, matched$sign, "*")

  c(
    phi_min = min(matched$phi),
    phi_mean = mean(matched$phi),
    ari = .ari(.partition_dominant(A), .partition_dominant(B)),
    rmsd = sqrt(mean((A - B_matched)^2)),
    max_loading = max(abs(B[, matched$unmatched]))
  )
}

.transition_table <- function(n = 0L) {
  data.frame(
    K_from = integer(n),
    K_to = integer(n),
    ELBO_gain = numeric(n),
    BIC_gain = numeric(n),
    phi_min = numeric(n),
    phi_mean = numeric(n),
    ari = numeric(n),
    rmsd = numeric(n),
    max_loading = numeric(n),
    stringsAsFactors = FALSE
  )
}

.build_transitions <- function(sweep, fits, bifactor = FALSE) {
  n <- length(sweep$K) - 1L
  if (n <= 0L) return(.transition_table())

  out <- .transition_table(n)
  out$K_from <- as.integer(sweep$K[-length(sweep$K)])
  out$K_to <- as.integer(sweep$K[-1L])
  out$ELBO_gain <- sweep$ELBO[-1L] - sweep$ELBO[-nrow(sweep)]
  out$BIC_gain <- sweep$BIC[-nrow(sweep)] - sweep$BIC[-1L]

  edge_converged <- vapply(seq_len(n), function(i) {
    isTRUE(sweep$converged[i]) && isTRUE(sweep$converged[i + 1L])
  }, logical(1))
  elbo_finite <- is.finite(sweep$ELBO[-nrow(sweep)]) &
    is.finite(sweep$ELBO[-1L])
  bic_finite <- is.finite(sweep$BIC[-nrow(sweep)]) &
    is.finite(sweep$BIC[-1L])
  out$ELBO_gain[!edge_converged | !elbo_finite] <- NA_real_
  out$BIC_gain[!edge_converged | !bic_finite] <- NA_real_

  for (i in seq_len(n)) {
    if (isTRUE(bifactor) || !edge_converged[i]) {
      stability <- .na_stability()
    } else {
      from <- fits[[as.character(out$K_from[i])]]
      to <- fits[[as.character(out$K_to[i])]]
      stability <- .transition_stability(from$Lam, to$Lam)
    }
    out[i, names(stability)] <- as.list(stability)
  }

  out
}
