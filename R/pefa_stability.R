## Adjacent-count stability machinery for pefa() sweeps (0.8.0).
##
## A transition compares the K-factor solution with the (K+1)-factor solution
## the sweep already fitted. Columns are matched one-to-one by *bottleneck*
## (max-min) Tucker congruence, not greedily in stored order: greedy matching
## is not permutation invariant and need not maximize the binding statistic.

## Documented implementation tolerances.
.pefa_tol <- list(
  loading = 1e-8,   # |loading| at or below this counts as zero
  norm    = 1e-8,   # column norm at or below this makes the column empty
  tie     = 1e-8    # two largest |loadings| within this are an ambiguous tie
)

## Sign-canonical orientation: the largest-absolute loading is made positive,
## with the smallest item index breaking an exact tie. Empty columns are
## returned as exact zeros so sub-tolerance noise cannot steer a match.
.canon_cols <- function(L, tol = .pefa_tol) {
  storage.mode(L) <- "double"
  for (j in seq_len(ncol(L))) {
    v <- L[, j]
    if (sqrt(sum(v^2)) <= tol$norm) { L[, j] <- 0; next }
    i <- which(abs(v) == max(abs(v)))[1L]
    if (v[i] < 0) L[, j] <- -v
  }
  L
}

## Tucker congruence (uncentered cosine); zero for an empty pair.
.tucker <- function(a, b) {
  na <- sqrt(sum(a^2)); nb <- sqrt(sum(b^2))
  if (na <= .pefa_tol$norm || nb <= .pefa_tol$norm) return(0)
  sum(a * b) / (na * nb)
}

## Exact lexicographic max-min assignment of A's columns into B's columns.
## Branch and bound over injections, ordered by descending congruence, pruning
## any branch whose running minimum cannot beat the incumbent. Exact for the
## K sizes a factor-number window can plausibly reach.
.assign_lexmaxmin <- function(P) {
  nA <- nrow(P); nB <- ncol(P)
  best <- list(idx = NULL, key = NULL)
  ## a candidate is better if its ascending-sorted congruence vector is
  ## lexicographically larger; ties fall through to the caller's index rule
  better <- function(key, ref) {
    if (is.null(ref)) return(TRUE)
    d <- key - ref
    nz <- which(abs(d) > 0)
    length(nz) > 0L && d[nz[1L]] > 0
  }
  ord <- lapply(seq_len(nA), function(i) order(P[i, ], decreasing = TRUE))
  rec <- function(i, used, vals) {
    if (i > nA) {
      key <- sort(vals)
      if (better(key, best$key)) best <<- list(idx = used, key = key)
      return(invisible(NULL))
    }
    ## bound: even the best remaining assignment cannot raise the minimum
    if (!is.null(best$key) && length(vals)) {
      if (min(vals) < best$key[1L]) return(invisible(NULL))
    }
    for (j in ord[[i]]) {
      if (j %in% used) next
      rec(i + 1L, c(used, j), c(vals, P[i, j]))
    }
    invisible(NULL)
  }
  rec(1L, integer(0), numeric(0))
  best$idx
}

## Match the smaller solution's columns into the larger's.
## Returns matched larger-column indices (in the caller's stored order),
## per-column congruences, and sign multipliers.
.match_columns <- function(A, B) {
  Ac <- .canon_cols(A); Bc <- .canon_cols(B)
  nA <- ncol(Ac); nB <- ncol(Bc)
  ## canonical column order removes stored-order dependence from tie-breaking
  key_of <- function(M) apply(M, 2, function(v) paste(sprintf("%.12e", v), collapse = "|"))
  oA <- order(key_of(Ac)); oB <- order(key_of(Bc))
  P <- matrix(0, nA, nB)
  for (i in seq_len(nA)) for (j in seq_len(nB))
    P[i, j] <- abs(.tucker(Ac[, oA[i]], Bc[, oB[j]]))
  sel <- .assign_lexmaxmin(P)
  ## map canonical positions back to the caller's stored order
  idx <- integer(nA); phi <- numeric(nA)
  for (i in seq_len(nA)) { idx[oA[i]] <- oB[sel[i]]; phi[oA[i]] <- P[i, sel[i]] }
  ## sign alignment uses the ORIGINAL columns; an exactly zero inner product
  ## is treated as +1 so sign(0) cannot blank a column
  sgn <- vapply(seq_len(nA), function(i) {
    ip <- sum(A[, i] * B[, idx[i]])
    if (sqrt(sum(A[, i]^2)) <= .pefa_tol$norm ||
        sqrt(sum(B[, idx[i]]^2)) <= .pefa_tol$norm) NA_real_
    else if (ip < 0) -1 else 1
  }, numeric(1))
  list(idx = idx, phi = phi, sign = sgn,
       surplus = setdiff(seq_len(nB), idx),
       empty_from = sum(apply(A, 2, function(v) sqrt(sum(v^2)) <= .pefa_tol$norm)),
       empty_to   = sum(apply(B, 2, function(v) sqrt(sum(v^2)) <= .pefa_tol$norm)))
}

## Dominant-factor partition. Items with every |loading| at or below tolerance,
## and items whose two largest |loadings| tie, become distinct singleton
## classes so unrelated fallback items cannot agree with one another.
.partition_dominant <- function(L, tol = .pefa_tol) {
  J <- nrow(L)
  lab <- integer(J); n_un <- 0L; n_amb <- 0L
  nxt <- ncol(L) + 1L
  for (i in seq_len(J)) {
    a <- abs(L[i, ])
    if (max(a) <= tol$loading) { lab[i] <- nxt; nxt <- nxt + 1L; n_un <- n_un + 1L; next }
    s <- sort(a, decreasing = TRUE)
    if (length(s) > 1L && (s[1L] - s[2L]) <= tol$tie) {
      lab[i] <- nxt; nxt <- nxt + 1L; n_amb <- n_amb + 1L; next
    }
    lab[i] <- which.max(a)
  }
  list(labels = lab, n_unassigned = n_un, n_ambiguous = n_amb)
}

## Adjusted Rand index; NA (never NaN) when the denominator degenerates.
.ari <- function(p, q) {
  tab <- table(p, q); n <- sum(tab)
  if (n < 2L) return(NA_real_)
  s1 <- sum(choose(tab, 2))
  s2 <- sum(choose(rowSums(tab), 2))
  s3 <- sum(choose(colSums(tab), 2))
  expct <- s2 * s3 / choose(n, 2)
  den <- (s2 + s3) / 2 - expct
  if (abs(den) <= .Machine$double.eps^0.5) return(NA_real_)
  (s1 - expct) / den
}

## One transition row's stability block. `A`/`B` are the K and K+1 loading
## matrices; gating (unsupported mode, non-converged endpoint) is applied by
## the caller, which passes status codes in `gate`.
.transition_stability <- function(A, B, gate = character(0)) {
  na_block <- list(phi_min = NA_real_, phi_mean = NA_real_, n_phi90 = NA_integer_,
                   ari = NA_real_, rmsd = NA_real_, surplus_max = NA_real_,
                   n_unassigned_from = NA_integer_, n_unassigned_to = NA_integer_,
                   n_ambiguous_from = NA_integer_, n_ambiguous_to = NA_integer_,
                   status = gate, match = NULL)
  if (length(gate)) return(na_block)

  m  <- .match_columns(A, B)
  Bm <- B[, m$idx, drop = FALSE]
  sg <- ifelse(is.na(m$sign), 1, m$sign)
  Bm <- sweep(Bm, 2, sg, "*")
  pa <- .partition_dominant(A); pb <- .partition_dominant(B)
  ari <- .ari(pa$labels, pb$labels)

  status <- character(0)
  if (m$empty_from > 0 || m$empty_to > 0)
    status <- c(status, sprintf("empty column: from = %d, to = %d",
                                m$empty_from, m$empty_to))
  if (is.na(ari)) status <- c(status, "ARI undefined: zero denominator")

  list(phi_min = min(m$phi), phi_mean = mean(m$phi),
       n_phi90 = sum(m$phi >= .90),
       ari = ari,
       rmsd = sqrt(mean((A - Bm)^2)),
       surplus_max = if (length(m$surplus)) max(abs(B[, m$surplus])) else NA_real_,
       n_unassigned_from = pa$n_unassigned, n_unassigned_to = pb$n_unassigned,
       n_ambiguous_from  = pa$n_ambiguous,  n_ambiguous_to  = pb$n_ambiguous,
       status = status, match = m)
}

## Assemble the whole $transitions table for a sweep.
.build_transitions <- function(sweep, fits, unsupported) {
  Ks <- sweep$K
  n  <- length(Ks) - 1L
  cols <- function(n) data.frame(
    K_from = integer(n), K_to = integer(n),
    Objective_gain = numeric(n), ELBO_gain = numeric(n), BIC_gain = numeric(n),
    phi_min = numeric(n), phi_mean = numeric(n), n_phi90 = integer(n),
    ari = numeric(n), rmsd = numeric(n), surplus_max = numeric(n),
    n_unassigned_from = integer(n), n_unassigned_to = integer(n),
    n_ambiguous_from = integer(n), n_ambiguous_to = integer(n),
    stability_status = character(n), stringsAsFactors = FALSE)
  if (n <= 0L) return(cols(0L))
  out <- cols(n)
  details <- vector("list", n)
  for (i in seq_len(n)) {
    kf <- Ks[i]; kt <- Ks[i + 1L]
    out$K_from[i] <- kf; out$K_to[i] <- kt
    out$Objective_gain[i] <- sweep$Objective[i + 1L] - sweep$Objective[i]
    out$ELBO_gain[i]      <- sweep$ELBO[i + 1L]      - sweep$ELBO[i]
    out$BIC_gain[i]       <- sweep$BIC[i]            - sweep$BIC[i + 1L]

    gate <- character(0)
    if (unsupported) gate <- c(gate, "unsupported: rotational indeterminacy")
    bad <- Ks[c(i, i + 1L)][!sweep$converged[c(i, i + 1L)]]
    if (length(bad))
      gate <- c(gate, sprintf("endpoint not converged: K = %s",
                              paste(bad, collapse = ", ")))
    A <- fits[[as.character(kf)]]$Lam; B <- fits[[as.character(kt)]]$Lam
    if (is.null(A) || is.null(B))
      gate <- c(gate, "fit not retained")
    st <- .transition_stability(A, B, gate)
    for (nm in c("phi_min", "phi_mean", "n_phi90", "ari", "rmsd", "surplus_max",
                 "n_unassigned_from", "n_unassigned_to",
                 "n_ambiguous_from", "n_ambiguous_to"))
      out[[nm]][i] <- st[[nm]]
    out$stability_status[i] <- if (length(st$status))
      paste(st$status, collapse = "; ") else "ok"
    details[[i]] <- st$match
  }
  attr(out, "detail") <- details
  out
}
