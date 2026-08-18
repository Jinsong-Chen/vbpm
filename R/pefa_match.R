## Structural matching for one ordered candidate pair (K_from < K_to).
##
## Section 2 of the 0.9.0 contract: pairwise sign alignment, an independent
## max-|phi| assignment with a two-level tie rule, the scored pair set
## P_kl = position-matched backbone union assigned exploratory, and the
## matched-pair metrics over it.  Every returned field follows its own
## dependency chain, so a missing input gives that one field its typed NA and
## blanks nothing else.  Only a malformed loading endpoint makes the row
## unavailable.  Endpoint convergence is deliberately not an input anywhere in
## this file, because suppressing a definable number would be exactly the
## endpoint-quality verdict that sections 0.1.2 and 2.10 reject.
##
## The caller supplies group-block matrices: in bifactor mode the leading
## general column has already been stripped, so every column index computed
## and reported below is a group-block index.

## Column sums of squares, written as an explicit per-column `sum(col^2)`.
## The eligibility screen `sum(col^2) >= stability_eps` is exact by contract,
## so nothing -- no rescaling, no rounding, no tolerance -- may come between
## the square and the comparison.
.pefa_ssl <- function(x) {
  vapply(seq_len(ncol(x)), function(k) sum(x[, k]^2), numeric(1))
}

## s(a, b) = 1 when a'b >= 0, else -1; 1 at a zero dot product, which is what
## keeps a zero-norm backbone column's aligned RMSD defined.  The dot product
## is taken on rescaled copies so that a large-magnitude endpoint cannot
## overflow before its sign is read; dividing by a positive scale leaves the
## sign untouched, and negating a whole column negates the sum exactly, which
## is what makes every metric below invariant to a column reflection.
.pefa_sign <- function(a, b) {
  scale_a <- if (length(a)) max(abs(a)) else 0
  scale_b <- if (length(b)) max(abs(b)) else 0
  if (!is.finite(scale_a) || !is.finite(scale_b) ||
      scale_a == 0 || scale_b == 0) {
    return(1)
  }
  if (sum((a / scale_a) * (b / scale_b)) >= 0) 1 else -1
}

## c(a, b) = |a'b| / (||a|| ||b||).  A zero-norm column returns NA_real_ as a
## value rather than being skipped, because section 2.6 forbids concealing an
## undefined scored congruence behind na.rm.
.pefa_congruence <- function(a, b) {
  scale_a <- if (length(a)) max(abs(a)) else 0
  scale_b <- if (length(b)) max(abs(b)) else 0
  if (!is.finite(scale_a) || !is.finite(scale_b) ||
      scale_a == 0 || scale_b == 0) {
    return(NA_real_)
  }
  unit_a <- a / scale_a
  unit_b <- b / scale_b
  value <- abs(sum(unit_a * unit_b)) /
    (sqrt(sum(unit_a^2)) * sqrt(sum(unit_b^2)))
  ## |cos| cannot exceed 1; only roundoff can push it past, and an unclamped
  ## 1 + 1e-16 would also make two exactly collinear targets compare unequal
  ## in the argmax, silently bypassing the tie rule.
  min(1, value)
}

## The largest power of two not exceeding a finite positive `x`.  Every
## rescaling below divides by this rather than by the raw magnitude, because
## dividing a double by an exact power of two only shifts its exponent: the
## scaled arithmetic then reproduces the unscaled result bit for bit wherever
## the unscaled one does not overflow, so removing the overflow moves no number
## that the ordinary fixtures pin.  log2()/floor() can only misjudge at the
## subnormal edge, where a zero or infinite scale would be worse than an
## inexact one, so the value itself is the fallback.
.pefa_pow2_scale <- function(x) {
  p <- 2^floor(log2(x))
  if (is.finite(p) && p > 0) p else x
}

## alignedRMSD(a, b) = ||a - s(a, b) b|| / sqrt(J), on a rescaled copy.
##
## Both columns are divided by one common scale BEFORE the subtraction, and the
## scale is restored only after the normalized RMSD is formed.  Rescaling the
## difference instead would be too late: a = (1e308, 1e308) against
## b = (-1e308, 1e308) has s = 1 and a perfectly representable RMSD of
## 1.41e308, but the raw difference overflows to Inf first and the whole column
## then reports NA for a number the double range holds.
.pefa_aligned_rmsd <- function(a, b, s) {
  J <- length(a)
  if (J == 0L) return(NA_real_)
  scale_ab <- max(abs(a), abs(b))
  if (!is.finite(scale_ab)) return(NA_real_)
  if (scale_ab == 0) return(0)
  scale_ab <- .pefa_pow2_scale(scale_ab)
  ## Each rescaled entry is now at most 2 in magnitude, so neither the
  ## difference nor its squares can overflow whatever the input magnitudes.
  d <- a / scale_ab - s * (b / scale_ab)
  scale_ab * sqrt(sum(d^2) / J)
}

## Adjusted Rand index between two item partitions (section 2.8), with the
## plan's exact contingency formula and roundoff clamp.
.pefa_ari <- function(p, q) {
  J <- length(p)
  if (J < 2L) return(NA_real_)
  c2 <- function(x) x * (x - 1) / 2
  tab <- table(p, q)
  A <- sum(c2(tab))
  B <- sum(c2(rowSums(tab)))
  C <- sum(c2(colSums(tab)))
  E <- B * C / c2(J)
  D <- (B + C) / 2 - E
  ## For valid partitions D = 0 only when both co-membership relations are the
  ## same degenerate partition -- all singletons, or one all-item cluster --
  ## which is agreement, not an undefined index.
  if (D == 0) return(1)
  pmin(1, pmax(-1, (A - E) / D))
}

## Every pair fact for one ordered candidate pair.  Returns the persistence row
## minus K_from/K_to/r, which the caller adds, plus the internal `.assign`
## per-column record that pair_detail() will expose in 0.9.1.
.pefa_pair_facts <- function(lam_from, lam_to, pip_from, pip_to,
                             reg_from, reg_to, K0, stability_eps,
                             loading_ok_from, loading_ok_to,
                             pip_ok_from, pip_ok_to) {
  K0 <- as.integer(K0)
  K_from <- ncol(lam_from)
  K_to <- ncol(lam_to)
  from_ok <- isTRUE(loading_ok_from)
  to_ok <- isTRUE(loading_ok_to)
  both_ok <- from_ok && to_ok

  backbone <- seq_len(K0)
  src_free <- if (K_from > K0) seq.int(K0 + 1L, K_from) else integer(0)
  tgt_free <- if (K_to > K0) seq.int(K0 + 1L, K_to) else integer(0)

  ## Eligibility counts require only their own loading endpoint (section 2.10),
  ## so a malformed target still leaves the source count available.
  eligible_src <- NULL
  if (from_ok) {
    ssl_src <- .pefa_ssl(lam_from[, src_free, drop = FALSE])
    eligible_src <- src_free[which(ssl_src >= stability_eps)]
  }
  eligible_tgt <- NULL
  if (to_ok) {
    ssl_tgt <- .pefa_ssl(lam_to[, tgt_free, drop = FALSE])
    eligible_tgt <- tgt_free[which(ssl_tgt >= stability_eps)]
  }
  source_ineligible_n <- if (from_ok) {
    length(src_free) - length(eligible_src)
  } else {
    NA_integer_
  }
  target_eligible_n <- if (to_ok) length(eligible_tgt) else NA_integer_

  ## The exact-zero boundary of $sweep$backbone_ssl_min at either endpoint.  It
  ## is a fact about the two position-matched columns being scored, so one
  ## inspectable endpoint is enough to establish TRUE, while FALSE needs both.
  backbone_zero <- function(lam) {
    ssl <- .pefa_ssl(lam[, backbone, drop = FALSE])
    any(is.finite(ssl) & ssl == 0)
  }
  backbone_degenerate <- if (K0 == 0L) {
    FALSE
  } else if ((from_ok && backbone_zero(lam_from)) ||
             (to_ok && backbone_zero(lam_to))) {
    TRUE
  } else if (both_ok) {
    FALSE
  } else {
    NA
  }

  ## Assignment (sections 2.3 and 2.4).  Each eligible exploratory source picks
  ## the eligible exploratory target of largest |phi| independently, so reuse is
  ## permitted, detected, and reported rather than prevented.  Ties go to the
  ## smaller aligned RMSD and then to the smaller target index, which makes
  ## repeat calls on identical input identical.
  assigned <- rep(NA_integer_, length(src_free))
  if (both_ok && length(eligible_src) && length(eligible_tgt)) {
    for (pos in match(eligible_src, src_free)) {
      a <- lam_from[, src_free[pos]]
      phi <- vapply(eligible_tgt,
                    function(j) .pefa_congruence(a, lam_to[, j]), numeric(1))
      usable <- which(!is.na(phi))
      if (!length(usable)) next
      best <- usable[phi[usable] == max(phi[usable])]
      if (length(best) > 1L) {
        keys <- vapply(best, function(tied) {
          b <- lam_to[, eligible_tgt[tied]]
          .pefa_aligned_rmsd(a, b, .pefa_sign(a, b))
        }, numeric(1))
        ## An undefined key cannot order the tie; the index rule then settles
        ## it, so the choice stays deterministic either way.
        if (!anyNA(keys)) best <- best[keys == min(keys)]
      }
      assigned[pos] <- eligible_tgt[best[1L]]
    }
  }

  ## P_kl (section 2.5): the K0 position-matched backbone pairs plus the
  ## exploratory pairs that could be formed.  Empty is a legitimate computed
  ## value, not an error.
  matched <- which(!is.na(assigned))
  pair_src <- if (both_ok) c(backbone, src_free[matched]) else integer(0)
  pair_tgt <- if (both_ok) c(backbone, assigned[matched]) else integer(0)
  m <- length(pair_src)
  n_pairs <- if (both_ok) m else NA_integer_

  pair_sign <- rep(NA_real_, m)
  pair_phi <- rep(NA_real_, m)
  pair_rmsd <- rep(NA_real_, m)
  if (m > 0L) {
    pair_sign <- vapply(seq_len(m), function(p) {
      .pefa_sign(lam_from[, pair_src[p]], lam_to[, pair_tgt[p]])
    }, numeric(1))
    pair_phi <- vapply(seq_len(m), function(p) {
      .pefa_congruence(lam_from[, pair_src[p]], lam_to[, pair_tgt[p]])
    }, numeric(1))
    pair_rmsd <- vapply(seq_len(m), function(p) {
      .pefa_aligned_rmsd(lam_from[, pair_src[p]], lam_to[, pair_tgt[p]],
                         pair_sign[p])
    }, numeric(1))
  }

  ## Matched-pair metrics (section 2.6).  A zero-norm scored backbone column
  ## takes phi_min out without touching the aligned quantities, which stay
  ## defined because s(a, b) = 1 at a zero dot product.
  phi_min <- NA_real_
  rmsd <- NA_real_
  rmsd_max <- NA_real_
  if (m > 0L) {
    if (!anyNA(pair_phi)) phi_min <- min(pair_phi)
    if (!anyNA(pair_rmsd)) {
      rmsd_max <- max(pair_rmsd)
      ## The sign-aligned RMSE pooled over every item-by-pair cell,
      ## sqrt(sum ||a - s b||^2 / (J * m)) -- equivalently the root mean of the
      ## squared column RMSDs, and not their average.  The sum of squares is
      ## formed on a common scale, because per-column RMSDs that are each
      ## finite can still overflow merely by being squared and added: two
      ## columns at 1.41e200 pool to 1.41e200, not to Inf.
      rmsd <- if (is.finite(rmsd_max) && rmsd_max > 0) {
        scale_r <- .pefa_pow2_scale(rmsd_max)
        scale_r * sqrt(sum((pair_rmsd / scale_r)^2) / m)
      } else {
        ## An all-zero pair set, whose exact 0 the fixtures pin, or a column
        ## RMSD that genuinely exceeded the double range: in both the unscaled
        ## form already carries the right answer.
        sqrt(sum(pair_rmsd^2) / m)
      }
    }
  }

  ## Collision and unmatched facts both need the assignment, hence both valid
  ## loadings; each is a recorded measurement, never a veto.
  collision <- NA
  n_collisions <- NA_integer_
  collision_targets <- NA_character_
  collision_multiplicities <- NA_character_
  unmatched_n <- NA_integer_
  unmatched_ssl <- NA_real_
  if (both_ok) {
    used <- assigned[matched]
    ## Reused targets, not reusing sources: three sources on one target is one
    ## collision.  The multiplicities follow the same sorted target order, so
    ## the two strings split positionally aligned.
    reused <- sort(unique(used[duplicated(used)]))
    n_collisions <- length(reused)
    collision <- n_collisions > 0L
    collision_targets <- paste(reused, collapse = ",")
    collision_multiplicities <- paste(
      vapply(reused, function(j) sum(used == j), integer(1)),
      collapse = ","
    )
    ## Section 2.7: every non-backbone target no source selected, screened-out
    ## columns included, sized by its largest column SSL rather than its sum.
    unmatched <- setdiff(tgt_free, used)
    unmatched_n <- length(unmatched)
    unmatched_ssl <- if (unmatched_n > 0L) {
      max(.pefa_ssl(lam_to[, unmatched, drop = FALSE]))
    } else {
      0
    }
  }

  ## ARI (section 2.8) deliberately uses its own independent partition column
  ## set at each endpoint -- all backbone columns plus every eligible
  ## non-backbone column there, including eligible targets no source selected.
  ari <- NA_real_
  if (both_ok) {
    part_cols_from <- c(backbone, eligible_src)
    part_cols_to <- c(backbone, eligible_tgt)
    if (length(part_cols_from) && length(part_cols_to)) {
      ## ties.method = "first" over an ascending column set is the plan's
      ## "exact ties go to the smaller original column index".
      part_from <- max.col(abs(lam_from[, part_cols_from, drop = FALSE]),
                           ties.method = "first")
      part_to <- max.col(abs(lam_to[, part_cols_to, drop = FALSE]),
                         ties.method = "first")
      ari <- .pefa_ari(part_from, part_to)
    }
  }

  ## PIP discrepancy (section 2.9) over the matched cells that are regularized
  ## at both endpoints.  Fixed cells and unmatched columns are excluded, and an
  ## unavailable PIP endpoint changes nothing loading-based.
  pip_rmsd <- NA_real_
  if (both_ok && isTRUE(pip_ok_from) && isTRUE(pip_ok_to) && m > 0L) {
    eligible_cells <- reg_from[, pair_src, drop = FALSE] &
      reg_to[, pair_tgt, drop = FALSE]
    if (any(eligible_cells)) {
      d <- pip_from[, pair_src, drop = FALSE][eligible_cells] -
        pip_to[, pair_tgt, drop = FALSE][eligible_cells]
      pip_rmsd <- sqrt(sum(d^2) / length(d))
    }
  }

  ## Matcher-input availability only.  This is the sole reason code, and no
  ## other fact -- nonconvergence, collision, a thin or degenerate backbone,
  ## ineligible sources, scarce targets -- may reach it.
  pair_status <- if (both_ok) "available" else "unavailable"
  pair_reason <- if (both_ok) NA_character_ else "malformed_loading"

  ## The per-column assignment record (section 3.4).  Backbone rows come first
  ## because the backbone occupies columns 1..K0, and they carry
  ## eligible_source = NA: the screen is defined only for exploratory columns.
  assign_target <- rep(NA_integer_, K_from)
  assign_sign <- rep(NA_real_, K_from)
  assign_phi <- rep(NA_real_, K_from)
  assign_rmsd <- rep(NA_real_, K_from)
  if (m > 0L) {
    assign_target[pair_src] <- pair_tgt
    assign_sign[pair_src] <- pair_sign
    assign_phi[pair_src] <- pair_phi
    assign_rmsd[pair_src] <- pair_rmsd
  }
  eligible_source <- rep(NA, K_from)
  if (from_ok && length(src_free)) {
    eligible_source[src_free] <- src_free %in% eligible_src
  }

  list(
    n_pairs = n_pairs,
    phi_min = phi_min,
    rmsd = rmsd,
    rmsd_max = rmsd_max,
    ari = ari,
    pip_rmsd = pip_rmsd,
    unmatched_n = unmatched_n,
    unmatched_ssl = unmatched_ssl,
    collision = collision,
    n_collisions = n_collisions,
    collision_targets = collision_targets,
    collision_multiplicities = collision_multiplicities,
    source_ineligible_n = source_ineligible_n,
    target_eligible_n = target_eligible_n,
    backbone_degenerate = backbone_degenerate,
    pair_status = pair_status,
    pair_reason = pair_reason,
    .assign = data.frame(
      source = seq_len(K_from),
      target = assign_target,
      sign = assign_sign,
      congruence = assign_phi,
      rmsd = assign_rmsd,
      backbone = seq_len(K_from) <= K0,
      eligible_source = eligible_source,
      stringsAsFactors = FALSE
    )
  )
}
