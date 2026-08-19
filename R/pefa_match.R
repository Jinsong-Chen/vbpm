## Sign-invariant Equation-17 correspondence and the seven pair measurements
## for one ordered candidate pair, K_from = k < l = K_to (note sections 4-5).
##
## The caller guarantees group-block matrices: in bifactor mode the leading
## general column has already been removed, so every column index computed and
## reported here is a group-block index.  Every entry is finite, because the
## constructor errors on a malformed candidate before reaching this file.
##
## The assignment rule is minimum sign-aligned squared distance.  Congruence is
## an exact-tie breaker and a reported shape diagnostic, never the rule itself.
## There is no eligibility screen of any kind, so every column participates and
## a comparison always holds exactly k scored pairs -- which is why no branch
## below needs an empty-set guard.

## The largest power of two not exceeding `x`, and 1 for a nonpositive or
## nonfinite `x` so a caller always receives a usable positive divisor.  Every
## rescaling below divides by this rather than by a raw magnitude, because
## dividing a double by an exact power of two only shifts its exponent: the
## scaled arithmetic reproduces the unscaled result bit for bit wherever the
## unscaled one neither overflows nor underflows, so removing either moves no
## number an ordinary fixture pins.  log2()/floor() can only misjudge at the
## subnormal edge, where a zero or infinite divisor would be worse than an
## inexact one, so the value itself is the fallback.
.pefa_pow2_scale <- function(x) {
  if (!is.finite(x) || x <= 0) return(1)
  p <- 2^floor(log2(x))
  if (is.finite(p) && p > 0) p else x
}

## s(a, b) = 1 when a'b >= 0, else -1.  The +1 convention at a'b = 0 is the
## note's arbitrary choice and is what keeps a zero-norm column's aligned
## distance defined.  The dot product is read on copies divided by one common
## power of two: an unscaled a'b can reach Inf + (-Inf) = NaN, and `NaN >= 0`
## is NA, which would turn the branch itself into an error rather than merely a
## wrong sign.  A positive divisor cannot change a sign, and negating a whole
## column negates the sum exactly, which is the root of reflection invariance.
.pefa_sign <- function(a, b) {
  scale <- .pefa_pow2_scale(max(abs(a), abs(b)))
  if (sum((a / scale) * (b / scale)) >= 0) 1 else -1
}

## d(a, b) = ||a - s b||^2 expressed on the basis `scale`, that is, the true
## aligned distance divided by scale^2.  Competing targets are ranked on this
## value directly and never on a restored one.
.pefa_sse <- function(a, b, s, scale) {
  d <- a / scale - s * (b / scale)
  sum(d^2)
}

## phi(a, b) = |a'b| / (||a|| ||b||).  A zero-norm column returns NA_real_ as a
## value rather than being skipped: note 5.1 forbids concealing an undefined
## scored congruence behind na.rm, and the raw 0/0 would be NaN rather than the
## documented NA flavor.  Congruence is scale free, so each column is
## normalized by its own peak; sum(unit^2) is then at least 1, so neither norm
## can underflow and neither product can overflow.
.pefa_congruence <- function(a, b) {
  scale_a <- max(abs(a))
  scale_b <- max(abs(b))
  if (scale_a == 0 || scale_b == 0) return(NA_real_)
  unit_a <- a / scale_a
  unit_b <- b / scale_b
  value <- abs(sum(unit_a * unit_b)) /
    (sqrt(sum(unit_a^2)) * sqrt(sum(unit_b^2)))
  ## |cos| cannot exceed 1; only roundoff can push it past, and an unclamped
  ## 1 + 1e-16 would make two exactly collinear targets compare unequal in the
  ## tie rule, pre-empting the index rule that should settle them.
  min(1, value)
}

## Adjusted Rand index between two item partitions (note 5.2).
.pefa_ari <- function(p, q) {
  J <- length(p)
  ## One item is one trivially identical partition at both endpoints.
  if (J < 2L) return(1)
  c2 <- function(x) x * (x - 1) / 2
  tab <- table(p, q)
  A <- sum(c2(tab))
  B <- sum(c2(rowSums(tab)))
  C <- sum(c2(colSums(tab)))
  E <- B * C / c2(J)
  D <- (B + C) / 2 - E
  if (D == 0) {
    ## A zero denominator forces both endpoints into the same degenerate
    ## partition -- all singletons, or one cluster -- so agreement is the only
    ## reachable answer and the 0 branch is purely defensive.  Coding each
    ## partition by order of first appearance compares them as partitions
    ## rather than as label vectors.
    same <- identical(as.integer(factor(p, levels = unique(p))),
                      as.integer(factor(q, levels = unique(q))))
    return(if (same) 1 else 0)
  }
  ## Roundoff, not mathematics, is the only way out of [-1, 1].
  min(1, max(-1, (A - E) / D))
}

## The seven structural facts for one ordered candidate pair.
.pefa_compare_pair <- function(lam_from, lam_to, pip_from, pip_to,
                               reg_from, reg_to, K0) {
  K0 <- as.integer(K0)
  J <- nrow(lam_from)
  k <- ncol(lam_from)
  l <- ncol(lam_to)
  backbone <- seq_len(K0)
  src_free <- if (k > K0) seq.int(K0 + 1L, k) else integer(0)
  tgt_free <- if (l > K0) seq.int(K0 + 1L, l) else integer(0)

  ## Assignment (note 4).  Backbone column h takes target h by position; each
  ## exploratory source then selects its own minimizer independently, so reuse
  ## is permitted, recorded, and never repaired into a one-to-one map.
  assigned <- integer(length(src_free))
  for (pos in seq_along(src_free)) {
    a <- lam_from[, src_free[pos]]
    signs <- vapply(tgt_free, function(j) .pefa_sign(a, lam_to[, j]),
                    numeric(1))
    ## One common power-of-two basis for this source and every target it may
    ## take, so the whole competing field is ranked on the same footing (note
    ## section 4).  Both ends of the double range manufacture a false exact tie
    ## otherwise: two astronomically distant targets share an Inf, and two
    ## astronomically near ones share a 0, and either hands the choice to the
    ## tie rule, which then resolves on congruence or index a distance had
    ## already settled.  Dividing by an exact power of two only shifts
    ## exponents, so a scaled SSE reproduces its unscaled value bit for bit
    ## wherever that value is representable, and one positive divisor common to
    ## the whole field is order preserving -- the mathematical ranking, exact
    ## ties included, is the unscaled one.
    scale <- .pefa_pow2_scale(max(abs(a), abs(lam_to[, tgt_free])))
    sse <- vapply(seq_along(tgt_free), function(t) {
      .pefa_sse(a, lam_to[, tgt_free[t]], signs[t], scale)
    }, numeric(1))
    best <- which(sse == min(sse))
    if (length(best) > 1L) {
      ## An exact distance tie goes to the larger finite absolute congruence,
      ## a finite congruence outranks an undefined one, and the smaller
      ## original target index settles what is left.  which.max() returns the
      ## first maximum over an ascending target set, which is that index rule.
      ## No tolerance widens a near-tie into a tie.
      phi <- vapply(best, function(t) {
        .pefa_congruence(a, lam_to[, tgt_free[t]])
      }, numeric(1))
      phi[is.na(phi)] <- -Inf
      best <- best[which.max(phi)]
    }
    assigned[pos] <- tgt_free[best[1L]]
  }

  ## Exactly k scored pairs: K0 position matches plus one per exploratory
  ## source.  Sources are distinct by construction; targets need not be.
  pair_src <- c(backbone, src_free)
  pair_tgt <- c(backbone, assigned)
  a_mat <- lam_from[, pair_src, drop = FALSE]
  b_mat <- lam_to[, pair_tgt, drop = FALSE]

  ## One common power-of-two basis for every scored cell, so a column whose
  ## aligned difference leaves the double range before the division by sqrt(J)
  ## still reports the representable RMSD instead of Inf.  Pooling and the
  ## columnwise maximum come off the same scaled sums, which is what keeps
  ## rmsd <= rmsd_max in the arithmetic and not only on paper.
  scale <- .pefa_pow2_scale(max(abs(a_mat), abs(b_mat)))
  signs <- vapply(seq_len(k), function(p) {
    .pefa_sign(a_mat[, p], b_mat[, p])
  }, numeric(1))
  sse <- vapply(seq_len(k), function(p) {
    .pefa_sse(a_mat[, p], b_mat[, p], signs[p], scale)
  }, numeric(1))
  rmsd <- scale * sqrt(sum(sse) / (J * k))
  rmsd_max <- scale * sqrt(max(sse) / J)

  ## phi_min runs over every scored pair, dead source columns included, so it
  ## is a weakest-link statistic: read a low value next to ssl() before reading
  ## it as movement in the retained structure.
  phi <- vapply(seq_len(k), function(p) {
    .pefa_congruence(a_mat[, p], b_mat[, p])
  }, numeric(1))
  phi_min <- if (anyNA(phi)) NA_real_ else min(phi)

  ## ARI uses every endpoint group column and is independent of the
  ## correspondence.  ties.method = "first" over an ascending column set is the
  ## note's "exact loading ties go to the lower column index".
  ari <- .pefa_ari(max.col(abs(lam_from), ties.method = "first"),
                   max.col(abs(lam_to), ties.method = "first"))

  ## PIPs are probabilities and are never sign-transformed.  Only cells
  ## regularized at both endpoints under the correspondence are comparable;
  ## fixed cells and unselected target columns drop out.
  sel <- reg_from[, pair_src, drop = FALSE] &
    reg_to[, pair_tgt, drop = FALSE]
  pip_rmsd <- NA_real_
  if (any(sel)) {
    d <- pip_from[, pair_src, drop = FALSE][sel] -
      pip_to[, pair_tgt, drop = FALSE][sel]
    pip_rmsd <- sqrt(sum(d^2) / length(d))
  }

  ## Every exploratory target no exploratory source selected, including the
  ## surplus a collision leaves behind.  The target has more exploratory
  ## columns than the source, so the set is nonempty and the NA_real_ branch
  ## exists only to keep max() off an empty vector.
  unmatched <- setdiff(tgt_free, assigned)
  unmatched_ssl <- if (length(unmatched)) {
    max(colSums(lam_to[, unmatched, drop = FALSE]^2))
  } else {
    NA_real_
  }
  ## Reused targets, not reusing sources; backbone positions never reuse.
  collision <- anyDuplicated(assigned) > 0L

  list(phi_min = phi_min, rmsd = rmsd, rmsd_max = rmsd_max, ari = ari,
       pip_rmsd = pip_rmsd, unmatched_ssl = unmatched_ssl,
       collision = collision)
}
