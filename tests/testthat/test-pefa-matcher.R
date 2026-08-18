## The 0.9.0 measurement kernel: sign-invariant Equation-17 correspondence and
## the seven structural facts (technical note sections 4-5, plan section 9.1).
##
## Every fixture calls the private comparator directly, so no model is fitted
## and every expected number below is hand-computable.  The comparator receives
## group-block matrices only: pefa() strips a bifactor general column before it
## gets here, so nothing in this file needs to know about bifactor mode.

## One ordered pair of candidate matrices.  PIPs and regularization masks
## default to values that keep the loading-based facts the subject of the test.
.cmp <- function(lam_from, lam_to, K0 = 0,
                 pip_from = NULL, pip_to = NULL,
                 reg_from = NULL, reg_to = NULL) {
  J <- nrow(lam_from)
  if (is.null(pip_from)) pip_from <- matrix(0.5, J, ncol(lam_from))
  if (is.null(pip_to)) pip_to <- matrix(0.5, J, ncol(lam_to))
  if (is.null(reg_from)) reg_from <- matrix(TRUE, J, ncol(lam_from))
  if (is.null(reg_to)) reg_to <- matrix(TRUE, J, ncol(lam_to))
  vbpm:::.pefa_compare_pair(lam_from, lam_to, pip_from, pip_to,
                            reg_from, reg_to, K0)
}

.facts <- c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
            "unmatched_ssl", "collision")

## Independent unsigned Equation 17: argmin_j sum((a_i - b_j)^2), with no sign
## alignment and no tie machinery.  On a fixture whose selected inner products
## are all nonnegative this must agree with the shipped rule (note invariant 1).
.eq17_unsigned <- function(lam_from, lam_to, K0) {
  J <- nrow(lam_from)
  k <- ncol(lam_from)
  l <- ncol(lam_to)
  src <- if (k > K0) seq.int(K0 + 1L, k) else integer(0)
  tgt <- if (l > K0) seq.int(K0 + 1L, l) else integer(0)
  assigned <- vapply(src, function(i) {
    d <- vapply(tgt, function(j) sum((lam_from[, i] - lam_to[, j])^2),
                numeric(1))
    tgt[which.min(d)]
  }, numeric(1))
  ps <- c(seq_len(K0), src)
  pt <- c(seq_len(K0), assigned)
  d <- lam_from[, ps, drop = FALSE] - lam_to[, pt, drop = FALSE]
  phi <- vapply(seq_len(k), function(p) {
    a <- lam_from[, ps[p]]
    b <- lam_to[, pt[p]]
    sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  }, numeric(1))
  list(assigned = assigned,
       rmsd = sqrt(sum(d^2) / (J * k)),
       rmsd_max = max(sqrt(colSums(d^2) / J)),
       phi_min = min(phi),
       unmatched_ssl = max(colSums(lam_to[, setdiff(tgt, assigned),
                                          drop = FALSE]^2)))
}


# ---- the defining behavior --------------------------------------------

test_that("minimum aligned SSE decides, not maximum congruence", {
  ## The release discriminator (note invariant 3).  `collinear` is a perfect
  ## rescaling of the source, so its congruence is exactly 1 and no target can
  ## beat it; `near_miss` is nearly the source itself, so its aligned squared
  ## distance is 100x smaller.  The two rules choose different targets.
  source_col <- c(1, 0)
  collinear <- c(2, 0)     # congruence 1,    aligned SSE 1
  near_miss <- c(1, 0.1)   # congruence < 1,  aligned SSE 0.01

  sse <- function(b) sum((source_col - b)^2)
  phi <- function(b) {
    abs(sum(source_col * b)) / sqrt(sum(source_col^2) * sum(b^2))
  }
  expect_gt(sse(collinear), sse(near_miss))     # SSE prefers near_miss
  expect_gt(phi(collinear), phi(near_miss))     # congruence prefers collinear

  f <- .cmp(cbind(source_col), cbind(collinear, near_miss))

  ## The SSE target wins.  Three independent readings pin the choice:
  ## the reported congruence is near_miss's, the pooled RMSD is near_miss's,
  ## and the leftover SSL is collinear's 2^2 rather than near_miss's 1.01.
  expect_equal(f$phi_min, 1 / sqrt(1.01))
  expect_equal(f$rmsd, sqrt(0.01 / 2))
  expect_equal(f$rmsd_max, sqrt(0.01 / 2))
  expect_equal(f$unmatched_ssl, 4)
  expect_false(f$collision)

  ## The comparator returns exactly the seven documented facts, in order.
  expect_named(f, .facts)
  expect_type(f$collision, "logical")
  expect_true(all(vapply(f[setdiff(.facts, "collision")], is.double,
                         logical(1))))
})

test_that("the rule reproduces unsigned Equation 17 when signs do not bind", {
  lam_from <- cbind(c(.9, .8, .1, 0, .1, 0),
                    c(.1, 0, .8, .7, 0, .1),
                    c(0, .1, .1, 0, .9, .8))
  lam_to <- cbind(c(.85, .75, .1, 0, .1, 0),
                  c(.05, 0, .85, .75, 0, .1),
                  c(0, .1, .1, 0, .85, .8),
                  c(.3, -.2, .1, .4, -.1, .2))
  ref <- .eq17_unsigned(lam_from, lam_to, K0 = 1)
  ## Every scored inner product is positive here, so s = +1 throughout and the
  ## shipped rule must equal the plain unsigned minimizer computed above.
  expect_identical(ref$assigned, c(2, 3))

  f <- .cmp(lam_from, lam_to, K0 = 1)
  expect_equal(f$rmsd, ref$rmsd)
  expect_equal(f$rmsd_max, ref$rmsd_max)
  expect_equal(f$phi_min, ref$phi_min)
  expect_equal(f$unmatched_ssl, ref$unmatched_ssl)
})


# ---- reflection invariance --------------------------------------------

test_that("whole-column reflection changes no fact", {
  lam_from <- cbind(c(.9, .8, .1, 0, .1, 0),
                    c(.1, 0, .8, .7, 0, .1),
                    c(0, .1, .1, 0, .9, .8))
  lam_to <- cbind(c(.85, .75, .1, 0, .1, 0),
                  c(.05, 0, .85, .75, 0, .1),
                  c(0, .1, .1, 0, .85, .8),
                  c(.3, -.2, .1, .4, -.1, .2))
  pip_from <- matrix(seq(0.05, 0.95, length.out = 18), 6, 3)
  pip_to <- matrix(seq(0.02, 0.98, length.out = 24), 6, 4)

  base <- .cmp(lam_from, lam_to, K0 = 1, pip_from = pip_from,
               pip_to = pip_to)
  expect_true(is.finite(base$phi_min))
  expect_gt(base$rmsd, 0)

  ## Every scored inner product is nonzero on this fixture, so the note
  ## requires bitwise -- not merely approximate -- invariance.  Source columns
  ## (backbone and exploratory), target columns (matched and unmatched), and
  ## then all four at once.
  for (cc in seq_len(ncol(lam_from))) {
    a <- lam_from
    a[, cc] <- -a[, cc]
    expect_identical(.cmp(a, lam_to, K0 = 1, pip_from = pip_from,
                          pip_to = pip_to), base)
  }
  for (cc in seq_len(ncol(lam_to))) {
    b <- lam_to
    b[, cc] <- -b[, cc]
    expect_identical(.cmp(lam_from, b, K0 = 1, pip_from = pip_from,
                          pip_to = pip_to), base)
  }
  a <- lam_from
  b <- lam_to
  a[, c(1L, 3L)] <- -a[, c(1L, 3L)]
  b[, c(2L, 4L)] <- -b[, c(2L, 4L)]
  expect_identical(.cmp(a, b, K0 = 1, pip_from = pip_from, pip_to = pip_to),
                   base)
})

test_that("reflection at a zero inner product holds to tolerance", {
  ## s = +1 at a'b = 0 is an arbitrary convention, so the note promises
  ## rounding-level rather than bitwise agreement here.  Both backbone pairs
  ## below are exactly orthogonal.
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  lam_to <- cbind(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 0, 1, 0))
  base <- .cmp(lam_from, lam_to, K0 = 2)
  expect_equal(base$phi_min, 0)

  for (cc in seq_len(ncol(lam_from))) {
    a <- lam_from
    a[, cc] <- -a[, cc]
    expect_equal(.cmp(a, lam_to, K0 = 2), base)
  }
  for (cc in seq_len(ncol(lam_to))) {
    b <- lam_to
    b[, cc] <- -b[, cc]
    expect_equal(.cmp(lam_from, b, K0 = 2), base)
  }
})


# ---- the tie rule -----------------------------------------------------

test_that("an exact distance tie breaks by congruence, then by index", {
  ## With ||a|| = 1 the aligned distance is 1 + ||b||^2 - 2 phi ||b||, so
  ## targets of different norm and different congruence can tie exactly.  All
  ## four targets below sit at aligned SSE exactly 1.
  a <- c(1, 0, 0, 0)
  zero <- c(0, 0, 0, 0)          # phi undefined
  half <- c(.5, .5, .5, .5)      # phi 0.5
  oneA <- c(1, 1, 0, 0)          # phi 1/sqrt(2)
  oneB <- c(1, 0, 1, 0)          # phi 1/sqrt(2), same norm, different pattern
  targets <- cbind(zero, half, oneA, oneB)
  expect_identical(apply(targets, 2L, function(b) sum((a - b)^2)),
                   c(zero = 1, half = 1, oneA = 1, oneB = 1))

  ## PIP column 4 differs from the source's, so pip_rmsd reads out whether the
  ## third or the fourth target was taken.
  pip_to <- cbind(rep(.5, 4), rep(.5, 4), rep(.5, 4), rep(.9, 4))
  f <- .cmp(cbind(a), targets, pip_to = pip_to)

  ## Largest finite congruence wins, so the two 1/sqrt(2) targets beat both
  ## the undefined zero column and the 0.5 column ...
  expect_equal(f$phi_min, 1 / sqrt(2))
  ## ... and the remaining tie goes to the lower original target index.
  expect_equal(f$pip_rmsd, 0)
  ## Leftovers are zero (0), half (1) and oneB (2): the max is oneB's.
  expect_equal(f$unmatched_ssl, 2)
})

test_that("a finite congruence outranks an undefined one in a tie", {
  a <- c(1, 0, 0, 0)
  zero <- c(0, 0, 0, 0)
  half <- c(.5, .5, .5, .5)
  expect_identical(sum((a - zero)^2), sum((a - half)^2))

  ## The zero column sits at the lower index, so only the congruence
  ## preference can move the choice off it.
  f <- .cmp(cbind(a), cbind(zero, half))
  expect_equal(f$phi_min, 0.5)
  expect_identical(f$unmatched_ssl, 0)

  ## When every tied congruence is undefined, the index rule settles it: the
  ## first zero column is taken and the second is the leftover.
  f0 <- .cmp(cbind(a), cbind(zero, zero))
  expect_identical(f0$phi_min, NA_real_)
  expect_identical(f0$unmatched_ssl, 0)
})

test_that("a near tie is not a tie", {
  ## Same fixture as the exact tie, except the high-congruence target is moved
  ## by 1e-9 so its aligned SSE is strictly larger.  No tolerance may widen
  ## that into a tie: the 0.5-congruence target must now win on distance.
  a <- c(1, 0, 0, 0)
  nudged <- c(1, 1 + 1e-9, 0, 0)
  half <- c(.5, .5, .5, .5)
  expect_gt(sum((a - nudged)^2), sum((a - half)^2))

  f <- .cmp(cbind(a), cbind(nudged, half))
  expect_equal(f$phi_min, 0.5)
  expect_equal(f$unmatched_ssl, sum(nudged^2))
})


# ---- no eligibility screen --------------------------------------------

test_that("weak, singleton and exact-zero columns all participate", {
  ## A salience screen of any size would drop the 1e-8 column, and phi_min
  ## would come back as 1 from the two clean pairs.  Every column is scored,
  ## so the dead column's arbitrary congruence is what is reported.
  weak <- c(1e-8, 0, 0, 0)
  lam_from <- cbind(c(.9, .9, 0, 0), weak, c(0, 0, 0, 1))
  lam_to <- cbind(c(.9, .9, 0, 0), c(0, 0, 1e-9, 0), c(0, 0, 0, 1),
                  c(.4, .4, .4, .4))

  f <- .cmp(lam_from, lam_to)
  ## The near-dead source matches the smallest-norm target, exactly as the
  ## note predicts, and drags phi_min to 0.
  expect_equal(f$phi_min, 0)
  ## Exactly k = 3 pairs are scored: pooled RMSD is the 1e-8 column's residual
  ## spread over J * k cells, which no screened-out column could produce.
  expect_equal(f$rmsd, sqrt((1e-16 + 1e-18) / 12))
  expect_equal(f$rmsd_max, sqrt((1e-16 + 1e-18) / 4))
  ## The only unselected target is the dense one.
  expect_equal(f$unmatched_ssl, 4 * 0.16)
  expect_false(f$collision)
})

test_that("a zero-norm source keeps RMSD finite and types phi_min as NA", {
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 0, 0, 0))
  lam_to <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
  f <- .cmp(lam_from, lam_to)

  ## NA_real_, not NaN: the undefined congruence is coerced rather than
  ## dropped with na.rm.
  expect_identical(f$phi_min, NA_real_)
  expect_false(is.nan(f$phi_min))
  ## Aligned distance is defined for a zero column, so both RMSDs stay finite.
  expect_equal(f$rmsd, sqrt(1 / 8))
  expect_equal(f$rmsd_max, sqrt(1 / 4))
  ## All three unit targets are equidistant from the dead column, so the index
  ## rule sends it onto the first one and raises a collision.
  expect_true(f$collision)
  expect_equal(f$unmatched_ssl, 1)
})


# ---- collision and leftovers ------------------------------------------

test_that("target reuse is recorded and never repaired", {
  ## Both sources are nearest to target 1, so the independent assignment
  ## reuses it instead of pushing the second source onto a runner-up.
  lam_from <- cbind(c(1, 0, 0, 0), c(.9, 0, 0, 0))
  lam_to <- cbind(c(1, 0, 0, 0), c(0, 3, 0, 0), c(0, 0, 2, 0))
  f <- .cmp(lam_from, lam_to)

  expect_true(f$collision)
  ## Both leftovers count, and the reported value is the largest SSL (9), not
  ## the sum (13) and not the smaller one (4).
  expect_equal(f$unmatched_ssl, 9)
  ## A collision blanks nothing: every numeric fact is still defined.
  expect_equal(f$phi_min, 1)
  expect_equal(f$rmsd, sqrt(0.01 / 8))
  expect_equal(f$rmsd_max, sqrt(0.01 / 4))
})

test_that("a clean one-to-one map reports collision FALSE", {
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  lam_to <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
  f <- .cmp(lam_from, lam_to)
  expect_false(f$collision)
  expect_identical(f$rmsd, 0)
  expect_identical(f$rmsd_max, 0)
  expect_equal(f$phi_min, 1)
  expect_equal(f$unmatched_ssl, 1)
})

test_that("an unused exact-zero target gives unmatched_ssl exactly zero", {
  f <- .cmp(cbind(c(1, 0)), cbind(c(1, 0), c(0, 0)))
  expect_identical(f$unmatched_ssl, 0)
  expect_false(f$collision)
})


# ---- the remaining metrics --------------------------------------------

test_that("ari compares max-loading partitions with lower-index ties", {
  ## Rows 3 and 4 are all-zero at the source and row 4 is all-zero at the
  ## target, so each is an exact loading tie resolved to the lowest column.
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  lam_to <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
  expect_identical(max.col(abs(lam_from), ties.method = "first"),
                   c(1L, 2L, 1L, 1L))
  expect_identical(max.col(abs(lam_to), ties.method = "first"),
                   c(1L, 2L, 3L, 1L))
  ## Contingency counts give A = 1, B = 3, C = 1, E = 1/2, D = 3/2.
  expect_equal(.cmp(lam_from, lam_to)$ari, (1 - 0.5) / 1.5)

  ## ARI reads every endpoint group column, not just the matched ones.
  ## Emptying the unselected third target leaves the assignment and both RMSDs
  ## untouched, yet moves the target partition and therefore the index.
  lam_dead <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 0))
  dead <- .cmp(lam_from, lam_dead)
  expect_identical(dead$rmsd, .cmp(lam_from, lam_to)$rmsd)
  expect_identical(dead$ari, 1)

  ## A negative index is reachable and is reported as such.
  lam_to2 <- cbind(c(1, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 1, 0))
  expect_equal(.cmp(lam_from, lam_to2)$ari, (1 - 1.5) / 1.5)
})

test_that("a one-item comparison returns the documented ari of 1", {
  f <- .cmp(matrix(1, 1, 1), matrix(c(1, 2), 1, 2))
  expect_identical(f$ari, 1)
})

test_that("pip_rmsd uses only cells regularized at both endpoints", {
  lam_from <- cbind(c(1, 1, 0, 0), c(0, 0, 1, 0))
  lam_to <- cbind(c(1, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
  ## Rows 1-2 of the backbone are fixed anchors at both endpoints, so only
  ## 2 backbone cells and 4 exploratory cells are jointly regularized.
  reg_from <- cbind(c(FALSE, FALSE, TRUE, TRUE), rep(TRUE, 4))
  reg_to <- cbind(c(FALSE, FALSE, TRUE, TRUE), rep(TRUE, 4), rep(TRUE, 4))
  pip_from <- cbind(c(.9, .9, .1, .2), c(.05, .05, .95, .10))
  pip_to <- cbind(c(.8, .8, .3, .2), c(.05, .15, .85, .10),
                  c(.5, .5, .5, .5))

  sel <- reg_from & reg_to[, c(1L, 2L)]
  expect_identical(sum(sel), 6L)
  d <- pip_from[sel] - pip_to[, c(1L, 2L)][sel]
  expect_equal(sum(d^2), 0.06)

  f <- .cmp(lam_from, lam_to, K0 = 1, pip_from = pip_from, pip_to = pip_to,
            reg_from = reg_from, reg_to = reg_to)
  expect_equal(f$pip_rmsd, sqrt(0.06 / 6))
  ## Column 3 of the target is unmatched, so its PIPs never enter.
  expect_equal(f$unmatched_ssl, 1)
})

test_that("pip_rmsd is NA when no cell is regularized at both endpoints", {
  lam_from <- cbind(c(1, 0), c(0, 1))
  lam_to <- cbind(c(1, 0), c(0, 1), c(1, 1))
  f <- .cmp(lam_from, lam_to,
            reg_from = matrix(FALSE, 2, 2), reg_to = matrix(TRUE, 2, 3))
  expect_identical(f$pip_rmsd, NA_real_)
  ## The loading facts stay defined.
  expect_identical(f$rmsd, 0)
})

test_that("pooled rmsd never exceeds the columnwise maximum", {
  fixtures <- list(
    list(cbind(c(1, 0, 0, 0), c(0, 1, 0, 0)),
         cbind(c(1, 0, 0, 0), c(0, 0, 0, 0), c(0, 0, 1, 0))),
    list(cbind(c(.9, .8, .1, 0), c(.1, 0, .8, .7)),
         cbind(c(.85, .75, .1, 0), c(.05, 0, .85, .75), c(.3, -.2, .1, .4))),
    list(cbind(c(1, 0, 0, 0), c(.9, 0, 0, 0)),
         cbind(c(1, 0, 0, 0), c(0, 3, 0, 0), c(0, 0, 2, 0))),
    list(cbind(c(1, 1, 1, 1)), cbind(c(1, 1, 1, 1), c(0, 0, 0, 1)))
  )
  for (fx in fixtures) {
    f <- .cmp(fx[[1L]], fx[[2L]])
    ## A mean cannot exceed a maximum in exact arithmetic, but two equal
    ## columns can invert the pair by one ulp, so this is a tolerance claim.
    expect_lte(f$rmsd, f$rmsd_max + 1e-12)
  }

  ## Unequal columns separate the two by a hand-computable margin: aligned
  ## squared distances 0 and 1 over J = 4 give sqrt(1/8) pooled and sqrt(1/4)
  ## as the weakest column.
  f <- .cmp(fixtures[[1L]][[1L]], fixtures[[1L]][[2L]])
  expect_equal(f$rmsd, sqrt(1 / 8))
  expect_equal(f$rmsd_max, sqrt(1 / 4))
  expect_lt(f$rmsd, f$rmsd_max)
})


# ---- backbone handling ------------------------------------------------

test_that("backbone columns pair by position, never by fit", {
  ## The target's first two columns are the source's in swapped order, so a
  ## free matcher would swap them and report RMSD 0.  Position pairing must
  ## instead report the orthogonal-pair distance.
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  lam_to <- cbind(c(0, 1, 0, 0), c(1, 0, 0, 0), c(0, 0, 1, 0))
  f <- .cmp(lam_from, lam_to, K0 = 2)

  expect_equal(f$rmsd, sqrt(0.5))
  expect_equal(f$rmsd_max, sqrt(0.5))
  expect_equal(f$phi_min, 0)
  ## k = K0 is a legal edge: there are no exploratory sources at all, so no
  ## reuse is possible and the single spare target is the leftover.
  expect_false(f$collision)
  expect_equal(f$unmatched_ssl, 1)
})

test_that("backbone position pairs are immune to a better exploratory target", {
  ## Column 3 of the target is an exact copy of backbone source column 1, yet
  ## the backbone still pairs 1 with 1.
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 0, 1, 0))
  lam_to <- cbind(c(.5, 0, 0, 0), c(0, 0, 1, 0), c(1, 0, 0, 0))
  f <- .cmp(lam_from, lam_to, K0 = 1)
  ## Backbone pair distance 0.25, exploratory pair distance 0.
  expect_equal(f$rmsd, sqrt(0.25 / 8))
  expect_equal(f$rmsd_max, sqrt(0.25 / 4))
  ## The exploratory source took target 2, so target 3 is the leftover.
  expect_equal(f$unmatched_ssl, 1)
})


# ---- direct, never chained --------------------------------------------

test_that("a direct comparison ignores what an intervening candidate matched", {
  ## L1 -> L2 sends the single source onto column 1 of L2; L2 -> L3 sends that
  ## column onto column 2 of L3.  Composing those two assignments would land
  ## on L3 column 2, but the direct L1 -> L3 comparison lands on column 1.
  L1 <- cbind(c(1, 0, 0))
  L2 <- cbind(c(1, 1, 0), c(0, 0, 1))
  L3 <- cbind(c(1, 0, 0), c(1, 1, 0), c(0, 0, 1))

  ## Leg 1: the source is closer to L2 column 1 than to column 2.
  expect_equal(.cmp(L1, L2)$rmsd, sqrt(1 / 3))
  ## Leg 2: that column is an exact copy of L3 column 2, so it goes there.
  expect_identical(.cmp(cbind(L2[, 1L]), L3)$rmsd, 0)

  direct <- .cmp(L1, L3)
  ## Direct: an exact copy of the source sits at L3 column 1.
  expect_identical(direct$rmsd, 0)
  expect_identical(direct$rmsd_max, 0)
  ## The leftovers are L3 columns 2 and 3, whose SSLs are 2 and 1.  Had the
  ## chained answer been used the leftovers would be columns 1 and 3, so
  ## unmatched_ssl would have been 1.
  expect_equal(direct$unmatched_ssl, 2)
  expect_equal(.cmp(L1, cbind(L3[, 2L]))$rmsd, sqrt(1 / 3))
})


# ---- numerical robustness ---------------------------------------------

test_that("extreme finite columns do not overflow into a manufactured tie", {
  ## Both unscaled aligned distances overflow to Inf, which without rescaling
  ## would tie and hand the choice to the tie rule.  On a common power-of-two
  ## basis the nearer target is still strictly nearer.
  a <- c(1e200, 0)
  near <- c(0, 1e200)
  far <- c(0, 2e200)
  expect_identical(sum((a - near)^2), Inf)
  expect_identical(sum((a - far)^2), Inf)

  f <- .cmp(cbind(a), cbind(near, far))
  ## sqrt(2) * 1e200 spread over J = 2 cells is exactly 1e200.
  expect_equal(f$rmsd, 1e200)
  expect_equal(f$rmsd_max, 1e200)
  expect_true(is.finite(f$rmsd))
  expect_false(f$collision)

  ## Reversing the target order reverses nothing: the choice is a distance,
  ## not an index accident.
  f2 <- .cmp(cbind(a), cbind(far, near))
  expect_equal(f2$rmsd, f$rmsd)
})

test_that("a huge but representable RMSD is formed instead of Inf", {
  ## The aligned sum of squares is 2e600 and cannot be held in a double, but
  ## the RMSD that comes out of it can, so scaling has to happen before the
  ## division rather than after it.
  a <- c(1e300, 0, 0, 0)
  targets <- cbind(c(0, 1e300, 0, 0), c(0, 0, 9e300, 0))
  f <- .cmp(cbind(a), targets)
  expect_true(is.finite(f$rmsd))
  expect_equal(f$rmsd, sqrt(2) / 2 * 1e300)
  expect_equal(f$rmsd_max, sqrt(2) / 2 * 1e300)
})
