## Section 2 of the 0.9.0 contract: pairwise sign alignment, the max-|phi|
## assignment and its tie rule, the scored pair set P_kl, and the
## metric-by-metric availability rules.  Every fixture calls the internal
## matcher directly, so the numbers below are hand-computable and no model is
## fitted.

## One ordered pair of loading matrices, with PIP/regularization defaults that
## keep the loading-based fields the subject of the test.
.pf <- function(lam_from, lam_to, K0, eps = 0.1,
                pip_from = NULL, pip_to = NULL,
                reg_from = NULL, reg_to = NULL,
                ok_from = TRUE, ok_to = TRUE,
                pip_ok_from = TRUE, pip_ok_to = TRUE) {
  J <- nrow(lam_from)
  if (is.null(pip_from)) pip_from <- matrix(0.5, J, ncol(lam_from))
  if (is.null(pip_to)) pip_to <- matrix(0.5, J, ncol(lam_to))
  if (is.null(reg_from)) reg_from <- matrix(TRUE, J, ncol(lam_from))
  if (is.null(reg_to)) reg_to <- matrix(TRUE, J, ncol(lam_to))
  vbpm:::.pefa_pair_facts(lam_from, lam_to, pip_from, pip_to,
                          reg_from, reg_to, K0, eps,
                          ok_from, ok_to, pip_ok_from, pip_ok_to)
}

## Everything the matcher returns except the internal assignment record.
.pf_metrics <- function(facts) facts[setdiff(names(facts), ".assign")]

.unit <- function(i, J = 8) {
  v <- numeric(J)
  v[i] <- 1
  v
}

test_that("reflecting a whole column changes no metric", {
  lam_from <- cbind(c(.9, .8, .1, 0, .1, 0),
                    c(.1, 0, .8, .7, 0, .1),
                    c(0, .1, .1, 0, .9, .8))
  lam_to <- cbind(c(.85, .75, .1, 0, .1, 0),
                  c(.05, 0, .85, .75, 0, .1),
                  c(0, .1, .1, 0, .85, .8),
                  c(.3, -.2, .1, .4, -.1, .2))
  pip_from <- matrix(seq(0.05, 0.95, length.out = 18), 6, 3)
  pip_to <- matrix(seq(0.02, 0.98, length.out = 24), 6, 4)

  base <- .pf_metrics(.pf(lam_from, lam_to, K0 = 1,
                          pip_from = pip_from, pip_to = pip_to))
  expect_identical(base$pair_status, "available")
  expect_true(is.finite(base$phi_min))

  ## Backbone at either endpoint, exploratory at either endpoint, and all four
  ## at once: s(a, b) absorbs every reflection, so nothing may move.
  flips <- list(from = 1L, from = 3L, to = 1L, to = 4L)
  for (i in seq_along(flips)) {
    a <- lam_from
    b <- lam_to
    if (names(flips)[i] == "from") {
      a[, flips[[i]]] <- -a[, flips[[i]]]
    } else {
      b[, flips[[i]]] <- -b[, flips[[i]]]
    }
    expect_identical(
      .pf_metrics(.pf(a, b, K0 = 1, pip_from = pip_from, pip_to = pip_to)),
      base
    )
  }

  a <- lam_from
  b <- lam_to
  a[, c(1L, 3L)] <- -a[, c(1L, 3L)]
  b[, c(1L, 4L)] <- -b[, c(1L, 4L)]
  expect_identical(
    .pf_metrics(.pf(a, b, K0 = 1, pip_from = pip_from, pip_to = pip_to)),
    base
  )
})

test_that("assignment maximizes absolute congruence, not minus SSE", {
  source_col <- c(1, 0)
  collinear <- c(2, 0)     # congruence 1, aligned SSE 1
  near_miss <- c(1, 0.1)   # congruence < 1, aligned SSE 0.01

  ## The two rules genuinely disagree on this fixture.
  expect_gt(sum((source_col - collinear)^2), sum((source_col - near_miss)^2))
  expect_gt(
    abs(sum(source_col * collinear)) /
      sqrt(sum(source_col^2) * sum(collinear^2)),
    abs(sum(source_col * near_miss)) /
      sqrt(sum(source_col^2) * sum(near_miss^2))
  )

  facts <- .pf(cbind(source_col), cbind(collinear, near_miss), K0 = 0,
               eps = 0.001)
  expect_identical(facts$.assign$target, 1L)
  expect_equal(facts$phi_min, 1, tolerance = 1e-12)
  ## The unselected near-miss column is the unmatched one.
  expect_identical(facts$unmatched_n, 1L)
  expect_equal(facts$unmatched_ssl, 1.01, tolerance = 1e-12)
})

test_that("ties resolve by aligned RMSD, then by target index", {
  source_col <- c(1, 1, 0, 0)
  scaled <- 2 * source_col   # congruence exactly 1, aligned RMSD sqrt(2)/2
  exact <- source_col        # congruence exactly 1, aligned RMSD 0
  duplicate <- source_col    # identical to `exact`, higher index

  targets <- cbind(scaled, exact, duplicate)
  facts <- .pf(cbind(source_col), targets, K0 = 0, eps = 0.001)
  ## The RMSD key beats the index rule: the exact copy sits at index 2.
  expect_identical(facts$.assign$target, 2L)
  expect_identical(facts$rmsd, 0)
  ## `exact` and `duplicate` tie on both keys, so the smaller index wins.
  expect_identical(
    .pf(cbind(source_col), targets, K0 = 0, eps = 0.001)$.assign$target,
    facts$.assign$target
  )
  expect_identical(.pf_metrics(.pf(cbind(source_col), targets, K0 = 0,
                                   eps = 0.001)),
                   .pf_metrics(facts))
})

test_that("P_kl carries the backbone and a k = K0 edge stays finite", {
  lam_from <- cbind(c(.9, .1, 0, 0), c(0, 0, .8, .2))
  lam_to <- cbind(c(.85, .05, 0, 0), c(0, 0, .75, .25),
                  c(.4, 0, .4, 0), c(0, .5, 0, .5))
  facts <- .pf(lam_from, lam_to, K0 = 2)

  expect_identical(facts$n_pairs, 2L)          # backbone only, complete
  expect_identical(facts$pair_status, "available")
  expect_true(is.na(facts$pair_reason))
  expect_true(is.finite(facts$phi_min))
  expect_true(is.finite(facts$rmsd_max))
  expect_false(any(is.infinite(c(facts$phi_min, facts$rmsd,
                                 facts$rmsd_max, facts$unmatched_ssl))))
  expect_identical(facts$n_collisions, 0L)
  expect_false(facts$collision)
  expect_identical(facts$source_ineligible_n, 0L)
  expect_identical(facts$unmatched_n, 2L)
  expect_identical(sum(!is.na(facts$.assign$target)), facts$n_pairs)
  expect_true(all(facts$.assign$backbone))
})

test_that("rmsd is the sign-aligned pooled RMSE bounded by rmsd_max", {
  ## Two position-matched backbone pairs: the first identical, the second
  ## differing in exactly one of four items.
  lam_from <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  lam_to <- cbind(c(1, 0, 0, 0), c(0, 1, 1, 0))
  facts <- .pf(lam_from, lam_to, K0 = 2)

  column_rmsd <- c(0, sqrt(1 / 4))
  expect_identical(facts$n_pairs, 2L)
  expect_equal(facts$rmsd_max, max(column_rmsd), tolerance = 1e-12)
  expect_equal(facts$rmsd, sqrt(sum(column_rmsd^2) / 2), tolerance = 1e-12)
  expect_equal(facts$rmsd, sqrt((0 + 1) / (4 * 2)), tolerance = 1e-12)
  expect_lte(facts$rmsd, facts$rmsd_max)
  ## Not the mean of the column RMSDs, which this fixture separates.
  expect_false(isTRUE(all.equal(facts$rmsd, mean(column_rmsd))))
  expect_equal(facts$phi_min, 1 / sqrt(2), tolerance = 1e-12)
})

test_that("ARI matches a hand-computed value and both degeneracies give one", {
  ari <- vbpm:::.pefa_ari
  ## p = (1,1,2,2), q = (1,2,2,2): A = 1, B = 2, C = 3, E = 1, D = 1.5.
  expect_equal(ari(c(1, 1, 2, 2), c(1, 2, 2, 2)), 0, tolerance = 1e-12)
  ## Both realizable D = 0 partitions are identities, so both return 1.
  expect_identical(ari(1:4, 1:4), 1)
  expect_identical(ari(rep(1, 4), rep(1, 4)), 1)
  ## Identity of co-membership, not of labels.
  expect_equal(ari(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1, tolerance = 1e-12)
  expect_identical(ari(1L, 1L), NA_real_)          # J < 2

  ## The roundoff clamp keeps every value inside [-1, 1].
  set.seed(11)
  sampled <- vapply(seq_len(200), function(i) {
    ari(sample(1:3, 12, TRUE), sample(1:4, 12, TRUE))
  }, numeric(1))
  expect_true(all(is.finite(sampled)))
  expect_true(all(sampled <= 1 & sampled >= -1))

  ## An unavailable endpoint gives the typed NA through the matcher.
  lam_from <- cbind(c(.9, .1, 0, 0), c(0, 0, .8, .2))
  lam_to <- cbind(c(.85, .05, 0, 0), c(0, 0, .75, .25), c(.4, 0, .4, 0))
  expect_identical(.pf(lam_from, lam_to, K0 = 1, ok_to = FALSE)$ari, NA_real_)
})

test_that("pip_rmsd averages only the jointly regularized matched cells", {
  ## Q0 has an anchor 1, an AZ fixed 0, and two regularized backbone cells.
  Q0 <- matrix(c(1L, 0L, -1L, -1L), 4, 1)
  reg_from <- cbind(Q0 == -1L, matrix(TRUE, 4, 1))
  reg_to <- cbind(Q0 == -1L, matrix(TRUE, 4, 2))
  lam_from <- cbind(c(.9, .8, .2, .1), c(0, 0, .7, .7))
  lam_to <- cbind(c(.85, .75, .2, .1), c(0, 0, .72, .68), c(.1, -.1, .1, -.1))
  pip_from <- matrix(0.5, 4, 2)
  pip_to <- cbind(rep(0.4, 4), rep(0.6, 4), rep(0.9, 4))

  facts <- .pf(lam_from, lam_to, K0 = 1, pip_from = pip_from,
               pip_to = pip_to, reg_from = reg_from, reg_to = reg_to)
  expect_identical(facts$.assign$target, c(1L, 2L))
  ## Two regularized backbone cells plus four exploratory cells.
  expect_identical(sum(reg_from[, c(1L, 2L)] & reg_to[, c(1L, 2L)]), 6L)
  expect_equal(facts$pip_rmsd, 0.1, tolerance = 1e-12)

  ## A fixed Q0 cell is excluded: moving it cannot move the number.
  moved_fixed <- pip_to
  moved_fixed[1:2, 1] <- 0
  expect_equal(
    .pf(lam_from, lam_to, K0 = 1, pip_from = pip_from, pip_to = moved_fixed,
        reg_from = reg_from, reg_to = reg_to)$pip_rmsd,
    0.1, tolerance = 1e-12
  )
  ## So is the unmatched target column.
  moved_unmatched <- pip_to
  moved_unmatched[, 3] <- 0
  expect_equal(
    .pf(lam_from, lam_to, K0 = 1, pip_from = pip_from,
        pip_to = moved_unmatched, reg_from = reg_from,
        reg_to = reg_to)$pip_rmsd,
    0.1, tolerance = 1e-12
  )
  ## An unavailable PIP endpoint blanks only pip_rmsd.
  degraded <- .pf(lam_from, lam_to, K0 = 1, pip_from = pip_from,
                  pip_to = pip_to, reg_from = reg_from, reg_to = reg_to,
                  pip_ok_to = FALSE)
  expect_identical(degraded$pip_rmsd, NA_real_)
  expect_true(is.finite(degraded$phi_min))
  expect_identical(degraded$pair_status, "available")
})

test_that("unmatched targets are sized by the largest column SSL", {
  lam_from <- cbind(.unit(1), .unit(2))
  lam_to <- cbind(.unit(1), .unit(2), .5 * .unit(4), .8 * .unit(4))
  facts <- .pf(lam_from, lam_to, K0 = 1)

  ssl_unmatched <- c(0.25, 0.64)
  expect_identical(facts$unmatched_n, 2L)
  expect_equal(facts$unmatched_ssl, max(ssl_unmatched), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(facts$unmatched_ssl, sum(ssl_unmatched))))

  ## An empty unmatched set is a known zero, not NA.
  full <- .pf(cbind(.unit(1), .unit(2)), cbind(.unit(1), .unit(2)), K0 = 1)
  expect_identical(full$unmatched_n, 0L)
  expect_identical(full$unmatched_ssl, 0)
})

test_that("collision counts reused targets and locates them", {
  ## Sources 2 and 3 both take group target 2; sources 4, 5 and 6 all take
  ## group target 4.  Targets 3, 5 and 6 are unselected, and target 6 also
  ## falls below the screen, so it must still be counted as unmatched.
  lam_from <- cbind(.unit(1), .unit(2), .unit(2), .unit(3), .unit(3),
                    .unit(3))
  lam_to <- cbind(.unit(1), .unit(2), .5 * .unit(4), .unit(3), .8 * .unit(4),
                  .1 * .unit(4))
  facts <- .pf(lam_from, lam_to, K0 = 1)

  expect_true(facts$collision)
  expect_identical(facts$n_collisions, 2L)          # two reused TARGETS
  expect_identical(facts$collision_targets, "2,4")
  expect_identical(facts$collision_multiplicities, "2,3")
  expect_identical(
    strsplit(facts$collision_multiplicities, ",", fixed = TRUE)[[1L]],
    c("2", "3")
  )
  expect_identical(facts$n_pairs, 6L)
  expect_identical(sum(!is.na(facts$.assign$target)), facts$n_pairs)
  expect_identical(facts$.assign$target, c(1L, 2L, 2L, 4L, 4L, 4L))
  expect_identical(facts$pair_status, "available")
  expect_true(is.finite(facts$phi_min))
  expect_identical(facts$target_eligible_n, 4L)
  expect_identical(facts$unmatched_n, 3L)           # includes the screened one
  expect_equal(facts$unmatched_ssl, 0.64, tolerance = 1e-12)

  ## Three sources on one target is one collision, not three.
  one_target <- .pf(cbind(.unit(1), .unit(2), .unit(2), .unit(2)),
                    cbind(.unit(1), .unit(2), .unit(5)), K0 = 1)
  expect_identical(one_target$n_collisions, 1L)
  expect_identical(one_target$collision_targets, "2")
  expect_identical(one_target$collision_multiplicities, "3")

  ## A clean computed edge stores two empty strings, never NA.
  clean <- .pf(cbind(.unit(1), .unit(2)), cbind(.unit(1), .unit(2), .unit(5)),
               K0 = 1)
  expect_false(clean$collision)
  expect_identical(clean$n_collisions, 0L)
  expect_identical(clean$collision_targets, "")
  expect_identical(clean$collision_multiplicities, "")
})

test_that("eligibility is exact at the stability screen", {
  ## SSL exactly 0.25, and SSL exactly one ulp below it.
  at_screen <- c(0.5, 0, 0, 0)
  below_screen <- c(sqrt(0.25 - 2^-54), 0, 0, 0)
  expect_identical(sum(at_screen^2), 0.25)
  expect_identical(sum(below_screen^2), 0.25 - 2^-54)
  expect_lt(sum(below_screen^2), 0.25)

  backbone <- c(1, 0, 0, 0)
  targets <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  expect_identical(
    .pf(cbind(backbone, at_screen), targets, K0 = 1,
        eps = 0.25)$source_ineligible_n, 0L)
  expect_identical(
    .pf(cbind(backbone, below_screen), targets, K0 = 1,
        eps = 0.25)$source_ineligible_n, 1L)
  expect_identical(
    .pf(cbind(backbone), cbind(backbone, at_screen), K0 = 1,
        eps = 0.25)$target_eligible_n, 1L)
  expect_identical(
    .pf(cbind(backbone), cbind(backbone, below_screen), K0 = 1,
        eps = 0.25)$target_eligible_n, 0L)

  ## The same boundary at the package default: SSL exactly 0.10 is eligible,
  ## and one ulp of separation is enough to exclude it.
  at_default <- c(sqrt(0.1), 0, 0, 0)
  expect_identical(sum(at_default^2), 0.1)
  expect_gt(0.1 + 2^-56, 0.1)
  expect_identical(
    .pf(cbind(backbone, at_default), targets, K0 = 1,
        eps = 0.1)$source_ineligible_n, 0L)
  expect_identical(
    .pf(cbind(backbone, at_default), targets, K0 = 1,
        eps = 0.1 + 2^-56)$source_ineligible_n, 1L)
})

test_that("a zero-norm backbone is metric-specific, not row-wide", {
  lam_from <- cbind(rep(0, 4), c(0, 0, 1, 0))
  lam_to <- cbind(c(.9, .1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
  facts <- .pf(lam_from, lam_to, K0 = 1)

  expect_identical(facts$pair_status, "available")
  expect_identical(facts$pair_reason, NA_character_)
  expect_true(facts$backbone_degenerate)
  expect_identical(facts$phi_min, NA_real_)   # one scored congruence undefined
  expect_true(is.finite(facts$rmsd))
  expect_true(is.finite(facts$rmsd_max))
  expect_true(is.finite(facts$ari))
  expect_true(is.finite(facts$pip_rmsd))
  expect_false(facts$collision)
  expect_identical(facts$n_collisions, 0L)
  expect_identical(facts$unmatched_n, 1L)
  expect_true(is.finite(facts$unmatched_ssl))
  expect_identical(facts$n_pairs, 2L)

  ## K0 = 0 has no backbone to degenerate.
  expect_false(.pf(cbind(c(1, 0, 0, 0)), cbind(c(1, 0, 0, 0), c(0, 1, 0, 0)),
                   K0 = 0)$backbone_degenerate)
})

test_that("an empty scored pair set is a computed value", {
  screened <- c(.05, 0, 0, 0)
  facts <- .pf(cbind(screened, screened),
               cbind(c(1, 0, 0, 0), c(0, 1, 0, 0)), K0 = 0)

  expect_identical(facts$n_pairs, 0L)
  expect_identical(facts$phi_min, NA_real_)
  expect_identical(facts$rmsd, NA_real_)
  expect_identical(facts$rmsd_max, NA_real_)
  expect_identical(facts$pip_rmsd, NA_real_)
  expect_identical(facts$ari, NA_real_)
  expect_false(facts$collision)
  expect_identical(facts$n_collisions, 0L)
  expect_identical(facts$collision_targets, "")
  expect_identical(facts$collision_multiplicities, "")
  expect_false(facts$backbone_degenerate)
  expect_identical(facts$source_ineligible_n, 2L)
  expect_identical(facts$unmatched_n, 2L)
  expect_equal(facts$unmatched_ssl, 1, tolerance = 1e-12)
  expect_identical(facts$pair_status, "available")
  expect_true(all(is.na(facts$.assign$target)))

  ## The same K0 = 0 shape with eligible columns returns ordinary metrics.
  ordinary <- .pf(cbind(c(1, 0, 0, 0)),
                  cbind(c(0, 1, 0, 0), c(1, 0, 0, 0)), K0 = 0)
  expect_identical(ordinary$n_pairs, 1L)
  expect_equal(ordinary$phi_min, 1, tolerance = 1e-12)
  expect_identical(ordinary$rmsd, 0)
  expect_identical(ordinary$.assign$target, 2L)
})

test_that("an ineligible source or scarce targets leave the pair available", {
  ## One eligible exploratory source, no eligible exploratory target.
  facts <- .pf(cbind(c(.9, .1, 0, 0), c(0, 0, 1, 0)),
               cbind(c(.85, .1, 0, 0), c(0, 0, .05, 0), c(0, 0, 0, .05)),
               K0 = 1)
  expect_identical(facts$pair_status, "available")
  expect_identical(facts$target_eligible_n, 0L)
  expect_identical(facts$source_ineligible_n, 0L)
  expect_identical(facts$n_pairs, 1L)              # backbone only
  expect_true(is.finite(facts$phi_min))
  expect_true(is.finite(facts$rmsd))
  expect_false(facts$collision)
  expect_identical(facts$unmatched_n, 2L)
  expect_identical(facts$.assign$target, c(1L, NA_integer_))
  expect_identical(facts$.assign$congruence[2L], NA_real_)
  expect_identical(facts$.assign$eligible_source, c(NA, TRUE))

  ## An ineligible source is dropped from P_kl and counted.
  scarce <- .pf(cbind(c(.9, .1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, .05)),
                cbind(c(.85, .1, 0, 0), c(0, 0, 1, 0)), K0 = 1)
  expect_identical(scarce$pair_status, "available")
  expect_identical(scarce$source_ineligible_n, 1L)
  expect_identical(scarce$n_pairs, 2L)
  expect_true(is.finite(scarce$phi_min))
  expect_identical(scarce$.assign$eligible_source, c(NA, TRUE, FALSE))
})

test_that("only a malformed loading endpoint makes a pair unavailable", {
  lam_from <- cbind(c(.9, .1, 0, 0), c(0, 0, .8, .2))
  lam_to <- cbind(c(.85, .05, 0, 0), c(0, 0, .75, .25), c(.4, 0, .4, 0))

  bad_source <- .pf(lam_from, lam_to, K0 = 1, ok_from = FALSE)
  expect_identical(bad_source$pair_status, "unavailable")
  expect_identical(bad_source$pair_reason, "malformed_loading")
  expect_identical(bad_source$source_ineligible_n, NA_integer_)
  expect_identical(bad_source$target_eligible_n, 2L)   # its own endpoint lives
  expect_identical(bad_source$n_pairs, NA_integer_)
  expect_identical(bad_source$collision, NA)
  expect_identical(bad_source$n_collisions, NA_integer_)
  expect_identical(bad_source$collision_targets, NA_character_)
  expect_identical(bad_source$unmatched_n, NA_integer_)
  expect_identical(bad_source$backbone_degenerate, NA)

  bad_target <- .pf(lam_from, lam_to, K0 = 1, ok_to = FALSE)
  expect_identical(bad_target$pair_status, "unavailable")
  expect_identical(bad_target$pair_reason, "malformed_loading")
  expect_identical(bad_target$source_ineligible_n, 0L)
  expect_identical(bad_target$target_eligible_n, NA_integer_)
  expect_identical(bad_target$n_pairs, NA_integer_)

  ## Nothing else reaches the status: every other fixture in this file is
  ## available, including collided, degenerate, empty and scarce ones.
  reasons <- vapply(
    list(
      .pf(lam_from, lam_to, K0 = 1),
      .pf(lam_from, lam_to, K0 = 1, pip_ok_from = FALSE),
      .pf(cbind(rep(0, 4), c(0, 0, 1, 0)), lam_to, K0 = 1),
      .pf(cbind(.unit(1, 4), .unit(2, 4), .unit(2, 4)),
          cbind(.unit(1, 4), .unit(2, 4), .unit(3, 4), .unit(4, 4)), K0 = 1)
    ),
    function(f) f$pair_status, character(1)
  )
  expect_identical(unique(reasons), "available")
})

test_that("the matcher never consults fit status", {
  ## The status computation may not read a convergence fact.  Deparsing the
  ## installed closures checks the shipped code even where R/ is unavailable.
  matcher <- c(".pefa_ssl", ".pefa_sign", ".pefa_congruence",
               ".pefa_aligned_rmsd", ".pefa_ari", ".pefa_pair_facts")
  code <- unlist(lapply(matcher, function(nm) {
    deparse(get(nm, envir = asNamespace("vbpm")))
  }), use.names = FALSE)
  expect_false(any(grepl("fit_status", code, fixed = TRUE)))
  expect_false(any(grepl("converged", code, fixed = TRUE)))
  expect_false(any(grepl("from_fit_ok", code, fixed = TRUE)))

  ## And the source file itself, whenever the tests run beside it.
  path <- test_path("..", "..", "R", "pefa_match.R")
  expect_false(file.exists(path) &&
                 any(grepl("fit_status", readLines(path), fixed = TRUE)))

  ## The matcher takes no convergence argument at all.
  expect_false(any(c("fit_ok_from", "fit_ok_to", "converged") %in%
                     names(formals(vbpm:::.pefa_pair_facts))))
})
