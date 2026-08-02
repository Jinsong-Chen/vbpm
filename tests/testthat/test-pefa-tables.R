## 0.8.0: the three-table pefa() object, bottleneck matching, stability
## indexes, the skip-gap selection policy, and the mode/convergence gates.

mk <- function(N = 300, K = 3, seed = 4) {
  set.seed(seed)
  J <- K * 6; grp <- rep(1:K, each = 6)
  Lam <- matrix(0, J, K); Lam[cbind(1:J, grp)] <- .7
  eta <- matrix(rnorm(N * K), N, K)
  Y <- eta %*% t(Lam) + matrix(rnorm(N * J, sd = sqrt(1 - .49)), N, J)
  Q0 <- matrix(-1L, J, 2)
  for (k in 1:2) { a <- which(grp == k)[1:2]; Q0[a, ] <- 0L; Q0[a, k] <- 1L }
  list(Y = Y, Q0 = Q0, grp = grp)
}

## ---- matching ------------------------------------------------------------

test_that("column matching is permutation and sign invariant", {
  set.seed(1)
  A <- matrix(rnorm(30), 10, 3)
  B <- cbind(A[, c(2, 3, 1)] + matrix(rnorm(30, sd = .01), 10, 3), rnorm(10))
  base <- vbpm:::.match_columns(A, B)
  for (p in list(c(2, 1, 3), c(3, 2, 1))) {
    q <- vbpm:::.match_columns(A[, p, drop = FALSE], B)
    expect_equal(sort(q$phi), sort(base$phi))          # same multiset
    expect_equal(min(q$phi), min(base$phi))            # same binding value
  }
  ## flipping a column's sign cannot change congruence
  A2 <- A; A2[, 2] <- -A2[, 2]
  expect_equal(sort(vbpm:::.match_columns(A2, B)$phi), sort(base$phi))
})

test_that("bottleneck matching beats greedy on the counterexample", {
  ## the plan's [[.95,.90],[.93,.10]] pattern, realized as loading vectors
  a1 <- c(1, 0, 0); a2 <- c(0, 1, 0)
  b1 <- c(1, .2, 0); b2 <- c(.9, 0, 0); b3 <- c(0, 0, 1)
  A <- cbind(a1, a2); B <- cbind(b1, b2, b3)
  m <- vbpm:::.match_columns(A, B)
  ## greedy in stored order would give a1 -> b1 and strand a2; max-min must not
  expect_gt(min(m$phi), 0)
  ## exhaustive optimum for comparison
  best <- max(vapply(list(c(1, 2), c(1, 3), c(2, 1), c(2, 3), c(3, 1), c(3, 2)),
                     function(ix) min(abs(vbpm:::.tucker(A[, 1], B[, ix[1]])),
                                      abs(vbpm:::.tucker(A[, 2], B[, ix[2]]))),
                     numeric(1)))
  expect_equal(min(m$phi), best)
})

test_that("empty columns match at zero congruence rather than being skipped", {
  A <- cbind(c(1, 0, 0), c(0, 0, 0))
  B <- cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  m <- vbpm:::.match_columns(A, B)
  expect_equal(min(m$phi), 0)            # not optimistic
  expect_true(is.na(m$sign[2]))          # sign of an empty column is not identified
  expect_identical(m$empty_from, 1L)
})

## ---- indexes -------------------------------------------------------------

test_that("ARI degenerates to NA, never NaN, and singletons cannot agree", {
  expect_true(is.na(vbpm:::.ari(rep(1, 5), rep(1, 5))))       # one class
  expect_true(is.na(vbpm:::.ari(1:5, 1:5)))                   # all singletons
  ## unassigned/ambiguous items become distinct singletons
  L <- rbind(c(.8, 0), c(0, .8), c(0, 0), c(.5, .5))
  p <- vbpm:::.partition_dominant(L)
  expect_identical(p$n_unassigned, 1L)
  expect_identical(p$n_ambiguous, 1L)
  expect_identical(anyDuplicated(p$labels[3:4]), 0L)          # not pooled
})

## ---- object shape --------------------------------------------------------

test_that("the sweep object has three tables indexed three ways", {
  d <- mk()
  p <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 4, v0 = .001, max_it = 300,
            verbose = FALSE)
  expect_false(any(grepl("_gain", names(p$sweep))))
  expect_identical(nrow(p$transitions), nrow(p$sweep) - 1L)
  expect_identical(p$transitions$K_from, 2:3)
  expect_identical(p$transitions$K_to, 3:4)
  expect_true(all(c("phi_min", "phi_mean", "n_phi90", "ari", "rmsd",
                    "surplus_max", "n_unassigned_from", "n_ambiguous_to",
                    "stability_status") %in% names(p$transitions)))
  ## gain orientation: larger is better for all three
  expect_equal(p$transitions$BIC_gain,
               p$sweep$BIC[-nrow(p$sweep)] - p$sweep$BIC[-1])
  ## no stored percentage anywhere
  expect_false(any(grepl("pct", c(names(p$sweep), names(p$transitions)))))
  expect_identical(p$object_schema_version, 2L)
})

test_that("a one-candidate window gives a typed zero-row transition table", {
  d <- mk()
  p <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 2, v0 = .001, max_it = 300,
            verbose = FALSE)
  expect_s3_class(p$transitions, "data.frame")
  expect_identical(nrow(p$transitions), 0L)
  expect_true(all(c("K_from", "phi_min", "stability_status") %in%
                  names(p$transitions)))
})

## ---- gates ---------------------------------------------------------------

test_that("orthogonal sweeps report gains but no loading-based stability", {
  d <- mk()
  p <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 3, bifactor = TRUE, v0 = .001,
            max_it = 300, verbose = FALSE)
  tr <- p$transitions
  expect_true(all(is.finite(tr$Objective_gain)))          # gains valid
  expect_true(all(is.na(tr$phi_min)))                     # stability gated
  expect_true(all(is.na(tr$ari)))
  expect_true(all(is.na(tr$n_unassigned_from)))
  expect_match(tr$stability_status[1], "rotational indeterminacy")
  expect_null(transition_detail(p, 2))                    # accessor respects it
})

test_that("keep_fits governs candidate slots only, and gates the accessor", {
  d <- mk()
  p <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 3, v0 = .001, max_it = 300,
            verbose = FALSE, keep_fits = "none")
  expect_true(all(vapply(p$fits, is.null, logical(1))))
  expect_s3_class(p$selected_fit, "vbpm_fit")
  ## every method still works with no candidate fits at all
  expect_output(print(p), "PEFA sweep")
  expect_output(print(summary(p)), "primary cut")
  ## stability was computed in-sweep, so the aggregate row survives
  expect_true(is.finite(p$transitions$phi_min[1]))
})

## ---- selection -----------------------------------------------------------

test_that("the selection table covers criterion x form x cut x scope", {
  d <- mk()
  p <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 4, v0 = .001, max_it = 300,
            verbose = FALSE, cuts = c(primary = 20, sensitivity = 10))
  s <- p$selection
  expect_identical(nrow(s), 3L * 2L * (1L + 2L))     # 3 criteria, 2 scopes, raw + 2 cuts
  expect_true(all(is.na(s$cut[s$form == "raw"])))
  expect_true(all(is.na(s$gain_max[s$form == "raw"])))
  expect_setequal(unique(s$eligibility_scope), c("converged+fit", "converged"))
  ## distinct cut names survive even at equal percentages
  p2 <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 3, v0 = .001, max_it = 300,
             verbose = FALSE, cuts = c(primary = 10, mirror = 10))
  expect_setequal(unique(p2$selection$cut[p2$selection$form == "gain"]),
                  c("primary", "mirror"))
})

test_that("cuts validation and the deprecated delta path", {
  d <- mk()
  expect_error(pefa(d$Q0, d$Y, 2, 3, cuts = c(primary = 10), delta = 10,
                    verbose = FALSE), "not both")
  expect_error(pefa(d$Q0, d$Y, 2, 3, cuts = c(a = 10), verbose = FALSE),
               "exactly one entry named 'primary'")
  expect_error(pefa(d$Q0, d$Y, 2, 3, cuts = c(primary = 200), verbose = FALSE),
               "\\[0, 100\\]")
  p <- pefa(d$Q0, d$Y, 2, 3, delta = 25, v0 = .001, max_it = 300,
            verbose = FALSE)
  expect_identical(p$cuts, c(primary = 25))
  expect_identical(p$delta, 25)                     # deprecated alias
  expect_identical(p$boundary, p$window_boundary)   # deprecated alias
})

test_that("an eligibility gap breaks the chain instead of being bridged", {
  d <- mk()
  p <- pefa(d$Q0, d$Y, Kmin = 2, Kmax = 5, v0 = .001, max_it = 300,
            verbose = FALSE)
  ## force an interior candidate to be ineligible, then re-run selection
  sw <- p$sweep; sw$converged[sw$K == 3] <- FALSE
  sw$eligible <- sw$converged & sw$pass_fit
  sel <- vbpm:::.pefa_selection(sw, p$transitions, c(primary = 10), 2,
                                p$objective_type)
  gainrows <- sel[sel$form == "gain" & sel$eligibility_scope == "converged+fit", ]
  ## no rule may use an edge touching K = 3
  usable <- with(p$transitions, K_from != 3 & K_to != 3)
  expect_true(all(gainrows$gain_max %in%
                  c(NA_real_, max(p$transitions$Objective_gain[usable]),
                    max(p$transitions$BIC_gain[usable]),
                    max(p$transitions$ELBO_gain[usable]))))
  expect_false(any(is.na(gainrows$selected_K)))
})
