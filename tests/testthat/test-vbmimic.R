## Tests for the MIMIC estimator and the sim_lvm generator.

## MIMIC design used throughout: K factors, each predicted by 3 of the P
## covariates, in disjoint blocks.
mimic_B <- function(K = 3, P = 9, val = .3) {
  B <- matrix(0, K, P)
  for (k in seq_len(K)) B[k, ((k - 1) * 3 + 1):(k * 3)] <- val
  B
}

## two anchors per factor -- the partially confirmatory idiom shared with vbfa()
anchor_QA <- function(J = 18, K = 3) {
  Q <- matrix(-1L, J, K); grp <- rep(seq_len(K), each = J / K)
  for (k in seq_len(K)) { a <- which(grp == k)[1:2]; Q[a, ] <- 0L; Q[a, k] <- 1L }
  Q
}

test_that("vbmimic recovers the structural pattern up to factor permutation", {
  B <- mimic_B()
  s <- sim_lvm(N = 600, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 11)
  Y <- s$dat[, 1:18]; X <- s$dat[, 19:27]
  f <- vbmimic(Y, X, anchor_QA(), matrix(-1L, 3, 9))

  expect_equal(f$flag, 1)
  expect_true(is.na(f$ELBO))                     # documented: no ELBO yet
  expect_equal(dim(f$A), c(18L, 3L))
  expect_equal(dim(f$B), c(3L, 9L))

  ## same covariate blocks, allowing for label switching
  truth <- apply(B != 0, 1, function(r) paste(which(r), collapse = ","))
  est   <- apply(abs(f$B) > .15, 1, function(r) paste(which(r), collapse = ","))
  expect_setequal(est, truth)
})

test_that("a fully exploratory Q_A AND Q_B still converges", {
  ## With nothing anchored in either part the solution is rotationally
  ## ambiguous: convergence is guaranteed to be clean, recovery is not
  ## (observed both recovering and merging factors, depending on
  ## standardization and seed). So only convergence is asserted here; the
  ## anchored test above carries the recovery claim, and the docs carry the
  ## "anchor at least one part" advice.
  B <- mimic_B()
  s <- sim_lvm(N = 600, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 11)
  f <- vbmimic(s$dat[, 1:18], s$dat[, 19:27],
               matrix(-1L, 18, 3), matrix(-1L, 3, 9))
  expect_equal(f$flag, 1)
})

test_that("vbmimic is deterministic (consumes no RNG)", {
  s <- sim_lvm(N = 300, K = 3, J = 18, P = 9, b = mimic_B(), phx = 0, rseed = 5)
  Y <- s$dat[, 1:18]; X <- s$dat[, 19:27]
  QA <- matrix(-1L, 18, 3); QB <- matrix(-1L, 3, 9)
  set.seed(1);   f1 <- vbmimic(Y, X, QA, QB)
  set.seed(999); f2 <- vbmimic(Y, X, QA, QB)
  expect_identical(f1$A, f2$A)
  expect_identical(f1$B, f2$B)
  expect_equal(f1$iter, f2$iter)
})

test_that("vbmimic validates its inputs", {
  s <- sim_lvm(N = 200, K = 3, J = 18, P = 9, b = mimic_B(), phx = 0, rseed = 2)
  Y <- s$dat[, 1:18]; X <- s$dat[, 19:27]
  QA <- matrix(-1L, 18, 3); QB <- matrix(-1L, 3, 9)

  Yna <- Y; Yna[1, 1] <- NA
  expect_error(vbmimic(Yna, X, QA, QB), "NA")
  Xna <- X; Xna[1, 1] <- NA
  expect_error(vbmimic(Y, Xna, QA, QB), "NA")
  expect_error(vbmimic(Y, X[1:10, ], QA, QB), "same number of rows")
  expect_error(vbmimic(Y, X, matrix(-1L, 17, 3), QB), "nrow\\(Q_A\\)")
  expect_error(vbmimic(Y, X, QA, matrix(-1L, 3, 8)), "ncol\\(Q_B\\)")
  expect_error(vbmimic(Y, X, QA, matrix(-1L, 2, 9)), "nrow\\(Q_B\\)")
  expect_error(vbmimic(Y, X, matrix(2L, 18, 3), QB), "Q_A entries")
  expect_error(vbmimic(Y, X, QA, matrix(2L, 3, 9)), "Q_B entries")
  expect_error(vbmimic(Y, X, QA, QB, v0 = -1), "v0 must be positive")
})

test_that("anchored entries of Q_B are estimated and fixed zeros stay zero", {
  B <- mimic_B()
  s <- sim_lvm(N = 500, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 4)
  QB <- matrix(0L, 3, 9)          # everything fixed to zero except one anchor
  QB[1, 1] <- 1L
  f <- vbmimic(s$dat[, 1:18], s$dat[, 19:27], matrix(-1L, 18, 3), QB)
  expect_true(all(f$B[QB == 0] == 0))
  expect_true(f$B[1, 1] != 0)
})

test_that("the v0 path runs and a scalar v0 is a single stage", {
  s <- sim_lvm(N = 300, K = 3, J = 18, P = 9, b = mimic_B(), phx = 0, rseed = 6)
  Y <- s$dat[, 1:18]; X <- s$dat[, 19:27]
  QA <- matrix(-1L, 18, 3); QB <- matrix(-1L, 3, 9)
  f_sc <- vbmimic(Y, X, QA, QB, v0 = 0.001)
  f_pa <- vbmimic(Y, X, QA, QB, v0 = c(0.01, 0.005, 0.001))
  expect_equal(f_sc$path$n_stage, 1L)
  expect_equal(f_pa$path$n_stage, 3L)
  expect_length(f_pa$path$stage_iters, 3)
  expect_equal(f_pa$flag, 1)
})

test_that("vb_fit refuses vbmimic fits", {
  s <- sim_lvm(N = 300, K = 3, J = 18, P = 9, b = mimic_B(), phx = 0, rseed = 8)
  f <- vbmimic(s$dat[, 1:18], s$dat[, 19:27],
              matrix(-1L, 18, 3), matrix(-1L, 3, 9))
  expect_error(vb_fit(f, s$dat[, 1:18], matrix(-1L, 18, 3)), "vbmimic")
})

## ---- sim_lvm ------------------------------------------------------------

test_that("sim_lvm confines missingness to the item block", {
  s <- sim_lvm(N = 200, K = 3, J = 18, P = 5, misp = .3, rseed = 9)
  expect_equal(ncol(s$dat), 23L)
  expect_gt(sum(is.na(s$dat[, 1:18])), 0)
  expect_equal(sum(is.na(s$dat[, 19:23])), 0)   # regression: LAWBL recycled here
})

test_that("phd correlates the structural disturbances", {
  B <- mimic_B(K = 5, P = 15)
  s0 <- sim_lvm(N = 800, K = 5, J = 30, P = 15, b = B, phx = 0, phd = 0, rseed = 1)
  s3 <- sim_lvm(N = 800, K = 5, J = 30, P = 15, b = B, phx = 0, phd = .3, rseed = 1)
  off <- function(m) mean(m[upper.tri(m)])
  expect_lt(abs(off(s0$PHI)), .10)
  expect_gt(off(s3$PHI), .10)
})

test_that("sim_lvm supports latent predictors and rejects bad coefficients", {
  s <- sim_lvm(N = 300, K = 4, J = 20, K1 = 2, b1 = .4, rseed = 3)
  expect_equal(dim(s$mb1), c(2L, 2L))
  ## dense b over-explains the factors -> negative disturbance variance
  expect_error(sim_lvm(N = 100, K = 5, J = 30, P = 15, b = .3),
               "disturbance variance")
  expect_error(sim_lvm(N = 100, K = 3, J = 18, K1 = 3), "K1")
})

test_that("sim_lvm's plain factor path is unchanged and still simulates LD", {
  s <- sim_lvm(N = 300, K = 3, J = 18, ecr = .3, rseed = 7)
  expect_equal(dim(s$dat), c(300L, 18L))
  expect_false(anyNA(s$dat))
  expect_gt(nrow(s$ofd_ind), 0)
})
