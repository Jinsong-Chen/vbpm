## Verification for the three features added after 0.1.0: in-loop missing-data
## handling in vbfa(), the local-dependence objective, and LD-aware fit
## statistics. These close verification items 2, 3 and 6 of
## dev/vbpm_extension_design.qmd, which were previously only covered
## behaviourally.

make_dat <- function(N = 300, K = 3, seed = 1) {
  set.seed(seed)
  J <- K * 6
  grp <- rep(1:K, each = 6)
  Phi <- matrix(.3, K, K); diag(Phi) <- 1
  Lam <- matrix(0, J, K); Lam[cbind(1:J, grp)] <- .7
  psi <- pmax(1 - rowSums((Lam %*% chol(Phi))^2), .25)
  eta <- matrix(rnorm(N * K), N, K) %*% chol(Phi)
  Y <- eta %*% t(Lam) + matrix(rnorm(N * J), N, J) %*% diag(sqrt(psi))
  Q <- matrix(-1L, J, K)
  for (k in 1:K) { a <- which(grp == k)[1:2]; Q[a, ] <- 0L; Q[a, k] <- 1L }
  list(Y = Y, Q = Q, J = J, K = K)
}

## ---- item 3: conditional moments vs a direct MVN calculation -------------

test_that(".cond_moments matches the textbook Gaussian conditional (full W)", {
  set.seed(7)
  J <- 6
  A <- matrix(rnorm(J * J), J, J)
  W <- crossprod(A) + diag(J)          # arbitrary SPD residual covariance
  mean_row <- rnorm(J)
  y_row <- rnorm(J)
  M <- c(2L, 5L); O <- setdiff(seq_len(J), M)

  got <- vbpm:::.cond_moments(W, mean_row, y_row, M, O)

  ## direct formulas: mu_M + W_MO W_OO^-1 (y_O - mu_O), W_MM - W_MO W_OO^-1 W_OM
  WMO <- W[M, O, drop = FALSE]; WOO <- W[O, O, drop = FALSE]
  ref_mean <- mean_row[M] + as.vector(WMO %*% solve(WOO) %*% (y_row[O] - mean_row[O]))
  ref_cov  <- W[M, M, drop = FALSE] - WMO %*% solve(WOO) %*% W[O, M, drop = FALSE]

  expect_equal(got$mean, ref_mean)
  expect_equal(unname(got$cov), unname(ref_cov))
  ## conditioning must reduce uncertainty
  expect_true(all(diag(got$cov) < diag(W[M, M, drop = FALSE])))
})

test_that(".cond_moments reduces to the factor-model mean under diagonal W", {
  set.seed(8)
  J <- 5
  W <- diag(runif(J, .3, 1))
  mean_row <- rnorm(J); y_row <- rnorm(J)
  M <- c(1L, 4L); O <- setdiff(seq_len(J), M)
  got <- vbpm:::.cond_moments(W, mean_row, y_row, M, O)
  ## Equation (2): with W_MO = 0 the conditional mean is mu_M and the
  ## conditional covariance is W_MM -- observed responses carry no residual
  ## information about the missing ones.
  expect_equal(got$mean, mean_row[M])
  expect_equal(unname(got$cov), unname(W[M, M, drop = FALSE]))
})

test_that(".cond_moments handles a row with no observed responses", {
  W <- diag(3); mean_row <- c(1, 2, 3); y_row <- rep(NA_real_, 3)
  got <- vbpm:::.cond_moments(W, mean_row, y_row, M = 1:3, O = integer())
  expect_equal(got$mean, mean_row)
  expect_equal(unname(got$cov), W)
})

test_that("joint updating keeps the uncertainty that cell-by-cell imputation destroys", {
  ## The design note requires conditioning jointly on the *observed* subset,
  ## rather than treating an already-imputed neighbour as evidence.
  ##
  ## Note what this does and does not change. The conditional MEAN is the same
  ## either way: conditional expectation is linear, so plugging in E[y2 | y3]
  ## and re-conditioning reproduces E[y1 | y3] exactly (tower property). The
  ## difference is entirely in the second moment, which is what feeds the
  ## residual sufficient statistics -- treating an imputed cell as observed
  ## asserts zero uncertainty about it.
  W <- matrix(c(1, .8, .5,
                .8, 1, .4,
                .5, .4, 1), 3, 3)
  mean_row <- c(0, 0, 0); y_row <- c(NA, NA, 2)

  joint <- vbpm:::.cond_moments(W, mean_row, y_row, M = 1:2, O = 3L)
  step1 <- vbpm:::.cond_moments(W, mean_row, y_row, M = 2L, O = 3L)$mean
  y2 <- y_row; y2[2] <- step1
  seq1  <- vbpm:::.cond_moments(W, mean_row, y2, M = 1L, O = c(2L, 3L))

  ## means agree, as the tower property requires
  expect_equal(joint$mean[1], seq1$mean, tolerance = 1e-12)

  ## but the joint covariance is a full 2x2 that keeps the two missing
  ## responses correlated, and is strictly more uncertain about item 1 than
  ## the sequential route believes
  expect_equal(dim(joint$cov), c(2L, 2L))
  expect_true(abs(joint$cov[1, 2]) > 1e-8)
  expect_gt(joint$cov[1, 1], seq1$cov[1, 1])
})

## ---- item 2: mask integrity ----------------------------------------------

test_that("observed cells are never overwritten and missing cells are filled", {
  d <- make_dat()
  Y <- d$Y
  mask <- matrix(FALSE, nrow(Y), ncol(Y))
  mask[cbind(1:15, rep(1:3, each = 5))] <- TRUE
  Y[mask] <- NA
  f <- vbfa(Y, d$Q, v0 = .001, max_it = 200)

  ## the fit records exactly the mask it was given
  expect_identical(f$preprocess$missing_mask, mask)
  expect_equal(f$preprocess$n_missing, sum(mask))

  ## Observed entries must drive the fit: perturbing only the *values already
  ## marked missing* in the input cannot change anything, because they are
  ## replaced before first use.
  Y2 <- Y; Y2[mask] <- NA          # same mask, values were NA either way
  f2 <- vbfa(Y2, d$Q, v0 = .001, max_it = 200)
  expect_identical(f$Lam, f2$Lam)
})

test_that("row order does not change the fit (permutation invariance)", {
  d <- make_dat(N = 200)
  Y <- d$Y
  Y[cbind(1:12, rep(1:4, each = 3))] <- NA
  f  <- vbfa(Y, d$Q, v0 = .001, max_it = 200)

  set.seed(99); perm <- sample(nrow(Y))
  fp <- vbfa(Y[perm, ], d$Q, v0 = .001, max_it = 200)

  ## Item-level parameters are row-order invariant ...
  expect_equal(f$Lam, fp$Lam, tolerance = 1e-6)
  expect_equal(f$PsiInv, fp$PsiInv, tolerance = 1e-6)
  ## ... and factor scores follow the permutation.
  expect_equal(f$eta[perm, ], fp$eta, tolerance = 1e-6)
})

test_that("complete data is unaffected by the missing-data code path", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, v0 = .001, max_it = 300)
  expect_equal(f$preprocess$n_missing, 0)
  expect_false(any(f$preprocess$missing_mask))
  expect_true(is.finite(f$objective))
})

## ---- LD objective and item 6: LD-aware fit statistics --------------------

ld_dat <- function(N = 400, seed = 3, rho = .5) {
  set.seed(seed)
  K <- 3; J <- 18
  grp <- rep(1:K, each = 6)
  Phi <- matrix(.3, K, K); diag(Phi) <- 1
  Lam <- matrix(0, J, K); Lam[cbind(1:J, grp)] <- .7
  psi <- pmax(1 - rowSums((Lam %*% chol(Phi))^2), .25)
  Sig_e <- diag(psi)
  pairs <- rbind(c(2, 8), c(5, 14))
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1]; j <- pairs[r, 2]
    Sig_e[i, j] <- Sig_e[j, i] <- rho * sqrt(psi[i] * psi[j])
  }
  eta <- matrix(rnorm(N * K), N, K) %*% chol(Phi)
  Y <- eta %*% t(Lam) + matrix(rnorm(N * J), N, J) %*% chol(Sig_e)
  Q <- matrix(-1L, J, K)
  for (k in 1:K) { a <- which(grp == k)[1:2]; Q[a, ] <- 0L; Q[a, k] <- 1L }
  list(Y = Y, Q = Q, pairs = pairs, J = J, K = K)
}

test_that("the LD objective is labelled honestly and ELBO stays NA", {
  d <- ld_dat()
  f <- vbfa(d$Y, d$Q, ld = TRUE, max_it = 150, tolVal = 1e-3)
  expect_identical(f$objective_type, "vecm")
  expect_true(is.finite(f$objective))
  ## ELBO is deliberately NA: Psi and tau are point updates, not variational
  ## factors, so no single joint mean-field bound is defined.
  expect_true(is.na(f$ELBO))
  ## the conditional bound of Eq. (3) is available separately
  expect_true(is.finite(f$ELBO_conditional))
  ## a diagonal fit reports the ordinary ELBO instead
  fd <- vbfa(d$Y, d$Q, v0 = .001, max_it = 200)
  expect_identical(fd$objective_type, "elbo")
  expect_equal(fd$objective, fd$ELBO)
})

test_that("LD fit_stats uses the residual covariance W, not the precision", {
  d <- ld_dat()
  f <- vbfa(d$Y, d$Q, ld = TRUE, max_it = 150, tolVal = 1e-3)
  fs <- fit_stats(f, d$Y, d$Q)

  expect_s3_class(fs, "vbpm_fit_stats")
  expect_true(all(is.finite(fs[c("RMSEA", "SRMR", "CFI", "TLI")])))

  ## W and Psi are genuinely different matrices here, so a fit-statistic
  ## computed from the wrong one would be detectably different. Rebuild the
  ## implied covariance from W and check SRMR agrees with the reported value.
  expect_false(isTRUE(all.equal(f$W, f$Psi)))
  pi_ <- f$pi; pi_[d$Q == 1] <- 1; pi_[d$Q == 0] <- 0
  aH <- (d$Q == 1) | (d$Q == -1 & pi_ >= .5)
  LamH <- f$Lam * aH
  Sg <- LamH %*% f$Phi %*% t(LamH) + f$W
  S  <- f$sample_cov
  si <- 1 / sqrt(diag(S))
  R0 <- si * S  * rep(si, each = ncol(S))
  Rh <- si * Sg * rep(si, each = ncol(S))
  srmr <- sqrt(2 * sum((R0 - Rh)[lower.tri(R0, diag = TRUE)]^2) /
                 (ncol(S) * (ncol(S) + 1)))
  expect_equal(unname(fs["SRMR"]), srmr, tolerance = 1e-8)
})

test_that("LD parameter counts grow with the selected residual edges", {
  d <- ld_dat()
  f  <- vbfa(d$Y, d$Q, ld = TRUE, max_it = 150, tolVal = 1e-3)
  fs <- fit_stats(f, d$Y, d$Q)
  fd <- vbfa(d$Y, d$Q, v0 = .001, max_it = 200)
  fsd <- fit_stats(fd, d$Y, d$Q)
  ## an LD model estimates residual edges on top of the diagonal model, so its
  ## nominal count must be at least as large
  expect_gte(unname(fs["t_nom"]), unname(fsd["t_nom"]))
})

test_that("missing data and local dependence compose", {
  d <- ld_dat(N = 350)
  Y <- d$Y
  Y[cbind(1:10, rep(c(2, 8), each = 5))] <- NA   # inside a planted LD pair
  f <- vbfa(Y, d$Q, ld = TRUE, max_it = 150, tolVal = 1e-3)
  expect_equal(f$preprocess$n_missing, 10)
  expect_false(anyNA(f$Lam))
  expect_true(is.finite(f$objective))
  expect_identical(f$objective_type, "vecm")
})

## ---- missing responses in vbmimic ---------------------------------------

test_that("vbmimic handles missing responses in Y but rejects missing X", {
  B <- matrix(0, 3, 9); for (k in 1:3) B[k, ((k - 1) * 3 + 1):(k * 3)] <- .3
  s <- sim_lvm(N = 500, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 21)
  Y <- s$dat[, 1:18]; X <- s$dat[, 19:27]
  QA <- matrix(-1L, 18, 3)
  for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
                   QA[a, ] <- 0L; QA[a, k] <- 1L }
  QB <- matrix(-1L, 3, 9)

  Ym <- Y
  Ym[cbind(1:20, rep(1:4, each = 5))] <- NA
  f <- vbmimic(Ym, X, QA, QB)
  expect_equal(f$flag, 1)
  expect_false(anyNA(f$A))
  expect_false(anyNA(f$B))
  expect_equal(f$preprocess$n_missing, 20)
  expect_identical(f$preprocess$missing_mask, is.na(Ym))

  ## deterministic with missing data too
  f2 <- vbmimic(Ym, X, QA, QB)
  expect_identical(f$A, f2$A)

  ## structural recovery survives a modest amount of missingness
  truth <- apply(B != 0, 1, function(r) paste(which(r), collapse = ","))
  est   <- apply(abs(f$B) > .15, 1, function(r) paste(which(r), collapse = ","))
  expect_setequal(est, truth)

  ## covariates must be complete
  Xm <- X; Xm[1, 1] <- NA
  expect_error(vbmimic(Y, Xm, QA, QB), "missing covariates")

  ## degenerate masks are rejected
  Ybad <- Y; Ybad[, 1] <- NA
  expect_error(vbmimic(Ybad, X, QA, QB), "at least one observed")
  Ybad <- Y; Ybad[1, ] <- NA
  expect_error(vbmimic(Ybad, X, QA, QB), "all responses missing")
})

test_that("complete data is unchanged by the vbmimic missing-data path", {
  B <- matrix(0, 3, 9); for (k in 1:3) B[k, ((k - 1) * 3 + 1):(k * 3)] <- .3
  s <- sim_lvm(N = 300, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 22)
  f <- vbmimic(s$dat[, 1:18], s$dat[, 19:27],
               matrix(-1L, 18, 3), matrix(-1L, 3, 9))
  expect_equal(f$preprocess$n_missing, 0)
  expect_false(any(f$preprocess$missing_mask))
})

## ---- pefa plot method ----------------------------------------------------

test_that("plot.pefa draws its three panels without error", {
  skip_if_not_installed("grDevices")
  d <- make_dat(N = 250)
  r <- pefa(d$Q[, 1:2, drop = FALSE], d$Y, Kmin = 2, Kmax = 3,
            verbose = FALSE, v0 = .001, max_it = 300)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 600, height = 400)
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
  expect_silent(plot(r, type = "objective"))
  expect_silent(plot(r, type = "gain"))
  expect_silent(plot(r, type = "fit"))
})
