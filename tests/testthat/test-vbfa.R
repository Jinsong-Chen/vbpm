## Core smoke tests for vbpm, ported from the analysis/vbfa_r2 suite.

## small oblique DGP with two anchors per factor
make_dat <- function(N = 400, K = 3, seed = 1, ld_pairs = NULL, ld_rho = .5) {
  set.seed(seed)
  J   <- K * 8
  grp <- rep(1:K, each = 8)
  Phi <- matrix(.3, K, K); diag(Phi) <- 1
  Lam <- matrix(0, J, K); Lam[cbind(1:J, grp)] <- .65
  cross <- seq(3, J, by = 5); Lam[cbind(cross, (grp[cross] %% K) + 1)] <- .35
  psi_var <- pmax(1 - rowSums((Lam %*% chol(Phi))^2), .25)
  Sig_e <- diag(psi_var)
  if (!is.null(ld_pairs)) for (r in seq_len(nrow(ld_pairs))) {
    i <- ld_pairs[r, 1]; j <- ld_pairs[r, 2]
    Sig_e[i, j] <- Sig_e[j, i] <- ld_rho * sqrt(psi_var[i] * psi_var[j])
  }
  eta <- matrix(rnorm(N * K), N, K) %*% chol(Phi)
  E   <- matrix(rnorm(N * J), N, J) %*% chol(Sig_e)
  Q <- matrix(-1L, J, K)
  for (k in 1:K) { a <- which(grp == k)[1:2]; Q[a, ] <- 0L; Q[a, k] <- 1L }
  list(Y = eta %*% t(Lam) + E, Q = Q, Lam = Lam, grp = grp)
}

congruence <- function(Le, Lt) {
  K <- ncol(Lt); cc <- abs(cor(Le, Lt)); perm <- integer(K); used <- integer(0)
  for (k in order(-apply(cc, 2, max))) {
    cand <- order(-cc[, k]); cand <- cand[!(cand %in% used)]
    perm[k] <- cand[1]; used <- c(used, cand[1])
  }
  vapply(1:K, function(k) abs(sum(Le[, perm[k]] * Lt[, k])) /
           sqrt(sum(Le[, perm[k]]^2) * sum(Lt[, k]^2)), numeric(1))
}

test_that("scalar v0 = 0.001 reproduces the fixed-spike estimator deterministically", {
  d <- make_dat()
  f1 <- vbfa(d$Y, d$Q, v0 = 0.001, max_it = 1500, convChk = FALSE)
  f2 <- vbfa(d$Y, d$Q, v0 = 0.001, max_it = 1500, convChk = FALSE)
  expect_identical(f1$Lam, f2$Lam)            # same seed -> identical
  expect_equal(f1$iter, f2$iter)
})

test_that("default dynamic path converges and recovers loadings", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, convChk = FALSE)
  expect_equal(f$flag, 1)
  expect_true(is.finite(f$ELBO))
  expect_gt(min(congruence(f$Lam, d$Lam)), 0.85)
})

test_that("orthogonal bifactor returns Phi = I", {
  d  <- make_dat()
  Qb <- cbind(1L, d$Q)
  f  <- vbfa(d$Y, Qb, orthogonal = TRUE, convChk = FALSE)
  expect_equal(f$flag, 1)
  expect_identical(f$Phi, diag(1, ncol(Qb), ncol(Qb)))
})

test_that("LD mode recovers planted residual pairs; ELBO is NA", {
  pairs <- rbind(c(4, 9), c(11, 20), c(15, 23))
  d <- make_dat(seed = 2, ld_pairs = pairs, ld_rho = .5)
  Qe <- matrix(-1L, 24, 24); diag(Qe) <- 0L
  f <- vbfa(d$Y, d$Q, Qe = Qe, convChk = FALSE, tolVal = 1e-3, max_it = 300)
  expect_true(is.na(f$ELBO))
  Poff <- abs(f$Psi); Poff[lower.tri(Poff, diag = TRUE)] <- 0
  top3 <- order(-Poff)[1:3]
  found <- apply(t(apply(cbind(row(Poff)[top3], col(Poff)[top3]), 1, sort)), 1,
                 paste, collapse = "-")
  planted <- apply(pairs, 1, paste, collapse = "-")
  expect_gte(sum(planted %in% found), 2)
})

test_that("NA input is rejected", {
  d <- make_dat(); Y <- d$Y; Y[1, 1] <- NA
  expect_error(vbfa(Y, d$Q), "NA")
})

test_that("vbfa_fit returns the documented statistics and rejects LD fits", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, convChk = FALSE)
  fs <- vbfa_fit(f, d$Y, d$Q)
  expect_named(fs, c("t_nom", "t", "t_S", "RMSEA", "SRMR", "CFI", "TLI",
                     "AIC", "BIC", "AIC_S", "BIC_S", "ELBO"))
  Qe <- matrix(-1L, 24, 24); diag(Qe) <- 0L
  fl <- vbfa(d$Y, d$Q, Qe = Qe, convChk = FALSE, max_it = 100)
  expect_error(vbfa_fit(fl, d$Y, d$Q), "LD")
})

test_that("pefa_vb selects a factor number within the window", {
  d <- make_dat()
  Q0 <- d$Q[, 1:2, drop = FALSE]
  r <- pefa_vb(Q0, d$Y, Kmin = 2, Kmax = 4, verbose = FALSE,
               max_it = 1500, v0 = 0.001)
  expect_true(r$selected_K %in% 2:4)
  expect_true(all(r$sweep$converged))
})

test_that("sim_fa produces data with planted local dependence", {
  s <- sim_fa(N = 300, K = 3, ipf = 6, lam = .7, lac = .3, ecr = .3, rseed = 7)
  expect_equal(dim(s$dat), c(300L, 18L))
  expect_false(anyNA(s$dat))
  expect_gt(nrow(s$ofd_ind), 0)
})
