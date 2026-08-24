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

test_that("vbfa(v0 = 0.001) reproduces the frozen r1.1 estimator", {
  ## THE regression test protecting the paper-1 artifact. The fixture stores
  ## both the data and the reference outputs of vbfa_r1.1.R (CS repo), so this
  ## cannot be invalidated by generator drift; regenerate only via
  ## dev/make_r11_fixture.R and only if the reference itself is re-frozen.
  ##
  ## Numerical agreement is asserted at a tolerance rather than bitwise. The
  ## estimator runs through compiled BLAS/LAPACK, whose SIMD and FMA choices
  ## differ between platforms, so the last one or two ulps are not portable --
  ## a reference frozen under one BLAS will not reproduce exactly under
  ## another. The tolerance is still eight orders of magnitude tighter than
  ## any change of estimator would produce, so the guard keeps its power.
  ## Bitwise identity is checked separately below, where it is meaningful.
  ref <- readRDS(test_path("fixtures", "vbfa_r11_ref.rds"))
  f <- vbfa(ref$Y, ref$Q, v0 = 0.001, max_it = 1500)
  expect_equal(unname(f$Lam), unname(ref$Lam), tolerance = 1e-8)
  expect_equal(unname(f$Phi), unname(ref$Phi), tolerance = 1e-8)
  expect_equal(unname(f$PsiInv), unname(ref$PsiInv), tolerance = 1e-8)
  expect_equal(f$ELBO, ref$ELBO, tolerance = 1e-8)
  ## The iteration count is discrete and must match exactly: a different
  ## stopping point means the path itself changed, not just its rounding.
  expect_identical(f$iter, ref$iter)     # 79, matching CS smoke test T1
})


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

test_that("LD mode recovers planted residual pairs and exposes VECM objective", {
  pairs <- rbind(c(4, 9), c(11, 20), c(15, 23))
  d <- make_dat(seed = 2, ld_pairs = pairs, ld_rho = .5)
  Qe <- matrix(-1L, 24, 24); diag(Qe) <- 0L
  f <- vbfa(d$Y, d$Q, ld = TRUE, Qe = Qe, convChk = FALSE,
            tolVal = 1e-3, max_it = 300)
  expect_true(is.na(f$ELBO))
  expect_identical(f$objective_type, "vecm")
  expect_true(is.finite(f$objective))
  expect_equal(sum(f$objective_terms), f$objective, tolerance = 1e-10)
  Poff <- abs(f$Psi); Poff[lower.tri(Poff, diag = TRUE)] <- 0
  top3 <- order(-Poff)[1:3]
  found <- apply(t(apply(cbind(row(Poff)[top3], col(Poff)[top3]), 1, sort)), 1,
                 paste, collapse = "-")
  planted <- apply(pairs, 1, paste, collapse = "-")
  expect_gte(sum(planted %in% found), 2)

  ## Sparsity, not just ranking. Measured on this scenario the package selects
  ## exactly the 3 planted edges while the V4 reference ld-vb_r.R selects 24.
  ## The cause is upstream of the LD machinery, which is identical in both:
  ## the reference attenuates the loadings (mean primary .475 against a true
  ## .65), so more common variance is left in the residual covariance and QUIC
  ## reads it as local dependence. The packaged branch takes its measurement
  ## updates from r1.1 instead, recovers the loadings (.640), and leaves little
  ## for QUIC to find. Guarding the false-positive rate is therefore the
  ## meaningful regression, not numerical agreement with the reference.
  expect_lte(sum(f$q_star[upper.tri(f$q_star)] >= .5), 6)
})

test_that("incomplete continuous data use deterministic joint moments", {
  d <- make_dat(N = 180); Y <- d$Y
  Y[cbind(1:20, rep(1:4, each = 5))] <- NA
  Y[31, 1:3] <- NA
  f1 <- vbfa(Y, d$Q, v0 = .001, max_it = 300, convChk = FALSE)
  f2 <- vbfa(Y, d$Q, v0 = .001, max_it = 300, convChk = FALSE)
  expect_false(anyNA(f1$Lam))
  expect_identical(f1$Lam, f2$Lam)
  expect_equal(f1$preprocess$n_missing, 23)
  expect_identical(f1$preprocess$missing_mask, is.na(Y))
  expect_true(is.finite(f1$objective))
  Y_bad <- Y; Y_bad[, 1] <- NA
  expect_error(vbfa(Y_bad, d$Q), "at least one observed")
  Y_bad <- Y; Y_bad[1, ] <- NA
  expect_error(vbfa(Y_bad, d$Q), "all responses missing")
})

test_that("the ld switch governs Qe", {
  d <- make_dat()
  Qe <- matrix(-1L, 24, 24); diag(Qe) <- 0L

  ## Qe without ld = TRUE is ignored, with a warning, and the fit is diagonal
  expect_warning(f <- vbfa(d$Y, d$Q, Qe = Qe, convChk = FALSE), "ignored")
  expect_false(f$ld)
  expect_false(is.null(f$PsiInv))
  expect_true(is.finite(f$ELBO))

  ## ld = TRUE with Qe = NULL means fully exploratory local dependence
  f2 <- vbfa(d$Y, d$Q, ld = TRUE, convChk = FALSE,
             tolVal = 1e-3, max_it = 100)
  expect_true(f2$ld)
  expect_null(f2$PsiInv)
  expect_true(is.na(f2$ELBO))
  expect_equal(dim(f2$Psi), c(24L, 24L))
})

test_that("fit_stats is generic and supports LD fits", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, convChk = FALSE)
  fs <- fit_stats(f)
  expect_named(fs, c("t_nom", "t", "t_S", "RMSEA", "SRMR", "CFI", "TLI",
                     "AIC", "BIC", "AIC_S", "BIC_S", "ELBO", "objective"))
  ## for a diagonal fit the objective IS the ELBO; under LD it is the VECM
  ## objective, where ELBO is NA. The README promises this element.
  expect_equal(unname(fs["objective"]), unname(fs["ELBO"]))
  fl <- vbfa(d$Y, d$Q, ld = TRUE, convChk = FALSE, max_it = 300, tolVal = 1e-3)
  fsl <- fit_stats(fl)
  expect_true(all(is.finite(fsl[c("RMSEA", "SRMR", "CFI", "TLI", "BIC")])))
  expect_identical(attr(fsl, "objective_type"), "vecm")
})

test_that("fit_stats reads orthogonal from the fit (the bifactor footgun)", {
  d  <- make_dat()
  Qb <- cbind(1L, d$Q)
  fb <- vbfa(d$Y, Qb, orthogonal = TRUE, convChk = FALSE)

  ## the default call now counts the bifactor correctly: no free factor
  ## correlations. Before the fix it silently added K(K-1)/2 of them.
  fs_default  <- fit_stats(fb, Y = d$Y, Q = Qb)
  fs_explicit <- fit_stats(fb, Y = d$Y, Q = Qb, orthogonal = TRUE)
  expect_identical(fs_default, fs_explicit)
  K <- ncol(Qb)
  expect_equal(unname(fs_default["t_nom"]),
               sum((Qb == 1) | (Qb == -1 & {p <- fb$pi; p[Qb == 1] <- 1
                                            p[Qb == 0] <- 0; p >= .5})) + 24)

  ## a contradictory value errors instead of miscounting
  expect_error(fit_stats(fb, Y = d$Y, Q = Qb, orthogonal = FALSE), "contradicts")
  fo <- vbfa(d$Y, d$Q, convChk = FALSE)
  expect_error(fit_stats(fo, Y = d$Y, Q = d$Q, orthogonal = TRUE), "contradicts")
})

test_that("pefa returns evidence over the window and selects nothing", {
  d <- make_dat()
  Q0 <- d$Q[, 1:2, drop = FALSE]
  r <- pefa(Q0, d$Y, Kmin = 2, Kmax = 4, verbose = FALSE,
               max_it = 1500, v0 = 0.001)
  expect_s3_class(r, "pefa")
  expect_false(inherits(r, "vbpm_sweep"))
  expect_true(all(r$sweep$converged))

  ## The evidence-only contract: nothing here decides a count.
  expect_null(r$selected_K)
  expect_null(r$boundary)
  expect_false(any(grepl("_gain", names(r$sweep))))

  expect_identical(nrow(r$transitions), nrow(r$sweep) - 1L)
  ## 0.9.0 stores all-pairs evidence as three upper-triangular W x W matrices
  ## rather than one long pair table.
  W <- nrow(r$sweep)
  expect_identical(names(r$persistence), c("phi", "rmsd", "collision"))
  for (m in r$persistence) expect_identical(dim(m), c(W, W))
  expect_true(all(vapply(r$loadings, is.matrix, logical(1))))
  expect_true(all(vapply(r$pips, is.matrix, logical(1))))

  s <- summary(r)
  expect_s3_class(s, "summary.pefa")
  expect_identical(s$ssl, ssl(r))
  expect_identical(s$settings, r$settings)
  expect_output(print(r), "PEFA sweep \\(ordinary\\)")
  ## No display names a chosen count, because the package owns no count rule.
  expect_false(any(grepl("selected", capture.output(print(r)),
                         ignore.case = TRUE)))
  expect_output(print(s), "Candidates:")
})

test_that("fits carry the vbpm_fit class; print and coef dispatch", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, convChk = FALSE)
  expect_s3_class(f, "vbfa")
  expect_s3_class(f, "vbpm_fit")
  expect_identical(f$Q, d$Q)                    # design stored in the fit
  out <- capture.output(print(f))
  expect_true(any(grepl("vbfa", out)))
  expect_true(any(grepl("converged", out)))
  expect_identical(coef(f), f$Lam)
  ## $ access is untouched by the class
  expect_identical(f$Lam, unclass(f)$Lam)
})

test_that("renamed output fields are present and old names are gone", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, convChk = FALSE)
  expect_true(all(c("Lam_var", "eta_cov", "Phi_inv_mean") %in% names(f)))
  expect_false(any(c("sigsq.q.Lam", "PHI.q.eta", "M.q.PHIiver", "plotDat")
                   %in% names(f)))
  expect_equal(dim(f$Lam_var), dim(f$Lam))
  expect_equal(dim(f$eta_cov), c(ncol(d$Q), ncol(d$Q)))
})

test_that("vbfa is quiet by default and messages when convChk = TRUE", {
  d <- make_dat(N = 150)
  ## max_it = 5 will not converge, so silence the expected warning and test
  ## only the message channel
  expect_no_message(suppressWarnings(
    vbfa(d$Y, d$Q, v0 = 0.001, max_it = 5)))
  expect_message(suppressWarnings(
    vbfa(d$Y, d$Q, v0 = 0.001, max_it = 5, convChk = TRUE)),
    "relative error")
})

test_that("sim_fa produces data with planted local dependence", {
  s <- sim_fa(N = 300, K = 3, ipf = 6, lam = .7, lac = .3, ecr = .3, rseed = 7)
  expect_equal(dim(s$dat), c(300L, 18L))
  expect_false(anyNA(s$dat))
  expect_gt(nrow(s$ofd_ind), 0)
})

test_that("bifactor mode equals the hand-built general column bit-for-bit", {
  d  <- make_dat()
  fb <- vbfa(d$Y, d$Q, bifactor = TRUE)
  fh <- vbfa(d$Y, cbind(1L, d$Q), orthogonal = TRUE)
  expect_identical(fb$Lam, fh$Lam)
  expect_identical(fb$ELBO, fh$ELBO)
  expect_identical(fb$pi, fh$pi)
  ## bookkeeping fields
  expect_true(fb$bifactor); expect_identical(fb$n_general, 1L)
  expect_identical(fb$Q, cbind(rep(1L, nrow(d$Q)), d$Q, deparse.level = 0))
  expect_true(fb$orthogonal)
  ## print shows the general/group split
  expect_output(print(fb), "1 general [+] 3 group factors")
})

test_that("bifactor = TRUE overrides orthogonal = FALSE, with a message", {
  d <- make_dat()
  expect_message(
    f <- suppressWarnings(          # max_it deliberately tiny; iteration
      vbfa(d$Y, d$Q, bifactor = TRUE, orthogonal = FALSE, max_it = 50)),
    "overridden")
  expect_true(f$orthogonal)
  ## no message when orthogonal is left at its default
  expect_no_message(suppressWarnings(
    vbfa(d$Y, d$Q, bifactor = TRUE, max_it = 50)))
})

test_that("bifactor mode validates general and rejects a general column in Q", {
  d <- make_dat()
  expect_error(vbfa(d$Y, cbind(1L, d$Q), bifactor = TRUE), "all-specified")
  expect_error(vbfa(d$Y, d$Q, bifactor = TRUE, general = 2), "-1, 0, or 1")
  expect_error(vbfa(d$Y, d$Q, bifactor = TRUE, general = 0), "at least one")
  ## bifactor-(S-1)-style: zero general loadings where general = 0
  g <- rep(1L, ncol(d$Y)); g[1:4] <- 0L
  f <- suppressWarnings(            # truncated run; only the zeros matter
    vbfa(d$Y, d$Q, bifactor = TRUE, general = g, max_it = 200))
  expect_true(all(f$Lam[1:4, 1] == 0))
})

test_that("pefa bifactor sweep counts group factors without storing K_total", {
  d  <- make_dat()
  Q0 <- d$Q[, 1:2, drop = FALSE]
  r  <- pefa(Q0, d$Y, Kmin = 2, Kmax = 3, bifactor = TRUE, v0 = .001,
             max_it = 500, verbose = FALSE)
  expect_true(isTRUE(r$settings$bifactor))
  expect_identical(r$sweep$K, 2:3)
  expect_false("K_total" %in% names(r$sweep))
  expected_dims <- stats::setNames(3:4, c("2", "3"))
  expect_identical(vapply(r$loadings, ncol, integer(1)), expected_dims)
  expect_identical(vapply(r$pips, ncol, integer(1)), expected_dims)
  ## The adjacent row carries all seven structural facts; $persistence
  ## projects three of them over every direct pair.
  expect_true(all(c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                    "unmatched_ssl", "collision") %in% names(r$transitions)))
  expect_type(r$persistence$collision, "logical")
  expect_true(is.finite(r$persistence$phi["2", "3"]))
  expect_true(is.finite(r$persistence$rmsd["2", "3"]))
  expect_output(print(r), "group factors [(][+]1 general")
  ## K is the group-factor count everywhere; no total-factor column is stored
  ## or displayed.
  expect_false(any(grepl("K_total", capture.output(print(summary(r))))))
})

test_that("the fit object stores each quantity exactly once", {
  d <- make_dat()
  f <- vbfa(d$Y, d$Q, v0 = .001, max_it = 500, convChk = FALSE)
  ## the mirror lists removed in 0.7.0 stay gone
  for (nm in c("coefficients", "design", "posterior", "settings"))
    expect_null(f[[nm]])
  ## complete data: no mask stored, but the count is still reported
  expect_null(f$preprocess$missing_mask)
  expect_equal(f$preprocess$n_missing, 0)
  ## incomplete data: the mask is stored
  Y <- d$Y; Y[1:5, 1] <- NA
  fm <- vbfa(Y, d$Q, v0 = .001, max_it = 500, convChk = FALSE)
  expect_identical(fm$preprocess$missing_mask, is.na(Y))
  expect_equal(fm$preprocess$n_missing, 5)
  ## no two top-level components hold identical large content
  big <- Filter(function(z) is.matrix(z) && length(z) > 100, f)
  if (length(big) > 1) {
    pairs <- utils::combn(seq_along(big), 2)
    dupes <- apply(pairs, 2, function(ij) identical(big[[ij[1]]], big[[ij[2]]]))
    expect_false(any(dupes))
  }
  ## pefa keeps only compact loading/PIP matrices, never full candidate fits
  r <- pefa(d$Q[, 1:2, drop = FALSE], d$Y, Kmin = 2, Kmax = 3, v0 = .001,
            max_it = 500, verbose = FALSE)
  expect_null(r[["fit", exact = TRUE]])
  expect_null(r[["selected_fit", exact = TRUE]])
  expect_null(r[["fits", exact = TRUE]])
  expect_true(all(vapply(r$loadings, is.matrix, logical(1))))
  expect_true(all(vapply(r$pips, is.matrix, logical(1))))
})
