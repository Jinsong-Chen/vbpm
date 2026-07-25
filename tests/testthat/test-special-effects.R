## special_effects(): SL check, higher-order parameters, Eq.-16 effect sizes.

ho_fit <- function(N = 800, bifactor = TRUE) {
  sim <- sim_fa(N = N, K = 3, ipf = 6, lam = .8, gamma = c(.8, .7, .6),
                rseed = 11)
  Q <- matrix(-1L, 18, 3)
  for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
                   Q[a, ] <- 0L; Q[a, k] <- 1L }
  if (bifactor) vbfa(sim$dat, Q, bifactor = TRUE, v0 = .001)
  else          vbfa(sim$dat, cbind(1L, Q), orthogonal = TRUE, v0 = .001)
}

test_that("special_effects reproduces the vignette's plain-R computations", {
  fit <- ho_fit()
  se  <- special_effects(fit)
  cl  <- rep(1:3, each = 6)

  ## independent reimplementation (the vignette's code, known cluster ids)
  rt  <- sapply(1:18, function(j) abs(fit$Lam[j, cl[j] + 1]) / abs(fit$Lam[j, 1]))
  cvs <- sapply(1:3, function(k) { r <- rt[cl == k]; sd(r) / mean(r) })
  rbar <- sapply(1:3, function(k) mean(rt[cl == k]))
  Dk  <- sapply(1:3, function(k) {
    col <- k + 1
    inc <- (fit$Q[, col] == 1) | (fit$Q[, col] == -1 & fit$pi[, col] > .5)
    sum(fit$Lam[inc, col]^2) / sum(inc) })

  expect_equal(se$groups$cv, cvs)
  expect_equal(se$mean_cv, sum(rep(cvs, each = 6)) / 18)
  expect_equal(unname(se$gamma), 1 / sqrt(1 + rbar^2))
  expect_equal(se$groups$D, Dk)
  expect_equal(se$groups$n_items, rep(6L, 3))
  ## lambda = |b_gen| / gamma of the item's group
  expect_equal(se$items$lambda,
               abs(fit$Lam[, 1]) / (1 / sqrt(1 + rbar^2))[cl])
  expect_length(se$unassigned, 0)
  expect_identical(se$source, "bifactor")
  expect_output(print(se), "Item-weighted mean CV")
})

test_that("legacy hand-built fits give identical numbers", {
  se_b <- special_effects(ho_fit(bifactor = TRUE))
  se_l <- special_effects(ho_fit(bifactor = FALSE))
  expect_identical(se_l$source, "legacy")
  expect_equal(se_b$groups, se_l$groups)
  expect_equal(se_b$mean_cv, se_l$mean_cv)
})

test_that("special_effects rejects non-bifactor fits with guidance", {
  sim <- sim_fa(N = 300, K = 2, ipf = 5, cpf = 0, lam = .7, rseed = 3)
  Q <- matrix(-1L, 10, 2)
  Q[1:2, ] <- 0L; Q[1:2, 1] <- 1L; Q[6:7, ] <- 0L; Q[6:7, 2] <- 1L
  f_obl <- vbfa(sim$dat, Q, v0 = .001, max_it = 500, convChk = FALSE)
  expect_error(special_effects(f_obl), "oblique")
  f_ort <- vbfa(sim$dat, Q, orthogonal = TRUE, v0 = .001, max_it = 500,
                convChk = FALSE)
  expect_error(special_effects(f_ort), "all-specified")
  expect_error(special_effects(ho_fit(), tau = 2), "tau")
})

test_that("membership flags: unassigned, crossloaded, no_general", {
  fit <- ho_fit()
  ## doctor a copy: item 3 loses its group loading entirely
  f2 <- fit
  f2$pi[3, 2:4] <- 0.01
  se <- special_effects(f2)
  expect_true(3 %in% se$unassigned)
  expect_equal(se$groups$n_items[1], 5L)
  expect_true(is.na(se$items$ratio[3]))
  ## item 4 active on two group factors -> assigned to larger |loading|
  f3 <- fit
  f3$pi[4, 3] <- 0.9
  se3 <- special_effects(f3)
  expect_true(4 %in% se3$crossloaded)
  expect_identical(se3$items$group[4], "group1")   # own loading is larger
  ## bifactor-(S-1)-style: general fixed to 0 for cluster 1's items
  g <- rep(1L, 18); g[1:6] <- 0L
  fs1 <- suppressWarnings(
    vbfa(sim_fa(N = 400, K = 3, ipf = 6, lam = .8, gamma = c(.8, .7, .6),
                rseed = 11)$dat,
         { Q <- matrix(-1L, 18, 3)
           for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
                            Q[a, ] <- 0L; Q[a, k] <- 1L }; Q },
         bifactor = TRUE, general = g, v0 = .001, max_it = 500))
  ses1 <- special_effects(fs1)
  expect_identical(ses1$no_general, 1:6)
  expect_true(all(is.na(ses1$items$ratio[1:6])))
  expect_true(all(is.na(ses1$items$lambda[1:6])))
  expect_output(print(ses1), "no general loading")
})
