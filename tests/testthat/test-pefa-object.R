## The 0.9.0 pefa() object contract: the seven top-level components, the three
## derived tables, the failure and validation rules, and the public surface
## that survived the lean redesign (technical note section 6, plan 9.2/9.3).
##
## Fits are deliberately tiny -- J = 8, N = 150, a single v0 -- because nothing
## here depends on estimation quality, only on shape, type and bookkeeping.

.mk <- function(N = 150, K = 2, ipf = 4, seed = 11) {
  set.seed(seed)
  J <- K * ipf
  grp <- rep(seq_len(K), each = ipf)
  Lam <- matrix(0, J, K)
  Lam[cbind(seq_len(J), grp)] <- .75
  Phi <- matrix(.3, K, K)
  diag(Phi) <- 1
  psi <- pmax(1 - rowSums((Lam %*% chol(Phi))^2), .25)
  eta <- matrix(rnorm(N * K), N, K) %*% chol(Phi)
  Y <- eta %*% t(Lam) + matrix(rnorm(N * J), N, J) %*% diag(sqrt(psi))
  Q0 <- matrix(-1L, J, K)
  for (k in seq_len(K)) {
    a <- which(grp == k)[1:2]
    Q0[a, ] <- 0L
    Q0[a, k] <- 1L
  }
  list(Y = Y, Q0 = Q0, J = as.integer(J), K0 = as.integer(K))
}

.fit3 <- function(d = .mk()) {
  pefa(d$Q0, d$Y, 2, 4, v0 = .001, max_it = 400, verbose = FALSE)
}

## The regularization mask pefa() rebuilds for candidate K from the stored
## backbone.  Reassembly and the direct-comparison checks both need it.
.reg_of <- function(Q0, K) {
  J <- nrow(Q0)
  K0 <- ncol(Q0)
  Q <- if (K > K0) cbind(Q0, matrix(-1L, J, K - K0)) else Q0
  matrix(Q == -1L, J, K)
}

## Group-block view: bifactor candidates carry a leading general column that
## every between-candidate comparison removes first.
.grp <- function(m, bifactor) if (bifactor) m[, -1L, drop = FALSE] else m

.recompare <- function(x, K_from, K_to) {
  bif <- isTRUE(x$settings$bifactor)
  a <- as.character(K_from)
  b <- as.character(K_to)
  vbpm:::.pefa_compare_pair(
    lam_from = .grp(x$loadings[[a]], bif), lam_to = .grp(x$loadings[[b]], bif),
    pip_from = .grp(x$pips[[a]], bif), pip_to = .grp(x$pips[[b]], bif),
    reg_from = .reg_of(x$Q0, K_from), reg_to = .reg_of(x$Q0, K_to),
    K0 = ncol(x$Q0))
}

.fact_names <- c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                 "unmatched_ssl", "collision")

## One synthetic pair result, so the derived-view builders can be exercised on
## values chosen by hand rather than produced by a fit.
.hand_fact <- function(phi_min, rmsd, rmsd_max, ari, pip_rmsd, unmatched_ssl,
                       collision) {
  list(phi_min = phi_min, rmsd = rmsd, rmsd_max = rmsd_max, ari = ari,
       pip_rmsd = pip_rmsd, unmatched_ssl = unmatched_ssl,
       collision = collision)
}


# ---- top-level shape --------------------------------------------------

test_that("the object holds exactly the seven documented components", {
  r <- .fit3()
  expect_identical(names(r),
                   c("sweep", "transitions", "persistence", "loadings",
                     "pips", "Q0", "settings"))
  expect_identical(class(r), "pefa")
  expect_type(r, "list")
  expect_length(r, 7L)
})

test_that("every component removed in 0.9.0 is absent", {
  r <- .fit3()
  gone <- c("provenance", "selected_K", "boundary", "profile", "evidence",
            "checkpoint", "resume", "lineage", "hashes", "hash", "build_id",
            "fit", "fits", "selected_fit", "status", "reason", "n_pairs",
            "stability_eps", "pair_detail", "pairs", "decisions")
  for (nm in gone) expect_null(r[[nm, exact = TRUE]])

  ## Neither a soft criterion nor a status taxonomy survives in any table.
  expect_false(any(c("t_S", "AIC_S", "BIC_S", "status", "reason") %in%
                     names(r$sweep)))
  expect_false(any(c("status", "reason", "n_pairs", "assignment",
                     "collision_location", "unmatched_n") %in%
                     names(r$transitions)))
  ## No speculative sweep superclass.
  expect_false(inherits(r, "vbpm_sweep"))
})

test_that("$settings holds exactly the resolved fit controls, in order", {
  d <- .mk()
  r <- .fit3(d)
  s <- r$settings
  expect_identical(names(s),
                   c("bifactor", "general", "v0", "max_it", "convChk",
                     "tolVal", "tau", "rank_adjust", "rank_max_J"))
  expect_identical(s$bifactor, FALSE)
  expect_null(s$general)
  expect_identical(s$v0, 0.001)
  expect_identical(s$max_it, 400L)
  expect_identical(s$convChk, FALSE)
  expect_identical(s$tolVal, 1e-4)
  expect_identical(s$tau, 0.5)
  expect_identical(s$rank_adjust, FALSE)
  expect_identical(s$rank_max_J, 100L)
  ## Unnamed scalars, and no analysis policy or verbosity value.
  expect_null(names(s$bifactor))
  expect_null(names(s$v0))
  expect_false(any(c("verbose", "cuts", "sustain", "delta", "eps",
                     "stability_eps") %in% names(s)))

  ## Bifactor mode resolves `general` to a length-J integer design.
  rb <- pefa(d$Q0, d$Y, 2, 3, bifactor = TRUE, v0 = .001, max_it = 300,
             verbose = FALSE)
  expect_identical(rb$settings$bifactor, TRUE)
  expect_identical(rb$settings$general, rep(1L, d$J))
})

test_that("$Q0, $loadings and $pips are the documented candidate primitives", {
  d <- .mk()
  r <- .fit3(d)
  expect_true(is.matrix(r$Q0))
  expect_identical(dim(r$Q0), c(d$J, 2L))
  expect_type(r$Q0, "integer")
  expect_true(all(r$Q0 %in% c(-1L, 0L, 1L)))
  expect_identical(colnames(r$Q0), c("F1", "F2"))

  keys <- as.character(2:4)
  expect_identical(names(r$loadings), keys)
  expect_identical(names(r$pips), keys)
  for (k in keys) {
    L <- r$loadings[[k]]
    P <- r$pips[[k]]
    expect_true(is.matrix(L) && is.double(L))
    expect_true(is.matrix(P) && is.double(P))
    expect_identical(dim(L), c(d$J, as.integer(k)))
    expect_identical(dim(P), dim(L))
    expect_identical(dimnames(P), dimnames(L))
    expect_identical(rownames(L), rownames(r$Q0))
    expect_true(all(is.finite(L)))
    expect_true(all(is.finite(P)) && all(P >= 0) && all(P <= 1))
  }
  ## Backbone names lead the factor dimnames; exploratory columns continue it.
  expect_identical(colnames(r$loadings[["4"]]), c("F1", "F2", "F3", "F4"))
})


# ---- $sweep -----------------------------------------------------------

test_that("$sweep has exactly the documented columns and types", {
  r <- .fit3()
  expect_s3_class(r$sweep, "data.frame")
  expect_identical(names(r$sweep),
                   c("K", "ELBO", "AIC", "BIC", "RMSEA", "SRMR", "CFI", "TLI",
                     "t", "iter", "converged"))
  expect_identical(nrow(r$sweep), 3L)
  expect_identical(r$sweep$K, 2:4)
  expect_type(r$sweep$K, "integer")
  expect_type(r$sweep$iter, "integer")
  expect_type(r$sweep$converged, "logical")
  for (nm in c("ELBO", "AIC", "BIC", "RMSEA", "SRMR", "CFI", "TLI", "t")) {
    expect_type(r$sweep[[nm]], "double")
  }
  ## Required candidate facts are finite; t is a whole nonnegative count that
  ## is nonetheless stored as a double.
  expect_true(all(is.finite(r$sweep$ELBO)))
  expect_true(all(is.finite(r$sweep$AIC)))
  expect_true(all(is.finite(r$sweep$BIC)))
  expect_true(all(is.finite(r$sweep$t)))
  expect_true(all(r$sweep$t >= 0))
  expect_identical(r$sweep$t, round(r$sweep$t))
  expect_false(anyNA(r$sweep$converged))
  expect_identical(rownames(r$sweep), as.character(1:3))
})


# ---- $transitions -----------------------------------------------------

test_that("$transitions has exactly the documented columns and types", {
  r <- .fit3()
  tr <- r$transitions
  expect_s3_class(tr, "data.frame")
  expect_identical(names(tr),
                   c("K_from", "K_to", "ELBO_gain_pct", "BIC_gain_pct",
                     .fact_names))
  expect_identical(ncol(tr), 11L)
  expect_identical(nrow(tr), nrow(r$sweep) - 1L)
  expect_identical(tr$K_from, 2:3)
  expect_identical(tr$K_to, 3:4)
  expect_type(tr$K_from, "integer")
  expect_type(tr$K_to, "integer")
  expect_type(tr$collision, "logical")
  for (nm in c("ELBO_gain_pct", "BIC_gain_pct",
               setdiff(.fact_names, "collision"))) {
    expect_type(tr[[nm]], "double")
  }
})

test_that("a one-candidate window keeps the typed zero-row schema", {
  d <- .mk()
  r <- pefa(d$Q0, d$Y, 2, 2, v0 = .001, max_it = 300, verbose = FALSE)
  expect_identical(nrow(r$sweep), 1L)
  tr <- r$transitions
  expect_identical(nrow(tr), 0L)
  expect_identical(names(tr),
                   c("K_from", "K_to", "ELBO_gain_pct", "BIC_gain_pct",
                     .fact_names))
  expect_type(tr$K_from, "integer")
  expect_type(tr$K_to, "integer")
  expect_type(tr$collision, "logical")
  for (nm in c("ELBO_gain_pct", "BIC_gain_pct",
               setdiff(.fact_names, "collision"))) {
    expect_type(tr[[nm]], "double")
  }
  ## The persistence triangles degenerate to a single typed NA cell.
  expect_identical(dim(r$persistence$phi), c(1L, 1L))
  expect_identical(r$persistence$phi[1L, 1L], NA_real_)
  expect_identical(r$persistence$collision[1L, 1L], NA)
})


# ---- $persistence -----------------------------------------------------

test_that("$persistence is three plain upper-triangular K-keyed matrices", {
  r <- .fit3()
  p <- r$persistence
  expect_type(p, "list")
  expect_identical(names(p), c("phi", "rmsd", "collision"))
  keys <- as.character(2:4)
  for (nm in names(p)) {
    m <- p[[nm]]
    expect_true(is.matrix(m))
    expect_identical(dim(m), c(3L, 3L))
    expect_identical(dimnames(m), list(keys, keys))
    ## Ordinary matrices only: no custom class, no smuggled attribute.
    expect_identical(names(attributes(m)), c("dim", "dimnames"))
    expect_identical(class(m), c("matrix", "array"))
    ## Diagonal and lower triangle are typed NA.
    expect_true(all(is.na(m[!upper.tri(m)])))
    expect_true(all(!is.na(m[upper.tri(m)])))
  }
  expect_type(p$phi, "double")
  expect_type(p$rmsd, "double")
  expect_type(p$collision, "logical")
  expect_identical(p$phi[2L, 1L], NA_real_)
  expect_identical(p$collision[3L, 3L], NA)
})

test_that("every adjacent persistence cell is its transition field", {
  r <- .fit3()
  tr <- r$transitions
  p <- r$persistence
  for (i in seq_len(nrow(tr))) {
    a <- as.character(tr$K_from[i])
    b <- as.character(tr$K_to[i])
    expect_identical(p$phi[a, b], tr$phi_min[i])
    ## The mapping is rmsd_max, NOT the pooled $transitions$rmsd.
    expect_identical(p$rmsd[a, b], tr$rmsd_max[i])
    expect_identical(p$collision[a, b], tr$collision[i])
  }
  ## The two RMSD definitions really do differ on this sweep, so the check
  ## above could not have passed by accident.
  expect_true(any(tr$rmsd != tr$rmsd_max))
  expect_false(isTRUE(all.equal(p$rmsd["2", "3"], tr$rmsd[1L])))
})

test_that("the derived views read one shared pair result, keyed correctly", {
  ## Hand-chosen facts with rmsd and rmsd_max deliberately unequal, so a
  ## crossed mapping cannot survive.
  Ks <- 2:4
  facts <- list(
    "2-3" = .hand_fact(0.91, 0.11, 0.21, 0.31, 0.41, 0.51, FALSE),
    "3-4" = .hand_fact(0.92, 0.12, 0.22, 0.32, 0.42, 0.52, TRUE),
    "2-4" = .hand_fact(0.93, 0.13, 0.23, 0.33, 0.43, 0.53, TRUE)
  )
  sweep <- data.frame(K = Ks, ELBO = c(-100, -50, -60), BIC = c(300, 280, 290))

  tr <- vbpm:::.pefa_transitions(sweep, facts)
  expect_identical(tr$K_from, 2:3)
  expect_identical(tr$phi_min, c(0.91, 0.92))
  expect_identical(tr$rmsd, c(0.11, 0.12))
  expect_identical(tr$rmsd_max, c(0.21, 0.22))
  expect_identical(tr$ari, c(0.31, 0.32))
  expect_identical(tr$pip_rmsd, c(0.41, 0.42))
  expect_identical(tr$unmatched_ssl, c(0.51, 0.52))
  expect_identical(tr$collision, c(FALSE, TRUE))

  p <- vbpm:::.pefa_persistence(Ks, facts)
  expect_identical(p$phi["2", "3"], 0.91)
  expect_identical(p$phi["2", "4"], 0.93)
  expect_identical(p$rmsd["2", "3"], 0.21)      # rmsd_max, not 0.11
  expect_identical(p$rmsd["3", "4"], 0.22)
  expect_identical(p$rmsd["2", "4"], 0.23)
  ## A collision on a non-adjacent pair reaches the mask on its own cell.
  expect_identical(p$collision["2", "4"], TRUE)
  expect_identical(p$collision["2", "3"], FALSE)
})

test_that("all-pairs cells are direct endpoint comparisons", {
  r <- .fit3()
  Ks <- r$sweep$K
  for (i in seq_along(Ks)) {
    for (j in seq_along(Ks)) {
      if (i >= j) next
      f <- .recompare(r, Ks[i], Ks[j])
      a <- as.character(Ks[i])
      b <- as.character(Ks[j])
      expect_identical(r$persistence$phi[a, b], f$phi_min)
      expect_identical(r$persistence$rmsd[a, b], f$rmsd_max)
      expect_identical(r$persistence$collision[a, b], f$collision)
    }
  }
})

test_that("a chained composition would give a different all-pairs answer", {
  ## Three synthetic candidates in which composing the adjacent assignments
  ## disagrees with the direct one: the single K = 1 source matches column 1
  ## of the K = 2 candidate, that column matches column 2 of the K = 3
  ## candidate, but the source itself matches column 1 of the K = 3 candidate.
  L <- list("1" = cbind(c(1, 0, 0)),
            "2" = cbind(c(1, 1, 0), c(0, 0, 1)),
            "3" = cbind(c(1, 0, 0), c(1, 1, 0), c(0, 0, 1)))
  P <- lapply(L, function(m) matrix(0.5, nrow(m), ncol(m)))
  Q0 <- matrix(0L, 3L, 0L)
  cmp <- function(a, b) {
    vbpm:::.pefa_compare_pair(L[[a]], L[[b]], P[[a]], P[[b]],
                              .reg_of(Q0, as.integer(a)),
                              .reg_of(Q0, as.integer(b)), 0L)
  }
  facts <- list("1-2" = cmp("1", "2"), "2-3" = cmp("2", "3"),
                "1-3" = cmp("1", "3"))
  p <- vbpm:::.pefa_persistence(1:3, facts)

  ## Direct: an exact copy of the source, so the weakest-column RMSD is 0.
  expect_identical(p$rmsd["1", "3"], 0)
  ## Chained would have gone through the K = 2 column and landed on the
  ## K = 3 column whose distance is sqrt(1/3).
  expect_equal(facts[["1-2"]]$rmsd_max, sqrt(1 / 3))
  expect_gt(facts[["1-2"]]$rmsd_max, p$rmsd["1", "3"])
})


# ---- gains ------------------------------------------------------------

test_that("gains are oriented toward K_to and scaled by their own maximum", {
  Ks <- 2:5
  facts <- stats::setNames(
    rep(list(.hand_fact(1, 0, 0, 1, 0, 0, FALSE)), 3L),
    c("2-3", "3-4", "4-5"))
  sweep <- data.frame(K = Ks,
                      ELBO = c(-100, -50, -60, -40),
                      BIC = c(300, 280, 290, 285))
  tr <- vbpm:::.pefa_transitions(sweep, facts)

  ## ELBO gain is ELBO_to - ELBO_from: 50, -10, 20, largest positive 50.
  expect_identical(tr$ELBO_gain_pct, c(100, -20, 40))
  ## BIC gain is BIC_from - BIC_to: 20, -10, 5, largest positive 20.
  expect_identical(tr$BIC_gain_pct, c(100, -50, 25))
})

test_that("a path with no positive gain is entirely NA", {
  Ks <- 2:4
  facts <- stats::setNames(
    rep(list(.hand_fact(1, 0, 0, 1, 0, 0, FALSE)), 2L), c("2-3", "3-4"))
  sweep <- data.frame(K = Ks, ELBO = c(-100, -150, -200),
                      BIC = c(300, 280, 260))
  tr <- vbpm:::.pefa_transitions(sweep, facts)
  expect_identical(tr$ELBO_gain_pct, c(NA_real_, NA_real_))
  ## The BIC path is unaffected: normalization is per path.
  expect_identical(tr$BIC_gain_pct, c(100, 100))

  ## A real all-nonpositive sweep reaches the same place.
  r <- .fit3()
  expect_true(all(diff(r$sweep$ELBO) < 0))
  expect_true(all(is.na(r$transitions$ELBO_gain_pct)))
  expect_true(all(is.na(r$transitions$BIC_gain_pct)))
})

test_that("a nonrepresentable adjacent difference is a controlled error", {
  facts <- list("2-3" = .hand_fact(1, 0, 0, 1, 0, 0, FALSE))
  sweep <- data.frame(K = 2:3, ELBO = c(-1e308, 1e308), BIC = c(1, 2))
  expect_error(vbpm:::.pefa_transitions(sweep, facts), "must be finite")
  sweep2 <- data.frame(K = 2:3, ELBO = c(1, 2), BIC = c(-1e308, 1e308))
  expect_error(vbpm:::.pefa_transitions(sweep2, facts), "must be finite")
})


# ---- reassembly -------------------------------------------------------

test_that("both derived views rebuild from the five source components", {
  r <- .fit3()
  local_mocked_bindings(
    vbfa = function(...) {
      stop("vbfa() must not be called when rebuilding a derived view.")
    },
    .package = "vbpm")

  Ks <- r$sweep$K
  K0 <- ncol(r$Q0)
  bif <- isTRUE(r$settings$bifactor)
  facts <- list()
  for (i in seq_along(Ks)) {
    for (j in seq_along(Ks)) {
      if (i >= j) next
      a <- as.character(Ks[i])
      b <- as.character(Ks[j])
      facts[[paste(Ks[i], Ks[j], sep = "-")]] <- vbpm:::.pefa_compare_pair(
        lam_from = .grp(r$loadings[[a]], bif),
        lam_to = .grp(r$loadings[[b]], bif),
        pip_from = .grp(r$pips[[a]], bif), pip_to = .grp(r$pips[[b]], bif),
        reg_from = .reg_of(r$Q0, Ks[i]), reg_to = .reg_of(r$Q0, Ks[j]),
        K0 = K0)
    }
  }
  expect_identical(vbpm:::.pefa_transitions(r$sweep, facts), r$transitions)
  expect_identical(vbpm:::.pefa_persistence(Ks, facts), r$persistence)
})


# ---- failure contract -------------------------------------------------

test_that("a negative or fractional t aborts the call naming the K", {
  ## t counts parameters, so these are malformed candidates rather than
  ## unusual ones (technical note 6.1, plan 9.2).
  d <- .mk()
  real <- vbpm::fit_stats
  for (bad in c(-5, 12.5)) {
    local({
      value <- bad
      local_mocked_bindings(
        fit_stats = function(object, ...) {
          out <- real(object, ...)
          if (ncol(object$Lam) == 3L) out[["t"]] <- value
          out
        },
        .package = "vbpm")
      expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200,
                        verbose = FALSE),
                   "PEFA candidate K = 3: t is not a nonnegative whole number")
    })
  }
})

test_that("a nonfinite required criterion aborts the call naming the K", {
  d <- .mk()
  real <- vbpm::fit_stats
  for (broken in c("AIC", "BIC", "t")) {
    local({
      field <- broken
      local_mocked_bindings(
        fit_stats = function(object, ...) {
          out <- real(object, ...)
          if (ncol(object$Lam) == 3L) out[[field]] <- Inf
          out
        },
        .package = "vbpm")
      expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200,
                        verbose = FALSE),
                   "PEFA candidate K = 3: AIC, BIC, or t is nonfinite")
    })
  }
})

test_that("a nonfinite ELBO aborts the call naming the K", {
  d <- .mk()
  real <- vbpm::vbfa
  local_mocked_bindings(
    vbfa = function(...) {
      f <- real(...)
      if (ncol(f$Lam) == 3L) f$ELBO <- NaN
      f
    },
    .package = "vbpm")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200,
                    verbose = FALSE),
               "PEFA candidate K = 3: ELBO is not one finite real number")
})

test_that("a malformed candidate matrix aborts the call naming the K", {
  ## fit_stats() is stubbed out here so the loading matrix itself is the only
  ## thing that can fail: the point is the candidate-normalization check, not
  ## how a fit statistic reacts to a broken matrix.
  d <- .mk()
  real <- vbpm::vbfa
  local_mocked_bindings(
    vbfa = function(...) {
      f <- real(...)
      if (ncol(f$Lam) == 3L) f$Lam[1L, 1L] <- Inf
      f
    },
    fit_stats = function(...) {
      c(AIC = 1, BIC = 2, t = 3, RMSEA = 0, SRMR = 0, CFI = 1, TLI = 1)
    },
    .package = "vbpm")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200,
                    verbose = FALSE),
               "PEFA candidate K = 3: the loading matrix has nonfinite")
})

test_that("a dimensionally inconsistent candidate matrix names its K", {
  d <- .mk()
  real <- vbpm::vbfa
  local_mocked_bindings(
    vbfa = function(...) {
      f <- real(...)
      if (ncol(f$Lam) == 2L) f$Lam <- f$Lam[, 1L, drop = FALSE]
      f
    },
    fit_stats = function(...) {
      c(AIC = 1, BIC = 2, t = 3, RMSEA = 0, SRMR = 0, CFI = 1, TLI = 1)
    },
    .package = "vbpm")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200,
                    verbose = FALSE),
               "PEFA candidate K = 2: the loading matrix is 8 by 1")
})

test_that("a PIP outside [0, 1] aborts the call naming the K", {
  d <- .mk()
  real <- vbpm::vbfa
  local_mocked_bindings(
    vbfa = function(...) {
      f <- real(...)
      if (ncol(f$Lam) == 2L) f$pi[1L, 1L] <- 1.5
      f
    },
    .package = "vbpm")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200,
                    verbose = FALSE),
               "PEFA candidate K = 2: the PIP matrix has values outside")
})

test_that("a vbfa() failure aborts the call naming the K", {
  d <- .mk()
  local_mocked_bindings(
    vbfa = function(Y, Q, ...) {
      if (ncol(Q) == 3L) stop("synthetic engine failure")
      stop("unexpected K")
    },
    .package = "vbpm")
  expect_error(pefa(d$Q0, d$Y, 3, 3, v0 = .001, max_it = 200,
                    verbose = FALSE),
               "PEFA candidate K = 3: vbfa\\(\\) failed: synthetic engine")
})

test_that("descriptive indices may be NA without erroring", {
  d <- .mk()
  real <- vbpm::fit_stats
  local_mocked_bindings(
    fit_stats = function(object, ...) {
      out <- real(object, ...)
      out[["RMSEA"]] <- NaN
      out[["TLI"]] <- NA_real_
      out
    },
    .package = "vbpm")
  r <- pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 200, verbose = FALSE)
  ## NaN and NA alike are normalized to the documented NA_real_ flavor.
  expect_identical(r$sweep$RMSEA, rep(NA_real_, 2L))
  expect_identical(r$sweep$TLI, rep(NA_real_, 2L))
  expect_false(any(is.nan(r$sweep$RMSEA)))
  expect_true(all(is.finite(r$sweep$SRMR)))
  expect_true(all(is.finite(r$sweep$AIC)))
})

test_that("nonconverged candidates are retained under one aggregate warning", {
  d <- .mk()
  seen <- character(0)
  r <- withCallingHandlers(
    pefa(d$Q0, d$Y, 2, 3, v0 = .001, max_it = 3, verbose = FALSE),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  ## Exactly one warning for the whole call, naming every affected K, and the
  ## per-candidate iteration-limit warnings are muffled.
  expect_length(seen, 1L)
  expect_match(seen, "did not converge at K = 2, 3")
  ## The candidates themselves survive with their numbers intact.
  expect_s3_class(r, "pefa")
  expect_identical(r$sweep$converged, c(FALSE, FALSE))
  expect_true(all(is.finite(r$sweep$ELBO)))
  expect_true(all(is.finite(r$sweep$BIC)))
  expect_identical(nrow(r$transitions), 1L)
  expect_identical(summary(r)$nonconverged_K, 2:3)
  expect_output(print(r), "nonconverged K: 2, 3")
})


# ---- argument validation ----------------------------------------------

test_that("the K window is validated before anything is fitted", {
  d <- .mk()
  expect_error(pefa(d$Q0, d$Y, 4, 2, verbose = FALSE), "non-decreasing")
  whole <- "must be a single whole number"
  expect_error(pefa(d$Q0, d$Y, 2.5, 4, verbose = FALSE),
               paste("Kmin", whole))
  expect_error(pefa(d$Q0, d$Y, 2, 3.5, verbose = FALSE),
               paste("Kmax", whole))
  expect_error(pefa(d$Q0, d$Y, Inf, 4, verbose = FALSE),
               paste("Kmin", whole))
  expect_error(pefa(d$Q0, d$Y, NA_integer_, 4, verbose = FALSE), "Kmin")
  expect_error(pefa(d$Q0, d$Y, NaN, 4, verbose = FALSE), "Kmin")
  expect_error(pefa(d$Q0, d$Y, 2 + 0i, 4, verbose = FALSE), "Kmin")
  expect_error(pefa(d$Q0, d$Y, "2", 4, verbose = FALSE), "Kmin")
  expect_error(pefa(d$Q0, d$Y, c(2, 3), 4, verbose = FALSE), "Kmin")
  expect_error(pefa(d$Q0, d$Y, 0, 4, verbose = FALSE), "Kmin")
  expect_error(pefa(d$Q0, d$Y, 2, 2^31, verbose = FALSE), "Kmax")
  ## The window must reach the backbone.
  expect_error(pefa(d$Q0, d$Y, 1, 3, verbose = FALSE),
               "must reach the backbone")
})

test_that("Q0 is validated", {
  d <- .mk()
  bad <- d$Q0
  bad[1L, 1L] <- 2L
  expect_error(pefa(bad, d$Y, 2, 3, verbose = FALSE), "-1, 0, or 1")
  bad2 <- d$Q0
  bad2[1L, 1L] <- NA_integer_
  expect_error(pefa(bad2, d$Y, 2, 3, verbose = FALSE), "missing values")
  expect_error(pefa(d$Q0[-1L, , drop = FALSE], d$Y, 2, 3, verbose = FALSE),
               "nrow\\(Q0\\) must equal ncol\\(Y\\)")
  ## Supplied item names must agree with the data rather than silently win.
  named <- d$Q0
  rownames(named) <- paste0("z", seq_len(nrow(named)))
  expect_error(pefa(named, d$Y, 2, 3, verbose = FALSE), "rownames\\(Q0\\)")
})

test_that("Y is validated", {
  d <- .mk()
  bad <- d$Y
  bad[1L, 1L] <- Inf
  expect_error(pefa(d$Q0, bad, 2, 3, verbose = FALSE), "not NaN or infinite")
  bad[1L, 1L] <- NaN
  expect_error(pefa(d$Q0, bad, 2, 3, verbose = FALSE), "not NaN or infinite")
  chr <- matrix(as.character(d$Y), nrow(d$Y), ncol(d$Y))
  expect_error(pefa(d$Q0, chr, 2, 3, verbose = FALSE),
               "nonempty numeric N by J matrix")
  ## An all-missing item is rejected downstream, and the K is still named.
  gap <- d$Y
  gap[, 1L] <- NA
  expect_error(pefa(d$Q0, gap, 2, 2, v0 = .001, verbose = FALSE),
               "PEFA candidate K = 2")
})

test_that("the remaining fit controls are validated", {
  d <- .mk()
  expect_error(pefa(d$Q0, d$Y, 2, 3, tau = 1.5, verbose = FALSE), "tau")
  expect_error(pefa(d$Q0, d$Y, 2, 3, tau = NA_real_, verbose = FALSE), "tau")
  expect_error(pefa(d$Q0, d$Y, 2, 3, tolVal = -1, verbose = FALSE), "tolVal")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = 0, verbose = FALSE), "v0")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = numeric(0), verbose = FALSE), "v0")
  expect_error(pefa(d$Q0, d$Y, 2, 3, v0 = c(.01, NA), verbose = FALSE), "v0")
  expect_error(pefa(d$Q0, d$Y, 2, 3, max_it = 0, verbose = FALSE), "max_it")
  expect_error(pefa(d$Q0, d$Y, 2, 3, rank_max_J = -1, verbose = FALSE),
               "rank_max_J")
  expect_error(pefa(d$Q0, d$Y, 2, 3, convChk = NA, verbose = FALSE),
               "convChk must be TRUE or FALSE")
  expect_error(pefa(d$Q0, d$Y, 2, 3, verbose = "yes"),
               "verbose must be TRUE or FALSE")
  expect_error(pefa(d$Q0, d$Y, 2, 3, bifactor = TRUE, general = 2,
                    verbose = FALSE), "-1, 0, or 1")
  expect_error(pefa(d$Q0, d$Y, 2, 3, bifactor = TRUE, general = 0,
                    verbose = FALSE), "at least one general loading")
  expect_error(pefa(d$Q0, d$Y, 2, 3, bifactor = TRUE,
                    general = rep(1L, 3L), verbose = FALSE), "length-J")
})


# ---- public surface ---------------------------------------------------

test_that("the formal surface matches the technical note exactly", {
  f <- formals(pefa)
  expect_identical(names(f)[1:4], c("Q0", "Y", "Kmin", "Kmax"))
  expect_identical(names(f),
                   c("Q0", "Y", "Kmin", "Kmax", "bifactor", "general", "v0",
                     "max_it", "convChk", "tolVal", "tau", "rank_adjust",
                     "rank_max_J", "verbose"))
  expect_identical(f$bifactor, FALSE)
  expect_identical(eval(f$general), 1)
  expect_identical(eval(f$v0), c(.01, .005, .002, .001))
  expect_identical(eval(f$max_it), 10000)
  expect_identical(f$convChk, FALSE)
  expect_identical(eval(f$tolVal), 1e-4)
  expect_identical(eval(f$tau), .50)
  expect_identical(f$rank_adjust, FALSE)
  expect_identical(eval(f$rank_max_J), 100)
  expect_identical(f$verbose, TRUE)

  ## No selector, cut, stability, checkpoint, save or dots argument survives.
  expect_false(any(c("...", "cuts", "sustain", "delta", "stability_eps",
                     "checkpoint", "resume", "save", "file", "select",
                     "ld", "Qe", "orthogonal") %in% names(f)))
})

test_that("every removed API is gone from the namespace and the exports", {
  ns <- asNamespace("vbpm")
  removed <- c("select_K_elbow", "verify_pefa", "pair_detail", "extend_pefa",
               "persistence", "persistence.pefa", "pefa_provenance",
               "pefa_evidence")
  for (nm in removed) {
    expect_false(exists(nm, envir = ns, inherits = FALSE),
                 label = paste0("namespace object ", nm))
    expect_false(nm %in% getNamespaceExports("vbpm"),
                 label = paste0("export ", nm))
  }
  ## The exported surface is exactly the delivered one.
  expect_identical(sort(getNamespaceExports("vbpm")),
                   sort(c("fit_stats", "pefa", "sim_fa", "sim_lvm",
                          "special_effects", "ssl", "vbfa", "vbmimic")))
  ## The retired matcher/provenance internals are gone too.
  for (nm in c(".pefa_pair_facts", ".pefa_aligned_rmsd", ".pefa_provenance",
               ".pefa_verify", ".pefa_hash", ".pefa_status")) {
    expect_false(exists(nm, envir = ns, inherits = FALSE), label = nm)
  }
})

test_that("ssl() is exactly colSums(L^2) and never screens", {
  d <- .mk()
  r <- .fit3(d)
  s <- ssl(r)
  expect_type(s, "list")
  expect_identical(names(s), names(r$loadings))
  for (k in names(s)) {
    expect_type(s[[k]], "double")
    expect_identical(s[[k]], colSums(r$loadings[[k]]^2))
    expect_identical(names(s[[k]]), colnames(r$loadings[[k]]))
  }
  ## Bifactor keeps the general column first and reports K + 1 entries.
  rb <- pefa(d$Q0, d$Y, 2, 3, bifactor = TRUE, v0 = .001, max_it = 300,
             verbose = FALSE)
  sb <- ssl(rb)
  expect_identical(lengths(sb), c("2" = 3L, "3" = 4L))
  expect_identical(names(sb[["3"]]), c("G", "F1", "F2", "F3"))
  expect_identical(sb[["3"]], colSums(rb$loadings[["3"]]^2))

  expect_error(ssl(unclass(r)), "must be a \"pefa\" object")
})


# ---- covered fitting cases --------------------------------------------

test_that("a fully exploratory window (K0 = 0) sweeps and compares", {
  d <- .mk()
  r <- pefa(matrix(0L, d$J, 0L), d$Y, 1, 3, v0 = .001, max_it = 300,
            verbose = FALSE)
  expect_identical(dim(r$Q0), c(d$J, 0L))
  expect_identical(r$sweep$K, 1:3)
  expect_identical(nrow(r$transitions), 2L)
  expect_true(all(is.finite(r$transitions$rmsd)))
  ## With no backbone every source column is exploratory, so a positive
  ## ELBO gain exists and the percentage path is populated.
  expect_true(any(is.finite(r$transitions$ELBO_gain_pct)))
  expect_identical(max(r$transitions$ELBO_gain_pct), 100)
  ## Every all-pairs cell is still a direct comparison.
  expect_identical(r$persistence$phi["1", "3"], .recompare(r, 1, 3)$phi_min)
})

test_that("a backbone-only source window is a legal edge", {
  ## Kmin == K0, so the K = 2 candidate has no exploratory source column and
  ## every pair leaving it is pure position matching.
  d <- .mk()
  r <- .fit3(d)
  expect_identical(ncol(r$Q0), 2L)
  expect_identical(min(r$sweep$K), 2L)
  expect_false(r$transitions$collision[1L])
  expect_true(is.finite(r$transitions$unmatched_ssl[1L]))
})

test_that("a bifactor sweep strips the general column without storing less", {
  d <- .mk()
  r <- pefa(d$Q0, d$Y, 2, 4, bifactor = TRUE, v0 = .001, max_it = 300,
            verbose = FALSE)
  expect_identical(r$sweep$K, 2:4)
  ## K counts group factors; the stored matrices keep the general column.
  expect_identical(vapply(r$loadings, ncol, integer(1)),
                   c("2" = 3L, "3" = 4L, "4" = 5L))
  expect_identical(vapply(r$pips, ncol, integer(1)),
                   c("2" = 3L, "3" = 4L, "4" = 5L))
  expect_identical(colnames(r$loadings[["3"]]), c("G", "F1", "F2", "F3"))
  ## Comparisons operate on the group block only, which is what makes them
  ## the same measurement as an ordinary sweep's.
  expect_identical(r$transitions$phi_min[1L], .recompare(r, 2, 3)$phi_min)
  expect_identical(r$persistence$rmsd["2", "4"], .recompare(r, 2, 4)$rmsd_max)
  ## Stripping is done on a copy: the stored matrices still have K + 1 columns.
  expect_identical(ncol(r$loadings[["4"]]), 5L)
})


# ---- display methods --------------------------------------------------

test_that("print, summary and print.summary read only stored fields", {
  d <- .mk()
  r <- .fit3(d)
  expect_output(print(r), "PEFA sweep \\(ordinary\\): K = 2:4, K0 = 2")
  expect_output(print(r), "sweep 3 x 11 \\| transitions 2 x 11")
  expect_output(print(r), "all candidates converged")
  invisible(capture.output(shown <- withVisible(print(r))))
  expect_false(shown$visible)

  s <- summary(r)
  expect_s3_class(s, "summary.pefa")
  expect_identical(names(s),
                   c("window", "sweep", "transitions", "persistence", "ssl",
                     "settings", "nonconverged_K"))
  expect_identical(s$window, c(Kmin = 2L, Kmax = 4L, K0 = 2L))
  expect_identical(s$sweep, r$sweep)
  expect_identical(s$transitions, r$transitions)
  expect_identical(s$persistence, r$persistence)
  expect_identical(s$ssl, ssl(r))
  expect_identical(s$settings, r$settings)
  expect_identical(s$nonconverged_K, integer(0))

  out <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(out, "Candidates:")
  expect_match(out, "SS loadings")
  expect_match(out, "Transitions:")
  expect_match(out, "Persistence")
  ## phi is displayed before rmsd.
  expect_lt(regexpr("phi \\(phi_min\\)", out),
            regexpr("rmsd \\(rmsd_max\\)", out))
  ## No method claims a count or a threshold.
  expect_false(grepl("selected|recommend|threshold|cut", out,
                     ignore.case = TRUE))
})

test_that("a collided persistence cell is marked for display only", {
  r <- .fit3()
  marked <- r
  marked$persistence$collision["2", "3"] <- TRUE
  out <- paste(capture.output(print(summary(marked))), collapse = "\n")
  expect_match(out, "\\*")
  expect_match(out, "target reuse")
  ## The stored triangles keep their storage modes: the marker is text only.
  expect_type(marked$persistence$phi, "double")
  expect_type(marked$persistence$collision, "logical")
  expect_identical(marked$persistence$phi, r$persistence$phi)

  ## With no collision anywhere the legend line is not printed.
  clean <- paste(capture.output(print(summary(r))), collapse = "\n")
  expect_false(grepl("target reuse", clean))
})

test_that("the compact displays run on bifactor and one-candidate objects", {
  d <- .mk()
  rb <- pefa(d$Q0, d$Y, 2, 3, bifactor = TRUE, v0 = .001, max_it = 300,
             verbose = FALSE)
  expect_output(print(rb), "bifactor")
  expect_output(print(rb), "group factors \\(\\+1 general\\)")
  expect_output(print(summary(rb)), "SS loadings .*general column first")

  r1 <- pefa(d$Q0, d$Y, 2, 2, v0 = .001, max_it = 300, verbose = FALSE)
  expect_output(print(r1), "sweep 1 x 11 \\| transitions 0 x 11")
  expect_output(print(summary(r1)), "single candidate")

  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 500, height = 400)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)
  ## A one-candidate gain view is an empty panel rather than an error.
  expect_silent(plot(r1, type = "gain", criterion = "ELBO"))
  expect_silent(plot(rb, type = "objective", criterion = "BIC"))
  expect_silent(plot(rb, type = "gain", criterion = "BIC", pct = TRUE))
  expect_silent(plot(rb, type = "fit"))
  ## The gain view has no AIC path, and no removed criterion is accepted.
  expect_error(plot(rb, type = "gain", criterion = "AIC"), "unavailable")
  expect_error(plot(rb, type = "gain", criterion = "BIC_S"), "unavailable")
  expect_error(plot(rb, type = "stability"), "arg")
})
