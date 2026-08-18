## The 0.9.0 evidence object: its exact shape, the three tables, the component
## status codebook, the one-warning rule, and the bifactor group-block offset.
## Fits are deliberately tiny and use a single spike (v0 = .001), and every
## table-level fact that can be checked without fitting is checked on a
## hand-built sweep instead.

.pefa_toy <- function(N = 60, J = 9, seed = 3) {
  set.seed(seed)
  K <- 3L
  group <- rep(seq_len(K), each = J / K)
  loading <- matrix(0, J, K)
  loading[cbind(seq_len(J), group)] <- .8
  scores <- matrix(rnorm(N * K), N, K)
  Y <- scores %*% t(loading) +
    matrix(rnorm(N * J, sd = sqrt(1 - .64)), N, J)
  Q0 <- matrix(-1L, J, 2)
  for (k in 1:2) {
    anchor <- which(group == k)[1:2]
    Q0[anchor, ] <- 0L
    Q0[anchor, k] <- 1L
  }
  list(Y = Y, Q0 = Q0)
}

.toy <- .pefa_toy()
## One healthy converged sweep, reused by the shape tests below.
.healthy <- pefa(.toy$Q0, .toy$Y, 2, 4, v0 = .001, max_it = 800,
                 verbose = FALSE)

.sweep_names <- c(
  "K", "ELBO", "t_nom", "t", "t_S", "AIC", "BIC", "AIC_S", "BIC_S",
  "RMSEA", "SRMR", "CFI", "TLI", "iter", "secs", "converged",
  "backbone_ssl_min", "fit_status", "fit_reason", "loading_status",
  "loading_reason", "pip_status", "pip_reason", "elbo_status", "elbo_reason",
  "stats_status", "stats_reason"
)
.transition_names <- c(
  "K_from", "K_to", "ELBO_gain", "AIC_gain", "BIC_gain", "AIC_S_gain",
  "BIC_S_gain", "ELBO_gain_pct", "AIC_gain_pct", "BIC_gain_pct",
  "AIC_S_gain_pct", "BIC_S_gain_pct"
)
.persistence_names <- c(
  "K_from", "K_to", "r", "n_pairs", "phi_min", "rmsd", "rmsd_max", "ari",
  "pip_rmsd", "unmatched_n", "unmatched_ssl", "collision", "n_collisions",
  "collision_targets", "collision_multiplicities", "source_ineligible_n",
  "target_eligible_n", "backbone_degenerate", "from_fit_ok", "to_fit_ok",
  "pair_status", "pair_reason"
)

## A sweep skeleton for the table builders, which read only these columns.
.stub_sweep <- function(Ks, converged = TRUE, loading = "ok", pip = "ok") {
  data.frame(K = as.integer(Ks),
             converged = rep_len(converged, length(Ks)),
             loading_status = rep_len(loading, length(Ks)),
             pip_status = rep_len(pip, length(Ks)),
             stringsAsFactors = FALSE)
}
.keyed <- function(mats, Ks) stats::setNames(mats, as.character(Ks))
.unit <- function(i, J = 6) {
  v <- numeric(J)
  v[i] <- 1
  v
}

test_that("pefa returns exactly the 0.9.0 evidence object", {
  expect_identical(
    names(formals(pefa)),
    c("Q0", "Y", "Kmin", "Kmax", "bifactor", "general", "v0", "max_it",
      "convChk", "tolVal", "tau", "stability_eps", "rank_adjust",
      "rank_max_J", "verbose")
  )
  expect_false(any(c("cuts", "sustain", "...") %in% names(formals(pefa))))

  expect_identical(class(.healthy), "pefa")
  expect_identical(names(.healthy), c("sweep", "transitions", "persistence",
                                      "loadings", "pips", "Q0", "settings",
                                      "provenance"))
  for (gone in c("selected_K", "boundary", "call", "window", "gain_status",
                 "fits", "ssl", "selection", "selected_fit")) {
    expect_false(gone %in% names(.healthy))
  }
  expect_identical(names(.healthy$loadings), as.character(.healthy$sweep$K))
  expect_identical(names(.healthy$pips), names(.healthy$loadings))
  expect_identical(dim(.healthy$Q0), c(ncol(.toy$Y), 2L))
  expect_identical(typeof(.healthy$Q0), "integer")
  for (key in names(.healthy$loadings)) {
    expect_identical(dim(.healthy$loadings[[key]]),
                     c(ncol(.toy$Y), as.integer(key)))
    expect_identical(dimnames(.healthy$pips[[key]]),
                     dimnames(.healthy$loadings[[key]]))
  }
})

test_that("the three tables keep their documented columns and types", {
  expect_named(.healthy$sweep, .sweep_names)
  expect_named(.healthy$transitions, .transition_names)
  expect_named(.healthy$persistence, .persistence_names)

  healthy_types <- vapply(.healthy$persistence, typeof, "")
  expect_identical(
    healthy_types,
    vapply(vbpm:::.pefa_empty_persistence(), typeof, "")
  )

  ## An all-unavailable table keeps the same signature: pair_reason stays
  ## character rather than collapsing to logical NA.
  Ks <- 2:3
  loadings <- .keyed(list(cbind(.unit(1), .unit(2)),
                          cbind(.unit(1), .unit(2), .unit(3))), Ks)
  pips <- .keyed(list(matrix(.5, 6, 2), matrix(.5, 6, 3)), Ks)
  unavailable <- vbpm:::.pefa_persistence_table(
    .stub_sweep(Ks, loading = "malformed_loading"), loadings, pips,
    matrix(-1L, 6, 1), K0 = 1, stability_eps = 0.1, bifactor = FALSE
  )
  expect_identical(vapply(unavailable, typeof, ""), healthy_types)
  expect_identical(unavailable$pair_status, "unavailable")
  expect_identical(unavailable$pair_reason, "malformed_loading")

  ## A one-candidate sweep returns typed zero-row companion tables.
  single <- pefa(.toy$Q0, .toy$Y, 2, 2, v0 = .001, max_it = 400,
                 verbose = FALSE)
  expect_identical(nrow(single$transitions), 0L)
  expect_identical(nrow(single$persistence), 0L)
  expect_named(single$transitions, .transition_names)
  expect_named(single$persistence, .persistence_names)
  expect_identical(vapply(single$persistence, typeof, ""), healthy_types)
  expect_identical(vapply(single$transitions, typeof, ""),
                   vapply(.healthy$transitions, typeof, ""))
})

test_that("persistence holds every ordered pair, ordered, with r", {
  z <- .healthy$persistence
  W <- nrow(.healthy$sweep)
  expect_identical(nrow(z), as.integer(W * (W - 1L) / 2L))
  expect_identical(z$K_from, c(2L, 2L, 3L))
  expect_identical(z$K_to, c(3L, 4L, 4L))
  expect_identical(z$r, z$K_to - z$K_from)
  expect_false(is.unsorted(z$K_from))
  for (k in unique(z$K_from)) {
    expect_false(is.unsorted(z$K_to[z$K_from == k]))
  }
  expect_identical(nrow(.healthy$transitions), W - 1L)
})

test_that("an r >= 2 edge is the direct endpoint comparison", {
  ## Source column u2 picks the blend w in K = 3, and w in turn picks u3 in
  ## K = 4 -- so a composed assignment would score u2 against u3 (phi = 0),
  ## while the direct comparison scores it against u2 (phi = 1).
  w <- 0.6 * .unit(2) + 0.8 * .unit(3)
  L2 <- cbind(.unit(1), .unit(2))
  L3 <- cbind(.unit(1), w, .unit(4))
  L4 <- cbind(.unit(1), .unit(2), .unit(3), .unit(5))
  Ks <- 2:4
  table <- vbpm:::.pefa_persistence_table(
    .stub_sweep(Ks), .keyed(list(L2, L3, L4), Ks),
    .keyed(list(matrix(.5, 6, 2), matrix(.5, 6, 3), matrix(.5, 6, 4)), Ks),
    matrix(-1L, 6, 1), K0 = 1, stability_eps = 0.1, bifactor = FALSE
  )
  direct <- table[table$K_from == 2L & table$K_to == 4L, ]
  expect_identical(direct$r, 2L)
  expect_equal(direct$phi_min, 1, tolerance = 1e-12)
  expect_equal(direct$rmsd_max, 0, tolerance = 1e-12)

  ## The two adjacent legs, and hence their composition, disagree with it.
  expect_equal(table$phi_min[table$K_from == 2L & table$K_to == 3L], 0.6,
               tolerance = 1e-12)
  expect_equal(table$phi_min[table$K_from == 3L & table$K_to == 4L], 0,
               tolerance = 1e-12)
})

test_that("gains carry the documented orientation and percent guards", {
  sweep <- data.frame(K = 2:5,
                      ELBO = c(0, 100, 120, 125),
                      AIC = c(500, 400, 390, 388),
                      BIC = c(600, 500, 495, 494),
                      AIC_S = c(510, 405, 400, 399),
                      BIC_S = c(610, 505, 503, 502))
  tr <- vbpm:::.pefa_transitions_table(sweep)
  expect_named(tr, .transition_names)
  ## ELBO larger-is-better; every IC gain is from - to, so positive favors
  ## the larger candidate on all five paths.
  expect_equal(tr$ELBO_gain, c(100, 20, 5))
  expect_equal(tr$AIC_gain, c(100, 10, 2))
  expect_equal(tr$BIC_gain, c(100, 5, 1))
  expect_equal(tr$AIC_S_gain, c(105, 5, 1))
  expect_equal(tr$BIC_S_gain, c(105, 2, 1))
  expect_equal(tr$ELBO_gain_pct, c(100, 20, 5))
  expect_equal(tr$BIC_gain_pct, c(100, 5, 1))

  ## One nonfinite candidate blanks only the edges it touches.
  broken <- sweep
  broken$ELBO <- c(0, 100, NA, 125)
  tr_broken <- vbpm:::.pefa_transitions_table(broken)
  expect_equal(tr_broken$ELBO_gain, c(100, NA, NA))
  expect_equal(tr_broken$ELBO_gain_pct, c(100, NA, NA))
  expect_equal(tr_broken$AIC_gain_pct, c(100, 10, 2))

  ## Delta_max <= 0 leaves the raw gains and drops the percentages.
  falling <- sweep
  falling$ELBO <- c(0, -10, -20, -30)
  tr_falling <- vbpm:::.pefa_transitions_table(falling)
  expect_equal(tr_falling$ELBO_gain, c(-10, -10, -10))
  expect_true(all(is.na(tr_falling$ELBO_gain_pct)))
  expect_equal(tr_falling$AIC_gain_pct, c(100, 10, 2))
})

test_that("an all-nonfinite gain path is NA and emits no warning", {
  sweep <- data.frame(K = 2:5,
                      ELBO = rep(NA_real_, 4),
                      AIC = c(500, 400, 390, 388),
                      BIC = c(600, 500, 495, 494),
                      AIC_S = c(510, 405, 400, 399),
                      BIC_S = c(610, 505, 503, 502))
  ## max(numeric(0)) would return -Inf and warn; the guard must prevent both.
  expect_silent(tr <- vbpm:::.pefa_transitions_table(sweep))
  emitted <- character(0)
  withCallingHandlers(
    vbpm:::.pefa_transitions_table(sweep),
    warning = function(w) {
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(emitted, 0L)
  expect_identical(tr$ELBO_gain, rep(NA_real_, 3))
  expect_identical(tr$ELBO_gain_pct, rep(NA_real_, 3))
  expect_identical(typeof(tr$ELBO_gain_pct), "double")
  expect_false(any(is.infinite(tr$ELBO_gain_pct)))
  expect_equal(tr$AIC_gain_pct, c(100, 10, 2))
})

test_that("gain percentages survive gains near the largest double", {
  xmax <- .Machine$double.xmax
  ## Each AIC gain is IC_from - IC_to, so these four finite criterion values
  ## produce gains of xmax, xmax / 2, and -xmax.
  sweep <- data.frame(K = 2:5,
                      ELBO = c(0, 100, 120, 125),
                      AIC = c(xmax, 0, -xmax / 2, xmax / 2),
                      BIC = c(600, 500, 495, 494),
                      AIC_S = c(510, 405, 400, 399),
                      BIC_S = c(610, 505, 503, 502))
  tr <- vbpm:::.pefa_transitions_table(sweep)
  expect_equal(tr$AIC_gain, c(xmax, xmax / 2, -xmax))
  ## 100 * gain overflows to Inf before the division and would blank all
  ## three cells; 100 * (gain / gain_max) keeps every ratio representable.
  expect_equal(tr$AIC_gain_pct, c(100, 50, -100))
  expect_false(any(is.na(tr$AIC_gain_pct)))
  expect_false(any(is.infinite(tr$AIC_gain_pct)))
  ## The other paths are untouched by the extreme one.
  expect_equal(tr$BIC_gain_pct, c(100, 5, 1))

  ## Percentages are scale-free: any positive rescaling of the gain vector
  ## leaves them alone, at either end of the representable range.
  expect_equal(vbpm:::.pefa_gain_pct(c(xmax, xmax / 2, -xmax)),
               c(100, 50, -100))
  for (s in c(1e-300, 1e-3, 1, 7, 1e300)) {
    expect_equal(vbpm:::.pefa_gain_pct(c(1, 0.5, -1) * s),
                 c(100, 50, -100))
  }
})

test_that("component statuses are independent", {
  dn <- list(paste0("item_", 1:4), c("F1", "F2"))
  stats_ok <- stats::setNames(
    c(10, 9, 8.5, 100, 110, 101, 111, .05, .04, .98, .97),
    c("t_nom", "t", "t_S", "AIC", "BIC", "AIC_S", "BIC_S",
      "RMSEA", "SRMR", "CFI", "TLI")
  )
  fit <- list(flag = 1, Lam = matrix(0.5, 4, 2), pi = matrix(0.4, 4, 2),
              ELBO = -100, iter = 12)

  healthy <- vbpm:::.pefa_candidate(fit, stats_ok, 2L, dn, 1L, FALSE, 0.5)$row
  expect_identical(healthy$converged, TRUE)
  expect_identical(
    unlist(healthy[c("fit_status", "loading_status", "pip_status",
                     "elbo_status", "stats_status")], use.names = FALSE),
    rep("ok", 5)
  )
  expect_identical(healthy$K, 2L)
  expect_identical(healthy$iter, 12L)
  expect_identical(typeof(healthy$secs), "double")

  ## Nonconvergence plus an invalid PIP: two unavailable components, three
  ## untouched, and the finite loading matrix stored as returned.
  combined_fit <- fit
  combined_fit$flag <- 0
  combined_fit$pi <- matrix(1.5, 4, 2)
  combined <- vbpm:::.pefa_candidate(combined_fit, stats_ok, 2L, dn, 1L,
                                     FALSE, 0.5)
  expect_identical(combined$row$converged, FALSE)
  expect_identical(combined$row$fit_status, "unavailable")
  expect_identical(combined$row$fit_reason, "nonconverged")
  expect_identical(combined$row$pip_status, "unavailable")
  expect_identical(combined$row$pip_reason, "out_of_range_pip")
  expect_identical(combined$row$loading_status, "ok")
  expect_identical(combined$row$elbo_status, "ok")
  expect_identical(combined$row$stats_status, "ok")
  expect_identical(combined$row$ELBO, -100)
  expect_identical(combined$loading, matrix(0.5, 4, 2, dimnames = dn))
  expect_true(all(is.na(combined$pip)))
  expect_identical(dim(combined$pip), c(4L, 2L))

  ## A fit error gives every component a controlled reason and expected-shape
  ## typed-NA matrices, with converged = NA rather than FALSE.
  err <- tryCatch(stop("boom"), error = identity)
  failed <- vbpm:::.pefa_candidate(err, err, 2L, dn, 1L, FALSE, 0.5)
  expect_identical(failed$row$converged, NA)
  expect_identical(failed$row$fit_reason, "fit_error")
  expect_identical(
    unlist(failed$row[c("loading_reason", "pip_reason", "elbo_reason",
                        "stats_reason")], use.names = FALSE),
    rep("upstream_fit_error", 4)
  )
  expect_identical(dimnames(failed$loading), dn)
  expect_true(all(is.na(failed$loading)))
  expect_identical(failed$row$backbone_ssl_min, NA_real_)
  expect_identical(failed$row$iter, NA_integer_)

  ## A missing convergence flag is malformed_fit, and stays NA.
  flagless <- fit
  flagless$flag <- NULL
  malformed <- vbpm:::.pefa_candidate(flagless, stats_ok, 2L, dn, 1L,
                                      FALSE, 0.5)$row
  expect_identical(malformed$converged, NA)
  expect_identical(malformed$fit_reason, "malformed_fit")
  expect_identical(malformed$loading_status, "ok")

  ## A malformed loading leaves the criteria alone and blanks the backbone SSL.
  wrong_shape <- fit
  wrong_shape$Lam <- matrix(0.5, 3, 2)
  shaped <- vbpm:::.pefa_candidate(wrong_shape, stats_ok, 2L, dn, 1L,
                                   FALSE, 0.5)$row
  expect_identical(shaped$loading_reason, "malformed_loading")
  expect_identical(shaped$backbone_ssl_min, NA_real_)
  expect_identical(shaped$pip_status, "ok")
  expect_identical(shaped$AIC, 100)

  ## A nonfinite cell is normalized in place after the schema passes.
  spoiled <- fit
  spoiled$Lam[1, 1] <- Inf
  nonfinite <- vbpm:::.pefa_candidate(spoiled, stats_ok, 2L, dn, 1L,
                                      FALSE, 0.5)
  expect_identical(nonfinite$row$loading_reason, "nonfinite_loading")
  expect_identical(dim(nonfinite$loading), c(4L, 2L))
  expect_true(is.na(nonfinite$loading[1, 1]))
  expect_identical(nonfinite$loading[2, 1], 0.5)

  ## Legitimate descriptive NAs coexist with stats_status "ok"; an invalid
  ## criterion does not.
  descriptive <- stats_ok
  descriptive[["RMSEA"]] <- NA_real_
  descriptive[["TLI"]] <- NaN
  soft <- vbpm:::.pefa_candidate(fit, descriptive, 2L, dn, 1L, FALSE, 0.5)$row
  expect_identical(soft$stats_status, "ok")
  expect_identical(soft$stats_reason, NA_character_)
  expect_identical(soft$RMSEA, NA_real_)
  expect_identical(soft$TLI, NA_real_)

  criterion <- stats_ok
  criterion[["BIC"]] <- Inf
  hard <- vbpm:::.pefa_candidate(fit, criterion, 2L, dn, 1L, FALSE, 0.5)$row
  expect_identical(hard$stats_status, "unavailable")
  expect_identical(hard$stats_reason, "invalid_criterion_stats")
  expect_identical(hard$BIC, NA_real_)
  expect_identical(hard$AIC, 100)          # the valid scalars survive
  expect_identical(hard$elbo_status, "ok")
  expect_identical(hard$loading_status, "ok")
})

test_that("a nonconverged endpoint keeps every defined metric", {
  warnings <- character(0)
  object <- withCallingHandlers(
    pefa(.toy$Q0, .toy$Y, 2, 3, v0 = .001, max_it = 3, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  ## The fit component records the nonconvergence; nothing else is erased.
  expect_identical(object$sweep$converged, c(FALSE, FALSE))
  expect_identical(object$sweep$fit_status, rep("unavailable", 2))
  expect_identical(object$sweep$fit_reason, rep("nonconverged", 2))
  expect_identical(object$sweep$loading_status, rep("ok", 2))
  expect_identical(object$sweep$pip_status, rep("ok", 2))
  expect_identical(object$sweep$elbo_status, rep("ok", 2))
  expect_identical(object$sweep$stats_status, rep("ok", 2))
  expect_true(all(is.finite(object$sweep$ELBO)))
  expect_true(all(is.finite(object$sweep$backbone_ssl_min)))
  for (key in names(object$loadings)) {
    expect_true(all(is.finite(object$loadings[[key]])))
    expect_true(all(is.finite(object$pips[[key]])))
  }

  ## The pair is available and fully measured, with the fact recorded.
  z <- object$persistence
  expect_identical(nrow(z), 1L)
  expect_identical(z$pair_status, "available")
  expect_identical(z$pair_reason, NA_character_)
  expect_identical(z$from_fit_ok, FALSE)
  expect_identical(z$to_fit_ok, FALSE)
  expect_identical(z$n_pairs, 2L)
  for (metric in c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                   "unmatched_ssl")) {
    expect_true(is.finite(z[[metric]]),
                info = paste("metric", metric, "was blanked"))
  }
  expect_false(is.na(z$collision))
  expect_identical(z$n_collisions, 0L)

  ## The matcher agrees with a direct call on the stored matrices, so nothing
  ## was suppressed on the way into the table.
  direct <- vbpm:::.pefa_pair_facts(
    object$loadings[["2"]], object$loadings[["3"]],
    object$pips[["2"]], object$pips[["3"]],
    matrix(object$Q0 == -1L, nrow(object$Q0), 2L),
    cbind(matrix(object$Q0 == -1L, nrow(object$Q0), 2L),
          matrix(TRUE, nrow(object$Q0), 1L)),
    K0 = 2L, stability_eps = object$settings$stability_eps,
    loading_ok_from = TRUE, loading_ok_to = TRUE,
    pip_ok_from = TRUE, pip_ok_to = TRUE
  )
  expect_equal(z$phi_min, direct$phi_min, tolerance = 1e-12)
  expect_equal(z$rmsd, direct$rmsd, tolerance = 1e-12)
  expect_equal(z$pip_rmsd, direct$pip_rmsd, tolerance = 1e-12)

  ## Exactly one aggregate warning, naming the affected K values by code.
  expect_length(warnings, 1L)
  expect_match(warnings, "nonconverged")
  expect_match(warnings, "K = 2, 3")

  ## The triangle keeps the value and flags the endpoint.
  triangle <- persistence(object, "phi_min")
  expect_true(is.finite(triangle["2", "3"]))
  expect_false(attr(triangle, "fit_ok")["2", "3"])
})

test_that("a converged sweep is silent and warns at most once", {
  expect_silent(pefa(.toy$Q0, .toy$Y, 2, 3, v0 = .001, max_it = 800,
                     verbose = FALSE))
  expect_identical(.healthy$sweep$converged, rep(TRUE, 3))
  expect_identical(.healthy$sweep$fit_status, rep("ok", 3))

  ## Three nonconverged candidates warn once, not three times: the
  ## per-candidate MFVB iteration-limit warning stays muffled.
  warnings <- character(0)
  withCallingHandlers(
    pefa(.toy$Q0, .toy$Y, 2, 4, v0 = .001, max_it = 3, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_false(any(grepl("MFVB", warnings, fixed = TRUE)))
})

test_that("settings records all ten controls and omits verbose", {
  expect_identical(
    names(.healthy$settings),
    c("bifactor", "general", "v0", "max_it", "convChk", "tolVal", "tau",
      "stability_eps", "rank_adjust", "rank_max_J")
  )
  expect_length(.healthy$settings, 10L)
  expect_false("verbose" %in% names(.healthy$settings))
  ## Every other formal of pefa() except the data arguments is recorded.
  recorded <- c(names(.healthy$settings), "Q0", "Y", "Kmin", "Kmax",
                "verbose")
  expect_true(all(names(formals(pefa)) %in% recorded))
  expect_identical(.healthy$settings$stability_eps, 0.1)
  expect_identical(.healthy$settings$v0, 0.001)
  expect_identical(.healthy$settings$max_it, 800L)
  expect_false(.healthy$settings$bifactor)
  ## K0 is recoverable from the object rather than stored twice.
  expect_identical(ncol(.healthy$Q0), 2L)
})

test_that("pefa validates the window and the stability screen", {
  Q0 <- .toy$Q0
  Y <- .toy$Y
  for (bad in list(0, -1, Inf, NA_real_, c(0.1, 0.2), "0.1")) {
    expect_error(
      pefa(Q0, Y, 2, 2, stability_eps = bad, verbose = FALSE),
      "stability_eps"
    )
  }
  expect_error(pefa(Q0, Y, 0, 2, verbose = FALSE), "Kmin")
  expect_error(pefa(Q0, Y, 2.5, 3, verbose = FALSE), "Kmin")
  expect_error(pefa(Q0, Y, 3, 2, verbose = FALSE), "must be non-decreasing")
  expect_error(pefa(Q0, Y, 1, 3, verbose = FALSE),
               "at least the backbone size")
  ## A fully exploratory backbone still starts at K = 1, never at zero.
  expect_error(
    pefa(matrix(integer(0), ncol(Y), 0), Y, 0, 2, verbose = FALSE),
    "Kmin"
  )
})

## Runs `expr`, returning the error message (or the value) plus every warning
## the call emitted, so a test can assert "controlled error and no warning".
.quiet_error <- function(expr) {
  emitted <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = conditionMessage),
    warning = function(w) {
      emitted <<- c(emitted, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(message = value, warnings = emitted)
}

test_that("integer-backed formals are range-checked before coercion", {
  Q0 <- .toy$Q0
  Y <- .toy$Y
  ## as.integer(3e10) is NA_integer_ plus a coercion warning, which would
  ## turn an argument mistake into a silent failed-fit evidence object.
  big <- 3e10
  oversized <- list(
    Kmin = function() pefa(Q0, Y, big, big + 1, verbose = FALSE),
    Kmax = function() pefa(Q0, Y, 2, big, verbose = FALSE),
    max_it = function() pefa(Q0, Y, 2, 2, max_it = big, verbose = FALSE),
    rank_max_J = function() pefa(Q0, Y, 2, 2, rank_max_J = big,
                                 verbose = FALSE)
  )
  for (nm in names(oversized)) {
    got <- .quiet_error(oversized[[nm]]())
    expect_true(is.character(got$message), info = nm)
    expect_match(got$message, nm, fixed = TRUE)
    expect_match(got$message, "whole number in", fixed = TRUE)
    expect_length(got$warnings, 0L)
  }

  ## Complex values are rejected outright rather than coerced.
  complexes <- list(
    Kmin = function() pefa(Q0, Y, complex(real = 2), 3, verbose = FALSE),
    Kmax = function() pefa(Q0, Y, 2, complex(real = 3), verbose = FALSE),
    max_it = function() pefa(Q0, Y, 2, 2, verbose = FALSE,
                             max_it = complex(real = 400, imaginary = 1)),
    rank_max_J = function() pefa(Q0, Y, 2, 2, verbose = FALSE,
                                 rank_max_J = complex(real = 100))
  )
  for (nm in names(complexes)) {
    got <- .quiet_error(complexes[[nm]]())
    expect_true(is.character(got$message), info = nm)
    expect_match(got$message, nm, fixed = TRUE)
    expect_length(got$warnings, 0L)
  }

  ## The range check is inclusive, so the largest representable integer is
  ## still a legal setting and reaches $settings unchanged.
  at_max <- pefa(Q0, Y, 2, 2, v0 = .001, max_it = 400,
                 rank_max_J = .Machine$integer.max, verbose = FALSE)
  expect_identical(at_max$settings$rank_max_J, .Machine$integer.max)
})

test_that("rank_max_J normalizes so 100 and 100L hash identically", {
  expect_identical(.healthy$settings$rank_max_J, 100L)
  expect_identical(typeof(.healthy$settings$rank_max_J), "integer")
  expect_null(names(.healthy$settings$rank_max_J))

  ## The default 100 is a double; 100L must not fork the provenance.
  typed <- pefa(.toy$Q0, .toy$Y, 2, 4, v0 = .001, max_it = 800,
                rank_max_J = 100L, verbose = FALSE)
  expect_identical(typed$settings, .healthy$settings)
  expect_identical(typed$provenance$settings_sha256,
                   .healthy$provenance$settings_sha256)
  expect_identical(typed$provenance$evidence_sha256,
                   .healthy$provenance$evidence_sha256)
  expect_identical(typed$provenance$lineage_sha256,
                   .healthy$provenance$lineage_sha256)

  ## Names on a stored control are input decoration, not evidence: they must
  ## not move any hash either.
  named <- pefa(.toy$Q0, .toy$Y, 2, 4, v0 = c(spike = .001),
                max_it = c(it = 800), convChk = c(cc = FALSE),
                tolVal = c(tv = 1e-4), tau = c(th = 0.5),
                stability_eps = c(se = 0.10), bifactor = c(bf = FALSE),
                rank_adjust = c(ra = FALSE), rank_max_J = c(cap = 100),
                verbose = FALSE)
  expect_identical(vapply(named$settings[-2L], function(v) is.null(names(v)),
                          logical(1)),
                   stats::setNames(rep(TRUE, 9L),
                                   names(named$settings)[-2L]))
  expect_identical(named$settings, .healthy$settings)
  expect_identical(named$provenance$settings_sha256,
                   .healthy$provenance$settings_sha256)
  expect_identical(named$provenance$evidence_sha256,
                   .healthy$provenance$evidence_sha256)
})

test_that("backbone_ssl_min tracks ssl() in both modes", {
  K0 <- ncol(.healthy$Q0)
  expected <- vapply(ssl(.healthy), function(v) min(v[seq_len(K0)]),
                     numeric(1))
  expect_equal(.healthy$sweep$backbone_ssl_min, unname(expected),
               tolerance = 1e-12)

  bifactor <- pefa(.toy$Q0, .toy$Y, 2, 3, bifactor = TRUE, v0 = .001,
                   max_it = 400, verbose = FALSE)
  ## Stored columns 2:(K0 + 1), never the general column.
  expected_bf <- vapply(ssl(bifactor), function(v) min(v[seq_len(K0) + 1L]),
                        numeric(1))
  expect_equal(bifactor$sweep$backbone_ssl_min, unname(expected_bf),
               tolerance = 1e-12)

  ## A general column smaller than every backbone column separates the two
  ## offsets, so a leak would be visible.
  lam <- cbind(c(0.1, 0, 0), c(1, 0, 0), c(0, 2, 0))
  expect_identical(
    vbpm:::.pefa_backbone_ssl_min(lam, K0 = 2L, bifactor = TRUE,
                                  loading_ok = TRUE),
    1
  )
  expect_equal(
    vbpm:::.pefa_backbone_ssl_min(lam, K0 = 2L, bifactor = FALSE,
                                  loading_ok = TRUE),
    0.01, tolerance = 1e-12
  )
  expect_identical(
    vbpm:::.pefa_backbone_ssl_min(lam, K0 = 2L, bifactor = TRUE,
                                  loading_ok = FALSE),
    NA_real_
  )

  ## K0 = 0 has no backbone to measure.
  free <- pefa(matrix(integer(0), ncol(.toy$Y), 0), .toy$Y, 1, 1, v0 = .001,
               max_it = 400, verbose = FALSE)
  expect_identical(free$sweep$backbone_ssl_min, NA_real_)
  expect_identical(dim(free$Q0), c(ncol(.toy$Y), 0L))
})

test_that("bifactor stores K + 1 columns and matches the group block", {
  object <- pefa(.toy$Q0, .toy$Y, 2, 3, bifactor = TRUE, v0 = .001,
                 max_it = 400, verbose = FALSE)
  expect_true(object$settings$bifactor)
  expect_identical(object$sweep$K, 2:3)
  expect_identical(unname(vapply(object$loadings, ncol, integer(1))), 3:4)
  expect_identical(unname(vapply(object$pips, ncol, integer(1))), 3:4)
  expect_identical(colnames(object$loadings[["3"]]),
                   c("G", "F1", "F2", "F3"))
  expect_identical(unname(lengths(ssl(object))), c(3L, 4L))
  ## The matcher works on the group block, so a k = K0 edge scores K0 pairs.
  expect_identical(object$persistence$n_pairs, 2L)
  expect_identical(object$persistence$pair_status, "available")
  expect_true(is.finite(object$persistence$phi_min))
})

test_that("bifactor collision targets are group-block indices", {
  ## The reused target sits at group position 2, i.e. stored position 3.
  general <- rep(0.5, 6)
  stored_from <- cbind(general, .unit(1), .unit(2), .unit(2))
  stored_to <- cbind(general, .unit(1), .unit(2), .unit(5), .unit(6))
  Ks <- 3:4
  table <- vbpm:::.pefa_persistence_table(
    .stub_sweep(Ks), .keyed(list(stored_from, stored_to), Ks),
    .keyed(list(matrix(.5, 6, 4), matrix(.5, 6, 5)), Ks),
    matrix(-1L, 6, 1), K0 = 1, stability_eps = 0.1, bifactor = TRUE
  )
  expect_identical(table$collision, TRUE)
  expect_identical(table$n_collisions, 1L)
  expect_identical(table$collision_targets, "2")   # never "3"
  expect_identical(table$collision_multiplicities, "2")
  expect_identical(table$n_pairs, 3L)
  expect_identical(table$unmatched_n, 2L)
})

test_that("bifactor validation runs on the stored K + 1 matrix", {
  ## The candidate checker is handed the stored dimnames, general column
  ## included, so a group-block-sized matrix is malformed rather than valid.
  dn <- list(paste0("item_", 1:4), c("G", "F1", "F2"))
  stats_ok <- stats::setNames(
    c(10, 9, 8.5, 100, 110, 101, 111, .05, .04, .98, .97),
    c("t_nom", "t", "t_S", "AIC", "BIC", "AIC_S", "BIC_S",
      "RMSEA", "SRMR", "CFI", "TLI")
  )
  stored <- list(flag = 1, Lam = cbind(c(.3, .3, .3, .3), c(1, 1, 0, 0),
                                       c(0, 0, 2, 2)),
                 pi = matrix(0.4, 4, 3), ELBO = -100, iter = 7)
  ok <- vbpm:::.pefa_candidate(stored, stats_ok, 2L, dn, 2L, TRUE, 0.1)
  expect_identical(ok$row$loading_status, "ok")
  expect_identical(dimnames(ok$loading), dn)
  ## Backbone SSL skips the general column: min(2, 8), not min(0.36, 2).
  expect_equal(ok$row$backbone_ssl_min, 2, tolerance = 1e-12)

  group_only <- stored
  group_only$Lam <- stored$Lam[, -1L, drop = FALSE]
  stripped <- vbpm:::.pefa_candidate(group_only, stats_ok, 2L, dn, 2L, TRUE,
                                     0.1)
  expect_identical(stripped$row$loading_status, "unavailable")
  expect_identical(stripped$row$loading_reason, "malformed_loading")
  expect_identical(dim(stripped$loading), c(4L, 3L))
})

test_that("an unavailable pair still carries its endpoint facts", {
  Ks <- 2:3
  loadings <- .keyed(list(cbind(.unit(1), .unit(2)),
                          cbind(.unit(1), .unit(2), .unit(3))), Ks)
  pips <- .keyed(list(matrix(.5, 6, 2), matrix(.5, 6, 3)), Ks)
  sweep <- .stub_sweep(Ks, converged = c(FALSE, TRUE))
  sweep$loading_status <- c("malformed_loading", "ok")
  row <- vbpm:::.pefa_persistence_table(
    sweep, loadings, pips, matrix(-1L, 6, 1), K0 = 1, stability_eps = 0.1,
    bifactor = FALSE
  )
  expect_identical(row$pair_status, "unavailable")
  expect_identical(row$pair_reason, "malformed_loading")
  ## Convergence does not depend on the loadings, so it is still reported.
  expect_identical(row$from_fit_ok, FALSE)
  expect_identical(row$to_fit_ok, TRUE)
  ## The target endpoint survived, so its eligibility count is available.
  expect_identical(row$target_eligible_n, 2L)
  expect_identical(row$source_ineligible_n, NA_integer_)
  expect_identical(row$K_from, 2L)
  expect_identical(row$r, 1L)
})

test_that("a decision profile composes from the stored facts alone", {
  ## Four edges: a nonconverged one, a collided one, one that fails the
  ## congruence gate, and a clean one.  The package grades none of them.
  z <- data.frame(
    K_from = c(3L, 3L, 3L, 3L), K_to = 4:7,
    phi_min = c(0.95, 0.93, 0.70, 0.97),
    collision = c(FALSE, TRUE, FALSE, FALSE),
    source_ineligible_n = rep(0L, 4L), target_eligible_n = rep(1L, 4L),
    from_fit_ok = rep(TRUE, 4L), to_fit_ok = c(FALSE, TRUE, TRUE, TRUE),
    pair_status = rep("available", 4L),
    stringsAsFactors = FALSE
  )
  K0 <- 2L
  phi_star <- 0.90

  profile <- function(z, endpoint_fit_policy, collision_veto) {
    endpoint_admissible <- switch(
      endpoint_fit_policy,
      require_converged = (z$from_fit_ok %in% TRUE) & (z$to_fit_ok %in% TRUE),
      allow_completed = !is.na(z$from_fit_ok) & !is.na(z$to_fit_ok)
    )
    known <- z$pair_status == "available" & is.finite(z$phi_min) &
      endpoint_admissible & !is.na(z$source_ineligible_n) &
      !is.na(z$target_eligible_n) & !is.na(z$collision)
    pass <- rep(NA, nrow(z))
    pass[known] <- z$source_ineligible_n[known] == 0L &
      z$target_eligible_n[known] >= (z$K_from[known] - K0) &
      (!collision_veto | !z$collision[known]) &
      z$phi_min[known] >= phi_star
    pass
  }

  expect_identical(profile(z, "require_converged", TRUE),
                   c(NA, FALSE, FALSE, TRUE))
  ## Only the endpoint policy changes: the nonconverged edge is now read from
  ## its retained metric.
  expect_identical(profile(z, "allow_completed", TRUE),
                   c(TRUE, FALSE, FALSE, TRUE))
  ## Only the collision policy changes: the collided edge is now evaluated.
  expect_identical(profile(z, "require_converged", FALSE),
                   c(NA, TRUE, FALSE, TRUE))

  ## An unavailable row and an available row with an NA metric are both
  ## unknown; a met-input condition failure is FALSE.
  z2 <- rbind(z, z[1L, ], z[1L, ])
  z2$to_fit_ok[5:6] <- TRUE
  z2$pair_status[5L] <- "unavailable"
  z2$phi_min[5L] <- NA_real_
  z2$phi_min[6L] <- NA_real_
  states <- profile(z2, "require_converged", TRUE)
  expect_identical(states[5:6], c(NA, NA))
  ## One failed required edge resolves the source as false even beside an NA;
  ## with no failure and at least one NA it stays unknown.
  resolve <- function(v) if (any(v %in% FALSE)) FALSE else if (anyNA(v)) NA
    else TRUE
  expect_identical(resolve(states), FALSE)
  expect_identical(resolve(states[c(1L, 4L)]), NA)
  expect_identical(resolve(states[4L]), TRUE)
})

test_that("print and summary announce no selected count", {
  printed <- capture.output(print(.healthy))
  expect_match(paste(printed, collapse = "\n"), "No factor count is selected")
  expect_match(paste(printed, collapse = "\n"), "select_K_elbow")
  expect_false(any(grepl("selected_K|boundary", printed)))
  expect_match(paste(printed, collapse = "\n"), "persistence 3 x 22")

  s <- summary(.healthy)
  expect_s3_class(s, "summary.pefa")
  expect_identical(s$window, c(Kmin = 2L, Kmax = 4L, K0 = 2L))
  expect_identical(s$persistence, .healthy$persistence)
  shown <- capture.output(print(s))
  expect_match(paste(shown, collapse = "\n"), "Candidates")
  expect_match(paste(shown, collapse = "\n"), "SS loadings")
  ## The persistence table is deferred to a row/column count.
  expect_match(paste(shown, collapse = "\n"),
               "Persistence: 3 ordered pairs x 22 fields")
  expect_false(any(grepl("pair_status", shown)))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(.healthy, type = "gain", criterion = "bic_s"))
  expect_invisible(plot(.healthy, type = "gain", criterion = "aic",
                        pct = TRUE))
  expect_invisible(plot(.healthy, type = "objective", criterion = "elbo"))
  expect_invisible(plot(.healthy, type = "fit"))
  expect_error(plot(.healthy, type = "gain", pct = NA), "pct")
})
