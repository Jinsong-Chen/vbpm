.transition_names_082 <- c(
  "K_from", "K_to", "ELBO_gain", "BIC_gain",
  "rmsd", "rmsd_max", "phi_min", "ari", "pip_sum", "pip_product",
  "max_unmatched_loading"
)

.make_sweep_082 <- function(loadings, pips, converged = NULL) {
  n <- length(loadings)
  if (is.null(converged)) converged <- rep(TRUE, n)
  out <- data.frame(
    K = seq_len(n) + 1L,
    ELBO = cumsum(c(100, rep(5, n - 1L))),
    BIC = cumsum(c(300, rep(-20, n - 1L))),
    converged = converged
  )
  out$loading <- I(loadings)
  out$pip <- I(pips)
  out
}

test_that("the 0.8.2 matcher fixes the backbone and deterministically matches free columns", {
  b <- c(.8, .8, .1, .1, .1, .1)
  f1 <- c(.1, .1, .8, .8, .1, .1)
  f2 <- c(.1, .1, .1, .1, .8, .8)
  tiny <- rep(.02, 6)
  lam1 <- cbind(backbone = b, first = f1, second = f2)
  lam2 <- cbind(backbone = -b, second = -f2, tiny = tiny, first = -f1)

  metrics <- vbpm:::.match_sweep_loadings(lam1, lam2, K0 = 1)
  expect_named(metrics, c(
    "rmsd", "rmsd_max", "phi_min", "ari", "pip_sum", "pip_product",
    "max_unmatched_loading"
  ))
  expect_equal(
    unname(metrics[c("rmsd", "rmsd_max", "phi_min", "ari")]),
    c(0, 0, 1, 1),
    tolerance = 1e-12
  )
  expect_true(all(is.na(metrics[c("pip_sum", "pip_product")])))
  expect_equal(unname(metrics["max_unmatched_loading"]), .02,
               tolerance = 1e-12)

  backbone <- c(.8, .8, 0, 0, 0, 0)
  free <- c(0, 0, .8, .8, 0, 0)
  decoy <- c(0, 0, 0, 0, .8, .8)
  fixed <- vbpm:::.match_sweep_loadings(
    cbind(backbone, free), cbind(decoy, free, backbone), K0 = 1
  )
  expect_identical(unname(fixed["phi_min"]), 0)
  expect_gt(fixed["rmsd"], 0)
})

test_that("free-column matching is injective, greedy, and has a stable tie rule", {
  x <- c(.8, .8, 0, 0)
  y <- c(0, 0, .8, .8)
  injective <- vbpm:::.match_sweep_loadings(cbind(x, x), cbind(x, y), K0 = 0)
  expect_identical(unname(injective["phi_min"]), 0)
  expect_gt(injective["rmsd"], 0)

  smaller <- cbind(c(.9, 0), c(.5, .5))
  larger <- cbind(c(1, 0), c(0, 2))
  forward <- vbpm:::.match_sweep_loadings(smaller, larger, K0 = 0)
  reverse <- vbpm:::.match_sweep_loadings(smaller[, 2:1], larger, K0 = 0)
  expect_lt(forward["rmsd"], reverse["rmsd"])
  expect_equal(unname(forward["rmsd"]), sqrt(2.51 / 4), tolerance = 1e-12)

  reference <- cbind(c(2, 0))
  candidate_a <- c(2, 1)
  candidate_b <- c(1, 0)
  first <- vbpm:::.match_sweep_loadings(
    reference, cbind(candidate_a, candidate_b), K0 = 0
  )
  reversed <- vbpm:::.match_sweep_loadings(
    reference, cbind(candidate_b, candidate_a), K0 = 0
  )
  expect_lt(first["phi_min"], reversed["phi_min"])
  expect_equal(unname(reversed["phi_min"]), 1, tolerance = 1e-12)
})

test_that("RMSD, ARI, screening, and unmatched size retain their distinct meanings", {
  core <- cbind(
    c(.8, .8, .8, 0, 0, 0),
    c(0, 0, 0, .8, .8, .8)
  )
  surplus <- c(0, 0, .9, .9, 0, 0)
  metrics <- vbpm:::.match_sweep_loadings(core, cbind(core, surplus), K0 = 1)
  expect_identical(unname(metrics["rmsd"]), 0)
  expect_identical(unname(metrics["phi_min"]), 1)
  expect_lt(metrics["ari"], 1)
  expect_equal(unname(metrics["max_unmatched_loading"]), .9)

  changed <- core
  changed[, 2] <- changed[, 2] + .4
  column_change <- vbpm:::.match_sweep_loadings(core, changed, K0 = 2)
  expect_gt(column_change["rmsd_max"], column_change["rmsd"])
  expect_equal(unname(column_change["rmsd_max"]), .4, tolerance = 1e-12)

  edge <- c(sqrt(.1), 0, 0, 0)
  retained <- vbpm:::.match_sweep_loadings(cbind(edge), cbind(edge),
                                           K0 = 0, eps = .1)
  excluded <- vbpm:::.match_sweep_loadings(cbind(edge / 2), cbind(edge / 2),
                                           K0 = 0, eps = .1)
  expect_identical(unname(retained["phi_min"]), 1)
  expect_true(is.na(excluded["rmsd"]) && is.na(excluded["phi_min"]))

  set.seed(1)
  rounding_edge <- rnorm(7)
  rounding_edge <- rounding_edge * sqrt(.1 / sum(rounding_edge^2))
  rounding_cutoff <- vbpm:::.match_sweep_loadings(
    cbind(rounding_edge), cbind(rounding_edge), K0 = 0, eps = .1
  )
  expect_gte(sum(rounding_edge^2), .1)
  expect_equal(unname(rounding_cutoff["phi_min"]), 1, tolerance = 1e-12)

  tiny <- rep(.02, 6)
  not_enough <- vbpm:::.match_sweep_loadings(
    cbind(core[, 1], core[, 2]), cbind(core[, 1], tiny), K0 = 0, eps = .1
  )
  expect_true(all(is.na(not_enough[c(
    "rmsd", "rmsd_max", "phi_min", "max_unmatched_loading"
  )])))
  expect_true(is.finite(not_enough["ari"]))
})

test_that("PIP reversals use matched cells and the inclusive .5 boundary", {
  lam1 <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  lam2 <- cbind(-lam1[, 1], rep(.01, 4), -lam1[, 2])
  pip1 <- matrix(.1, 4, 2)
  pip1[1, 1] <- .9
  pip1[3, 2] <- .5
  pip2 <- matrix(.9, 4, 3)
  pip2[, 1] <- pip1[, 1]
  pip2[, 3] <- pip1[, 2]
  pip2[1, 1] <- .1
  pip2[2, 1] <- .9

  metrics <- vbpm:::.match_sweep_loadings(
    lam1, lam2, K0 = 1, pip1 = pip1, pip2 = pip2
  )
  expect_equal(
    unname(metrics[c("pip_sum", "pip_product")]),
    c(2 / 8, 1 / 64),
    tolerance = 1e-12
  )
  expect_lte(metrics["pip_product"], metrics["pip_sum"]^2 / 4)

  boundary <- vbpm:::.match_sweep_loadings(
    cbind(c(.8, .8)), cbind(c(.8, .8)), K0 = 1,
    pip1 = cbind(c(.5, .49)), pip2 = cbind(c(.49, .5))
  )
  expect_equal(unname(boundary[c("pip_sum", "pip_product")]), c(1, .25))

  one_way <- vbpm:::.match_sweep_loadings(
    cbind(c(.8, .8)), cbind(c(.8, .8)), K0 = 1,
    pip1 = cbind(c(.9, .1)), pip2 = cbind(c(.1, .1))
  )
  expect_equal(unname(one_way[c("pip_sum", "pip_product")]), c(.5, 0))
})

test_that("orientation and scaled arithmetic preserve reviewed edge behavior", {
  zero_mean <- vbpm:::.match_sweep_loadings(
    cbind(c(1, -1)), cbind(c(-1, 1)), K0 = 1
  )
  expect_equal(unname(zero_mean["phi_min"]), -1, tolerance = 1e-12)

  huge <- cbind(c(1e200, 0))
  candidates <- cbind(c(0, 1e200), c(9e199, 0))
  huge_match <- vbpm:::.match_sweep_loadings(huge, candidates, K0 = 0)
  expect_equal(unname(huge_match["phi_min"]), 1, tolerance = 1e-12)
  expect_true(is.finite(huge_match["rmsd"]))

  expect_error(
    vbpm:::.match_sweep_loadings(matrix(1, 2, 2), matrix(1, 2, 1), K0 = 0),
    "no greater"
  )
  expect_error(
    vbpm:::.match_sweep_loadings(huge, huge, K0 = 0, eps = 0),
    "positive"
  )
  expect_error(
    vbpm:::.match_sweep_loadings(huge, huge, K0 = 0,
                                 pip1 = matrix(.5, 2, 1)),
    "both"
  )
})

test_that("transition tables have the exact 0.8.2 schema and typed empty columns", {
  empty <- vbpm:::.transition_table()
  expect_s3_class(empty, "data.frame")
  expect_named(empty, .transition_names_082)
  expect_identical(nrow(empty), 0L)
  expect_type(empty$K_from, "integer")
  expect_type(empty$K_to, "integer")
  expect_true(all(vapply(empty[-(1:2)], is.double, logical(1))))

  one <- matrix(c(.8, .8, 0, 0), 4, 1)
  sweep <- .make_sweep_082(list(one), list(matrix(.5, 4, 1)))
  built <- vbpm:::.build_transitions(sweep, K0 = 1, stability_eps = .1)
  expect_named(built, .transition_names_082)
  expect_identical(nrow(built), 0L)
})

test_that("the builder preserves gains while gating bifactor and failed endpoints", {
  first <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  second <- cbind(first, c(.7, 0, 0, .7))
  pips <- list(matrix(.5, 4, 2), matrix(.5, 4, 3))
  sweep <- .make_sweep_082(list(first, second), pips)

  ordinary <- vbpm:::.build_transitions(sweep, K0 = 1,
                                        stability_eps = .1)
  expect_named(ordinary, .transition_names_082)
  expect_identical(ordinary$K_from, 2L)
  expect_identical(ordinary$K_to, 3L)
  expect_equal(ordinary$ELBO_gain, 5)
  expect_equal(ordinary$BIC_gain, 20)
  expect_true(all(is.finite(unlist(ordinary[c(
    "rmsd", "rmsd_max", "phi_min", "pip_sum", "pip_product",
    "max_unmatched_loading"
  )], use.names = FALSE))))

  failed <- sweep
  failed$converged[2] <- FALSE
  expect_warning(
    failed_edge <- vbpm:::.build_transitions(
      failed, K0 = 1, stability_eps = .1
    ),
    NA
  )
  expect_equal(failed_edge$ELBO_gain, 5)
  expect_equal(failed_edge$BIC_gain, 20)
  expect_true(all(is.na(failed_edge[.transition_names_082[-(1:4)]])))

  expect_warning(
    bifactor_edge <- vbpm:::.build_transitions(
      sweep, K0 = 1, stability_eps = .1, bifactor = TRUE
    ),
    NA
  )
  expect_equal(bifactor_edge$ELBO_gain, 5)
  expect_equal(bifactor_edge$BIC_gain, 20)
  expect_true(all(is.na(bifactor_edge[.transition_names_082[-(1:4)]])))
})

test_that("malformed converged endpoints warn once without losing finite gains", {
  first <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  second <- cbind(first, c(.7, 0, 0, .7))
  third <- cbind(second, c(0, .7, .7, 0))
  malformed <- second
  malformed[1, 1] <- Inf
  loadings <- list(first, malformed, third)
  pips <- lapply(loadings, function(x) matrix(.5, nrow(x), ncol(x)))
  sweep <- .make_sweep_082(loadings, pips)

  expect_error(
    vbpm:::.match_sweep_loadings(
      first, malformed, K0 = 1, eps = .1,
      pip1 = pips[[1]], pip2 = pips[[2]]
    ),
    "finite numeric matrices"
  )

  warnings <- character(0)
  transitions <- withCallingHandlers(
    vbpm:::.build_transitions(sweep, K0 = 1, stability_eps = .1),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_match(warnings, "malformed or nonfinite")
  expect_match(warnings, "2->3, 3->4", fixed = TRUE)
  expect_equal(transitions$ELBO_gain, c(5, 5))
  expect_equal(transitions$BIC_gain, c(20, 20))
  expect_true(all(is.na(transitions[.transition_names_082[-(1:4)]])))
})
