.transition_names_083 <- c(
  "K_from", "K_to", "ELBO_gain", "BIC_gain",
  "rmsd", "rmsd_max", "phi_min", "ari", "pip_rmsd", "collision",
  "unmatched_max"
)

.make_sweep_083 <- function(loadings, pips, converged = NULL) {
  n <- length(loadings)
  if (is.null(converged)) converged <- rep(TRUE, n)
  sweep <- data.frame(
    K = seq_len(n) + 1L,
    ELBO = cumsum(c(100, rep(5, n - 1L))),
    BIC = cumsum(c(300, rep(-20, n - 1L))),
    converged = converged
  )
  keys <- as.character(sweep$K)
  list(
    sweep = sweep,
    loadings = stats::setNames(loadings, keys),
    pips = stats::setNames(pips, keys)
  )
}

test_that("the 0.8.3 matcher fixes the backbone and deterministically matches free columns", {
  b <- c(.8, .8, .1, .1, .1, .1)
  f1 <- c(.1, .1, .8, .8, .1, .1)
  f2 <- c(.1, .1, .1, .1, .8, .8)
  tiny <- rep(.02, 6)
  lam1 <- cbind(backbone = b, first = f1, second = f2)
  lam2 <- cbind(backbone = -b, second = -f2, tiny = tiny, first = -f1)

  metrics <- vbpm:::.match_sweep_loadings(lam1, lam2, K0 = 1)
  expect_named(metrics, c(
    "rmsd", "rmsd_max", "phi_min", "ari", "pip_rmsd", "collision",
    "unmatched_max"
  ))
  expect_equal(
    unname(metrics[c("rmsd", "rmsd_max", "phi_min", "ari")]),
    c(0, 0, 1, 1),
    tolerance = 1e-12
  )
  expect_true(is.na(metrics["pip_rmsd"]))
  expect_equal(unname(metrics["unmatched_max"]), .02,
               tolerance = 1e-12)
  expect_false(as.logical(metrics["collision"]))

  backbone <- c(.8, .8, 0, 0, 0, 0)
  free <- c(0, 0, .8, .8, 0, 0)
  decoy <- c(0, 0, 0, 0, .8, .8)
  fixed <- vbpm:::.match_sweep_loadings(
    cbind(backbone, free), cbind(decoy, free, backbone), K0 = 1
  )
  expect_identical(unname(fixed["phi_min"]), 0)
  expect_gt(fixed["rmsd"], 0)
})

test_that("free-column matching is independent and has a stable tie rule", {
  x <- c(.8, .8, 0, 0)
  y <- c(0, 0, .8, .8)
  reused <- vbpm:::.match_sweep_loadings(cbind(x, x), cbind(x, y), K0 = 0)
  expect_equal(unname(reused["phi_min"]), 1, tolerance = 1e-12)
  expect_identical(unname(reused["rmsd"]), 0)
  expect_true(as.logical(reused["collision"]))

  distinct <- vbpm:::.match_sweep_loadings(cbind(x, y), cbind(y, x), K0 = 0)
  expect_false(as.logical(distinct["collision"]))

  smaller <- cbind(c(.9, 0), c(.5, .5))
  larger <- cbind(c(1, 0), c(0, 2))
  forward <- vbpm:::.match_sweep_loadings(smaller, larger, K0 = 0)
  reverse <- vbpm:::.match_sweep_loadings(smaller[, 2:1], larger, K0 = 0)
  expect_equal(unname(forward["rmsd"]), unname(reverse["rmsd"]),
               tolerance = 1e-12)
  expect_equal(unname(forward["rmsd"]), sqrt(.51 / 4), tolerance = 1e-12)
  expect_true(as.logical(forward["collision"]))

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
  expect_equal(unname(metrics["unmatched_max"]), .9)

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
  reused_candidate <- vbpm:::.match_sweep_loadings(
    cbind(core[, 1], core[, 2]), cbind(core[, 1], tiny), K0 = 0, eps = .1
  )
  expect_true(all(is.finite(reused_candidate[c(
    "rmsd", "rmsd_max", "phi_min", "unmatched_max"
  )])))
  expect_true(as.logical(reused_candidate["collision"]))

  unavailable <- vbpm:::.match_sweep_loadings(
    cbind(core[, 1], core[, 2]), cbind(tiny, tiny), K0 = 0, eps = .1
  )
  expect_true(all(is.na(unavailable[c(
    "rmsd", "rmsd_max", "phi_min", "unmatched_max", "collision"
  )])))
  expect_true(is.na(unavailable["ari"]))
})

test_that("PIP RMSD is continuous and uses only mutually regularized matched cells", {
  lam1 <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  lam2 <- cbind(-lam1[, 1], rep(.01, 4), -lam1[, 2])
  pip1 <- matrix(.1, 4, 2)
  pip2 <- matrix(.9, 4, 3)       # unmatched column is irrelevant
  pip2[, 1] <- pip1[, 1]
  pip2[, 3] <- pip1[, 2]
  pip1[3, 2] <- .49
  pip2[3, 3] <- .09              # continuous change without a .5 crossing
  pip1[4, 2] <- .3
  pip2[4, 3] <- .1
  pip2[1, 1] <- 1                # excluded fixed-cell change

  regularized1 <- matrix(FALSE, 4, 2)
  regularized2 <- matrix(FALSE, 4, 3)
  regularized1[3:4, 2] <- TRUE
  regularized2[3:4, 3] <- TRUE
  regularized1[2, 1] <- TRUE     # only one endpoint TRUE: excluded

  metrics <- vbpm:::.match_sweep_loadings(
    lam1, lam2, K0 = 1, pip1 = pip1, pip2 = pip2,
    regularized1 = regularized1, regularized2 = regularized2
  )
  expect_equal(unname(metrics["pip_rmsd"]), sqrt((.4^2 + .2^2) / 2),
               tolerance = 1e-12)

  none <- vbpm:::.match_sweep_loadings(
    cbind(c(.8, .8)), cbind(c(.8, .8)), K0 = 1,
    pip1 = cbind(c(.9, .1)), pip2 = cbind(c(.1, .9)),
    regularized1 = matrix(FALSE, 2, 1),
    regularized2 = matrix(FALSE, 2, 1)
  )
  expect_true(is.na(none["pip_rmsd"]))
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

test_that("transition tables have the exact 0.8.3 schema and typed empty columns", {
  empty <- vbpm:::.transition_table()
  expect_s3_class(empty, "data.frame")
  expect_named(empty, .transition_names_083)
  expect_identical(nrow(empty), 0L)
  expect_type(empty$K_from, "integer")
  expect_type(empty$K_to, "integer")
  expect_true(all(vapply(
    empty[setdiff(names(empty), c("K_from", "K_to", "collision"))],
    is.double, logical(1)
  )))
  expect_type(empty$collision, "logical")

  one <- matrix(c(.8, .8, 0, 0), 4, 1)
  object <- .make_sweep_083(list(one), list(matrix(.5, 4, 1)))
  built <- vbpm:::.build_transitions(
    object$sweep, object$loadings, object$pips,
    Q0 = matrix(-1, 4, 1), stability_eps = .1
  )
  expect_named(built, .transition_names_083)
  expect_identical(nrow(built), 0L)

  one_groups <- cbind(one, rev(one))
  bifactor_one <- vbpm:::.build_transitions(
    object$sweep,
    stats::setNames(list(cbind(general = rep(.7, 4), one_groups)), "2"),
    stats::setNames(list(matrix(.5, 4, 3)), "2"),
    Q0 = matrix(-1, 4, 1), stability_eps = .1, bifactor = TRUE
  )
  expect_named(bifactor_one, .transition_names_083)
  expect_identical(nrow(bifactor_one), 0L)
  expect_type(bifactor_one$collision, "logical")

  empty_sweep <- object$sweep[FALSE, , drop = FALSE]
  bifactor_zero <- vbpm:::.build_transitions(
    empty_sweep, list(), list(), Q0 = matrix(-1, 4, 1),
    stability_eps = .1, bifactor = TRUE
  )
  expect_named(bifactor_zero, .transition_names_083)
  expect_identical(nrow(bifactor_zero), 0L)
})

test_that("the builder computes bifactor group diagnostics and gates failures", {
  first <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  second <- cbind(first, c(.7, 0, 0, .7))
  pips <- list(matrix(.5, 4, 2), matrix(.5, 4, 3))
  object <- .make_sweep_083(list(first, second), pips)
  Q0 <- matrix(c(1, -1, 0, -1), 4, 1)

  ordinary <- vbpm:::.build_transitions(
    object$sweep, object$loadings, object$pips, Q0 = Q0,
    stability_eps = .1
  )
  expect_named(ordinary, .transition_names_083)
  expect_identical(ordinary$K_from, 2L)
  expect_identical(ordinary$K_to, 3L)
  expect_equal(ordinary$ELBO_gain, 5)
  expect_equal(ordinary$BIC_gain, 20)
  expect_true(all(is.finite(unlist(ordinary[c(
    "rmsd", "rmsd_max", "phi_min", "pip_rmsd",
    "unmatched_max"
  )], use.names = FALSE))))
  expect_false(ordinary$collision)

  failed <- object$sweep
  failed$converged[2] <- FALSE
  expect_warning(
    failed_edge <- vbpm:::.build_transitions(
      failed, object$loadings, object$pips, Q0 = Q0, stability_eps = .1
    ),
    NA
  )
  expect_equal(failed_edge$ELBO_gain, 5)
  expect_equal(failed_edge$BIC_gain, 20)
  expect_true(all(is.na(failed_edge[.transition_names_083[-(1:4)]])))

  general1 <- c(.95, .8, .7, .6)
  general2 <- c(-.1, -.2, -.3, -.4)
  bifactor_loadings <- stats::setNames(list(
    cbind(general1, first), cbind(general2, second)
  ), names(object$loadings))
  bifactor_pips <- stats::setNames(list(
    cbind(c(.01, .02, .03, .04), pips[[1L]]),
    cbind(c(.99, .98, .97, .96), pips[[2L]])
  ), names(object$pips))
  bifactor_edge <- vbpm:::.build_transitions(
    object$sweep, bifactor_loadings, bifactor_pips, Q0 = Q0,
    stability_eps = .1, bifactor = TRUE
  )
  expect_named(bifactor_edge, .transition_names_083)
  expect_equal(bifactor_edge$ELBO_gain, 5)
  expect_equal(bifactor_edge$BIC_gain, 20)
  expect_equal(
    bifactor_edge[vbpm:::.transition_numeric_metric_names],
    ordinary[vbpm:::.transition_numeric_metric_names],
    tolerance = 1e-12
  )
  expect_identical(bifactor_edge$collision, ordinary$collision)

  failed_bifactor <- vbpm:::.build_transitions(
    failed, bifactor_loadings, bifactor_pips, Q0 = Q0,
    stability_eps = .1, bifactor = TRUE
  )
  expect_true(all(is.na(
    failed_bifactor[.transition_names_083[-(1:4)]]
  )))

  malformed_bifactor_loadings <- bifactor_loadings
  malformed_bifactor_loadings[["3"]] <- second
  expect_warning(
    malformed_bifactor <- vbpm:::.build_transitions(
      object$sweep, malformed_bifactor_loadings, bifactor_pips, Q0 = Q0,
      stability_eps = .1, bifactor = TRUE
    ),
    "malformed or nonfinite"
  )
  expect_equal(malformed_bifactor$ELBO_gain, 5)
  expect_equal(malformed_bifactor$BIC_gain, 20)
  expect_true(all(is.na(
    malformed_bifactor[.transition_names_083[-(1:4)]]
  )))

  malformed_bifactor_pips <- bifactor_pips
  malformed_bifactor_pips[["3"]] <- pips[[2L]]
  expect_warning(
    malformed_pip <- vbpm:::.build_transitions(
      object$sweep, bifactor_loadings, malformed_bifactor_pips, Q0 = Q0,
      stability_eps = .1, bifactor = TRUE
    ),
    "malformed or nonfinite"
  )
  expect_true(is.na(malformed_pip$pip_rmsd))
  expect_equal(
    malformed_pip[vbpm:::.transition_numeric_metric_names[1:4]],
    ordinary[vbpm:::.transition_numeric_metric_names[1:4]],
    tolerance = 1e-12
  )
  expect_identical(malformed_pip$collision, ordinary$collision)
})

test_that("bifactor collision is determined only by the group block", {
  x <- c(.8, .8, 0, 0)
  y <- c(0, 0, .8, .8)
  z <- c(.4, 0, .4, 0)
  group1 <- cbind(x, x)
  group2 <- cbind(x, y, z)
  pips1 <- matrix(.5, 4, 2)
  pips2 <- matrix(.5, 4, 3)
  dimnames(pips1) <- dimnames(group1)
  dimnames(pips2) <- dimnames(group2)
  object <- .make_sweep_083(list(group1, group2), list(pips1, pips2))
  Q0 <- matrix(integer(0), 4, 0)

  ordinary <- vbpm:::.build_transitions(
    object$sweep, object$loadings, object$pips, Q0 = Q0,
    stability_eps = .1
  )
  bifactor <- vbpm:::.build_transitions(
    object$sweep,
    stats::setNames(list(
      cbind(general = c(.9, .8, .7, .6), group1),
      cbind(general = c(-.4, -.3, -.2, -.1), group2)
    ), names(object$loadings)),
    stats::setNames(list(
      cbind(general = c(.01, .02, .03, .04), pips1),
      cbind(general = c(.99, .98, .97, .96), pips2)
    ), names(object$pips)),
    Q0 = Q0, stability_eps = .1, bifactor = TRUE
  )

  expect_true(ordinary$collision)
  expect_true(bifactor$collision)
  expect_equal(
    bifactor[vbpm:::.transition_numeric_metric_names],
    ordinary[vbpm:::.transition_numeric_metric_names],
    tolerance = 1e-12
  )
})

test_that("the builder uses K keys and excludes fixed Q0 cells from PIP RMSD", {
  first <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  second <- cbind(first, c(.7, 0, 0, .7))
  pip1 <- matrix(.1, 4, 2)
  pip2 <- cbind(pip1, rep(.9, 4))
  pip2[1, 1] <- .9              # fixed by Q0 and therefore excluded
  pip2[2, 1] <- .3              # regularized backbone difference .2
  pip2[3, 2] <- .5              # added-factor difference .4
  object <- .make_sweep_083(list(first, second), list(pip1, pip2))
  Q0 <- matrix(c(1, -1, 0, -1), 4, 1)

  ## Reverse storage order to ensure K names, rather than positions, are used.
  transitions <- vbpm:::.build_transitions(
    object$sweep, rev(object$loadings), rev(object$pips), Q0 = Q0,
    stability_eps = .1
  )
  expect_equal(
    transitions$pip_rmsd,
    sqrt((.2^2 + .4^2) / 6),
    tolerance = 1e-12
  )
})

test_that("malformed PIPs leave loading diagnostics intact", {
  first <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  second <- cbind(first, c(.7, 0, 0, .7))
  object <- .make_sweep_083(
    list(first, second),
    list(matrix(.5, 4, 2), matrix(.5, 4, 3))
  )
  object$pips[["3"]][1, 1] <- Inf

  expect_warning(
    transitions <- vbpm:::.build_transitions(
      object$sweep, object$loadings, object$pips,
      Q0 = matrix(-1, 4, 1), stability_eps = .1
    ),
    "malformed or nonfinite"
  )
  expect_true(is.na(transitions$pip_rmsd))
  expect_true(all(is.finite(unlist(transitions[c(
    "rmsd", "rmsd_max", "phi_min", "unmatched_max"
  )], use.names = FALSE))))
  expect_false(transitions$collision)
  expect_equal(transitions$ELBO_gain, 5)
  expect_equal(transitions$BIC_gain, 20)
})

test_that("malformed converged endpoints warn once without losing finite gains", {
  first <- cbind(c(.8, .8, 0, 0), c(0, 0, .8, .8))
  second <- cbind(first, c(.7, 0, 0, .7))
  third <- cbind(second, c(0, .7, .7, 0))
  malformed <- second
  malformed[1, 1] <- Inf
  loadings <- list(first, malformed, third)
  pips <- lapply(loadings, function(x) matrix(.5, nrow(x), ncol(x)))
  object <- .make_sweep_083(loadings, pips)

  expect_error(
    vbpm:::.match_sweep_loadings(
      first, malformed, K0 = 1, eps = .1,
      pip1 = pips[[1]], pip2 = pips[[2]],
      regularized1 = matrix(TRUE, 4, 2),
      regularized2 = matrix(TRUE, 4, 3)
    ),
    "finite numeric matrices"
  )

  warnings <- character(0)
  transitions <- withCallingHandlers(
    vbpm:::.build_transitions(
      object$sweep, object$loadings, object$pips,
      Q0 = matrix(-1, 4, 1), stability_eps = .1
    ),
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
  expect_true(all(is.na(transitions[.transition_names_083[-(1:4)]])))
})
