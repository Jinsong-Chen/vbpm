## Compact PEFA edge table and adjacent-count summaries (0.8.1).

.make_pefa_data <- function(N = 300, K = 3, seed = 4) {
  set.seed(seed)
  J <- K * 6
  group <- rep(seq_len(K), each = 6)
  loading <- matrix(0, J, K)
  loading[cbind(seq_len(J), group)] <- .7
  score <- matrix(rnorm(N * K), N, K)
  Y <- score %*% t(loading) +
    matrix(rnorm(N * J, sd = sqrt(1 - .49)), N, J)
  Q0 <- matrix(-1L, J, 2)
  for (k in 1:2) {
    anchor <- which(group == k)[1:2]
    Q0[anchor, ] <- 0L
    Q0[anchor, k] <- 1L
  }
  list(Y = Y, Q0 = Q0)
}

.simple_loadings <- function() {
  A <- rbind(
    c(.8, .1), c(.7, .1),
    c(.1, .8), c(.1, .7)
  )
  surplus <- c(.05, .06, .04, .03)
  B <- cbind(-A[, 2], A[, 1], surplus)
  C <- cbind(B, c(.02, .03, .01, .02))
  list(A = A, B = B, C = C, surplus = surplus)
}

.transition_names <- c(
  "K_from", "K_to", "ELBO_gain", "BIC_gain",
  "phi_min", "phi_mean", "ari", "rmsd", "max_loading"
)

.selection_names <- c(
  "criterion", "form", "cut_name", "cut", "gain_max", "threshold",
  "selected_K", "boundary"
)

test_that("the elbow rule requires the full sustained run", {
  expect_identical(
    select_K_elbow(1:5, cumsum(c(0, 100, 5, 4, 3)), delta = 10, sustain = 2),
    2L
  )
  ## A lone final drop is not a shortened right-edge run.
  expect_identical(
    select_K_elbow(1:5, cumsum(c(0, 100, 50, 20, 5)), delta = 10, sustain = 2),
    5L
  )
  expect_identical(select_K_elbow(2L, 10), 2L)
  expect_identical(select_K_elbow(2:4, c(10, 9, 8)), 2L)
  ## Equality does not count: the paper's comparison is strictly below.
  expect_identical(
    select_K_elbow(1:4, cumsum(c(0, 100, 10, 5)),
                   delta = 10, sustain = 1),
    3L
  )
})

test_that("pefa has the exact compact public signature", {
  expect_identical(
    names(formals(pefa)),
    c("Q0", "Y", "Kmin", "Kmax", "cuts", "sustain", "bifactor",
      "general", "max_it", "tau", "verbose", "...")
  )
})

test_that("selection contains only ELBO and BIC raw/gain rules", {
  sweep <- data.frame(
    K = 2:6,
    ELBO = c(0, 100, 105, 109, 112),
    BIC = c(200, 100, 95, 91, 88),
    converged = rep(TRUE, 5)
  )
  transitions <- vbpm:::.transition_table(4L)
  transitions$K_from <- 2:5
  transitions$K_to <- 3:6
  transitions$ELBO_gain <- c(100, 5, 4, 3)
  transitions$BIC_gain <- c(100, 5, 4, 3)

  selection <- vbpm:::.pefa_selection(
    sweep, transitions, c(primary = 10, sensitivity = 20), sustain = 2
  )
  expect_named(selection, .selection_names)
  expect_identical(nrow(selection), 6L)
  expect_setequal(unique(selection$criterion), c("elbo", "bic"))
  expect_true(all(is.na(selection$cut_name[selection$form == "raw"])))
  expect_true(all(is.na(selection$cut[selection$form == "raw"])))

  gain <- selection[selection$form == "gain", ]
  expect_equal(gain$gain_max, rep(100, 4))
  expect_equal(gain$threshold, c(10, 20, 10, 20))
  expect_equal(gain$selected_K, rep(3L, 4))
  expect_true(all(gain$boundary == "interior"))

  raw <- selection[selection$form == "raw", ]
  expect_equal(raw$selected_K, c(6L, 6L))
  expect_true(all(raw$boundary == "upper"))
})

test_that("a broken path makes that criterion unavailable", {
  sweep <- data.frame(
    K = 2:5,
    ELBO = c(0, 100, 105, 109),
    BIC = c(200, 100, 95, 91),
    converged = c(TRUE, FALSE, TRUE, TRUE)
  )
  transitions <- vbpm:::.transition_table(3L)
  transitions$K_from <- 2:4
  transitions$K_to <- 3:5
  transitions$ELBO_gain <- c(NA, NA, 4)
  transitions$BIC_gain <- c(NA, NA, 4)

  expect_warning(
    selection <- vbpm:::.pefa_selection(
      sweep, transitions, c(primary = 10), sustain = 2
    ),
    "ELBO and BIC"
  )
  expect_true(all(is.na(selection$selected_K)))
  expect_true(all(selection$boundary == "none"))

  ## Non-finite ELBO does not disable an otherwise complete BIC path.
  sweep$converged[] <- TRUE
  sweep$ELBO[2] <- Inf
  transitions$BIC_gain <- c(100, 5, 4)
  expect_warning(
    selection <- vbpm:::.pefa_selection(
      sweep, transitions, c(primary = 10), sustain = 2
    ),
    "ELBO"
  )
  expect_true(all(is.na(selection$selected_K[selection$criterion == "elbo"])))
  expect_true(all(!is.na(selection$selected_K[selection$criterion == "bic"])))
})

test_that("named cuts are validated and normalized", {
  expect_identical(vbpm:::.pefa_cuts(c(` primary ` = 10)), c(primary = 10))
  expect_error(vbpm:::.pefa_cuts(10), "names")
  expect_error(vbpm:::.pefa_cuts(c(primary = 10, primary = 20)), "unique")
  expect_error(vbpm:::.pefa_cuts(c(other = 10)), "exactly one")
  expect_error(vbpm:::.pefa_cuts(c(primary = 101)), "\\[0, 100\\]")
  bad_names <- c(10, 20)
  names(bad_names) <- c("primary", NA_character_)
  expect_error(vbpm:::.pefa_cuts(bad_names), "unique, non-empty")
})

test_that("hard-selected CFI is one when its denominator is zero", {
  J <- 3L
  Q <- matrix(0L, J, 1L)
  fit <- structure(
    list(
      Lam = matrix(0, J, 1L), pi = matrix(0, J, 1L),
      Phi = matrix(1), PsiInv = rep(1, J), ld = FALSE,
      orthogonal = FALSE, sample_cov = diag(J), nobs = 100L,
      ELBO = 0, objective = 0, objective_type = "elbo", Q = Q
    ),
    class = c("vbfa", "vbpm_fit")
  )
  stats <- fit_stats(fit, Q = Q, rank_adjust = FALSE)
  expect_identical(unname(stats[["CFI"]]), 1)
})

test_that("best-unused alignment is sign invariant", {
  loading <- .simple_loadings()
  matched <- vbpm:::.match_columns(loading$A, loading$B)

  expect_identical(matched$idx, c(2L, 1L))
  expect_equal(matched$phi, c(1, 1))
  expect_identical(matched$sign, c(1, -1))
  expect_identical(matched$unmatched, 3L)

  flipped <- loading$A
  flipped[, 1] <- -flipped[, 1]
  expect_equal(
    vbpm:::.transition_stability(flipped, loading$B),
    vbpm:::.transition_stability(loading$A, loading$B)
  )
})

test_that("adjacent stability returns the five retained summaries", {
  loading <- .simple_loadings()
  stability <- vbpm:::.transition_stability(loading$A, loading$B)

  expect_named(
    stability,
    c("phi_min", "phi_mean", "ari", "rmsd", "max_loading")
  )
  expect_equal(unname(stability[c("phi_min", "phi_mean")]), c(1, 1))
  expect_equal(unname(stability["ari"]), 1)
  expect_equal(unname(stability["rmsd"]), 0)
  expect_equal(unname(stability["max_loading"]), max(abs(loading$surplus)))
})

test_that("zero-norm endpoints and degenerate ARI return typed NA", {
  loading <- .simple_loadings()
  B0 <- loading$B
  B0[, 3] <- 0

  expect_true(all(is.na(vbpm:::.transition_stability(loading$A, B0))))
  expect_true(is.na(vbpm:::.ari(rep(1, 5), rep(1, 5))))
  expect_false(is.nan(vbpm:::.ari(rep(1, 5), rep(1, 5))))
})

test_that("transition construction has the exact compact schema", {
  loading <- .simple_loadings()
  sweep <- data.frame(
    K = 2:4,
    ELBO = c(10, 15, 16),
    BIC = c(100, 80, 78),
    converged = rep(TRUE, 3)
  )
  fits <- list(
    `2` = list(Lam = loading$A),
    `3` = list(Lam = loading$B),
    `4` = list(Lam = loading$C)
  )

  transitions <- vbpm:::.build_transitions(sweep, fits)
  expect_named(transitions, .transition_names)
  expect_identical(transitions$K_from, 2:3)
  expect_identical(transitions$K_to, 3:4)
  expect_equal(transitions$ELBO_gain, c(5, 1))
  expect_equal(transitions$BIC_gain, c(20, 2))
  expect_true(all(is.finite(transitions$phi_min)))
  expect_null(attr(transitions, "detail"))
})

test_that("one candidate gives a typed zero-row transition table", {
  sweep <- data.frame(K = 2L, ELBO = 10, BIC = 100, converged = TRUE)
  transitions <- vbpm:::.build_transitions(sweep, list(`2` = list()))

  expect_s3_class(transitions, "data.frame")
  expect_named(transitions, .transition_names)
  expect_identical(nrow(transitions), 0L)
  expect_type(transitions$K_from, "integer")
  expect_true(all(vapply(transitions[-(1:2)], is.double, logical(1))))

  selection <- vbpm:::.pefa_selection(
    sweep, transitions, c(primary = 10), sustain = 2
  )
  expect_equal(selection$selected_K, rep(2L, 4))
  expect_true(all(selection$boundary == "single"))
})

test_that("selected_fit returns NULL when primary selection is unavailable", {
  object <- structure(
    list(selected_K = NA_integer_, fits = list(`2` = list())),
    class = "pefa"
  )
  expect_null(selected_fit(object))
})

test_that("failed and non-finite endpoints produce NA gains", {
  loading <- .simple_loadings()
  sweep <- data.frame(
    K = 2:4,
    ELBO = c(10, 15, 16),
    BIC = c(100, 80, 78),
    converged = rep(TRUE, 3)
  )
  fits <- list(
    `2` = list(Lam = loading$A),
    `3` = list(Lam = loading$B),
    `4` = list(Lam = loading$C)
  )

  failed <- sweep
  failed$converged[2] <- FALSE
  transitions <- vbpm:::.build_transitions(failed, fits)
  expect_true(all(is.na(transitions$ELBO_gain)))
  expect_true(all(is.na(transitions$BIC_gain)))
  expect_true(all(is.na(transitions[c(
    "phi_min", "phi_mean", "ari", "rmsd", "max_loading"
  )])))

  nonfinite <- sweep
  nonfinite$ELBO[2] <- Inf
  transitions <- vbpm:::.build_transitions(nonfinite, fits)
  expect_true(all(is.na(transitions$ELBO_gain)))
  expect_true(all(is.finite(transitions$BIC_gain)))
})

test_that("bifactor edges retain gains but not stability summaries", {
  loading <- .simple_loadings()
  sweep <- data.frame(
    K = 2:3,
    ELBO = c(10, 15),
    BIC = c(100, 80),
    converged = c(TRUE, TRUE)
  )
  fits <- list(`2` = list(Lam = loading$A), `3` = list(Lam = loading$B))
  transitions <- vbpm:::.build_transitions(sweep, fits, bifactor = TRUE)

  expect_equal(transitions$ELBO_gain, 5)
  expect_equal(transitions$BIC_gain, 20)
  expect_true(all(is.na(transitions[c(
    "phi_min", "phi_mean", "ari", "rmsd", "max_loading"
  )])))
})

test_that("pefa exposes only the retained transition columns", {
  data <- .make_pefa_data()
  object <- pefa(
    data$Q0, data$Y, Kmin = 2, Kmax = 3,
    v0 = .001, max_it = 300, verbose = FALSE
  )

  expect_named(object$transitions, .transition_names)
  expect_identical(nrow(object$transitions), nrow(object$sweep) - 1L)
  expect_false(any(grepl("_gain", names(object$sweep))))
  expect_null(attr(object$transitions, "detail"))

  nominal_t <- vapply(
    object$fits,
    function(fit) unname(fit_stats(
      fit, tau = object$tau, rank_adjust = FALSE
    )[["t"]]),
    numeric(1)
  )
  expect_equal(object$sweep$t, unname(nominal_t))
})

test_that("a failed sweep emits one aggregate selection warning", {
  data <- .make_pefa_data(N = 120)
  warnings <- character(0)
  object <- withCallingHandlers(
    pefa(data$Q0, data$Y, Kmin = 2, Kmax = 3,
         v0 = .001, max_it = 1, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_match(warnings, "selection is unavailable")
  expect_true(all(is.na(object$selection$selected_K)))
})
