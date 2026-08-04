## Compact two-table PEFA contract (0.8.2).

.make_pefa_data <- function(N = 180, K = 3, seed = 4) {
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

.transition_names <- c(
  "K_from", "K_to", "ELBO_gain", "BIC_gain",
  "rmsd", "rmsd_max", "phi_min", "ari", "pip_sum", "pip_product",
  "max_unmatched_loading"
)

.sweep_names <- c(
  "K", "ELBO", "AIC", "BIC", "RMSEA", "SRMR", "CFI", "TLI",
  "t", "iter", "secs", "converged", "loading", "pip"
)

test_that("the elbow rule requires the full sustained run", {
  expect_identical(
    select_K_elbow(1:5, cumsum(c(0, 100, 5, 4, 3)), delta = 10, sustain = 2),
    2L
  )
  expect_identical(
    select_K_elbow(1:5, cumsum(c(0, 100, 50, 20, 5)), delta = 10, sustain = 2),
    5L
  )
  expect_identical(select_K_elbow(2L, 10), 2L)
  expect_identical(select_K_elbow(2:4, c(10, 9, 8)), 2L)
  expect_identical(
    select_K_elbow(1:4, cumsum(c(0, 100, 10, 5)),
                   delta = 10, sustain = 1),
    3L
  )
})

test_that("pefa has the exact 0.8.2 public signature", {
  expect_identical(
    names(formals(pefa)),
    c("Q0", "Y", "Kmin", "Kmax", "cuts", "sustain", "bifactor",
      "general", "max_it", "tau", "stability_eps", "rank_adjust",
      "rank_max_J", "verbose", "...")
  )
})

test_that("named cuts and compact gain selections are validated", {
  expect_identical(vbpm:::.pefa_cuts(c(` primary ` = 20)), c(primary = 20))
  expect_error(vbpm:::.pefa_cuts(10), "names")
  expect_error(vbpm:::.pefa_cuts(c(primary = 10, primary = 20)), "unique")
  expect_error(vbpm:::.pefa_cuts(c(other = 10)), "exactly one")
  expect_error(vbpm:::.pefa_cuts(c(primary = 101)), "\\[0, 100\\]")

  cuts <- c(primary = 20, sensitivity = 10)
  selected <- vbpm:::.sweep_gain_selections(
    2:6, c(100, 5, 4, 3), cuts, sustain = 2, complete = TRUE
  )
  expect_identical(names(selected$selected_K), names(cuts))
  expect_identical(names(selected$boundary), names(cuts))
  expect_identical(selected$selected_K, c(primary = 3L, sensitivity = 3L))
  expect_true(all(selected$boundary == "interior"))

  unavailable <- vbpm:::.sweep_gain_selections(
    2:6, c(NA_real_, 5, 4, 3), cuts, sustain = 2, complete = FALSE
  )
  expect_true(all(is.na(unavailable$selected_K)))
  expect_true(all(unavailable$boundary == "none"))

  single <- vbpm:::.sweep_gain_selections(
    2L, numeric(0), c(primary = 10), sustain = 2, complete = TRUE
  )
  expect_identical(single$selected_K, c(primary = 2L))
  expect_identical(single$boundary, c(primary = "single"))
})

test_that("the internal runner is model-neutral and preserves typed failures", {
  fit_one <- function(K) {
    if (K == 3L) stop("synthetic failure")
    list(K = K)
  }
  payload <- function(K, elapsed, converged) {
    list(
      row = data.frame(
        K = K, ELBO = if (converged) K else NA_real_,
        BIC = if (converged) -K else NA_real_, converged = converged
      ),
      loading = matrix(if (converged) K else NA_real_, 2, K),
      pip = matrix(if (converged) .5 else NA_real_, 2, K)
    )
  }
  sweep <- vbpm:::.run_sweep(
    2:3,
    fit_one = fit_one,
    extract_one = function(fit, K) {
      Sys.sleep(.05)
      payload(K, NA_real_, TRUE)
    },
    failure_one = function(K, elapsed, error) payload(K, elapsed, FALSE)
  )
  expect_identical(sweep$K, 2:3)
  expect_identical(sweep$converged, c(TRUE, FALSE))
  expect_true(is.matrix(sweep$loading[[1L]]))
  expect_true(all(is.na(sweep$pip[[2L]])))
  expect_gte(sweep$secs[[1L]], .02)
  extraction_failed <- vbpm:::.run_sweep(
    2L,
    fit_one = function(K) list(K = K),
    extract_one = function(fit, K) stop("synthetic extraction failure"),
    failure_one = function(K, elapsed, error) payload(K, elapsed, FALSE)
  )
  expect_false(extraction_failed$converged[[1L]])
  expect_true(is.finite(extraction_failed$secs[[1L]]))
  expect_error(
    vbpm:::.run_sweep(
      c(2L, 4L), fit_one = fit_one,
      extract_one = function(fit, K) payload(K, NA_real_, TRUE),
      failure_one = function(K, elapsed, error) payload(K, elapsed, FALSE)
    ),
    "consecutive increasing"
  )
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
  stats <- fit_stats(fit, Q = Q)
  expect_identical(unname(stats[["CFI"]]), 1)
})

test_that("pefa returns the exact compact two-table object", {
  dat <- .make_pefa_data()
  object <- pefa(
    dat$Q0, dat$Y, Kmin = 2, Kmax = 3,
    cuts = c(primary = 20, sensitivity = 10),
    v0 = .001, max_it = 300, verbose = FALSE
  )

  expect_s3_class(object, "pefa")
  expect_false(inherits(object, "vbpm_sweep"))
  expect_named(object, c(
    "call", "model", "objective_type", "Q0", "selected_K", "boundary",
    "sweep", "transitions", "window", "cuts", "sustain", "bifactor",
    "general", "tau", "stability_eps", "rank_adjust", "rank_max_J"
  ))
  expect_identical(object$model, "vbfa")
  expect_identical(object$objective_type, "elbo")
  expect_identical(object$Q0, dat$Q0)
  expect_named(object$sweep, .sweep_names)
  expect_named(object$transitions, .transition_names)
  expect_identical(object$window, c(Kmin = 2L, Kmax = 3L, K0 = 2L))
  expect_identical(names(object$selected_K), names(object$cuts))
  expect_identical(names(object$boundary), names(object$cuts))
  expect_null(object$fits)
  expect_null(object$selection)
  expect_null(object$selected_fit)
  expect_false(exists("selected_fit", envir = asNamespace("vbpm"),
                      inherits = FALSE))

  for (i in seq_len(nrow(object$sweep))) {
    expect_true(is.matrix(object$sweep$loading[[i]]))
    expect_true(is.matrix(object$sweep$pip[[i]]))
    expect_identical(dim(object$sweep$loading[[i]]),
                     c(ncol(dat$Y), object$sweep$K[i]))
    expect_identical(dim(object$sweep$pip[[i]]),
                     dim(object$sweep$loading[[i]]))
  }
  expect_false(any(grepl("_gain", names(object$sweep))))
  expect_identical(nrow(object$transitions), nrow(object$sweep) - 1L)

  summary <- summary(object)
  expect_s3_class(summary, "summary.pefa")
  expect_false(any(c("loading", "pip") %in% names(summary$sweep)))
  expect_identical(summary$selected_K, object$selected_K)
  expect_identical(summary$boundary, object$boundary)
  expect_output(print(object), "ELBO gain \\[primary = 20%\\]")
  expect_output(print(summary), "Candidates")
})

test_that("pefa validates new controls before fitting", {
  dat <- .make_pefa_data(N = 80)
  expect_error(
    pefa(dat$Q0, dat$Y, 2, 2, stability_eps = 0, verbose = FALSE),
    "stability_eps"
  )
  expect_error(
    pefa(dat$Q0, dat$Y, 2, 2, rank_adjust = 1, verbose = FALSE),
    "rank_adjust"
  )
  expect_error(
    pefa(dat$Q0, dat$Y, 2, 2, rank_max_J = 1.5, verbose = FALSE),
    "rank_max_J"
  )
  expect_error(
    pefa(dat$Q0, dat$Y, 2, 2, rank_adjust = TRUE, rank_max_J = 2,
         verbose = FALSE),
    "J <= rank_max_J"
  )
  expect_error(
    pefa(dat$Q0, dat$Y, 2, 2, bifactor = TRUE, general = c(1, 0),
         verbose = FALSE),
    "scalar or length-J"
  )
})

test_that("print gives actionable advice for every named open boundary", {
  object <- structure(
    list(
      bifactor = FALSE,
      window = c(Kmin = 2L, Kmax = 4L, K0 = 2L),
      cuts = c(primary = 10),
      selected_K = c(primary = 4L),
      boundary = c(primary = "upper"),
      sweep = data.frame(K = 2:4, converged = TRUE)
    ),
    class = "pefa"
  )
  expect_output(print(object), "larger Kmax")

  object$selected_K[["primary"]] <- 2L
  object$boundary[["primary"]] <- "lower"
  expect_output(print(object), "reduce or revise Q0")

  object$window <- c(Kmin = 3L, Kmax = 4L, K0 = 2L)
  object$selected_K[["primary"]] <- 3L
  expect_output(print(object), "smaller Kmin")

  object$window <- c(Kmin = 2L, Kmax = 4L, K0 = 2L)
  object$cuts <- c(primary = 20, sensitivity = 10)
  object$selected_K <- c(primary = 3L, sensitivity = 4L)
  object$boundary <- c(primary = "interior", sensitivity = "upper")
  expect_output(print(object), "sensitivity.*larger Kmax")
})

test_that("pefa stores fit statistics under the requested rank policy", {
  set.seed(9)
  Y <- matrix(rnorm(80 * 4), 80, 4)
  ## One indicator on one factor makes its loading and uniqueness redundant
  ## in the implied covariance, so the explicit Jacobian policy differs.
  Q0 <- matrix(c(1L, 0L, 0L, 0L), 4, 1)
  fit <- vbfa(Y, Q0, v0 = .001, max_it = 1500)

  nominal <- pefa(
    Q0, Y, 1, 1, v0 = .001, max_it = 1500,
    rank_adjust = FALSE, rank_max_J = 3e9, verbose = FALSE
  )
  adjusted <- pefa(
    Q0, Y, 1, 1, v0 = .001, max_it = 1500,
    rank_adjust = TRUE, verbose = FALSE
  )
  direct_nominal <- fit_stats(
    fit, Y = Y, Q = fit$Q, rank_adjust = FALSE
  )
  direct_adjusted <- fit_stats(
    fit, Y = Y, Q = fit$Q, rank_adjust = TRUE
  )
  fields <- c("AIC", "BIC", "RMSEA", "SRMR", "CFI", "TLI", "t")

  expect_equal(
    unname(unlist(nominal$sweep[1L, fields], use.names = TRUE)),
    unname(direct_nominal[fields])
  )
  expect_equal(
    unname(unlist(adjusted$sweep[1L, fields], use.names = TRUE)),
    unname(direct_adjusted[fields])
  )
  expect_false(nominal$rank_adjust)
  expect_true(adjusted$rank_adjust)
  expect_identical(nominal$rank_max_J, 3e9)
  expect_lt(adjusted$sweep$t, nominal$sweep$t)
})

test_that("a failed sweep emits one aggregate selection warning", {
  dat <- .make_pefa_data(N = 100)
  warnings <- character(0)
  object <- withCallingHandlers(
    pefa(dat$Q0, dat$Y, Kmin = 2, Kmax = 3,
         v0 = .001, max_it = 1, verbose = FALSE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_match(warnings, "selection is unavailable")
  expect_true(all(is.na(object$selected_K)))
  expect_true(all(object$boundary == "none"))
})

test_that("bifactor sweeps store K plus general matrices and gate stability", {
  dat <- .make_pefa_data(N = 120)
  object <- suppressWarnings(pefa(
    dat$Q0, dat$Y, 2, 3, bifactor = TRUE,
    cuts = c(primary = 20, sensitivity = 10),
    v0 = .001, max_it = 300, verbose = FALSE
  ))
  expect_true(object$bifactor)
  expect_identical(object$general, rep(1L, ncol(dat$Y)))
  expect_identical(object$sweep$K, 2:3)
  expect_identical(
    vapply(object$sweep$loading, ncol, integer(1)), 3:4
  )
  expect_true(all(is.na(object$transitions[vbpm:::.transition_metric_names])))
  expect_output(print(object), "group factors [(][+] 1 general")
  expect_output(print(summary(object)), "K_total")
})
