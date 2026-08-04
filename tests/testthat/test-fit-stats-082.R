## Focused regression tests for the 0.8.2 parameter-count policy.

make_rank_fit_082 <- function(J = 3L) {
  Q <- matrix(0L, J, 1L)
  Q[1L, 1L] <- 1L
  structure(
    list(
      Q = Q,
      Lam = matrix(c(0.7, rep(0, J - 1L)), J, 1L),
      pi = matrix(c(1, rep(0, J - 1L)), J, 1L),
      orthogonal = FALSE,
      nobs = 120L,
      ld = FALSE,
      PsiInv = rep(2, J),
      Phi = matrix(1, 1L, 1L),
      sample_cov = diag(J),
      ELBO = -20,
      objective = -20,
      objective_type = "elbo"
    ),
    class = c("vbfa", "vbpm_fit")
  )
}

test_that("nominal counting is the deterministic default", {
  fit <- make_rank_fit_082()

  ## Shadow requireNamespace() only for a copy of the S3 method. If the
  ## nominal branch inspects the optional package, this sentinel fails.
  method <- getS3method("fit_stats", "vbfa")
  method_env <- new.env(parent = environment(method))
  method_env$requireNamespace <- function(...) {
    stop("nominal path touched the optional namespace", call. = FALSE)
  }
  environment(method) <- method_env

  stats <- method(fit)
  expect_identical(unname(stats[["t"]]), unname(stats[["t_nom"]]))
})

test_that("rank controls reject malformed values", {
  fit <- make_rank_fit_082()

  bad_adjust <- list(NULL, NA, 1, "TRUE", logical(), c(TRUE, FALSE))
  for (value in bad_adjust) {
    expect_error(
      fit_stats(fit, rank_adjust = value),
      "rank_adjust.*one non-missing logical"
    )
  }

  bad_max <- list(NULL, NA_real_, NaN, Inf, -Inf, 0, -1, 1.5, TRUE,
                  "100", numeric(), c(50, 100))
  for (value in bad_max) {
    expect_error(
      fit_stats(fit, rank_max_J = value),
      "rank_max_J.*finite positive whole number"
    )
  }

  expect_silent(fit_stats(fit, rank_adjust = FALSE, rank_max_J = 1L))
})

test_that("an unavailable numDeriv package is an explicit opt-in error", {
  fit <- make_rank_fit_082()

  ## Give a private copy of the method a deterministic unavailable-package
  ## result without changing the package namespace or the machine libraries.
  method <- getS3method("fit_stats", "vbfa")
  method_env <- new.env(parent = environment(method))
  method_env$requireNamespace <- function(package, quietly = FALSE) FALSE
  environment(method) <- method_env

  expect_error(
    method(fit, rank_adjust = TRUE),
    "rank_adjust = TRUE.*requires.*numDeriv"
  )
})

test_that("large explicit rank requests do not silently fall back", {
  fit <- make_rank_fit_082(J = 4L)

  expect_error(
    fit_stats(fit, rank_adjust = TRUE, rank_max_J = 3L),
    "J <= rank_max_J.*Raise.*rank_max_J.*rank_adjust = FALSE"
  )
})

test_that("explicit rank adjustment computes the Jacobian rank", {
  skip_if_not_installed("numDeriv")
  fit <- make_rank_fit_082()

  nominal <- fit_stats(fit)
  adjusted <- fit_stats(fit, rank_adjust = TRUE)

  ## With one indicator, its loading and uniqueness contribute to the same
  ## variance moment, so this deliberately underidentified example loses one
  ## degree of Jacobian rank.
  expect_identical(unname(nominal[["t"]]), 4)
  expect_identical(unname(adjusted[["t_nom"]]), 4)
  expect_identical(unname(adjusted[["t"]]), 3)
})
