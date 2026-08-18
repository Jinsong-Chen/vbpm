## persistence(), its print method, and ssl().  The accessors read stored
## fields only, so most fixtures here are hand-built evidence objects: if a
## number appears in the triangle that is not in $persistence, the pivot
## invented it.

.acc_persistence <- function() {
  data.frame(
    K_from = c(2L, 2L, 3L), K_to = c(3L, 4L, 4L), r = c(1L, 2L, 1L),
    n_pairs = c(2L, 1L, 3L),
    phi_min = c(0.912, NA_real_, 0.4567),
    rmsd = c(0.10, 0.20, 0.30), rmsd_max = c(0.20, 0.30, 0.40),
    ari = c(1, 0.5, 0.25), pip_rmsd = c(0.01, 0.02, 0.03),
    unmatched_n = c(1L, 2L, 1L), unmatched_ssl = c(0.5, 0.6, 0.7),
    collision = c(TRUE, FALSE, FALSE), n_collisions = c(2L, 0L, 0L),
    collision_targets = c("2,3", "", ""),
    collision_multiplicities = c("2,2", "", ""),
    source_ineligible_n = c(0L, 1L, 0L), target_eligible_n = c(2L, 2L, 1L),
    backbone_degenerate = c(FALSE, FALSE, FALSE),
    from_fit_ok = c(TRUE, TRUE, FALSE), to_fit_ok = c(TRUE, TRUE, TRUE),
    pair_status = rep("available", 3L), pair_reason = rep(NA_character_, 3L),
    stringsAsFactors = FALSE
  )
}

.acc_object <- function(pairs = .acc_persistence()) {
  structure(
    list(
      sweep = data.frame(K = 2:4),
      persistence = pairs,
      loadings = list(`2` = cbind(c(1, 0, 0), c(0, 2, 0)),
                      `3` = cbind(c(1, 0, 0), c(0, 2, 0), c(0, 0, 3)),
                      `4` = cbind(c(1, 0, 0), c(0, 2, 0), c(0, 0, 3),
                                  c(0, 0, 0))),
      Q0 = matrix(-1L, 3, 1),
      settings = list(bifactor = FALSE)
    ),
    class = "pefa"
  )
}

## The same evidence with nothing to mark: no collision, no reduced pair set,
## no nonconverged endpoint, and no unavailable cell.
.acc_clean <- function() {
  pairs <- .acc_persistence()
  pairs$collision <- FALSE
  pairs$n_collisions <- 0L
  pairs$n_pairs <- c(2L, 2L, 3L)
  pairs$from_fit_ok <- TRUE
  pairs$phi_min <- c(0.9, 0.8, 0.7)
  .acc_object(pairs)
}

## One small real sweep, used only where a fitted object is required.
.acc_fit <- function(bifactor = FALSE) {
  set.seed(5)
  J <- 6L
  N <- 50L
  group <- rep(1:2, each = 3)
  loading <- matrix(0, J, 2)
  loading[cbind(seq_len(J), group)] <- .8
  Y <- matrix(rnorm(N * 2), N, 2) %*% t(loading) +
    matrix(rnorm(N * J, sd = sqrt(1 - .64)), N, J)
  Q0 <- matrix(-1L, J, 1)
  Q0[1:2, 1] <- 1L
  pefa(Q0, Y, 1, 2, bifactor = bifactor, v0 = .001, max_it = 600,
       verbose = FALSE)
}

test_that("persistence() pivots the stored table into a W x W triangle", {
  x <- .acc_object()
  m <- persistence(x, "phi_min")

  expect_identical(class(m), c("vbpm_persistence", "matrix", "array"))
  expect_identical(dim(m), c(3L, 3L))
  expect_identical(dimnames(m), list(c("2", "3", "4"), c("2", "3", "4")))
  expect_identical(storage.mode(m), "double")
  expect_equal(m[["2", "3"]], 0.912, tolerance = 1e-12)
  expect_identical(m[["2", "4"]], NA_real_)
  expect_equal(m[["3", "4"]], 0.4567, tolerance = 1e-12)
  ## The diagonal and the lower half are unusable, not zero.
  expect_true(all(is.na(m[!upper.tri(m)])))

  ## The count metric keeps integer storage so a count is never printed as a
  ## fraction of one.
  counts <- persistence(x, "n_collisions")
  expect_identical(storage.mode(counts), "integer")
  expect_identical(counts[["2", "3"]], 2L)
  expect_identical(counts[["2", "4"]], 0L)
  expect_true(all(is.na(counts[!upper.tri(counts)])))

  ## Every documented metric pivots the stored column of the same name, in
  ## the (2,3), (2,4), (3,4) order the triangle reads out.
  for (metric in c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                   "unmatched_ssl", "n_collisions")) {
    pivoted <- persistence(x, metric)
    expect_equal(pivoted[upper.tri(pivoted)],
                 x$persistence[[metric]],
                 tolerance = 1e-12, info = metric)
  }
})

test_that("the three triangle attributes carry the pair facts", {
  x <- .acc_object()
  x$persistence$collision[2L] <- NA
  x$persistence$from_fit_ok[2L] <- NA
  x$persistence$n_pairs[2L] <- NA_integer_
  m <- persistence(x, "rmsd")

  for (name in c("collision", "reduced", "fit_ok")) {
    a <- attr(m, name, exact = TRUE)
    expect_true(is.logical(a), info = name)
    expect_identical(dim(a), dim(m), info = name)
    expect_identical(dimnames(a), dimnames(m), info = name)
    expect_true(all(is.na(a[!upper.tri(a)])), info = name)
    ## NA exactly where the underlying factual field is unavailable.
    expect_true(is.na(a[["2", "4"]]), info = name)
  }
  expect_identical(attr(m, "collision")[["2", "3"]], TRUE)
  expect_identical(attr(m, "reduced")[["2", "3"]], FALSE)   # n_pairs 2, k 2
  expect_identical(attr(m, "reduced")[["3", "4"]], FALSE)   # n_pairs 3, k 3
  expect_identical(attr(m, "fit_ok")[["3", "4"]], FALSE)
  expect_identical(attr(m, "fit_ok")[["2", "3"]], TRUE)

  ## A reduced pair set is n_pairs < K_from.
  y <- .acc_object()
  y$persistence$n_pairs <- c(1L, 1L, 3L)
  expect_identical(attr(persistence(y, "rmsd"), "reduced")[["2", "3"]], TRUE)
})

test_that("the metric enumeration is closed and names the right table", {
  x <- .acc_object()
  expect_identical(
    eval(formals(persistence)$metric),
    c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd", "unmatched_ssl",
      "n_collisions")
  )
  for (count in c("unmatched_n", "source_ineligible_n", "target_eligible_n",
                  "n_pairs")) {
    expect_error(persistence(x, count), "x\\$persistence", info = count)
    expect_error(persistence(x, count), "matching-process count", info = count)
  }
  expect_error(persistence(x, "backbone_ssl_min"), "x\\$sweep")
  expect_error(persistence(x, "backbone_ssl_min"), "per-candidate quantity")
  expect_error(persistence(x, "collision"), "should be one of")
  expect_error(persistence(list(a = 1), "rmsd"), "must be a \"pefa\" object")
  ## The default metric is the first documented one.
  expect_identical(persistence(x), persistence(x, "phi_min"))
})

test_that("print.vbpm_persistence marks the facts and legends only those", {
  x <- .acc_object()
  shown <- capture.output(print(persistence(x, "phi_min")))

  expect_match(paste(shown, collapse = "\n"), "0\\.91\\*")
  expect_match(paste(shown, collapse = "\n"), "\\?~")
  expect_match(paste(shown, collapse = "\n"), "0\\.46!")
  expect_identical(
    shown[length(shown)],
    paste("* = collision, ~ = reduced pair set,",
          "! = nonconverged endpoint, ? = unavailable")
  )
  ## The lower half and the diagonal print blank, not "NA".
  expect_false(any(grepl("NA", shown, fixed = TRUE)))
  expect_match(shown[length(shown) - 1L], "^4\\s*$")

  ## digits is honored and the markers travel with the reformatted value.
  wider <- capture.output(print(persistence(x, "phi_min"), digits = 3))
  expect_match(paste(wider, collapse = "\n"), "0\\.912\\*")
  expect_error(print(persistence(x, "phi_min"), digits = -1), "digits")

  ## A count prints as a whole number.
  counts <- capture.output(print(persistence(x, "n_collisions")))
  expect_match(paste(counts, collapse = "\n"), "2\\*")
  expect_false(any(grepl("2.00", counts, fixed = TRUE)))

  capture.output(expect_invisible(print(persistence(x, "rmsd"))))
})

test_that("a clean triangle prints no legend and no markers", {
  clean <- capture.output(print(persistence(.acc_clean(), "phi_min")))
  expect_match(paste(clean, collapse = "\n"), "0\\.90")
  expect_false(any(grepl("[*~!?]", clean)))
  expect_length(clean, 4L)          # header plus three rows, no legend line

  ## Restoring one fact restores exactly one legend entry.
  partial <- .acc_clean()
  partial$persistence$collision[1L] <- TRUE
  legend <- capture.output(print(persistence(partial, "phi_min")))
  expect_identical(legend[length(legend)], "* = collision")
})

test_that("the deferred accessors are not public surface in 0.9.0", {
  exports <- getNamespaceExports("vbpm")
  expect_false("pair_detail" %in% exports)
  expect_false("extend_pefa" %in% exports)
  expect_true(all(c("pefa", "persistence", "ssl", "select_K_elbow") %in%
                    exports))
  ## The per-column assignment record still exists as an internal helper, so
  ## deferring the accessor costs no evidence.
  facts <- vbpm:::.pefa_pair_facts(
    cbind(c(1, 0, 0), c(0, 1, 0)), cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)),
    matrix(0.5, 3, 2), matrix(0.5, 3, 3),
    matrix(TRUE, 3, 2), matrix(TRUE, 3, 3),
    K0 = 1, stability_eps = 0.1,
    loading_ok_from = TRUE, loading_ok_to = TRUE,
    pip_ok_from = TRUE, pip_ok_to = TRUE
  )
  expect_named(facts$.assign, c("source", "target", "sign", "congruence",
                                "rmsd", "backbone", "eligible_source"))
  expect_identical(nrow(facts$.assign), 2L)
})

test_that("ssl() is the one definition of the column sizes", {
  x <- .acc_object()
  expect_identical(ssl(x), list(`2` = c(1, 4), `3` = c(1, 4, 9),
                                `4` = c(1, 4, 9, 0)))
  expect_identical(names(ssl(x)), names(x$loadings))
  expect_null(names(ssl(x)[["2"]]))
  expect_error(ssl(list()), "must be a \"pefa\" object")

  ## summary() calls it rather than recomputing, and the object stores no copy.
  toy <- .acc_fit()
  expect_identical(summary(toy)$ssl, ssl(toy))
  expect_false("ssl" %in% names(toy))
  expect_identical(unname(lengths(ssl(toy))), toy$sweep$K)
  expect_equal(ssl(toy)[["2"]], unname(colSums(toy$loadings[["2"]]^2)),
               tolerance = 1e-12)

  ## Bifactor keeps all K + 1 columns, general factor first.
  bifactor <- .acc_fit(bifactor = TRUE)
  expect_identical(unname(lengths(ssl(bifactor))), bifactor$sweep$K + 1L)
  expect_identical(summary(bifactor)$ssl, ssl(bifactor))
  expect_equal(ssl(bifactor)[["2"]],
               unname(colSums(bifactor$loadings[["2"]]^2)),
               tolerance = 1e-12)
})
