## select_K_elbow(): the only packaged count calculator, and one pefa() never
## invokes.  Its whole contract is the rule plus the NA cases, so both are
## pinned here against hand-computed paths.

test_that("the convenience defaults are 20 and 1", {
  expect_identical(formals(select_K_elbow)$delta, 20)
  expect_identical(formals(select_K_elbow)$sustain, 1)
  ## Gains 100, 5, 3 against a 20% cut of 100: K = 3 is the first crossing.
  path <- cumsum(c(0, 100, 5, 3))
  expect_identical(select_K_elbow(2:5, path), 3L)
  expect_identical(select_K_elbow(2:5, path, delta = 20, sustain = 1), 3L)
})

test_that("the rule returns the smallest qualifying count", {
  ## Gains 100, 5, 3: candidates 3 and 4 both qualify at 20%, so the smaller
  ## one is returned.
  path <- cumsum(c(0, 100, 5, 3))
  expect_identical(select_K_elbow(2:5, path, delta = 20, sustain = 1), 3L)
  expect_identical(select_K_elbow(2:5, path, delta = 20, sustain = 2), 3L)
  ## A larger cut moves the reading up the same path; a smaller one moves it
  ## down, and neither is a package preference.
  expect_identical(select_K_elbow(2:5, path, delta = 4, sustain = 1), 4L)
  expect_identical(select_K_elbow(2:5, path, delta = 100, sustain = 1), 3L)
})

test_that("the comparison is strict and the look-ahead is complete", {
  ## Gains 100, 50, 20, 5 at a 20% cut of 100: the 20 gain is NOT below 20.
  path <- cumsum(c(0, 100, 50, 20, 5))
  expect_identical(select_K_elbow(2:6, path, delta = 20, sustain = 1), 5L)

  ## sustain = 2 needs two complete gains after the candidate; the window
  ## here contains only one, so nothing qualifies.
  expect_identical(select_K_elbow(4:6, c(0, 100, 105), delta = 20,
                                  sustain = 2), NA_integer_)
  ## The same path with room for the second gain does qualify.
  expect_identical(select_K_elbow(4:7, c(0, 100, 105, 107), delta = 20,
                                  sustain = 2), 5L)
  ## A jagged path separates the two sustain settings.
  jagged <- c(-5400, -5100, -5095, -4900, -4897, -4895)
  expect_identical(select_K_elbow(2:7, jagged, delta = 20, sustain = 1), 3L)
  expect_identical(select_K_elbow(2:7, jagged, delta = 20, sustain = 2), 5L)
})

test_that("every no-answer case returns typed NA_integer_", {
  ## A one-candidate window holds no gain at all.
  expect_identical(select_K_elbow(3L, -4900), NA_integer_)
  expect_identical(select_K_elbow(3L, -4900, delta = 0, sustain = 5),
                   NA_integer_)
  ## Delta_max <= 0: a percentage of a nonpositive largest gain is no cut.
  expect_identical(select_K_elbow(2:4, c(10, 9, 8)), NA_integer_)
  expect_identical(select_K_elbow(2:4, c(10, 10, 10)), NA_integer_)
  ## No start index has a complete run of qualifying gains after it.
  expect_identical(select_K_elbow(2:5, cumsum(c(0, 100, 90, 95)),
                                  delta = 20, sustain = 1), NA_integer_)
  ## delta = 0 makes the strict inequality unsatisfiable for a positive gain.
  expect_identical(select_K_elbow(2:5, cumsum(c(0, 100, 5, 3)), delta = 0,
                                  sustain = 1), NA_integer_)
})

test_that("valid input is sorted by K before the gains are computed", {
  K <- c(5L, 2L, 4L, 3L)
  score <- c(105, 0, 100, 95)
  expect_identical(select_K_elbow(K, score, delta = 20, sustain = 1),
                   select_K_elbow(sort(K), score[order(K)], delta = 20,
                                  sustain = 1))
  expect_identical(select_K_elbow(K, score, delta = 20, sustain = 1), 3L)
  ## The score is oriented larger-is-better, so an IC path is negated.
  bic <- c(10900, 10500, 10310, 10395, 10480)
  expect_identical(select_K_elbow(2:6, -bic, delta = 20, sustain = 1), 4L)
})

test_that("malformed paths are errors, never NA", {
  expect_error(select_K_elbow(2:4, c(1, 2)), "same length")
  expect_error(select_K_elbow(numeric(0), numeric(0)), "must not be empty")
  expect_error(select_K_elbow(c(2, 2, 3), c(1, 2, 3)), "duplicates")
  expect_error(select_K_elbow(c(2, 4, 6), c(1, 2, 3)), "consecutive")
  expect_error(select_K_elbow(c(2, 3, 4.5), c(1, 2, 3)), "whole numbers")
  expect_error(select_K_elbow(c(2, 3, Inf), c(1, 2, 3)), "K must be finite")
  expect_error(select_K_elbow(2:4, c(1, 2, NA)), "score must be finite")
  expect_error(select_K_elbow(2:4, c(1, 2, Inf)), "score must be finite")
  expect_error(select_K_elbow(matrix(2:5, 2), matrix(1:4, 2)),
               "numeric vectors")
  expect_error(select_K_elbow(as.character(2:4), c(1, 2, 3)),
               "numeric vectors")

  path <- cumsum(c(0, 100, 5, 3))
  expect_error(select_K_elbow(2:5, path, delta = -1), "\\[0, 100\\]")
  expect_error(select_K_elbow(2:5, path, delta = 101), "\\[0, 100\\]")
  expect_error(select_K_elbow(2:5, path, delta = c(10, 20)), "\\[0, 100\\]")
  expect_error(select_K_elbow(2:5, path, delta = NA_real_), "\\[0, 100\\]")
  expect_error(select_K_elbow(2:5, path, sustain = 0), "positive whole")
  expect_error(select_K_elbow(2:5, path, sustain = -1), "positive whole")
  expect_error(select_K_elbow(2:5, path, sustain = 1.5), "positive whole")
  expect_error(select_K_elbow(2:5, path, sustain = c(1, 2)), "positive whole")
})

test_that("oversized integer-backed input errors before any coercion", {
  path <- cumsum(c(0, 100, 5, 3))
  ## as.integer(3e10) is NA_integer_ plus a coercion warning.  Before the range
  ## check that NA reached the look-ahead and surfaced as "missing value where
  ## TRUE/FALSE needed" -- a downstream failure, not an argument error.
  expect_error(select_K_elbow(2:5, path, sustain = 3e10),
               "positive whole number")
  expect_error(select_K_elbow(2:5, path,
                              sustain = .Machine$integer.max + 1),
               "positive whole number")
  ## Controlled error and no warning: the coercion never happens.
  for (bad in list(3e10, .Machine$integer.max + 1)) {
    expect_silent(tryCatch(select_K_elbow(2:5, path, sustain = bad),
                           error = function(e) NULL))
  }
  ## The representable boundary is admissible; it simply finds no answer,
  ## because no start index has that many gains after it.
  expect_identical(select_K_elbow(2:5, path,
                                  sustain = .Machine$integer.max),
                   NA_integer_)

  ## K carries the same pre-coercion rule, and likewise never warns.
  expect_error(select_K_elbow(3e10 + 0:3, path), "representable as integers")
  expect_silent(tryCatch(select_K_elbow(3e10 + 0:3, path),
                         error = function(e) NULL))
  expect_error(select_K_elbow(c(2, 3, .Machine$integer.max + 1), path[1:3]),
               "representable as integers")
})

test_that("complex and nonpositive input are rejected by name", {
  path <- cumsum(c(0, 100, 5, 3))
  ## is.numeric() already excludes complex, but the contract asks for the
  ## rejection to be explicit rather than incidental.
  expect_error(select_K_elbow(as.complex(2:5), path), "not complex")
  expect_error(select_K_elbow(2:5, as.complex(path)), "not complex")
  expect_error(select_K_elbow(2:5, path, delta = 20 + 0i), "not complex")
  expect_error(select_K_elbow(2:5, path, sustain = 1 + 0i), "not complex")
  for (bad in list(as.complex(2:5), 2:5)) {
    expect_silent(tryCatch(select_K_elbow(bad, path, sustain = 1 + 0i),
                           error = function(e) NULL))
  }

  ## A factor count is positive; pefa() never sweeps below 1.
  expect_error(select_K_elbow(c(0, 1, 2, 3), path), "positive whole numbers")
  expect_error(select_K_elbow(-1:2, path), "positive whole numbers")
  expect_identical(select_K_elbow(1:4, path), 2L)
})

test_that("pefa() never calls the calculator", {
  code <- deparse(pefa)
  expect_false(any(grepl("select_K_elbow", code, fixed = TRUE)))
})
