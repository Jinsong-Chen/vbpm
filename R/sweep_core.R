## Internal model-neutral runner for consecutive factor-count candidates.
## Public adapters own validation, model construction, field names, warnings,
## and selection semantics.

.run_sweep <- function(K, fit_one, extract_one, failure_one,
                       progress_one = NULL) {
  if (!is.integer(K) || length(K) < 1L || anyNA(K) ||
      (length(K) > 1L && any(diff(K) != 1L))) {
    stop("Internal sweep K must be consecutive increasing integers.",
         call. = FALSE)
  }
  stopifnot(is.function(fit_one), is.function(extract_one),
            is.function(failure_one))

  rows <- vector("list", length(K))
  loading <- vector("list", length(K))
  pip <- vector("list", length(K))

  for (i in seq_along(K)) {
    k <- K[i]
    started <- proc.time()[["elapsed"]]
    elapsed <- NA_real_
    fit <- tryCatch(fit_one(k), error = identity)

    value <- if (inherits(fit, "error")) {
      elapsed <- max(0, proc.time()[["elapsed"]] - started)
      failure_one(k, elapsed, fit)
    } else {
      tryCatch(
        {
          extracted <- extract_one(fit, k)
          elapsed <- max(0, proc.time()[["elapsed"]] - started)
          extracted
        },
        error = function(e) {
          elapsed <<- max(0, proc.time()[["elapsed"]] - started)
          failure_one(k, elapsed, e)
        }
      )
    }

    valid_K <- is.list(value) && is.data.frame(value$row) &&
      nrow(value$row) == 1L && "K" %in% names(value$row) &&
      is.numeric(value$row$K) && !is.complex(value$row$K) &&
      length(value$row$K) == 1L && !is.na(value$row$K) &&
      is.finite(value$row$K) && value$row$K == floor(value$row$K) &&
      identical(as.integer(value$row$K), k)
    valid_matrices <- is.list(value) && is.matrix(value$loading) &&
      is.matrix(value$pip) && identical(dim(value$loading), dim(value$pip))
    if (!valid_K || !valid_matrices) {
      stop("Internal sweep adapter returned a malformed candidate payload.",
           call. = FALSE)
    }
    value$row$secs <- as.numeric(elapsed)
    rows[[i]] <- value$row
    loading[[i]] <- value$loading
    pip[[i]] <- value$pip
    if (is.function(progress_one)) progress_one(value$row)
  }

  sweep <- do.call(rbind, rows)
  rownames(sweep) <- NULL
  candidate_names <- as.character(sweep$K)
  names(loading) <- candidate_names
  names(pip) <- candidate_names

  list(sweep = sweep, loadings = loading, pips = pip)
}
