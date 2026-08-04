## Internal selection helpers shared by factor-count sweep adapters.

## Validate and normalize named percentage cuts.
.pefa_cuts <- function(cuts) {
  if (!is.numeric(cuts) || !length(cuts) || anyNA(cuts) ||
      any(!is.finite(cuts))) {
    stop("cuts must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (any(cuts < 0 | cuts > 100))
    stop("cuts must lie in [0, 100].", call. = FALSE)

  nm <- names(cuts)
  if (is.null(nm) || anyNA(nm))
    stop("cuts must have unique, non-empty names.", call. = FALSE)
  nm <- trimws(nm)
  if (any(!nzchar(nm)) || anyDuplicated(nm))
    stop("cuts must have unique, non-empty trimmed names.", call. = FALSE)
  if (sum(nm == "primary") != 1L)
    stop("cuts must contain exactly one entry named 'primary'.", call. = FALSE)
  names(cuts) <- nm
  cuts
}

.pefa_sustain <- function(sustain) {
  if (length(sustain) != 1L || !is.numeric(sustain) ||
      is.na(sustain) || !is.finite(sustain) ||
      sustain < 1 || sustain != round(sustain)) {
    stop("sustain must be a single positive whole number.", call. = FALSE)
  }
  as.integer(sustain)
}

## Shared implementation of the sustained-drop rule. The caller supplies a
## complete finite path and larger-is-better adjacent gains.
.scan_gain <- function(K, gain, cut, sustain) {
  if (length(gain) != max(0L, length(K) - 1L) || any(!is.finite(gain)))
    stop("gain must be a complete finite adjacent path.", call. = FALSE)
  if (length(K) == 1L) {
    return(list(selected_K = as.integer(K), gain_max = NA_real_,
                threshold = NA_real_))
  }

  gmax <- max(gain)
  threshold <- unname(cut) / 100 * gmax
  if (!any(gain > 0)) {
    return(list(selected_K = as.integer(K[1L]), gain_max = gmax,
                threshold = threshold))
  }

  selected <- as.integer(K[length(K)])
  last_start <- length(gain) - sustain + 1L
  if (last_start >= 1L) {
    for (i in seq_len(last_start)) {
      ## Strict inequality and full look-ahead are part of the published rule.
      if (all(gain[i:(i + sustain - 1L)] < threshold)) {
        selected <- as.integer(K[i])
        break
      }
    }
  }
  list(selected_K = selected, gain_max = gmax, threshold = threshold)
}

.pefa_boundary <- function(k, K) {
  if (length(k) != 1L || is.na(k)) return("none")
  if (length(K) == 1L) return("single")
  if (k == K[1L]) return("lower")
  if (k == K[length(K)]) return("upper")
  "interior"
}

## Return same-named compact selections and boundaries. Warning policy belongs
## to the public adapter so one incomplete path produces one user-facing warning.
.sweep_gain_selections <- function(K, gain, cuts, sustain, complete = TRUE) {
  nm <- names(cuts)
  if (!isTRUE(complete)) {
    selected <- stats::setNames(rep(NA_integer_, length(cuts)), nm)
    boundary <- stats::setNames(rep("none", length(cuts)), nm)
    return(list(selected_K = selected, boundary = boundary))
  }

  selected <- vapply(cuts, function(cut) {
    .scan_gain(K, gain, cut, sustain)$selected_K
  }, integer(1))
  selected <- stats::setNames(as.integer(selected), nm)
  boundary <- stats::setNames(
    vapply(selected, .pefa_boundary, character(1), K = K), nm
  )
  list(selected_K = selected, boundary = boundary)
}
