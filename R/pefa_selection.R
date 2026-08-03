## Selection rules for PEFA sweeps.

## Validate and normalize named sensitivity cuts.
.pefa_cuts <- function(cuts) {
  if (!is.numeric(cuts) || !length(cuts) || anyNA(cuts) ||
      any(!is.finite(cuts)))
    stop("cuts must be a non-empty finite numeric vector.", call. = FALSE)
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
  if (length(sustain) != 1L || !is.numeric(sustain) || !is.finite(sustain) ||
      sustain < 1 || sustain != round(sustain))
    stop("sustain must be a single positive whole number.", call. = FALSE)
  as.integer(sustain)
}

## Shared implementation of the paper's sustained-drop rule. The caller
## supplies a complete finite path and larger-is-better adjacent gains.
.scan_gain <- function(K, gain, cut, sustain) {
  if (length(gain) != max(0L, length(K) - 1L) || any(!is.finite(gain)))
    stop("gain must be a complete finite adjacent path.", call. = FALSE)
  if (length(K) == 1L)
    return(list(selected_K = as.integer(K), gain_max = NA_real_,
                threshold = NA_real_))

  gmax <- max(gain)
  threshold <- unname(cut) / 100 * gmax
  if (!any(gain > 0))
    return(list(selected_K = as.integer(K[1L]), gain_max = gmax,
                threshold = threshold))

  selected <- as.integer(K[length(K)])
  last_start <- length(gain) - sustain + 1L
  if (last_start >= 1L) {
    for (i in seq_len(last_start)) {
      ## Strict inequality and full look-ahead are both part of the rule.
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

## Evaluate ELBO and BIC raw/gain rules on one complete candidate path.
.pefa_selection <- function(sweep, transitions, cuts, sustain) {
  K <- sweep$K
  specs <- list(
    elbo = list(node = "ELBO", gain = "ELBO_gain", raw = which.max),
    bic  = list(node = "BIC",  gain = "BIC_gain",  raw = which.min)
  )
  rows <- list()
  unavailable <- character(0)

  for (criterion in names(specs)) {
    spec <- specs[[criterion]]
    node <- sweep[[spec$node]]
    gain <- transitions[[spec$gain]]
    complete <- isTRUE(all(sweep$converged & is.finite(node))) &&
      length(gain) == max(0L, length(K) - 1L) &&
      !any(!is.finite(gain))

    raw_K <- if (complete) as.integer(K[spec$raw(node)]) else NA_integer_
    rows[[length(rows) + 1L]] <- data.frame(
      criterion = criterion,
      form = "raw",
      cut_name = NA_character_,
      cut = NA_real_,
      gain_max = NA_real_,
      threshold = NA_real_,
      selected_K = raw_K,
      boundary = .pefa_boundary(raw_K, K),
      stringsAsFactors = FALSE
    )

    if (!complete)
      unavailable <- c(unavailable, toupper(criterion))

    for (cut_name in names(cuts)) {
      ans <- if (complete) {
        .scan_gain(K, gain, cuts[[cut_name]], sustain)
      } else {
        list(selected_K = NA_integer_, gain_max = NA_real_,
             threshold = NA_real_)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        criterion = criterion,
        form = "gain",
        cut_name = cut_name,
        cut = unname(cuts[[cut_name]]),
        gain_max = ans$gain_max,
        threshold = ans$threshold,
        selected_K = ans$selected_K,
        boundary = .pefa_boundary(ans$selected_K, K),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(unavailable))
    warning("PEFA ", paste(unique(unavailable), collapse = " and "),
            " selection is unavailable because its candidate path contains ",
            "a failed or non-finite fit; repair or rerun the full sweep.",
            call. = FALSE)
  do.call(rbind, rows)
}
