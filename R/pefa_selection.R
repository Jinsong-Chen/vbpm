## Selection rules for pefa() sweeps (0.8.0).
##
## Skip-gap policy: a gain edge is usable only when both endpoints are eligible
## in the scope being evaluated and its raw gain is finite. Non-adjacent gains
## are never constructed, so an ineligible interior candidate breaks the chain
## instead of being bridged: differencing across a gap would report a
## two-factor jump as a one-factor marginal gain.

## Normalize the `cuts`/`delta` pair into a named numeric vector.
.pefa_cuts <- function(cuts, delta, cuts_given, delta_given) {
  if (cuts_given && delta_given)
    stop("Supply either `cuts` or the deprecated `delta`, not both.",
         call. = FALSE)
  if (delta_given && !is.null(delta)) cuts <- c(primary = delta)
  if (!is.numeric(cuts) || anyNA(cuts) || any(!is.finite(cuts)))
    stop("cuts must be a finite numeric vector.", call. = FALSE)
  if (any(cuts < 0 | cuts > 100))
    stop("cuts must lie in [0, 100].", call. = FALSE)
  nm <- names(cuts)
  if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("cuts must have unique, non-empty names.", call. = FALSE)
  if (sum(nm == "primary") != 1L)
    stop("cuts must contain exactly one entry named 'primary'.", call. = FALSE)
  cuts
}

## One rule: scan usable adjacent edges, resetting the sustained-run counter
## wherever the chain is broken.
.scan_gain <- function(Kf, gain, usable, thr, sustain) {
  n <- length(gain)
  run <- 0L; start <- NA_integer_
  for (i in seq_len(n)) {
    if (!usable[i]) { run <- 0L; start <- NA_integer_; next }
    if (gain[i] < thr) {
      if (run == 0L) start <- Kf[i]
      run <- run + 1L
      if (run >= sustain) return(start)
    } else { run <- 0L; start <- NA_integer_ }
  }
  NA_integer_
}

## Evaluate every criterion x form x cut x scope combination.
.pefa_selection <- function(sweep, transitions, cuts, sustain, types) {
  scopes <- list(`converged+fit` = sweep$eligible, converged = sweep$converged)
  crits  <- list(objective = "Objective", elbo = "ELBO", bic = "BIC")
  rows <- list()
  unavailable <- character(0)

  for (sc in names(scopes)) {
    ok <- scopes[[sc]]
    for (cr in names(crits)) {
      node_col <- crits[[cr]]
      node_val <- sweep[[node_col]]
      ## orient so larger is better at the node level
      node_score <- if (cr == "bic") -node_val else node_val
      fin <- ok & is.finite(node_score)
      gain_col <- paste0(node_col, "_gain")
      gain <- transitions[[gain_col]]
      usable <- logical(nrow(transitions))
      if (nrow(transitions))
        usable <- is.finite(gain) &
          fin[match(transitions$K_from, sweep$K)] &
          fin[match(transitions$K_to,   sweep$K)]
      gmax <- if (any(usable)) max(gain[usable]) else NA_real_

      ## ---- raw row -----------------------------------------------------
      raw_K <- if (any(fin)) sweep$K[fin][which.max(node_score[fin])] else NA_integer_
      rows[[length(rows) + 1L]] <- data.frame(
        criterion = cr, form = "raw", cut = NA_character_,
        cut_value = NA_real_, eligibility_scope = sc,
        gain_max = NA_real_, threshold = NA_real_,
        selected_K = as.integer(raw_K), stringsAsFactors = FALSE)
      if (!any(fin)) unavailable <- c(unavailable, sprintf("%s/raw/%s", cr, sc))

      ## ---- gain rows, one per named cut --------------------------------
      for (cn in names(cuts)) {
        sel <- NA_integer_; thr <- NA_real_
        if (!any(fin)) {
          unavailable <- c(unavailable, sprintf("%s/gain[%s]/%s", cr, cn, sc))
        } else if (!any(usable)) {
          sel <- max(sweep$K[fin])                       # no usable edge
        } else if (!is.finite(gmax) || gmax <= 0) {
          sel <- min(sweep$K[fin])                       # nothing worth adding
        } else {
          thr <- (cuts[[cn]] / 100) * gmax
          sel <- .scan_gain(transitions$K_from, gain, usable, thr, sustain)
          if (is.na(sel)) sel <- max(sweep$K[fin])       # never sustained
        }
        rows[[length(rows) + 1L]] <- data.frame(
          criterion = cr, form = "gain", cut = cn,
          cut_value = unname(cuts[[cn]]), eligibility_scope = sc,
          gain_max = if (any(usable)) gmax else NA_real_,
          threshold = thr, selected_K = as.integer(sel),
          stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)

  ## boundary readouts: against the requested window, and against the nodes
  ## actually usable for that row
  Kmin <- min(sweep$K); Kmax <- max(sweep$K)
  wb <- function(k) if (is.na(k)) "none" else if (k == Kmin) "lower" else
                    if (k == Kmax) "upper" else "interior"
  out$window_boundary <- vapply(out$selected_K, wb, character(1))
  out$support_boundary <- vapply(seq_len(nrow(out)), function(i) {
    k <- out$selected_K[i]
    if (is.na(k)) return("none")
    ok <- if (out$eligibility_scope[i] == "converged") sweep$converged else sweep$eligible
    cr <- out$criterion[i]
    v  <- sweep[[crits[[cr]]]]
    sup <- sweep$K[ok & is.finite(v)]
    if (!length(sup)) "none"
    else if (length(sup) == 1L) "single"
    else if (k == min(sup)) "lower"
    else if (k == max(sup)) "upper"
    else "interior"
  }, character(1))
  attr(out, "unavailable") <- unique(unavailable)
  out
}

## The one row that drives the top-level result.
.pefa_primary_row <- function(selection) {
  which(selection$criterion == "objective" & selection$form == "gain" &
        selection$cut == "primary" &
        selection$eligibility_scope == "converged+fit")[1L]
}
