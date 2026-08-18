## Partially exploratory factor analysis: an evidence-only factor-count sweep.
##
## pefa() returns evidence and decides nothing.  It selects no factor count,
## applies no threshold, vetoes no comparison, and grades nothing.  Every
## quantity needed to apply a downstream count rule, persistence rule,
## threshold, horizon, or collision policy is present in the returned object;
## turning those facts into a conclusion is the assessor's job.
##
## Component evidence is protected independently: a statistics failure never
## erases a valid loading matrix, a malformed PIP never erases a valid ELBO,
## and a fit that stopped at max_it still contributes its returned matrices.
## Convergence is recorded, never used to suppress a number that is defined.


# ---- canonical naming -------------------------------------------------

## Item names come from the data when usable, otherwise from a fixed pattern.
## A supplied rownames(Q0) must agree exactly rather than silently winning.
.pefa_item_names <- function(Y, Q0, J) {
  nm <- colnames(Y)
  nm <- if (!is.null(nm) && length(nm) == J && !anyNA(nm) &&
            all(nzchar(nm)) && !anyDuplicated(nm)) {
    as.character(nm)
  } else {
    paste0("item_", seq_len(J))
  }
  qn <- rownames(Q0)
  if (!is.null(qn) && !identical(as.character(qn), nm)) {
    stop("rownames(Q0) must match the item names taken from colnames(Y).",
         call. = FALSE)
  }
  nm
}

.pefa_backbone_names <- function(Q0, K0) {
  nm <- colnames(Q0)
  if (!is.null(nm) && length(nm) == K0 && !anyNA(nm) &&
      all(nzchar(nm)) && !anyDuplicated(nm)) {
    return(as.character(nm))
  }
  if (K0 == 0L) character(0) else paste0("F", seq_len(K0))
}

## Group-block factor names for candidate K: the backbone, then the appended
## exploratory columns continuing the F-numbering from K0 + 1.
.pefa_factor_names <- function(backbone_names, K0, K) {
  extra <- if (K > K0) paste0("F", seq.int(K0 + 1L, K)) else character(0)
  c(backbone_names, extra)
}

.pefa_dimnames <- function(item_names, backbone_names, K0, K, bifactor) {
  cols <- .pefa_factor_names(backbone_names, K0, K)
  if (isTRUE(bifactor)) cols <- c("G", cols)
  list(item_names, cols)
}


# ---- component status codebook (section 3.1.1) ------------------------

.pefa_na_matrix <- function(dn) {
  matrix(NA_real_, length(dn[[1L]]), length(dn[[2L]]), dimnames = dn)
}

.pefa_ok <- function(value) {
  list(status = "ok", reason = NA_character_, value = value)
}

.pefa_bad <- function(reason, value) {
  list(status = "unavailable", reason = reason, value = value)
}

## A recognized finite scalar flag that does not indicate convergence is
## `nonconverged`; anything unrecognizable is `malformed_fit`.  The two are
## different evidence and must not collapse into one code.
.pefa_check_fit <- function(fit) {
  if (inherits(fit, "error")) {
    return(list(status = "unavailable", reason = "fit_error",
                converged = NA))
  }
  flag <- if (is.list(fit)) fit[["flag"]] else NULL
  ok_flag <- !is.null(flag) && length(flag) == 1L && is.numeric(flag) &&
    !is.complex(flag) && !is.na(flag) && is.finite(flag)
  if (!ok_flag) {
    return(list(status = "unavailable", reason = "malformed_fit",
                converged = NA))
  }
  if (flag == 1) {
    list(status = "ok", reason = NA_character_, converged = TRUE)
  } else {
    list(status = "unavailable", reason = "nonconverged", converged = FALSE)
  }
}

## Schema first, then finiteness: a wrong-shaped object cannot have its cells
## inspected, so `nonfinite_*` applies only once the shape is known good.
.pefa_check_matrix <- function(x, dn, kind, fit_failed) {
  na_out <- .pefa_na_matrix(dn)
  codes <- if (kind == "loading") {
    c(missing = "missing_loading", malformed = "malformed_loading",
      nonfinite = "nonfinite_loading")
  } else {
    c(missing = "missing_pip", malformed = "malformed_pip",
      nonfinite = "nonfinite_pip")
  }
  if (isTRUE(fit_failed)) return(.pefa_bad("upstream_fit_error", na_out))
  if (is.null(x)) return(.pefa_bad(codes[["missing"]], na_out))

  shape_ok <- is.matrix(x) && is.numeric(x) && !is.complex(x) &&
    identical(dim(x), c(length(dn[[1L]]), length(dn[[2L]])))
  if (!shape_ok) return(.pefa_bad(codes[["malformed"]], na_out))

  value <- matrix(as.double(x), nrow(x), ncol(x), dimnames = dn)
  bad <- !is.finite(value)
  if (any(bad)) {
    value[bad] <- NA_real_
    return(.pefa_bad(codes[["nonfinite"]], value))
  }
  if (kind == "pip") {
    out_of_range <- value < 0 | value > 1
    if (any(out_of_range)) {
      value[out_of_range] <- NA_real_
      return(.pefa_bad("out_of_range_pip", value))
    }
  }
  .pefa_ok(value)
}

.pefa_check_elbo <- function(fit, fit_failed) {
  if (isTRUE(fit_failed)) return(.pefa_bad("upstream_fit_error", NA_real_))
  value <- if (is.list(fit)) fit[["ELBO"]] else NULL
  if (is.null(value)) return(.pefa_bad("missing_elbo", NA_real_))
  if (length(value) != 1L || !is.numeric(value) || is.complex(value)) {
    return(.pefa_bad("malformed_elbo", NA_real_))
  }
  value <- as.double(value)
  if (!is.finite(value)) return(.pefa_bad("nonfinite_elbo", NA_real_))
  .pefa_ok(value)
}

.pefa_stats_fields <- c("t_nom", "t", "t_S", "AIC", "BIC", "AIC_S", "BIC_S",
                        "RMSEA", "SRMR", "CFI", "TLI")
.pefa_criterion_fields <- c("t_nom", "t", "t_S", "AIC", "BIC", "AIC_S", "BIC_S")

.pefa_na_stats <- function() {
  stats::setNames(rep(NA_real_, length(.pefa_stats_fields)),
                  .pefa_stats_fields)
}

## RMSEA/SRMR/CFI/TLI are descriptive: a legitimate NA_real_ (undefined degrees
## of freedom, typically) does not invalidate the criterion path.  Only the
## counts and information criteria can set `invalid_criterion_stats`.
.pefa_check_stats <- function(fs, fit_failed) {
  na_out <- .pefa_na_stats()
  if (isTRUE(fit_failed)) return(.pefa_bad("upstream_fit_error", na_out))
  if (inherits(fs, "error")) return(.pefa_bad("stats_error", na_out))

  nm <- names(fs)
  usable <- is.numeric(fs) && !is.complex(fs) && !is.null(nm) &&
    all(.pefa_stats_fields %in% nm) && !anyDuplicated(nm[nm %in% .pefa_stats_fields])
  if (!usable) return(.pefa_bad("malformed_stats", na_out))

  value <- vapply(.pefa_stats_fields,
                  function(f) as.double(unname(fs[[f]])), numeric(1))
  descriptive <- setdiff(.pefa_stats_fields, .pefa_criterion_fields)
  value[descriptive][!is.finite(value[descriptive])] <- NA_real_

  crit <- value[.pefa_criterion_fields]
  whole_ok <- function(x) is.finite(x) && x >= 0 && x == floor(x)
  invalid <- any(!is.finite(crit)) ||
    !whole_ok(value[["t_nom"]]) || !whole_ok(value[["t"]]) ||
    !is.finite(value[["t_S"]]) || value[["t_S"]] < 0
  if (invalid) {
    value[.pefa_criterion_fields][!is.finite(crit)] <- NA_real_
    return(.pefa_bad("invalid_criterion_stats", value))
  }
  .pefa_ok(value)
}

## Smallest SSL among the candidate's K0 position-matched backbone columns.
## In bifactor mode those are stored columns 2:(K0 + 1); the general column is
## never part of the backbone scale.
.pefa_backbone_ssl_min <- function(lam, K0, bifactor, loading_ok) {
  if (K0 == 0L || !isTRUE(loading_ok)) return(NA_real_)
  offset <- as.integer(isTRUE(bifactor))
  cols <- seq_len(K0) + offset
  min(colSums(lam[, cols, drop = FALSE]^2))
}


# ---- one candidate ----------------------------------------------------

.pefa_candidate <- function(fit, fs, K, dn, K0, bifactor, elapsed) {
  fit_c <- .pefa_check_fit(fit)
  fit_failed <- identical(fit_c$reason, "fit_error")

  lam_c <- .pefa_check_matrix(if (is.list(fit)) fit[["Lam"]] else NULL,
                              dn, "loading", fit_failed)
  pip_c <- .pefa_check_matrix(if (is.list(fit)) fit[["pi"]] else NULL,
                              dn, "pip", fit_failed)
  elbo_c <- .pefa_check_elbo(fit, fit_failed)
  stats_c <- .pefa_check_stats(fs, fit_failed)

  iter <- if (!fit_failed && is.list(fit)) fit[["iter"]] else NULL
  iter <- if (!is.null(iter) && length(iter) == 1L && is.numeric(iter) &&
              !is.na(iter) && is.finite(iter)) as.integer(iter) else NA_integer_

  sv <- stats_c$value
  row <- data.frame(
    K = as.integer(K),
    ELBO = elbo_c$value,
    t_nom = sv[["t_nom"]], t = sv[["t"]], t_S = sv[["t_S"]],
    AIC = sv[["AIC"]], BIC = sv[["BIC"]],
    AIC_S = sv[["AIC_S"]], BIC_S = sv[["BIC_S"]],
    RMSEA = sv[["RMSEA"]], SRMR = sv[["SRMR"]],
    CFI = sv[["CFI"]], TLI = sv[["TLI"]],
    iter = iter,
    secs = as.double(elapsed),
    converged = fit_c$converged,
    backbone_ssl_min = .pefa_backbone_ssl_min(
      lam_c$value, K0, bifactor, identical(lam_c$status, "ok")
    ),
    fit_status = fit_c$status, fit_reason = fit_c$reason,
    loading_status = lam_c$status, loading_reason = lam_c$reason,
    pip_status = pip_c$status, pip_reason = pip_c$reason,
    elbo_status = elbo_c$status, elbo_reason = elbo_c$reason,
    stats_status = stats_c$status, stats_reason = stats_c$reason,
    stringsAsFactors = FALSE
  )
  list(row = row, loading = lam_c$value, pip = pip_c$value)
}


# ---- table builders ---------------------------------------------------

## Percent-of-largest-gain, guarded so that one bad candidate does not blank
## the column and an all-nonfinite path does not reach max(numeric(0)).
.pefa_gain_pct <- function(gain) {
  finite <- gain[is.finite(gain)]
  if (!length(finite)) return(rep(NA_real_, length(gain)))
  gmax <- max(finite)
  if (!(gmax > 0)) return(rep(NA_real_, length(gain)))
  ## Divide first.  For gains near the largest finite double, 100 * gain
  ## overflows to Inf before the division ever happens, blanking a ratio that
  ## is perfectly representable; gain / gmax is bounded by 1 in magnitude.
  out <- 100 * (gain / gmax)
  out[!is.finite(out)] <- NA_real_
  out
}

.pefa_transitions_table <- function(sweep) {
  n <- nrow(sweep) - 1L
  if (n <= 0L) {
    empty <- data.frame(K_from = integer(0), K_to = integer(0),
                        stringsAsFactors = FALSE)
    for (nm in c("ELBO_gain", "AIC_gain", "BIC_gain", "AIC_S_gain",
                 "BIC_S_gain", "ELBO_gain_pct", "AIC_gain_pct",
                 "BIC_gain_pct", "AIC_S_gain_pct", "BIC_S_gain_pct")) {
      empty[[nm]] <- numeric(0)
    }
    return(empty)
  }
  from <- seq_len(n)
  to <- from + 1L
  out <- data.frame(K_from = as.integer(sweep$K[from]),
                    K_to = as.integer(sweep$K[to]),
                    stringsAsFactors = FALSE)

  ## ELBO is larger-is-better; each IC gain is from - to, so positive always
  ## favors the larger candidate on every retained path.
  gains <- list(
    ELBO = as.double(sweep$ELBO[to] - sweep$ELBO[from]),
    AIC = as.double(sweep$AIC[from] - sweep$AIC[to]),
    BIC = as.double(sweep$BIC[from] - sweep$BIC[to]),
    AIC_S = as.double(sweep$AIC_S[from] - sweep$AIC_S[to]),
    BIC_S = as.double(sweep$BIC_S[from] - sweep$BIC_S[to])
  )
  for (nm in names(gains)) {
    g <- gains[[nm]]
    g[!is.finite(g)] <- NA_real_
    out[[paste0(nm, "_gain")]] <- g
  }
  for (nm in names(gains)) {
    out[[paste0(nm, "_gain_pct")]] <- .pefa_gain_pct(out[[paste0(nm, "_gain")]])
  }
  out
}

.pefa_persistence_cols <- c(
  "K_from", "K_to", "r", "n_pairs", "phi_min", "rmsd", "rmsd_max", "ari",
  "pip_rmsd", "unmatched_n", "unmatched_ssl", "collision", "n_collisions",
  "collision_targets", "collision_multiplicities", "source_ineligible_n",
  "target_eligible_n", "backbone_degenerate", "from_fit_ok", "to_fit_ok",
  "pair_status", "pair_reason"
)

.pefa_empty_persistence <- function() {
  data.frame(
    K_from = integer(0), K_to = integer(0), r = integer(0),
    n_pairs = integer(0), phi_min = numeric(0), rmsd = numeric(0),
    rmsd_max = numeric(0), ari = numeric(0), pip_rmsd = numeric(0),
    unmatched_n = integer(0), unmatched_ssl = numeric(0),
    collision = logical(0), n_collisions = integer(0),
    collision_targets = character(0),
    collision_multiplicities = character(0),
    source_ineligible_n = integer(0), target_eligible_n = integer(0),
    backbone_degenerate = logical(0), from_fit_ok = logical(0),
    to_fit_ok = logical(0), pair_status = character(0),
    pair_reason = character(0),
    stringsAsFactors = FALSE
  )
}

## Every ordered pair k < l, computed directly from its two endpoints and never
## by chaining adjacent assignments.
.pefa_persistence_table <- function(sweep, loadings, pips, Q0, K0,
                                    stability_eps, bifactor) {
  Ks <- as.integer(sweep$K)
  W <- length(Ks)
  if (W < 2L) return(.pefa_empty_persistence())

  offset <- as.integer(isTRUE(bifactor))
  group_block <- function(m) {
    if (offset == 0L) m else m[, -1L, drop = FALSE]
  }
  reg_mask <- function(K) {
    Q <- if (K > K0) cbind(Q0, matrix(-1L, nrow(Q0), K - K0)) else Q0
    matrix(Q == -1L, nrow(Q0), K)
  }
  ok <- function(field, K) identical(sweep[[field]][match(K, Ks)], "ok")
  fit_ok <- function(K) {
    value <- sweep$converged[match(K, Ks)]
    if (is.na(value)) NA else isTRUE(value)
  }

  pairs <- expand.grid(i = seq_len(W), j = seq_len(W))
  pairs <- pairs[pairs$i < pairs$j, , drop = FALSE]
  pairs <- pairs[order(pairs$i, pairs$j), , drop = FALSE]

  rows <- lapply(seq_len(nrow(pairs)), function(p) {
    k <- Ks[pairs$i[p]]
    l <- Ks[pairs$j[p]]
    facts <- .pefa_pair_facts(
      lam_from = group_block(loadings[[as.character(k)]]),
      lam_to = group_block(loadings[[as.character(l)]]),
      pip_from = group_block(pips[[as.character(k)]]),
      pip_to = group_block(pips[[as.character(l)]]),
      reg_from = reg_mask(k), reg_to = reg_mask(l),
      K0 = K0, stability_eps = stability_eps,
      loading_ok_from = ok("loading_status", k),
      loading_ok_to = ok("loading_status", l),
      pip_ok_from = ok("pip_status", k),
      pip_ok_to = ok("pip_status", l)
    )
    facts[[".assign"]] <- NULL
    data.frame(
      K_from = as.integer(k), K_to = as.integer(l),
      r = as.integer(l - k),
      facts,
      from_fit_ok = fit_ok(k), to_fit_ok = fit_ok(l),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, .pefa_persistence_cols, drop = FALSE]
}


# ---- constructor ------------------------------------------------------

#' Partially exploratory factor analysis over a factor-count window
#'
#' Fits [vbfa()] over a consecutive factor-count window while holding the
#' partially specified backbone `Q0` fixed and adding exploratory columns.
#'
#' `pefa()` **returns evidence and decides nothing.** It selects no factor
#' count, applies no threshold, vetoes no comparison, and grades nothing. The
#' object carries every quantity needed to apply a downstream count rule,
#' persistence rule, threshold, horizon, or collision policy; composing those
#' facts into a conclusion is the assessor's job. [select_K_elbow()] is a
#' criterion-generic calculator the user may invoke, never something `pefa()`
#' calls.
#'
#' Ordinary candidates are oblique with diagonal residuals. The optional
#' bifactor sweep is a research extension: its reported `K` is the number of
#' group factors and [vbfa()] adds one general factor, so a candidate stores
#' `K + 1` columns. Structural comparisons strip that labelled general column
#' and apply the matcher to the group block, while `$loadings`, `$pips`, and
#' [ssl()] retain all `K + 1` columns with the general factor first.
#'
#' @section Component evidence:
#' Each candidate carries five independent status/reason pairs — fit, loading,
#' PIP, ELBO, and statistics. A status of `"unavailable"` records that the
#' component's stated condition was not met, not that its payload is missing:
#' a fit that stopped at `max_it` has `fit_status = "unavailable"` with
#' `fit_reason = "nonconverged"` while its finite matrices are stored
#' unchanged. No component is erased because another is unavailable, and
#' convergence never suppresses a structural number that is mathematically
#' defined. `$persistence` records `from_fit_ok`/`to_fit_ok` rather than
#' gating on them; a consumer wanting converged endpoints filters on those
#' columns.
#'
#' @section The evidence object:
#' The returned list holds exactly eight members, in this order: `sweep`,
#' `transitions`, `persistence`, `loadings`, `pips`, `Q0`, `settings`, and
#' `provenance`. There is no `selected_K`, no `boundary`, no stored call and
#' no stored window: the window is `range(x$sweep$K)` and the backbone size
#' is `K0 = ncol(x$Q0)`.
#'
#' \describe{
#'   \item{`$loadings`, `$pips`}{Named lists keyed by `as.character(K)`. Each
#'     element is a double `J` by `C_K` matrix, where
#'     `C_K = K + as.integer(bifactor)`. A missing or malformed matrix
#'     becomes an all-`NA_real_` matrix of that expected shape. A correctly
#'     shaped matrix with nonfinite or out-of-range cells keeps its shape
#'     with only the invalid cells normalized to `NA_real_`, and its
#'     component status stays `"unavailable"`, so a partial endpoint is never
#'     mistaken for a complete one. The numeric cells of a valid finite
#'     matrix returned by a nonconverged fit are retained unchanged.}
#'   \item{`$Q0`}{The canonical integer backbone, `J` by `K0`, carrying the
#'     same item and backbone-factor names.}
#' }
#'
#' Canonical names are part of that payload. Item names are `colnames(Y)`
#' when those are non-missing, nonempty and unique, otherwise `item_1`, ...,
#' `item_J`; a non-`NULL` `rownames(Q0)` must agree exactly or `pefa()`
#' errors. Backbone factor names are nonempty unique `colnames(Q0)` when
#' supplied, otherwise `F1`, ..., `FK0`; appended exploratory factors
#' continue as `F(K0 + 1)`, ..., `FK`. In bifactor mode the general column is
#' named `G` and precedes the group columns. Loading and PIP matrices always
#' share dimnames.
#'
#' @section The candidate sweep:
#' `$sweep` is a data frame with one row per candidate `K`, in increasing
#' order, and these 27 columns:
#'
#' \describe{
#'   \item{`K`}{integer; the candidate factor count. In bifactor mode this is
#'     the number of *group* factors, and the stored matrices have `K + 1`
#'     columns.}
#'   \item{`ELBO`}{double; the variational objective, larger is better.
#'     `NA_real_` whenever `elbo_status` is not `"ok"`.}
#'   \item{`t_nom`, `t`, `t_S`}{double parameter counts from [fit_stats()]:
#'     nominal, rank-adjusted hard, and soft/effective.}
#'   \item{`AIC`, `BIC`, `AIC_S`, `BIC_S`}{double information criteria under
#'     hard (`tau`) and soft selection; smaller is better.}
#'   \item{`RMSEA`, `SRMR`, `CFI`, `TLI`}{double descriptive fit indices. A
#'     mathematically legitimate `NA_real_` — undefined degrees of freedom,
#'     typically — does not by itself make `stats_status` unavailable, and
#'     any `NaN` or infinite index is normalized to `NA_real_`.}
#'   \item{`iter`}{integer iterations used, or `NA_integer_` when the engine
#'     did not report a usable count.}
#'   \item{`secs`}{double, finite and non-negative; wall-clock seconds for
#'     that candidate. Machine-dependent, and therefore the one field
#'     excluded from `evidence_sha256`.}
#'   \item{`converged`}{logical; `TRUE` for a recognized converged fit,
#'     `FALSE` for a recognized completed nonconverged fit, and `NA` after a
#'     fit error or an unrecognizable convergence flag. It is never
#'     reconstructed from matrix or criterion availability.}
#'   \item{`backbone_ssl_min`}{double; the smallest SSL among that
#'     candidate's `K0` position-matched backbone columns — in bifactor mode
#'     the group block's first `K0`, i.e. stored columns `2:(K0 + 1)`, never
#'     the general column. `NA_real_` when `K0 = 0` or
#'     `loading_status != "ok"`. It gives a continuous scale for judging thin
#'     anchors instead of a single cliff.}
#'   \item{`fit_status`, `loading_status`, `pip_status`, `elbo_status`,
#'     `stats_status`}{character, each exactly `"ok"` or `"unavailable"`.}
#'   \item{`fit_reason`, `loading_reason`, `pip_reason`, `elbo_reason`,
#'     `stats_reason`}{character; `NA_character_` if and only if the matching
#'     status is `"ok"`, otherwise one code from the vocabulary below.}
#' }
#'
#' @section Component reason-code vocabulary:
#' Within a component the reasons are ordered: when more than one defect
#' applies, the first listed reason wins. `upstream_fit_error` means the fit
#' itself errored, so that component had no input to inspect.
#'
#' \describe{
#'   \item{fit}{`fit_error`, `malformed_fit`, `nonconverged`. `"ok"` requires
#'     that the engine returned an object whose recognized finite scalar
#'     convergence flag (`flag == 1` for [vbfa()]) normalizes to `TRUE`. A
#'     returned fit lacking such a flag is `malformed_fit`; a recognized flag
#'     that does not indicate convergence is `nonconverged`. The two are
#'     different evidence and never collapse into one code.}
#'   \item{loading}{`upstream_fit_error`, `missing_loading`,
#'     `malformed_loading`, `nonfinite_loading`. `"ok"` requires a real
#'     numeric matrix of the exact expected dimensions with every cell
#'     finite. `malformed_loading` covers wrong type, complex storage, or
#'     wrong dimensions; `nonfinite_loading` applies only once the schema is
#'     known good.}
#'   \item{PIP}{`upstream_fit_error`, `missing_pip`, `malformed_pip`,
#'     `nonfinite_pip`, `out_of_range_pip`. `"ok"` requires a real numeric
#'     matrix of the exact loading dimensions, every cell finite, and every
#'     value in the closed interval `[0, 1]`; a finite value outside it gives
#'     `out_of_range_pip`.}
#'   \item{ELBO}{`upstream_fit_error`, `missing_elbo`, `malformed_elbo`,
#'     `nonfinite_elbo`. `"ok"` requires exactly one real finite numeric
#'     scalar; `malformed_elbo` means anything else.}
#'   \item{statistics}{`upstream_fit_error`, `stats_error`,
#'     `malformed_stats`, `invalid_criterion_stats`. `"ok"` requires that the
#'     protected [fit_stats()] extraction returned the required named real
#'     scalars and that every parameter-count and hard/soft information
#'     criterion is valid. `malformed_stats` means a non-real-numeric result
#'     or absent/duplicated required names. `invalid_criterion_stats` means
#'     any of `t_nom`, `t`, `t_S`, `AIC`, `BIC`, `AIC_S`, `BIC_S` is
#'     nonfinite, `t_nom` or `t` is negative or non-whole-valued, or `t_S` is
#'     negative; independently valid scalars are kept and only the invalid
#'     ones are normalized.}
#' }
#'
#' A valid-shaped engine matrix is read in the engine's positional order and
#' given the canonical dimnames; returned dimname metadata is not a second
#' schema. `summary()` compacts these ten columns into one `status` label
#' using the frozen display precedence
#' `fit -> loading -> pip -> elbo -> stats`, which never changes a stored
#' status or payload.
#'
#' @section Criterion transitions:
#' `$transitions` is a data frame with one row per **adjacent** pair, so
#' `W - 1` rows for a `W`-candidate window, and 12 columns: integer `K_from`
#' and `K_to`, then the five double raw gains `ELBO_gain`, `AIC_gain`,
#' `BIC_gain`, `AIC_S_gain`, `BIC_S_gain` and their five double percentages
#' `ELBO_gain_pct`, `AIC_gain_pct`, `BIC_gain_pct`, `AIC_S_gain_pct`,
#' `BIC_S_gain_pct`.
#'
#' The table is criterion-path only; it holds no structural column, and the
#' underlying counts `t_nom`, `t`, and `t_S` live on `$sweep`. ELBO is
#' oriented larger-is-better, and every information-criterion gain is
#' `IC_from - IC_to`, so a positive value favors the larger candidate on
#' every path. Raw gains are kept alongside the percentages because they are
#' the primitive from which any cut is computed and they carry the absolute
#' scale that percentages destroy.
#'
#' Percentages are `100 * (gain / gain_max)`, where `gain_max` is the maximum
#' over the **finite** gains of that path, and are `NA_real_` when
#' `gain_max <= 0` or no gain on the path is finite. One nonfinite gain
#' blanks only the edges it touches, never the whole column. The division is
#' grouped so that gains near the largest finite double do not overflow
#' before being scaled.
#'
#' @section Structural persistence:
#' `$persistence` is a data frame with one row for **every ordered pair**
#' `k < l` in the window — `W(W-1)/2` rows — ordered by `K_from` ascending
#' then `K_to` ascending, and **including rows whose `pair_status` is
#' `"unavailable"`**: a missing row and an unavailable one are different
#' evidence, and a downstream join needs the row to exist. Every pair is
#' computed directly from its two endpoints, never by chaining adjacent
#' assignments, so the rows at horizon `1` are the adjacent comparisons and the
#' longer horizons are genuine long-horizon measurements.
#'
#' The matcher aligns columns pairwise by sign, `s(a, b) = sign(a'b)`, and
#' scores the pair set `P_kl` made of the `K0` position-matched backbone
#' pairs plus one assigned pair for each eligible exploratory source column.
#' A non-backbone column is **eligible** exactly when `sum(col^2)` is at
#' least `stability_eps` (equality is eligible, compared exactly). The
#' backbone is **never screened**: its columns keep their prespecified
#' positions in both candidates and are always scored, with their thinness
#' reported through `$sweep$backbone_ssl_min` and the pair-level
#' `backbone_degenerate` rather than through exclusion. Each eligible source
#' independently takes the eligible target of largest absolute congruence,
#' ties broken by smaller aligned RMSD then smaller target index; reuse of a
#' target is permitted, detected, and reported, never prevented.
#'
#' \describe{
#'   \item{`K_from`, `K_to`, `r`}{integer; the source count `k`, the target
#'     count `l`, and the horizon `r`, equal to `K_to - K_from`.}
#'   \item{`n_pairs`}{integer size of `P_kl`: `K0` plus the eligible
#'     exploratory sources that could be assigned. `NA_integer_` when the
#'     assignment cannot be computed. It governs the matched-pair metrics
#'     exactly, so a cell with `n_pairs < K_from` did not score every column
#'     of the smaller candidate.}
#'   \item{`phi_min`}{double; the minimum absolute congruence over `P_kl`.
#'     Non-negative by construction.}
#'   \item{`rmsd`}{double; the sign-aligned RMSE pooled over every
#'     item-by-matched-pair cell, `sqrt(sum ||a - s b||^2 / (J * m))` for
#'     `m = n_pairs`. Not an unsigned norm, not an average of column RMSDs,
#'     and not exploratory-only.}
#'   \item{`rmsd_max`}{double; the largest aligned column RMSD over `P_kl`.}
#'   \item{`ari`}{double adjusted Rand index between two **independently**
#'     formed partitions, one per endpoint, each assigning every item to the
#'     column of largest absolute loading among all `K0` backbone columns
#'     plus every eligible non-backbone column at that endpoint — including
#'     eligible targets no source selected. It is deliberately not the
#'     matched subset, and it is clamped to `[-1, 1]`.}
#'   \item{`pip_rmsd`}{double RMSD between the two PIP matrices over the
#'     matched item-factor cells that are unspecified (`Q == -1`) at **both**
#'     endpoints. Fixed cells and unmatched columns are excluded.}
#'   \item{`unmatched_n`, `unmatched_ssl`}{integer count and double maximum
#'     column SSL over the non-backbone target columns that no source
#'     selected, *including columns below the screen*. `unmatched_ssl` is
#'     `0` when that set is empty.}
#'   \item{`collision`, `n_collisions`}{logical flag and integer count of
#'     **reused targets** — three sources landing on one target is one
#'     collision. `collision` is exactly `n_collisions > 0`. Collision is
#'     recorded evidence; it blanks no metric and vetoes no comparison.}
#'   \item{`collision_targets`, `collision_multiplicities`}{character; the
#'     reused **group-block** target-column indices sorted increasingly and
#'     joined by commas, and the reuse counts in that same order, formatted
#'     identically. The two `strsplit(x, ",", fixed = TRUE)` results are
#'     therefore positionally aligned. Both are `""` for a computed clean
#'     edge and `NA_character_` only when the assignment itself could not be
#'     computed.}
#'   \item{`source_ineligible_n`}{integer; exploratory source columns
#'     screened out at `K_from`. Requires only that endpoint's loadings.}
#'   \item{`target_eligible_n`}{integer; eligible exploratory target columns
#'     at `K_to`. Requires only that endpoint's loadings. A value of `0` with
#'     eligible sources present means no exploratory assignment could be
#'     formed and the row retains the backbone-only measurement.}
#'   \item{`backbone_degenerate`}{logical; `FALSE` when `K0 = 0`, `TRUE` when
#'     any inspectable valid endpoint has a finite zero-norm backbone column,
#'     `FALSE` when both valid endpoints establish that none exists, and `NA`
#'     otherwise. It is the exact-zero boundary of `backbone_ssl_min`, and a
#'     `TRUE` leaves `pair_status = "available"`.}
#'   \item{`from_fit_ok`, `to_fit_ok`}{logical endpoint-convergence facts
#'     copied from `$sweep$converged`, hence `NA` after a fit error or a
#'     malformed flag. They are **recorded, not applied**: the matcher does
#'     not consult `fit_status`, and a nonconverged endpoint that returned a
#'     finite well-formed loading matrix yields the same computation. Filter
#'     on these columns to obtain a converged-endpoints-only reading.}
#'   \item{`pair_status`, `pair_reason`}{character; see below.}
#' }
#'
#' @section Pair status and metric-specific computability:
#' `pair_status` takes exactly two values and has exactly two reason codes:
#'
#' \describe{
#'   \item{`"available"` with `pair_reason = NA_character_`}{Both endpoints
#'     have `loading_status == "ok"`, so the matching inputs are usable.}
#'   \item{`"unavailable"` with `pair_reason = "malformed_loading"`}{At least
#'     one endpoint has `loading_status != "ok"`. This is the sole
#'     unavailability reason; nonconvergence, collision, a thin or degenerate
#'     backbone, ineligible sources, and scarce targets never reach it.}
#' }
#'
#' The row-level status says whether matching input is usable, not whether
#' every optional scalar is defined. Each field then follows its own
#' dependency, and a missing dependency gives that field its typed `NA`
#' without blanking an independently defined neighbour:
#'
#' \itemize{
#'   \item endpoint convergence facts are always logical and are `NA` only
#'     for a fit error or a malformed convergence flag;
#'   \item `source_ineligible_n` and `target_eligible_n` need only their own
#'     valid loading endpoint;
#'   \item the assignment, `n_pairs`, the collision fields, and the unmatched
#'     fields need both valid loadings;
#'   \item `phi_min` additionally needs a nonempty `P_kl` and a nonzero norm
#'     for every scored pair, so a zero-norm scored backbone column makes it
#'     `NA_real_` while the aligned `rmsd` and `rmsd_max` stay defined;
#'   \item `rmsd` and `rmsd_max` need a nonempty `P_kl`;
#'   \item `ari` needs both valid loadings and a retained partition column at
#'     each endpoint;
#'   \item `pip_rmsd` additionally needs both PIP endpoints valid and at
#'     least one jointly unspecified matched cell.
#' }
#'
#' For a computed empty `P_kl` the row stores `n_pairs = 0L`, `NA_real_` for
#' `phi_min`, `rmsd`, `rmsd_max` and `pip_rmsd`, `collision = FALSE`,
#' `n_collisions = 0L`, `""` for both collision strings, and the
#' independently computable unmatched facts — and still reports
#' `pair_status = "available"` when both loading matrices are valid.
#'
#' In bifactor mode the labelled general column is stripped before any of
#' this, so every persistence index and count refers to the group block;
#' dimension and finiteness validation happens on the stored `K + 1` matrix
#' beforehand.
#'
#' @section Run settings:
#' `$settings` is a named list of exactly the ten controls that can change a
#' returned number, in this order: `bifactor`, `general`, `v0`, `max_it`,
#' `convChk`, `tolVal`, `tau`, `stability_eps`, `rank_adjust`, and
#' `rank_max_J`. `general` is `NULL` unless `bifactor = TRUE`, in which case
#' it is the resolved length-`J` integer design. `pefa()` accepts no
#' pass-through control absent from this list, so a run is interpretable from
#' the object alone.
#'
#' `verbose` is the one formal deliberately **omitted**: it changes console
#' output and no returned value, so recording it would make two otherwise
#' identical runs hash differently for a purely cosmetic reason.
#'
#' Every integer-backed formal — `Kmin`, `Kmax`, `max_it`, `rank_max_J` — is
#' validated as real, finite, whole-valued and within
#' `.Machine$integer.max` **before** coercion; a complex or out-of-range
#' value is a controlled argument error, never a coercion warning or an
#' accidental `NA_integer_`. The stored settings are then normalized:
#' `rank_max_J` is one unnamed integer, the logical controls are unnamed
#' logicals, and the numeric controls are unnamed doubles. `100` and `100L`
#' therefore produce identical `$settings` and identical hashes, as do named
#' and unnamed spellings of the same value.
#'
#' @section Provenance:
#' `$provenance` is a named list, in this order: `evidence_schema_id`,
#' `engine_semantics_id`, `package_version`, `package_build_id`,
#' `standardized_data_sha256`, `Q0_sha256`, `settings_sha256`,
#' `evidence_sha256`, `lineage_sha256`, `parent_evidence_sha256`, `reused_K`,
#' and `newly_fitted_K`. The data hash covers the engine-scale standardized
#' matrix and its missingness mask, not the raw input's units.
#' `evidence_sha256` covers the three tables, `$loadings`, `$pips`, `$Q0`,
#' `$settings`, the two identity strings, the package version and the three
#' input hashes, and excludes the machine-dependent `$sweep$secs`, itself,
#' and the lineage fields. A root object has
#' `parent_evidence_sha256 = NA_character_`, `reused_K = integer(0)`, and
#' `newly_fitted_K = Kmin:Kmax`. Provenance carries no count rule, criterion
#' preference, cut, threshold, horizon, or collision policy.
#'
#' @section Standalone fit_stats provenance:
#' Inside `pefa()` the resolved `Q`, the data, `tau`, the orthogonality
#' implied by `bifactor`, `rank_adjust`, and the integer-normalized
#' `rank_max_J` are all recorded by `$settings` and the provenance hashes.
#' Called on its own, [fit_stats()] remains a **numeric calculator**: it
#' returns a named numeric vector and deliberately attaches no partial record
#' of its invocation controls, because a partial record would imply
#' self-certification while omitting other result-changing inputs.
#'
#' A provenance-sensitive caller therefore owns that record. Beside the
#' returned vector it stores the resolved `Q` actually used, the data
#' identity, `tau`, whether the fit was orthogonal, `rank_adjust`, and the
#' integer-normalized `rank_max_J`, and its consumer gates inspect that
#' caller-owned record. Do not expect — or add — provenance attributes on the
#' vector [fit_stats()] returns.
#'
#' @param Q0 A `J` by `K0` backbone matrix with entries `1`, `0`, or `-1`.
#'   `K0 = 0` (fully exploratory) is supported.
#' @param Y An `N` by `J` data matrix.
#' @param Kmin,Kmax Inclusive factor-count window, with `Kmin >= ncol(Q0)`.
#' @param bifactor Fit an orthogonal bifactor parameterization.
#' @param general Scalar or length-`J` general-factor design in `{-1, 0, 1}`,
#'   used only when `bifactor = TRUE`.
#' @param v0 Spike variance path passed to [vbfa()].
#' @param max_it,convChk,tolVal Optimizer controls passed to [vbfa()].
#' @param tau Loading-inclusion threshold used by [fit_stats()].
#' @param stability_eps Squared-norm screen for exploratory columns. Must be
#'   strictly positive: at zero a zero-norm column would clear the screen and
#'   produce an undefined congruence.
#' @param rank_adjust,rank_max_J Rank-adjustment controls for [fit_stats()].
#' @param verbose Print one progress line per candidate.
#'
#' @return An object of class `"pefa"`: a list with `sweep`, `transitions`,
#'   `persistence`, `loadings`, `pips`, `Q0`, `settings`, and `provenance`,
#'   in that order. See [persistence()] and [ssl()] for accessors.
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' Chen, J. (2023). Fully and partially exploratory factor analysis with
#' bi-level Bayesian regularization. *Behavior Research Methods*, 55(4),
#' 2125--2142. \doi{10.3758/s13428-022-01884-7}
#'
#' @examples
#' \donttest{
#' sim <- sim_fa(N = 400, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' Q0 <- matrix(-1L, ncol(sim$dat), 2)
#' groups <- rep(1:3, each = 6)
#' for (k in 1:2) Q0[which(groups == k)[1:2], k] <- 1L
#'
#' fit <- pefa(Q0, sim$dat, 2, 4, verbose = FALSE)
#' fit$sweep[, c("K", "ELBO", "BIC", "converged")]
#' fit$transitions[, c("K_from", "K_to", "ELBO_gain", "ELBO_gain_pct")]
#' persistence(fit, "phi_min")
#'
#' ## The count rule is the user's, not the package's.
#' select_K_elbow(fit$sweep$K, fit$sweep$ELBO, delta = 20, sustain = 1)
#' }
#'
#' @seealso [vbfa()], [fit_stats()], [select_K_elbow()], [persistence()],
#'   [ssl()]
#' @export
pefa <- function(Q0, Y, Kmin, Kmax, bifactor = FALSE, general = 1,
                 v0 = c(0.01, 0.005, 0.002, 0.001), max_it = 10000,
                 convChk = FALSE, tolVal = 1e-4, tau = 0.50,
                 stability_eps = 0.10, rank_adjust = FALSE, rank_max_J = 100,
                 verbose = TRUE) {
  .validate_rank_controls(rank_adjust, rank_max_J)

  Y <- as.matrix(Y)
  Q0 <- as.matrix(Q0)
  J <- ncol(Y)
  K0 <- ncol(Q0)
  if (!length(dim(Y)) || nrow(Y) < 1L || J < 1L) {
    stop("Y must be a non-empty N by J matrix.", call. = FALSE)
  }
  if (nrow(Q0) != J) stop("nrow(Q0) must equal ncol(Y).", call. = FALSE)
  if (K0 > 0L && !all(Q0 %in% c(-1, 0, 1))) {
    stop("Q0 entries must be in {-1, 0, 1}.", call. = FALSE)
  }

  ## Every integer-backed formal is checked as real, finite, whole, and inside
  ## the representable R-integer range BEFORE as.integer() sees it: past that
  ## range as.integer() returns NA_integer_ with a coercion warning, which
  ## would turn an argument mistake into a silent failed-fit evidence object.
  ## is.complex() is stated explicitly even though is.numeric() already
  ## excludes complex, because rejecting complex is a contract requirement.
  whole_scalar <- function(value, nm, min_value = 1) {
    ok <- length(value) == 1L && is.numeric(value) && !is.complex(value) &&
      !is.na(value) && is.finite(value) && value >= min_value &&
      value == round(value) && value <= .Machine$integer.max
    if (!ok) {
      stop(nm, " must be a single whole number in [", min_value, ", ",
           .Machine$integer.max, "].", call. = FALSE)
    }
    ## as.integer() also drops names and every other attribute, so a named or
    ## double-typed input cannot change the settings hash for the same run.
    as.integer(value)
  }
  Kmin <- whole_scalar(Kmin, "Kmin")
  Kmax <- whole_scalar(Kmax, "Kmax")
  max_it <- whole_scalar(max_it, "max_it")
  ## Normalized, not merely validated: 100 and 100L must reach $settings and
  ## fit_stats() as one and the same unnamed integer.
  rank_max_J <- whole_scalar(rank_max_J, "rank_max_J")
  if (Kmax < Kmin) {
    stop(sprintf("Kmax (%d) < Kmin (%d): the window must be non-decreasing.",
                 Kmax, Kmin), call. = FALSE)
  }
  if (Kmin < K0) {
    stop(sprintf(
      "Kmin (%d) < K0 (%d): the window must be at least the backbone size.",
      Kmin, K0), call. = FALSE)
  }

  unit_scalar <- function(value, nm, lo, hi) {
    if (length(value) != 1L || !is.numeric(value) || is.complex(value) ||
        is.na(value) || !is.finite(value) || value < lo || value > hi) {
      stop(nm, " must be a single number in [", lo, ", ", hi, "].",
           call. = FALSE)
    }
    as.double(value)
  }
  tau <- unit_scalar(tau, "tau", 0, 1)
  tolVal <- unit_scalar(tolVal, "tolVal", 0, Inf)

  ## Strictly positive: at stability_eps = 0 a zero-norm exploratory column
  ## would clear the screen and give a 0/0 congruence that no pair_reason
  ## covers, since backbone degeneracy is a separate, backbone-only fact.
  if (length(stability_eps) != 1L || !is.numeric(stability_eps) ||
      is.complex(stability_eps) || is.na(stability_eps) ||
      !is.finite(stability_eps) || stability_eps <= 0) {
    stop("stability_eps must be one finite, real, strictly positive number.",
         call. = FALSE)
  }
  stability_eps <- as.double(stability_eps)

  if (!is.numeric(v0) || is.complex(v0) || !length(v0) || anyNA(v0) ||
      any(!is.finite(v0)) || any(v0 <= 0)) {
    stop("v0 must be a non-empty vector of finite positive numbers.",
         call. = FALSE)
  }
  v0 <- as.double(v0)

  flag_scalar <- function(value, nm) {
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    }
    ## as.logical() strips names and any other attribute for the same reason
    ## as whole_scalar(): a stored setting must not carry input decoration.
    as.logical(value)
  }
  bifactor <- flag_scalar(bifactor, "bifactor")
  convChk <- flag_scalar(convChk, "convChk")
  verbose <- flag_scalar(verbose, "verbose")
  ## Already validated by .validate_rank_controls(); this call only normalizes.
  rank_adjust <- flag_scalar(rank_adjust, "rank_adjust")

  if (isTRUE(rank_adjust)) {
    if (J > rank_max_J) {
      stop("rank_adjust = TRUE requires J <= rank_max_J (", J, " > ",
           rank_max_J, "). Raise rank_max_J or use rank_adjust = FALSE.",
           call. = FALSE)
    }
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("rank_adjust = TRUE requires the optional numDeriv package. ",
           "Install numDeriv or use rank_adjust = FALSE.", call. = FALSE)
    }
  }

  general_design <- NULL
  if (isTRUE(bifactor)) {
    if (!is.numeric(general) || is.complex(general) || anyNA(general) ||
        !length(general) || !all(general %in% c(-1, 0, 1)) ||
        !(length(general) %in% c(1L, J))) {
      stop("general must be a scalar or length-J vector of -1, 0, or 1.",
           call. = FALSE)
    }
    general_design <- rep_len(as.integer(general), J)
    if (!any(general_design == 1L)) {
      stop("general must specify (1) at least one general loading.",
           call. = FALSE)
    }
    if (K0 > 0L && any(apply(Q0 == 1, 2L, all))) {
      stop("With bifactor = TRUE, Q0 must contain only the group backbone; ",
           "the general column is added internally.", call. = FALSE)
    }
  }

  item_names <- .pefa_item_names(Y, Q0, J)
  backbone_names <- .pefa_backbone_names(Q0, K0)
  Q0_canon <- matrix(as.integer(Q0), J, K0,
                     dimnames = list(item_names, backbone_names))

  Ks <- seq.int(Kmin, Kmax)
  settings <- list(
    bifactor = bifactor, general = general_design, v0 = v0,
    max_it = max_it, convChk = convChk, tolVal = tolVal, tau = tau,
    stability_eps = stability_eps, rank_adjust = rank_adjust,
    rank_max_J = rank_max_J
  )

  rows <- vector("list", length(Ks))
  loadings <- vector("list", length(Ks))
  pips <- vector("list", length(Ks))

  for (i in seq_along(Ks)) {
    K <- Ks[i]
    dn <- .pefa_dimnames(item_names, backbone_names, K0, K, bifactor)
    Q <- if (K > K0) cbind(Q0_canon, matrix(-1L, J, K - K0)) else Q0_canon
    started <- proc.time()[["elapsed"]]

    args <- list(Y = Y, Q = Q, v0 = v0, max_it = max_it,
                 convChk = convChk, tolVal = tolVal)
    if (isTRUE(bifactor)) {
      args <- c(args, list(bifactor = TRUE, general = general_design))
    }
    ## The per-candidate iteration-limit warning is muffled here because
    ## nonconvergence is a recorded status, not a repeated console event; the
    ## sweep discloses it once, at the end, naming every affected K.
    fit <- tryCatch(
      withCallingHandlers(
        do.call(vbfa, args),
        warning = function(w) {
          if (startsWith(conditionMessage(w),
                         "Maximum number of MFVB iterations exceeded")) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = identity
    )

    ## Statistics extraction is protected separately so that a failure there
    ## cannot erase valid fit, loading, PIP, or ELBO evidence.
    fs <- if (inherits(fit, "error")) fit else tryCatch(
      fit_stats(fit, Y = Y, Q = fit$Q, tau = tau,
                rank_adjust = rank_adjust, rank_max_J = rank_max_J),
      error = identity
    )
    elapsed <- max(0, proc.time()[["elapsed"]] - started)

    candidate <- .pefa_candidate(fit, fs, K, dn, K0, bifactor, elapsed)
    rows[[i]] <- candidate$row
    loadings[[i]] <- candidate$loading
    pips[[i]] <- candidate$pip

    if (isTRUE(verbose)) {
      cat(sprintf("K=%2d  ELBO=%-12s BIC=%-10s CFI=%-7s conv=%-5s %.1fs\n",
                  K, format(candidate$row$ELBO, digits = 7),
                  format(candidate$row$BIC, digits = 7),
                  format(candidate$row$CFI, digits = 3),
                  as.character(candidate$row$converged), elapsed))
    }
  }

  sweep <- do.call(rbind, rows)
  rownames(sweep) <- NULL
  names(loadings) <- as.character(Ks)
  names(pips) <- as.character(Ks)

  evidence <- list(
    sweep = sweep,
    transitions = .pefa_transitions_table(sweep),
    persistence = .pefa_persistence_table(
      sweep, loadings, pips, Q0_canon, K0, stability_eps, bifactor
    ),
    loadings = loadings,
    pips = pips,
    Q0 = Q0_canon,
    settings = settings
  )

  out <- c(evidence, list(provenance = .pefa_provenance(
    evidence = evidence,
    standardized_data_sha256 = .pefa_standardized_hash(Y),
    Q0_sha256 = .pefa_sha256(Q0_canon),
    settings_sha256 = .pefa_sha256(settings),
    package_version = as.character(utils::packageVersion("vbpm")),
    newly_fitted_K = as.integer(Ks)
  )))
  class(out) <- "pefa"

  .pefa_status_warning(sweep)
  out
}

## One aggregate warning per call.  Reason codes only, never raw condition
## text, and never one warning per candidate.
.pefa_status_warning <- function(sweep) {
  components <- c("fit", "loading", "pip", "elbo", "stats")
  hits <- lapply(components, function(cm) {
    bad <- sweep[[paste0(cm, "_status")]] != "ok"
    if (!any(bad)) return(NULL)
    sprintf("%s (%s) at K = %s", cm,
            paste(unique(sweep[[paste0(cm, "_reason")]][bad]), collapse = "/"),
            paste(sweep$K[bad], collapse = ", "))
  })
  hits <- unlist(hits, use.names = FALSE)
  if (!length(hits)) return(invisible(NULL))
  warning("Some PEFA candidate components are unavailable: ",
          paste(hits, collapse = "; "),
          ". Every independently valid payload was retained; see $sweep.",
          call. = FALSE)
  invisible(NULL)
}


# ---- methods ----------------------------------------------------------

#' Methods for a PEFA factor-count sweep
#'
#' @param object,x A result returned by [pefa()].
#' @param ... Further arguments; for `plot()`, graphical arguments passed to
#'   the first base plot in the selected display.
#' @param type Plot type: criterion trajectories (`"objective"`), one gain
#'   trajectory (`"gain"`), or descriptive absolute-fit trajectories
#'   (`"fit"`).
#' @param criterion Criterion path for `type = "objective"` and
#'   `type = "gain"`. Ignored when `type = "fit"`.
#' @param pct Plot percent-of-largest gains instead of raw gains.
#' @param digits Number of decimals used when printing.
#'
#' @return `summary()` returns a compact `summary.pefa` object. Print and plot
#'   methods return their input invisibly. No method announces a selected
#'   factor count, because the package selects none.
#'
#' @name pefa-methods
NULL

.pefa_window <- function(x) {
  c(Kmin = min(x$sweep$K), Kmax = max(x$sweep$K), K0 = ncol(x$Q0))
}

.pefa_status_counts <- function(sweep) {
  components <- c("fit", "loading", "pip", "elbo", "stats")
  vapply(components,
         function(cm) sum(sweep[[paste0(cm, "_status")]] != "ok"),
         integer(1))
}

## Display precedence only.  It never changes a component status or payload.
.pefa_first_failure <- function(sweep) {
  components <- c("fit", "loading", "pip", "elbo", "stats")
  vapply(seq_len(nrow(sweep)), function(i) {
    for (cm in components) {
      if (sweep[[paste0(cm, "_status")]][i] != "ok") return(cm)
    }
    "ok"
  }, character(1))
}

#' @rdname pefa-methods
#' @export
print.pefa <- function(x, ...) {
  w <- .pefa_window(x)
  mode <- if (isTRUE(x$settings$bifactor)) {
    " group factors (+ 1 general)"
  } else ""
  cat(sprintf("PEFA evidence: K = %d:%d%s, K0 = %d; %d candidates\n",
              w[["Kmin"]], w[["Kmax"]], mode, w[["K0"]], nrow(x$sweep)))

  counts <- .pefa_status_counts(x$sweep)
  if (any(counts > 0L)) {
    cat("  unavailable components: ",
        paste(sprintf("%s=%d", names(counts)[counts > 0L],
                      counts[counts > 0L]), collapse = ", "), "\n", sep = "")
  } else {
    cat("  all components ok\n")
  }
  cat(sprintf("  sweep %d x %d | transitions %d x %d | persistence %d x %d\n",
              nrow(x$sweep), ncol(x$sweep),
              nrow(x$transitions), ncol(x$transitions),
              nrow(x$persistence), ncol(x$persistence)))
  cat("  No factor count is selected. Apply your own rule with",
      "select_K_elbow()\n  or read the criterion paths directly;",
      "see vignette(\"pefa\") and vignette(\"bifactor\").\n")
  invisible(x)
}

#' @rdname pefa-methods
#' @export
summary.pefa <- function(object, ...) {
  out <- list(
    window = .pefa_window(object),
    sweep = object$sweep,
    transitions = object$transitions,
    persistence = object$persistence,
    ssl = ssl(object),
    settings = object$settings,
    provenance = object$provenance,
    status_counts = .pefa_status_counts(object$sweep)
  )
  class(out) <- "summary.pefa"
  out
}

#' @rdname pefa-methods
#' @export
print.summary.pefa <- function(x, digits = 3, ...) {
  mode <- if (isTRUE(x$settings$bifactor)) {
    " group factors (+ 1 general)"
  } else ""
  cat(sprintf("PEFA evidence summary: K = %d:%d%s, K0 = %d\n",
              x$window[["Kmin"]], x$window[["Kmax"]], mode,
              x$window[["K0"]]))

  cat("\nCandidates:\n")
  keep <- c("K", "ELBO", "AIC", "BIC", "RMSEA", "SRMR", "CFI", "TLI", "t",
            "iter", "converged")
  candidates <- x$sweep[, intersect(keep, names(x$sweep)), drop = FALSE]
  if (isTRUE(x$settings$bifactor)) {
    candidates <- cbind(K = candidates$K, K_total = candidates$K + 1L,
                        candidates[, setdiff(names(candidates), "K"),
                                   drop = FALSE])
  }
  for (nm in intersect(c("ELBO", "AIC", "BIC"), names(candidates))) {
    candidates[[nm]] <- round(candidates[[nm]], 1L)
  }
  for (nm in intersect(c("RMSEA", "SRMR", "CFI", "TLI"), names(candidates))) {
    candidates[[nm]] <- round(candidates[[nm]], digits)
  }
  candidates$status <- .pefa_first_failure(x$sweep)
  print(candidates, row.names = FALSE)

  cat("\nSS loadings (colSums(Lam^2)",
      if (isTRUE(x$settings$bifactor)) "; general column first" else "",
      "):\n", sep = "")
  for (nm in names(x$ssl)) {
    cat(sprintf("  K=%s: %s\n", nm,
                paste(sprintf("%.2f", x$ssl[[nm]]), collapse = " ")))
  }

  if (nrow(x$transitions)) {
    cat("\nTransitions:\n")
    tr <- x$transitions
    num <- vapply(tr, is.numeric, logical(1))
    tr[num] <- lapply(tr[num], round, digits)
    print(tr, row.names = FALSE)
  }

  cat(sprintf(
    "\nPersistence: %d ordered pairs x %d fields; use persistence(x, metric)\n",
    nrow(x$persistence), ncol(x$persistence)
  ))
  invisible(x)
}

#' @rdname pefa-methods
#' @importFrom graphics abline axis legend lines mtext par plot
#' @export
plot.pefa <- function(x, type = c("objective", "gain", "fit"),
                      criterion = c("elbo", "aic", "bic", "aic_s", "bic_s"),
                      pct = FALSE, ...) {
  type <- match.arg(type)
  criterion <- match.arg(criterion)
  if (length(pct) != 1L || !is.logical(pct) || is.na(pct)) {
    stop("pct must be TRUE or FALSE.", call. = FALSE)
  }
  d <- x$sweep
  K <- d$K
  col1 <- "#1f77b4"
  col2 <- "#d62728"
  xlab <- if (isTRUE(x$settings$bifactor)) {
    "Number of group factors (K)"
  } else "Number of factors (K)"
  path <- c(elbo = "ELBO", aic = "AIC", bic = "BIC",
            aic_s = "AIC_S", bic_s = "BIC_S")[[criterion]]

  padded_range <- function(values, top = 0.20) {
    values <- values[is.finite(values)]
    if (!length(values)) return(c(0, 1))
    r <- range(values)
    span <- diff(r)
    if (span == 0) span <- abs(r[1L]) * 0.01 + 1
    c(r[1L] - 0.05 * span, r[2L] + top * span)
  }

  if (type == "objective") {
    better <- if (criterion == "elbo") "higher is better" else
      "lower is better"
    plot(K, d[[path]], type = "b", pch = 16, col = col1,
         ylim = padded_range(d[[path]]), xaxt = "n", xlab = xlab,
         ylab = paste0(path, " (", better, ")"),
         main = paste("PEFA", path, "path"), ...)
    axis(1, at = K)
  } else if (type == "gain") {
    ## No threshold line and no selection marker: the package chooses none.
    column <- paste0(path, "_gain", if (isTRUE(pct)) "_pct" else "")
    gains <- x$transitions[[column]]
    ylab <- if (isTRUE(pct)) "Gain (% of largest)" else "Marginal gain"
    if (nrow(x$transitions)) {
      plot(x$transitions$K_to, gains, type = "b", pch = 16, col = col1,
           xlim = range(K), ylim = padded_range(c(gains, 0)), xaxt = "n",
           xlab = xlab, ylab = ylab,
           main = paste(path, "gain sweep"), ...)
    } else {
      plot(K, rep(NA_real_, length(K)), type = "n", xlim = range(K),
           ylim = c(0, 1), xaxt = "n", xlab = xlab, ylab = ylab,
           main = paste(path, "gain sweep"), ...)
    }
    axis(1, at = K)
    abline(h = 0, col = "grey85")
  } else {
    op <- par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1) + 0.1,
              oma = c(0, 0, 2, 0))
    on.exit(par(op))
    plot(K, d$RMSEA, type = "b", pch = 16, col = col1,
         ylim = padded_range(c(d$RMSEA, d$SRMR)), xaxt = "n", xlab = xlab,
         ylab = "Index (lower is better)", main = "Absolute misfit", ...)
    axis(1, at = K)
    lines(K, d$SRMR, type = "b", pch = 17, lty = 2, col = col2)
    legend("top", c("RMSEA", "SRMR"), col = c(col1, col2), lty = c(1, 2),
           pch = c(16, 17), bty = "n", cex = 0.9)

    plot(K, d$CFI, type = "b", pch = 16, col = col1,
         ylim = padded_range(c(d$CFI, d$TLI)), xaxt = "n", xlab = xlab,
         ylab = "Index (higher is better)", main = "Incremental fit")
    axis(1, at = K)
    lines(K, d$TLI, type = "b", pch = 17, lty = 2, col = col2)
    legend("bottom", c("CFI", "TLI"), col = c(col1, col2), lty = c(1, 2),
           pch = c(16, 17), bty = "n", cex = 0.9)
    mtext("PEFA descriptive fit trajectories", outer = TRUE, font = 2)
  }
  invisible(x)
}
