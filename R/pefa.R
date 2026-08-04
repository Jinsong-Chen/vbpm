## Partially exploratory factor analysis: an in-memory factor-count sweep.

#' Gain rule for choosing the number of factors
#'
#' Selects a factor count by comparing each marginal gain with a percentage of
#' the largest gain on the same fixed criterion path. The rule is invariant to
#' positive affine transformations of that path; candidates must otherwise use
#' the same data, priors, tuning, and optimization setup.
#'
#' A candidate is selected only when all `sustain` gains following it are
#' strictly below the threshold. The complete look-ahead is required: a short
#' suffix at the right edge does not qualify.
#'
#' @param K Consecutive candidate factor counts (any order).
#' @param score A finite criterion for which larger is better. Pass the negative
#'   of an information criterion for which smaller is better.
#' @param delta Threshold as a percentage of the largest marginal gain.
#' @param sustain Number of consecutive sub-threshold gains required.
#'
#' @return One factor count from `K`.
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' @examples
#' elbo <- c(-5400, -5100, -4900, -4895, -4892)
#' select_K_elbow(2:6, elbo)
#'
#' bic <- c(10900, 10500, 10310, 10395, 10480)
#' select_K_elbow(2:6, -bic)
#'
#' jagged <- c(-5400, -5100, -5095, -4900, -4897, -4895)
#' select_K_elbow(2:7, jagged, sustain = 1)
#' select_K_elbow(2:7, jagged, sustain = 2)
#'
#' @seealso [pefa()]
#' @export
select_K_elbow <- function(K, score, delta = 10, sustain = 2) {
  if (length(K) != length(score)) {
    stop("K and score must have the same length (got ", length(K), " and ",
         length(score), ").", call. = FALSE)
  }
  if (!length(K)) stop("K must not be empty.", call. = FALSE)
  if (!is.numeric(K) || anyNA(K) || any(!is.finite(K)))
    stop("K must be finite and numeric.", call. = FALSE)
  if (any(K != round(K)))
    stop("K must contain whole numbers.", call. = FALSE)
  if (anyDuplicated(K))
    stop("K must not contain duplicates.", call. = FALSE)
  if (!is.numeric(score) || anyNA(score) || any(!is.finite(score)))
    stop("score must be finite and numeric.", call. = FALSE)
  if (length(delta) != 1L || !is.numeric(delta) || is.na(delta) ||
      !is.finite(delta) || delta < 0 || delta > 100) {
    stop("delta must be a single number in [0, 100].", call. = FALSE)
  }
  sustain <- .pefa_sustain(sustain)

  o <- order(K)
  K <- as.integer(K[o])
  score <- score[o]
  if (length(K) > 1L && any(diff(K) != 1L))
    stop("K must contain consecutive whole numbers.", call. = FALSE)
  .scan_gain(K, diff(score), delta, sustain)$selected_K
}

#' Partially exploratory factor analysis with gain-rule selection
#'
#' Fits [vbfa()] over a consecutive factor-count window while holding the
#' partially specified backbone `Q0` fixed and adding exploratory columns.
#' Automatic selection uses the named ELBO-gain cuts only. AIC, BIC, absolute
#' fit indices, BIC gains, and adjacent structural summaries are descriptive.
#'
#' Ordinary `pefa()` candidates are oblique and have diagonal residuals. The
#' optional bifactor sweep is retained as a research extension: its reported
#' `K` is the number of group factors and [vbfa()] adds one general factor.
#' For directly comparable adjacent-count summaries, bifactor transitions
#' remove that labelled general column and apply the ordinary matcher to the
#' remaining group-loading and group-PIP blocks. Whether the bifactor extension
#' recovers a useful group-factor count is an empirical question rather than a
#' guarantee of the paper-defined ordinary workflow.
#'
#' @param Q0 A `J` by `K0` backbone matrix with entries `1`, `0`, or `-1`.
#' @param Y An `N` by `J` data matrix.
#' @param Kmin,Kmax Inclusive factor-count window, with `Kmin >= ncol(Q0)`.
#' @param cuts A named numeric vector of thresholds in `[0, 100]`. Names must
#'   be unique after trimming and exactly one must be `primary`.
#' @param sustain Number of consecutive sub-threshold gains required.
#' @param bifactor Logical; fit the labelled bifactor extension. In this mode
#'   `Kmin` and `Kmax` count group factors.
#' @param general General-column design passed to [vbfa()] in bifactor mode: a
#'   scalar or length-`J` vector in `-1`, `0`, `1`.
#' @param max_it Maximum iterations per [vbfa()] path stage.
#' @param tau Hard-selection threshold passed to [fit_stats()].
#' @param stability_eps Positive squared-L2-norm cutoff for retaining a
#'   non-backbone column in adjacent loading comparisons. It is not a loading
#'   magnitude, PIP, or delivery threshold.
#' @param rank_adjust Logical; request Jacobian-rank parameter counts for every
#'   candidate. The deterministic default uses nominal counts.
#' @param rank_max_J Positive whole-number size guard for explicit Jacobian-rank
#'   counting.
#' @param verbose Logical; print per-candidate progress and named results.
#' @param ... Named estimation controls passed to [vbfa()]. PEFA controls its
#'   own model class: `ld`, `Qe`, `ld_control`, and `orthogonal` are rejected.
#'
#' @return An object of class `pefa` with two data frames and two matrix lists:
#'   \describe{
#'     \item{`$sweep`}{Candidate properties: `K`, `ELBO`, `AIC`, `BIC`,
#'       `RMSEA`, `SRMR`, `CFI`, `TLI`, `t`, `iter`, `secs`, and `converged`.}
#'     \item{`$transitions`}{Adjacent edges: `K_from`, `K_to`, `ELBO_gain`,
#'       `BIC_gain`, `rmsd`, `rmsd_max`, `phi_min`, `ari`, `pip_rmsd`, the
#'       logical `collision` flag, and `unmatched_max`.}
#'     \item{`$loadings`, `$pips`}{Named lists containing the natural loading
#'       and PIP matrix for every candidate `K`.}
#'   }
#'   In order, the other top-level fields are named `selected_K` and `boundary`
#'   vectors, the backbone `Q0`, and `settings`. The settings are `cuts`,
#'   `sustain`, `bifactor`, `general`, `tau`, `stability_eps`, and
#'   `rank_adjust`. Full candidate fits are not retained.
#'
#' @section Adjacent loading and PIP summaries:
#' Every compared loading column is oriented to have a nonnegative mean. In a
#' bifactor sweep the labelled general column is first removed, so the same
#' definitions below operate on group factors in both modes. The first
#' `ncol(Q0)` backbone columns retain their positions and bypass screening;
#' remaining columns are screened by `stability_eps`, then each smaller-fit
#' column is independently matched to its minimum-SSE larger-fit column. `rmsd`
#' pools common cells; `rmsd_max` is the largest
#' column-wise RMSD; and `phi_min` is signed Tucker congruence. For bipolar or
#' contrast columns with means at or near zero, reflection orientation is weak,
#' so interpret `phi_min` beside both RMSDs. `ari` compares dominant-loading
#' item partitions using all retained columns, including an unmatched larger-
#' fit column. `collision` is `TRUE` when two or more retained non-backbone
#' columns in the smaller fit select the same retained larger-fit column; this
#' reuse is permitted by the independent matching rule.
#'
#' `pip_rmsd` is the continuous PIP RMSD over matched cells that are
#' unspecified/regularized in both adjacent candidates; fixed backbone cells
#' are excluded. `unmatched_max` is the largest absolute loading outside
#' the completely matched non-backbone core, including screened-out columns; it
#' is descriptive and is not an added-factor or salience rule. The complete
#' bifactor loading and PIP matrices, including the general column, remain in
#' `$loadings` and `$pips` for candidate-level inspection.
#'
#' @section Complete paths and boundaries:
#' ELBO-gain selection requires every candidate to have converged and the whole
#' ELBO-gain path to be finite. An incomplete path returns same-named
#' `NA_integer_` selections, `"none"` boundaries, and one warning; candidates
#' on opposite sides of a gap are never joined. An upper result calls for a
#' rerun with larger `Kmax`. At a lower result, reduce `Kmin` when `Kmin > K0`;
#' when `Kmin == K0`, reduce or revise `Q0`. A one-candidate window has boundary
#' `"single"`.
#'
#' @section Comparability and refitting:
#' Positive affine transformations of one fixed criterion path leave the gain
#' rule unchanged. This does not make sweeps with different data,
#' preprocessing, priors, hyperparameters, or optimization settings
#' comparable. The object deliberately has no checkpoint/resume behavior. To
#' inspect `Phi`, residuals, scores, or other complete-fit quantities, rebuild
#' the selected design from `Q0` and refit it with [vbfa()].
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
#' for (k in 1:2) {
#'   anchors <- which(groups == k)[1:2]
#'   Q0[anchors, k] <- 1L
#' }
#' result <- pefa(Q0, sim$dat, 2, 4, verbose = FALSE)
#' result$selected_K
#' result$transitions
#' K_selected <- result$selected_K[["primary"]]
#' result$loadings[[as.character(K_selected)]]
#' }
#'
#' @seealso [vbfa()], [fit_stats()], [select_K_elbow()]
#' @export
pefa <- function(Q0, Y, Kmin, Kmax,
                 cuts = c(primary = 10), sustain = 2,
                 bifactor = FALSE, general = 1,
                 max_it = 10000, tau = 0.50,
                 stability_eps = 0.1,
                 rank_adjust = FALSE, rank_max_J = 100,
                 verbose = TRUE, ...) {
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names) || anyNA(dot_names) ||
        any(!nzchar(trimws(dot_names)))) {
      stop("Every argument in ... must be named.", call. = FALSE)
    }
    if (anyDuplicated(dot_names))
      stop("Arguments in ... must have unique names.", call. = FALSE)

    forbidden <- c("ld", "Qe", "ld_control", "orthogonal")
    changes_mode <- vapply(dot_names, function(nm) {
      any(startsWith(forbidden, nm))
    }, logical(1))
    if (any(changes_mode)) {
      stop("pefa() does not allow mode controls in ...: ",
           paste(dot_names[changes_mode], collapse = ", "),
           ". Use vbfa() directly for orthogonal or local-dependence fits.",
           call. = FALSE)
    }

    controlled <- c("Y", "Q", "bifactor", "general", "max_it")
    collisions <- vapply(dot_names, function(nm) {
      any(startsWith(controlled, nm))
    }, logical(1))
    if (any(collisions)) {
      stop("pefa() controls these vbfa() arguments: ",
           paste(dot_names[collisions], collapse = ", "), ".", call. = FALSE)
    }
  }

  cuts <- .pefa_cuts(cuts)
  sustain <- .pefa_sustain(sustain)
  .validate_rank_controls(rank_adjust, rank_max_J)

  Y <- as.matrix(Y)
  Q0 <- as.matrix(Q0)
  J <- ncol(Y)
  K0 <- ncol(Q0)
  if (!length(dim(Y)) || nrow(Y) < 1L || J < 1L)
    stop("Y must be a non-empty N by J matrix.", call. = FALSE)
  if (nrow(Q0) != J)
    stop("nrow(Q0) must equal ncol(Y).", call. = FALSE)
  if (!all(Q0 %in% c(-1, 0, 1)))
    stop("Q0 entries must be in {-1, 0, 1}.", call. = FALSE)

  for (nm in c("Kmin", "Kmax")) {
    value <- get(nm)
    if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
        !is.finite(value) || value < 1 || value != round(value)) {
      stop(nm, " must be a single positive whole number.", call. = FALSE)
    }
  }
  Kmin <- as.integer(Kmin)
  Kmax <- as.integer(Kmax)
  if (Kmax < Kmin) {
    stop(sprintf("Kmax (%d) < Kmin (%d): the window must be non-decreasing.",
                 Kmax, Kmin), call. = FALSE)
  }
  if (Kmin < K0) {
    stop(sprintf(
      "Kmin (%d) < K0 (%d): the window must be at least the backbone size.",
      Kmin, K0), call. = FALSE)
  }
  if (length(max_it) != 1L || !is.numeric(max_it) || is.na(max_it) ||
      !is.finite(max_it) || max_it < 1 || max_it != round(max_it)) {
    stop("max_it must be a single positive whole number.", call. = FALSE)
  }
  max_it <- as.integer(max_it)
  if (length(tau) != 1L || !is.numeric(tau) || is.na(tau) ||
      !is.finite(tau) || tau < 0 || tau > 1) {
    stop("tau must be a single number in [0, 1].", call. = FALSE)
  }
  if (length(stability_eps) != 1L || !is.numeric(stability_eps) ||
      is.complex(stability_eps) || is.na(stability_eps) ||
      !is.finite(stability_eps) || stability_eps <= 0) {
    stop("stability_eps must be one finite, real, positive number.",
         call. = FALSE)
  }
  if (length(bifactor) != 1L || !is.logical(bifactor) || is.na(bifactor))
    stop("bifactor must be TRUE or FALSE.", call. = FALSE)
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose))
    stop("verbose must be TRUE or FALSE.", call. = FALSE)

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
      stop("general must be a scalar or length-J vector containing only -1, 0, or 1.",
           call. = FALSE)
    }
    general_design <- rep_len(as.integer(general), J)
    if (!any(general_design == 1L))
      stop("general must specify (1) at least one general loading.", call. = FALSE)
    if (K0 > 0L && any(apply(Q0 == 1, 2L, all))) {
      stop("With bifactor = TRUE, Q0 must contain only the group backbone; ",
           "the general column is added internally.", call. = FALSE)
    }
  }

  Ks <- seq.int(Kmin, Kmax)
  pad <- function(K) {
    if (K == K0) Q0 else cbind(Q0, matrix(-1L, J, K - K0))
  }
  failed_payload <- function(K, elapsed, error) {
    total_K <- K + as.integer(bifactor)
    list(
      row = data.frame(
        K = as.integer(K), ELBO = NA_real_, AIC = NA_real_, BIC = NA_real_,
        RMSEA = NA_real_, SRMR = NA_real_, CFI = NA_real_, TLI = NA_real_,
        t = NA_real_, iter = NA_integer_, secs = as.numeric(elapsed),
        converged = FALSE, stringsAsFactors = FALSE
      ),
      loading = matrix(NA_real_, J, total_K),
      pip = matrix(NA_real_, J, total_K)
    )
  }
  fit_one <- function(K) {
    args <- list(Y = Y, Q = pad(K), max_it = max_it)
    if (isTRUE(bifactor)) {
      args <- c(args, list(bifactor = TRUE, general = general_design))
    }
    withCallingHandlers(
      do.call(vbfa, c(args, dots)),
      warning = function(w) {
        if (startsWith(conditionMessage(w),
                       "Maximum number of MFVB iterations exceeded")) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  extract_one <- function(fit, K) {
    fs <- fit_stats(
      fit, Y = Y, Q = fit$Q, tau = tau,
      rank_adjust = rank_adjust, rank_max_J = rank_max_J
    )
    list(
      row = data.frame(
        K = as.integer(K),
        ELBO = { e <- as.numeric(fit$ELBO); if (length(e)) e[1L] else NA_real_ },
        AIC = unname(fs[["AIC"]]), BIC = unname(fs[["BIC"]]),
        RMSEA = unname(fs[["RMSEA"]]), SRMR = unname(fs[["SRMR"]]),
        CFI = unname(fs[["CFI"]]), TLI = unname(fs[["TLI"]]),
        t = unname(fs[["t"]]), iter = as.integer(fit$iter),
        secs = NA_real_, converged = isTRUE(fit$converged),
        stringsAsFactors = FALSE
      ),
      loading = as.matrix(fit$Lam),
      pip = as.matrix(fit$pi)
    )
  }
  progress_one <- if (isTRUE(verbose)) function(row) {
    cat(sprintf("K=%2d  ELBO=%.1f  BIC=%.0f  CFI=%.3f  conv=%d  %.1fs\n",
                row$K, row$ELBO, row$BIC, row$CFI,
                as.integer(row$converged), row$secs))
  } else NULL

  run <- .run_sweep(
    Ks, fit_one = fit_one, extract_one = extract_one,
    failure_one = failed_payload, progress_one = progress_one
  )
  sweep <- run$sweep
  loadings <- run$loadings
  pips <- run$pips
  transitions <- .build_transitions(
    sweep, loadings, pips, Q0 = Q0,
    stability_eps = stability_eps, bifactor = bifactor
  )

  complete <- isTRUE(all(sweep$converged & is.finite(sweep$ELBO))) &&
    nrow(transitions) == max(0L, length(Ks) - 1L) &&
    all(is.finite(transitions$ELBO_gain))
  selection <- .sweep_gain_selections(
    Ks, transitions$ELBO_gain, cuts, sustain, complete = complete
  )
  if (!complete) {
    failed_K <- Ks[!sweep$converged | !is.finite(sweep$ELBO)]
    detail <- if (length(failed_K)) {
      paste0(" (K = ", paste(failed_K, collapse = ", "), ")")
    } else ""
    warning("PEFA ELBO-gain selection is unavailable because its candidate ",
            "path contains a failed or non-finite fit", detail,
            "; repair or rerun the complete sweep.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    for (nm in names(cuts)) {
      shown <- if (is.na(selection$selected_K[[nm]])) "NA" else
        as.character(selection$selected_K[[nm]])
      cat(sprintf("PEFA ELBO gain [%s = %g%%]: K = %s (%s)\n",
                  nm, cuts[[nm]], shown, selection$boundary[[nm]]))
    }
  }

  out <- list(
    sweep = sweep,
    transitions = transitions,
    loadings = loadings,
    pips = pips,
    selected_K = selection$selected_K,
    boundary = selection$boundary,
    Q0 = Q0,
    settings = list(
      cuts = cuts,
      sustain = sustain,
      bifactor = bifactor,
      general = general_design,
      tau = tau,
      stability_eps = stability_eps,
      rank_adjust = rank_adjust
    )
  )
  class(out) <- "pefa"
  out
}

#' Methods for a PEFA factor-number sweep
#'
#' Methods consume the compact sweep and transition tables directly. Loading
#' and PIP matrices remain available in separate named lists on the original
#' object but are omitted from summaries and console candidate tables.
#'
#' @param object,x A result returned by [pefa()].
#' @param ... Further arguments. For `plot()`, graphical arguments passed to
#'   the first base plot in the selected display.
#' @param type Plot type: criterion trajectories (`"objective"`), one gain
#'   trajectory (`"gain"`), or descriptive absolute-fit trajectories
#'   (`"fit"`).
#' @param criterion Criterion used when `type = "gain"`: `"elbo"` (default)
#'   or `"bic"`. Only the ELBO panel draws automatic cut thresholds and the
#'   primary selection; the BIC-gain panel is descriptive.
#' @param digits Number of decimals used when printing fit and stability
#'   summaries.
#'
#' @return `summary()` returns a compact `summary.pefa` object without loading
#'   or PIP matrices. It carries `ssl`, the per-candidate sums of squared
#'   loadings (`colSums(Lam^2)`), derived from `$loadings` at summary time
#'   rather than stored on the sweep object; in bifactor mode each vector
#'   covers all `K + 1` columns with the general factor first. SS loadings
#'   are descriptive column sizes, not Gram eigenvalues, and take part in no
#'   selection rule. Print and plot methods return their input invisibly.
#'
#' @name pefa-methods
NULL

.print_pefa_rules <- function(x, indent = "  ") {
  cuts <- x$settings$cuts
  for (nm in names(cuts)) {
    selected <- x$selected_K[[nm]]
    shown <- if (is.na(selected)) "NA" else as.character(selected)
    cat(sprintf("%sELBO gain [%s = %g%%] -> K = %-3s (%s)\n",
                indent, nm, cuts[[nm]], shown, x$boundary[[nm]]))
  }
  invisible(NULL)
}

.pefa_window <- function(x) {
  c(Kmin = min(x$sweep$K), Kmax = max(x$sweep$K), K0 = ncol(x$Q0))
}

.pefa_boundary_advice <- function(x) {
  window <- .pefa_window(x)
  upper <- names(x$boundary)[x$boundary == "upper"]
  lower <- names(x$boundary)[x$boundary == "lower"]
  if (length(upper)) {
    cat(sprintf(
      "  Upper-window selection [%s]: rerun the complete sweep with a larger Kmax.\n",
      paste(upper, collapse = ", ")
    ))
  }
  if (length(lower)) {
    minimum_K <- max(1L, window[["K0"]])
    if (window[["Kmin"]] > minimum_K) {
      cat(sprintf(
        "  Lower-window selection [%s]: rerun the complete sweep with a smaller Kmin.\n",
        paste(lower, collapse = ", ")
      ))
    } else if (window[["K0"]] == 0L) {
      cat(sprintf(
        "  Lower selection [%s] is already at the minimum supported K; revise the search or model.\n",
        paste(lower, collapse = ", ")
      ))
    } else {
      cat(sprintf(
        "  Lower selection [%s] equals the backbone size; reduce or revise Q0 to search lower.\n",
        paste(lower, collapse = ", ")
      ))
    }
  }
  invisible(NULL)
}

#' @rdname pefa-methods
#' @export
print.pefa <- function(x, ...) {
  window <- .pefa_window(x)
  mode <- if (isTRUE(x$settings$bifactor))
    " group factors (+ 1 general)" else ""
  cat(sprintf("PEFA sweep: K = %d:%d%s; %d of %d converged\n",
              window[["Kmin"]], window[["Kmax"]], mode,
              sum(x$sweep$converged), nrow(x$sweep)))
  .print_pefa_rules(x)
  .pefa_boundary_advice(x)
  invisible(x)
}

#' @rdname pefa-methods
#' @export
summary.pefa <- function(object, ...) {
  out <- list(
    window = .pefa_window(object),
    sweep = object$sweep,
    transitions = object$transitions,
    ## Derived at summary time from $loadings; never stored on the object.
    ## In bifactor mode each vector covers all K + 1 columns, general first.
    ssl = lapply(object$loadings, function(m) colSums(m^2)),
    selected_K = object$selected_K,
    boundary = object$boundary,
    settings = object$settings,
    convergence = c(converged = sum(object$sweep$converged),
                    total = nrow(object$sweep))
  )
  class(out) <- "summary.pefa"
  out
}

#' @rdname pefa-methods
#' @export
print.summary.pefa <- function(x, digits = 3, ...) {
  mode <- if (isTRUE(x$settings$bifactor))
    " group factors (+ 1 general)" else ""
  cat(sprintf("PEFA sweep summary: K = %d:%d%s; sustain = %d\n",
              x$window[["Kmin"]], x$window[["Kmax"]], mode,
              x$settings$sustain))
  .print_pefa_rules(x)

  cat("\nCandidates:\n")
  candidates <- x$sweep
  if (isTRUE(x$settings$bifactor)) {
    candidates <- cbind(
      K = candidates$K, K_total = candidates$K + 1L,
      candidates[, setdiff(names(candidates), "K"), drop = FALSE]
    )
  }
  for (nm in intersect(c("ELBO", "AIC", "BIC"), names(candidates)))
    candidates[[nm]] <- round(candidates[[nm]], 1L)
  for (nm in intersect(c("RMSEA", "SRMR", "CFI", "TLI"), names(candidates)))
    candidates[[nm]] <- round(candidates[[nm]], digits)
  selected <- x$selected_K[["primary"]]
  candidates$selected <- ifelse(
    !is.na(selected) & candidates$K == selected, "*", ""
  )
  print(candidates, row.names = FALSE)

  cat("\nSS loadings (colSums(Lam^2)",
      if (isTRUE(x$settings$bifactor)) "; general column first",
      "):\n", sep = "")
  for (nm in names(x$ssl)) {
    cat(sprintf("  K=%s: %s\n", nm,
                paste(sprintf("%.2f", x$ssl[[nm]]), collapse = " ")))
  }

  if (nrow(x$transitions)) {
    cat("\nTransitions:\n")
    transitions <- x$transitions
    numeric_cols <- vapply(transitions, is.numeric, logical(1))
    transitions[numeric_cols] <- lapply(transitions[numeric_cols], round, digits)
    print(transitions, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname pefa-methods
#' @importFrom graphics abline axis legend lines mtext par plot
#' @export
plot.pefa <- function(x, type = c("objective", "gain", "fit"),
                      criterion = c("elbo", "bic"), ...) {
  type <- match.arg(type)
  criterion <- match.arg(criterion)
  d <- x$sweep
  K <- d$K
  selected <- x$selected_K[["primary"]]
  col1 <- "#1f77b4"
  col2 <- "#d62728"
  factor_xlab <- if (isTRUE(x$settings$bifactor))
    "Number of group factors (K)" else "Number of factors (K)"

  padded_range <- function(values, top = 0.20) {
    values <- values[is.finite(values)]
    if (!length(values)) return(c(0, 1))
    r <- range(values)
    span <- diff(r)
    if (span == 0) span <- abs(r[1L]) * 0.01 + 1
    c(r[1L] - 0.05 * span, r[2L] + top * span)
  }
  mark_selected <- function(value = selected) {
    if (length(value) == 1L && !is.na(value))
      abline(v = value, lty = 3, col = "grey40")
  }

  if (type == "objective") {
    op <- par(mar = c(5, 4.2, 4, 4.2) + 0.1)
    on.exit(par(op))
    plot(K, d$ELBO, type = "b", pch = 16, col = col1,
         ylim = padded_range(d$ELBO), xaxt = "n",
         xlab = factor_xlab,
         ylab = "ELBO (higher is better)", main = "PEFA criteria", ...)
    axis(1, at = K)
    mark_selected()
    par(new = TRUE)
    plot(K, d$BIC, type = "b", pch = 17, lty = 2, col = col2,
         ylim = padded_range(d$BIC), axes = FALSE, xlab = "", ylab = "")
    axis(4, col = col2, col.axis = col2)
    mtext("BIC (lower is better)", side = 4, line = 2.8, col = col2)
    legend("top", c("ELBO (left axis)", "BIC (right axis)"),
           col = c(col1, col2), lty = c(1, 2), pch = c(16, 17),
           bty = "n", cex = 0.9)
  } else if (type == "gain") {
    transition_col <- if (criterion == "elbo") "ELBO_gain" else "BIC_gain"
    gains <- x$transitions[[transition_col]]
    thresholds <- numeric(0)
    if (criterion == "elbo" && length(gains) && all(is.finite(gains)))
      thresholds <- unname(x$settings$cuts / 100 * max(gains))
    y_values <- c(gains, thresholds, 0)
    if (nrow(x$transitions)) {
      plot(x$transitions$K_to, gains, type = "b", pch = 16, col = col1,
           xlim = range(K), ylim = padded_range(y_values), xaxt = "n",
           xlab = factor_xlab, ylab = "Marginal gain",
           main = paste(toupper(criterion), "gain sweep"), ...)
    } else {
      plot(K, rep(NA_real_, length(K)), type = "n",
           xlim = range(K), ylim = padded_range(y_values), xaxt = "n",
           xlab = factor_xlab, ylab = "Marginal gain",
           main = paste(toupper(criterion), "gain sweep"), ...)
    }
    axis(1, at = K)
    abline(h = 0, col = "grey85")
    if (length(thresholds))
      for (value in thresholds) abline(h = value, lty = 3, col = "grey40")
    if (criterion == "elbo") mark_selected()

    labels <- toupper(paste(criterion, "gain"))
    lty <- 1
    pch <- 16
    cols <- col1
    if (length(thresholds)) {
      labels <- c(labels,
                  sprintf("%s = %g", names(x$settings$cuts), thresholds))
      lty <- c(lty, rep(3, length(thresholds)))
      pch <- c(pch, rep(NA, length(thresholds)))
      cols <- c(cols, rep("grey40", length(thresholds)))
    }
    legend("top", labels, col = cols, lty = lty, pch = pch,
           bty = "n", cex = 0.9)
  } else {
    op <- par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1) + 0.1,
              oma = c(0, 0, 2, 0))
    on.exit(par(op))
    plot(K, d$RMSEA, type = "b", pch = 16, col = col1,
         ylim = padded_range(c(d$RMSEA, d$SRMR)), xaxt = "n",
         xlab = factor_xlab, ylab = "Index (lower is better)",
         main = "Absolute misfit", ...)
    axis(1, at = K)
    lines(K, d$SRMR, type = "b", pch = 17, lty = 2, col = col2)
    mark_selected()
    legend("top", c("RMSEA", "SRMR"), col = c(col1, col2),
           lty = c(1, 2), pch = c(16, 17), bty = "n", cex = 0.9)

    plot(K, d$CFI, type = "b", pch = 16, col = col1,
         ylim = padded_range(c(d$CFI, d$TLI)), xaxt = "n",
         xlab = factor_xlab, ylab = "Index (higher is better)",
         main = "Incremental fit")
    axis(1, at = K)
    lines(K, d$TLI, type = "b", pch = 17, lty = 2, col = col2)
    mark_selected()
    legend("bottom", c("CFI", "TLI"), col = c(col1, col2),
           lty = c(1, 2), pch = c(16, 17), bty = "n", cex = 0.9)
    mtext("PEFA descriptive fit trajectories", outer = TRUE, font = 2)
  }
  invisible(x)
}
