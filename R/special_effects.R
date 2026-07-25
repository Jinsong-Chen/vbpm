## Special-effect summaries for orthogonal bifactor fits: the Schmid-Leiman
## proportionality check, the (approximate) higher-order parameterization,
## and testlet/special effect sizes. Pure post-processing of a vbfa fit;
## no re-estimation.

#' Special-effect summaries of a bifactor fit
#'
#' Post-processes an orthogonal bifactor [vbfa()] fit into the three
#' summaries used throughout the `bifactor` vignette: (i) the
#' **Schmid-Leiman (SL) proportionality check** -- the within-group
#' coefficient of variation (CV) of the group-to-general loading ratio, per
#' group factor and item-weighted overall; (ii) the **(approximate)
#' higher-order parameterization** -- second-order loadings
#' \eqn{\hat\gamma_k = 1/\sqrt{1 + \bar r_k^2}} from each group's mean ratio
#' \eqn{\bar r_k}, and first-order loadings
#' \eqn{\hat\lambda_j = |b^{gen}_j| / \hat\gamma_{k(j)}}; and (iii)
#' **testlet/special effect sizes** by Eq. 16 of Zhang and Chen (2024): the
#' average squared loading on each group factor, counting specified loadings
#' always and unspecified ones only when their posterior inclusion
#' probability exceeds `tau`.
#'
#' All output is **descriptive**. Near-zero CVs signal SL proportionality --
#' the fitted bifactor is then (approximately) a higher-order model, and by
#' the same token a testlet model, so the higher-order parameters are a
#' faithful re-expression. Values near `.1` are conventionally read as
#' approximately higher-order, but the cutoff is suggestive, not a test; the
#' higher-order parameters and effect sizes are reported regardless, and
#' with clearly non-constant ratios they are a deliberately lossy summary of
#' a genuinely richer bifactor. The ratio inherits the sampling softness of
#' the general/group split, so CVs are noisy below (roughly) a few thousand
#' observations -- see the `bifactor` vignette for demonstrations, caveats,
#' and the covariance geometry behind all three summaries.
#'
#' **Which fits are accepted.** A fit with `bifactor = TRUE` is used as-is
#' (general column first, group columns after). A legacy hand-built fit is
#' also accepted when it has `orthogonal = TRUE` and exactly one
#' all-specified column in `Q`, which is taken as the general factor. Group
#' membership of each item is the group column where it is anchored
#' (`Q == 1`) or, for unspecified rows, its single active group column
#' (inclusion probability above `tau`). Items with no active group loading
#' are excluded from the ratio summaries and reported; items active on
#' several group columns are assigned to the anchored one (else the largest
#' absolute loading) and reported. Ratios use absolute loadings; items whose
#' general loading is fixed at zero (bifactor-(S-1)-style designs) carry no
#' ratio and are reported.
#'
#' @param object A fitted [vbfa()] object with a bifactor structure (see
#'   Details).
#' @param tau Hard-selection threshold on the posterior inclusion
#'   probability (default `0.5`), used both for group-membership inference
#'   and for the Eq.-16 effect sizes.
#'
#' @return An object of class `special_effects`: a list with
#'   \item{groups}{A data frame with one row per group factor: `n_items`
#'     (assigned items), `mean_ratio`, `sd_ratio`, `cv`, `gamma`
#'     (second-order loading), and `D` (Eq.-16 effect size).}
#'   \item{mean_cv}{The item-weighted mean CV across group factors.}
#'   \item{items}{A data frame with one row per item: assigned `group`
#'     (`NA` if unassigned), `b_gen`, `b_grp`, `ratio`, and `lambda` (the
#'     first-order loading implied by the higher-order parameterization).}
#'   \item{gamma}{Named vector of second-order loadings (same as in
#'     `groups`).}
#'   \item{unassigned, crossloaded, no_general}{Integer vectors of item
#'     indices flagged by the membership rules above.}
#'   \item{tau, source, call}{Bookkeeping (`source` is `"bifactor"` for
#'     native fits, `"legacy"` for hand-built ones).}
#'
#' @references
#' Schmid, J., & Leiman, J. M. (1957). The development of hierarchical
#' factor solutions. *Psychometrika*, 22(1), 53-61.
#'
#' Yung, Y.-F., Thissen, D., & McLeod, L. D. (1999). On the relationship
#' between the higher-order factor model and the hierarchical factor model.
#' *Psychometrika*, 64(2), 113-128.
#'
#' Zhang, Y., & Chen, J. (2024). Accommodating and extending various models
#' for special effects within the generalized partially confirmatory factor
#' analysis framework. *Applied Psychological Measurement*, 48(4-5),
#' 208-229. \doi{10.1177/01466216241261704}
#'
#' @examples
#' ## higher-order data: the ratios are near-constant and gamma is recovered
#' sim <- sim_fa(N = 500, K = 3, ipf = 6, lam = .8, gamma = c(.8, .7, .6),
#'               rseed = 11)
#' Q <- matrix(-1L, 18, 3)
#' for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
#'                  Q[a, ] <- 0L; Q[a, k] <- 1L }
#' fit <- vbfa(sim$dat, Q, bifactor = TRUE, v0 = .001)
#' special_effects(fit)
#'
#' @seealso [vbfa()], and the `bifactor` vignette for the full workflow.
#' @export
special_effects <- function(object, tau = 0.5) {
  if (!inherits(object, "vbfa"))
    stop("special_effects() needs a vbfa fit.", call. = FALSE)
  if (length(tau) != 1L || !is.finite(tau) || tau < 0 || tau > 1)
    stop("tau must be a single number in [0, 1].", call. = FALSE)
  Q   <- object$Q
  Lam <- object$Lam
  J   <- nrow(Q); P <- ncol(Q)

  ## ---- locate the general and group columns ---------------------------
  if (isTRUE(object$bifactor)) {
    gcol <- 1L
    src  <- "bifactor"
  } else {
    if (!isTRUE(object$orthogonal))
      stop("Not a bifactor fit: the factors are oblique. Fit with ",
           "vbfa(..., bifactor = TRUE).", call. = FALSE)
    full <- which(apply(Q == 1, 2, all))
    if (length(full) != 1L)
      stop("Not recognizably a bifactor fit: expected exactly one ",
           "all-specified column in Q (found ", length(full), "). Fit with ",
           "vbfa(..., bifactor = TRUE) to declare the structure.",
           call. = FALSE)
    gcol <- full
    src  <- "legacy"
  }
  gc <- setdiff(seq_len(P), gcol)
  if (length(gc) < 1L)
    stop("No group factors found alongside the general column.", call. = FALSE)

  ## effective inclusion probabilities (specified -> 1, fixed zero -> 0)
  pi_eff <- object$pi
  pi_eff[Q == 1] <- 1; pi_eff[Q == 0] <- 0

  ## ---- group membership per item ---------------------------------------
  grp <- rep(NA_integer_, J)          # index into gc
  crossloaded <- integer(0)
  for (j in seq_len(J)) {
    act <- which(pi_eff[j, gc] > tau)
    if (length(act) == 1L) grp[j] <- act
    else if (length(act) > 1L) {
      crossloaded <- c(crossloaded, j)
      anc <- act[Q[j, gc[act]] == 1]
      grp[j] <- if (length(anc)) anc[1] else act[which.max(abs(Lam[j, gc[act]]))]
    }
  }
  unassigned <- which(is.na(grp))

  ## ---- ratios, per-group summaries, higher-order parameters ------------
  b_gen <- Lam[, gcol]
  b_grp <- ifelse(is.na(grp), NA_real_, Lam[cbind(seq_len(J), gc[grp])])
  no_general <- which(Q[, gcol] == 0)          # (S-1)-style: no ratio defined
  ratio <- abs(b_grp) / abs(b_gen)
  ratio[no_general] <- NA_real_

  K <- length(gc)
  gnames <- colnames(Q)[gc]
  if (is.null(gnames)) gnames <- paste0("group", seq_len(K))
  n_it <- mr <- sr <- cv <- gam <- D <- rep(NA_real_, K)
  for (k in seq_len(K)) {
    r <- ratio[!is.na(grp) & grp == k & !is.na(ratio)]
    n_it[k] <- length(r)
    if (length(r) >= 1L) mr[k] <- mean(r)
    if (length(r) >= 2L) { sr[k] <- stats::sd(r); cv[k] <- sr[k] / mr[k] }
    if (!is.na(mr[k])) gam[k] <- 1 / sqrt(1 + mr[k]^2)
    ## Eq. 16 (Zhang & Chen, 2024): average squared loading on the group
    ## factor; specified loadings always count, unspecified ones only when
    ## their inclusion probability exceeds tau.
    col <- gc[k]
    inc <- (Q[, col] == 1) | (Q[, col] == -1 & object$pi[, col] > tau)
    if (any(inc)) D[k] <- sum(Lam[inc, col]^2) / sum(inc)
  }
  ok <- !is.na(cv)
  mean_cv <- if (any(ok)) sum(n_it[ok] * cv[ok]) / sum(n_it[ok]) else NA_real_
  lambda <- abs(b_gen) / gam[grp]
  lambda[no_general] <- NA_real_

  out <- list(
    groups = data.frame(group = gnames, n_items = as.integer(n_it),
                        mean_ratio = mr, sd_ratio = sr, cv = cv,
                        gamma = gam, D = D, row.names = NULL),
    mean_cv = mean_cv,
    items = data.frame(item = seq_len(J),
                       group = ifelse(is.na(grp), NA_character_, gnames[grp]),
                       b_gen = b_gen, b_grp = b_grp, ratio = ratio,
                       lambda = lambda, row.names = NULL),
    gamma = stats::setNames(gam, gnames),
    unassigned = unassigned, crossloaded = crossloaded,
    no_general = no_general,
    tau = tau, source = src, call = match.call())
  class(out) <- "special_effects"
  out
}

#' @rdname special_effects
#' @param x A `special_effects` object.
#' @param digits Number of decimals in the printed table.
#' @param ... Unused.
#' @export
print.special_effects <- function(x, digits = 3, ...) {
  cat(sprintf("Special effects of a bifactor fit (%d group factors%s; tau = %g)\n",
              nrow(x$groups),
              if (identical(x$source, "legacy")) ", hand-built design" else "",
              x$tau))
  tab <- x$groups
  tab$mean_ratio <- round(tab$mean_ratio, digits)
  tab$sd_ratio   <- round(tab$sd_ratio, digits)
  tab$cv         <- round(tab$cv, digits)
  tab$gamma      <- round(tab$gamma, digits)
  tab$D          <- round(tab$D, digits)
  print(tab, row.names = FALSE)
  cat(sprintf("Item-weighted mean CV of the group/general loading ratio: %s\n",
              format(round(x$mean_cv, digits))))
  cat("  (near-zero CV = Schmid-Leiman proportional = the fit re-expresses a\n",
      "  higher-order/testlet model, and gamma/lambda are faithful; values\n",
      "  under ~.1 are conventionally read that way -- suggestive, not a test.\n",
      "  With clearly non-constant ratios the higher-order parameters are a\n",
      "  lossy summary of a genuinely richer bifactor.)\n", sep = "")
  if (length(x$no_general))
    cat(sprintf("Note: %d item(s) carry no general loading (general design 0): %s\n",
                length(x$no_general), paste(x$no_general, collapse = ", ")))
  if (length(x$unassigned))
    cat(sprintf("Note: %d item(s) have no active group loading and are excluded from the ratio summaries: %s\n",
                length(x$unassigned), paste(x$unassigned, collapse = ", ")))
  if (length(x$crossloaded))
    cat(sprintf("Note: %d item(s) are active on more than one group factor (assigned to anchor/largest): %s\n",
                length(x$crossloaded), paste(x$crossloaded, collapse = ", ")))
  invisible(x)
}
