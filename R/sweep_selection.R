## Factor-count calculators for pefa() sweeps.
##
## 0.9.0 keeps exactly one, and pefa() never invokes it: the package returns
## criterion paths and the analyst applies a rule. The 0.8.x internals that
## selected inside pefa() (.pefa_cuts, .pefa_sustain, .scan_gain,
## .pefa_boundary, .sweep_gain_selections) went with $selected_K.

#' Sustained-drop elbow rule for a factor count
#'
#' Applies the sustained-drop elbow rule to one criterion path: a candidate
#' qualifies when the marginal gains that follow it have all fallen below a
#' percentage of the largest gain on the same path. The rule is invariant to
#' positive affine transformations of that path, so the candidates must
#' otherwise share the same data, priors, tuning, and optimization setup for
#' the comparison to carry meaning.
#'
#' `score` is oriented larger-is-better, so a single call serves an ELBO path,
#' `-AIC`, or `-BIC`. A candidate qualifies only when **all** `sustain` gains
#' following it are **strictly** below the threshold. The look-ahead is
#' complete: a candidate with fewer than `sustain` gains after it cannot
#' qualify, so a short flat suffix at the right edge of the window never
#' selects. When nothing qualifies the answer is `NA_integer_`. There is no
#' boundary fallback, and the returned count is therefore never a consolation
#' value read off either end of the window.
#'
#' [pefa()] does not call this function. It returns the criterion paths and
#' leaves the rule, and the responsibility for the rule, with the analyst; the
#' defaults below are conveniences of a criterion-generic calculator, not a
#' claim that one setting is universally preferable.
#'
#' @param K Candidate factor counts: a finite, unique, integer-valued numeric
#'   vector that is consecutive once sorted. Any input order is accepted, since
#'   `K` and `score` are sorted by `K` before gains are computed.
#' @param score The criterion at each candidate in `K`, oriented so that larger
#'   is better; pass the negative of an information criterion. Must be finite,
#'   the same length as `K`, and on a scale whose adjacent differences remain
#'   representable.
#' @param delta Threshold as a percentage of the largest marginal gain on the
#'   path: one finite number in `[0, 100]`.
#' @param sustain Number of consecutive sub-threshold gains required: one
#'   positive whole number.
#'
#' @return The smallest qualifying count, as a length-one integer, or
#'   `NA_integer_` when no count qualifies. There are three `NA_integer_`
#'   cases: a one-candidate path, which holds no gain for any rule to read; a
#'   path whose largest gain is not positive (`Delta_max <= 0`), where a
#'   percentage of it is not a usable threshold; and a path on which no start
#'   index is followed by a complete run of `sustain` sub-threshold gains.
#'   Invalid input is an error, never `NA`.
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' @examples
#' ## An ELBO path that flattens: gains 300, 200, 5, 3 against a 20% cut of 300.
#' elbo <- c(-5400, -5100, -4900, -4895, -4892)
#' select_K_elbow(2:6, elbo, delta = 20, sustain = 1)
#'
#' ## The same call reads a BIC path once it is negated to larger-is-better.
#' bic <- c(10900, 10500, 10310, 10395, 10480)
#' select_K_elbow(2:6, -bic, delta = 20, sustain = 1)
#'
#' ## Requiring the drop to persist changes the reading on a jagged path.
#' jagged <- c(-5400, -5100, -5095, -4900, -4897, -4895)
#' select_K_elbow(2:7, jagged, delta = 20, sustain = 1)
#' select_K_elbow(2:7, jagged, delta = 20, sustain = 2)
#'
#' ## No positive gain anywhere: NA_integer_, not the left edge.
#' select_K_elbow(2:4, c(10, 9, 8), delta = 20, sustain = 1)
#'
#' ## Complete look-ahead: K = 5 would need a gain from K = 6 to K = 7, which
#' ## this window does not contain, so the path yields NA_integer_.
#' select_K_elbow(4:6, c(0, 100, 105), delta = 20, sustain = 2)
#'
#' ## A one-candidate window contains no gain at all.
#' select_K_elbow(3L, -4900, delta = 20, sustain = 1)
#'
#' @seealso [pefa()]
#' @export
select_K_elbow <- function(K, score, delta = 20, sustain = 1) {
  ## A matrix would silently change what diff() means, so paths must be plain
  ## numeric vectors.
  path_vector <- function(x) is.numeric(x) && is.null(dim(x))
  if (!path_vector(K) || !path_vector(score))
    stop("K and score must be numeric vectors.", call. = FALSE)
  if (length(K) != length(score)) {
    stop("K and score must have the same length (got ", length(K), " and ",
         length(score), ").", call. = FALSE)
  }
  if (!length(K)) stop("K must not be empty.", call. = FALSE)
  if (!all(is.finite(K))) stop("K must be finite.", call. = FALSE)
  if (any(K != round(K)) || any(abs(K) > .Machine$integer.max)) {
    stop("K must contain whole numbers representable as integers.",
         call. = FALSE)
  }
  if (anyDuplicated(K)) stop("K must not contain duplicates.", call. = FALSE)
  if (!all(is.finite(score))) stop("score must be finite.", call. = FALSE)
  if (length(delta) != 1L || !is.numeric(delta) || !is.finite(delta) ||
      delta < 0 || delta > 100) {
    stop("delta must be a single number in [0, 100].", call. = FALSE)
  }
  if (length(sustain) != 1L || !is.numeric(sustain) || !is.finite(sustain) ||
      sustain < 1 || sustain != round(sustain)) {
    stop("sustain must be a single positive whole number.", call. = FALSE)
  }
  sustain <- as.integer(sustain)

  ## Sort first: the rule is defined on the ascending path, and a gap in it is
  ## rejected rather than bridged, because a bridged gain is not a gain.
  o <- order(K)
  K <- as.integer(K[o])
  ## Double storage keeps an integer criterion path from overflowing in diff().
  score <- as.double(score[o])
  if (any(diff(K) != 1L))
    stop("K must contain consecutive whole numbers.", call. = FALSE)

  ## A one-candidate window holds no gain, so no rule can fire on it.
  if (length(K) == 1L) return(NA_integer_)

  gain <- diff(score)
  if (!all(is.finite(gain))) {
    ## Finite scores can still overflow the double range when differenced.
    stop("Adjacent score differences must be finite.", call. = FALSE)
  }

  gain_max <- max(gain)
  ## A percentage of a nonpositive largest gain is not a usable threshold.
  if (gain_max <= 0) return(NA_integer_)
  threshold <- delta / 100 * gain_max

  ## Complete look-ahead: the last admissible start is the one that still has
  ## `sustain` gains after it, so a short suffix at the right edge cannot
  ## qualify however flat it is.
  last_start <- length(gain) - sustain + 1L
  if (last_start < 1L) return(NA_integer_)
  for (i in seq_len(last_start)) {
    ## Strict inequality is part of the published rule; the first hit is the
    ## smallest qualifying count.
    if (all(gain[i:(i + sustain - 1L)] < threshold)) return(K[i])
  }
  NA_integer_
}
