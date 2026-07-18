## S3 methods for fitted vbpm models. The class was introduced at 0.3.0; see
## the vbpm_fit help page for the full rationale and what the class does and
## does not change for existing code.

#' The vbpm_fit class: what it is and why it exists
#'
#' Every estimator in `vbpm` returns its results with an S3 class attached:
#' [vbfa()] returns `c("vbfa", "vbpm_fit")` and [vbmimic()] returns
#' `c("vbmimic", "vbpm_fit")`. This page explains what that means in
#' practice, because "S3 class" sounds like more machinery than it is.
#'
#' @section What an S3 class actually is:
#' An S3 class is nothing more than a **label attached to an ordinary R
#' object**. A `vbpm_fit` is still a plain named list — the class attribute
#' changes none of its contents and none of the ways you already use it:
#'
#' * `fit$Lam`, `fit$eta`, `fit$ELBO`, ... work exactly as before;
#' * `names(fit)` lists every component;
#' * `str(fit)` shows the full structure;
#' * `unclass(fit)` strips the label and gives back the bare list, should you
#'   ever want it.
#'
#' The label does exactly two things. First, generic functions such as
#' [print()] can now **dispatch**: typing `fit` at the console no longer
#' floods the screen with every matrix in the object, but shows a compact
#' summary (model type, dimensions, convergence, ELBO, active loadings).
#' Second, functions receiving a fit can **recognize what it is**: this is
#' how [vb_fit()] knows a bifactor fit from an oblique one without being told
#' (the fit carries its own `orthogonal` and `ld` settings), and how it
#' refuses a [vbmimic()] fit rather than computing meaningless statistics
#' on it.
#'
#' @section Why the family shares a parent class:
#' The specific class (`"vbfa"`, `"vbmimic"`) comes first, the family class
#' `"vbpm_fit"` second. Methods written for `vbpm_fit` therefore apply to
#' every model in the family unless a model overrides them — one `print`
#' method serves all estimators, and future family members (e.g. further
#' variational models) inherit sensible behaviour on day one.
#'
#' @section What is deliberately NOT hidden:
#' Some R packages wrap results in opaque objects whose internals are
#' discouraged territory. `vbpm` does the opposite: the list components are
#' documented in each estimator's help page (`?vbfa`, `?vbmimic`) and are
#' part of the public API. The class adds convenience on top; it takes
#' nothing away.
#'
#' @name vbpm_fit
#' @aliases vbpm_fit-class
#' @seealso [vbfa()], [vbmimic()], [print.vbpm_fit()]
NULL


#' Compact display of a fitted vbpm model
#'
#' Prints a short summary — model type, dimensions, convergence, timing, and
#' selection results — instead of dumping the full list. Access the parts as
#' usual (`fit$Lam` etc.; see [vbpm_fit]); `str(fit)` shows everything.
#'
#' @param x A fit returned by [vbfa()] or [vbmimic()].
#' @param tau Threshold on the posterior inclusion probability used only for
#'   the "active loadings" line (default `0.5`, matching [vb_fit()]).
#' @param ... Ignored (present for compatibility with the [print()] generic).
#'
#' @return `x`, invisibly (the R convention for print methods, so `print(fit)`
#'   can be used in pipes without losing the object).
#'
#' @examples
#' sim <- sim_fa(N = 300, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' Q <- matrix(-1L, ncol(sim$dat), 3)
#' for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
#'                  Q[a, ] <- 0L; Q[a, k] <- 1L }
#' fit <- vbfa(sim$dat, Q)
#' fit               # dispatches here: compact summary, not a list dump
#' fit$Lam[1:3, ]    # the object is still an ordinary list underneath
#'
#' @export
print.vbpm_fit <- function(x, tau = 0.5, ...) {
  is_mimic <- inherits(x, "vbmimic")

  if (is_mimic) {
    J <- nrow(x$A); K <- ncol(x$A); P <- ncol(x$B); N <- nrow(x$eta)
    cat("vbmimic: regularized VB MIMIC model\n")
    cat(sprintf("  %d items, %d factors, %d covariates, N = %d\n",
                J, K, P, N))
  } else {
    J <- nrow(x$Lam); K <- ncol(x$Lam); N <- nrow(x$eta)
    cat("vbfa: VB partially confirmatory factor analysis\n")
    cat(sprintf("  %d items, %d factors, N = %d  (%s%s)\n", J, K, N,
                if (isTRUE(x$orthogonal)) "orthogonal/bifactor" else "oblique",
                if (isTRUE(x$ld)) ", local dependence" else ""))
  }

  n_stage <- x$path$n_stage
  cat(sprintf("  %s in %d iterations%s (%.2f secs)\n",
              if (isTRUE(x$flag == 1)) "converged" else "NOT converged",
              x$iter,
              if (!is.null(n_stage) && n_stage > 1)
                sprintf(" over %d v0 stages", n_stage) else "",
              as.numeric(x$time, units = "secs")))

  if (!is.null(x$ELBO))
    cat("  ELBO:", if (is.na(x$ELBO)) "NA (not available for this model)"
        else format(x$ELBO, digits = 8), "\n")

  if (is_mimic) {
    unsA <- x$Q_A == -1
    unsB <- x$Q_B == -1
    if (any(unsA))
      cat(sprintf("  active unspecified loadings    (PIP >= %.2g): %d of %d\n",
                  tau, sum(x$pi_A[unsA] >= tau), sum(unsA)))
    if (any(unsB))
      cat(sprintf("  active structural coefficients (PIP >= %.2g): %d of %d\n",
                  tau, sum(x$pi_B[unsB] >= tau), sum(unsB)))
  } else {
    uns <- x$Q == -1
    if (any(uns))
      cat(sprintf("  active unspecified loadings (PIP >= %.2g): %d of %d\n",
                  tau, sum(x$pi[uns] >= tau), sum(uns)))
  }

  cat("\nComponents: ", paste(names(x), collapse = ", "), "\n", sep = "")
  cat("Access them with $ as usual; see ?vbpm_fit.\n")
  invisible(x)
}


#' Extract coefficients from a fitted vbpm model
#'
#' For a [vbfa()] fit, returns the loading matrix. For a [vbmimic()] fit,
#' returns a list with the measurement loadings `A` and the structural
#' coefficients `B`.
#'
#' @param object A fit returned by [vbfa()] or [vbmimic()].
#' @param ... Ignored.
#' @return A matrix (vbfa) or a list of two matrices (vbmimic).
#' @examples
#' sim <- sim_fa(N = 300, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' Q <- matrix(-1L, ncol(sim$dat), 3)
#' for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
#'                  Q[a, ] <- 0L; Q[a, k] <- 1L }
#' fit <- vbfa(sim$dat, Q)
#' round(coef(fit), 2)
#' @export
coef.vbpm_fit <- function(object, ...) {
  if (inherits(object, "vbmimic")) list(A = object$A, B = object$B)
  else object$Lam
}
