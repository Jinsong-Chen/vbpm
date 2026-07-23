## SEM-like fit statistics for vbpm models. Hard selection (posterior
## inclusion probability >= tau) is the default for every statistic; soft
## variants carry an "_S" suffix.
##
## fit_stats() is the S3 generic; vb_fit() at the foot of this file is the
## pre-0.4.0 name, kept as a thin compatibility wrapper. See the roxygen
## documentation of fit_stats() below.

#' Fit statistics for vbpm models
#'
#' Computes structural-equation-model fit statistics for a fitted [vbfa()]
#' object under hard selection (an unspecified loading counts as active when its
#' posterior inclusion probability is at least `tau`). Absolute indices
#' (RMSEA, SRMR, CFI, TLI) and information criteria (AIC, BIC) are returned,
#' together with the model's ELBO.
#'
#' @param object,fit A fitted `vbpm_fit` object. `fit` is the historical
#'   argument name retained by [vb_fit()].
#' @param Y The `N x J` data matrix the model was fit to.
#' @param Q The `J x K` loading design matrix used in the fit.
#' @param tau Hard-selection threshold on the posterior inclusion probability
#'   (default `0.5`).
#' @param gamma Reserved (extended BIC tuning); currently unused in the
#'   returned vector.
#' @param orthogonal Logical; whether the fit is an orthogonal (bifactor)
#'   model, in which case factor correlations are not free parameters. The
#'   default `NULL` reads the setting from the fit object itself, which is the
#'   recommended use. Supplying a value that contradicts the fit raises an
#'   error rather than silently miscounting parameters.
#' @param rank_adjust Logical; if `TRUE`, the effective parameter count is the
#'   rank of the Jacobian of the model-implied covariance (needs the `numDeriv`
#'   package). Automatically skipped for large `J`.
#' @param rank_max_J Skip the Jacobian rank check when `J` exceeds this value
#'   (default `100`); the nominal parameter count is used instead.
#' @param ... Further arguments passed to methods. Unused by the methods
#'   supplied with the package.
#'
#' @return An object of class `vbpm_fit_stats` (a named numeric vector) with
#'   `t_nom`, `t`, `t_S` (parameter counts), `RMSEA`, `SRMR`, `CFI`, `TLI`,
#'   `AIC`, `BIC` (hard selection), `AIC_S`, `BIC_S` (soft selection), `ELBO`,
#'   and `objective`. The last is the model's comparability score: equal to
#'   `ELBO` for a diagonal fit, and the terminal VECM objective under local
#'   dependence, where `ELBO` is `NA`. It is comparable only across fits
#'   sharing the `objective_type` attribute, which is also attached.
#'
#'   The `vbmimic` method returns the same shape with `n_active_coef` in place
#'   of a full soft count and typed `NA` elsewhere; see [vbmimic()].
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' Jin, Y., & Chen, J. (2025). Regularized variational approximation for
#' partially confirmatory factor analysis. *Structural Equation Modeling: A
#' Multidisciplinary Journal*, 32(3), 437-449.
#' \doi{10.1080/10705511.2024.2432612}
#'
#' @examples
#' sim <- sim_fa(N = 300, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' Q <- matrix(-1L, ncol(sim$dat), 3)
#' for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]
#'                  Q[a, ] <- 0L; Q[a, k] <- 1L }
#' fit <- vbfa(sim$dat, Q)
#' round(fit_stats(fit), 3)
#'
#' \donttest{
#' ## bifactor: orthogonal is read from the fit object, so the same call works
#' Qb <- cbind(1L, Q)
#' fb <- vbfa(sim$dat, Qb, orthogonal = TRUE)
#' round(fit_stats(fb), 3)
#' }
#'
#' @seealso [vbfa()], [pefa()]
#' @export
fit_stats <- function(object, ...) UseMethod("fit_stats")

#' @rdname fit_stats
#' @export
## `...` is required for S3 consistency with the generic fit_stats(object, ...)
fit_stats.vbfa <- function(object, Y = NULL, Q = object$Q, tau = 0.50, gamma = 0.5,
                     orthogonal = NULL, rank_adjust = TRUE, rank_max_J = 100,
                     ...) {
    fit <- object

    ## resolve `orthogonal` from the fit object; this used to default to FALSE,
    ## which silently miscounted K(K-1)/2 factor correlations for bifactor fits
    if (is.null(orthogonal)) {
        orthogonal <- isTRUE(fit$orthogonal)
    } else if (!is.null(fit$orthogonal) &&
               !identical(isTRUE(orthogonal), isTRUE(fit$orthogonal))) {
        stop("orthogonal = ", orthogonal, " contradicts the fit, which was ",
             "estimated with orthogonal = ", isTRUE(fit$orthogonal),
             ". Omit the argument to use the fit's own setting.",
             call. = FALSE)
    }

    N <- if (is.null(Y)) fit$nobs else nrow(Y)
    J <- nrow(fit$Lam); K <- ncol(Q); Ns <- N - 1
    m  <- J * (J + 1) / 2                    # unique covariance moments

    ## ---- plug-in variational posterior summaries ----------------------
    Lam_h <- fit$Lam                    # Lambda_hat
    is_ld <- isTRUE(fit$ld)
    Wres  <- if (is_ld) fit$W else diag(1 / fit$PsiInv, J)
    if (is.null(Wres)) stop("Residual covariance is unavailable.", call. = FALSE)
    Phi   <- fit$Phi

    pi_ <- fit$pi;  pi_[Q == 1] <- 1;  pi_[Q == 0] <- 0   # pi_hat_{jk}
    aH  <- (Q == 1) | (Q == -1 & pi_ >= tau)             # hard-selection set
    Lam_H <- Lam_h * aH                                  # hard-selected
    Lam_S <- Lam_h * pi_                                 # soft-selected

    ## ---- model-implied covariances:  Sigma = Lam Phi Lam' + Psi -------
    sigma <- function(L) L %*% Phi %*% t(L) + Wres
    SgH <- sigma(Lam_H); SgS <- sigma(Lam_S)
    S <- if (!is.null(fit$sample_cov)) fit$sample_cov else stats::cor(Y)
    Sg0 <- diag(diag(S))

    ## ---- discrepancy F + plug-in log-likelihood -----------------------
    ## stable log-determinant (det() underflows to 0 for large, collinear S,
    ## which made log(det(S)) = -Inf and RMSEA/CFI degenerate at large J).
    logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
    F_ML <- function(Sg) {                   # log|Sg| + tr(S Sg^-1) - log|S| - J
        logdet(Sg) + sum(diag(S %*% solve(Sg))) - logdet(S) - J
    }
    ll <- function(Sg) -0.5 * Ns * (J * log(2 * pi) + logdet(Sg) +
                                    sum(diag(S %*% solve(Sg))))
    T_H <- Ns * F_ML(SgH);  T_0 <- Ns * F_ML(Sg0)

    ## ---- parameter counts and degrees of freedom ----------------------
    p_H   <- sum(aH)                          # required + active-unspecified
    n_phi <- if (orthogonal) 0L else K * (K - 1) / 2  # free factor correlations
    edge_H <- matrix(FALSE, J, J)
    edge_S <- 0
    if (is_ld) {
        qe_upper <- upper.tri(fit$Qe)
        edge_H <- qe_upper & (fit$Qe == 1 | (fit$Qe == -1 & fit$q_star >= tau))
        edge_S <- sum(fit$Qe[qe_upper] == 1) +
                  sum(fit$q_star[qe_upper & fit$Qe == -1])
    }
    n_edge <- sum(edge_H)
    t_nom <- p_H + J + n_phi + n_edge

    ## rank-adjusted count: t = rank( d vech(Sigma) / d theta' )
    resid0 <- if (is_ld) c(diag(fit$Psi), fit$Psi[edge_H]) else diag(Wres)
    th0 <- c(Lam_H[aH], if (n_phi > 0) Phi[lower.tri(Phi)], resid0)
    rebuild_vech <- function(th) {
        L <- matrix(0, J, K); L[aH] <- th[seq_len(p_H)]
        P <- diag(K)
        if (n_phi > 0) {
            P[lower.tri(P)] <- th[p_H + seq_len(n_phi)]
            P[upper.tri(P)] <- t(P)[upper.tri(P)]
        }
        pos <- p_H + n_phi
        rd <- th[pos + seq_len(J)]
        if (is_ld) {
            Prec <- diag(rd, J)
            if (n_edge > 0) {
                ev <- th[pos + J + seq_len(n_edge)]
                Prec[edge_H] <- ev
                Prec[t(edge_H)] <- ev
            }
            Rcov <- solve(Prec)
        } else Rcov <- diag(rd)
        Sg <- L %*% P %*% t(L) + Rcov
        Sg[lower.tri(Sg, diag = TRUE)]
    }
    ## Rank-adjusted t needs a numDeriv Jacobian of vech(Sigma) (J(J+1)/2 rows),
    ## whose cost explodes with J; auto-skip for large J and fall back to t_nom.
    do_rank <- isTRUE(rank_adjust) && (J <= rank_max_J) &&
        requireNamespace("numDeriv", quietly = TRUE)
    if (do_rank) {
        Delta  <- numDeriv::jacobian(rebuild_vech, th0)
        sv     <- svd(Delta, nu = 0, nv = 0)$d
        t_H    <- sum(sv > 1e-6 * max(sv))    # rank count; == t_nom when identified
    } else {
        t_H    <- t_nom                       # nominal count (rank check skipped)
    }
    t_S <- J + n_phi + sum(Q == 1) + sum(pi_[Q == -1]) + edge_S

    df_H <- m - t_H
    df_0 <- J * (J - 1) / 2

    ## ---- absolute / incremental fit (hard-selected) -------------------
    RMSEA <- if (df_H > 0)
                 sqrt(max(0, (T_H - df_H) / (df_H * Ns)))
             else NA_real_
    si <- 1 / sqrt(diag(S))
    R0 <- si * S   * rep(si, each = J)
    Rh <- si * SgH * rep(si, each = J)
    SRMR <- sqrt(2 * sum((R0 - Rh)[lower.tri(R0, diag = TRUE)]^2) /
                 (J * (J + 1)))
    denom <- max(T_H - df_H, T_0 - df_0, 0)
    CFI <- if (denom > 0) 1 - max(T_H - df_H, 0) / denom else NA_real_
    TLI <- if (df_H > 0) (T_0 / df_0 - T_H / df_H) / (T_0 / df_0 - 1) else NA_real_

    ## ---- relative fit: information criteria ---------------------------
    llH <- ll(SgH); llS <- ll(SgS)
    AIC   <- -2 * llH + 2 * t_H
    BIC   <- -2 * llH + log(Ns) * t_H
    AIC_S <- -2 * llS + 2 * t_S
    BIC_S <- -2 * llS + log(Ns) * t_S

    out <- c(t_nom = t_nom, t = t_H, t_S = t_S,
      RMSEA = RMSEA, SRMR = SRMR,
      CFI = CFI, TLI = TLI,
      AIC = AIC, BIC = BIC,
      AIC_S = AIC_S, BIC_S = BIC_S,
      ELBO = { e <- as.numeric(fit$ELBO); if (length(e) && is.finite(e)) e else NA_real_ },
      ## The model's own comparability score. For a diagonal fit this equals
      ## ELBO; under LD it is the terminal VECM objective, where ELBO is NA.
      ## Comparable only across fits sharing the objective_type attribute.
      objective = { o <- as.numeric(fit$objective)
                    if (length(o) && is.finite(o)) o else NA_real_ })
    class(out) <- c("vbpm_fit_stats", "numeric")
    attr(out, "model") <- "vbfa"
    attr(out, "objective_type") <- if (!is.null(fit$objective_type))
        fit$objective_type else if (is_ld) "vecm" else "elbo"
    out
}

#' @rdname fit_stats
#' @export
fit_stats.vbmimic <- function(object, Y = NULL, Q = NULL, tau = 0.5, ...) {
    ## Soft-selected COEFFICIENTS only: the measurement loadings plus the
    ## structural coefficients. This is deliberately not reported as `t_S`,
    ## which elsewhere means a full soft parameter count -- it omits the
    ## residual precisions, the factor covariance, and the inclusion rates.
    ## It must not be fed to an information criterion.
    nA <- sum(object$Q_A == 1) + sum(object$pi_A[object$Q_A == -1])
    nB <- sum(object$Q_B == 1) + sum(object$pi_B[object$Q_B == -1])
    out <- c(t_nom = NA_real_, t = NA_real_, t_S = NA_real_,
             n_active_coef = nA + nB,
             RMSEA = NA_real_, SRMR = NA_real_, CFI = NA_real_, TLI = NA_real_,
             AIC = NA_real_, BIC = NA_real_, AIC_S = NA_real_, BIC_S = NA_real_,
             ELBO = NA_real_, objective = NA_real_)
    class(out) <- c("vbpm_fit_stats", "numeric")
    attr(out, "model") <- "vbmimic"
    attr(out, "objective_type") <- "none"
    attr(out, "note") <- paste(
        "SEM-like fit statistics require the MIMIC objective and joint",
        "covariance derivation, neither of which is available yet.",
        "n_active_coef counts soft-selected measurement and structural",
        "coefficients only and is not a full parameter count.")
    out
}

#' @rdname fit_stats
#' @export
fit_stats.vbpm_fit <- function(object, ...)
    stop("No fit_stats() method is defined for model class: ", class(object)[1],
         call. = FALSE)

# Historical API retained for v0.4.0 migration.
#' @rdname fit_stats
#' @export
vb_fit <- function(fit, Y, Q = fit$Q, ...) fit_stats(fit, Y = Y, Q = Q, ...)
