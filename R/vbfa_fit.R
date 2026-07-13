###############################################################################
##  vbfa_fit.R  --  SEM-like fit statistics for PCFA-VB / PEFA  (vbfa_r2 stack)
##
##  Updated companion of analysis/vbfa_fit.R for the vbfa_r2 estimator
##  (2026-07-13; the analysis/ original stays untouched for the frozen P1
##  scripts). Notation follows Jin & Chen (2025), Regularized Variational
##  Approximation for Partially Confirmatory Factor Analysis, SEM 32(3).
##
##  CHANGES vs analysis/vbfa_fit.R:
##   * HARD selection is the default for every statistic — no _H suffix:
##       t_nom t t_S RMSEA SRMR CFI TLI AIC BIC AIC_S BIC_S ELBO
##     (_S marks the soft-selection variants where one exists; everything
##     unsuffixed is hard-selected at threshold tau.)
##   * LD guard: an LD fit from vbfa_r2 (Qe supplied) has PsiInv = NULL and
##     is REJECTED with a clear error. r2.1 upgrade: use the full residual
##     covariance W (already returned by QUIC) in the implied covariance
##     Sigma = Lam Phi Lam' + W, and add the active residual off-diagonals
##     (q_star >= tau) to the parameter counts.
##
##  Usage:   fit <- vbfa(Y, Q);  idx <- vbfa_fit(fit, Y, Q)
###############################################################################

#' SEM-like fit statistics for a vbfa model
#'
#' Computes structural-equation-model fit statistics for a fitted [vbfa()]
#' object under hard selection (an unspecified loading counts as active when its
#' posterior inclusion probability is at least `tau`). Absolute indices
#' (RMSEA, SRMR, CFI, TLI) and information criteria (AIC, BIC) are returned,
#' together with the model's ELBO.
#'
#' @param fit A fitted object returned by [vbfa()] (diagonal-residual model;
#'   local-dependence fits are not yet supported and raise an error).
#' @param Y The `N x J` data matrix the model was fit to.
#' @param Q The `J x K` loading design matrix used in the fit.
#' @param tau Hard-selection threshold on the posterior inclusion probability
#'   (default `0.5`).
#' @param gamma Reserved (extended BIC tuning); currently unused in the
#'   returned vector.
#' @param orthogonal Logical; must match the value used in [vbfa()]. When
#'   `TRUE`, factor correlations are not free parameters.
#' @param rank_adjust Logical; if `TRUE`, the effective parameter count is the
#'   rank of the Jacobian of the model-implied covariance (needs the `numDeriv`
#'   package). Automatically skipped for large `J`.
#' @param rank_max_J Skip the Jacobian rank check when `J` exceeds this value
#'   (default `100`); the nominal parameter count is used instead.
#'
#' @return A named numeric vector: `t_nom`, `t`, `t_S` (parameter counts),
#'   `RMSEA`, `SRMR`, `CFI`, `TLI`, `AIC`, `BIC` (hard selection), `AIC_S`,
#'   `BIC_S` (soft selection), and `ELBO`.
#'
#' @references
#' Chen, J., & Jin, Y. (2026). Recovering latent structures after variational
#' Bayesian variable selection: Fit assessment and factor-number selection in
#' partially exploratory factor analysis. *arXiv preprint* arXiv:2607.07159.
#'
#' Jin, Y., & Chen, J. (2025). Regularized variational approximation for
#' partially confirmatory factor analysis. *Structural Equation Modeling: A
#' Multidisciplinary Journal*, 32(3), 437–449.
#' \doi{10.1080/10705511.2024.2432612}
#'
#' @seealso [vbfa()], [pefa_vb()]
#' @export
vbfa_fit <- function(fit, Y, Q, tau = 0.50, gamma = 0.5,
                     orthogonal = FALSE, rank_adjust = TRUE, rank_max_J = 100) {
    ## LD fits (vbfa_r2 with Qe) are not supported yet — see header (r2.1).
    if (is.null(fit$PsiInv))
        stop("vbfa_fit does not support LD fits yet (fit$PsiInv is NULL). ",
             "r2.1 will use W in the implied covariance and count residual ",
             "edges; until then compute fit statistics on a non-LD fit.",
             call. = FALSE)

    N  <- nrow(Y); J <- ncol(Y); K <- ncol(Q); Ns <- N - 1
    m  <- J * (J + 1) / 2                    # unique covariance moments

    ## ---- plug-in variational posterior summaries ----------------------
    Lam_h <- fit$Lam                    # Lambda_hat
    Psi   <- 1 / fit$PsiInv             # error variances (diag Psi_hat)
    Phi   <- fit$Phi

    pi_ <- fit$pi;  pi_[Q == 1] <- 1;  pi_[Q == 0] <- 0   # pi_hat_{jk}
    aH  <- (Q == 1) | (Q == -1 & pi_ >= tau)             # hard-selection set
    Lam_H <- Lam_h * aH                                  # hard-selected
    Lam_S <- Lam_h * pi_                                 # soft-selected

    ## ---- model-implied covariances:  Sigma = Lam Phi Lam' + Psi -------
    sigma <- function(L) L %*% Phi %*% t(L) + diag(Psi)
    SgH   <- sigma(Lam_H);  SgS <- sigma(Lam_S)
    Sg0   <- diag(diag(cor(Y)))              # independence baseline diag(S)
    S     <- cor(Y)                          # vbfa fits scale(Y), so S = R

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
    t_nom <- p_H + J + n_phi                  # nominal count (J = error variances)

    ## rank-adjusted count: t = rank( d vech(Sigma) / d theta' )
    th0 <- if (orthogonal) c(Lam_H[aH], Psi) else c(Lam_H[aH], Phi[lower.tri(Phi)], Psi)
    rebuild_vech <- function(th) {
        L <- matrix(0, J, K); L[aH] <- th[seq_len(p_H)]
        P <- diag(K)
        if (n_phi > 0) {
            P[lower.tri(P)] <- th[p_H + seq_len(n_phi)]
            P[upper.tri(P)] <- t(P)[upper.tri(P)]
        }
        Pp <- th[p_H + n_phi + seq_len(J)]
        Sg <- L %*% P %*% t(L) + diag(Pp)
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
    t_S    <- J + n_phi + sum(Q == 1) + sum(pi_[Q == -1])  # soft count

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

    c(t_nom = t_nom, t = t_H, t_S = t_S,
      RMSEA = RMSEA, SRMR = SRMR,
      CFI = CFI, TLI = TLI,
      AIC = AIC, BIC = BIC,
      AIC_S = AIC_S, BIC_S = BIC_S,
      ELBO = { e <- as.numeric(fit$ELBO); if (is.finite(e)) e else NA_real_ })
}
