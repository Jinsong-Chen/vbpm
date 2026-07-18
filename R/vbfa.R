## Internal: log-normalizing constant of the inverse generalized Wishart,
## used in the oblique ELBO.
logC_igw <- function(xi, Lam, P) {
  0.5 * (xi - P + 1) * log(det(Lam)) -
    0.5 * P * (xi - P + 1) * log(2) -
    0.25 * P * (P - 1) * log(pi) -
    sum(lgamma((xi - P - (1:P)) / 2 + 1))
}


#' Variational Bayes partially confirmatory factor analysis
#'
#' Fits a partially confirmatory factor-analytic measurement model by
#' regularized mean-field variational Bayes. Loadings on *specified* factors
#' are anchored by the design matrix `Q`; loadings on *unspecified* factors are
#' selected from the data through continuous spike-and-slab priors. Optional
#' features are a warm-started dynamic regularization path, an orthogonal
#' bifactor parameterization, and sparse residual (local-dependence)
#' estimation.
#'
#' @param Y Numeric `N x J` data matrix (observations by items). Standardized
#'   internally. Missing values are not supported and raise an error.
#' @param Q Integer `J x K` design matrix for the loadings: `1` = specified
#'   (anchored) loading, `0` = fixed zero, `-1` = unspecified (estimated by
#'   spike-and-slab). `K` is the number of factors.
#' @param ld Logical switch for local-dependence (sparse residual) estimation.
#'   `FALSE` (default) fits diagonal residuals and ignores `Qe` (with a warning
#'   if one is supplied). `TRUE` estimates a sparse residual precision by a
#'   graphical spike-and-slab prior solved by QUIC, with the search space given
#'   by `Qe`.
#' @param Qe Residual design matrix (`J x J`), read only when `ld = TRUE`.
#'   Entries `1` (freely estimated residual dependence), `-1` (uncertain;
#'   selected by spike-and-slab), `0` (shrunk to zero). Symmetric; its diagonal
#'   is ignored. The default `NULL` means fully exploratory local dependence,
#'   i.e. `Qe = matrix(-1, J, J)`.
#' @param orthogonal Logical. If `TRUE`, the factor correlation is fixed at the
#'   identity (an orthogonal bifactor model: one general plus orthogonal group
#'   factors). Default `FALSE` (oblique; factor correlation estimated).
#' @param v0 Loading-spike variance. A decreasing **vector** gives a
#'   warm-started dynamic regularization path (each stage runs to convergence
#'   and warm-starts the next); a **scalar** gives a single fixed spike. The
#'   default path ends at `0.001`; a scalar `v0 = 0.001` reproduces the
#'   fixed-spike estimator exactly.
#'
#'   For a **gentler path that aids convergence**, pass a finer, longer,
#'   log-spaced descent starting from a larger `v0`, e.g.
#'   `v0 = exp(seq(log(0.05), log(0.001), length.out = 10))`. Each stage is
#'   then a small perturbation of the previous one, so its warm start is a good
#'   initialization and it converges quickly and stably. Prefer this when a fit
#'   fails to converge within `max_it`, when the loading solution looks
#'   unstable, or for near-singular / strongly regularized models: it navigates
#'   the non-convex spike-and-slab landscape better than a single hard fit, and
#'   is usually preferable to (and often faster than) simply raising `max_it`.
#'   Keep the endpoint at your target; log spacing suits the multiplicative
#'   spike.
#' @param max_it Maximum variational iterations **per stage**. Default `5000`.
#' @param convChk Logical; print per-iteration relative error. Default `TRUE`.
#' @param tolVal Convergence tolerance on the maximum relative error.
#' @param ld_control A named list of local-dependence controls, read only when
#'   `ld = TRUE`: `xi0` (spike-penalty path, units of `N`; default
#'   `seq(0.1, 1, length.out = 5)`), `xi1` (slab penalty, units of `N`;
#'   default `0.01`), `diag_penalty` (`1` penalizes the precision diagonal with
#'   `xi1`, `0` leaves it unpenalized), `quic_eps`, `quic_max_it`, and the Beta
#'   prior parameters `a1`, `b1` on the local-dependence proportion.
#'
#' @return A list with elements including `Lam` (posterior mean loadings),
#'   `eta` (word/observation factor scores), `Phi` (factor correlation; the
#'   identity when `orthogonal = TRUE`), `pi` (posterior inclusion
#'   probabilities), `PsiInv` (residual precisions; `NULL` under local
#'   dependence), `ELBO` (evidence lower bound; `NA` under local dependence),
#'   `iter`, `flag` (`1` if converged), `orthogonal` and `ld` (the settings the
#'   model was fit with), and — when `ld = TRUE` — the residual precision
#'   `Psi`, its inverse `W`, the residual-edge activeness `q_star`, and path
#'   metadata.
#'
#' @section Notes:
#' The diagonal-residual estimator is **deterministic** (fixed initialization,
#' no random numbers), so results are reproducible without a seed; there is no
#' `rseed` argument. (The local-dependence branch's QUIC solver shuffles its
#' active set with the RNG, but that is a convex solve, so the shuffle order
#' does not change the solution.)
#'
#' Missing-data support and the local-dependence ELBO are planned for a future
#' release; see `NEWS`.
#'
#' @references
#' Chen, J., Guo, Z., Zhang, L., & Pan, J. (2021). A partially confirmatory
#' approach to scale development with the Bayesian Lasso. *Psychological
#' Methods*, 26(2), 210–235. \doi{10.1037/met0000293}
#'
#' Jin, Y., & Chen, J. (2025). Regularized variational approximation for
#' partially confirmatory factor analysis. *Structural Equation Modeling: A
#' Multidisciplinary Journal*, 32(3), 437–449.
#' \doi{10.1080/10705511.2024.2432612}
#'
#' Ročková, V., & George, E. I. (2018). The spike-and-slab LASSO.
#' *Journal of the American Statistical Association*, 113(521), 431–444.
#'
#' @examples
#' ## partially confirmatory: two anchors per factor, the rest exploratory
#' sim <- sim_fa(N = 300, K = 3, ipf = 6, lam = .7, lac = .3, rseed = 1)
#' Q <- matrix(-1L, ncol(sim$dat), 3)
#' for (k in 1:3) { a <- which(rep(1:3, each = 6) == k)[1:2]; Q[a, ] <- 0; Q[a, k] <- 1 }
#' fit <- vbfa(sim$dat, Q, convChk = FALSE)
#' fit$flag                    # 1 if converged
#' round(fit$Lam, 2)           # loadings
#' round(fit$pi, 2)            # posterior inclusion probabilities
#'
#' \donttest{
#' ## orthogonal bifactor: one general column plus the group factors
#' Qb  <- cbind(1L, Q)
#' fb  <- vbfa(sim$dat, Qb, orthogonal = TRUE, convChk = FALSE)
#' fb$Phi                      # identity by construction
#'
#' ## local dependence: simulate correlated residuals, then set ld = TRUE
#' simLD <- sim_fa(N = 500, K = 3, ipf = 6, lam = .7, lac = .3, ecr = .3,
#'                 rseed = 2)
#' fLD <- vbfa(simLD$dat, Q, ld = TRUE, convChk = FALSE, max_it = 300,
#'             tolVal = 1e-3)
#' ## largest recovered residual edges vs the planted ones
#' Poff <- abs(fLD$Psi); Poff[lower.tri(Poff, diag = TRUE)] <- 0
#' which(Poff >= sort(Poff, decreasing = TRUE)[3], arr.ind = TRUE)
#' simLD$ofd_ind
#'
#' ## empirical illustration: NLSY 1997 (27 mixed-type items, 3 factors).
#' ## vbfa() does not support NA yet, so complete cases are used, and the
#' ## polytomous items are treated as continuous for this illustration.
#' data(nlsy27)
#' Yn <- as.matrix(nlsy27$dat[stats::complete.cases(nlsy27$dat), ])
#' fn <- vbfa(Yn, nlsy27$Q, convChk = FALSE)
#' round(fn$Lam, 2)
#' }
#'
#' @export
vbfa <- function(Y, Q, ld = FALSE, Qe = NULL, orthogonal = FALSE,
                 v0 = c(0.01, 0.005, 0.002, 0.001),
                 max_it = 5000, convChk = TRUE, tolVal = 1e-4,
                 ld_control = list()) {

  ## ---- input checks --------------------------------------------------
  ## missing data are not yet supported (a future release will impute inside
  ## the variational loop; see NEWS)
  if (anyNA(Y))
    stop("Y contains NA. vbfa() does not support missing data; ",
         "impute beforehand.", call. = FALSE)
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  if (nrow(Q) != ncol(Y))
    stop("The numbers of items in data and Q are unequal.", call. = FALSE)
  if (!all(Q %in% c(-1, 0, 1)))
    stop("Q entries must be in {-1, 0, 1}.", call. = FALSE)

  ld <- isTRUE(ld)
  if (!ld) {
    if (!is.null(Qe)) {
      warning("Qe is ignored because ld = FALSE; set ld = TRUE to enable ",
              "local-dependence estimation.", call. = FALSE)
      Qe <- NULL
    }
  } else {
    if (is.null(Qe))
      Qe <- matrix(-1, ncol(Y), ncol(Y))   # fully exploratory local dependence
    Qe <- as.matrix(Qe)
    if (nrow(Qe) != ncol(Y) || ncol(Qe) != ncol(Y))
      stop("Qe must be J x J (J = ncol(Y)).", call. = FALSE)
    if (!all(Qe %in% c(-1, 0, 1)))
      stop("Qe entries must be in {-1, 0, 1}.", call. = FALSE)
    if (!isSymmetric(unname(Qe)))
      stop("Qe must be symmetric (its diagonal is ignored).", call. = FALSE)
    quic_fun <- .quic_rcpp   # package-internal QUIC solver (src/quic_impl.cpp)
  }

  ## ---- ld_control defaults -------------------------------------------
  J <- ncol(Y)
  N <- nrow(Y)
  P <- ncol(Q)
  ctrl <- list(xi0 = seq(0.1, 1, length.out = 5),  # units of N
               xi1 = 0.01,                          # units of N
               diag_penalty = 1,
               quic_eps = 1e-3, quic_max_it = 1000,
               a1 = 1, b1 = J)
  if (length(ld_control)) {
    unknown <- setdiff(names(ld_control), names(ctrl))
    if (length(unknown))
      stop("Unknown ld_control entries: ", paste(unknown, collapse = ", "),
           call. = FALSE)
    ctrl[names(ld_control)] <- ld_control
  }

  ## ---- dynamic-path staging -------------------------------------------
  ## stages = max(length(v0), length(xi0)); shorter held at its last value.
  v0_seq <- as.numeric(v0)
  if (any(v0_seq <= 0)) stop("v0 must be positive.", call. = FALSE)
  if (ld) {
    xi0_seq <- as.numeric(ctrl$xi0)
    n_stage <- max(length(v0_seq), length(xi0_seq))
    v0_seq  <- c(v0_seq,  rep(v0_seq[length(v0_seq)],
                              n_stage - length(v0_seq)))
    xi0_seq <- c(xi0_seq, rep(xi0_seq[length(xi0_seq)],
                              n_stage - length(xi0_seq)))
  } else {
    n_stage <- length(v0_seq)
    xi0_seq <- rep(NA_real_, n_stage)
  }

  ## ---- RNG hygiene ---------------------------------------------------
  ## The diagonal-residual estimator consumes no random numbers (deterministic).
  ## The only RNG use is the QUIC active-set shuffle in the local-dependence
  ## branch -- a convex solve, so the order does not change the solution. We
  ## save and restore the caller's RNG state so a fit leaves it untouched.
  if (exists(".Random.seed", .GlobalEnv)) {
    oldseed <- .GlobalEnv$.Random.seed
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv))
  }

  ## ---- shared setup and initialization -------------------------------
  Y <- scale(Y)          # standardize (both modes)

  flag <- 0
  xi.PHI  <- 2 * P + 1
  Lam.PHI <- diag(1, P, P)

  mu.q.Lam <- mu.q.Lam_sq <- matrix(0.8, J, P)
  mu.q.Lam[Q == -1] <- mu.q.Lam_sq[Q == -1] <- 0.2
  mu.q.Lam[Q == 0]  <- mu.q.Lam_sq[Q == 0]  <- 0
  sigsq.q.Lam <- matrix(1, J, P)

  mu.q.eta <- mu.q.eta_sq <- matrix(1, N, P)
  mu.q.d <- matrix(1, J, P)
  v1 <- 1
  inv_v1 <- 1 / v1
  h_0 <- 1
  gamma.q.d <- rep(1, J)

  PHI.q.eta <- M.q.LamTLam <- matrix(1, P, P)
  M.q.PHIiver <- diag(1, P, P)   # orthogonal: stays I throughout (Phi = I)

  xi.q.PHI <- N + xi.PHI

  a <- 1
  b <- 0.1
  a.q.psi <- rep(N / 2 + a, J)
  b.q.psi <- rep(b, J)
  mu.q.psiInv <- rep(1, J)
  plotDat <- NULL

  a0 <- 1
  b0 <- 1

  q <- Q
  rho <- 0.5

  mu_tem2 <- matrix(0, J, P)

  ## ---- LD-specific state ------------------------------------------------
  if (ld) {
    Psi <- W <- diag(1, J)       # full residual PRECISION (Psi) and its
    tau <- 0.5                   # inverse W; tau = LD proportion
    q_star <- xi_star <- matrix(0.5, J, J)
    xi1_abs <- ctrl$xi1 * N
    a1 <- ctrl$a1
    b1 <- ctrl$b1
    upper_tri_idx  <- which(upper.tri(Psi), arr.ind = TRUE)
    uncertain_upper <- upper.tri(Qe) & Qe == -1
  }

  ## ---- staged VB loop ----------------------------------------------------
  itNum <- 0
  stage_iters <- integer(n_stage)
  MFVBstartTime <- Sys.time()

  for (stage in seq_len(n_stage)) {
    v0 <- v0_seq[stage]
    inv_v0 <- 1 / v0
    if (ld) xi0_abs <- xi0_seq[stage] * N
    converged <- FALSE
    parVecHatPrev <- NULL
    stage_it <- 0
    flag <- 0        # reset per stage; reflects the FINAL stage on return
    if (convChk && n_stage > 1)
      cat("Stage", stage, "/", n_stage, ": v0 =", signif(v0, 4),
          if (ld) paste0("; xi0 = ", signif(xi0_seq[stage], 4), " * N"),
          "\n")

    while (!converged) {
      itNum <- itNum + 1
      stage_it <- stage_it + 1

      ## -- eta update ---------------------------------------------------
      if (!ld) {
        M.q.LamTpsiInvLam <- t(mu.q.Lam) %*% diag(mu.q.psiInv, J, J) %*% mu.q.Lam
        diag(M.q.LamTpsiInvLam) <- apply(mu.q.Lam_sq * mu.q.psiInv, 2, sum)
        PHI.q.eta <- solve(M.q.LamTpsiInvLam + M.q.PHIiver)
        mu.q.eta <- t(PHI.q.eta %*% t(Y %*% (mu.q.Lam * mu.q.psiInv)))
      } else {
        M.q.LamTpsiInvLam <- t(mu.q.Lam) %*% Psi %*% mu.q.Lam
        diag(M.q.LamTpsiInvLam) <- diag(M.q.LamTpsiInvLam) +
          apply(diag(Psi) * (mu.q.Lam_sq - mu.q.Lam^2), 2, sum)
        PHI.q.eta <- solve(M.q.LamTpsiInvLam + M.q.PHIiver)
        mu.q.eta <- t(PHI.q.eta %*% t(Y %*% Psi %*% mu.q.Lam))
      }

      ## standardize factor variance to 1 via the posterior second moment
      ## E_q[eta_p^2] = mean(eta_p^2) + Var_q(eta_p)   (both modes)
      scale_k <- sqrt(colMeans(mu.q.eta^2) + diag(PHI.q.eta))
      mu.q.eta <- t(t(mu.q.eta) / scale_k)
      PHI.q.eta <- PHI.q.eta / outer(scale_k, scale_k)
      mu.q.eta_sq <- sweep(mu.q.eta^2, 2, diag(PHI.q.eta), '+')

      E_eta_eta <- t(mu.q.eta) %*% mu.q.eta + N * PHI.q.eta

      if (!orthogonal) {               # oblique: inverse-Wishart update
        Lam.q.PHI <- E_eta_eta + Lam.PHI
        M.q.PHIiver <- (xi.q.PHI - P + 1) * solve(Lam.q.PHI)
      }                                # orthogonal: M.q.PHIiver stays I

      etaY <- t(mu.q.eta) %*% Y

      ## -- loading (and, non-LD, uniqueness) updates ----------------------
      if (!ld) {
        for (j in 1:J) {
          E_lam_lam_j <- tcrossprod(mu.q.Lam[j, ])
          diag(E_lam_lam_j) <- mu.q.Lam_sq[j, ]
          RSS_j <- sum(Y[, j]^2) -
            2 * sum(etaY[, j] * mu.q.Lam[j, ]) +
            sum(E_eta_eta * E_lam_lam_j)
          b.q.psi[j] <- RSS_j / 2 + b
          mu.q.psiInv[j] <- a.q.psi[j] / b.q.psi[j]
          for (p in 1:P) {
            cross_part <- sum(E_eta_eta[p, -p] * mu.q.Lam[j, -p])
            if (Q[j, p] == 1) {
              sigsq.q.Lam[j, p] <- (E_eta_eta[p, p] * mu.q.psiInv[j] + h_0)^(-1)
              sigsq.q.Lam[j, p] <- max(sigsq.q.Lam[j, p], 1e-6)
              sigsq.q.Lam[j, p] <- min(sigsq.q.Lam[j, p], 1e6)
              mu.q.Lam[j, p] <- sigsq.q.Lam[j, p] * mu.q.psiInv[j] *
                (etaY[p, j] - cross_part)
              mu.q.Lam_sq[j, p] <- mu.q.Lam[j, p]^2 + sigsq.q.Lam[j, p]
              q[j, p] <- 1
            } else if (Q[j, p] == -1) {
              tem <- 1/2 * log(inv_v1 / inv_v0) + log(rho / (1 - rho)) +
                1/2 * mu.q.Lam_sq[j, p] * (inv_v0 - inv_v1)
              tem <- pmin(tem, 100)
              q[j, p] <- exp(tem) / (1 + exp(tem))
              mu_tem2[j, p] <- (1 - q[j, p]) * inv_v0 + q[j, p] * inv_v1
              sigsq.q.Lam[j, p] <- (E_eta_eta[p, p] * mu.q.psiInv[j] +
                                      mu_tem2[j, p])^(-1)
              sigsq.q.Lam[j, p] <- max(sigsq.q.Lam[j, p], 1e-6)
              sigsq.q.Lam[j, p] <- min(sigsq.q.Lam[j, p], 1e6)
              mu.q.Lam[j, p] <- sigsq.q.Lam[j, p] * mu.q.psiInv[j] *
                (etaY[p, j] - cross_part)
              mu.q.Lam_sq[j, p] <- mu.q.Lam[j, p]^2 + sigsq.q.Lam[j, p]
            }
          }
        }
      } else {
        ## LD: loading update gains a cross-item term (tem2) because the
        ## residual precision Psi is non-diagonal
        for (j in 1:J) {
          for (p in 1:P) {
            if (Q[j, p] == 0) next
            if (Q[j, p] == 1) {
              sigsq.q.Lam[j, p] <- (E_eta_eta[p, p] * Psi[j, j] + h_0)^(-1)
              sigsq.q.Lam[j, p] <- max(sigsq.q.Lam[j, p], 1e-6)
              sigsq.q.Lam[j, p] <- min(sigsq.q.Lam[j, p], 1e6)
              tem1 <- Psi[j, j] * (etaY[p, j] -
                                     sum(E_eta_eta[p, -p] * mu.q.Lam[j, -p]))
              tem2_vec <- etaY[p, -j] -
                as.vector(mu.q.Lam[-j, , drop = FALSE] %*% E_eta_eta[p, ])
              tem2 <- Psi[j, -j] %*% matrix(tem2_vec, ncol = 1)
              mu.q.Lam[j, p] <- sigsq.q.Lam[j, p] * (tem1 + tem2)
              mu.q.Lam_sq[j, p] <- mu.q.Lam[j, p]^2 + sigsq.q.Lam[j, p]
              q[j, p] <- 1
            } else if (Q[j, p] == -1) {
              tem <- 1/2 * log(inv_v1 / inv_v0) + log(rho / (1 - rho)) +
                1/2 * mu.q.Lam_sq[j, p] * (inv_v0 - inv_v1)
              tem <- pmin(pmax(tem, -100), 100)
              q[j, p] <- 1 / (1 + exp(-tem))
              mu_tem2[j, p] <- (1 - q[j, p]) * inv_v0 + q[j, p] * inv_v1
              sigsq.q.Lam[j, p] <- (E_eta_eta[p, p] * Psi[j, j] +
                                      mu_tem2[j, p])^(-1)
              sigsq.q.Lam[j, p] <- max(sigsq.q.Lam[j, p], 1e-6)
              sigsq.q.Lam[j, p] <- min(sigsq.q.Lam[j, p], 1e6)
              tem1 <- Psi[j, j] * (etaY[p, j] -
                                     sum(E_eta_eta[p, -p] * mu.q.Lam[j, -p]))
              tem2_vec <- etaY[p, -j] -
                as.vector(mu.q.Lam[-j, , drop = FALSE] %*% E_eta_eta[p, ])
              tem2 <- Psi[j, -j] %*% matrix(tem2_vec, ncol = 1)
              mu.q.Lam[j, p] <- sigsq.q.Lam[j, p] * (tem1 + tem2)
              mu.q.Lam_sq[j, p] <- mu.q.Lam[j, p]^2 + sigsq.q.Lam[j, p]
            }
          }
        }
      }

      ## -- rho update (shared) --------------------------------------------
      a0_q_rho <- a0 + sum(q[Q == -1])
      b0_q_rho <- b0 + sum(Q == -1) - sum(q[Q == -1])
      rho <- a0_q_rho / (sum(Q == -1) + a0 + b0)

      ## -- LD block: q_star / tau / S / QUIC -------------------------------
      if (ld) {
        for (i in seq_len(nrow(upper_tri_idx))) {
          jr <- upper_tri_idx[i, 1]
          jc <- upper_tri_idx[i, 2]
          if (Qe[jr, jc] == -1) {
            tem <- plogis(log(tau / (1 - tau)) + log(xi1_abs / xi0_abs) +
                            (xi0_abs - xi1_abs) * abs(Psi[jr, jc]))
          } else if (Qe[jr, jc] == 1) {
            tem <- 1
          } else {
            tem <- 0
          }
          q_star[jr, jc] <- tem
          q_star[jc, jr] <- tem
        }

        xi_star <- xi1_abs * q_star + xi0_abs * (1 - q_star)
        if (ctrl$diag_penalty == 0) {
          diag(xi_star) <- 0
        } else {
          diag(xi_star) <- xi1_abs
        }

        if (sum(uncertain_upper) > 0) {
          tau <- (a1 - 1 + sum(q_star[uncertain_upper])) /
            (a1 + b1 - 2 + sum(uncertain_upper))
        }

        xi_star <- xi_star / N

        ## expected residual second-moment matrix S
        tem2m <- mu.q.Lam %*% E_eta_eta %*% t(mu.q.Lam)
        diag(tem2m) <- diag(tem2m) +
          rowSums(sigsq.q.Lam * matrix(diag(E_eta_eta), J, P, byrow = TRUE))
        cross_term <- t(Y) %*% mu.q.eta %*% t(mu.q.Lam)
        S <- (tem2m + t(Y) %*% Y - cross_term - t(cross_term)) / N

        res_quic <- quic_fun(J, S, xi_star,
                             eps = ctrl$quic_eps,
                             max_iter_quic = ctrl$quic_max_it)
        Psi <- res_quic[, , 1]
        W   <- res_quic[, , 2]
      }

      ## -- convergence check -----------------------------------------------
      ## (the two expressions coincide when Psi is diagonal; the non-LD form
      ## is kept verbatim for bit-exact reproduction of the fixed-spike fit)
      if (!ld) {
        parVecHat <- c(t(t(abs(Y - mu.q.eta %*% t(mu.q.Lam))) * mu.q.psiInv))
      } else {
        parVecHat <- c(abs(Y - mu.q.eta %*% t(mu.q.Lam)) %*% Psi)
      }

      if (!is.null(parVecHatPrev)) {
        relErrTemp <- abs(parVecHat - parVecHatPrev)
        relErrTemp[is.nan(relErrTemp)] <- 0
        relErr <- max(relErrTemp)
        plotDat <- append(plotDat, relErr)
        if (convChk == TRUE)
          cat("Iteration number:", itNum, ";",
              "MFVB relative error:", signif(relErr, 5), "\n")
        if (relErr < tolVal) {
          converged <- TRUE
          flag <- 1
        }
      }
      if (stage_it >= max_it) {
        converged <- TRUE
        warning("Maximum number of MFVB iterations exceeded",
                if (n_stage > 1) paste0(" at stage ", stage), ".\n")
      }
      parVecHatPrev <- parVecHat
    }
    stage_iters[stage] <- stage_it
  }

  ## ---- ELBO ---------------------------------------------------------------
  if (!ld) {
    q_uns <- q[Q == -1]
    q_uns <- pmin(pmax(q_uns, 1e-12), 1 - 1e-12)
    if (orthogonal) {
      eta_phi_term <- -0.5 * sum(diag(E_eta_eta)) +
        0.5 * N * log(det(PHI.q.eta)) +
        0.5 * N * P
    } else {
      eta_phi_term <- 0.5 * N * log(det(PHI.q.eta)) +
        0.5 * N * P +
        logC_igw(xi.PHI, Lam.PHI, P) -
        logC_igw(xi.q.PHI, Lam.q.PHI, P)
    }
    ELBO <- (
      -0.5 * N * J * log(2 * pi) +
        sum(a * log(b) - lgamma(a) -
              a.q.psi * log(b.q.psi) + lgamma(a.q.psi)) +
        eta_phi_term +
        0.5 * sum(
          log(sigsq.q.Lam[Q == 1] / h_0) + 1 -
            mu.q.Lam_sq[Q == 1] / h_0
        ) +
        0.5 * sum(
          log(sigsq.q.Lam[Q == -1]) + 1 -
            (1 - q_uns) * log(v0) -
            q_uns * log(v1) -
            mu.q.Lam_sq[Q == -1] *
            ((1 - q_uns) / v0 + q_uns / v1)
        ) +
        lbeta(a0_q_rho, b0_q_rho) -
        lbeta(a0, b0) -
        sum(q_uns * log(q_uns) + (1 - q_uns) * log(1 - q_uns))
    )
  } else {
    ## the ELBO is not yet available under local dependence; a future release
    ## will assemble it from the graphical spike-and-slab prior over Psi's
    ## off-diagonals, log det Psi (via W from QUIC), the q_star entropy, and the
    ## Beta(tau) terms.
    ELBO <- NA_real_
  }

  ## ---- sign check ----------------------------------------------------
  mean_mu.q.Lam <- apply(mu.q.Lam, 2, mean)
  sign_flag <- which(mean_mu.q.Lam < 0)
  if (length(sign_flag) > 0) {
    mu.q.Lam[, sign_flag] <- -mu.q.Lam[, sign_flag]
    mu.q.eta[, sign_flag] <- -mu.q.eta[, sign_flag]
    sign_vec <- rep(1, ncol(mu.q.Lam))
    sign_vec[sign_flag] <- -1
    PHI.q.eta <- PHI.q.eta * outer(sign_vec, sign_vec)
    M.q.PHIiver <- M.q.PHIiver * outer(sign_vec, sign_vec)
    if (!orthogonal) Lam.q.PHI <- Lam.q.PHI * outer(sign_vec, sign_vec)
  }

  MFVBendTime <- Sys.time()
  MFVBtime <- MFVBendTime - MFVBstartTime

  out <- list(sigsq.q.Lam = sigsq.q.Lam, iter = itNum,
              PHI.q.eta = PHI.q.eta, eta = mu.q.eta,
              Phi = if (orthogonal) diag(1, P, P) else cov2cor(Lam.q.PHI),
              M.q.PHIiver = M.q.PHIiver, time = MFVBtime,
              PsiInv = if (!ld) mu.q.psiInv else NULL,
              plotDat = plotDat, flag = flag, Lam = mu.q.Lam, pi = q,
              ELBO = ELBO, rho = rho,
              ## settings the model was fit with, so downstream consumers
              ## (vb_fit) need not be told twice
              orthogonal = orthogonal, ld = ld,
              ## LD extras (NULL when Qe = NULL)
              Psi = if (ld) Psi else NULL,
              W = if (ld) W else NULL,
              q_star = if (ld) q_star else NULL,
              xi_star = if (ld) xi_star else NULL,
              tau = if (ld) tau else NULL,
              ## dynamic-path metadata
              path = list(v0 = v0_seq,
                          xi0 = if (ld) xi0_seq else NULL,
                          n_stage = n_stage, stage_iters = stage_iters))
  return(out)
}
