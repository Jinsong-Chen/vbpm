## Regularized variational Bayes for the extended MIMIC model.
##
## Ported from the estimator accompanying Jin & Chen (2025, Multivariate
## Behavioral Research). The variational updates are unchanged from that
## implementation; the packaging changes are input validation, namespace
## imports, a staged (warm-started) v0 path, a convergence flag, an optional
## standardization step, and vbfa()'s sign convention on the output. The
## defaults (scalar v0 = 0.001, standardize = FALSE) reproduce the published
## estimator exactly.

#' Regularized variational Bayes for extended MIMIC models
#'
#' Fits a partially confirmatory multiple-indicators multiple-causes (MIMIC)
#' model by a variational Bayes EM algorithm, with continuous spike-and-slab
#' regularization on **both** the measurement and the structural part. Items
#' load on factors through a measurement design matrix `Q_A`, while covariates
#' predict the factors through a structural design matrix `Q_B`; entries of
#' either that are left unspecified are selected from the data.
#'
#' This is the MIMIC member of the `vbpm` family and is the companion of
#' [vbfa()], which fits the measurement part alone. The two share the `-1/0/1`
#' design-matrix convention and the spike-and-slab formulation.
#'
#' @param Y Numeric `N x J` data matrix (observations by items). Missing values
#'   in `Y` are supported: because the residual covariance is diagonal, a
#'   missing response is replaced in-loop by its conditional expectation
#'   `eta_i A'`, with the conditional variance carried into the residual sum
#'   of squares. Items or rows with no observed response are rejected.
#' @param X Numeric `N x P` matrix of observed covariates (predictors of the
#'   factors), with the same number of rows as `Y`. Missing covariates are
#'   **not** supported: `X` is conditioned on rather than modelled, so
#'   imputing it would require a distributional assumption the MIMIC model
#'   does not make.
#' @param Q_A Integer `J x K` measurement design matrix: `1` = specified
#'   (anchored) loading, `0` = fixed zero, `-1` = unspecified (selected by
#'   spike-and-slab). `K` is the number of factors.
#' @param Q_B Integer `K x P` structural design matrix using the same coding,
#'   for the regression of the `K` factors on the `P` covariates.
#' @param v0 Spike variance. A scalar gives a single fixed spike; a decreasing
#'   vector gives a warm-started regularization path, each stage initialized at
#'   the previous stage's solution (Rockova & George, 2018). The same schedule
#'   is applied to the measurement and structural parts. `v0 = 0.001` (the
#'   default) is the fixed spike of the published estimator.
#' @param standardize Logical; standardize `Y` and `X` before fitting. The
#'   default `FALSE` is the behaviour of the published estimator. Note that
#'   [vbfa()] always standardizes internally; set `standardize = TRUE` for the
#'   analogous behaviour here.
#' @param max_it Maximum number of variational iterations *per stage*.
#' @param convChk Logical; print per-iteration convergence information via
#'   [message()]. Default `FALSE` (quiet).
#' @param tolVal Convergence tolerance on the maximum absolute change in the
#'   weighted residual vector.
#'
#' @return An object of class `c("vbmimic", "vbpm_fit")` -- a named list (see
#'   [vbpm_fit]) with components:
#'   \describe{
#'     \item{A}{`J x K` posterior mean measurement loadings.}
#'     \item{B}{`K x P` posterior mean structural coefficients.}
#'     \item{pi_A, pi_B}{Posterior inclusion probabilities for the unspecified
#'       entries of `Q_A` and `Q_B` (specified entries are held at 1, fixed
#'       zeros at 0).}
#'     \item{Q_A, Q_B}{The design matrices the model was fit with.}
#'     \item{eta}{`N x K` posterior mean factor scores.}
#'     \item{Phi}{`K x K` factor correlation matrix implied by `Sig`.}
#'     \item{Sig, U}{Factor covariance matrix and its inverse.}
#'     \item{V}{Length-`J` residual precisions.}
#'     \item{rho, theta}{Inclusion rates for the measurement and structural
#'       parts.}
#'     \item{A_var, B_var}{Posterior variances of the loadings and structural
#'       coefficients.}
#'     \item{iter, flag, time}{Total iterations, `1` if the final stage
#'       converged within `max_it` (else `0`), and elapsed time.}
#'     \item{ELBO}{`NA_real_`; see Note.}
#'     \item{standardize}{The `standardize` setting the model was fit with.}
#'     \item{path}{The `v0` schedule and per-stage iteration counts.}
#'     \item{preprocess}{`n_missing`, the response type, and `missing_mask` --
#'       stored only when the data contain missing values, `NULL` otherwise.}
#'   }
#'
#'   Every quantity appears exactly once: as of 0.7.0 the fit no longer
#'   carries `$coefficients`, `$design`, `$posterior`, or `$settings`, which
#'   duplicated the elements above.
#'
#' @note The evidence lower bound is **not** currently returned. The published
#'   algorithm converges on the weighted residual vector rather than on the
#'   bound, and the bound's implementable closed form has not been assembled.
#'   [fit_stats()] does accept a `vbmimic` fit, but returns a deliberately
#'   limited result: `n_active_coef` (soft-selected measurement and structural
#'   coefficients) plus typed `NA` for every covariance-based index and
#'   parameter count, with an explanatory `note` attribute. [pefa()] does not
#'   accept `vbmimic` fits at all, because factor-number selection needs a
#'   comparable objective.
#'
#' @references
#' Jin, Y., & Chen, J. (2025). Regularized variational Bayesian approximations
#' for variable selection in extended multiple-indicators multiple-causes
#' models. *Multivariate Behavioral Research*.
#' \doi{10.1080/00273171.2025.2483253}
#'
#' Rockova, V., & George, E. I. (2018). The spike-and-slab LASSO.
#' *Journal of the American Statistical Association*, 113(521), 431-444.
#'
#' @seealso [vbfa()] for the measurement model alone, [sim_lvm()] to simulate
#'   MIMIC data.
#'
#' @examples
#' ## 3 factors, 18 items, 9 covariates in disjoint blocks of 3
#' B <- matrix(0, 3, 9)
#' for (k in 1:3) B[k, ((k - 1) * 3 + 1):(k * 3)] <- .3
#' sim <- sim_lvm(N = 300, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 1)
#' Y <- sim$dat[, 1:18]
#' X <- sim$dat[, 19:27]
#'
#' ## two anchored items per factor; the structural part left exploratory.
#' ## Anchoring at least one part is advisable: with both Q_A and Q_B fully
#' ## exploratory the solution can be rotationally ambiguous.
#' Q_A <- matrix(-1L, 18, 3)
#' for (k in 1:3) {
#'   a <- which(rep(1:3, each = 6) == k)[1:2]
#'   Q_A[a, ] <- 0L; Q_A[a, k] <- 1L
#' }
#' Q_B <- matrix(-1L, 3, 9)
#'
#' fit <- vbmimic(Y, X, Q_A, Q_B)
#' fit$flag           # 1 if converged
#' round(fit$B, 2)    # structural coefficients
#'
#' @export
vbmimic <- function(Y, X, Q_A, Q_B, v0 = 0.001, standardize = FALSE,
                    max_it = 1000, convChk = FALSE, tolVal = 1e-5) {

  call <- match.call()

  ## ---- input checks --------------------------------------------------
  ## Missing responses in Y are supported (see below). Missing covariates are
  ## not: X is conditioned on rather than modelled, so imputing it would
  ## require assuming a distribution for the covariates that the MIMIC model
  ## does not specify.
  if (anyNA(X))
    stop("X contains NA. vbmimic() does not support missing covariates: X is ",
         "conditioned on, not modelled. Impute or drop them beforehand.",
         call. = FALSE)

  Y   <- as.matrix(Y)
  X   <- as.matrix(X)
  Q_A <- as.matrix(Q_A)
  Q_B <- as.matrix(Q_B)

  if (nrow(X) != nrow(Y))
    stop("Y and X must have the same number of rows (observations).",
         call. = FALSE)
  if (nrow(Q_A) != ncol(Y))
    stop("nrow(Q_A) must equal ncol(Y) (the number of items).", call. = FALSE)
  if (ncol(Q_B) != ncol(X))
    stop("ncol(Q_B) must equal ncol(X) (the number of covariates).",
         call. = FALSE)
  if (nrow(Q_B) != ncol(Q_A))
    stop("nrow(Q_B) must equal ncol(Q_A) (the number of factors).",
         call. = FALSE)
  if (!all(Q_A %in% c(-1, 0, 1)))
    stop("Q_A entries must be in {-1, 0, 1}.", call. = FALSE)
  if (!all(Q_B %in% c(-1, 0, 1)))
    stop("Q_B entries must be in {-1, 0, 1}.", call. = FALSE)

  v0_seq <- as.numeric(v0)
  if (any(v0_seq <= 0)) stop("v0 must be positive.", call. = FALSE)
  n_stage <- length(v0_seq)

  J <- ncol(Y)
  N <- nrow(Y)
  K <- ncol(Q_A)
  P <- ncol(X)
  if (K < 2)
    stop("vbmimic() requires at least two factors (ncol(Q_A) >= 2).",
         call. = FALSE)

  ## ---- missing responses ----------------------------------------------
  ## The residual covariance of Y is diagonal here (V holds residual
  ## precisions), so the joint Gaussian conditional of vbfa() collapses to its
  ## diagonal case: E[y_ij | y_i,obs] = eta_i A_j'. The covariate information
  ## is already carried by eta_i, whose own update uses X and B, so no extra
  ## conditioning on X is required. The conditional variance of an imputed
  ## cell is 1 / V_j and is added to the residual sum of squares below.
  missing_mask <- is.na(Y)
  has_missing  <- any(missing_mask)
  if (has_missing) {
    if (any(colSums(!missing_mask) == 0L))
      stop("Every item needs at least one observed response.", call. = FALSE)
    if (any(rowSums(!missing_mask) == 0L))
      stop("Rows with all responses missing are not supported.", call. = FALSE)
  }
  miss_count <- colSums(missing_mask)

  if (standardize) {
    ## centre and scale from observed entries only
    center <- colMeans(Y, na.rm = TRUE)
    scale_ <- apply(Y, 2, stats::sd, na.rm = TRUE)
    Y <- sweep(sweep(Y, 2, center, "-"), 2, scale_, "/")
    X <- scale(X)
  }
  Y[missing_mask] <- 0   # replaced by the conditional mean in iteration 1

  ## ---- initialization -------------------------------------------------
  ## Deterministic: the estimator consumes no random numbers, so there is no
  ## seed argument (cf. vbfa()).
  mu.q.A <- mu.q.A_sq <- matrix(0.8, J, K)
  mu.q.A[Q_A == -1] <- mu.q.A_sq[Q_A == -1] <- 0.2
  mu.q.A[Q_A == 0]  <- mu.q.A_sq[Q_A == 0]  <- 0
  sigsq.q.A <- matrix(1, J, K)
  Alpha <- matrix(0.5, J, K)
  Alpha[Q_A == 0] <- 0
  Alpha[Q_A == 1] <- 1
  RHO <- 0.5

  mu.q.B <- mu.q.B_sq <- matrix(0.8, K, P)
  mu.q.B[Q_B == -1] <- mu.q.B_sq[Q_B == -1] <- 0.2
  mu.q.B[Q_B == 0]  <- mu.q.B_sq[Q_B == 0]  <- 0
  sigsq.q.B <- matrix(1, K, P)
  Beta <- matrix(0.5, K, P)
  Beta[Q_B == 0] <- 0
  Beta[Q_B == 1] <- 1
  Theta <- 0.5

  mu.F <- mu.F_sq <- matrix(1, N, K)
  sigsq.F <- matrix(1, N, K)
  V <- rep(1, J)
  U <- Sig <- diag(1, K, K)

  invH_a <- invH_b <- 1
  a_v <- 1
  b_v <- 0.1
  inv_va1 <- inv_vb1 <- 1      # slab precision (v1 = 1)
  a0 <- b0 <- 1

  ## ---- staged VB loop ---------------------------------------------------
  itNum <- 0
  flag <- 0
  stage_iters <- integer(n_stage)
  startTime <- Sys.time()

  for (stage in seq_len(n_stage)) {
    ## spike precision for this stage, shared by the measurement and
    ## structural parts
    inv_va0 <- inv_vb0 <- 1 / v0_seq[stage]
    if (convChk && n_stage > 1)
      message("Stage ", stage, "/", n_stage, ": v0 = ", signif(v0_seq[stage], 4))

    converged <- FALSE
    parVecHatPrev <- NULL
    stage_it <- 0
    flag <- 0        # reset per stage; reflects the FINAL stage on return

    while (!converged) {
      itNum <- itNum + 1
      stage_it <- stage_it + 1

      ## ---- factors ------------------------------------------------------
      for (k in 1:K) {
        sigsq.F[, k] <- 1 / (U[k, k] + sum(mu.q.A_sq[, k] * V))
        tem1 <- (Y - mu.F[, -k, drop = FALSE] %*% t(mu.q.A[, -k, drop = FALSE])) %*%
          matrix(mu.q.A[, k] * V, nrow = J, ncol = 1)
        tem2 <- X %*% matrix(mu.q.B[k, ], nrow = P, ncol = 1) * U[k, k]
        tem3 <- apply(t(t(mu.F[, -k, drop = FALSE]) * U[k, -k]), 1, sum)
        tem4 <- apply(t(t(X %*% t(mu.q.B[-k, , drop = FALSE])) * U[k, -k]), 1, sum)
        mu.F[, k] <- sigsq.F[, k] * (tem1 + tem2 - tem3 + tem4)
      }

      ## rescale factors to unit variance (identification)
      sdF <- apply(mu.F, 2, sd)
      mu.F <- t(t(mu.F) / sdF)
      mu.F_sq <- mu.F^2 + sigsq.F

      ## ---- missing responses: conditional mean and variance --------------
      ## Diagonal residual covariance => Equation (2) of the design note.
      ## V from the start of the sweep is used for the variance correction,
      ## as in vbfa()'s coordinate-ascent ordering.
      if (has_missing) {
        mean_now <- mu.F %*% t(mu.q.A)
        Y[missing_mask] <- mean_now[missing_mask]
        miss_var <- miss_count / V           # length J; 1/V_j per missing cell
      }

      ## ---- residual precisions and measurement loadings ------------------
      for (j in 1:J) {
        ## expressions kept in the published associativity so the fit is
        ## bit-identical to the reference implementation; do not "simplify"
        ## without re-running the equivalence test
        tem1 <- matrix(Y[, j], nrow = 1, ncol = N) %*%
          matrix(Y[, j], nrow = N, ncol = 1)
        tem2 <- Y[, j] %*% (mu.F %*% matrix(mu.q.A[j, ], nrow = K, ncol = 1))
        tem3 <- (mu.q.A[j, ]) %*% t(mu.F) %*% mu.F %*%
          matrix(mu.q.A[j, ], nrow = K, ncol = 1)
        sub <- sum(t(mu.F^2) * mu.q.A[j, ]^2)
        add <- sum(t(mu.F_sq) * (mu.q.A_sq[j, ]))
        tem4 <- (tem1 - 2 * tem2 + tem3 - sub + add)
        ## E[y_ij^2] for an imputed cell is (conditional mean)^2 + 1/V_j; the
        ## squared mean is already in tem1, so only the variance is added.
        if (has_missing) tem4 <- tem4 + miss_var[j]
        V[j] <- (N / 2 + a_v) / (tem4 / 2 + b_v)

        for (k in 1:K) {
          if (Q_A[j, k] == -1) {
            tem <- 1 / 2 * log(inv_va1 / inv_va0) + log(RHO / (1 - RHO)) +
              1 / 2 * mu.q.A_sq[j, k] * (inv_va0 - inv_va1)
            tem <- pmin(tem, 100)
            Alpha[j, k] <- exp(tem) / (1 + exp(tem))
            tem2b <- (1 - Alpha[j, k]) * inv_va0 + Alpha[j, k] * inv_va1
            sigsq.q.A[j, k] <- 1 / (V[j] * sum(mu.F_sq[, k]) + tem2b)
            sigsq.q.A[j, k] <- max(sigsq.q.A[j, k], 1e-6)
            sigsq.q.A[j, k] <- min(sigsq.q.A[j, k], 1e6)
            tem1b <- Y[, j] - (mu.q.A[j, -k]) %*% t(mu.F[, -k, drop = FALSE])
            mu.q.A[j, k] <- sigsq.q.A[j, k] * V[j] *
              (matrix(mu.F[, k], nrow = 1, ncol = N) %*%
                 matrix(tem1b, nrow = N, ncol = 1))
            mu.q.A_sq[j, k] <- mu.q.A[j, k]^2 + sigsq.q.A[j, k]
          } else if (Q_A[j, k] == 1) {
            sigsq.q.A[j, k] <- 1 / (V[j] * sum(mu.F_sq[, k]) + invH_a)
            sigsq.q.A[j, k] <- max(sigsq.q.A[j, k], 1e-6)
            sigsq.q.A[j, k] <- min(sigsq.q.A[j, k], 1e6)
            tem1b <- Y[, j] - (mu.q.A[j, -k]) %*% t(mu.F[, -k, drop = FALSE])
            mu.q.A[j, k] <- sigsq.q.A[j, k] * V[j] *
              (matrix(mu.F[, k], nrow = 1, ncol = N) %*%
                 matrix(tem1b, nrow = N, ncol = 1))
            mu.q.A_sq[j, k] <- mu.q.A[j, k]^2 + sigsq.q.A[j, k]
          }
        }
      }
      RHO <- (a0 + sum(Alpha[Q_A == -1])) / (sum(Q_A == -1) + a0 + b0)

      ## ---- structural coefficients ---------------------------------------
      for (k in 1:K) {
        for (p in 1:P) {
          if (Q_B[k, p] == -1 || Q_B[k, p] == 1) {
            if (Q_B[k, p] == -1) {
              tem <- 1 / 2 * log(inv_vb1 / inv_vb0) + log(Theta / (1 - Theta)) +
                1 / 2 * mu.q.B_sq[k, p] * (inv_vb0 - inv_vb1)
              tem <- pmin(tem, 100)
              Beta[k, p] <- exp(tem) / (1 + exp(tem))
              pen <- (1 - Beta[k, p]) * inv_vb0 + Beta[k, p] * inv_vb1
            } else {
              pen <- invH_b
            }
            xtx <- matrix(X[, p], nrow = 1, ncol = N) %*%
              matrix(X[, p], nrow = N, ncol = 1)
            sigsq.q.B[k, p] <- 1 / (U[k, k] * xtx + pen)
            sigsq.q.B[k, p] <- max(sigsq.q.B[k, p], 1e-6)
            sigsq.q.B[k, p] <- min(sigsq.q.B[k, p], 1e6)
            tem1 <- sum(matrix(X[, p], nrow = 1, ncol = N) %*% mu.F * U[k, ])
            tem2 <- sum((mu.q.B[, -p, drop = FALSE]) %*%
                          t(X[, -p, drop = FALSE]) %*%
                          matrix(X[, p], nrow = N, ncol = 1) * U[k, ])
            tem3 <- sum(U[k, -k] * mu.q.B[-k, p]) * xtx
            mu.q.B[k, p] <- sigsq.q.B[k, p] * (tem1 - tem2 - tem3)
            mu.q.B_sq[k, p] <- mu.q.B[k, p]^2 + sigsq.q.B[k, p]
          }
        }
      }
      Theta <- (a0 + sum(Beta[Q_B == -1])) / (sum(Q_B == -1) + a0 + b0)

      ## ---- factor covariance ----------------------------------------------
      for (k in 1:K) {
        rk <- mu.F[, k] - X %*% matrix(mu.q.B[k, ], nrow = P, ncol = 1)
        tem1 <- matrix(rk, nrow = 1, ncol = N) %*% matrix(rk, nrow = N, ncol = 1)
        tem2 <- sum(mu.F_sq[, k] - mu.F[, k]^2) +
          apply(X^2, 2, sum) %*% matrix(mu.q.B_sq[k, ] - mu.q.B[k, ]^2,
                                        nrow = P, ncol = 1)
        Sig[k, k] <- (tem1 + tem2) / N
      }
      for (k in 1:(K - 1)) {
        for (t in (k + 1):K) {
          rk <- mu.F[, k] - X %*% matrix(mu.q.B[k, ], nrow = P, ncol = 1)
          rt <- mu.F[, t] - X %*% matrix(mu.q.B[t, ], nrow = P, ncol = 1)
          Sig[k, t] <- Sig[t, k] <-
            matrix(rk, nrow = 1, ncol = N) %*% matrix(rt, nrow = N, ncol = 1) / N
        }
      }
      U <- solve(Sig)

      ## ---- convergence -----------------------------------------------------
      ## weighted absolute residuals, in the published form. `%*% diag(V)`
      ## builds a J x J matrix and is wasteful for large J, but t(t(R) * V) --
      ## although mathematically identical -- sums in a different order and
      ## breaks bit-equivalence with the reference. Kept deliberately.
      parVecHat <- c(t(t(abs(Y - mu.F %*% t(mu.q.A)) %*% diag(V))))

      if (stage_it > 1) {
        relErrTemp <- abs(parVecHat - parVecHatPrev)
        relErrTemp[is.nan(relErrTemp)] <- 0
        relErr <- max(relErrTemp)
        if (convChk)
          message("Iteration ", itNum, "; MFVB relative error: ",
                  signif(relErr, 5))
        if (relErr < tolVal) {
          converged <- TRUE
          flag <- 1
        }
      }

      if (stage_it >= max_it) {
        converged <- TRUE
        if (flag != 1)
          warning("Maximum number of MFVB iterations exceeded in stage ",
                  stage, ".", call. = FALSE)
      }
      parVecHatPrev <- parVecHat
    }
    stage_iters[stage] <- stage_it
  }

  ## ---- sign convention (matches vbfa()) ---------------------------------
  ## Factor signs are arbitrary; flip any factor whose mean loading is
  ## negative and propagate to the scores and the structural coefficients.
  mean_A <- apply(mu.q.A, 2, mean)
  sign_flag <- which(mean_A < 0)
  if (length(sign_flag) > 0) {
    sign_vec <- rep(1, K)
    sign_vec[sign_flag] <- -1
    mu.q.A <- t(t(mu.q.A) * sign_vec)
    mu.F   <- t(t(mu.F) * sign_vec)
    mu.q.B <- mu.q.B * sign_vec        # rows are factors
    Sig    <- Sig * outer(sign_vec, sign_vec)
    U      <- U   * outer(sign_vec, sign_vec)
  }

  endTime <- Sys.time()

  pi_A <- Alpha
  pi_B <- Beta

  out <- list(model = "vbmimic", call = call,
              nobs = N, nitem = J, nfactor = K,
              converged = isTRUE(flag == 1),
              A = mu.q.A, B = mu.q.B,
              pi_A = pi_A, pi_B = pi_B,
              Q_A = Q_A, Q_B = Q_B,
              eta = mu.F, Phi = cov2cor(Sig), Sig = Sig, U = U, V = V,
              rho = RHO, theta = Theta,
              A_var = sigsq.q.A, B_var = sigsq.q.B,
              iter = itNum, flag = flag,
              time = difftime(endTime, startTime, units = "secs"),
              ELBO = NA_real_, objective = NA_real_, objective_type = "none",
              standardize = standardize,
              ## stored only when something is missing; see vbfa()
              preprocess = list(missing_mask = if (has_missing) missing_mask else NULL,
                                n_missing = sum(missing_mask),
                                response_type = "continuous"),
              path = list(v0 = v0_seq, n_stage = n_stage,
                          stage_iters = stage_iters))
  class(out) <- c("vbmimic", "vbpm_fit")
  out
}
