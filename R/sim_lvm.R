## Latent-variable-model data generator, absorbed from LAWBL::sim_lvm (Chen,
## CRAN 1.5.0) as part of retiring LAWBL into vbpm. Changes on absorption:
##   * `phd` added -- correlated structural disturbances, so the MIMIC generator
##     of Jin & Chen (2025) is covered (its factor residual covariance U is
##     unrestricted; LAWBL's was always diagonal).
##   * missingness bug fixed -- the N x J mask was applied to an N x (J + P)
##     matrix, so R recycled it and NAs bled into the covariate block.
##   * `P > 0` no longer forces the covariates to share `phi` with the factors;
##     `phx` controls them separately (default = `phi` for back-compatibility).
##   * input validation, namespace imports, no partial-argument reliance.
## The generating equations are otherwise unchanged.

#' Simulate data from a latent variable model
#'
#' Generates factor-analytic or item-response data with optional observed and
#' latent predictors (MIMIC / structural models), cross-loadings, local
#' dependence, categorical or mixed response formats, and missingness.
#'
#' `sim_lvm()` and [sim_fa()] are complementary rather than interchangeable:
#'
#' \describe{
#'   \item{`sim_lvm()`}{is the *model-side* generator. Use it when the data
#'     structure is driven by a supplied loading matrix (`mla`), or when the
#'     design involves **predictors** -- observed covariates (`P`, `b`) or
#'     latent predictors (`K1`, `ph1`, `b1`) -- or **mixed response formats**
#'     per item (`ilvl`). This is the generator for [vbmimic()].}
#'   \item{[sim_fa()]}{is the *structure-side* generator. Use it when the
#'     question is about the loading pattern itself: items per factor (`ipf`),
#'     alternating-sign cross-loadings (`alt_sign`), or minor factors
#'     (`minor`, definition-based per Auerswald & Moshagen, 2019). This is the
#'     generator for [vbfa()].}
#' }
#'
#' Note one deliberate difference in argument naming, inherited from the two
#' functions' separate histories: `sim_lvm()` uses `ph12` for the factor 1-2
#' correlation and `ph1` for correlations **among latent predictors**, whereas
#' [sim_fa()] uses `ph1` for the factor 1-2 correlation. They are not
#' interchangeable; check which function you are calling.
#'
#' @param N Sample size.
#' @param mla Population loading matrix. If supplied, `K`, `J` and `cpf` are
#'   taken from it.
#' @param K Number of factors (used when `mla = NULL`). Must divide `J`.
#' @param J Number of items (used when `mla = NULL`).
#' @param cpf Number of cross-loadings per factor (used when `mla = NULL`).
#' @param lam Value of the primary loadings.
#' @param lac Value of the cross-loadings.
#' @param phi Homogeneous correlation between any two factors.
#' @param ph12 Correlation between factors 1 and 2, when it differs from `phi`;
#'   `-1` (default) to use `phi`.
#' @param phx Correlation among the observed covariates. Defaults to `phi`,
#'   which is LAWBL's behaviour; set `phx = 0` for uncorrelated covariates.
#' @param ecr Residual correlation (local dependence).
#' @param P Number of observed predictors (MIMIC model); `0` for none.
#' @param b Coefficients of the observed predictors: a scalar, or a
#'   `(K - K1) x P` matrix.
#' @param K1 Number of latent predictors (MIMIC / structural model).
#' @param ph1 Correlation among the latent predictors.
#' @param b1 Coefficients of the latent predictors: a scalar, or a
#'   `(K - K1) x K1` matrix.
#' @param phd Correlation among the structural disturbances (the residuals of
#'   the regression of factors on predictors). `0` (default) gives independent
#'   disturbances, which is LAWBL's behaviour; a positive value reproduces
#'   designs with a freely correlated factor residual covariance. Ignored when
#'   `P = 0` and `K1 = 0`.
#' @param ilvl Levels of every item, length `J + P`; values below 2 mean
#'   continuous. `NULL` to use `cati`/`noc` instead.
#' @param cati Polytomous items by index; `NULL` for none, `-1` for all. Used
#'   only when `ilvl = NULL`.
#' @param noc Number of categories for the polytomous items.
#' @param misp Proportion of missingness, applied to the `J` items only (never
#'   to the covariate block).
#' @param ome_out Return the factor scores.
#' @param necw Number of within-factor local dependencies.
#' @param necb Number of between-factor local dependencies.
#' @param add_ind Indices of items loading on an additional minor factor.
#' @param add_la Loading value on the additional minor factor.
#' @param add_phi Correlation between the additional minor factor and the
#'   others.
#' @param zero_it Number of surplus items with zero loadings.
#' @param rseed Random seed.
#' @param digits Significant digits used while the function is running.
#'
#' @return A list with the data (`dat`; items first, then any covariates),
#'   the loading matrix (`MLA`), factor correlations (`PHI`), residual
#'   covariance (`PSX`), the indices of the local-dependence pairs
#'   (`ofd_ind`), and, when predictors are present, the coefficient matrices
#'   (`mb`, `mb1`) and covariate correlation (`PHX`).
#'
#' @references
#' Chen, J. (2022). *LAWBL: Latent (variable) analysis with Bayesian
#' learning* (R package version 1.5.0).
#'
#' Jin, Y., & Chen, J. (2025). Regularized variational Bayesian approximations
#' for variable selection in extended multiple-indicators multiple-causes
#' models. *Multivariate Behavioral Research*.
#' \doi{10.1080/00273171.2025.2483253}
#'
#' @seealso [sim_fa()] for the structure-side generator, [vbmimic()] and
#'   [vbfa()] for the estimators.
#'
#' @examples
#' ## MIMIC data: 18 items on 3 factors, 9 covariates in disjoint blocks of 3.
#' ## Supply `b` as a matrix for a sparse structural design; a scalar `b` makes
#' ## every covariate predict every factor, which quickly leaves no disturbance
#' ## variance (the function stops if it does).
#' B <- matrix(0, 3, 9)
#' for (k in 1:3) B[k, ((k - 1) * 3 + 1):(k * 3)] <- .3
#'
#' s <- sim_lvm(N = 500, K = 3, J = 18, P = 9, b = B, phx = 0, rseed = 1)
#' dim(s$dat)      # 500 x 27 -- items first, covariates last
#' s$mb            # the generating structural coefficients
#'
#' ## correlated disturbances (unrestricted factor residual covariance)
#' s2 <- sim_lvm(N = 500, K = 3, J = 18, P = 9, b = B, phx = 0, phd = .3,
#'               rseed = 1)
#' round(s2$PHI, 2)
#'
#' @importFrom MASS mvrnorm
#' @export
sim_lvm <- function(N = 1000, mla = NULL, K = 3, J = 18, cpf = 0,
                    lam = 0.7, lac = 0.3, phi = 0.3, ph12 = -1, phx = NULL,
                    ecr = 0, P = 0, b = 0.3, K1 = 0, ph1 = 0.2, b1 = 0.3,
                    phd = 0, ilvl = NULL, cati = NULL, noc = c(4), misp = 0,
                    ome_out = FALSE, necw = K, necb = K, add_ind = c(),
                    add_la = 0.5, add_phi = 0, zero_it = 0, rseed = 333,
                    digits = 4) {

  if (is.null(mla) && J %% K != 0)
    stop("J should be a multiple of K.", call. = FALSE)
  if (K1 > 0 && K1 >= K)
    stop("K1 (latent predictors) must be smaller than K.", call. = FALSE)
  if (phd < 0 || phd >= 1)
    stop("phd must be in [0, 1).", call. = FALSE)
  if (is.null(phx)) phx <- phi

  if (exists(".Random.seed", .GlobalEnv))
    oldseed <- .GlobalEnv$.Random.seed
  else oldseed <- NULL
  set.seed(rseed)

  oo <- options()
  on.exit(options(oo))
  options(digits = digits)

  if (!is.null(mla)) {
    mla <- as.matrix(mla)
    K   <- ncol(mla)
    J   <- nrow(mla)
    cpf <- 0
  }

  x <- NULL

  if (P == 0 && K1 == 0) {
    ## ---- plain factor model ------------------------------------------
    PHI <- matrix(phi, K, K)
    diag(PHI) <- 1
    if (ph12 > -1 && ph12 < 1) PHI[1, 2] <- PHI[2, 1] <- ph12
    PHX <- NULL
    eta <- mvrnorm(N, rep(0, K), PHI)
  } else {
    ## ---- MIMIC / structural model -------------------------------------
    mb <- mb1 <- PHX <- NULL
    tmp <- eta <- 0
    Kd  <- K - K1                 # factors that are regressed on predictors

    if (P > 0) {
      PHX <- matrix(phx, P, P)
      diag(PHX) <- 1
      x <- mvrnorm(N, rep(0, P), PHX)
      mb <- if (length(b) == 1) matrix(b, Kd, P) else as.matrix(b)
      if (nrow(mb) != Kd || ncol(mb) != P)
        stop("b must be a scalar or a (K - K1) x P matrix.", call. = FALSE)
      tmp <- mb %*% PHX %*% t(mb)
      eta <- t(mb %*% t(x))
    }

    if (K1 > 0) {
      PH1 <- matrix(ph1, K1, K1)
      diag(PH1) <- 1
      ome <- mvrnorm(N, rep(0, K1), PH1)
      mb1 <- if (length(b1) == 1) matrix(b1, Kd, K1) else as.matrix(b1)
      if (nrow(mb1) != Kd || ncol(mb1) != K1)
        stop("b1 must be a scalar or a (K - K1) x K1 matrix.", call. = FALSE)
      tmp <- tmp + mb1 %*% PH1 %*% t(mb1)
      eta <- eta + t(mb1 %*% t(ome))
    }

    ## structural disturbances: variances make each factor unit-variance;
    ## phd > 0 correlates them (an unrestricted factor residual covariance).
    dvar <- 1 - diag(tmp)
    if (any(dvar <= 0))
      stop("Predictor coefficients imply non-positive disturbance variance; ",
           "reduce b/b1.", call. = FALSE)
    if (Kd > 1) {
      DIS <- diag(dvar)
      if (phd > 0) {
        sdv <- sqrt(dvar)
        DIS <- phd * outer(sdv, sdv)
        diag(DIS) <- dvar
      }
    } else {
      DIS <- dvar
    }
    eta <- eta + mvrnorm(N, rep(0, Kd), DIS)
    if (K1 > 0) eta <- cbind(eta, ome)
  }

  ## ---- loading matrix ---------------------------------------------------
  ipf <- J / K
  if (is.null(mla)) {
    if (K > 1) {
      lam0 <- c(rep(lam, ipf), rep(0, ipf - cpf), rep(lac, cpf),
                rep(0, (K - 2) * ipf))
    } else {
      lam0 <- lam
    }
    lam1 <- matrix(lam0, ipf, K)
    if (K == 1) {
      mla <- lam1
    } else {
      mla <- c()
      for (i in K:1) {
        ind <- (c(1:K) + i) %% K
        ind[ind == 0] <- K
        mla <- rbind(mla, lam1[, ind])
      }
    }
    if (P == 0 && K1 == 0) {
      Kp1 <- K
      if (!is.null(add_ind)) {
        Kp1 <- K + 1
        mla <- cbind(mla, 0)
        mla[add_ind, Kp1] <- add_la
        PHI <- cbind(PHI, add_phi)
        PHI <- rbind(PHI, add_phi)
        PHI[Kp1, Kp1] <- 1
      }
      if (zero_it > 0) {
        J <- J + zero_it
        mla <- rbind(mla, matrix(0, zero_it, Kp1))
      }
      eta <- mvrnorm(N, rep(0, Kp1), PHI)
    }
  }

  ## ---- residual covariance (local dependence) ---------------------------
  ecm <- matrix(0, J, J)
  if (ecr > 0) {
    li <- c(1:K) * ipf - cpf
    for (i in 1:necw) {
      if (necw > 0) ecm[li[i], li[i] - 1] <- ecm[li[i] - 1, li[i]] <- ecr
    }
    for (i in 1:necb) {
      i1 <- i %% K + 1
      if (necb > 0)
        ecm[li[i] - 2, li[i1] - 3] <- ecm[li[i1] - 3, li[i] - 2] <- ecr
    }
  }

  PHI <- cor(eta)
  evr <- 1 - diag(mla %*% PHI %*% t(mla))
  diag(ecm) <- evr
  er <- mvrnorm(N, rep(0, J), ecm)
  y  <- t(mla %*% t(eta)) + er
  if (P > 0) y <- cbind(y, x)

  scale <- apply(y, 2, sd)
  eigen <- diag(crossprod(mla))
  out <- list(N = N, PHI = PHI, MLA = mla, Eigen = eigen, PSX = ecm,
              scale = scale)
  if (ome_out) out$OME <- eta

  pos  <- lower.tri(ecm)
  ind  <- which(pos, arr.ind = TRUE)
  rind <- which(ecm[pos] > 0)
  out$ofd_ind <- ind[rind, ]

  ## ---- categorization ----------------------------------------------------
  UL <- 3
  LL <- -3
  if (is.null(ilvl)) {
    Jp <- length(cati)
    if (Jp > 0) {
      if (Jp == 1 && cati == -1) {
        cati <- c(1:J)
        Jp <- J
      }
      yc  <- y[, cati]
      len <- length(noc)
      if (Jp %% len != 0)
        stop("cati is not divisable by length(noc).", call. = FALSE)
      for (i in 1:len) {
        M   <- noc[i]
        stp <- LL + c(1:(M - 1)) * (UL - LL) / M
        i0  <- (i - 1) * Jp / len + 1
        i1  <- i * Jp / len
        for (j in i0:i1) {
          tmp2 <- yc[, j] > matrix(stp, N, (M - 1), byrow = TRUE)
          yc[, j] <- rowSums(tmp2) + 1
        }
      }
      y[, cati] <- yc
      out$noc  <- noc
      out$cati <- cati
    }
  } else {
    if (length(ilvl) != J + P)
      stop("number of item levels should be J + P.", call. = FALSE)
    for (i in 1:(J + P)) {
      M <- ilvl[i]
      if (M >= 2) {
        M    <- round(M)
        stp  <- LL + c(1:(M - 1)) * (UL - LL) / M
        tmp2 <- y[, i] > matrix(stp, N, (M - 1), byrow = TRUE)
        y[, i] <- rowSums(tmp2) + 1
      }
    }
  }

  ## ---- missingness -------------------------------------------------------
  ## applied to the J items only. LAWBL built an N x J mask and indexed an
  ## N x (J + P) matrix with it, so the logical vector was recycled and NAs
  ## appeared in the covariate block; the explicit column restriction below
  ## is the fix.
  if (misp > 0) {
    mind <- matrix(rbinom(N * J, 1, misp), N, J)
    yi   <- y[, 1:J, drop = FALSE]
    yi[mind == 1] <- NA
    y[, 1:J] <- yi
  }
  out$misp <- misp
  out$dat  <- y

  if (P != 0 || K1 != 0) {
    out$mb  <- mb
    out$mb1 <- mb1
    out$PHX <- PHX
  }

  if (!is.null(oldseed))
    .GlobalEnv$.Random.seed <- oldseed
  else rm(".Random.seed", envir = .GlobalEnv)

  out
}
