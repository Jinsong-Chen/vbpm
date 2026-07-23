#' @title Simulating data with Latent Variable Modeling
#'
#' @description \code{sim_fa} can simulate data based on factor analysis or
#' item response models with different response formats (continuous or categorical),
#' loading patterns and residual covariance (local dependence) structures.
#'
#' \code{sim_fa()} and [sim_lvm()] are complementary rather than
#' interchangeable: \code{sim_fa()} is the *structure-side* generator, driven
#' by the loading pattern itself -- items per factor (`ipf`), alternating-sign
#' cross-loadings (`alt_sign`), or minor factors (`minor`, definition-based
#' per Auerswald & Moshagen, 2019) -- and is the generator for [vbfa()].
#' [sim_lvm()] is the *model-side* generator, adding observed/latent
#' predictors (`P`/`b`, `K1`/`ph1`/`b1`) and mixed response formats (`ilvl`),
#' and is the generator for [vbmimic()]. Both generators also accept a
#' population loading matrix directly (`mla`) with identical semantics.
#'
#' @name sim_fa
#'
#' @param N Sample size.
#'
#' @param mla Population loading matrix (J x K). If supplied, the pattern
#' arguments (\code{ipf}, \code{cpf}, \code{lam}, \code{lac}, \code{alt_sign},
#' \code{add_*}, \code{minor}, \code{zero_it}) are ignored and \code{K},
#' \code{J} are taken from it -- the same convention as [sim_lvm()].
#'
#' @param K Number of factors.
#'
#' @param ipf Items per factor.
#'
#' @param cpf Cross-loadings per factor.
#'
#' @param lam Number of formal iterations for posterior sampling.
#'
#' @param lac Number of iterations to update the sampling information.
#'
#' @param phi Homogeneous correlations between any two factors.
#'
#' @param ph12 Correlation between factor 1 and 2 (if it's different from
#' \code{phi}). Shared with [sim_lvm()], which uses \code{ph1} for a
#' different purpose (correlation among latent predictors); see [sim_lvm()]
#' for that distinction.
#'
#' @param ecr Residual correlation (local dependence).
#'
#' @param ome_out Output factor score or not.
#'
#' @param cati The set of categorical (polytomous) items in sequence number (i.e., 1 to \eqn{J});
#' \code{NULL} for no and -1 for all (default is \code{NULL}).
#'
#' @param noc Number of categories for categorical items
#'
#' @param misp Proportion of missingness.
#'
#' @param rseed An integer for the random seed.
#'
#' @param necb Number of between-factor local dependence.
#'
#' @param necw Number of within-factor local dependence.
#'
#' @param add_ind (Additional) minor factor with cross-loadings.
#'
#' @param add_la Value of cross-loadings on (Additional) minor factor.
#'
#' @param add_phi Correlations between (Additional) minor factor and other factors.
#'
#' @param alt_sign Logical; alternate the sign of the cross-loadings across the
#' cross-loading items (default \code{TRUE}).
#'
#' @param add1_ind A second (additional) minor factor with cross-loadings,
#' given as item indices; \code{c()} for none.
#'
#' @param add1_la Value of the cross-loadings on the second minor factor.
#'
#' @param add1_phi Correlations between the second minor factor and other
#' factors.
#'
#' @param minor Optional character vector requesting orthogonal minor factors
#' loading on every indicator, each \code{"weak"} or \code{"moderate"}
#' (Auerswald & Moshagen, 2019); \code{NULL} for none.
#'
#' @param zero_it Surplus items with zero loading.
#'
#' @param digits Number of significant digits to print when printing numeric values.
#'
#' @return An object of class \code{list} containing the data, loading, and factorial correlation matrix.
#'
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
#' @examples
#'
#' # for continuous data with cross-loadings and local dependence effect .3
#' out <- sim_fa(N=1000,K=3,ipf=6,lam = .7, lac=.3,ecr=.3)
#' summary(out$dat)
#' out$MLA
#' out$ofd_ind
#'
#' # for categorical data with cross-loadings .4 and 10% missingness
#' out <- sim_fa(N=1000,K=3,ipf=6,lam = .7, lac=.4,cati=-1,noc=4,misp=.1)
#' summary(out$dat)
#' out$MLA
#' out$ofd_ind
#'
#' ## matrix-driven: supply the population loading matrix directly
#' mla <- matrix(0, 12, 2)
#' mla[1:6, 1] <- .7; mla[7:12, 2] <- c(.8, .7, .6, .7, .8, .6)
#' s <- sim_fa(N = 300, mla = mla, phi = .3, rseed = 4)
#' s$MLA                      # the generating matrix, as supplied
#'
#' @importFrom MASS mvrnorm
#' @export
sim_fa <- function(N = 1000, mla = NULL, K = 3, ipf = 8, cpf = 2, lam = 0.7, lac = 0.3, alt_sign = TRUE, phi = 0.5, ph12 = -1,
        ecr = .0, ome_out = FALSE,cati = NULL, noc = c(4), misp = 0, rseed = 333,
        necw=K,necb=K,add_ind=c(),add_la=0,add_phi=0,
        add1_ind=c(), add1_la=0, add1_phi=0, zero_it=0, minor=NULL, digits = 4) {

    ## `minor`: an optional character vector requesting orthogonal minor factors,
    ## each loading on EVERY indicator, following Auerswald & Moshagen (2019,
    ## Psychological Methods) / Lorenzo-Seva et al. (2011):
    ##   "weak"     -> loadings ~ U(-.10, .10)            (~0.33% explained var)
    ##   "moderate" -> loadings ~ +/- U(.09, .11)          (~1% explained var)
    ## e.g. minor = c("weak","weak","moderate","moderate") adds 2 weak + 2 moderate
    ## minor factors. Minor factors are uncorrelated with the major factors and
    ## with one another. This is the modern, definition-based alternative to the
    ## ad hoc add_ind/add1_ind mechanism above.

    if (exists(".Random.seed", .GlobalEnv))
        oldseed <- .GlobalEnv$.Random.seed else oldseed <- NULL
        set.seed(rseed)

    oo <- options()       # code line i
    on.exit(options(oo))  # code line i+1
    # old_digits <- getOption("digits")
    options(digits = digits)

    use_mla <- !is.null(mla)
    if (use_mla) {
        ## ---- matrix-driven: supplied mla IS the population loading matrix.
        ## Pattern arguments (ipf, cpf, lam, lac, alt_sign, add_*, minor,
        ## zero_it) are ignored; K, J, ipf, cpf follow sim_lvm's convention.
        ## Determined here (before PHI is built) so PHI comes out K x K.
        mla <- as.matrix(mla)
        K   <- ncol(mla)
        J   <- nrow(mla)
        K1  <- K
        ipf <- J / K
        cpf <- 0
        if (ecr > 0 && J %% K != 0)
            stop("Local dependence (ecr > 0) with a supplied mla assumes J is divisible by K; use ecr = 0 with an unbalanced mla.", call. = FALSE)
    }

    PHI <- matrix(phi, K, K)
    diag(PHI) <- 1
    if (ph12 > -1 && ph12 < 1) PHI[1,2] <- PHI[2,1]<-ph12

    if (!use_mla) {
    if (alt_sign) {
      cross_vals <- rep(c(lac, -lac), length.out = cpf)
    } else {
      cross_vals <- rep(lac, cpf)
    }
    lam0 <- c(rep(lam, ipf), rep(0, ipf - cpf), cross_vals, rep(0, (K - 2) * ipf))
    # lam0 <- c(rep(lam, ipf), rep(0, ipf - cpf), rep(lac, cpf), rep(0, (K - 2) * ipf))

    lam1 <- matrix(lam0, ipf, K)

    mla <- c()
    for (i in K:1) {
        # i <- 1
        ind <- (c(1:K) + i)%%K
        ind[ind == 0] <- K
        mla <- rbind(mla, lam1[, ind])
    }

    J <- K * ipf
    K1<-K
    if(length(add_ind)>0){
        # len<-length(add_ind)
        K1<-K+1
        mla<-cbind(mla,0)
        mla[add_ind,K1]<-add_la
        PHI<-cbind(PHI,add_phi)
        PHI<-rbind(PHI,add_phi)
        PHI[K1,K1]<-1
    }

    if(length(add1_ind)>0){
      # len<-length(add_ind)
      K1<-K1+1
      mla<-cbind(mla,0)
      mla[add1_ind,K1]<-add1_la
      PHI<-cbind(PHI,add1_phi)
      PHI<-rbind(PHI,add1_phi)
      PHI[K1,K1]<-1
    }

    ## modern minor-factor mechanism (definition-based; see signature comment)
    if (length(minor) > 0) {
        Jm <- nrow(mla)                       # current number of items (= K*ipf)
        for (mt in minor) {
            load_col <- switch(mt,
                weak     = runif(Jm, -0.10, 0.10),
                moderate = sample(c(-1, 1), Jm, replace = TRUE) * runif(Jm, 0.09, 0.11),
                stop("minor must contain only 'weak' or 'moderate'"))
            K1  <- K1 + 1
            mla <- cbind(mla, load_col)       # minor factor loads on every item
            PHI <- cbind(PHI, 0); PHI <- rbind(PHI, 0); PHI[K1, K1] <- 1  # orthogonal
        }
    }

    if (zero_it>0){
        J <- J + zero_it
        mla <- rbind(mla,matrix(0,zero_it,K1))
    }
    }

    ecm <- matrix(0, J, J)
    if (ecr > 0) {
        iecb <- iecw <- NULL
        li <- c(1:K) * ipf - cpf
        for (i in 1:necw) {
            if(necw>0) ecm[li[i], li[i] - 1] <- ecm[li[i] - 1,li[i]] <- ecr
        }
        for (i in 1:necb) {
            i1 <- i%%K + 1
            if(necb>0) ecm[li[i] - 2, li[i1] - 3] <- ecm[li[i1] - 3,li[i] - 2]<- ecr
        }
    }

    # iecb <- iecw <- NULL
    # li <- c(1:K) * ipf - cpf
    # for (i in 1:necw) {
    #     if(necw>0) iecw <- rbind(iecw, c(li[i], li[i] - 1))
    #     # i1 <- i%%K + 1
    #     # iecb <- rbind(iecb, c(li[i] - 2, li[i1] - 3))
    # }
    # for (i in 1:necb) {
    #     # iecw <- rbind(iecw, c(li[i], li[i] - 1))
    #     i1 <- i%%K + 1
    #     if(necb>0) iecb <- rbind(iecb, c(li[i] - 2, li[i1] - 3))
    # }

    # ecm <- matrix(0, J, J)
    #
    # if (ecr > 0) {
    #     ecm[iecw] <- ecm[iecw[, 2:1]] <- ecr
    #     ecm[iecb] <- ecm[iecb[, 2:1]] <- ecr
    # }

    evr <- rep(0, J)
    for (j in 1:J) {
        evr[j] = 1 - t(mla[j, ]) %*% PHI %*% mla[j, ]
    }  #end j

    diag(ecm) <- evr
    er <- mvrnorm(N, rep(0, J), ecm)
    eta <- mvrnorm(N, rep(0, K1), PHI)
    y <- t(mla %*% t(eta)) + er
    scale <- apply(y, 2, sd)
    eigen <- diag(crossprod(mla))

    out <- list(N = N, PHI = PHI, MLA = mla, Eigen = eigen, PSX = ecm, scale = scale)
    if (ome_out)
        out$OME <- eta

    pos <- lower.tri(ecm)
    ind <- which(pos, arr.ind = T)
    rind <- which(ecm[pos] > 0)
    out$ofd_ind <- ind[rind, ]

    Jp <- length(cati)
    if (Jp > 0)
        {
            if (Jp == 1 && cati == -1) {
                cati <- c(1:J)
                Jp <- J
            }
            yc <- y[, cati]
            len <- length(noc)
            if (Jp%%len != 0)
                stop("cati is not divisable by length(noc).", call. = F)
            for (i in 1:len) {
                M <- noc[i]  #categories
                UL <- 3
                LL <- -3
                stp <- LL + c(1:(M - 1)) * (UL - LL)/M  #M-1 step points
                i0 <- (i - 1) * Jp/len + 1
                i1 <- i * Jp/len
                for (j in i0:i1) {
                  tmp <- yc[, j] > matrix(stp, N, (M - 1), byrow = T)
                  yc[, j] <- rowSums(tmp) + 1  # value starting from 1
                }
            }  #end len
            y[, cati] <- yc
            out$noc <- noc
            out$cati <- cati
        }  #end Jp

    mind <- matrix(rbinom(N * J, 1, misp), N, J)
    y[mind == 1] <- NA
    out$misp <- misp
    out$dat <- y

    if (!is.null(oldseed))
        .GlobalEnv$.Random.seed <- oldseed else rm(".Random.seed", envir = .GlobalEnv)

    return(out)
}
