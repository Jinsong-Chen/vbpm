## ---------------------------------------------------------------------------
## data-raw/make_r8_fixture.R
##
## Rebuilds inst/extdata/vbpm_r8_fixture.rds, the R8 integration fixture of
## plan section 6.1.  Runnable standalone from the package root (or from
## data-raw/):
##
##   Rscript data-raw/make_r8_fixture.R
##
## The file holds two things: `$members`, a named list of low-level SWEEPS,
## and `$integrity`, the one complete sealed `pefa` object of section 8.1.
## Every loading and PIP matrix in either is constructed by hand or from a
## seeded rnorm jitter -- no fitted number is ever stored, so nothing frozen
## here depends on a BLAS.  Each member keeps its matrices together with the
## FROZEN expected $transitions, $persistence, and pivoted persistence()
## matrices, so tests/testthat/test-r8-fixture.R re-runs only the matcher and
## the two table builders.
##
## Two independent checks run before anything is frozen, because a fixture that
## freezes a bug is worse than no fixture:
##
##   1. `check_facts()` compares every discrete pair fact of every edge against
##      a hand-derived table written out below with the member that produced
##      it.
##   2. `check_reference()` recomputes all seventeen matcher fields with the
##      naive reference implementation in this file -- direct formulas, no
##      rescaling, explicit loops -- and requires agreement.
##
## The integrity object additionally has to pass the installed verifier both
## ways before it is written; see the long note above its section.
##
## data-raw is .Rbuildignore'd; nothing here ships in the package.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(vbpm))

FIXTURE_VERSION <- "vbpm_r8_fixture_1"
EVIDENCE_SCHEMA_ID <- "vbpm_pefa_evidence_1"

METRICS <- c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
             "unmatched_ssl", "n_collisions")

DOUBLE_PAIR_COLS <- c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                      "unmatched_ssl")


# ---- output location --------------------------------------------------

find_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (p in c(here, dirname(here))) {
    if (file.exists(file.path(p, "DESCRIPTION")) &&
        dir.exists(file.path(p, "R"))) {
      return(p)
    }
  }
  stop("Run this script from the vbpm package root or from data-raw/.",
       call. = FALSE)
}


# ---- canonical shapes -------------------------------------------------

J <- 8L
ITEMS <- paste0("item_", seq_len(J))

factor_names <- function(K0, K, bifactor) {
  nm <- if (K0 > 0L) paste0("F", seq_len(K0)) else character(0)
  if (K > K0) nm <- c(nm, paste0("F", seq.int(K0 + 1L, K)))
  if (isTRUE(bifactor)) nm <- c("G", nm)
  nm
}

## Assemble one stored loading matrix from its columns.  In bifactor mode the
## first column supplied is the general factor, exactly as vbfa() stores it.
lam_of <- function(cols, K0, bifactor = FALSE) {
  m <- do.call(cbind, cols)
  storage.mode(m) <- "double"
  K <- ncol(m) - as.integer(isTRUE(bifactor))
  dimnames(m) <- list(ITEMS, factor_names(K0, K, bifactor))
  m
}

## Deterministic jitter: set.seed() inside, so the result never depends on the
## order in which the members are built.
jitter_lam <- function(m, seed, sd = 0.01) {
  set.seed(seed)
  out <- m + matrix(stats::rnorm(length(m), sd = sd), nrow(m), ncol(m))
  dimnames(out) <- dimnames(m)
  out
}

## Full-width design for one candidate: the backbone, the appended
## regularized exploratory columns, and (bifactor) a leading general column.
design_of <- function(Q0, K0, K, bifactor) {
  Q <- if (K > K0) cbind(Q0, matrix(-1L, nrow(Q0), K - K0)) else Q0
  if (isTRUE(bifactor)) Q <- cbind(matrix(-1L, nrow(Q0), 1L), Q)
  matrix(as.integer(Q), nrow(Q0), K + as.integer(isTRUE(bifactor)))
}

## A plausible PIP payload: monotone in |loading|, inside [0, 1], with the
## Q0-fixed cells pinned so that the matcher's "regularized at both endpoints"
## rule has something real to exclude.
pip_of <- function(lam, Q) {
  a <- abs(lam)
  p <- 0.02 + 0.94 * (a / max(a))
  p[Q == 1L] <- 1
  p[Q == 0L] <- 0
  p <- pmin(1, pmax(0, p))
  dim(p) <- dim(lam)
  dimnames(p) <- dimnames(lam)
  p
}

## Independent restatement of section 3.1: the smallest SSL among the K0
## position-matched backbone columns, which in bifactor mode are the stored
## columns 2:(K0 + 1).
backbone_ssl_min_ref <- function(lam, K0, bifactor, loading_ok) {
  if (K0 == 0L || !isTRUE(loading_ok)) return(NA_real_)
  cols <- seq_len(K0) + as.integer(isTRUE(bifactor))
  min(vapply(cols, function(j) sum(lam[, j]^2), numeric(1)))
}


# ---- synthetic candidate table ----------------------------------------

## Criterion paths with diminishing gains on every retained path, so that
## $transitions has something to normalize and select_K_elbow() has a shape.
crit_paths <- function(n) {
  take <- function(v) v[seq_len(n)]
  ELBO <- take(c(-4000, -3700, -3610, -3585, -3578))
  AIC <- take(c(8100, 7560, 7420, 7390, 7385))
  BIC <- AIC + take(c(60, 84, 108, 132, 156))
  list(
    ELBO = as.double(ELBO),
    t_nom = take(c(20, 28, 36, 44, 52)),
    t_S = take(c(18.5, 26.0, 33.5, 41.0, 48.5)),
    AIC = as.double(AIC), BIC = as.double(BIC),
    AIC_S = AIC - take(c(4, 5, 6, 7, 8)),
    BIC_S = BIC - take(c(12, 15, 18, 21, 24)),
    RMSEA = take(c(.085, .062, .048, .041, .039)),
    SRMR = take(c(.070, .055, .043, .038, .036)),
    CFI = take(c(.900, .941, .968, .975, .977)),
    TLI = take(c(.880, .928, .958, .966, .969)),
    iter = take(c(140L, 210L, 265L, 320L, 380L))
  )
}

build_sweep <- function(Ks, loadings, K0, bifactor,
                        converged = NULL, loading_reason = NULL,
                        pip_reason = NULL, elbo_reason = NULL,
                        stats_reason = NULL) {
  n <- length(Ks)
  blank <- rep(NA_character_, n)
  if (is.null(converged)) converged <- rep(TRUE, n)
  if (is.null(loading_reason)) loading_reason <- blank
  if (is.null(pip_reason)) pip_reason <- blank
  if (is.null(elbo_reason)) elbo_reason <- blank
  if (is.null(stats_reason)) stats_reason <- blank

  cp <- crit_paths(n)
  st <- function(reason) ifelse(is.na(reason), "ok", "unavailable")
  ## Section 3.1.1: the fit component is "ok" exactly when the recognized
  ## convergence flag normalizes to TRUE; a recognized non-converging flag is
  ## `nonconverged`, an unrecognizable one is `malformed_fit`.
  fit_status <- ifelse(converged %in% TRUE, "ok", "unavailable")
  fit_reason <- ifelse(converged %in% TRUE, NA_character_,
                       ifelse(is.na(converged), "malformed_fit",
                              "nonconverged"))
  ELBO <- cp$ELBO
  ELBO[!is.na(elbo_reason)] <- NA_real_
  loading_ok <- is.na(loading_reason)

  data.frame(
    K = as.integer(Ks), ELBO = ELBO,
    t_nom = as.double(cp$t_nom), t = as.double(cp$t_nom),
    t_S = as.double(cp$t_S),
    AIC = cp$AIC, BIC = cp$BIC,
    AIC_S = as.double(cp$AIC_S), BIC_S = as.double(cp$BIC_S),
    RMSEA = as.double(cp$RMSEA), SRMR = as.double(cp$SRMR),
    CFI = as.double(cp$CFI), TLI = as.double(cp$TLI),
    iter = as.integer(cp$iter), secs = rep(0, n),
    converged = converged,
    backbone_ssl_min = vapply(
      seq_len(n),
      function(i) backbone_ssl_min_ref(loadings[[i]], K0, bifactor,
                                       loading_ok[i]),
      numeric(1)
    ),
    fit_status = fit_status, fit_reason = fit_reason,
    loading_status = st(loading_reason), loading_reason = loading_reason,
    pip_status = st(pip_reason), pip_reason = pip_reason,
    elbo_status = st(elbo_reason), elbo_reason = elbo_reason,
    stats_status = st(stats_reason), stats_reason = stats_reason,
    stringsAsFactors = FALSE
  )
}


# ---- naive reference matcher (second implementation) ------------------

## Section 2, written the slow obvious way: no rescaling before the dot
## product, explicit loops for the pooled RMSD and the PIP discrepancy, and a
## contingency table built from factors rather than from max.col().  Its only
## job is to disagree with R/pefa_match.R if either of them is wrong.
ref_pair_facts <- function(lam_from, lam_to, pip_from, pip_to,
                           reg_from, reg_to, K0, eps,
                           from_ok, to_ok, pip_ok_from, pip_ok_to) {
  Jn <- nrow(lam_from)
  Kf <- ncol(lam_from)
  Kt <- ncol(lam_to)
  both <- from_ok && to_ok

  nrm <- function(v) sqrt(sum(v * v))
  cg <- function(a, b) {
    if (nrm(a) == 0 || nrm(b) == 0) return(NA_real_)
    min(1, abs(sum(a * b)) / (nrm(a) * nrm(b)))
  }
  sg <- function(a, b) if (sum(a * b) >= 0) 1 else -1
  ar <- function(a, b, s) sqrt(sum((a - s * b)^2) / Jn)
  ssl <- function(m, cols) vapply(cols, function(j) sum(m[, j]^2), numeric(1))

  bb <- seq_len(K0)
  sf <- if (Kf > K0) seq.int(K0 + 1L, Kf) else integer(0)
  tf <- if (Kt > K0) seq.int(K0 + 1L, Kt) else integer(0)

  es <- if (from_ok) sf[ssl(lam_from, sf) >= eps] else integer(0)
  et <- if (to_ok) tf[ssl(lam_to, tf) >= eps] else integer(0)
  source_ineligible_n <- if (from_ok) length(sf) - length(es) else NA_integer_
  target_eligible_n <- if (to_ok) length(et) else NA_integer_

  zero_backbone <- function(m) {
    if (K0 == 0L) return(FALSE)
    s <- ssl(m, bb)
    any(is.finite(s) & s == 0)
  }
  backbone_degenerate <- if (K0 == 0L) {
    FALSE
  } else if ((from_ok && zero_backbone(lam_from)) ||
             (to_ok && zero_backbone(lam_to))) {
    TRUE
  } else if (both) {
    FALSE
  } else {
    NA
  }

  asg <- rep(NA_integer_, Kf)
  if (both && length(es) && length(et)) {
    for (i in es) {
      a <- lam_from[, i]
      ph <- vapply(et, function(j) cg(a, lam_to[, j]), numeric(1))
      usable <- which(!is.na(ph))
      if (!length(usable)) next
      best <- usable[ph[usable] == max(ph[usable])]
      if (length(best) > 1L) {
        keys <- vapply(best, function(p) {
          b <- lam_to[, et[p]]
          ar(a, b, sg(a, b))
        }, numeric(1))
        if (!anyNA(keys)) best <- best[keys == min(keys)]
      }
      asg[i] <- et[best[1L]]
    }
  }

  ps <- if (both) c(bb, which(!is.na(asg))) else integer(0)
  pt <- if (both) c(bb, asg[!is.na(asg)]) else integer(0)
  m <- length(ps)
  n_pairs <- if (both) m else NA_integer_

  phi_min <- NA_real_
  rmsd <- NA_real_
  rmsd_max <- NA_real_
  if (m > 0L) {
    phis <- numeric(m)
    rms <- numeric(m)
    for (p in seq_len(m)) {
      a <- lam_from[, ps[p]]
      b <- lam_to[, pt[p]]
      phis[p] <- cg(a, b)
      rms[p] <- ar(a, b, sg(a, b))
    }
    if (!anyNA(phis)) phi_min <- min(phis)
    if (!anyNA(rms)) {
      rmsd_max <- max(rms)
      total <- 0
      for (p in seq_len(m)) {
        a <- lam_from[, ps[p]]
        b <- lam_to[, pt[p]]
        total <- total + sum((a - sg(a, b) * b)^2)
      }
      rmsd <- sqrt(total / (Jn * m))
    }
  }

  collision <- NA
  n_collisions <- NA_integer_
  collision_targets <- NA_character_
  collision_multiplicities <- NA_character_
  unmatched_n <- NA_integer_
  unmatched_ssl <- NA_real_
  if (both) {
    used <- asg[!is.na(asg)]
    counts <- table(used)
    reused <- sort(as.integer(names(counts)[counts > 1L]))
    n_collisions <- length(reused)
    collision <- n_collisions > 0L
    collision_targets <- paste(reused, collapse = ",")
    collision_multiplicities <- paste(
      vapply(reused, function(j) sum(used == j), integer(1)), collapse = ","
    )
    unm <- tf[!(tf %in% used)]
    unmatched_n <- length(unm)
    unmatched_ssl <- if (unmatched_n > 0L) max(ssl(lam_to, unm)) else 0
  }

  ari <- NA_real_
  if (both) {
    pcf <- c(bb, es)
    pct <- c(bb, et)
    if (length(pcf) && length(pct) && Jn >= 2L) {
      lab <- function(m2, cols) {
        vapply(seq_len(Jn),
               function(i) which.max(abs(m2[i, cols])), integer(1))
      }
      p <- factor(lab(lam_from, pcf), levels = seq_along(pcf))
      q <- factor(lab(lam_to, pct), levels = seq_along(pct))
      tb <- table(p, q)
      c2 <- function(x) x * (x - 1) / 2
      A <- sum(c2(tb))
      B <- sum(c2(rowSums(tb)))
      Cc <- sum(c2(colSums(tb)))
      E <- B * Cc / c2(Jn)
      D <- (B + Cc) / 2 - E
      ari <- if (D == 0) 1 else min(1, max(-1, (A - E) / D))
    }
  }

  pip_rmsd <- NA_real_
  if (both && isTRUE(pip_ok_from) && isTRUE(pip_ok_to) && m > 0L) {
    total <- 0
    cnt <- 0L
    for (p in seq_len(m)) {
      for (i in seq_len(Jn)) {
        if (reg_from[i, ps[p]] && reg_to[i, pt[p]]) {
          total <- total + (pip_from[i, ps[p]] - pip_to[i, pt[p]])^2
          cnt <- cnt + 1L
        }
      }
    }
    if (cnt > 0L) pip_rmsd <- sqrt(total / cnt)
  }

  list(
    n_pairs = n_pairs, phi_min = phi_min, rmsd = rmsd, rmsd_max = rmsd_max,
    ari = ari, pip_rmsd = pip_rmsd, unmatched_n = unmatched_n,
    unmatched_ssl = unmatched_ssl, collision = collision,
    n_collisions = n_collisions, collision_targets = collision_targets,
    collision_multiplicities = collision_multiplicities,
    source_ineligible_n = source_ineligible_n,
    target_eligible_n = target_eligible_n,
    backbone_degenerate = backbone_degenerate,
    pair_status = if (both) "available" else "unavailable",
    pair_reason = if (both) NA_character_ else "malformed_loading"
  )
}


# ---- verification -----------------------------------------------------

bail <- function(...) stop(paste0(...), call. = FALSE)

fmt <- function(v) paste(format(v, trim = TRUE), collapse = " | ")

## Hand-derived discrete facts, one row per ordered pair.
check_facts <- function(name, pe, want) {
  if (!identical(paste(pe$K_from, pe$K_to), paste(want$K_from, want$K_to))) {
    bail(name, ": edge list disagrees with the hand table.")
  }
  for (nm in setdiff(names(want), c("K_from", "K_to"))) {
    if (!identical(pe[[nm]], want[[nm]])) {
      bail(name, ": hand-computed `", nm, "` disagrees with the builder.\n",
           "  built: ", fmt(pe[[nm]]), "\n  hand:  ", fmt(want[[nm]]))
    }
  }
  invisible(TRUE)
}

## Every field of every edge against the naive reference implementation.
check_reference <- function(name, sweep, loadings, pips, Q0, K0, eps,
                            bifactor, pe) {
  Ks <- as.integer(sweep$K)
  strip <- function(m) if (isTRUE(bifactor)) m[, -1L, drop = FALSE] else m
  reg <- function(K) {
    Q <- if (K > K0) cbind(Q0, matrix(-1L, nrow(Q0), K - K0)) else Q0
    matrix(Q == -1L, nrow(Q0), K)
  }
  status_ok <- function(field, K) {
    identical(sweep[[field]][match(K, Ks)], "ok")
  }
  row <- 0L
  for (i in seq_along(Ks)) {
    for (j in seq_along(Ks)) {
      if (i >= j) next
      row <- row + 1L
      k <- Ks[i]
      l <- Ks[j]
      if (!identical(c(pe$K_from[row], pe$K_to[row]), c(k, l))) {
        bail(name, ": row ", row, " is not the (", k, ", ", l, ") edge.")
      }
      ref <- ref_pair_facts(
        strip(loadings[[as.character(k)]]), strip(loadings[[as.character(l)]]),
        strip(pips[[as.character(k)]]), strip(pips[[as.character(l)]]),
        reg(k), reg(l), K0, eps,
        status_ok("loading_status", k), status_ok("loading_status", l),
        status_ok("pip_status", k), status_ok("pip_status", l)
      )
      for (nm in names(ref)) {
        got <- pe[[nm]][row]
        exp <- ref[[nm]]
        agree <- if (nm %in% DOUBLE_PAIR_COLS) {
          identical(is.na(got), is.na(exp)) &&
            (is.na(got) || abs(got - exp) <= 1e-12 * max(1, abs(exp)))
        } else {
          identical(got, exp)
        }
        if (!agree) {
          bail(name, ": edge (", k, ", ", l, ") field `", nm,
               "` disagrees with the reference matcher.\n",
               "  package:   ", fmt(got), "\n  reference: ", fmt(exp))
        }
      }
    }
  }
  invisible(TRUE)
}


# ---- member assembly --------------------------------------------------

as_stub <- function(sweep, tr, pe, loadings, pips, Q0, bifactor, eps) {
  structure(
    list(sweep = sweep, transitions = tr, persistence = pe,
         loadings = loadings, pips = pips, Q0 = Q0,
         settings = list(bifactor = bifactor, stability_eps = eps),
         provenance = list()),
    class = "pefa"
  )
}

finish <- function(description, features, Q0, K0, bifactor, eps, Ks,
                   loadings, pips, sweep, want) {
  names(loadings) <- as.character(Ks)
  names(pips) <- as.character(Ks)
  tr <- vbpm:::.pefa_transitions_table(sweep)
  pe <- vbpm:::.pefa_persistence_table(sweep, loadings, pips, Q0, K0, eps,
                                       bifactor)
  check_facts(description, pe, want)
  check_reference(description, sweep, loadings, pips, Q0, K0, eps, bifactor,
                  pe)
  stub <- as_stub(sweep, tr, pe, loadings, pips, Q0, bifactor, eps)
  pivots <- lapply(METRICS, function(mt) vbpm::persistence(stub, mt))
  names(pivots) <- METRICS
  list(
    description = description,
    features = features,
    K0 = as.integer(K0),
    bifactor = bifactor,
    stability_eps = eps,
    Q0 = Q0,
    sweep = sweep,
    loadings = loadings,
    pips = pips,
    expected = list(transitions = tr, persistence = pe, pivots = pivots)
  )
}

## Build the PIP list that goes with a loading list, using each candidate's
## own full-width design.
pips_for <- function(loadings, Ks, Q0, K0, bifactor) {
  lapply(seq_along(Ks), function(i) {
    pip_of(loadings[[i]], design_of(Q0, K0, Ks[i], bifactor))
  })
}


# ---- shared column shapes ---------------------------------------------

Q0_2 <- matrix(-1L, J, 2L, dimnames = list(ITEMS, c("F1", "F2")))
Q0_2[1:2, 1] <- 1L
Q0_2[1:2, 2] <- 0L
Q0_2[4:5, 2] <- 1L
Q0_2[4:5, 1] <- 0L

Q0_1 <- matrix(-1L, J, 1L, dimnames = list(ITEMS, "F1"))
Q0_1[1:2, 1] <- 1L

Q0_0 <- matrix(integer(0), J, 0L, dimnames = list(ITEMS, character(0)))

b1 <- c(.85, .80, .75, 0, 0, 0, 0, 0)
b2 <- c(0, 0, 0, .82, .78, .74, 0, 0)
zero_col <- rep(0, J)

e_a <- c(0, 0, 0, 0, 0, 0, .85, .60)   # SSL 1.0825
s3 <- c(0, 0, 0, 0, 0, 0, .90, .30)    # SSL 0.9000
s4 <- c(0, 0, 0, 0, 0, 0, .75, .45)    # SSL 0.7650
t3 <- c(0, 0, 0, 0, 0, .60, 0, 0)      # SSL 0.3600
t4 <- c(0, 0, 0, 0, 0, 0, .85, .35)    # SSL 0.8450
t5 <- c(.50, 0, 0, 0, 0, 0, 0, 0)      # SSL 0.2500
w1 <- c(0, 0, 0, 0, 0, 0, .15, .10)    # SSL 0.0325  (below the 0.10 screen)
w2 <- c(.12, 0, 0, 0, 0, 0, 0, .08)    # SSL 0.0208  (below the 0.10 screen)

EPS <- 0.10

members <- list()


# ---- 1. all_ok --------------------------------------------------------
#
# K0 = 2 over K = 2:5, every component "ok".  The k = K0 edges are (2,3),
# (2,4) and (2,5): backbone-only comparisons with n_pairs == K0 and no
# reduced-pair-set marker.  The (4,5) edge is the collision: sources F3 and F4
# both sit on items 7-8, target F4 is the only item-7/8 column at K = 5, so
# both choose it.

lam <- list(
  jitter_lam(lam_of(list(b1, b2), 2L), 1101),
  jitter_lam(lam_of(list(b1, b2, e_a), 2L), 1102),
  jitter_lam(lam_of(list(b1, b2, s3, s4), 2L), 1103),
  jitter_lam(lam_of(list(b1, b2, t3, t4, t5), 2L), 1104)
)
Ks <- 2:5
members$all_ok <- finish(
  description = "all components ok; k = K0 edges and a collision edge",
  features = c("all_ok", "k_equals_K0_edge", "collision_edge"),
  Q0 = Q0_2, K0 = 2L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_2, 2L, FALSE),
  sweep = build_sweep(Ks, lam, 2L, FALSE),
  want = data.frame(
    K_from = c(2L, 2L, 2L, 3L, 3L, 4L),
    K_to = c(3L, 4L, 5L, 4L, 5L, 5L),
    r = c(1L, 2L, 3L, 1L, 2L, 1L),
    n_pairs = c(2L, 2L, 2L, 3L, 3L, 4L),
    unmatched_n = c(1L, 2L, 3L, 1L, 2L, 2L),
    collision = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    n_collisions = c(0L, 0L, 0L, 0L, 0L, 1L),
    collision_targets = c("", "", "", "", "", "4"),
    collision_multiplicities = c("", "", "", "", "", "2"),
    source_ineligible_n = c(0L, 0L, 0L, 0L, 0L, 0L),
    target_eligible_n = c(1L, 2L, 3L, 2L, 3L, 3L),
    backbone_degenerate = rep(FALSE, 6),
    from_fit_ok = rep(TRUE, 6), to_fit_ok = rep(TRUE, 6),
    pair_status = rep("available", 6),
    pair_reason = rep(NA_character_, 6),
    stringsAsFactors = FALSE
  )
)


# ---- 2. nonconverged_endpoint -----------------------------------------
#
# K = 3 stopped at max_it but returned a finite well-formed loading matrix.
# fit_status is unavailable/nonconverged while every other component stays ok,
# and both edges touching it stay fully computed, flagged only through
# from_fit_ok / to_fit_ok.

lam <- list(
  jitter_lam(lam_of(list(b1, b2), 2L), 1201),
  jitter_lam(lam_of(list(b1, b2, e_a), 2L), 1202),
  jitter_lam(lam_of(list(b1, b2, s3, s4), 2L), 1203)
)
Ks <- 2:4
members$nonconverged_endpoint <- finish(
  description = "a nonconverged candidate whose pair metrics stay computed",
  features = c("nonconverged_endpoint", "fully_computed_pairs"),
  Q0 = Q0_2, K0 = 2L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_2, 2L, FALSE),
  sweep = build_sweep(Ks, lam, 2L, FALSE,
                      converged = c(TRUE, FALSE, TRUE)),
  want = data.frame(
    K_from = c(2L, 2L, 3L), K_to = c(3L, 4L, 4L), r = c(1L, 2L, 1L),
    n_pairs = c(2L, 2L, 3L),
    unmatched_n = c(1L, 2L, 1L),
    collision = c(FALSE, FALSE, FALSE),
    n_collisions = c(0L, 0L, 0L),
    collision_targets = c("", "", ""),
    collision_multiplicities = c("", "", ""),
    source_ineligible_n = c(0L, 0L, 0L),
    target_eligible_n = c(1L, 2L, 2L),
    backbone_degenerate = rep(FALSE, 3),
    from_fit_ok = c(TRUE, TRUE, FALSE),
    to_fit_ok = c(FALSE, TRUE, TRUE),
    pair_status = rep("available", 3),
    pair_reason = rep(NA_character_, 3),
    stringsAsFactors = FALSE
  )
)


# ---- 3. malformed_loading ---------------------------------------------
#
# K = 4 returned a wrongly shaped Lam, so its stored payload is the
# expected-shape all-NA matrix and loading_status is
# unavailable/malformed_loading -- while its fit converged and its PIP is
# valid.  Its ELBO is separately unavailable, which makes both ELBO gains
# nonfinite and exercises the empty-max() guard in .pefa_gain_pct().

lam_ok3 <- jitter_lam(lam_of(list(b1, b2, e_a), 2L), 1301)
lam_int4 <- jitter_lam(lam_of(list(b1, b2, s3, s4), 2L), 1302)
lam_ok5 <- jitter_lam(lam_of(list(b1, b2, t3, t4, t5), 2L), 1303)
lam_bad4 <- matrix(NA_real_, J, 4L, dimnames = dimnames(lam_int4))
lam <- list(lam_ok3, lam_bad4, lam_ok5)
Ks <- 3:5
pips_ml <- list(
  pip_of(lam_ok3, design_of(Q0_2, 2L, 3L, FALSE)),
  pip_of(lam_int4, design_of(Q0_2, 2L, 4L, FALSE)),
  pip_of(lam_ok5, design_of(Q0_2, 2L, 5L, FALSE))
)
members$malformed_loading <- finish(
  description = "a malformed loading endpoint with an independently valid PIP",
  features = c("malformed_loading_endpoint", "unavailable_pair",
               "all_nonfinite_elbo_gain_path"),
  Q0 = Q0_2, K0 = 2L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_ml,
  sweep = build_sweep(
    Ks, lam, 2L, FALSE,
    loading_reason = c(NA, "malformed_loading", NA),
    elbo_reason = c(NA, "nonfinite_elbo", NA)
  ),
  want = data.frame(
    K_from = c(3L, 3L, 4L), K_to = c(4L, 5L, 5L), r = c(1L, 2L, 1L),
    n_pairs = c(NA_integer_, 3L, NA_integer_),
    unmatched_n = c(NA_integer_, 2L, NA_integer_),
    collision = c(NA, FALSE, NA),
    n_collisions = c(NA_integer_, 0L, NA_integer_),
    collision_targets = c(NA_character_, "", NA_character_),
    collision_multiplicities = c(NA_character_, "", NA_character_),
    source_ineligible_n = c(0L, 0L, NA_integer_),
    target_eligible_n = c(NA_integer_, 3L, 3L),
    backbone_degenerate = c(NA, FALSE, NA),
    from_fit_ok = rep(TRUE, 3), to_fit_ok = rep(TRUE, 3),
    pair_status = c("unavailable", "available", "unavailable"),
    pair_reason = c("malformed_loading", NA_character_, "malformed_loading"),
    stringsAsFactors = FALSE
  )
)


# ---- 4. degenerate_backbone -------------------------------------------
#
# The second backbone column of K = 3 has exact zero norm.  Both edges that
# touch it keep pair_status "available" and lose only the congruence, while
# the (2,4) edge in the same sweep keeps a finite phi_min for contrast.

lam3_deg <- jitter_lam(lam_of(list(b1, b2, e_a), 2L), 1402)
lam3_deg[, 2L] <- 0
lam <- list(
  jitter_lam(lam_of(list(b1, b2), 2L), 1401),
  lam3_deg,
  jitter_lam(lam_of(list(b1, b2, s3, s4), 2L), 1403)
)
Ks <- 2:4
members$degenerate_backbone <- finish(
  description = "a finite zero-norm backbone column at one endpoint",
  features = c("degenerate_backbone", "metric_specific_availability"),
  Q0 = Q0_2, K0 = 2L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_2, 2L, FALSE),
  sweep = build_sweep(Ks, lam, 2L, FALSE),
  want = data.frame(
    K_from = c(2L, 2L, 3L), K_to = c(3L, 4L, 4L), r = c(1L, 2L, 1L),
    n_pairs = c(2L, 2L, 3L),
    unmatched_n = c(1L, 2L, 1L),
    collision = rep(FALSE, 3),
    n_collisions = rep(0L, 3),
    collision_targets = rep("", 3),
    collision_multiplicities = rep("", 3),
    source_ineligible_n = c(0L, 0L, 0L),
    target_eligible_n = c(1L, 2L, 2L),
    backbone_degenerate = c(TRUE, FALSE, TRUE),
    from_fit_ok = rep(TRUE, 3), to_fit_ok = rep(TRUE, 3),
    pair_status = rep("available", 3),
    pair_reason = rep(NA_character_, 3),
    stringsAsFactors = FALSE
  )
)


# ---- 5. bifactor_thin_backbone ----------------------------------------
#
# K0 = 1 group backbone under bifactor = TRUE, so stored matrices carry K + 1
# columns with G first.  The backbone column is small but nonzero (SSL ~ .06,
# below the .10 screen it is never subject to) and the general column is
# smaller still (SSL ~ .02), so backbone_ssl_min taken over stored columns
# 1:K0 instead of 2:(K0 + 1) would return a different number.  The (3,4) edge
# reuses GROUP target 2, which is STORED column 3.

g_col <- c(.08, .08, .08, 0, 0, 0, 0, 0)   # SSL 0.0192, the general factor
g1 <- c(.20, .10, .10, 0, 0, 0, 0, 0)      # SSL 0.0600, the thin backbone
g2 <- c(0, 0, 0, .90, .30, 0, 0, 0)        # SSL 0.9000
g3 <- c(0, 0, 0, .80, .40, 0, 0, 0)        # SSL 0.8000
u2 <- c(0, 0, 0, .85, .35, 0, 0, 0)        # SSL 0.8450, reused target
u3 <- c(0, 0, 0, 0, 0, .75, 0, 0)          # SSL 0.5625
u4 <- c(0, 0, 0, 0, 0, 0, .60, .25)        # SSL 0.4225

lam <- list(
  jitter_lam(lam_of(list(g_col, g1), 1L, TRUE), 1501, sd = 0.004),
  jitter_lam(lam_of(list(g_col, g1, g2), 1L, TRUE), 1502, sd = 0.004),
  jitter_lam(lam_of(list(g_col, g1, g2, g3), 1L, TRUE), 1503, sd = 0.004),
  jitter_lam(lam_of(list(g_col, g1, u2, u3, u4), 1L, TRUE), 1504, sd = 0.004)
)
Ks <- 1:4
members$bifactor_thin_backbone <- finish(
  description = paste("bifactor sweep with a small-but-nonzero backbone and a",
                      "group-block collision"),
  features = c("bifactor", "thin_backbone", "backbone_ssl_min_offset",
               "k_equals_K0_edge", "collision_edge"),
  Q0 = Q0_1, K0 = 1L, bifactor = TRUE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_1, 1L, TRUE),
  sweep = build_sweep(Ks, lam, 1L, TRUE),
  want = data.frame(
    K_from = c(1L, 1L, 1L, 2L, 2L, 3L),
    K_to = c(2L, 3L, 4L, 3L, 4L, 4L),
    r = c(1L, 2L, 3L, 1L, 2L, 1L),
    n_pairs = c(1L, 1L, 1L, 2L, 2L, 3L),
    unmatched_n = c(1L, 2L, 3L, 1L, 2L, 2L),
    collision = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    n_collisions = c(0L, 0L, 0L, 0L, 0L, 1L),
    collision_targets = c("", "", "", "", "", "2"),
    collision_multiplicities = c("", "", "", "", "", "2"),
    source_ineligible_n = rep(0L, 6),
    target_eligible_n = c(1L, 2L, 3L, 2L, 3L, 3L),
    backbone_degenerate = rep(FALSE, 6),
    from_fit_ok = rep(TRUE, 6), to_fit_ok = rep(TRUE, 6),
    pair_status = rep("available", 6),
    pair_reason = rep(NA_character_, 6),
    stringsAsFactors = FALSE
  )
)


# ---- 6. ineligible_source ---------------------------------------------
#
# The single exploratory column of K = 3 sits below the screen.  As a source
# it is dropped from P_kl (source_ineligible_n = 1, a reduced pair set); as a
# target it still counts in the unmatched set, which section 2.7 says includes
# screened-out columns.

lam <- list(
  jitter_lam(lam_of(list(b1, b2), 2L), 1601),
  jitter_lam(lam_of(list(b1, b2, w1), 2L), 1602, sd = 0.004),
  jitter_lam(lam_of(list(b1, b2, s3, s4), 2L), 1603)
)
Ks <- 2:4
members$ineligible_source <- finish(
  description = "an exploratory source below the SSL screen",
  features = c("ineligible_source", "reduced_pair_set",
               "screened_target_counted_unmatched"),
  Q0 = Q0_2, K0 = 2L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_2, 2L, FALSE),
  sweep = build_sweep(Ks, lam, 2L, FALSE),
  want = data.frame(
    K_from = c(2L, 2L, 3L), K_to = c(3L, 4L, 4L), r = c(1L, 2L, 1L),
    n_pairs = c(2L, 2L, 2L),
    unmatched_n = c(1L, 2L, 2L),
    collision = rep(FALSE, 3),
    n_collisions = rep(0L, 3),
    collision_targets = rep("", 3),
    collision_multiplicities = rep("", 3),
    source_ineligible_n = c(0L, 0L, 1L),
    target_eligible_n = c(0L, 2L, 2L),
    backbone_degenerate = rep(FALSE, 3),
    from_fit_ok = rep(TRUE, 3), to_fit_ok = rep(TRUE, 3),
    pair_status = rep("available", 3),
    pair_reason = rep(NA_character_, 3),
    stringsAsFactors = FALSE
  )
)


# ---- 7. scarce_targets ------------------------------------------------
#
# K = 4 has two exploratory columns and both are below the screen, so an
# eligible source at K = 3 has nowhere to go: the (3,4) edge keeps the
# backbone-only measurement plus target_eligible_n = 0.

lam <- list(
  jitter_lam(lam_of(list(b1, b2, e_a), 2L), 1701),
  jitter_lam(lam_of(list(b1, b2, w1, w2), 2L), 1702, sd = 0.004),
  jitter_lam(lam_of(list(b1, b2, t3, t4, t5), 2L), 1703)
)
Ks <- 3:5
members$scarce_targets <- finish(
  description = "an eligible source with zero eligible exploratory targets",
  features = c("scarce_targets", "reduced_pair_set", "ineligible_source"),
  Q0 = Q0_2, K0 = 2L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_2, 2L, FALSE),
  sweep = build_sweep(Ks, lam, 2L, FALSE),
  want = data.frame(
    K_from = c(3L, 3L, 4L), K_to = c(4L, 5L, 5L), r = c(1L, 2L, 1L),
    n_pairs = c(2L, 3L, 2L),
    unmatched_n = c(2L, 2L, 3L),
    collision = rep(FALSE, 3),
    n_collisions = rep(0L, 3),
    collision_targets = rep("", 3),
    collision_multiplicities = rep("", 3),
    source_ineligible_n = c(0L, 0L, 2L),
    target_eligible_n = c(0L, 3L, 3L),
    backbone_degenerate = rep(FALSE, 3),
    from_fit_ok = rep(TRUE, 3), to_fit_ok = rep(TRUE, 3),
    pair_status = rep("available", 3),
    pair_reason = rep(NA_character_, 3),
    stringsAsFactors = FALSE
  )
)


# ---- 8. k0_empty_pairs ------------------------------------------------
#
# Fully exploratory with the only source column below the screen: P_kl is
# empty, so the four matched/PIP metrics are typed NA while collision is a
# known clean zero and the unmatched evidence is still computed.

v1 <- c(.80, .70, 0, 0, 0, 0, 0, 0)   # SSL 1.1300
v2 <- c(0, 0, 0, .75, .65, 0, 0, 0)   # SSL 0.9850
lam <- list(
  jitter_lam(lam_of(list(w1), 0L), 1801, sd = 0.004),
  jitter_lam(lam_of(list(v1, v2), 0L), 1802)
)
Ks <- 1:2
members$k0_empty_pairs <- finish(
  description = "K0 = 0 with an empty scored pair set",
  features = c("k0_zero", "empty_pair_set", "reduced_pair_set"),
  Q0 = Q0_0, K0 = 0L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_0, 0L, FALSE),
  sweep = build_sweep(Ks, lam, 0L, FALSE),
  want = data.frame(
    K_from = 1L, K_to = 2L, r = 1L,
    n_pairs = 0L,
    unmatched_n = 2L,
    collision = FALSE,
    n_collisions = 0L,
    collision_targets = "",
    collision_multiplicities = "",
    source_ineligible_n = 1L,
    target_eligible_n = 2L,
    backbone_degenerate = FALSE,
    from_fit_ok = TRUE, to_fit_ok = TRUE,
    pair_status = "available",
    pair_reason = NA_character_,
    stringsAsFactors = FALSE
  )
)


# ---- 9. k0_nonempty_pairs ---------------------------------------------
#
# The same fully exploratory shape with every column eligible, so the matched
# metrics are ordinary finite numbers.

v1a <- c(.80, .70, .20, 0, 0, 0, 0, 0)
v1b <- c(.78, .72, .15, 0, 0, 0, 0, 0)
v1c <- c(.76, .70, .18, 0, 0, 0, 0, 0)
v2c <- c(0, 0, 0, .72, .68, 0, 0, 0)
v3 <- c(0, 0, 0, 0, 0, .70, .30, 0)
lam <- list(
  jitter_lam(lam_of(list(v1a), 0L), 1901),
  jitter_lam(lam_of(list(v1b, v2), 0L), 1902),
  jitter_lam(lam_of(list(v1c, v2c, v3), 0L), 1903)
)
Ks <- 1:3
members$k0_nonempty_pairs <- finish(
  description = "K0 = 0 with a nonempty scored pair set",
  features = c("k0_zero", "nonempty_pair_set"),
  Q0 = Q0_0, K0 = 0L, bifactor = FALSE, eps = EPS, Ks = Ks,
  loadings = lam, pips = pips_for(lam, Ks, Q0_0, 0L, FALSE),
  sweep = build_sweep(Ks, lam, 0L, FALSE),
  want = data.frame(
    K_from = c(1L, 1L, 2L), K_to = c(2L, 3L, 3L), r = c(1L, 2L, 1L),
    n_pairs = c(1L, 1L, 2L),
    unmatched_n = c(1L, 2L, 1L),
    collision = rep(FALSE, 3),
    n_collisions = rep(0L, 3),
    collision_targets = rep("", 3),
    collision_multiplicities = rep("", 3),
    source_ineligible_n = rep(0L, 3),
    target_eligible_n = c(2L, 3L, 3L),
    backbone_degenerate = rep(FALSE, 3),
    from_fit_ok = rep(TRUE, 3), to_fit_ok = rep(TRUE, 3),
    pair_status = rep("available", 3),
    pair_reason = rep(NA_character_, 3),
    stringsAsFactors = FALSE
  )
)


# ---- 10. the discriminating integrity object --------------------------
#
# Release-plan section 8.1 and revision-contract section 6.1: one complete,
# classed, current-schema `pefa` object whose two-candidate geometry decides
# between the 0.9.0 matcher rule and the rule 0.8.3 retracted, plus the exact
# deterministic `Y` behind its standardized-data fingerprint.
#
# The geometry is the plan's, unchanged.  On backbone (.7, .7, .7, 0, 0, 0)
# with exploratory source a = (.6, .6, 0, 0, 0, 0) and the two targets
# t1 = 2a and t2 = (.55, .55, .30, 0, 0, 0), the raw squared distances are
# .72 and .095 while the absolute congruences are 1 and .9330078.  The two
# rules therefore pick different targets by construction: 0.9.0 assigns by
# max |phi| and reports phi_min = 1; the retracted minimum-raw-distance rule
# assigns t2 and reports phi_min = 0.93300782.
#
# ---------------------------------------------------------------------------
# WHY THIS OBJECT IS CONSTRUCTED AND RE-SEALED RATHER THAN FITTED
# ---------------------------------------------------------------------------
# The soundest route would be to obtain the object from pefa() itself, so that
# its provenance is a byproduct of a real sweep.  That route cannot produce
# this fixture: the discriminator is a statement about four EXACT loading
# columns, and a variational fit returns whatever the data and the
# spike-and-slab path give it.  `conformance` below runs a real default
# pefa() on this very Q0 and Y and records how far its K = 3 estimate lands
# from the target geometry -- the fitted columns are nowhere near collinear
# with 2a, so no seed, v0 path, or iteration budget turns a fit into this
# case.  Screening seeds until a fit happened to look close would also be the
# worst of both worlds: still not the plan's numbers, and no longer a
# hand-checkable case.
#
# So the seven scientific components are constructed, exactly as the low-level
# members above are, and the object is then SEALED WITH THE PACKAGE'S OWN
# .pefa_provenance() over the same payload pefa() seals.  Nothing is
# hand-written into a hash field:
#
#   * $provenance is the return value of .pefa_provenance(), so the payload
#     assembly, the `secs` exclusion, the canonicalization, and both the
#     evidence and lineage hashes are the sealing code's, not a copy;
#   * standardized_data_sha256 is .pefa_standardized_hash(Y_int) on a
#     real 60 x 6 matrix that is stored beside the object, so
#     verify_pefa(object, Y) checks a fingerprint that genuinely came from
#     that matrix rather than one written to match;
#   * Q0_sha256 and settings_sha256 are .pefa_sha256() of the stored $Q0 and
#     $settings;
#   * $transitions and $persistence are built by the package's own builders
#     from the stored $sweep, $loadings, $pips, $Q0 and $settings, which is
#     exactly what verify_pefa() rebuilds and compares against.
#
# Every hash is therefore genuinely correct FOR THE STORED CONTENT, which is
# the whole of what verify_pefa() asserts.  What the object does not claim,
# and what no field of it states, is that a vbfa() run produced these
# columns; the honest description of it is the one carried in
# $integrity$route below and repeated in the tests.  Y is what section 6.1
# says it is -- the exact deterministic input whose internal standardization
# produced the stored data hash -- and it is drawn from this same backbone
# and source, so it is the right data for the object to be bound to.
#
# Determinism: every scientific number here is a literal or comes from a
# set.seed()ed rnorm, and the conformance fit contributes no stored byte -- it
# compares only names, storage modes, and the settings list, none of which any
# BLAS can move.  Two runs of this script against one installation therefore
# write the identical file, which was checked.  The two fields that are
# deliberately installation-dependent, $provenance$package_version and
# $package_build_id, move with the installed DESCRIPTION as they are meant to;
# the first is inside evidence_sha256 and the second outside both hashes, so
# regenerating after a version bump reseals consistently and the shipped
# object keeps verifying either way.

J6 <- 6L
ITEMS6 <- paste0("item_", seq_len(J6))

int_backbone <- c(.70, .70, .70, 0, 0, 0)   # SSL 1.4700
int_source   <- c(.60, .60, 0, 0, 0, 0)     # SSL 0.7200
int_target_1 <- 2 * int_source              # SSL 2.8800, collinear with a
int_target_2 <- c(.55, .55, .30, 0, 0, 0)   # SSL 0.6950, the near miss

## The plan's four numbers, restated here and checked rather than trusted.
INT_RAW_SQ_DISTANCE <- c(target_1 = 0.72, target_2 = 0.095)
INT_ABS_CONGRUENCE <- c(target_1 = 1, target_2 = 0.9330078)
INT_RETRACTED_PHI_MIN <- 0.93300782

int_raw <- function(x, y) sum((x - y)^2)
int_cong <- function(x, y) abs(sum(x * y)) / sqrt(sum(x^2) * sum(y^2))

if (!isTRUE(all.equal(unname(INT_RAW_SQ_DISTANCE),
                      c(int_raw(int_source, int_target_1),
                        int_raw(int_source, int_target_2)),
                      tolerance = 1e-12))) {
  bail("integrity: the raw squared distances are not the plan's .72/.095.")
}
if (!isTRUE(all.equal(unname(INT_ABS_CONGRUENCE),
                      c(int_cong(int_source, int_target_1),
                        int_cong(int_source, int_target_2)),
                      tolerance = 1e-7))) {
  bail("integrity: the congruences are not the plan's 1/.9330078.")
}
## The retracted rule's answer is the near miss, to the digits the plan
## quotes; the two rules must be numerically far apart, not merely different.
if (abs(int_cong(int_source, int_target_2) - INT_RETRACTED_PHI_MIN) > 1e-8) {
  bail("integrity: 0.93300782 is not congruence(a, t2).")
}

## K0 = 1 with two anchors, one regularized backbone cell, and three
## administrative zeros, so the matched backbone pair contributes real PIP
## cells and the fixed cells are excluded from pip_rmsd.
Q0_int <- matrix(c(1L, 1L, -1L, 0L, 0L, 0L), J6, 1L,
                 dimnames = list(ITEMS6, "F1"))
K0_INT <- 1L

int_lam <- function(cols, K) {
  m <- matrix(unlist(cols), J6, length(cols))
  storage.mode(m) <- "double"
  dimnames(m) <- list(ITEMS6, paste0("F", seq_len(K)))
  m
}
int_loadings <- list(
  `2` = int_lam(list(int_backbone, int_source), 2L),
  `3` = int_lam(list(int_backbone, int_target_1, int_target_2), 3L)
)
int_pips <- list(
  `2` = pip_of(int_loadings[["2"]], design_of(Q0_int, K0_INT, 2L, FALSE)),
  `3` = pip_of(int_loadings[["3"]], design_of(Q0_int, K0_INT, 3L, FALSE))
)

## A real 60 x 6 matrix drawn from the same backbone and source the object's
## K = 2 candidate describes.  It is the object's data of record: its
## standardization is what produced $provenance$standardized_data_sha256.
set.seed(2001)
N_INT <- 60L
int_factors <- matrix(stats::rnorm(N_INT * 2L), N_INT, 2L)
int_resid_sd <- sqrt(pmax(0.10, 1 - (int_backbone^2 + int_source^2)))
Y_int <- int_factors %*% rbind(int_backbone, int_source) +
  matrix(stats::rnorm(N_INT * J6), N_INT, J6) * rep(int_resid_sd, each = N_INT)
colnames(Y_int) <- ITEMS6

## The ten-field settings list pefa() stores for a default non-bifactor call.
## `conformance` below proves this is identical to a real one rather than
## merely plausible.
settings_int <- list(
  bifactor = FALSE, general = NULL, v0 = c(0.01, 0.005, 0.002, 0.001),
  max_it = 10000L, convChk = FALSE, tolVal = 1e-4, tau = 0.50,
  stability_eps = EPS, rank_adjust = FALSE, rank_max_J = 100L
)

sweep_int <- build_sweep(2:3, int_loadings, K0_INT, FALSE)

int_persistence <- vbpm:::.pefa_persistence_table(
  sweep_int, int_loadings, int_pips, Q0_int, K0_INT, EPS, FALSE
)

## The same two checks every member gets: the hand table, then the naive
## reference matcher.
check_facts(
  "integrity", int_persistence,
  data.frame(
    K_from = 2L, K_to = 3L, r = 1L,
    n_pairs = 2L,              # the backbone pair plus a -> t1
    unmatched_n = 1L,          # t2 is the target no source took
    collision = FALSE,
    n_collisions = 0L,
    collision_targets = "",
    collision_multiplicities = "",
    source_ineligible_n = 0L,
    target_eligible_n = 2L,
    backbone_degenerate = FALSE,
    from_fit_ok = TRUE, to_fit_ok = TRUE,
    pair_status = "available",
    pair_reason = NA_character_,
    stringsAsFactors = FALSE
  )
)
check_reference("integrity", sweep_int, int_loadings, int_pips, Q0_int,
                K0_INT, EPS, FALSE, int_persistence)

## The discrimination itself, on the built table.
##
## phi_min is 1 to within two ulps rather than bit-exactly 1: the congruence
## of a with 2a is formed as |a'b| / (||a|| ||b||) on rescaled copies, and
## sqrt(2) * sqrt(2) is 2 + 4.4e-16, so the ratio lands at 1 - 2.2e-16 and the
## clamp at 1 never fires.  The plan's tolerance for this cell is 1e-10, and
## the retracted answer is 0.067 away, so nothing here rests on the last bit.
INT_PHI_MIN <- int_persistence$phi_min[1L]
if (!isTRUE(all.equal(INT_PHI_MIN, 1, tolerance = 1e-12))) {
  bail("integrity: phi_min is ", format(INT_PHI_MIN, digits = 17),
       ", not the max-congruence answer 1.")
}
if (abs(INT_PHI_MIN - INT_RETRACTED_PHI_MIN) < 0.06) {
  bail("integrity: phi_min is the retracted minimum-raw-distance answer.")
}
## The assignment, the unmatched target, and its size all follow the same
## choice, so the row discriminates in three independent places.
if (!identical(int_persistence$unmatched_n[1L], 1L) ||
    !isTRUE(all.equal(int_persistence$unmatched_ssl[1L], sum(int_target_2^2),
                      tolerance = 1e-12))) {
  bail("integrity: the unmatched target is not t2.")
}

## The low-level case, called exactly as the builder calls it, so the raw
## matcher is covered beside the sealed object.
int_reg <- function(K) {
  Q <- if (K > K0_INT) {
    cbind(Q0_int, matrix(-1L, J6, K - K0_INT))
  } else {
    Q0_int
  }
  matrix(Q == -1L, J6, K)
}
int_facts <- vbpm:::.pefa_pair_facts(
  lam_from = int_loadings[["2"]], lam_to = int_loadings[["3"]],
  pip_from = int_pips[["2"]], pip_to = int_pips[["3"]],
  reg_from = int_reg(2L), reg_to = int_reg(3L),
  K0 = K0_INT, stability_eps = EPS,
  loading_ok_from = TRUE, loading_ok_to = TRUE,
  pip_ok_from = TRUE, pip_ok_to = TRUE
)
if (!identical(int_facts$.assign$target, c(1L, 2L))) {
  bail("integrity: the source did not take target 1 (stored column 2).")
}
if (!identical(int_facts$phi_min, INT_PHI_MIN)) {
  bail("integrity: the direct matcher call disagrees with the built table.")
}

## Seal.  Everything below is the package's own hashing over the payload
## above; nothing is written by hand.
int_evidence <- list(
  sweep = sweep_int,
  transitions = vbpm:::.pefa_transitions_table(sweep_int),
  persistence = int_persistence,
  loadings = int_loadings,
  pips = int_pips,
  Q0 = Q0_int,
  settings = settings_int
)
int_object <- c(int_evidence, list(provenance = vbpm:::.pefa_provenance(
  evidence = int_evidence,
  standardized_data_sha256 = vbpm:::.pefa_standardized_hash(Y_int),
  Q0_sha256 = vbpm:::.pefa_sha256(Q0_int),
  settings_sha256 = vbpm:::.pefa_sha256(settings_int),
  package_version = as.character(utils::packageVersion("vbpm")),
  newly_fitted_K = 2:3
)))
class(int_object) <- "pefa"

## The acceptance gate, run here so a bad object never reaches the .rds.
if (!isTRUE(vbpm::verify_pefa(int_object))) {
  bail("integrity: verify_pefa(object) did not return TRUE.")
}
if (!isTRUE(vbpm::verify_pefa(int_object, Y_int))) {
  bail("integrity: verify_pefa(object, Y) did not return TRUE.")
}
if (isTRUE(tryCatch(vbpm::verify_pefa(int_object, Y_int + 1),
                    error = function(e) FALSE))) {
  bail("integrity: a wrong Y passed the data-hash check.")
}

## The negative copy: the persistence row rewritten with the retracted rule's
## answer.  It must fail, and fail as a persistence/evidence mismatch rather
## than as a shape or schema complaint.
int_retracted <- int_object
int_retracted$persistence$phi_min[1L] <- INT_RETRACTED_PHI_MIN
int_retracted_message <- tryCatch(vbpm::verify_pefa(int_retracted),
                                  error = conditionMessage)
if (!is.character(int_retracted_message) ||
    !grepl("`$persistence` disagrees with the table rebuilt",
           int_retracted_message, fixed = TRUE) ||
    !grepl("`$provenance$evidence_sha256` does not match",
           int_retracted_message, fixed = TRUE)) {
  bail("integrity: the retracted copy did not fail with a persistence/",
       "evidence mismatch.\n  got: ",
       paste(int_retracted_message, collapse = " "))
}

## Conformance: a real default pefa() call on this Q0 and Y, used only to
## check the hand-built schema against a live one and to record why the fit
## route was rejected.  Its numbers are never stored.
conformance <- suppressWarnings(
  pefa(Q0_int, Y_int, 2, 3, verbose = FALSE)
)
if (!identical(conformance$settings, settings_int)) {
  bail("integrity: $settings is not what a real default pefa() call stores.")
}
if (!identical(names(conformance$sweep), names(sweep_int)) ||
    !identical(vapply(conformance$sweep, typeof, ""),
               vapply(sweep_int, typeof, ""))) {
  bail("integrity: $sweep is not the current-schema candidate table.")
}
if (!identical(names(conformance), names(int_object)) ||
    !identical(class(conformance), class(int_object)) ||
    !identical(names(conformance$provenance), names(int_object$provenance)) ||
    !identical(dimnames(conformance$loadings[["3"]]),
               dimnames(int_loadings[["3"]])) ||
    !identical(conformance$Q0, Q0_int)) {
  bail("integrity: the object shape is not the one pefa() issues.")
}
## And the reason for the route: the fit does not land on the geometry.
int_fit_gap <- max(abs(conformance$loadings[["3"]] - int_loadings[["3"]]))
if (!(int_fit_gap > 0.1)) {
  bail("integrity: the fitted candidate now reproduces the constructed ",
       "geometry; rebuild this case through pefa() instead.")
}

integrity <- list(
  description = paste(
    "The section 8.1 discriminator as a complete sealed pefa object:",
    "backbone (.7,.7,.7,0,0,0), exploratory source (.6,.6,0,0,0,0), targets",
    "2*source and (.55,.55,.30,0,0,0), K0 = 1, stability_eps = .10.",
    "The 0.9.0 max-absolute-congruence rule reports phi_min = 1; the",
    "retracted minimum-raw-distance rule reports 0.93300782."
  ),
  route = paste(
    "Constructed and re-sealed, not fitted. The seven scientific components",
    "are built by hand (a variational fit cannot return these exact loading",
    "columns); $transitions and $persistence come from the package's own",
    "builders and $provenance from .pefa_provenance(), so every stored hash",
    "is genuinely correct for the stored content and $Y is the real matrix",
    "whose .pefa_standardized_hash() produced standardized_data_sha256.",
    "See data-raw/make_r8_fixture.R for the full argument."
  ),
  features = c("discriminating_max_congruence", "sealed_pefa_object",
               "low_level_matcher_case", "retracted_min_distance_copy"),
  K0 = K0_INT,
  stability_eps = EPS,
  columns = list(backbone = int_backbone, source = int_source,
                 target_1 = int_target_1, target_2 = int_target_2),
  rules = list(
    raw_sq_distance = INT_RAW_SQ_DISTANCE,
    abs_congruence = INT_ABS_CONGRUENCE,
    max_congruence_target = 2L,      # stored column of t1
    min_distance_target = 3L,        # stored column of t2
    phi_min = INT_PHI_MIN,
    retracted_phi_min = INT_RETRACTED_PHI_MIN
  ),
  object = int_object,
  Y = Y_int,
  expected = list(
    persistence = int_persistence,
    transitions = int_evidence$transitions,
    phi_min = INT_PHI_MIN,
    n_pairs = 2L,
    unmatched_n = 1L,
    unmatched_ssl = int_persistence$unmatched_ssl[1L],
    assign_target = int_facts$.assign$target
  ),
  matcher_case = list(
    lam_from = int_loadings[["2"]], lam_to = int_loadings[["3"]],
    pip_from = int_pips[["2"]], pip_to = int_pips[["3"]],
    reg_from = int_reg(2L), reg_to = int_reg(3L),
    K0 = K0_INT, stability_eps = EPS,
    expected = int_facts[setdiff(names(int_facts), ".assign")],
    expected_assign = int_facts$.assign
  ),
  retracted_copy = int_retracted,
  retracted_copy_note = paste(
    "DELIBERATELY CORRUPT. Its persistence phi_min carries the retracted",
    "minimum-raw-distance answer, so verify_pefa() must reject it with a",
    "persistence/evidence mismatch. Never use it as evidence."
  )
)


# ---- write ------------------------------------------------------------

fixture <- list(
  fixture_version = FIXTURE_VERSION,
  evidence_schema_id = EVIDENCE_SCHEMA_ID,
  metrics = METRICS,
  double_pair_cols = DOUBLE_PAIR_COLS,
  members = members,
  integrity = integrity
)

root <- find_root()
dir.create(file.path(root, "inst", "extdata"), recursive = TRUE,
           showWarnings = FALSE)
out_file <- file.path(root, "inst", "extdata", "vbpm_r8_fixture.rds")
saveRDS(fixture, out_file, version = 3L, compress = "xz")

cat("wrote ", out_file, "\n", sep = "")
cat("members: ", paste(names(members), collapse = ", "), "\n", sep = "")
cat("persistence rows: ",
    sum(vapply(members, function(m) nrow(m$expected$persistence), integer(1))),
    "\n", sep = "")
cat("integrity: phi_min = ", format(INT_PHI_MIN, digits = 17),
    " (retracted rule: ", format(INT_RETRACTED_PHI_MIN, digits = 9), ")\n",
    sep = "")
cat("integrity: verify_pefa() passes with and without Y; the retracted copy ",
    "fails.\n", sep = "")
cat("integrity: fitted-vs-constructed loading gap = ",
    format(int_fit_gap, digits = 4), "\n", sep = "")
