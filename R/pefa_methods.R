## Description and display for the object returned by pefa().  Everything in
## this file reads stored fields: nothing here refits, rematches columns, or
## applies an analysis rule, so a display can never disagree with $transitions
## or $persistence.

#' Sums of squared loadings for a PEFA sweep
#'
#' Descriptive column sizes for every candidate in a [pefa()] sweep:
#' `colSums(Lam^2)` on the stored loading matrix of each candidate. They are
#' not Gram eigenvalues, they take part in no correspondence or selection
#' rule, and no cutoff is applied to them. `summary()` calls this function, so
#' there is one definition.
#'
#' A low `phi_min` in `$transitions` or `$persistence$phi` is a weakest-link
#' reading over every source column, including one the fit has driven to
#' nearly zero. Such a column matches the smallest-norm target and carries an
#' essentially arbitrary congruence. `ssl()` is what separates "this candidate
#' contains a column with nothing in it" from "the retained structure moved".
#'
#' @param x A result returned by [pefa()].
#'
#' @return A list with one element per candidate, named as `x$loadings` is
#'   (the candidate `K` as a character string). Each element is the double
#'   vector `colSums(L^2)` of that candidate's loading matrix and keeps the
#'   stored loading-column names when the matrix has them. In bifactor mode
#'   each vector covers all `K + 1` stored columns with the general factor
#'   first, so the group-block scale used by every between-candidate
#'   comparison is `ssl(x)[[k]][-1]`.
#'
#' @seealso [pefa()]
#'
#' @export
ssl <- function(x) {
  if (!inherits(x, "pefa")) {
    stop("`x` must be a \"pefa\" object returned by pefa().", call. = FALSE)
  }
  loadings <- x[["loadings"]]
  if (!is.list(loadings)) {
    stop("`x$loadings` must be the named list of candidate loading matrices.",
         call. = FALSE)
  }
  lapply(loadings, function(m) colSums(m^2))
}


# ---- private display helpers ------------------------------------------

.pefa_disp_window <- function(x) {
  K <- x[["sweep"]][["K"]]
  c(Kmin = as.integer(min(K)), Kmax = as.integer(max(K)),
    K0 = as.integer(ncol(x[["Q0"]])))
}

## which() rather than a logical index so a defensive NA cannot silently
## widen the reported set.
.pefa_disp_nonconverged <- function(sweep) {
  as.integer(sweep[["K"]][which(!sweep[["converged"]])])
}

.pefa_disp_mode <- function(settings) {
  if (isTRUE(settings[["bifactor"]])) "bifactor" else "ordinary"
}

.pefa_disp_general <- function(settings) {
  if (isTRUE(settings[["bifactor"]])) " group factors (+1 general)" else ""
}

## Display copy only.  The stored numeric and logical triangles keep their
## own storage modes, so a "*" marker can never leak back into the object.
.pefa_disp_tri <- function(values, collision, digits) {
  out <- matrix("", nrow(values), ncol(values),
                dimnames = dimnames(values))
  upper <- upper.tri(values)
  if (any(upper)) {
    flag <- collision[upper]
    out[upper] <- paste0(formatC(values[upper], format = "f",
                                 digits = digits),
                         ifelse(!is.na(flag) & flag, "*", ""))
  }
  out
}

.pefa_disp_persistence <- function(persistence, digits) {
  phi <- persistence[["phi"]]
  rmsd <- persistence[["rmsd"]]
  collision <- persistence[["collision"]]
  cat("\nPersistence (direct K_from x K_to comparisons):\n")
  cat("  phi (phi_min):\n")
  print(.pefa_disp_tri(phi, collision, digits), quote = FALSE, right = TRUE)
  cat("  rmsd (rmsd_max):\n")
  print(.pefa_disp_tri(rmsd, collision, digits), quote = FALSE, right = TRUE)
  ## One legend line, and only when the reader can actually see a marker.
  if (any(upper.tri(collision) & !is.na(collision) & collision)) {
    cat("  * = target reuse (collision)\n")
  }
  invisible(NULL)
}

## A single K gives range() a degenerate interval, which is not a usable
## plotting window.
.pefa_disp_xlim <- function(K) {
  if (length(K) == 1L) return(c(K - 0.5, K + 0.5))
  range(K)
}

.pefa_disp_range <- function(values, top = 0.20) {
  values <- values[is.finite(values)]
  if (!length(values)) return(c(0, 1))
  r <- range(values)
  span <- diff(r)
  if (span == 0) span <- abs(r[1L]) * 0.01 + 1
  c(r[1L] - 0.05 * span, r[2L] + top * span)
}

## The formal default carries every objective criterion, so a gain view has
## to narrow an untouched default rather than reject it.
.pefa_disp_criterion <- function(criterion, type) {
  allowed <- if (identical(type, "gain")) {
    c("ELBO", "BIC")
  } else c("ELBO", "AIC", "BIC")
  if (!is.character(criterion) || !length(criterion) || anyNA(criterion)) {
    stop("criterion must be a criterion name.", call. = FALSE)
  }
  chosen <- toupper(criterion)
  if (length(chosen) > 1L) {
    if (!identical(chosen, c("ELBO", "AIC", "BIC"))) {
      stop("criterion must be a single criterion name.", call. = FALSE)
    }
    return(allowed[[1L]])
  }
  if (!(chosen %in% allowed)) {
    stop(sprintf(
      "criterion = \"%s\" is unavailable for type = \"%s\"; use one of %s.",
      criterion, type, paste0("\"", allowed, "\"", collapse = ", ")),
      call. = FALSE)
  }
  chosen
}


# ---- methods ----------------------------------------------------------

#' Methods for a PEFA factor-count sweep
#'
#' Descriptive readers for the object returned by [pefa()]. They report what
#' was fitted and measured and nothing else: no method selects, recommends, or
#' marks a factor count, and none draws a threshold, because the package owns
#' no count rule. Count rules, persistence thresholds, and collision policies
#' belong to the analysis that consumes this object.
#'
#' @param object,x A result returned by [pefa()]; for
#'   `print.summary.pefa()`, the object returned by `summary()`.
#' @param ... Further arguments; for `plot()`, graphical arguments passed to
#'   the first base plot in the selected display.
#'
#' @return `print()` and `plot()` return their input invisibly.
#'   `summary()` returns a list of class `"summary.pefa"` with exactly the
#'   components `window` (named integer `Kmin`, `Kmax`, `K0`), `sweep`,
#'   `transitions`, `persistence`, `ssl`, `settings`, and `nonconverged_K`
#'   (an integer vector of the K values whose candidate did not converge,
#'   `integer(0)` when all converged).
#'
#' @seealso [pefa()], [ssl()]
#'
#' @rdname pefa-methods
#' @export
print.pefa <- function(x, ...) {
  w <- .pefa_disp_window(x)
  settings <- x[["settings"]]
  cat(sprintf("PEFA sweep (%s): K = %d:%d%s, K0 = %d\n",
              .pefa_disp_mode(settings), w[["Kmin"]], w[["Kmax"]],
              .pefa_disp_general(settings), w[["K0"]]))
  phi <- x[["persistence"]][["phi"]]
  cat(sprintf(
    "  sweep %d x %d | transitions %d x %d | persistence 3 x %d x %d\n",
    nrow(x[["sweep"]]), ncol(x[["sweep"]]),
    nrow(x[["transitions"]]), ncol(x[["transitions"]]),
    nrow(phi), ncol(phi)))
  bad <- .pefa_disp_nonconverged(x[["sweep"]])
  if (length(bad)) {
    cat("  nonconverged K: ", paste(bad, collapse = ", "), "\n", sep = "")
  } else {
    cat("  all candidates converged\n")
  }
  invisible(x)
}

#' @rdname pefa-methods
#' @export
summary.pefa <- function(object, ...) {
  out <- list(
    window = .pefa_disp_window(object),
    sweep = object[["sweep"]],
    transitions = object[["transitions"]],
    persistence = object[["persistence"]],
    ssl = ssl(object),
    settings = object[["settings"]],
    nonconverged_K = .pefa_disp_nonconverged(object[["sweep"]])
  )
  class(out) <- "summary.pefa"
  out
}

#' @param digits Number of decimals used when printing measured quantities.
#'
#' @rdname pefa-methods
#' @export
print.summary.pefa <- function(x, digits = 3, ...) {
  settings <- x[["settings"]]
  w <- x[["window"]]
  cat(sprintf("PEFA sweep summary (%s): K = %d:%d%s, K0 = %d\n",
              .pefa_disp_mode(settings), w[["Kmin"]], w[["Kmax"]],
              .pefa_disp_general(settings), w[["K0"]]))

  cat("\nCandidates:\n")
  candidates <- x[["sweep"]]
  ## Criteria and indices live on different scales; rounding each to its own
  ## readable width keeps the table narrow without hiding a real difference.
  for (nm in intersect(c("ELBO", "AIC", "BIC"), names(candidates))) {
    candidates[[nm]] <- round(candidates[[nm]], 1L)
  }
  for (nm in intersect(c("RMSEA", "SRMR", "CFI", "TLI"),
                       names(candidates))) {
    candidates[[nm]] <- round(candidates[[nm]], digits)
  }
  print(candidates, row.names = FALSE)

  cat("\nSS loadings (colSums(Lam^2)",
      if (isTRUE(settings[["bifactor"]])) "; general column first" else "",
      "):\n", sep = "")
  for (nm in names(x[["ssl"]])) {
    cat(sprintf("  K=%s: %s\n", nm,
                paste(sprintf("%.2f", x[["ssl"]][[nm]]), collapse = " ")))
  }

  cat("\nTransitions:\n")
  transitions <- x[["transitions"]]
  if (nrow(transitions)) {
    num <- vapply(transitions, is.numeric, logical(1))
    transitions[num] <- lapply(transitions[num], round, digits)
    print(transitions, row.names = FALSE)
  } else {
    cat("  none: the window holds a single candidate\n")
  }

  .pefa_disp_persistence(x[["persistence"]], digits)
  invisible(x)
}

#' @param type Display: criterion trajectories (`"objective"`), one gain
#'   trajectory (`"gain"`), or descriptive absolute and incremental fit
#'   trajectories (`"fit"`).
#' @param criterion Criterion path. `"ELBO"`, `"AIC"`, or `"BIC"` for
#'   `type = "objective"`; `"ELBO"` or `"BIC"` for `type = "gain"`, which has
#'   no AIC path. Ignored when `type = "fit"`.
#' @param pct Gain view only. `FALSE` (the default) computes the raw oriented
#'   gains `ELBO_to - ELBO_from` or `BIC_from - BIC_to` from `$sweep`;
#'   `TRUE` reads the stored `ELBO_gain_pct` or `BIC_gain_pct` column of
#'   `$transitions`, which is `NA` for a path with no positive gain.
#'
#' @rdname pefa-methods
#' @importFrom graphics axis legend lines mtext par plot
#' @export
plot.pefa <- function(x, type = c("objective", "gain", "fit"),
                      criterion = c("ELBO", "AIC", "BIC"),
                      pct = FALSE, ...) {
  type <- match.arg(type)
  ## The fit view has no criterion path, so it ignores the argument outright
  ## rather than rejecting a value it would not have read.
  if (type != "fit") criterion <- .pefa_disp_criterion(criterion, type)
  if (length(pct) != 1L || !is.logical(pct) || is.na(pct)) {
    stop("pct must be TRUE or FALSE.", call. = FALSE)
  }
  d <- x[["sweep"]]
  K <- d[["K"]]
  col1 <- "#1f77b4"
  col2 <- "#d62728"
  xlim <- .pefa_disp_xlim(K)
  xlab <- if (isTRUE(x[["settings"]][["bifactor"]])) {
    "Number of group factors (K)"
  } else "Number of factors (K)"

  if (type == "objective") {
    better <- if (criterion == "ELBO") {
      "higher is better"
    } else "lower is better"
    plot(K, d[[criterion]], type = "b", pch = 16, col = col1, xlim = xlim,
         ylim = .pefa_disp_range(d[[criterion]]), xaxt = "n", xlab = xlab,
         ylab = paste0(criterion, " (", better, ")"),
         main = paste("PEFA", criterion, "path"), ...)
    axis(1, at = K)
  } else if (type == "gain") {
    ylab <- if (pct) {
      "Gain (% of largest positive gain)"
    } else "Oriented gain (positive favors K_to)"
    main <- paste("PEFA", criterion, if (pct) "gain (%)" else "gain", "path")
    if (nrow(d) < 2L) {
      ## A one-candidate window has no adjacent pair.  That is a legal
      ## window, so the gain view stays an empty panel instead of erroring.
      plot(K, rep(NA_real_, length(K)), type = "n", xlim = xlim,
           ylim = c(0, 1), xaxt = "n", xlab = xlab, ylab = ylab,
           main = main, ...)
      axis(1, at = K)
      mtext("no adjacent pair in a one-candidate window", side = 3,
            line = 0.2, cex = 0.9)
    } else {
      transitions <- x[["transitions"]]
      gains <- if (pct) {
        transitions[[paste0(criterion, "_gain_pct")]]
      } else if (criterion == "ELBO") {
        diff(d[["ELBO"]])
      } else -diff(d[["BIC"]])
      at <- if (pct) transitions[["K_to"]] else K[-1L]
      ## No reference line is drawn anywhere in this method, so zero is held
      ## inside the window instead: the sign of a gain is what a reader needs.
      plot(at, gains, type = "b", pch = 16, col = col1, xlim = xlim,
           ylim = .pefa_disp_range(c(gains, 0)), xaxt = "n", xlab = xlab,
           ylab = ylab, main = main, ...)
      axis(1, at = K)
    }
  } else {
    op <- par(mfrow = c(1, 2), mar = c(5, 4.2, 3, 1) + 0.1,
              oma = c(0, 0, 2, 0))
    on.exit(par(op))
    plot(K, d[["RMSEA"]], type = "b", pch = 16, col = col1, xlim = xlim,
         ylim = .pefa_disp_range(c(d[["RMSEA"]], d[["SRMR"]])), xaxt = "n",
         xlab = xlab, ylab = "Index (lower is better)",
         main = "Absolute misfit", ...)
    axis(1, at = K)
    lines(K, d[["SRMR"]], type = "b", pch = 17, lty = 2, col = col2)
    legend("top", c("RMSEA", "SRMR"), col = c(col1, col2), lty = c(1, 2),
           pch = c(16, 17), bty = "n", cex = 0.9)

    plot(K, d[["CFI"]], type = "b", pch = 16, col = col1, xlim = xlim,
         ylim = .pefa_disp_range(c(d[["CFI"]], d[["TLI"]])), xaxt = "n",
         xlab = xlab, ylab = "Index (higher is better)",
         main = "Incremental fit")
    axis(1, at = K)
    lines(K, d[["TLI"]], type = "b", pch = 17, lty = 2, col = col2)
    legend("bottom", c("CFI", "TLI"), col = c(col1, col2), lty = c(1, 2),
           pch = c(16, 17), bty = "n", cex = 0.9)
    mtext("PEFA descriptive fit trajectories", outer = TRUE, font = 2)
  }
  invisible(x)
}
