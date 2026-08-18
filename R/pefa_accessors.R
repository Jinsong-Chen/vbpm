## Accessors for the evidence object returned by pefa().  Both read stored
## fields only: switching the metric must not be able to disagree with
## $persistence, so nothing here recomputes a matcher quantity.

## Fields a reader plausibly asks the triangle for and will not find there.
## match.arg() would only list the seven legal names, which does not say where
## the thing they asked for actually lives.
.pefa_persistence_elsewhere <- c(
  unmatched_n = "persistence",
  source_ineligible_n = "persistence",
  target_eligible_n = "persistence",
  n_pairs = "persistence",
  backbone_ssl_min = "sweep"
)

#' Persistence triangle for a PEFA sweep
#'
#' Pivots the long-format `$persistence` table of a [pefa()] object into an
#' upper-triangular `W x W` matrix holding one stored metric: rows are the
#' source count `k`, columns the target count `l`, and both dimnames are
#' `as.character(x$sweep$K)`. The diagonal and the lower half are `NA`.
#' Nothing is recomputed, so switching `metric` is free and cannot disagree
#' with the table.
#'
#' @param x A result returned by [pefa()]; for the print method, the matrix
#'   `persistence()` returned.
#' @param metric Name of the stored pair metric to pivot. One of `"phi_min"`,
#'   `"rmsd"`, `"rmsd_max"`, `"ari"`, `"pip_rmsd"`, `"unmatched_ssl"`, or
#'   `"n_collisions"`. The enumeration is closed; see Details.
#' @param digits Number of decimals used when printing the triangle.
#' @param ... Further arguments. Unused.
#'
#' @details
#' A cell is `NA` when the selected metric itself is unavailable, whether
#' because the matcher input was unavailable or because that metric's own
#' narrower dependencies were unmet. Collided, nonconverged-endpoint, and
#' reduced-pair-set cells retain any computed value.
#'
#' The `metric` enumeration is **closed** and deliberately excludes
#' `unmatched_n`, `source_ineligible_n`, `target_eligible_n`, and `n_pairs`:
#' those are matching-process counts rather than structural measurements, and
#' a `W x W` triangle is the wrong shape for reading them. They stay in
#' `$persistence`, which is already long-format and is the natural place to
#' filter or join on them. `backbone_ssl_min` is likewise absent because it is
#' a per-candidate quantity on `$sweep`, not a pair metric -- a triangle of it
#' would just repeat one vector twice. Asking for any of these five is an
#' error that names the table the field lives in.
#'
#' Cells are **not on a common pair-set scale.** `phi_min` in row `k` is a
#' minimum over at most `k` pairs, so rows are not directly comparable, and
#' within a row a `~` cell scored fewer columns than its neighbours. The
#' `k = K0` row is backbone-only throughout, and by construction carries
#' **no** `~`: its pair set is complete for that source
#' (`n_pairs == K_from == K0`), it is simply a small source. A high value
#' there says the anchors did not move, not that an exploratory structure
#' persisted. `n_pairs` gives the exact count.
#'
#' The printed `?` marks **selected-metric unavailability only**. It is not
#' the complete tri-state edge state, which a decision profile computes from
#' `pair_status`, the chosen metric's own `NA`, and its declared endpoint,
#' collision, and eligibility conditions.
#'
#' The three markers get equal display weight on purpose: each is a fact that
#' changes how the printed number should be read, and singling one out for a
#' glyph would reintroduce a package preference through the print method.
#'
#' @return `persistence()` returns a `W x W` matrix of class
#'   `c("vbpm_persistence", "matrix", "array")`, integer for `"n_collisions"`
#'   and double otherwise. It carries three logical matrices of the same `dim`
#'   and `dimnames` as attributes: `"collision"` (the stored collision flag),
#'   `"reduced"` (`n_pairs < K_from`), and `"fit_ok"`
#'   (`from_fit_ok & to_fit_ok`). Each is `NA` only where its own underlying
#'   factual field is unavailable. `print()` returns its input invisibly.
#'
#' @seealso [pefa()], [ssl()]
#'
#' @export
persistence <- function(x, metric = c("phi_min", "rmsd", "rmsd_max", "ari",
                                      "pip_rmsd", "unmatched_ssl",
                                      "n_collisions")) {
  if (!inherits(x, "pefa")) {
    stop("`x` must be a \"pefa\" object returned by pefa().", call. = FALSE)
  }
  if (is.character(metric) && length(metric) == 1L && !is.na(metric) &&
      metric %in% names(.pefa_persistence_elsewhere)) {
    field <- .pefa_persistence_elsewhere[[metric]]
    stop(sprintf(
      paste0("metric = \"%s\" is not a persistence-triangle metric: it is %s. ",
             "Read it from `x$%s`, e.g. `x$%s[[\"%s\"]]`."),
      metric,
      if (identical(field, "sweep")) {
        "a per-candidate quantity"
      } else {
        "a matching-process count"
      },
      field, field, metric
    ), call. = FALSE)
  }
  metric <- match.arg(metric)

  sweep <- x[["sweep"]]
  pairs <- x[["persistence"]]
  if (!is.data.frame(sweep) || is.null(sweep[["K"]])) {
    stop("`x$sweep` must be a data frame with a `K` column.", call. = FALSE)
  }
  if (!is.data.frame(pairs)) {
    stop("`x$persistence` must be a data frame.", call. = FALSE)
  }
  needed <- c("K_from", "K_to", metric, "collision", "n_pairs",
              "from_fit_ok", "to_fit_ok")
  absent <- needed[!(needed %in% names(pairs))]
  if (length(absent)) {
    stop(sprintf("`x$persistence` is missing the column%s %s.",
                 if (length(absent) == 1L) "" else "s",
                 paste0("`", absent, "`", collapse = ", ")),
         call. = FALSE)
  }

  labels <- as.character(sweep[["K"]])
  W <- length(labels)
  ## An integer matrix cannot hold NA_real_, so the unusable half takes the NA
  ## flavor of the storage mode the metric itself demands.
  count_metric <- identical(metric, "n_collisions")
  out <- matrix(if (count_metric) NA_integer_ else NA_real_,
                nrow = W, ncol = W, dimnames = list(labels, labels))

  ## Only the strict upper triangle is addressable; a stored row naming a
  ## candidate outside the sweep window has no cell and is dropped rather than
  ## silently recycled into one.
  from <- match(as.character(pairs[["K_from"]]), labels)
  to <- match(as.character(pairs[["K_to"]]), labels)
  keep <- !is.na(from) & !is.na(to) & from < to
  cell <- cbind(from[keep], to[keep])
  filled <- any(keep)

  values <- pairs[[metric]]
  values <- if (count_metric) as.integer(values) else as.numeric(values)
  if (filled) out[cell] <- values[keep]

  triangle <- function(v) {
    m <- matrix(NA, nrow = W, ncol = W, dimnames = list(labels, labels))
    if (filled) m[cell] <- v[keep]
    m
  }
  ## `&` already gives FALSE where one endpoint is a known nonconverged fit and
  ## the other is unknown, so the pair fact survives a single missing flag.
  attr(out, "collision") <- triangle(as.logical(pairs[["collision"]]))
  attr(out, "reduced") <- triangle(pairs[["n_pairs"]] < pairs[["K_from"]])
  attr(out, "fit_ok") <- triangle(as.logical(pairs[["from_fit_ok"]]) &
                                    as.logical(pairs[["to_fit_ok"]]))
  class(out) <- c("vbpm_persistence", "matrix", "array")
  out
}

#' @rdname persistence
#' @export
print.vbpm_persistence <- function(x, digits = 2, ...) {
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0) {
    stop("`digits` must be one finite nonnegative number.", call. = FALSE)
  }
  values <- unclass(x)
  if (!is.matrix(values)) {
    stop("`x` must be a matrix returned by persistence().", call. = FALSE)
  }
  W <- nrow(values)

  flag <- function(name) {
    m <- attr(x, name, exact = TRUE)
    if (is.logical(m) && identical(dim(m), dim(values))) {
      m
    } else {
      matrix(NA, nrow = W, ncol = ncol(values))
    }
  }
  collision <- flag("collision")
  reduced <- flag("reduced")
  fit_ok <- flag("fit_ok")

  ## An unavailable fact marks nothing; only an established TRUE does.
  established <- function(v) !is.na(v) & v

  cells <- rep("", length(values))
  usable <- which(upper.tri(values))
  shown <- as.vector(values)[usable]
  ## Fixed-decimal formatting keeps the leading zero on a metric like 0.91.
  ## An integer-stored metric is a count, though, and "1.00 collisions" is
  ## just wrong, so counts print as whole numbers.
  text <- if (is.integer(values)) {
    formatC(shown, format = "d")
  } else {
    formatC(shown, format = "f", digits = digits)
  }
  text[is.na(shown)] <- "?"
  ## Markers describe the pair, not the number, so they are appended to an
  ## unavailable cell too whenever their own fact is established.
  collided <- established(collision[usable])
  reduced_pairs <- established(reduced[usable])
  nonconverged <- established(!fit_ok[usable])
  cells[usable] <- paste0(text,
                          ifelse(collided, "*", ""),
                          ifelse(reduced_pairs, "~", ""),
                          ifelse(nonconverged, "!", ""))

  grid <- matrix(cells, nrow = W, ncol = ncol(values),
                 dimnames = dimnames(values))
  print(grid, quote = FALSE, right = TRUE)

  legend <- character(0)
  if (any(collided)) legend <- c(legend, "* = collision")
  if (any(reduced_pairs)) legend <- c(legend, "~ = reduced pair set")
  if (any(nonconverged)) legend <- c(legend, "! = nonconverged endpoint")
  if (any(is.na(shown))) legend <- c(legend, "? = unavailable")
  if (length(legend)) cat(paste(legend, collapse = ", "), "\n", sep = "")

  invisible(x)
}

#' Sums of squared loadings for a PEFA sweep
#'
#' Descriptive column sizes for every candidate in a [pefa()] sweep:
#' `colSums(Lam^2)` on the stored loading matrix of each candidate. They are
#' not Gram eigenvalues and take part in no selection rule. `summary()` calls
#' this function, so there is one definition.
#'
#' @param x A result returned by [pefa()].
#'
#' @return A list with one element per candidate, named as `x$loadings` is
#'   (the candidate `K` as a character string). Element vectors are unnamed.
#'   In bifactor mode each vector covers all `K + 1` stored columns with the
#'   general factor first, so the group-block comparison scale is
#'   `ssl(x)[[k]][-1]`.
#'
#' @seealso [pefa()], [persistence()]
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
  out <- lapply(loadings, function(m) unname(colSums(m^2)))
  names(out) <- names(loadings)
  out
}
