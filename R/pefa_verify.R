## Public integrity verification for pefa() evidence objects (plan section
## 3.4).  This is the one new audit API in 0.9.0: a consuming project must be
## able to check a package-issued evidence object without reproducing package
## internals, and without gaining an entry point that hands back rederived
## evidence.  Everything rebuilt or rehashed here therefore stays inside the
## function -- the only outward signal is invisible(TRUE) or a controlled
## error naming the component that failed.
##
## Nothing in this file hashes anything itself.  Every fingerprint is
## recomputed by calling the same canonical helpers in R/pefa_provenance.R
## that pefa() calls when it seals the object, so a verification can never
## drift from the sealing code by one payload field, one dropped `secs`, or
## one canonicalization rule.

## Top-level components of a current-schema evidence object, in construction
## order.  The shape is checked directly and exactly because evidence_sha256
## cannot see it: the hashed payload is assembled by name, so an appended
## component -- a smuggled `$selected_K`, say -- would leave every hash intact.
.pefa_verify_components <- c("sweep", "transitions", "persistence", "loadings",
                             "pips", "Q0", "settings", "provenance")

## The section 3.6 provenance fields, in order.  This restates the schema
## rather than deriving it from .pefa_provenance(), on purpose: an independent
## restatement is what makes the verifier a check instead of a mirror. Any
## change to the field list is a schema change and must bump
## .pefa_evidence_schema_id, which is gated above this check, so the two can
## only ever disagree during a botched edit -- exactly when a mirror would
## stay silent.
.pefa_verify_provenance_fields <- c(
  "evidence_schema_id", "engine_semantics_id", "package_version",
  "package_build_id", "standardized_data_sha256", "Q0_sha256",
  "settings_sha256", "evidence_sha256", "lineage_sha256",
  "parent_evidence_sha256", "reused_K", "newly_fitted_K"
)

## Every refusal leaves through here, so a failure always names the function
## and never leaks the internal call that raised it.
.pefa_verify_stop <- function(...) {
  stop("verify_pefa(): ", ..., call. = FALSE)
}

.pefa_verify_names <- function(nm) {
  paste0("`$", nm, "`", collapse = ", ")
}

## Render a stored identifier for a message without pasting an arbitrary
## object into it: a tampered field can hold anything at all.
.pefa_verify_label <- function(value) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    return(paste0("\"", value, "\""))
  }
  paste0("<", class(value)[1L], " of length ", length(value), ">")
}

.pefa_verify_hash_field <- function(prov, field) {
  value <- prov[[field]]
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    .pefa_verify_stop("`$provenance$", field, "` must be one nonmissing, ",
                      "nonempty character string.")
  }
  value
}

## Compare a stored evidence table against the one the installed builders
## produce.  Returns NULL when they agree and a one-clause description of the
## first disagreement otherwise.
##
## The policy is the one section 6.1 fixes for evidence tables and ?verify_pefa
## states for users: identical() for everything that is not a double column,
## and all.equal(tolerance = 1e-10) -- the comparison testthat::expect_equal()
## performs -- for the double cells.  The NA *and* NaN patterns are compared
## with identical() first, because all.equal() treats an NA_real_ as matching a
## NaN and an integrity checker must not.
.pefa_verify_table <- function(stored, rebuilt) {
  if (!is.data.frame(stored)) return("the stored value is not a data frame")
  if (!identical(class(stored), class(rebuilt))) return("the class differs")
  if (!identical(names(stored), names(rebuilt))) {
    return("the column names or their order differ")
  }
  if (!identical(attr(stored, "row.names"), attr(rebuilt, "row.names"))) {
    return("the row count or row names differ")
  }
  if (!identical(vapply(stored, typeof, ""), vapply(rebuilt, typeof, ""))) {
    return("a column storage mode differs")
  }
  for (nm in names(rebuilt)) {
    a <- stored[[nm]]
    b <- rebuilt[[nm]]
    if (identical(typeof(b), "double")) {
      if (!identical(attributes(a), attributes(b))) {
        return(sprintf("column `%s` carries different attributes", nm))
      }
      if (!identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b))) {
        return(sprintf("column `%s` has a different NA pattern", nm))
      }
      cmp <- all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE)
      if (!isTRUE(cmp)) {
        return(sprintf("column `%s` differs beyond tolerance 1e-10 (%s)",
                       nm, paste(cmp, collapse = "; ")))
      }
    } else if (!identical(a, b)) {
      return(sprintf("column `%s` is not identical to the rebuilt column", nm))
    }
  }
  NULL
}

#' Verify the integrity of a PEFA evidence object
#'
#' The sole public integrity checker for objects returned by [pefa()]. It
#' accepts the current evidence schema and engine semantics only, rebuilds
#' `$transitions` and `$persistence` from the stored primitive components with
#' the installed package, and recomputes the Q0, settings, evidence, and
#' lineage fingerprints. Supplying the original data additionally verifies the
#' standardized-data fingerprint. It returns `invisible(TRUE)` when every
#' applicable check passes and otherwise raises an error naming the component
#' that failed.
#'
#' @param x An object of class `"pefa"` returned by [pefa()].
#' @param Y Optional. The original `N` by `J` data matrix the object was
#'   fitted to. When supplied, `standardized_data_sha256` is recomputed from
#'   it through the same internal standardization [pefa()] applies before
#'   fitting. When `NULL` (the default) the stored data fingerprint is checked
#'   only as part of the evidence payload; see Details for exactly what that
#'   does and does not prove.
#'
#' @details
#' `verify_pefa()` answers one question: is this object still exactly the
#' evidence the installed package sealed? It runs the checks below and reports
#' every one that fails, so a single call locates the damage rather than
#' stopping at the first symptom.
#'
#' **Compatibility.** `$provenance$evidence_schema_id` and
#' `$provenance$engine_semantics_id` must equal the identifiers frozen in this
#' installation. An object written under any other schema or engine is refused
#' outright rather than checked, because neither its payload nor its hash
#' canonicalization is this package's. `package_version` and
#' `package_build_id` are not compatibility gates and are not required to
#' match the running installation.
#'
#' **Rebuilt tables.** `$transitions` is rebuilt from `$sweep`, and
#' `$persistence` from `$sweep`, `$loadings`, `$pips`, `$Q0`, and the
#' `stability_eps` and `bifactor` settings, using the installed table builders
#' and matcher. [vbfa()] is never called and nothing is refitted. The rebuilt
#' tables are used for comparison only and are never returned.
#'
#' **Fingerprints.** `Q0_sha256` is recomputed from `$Q0`, `settings_sha256`
#' from `$settings`, `evidence_sha256` from the stored evidence payload, and
#' `lineage_sha256` from the stored evidence hash together with the
#' parent/reused/newly-fitted declaration. Elapsed `$sweep$secs` lies outside
#' the evidence payload by design, so the same analysis timed on a slower
#' machine is not a mutation. Because `lineage_sha256` is a function of
#' `evidence_sha256`, it is checked only once that hash agrees.
#'
#' **Comparison policy.** Stored and rebuilt tables are compared under the
#' policy the package fixes for its evidence tables, stated here so that
#' downstream reporting can rely on it: `identical()` for integers, logicals,
#' status/reason character vectors, `r`, dimnames, storage modes, and the NA
#' pattern; `expect_equal(tolerance = 1e-10)` for double metric cells. The
#' implementation applies base R's `all.equal(tolerance = 1e-10)` to the
#' double columns, which is the comparison `testthat::expect_equal()`
#' performs, and `identical()` to column names and order, row names,
#' per-column storage modes, every non-double column, and the NA and NaN
#' patterns of the double ones. Fingerprints are never compared at a
#' tolerance: they are compared with `identical()`.
#'
#' **What `Y = NULL` does and does not prove.**
#' `standardized_data_sha256` is one of the fields covered by
#' `evidence_sha256`. With `Y = NULL` the stored data fingerprint is therefore
#' verified to be the value that was sealed with the rest of the payload:
#' nobody has swapped the recorded data identity while leaving the evidence
#' intact. That is all it proves. No preimage is reconstructed, no data set is
#' examined, and no claim is made about which data the object was fitted to.
#' Supplying `Y` adds exactly that claim and only that claim: the supplied
#' matrix, standardized the way [pefa()] standardizes its input, hashes to the
#' stored value, so it is the fitted input in the same row and column order
#' and with the same missingness pattern.
#'
#' The function makes no scientific decision. It does not judge convergence,
#' persistence, or factor count, and a passing result says nothing about
#' whether the evidence supports any conclusion.
#'
#' @return `invisible(TRUE)` when every applicable check passes. There is no
#'   other return value: the rebuilt tables and recomputed fingerprints stay
#'   inside the call, since 0.9.0 exports nothing that returns rederived
#'   evidence. Any failure is a controlled error naming the component, never a
#'   `FALSE`. `x` is never modified.
#'
#' @examples
#' ## The shipped R8 integration fixture carries a complete, classed evidence
#' ## object together with the exact data that produced its standardized-data
#' ## fingerprint, so both forms of the check run without fitting anything.
#' path <- system.file("extdata", "vbpm_r8_fixture.rds", package = "vbpm")
#' if (nzchar(path) && file.exists(path)) {
#'   integrity <- readRDS(path)[["integrity"]]
#'   if (!is.null(integrity)) {
#'     verify_pefa(integrity[["object"]])                    # TRUE, invisibly
#'     verify_pefa(integrity[["object"]], integrity[["Y"]])  # also checks Y
#'   }
#' }
#'
#' @seealso [pefa()], [persistence()], [ssl()]
#'
#' @export
verify_pefa <- function(x, Y = NULL) {
  ## ---- object shape --------------------------------------------------
  if (!inherits(x, "pefa") || !is.list(x)) {
    .pefa_verify_stop("`x` must be a \"pefa\" object returned by pefa().")
  }
  if (!identical(names(x), .pefa_verify_components)) {
    absent <- .pefa_verify_components[!.pefa_verify_components %in% names(x)]
    extra <- names(x)[!names(x) %in% .pefa_verify_components]
    .pefa_verify_stop(
      "`x` does not have the current evidence-object shape",
      if (length(absent)) {
        paste0("; missing component(s) ", .pefa_verify_names(absent))
      } else "",
      if (length(extra)) {
        paste0("; unexpected component(s) ", .pefa_verify_names(extra))
      } else "",
      if (!length(absent) && !length(extra)) {
        "; the components are out of order"
      } else "",
      ".",
      if ("provenance" %in% absent) {
        paste0(" An object carrying no `$provenance` predates evidence ",
               "schema \"", .pefa_evidence_schema_id, "\".")
      } else ""
    )
  }

  ## ---- compatibility -------------------------------------------------
  prov <- x[["provenance"]]
  if (!is.list(prov) ||
      !identical(names(prov), .pefa_verify_provenance_fields)) {
    .pefa_verify_stop(
      "`$provenance` does not carry the twelve fields of evidence schema \"",
      .pefa_evidence_schema_id, "\" in order, so the object cannot be ",
      "verified.")
  }
  if (!identical(prov[["evidence_schema_id"]], .pefa_evidence_schema_id)) {
    .pefa_verify_stop(
      "`$provenance$evidence_schema_id` is ",
      .pefa_verify_label(prov[["evidence_schema_id"]]),
      "; this installation verifies only \"", .pefa_evidence_schema_id,
      "\" objects.")
  }
  if (!identical(prov[["engine_semantics_id"]], .pefa_engine_semantics_id)) {
    .pefa_verify_stop(
      "`$provenance$engine_semantics_id` is ",
      .pefa_verify_label(prov[["engine_semantics_id"]]),
      "; this installation verifies only \"", .pefa_engine_semantics_id,
      "\" objects.")
  }

  ## ---- provenance field types ----------------------------------------
  ## Checked before anything is recomputed so that a damaged record is named
  ## as such, rather than surfacing as a hash mismatch it did not cause.
  for (field in c("package_version", "standardized_data_sha256", "Q0_sha256",
                  "settings_sha256", "evidence_sha256", "lineage_sha256")) {
    .pefa_verify_hash_field(prov, field)
  }
  parent <- prov[["parent_evidence_sha256"]]
  if (!is.character(parent) || length(parent) != 1L ||
      (!is.na(parent) && !nzchar(parent))) {
    .pefa_verify_stop("`$provenance$parent_evidence_sha256` must be one ",
                      "character string, or NA_character_ for a root object.")
  }
  for (field in c("reused_K", "newly_fitted_K")) {
    value <- prov[[field]]
    if (!is.integer(value) || anyNA(value)) {
      .pefa_verify_stop("`$provenance$", field,
                        "` must be an integer vector holding no NA.")
    }
  }

  ## ---- components the rebuild consumes -------------------------------
  sweep <- x[["sweep"]]
  if (!is.data.frame(sweep) || is.null(sweep[["K"]])) {
    .pefa_verify_stop("`$sweep` must be a data frame with a `K` column.")
  }
  for (field in c("transitions", "persistence")) {
    if (!is.data.frame(x[[field]])) {
      .pefa_verify_stop("`$", field, "` must be a data frame.")
    }
  }
  for (field in c("loadings", "pips")) {
    if (!is.list(x[[field]]) || is.null(names(x[[field]]))) {
      .pefa_verify_stop("`$", field, "` must be a list keyed by `K`.")
    }
  }
  Q0 <- x[["Q0"]]
  if (!is.matrix(Q0)) {
    .pefa_verify_stop("`$Q0` must be a matrix.")
  }
  settings <- x[["settings"]]
  if (!is.list(settings)) {
    .pefa_verify_stop("`$settings` must be a list.")
  }
  needed <- c("bifactor", "stability_eps")
  absent <- needed[!needed %in% names(settings)]
  if (length(absent)) {
    .pefa_verify_stop("`$settings` is missing ",
                      paste0("`", absent, "`", collapse = ", "),
                      ", which the rebuild needs.")
  }

  ## An unusable `Y` is a call problem, not an integrity finding, so it stops
  ## here instead of joining the report below.
  data_hash <- NULL
  if (!is.null(Y)) {
    data_hash <- tryCatch(.pefa_standardized_hash(Y), error = identity)
    if (inherits(data_hash, "error")) {
      .pefa_verify_stop("`Y` cannot be standardized the way pefa() ",
                        "standardizes its input: ",
                        conditionMessage(data_hash))
    }
  }

  ## ---- integrity checks ----------------------------------------------
  ## Collected rather than raised one at a time: a reader repairing an archive
  ## needs to see every component that moved, and the more specific findings
  ## above narrow the less specific ones below.
  failures <- character(0)
  add <- function(...) failures <<- c(failures, paste0(...))

  if (!identical(.pefa_sha256(Q0), prov[["Q0_sha256"]])) {
    add("`$Q0` does not reproduce `$provenance$Q0_sha256`.")
  }
  if (!identical(.pefa_sha256(settings), prov[["settings_sha256"]])) {
    add("`$settings` does not reproduce `$provenance$settings_sha256`.")
  }

  rebuilt <- tryCatch(.pefa_transitions_table(sweep), error = identity)
  if (inherits(rebuilt, "error")) {
    add("`$transitions` could not be rebuilt from `$sweep`: ",
        conditionMessage(rebuilt), ".")
  } else {
    difference <- .pefa_verify_table(x[["transitions"]], rebuilt)
    if (!is.null(difference)) {
      add("`$transitions` disagrees with the table rebuilt from `$sweep`: ",
          difference, ".")
    }
  }

  rebuilt <- tryCatch(
    .pefa_persistence_table(sweep, x[["loadings"]], x[["pips"]], Q0,
                            ncol(Q0), settings[["stability_eps"]],
                            settings[["bifactor"]]),
    error = identity
  )
  if (inherits(rebuilt, "error")) {
    add("`$persistence` could not be rebuilt from `$sweep`, `$loadings`, ",
        "`$pips`, `$Q0`, and `$settings`: ", conditionMessage(rebuilt), ".")
  } else {
    difference <- .pefa_verify_table(x[["persistence"]], rebuilt)
    if (!is.null(difference)) {
      add("`$persistence` disagrees with the table rebuilt from `$sweep`, ",
          "`$loadings`, `$pips`, `$Q0`, and `$settings`: ", difference, ".")
    }
  }

  ## The recomputation runs through .pefa_provenance() itself, so the payload,
  ## its field order, the `secs` exclusion, and the canonicalization are the
  ## sealing code's rather than a copy of it.
  recomputed <- tryCatch(
    .pefa_provenance(
      evidence = list(sweep = sweep, transitions = x[["transitions"]],
                      persistence = x[["persistence"]],
                      loadings = x[["loadings"]], pips = x[["pips"]],
                      Q0 = Q0, settings = settings),
      standardized_data_sha256 = prov[["standardized_data_sha256"]],
      Q0_sha256 = prov[["Q0_sha256"]],
      settings_sha256 = prov[["settings_sha256"]],
      package_version = prov[["package_version"]],
      newly_fitted_K = prov[["newly_fitted_K"]],
      reused_K = prov[["reused_K"]],
      parent_evidence_sha256 = prov[["parent_evidence_sha256"]]
    ),
    error = identity
  )
  if (inherits(recomputed, "error")) {
    add("`$provenance` could not be recomputed from the stored payload: ",
        conditionMessage(recomputed), ".")
  } else if (!identical(recomputed[["evidence_sha256"]],
                        prov[["evidence_sha256"]])) {
    add("`$provenance$evidence_sha256` does not match the hash recomputed ",
        "from the stored payload: one of `$sweep` (elapsed `secs` excluded), ",
        "`$transitions`, `$persistence`, `$loadings`, `$pips`, `$Q0`, ",
        "`$settings`, `$provenance$package_version`, and ",
        "`$provenance$standardized_data_sha256` has changed since the object ",
        "was sealed. `$provenance$lineage_sha256` cannot be checked ",
        "independently until it agrees.")
  } else if (!identical(recomputed[["lineage_sha256"]],
                        prov[["lineage_sha256"]])) {
    add("`$provenance$lineage_sha256` does not match the hash recomputed ",
        "from the stored declaration: `$provenance$parent_evidence_sha256`, ",
        "`$provenance$reused_K`, or `$provenance$newly_fitted_K` has changed.")
  }

  if (!is.null(data_hash) &&
      !identical(data_hash, prov[["standardized_data_sha256"]])) {
    add("`Y` does not reproduce `$provenance$standardized_data_sha256`: the ",
        "supplied matrix is not the input this object was fitted to, or is ",
        "not in the same row and column order.")
  }

  if (length(failures)) {
    .pefa_verify_stop(
      "the evidence object failed ", length(failures),
      if (length(failures) == 1L) " integrity check:\n" else
        " integrity checks:\n",
      paste0("  - ", failures, collapse = "\n")
    )
  }
  invisible(TRUE)
}
