## Immutable provenance for pefa() evidence objects (plan section 3.6).
##
## This is the only file in the package that may use `digest`, because it is the
## only file that hashes anything: canonicalization lives in one place so a
## future change to the byte stream is a single edit and a schema-ID bump.
##
## The two frozen identifiers below are the compatibility contract. Bump
## .pefa_evidence_schema_id for any change to public field names, order, types,
## normalization payloads, or hash canonicalization; bump
## .pefa_engine_semantics_id for any change to standardization, fitting,
## assignment, scored-pair membership, or metric definitions, even when the
## object shape is untouched.

.pefa_evidence_schema_id <- "vbpm_pefa_evidence_1"

.pefa_engine_semantics_id <- "vbpm_pefa_engine_0.9"

## Force every atomic vector in a payload, attributes included, into its
## materialized form. Serialization version 3 writes an ALTREP object in its
## compact representation, so `2:3` and `as.integer(c(2, 3))` -- which are
## `identical()` and mean the same thing -- produce different bytes, as do
## `as.character(1:3)` against a literal character vector and a data frame's
## automatic row names against explicitly assigned ones. Hashing the
## representation instead of the value would report a mutation whenever an
## equivalent payload happened to be built a different way, so the value is
## materialized first and only then serialized.
##
## Recursion is by attribute and list element rather than rapply() because
## attributes carry ALTREP too, and the walk must not convert a pairlist,
## language object, function, environment, or S4 object into something else.
##
## Attributes are also reattached in name order. They are serialized in storage
## order, but `identical()` treats them as a set, so two objects that R calls
## equal -- a data frame whose class was set before its row names against one
## where it was set after -- would otherwise hash apart. Radix order sorts in
## the C locale, so the result does not depend on the session's collation, and
## `dim` sorts before `dimnames`, which is the one ordering constraint R
## enforces on assignment.
##
## Value normalization happens here too, after materialization and before any
## byte ever leaves: IEEE -0 and +0 are `identical()` and `==` yet serialize
## differently, and the same text carried as latin1 rather than UTF-8 is again
## `identical()` (R translates when it compares) yet serializes to different
## bytes under a different declared encoding. Either one would make an
## unchanged scientific payload look mutated. Names and dimnames need no
## special case: both are attributes, so the recursion below feeds them back
## through this same branch at every level.
.pefa_canonical <- function(x) {
  if (isS4(x) || is.function(x) || is.environment(x)) return(x)
  a <- attributes(x)
  if (!is.null(a)) attributes(x) <- NULL
  if (typeof(x) == "list") {
    x <- lapply(x, .pefa_canonical)
  } else if (is.atomic(x) && length(x) > 0L) {
    x[1L] <- x[1L]   # subassignment expands any compact or deferred form
    if (is.double(x)) {
      ## which() drops the NA that NA and NaN return from the comparison, and
      ## +-Inf compares FALSE, so only true signed zeros are rewritten.
      zeros <- which(x == 0)
      if (length(zeros)) x[zeros] <- 0
    } else if (is.character(x)) {
      x <- enc2utf8(x)
    }
  }
  if (!is.null(a)) {
    a <- lapply(a, .pefa_canonical)
    names(a) <- enc2utf8(names(a))
    attributes(x) <- a[order(names(a), method = "radix")]
  }
  x
}

## serialize() prefixes every stream with a header that records the WRITING R
## version (bytes 7-10) and, in version 3, the session's native encoding. Those
## bytes are not part of the object, so hashing them would make an unchanged
## scientific payload hash differently after an R upgrade or on a non-UTF-8
## locale -- defeating the one job evidence_sha256 has, and silently
## invalidating the shipped fixture and any 0.9.1 reuse check. Drop the header
## and hash only the payload.
##
## Version-3 layout: 2 magic bytes, then three 4-byte big-endian integers
## (format version, writer version, minimum reader version), then a 4-byte
## encoding length and that many encoding bytes.
.pefa_payload <- function(x) {
  raw_bytes <- serialize(x, NULL, xdr = TRUE, version = 3L)
  enc_len <- sum(as.integer(raw_bytes[15:18]) * 256^(3:0))
  raw_bytes[-seq_len(18L + enc_len)]
}

## Single canonicalization point for every hash in the package. `xdr = TRUE`
## and `version = 3L` are pinned inside .pefa_payload() so the byte stream
## follows neither the host's endianness nor the session's default
## serialization version; `serialize = FALSE` then tells digest() to hash the
## raw vector it was handed instead of re-serializing it under its own options.
.pefa_sha256 <- function(x) {
  digest::digest(.pefa_payload(.pefa_canonical(x)),
                 algo = "sha256", serialize = FALSE)
}

## Reduce a matrix to values, shape, and names, discarding every other
## attribute -- a class, a caller's bookkeeping tag, or the scaled:center and
## scaled:scale that base scale() leaves behind. Those last two are functions of
## the data the payload already carries, and keeping them would make the
## complete-data branch below hash differently from the incomplete-data branch
## on the same numbers.
.pefa_bare_matrix <- function(m) {
  d <- dim(m)
  dn <- dimnames(m)
  attributes(m) <- NULL
  if (!is.null(d)) dim(m) <- d
  if (!is.null(dn)) dimnames(m) <- dn
  m
}

## Hash the engine-scale data rather than the raw input, so the hash binds what
## vbfa() actually fits: standardized values plus the missingness pattern, in
## the caller's row/column order and under the caller's dimnames. The body below
## mirrors vbfa()'s standardization block line for line, including its split
## between base scale() for complete data -- kept to preserve the historical
## complete-data arithmetic bit for bit -- and the explicit two-sweep()
## transformation otherwise. Any drift here silently breaks reuse verification.
.pefa_standardized_hash <- function(Y) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) && !is.logical(Y)) {
    stop("`Y` must be a real numeric matrix.", call. = FALSE)
  }
  Y <- .pefa_bare_matrix(Y)
  storage.mode(Y) <- "double"
  if (nrow(Y) == 0L || ncol(Y) == 0L) {
    stop("`Y` must have at least one row and one column.", call. = FALSE)
  }

  missing_mask <- is.na(Y)
  has_missing <- any(missing_mask)
  center <- colMeans(Y, na.rm = TRUE)
  scale_ <- apply(Y, 2, stats::sd, na.rm = TRUE)
  if (any(!is.finite(scale_) | scale_ <= 0)) {
    stop("Every item must have a positive observed-data standard deviation.",
         call. = FALSE)
  }
  Z <- if (has_missing) sweep(sweep(Y, 2, center, "-"), 2, scale_, "/") else
       scale(Y)
  Z[missing_mask] <- 0   # observed-scale mean, as in vbfa()

  ## Shape and names are hashed twice over, once inside the matrix and once as
  ## their own entries. The redundancy is deliberate: it fixes the payload
  ## against a future change to how a matrix carries its dimnames.
  .pefa_sha256(list(matrix = .pefa_bare_matrix(Z),
                    missing_mask = .pefa_bare_matrix(missing_mask),
                    dim = as.integer(dim(Z)),
                    dimnames = dimnames(Z)))
}

## Read one installed-DESCRIPTION field, returning NA_character_ for absent,
## missing, empty, or non-scalar values. packageDescription() returns a logical
## NA (with a warning) when the package is not installed, so the is.list() gate
## is load-bearing.
.pefa_desc_field <- function(desc, field) {
  if (!is.list(desc) || !field %in% names(desc)) return(NA_character_)
  value <- desc[[field]]
  if (!is.character(value) && !is.numeric(value)) return(NA_character_)
  if (length(value) != 1L || is.na(value)) return(NA_character_)
  value <- trimws(as.character(value))
  if (!nzchar(value)) return(NA_character_)
  value
}

## Informational build stamp; compatibility is gated by the schema and engine
## identifiers, never by this. Must survive a devtools/pkgload load where the
## installed metadata does not exist at all, hence the fallback chain rather
## than a hard read.
.pefa_build_id <- function(pkg = "vbpm") {
  desc <- tryCatch(suppressWarnings(utils::packageDescription(pkg)),
                   error = function(e) NULL)
  version <- .pefa_desc_field(desc, "Version")
  if (is.na(version)) {
    version <- tryCatch(as.character(utils::packageVersion(pkg)),
                        error = function(e) NA_character_)
  }
  if (is.na(version)) version <- "unknown"

  sha <- .pefa_desc_field(desc, "RemoteSha")
  if (is.na(sha)) sha <- .pefa_desc_field(desc, "GithubSHA1")
  if (!is.na(sha)) return(paste0("git:", tolower(sha)))

  repository <- .pefa_desc_field(desc, "Repository")
  if (!is.na(repository) && identical(toupper(repository), "CRAN")) {
    return(paste0("cran:", version))
  }
  paste0("source:", version)
}

.pefa_chr_scalar <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", what, "` must be one nonmissing, nonempty character string.",
         call. = FALSE)
  }
  x
}

## A lineage declaration is a set of K values, so it is stored as a bare integer
## vector: as.integer() also strips names, which would otherwise leak into the
## lineage hash and make two identical declarations disagree.
.pefa_k_vector <- function(x, what) {
  if (!is.numeric(x) || is.complex(x) || anyNA(x) || !all(is.finite(x)) ||
      !all(x == floor(x)) || !all(abs(x) <= .Machine$integer.max)) {
    stop("`", what, "` must be a vector of finite whole numbers.",
         call. = FALSE)
  }
  as.integer(x)
}

## A root object declares no parent. Plain NA is normalized to NA_character_
## because the two serialize to different bytes and would give a root object two
## possible lineage hashes.
.pefa_parent_hash <- function(x) {
  if (length(x) == 1L && is.na(x)) return(NA_character_)
  .pefa_chr_scalar(x, "parent_evidence_sha256")
}

## $sweep$secs is elapsed time: the one sweep field that a byte-identical rerun
## on a different machine may legitimately change. Drop the element rather than
## normalizing it, so evidence identity never depends on machine speed. Written
## against the list interface, which the sweep data frame satisfies.
.pefa_drop_secs <- function(sweep_table) {
  if (is.list(sweep_table) && "secs" %in% names(sweep_table)) {
    sweep_table[["secs"]] <- NULL
  }
  sweep_table
}

## Assemble the twelve-field provenance list of plan section 3.6, in that exact
## order. `evidence` is the assembled list(sweep, transitions, persistence,
## loadings, pips, Q0, settings).
##
## evidence_sha256 covers the scientific payload and the compatibility identity
## and nothing else: it excludes $sweep$secs, itself, and the three lineage
## fields, so a fresh root and a verified child that fit the same numbers share
## it. lineage_sha256 then protects the parent/reused/new-K declaration
## separately, excluding only itself. package_build_id is deliberately outside
## both: it is platform metadata, and hashing it would make an identical
## analysis look different when installed from a different source.
.pefa_provenance <- function(evidence, standardized_data_sha256, Q0_sha256,
                             settings_sha256, package_version, newly_fitted_K,
                             reused_K = integer(0),
                             parent_evidence_sha256 = NA_character_) {
  if (!is.list(evidence)) {
    stop("`evidence` must be the assembled evidence list.", call. = FALSE)
  }
  required <- c("sweep", "transitions", "persistence", "loadings", "pips",
                "Q0", "settings")
  absent <- required[!required %in% names(evidence)]
  if (length(absent)) {
    ## An absent field would hash as NULL and collide with a present-but-NULL
    ## one, so refuse rather than stamp a hash over an incomplete payload.
    stop("`evidence` is missing required field(s): ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }

  standardized_data_sha256 <- .pefa_chr_scalar(standardized_data_sha256,
                                                "standardized_data_sha256")
  Q0_sha256 <- .pefa_chr_scalar(Q0_sha256, "Q0_sha256")
  settings_sha256 <- .pefa_chr_scalar(settings_sha256, "settings_sha256")
  package_version <- .pefa_chr_scalar(package_version, "package_version")
  newly_fitted_K <- .pefa_k_vector(newly_fitted_K, "newly_fitted_K")
  reused_K <- .pefa_k_vector(reused_K, "reused_K")
  parent_evidence_sha256 <- .pefa_parent_hash(parent_evidence_sha256)

  evidence_payload <- list(
    sweep                    = .pefa_drop_secs(evidence[["sweep"]]),
    transitions              = evidence[["transitions"]],
    persistence              = evidence[["persistence"]],
    loadings                 = evidence[["loadings"]],
    pips                     = evidence[["pips"]],
    Q0                       = evidence[["Q0"]],
    settings                 = evidence[["settings"]],
    evidence_schema_id       = .pefa_evidence_schema_id,
    engine_semantics_id      = .pefa_engine_semantics_id,
    package_version          = package_version,
    standardized_data_sha256 = standardized_data_sha256,
    Q0_sha256                = Q0_sha256,
    settings_sha256          = settings_sha256
  )
  evidence_sha256 <- .pefa_sha256(evidence_payload)

  lineage_sha256 <- .pefa_sha256(list(
    evidence_sha256        = evidence_sha256,
    parent_evidence_sha256 = parent_evidence_sha256,
    reused_K               = reused_K,
    newly_fitted_K         = newly_fitted_K
  ))

  list(
    evidence_schema_id       = .pefa_evidence_schema_id,
    engine_semantics_id      = .pefa_engine_semantics_id,
    package_version          = package_version,
    package_build_id         = .pefa_build_id(),
    standardized_data_sha256 = standardized_data_sha256,
    Q0_sha256                = Q0_sha256,
    settings_sha256          = settings_sha256,
    evidence_sha256          = evidence_sha256,
    lineage_sha256           = lineage_sha256,
    parent_evidence_sha256   = parent_evidence_sha256,
    reused_K                 = reused_K,
    newly_fitted_K           = newly_fitted_K
  )
}
