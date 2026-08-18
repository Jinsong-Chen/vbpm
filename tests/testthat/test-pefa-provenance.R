## Immutable provenance (section 3.6): twelve fields in a fixed order, hashes
## that bind the scientific payload and the compatibility identity, and a
## lineage declaration protected separately from it.

.prov_data <- function(N = 50, J = 6, seed = 8) {
  set.seed(seed)
  group <- rep(1:2, each = J / 2)
  loading <- matrix(0, J, 2)
  loading[cbind(seq_len(J), group)] <- .8
  Y <- matrix(rnorm(N * 2), N, 2) %*% t(loading) +
    matrix(rnorm(N * J, sd = sqrt(1 - .64)), N, J)
  Q0 <- matrix(-1L, J, 1)
  Q0[1:2, 1] <- 1L
  list(Y = Y, Q0 = Q0)
}

.prov_fit <- function(dat, ...) {
  pefa(dat$Q0, dat$Y, 1, 2, v0 = .001, max_it = 500, verbose = FALSE, ...)
}

.prov_evidence <- function(x) {
  x[c("sweep", "transitions", "persistence", "loadings", "pips", "Q0",
      "settings")]
}

.prov_dat <- .prov_data()
.prov_first <- .prov_fit(.prov_dat)

test_that("provenance holds the twelve fields in order with root lineage", {
  p <- .prov_first$provenance
  expect_type(p, "list")
  expect_identical(
    names(p),
    c("evidence_schema_id", "engine_semantics_id", "package_version",
      "package_build_id", "standardized_data_sha256", "Q0_sha256",
      "settings_sha256", "evidence_sha256", "lineage_sha256",
      "parent_evidence_sha256", "reused_K", "newly_fitted_K")
  )
  expect_identical(p$evidence_schema_id, "vbpm_pefa_evidence_1")
  expect_identical(p$engine_semantics_id, "vbpm_pefa_engine_0.9")
  expect_identical(p$package_version,
                   as.character(utils::packageVersion("vbpm")))
  expect_match(p$package_build_id, "^(git:|cran:|source:)")

  ## The first nine fields are nonmissing character scalars.
  for (field in names(p)[1:9]) {
    value <- p[[field]]
    expect_type(value, "character")
    expect_length(value, 1L)
    expect_false(is.na(value), info = field)
    expect_true(nzchar(value), info = field)
  }
  for (field in c("standardized_data_sha256", "Q0_sha256", "settings_sha256",
                  "evidence_sha256", "lineage_sha256")) {
    expect_match(p[[field]], "^[0-9a-f]{64}$", info = field)
  }

  ## A root object declares no parent and reuses nothing.
  expect_identical(p$parent_evidence_sha256, NA_character_)
  expect_identical(p$reused_K, integer(0))
  expect_identical(p$newly_fitted_K, 1:2)
  expect_identical(typeof(p$newly_fitted_K), "integer")

  ## No decision vocabulary may appear in provenance.
  expect_false(any(grepl(
    "profile|criterion|threshold|horizon|layer|cut|sustain|selected",
    names(p)
  )))
})

test_that("the evidence hash is reproducible and excludes elapsed time", {
  second <- .prov_fit(.prov_dat)
  a <- .prov_first$provenance
  b <- second$provenance

  ## The scientific payload is identical, so every hash agrees.
  expect_identical(b$evidence_sha256, a$evidence_sha256)
  expect_identical(b$lineage_sha256, a$lineage_sha256)
  expect_identical(b$standardized_data_sha256, a$standardized_data_sha256)
  expect_identical(b$Q0_sha256, a$Q0_sha256)
  expect_identical(b$settings_sha256, a$settings_sha256)
  ## ... even though `secs` is a runtime measurement rather than evidence.
  timeless <- setdiff(names(.prov_first$sweep), "secs")
  expect_identical(second$sweep[, timeless], .prov_first$sweep[, timeless])
  expect_true(all(is.finite(second$sweep$secs)))

  ## Recomputing provenance from the stored evidence reproduces the object.
  recomputed <- vbpm:::.pefa_provenance(
    evidence = .prov_evidence(.prov_first),
    standardized_data_sha256 = vbpm:::.pefa_standardized_hash(.prov_dat$Y),
    Q0_sha256 = vbpm:::.pefa_sha256(.prov_first$Q0),
    settings_sha256 = vbpm:::.pefa_sha256(.prov_first$settings),
    package_version = as.character(utils::packageVersion("vbpm")),
    newly_fitted_K = 1:2
  )
  expect_identical(recomputed, a)

  ## Perturbing only the elapsed time leaves the evidence hash alone.
  slower <- .prov_evidence(.prov_first)
  slower$sweep$secs <- slower$sweep$secs + 1000
  expect_identical(
    vbpm:::.pefa_provenance(
      evidence = slower,
      standardized_data_sha256 = a$standardized_data_sha256,
      Q0_sha256 = a$Q0_sha256, settings_sha256 = a$settings_sha256,
      package_version = a$package_version, newly_fitted_K = 1:2
    )$evidence_sha256,
    a$evidence_sha256
  )

  ## Any scientific payload mutation is detected.
  for (mutate in list(
    function(e) { e$sweep$ELBO[1L] <- 0; e },
    function(e) { e$persistence$phi_min[1L] <- -1; e },
    function(e) { e$transitions$ELBO_gain[1L] <- -999; e },
    function(e) { e$loadings[["1"]][1L, 1L] <- 42; e },
    function(e) { e$pips[["1"]][1L, 1L] <- 42; e }
  )) {
    expect_false(identical(
      vbpm:::.pefa_provenance(
        evidence = mutate(.prov_evidence(.prov_first)),
        standardized_data_sha256 = a$standardized_data_sha256,
        Q0_sha256 = a$Q0_sha256, settings_sha256 = a$settings_sha256,
        package_version = a$package_version, newly_fitted_K = 1:2
      )$evidence_sha256,
      a$evidence_sha256
    ))
  }
})

test_that("a changed setting moves its own hash and the evidence hash", {
  changed <- .prov_fit(.prov_dat, stability_eps = 0.2)
  a <- .prov_first$provenance
  b <- changed$provenance

  expect_false(identical(b$settings_sha256, a$settings_sha256))
  expect_false(identical(b$evidence_sha256, a$evidence_sha256))
  ## The data and backbone are untouched, so their hashes must not move.
  expect_identical(b$standardized_data_sha256, a$standardized_data_sha256)
  expect_identical(b$Q0_sha256, a$Q0_sha256)
  expect_identical(b$evidence_schema_id, a$evidence_schema_id)
  expect_identical(b$engine_semantics_id, a$engine_semantics_id)
})

test_that("the input hashes bind the data, its order, and the backbone", {
  Y <- .prov_dat$Y
  base <- vbpm:::.pefa_standardized_hash(Y)
  expect_identical(vbpm:::.pefa_standardized_hash(Y), base)

  moved <- Y
  moved[1L, 1L] <- moved[1L, 1L] + 1e-6
  expect_false(identical(vbpm:::.pefa_standardized_hash(moved), base))
  expect_false(identical(vbpm:::.pefa_standardized_hash(Y[, c(2:6, 1L)]),
                         base))
  named <- Y
  colnames(named) <- paste0("v", seq_len(ncol(named)))
  expect_false(identical(vbpm:::.pefa_standardized_hash(named), base))
  ## Missingness is part of the payload, not silently dropped.
  gappy <- Y
  gappy[1L, 1L] <- NA_real_
  expect_false(identical(vbpm:::.pefa_standardized_hash(gappy), base))

  Q0 <- .prov_first$Q0
  expect_identical(vbpm:::.pefa_sha256(Q0), .prov_first$provenance$Q0_sha256)
  other <- Q0
  other[3L, 1L] <- 0L
  expect_false(identical(vbpm:::.pefa_sha256(other),
                         .prov_first$provenance$Q0_sha256))
})

test_that("signed zero and source encoding cannot move a hash", {
  h <- vbpm:::.pefa_sha256

  ## -0 and 0 are `identical()` and `==` but serialize to different bytes, so
  ## before normalization an unchanged payload carrying one looked mutated.
  expect_true(identical(0, -0))
  expect_identical(h(0), h(-0))
  expect_identical(h(c(1, 0, 3)), h(c(1, -0, 3)))
  expect_identical(h(matrix(c(0, 1, 2, 3), 2)), h(matrix(c(-0, 1, 2, 3), 2)))

  ## NA, NaN and both infinities pass through the normalization untouched.
  canonical <- vbpm:::.pefa_canonical(c(NA, NaN, Inf, -Inf, 1, -0))
  expect_true(is.na(canonical[1L]) && !is.nan(canonical[1L]))
  expect_true(is.nan(canonical[2L]))
  expect_identical(canonical[3L], Inf)
  expect_identical(canonical[4L], -Inf)
  expect_identical(canonical[5L], 1)
  expect_identical(1 / canonical[6L], Inf)   # the -0 really did become +0
  expect_false(identical(h(NaN), h(NA_real_)))

  ## The same text in latin1 and in UTF-8 is one value to R -- comparison
  ## translates -- but serialize() writes the declared encoding and the source
  ## bytes.  Values, names, and dimnames must all normalize.
  latin <- "caf\xe9"
  Encoding(latin) <- "latin1"
  utf8 <- enc2utf8(latin)
  expect_true(identical(latin, utf8))
  expect_false(identical(Encoding(latin), Encoding(utf8)))
  expect_identical(h(latin), h(utf8))

  named_a <- c(1, 2)
  named_b <- c(1, 2)
  names(named_a) <- c(latin, "b")
  names(named_b) <- c(utf8, "b")
  expect_true(identical(named_a, named_b))
  expect_identical(h(named_a), h(named_b))

  dim_a <- matrix(1:4, 2, dimnames = list(c(latin, "r2"), c("c1", latin)))
  dim_b <- matrix(1:4, 2, dimnames = list(c(utf8, "r2"), c("c1", utf8)))
  expect_true(identical(dim_a, dim_b))
  expect_identical(h(dim_a), h(dim_b))

  ## Nested one level deeper, which is how the evidence payload carries its
  ## matrices, and with a signed zero in the same object.
  expect_identical(h(list(a = list(b = dim_a, s = latin, z = -0))),
                   h(list(a = list(b = dim_b, s = utf8, z = 0))))

  ## The normalization is not a blunt instrument: genuinely different values,
  ## names, and dimnames still hash apart.
  expect_false(identical(h(latin), h("cafe")))
  expect_false(identical(h(c(1, 0, 3)), h(c(1, 1e-300, 3))))
  renamed <- named_a
  names(renamed) <- c("a", "b")
  expect_false(identical(h(named_a), h(renamed)))
})

test_that("the lineage hash protects only the lineage declaration", {
  a <- .prov_first$provenance
  child <- vbpm:::.pefa_provenance(
    evidence = .prov_evidence(.prov_first),
    standardized_data_sha256 = a$standardized_data_sha256,
    Q0_sha256 = a$Q0_sha256, settings_sha256 = a$settings_sha256,
    package_version = a$package_version,
    newly_fitted_K = 2L, reused_K = 1L,
    parent_evidence_sha256 = a$evidence_sha256
  )
  ## The same science, a different declaration.
  expect_identical(child$evidence_sha256, a$evidence_sha256)
  expect_false(identical(child$lineage_sha256, a$lineage_sha256))
  expect_identical(child$reused_K, 1L)
  expect_identical(child$newly_fitted_K, 2L)
  expect_identical(child$parent_evidence_sha256, a$evidence_sha256)

  ## An incomplete evidence payload is refused rather than stamped.
  expect_error(
    vbpm:::.pefa_provenance(
      evidence = .prov_evidence(.prov_first)[1:3],
      standardized_data_sha256 = a$standardized_data_sha256,
      Q0_sha256 = a$Q0_sha256, settings_sha256 = a$settings_sha256,
      package_version = a$package_version, newly_fitted_K = 1:2
    ),
    "missing required field"
  )
})
