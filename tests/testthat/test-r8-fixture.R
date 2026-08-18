## The R8 integration fixture (plan section 6.1 and release-plan section 8.1).
##
## inst/extdata/vbpm_r8_fixture.rds holds `$members`, a named list of sweeps
## whose loading and PIP matrices were built by hand, never fitted, and
## `$integrity`, one complete sealed evidence object carrying the
## discriminating max-congruence case.  Re-running the matcher and the two
## table builders on those stored matrices therefore exercises only
## R/pefa_match.R, R/pefa.R's builders, the persistence() pivot, and
## verify_pefa() -- vbfa() is never called, so the frozen expectations are
## BLAS-stable.
##
## Comparison policy, as stated in section 6.1 and in ?pefa: identical() for
## integers, logicals, status/reason character vectors, `r`, dimnames, storage
## modes, and the NA pattern; expect_equal(tolerance = 1e-10) for the double
## metric cells.
##
## data-raw/make_r8_fixture.R rebuilds the file deterministically.

r8_fixture <- function() {
  path <- system.file("extdata", "vbpm_r8_fixture.rds", package = "vbpm")
  if (!nzchar(path) || !file.exists(path)) {
    skip("inst/extdata/vbpm_r8_fixture.rds is not installed.")
  }
  readRDS(path)
}

## Rebuild both tables from a member's stored matrices.
r8_rebuild <- function(member) {
  list(
    transitions = vbpm:::.pefa_transitions_table(member$sweep),
    persistence = vbpm:::.pefa_persistence_table(
      member$sweep, member$loadings, member$pips, member$Q0, member$K0,
      member$stability_eps, member$bifactor
    )
  )
}

## persistence() needs a "pefa" object; only $sweep$K and $persistence are read.
r8_stub <- function(member, persistence_table) {
  structure(
    list(sweep = member$sweep, transitions = member$expected$transitions,
         persistence = persistence_table, loadings = member$loadings,
         pips = member$pips, Q0 = member$Q0,
         settings = list(bifactor = member$bifactor,
                         stability_eps = member$stability_eps),
         provenance = list()),
    class = "pefa"
  )
}

r8_row <- function(pe, from, to) {
  hit <- which(pe$K_from == from & pe$K_to == to)
  expect_length(hit, 1L)
  pe[hit, , drop = FALSE]
}

## identical() everywhere except the double metric cells, which get 1e-10.
expect_frozen_table <- function(observed, expected, double_cols, label) {
  expect_identical(names(observed), names(expected), info = label)
  expect_identical(vapply(observed, typeof, ""), vapply(expected, typeof, ""),
                   info = label)
  expect_identical(nrow(observed), nrow(expected), info = label)
  expect_identical(rownames(observed), rownames(expected), info = label)
  for (nm in names(expected)) {
    tag <- paste0(label, "$", nm)
    if (nm %in% double_cols) {
      ## The NA pattern is checked exactly even where the value is not.
      expect_identical(is.na(observed[[nm]]), is.na(expected[[nm]]),
                       info = tag)
      expect_equal(observed[[nm]], expected[[nm]], tolerance = 1e-10,
                   info = tag)
    } else {
      ## identical() already pins the NA pattern for these.
      expect_identical(observed[[nm]], expected[[nm]], info = tag)
    }
  }
}

expect_frozen_pivot <- function(observed, expected, label) {
  expect_identical(class(observed), class(expected), info = label)
  expect_identical(dim(observed), dim(expected), info = label)
  expect_identical(dimnames(observed), dimnames(expected), info = label)
  expect_identical(typeof(observed), typeof(expected), info = label)
  got <- unclass(observed)
  want <- unclass(expected)
  attributes(got) <- attributes(got)[c("dim", "dimnames")]
  attributes(want) <- attributes(want)[c("dim", "dimnames")]
  expect_identical(is.na(got), is.na(want), info = label)
  if (typeof(observed) == "integer") {
    expect_identical(got, want, info = label)
  } else {
    expect_equal(got, want, tolerance = 1e-10, info = label)
  }
  for (a in c("collision", "reduced", "fit_ok")) {
    expect_identical(attr(observed, a, exact = TRUE),
                     attr(expected, a, exact = TRUE),
                     info = paste0(label, " attr:", a))
  }
}


test_that("the R8 fixture ships every member section 6.1 requires", {
  fx <- r8_fixture()

  expect_identical(fx$fixture_version, "vbpm_r8_fixture_1")
  expect_identical(fx$evidence_schema_id, "vbpm_pefa_evidence_1")
  ## The file holds the low-level members and the one sealed integrity case,
  ## and nothing else: an unannounced component would travel to R8 unread.
  expect_identical(names(fx),
                   c("fixture_version", "evidence_schema_id", "metrics",
                     "double_pair_cols", "members", "integrity"))
  expect_identical(fx$metrics,
                   c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
                     "unmatched_ssl", "n_collisions"))
  ## The pivot enumeration the fixture freezes must be persistence()'s own.
  expect_identical(fx$metrics, eval(formals(persistence)$metric))

  expect_identical(
    names(fx$members),
    c("all_ok", "nonconverged_endpoint", "malformed_loading",
      "degenerate_backbone", "bifactor_thin_backbone", "ineligible_source",
      "scarce_targets", "k0_empty_pairs", "k0_nonempty_pairs")
  )

  ## Every feature section 6.1 enumerates is claimed by some member.
  features <- unlist(lapply(fx$members, `[[`, "features"), use.names = FALSE)
  expect_true(all(c("k_equals_K0_edge", "collision_edge",
                    "nonconverged_endpoint", "malformed_loading_endpoint",
                    "degenerate_backbone", "bifactor", "thin_backbone",
                    "backbone_ssl_min_offset", "ineligible_source",
                    "scarce_targets", "empty_pair_set", "nonempty_pair_set")
                  %in% features))

  ## The members must be evidence, not fits: nothing carries a fit object.
  for (nm in names(fx$members)) {
    m <- fx$members[[nm]]
    expect_identical(
      sort(names(m)),
      sort(c("description", "features", "K0", "bifactor", "stability_eps",
             "Q0", "sweep", "loadings", "pips", "expected"))
    )
    expect_identical(sort(names(m$expected)),
                     c("persistence", "pivots", "transitions"))
    expect_identical(names(m$loadings), as.character(m$sweep$K))
    expect_identical(names(m$pips), names(m$loadings))
    expect_true(is.integer(m$Q0))
  }
})


test_that("the matcher reproduces every frozen $persistence table", {
  fx <- r8_fixture()
  double_cols <- fx$double_pair_cols
  for (nm in names(fx$members)) {
    member <- fx$members[[nm]]
    pe <- r8_rebuild(member)$persistence
    expect_frozen_table(pe, member$expected$persistence, double_cols,
                        paste0(nm, ":persistence"))
    ## Section 3.3: every ordered pair present, ordered, with r = K_to - K_from.
    W <- nrow(member$sweep)
    expect_identical(nrow(pe), as.integer(W * (W - 1L) / 2L))
    expect_identical(pe$r, pe$K_to - pe$K_from)
    expect_false(is.unsorted(pe$K_from))
    expect_identical(order(pe$K_from, pe$K_to), seq_len(nrow(pe)))
  }
})


test_that("the builder reproduces every frozen $transitions table", {
  fx <- r8_fixture()
  double_cols <- setdiff(names(fx$members[[1]]$expected$transitions),
                         c("K_from", "K_to"))
  for (nm in names(fx$members)) {
    member <- fx$members[[nm]]
    tr <- r8_rebuild(member)$transitions
    expect_frozen_table(tr, member$expected$transitions, double_cols,
                        paste0(nm, ":transitions"))
    expect_identical(nrow(tr), nrow(member$sweep) - 1L)
  }
})


test_that("persistence() reproduces every frozen pivoted matrix", {
  fx <- r8_fixture()
  for (nm in names(fx$members)) {
    member <- fx$members[[nm]]
    stub <- r8_stub(member, r8_rebuild(member)$persistence)
    for (metric in fx$metrics) {
      expect_frozen_pivot(persistence(stub, metric),
                          member$expected$pivots[[metric]],
                          paste0(nm, ":", metric))
    }
    labels <- as.character(member$sweep$K)
    pivot <- member$expected$pivots$phi_min
    expect_identical(dimnames(pivot), list(labels, labels))
    expect_true(all(is.na(pivot[lower.tri(pivot, diag = TRUE)])))
  }
})


test_that("the all-ok member carries a k = K0 edge and a collision edge", {
  fx <- r8_fixture()
  member <- fx$members$all_ok
  pe <- r8_rebuild(member)$persistence

  for (cm in c("fit", "loading", "pip", "elbo", "stats")) {
    expect_true(all(member$sweep[[paste0(cm, "_status")]] == "ok"))
    expect_true(all(is.na(member$sweep[[paste0(cm, "_reason")]])))
  }

  ## k = K0 > 0: the backbone-only comparison, complete for its source.
  edge <- r8_row(pe, member$K0, member$K0 + 1L)
  expect_identical(edge$K_from, member$K0)
  expect_identical(edge$n_pairs, member$K0)
  expect_identical(edge$pair_status, "available")
  expect_true(is.finite(edge$phi_min))
  expect_true(is.finite(edge$rmsd_max))
  expect_identical(edge$n_collisions, 0L)
  expect_false(edge$n_pairs < edge$K_from)   # never a reduced pair set

  ## The collision edge really collides, and keeps every metric.
  hit <- pe[!is.na(pe$n_collisions) & pe$n_collisions > 0L, , drop = FALSE]
  expect_identical(nrow(hit), 1L)
  expect_identical(c(hit$K_from, hit$K_to), c(4L, 5L))
  expect_true(hit$collision)
  expect_identical(hit$n_collisions, 1L)
  expect_identical(hit$collision_targets, "4")
  expect_identical(hit$collision_multiplicities, "2")
  expect_identical(hit$pair_status, "available")
  expect_true(all(is.finite(c(hit$phi_min, hit$rmsd, hit$rmsd_max, hit$ari,
                              hit$pip_rmsd, hit$unmatched_ssl))))
  expect_lte(hit$rmsd, hit$rmsd_max)
  ## Two sources on one target is one collision, not two.
  expect_identical(hit$n_pairs, 4L)
  expect_identical(hit$unmatched_n, 2L)

  ## The pivot's "collision" attribute marks exactly that cell.
  flag <- attr(member$expected$pivots$phi_min, "collision", exact = TRUE)
  expect_identical(which(flag %in% TRUE), which(flag == TRUE))
  expect_identical(sum(flag, na.rm = TRUE), 1L)
  expect_true(flag["4", "5"])
})


test_that("the nonconverged member keeps every metric and only flags it", {
  fx <- r8_fixture()
  member <- fx$members$nonconverged_endpoint
  sweep <- member$sweep
  pe <- r8_rebuild(member)$persistence

  bad <- which(sweep$converged %in% FALSE)
  expect_identical(bad, 2L)
  expect_identical(sweep$fit_status[bad], "unavailable")
  expect_identical(sweep$fit_reason[bad], "nonconverged")
  ## Components are independent: the payload survives the fit status.
  for (cm in c("loading", "pip", "elbo", "stats")) {
    expect_identical(sweep[[paste0(cm, "_status")]][bad], "ok")
  }
  expect_true(all(is.finite(member$loadings[[bad]])))
  expect_true(is.finite(sweep$backbone_ssl_min[bad]))

  ## Every pair, including the two touching K = 3, is fully computed.
  expect_true(all(pe$pair_status == "available"))
  expect_true(all(is.na(pe$pair_reason)))
  metrics <- c("phi_min", "rmsd", "rmsd_max", "ari", "pip_rmsd",
               "unmatched_ssl")
  for (nm in metrics) expect_true(all(is.finite(pe[[nm]])), info = nm)
  expect_false(anyNA(pe$n_pairs))
  expect_false(anyNA(pe$collision))

  ## The only trace of nonconvergence is from_fit_ok / to_fit_ok.
  expect_identical(pe$from_fit_ok, c(TRUE, TRUE, FALSE))
  expect_identical(pe$to_fit_ok, c(FALSE, TRUE, TRUE))
  fit_ok <- attr(member$expected$pivots$phi_min, "fit_ok", exact = TRUE)
  expect_false(fit_ok["2", "3"])
  expect_true(fit_ok["2", "4"])
  expect_false(fit_ok["3", "4"])
})


test_that("the malformed-loading member is the only unavailable pair", {
  fx <- r8_fixture()
  member <- fx$members$malformed_loading
  sweep <- member$sweep
  pe <- r8_rebuild(member)$persistence

  bad <- which(sweep$loading_status != "ok")
  expect_identical(bad, 2L)
  expect_identical(sweep$loading_reason[bad], "malformed_loading")
  expect_identical(sweep$backbone_ssl_min[bad], NA_real_)
  expect_true(all(is.na(member$loadings[[bad]])))
  ## Its PIP and its fit are independently fine, so no other status leaks in.
  expect_identical(sweep$pip_status[bad], "ok")
  expect_identical(sweep$fit_status[bad], "ok")
  expect_true(sweep$converged[bad])

  expect_identical(pe$pair_status,
                   c("unavailable", "available", "unavailable"))
  expect_identical(pe$pair_reason,
                   c("malformed_loading", NA_character_, "malformed_loading"))
  ## malformed_loading is the sole reason code in 0.9.0.
  expect_identical(unique(stats::na.omit(pe$pair_reason)), "malformed_loading")
  ## Convergence never drives availability.
  expect_true(all(pe$from_fit_ok) && all(pe$to_fit_ok))

  ## Each eligibility count survives wherever its own endpoint survives.
  malformed_target <- r8_row(pe, 3L, 4L)
  expect_identical(malformed_target$source_ineligible_n, 0L)
  expect_identical(malformed_target$target_eligible_n, NA_integer_)
  malformed_source <- r8_row(pe, 4L, 5L)
  expect_identical(malformed_source$source_ineligible_n, NA_integer_)
  expect_identical(malformed_source$target_eligible_n, 3L)
  for (row in list(malformed_target, malformed_source)) {
    expect_identical(row$n_pairs, NA_integer_)
    expect_identical(row$collision, NA)
    expect_identical(row$collision_targets, NA_character_)
    expect_identical(row$backbone_degenerate, NA)
  }

  ## Its ELBO is separately unavailable, so that whole gain path is nonfinite:
  ## the percent column is all NA and no max(numeric(0)) warning escapes.
  expect_identical(sweep$elbo_status[bad], "unavailable")
  tr <- expect_no_warning(vbpm:::.pefa_transitions_table(sweep))
  expect_true(all(is.na(tr$ELBO_gain)))
  expect_true(all(is.na(tr$ELBO_gain_pct)))
  expect_true(all(is.finite(tr$BIC_gain_pct)))
})


test_that("the degenerate backbone loses only the congruence", {
  fx <- r8_fixture()
  member <- fx$members$degenerate_backbone
  sweep <- member$sweep
  pe <- r8_rebuild(member)$persistence

  deg <- 2L
  expect_identical(sweep$backbone_ssl_min[deg], 0)
  expect_identical(sweep$loading_status[deg], "ok")   # finite, well shaped
  expect_true(all(is.finite(member$loadings[[deg]])))
  expect_identical(sum(member$loadings[[deg]][, 2L]^2), 0)

  expect_identical(pe$backbone_degenerate, c(TRUE, FALSE, TRUE))
  touching <- pe[pe$backbone_degenerate, , drop = FALSE]
  expect_identical(touching$pair_status, rep("available", 2L))
  expect_identical(touching$pair_reason, rep(NA_character_, 2L))
  expect_identical(touching$phi_min, rep(NA_real_, 2L))
  ## Aligned RMSD, ARI, PIP, collision and unmatched evidence are retained.
  for (nm in c("rmsd", "rmsd_max", "ari", "pip_rmsd", "unmatched_ssl")) {
    expect_true(all(is.finite(touching[[nm]])), info = nm)
  }
  expect_identical(touching$collision, c(FALSE, FALSE))
  expect_identical(touching$n_collisions, c(0L, 0L))

  clean <- r8_row(pe, 2L, 4L)
  expect_false(clean$backbone_degenerate)
  expect_true(is.finite(clean$phi_min))

  ## The printed triangle marks the lost congruence with "?" and keeps its
  ## computed neighbour.
  out <- utils::capture.output(print(member$expected$pivots$phi_min))
  expect_true(any(grepl("?", out, fixed = TRUE)))
  expect_true(any(grepl("? = unavailable", out, fixed = TRUE)))
})


test_that("the bifactor member offsets backbone_ssl_min past the general", {
  fx <- r8_fixture()
  member <- fx$members$bifactor_thin_backbone
  sweep <- member$sweep
  pe <- r8_rebuild(member)$persistence

  expect_true(member$bifactor)
  expect_identical(member$K0, 1L)
  scales <- ssl(r8_stub(member, pe))
  expect_identical(names(scales), names(member$loadings))

  for (i in seq_along(sweep$K)) {
    s <- scales[[i]]
    ## K + 1 stored columns, general first.
    expect_length(s, sweep$K[i] + 1L)
    right <- min(s[seq_len(member$K0) + 1L])
    wrong <- min(s[seq_len(member$K0)])
    expect_equal(sweep$backbone_ssl_min[i], right, tolerance = 1e-10)
    ## A leaked general column would have produced `wrong`.
    expect_false(isTRUE(all.equal(right, wrong)))
    expect_lt(wrong, right)
    expect_identical(
      vbpm:::.pefa_backbone_ssl_min(member$loadings[[i]], member$K0,
                                    TRUE, TRUE),
      sweep$backbone_ssl_min[i]
    )
    ## Small but nonzero: below the screen the backbone is never subject to.
    expect_lt(right, member$stability_eps)
    expect_gt(right, 0)
  }

  ## The thin backbone stays paired: not degenerate, congruence still defined.
  expect_true(all(!pe$backbone_degenerate))
  expect_true(all(is.finite(pe$phi_min)))
  backbone_only <- r8_row(pe, member$K0, member$K0 + 1L)
  expect_identical(backbone_only$n_pairs, member$K0)
  expect_true(is.finite(backbone_only$phi_min))

  ## Group-block indices: the reused target sits at group 2 = stored 3.
  hit <- r8_row(pe, 3L, 4L)
  expect_true(hit$collision)
  expect_identical(hit$n_collisions, 1L)
  expect_identical(hit$collision_targets, "2")
  expect_identical(hit$collision_multiplicities, "2")
  reused <- as.integer(strsplit(hit$collision_targets, ",", fixed = TRUE)[[1]])
  expect_identical(reused, 2L)
  expect_identical(ncol(member$loadings[["4"]]), 5L)  # G + 4 group columns
  expect_identical(colnames(member$loadings[["4"]])[1L], "G")
  ## Group column 2 is stored column 3 and is the second group factor; a
  ## leaked general-column offset would have reported group index 3 instead.
  stored <- reused + 1L
  expect_identical(stored, 3L)
  expect_identical(colnames(member$loadings[["4"]])[stored], "F2")
})


test_that("an ineligible source reduces the pair set without failing", {
  fx <- r8_fixture()
  member <- fx$members$ineligible_source
  pe <- r8_rebuild(member)$persistence

  hit <- r8_row(pe, 3L, 4L)
  expect_identical(hit$source_ineligible_n, 1L)
  expect_identical(hit$pair_status, "available")
  expect_identical(hit$pair_reason, NA_character_)
  expect_identical(hit$n_pairs, member$K0)
  expect_lt(hit$n_pairs, hit$K_from)              # a reduced pair set
  expect_true(all(is.finite(c(hit$phi_min, hit$rmsd, hit$rmsd_max, hit$ari,
                              hit$pip_rmsd, hit$unmatched_ssl))))
  expect_false(hit$collision)

  ## Section 2.7: a screened-out column is still an unmatched target.
  screened <- r8_row(pe, 2L, 3L)
  expect_identical(screened$target_eligible_n, 0L)
  expect_identical(screened$unmatched_n, 1L)
  expect_lt(screened$unmatched_ssl, member$stability_eps)

  reduced <- attr(member$expected$pivots$phi_min, "reduced", exact = TRUE)
  expect_true(reduced["3", "4"])
  expect_false(reduced["2", "3"])
  out <- utils::capture.output(print(member$expected$pivots$phi_min))
  expect_true(any(grepl("~ = reduced pair set", out, fixed = TRUE)))
})


test_that("scarce eligible targets leave a backbone-only measurement", {
  fx <- r8_fixture()
  member <- fx$members$scarce_targets
  pe <- r8_rebuild(member)$persistence

  hit <- r8_row(pe, 3L, 4L)
  expect_identical(hit$target_eligible_n, 0L)
  expect_identical(hit$source_ineligible_n, 0L)   # the source was eligible
  expect_identical(hit$n_pairs, member$K0)
  expect_identical(hit$pair_status, "available")
  expect_false(hit$collision)
  expect_identical(hit$n_collisions, 0L)
  expect_identical(hit$collision_targets, "")
  expect_identical(hit$collision_multiplicities, "")
  expect_true(all(is.finite(c(hit$phi_min, hit$rmsd, hit$rmsd_max, hit$ari,
                              hit$pip_rmsd))))
  ## Both screened-out targets still count as unmatched.
  expect_identical(hit$unmatched_n, 2L)
  expect_true(is.finite(hit$unmatched_ssl))

  ## The internal per-column record holds typed NA for the unplaceable source.
  facts <- vbpm:::.pefa_pair_facts(
    lam_from = member$loadings[["3"]], lam_to = member$loadings[["4"]],
    pip_from = member$pips[["3"]], pip_to = member$pips[["4"]],
    reg_from = matrix(TRUE, nrow(member$Q0), 3L),
    reg_to = matrix(TRUE, nrow(member$Q0), 4L),
    K0 = member$K0, stability_eps = member$stability_eps,
    loading_ok_from = TRUE, loading_ok_to = TRUE,
    pip_ok_from = TRUE, pip_ok_to = TRUE
  )
  expect_identical(facts$n_pairs, hit$n_pairs)
  expect_identical(sum(!is.na(facts$.assign$target)), facts$n_pairs)
  expect_identical(facts$.assign$target[3L], NA_integer_)
  expect_true(facts$.assign$eligible_source[3L])
})


test_that("K0 = 0 gives both an empty and a nonempty scored pair set", {
  fx <- r8_fixture()

  empty <- fx$members$k0_empty_pairs
  pe <- r8_rebuild(empty)$persistence
  expect_identical(empty$K0, 0L)
  expect_identical(ncol(empty$Q0), 0L)
  expect_true(all(is.na(empty$sweep$backbone_ssl_min)))
  hit <- r8_row(pe, 1L, 2L)
  expect_identical(hit$n_pairs, 0L)
  for (nm in c("phi_min", "rmsd", "rmsd_max", "pip_rmsd", "ari")) {
    expect_identical(hit[[nm]], NA_real_, info = nm)
  }
  expect_identical(hit$collision, FALSE)          # a known clean zero
  expect_identical(hit$n_collisions, 0L)
  expect_identical(hit$collision_targets, "")
  expect_identical(hit$collision_multiplicities, "")
  expect_identical(hit$backbone_degenerate, FALSE)
  expect_identical(hit$source_ineligible_n, 1L)
  expect_identical(hit$unmatched_n, 2L)           # still computed
  expect_true(is.finite(hit$unmatched_ssl))
  expect_identical(hit$pair_status, "available")

  nonempty <- fx$members$k0_nonempty_pairs
  pe <- r8_rebuild(nonempty)$persistence
  expect_identical(nonempty$K0, 0L)
  expect_true(all(pe$n_pairs > 0L))
  expect_true(all(is.finite(c(pe$phi_min, pe$rmsd, pe$rmsd_max, pe$ari,
                              pe$pip_rmsd))))
  expect_true(all(!pe$backbone_degenerate))
  expect_true(all(pe$pair_status == "available"))
})


test_that("the frozen tables are stable across repeat matcher runs", {
  fx <- r8_fixture()
  for (nm in names(fx$members)) {
    member <- fx$members[[nm]]
    first <- r8_rebuild(member)$persistence
    second <- r8_rebuild(member)$persistence
    expect_identical(first, second, info = nm)
  }
})


## ---------------------------------------------------------------------------
## The section 8.1 integrity case.
##
## Everything above re-runs the matcher on hand-built matrices.  The integrity
## case is the other half of the fixture: one complete, classed, sealed pefa
## object that the INSTALLED package must accept, whose two-candidate geometry
## decides between the 0.9.0 max-absolute-congruence rule and the rule 0.8.3
## retracted.  The object is constructed and re-sealed rather than fitted --
## no variational fit returns four exact loading columns -- but every hash it
## carries was computed by the package's own sealing code over the stored
## content, and $Y is the real matrix whose standardization produced the
## stored data fingerprint.  data-raw/make_r8_fixture.R argues the route in
## full and re-runs this gate before it writes the file.
##
## Reached through ::: so the tests resolve whether or not NAMESPACE has been
## regenerated from the @export tag; test-pefa-verify.R asserts the export
## itself.
.r8_verify <- vbpm:::verify_pefa

r8_integrity <- function() {
  integrity <- r8_fixture()$integrity
  if (is.null(integrity)) {
    skip("The installed fixture carries no integrity case.")
  }
  integrity
}


test_that("the fixture ships the section 8.1 integrity case whole", {
  it <- r8_integrity()

  expect_identical(
    names(it),
    c("description", "route", "features", "K0", "stability_eps", "columns",
      "rules", "object", "Y", "expected", "matcher_case", "retracted_copy",
      "retracted_copy_note")
  )
  ## The fixture states its own construction: a reader must not have to infer
  ## from the numbers that the candidates were built rather than fitted.
  expect_match(it$route, "Constructed and re-sealed, not fitted", fixed = TRUE)
  expect_match(it$retracted_copy_note, "DELIBERATELY CORRUPT", fixed = TRUE)

  ## A complete current-schema object, in construction order, with the twelve
  ## provenance fields and the ten settings.
  obj <- it$object
  expect_s3_class(obj, "pefa")
  expect_identical(names(obj),
                   c("sweep", "transitions", "persistence", "loadings",
                     "pips", "Q0", "settings", "provenance"))
  expect_identical(
    names(obj$provenance),
    c("evidence_schema_id", "engine_semantics_id", "package_version",
      "package_build_id", "standardized_data_sha256", "Q0_sha256",
      "settings_sha256", "evidence_sha256", "lineage_sha256",
      "parent_evidence_sha256", "reused_K", "newly_fitted_K")
  )
  expect_identical(obj$provenance$evidence_schema_id, "vbpm_pefa_evidence_1")
  expect_identical(obj$provenance$engine_semantics_id, "vbpm_pefa_engine_0.9")
  expect_identical(obj$provenance$newly_fitted_K, 2:3)
  expect_identical(obj$provenance$reused_K, integer(0))
  expect_identical(obj$provenance$parent_evidence_sha256, NA_character_)
  expect_identical(
    names(obj$settings),
    c("bifactor", "general", "v0", "max_it", "convChk", "tolVal", "tau",
      "stability_eps", "rank_adjust", "rank_max_J")
  )
  expect_identical(obj$settings$stability_eps, 0.10)
  expect_false(obj$settings$bifactor)

  ## Canonical dimensions: J = 6 items, a one-column backbone, a two-candidate
  ## window, and the single ordered pair that window admits.
  expect_identical(dim(obj$Q0), c(6L, 1L))
  expect_identical(it$K0, ncol(obj$Q0))
  expect_true(is.integer(obj$Q0))
  expect_identical(obj$sweep$K, 2:3)
  expect_identical(names(obj$loadings), c("2", "3"))
  expect_identical(names(obj$pips), names(obj$loadings))
  expect_identical(dim(obj$loadings[["2"]]), c(6L, 2L))
  expect_identical(dim(obj$loadings[["3"]]), c(6L, 3L))
  expect_identical(dimnames(obj$loadings[["3"]]),
                   list(paste0("item_", 1:6), c("F1", "F2", "F3")))
  expect_identical(nrow(obj$persistence), 1L)
  expect_identical(nrow(obj$transitions), 1L)

  ## Y is a real matrix in the object's item order, not a placeholder.
  expect_true(is.matrix(it$Y) && is.numeric(it$Y))
  expect_identical(ncol(it$Y), nrow(obj$Q0))
  expect_identical(colnames(it$Y), rownames(obj$Q0))
  expect_false(anyNA(it$Y))

  ## The fixture's self-description is the object's own geometry, so the two
  ## can never drift into disagreeing about which case this is.
  expect_equal(unname(obj$loadings[["2"]][, 1]), it$columns$backbone,
               tolerance = 1e-15)
  expect_equal(unname(obj$loadings[["2"]][, 2]), it$columns$source,
               tolerance = 1e-15)
  expect_equal(unname(obj$loadings[["3"]][, 1]), it$columns$backbone,
               tolerance = 1e-15)
  expect_equal(unname(obj$loadings[["3"]][, 2]), it$columns$target_1,
               tolerance = 1e-15)
  expect_equal(unname(obj$loadings[["3"]][, 3]), it$columns$target_2,
               tolerance = 1e-15)
})


test_that("the integrity geometry discriminates the two rules", {
  it <- r8_integrity()
  a <- it$columns$source
  t1 <- it$columns$target_1
  t2 <- it$columns$target_2

  ## Plan section 8.1's four numbers, recomputed from the shipped columns.
  cong <- function(x, y) abs(sum(x * y)) / sqrt(sum(x^2) * sum(y^2))
  expect_equal(sum((a - t1)^2), 0.72, tolerance = 1e-12)
  expect_equal(sum((a - t2)^2), 0.095, tolerance = 1e-12)
  expect_equal(cong(a, t1), 1, tolerance = 1e-12)
  expect_equal(cong(a, t2), 0.9330078, tolerance = 1e-7)

  ## The disagreement is the point: the nearer target is the less congruent
  ## one, and every column clears the eligibility screen.
  expect_lt(sum((a - t2)^2), sum((a - t1)^2))
  expect_gt(cong(a, t1), cong(a, t2))
  expect_gte(sum(a^2), it$stability_eps)
  expect_gte(sum(t1^2), it$stability_eps)
  expect_gte(sum(t2^2), it$stability_eps)

  ## The frozen answers are the two rules' answers, far apart.
  expect_equal(it$rules$phi_min, 1, tolerance = 1e-12)
  expect_identical(it$rules$retracted_phi_min, 0.93300782)
  expect_equal(it$rules$retracted_phi_min, cong(a, t2), tolerance = 1e-8)
  expect_gt(abs(it$rules$phi_min - it$rules$retracted_phi_min), 0.06)
})


test_that("the installed package recomputes phi_min = 1, not 0.9330078", {
  fx <- r8_fixture()
  it <- r8_integrity()
  obj <- it$object

  ## Recomputed with the installed builders from the stored primitives --
  ## never read out of the frozen table.
  rebuilt <- vbpm:::.pefa_persistence_table(
    obj$sweep, obj$loadings, obj$pips, obj$Q0, ncol(obj$Q0),
    obj$settings$stability_eps, obj$settings$bifactor
  )
  expect_frozen_table(rebuilt, obj$persistence, fx$double_pair_cols,
                      "integrity:sealed")
  expect_frozen_table(rebuilt, it$expected$persistence, fx$double_pair_cols,
                      "integrity:frozen")

  row <- r8_row(rebuilt, 2L, 3L)
  expect_equal(row$phi_min, 1, tolerance = 1e-12)
  expect_gt(abs(row$phi_min - it$rules$retracted_phi_min), 0.06)
  expect_false(isTRUE(all.equal(row$phi_min, it$rules$retracted_phi_min,
                                tolerance = 1e-6)))

  ## The same choice shows in three independent places on the row: how many
  ## pairs were scored, which target was left over, and how big it is.
  expect_identical(row$n_pairs, 2L)
  expect_identical(row$unmatched_n, 1L)
  expect_equal(row$unmatched_ssl, sum(it$columns$target_2^2),
               tolerance = 1e-12)
  expect_false(isTRUE(all.equal(row$unmatched_ssl, sum(it$columns$target_1^2),
                                tolerance = 1e-6)))
  expect_identical(row$pair_status, "available")
  expect_identical(row$collision, FALSE)
  expect_false(row$backbone_degenerate)

  ## The triangle a reader prints carries the same number.
  triangle <- persistence(obj, "phi_min")
  expect_identical(dimnames(triangle), list(c("2", "3"), c("2", "3")))
  expect_equal(triangle[["2", "3"]], 1, tolerance = 1e-12)
})


test_that("the low-level discriminating call reproduces the sealed row", {
  it <- r8_integrity()
  case <- it$matcher_case

  direct <- function() {
    vbpm:::.pefa_pair_facts(
      lam_from = case$lam_from, lam_to = case$lam_to,
      pip_from = case$pip_from, pip_to = case$pip_to,
      reg_from = case$reg_from, reg_to = case$reg_to,
      K0 = case$K0, stability_eps = case$stability_eps,
      loading_ok_from = TRUE, loading_ok_to = TRUE,
      pip_ok_from = TRUE, pip_ok_to = TRUE
    )
  }
  facts <- direct()

  ## Every field, frozen: the raw matcher is covered beside the sealed object.
  expect_identical(facts[setdiff(names(facts), ".assign")], case$expected)
  expect_identical(facts$.assign, case$expected_assign)

  ## The assignment is the whole discriminator. Source column 2 takes target
  ## column 2 (t1, congruence 1); the retracted rule would have taken the
  ## nearer column 3.
  expect_identical(facts$.assign$target, c(1L, 2L))
  expect_identical(facts$.assign$target[2L], it$rules$max_congruence_target)
  expect_false(identical(facts$.assign$target[2L],
                         it$rules$min_distance_target))
  expect_equal(facts$.assign$congruence[2L], 1, tolerance = 1e-12)
  expect_equal(facts$phi_min, 1, tolerance = 1e-12)

  ## And it is the same number the sealed table stores.
  expect_identical(facts$phi_min, it$object$persistence$phi_min[1L])
  expect_identical(facts$phi_min, it$expected$phi_min)
  expect_identical(facts$unmatched_n, it$expected$unmatched_n)
  expect_identical(facts$n_pairs, it$expected$n_pairs)

  ## Repeat calls are identical: the tie rule leaves nothing to chance.
  expect_identical(direct(), facts)
})


test_that("the sealed integrity object verifies against the installation", {
  it <- r8_integrity()
  obj <- it$object

  expect_true(.r8_verify(obj))
  expect_true(.r8_verify(obj, it$Y))
  expect_invisible(.r8_verify(obj))

  ## The data hash is checked, not assumed: a perturbed or permuted Y is
  ## refused, so the passing call above is evidence about this matrix.
  expect_error(.r8_verify(obj, it$Y + 1),
               "`Y` does not reproduce `$provenance$standardized_data_sha256`",
               fixed = TRUE)
  expect_error(.r8_verify(obj, it$Y[, c(2:6, 1)]),
               "`Y` does not reproduce", fixed = TRUE)

  ## Verification is inert, so a fixture read twice is the same fixture.
  before <- serialize(obj, NULL)
  .r8_verify(obj)
  .r8_verify(obj, it$Y)
  expect_identical(serialize(obj, NULL), before)
})


test_that("the retracted-rule copy fails on persistence and evidence", {
  it <- r8_integrity()

  ## The shipped negative copy, and one made here from the good object: the
  ## first proves what the artifact carries, the second that the failure is a
  ## property of the value rather than of how the copy was stored.
  local_copy <- it$object
  local_copy$persistence$phi_min[1L] <- it$rules$retracted_phi_min

  for (bad in list(it$retracted_copy, local_copy)) {
    expect_identical(bad$persistence$phi_min[1L], 0.93300782)
    message <- tryCatch(.r8_verify(bad), error = conditionMessage)
    expect_true(is.character(message))
    ## Controlled: the rebuilt table disagrees, and the sealed payload hash
    ## disagrees with it. Two findings, both named, and no other complaint.
    expect_match(message, "failed 2 integrity checks", fixed = TRUE)
    expect_match(message, "`$persistence` disagrees with the table rebuilt",
                 fixed = TRUE)
    expect_match(message, "column `phi_min` differs beyond tolerance",
                 fixed = TRUE)
    expect_match(message, "`$provenance$evidence_sha256` does not match",
                 fixed = TRUE)
    ## Not a shape, schema, or engine complaint: this is a valid object
    ## carrying one wrong number.
    expect_false(grepl("evidence_schema_id", message, fixed = TRUE))
    expect_false(grepl("engine_semantics_id", message, fixed = TRUE))
    expect_false(grepl("does not have the current evidence-object shape",
                       message, fixed = TRUE))
  }

  ## The good object is untouched by any of it.
  expect_true(.r8_verify(it$object))
  expect_equal(it$object$persistence$phi_min[1L], 1, tolerance = 1e-12)
})
