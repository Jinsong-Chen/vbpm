## The R8 integration fixture (plan section 6.1).
##
## inst/extdata/vbpm_r8_fixture.rds is a named list of sweeps whose loading and
## PIP matrices were built by hand, never fitted.  Re-running the matcher and
## the two table builders on those stored matrices therefore exercises only
## R/pefa_match.R, R/pefa.R's builders, and the persistence() pivot -- vbfa()
## is never called, so the frozen expectations are BLAS-stable.
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
